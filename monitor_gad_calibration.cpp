#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TStyle.h"
#include "TColor.h"
#include "TApplication.h"
#include "TMath.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TError.h"
#include "Math/MinimizerOptions.h"
#include "TVirtualFitter.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <thread>
#include <future>
#include <cassert>
#include <sys/stat.h>  // dirname and basename
#include <sys/types.h> // for stat() test to see if file or folder

class Plotter{
	public:
	
	// members
	TChain* c_dark = nullptr;
	TChain* c_275_A = nullptr;
	TChain* c_275_B = nullptr;
	TChain* c_R_G_B = nullptr;
	TChain* c_White_385 = nullptr;
	
	std::vector<double> values;
	std::vector<double> wavelengths;
	std::vector<double> errors;
	std::vector<double> values_dark;
	std::vector<double> errors_dark;
	std::vector<double> wavelength_errors;  // make this up based on 1/2 bin width.
	Short_t yr, mon, dy, hr, mn, sc;
	
	// think these need to remain in scope, for some reason
	std::vector<double>* values_p = &values;
	std::vector<double>* errors_p = &errors;
	std::vector<double>* values_dark_p = &values_dark;
	std::vector<double>* errors_dark_p = &errors_dark;
	std::vector<double>* wavelengths_p = &wavelengths;
	
	std::map<std::string, TChain*> chains;
	std::map<std::string, int> files;
	std::map<std::string,int> offsets;
	int n_datapoints;
	int num_meas;
	int darksperloop=12;
	std::map<std::string, std::pair<double,double>> calibration_data;
	std::string purerefset;
	bool do_calibration;
	bool generate_pure;  // MakePure will not overwrite unless told to do so
	
	// for finding input data, used by LoadNewData
	std::string dirname;
	std::string expectedprefix;
	int first_entry=0;
	int max_entries=-1;
	
	// methods
	bool LoadConfig(std::string filename);
	int LoadNewData();
	int SetBranchAddresses();
	int GetNextDarkEntry(std::string name, int ledon_entry_num, bool check=true);
	
	bool CheckPath(std::string path, std::string& type);
	TGraphErrors* MakePure(std::string name, bool overwrite=false);
	bool LoadConcentrations(std::string filename);
	std::vector<double> GetCalibCurve(std::string led, std::string fitmethod);
	TF1* PureScaledPlusExtras(int led_num);
	std::pair<double,double> CalculateError(TF1* abs_func, double peak1_pos, double peak2_pos);
	bool DoRawFit(TGraphErrors* g_abs, std::pair<double,double>& peak_posns, std::pair<double,double>& peak_heights, std::pair<double,double>& peak_errs, std::string pngname="", bool draw=false);
	bool DoSimpleFit(TGraph* abs_graph, std::pair<double,double>& simple_posns, std::pair<double,double>& simple_peaks, std::pair<double,double>& simple_errs, std::string pngname="", bool plot=false);
	bool DoComplexFit(TGraphErrors* g_abs, std::pair<double,double>& peak_posns, std::pair<double,double>& peak_heights, std::pair<double,double>& peak_errs, std::pair<double,double> initial_peak_heights, std::string pngname="", bool draw=false);
	
	TF1* GetAbsFunc();
	std::map<std::string, TMultiGraph*> FitCalibrationData(std::string name, bool make_cal_curve);
	
	void SetTimeAxis(TH1* hist, long t0_seconds);
	std::vector<std::string> GetListOfFiles(std::string inputdir);
	std::string GetStdoutFromCommand(std::string cmd, int bufsize=500);
	
	int ok; // general purpose use
	
	int sample_273=0;
	int sample_276=0;
	
	int Initialise(std::string configfile);
	std::map<std::string, TMultiGraph*> Execute();
	std::map<std::string,int> dark_sub_pures_byname; // map led name to index in dark_sub_pures
};
std::vector<TGraphErrors*> dark_sub_pures; // global

bool Plotter::LoadConfig(std::string filename){
	
	// read config file keys into a std::map
	std::map<std::string,std::string> configs;
	std::cout<<"loading configuration from file "<<filename<<std::endl;
	std::ifstream config_file(filename.c_str());
	if(not config_file.is_open()){
		std::cerr<<"couldn't open config file "<<filename<<std::endl;
		return false;
	}
	
	std::string line;
	//std::cout<<"getting lines"<<std::endl;
	while(getline(config_file, line)){
		//std::cout<<"processing line "<<line<<std::endl;
		if(line.empty()) continue;
		if(line[0]=='#') continue;
		std::stringstream ss(line);
		std::string key, value;
		//std::cout<<"parsing line"<<std::endl;
		if(!(ss >> key >> value)){
			std::cerr<<"Failed to parse line "<<line<<"!"<<std::endl;
			break;
		}
		configs.emplace(std::pair<std::string,std::string>{key,value});
	}
	config_file.close();
	
	// scan the map for key variables, check everything is present
	// and set the corresponding member variables
	if(configs.count("concentrations_file")!=0){
		// if specified, this implies we're generating a calibation curve
		do_calibration=true;
		// read in the file that specifies the concentration at each measurement
		bool ok = LoadConcentrations(configs.at("concentrations_file"));
		if(not ok) return false;
		// we can use the first entry in the TChain to make a pure-water
		// LED reference trace, if one doesn't exist already
		generate_pure=true;
	} else {
		// if not specified, this implies we're making a stability plot
		do_calibration=false;
		// in this case we must have a pure reference file already
		// since it is assumed that our first measurement is not pure water
		generate_pure=false;
	}
	
	// in either case we need a filename for the pure water reference (to read or create)
	if(configs.count("pure_ref_name")==0){
		std::cerr<<"no pure reference pure_ref_name given in config file!"<<std::endl;
		return false;
	} else {
		// this string is used to form the name of the pure reference file
		// used (or generated, if in calibration mode) by MakePure
		purerefset = configs.at("pure_ref_name");
	}
	
	// get the directory of input files
	if(configs.count("input_dir")==0){
		std::cerr<<"no input_dir specified!"<<std::endl;
		return false;
	} else {
		dirname = configs.at("input_dir");
	}
	// and the file prefix to identify a subset of files in this directory
	if(configs.count("input_prefix")==0){
		std::cerr<<"no input_prefix specified!"<<std::endl;
		return false;
	} else {
		expectedprefix=configs.at("input_prefix");
	}
	
	// in case we want to analyse a subset of measurements
	// in the resulting TChain from these files
	if(configs.count("first_entry")!=0){
		first_entry=stoi(configs.at("first_entry"));
	}
	if(configs.count("max_entries")!=0){
		max_entries=stoi(configs.at("max_entries"));
	}
	
	return true;
}


int main(int argc, const char* argv[]){
	
	if(argc==0){
		std::cout<<"usage: "<<argv[0]<<" <configfile>"<<std::endl;
		return 0;
	}
	std::string configfile = argv[1];
	
	TApplication myapp("rootTApp",0,0);
	gErrorIgnoreLevel = kWarning; // decrease ROOT verbosity, it never has anything useful to say.
	
	// the complex gaussian fit "fails" (does a good job but reports invalid status)
	// this seems to be because it doesn't think it's converged. not sure why.. the fit is good.
	// following stolen from DEAPFit stuff
	// from https://root-forum.cern.ch/t/speeding-up-fitting-to-a-landau-distribution/25140/2
	// see also https://root-forum.cern.ch/t/tolerance-for-th1-fit-interface/9312/4
	//ROOT::Math::MinimizerOptions::SetDefaultStrategy(0);
	//ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(100);
	//ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100);
	//ROOT::Math::MinimizerOptions::SetDefaultTolerance(0.1);
	//ROOT::Math::MinimizerOptions::SetDefaultPrecision(1E-6);  // this sets minimiser EPS - i.e numerical accuracy
	//TVirtualFitter::SetPrecision(0.5);                        // this is equivalent to SetDefaultTolerance
	
	TCanvas* c1 = new TCanvas("c1","c1",1024,800);  // do not change the name of this canvas, closing it will terminate the program.
	
	// create plotter class
	Plotter myplotter;
	
	// read configs, create toolchains and set branch addresses
	myplotter.Initialise(configfile);
	
	// each execute call will return a map of various different plots
	std::map<std::string,TMultiGraph*> plots;
	// we'll need a canvas for each
	std::map<std::string,TCanvas*> canvases;
	
	// keep looking for new files until user closes main canvas
	int loopi=0;
	int last_num_files=0;
	bool watch=false;  // should we continually watch for new files, or just process once
	while(c1!=nullptr){
		std::cout<<"."<<std::flush;
		
		if((watch && (loopi%10000)==0) || loopi==0){
			// look for new files
			std::cout<<"\nchecking for new data"<<std::endl;
			myplotter.LoadNewData();
			std::cout<<"found "<<myplotter.files.size()<<" files"<<std::endl;
			
			// if new files found...
			if(myplotter.files.size()!=last_num_files){
				//std::cout<<"new data: deleting graphs"<<std::endl;
				// delete the old plots
				for(std::pair<const std::string, TMultiGraph*> aplot : plots){
					std::cout<<"deleting..."<<aplot.first<<std::endl;
					delete aplot.second;
					std::cout<<"...deleted"<<std::endl;
				}
				
				// re-create with the full dataset
				//std::cout<<"making new graphs"<<std::endl;
				plots = myplotter.Execute();
				
				std::cout<<"drawing new plots"<<std::endl;
				for(std::pair<const std::string, TMultiGraph*> aplot : plots){
					std::string canvasname = std::string("c_")+aplot.first;
					TCanvas* next_canv = (TCanvas*)gROOT->FindObject(canvasname.c_str());
					if(next_canv==nullptr){
						next_canv = new TCanvas(canvasname.c_str(),canvasname.c_str(),1024,800);
					}
					if(canvases.count(canvasname)!=0 && canvases.at(canvasname)!=next_canv){
						delete canvases.at(canvasname);
					}
					canvases[canvasname] = next_canv;
					next_canv->cd();
					aplot.second->Draw("AXP");
					next_canv->Modified();
					next_canv->Update();
					gSystem->ProcessEvents();
					std::string canvasfilename = canvasname+".root";
					next_canv->SaveAs(canvasfilename.c_str());
				}
				
				last_num_files = myplotter.files.size();
			}
		}
		
		gSystem->ProcessEvents();
		std::this_thread::sleep_for(std::chrono::milliseconds(100));
		c1 = (TCanvas*)gROOT->FindObject("c1");
		
		++loopi;
	}
	
	for(std::pair<const std::string, TMultiGraph*> aplot : plots){
		delete aplot.second;
		std::string canvasname = std::string("c_")+aplot.first;
		TCanvas* canv = (TCanvas*)gROOT->FindObject(canvasname.c_str());
		if(canv!=nullptr){
			delete canv;
		}
	}
	
	return 0;
}

int Plotter::LoadNewData(){
	
	// scan folder for new files
	std::cout<<"scanning "<<dirname<<" for new files beginning with prefix "<<expectedprefix<<std::endl;
	
	// for some reason ROOT doesn't like bind-mounted directories in containers???
	// instead get list of files manually
	std::vector<std::string> filestlist = GetListOfFiles(dirname);
	
	// Unfortunately if you call TChain::Add with the same file twice
	// (or patterns that match the same file twice),
	// it'll get re-added and you'll have duplicate entries.
	// so we manually keep track of what's in the chain
	// and re-scan the directory, adding only what's new
	std::map<std::string,bool> newfiles;
	if(filestlist.size()){
		//std::cout<<"got a file list"<<std::endl;
		for(auto&& fname : filestlist){
			//std::cout<<"checking file "<<fname<<std::endl;
			std::string prefix="[N/A]";
			if(fname.length()>expectedprefix.length()) prefix = fname.substr(0,expectedprefix.length());
			bool ourfile = (prefix==expectedprefix);
			//std::cout<<"prefix is "<<prefix<<", ourfile is "<<ourfile<<std::endl;
			if( ourfile && files.count(fname)==0){
				std::string thisfname = dirname + "/" + fname;
				newfiles.emplace(thisfname,true);
				
				// note this file is now added
				files.emplace(fname,1);
			}
			
			// update num measurements
			num_meas = c_275_A->GetEntries();
			//std::cout<<"had "<<num_meas<<" measurements"<<std::endl;
		}
	}
	// add the new files, which should be sorted by name now.
	for(auto&& afile : newfiles){
		std::string thisfname = afile.first;
		// add the new file
		//std::cout<<"adding "<<thisfname<<" to tchains"<<std::endl;
		for(auto&& c : chains){
			c.second->Add(thisfname.c_str());
			c_dark->Add(thisfname.c_str());
		}
	}
	
	// get number of measurements
	num_meas = c_275_A->GetEntries();
	//std::cout<<"had "<<num_meas<<" measurements"<<std::endl;
	if(num_meas==0) return 0;
	
	// get num datapoints
	//std::cout<<"getting n datapoints"<<std::endl;
	c_dark->GetEntry(0);
	n_datapoints = values_dark.size();
	for(int i=0; i<n_datapoints; ++i){
		if(wavelengths.at(i)>272.8 && sample_273==0){
			sample_273 = i;
		} else if(wavelengths.at(i)>275.3){
			sample_276 = i;
			break;
		}
	}
	// fill wavevelength errors using bin widths
	wavelength_errors.resize(n_datapoints);
	std::fill(wavelength_errors.begin(), wavelength_errors.end(), wavelengths.at(1)-wavelengths.at(0));
	
	return num_meas;
}


int Plotter::SetBranchAddresses(){
	
	std::cout<<"setting led chain branch addresses"<<std::endl;
	
	for(auto&& c : chains){
		c.second->SetBranchAddress("value", &values_p);
		c.second->SetBranchAddress("error", &errors_p);
		c.second->SetBranchAddress("year",&yr);
		c.second->SetBranchAddress("month",&mon);
		c.second->SetBranchAddress("day",&dy);
		c.second->SetBranchAddress("hour",&hr);
		c.second->SetBranchAddress("min",&mn);
		c.second->SetBranchAddress("sec",&sc);
	}
	
	std::cout<<"setting dark chain branch addresses"<<std::endl;
	c_dark->SetBranchAddress("value", &values_dark_p);
	c_dark->SetBranchAddress("wavelength", &wavelengths_p);
	c_dark->SetBranchAddress("error", &errors_dark_p);
	c_dark->SetBranchAddress("year",&yr);
	c_dark->SetBranchAddress("month",&mon);
	c_dark->SetBranchAddress("day",&dy);
	c_dark->SetBranchAddress("hour",&hr);
	c_dark->SetBranchAddress("min",&mn);
	c_dark->SetBranchAddress("sec",&sc);
	
	chains.emplace("dark",c_dark);
	
	return 1;
}

bool Plotter::LoadConcentrations(std::string filename){
	
	std::cout<<"loading concentrations from file "<<filename<<std::endl;
	std::ifstream calib_data_file(filename.c_str());
	if(not calib_data_file.is_open()){
		std::cerr<<"couldn't open calibration data file "<<filename<<std::endl;
		return false;
	}
	
	calibration_data.clear();
	std::string line;
	//std::cout<<"getting lines"<<std::endl;
	while(getline(calib_data_file, line)){
		//std::cout<<"processing line "<<line<<std::endl;
		if(line.empty()) continue;
		if(line[0]=='#') continue;
		std::stringstream ss(line);
		std::string fname;
		double conc, concerr;
		//std::cout<<"parsing line"<<std::endl;
		if(!(ss >> fname >> conc >> concerr)){
			std::cerr<<"Failed to parse line "<<line<<"!"<<std::endl;
			break;
		}
		calibration_data.emplace(fname, std::pair<double,double>{conc,concerr});
		std::cout<<"file "<<fname<<" has concentration "<<conc<<std::endl;
	}
	calib_data_file.close();
	std::cout<<"loaded "<<calibration_data.size()<<" concentrations"<<std::endl;
	
	return true;
}

int Plotter::GetNextDarkEntry(std::string name, int ledon_entry_num, bool check){
	
	/*
	 Following code taken from the LoadOldFiles Tool:
	 In order to pull the correct dark trace for dark subtraction, rather than hard-coding the number of measurements,
	 pull the timestamps of the led-on measurement, then scan over the dark TTree until we find the dark entry
	 with the closest preceding timestamp.
	*/
	
	// get led timestamp
	TChain* c_led = chains.at(name);
	c_led->GetEntry(ledon_entry_num);
	struct tm ledtime;
	ledtime.tm_year = yr - 1900;
	ledtime.tm_mon = mon - 1;
	ledtime.tm_mday = dy;
	ledtime.tm_hour = hr;
	ledtime.tm_min = mn;
	ledtime.tm_sec = sc;
	time_t ledtime_t = mktime(&ledtime);
	
	// get the dark tree
	TChain* c_dark = chains.at("dark");
	// loop over dark entries
	int darkentry=-1;
	for(int i=0; i<c_dark->GetEntries(); ++i){
		c_dark->GetEntry(i);
		
		struct tm darktime;
		darktime.tm_year = yr - 1900;
		darktime.tm_mon = mon - 1;
		darktime.tm_mday = dy;
		darktime.tm_hour = hr;
		darktime.tm_min = mn;
		darktime.tm_sec = sc;
		time_t darktime_t = mktime(&darktime);
		
		// difftime does [ time_a - time_b ]
		double numsecs = difftime(ledtime_t, darktime_t);
		// do we use the closest in time dark, or the last dark before the led-on?
		// what if someone decides to do the dark measurements after the led?
		if(numsecs<0 && darkentry>=0) break;
		darkentry=i;
	}
	
	if(check){
		// compare against old method to check if it's good.
		int offset = offsets.at(name);
		int dark_entry = offset + ledon_entry_num*darksperloop;
		if(darkentry!=dark_entry){
			std::cerr<<"WARNING: closest dark entry by timestamp is "<<darkentry
				     <<", but offsets suggests entry "<<dark_entry<<std::endl;
			//assert(false);
			//exit(1);
			// not triggered: seems ok!
		}
	}
	
	return darkentry;
}

int Plotter::Initialise(std::string configfile){
	
	std::cout<<"making chains"<<std::endl;
	c_dark = new TChain("Dark");
	c_275_A = new TChain("275_A");
	c_275_B = new TChain("275_B");
	c_R_G_B = new TChain("R_G_B");
	c_White_385 = new TChain("White_385");
	
	chains = std::map<std::string, TChain*>{
		{"275_A", c_275_A},
		{"275_B", c_275_B},
		{"R_G_B", c_R_G_B},
		{"White_385", c_White_385}
	};
	
	SetBranchAddresses();

	bool ok = LoadConfig(configfile);
	
	return ok;
}

std::map<std::string,TMultiGraph*> Plotter::Execute(){
	
	// step 1: analyse calibration data for polynomial coffiecients
	std::cout<<"calling FitCalibData with do_calibration: "<<do_calibration<<std::endl;
	std::map<std::string, TMultiGraph*> plots_A = FitCalibrationData("275_A",do_calibration);
	std::map<std::string, TMultiGraph*> plots_B = FitCalibrationData("275_B",do_calibration);
	
	// combine plots
	//plots_A.insert(plots_B.begin(),plots_B.end());
	for(std::pair<const std::string, TMultiGraph*>& aplot : plots_A){
		std::string plotname = aplot.first;
		// edit the LED name to make the corresponding key in the plots_B map
		plotname.back()='B';
		std::cout<<"stripped plotname is "<<plotname<<std::endl;
		TMultiGraph* as_plots = aplot.second;
		int num_plots = as_plots->GetListOfGraphs()->GetEntries();
		TMultiGraph* bs_plots = plots_B.at(plotname);
		as_plots->Add(bs_plots);
		for(int i=0; i<num_plots; ++i){
			TGraphErrors* g = (TGraphErrors*)as_plots->GetListOfGraphs()->At(i);
			g->SetMarkerStyle(2);
		}
		for(int i=num_plots; i<(2*num_plots); ++i){
			TGraphErrors* g = (TGraphErrors*)as_plots->GetListOfGraphs()->At(i);
			g->SetMarkerStyle(5);
		}
	}
	
	/* actually, we can't parallelise these: they both use the same plotter object
	   so will all be sharing member variables, which isn't supported. yet.
	// XXX note that calling a function with async requires all arguments to be given -
	// even those with default values given in the declaration!
	std::vector<std::future<TMultiGraph*>> promises;
	promises.push_back(std::move(std::async(&Plotter::ExtractAbsorbance, this, "275_A","raw",true)));
	mg_all->Add(promises.front().get());
	*/
	
	return plots_A;
}

/////////////////////////////////

double PureFuncv2(double* x, double* par){
	// ok this function needs access to the appropriate pure graph for the correct LED
	// But to be invoked by a TF1, this function can only take an array of doubles
	// as arguments. So the last argument, par[8], is an LED index, which is used to
	// look up the appropriate pure reference TGraphErrors* from a global map.
	
	// par [0] = y-scaling
	// par [1] = x-scaling
	// par [2] = x-offset
	// par [3] = y-offset
	// par [4] = linear baseline offset
	// par [5] = shoulder gaussian scaling
	// par [6] = shoulder gaussian centre, restricted to > 282nm (RH shoulder)
	// par [7] = shoulder gaussian spread
	// par [8] = LED index for retrieving appropriate pure
	
	TGraphErrors* dark_subbed_pure = (TGraphErrors*)(dark_sub_pures.at(par[8]));
	
	double purepart = (par[0] * dark_subbed_pure->Eval((par[1]*x[0])+par[2]));
	double linpart = (par[4] * x[0]) + par[3];
	double shoulderpart = par[5]*exp(-0.5*TMath::Sq((x[0]-282.-abs(par[6]))/par[7]));
	double retval = purepart + linpart + shoulderpart;
	
	return retval;
}

TF1* Plotter::PureScaledPlusExtras(int led_num){
	
	// construct functional fit. We'll scale and add a linear background.
	// limit the pure function to a region in which we have light - no point fitting outside this region
	static int purever=0;
	purever++;
	std::string name="purev2_fct"+std::to_string(purever);
	const int wave_min = 260, wave_max = 300, numb_of_fitting_parameters = 9;
	TF1* pure = new TF1(name.c_str(), PureFuncv2, wave_min, wave_max, numb_of_fitting_parameters);
	// set default parameters
	pure->SetParameters(1.,1.,0.,0.,0.,0.,0.,10.);
	pure->FixParameter(8,led_num);
	return pure;
	
}

TF1* Plotter::GetAbsFunc(){
	// to really have a reasonable fit to data we in fact need *four* gaussians.
	TF1* twinpeaks = new TF1("twinpeaks"," [0]*exp(-0.5*TMath::Sq((x[0]-[1])/[2]))"
	                                     "+[3]*[0]*exp(-0.5*TMath::Sq((x[0]-[4])/[5]))"
	                                     "+[6]*[0]*exp(-0.5*TMath::Sq((x[0]-[7])/[8]))"
	                                     "+[0]*[9]*exp(-0.5*TMath::Sq((x[0]-[10])/[11]))",
	                                     270,277);
	
	// we may wish to make the amplitudes of all peaks relative to the amplitude of peak 1
	// by limiting the amplitudes 0-1, this ensures peak 1 is the largest, which helps
	// prevent misfits at low concentratiosn
	
	// par [0] = peak 1 (LH) amplitude
	// par [1] = peak 1 (LH) position
	// par [2] = peak 1 (LH) width
	
	// par [3] = amplitude of peak 1 RH shoulder (relative to peak 1)
	// par [4] = peak 1 RH shoulder position
	// par [5] = peak 1 RH shoulder width
	
	// par [6] = amplitude of peak 1 LH shoulder (relative to peak 1)
	// par [7] = peak 1 LH shoulder position
	// par [8] = peak 1 LH shoulder width
	
	// par [9]  = amplitude of peak 2
	// par [10] = peak 2 position
	// par [11] = peak 2 width
	twinpeaks->SetParName(0,"peak 1 amp");
	twinpeaks->SetParName(1,"peak 1 pos");
	twinpeaks->SetParName(2,"peak 1 wid");
	
	twinpeaks->SetParName(3,"RH shoulder amp");
	twinpeaks->SetParName(4,"RH shoulder pos");
	twinpeaks->SetParName(5,"RH shoulder wid");
	
	twinpeaks->SetParName(6,"LH shoulder amp");
	twinpeaks->SetParName(7,"LH shoulder pos");
	twinpeaks->SetParName(8,"LH shoulder wid");
	
	twinpeaks->SetParName(9,"peak 2 amp");
	twinpeaks->SetParName(10,"peak 2 pos");
	twinpeaks->SetParName(11,"peak 2 wid");
	
	// set initial positions of the absorption peaks
	twinpeaks->SetParameter(1,272.9);             // peak 1
	twinpeaks->SetParameter(4,274.);              // peak 1 RH shoulder
	twinpeaks->SetParameter(7,271.8);             // peak 1 LH shoulder
	twinpeaks->SetParameter(10,275.66);           // peak 2
	// restrict ranges to ensure each gaussian fits its appropriate peak
	twinpeaks->SetParLimits(1,272.5,273.5);       // peak 1
	twinpeaks->SetParLimits(4,274.,275.);         // peak 1 RH shoulder
	twinpeaks->SetParLimits(7,270,272.5);         // peak 1 LH shoulder
	twinpeaks->SetParLimits(10,275.,276.);        // peak 2
	
	// initial widths based on a good fit
	twinpeaks->SetParameter(2,0.5);
	twinpeaks->SetParameter(5,0.5);
	twinpeaks->SetParameter(8,0.6);
	twinpeaks->SetParameter(11,0.55);
	// restrict all widths to prevent crazy values
	twinpeaks->SetParLimits(2,0.3,1.);
	twinpeaks->SetParLimits(5,0.3,1.);
	twinpeaks->SetParLimits(8,0.5,1.);
	twinpeaks->SetParLimits(11,0.3,1.);
	
	// initial amplitudes from a good fit - amplitudes of peaks 1 and 2
	// will be overrriden based on results from one of the simpler methods
	twinpeaks->SetParameter(0,0.46);
	twinpeaks->SetParameter(3,0.5);    // absolute amplitude ~0.23
	twinpeaks->SetParameter(6,0.2);    // absolute amplitude ~0.09
	twinpeaks->SetParameter(9,0.6);    // absolute amplitude ~0.27
	// set limits. peak 1 is 0-1 by definition of absorption
	// others are 0-1 since they're all smaller than peak 1
	twinpeaks->SetParLimits(0,0.,1.);  // peak 1
	twinpeaks->SetParLimits(3,0.,1.);  // peak 1
	twinpeaks->SetParLimits(6,0.,1.);  // peak 1
	twinpeaks->SetParLimits(9,0.,1.);  // peak 1
	
	
	return twinpeaks;
}

std::pair<double,double> Plotter::CalculateError(TF1* abs_func, double peak1_pos, double peak2_pos){
	
	// ok so we have multiple (4) contributions to the function at each peak.
	// to estimate the error on the height of each peak, i'm gonna add in quadrature
	// the error from each contribution to that peak. The error from each contribution is
	// the error on that gaussian's amplitude, but scaled by it's contribution
	// to the function at that peak's position. e.g. if gaussian 4
	// (the peak2 RH shoulder) does not contribute to the amplitude at
	// peak 1, then its error would be scaled to 0. The contribution from the gaussian
	// centered on peak 1 would be scaled to 1, whereas the contribution to error
	// from the shoulders (peaks 2 and 3) would be something between.
	
	double totalerror1=0;
	double totalerror2=0;
	
	// we have 4 gaussians
	for(int i=0; i<4; ++i){
		
		// error on amplitude of this gaus
		double error_centre = abs_func->GetParError(i*3);
		// note that amplitudes of gaussians 1-3 are relative to that of gaussian 0
		// so an error of 0.5 on gaussian 1 is actually only an error of 0.05
		// in the case that gaussian 0 peak height is 0.1
		if(i!=0) error_centre *= abs_func->GetParameter(0);
		
		// get the relative contribution to the function from this gaussian
		// first get the height at it's centre...
		double amp_centre;
		if(i==0){
			amp_centre = abs_func->GetParameter(0);
		} else {
			amp_centre = abs_func->GetParameter(i*3)*abs_func->GetParameter(0);
		}
		// now make a temporary gaus to obtain its height at peak 1
		TF1 agaus("agaus","[0]*exp(-0.5*TMath::Sq((x-[1])/[2]))",270., 300.);
		// copy over its parameters...
		for(int j=0; j<3; ++j){
			if(j==0 && i!=0){
				agaus.SetParameter(j,abs_func->GetParameter(i*3+j)*abs_func->GetParameter(0));
			} else {
				agaus.SetParameter(j,abs_func->GetParameter(j+i*3));
			}
		}
		// evaluate at peak 1
		double amp_here = agaus.Eval(peak1_pos);
		double relative_amp = amp_here/amp_centre;
		// scale the error (on the ampltiude at the center) by the relative height
		// at peak 1, compared to at the centre
		double error_here = error_centre*relative_amp;
		
		/*
		std::cout<<"complex peak fit gaus "<<i<<" has ampltiude "
		         <<agaus.GetParameter(0)<<"+-"<<error_centre<<" scaled by "
		         <<relative_amp<<" to "<<error_here<<" at peak 1";
		*/
		
		totalerror1+= TMath::Sq(error_here);
		
		// repeat for contribution to peak2
		amp_here = agaus.Eval(peak2_pos);
		relative_amp = amp_here/amp_centre;
		// scale the error on the amplitude at the centre
		error_here = error_centre*relative_amp;
		
		//std::cout<<" and scaled by "<<relative_amp<<" to "<<error_here<<" at peak 2"<<std::endl;
		
		totalerror2+= TMath::Sq(error_here);
	}
	
	// take sqrt of totalerror1 and totalerror2 each to get quadrature sum of contributions...
	// but then take square to add the errors on the two peak amplitudes in quadrature
	//std::cout<<"total error at peak 1 "<<sqrt(totalerror1)<<" and at peak 2 "<<sqrt(totalerror2)<<std::endl;
	
	// FIXME this is WRONG because it doesn't account for correlation between the components
	// ALSO when fixing it: for compatibility with other formulas, we should return two values
	// that add in quadrature to give the total error on the difference. :/ (or generalise the calling code).
	
	return std::pair<double,double>{sqrt(totalerror1),sqrt(totalerror2)};
}

std::map<std::string, TMultiGraph*> Plotter::FitCalibrationData(std::string name, bool make_cal_curve){
	
	// loop over the data, do dark subtraction,
	// fit the dark-subtracted data, extract the absorption graph,
	// fit the absorption graph and pull the peak heights and difference.
	
	// depending on whether we are performing a calibration or making a stability plot,
	// either use the peak height difference to look up concentration,
	// or look up the concentration and add a point to the calibration curve
	
	// keep vectors of the key parameters for making plots
	// absorption peak height diffs
	std::vector<double> raw_peakdiffs, raw_peakdiff_errs;
	std::vector<double> simple_peakdiffs, simple_peakdiff_errs;
	std::vector<double> complex_peakdiffs, complex_peakdiff_errs;
	// true concentrations if making a calibration curve
	std::vector<double> true_concentrations, true_concentration_errs;
	// calculated concentrations if using pre-existing calibration curve
	std::vector<double> raw_concentrations, raw_concentration_errs;
	std::vector<double> simple_concentrations, simple_concentration_errs;
	std::vector<double> complex_concentrations, complex_concentration_errs;
	
	//std::cout<<"getting pure"<<std::endl;
	int pure_index=-1;
	if(dark_sub_pures_byname.count(name)==0){
		// do not currently have this pure reference trace in memory;
		// get it from file, making the file from the first entry
		// in the chains if it doesn't exist.
		TGraphErrors* dark_subbed_pure = MakePure(name);
		if(dark_subbed_pure==nullptr){
			return std::map<std::string,TMultiGraph*>{};
		}
		pure_index = dark_sub_pures.size();
		dark_sub_pures_byname.emplace(std::pair<std::string,int>{name,pure_index});
		dark_sub_pures.push_back(dark_subbed_pure);
	} else {
		// already have this pure reference in memory
		pure_index = dark_sub_pures_byname.at(name);
	}
	// construct the pure water function based on the pure water trace.
	TF1* pureplusextras = PureScaledPlusExtras(pure_index);
	pureplusextras->SetLineColor(kBlack);
	pureplusextras->SetLineWidth(1);
	
	//std::cout<<"getting chains"<<std::endl;
	TChain* c_led = chains.at(name);
	TChain* c_dark = chains.at("dark");
	int n_entries = c_led->GetEntries(); //53
	
	std::string title = std::string("g_") + name;
	
	// functions that map peak height diff to concentration
	// we will either need these as we go (for making concentration stability plots)
	// or we'll make them by fitting the data at the end.
	TF1* calib_func_raw = new TF1("calib_func_raw", "pol6", 0, 0.25);
	// line properties have to be set before calling Fit to apply to the saved function
	calib_func_raw->SetLineWidth(1);
	calib_func_raw->SetLineColor(kSpring-5);
	TF1* calib_func_simple = new TF1("calib_func_simple", "pol6", 0, 0.25);
	calib_func_simple->SetLineColor(kRed);
	calib_func_simple->SetLineWidth(1);
	TF1* calib_func_complex = new TF1("calib_func_complex", "pol6", 0, 0.25);
	calib_func_complex->SetLineWidth(1);
	calib_func_complex->SetLineColor(kBlue);
	
	if(!make_cal_curve){
		// If we're making a stability plot of concentration look up
		// the calibration curve parameters based on a past calibration
		std::vector<double> calib_coefficients_raw = GetCalibCurve(name, "raw");
		calib_func_raw->SetParameters(calib_coefficients_raw.data());
		std::vector<double> calib_coefficients_simple = GetCalibCurve(name, "simple");
		calib_func_simple->SetParameters(calib_coefficients_simple.data());
		std::vector<double> calib_coefficients_complex = GetCalibCurve(name, "complex");
		calib_func_complex->SetParameters(calib_coefficients_complex.data());
	}
	
	// a temporary canvas for showing the fits as we go.
	TCanvas* c_temp = new TCanvas("c_temp","c_temp",1280,800);
	
	int n_measurements=0;
	double lastconc=-99;
	
	if(max_entries>0) n_entries = std::min(n_entries,first_entry+max_entries);
	std::cout<<"looping over entries "<<first_entry<<" - "<<n_entries<<std::endl;
	for(int i=first_entry; i<n_entries; i++){
		
		std::cout<<"measurement "<<i<<": ";
		
		// get LED on trace
		int nbytes = c_led->GetEntry(i);
		if(nbytes<=0){
			std::cerr<<"COULDN'T LOAD ENTRY "<<i<<" FROM LED CHAIN"<<std::endl;
			exit(-1);
		}
		
		// get nearest dark trace
		int dark_entry = GetNextDarkEntry(name, i, false);
		c_dark->GetEntry(dark_entry);
		
		// subtract dark
		for(int j=0; j<n_datapoints; ++j){
			values.at(j) -= values_dark.at(j);
		}
		
		// split into in- and -out of absorbance bands
		std::vector<double> inband_values;
		std::vector<double> inband_wls;
		std::vector<double> inband_errors;
		std::vector<double> sideband_values;
		std::vector<double> sideband_wls;
		std::vector<double> other_values;
		std::vector<double> other_wls;
		sample_273=0;
		sample_276=0;
		for(int j=0; j<n_datapoints; ++j){
			if(wavelengths.at(j)>270 && wavelengths.at(j)<280){
				// in band
				inband_values.push_back(values.at(j));
				inband_wls.push_back(wavelengths.at(j));
				inband_errors.push_back(errors.at(j));
			} else if(wavelengths.at(j)>260 && wavelengths.at(j)<300){
				// side band (lobes)
				sideband_values.push_back(values.at(j));
				sideband_wls.push_back(wavelengths.at(j));
			} else {
				other_values.push_back(values.at(j));
				other_wls.push_back(wavelengths.at(j));
			}
			
			// make a note of the index of the 273 and 276 peaks
			// within the raw in-band absorption graph arrays
			if(wavelengths.at(j)>272.8 && sample_273==0){
				sample_273 = inband_values.size()-1;
			} else if(wavelengths.at(j)>275.3 && sample_276==0){
				sample_276 = inband_values.size()-1;
			}
		}
		
		// make TGraphErrors from data
		TGraph* g_inband = new TGraph(inband_values.size(), inband_wls.data(), inband_values.data());
		TGraph* g_sideband = new TGraph(sideband_values.size(), sideband_wls.data(), sideband_values.data());
		TGraph* g_other = new TGraph(other_values.size(), other_wls.data(), other_values.data());
		
		std::string thistitle = title+"_"+std::to_string(i);
		std::string intitle = thistitle+"_in";
		std::string sidetitle = thistitle+"_side";
		std::string othertitle = thistitle+"_other";
		g_inband->SetTitle(thistitle.c_str());
		g_inband->SetName(thistitle.c_str());
		g_sideband->SetTitle(sidetitle.c_str());
		g_sideband->SetName(sidetitle.c_str());
		g_other->SetTitle(othertitle.c_str());
		g_other->SetName(othertitle.c_str());
		g_inband->SetMarkerColor(kRed);
		g_inband->SetMarkerStyle(2);
		g_sideband->SetMarkerColor(kBlue);
		g_sideband->SetMarkerStyle(5);
		
		// fit with pure scaled
		std::cout<<"fitting pure"<<std::endl;
		g_sideband->Fit(pureplusextras,"RNQ");
		// save an image to check the fit
		std::string purename = "images/purefit_"+name+"_"+std::to_string(i)+".png";
		/*if(!draw)*/ gROOT->SetBatch(true);
		c_temp->cd();
		g_sideband->Draw("AP");
		g_sideband->GetYaxis()->SetRangeUser(0,1.1*(*std::max_element(g_inband->GetY(),g_inband->GetY()+g_inband->GetN())));
		g_inband->Draw("same P");
		pureplusextras->Draw("same");
		
		c_temp->SaveAs(purename.c_str());
		gROOT->SetBatch(false);

		
		// calculate absorbance from ratio of fit to data in absorption region
		std::cout<<"generating absorption plot"<<std::endl;
		TGraphErrors* g_abs = new TGraphErrors(inband_values.size());
		std::string restitle = "g_abs_"+std::to_string(i);
		for(int k=0; k<inband_values.size(); ++k){
			double wlval, dataval;
			g_inband->GetPoint(k, wlval, dataval);
			double fitval = pureplusextras->Eval(wlval);
			g_abs->SetPoint(k, wlval, log10(fitval/dataval));
			//std::cout<<"setting absorption point "<<k<<" to log10("<<fitval<<"/"<<dataval
			//         <<") = "<<log10(fitval/dataval)<<std::endl;
			
			// error on ratio is sqrt((Δa/a)²+(Δb/b)²)
			// https://www.statisticshowto.com/error-propagation/
			// error on logY(X) is (ΔX/X)*(1/ln(Y))
			// https://physics.stackexchange.com/questions/95254/the-error-of-the-natural-logarithm/95278#95278
			// at least in the regime where ΔX<<X. if outside this regime the non-linear
			// nature of logs means the errors will be asymmetric. In this case one can use
			// as a rough guide ΔY = log(X-ΔX) -> log(X+ΔX).
			double error_on_data = inband_errors.at(k);
			// uhhh what's the error on the extrapolated fit value - i.e. transmitted intensity.
			// no idea. Let's assume it's of the same order as the data??? FIXME
			double error_on_fitval = inband_errors.at(k);
			double err_on_ratio = sqrt(TMath::Sq(error_on_data/dataval)
			                          +TMath::Sq(error_on_fitval/fitval));
			double err_on_ratio_over_ratio = err_on_ratio / (fitval/dataval);
			/*
			std::cout<<"measured intensity: "<<dataval<<"+- "<<error_on_data<<"\n"
			         <<", fractional error: "<<(error_on_data/dataval)<<"\n"
			         <<", transmitted intensity: "<<fitval<<", +- "<<error_on_fitval<<"\n"
			         <<", fractional error: "<<(error_on_fitval/fitval)<<"\n"
			         <<", total error on ratio: "<<err_on_ratio
			         <<", ratio of measured / transmitted: "<<(fitval/dataval)
			         <<", ratio of error to value: "<<err_on_ratio_over_ratio
			         <<std::endl;
			*/
			g_abs->SetPointError(k, wavelength_errors.at(0)/2., err_on_ratio_over_ratio*(1./log(10.)));
		}
		
		// Extract the height of the peaks.
		// as the simplest estimate of peak heights, just take the graph value at the peaks.
		std::pair<double,double> peak_posns_raw;
		std::pair<double,double> peak_heights_raw;
		std::pair<double,double> peak_errs_raw;
		std::string rawname = "images/raw_"+name+"_"+std::to_string(i)+".png";
		std::cout<<"fitting raw"<<std::endl;
		ok = DoRawFit(g_abs, peak_posns_raw, peak_heights_raw, peak_errs_raw, rawname);
		
		// A better estimate fits the peaks with gaussians to find a better estimate of maximum.
		std::pair<double,double> peak_posns_simple;
		std::pair<double,double> peak_heights_simple;
		std::pair<double,double> peak_errs_simple;
		std::string simplename = "images/simple_"+name+"_"+std::to_string(i)+".png";
		std::cout<<"fitting simple"<<std::endl;
		ok = DoSimpleFit(g_abs, peak_posns_simple, peak_heights_simple, peak_errs_simple, simplename);
		
		// as the most complex method, we fit the whole absorbance region with 4 overlapping gaussians
		std::pair<double,double> peak_posns_complex;
		std::pair<double,double> peak_heights_complex;
		std::pair<double,double> peak_errs_complex;
		std::string complexname = "images/complex_"+name+"_"+std::to_string(i)+".png";
		std::cout<<"fitting complex"<<std::endl;
		ok = DoComplexFit(g_abs, peak_posns_complex, peak_heights_complex, peak_errs_complex, peak_heights_raw, complexname);
		
		double raw_peakdiff = peak_heights_raw.first - peak_heights_raw.second;
		double raw_peakdiff_err = sqrt(TMath::Sq(peak_errs_raw.first)+TMath::Sq(peak_errs_raw.second));
		
		double simple_peakdiff = peak_heights_simple.first - peak_heights_simple.second;
		double simple_peakdiff_err = sqrt(TMath::Sq(peak_errs_simple.first)+TMath::Sq(peak_errs_simple.second));
		
		double complex_peakdiff = peak_heights_complex.first - peak_heights_complex.second;
		double complex_peakdiff_err = sqrt(TMath::Sq(peak_errs_complex.first)+TMath::Sq(peak_errs_complex.second));
		
		
		std::cout<<"raw diff: "<<peak_heights_raw.first<<" - "<<peak_heights_raw.second
		         <<" = "<<raw_peakdiff<<std::endl;
		std::cout<<"simple diff: "<<peak_heights_simple.first<<" - "<<peak_heights_simple.second
		         <<" = "<<simple_peakdiff<<std::endl;
		std::cout<<"complex diff: "<<peak_heights_complex.first<<" - "<<peak_heights_complex.second
		         <<" = "<<complex_peakdiff<<std::endl;
		
		
		// one last thing: how do we handle bad fit values (negative peak height differences?)
		// XXX coerce or skip, that is the question....
		raw_peakdiff=std::max(0.,raw_peakdiff);
		simple_peakdiff=std::max(0.,simple_peakdiff);
		complex_peakdiff=std::max(0.,complex_peakdiff);
		++n_measurements; // in case we choose to skip.
		
		// add them to the vectors
		raw_peakdiffs.push_back(raw_peakdiff);
		simple_peakdiffs.push_back(simple_peakdiff);
		complex_peakdiffs.push_back(complex_peakdiff);
		
		// add the errors FIXME for now since errors aren't correct, just use nominal values...
		raw_peakdiff_errs.push_back(0.03);      //raw_peakdiff_err);
		simple_peakdiff_errs.push_back(0.03);   //simple_peakdiff_err);
		complex_peakdiff_errs.push_back(0.03);  //complex_peakdiff_err);
		
		// Map difference in peak height to concentration.
		if(make_cal_curve){
			// if we're making a calibration curve then the true concentration is known,
			// so look it up from the reference file.
			
			// first get name of the file this measurement is from
			std::string current_file = c_led->GetTree()->GetCurrentFile()->GetName();
			//std::cout<<"filename: "<<current_file<<std::endl;
			
			// strip off path
			if(current_file.find("/")!=std::string::npos){
				current_file = current_file.substr(current_file.find_last_of("/")+1,std::string::npos);
			}
			
			// sanity check that this filename is in the concentration reference map
			if(calibration_data.count(current_file)==0){
				std::cerr<<"Couldn't find file "<<current_file<<" in calibration data map!"<<std::endl;
				for(auto&& am : calibration_data){
					std::cerr<<am.first<<", ";
				}
				std::cerr<<std::endl;
				return std::map<std::string,TMultiGraph*>{};
			}
			
			// get the concentration and error
			double current_conc = calibration_data.at(current_file).first;
			std::cout<<"concentration for measurement "<<i<<" is "<<current_conc<<std::endl;
			double current_concerr = calibration_data.at(current_file).second;
			//if(lastconc==current_conc) break; // HACK - stop the loop once the concentration stops changing
			// XXX double comparison; this isn't reliable and can result in breaking early.
			lastconc=current_conc;
			
			// add them to the vectors
			true_concentrations.push_back(current_conc);
			true_concentration_errs.push_back(current_concerr);
			
		} else {
			
			// if we're making a stability plot we already have a calibration curve available
			// so estimate the concentration based on the peak height difference
			
			// each fitting method has its own calibration function, which we need to
			// solve for concentration (x) from absorbance (y), with 0.01 < x < 0.21
			//std::cout<<"converting to concentration"<<std::endl;
			std::cout<<"calculating raw concentration"<<std::endl;
			double raw_conc = calib_func_raw->GetX(raw_peakdiff, 0.001, 0.22);			
			std::cout<<"calculating simple concentration"<<std::endl;
			double simple_conc = calib_func_simple->GetX(simple_peakdiff, 0.001, 0.22);
			std::cout<<"calculating complex concentration"<<std::endl;
			double complex_conc = calib_func_complex->GetX(complex_peakdiff, 0.001, 0.22);
			raw_concentrations.push_back(raw_conc);
			simple_concentrations.push_back(simple_conc);
			complex_concentrations.push_back(complex_conc);
			
			std::cout<<"raw conc for peak diff "<<raw_peakdiff<<" = "<<raw_conc<<std::endl;
			std::cout<<"simple conc for peak diff "<<simple_peakdiff<<" = "<<simple_conc<<std::endl;
			std::cout<<"complex conc for peak diff "<<complex_peakdiff<<" = "<<complex_conc<<std::endl;
			
			// error on concentration is error on height times gradient at that point
			double raw_conc_err = raw_peakdiff_err * calib_func_raw->Derivative(raw_peakdiff);
			double simple_conc_err = simple_peakdiff_err * calib_func_simple->Derivative(simple_peakdiff);
			double complex_conc_err = complex_peakdiff_err * calib_func_complex->Derivative(complex_peakdiff);
			raw_concentration_errs.push_back(raw_conc_err);
			simple_concentration_errs.push_back(simple_conc_err);
			complex_concentration_errs.push_back(complex_conc_err);
			
		}
		
		delete g_inband;
		delete g_sideband;
		delete g_other;
		delete g_abs;
		
	}
	
	// general container for returned results.
	std::map<std::string,TMultiGraph*> output_plots;
	
	// make a generic numberline for plots against measurement number
	std::cout<<"n_measurements is "<<n_measurements<<", true concs size is "<<true_concentrations.size()<<std::endl;
	std::vector<double> numberline(n_measurements);
	std::iota(numberline.begin(), numberline.end(), 0);
	std::vector<double> zeros(n_measurements); // errors on measurement number: 0
	std::fill(zeros.begin(), zeros.end(), 0);
	
	// plot peak height difference vs measurement
	// generally useful whether we're doing calibration or stability monitoring
	TGraphErrors* g_peakdiff_raw = new TGraphErrors(numberline.size(), numberline.data(), raw_peakdiffs.data(), zeros.data(), raw_peakdiff_errs.data());
	title = "g_peakdiff_raw_"+name;
	g_peakdiff_raw->SetName(title.c_str());
	g_peakdiff_raw->SetTitle(title.c_str());
	g_peakdiff_raw->SetLineColor(kSpring-5);
	g_peakdiff_raw->SetMarkerColor(kSpring-5);
	
	TGraphErrors* g_peakdiff_simple = new TGraphErrors(numberline.size(), numberline.data(), simple_peakdiffs.data(), zeros.data(), simple_peakdiff_errs.data());
	title = "g_peakdiff_simple_"+name;
	g_peakdiff_simple->SetName(title.c_str());
	g_peakdiff_simple->SetTitle(title.c_str());
	g_peakdiff_simple->SetLineColor(kRed);
	g_peakdiff_simple->SetMarkerColor(kRed);
	
	TGraphErrors* g_peakdiff_complex = new TGraphErrors(numberline.size(), numberline.data(), complex_peakdiffs.data(), zeros.data(), complex_peakdiff_errs.data());
	title = "g_peakdiff_complex_"+name;
	g_peakdiff_complex->SetName(title.c_str());
	g_peakdiff_complex->SetTitle(title.c_str());
	g_peakdiff_complex->SetLineColor(kBlue);
	g_peakdiff_complex->SetMarkerColor(kBlue);
	
	// save results
	title = "mg_peakdiffs_"+name;
	TMultiGraph* mg_peakdiffs = new TMultiGraph(title.c_str(),title.c_str());
	mg_peakdiffs->Add(g_peakdiff_raw);
	mg_peakdiffs->Add(g_peakdiff_simple);
	mg_peakdiffs->Add(g_peakdiff_complex);
	output_plots.emplace(std::pair<std::string,TMultiGraph*>{title,mg_peakdiffs});
	title = title+".root";
	mg_peakdiffs->SaveAs(title.c_str());
	
	if(make_cal_curve){
		
		// Make the calibration curves
		// ===========================
		
		// 1. Raw values method
		TGraphErrors* calib_curve_raw = new TGraphErrors(n_measurements, true_concentrations.data(), raw_peakdiffs.data(), true_concentration_errs.data(), raw_peakdiff_errs.data());
		
		title+"g_calib_raw_"+name;
		calib_curve_raw->SetName(title.c_str());
		calib_curve_raw->SetTitle(title.c_str());
		calib_curve_raw->SetMarkerColor(kSpring-5);
		calib_curve_raw->SetLineColor(kSpring-5);
		calib_curve_raw->SetLineWidth(0);
		if(name=="275_A"){
			calib_curve_raw->SetMarkerStyle(20);
		} else if(name=="275_B"){
			calib_curve_raw->SetMarkerStyle(30);
		}
		
		// 2. Simple 2-gaussian method
		TGraphErrors* calib_curve_simple = new TGraphErrors(n_measurements, true_concentrations.data(), simple_peakdiffs.data(), true_concentration_errs.data(), simple_peakdiff_errs.data());
		title = "g_calib_simple_"+name;
		calib_curve_simple->SetName(title.c_str());
		calib_curve_simple->SetTitle(title.c_str());
		calib_curve_simple->SetMarkerColor(kRed);
		calib_curve_simple->SetLineColor(kRed);
		calib_curve_simple->SetLineWidth(0);
		if(name=="275_A"){
			calib_curve_simple->SetMarkerStyle(20);
		} else if(name=="275_B"){
			calib_curve_simple->SetMarkerStyle(30); // 30=☆, 34=+
		}
		
		// 3. complex 4-gaussian method
		TGraphErrors* calib_curve_complex = new TGraphErrors(n_measurements, true_concentrations.data(), complex_peakdiffs.data(), true_concentration_errs.data(), complex_peakdiff_errs.data());
		title = "g_calib_complex_"+name;
		calib_curve_complex->SetName(title.c_str());
		calib_curve_complex->SetTitle(title.c_str());
		calib_curve_complex->SetMarkerColor(kBlue);
		calib_curve_complex->SetLineColor(kBlue);
		calib_curve_complex->SetLineWidth(0);
		if(name=="275_A"){
			calib_curve_complex->SetMarkerStyle(20);
		} else if(name=="275_B"){
			calib_curve_complex->SetMarkerStyle(30);
		}
		
		// Fit the calibration curves
		// ==========================
		// 1. Raw Values method
		
		calib_curve_raw->Fit("calib_func_raw","Q","",0,lastconc); // "N"
		
		// save fit parameters in a simple text file
		std::string rawcoeffs_filename = std::string("calib_coeffs_")+name+"_raw.txt";
		// open file for writing, discard any existing content
		std::ofstream rawcoeffs_file(rawcoeffs_filename, std::ofstream::out | std::ofstream::trunc);
		if(!rawcoeffs_file.is_open()){
			std::cerr<<"Failed to open file "<<rawcoeffs_filename
				     <<" for writing fit parameters"<<std::endl;
		}
		std::cout<<"raw fit curve pars for led "<<name<<" are {";
		for(int i=0; i<calib_func_raw->GetNpar(); ++i){
			if(i>0) std::cout<<", ";
			std::cout<<calib_func_raw->GetParameter(i);
			if(rawcoeffs_file.is_open()) rawcoeffs_file << calib_func_raw->GetParameter(i) << "\n";
		}
		std::cout<<"}"<<std::endl;
		rawcoeffs_file.close();
		
		// 2. Simple Fit method
		calib_curve_simple->Fit("calib_func_simple","Q","",0,lastconc); // "N"
		
		std::string simplecoeffs_filename = std::string("calib_coeffs_")+name+"_simple.txt";
		// open file for writing, discard any existing content
		std::ofstream simplecoeffs_file(simplecoeffs_filename, std::ofstream::out | std::ofstream::trunc);
		if(!simplecoeffs_file.is_open()){
			std::cerr<<"Failed to open file "<<simplecoeffs_filename
				     <<" for writing fit parameters"<<std::endl;
		}
		std::cout<<"simple fit curve pars for led "<<name<<" are {";
		for(int i=0; i<calib_func_simple->GetNpar(); ++i){
			if(i>0) std::cout<<", ";
			std::cout<<calib_func_simple->GetParameter(i);
			if(simplecoeffs_file.is_open()) simplecoeffs_file << calib_func_simple->GetParameter(i) << "\n";
		}
		std::cout<<"}"<<std::endl;
		simplecoeffs_file.close();
		
		// 3. complex fit
		calib_curve_complex->Fit("calib_func_complex","Q","",0,lastconc); // "N"
		
		std::string complexcoeffs_filename = std::string("calib_coeffs_")+name+"_complex.txt";
		// open file for writing, discard any existing content
		std::ofstream complexcoeffs_file(complexcoeffs_filename, std::ofstream::out | std::ofstream::trunc);
		if(!complexcoeffs_file.is_open()){
			std::cerr<<"Failed to open file "<<complexcoeffs_filename
				     <<" for writing fit parameters"<<std::endl;
		}
		std::cout<<"complex fit curve pars for led "<<name<<" are {";
		for(int i=0; i<calib_func_complex->GetNpar(); ++i){
			if(i>0) std::cout<<", ";
			std::cout<<calib_func_complex->GetParameter(i);
			if(complexcoeffs_file.is_open()) complexcoeffs_file << calib_func_complex->GetParameter(i) << "\n";
		}
		std::cout<<"}"<<std::endl;
		complexcoeffs_file.close();
		
		/*
		TCanvas cc("cc","cc",1024,800);
		calib_curve_simple->Draw("AP");
		calib_func_simple->Draw("same");
		calib_curve_complex->Draw("same P");
		calib_func_complex->Draw("same");
		calib_curve_raw->Draw("same P");
		calib_func_raw->Draw("same");
		cc.Modified();
		cc.Update();
		*/
		
		// save results
		// ============
		title = "mg_calibs_"+name;
		TMultiGraph* mg_calibs = new TMultiGraph(title.c_str(), title.c_str());
		mg_calibs->Add(calib_curve_simple);
		mg_calibs->Add(calib_curve_complex);
		mg_calibs->Add(calib_curve_raw);
		output_plots.emplace(std::pair<std::string,TMultiGraph*>{title,mg_calibs});
		title = title+".root";
		mg_calibs->SaveAs(title.c_str());
		
		// make a plot of true concentration vs measurement number for reference.
		TGraphErrors* g_true_concentrations = new TGraphErrors(numberline.size(), numberline.data(), true_concentrations.data(), zeros.data(), true_concentration_errs.data());
		title = "g_true_concs_"+name;
		g_true_concentrations->SetName(title.c_str());
		g_true_concentrations->SetTitle(title.c_str());
		title = std::string("m")+title;
		TMultiGraph* mg_true_concs = new TMultiGraph(title.c_str(),title.c_str());
		mg_true_concs->Add(g_true_concentrations);
		output_plots.emplace(std::pair<std::string,TMultiGraph*>{title,mg_true_concs});
		title = title+".root";
		mg_true_concs->SaveAs(title.c_str());
		
	} else {
		
		// make stability with derived concentration if not doing calibration
		TGraphErrors* g_conc_raw = new TGraphErrors(numberline.size(), numberline.data(), raw_concentrations.data(), zeros.data(), raw_concentration_errs.data());
		title = "g_conc_raw_"+name;
		g_conc_raw->SetName(title.c_str());
		g_conc_raw->SetTitle(title.c_str());
		g_conc_raw->SetLineColor(kSpring-5);
		g_conc_raw->SetMarkerColor(kSpring-5);
		
		TGraphErrors* g_conc_simple = new TGraphErrors(numberline.size(), numberline.data(), simple_concentrations.data(), zeros.data(), simple_concentration_errs.data());
		title="g_conc_simple_"+name;
		g_conc_simple->SetName(title.c_str());
		g_conc_simple->SetTitle(title.c_str());
		g_conc_simple->SetLineColor(kRed);
		g_conc_simple->SetMarkerColor(kRed);
		
		TGraphErrors* g_conc_complex = new TGraphErrors(numberline.size(), numberline.data(), complex_concentrations.data(), zeros.data(), complex_concentration_errs.data());
		title="g_conc_complex_"+name;
		g_conc_complex->SetName(title.c_str());
		g_conc_complex->SetTitle(title.c_str());
		g_conc_complex->SetLineColor(kBlue);
		g_conc_complex->SetMarkerColor(kBlue);
		
		// save results
		title = "mg_concs_"+name;
		TMultiGraph* mg_concs = new TMultiGraph(title.c_str(),title.c_str());
		mg_concs->Add(g_conc_raw);
		mg_concs->Add(g_conc_simple);
		mg_concs->Add(g_conc_complex);
		output_plots.emplace(std::pair<std::string,TMultiGraph*>{title,mg_concs});
		title=title+".root";
		mg_concs->SaveAs(title.c_str());
		
	}
	
	/*
	gSystem->ProcessEvents();
	while(gROOT->FindObject("c_temp")!=nullptr){
		gSystem->ProcessEvents();
		std::this_thread::sleep_for(std::chrono::milliseconds(100));
	}
	*/
	
	delete c_temp;
	delete calib_func_raw;
	delete calib_func_simple;
	delete calib_func_complex;
	
	//std::cout<<"returning multigraph"<<std::endl;
	return output_plots;
	
}

void Plotter::SetTimeAxis(TH1* hist, long t0_seconds){
	// start_time is in seconds: plotted times should be in seconds relative to this
	gStyle->SetTimeOffset(t0_seconds);
	// then you set the x-axis to time format
	hist->GetXaxis()->SetTimeDisplay(1);
	hist->GetXaxis()->SetTimeFormat("%m/%d/%y %H:%M");
	hist->GetXaxis()->SetLabelSize(0.03);
}

bool Plotter::CheckPath(std::string path, std::string& type){
	struct stat s;
	if(stat(path.c_str(),&s)==0){
		if(s.st_mode & S_IFDIR){        // mask to extract if it's a directory?? how does this work?
			type="d";  //it's a directory
			return true;
		} else if(s.st_mode & S_IFREG){ // mask to check if it's a file??
			type="f"; //it's a file
			return true;
		} else {
			// exists, but neither file nor directory?
			type="???";
			return false;
			//assert(false&&"Check input path: stat says it's neither file nor directory..?");
		}
	} else {
		// does not exist - could be a pattern, e.g. "/path/to/rootfiles_*.root"
		type="none";
		return false;
	}
	return false;
}

TGraphErrors* Plotter::MakePure(std::string name, bool overwrite){
	
	std::string purefile = std::string("../GDConcMeasure/pureDarkSubtracted_") + name + "_" + purerefset + ".root";
	if(!overwrite){
		// check if the pure reference file already exists
		std::string type="";
		bool exists = CheckPath(purefile, type);
		if(exists && type=="f"){
			std::cout<<"Pure file "<<purefile<<" already exists, using it"<<std::endl;
			TFile* _file0 = TFile::Open(purefile.c_str());
			TGraphErrors* dark_subtracted_pure = (TGraphErrors*)_file0->Get("Graph");
			return dark_subtracted_pure;
		}
		else if(exists && type!="f"){
			std::cerr<<"Pure file "<<purefile<<" exists but is not a standard file. Please investigate"<<std::endl;
			return nullptr;
		} else if(!exists && !generate_pure){
			std::cerr<<"Pure file "<<purefile<<" does not exist"<<std::endl;
			return nullptr;
		}
	}
	
	std::cout<<"Pure file "<<purefile<<" does not exist or overwriting it"<<std::endl;
	
	TChain* c_led = chains.at(name);
	TChain* c_dark = chains.at("dark");
	
	// get LED on trace
	int nbytes = c_led->GetEntry(0);
	if(nbytes<=0){
		std::cerr<<"COULDN'T LOAD ENTRY 0 FROM LED CHAIN"<<std::endl;
		return nullptr;
	}
	
	/*
	// sanity check: look up the concentration of this measurement.
	// should be 0 if it's entry 0 in the chain, but we should make sure just in case.
	std::string current_file = c_led->GetTree()->GetCurrentFile()->GetName();
	std::cout<<"making reference pure from file: "<<current_file<<std::endl;
	// then strip off preceding path
	if(current_file.find("/")!=std::string::npos){
		current_file = current_file.substr(current_file.find_last_of("/")+1,std::string::npos);
	}
	// then use filename as a key to look up the concentration and error
	if(calibration_data.count(current_file)==0){
		std::cerr<<"Couldn't find file "<<current_file<<" in calibration data map!"<<std::endl;
		for(auto&& am : calibration_data){
			std::cerr<<am.first<<", ";
		}
		std::cout<<std::endl;
		return false;
	}
	// if concentration for this file is not 0, refuse to make the pure from it
	double current_conc = calibration_data.at(current_file).first;
	if(current_conc!=0){
		std::cerr<<"Recorded concentration is not 0!"<<std::endl;
		return false;
	}
	*/
	
	// get nearest dark trace
	int dark_entry = GetNextDarkEntry(name, 0, false);
	c_dark->GetEntry(dark_entry);
	
	// subtract dark
	for(int j=0; j<n_datapoints; ++j){
		values.at(j) -= values_dark.at(j);
	}

	// build errors
	std::vector<double> xerrs, yerrs;
	for(int j=0; j<n_datapoints; ++j){
		xerrs.push_back(0.5*(wavelengths.at(1)-wavelengths.at(0)));
		double yerr = sqrt(errors.at(j)*errors.at(j) + errors_dark.at(j)*errors_dark.at(j));
		yerrs.push_back(yerr);
	}
	
	// make a TGraph from the data
	TGraphErrors* graph = new TGraphErrors(n_datapoints, wavelengths.data(), values.data(), xerrs.data(), yerrs.data());
	graph->SetName("Graph");
	graph->SetTitle("Graph");
	graph->SaveAs(purefile.c_str());
	
	return graph;
}

std::vector<double> Plotter::GetCalibCurve(std::string led, std::string fitmethod){
	
	std::string filename = std::string("calib_coeffs_")+led+"_"+fitmethod+".txt";
	std::cout<<"loading calibration parameters from file "<<filename<<std::endl;
	std::ifstream coeffs_file(filename.c_str());
	if(not coeffs_file.is_open()){
		std::cerr<<"couldn't open coefficients file "<<filename<<std::endl;
		return std::vector<double>{};
	}
	
	std::vector<double> coefficients;
	std::string line;
	while(getline(coeffs_file, line)){
		if(line.empty()) continue;
		if(line[0]=='#') continue;
		std::stringstream ss(line);
		double next_coeff;
		if(!(ss >> next_coeff)){
			std::cerr<<"Failed to parse line "<<line<<"!"<<std::endl;
			break;
		}
		coefficients.push_back(next_coeff);
	}
	coeffs_file.close();
	return coefficients;
	
}

/*  for posterity, for the time being, the old hard-coded calibration curves
			if(calib_ver=="old" && (name=="275_A" || name=="275_B")){
				std::vector<double> calib_coefficients{ -2.2420182e-05,
					                                 2.4347342,
					                                -10.671675,
					                                 25.117418,
					                                -15.640706,
					                                -35.283659,
					                                 67.871408 };
			} else if(calib_ver=="raw" && name=="275_A"){
				//std::vector<double> calib_coefficients{0.00517657, 2.86027, -7.36815, -19.955, 78.0564, 491.547, -1766.58}; // jul08
				std::vector<double> calib_coefficients{0.00141186, 2.85759, -7.43248, -18.2004, 78.8485, 472.124, -1814.29}; // jul12
				//std::vector<double> calib_coefficients{ 0.00396146,
				//	                                2.66054,
				//	                                -12.2432,
				//	                                42.909,
				//	                                -15.7365,
				//	                                -565.628,
				//	                                1362.37 };
				
				calib_curve.SetParameters(calib_coefficients.data());
			} else if(calib_ver=="simple" && name=="275_A"){
				//std::vector<double> calib_coefficients{0.00100964, 3.03952, -19.4798, 77.005, 27.6594, -1608.27, 4177.23}; // jul08
				std::vector<double> calib_coefficients{0.000901656, 2.7328, -11.972, 19.2158, 70.2304, -429.644, 743.848}; // jul12
				//std::vector<double> calib_coefficients{ 0.00669663,
				//	                                2.30659,
				//	                                -6.70807,
				//	                                -8.89578,
				//	                                65.2506,
				//	                                256.543,
				//	                                -1151.12 };
				
				calib_curve.SetParameters(calib_coefficients.data());
			} else if(calib_ver=="complex" && name=="275_A"){
				//std::vector<double> calib_coefficients{0.00532017, 2.62183, -8.00979, -18.5983, 91.0425, 516.762, -1989.37}; // jul08
				std::vector<double> calib_coefficients{0.00256741, 2.53608, -5.59115, -34.2931, 86.1831, 887.119, -2934}; // jul12
				//std::vector<double> calib_coefficients{ 0.00668535,
				//	                                2.02915,
				//	                                16.099,
				//	                                -501.567,
				//	                                4527.85,
				//	                                -17774,
				//	                                25663.8 };
				
				calib_curve.SetParameters(calib_coefficients.data());
			} else if(calib_ver=="raw" && name=="275_B"){
				//std::vector<double> calib_coefficients{0.0121442, 2.90637, -7.30414, -19.9099, 77.1875, 485.482, -1764.29}; // jul08
				std::vector<double> calib_coefficients{0.00151493, 2.67998, -6.92564, -16.6946, 75.3407, 444.124, -1731.61}; // jul12
				//std::vector<double> calib_coefficients{ 0.0162978,
				//	                                2.57888,
				//	                                -10.2983,
				//	                                26.7082,
				//	                                14.5267,
				//	                                -414.93,
				//	                                875.207 };
				
				calib_curve.SetParameters(calib_coefficients.data());
			} else if(calib_ver=="simple" && name=="275_B"){
				//std::vector<double> calib_coefficients{0.00958699, 3.01173, -17.3047, 56.2216, 55.7341, -1232.64, 2979.2}; // jul08
				std::vector<double> calib_coefficients{0.000557654, 2.58452, -12.3119, 31.4942, 55.7869, -735.182, 1678.52}; // jul12
				/std::vector<double> calib_coefficients{ 0.0194698,
				//	                                2.28874,
				//	                                -6.76499,
				//	                                -9.3567,
				//	                                64.0635,
				//	                                261.506,
				//	                                -1119.84 };
				
				calib_curve.SetParameters(calib_coefficients.data());
			} else if(calib_ver=="complex" && name=="275_B"){
				//std::vector<double> calib_coefficients{0.00603398, 2.88702, -9.77863, -18.9483, 116.117, 594.401, -2502.56}; // jul08
				std::vector<double> calib_coefficients{0.00197077, 2.35596, -5.21091, -30.7562, 100.809, 708.71, -2614.03};  // jul12
				//std::vector<double> calib_coefficients{ 0.0204775,
				//	                                1.79824,
				//	                                19.8836,
				//	                                -516.837,
				//	                                4389.5,
				//	                                -16604.5,
				//	                                23385.5 };
				
				calib_curve.SetParameters(calib_coefficients.data());
			}
*/

bool Plotter::DoRawFit(TGraphErrors* g_abs, std::pair<double,double>& peak_posns, std::pair<double,double>& peak_heights, std::pair<double,double>& peak_errs, std::string pngname, bool draw){
	g_abs->GetPoint(sample_273, peak_posns.first, peak_heights.first);
	g_abs->GetPoint(sample_276, peak_posns.second, peak_heights.second);
	peak_errs.first = errors.at(sample_273);
	peak_errs.second = errors.at(sample_276);
	//std::cout<<"raw peaks are: "<<peak_heights.first<<", "<<peak_heights.second<<std::endl;
	
	bool ok = true;
	
	// check if extracted peak amplitudes are positive
	if(peak_heights.first<0 || peak_heights.second<0){
		// this tends to only happen for the first value, and results in a point at (0,0)
		// plotting this point then tends to lie out of the trend, so may throw the fit off.
		// instead, just skip this concentration...?
		//std::cerr<<"skipping concentration measurement "<<i<<" as one of the peaks is negative"<<std::endl;
		//continue;
		// or alternatively coerce to 0
		//peak_heights.first=std::max(0.,peak_heights.first);
		//peak_heights.second=std::max(0.,peak_heights.second);
		
		ok = false;
	}
	
	if(pngname!=""){
		if(!draw) gROOT->SetBatch(true);
		TCanvas* c_temp=(TCanvas*)gROOT->FindObject("c_temp");
		c_temp->cd();
		g_abs->Draw("AP");
		c_temp->SaveAs(pngname.c_str());
		gROOT->SetBatch(false);
	}
	
	return true;
}

bool Plotter::DoSimpleFit(TGraph* abs_graph, std::pair<double,double>& peak_posns, std::pair<double,double>& peak_heights, std::pair<double,double>& peak_errs, std::string pngname, bool plot){
	// given the absorption graph, fit the two main peaks at 273 and 275nm
	// with gaussians, but only over a narrow region around the peak centre,
	// where a gaussian approximation is reasonable
	TF1 gaus1("gaus1","gaus",272.5,273.5);   // very narrow.... too narrow?
	TF1 gaus2("gaus2","gaus",275.2,276.2);   // very narrow.... too narrow?
	
	double peakval = *std::max_element(abs_graph->GetY(),abs_graph->GetY()+abs_graph->GetN());
	
	// defaults and limits
	//gaus1.SetParameter(0,0.5);
	gaus1.SetParameter(0,peakval);
	gaus1.SetParameter(1,273);
	gaus1.SetParameter(2,0.6);
	gaus1.SetParLimits(0,0.,1.);
	gaus1.SetParLimits(1,272.75,273.25);
	gaus1.SetParLimits(2,0.3,1.);
	
	//gaus2.SetParameter(0,0.3);
	gaus2.SetParameter(0,peakval*0.1);
	gaus2.SetParameter(1,275.65);
	gaus2.SetParameter(2,0.55);
	gaus2.SetParLimits(0,0.,1.);
	gaus2.SetParLimits(1,275.3,276.0);
	gaus2.SetParLimits(2,0.3,1.);
	
	// fit the two gaussians
	TFitResultPtr gfptr1 = abs_graph->Fit("gaus1","NRQS");
	TFitResultPtr gfptr2 = abs_graph->Fit("gaus2","NRQS");
	
	double gausamp1 = gaus1.GetParameter(0);
	double gausamp2 = gaus2.GetParameter(0);
	
	double gaus1amperr = gaus1.GetParError(0);
	double gaus2amperr = gaus2.GetParError(0);
	//std::cout<<"simple peak fit gaus 1 has ampltiude "<<gausamp1<<"+-"<<gaus1amperr<<std::endl;
	//std::cout<<"simple peak fit gaus 2 has ampltiude "<<gausamp2<<"+-"<<gaus2amperr<<std::endl;
	
	double gaus1pos = gaus1.GetParameter(1);
	double gaus2pos = gaus2.GetParameter(1);
	
	peak_posns = std::pair<double,double>{gaus1pos,gaus2pos};
	peak_heights = std::pair<double,double>{gausamp1,gausamp2};
	peak_errs = std::pair<double,double>{gaus1amperr,gaus2amperr};
	
	bool ok = 1;
	// note that these are not necessarily reliable indicators of whether the fit was ok or not
	if(gfptr1->IsEmpty() || !gfptr1->IsValid() || gfptr1->Status()!=0){
		std::cerr<<"gaus1 fit failed"<<std::endl;
		ok = 0;
	}
	if(gfptr2->IsEmpty() || !gfptr2->IsValid() || gfptr2->Status()!=0){
		std::cerr<<"gaus2 fit failed"<<std::endl;
		ok = 0;
	}
	
	//std::cout<<"simple peaks are: "<<peak_heights.first<<", "<<peak_heights.second<<std::endl;
	
	if(pngname!=""){
		if(!plot) gROOT->SetBatch(true);
		TCanvas* c_temp=(TCanvas*)gROOT->FindObject("c_temp");
		c_temp->cd();
		abs_graph->Draw("AP");
		gaus1.Draw("same");
		gaus2.Draw("same");
		c_temp->SaveAs(pngname.c_str());
		gROOT->SetBatch(false);
	}
	
	return ok;
}

bool Plotter::DoComplexFit(TGraphErrors* g_abs, std::pair<double,double>& peak_posns, std::pair<double,double>& peak_heights, std::pair<double,double>& peak_errs, std::pair<double,double> initial_peak_heights, std::string pngname, bool draw){
	
	TF1* abs_func = GetAbsFunc();
	// Seed the initial values based on the raw fit.
	// The raw heights will be over-estimates, since in the complex fit the two main peaks
	// are fit with additional shoulder gaussians, so the main component amplitude is less,
	// but it should be close enough for a starting value.
	// in case of a bad initial value from low concentrations, coerce to at least 0.
	initial_peak_heights.first = std::max(0.,initial_peak_heights.first);
	initial_peak_heights.second = std::max(0.,initial_peak_heights.second);
	abs_func->SetParameter("peak 1 amp",initial_peak_heights.first);
	abs_func->SetParameter("peak 2 amp",initial_peak_heights.second/initial_peak_heights.first);
	TFitResultPtr frptr = g_abs->Fit(abs_func,"RQNS");
	
	bool ok = true;
	// note these are not reliable indicators that the fit is bad...
	if( frptr->IsEmpty() || !frptr->IsValid() || frptr->Status()!=0){
			// fit failed; skip value?
			std::cout<<"FIT INVALID!"<<std::endl;
			std::cout<<"empty: "<<frptr->IsEmpty()<<", valid: "<<frptr->IsValid()
			         <<", status: "<<frptr->Status()<<std::endl;
			//frptr->Print();
			//continue;
			ok = false;
	}
	
	// extract the height of the peaks by finding the maximum of the curve
	// within a region around the peak central location. Note we CANNOT
	// just take the gaussian amplitude because there are multiple overlapping
	// contributions to the peak here. We could just evaluate the function at
	// the known peak position, but this allows more flexibility.
	peak_posns.first = abs_func->GetMaximumX(272.5,273.5);
	peak_posns.second = abs_func->GetMaximumX(275.,276.);
	// extract the peak heights by evaluating the function at that x
	peak_heights.first   = abs_func->Eval(peak_posns.first);
	peak_heights.second  = abs_func->Eval(peak_posns.second);
	
	// this is a pretty robust indicator that the fit is bad!
	if(TMath::IsNaN(peak_heights.first)||TMath::IsNaN(peak_heights.second)){
		std::cout<<"NaN peak height!"<<std::endl;
		//continue;
		if(TMath::IsNaN(peak_heights.first)) peak_heights.first=0;
		if(TMath::IsNaN(peak_heights.second)) peak_heights.second=0;
		ok = false;
	}
	//std::cout<<"complex peaks are: "<<peak_heights.first<<", "<<peak_heights.second<<std::endl;
	
	// calculate error on peak heights
	peak_errs = CalculateError(abs_func, peak_posns.first, peak_posns.second);
	
	// if given an output name, save an image of the fit
	if(pngname!=""){
		if(!draw) gROOT->SetBatch(true);
		TCanvas* c_temp=(TCanvas*)gROOT->FindObject("c_temp");
		c_temp->cd();
		g_abs->Draw("AL");
		abs_func->Draw("same");
		c_temp->SaveAs(pngname.c_str());
		gROOT->SetBatch(false);
	}
	
	/*
	if(draw){
		// draw the absorption curve and the complex fit.
		TCanvas cnew("cnew","cnew",1024,800);
		g_abs->Draw("ALP");
		abs_func->Draw("same");
		// these have to go *after* a draw or the graph has no axes
		//g_abs->GetYaxis()->SetRangeUser(-0.1,0.6);
		g_abs->GetXaxis()->SetRangeUser(240,320);
		cnew.Modified();
		cnew.Update();
		gSystem->ProcessEvents();
		std::cout<<"waiting for user to close canvas"<<std::endl;
		while(gROOT->FindObject("cnew")!=nullptr){
			gSystem->ProcessEvents();
			std::this_thread::sleep_for(std::chrono::milliseconds(100));
		}
		std::this_thread::sleep_for(std::chrono::milliseconds(100));
	}
	*/
	
	delete abs_func;
	
	return ok;
	
}

std::vector<std::string> Plotter::GetListOfFiles(std::string inputdir){
	// convert to absolute path if required
	std::string absdir = GetStdoutFromCommand(std::string("readlink -f ")+inputdir);
	absdir.erase(absdir.find_last_not_of(" \t\n\015\014\013")+1);  // strip trailing whitespace
	//std::cout<<"input dir is '"<<inputdir<<"', absolute path is '"<<absdir<<"'"<<std::endl;
	std::string lscommand = "find " + absdir + " -iname '*.root' 2>/dev/null";
	//std::cout<<"lscommand is '"<<lscommand<<"'"<<std::endl;
	std::string fileliststring = GetStdoutFromCommand(lscommand);
	//std::cout<<"fileliststring is '"<<fileliststring<<"'"<<std::endl;
	
	std::stringstream ssl;
	ssl << fileliststring;
	std::string nextfilestring;
	std::vector<std::string> paths;
	
	while(getline(ssl,nextfilestring)){
		std::size_t last_char_loc = nextfilestring.find_last_of("/\\");
		std::string fname = nextfilestring.substr(last_char_loc+1);
		paths.push_back(fname);
	}
	return paths;
}

std::string Plotter::GetStdoutFromCommand(std::string cmd, int bufsize){
	/*
	  credit: Jeremy Morgan, source:
	  https://www.jeremymorgan.com/tutorials/c-programming/how-to-capture-the-output-of-a-linux-command-in-c/
	*/
	std::string data;
	FILE * stream;
	char* buffer = new char[bufsize];
	cmd.append(" 2>&1");
	
	stream = popen(cmd.c_str(), "r");
	if(stream){
		while(!feof(stream)){
			if (fgets(buffer, bufsize, stream) != NULL) data.append(buffer);
		}
		pclose(stream);
	}
	delete[] buffer;
	return data;
}
