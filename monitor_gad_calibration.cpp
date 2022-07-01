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

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <thread>
#include <cassert>


class Plotter{
	public:
	
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
	int refentry=0;
	int darksperloop=12;
	
	std::map<std::string, std::pair<double,double>> calibration_data;
	
	int LoadNewData();
	int SetBranchAddresses();
	int GetNextDarkEntry(std::string name, int ledon_entry_num, bool check=true);
	
	bool LoadConcentrations(std::string filename);
	TF1* PureScaledPlusExtras();
	TMultiGraph* ExtractAbsorbance(std::string name, std::string calib_version="old");
	std::pair<double,double> CalculateError(TF1* abs_func, double peak1_pos, double peak2_pos);
	int FitTwoGaussians(TGraph* abs_graph, std::pair<double,double>& simple_peaks, std::pair<double,double>& simple_errs, std::pair<double,double>& simple_posns, bool plot=false);
	TF1* GetAbsFunc();
	TMultiGraph* FitCalibrationData(std::string name, std::string dataset);
	
	void SetTimeAxis(TH1* hist, long t0_seconds);
	
	int sample_273=0;
	int sample_276=0;
	
	//TF1* purefunc=nullptr;
	std::string purefile="../GDConcMeasure/pureDarkSubtracted.root";
	
	int Initialise();
	TMultiGraph* Execute();
};

int main(){
	
	TApplication myapp("rootTApp",0,0);
	TCanvas* c1 = new TCanvas("c1","c1",1024,800);  // do not change the name of this canvas, closing it will terminate the program.
	
	// create plotter class
	Plotter myplotter;
	
	// create toolchains and set branch addresses
	myplotter.Initialise();
	int last_num_files=0;
	
	// multigraph of the calibration curves so far
	TMultiGraph* mg_cals=nullptr;
	
	// keep looking for new files until user closes main canvas
	while(c1!=nullptr){
		
		// look for new files
		std::cout<<"loading new data"<<std::endl;
		myplotter.LoadNewData();
		
		// if new files found...
		if(myplotter.files.size()!=last_num_files){
			std::cout<<"new data: deleting graphs"<<std::endl;
			// delete the old calibration curves
			if(mg_cals) delete mg_cals;
			
			// re-create with the full dataset
			std::cout<<"making new graphs"<<std::endl;
			mg_cals = myplotter.Execute();
			std::cout<<"drawing"<<std::endl;
			mg_cals->Draw("AL");
			c1->Modified();
			c1->Update();
			
			last_num_files = myplotter.files.size();
		}
		
		gSystem->ProcessEvents();
		std::this_thread::sleep_for(std::chrono::milliseconds(100));
		c1 = (TCanvas*)gROOT->FindObject("c1");
	}
	
	delete mg_cals;
	
	return 0;
}

int Plotter::LoadNewData(){
	
	// scan folder for new files
	const char* dirname = "../GDConcMeasure/data/2022/06/rename";
//	std::cout<<"scanning "<<dirname<<" for new files"<<std::endl;
	TSystemDirectory dir(dirname, dirname);
	TList* filestlist = dir.GetListOfFiles();
	
	// Unfortunately if you call TChain::Add with the same file twice
	// (or patterns that match the same file twice),
	// it'll get re-added and you'll have duplicate entries.
	// so we manually keep track of what's in the chain
	// and re-scan the directory, adding only what's new
	if(filestlist){
//		std::cout<<"got a file list"<<std::endl;
		TSystemFile* file;
		TIter next(filestlist);
		while( file = (TSystemFile*)next() ) {
			TString fname = file->GetName();
//			std::cout<<"checking file "<<fname<<std::endl;
			// TODO check name matches *01Jun22Calib_[0-9]*
			if( !file->IsDirectory() && (fname.EndsWith(".root")) && (files.count(fname.Data())==0) ){
				std::string thisfname = std::string(dirname) + "/" + std::string(fname.Data());
//				std::cout<<"adding "<<thisfname<<" to tchains"<<std::endl;
				// add the new file
				for(auto&& c : chains){
					c.second->Add(thisfname.c_str());
					c_dark->Add(thisfname.c_str());
				}
				
				// note this file is now added
				files.emplace(fname.Data(),1);
			}
			
			// update num measurements
			num_meas = c_275_A->GetEntries();
//			std::cout<<"had "<<num_meas<<" measurements"<<std::endl;
		}
	}
	
	// get number of measurements
	num_meas = c_275_A->GetEntries();
	std::cout<<"had "<<num_meas<<" measurements"<<std::endl;
	if(num_meas==0) return 0;
	
	// get num datapoints
	std::cout<<"getting n datapoints"<<std::endl;
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
	
	std::ifstream calib_data_file(filename.c_str());
	if(not calib_data_file.is_open()){
		std::cerr<<"couldn't open calibration data file "<<filename<<std::endl;
		return false;
	}
	
	calibration_data.clear();
	std::string line;
	while(getline(calib_data_file, line)){
		if(line.empty()) continue;
		if(line[0]=='#') continue;
		std::stringstream ss(line);
		std::string fname, conc, concerr;
		if(!(ss >> fname >> conc >> concerr)){
			std::cerr<<"Failed to parse line "<<line<<"!"<<std::endl;
			break;
		}
		calibration_data.emplace(fname, std::pair<double,double>{strtod(conc.c_str(), nullptr),strtod(concerr.c_str(), nullptr)});
	}
	calib_data_file.close();
	
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

int Plotter::Initialise(){
	
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
	
	return 0;
}

TMultiGraph* Plotter::Execute(){
	
	// step 1: analyse calibration data for polynomial coffiecients,
	// take the results printed out and populate the coefficients list in ExtractAbsorbance
	std::string dataset="jun01_calibration_info.txt";
	TMultiGraph* mg_cals_A = FitCalibrationData("275_A", dataset);
	TMultiGraph* mg_cals_B = FitCalibrationData("275_B", dataset);
	
	TMultiGraph* mg_cals = new TMultiGraph("mg_cals","mg_cals");
	mg_cals->Add(mg_cals_A);
	mg_cals->Add(mg_cals_B);
	
	/*
	mg_cals->Draw("AL");
	
	// step 2: analyse data to extract plots of concentration stability
	// we can do it for each of the calibration fit methods
	TMultiGraph* mg_all = new TMultiGraph("mg_all","mg_all");
	mg_all->Add(ExtractAbsorbance("275_A", "raw"));
	mg_all->Add(ExtractAbsorbance("275_A", "simple"));
	mg_all->Add(ExtractAbsorbance("275_A", "complex"));
	mg_all->Add(ExtractAbsorbance("275_B", "raw"));
	mg_all->Add(ExtractAbsorbance("275_B", "simple"));
	mg_all->Add(ExtractAbsorbance("275_B", "complex"));
	
	// for comparison, the true concentrations for input
	bool ok = LoadCalibrationData("feb24_calibration_info.txt");
	std::vector<double> calib_concs;
	for(auto&& acalibmeas : calibration_data){
		calib_concs.push_back(acalibmeas.second.first);
	}
	std::cout<<"we had "<<calib_concs.size()<<" calibration measurements and "
	         <<num_meas<<" data measurements"<<std::endl;
	std::vector<double> numberline(num_meas);
	std::iota(numberline.begin(),numberline.end(),0);
	TGraph* calgraph = new TGraph(calib_concs.size(), numberline.data(), calib_concs.data());
	calgraph->SetLineColor(kMagenta);
	calgraph->SetMarkerColor(kMagenta);
	calgraph->SetMarkerStyle(30);
	mg_all->Add(calgraph);
	
	c1->cd();
	mg_all->Draw("ALP");
	*/
	
	
	return mg_cals;
}

/////////////////////////////////

TGraphErrors* dark_subtracted_pure =nullptr;
double PureFuncv2(double* x, double* par){
	if(dark_subtracted_pure==nullptr){
		TFile *_file0 = TFile::Open("../GDConcMeasure/pureDarkSubtracted.root");
		dark_subtracted_pure = (TGraphErrors*)_file0->Get("Graph");
	}
	// par [0] = y-scaling
	// par [1] = x-scaling
	// par [2] = x-offset
	// par [3] = y-offset
	// par [4] = linear baseline offset
	// par [5] = shoulder gaussian scaling
	// par [6] = shoulder gaussian centre, restricted to > 282nm (RH shoulder)
	// par [7] = shoulder gaussian spread
	double purepart = (par[0] * dark_subtracted_pure->Eval((par[1]*x[0])+par[2]));
	double linpart = (par[4] * x[0]) + par[3];
	double shoulderpart = par[5]*exp(-0.5*TMath::Sq((x[0]-282.-abs(par[6]))/par[7]));
	double retval = purepart + linpart + shoulderpart;
	
	return retval;
}

TF1* Plotter::PureScaledPlusExtras(){
	
	// construct functional fit. We'll scale and add a linear background.
	// limit the pure function to a region in which we have light - no point fitting outside this region
	static int purever=0;
	purever++;
	std::string name="purev2_fct"+std::to_string(purever);
	const int wave_min = 260, wave_max = 300, numb_of_fitting_parameters = 8;
	TF1* pure = new TF1(name.c_str(), PureFuncv2, wave_min, wave_max, numb_of_fitting_parameters);
	// set default parameters
	pure->SetParameters(1.,1.,0.,0.,0.,0.,0.,10.);
	
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
	
	return std::pair<double,double>{sqrt(totalerror1),sqrt(totalerror2)};
}

int Plotter::FitTwoGaussians(TGraph* abs_graph, std::pair<double,double>& simple_peaks, std::pair<double,double>& simple_errs, std::pair<double,double>& peak_posns, bool plot){
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
	std::cout<<"simple peak fit gaus 1 has ampltiude "<<gausamp1<<"+-"<<gaus1amperr<<std::endl;
	std::cout<<"simple peak fit gaus 2 has ampltiude "<<gausamp2<<"+-"<<gaus2amperr<<std::endl;
	
	double gaus1pos = gaus1.GetParameter(1);
	double gaus2pos = gaus2.GetParameter(1);
	
	simple_peaks = std::pair<double,double>{gausamp1,gausamp2};
	simple_errs = std::pair<double,double>{gaus1amperr,gaus2amperr};
	peak_posns = std::pair<double,double>{gaus1pos,gaus2pos};
	
	int ok = 1;
	if(gfptr1->IsEmpty() || !gfptr1->IsValid() || gfptr1->Status()!=0){
		std::cerr<<"gaus1 fit failed"<<std::endl;
		ok = 0;
	}
	if(gfptr2->IsEmpty() || !gfptr2->IsValid() || gfptr2->Status()!=0){
		std::cerr<<"gaus2 fit failed"<<std::endl;
		ok = 0;
	}
	return ok;
}

TMultiGraph* Plotter::ExtractAbsorbance(std::string name, std::string calib_ver){
	
	
	TF1* pureplusextras = PureScaledPlusExtras();
	
	TChain* c_led = chains.at(name);
	TChain* c_dark = chains.at("dark");
	int offset = offsets.at(name);
	int n_entries = c_led->GetEntries();
	
	std::string title = name+"_g";
	TMultiGraph* mg = new TMultiGraph(title.c_str(), title.c_str());
	
	// FIXME TODO
	//TGraphErrors* g_concentrations = new TGraphErrors(n_entries);
	TGraph* g_concentrations = new TGraph(n_entries);
	if(calib_ver=="raw"){
		g_concentrations->SetMarkerColor(kSpring-5);
		g_concentrations->SetLineColor(kSpring-5);
	} else if(calib_ver=="simple"){
		g_concentrations->SetMarkerColor(kRed);
		g_concentrations->SetLineColor(kRed);
	} else if(calib_ver=="complex"){
		g_concentrations->SetMarkerColor(kBlue);
		g_concentrations->SetLineColor(kBlue);
	}
	if(name=="275_A"){
		g_concentrations->SetMarkerStyle(20);
	} else if(name=="275_B"){
		g_concentrations->SetMarkerStyle(34);
	}
	
	TGraph g_relative_heights(n_entries);
	TGraph g_relative_ratios(n_entries);
	
	int n_concentration_vals=-1;
	for(int i=0; i<n_entries; i++){
		
		std::cout<<i<<"\n";
		
		// get LED on trace
		int nbytes = c_led->GetEntry(i);
		if(nbytes<=0){
			std::cerr<<"COULDN'T LOAD ENTRY "<<i<<" FROM LED CHAIN"<<std::endl;
			exit(-1);
		}
		
		// get nearest dark trace
		int dark_entry = offset + i*darksperloop;
		//int dark_entry = GetNextDarkEntry(name, i, false);
		c_dark->GetEntry(dark_entry);
		
		// subtract dark
		for(int j=0; j<n_datapoints; ++j){
			values.at(j) -= values_dark.at(j);
		}
		
		// split into in- and -out of absorbance bands
		std::vector<double> inband_values;
		std::vector<double> inband_wls;
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
		
		// fit with pure scaled
		g_sideband->Fit(pureplusextras,"RNQ");
		
		// calculate absorbance from ratio of fit to data in absorption region
		TGraph* g_abs = new TGraph(inband_values.size());
		std::string restitle = "g_abs_"+std::to_string(i);
		for(int k=0; k<inband_values.size(); ++k){
			double wlval, dataval;
			g_inband->GetPoint(k, wlval, dataval);
			double fitval = pureplusextras->Eval(wlval);
			g_abs->SetPoint(k, wlval, log10(fitval/dataval));
		}
		
		// extract the peak heights and errors
		std::pair<double,double> peak_posns;
		std::pair<double,double> peak_heights;
		std::pair<double,double> peak_errs;
		
		// raw peaks: just the graph value at the peak wavelengths
		if(calib_ver=="raw" || calib_ver=="complex"){
			peak_posns.first = wavelengths.at(sample_273);
			peak_posns.second = wavelengths.at(sample_276);
			double wlval;
			g_abs->GetPoint(sample_273, wlval, peak_heights.first);
			g_abs->GetPoint(sample_276, wlval, peak_heights.second);
			
			// if peaks are negative, coerce to 0
			peak_heights.first=std::max(0.,peak_heights.first);
			peak_heights.second=std::max(0.,peak_heights.second);
			
			peak_errs = std::pair<double,double>{0.1,0.1}; // FIXME
		}
		
		// simple peaks: fit the peaks with gaussians
		if(calib_ver=="simple"){
			int ok = FitTwoGaussians(g_abs, peak_heights, peak_errs, peak_posns);
		}
		
		// complex peaks: fit entire region with a combination of 4 gaussians
		TF1* abs_func=nullptr;
		if(calib_ver=="complex"){
			
			abs_func = GetAbsFunc();  // we own the returned TF1
			abs_func->SetParameter("peak 1 amp",peak_heights.first);
			abs_func->SetParameter("peak 2 amp",peak_heights.second/peak_heights.first);
			TFitResultPtr frptr = g_abs->Fit(abs_func,"RQNS");
			
			if(frptr->IsEmpty() || !frptr->IsValid() || frptr->Status()!=0){
				// fit failed; skip value?
				//continue;
				peak_heights.first = 0;
				peak_heights.second = 0;
			} else {
				// get the height of the peaks as the maximum of the curve around the peak location
				peak_posns.first = abs_func->GetMaximumX(272.5,273.5);
				peak_posns.second = abs_func->GetMaximumX(275.,276.);
				peak_heights.first  = abs_func->Eval(peak_posns.first);
				peak_heights.second  = abs_func->Eval(peak_posns.second);
				if(TMath::IsNaN(peak_heights.first)||TMath::IsNaN(peak_heights.second)){
					//continue;
					peak_heights.first = 0;
					peak_heights.second = 0;
				}
				
				std::cout<<"complex peaks at: "<<peak_posns.first<<":"<<peak_posns.second
				         <<" are "<<peak_heights.first<<" and "<<peak_heights.first
				         <<" with diff "<<peak_heights.first-peak_heights.second
				         <<" and ratio "<<peak_heights.first/peak_heights.second<<std::endl;
				
				peak_errs = CalculateError(abs_func, peak_posns.first, peak_posns.second);
			}
		}
		++n_concentration_vals;
		
		g_relative_heights.SetPoint(n_concentration_vals,n_concentration_vals,
		                            peak_heights.first-peak_heights.second);
		double peak_ratio = (peak_heights.second>0) ? peak_heights.first/peak_heights.second : 1;
		g_relative_ratios.SetPoint(n_concentration_vals,n_concentration_vals,peak_ratio);
		
		// total error on difference in peak heights
		double peakheightdifferr = sqrt(TMath::Sq(peak_errs.first)+TMath::Sq(peak_errs.second));
		
		// convert into a concentration value. We need the calibration curve:
		TF1 calib_curve("calib", "pol6", 0, 0.4);
		if(calib_ver=="old" && (name=="275_A" || name=="275_B")){
			// did we have calibration coefficients for 275_A vs 275_B?
			std::vector<double> calib_coefficients{ -2.2420182e-05,
			                                         2.4347342,
			                                        -10.671675,
			                                         25.117418,
			                                        -15.640706,
			                                        -35.283659,
			                                         67.871408 };
			calib_curve.SetParameters(calib_coefficients.data());
		} else if(calib_ver=="raw" && name=="275_A"){
			std::vector<double> calib_coefficients{ 0.00396146,
			                                        2.66054,
			                                        -12.2432,
			                                        42.909,
			                                        -15.7365,
			                                        -565.628,
			                                        1362.37 };
			calib_curve.SetParameters(calib_coefficients.data());
		} else if(calib_ver=="simple" && name=="275_A"){
			std::vector<double> calib_coefficients{ 0.00669663,
			                                        2.30659,
			                                        -6.70807,
			                                        -8.89578,
			                                        65.2506,
			                                        256.543,
			                                        -1151.12 };
			calib_curve.SetParameters(calib_coefficients.data());
		} else if(calib_ver=="complex" && name=="275_A"){
			std::vector<double> calib_coefficients{ 0.00668535,
			                                        2.02915,
			                                        16.099,
			                                        -501.567,
			                                        4527.85,
			                                        -17774,
			                                        25663.8 };
			calib_curve.SetParameters(calib_coefficients.data());
		} else if(calib_ver=="raw" && name=="275_B"){
			std::vector<double> calib_coefficients{ 0.0162978,
			                                        2.57888,
			                                        -10.2983,
			                                        26.7082,
			                                        14.5267,
			                                        -414.93,
			                                        875.207 };
			calib_curve.SetParameters(calib_coefficients.data());
		} else if(calib_ver=="simple" && name=="275_B"){
			std::vector<double> calib_coefficients{ 0.0194698,
			                                        2.28874,
			                                        -6.76499,
			                                        -9.3567,
			                                        64.0635,
			                                        261.506,
			                                        -1119.84 };
			calib_curve.SetParameters(calib_coefficients.data());
		} else if(calib_ver=="complex" && name=="275_B"){
			std::vector<double> calib_coefficients{ 0.0204775,
			                                        1.79824,
			                                        19.8836,
			                                        -516.837,
			                                        4389.5,
			                                        -16604.5,
			                                        23385.5 };
			calib_curve.SetParameters(calib_coefficients.data());
		}
		
		// solve for concentration (x) from absorbance (y), with 0.01 < x < 0.21
		double conc = calib_curve.GetX(peak_heights.first - peak_heights.second, 0.001, 0.25);
		g_concentrations->SetPoint(n_concentration_vals,n_concentration_vals,conc);
		// error on concentration is error on height times gradient at that point FIXME TODO
//		double conc_err = peakheightdifferr * calib_curve.Derivative(peak_heights.first - peak_heights.second);
//		g_concentrations->SetPointError(i,wavelength_errors.at(i),conc_err);
		
		/*
		TCanvas cnew("cnew","cnew",1024,800);
		std::string mgtitle = "mgg_"+std::to_string(i);
		TMultiGraph* mgg = new TMultiGraph(mgtitle.c_str(),mgtitle.c_str());
		mgg->Add(g_abs);
		mgg->Draw("ALP");
		abs_func->Draw("same");
		// these have to go *after* a draw or the graph has no axes
		mgg->GetYaxis()->SetRangeUser(-0.1,0.6);
		mgg->GetXaxis()->SetRangeUser(240,320);
		cnew.Modified();
		cnew.Update();
		gSystem->ProcessEvents();
		
		std::cout<<"waiting for user to close canvas"<<std::endl;
		while(gROOT->FindObject("cnew")!=nullptr){
			gSystem->ProcessEvents();
			std::this_thread::sleep_for(std::chrono::milliseconds(100));
		}
		std::this_thread::sleep_for(std::chrono::milliseconds(100));
		*/
		
		delete g_abs;       // use this or above
		
		//delete mgg;
		delete g_inband;
		delete g_sideband;
		delete g_other;
		delete abs_func;
	}
	g_concentrations->Set(n_concentration_vals);
	
	/*
	TCanvas cc("cc","cc",1024,800);
	g_concentrations->Draw("ALP");
	cc.Modified();
	cc.Update();
	
	TCanvas ccc("ccc","ccc",1024,800);
	g_relative_heights.Draw("ALP");
	ccc.Modified();
	ccc.Update();
	
	TCanvas c4("c4","c4",1024,800);
	g_relative_ratios.Draw("ALP");
	g_relative_ratios.GetYaxis()->SetRangeUser(0,2);
	c4.Modified();
	c4.Update();
	
	gSystem->ProcessEvents();
	while(gROOT->FindObject("ccc")!=nullptr){
		gSystem->ProcessEvents();
		std::this_thread::sleep_for(std::chrono::milliseconds(100));
	}
	*/
	
	std::cout<<"returning multigraph"<<std::endl;
	mg->Add(g_concentrations);
	return mg;
	
}

//####################################

TMultiGraph* Plotter::FitCalibrationData(std::string name, std::string dataset){
	
	// same basic procedure, loop over the data, do dark subtraction,
	// fit the dark-subtracted data, extract the absorption graph,
	// fit the absorption graph and pull the peak heights and difference.
	// only difference is this time rather than using peak height difference
	// to look up concentration, we look up the known concentration,
	// and add a point mapping peak height to concentration.
	
	// load calibration info that specifies our concentrations
	bool ok = LoadConcentrations(dataset);
	
	TF1* pureplusextras = PureScaledPlusExtras();
	
	TChain* c_led = chains.at(name);
	TChain* c_dark = chains.at("dark");
	int n_entries = c_led->GetEntries();
	
	std::string title = name+"_g";
	TGraphErrors* calib_curve_raw = new TGraphErrors(n_entries);
	std::string title0 = title+"_raw";
	calib_curve_raw->SetName(title0.c_str());
	calib_curve_raw->SetTitle(title0.c_str());
	TGraphErrors* calib_curve_simple = new TGraphErrors(n_entries);
	std::string title1 = title+"_simple";
	calib_curve_simple->SetName(title1.c_str());
	calib_curve_simple->SetTitle(title1.c_str());
	TGraphErrors* calib_curve_complex = new TGraphErrors(n_entries);
	std::string title2 = title+"_complex";
	calib_curve_complex->SetName(title2.c_str());
	calib_curve_complex->SetTitle(title2.c_str());
	TMultiGraph* mg_calibs = new TMultiGraph(title.c_str(), title.c_str());
	mg_calibs->Add(calib_curve_simple);
	mg_calibs->Add(calib_curve_complex);
	mg_calibs->Add(calib_curve_raw);
	
	// as a back-check we'll also plot the concentrations vs measurement
	// and peak height difference vs measurement indepdently
	// (as well as the usual one against the other)
	std::vector<double> concentrations, concentration_errs;
	
	
	int next_graph_point=0;
	double lastconc=-99;
	//std::cout<<"looping over entries"<<std::endl;
	for(int i=0; i<n_entries; i++){
		
		std::cout<<i<<": ";
		
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
		
		// fit with pure scaled
		g_sideband->Fit(pureplusextras,"RNQ");
		
		std::cout<<"YYY"<<std::endl;
		// calculate absorbance from ratio of fit to data in absorption region
		TGraphErrors* g_abs = new TGraphErrors(inband_values.size());
		std::string restitle = "g_abs_"+std::to_string(i);
		for(int k=0; k<inband_values.size(); ++k){
			double wlval, dataval;
			g_inband->GetPoint(k, wlval, dataval);
			double fitval = pureplusextras->Eval(wlval);
			std::cout<<"setting point to log10("<<fitval<<"/"<<dataval<<") = "<<log10(fitval/dataval)<<std::endl;
			g_abs->SetPoint(k, wlval, log10(fitval/dataval));
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
//			/*
			std::cout<<"data val: "<<inband_values.at(k)<<"+- "<<error_on_data<<"\n"
			         <<", fractional error: "<<(error_on_data/dataval)<<"\n"
			         <<", fit val: "<<fitval<<", +- "<<error_on_fitval<<"\n"
			         <<", fractional error: "<<(error_on_fitval/fitval)<<"\n"
			         <<", total error on ratio: "<<err_on_ratio
			         <<", ratio of abs / transmitted: "<<(fitval/dataval)
			         <<", ratio of error to value: "<<err_on_ratio_over_ratio<<std::endl;
//			*/
			g_abs->SetPointError(k, wavelength_errors.at(0)/2., err_on_ratio_over_ratio*(1./log(10.)));
		}
		
		// as the simplest estimate of peak heights, just take the graph value at the peaks.
		std::pair<double,double> raw_peaks;
		double wlval;
		g_abs->GetPoint(sample_273, wlval, raw_peaks.first);
		g_abs->GetPoint(sample_276, wlval, raw_peaks.second);
		std::cout<<"raw_peaks vals are: "<<raw_peaks.first<<", "<<raw_peaks.second<<std::endl;
		
		// if peaks are negative, we could coerce to 0...
		//raw_peaks.first=std::max(0.,raw_peaks.first);
		//raw_peaks.second=std::max(0.,raw_peaks.second);
		// but this tends to only happen for the first value, and results in a point at (0,0)
		// this point then seems to be way out of the trend, so probably throws the fit off.
		// instead, just skip this concentration
		if(raw_peaks.first<0 || raw_peaks.second<0){
			std::cerr<<"skipping concentration measurement "<<i<<" as one of the peaks is negative"<<std::endl;
			continue;
		}
		
		std::cout<<"ZZZ"<<std::endl;
		// just taking the data value at a specific wavelength may be a bit naff
		// as the sampling is sparse, and we may slightly miss the absorption peak.
		// We may do better by interpolating the peak to find a better estimate of maximum.
		// We can do this by fitting the peaks with gaussians, but since these are peaks
		// on a non-uniform background, we can only fit within a narrow region close to the peak.
		std::pair<double,double> simple_peaks, simple_errs, simple_posns;
		int ok = FitTwoGaussians(g_abs, simple_peaks, simple_errs, simple_posns);
		// total error from adding in quadrature
		double simplerr = sqrt(TMath::Sq(simple_errs.first)+TMath::Sq(simple_errs.second));
		
		// fit it with a combination of 4 gaussians
		TF1* abs_func = GetAbsFunc();
		// The fit has a tendency to screw up, so seed the initial values based on the raw fit.
		// These initial values will be over-estimates, since the complex fit peak heights also
		// have contributions from the shoulder gaussians, so peak1 component amplitude is less,
		abs_func->SetParameter("peak 1 amp",raw_peaks.first);
		abs_func->SetParameter("peak 2 amp",raw_peaks.second/raw_peaks.first);
		// also set the other components relative to this - no longer needed, par definitions are already relative
		//abs_func->SetParameter("RH shoulder amp",raw_peaks.first*0.5);
		//abs_func->SetParameter("LH shoulder amp",raw_peaks.second*0.2);
		g_abs->Fit(abs_func,"RQN");
		
		/*
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
		*/
		
		// extract the height of the peaks by finding the maximum of the curve
		// within a region around the peak central location. Note we CANNOT
		// just take the gaussian amplitude because there are multiple overlapping
		// contributions to the peak here. We could just evaluate the function at
		// the known peak position, but this allows more flexibility.
		double peak1_pos = abs_func->GetMaximumX(272.5,273.5);
		double peak2_pos = abs_func->GetMaximumX(275.,276.);
		// extract the peak heights by evaluating the function at that x
		double peak1_height  = abs_func->Eval(peak1_pos);
		double peak2_height  = abs_func->Eval(peak2_pos);
		// calculate error on peak heights
		std::pair<double,double> complexerrp = CalculateError(abs_func, peak1_pos, peak2_pos);
		double complexerr = sqrt(TMath::Sq(complexerrp.first)+TMath::Sq(complexerrp.second));
		
		// now look up the concentration
		// get name of the file this measuerement is from
		std::string current_file = c_led->GetTree()->GetCurrentFile()->GetName();
		std::cout<<"filename: "<<current_file<<std::endl;
		// strip off path
		if(current_file.find("/")!=std::string::npos){
			current_file = current_file.substr(current_file.find_last_of("/")+1,std::string::npos);
		}
		
		// look up its concentration (and error)
		if(calibration_data.count(current_file)==0){
			std::cerr<<"Couldn't find file "<<current_file<<" in calibration data map!"<<std::endl;
			for(auto&& am : calibration_data){
				std::cerr<<am.first<<", ";
			}
			std::cout<<std::endl;
			return nullptr;
		}
		double current_conc = calibration_data.at(current_file).first;
		double current_concerr = calibration_data.at(current_file).second;
		if(lastconc==current_conc) continue; /// FIXME HACK because i messed up, skip this file
		lastconc=current_conc;
		
		concentrations.push_back(current_conc);
		concentration_errs.push_back(current_concerr);
		
		// add to the calibration data curve
		calib_curve_raw->SetPoint(next_graph_point,current_conc,raw_peaks.first-raw_peaks.second);
		calib_curve_simple->SetPoint(next_graph_point,current_conc,simple_peaks.first-simple_peaks.second);
		calib_curve_complex->SetPoint(next_graph_point,current_conc,peak1_height-peak2_height);
		// FIXME calculate errors, error on peak 2 needs to account for error on peak 1?
		// error on raw result needs to come from error from spectrometer
		calib_curve_raw->SetPointError(next_graph_point,current_concerr,0.1); // XXX ??? TODO
		calib_curve_simple->SetPointError(next_graph_point,current_concerr, simplerr);
		calib_curve_complex->SetPointError(next_graph_point,current_concerr, complexerr);
		++next_graph_point;
		/*
		std::cout<<"raw diff: "<<raw_peaks.first<<" - "<<raw_peaks.second<<" = "<<raw_peaks.first-raw_peaks.second<<std::endl;
		std::cout<<"simple diff: "<<simple_peaks.first<<" - "<<simple_peaks.second<<" = "<<simple_peaks.first-simple_peaks.second<<std::endl;
		std::cout<<"complex diff: "<<peak1_height<<" - "<<peak2_height<<" = "<<peak1_height-peak2_height<<std::endl;
		std::cout<<"raw err: "<<0<<", simple err: "<<simplerr<<", complexerr: "<<complexerr<<std::endl;
		*/
		
		//delete mgg;
		delete g_inband;
		delete g_sideband;
		delete g_other;
		delete g_abs;
		delete abs_func;
		
	}
	
	std::cout<<"000"<<std::endl;
	
	// curtail the TGraphs in case we didn't actually set all the points due to bad fits
	calib_curve_raw->Set(next_graph_point);
	calib_curve_simple->Set(next_graph_point);
	calib_curve_complex->Set(next_graph_point);
	
	// Fit the calibration curves
	TF1* calib_func_raw = new TF1("calib_func_raw", "pol6", 0, 0.4);
	calib_curve_raw->Fit("calib_func_raw","RQ");
	
	TF1* calib_func_simple = new TF1("calib_func_simple", "pol6", 0, 0.4);
	calib_curve_simple->Fit("calib_func_simple","RQ");
	
	TF1* calib_func_complex = new TF1("calib_func_complex", "pol6", 0, 0.4);
	calib_curve_complex->Fit("calib_func_complex","RQ");
	
	calib_curve_simple->SetLineColor(kRed);
	calib_curve_simple->SetLineWidth(0);
	if(name=="275_A"){
		calib_curve_simple->SetMarkerStyle(20);
	} else if(name=="275_B"){
		calib_curve_simple->SetMarkerStyle(30);
	}
	calib_curve_simple->SetMarkerColor(kRed);
	calib_func_simple->SetLineWidth(1);
	calib_func_simple->SetLineColor(kRed);
	
	calib_curve_complex->SetLineWidth(0);
	calib_curve_complex->SetLineColor(kBlue);
	if(name=="275_A"){
		calib_curve_complex->SetMarkerStyle(20);
	} else if(name=="275_B"){
		calib_curve_complex->SetMarkerStyle(30);
	}
	calib_curve_complex->SetMarkerColor(kBlue);
	calib_func_complex->SetLineWidth(1);
	calib_func_complex->SetLineColor(kBlue);
	
	calib_curve_raw->SetLineColor(kSpring-5);
	calib_curve_raw->SetLineWidth(0);
	if(name=="275_A"){
		calib_curve_raw->SetMarkerStyle(20);
	} else if(name=="275_B"){
		calib_curve_raw->SetMarkerStyle(30);
	}
	calib_curve_raw->SetMarkerColor(kSpring-5);
	calib_func_raw->SetLineWidth(1);
	calib_func_raw->SetLineColor(kSpring-5);
	
	
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
	
	TCanvas c4("c4","c4",1024,800);
	std::vector<double> numberline(concentrations.size());
	std::iota(numberline.begin(), numberline.end(), 0);
	std::vector<double> zeros(concentrations.size()); // errors on measurement number: 0
	std::fill(zeros.begin(), zeros.end(), 0);
	TGraphErrors g4(concentrations.size(), numberline.data(), concentrations.data(), zeros.data(), concentration_errs.data());
	g4.Draw("ALP");
	
	//TCanvas c5("c5","c5",1024,800);
	// to try to compare shapes, scale to the same max value then plot on the same canvas
	TGraphErrors g5(concentrations.size(), numberline.data(), calib_curve_raw->GetY(), zeros.data(), calib_curve_raw->GetEY());
	g5.SetLineColor(kRed);
	g5.SetMarkerColor(kRed);
	g5.Draw("LP same");
	
	gSystem->ProcessEvents();
	while(gROOT->FindObject("cc")!=nullptr){
		gSystem->ProcessEvents();
		std::this_thread::sleep_for(std::chrono::milliseconds(100));
	}
	
	std::cout<<"returning multigraph"<<std::endl;
	return mg_calibs;
	
}


void Plotter::SetTimeAxis(TH1* hist, long t0_seconds){
	// start_time is in seconds: plotted times should be in seconds relative to this
	gStyle->SetTimeOffset(t0_seconds);
	// then you set the x-axis to time format
	hist->GetXaxis()->SetTimeDisplay(1);
	hist->GetXaxis()->SetTimeFormat("%m/%d/%y %H:%M");
	hist->GetXaxis()->SetLabelSize(0.03);
}


