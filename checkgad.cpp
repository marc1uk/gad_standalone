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
#include "TObjectTable.h"
#include "TMath.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TMinuit.h"
#include "TMatrixD.h"
#include "TGaxis.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <thread>
#include <cassert>
#include <sys/stat.h>  // dirname and basename
#include <sys/types.h> // for stat() test to see if file or folder


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
	std::map<std::string,int> offsets;
	int n_datapoints;
	int num_meas;
	int n_multigraph_lines = 200;
	int n_colours;
	int refentry=0;
	int darksperloop=12;
	
	std::map<std::string, std::pair<double,double>> calibration_data;
	
	int LoadData();
	int LoadCalibrationData();
	int LoadCalibrationData2();
	int LoadCalibrationData3();
	int LoadCalibrationData4();
	int LoadCalibrationData5();
	int LoadCalibrationData6();
	int LoadCalibrationData7();
	int SetBranchAddresses();
	int GetNextDarkEntry(std::string name, int ledon_entry_num, bool check=true);
	TMultiGraph* MakeMultigraph(std::string name, bool relative_to_first, bool darksub, char relop='/');
	TGraph* MakeTrendGraph(std::string name, bool relative_to_first, bool darksub, bool darkval, int trendsample);
	
	TGraphErrors* LoadPureFile(std::string purefilename);
	bool LoadConcentrations(std::string filename);
	TF1* PureScaledPlusLinear(TGraphErrors* dark_subtracted_pure);
	TF1* PureScaledPlusLinearv2();
	TF1* PureScaledPlusExtras();
	TMultiGraph* FitPurePlusLinear(std::string name);
	TMultiGraph* ExtractAbsorbance(std::string name, std::string calib_version="old");
	std::pair<double,double> CalculateError(TF1* abs_func, double peak1_pos, double peak2_pos);
	int FitTwoGaussians(TGraph* abs_graph, std::pair<double,double>& simple_peaks, std::pair<double,double>& simple_errs, std::pair<double,double>& simple_posns, bool plot=false);
	TF1* GetAbsFunc();
	TMultiGraph* FitCalibrationData(std::string name, std::string dataset);
	
	void SetTimeAxis(TH1* hist, long t0_seconds);
	bool CheckPath(std::string path, std::string& type);
	bool MakePure(std::string name, bool overwrite);
	
	int sample_273=0;
	int sample_276=0;
	
	//TF1* purefunc=nullptr;
	
	TGraph* trend_A=nullptr;
	TGraph* trend_B=nullptr;
	
	int Execute();
};
std::string purefile="../GDConcMeasure/pureDarkSubtracted.root";

int main(){
	
	TApplication myapp("rootTApp",0,0);
	
	Plotter myplotter;
	myplotter.Execute();
	
	TCanvas* c1 = (TCanvas*)gROOT->FindObject("c1");
	while(c1!=nullptr){
		gSystem->ProcessEvents();
		std::this_thread::sleep_for(std::chrono::milliseconds(100));
		c1 = (TCanvas*)gROOT->FindObject("c1");
	}
	
	return 0;
}

int Plotter::LoadData(){
	

	std::cout<<"adding run july data to chains"<<std::endl;
	for(auto&& c : chains){
		c.second->Add("../GDConcMeasure/data/2022/07/*08*");
		c_dark->Add("../GDConcMeasure/data/2022/07/*08*");
	}

	std::cout<<"adding run 66 april data to chains"<<std::endl;
	for(auto&& c : chains){
		c.second->Add("../GDConcMeasure/data/2022/04/66_*");
		c_dark->Add("../GDConcMeasure/data/2022/04/66_*");
	}
	std::cout<<"had "<<c_275_A->GetEntries()<<" april entries"<<std::endl;
	
	std::cout<<"adding run 66 may data to chains"<<std::endl;
	for(auto&& c : chains){
		c.second->Add("../GDConcMeasure/data/2022/05/66_*");
		c_dark->Add("../GDConcMeasure/data/2022/05/66_*");
	}
	std::cout<<"had "<<c_275_A->GetEntries()<<" april+may entries"<<std::endl;
	
	std::cout<<"adding second may dataset to led chains"<<std::endl;
	for(auto&& c : chains){
		c.second->Add("../GDConcMeasure/data/2022/05/stability_*");
	}
	
	// the Dark tree changed its name from 'dark' to 'Dark',
	// so we can't add entries from a TTree called 'Dark' to a chain
	// with name 'dark' with the normal method - it woudn't pick the tree up.
	// allegedly we can override the name of the tree in a file being added like this:
	//c_dark->Add("/home/pi/GDConcMeasure/data/2022/05/stability_*?#Dark");
	// but this doesn't seem to work when combined with wildcarding.
	// Instead we need to iterate over the files, and add them to the TChain
	// one-by-one, overriding the name of the TTree they contain.
	// One way we can build the list of matching files is with a dummy TChain.
	std::cout<<"adding second may dataset to dark chains"<<std::endl;
	TChain dummyc("");
	dummyc.Add("../GDConcMeasure/data/2022/05/stability_*");
	TIter nextIt(dummyc.GetListOfFiles());
	while(TObject* nextObj = nextIt()){
		std::string fname = nextObj->GetTitle();
		fname += "?#Dark";
		c_dark->Add(fname.c_str());
	}
	
	// get number of measurements
	num_meas = c_275_A->GetEntries();
	std::cout<<"had "<<num_meas<<" april+may+untagged entries"<<std::endl;
	
	return num_meas;
}

int Plotter::LoadCalibrationData(){
	
	for(auto&& c : chains){
		c.second->Add("../GDConcMeasure/Feb24thCalib/Feb24thCalib_*.root");
	}
	c_dark->Add("../GDConcMeasure/Feb24thCalib/Feb24thCalib_*.root");
	
	// get number of measurements
	num_meas = c_275_A->GetEntries();
	std::cout<<"had "<<num_meas<<" april+may+untagged entries"<<std::endl;
	
	return num_meas;
}

int Plotter::LoadCalibrationData2(){
	
	for(auto&& c : chains){
		c.second->Add("../GDConcMeasure/data/2022/05/[0-9]*26May22Calib_[0-9]*");
	}
	c_dark->Add("../GDConcMeasure/data/2022/05/[0-9]*26May22Calib_[0-9]*");
	
	// get number of measurements
	num_meas = c_275_A->GetEntries();
	std::cout<<"had "<<num_meas<<" may calibration entries"<<std::endl;
	
	return num_meas;
}

int Plotter::LoadCalibrationData3(){
	
	for(auto&& c : chains){
		c.second->Add("../GDConcMeasure/data/2022/06/rename/[0-9]*01Jun22Calib_*");
	}
	c_dark->Add("../GDConcMeasure/data/2022/06/rename/[0-9]*01Jun22Calib_*");
	
	// get number of measurements
	num_meas = c_275_A->GetEntries();
	std::cout<<"had "<<num_meas<<" jun calibration entries"<<std::endl;
	
	return num_meas;
}

int Plotter::LoadCalibrationData4(){
	
	for(auto&& c : chains){
		c.second->Add("../GDConcMeasure/data/2022/06/jun06calib/*");
		//c.second->Add("../GDConcMeasure/data/2022/06/[0-9]*06Jun22Calib_*");
	}
	c_dark->Add("../GDConcMeasure/data/2022/06/jun06calib/*");
	//c_dark->Add("../GDConcMeasure/data/2022/06/[0-9]*06Jun22Calib_*");
	
	// get number of measurements
	num_meas = c_275_A->GetEntries();
	std::cout<<"had "<<num_meas<<" jun calibration entries"<<std::endl;
	
	return num_meas;
}

int Plotter::LoadCalibrationData5(){
	
	for(auto&& c : chains){
		c.second->Add("../GDConcMeasure/data/2022/06/[0-9]*08Jun22Calib_*");
	}
	c_dark->Add("../GDConcMeasure/data/2022/06/[0-9]*08Jun22Calib_*");
	
	// get number of measurements
	num_meas = c_275_A->GetEntries();
	std::cout<<"had "<<num_meas<<" jun calibration entries"<<std::endl;
	
	return num_meas;
}

int Plotter::LoadCalibrationData6(){
	
	for(auto&& c : chains){
		c.second->Add("../GDConcMeasure/data/2022/07/[0-9]*08July22Calib_*");
	}
	c_dark->Add("../GDConcMeasure/data/2022/07/[0-9]*08July22Calib_*");
	
	// get number of measurements
	num_meas = c_275_A->GetEntries();
	std::cout<<"had "<<num_meas<<" july calibration entries"<<std::endl;
	
	return num_meas;
}

int Plotter::LoadCalibrationData7(){
	
	for(auto&& c : chains){
		c.second->Add("../GDConcMeasure/data/2022/07/[0-9]*12July22Calib_*");
	}
	c_dark->Add("../GDConcMeasure/data/2022/07/[0-9]*12July22Calib_*");
	
	// get number of measurements
	num_meas = c_275_A->GetEntries();
	std::cout<<"had "<<num_meas<<" july calibration entries"<<std::endl;
	
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
	
	std::cout<<"getting c_275_A entry 0"<<std::endl;
	int nbytes = c_275_A->GetEntry(0);
	if(nbytes<=0){
		std::cerr<<"failure getting entry 0!"<<std::endl;
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
	
	return 1;
}

int Plotter::Execute(){
	
	std::cout<<"making chains"<<std::endl;
	std::string dataset="julcal2";
	if(dataset=="febstab"){
		// data
		c_dark = new TChain("dark");
		c_275_A = new TChain("275_A");
		c_275_B = new TChain("275_B");
		c_R_G_B = new TChain("R_G_B");
		c_White_385 = new TChain("White_385");
	} else if(dataset=="febcal"){
		// calibration
		c_dark = new TChain("dark");
		c_275_A = new TChain("260");
		c_275_B = new TChain("275");
		c_R_G_B = new TChain("R_G_B");
		c_White_385 = new TChain("White_385");
	} else if(dataset=="maycal" || dataset=="juncal" || dataset=="juncal2" || dataset=="juncal3" || dataset=="julcal" || dataset=="julcal2" ){
		// calibration
		c_dark = new TChain("Dark");
		c_275_A = new TChain("275_A");
		c_275_B = new TChain("275_B");
		c_R_G_B = new TChain("R_G_B");
		c_White_385 = new TChain("White_385");
	}
	
	chains = std::map<std::string, TChain*>{
		{"275_A", c_275_A},
		{"275_B", c_275_B},
		{"R_G_B", c_R_G_B},
		{"White_385", c_White_385}
	};
	
	if(dataset=="febstab"){
		LoadData();
	} else if(dataset=="febcal"){
		LoadCalibrationData();
	} else if(dataset=="maycal"){
		LoadCalibrationData2();
	} else if(dataset=="juncal"){
		LoadCalibrationData3();
	} else if(dataset=="juncal2"){
		LoadCalibrationData4();
	} else if(dataset=="juncal3"){
		LoadCalibrationData5();
	} else if(dataset=="julcal"){
		LoadCalibrationData6();
	} else if(dataset=="julcal2"){
		LoadCalibrationData7();
	}
	SetBranchAddresses();
	
	gStyle->SetPalette(kBird);
	n_colours = TColor::GetNumberOfColors();
	
	TCanvas* c1 = new TCanvas("c1","c1",1024,800);  // do not change the name of this canvas, closing it will terminate the program.
	
	/* based on the measurement command list:
		start
		power on
		wait 900
		pump on
		start_loop
		valve outlet open
		valve inlet open
		measure Dark
		measure Dark
		measure Dark
		measure Dark
		valve inlet close
		valve outlet close
		wait 150
		measure Dark
		measure Dark
		measure Dark
		measure Dark
		measure Dark
		measure 275_B
		measure Dark
		measure 275_A
		measure Dark
		measure R G B
		measure Dark
		measure White 385
		save trash
		wait 600
		loop
		quit
	we take 9 darks before the first 275_B measurement,
	so the corressponding dark for this trace would be dark entry 8
	similarly for the 275_B measurement the corresponding dark entry is 9
	10 for R_G_B and 11 for White_385
	we take 12 darks in total each loop, so add 12 to each index per loop
	*/
	offsets = std::map<std::string,int>{
		{"275_B",8},
		{"275_A",9},
		{"R_G_B",10},
		{"White_385",11}
	};
	
	std::vector<double> numberline(num_meas);
	std::iota(numberline.begin(),numberline.end(),0);
	
	/*
	// plot all traces over time, coloured by time
	// MakeMultigraph(std::string name, bool relative_to_first, bool darksub)
	TMultiGraph* mg_275_A = MakeMultigraph("275_B", false, true);
	//TCanvas* c1 = new TCanvas("c1","c1",1024,800);  // do not change the name of this canvas, closing it will terminate the program.
	c1->cd();
	mg_275_A->Draw("AL");
	c1->Modified();
	c1->Update();
	
	// plot all traces over time, relative to the first measurement
	TMultiGraph* mg_275_A_ref = MakeMultigraph("275_B", true, true);
	TCanvas* c2 = new TCanvas("c2","c2",1024,800);
	mg_275_A_ref->Draw("AL");
	c2->Modified();
	c2->Update();
	
	// plot trend of the 275_A trace at 273nm over time
	// MakeTrendGraph(std::string name, bool relative_to_first, bool darksub, bool darkval, int trendsample)
	TGraph* trend_273_rel = MakeTrendGraph("275_B", true, true, false, sample_273);
	trend_273_rel->SetLineColor(kRed);
	// same trend at 276nm over time
	TGraph* trend_276_rel = MakeTrendGraph("275_B", true, true, false, sample_276);
	trend_276_rel->SetLineColor(kBlue);
	// plot the two
	TMultiGraph* mgtrend = new TMultiGraph("273_276_trend","273_276_trend");
	mgtrend->Add(trend_273_rel);
	mgtrend->Add(trend_276_rel);
	TCanvas* c3 = new TCanvas("c3","c3",1024,800);
	mgtrend->Draw("AL");
	c3->Modified();
	c3->Update();
	
	// plot the difference between ratios at 273 vs 276 over time:
	// this shows whether the relative amplitude at 273 changes vs the relative amplitude at 276
	// after correcting for LED power - i.e., whether the shape of the spectrum changes
	//TGraph* trend_273 = MakeTrendGraph("275_A", false, true, false, sample_273);
	//TGraph* trend_276 = MakeTrendGraph("275_A", false, true, false, sample_276);
	std::vector<double> peak_diffs(num_meas);
	for(int i=0; i<num_meas; ++i){
		double xval_a, yval_a, xval_b, yval_b;
		trend_273_rel->GetPoint(i,xval_a, yval_a);
		trend_276_rel->GetPoint(i,xval_b, yval_b);
		peak_diffs.at(i) = yval_b - yval_a;
	}
	TGraph* diff_trend = new TGraph(num_meas, numberline.data(), peak_diffs.data());
	diff_trend->SetTitle("diff_trend");
	diff_trend->SetName("diff_trend");
	TCanvas* c4 = new TCanvas("c4","c4",1024,800);
	diff_trend->Draw("AL");
	c4->Modified();
	c4->Update();
	
	// plot the trend in dark data over time. this may be a proxy for spectrometer temperature, possibly.
	// plot it for both LEDs, to see if their trends are in sync, suggesting temperatures in sync
	// MakeTrendGraph(std::string name, bool relative_to_first, bool darksub, bool darkval, int trendsample)
	TGraph* trend_dark_A_273 = MakeTrendGraph("275_A", false, false, true, sample_273);
	TGraph* trend_dark_A_276 = MakeTrendGraph("275_A", false, false, true, sample_276);
	TGraph* trend_dark_B_273 = MakeTrendGraph("275_B", false, false, true, sample_273);
	trend_dark_A_273->SetLineColor(kRed);
	trend_dark_A_276->SetLineColor(kMagenta);
	trend_dark_B_273->SetLineColor(kBlue);
	TMultiGraph* mg_darks = new TMultiGraph("mg_darks","mg_darks");
	mg_darks->Add(trend_dark_A_273);
	mg_darks->Add(trend_dark_A_276);
	mg_darks->Add(trend_dark_B_273);
	TCanvas* c5 = new TCanvas("c5","c5",1024,800);
	mg_darks->Draw("AL");
	c5->Modified();
	c5->Update();
	
	std::vector<double> darkdiff273_276(num_meas);
	for(int i=0; i<num_meas; ++i){
		double tmpx273, tmpy273, tmpx276, tmpy276;
		trend_dark_A_273->GetPoint(i,tmpx273, tmpy273);
		trend_dark_A_276->GetPoint(i,tmpx276, tmpy276);
		darkdiff273_276.at(i) = tmpy273 - tmpy276;
	}
	TGraph* darkdiff273_276g = new TGraph(num_meas, numberline.data(), darkdiff273_276.data());
	TCanvas* c6 = new TCanvas("c6","c6",1024,800);
	darkdiff273_276g->Draw("AL");
	c6->Modified();
	c6->Update();
	
	// plot the dark value (temperature) vs peak difference (concentration)
	TGraph* gversus = new TGraph(num_meas, trend_dark_A_273->GetY(), peak_diffs.data());
	gversus->SetTitle("[dark rate] vs [LHS-RHS peak variance diff]");
	gversus->SetName("[dark rate] vs [LHS-RHS peak variance diff]");
	TCanvas* c7 = new TCanvas("c7","c7",1024,800);
	gversus->Draw("AP*");
	c7->Modified();
	c7->Update();
	
	///////////////////////////////////
	// plot the variation in amplitude ratio between the first and last peak against wavelength
	//MakeMultigraph(std::string name, bool relative_to_first, bool darksub, char relop)
	TMultiGraph* mg_275_A_dif = MakeMultigraph("275_A", true, true,'/');
	TMultiGraph* aaa = new TMultiGraph("aaa","aaa");
	TF1 afunc("afunc","[0] + [1]*x",270,330);
	std::vector<double> grads(num_meas);
	for(int i=0; i<num_meas; ++i){
		// we actually only want the last tgraph
		TGraph* diffgraph = (TGraph*)mg_275_A_dif->GetListOfGraphs()->At(i); // num_meas-1
		// as a reference point we'll use 273
		double* diffvals = diffgraph->GetY();
		double diffat273 = diffvals[sample_273];
		std::vector<double> diffratios(n_datapoints);
		for(int i=0; i<n_datapoints; ++i){
			diffratios.at(i) = diffvals[i] / diffat273;
		}
		TGraph* gbackg = new TGraph(n_datapoints, wavelengths.data(), diffratios.data());
		int colour_i = int(double(i)*(double(n_colours)/double(num_meas)));
		gbackg->SetLineColor(TColor::GetColorPalette(colour_i));
		aaa->Add(gbackg);
		gbackg->Fit("afunc","R");
		double grad = afunc.GetParameter(1);
		gbackg->GetListOfFunctions()->Clear();
		grads.at(i) = grad;
	}
	TCanvas* c8 = new TCanvas("c8","c8",1024,800);
	//gbackg->Draw("AL");
	TGraph* ggrads = new TGraph(num_meas, numberline.data(), grads.data());
	ggrads->Draw("AL");
	c8->Modified();
	c8->Update();
	
	TCanvas* c9 = new TCanvas("c9","c9",1024,800);
	aaa->Draw("AL");
	c9->Modified();
	c9->Update();
	
	TMultiGraph* mgresid = FitPurePlusLinear("275_A");
	mgresid->Draw("AL");
	*/
	
	// step 1: analyse calibration data for polynomial coffiecients,
	// take the results printed out and populate the coefficients list
	// in ExtractAbsorbance
	TMultiGraph* mg_cals = FitCalibrationData("275_B", dataset);  // run for both 275_A and 275_B
	if(mg_cals==nullptr) return 0;
	mg_cals->Draw("AL");
	
	/*
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
	TGraph* calgraph = new TGraph(calib_concs.size(), numberline.data(), calib_concs.data());
	calgraph->SetLineColor(kMagenta);
	calgraph->SetMarkerColor(kMagenta);
	calgraph->SetMarkerStyle(30);
	mg_all->Add(calgraph);
	
	c1->cd();
	mg_all->Draw("ALP");
	*/
	
	
	return 0;
}

TGraphErrors* Plotter::LoadPureFile(std::string purefilename){
	// how the heck does ownership and scope happen here?
	// pure_file is local, presumably dark_subtracted_pure is therefore also local
	// which would mean purefunc is calling Eval on a TGraphErrors that goes out of scope?
	// As PepeLePew says: https://root-forum.cern.ch/t/object-ownership-by-tfile/13772
	// the TGraphErrors does NOT go out of scope, so presumably we are now responsible for it...
	
	//std::cout<<"TObjectTable before loading pure file:\n";
	//gObjectTable->Print();
	
	TFile pure_file(purefilename.c_str(), "READ");
	if (! pure_file.IsOpen()){
		std::cerr<<"failed to open "<<purefilename<<std::endl;
		return nullptr;
	}
	
	// retrieve the pure water trace (LED spectrum) as a TGraph
	TGraphErrors* dark_subtracted_pure = (TGraphErrors*) pure_file.Get("Graph");
	if(dark_subtracted_pure==nullptr){
		std::cerr<<"failed to get 'Graph' from pure file "<<purefilename<<std::endl;
		return nullptr;
	}
	
	// minor optimization; if the data in the TGraph is sorted by x value
	//(which it should be for us), then setting the following option can speed up Eval() calls
	dark_subtracted_pure->SetBit(TGraph::kIsSortedX); // XXX note only available on newer ROOT
	
	//std::cout<<"TObjectTable before returning from loadpurefile:\n";
	//gObjectTable->Print();
	
	return dark_subtracted_pure;
}

TF1* Plotter::PureScaledPlusLinear(TGraphErrors* dark_subtracted_pure){
	
	// construct functional fit. We'll scale and add a linear background.
	// limit the pure function to a region in which we have light - no point fitting outside this region
	const int wave_min = 260, wave_max = 300, numb_of_fitting_parameters = 3;
	TF1* pure = new TF1("pure_fct",
	                       // TF1 wraps a c++ lambda function that evaluates
	                       // the reference TGraph at a given wavelength
	                      [dark_subtracted_pure](double* x, double* par){
	                      return (par[0] * dark_subtracted_pure->Eval(x[0])) + (par[1] * x[0]) + par[2];},
	                      wave_min, wave_max, numb_of_fitting_parameters);
	
	return pure;
	
}

TGraphErrors* dark_subtracted_pure =nullptr;
// my local installation (v6.06) does not like using lambdas - use a global function instead
double PureFunc(double* x, double* par){
	if(dark_subtracted_pure==nullptr){
		TFile *_file0 = TFile::Open(purefile.c_str());
		dark_subtracted_pure = (TGraphErrors*)_file0->Get("Graph");
	}
	double purepart = (par[0] * dark_subtracted_pure->Eval(x[0]));
	double linpart = (par[1] * x[0]) + par[2];
	double retval = purepart + linpart;
	return retval;
}

TF1* Plotter::PureScaledPlusLinearv2(){
	
	// construct functional fit. We'll scale and add a linear background.
	// limit the pure function to a region in which we have light - no point fitting outside this region
	static int purever=0;
	purever++;
	std::string name="pure_fct"+std::to_string(purever);
	const int wave_min = 260, wave_max = 300, numb_of_fitting_parameters = 3;
	TF1* pure = new TF1(name.c_str(), PureFunc, wave_min, wave_max, numb_of_fitting_parameters);
	// set default parameters
	pure->SetParameters(1.,0.,0.);
	
	return pure;
	
}

double PureFuncv2(double* x, double* par){
	if(dark_subtracted_pure==nullptr){
		TFile *_file0 = TFile::Open(purefile.c_str());
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

TMultiGraph* Plotter::FitPurePlusLinear(std::string name){
	
	/*
	std::cout<<"loading pure file and getting TGraphErrors"<<std::endl;
	TGraphErrors* puretgraph = LoadPureFile(purefile);
	
	//std::cout<<"TObjectTable having returned from loadpurefile:\n";
	//gObjectTable->Print();
	
	std::cout<<"making pure plus linear TF1"<<std::endl;
	TF1* purepluslinear = PureScaledPlusLinear(puretgraph);
	*/
	
	TF1* purepluslinear = PureScaledPlusExtras();
//	TF1* purepluslinear = PureScaledPlusLinearv2();
	std::cout<<"makde tf1"<<std::endl;
	
	TChain* c_led = chains.at(name);
	std::cout<<"led chain "<<name<<" is at "<<c_led<<std::endl;
	TChain* c_dark = chains.at("dark");
	int offset = offsets.at(name);
	int n_entries = c_led->GetEntries();
	std::cout<<"led chain has "<<n_entries<<" entries"<<std::endl;
	
	std::string title = name+"_g";
	TMultiGraph* mg = new TMultiGraph(title.c_str(), title.c_str());
	
	TCanvas cnew("cnew","cnew",1024,800);
	
	std::cout<<"looping over "<<n_entries<<" entries"<<std::endl;
	std::vector<int> interesting_entries{0,1};
	for(int i=2; i>0; --i) interesting_entries.push_back(n_entries-i);
	
//	for(int& i : interesting_entries){
	for(int i=0; i<n_entries; i++){
		
		// get LED on trace
		std::cout<<"getting entry "<<i<<std::endl;
		int nbytes = c_led->GetEntry(i);
		if(nbytes<=0){
			std::cerr<<"COULDN'T LOAD ENTRY "<<i<<" FROM LED CHAIN"<<std::endl;
			exit(-1);
		}
		
		// get nearest dark trace
		std::cout<<"looking up matching dark entry"<<std::endl;
		int dark_entry = offset + i*darksperloop;
		//int dark_entry = GetNextDarkEntry(name, i);
		std::cout<<"getting dark entry "<<dark_entry<<std::endl;
		c_dark->GetEntry(dark_entry);
		
		// subtract dark
		std::cout<<"doing dark subtraction"<<std::endl;
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
		}
		
		// make TGraphErrors from data
		std::cout<<"making tgrapherrors"<<std::endl;
		TGraph* g_inband = new TGraph(inband_values.size(), inband_wls.data(), inband_values.data());
		TGraph* g_sideband = new TGraph(sideband_values.size(), sideband_wls.data(), sideband_values.data());
		TGraph* g_other = new TGraph(other_values.size(), other_wls.data(), other_values.data());
		
		std::string thistitle = title+"_"+std::to_string(i);
		std::string intitle = thistitle+"_in";
		std::string sidetitle = thistitle+"_side";
		std::string othertitle = thistitle+"_other";
		g_inband->SetTitle(thistitle.c_str());
		g_inband->SetName(thistitle.c_str());
		g_inband->SetLineWidth(1);
		g_inband->SetMarkerStyle(7);
		g_inband->SetMarkerColor(kSpring-5);  // kind of lightish green
		g_sideband->SetTitle(sidetitle.c_str());
		g_sideband->SetName(sidetitle.c_str());
		g_sideband->SetLineWidth(0);
		g_sideband->SetMarkerStyle(7);
		g_sideband->SetMarkerColor(kRed);
		g_other->SetTitle(othertitle.c_str());
		g_other->SetName(othertitle.c_str());
		g_other->SetLineWidth(0);
		g_other->SetMarkerStyle(7);
		g_other->SetMarkerColor(kMagenta);
		
		// this has to come before the fit for it to register
		purepluslinear->SetLineColor(kRed);
		purepluslinear->SetLineWidth(4);
		
		// fit with pure scaled
		std::cout<<"fitting with purepluslinear"<<std::endl;
		g_sideband->Fit(purepluslinear,"RNQ");
		
		// split into components
		std::cout<<"cloning"<<std::endl;
		std::string purename = name+"_pure_"+std::to_string(i);
		TF1* purepart = PureScaledPlusExtras();
//		TF1* purepart = PureScaledPlusLinearv2();
		purepart->SetParameter(0,purepluslinear->GetParameter(0));
		purepart->SetParameter(1,0.);
		purepart->SetParameter(2,0.);
		purepart->SetLineColor(kBlue);
		purepart->SetLineWidth(2);
		
		std::cout<<"making linear component"<<std::endl;
		std::string linename = name+"_lin_"+std::to_string(i);
		TF1 linearpart(linename.c_str(),"[0]+[1]*x",260,300);
		linearpart.SetParameter(0,purepluslinear->GetParameter(2));
		linearpart.SetParameter(1,purepluslinear->GetParameter(1));
		linearpart.SetLineColor(kGreen-1);
		
		std::cout<<"drawing"<<std::endl;
		if(i>0) cnew.Clear();
		//TCanvas cnew2("cnew2","cnew2",1024,800);  // XXX uncomment to wait until user closes
		
		std::string mgtitle = "mgg_"+std::to_string(i);
		TMultiGraph* mgg = new TMultiGraph(mgtitle.c_str(),mgtitle.c_str());
		mgg->Add(g_inband);
		mgg->Add(g_sideband);
		mgg->Add(g_other);
		mgg->Draw("ALP");
		// these have to go *after* a draw or the graph has no axes
		mgg->GetYaxis()->SetRangeUser(-100,2000);
		mgg->GetXaxis()->SetRangeUser(240,320);
		purepluslinear->Draw("same");
//		purepart->Draw("same");
//		linearpart.Draw("same");
		cnew.Modified();
		cnew.Update();
//		cnew2.Modified();
//		cnew2.Update();
		gSystem->ProcessEvents();
		
		/*
		// XXX uncomment cnew2 to use
		std::cout<<"waiting for user to close canvas"<<std::endl;
		while(gROOT->FindObject("cnew2")!=nullptr){
			gSystem->ProcessEvents();
			std::this_thread::sleep_for(std::chrono::milliseconds(100));
		}
		*/
		std::this_thread::sleep_for(std::chrono::milliseconds(100));
		
		// make a plot of the fit residuals
		TGraph* g_resid = new TGraph(sideband_values.size());
		for(int k=0; k<sideband_values.size(); ++k){
			double wlval, dataval;
			g_sideband->GetPoint(k, wlval, dataval);
			double fitval = purepluslinear->Eval(wlval);
			g_resid->SetPoint(k,wlval,fitval-dataval);
		}
		std::string restitle = "g_resid_"+std::to_string(i);
		g_resid->SetTitle(restitle.c_str());
		g_resid->SetName(restitle.c_str());
		// coerce to the nearest colour index in the palette
		int colour_i = int(double(i)*(double(n_colours)/double(n_entries)));
		g_resid->SetLineColor(TColor::GetColorPalette(colour_i));
		
		std::cout<<"adding to multigraph"<<std::endl;
		mg->Add(g_resid);
		
		std::cout<<"cleaning up temp functions"<<std::endl;
		delete purepart;
		delete mgg;
		
	}
	
	std::cout<<"returning multigraph"<<std::endl;
	return mg;
	
}

/////////////////////////////////

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
		
		std::cout<<"complex peak fit gaus "<<i<<" has ampltiude "
		         <<agaus.GetParameter(0)<<"+-"<<error_centre<<" scaled by "
		         <<relative_amp<<" to "<<error_here<<" at peak 1";
		
		totalerror1+= TMath::Sq(error_here);
		
		// repeat for contribution to peak2
		amp_here = agaus.Eval(peak2_pos);
		relative_amp = amp_here/amp_centre;
		// scale the error on the amplitude at the centre
		error_here = error_centre*relative_amp;
		
		std::cout<<" and scaled by "<<relative_amp<<" to "<<error_here<<" at peak 2"<<std::endl;
		
		totalerror2+= TMath::Sq(error_here);
	}
	
	// take sqrt of totalerror1 and totalerror2 each to get quadrature sum of contributions...
	// but then take square to add the errors on the two peak amplitudes in quadrature
	std::cout<<"total error at peak 1 "<<sqrt(totalerror1)<<" and at peak 2 "<<sqrt(totalerror2)<<std::endl;
	
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
	
	if(plot){
		TCanvas ctemp("ctemp","ctemp",1024,800);
		
		//std::cout<<"peak 1 amp: "<<gausamp1<<", gaus 2 amp: "<<gausamp2<<std::endl;
		abs_graph->Draw("ALP");
		gaus1.Draw("same");
		gaus2.Draw("same");
		ctemp.Modified();
		ctemp.Update();
		
		/*
		gfptr1->NormalizeErrors();
		gfptr1->Print();
		// print correlation matrix
		gMinuit->mnmatu(1);  /// ??? 2x2 matrix, what about 3x3?
		// print covariance matrix
		TMatrixD matrix0(3,3);
		gMinuit->mnemat(matrix0.GetMatrixArray(),3);
		matrix0.Print();
		
		// draw a 1-sigma error plot
		gMinuit->SetErrorDef(1); // 1-sigma. Argument is N² for N-sigma.
		TGraph *gr0 = (TGraph *)gMinuit->Contour(30,0,2);  // 30 points, param 0 vs 2
		//TGraph *gr0=nullptr;
		TCanvas ctemp2("ctemp2","ctemp2",1024,800);
		if(gr0){
			gr0->SetLineColor(kRed);
			gr0->Draw("alp");
			gr0->GetXaxis()->SetTitle("parameter 0 (amp)");
			gr0->GetYaxis()->SetTitle("parameter 2 (width)");
			gr0->SetTitle("1-sigma uncertainties on fit parameters");
			ctemp2.Modified();
			ctemp2.Update();
		}
		*/
		
		gSystem->ProcessEvents();
		while(gROOT->FindObject("ctemp")!=nullptr){
			gSystem->ProcessEvents();
			std::this_thread::sleep_for(std::chrono::milliseconds(100));
		}
	}
	
	
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
			
			if( false /*frptr->IsEmpty() || !frptr->IsValid() || frptr->Status()!=0*/){
				// fit failed; skip value?
				continue;
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
		
		// coerce to the nearest colour index in the palette
		//int colour_i = int(double(i)*(double(n_colours)/double(n_entries)));
		//g_abs->SetLineColor(TColor::GetColorPalette(colour_i));
		//mg->Add(g_abs);   // use this or below
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
	if(dataset=="febcal"){
		bool ok = LoadConcentrations("feb24_calibration_info.txt");
	} else if(dataset=="maycal"){
		bool ok = LoadConcentrations("may26_calibration_info_fix2.txt");
	} else if(dataset=="juncal"){
		bool ok = LoadConcentrations("jun01_calibration_info.txt"); //n.b. theoretical concs FIXME
	} else if(dataset=="juncal2"){
		bool ok = LoadConcentrations("jun06_calibration_info.txt");
	} else if(dataset=="juncal3"){
		bool ok = LoadConcentrations("jun08_calibration_info.txt");
	} else if(dataset=="julcal"){
		bool ok = LoadConcentrations("jul08_calibration_info.txt");
	} else if(dataset=="julcal2"){
		bool ok = LoadConcentrations("jul12_calibration_info.txt");
	}
	
	purefile = std::string("../GDConcMeasure/pureDarkSubtracted_") + name + "_" + dataset + ".root";
	bool ok = MakePure(name, false);
	if(not ok){
		return nullptr;
	}
	
	TF1* pureplusextras = PureScaledPlusExtras();
	
	TChain* c_led = chains.at(name);
	TChain* c_dark = chains.at("dark");
	int offset = offsets.at(name);
	int n_entries = c_led->GetEntries();
	n_entries = 53;
	
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
	std::vector<double> concentrations, concentration_errs, measnums;
	
	int next_graph_point=0;
	double lastconc=-99;
	for(int i=0; i<n_entries; i++){
		
		std::cout<<i<<": ";
		
		// get LED on trace
		int nbytes = c_led->GetEntry(i);
		if(nbytes<=0){
			std::cerr<<"COULDN'T LOAD ENTRY "<<i<<" FROM LED CHAIN"<<std::endl;
			exit(-1);
		}
		
		// get nearest dark trace
		//int dark_entry = offset + i*darksperloop;
		int dark_entry = GetNextDarkEntry(name, i, false);
		c_dark->GetEntry(dark_entry);
		
		// subtract dark
		for(int j=0; j<n_datapoints; ++j){
			values.at(j) -= values_dark.at(j);
		}
		
		// look up the concentration of this measurement.
		// first get name of the file this measuerement is from
		std::string current_file = c_led->GetTree()->GetCurrentFile()->GetName();
		std::cout<<"filename: "<<current_file<<std::endl;
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
			break;
			return nullptr;
		}
		double current_conc = calibration_data.at(current_file).first;
		double current_concerr = calibration_data.at(current_file).second;
		//if(lastconc==current_conc) continue; /// hack because i messed up, skip this file
		//lastconc=current_conc;
		
		concentrations.push_back(current_conc);
		concentration_errs.push_back(current_concerr);
		
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
		
		// calculate absorbance from ratio of fit to data in absorption region
		TGraphErrors* g_abs = new TGraphErrors(inband_values.size());
		std::string restitle = "g_abs_"+std::to_string(i);
		for(int k=0; k<inband_values.size(); ++k){
			double wlval, dataval;
			g_inband->GetPoint(k, wlval, dataval);
			double fitval = pureplusextras->Eval(wlval);
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
			/*
			std::cout<<"data val: "<<inband_values.at(k)<<"+- "<<error_on_data<<"\n"
			         <<", fractional error: "<<(error_on_data/dataval)<<"\n"
			         <<", fit val: "<<fitval<<", +- "<<error_on_fitval<<"\n"
			         <<", fractional error: "<<(error_on_fitval/fitval)<<"\n"
			         <<", total error on ratio: "<<err_on_ratio
			         <<", ratio of abs / transmitted: "<<(fitval/dataval)
			         <<", ratio of error to value: "<<err_on_ratio_over_ratio<<std::endl;
			*/
			g_abs->SetPointError(k, wavelength_errors.at(0)/2., err_on_ratio_over_ratio*(1./log(10.)));
		}
		
		// as the simplest estimate of peak heights, just take the graph value at the peaks.
		std::pair<double,double> raw_peaks;
		double wlval;
		g_abs->GetPoint(sample_273, wlval, raw_peaks.first);
		g_abs->GetPoint(sample_276, wlval, raw_peaks.second);
		
		// if peaks are negative, we could coerce to 0...
		// but this tends to only happen for the first value, and results in a point at (0,0)
		// this point then seems to be way out of the trend, so probably throws the fit off.
		// instead, just skip this concentration
		if(raw_peaks.first<0 || raw_peaks.second<0){
			std::cerr<<"skipping concentration measurement "<<i<<" as one of the peaks is negative"<<std::endl;
			raw_peaks.first=std::max(0.,raw_peaks.first);
			raw_peaks.second=std::max(0.,raw_peaks.second);
			continue;
		}
		
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
		std::cout<<"peak1_pos="<<peak1_pos<<", 2 pos="<<peak2_pos<<std::endl;
		// extract the peak heights by evaluating the function at that x
		double peak1_height  = abs_func->Eval(peak1_pos);
		double peak2_height  = abs_func->Eval(peak2_pos);
		if(TMath::IsNaN(peak1_height) || TMath::IsNaN(peak2_height)){
			peak1_height=0;
			peak2_height=0;
		}
		
		// calculate error on peak heights
		std::pair<double,double> complexerrp = CalculateError(abs_func, peak1_pos, peak2_pos);
		double complexerr = sqrt(TMath::Sq(complexerrp.first)+TMath::Sq(complexerrp.second));
		
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
		
		// extract measurement number from filename since they're not contiguous in a TChain
		// unless the filenames are alphabetical
		// 00113_06Jun22Calib_99.root
		std::string thismeasnum = current_file.substr(current_file.find_last_of("_")+1,current_file.length()-24);
		int this_meas_num = std::stoi(thismeasnum);
		std::cout<<"file "<<current_file<<" is measurement "<<this_meas_num<<std::endl;
		measnums.push_back(this_meas_num);
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
	
	// curtail the TGraphs in case we didn't actually set all the points due to bad fits
	calib_curve_raw->Set(next_graph_point);
	calib_curve_simple->Set(next_graph_point);
	calib_curve_complex->Set(next_graph_point);
	
	// Fit the calibration curves
	TF1* calib_func_raw = new TF1("calib_func_raw", "pol6", 0, 0.4);
	calib_curve_raw->Fit("calib_func_raw","RQN");
	
	TF1* calib_func_simple = new TF1("calib_func_simple", "pol6", 0, 0.4);
	calib_curve_simple->Fit("calib_func_simple","RQN");
	
	TF1* calib_func_complex = new TF1("calib_func_complex", "pol6", 0, 0.4);
	calib_curve_complex->Fit("calib_func_complex","RQN");
	
	// get fit parameters
	int npar=calib_func_simple->GetNpar();
	std::cout<<"npar: "<<npar<<std::endl;
	std::vector<double> calib_coefficients_raw(calib_func_raw->GetParameters(),calib_func_raw->GetParameters()+npar);
	std::vector<double> calib_coefficients_simple(calib_func_simple->GetParameters(),calib_func_simple->GetParameters()+npar);
	std::vector<double> calib_coefficients_complex(calib_func_complex->GetParameters(),calib_func_complex->GetParameters()+npar);
	for(int i=0; i<npar; ++i){
		std::cout<<"calibration coefficient "<<i<<":\n"
		         <<"raw: "<<calib_coefficients_raw.at(i)
		         <<", simple: "<<calib_coefficients_simple.at(i)
		         <<", complex: "<<calib_coefficients_complex.at(i)<<std::endl;
	}
	std::cout<<"raw fit chi2/NDF = "<<calib_func_raw->GetChisquare()<<"/"<<calib_func_raw->GetNDF()
	         <<" = "<<(calib_func_raw->GetChisquare()/calib_func_raw->GetNDF())<<std::endl;
	std::cout<<"simple fit chi2/NDF = "<<calib_func_simple->GetChisquare()<<"/"<<calib_func_simple->GetNDF()
	         <<" = "<<(calib_func_simple->GetChisquare()/calib_func_simple->GetNDF())<<std::endl;
	std::cout<<"complex fit chi2/NDF = "<<calib_func_complex->GetChisquare()<<"/"<<calib_func_complex->GetNDF()
	         <<" = "<<(calib_func_complex->GetChisquare()/calib_func_complex->GetNDF())<<std::endl;
	// LED 275_A
	//raw fit chi2/NDF = 0.0168018/69 = 0.000243504
	//simple fit chi2/NDF = 0.000357169/69 = 5.17636e-06
	//complex fit chi2/NDF = 11.0688/69 = 0.160417
	// LED 275_B
	//raw fit chi2/NDF = 0.0113123/69 = 0.000163946
	//simple fit chi2/NDF = 0.00022226/69 = 3.22115e-06
	//complex fit chi2/NDF = 4.26198/69 = 0.0617678
	
	TCanvas cc("cc","cc",1024,800);
	calib_curve_simple->SetLineColor(kRed);
	calib_curve_simple->SetLineWidth(0);
	calib_curve_simple->SetMarkerStyle(20);
	calib_curve_simple->SetMarkerColor(kRed);
	calib_curve_simple->Draw("AP");
	calib_func_simple->SetLineWidth(1);
	calib_func_simple->SetLineColor(kRed);
	calib_func_simple->Draw("same");
	
	calib_curve_complex->SetLineWidth(0);
	calib_curve_complex->SetLineColor(kBlue);
	calib_curve_complex->SetMarkerStyle(20);
	calib_curve_complex->SetMarkerColor(kBlue);
	calib_curve_complex->Draw("same P");
	calib_func_complex->SetLineWidth(1);
	calib_func_complex->SetLineColor(kBlue);
	calib_func_complex->Draw("same");
	
	calib_curve_raw->SetLineColor(kSpring-5);
	calib_curve_raw->SetLineWidth(0);
	calib_curve_raw->SetMarkerStyle(20);
	calib_curve_raw->SetMarkerColor(kSpring-5);
	calib_curve_raw->Draw("same P");
	calib_func_raw->SetLineWidth(1);
	calib_func_raw->SetLineColor(kSpring-5);
	calib_func_raw->Draw("same");
	
	cc.Modified();
	cc.Update();
	
	/*
	TGraph residuals_raw(next_graph_point);
	TGraph residuals_simple(next_graph_point);
	TGraph residuals_complex(next_graph_point);
	for(int i=0; i<next_graph_point; ++i){
		double conctmp, difftmp;
		calib_curve_raw->GetPoint(i, conctmp, difftmp);
		residuals_raw.SetPoint(i, conctmp, difftmp - calib_func_raw->Eval(conctmp));
		calib_curve_simple->GetPoint(i, conctmp, difftmp);
		residuals_simple.SetPoint(i, conctmp, difftmp - calib_func_simple->Eval(conctmp));
		calib_curve_complex->GetPoint(i, conctmp, difftmp);
		residuals_complex.SetPoint(i, conctmp, difftmp - calib_func_complex->Eval(conctmp));
	}
	TCanvas ccc("ccc","ccc",1024,800);
	residuals_simple.SetLineColor(kRed);
	residuals_complex.SetLineColor(kBlue);
	residuals_raw.SetLineColor(kSpring-5);
	residuals_raw.Draw("AL");
	residuals_simple.Draw("same L");
	residuals_complex.Draw("same L");
	ccc.Modified();
	ccc.Update();
	*/
	
	TCanvas c4("c4","c4",1024,800);
	std::vector<double> numberline(concentrations.size());
	std::iota(numberline.begin(), numberline.end(), 0);
	std::vector<double> zeros(concentrations.size()); // errors on measurement number: 0
	std::fill(zeros.begin(), zeros.end(), 0);
	TGraphErrors g4(concentrations.size(), measnums.data(), concentrations.data(), zeros.data(), concentration_errs.data());
	g4.Draw("ALP");
	c4.Modified();
	c4.Update();
	
	//TCanvas c5("c5","c5",1024,800);
	// to try to compare shapes, scale to the same max value then plot on the same canvas
	TGraph g5(concentrations.size(), measnums.data(), calib_curve_raw->GetY()); //, zeros.data(), calib_curve_raw->GetEY());
	g5.SetLineColor(kRed);
	g5.SetMarkerColor(kRed);
	g5.SetMarkerStyle(22);
	//scale to the pad coordinates
	Float_t rightmax = (*std::max_element(g5.GetY(),g5.GetY()+g5.GetN()));
	Float_t rightmin = (*std::min_element(g5.GetY(),g5.GetY()+g5.GetN()));
	Float_t scale = gPad->GetUymax()/(1.1*rightmax);  // Uymax is 1.1 * max
	std::cout<<"scaling is "<<scale<<std::endl;
	for(int i=0; i<g5.GetN(); ++i) g5.GetY()[i] *= scale;
	g5.Draw("LP same");
	//draw associated axis on the right side
	TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
		 gPad->GetUxmax(), gPad->GetUymax(),rightmin,rightmax,510,"+L");
	axis->SetLineColor(kRed);
	axis->SetLabelColor(kRed);
	axis->SetLabelSize(0.03);
	axis->Draw();
	
	gSystem->ProcessEvents();
	while(gROOT->FindObject("cc")!=nullptr){
		gSystem->ProcessEvents();
		std::this_thread::sleep_for(std::chrono::milliseconds(100));
	}
	
	std::cout<<"returning multigraph"<<std::endl;
	return mg_calibs;
	
}


//////////////////////////////////////

TMultiGraph* Plotter::MakeMultigraph(std::string name, bool relative_to_first, bool darksub, char relop){
	
	TChain* c_led = chains.at(name);
	TChain* c_dark = chains.at("dark");
	int offset = offsets.at(name);
	int n_entries = c_led->GetEntries();
	
	int increment = 1;
	// we may wish to limit the number of graphs we plot while covering the full time span
	if(n_entries>n_multigraph_lines){
		// in this case plot every n-th entry
		increment = int(double(n_entries)/double(n_multigraph_lines));
	}
	
	std::vector<double> ref;
	if(relative_to_first){
		c_led->GetEntry(refentry);
		int dark_entry = offset+(refentry*darksperloop);
		c_dark->GetEntry(dark_entry);
		ref = values;
		
		if(darksub){
			for(int i=0; i<n_datapoints; ++i){
				ref.at(i) -= values_dark.at(i);
			}
		}
	}
	
	std::string title = name+"_g";
	if(relative_to_first) title = title+"_var";
	TMultiGraph* mg = new TMultiGraph(title.c_str(), title.c_str());
	
	
	for(int i=0; i<n_entries; i+=increment){
		
		// get LED on trace
		c_led->GetEntry(i);
		
		if(darksub){
			// get nearest dark trace
			//int dark_entry = offset + i*darksperloop;
			int dark_entry = GetNextDarkEntry(name, i);
			c_dark->GetEntry(dark_entry);
			
			// subtract dark
			for(int j=0; j<n_datapoints; ++j){
				values.at(j) = values.at(j) - values_dark.at(j);
			}
		}
		
		// if plotting relative to first trace subtract reference values
		if(relative_to_first){
			for(int j=0; j<n_datapoints; ++j){
				if(relop=='/') values.at(j) /= ref.at(j);
				if(relop=='-') values.at(j) -= ref.at(j);
			}
		}
		
		// build dark subtracted TGraph and add to the MultiGraph
		TGraph* g = new TGraph(n_datapoints, wavelengths.data(), values.data());
		std::string thistitle = title+"_"+std::to_string(i);
		g->SetTitle(thistitle.c_str());
		g->SetName(thistitle.c_str());
		// coerce to the nearest colour index in the palette
		int colour_i = int(double(i)*(double(n_colours)/double(n_entries)));
		g->SetLineColor(TColor::GetColorPalette(colour_i));
		mg->Add(g);
		
	}
	
	return mg;
}


TGraph* Plotter::MakeTrendGraph(std::string name, bool relative_to_first, bool darksub, bool darkval, int trendsample){
	
	if(darkval) darksub=false;
	
	TChain* c_led = chains.at(name);
	TChain* c_dark = chains.at("dark");
	int offset = offsets.at(name);
	
	double refval;
	if(relative_to_first){
		c_led->GetEntry(refentry);
		int dark_entry = GetNextDarkEntry(name, refentry);
		c_dark->GetEntry(dark_entry);
		
		refval = (darkval) ? values_dark.at(trendsample) : values.at(trendsample);
		
		if(darksub){
			refval -= values_dark.at(trendsample);
		}
	}
	
	TGraph* trend = new TGraph();
	std::string title = name+"_trend";
	trend->SetTitle(title.c_str());
	trend->SetName(title.c_str());
	
	int n_entries = c_led->GetEntries();
	for(int i=0; i<n_entries; ++i){
		
		// get LED on trace
		c_led->GetEntry(i);
		
		if(darkval || darksub){
			// get nearest dark trace
			int dark_entry = i;
			if(darksub) dark_entry = GetNextDarkEntry(name, i);
			c_dark->GetEntry(dark_entry);
		}
		
		double value = (darkval) ? values_dark.at(trendsample) : values.at(trendsample);
		
		if(darksub){
			// subtract dark
			value -= values_dark.at(trendsample);
		}
		
		// if plotting relative to first trace subtract reference values
		if(relative_to_first){
			value /= refval;
		}
		
		trend->SetPoint(i,i,value);
		
	}
	
	return trend;
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

void Plotter::SetTimeAxis(TH1* hist, long t0_seconds){
	// start_time is in seconds: plotted times should be in seconds relative to this
	gStyle->SetTimeOffset(t0_seconds);
	// then you set the x-axis to time format
	hist->GetXaxis()->SetTimeDisplay(1);
	hist->GetXaxis()->SetTimeFormat("%m/%d/%y %H:%M");
	hist->GetXaxis()->SetLabelSize(0.03);
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

bool Plotter::MakePure(std::string name, bool overwrite){
	
	if(!overwrite){
		// check if the pure reference file already exists
		std::string type="";
		bool exists = CheckPath(purefile, type);
		if(exists && type=="f"){
			std::cout<<"Pure file "<<purefile<<" already exists, using it"<<std::endl;
			return true;
		}
		else if(exists && type!="f"){
			std::cerr<<"Pure file "<<purefile<<" exists but is not a standard file. Please investigate"<<std::endl;
			return false;
		}
	}
	
	std::cout<<"Pure file "<<purefile<<" does not exist or overwriting it"<<std::endl;
	
	TChain* c_led = chains.at(name);
	TChain* c_dark = chains.at("dark");
	
	// get LED on trace
	int nbytes = c_led->GetEntry(0);
	if(nbytes<=0){
		std::cerr<<"COULDN'T LOAD ENTRY 0 FROM LED CHAIN"<<std::endl;
		return false;
	}
	
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
	
	// get nearest dark trace
	int dark_entry = GetNextDarkEntry(name, 0, false);
	c_dark->GetEntry(dark_entry);
	
	// subtract dark
	for(int j=0; j<n_datapoints; ++j){
		values.at(j) -= values_dark.at(j);
	}
	
	// make a TGraph from the data
	TGraph graph(n_datapoints, wavelengths.data(), values.data());
	graph.SetName("Graph");
	graph.SetTitle("Graph");
	graph.SaveAs(purefile.c_str());
	
	return true;
}

