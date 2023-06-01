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
#include <sstream>
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
	int n_datapoints;
	int num_meas;
	int n_colours;
	int refentry=0;
	
	std::map<std::string, std::pair<double,double>> calibration_data;
	
	// initialisation
	int LoadData(std::string filepattern);
	std::string GetStdoutFromCommand(std::string cmd, int bufsize=500);
	int SetBranchAddresses();
	bool LoadConcentrations(std::string filename);
	bool MakePure(std::string name, bool overwrite);
	bool CheckPath(std::string path, std::string& type);
	
	// main processing functions
	int Execute(std::string concs_file, std::string datafiles, std::string datasetname);
	TMultiGraph* FitStabilityData(std::string name, std::string dataset, std::string calib_info_file="");
	TMultiGraph* FitCalibrationData(std::string name, std::string purefile, std::string calib_info_file);
	
	// Fitting functions
	std::pair<double,double> CalculateError(TF1* abs_func, double peak1_pos, double peak2_pos);
	int FitTwoGaussians(TGraph* abs_graph, std::pair<double,double>& simple_peaks, std::pair<double,double>& simple_errs, std::pair<double,double>& simple_posns, bool plot=false);
	TF1* PureScaledPlusExtras();
	TF1* GetLaungaus(); // tried it, not great.
	TF1* GetFourGausFit();
	
	// misc
	bool GetMetrics(std::string ledname, int entrynum, TF1* purefit, std::map<std::string,std::pair<double,double>>& metric_and_err);
	int GetNextDarkEntry(std::string name, int ledon_entry_num);
	
	int sample_273=0;
	int sample_276=0;
	
	TGraph* trend_A=nullptr;
	TGraph* trend_B=nullptr;
	
};
std::string purefile;

// ================
// Main Application
// ================

int main(){
	
	TApplication myapp("rootTApp",0,0);
	
	// the known concentrations for each measurement
	std::string concs_file = "april_2023_calib/apr_2023_calibration_info.txt";
	// pattern of data files to analyse
	//std::string datafiles = "april_2023_calib/*.root";
	std::string datafiles = "../GDConcMeasure/data/2023/06/00185_25_Apr_2023_EGADS_Run_031*";
	// 'datasetname' is unique part of the generated pure reference file name
	// the generated pure reference file is named according to:
	// purefile = "../GDConcMeasure/pureDarkSubtracted_"+ ledname + "_" + dataset + ".root";
	// since it depends on the ledname, we can't pass the generated filename in directly
	std::string datasetname="apr2023cal";
	
	Plotter myplotter;
	myplotter.Execute(concs_file, datafiles, datasetname);
	
	return 0;
}

int Plotter::Execute(std::string concs_file, std::string datafiles, std::string datasetname){
	
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
	
	std::cout<<"loading data"<<std::endl;
	LoadData(datafiles);
	std::cout<<"setting branch addresses"<<std::endl;
	SetBranchAddresses();
	
	gStyle->SetPalette(kBird);
	n_colours = TColor::GetNumberOfColors();
	TCanvas* c1 = nullptr;
	
	// ============================================================================
	
	// step 1: analyse calibration data for polynomial coffiecients
	// ------------------------------------------------------------
	/*
	std::cout<<"calling FitCalibrationData"<<std::endl;
	TMultiGraph* mg_calib = new TMultiGraph("mg_data","mg_data");
	mg_calib->Add(FitCalibrationData("275_A", datasetname, concs_file));
	mg_calib->Add(FitCalibrationData("275_B", datasetname, concs_file));
	
	// display output for user
	if(gROOT->FindObject("c1")==nullptr) c1 = new TCanvas("c1","c1",1024,800);
	c1->cd();
	mg_calib->Draw("AP");
	*/
	
	// step 2: analyse data to extract plots of concentration stability
	// ----------------------------------------------------------------
	// we can do it for each of the calibration fit methods
	TMultiGraph* mg_data = new TMultiGraph("mg_data","mg_data");
	mg_data->Add(FitStabilityData("275_A", datasetname/*, concs_file*/));  // concs_file optional
	mg_data->Add(FitStabilityData("275_B", datasetname/*, concs_file*/));  // for validation
	
	///*
	// for comparison, plot true concentrations vs extracted concentrations
	// (note: don't do this if passing concs_file to FitStabilityData)
	bool ok = LoadConcentrations(concs_file);
	std::vector<double> calib_concs;
	for(auto&& acalibmeas : calibration_data){
		calib_concs.push_back(acalibmeas.second.first);
	}
	std::sort(calib_concs.begin(),calib_concs.end());
	std::vector<double> numberline(calib_concs.size());
	std::iota(numberline.begin(),numberline.end(),0);
	TGraph* g_cal = new TGraph(calib_concs.size(), numberline.data(), calib_concs.data());
	g_cal->SetTitle("true conc");
	g_cal->SetLineColor(kMagenta);
	g_cal->SetMarkerColor(kMagenta);
	g_cal->SetMarkerStyle(30);
	mg_data->Add(g_cal);
	//*/
	
	/*
	// if doing a validation plot of measured conc vs true conc, add a line of y=x
	bool ok = LoadConcentrations(concs_file);
	std::vector<double> calib_concs;
	for(auto&& acalibmeas : calibration_data){
		calib_concs.push_back(acalibmeas.second.first);
	}
	std::sort(calib_concs.begin(),calib_concs.end());
	std::vector<double> numberline(calibration_data.size());
	std::iota(numberline.begin(),numberline.end(),0);
	for(auto& elem : numberline) elem *= calib_concs.back()/numberline.size();
	TGraph* g_diag = new TGraph(numberline.size(), numberline.data(), numberline.data());
	g_diag->SetLineColor(kBlack);
	g_diag->SetMarkerStyle(0);
	g_diag->SetLineStyle(9);
	mg_data->Add(g_diag);
	*/
	
	// display output for user
	if(gROOT->FindObject("c1")==nullptr) c1 = new TCanvas("c1","c1",1024,800);
	c1->cd();
	mg_data->Draw("ALP");
	
	// ============================================================================
	
	while(gROOT->FindObject("c1")!=nullptr){
		gSystem->ProcessEvents();
		std::this_thread::sleep_for(std::chrono::milliseconds(100));
	}
	
	
	return 0;
}

//==================

// Main Function Part 1 - Extracting the Calibration Curve
// -------------------------------------------------------

TMultiGraph* Plotter::FitCalibrationData(std::string name, std::string dataset, std::string concs_file){
	
	// loop over the data, do dark subtraction, fit the dark-subtracted data, extract the absorption graph.
	// Fit the absorption graph and pull peak heights and difference.
	// Finally, look up the known concentration and add a point mapping peak height diff to concentration.
	
	// load calibration info that specifies our concentrations
	std::cout<<"loading true concentrations"<<std::endl;
	bool ok = LoadConcentrations(concs_file);
	if(!ok) return nullptr;
	
	// build the name of the pure file based on the led name and an arbitrary dataset name.
	// MakePure will make it if required, or do nothing if it already exists.
	// because it depends on the led name, we can't pass the filename in directly.
	// note that the 'purefile' variable is a global, as it will be opened by the global function
	// 'PureFuncv2', which returns a TF1 functional fit based on the pure reference trace in 'purefile'.
	purefile="../GDConcMeasure/pureDarkSubtracted_"+ name + "_" + dataset + ".root";
	ok = MakePure(name, false);
	if(not ok){
		return nullptr;
	}
	
	std::cout<<"getting reference pure"<<std::endl;
	TF1* purefit = PureScaledPlusExtras();
	//TF1* purefit = GetLaungaus();
	
	TChain* c_led = chains.at(name);
	int n_entries = c_led->GetEntries();
	//n_entries = 53;
	
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
	
	// as a back-check we'll also plot the concentrations vs measurement
	// and peak height difference vs measurement indepdently
	// (as well as the usual one against the other)
	std::vector<double> concentrations, concentration_errs, measnums;
	
	int next_graph_point=0;
	double lastconc=-99;
	std::cout<<"looping over calibration entries"<<std::endl;
	for(int i=0; i<n_entries; i++){
		
		std::cout<<i<<": ";
		
		std::map<std::string, std::pair<double,double>> metric_and_err;
		ok = GetMetrics(name, i, purefit, metric_and_err);
		//std::cout<<"GetMetrics returned "<<ok<<std::endl;
		// XXX should we 'continue' if false?
		
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
		//std::cout<<"getting concentration"<<std::endl;
		double current_conc = calibration_data.at(current_file).first;
		double current_concerr = calibration_data.at(current_file).second;
		//if(lastconc==current_conc) continue; /// hack because i messed up, skip this file
		//lastconc=current_conc;
		
		concentrations.push_back(current_conc);
		concentration_errs.push_back(current_concerr);
		
		//std::cout<<"adding point to calibration curve"<<std::endl;
		// add to the calibration data curve
		calib_curve_raw->SetPoint(next_graph_point,current_conc,metric_and_err.at("raw").first);
		calib_curve_simple->SetPoint(next_graph_point,current_conc,metric_and_err.at("simple").first);
		calib_curve_complex->SetPoint(next_graph_point,current_conc,metric_and_err.at("complex").first);
		
		calib_curve_raw->SetPointError(next_graph_point,current_concerr,metric_and_err.at("raw").second);
		calib_curve_simple->SetPointError(next_graph_point,current_concerr, metric_and_err.at("simple").second);
		calib_curve_complex->SetPointError(next_graph_point,current_concerr, metric_and_err.at("complex").second);
		++next_graph_point;
		
		// extract measurement number from filename since they're not contiguous in a TChain
		// unless the filenames are alphabetical
		// 00113_06Jun22Calib_99.root
		std::string thismeasnum = current_file.substr(current_file.find_last_of("_")+1,current_file.length()-24);
		int this_meas_num = std::stoi(thismeasnum);
		//std::cout<<"file "<<current_file<<" is measurement "<<this_meas_num<<std::endl;
		measnums.push_back(this_meas_num);
		/*
		std::cout<<"raw diff: "<<raw_peaks.first<<" - "<<raw_peaks.second<<" = "<<raw_peaks.first-raw_peaks.second<<std::endl;
		std::cout<<"simple diff: "<<simple_peaks.first<<" - "<<simple_peaks.second<<" = "<<simple_peaks.first-simple_peaks.second<<std::endl;
		std::cout<<"complex diff: "<<peak1_height<<" - "<<peak2_height<<" = "<<peak1_height-peak2_height<<std::endl;
		std::cout<<"raw err: "<<0<<", simple err: "<<simplerr<<", complexerr: "<<complexerr<<std::endl;
		*/
		
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
	
	// print for info
	/*
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
	*/
	
	// save to file
	std::string parfilename = "../GDConcMeasure/fitParams_"+ name + "_" + dataset + ".root";
	TFile* fparams = new TFile(parfilename.c_str(),"RECREATE");
	TTree* tparams = new TTree("params","params");
	//std::vector<double>* calib_coefficients_rawp;
	//std::vector<double>* calib_coefficients_simplep;
	//std::vector<double>* calib_coefficients_complexp;
	tparams->Branch("raw", &calib_coefficients_raw);
	tparams->Branch("simple", &calib_coefficients_simple);
	tparams->Branch("complex", &calib_coefficients_complex);
	tparams->Fill();
	fparams->Write();
	fparams->Close();
	delete fparams;
	
	TCanvas* cc = nullptr;
	//cc = new TCanvas("cc","cc",1024,800);
	// drawing the tgraphs and functions here seems to kill them when we close the canvas
	// so that when we then try to draw both LEDs together in Execute, it segs. :(
	if(name=="275_A"){
		calib_curve_raw->SetMarkerStyle(20);
		calib_curve_simple->SetMarkerStyle(20);
		calib_curve_complex->SetMarkerStyle(20);
	} else if(name=="275_B"){
		calib_curve_raw->SetMarkerStyle(34);
		calib_curve_simple->SetMarkerStyle(34);
		calib_curve_complex->SetMarkerStyle(34);
	}
	calib_curve_simple->SetLineColor(kRed);
	calib_curve_simple->SetLineWidth(0);
	calib_curve_simple->SetMarkerColor(kRed);
	if(cc) calib_curve_simple->Draw("AP");
	calib_func_simple->SetLineWidth(1);
	calib_func_simple->SetLineColor(kRed);
	if(cc) calib_func_simple->Draw("same");
	
	calib_curve_complex->SetLineWidth(0);
	calib_curve_complex->SetLineColor(kBlue);
	calib_curve_complex->SetMarkerColor(kBlue);
	if(cc) calib_curve_complex->Draw("same P");
	calib_func_complex->SetLineWidth(1);
	calib_func_complex->SetLineColor(kBlue);
	if(cc) calib_func_complex->Draw("same");
	
	calib_curve_raw->SetLineColor(kSpring-5);
	calib_curve_raw->SetLineWidth(0);
	calib_curve_raw->SetMarkerColor(kSpring-5);
	if(cc) calib_curve_raw->Draw("same P");
	calib_func_raw->SetLineWidth(1);
	calib_func_raw->SetLineColor(kSpring-5);
	if(cc) calib_func_raw->Draw("same");
	
	if(cc){
		cc->Modified();
		cc->Update();
	}
	
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
	
	/*
	// draw concentration vs measurement number
	// and metric vs measurement number
	TCanvas c4("c4","c4",1024,800);
	std::vector<double> numberline(concentrations.size());
	std::iota(numberline.begin(), numberline.end(), 0);
	std::vector<double> zeros(concentrations.size()); // errors on measurement number: 0
	std::fill(zeros.begin(), zeros.end(), 0);
	TGraphErrors g4(concentrations.size(), measnums.data(), concentrations.data(), zeros.data(), concentration_errs.data());
	g4.SetTitle("true_conc;measurement num;concentration [%]");
	g4.Draw("ALP");
	c4.Modified();
	c4.Update();
	
	//TCanvas c5("c5","c5",1024,800);
	// to try to compare shapes, scale to the same max value then plot on the same canvas
	TGraph g5(concentrations.size(), measnums.data(), calib_curve_raw->GetY()); //, zeros.data(), calib_curve_raw->GetEY());
	g5.SetTitle("metric_raw;measurement num;raw metric");
	g5.SetLineColor(kRed);
	g5.SetMarkerColor(kRed);
	g5.SetMarkerStyle(22);
	//scale to the pad coordinates
	Float_t rightmax = (*std::max_element(g5.GetY(),g5.GetY()+g5.GetN()));
	Float_t rightmin = (*std::min_element(g5.GetY(),g5.GetY()+g5.GetN()));
	Float_t scale = gPad->GetUymax()/(1.1*rightmax);  // Uymax is 1.1 * max
	//std::cout<<"scaling is "<<scale<<std::endl;
	for(int i=0; i<g5.GetN(); ++i) g5.GetY()[i] *= scale;
	g5.Draw("LP same");
	//draw associated axis on the right side
	TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
		 gPad->GetUxmax(), gPad->GetUymax(),rightmin,rightmax,510,"+L");
	axis->SetLineColor(kRed);
	axis->SetLabelColor(kRed);
	axis->SetLabelSize(0.03);
	axis->Draw();
	*/
	
	// hold to allow user to inspect, including any other local graphs
	std::cout<<"waiting for user to close canvas cc"<<std::endl;
	gSystem->ProcessEvents();
	while(gROOT->FindObject("cc")!=nullptr){
		gSystem->ProcessEvents();
		std::this_thread::sleep_for(std::chrono::milliseconds(100));
	}
	std::cout<<"Done!"<<std::endl;
	
	TMultiGraph* mg_calibs = new TMultiGraph(title.c_str(), title.c_str());
	mg_calibs->Add(calib_curve_simple);
	mg_calibs->Add(calib_curve_complex);
	mg_calibs->Add(calib_curve_raw);
	return mg_calibs;
	
}

// Main Function Part 2 - Fitting Stability Data
// ---------------------------------------------

TMultiGraph* Plotter::FitStabilityData(std::string name, std::string dataset, std::string concs_file){
	
	if(concs_file!=""){
		// validating calibration data: plot measured conc vs true conc
		// load calibration info that specifies our concentrations
		std::cout<<"loading true concentrations"<<std::endl;
		bool ok = LoadConcentrations(concs_file);
		if(!ok) return nullptr;
	}
	
	// get the pure function.
	purefile="../GDConcMeasure/pureDarkSubtracted_"+ name + "_" + dataset + ".root";
	std::cout<<"getting pure "<<purefile<<std::endl;
	TF1* purefit = PureScaledPlusExtras();
//	TF1* purefit = GetLaungaus();
	
	TChain* c_led = chains.at(name);
	TChain* c_dark = chains.at("dark");
	int n_entries = c_led->GetEntries();
	
	std::string title = name+"_g";
	TMultiGraph* mg = new TMultiGraph(title.c_str(), title.c_str());
	
	// FIXME TODO make them TGraphErrors
	TGraph* g_concentrations_raw = new TGraph(n_entries);
	g_concentrations_raw->SetMarkerColor(kSpring-5);
	g_concentrations_raw->SetLineColor(kSpring-5);
	TGraph* g_concentrations_simple = new TGraph(n_entries);
	g_concentrations_simple->SetMarkerColor(kRed);
	g_concentrations_simple->SetLineColor(kRed);
	TGraph* g_concentrations_complex = new TGraph(n_entries);
	g_concentrations_complex->SetMarkerColor(kBlue);
	g_concentrations_complex->SetLineColor(kBlue);
	if(name=="275_A"){
		g_concentrations_raw->SetMarkerStyle(20);
		g_concentrations_simple->SetMarkerStyle(20);
		g_concentrations_complex->SetMarkerStyle(20);
	} else if(name=="275_B"){
		g_concentrations_raw->SetMarkerStyle(34);
		g_concentrations_simple->SetMarkerStyle(34);
		g_concentrations_complex->SetMarkerStyle(34);
	}
	mg->Add(g_concentrations_raw);
	mg->Add(g_concentrations_simple);
	mg->Add(g_concentrations_complex);
	
	// for converting metric to concentration
	std::string parfilename = "../GDConcMeasure/fitParams_"+ name + "_" + dataset + ".root";
	std::cout<<"getting calibration params from "<<parfilename<<std::endl;
	TFile* fparams = TFile::Open(parfilename.c_str(),"READ");
	TTree* tparams = (TTree*)fparams->Get("params");
	std::vector<double> calib_coefficients_raw;
	std::vector<double> calib_coefficients_simple;
	std::vector<double> calib_coefficients_complex;
	std::vector<double>* calib_coefficients_rawp = &calib_coefficients_raw;
	std::vector<double>* calib_coefficients_simplep = &calib_coefficients_simple;
	std::vector<double>* calib_coefficients_complexp = &calib_coefficients_complex;
	tparams->SetBranchAddress("raw", &calib_coefficients_rawp);
	tparams->SetBranchAddress("simple", &calib_coefficients_simplep);
	tparams->SetBranchAddress("complex", &calib_coefficients_complexp);
	tparams->GetEntry(0);
	tparams->ResetBranchAddresses();
	fparams->Close();
	TF1 calib_curve_raw("calib_raw", "pol6", 0, 0.4);
	TF1 calib_curve_simple("calib_simple", "pol6", 0, 0.4);
	TF1 calib_curve_complex("calib_complex", "pol6", 0, 0.4);
	calib_curve_raw.SetParameters(calib_coefficients_raw.data());
	calib_curve_simple.SetParameters(calib_coefficients_simple.data());
	calib_curve_complex.SetParameters(calib_coefficients_complex.data());
	
	int n_concentration_vals=-1;
	std::cout<<"analysing "<<n_entries<<" entries"<<std::endl;
	for(int i=0; i<n_entries; i++){
		
		std::cout<<i<<": ";
		std::map<std::string, std::pair<double,double>> metric_and_err;
		bool ok = GetMetrics(name, i, purefit, metric_and_err);
		//if(!ok) continue;
		
		std::string current_file = c_led->GetTree()->GetCurrentFile()->GetName();
		std::cout<<"filename: "<<current_file<<std::endl;
		
		std::cout<<"metrics:\nraw: "<<metric_and_err.at("raw").first
		         <<", simple: "<<metric_and_err.at("simple").first
		         <<", complex: "<<metric_and_err.at("complex").first<<std::endl;
		
		// solve for concentration (x) from absorbance (y), with 0.01 < x < 0.21
		// XXX note that the pol6 is not monotonic, so we MUST constrain the range
		// to find a solution in the first monotonic region.
		// (TODO perhaps we could find the first point of inflexion and set that as the upper limit?)
		// (this is only really a problem when re-analysing calibration data as we go right
		// to the edge of the calibrated range - data won't do that)
		double conc_raw = calib_curve_raw.GetX(metric_and_err.at("raw").first, 0.001, 0.23);
		double conc_simple = calib_curve_simple.GetX(metric_and_err.at("simple").first, 0.001, 0.23);
		double conc_complex = calib_curve_complex.GetX(metric_and_err.at("complex").first, 0.001, 0.23);
	
		std::cout<<"concentrations:\nraw: "<<conc_raw
		         <<", simple: "<<conc_simple
		         <<", complex: "<<conc_complex<<std::endl;
		
		++n_concentration_vals;
		double x_val = n_concentration_vals;
		
		// if doing validation of calibration data, plot measured conc vs true conc
		// this should be a straight line along y=x, with small fluctuations 
		// from fitting the calibration datapoints with a pol6
		if(concs_file!=""){
			// strip off preceding path to get key in concentrations map
			if(current_file.find("/")!=std::string::npos){
				current_file = current_file.substr(current_file.find_last_of("/")+1,std::string::npos);
			}
			if(calibration_data.count(current_file)==0){
				std::cerr<<"Couldn't find file "<<current_file<<" in calibration data map!"<<std::endl;
				for(auto&& am : calibration_data){
					std::cerr<<am.first<<", ";
				}
				std::cout<<std::endl;
				break;
				return nullptr;
			}
			x_val = calibration_data.at(current_file).first;
			std::cout<<"true conc: "<<x_val<<std::endl;
		}
		
		g_concentrations_raw->SetPoint(n_concentration_vals,x_val,conc_raw);
		g_concentrations_simple->SetPoint(n_concentration_vals,x_val,conc_simple);
		g_concentrations_complex->SetPoint(n_concentration_vals,x_val,conc_complex);
		
		/* FIXME TODO
		// error on concentration is error on height times gradient at that point
		double conc_err_raw = metric_and_err.at("raw").second * calib_curve_raw.Derivative(metric_and_err.at("raw").first);
		double conc_err_simple = metric_and_err.at("simple").second * calib_curve_simple.Derivative(metric_and_err.at("simple").first);
		double conc_err_complex = metric_and_err.at("complex").second * calib_curve_complex.Derivative(metric_and_err.at("complex").first);
		g_concentrations_raw->SetPointError(i,wavelength_errors.at(i),conc_err_raw);
		g_concentrations_simple->SetPointError(i,wavelength_errors.at(i),conc_err_simple);
		g_concentrations_complex->SetPointError(i,wavelength_errors.at(i),conc_err_complex);
		*/
		
	}
	
	g_concentrations_raw->Set(n_concentration_vals);
	g_concentrations_simple->Set(n_concentration_vals);
	g_concentrations_complex->Set(n_concentration_vals);
	
	std::cout<<"returning multigraph"<<std::endl;
	return mg;
	
}

// ============
// Fitting Data
// ============

bool Plotter::GetMetrics(std::string ledname, int entrynum, TF1* purefit, std::map<std::string,std::pair<double,double>>& metric_and_err){
	
	TChain* c_led = chains.at(ledname);
	TChain* c_dark = chains.at("dark");
	
	// get LED on trace
	int nbytes = c_led->GetEntry(entrynum);
	if(nbytes<=0){
		std::cerr<<"COULDN'T LOAD ENTRY "<<entrynum<<" FROM LED CHAIN"<<std::endl;
		exit(-1);
	}
	
	// get nearest dark trace
	int dark_entry = GetNextDarkEntry(ledname, entrynum);
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
	
	std::string title = ledname+"_g";
	std::string thistitle = title+"_"+std::to_string(entrynum);
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
	std::cout<<"fitting sideband"<<std::endl;
	g_sideband->Fit(purefit,"RNMQ");
	///*
	TCanvas ccnew("ccnew","ccnew",1024,800);
	g_sideband->SetMarkerColor(kRed);
	g_sideband->SetMarkerStyle(2);
	g_inband->SetMarkerColor(kBlue);
	g_inband->SetMarkerStyle(2);
	purefit->Draw();
	g_inband->Draw("P same");
	g_sideband->Draw("P same");
	// these have to go *after* a draw or the graph has no axes
	// do not seem to work, however...
	//g_inband->GetYaxis()->SetRangeUser(-0.1,0.6);
	//g_inband->GetXaxis()->SetRangeUser(240,320);
	ccnew.Modified();
	ccnew.Update();
	gSystem->ProcessEvents();
	std::cout<<"waiting for user to close canvas"<<std::endl;
	while(gROOT->FindObject("ccnew")!=nullptr){
		gSystem->ProcessEvents();
		std::this_thread::sleep_for(std::chrono::milliseconds(100));
	}
	std::this_thread::sleep_for(std::chrono::milliseconds(100));
	//*/
	
	std::cout<<"generating absorption graph"<<std::endl;
	// calculate absorbance from ratio of fit to data in absorption region
	TGraphErrors* g_abs = new TGraphErrors(inband_values.size());
	std::string restitle = "g_abs_"+std::to_string(entrynum);
	for(int k=0; k<inband_values.size(); ++k){
		double wlval, dataval;
		g_inband->GetPoint(k, wlval, dataval);
		double fitval = purefit->Eval(wlval);
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
	std::cout<<"getting raw peaks"<<std::endl;
	std::pair<double,double> raw_peaks;
	double wlval;
	g_abs->GetPoint(sample_273, wlval, raw_peaks.first);
	g_abs->GetPoint(sample_276, wlval, raw_peaks.second);
	
	// if peaks are negative, we could coerce to 0...
	// but this tends to only happen for the first value, and results in a point at (0,0)
	// this point then seems to be way out of the trend, so probably throws the fit off.
	// instead, just skip this concentration
	if(raw_peaks.first<0 || raw_peaks.second<0){
		std::cerr<<"skipping concentration measurement "<<entrynum
		         <<" as one of the peaks is negative"<<std::endl;
		raw_peaks.first=std::max(0.,raw_peaks.first);
		raw_peaks.second=std::max(0.,raw_peaks.second);
		delete g_inband;
		delete g_sideband;
		delete g_other;
		delete g_abs;
		return false;
	}
	
	// just taking the data value at a specific wavelength may be a bit naff
	// as the sampling is sparse, and we may slightly miss the absorption peak.
	// We may do better by interpolating the peak to find a better estimate of maximum.
	// We can do this by fitting the peaks with gaussians, but since these are peaks
	// on a non-uniform background, we can only fit within a narrow region close to the peak.
	std::cout<<"getting simple peaks"<<std::endl;
	std::pair<double,double> simple_peaks, simple_errs, simple_posns;
	int ok = FitTwoGaussians(g_abs, simple_peaks, simple_errs, simple_posns);
	// total error from adding in quadrature
	double simplerr = sqrt(TMath::Sq(simple_errs.first)+TMath::Sq(simple_errs.second));
	
	// fit it with a combination of 4 gaussians
	std::cout<<"getting complex peaks"<<std::endl;
	TF1* abs_func = GetFourGausFit();  // we own the returned TF1
	// The fit has a tendency to screw up, so seed the initial values based on the raw fit.
	// These initial values will be over-estimates, since the complex fit peak heights also
	// have contributions from the shoulder gaussians, so peak1 component amplitude is less,
	abs_func->SetParameter("peak 1 amp",raw_peaks.first);
	abs_func->SetParameter("peak 2 amp",raw_peaks.second/raw_peaks.first);
	// also set the other components relative to this - no longer needed, par definitions are already relative
	//abs_func->SetParameter("RH shoulder amp",raw_peaks.first*0.5);
	//abs_func->SetParameter("LH shoulder amp",raw_peaks.second*0.2);
	TFitResultPtr frptr = g_abs->Fit(abs_func,"RQNS");
	
	if( false /*frptr->IsEmpty() || !frptr->IsValid() || frptr->Status()!=0*/){
		// fit failed; skip value? nah, this happens a lot.
		delete g_inband;
		delete g_sideband;
		delete g_other;
		delete g_abs;
		delete abs_func;
		return false;
	}
	
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
	//std::cout<<"peak1_pos="<<peak1_pos<<", 2 pos="<<peak2_pos<<std::endl;
	// extract the peak heights by evaluating the function at that x
	double peak1_height  = abs_func->Eval(peak1_pos);
	double peak2_height  = abs_func->Eval(peak2_pos);
	if(TMath::IsNaN(peak1_height) || TMath::IsNaN(peak2_height)){
		//continue;
		peak1_height=0;
		peak2_height=0;
	}
	
	// calculate error on peak heights
	std::pair<double,double> complexerrp = CalculateError(abs_func, peak1_pos, peak2_pos);
	double complexerr = sqrt(TMath::Sq(complexerrp.first)+TMath::Sq(complexerrp.second));
	
	// FIXME better errors, error on peak 2 needs to account for error on peak 1?
	// error on raw result needs to come from error from spectrometer
	
	metric_and_err["raw"]=std::pair<double,double>{raw_peaks.first-raw_peaks.second,0.1};
	metric_and_err["simple"]=std::pair<double,double>{simple_peaks.first-simple_peaks.second,simplerr};
	metric_and_err["complex"]=std::pair<double,double>{peak1_height-peak2_height,complexerr};
	
	// cleanup
	delete g_inband;
	delete g_sideband;
	delete g_other;
	delete g_abs;
	delete abs_func;
	
	return true;
}

// ========================
// Initialisation Functions
// ========================


int Plotter::LoadData(std::string filepattern){
	
	std::string lscommand = "ls -v " + filepattern + " 2>/dev/null";
	//std::cout<<"lscommand is '"<<lscommand<<"'"<<std::endl;
	std::string fileliststring = GetStdoutFromCommand(lscommand);
	fileliststring.erase(fileliststring.find_last_not_of(" \t\n\015\014\013")+1);  // strip trailing whitespace
	
	std::stringstream ssl;
	ssl << fileliststring;
	std::string nextfilestring;
	std::vector<std::string> paths;
	
	while(getline(ssl,nextfilestring)){
		std::cout<<"adding "<<nextfilestring<<std::endl;
		for(auto&& c : chains){
			c.second->Add(nextfilestring.c_str());
		}
		c_dark->Add(nextfilestring.c_str());
	}
	
	// get number of measurements
	num_meas = c_275_A->GetEntries();
	std::cout<<"had "<<num_meas<<" data files"<<std::endl;
	
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

// =================
// Fitting Functions
// =================

// Fitting the Pure
// ----------------

// this is the c-function underlying our TF1 pure functional fit
// my local installation (v6.06) does not like using lambdas - use a global function instead
TGraphErrors* dark_subtracted_pure =nullptr;
double PureFuncv2(double* x, double* par){
	if(dark_subtracted_pure==nullptr){
		TFile *_file0 = TFile::Open(purefile.c_str());
		dark_subtracted_pure = (TGraphErrors*)_file0->Get("Graph");
		// minor optimization; if the data in the TGraph is sorted by x value
		//(which it should be for us), then setting the following option can speed up Eval() calls
		//dark_subtracted_pure->SetBit(TGraph::kIsSortedX); // XXX note only available on newer ROOT
	}
	// par [0] = y-scaling
	// par [1] = x-scaling
	// par [2] = x-offset
	// par [3] = y-offset
	// par [4] = linear baseline offset
	// par [5] = shoulder gaussian scaling
	// par [6] = shoulder gaussian centre, restricted to > 282nm (RH shoulder)
	// par [7] = shoulder gaussian spread
	double purepart = par[0] * dark_subtracted_pure->Eval((par[1]*(x[0]-276))+276+par[2]);
	double linpart = (par[4] * (x[0]-276)) + par[3];
	double shoulderpart = par[5]*exp(-0.5*TMath::Sq((x[0]-282.-abs(par[6]))/par[7]));
	double retval = purepart + linpart + shoulderpart;
	
	return retval;
}

// this concstruct the TF1 functional fit, based on the above c-function
TF1* Plotter::PureScaledPlusExtras(){
	
	if(dark_subtracted_pure) delete dark_subtracted_pure;
	dark_subtracted_pure=nullptr;
	
	// construct functional fit. We'll scale and add a linear background.
	// limit the pure function to a region in which we have light - no point fitting outside this region
	static int purever=0;
	purever++;
	std::string name="purev2_fct"+std::to_string(purever);
	const int wave_min = 260, wave_max = 300, numb_of_fitting_parameters = 8;
	TF1* pure = new TF1(name.c_str(), PureFuncv2, wave_min, wave_max, numb_of_fitting_parameters);
	// set default parameters
	//pure->SetParameters(1.,1.,0.,0.,0.,0.,0.,10.);
	pure->SetParameters(17.,0.8,1.8,-0.9,-58,-10453,4.,9.);
	
	// set parameter limits
	pure->SetParLimits(0,0,30);     // stretch y
	pure->SetParLimits(1,0.8,1.2);  // stretch x
	pure->SetParLimits(2,-20,20);   // x offset
	pure->SetParLimits(3,-10,10);   // y offset
	pure->SetParLimits(4,-500,10);   // linear baseline gradient
	pure->SetParLimits(5,-15000,500); // shoulder gaussian amplitude XXX reduce
	pure->SetParLimits(6,0,60);     // shoulder gaussian position - 282
	pure->SetParLimits(7,0,50);     // shoulder gaussian width
	
	pure->SetParName(0,"y scaling");
	pure->SetParName(1,"x scaling");
	pure->SetParName(2,"x offset");
	pure->SetParName(3,"y offset");
	pure->SetParName(4,"linear gradient");
	pure->SetParName(5,"shoulder amp");
	pure->SetParName(6,"shoulder pos");
	pure->SetParName(7,"shoulder width");
	
	return pure;
	
}

// Fitting the Absorbance
// ----------------------

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
	//std::cout<<"simple peak fit gaus 1 has ampltiude "<<gausamp1<<"+-"<<gaus1amperr<<std::endl;
	//std::cout<<"simple peak fit gaus 2 has ampltiude "<<gausamp2<<"+-"<<gaus2amperr<<std::endl;
	
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

// Fitting the Absorbance v2
// -------------------------

TF1* Plotter::GetFourGausFit(){
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
		
		//std::cout<<"complex peak fit gaus "<<i<<" has ampltiude "
		//         <<agaus.GetParameter(0)<<"+-"<<error_centre<<" scaled by "
		//         <<relative_amp<<" to "<<error_here<<" at peak 1";
		
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

// ================
// Helper Functions
// ================

int Plotter::GetNextDarkEntry(std::string name, int ledon_entry_num){
	
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
	
	return darkentry;
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
	std::cout<<"had "<<calibration_data.size()<<" calibration concentrations"<<std::endl;
	
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
	int dark_entry = GetNextDarkEntry(name, 0);
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

// ==================

// Langaus
// -------
double combinedfunc(double* x, double* par);
double landau_pdf(double x, double xi=1, double x0=0);
double Landau(double x, double mu, double sigma=1, bool norm=false);
double Gaus(double x, double mean=0, double sigma=1, bool norm=false);
double langaufun(double *x, double *par);

TF1* Plotter::GetLaungaus(){
	double wave_min = 260, wave_max = 300;
	TF1* fit = new TF1("fit",langaufun,wave_min,wave_max,4);
	fit->SetParName(0,"landau width");
	fit->SetParName(1,"landau centre");
	fit->SetParName(2,"integral of landaugaus");
	fit->SetParName(3,"width of centre gaus component");
	
	// these two are fairly straightforward
	fit->SetParLimits(1,260,280);    // centre of peak
	fit->SetParLimits(2,0.5E6,5E6);  // amplitude scaling
	
	// these two are the tricky ones that control shape: they both change the width
	// and asymmetry, but they seem more or less to do the same thing...
	// and still not achieve what's wanted. :(
	fit->SetParLimits(0,0.5,7.);     // below 0.5 this breaks down
	fit->SetParLimits(3,1,7);        // 
	
	// set initial values
	std::vector<double> initvals{1., 273., 1.0E6, 5.};
	fit->SetParameters(initvals.data());
	return fit;
}

// copied from $ROOTSYS/tutorials/fit/langaus.C
// dependant functions from TMath are copied below.
double langaufun(double *x, double *par)
{
	
	//Fit parameters:
	//par[0]=Width (scale) parameter of Landau density
	//par[1]=Most Probable (MP, location) parameter of Landau density
	//par[2]=Total area (integral -inf to inf, normalization constant)
	//par[3]=Width (sigma) of convoluted Gaussian function
	//
	//In the Landau distribution (represented by the CERNLIB approximation),
	//the maximum is located at x=-0.22278298 with the location parameter=0.
	//This shift is corrected within this function, so that the actual
	//maximum is identical to the MP parameter.
	
	// Numeric constants
	double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
	double mpshift  = -0.22278298;       // Landau maximum location
	
	// Control constants
	double np = 100.0;      // number of convolution steps
	double sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
	
	// Variables
	double xx;
	double mpc;
	double fland;
	double sum = 0.0;
	double xlow,xupp;
	double step;
	double i;
	
	
	// MP shift correction
	mpc = par[1] - mpshift * par[0];
	
	// Range of convolution integral
	xlow = x[0] - sc * par[3];
	xupp = x[0] + sc * par[3];
	
	step = (xupp-xlow) / np;
	
	// Convolution integral of Landau and Gaussian by sum
	for(i=1.0; i<=np/2; i++) {
		xx = xlow + (i-.5) * step;
		fland = Landau(xx,mpc,par[0]) / par[0];
		sum += fland * Gaus(x[0],xx,par[3]);
		
		xx = xupp - (i-.5) * step;
		fland = Landau(xx,mpc,par[0]) / par[0];
		sum += fland * Gaus(x[0],xx,par[3]);
	}
	
	return (par[2] * step * sum * invsq2pi / par[3]);
}

// copied from https://root.cern.ch/doc/master/TMath_8cxx_source.html#l00448
double Gaus(double x, double mean, double sigma, bool norm)
{
	if (sigma == 0) return 1.e30;
	double arg = (x-mean)/sigma;
	// for |arg| > 39  result is zero in double precision
	if (arg < -39.0 || arg > 39.0) return 0.0;
	double res = std::exp(-0.5*arg*arg);
	if (!norm) return res;
	return res/(2.50662827463100024*sigma); //sqrt(2*Pi)=2.50662827463100024
}

// copied from https://root.cern.ch/doc/master/TMath_8cxx_source.html#l00448
double Landau(double x, double mu, double sigma, bool norm)
{
	if (sigma <= 0) return 0;
	double den = landau_pdf( (x-mu)/sigma );
	if (!norm) return den;
	return den/sigma;
}

// copied from https://root.cern.ch/doc/master/PdfFuncMathCore_8cxx_source.html#l00191 
double landau_pdf(double x, double xi, double x0)
{
	// LANDAU pdf : algorithm from CERNLIB G110 denlan
	// same algorithm is used in GSL
	
	static double p1[5] = {0.4259894875,-0.1249762550, 0.03984243700, -0.006298287635,   0.001511162253};
	static double q1[5] = {1.0         ,-0.3388260629, 0.09594393323, -0.01608042283,    0.003778942063};
	
	static double p2[5] = {0.1788541609, 0.1173957403, 0.01488850518, -0.001394989411,   0.0001283617211};
	static double q2[5] = {1.0         , 0.7428795082, 0.3153932961,   0.06694219548,    0.008790609714};
	
	static double p3[5] = {0.1788544503, 0.09359161662,0.006325387654, 0.00006611667319,-0.000002031049101};
	static double q3[5] = {1.0         , 0.6097809921, 0.2560616665,   0.04746722384,    0.006957301675};
	
	static double p4[5] = {0.9874054407, 118.6723273,  849.2794360,   -743.7792444,      427.0262186};
	static double q4[5] = {1.0         , 106.8615961,  337.6496214,    2016.712389,      1597.063511};
	
	static double p5[5] = {1.003675074,  167.5702434,  4789.711289,    21217.86767,     -22324.94910};
	static double q5[5] = {1.0         , 156.9424537,  3745.310488,    9834.698876,      66924.28357};
	
	static double p6[5] = {1.000827619,  664.9143136,  62972.92665,    475554.6998,     -5743609.109};
	static double q6[5] = {1.0         , 651.4101098,  56974.73333,    165917.4725,     -2815759.939};
	
	static double a1[3] = {0.04166666667,-0.01996527778, 0.02709538966};
	
	static double a2[2] = {-1.845568670,-4.284640743};
	
	if (xi <= 0) return 0;
	double v = (x - x0)/xi;
	double u, ue, us, denlan;
	if (v < -5.5) {
		u   = std::exp(v+1.0);
		if (u < 1e-10) return 0.0;
		ue  = std::exp(-1/u);
		us  = std::sqrt(u);
		denlan = 0.3989422803*(ue/us)*(1+(a1[0]+(a1[1]+a1[2]*u)*u)*u);
	} else if(v < -1) {
		u   = std::exp(-v-1);
		denlan = std::exp(-u)*std::sqrt(u)*
		 (p1[0]+(p1[1]+(p1[2]+(p1[3]+p1[4]*v)*v)*v)*v)/
		 (q1[0]+(q1[1]+(q1[2]+(q1[3]+q1[4]*v)*v)*v)*v);
	} else if(v < 1) {
		denlan = (p2[0]+(p2[1]+(p2[2]+(p2[3]+p2[4]*v)*v)*v)*v)/
		 (q2[0]+(q2[1]+(q2[2]+(q2[3]+q2[4]*v)*v)*v)*v);
	} else if(v < 5) {
		denlan = (p3[0]+(p3[1]+(p3[2]+(p3[3]+p3[4]*v)*v)*v)*v)/
		 (q3[0]+(q3[1]+(q3[2]+(q3[3]+q3[4]*v)*v)*v)*v);
	} else if(v < 12) {
		u   = 1/v;
		denlan = u*u*(p4[0]+(p4[1]+(p4[2]+(p4[3]+p4[4]*u)*u)*u)*u)/
		 (q4[0]+(q4[1]+(q4[2]+(q4[3]+q4[4]*u)*u)*u)*u);
	} else if(v < 50) {
		u   = 1/v;
		denlan = u*u*(p5[0]+(p5[1]+(p5[2]+(p5[3]+p5[4]*u)*u)*u)*u)/
		 (q5[0]+(q5[1]+(q5[2]+(q5[3]+q5[4]*u)*u)*u)*u);
	} else if(v < 300) {
		u   = 1/v;
		denlan = u*u*(p6[0]+(p6[1]+(p6[2]+(p6[3]+p6[4]*u)*u)*u)*u)/
		 (q6[0]+(q6[1]+(q6[2]+(q6[3]+q6[4]*u)*u)*u)*u);
	} else {
		u   = 1/(v-v*std::log(v)/(v+1));
		denlan = u*u*(1+(a2[0]+a2[1]*u)*u);
	}
	return denlan/xi;
}

