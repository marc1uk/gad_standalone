#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TApplication.h"
#include "TSystem.h"
#include "TCanvas.h"
#include <thread>
#include <chrono>

int GetDarkEntryNum(TTree* ledtree, TTree* darktree);
TGraph* MakeTGraph(TTree* tl, TTree* td, int darkentry);
void NormaliseTGraph(TGraph* g);
TGraph* GetFitGraph(TF1* pure_fit, TGraph* g, std::vector<double>& pars, bool reinit, bool refit=false);
TF1* PureScaledPlusExtras();

TGraph* dark_subtracted_pure=nullptr;
std::string purefile="";
std::string pureref_275A = "../GDConcMeasure/pureDarkSubtracted_275_A_18mar2023.root";
std::string pureref_275B = "../GDConcMeasure/pureDarkSubtracted_275_B_18mar2023.root";

int main(){
	
	TApplication tapp("app",0,0,0);
	
	TFile* _file0 = TFile::Open("../GDConcMeasure/data/2023/03/00150_18Mar2023Calibration_00028.root");
	TFile* _file1 = TFile::Open("../GDConcMeasure/data/2023/03/00150_18Mar2023Calibration_00029.root");
	TFile* _file2 = TFile::Open("../GDConcMeasure/data/2023/03/00150_18Mar2023Calibration_00030.root");
	
	TTree* t1a=(TTree*)_file0->Get("275_A");
	TTree* t1b=(TTree*)_file0->Get("275_B");
	TTree* t1d=(TTree*)_file0->Get("Dark");
	
	TTree* t2a=(TTree*)_file1->Get("275_A");
	TTree* t2b=(TTree*)_file1->Get("275_B");
	TTree* t2d=(TTree*)_file1->Get("Dark");
	
	TTree* t3a=(TTree*)_file2->Get("275_A");
	TTree* t3b=(TTree*)_file2->Get("275_B");
	TTree* t3d=(TTree*)_file2->Get("Dark");
	
	// will all be the same so only need to do it for one
	int darkentry_A = GetDarkEntryNum(t1a, t1d);
	int darkentry_B = GetDarkEntryNum(t1b, t1d);
	
	// generate dark subtracted tgraphs
	TGraph* g1a = MakeTGraph(t1a, t1d, darkentry_A);
	TGraph* g1b = MakeTGraph(t1b, t1d, darkentry_B);
	TGraph* g2a = MakeTGraph(t2a, t2d, darkentry_A);
	TGraph* g2b = MakeTGraph(t2b, t2d, darkentry_B);
	TGraph* g3a = MakeTGraph(t3a, t3d, darkentry_A);
	TGraph* g3b = MakeTGraph(t3b, t3d, darkentry_B);
	
	// name them
	g1a->SetName("275A_1"); g1a->SetTitle("275A_1");
	g1b->SetName("275B_1"); g1b->SetTitle("275B_1");
	g2a->SetName("275A_2"); g2a->SetTitle("275A_2");
	g2b->SetName("275B_2"); g2b->SetTitle("275B_2");
	g3a->SetName("275A_3"); g3a->SetTitle("275A_3");
	g3b->SetName("275B_3"); g3b->SetTitle("275B_3");
	
	// colour them
	g1a->SetMarkerColor(kRed);
	g1b->SetMarkerColor(kRed);
	g2a->SetMarkerColor(kMagenta);
	g2b->SetMarkerColor(kMagenta);
	g3a->SetMarkerColor(kBlue);
	g3b->SetMarkerColor(kBlue);
	
	// distinguish the plots
	g1b->SetMarkerStyle(2);
	g2b->SetMarkerStyle(2);
	g3b->SetMarkerStyle(2);
	
	g1a->SetMarkerStyle(5);
	g2a->SetMarkerStyle(5);
	g3a->SetMarkerStyle(5);
	
	g1a->SetLineWidth(0);
	g1b->SetLineWidth(0);
	g2a->SetLineWidth(0);
	g2b->SetLineWidth(0);
	g3a->SetLineWidth(0);
	g3b->SetLineWidth(0);
	
	// scale them to compare shape only
	/*
	NormaliseTGraph(g1a);
	NormaliseTGraph(g1b);
	NormaliseTGraph(g2a);
	NormaliseTGraph(g2b);
	NormaliseTGraph(g3a);
	NormaliseTGraph(g3b);
	*/
	
	TMultiGraph mg("mg","mg");
	mg.Add(g1a);
	mg.Add(g1b);
	//mg.Add(g2a);
	//mg.Add(g2b);
	mg.Add(g3a);
	mg.Add(g3b);
	
	// get the fit parameters
	TFile* f_fp = TFile::Open("monitor_purefits.root");
	TTree* tpurea = (TTree*)f_fp->Get("LED275A");
	TTree* tpureb = (TTree*)f_fp->Get("LED275B");
	std::vector<double> fit_pars;
	std::vector<double>* fit_parsp = &fit_pars;
	tpurea->SetBranchAddress("purefits",&fit_parsp);
	tpureb->SetBranchAddress("purefits",&fit_parsp);
	
	// Generate TGraphs from the pure fits and overlay them to check it was ok
	TF1* pure_fit = PureScaledPlusExtras();
	
	purefile=pureref_275A;
	tpurea->GetEntry(28); // 1
	TGraph* g1af = GetFitGraph(pure_fit, g1a, fit_pars, true);
	tpurea->GetEntry(29); // 2
	TGraph* g2af = GetFitGraph(pure_fit, g2a, fit_pars, false);
	tpurea->GetEntry(30); // 3
	TGraph* g3af = GetFitGraph(pure_fit, g3a, fit_pars, false);
	purefile = pureref_275B;
	tpureb->GetEntry(28); // 1
	TGraph* g1bf = GetFitGraph(pure_fit, g1b, fit_pars, true);
	tpureb->GetEntry(29); // 2
	TGraph* g2bf = GetFitGraph(pure_fit, g2b, fit_pars, false);
	tpureb->GetEntry(30); // 2
	TGraph* g3bf = GetFitGraph(pure_fit, g3b, fit_pars, false);
	
	// name them
	g1af->SetName("275A_1_f"); g1af->SetTitle("275A_1_f");
	g1bf->SetName("275B_1_f"); g1bf->SetTitle("275B_1_f");
	g2af->SetName("275A_2_f"); g2af->SetTitle("275A_2_f");
	g2bf->SetName("275B_2_f"); g2bf->SetTitle("275B_2_f");
	g3af->SetName("275A_3_f"); g3af->SetTitle("275A_3_f");
	g3bf->SetName("275B_3_f"); g3bf->SetTitle("275B_3_f");
	
	// configure line styles
	g1af->SetLineColor(kRed);
	g1bf->SetLineColor(kRed);
	g2af->SetLineColor(kMagenta);
	g2bf->SetLineColor(kMagenta);
	g3af->SetLineColor(kBlue);
	g3bf->SetLineColor(kBlue);
	
	g1af->SetLineStyle(2);
	g2af->SetLineStyle(2);
	g3af->SetLineStyle(2);
	
	// collect into a mutigraph
	mg.Add(g1af);
	mg.Add(g1bf);
	//mg.Add(g2af);
	//mg.Add(g2bf);
	mg.Add(g3af);
	mg.Add(g3bf);
	
	TCanvas c1("c1","c1");
	mg.Draw("ALP");
	c1.BuildLegend();
	c1.Modified();
	c1.Update();
	
	// wait till user closes tgraph
	while(gROOT->FindObject("c1")){
		gSystem->ProcessEvents();
		std::this_thread::sleep_for(std::chrono::milliseconds(500));
	}
	
	return 0;
}


TGraph* MakeTGraph(TTree* tl, TTree* td, int darkentry){
	
	std::vector<double> wls;
	std::vector<double>* wlsp = &wls;
	
	std::vector<double> vals;
	std::vector<double>* valsp = &vals;
	
	std::vector<double> darks;
	std::vector<double>* darksp = &darks;
	
	tl->SetBranchAddress("wavelength",&wlsp);
	tl->SetBranchAddress("value",&valsp);
	td->SetBranchAddress("value",&darksp);
	
	tl->GetEntry(0);
	td->GetEntry(darkentry);
	std::vector<double> darksub_vals(wls.size());
	for(int i=0; i<wls.size(); ++i){
		darksub_vals.at(i) = vals.at(i) - darks.at(i);
	}
	TGraph* g = new TGraph(wls.size(), wls.data(), darksub_vals.data());
	return g;
}

void NormaliseTGraph(TGraph* g){
	double max = *std::max_element(g->GetY(), g->GetY()+g->GetN());
	for(int i=0; i<g->GetN(); ++i){
		g->GetY()[i] /= max;
	}
	return;
}

int GetDarkEntryNum(TTree* ledtree, TTree* darktree){
	
	// find the dark entry number representing the last dark before the given light
	
	Short_t yr, mon, dy, hr, mn, sc;
	ledtree->SetBranchAddress("year",&yr);
	ledtree->SetBranchAddress("month",&mon);
	ledtree->SetBranchAddress("day",&dy);
	ledtree->SetBranchAddress("hour",&hr);
	ledtree->SetBranchAddress("min",&mn);
	ledtree->SetBranchAddress("sec",&sc);
	
	// get led timestamp - only one entry
	ledtree->GetEntry(0);
	struct tm ledtime;
	ledtime.tm_year = yr - 1900;
	ledtime.tm_mon = mon - 1;
	ledtime.tm_mday = dy;
	ledtime.tm_hour = hr;
	ledtime.tm_min = mn;
	ledtime.tm_sec = sc;
	time_t ledtime_t = mktime(&ledtime);
	
	// set up dark tree
	darktree->SetBranchAddress("year",&yr);
	darktree->SetBranchAddress("month",&mon);
	darktree->SetBranchAddress("day",&dy);
	darktree->SetBranchAddress("hour",&hr);
	darktree->SetBranchAddress("min",&mn);
	darktree->SetBranchAddress("sec",&sc);
	
	// loop over dark entries
	int darkentry=-1;
	for(int i=0; i<darktree->GetEntries(); ++i){
		darktree->GetEntry(i);
		
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
		if(numsecs<0 && darkentry>=0) break;
		darkentry=i;
	}
	
	return darkentry;
}

double PureFuncv2(double* x, double* par){
	if(dark_subtracted_pure==nullptr){
		TFile *_file0 = TFile::Open(purefile.c_str());
		dark_subtracted_pure = (TGraph*)_file0->Get("Graph");
	}
	// par [0] = y-scaling
	// par [1] = x-scaling
	// par [2] = x-offset
	// par [3] = y-offset
	// par [4] = linear baseline offset
	// par [5] = shoulder gaussian scaling
	// par [6] = shoulder gaussian centre, restricted to > 282nm (RH shoulder)
	// par [7] = shoulder gaussian spread
	double purepart = (par[0] * dark_subtracted_pure->Eval((par[1]*(x[0]-276))+276+par[2]));
	double linpart = (par[4] * (x[0]-276)) + par[3];
	double shoulderpart = par[5]*exp(-0.5*TMath::Sq((x[0]-282.-abs(par[6]))/par[7]));
	double retval = purepart + linpart + shoulderpart;

	return retval;
}

TF1* PureScaledPlusExtras(){
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

TGraph* GetFitGraph(TF1* pure_fit, TGraph* g, std::vector<double>& pars, bool reinit, bool refit){
	// we need to delete and re-set the pure reference when changing LEDs
	if(reinit){
		if(dark_subtracted_pure) delete dark_subtracted_pure;
		dark_subtracted_pure = nullptr;
	}
	if(refit) g->Fit(pure_fit,"RNQS");
	else pure_fit->SetParameters(pars.data());
	std::vector<double> fitvals(g->GetN());
	for(int i=0; i<g->GetN(); ++i){
		fitvals.at(i)  = pure_fit->Eval(g->GetX()[i]);
	}
	TGraph* gout = new TGraph(g->GetN(), g->GetX(), fitvals.data());
	return gout;
}
