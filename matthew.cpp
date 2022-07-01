#include <iostream>
#include <filesystem>
#include <set>
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TROOT.h"
#include "TF1.h"
int main(int argc, char* argv[]){
  TCanvas* dark_sub_plot = new TCanvas("dark_sub_plot", "pure traces", 2000, 2000);
  TCanvas* ratios_plot = new TCanvas("ratios_plot", "pure ratios", 2000, 2000);
  TMultiGraph* pure_traces = new TMultiGraph();
  TMultiGraph* ratios = new TMultiGraph();
  TGraph* grad = new TGraph();
  TGraph* offs = new TGraph();
  int counter = 0;
  std::vector<double>* led_value_init;
  std::vector<double>* dark_value_init;
  std::set<std::filesystem::path> sorted_names;
  for (const auto& dir_entry : std::filesystem::directory_iterator(argv[2])){
    if (dir_entry.path().extension() == ".root"){
      sorted_names.insert(dir_entry.path());
    }
  }
  for (const auto& dir_entry : sorted_names){
    
    std::cout << "opening: " << dir_entry << "\n";
    TFile c_file(std::string(dir_entry).c_str());
    TTree* led_tree = (TTree*) c_file.Get(argv[1]);
    //    TTree* dark_tree = (TTree*) c_file.Get("Dark");
    std::vector<double>* led_value = 0;
    std::vector<double>* dark_value = 0;
    std::vector<double>* wavelength = 0;
    led_tree->SetBranchAddress("value", &led_value);
    // dark_tree->SetBranchAddress("value", &dark_value);    
    led_tree->SetBranchAddress("wavelength", &wavelength);
    
    led_tree->GetEntry(0);
    // dark_tree->GetEntry(atoi(argv[3]));
    if (counter == 0){
      led_value_init = led_value;
      // dark_value_init = dark_value;
    }
    
    
    int N = led_value->size();
    
    TGraph* dark_sub = new TGraph(N);
    TGraph* ratio = new TGraph(N);
    for (auto i = 0; i < N; ++i){
      dark_sub->SetPoint(i, wavelength->at(i), led_value->at(i) - led_value->at(500) ); // - dark_value->at(i));
      if (wavelength->at(i) < 1000){
    ratio->SetPoint(i, wavelength->at(i), ( led_value->at(i) - led_value->at(500) /* - dark_value->at(i)*/ ) / (led_value_init->at(i) - led_value_init->at(500) /* - dark_value_init->at(i)*/ ) );
      }
      else { ratio->SetPoint(i, wavelength->at(i), 0);}
      
    }
    TF1* lin = new TF1("lin", "pol1", 270, 290);
    
    ratio->Fit(lin, "+R");
    grad->AddPoint(counter, lin->GetParameter(1));
    offs->AddPoint(counter, lin->GetParameter(0));    
    // dark_sub->SetLineColor();
    //    ratios->SetLineColor(
    pure_traces->Add(dark_sub);
    ratios->Add(ratio);
    ++counter;
  }
  grad->SaveAs( (std::string(argv[2]) + "/plots/fit_grad_" + std::string(argv[1]) + ".root").c_str() );
  offs->SaveAs( (std::string(argv[2]) + "/plots/fit_offs_" + std::string(argv[1]) + ".root").c_str() );
  
  dark_sub_plot->cd();
  pure_traces->Draw("AP PMC PLC PFC L");
  dark_sub_plot->SaveAs( (std::string(argv[2]) + "/plots/pure_traces" + std::string(argv[1]) + ".root").c_str() );
  
  ratios_plot->cd();
  ratios->Draw("AP PMC PLC PFC L");
  ratios_plot->SaveAs( (std::string(argv[2]) + "/plots/ratio" + std::string(argv[1]) + ".root").c_str() );
  
  return 0;
}
