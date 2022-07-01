{
	//#include <algorithm>
	
	/*
	TFile* may = TFile::Open("may.root");
	TCanvas* cmay = (TCanvas*)may->Get("cc");
	TGraphErrors* may_raw = (TGraphErrors*)cmay->GetListOfPrimitives()->At(5);
	TGraphErrors* may_simple = (TGraphErrors*)cmay->GetListOfPrimitives()->At(1);
	TGraphErrors* may_complex = (TGraphErrors*)cmay->GetListOfPrimitives()->At(3);
	
	TFile* june = TFile::Open("june.root");
	TCanvas* cjune = (TCanvas*)june->Get("cc");
	TGraphErrors* june_raw = (TGraphErrors*)cjune->GetListOfPrimitives()->At(5);
	TGraphErrors* june_simple = (TGraphErrors*)cjune->GetListOfPrimitives()->At(1);
	TGraphErrors* june_complex = (TGraphErrors*)cjune->GetListOfPrimitives()->At(3);
	
	june_raw->SetMarkerStyle(30);
	june_simple->SetMarkerStyle(30);
	june_complex->SetMarkerStyle(30);
	*/
	
	/*
	TFile* june2_A = TFile::Open("june2_A.root");
	TCanvas* cjune2_A = (TCanvas*)june2_A->Get("cc");
	TGraphErrors* june_raw2_A = (TGraphErrors*)cjune2_A->GetListOfPrimitives()->At(5);
	TGraphErrors* june_simple2_A = (TGraphErrors*)cjune2_A->GetListOfPrimitives()->At(1);
	TGraphErrors* june_complex2_A = (TGraphErrors*)cjune2_A->GetListOfPrimitives()->At(3);
	
	june_raw2_A->SetMarkerStyle(30);
	june_simple2_A->SetMarkerStyle(30);
	june_complex2_A->SetMarkerStyle(30);
	
	TFile* june2_B = TFile::Open("june2_B.root");
	TCanvas* cjune2_B = (TCanvas*)june2_B->Get("cc");
	TGraphErrors* june_raw2_B = (TGraphErrors*)cjune2_B->GetListOfPrimitives()->At(5);
	TGraphErrors* june_simple2_B = (TGraphErrors*)cjune2_B->GetListOfPrimitives()->At(1);
	TGraphErrors* june_complex2_B = (TGraphErrors*)cjune2_B->GetListOfPrimitives()->At(3);
	
	june_raw2_B->SetMarkerStyle(30);
	june_raw2_B->SetMarkerColor(kMagenta);
	june_simple2_B->SetMarkerStyle(30);
	june_simple2_B->SetMarkerColor(kAzure+1);
	june_complex2_B->SetMarkerStyle(30);
	june_complex2_B->SetMarkerColor(kViolet+5);
	*/
	
	TFile* june3_A = TFile::Open("june3_A.root");
	TCanvas* cjune3_A = (TCanvas*)june3_A->Get("cc");
	TGraphErrors* june_raw3_A = (TGraphErrors*)cjune3_A->GetListOfPrimitives()->At(5);
	TGraphErrors* june_simple3_A = (TGraphErrors*)cjune3_A->GetListOfPrimitives()->At(1);
	TGraphErrors* june_complex3_A = (TGraphErrors*)cjune3_A->GetListOfPrimitives()->At(3);
	
	june_raw3_A->SetMarkerStyle(30);
	june_simple3_A->SetMarkerStyle(30);
	june_complex3_A->SetMarkerStyle(30);
	
	TFile* june3_B = TFile::Open("june3_B.root");
	TCanvas* cjune3_B = (TCanvas*)june3_B->Get("cc");
	TGraphErrors* june_raw3_B = (TGraphErrors*)cjune3_B->GetListOfPrimitives()->At(5);
	TGraphErrors* june_simple3_B = (TGraphErrors*)cjune3_B->GetListOfPrimitives()->At(1);
	TGraphErrors* june_complex3_B = (TGraphErrors*)cjune3_B->GetListOfPrimitives()->At(3);
	
	june_raw3_B->SetMarkerStyle(28);
	june_raw3_B->SetMarkerColor(kMagenta);
	june_simple3_B->SetMarkerStyle(28);
	june_simple3_B->SetMarkerColor(kAzure+1);
	june_complex3_B->SetMarkerStyle(28);
	june_complex3_B->SetMarkerColor(kViolet+5);
	
	/*
	TFile* feb_A = TFile::Open("feb_A.root");
	TCanvas* cfeb_A = (TCanvas*)feb_A->Get("cc");
	TGraphErrors* feb_raw_A = (TGraphErrors*)cfeb_A->GetListOfPrimitives()->At(5);
	TGraphErrors* feb_simple_A = (TGraphErrors*)cfeb_A->GetListOfPrimitives()->At(1);
	TGraphErrors* feb_complex_A = (TGraphErrors*)cfeb_A->GetListOfPrimitives()->At(3);
	
	feb_raw_A->SetMarkerStyle(3);
	feb_simple_A->SetMarkerStyle(3);
	feb_complex_A->SetMarkerStyle(3);
	
	TFile* feb_B = TFile::Open("feb_B.root");
	TCanvas* cfeb_B = (TCanvas*)feb_B->Get("cc");
	TGraphErrors* feb_raw_B = (TGraphErrors*)cfeb_B->GetListOfPrimitives()->At(5);
	TGraphErrors* feb_simple_B = (TGraphErrors*)cfeb_B->GetListOfPrimitives()->At(1);
	TGraphErrors* feb_complex_B = (TGraphErrors*)cfeb_B->GetListOfPrimitives()->At(3);
	
	feb_raw_B->SetMarkerStyle(3);
	feb_raw_B->SetMarkerColor(kMagenta);
	feb_simple_B->SetMarkerStyle(3);
	feb_simple_B->SetMarkerColor(kAzure+1);
	feb_complex_B->SetMarkerStyle(3);
	feb_complex_B->SetMarkerColor(kViolet+5);
	*/
	
	// fit and print curves
	std::map<std::string, TGraphErrors*> curves{{"june3_A_raw",june_raw3_A},
	                                            {"june3_A_simple",june_simple3_A},
	                                            {"june3_A_complex",june_complex3_A},
	                                            {"june3_B_raw",june_raw3_B},
	                                            {"june3_B_simple",june_simple3_B},
	                                            {"june3_B_complex",june_complex3_B}};
	TCanvas c0;
	for(auto&& acurve : curves){
		c0.Clear();
		acurve.second->Draw("AP");
		c0.Modified();
		c0.Update();
		gSystem->ProcessEvents();
		gPad->WaitPrimitive();
		std::cout<<"fitting "<<acurve.first<<" from 0 to "<<acurve.second->GetX()[acurve.second->GetN()-1]<<std::endl;
		TF1* nextfunc = new TF1(acurve.first.c_str(),"pol6",0,acurve.second->GetX()[acurve.second->GetN()-1]);
		TFitResultPtr fptr = acurve.second->Fit(acurve.first.c_str(),"RSQ");
		std::cout<<"Fit status was: "<<fptr->Status()<<std::endl;
		c0.Modified();
		c0.Update();
		gSystem->ProcessEvents();
		gPad->WaitPrimitive();
		std::cout<<"curve "<<acurve.first<<" fit parameters: [";
		for(int i=0; i<nextfunc->GetNpar(); ++i){
			if(i>0) std::cout<<", ";
			std::cout<<nextfunc->GetParameters()[i];
		}
		std::cout<<"]"<<std::endl;
	}
	
	// draw all curves
	TCanvas c1; c1.cd();
	TMultiGraph mg("mg","mg");
	/*
	mg.Add(may_raw);
	mg.Add(may_simple);
	mg.Add(may_complex);
	mg.Add(june_raw);
	mg.Add(june_simple);
	mg.Add(june_complex);
	*/
	
	/*
	mg.Add(june_raw2_A);
	mg.Add(june_simple2_A);
	mg.Add(june_complex2_A);
	*/
	
	/*
	mg.Add(june_raw2_B);
	mg.Add(june_simple2_B);
	mg.Add(june_complex2_B);
	*/
	
	/*
	mg.Add(feb_raw_A);
	mg.Add(feb_simple_A);
	mg.Add(feb_complex_A);
	*/
	
	/*
	mg.Add(feb_raw_B);
	mg.Add(feb_simple_B);
	mg.Add(feb_complex_B);
	*/
	
	mg.Add(june_raw3_A);
	mg.Add(june_simple3_A);
	mg.Add(june_complex3_A);
	
	mg.Add(june_raw3_B);
	mg.Add(june_simple3_B);
	mg.Add(june_complex3_B);
	
	
	mg.Draw("ALP");
	mg.GetYaxis()->SetRangeUser(-0.01,0.3);
	c1.Modified();
	c1.Update();
	gSystem->ProcessEvents();
	
	return 0;
	
//	std::cout<<"june2_conc, june2_A_peaks, june2_B_peaks, feb_concs, feb_A_peaks, feb_B_peaks"<<"\n";
//	for(int i=0; i<feb_raw_B->GetN(); ++i){
//		if(i<june_raw2_A->GetN()){
//			std::cout<<june_raw2_A->GetX()[i]<<","<<june_raw2_A->GetY()[i]<<","<<june_raw2_B->GetY()[i]
//			         <<","<<feb_raw_A->GetX()[i]<<","<<feb_raw_A->GetY()[i]<<","<<feb_raw_B->GetY()[i]<<"\n";
//		} else {
//			std::cout<<",,,"
//			         <<","<<feb_raw_A->GetX()[i]<<","<<feb_raw_A->GetY()[i]<<","<<feb_raw_B->GetY()[i]<<"\n";
//		}
//	}
//	std::cout<<std::endl;
//	
//	//TGraph* gA = june_raw2_B;
//	TGraph* gA = feb_raw_B;
//	TGraph* gB = june_raw3_B;
//	
//	TGraph* adata = new TGraph(gA->GetN(),gA->GetX(),gA->GetY());
//	TGraph* bdata = new TGraph(gB->GetN(),gB->GetX(),gB->GetY());
//	int npoints = gB->GetN();
//	TGraph* ratio = new TGraph(npoints);
//	TF1* junefit = new TF1("afit","pol3",0,adata->GetX()[adata->GetN()-1]);
//	junefit->SetLineColor(kRed);
//	adata->Fit("afit","RQ");
//	TF1* febfit = new TF1("bfit","pol3",0,bdata->GetX()[gB->GetN()-1]);
//	febfit->SetLineColor(kBlue);
//	bdata->Fit("bfit","RQ");
//	int nextpt=0;
//	for(int i=0; i<npoints; ++i){
//		double nextdiff = (double(i)/double(npoints))*adata->GetY()[npoints-1];
//		double aconc = junefit->GetX(nextdiff,0.0,adata->GetY()[npoints-1]);
//		double bconc = febfit->GetX(nextdiff,0.0,bdata->GetY()[npoints-1]);
//		//std::cout<<"ratio = "<<aconc<<"/"<<bconc<<" = "<<aconc/bconc<<" at "<<nextdiff<<std::endl;
//		aconc = std::max(aconc,0.);
//		bconc = std::max(bconc,0.);
//		if(bconc==0 || aconc==0) continue;
//		ratio->SetPoint(nextpt,nextdiff,aconc/bconc);
//		nextpt++;
//	}
//	ratio->Set(nextpt);
//	TCanvas c2("c2","c2",1024,800);
//	ratio->Draw("A*");
//	
//	TCanvas c3("c3","c3",1024,800);
//	adata->Draw("A*");
//	bdata->Draw("same *");
	
	
}
