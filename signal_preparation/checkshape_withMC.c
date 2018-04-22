using namespace std;
using namespace RooFit;

void docheckshape(int cat, float mh, float sig,TString ch){
	int debug=0;	

	TFile *fnew = new TFile(Form("workspace125_onshell/hzz%s_13TeV.input_func_0jet_cat%d.root",ch.Data(),cat));
	if(debug) cout<<"Done with fnew"<<endl;
	TFile *fold = new TFile(Form("withVHreso/hzz%s_13TeV.input_func_0jet_cat%d.root",ch.Data(),cat));
	if(debug)  cout<<"Done with fold"<<endl;
	TFile *fmc = new TFile(Form("~/work/public/CMSSW_8_0_26_patch1/src/MC_check_from_Heshy/ggTo%s_0PMH125_MCFM701_xrd.root",ch.Data()));
        if(debug)  cout<<"Done with fmc"<<endl;

	RooWorkspace* wold = (RooWorkspace*) fold->Get("w");
	RooWorkspace* wnew = (RooWorkspace*) fnew->Get("w");
	RooWorkspace* wmc = (RooWorkspace*) fmc->Get("w");

	RooRealVar* mean = (RooRealVar*) wold->var("mean_pole");
	RooRealVar* sigma = (RooRealVar*) wold->var("sigma_pole");
        RooRealVar* mean2 = (RooRealVar*) wnew->var("mean_pole");
        RooRealVar* sigma2 = (RooRealVar*) wnew->var("sigma_pole");
	RooRealVar* mreco = (RooRealVar*) wold->var("mreco");
	RooRealVar* mreco2 = (RooRealVar*) wnew->var("mreco");
	
	mean2->setVal(mh);
	sigma2->setVal(sig);
//        mean2->setVal(130);
//        sigma2->setVal(0.004);

	RooPlot *plot = (RooPlot*) mreco2->frame();
	int nbin = 800;	
	RooExtendPdf* gghold = (RooExtendPdf*) wnew->obj(Form("pdf_x%s%d",ch.Data(),cat));
        if(debug)  cout<<"Done with gghold"<<endl;
	TH1F* hold = (TH1F*) gghold->createHistogram("mreco",nbin);
	double normold = hold->Integral();
	hold->Scale(1/normold);
	hold->SetLineColor(1);
	hold->SetLineWidth(2);
	hold->SetAxisRange(105,140);
	RooRelBW1* gghnew = (RooRelBW1*) wnew->obj(Form("BR_ggHNominal_%s_iCat%d_0jet",ch.Data(),cat));
        if(debug)  cout<<"Done with gghnew"<<endl;
	TH1F* hnew = (TH1F*) gghnew->createHistogram("mreco", nbin);
	double normnew = hnew->Integral();
        hnew->Scale(1/normnew);
	hnew->SetLineColor(2);
	hnew->SetLineWidth(2);
        hnew->SetAxisRange(105,140);
	TH1F* hmc = (TH1F*) wmc->obj("h");
	double normmc = hmc->Integral();
        hmc->Scale(1/normmc);	
        if(debug)  cout<<"Done with hmc"<<endl;

	hmc->SetLineColor(4);
        hmc->SetLineWidth(2);
	
	hratio = (TH1F*) hnew->Clone();
	cout<<hnew->GetXaxis()->GetNbins()<<endl;
	cout<<hmc->GetXaxis()->GetNbins()<<endl;
	hratio->Divide(hmc);
	hratio->SetMarkerStyle(5);
	hratio->SetLineColor(1);
	hratio->SetAxisRange(105,140);

	
	TCanvas* c1  = new TCanvas("","",1000,800);
	c1->Divide(1,2);
//	gghold->plotOn(plot,RooFit::LineColor(1));
//	gghnew->plotOn(plot,RooFit::LineColor(2));
	c1->cd(1);
	hnew->Draw("hist");
//	cout<<hold<<endl;
//	hold->DrawNormalized("same hist");
	hmc->Draw("same hist");
//	plot->Draw();
	c1->cd(2);	
	hratio->Draw("p");
	gPad->SetLogy();
	gPad->Update();
	c1->SaveAs(Form("~/www/checkbrshape/check_mreco_%s%d_%.1f_%.1f_withMC_withratio.png",ch.Data(),cat,mh,sig));
	c1->SaveAs(Form("~/www/checkbrshape/check_mreco_%s%d_%.1f_%.1f_withMC_withratio.pdf",ch.Data(),cat,mh,sig));
}


void checkshape_withMC(){
	vector<TString> chans={"4e","4mu","2e2mu"};
	float mh = 125.; float sig = 0.004;
	for(int cat = 0; cat<1; cat++){
	  for(int ch=0; ch<3; ch++){
  	      docheckshape(cat,mh,sig,chans[ch]);
    }
  }
}
