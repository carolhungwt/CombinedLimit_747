using namespace std;
using namespace RooFit;

void docheckshape(int cat, float mh, float sig){
	TFile *fnew = new TFile(Form("workspace125_onshell/hzz4e_13TeV.input_func_0jet_cat%d.root",cat));
	TFile *fold = new TFile(Form("withVHreso/hzz4e_13TeV.input_func_0jet_cat%d.root",cat));
	RooWorkspace* wold = (RooWorkspace*) fold->Get("w");
	RooWorkspace* wnew = (RooWorkspace*) fnew->Get("w");
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
	
	RooExtendPdf* gghold = (RooExtendPdf*) wnew->obj(Form("pdf_x4e%d",cat));
	TH1F* hold = (TH1F*) gghold->createHistogram("mreco",40);
	hold->SetLineColor(1);
	hold->SetLineWidth(2);
	RooBreitWigner* gghnew = (RooBreitWigner*) wnew->obj(Form("BR_ggHNominal_4e_iCat%d_0jet",cat));
	TH1F* hnew = (TH1F*) gghnew->createHistogram("mreco", 40);
	hnew->SetLineColor(2);
	hnew->SetLineWidth(2);
	
	TCanvas* c1  = new TCanvas("","",500,400);
//	gghold->plotOn(plot,RooFit::LineColor(1));
//	gghnew->plotOn(plot,RooFit::LineColor(2));
	hnew->DrawNormalized("hist");
//	cout<<hold<<endl;
	hold->DrawNormalized("same hist");
//	plot->Draw();
	c1->SaveAs(Form("~/www/checkbrshape/check_mreco_4e%d_%.1f_%.5f.png",cat,mh,sig));
	c1->SaveAs(Form("~/www/checkbrshape/check_mreco_4e%d_%.1f_%.5f.pdf",cat,mh,sig));
}


void checkshape(){
	for(int cat = 0; cat<1; cat++){
	  for(float mh=120; mh<160; mh+=5.){
	    for(float sig=0.004; sig<6; sig+=1.){
  	      docheckshape(cat,mh,sig);
      }
    }
  }
}
