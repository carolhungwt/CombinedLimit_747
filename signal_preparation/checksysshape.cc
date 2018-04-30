#include "HiggsAnalysis/CombinedLimit/interface/RooNCSpline_2D_fast.h"
using namespace RooFit;
using namespace std;


void makeplot(TString channel, int cat, TString vbfcate){
	TFile *f = new TFile(Form("%s_rpdfWS_cat%d_%s_0jet_.root",channel.Data(),cat,vbfcate.Data()),"read");
	RooWorkspace *w = (RooWorkspace*) f->Get("w");
	RooRealVar* mh = (RooRealVar*) w->var("mean_pole");
	RooRealVar* sigma = (RooRealVar*) w->var("sigma_pole");
	RooRealVar* m4l = (RooRealVar*) w->var("mreco");

	mh->setVal(125.);
	sigma->setVal(0.004);

	TCanvas *c1 = new TCanvas("","",1000,800);

	TString pdfname = "ggH";
	if(vbfcate=="vbf")  pdfname="qqH";
	else if(vbfcate=="vh")  pdfname="VH";

	RooPlot *rooplot = m4l->frame();
	RooAbsReal* nominal = (RooAbsReal*) w->obj(Form("r1_%s%d_%s",channel.Data(),cat,pdfname.Data()));
//	nominal->plotOn(rooplot,LineColor(1));
	RooAbsReal* resdn = (RooAbsReal*) w->obj(Form("r1_%s%d_%s_ResDown",channel.Data(),cat,pdfname.Data()));
//	resdn->plotOn(rooplot,LineColor(2));
	RooAbsReal* resup = (RooAbsReal*) w->obj(Form("r1_%s%d_%s_ResUp",channel.Data(),cat,pdfname.Data()));
//	resup->plotOn(rooplot,LineColor(6));
	RooAbsReal* scaledn = (RooAbsReal*) w->obj(Form("r1_%s%d_%s_ScaleDown",channel.Data(),cat,pdfname.Data()));
	scaledn->plotOn(rooplot,LineColor(4));
	RooAbsReal* scaleup = (RooAbsReal*) w->obj(Form("r1_%s%d_%s_ScaleUp",channel.Data(),cat,pdfname.Data()));
	scaleup->plotOn(rooplot,LineColor(9));
	
	rooplot->Draw();
	c1->SaveAs("~/www/check_sys_"+pdfname+".png");
	c1->SaveAs("~/www/check_sys_"+pdfname+".pdf");
}

void checksysshape(){
	makeplot("4e", 1, "ggH");
	makeplot("4e", 1, "vbf");
	makeplot("4e", 1, "vh");
}
