//#include "loadfirst.h"

const char* prodtag(const char* prod){
	if(prod=="ggH")	return "";
	else if(prod=="qqH") return "_vbf";
	else if(prod=="VH")	return "_vh";
	else return "";
}

void check1Dspline(const char* chan = "4e", const char* prod = "ggH",int mh = 125, float sigma=0.004, int _debug=1){
	//const char* newrootfolder = "/afs/cern.ch/user/w/wahung/work/public/CombineLimitDbkgkin/CMSSW_7_4_7/src/HiggsAnalysis/CombinedLimit/signal_preparation/workspace125_onshell";
	//const char* rootbase = Form("hzz%s_13TeV.input_func%s_1jet_cat0",chan,prodtag(prod));
	//const char* f1dname = Form("%s/%s.root",newrootfolder,rootbase); 
	//
	
	const char* f1dname = Form("../spline_prod/spline_WS_withCat_quad9/%s_spline_WS_cat0.root",chan); 
	if(_debug)  cout<<f1dname<<endl;
	TFile *f1d = new TFile(f1dname,"read");
	RooWorkspace* w1d = (RooWorkspace*) f1d->Get("w");
//	RooRealFlooredSumPdf* sp1d = (RooRealFlooredSumPdf*) w1d->pdf(prod);
	RooRealFlooredSumPdf* sp1d = (RooRealFlooredSumPdf*) w1d->pdf(Form("r1_%s_untagged",prod));
	if(_debug)  cout<<"Done with reading ggH"<<endl;

	const char* f2dname = Form("/eos/user/w/wahung/Mass_Width_Measurement/171217_refit/mass/spline_WS_withCat_quad9/%s_spline_WS_cat0.root",chan);
	TFile *f2d = new TFile(f2dname,"read");
	RooWorkspace *w2d = (RooWorkspace*) f2d->Get("w");
	w2d->exportToCint("w2d");
	RooRealFlooredSumPdf* sp2d = (RooRealFlooredSumPdf*) w2d->obj("r1_ggH");
	if(_debug)  cout<<"Done with reading r1_ggH"<<endl;
	w1d->var("sigma_pole")->setVal(0.00407);
	w1d->var("sigma_pole")->setConstant(true);
	if(_debug)  cout<<"Done with setting sigma"<<endl;
	w1d->var("mean_pole")->setVal(125);
	w1d->var("mean_pole")->setConstant(true);
        if(_debug)  cout<<"Done with setting mh"<<endl;	
	RooRealVar* m4l = (RooRealVar*) w1d->var("mreco");
        if(_debug)  cout<<"Done with getting mreco"<<endl;	

	w2d->var("sigma_pole")->setVal(0.00407);
        w2d->var("sigma_pole")->setConstant(true);
	w2d->var("mean_pole")->setVal(125);
        w2d->var("mean_pole")->setConstant(true);

	cout<<sp1d<<" "<<sp2d<<endl;
	RooPlot* plot = (RooPlot*) m4l->frame();	
	sp1d->plotOn(plot,LineColor(kRed),LineWidth(2));
	cout<<sp1d<<" "<<sp2d<<endl;
	sp2d->plotOn(plot,LineColor(kBlack),LineWidth(2));
        if(_debug)  cout<<"Done plotting"<<endl;
	TCanvas *c = new TCanvas("c","c");
	c->cd();
	plot->Draw();
	c->SaveAs(Form("~/www/checkspline_mh%d_sigma%0.4f_%s_%s.png",mh,sigma,chan,prod));


}

