#include "../interface/RooNCSplinePdf_1D_fast.h"
#include "TROOT.h"
#include "../interface/RooNCSplineFactory_1D.h"
#include "../interface/RooRealFlooredSumPdf.h"
#include "TFile.h"
#include "TH2F.h"
#include "RooProdPdf.h"
#include "RooConstVar.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "../interface/FastTemplateFunc.h"
#include "../interface/FastTemplate.h"
#include "TCanvas.h"
using namespace RooFit;
using std::vector;
using std::endl;
using std::cout;
using std::pair;
using namespace std;


vector<double> set_grid(RooRealVar* var);
RooNCSpline_1D_fast* make_1D_spPDF(RooRealVar& projvar, vector<double> projval, int ndim, vector<double> fixedVals,RooWorkspace* w, TString tag, int prod_cate);
RooRealFlooredSumPdf* getTmpPdf(TH2F* sigtmp, RooArgList inObsList);
void createWS(RooProdPdf* pdf,char* rootname,vector<RooRealVar*> projvars);
vector<double> convert_arr_to_vec(double vals[], int size);
vector<int> gen_ndim(vector<vector<double>> vec_vec);


void do_make_rpdfWS_withCat(TString outpath, TString tag, int cat,int quad, int prod_cate){
	clock_t tStart = clock();
	char* inputDir = "../signal_preparation";
	char* tag2 = (char*) tag.Data();
	TString rootprodtag="";
	if(prod_cate==0)	rootprodtag="";
	else if(prod_cate==1)	rootprodtag="vbf_";
	else if(prod_cate==2)	rootprodtag="vh_";
	TFile *fm = new TFile(Form("%s/workspace125_onshell/hzz%s_13TeV.input_func_%s0jet_cat%d.root",inputDir,tag2,rootprodtag.Data(),cat));
	RooWorkspace* ws = (RooWorkspace*) fm->Get("w");


	// Get all variables from workspace
	std::vector<RooRealVar*> projvars;
	RooRealVar *mean_pole = ws->var("mean_pole");
	mean_pole->setRange(122.,128.);
	RooRealVar *mreco = ws->var("mreco");
	mreco->setRange(105.,140.);
	RooRealVar *sigma_pole = ws->var("sigma_pole");
	sigma_pole->setRange(0.,5.);
	RooRealVar *r = ws->var("r");
	r->setVal(1);
	RooRealVar* dbkg = new RooRealVar("dbkgkin","dbkgkin",0.,0.,1.);

	// 0: mH, 1: m4l, 2: gammaH
	projvars.push_back(mean_pole);
	projvars.push_back(mreco);
	projvars.push_back(sigma_pole);


// Construct 1D spline 
	// set grid markers for spline construction
	vector<vector<double>> projvals;
	for(unsigned int i=0; i<projvars.size(); i++){
		projvals.push_back(set_grid(projvars[i]));
	}
	// get number of grid points
	vector<int> ndim = gen_ndim(projvals);
	int totalpts = 0;
	for(unsigned int i=0; i<ndim.size(); i++)
	{
		totalpts+=ndim[i];
	}

	TString prodtag="ggH";
	if(prod_cate==1)	prodtag="qqH";

	// Extract Template
	TFile *ftmp = new TFile("../TemplateBuilder/Normalized_spline_tmpl.root","read");
	char* tmpname = Form("ggzz_sig_%s",tag.Data());
	if(prod_cate==1)	tmpname = Form("vbf_sig_%s",tag.Data());
	TH2F* sigtmp = (TH2F*) ftmp->Get(tmpname);
	sigtmp->SetName(Form("%s_sig_%s",prodtag.Data(),tag.Data()));

	// push dbkg in the vector after getting the grids
	projvars.push_back(dbkg);
	RooArgList inObsList(*mreco,*dbkg);

	// Convert tmp into RooRealFlooredSumPdf
	RooRealFlooredSumPdf* tmppdf = getTmpPdf(sigtmp, inObsList);

	vector<double> fixedVals(2);
	for(int imh=0; imh<ndim[0]; imh++){
		for(int isigma=0; isigma<ndim[2]; isigma++){
			fixedVals[0]=projvals[0][imh];
			fixedVals[1]=projvals[2][isigma];
			//Make 1D spline for signal
			RooNCSpline_1D_fast* r1 = make_1D_spPDF(*projvars[1], projvals[1], ndim[1],fixedVals, ws,"r1",prod_cate);
			r1->Print();
			cout<<r1<<endl;
//			RooAbsReal*	r1 = (RooAbsReal*) r1fac->getFunc();
			cout<<"Done making spline"<<endl;
			RooConstVar* one = new RooConstVar("one","",1.);
			RooRealFlooredSumPdf* pdfsig = new RooRealFlooredSumPdf("pdfr1","",RooArgList(*r1),RooArgList(*one));
			char* outputroot = Form("%s/1DSpline_mH%.2f_sigma%.2f.root",outpath.Data(),fixedVals[0],fixedVals[1]);
			//Create WS and store in root file
			RooProdPdf *finalpdf = new RooProdPdf("signal","",*pdfsig,*tmppdf);
			cout<<finalpdf->GetName()<<endl;
			cout<<"Done making signal pdf"<<endl;
	
//			createWS(finalpdf,outputroot, projvars);
	TFile ftest(outputroot,"recreate");
	cout<<"Done creating root file"<<endl;
	ftest.cd();
	RooWorkspace* w = new RooWorkspace("w", "");
	w->addClassDeclImportDir("../interface/");
	w->importClassCode(RooNCSplineCore::Class(),true);
	w->importClassCode(RooNCSpline_1D_fast::Class(), true);
	w->importClassCode("FastHisto2D_f",true);
	w->importClassCode(FastHisto2DFunc_f::Class(),true);
	w->importClassCode(RooRealFlooredSumPdf::Class(),true);
	cout<<"Done importing classes"<<endl;
	w->cd();
	w->import(*r1);cout<<"done with r1"<<endl;
	w->import(*tmppdf); cout<<"done with tmppdf"<<endl;
	w->import(*finalpdf, RooFit::RecycleConflictNodes());
	cout<<"Done importing pdf"<<endl;
	for(unsigned int i=0; i<projvars.size();i++){
		w->import(*projvars[i],RooFit::RecycleConflictNodes());
	}	

	w->Print();
	ftest.WriteTObject(w);
	ftest.Close();



		}
	}

	fm->Close();
	return;
}

RooRealFlooredSumPdf* getTmpPdf(TH2F* sigtmp, RooArgList inObsList){

	const char* tmpname = sigtmp->GetName();
	FastHisto2D_f *hist2df = new FastHisto2D_f(*sigtmp,true);
	FastHisto2DFunc_f *hist2dfuncf= new FastHisto2DFunc_f(Form("Histo2DFunc_%s",tmpname),"",inObsList,*hist2df);
	RooConstVar* one = new RooConstVar("one","",1);
	RooRealFlooredSumPdf* tmppdf = new RooRealFlooredSumPdf(Form("RooRealFlooredSumPdf_%s",tmpname),"",RooArgList(*hist2dfuncf),RooArgList(*one));
	return tmppdf;
}

RooNCSpline_1D_fast* make_1D_spPDF(RooRealVar& projvar, vector<double> projval, int ndim, vector<double> fixedVals,RooWorkspace* w, TString tag, int prod_cate){
	// 1D spline only deals with bkg
	int r=0;
	RooAbsPdf* pdf;
	if(prod_cate==0){
		w->var("r")->setVal(r);
		pdf = w->pdf("ggH"); 
	}
	else if(prod_cate==1){
		w->var("r")->setVal(r);
		w->var("rvbf_ggh")->setVal(1);
		pdf = w->pdf("qqH");
	}
	else {
		w->var("r")->setVal(r);
		w->var("rvbf_ggh")->setVal(1);
		pdf = w->pdf("VH");
	}

	double temp = 0.;
	RooNCSplineFactory_1D* spFac =  new RooNCSplineFactory_1D(projvar);
	double fcn =0.;
//	vector<pair<double,double>> points;
	vector<float> fcns;
	vector<float> xvals;	
	// Fixing parameter
	TString fixedVars[] = {"mean_pole","sigma_pole"};
	for(int i=0; i<2; i++){
		w->var(fixedVars[i])->setConstant(false);
		w->var(fixedVars[i])->setVal(fixedVals[i]);
		w->var(fixedVars[i])->setConstant(true);
	}

	for(int ix=0; ix<ndim; ix++){
		temp = projval.at(ix);
		w->var("mreco")->setConstant(false);
		w->var("mreco")->setVal(temp);
		w->var("mreco")->setConstant(true);
		fcn = pdf->getVal();
//		typedef pair<float,float> tmppair;		const tmppair p(temp,fcn);
		fcns.push_back(fcn);	xvals.push_back(temp);
	}
	cout<<"Done extracting points"<<endl;
	// plot sth
	RooRealVar xvar("mreco","",125,105,140);
//	RooPlot* plot = xvar.frame();
//	TCanvas* c1;
//	c1->cd();	
	spFac->setPoints(xvals,fcns);
	cout<<"Done setting points"<<endl;
	RooNCSpline_1D_fast* spPDF = spFac->getFunc();
	cout<<"done converting to RooNCSpline_1D_fast"<<endl;
//	spPDF->SetNameTitle(tag, tag);
//	spPDF->plotOn(plot);
//	plot->Draw();
//	c1->SaveAs(Form("~/www/All1DSpline/check_r0_mH%.2f_sigma%.2f.png",fixedVals[0],fixedVals[1]));
	return spPDF;
}



void createWS(RooProdPdf *pdf,char* rootname,vector<RooRealVar*> projvars){
	cout<<rootname<<endl;
	TFile ftest(rootname,"recreate");
	cout<<"Done creating root file"<<endl;
	ftest.cd();
	RooWorkspace* w = new RooWorkspace("w", "");
	w->addClassDeclImportDir("../interface/");
	w->importClassCode(RooNCSplineCore::Class(),true);
	w->importClassCode(RooNCSpline_1D_fast::Class(), true);
//	w->importClassCode(FastHisto2D_f::Class(),true);
	w->importClassCode(FastHisto2DFunc_f::Class(),true);
	w->importClassCode(RooRealFlooredSumPdf::Class(),true);
	cout<<"Done importing classes"<<endl;
	w->import(*pdf, RooFit::RecycleConflictNodes());
	cout<<"Done importing pdf"<<endl;
	for(unsigned int i=0; i<projvars.size();i++){
		w->import(*projvars[i],RooFit::RecycleConflictNodes());
	}	

	w->Print();
	ftest.WriteTObject(w);
	ftest.Close();

}

vector<int> gen_ndim(vector<vector<double>> vec_vec){
	vector<int> ndim;
	int layer1size = vec_vec.size();
	for(int i=0; i<layer1size; i++){
		ndim.push_back(vec_vec[i].size());
		cout<<vec_vec[i].size()<<" ";
	}
	cout<<endl;
	return ndim;
}

vector<double> convert_arr_to_vec(double vals[], int size){
//	int size = sizeof(*vals)/sizeof(double);
	vector<double> vecs;
	for(int i=0; i<size; i++){
		vecs.push_back(vals[i]);
		cout<<vals[i]<<" ";
	}	cout<<endl;
	return vecs;
}

vector<double> set_grid(RooRealVar* var){
	vector<double> gridMarkers;
	TString varName = var->GetName();
	double inc = 0.2;
	int gridnum = 10;
	if (varName=="mean_pole"){
		inc = 1.;
		double min_mH = 125.;
		double max_mH = 128.;
		gridnum = (max_mH-min_mH)/inc+1;
		for(int i = 0; i<gridnum; i++){
			gridMarkers.push_back(min_mH+inc*i);
		}	
	}	else if(varName == "mreco"){
		double temp_vals[]={105.,106., 107.,108., 109., 110.,111.,112., 113., 114., 115.,116., 117.,118., 119., 120.,121., 122., 123., 124., 125., 126., 127., 128., 129.5, 131., 132.5, 134., 135.5, 137., 138.5, 140.};
		int size = sizeof(temp_vals)/sizeof(double);
		gridMarkers = convert_arr_to_vec(temp_vals, size);
	}	else if(varName == "sigma_pole"){
		double temp_vals[]={0.004};//,0.001, 0.002,.00407,0.007,0.01,0.05,0.08,0.1};//,0.15, 0.2,0.3, 0.5,0.7,0.9, 1.,1.15,1.25,1.4, 1.5,1.6,1.75, 2.,2.1,2.25,2.4,2.5,2.6,2.75,2.9,3.,3.1,3.25,3.4,3.5,3.6,3.75,3.9,4.,4.25,4.5,4.75, 5.};//,6.,7.,8.5, 10.};
		int size = sizeof(temp_vals)/sizeof(double);
//		double *temp_vals=temp;
		gridMarkers = convert_arr_to_vec(temp_vals, size);
	}
	return gridMarkers;
}

void make_rpdfWS_withCat_1D(TString outdir=".", TString tag="4e",int quad=9){
  gROOT->ProcessLine("gSystem->Load(\"make_rpdfWS_withCat_1D_c.so\")");
  gROOT->ProcessLine("gSystem->AddIncludePath(\"-I$ROOFITSYS/include/\")");
  gROOT->ProcessLine("gSystem->Load(\"libRooFit\")");
  gROOT->ProcessLine("gSystem->Load(\"libHiggsAnalysisCombinedLimit.so\")");

  	int nprods=2;
	TString prodtags[]={"ggH","qqH"};

	for(int cat =0; cat<1; cat++){
 		for(int prod_cate=0; prod_cate<nprods; prod_cate++){
 			TString outpath = outdir+"/1DSpline/"+tag+"/Dmass_Category"+to_string(cat)+"/"+prodtags[prod_cate];
			cout<<outpath<<endl;
 			system(Form("mkdir -p %s",outpath.Data()));
	 		do_make_rpdfWS_withCat(outpath,tag,cat,quad,prod_cate);
		
    }	
  }
}
