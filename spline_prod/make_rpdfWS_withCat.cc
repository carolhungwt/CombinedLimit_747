#include <iostream>
#include <utility>
#include <fstream>
#include <cmath>
#include <string>
#include <algorithm>
#include "TSystem.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include "TGaxis.h"
#include "TString.h"
#include "TChain.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussModel.h"
#include "RooRealIntegral.h"
#include "RooDecay.h"
#include "RooBMixDecay.h"
#include "RooCategory.h"
#include "RooBinning.h"
#include "RooRealSumPdf.h"
#include "RooPlot.h"
#include "RooNumIntConfig.h"
#include "RooWorkspace.h"
#include <ZZMatrixElement/MELA/interface/Mela.h>
#include <ZZMatrixElement/MELA/interface/ScalarPdfFactory_HVV.h>
#include "HiggsAnalysis/CombinedLimit/interface/RooNCSplineFactory_1D.h"
#include "HiggsAnalysis/CombinedLimit/interface/Width_conv.h"
#include "HiggsAnalysis/CombinedLimit/interface/SplinePdf.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooNCSplineCore.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooNCSplineFactory_2D.h"
#include "/cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/gcc/4.9.1-cms/include/c++/4.9.1/bits/stl_pair.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooNCSplineFactory_3D.h"
#include <time.h>
#include <iomanip>
#include <ctime>
using namespace RooFit;
using std::vector;
using std::endl;
using std::cout;

RooNCSpline_3D_fast* make_3D_spPDF(std::vector<RooRealVar*> projvars, vector<vector<double>> projvals, vector<int> ndim, RooWorkspace* w,  TString tag,int r,int cate_vbf,int prod_cate);
RooNCSpline_2D_fast* make_2D_spPDF(std::vector<RooRealVar*> projvars, vector<vector<double>> projvals, vector<int> ndim, RooWorkspace* w,  TString tag ,int r,int cate_vbf,int prod_cate);
RooNCSpline_1D_fast* make_1D_spPDF(RooRealVar projvar, vector<double> projval, int ndim, RooWorkspace* w, TString tag, int cate_vbf, int prod_cate);
void setVars(RooWorkspace* w, vector<float> Vals);
void setVar(RooWorkspace* w, const int var, const float Val);
vector<double> set_grid(RooRealVar* var);
vector<int> gen_ndim(vector<vector<double>> vec_vec);
vector<double> convert_arr_to_vec(double vals[],int size);
RooWorkspace* createWS(vector<RooAbsReal*> pdfs,TString rootname, vector<RooRealVar*> projvars);
RooPlot* quickplot(RooWorkspace* wo, vector<RooAbsPdf*> spPDFs, vector<RooRealVar*> projvars, const unsigned int chosen, vector<float> varset, int rval);
void plotOnC(TCanvas* c, RooPlot* plot, const char *tag,int rval);

void do_make_rpdfWS_withCat(TString workdir, TString tag, int cat,int quad, int cate_vbf, int prod_cate){
	
	// For timing the script
	clock_t tStart = clock();

	char* inputDir = "../signal_preparation";
	char* tag2 = (char*) tag.Data();
	TString tag1="";
	if(cate_vbf==0)		tag1="0jet_";
	else if(cate_vbf==1)	tag1="1jet_";
	else if(cate_vbf==2)	tag1="2jet_";
	TString prodtag="";
	if(prod_cate==0)	prodtag="";
	else if(prod_cate==1)	prodtag="vbf_";
	else if(prod_cate==2)	prodtag="vh_";
	TFile *fm = new TFile(Form("%s/workspace125_onshell/hzz%s_13TeV.input_func_%s%scat%d.root",inputDir,tag2,prodtag.Data(),tag1.Data(),cat));
	RooWorkspace* ws = (RooWorkspace*) fm->Get("w");


	// Get all variables from workspace
	std::vector<RooRealVar*> projvars;
	RooRealVar *m12 = ws->var("mean_pole");
	m12->setRange(122.,128.);
	RooRealVar *m1 = ws->var("mreco");
	m1->setRange(105.,140.);
	RooRealVar *m2 = ws->var("sigma_pole");
	m2->setRange(0.,5.);
	RooRealVar *r = ws->var("r");
	r->setVal(1);

	// 0: mH, 1: m4l, 2: gammaH
	projvars.push_back(m12);
	projvars.push_back(m1);
	projvars.push_back(m2);


// Construct 3D spline 
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
	//Make 2D spline for signal
	RooNCSplineCore* pdfr1 = make_3D_spPDF(projvars, projvals, ndim, ws,"r1",1,cate_vbf,prod_cate);


	vector<RooAbsReal*> vec_pdfs;
	TString curdir = ".";
	vec_pdfs.push_back(pdfr1);
	if(prod_cate==0)	prodtag="ggH_";
	RooWorkspace* rpdfWS = createWS(vec_pdfs,workdir+"/rpdfWS_withCat/"+tag+"_rpdfWS_cat"+to_string(cat)+"_"+prodtag+tag1+".root", projvars);

	rpdfWS->Print();


	return;
	// Timing the script
	string Timefile = Form("%s/Timefile.txt",workdir.Data());
	auto t = std::time(nullptr);
	auto tm = *std::localtime(&t);
}


RooNCSpline_3D_fast* make_3D_spPDF(std::vector<RooRealVar*> projvars, vector<vector<double>> projvals, vector<int> ndim, RooWorkspace* w,  TString tag ,int r,int cate_vbf,int prod_cate){
	double* dim = new double[projvars.size()];
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
	
	RooNCSplineFactory_3D* spFac =  new RooNCSplineFactory_3D(*(projvars.at(0)), *(projvars.at(1)), *(projvars.at(2)));
	double fcn = 0.;
	ofstream out("crosscheck"+tag+".txt");
	vector<splineQuadruplet_t> points;
		for(int ix=0; ix<ndim[0]; ix++){
			dim[0]=projvals.at(0).at(ix);
			w->var("mean_pole")->setConstant(false);
                        w->var("mean_pole")->setVal(dim[0]);
                        w->var("mean_pole")->setConstant(true);
			for (int iy=0; iy<ndim[1]; iy++){
       			dim[1]=projvals.at(1).at(iy);
       			w->var("mreco")->setConstant(false);
        		w->var("mreco")->setVal(dim[1]);
        		w->var("mreco")->setConstant(true);
				for (int iz=0; iz<ndim[2]; iz++){
    			    dim[2]=projvals.at(2).at(iz);
		        	w->var("sigma_pole")->setConstant(false);
          			w->var("sigma_pole")->setVal(dim[2]);
          			w->var("sigma_pole")->setConstant(true);	

    			    double fcn = pdf->getVal();
    			    cout << dim[0]<<" "<<dim[1]<<" "<<dim[2]<< " "<<fcn<<" \n";
					points.push_back(splineQuadruplet_t(dim[0], dim[1], dim[2], fcn));
				}
			}
		}
	spFac -> setPoints(points);
	RooNCSpline_3D_fast* spPDF = spFac->getFunc();
	spPDF->SetNameTitle(tag,tag);
	return spPDF;
}


void plotOnC(TCanvas* c, RooPlot* plot, const char *tag,int rval){
	c->SetGrid(0,1);
	int cat=0;
	c->cd();
	plot->Draw();
	c->SaveAs(Form("~/www/spline_test/r%d_test/test_spline_signal_%s_cat%d.png",rval,tag,cat));
}

RooPlot* quickplot(RooWorkspace* wo, vector<RooAbsPdf*> spPDFs, vector<RooRealVar*> projvars, const unsigned int chosen, vector<float> varset, int rval){
	RooPlot* plot = projvars[chosen]->frame(24);
	plot->GetXaxis()->CenterTitle();
    plot->GetYaxis()->SetTitleOffset(1.2);
    plot->GetYaxis()->CenterTitle();
    plot->GetXaxis()->SetTitle("M_{4l} (GeV)");//projvars.at(1)->GetTitle());
    plot->GetYaxis()->SetTitle("H4l amplitude");
    if(rval==0) plot->GetYaxis()->SetRangeUser(0,0.1);
    if(rval==1) plot->GetYaxis()->SetRangeUser(0,0.3);
    if(rval==2) plot->GetYaxis()->SetRangeUser(0,0.3);
    plot->GetXaxis()->SetNdivisions(-505);
    plot->SetTitle(Form("mH = %0.2f GeV #Gamma_{H} = %0.2f GeV",varset[0],varset[2]));
    for(unsigned int i=0; i<projvars.size(); i++){
    	if (i==chosen)	continue;
        setVar(wo,i,varset[i]);
	projvars[i]->setConstant(false);
	projvars[i]->setVal(varset[i]);
	projvars[i]->setConstant(true);
    }
    setVar(wo,3,rval);
    wo->pdf("ggH")->plotOn(plot,LineColor(1), LineWidth(2), LineStyle(1));
    plot->SetTitle(Form("mH = %0.2f GeV #Gamma_{H} = %0.2f GeV",varset[0],varset[2]));
    spPDFs[rval]->plotOn(plot,LineColor(2), LineWidth(2), LineStyle(2));

//	varset[0]+=1.;
//	setVar(wo,0,varset[0]);
//     spPDFs[rval]->plotOn(plot,LineColor(1), LineWidth(2), LineStyle(1));
	return plot;
}


RooWorkspace* createWS(vector<RooAbsReal*> pdfs,TString rootname,  vector<RooRealVar*> projvars){
	TFile ftest(rootname.Data(),"recreate");

	ftest.cd();
	RooWorkspace* w = new RooWorkspace("w", "");
	w->addClassDeclImportDir("../interface/");
	w->importClassCode(RooNCSplineCore::Class(),true);
	w->importClassCode(RooNCSpline_1D_fast::Class(), true);
	w->importClassCode(RooNCSpline_2D_fast::Class(), true);		
	w->importClassCode(RooNCSpline_3D_fast::Class(), true);
	int j=0;
	for(unsigned int i=0; i<pdfs.size(); i++){	
		cout<<pdfs[j]->GetName()<<endl;
		cout<<pdfs[j]<<endl;
		w->import(*pdfs[j], RooFit::RecycleConflictNodes());
		j++;
		cout<<j<<endl;
	}
	for(unsigned int i=0; i<projvars.size();i++){
		w->import(*projvars[i],RooFit::RecycleConflictNodes());
	}	

//	return 0;
	ftest.WriteTObject(w);
	ftest.Close();
	return w;

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
		inc = 0.05;
		double min_mH = 122.;
		double max_mH = 128.;
		gridnum = (max_mH-min_mH)/inc+1;
		for(int i = 0; i<gridnum; i++){
			gridMarkers.push_back(min_mH+inc*i);
		}	
	}	else if(varName == "mreco"){
		double temp_vals[]={105.,106., 107.,108., 109., 110.,111.,112., 113., 114., 115.,116., 117.,118., 119., 120.,121., 122., 123., 124., 125., 126., 127., 128., 129.5, 131., 132.5, 134., 135.5, 137., 138.5, 140.};
//		inc = 0.5;
//		double min_m4l = 105.;
//		double max_m4l = 140.;
//		gridnum = (max_m4l-min_m4l)/inc+1;
//		double temp_vals[gridnum];
//		for(int i=0; i<gridnum; i++){
//			gridMarkers.push_back(min_m4l+inc*i);
//		}
		int size = sizeof(temp_vals)/sizeof(double);
		gridMarkers = convert_arr_to_vec(temp_vals, size);
	}	else if(varName == "sigma_pole"){
		double temp_vals[]={0.0001,0.001, 0.002,.00407,0.007,0.01,0.05,0.08,0.1,0.15, 0.2,0.3, 0.5,0.7,0.9, 1.,1.25,1.4, 1.5,1.6,1.75,1.9, 2.,2.1,2.25,2.4,2.5,2.6,2.75,2.9,3.,3.1,3.25,3.4,3.5,3.6,3.75,3.9,4.,4.25,4.5,4.75, 5.};//,6.,7.,8.5, 10.};
		int size = sizeof(temp_vals)/sizeof(double);
//		double *temp_vals=temp;
		gridMarkers = convert_arr_to_vec(temp_vals, size);
	}
	return gridMarkers;
}


void setVars(RooWorkspace* w, const vector<float> Vals){
	TString curvar="";
	for(unsigned int i=0; i<Vals.size(); i++){
		if(i==0)	curvar="mean_pole";
		if(i==1)	curvar="mreco";
		if(i==2)	curvar="sigma_pole";
		w->var(curvar)->setConstant(false);
		w->var(curvar)->setVal(Vals[i]);
		w->var(curvar)->setConstant(true);
	}
}

void setVar(RooWorkspace* w, const int var, const float Val){
        TString curvar = "";
        if(var==0)  curvar="mean_pole";
        if(var==1)  curvar="mreco";
        if(var==2)  curvar="sigma_pole";
        if(var==3)  curvar="r";
        w->var(curvar)->setConstant(false);
        w->var(curvar)->setVal(Val);
        w->var(curvar)->setConstant(true);
}
void make_rpdfWS_withCat(TString workdir=".", TString tag="4e", int cat=0,int quad=9){
  gROOT->ProcessLine("gSystem->AddIncludePath(\"-I$ROOFITSYS/include/\")");
  gROOT->ProcessLine("gSystem->Load(\"libRooFit\")");
  gROOT->ProcessLine("gSystem->Load(\"libHiggsAnalysisCombinedLimit.so\")");
for(int cate_vbf=0; cate_vbf<3; cate_vbf++){
   for(int prod_cate=0; prod_cate<3; prod_cate++){
         do_make_rpdfWS_withCat(workdir,tag, cat,quad, cate_vbf, prod_cate);
        }
    }
}
