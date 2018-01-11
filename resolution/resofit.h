#include "TDirectory.h"
#include "TPad.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "THStack.h"
#include "TH2.h"
#include "TSpline.h"
#include "TF1.h"
#include "TLine.h"
#include "TCut.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"
//#include "HZZ2L2QRooPdfs.h"
//#include "HZZ2L2QRooPdfs.cc"
#include "HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooDoubleCBFast.h"
#include "HiggsAnalysis/CombinedLimit/interface/Category.h"
#include "HiggsAnalysis/CombinedLimit/interface/Discriminants.h"
//#include "ZZAnalysis/AnalysisStep/interface/Category.h"
//#include "Math/GenVector/LorentzVector.h"
//#include "Math/GenVector/PtEtaPhiM4D.h"
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <string>

//#include "Math/GenVector/LorentzVector.h"
//#include "Math/GenVector/PtEtaPhiM4D.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooWorkspace.h"
#include "RooHist.h"
#include "TCanvas.h"

#include "TLegend.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "TPaveText.h"
#include "RooAddPdf.h"
#include "RooBreitWigner.h"
#include "RooFitResult.h"
#include "RooFFTConvPdf.h"
#include "RooAddition.h"
#include "RooMinuit.h"
#include "Math/MinimizerOptions.h"

#include "RooAbsCollection.h"


using namespace RooFit;
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>


  float m4l,genM,Z1Mass,Z2Mass;
  std::vector<float> *wt=0;
  Short_t z1flav, z2flav,zzsel;
  float weight;
  char tempmass[100];
  char tempmass2[100];
  double width[100];
  double xMinF=-2000;
  double xMaxF=1500;
  double xMin[100];
  double xMax[100];
  int total_ent=0;
  int fe_ent=0;
  int fmu_ent=0;
  int tetmu_ent=0;


  RooRealVar x("reso","m_{reco}-m_{true}",0.,xMinF,xMaxF,"GeV");
  RooRealVar w("myW","myW",1.0,-2000.,1500.);
  RooCategory massrc("massrc","massrc");
  RooDataSet* dataset;
  RooDataSet* dataset_sub[100];
  vector<TString> files;
  char inputfile[1000];

  RooDoubleCBFast* DCBall[100];
  RooConstVar* xmgen[100];
  RooRealVar* mean_ind[100];
  RooRealVar* sigma_ind[100];
  RooRealVar* a1_ind[100];
  RooRealVar* a2_ind[100];
  RooRealVar* n1_ind[100];
  RooRealVar* n2_ind[100];
  RooFitResult* fitres[100];

//  string inputDir =;
//  TString inputDir = "/afs/cern.ch/user/w/wahung/work/public/CombineLimitDbkgkin/resolution_pkg/resolution_Quad9/Quad9_tree";

//  TString TreeinputDir="/eos/user/w/wahung/Mass_Width_Measurement/Data_Sample/";

  const int maxCat = 9;
//  int inputfiles[]={125};
//  int massBin[]={125};
  int massBin[]={115,120,124,125,126,130,135,140};//,145,150,155,160,165,170,175,180,190,200};
  int inputfiles[]={115,120,124,125,126,130,135,140};//,145,150,155,160,165,170,175,180,190,200};
  int maxMassBin=sizeof(massBin)/sizeof(*massBin);;
  int Nfiles=sizeof(inputfiles)/sizeof(*inputfiles);

