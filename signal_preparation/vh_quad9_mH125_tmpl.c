/** Things to remember later
 * mzz range defined by ggZZ workspace, create 2 versions for onshell and offshell
 * give different ggzz RooKeysPdf name for 2e2mu and 4e
 *
 ****/
using namespace std;
using namespace RooFit;

void dosomething(TString chan ="2e2mu",unsigned int cate_vbf=1,bool onshell=true, int quad=3, int iCat=0){
	int nparampar=2;
	double lumi = 35.9;
	double ggzz_xsec = 0.488;
	double x_xsec = 0.312*1.108; 

	if(onshell){
		if(chan=="2e2mu"){
		ggzz_xsec = 0.1009*0;
		x_xsec = 0.275*1.108*0.122929291; 
		}
		else if(chan=="4mu"){
			ggzz_xsec = 0.0473*0;
			x_xsec = 0.149*1.108*0.121647465;
		}
		else if(chan=="4e"){
			ggzz_xsec = 0.0473*0;
			x_xsec = 0.149*1.108*0.122419844;
		}
	}
		else{
			if(chan!="2e2mu"){
			ggzz_xsec=0.239;
			x_xsec = 0.167*1.108;; 
			}
		}

	gStyle->SetOptStat(0);
<dcbPara_4e_2nd>

<dcbPara_4mu_2nd>

<dcbPara_2e2mu_2nd>

	double eff[3][11]={
//VBF->0p_4e_param_0:
-4.457080E+00, 4.583149E+00, -3.220522E+02, 2.772364E+02, 2.590067E+00, 1.438727E-03, -5.122000E-07, 1.369010E-01, 1.296028E+02, 6.672205E+01, 4.702222E-11,
//VBF->0p_4mu_param_0:
-4.464713E+00, 4.575797E+00, -7.792695E+03, 4.494984E+03, 1.001845E+01, -1.132429E-03, -6.647936E-08, -4.162270E-01, 5.929335E+01, 4.422993E+01, 4.317671E-11,
//VBF->0p_2e2mu_param_0:
-4.487157E+00, 4.555744E+00, -3.852006E+02, 3.212818E+02, 6.250167E+00, 3.400357E-03, -1.287049E-06, 6.712138E-01, 6.316531E+01, 1.033601E+02, 1.402357E-10,
};

 double eff_adj[3][9][2]={
//Channel 4e
{{0.142951, -0.000237},{0.140153, -0.000233},{0.066292, 0.000366},{0.034540, 0.000634},{0.086560, 0.000192},{0.109414, -0.000006},{0.096713, 0.000119},{0.152212, -0.000342},{0.171165, -0.000493}},
//Channel 4mu
{{0.031283, 0.000614},{0.086497, 0.000176},{0.102232, 0.000075},{0.047921, 0.000503},{0.035003, 0.000641},{0.055539, 0.000465},{0.136926, -0.000230},{0.256384, -0.001143},{0.248215, -0.001100}},
////Channel 2e2mu
{{0.033567, 0.000617},{0.066665, 0.000353},{0.085146, 0.000204},{0.074805, 0.000285},{0.020880, 0.000707},{0.059935, 0.000412},{0.107745, 0.000028},{0.226598, -0.000925},{0.324660, -0.001680}},
};

  double dcbPara_2nd[6][3];
  double* effsig;
  double* effadj;

  if (chan=="4e"){
    for (int p=0; p<6; p++){ for (int k=0; k<nparampar; k++) dcbPara_2nd[p][k]=dcbPara_4e_2nd[p][iCat][k]; }
    effsig=eff[0];
    effadj = eff_adj[0][iCat];
  }
  if (chan=="4mu"){
    for (int p=0; p<6; p++){ for (int k=0; k<nparampar; k++) dcbPara_2nd[p][k]=dcbPara_4mu_2nd[p][iCat][k]; }
    effsig=eff[1];
    effadj = eff_adj[1][iCat];
  }
  if (chan=="2e2mu"){
    for (int p=0; p<6; p++){ for (int k=0; k<nparampar; k++) dcbPara_2nd[p][k]=dcbPara_2e2mu_2nd[p][iCat][k]; }
    effsig=eff[2];
    effadj = eff_adj[2][iCat];
  }



	double lowarr[2]={100.5,115};
	double higharr[2]={1500.5,140};
	const int nbinsarr[2]={1500,100};

	double recolowarr[2]={104,105};
	double recohigharr[2]={1604.,140.};
	const int reconbinsarr[2]={750,100};

	const double low= lowarr[onshell];
	const double high=higharr[onshell];
	const int nbins= nbinsarr[onshell]; 

	const double low_reco=recolowarr[onshell];
	const double high_reco=recohigharr[onshell];
	const int nbins_reco=reconbinsarr[onshell];

	cout << low<<"\t"<<high<<endl;
	cout << low_reco<<"\t"<<high_reco<<endl;

	TString ap = (onshell? "_onshell":"");
	TFile *fpdfbkg = new TFile("vbfpdfs"+ap+".root");
	RooWorkspace *wbkg =( RooWorkspace*)fpdfbkg->Get ("w"); 

	RooRealVar* mzz = new RooRealVar("ZZMass","M_{ZZ} (GeV)",125,low_reco,high_reco);
	//RooRealVar* mzz = wbkg->var("ZZMass"); 
	RooRealVar* mreco= new RooRealVar("mreco","M_{ZZ} (GeV)",125,low_reco,high_reco);
	RooRealVar* mdiff= new RooRealVar("mdiff","M_{ZZ} (GeV)",125,low_reco,high_reco);

	RooRealVar *r= new RooRealVar("r","signal strength",1.,0.000,1000);
	RooRealVar *rvbf_ggh= new RooRealVar("rvbf_ggh","rvb_ggh",1.,0.000,1000);
	RooFormulaVar *rvbf= new RooFormulaVar("rvbf","@0*@1",RooArgSet(*r,*rvbf_ggh));

	RooRealVar* mean = new RooRealVar("mean_pole","mean_pole",125,100,180);
	RooRealVar* sigma= new RooRealVar("sigma_pole","sigma_pole",0.00418,0.,10.);

	RooConstVar* mean_125= new RooConstVar("mean_125","mean_125",125);
	RooConstVar* sigma_125= new RooConstVar("sigma_125","sigma_125",0.00407);

	RooPlot* frame= mreco->frame(low_reco,high_reco) ;
	//	RooPlot* frame= mreco->frame(150,250) ;
	RooPlot* frame_mzz= mzz->frame(Title("Z mass")) ;
	RooPlot* frame_width= sigma->frame(Title("width")) ;
	RooPlot* frame_mean= mean->frame(Title("mean")) ;

	TFile *flo=new TFile("width_new_spline.root","read");
	TString chn = "2e2mu";
	if(chan!="2e2mu")
		chn="4e";
	TSpline3 *lo=(TSpline3*) flo->Get("br_"+chn);

	SplinePdf par2_int ("vbfpar2_int"+chan+Form("%d_%djet",iCat,cate_vbf),"",*mzz,*mean,*sigma,*lo);
	RooRealVar m_gauss ("m_gauss","",125);
	RooRealVar w_gauss ("w_gauss","",0.004);
	RooGaussian gauss("gauss","",*mzz,m_gauss,w_gauss);

	TString pdfn = "2e2mu";
	if(chan!="2e2mu")
		pdfn = "4e";
	RooAbsPdf *pdfbkg = (RooAbsPdf*) wbkg->pdf("vbfpdfbkg_"+pdfn);

	float zero_bkg = 0;
	RooConstVar *ggzznorm= new RooConstVar("vhbkgnorm"+chan+Form("%d_%djet",iCat,cate_vbf),"",lumi*ggzz_xsec*zero_bkg);
	RooExtendPdf pdf_ggzz("vhpdf_bkg"+chan+Form("%d_%djet",iCat,cate_vbf),"vhpdf_bkg"+chan+Form("%d_%djet",iCat,cate_vbf),*pdfbkg,*ggzznorm);
//	pdf_ggzz.plotOn(frame_mzz);
//	frame_mzz->Draw();
//	return;

	RooConstVar *xnorm= new RooConstVar("vhxnorm"+chan+Form("%d_%djet",iCat,cate_vbf),"",lumi*x_xsec);
	RooExtendPdf pdf_x("vhpdf_x"+chan+Form("%d_%djet",iCat,cate_vbf),"vhpdf_x"+chan+Form("%d_%djet",iCat,cate_vbf),par2_int,*xnorm);


  TFile *fphase_noweight=new TFile("fgraph_vbf_phase.root");
  TGraph *cosfunc = (TGraph*)fphase_noweight->Get("cos");
  TGraph *sinfunc = (TGraph*)fphase_noweight->Get("sin");


	TFile *fbkg = new TFile("./bkg_reg_eff.root");
	TGraph *eff_bkg =  (TGraph*)fbkg->Get("VBF_reg_"+chan);


	TGraph *effxkf_sig= new TGraph(nbins);
	TGraph *effxkf_bkg= new TGraph(nbins);

  TFile *fxsec=new TFile("xsec.root");
  TGraph *vbfxs= (TGraph*)fxsec->Get("vbf");
  TGraph *whxs= (TGraph*)fxsec->Get("wh");
  TGraph *zhxs= (TGraph*)fxsec->Get("zh");


	for(int i =0;i<nbins;i++){
		double cva = low+ i*(high-low)/double(nbins);
		double effval_sig = (effsig[0]+effsig[1]*TMath::Erf( (cva-effsig[2])/effsig[3] ))*(effsig[4]+effsig[5]*cva+effsig[6]*cva*cva+effsig[10]*cva*cva*cva)+effsig[7]*TMath::Gaus(cva,effsig[8],effsig[9]);
    double m4l = cva;
    //if (m4l > 900.) m4l = 900.;
    //double effcate = 3.844411e-01 + 2.755289e-04*m4l + -5.098494e-07*m4l*m4l + 2.491489e-10*m4l*m4l*m4l;
	//	double effcate_sig = 0.476383 - 4.94232e-05*m4l;
		double effcate_sig = 0.6325;
		double effcate_bkg = 0;
		double vhlep_vbf2j_sig = 0.0984;

		double vh_vhlep_sig = 1-effcate_sig-vhlep_vbf2j_sig;
//		double vhlepfrac=1;
//		if(cva<2000)
		double vhlepfrac = 0.3*(1-vbfxs->Eval(cva)*0.955/(whxs->Eval(cva)*0.654+zhxs->Eval(cva)*0.669+vbfxs->Eval(cva)*0.955))/0.7;
		if(cate_vbf==0){	// untagged
				effcate_sig= effcate_sig*vhlepfrac;
		}
   	else if(cate_vbf==1){	// vbf2j tagged 
       effcate_sig= vhlepfrac* vhlep_vbf2j_sig;
		}
	else if (cate_vbf==2){	// vhlep tagged
		effcate_sig = vhlepfrac*vh_vhlep_sig;

	}

		double va_bkg= eff_bkg->Eval(cva)*effcate_bkg; 	
		double va_sig= effval_sig*effcate_sig; 	
//		if(chan=="2e2mu"){
//			va_bkg*=0.95;
//			va_sig*=0.95;
//		}
		effxkf_sig->SetPoint(effxkf_sig->GetN(),cva,va_sig);
		effxkf_bkg->SetPoint(effxkf_bkg->GetN(),cva,va_bkg);
	}
	effxkf_sig->SetName("vhsigeffxkf"+chan+Form("%d_%djet",iCat,cate_vbf));
	effxkf_bkg->SetName("vhbkgeffxkf"+chan+Form("%d_%djet",iCat,cate_vbf));

	//effxkf_sig->Draw("al");
	//effxkf_bkg->Draw("lsame");
	//return;

	mean->setRange(100,1500);
	sigma->setRange(0.00005,100.);
	mean->setVal(125);
	sigma->setVal(0.004165);

	TString formu_2nd = "@1+@0*@2";	

	RooArgList formuList_a1;
	RooArgList formuList_a2;
	RooArgList formuList_mean;
	RooArgList formuList_n1;
	RooArgList formuList_n2;
	RooArgList formuList_sigma;
	formuList_a1.add(*mzz);
	formuList_a2.add(*mzz);
	formuList_mean.add(*mzz);
	formuList_n1.add(*mzz);
	formuList_n2.add(*mzz);
	formuList_sigma.add(*mzz);

	RooConstVar* a1_p0_0_2nd[11] ;
	RooConstVar* a2_p0_0_2nd[11] ;
	RooConstVar* mean_p0_0_2nd[11] ;
	RooConstVar* n1_p0_0_2nd[11] ;
	RooConstVar* n2_p0_0_2nd[11] ;
	RooConstVar* sigma_p0_0_2nd[11] ;
  for (int i=0; i<nparampar; i++){
    a1_p0_0_2nd[i]= new RooConstVar(Form("%s_%d_a1_p0_0_2nd_iCat%d_%djet", chan.Data(), i, iCat, cate_vbf), Form("%s_%d_a1_p0_0_2nd_iCat%d_%djet", chan.Data(), i, iCat,cate_vbf), dcbPara_2nd[0][i]);
    a2_p0_0_2nd[i]= new RooConstVar(Form("%s_%d_a2_p0_0_2nd_iCat%d_%djet", chan.Data(), i, iCat,cate_vbf), Form("%s_%d_a2_p0_0_2nd_iCat%d_%djet", chan.Data(), i, iCat, cate_vbf), dcbPara_2nd[1][i]);
    mean_p0_0_2nd[i]= new RooConstVar(Form("%s_%d_mean_p0_0_2nd_iCat%d_%djet", chan.Data(), i, iCat, cate_vbf), Form("%s_%d_mean_p0_0_2nd_iCat%d_%djet", chan.Data(), i, iCat, cate_vbf), dcbPara_2nd[2][i]);
    n1_p0_0_2nd[i]= new RooConstVar(Form("%s_%d_n1_p0_0_2nd_iCat%d_%djet", chan.Data(), i, iCat, cate_vbf), Form("%s_%d_n1_p0_0_2nd_iCat%d_%djet", chan.Data(), i, iCat, cate_vbf), dcbPara_2nd[3][i]);
    n2_p0_0_2nd[i]= new RooConstVar(Form("%s_%d_n2_p0_0_2nd_iCat%d_%djet", chan.Data(), i, iCat, cate_vbf), Form("%s_%d_n2_p0_0_2nd_iCat%d_%djet", chan.Data(), i, iCat, cate_vbf), dcbPara_2nd[4][i]);
    sigma_p0_0_2nd[i]= new RooConstVar(Form("%s_%d_sigma_p0_0_2nd_iCat%d_%djet", chan.Data(), i, iCat, cate_vbf), Form("%s_%d_sigma_p0_0_2nd_iCat%d_%djet", chan.Data(), i, iCat, cate_vbf), dcbPara_2nd[5][i]);
std::cout<<"***********************************"<<endl;
if(i==0) std::cout<<(a1_p0_0_2nd[i])->getValV()<<endl;
std::cout<<"***********************************"<<endl;
    formuList_a1.add(*a1_p0_0_2nd[i]);
    formuList_a2.add(*a2_p0_0_2nd[i]);
    formuList_mean.add(*mean_p0_0_2nd[i]);
    formuList_n1.add(*n1_p0_0_2nd[i]);
    formuList_n2.add(*n2_p0_0_2nd[i]);
    formuList_sigma.add(*sigma_p0_0_2nd[i]);
  }
	RooWorkspace w("w");
  RooFormulaVar* a1_p0_2nd= new RooFormulaVar("a1_p0_2nd"+chan+Form("_iCat%d_%djet",iCat, cate_vbf), "a1_p0_2nd"+chan+Form("_iCat%d_%djet",iCat, cate_vbf), formu_2nd, formuList_a1);
  RooFormulaVar* a2_p0_2nd= new RooFormulaVar("a2_p0_2nd"+chan+Form("_iCat%d_%djet",iCat, cate_vbf), "a2_p0_2nd"+chan+Form("_iCat%d_%djet",iCat, cate_vbf), formu_2nd, formuList_a2);
  RooFormulaVar* mean_p0_2nd= new RooFormulaVar("mean_p0_2nd"+chan+Form("_iCat%d_%djet",iCat, cate_vbf), "mean_p0_2nd"+chan+Form("_iCat%d_%djet",iCat, cate_vbf), formu_2nd, formuList_mean);
  RooFormulaVar* n1_p0_2nd= new RooFormulaVar("n1_p0_2nd"+chan+Form("_iCat%d_%djet",iCat, cate_vbf), "n1_p0_2nd"+chan+Form("_iCat%d_%djet",iCat, cate_vbf), formu_2nd, formuList_n1);
  RooFormulaVar* n2_p0_2nd= new RooFormulaVar("n2_p0_2nd"+chan+Form("_iCat%d_%djet",iCat, cate_vbf), "n2_p0_2nd"+chan+Form("_iCat%d_%djet",iCat, cate_vbf), formu_2nd, formuList_n2);
  RooFormulaVar* sigma_p0_2nd= new RooFormulaVar("sigma_p0_2nd"+chan+Form("_iCat%d_%djet",iCat, cate_vbf), "sigma_p0_2nd"+chan+Form("_iCat%d_%djet",iCat, cate_vbf), formu_2nd, formuList_sigma);


  RooFormulaVar* sigma_p0_up = new RooFormulaVar("sigma_p0_up", "", "@0+0.2*@0", *sigma_p0_2nd);
  RooFormulaVar* sigma_p0_dn = new RooFormulaVar("sigma_p0_dn", "", "@0-0.2*@0", *sigma_p0_2nd);
  RooDoubleCB dcrReso(Form("dcrReso_%s_iCat%d_%djet",chan.Data(), iCat, cate_vbf), "Double Crystal ball ", *mreco, *mean, *mean_p0_2nd, *sigma_p0_2nd, *a1_p0_2nd, *n1_p0_2nd, *a2_p0_2nd, *n2_p0_2nd);
  RooDoubleCB dcrReso_up(Form("dcrReso_%s_iCat%d_%djet_up",chan.Data(), iCat, cate_vbf), "dcb up", *mreco, *mean, *mean_p0_2nd, *sigma_p0_up, *a1_p0_2nd, *n1_p0_2nd, *a2_p0_2nd, *n2_p0_2nd);
  RooDoubleCB dcrReso_dn(Form("dcrReso_%s_iCat%d_%djet_dn",chan.Data(),iCat, cate_vbf), "dcb dn", *mreco, *mean, *mean_p0_2nd, *sigma_p0_dn, *a1_p0_2nd, *n1_p0_2nd, *a2_p0_2nd, *n2_p0_2nd);

	RooAbsReal *final_integral; 

//if(onshell){
	Width_conv convpdf_spline(Form("VH_%s_iCat%d_%djet",chan.Data(),iCat, cate_vbf), Form("VH_%s_iCat%d_%djet",chan.Data(),iCat, cate_vbf),*mreco, *mean, *sigma, *rvbf, RooArgList(pdf_x, pdf_ggzz,dcrReso),*cosfunc, *sinfunc, *effxkf_sig,*effxkf_bkg); 
	convpdf_spline.SetNameTitle("VH","VH");
//	convpdf_spline.plotOn(frame);
	w.import(convpdf_spline,RecycleConflictNodes());
	Width_conv convpdf_spline_up("VH_Res"+chan+"Up", "VH_Res"+chan+"Up",*mreco, *mean, *sigma, *r, RooArgList(pdf_x, pdf_ggzz,dcrReso_up),*cosfunc, *sinfunc, *effxkf_sig,*effxkf_bkg); 
	Width_conv convpdf_spline_dn("VH_Res"+chan+"Down", "VH_Res"+chan+"Down",*mreco, *mean, *sigma, *r, RooArgList(pdf_x, pdf_ggzz,dcrReso_dn),*cosfunc, *sinfunc, *effxkf_sig,*effxkf_bkg); 
	w.import(convpdf_spline_up,RecycleConflictNodes());
	w.import(convpdf_spline_dn,RecycleConflictNodes());
	final_integral = convpdf_spline.createIntegral(*mreco);
//}

	mean->setVal(125);
	mreco->setVal(126);

	r->setVal(0);
//	double bexp =  final_integral->getVal();
	double bexp = 0.;
	cout << bexp<<endl;

	RooConstVar *bkg_integral= new RooConstVar("bkg_integral"+chan+Form("%d_%djet",iCat, cate_vbf),"",bexp);

	mean->setVal(125);
	//ROOT::Math::Interpolator inter(200, ROOT::Math::Interpolation::kCSPLINE);
	sigma->setRange(0.00005,100.);

	TH2F *hint ; 
	TH2F *hsig ;
	if(!onshell){	
		hint= new TH2F("hint","",101,119.95,130.05,101,-0.0005,0.1005);
		hsig= new TH2F("hsig","",101,119.95,130.05,101,-0.0005,0.1005);
	}
	else{
		hint= new TH2F("hint","",101,119.95,130.05,101,-0.025,5.025);
		hsig= new TH2F("hsig","",101,119.95,130.05,101,-0.025,5.025);
		//hint= new TH2F("hint","",101,122.475,127.525,101,-0.025,5.025);
		//hsig= new TH2F("hsig","",101,122.475,127.525,101,-0.025,5.025);
	}
	//		double xi[201]; 
	//		double yi[201] ;

	for(int i = 0;i <101;i++){
		if(i%10==0)
			cout<<i<<endl;
		for(int j = 0;j <101;j++){
			double mv=hint->GetXaxis()->GetBinCenter(i+1);
			double sv=hint->GetYaxis()->GetBinCenter(j+1);
			sigma->setVal(sv);
			mean->setVal(mv);
			//sigma->setVal(0.001*(i+1));
			//mean->setVal(100+0.5*(j+1));
			r->setVal(1);
			double sbi =  final_integral->getVal();
			r->setVal(4);
			double sbi2 =  final_integral->getVal();
			double sexp = ((sbi2-sbi*2)+bexp)/2.;
			double iexp = sbi -sexp -bexp; 
			hint->Fill(mean->getVal(),sigma->getVal(),0);
			hsig->Fill(mean->getVal(),sigma->getVal(),sexp);
		}
	}
	RooDataHist* hinthist= new RooDataHist("vhhinthist"+chan+Form("%d_%djet",iCat,cate_vbf),"vhhinthist"+chan+Form("%d_%djet",iCat,cate_vbf),RooArgSet(*mean,*sigma),hint);
	RooHistFunc *hintfunc = new RooHistFunc("vhhintfunc"+chan+Form("%d_%djet",iCat,cate_vbf),"",RooArgSet(*mean,*sigma),*hinthist);
	Width_integral inter_intergral ("vhint_integral"+chan+Form("%d_%djet",iCat,cate_vbf),"",*mean,*sigma,RooArgList(*hintfunc));

	RooDataHist* hsighist= new RooDataHist("vhhsighist"+chan+Form("%d_%djet",iCat,cate_vbf),"hsighist"+chan+Form("%d_%djet",iCat,cate_vbf),RooArgSet(*mean,*sigma),hsig);
	RooHistFunc *hsigfunc = new RooHistFunc("vhhsigfunc"+chan+Form("%d_%djet",iCat,cate_vbf),"",RooArgSet(*mean,*sigma),*hsighist);
	Width_integral sig_intergral ("vhsig_integral"+chan+Form("%d_%djet",iCat,cate_vbf),"",*mean,*sigma,RooArgList(*hsigfunc));
	RooFormulaVar overall_intergral("VH_norm","","@0*@2+ @1 + sqrt(@2)*@3",RooArgList(sig_intergral, *bkg_integral, *rvbf,inter_intergral));

//****************** Calculate yield for different cats ******************
	fstream yield("vh_yield.txt", ios::out|ios::app);
	mean->setVal(125);
	sigma->setVal(0.004);
	yield<<"Channel: "<<chan<<"\t Category: "<<iCat<<"\t Jet Category: "<<cate_vbf<<endl;
	yield<<"Signal \t Background \t Interference \t"<<endl;
	yield<< sig_intergral.getVal()<<"\t"<< bkg_integral->getVal()<<"\t"<<inter_intergral.getVal()<<endl;

	cout<< sig_intergral.getVal()<<"\t"<< bkg_integral->getVal()<<"\t"<<inter_intergral.getVal()<<endl;

	Width_conv_integral convpdf_spline_integral("VH", "VH",*mreco, *mean, *sigma, *rvbf, RooArgList(pdf_x, pdf_ggzz,dcrReso,*hsigfunc,*hintfunc),*cosfunc, *sinfunc, *effxkf_sig,*effxkf_bkg); 
	Width_conv_integral convpdf_spline_integral_up("VH_Res"+chan+"Up", "VH_Res"+chan+"Up",*mreco, *mean, *sigma, *rvbf, RooArgList(pdf_x, pdf_ggzz,dcrReso_up,*hsigfunc,*hintfunc),*cosfunc, *sinfunc, *effxkf_sig,*effxkf_bkg); 
	Width_conv_integral convpdf_spline_integral_dn("VH_Res"+chan+"Down", "VH_Res"+chan+"Down",*mreco, *mean, *sigma, *rvbf, RooArgList(pdf_x, pdf_ggzz,dcrReso_dn,*hsigfunc,*hintfunc),*cosfunc, *sinfunc, *effxkf_sig,*effxkf_bkg); 

	w.import(convpdf_spline_integral,RecycleConflictNodes());
	w.import(convpdf_spline_integral_up,RecycleConflictNodes());
	w.import(convpdf_spline_integral_dn,RecycleConflictNodes());
	overall_intergral.plotOn(frame_width);
	frame_width->Draw();

	mzz->setConstant(0);
	mean->setConstant(0);
	sigma->setConstant(0);
	mreco->setConstant(0);

	mreco->setRange(low_reco,high_reco);
	mean->setVal(125);
	sigma->setVal(0.004);
	r->setVal(1);

	convpdf_spline.plotOn(frame,LineColor(2));
	frame->Draw();

	w.import(overall_intergral,RecycleConflictNodes());

	mreco->Print("v");
	cout << mreco->getBins();

	TString filename = "workspace125/hzz4l_"+chan+Form("_13TeV.input_func_vh_%djet_cat%d.root",cate_vbf, iCat);
	if(onshell)
		filename = "workspace125_onshell/hzz"+chan+Form("_13TeV.input_func_vh_%djet_cat%d.root",cate_vbf, iCat);
	TFile *f=new TFile(filename,"recreate");	
	f->cd();
	w.Write();
	f->Close();
	return;
	
}

void vh_quad9_mH125(TString chan="4e",unsigned int cate_vbf=1, bool onshell=1,int quad=3){
	gROOT->ProcessLine("gSystem->AddIncludePath(\"-I$ROOFITSYS/include/\")");
	gROOT->ProcessLine("gSystem->Load(\"libRooFit\")");
	gROOT->ProcessLine("gSystem->Load(\"libHiggsAnalysisCombinedLimit.so\")");
for (int iCat=0; iCat<quad; iCat++){
	dosomething(chan,cate_vbf, onshell,quad,iCat);
	}
}
