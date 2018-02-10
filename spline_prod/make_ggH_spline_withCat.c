using namespace RooFit;
using namespace std;

//RooFormulaVar* extract_Formula(std::vector<RooRealVar*> projvars, RooWorkspace* wo, char *tag= "bkg");
RooPlot* quickplot(RooWorkspace* w, std::vector<RooAbsPdf*> pdfs, vector<RooRealVar*> projvars, unsigned int chosen, vector<float> varset, vector<unsigned int> ndim);
void plotOnC(TCanvas* c, RooPlot* plot, TString tag);
RooArgList extract_pdf(RooWorkspace* w);
void setVar(RooWorkspace* w, const int var, const float Val);
RooFormulaVar* writeFormula(RooWorkspace* w, RooRealVar* r, TString tag,int cat, int cate_vbf, int prod_cate);

void do_make_ggH_spline_withCat(TString tag, int cat,int quad, TString workdir="./rpdfWS_withCat"){

	RooWorkspace* wo, *ws;
	TFile* fm, *f, *fr1;
	vector<RooAbsReal*> vec_pdfs;
	RooAbsReal* temppdfr1;
	RooFormulaVar* tempcoeff;

	RooFormulaVar *ggH_norm;
	TString pdfname, tmplname;
	RooRealVar *r= new RooRealVar("r","signal strength",1.,0.000,1000);
	RooConstVar *constone = new RooConstVar("const1","const1",1.);
	TH2F* th2ftmpl;
	RooRealVar *dbkg = new RooRealVar("dbkg_kin","dbkg_kin",0.,1.);
	dbkg->setBins(10);

	RooArgSet argset;
    	fm = new TFile("../signal_preparation/all_sig_integral.root","read");
	ws = (RooWorkspace*) fm->Get("w");
	ws->exportToCint("ws");
	//ws->Print();
	TFile ftest(Form("spline_WS_withCat_quad%d/%s_spline_WS_cat%d.root",quad,tag.Data(),cat),"recreate");
	RooWorkspace* newWS = new RooWorkspace("w", "");
	newWS->addClassDeclImportDir("../interface/");
	newWS->importClassCode(RooNCSplineCore::Class(),true);
	newWS->importClassCode(RooNCSpline_2D_fast::Class(), true);
	newWS->importClassCode(RooNCSpline_3D_fast::Class(), true);
	TFile* ftmpl = new TFile("../TemplateBuilder/Normalized_spline_tmpl.root");
	RooWorkspace* wtmpl = (RooWorkspace*) ftmpl->Get("w");	
 	for(int cate_vbf=0; cate_vbf<3; cate_vbf++){
 		for(int prod_cate=0; prod_cate<3; prod_cate++){

	char* inputDir = "../signal_preparation";
 	
	cout<< inputDir<<endl;
        char* tag2 = (char*) tag.Data();
        TString tag1="";
        if(cate_vbf==0)         tag1="0jet_";
        else if(cate_vbf==1)    tag1="1jet_";
        else if(cate_vbf==2)    tag1="2jet_";
        TString prodtag="";
        if(prod_cate==0)        prodtag="";
        else if(prod_cate==1)   prodtag="vbf_";
        else if(prod_cate==2)   prodtag="vh_";

//	TString workdir="./rpdfWS_withCat";

	fr1 = new TFile(Form("%s/%s_rpdfWS_cat%d_ggH_%s.root",workdir.Data(), tag.Data(), cat,tag1.Data()));

	RooWorkspace* wr1 = (RooWorkspace*) fr1->Get("w");
	
	if(prod_cate==0)        prodtag="ggH_";
	f = new TFile(Form("%s/%s_rpdfWS_cat%d_%s%s.root",workdir.Data(), tag.Data(), cat,prodtag.Data(),tag1.Data()));
	//4e_rpdfWS_cat0.root
	wo = (RooWorkspace*) f->Get("w");
	cout<<wr1<<endl;


	
    if(cate_vbf==0)         tag1="_untagged";
    else if(cate_vbf==1)    tag1="_vbf2j";
    else if(cate_vbf==2)    tag1="_vhhad";
	if(prod_cate==0)		pdfname="ggH";
	else if (prod_cate==1)		pdfname="qqH";
	else if (prod_cate==2)		pdfname="VH";
	pdfname += tag1;
	RooRealVar* mreco = (RooRealVar*) wo->var("mreco");
	argset.add(*mreco);
	argset.add(*dbkg);
	temppdfr1 = (RooAbsReal*) wo->obj("r1");	
	temppdfr1->SetName(Form("r1_%s_cat%d_%s%s",tag.Data(),cat,pdfname.Data(),tag.Data()));
	tempcoeff = writeFormula(ws,r, tag,cat,cate_vbf,prod_cate);
	
	//ggH_norm = (RooFormulaVar*) ws->obj(Form("%s%s_norm",pdfname.Data(),tag1.Data()));

	RooRealFlooredSumPdf *temppdf = new RooRealFlooredSumPdf(Form("r1_%s_cat%d_%s%s",tag.Data(),cat,pdfname.Data(),tag.Data()),"",RooArgList(*temppdfr1),RooArgList(*constone));	

//**************** Extract template from root file*************************

	if(prod_cate==0)		tmplname="ggzz";
	else				tmplname="vbf";
	cout<<Form("%s_sig_%s",tmplname.Data(),tag.Data())<<endl;	
	
	TH2F* th2ftmpl = (TH2F*) ftmpl->Get(Form("%s_sig_%s",tmplname.Data(),tag.Data()));
	RooDataHist *th2fdatahist = new RooDataHist(Form("datahist_%s_sig_%s",tmplname.Data(),tag.Data()),"",argset,th2ftmpl);
	RooHistPdf *th2fhistpdf = new RooHistPdf(Form("histpdf_%s_sig_%s",tmplname.Data(),tag.Data()),"",argset,*th2fdatahist);
//	RooProdPdf* tempprodpdf = new RooProdPdf("tempprodpdf","tempprodpdf",*temppdf, *tmplfuncpdf);

	temppdf->SetNameTitle(Form("r1_%s",pdfname),Form("r1_%s",pdfname));	
	
//	ggH = new RooRealSumPdf(pdfname.Data(), pdfname.Data(), RooArgList(*tempprodpdf), RooArgList(*tempcoeff)); 
	RooProdPdf *ggH = new RooProdPdf(Form("%s_prodpdf",pdfname.Data()),"", *temppdf, Conditional(*th2fhistpdf,*dbkg));
//	RooRealFlooredSumPdf *ggH = new RooRealFlooredSumPdf(pdfname.Data(),"",RooArgList(*ggH_prod),RooArgList(*ggH_norm));
		
	ggH->SetNameTitle(pdfname.Data(), pdfname.Data());
	ggH_norm = new RooFormulaVar(Form("%s_norm",pdfname.Data()),Form("%s_norm",pdfname.Data()),"@0",RooArgList(*tempcoeff));	

	ftest.cd();
	newWS->import(*ggH,RooFit::RecycleConflictNodes());
	
	cout<<ggH->GetName()<<endl;
	cout<<ggH_norm->GetName()<<endl;
	newWS->import(*ggH_norm,RooFit::RecycleConflictNodes());
	
	f->Close();

	}
}

	ftest.WriteTObject(newWS);
	ftest.Close();
	ftmpl->Close();
	fm->Close();
	newWS->Print();
	return;


}

RooFormulaVar* writeFormula(RooWorkspace* w, RooRealVar* r, TString tag,int cat, int cate_vbf, int prod_cate){
	// The following is how the coefficient is derived
	//
	//f0 = bkg_integral2e2mu0Nominal;
	//f1 = bkg_integral2e2mu0Nominal+sig_integral2e2mu0Nominal+int_integral2e2mu0Nominal;
	//f2 = bkg_integral2e2mu0Nominal+4*sig_integral2e2mu0Nominal+2*int_integral2e2mu0Nominal;

	//f0*r0 = b;
	//f1*r1 = b+s+I;
	//f2*r2 = b+4s+2I;

	//s = ( r2*f2 - 2*r1*f1 +   r0*f0)/2.;
	//I = (-r2*f2 + 4*r1*f1 - 3*r0*f0)/2.;

	//Pdf = s *r + b+ I*sqrt(r);

	//Pdf = r0 * f0 * (r/2. - 3./2.*sqrt(r)+1) + r1*f1 *(-r+2*sqrt(r)) + r2*f2*(r/2.-sqrt(r)/2.); 
	

//	RooRealVar* r = w->var("r");	
	//cout<<prod_cate<<endl;
	RooAbsReal* bkg_integral2e2mu0Nominal;
//	RooConstVar* bkg_integral2e2mu0Nominal;
	RooAbsReal* sig_integral2e2mu0Nominal;
	RooAbsReal* int_integral2e2mu0Nominal;
	RooFormulaVar* rvbf;
	if(prod_cate==0){
//		bkg_integral2e2mu0Nominal= (RooAbsReal*) w->obj(Form("bkg_integral%s%d_%djetNominal",tag.Data(),cat,cate_vbf)); 
		sig_integral2e2mu0Nominal= (RooAbsReal*) w->obj(Form("sig_integral%s%d_%djetNominal",tag.Data(),cat,cate_vbf)); 
//		int_integral2e2mu0Nominal= (RooAbsReal*) w->obj(Form("int_integral%s%d_%djetNominal",tag.Data(),cat,cate_vbf)); 
		r = (RooRealVar*) w->var("r");
	}
	else if (prod_cate==1){
		cout <<Form("vbfsig_integral%s%d_%djet",tag.Data(),cat,cate_vbf)<<endl;
                sig_integral2e2mu0Nominal= (RooAbsReal*) w->obj(Form("vbfsig_integral%s%d_%djet",tag.Data(),cat,cate_vbf)); 
//vbfsig_integral4e3_0jet
    //            int_integral2e2mu0Nominal= (RooAbsReal*) w->obj(Form("vbfint_integral%s%d_%djet",tag.Data(),cat,cate_vbf)); 	
		rvbf = (RooFormulaVar*) w->obj("rvbf");
	}
	else if (prod_cate==2){
//		bkg_integral2e2mu0Nominal= new RooConstVar(Form("bkg_integral%s%d_%djet",tag.Data(),cat,cate_vbf),Form("bkg_integral%s%d_%djet",tag.Data(),cat,cate_vbf),0);
     //           bkg_integral2e2mu0Nominal= (RooAbsReal*) w->obj(Form("bkg_integral%s%d_%djet",tag.Data(),cat,cate_vbf)); 
                sig_integral2e2mu0Nominal= (RooAbsReal*) w->obj(Form("vhsig_integral%s%d_%djet",tag.Data(),cat,cate_vbf)); 
       //         int_integral2e2mu0Nominal= (RooAbsReal*) w->obj(Form("vhint_integral%s%d_%djet",tag.Data(),cat,cate_vbf)); 
		rvbf = (RooFormulaVar*) w->obj("rvbf");
	}
	//cout<<"bkg_integral4e0Nominal:   "<<bkg_integral2e2mu0Nominal->getVal()<<endl;
	RooFormulaVar* f1= new RooFormulaVar("f1","f1", "@0", RooArgList(*sig_integral2e2mu0Nominal));

	f1->SetNameTitle(Form("f1%s_cat%d_vbf%d_prod%d",tag.Data(),cat,cate_vbf,prod_cate),"f1"+tag);

	RooFormulaVar* f1_final;


	if(prod_cate==0){
	f1_final= new RooFormulaVar("f1_final","f1_final", "@1*@0", RooArgList(*f1,*r));
	f1_final->SetNameTitle(Form("f1_final%s_cat%d_vbf%d_prod%d",tag.Data(),cat,cate_vbf,prod_cate),Form("f1_final%s_cat%d",tag.Data(),cat));
}	else{
        f1_final= new RooFormulaVar("f1_final","f1_final", "@0*@1", RooArgList(*f1,*rvbf));
        f1_final->SetNameTitle(Form("f1_final%s_cat%d_vbf%d_prod%d",tag.Data(),cat,cate_vbf,prod_cate),Form("f1_final%s_cat%d",tag.Data(),cat));
}
	return f1_final;

//	form_list.add(*f1_final);
//	return form_list;
}


RooArgList extract_pdf(RooWorkspace* w){
	RooArgList pdf_list;
	for(int rval=0; rval<3; rval++){
		TString pdfname = Form("r%d",rval);
		RooAbsPdf* pdf = w->pdf(pdfname);
		pdf_list.add(*pdf);
	}
	return pdf_list;
}

void plotOnC(TCanvas* c, RooPlot* plot, TString tag){
	c->SetGrid(0,1);
	c->cd();
	plot->Draw();
	c->SaveAs("~/www/test_spline_signal_"+tag+".png");
}

RooPlot* quickplot(RooWorkspace* w, std::vector<RooAbsPdf*> pdfs, vector<RooRealVar*> projvars, unsigned int chosen, vector<float> varset, vector<unsigned int> ndim){
	RooPlot* plot = projvars[chosen]->frame(ndim[chosen]-1);
	plot->GetXaxis()->CenterTitle();
    plot->GetYaxis()->SetTitleOffset(1.2);
    plot->GetYaxis()->CenterTitle();
    plot->GetXaxis()->SetTitle("M_{4l} (GeV)");//projvars.at(1)->GetTitle());
    plot->GetYaxis()->SetTitle("H4l amplitude");
    //plot->GetYaxis()->SetRangeUser(0,0.3);
    plot->GetXaxis()->SetNdivisions(-505);
    for(unsigned int i=0; i<projvars.size(); i++){
    	if (i==chosen)	continue;
    	setVar(w,i,varset[i]);
    }
    for(unsigned int k=0; k<pdfs.size(); k++){
	    pdfs.at(k)->plotOn(plot,LineColor(k), LineWidth(2), LineStyle(2));
	}
	return plot;
}



void setVar(RooWorkspace* w, const int var, const float Val){
		TString curvar = "";
		if(var==0)	curvar="mean_pole";
		if(var==1)	curvar="mreco";
		if(var==2)	curvar="sigma_pole";
		w->var(curvar)->setConstant(false);
		w->var(curvar)->setVal(Val);
		w->var(curvar)->setConstant(true);
}

void make_ggH_spline_withCat(TString tag, int cat,int quad, TString workdir="./rpdfWS_withCat"){
  
  do_make_ggH_spline_withCat( tag, cat, quad, workdir);

}


