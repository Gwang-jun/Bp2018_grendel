#include "uti.h"
#include "parameters.h"
#include "TLegendEntry.h"
#include "../Systematics/systematics.C"
#include <iostream>
#include <fstream>
using namespace std;

float tpadr = 0.7;

int _nBins = nBinsCent;
double *_ptBins = ptBinsCent;

bool plotFONLL = false;
int addpbpb = 0;
TString outplotf = "plotCrossSectionCent";

void CrossSectionRatioCent(TString inputFONLL="", TString input="", TString efficiency="",TString outputplot="",int usePbPb=1,TString label="PbPb", int doDataCor = 0,double lumi=1.,Float_t centMin=0.,Float_t centMax=90.)
{
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
	gStyle->SetEndErrorSize(0);
	gStyle->SetMarkerStyle(20);

	std::cout<<"step1"<<std::endl;

	TFile* fileReference = new TFile(inputFONLL.Data());  
	TGraphAsymmErrors* gaeBplusReference = (TGraphAsymmErrors*)fileReference->Get("gaeSigmaBplus");

	if (!(usePbPb==1||usePbPb==0)) std::cout<<"ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!, you are using a non valid isPbPb option"<<std::endl;
	bool isPbPb=(bool)(usePbPb);

	float tpadpos = 1-tpadr;
	if(!isPbPb) tpadr = 1;

	//float cRightMargin = 0.1451613;
	//float cLeftMargin = 0.05443548;

	float cRightMargin = 0.1;
	float cLeftMargin = 0.15;
	float cTopMargin = 0.08474576*tpadr;
	float cBottomMargin = 0.3*(1-tpadr);

	double mod_ptBins[_nBins];
	for(int i=0;i<_nBins;i++){
	  mod_ptBins[i] = npart[_nBins-1-i];
	  cout<<mod_ptBins[i]<<endl;
	}

	double apt[_nBins];
	//bin half width
	double aptl[_nBins];
	//number of every rebined bin
	double bin_num[_nBins];

	for (int ibin=0; ibin<_nBins; ibin++){
	  apt[ibin]  = mod_ptBins[ibin];
	  aptl[ibin] = 2.5;
	}

	double pti = 0.;
	double pte = 300.;

	//TFile* file = new TFile(input.Data());
	//TFile* fileeff = new TFile(efficiency.Data());


        TFile* file1 = new TFile("unbinnedfiles/yields_Bp_Fid_pt650_Cent0-90.root");
        TFile* fileeff1 = new TFile("ROOTfiles/obsolete/MCstudiesPbPb2D_newpt_Cent0-90.root");
	TFile* fileinveff1 = new TFile("ROOTfiles/MCstudiesPbPbAverage_Fid2D_pt650_Cent0-90.root");
        TH1F* hPtSigmadiff1 = (TH1F*)file1->Get("hPt");
        TH1F* hEff1 = (TH1F*)fileeff1->Get("hEff");
	TH1F* hinvEff1 = (TH1F*)fileinveff1->Get("invEffave");
	//hPtSigmadiff1->Divide(hEff1);
	hPtSigmadiff1->Multiply(hinvEff1);
        double ptwidth[nBins];
        double coryield1 = 0.0;
        double coryieldErr1 = 0.0;

	/*
        TFile* file2 = new TFile("unbinnedfiles/yields_Bp_Fid_full_Cent30-90.root");
        TFile* fileeff2 = new TFile("unbinnedfiles/MCstudiesPbPb_newpt_Cent30-90.root");
	TFile* fileinveff2 = new TFile("ROOTfiles/MCstudiesPbPbAverage_Fid2D_inclusive_Cent30-90.root");
        TH1F* hPtSigmadiff2 = (TH1F*)file2->Get("hPt");
        TH1F* hEff2 = (TH1F*)fileeff2->Get("hEff");
	TH1F* hinvEff2 = (TH1F*)fileinveff2->Get("invEffave");
        //hPtSigmadiff2->Divide(hEff2);
	hPtSigmadiff2->Multiply(hinvEff2);
        double coryield2 = 0.0;
        double coryieldErr2 = 0.0;
	*/

	hPtSigmadiff1->Scale(44.);
	//hPtSigmadiff1->Scale(1./(2*TAA[0]*BRchain));//6.274 15.41 1.705
	//hPtSigmadiff1->Scale(1./(2*BRchain));//6.274 15.41 1.705

        for(int k=0;k<nBins;k++) //must use differential pt bins that are used to calculate inclusive pt corrected yield
          {
            ptwidth[k] = (ptBins[k+1]-ptBins[k])/(ptBins[nBins]-ptBins[0]);
	    std::cout<<"ptbin "<<k<<" "<<ptwidth[k]*hPtSigmadiff1->GetBinContent(k+1)<<" pm "<<ptwidth[k]*hPtSigmadiff1->GetBinError(k+1)<<std::endl;
	    coryield1 += ptwidth[k]*hPtSigmadiff1->GetBinContent(k+1);
            coryieldErr1 += ptwidth[k]*ptwidth[k]*hPtSigmadiff1->GetBinError(k+1)*hPtSigmadiff1->GetBinError(k+1);
            //coryield2 += ptwidth[k]*hPtSigmadiff2->GetBinContent(k+1);
            //coryieldErr2 += ptwidth[k]*ptwidth[k]*hPtSigmadiff2->GetBinError(k+1)*hPtSigmadiff2->GetBinError(k+1);
          }
	

        TH1D* hPtSigmahibin = new TH1D("hPtSigma","",nBinsCent,ptBinsCent);
        //hPtSigmahibin->SetBinContent(1,hPtSigmadiff1->GetBinContent(1));
        //hPtSigmahibin->SetBinError(1,hPtSigmadiff1->GetBinError(1));
        //hPtSigmahibin->SetBinContent(2,hPtSigmadiff1->GetBinContent(2));
        //hPtSigmahibin->SetBinError(2,hPtSigmadiff1->GetBinError(2));
	hPtSigmahibin->SetBinContent(1,coryield1);
        hPtSigmahibin->SetBinError(1,sqrt(coryieldErr1));
        //hPtSigmahibin->SetBinContent(2,coryield2);
        //hPtSigmahibin->SetBinError(2,sqrt(coryieldErr2));

        //hPtSigmahibin->Scale(1./(2*BRchain));
        for(int k=0;k<nBinsCent;k++)
          {
            //hPtSigmahibin->SetBinContent(k+1,hPtSigmahibin->GetBinContent(k+1)/TAA[k]);
            //hPtSigmahibin->SetBinError(k+1,hPtSigmahibin->GetBinError(k+1)/TAA[k]);
          }

        for(int k=0;k<nBinsCent;k++){
	printf("Cent bin %.0f-%.0f Corrected Yield: %f Err: %f\n", ptBinsCent[k], ptBinsCent[k+1], hPtSigmahibin->GetBinContent(k+1), hPtSigmahibin->GetBinError(k+1));}

	Double_t xr[_nBins], xrlow[_nBins], xrhigh[_nBins], ycross[_nBins],ycrossstat[_nBins],ycrosssysthigh[_nBins],ycrosssystlow[_nBins], yFONLL[_nBins];
	Double_t yratiocrossFONLL[_nBins], yratiocrossFONLLstat[_nBins], yratiocrossFONLLsysthigh[_nBins], yratiocrossFONLLsystlow[_nBins];
	Double_t yFONLLrelunclow[_nBins], yFONLLrelunchigh[_nBins], yunity[_nBins];

	Double_t ycrosshibin[_nBins], ycrossstathibin[_nBins], ycrosssysthighhibin[_nBins], ycrosssystlowhibin[_nBins];

	std::cout<<"step2"<<std::endl;

	for(int i=0;i<_nBins;i++)
	{
      	  //gaeBplusReference->GetPoint(i,xr[i],yFONLL[i]);
	  //xrlow[i] = gaeBplusReference->GetErrorXlow(i);
	  //xrhigh[i] = gaeBplusReference->GetErrorXhigh(i);
	  ycrosshibin[i] = hPtSigmahibin->GetBinContent(i+1);
	  ycrossstathibin[i] = hPtSigmahibin->GetBinError(i+1);
	  double systematic=0.;
	  
	  if (!isPbPb) systematic=0.01*systematicsPP(xr[i],0.);
	  else  systematic=0.01*systematicsPbPb(30.,0,_ptBins[i]/2.,_ptBins[i+1]/2.,0.);     
	  
	  //printf("bin: %d\n", i);
	  //printf("xsec value: %f\n", ycross[i]);
	  //printf("stat err: %f\n", ycrossstat[i]/ycross[i]);
	  //printf("sys  err: %f\n", systematic);
	  
	  ycrosssysthighhibin[i]= hPtSigmahibin->GetBinContent(i+1)*systematic;
	  ycrosssystlowhibin[i]= hPtSigmahibin->GetBinContent(i+1)*systematic;
      	  //yratiocrossFONLL[i] = ycross[i]/yFONLL[i];
	  //yratiocrossFONLLstat[i] = ycrossstat[i]/yFONLL[i];
	  //yratiocrossFONLLsysthigh[i] = ycrosssysthigh[i]/yFONLL[i];
	  //yratiocrossFONLLsystlow[i] = ycrosssystlow[i]/yFONLL[i];
	  //yFONLLrelunclow[i] = gaeBplusReference->GetErrorYlow(i)/yFONLL[i];
	  //yFONLLrelunchigh[i] = gaeBplusReference->GetErrorYhigh(i)/yFONLL[i];
	  //yunity[i] = yFONLL[i]/yFONLL[i];
	  
	  
	  //ofstream foutResults(Form("ResultFile/Cross_%d.tex",i));
	  //foutResults	<< "Mean = "  <<  ycross[i] << endl;
	  //foutResults	<< "StatError = "  << ycrossstat[i]/ycross[i] * 100 << "%" << endl;
	  //foutResults	<< "SystError = "  << systematic * 100 << "%" << endl;	
	  

	}

	TH1F* hPtSigmaref = new TH1F("hPtSigmaref","",_nBins,mod_ptBins);
	double hPtSigmarefcent[_nBins], hPtSigmareferror[_nBins];

        for(int i=0;i<_nBins;i++)
          {
	    hPtSigmarefcent[i]=hPtSigmahibin->GetBinContent(_nBins-i);
	    hPtSigmareferror[i]=hPtSigmahibin->GetBinError(_nBins-i);
            hPtSigmaref->SetBinContent(i+1,hPtSigmarefcent[i]);
            hPtSigmaref->SetBinError(i+1,hPtSigmareferror[i]);
	    ycross[i]=ycrosshibin[_nBins-i-1];
	    ycrossstat[i]=ycrossstathibin[_nBins-i-1];
	    ycrosssysthigh[i]=ycrosssysthighhibin[_nBins-i-1];
	    ycrosssystlow[i]=ycrosssystlowhibin[_nBins-i-1];   
          }
	
        TGraphAsymmErrors* hPtSigma = new TGraphAsymmErrors(_nBins,apt,hPtSigmarefcent,aptl,aptl,hPtSigmareferror,hPtSigmareferror);
        hPtSigma->SetName("hPtSigma");

	//TGraphAsymmErrors* gaeCrossSyst = new TGraphAsymmErrors(_nBins,xr,ycross,xrlow,xrhigh,ycrosssystlow,ycrosssysthigh);
	TGraphAsymmErrors* gaeCrossSyst = new TGraphAsymmErrors(_nBins,apt,ycross,aptl,aptl,ycrosssystlow,ycrosssysthigh);
	gaeCrossSyst->SetName("gaeCrossSyst");
	gaeCrossSyst->SetMarkerStyle(20);
	gaeCrossSyst->SetMarkerSize(0.8/tpadr);

	/*
	TGraphAsymmErrors* gaeRatioCrossFONLLstat = new TGraphAsymmErrors(_nBins,xr,yratiocrossFONLL,xrlow,xrhigh,yratiocrossFONLLstat,yratiocrossFONLLstat);
	gaeRatioCrossFONLLstat->SetName("gaeRatioCrossFONLLstat");
	gaeRatioCrossFONLLstat->SetMarkerStyle(20);
	gaeRatioCrossFONLLstat->SetMarkerSize(0.8/tpadr);

	TGraphAsymmErrors* gaeRatioCrossFONLLsyst= new TGraphAsymmErrors(_nBins,xr,yratiocrossFONLL,xrlow,xrhigh,yratiocrossFONLLsystlow,yratiocrossFONLLsysthigh);
	gaeRatioCrossFONLLsyst->SetName("gaeRatioCrossFONLLsyst");
	gaeRatioCrossFONLLsyst->SetLineWidth(2);
	gaeRatioCrossFONLLsyst->SetLineColor(1);
	gaeRatioCrossFONLLsyst->SetFillColor(5);
	gaeRatioCrossFONLLsyst->SetFillStyle(0);

	TGraphAsymmErrors* gaeRatioCrossFONLLunity = new TGraphAsymmErrors(_nBins,xr,yunity,xrlow,xrhigh,yFONLLrelunclow,yFONLLrelunchigh);
	gaeRatioCrossFONLLunity->SetName("gaeRatioCrossFONLLunity");
	gaeRatioCrossFONLLunity->SetLineWidth(2);
	gaeRatioCrossFONLLunity->SetLineColor(kOrange);
	gaeRatioCrossFONLLunity->SetFillColor(kOrange);
	gaeRatioCrossFONLLunity->SetFillStyle(1001);
	*/

	TCanvas* cSigma = new TCanvas("cSigma","",600,600);
	cSigma->SetFrameBorderMode(0);
	cSigma->SetFrameBorderMode(0);
	cSigma->Range(-1.989924,-0.2917772,25.49622,2.212202);
	cSigma->SetFillColor(0);
	cSigma->SetBorderMode(0);
	cSigma->SetBorderSize(2);
	cSigma->SetRightMargin(cRightMargin);
	cSigma->SetLeftMargin(cLeftMargin);
	cSigma->SetTopMargin(cTopMargin);
	cSigma->SetBottomMargin(cBottomMargin);
	cSigma->SetFrameBorderMode(0);
	cSigma->SetFrameBorderMode(0);
	cSigma->SetLogy();
	cSigma->cd();
	TPad* pSigma = new TPad("pSigma","",0.,tpadpos,1.,1.);
	pSigma->SetFillColor(0);
	pSigma->SetBorderMode(0);
	pSigma->SetBorderSize(2);
	pSigma->SetRightMargin(cRightMargin);
	pSigma->SetLeftMargin(cLeftMargin);
	pSigma->SetTopMargin(cTopMargin/tpadr);
	pSigma->SetBottomMargin(0);
	pSigma->SetLogy();
	if(plotFONLL){
		pSigma->Draw();
		pSigma->cd();
	}

	Float_t yaxisMin=1.e+5,yaxisMax=1.e+7;
	if(plotFONLL) yaxisMin=1.e+3;

	//TH2F* hemptySigma=new TH2F("hemptySigma","",50,_ptBins[0]-5.,_ptBins[_nBins]+5.,10.,yaxisMin,yaxisMax);  
	TH2F* hemptySigma=new TH2F("hemptySigma","",50,pti,pte,10.,yaxisMin,yaxisMax);  
	hemptySigma->GetXaxis()->CenterTitle();
	hemptySigma->GetYaxis()->CenterTitle();
	//hemptySigma->GetYaxis()->SetTitle("#frac{d#sigma_{pp}}{dp_{T}} ( pb GeV^{-1}c)");
	//if(isPbPb) hemptySigma->GetYaxis()->SetTitle("#frac{1}{T_{AA}} #frac{dN}{dp_{T}} ( pb GeV^{-1}c)");
	hemptySigma->GetYaxis()->SetTitle("#frac{1}{T_{AA}} #frac{dN}{dp_{T}} (mb*GeV^{-1}c)");
	hemptySigma->GetXaxis()->SetTitle("<N_{part}>");
        hemptySigma->GetXaxis()->SetTitleOffset(1.);
        hemptySigma->GetYaxis()->SetTitleOffset(1./tpadr);
        hemptySigma->GetXaxis()->SetTitleSize(0.12*tpadpos);
        hemptySigma->GetYaxis()->SetTitleSize(0.06*tpadr);
        hemptySigma->GetXaxis()->SetTitleFont(42);
        hemptySigma->GetYaxis()->SetTitleFont(42);
	hemptySigma->GetXaxis()->SetLabelFont(42);
        hemptySigma->GetYaxis()->SetLabelFont(42);
        hemptySigma->GetXaxis()->SetLabelSize(0.12*tpadpos);
        hemptySigma->GetYaxis()->SetLabelSize(0.06*tpadr);
        hemptySigma->GetXaxis()->SetLabelOffset(0.005*tpadpos);
        hemptySigma->SetMaximum(2);
        hemptySigma->SetMinimum(0.);
        hemptySigma->Draw();

        gaeBplusReference->SetFillColor(kOrange);
        gaeBplusReference->SetFillStyle(1001);
        gaeBplusReference->SetLineWidth(3);
        gaeBplusReference->SetLineColor(kOrange);
        //gaeBplusReference->Draw("5same");
        //if(!isPbPb)gaeBplusReference->Draw("5same");
        hPtSigma->SetLineColor(1);
        hPtSigma->SetLineWidth(2);
        hPtSigma->SetMarkerStyle(20);
        hPtSigma->SetMarkerSize(1.2*tpadr);
        hPtSigma->Draw("epsame");
        gaeCrossSyst->SetFillColor(1);
        gaeCrossSyst->SetFillStyle(0);
        gaeCrossSyst->SetLineWidth(2);
        gaeCrossSyst->SetLineColor(1);
        gaeCrossSyst->Draw("5same");

        //TLatex* texCms = new TLatex(0.16,0.95, "#scale[1.25]{CMS}");
        TLatex* texCms = new TLatex(0.52,1-(1-0.88)*tpadr, "CMS");
        texCms->SetNDC();
        texCms->SetTextAlign(13);
        texCms->SetTextSize(0.08*tpadr);
        texCms->SetTextFont(62);
        texCms->Draw();

        //TLatex* texPrel = new TLatex(0.20,0.95, "#scale[1.25]Preliminary");
        TLatex* texPrel = new TLatex(0.30,0.89, "Preliminary");
        texPrel->SetNDC();
        texPrel->SetTextAlign(13);
        texPrel->SetTextSize(0.050);
        texPrel->SetTextFont(52);
        //texPrel->Draw();

        TString text;
        if (label=="PbPb") { text="1.5 nb^{-1} (PbPb 5.02 TeV)";}
        else {text="28.0 pb^{-1} (pp 5.02 TeV)";}
        TLatex* texlumi = new TLatex(0.90,1-(1-0.936)*tpadr,text.Data());
        texlumi->SetNDC();
        texlumi->SetTextAlign(31);
        texlumi->SetTextFont(42);
        texlumi->SetTextSize(0.050*tpadr);
        texlumi->SetLineWidth(2);
        texlumi->Draw();

        TLatex* texCol = new TLatex(0.94,0.95, Form("%s #sqrt{s_{NN}} = 5.02 TeV",label.Data()));
        texCol->SetNDC();
        texCol->SetTextAlign(32);
        texCol->SetTextSize(0.04*tpadr);
        texCol->SetTextFont(42);
        //texCol->Draw();

        TString texper="%";
        TLatex* texCent = new TLatex(0.53,0.600,Form("Cent. %.0f-%.0f%s",centMin,centMax,texper.Data()));
        texCent->SetNDC();
        texCent->SetTextFont(42);
        texCent->SetTextSize(0.04);
        //if(isPbPb) texCent->Draw();

        TLatex* texY = new TLatex(0.53,1-(1-0.65)*tpadr,"|y| < 2.4");
        texY->SetNDC();
        texY->SetTextFont(42);
        texY->SetTextSize(0.05*tpadr);
        texY->SetLineWidth(2);
        texY->Draw();

        //TLatex* texB = new TLatex(0.77,0.82,"B^{#plus}+B^{#minus}");
        TLatex* texB = new TLatex(0.70,1-(1-0.82)*tpadr,"B^{+}");
        texB->SetNDC();
        texB->SetTextFont(62);
        texB->SetTextSize(0.09*tpadr);
        texB->SetLineWidth(2);
        texB->Draw();

        TLatex* texGlobal = new TLatex(0.53,0.59,Form("Global uncert. %.1f%s",normalizationUncertaintyForPP(),texper.Data()));
        if(isPbPb)  texGlobal = new TLatex(0.53,1-(1-0.59)*tpadr,Form("Global uncert. #plus%.1f, #minus%.1f%s",normalizationUncertaintyForPbPb(1),normalizationUncertaintyForPbPb(0),texper.Data()));
        texGlobal->SetNDC();
        texGlobal->SetTextFont(42);
        texGlobal->SetTextSize(0.05*tpadr);
        texGlobal->SetLineWidth(2);
        if(!addpbpb) texGlobal->Draw();

        TLegend* leg_CS = new TLegend(0.52,1-(1-0.70)*tpadr,0.85,1-(1-0.80)*tpadr);
        leg_CS->SetBorderSize(0);
        leg_CS->SetFillStyle(0);
        leg_CS->SetTextSize(0.05*tpadr);
        leg_CS->AddEntry(hPtSigma,"Data","pf");
        //leg_CS->AddEntry(gaeBplusReference,"FONLL pp ref.","f");//PAS
        //if(!isPbPb) leg_CS->AddEntry(gaeBplusReference,"FONLL","f");//paper
        leg_CS->Draw("same");

	/*
	if(!addpbpb){
		leg_CS->AddEntry(hPtSigma,"Data","pf");
		if(plotFONLL) leg_CS->AddEntry(gaeBplusReference,"FONLL","f");//paper
		leg_CS->Draw("same");
	//	texGlobal->Draw();
	//	texGlobal_num->Draw();
	}
	else{
		TFile* filepbpb = new TFile("ROOTfiles/CrossSectionPbPb_CENT.root");
		TGraphAsymmErrors* gaeCrossSyst_PbPb = (TGraphAsymmErrors*)filepbpb->Get("gaeCrossSyst");
		gaeCrossSyst_PbPb->SetLineColor(2);
		gaeCrossSyst_PbPb->Draw("5same");
		TH1F* hPtSigma_PbPb = (TH1F*)filepbpb->Get("hPtSigma");
		hPtSigma_PbPb->SetLineColor(2);
		hPtSigma_PbPb->SetLineWidth(2);
		hPtSigma_PbPb->SetMarkerColor(2);
		hPtSigma_PbPb->SetMarkerStyle(21);
		hPtSigma_PbPb->SetMarkerSize(1.2/tpadr);
		hPtSigma_PbPb->Draw("epsame");
		leg_CS->SetX1(0.42);
		leg_CS->SetX2(0.80);
		leg_CS->SetY2(1-(1-0.88)/tpadr);
		leg_CS->AddEntry(hPtSigma,"pp","pf");
		leg_CS->AddEntry(hPtSigma_PbPb,"PbPb Cent. 0-100%","pf");
		if(plotFONLL) leg_CS->AddEntry(gaeBplusReference,"FONLL pp ref.","f");//PAS
		leg_CS->Draw("same");
		hemptySigma->GetYaxis()->SetTitle("#frac{d#sigma_{pp}}{dp_{T}} or #frac{1}{T_{AA}} #frac{dN_{PbPb}}{dp_{T}} ( pb GeV^{-1}c)");
		hemptySigma->GetYaxis()->SetTitleSize(0.045/tpadr);
		hemptySigma->GetYaxis()->SetTitleOffset(1.7*tpadr);
		uncGlobal_pp = Form("pp:     #pm%.1f%s",normalizationUncertaintyForPP(),texper.Data());
		uncGlobal_PbPb = Form("PbPb: #plus%.1f, #minus%.1f%s",normalizationUncertaintyForPbPb(1),normalizationUncertaintyForPbPb(0),texper.Data());
		texGlobal = new TLatex(0.53,0.594,Form("Global uncertainty"));
		TLatex* texGlobal_pp = new TLatex(0.53,1-(1-0.54)/tpadr,uncGlobal_pp);
		TLatex* texGlobal_PbPb = new TLatex(0.53,1-(1-0.48)/tpadr,uncGlobal_PbPb);
		texGlobal->SetNDC();
		texGlobal->SetTextFont(42);
		texGlobal->SetTextSize(0.05*tpadr);
		texGlobal->SetLineWidth(2);
		//texGlobal->Draw();
		texGlobal_pp->SetNDC();
		texGlobal_pp->SetTextFont(42);
		texGlobal_pp->SetTextSize(0.05*tpadr);
		texGlobal_pp->SetLineWidth(2);
		texGlobal_pp->Draw();
		texGlobal_PbPb->SetNDC();
		texGlobal_PbPb->SetTextFont(42);
		texGlobal_PbPb->SetTextSize(0.05*tpadr);
		texGlobal_PbPb->SetLineWidth(2);
		//texGlobal_PbPb->Draw();
	}
	*/

	cSigma->cd();
	TPad* pRatio = new TPad("pRatio","",0.,0.,1.,tpadpos);
	pRatio->SetRightMargin(cRightMargin);
	pRatio->SetLeftMargin(cLeftMargin);
	pRatio->SetTopMargin(0);
	pRatio->SetBottomMargin(cBottomMargin/(1-tpadr));//0.25

	TH2F* hemptyRatio=new TH2F("hemptyRatio","",50,_ptBins[0]-5.,_ptBins[_nBins]+5.,10.,0.2,1.8);//paper
	hemptyRatio->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	hemptyRatio->GetXaxis()->CenterTitle();
	//hemptyRatio->GetYaxis()->CenterTitle();
	//hemptyRatio->GetYaxis()->SetTitle("Data/FONLL");
	hemptyRatio->GetYaxis()->SetTitle("#frac{Data}{FONLL}");
	hemptyRatio->GetXaxis()->SetTitleOffset(1.);
	hemptyRatio->GetYaxis()->SetTitleOffset(1.2);
	hemptyRatio->GetYaxis()->SetTitleOffset(1.4*(1-tpadr));
	hemptyRatio->GetXaxis()->SetTitleSize(0.055/(1-tpadr));
	hemptyRatio->GetYaxis()->SetTitleSize(0.055/(1-tpadr));
	hemptyRatio->GetXaxis()->SetTitleFont(42);
	hemptyRatio->GetYaxis()->SetTitleFont(42);
	hemptyRatio->GetXaxis()->SetLabelFont(42);
	hemptyRatio->GetYaxis()->SetLabelFont(42);
	hemptyRatio->GetXaxis()->SetLabelSize(0.048/(1-tpadr));
	hemptyRatio->GetYaxis()->SetLabelSize(0.048/(1-tpadr));
	hemptyRatio->GetYaxis()->SetNdivisions(505);
	hemptyRatio->GetXaxis()->SetTickLength(0.03/tpadpos);

	TLine* l = new TLine(_ptBins[0]-5.,1,_ptBins[_nBins]+5.,1);
	l->SetLineWidth(1);
	l->SetLineStyle(2);

	if(plotFONLL){
		pRatio->Draw();
		pRatio->cd();
		hemptyRatio->Draw();
		l->Draw("same");
       		//gaeRatioCrossFONLLunity->Draw("5same");
		//gaeRatioCrossFONLLstat->Draw("epsame");
		//gaeRatioCrossFONLLsyst->Draw("5same");
	}

	if(addpbpb){
		hemptyRatio->GetYaxis()->SetTitle("pp / FONLL");
	}

	TString _postfix = "";
	if(doDataCor==1) _postfix += "_EFFCOR";
	if(addpbpb) _postfix += "_AddPbPb";
	if(!isPbPb) cSigma->SaveAs(Form("%s/canvasSigmaBplusRatio%s%s.pdf",outplotf.Data(),label.Data(),_postfix.Data()));
	else cSigma->SaveAs(Form("%s/canvasSigmaBplusRatio%s_%.0f_%.0f%s.pdf",outplotf.Data(),label.Data(),centMin,centMax,_postfix.Data()));
	if(!isPbPb) cSigma->SaveAs(Form("%s/canvasSigmaBplusRatio%s%s.png",outplotf.Data(),label.Data(),_postfix.Data()));
	else cSigma->SaveAs(Form("%s/canvasSigmaBplusRatio%s_%.0f_%.0f%s.png",outplotf.Data(),label.Data(),centMin,centMax,_postfix.Data()));
	if(!isPbPb) cSigma->SaveAs(Form("%s/canvasSigmaBplusRatio%s%s.C",outplotf.Data(),label.Data(),_postfix.Data()));
	else cSigma->SaveAs(Form("%s/canvasSigmaBplusRatio%s_%.0f_%.0f%s.C",outplotf.Data(),label.Data(),centMin,centMax,_postfix.Data()));

	TCanvas* cEff = new TCanvas("cEff","",550,500);
	TH2F* hemptyEff=new TH2F("hemptyEff","",50,0.,110.,10.,0,1.);  
	hemptyEff->GetXaxis()->CenterTitle();
	hemptyEff->GetYaxis()->CenterTitle();
	hemptyEff->GetYaxis()->SetTitle("#alpha x #epsilon_{reco} x #epsilon_{sel} ");
	hemptyEff->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	hemptyEff->GetXaxis()->SetTitleOffset(1.);
	hemptyEff->GetYaxis()->SetTitleOffset(1.4);
	hemptyEff->GetXaxis()->SetTitleSize(0.045);
	hemptyEff->GetYaxis()->SetTitleSize(0.045);
	hemptyEff->GetXaxis()->SetTitleFont(42);
	hemptyEff->GetYaxis()->SetTitleFont(42);
	hemptyEff->GetXaxis()->SetLabelFont(42);
	hemptyEff->GetYaxis()->SetLabelFont(42);
	hemptyEff->GetXaxis()->SetLabelSize(0.04);
	hemptyEff->GetYaxis()->SetLabelSize(0.04);  
	hemptyEff->SetMaximum(2);
	hemptyEff->SetMinimum(0.);
	hemptyEff->Draw();
	cEff->cd();
	hemptyEff->Draw();
	//hEff->SetLineWidth(2);
	//hEff->SetLineColor(1);
	//hEff->SetMarkerStyle(20);
	//hEff->SetMarkerSize(1.2);
	//hEff->Draw("same");
	if(!isPbPb) cEff->SaveAs(Form("%s/efficiency%s%s.pdf",outplotf.Data(),label.Data(),_postfix.Data()));
	else cEff->SaveAs(Form("%s/efficiency%s_%.0f_%.0f%s.pdf",outplotf.Data(),label.Data(),centMin,centMax,_postfix.Data()));
	if(!isPbPb) cEff->SaveAs(Form("%s/efficiency%s%s.png",outplotf.Data(),label.Data(),_postfix.Data()));
	else cEff->SaveAs(Form("%s/efficiency%s_%.0f_%.0f%s.png",outplotf.Data(),label.Data(),centMin,centMax,_postfix.Data()));

	TFile *outputfile=new TFile(outputplot.Data(),"recreate");
	outputfile->cd();
	//gaeCrossSyst->Write();
	//gaeRatioCrossFONLLstat->Write();
	//gaeBplusReference->Write();
	//hPtSigma->Write();
	hPtSigmahibin->Write();
	//gaeRatioCrossFONLLstat->Write();
	//gaeRatioCrossFONLLsyst->Write();
	//gaeRatioCrossFONLLunity->Write();
	//hEff->Write();

}

int main(int argc, char *argv[])
{
  if(argc==11)
    {
      CrossSectionRatioCent(argv[1], argv[2], argv[3],argv[4],atoi(argv[5]),argv[6],atoi(argv[7]),atof(argv[8]),atof(argv[9]),atof(argv[10]));
      return 0;
    }
  else if(argc==9)
    {
      CrossSectionRatioCent(argv[1], argv[2], argv[3],argv[4],atoi(argv[5]),argv[6],atoi(argv[7]),atof(argv[8]));
      return 0;
    }
  else
    {
      std::cout << "Wrong number of inputs" << std::endl;
      return 1;
    }

}
