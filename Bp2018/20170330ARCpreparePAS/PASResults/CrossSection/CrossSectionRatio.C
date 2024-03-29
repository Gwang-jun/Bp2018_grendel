#include "uti.h"
#include "parameters.h"
#include "TLegendEntry.h"
#include "../Systematics/systematics.C"

void CrossSectionRatio(TString inputFONLL="ROOTfiles/output_inclusiveDd0meson_5TeV_y1.root", TString input="ROOTfiles/hPtSpectrumBplusPP.root", TString efficiency="test.root",TString outputplot="myplot.root",int usePbPb=1,TString label="PbPb",int doDataCor = 0,double lumi=1.,Float_t centMin=0., Float_t centMax=100.)
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

	TFile* file = new TFile(input.Data());  
	TFile* fileeff = new TFile(efficiency.Data());
	TH1F* hEff = (TH1F*)fileeff->Get("hEff");

	TH1F* hPtSigma = (TH1F*)file->Get("hPt");
	if(doDataCor != 1) hPtSigma->Divide(hEff);
	hPtSigma->Scale(1./(2*lumi*BRchain));
	hPtSigma->SetName("hPtSigma");

	Double_t xr[nBins], xrlow[nBins], xrhigh[nBins], ycross[nBins],ycrossstat[nBins],ycrosssysthigh[nBins],ycrosssystlow[nBins], yFONLL[nBins];
	Double_t yratiocrossFONLL[nBins], yratiocrossFONLLstat[nBins], yratiocrossFONLLsysthigh[nBins], yratiocrossFONLLsystlow[nBins];
	Double_t yFONLLrelunclow[nBins], yFONLLrelunchigh[nBins], yunity[nBins];

	std::cout<<"step2"<<std::endl;

	for(int i=0;i<nBins;i++)
	{
		gaeBplusReference->GetPoint(i,xr[i],yFONLL[i]);
		xrlow[i] = gaeBplusReference->GetErrorXlow(i);
		xrhigh[i] = gaeBplusReference->GetErrorXhigh(i);
		ycross[i] = hPtSigma->GetBinContent(i+1);
		ycrossstat[i] = hPtSigma->GetBinError(i+1);
		double systematic=0.;
		if (!isPbPb) systematic=0.01*systematicsPP(xr[i],0.);
		else  systematic=0.01*systematicsPbPb(xr[i],centMin,centMax,0.);     

		ycrosssysthigh[i]= hPtSigma->GetBinContent(i+1)*systematic;
		ycrosssystlow[i]= hPtSigma->GetBinContent(i+1)*systematic;
		yratiocrossFONLL[i] = ycross[i]/yFONLL[i];
		yratiocrossFONLLstat[i] = ycrossstat[i]/yFONLL[i];
		yratiocrossFONLLsysthigh[i] = ycrosssysthigh[i]/yFONLL[i];
		yratiocrossFONLLsystlow[i] = ycrosssystlow[i]/yFONLL[i];
		yFONLLrelunclow[i] = gaeBplusReference->GetErrorYlow(i)/yFONLL[i];
		yFONLLrelunchigh[i] = gaeBplusReference->GetErrorYhigh(i)/yFONLL[i];
		yunity[i] = yFONLL[i]/yFONLL[i];
	}

	std::cout<<"step3"<<std::endl;


	TGraphAsymmErrors* gaeCrossSyst = new TGraphAsymmErrors(nBins,xr,ycross,xrlow,xrhigh,ycrosssystlow,ycrosssysthigh);
	gaeCrossSyst->SetName("gaeCrossSyst");
	gaeCrossSyst->SetMarkerStyle(20);
	gaeCrossSyst->SetMarkerSize(0.8);

	TGraphAsymmErrors* gaeRatioCrossFONLLstat = new TGraphAsymmErrors(nBins,xr,yratiocrossFONLL,xrlow,xrhigh,yratiocrossFONLLstat,yratiocrossFONLLstat);
	gaeRatioCrossFONLLstat->SetName("gaeRatioCrossFONLLstat");
	gaeRatioCrossFONLLstat->SetMarkerStyle(20);
	gaeRatioCrossFONLLstat->SetMarkerSize(0.8);

	TGraphAsymmErrors* gaeRatioCrossFONLLsyst= new TGraphAsymmErrors(nBins,xr,yratiocrossFONLL,xrlow,xrhigh,yratiocrossFONLLsystlow,yratiocrossFONLLsysthigh);
	gaeRatioCrossFONLLsyst->SetName("gaeRatioCrossFONLLsyst");
	gaeRatioCrossFONLLsyst->SetLineWidth(2);
	gaeRatioCrossFONLLsyst->SetLineColor(1);
	gaeRatioCrossFONLLsyst->SetFillColor(5);//Color(2);
	gaeRatioCrossFONLLsyst->SetFillStyle(0);

	TGraphAsymmErrors* gaeRatioCrossFONLLunity = new TGraphAsymmErrors(nBins,xr,yunity,xrlow,xrhigh,yFONLLrelunclow,yFONLLrelunchigh);
	gaeRatioCrossFONLLunity->SetName("gaeRatioCrossFONLLunity");
	gaeRatioCrossFONLLunity->SetLineWidth(2);
	//gaeRatioCrossFONLLunity->SetLineColor(5);//Color(2);
	//gaeRatioCrossFONLLunity->SetFillColor(5);//Color(2);
	gaeRatioCrossFONLLunity->SetLineColor(kOrange);//Color(2);
	gaeRatioCrossFONLLunity->SetFillColor(kOrange);//Color(2);
	//gaeRatioCrossFONLLunity->SetFillStyle(3001);
	gaeRatioCrossFONLLunity->SetFillStyle(1001);

	TCanvas* cSigma = new TCanvas("cSigma","",600,750);
	cSigma->SetFrameBorderMode(0);
	cSigma->SetFrameBorderMode(0);
	cSigma->Range(-1.989924,-0.2917772,25.49622,2.212202);
	cSigma->SetFillColor(0);
	cSigma->SetBorderMode(0);
	cSigma->SetBorderSize(2);
	cSigma->SetLeftMargin(0.1451613);
	cSigma->SetRightMargin(0.05443548);
	cSigma->SetTopMargin(0.08474576);
	cSigma->SetBottomMargin(0.1165254);
	cSigma->SetFrameBorderMode(0);
	cSigma->SetFrameBorderMode(0);
	cSigma->cd();
	TPad* pSigma = new TPad("pSigma","",0,0.25,1,1);
	pSigma->SetFillColor(0);
	pSigma->SetBorderMode(0);
	pSigma->SetBorderSize(2);
	pSigma->SetLeftMargin(0.1451613);
	pSigma->SetRightMargin(0.05443548);
	pSigma->SetTopMargin(0.08474576);
	pSigma->SetBottomMargin(0);
	pSigma->SetLogy();
	pSigma->Draw();
	pSigma->cd();


	Float_t yaxisMin=1.1,yaxisMax=1.e+9;//PAS
//	Float_t yaxisMin=1.e3,yaxisMax=1.e+7;//paper 20170224
	if(label=="PPMB"||label=="PbPbMB")
	{
		yaxisMin=1.e+4;
		yaxisMax=1.e+13;
	}
	TH2F* hemptySigma=new TH2F("hemptySigma","",50,ptBins[0]-5.,ptBins[nBins]+5.,10.,yaxisMin,yaxisMax);  
	hemptySigma->GetXaxis()->CenterTitle();
	hemptySigma->GetYaxis()->CenterTitle();
	hemptySigma->GetYaxis()->SetTitle("d#sigma / dp_{T}( pb GeV^{-1}c)");
	if(isPbPb) hemptySigma->GetYaxis()->SetTitle("1/T_{AA} * dN / dp_{T}( pb GeV^{-1}c)");
	hemptySigma->GetXaxis()->SetTitleOffset(1.);
	hemptySigma->GetYaxis()->SetTitleOffset(1.3);
	hemptySigma->GetXaxis()->SetTitleSize(0.045);
	hemptySigma->GetYaxis()->SetTitleSize(0.045);
	hemptySigma->GetXaxis()->SetTitleFont(42);
	hemptySigma->GetYaxis()->SetTitleFont(42);
	hemptySigma->GetXaxis()->SetLabelFont(42);
	hemptySigma->GetYaxis()->SetLabelFont(42);
	hemptySigma->GetXaxis()->SetLabelSize(0.04);
	hemptySigma->GetYaxis()->SetLabelSize(0.04);  
	hemptySigma->SetMaximum(2);
	hemptySigma->SetMinimum(0.);
	hemptySigma->Draw();
	//gaeBplusReference->SetFillColor(5);//Color(2);
	gaeBplusReference->SetFillColor(kOrange);
	//gaeBplusReference->SetFillStyle(3001); 
	gaeBplusReference->SetFillStyle(1001); 
	gaeBplusReference->SetLineWidth(3);
	//gaeBplusReference->SetLineColor(5);//Color(2);
	gaeBplusReference->SetLineColor(kOrange);
	gaeBplusReference->Draw("5same");
	hPtSigma->SetLineColor(1);
	hPtSigma->SetLineWidth(2);
	hPtSigma->SetMarkerStyle(20);
	hPtSigma->SetMarkerSize(1.2);
	hPtSigma->Draw("epsame"); 
	gaeCrossSyst->SetFillColor(1);
	gaeCrossSyst->SetFillStyle(0); 
	gaeCrossSyst->SetLineWidth(2);
	gaeCrossSyst->SetLineColor(1);
	gaeCrossSyst->Draw("5same");  

	//TLatex* texCms = new TLatex(0.16,0.95, "#scale[1.25]{CMS}");
	TLatex* texCms = new TLatex(0.18,0.87, "CMS");
	texCms->SetNDC();
	texCms->SetTextAlign(13);
	texCms->SetTextSize(0.07);
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
	if (label=="PbPb") { text="350.68 #mub^{-1} (PbPb 5.02 TeV)";}
	else {text="28.0 pb^{-1} (pp 5.02 TeV)";}
	TLatex* texlumi = new TLatex(0.90,0.936,text.Data());
	texlumi->SetNDC();
	texlumi->SetTextAlign(31);
	texlumi->SetTextFont(42);
	texlumi->SetTextSize(0.050);
	texlumi->SetLineWidth(2);
	texlumi->Draw();

	TLatex* texCol = new TLatex(0.94,0.95, Form("%s #sqrt{s_{NN}} = 5.02 TeV",label.Data()));
	texCol->SetNDC();
	texCol->SetTextAlign(32);
	texCol->SetTextSize(0.04);
	texCol->SetTextFont(42);
	//texCol->Draw();

	TString texper="%";
	TLatex* texCent = new TLatex(0.53,0.815,Form("Centrality %.0f - %.0f%s",centMin,centMax,texper.Data()));
	texCent->SetNDC();
	texCent->SetTextFont(42);
	texCent->SetTextSize(0.04);
	if(isPbPb) texCent->Draw();

	TLatex* texY = new TLatex(0.53,0.65,"|y| < 2.4");
	texY->SetNDC();
	texY->SetTextFont(42);
	texY->SetTextSize(0.05);
	texY->SetLineWidth(2);
	texY->Draw();

	TLatex* texB = new TLatex(0.82,0.82,"B^{+}");
	texB->SetNDC();
	texB->SetTextFont(42);
	texB->SetTextSize(0.07);
	texB->SetLineWidth(2);
	texB->Draw();

	//TLatex* texGlobal = new TLatex(0.53,0.59,Form("Global uncert. %.1f%s",normalizationUncertaintyForPP(),texper.Data()));
	TLatex* texGlobal = new TLatex(0.53,0.59,Form("L_{pp} + B uncert. %.1f%s",normalizationUncertaintyForPP(),texper.Data()));
	texGlobal->SetNDC();
	texGlobal->SetTextFont(42);
	texGlobal->SetTextSize(0.05);
	texGlobal->SetLineWidth(2);
	texGlobal->Draw();

	TLegend* leg_CS = new TLegend(0.52,0.70,0.85,0.80);
	leg_CS->SetBorderSize(0);
	leg_CS->SetFillStyle(0);
	leg_CS->SetTextSize(0.05);
	leg_CS->AddEntry(hPtSigma,"Data","pf");
	//leg_CS->AddEntry(gaeBplusReference,"FONLL pp ref.","f");//PAS
	leg_CS->AddEntry(gaeBplusReference,"FONLL","f");//paper
	leg_CS->Draw("same");

	/*
	   TLatex * tlatexlumi=new TLatex(0.671371,0.7801268,"L = 9.97 pb^{-1}");
	   tlatexlumi->SetNDC();
	   tlatexlumi->SetTextColor(1);
	   tlatexlumi->SetTextFont(42);
	   tlatexlumi->SetTextSize(0.045);
	   */

	cSigma->cd();
	TPad* pRatio = new TPad("pRatio","",0,0,1,0.25);
	pRatio->SetLeftMargin(0.1451613);
	pRatio->SetRightMargin(0.05443548);
	pRatio->SetTopMargin(0);
	pRatio->SetBottomMargin(0.30);//0.25
	pRatio->Draw();
	pRatio->cd();

	TH2F* hemptyRatio=new TH2F("hemptyRatio","",50,ptBins[0]-5.,ptBins[nBins]+5.,10.,0.,3.1);//PAS
	//TH2F* hemptyRatio=new TH2F("hemptyRatio","",50,ptBins[0]-5.,ptBins[nBins]+5.,10.,0.,2.1);//
	//TH2F* hemptyRatio=new TH2F("hemptyRatio","",50,ptBins[0]-5.,ptBins[nBins]+5.,10.,0.2,1.8);//paper
	hemptyRatio->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	hemptyRatio->GetXaxis()->CenterTitle();
	hemptyRatio->GetYaxis()->CenterTitle();
	hemptyRatio->GetYaxis()->SetTitle("Data / FONLL");
	hemptyRatio->GetXaxis()->SetTitleOffset(0.9);
	hemptyRatio->GetYaxis()->SetTitleOffset(0.5);
	hemptyRatio->GetXaxis()->SetTitleSize(0.12);
	hemptyRatio->GetYaxis()->SetTitleSize(0.12);
	hemptyRatio->GetXaxis()->SetTitleFont(42);
	hemptyRatio->GetYaxis()->SetTitleFont(42);
	hemptyRatio->GetXaxis()->SetLabelFont(42);
	hemptyRatio->GetYaxis()->SetLabelFont(42);
	hemptyRatio->GetXaxis()->SetLabelSize(0.12);//0.1
	hemptyRatio->GetYaxis()->SetLabelSize(0.12);//0.1  
	hemptyRatio->GetYaxis()->SetNdivisions(505);
	hemptyRatio->Draw();

	TLine* l = new TLine(ptBins[0]-5.,1,ptBins[nBins]+5.,1);
	l->SetLineWidth(1);
	l->SetLineStyle(2);
	gaeRatioCrossFONLLunity->Draw("5same");
	gaeRatioCrossFONLLstat->Draw("epsame");
	gaeRatioCrossFONLLsyst->Draw("5same");
	l->Draw("same");

	TString _postfix = "";
	if(doDataCor==1) _postfix = "_EFFCOR";
	if(!isPbPb) cSigma->SaveAs(Form("plotCrossSection/canvasSigmaBplusRatio%s%s.pdf",label.Data(),_postfix.Data()));
	else cSigma->SaveAs(Form("plotCrossSection/canvasSigmaBplusRatio%s_%.0f_%.0f%s.pdf",label.Data(),centMin,centMax,_postfix.Data()));
	if(!isPbPb) cSigma->SaveAs(Form("plotCrossSection/canvasSigmaBplusRatio%s%s.png",label.Data(),_postfix.Data()));
	else cSigma->SaveAs(Form("plotCrossSection/canvasSigmaBplusRatio%s_%.0f_%.0f%s.pgn",label.Data(),centMin,centMax,_postfix.Data()));


	TCanvas* cEff = new TCanvas("cEff","",550,500);

	TH2F* hemptyEff=new TH2F("hemptyEff","",50,0.,110.,10.,0,1.);  
	hemptyEff->GetXaxis()->CenterTitle();
	hemptyEff->GetYaxis()->CenterTitle();
	hemptyEff->GetYaxis()->SetTitle("#alpha x #epsilon_{reco} x #epsilon_{sel} ");
	hemptyEff->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	hemptyEff->GetXaxis()->SetTitleOffset(0.9);
	hemptyEff->GetYaxis()->SetTitleOffset(1.05);
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
	hEff->SetLineWidth(2);
	hEff->SetLineColor(1);
	hEff->SetMarkerStyle(20);
	hEff->SetMarkerSize(1.2);
	hEff->Draw("same");


	//  TString text;
	TString sample;
	/*
	   if (label=="PbPb") { text="CMS Preliminary     PbPb #sqrt{s}= 5.02 TeV"; sample="Pythia8+Hydjet MC simulation, prompt D^{0}";}
	   else {text="CMS Preliminary     pp #sqrt{s}= 5.02 TeV"; sample="Pythia8 MC simulation, prompt D^{0}";}

	   TLatex * tlatexeff=new TLatex(0.1612903,0.8525793,text.Data());
	   tlatexeff->SetNDC();
	   tlatexeff->SetTextColor(1);
	   tlatexeff->SetTextFont(42);
	   tlatexeff->SetTextSize(0.04);
	   tlatexeff->Draw();
	   */

	TLatex * tlatexeff2=new TLatex(0.1612903,0.7925793,sample.Data());
	tlatexeff2->SetNDC();
	tlatexeff2->SetTextColor(1);
	tlatexeff2->SetTextFont(42);
	tlatexeff2->SetTextSize(0.04);
	tlatexeff2->Draw();

	if(!isPbPb) cEff->SaveAs(Form("plotOthers/efficiency%s%s.pdf",label.Data(),_postfix.Data()));
	else cEff->SaveAs(Form("plotOthers/efficiency%s_%.0f_%.0f%s.pdf",label.Data(),centMin,centMax,_postfix.Data()));
	if(!isPbPb) cEff->SaveAs(Form("plotOthers/efficiency%s%s.png",label.Data(),_postfix.Data()));
	else cEff->SaveAs(Form("plotOthers/efficiency%s_%.0f_%.0f%s.png",label.Data(),centMin,centMax,_postfix.Data()));


	TFile *outputfile=new TFile(outputplot.Data(),"recreate");
	outputfile->cd();
	gaeCrossSyst->Write();
	gaeRatioCrossFONLLstat->Write();
	gaeBplusReference->Write();
	hPtSigma->Write();
	gaeRatioCrossFONLLstat->Write();
	gaeRatioCrossFONLLsyst->Write();
	gaeRatioCrossFONLLunity->Write();
	hEff->Write();
}


int main(int argc, char *argv[])
{
	if(argc==11)
	{
		CrossSectionRatio(argv[1], argv[2], argv[3],argv[4],atoi(argv[5]),argv[6],atoi(argv[7]),atof(argv[8]),atof(argv[9]),atof(argv[10]));
		return 0;
	}
	else if(argc==9)
	{
		CrossSectionRatio(argv[1], argv[2], argv[3],argv[4],atoi(argv[5]),argv[6],atoi(argv[7]),atof(argv[8]));
		return 0;
	}
	else
	{
		std::cout << "Wrong number of inputs" << std::endl;
		return 1;
	}

}


