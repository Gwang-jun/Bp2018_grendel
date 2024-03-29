#include "uti.h"
#include "parameters.h"
#include "TF1.h"
#include <TFitResultPtr.h>

int plot2D=0;
int plotClosure=1;
TString closureplotname = "plotClosure/Closure_Bplusbin";

Double_t setparam0=100.;
Double_t setparam1=5.28;
Double_t setparam2=0.05;
Double_t setparam3=0.03;
Double_t fixparam1=5.279;

Double_t minhisto=5.0;
Double_t maxhisto=6.0;
Double_t nbinsmasshisto=50;
Double_t binwidthmass=(maxhisto-minhisto)/nbinsmasshisto;

TString weight;
TString weightgen;
TString weightdata;
TString weightdataBgenpt;
TString seldata;
TString selmc;
TString selmceff;
TString selmcgen;
TString collisionsystem;
Float_t hiBinMin,hiBinMax,centMin,centMax;
double _ErrCor=1;

int _nBins = nBins;
double *_ptBins = ptBins;

Double_t yield;
Double_t yieldErr;

void fitB(int usePbPb=0, TString inputdata="" , TString inputmc="", TString trgselection="1",  TString cut="", TString cutmcgen="", int isMC=0, Double_t luminosity=1., int doweight=0, TString collsyst="PbPb", TString outputfile="", TString npfit="0", int doDataCor = 0, Float_t centmin=0., Float_t centmax=90.)
{
  collisionsystem=collsyst;
  if(collisionsystem=="ppInc"||collisionsystem=="PbPbInc"){
    _nBins = nBinsInc;
    _ptBins = ptBinsInc;
  }
  
  hiBinMin = centmin*2;
  hiBinMax = centmax*2;
  centMin = centmin;
  centMax = centmax;

  double ErrorOnSigma(double width, double errwidth, double smear, double errsmearing);

  if (!(usePbPb==1||usePbPb==0)) std::cout<<"ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!, you are using a non valid isPbPb option"<<std::endl;
  bool isPbPb=(bool)(usePbPb);

  if(!isPbPb)
    {
      seldata = Form("%s&&%s",trgselection.Data(),cut.Data());
      selmceff = Form("%s&&%s",trgselection.Data(),cut.Data());
      selmcgen = Form("%s",cutmcgen.Data());
      selmc = Form("%s&&%s",trgselection.Data(),cut.Data());
    }
  else
    {
      seldata = Form("%s&&%s&&hiBin>=%f&&hiBin<=%f",trgselection.Data(),cut.Data(),hiBinMin,hiBinMax);
      selmceff = Form("%s&&%s&&hiBin>=%f&&hiBin<=%f",trgselection.Data(),cut.Data(),hiBinMin,hiBinMax);
      selmcgen = Form("%s&&hiBin>=%f&&hiBin<=%f",cutmcgen.Data(),hiBinMin,hiBinMax);
      selmc = Form("%s&&%s&&hiBin>=%f&&hiBin<=%f",trgselection.Data(),cut.Data(),hiBinMin,hiBinMax);
    }

  if(plot2D==1)
    {
      gStyle->SetPalette(55,0);
      gStyle->SetOptStat(0);
    }
  if(plotClosure==1)
    {
      gStyle->SetOptStat(0);
      gStyle->SetTextSize(0.05);
      gStyle->SetTextFont(42);
      gStyle->SetPadRightMargin(0.043);
      gStyle->SetPadLeftMargin(0.18);
      gStyle->SetPadTopMargin(0.1);
      gStyle->SetPadBottomMargin(0.145);
      gStyle->SetTitleX(.0f);
    }
  if(plot2D==0 && plotClosure==0)
    {
      gStyle->SetTextSize(0.05);
      gStyle->SetTextFont(42);
      gStyle->SetPadRightMargin(0.043);
      gStyle->SetPadLeftMargin(0.18);
      gStyle->SetPadTopMargin(0.1);
      gStyle->SetPadBottomMargin(0.145);
      gStyle->SetTitleX(.0f);
    }

void clean0 (TH1D* h);
void getNPFnPar(TString npfname, float par[]);
//TF1* fit (TTree* nt, TTree* ntMC, double ptmin, double ptmax, int isMC,bool, TF1* &total,Float_t centmin, Float_t centmax, float NPpar[]);
TF1* fit (TTree* nt, TTree* ntMC, double ptmin, double ptmax, int isMC,bool, TF1* &total,Float_t centmin, Float_t centmax, TString npfit);
//float NPpar[2];
//getNPFnPar(npfile, NPpar);
//std::cout<<"NP parameter 0: "<<NPpar[0]<<std::endl;
//std::cout<<"NP parameter 1: "<<NPpar[1]<<std::endl;
 
//weightdata="TMath::Exp(3.11695e-08-5.16020e-02*Bpt+2.69860e-03*Bpt*Bpt-3.06583e-05*Bpt*Bpt*Bpt+4.61374e+01/Bpt)";
//weightdata="(3.91249e+00+TMath::Exp(-7.63499e-01*(Bpt-1.55399e+01))+TMath::Exp(-1.54551e-01*(Bpt-2.98508e+01)))";//Cent0-90%
//weightdata="(3.59596e+00+TMath::Exp(-7.33720e-02*(Bpt-3.70229e+01))+TMath::Exp(-9.03380e-01*(Bpt-1.44216e+01))+TMath::Exp(-2.85942e-01*(Bpt-2.22367e+01)))";//Cent0-30%
//weightdata="(2.66429e+00+TMath::Exp(-9.09252e-01*(Bpt-1.30721e+01))+TMath::Exp(-9.62700e-02*(Bpt-2.85406e+01))+TMath::Exp(-3.34052e-01*(Bpt-1.93575e+01)))";//Cent30-90%
weightdata="1";
weightdataBgenpt="1";
 if(!isPbPb)
   {
     weightgen="pthatweight*(0.599212+-0.020703*Gpt+0.003143*Gpt*Gpt+-0.000034*Gpt*Gpt*Gpt)";
     weight="pthatweight*(0.599212+-0.020703*Bgenpt+0.003143*Bgenpt*Bgenpt+-0.000034*Bgenpt*Bgenpt*Bgenpt)";
   }
 else
   {
     //weightgen="pthatweight*((3.506006+0.963473*Gpt+-0.258731*Gpt*Gpt)*TMath::Exp(-0.386065*Gpt)+1.139897)";
     //weight="pthatweight*Ncoll*(TMath::Gaus(PVz,0.427450,4.873825)/(sqrt(2*3.14159)*4.873825))/(TMath::Gaus(PVz,0.909938,4.970989)/(sqrt(2*3.14159)*4.970989))*((3.506006+0.963473*Bgenpt+-0.258731*Bgenpt*Bgenpt)*TMath::Exp(-0.386065*Bgenpt)+1.139897)";
     weightgen="pthatweight*(3.00448277-0.35865276*Gpt+0.01997413*Gpt*Gpt-0.00042585*Gpt*Gpt*Gpt+0.00000315*Gpt*Gpt*Gpt*Gpt)";
     //weight="pthatweight*Ncoll*(TMath::Gaus(PVz,0.427450,4.873825)/(sqrt(2*3.14159)*4.873825))/(TMath::Gaus(PVz,0.909938,4.970989)/(sqrt(2*3.14159)*4.970989))*(3.00448277-0.35865276*Bgenpt+0.01997413*Bgenpt*Bgenpt-0.00042585*Bgenpt*Bgenpt*Bgenpt+0.00000315*Bgenpt*Bgenpt*Bgenpt*Bgenpt)";
     //weight="pthatweight*Ncoll*(TMath::Gaus(PVz,0.427450,4.873825)/(sqrt(2*3.14159)*4.873825))/(TMath::Gaus(PVz,0.909938,4.970989)/(sqrt(2*3.14159)*4.970989))";
     weight="1";
   }

std::cout<<"we are using weight="<<weight<<std::endl;
std::cout<<"we are using weightdata="<<weightdata<<std::endl;
std::cout<<"we are using centrality="<<centmin<<"-"<<centmax<<"%"<<std::endl;

TFile* inf = new TFile(inputdata.Data());
TFile* infMC = new TFile(inputmc.Data());


//For 2018 PbPb data, MC
TTree* nt = (TTree*)inf->Get("Bfinder/ntKp");
nt->AddFriend("hltanalysis/HltTree");
nt->AddFriend("hiEvtAnalyzer/HiTree");
nt->AddFriend("skimanalysis/HltTree");
nt->AddFriend("BDT");
TTree* ntGen = (TTree*)infMC->Get("Bfinder/ntGen");
ntGen->AddFriend("hltanalysis/HltTree");
ntGen->AddFriend("hiEvtAnalyzer/HiTree");
ntGen->AddFriend("Bfinder/ntKp"); //call PVz
ntGen->AddFriend("skimanalysis/HltTree");
ntGen->AddFriend("BDT");
TTree* ntMC = (TTree*)infMC->Get("Bfinder/ntKp");
ntMC->AddFriend("hltanalysis/HltTree");
ntMC->AddFriend("hiEvtAnalyzer/HiTree");
ntMC->AddFriend("Bfinder/ntGen"); //call Bgen
ntMC->AddFriend("skimanalysis/HltTree");
ntMC->AddFriend("BDT");

/*
//For 2015 PbPb, pp data, MC
TTree* nt = (TTree*)inf->Get("ntKp");
nt->AddFriend("ntHlt");
nt->AddFriend("ntHi");
nt->AddFriend("ntSkim");
nt->AddFriend("mvaTree");
TTree* ntGen = (TTree*)infMC->Get("ntGen");
ntGen->AddFriend("ntHlt");
ntGen->AddFriend("ntHi");
TTree* ntMC = (TTree*)infMC->Get("ntKp");
ntMC->AddFriend("ntHlt");
ntMC->AddFriend("ntHi");
ntMC->AddFriend("ntSkim");
ntMC->AddFriend("ntGen");
ntMC->AddFriend("mvaTree");
*/
 
TH1D* hMean = new TH1D("hMean","",_nBins,_ptBins);                       
TH1D* hSigmaGaus1 = new TH1D("hSigmaGaus1","",_nBins,_ptBins); 
TH1D* hSigmaGaus2 = new TH1D("hSigmaGaus2","",_nBins,_ptBins); 
TF1 * totalmass;

 if(plot2D==1)
   {

     double LowBinWidth = 0.5;
     int NLowBin = 5/LowBinWidth;
     double HighBinWidth = 1;
     int NHighBin = 50/HighBinWidth;
     const int BptBin = NHighBin + NLowBin;
     double BptBinning[BptBin + 1];
     for(int i = 0; i < NLowBin; i++){
       BptBinning[i] = 5 + i * LowBinWidth;
     }
     for(int i = 0; i <  NHighBin+1; i++){
       BptBinning[i+NLowBin] = 10 + i * HighBinWidth;
     }

     TH2D* h2D = new TH2D("h2D","",BptBin,BptBinning,nBinsY,ptBinsY);
     nt->Project("h2D","abs(By):Bpt",Form("%s&&Bpt>%f&&Bpt<%f&&((Bpt>5&&Bpt<10&&TMath::Abs(By)>1.5)||(Bpt>10))",seldata.Data(),ptBins[0],ptBins[nBins]));
     h2D->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
     h2D->GetYaxis()->SetTitle("B^{+} |y|");
     h2D->GetYaxis()->SetTitleOffset(1.5);
     h2D->GetXaxis()->CenterTitle();
     h2D->GetYaxis()->CenterTitle();
     TCanvas* c2D = new TCanvas("","",600,600);
     c2D->cd();
     h2D->Draw("COLZ");
     c2D->SaveAs(Form("plotAverageEff/hMass2D_data_hyperfine_Cent%.0f-%.0f.png",centmin,centmax));
     c2D->SaveAs(Form("plotAverageEff/hMass2D_data_hyperfine_Cent%.0f-%.0f.pdf",centmin,centmax));

     /*
     TH2D* h2D = new TH2D("h2D","",45,5,50,60,0,180);
     nt->Project("h2D","hiBin:Bpt",Form("(%s&&Bpt>%f&&Bpt<%f)",seldata.Data(),ptBins[0],ptBins[nBins]));
     h2D->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
     h2D->GetYaxis()->SetTitle("B^{+} hiBin (Cent.*2)");
     h2D->GetYaxis()->SetTitleOffset(1.5);
     h2D->GetXaxis()->CenterTitle();
     h2D->GetYaxis()->CenterTitle();
     TCanvas* c2D = new TCanvas("","",600,600);
     c2D->cd();
     h2D->Draw("COLZ");
     c2D->SaveAs(Form("plotAverageEff/hMass2D_hiBin_Cent%.0f-%.0f.png",centmin,centmax));
     c2D->SaveAs(Form("plotAverageEff/hMass2D_hiBin_Cent%.0f-%.0f.pdf",centmin,centmax));
     */
   }
 
TString outputf;
outputf = Form("%s",outputfile.Data());
TFile* outf = new TFile(outputf.Data(),"recreate");
 outf->cd();

TH1D* hPt = new TH1D("hPt","",_nBins,_ptBins);
TH1D* hPtRecoTruth = new TH1D("hPtRecoTruth","",_nBins,_ptBins);
TH1D* hPtMC = new TH1D("hPtMC","",_nBins,_ptBins);
TH1D* hPtGen = new TH1D("hPtGen","",_nBins,_ptBins);

 for(int i=0;i<_nBins;i++)
   {
     //TF1* f = fit(nt,ntMC,_ptBins[i],_ptBins[i+1],isMC,isPbPb, totalmass,centmin, centmax, NPpar);     
     //if(i==(_nBins-1)){nbinsmasshisto=25;binwidthmass=(maxhisto-minhisto)/nbinsmasshisto;}
     
     TF1* f = fit(nt,ntMC,_ptBins[i],_ptBins[i+1],isMC,isPbPb,totalmass,centmin,centmax,npfit);
     hMean->SetBinContent(i+1,f->GetParameter(1));
     hMean->SetBinError(i+1,f->GetParError(1));     
     yieldErr = yieldErr*_ErrCor;
     hPt->SetBinContent(i+1,yield/(_ptBins[i+1]-_ptBins[i]));
     hPt->SetBinError(i+1,yieldErr/(_ptBins[i+1]-_ptBins[i]));
   }  
 
 ntMC->Project("hPtMC","Bpt",TCut(weight)*(TCut(selmceff.Data())&&"(Bgen==23333)"));
 divideBinWidth(hPtMC);
 ntMC->Project("hPtRecoTruth","Bpt",TCut(selmceff.Data())&&"(Bgen==23333)");
 divideBinWidth(hPtRecoTruth);
 ntGen->Project("hPtGen","Gpt",TCut(weightgen)*(TCut(selmcgen.Data())));
 divideBinWidth(hPtGen);

 TCanvas* cPt =  new TCanvas("cPt","",600,600);
 cPt->cd();
 cPt->SetLogy();
 hPt->SetXTitle("B^{+} p_{T} (GeV/c)");
 hPt->SetYTitle("Uncorrected dN(B^{+})/dp_{T}");
 hPt->Sumw2();
 hPt->Draw();
 if(isMC==1)
   {
     hPtMC->Draw("same hist");
     TLegend* legPt = myLegend(0.55,0.80,0.90,0.94);
     legPt->AddEntry(hPt,"Signal extraction","pl");
     legPt->AddEntry(hPtMC,"Matched reco","lf");
     legPt->Draw("same");  
   }
 hPtMC->Sumw2();
 
 TH1D* hEff = (TH1D*)hPtMC->Clone("hEff");
 hEff->SetTitle(";B^{+} p_{T} (GeV/c);Efficiency");
 hEff->Sumw2();
 hEff->Divide(hPtGen);
 
 TCanvas* cEff = new TCanvas("cEff","",600,600);
 hEff->Draw();
 
 TH1D* hPtCor = (TH1D*)hPt->Clone("hPtCor");
 hPtCor->SetTitle(";B^{+} p_{T} (GeV/c);Corrected dN(B^{+})/dp_{T}");
 hPtCor->Divide(hEff);
 TCanvas* cPtCor=  new TCanvas("cCorResult","",600,600);
 cPtCor->cd();
 cPtCor->SetLogy();
 hPtCor->Draw();
 if(isMC==1)
   {
     hPtGen->Draw("same hist");
     TLegend* legPtCor = myLegend(0.55,0.80,0.90,0.94);
     legPtCor->AddEntry(hPtCor,"Corrected signal","pl");
     legPtCor->AddEntry(hPtGen,"Generated B^{+}","lf");
     //legPtCor->Draw("same");
     hPtCor->SetLineColor(kRed);
     hPtGen->SetLineColor(kBlue);
     //hPtCor->Draw("same");
     //hPtGen->Draw("same");
   }
 
 TH1D* hPtSigma= (TH1D*)hPtCor->Clone("hPtSigma");
 hPtSigma->SetTitle(";B^{+} p_{T} (GeV/c);d#sigma(B^{+})/dp_{T} (pb/GeV)");
 hPtSigma->Scale(1./(2*luminosity*BRchain));
 
 TCanvas* cPtSigma= new TCanvas("cPtSigma","",600,600);
 cPtSigma->cd();
 cPtSigma->SetLogy();
 hPtSigma->Draw();

 TCanvas* cClosure = new TCanvas("cClosure","",600,600);
 cClosure->cd();
 TH1D* hClosure = (TH1D*)hPt->Clone("hClosure");
 hClosure->Divide(hPtMC);
 TH2F* hemptyClosure=new TH2F("hemptyClosure","",50,ptBins[0]-5,ptBins[nBins]+5,8,0.8,1.2);
 hemptyClosure->GetXaxis()->CenterTitle();
 hemptyClosure->GetYaxis()->CenterTitle();
 hemptyClosure->GetXaxis()->SetTitle("p_{T} (GeV/c)");
 hemptyClosure->GetYaxis()->SetTitle("MC corrected yield/Gen yield");
 hemptyClosure->GetXaxis()->SetTitleOffset(0.9);
 hemptyClosure->GetYaxis()->SetTitleOffset(1.3);
 hemptyClosure->GetXaxis()->SetTitleSize(0.04);
 hemptyClosure->GetYaxis()->SetTitleSize(0.04);
 hemptyClosure->GetXaxis()->SetTitleFont(42);
 hemptyClosure->GetYaxis()->SetTitleFont(42);
 hemptyClosure->GetXaxis()->SetLabelFont(42);
 hemptyClosure->GetYaxis()->SetLabelFont(42);
 hemptyClosure->GetXaxis()->SetLabelSize(0.035);
 hemptyClosure->GetYaxis()->SetLabelSize(0.035);
 hemptyClosure->Draw();
 hClosure->Draw("same");

 cClosure->SaveAs(Form("%s_Cent%.0f-%.0f.png",closureplotname.Data(),centmin,centmax));
 cClosure->SaveAs(Form("%s_Cent%.0f-%.0f.pdf",closureplotname.Data(),centmin,centmax));

 for(int i=0;i<nBins;i++)
 {
   std::cout<<i<<"-th bin Closure deviation: "<<100*(hClosure->GetBinContent(i+1)-1.0)<<"%"<<std::endl;
 }

 hPt->Write();
 hEff->Write();
 hMean->Write();
 hPtGen->Write();
 hPtMC->Write();
 hPtCor->Write();
 hPtSigma->Write();
 hClosure->Write();
 outf->Close();
 
}

void clean0(TH1D* h)
{
  for (int i=1;i<=h->GetNbinsX();i++)
    {
      if(h->GetBinContent(i)==0) h->SetBinError(i,1);
    }
}

void getNPFnPar(TString npfname, float par[]){
  TFile* npf = new TFile(npfname.Data());
  TF1* f = (TF1*)npf->Get("f1");
  par[0] = f->GetParameter(1);
  par[1] = f->GetParameter(2);
}

//TF1 *fit(TTree *nt, TTree *ntMC, Double_t ptmin, Double_t ptmax, int isMC,bool isPbPb,TF1* &total,Float_t centmin, Float_t centmax, float NPpar[])
  TF1 *fit(TTree *nt, TTree *ntMC, Double_t ptmin, Double_t ptmax, int isMC,bool isPbPb,TF1* &total,Float_t centmin, Float_t centmax, TString npfit)
{
  static Int_t count=0;
  count++;
  TCanvas* c= new TCanvas(Form("c%d",count),"",600,600);
  //c->cd();

  /*
  TPad* pFit = new TPad("pFit","",0,0.3,1,1);
  pFit->SetBottomMargin(0);
  pFit->Draw();
  pFit->cd();
  */

  /*
  h->Draw("e");
  background->Draw("same");
  mass->Draw("same");
  f->Draw("same");
  texCms->Draw();
  texCol->Draw();
  texPt->Draw();
  texY->Draw();
  texYield->Draw();
  c->cd();
  */

  TH1D* h = new TH1D(Form("h-%d",count),"",nbinsmasshisto,minhisto,maxhisto);
  TH1D* hMCSignal = new TH1D(Form("hMCSignal-%d",count),"",nbinsmasshisto,minhisto,maxhisto);
  TString iNP = npfit;
  TF1* f = new TF1(Form("f%d",count),"[0]*([7]*TMath::Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*TMath::Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[5]*("+iNP+")");
  f->SetNpx(5000);
  f->SetLineWidth(5);
  
  if(isMC==1) ntMC->Project(Form("h-%d",count),"Bmass",Form("%s*(%s&&Bgen==23333&&Bpt>%f&&Bpt<%f)*(%s)",weight.Data(),seldata.Data(),ptmin,ptmax,weightdata.Data())); //Closure
  //if(isMC==1) ntMC->Project(Form("h-%d",count),"Bmass",Form("(%s&&Bpt>%f&&Bpt<%f)",seldata.Data(),ptmin,ptmax));
  else nt->Project(Form("h-%d",count),"Bmass",Form("(%s&&Bpt>%f&&Bpt<%f)*(%s)",seldata.Data(),ptmin,ptmax,weightdata.Data()));   
  ntMC->Project(Form("hMCSignal-%d",count),"Bmass",Form("(%s&&Bgen==23333&&Bpt>%f&&Bpt<%f)*(%s)",selmc.Data(),ptmin,ptmax,weightdataBgenpt.Data()));

  clean0(h);

  f->SetParLimits(0,0,1e4);//1e7
  f->SetParLimits(7,0,1);
  f->SetParLimits(1,5.25,5.30);  
  f->SetParLimits(2,0.01,0.05);
  f->SetParLimits(8,0.01,0.05);
  //f->SetParLimits(4,-1e5,1e5);
  f->SetParLimits(4,-5,5);
  //f->SetParLimits(5,0,1e4);
  f->SetParLimits(5,0,1);
  
  //Do the signal fit first
  
  f->SetParameter(0,setparam0);
  f->SetParameter(1,setparam1);
  f->SetParameter(2,setparam2);
  f->SetParameter(8,setparam3);
  //f->FixParameter(1,fixparam1);

  /*
  f->FixParameter(3,0);
  f->FixParameter(4,0);
  f->FixParameter(5,0);

  h->GetEntries();
  
  hMCSignal->Fit(Form("f%d",count),"q","",minhisto,maxhisto);
  hMCSignal->Fit(Form("f%d",count),"q","",minhisto,maxhisto);
  //f->ReleaseParameter(1);
  hMCSignal->Fit(Form("f%d",count),"L q","",minhisto,maxhisto);
  hMCSignal->Fit(Form("f%d",count),"L q","",minhisto,maxhisto);
  hMCSignal->Fit(Form("f%d",count),"L q","",minhisto,maxhisto);
  hMCSignal->Fit(Form("f%d",count),"L m","",minhisto,maxhisto);
  
  //Fix the signal fit and do the background fits
  
  f->FixParameter(1,f->GetParameter(1));
  f->FixParameter(2,f->GetParameter(2));
  f->FixParameter(7,f->GetParameter(7));
  f->FixParameter(8,f->GetParameter(8));
  
  if(isMC==0)
    {
      f->ReleaseParameter(3);
      f->ReleaseParameter(4);
      f->ReleaseParameter(5);
      f->SetParLimits(5,0,1000);
    }
  
  printf("Fixed para.:\n");
  printf("%f, %f, %f\n", f->GetParameter(2), f->GetParameter(7), f->GetParameter(8));
  */

  h->Fit(Form("f%d",count),"q","",minhisto,maxhisto);
  h->Fit(Form("f%d",count),"q","",minhisto,maxhisto);
  //f->ReleaseParameter(1);
  h->Fit(Form("f%d",count),"L q","",minhisto,maxhisto);
  h->Fit(Form("f%d",count),"L q","",minhisto,maxhisto);
  h->Fit(Form("f%d",count),"L q","",minhisto,maxhisto);
  h->Fit(Form("f%d",count),"L m S","",minhisto,maxhisto);

  TF1 *background = new TF1(Form("background%d",count),"[0]+[1]*x");
  background->SetParameter(0,f->GetParameter(3));
  background->SetParameter(1,f->GetParameter(4));
  background->SetLineColor(4);
  background->SetRange(minhisto,maxhisto);
  //background->SetLineStyle(2);//PAS
  background->SetLineStyle(7);//paper
  background->SetLineWidth(5);
  
  TF1 *Bkpi = new TF1(Form("fBkpi%d",count),"[0]*("+iNP+")");
  Bkpi->SetParameter(0,f->GetParameter(5));
  Bkpi->SetRange(minhisto,maxhisto);
  Bkpi->SetLineStyle(7);
  //Bkpi->SetFillStyle(3004);//PAS
  Bkpi->SetFillStyle(3005);//paper
  //Bkpi->SetLineColor(kGreen+1);//PAS
  //Bkpi->SetFillColor(kGreen+1);//PAS
  Bkpi->SetLineColor(kGreen+4);//paper
  Bkpi->SetFillColor(kGreen+4);//paper
  Bkpi->SetLineWidth(5);
  
  TF1 *mass = new TF1(Form("fmass%d",count),"[0]*([3]*TMath::Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[3])*TMath::Gaus(x,[1],[4])/(sqrt(2*3.14159)*[4]))");
  mass->SetParameters(f->GetParameter(0),f->GetParameter(1),f->GetParameter(2),f->GetParameter(7),f->GetParameter(8));
  mass->SetParError(0,f->GetParError(0));
  mass->SetParError(1,f->GetParError(1));
  mass->SetParError(2,f->GetParError(2));
  mass->SetParError(7,f->GetParError(7));
  mass->SetParError(8,f->GetParError(8));
  //mass->SetLineColor(2);//PAS
  //mass->SetFillColor(2);//PAS
  mass->SetFillColor(kOrange-3);//paper
  mass->SetLineColor(kOrange-3);//paper
  //mass->SetFillStyle(3004);//PAS
  mass->SetFillStyle(3002);//paper
  mass->SetLineWidth(5);
  //mass->SetLineStyle(2);//PAS
  mass->SetLineStyle(7);//paper
  
  //h->SetXTitle("m_{#mu#muK} (GeV/c^{2})");
  h->SetXTitle("m_{B} (GeV/c^{2})");
  
  double yunit = 1000./nbinsmasshisto;
  h->SetYTitle(Form("Events / (%.0f MeV/c^{2})",yunit));//basically, the unit is equivalent to binwidthmass in MeV
  
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  h->SetAxisRange(0,h->GetMaximum()*1.4*1.2,"Y");
  h->GetXaxis()->SetTitleOffset(0.9);
  h->GetYaxis()->SetTitleOffset(1.3);
  //h->GetXaxis()->SetLabelOffset(0.007);
  //h->GetYaxis()->SetLabelOffset(0.007);
  h->GetXaxis()->SetTitleSize(0.07);
  h->GetYaxis()->SetTitleSize(0.07);
  h->GetXaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleFont(42);
  h->GetXaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelSize(0.06);
  h->GetYaxis()->SetLabelSize(0.06);
  h->SetMarkerSize(1.55);
  h->SetMarkerStyle(20);
  h->SetLineColor(1);
  h->SetLineWidth(5);
  h->SetStats(0);
  h->GetXaxis()->SetNdivisions(-50205);
  h->Draw("e");
  Bkpi->SetRange(minhisto,maxhisto);
  if(isMC==0) {Bkpi->Draw("same"); background->Draw("same");}
  mass->SetRange(minhisto,maxhisto);
  mass->Draw("same");
  f->Draw("same");
  c->RedrawAxis();

  yield = mass->Integral(minhisto,maxhisto)/binwidthmass;
  yieldErr = mass->Integral(minhisto,maxhisto)/binwidthmass*mass->GetParError(0)/mass->GetParameter(0);
  printf("p_t bin %.0f-%.0f     yield: %f     yieldErr: %f\n", ptmin, ptmax, yield, yieldErr);
  
  Double_t Signal = mass->Integral(5.19932,5.35932)/binwidthmass; //B+ mass_pdg=5.27932GeV, signal region = pm 0.08GeV
  Double_t Bkg_comb = background->Integral(5.19932,5.35932)/binwidthmass;
  Double_t Bkg_nonprompt = Bkpi->Integral(5.19932,5.35932)/binwidthmass;
  
  printf("p_t bin %.0f-%.0f S: %f B(comb): %f B(comb+np): %f\n", ptmin, ptmax, Signal, Bkg_comb, Bkg_comb+Bkg_nonprompt);
  printf("p_t bin %.0f-%.0f sig(comb): %f sig(comb+np): %f\n", ptmin, ptmax, Signal/sqrt(Signal+Bkg_comb), Signal/sqrt(Signal+Bkg_comb+Bkg_nonprompt));

  TLegend *leg = new TLegend(0.55,0.55,0.875,0.775,NULL,"brNDC");//paper
  leg->SetBorderSize(0);
  leg->SetTextSize(0.045);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  if(isMC==0) leg->AddEntry(h,"Data","pl");
  if(isMC==1) leg->AddEntry(h,"MC","pl");
  leg->AddEntry(f,"Fit","l");
  leg->AddEntry(mass,"Signal","f");
  if(isMC==0) {leg->AddEntry(background,"Combinatorial","l"); leg->AddEntry(Bkpi,"B #rightarrow J/#psi X","f");}
  leg->Draw("same");
  
  TLatex* texYield = new TLatex(0.55,0.51,Form("Yield:%.2f#pm%.2f", yield, yieldErr));
  texYield->SetNDC();
  texYield->SetTextAlign(12);
  texYield->SetTextSize(0.035);
  texYield->SetTextFont(42);
  texYield->Draw();
  
  TLatex* texChi = new TLatex(0.55,0.475,Form("#chi^{2}/nDOF:%.2f/%d=%.2f", f->GetChisquare(), f->GetNDF(), f->GetChisquare()/f->GetNDF()));
  texChi->SetNDC();
  texChi->SetTextAlign(12);
  texChi->SetTextSize(0.035);
  texChi->SetTextFont(42);
  texChi->Draw();
  printf("NDF: %d, chi2: %f, prob: %f\n", f->GetNDF(), f->GetChisquare(), f->GetProb());
  
  TLatex* texcms = new TLatex(0.22,0.87,"CMS");
  texcms->SetNDC();
  texcms->SetTextAlign(13);
  texcms->SetTextFont(62);
  texcms->SetTextSize(0.08);
  texcms->SetLineWidth(2);
  texcms->Draw();
  
  TLatex* texB = new TLatex(0.22,0.73,"B^{+}+B^{-}");
  texB->SetNDC();
  texB->SetTextFont(42);
  texB->SetTextSize(0.07);
  texB->SetLineWidth(2);
  texB->Draw();
  
  TLatex* texCol;
  if(collisionsystem=="pp"||collisionsystem=="PP"||collisionsystem=="ppInc"||collisionsystem=="PbPbInc") texCol= new TLatex(0.945,0.94, Form("28.0 pb^{-1} (%s 5.02 TeV)","pp"));
  else texCol= new TLatex(0.94,0.93, Form("1.5 nb^{-1} (%s 5.02 TeV)","PbPb"));
  texCol->SetNDC();
  texCol->SetTextAlign(32);
  texCol->SetTextSize(0.06);
  texCol->SetTextFont(42);
  texCol->Draw();
  
  TLatex* tex;
  tex = new TLatex(0.49,0.845,Form("%.0f<p_{T}<%.0f GeV/c",ptmin,ptmax));
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.045);
  tex->SetLineWidth(2);
  tex->Draw();
  
  TString texper="%";
  if(isPbPb){
  tex = new TLatex(0.488,0.80,Form("Cent. %.0f-%.0f%s  |y|<2.4",centmin,centmax,texper.Data()));
  tex->SetNDC();
  tex->SetTextColor(1);
  tex->SetTextFont(42);
  tex->SetTextSize(0.045);
  tex->SetLineWidth(2);
  tex->Draw();
  }  

  /*
  TCanvas* c100 = new TCanvas("","",600,600);
  c100->cd();

  TH1D* hPull = new TH1D("","",nbinsmasshisto,minhisto,maxhisto);
  for(int i=0;i<nbinsmasshisto;i++)
    {
      hPull->SetBinContent(i+1,(h->GetBinContent(i+1)-f->Eval(h->GetBinCenter(i+1)))/h->GetBinError(i+1));
      hPull->SetBinError(i+1,0);
    }
  
  gStyle->SetOptStat(0);
  
  hPull->SetMinimum(-4.);
  hPull->SetMaximum(4.);
  hPull->SetYTitle("Pull");
  hPull->SetMarkerSize(1);
  hPull->GetXaxis()->SetTitleOffset(1.);
  hPull->GetYaxis()->SetTitleOffset(0.65);
  hPull->GetXaxis()->SetLabelOffset(0.007);
  hPull->GetYaxis()->SetLabelOffset(0.007);
  hPull->GetXaxis()->SetTitleSize(0.05);
  hPull->GetYaxis()->SetTitleSize(0.05);
  hPull->GetXaxis()->SetLabelSize(0.03);
  hPull->GetYaxis()->SetLabelSize(0.03);
  hPull->GetYaxis()->SetNdivisions(504);
  
  hPull->Draw("hist");
  
  if(isPbPb && isMC==0)
  {
    c100->SaveAs(Form("plotPull/Pull_AN_data_PbPb_pt%.0f-%.0f_cent%.0f-%.0f.png",ptmin,ptmax,centmin,centmax));
    c100->SaveAs(Form("plotPull/Pull_AN_data_PbPb_pt%.0f-%.0f_cent%.0f-%.0f.pdf",ptmin,ptmax,centmin,centmax));
  }
  if(isPbPb && isMC==1)
    {
      c100->SaveAs(Form("plotPull/Pull_AN_mc_PbPb_pt%.0f-%.0f_cent%.0f-%.0f.png",ptmin,ptmax,centmin,centmax));
      c100->SaveAs(Form("plotPull/Pull_AN_mc_PbPb_pt%.0f-%.0f_cent%.0f-%.0f.pdf",ptmin,ptmax,centmin,centmax));
    }
  if(!isPbPb && isMC==0)
    {
      c100->SaveAs(Form("plotPull/Pull_AN_data_pp_pt%.0f-%.0f.png",ptmin,ptmax));
      c100->SaveAs(Form("plotPull/Pull_AN_data_pp_pt%.0f-%.0f.pdf",ptmin,ptmax));
    }
  if(!isPbPb && isMC==1)
    {
      c100->SaveAs(Form("plotPull/Pull_AN_mc_pp_pt%.0f-%.0f.png",ptmin,ptmax));
      c100->SaveAs(Form("plotPull/Pull_AN_mc_pp_pt%.0f-%.0f.pdf",ptmin,ptmax));
    }

  c->cd();
  */

  /*
  TLine* lPull = new TLine(5.0,0,6.0,0);
  lPull->SetLineWidth(1);
  lPull->SetLineStyle(7);
  lPull->SetLineColor(1);

  TPad* pPull = new TPad("pPull","",0,0,1,0.3);
  pPull->SetTopMargin(0);
  pPull->SetBottomMargin(0.3);
  pPull->Draw();
  pPull->cd();
  hPull->Draw("p");
  lPull->Draw();
  c->cd();
  */

  total=f;
  
  TF1* t = (TF1*)h->GetFunction(Form("f%d",count))->Clone();
  h->GetFunction(Form("f%d",count))->Delete();
  t->Draw("same");
  h->Draw("e same");
  h->Write();
  hMCSignal->Write();
  
  TString postfix = "%";
  //if(weightdata!="1") postfix = "EFFCOR";
  if(isPbPb && isMC==0)
  {
    c->SaveAs(Form("plotFits/weighted_AN_data_PbPb_pt%.0f-%.0f_cent%.0f-%.0f.png",ptmin,ptmax,centmin,centmax));
    c->SaveAs(Form("plotFits/weighted_AN_data_PbPb_pt%.0f-%.0f_cent%.0f-%.0f.pdf",ptmin,ptmax,centmin,centmax));
  }
  if(isPbPb && isMC==1)
    {
      c->SaveAs(Form("plotFits/AN_mc_PbPb_pt%.0f-%.0f_cent%.0f-%.0f.png",ptmin,ptmax,centmin,centmax));
      c->SaveAs(Form("plotFits/AN_mc_PbPb_pt%.0f-%.0f_cent%.0f-%.0f.pdf",ptmin,ptmax,centmin,centmax));
    }
  if(!isPbPb && isMC==0)
    {
      c->SaveAs(Form("plotFits/AN_data_pp_pt%.0f-%.0f.png",ptmin,ptmax));
      c->SaveAs(Form("plotFits/AN_data_pp_pt%.0f-%.0f.pdf",ptmin,ptmax));
    }
  if(!isPbPb && isMC==1)
    {
      c->SaveAs(Form("plotFits/AN_mc_pp_pt%.0f-%.0f.png",ptmin,ptmax));
      c->SaveAs(Form("plotFits/AN_mc_pp_pt%.0f-%.0f.pdf",ptmin,ptmax));
    }
  return mass;
}

int main(int argc, char *argv[])
{
  if(argc==16)
    {
      fitB(atoi(argv[1]), argv[2], argv[3], argv[4], argv[5], argv[6], atoi(argv[7]), atof(argv[8]), atoi(argv[9]), argv[10], argv[11], argv[12], atoi(argv[13]), atof(argv[14]), atof(argv[15]));
      return 0;
    }
  else if(argc==14)
    {
      fitB(atoi(argv[1]), argv[2], argv[3], argv[4], argv[5], argv[6], atoi(argv[7]), atof(argv[8]), atoi(argv[9]), argv[10], argv[11], argv[12], atoi(argv[13]));
      return 0;
    }
  else
    {
      std::cout << "Wrong number of inputs" << std::endl;
      return 1;
    }
}

double ErrorOnSigma(double width, double errwidth, double smear, double errsmearing){
  double squarederroronsigma=(1+smear)*(1+smear)*errwidth*errwidth+width*width*errsmearing*errsmearing;
  double erroronsigma=TMath::Sqrt(squarederroronsigma);
  return erroronsigma;
}
