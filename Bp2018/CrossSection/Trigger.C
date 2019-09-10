#include "uti.h"
#include "parameters.h"
#include "TF1.h"

TString weight;
TString weightgen;
TString weightdata;
TString seldata;
TString selmc;
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

  bool isPbPb=(bool)(usePbPb);
  
  if(!isPbPb)
    {
      seldata = Form("%s&&%s",trgselection.Data(),cut.Data());
      selmc = Form("%s&&%s",trgselection.Data(),cut.Data());
    }
  else
    {
      seldata = Form("%s&&hiBin>=%f&&hiBin<=%f",cut.Data(),hiBinMin,hiBinMax);
      selmc = Form("%s&&hiBin>=%f&&hiBin<=%f",cut.Data(),hiBinMin,hiBinMax);
    }

  gStyle->SetOptStat(0);
  /*
  gStyle->SetTextSize(0.05);
  gStyle->SetTextFont(42);
  gStyle->SetPadRightMargin(0.043);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.145);
  gStyle->SetTitleX(.0f);
  */  

  weightdata="1";
  if(!isPbPb)
    {
      weightgen="pthatweight*(0.599212+-0.020703*Gpt+0.003143*Gpt*Gpt+-0.000034*Gpt*Gpt*Gpt)";
      weight="pthatweight*(0.599212+-0.020703*Bgenpt+0.003143*Bgenpt*Bgenpt+-0.000034*Bgenpt*Bgenpt*Bgenpt)";
    }
 else
   {
     weightgen="pthatweight*((3.506006+0.963473*Gpt+-0.258731*Gpt*Gpt)*TMath::Exp(-0.386065*Gpt)+1.139897)";
     weight="pthatweight*Ncoll*(TMath::Gaus(PVz,0.427450,4.873825)/(sqrt(2*3.14159)*4.873825))/(TMath::Gaus(PVz,0.909938,4.970989)/(sqrt(2*3.14159)*4.970989))*((3.506006+0.963473*Bgenpt+-0.258731*Bgenpt*Bgenpt)*TMath::Exp(-0.386065*Bgenpt)+1.139897)";
     //weightgen="pthatweight*((2.907795+-0.436572*Gpt+0.006372*Gpt*Gpt)*TMath::Exp(-0.157563*Gpt)+1.01308)";
     //weight="pthatweight*Ncoll*(TMath::Gaus(PVz,0.427450,4.873825)/(sqrt(2*3.14159)*4.873825))/(TMath::Gaus(PVz,0.909938,4.970989)/(sqrt(2*3.14159)*4.970989))*((2.907795+-0.436572*Bgenpt+0.006372*Bgenpt*Bgenpt)*TMath::Exp(-0.157563*Bgenpt)+1.01308)";
     //weightgen="pthatweight*(0.889175+0.000791*Gpt+0.000015*Gpt*Gpt)";
     //weight="pthatweight*Ncoll*(TMath::Gaus(PVz,0.427450,4.873825)/(sqrt(2*3.14159)*4.873825))/(TMath::Gaus(PVz,0.909938,4.970989)/(sqrt(2*3.14159)*4.970989))*(0.889175+0.000791*Bgenpt+0.000015*Bgenpt*Bgenpt)";
     //weightgen="pthatweight*(0.094376+0.028350*Gpt+-0.000225*Gpt*Gpt+5.369348/Gpt)";
     //weight="pthatweight*Ncoll*(TMath::Gaus(PVz,0.427450,4.873825)/(sqrt(2*3.14159)*4.873825))/(TMath::Gaus(PVz,0.909938,4.970989)/(sqrt(2*3.14159)*4.970989))*(0.094376+0.028350*Bgenpt+-0.000225*Bgenpt*Bgenpt+5.369348/Bgenpt)";
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
  
  TH1D* hTrg = new TH1D("hTrg","",201,0,200);
  
  int n, m;
  if(isMC!=0)
    {
      for(int i=0;i<=200;i++)
	{
	  n=ntMC->GetEntries(Form("%s&&hiBin<=%d",selmc.Data(),i));
	  m=ntMC->GetEntries(Form("%s&&%s&&hiBin<=%d",selmc.Data(),trgselection.Data(),i));	  
	  hTrg->SetBinContent(i,m/n);
	  hTrg->SetBinError(i,0);
	}
    }

  if(isMC==0)
    {
      for(int i=0;i<=200;i++)
	{
	  n=nt->GetEntries(Form("%s&&hiBin<=%d",seldata.Data(),i));
	  m=nt->GetEntries(Form("%s&&%s&&hiBin<=%d",seldata.Data(),trgselection.Data(),i));	  
	  hTrg->SetBinContent(i,m/n);
	  hTrg->SetBinError(i,0);        
	} 
    }
  
  TString outputf;
  outputf = Form("%s",outputfile.Data());
  TFile* outf = new TFile(outputf.Data(),"recreate");
  outf->cd();
      
  TCanvas* c =  new TCanvas("c","",600,600);
  hTrg->Draw("p");
  c->SaveAs("Trigger.png");

  hTrg->Write();
  outf->Close();
 
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
