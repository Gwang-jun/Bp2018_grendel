#include "TROOT.h"
#include "TH1.h"
#include "TTree.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TSystem.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom.h"
#include <iostream>
#include <fstream>
#include "uti.h"
#include "parameters.h"

TString CommonCut = "Btrk1Pt>0.9 && Bpt>5.0 && (BsvpvDistance/BsvpvDisErr)>2.0 && Bchi2cl>0.05 && TMath::Abs(Btrk1Eta)<2.4 && TMath::Abs(By)<2.4 && TMath::Abs(PVz)<15 && Bmass>5 && Bmass<6 && TMath::Abs(Bmumumass-3.096900)<0.15 && Bmu1SoftMuID && Bmu2SoftMuID && Bmu1isAcc && Bmu2isAcc && Bmu1isTriggered && Bmu2isTriggered && (Btrk1PixelHit+Btrk1StripHit)>=11 && (Btrk1Chi2ndf/(Btrk1nStripLayer+Btrk1nPixelLayer))<0.18 && TMath::Abs(Btrk1PtErr/Btrk1Pt)<0.1 && hiBin>=0 && hiBin <=180";
TString SignalCut = "$CommonCut && (abs(Bmass-5.27932)<0.08)";
TString SideBandCut = "$CommonCut && (abs(Bmass-5.27932)>0.15 && abs(Bmass-5.27932)<0.25)";
TString GenCut = "$CommonCut && (abs(Bmass-5.27932)<0.08 && Bgen==23333)";
TString weight="pthatweight*Ncoll*(1.034350*TMath::Exp(-0.000844*(PVz+3.502992)*(PVz+3.502992)))*(0.715021+0.039896*Bgenpt-0.000834*Bgenpt*Bgenpt+0.000006*Bgenpt*Bgenpt*Bgenpt)";


//#include "his.h"                                                                                                                             
 
using namespace std;

using std::cout;
using std::endl;

void ReweightBpt(){

  gStyle->SetOptStat(0);

  TFile *finData = new TFile("/afs/cern.ch/work/g/gwkim/private/samples/crab_Bfinder_20190513_HIDoubleMuon__PsiPeri__HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_skimhltBsize_ntKp.root");
  TFile *finMC = new TFile("/afs/cern.ch/work/g/gwkim/private/samples/crab_Bfinder_20190624_Hydjet_Pythia8_Official_BuToJpsiK_1033p1_pt3tkpt0p7dls2_allpthat_pthatweight.root");

  TTree * tData = (TTree*) finData->Get("Bfinder/ntKp");
  tData->AddFriend("hltanalysis/HltTree");
  tData->AddFriend("hiEvtAnalyzer/HiTree");
  TTree * tMC = (TTree*) finMC->Get("Bfinder/ntKp");
  tMC->AddFriend("hltanalysis/HltTree");
  tMC->AddFriend("hiEvtAnalyzer/HiTree");
  tMC->AddFriend("Bfinder/ntGen");
  TH1D * BptData = new TH1D("BptData","BptData",nBins,ptBins);
  TH1D * BptMC = new TH1D("BptMC","BptMC",nBins,ptBins);

  tData->Project("BptData","Bpt");
  BptData->GetXaxis()->SetTitle("Bpt (GeV)");
  BptData->GetYaxis()->SetTitle("Normalized Counts");
  BptData->SetTitle("Bpt Distribution in B^{+} Data");
  BptData->SetMarkerStyle(20);
  BptData->SetMarkerColor(kBlack);
  BptData->SetMarkerSize(1);
  BptData->Sumw2();
  divideBinWidth(BptData);
  BptData->Scale(1.0/BptData->Integral("width"));

  tMC->Project("BptMC","Bpt",TCut(weight)*"(Bgen==23333)");
  BptMC->GetXaxis()->SetTitle("Bpt (GeV)");
  BptMC->GetYaxis()->SetTitle("Normalized Counts");
  BptMC->SetTitle("Bpt Distribution in B^{+} MC");
  BptMC->SetMarkerStyle(20);
  BptMC->SetMarkerColor(kBlack);
  BptMC->SetMarkerSize(1);
  BptMC->Sumw2();
  divideBinWidth(BptMC);
  BptMC->Scale(1.0/BptMC->Integral("width"));

  TH1D * Ratio = (TH1D*) BptData->Clone("Ratio");
  Ratio->Sumw2();
  Ratio->Divide(BptMC);
  Ratio->GetXaxis()->SetTitle("Bpt (GeV)");
  Ratio->GetYaxis()->SetTitle("Ratio of the Normalized Counts");
  Ratio->SetTitle("Ratio of the Bpt Distributions in B^{+}");
  Ratio->SetMarkerStyle(20);
  Ratio->SetMarkerColor(kBlack);
  Ratio->SetMarkerSize(1);
  Ratio->Scale(1.0/Ratio->Integral());

  TString FuncData;
  TString FuncMC;
  TString FuncRatio;

  TF1 * f1 = new TF1("f1","[0]+[1]*x+[2]*x*x+[3]*x*x*x",nBins,ptBins);
  f1->SetParLimits(0,0,10);
  f1->SetParLimits(1,-10,10);
  f1->SetParLimits(2,-10,10);
  f1->SetParLimits(3,-10,10);

  TCanvas* c = new TCanvas("c","",600,600);
  c->cd();

  BptData->Fit(f1,"R");
  double p0Data = f1->GetParameter(0);
  double p1Data = f1->GetParameter(1);
  double p2Data = f1->GetParameter(2);
  FuncData = Form("%f*TMath::Exp(-%f*(Bpt-%f)*(Bpt-%f))",f1->GetParameter(0),f1->GetParameter(1),f1->GetParameter(2),f1->GetParameter(2));
  cout << "Funct Data = " << FuncData.Data() << endl;
  BptData->Draw("ep");
  c->SaveAs("BptDataFit.png");

  BptMC->Fit(f1,"R");
  double p0MC = f1->GetParameter(0);
  double p1MC = f1->GetParameter(1);
  double p2MC = f1->GetParameter(2);
  FuncMC = Form("%f*TMath::Exp(-%f*(Bpt-%f)*(Bpt-%f))",f1->GetParameter(0),f1->GetParameter(1),f1->GetParameter(2),f1->GetParameter(2));
  cout << "Funct MC = " << FuncMC.Data() << endl;
  BptMC->Draw("ep");
  c->SaveAs("BptMCFit.png");

  Ratio->Fit(f1,"R");
  double p0Ratio = f1->GetParameter(0);
  double p1Ratio = f1->GetParameter(1);
  double p2Ratio = f1->GetParameter(2);
  FuncRatio = Form("%f*TMath::Exp(-%f*(Bpt-%f)*(Bpt-%f))",f1->GetParameter(0),f1->GetParameter(1),f1->GetParameter(2),f1->GetParameter(2));
  cout << "Funct Ratio = " << FuncRatio.Data() << endl;
  Ratio->Draw("ep");
  c->SaveAs("BptRatioFit.png");
}

int main()
{
  ReweightBpt();
  return 0;
}
