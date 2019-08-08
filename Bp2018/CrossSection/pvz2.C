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

const int NPVz = 60;
const double PVzMin = -30;
const double PVzMax = 30;

//#include "his.h"                                                                                                                             
 
using namespace std;

using std::cout;
using std::endl;


void ReweightPVz(){

  gStyle->SetOptStat(0);

  TFile *finData = new TFile("/afs/cern.ch/work/g/gwkim/private/samples/crab_Bfinder_20190513_HIDoubleMuon__PsiPeri__HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_skimhltBsize_ntKp.root");
  TFile *finMC = new TFile("/afs/cern.ch/work/g/gwkim/private/samples/crab_Bfinder_20190624_Hydjet_Pythia8_Official_BuToJpsiK_1033p1_pt3tkpt0p7dls2_allpthat_pthatweight.root");
  TFile *finMCprivate = new TFile("/afs/cern.ch/work/g/gwkim/private/samples/crab_Bfinder_20190520_Hydjet_Pythia8_BuToJpsiK_1033p1_pt3tkpt0p7dls2_v2_addSamplePthat_pthatweight.root");

  TTree * tData = (TTree*) finData->Get("Bfinder/ntKp");
  TTree * tMC = (TTree*) finMC->Get("Bfinder/ntKp");
  TTree * tMCprivate = (TTree*) finMCprivate->Get("Bfinder/ntKp");
  TH1D * PVzData = new TH1D("PVzData","PVzData",NPVz,PVzMin,PVzMax);
  TH1D * PVzMC = new TH1D("PVzMC","PVzMC",NPVz,PVzMin,PVzMax);
  TH1D * PVzMCprivate = new TH1D("PVzMCprivate","PVzMCprivate",NPVz,PVzMin,PVzMax);

  tData->Project("PVzData","PVz");
  PVzData->GetXaxis()->SetTitle("PVz (mm)");
  PVzData->GetYaxis()->SetTitle("Normalized Counts");
  PVzData->GetYaxis()->SetTitleOffset(1.1);
  PVzData->SetTitle("PVz Distribution of B^{+} Data");
  PVzData->SetMarkerStyle(20);
  PVzData->SetMarkerColor(kBlack);
  PVzData->SetMarkerSize(1);
  PVzData->Sumw2();
  PVzData->Scale(1.0/PVzData->Integral());

  //tMC->Project("PVzMC","PVz-0.321904");
  tMC->Project("PVzMC","PVz");
  PVzMC->GetXaxis()->SetTitle("PVz (mm)");
  PVzMC->GetYaxis()->SetTitle("Normalized Counts");
  PVzMC->GetYaxis()->SetTitleOffset(1.1);
  PVzMC->SetTitle("PVz Distribution of B^{+} MC");
  PVzMC->SetMarkerStyle(20);
  PVzMC->SetMarkerColor(kBlack);
  PVzMC->SetMarkerSize(1);
  PVzMC->Sumw2();
  PVzMC->Scale(1.0/PVzMC->Integral());

  tMCprivate->Project("PVzMCprivate","PVz");
  PVzMCprivate->GetXaxis()->SetTitle("PVz (mm)");
  PVzMCprivate->GetYaxis()->SetTitle("Normalized Counts");
  PVzMCprivate->GetYaxis()->SetTitleOffset(1.1);
  PVzMCprivate->SetTitle("PVz Distribution of B^{+} private MC");
  PVzMCprivate->SetMarkerStyle(20);
  PVzMCprivate->SetMarkerColor(kRed);
  PVzMCprivate->SetMarkerSize(1);
  PVzMCprivate->Sumw2();
  PVzMCprivate->Scale(1.0/PVzMCprivate->Integral());

  TCanvas* c0 = new TCanvas("c0","",600,600);
  c0->cd();
  PVzMC->Draw("ep");
  PVzMCprivate->Draw("ep same");
  c0->SaveAs("PVzOffPriComparison.png");

  TH1D * Ratio = (TH1D*) PVzData->Clone("Ratio");
  Ratio->Sumw2();
  Ratio->Divide(PVzMC);
  Ratio->GetXaxis()->SetTitle("PVz (mm)");
  Ratio->GetYaxis()->SetTitle("Ratio of the Normalized Counts");
  Ratio->GetYaxis()->SetTitleOffset(1.1);
  Ratio->SetTitle("Ratio of the PVz Distributions of B^{+}");
  Ratio->SetMarkerStyle(20);
  Ratio->SetMarkerColor(kBlack);
  Ratio->SetMarkerSize(1);
  Ratio->Scale(1.0/Ratio->Integral());

  TString FuncData;
  TString FuncMC;
  TString FuncRatio;

  TF1 * f1 = new TF1("f1","TMath::Gaus(x,[0],[1])/(sqrt(2*3.14159)*[1])",PVzMin,PVzMax);
  f1->SetParLimits(0,-15,15);
  f1->SetParLimits(1,3,30);

  TCanvas* c = new TCanvas("c","",600,600);
  c->cd();

  PVzData->Fit(f1,"R");
  double p0Data = f1->GetParameter(0);
  double p1Data = f1->GetParameter(1);
  FuncData = Form("TMath::Gaus(PVz,%f,%f)/(sqrt(2*3.14159)*%f)",f1->GetParameter(0),f1->GetParameter(1),f1->GetParameter(1));
  cout << "Funct Data = " << FuncData.Data() << endl;
  PVzData->Draw("ep");
  c->SaveAs("PVzDataFit.png");

  PVzMC->Fit(f1,"R");
  double p0MC = f1->GetParameter(0);
  double p1MC = f1->GetParameter(1);
  FuncMC = Form("TMath::Gaus(PVz,%f,%f)/(sqrt(2*3.14159)*%f)",f1->GetParameter(0),f1->GetParameter(1),f1->GetParameter(1));
  cout << "Funct MC = " << FuncMC.Data() << endl;
  PVzMC->Draw("ep");
  c->SaveAs("PVzMCFit.png");

  TF1 * f2 = new TF1("f2","[2]*TMath::Gaus(x,[0],[1])/(sqrt(2*3.14159)*[1])",PVzMin,PVzMax);
  f2->SetParLimits(0,-15,15);
  f2->SetParLimits(1,3,10);
  f2->SetParLimits(2,0,100);

  Ratio->Fit(f2,"R");
  double p0Ratio = f2->GetParameter(0);
  double p1Ratio = f2->GetParameter(1);
  double p2Ratio = f2->GetParameter(2);
  FuncRatio = Form("%f*TMath::Gaus(PVz,%f,%f)/(sqrt(2*3.14159)*%f)",f2->GetParameter(2),f2->GetParameter(0),f2->GetParameter(1),f2->GetParameter(1));
  cout << "Funct Ratio = " << FuncRatio.Data() << endl;
  Ratio->Draw("ep");
  c->SaveAs("PVzRatioFit.png");
}

int main()
{
  ReweightPVz();
  return 0;
}
