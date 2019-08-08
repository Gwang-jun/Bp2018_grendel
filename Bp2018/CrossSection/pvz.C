#include "uti.h"

TString inputdataname = "/mnt/T2_US_MIT/submit-hi2/scratch/gwangjun/crab_Bfinder_20190513_HIDoubleMuonPsiPeri_HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_327123_327564_skimhltBsize_ntKp.root";
TString inputmcname = "/mnt/T2_US_MIT/submit-hi2/scratch/gwangjun/crab_Bfinder_20190520_Hydjet_Pythia8_BuToJpsiK_1033p1_pt3tkpt0p7dls2_v2_addSamplePthat_pthatweight.root";
TString CommonCut = "Btrk1Pt>0.9 && Bpt>5.0 && (BsvpvDistance/BsvpvDisErr)>2.0 && Bchi2cl>0.05 && TMath::Abs(Btrk1Eta)<2.4 && TMath::Abs(By)<2.4 && TMath::Abs(PVz)<15 && Bmass>5 && Bmass<6 && TMath::Abs(Bmumumass-3.096900)<0.15 && Bmu1SoftMuID && Bmu2SoftMuID && Bmu1isAcc && Bmu2isAcc && Bmu1isTriggered && Bmu2isTriggered && (Btrk1PixelHit+Btrk1StripHit)>=11 && (Btrk1Chi2ndf/(Btrk1nStripLayer+Btrk1nPixelLayer))<0.18 && TMath::Abs(Btrk1PtErr/Btrk1Pt)<0.1 && hiBin>=0 && hiBin <=180";
TString SignalCut = "$CommonCut && (abs(Bmass-5.27932)<0.08)";
TString SideBandCut = "$CommonCut && (abs(Bmass-5.27932)>0.15 && abs(Bmass-5.27932)<0.25)";
TString GenCut = "$CommonCut && (abs(Bmass-5.27932)<0.08 && Bgen==23333)";
TString weight = "pthatweight*Ncoll";

void pvz()
{
  TFile* inputdata = new TFile(inputdataname.Data());
  TFile* inputmc = new TFile(inputmcname.Data());
  
  TTree* data = (TTree*) inputdata->Get("Bfinder/ntKp");
  data->AddFriend("hltanalysis/HltTree");
  data->AddFriend("hiEvtAnalyzer/HiTree");
  TTree* mc = (TTree*) inputmc->Get("Bfinder/ntKp");
  mc->AddFriend("hltanalysis/HltTree");
  mc->AddFriend("hiEvtAnalyzer/HiTree");
  mc->AddFriend("Bfinder/ntGen");
  
  TH1F* pvz_data = new TH1F("pvz_data","",60,-30,30);
  TH1F* pvz_mc = new TH1F("pvz_mc","",60,-30,30);
  
  data->Project("pvz_data","PVz");
  mc->Project("pvz_mc","PVz - 0.321904");
  
  pvz_data->Scale(1.0/pvz_data->Integral());
  pvz_mc->Scale(1.0/pvz_mc->Integral());
  
  TF1* f_data = new TF1("f_data","[0]*TMath::Exp(-[1]*(x-[2])*(x-[2]))",-30,30);
  f_data->SetParameters(1,-1,0);
  f_data->SetParLimits(0,0,10);
  f_data->SetParLimits(1,0,10);
  f_data->SetParLimits(2,-10,10);
  
  pvz_data->Fit("f_data");
  
  TF1* f_mc = new TF1("f_mc","[0]*TMath::Exp(-[1]*(x-[2])*(x-[2]))",-30,30);
  f_mc->SetParameters(1,-1,0);
  f_mc->SetParLimits(0,0,10);
  f_mc->SetParLimits(1,0,10);
  f_mc->SetParLimits(2,-10,10);
  
  pvz_mc->Fit("f_mc");
  
  return;
  
}

int main()
{
  pvz();
  return 0;
}
