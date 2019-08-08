#include "uti.h"

const int nBins=7;
double ptBins[nBins+1] = {5., 7., 10., 15., 20., 30., 50., 100.};
TString inputdataname = "/mnt/T2_US_MIT/submit-hi2/scratch/gwangjun/crab_Bfinder_20190513_HIDoubleMuonPsiPeri_HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_327123_327564_skimhltBsize_ntKp.root";
TString inputmcname = "/mnt/T2_US_MIT/submit-hi2/scratch/gwangjun/crab_Bfinder_20190520_Hydjet_Pythia8_BuToJpsiK_1033p1_pt3tkpt0p7dls2_v2_addSamplePthat_pthatweight.root";
TString CommonCut = "Btrk1Pt>0.9 && Bpt>5.0 && (BsvpvDistance/BsvpvDisErr)>2.0 && Bchi2cl>0.05 && TMath::Abs(Btrk1Eta)<2.4 && TMath::Abs(By)<2.4 && TMath::Abs(PVz)<15 && Bmass>5 && Bmass<6 && TMath::Abs(Bmumumass-3.096900)<0.15 && Bmu1SoftMuID && Bmu2SoftMuID && Bmu1isAcc && Bmu2isAcc && Bmu1isTriggered && Bmu2isTriggered && (Btrk1PixelHit+Btrk1StripHit)>=11 && (Btrk1Chi2ndf/(Btrk1nStripLayer+Btrk1nPixelLayer))<0.18 && TMath::Abs(Btrk1PtErr/Btrk1Pt)<0.1 && hiBin>=0 && hiBin <=180 && (abs(Bmass-5.27932)<0.08)";
//TString SignalCut = "$CommonCut && (abs(Bmass-5.27932)<0.08)";
//TString SideBandCut = "$CommonCut && (abs(Bmass-5.27932)>0.15 && abs(Bmass-5.27932)<0.25)";
//TString GenCut = "$CommonCut && (abs(Bmass-5.27932)<0.08 && Bgen==23333)";
//TString weight = "pthatweight*Ncoll";
TString weight="pthatweight*Ncoll*(1.034350*TMath::Exp(-0.000844*(PVz+3.502992)*(PVz+3.502992)))*(0.715021+0.039896*Bgenpt-0.000834*Bgenpt*Bgenpt+0.000006*Bgenpt*Bgenpt*Bgenpt)";

void pvz()
{
  gStyle->SetOptStat(0);

  TFile* inputdata = new TFile(inputdataname.Data());
  TFile* inputmc = new TFile(inputmcname.Data());
  
  TTree* data = (TTree*) inputdata->Get("Bfinder/ntKp");
  data->AddFriend("hltanalysis/HltTree");
  data->AddFriend("hiEvtAnalyzer/HiTree");
  TTree* mc = (TTree*) inputmc->Get("Bfinder/ntKp");
  mc->AddFriend("hltanalysis/HltTree");
  mc->AddFriend("hiEvtAnalyzer/HiTree");
  mc->AddFriend("Bfinder/ntGen");
  
  TH1F* Bpt_data = new TH1F("Bpt_data","",nBins,ptBins);
  TH1F* Bpt_mc = new TH1F("Bpt_mc","",nBins,ptBins);
  
  data->Project("Bpt_data","Bpt",TCut(CommonCut));
  mc->Project("Bpt_mc","Bpt",TCut(CommonCut)*TCut(weight)&&"(Bgen==23333)");

  divideBinWidth(Bpt_data);
  Bpt_data->Scale(1.0/Bpt_data->Integral("width"));  
  divideBinWidth(Bpt_mc);
  Bpt_mc->Scale(1.0/Bpt_mc->Integral("width"));

  TH1F* Bpt_ratio = (TH1F*) Bpt_data->Clone("Bpt_ratio");
  Bpt_ratio->Sumw2();
  Bpt_ratio->Divide(Bpt_mc);

  TCanvas* c1 = new TCanvas("","",600,600);
  c1->cd();
  Bpt_ratio->Draw();
  
  c1->SaveAs("Bpt_ratio.png");
  
  return;
  
}

int main()
{
  pvz();
  return 0;
}
