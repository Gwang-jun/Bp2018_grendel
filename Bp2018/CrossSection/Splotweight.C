#include "uti.h"
#include "ntuple.h"

TString weightname = "/raid5/data/gwangjun/weights.root";
TString inputmcname = "/raid5/data/gwangjun/crab_Bfinder_20190624_Hydjet_Pythia8_Official_BuToJpsiK_1033p1_pt3tkpt0p7dls2_allpthat_pthatweight_MuonJSON_BDT.root";
TString outputmcname = "/raid5/data/gwangjun/MC_splot_BDT_trk.root";

void Splotweight()
{
  TFile* inputweight = new TFile(weightname.Data());
  TH1F* weights_BDT_pt_5_7 = (TH1F*) inputweight->Get("weights_BDT_pt_5_7");
  TH1F* weights_BDT_pt_7_10 = (TH1F*) inputweight->Get("weights_BDT_pt_7_10");
  TH1F* weights_BDT_pt_10_15 = (TH1F*) inputweight->Get("weights_BDT_pt_10_15");
  TH1F* weights_BDT_pt_15_20 = (TH1F*) inputweight->Get("weights_BDT_pt_15_20");
  TH1F* weights_BDT_pt_20_30 = (TH1F*) inputweight->Get("weights_BDT_pt_20_30");
  TH1F* weights_BDT_pt_30_40 = (TH1F*) inputweight->Get("weights_BDT_pt_30_40");
  TH1F* weights_BDT_pt_40_50 = (TH1F*) inputweight->Get("weights_BDT_pt_40_50");
  TH1F* weights_BDT_pt_50_60 = (TH1F*) inputweight->Get("weights_BDT_pt_50_60");

  TH1F* weights_Btrk1Eta = (TH1F*) inputweight->Get("weights_Btrk1eta");
  TH1F* weights_Btrk1Y = (TH1F*) inputweight->Get("weights_Btrk1Y");
  TH1F* weights_Btrk1Pt = (TH1F*) inputweight->Get("weights_Btrk1pt");
  TH1F* weights_Btrk1Dz1 = (TH1F*) inputweight->Get("weights_Btrk1Dz1");
  TH1F* weights_Btrk1DzError1 = (TH1F*) inputweight->Get("weights_Btrk1DzError1");
  TH1F* weights_Btrk1Dxy1 = (TH1F*) inputweight->Get("weights_Btrk1Dxy1");
  TH1F* weights_Btrk1DxyError1 = (TH1F*) inputweight->Get("weights_Btrk1DxyError1");
  

  TFile* inputMC = new TFile(inputmcname.Data());
  TTree* ntMC = (TTree*)inputMC->Get("Bfinder/ntKp");
  ntMC->AddFriend("hltanalysis/HltTree");
  ntMC->AddFriend("hiEvtAnalyzer/HiTree");
  ntMC->AddFriend("Bfinder/ntGen");
  ntMC->AddFriend("skimanalysis/HltTree");
  ntMC->AddFriend("BDT");

  //setbranchaddress(inputMC,ntMC);

  Int_t Bsize;
  Float_t Bpt[MAX_XB];
  Double_t BDT_5_7[MAX_XB];
  Double_t BDT_7_10[MAX_XB];
  Double_t BDT_10_15[MAX_XB];
  Double_t BDT_15_20[MAX_XB];
  Double_t BDT_20_30[MAX_XB];
  Double_t BDT_30_40[MAX_XB];
  Double_t BDT_40_50[MAX_XB];
  Double_t BDT_50_60[MAX_XB];

  Float_t Btrk1Eta[MAX_XB];
  Float_t Btrk1Y[MAX_XB];
  Float_t Btrk1Pt[MAX_XB];
  Float_t Btrk1Dz1[MAX_XB];
  Float_t Btrk1DzError1[MAX_XB];
  Float_t Btrk1Dxy1[MAX_XB];
  Float_t Btrk1DxyError1[MAX_XB];

  ntMC->SetBranchAddress("Bsize",&Bsize);
  ntMC->SetBranchAddress("Bpt",Bpt);
  ntMC->SetBranchAddress("BDT_5_7",BDT_5_7);
  ntMC->SetBranchAddress("BDT_7_10",BDT_7_10);
  ntMC->SetBranchAddress("BDT_10_15",BDT_10_15);
  ntMC->SetBranchAddress("BDT_15_20",BDT_15_20);
  ntMC->SetBranchAddress("BDT_20_30",BDT_20_30);
  ntMC->SetBranchAddress("BDT_30_40",BDT_30_40);
  ntMC->SetBranchAddress("BDT_40_50",BDT_40_50);
  ntMC->SetBranchAddress("BDT_50_60",BDT_50_60);

  ntMC->SetBranchAddress("Btrk1Eta",Btrk1Eta);
  ntMC->SetBranchAddress("Btrk1Y",Btrk1Y);
  ntMC->SetBranchAddress("Btrk1Pt",Btrk1Pt);
  ntMC->SetBranchAddress("Btrk1Dz1",Btrk1Dz1);
  ntMC->SetBranchAddress("Btrk1DzError1",Btrk1DzError1);
  ntMC->SetBranchAddress("Btrk1Dxy1",Btrk1Dxy1);
  ntMC->SetBranchAddress("Btrk1DxyError1",Btrk1DxyError1);


  TFile* outputMC = new TFile(outputmcname.Data(),"recreate");
  outputMC->cd();
  
  TTree* splot = new TTree("splot","splot");

  Int_t BsizeNew;
  Float_t splotweight[MAX_XB];
  Float_t Btrk1Eta_weight[MAX_XB];
  Float_t Btrk1Y_weight[MAX_XB];
  Float_t Btrk1Pt_weight[MAX_XB];
  Float_t Btrk1Dz1_weight[MAX_XB];
  Float_t Btrk1DzError1_weight[MAX_XB];
  Float_t Btrk1Dxy1_weight[MAX_XB];
  Float_t Btrk1DxyError1_weight[MAX_XB];

  splot->Branch("Bsize",&BsizeNew,"Bsize/I");
  splot->Branch("splotweight",splotweight,"splotweight[Bsize]/F");
  splot->Branch("Btrk1Eta_weight",Btrk1Eta_weight,"Btrk1Eta_weight[Bsize]/F");
  splot->Branch("Btrk1Y_weight",Btrk1Y_weight,"Btrk1Y_weight[Bsize]/F");
  splot->Branch("Btrk1Pt_weight",Btrk1Pt_weight,"Btrk1Pt_weight[Bsize]/F");
  splot->Branch("Btrk1Dz1_weight",Btrk1Dz1_weight,"Btrk1Dz1_weight[Bsize]/F");
  splot->Branch("Btrk1DzError1_weight",Btrk1DzError1_weight,"Btrk1DzError1_weight[Bsize]/F");
  splot->Branch("Btrk1Dxy1_weight",Btrk1Dxy1_weight,"Btrk1Dxy1_weight[Bsize]/F");
  splot->Branch("Btrk1DxyError1_weight",Btrk1DxyError1_weight,"Btrk1DxyError1_weight[Bsize]/F");

  //for(int evt=0;evt<10000;evt++)
  for(int evt=0;evt<ntMC->GetEntries();evt++)
    {
      if(evt%100000==0) std::cout<<std::setiosflags(std::ios::left)<<"  [ \033[1;36m"<<std::setw(10)<<evt<<"\033[0m"<<" / "<<std::setw(10)<<ntMC->GetEntries()<<" ] "<<"\033[1;36m"<<(int)(100.*evt/ntMC->GetEntries())<<"%\033[0m"<<"\r"<<std::flush;

      //if(evt%100000==0) std::cout<<std::setiosflags(std::ios::left)<<"  [ \033[1;36m"<<std::setw(10)<<evt<<"\033[0m"<<" / "<<std::setw(10)<<ntMC->GetEntries()<<" ] "<<"\033[1;36m"<<Form("%.0f",100.*evt/ntMC->GetEntries())<<"%\033[0m"<<"\r"<<std::flush;

      ntMC->GetEntry(evt);
      BsizeNew = Bsize;
      for(int i=0;i<Bsize;i++)
	{
	  if(Bpt[i]>=ptBins[0] && Bpt[i]<ptBins[1]) splotweight[i] = weights_BDT_pt_5_7->GetBinContent(weights_BDT_pt_5_7->GetXaxis()->FindBin(BDT_5_7[i]));
	  if(Bpt[i]>=ptBins[1] && Bpt[i]<ptBins[2]) splotweight[i] = weights_BDT_pt_7_10->GetBinContent(weights_BDT_pt_7_10->GetXaxis()->FindBin(BDT_7_10[i]));
	  if(Bpt[i]>=ptBins[2] && Bpt[i]<ptBins[3]) splotweight[i] = weights_BDT_pt_10_15->GetBinContent(weights_BDT_pt_10_15->GetXaxis()->FindBin(BDT_10_15[i]));
	  if(Bpt[i]>=ptBins[3] && Bpt[i]<ptBins[4]) splotweight[i] = weights_BDT_pt_15_20->GetBinContent(weights_BDT_pt_15_20->GetXaxis()->FindBin(BDT_15_20[i]));
	  if(Bpt[i]>=ptBins[4] && Bpt[i]<ptBins[5]) splotweight[i] = weights_BDT_pt_20_30->GetBinContent(weights_BDT_pt_20_30->GetXaxis()->FindBin(BDT_20_30[i]));
	  if(Bpt[i]>=ptBins[5] && Bpt[i]<ptBins[6]) splotweight[i] = weights_BDT_pt_30_40->GetBinContent(weights_BDT_pt_30_40->GetXaxis()->FindBin(BDT_30_40[i]));
	  if(Bpt[i]>=ptBins[6] && Bpt[i]<ptBins[7]) splotweight[i] = weights_BDT_pt_40_50->GetBinContent(weights_BDT_pt_40_50->GetXaxis()->FindBin(BDT_40_50[i]));
	  if(Bpt[i]>=ptBins[7] && Bpt[i]<ptBins[8]) splotweight[i] = weights_BDT_pt_50_60->GetBinContent(weights_BDT_pt_50_60->GetXaxis()->FindBin(BDT_50_60[i]));
      
	  Btrk1Eta_weight[i] = weights_Btrk1Eta->GetBinContent(weights_Btrk1Eta->GetXaxis()->FindBin(Btrk1Eta[i]));
	  Btrk1Y_weight[i] = weights_Btrk1Y->GetBinContent(weights_Btrk1Y->GetXaxis()->FindBin(Btrk1Y[i]));
	  Btrk1Pt_weight[i] = weights_Btrk1Pt->GetBinContent(weights_Btrk1Pt->GetXaxis()->FindBin(Btrk1Pt[i]));
	  Btrk1Dz1_weight[i] = weights_Btrk1Dz1->GetBinContent(weights_Btrk1Dz1->GetXaxis()->FindBin(Btrk1Dz1[i]));
	  Btrk1DzError1_weight[i] = weights_Btrk1DzError1->GetBinContent(weights_Btrk1DzError1->GetXaxis()->FindBin(Btrk1DzError1[i]));
	  Btrk1Dxy1_weight[i] = weights_Btrk1Dxy1->GetBinContent(weights_Btrk1Dxy1->GetXaxis()->FindBin(Btrk1Dxy1[i]));
	  Btrk1DxyError1_weight[i] = weights_Btrk1DxyError1->GetBinContent(weights_Btrk1DxyError1->GetXaxis()->FindBin(Btrk1DxyError1[i]));
	}

      splot->Fill();
    }

  outputMC->Write();
  outputMC->Close();

  return;

}

int main()
{
  Splotweight();
  return 0;
}
