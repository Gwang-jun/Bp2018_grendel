#include "uti.h"
#include "TF1.h"
#include <TFitResultPtr.h>
#include "ntuple2.h"

//TString inputdataname = "/raid5/data/gwangjun/selected_data_ntKp_PbPb_2018_new_train_BDT.root";
TString inputdataname = "/raid5/data/gwangjun/selected_data_ntKp_PbPb_2018_MuonJSON_BDT.root";

TString seldata;
TString selmc;
TString selmceff;
TString selmcgen;
TString collisionsystem;
Float_t hiBinMin,hiBinMax,centMin,centMax;

bool useTNP = 1;

void AverageEff(int usePbPb=0, TString inputdata="" , TString inputmc="", TString trgselection="1",  TString cut="", TString cutmcgen="", int isMC=0, Double_t luminosity=1., int doweight=0, TString collsyst="PbPb", TString outputfile="", TString npfit="0", int doDataCor = 0, Float_t centmin=0., Float_t centmax=90.)
{
  //TString inputeffname = Form("ROOTfiles/MCstudiesPbPb2D_Bptweight_Cent%.0f-%.0f.root",centmin,centmax);
  //TString outputeffname = Form("ROOTfiles/MCstudiesPbPbAverage_Bptweight_Cent%.0f-%.0f.root",centmin,centmax);
  TString inputeffname = Form("ROOTfiles/MCstudiesPbPb_2D_Cent%.0f-%.0f.root",centmin,centmax);
  TString outputeffname = Form("ROOTfiles/MCstudiesPbPbAverage_2D_pt650_Cent%.0f-%.0f.root",centmin,centmax);
  TString inputTNPname = Form("ROOTfiles/TNP2D_Bplus_Cent%.0f-%.0f.root",centmin,centmax);

  hiBinMin = centmin*2;
  hiBinMax = centmax*2;
  centMin = centmin;
  centMax = centmax;

  if (!(usePbPb==1||usePbPb==0)) std::cout<<"ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!, you are using a non valid isPbPb option"<<std::endl;
  bool isPbPb=(bool)(usePbPb);

  gStyle->SetTextSize(0.05);
  gStyle->SetTextFont(42);
  gStyle->SetPadRightMargin(0.043);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.145);
  gStyle->SetTitleX(.0f);
  
  //TFile* inf = new TFile(inputdata.Data());
  TFile* inf = new TFile(inputdataname.Data());
  TTree* nt = (TTree*)inf->Get("ntKp");
  //nt->AddFriend("hltanalysis/HltTree");
  //nt->AddFriend("hiEvtAnalyzer/HiTree");
  //nt->AddFriend("skimanalysis/HltTree");
  //nt->AddFriend("BDT");
  
  setbranchaddress(inf,nt);
  
  TFile* Eff2Dfile = new TFile(inputeffname.Data());
  TH2D* hEff2D = (TH2D*)Eff2Dfile->Get("hEff2D");
  TH2D* invEff2D = (TH2D*)Eff2Dfile->Get("invEff2D");
  TFile* TNP2D = new TFile(inputTNPname.Data());
  TH2D* tnp_scale = (TH2D*)TNP2D->Get("tnp_scale");
  TH2D* tnp_total_d = (TH2D*)TNP2D->Get("tnp_total_d");
  TH2D* tnp_total_u = (TH2D*)TNP2D->Get("tnp_total_u");

  int nevts[nBins];
  float sum[nBins], sumerr[nBins], sum_d[nBins], sumerr_d[nBins], sum_u[nBins], sumerr_u[nBins];
  for(int i=0;i<nBins;i++)
    {
      nevts[i]=0;
      sum[i]=0.0;
      sumerr[i]=0.0;
      sum_d[i]=0.0;
      sumerr_d[i]=0.0;
      sum_u[i]=0.0;
      sumerr_u[i]=0.0;
    }

  for(int i=0;i<nt->GetEntries();i++)
    {
      nt->GetEntry(i);
      for(int k=0;k<nBins;k++)
	{
	  if(!(TMath::Abs(Bmass-5.27932)<0.08)) continue;
	  if(!(Bpt>ptBins[k] && Bpt<ptBins[k+1] && HiBin>=centMin && HiBin<centMax)) continue;
	  if(hEff2D->GetBinContent(invEff2D->GetBin(findBptbin(Bpt),findBybin(By)))==0) continue;
	  nevts[k]++;

	  if(!useTNP)
	    {
	      sum[k]+=invEff2D->GetBinContent(invEff2D->GetBin(findBptbin(Bpt),findBybin(By)));
	      sumerr[k]+=invEff2D->GetBinError(invEff2D->GetBin(findBptbin(Bpt),findBybin(By)))*invEff2D->GetBinError(invEff2D->GetBin(findBptbin(Bpt),findBybin(By)));
	    }

	  if(useTNP)
	    {
	      sum[k]+=invEff2D->GetBinContent(invEff2D->GetBin(findBptbin(Bpt),findBybin(By)))/tnp_scale->GetBinContent(tnp_scale->GetBin(findBptbin(Bpt),findBybin(By)));
	      sumerr[k]+=(invEff2D->GetBinError(invEff2D->GetBin(findBptbin(Bpt),findBybin(By)))/tnp_scale->GetBinContent(tnp_scale->GetBin(findBptbin(Bpt),findBybin(By))))*(invEff2D->GetBinError(invEff2D->GetBin(findBptbin(Bpt),findBybin(By)))/tnp_scale->GetBinContent(tnp_scale->GetBin(findBptbin(Bpt),findBybin(By))));
	      
	      sum_d[k]+=invEff2D->GetBinContent(invEff2D->GetBin(findBptbin(Bpt),findBybin(By)))/(tnp_scale->GetBinContent(tnp_scale->GetBin(findBptbin(Bpt),findBybin(By)))*(1.0-tnp_total_d->GetBinContent(tnp_total_d->GetBin(findBptbin(Bpt),findBybin(By)))));
	      sumerr_d[k]+=(invEff2D->GetBinError(invEff2D->GetBin(findBptbin(Bpt),findBybin(By)))/(tnp_scale->GetBinContent(tnp_scale->GetBin(findBptbin(Bpt),findBybin(By)))*(1.0-tnp_total_d->GetBinContent(tnp_total_d->GetBin(findBptbin(Bpt),findBybin(By))))))*(invEff2D->GetBinError(invEff2D->GetBin(findBptbin(Bpt),findBybin(By)))/(tnp_scale->GetBinContent(tnp_scale->GetBin(findBptbin(Bpt),findBybin(By)))*(1.0-tnp_total_d->GetBinContent(tnp_total_d->GetBin(findBptbin(Bpt),findBybin(By))))));
	      
	      sum_u[k]+=invEff2D->GetBinContent(invEff2D->GetBin(findBptbin(Bpt),findBybin(By)))/(tnp_scale->GetBinContent(tnp_scale->GetBin(findBptbin(Bpt),findBybin(By)))*(1.0+tnp_total_u->GetBinContent(tnp_total_u->GetBin(findBptbin(Bpt),findBybin(By)))));
	      sumerr_u[k]+=(invEff2D->GetBinError(invEff2D->GetBin(findBptbin(Bpt),findBybin(By)))/(tnp_scale->GetBinContent(tnp_scale->GetBin(findBptbin(Bpt),findBybin(By)))*(1.0+tnp_total_u->GetBinContent(tnp_total_u->GetBin(findBptbin(Bpt),findBybin(By))))))*(invEff2D->GetBinError(invEff2D->GetBin(findBptbin(Bpt),findBybin(By)))/(tnp_scale->GetBinContent(tnp_scale->GetBin(findBptbin(Bpt),findBybin(By)))*(1.0+tnp_total_u->GetBinContent(tnp_total_u->GetBin(findBptbin(Bpt),findBybin(By))))));
	    }

	}
   }


  /*
  TFile* Eff1Dfile = new TFile(inputeffname.Data());
  TH1D* hEff1D = (TH1D*)Eff1Dfile->Get("hEff");
  TH1D* invEff1D = (TH1D*)Eff1Dfile->Get("invEff");

  int nevts[nBins];
  float sum[nBins], sumerr[nBins];
  for(int i=0;i<nBins;i++)
    {
      nevts[i]=0;
      sum[i]=0.0;
      sumerr[i]=0.0;
    }


  for(int i=0;i<nt->GetEntries();i++)
    {
      nt->GetEntry(i);
      for(int k=0;k<nBins;k++)
	{
	  if(!(TMath::Abs(Bmass-5.27932)<0.08)) continue;
	  if(!(Bpt>ptBins[k] && Bpt<ptBins[k+1] && HiBin>=centmin && HiBin<centmax)) continue;
	  if(hEff1D->GetBinContent(invEff1D->GetBin(findBptbin(Bpt)))==0) continue;
	  nevts[k]++;
	  sum[k]+=invEff1D->GetBinContent(invEff1D->GetBin(findBptbin(Bpt)));
	  sumerr[k]+=invEff1D->GetBinError(invEff1D->GetBin(findBptbin(Bpt)))*invEff1D->GetBinError(invEff1D->GetBin(findBptbin(Bpt)));
	}
   }
  */

  /*
  for(int i=0;i<nt->GetEntries();i++)
    {
      nt->GetEntry(i);
      for(int k=0;k<nBins;k++)
	{
	  if(!(TMath::Abs(Bmass-5.27932)<0.08)) continue;
	  if(!(Bpt>ptBins[k] && Bpt<ptBins[k+1] && HiBin>=centmin && HiBin<centmax)) continue;
	  if(hEff1D->GetBinContent(invEff1D->GetBin(findCentbin(HiBin)))==0) continue;
	  nevts[k]++;
	  sum[k]+=invEff1D->GetBinContent(invEff1D->GetBin(findCentbin(HiBin)));
	  sumerr[k]+=invEff1D->GetBinError(invEff1D->GetBin(findCentbin(HiBin)))*invEff1D->GetBinError(invEff1D->GetBin(findCentbin(HiBin)));
	}
   }
  */
  
  TH1F* invEffave = new TH1F("invEffave","",nBins,ptBins);
  TH1F* invEffave_d = new TH1F("invEffave_d","",nBins,ptBins);
  TH1F* invEffave_u = new TH1F("invEffave_u","",nBins,ptBins);

  for(int i=0;i<nBins;i++)
    {
      sum[i]=sum[i]/nevts[i];
      sumerr[i]=sqrt(sumerr[i])/nevts[i];
      invEffave->SetBinContent(i+1,sum[i]);
      invEffave->SetBinError(i+1,sumerr[i]);

      if(useTNP)
	{
	  sum_d[i]=sum_d[i]/nevts[i];
	  sumerr_d[i]=sqrt(sumerr_d[i])/nevts[i];
	  sum_u[i]=sum_u[i]/nevts[i];
	  sumerr_u[i]=sqrt(sumerr_u[i])/nevts[i];

	  invEffave_d->SetBinContent(i+1,sum_u[i]);
	  invEffave_d->SetBinError(i+1,sumerr_u[i]);
	  invEffave_u->SetBinContent(i+1,sum_d[i]);
	  invEffave_u->SetBinError(i+1,sumerr_d[i]);
	}

      std::cout<<"ptbin "<<ptBins[i]<<"-"<<ptBins[i+1]<<" Nevts: "<<nevts[i]<<", Average 1/(acc*eff) (Cent): "<<sum[i]<<" pm "<<sumerr[i]<<std::endl;
      if(useTNP)
	{
	  std::cout<<"ptbin "<<ptBins[i]<<"-"<<ptBins[i+1]<<" Nevts: "<<nevts[i]<<", Average 1/(acc*eff) (Low): "<<sum_u[i]<<" pm "<<sumerr_u[i]<<", TNP syst (Low): "<<((sum[i]-sum_u[i])/sum[i])*100<<"%"<<std::endl;
	  std::cout<<"ptbin "<<ptBins[i]<<"-"<<ptBins[i+1]<<" Nevts: "<<nevts[i]<<", Average 1/(acc*eff) (High): "<<sum_d[i]<<" pm "<<sumerr_d[i]<<", TNP syst (High): "<<((sum_d[i]-sum[i])/sum[i])*100<<"%"<<std::endl;
	}
    }

  TFile* outf = new TFile(outputeffname.Data(),"recreate");
  outf->cd();
  invEffave->Write();
  if(useTNP)
    {
      invEffave_d->Write();
      invEffave_u->Write();
    }
  outf->Close(); 
}

int main(int argc, char *argv[])
{
  if(argc==16)
    {
      AverageEff(atoi(argv[1]), argv[2], argv[3], argv[4], argv[5], argv[6], atoi(argv[7]), atof(argv[8]), atoi(argv[9]), argv[10], argv[11], argv[12], atoi(argv[13]), atof(argv[14]), atof(argv[15]));
      return 0;
    }
  else if(argc==14)
    {
      AverageEff(atoi(argv[1]), argv[2], argv[3], argv[4], argv[5], argv[6], atoi(argv[7]), atof(argv[8]), atoi(argv[9]), argv[10], argv[11], argv[12], atoi(argv[13]));
      return 0;
    }
  else
    {
      std::cout << "Wrong number of inputs" << std::endl;
      return 1;
    }
}
