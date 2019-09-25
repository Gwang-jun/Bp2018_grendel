#include "uti.h"
#include "TF1.h"
#include <TFitResultPtr.h>
#include "ntuple2.h"

TString inputdataname = "/raid5/data/gwangjun/selected_data_ntKp_PbPb_2018_new_train_BDT.root";
TString inputeffname = "ROOTfiles/MCstudiesPbPb.root";
TString outputeffname = "ROOTfiles/MCstudiesPbPbave.root";

TString seldata;
TString selmc;
TString selmceff;
TString selmcgen;
TString collisionsystem;
Float_t hiBinMin,hiBinMax,centMin,centMax;
double _ErrCor=1;

int _nBins = nBins;
double *_ptBins = ptBins;

void AverageEff(int usePbPb=0, TString inputdata="" , TString inputmc="", TString trgselection="1",  TString cut="", TString cutmcgen="", int isMC=0, Double_t luminosity=1., int doweight=0, TString collsyst="PbPb", TString outputfile="", TString npfit="0", int doDataCor = 0, Float_t centmin=0., Float_t centmax=90.)
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

  int nevts[nBins];
  double sum[nBins], sumerr[nBins];
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
	  if(!(Bpt>ptBins[k] && Bpt<ptBins[k+1])) continue;
	  if(hEff2D->GetBinContent(invEff2D->GetBin(findBptbin(Bpt),findBybin(By),0))==0) continue;
	  nevts[k]++;
	  sum[k]+=invEff2D->GetBinContent(invEff2D->GetBin(findBptbin(Bpt),findBybin(By),0));
	  sumerr[k]+=invEff2D->GetBinError(invEff2D->GetBin(findBptbin(Bpt),findBybin(By),0))*invEff2D->GetBinError(invEff2D->GetBin(findBptbin(Bpt),findBybin(By),0));
	}
   }

  /*
  for(int i=0;i<nt->GetEntries();i++)
    {
      nt->GetEntry(i);
      for(int j=0;j<Bsize;j++)
       {
	 if(!(pprimaryVertexFilter && phfCoincFilter2Th4 && pclusterCompatibilityFilter && Btrk1Pt[j]>0.9 && Bpt[j]>5.0 && (BsvpvDistance[j]/BsvpvDisErr[j])>2.0 && Bchi2cl[j]>0.05 && TMath::Abs(Btrk1Eta[j])<2.4 && TMath::Abs(By[j])<2.4 && TMath::Abs(PVz)<15 && Bmass[j]>5 && Bmass[j]<6 && TMath::Abs(Bmumumass[j]-3.096900)<0.15 && Bmu1SoftMuID[j] && Bmu2SoftMuID[j] && ((TMath::Abs(Bmu1eta[j])<1.2 && Bmu1pt[j]>3.5) || (TMath::Abs(Bmu1eta[j])>1.2 && TMath::Abs(Bmu1eta[j])<2.1 && Bmu1pt[j]>5.47-1.89*TMath::Abs(Bmu1eta[j])) || (TMath::Abs(Bmu1eta[j])>2.1 && TMath::Abs(Bmu1eta[j])<2.4 && Bmu1pt[j]>1.5)) && ((TMath::Abs(Bmu2eta[j])<1.2 && Bmu2pt[j]>3.5) || (TMath::Abs(Bmu2eta[j])>1.2 && TMath::Abs(Bmu2eta[j])<2.1 && Bmu2pt[j]>5.47-1.89*TMath::Abs(Bmu2eta[j])) || (TMath::Abs(Bmu2eta[j])>2.1 && TMath::Abs(Bmu2eta[j])<2.4 && Bmu2pt[j]>1.5)) && Bmu1isTriggered[j] && Bmu2isTriggered[j] && (Btrk1PixelHit[j]+Btrk1StripHit[j])>=11 && (Btrk1Chi2ndf[j]/(Btrk1nStripLayer[j]+Btrk1nPixelLayer[j]))<0.18 && TMath::Abs(Btrk1PtErr[j]/Btrk1Pt[j])<0.1 && ((Bpt[j]>5 && Bpt[j]<7 && (BsvpvDistance[j]/BsvpvDisErr[j])>12 && cos(Bdtheta[j])>0.95 && BDT_5_7[j]>0.07) || (Bpt[j]>7 && Bpt[j]<10 && BDT_7_10[j]>0.08) || (Bpt[j]>10 && Bpt[j]<15 && BDT_10_15[j]>0.09) || (Bpt[j]>15 && Bpt[j]<20 && BDT_15_20[j]>0.09) || (Bpt[j]>20 && Bpt[j]<30 && BDT_20_30[j]>0.07) || (Bpt[j]>30 && Bpt[j]<50 && BDT_30_50[j]>0.12) || (Bpt[j]>50 && Bpt[j]<100 && BDT_50_100[j]>0.24)))) continue;
	 for(int k=0;k<nBins;k++)
	   {
	     if(!(Bpt[j]>ptBins[k] && Bpt[j]<ptBins[k+1])) continue;
	     if(hEff2D->GetBinContent(invEff2D->GetBin(findBptbin(Bpt[j]),findBybin(By[j]),0))==0) continue;
	     nevts[k]++;
	     sum[k]+=invEff2D->GetBinContent(invEff2D->GetBin(findBptbin(Bpt[j]),findBybin(By[j]),0));
	     sumerr[k]+=invEff2D->GetBinError(invEff2D->GetBin(findBptbin(Bpt[j]),findBybin(By[j]),0))*invEff2D->GetBinError(invEff2D->GetBin(findBptbin(Bpt[j]),findBybin(By[j]),0));
	   }
       }
   }
  */

  TH1D* invEffave = new TH1D("invEffave","",nBins,ptBins);
  for(int i=0;i<nBins;i++)
    {
      sum[i]=sum[i]/nevts[i];
      sumerr[i]=sqrt(sumerr[i])/nevts[i];
      invEffave->SetBinContent(i+1,sum[i]);
      invEffave->SetBinError(i+1,sumerr[i]);
      std::cout<<"ptbin "<<i<<" Average 1/(acc*eff): "<<sum[i]<<" pm "<<sumerr[i]<<std::endl;      
    }

  TFile* outf = new TFile(outputeffname.Data(),"recreate");
  outf->cd();
  invEffave->Write();
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

double ErrorOnSigma(double width, double errwidth, double smear, double errsmearing){
  double squarederroronsigma=(1+smear)*(1+smear)*errwidth*errwidth+width*width*errsmearing*errsmearing;
  double erroronsigma=TMath::Sqrt(squarederroronsigma);
  return erroronsigma;
}
