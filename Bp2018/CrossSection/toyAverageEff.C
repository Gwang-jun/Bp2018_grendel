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

using namespace std;

using std::cout;
using std::endl;

void toyAverageEff(int CentMin, int CentMax)
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  
  const int yBinN = 5;
  double yBinning[yBinN+1] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.4};
  
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
  
  TString FileName = Form("ROOTfiles/MCstudiesPbPb_Fid2D_Cent%d-%d.root",CentMin,CentMax);
  TFile * fin = new TFile(FileName.Data());
  fin->cd();
  TH2D * EffBptByInv = (TH2D * ) fin->Get("invEff2D");
  TH2D * EffBptBy = (TH2D * ) fin->Get("hEff2D");
  TH2D * AccBptByInv = (TH2D * ) fin->Get("invAcc2D");
  TH2D * AccBptBy = (TH2D * ) fin->Get("hAcc2D");
  TH2D * EffonlyBptByInv = (TH2D * ) fin->Get("invEffonly2D");
  TH2D * EffonlyBptBy = (TH2D * ) fin->Get("hEffonly2D");

  int NTrials = 10000;
  
  double EffInv2DMapCenter[BptBin][yBinN];
  double EffInv2DMapError[BptBin][yBinN];
  double AccInv2DMapCenter[BptBin][yBinN];
  double AccInv2DMapError[BptBin][yBinN];
  double EffonlyInv2DMapCenter[BptBin][yBinN];
  double EffonlyInv2DMapError[BptBin][yBinN];
  
  //loop through the 2D map histogram//  
  
  for(int i = 0; i < BptBin; i++){    
    for(int j = 0; j < yBinN; j++){      
      EffInv2DMapCenter[i][j] = EffBptByInv->GetBinContent(i+1,j+1);
      EffInv2DMapError[i][j] = EffBptByInv->GetBinError(i+1,j+1);		
      AccInv2DMapCenter[i][j] = AccBptByInv->GetBinContent(i+1,j+1);
      AccInv2DMapError[i][j] = AccBptByInv->GetBinError(i+1,j+1);		
      EffonlyInv2DMapCenter[i][j] = EffonlyBptByInv->GetBinContent(i+1,j+1);
      EffonlyInv2DMapError[i][j] = EffonlyBptByInv->GetBinError(i+1,j+1);		

      cout<<"i = "<<i<<" j = "<<j<<" Acc*Eff = "<<EffInv2DMapCenter[i][j]<<" Acc*Eff Err = "<<EffInv2DMapError[i][j]<<" ("<<100*EffInv2DMapError[i][j]/EffInv2DMapCenter[i][j]<<"%)"<<endl;
      cout<<"i = "<<i<<" j = "<<j<<" Acc = "<<AccInv2DMapCenter[i][j]<<" Acc Err = "<<AccInv2DMapError[i][j]<<" ("<<100*AccInv2DMapError[i][j]/AccInv2DMapCenter[i][j]<<"%)"<<endl;
      cout<<"i = "<<i<<" j = "<<j<<" Eff = "<<EffonlyInv2DMapCenter[i][j]<<" Eff Err = "<<EffonlyInv2DMapError[i][j]<<" ("<<100*EffonlyInv2DMapError[i][j]/EffonlyInv2DMapCenter[i][j]<<"%)"<<endl;
      cout<<endl;
    }   
  }
  

  TFile * fout = new TFile(Form("ROOTfiles/GenStatSyst_Cent%d-%d_AccEff.root",CentMin,CentMax),"RECREATE");
  fout->cd();
  
  TH2D * EffBptByInvTrial;
  double GeneratedEffInv;
  TH2D * AccBptByInvTrial;
  double GeneratedAccInv;
  TH2D * EffonlyBptByInvTrial;
  double GeneratedEffonlyInv;
  
  for(int iTrial=0;iTrial < NTrials; iTrial++){    
    //cout << "Now working on Trial  " << iTrial <<  "  out of  " << NTrials  << endl;    
    EffBptByInvTrial = new TH2D(Form("EffBptByInvTrial%d",iTrial),"",BptBin,BptBinning,yBinN,yBinning);    
    AccBptByInvTrial = new TH2D(Form("AccBptByInvTrial%d",iTrial),"",BptBin,BptBinning,yBinN,yBinning);    
    EffonlyBptByInvTrial = new TH2D(Form("EffonlyBptByInvTrial%d",iTrial),"",BptBin,BptBinning,yBinN,yBinning);    
    for(int i = 0; i < BptBin; i++)
      {      
	for(int j = 0; j < yBinN; j++)
	  {
	    if(EffInv2DMapCenter[i][j]!=0)
	      {
		GeneratedEffInv = gRandom->Gaus(EffInv2DMapCenter[i][j],EffInv2DMapError[i][j]);
		EffBptByInvTrial->SetBinContent(i+1,j+1,GeneratedEffInv);
	      }
	    if(AccInv2DMapCenter[i][j]!=0)
	      {
		GeneratedAccInv = gRandom->Gaus(AccInv2DMapCenter[i][j],AccInv2DMapError[i][j]);
		AccBptByInvTrial->SetBinContent(i+1,j+1,GeneratedAccInv);
	      }
	    if(EffonlyInv2DMapCenter[i][j]!=0)
	      {
		GeneratedEffonlyInv = gRandom->Gaus(EffonlyInv2DMapCenter[i][j],EffonlyInv2DMapError[i][j]);
		EffonlyBptByInvTrial->SetBinContent(i+1,j+1,GeneratedEffonlyInv);
	      }
	  }
      }
    EffBptByInvTrial->Write();
    AccBptByInvTrial->Write();
    EffonlyBptByInvTrial->Write();
  }
  fout->Close();

}

int main(int argc, char* argv[])
{
  toyAverageEff(atoi(argv[1]),atoi(argv[2]));
  return 0;
}
