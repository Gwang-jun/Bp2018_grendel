#include <TTree.h>
#include <TMath.h>
#include <iostream>
#include "parameters.h"

#ifndef MAX_XB
#define MAX_XB      20000
#endif

//float  pthatweight;
//float  Ncoll;
float  PVz;
int    HiBin;
int    hiBin;
int    HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1;
int    pprimaryVertexFilter;
int    phfCoincFilter2Th4;
int    pclusterCompatibilityFilter;
double BDT_5_7[MAX_XB];
double BDT_7_10[MAX_XB];
double BDT_10_15[MAX_XB];
double BDT_15_20[MAX_XB];
double BDT_20_30[MAX_XB];
double BDT_30_40[MAX_XB];
double BDT_40_50[MAX_XB];
double BDT_50_60[MAX_XB];

int    Bsize;
float  Bgen[MAX_XB];
//float  Bpt[MAX_XB];
float  Bpt;
float  Bgenpt[MAX_XB];
float  Balpha[MAX_XB];
float  Blxy[MAX_XB];
float  Btrk1Pt[MAX_XB];
float  Btrk2Pt[MAX_XB];
//float  Bmass[MAX_XB];
float  Bmass;
float  Btktkmass[MAX_XB];
float  Bujpt[MAX_XB];
float  Bchi2cl[MAX_XB];
float  Bdtheta[MAX_XB];
float  Bujphi[MAX_XB];
float  Btrk1Phi[MAX_XB];
float  Btrk2Phi[MAX_XB];
float  BsvpvDistance[MAX_XB];
float  BsvpvDisErr[MAX_XB];
float  BsvpvDistance_2D[MAX_XB];
float  BsvpvDisErr_2D[MAX_XB];
float  Bujeta[MAX_XB];
float  Btrk1Eta[MAX_XB];
float  Btrk2Eta[MAX_XB];
float  Btrk1Dxy1[MAX_XB];
float  Btrk2Dxy1[MAX_XB];
float  Btrk1DxyError1[MAX_XB];
float  Btrk2DxyError1[MAX_XB];
float  Btrk1Chi2ndf[MAX_XB];
float  Btrk2Chi2ndf[MAX_XB];
float  Btrk1nStripLayer[MAX_XB];
float  Btrk2nStripLayer[MAX_XB];
float  Btrk1nPixelLayer[MAX_XB];
float  Btrk2nPixelLayer[MAX_XB];
//float  By[MAX_XB];
float  By;

bool   Bmu1TMOneStationTight[MAX_XB];
int    Bmu1InPixelLayer[MAX_XB];
int    Bmu1InStripLayer[MAX_XB];
float  Bmu1dxyPV[MAX_XB];
float  Bmu1dzPV[MAX_XB];
bool   Bmu1isGlobalMuon[MAX_XB];
float  Bmu1eta[MAX_XB];
float  Bmu1pt[MAX_XB];
bool   Bmu2TMOneStationTight[MAX_XB];
int    Bmu2InPixelLayer[MAX_XB];
int    Bmu2InStripLayer[MAX_XB];
float  Bmu2dxyPV[MAX_XB];
float  Bmu2dzPV[MAX_XB];
bool   Bmu2isGlobalMuon[MAX_XB];
float  Bmu2eta[MAX_XB];
float  Bmu2pt[MAX_XB];
float  Bmumumass[MAX_XB];
bool   Btrk1highPurity[MAX_XB];
float  Btrk1PixelHit[MAX_XB];
float  Btrk1StripHit[MAX_XB];
float  Btrk1PtErr[MAX_XB];
bool   Btrk2highPurity[MAX_XB];
float  Btrk2PixelHit[MAX_XB];
float  Btrk2StripHit[MAX_XB];
float  Btrk2PtErr[MAX_XB];

bool   Bmu1isTriggered[MAX_XB];
bool   Bmu2isTriggered[MAX_XB];
bool   Bmu1SoftMuID[MAX_XB];
bool   Bmu2SoftMuID[MAX_XB];
bool   Bmu1isAcc[MAX_XB];
bool   Bmu2isAcc[MAX_XB];

void setbranchaddress(TFile* ffile,TTree* fnt)
{
  fnt->SetBranchAddress("HiBin", &HiBin);
  fnt->SetBranchAddress("hiBin", &hiBin);
  //fnt->SetBranchAddress("Bpt", Bpt);
  fnt->SetBranchAddress("Bpt", &Bpt);
  //fnt->SetBranchAddress("Bmass", Bmass);
  fnt->SetBranchAddress("Bmass", &Bmass);  
  //fnt->SetBranchAddress("By", By);
  fnt->SetBranchAddress("By", &By);
}

int findBptbin(float Bpt)
{
  int n;
  for(int i=0;i<10;i++)
    {
      if(Bpt>5.0+0.5*i && Bpt<5.0+0.5*(i+1))
	{
	  n=i+1;
	  break;
	}
    }
  for(int j=0;j<50;j++)
    {
      if(Bpt>10.0+j && Bpt<10.0+(j+1))
	{
	  n=j+11;
	  break;
	}
    }

  /*
  int n;
  for(int i=0;i<nBinsFine;i++)
    {
      if(Bpt>5.0+45.0/nBinsFine*i && Bpt<5.0+45.0/nBinsFine*(i+1))
	{
	  n=i+1;
	  break;
	}
    }
  */

  return n;
}

int findBybin(float By)
{
  int n;
  if(abs(By)>0.0 && abs(By)<0.5) n=1;
  if(abs(By)>0.5 && abs(By)<1.0) n=2;
  if(abs(By)>1.0 && abs(By)<1.5) n=3;
  if(abs(By)>1.5 && abs(By)<2.0) n=4;
  if(abs(By)>2.0 && abs(By)<2.4) n=5;
  
  /*
  for(int i=0;i<nBinsYFine;i++)
    {
      if(By>-2.4+4.8/nBinsYFine*i && By<-2.4+4.8/nBinsYFine*(i+1))
	{
	  n=i+1;
	  break;
	}
    }
  */
  
  return n;
}

int findCentbin(float HiBin)
{
  int n;
  for(int i=0;i<40;i++)
    {
      if(HiBin>=60+3*i && HiBin<60+3*(i+1))
	{
	  n=i+1;
	  break;
	}
    }
  return n;
}
