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
double BDT_30_50[MAX_XB];
double BDT_50_100[MAX_XB];

int    Bsize;
float  Bgen[MAX_XB];
//float  Bpt[MAX_XB];
float  Bpt;
float  Bgenpt[MAX_XB];
float  Balpha[MAX_XB];
float  Blxy[MAX_XB];
float  Btrk1Pt[MAX_XB];
float  Btrk2Pt[MAX_XB];
float  Bmass[MAX_XB];
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
  //fnt->SetBranchAddress("pthatweight", &pthatweight); 
  //fnt->SetBranchAddress("Ncoll", &Ncoll);
  fnt->SetBranchAddress("PVz", &PVz);
  //fnt->SetBranchAddress("HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1", &HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1); 
  //fnt->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter); 
  //fnt->SetBranchAddress("phfCoincFilter2Th4", &phfCoincFilter2Th4); 
  //fnt->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter);       
  //fnt->SetBranchAddress("BDT_5_7", BDT_5_7);
  //fnt->SetBranchAddress("BDT_7_10", BDT_7_10);
  //fnt->SetBranchAddress("BDT_10_15", BDT_10_15);
  //fnt->SetBranchAddress("BDT_15_20", BDT_15_20);
  //fnt->SetBranchAddress("BDT_20_30", BDT_20_30);
  //fnt->SetBranchAddress("BDT_30_50", BDT_30_50);
  //fnt->SetBranchAddress("BDT_50_100", BDT_50_100);
  //fnt->SetBranchAddress("hiBin", &hiBin);
  //fnt->SetBranchAddress("Bsize", &Bsize);
  //fnt->SetBranchAddress("Bgen", Bgen);
  //fnt->SetBranchAddress("Bpt", Bpt);
  fnt->SetBranchAddress("Bpt", &Bpt);
  fnt->SetBranchAddress("Bgenpt", Bgenpt);
  fnt->SetBranchAddress("Balpha", Balpha);
  //fnt->SetBranchAddress("Blxy", Blxy);
  //fnt->SetBranchAddress("Btrk1Pt", Btrk1Pt);
  //fnt->SetBranchAddress("Btrk2Pt", Btrk2Pt);
  fnt->SetBranchAddress("Bmass", Bmass);
  //fnt->SetBranchAddress("Btktkmass", Btktkmass);
  //fnt->SetBranchAddress("Bujpt", Bujpt);
  fnt->SetBranchAddress("Bchi2cl", Bchi2cl);
  //fnt->SetBranchAddress("Bdtheta", Bdtheta);
  //fnt->SetBranchAddress("Bujphi", Bujphi);
  //fnt->SetBranchAddress("Btrk1Phi", Btrk1Phi);
  //fnt->SetBranchAddress("Btrk2Phi", Btrk2Phi);
  fnt->SetBranchAddress("BsvpvDistance", BsvpvDistance);
  //fnt->SetBranchAddress("BsvpvDisErr", BsvpvDisErr);
  //fnt->SetBranchAddress("BsvpvDistance_2D", BsvpvDistance_2D);
  //fnt->SetBranchAddress("BsvpvDisErr_2D", BsvpvDisErr_2D);
  //fnt->SetBranchAddress("Bujeta", Bujeta);
  //fnt->SetBranchAddress("Btrk1Eta", Btrk1Eta);
  //fnt->SetBranchAddress("Btrk2Eta", Btrk2Eta);
  //fnt->SetBranchAddress("Btrk1Dxy1", Btrk1Dxy1);
  //fnt->SetBranchAddress("Btrk2Dxy1", Btrk2Dxy1);
  //fnt->SetBranchAddress("Btrk1DxyError1", Btrk1DxyError1);
  //fnt->SetBranchAddress("Btrk2DxyError1", Btrk2DxyError1);
  //fnt->SetBranchAddress("Btrk1Chi2ndf", Btrk1Chi2ndf);
  //fnt->SetBranchAddress("Btrk2Chi2ndf", Btrk2Chi2ndf);
  //fnt->SetBranchAddress("Btrk1nStripLayer", Btrk1nStripLayer);
  //fnt->SetBranchAddress("Btrk2nStripLayer", Btrk2nStripLayer);
  //fnt->SetBranchAddress("Btrk1nPixelLayer", Btrk1nPixelLayer);
  //fnt->SetBranchAddress("Btrk2nPixelLayer", Btrk2nPixelLayer);
  //fnt->SetBranchAddress("Bmu1TMOneStationTight", Bmu1TMOneStationTight);
  //fnt->SetBranchAddress("Bmu1SoftMuID", Bmu1SoftMuID);
  //fnt->SetBranchAddress("Bmu1isAcc", Bmu1isAcc);
  //fnt->SetBranchAddress("Bmu1isTriggered", Bmu1isTriggered);
  //fnt->SetBranchAddress("Bmu1InPixelLayer", Bmu1InPixelLayer);
  //fnt->SetBranchAddress("Bmu1InStripLayer", Bmu1InStripLayer);
  //fnt->SetBranchAddress("Bmu1dxyPV", Bmu1dxyPV);
  //fnt->SetBranchAddress("Bmu1dzPV", Bmu1dzPV);
  //fnt->SetBranchAddress("Bmu1isGlobalMuon", Bmu1isGlobalMuon);
  //fnt->SetBranchAddress("Bmu1eta", Bmu1eta);
  //fnt->SetBranchAddress("Bmu1pt", Bmu1pt);
  //fnt->SetBranchAddress("Bmu2TMOneStationTight", Bmu2TMOneStationTight);
  //fnt->SetBranchAddress("Bmu2SoftMuID", Bmu2SoftMuID);
  //fnt->SetBranchAddress("Bmu2isAcc", Bmu2isAcc);
  //fnt->SetBranchAddress("Bmu2isTriggered", Bmu2isTriggered);
  //fnt->SetBranchAddress("Bmu2InPixelLayer", Bmu2InPixelLayer);
  //fnt->SetBranchAddress("Bmu2InStripLayer", Bmu2InStripLayer);
  //fnt->SetBranchAddress("Bmu2dxyPV", Bmu2dxyPV);
  //fnt->SetBranchAddress("Bmu2dzPV", Bmu2dzPV);
  //fnt->SetBranchAddress("Bmu2isGlobalMuon", Bmu2isGlobalMuon);
  //fnt->SetBranchAddress("Bmu2eta", Bmu2eta);
  //fnt->SetBranchAddress("Bmu2pt", Bmu2pt);
  //fnt->SetBranchAddress("Bmumumass", Bmumumass);
  //fnt->SetBranchAddress("Btrk1highPurity", Btrk1highPurity);
  //fnt->SetBranchAddress("Btrk1PixelHit", Btrk1PixelHit);
  //fnt->SetBranchAddress("Btrk1StripHit", Btrk1StripHit);
  //fnt->SetBranchAddress("Btrk1PtErr", Btrk1PtErr);
  //fnt->SetBranchAddress("Btrk2highPurity", Btrk2highPurity);
  //fnt->SetBranchAddress("Btrk2PixelHit", Btrk2PixelHit);
  //fnt->SetBranchAddress("Btrk2StripHit", Btrk2StripHit);
  //fnt->SetBranchAddress("Btrk2PtErr", Btrk2PtErr);
  //fnt->SetBranchAddress("By", By);
  fnt->SetBranchAddress("By", &By);
  //fnt->SetBranchAddress("Btktkmass", Btktkmass);
}

int findBptbin(float Bpt)
{
  int n;
  for(int i=0;i<nBinsFine;i++)
    {
      if(Bpt>=5+45.0/nBinsFine*i && Bpt<5+45.0/nBinsFine*(i+1))
	{
	  n=i;
	  break;
	}
    }
  return n;
}

int findBybin(float By)
{
  int n;
  for(int i=0;i<nBinsYFine;i++)
    {
      if(By>=-2.4+4.8/nBinsYFine*i && By<-2.4+4.8/nBinsYFine*(i+1))
	{
	  n=i;
	  break;
	}
    }
  return n;
}
