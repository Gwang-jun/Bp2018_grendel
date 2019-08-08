#include "TTree.h"
#define MAX_XB       20000
#define MAX_GEN      6000
//pp
int     pBeamScrapingFilter;
int     pPAprimaryVertexFilter;
int     HLT_HIL1DoubleMu0ForPPRef_v1;
//pbpb
int     pclusterCompatibilityFilter;
int     pprimaryVertexFilter;
int     phfCoincFilter2Th4;
int     HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1;
float     Ncoll;
//common
float   pthatweight;
int     Gsize;
float   Gpt[MAX_GEN];
float   Gy[MAX_GEN];
int   GpdgId[MAX_GEN];   
int   GisSignal[MAX_GEN];
int     Bsize;
float   Bgen[MAX_XB];
float   Bgenpt[MAX_XB];
float   PVz;
float   Bpt[MAX_XB];
float   By[MAX_XB];
float   Bmumumass[MAX_XB];
float   Bmass[MAX_XB];
float   Bmu1pt[MAX_XB];
float   Bmu2pt[MAX_XB];
float   Bmu1p[MAX_XB];
float   Bmu2p[MAX_XB];
float   Bmu1eta[MAX_XB];
float   Bmu2eta[MAX_XB];
bool    Bmu1TMOneStationTight[MAX_XB];
bool    Bmu2TMOneStationTight[MAX_XB];
int     Bmu1InPixelLayer[MAX_XB];
int     Bmu2InPixelLayer[MAX_XB];
int     Bmu1InStripLayer[MAX_XB];
int     Bmu2InStripLayer[MAX_XB];
float   Bmu1dxyPV[MAX_XB];
float   Bmu2dxyPV[MAX_XB];
float   Bmu1dzPV[MAX_XB];
float   Bmu2dzPV[MAX_XB];
bool    Bmu1isGlobalMuon[MAX_XB];
bool    Bmu2isGlobalMuon[MAX_XB];
float   Bmu1TrgMatchFilterE[MAX_XB];
float   Bmu2TrgMatchFilterE[MAX_XB];
bool    Btrk1highPurity[MAX_XB];
bool    Btrk2highPurity[MAX_XB];
float   Btrk1Eta[MAX_XB];
float   Btrk1Pt[MAX_XB];
float   Btrk1PtErr[MAX_XB];
float   Bchi2cl[MAX_XB];
float   BsvpvDistance[MAX_XB];
float   BsvpvDisErr[MAX_XB];
float   Bdtheta[MAX_XB];
bool   Bmu1SoftMuID[MAX_XB];
bool   Bmu2SoftMuID[MAX_XB];
bool   Bmu1isAcc[MAX_XB];
bool   Bmu2isAcc[MAX_XB];
bool   Bmu1isTriggered[MAX_XB];
bool   Bmu2isTriggered[MAX_XB];
float   Btrk1PixelHit[MAX_XB];
float   Btrk1StripHit[MAX_XB];
float   Btrk1Chi2ndf[MAX_XB];
float   Btrk1nStripLayer[MAX_XB];
float   Btrk1nPixelLayer[MAX_XB];
float   Btrk1Dxy1[MAX_XB];
float   Btrk1DxyError1[MAX_XB];
float   Btrk1Dz1[MAX_XB];
float   Btrk1DzError1[MAX_XB];

void setAddressTree(TTree* ntKp, TTree* ntHlt, TTree* ntSkim, TTree* ntHi, TTree* ntGen, bool ispp){
	if(ispp){
		ntSkim->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter);
		ntSkim->SetBranchAddress("pPAprimaryVertexFilter",&pPAprimaryVertexFilter);
		ntHlt->SetBranchAddress("HLT_HIL1DoubleMu0ForPPRef_v1",&HLT_HIL1DoubleMu0ForPPRef_v1);
	}
	else{
		ntSkim->SetBranchAddress("pclusterCompatibilityFilter",&pclusterCompatibilityFilter);
		ntSkim->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter);
		ntSkim->SetBranchAddress("phfCoincFilter2Th4",&phfCoincFilter2Th4);
		ntHlt->SetBranchAddress("HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1",&HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1);
		ntHi->SetBranchAddress("Ncoll",&Ncoll);
	}
	ntHi->SetBranchAddress("pthatweight",&pthatweight);
	ntGen->SetBranchAddress("Gsize",&Gsize);
	ntGen->SetBranchAddress("Gpt",Gpt);
	ntGen->SetBranchAddress("Gy",Gy);
	ntGen->SetBranchAddress("GpdgId",GpdgId);
	ntGen->SetBranchAddress("GisSignal",GisSignal);
	ntKp->SetBranchAddress("Bsize",&Bsize);
	ntKp->SetBranchAddress("Bgen",Bgen);
	ntKp->SetBranchAddress("Bgenpt",Bgenpt);
	ntKp->SetBranchAddress("PVz",&PVz);
	ntKp->SetBranchAddress("Bpt",Bpt);
	ntKp->SetBranchAddress("By",By);
	ntKp->SetBranchAddress("Bmumumass",Bmumumass);
	ntKp->SetBranchAddress("Bmass",Bmass);
	ntKp->SetBranchAddress("Bmu1pt",Bmu1pt);
	ntKp->SetBranchAddress("Bmu2pt",Bmu2pt);
	ntKp->SetBranchAddress("Bmu1p",Bmu1p);
	ntKp->SetBranchAddress("Bmu1p",Bmu1p);
	ntKp->SetBranchAddress("Bmu1eta",Bmu1eta);
	ntKp->SetBranchAddress("Bmu2eta",Bmu2eta);
	ntKp->SetBranchAddress("Bmu1TMOneStationTight",Bmu1TMOneStationTight);
	ntKp->SetBranchAddress("Bmu2TMOneStationTight",Bmu2TMOneStationTight);
	ntKp->SetBranchAddress("Bmu1InPixelLayer",Bmu1InPixelLayer);
	ntKp->SetBranchAddress("Bmu2InPixelLayer",Bmu2InPixelLayer);
	ntKp->SetBranchAddress("Bmu1InStripLayer",Bmu1InStripLayer);
	ntKp->SetBranchAddress("Bmu2InStripLayer",Bmu2InStripLayer);
	ntKp->SetBranchAddress("Bmu1dxyPV",Bmu1dxyPV);
	ntKp->SetBranchAddress("Bmu2dxyPV",Bmu2dxyPV);
	ntKp->SetBranchAddress("Bmu1dzPV",Bmu1dzPV);
	ntKp->SetBranchAddress("Bmu2dzPV",Bmu2dzPV);
	ntKp->SetBranchAddress("Bmu1isGlobalMuon",Bmu1isGlobalMuon);
	ntKp->SetBranchAddress("Bmu2isGlobalMuon",Bmu2isGlobalMuon);
	ntKp->SetBranchAddress("Bmu1TrgMatchFilterE",Bmu1TrgMatchFilterE);
	ntKp->SetBranchAddress("Bmu2TrgMatchFilterE",Bmu2TrgMatchFilterE);
	ntKp->SetBranchAddress("Btrk1highPurity",Btrk1highPurity);
	ntKp->SetBranchAddress("Btrk2highPurity",Btrk2highPurity);
	ntKp->SetBranchAddress("Btrk1Eta",Btrk1Eta);
	ntKp->SetBranchAddress("Btrk1Pt",Btrk1Pt);
	ntKp->SetBranchAddress("Btrk1PtErr",Btrk1Pt);
	ntKp->SetBranchAddress("Bchi2cl",Bchi2cl);
	ntKp->SetBranchAddress("BsvpvDistance",BsvpvDistance);
	ntKp->SetBranchAddress("BsvpvDisErr",BsvpvDisErr);
	ntKp->SetBranchAddress("Bdtheta",Bdtheta);
	ntKp->SetBranchAddress("Bmu1SoftMuID",Bmu1SoftMuID);
	ntKp->SetBranchAddress("Bmu2SoftMuID",Bmu2SoftMuID);
	ntKp->SetBranchAddress("Bmu1isAcc",Bmu1isAcc);
	ntKp->SetBranchAddress("Bmu2isAcc",Bmu2isAcc);
	ntKp->SetBranchAddress("Bmu1isTriggered",Bmu1isTriggered);
	ntKp->SetBranchAddress("Bmu2isTriggered",Bmu2isTriggered);
	ntKp->SetBranchAddress("Btrk1PixelHit",Btrk1PixelHit);
	ntKp->SetBranchAddress("Btrk1StripHit",Btrk1StripHit);
	ntKp->SetBranchAddress("Btrk1Chi2ndf",Btrk1Chi2ndf);
	ntKp->SetBranchAddress("Btrk1nStripLayer",Btrk1nStripLayer);
	ntKp->SetBranchAddress("Btrk1nPixelLayer",Btrk1nPixelLayer);
	ntKp->SetBranchAddress("Btrk1Dxy1",Btrk1Dxy1);
	ntKp->SetBranchAddress("Btrk1DxyError1",Btrk1DxyError1);
	ntKp->SetBranchAddress("Btrk1Dz1",Btrk1Dz1);
	ntKp->SetBranchAddress("Btrk1DzError1",Btrk1DzError1);

}
