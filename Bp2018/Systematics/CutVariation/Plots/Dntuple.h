using namespace std;

#ifndef _DNTUPLE_H_
#define _DNTUPLE_H_
#include "format.h"
#include "TMVA_BDTG_pp_pt_1p0_2p0.class.C"
#include "uti.h"
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <vector>
#include <iostream>
class DntupleBranches
{//{{{
	public:
		//EvtInfo
		int     RunNo;
		int     EvtNo;
		int     LumiNo;
		int     Dsize;
		float   PVx;
		float   PVy;
		float   PVz;
		float   PVxE;
		float   PVyE;
		float   PVzE;
		float   PVnchi2;
		float   PVchi2;
		float   BSx;
		float   BSy;
		float   BSz;
		float   BSxErr;
		float   BSyErr;
		float   BSzErr;
		float   BSdxdz;
		float   BSdydz;
		float   BSdxdzErr;
		float   BSdydzErr;
		float   BSWidthX;
		float   BSWidthXErr;
		float   BSWidthY;
		float   BSWidthYErr;
		//DInfo
		int     Dindex[MAX_XB];
		int     Dtype[MAX_XB];
		float   Dmass[MAX_XB];
		float   Dpt[MAX_XB];
		float   Deta[MAX_XB];
		float   Dphi[MAX_XB];
		float   Dy[MAX_XB];
		float   DvtxX[MAX_XB];
		float   DvtxY[MAX_XB];
		float   DvtxZ[MAX_XB];
		float   Dd0[MAX_XB];
		float   Dd0Err[MAX_XB];
		float   Ddxyz[MAX_XB];
		float   DdxyzErr[MAX_XB];
		float   Dchi2ndf[MAX_XB];
		float   Dchi2cl[MAX_XB];
		float   Ddtheta[MAX_XB];
		float   Dlxy[MAX_XB];
		float   Dalpha[MAX_XB];
		float   DsvpvDistance[MAX_XB];
		float   DsvpvDisErr[MAX_XB];
		float   DsvpvDistance_2D[MAX_XB];
		float   DsvpvDisErr_2D[MAX_XB];
		float   Ddca[MAX_XB];
		float   DlxyBS[MAX_XB];
		float   DlxyBSErr[MAX_XB];
		float   DMaxDoca[MAX_XB];
		float   DMaxTkPt[MAX_XB];
		float   DMinTkPt[MAX_XB];

		//DInfo.trkInfo
		int     Dtrk1Idx[MAX_XB];
		int     Dtrk2Idx[MAX_XB];
		int     Dtrk3Idx[MAX_XB];
		int     Dtrk4Idx[MAX_XB];
		float   Dtrk1Pt[MAX_XB];
		float   Dtrk2Pt[MAX_XB];
		float   Dtrk3Pt[MAX_XB];
		float   Dtrk4Pt[MAX_XB];
		float   Dtrk1Eta[MAX_XB];
		float   Dtrk2Eta[MAX_XB];
		float   Dtrk3Eta[MAX_XB];
		float   Dtrk4Eta[MAX_XB];
		float   Dtrk1Phi[MAX_XB];
		float   Dtrk2Phi[MAX_XB];
		float   Dtrk3Phi[MAX_XB];
		float   Dtrk4Phi[MAX_XB];
		float   Dtrk1P[MAX_XB];
		float   Dtrk2P[MAX_XB];
		float   Dtrk3P[MAX_XB];
		float   Dtrk4P[MAX_XB];
		float   Dtrk1PtErr[MAX_XB];
		float   Dtrk2PtErr[MAX_XB];
		float   Dtrk3PtErr[MAX_XB];
		float   Dtrk4PtErr[MAX_XB];
		float   Dtrk1EtaErr[MAX_XB];
		float   Dtrk2EtaErr[MAX_XB];
		float   Dtrk3EtaErr[MAX_XB];
		float   Dtrk4EtaErr[MAX_XB];
		float   Dtrk1PhiErr[MAX_XB];
		float   Dtrk2PhiErr[MAX_XB];
		float   Dtrk3PhiErr[MAX_XB];
		float   Dtrk4PhiErr[MAX_XB];
		float   Dtrk1Y[MAX_XB];
		float   Dtrk2Y[MAX_XB];
		float   Dtrk3Y[MAX_XB];
		float   Dtrk4Y[MAX_XB];
		float   Dtrk1Dz[MAX_XB];
		float   Dtrk2Dz[MAX_XB];
		float   Dtrk3Dz[MAX_XB];
		float   Dtrk4Dz[MAX_XB];
		float   Dtrk1Dxy[MAX_XB];
		float   Dtrk2Dxy[MAX_XB];
		float   Dtrk3Dxy[MAX_XB];
		float   Dtrk4Dxy[MAX_XB];
		float   Dtrk1D0[MAX_XB];
		float   Dtrk2D0[MAX_XB];
		float   Dtrk3D0[MAX_XB];
		float   Dtrk4D0[MAX_XB];
		float   Dtrk1D0Err[MAX_XB];
		float   Dtrk2D0Err[MAX_XB];
		float   Dtrk3D0Err[MAX_XB];
		float   Dtrk4D0Err[MAX_XB];
		float   Dtrk1PixelHit[MAX_XB];
		float   Dtrk2PixelHit[MAX_XB];
		float   Dtrk3PixelHit[MAX_XB];
		float   Dtrk4PixelHit[MAX_XB];
		float   Dtrk1StripHit[MAX_XB];
		float   Dtrk2StripHit[MAX_XB];
		float   Dtrk3StripHit[MAX_XB];
		float   Dtrk4StripHit[MAX_XB];
		float   Dtrk1nStripLayer[MAX_XB];
		float   Dtrk2nStripLayer[MAX_XB];
		float   Dtrk3nStripLayer[MAX_XB];
		float   Dtrk4nStripLayer[MAX_XB];
		float   Dtrk1nPixelLayer[MAX_XB];
		float   Dtrk2nPixelLayer[MAX_XB];
		float   Dtrk3nPixelLayer[MAX_XB];
		float   Dtrk4nPixelLayer[MAX_XB];
		float   Dtrk1Chi2ndf[MAX_XB];
		float   Dtrk2Chi2ndf[MAX_XB];
		float   Dtrk3Chi2ndf[MAX_XB];
		float   Dtrk4Chi2ndf[MAX_XB];
		float   Dtrk1MassHypo[MAX_XB];
		float   Dtrk2MassHypo[MAX_XB];
		float   Dtrk3MassHypo[MAX_XB];
		float   Dtrk4MassHypo[MAX_XB];
		float   Dtrk1MVAVal[MAX_XB];
		float   Dtrk2MVAVal[MAX_XB];
		float   Dtrk3MVAVal[MAX_XB];
		float   Dtrk4MVAVal[MAX_XB];
		int     Dtrk1Algo[MAX_XB];
		int     Dtrk2Algo[MAX_XB];
		int     Dtrk3Algo[MAX_XB];
		int     Dtrk4Algo[MAX_XB];
		int     Dtrk1originalAlgo[MAX_XB];
		int     Dtrk2originalAlgo[MAX_XB];
		int     Dtrk3originalAlgo[MAX_XB];
		int     Dtrk4originalAlgo[MAX_XB];
		bool    Dtrk1highPurity[MAX_XB];
		bool    Dtrk2highPurity[MAX_XB];
		bool    Dtrk3highPurity[MAX_XB];
		bool    Dtrk4highPurity[MAX_XB];
		int     Dtrk1Quality[MAX_XB];
		int     Dtrk2Quality[MAX_XB];
		int     Dtrk3Quality[MAX_XB];
		int     Dtrk4Quality[MAX_XB];
		float   Dtrk1dedx[MAX_XB];
		float   Dtrk2dedx[MAX_XB];
		float   Dtrk3dedx[MAX_XB];
		float   Dtrk4dedx[MAX_XB];
		float   Dtrk1thetastar[MAX_XB];
		float   Dtrk2thetastar[MAX_XB];
		float   Dtrk3thetastar[MAX_XB];
		float   Dtrk4thetastar[MAX_XB];
		float   BDTG_pt_1_2[MAX_XB];

		//DInfo.tktkResInfo
		float   DtktkResmass[MAX_XB];
		float   DtktkRespt[MAX_XB];
		float   DtktkReseta[MAX_XB];
		float   DtktkResphi[MAX_XB];
		float   DtktkRes_chi2ndf[MAX_XB];
		float   DtktkRes_chi2cl[MAX_XB];
		float   DtktkRes_alpha[MAX_XB];
		float   DtktkRes_alphaToSV[MAX_XB];
		float   DtktkRes_svpvDistance[MAX_XB];
		float   DtktkRes_svpvDisErr[MAX_XB];
		float   DtktkRes_svpvDistanceToSV[MAX_XB];
		float   DtktkRes_svpvDisErrToSV[MAX_XB];
		float   DtktkRes_dca[MAX_XB];
		float   DtktkRes_dcaToSV[MAX_XB];
		float   DtktkRes_lxyBS[MAX_XB];
		float   DtktkRes_lxyBSErr[MAX_XB];
		float   DtktkRes_angleToTrk1[MAX_XB];
		float   DtktkRes_ptAsymToTrk1[MAX_XB];
		float   DtktkRes_unfitter_ptAsymToTrk1[MAX_XB];

		float   DRestrk1Pt[MAX_XB];
		float   DRestrk1Eta[MAX_XB];
		float   DRestrk1Phi[MAX_XB];
		float   DRestrk1Y[MAX_XB];
		float   DRestrk1P[MAX_XB];
		float   DRestrk1Dz[MAX_XB];
		float   DRestrk1Dxy[MAX_XB];
		float   DRestrk1D0[MAX_XB];
		float   DRestrk1D0Err[MAX_XB];
		int     DRestrk1originalAlgo[MAX_XB];
		float   DRestrk2Pt[MAX_XB];
		float   DRestrk2Eta[MAX_XB];
		float   DRestrk2Phi[MAX_XB];
		float   DRestrk2Y[MAX_XB];
		float   DRestrk2P[MAX_XB];
		float   DRestrk2Dz[MAX_XB];
		float   DRestrk2Dxy[MAX_XB];
		float   DRestrk2D0[MAX_XB];
		float   DRestrk2D0Err[MAX_XB];
		int     DRestrk2originalAlgo[MAX_XB];
		float   DRestrk3Pt[MAX_XB];
		float   DRestrk3Eta[MAX_XB];
		float   DRestrk3Phi[MAX_XB];
		float   DRestrk3Y[MAX_XB];
		float   DRestrk3P[MAX_XB];
		float   DRestrk3Dz[MAX_XB];
		float   DRestrk3Dxy[MAX_XB];
		float   DRestrk3D0[MAX_XB];
		float   DRestrk3D0Err[MAX_XB];
		int     DRestrk3originalAlgo[MAX_XB];
		float   DRestrk4Pt[MAX_XB];
		float   DRestrk4Eta[MAX_XB];
		float   DRestrk4Phi[MAX_XB];
		float   DRestrk4P[MAX_XB];
		float   DRestrk4Y[MAX_XB];
		float   DRestrk4Dz[MAX_XB];
		float   DRestrk4Dxy[MAX_XB];
		float   DRestrk4D0[MAX_XB];
		float   DRestrk4D0Err[MAX_XB];
		int     DRestrk4originalAlgo[MAX_XB];
		float   DRestrk1PtErr[MAX_XB];
		float   DRestrk2PtErr[MAX_XB];
		float   DRestrk3PtErr[MAX_XB];
		float   DRestrk4PtErr[MAX_XB];
		bool    DRestrk1highPurity[MAX_XB];
		bool    DRestrk2highPurity[MAX_XB];
		bool    DRestrk3highPurity[MAX_XB];
		bool    DRestrk4highPurity[MAX_XB];
		int     DRestrk1Quality[MAX_XB];
		int     DRestrk2Quality[MAX_XB];
		int     DRestrk3Quality[MAX_XB];
		int     DRestrk4Quality[MAX_XB];
		float   DRestrk1dedx[MAX_XB];
		float   DRestrk2dedx[MAX_XB];
		float   DRestrk3dedx[MAX_XB];
		float   DRestrk4dedx[MAX_XB];
		float   DRestrk1thetastar[MAX_XB];
		float   DRestrk2thetastar[MAX_XB];
		float   DRestrk3thetastar[MAX_XB];
		float   DRestrk4thetastar[MAX_XB];
		//DInfo.genInfo
		float   Dgen[MAX_XB];
		int     DgennDa[MAX_XB];
		int     DgenIndex[MAX_XB];
		float   Dgenpt[MAX_XB];
		float   Dgeneta[MAX_XB];
		float   Dgenphi[MAX_XB];
		float   Dgeny[MAX_XB];
		int     DgencollisionId[MAX_XB];
		float   DgenBAncestorpt[MAX_XB];
		int     DgenBAncestorpdgId[MAX_XB];
		float   DgenprodvtxX[MAX_XB];
		float   DgenprodvtxY[MAX_XB];
		float   DgenprodvtxZ[MAX_XB];
		float   DgendecayvtxX[MAX_XB];
		float   DgendecayvtxY[MAX_XB];
		float   DgendecayvtxZ[MAX_XB];
		int     DgenfromgenPV[MAX_XB];

		void buildDBranch(TTree* dnt, bool D0kpimode=false, bool detailMode=true)
		{
			//EvtInfo
			dnt->Branch("RunNo",&RunNo);
			dnt->Branch("EvtNo",&EvtNo);
			dnt->Branch("LumiNo",&LumiNo);
			dnt->Branch("Dsize",&Dsize);
			dnt->Branch("PVx",&PVx);
			dnt->Branch("PVy",&PVy);
			dnt->Branch("PVz",&PVz);
			dnt->Branch("PVnchi2",&PVnchi2);
			dnt->Branch("BSx",&BSx);
			dnt->Branch("BSy",&BSy);
			dnt->Branch("BSz",&BSz);
			if(detailMode)
			{
				dnt->Branch("PVxE",&PVxE);
				dnt->Branch("PVyE",&PVyE);
				dnt->Branch("PVzE",&PVzE);
				dnt->Branch("BSxErr",&BSxErr);
				dnt->Branch("BSyErr",&BSyErr);
				dnt->Branch("BSzErr",&BSzErr);
				dnt->Branch("BSdxdz",&BSdxdz);
				dnt->Branch("BSdydz",&BSdydz);
				dnt->Branch("BSdxdzErr",&BSdxdzErr);
				dnt->Branch("BSdydzErr",&BSdydzErr);
				dnt->Branch("BSWidthX",&BSWidthX);
				dnt->Branch("BSWidthXErr",&BSWidthXErr);
				dnt->Branch("BSWidthY",&BSWidthY);
				dnt->Branch("BSWidthYErr",&BSWidthYErr);        
			}
			//DInfo
			dnt->Branch("Dindex",Dindex,"Dindex[Dsize]/I");
			dnt->Branch("Dtype",Dtype,"Dtype[Dsize]/I");
			dnt->Branch("Dmass",Dmass,"Dmass[Dsize]/F");
			dnt->Branch("Dpt",Dpt,"Dpt[Dsize]/F");
			dnt->Branch("Deta",Deta,"Deta[Dsize]/F");
			dnt->Branch("Dphi",Dphi,"Dphi[Dsize]/F");
			dnt->Branch("Dy",Dy,"Dy[Dsize]/F");
			dnt->Branch("DvtxX",DvtxX,"DvtxX[Dsize]/F");
			dnt->Branch("DvtxY",DvtxY,"DvtxY[Dsize]/F");
			dnt->Branch("DvtxZ",DvtxZ,"DvtxZ[Dsize]/F");
			dnt->Branch("Dd0",Dd0,"Dd0[Dsize]/F");
			dnt->Branch("Dd0Err",Dd0Err,"Dd0Err[Dsize]/F");
			dnt->Branch("Ddxyz",Ddxyz,"Ddxyz[Dsize]/F");
			dnt->Branch("DdxyzErr",DdxyzErr,"DdxyzErr[Dsize]/F");
			dnt->Branch("Dchi2ndf",Dchi2ndf,"Dchi2ndf[Dsize]/F");
			dnt->Branch("Dchi2cl",Dchi2cl,"Dchi2cl[Dsize]/F");
			dnt->Branch("Ddtheta",Ddtheta,"Ddtheta[Dsize]/F");
			dnt->Branch("Dlxy",Dlxy,"Dlxy[Dsize]/F");
			dnt->Branch("Dalpha",Dalpha,"Dalpha[Dsize]/F");
			dnt->Branch("DsvpvDistance",DsvpvDistance,"DsvpvDistance[Dsize]/F");
			dnt->Branch("DsvpvDisErr",DsvpvDisErr,"DsvpvDisErr[Dsize]/F");
			dnt->Branch("DsvpvDistance_2D",DsvpvDistance_2D,"DsvpvDistance_2D[Dsize]/F");
			dnt->Branch("DsvpvDisErr_2D",DsvpvDisErr_2D,"DsvpvDisErr_2D[Dsize]/F");
			dnt->Branch("Ddca",Ddca,"Ddca[Dsize]/F");
			dnt->Branch("DlxyBS",DlxyBS,"DlxyBS[Dsize]/F");
			dnt->Branch("DlxyBSErr",DlxyBSErr,"DlxyBSErr[Dsize]/F");
			dnt->Branch("DMaxDoca",DMaxDoca,"DMaxDoca[Dsize]/F");
			dnt->Branch("DMaxTkPt",DMaxTkPt,"DMaxTkPt[Dsize]/F");
			dnt->Branch("DMinTkPt",DMinTkPt,"DMinTkPt[Dsize]/F");

			//DInfo.trkInfo
			dnt->Branch("Dtrk1Pt",Dtrk1Pt,"Dtrk1Pt[Dsize]/F");
			dnt->Branch("Dtrk2Pt",Dtrk2Pt,"Dtrk2Pt[Dsize]/F");
			dnt->Branch("Dtrk1Eta",Dtrk1Eta,"Dtrk1Eta[Dsize]/F");
			dnt->Branch("Dtrk2Eta",Dtrk2Eta,"Dtrk2Eta[Dsize]/F");
			dnt->Branch("Dtrk1Phi",Dtrk1Phi,"Dtrk1Phi[Dsize]/F");
			dnt->Branch("Dtrk2Phi",Dtrk2Phi,"Dtrk2Phi[Dsize]/F");
			dnt->Branch("Dtrk1P",Dtrk1P,"Dtrk1P[Dsize]/F");
			dnt->Branch("Dtrk2P",Dtrk2P,"Dtrk2P[Dsize]/F");
			dnt->Branch("Dtrk1PtErr",Dtrk1PtErr,"Dtrk1PtErr[Dsize]/F");
			dnt->Branch("Dtrk2PtErr",Dtrk2PtErr,"Dtrk2PtErr[Dsize]/F");
			dnt->Branch("Dtrk1Dz",Dtrk1Dz,"Dtrk1Dz[Dsize]/F");
			dnt->Branch("Dtrk2Dz",Dtrk2Dz,"Dtrk2Dz[Dsize]/F");
			dnt->Branch("Dtrk1Dxy",Dtrk1Dxy,"Dtrk1Dxy[Dsize]/F");
			dnt->Branch("Dtrk2Dxy",Dtrk2Dxy,"Dtrk2Dxy[Dsize]/F");
			dnt->Branch("Dtrk1PixelHit",Dtrk1PixelHit,"Dtrk1PixelHit[Dsize]/F");
			dnt->Branch("Dtrk2PixelHit",Dtrk2PixelHit,"Dtrk2PixelHit[Dsize]/F");
			dnt->Branch("Dtrk1StripHit",Dtrk1StripHit,"Dtrk1StripHit[Dsize]/F");
			dnt->Branch("Dtrk2StripHit",Dtrk2StripHit,"Dtrk2StripHit[Dsize]/F");
			dnt->Branch("Dtrk1nStripLayer",Dtrk1nStripLayer,"Dtrk1nStripLayer[Dsize]/F");
			dnt->Branch("Dtrk2nStripLayer",Dtrk2nStripLayer,"Dtrk2nStripLayer[Dsize]/F");
			dnt->Branch("Dtrk1nPixelLayer",Dtrk1nPixelLayer,"Dtrk1nPixelLayer[Dsize]/F");
			dnt->Branch("Dtrk2nPixelLayer",Dtrk2nPixelLayer,"Dtrk2nPixelLayer[Dsize]/F");
			dnt->Branch("Dtrk1Chi2ndf",Dtrk1Chi2ndf,"Dtrk1Chi2ndf[Dsize]/F");
			dnt->Branch("Dtrk2Chi2ndf",Dtrk2Chi2ndf,"Dtrk2Chi2ndf[Dsize]/F");
			dnt->Branch("Dtrk1MassHypo",Dtrk1MassHypo,"Dtrk1MassHypo[Dsize]/F");
			dnt->Branch("Dtrk2MassHypo",Dtrk2MassHypo,"Dtrk2MassHypo[Dsize]/F");
			dnt->Branch("Dtrk1Algo",Dtrk1Algo,"Dtrk1Algo[Dsize]/I");
			dnt->Branch("Dtrk2Algo",Dtrk2Algo,"Dtrk2Algo[Dsize]/I");
			dnt->Branch("Dtrk1originalAlgo",Dtrk1originalAlgo,"Dtrk1originalAlgo[Dsize]/I");
			dnt->Branch("Dtrk2originalAlgo",Dtrk2originalAlgo,"Dtrk2originalAlgo[Dsize]/I");
			dnt->Branch("Dtrk1highPurity",Dtrk1highPurity,"Dtrk1highPurity[Dsize]/O");
			dnt->Branch("Dtrk2highPurity",Dtrk2highPurity,"Dtrk2highPurity[Dsize]/O");
			dnt->Branch("Dtrk1dedx",Dtrk1dedx,"Dtrk1dedx[Dsize]/F");
			dnt->Branch("Dtrk2dedx",Dtrk2dedx,"Dtrk2dedx[Dsize]/F");
			dnt->Branch("Dtrk1thetastar",Dtrk1thetastar,"Dtrk1thetastar[Dsize]/F");
			dnt->Branch("Dtrk2thetastar",Dtrk2thetastar,"Dtrk2thetastar[Dsize]/F");
			dnt->Branch("BDTG_pt_1_2",BDTG_pt_1_2,"BDTG_pt_1_2[Dsize]/F");
			if(!D0kpimode)
			{
				dnt->Branch("Dtrk3Pt",Dtrk3Pt,"Dtrk3Pt[Dsize]/F");
				dnt->Branch("Dtrk4Pt",Dtrk4Pt,"Dtrk4Pt[Dsize]/F");
				dnt->Branch("Dtrk3Eta",Dtrk3Eta,"Dtrk3Eta[Dsize]/F");
				dnt->Branch("Dtrk4Eta",Dtrk4Eta,"Dtrk4Eta[Dsize]/F");
				dnt->Branch("Dtrk3Phi",Dtrk3Phi,"Dtrk3Phi[Dsize]/F");
				dnt->Branch("Dtrk4Phi",Dtrk4Phi,"Dtrk4Phi[Dsize]/F");
				dnt->Branch("Dtrk3P",Dtrk3P,"Dtrk3P[Dsize]/F");
				dnt->Branch("Dtrk4P",Dtrk4P,"Dtrk4P[Dsize]/F");
				dnt->Branch("Dtrk3PtErr",Dtrk3PtErr,"Dtrk3PtErr[Dsize]/F");
				dnt->Branch("Dtrk4PtErr",Dtrk4PtErr,"Dtrk4PtErr[Dsize]/F");
				dnt->Branch("Dtrk3Dz",Dtrk3Dz,"Dtrk3Dz[Dsize]/F");
				dnt->Branch("Dtrk4Dz",Dtrk4Dz,"Dtrk4Dz[Dsize]/F");
				dnt->Branch("Dtrk3Dxy",Dtrk3Dxy,"Dtrk3Dxy[Dsize]/F");
				dnt->Branch("Dtrk4Dxy",Dtrk4Dxy,"Dtrk4Dxy[Dsize]/F");
				dnt->Branch("Dtrk3PixelHit",Dtrk3PixelHit,"Dtrk3PixelHit[Dsize]/F");
				dnt->Branch("Dtrk4PixelHit",Dtrk4PixelHit,"Dtrk4PixelHit[Dsize]/F");
				dnt->Branch("Dtrk3StripHit",Dtrk3StripHit,"Dtrk3StripHit[Dsize]/F");
				dnt->Branch("Dtrk4StripHit",Dtrk4StripHit,"Dtrk4StripHit[Dsize]/F");
				dnt->Branch("Dtrk3nStripLayer",Dtrk3nStripLayer,"Dtrk3nStripLayer[Dsize]/F");
				dnt->Branch("Dtrk4nStripLayer",Dtrk4nStripLayer,"Dtrk4nStripLayer[Dsize]/F");
				dnt->Branch("Dtrk3nPixelLayer",Dtrk3nPixelLayer,"Dtrk3nPixelLayer[Dsize]/F");
				dnt->Branch("Dtrk4nPixelLayer",Dtrk4nPixelLayer,"Dtrk4nPixelLayer[Dsize]/F");
				dnt->Branch("Dtrk3Chi2ndf",Dtrk3Chi2ndf,"Dtrk3Chi2ndf[Dsize]/F");
				dnt->Branch("Dtrk4Chi2ndf",Dtrk4Chi2ndf,"Dtrk4Chi2ndf[Dsize]/F");
				dnt->Branch("Dtrk3MassHypo",Dtrk3MassHypo,"Dtrk3MassHypo[Dsize]/F");
				dnt->Branch("Dtrk4MassHypo",Dtrk4MassHypo,"Dtrk4MassHypo[Dsize]/F");
				dnt->Branch("Dtrk3Algo",Dtrk3Algo,"Dtrk3Algo[Dsize]/I");
				dnt->Branch("Dtrk4Algo",Dtrk4Algo,"Dtrk4Algo[Dsize]/I");
				dnt->Branch("Dtrk3originalAlgo",Dtrk3originalAlgo,"Dtrk3originalAlgo[Dsize]/I");
				dnt->Branch("Dtrk4originalAlgo",Dtrk4originalAlgo,"Dtrk4originalAlgo[Dsize]/I");
				dnt->Branch("Dtrk3highPurity",Dtrk3highPurity,"Dtrk3highPurity[Dsize]/O");
				dnt->Branch("Dtrk4highPurity",Dtrk4highPurity,"Dtrk4highPurity[Dsize]/O");
				dnt->Branch("Dtrk3dedx",Dtrk3dedx,"Dtrk3dedx[Dsize]/F");
				dnt->Branch("Dtrk4dedx",Dtrk4dedx,"Dtrk4dedx[Dsize]/F");
				dnt->Branch("Dtrk3thetastar",Dtrk3thetastar,"Dtrk3thetastar[Dsize]/F");
				dnt->Branch("Dtrk4thetastar",Dtrk4thetastar,"Dtrk4thetastar[Dsize]/F");
			}
			if(detailMode)
			{
				dnt->Branch("Dtrk1Idx",Dtrk1Idx,"Dtrk1Idx[Dsize]/I");
				dnt->Branch("Dtrk2Idx",Dtrk2Idx,"Dtrk2Idx[Dsize]/I");
				dnt->Branch("Dtrk1EtaErr",Dtrk1EtaErr,"Dtrk1EtaErr[Dsize]/F");
				dnt->Branch("Dtrk2EtaErr",Dtrk2EtaErr,"Dtrk2EtaErr[Dsize]/F");
				dnt->Branch("Dtrk1PhiErr",Dtrk1PhiErr,"Dtrk1PhiErr[Dsize]/F");
				dnt->Branch("Dtrk2PhiErr",Dtrk2PhiErr,"Dtrk2PhiErr[Dsize]/F");
				dnt->Branch("Dtrk1Y",Dtrk1Y,"Dtrk1Y[Dsize]/F");
				dnt->Branch("Dtrk2Y",Dtrk2Y,"Dtrk2Y[Dsize]/F");
				dnt->Branch("Dtrk1D0",Dtrk1D0,"Dtrk1D0[Dsize]/F");
				dnt->Branch("Dtrk2D0",Dtrk2D0,"Dtrk2D0[Dsize]/F"); 
				dnt->Branch("Dtrk1D0Err",Dtrk1D0Err,"Dtrk1D0Err[Dsize]/F");
				dnt->Branch("Dtrk2D0Err",Dtrk2D0Err,"Dtrk2D0Err[Dsize]/F"); 
				dnt->Branch("Dtrk1MVAVal",Dtrk1MVAVal,"Dtrk1MVAVal[Dsize]/F");
				dnt->Branch("Dtrk2MVAVal",Dtrk2MVAVal,"Dtrk2MVAVal[Dsize]/F");
				dnt->Branch("Dtrk1Quality",Dtrk1Quality,"Dtrk1Quality[Dsize]/I");
				dnt->Branch("Dtrk2Quality",Dtrk2Quality,"Dtrk2Quality[Dsize]/I");
				if(!D0kpimode)
				{
					dnt->Branch("Dtrk3Idx",Dtrk3Idx,"Dtrk3Idx[Dsize]/I");
					dnt->Branch("Dtrk4Idx",Dtrk4Idx,"Dtrk4Idx[Dsize]/I");
					dnt->Branch("Dtrk3EtaErr",Dtrk3EtaErr,"Dtrk3EtaErr[Dsize]/F");
					dnt->Branch("Dtrk4EtaErr",Dtrk4EtaErr,"Dtrk4EtaErr[Dsize]/F");
					dnt->Branch("Dtrk3PhiErr",Dtrk3PhiErr,"Dtrk3PhiErr[Dsize]/F");
					dnt->Branch("Dtrk4PhiErr",Dtrk4PhiErr,"Dtrk4PhiErr[Dsize]/F");
					dnt->Branch("Dtrk3Y",Dtrk3Y,"Dtrk3Y[Dsize]/F");
					dnt->Branch("Dtrk4Y",Dtrk4Y,"Dtrk4Y[Dsize]/F");
					dnt->Branch("Dtrk3D0",Dtrk3D0,"Dtrk3D0[Dsize]/F");
					dnt->Branch("Dtrk4D0",Dtrk4D0,"Dtrk4D0[Dsize]/F");
					dnt->Branch("Dtrk3D0Err",Dtrk3D0Err,"Dtrk3D0Err[Dsize]/F");
					dnt->Branch("Dtrk4D0Err",Dtrk4D0Err,"Dtrk4D0Err[Dsize]/F");
					dnt->Branch("Dtrk3MVAVal",Dtrk3MVAVal,"Dtrk3MVAVal[Dsize]/F");
					dnt->Branch("Dtrk4MVAVal",Dtrk4MVAVal,"Dtrk4MVAVal[Dsize]/F");
					dnt->Branch("Dtrk3Quality",Dtrk3Quality,"Dtrk3Quality[Dsize]/I");
					dnt->Branch("Dtrk4Quality",Dtrk4Quality,"Dtrk4Quality[Dsize]/I");
				}
			}
			//DInfo.tktkResInfo
			if(!D0kpimode)
			{
				dnt->Branch("DtktkResmass",DtktkResmass,"DtktkResmass[Dsize]/F");
				dnt->Branch("DtktkRespt",DtktkRespt,"DtktkRespt[Dsize]/F");
				dnt->Branch("DtktkReseta",DtktkReseta,"DtktkReseta[Dsize]/F");
				dnt->Branch("DtktkResphi",DtktkResphi,"DtktkResphi[Dsize]/F");
				dnt->Branch("DtktkRes_chi2ndf",DtktkRes_chi2ndf,"DtktkRes_chi2ndf[Dsize]/F");
				dnt->Branch("DtktkRes_chi2cl",DtktkRes_chi2cl,"DtktkRes_chi2cl[Dsize]/F");
				dnt->Branch("DtktkRes_alpha",DtktkRes_alpha,"DtktkRes_alpha[Dsize]/F");
				dnt->Branch("DtktkRes_alphaToSV",DtktkRes_alphaToSV,"DtktkRes_alphaToSV[Dsize]/F");
				dnt->Branch("DtktkRes_svpvDistance",DtktkRes_svpvDistance,"DtktkRes_svpvDistance[Dsize]/F");
				dnt->Branch("DtktkRes_svpvDisErr",DtktkRes_svpvDisErr,"DtktkRes_svpvDisErr[Dsize]/F");
				dnt->Branch("DtktkRes_svpvDistanceToSV",DtktkRes_svpvDistanceToSV,"DtktkRes_svpvDistanceToSV[Dsize]/F");
				dnt->Branch("DtktkRes_svpvDisErrToSV",DtktkRes_svpvDisErrToSV,"DtktkRes_svpvDisErrToSV[Dsize]/F");
				dnt->Branch("DtktkRes_dca",DtktkRes_dca,"DtktkRes_dca[Dsize]/F");
				dnt->Branch("DtktkRes_dcaToSV",DtktkRes_dcaToSV,"DtktkRes_dcaToSV[Dsize]/F");
				dnt->Branch("DtktkRes_lxyBS",DtktkRes_lxyBS,"DtktkRes_lxyBS[Dsize]/F");
				dnt->Branch("DtktkRes_lxyBSErr",DtktkRes_lxyBSErr,"DtktkRes_lxyBSErr[Dsize]/F");
				dnt->Branch("DtktkRes_angleToTrk1",DtktkRes_angleToTrk1,"DtktkRes_angleToTrk1[Dsize]/F");
				dnt->Branch("DtktkRes_ptAsymToTrk1",DtktkRes_ptAsymToTrk1,"DtktkRes_ptAsymToTrk1[Dsize]/F");
				dnt->Branch("DtktkRes_unfitter_ptAsymToTrk1",DtktkRes_unfitter_ptAsymToTrk1,"DtktkRes_unfitter_ptAsymToTrk1[Dsize]/F");

				dnt->Branch("DRestrk1Pt",DRestrk1Pt,"DRestrk1Pt[Dsize]/F");
				dnt->Branch("DRestrk1Eta",DRestrk1Eta,"DRestrk1Eta[Dsize]/F");
				dnt->Branch("DRestrk1Phi",DRestrk1Phi,"DRestrk1Phi[Dsize]/F");
				dnt->Branch("DRestrk1P",DRestrk1P,"DRestrk1P[Dsize]/F");
				dnt->Branch("DRestrk1Dz",DRestrk1Dz,"DRestrk1Dz[Dsize]/F");
				dnt->Branch("DRestrk1Dxy",DRestrk1Dxy,"DRestrk1Dxy[Dsize]/F");
				dnt->Branch("DRestrk1originalAlgo",DRestrk1originalAlgo,"DRestrk1originalAlgo[Dsize]/I");
				dnt->Branch("DRestrk2Pt",DRestrk2Pt,"DRestrk2Pt[Dsize]/F");
				dnt->Branch("DRestrk2Eta",DRestrk2Eta,"DRestrk2Eta[Dsize]/F");
				dnt->Branch("DRestrk2Phi",DRestrk2Phi,"DRestrk2Phi[Dsize]/F");
				dnt->Branch("DRestrk2P",DRestrk2P,"DRestrk2P[Dsize]/F");
				dnt->Branch("DRestrk2Dz",DRestrk2Dz,"DRestrk2Dz[Dsize]/F");
				dnt->Branch("DRestrk2Dxy",DRestrk2Dxy,"DRestrk2Dxy[Dsize]/F");
				dnt->Branch("DRestrk2originalAlgo",DRestrk2originalAlgo,"DRestrk2originalAlgo[Dsize]/I");
				dnt->Branch("DRestrk3Pt",DRestrk3Pt,"DRestrk3Pt[Dsize]/F");
				dnt->Branch("DRestrk3Eta",DRestrk3Eta,"DRestrk3Eta[Dsize]/F");
				dnt->Branch("DRestrk3Phi",DRestrk3Phi,"DRestrk3Phi[Dsize]/F");
				dnt->Branch("DRestrk3P",DRestrk3P,"DRestrk3P[Dsize]/F");
				dnt->Branch("DRestrk3Dz",DRestrk3Dz,"DRestrk3Dz[Dsize]/F");
				dnt->Branch("DRestrk3Dxy",DRestrk3Dxy,"DRestrk3Dxy[Dsize]/F");
				dnt->Branch("DRestrk3originalAlgo",DRestrk3originalAlgo,"DRestrk3originalAlgo[Dsize]/I");
				dnt->Branch("DRestrk4Pt",DRestrk4Pt,"DRestrk4Pt[Dsize]/F");
				dnt->Branch("DRestrk4Eta",DRestrk4Eta,"DRestrk4Eta[Dsize]/F");
				dnt->Branch("DRestrk4Phi",DRestrk4Phi,"DRestrk4Phi[Dsize]/F");
				dnt->Branch("DRestrk4P",DRestrk4P,"DRestrk4P[Dsize]/F");
				dnt->Branch("DRestrk4Dz",DRestrk4Dz,"DRestrk4Dz[Dsize]/F");
				dnt->Branch("DRestrk4Dxy",DRestrk4Dxy,"DRestrk4Dxy[Dsize]/F");
				dnt->Branch("DRestrk4originalAlgo",DRestrk4originalAlgo,"DRestrk4originalAlgo[Dsize]/I");
				dnt->Branch("DRestrk1PtErr",DRestrk1PtErr,"DRestrk1PtErr[Dsize]/F");
				dnt->Branch("DRestrk2PtErr",DRestrk2PtErr,"DRestrk2PtErr[Dsize]/F");
				dnt->Branch("DRestrk3PtErr",DRestrk3PtErr,"DRestrk3PtErr[Dsize]/F");
				dnt->Branch("DRestrk4PtErr",DRestrk4PtErr,"DRestrk4PtErr[Dsize]/F");
				dnt->Branch("DRestrk1highPurity",DRestrk1highPurity,"DRestrk1highPurity[Dsize]/O");
				dnt->Branch("DRestrk2highPurity",DRestrk2highPurity,"DRestrk2highPurity[Dsize]/O");
				dnt->Branch("DRestrk3highPurity",DRestrk3highPurity,"DRestrk3highPurity[Dsize]/O");
				dnt->Branch("DRestrk4highPurity",DRestrk4highPurity,"DRestrk4highPurity[Dsize]/O");
				dnt->Branch("DRestrk1dedx",DRestrk1dedx,"DRestrk1dedx[Dsize]/F");
				dnt->Branch("DRestrk2dedx",DRestrk2dedx,"DRestrk2dedx[Dsize]/F");
				dnt->Branch("DRestrk3dedx",DRestrk3dedx,"DRestrk3dedx[Dsize]/F");
				dnt->Branch("DRestrk4dedx",DRestrk4dedx,"DRestrk4dedx[Dsize]/F");
				dnt->Branch("DRestrk1thetastar",DRestrk1thetastar,"DRestrk1thetastar[Dsize]/F");
				dnt->Branch("DRestrk2thetastar",DRestrk2thetastar,"DRestrk2thetastar[Dsize]/F");
				dnt->Branch("DRestrk3thetastar",DRestrk3thetastar,"DRestrk3thetastar[Dsize]/F");
				dnt->Branch("DRestrk4thetastar",DRestrk4thetastar,"DRestrk4thetastar[Dsize]/F");
			}
			if(detailMode)
			{
				if(!D0kpimode)
				{
					dnt->Branch("DRestrk1Y",DRestrk1Y,"DRestrk1Y[Dsize]/F");
					dnt->Branch("DRestrk1D0",DRestrk1D0,"DRestrk1D0[Dsize]/F");
					dnt->Branch("DRestrk1D0Err",DRestrk1D0Err,"DRestrk1D0Err[Dsize]/F");
					dnt->Branch("DRestrk2Y",DRestrk2Y,"DRestrk2Y[Dsize]/F");
					dnt->Branch("DRestrk2D0",DRestrk2D0,"DRestrk2D0[Dsize]/F");
					dnt->Branch("DRestrk2D0Err",DRestrk2D0Err,"DRestrk2D0Err[Dsize]/F");
					dnt->Branch("DRestrk3Y",DRestrk3Y,"DRestrk3Y[Dsize]/F");
					dnt->Branch("DRestrk3D0",DRestrk3D0,"DRestrk3D0[Dsize]/F");
					dnt->Branch("DRestrk3D0Err",DRestrk3D0Err,"DRestrk3D0Err[Dsize]/F");
					dnt->Branch("DRestrk4Y",DRestrk4Y,"DRestrk4Y[Dsize]/F");
					dnt->Branch("DRestrk4D0",DRestrk4D0,"DRestrk4D0[Dsize]/F");
					dnt->Branch("DRestrk4D0Err",DRestrk4D0Err,"DRestrk4D0Err[Dsize]/F");
					dnt->Branch("DRestrk1Quality",DRestrk1Quality,"DRestrk1Quality[Dsize]/I");
					dnt->Branch("DRestrk2Quality",DRestrk2Quality,"DRestrk2Quality[Dsize]/I");
					dnt->Branch("DRestrk3Quality",DRestrk3Quality,"DRestrk3Quality[Dsize]/I");
					dnt->Branch("DRestrk4Quality",DRestrk4Quality,"DRestrk4Quality[Dsize]/I");
				}
			}
			//DInfo.genInfo
			dnt->Branch("Dgen",Dgen,"Dgen[Dsize]/F");
			dnt->Branch("DgenIndex",DgenIndex,"DgenIndex[Dsize]/I");
			dnt->Branch("DgennDa",DgennDa,"DgennDa[Dsize]/I");
			dnt->Branch("Dgenpt",Dgenpt,"Dgenpt[Dsize]/F");
			dnt->Branch("Dgeneta",Dgeneta,"Dgeneta[Dsize]/F");
			dnt->Branch("Dgenphi",Dgenphi,"Dgenphi[Dsize]/F");
			dnt->Branch("Dgeny",Dgeny,"Dgeny[Dsize]/F");
			dnt->Branch("DgencollisionId",DgencollisionId,"DgencollisionId[Dsize]/I");
			dnt->Branch("DgenBAncestorpt",DgenBAncestorpt,"DgenBAncestorpt[Dsize]/F");
			dnt->Branch("DgenBAncestorpdgId",DgenBAncestorpdgId,"DgenBAncestorpdgId[Dsize]/I");
			dnt->Branch("DgenprodvtxX",DgenprodvtxX,"DgenprodvtxX[Dsize]/F");
			dnt->Branch("DgenprodvtxY",DgenprodvtxY,"DgenprodvtxY[Dsize]/F");
			dnt->Branch("DgenprodvtxZ",DgenprodvtxZ,"DgenprodvtxZ[Dsize]/F");
			dnt->Branch("DgendecayvtxX",DgendecayvtxX,"DgendecayvtxX[Dsize]/F");
			dnt->Branch("DgendecayvtxY",DgendecayvtxY,"DgendecayvtxY[Dsize]/F");
			dnt->Branch("DgendecayvtxZ",DgendecayvtxZ,"DgendecayvtxZ[Dsize]/F");
			dnt->Branch("DgenfromgenPV",DgenfromgenPV,"DgenfromgenPV[Dsize]/I");
		}

		//GenInfo
		float   GPVx;
		float   GPVy;
		float   GPVz;
		int     Gsize;
		float   Gy[MAX_GEN];
		float   Geta[MAX_GEN];
		float   Gphi[MAX_GEN];
		float   Gpt[MAX_GEN];
		int     GpdgId[MAX_GEN];
		int     GcollisionId[MAX_GEN];
		int     GisSignal[MAX_GEN];
		float   GprodvtxX[MAX_GEN];//gen production vertex
		float   GprodvtxY[MAX_GEN];
		float   GprodvtxZ[MAX_GEN];
		float   GdecayvtxX[MAX_GEN];//gen decay vertex
		float   GdecayvtxY[MAX_GEN];
		float   GdecayvtxZ[MAX_GEN];
		float   GBAncestorpt[MAX_GEN];
		int     GBAncestorpdgId[MAX_GEN];
		int     GfromgenPV[MAX_GEN];
		float   Gtk1pt[MAX_GEN];
		float   Gtk1eta[MAX_GEN];
		float   Gtk1y[MAX_GEN];
		float   Gtk1phi[MAX_GEN];
		float   Gtk2pt[MAX_GEN];
		float   Gtk2eta[MAX_GEN];
		float   Gtk2y[MAX_GEN];
		float   Gtk2phi[MAX_GEN];
		float   Gtk3pt[MAX_GEN];
		float   Gtk3eta[MAX_GEN];
		float   Gtk3y[MAX_GEN];
		float   Gtk3phi[MAX_GEN];
		float   Gtk4pt[MAX_GEN];
		float   Gtk4eta[MAX_GEN];
		float   Gtk4y[MAX_GEN];
		float   Gtk4phi[MAX_GEN];
		float   GRestk1pt[MAX_GEN];
		float   GRestk1eta[MAX_GEN];
		float   GRestk1y[MAX_GEN];
		float   GRestk1phi[MAX_GEN];
		float   GRestk2pt[MAX_GEN];
		float   GRestk2eta[MAX_GEN];
		float   GRestk2y[MAX_GEN];
		float   GRestk2phi[MAX_GEN];
		float   GRestk3pt[MAX_GEN];
		float   GRestk3eta[MAX_GEN];
		float   GRestk3y[MAX_GEN];
		float   GRestk3phi[MAX_GEN];
		float   GRestk4pt[MAX_GEN];
		float   GRestk4eta[MAX_GEN];
		float   GRestk4y[MAX_GEN];
		float   GRestk4phi[MAX_GEN];  

		void buildGenBranch(TTree* nt)
		{
			nt->Branch("GPVx",&GPVx);
			nt->Branch("GPVy",&GPVy);
			nt->Branch("GPVz",&GPVz);
			nt->Branch("Gsize",&Gsize);
			nt->Branch("Gy",Gy,"Gy[Gsize]/F");
			nt->Branch("Geta",Geta,"Geta[Gsize]/F");
			nt->Branch("Gphi",Gphi,"Gphi[Gsize]/F");
			nt->Branch("Gpt",Gpt,"Gpt[Gsize]/F");
			nt->Branch("GpdgId",GpdgId,"GpdgId[Gsize]/I");
			nt->Branch("GcollisionId",GcollisionId,"GcollisionId[Gsize]/I");
			nt->Branch("GisSignal",GisSignal,"GisSignal[Gsize]/I");
			nt->Branch("GBAncestorpt",GBAncestorpt,"GBAncestorpt[Gsize]/F");
			nt->Branch("GBAncestorpdgId",GBAncestorpdgId,"GBAncestorpdgId[Gsize]/I");
			nt->Branch("GfromgenPV",GfromgenPV,"GfromgenPV[Gsize]/I");
			nt->Branch("GprodvtxX",GprodvtxX,"GprodvtxX[Gsize]/F");
			nt->Branch("GprodvtxY",GprodvtxY,"GprodvtxY[Gsize]/F");
			nt->Branch("GprodvtxZ",GprodvtxZ,"GprodvtxZ[Gsize]/F");
			nt->Branch("GdecayvtxX",GdecayvtxX,"GdecayvtxX[Gsize]/F");
			nt->Branch("GdecayvtxY",GdecayvtxY,"GdecayvtxY[Gsize]/F");
			nt->Branch("GdecayvtxZ",GdecayvtxZ,"GdecayvtxZ[Gsize]/F");
			nt->Branch("Gtk1pt",Gtk1pt,"Gtk1pt[Gsize]/F");
			nt->Branch("Gtk1eta",Gtk1eta,"Gtk1eta[Gsize]/F");
			nt->Branch("Gtk1y",Gtk1y,"Gtk1y[Gsize]/F");
			nt->Branch("Gtk1phi",Gtk1phi,"Gtk1phi[Gsize]/F");
			nt->Branch("Gtk2pt",Gtk2pt,"Gtk2pt[Gsize]/F");
			nt->Branch("Gtk2eta",Gtk2eta,"Gtk2eta[Gsize]/F");
			nt->Branch("Gtk2y",Gtk2y,"Gtk2y[Gsize]/F");
			nt->Branch("Gtk2phi",Gtk2phi,"Gtk2phi[Gsize]/F");
			nt->Branch("Gtk3pt",Gtk3pt,"Gtk3pt[Gsize]/F");
			nt->Branch("Gtk3eta",Gtk3eta,"Gtk3eta[Gsize]/F");
			nt->Branch("Gtk3y",Gtk3y,"Gtk3y[Gsize]/F");
			nt->Branch("Gtk3phi",Gtk3phi,"Gtk3phi[Gsize]/F");
			nt->Branch("Gtk4pt",Gtk4pt,"Gtk4pt[Gsize]/F");
			nt->Branch("Gtk4eta",Gtk4eta,"Gtk4eta[Gsize]/F");
			nt->Branch("Gtk4y",Gtk4y,"Gtk4y[Gsize]/F");
			nt->Branch("Gtk4phi",Gtk4phi,"Gtk4phi[Gsize]/F");
			nt->Branch("GRestk1pt",GRestk1pt,"GRestk1pt[Gsize]/F");
			nt->Branch("GRestk1eta",GRestk1eta,"GRestk1eta[Gsize]/F");
			nt->Branch("GRestk1y",GRestk1y,"GRestk1y[Gsize]/F");
			nt->Branch("GRestk1phi",GRestk1phi,"GRestk1phi[Gsize]/F");
			nt->Branch("GRestk2pt",GRestk2pt,"GRestk2pt[Gsize]/F");
			nt->Branch("GRestk2eta",GRestk2eta,"GRestk2eta[Gsize]/F");
			nt->Branch("GRestk2y",GRestk2y,"GRestk2y[Gsize]/F");
			nt->Branch("GRestk2phi",GRestk2phi,"GRestk2phi[Gsize]/F");
			nt->Branch("GRestk3pt",GRestk3pt,"GRestk3pt[Gsize]/F");
			nt->Branch("GRestk3eta",GRestk3eta,"GRestk3eta[Gsize]/F");
			nt->Branch("GRestk3y",GRestk3y,"GRestk3y[Gsize]/F");
			nt->Branch("GRestk3phi",GRestk3phi,"GRestk3phi[Gsize]/F");
			nt->Branch("GRestk4pt",GRestk4pt,"GRestk4pt[Gsize]/F");
			nt->Branch("GRestk4eta",GRestk4eta,"GRestk4eta[Gsize]/F");
			nt->Branch("GRestk4y",GRestk4y,"GRestk4y[Gsize]/F");
			nt->Branch("GRestk4phi",GRestk4phi,"GRestk4phi[Gsize]/F");
		}

		void makeDNtuple(int isDchannel[], int Dtypesize[], bool REAL, bool fillZeroCandEvt, bool skim, EvtInfoBranches *EvtInfo, VtxInfoBranches *VtxInfo, TrackInfoBranches *TrackInfo, DInfoBranches *DInfo, GenInfoBranches *GenInfo, TTree* ntD1, TTree* ntD2, TTree* ntD3, TTree* ntD4, TTree* ntD5, TTree* ntD6, TTree* ntD7)
		{//{{{
			std::vector<std::string> xtheInputVars;
			string xa1="DsvpvDistance/DsvpvDisErr";
			string xa2="Dchi2cl";
			string xa3 = "DsvpvDistance*TMath::Sin(Dalpha)";
			string xa4 = "TMath::Cos(Dtrk1thetastar)";
			string xa5="Dtrk1dedx";
			string xa6="Dtrk2dedx";
			string xa7 = "Dtrk1P";
			string xa8="Dtrk2P";
			string xa9 = "Dtrk1Pt";
			string xa10="Dtrk2Pt";

			vector<double> xinputValues;

			xtheInputVars.push_back(xa1);
			xtheInputVars.push_back(xa2);
			xtheInputVars.push_back(xa3);
			xtheInputVars.push_back(xa4); 
			xtheInputVars.push_back(xa5); 
			xtheInputVars.push_back(xa6); 
			xtheInputVars.push_back(xa7);
			xtheInputVars.push_back(xa8);
			xtheInputVars.push_back(xa9);
			xtheInputVars.push_back(xa10);

			TVector3* bP = new TVector3;
			TVector3* bVtx = new TVector3;
			TLorentzVector* b4P = new TLorentzVector;
			TVector3* boost = new TVector3();
			TVector3* D3Vec = new TVector3();
			fillTreeEvt(EvtInfo);
			bool zeroCand = true;
			for(int t=0;t<14;t++)
			{
				if(t%2==0)
				{
					Dsize = 0;
				}
				if(isDchannel[t]==1)
				{
					for(int j=0;j<DInfo->size;j++)
					{
						// b4P->SetPtEtaPhiM(DInfo->pt[j],DInfo->eta[j],DInfo->phi[j],DInfo->mass[j]);
						//
						//

						b4P->SetPtEtaPhiM(DInfo->pt[j],DInfo->eta[j],DInfo->phi[j],DInfo->mass[j]);
						double xDpt = DInfo->pt[j];
					//cout << "xDpt = " << xDpt << endl;
							double  xMVA = -999;
						if(xDpt < 2 && xDpt > 1)
						{
						bP->SetPtEtaPhi(DInfo->pt[j],DInfo->eta[j]*0,DInfo->phi[j]);
						bVtx->SetXYZ(DInfo->vtxX[j]-EvtInfo->PVx,
								DInfo->vtxY[j]-EvtInfo->PVy,
								DInfo->vtxZ[j]*0-EvtInfo->PVz*0);
						b4P->SetPtEtaPhiM(DInfo->pt[j],DInfo->eta[j],DInfo->phi[j],DInfo->mass[j]);
						boost->SetXYZ(b4P->BoostVector().X(), b4P->BoostVector().Y(), b4P->BoostVector().Z());
						D3Vec->SetXYZ(b4P->Vect().X(), b4P->Vect().Y(), b4P->Vect().Z());

						//          b4P->SetPtEtaPhiM(DInfo->pt[j],DInfo->eta[j],DInfo->phi[j],DInfo->mass[j]);
						//Calculate MVA Value//
						double trk1mass = DInfo->rftk1_MassHypo[j];
						double   trk2mass = DInfo->rftk2_MassHypo[j];
						double xDtrk1Idx = DInfo->rftk1_index[j];
						double xDtrk1Pt = TrackInfo->pt[DInfo->rftk1_index[j]];
						double xDtrk1Eta = TrackInfo->eta[DInfo->rftk1_index[j]];
						double xDtrk1Phi = TrackInfo->phi[DInfo->rftk1_index[j]];
						double xDtrk1PtErr = TrackInfo->ptErr[DInfo->rftk1_index[j]];
						double xDtrk1EtaErr = TrackInfo->etaErr[DInfo->rftk1_index[j]];
						double xDtrk1PhiErr = TrackInfo->phiErr[DInfo->rftk1_index[j]];
						b4P->SetPtEtaPhiM(TrackInfo->pt[DInfo->rftk1_index[j]],TrackInfo->eta[DInfo->rftk1_index[j]],TrackInfo->phi[DInfo->rftk1_index[j]],trk1mass);
						double xDtrk1Y = b4P->Rapidity();
						double xDtrk1P = b4P->P();
						b4P->Boost(-*boost);
						double xDtrk1thetastar = b4P->Angle(*D3Vec);
						double xDtrk1Dz = TrackInfo->dzPV[DInfo->rftk1_index[j]];
						double xDtrk1Dxy = TrackInfo->dxyPV[DInfo->rftk1_index[j]];
						double xDtrk1D0 = TrackInfo->d0[DInfo->rftk1_index[j]];
						double xDtrk1D0Err = TrackInfo->d0error[DInfo->rftk1_index[j]];
						double xDtrk1PixelHit = TrackInfo->pixelhit[DInfo->rftk1_index[j]];
						double xDtrk1StripHit = TrackInfo->striphit[DInfo->rftk1_index[j]];
						double xDtrk1nPixelLayer = TrackInfo->nPixelLayer[DInfo->rftk1_index[j]];
						double xDtrk1nStripLayer = TrackInfo->nStripLayer[DInfo->rftk1_index[j]];
						double xDtrk1MassHypo = DInfo->rftk1_MassHypo[j]*TrackInfo->charge[DInfo->rftk1_index[j]];
						double xDtrk1Chi2ndf = TrackInfo->chi2[DInfo->rftk1_index[j]]/TrackInfo->ndf[DInfo->rftk1_index[j]];
						double xDtrk1MVAVal = TrackInfo->trkMVAVal[DInfo->rftk1_index[j]];
						double xDtrk1Algo = TrackInfo->trkAlgo[DInfo->rftk1_index[j]];
						double xDtrk1originalAlgo = TrackInfo->originalTrkAlgo[DInfo->rftk1_index[j]];
						double xDtrk1highPurity = TrackInfo->highPurity[DInfo->rftk1_index[j]];
						double xDtrk1Quality = TrackInfo->trackQuality[DInfo->rftk1_index[j]];
						double xDtrk1dedx = TrackInfo->dedx[DInfo->rftk1_index[j]];
						double xDsvpvDistance= DInfo->svpvDistance[j];
						double xDtrk2Idx = DInfo->rftk2_index[j];
						double xDtrk2Pt = TrackInfo->pt[DInfo->rftk2_index[j]];
						double xDtrk2Eta = TrackInfo->eta[DInfo->rftk2_index[j]];
						double xDtrk2Phi = TrackInfo->phi[DInfo->rftk2_index[j]];
						double xDtrk2PtErr = TrackInfo->ptErr[DInfo->rftk2_index[j]];
						double xDtrk2EtaErr = TrackInfo->etaErr[DInfo->rftk2_index[j]];
						double xDtrk2PhiErr = TrackInfo->phiErr[DInfo->rftk2_index[j]];
						b4P->SetPtEtaPhiM(TrackInfo->pt[DInfo->rftk2_index[j]],TrackInfo->eta[DInfo->rftk2_index[j]],TrackInfo->phi[DInfo->rftk2_index[j]],trk2mass);
						double xDtrk2Y = b4P->Rapidity();
						double xDtrk2P = b4P->P();
						b4P->Boost(-*boost);
						double xDtrk2thetastar = b4P->Angle(*D3Vec);
						double xDtrk2Dz = TrackInfo->dzPV[DInfo->rftk2_index[j]];
						double xDtrk2Dxy = TrackInfo->dxyPV[DInfo->rftk2_index[j]];
						double xDtrk2D0 = TrackInfo->d0[DInfo->rftk2_index[j]];
						double xDtrk2D0Err = TrackInfo->d0error[DInfo->rftk2_index[j]];
						double xDtrk2PixelHit = TrackInfo->pixelhit[DInfo->rftk2_index[j]];
						double xDtrk2StripHit = TrackInfo->striphit[DInfo->rftk2_index[j]];
						double xDtrk2nPixelLayer = TrackInfo->nPixelLayer[DInfo->rftk2_index[j]];
						double xDtrk2nStripLayer = TrackInfo->nStripLayer[DInfo->rftk2_index[j]];
						double xDtrk2Chi2ndf = TrackInfo->chi2[DInfo->rftk2_index[j]]/TrackInfo->ndf[DInfo->rftk2_index[j]];
						double xDtrk2MassHypo = DInfo->rftk2_MassHypo[j]*TrackInfo->charge[DInfo->rftk2_index[j]];
						double xDtrk2MVAVal = TrackInfo->trkMVAVal[DInfo->rftk2_index[j]];
						double xDtrk2Algo = TrackInfo->trkAlgo[DInfo->rftk2_index[j]];
						double xDtrk2originalAlgo = TrackInfo->originalTrkAlgo[DInfo->rftk2_index[j]];
						double xDtrk2highPurity = TrackInfo->highPurity[DInfo->rftk2_index[j]];
						double xDtrk2Quality = TrackInfo->trackQuality[DInfo->rftk2_index[j]];
						double xDtrk2dedx = TrackInfo->dedx[DInfo->rftk2_index[j]];
						double xDalpha = DInfo->alpha[j];

						double xr2lxyBS = (DInfo->vtxX[j]-EvtInfo->BSx+(DInfo->vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdxdz) * (DInfo->vtxX[j]-EvtInfo->BSx+(DInfo->vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdxdz)
							+ (DInfo->vtxY[j]-EvtInfo->BSy+(DInfo->vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdydz) * (DInfo->vtxY[j]-EvtInfo->BSy+(DInfo->vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdydz);
						float xxlxyBS = DInfo->vtxX[j]-EvtInfo->BSx + (DInfo->vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdxdz;
						float xylxyBS = DInfo->vtxY[j]-EvtInfo->BSy + (DInfo->vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdydz;
						double xDlxyBS = TMath::Sqrt(xr2lxyBS);
						double xDlxyBSErr = (1./xr2lxyBS) * ((xxlxyBS*xxlxyBS)*DInfo->vtxXErr[j] + (2*xxlxyBS*xylxyBS)*DInfo->vtxYXErr[j] + (xylxyBS*xylxyBS)*DInfo->vtxYErr[j]);
						double xDchi2cl = TMath::Prob(DInfo->vtxchi2[j],DInfo->vtxdof[j]);
						double xDsvpvDisErr = DInfo->svpvDisErr[j];


						//MVA Calculations//
						xinputValues.clear();
						xinputValues.push_back(xDsvpvDistance/xDsvpvDisErr );
						xinputValues.push_back( xDchi2cl );
						xinputValues.push_back( xDsvpvDistance * TMath::Sin( xDalpha) );
						xinputValues.push_back(TMath::Cos( xDtrk1thetastar));
						xinputValues.push_back(xDtrk1dedx);
						xinputValues.push_back(xDtrk2dedx);
						xinputValues.push_back(xDtrk1P);
						xinputValues.push_back(xDtrk2P);
						xinputValues.push_back(xDtrk1Pt);
						xinputValues.push_back(xDtrk2Pt);
					/*
						   cout << "V1 Before " <<  xDsvpvDistance/xDsvpvDisErr  << endl;
						   cout << "V2 Before = " <<  xDchi2cl << endl;
						   cout << "V3 Before = " <<  xDsvpvDistance * TMath::Sin( xDalpha) << endl;
						   cout << "V4 Before = " <<  TMath::Cos( xDtrk1thetastar) << endl;
						   cout << "V5 Before = " <<xDtrk1dedx << endl;
						   cout << "V6 Before = " <<  xDtrk2dedx << endl;
						   cout << "V7 Before = " << xDtrk1P << endl;
						   cout << "V8 Before = " <<  xDtrk2P << endl;
										   cout << "V9 Before = " << xDtrk1Pt << endl;
						   cout << "V10 Before = " <<  xDtrk2Pt << endl;
*/
						//       b4P->SetPtEtaPhiM(TrackInfo->pt[DInfo->rftk2_index[j]],TrackInfo->eta[DInfo->rftk2_index[j]],TrackInfo->phi[DInfo->rftk2_index[j]],trk2mass);

						//       b4P->SetPtEtaPhiM(DInfo->pt[j],DInfo->eta[j],DInfo->phi[j],DInfo->mass[j]);


						/*if(xDpt < 0.5)  
						  {	
						  ReadBDTG mva(xtheInputVars);
						  xMVA=mva.GetMvaValue(xinputValues);    
						//if( t < 2)	cout << "Pt < 0.5 MVABefore MVA = " <<  xMVA << "  Channel = " << t <<  endl;

						//cout << "Dpt = " <<  DInfo->pt[j] << endl;
						}

						if(xDpt < 1.0 && xDpt > 0.5)   
						{
						ReadBDTG2 mva(xtheInputVars);
						xMVA2=mva.GetMvaValue(xinputValues);   
						//if( t < 2)		  cout << "0.5 < Pt < 1.0 MVABefore MVA = " <<  xMVA2 <<  "  Channel = " << t << endl;

						//			cout << "Dpt = " <<  DInfo->pt[j] << endl;

						}
						if(xDpt < 2.0 && xDpt > 1.0)   
						{
						ReadBDTG3 mva(xtheInputVars);
						xMVA3=mva.GetMvaValue(xinputValues);   
						//		if( t < 2)		  cout << "1.0 < Pt < 2.0 MVABefore MVA = " <<  xMVA3 <<  "  Channel = " << t << endl;
						//cout << "Dpt = " <<  DInfo->pt[j] << endl;

						}
						*/
							
						ReadBDTG mva(xtheInputVars);
							xMVA=mva.GetMvaValue(xinputValues);  
							//cout << "BDTG Before = " << xMVA << endl;
						}
						//		if( t < 2)		  cout << "1.0 < Pt < 2.0 MVABefore MVA = " <<  xMVA3 <<  "  Channel = " << t << endl;
						//cout << "Dpt = " <<  DInfo->pt[j] << endl;
				//		cout << "BDTG Before = " << xMVA << endl;


						b4P->SetPtEtaPhiM(DInfo->pt[j],DInfo->eta[j],DInfo->phi[j],DInfo->mass[j]);
						if(skim)
						{
							if(t == 0 || t == 1){
								/*
								   if(DInfo->pt[j] < 0.5) 
								   {
								   if(xMVA < 0.27) continue;
								   }

								   if(DInfo->pt[j] < 1.0 && DInfo->pt[j] > 0.5 ) 
								   {
								   if(xMVA2 < 0.10) continue;
								   }

								   if(DInfo->pt[j] < 2.0 && DInfo->pt[j] > 1.0 ) 
								   {
								   if(xMVA3 < 0.50) continue;
								   }
								   */

								//	cout << "PT After =" << DInfo->pt[j]<< endl;
								if(DInfo->pt[j] < 1.0) continue;

								if(DInfo->pt[j] < 2.0 && DInfo->pt[j] > 1.0) 
								{
									if(xMVA < 0.90) continue;
									if(fabs(b4P->Rapidity())>2.0) continue;
								}

								if(DInfo->pt[j] > 2.0 &&  DInfo->pt[j] < 3.0){ 
									if(fabs(b4P->Rapidity())>2.0) continue;
									if(fabs(TrackInfo->eta[DInfo->rftk1_index[j]])>2.0) continue;
									if(fabs(TrackInfo->eta[DInfo->rftk2_index[j]])>2.0) continue;
									if(TrackInfo->pt[DInfo->rftk1_index[j]]<0.2) continue;
									if(TrackInfo->pt[DInfo->rftk2_index[j]]<0.2) continue;
									if(TrackInfo->ptErr[DInfo->rftk1_index[j]]/TrackInfo->pt[DInfo->rftk1_index[j]] > 0.5) continue;
									if(TrackInfo->ptErr[DInfo->rftk2_index[j]]/TrackInfo->pt[DInfo->rftk2_index[j]] > 0.5) continue;
									if((DInfo->svpvDistance[j]/DInfo->svpvDisErr[j])< 3.91244) continue;
									if((DInfo->svpvDistance[j]* TMath::Sin(DInfo->alpha[j])) > 0.0111462) continue;
									if(TMath::Prob(DInfo->vtxchi2[j],DInfo->vtxdof[j]) < 0.272068) continue;
								}



								if(DInfo->pt[j] > 3.0 &&  DInfo->pt[j] < 4.0 ){ 
											if(fabs(b4P->Rapidity())>2.0) continue;
									if(fabs(TrackInfo->eta[DInfo->rftk1_index[j]])>2.0) continue;
									if(fabs(TrackInfo->eta[DInfo->rftk2_index[j]])>2.0) continue;
									if(TrackInfo->pt[DInfo->rftk1_index[j]]<0.2) continue;
									if(TrackInfo->pt[DInfo->rftk2_index[j]]<0.2) continue;
									if(TrackInfo->ptErr[DInfo->rftk1_index[j]]/TrackInfo->pt[DInfo->rftk1_index[j]] > 0.5) continue;
									if(TrackInfo->ptErr[DInfo->rftk2_index[j]]/TrackInfo->pt[DInfo->rftk2_index[j]] > 0.5) continue;
									if((DInfo->svpvDistance[j]/DInfo->svpvDisErr[j])<5.62064) continue;
									if((DInfo->svpvDistance[j]* TMath::Sin(DInfo->alpha[j])) > 0.00795336) continue;
									if(TMath::Prob(DInfo->vtxchi2[j],DInfo->vtxdof[j]) < 0.0697892) continue;

								}


								if(DInfo->pt[j] > 4.0 && DInfo->pt[j] < 5.0  ){ 
												if(fabs(b4P->Rapidity())>2.0) continue;
									if(fabs(TrackInfo->eta[DInfo->rftk1_index[j]])>2.0) continue;
									if(fabs(TrackInfo->eta[DInfo->rftk2_index[j]])>2.0) continue;
									if(TrackInfo->pt[DInfo->rftk1_index[j]]<0.2) continue;
									if(TrackInfo->pt[DInfo->rftk2_index[j]]<0.2) continue;
									if(TrackInfo->ptErr[DInfo->rftk1_index[j]]/TrackInfo->pt[DInfo->rftk1_index[j]] > 0.5) continue;
									if(TrackInfo->ptErr[DInfo->rftk2_index[j]]/TrackInfo->pt[DInfo->rftk2_index[j]] > 0.5) continue;
									if((DInfo->svpvDistance[j]/DInfo->svpvDisErr[j])<4.09504) continue;
									if((DInfo->svpvDistance[j]* TMath::Sin(DInfo->alpha[j])) > 0.00926084) continue;
									if(TMath::Prob(DInfo->vtxchi2[j],DInfo->vtxdof[j]) < 0.212428) continue;


								}


									if(DInfo->pt[j] > 5.0 && DInfo->pt[j] < 6.0  ){ 
									if(fabs(b4P->Rapidity())>2.0) continue;
									if(fabs(TrackInfo->eta[DInfo->rftk1_index[j]])>2.0) continue;
									if(fabs(TrackInfo->eta[DInfo->rftk2_index[j]])>2.0) continue;
									if(TrackInfo->pt[DInfo->rftk1_index[j]]<0.2) continue;
									if(TrackInfo->pt[DInfo->rftk2_index[j]]<0.2) continue;
									if(TrackInfo->ptErr[DInfo->rftk1_index[j]]/TrackInfo->pt[DInfo->rftk1_index[j]] > 0.5) continue;
									if(TrackInfo->ptErr[DInfo->rftk2_index[j]]/TrackInfo->pt[DInfo->rftk2_index[j]] > 0.5) continue;
									if((DInfo->svpvDistance[j]/DInfo->svpvDisErr[j])<3.47365) continue;
									if((DInfo->svpvDistance[j]* TMath::Sin(DInfo->alpha[j])) > 0.00884528) continue;
									if(TMath::Prob(DInfo->vtxchi2[j],DInfo->vtxdof[j]) < 0.0602953) continue;


								}


												if(DInfo->pt[j] > 6.0 && DInfo->pt[j] < 8.0  ){ 
									if(fabs(b4P->Rapidity())>2.0) continue;
									if(fabs(TrackInfo->eta[DInfo->rftk1_index[j]])>2.0) continue;
									if(fabs(TrackInfo->eta[DInfo->rftk2_index[j]])>2.0) continue;
									if(TrackInfo->pt[DInfo->rftk1_index[j]]<0.2) continue;
									if(TrackInfo->pt[DInfo->rftk2_index[j]]<0.2) continue;
									if(TrackInfo->ptErr[DInfo->rftk1_index[j]]/TrackInfo->pt[DInfo->rftk1_index[j]] > 0.5) continue;
									if(TrackInfo->ptErr[DInfo->rftk2_index[j]]/TrackInfo->pt[DInfo->rftk2_index[j]] > 0.5) continue;
									if((DInfo->svpvDistance[j]/DInfo->svpvDisErr[j])<3.8902) continue;
									if((DInfo->svpvDistance[j]* TMath::Sin(DInfo->alpha[j])) > 0.011302) continue;
									if(TMath::Prob(DInfo->vtxchi2[j],DInfo->vtxdof[j]) < 0.0525097) continue;


								}

																				if(DInfo->pt[j] > 8.0 && DInfo->pt[j] < 10.0  ){ 
									if(fabs(b4P->Rapidity())>2.0) continue;
									if(fabs(TrackInfo->eta[DInfo->rftk1_index[j]])>2.0) continue;
									if(fabs(TrackInfo->eta[DInfo->rftk2_index[j]])>2.0) continue;
									if(TrackInfo->pt[DInfo->rftk1_index[j]]<0.2) continue;
									if(TrackInfo->pt[DInfo->rftk2_index[j]]<0.2) continue;
									if(TrackInfo->ptErr[DInfo->rftk1_index[j]]/TrackInfo->pt[DInfo->rftk1_index[j]] > 0.5) continue;
									if(TrackInfo->ptErr[DInfo->rftk2_index[j]]/TrackInfo->pt[DInfo->rftk2_index[j]] > 0.5) continue;
									if((DInfo->svpvDistance[j]/DInfo->svpvDisErr[j])<2.94532) continue;
									if((DInfo->svpvDistance[j]* TMath::Sin(DInfo->alpha[j])) > 0.0137287) continue;
									if(TMath::Prob(DInfo->vtxchi2[j],DInfo->vtxdof[j]) < 0.0538837) continue;


								}

										if(DInfo->pt[j] > 10.0 && DInfo->pt[j] < 12.5  ){ 
									if(fabs(b4P->Rapidity())>2.0) continue;
									if(fabs(TrackInfo->eta[DInfo->rftk1_index[j]])>2.0) continue;
									if(fabs(TrackInfo->eta[DInfo->rftk2_index[j]])>2.0) continue;
									if(TrackInfo->pt[DInfo->rftk1_index[j]]<0.2) continue;
									if(TrackInfo->pt[DInfo->rftk2_index[j]]<0.2) continue;
									if(TrackInfo->ptErr[DInfo->rftk1_index[j]]/TrackInfo->pt[DInfo->rftk1_index[j]] > 0.5) continue;
									if(TrackInfo->ptErr[DInfo->rftk2_index[j]]/TrackInfo->pt[DInfo->rftk2_index[j]] > 0.5) continue;
									if((DInfo->svpvDistance[j]/DInfo->svpvDisErr[j])<3.3494) continue;
									if((DInfo->svpvDistance[j]* TMath::Sin(DInfo->alpha[j])) > 0.0103252) continue;
									if(TMath::Prob(DInfo->vtxchi2[j],DInfo->vtxdof[j]) < 0.0787636) continue;


								}


												if(DInfo->pt[j] > 12.5 && DInfo->pt[j] < 15.0  ){ 
									if(fabs(b4P->Rapidity())>2.0) continue;
									if(fabs(TrackInfo->eta[DInfo->rftk1_index[j]])>2.0) continue;
									if(fabs(TrackInfo->eta[DInfo->rftk2_index[j]])>2.0) continue;
									if(TrackInfo->pt[DInfo->rftk1_index[j]]<0.2) continue;
									if(TrackInfo->pt[DInfo->rftk2_index[j]]<0.2) continue;
									if(TrackInfo->ptErr[DInfo->rftk1_index[j]]/TrackInfo->pt[DInfo->rftk1_index[j]] > 0.5) continue;
									if(TrackInfo->ptErr[DInfo->rftk2_index[j]]/TrackInfo->pt[DInfo->rftk2_index[j]] > 0.5) continue;
									if((DInfo->svpvDistance[j]/DInfo->svpvDisErr[j])<3.49322) continue;
									if((DInfo->svpvDistance[j]* TMath::Sin(DInfo->alpha[j])) > 0.0104016) continue;
									if(TMath::Prob(DInfo->vtxchi2[j],DInfo->vtxdof[j]) < 0.0554402) continue;


								}

									if(DInfo->pt[j] > 12.5 && DInfo->pt[j] < 15.0  ){ 
									if(fabs(b4P->Rapidity())>2.0) continue;
									if(fabs(TrackInfo->eta[DInfo->rftk1_index[j]])>2.0) continue;
									if(fabs(TrackInfo->eta[DInfo->rftk2_index[j]])>2.0) continue;
									if(TrackInfo->pt[DInfo->rftk1_index[j]]<0.2) continue;
									if(TrackInfo->pt[DInfo->rftk2_index[j]]<0.2) continue;
									if(TrackInfo->ptErr[DInfo->rftk1_index[j]]/TrackInfo->pt[DInfo->rftk1_index[j]] > 0.5) continue;
									if(TrackInfo->ptErr[DInfo->rftk2_index[j]]/TrackInfo->pt[DInfo->rftk2_index[j]] > 0.5) continue;
									if((DInfo->svpvDistance[j]/DInfo->svpvDisErr[j])<3.42336) continue;
									if((DInfo->svpvDistance[j]* TMath::Sin(DInfo->alpha[j])) > 0.00975216) continue;
									if(TMath::Prob(DInfo->vtxchi2[j],DInfo->vtxdof[j]) < 0.0546604) continue;


								}




							}

						}

						if(DInfo->type[j]==(t+1))
						{
							fillDTree(bP,bVtx,b4P,boost,D3Vec,j,Dtypesize[t/2],REAL,EvtInfo,VtxInfo,TrackInfo,DInfo,GenInfo);
							Dtypesize[t/2]++;
						}
					}
					if(Dtypesize[t/2] != 0) zeroCand = false;
					else continue;
					if(t==1)       ntD1->Fill();
					else if(t==3)  ntD2->Fill();
					else if(t==5)  ntD3->Fill();
					else if(t==7)  ntD4->Fill();
					else if(t==9)  ntD5->Fill();
					else if(t==11) ntD6->Fill();
					else if(t==13) ntD7->Fill();
				}
			}

			Dsize = 0;
			for(int t = 1; t < 14; t+=2)
			{
				if(isDchannel[t]==1 && Dtypesize[t/2]==0)
				{
					if(!zeroCand || fillZeroCandEvt)
					{
						if(t==1)       ntD1->Fill();
						else if(t==3)  ntD2->Fill();
						else if(t==5)  ntD3->Fill();
						else if(t==7)  ntD4->Fill();
						else if(t==9)  ntD5->Fill();
						else if(t==11) ntD6->Fill();
						else if(t==13) ntD7->Fill();
					}
				}
			}

		}//}}}

		void fillDGenTree(TTree* ntGen, GenInfoBranches *GenInfo, bool gskim=true)
		{
			GPVx = GenInfo->genPVx;
			GPVy = GenInfo->genPVy;
			GPVz = GenInfo->genPVz;
			TLorentzVector* bGen = new TLorentzVector;
			int gt=0,sigtype=0;
			int gsize=0;
			int BAncestorindex=-99;
			Gsize = 0;
			for(int j=0;j<GenInfo->size;j++)
			{
				if(TMath::Abs(GenInfo->pdgId[j])!=DZERO_PDGID&&
						TMath::Abs(GenInfo->pdgId[j])!=DPLUS_PDGID&&
						TMath::Abs(GenInfo->pdgId[j])!=DSUBS_PDGID&&
						TMath::Abs(GenInfo->pdgId[j])!=DSTAR_PDGID&& 
						TMath::Abs(GenInfo->pdgId[j])!=BPLUS_PDGID&&gskim) continue;
				Gsize = gsize+1;
				Gpt[gsize] = GenInfo->pt[j];
				Geta[gsize] = GenInfo->eta[j];
				Gphi[gsize] = GenInfo->phi[j];
				GpdgId[gsize] = GenInfo->pdgId[j];
				GcollisionId[gsize] = GenInfo->collisionId[j];
				GprodvtxX[gsize] = GenInfo->vtxX[j];
				GprodvtxY[gsize] = GenInfo->vtxY[j];
				GprodvtxZ[gsize] = GenInfo->vtxZ[j];
				if( fabs(GprodvtxX[gsize]-GPVx) < 0.001 && fabs(GprodvtxY[gsize]-GPVy) < 0.001 && fabs(GprodvtxZ[gsize]-GPVz) < 0.001 )
					GfromgenPV[gsize] = 1;
				else
					GfromgenPV[gsize] = -1;
				bGen->SetPtEtaPhiM(GenInfo->pt[j],GenInfo->eta[j],GenInfo->phi[j],GenInfo->mass[j]);
				Gy[gsize] = bGen->Rapidity();
				sigtype=0;
				for(gt=1;gt<15;gt++)
				{
					if(isDsignalGen(gt,j,GenInfo))
					{
						sigtype=gt;
						break;
					}
				}
				GisSignal[gsize] = sigtype;
				GBAncestorpt[gsize] = -99.;
				GBAncestorpdgId[gsize] = 0;
				BAncestorindex = findBAncestor(j, GenInfo);
				if(BAncestorindex>=0)
				{
					GBAncestorpt[gsize] = GenInfo->pt[BAncestorindex];
					GBAncestorpdgId[gsize] = GenInfo->pdgId[BAncestorindex];
				}
				GdecayvtxX[gsize] = -999;
				GdecayvtxY[gsize] = -999;
				GdecayvtxZ[gsize] = -999;
				Gtk1pt[gsize] = -1;
				Gtk1eta[gsize] = -20;
				Gtk1phi[gsize] = -20;
				Gtk1y[gsize] = -1;
				Gtk2pt[gsize] = -1;
				Gtk2eta[gsize] = -20;
				Gtk2phi[gsize] = -20;
				Gtk2y[gsize] = -1;
				Gtk3pt[gsize] = -1;
				Gtk3eta[gsize] = -20;
				Gtk3phi[gsize] = -20;
				Gtk3y[gsize] = -1;
				Gtk4pt[gsize] = -1;
				Gtk4eta[gsize] = -20;
				Gtk4phi[gsize] = -20;
				Gtk4y[gsize] = -1;
				GRestk1pt[gsize] = -1;
				GRestk1eta[gsize] = -20;
				GRestk1phi[gsize] = -20;
				GRestk1y[gsize] = -1;
				GRestk2pt[gsize] = -1;
				GRestk2eta[gsize] = -20;
				GRestk2phi[gsize] = -20;
				GRestk2y[gsize] = -1;
				GRestk3pt[gsize] = -1;
				GRestk3eta[gsize] = -20;
				GRestk3phi[gsize] = -20;
				GRestk3y[gsize] = -1;
				GRestk4pt[gsize] = -1;
				GRestk4eta[gsize] = -20;
				GRestk4phi[gsize] = -20;
				GRestk4y[gsize] = -1;
				if(GisSignal[gsize]==1||GisSignal[gsize]==2||GisSignal[gsize]==3||GisSignal[gsize]==4||GisSignal[gsize]==5||GisSignal[gsize]==6)
				{
					GdecayvtxX[gsize] = GenInfo->vtxX[GenInfo->da1[j]];//all daughers should be from the same vertex, can be double checked here
					GdecayvtxY[gsize] = GenInfo->vtxY[GenInfo->da1[j]];
					GdecayvtxZ[gsize] = GenInfo->vtxZ[GenInfo->da1[j]];
					Gtk1pt[gsize] = GenInfo->pt[GenInfo->da1[j]];
					Gtk1eta[gsize] = GenInfo->eta[GenInfo->da1[j]];
					Gtk1phi[gsize] = GenInfo->phi[GenInfo->da1[j]];
					bGen->SetPtEtaPhiM(GenInfo->pt[GenInfo->da1[j]],GenInfo->eta[GenInfo->da1[j]],GenInfo->phi[GenInfo->da1[j]],GenInfo->mass[GenInfo->da1[j]]);
					Gtk1y[gsize] = bGen->Rapidity();
					Gtk2pt[gsize] = GenInfo->pt[GenInfo->da2[j]];
					Gtk2eta[gsize] = GenInfo->eta[GenInfo->da2[j]];
					Gtk2phi[gsize] = GenInfo->phi[GenInfo->da2[j]];
					bGen->SetPtEtaPhiM(GenInfo->pt[GenInfo->da2[j]],GenInfo->eta[GenInfo->da2[j]],GenInfo->phi[GenInfo->da2[j]],GenInfo->mass[GenInfo->da2[j]]);
					Gtk2y[gsize] = bGen->Rapidity();
					if(GisSignal[gsize]==3||GisSignal[gsize]==4||GisSignal[gsize]==5||GisSignal[gsize]==6)
					{
						Gtk3pt[gsize] = GenInfo->pt[GenInfo->da3[j]];
						Gtk3eta[gsize] = GenInfo->eta[GenInfo->da3[j]];
						Gtk3phi[gsize] = GenInfo->phi[GenInfo->da3[j]];
						bGen->SetPtEtaPhiM(GenInfo->pt[GenInfo->da3[j]],GenInfo->eta[GenInfo->da3[j]],GenInfo->phi[GenInfo->da3[j]],GenInfo->mass[GenInfo->da3[j]]);
						Gtk3y[gsize] = bGen->Rapidity();
						if(GisSignal[gsize]==5||GisSignal[gsize]==6)
						{
							Gtk4pt[gsize] = GenInfo->pt[GenInfo->da4[j]];
							Gtk4eta[gsize] = GenInfo->eta[GenInfo->da4[j]];
							Gtk4phi[gsize] = GenInfo->phi[GenInfo->da4[j]];
							bGen->SetPtEtaPhiM(GenInfo->pt[GenInfo->da4[j]],GenInfo->eta[GenInfo->da4[j]],GenInfo->phi[GenInfo->da4[j]],GenInfo->mass[GenInfo->da4[j]]);
							Gtk4y[gsize] = bGen->Rapidity();
						}
					}
				}
				if(GisSignal[gsize]==7||GisSignal[gsize]==8||GisSignal[gsize]==9||GisSignal[gsize]==10||GisSignal[gsize]==11||GisSignal[gsize]==12||GisSignal[gsize]==13||GisSignal[gsize]==14)
				{
					GdecayvtxX[gsize] = GenInfo->vtxX[GenInfo->da1[j]];
					GdecayvtxY[gsize] = GenInfo->vtxY[GenInfo->da1[j]];
					GdecayvtxZ[gsize] = GenInfo->vtxZ[GenInfo->da1[j]];
					Gtk1pt[gsize] = GenInfo->pt[GenInfo->da2[j]];
					Gtk1eta[gsize] = GenInfo->eta[GenInfo->da2[j]];
					Gtk1phi[gsize] = GenInfo->phi[GenInfo->da2[j]];
					bGen->SetPtEtaPhiM(GenInfo->pt[GenInfo->da2[j]],GenInfo->eta[GenInfo->da2[j]],GenInfo->phi[GenInfo->da2[j]],GenInfo->mass[GenInfo->da2[j]]);
					Gtk1y[gsize] = bGen->Rapidity();
					GRestk1pt[gsize] = GenInfo->pt[GenInfo->da1[GenInfo->da1[j]]];
					GRestk1eta[gsize] = GenInfo->eta[GenInfo->da1[GenInfo->da1[j]]];
					GRestk1phi[gsize] = GenInfo->phi[GenInfo->da1[GenInfo->da1[j]]];
					bGen->SetPtEtaPhiM(GenInfo->pt[GenInfo->da1[GenInfo->da1[j]]],GenInfo->eta[GenInfo->da1[GenInfo->da1[j]]],GenInfo->phi[GenInfo->da1[GenInfo->da1[j]]],GenInfo->mass[GenInfo->da1[GenInfo->da1[j]]]);
					GRestk1y[gsize] = bGen->Rapidity();
					GRestk2pt[gsize] = GenInfo->pt[GenInfo->da2[GenInfo->da1[j]]];
					GRestk2eta[gsize] = GenInfo->eta[GenInfo->da2[GenInfo->da1[j]]];
					GRestk2phi[gsize] = GenInfo->phi[GenInfo->da2[GenInfo->da1[j]]];
					bGen->SetPtEtaPhiM(GenInfo->pt[GenInfo->da2[GenInfo->da1[j]]],GenInfo->eta[GenInfo->da2[GenInfo->da1[j]]],GenInfo->phi[GenInfo->da2[GenInfo->da1[j]]],GenInfo->mass[GenInfo->da2[GenInfo->da1[j]]]);
					GRestk2y[gsize] = bGen->Rapidity();
					if(GisSignal[gsize]==11||GisSignal[gsize]==12)
					{
						GRestk3pt[gsize] = GenInfo->pt[GenInfo->da3[GenInfo->da1[j]]];
						GRestk3eta[gsize] = GenInfo->eta[GenInfo->da3[GenInfo->da1[j]]];
						GRestk3phi[gsize] = GenInfo->phi[GenInfo->da3[GenInfo->da1[j]]];
						bGen->SetPtEtaPhiM(GenInfo->pt[GenInfo->da3[GenInfo->da1[j]]],GenInfo->eta[GenInfo->da3[GenInfo->da1[j]]],GenInfo->phi[GenInfo->da3[GenInfo->da1[j]]],GenInfo->mass[GenInfo->da3[GenInfo->da1[j]]]);
						GRestk3y[gsize] = bGen->Rapidity();
						GRestk4pt[gsize] = GenInfo->pt[GenInfo->da4[GenInfo->da1[j]]];
						GRestk4eta[gsize] = GenInfo->eta[GenInfo->da4[GenInfo->da1[j]]];
						GRestk4phi[gsize] = GenInfo->phi[GenInfo->da4[GenInfo->da1[j]]];
						bGen->SetPtEtaPhiM(GenInfo->pt[GenInfo->da4[GenInfo->da1[j]]],GenInfo->eta[GenInfo->da4[GenInfo->da1[j]]],GenInfo->phi[GenInfo->da4[GenInfo->da1[j]]],GenInfo->mass[GenInfo->da4[GenInfo->da1[j]]]);
						GRestk4y[gsize] = bGen->Rapidity();
					}
				}
				gsize++;
			}
			ntGen->Fill();
		}
		/*
		   double findMass(int particlePdgId)
		   {
		   if(TMath::Abs(particlePdgId)==211) return PION_MASS;
		   if(TMath::Abs(particlePdgId)==321) return KAON_MASS;
		   else
		   {
		   cout<<"ERROR: find particle mass falied >> Particle pdgId: "<<particlePdgId<<endl;
		   return 0;
		   }
		   }
		   */
		int findPdgid(Double_t tkmass)
		{
			if(TMath::Abs(tkmass-KAON_MASS)<0.1) return KAON_PDGID;
			else if(TMath::Abs(tkmass-PION_MASS)<0.1) return PION_PDGID;
			else
			{
				cout<<"ERROR: find particle pdgid falied >> Particle mass: "<<tkmass<<endl;
				return 0;
			}
		}

		void fillTreeEvt(EvtInfoBranches *EvtInfo)
		{
			//Event Info
			RunNo = EvtInfo->RunNo;
			EvtNo = EvtInfo->EvtNo;
			LumiNo = EvtInfo->LumiNo;
			PVx = EvtInfo->PVx;
			PVy = EvtInfo->PVy;
			PVz = EvtInfo->PVz;
			PVxE = EvtInfo->PVxE;
			PVyE = EvtInfo->PVyE;
			PVzE = EvtInfo->PVzE;
			PVnchi2 = EvtInfo->PVnchi2;
			PVchi2 = EvtInfo->PVchi2;
			BSx = EvtInfo->BSx;
			BSy = EvtInfo->BSy;
			BSz = EvtInfo->BSz;
			BSxErr = EvtInfo->BSxErr;
			BSyErr = EvtInfo->BSyErr;
			BSzErr = EvtInfo->BSzErr;
			BSdxdz = EvtInfo->BSdxdz;
			BSdydz = EvtInfo->BSdydz;
			BSdxdzErr = EvtInfo->BSdxdzErr;
			BSdydzErr = EvtInfo->BSdydzErr;
			BSWidthX = EvtInfo->BSWidthX;
			BSWidthXErr = EvtInfo->BSWidthXErr;
			BSWidthY = EvtInfo->BSWidthY;
			BSWidthYErr = EvtInfo->BSWidthYErr;
		}

		void fillDTree(TVector3* bP, TVector3* bVtx, TLorentzVector* b4P, TVector3* boost, TVector3* D3Vec, int j, int typesize, bool REAL, EvtInfoBranches *EvtInfo, VtxInfoBranches *VtxInfo, TrackInfoBranches *TrackInfo, DInfoBranches *DInfo, GenInfoBranches *GenInfo)
		{
			//EvtInfo
			Dsize = typesize+1;

			//DInfo
			bP->SetPtEtaPhi(DInfo->pt[j],DInfo->eta[j]*0,DInfo->phi[j]);
			bVtx->SetXYZ(DInfo->vtxX[j]-EvtInfo->PVx,
					DInfo->vtxY[j]-EvtInfo->PVy,
					DInfo->vtxZ[j]*0-EvtInfo->PVz*0);
			b4P->SetPtEtaPhiM(DInfo->pt[j],DInfo->eta[j],DInfo->phi[j],DInfo->mass[j]);
			boost->SetXYZ(b4P->BoostVector().X(), b4P->BoostVector().Y(), b4P->BoostVector().Z());
			D3Vec->SetXYZ(b4P->Vect().X(), b4P->Vect().Y(), b4P->Vect().Z());
			Dindex[typesize] = typesize;
			Dtype[typesize] = DInfo->type[j];
			Dmass[typesize] = DInfo->mass[j];
			Dpt[typesize] = DInfo->pt[j];
			Deta[typesize] = DInfo->eta[j];
			Dphi[typesize] = DInfo->phi[j];
			Dy[typesize] = b4P->Rapidity();
			DvtxX[typesize] = DInfo->vtxX[j] - EvtInfo->PVx;
			DvtxY[typesize] = DInfo->vtxY[j] - EvtInfo->PVy;
			DvtxZ[typesize] = DInfo->vtxZ[j] - EvtInfo->PVz;
			Dd0[typesize] = TMath::Sqrt((DInfo->vtxX[j]-EvtInfo->PVx)*(DInfo->vtxX[j]-EvtInfo->PVx)+(DInfo->vtxY[j]-EvtInfo->PVy)*(DInfo->vtxY[j]-EvtInfo->PVy));
			Dd0Err[typesize] = TMath::Sqrt(DInfo->vtxXErr[j]*DInfo->vtxXErr[j]+DInfo->vtxYErr[j]*DInfo->vtxYErr[j]);
			Ddxyz[typesize] = TMath::Sqrt((DInfo->vtxX[j]-EvtInfo->PVx)*(DInfo->vtxX[j]-EvtInfo->PVx)+(DInfo->vtxY[j]-EvtInfo->PVy)*(DInfo->vtxY[j]-EvtInfo->PVy)+(DInfo->vtxZ[j]-EvtInfo->PVz)*(DInfo->vtxZ[j]-EvtInfo->PVz));
			DdxyzErr[typesize] = TMath::Sqrt(DInfo->vtxXErr[j]*DInfo->vtxXErr[j]+DInfo->vtxYErr[j]*DInfo->vtxYErr[j]+DInfo->vtxZErr[j]*DInfo->vtxZErr[j]);
			Dchi2ndf[typesize] = DInfo->vtxchi2[j]/DInfo->vtxdof[j];
			Dchi2cl[typesize] = TMath::Prob(DInfo->vtxchi2[j],DInfo->vtxdof[j]);
			Ddtheta[typesize] = bP->Angle(*bVtx);
			Dlxy[typesize] = ((DInfo->vtxX[j]-EvtInfo->PVx)*b4P->Px() + (DInfo->vtxY[j]-EvtInfo->PVy)*b4P->Py())/DInfo->pt[j];
			Dalpha[typesize] = DInfo->alpha[j];
			DsvpvDistance[typesize] = DInfo->svpvDistance[j];
			DsvpvDisErr[typesize] = DInfo->svpvDisErr[j];
			DsvpvDistance_2D[typesize] = DInfo->svpvDistance_2D[j];
			DsvpvDisErr_2D[typesize] = DInfo->svpvDisErr_2D[j];
			Ddca[typesize] = DInfo->svpvDistance[j]*TMath::Sin(DInfo->alpha[j]);

			float r2lxyBS = (DInfo->vtxX[j]-EvtInfo->BSx+(DInfo->vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdxdz) * (DInfo->vtxX[j]-EvtInfo->BSx+(DInfo->vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdxdz)
				+ (DInfo->vtxY[j]-EvtInfo->BSy+(DInfo->vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdydz) * (DInfo->vtxY[j]-EvtInfo->BSy+(DInfo->vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdydz);
			float xlxyBS = DInfo->vtxX[j]-EvtInfo->BSx + (DInfo->vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdxdz;
			float ylxyBS = DInfo->vtxY[j]-EvtInfo->BSy + (DInfo->vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdydz;
			DlxyBS[typesize] = TMath::Sqrt(r2lxyBS);
			DlxyBSErr[typesize] = (1./r2lxyBS) * ((xlxyBS*xlxyBS)*DInfo->vtxXErr[j] + (2*xlxyBS*ylxyBS)*DInfo->vtxYXErr[j] + (ylxyBS*ylxyBS)*DInfo->vtxYErr[j]);
			DMaxDoca[typesize] = DInfo->MaxDoca[j];

			//
			DtktkResmass[typesize] = -1;
			DtktkRespt[typesize] = -1;
			DtktkReseta[typesize] = -20;
			DtktkResphi[typesize] = -20;
			DtktkRes_chi2ndf[typesize] = -1;
			DtktkRes_chi2cl[typesize] = -1;
			DtktkRes_alpha[typesize] = -1;
			DtktkRes_alphaToSV[typesize] = -1;
			DtktkRes_svpvDistance[typesize] = -1;
			DtktkRes_svpvDisErr[typesize] = -1;
			DtktkRes_svpvDistanceToSV[typesize] = -1;
			DtktkRes_svpvDisErrToSV[typesize] = -1;
			DtktkRes_dca[typesize] = -1;
			DtktkRes_dcaToSV[typesize] = -1;
			DtktkRes_angleToTrk1[typesize] = -1;
			DtktkRes_ptAsymToTrk1[typesize] = -1;
			DtktkRes_unfitter_ptAsymToTrk1[typesize] = -1;
			DtktkRes_lxyBS[typesize] = -1;
			DtktkRes_lxyBSErr[typesize] = -1;

			DRestrk1Pt[typesize] = -1;
			DRestrk1Eta[typesize] = -20;
			DRestrk1Phi[typesize] = -20;
			DRestrk1P[typesize] = -1;
			DRestrk1Y[typesize] = -1;
			DRestrk1Dz[typesize] = -1;
			DRestrk1Dxy[typesize] = -1;
			DRestrk1D0[typesize] = -1;
			DRestrk1D0Err[typesize] = -1;
			DRestrk1originalAlgo[typesize] = 0;
			DRestrk2Pt[typesize] = -1;
			DRestrk2Eta[typesize] = -20;
			DRestrk2Phi[typesize] = -20;
			DRestrk2P[typesize] = -1;
			DRestrk2Y[typesize] = -1;
			DRestrk2Dz[typesize] = -1;
			DRestrk2Dxy[typesize] = -1;
			DRestrk2D0[typesize] = -1;
			DRestrk2D0Err[typesize] = -1;
			DRestrk2originalAlgo[typesize] = 0;
			DRestrk3Pt[typesize] = -1;
			DRestrk3Eta[typesize] = -20;
			DRestrk3Phi[typesize] = -20;
			DRestrk3P[typesize] = -1;
			DRestrk3Y[typesize] = -1;
			DRestrk3Dz[typesize] = -1;
			DRestrk3Dxy[typesize] = -1;
			DRestrk3D0[typesize] = -1;
			DRestrk3D0Err[typesize] = -1;
			DRestrk3originalAlgo[typesize] = 0;
			DRestrk4Pt[typesize] = -1;
			DRestrk4Eta[typesize] = -20;
			DRestrk4Phi[typesize] = -20;
			DRestrk4P[typesize] = -1;
			DRestrk4Y[typesize] = -1;
			DRestrk4Dz[typesize] = -1;
			DRestrk4Dxy[typesize] = -1;
			DRestrk4D0[typesize] = -1;
			DRestrk4D0Err[typesize] = -1;
			DRestrk4originalAlgo[typesize] = 0;
			DRestrk1PtErr[typesize] = -1;
			DRestrk2PtErr[typesize] = -1;
			DRestrk3PtErr[typesize] = -1;
			DRestrk4PtErr[typesize] = -1;
			DRestrk1highPurity[typesize] = -1;
			DRestrk2highPurity[typesize] = -1;
			DRestrk3highPurity[typesize] = -1;
			DRestrk4highPurity[typesize] = -1;
			DRestrk1Quality[typesize] = -1;
			DRestrk2Quality[typesize] = -1;
			DRestrk3Quality[typesize] = -1;
			DRestrk4Quality[typesize] = -1;
			DRestrk1dedx[typesize] = -20;
			DRestrk2dedx[typesize] = -20;
			DRestrk3dedx[typesize] = -20;
			DRestrk4dedx[typesize] = -20;
			DRestrk1thetastar[typesize] = -20;
			DRestrk2thetastar[typesize] = -20;
			DRestrk3thetastar[typesize] = -20;
			DRestrk4thetastar[typesize] = -20;
			BDTG_pt_1_2[typesize] = -99;	

			//MVA Part (BDT)
			//

			std::vector<std::string> theInputVars;
			std::string a1 = "DsvpvDistance/DsvpvDisErr";
			std::string a2 = "Dchi2cl";

			//  std::string a4 = "MVA";
			std::string a3 = "DsvpvDistance*TMath::Sin(Dalpha)";
			//	std::string a3 = "DsvpvDistance*TMath::Sin(Dalpha)";

			std::string a4 = "TMath::Cos(Dtrk1thetastar)";

			std::string a5 = "Dtrk1dedx";
			std::string a6 = "Dtrk2dedx";
			std::string a7 = "Dtrk1P";
			std::string a8 = "Dtrk2P";
			std::string a9 = "Dtrk1Pt";
			std::string a10 = "Dtrk2Pt";

	
			theInputVars.push_back(a1);
			theInputVars.push_back(a2);
			theInputVars.push_back(a3);
			theInputVars.push_back(a4); 
			theInputVars.push_back(a5); 
			theInputVars.push_back(a6); 
			theInputVars.push_back(a7);
			theInputVars.push_back(a8);
						theInputVars.push_back(a9);
			theInputVars.push_back(a10);

			vector<double> inputValues;


			//DInfo.trkInfo
			float trk1mass,trk2mass,trk3mass,trk4mass;
			if(DInfo->type[j]==1||DInfo->type[j]==2||DInfo->type[j]==3||DInfo->type[j]==4||DInfo->type[j]==5||DInfo->type[j]==6)
			{
				Dtrk1Idx[typesize] = DInfo->rftk1_index[j];
				Dtrk1Pt[typesize] = TrackInfo->pt[DInfo->rftk1_index[j]];
				Dtrk1Eta[typesize] = TrackInfo->eta[DInfo->rftk1_index[j]];
				Dtrk1Phi[typesize] = TrackInfo->phi[DInfo->rftk1_index[j]];
				Dtrk1PtErr[typesize] = TrackInfo->ptErr[DInfo->rftk1_index[j]];
				Dtrk1EtaErr[typesize] = TrackInfo->etaErr[DInfo->rftk1_index[j]];
				Dtrk1PhiErr[typesize] = TrackInfo->phiErr[DInfo->rftk1_index[j]];
				trk1mass = DInfo->rftk1_MassHypo[j];
				b4P->SetPtEtaPhiM(TrackInfo->pt[DInfo->rftk1_index[j]],TrackInfo->eta[DInfo->rftk1_index[j]],TrackInfo->phi[DInfo->rftk1_index[j]],trk1mass);
				Dtrk1Y[typesize] = b4P->Rapidity();
				Dtrk1P[typesize] = b4P->P();
				b4P->Boost(-*boost);
				Dtrk1thetastar[typesize] = b4P->Angle(*D3Vec);
				Dtrk1Dz[typesize] = TrackInfo->dzPV[DInfo->rftk1_index[j]];
				Dtrk1Dxy[typesize] = TrackInfo->dxyPV[DInfo->rftk1_index[j]];
				Dtrk1D0[typesize] = TrackInfo->d0[DInfo->rftk1_index[j]];
				Dtrk1D0Err[typesize] = TrackInfo->d0error[DInfo->rftk1_index[j]];
				Dtrk1PixelHit[typesize] = TrackInfo->pixelhit[DInfo->rftk1_index[j]];
				Dtrk1StripHit[typesize] = TrackInfo->striphit[DInfo->rftk1_index[j]];
				Dtrk1nPixelLayer[typesize] = TrackInfo->nPixelLayer[DInfo->rftk1_index[j]];
				Dtrk1nStripLayer[typesize] = TrackInfo->nStripLayer[DInfo->rftk1_index[j]];
				Dtrk1MassHypo[typesize] = DInfo->rftk1_MassHypo[j]*TrackInfo->charge[DInfo->rftk1_index[j]];
				Dtrk1Chi2ndf[typesize] = TrackInfo->chi2[DInfo->rftk1_index[j]]/TrackInfo->ndf[DInfo->rftk1_index[j]];
				Dtrk1MVAVal[typesize] = TrackInfo->trkMVAVal[DInfo->rftk1_index[j]];
				Dtrk1Algo[typesize] = TrackInfo->trkAlgo[DInfo->rftk1_index[j]];
				Dtrk1originalAlgo[typesize] = TrackInfo->originalTrkAlgo[DInfo->rftk1_index[j]];
				Dtrk1highPurity[typesize] = TrackInfo->highPurity[DInfo->rftk1_index[j]];
				Dtrk1Quality[typesize] = TrackInfo->trackQuality[DInfo->rftk1_index[j]];
				Dtrk1dedx[typesize] = TrackInfo->dedx[DInfo->rftk1_index[j]];

				Dtrk2Idx[typesize] = DInfo->rftk2_index[j];
				Dtrk2Pt[typesize] = TrackInfo->pt[DInfo->rftk2_index[j]];
				Dtrk2Eta[typesize] = TrackInfo->eta[DInfo->rftk2_index[j]];
				Dtrk2Phi[typesize] = TrackInfo->phi[DInfo->rftk2_index[j]];
				Dtrk2PtErr[typesize] = TrackInfo->ptErr[DInfo->rftk2_index[j]];
				Dtrk2EtaErr[typesize] = TrackInfo->etaErr[DInfo->rftk2_index[j]];
				Dtrk2PhiErr[typesize] = TrackInfo->phiErr[DInfo->rftk2_index[j]];
				trk2mass = DInfo->rftk2_MassHypo[j];
				b4P->SetPtEtaPhiM(TrackInfo->pt[DInfo->rftk2_index[j]],TrackInfo->eta[DInfo->rftk2_index[j]],TrackInfo->phi[DInfo->rftk2_index[j]],trk2mass);
				Dtrk2Y[typesize] = b4P->Rapidity();
				Dtrk2P[typesize] = b4P->P();
				b4P->Boost(-*boost);
				Dtrk2thetastar[typesize] = b4P->Angle(*D3Vec);
				Dtrk2Dz[typesize] = TrackInfo->dzPV[DInfo->rftk2_index[j]];
				Dtrk2Dxy[typesize] = TrackInfo->dxyPV[DInfo->rftk2_index[j]];
				Dtrk2D0[typesize] = TrackInfo->d0[DInfo->rftk2_index[j]];
				Dtrk2D0Err[typesize] = TrackInfo->d0error[DInfo->rftk2_index[j]];
				Dtrk2PixelHit[typesize] = TrackInfo->pixelhit[DInfo->rftk2_index[j]];
				Dtrk2StripHit[typesize] = TrackInfo->striphit[DInfo->rftk2_index[j]];
				Dtrk2nPixelLayer[typesize] = TrackInfo->nPixelLayer[DInfo->rftk2_index[j]];
				Dtrk2nStripLayer[typesize] = TrackInfo->nStripLayer[DInfo->rftk2_index[j]];
				Dtrk2Chi2ndf[typesize] = TrackInfo->chi2[DInfo->rftk2_index[j]]/TrackInfo->ndf[DInfo->rftk2_index[j]];
				Dtrk2MassHypo[typesize] = DInfo->rftk2_MassHypo[j]*TrackInfo->charge[DInfo->rftk2_index[j]];
				Dtrk2MVAVal[typesize] = TrackInfo->trkMVAVal[DInfo->rftk2_index[j]];
				Dtrk2Algo[typesize] = TrackInfo->trkAlgo[DInfo->rftk2_index[j]];
				Dtrk2originalAlgo[typesize] = TrackInfo->originalTrkAlgo[DInfo->rftk2_index[j]];
				Dtrk2highPurity[typesize] = TrackInfo->highPurity[DInfo->rftk2_index[j]];
				Dtrk2Quality[typesize] = TrackInfo->trackQuality[DInfo->rftk2_index[j]];
				Dtrk2dedx[typesize] = TrackInfo->dedx[DInfo->rftk2_index[j]];
				inputValues.clear();
				
				if(	Dpt[typesize] > 1 && 	Dpt[typesize] < 2)
				{
		//	cout << "typesize = " << typesize << endl;
					//    inputValues.push_back(Dalpha[j]);
					inputValues.push_back(DsvpvDistance[typesize]/DsvpvDisErr[typesize]);
				inputValues.push_back(Dchi2cl[typesize]);
				//		  inputValues.push_back(MVA[j]);
				inputValues.push_back(DsvpvDistance[typesize]*TMath::Sin(Dalpha[typesize]));
				inputValues.push_back(TMath::Cos(Dtrk1thetastar[typesize]));
				inputValues.push_back(Dtrk1dedx[typesize]);
				inputValues.push_back(Dtrk2dedx[typesize]);
				inputValues.push_back(Dtrk1P[typesize]);
				inputValues.push_back(Dtrk2P[typesize]);
				inputValues.push_back(Dtrk1Pt[typesize]);
				inputValues.push_back(Dtrk2Pt[typesize]);
		/*
					   cout << "V1 After " <<  DsvpvDistance[typesize]/DsvpvDisErr[typesize]  << endl;
						   cout << "V2 After = " <<  Dchi2cl[typesize] << endl;
						   cout << "V3 After = " <<  DsvpvDistance[typesize] * TMath::Sin( Dalpha[typesize]) << endl;
						   cout << "V4 After = " <<  TMath::Cos(Dtrk1thetastar[typesize]) << endl;
						   cout << "V5 After = " <<Dtrk1dedx[typesize] << endl;
						   cout << "V6 After = " <<  Dtrk2dedx[typesize] << endl;
						   cout << "V7 After = " << Dtrk1P[typesize] << endl;
						   cout << "V8 After = " <<  Dtrk2P[typesize] << endl;
						     cout << "V9 After = " << Dtrk1Pt[typesize] << endl;
						   cout << "V10 After = " <<  Dtrk2Pt[typesize] << endl;
*/

				ReadBDTG mva(theInputVars);
							BDTG_pt_1_2[typesize]=mva.GetMvaValue(inputValues);
		//		cout << "BDTG_pt_1_2 After = " << BDTG_pt_1_2[typesize] << endl;
				}
				//			cout << "Dpt After =" << Dpt[typesize] << endl;
			//cout << "BDTG_pt_1_2 After = " << BDTG_pt_1_2[typesize] << endl;


				if(DInfo->type[j]==1||DInfo->type[j]==2)
				{
					Dtrk3Idx[typesize] = -1;
					Dtrk3Pt[typesize] = -1;
					Dtrk3Eta[typesize] = -20;
					Dtrk3Phi[typesize] = -20;
					Dtrk3P[typesize] = -1;
					Dtrk3PtErr[typesize] = 0;
					Dtrk3EtaErr[typesize] = 0;
					Dtrk3PhiErr[typesize] = 0;
					Dtrk3Y[typesize] = -1;
					Dtrk3Dz[typesize] = -1;
					Dtrk3Dxy[typesize] = -1;
					Dtrk3D0[typesize] = -1;
					Dtrk3D0Err[typesize] = -1;
					Dtrk3PixelHit[typesize] = -1;
					Dtrk3StripHit[typesize] = -1;
					Dtrk3nPixelLayer[typesize] = -1;
					Dtrk3nStripLayer[typesize] = -1;
					Dtrk3Chi2ndf[typesize] = -1;
					Dtrk3MassHypo[typesize] = 0;
					Dtrk3MVAVal[typesize] = -100;
					Dtrk3Algo[typesize] = 0;
					Dtrk3originalAlgo[typesize] = 0;
					Dtrk3Quality[typesize] = 0;
					Dtrk3highPurity[typesize] = false;
					Dtrk3dedx[typesize] = -20;
					Dtrk3thetastar[typesize] = -20;
					Dtrk4Idx[typesize] = -1;
					Dtrk4Pt[typesize] = -1;
					Dtrk4Eta[typesize] = -20;
					Dtrk4Phi[typesize] = -20;
					Dtrk4P[typesize] = -1;
					Dtrk4PtErr[typesize] = 0;
					Dtrk4EtaErr[typesize] = 0;
					Dtrk4PhiErr[typesize] = 0;
					Dtrk4Y[typesize] = -1;
					Dtrk4Dz[typesize] = -1;
					Dtrk4Dxy[typesize] = -1;
					Dtrk4D0[typesize] = -1;
					Dtrk4D0Err[typesize] = -1;
					Dtrk4PixelHit[typesize] = -1;
					Dtrk4StripHit[typesize] = -1;
					Dtrk4nPixelLayer[typesize] = -1;
					Dtrk4nStripLayer[typesize] = -1;
					Dtrk4Chi2ndf[typesize] = -1;
					Dtrk4MassHypo[typesize] = 0;
					Dtrk4MVAVal[typesize] = -100;
					Dtrk4Algo[typesize] = 0;
					Dtrk4originalAlgo[typesize] = 0;
					Dtrk4Quality[typesize] = 0;
					Dtrk4highPurity[typesize] = false;
					Dtrk4dedx[typesize] = -20;
					Dtrk4thetastar[typesize] = -20;
				}
				else if(DInfo->type[j]==3||DInfo->type[j]==4)
				{
					Dtrk3Idx[typesize] = DInfo->rftk3_index[j];
					Dtrk3Pt[typesize] = TrackInfo->pt[DInfo->rftk3_index[j]];
					Dtrk3Eta[typesize] = TrackInfo->eta[DInfo->rftk3_index[j]];
					Dtrk3Phi[typesize] = TrackInfo->phi[DInfo->rftk3_index[j]];
					Dtrk3PtErr[typesize] = TrackInfo->ptErr[DInfo->rftk3_index[j]];
					Dtrk3EtaErr[typesize] = TrackInfo->etaErr[DInfo->rftk3_index[j]];
					Dtrk3PhiErr[typesize] = TrackInfo->phiErr[DInfo->rftk3_index[j]];
					trk3mass = DInfo->rftk3_MassHypo[j];
					b4P->SetPtEtaPhiM(TrackInfo->pt[DInfo->rftk3_index[j]],TrackInfo->eta[DInfo->rftk3_index[j]],TrackInfo->phi[DInfo->rftk3_index[j]],trk3mass);
					Dtrk3Y[typesize] = b4P->Rapidity();
					Dtrk3P[typesize] = b4P->P();
					b4P->Boost(-*boost);
					Dtrk3thetastar[typesize] = b4P->Angle(*D3Vec);
					Dtrk3Dz[typesize] = TrackInfo->dzPV[DInfo->rftk3_index[j]];
					Dtrk3Dxy[typesize] = TrackInfo->dxyPV[DInfo->rftk3_index[j]];
					Dtrk3D0[typesize] = TrackInfo->d0[DInfo->rftk3_index[j]];
					Dtrk3D0Err[typesize] = TrackInfo->d0error[DInfo->rftk3_index[j]];
					Dtrk3PixelHit[typesize] = TrackInfo->pixelhit[DInfo->rftk3_index[j]];
					Dtrk3StripHit[typesize] = TrackInfo->striphit[DInfo->rftk3_index[j]];
					Dtrk3nPixelLayer[typesize] = TrackInfo->nPixelLayer[DInfo->rftk3_index[j]];
					Dtrk3nStripLayer[typesize] = TrackInfo->nStripLayer[DInfo->rftk3_index[j]];
					Dtrk3Chi2ndf[typesize] = TrackInfo->chi2[DInfo->rftk3_index[j]]/TrackInfo->ndf[DInfo->rftk3_index[j]];
					Dtrk3MassHypo[typesize] = DInfo->rftk3_MassHypo[j]*TrackInfo->charge[DInfo->rftk3_index[j]];
					Dtrk3MVAVal[typesize] = TrackInfo->trkMVAVal[DInfo->rftk3_index[j]];
					Dtrk3Algo[typesize] = TrackInfo->trkAlgo[DInfo->rftk3_index[j]];
					Dtrk3originalAlgo[typesize] = TrackInfo->originalTrkAlgo[DInfo->rftk3_index[j]];
					Dtrk3highPurity[typesize] = TrackInfo->highPurity[DInfo->rftk3_index[j]];
					Dtrk3Quality[typesize] = TrackInfo->trackQuality[DInfo->rftk3_index[j]];
					Dtrk3dedx[typesize] = TrackInfo->dedx[DInfo->rftk3_index[j]];
					Dtrk4Idx[typesize] = -1;
					Dtrk4Pt[typesize] = -1;
					Dtrk4Eta[typesize] = -20;
					Dtrk4Phi[typesize] = -20;
					Dtrk4P[typesize] = -1;
					Dtrk4PtErr[typesize] = 0;
					Dtrk4EtaErr[typesize] = 0;
					Dtrk4PhiErr[typesize] = 0;
					Dtrk4Y[typesize] = -1;
					Dtrk4Dz[typesize] = -1;
					Dtrk4Dxy[typesize] = -1;
					Dtrk4D0[typesize] = -1;
					Dtrk4D0Err[typesize] = -1;
					Dtrk4PixelHit[typesize] = -1;
					Dtrk4StripHit[typesize] = -1;
					Dtrk4nPixelLayer[typesize] = -1;
					Dtrk4nStripLayer[typesize] = -1;
					Dtrk4Chi2ndf[typesize] = -1;
					Dtrk4MassHypo[typesize] = 0;
					Dtrk4MVAVal[typesize] = -100;
					Dtrk4Algo[typesize] = 0;
					Dtrk4originalAlgo[typesize] = 0;
					Dtrk4Quality[typesize] = 0;
					Dtrk4highPurity[typesize] = false;
					Dtrk4dedx[typesize] = -20;
					Dtrk4thetastar[typesize] = -20;
				}
				else if(DInfo->type[j]==5||DInfo->type[j]==6)
				{
					Dtrk3Idx[typesize] = DInfo->rftk3_index[j];
					Dtrk4Idx[typesize] = DInfo->rftk4_index[j];
					Dtrk3Pt[typesize] = TrackInfo->pt[DInfo->rftk3_index[j]];
					Dtrk4Pt[typesize] = TrackInfo->pt[DInfo->rftk4_index[j]];
					Dtrk3Eta[typesize] = TrackInfo->eta[DInfo->rftk3_index[j]];
					Dtrk4Eta[typesize] = TrackInfo->eta[DInfo->rftk4_index[j]];
					Dtrk3Phi[typesize] = TrackInfo->phi[DInfo->rftk3_index[j]];
					Dtrk4Phi[typesize] = TrackInfo->phi[DInfo->rftk4_index[j]];
					Dtrk3PtErr[typesize] = TrackInfo->ptErr[DInfo->rftk3_index[j]];
					Dtrk4PtErr[typesize] = TrackInfo->ptErr[DInfo->rftk4_index[j]];
					Dtrk3EtaErr[typesize] = TrackInfo->etaErr[DInfo->rftk3_index[j]];
					Dtrk4EtaErr[typesize] = TrackInfo->etaErr[DInfo->rftk4_index[j]];
					Dtrk3PhiErr[typesize] = TrackInfo->phiErr[DInfo->rftk3_index[j]];
					Dtrk4PhiErr[typesize] = TrackInfo->phiErr[DInfo->rftk4_index[j]];
					trk3mass = DInfo->rftk3_MassHypo[j];
					trk4mass = DInfo->rftk4_MassHypo[j];
					b4P->SetPtEtaPhiM(TrackInfo->pt[DInfo->rftk3_index[j]],TrackInfo->eta[DInfo->rftk3_index[j]],TrackInfo->phi[DInfo->rftk3_index[j]],trk3mass);
					Dtrk3Y[typesize] = b4P->Rapidity();
					Dtrk3P[typesize] = b4P->P();
					b4P->Boost(-*boost);
					Dtrk3thetastar[typesize] = b4P->Angle(*D3Vec);
					b4P->SetPtEtaPhiM(TrackInfo->pt[DInfo->rftk4_index[j]],TrackInfo->eta[DInfo->rftk4_index[j]],TrackInfo->phi[DInfo->rftk4_index[j]],trk4mass);
					Dtrk4Y[typesize] = b4P->Rapidity();
					Dtrk4P[typesize] = b4P->P();
					b4P->Boost(-*boost);
					Dtrk4thetastar[typesize] = b4P->Angle(*D3Vec);
					Dtrk3Dz[typesize] = TrackInfo->dzPV[DInfo->rftk3_index[j]];
					Dtrk4Dz[typesize] = TrackInfo->dzPV[DInfo->rftk4_index[j]];
					Dtrk3Dxy[typesize] = TrackInfo->dxyPV[DInfo->rftk3_index[j]];
					Dtrk4Dxy[typesize] = TrackInfo->dxyPV[DInfo->rftk4_index[j]];
					Dtrk3D0[typesize] = TrackInfo->d0[DInfo->rftk3_index[j]];
					Dtrk4D0[typesize] = TrackInfo->d0[DInfo->rftk4_index[j]];
					Dtrk3D0Err[typesize] = TrackInfo->d0error[DInfo->rftk3_index[j]];
					Dtrk4D0Err[typesize] = TrackInfo->d0error[DInfo->rftk4_index[j]];
					Dtrk3PixelHit[typesize] = TrackInfo->pixelhit[DInfo->rftk3_index[j]];
					Dtrk4PixelHit[typesize] = TrackInfo->pixelhit[DInfo->rftk4_index[j]];
					Dtrk3StripHit[typesize] = TrackInfo->striphit[DInfo->rftk3_index[j]];
					Dtrk4StripHit[typesize] = TrackInfo->striphit[DInfo->rftk4_index[j]];
					Dtrk3nPixelLayer[typesize] = TrackInfo->nPixelLayer[DInfo->rftk3_index[j]];
					Dtrk4nPixelLayer[typesize] = TrackInfo->nPixelLayer[DInfo->rftk4_index[j]];
					Dtrk3nStripLayer[typesize] = TrackInfo->nStripLayer[DInfo->rftk3_index[j]];
					Dtrk4nStripLayer[typesize] = TrackInfo->nStripLayer[DInfo->rftk4_index[j]];
					Dtrk3Chi2ndf[typesize] = TrackInfo->chi2[DInfo->rftk3_index[j]]/TrackInfo->ndf[DInfo->rftk3_index[j]];
					Dtrk4Chi2ndf[typesize] = TrackInfo->chi2[DInfo->rftk4_index[j]]/TrackInfo->ndf[DInfo->rftk4_index[j]];
					Dtrk3MassHypo[typesize] = DInfo->rftk3_MassHypo[j]*TrackInfo->charge[DInfo->rftk3_index[j]];
					Dtrk4MassHypo[typesize] = DInfo->rftk4_MassHypo[j]*TrackInfo->charge[DInfo->rftk4_index[j]];
					Dtrk3MVAVal[typesize] = TrackInfo->trkMVAVal[DInfo->rftk3_index[j]];
					Dtrk4MVAVal[typesize] = TrackInfo->trkMVAVal[DInfo->rftk4_index[j]];
					Dtrk3Algo[typesize] = TrackInfo->trkAlgo[DInfo->rftk3_index[j]];
					Dtrk4Algo[typesize] = TrackInfo->trkAlgo[DInfo->rftk4_index[j]];
					Dtrk3originalAlgo[typesize] = TrackInfo->originalTrkAlgo[DInfo->rftk3_index[j]];
					Dtrk4originalAlgo[typesize] = TrackInfo->originalTrkAlgo[DInfo->rftk4_index[j]];
					Dtrk3highPurity[typesize] = TrackInfo->highPurity[DInfo->rftk3_index[j]];
					Dtrk4highPurity[typesize] = TrackInfo->highPurity[DInfo->rftk4_index[j]];
					Dtrk3Quality[typesize] = TrackInfo->trackQuality[DInfo->rftk3_index[j]];
					Dtrk4Quality[typesize] = TrackInfo->trackQuality[DInfo->rftk4_index[j]];
					Dtrk3dedx[typesize] = TrackInfo->dedx[DInfo->rftk3_index[j]];
					Dtrk4dedx[typesize] = TrackInfo->dedx[DInfo->rftk4_index[j]];
				}
			}
				else if(DInfo->type[j]==7||DInfo->type[j]==8||DInfo->type[j]==9||DInfo->type[j]==10||DInfo->type[j]==11||DInfo->type[j]==12||DInfo->type[j]==13||DInfo->type[j]==14)
				{
					Dtrk1Idx[typesize] = DInfo->rftk2_index[j];
					Dtrk1Pt[typesize] = TrackInfo->pt[DInfo->rftk2_index[j]];
					Dtrk1Eta[typesize] = TrackInfo->eta[DInfo->rftk2_index[j]];
					Dtrk1Phi[typesize] = TrackInfo->phi[DInfo->rftk2_index[j]];
					Dtrk1PtErr[typesize] = TrackInfo->ptErr[DInfo->rftk2_index[j]];
					Dtrk1EtaErr[typesize] = TrackInfo->etaErr[DInfo->rftk2_index[j]];
					Dtrk1PhiErr[typesize] = TrackInfo->phiErr[DInfo->rftk2_index[j]];
					b4P->SetPtEtaPhiM(TrackInfo->pt[DInfo->rftk2_index[j]],TrackInfo->eta[DInfo->rftk2_index[j]],TrackInfo->phi[DInfo->rftk2_index[j]],PION_MASS);
					Dtrk1Y[typesize] = b4P->Rapidity();
					Dtrk1P[typesize] = b4P->P();
					b4P->Boost(-*boost);
					Dtrk1thetastar[typesize] = b4P->Angle(*D3Vec);
					Dtrk1Dz[typesize] = TrackInfo->dzPV[DInfo->rftk2_index[j]];
					Dtrk1Dxy[typesize] = TrackInfo->dxyPV[DInfo->rftk2_index[j]];
					Dtrk1D0[typesize] = TrackInfo->d0[DInfo->rftk2_index[j]];
					Dtrk1D0Err[typesize] = TrackInfo->d0error[DInfo->rftk2_index[j]];
					Dtrk1PixelHit[typesize] = TrackInfo->pixelhit[DInfo->rftk2_index[j]];
					Dtrk1StripHit[typesize] = TrackInfo->striphit[DInfo->rftk2_index[j]];
					Dtrk1nPixelLayer[typesize] = TrackInfo->nPixelLayer[DInfo->rftk2_index[j]];
					Dtrk1nStripLayer[typesize] = TrackInfo->nStripLayer[DInfo->rftk2_index[j]];
					Dtrk1Chi2ndf[typesize] = TrackInfo->chi2[DInfo->rftk2_index[j]]/TrackInfo->ndf[DInfo->rftk2_index[j]];
					Dtrk1MassHypo[typesize] = DInfo->rftk2_MassHypo[j]*TrackInfo->charge[DInfo->rftk2_index[j]];
					Dtrk1MVAVal[typesize] = TrackInfo->trkMVAVal[DInfo->rftk2_index[j]];
					Dtrk1Algo[typesize] = TrackInfo->trkAlgo[DInfo->rftk2_index[j]];
					Dtrk1originalAlgo[typesize] = TrackInfo->originalTrkAlgo[DInfo->rftk2_index[j]];
					Dtrk1highPurity[typesize] = TrackInfo->highPurity[DInfo->rftk2_index[j]];
					Dtrk1Quality[typesize] = TrackInfo->trackQuality[DInfo->rftk2_index[j]];
					Dtrk1dedx[typesize] = TrackInfo->dedx[DInfo->rftk2_index[j]];

					Dtrk2Idx[typesize] = -1;
					Dtrk2Pt[typesize] = -1;
					Dtrk2Eta[typesize] = -20;
					Dtrk2Phi[typesize] = -20;
					Dtrk2P[typesize] = -1;
					Dtrk2PtErr[typesize] = 0;
					Dtrk2EtaErr[typesize] = 0;
					Dtrk2PhiErr[typesize] = 0;
					Dtrk2Y[typesize] = -1;
					Dtrk2Dz[typesize] = -1;
					Dtrk2Dxy[typesize] = -1;
					Dtrk2D0[typesize] = -1;
					Dtrk2D0Err[typesize] = -1;
					Dtrk2PixelHit[typesize] = -1;
					Dtrk2StripHit[typesize] = -1;
					Dtrk2nPixelLayer[typesize] = -1;
					Dtrk2nStripLayer[typesize] = -1;
					Dtrk2Chi2ndf[typesize] = -1;
					Dtrk2MassHypo[typesize] = 0;
					Dtrk2MVAVal[typesize] = -100;
					Dtrk2Algo[typesize] = 0;
					Dtrk2originalAlgo[typesize] = 0;
					Dtrk2Quality[typesize] = 0;
					Dtrk2highPurity[typesize] = false;
					Dtrk2dedx[typesize] = -20;
					Dtrk2thetastar[typesize] = -20;
					Dtrk3Idx[typesize] = -1;
					Dtrk3Pt[typesize] = -1;
					Dtrk3Eta[typesize] = -20;
					Dtrk3Phi[typesize] = -20;
					Dtrk3P[typesize] = -1;
					Dtrk3PtErr[typesize] = 0;
					Dtrk3EtaErr[typesize] = 0;
					Dtrk3PhiErr[typesize] = 0;
					Dtrk3Y[typesize] = -1;
					Dtrk3Dz[typesize] = -1;
					Dtrk3Dxy[typesize] = -1;
					Dtrk3D0[typesize] = -1;
					Dtrk3D0Err[typesize] = -1;
					Dtrk3PixelHit[typesize] = -1;
					Dtrk3StripHit[typesize] = -1;
					Dtrk3nPixelLayer[typesize] = -1;
					Dtrk3nStripLayer[typesize] = -1;
					Dtrk3Chi2ndf[typesize] = -1;
					Dtrk3MassHypo[typesize] = 0;
					Dtrk3MVAVal[typesize] = -100;
					Dtrk3Algo[typesize] = 0;
					Dtrk3originalAlgo[typesize] = 0;
					Dtrk3Quality[typesize] = 0;
					Dtrk3highPurity[typesize] = false;
					Dtrk3dedx[typesize] = -20;
					Dtrk3thetastar[typesize] = -20;
					Dtrk4Idx[typesize] = -1;
					Dtrk4Pt[typesize] = -1;
					Dtrk4Eta[typesize] = -20;
					Dtrk4Phi[typesize] = -20;
					Dtrk4P[typesize] = -1;
					Dtrk4PtErr[typesize] = 0;
					Dtrk4EtaErr[typesize] = 0;
					Dtrk4PhiErr[typesize] = 0;
					Dtrk4Y[typesize] = -1;
					Dtrk4Dz[typesize] = -1;
					Dtrk4Dxy[typesize] = -1;
					Dtrk4D0[typesize] = -1;
					Dtrk4D0Err[typesize] = -1;
					Dtrk4PixelHit[typesize] = -1;
					Dtrk4StripHit[typesize] = -1;
					Dtrk4nPixelLayer[typesize] = -1;
					Dtrk4nStripLayer[typesize] = -1;
					Dtrk4Chi2ndf[typesize] = -1;
					Dtrk4MassHypo[typesize] = 0;
					Dtrk4MVAVal[typesize] = -100;
					Dtrk4Algo[typesize] = 0;
					Dtrk4originalAlgo[typesize] = 0;
					Dtrk4Quality[typesize] = 0;
					Dtrk4highPurity[typesize] = false;
					Dtrk4dedx[typesize] = -20;
					Dtrk4thetastar[typesize] = -20;

					DRestrk1Pt[typesize] = TrackInfo->pt[DInfo->tktkRes_rftk1_index[j]];
					DRestrk1Eta[typesize] = TrackInfo->eta[DInfo->tktkRes_rftk1_index[j]];
					DRestrk1Phi[typesize] = TrackInfo->phi[DInfo->tktkRes_rftk1_index[j]];
					b4P->SetPtEtaPhiM(TrackInfo->pt[DInfo->tktkRes_rftk1_index[j]],TrackInfo->eta[DInfo->tktkRes_rftk1_index[j]],TrackInfo->phi[DInfo->tktkRes_rftk1_index[j]],DInfo->tktkRes_rftk1_mass[j]);
					DRestrk1Y[typesize] = b4P->Rapidity();
					DRestrk1P[typesize] = b4P->P();
					b4P->Boost(-*boost);
					DRestrk1thetastar[typesize] = b4P->Angle(*D3Vec);
					DRestrk1Dz[typesize] = TrackInfo->dzPV[DInfo->tktkRes_rftk1_index[j]];
					DRestrk1Dxy[typesize] = TrackInfo->dxyPV[DInfo->tktkRes_rftk1_index[j]];
					DRestrk1D0[typesize] = TrackInfo->d0[DInfo->tktkRes_rftk1_index[j]];
					DRestrk1D0Err[typesize] = TrackInfo->d0error[DInfo->tktkRes_rftk1_index[j]];
					DRestrk1originalAlgo[typesize] = TrackInfo->originalTrkAlgo[DInfo->tktkRes_rftk1_index[j]];
					DRestrk2Pt[typesize] = TrackInfo->pt[DInfo->tktkRes_rftk2_index[j]];
					DRestrk2Eta[typesize] = TrackInfo->eta[DInfo->tktkRes_rftk2_index[j]];
					DRestrk2Phi[typesize] = TrackInfo->phi[DInfo->tktkRes_rftk2_index[j]];
					b4P->SetPtEtaPhiM(TrackInfo->pt[DInfo->tktkRes_rftk2_index[j]],TrackInfo->eta[DInfo->tktkRes_rftk2_index[j]],TrackInfo->phi[DInfo->tktkRes_rftk2_index[j]],DInfo->tktkRes_rftk2_mass[j]);
					DRestrk2Y[typesize] = b4P->Rapidity();
					DRestrk2P[typesize] = b4P->P();
					b4P->Boost(-*boost);
					DRestrk2thetastar[typesize] = b4P->Angle(*D3Vec);
					DRestrk2Dz[typesize] = TrackInfo->dzPV[DInfo->tktkRes_rftk2_index[j]];
					DRestrk2Dxy[typesize] = TrackInfo->dxyPV[DInfo->tktkRes_rftk2_index[j]];
					DRestrk2D0[typesize] = TrackInfo->d0[DInfo->tktkRes_rftk2_index[j]];
					DRestrk2D0Err[typesize] = TrackInfo->d0error[DInfo->tktkRes_rftk2_index[j]];
					DRestrk2originalAlgo[typesize] = TrackInfo->originalTrkAlgo[DInfo->tktkRes_rftk2_index[j]];
					DRestrk3Pt[typesize] = -1;
					DRestrk3Eta[typesize] = -20;
					DRestrk3Phi[typesize] = -20;
					DRestrk3P[typesize] = -1;
					DRestrk3Y[typesize] = -1;
					DRestrk3Dz[typesize] = -1;
					DRestrk3Dxy[typesize] = -1;
					DRestrk3D0[typesize] = -1;
					DRestrk3D0Err[typesize] = -1;
					DRestrk3originalAlgo[typesize] = 0;
					DRestrk4Pt[typesize] = -1;
					DRestrk4Eta[typesize] = -20;
					DRestrk4Phi[typesize] = -20;
					DRestrk4P[typesize] = -1;
					DRestrk4Y[typesize] = -1;
					DRestrk4Dz[typesize] = -1;
					DRestrk4Dxy[typesize] = -1;
					DRestrk4D0[typesize] = -1;
					DRestrk4D0Err[typesize] = -1;
					DRestrk4originalAlgo[typesize] = 0;
					DRestrk1PtErr[typesize] = TrackInfo->ptErr[DInfo->tktkRes_rftk1_index[j]];
					DRestrk2PtErr[typesize] = TrackInfo->ptErr[DInfo->tktkRes_rftk2_index[j]];
					DRestrk3PtErr[typesize] = -1;
					DRestrk4PtErr[typesize] = -1;
					DRestrk1highPurity[typesize] = TrackInfo->highPurity[DInfo->tktkRes_rftk1_index[j]];
					DRestrk2highPurity[typesize] = TrackInfo->highPurity[DInfo->tktkRes_rftk2_index[j]];
					DRestrk3highPurity[typesize] = -1;
					DRestrk4highPurity[typesize] = -1;
					DRestrk1Quality[typesize] = TrackInfo->trackQuality[DInfo->tktkRes_rftk1_index[j]];
					DRestrk2Quality[typesize] = TrackInfo->trackQuality[DInfo->tktkRes_rftk2_index[j]];
					DRestrk3Quality[typesize] = -1;
					DRestrk4Quality[typesize] = -1;
					DRestrk1dedx[typesize] = TrackInfo->dedx[DInfo->tktkRes_rftk1_index[j]];
					DRestrk2dedx[typesize] = TrackInfo->dedx[DInfo->tktkRes_rftk2_index[j]];
					DRestrk3dedx[typesize] = -20;
					DRestrk4dedx[typesize] = -20;
					DRestrk3thetastar[typesize] = -20;
					DRestrk4thetastar[typesize] = -20;

					DtktkResmass[typesize] = DInfo->tktkRes_mass[j];
					DtktkRespt[typesize] = DInfo->tktkRes_pt[j];
					DtktkReseta[typesize] = DInfo->tktkRes_eta[j];
					DtktkResphi[typesize] = DInfo->tktkRes_phi[j];

					DtktkRes_chi2ndf[typesize] = DInfo->tktkRes_vtxchi2[j]/DInfo->tktkRes_vtxdof[j];
					DtktkRes_chi2cl[typesize] = TMath::Prob(DInfo->tktkRes_vtxchi2[j], DInfo->tktkRes_vtxdof[j]);
					DtktkRes_alpha[typesize] = DInfo->tktkRes_alpha[j];
					TVector3 *DisSvResVtx = new TVector3;
					DisSvResVtx->SetXYZ(DInfo->tktkRes_vtxX[j]-DInfo->vtxX[j],
							DInfo->tktkRes_vtxY[j]-DInfo->vtxY[j],
							DInfo->tktkRes_vtxZ[j]-DInfo->vtxZ[j]);
					TLorentzVector *tktkRes4Vec = new TLorentzVector;
					tktkRes4Vec->SetPtEtaPhiM(DInfo->tktkRes_pt[j], DInfo->tktkRes_eta[j], DInfo->tktkRes_phi[j], DInfo->tktkRes_mass[j]);
					DtktkRes_alphaToSV[typesize] = tktkRes4Vec->Angle(*DisSvResVtx);
					//DtktkRes_alphaToSV[typesize] = DInfo->tktkRes_alphaToSV[j]; // update after moving to new Dfinder
					DtktkRes_svpvDistance[typesize] = DInfo->tktkRes_svpvDistance[j];
					DtktkRes_svpvDisErr[typesize] = DInfo->tktkRes_svpvDisErr[j];
					DtktkRes_svpvDistanceToSV[typesize] = DisSvResVtx->Mag();
					//DtktkRes_svpvDistanceToSV[typesize] = DInfo->tktkRes_svpvDistanceToSV[j]; // update after moving to new Dfinder
					DtktkRes_svpvDisErrToSV[typesize] = 1;
					//DtktkRes_svpvDisErrToSV[typesize] = DInfo->tktkRes_svpvDisErrToSV[j]; // update after moving to new Dfinder
					DtktkRes_dca[typesize] = DInfo->tktkRes_svpvDistance[j]*TMath::Sin(DInfo->tktkRes_alpha[j]);
					DtktkRes_dcaToSV[typesize] = DtktkRes_svpvDistanceToSV[typesize]*TMath::Sin(DtktkRes_alphaToSV[typesize]);
					//DtktkRes_dcaToSV[typesize] = DInfo->tktkRes_svpvDistanceToSV[j]*TMath::Sin(DInfo->tktkRes_alphaToSV[j]); // update after moving to new Dfinder

					r2lxyBS = (DInfo->tktkRes_vtxX[j]-EvtInfo->BSx+(DInfo->tktkRes_vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdxdz) * (DInfo->tktkRes_vtxX[j]-EvtInfo->BSx+(DInfo->tktkRes_vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdxdz)
						+ (DInfo->tktkRes_vtxY[j]-EvtInfo->BSy+(DInfo->tktkRes_vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdydz) * (DInfo->tktkRes_vtxY[j]-EvtInfo->BSy+(DInfo->tktkRes_vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdydz);
					xlxyBS = DInfo->tktkRes_vtxX[j]-EvtInfo->BSx + (DInfo->tktkRes_vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdxdz;
					ylxyBS = DInfo->tktkRes_vtxY[j]-EvtInfo->BSy + (DInfo->tktkRes_vtxZ[j]-EvtInfo->BSz)*EvtInfo->BSdydz;
					DtktkRes_lxyBS[typesize] = TMath::Sqrt(r2lxyBS);
					DtktkRes_lxyBSErr[typesize] = (1./r2lxyBS) * ((xlxyBS*xlxyBS)*DInfo->tktkRes_vtxXErr[j] + (2*xlxyBS*ylxyBS)*DInfo->tktkRes_vtxYXErr[j] + (ylxyBS*ylxyBS)*DInfo->tktkRes_vtxYErr[j]);

					TLorentzVector *trk14Vec = new TLorentzVector;
					trk14Vec->SetPtEtaPhiM(TrackInfo->pt[DInfo->rftk2_index[j]], TrackInfo->eta[DInfo->rftk2_index[j]], TrackInfo->phi[DInfo->rftk2_index[j]], PION_MASS);
					DtktkRes_angleToTrk1[typesize] = tktkRes4Vec->Angle(trk14Vec->Vect());

					DtktkRes_ptAsymToTrk1[typesize] = (DInfo->tktkRes_pt[j]-TrackInfo->pt[DInfo->rftk2_index[j]])/(DInfo->tktkRes_pt[j]+TrackInfo->pt[DInfo->rftk2_index[j]]);
					DtktkRes_unfitter_ptAsymToTrk1[typesize] = (DInfo->tktkRes_unfitted_pt[j]-TrackInfo->pt[DInfo->rftk2_index[j]])/(DInfo->tktkRes_unfitted_pt[j]+TrackInfo->pt[DInfo->rftk2_index[j]]);

					if(DInfo->type[j]==11||DInfo->type[j]==12)
					{
						DRestrk3Pt[typesize] = TrackInfo->pt[DInfo->tktkRes_rftk3_index[j]];
						DRestrk3Eta[typesize] = TrackInfo->eta[DInfo->tktkRes_rftk3_index[j]];
						DRestrk3Phi[typesize] = TrackInfo->phi[DInfo->tktkRes_rftk3_index[j]];
						b4P->SetPtEtaPhiM(TrackInfo->pt[DInfo->tktkRes_rftk3_index[j]],TrackInfo->eta[DInfo->tktkRes_rftk3_index[j]],TrackInfo->phi[DInfo->tktkRes_rftk3_index[j]],DInfo->tktkRes_rftk3_mass[j]);
						DRestrk3Y[typesize] = b4P->Rapidity();
						DRestrk3P[typesize] = b4P->P();
						b4P->Boost(-*boost);
						DRestrk3thetastar[typesize] = b4P->Angle(*D3Vec);
						DRestrk3Dz[typesize] = TrackInfo->dzPV[DInfo->tktkRes_rftk3_index[j]];
						DRestrk3Dxy[typesize] = TrackInfo->dxyPV[DInfo->tktkRes_rftk3_index[j]];
						DRestrk3D0[typesize] = TrackInfo->d0[DInfo->tktkRes_rftk3_index[j]];
						DRestrk3D0Err[typesize] = TrackInfo->d0error[DInfo->tktkRes_rftk3_index[j]];
						DRestrk3originalAlgo[typesize] = TrackInfo->originalTrkAlgo[DInfo->tktkRes_rftk3_index[j]];
						DRestrk4Pt[typesize] = TrackInfo->pt[DInfo->tktkRes_rftk4_index[j]];
						DRestrk4Eta[typesize] = TrackInfo->eta[DInfo->tktkRes_rftk4_index[j]];
						DRestrk4Phi[typesize] = TrackInfo->phi[DInfo->tktkRes_rftk4_index[j]];
						b4P->SetPtEtaPhiM(TrackInfo->pt[DInfo->tktkRes_rftk4_index[j]],TrackInfo->eta[DInfo->tktkRes_rftk4_index[j]],TrackInfo->phi[DInfo->tktkRes_rftk4_index[j]],DInfo->tktkRes_rftk4_mass[j]);
						DRestrk4Y[typesize] = b4P->Rapidity();
						DRestrk4P[typesize] = b4P->P();
						b4P->Boost(-*boost);
						DRestrk4thetastar[typesize] = b4P->Angle(*D3Vec);
						DRestrk4Dz[typesize] = TrackInfo->dzPV[DInfo->tktkRes_rftk4_index[j]];
						DRestrk4Dxy[typesize] = TrackInfo->dxyPV[DInfo->tktkRes_rftk4_index[j]];
						DRestrk4D0[typesize] = TrackInfo->d0[DInfo->tktkRes_rftk4_index[j]];
						DRestrk4D0Err[typesize] = TrackInfo->d0error[DInfo->tktkRes_rftk4_index[j]];
						DRestrk4originalAlgo[typesize] = TrackInfo->originalTrkAlgo[DInfo->tktkRes_rftk4_index[j]];
						DRestrk3PtErr[typesize] = TrackInfo->ptErr[DInfo->tktkRes_rftk3_index[j]];
						DRestrk4PtErr[typesize] = TrackInfo->ptErr[DInfo->tktkRes_rftk4_index[j]];
						DRestrk3highPurity[typesize] = TrackInfo->highPurity[DInfo->tktkRes_rftk3_index[j]];
						DRestrk4highPurity[typesize] = TrackInfo->highPurity[DInfo->tktkRes_rftk4_index[j]];
						DRestrk3Quality[typesize] = TrackInfo->trackQuality[DInfo->tktkRes_rftk3_index[j]];
						DRestrk4Quality[typesize] = TrackInfo->trackQuality[DInfo->tktkRes_rftk4_index[j]];
						DRestrk3dedx[typesize] = TrackInfo->dedx[DInfo->tktkRes_rftk3_index[j]];
						DRestrk4dedx[typesize] = TrackInfo->dedx[DInfo->tktkRes_rftk4_index[j]];
					}
				}
				DMaxTkPt[typesize] = max(Dtrk1Pt[typesize], max(Dtrk2Pt[typesize], max(Dtrk3Pt[typesize], max(Dtrk4Pt[typesize], max(DRestrk1Pt[typesize], max(DRestrk2Pt[typesize], max(DRestrk3Pt[typesize], DRestrk4Pt[typesize])))))));
				DMinTkPt[typesize] = max(1/Dtrk1Pt[typesize], max(1/Dtrk2Pt[typesize], max(1/Dtrk3Pt[typesize], max(1/Dtrk4Pt[typesize], max(1/DRestrk1Pt[typesize], max(1/DRestrk2Pt[typesize], max(1/DRestrk3Pt[typesize], 1/DRestrk4Pt[typesize])))))));
				DMinTkPt[typesize] = 1/DMinTkPt[typesize];

				int DpdgId=0,RpdgId=0;
				int dGenIdxRes = -1;
				if(DInfo->type[j]==1||DInfo->type[j]==2||DInfo->type[j]==5||DInfo->type[j]==6) DpdgId=DZERO_PDGID;
				else if(DInfo->type[j]==3||DInfo->type[j]==4) DpdgId=DPLUS_PDGID;
				else if(DInfo->type[j]==7||DInfo->type[j]==8) DpdgId=DSUBS_PDGID;
				else if(DInfo->type[j]==9||DInfo->type[j]==10||DInfo->type[j]==11||DInfo->type[j]==12) DpdgId=DSTAR_PDGID;
				else if(DInfo->type[j]==13||DInfo->type[j]==14) DpdgId=BPLUS_PDGID;
				if(DInfo->type[j]==7||DInfo->type[j]==8) RpdgId=PHI_PDGID;
				else if(DInfo->type[j]==9||DInfo->type[j]==10||DInfo->type[j]==11||DInfo->type[j]==12||DInfo->type[j]==13||DInfo->type[j]==14) RpdgId=DZERO_PDGID;
				Dgen[typesize] = 0;//gen init
				DgenIndex[typesize] = -1;
				DgennDa[typesize] = -1;
				Dgenpt[typesize] = -1;
				Dgeneta[typesize] = -20;
				Dgenphi[typesize] = -20;
				Dgeny[typesize] = -1;
				DgenBAncestorpt[typesize] = -99;
				DgenBAncestorpdgId[typesize] = 0;
				DgencollisionId[typesize] = -99;
				DgenprodvtxX[typesize] = -999;
				DgenprodvtxY[typesize] = -999;
				DgenprodvtxZ[typesize] = -999;
				DgendecayvtxX[typesize] = -999;
				DgendecayvtxY[typesize] = -999;
				DgendecayvtxZ[typesize] = -999;
				DgenfromgenPV[typesize] = -999;
				if(!REAL)
				{
					if(DInfo->type[j]==1||DInfo->type[j]==2)
					{
						if(DInfo->rftk1_index[j]>-1 && DInfo->rftk2_index[j]>-1)
						{
							if(TrackInfo->geninfo_index[DInfo->rftk1_index[j]]>-1 && 
									TrackInfo->geninfo_index[DInfo->rftk2_index[j]]>-1)
							{
								if(GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]]>-1 && 
										GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]]>-1)
								{
									if(TMath::Abs(GenInfo->pdgId[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]]])==DpdgId && 
											GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]]==GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])
									{
										if(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==DInfo->rftk1_MassHypo[j] && 
												TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==DInfo->rftk2_MassHypo[j])
										{
											Dgen[typesize] = 23333;
											dGenIdxRes = GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]];
										}
										else if(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==DInfo->rftk1_MassHypo[j] && 
												TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==DInfo->rftk2_MassHypo[j])
										{
											Dgen[typesize] = 23344;
											dGenIdxRes = GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]];
										}
									}
								}
							}
						}
					}
					else if(DInfo->type[j]==3||DInfo->type[j]==4)
					{
						if(DInfo->rftk1_index[j]>-1 && DInfo->rftk2_index[j]>-1 && DInfo->rftk3_index[j]>-1)
						{
							if(TrackInfo->geninfo_index[DInfo->rftk1_index[j]]>-1 && 
									TrackInfo->geninfo_index[DInfo->rftk2_index[j]]>-1 && 
									TrackInfo->geninfo_index[DInfo->rftk3_index[j]]>-1)
							{
								if(GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]]>-1 && 
										GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]]>-1 && 
										GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]]>-1)
								{
									if(TMath::Abs(GenInfo->pdgId[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]]])==DpdgId && 
											GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]]==GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]] &&
											GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]]==GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])
									{
										if(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==DInfo->rftk1_MassHypo[j] && 
												TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==DInfo->rftk2_MassHypo[j] &&
												TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])==DInfo->rftk3_MassHypo[j])
										{
											Dgen[typesize] = 23333;
											dGenIdxRes = GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]];
										}
										else if((TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==DInfo->rftk2_MassHypo[j]&&
													TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==DInfo->rftk1_MassHypo[j]&&
													TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])==DInfo->rftk3_MassHypo[j]&&
													TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])==PION_PDGID) ||
												(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==DInfo->rftk3_MassHypo[j]&&
												 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])==DInfo->rftk1_MassHypo[j]&&
												 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==DInfo->rftk2_MassHypo[j]&&
												 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==PION_PDGID) ||
												(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==DInfo->rftk3_MassHypo[j]&&
												 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])==DInfo->rftk2_MassHypo[j]&&
												 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==DInfo->rftk1_MassHypo[j]&&
												 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==PION_PDGID))
										{
											Dgen[typesize] = 23344;
											dGenIdxRes = GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]];
										}
									}
								}
							}
						}
					}
					else if(DInfo->type[j]==5||DInfo->type[j]==6)
					{
						if(DInfo->rftk1_index[j]>-1 && DInfo->rftk2_index[j]>-1 && DInfo->rftk3_index[j]>-1 && DInfo->rftk4_index[j]>-1)
						{
							if(TrackInfo->geninfo_index[DInfo->rftk1_index[j]]>-1 && 
									TrackInfo->geninfo_index[DInfo->rftk2_index[j]]>-1 && 
									TrackInfo->geninfo_index[DInfo->rftk3_index[j]]>-1 && 
									TrackInfo->geninfo_index[DInfo->rftk4_index[j]]>-1)
							{
								if(GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]]>-1 && 
										GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]]>-1 && 
										GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]]>-1 && 
										GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk4_index[j]]]>-1)
								{
									if(TMath::Abs(GenInfo->pdgId[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]]])==DpdgId && 
											GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]]==GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]] &&
											GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]]==GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]] &&
											GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]]==GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk4_index[j]]])
									{
										if(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==DInfo->rftk1_MassHypo[j] && 
												TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==DInfo->rftk2_MassHypo[j] &&
												TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])==DInfo->rftk3_MassHypo[j] &&
												TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk4_index[j]]])==DInfo->rftk4_MassHypo[j])
										{
											Dgen[typesize] = 23333;
											dGenIdxRes = GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]];
										}
										else if((TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==DInfo->rftk2_MassHypo[j] &&
													TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==DInfo->rftk1_MassHypo[j] &&
													TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])==DInfo->rftk3_MassHypo[j] &&
													TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk4_index[j]]])==DInfo->rftk4_MassHypo[j] &&
													TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])==PION_PDGID &&
													TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk4_index[j]]])==PION_PDGID) ||
												(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==DInfo->rftk3_MassHypo[j] &&
												 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])==DInfo->rftk1_MassHypo[j] &&
												 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==DInfo->rftk2_MassHypo[j] &&
												 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk4_index[j]]])==DInfo->rftk4_MassHypo[j] &&
												 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==PION_PDGID &&
												 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk4_index[j]]])==PION_PDGID) ||
												(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==DInfo->rftk4_MassHypo[j] &&
												 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk4_index[j]]])==DInfo->rftk1_MassHypo[j] &&
												 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==DInfo->rftk2_MassHypo[j] &&
												 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])==DInfo->rftk3_MassHypo[j] &&
												 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==PION_PDGID &&
												 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])==PION_PDGID) ||
												(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==DInfo->rftk3_MassHypo[j] &&
												 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])==DInfo->rftk2_MassHypo[j] &&
												 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==DInfo->rftk1_MassHypo[j] &&
												 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk4_index[j]]])==DInfo->rftk4_MassHypo[j] &&
												 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==PION_PDGID &&
												 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk4_index[j]]])==PION_PDGID) ||
												(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==DInfo->rftk4_MassHypo[j] &&
												 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk4_index[j]]])==DInfo->rftk2_MassHypo[j] &&
												 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==DInfo->rftk1_MassHypo[j] &&
												 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])==DInfo->rftk3_MassHypo[j] &&
												 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==PION_PDGID &&
												 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])==PION_PDGID) ||
												(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk3_index[j]]])==DInfo->rftk4_MassHypo[j] &&
												 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk4_index[j]]])==DInfo->rftk3_MassHypo[j] &&
												 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==DInfo->rftk1_MassHypo[j] &&
												 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==DInfo->rftk2_MassHypo[j] &&
												 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]])==PION_PDGID &&
												 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==PION_PDGID))
												 {
													 Dgen[typesize] = 23344;
													 dGenIdxRes = GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk1_index[j]]];
												 }
									}
								}
							}
						}
					}
					else if(DInfo->type[j]==7||DInfo->type[j]==8||DInfo->type[j]==9||DInfo->type[j]==10||DInfo->type[j]==13||DInfo->type[j]==14)
					{
						if(DInfo->tktkRes_rftk1_index[j]>-1 && DInfo->tktkRes_rftk2_index[j]>-1)
						{
							if(TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]>-1 && 
									TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]>-1)
							{
								if(GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]]>-1 && 
										GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]]>-1)
								{
									if(TMath::Abs(GenInfo->pdgId[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]]])==RpdgId && 
											GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]]==GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]])
									{
										if(GenInfo->mo1[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]]]>-1 &&
												GenInfo->mo1[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]]]>-1)
										{
											if(TMath::Abs(GenInfo->pdgId[GenInfo->mo1[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]]]])==DpdgId)
											{
												if(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]])==findPdgid(DInfo->tktkRes_rftk1_mass[j]) && 
														TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]])==findPdgid(DInfo->tktkRes_rftk2_mass[j]))
												{
													Dgen[typesize] = 3333;
													dGenIdxRes = GenInfo->mo1[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]]];
												}
												if((DInfo->type[j]==9||DInfo->type[j]==10||DInfo->type[j]==3||DInfo->type[j]==14) &&
														TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]])==findPdgid(DInfo->tktkRes_rftk1_mass[j]) && 
														TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]])==findPdgid(DInfo->tktkRes_rftk2_mass[j]))
												{
													Dgen[typesize] = 3344;
													dGenIdxRes = GenInfo->mo1[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]]];
												}
											}
										}
									}
								}
							}
						}
						if(DInfo->rftk2_index[j]>-1)
						{
							if(TrackInfo->geninfo_index[DInfo->rftk2_index[j]]>-1)
							{
								if(GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]]>-1)
								{
									if(TMath::Abs(GenInfo->pdgId[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]]])==DpdgId &&
											GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]]==dGenIdxRes)
									{
										if(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==PION_PDGID)
										{
											Dgen[typesize]+=20000;
										}
									}
								}
							}
						}
					}
					else if(DInfo->type[j]==11||DInfo->type[j]==12)
					{
						if(DInfo->tktkRes_rftk1_index[j]>-1 && DInfo->tktkRes_rftk2_index[j]>-1 && DInfo->tktkRes_rftk3_index[j]>-1 && DInfo->tktkRes_rftk4_index[j]>-1)
						{
							if(TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]>-1 && 
									TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]>-1 && 
									TrackInfo->geninfo_index[DInfo->tktkRes_rftk3_index[j]]>-1 && 
									TrackInfo->geninfo_index[DInfo->tktkRes_rftk4_index[j]]>-1)
							{
								if(GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]]>-1 && 
										GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]]>-1 &&
										GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk3_index[j]]]>-1 &&
										GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk4_index[j]]]>-1)
								{
									if(TMath::Abs(GenInfo->pdgId[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]]])==RpdgId && 
											GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]]==GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]] &&
											GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]]==GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk3_index[j]]] &&
											GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]]==GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk4_index[j]]])
									{
										if(GenInfo->mo1[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]]]>-1)
										{
											if(TMath::Abs(GenInfo->pdgId[GenInfo->mo1[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]]]])==DpdgId)
											{
												if(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]])==findPdgid(DInfo->tktkRes_rftk1_mass[j]) && 
														TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]])==findPdgid(DInfo->tktkRes_rftk2_mass[j]) &&
														TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk3_index[j]]])==findPdgid(DInfo->tktkRes_rftk3_mass[j]) &&
														TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk4_index[j]]])==findPdgid(DInfo->tktkRes_rftk4_mass[j]))
												{
													Dgen[typesize] = 3333;
													dGenIdxRes = GenInfo->mo1[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]]];
												}
												else if((TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]])==findPdgid(DInfo->tktkRes_rftk2_mass[j]) && 
															TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]])==findPdgid(DInfo->tktkRes_rftk1_mass[j]) &&
															TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk3_index[j]]])==findPdgid(DInfo->tktkRes_rftk3_mass[j]) &&
															TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk4_index[j]]])==findPdgid(DInfo->tktkRes_rftk4_mass[j]) &&
															TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk3_index[j]]])==PION_PDGID &&
															TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk4_index[j]]])==PION_PDGID) ||
														(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]])==findPdgid(DInfo->tktkRes_rftk3_mass[j]) && 
														 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk3_index[j]]])==findPdgid(DInfo->tktkRes_rftk1_mass[j]) &&
														 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]])==findPdgid(DInfo->tktkRes_rftk2_mass[j]) &&
														 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk4_index[j]]])==findPdgid(DInfo->tktkRes_rftk4_mass[j]) &&
														 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]])==PION_PDGID &&
														 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk4_index[j]]])==PION_PDGID) ||
														(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]])==findPdgid(DInfo->tktkRes_rftk4_mass[j]) && 
														 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk4_index[j]]])==findPdgid(DInfo->tktkRes_rftk1_mass[j]) &&
														 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]])==findPdgid(DInfo->tktkRes_rftk2_mass[j]) &&
														 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk3_index[j]]])==findPdgid(DInfo->tktkRes_rftk3_mass[j]) &&
														 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]])==PION_PDGID &&
														 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk3_index[j]]])==PION_PDGID) ||
														(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]])==findPdgid(DInfo->tktkRes_rftk3_mass[j]) && 
														 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk3_index[j]]])==findPdgid(DInfo->tktkRes_rftk2_mass[j]) &&
														 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]])==findPdgid(DInfo->tktkRes_rftk1_mass[j]) &&
														 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk4_index[j]]])==findPdgid(DInfo->tktkRes_rftk4_mass[j]) &&
														 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]])==PION_PDGID &&
														 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk4_index[j]]])==PION_PDGID) ||
														(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]])==findPdgid(DInfo->tktkRes_rftk4_mass[j]) && 
														 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk4_index[j]]])==findPdgid(DInfo->tktkRes_rftk2_mass[j]) &&
														 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]])==findPdgid(DInfo->tktkRes_rftk1_mass[j]) &&
														 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk3_index[j]]])==findPdgid(DInfo->tktkRes_rftk3_mass[j]) &&
														 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]])==PION_PDGID &&
														 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk3_index[j]]])==PION_PDGID) ||
														(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk3_index[j]]])==findPdgid(DInfo->tktkRes_rftk4_mass[j]) && 
														 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk4_index[j]]])==findPdgid(DInfo->tktkRes_rftk3_mass[j]) &&
														 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]])==findPdgid(DInfo->tktkRes_rftk1_mass[j]) &&
														 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]])==findPdgid(DInfo->tktkRes_rftk2_mass[j]) &&
														 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]])==PION_PDGID &&
														 TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->tktkRes_rftk2_index[j]]])==PION_PDGID))
														 {
															 Dgen[typesize] = 3344;
															 dGenIdxRes = GenInfo->mo1[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->tktkRes_rftk1_index[j]]]];
														 }
											}
										}
									}
								}
							}
						}
						if(DInfo->rftk2_index[j]>-1)
						{
							if(TrackInfo->geninfo_index[DInfo->rftk2_index[j]]>-1)
							{
								if(GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]]>-1)
								{
									if(TMath::Abs(GenInfo->pdgId[GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]]])==DpdgId &&
											GenInfo->mo1[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]]==dGenIdxRes)
									{
										if(TMath::Abs(GenInfo->pdgId[TrackInfo->geninfo_index[DInfo->rftk2_index[j]]])==PION_PDGID)
										{
											Dgen[typesize]+=20000;
										}
									}
								}
							}
						}
					}

					if(Dgen[typesize]==23333||Dgen[typesize]==23344)
					{
						if(dGenIdxRes<0) cout<<"ERROR: Gen-Matched D index is -1"<<endl;
						else
						{
							DgenIndex[typesize] = dGenIdxRes;
							if((DInfo->type[j]==1||DInfo->type[j]==2)&&GenInfo->nDa[DgenIndex[typesize]]>2) Dgen[typesize]=41000;
							DgennDa[typesize] = GenInfo->nDa[DgenIndex[typesize]];
							Dgenpt[typesize] = GenInfo->pt[DgenIndex[typesize]];
							Dgeneta[typesize] = GenInfo->eta[DgenIndex[typesize]];
							Dgenphi[typesize] = GenInfo->phi[DgenIndex[typesize]];
							DgencollisionId[typesize] = GenInfo->collisionId[DgenIndex[typesize]];
							b4P->SetXYZM(GenInfo->pt[DgenIndex[typesize]]*cos(GenInfo->phi[DgenIndex[typesize]]),
									GenInfo->pt[DgenIndex[typesize]]*sin(GenInfo->phi[DgenIndex[typesize]]),
									GenInfo->pt[DgenIndex[typesize]]*sinh(GenInfo->eta[DgenIndex[typesize]]),
									GenInfo->mass[DgenIndex[typesize]]);
							Dgeny[typesize] = b4P->Rapidity();
							DgenprodvtxX[typesize] = GenInfo->vtxX[DgenIndex[typesize]];
							DgenprodvtxY[typesize] = GenInfo->vtxY[DgenIndex[typesize]];
							DgenprodvtxZ[typesize] = GenInfo->vtxZ[DgenIndex[typesize]];
							DgendecayvtxX[typesize] = GenInfo->vtxX[GenInfo->da1[DgenIndex[typesize]]]; //production vertex of first daughter
							DgendecayvtxY[typesize] = GenInfo->vtxY[GenInfo->da1[DgenIndex[typesize]]];
							DgendecayvtxZ[typesize] = GenInfo->vtxZ[GenInfo->da1[DgenIndex[typesize]]];
							//decide if from gen PV or not
							if( fabs(DgenprodvtxX[typesize] - GenInfo->genPVx) < 0.001 && fabs(DgenprodvtxY[typesize] - GenInfo->genPVy) < 0.001 && fabs(DgenprodvtxZ[typesize] - GenInfo->genPVz) < 0.001 )
								DgenfromgenPV[typesize] = 1;
							else
								DgenfromgenPV[typesize] = -1;
							int DgenBAncestorindex = findBAncestor(DgenIndex[typesize], GenInfo);
							if( DgenBAncestorindex >= 0 )
							{
								DgenBAncestorpt[typesize] = GenInfo->pt[DgenBAncestorindex];
								DgenBAncestorpdgId[typesize] = GenInfo->pdgId[DgenBAncestorindex];
							}
						}
					}
				}//if(!real)
			}//fillDtree

			int findBAncestor(int j, GenInfoBranches *GenInfo)
			{
				int BAncestorindex = -999;
				if( GenInfo->nMo[j] == 0 ) return BAncestorindex;

				int motherindex = -999;
				int daughterindex = j;
				int igeneration = 0;//to control how many generations loop up, 50 is big enough. To avoid infinite loop
				//just work for 1 mom case yet, but more than 99.9% particle just have one mom
				//and in B->D decay chain, all particles should just have one mom (they are not from collision). Should be checked
				while( GenInfo->nMo[daughterindex] == 1 && BAncestorindex < 0 && igeneration < 50 )
				{
					motherindex = GenInfo->mo1[daughterindex];
					if( motherindex < 0 ) break;
					if( ( TMath::Abs( GenInfo->pdgId[motherindex] ) > 500  && TMath::Abs( GenInfo->pdgId[motherindex] ) < 600 ) || ( TMath::Abs( GenInfo->pdgId[motherindex] ) > 5000  && TMath::Abs( GenInfo->pdgId[motherindex] ) < 6000 ) )
						BAncestorindex = motherindex;
					igeneration++;
					daughterindex = motherindex;
				}

				return BAncestorindex;
			}

			bool isDsignalGen(int dmesontype, int j, GenInfoBranches *GenInfo)
			{
				bool flag=false;
				if(dmesontype==1||dmesontype==2)
				{
					if(TMath::Abs(GenInfo->pdgId[j])==DZERO_PDGID&&GenInfo->nDa[j]==2&&GenInfo->da1[j]!=-1&&GenInfo->da2[j]!=-1)
					{
						if((((GenInfo->pdgId[GenInfo->da1[j]]==KAON_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[j]])==PION_PDGID)||
										(GenInfo->pdgId[GenInfo->da2[j]]==KAON_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[j]])==PION_PDGID))&&dmesontype==1) ||
								(((GenInfo->pdgId[GenInfo->da1[j]]==(0-KAON_PDGID)&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[j]])==PION_PDGID)||
								  (GenInfo->pdgId[GenInfo->da2[j]]==(0-KAON_PDGID)&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[j]])==PION_PDGID))&&dmesontype==2))
						{
							flag=true;
						}
					}
				}
				if(dmesontype==3||dmesontype==4)
				{
					if(TMath::Abs(GenInfo->pdgId[j])==DPLUS_PDGID&&GenInfo->nDa[j]==3&&GenInfo->da1[j]!=-1&&GenInfo->da2[j]!=-1&&GenInfo->da3[j]!=-1)
					{
						if((((GenInfo->pdgId[GenInfo->da1[j]]==KAON_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da3[j]])==PION_PDGID)||
										(GenInfo->pdgId[GenInfo->da2[j]]==KAON_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da3[j]])==PION_PDGID)||
										(GenInfo->pdgId[GenInfo->da3[j]]==KAON_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[j]])==PION_PDGID))&&dmesontype==4) ||
								(((GenInfo->pdgId[GenInfo->da1[j]]==(0-KAON_PDGID)&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da3[j]])==PION_PDGID)||
								  (GenInfo->pdgId[GenInfo->da2[j]]==(0-KAON_PDGID)&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da3[j]])==PION_PDGID)||
								  (GenInfo->pdgId[GenInfo->da3[j]]==(0-KAON_PDGID)&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[j]])==PION_PDGID))&&dmesontype==3))
						{
							flag=true;
						}
					}
				}
				if(dmesontype==5||dmesontype==6)
				{
					if(TMath::Abs(GenInfo->pdgId[j])==DZERO_PDGID&&GenInfo->nDa[j]==4&&GenInfo->da1[j]!=-1&&GenInfo->da2[j]!=-1&&GenInfo->da3[j]!=-1&&GenInfo->da4[j]!=-1)
					{
						if((((GenInfo->pdgId[GenInfo->da1[j]]==KAON_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da3[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da4[j]])==PION_PDGID)||
										(GenInfo->pdgId[GenInfo->da2[j]]==KAON_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da3[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da4[j]])==PION_PDGID)||
										(GenInfo->pdgId[GenInfo->da3[j]]==KAON_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da4[j]])==PION_PDGID)||
										(GenInfo->pdgId[GenInfo->da4[j]]==KAON_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da3[j]])==PION_PDGID))&&dmesontype==6) ||
								(((GenInfo->pdgId[GenInfo->da1[j]]==(0-KAON_PDGID)&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da3[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da4[j]])==PION_PDGID)||
								  (GenInfo->pdgId[GenInfo->da2[j]]==(0-KAON_PDGID)&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da3[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da4[j]])==PION_PDGID)||
								  (GenInfo->pdgId[GenInfo->da3[j]]==(0-KAON_PDGID)&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da4[j]])==PION_PDGID)||
								  (GenInfo->pdgId[GenInfo->da4[j]]==(0-KAON_PDGID)&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[j]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da3[j]])==PION_PDGID))&&dmesontype==5))
						{
							flag=true;
						}
					}
				}
				if(dmesontype==7||dmesontype==8)
				{
					if(TMath::Abs(GenInfo->pdgId[j])==DSUBS_PDGID&&GenInfo->nDa[j]==2&&GenInfo->da1[j]!=-1&&GenInfo->da2[j]!=-1)
					{
						if(TMath::Abs(GenInfo->pdgId[GenInfo->da1[j]])==PHI_PDGID)
						{
							if(GenInfo->nDa[GenInfo->da1[j]]==2&&GenInfo->da1[GenInfo->da1[j]]!=-1&&GenInfo->da2[GenInfo->da1[j]]!=-1)
							{
								if(TMath::Abs(GenInfo->pdgId[GenInfo->da1[GenInfo->da1[j]]])==KAON_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[GenInfo->da1[j]]])==KAON_PDGID)
								{
									if((GenInfo->pdgId[GenInfo->da2[j]]==PION_PDGID&&dmesontype==7) ||
											(GenInfo->pdgId[GenInfo->da2[j]]==(0-PION_PDGID)&&dmesontype==8))
										flag=true;                      
								}
							}
						}
					}
				}
				if(dmesontype==9||dmesontype==10)
				{
					if(TMath::Abs(GenInfo->pdgId[j])==DSTAR_PDGID&&GenInfo->nDa[j]==2&&GenInfo->da1[j]!=-1&&GenInfo->da2[j]!=-1)
					{
						if(TMath::Abs(GenInfo->pdgId[GenInfo->da1[j]])==DZERO_PDGID)
						{
							if(GenInfo->nDa[GenInfo->da1[j]]==2&&GenInfo->da1[GenInfo->da1[j]]!=-1&&GenInfo->da2[GenInfo->da1[j]]!=-1)
							{
								if((((GenInfo->pdgId[GenInfo->da2[GenInfo->da1[j]]]==KAON_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[GenInfo->da1[j]]])==PION_PDGID)||
												(GenInfo->pdgId[GenInfo->da1[GenInfo->da1[j]]]==KAON_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[GenInfo->da1[j]]])==PION_PDGID))&&GenInfo->pdgId[GenInfo->da2[j]]==(0-PION_PDGID)&&dmesontype==10) ||
										(((GenInfo->pdgId[GenInfo->da2[GenInfo->da1[j]]]==(0-KAON_PDGID)&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[GenInfo->da1[j]]])==PION_PDGID)||
										  (GenInfo->pdgId[GenInfo->da1[GenInfo->da1[j]]]==(0-KAON_PDGID)&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[GenInfo->da1[j]]])==PION_PDGID))&&GenInfo->pdgId[GenInfo->da2[j]]==PION_PDGID&&dmesontype==9))
								{
									flag=true;                      
								}
							}
						}
					}
				}
				if(dmesontype==11||dmesontype==12)
				{
					if(TMath::Abs(GenInfo->pdgId[j])==DSTAR_PDGID&&GenInfo->nDa[j]==2&&GenInfo->da1[j]!=-1&&GenInfo->da2[j]!=-1)
					{
						if(TMath::Abs(GenInfo->pdgId[GenInfo->da1[j]])==DZERO_PDGID)
						{
							if(GenInfo->nDa[GenInfo->da1[j]]==4&&GenInfo->da1[GenInfo->da1[j]]!=-1&&GenInfo->da2[GenInfo->da1[j]]!=-1&&GenInfo->da3[GenInfo->da1[j]]!=-1&&GenInfo->da4[GenInfo->da1[j]]!=-1)
							{
								if((((GenInfo->pdgId[GenInfo->da1[GenInfo->da1[j]]]==KAON_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[GenInfo->da1[j]]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da3[GenInfo->da1[j]]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da4[GenInfo->da1[j]]])==PION_PDGID)||
												(GenInfo->pdgId[GenInfo->da2[GenInfo->da1[j]]]==KAON_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[GenInfo->da1[j]]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da3[GenInfo->da1[j]]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da4[GenInfo->da1[j]]])==PION_PDGID)||
												(GenInfo->pdgId[GenInfo->da3[GenInfo->da1[j]]]==KAON_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[GenInfo->da1[j]]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[GenInfo->da1[j]]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da4[GenInfo->da1[j]]])==PION_PDGID)||
												(GenInfo->pdgId[GenInfo->da4[GenInfo->da1[j]]]==KAON_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[GenInfo->da1[j]]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[GenInfo->da1[j]]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da3[GenInfo->da1[j]]])==PION_PDGID))&&GenInfo->pdgId[GenInfo->da2[j]]==(0-PION_PDGID)&&dmesontype==12) ||
										(((GenInfo->pdgId[GenInfo->da1[GenInfo->da1[j]]]==(0-KAON_PDGID)&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[GenInfo->da1[j]]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da3[GenInfo->da1[j]]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da4[GenInfo->da1[j]]])==PION_PDGID)||
										  (GenInfo->pdgId[GenInfo->da2[GenInfo->da1[j]]]==(0-KAON_PDGID)&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[GenInfo->da1[j]]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da3[GenInfo->da1[j]]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da4[GenInfo->da1[j]]])==PION_PDGID)||
										  (GenInfo->pdgId[GenInfo->da3[GenInfo->da1[j]]]==(0-KAON_PDGID)&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[GenInfo->da1[j]]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[GenInfo->da1[j]]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da4[GenInfo->da1[j]]])==PION_PDGID)||
										  (GenInfo->pdgId[GenInfo->da4[GenInfo->da1[j]]]==(0-KAON_PDGID)&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[GenInfo->da1[j]]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[GenInfo->da1[j]]])==PION_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da3[GenInfo->da1[j]]])==PION_PDGID))&&GenInfo->pdgId[GenInfo->da2[j]]==PION_PDGID&&dmesontype==11))
								{
									flag=true;                      
								}
							}
						}
					}
				}
				if(dmesontype==13||dmesontype==14)
				{
					if(TMath::Abs(GenInfo->pdgId[j])==BPLUS_PDGID&&GenInfo->nDa[j]==2&&GenInfo->da1[j]!=-1&&GenInfo->da2[j]!=-1)
					{
						if(TMath::Abs(GenInfo->pdgId[GenInfo->da1[j]])==DZERO_PDGID)
						{
							if(GenInfo->nDa[GenInfo->da1[j]]==2&&GenInfo->da1[GenInfo->da1[j]]!=-1&&GenInfo->da2[GenInfo->da1[j]]!=-1)
							{
								if((((GenInfo->pdgId[GenInfo->da2[GenInfo->da1[j]]]==KAON_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[GenInfo->da1[j]]])==PION_PDGID)||
												(GenInfo->pdgId[GenInfo->da1[GenInfo->da1[j]]]==KAON_PDGID&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[GenInfo->da1[j]]])==PION_PDGID))&&GenInfo->pdgId[GenInfo->da2[j]]==PION_PDGID&&dmesontype==13) ||
										(((GenInfo->pdgId[GenInfo->da2[GenInfo->da1[j]]]==(0-KAON_PDGID)&&TMath::Abs(GenInfo->pdgId[GenInfo->da1[GenInfo->da1[j]]])==PION_PDGID)||
										  (GenInfo->pdgId[GenInfo->da1[GenInfo->da1[j]]]==(0-KAON_PDGID)&&TMath::Abs(GenInfo->pdgId[GenInfo->da2[GenInfo->da1[j]]])==PION_PDGID))&&GenInfo->pdgId[GenInfo->da2[j]]==(0-PION_PDGID)&&dmesontype==14))
								{
									flag=true;                      
								}
							}
						}
					}
				}

				return flag;
			}

		};

#endif
