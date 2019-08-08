#include <iostream>
#include "parameters.h"
#include "branchStruct.h"
#include "tnp_weight_lowptPbPb.h"
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>
using namespace std;

bool ispp = 0;
TString inputmc;

int _nBins = nBins;
double *_ptBins = ptBins;

void converter(){
	string label = "";
	if(ispp){
		label = "pp"; 
		inputmc = "/afs/cern.ch/work/g/gwkim/private/samples/Bntuple20160816_Bpt7svpv5p5Bpt10svpv3p5_BfinderMC_pp_Pythia8_BuToJpsiK_TuneCUEP8M1_20160816_bPt0jpsiPt0tkPt0p5_Bp_pthatweight_JingBDT.root";
    }
	else{
		label = "pbpb";
		inputmc = "/afs/cern.ch/work/g/gwkim/private/samples/crab_Bfinder_20190624_Hydjet_Pythia8_Official_BuToJpsiK_1033p1_pt3tkpt0p7dls2_allpthat_pthatweight.root";
	}
    TFile* infMC = new TFile(inputmc.Data());
    TTree* ntKp = (TTree*)infMC->Get("Bfinder/ntKp");
    TTree* ntHlt = (TTree*)infMC->Get("hltanalysis/HltTree");
    TTree* ntSkim = (TTree*)infMC->Get("skimanalysis/HltTree");
    TTree* ntHi = (TTree*)infMC->Get("hiEvtAnalyzer/HiTree");
    TTree* ntGen = (TTree*)infMC->Get("Bfinder/ntGen");
    ntKp->AddFriend("hltanalysis/HltTree");
    ntKp->AddFriend("hiEvtAnalyzer/HiTree");
    ntKp->AddFriend("skimanalysis/HltTree");
    setAddressTree(ntKp, ntHlt, ntSkim, ntHi, ntGen, ispp);
    int nevents_total = ntKp->GetEntries();

    TFile *f= new TFile(Form("results_%s.root", label.c_str()), "recreate"); 
    TH1D* hGen = new TH1D("hGen","",_nBins,_ptBins);
    TH1D* hNominal = new TH1D("hNominal","",_nBins,_ptBins);
    TH1D* hScale = new TH1D("hScale","",_nBins,_ptBins);
    TH1D* hStatTrg[100];
    TH1D* hStatMuid[100];
    //TH1D* hStatSTA[100];
        for(int i = 1; i <= 100; i++){
		hStatTrg[i-1] = new TH1D(Form("hStatTrg_%d", i),"",_nBins,_ptBins);
		hStatMuid[i-1] = new TH1D(Form("hStatMuid_%d", i),"",_nBins,_ptBins);
		//hStatSTA[i-1] = new TH1D(Form("hStatSTA_%d", i),"",_nBins,_ptBins);
	}
	TH1D* hSysHiTrg  = new TH1D("hSysHiTrg","",_nBins,_ptBins);
	TH1D* hSysHiTrk  = new TH1D("hSysHiTrk","",_nBins,_ptBins);
	TH1D* hSysHiMuid = new TH1D("hSysHiMuid","",_nBins,_ptBins);
	//TH1D* hSysHiSTA  = new TH1D("hSysHiSTA","",_nBins,_ptBins);
	TH1D* hSysLoTrg  = new TH1D("hSysLoTrg","",_nBins,_ptBins);
	TH1D* hSysLoTrk  = new TH1D("hSysLoTrk","",_nBins,_ptBins);
	TH1D* hSysLoMuid = new TH1D("hSysLoMuid","",_nBins,_ptBins);
	//TH1D* hSysLoSTA  = new TH1D("hSysLoSTA","",_nBins,_ptBins);
	TH1D* hSysBinnedTrg  = new TH1D("hSysBinnedTrg","",_nBins,_ptBins);
	TH1D* hStatUnc = new TH1D("hStatUnc","",_nBins,_ptBins);
	TH1D* hSysUnc = new TH1D("hSysUnc","",_nBins,_ptBins);
	TH1D* hAllUnc = new TH1D("hAllUnc","",_nBins,_ptBins);
	TH1D* hSF = new TH1D("hSF","",_nBins,_ptBins);
    for(int entry=0; entry<nevents_total; entry++){
//    for(int entry=111; entry<112; entry++){
	    //if ((entry%10000) == 0) printf("Loading event #%d of %d.\n",entry,nevents_total);
    	ntKp->GetEntry(entry);
    	ntHlt->GetEntry(entry);
    	ntSkim->GetEntry(entry);
    	ntHi->GetEntry(entry);
    	ntGen->GetEntry(entry);
	    for(int g=0; g<Gsize; g++){
        }
	    for(int b=0; b<Bsize; b++){
	                 if(ispp){
				if(HLT_HIL1DoubleMu0ForPPRef_v1)
				if(abs(PVz)<15
				&&pBeamScrapingFilter
				&&pPAprimaryVertexFilter
				&&TMath::Abs(By[b])<2.4
				&&TMath::Abs(Bmumumass[b]-3.096916)<0.15
				&&Bmass[b]>5&&Bmass[b]<6
				&& ((abs(Bmu1eta[b])<1.2 && Bmu1pt[b]>3.5) || (abs(Bmu1eta[b])>1.2 && abs(Bmu1eta[b])<2.1 && Bmu1pt[b]>(5.77-1.8*abs(Bmu1eta[b]))) || (abs(Bmu1eta[b])>2.1 && abs(Bmu1eta[b])<2.4 && Bmu1pt[b]>1.8)) 
				&& ((abs(Bmu2eta[b])<1.2 && Bmu2pt[b]>3.5) || (abs(Bmu2eta[b])>1.2 && abs(Bmu2eta[b])<2.1 && Bmu2pt[b]>(5.77-1.8*abs(Bmu2eta[b]))) || (abs(Bmu2eta[b])>2.1 && abs(Bmu2eta[b])<2.4 && Bmu2pt[b]>1.8)) 
				&& Bmu1TMOneStationTight[b] && Bmu2TMOneStationTight[b] 
				&& Bmu1InPixelLayer[b] > 0 && (Bmu1InPixelLayer[b]+Bmu1InStripLayer[b]) > 5 && Bmu2InPixelLayer[b]> 0 && (Bmu2InPixelLayer[b]+Bmu2InStripLayer[b]) > 5 
				&& Bmu1dxyPV[b]< 0.3 && Bmu2dxyPV[b]< 0.3 && Bmu1dzPV[b]<20 && Bmu2dzPV[b]<20 
				&& Bmu1isGlobalMuon[b] && Bmu2isGlobalMuon[b] 
				&& Bmu1TrgMatchFilterE[b]>0 && Bmu2TrgMatchFilterE[b]>0 
				&& Btrk1highPurity[b] && abs(Btrk1Eta[b])<2.4 && Btrk1Pt[b]>0.5 && Bchi2cl[b]>0.005 
				&& ((Bpt[b]<10 && (BsvpvDistance[b]/BsvpvDisErr[b])>5.5) || (Bpt[b]>10 && (BsvpvDistance[b]/BsvpvDisErr[b])>3.5)) 
				&& ((Bpt[b]>7 && Bpt[b]<10 && Btrk1Pt[b]>0.748 && Bchi2cl[b]>0.0322 && (BsvpvDistance[b]/BsvpvDisErr[b])>0 && cos(Bdtheta[b])>0.989 && abs(Btrk1Eta[b])<2.4) 
				|| (Bpt[b]>10 && Bpt[b]<15 && Btrk1Pt[b]>0.878 && Bchi2cl[b]>0.0014 && (BsvpvDistance[b]/BsvpvDisErr[b])>0 && cos(Bdtheta[b])>0.971 && abs(Btrk1Eta[b])<2.39) 
				|| (Bpt[b]>15 && Bpt[b]<20 && Btrk1Pt[b]>0.818 && Bchi2cl[b]>0.00838 && (BsvpvDistance[b]/BsvpvDisErr[b])>0 && cos(Bdtheta[b])>0.959 && abs(Btrk1Eta[b])<2.42) 
				|| (Bpt[b]>20 && Bpt[b]<30 && Btrk1Pt[b]>0.84 && Bchi2cl[b]>0.014 && (BsvpvDistance[b]/BsvpvDisErr[b])>0 && cos(Bdtheta[b])>0.603 && abs(Btrk1Eta[b])<2.38) 
				|| (Bpt[b]>30 && Bpt[b]<50 && Btrk1Pt[b]>1.06 && Bchi2cl[b]>0.0155 && (BsvpvDistance[b]/BsvpvDisErr[b])>3.53 && cos(Bdtheta[b])>0.995 && abs(Btrk1Eta[b])<2.37))
				) 
				if(Bgen[b]==23333){
					double _weight = pthatweight*((pow(10,-0.094152+0.008102*Bgenpt[b]+Bgenpt[b]*Bgenpt[b]*0.000171+Bgenpt[b]*Bgenpt[b]*Bgenpt[b]*-0.000005+Bgenpt[b]*Bgenpt[b]*Bgenpt[b]*Bgenpt[b]*-0.000000+Bgenpt[b]*Bgenpt[b]*Bgenpt[b]*Bgenpt[b]*Bgenpt[b]*0.000000)));
					hNominal->Fill(Bpt[b], _weight);
					////passed and fill
					//double scale = tnp_weight_trg_pp(Bmu1pt[b], Bmu1eta[b], 0)*tnp_weight_trk_pp(0)*tnp_weight_trg_pp(Bmu2pt[b], Bmu2eta[b], 0)*tnp_weight_trk_pp(0);
					//hScale->Fill(Bpt[b], _weight*scale);
					//for(int i = 1; i <= 100; i ++){
					//double statTrg = tnp_weight_trg_pp(Bmu1pt[b], Bmu1eta[b], i)*tnp_weight_trk_pp(0)*tnp_weight_trg_pp(Bmu2pt[b], Bmu2eta[b], i)*tnp_weight_trk_pp(0);
					//double statMuid = scale*tnp_weight_muid_pp(Bmu1pt[b], Bmu1eta[b], i)*tnp_weight_muid_pp(Bmu2pt[b], Bmu2eta[b], i);
					//double statSTA = scale*tnp_weight_sta_pp(Bmu1pt[b], Bmu1eta[b], i)*tnp_weight_sta_pp(Bmu2pt[b], Bmu2eta[b], i);
					//hStatTrg[i-1]->Fill(Bpt[b], _weight*statTrg);
					//hStatMuid[i-1]->Fill(Bpt[b], _weight*statMuid);
					//hStatSTA[i-1]->Fill(Bpt[b], _weight*statSTA);
					//}
				//double sysHiTrg = tnp_weight_trg_pp(Bmu1pt[b], Bmu1eta[b], -1)*tnp_weight_trk_pp(0)*tnp_weight_trg_pp(Bmu2pt[b], Bmu2eta[b], -1)*tnp_weight_trk_pp(0);
				//double sysLoTrg = tnp_weight_trg_pp(Bmu1pt[b], Bmu1eta[b], -2)*tnp_weight_trk_pp(0)*tnp_weight_trg_pp(Bmu2pt[b], Bmu2eta[b], -2)*tnp_weight_trk_pp(0);
				//double sysHiTrk = tnp_weight_trg_pp(Bmu1pt[b], Bmu1eta[b], 0)*tnp_weight_trk_pp(-1)*tnp_weight_trg_pp(Bmu2pt[b], Bmu2eta[b], 0)*tnp_weight_trk_pp(-1);
				//double sysLoTrk = tnp_weight_trg_pp(Bmu1pt[b], Bmu1eta[b], 0)*tnp_weight_trk_pp(-2)*tnp_weight_trg_pp(Bmu2pt[b], Bmu2eta[b], 0)*tnp_weight_trk_pp(-2);
				//double sysHiMuid = scale*tnp_weight_muid_pp(Bmu1pt[b], Bmu1eta[b], -1)*tnp_weight_muid_pp(Bmu2pt[b], Bmu2eta[b], -1);
				//double sysLoMuid = scale*tnp_weight_muid_pp(Bmu1pt[b], Bmu1eta[b], -2)*tnp_weight_muid_pp(Bmu2pt[b], Bmu2eta[b], -2);
				//double sysHiSTA = scale*tnp_weight_sta_pp(Bmu1pt[b], Bmu1eta[b], -1)*tnp_weight_sta_pp(Bmu2pt[b], Bmu2eta[b], -1);
				//double sysLoSTA = scale*tnp_weight_sta_pp(Bmu1pt[b], Bmu1eta[b], -2)*tnp_weight_sta_pp(Bmu2pt[b], Bmu2eta[b], -2);
				//double sysbinnedTrg = tnp_weight_trg_pp(Bmu1pt[b], Bmu1eta[b], -10)*tnp_weight_trk_pp(0)*tnp_weight_trg_pp(Bmu2pt[b], Bmu2eta[b], -10)*tnp_weight_trk_pp(0);
				//hSysHiTrg ->Fill(Bpt[b], _weight*sysHiTrg);
				//hSysHiTrk ->Fill(Bpt[b], _weight*sysHiTrk);
				//hSysHiMuid->Fill(Bpt[b], _weight*sysHiMuid);
				//hSysHiSTA ->Fill(Bpt[b], _weight*sysHiSTA);
				//hSysLoTrg ->Fill(Bpt[b], _weight*sysLoTrg);
				//hSysLoTrk ->Fill(Bpt[b], _weight*sysLoTrk);
				//hSysLoMuid->Fill(Bpt[b], _weight*sysLoMuid);
				//hSysLoSTA ->Fill(Bpt[b], _weight*sysLoSTA);
				//hSysBinnedTrg ->Fill(Bpt[b], _weight*sysbinnedTrg);
				}
	    }
			else{
				if(HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1)
				if(pclusterCompatibilityFilter
				&&pprimaryVertexFilter
				&&phfCoincFilter2Th4
				&&abs(PVz)<15
				&&TMath::Abs(By[b])<2.4
				&&TMath::Abs(Bmumumass[b]-3.096900)<0.15
				&&Bmass[b]>5&&Bmass[b]<6
				&& ((abs(Bmu1eta[b])<1.2 && Bmu1pt[b]>=3.5) || (abs(Bmu1eta[b])>=1.2 && abs(Bmu1eta[b])<2.1 && Bmu1pt[b]>=(5.77-1.89*abs(Bmu1eta[b]))) || (abs(Bmu1eta[b])>=2.1 && Bmu1pt[b]>=1.8)) 
				&& ((abs(Bmu2eta[b])<1.2 && Bmu2pt[b]>=3.5) || (abs(Bmu2eta[b])>=1.2 && abs(Bmu2eta[b])<2.1 && Bmu2pt[b]>=(5.77-1.89*abs(Bmu2eta[b]))) || (abs(Bmu2eta[b])>=2.1 && Bmu2pt[b]>=1.8)) 
				   && Bmu1SoftMuID[b] && Bmu2SoftMuID[b] && Bmu1isAcc[b] && Bmu2isAcc[b] && Bmu1isTriggered[b] && Bmu2isTriggered[b] && (Btrk1PixelHit[b]+Btrk1StripHit[b])>=11 && (Btrk1Chi2ndf[b]/(Btrk1nStripLayer[b]+Btrk1nPixelLayer[b]))<0.18 && TMath::Abs(Btrk1PtErr[b]/Btrk1Pt[b])<0.1
				&& (BsvpvDistance[b]/BsvpvDisErr[b])>2.0 && abs(Btrk1Eta[b])<2.4 && Btrk1Pt[b]>0.9 && Bchi2cl[b]>0.005 
				&& ((Bpt[b]>5 && Bpt[b]<7 && (BsvpvDistance[b]/BsvpvDisErr[b])>13.028 && cos(Bdtheta[b])>-0.785 && TMath::Abs(Btrk1Dxy1[b]/Btrk1DxyError1[b])>3.664 && Btrk1Pt[b]>1.104 && Bchi2cl[b]>0.224) || (Bpt[b]>7 && Bpt[b]<10 && (BsvpvDistance[b]/BsvpvDisErr[b])>6.151 && cos(Bdtheta[b])>-0.279 && TMath::Abs(Btrk1Dxy1[b]/Btrk1DxyError1[b])>3.379 && TMath::Abs(Btrk1Dz1[b]/Btrk1DzError1[b])>0.628 && Btrk1Pt[b]>1.185) || (Bpt[b]>10 && Bpt[b]<15 && (BsvpvDistance[b]/BsvpvDisErr[b])>9.641 && cos(Bdtheta[b])>-0.510 && TMath::Abs(Btrk1Dxy1[b]/Btrk1DxyError1[b])>3.464 && Btrk1Pt[b]>1.287 && Bchi2cl[b]>0.185) || (Bpt[b]>15 && Bpt[b]<20 && (BsvpvDistance[b]/BsvpvDisErr[b])>6.520 && cos(Bdtheta[b])>0.971 && Btrk1Pt[b]>1.837 && Bchi2cl[b]>0.089) || (Bpt[b]>20 && Bpt[b]<30 && (BsvpvDistance[b]/BsvpvDisErr[b])>4.171 && cos(Bdtheta[b])>0.998 && Btrk1Pt[b]>1.692) || (Bpt[b]>30 && Bpt[b]<50 && (BsvpvDistance[b]/BsvpvDisErr[b])>3.850 && cos(Bdtheta[b])>0.571 && Btrk1Pt[b]>1.723) || (Bpt[b]>50 && Bpt[b]<100 && cos(Bdtheta[b])>0.743 && TMath::Abs(Btrk1Dxy1[b]/Btrk1DxyError1[b])>0.203 && Btrk1Pt[b]>2.973 && Bchi2cl[b]>0.062))
				)	
                if(Bgen[b]==23333){
		  double _weight = pthatweight*Ncoll*(1.034350*TMath::Exp(-0.000844*(PVz+3.502992)*(PVz+3.502992)))*(0.715021+0.039896*Bgenpt[b]-0.000834*Bgenpt[b]*Bgenpt[b]+0.000006*Bgenpt[b]*Bgenpt[b]*Bgenpt[b]);
					hNominal->Fill(Bpt[b], _weight);
					//passed and fill
					double scale = tnp_weight_trg_pbpb(Bmu1pt[b], Bmu1eta[b], 0)*tnp_weight_trk_pbpb(0)*tnp_weight_trg_pbpb(Bmu2pt[b], Bmu2eta[b], 0)*tnp_weight_trk_pbpb(0);
					hScale->Fill(Bpt[b], _weight*scale);
					for(int i = 1; i <= 100; i ++){
						double statTrg = tnp_weight_trg_pbpb(Bmu1pt[b], Bmu1eta[b], i)*tnp_weight_trk_pbpb(0)*tnp_weight_trg_pbpb(Bmu2pt[b], Bmu2eta[b], i)*tnp_weight_trk_pbpb(0);
						double statMuid = scale*tnp_weight_muid_pbpb(Bmu1pt[b], Bmu1eta[b], i)*tnp_weight_muid_pbpb(Bmu2pt[b], Bmu2eta[b], i);
						//double statSTA = scale*tnp_weight_sta_pbpb(Bmu1pt[b], Bmu1eta[b], i)*tnp_weight_sta_pbpb(Bmu2pt[b], Bmu2eta[b], i);
						hStatTrg[i-1]->Fill(Bpt[b], _weight*statTrg);
						hStatMuid[i-1]->Fill(Bpt[b], _weight*statMuid);
						//hStatSTA[i-1]->Fill(Bpt[b], _weight*statSTA);
					}
					double sysHiTrg = tnp_weight_trg_pbpb(Bmu1pt[b], Bmu1eta[b], -1)*tnp_weight_trk_pbpb(0)*tnp_weight_trg_pbpb(Bmu2pt[b], Bmu2eta[b], -1)*tnp_weight_trk_pbpb(0);
					double sysLoTrg = tnp_weight_trg_pbpb(Bmu1pt[b], Bmu1eta[b], -2)*tnp_weight_trk_pbpb(0)*tnp_weight_trg_pbpb(Bmu2pt[b], Bmu2eta[b], -2)*tnp_weight_trk_pbpb(0);
					double sysHiTrk = tnp_weight_trg_pbpb(Bmu1pt[b], Bmu1eta[b], 0)*tnp_weight_trk_pbpb(-1)*tnp_weight_trg_pbpb(Bmu2pt[b], Bmu2eta[b], 0)*tnp_weight_trk_pbpb(-1);
					double sysLoTrk = tnp_weight_trg_pbpb(Bmu1pt[b], Bmu1eta[b], 0)*tnp_weight_trk_pbpb(-2)*tnp_weight_trg_pbpb(Bmu2pt[b], Bmu2eta[b], 0)*tnp_weight_trk_pbpb(-2);
                    double sysHiMuid = scale*tnp_weight_muid_pbpb(Bmu1pt[b], Bmu1eta[b], -1)*tnp_weight_muid_pbpb(Bmu2pt[b], Bmu2eta[b], -1);
                    double sysLoMuid = scale*tnp_weight_muid_pbpb(Bmu1pt[b], Bmu1eta[b], -2)*tnp_weight_muid_pbpb(Bmu2pt[b], Bmu2eta[b], -2);
                    //double sysHiSTA = scale*tnp_weight_sta_pbpb(Bmu1pt[b], Bmu1eta[b], -1)*tnp_weight_sta_pbpb(Bmu2pt[b], Bmu2eta[b], -1);
                    //double sysLoSTA = scale*tnp_weight_sta_pbpb(Bmu1pt[b], Bmu1eta[b], -2)*tnp_weight_sta_pbpb(Bmu2pt[b], Bmu2eta[b], -2);
					double sysbinnedTrg = tnp_weight_trg_pbpb(Bmu1pt[b], Bmu1eta[b], -10)*tnp_weight_trk_pbpb(0)*tnp_weight_trg_pbpb(Bmu2pt[b], Bmu2eta[b], -10)*tnp_weight_trk_pbpb(0);
					hSysHiTrg ->Fill(Bpt[b], _weight*sysHiTrg);
					hSysHiTrk ->Fill(Bpt[b], _weight*sysHiTrk);
					hSysHiMuid->Fill(Bpt[b], _weight*sysHiMuid);
					//hSysHiSTA ->Fill(Bpt[b], _weight*sysHiSTA);
					hSysLoTrg ->Fill(Bpt[b], _weight*sysLoTrg);
					hSysLoTrk ->Fill(Bpt[b], _weight*sysLoTrk);
					hSysLoMuid->Fill(Bpt[b], _weight*sysLoMuid);
					//hSysLoSTA ->Fill(Bpt[b], _weight*sysLoSTA);
					hSysBinnedTrg ->Fill(Bpt[b], _weight*sysbinnedTrg);
				}
			}
        }//Bsize loop
	}//event loop

	cout<<"Stat. unc."<<endl;
	cout<<"vStatTrg"<<endl;
	for(int j = 0; j < 100; j++){
		for(int i = 0; i < _nBins; i++){	printf("%.2f, ", hStatTrg[j]->GetBinContent(i+1));}printf("\n");
	}
	cout<<"vStatMuid"<<endl;
	for(int j = 0; j < 100; j++){
		for(int i = 0; i < _nBins; i++){	printf("%.2f, ", hStatMuid[j]->GetBinContent(i+1));}printf("\n");
	}
	//cout<<"vStatSTA"<<endl;
	//for(int j = 0; j < 100; j++){
	//	for(int i = 0; i < _nBins; i++){	printf("%.2f, ", hStatSTA[j]->GetBinContent(i+1));}printf("\n");
	//}

	//Calculate stat. unc.
	printf("vStatTrg, vStatMuid, StatAll\n");
    for(int i = 0; i < _nBins; i++){
		double vStatTrg = 0;
		double vStatMuid = 0;
		//double vStatSTA = 0;
		for(int j = 0; j < 100; j++){
			double _stattrg = fabs(hStatTrg[j]->GetBinContent(i+1)-hScale->GetBinContent(i+1))/hScale->GetBinContent(i+1);
			vStatTrg += _stattrg*_stattrg;
//			printf("%.4f, ", _stattrg);
		}
//		cout<<endl;
		for(int j = 0; j < 100; j++){
			double _statmuid = fabs(hStatMuid[j]->GetBinContent(i+1)-hScale->GetBinContent(i+1))/hScale->GetBinContent(i+1);
			vStatMuid += _statmuid*_statmuid;
//			printf("%.4f, ", _statmuid);
		}
//		cout<<endl;
		//for(int j = 0; j < 100; j++){
		//double _statsta = fabs(hStatSTA[j]->GetBinContent(i+1)-hScale->GetBinContent(i+1))/hScale->GetBinContent(i+1);
		//	vStatSTA += _statsta*_statsta;
//			printf("%.4f, ", _statsta);
		//}
//		cout<<endl;
		vStatTrg = sqrt(vStatTrg/100.);
		vStatMuid = sqrt(vStatMuid/100.);
		//vStatSTA = sqrt(vStatSTA/100.);
		double StatAll = vStatTrg*vStatTrg + vStatMuid*vStatMuid;
		printf("%.4f, %.4f, %.4f\n", vStatTrg, vStatMuid, StatAll);
		hStatUnc->SetBinContent(i+1, StatAll);
	}
	
	cout<<"Sys. unc."<<endl;
	cout<<"vNominal"<<endl;
	for(int i = 0; i < _nBins; i++){	printf("%.2f, ", hNominal->GetBinContent(i+1));}printf("\n");
	cout<<"vScale"<<endl;
	for(int i = 0; i < _nBins; i++){	printf("%.2f, ", hScale->GetBinContent(i+1));}printf("\n");
	cout<<"vSysHiTrg"<<endl;
	for(int i = 0; i < _nBins; i++){	printf("%.2f, ", hSysHiTrg->GetBinContent(i+1));}printf("\n");
	cout<<"vSysLoTrg"<<endl;
	for(int i = 0; i < _nBins; i++){	printf("%.2f, ", hSysLoTrg->GetBinContent(i+1));}printf("\n");
	cout<<"vSysHiTrk"<<endl;
	for(int i = 0; i < _nBins; i++){	printf("%.2f, ", hSysHiTrk->GetBinContent(i+1));}printf("\n");
	cout<<"vSysLoTrk"<<endl;
	for(int i = 0; i < _nBins; i++){	printf("%.2f, ", hSysLoTrk->GetBinContent(i+1));}printf("\n");
	cout<<"vSysHiMuid"<<endl;
	for(int i = 0; i < _nBins; i++){	printf("%.2f, ", hSysHiMuid->GetBinContent(i+1));}printf("\n");
	cout<<"vSysLoMuid"<<endl;
	for(int i = 0; i < _nBins; i++){	printf("%.2f, ", hSysLoMuid->GetBinContent(i+1));}printf("\n");
//cout<<"vSysHiSTA"<<endl;
//for(int i = 0; i < _nBins; i++){	printf("%.2f, ", hSysHiSTA->GetBinContent(i+1));}printf("\n");
//cout<<"vSysLoSTA"<<endl;
//for(int i = 0; i < _nBins; i++){	printf("%.2f, ", hSysLoSTA->GetBinContent(i+1));}printf("\n");
	cout<<"vSysBinnedTrg"<<endl;
	for(int i = 0; i < _nBins; i++){	printf("%.2f, ", hSysBinnedTrg->GetBinContent(i+1));}printf("\n");
	
	//Calculate sys. unc.
	printf("vSysTrg, vSysTrk, vSysMuid, vSysBinnedTrg, Sysall\n");
	for(int i = 0; i < _nBins; i++){
		double vSysTrg = max(fabs(hSysHiTrg->GetBinContent(i+1)-hScale->GetBinContent(i+1)), fabs(hSysLoTrg->GetBinContent(i+1)-hScale->GetBinContent(i+1)))/hScale->GetBinContent(i+1);
		double vSysTrk = max(fabs(hSysHiTrk->GetBinContent(i+1)-hScale->GetBinContent(i+1)), fabs(hSysLoTrk->GetBinContent(i+1)-hScale->GetBinContent(i+1)))/hScale->GetBinContent(i+1);
		double vSysMuid = max(fabs(hSysHiMuid->GetBinContent(i+1)-hScale->GetBinContent(i+1)), fabs(hSysLoMuid->GetBinContent(i+1)-hScale->GetBinContent(i+1)))/hScale->GetBinContent(i+1);
		//double vSysSTA = max(fabs(hSysHiSTA->GetBinContent(i+1)-hScale->GetBinContent(i+1)), fabs(hSysLoSTA->GetBinContent(i+1)-hScale->GetBinContent(i+1)))/hScale->GetBinContent(i+1);
		double vSysBinnedTrg = fabs(hSysBinnedTrg->GetBinContent(i+1)-hScale->GetBinContent(i+1))/hScale->GetBinContent(i+1);
		double Sysall = vSysTrg*vSysTrg + vSysTrk*vSysTrk + vSysMuid*vSysMuid + vSysBinnedTrg*vSysBinnedTrg;
		Sysall = sqrt(Sysall);
		printf("%.4f, %.4f, %.4f, %.4f, %.4f\n", vSysTrg, vSysTrk, vSysMuid, vSysBinnedTrg, Sysall);
		hSysUnc->SetBinContent(i+1, Sysall);
	}
	printf("sys + stat uncertainty:\n");
	printf("double tnpUnc_%s[%d] = {", label.c_str(), _nBins);
	for(int i = 0; i < _nBins; i++){
		double sys = hSysUnc->GetBinContent(i+1);
		double stat = hStatUnc->GetBinContent(i+1);
		double all = sqrt(sys*sys + stat*stat);
        hAllUnc->SetBinContent(i+1, all);	
		printf("%4f, ", all*100.);
	}printf("};\n");
	printf("scaleing:\n");
	for(int i = 0; i < _nBins; i++){	
		double SF = hScale->GetBinContent(i+1)/hNominal->GetBinContent(i+1);
	    hSF->SetBinContent(i+1, SF);
		printf("%.4f, ", SF);
	}printf("\n");
	printf("double sf_%s[%d] = {", label.c_str(), _nBins);
	for(int i = 0; i < _nBins; i++){	
		printf("%.2f/%.2f, ", hScale->GetBinContent(i+1), hNominal->GetBinContent(i+1));
	}printf("};\n");
}

int main()
{
  converter();
  return 0;
}
