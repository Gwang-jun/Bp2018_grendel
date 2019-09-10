#!/bin/bash
#source clean.sh

DOTRIGGERSTUDY=1
ISMCPbPb=1

CENTPbPbMIN=0
CENTPbPbMAX=100

# 2018 PbPb
INPUTMCPbPb="/raid5/data/gwangjun/crab_Bfinder_20190624_Hydjet_Pythia8_Official_BuToJpsiK_1033p1_pt3tkpt0p7dls2_allpthat_pthatweight_BDT.root"
INPUTDATAPbPb="/raid5/data/gwangjun/crab_Bfinder_20190513_HIDoubleMuon__PsiPeri__HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_skimhltBsize_ntKp_BDT.root"

LABELPbPb="PbPb"

#LUMIPbPb=56.564165324 #2018 PbPb 0-100% (SUBJECT TO CHANGE!)
LUMIPbPb=62.546428573 #2018 PbPb 0-90% (SUBJECT TO CHANGE!)

NMBEVT=9984367691.260341 #Number of MB events (SUBJECT TO CHANGE!)
#NMBEVT=2329685794.627413

ISDOWEIGHTPbPb=1
SELGENPbPb="TMath::Abs(Gy)<2.4 && TMath::Abs(GpdgId)==521 && GisSignal==1 && GcollisionId==0"
SELGENPbPbACCPbPb="TMath::Abs(Gy)<2.4 && TMath::Abs(GpdgId)==521 && GisSignal==1 && GcollisionId==0 && Gtk1pt>0.7 && TMath::Abs(Gtk1eta)<2.4 && ((TMath::Abs(Gmu1eta)<1.2 && Gmu1pt>3.5) || (TMath::Abs(Gmu1eta)>1.2 && TMath::Abs(Gmu1eta)<2.1 && Gmu1pt>5.47-1.89*TMath::Abs(Gmu1eta)) || (TMath::Abs(Gmu1eta)>2.1 && TMath::Abs(Gmu2eta)<2.4 && Gmu1pt>1.5)) && ((TMath::Abs(Gmu2eta)<1.2 && Gmu2pt>3.5) || (TMath::Abs(Gmu2eta)>1.2 && TMath::Abs(Gmu2eta)<2.1 && Gmu2pt>5.47-1.89*TMath::Abs(Gmu2eta)) || (TMath::Abs(Gmu2eta)>2.1 && TMath::Abs(Gmu2eta)<2.4 && Gmu2pt>1.5))"
RECOONLYPbPb="pprimaryVertexFilter && phfCoincFilter2Th4 && pclusterCompatibilityFilter && Btrk1Pt>0.9 && Bpt>5.0 && (BsvpvDistance/BsvpvDisErr)>2.0 && Bchi2cl>0.05 && TMath::Abs(Btrk1Eta)<2.4 && TMath::Abs(By)<2.4 && TMath::Abs(PVz)<15 && Bmass>5 && Bmass<6 && TMath::Abs(Bmumumass-3.096900)<0.15 && Bmu1SoftMuID && Bmu2SoftMuID && ((TMath::Abs(Bmu1eta)<1.2 && Bmu1pt>3.5) || (TMath::Abs(Bmu1eta)>1.2 && TMath::Abs(Bmu1eta)<2.1 && Bmu1pt>5.47-1.89*TMath::Abs(Bmu1eta)) || (TMath::Abs(Bmu1eta)>2.1 && TMath::Abs(Bmu1eta)<2.4 && Bmu1pt>1.5)) && ((TMath::Abs(Bmu2eta)<1.2 && Bmu2pt>3.5) || (TMath::Abs(Bmu2eta)>1.2 && TMath::Abs(Bmu2eta)<2.1 && Bmu2pt>5.47-1.89*TMath::Abs(Bmu2eta)) || (TMath::Abs(Bmu2eta)>2.1 && TMath::Abs(Bmu2eta)<2.4 && Bmu2pt>1.5)) && Bmu1isTriggered && Bmu2isTriggered && (Btrk1PixelHit+Btrk1StripHit)>=11 && (Btrk1Chi2ndf/(Btrk1nStripLayer+Btrk1nPixelLayer))<0.18 && TMath::Abs(Btrk1PtErr/Btrk1Pt)<0.1"

#GA
BASECUTPbPb="pprimaryVertexFilter && phfCoincFilter2Th4 && pclusterCompatibilityFilter && Btrk1Pt>0.9 && Bpt>5.0 && (BsvpvDistance/BsvpvDisErr)>2.0 && Bchi2cl>0.05 && TMath::Abs(Btrk1Eta)<2.4 && TMath::Abs(By)<2.4 && TMath::Abs(PVz)<15 && Bmass>5 && Bmass<6 && TMath::Abs(Bmumumass-3.096900)<0.15 && Bmu1SoftMuID && Bmu2SoftMuID && ((TMath::Abs(Bmu1eta)<1.2 && Bmu1pt>3.5) || (TMath::Abs(Bmu1eta)>1.2 && TMath::Abs(Bmu1eta)<2.1 && Bmu1pt>5.47-1.89*TMath::Abs(Bmu1eta)) || (TMath::Abs(Bmu1eta)>2.1 && TMath::Abs(Bmu1eta)<2.4 && Bmu1pt>1.5)) && ((TMath::Abs(Bmu2eta)<1.2 && Bmu2pt>3.5) || (TMath::Abs(Bmu2eta)>1.2 && TMath::Abs(Bmu2eta)<2.1 && Bmu2pt>5.47-1.89*TMath::Abs(Bmu2eta)) || (TMath::Abs(Bmu2eta)>2.1 && TMath::Abs(Bmu2eta)<2.4 && Bmu2pt>1.5)) && Bmu1isTriggered && Bmu2isTriggered && (Btrk1PixelHit+Btrk1StripHit)>=11 && (Btrk1Chi2ndf/(Btrk1nStripLayer+Btrk1nPixelLayer))<0.18 && TMath::Abs(Btrk1PtErr/Btrk1Pt)<0.1" ## pre-filter
CUTPbPb="$BASECUTPbPb && ((Bpt>5 && Bpt<7 && BDT_5_7>-0.01) || (Bpt>7 && Bpt<10 && BDT_7_10>0.08) || (Bpt>10 && Bpt<15 && BDT_10_15>0.09) || (Bpt>15 && Bpt<20 && BDT_15_20>0.09) || (Bpt>20 && Bpt<30 && BDT_20_30>0.07) || (Bpt>30 && Bpt<50 && BDT_30_50>0.12) || (Bpt>50 && Bpt<100 && BDT_50_100>0.24))"

TRGPbPb="HLT_HIL1DoubleMuOpen_Centrality_50_100_v1"
TRGPbPbMC="HLT_HIL1DoubleMuOpen_Centrality_50_100_v1"
OUTPUTFILEPbPb="ROOTfiles/Trigger.root"
NPFIT_PbPb="701.019629*TMath::Erf((x-5.140349)/-0.035471)+701.019629+16.946432*TMath::Gaus(x,5.343914,0.040000)/(sqrt(2*3.14159)*0.040000)"

if [ $DOTRIGGERSTUDY -eq 1 ]; then      
g++ Trigger.C $(root-config --cflags --libs) -g -o Trigger.exe 
./Trigger.exe 1 "$INPUTDATAPbPb"  "$INPUTMCPbPb"  "$TRGPbPb" "$CUTPbPb"   "$SELGENPbPb"   "$ISMCPbPb"   1   "$ISDOWEIGHTPbPb"   "$LABELPbPb"  "$OUTPUTFILEPbPb" "$NPFIT_PbPb" 0 "$CENTPbPbMIN" "$CENTPbPbMAX"
rm Trigger.exe
fi 
