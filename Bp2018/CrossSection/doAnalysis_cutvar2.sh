#!/bin/bash
#source clean.sh
CENTPbPbMIN=0
CENTPbPbMAX=90
#Central point of the analysis

DOANALYSISPP_CUTVAR=0
DOANALYSISPbPb_CUTVAR=1

#LABEL_VARNUM=0

# 2015 PP
INPUTMCPP="/raid5/data/BntupleRun2018/MC/Bntuple20160816_Bpt7svpv5p5Bpt10svpv3p5_BfinderMC_pp_Pythia8_BuToJpsiK_TuneCUEP8M1_20160816_bPt0jpsiPt0tkPt0p5_Bp_pthatweight_JingBDT.root"
INPUTDATAPP="/raid5/data/BntupleRun2018/Data/Bntuple20160816_Bpt7svpv5p5Bpt10svpv3p5_BfinderData_pp_20160816_bPt0jpsiPt0tkPt0p5_Bp_JingBDT.root"

# 2018 PbPb
INPUTMCPbPb="/raid5/data/gwangjun/crab_Bfinder_20190624_Hydjet_Pythia8_Official_BuToJpsiK_1033p1_pt3tkpt0p7dls2_allpthat_pthatweight.root"
INPUTDATAPbPb="/raid5/data/gwangjun/crab_Bfinder_20190513_HIDoubleMuon__PsiPeri__HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_skimhltBsize_ntKp.root"

LABELPP="pp"
LUMIPP=27.748

ISMCPP=0
ISDOWEIGHTPP=0
SELGENPP="TMath::Abs(Gy)<2.4 && abs(GpdgId)==521 && GisSignal==1 && GcollisionId==0"
SELGENPPACCPP="TMath::Abs(Gy)<2.4 && abs(GpdgId)==521 && GisSignal==1 && GcollisionId==0 && TMath::Abs(Gy)<2.4 && ((TMath::Abs(Gmu1eta)<1.2 && Gmu1pt>=3.5) || (TMath::Abs(Gmu1eta)>=1.2 && TMath::Abs(Gmu1eta)<2.1 && Gmu1pt>=(5.77-1.89*TMath::Abs(Gmu1eta))) || (TMath::Abs(Gmu1eta)>=2.1 && Gmu1pt>=1.8)) && ((TMath::Abs(Gmu2eta)<1.2 && Gmu2pt>3.5) || (TMath::Abs(Gmu2eta)>=1.2 && TMath::Abs(Gmu2eta)<2.1 && Gmu2pt>=(5.77-1.89*TMath::Abs(Gmu2eta))) || (TMath::Abs(Gmu2eta)>=2.1 && Gmu2pt>=1.8)) && Gtk1pt>0.9 && TMath::Abs(Gtk1eta)<2.4"
#pp GA
BASECUTPP="abs(PVz)<15&&pBeamScrapingFilter&&pPAprimaryVertexFilter&&TMath::Abs(By)<2.4&&TMath::Abs(Bmumumass-3.096900)<0.15&&Bmass>5&&Bmass<6&& ((abs(Bmu1eta)<1.2 && Bmu1pt>3.5) || (abs(Bmu1eta)>1.2 && abs(Bmu1eta)<2.1 && Bmu1pt>(5.77-1.8*abs(Bmu1eta))) || (abs(Bmu1eta)>2.1 && abs(Bmu1eta)<2.4 && Bmu1pt>1.8)) && ((abs(Bmu2eta)<1.2 && Bmu2pt>3.5) || (abs(Bmu2eta)>1.2 && abs(Bmu2eta)<2.1 && Bmu2pt>(5.77-1.8*abs(Bmu2eta))) || (abs(Bmu2eta)>2.1 && abs(Bmu2eta)<2.4 && Bmu2pt>1.8)) && Bmu1TMOneStationTight && Bmu2TMOneStationTight && Bmu1InPixelLayer > 0 && (Bmu1InPixelLayer+Bmu1InStripLayer) > 5 && Bmu2InPixelLayer > 0 && (Bmu2InPixelLayer+Bmu2InStripLayer) > 5 && Bmu1dxyPV< 0.3 && Bmu2dxyPV< 0.3 && Bmu1dzPV<20 && Bmu2dzPV<20 && Bmu1isGlobalMuon && Bmu2isGlobalMuon && Bmu1TrgMatchFilterE>0 && Bmu2TrgMatchFilterE>0 && Btrk1highPurity && TMath::Abs(Btrk1PtErr/Btrk1Pt)<0.3 && abs(Btrk1Eta)<2.4 && Btrk1Pt>0.5 && Bchi2cl>0.005 && (BsvpvDistance/BsvpvDisErr)>3.5"
CUTPP="$BASECUTPP && ((Bpt>5 && Bpt<10 && Btrk1Pt>0.748 && Bchi2cl>0.0322 && (BsvpvDistance/BsvpvDisErr)>0 && cos(Bdtheta)>0.989 && abs(Btrk1Eta)<2.4)|| (Bpt>10 && Bpt<15 && Btrk1Pt>0.878 && Bchi2cl>0.0014 && (BsvpvDistance/BsvpvDisErr)>0 && cos(Bdtheta)>0.971 && abs(Btrk1Eta)<2.39) || (Bpt>15 && Bpt<20 && Btrk1Pt>0.818 && Bchi2cl>0.00838 && (BsvpvDistance/BsvpvDisErr)>0 && cos(Bdtheta)>0.959 && abs(Btrk1Eta)<2.42) || (Bpt>20 && Bpt<30 && Btrk1Pt>0.84 && Bchi2cl>0.014 && (BsvpvDistance/BsvpvDisErr)>0 && cos(Bdtheta)>0.603 && abs(Btrk1Eta)<2.38) || (Bpt>30 && Bpt<100 && Btrk1Pt>1.06 && Bchi2cl>0.0155 && (BsvpvDistance/BsvpvDisErr)>3.53 && cos(Bdtheta)>0.995 && abs(Btrk1Eta)<2.37))"

RECOONLYPP=$CUTPP
TRGPP="(HLT_HIL1DoubleMu0_v1)"
TRGPPMC="(HLT_HIL1DoubleMu0ForPPRef_v1)"
#TRGPPMC="(1)"
OUTPUTFILEPPSAVEHIST="ROOTfiles/hPtSpectrumSaveHistBplusPP.root"
OUTPUTFILEPP="ROOTfiles/hPtSpectrumBplusPP.root"
PREFIXPP="ROOTfiles/hPtSpectrumBplusPP" #IF YOU CHANGE OUTPUTFILEPP, CHANGE PREFIX PP ACCORDINGLY
OUTPUTFILEMCSTUDYPP="ROOTfiles/MCstudiesPP.root"
OUTPUTFILEPlotPP="ROOTfiles/CrossSectionPP.root"
OUTPUTFILEPPDATA="ROOTfiles/data_pp.root"
OUTPUTFILEPPMC="ROOTfiles/mc_pp.root"

NPFIT_PP="3.12764e1*TMath::Gaus(x,5.33166,3.64663e-2)*(x<5.33166)+(x>=5.33166)*3.12764e1*TMath::Gaus(x,5.33166,1.5204e-1)+2.11124e2*TMath::Erf(-(x-5.14397)/6.43194e-2) + 2.11124e2"

if [ $DOANALYSISPP_CUTVAR -eq 1 ]; then      
g++ fitB_cutvar.C $(root-config --cflags --libs) -g -o fitB_cutvar.exe 
./fitB_cutvar.exe 0 "$INPUTDATAPP"  "$INPUTMCPP"  "$TRGPP" "$CUTPP"   "$SELGENPP"   "$ISMCPP"   1   "$ISDOWEIGHTPP"   "$LABELPP"  "$OUTPUTFILEPP" "$NPFIT_PP" 0 "$LABEL_VARNUM" 
rm fitB_cutvar.exe
fi 

LABELPbPb="PbPb"
LUMIPbPb=62.546428573 #2018 PbPb 0-90% (SUBJECT TO CHANGE!)

NMBEVT=9984367691.260341 #Number of MB events (SUBJECT TO CHANGE!)
ISMCPbPb=0
ISDOWEIGHTPbPb=1
SELGENPbPb="TMath::Abs(Gy)<2.4 && abs(GpdgId)==521 && GisSignal==1 && GcollisionId==0"
SELGENPbPbACCPbPb="TMath::Abs(Gy)<2.4 && abs(GpdgId)==521 && GisSignal==1 && GcollisionId==0 && TMath::Abs(Gy)<2.4 && ((TMath::Abs(Gmu1eta)<1.2 && Gmu1pt>=3.5) || (TMath::Abs(Gmu1eta)>=1.2 && TMath::Abs(Gmu1eta)<2.1 && Gmu1pt>=(5.77-1.89*TMath::Abs(Gmu1eta))) || (TMath::Abs(Gmu1eta)>=2.1 && Gmu1pt>=1.8)) && ((TMath::Abs(Gmu2eta)<1.2 && Gmu2pt>=3.5) || (TMath::Abs(Gmu2eta)>=1.2 && TMath::Abs(Gmu2eta)<2.1 && Gmu2pt>=(5.77-1.89*TMath::Abs(Gmu2eta))) || (TMath::Abs(Gmu2eta)>=2.1 && Gmu2pt>=1.8)) && Gtk1pt>0.9 && TMath::Abs(Gtk1eta)<2.4"
RECOONLYPbPb="Btrk1Pt>0.7 && Bpt>3.0 && (BsvpvDistance/BsvpvDisErr)>2.0 && Bchi2cl>0.05"

#GA
BASECUTPbPb="Btrk1Pt>0.9 && Bpt>5.0 && (BsvpvDistance/BsvpvDisErr)>2.0 && Bchi2cl>0.05 && TMath::Abs(Btrk1Eta)<2.4 && TMath::Abs(By)<2.4 && TMath::Abs(PVz)<15 && Bmass>5 && Bmass<6 && TMath::Abs(Bmumumass-3.096900)<0.15 && Bmu1SoftMuID && Bmu2SoftMuID && (((TMath::Abs(Bmu1eta[j])<1.2 && Bmu1pt[j])>3.5) || (TMath::Abs(Bmu1eta[j])>1.2 && TMath::Abs(Bmu1eta[j])<2.1 && Bmu1pt[j]>5.47-1.89*TMath::Abs(Bmu1eta[j])) || (TMath::Abs(Bmu1eta[j])>2.1 && TMath::Abs(Bmu1eta[j])<2.4 && Bmu1pt[j]>1.5)) && (((TMath::Abs(Bmu2eta[j])<1.2 && Bmu2pt[j])>3.5) || (TMath::Abs(Bmu2eta[j])>1.2 && TMath::Abs(Bmu2eta[j])<2.1 && Bmu2pt[j]>5.47-1.89*TMath::Abs(Bmu2eta[j])) || (TMath::Abs(Bmu2eta[j])>2.1 && TMath::Abs(Bmu2eta[j])<2.4 && Bmu2pt[j]>1.5)) && Bmu1isTriggered && Bmu2isTriggered && (Btrk1PixelHit+Btrk1StripHit)>=11 && (Btrk1Chi2ndf/(Btrk1nStripLayer+Btrk1nPixelLayer))<0.18 && TMath::Abs(Btrk1PtErr/Btrk1Pt)<0.1" ## pre-filter
CUTPbPb="$BASECUTPbPb"

TRGPbPb="(HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1)"
TRGPbPbMC="(HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1)"
OUTPUTFILEPbPbSAVEHIST="ROOTfiles/hPtSpectrumSaveHistBplusPbPb.root"
OUTPUTFILEPbPb="ROOTfiles/hPtSpectrumBplusPbPb.root"
PREFIXPbPb="ROOTfiles/hPtSpectrumBplusPbPb" #IF YOU CHANGE OUTPUTFILEPbPb, CHANGE PREFIXPbPb ACCORDINGLY
OUTPUTFILEMCSTUDYPbPb="ROOTfiles/MCstudiesPbPb.root"
OUTPUTFILEPlotPbPb="ROOTfiles/CrossSectionPbPb.root"
OUTPUTFILEPbPbDATA="ROOTfiles/data_PbPb.root"
OUTPUTFILEPbPbMC="ROOTfiles/mc_PbPb.root"

#NPFIT_PbPb="787.503468*TMath::Erf(-(x-5.105477)/0.080109)+787.503468+61.466850*(0.365515*TMath::Gaus(x,5.340806,0.056890)/(sqrt(2*3.14159)*0.056890)+(1-0.365515)*TMath::Gaus(x,5.340806,0.279220)/(sqrt(2*3.14159)*0.279220))+43.437777*TMath::Gaus(x,5.105719,0.025197)/(sqrt(2*3.14159)*0.025197)"
NPFIT_PbPb="701.019629*TMath::Erf((x-5.140349)/-0.035471)+701.019629+16.946432*TMath::Gaus(x,5.343914,0.040000)/(sqrt(2*3.14159)*0.040000)"

if [ $DOANALYSISPbPb_CUTVAR -eq 1 ]; then      
g++ fitB_cutvar2.C $(root-config --cflags --libs) -g -o fitB_cutvar2.exe 
#./fitB_cutvar2.exe 1 "$INPUTDATAPbPb"  "$INPUTMCPbPb"  "$TRGPbPb" "$CUTPbPb"   "$SELGENPbPb"   "$ISMCPbPb"   1   "$ISDOWEIGHTPbPb"   "$LABELPbPb"  "$OUTPUTFILEPbPb" "$NPFIT_PbPb" 0 "$CENTPbPbMIN" "$CENTPbPbMAX" 0
#./fitB_cutvar2.exe 1 "$INPUTDATAPbPb"  "$INPUTMCPbPb"  "$TRGPbPb" "$CUTPbPb"   "$SELGENPbPb"   "$ISMCPbPb"   1   "$ISDOWEIGHTPbPb"   "$LABELPbPb"  "$OUTPUTFILEPbPb" "$NPFIT_PbPb" 0 "$CENTPbPbMIN" "$CENTPbPbMAX" 1
#./fitB_cutvar2.exe 1 "$INPUTDATAPbPb"  "$INPUTMCPbPb"  "$TRGPbPb" "$CUTPbPb"   "$SELGENPbPb"   "$ISMCPbPb"   1   "$ISDOWEIGHTPbPb"   "$LABELPbPb"  "$OUTPUTFILEPbPb" "$NPFIT_PbPb" 0 "$CENTPbPbMIN" "$CENTPbPbMAX" 2
./fitB_cutvar2.exe 1 "$INPUTDATAPbPb"  "$INPUTMCPbPb"  "$TRGPbPb" "$CUTPbPb"   "$SELGENPbPb"   "$ISMCPbPb"   1   "$ISDOWEIGHTPbPb"   "$LABELPbPb"  "$OUTPUTFILEPbPb" "$NPFIT_PbPb" 0 "$CENTPbPbMIN" "$CENTPbPbMAX" 3
#./fitB_cutvar2.exe 1 "$INPUTDATAPbPb"  "$INPUTMCPbPb"  "$TRGPbPb" "$CUTPbPb"   "$SELGENPbPb"   "$ISMCPbPb"   1   "$ISDOWEIGHTPbPb"   "$LABELPbPb"  "$OUTPUTFILEPbPb" "$NPFIT_PbPb" 0 "$CENTPbPbMIN" "$CENTPbPbMAX" 4
rm fitB_cutvar2.exe
fi
