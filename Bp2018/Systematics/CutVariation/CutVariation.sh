
doppCutVariation=0
doPbPbCutVariation=1

ppData=""
ppMC=""

PbPbData="/mnt/T2_US_MIT/submit-hi2/scratch/gwangjun/crab_Bfinder_20190513_HIDoubleMuonPsiPeri_HIRun2018A_04Apr2019_v1_1033p1_GoldenJSON_327123_327564_skimhltBsize_ntKp.root"
PbPbMC="/mnt/T2_US_MIT/submit-hi2/scratch/gwangjun/crab_Bfinder_20190520_Hydjet_Pythia8_BuToJpsiK_1033p1_pt3tkpt0p7dls2_v2_addSamplePthat_pthatweight.root"

BASECUTPbPb="Btrk1Pt>0.9 && Bpt>5.0 && (BsvpvDistance/BsvpvDisErr)>2.0 && Bchi2cl>0.05 && TMath::Abs(Btrk1Eta)<2.4 && TMath::Abs(By)<2.4 && TMath::Abs(PVz)<15 && Bmass>5 && Bmass<6 && TMath::Abs(Bmumumass-3.096900)<0.15 && Bmu1SoftMuID && Bmu2SoftMuID && Bmu1isAcc && Bmu2isAcc && Bmu1isTriggered && Bmu2isTriggered && (Btrk1PixelHit+Btrk1StripHit)>=11 && (Btrk1Chi2ndf/(Btrk1nStripLayer+Btrk1nPixelLayer))<0.18 && TMath::Abs(Btrk1PtErr/Btrk1Pt)<0.1"
CUTPbPb=${BASECUTPbPb}

if [ $doppCutVariation -eq 1 ]; then
g++ CutVariation.C $(root-config --cflags --libs) -g -o CutVariation.exe 
./CutVariation.exe "pp" "$ppData"  "$ppMC" "$CUTPP"
rm CutVariation.exe
fi

if [ $doPbPbCutVariation -eq 1 ]; then
g++ CutVariation.C $(root-config --cflags --libs) -g -o CutVariation.exe 
./CutVariation.exe "PbPb" "$PbPbData"  "$PbPbMC" "$CUTPbPb"
rm CutVariation.exe
fi
