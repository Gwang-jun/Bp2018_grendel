#include "uti.h"
#include "parameters.h"
#include "TF1.h"

Double_t setparam0=100.;
Double_t setparam1=5.28;
Double_t setparam2=0.05;
Double_t setparam3=0.03;
Double_t fixparam1=5.279;

Double_t minhisto=5.0;
Double_t maxhisto=6.0;
Double_t nbinsmasshisto=50;
Double_t binwidthmass=(maxhisto-minhisto)/nbinsmasshisto;

Double_t hiBinMin = 0.0; Double_t hiBinMax = 200.0;

//TString trgselection="(HLT_HIL1DoubleMu0_v1 || HLT_HIL1DoubleMu0_part1_v1 || HLT_HIL1DoubleMu0_part2_v1 || HLT_HIL1DoubleMu0_part3_v1)";
TString cut="Btrk1Pt>1.0 && Bpt>5.0 && (BsvpvDistance/BsvpvDisErr)>3.5 && Bchi2cl>0.05 && TMath::Abs(Btrk1Eta)<2.4 && TMath::Abs(By)<2.4";
//TString cut="Btrk1Pt>1.0 && Bpt>5.0 && (BsvpvDistance/BsvpvDisErr)>3.5 && Bchi2cl>0.05 && TMath::Abs(Btrk1Eta)<2.4 && TMath::Abs(By)<2.4 && TMath::Abs(PVz)<15 && Bmass>5 && Bmass<6 && Bmu1isGlobalMuon && Bmu2isGlobalMuon && Bmu1dxyPV<0.3 && Bmu2dxyPV<0.3 && Bmu1dzPV<20 && Bmu2dzPV<20 && Bmu1InPixelLayer>0 && (Bmu1InPixelLayer+Bmu1InStripLayer)>5 && Bmu2InPixelLayer>0 && (Bmu2InPixelLayer+Bmu2InStripLayer)>5 && ((TMath::Abs(Bmu1eta)<1.2 && Bmu1pt>3.5) || (TMath::Abs(Bmu1eta)>1.2 && TMath::Abs(Bmu1eta)<2.1 && Bmu1pt>(5.77-1.8*TMath::Abs(Bmu1eta))) || (TMath::Abs(Bmu1eta)>2.1 && TMath::Abs(Bmu1eta)<2.4 && Bmu1pt>1.8)) && TMath::Abs(Bmumumass-3.096900)<0.15 && Bmu1TMOneStationTight && Bmu2TMOneStationTight && Btrk1highPurity && (Btrk1PixelHit+Btrk1StripHit)>=11 && (Btrk1Chi2ndf/(Btrk1nStripLayer+Btrk1nPixelLayer))<0.18 && TMath::Abs(Btrk1PtErr/Btrk1Pt)<0.1";
TString Nonprompt_fit="3.39711e1*TMath::Gaus(x,5.35002e0,3.87952e-2)*(x<5.35002e0)+(x>=5.35002e0)*3.39711e1*TMath::Gaus(x,5.35002e0,1.14232e-1)+2.16376e2*TMath::Erf(-(x-5.14189)/8.66243e-2) + 2.16376e2";

// 2018 PbPb data(bkg)
TString inputbkg="/mnt/hadoop/cms/store/user/jwang/BntupleRun2018/crab_Bfinder_20181220_HIDoubleMuon_HIRun2018A_PromptReco_v1v2_1031_NoJSON_skimhltBsize_ntKp.root";
TString outputfile="bkgprime.root";

TF1 *fit(TTree *nt, Double_t ptmin, Double_t ptmax)
{
  static Int_t count=0;
  count++;
  TCanvas* c= new TCanvas(Form("c%d",count),"",600,600);
  TH1D* h = new TH1D(Form("h-%d",count),"",nbinsmasshisto,minhisto,maxhisto);

  TString seldata = Form("%s",cut.Data());
  TString iNP = Form("%s",Nonprompt_fit.Data());

  TF1* f = new TF1(Form("f%d",count),"[0]*([7]*TMath::Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*TMath::Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[5]*("+iNP+")");
  f->SetNpx(5000);
  f->SetLineWidth(5);
  nt->Project(Form("h-%d",count),"Bmass",Form("(%s&&Bpt>%f&&Bpt<%f)",seldata.Data(),ptmin,ptmax));
  // clean0(h);

  //h->Draw();
  f->SetParLimits(4,-1000,1000);
  f->SetParLimits(2,0.01,1000);
  f->SetParLimits(8,0.01,1000);
  f->SetParLimits(7,0,1);
  f->SetParLimits(5,0,1000);

  f->SetParameter(0,setparam0);
  f->SetParameter(1,setparam1);
  f->SetParameter(2,setparam2);
  f->SetParameter(8,setparam3);
  f->FixParameter(1,fixparam1);
  f->FixParameter(5,0);
  //h->GetEntries();

  f->FixParameter(1,f->GetParameter(1));
  f->FixParameter(2,f->GetParameter(2));
  f->FixParameter(7,f->GetParameter(7));
  f->FixParameter(8,f->GetParameter(8));
  f->ReleaseParameter(5);
  printf("Fixed para.:\n");
  printf("%f, %f\n", f->GetParameter(3), f->GetParameter(4));
  h->Fit(Form("f%d",count),"q","",minhisto,maxhisto);
  h->Fit(Form("f%d",count),"q","",minhisto,maxhisto);
  f->ReleaseParameter(1);
  h->Fit(Form("f%d",count),"L q","",minhisto,maxhisto);
  h->Fit(Form("f%d",count),"L q","",minhisto,maxhisto);
  h->Fit(Form("f%d",count),"L q","",minhisto,maxhisto);
  h->Fit(Form("f%d",count),"L m","",minhisto,maxhisto);

  TF1 *background = new TF1(Form("background%d",count),"[0]+[1]*x");
  background->SetParameter(0,f->GetParameter(3));
  background->SetParameter(1,f->GetParameter(4));
  background->SetLineColor(4);
  background->SetRange(minhisto,maxhisto);
  //background->SetLineStyle(2);//PAS
  background->SetLineStyle(7);//paper
  //background->SetLineWidth(5);
  background->SetLineWidth(9);

  TF1 *Bkpi = new TF1(Form("fBkpi%d",count),"[0]*("+iNP+")");
  Bkpi->SetParameter(0,f->GetParameter(5));
  Bkpi->SetRange(minhisto,maxhisto);
  Bkpi->SetLineStyle(1);
  //Bkpi->SetFillStyle(3004);//PAS
  Bkpi->SetFillStyle(3005);//paper
  //Bkpi->SetLineColor(kGreen+1);//PAS
  //Bkpi->SetFillColor(kGreen+1);//PAS
  Bkpi->SetLineColor(kGreen+4);//paper
  Bkpi->SetFillColor(kGreen+4);//paper
  //Bkpi->SetLineWidth(5);
  Bkpi->SetLineWidth(9);
  
  Double_t Bkg_comb = 0.8*(background->Integral(4.97932,5.07932)+background->Integral(5.47932,5.57932))/binwidthmass; //B+ mass_pdg=5.27932GeV, bkg region = pm 0.2~0.3GeV
  Double_t Bkg_np = 0.8*(Bkpi->Integral(4.97932,5.07932)+Bkpi->Integral(5.47932,5.57932))/binwidthmass;
  printf("p_t bin %.0f-%.0f     bkgprime(comb): %f     bkgprime(comb+np): %f\n",ptmin,ptmax,Bkg_comb,Bkg_comb+Bkg_np);

  return background;
}

int main()
{
  
  TFile* inf = new TFile(inputbkg.Data());

  //For 2018 PbPb bkg(data)                                                                                                                    
  TTree* nt = (TTree*)inf->Get("Bfinder/ntKp");
  nt->AddFriend("hltanalysis/HltTree");
  nt->AddFriend("hiEvtAnalyzer/HiTree");
  nt->AddFriend("skimanalysis/HltTree");

  TString outputf;
  outputf = Form("%s",outputfile.Data());
  TFile* outf = new TFile(outputf.Data(),"recreate");
  outf->cd();
  TH1D* bkgprime = new TH1D("bkgprime","",nBins,ptBins);
  TCanvas *c = new TCanvas("c", "bkgprime");
  
  double bkgref[nBins];
  double bkgl, bkgh;

for(int i=0;i<nBins;i++)
{
    // if(i==(nBins-1)){nbinsmasshisto=20;binwidthmass=(maxhisto-minhisto)/nbinsmasshisto;}
    
    TF1* f = fit(nt,ptBins[i],ptBins[i+1]);
    //bkgref[i] = bkgl+bkgh;
    //bkgprime->Fill(bkgref[i]);
    //printf("p_t bin %.0f-%.0f    bkgprime: %f     bkg: %f\n", ptBins[i], ptBins[i+1], bkgref[i], 0.8*bkgref[i]);
}
//bkgprime->Write();

 return 1;

}
