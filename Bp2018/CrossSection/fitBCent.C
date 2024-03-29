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

TString weight;
TString weightgen;
TString weightdata;
TString seldata, seldataref;
TString selmc, selmcref;
TString selmceff, selmceffref;
TString selmcgen, selmcgenref;
TString collisionsystem;
Float_t hiBinMin,hiBinMax,centMin,centMax;
double _ErrCor=1;

int _nBins = nBinsCent;
double *_ptBins = ptBinsCent;

Double_t yield;
Double_t yieldErr;

void fitBCent(int usePbPb=0, TString inputdata="" , TString inputmc="", TString trgselection="1",  TString cut="", TString cutmcgen="", int isMC=0, Double_t luminosity=1., int doweight=0, TString collsyst="PbPb", TString outputfile="", TString npfit="0", int doDataCor = 0, Float_t centmin=0., Float_t centmax=100.)
{
  collisionsystem=collsyst;  
  
  double ErrorOnSigma(double width, double errwidth, double smear, double errsmearing);
  
  if (!(usePbPb==1||usePbPb==0)) std::cout<<"ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!, you are using a non valid isPbPb option"<<std::endl;
  bool isPbPb=(bool)(usePbPb);

  seldataref = Form("%s&&%s",trgselection.Data(),cut.Data());
  selmceffref = Form("%s&&%s",trgselection.Data(),cut.Data());
  selmcgenref = Form("%s",cutmcgen.Data());
  selmcref = Form("%s&&%s",trgselection.Data(),cut.Data());

  gStyle->SetTextSize(0.05);
  gStyle->SetTextFont(42);
  gStyle->SetPadRightMargin(0.043);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.145);
  gStyle->SetTitleX(.0f);

  void clean0 (TH1D* h);
  void getNPFnPar(TString npfname, float par[]);
  //TF1* fit (TTree* nt, TTree* ntMC, double ptmin, double ptmax, int isMC,bool, TF1* &total,Float_t centmin, Float_t centmax, float NPpar[]);
  TF1* fit (TTree* nt, TTree* ntMC, double ptmin, double ptmax, int isMC,bool, TF1* &total,Float_t centMin, Float_t centMax, TString npfit);
  //float NPpar[2];
  //getNPFnPar(npfile, NPpar);
  //std::cout<<"NP parameter 0: "<<NPpar[0]<<std::endl;
  //std::cout<<"NP parameter 1: "<<NPpar[1]<<std::endl;

  weightdata="1";
  if(!isPbPb)
    {
     weightgen="pthatweight*(0.599212+-0.020703*Gpt+0.003143*Gpt*Gpt+-0.000034*Gpt*Gpt*Gpt)";
     weight="pthatweight*(0.599212+-0.020703*Bgenpt+0.003143*Bgenpt*Bgenpt+-0.000034*Bgenpt*Bgenpt*Bgenpt)";
    }
  else
    {
      //weightgen="pthatweight";
      //weight="pthatweight*Ncoll*(TMath::Gaus(PVz,0.427450,4.873825)/(sqrt(2*3.14159)*4.873825))/(TMath::Gaus(PVz,0.909938,4.970989)/(sqrt(2*3.14159)*4.970989))";
      weightgen="pthatweight*((3.506006+0.963473*Gpt+-0.258731*Gpt*Gpt)*TMath::Exp(-0.386065*Gpt)+1.139897)";
      weight="pthatweight*Ncoll*(TMath::Gaus(PVz,0.427450,4.873825)/(sqrt(2*3.14159)*4.873825))/(TMath::Gaus(PVz,0.909938,4.970989)/(sqrt(2*3.14159)*4.970989))*((3.506006+0.963473*Bgenpt+-0.258731*Bgenpt*Bgenpt)*TMath::Exp(-0.386065*Bgenpt)+1.139897)";
      //weightgen="pthatweight*((2.907795+-0.436572*Gpt+0.006372*Gpt*Gpt)*TMath::Exp(-0.157563*Gpt)+1.01308)";
      //weight="pthatweight*Ncoll*(TMath::Gaus(PVz,0.427450,4.873825)/(sqrt(2*3.14159)*4.873825))/(TMath::Gaus(PVz,0.909938,4.970989)/(sqrt(2*3.14159)*4.970989))*((2.907795+-0.436572*Bgenpt+0.006372*Bgenpt*Bgenpt)*TMath::Exp(-0.157563*Bgenpt)+1.01308)";
      //weightgen="pthatweight*(0.889175+0.000791*Gpt+0.000015*Gpt*Gpt)";
      //weight="pthatweight*Ncoll*(TMath::Gaus(PVz,0.427450,4.873825)/(sqrt(2*3.14159)*4.873825))/(TMath::Gaus(PVz,0.909938,4.970989)/(sqrt(2*3.14159)*4.970989))*(0.889175+0.000791*Bgenpt+0.000015*Bgenpt*Bgenpt)";
      //weightgen="pthatweight*Ncoll*(1.034350*TMath::Exp(-0.000844*(PVz+3.502992)*(PVz+3.502992)))*(0.715021+0.039896*Gpt-0.000834*Gpt*Gpt+0.000006*Gpt*Gpt*Gpt)"; // MC Gpt
      //weight="pthatweight*Ncoll*(1.034350*TMath::Exp(-0.000844*(PVz+3.502992)*(PVz+3.502992)))*(0.715021+0.039896*Bgenpt-0.000834*Bgenpt*Bgenpt+0.000006*Bgenpt*Bgenpt*Bgenpt)"; // MC Gpt
    }

  std::cout<<"we are using weight="<<weight<<std::endl;
  std::cout<<"we are using weightdata="<<weightdata<<std::endl;
  std::cout<<"we are using total centrality="<<centmin<<"-"<<centmax<<"%"<<std::endl;
  
  TFile* inf = new TFile(inputdata.Data());
  TFile* infMC = new TFile(inputmc.Data());

  //For 2018 PbPb data, MC
  TTree* nt = (TTree*)inf->Get("Bfinder/ntKp");
  nt->AddFriend("hltanalysis/HltTree");
  nt->AddFriend("hiEvtAnalyzer/HiTree");
  nt->AddFriend("skimanalysis/HltTree");
  nt->AddFriend("BDT");
  TTree* ntGen = (TTree*)infMC->Get("Bfinder/ntGen");
  ntGen->AddFriend("hltanalysis/HltTree");
  ntGen->AddFriend("hiEvtAnalyzer/HiTree");
  ntGen->AddFriend("Bfinder/ntKp"); //call PVz
  ntGen->AddFriend("skimanalysis/HltTree");
  ntGen->AddFriend("BDT");
  TTree* ntMC = (TTree*)infMC->Get("Bfinder/ntKp");
  ntMC->AddFriend("hltanalysis/HltTree");
  ntMC->AddFriend("hiEvtAnalyzer/HiTree");
  ntMC->AddFriend("Bfinder/ntGen"); //call Bgen
  ntMC->AddFriend("skimanalysis/HltTree");
  ntMC->AddFriend("BDT");

  /*
  //For 2015 PbPb, pp data, MC
  TTree* nt = (TTree*)inf->Get("ntKp");
  nt->AddFriend("ntHlt");
  nt->AddFriend("ntHi");
  nt->AddFriend("ntSkim");
  nt->AddFriend("mvaTree");
  TTree* ntGen = (TTree*)infMC->Get("ntGen");
  ntGen->AddFriend("ntHlt");
  ntGen->AddFriend("ntHi");
  TTree* ntMC = (TTree*)infMC->Get("ntKp");
  ntMC->AddFriend("ntHlt");
  ntMC->AddFriend("ntHi");
  ntMC->AddFriend("ntSkim");
  ntMC->AddFriend("mvaTree");
  ntMC->AddFriend("ntGen");
  */

  TH1D* hMean = new TH1D("hMean","",_nBins,_ptBins);                       
  TH1D* hSigmaGaus1 = new TH1D("hSigmaGaus1","",_nBins,_ptBins); 
  TH1D* hSigmaGaus2 = new TH1D("hSigmaGaus2","",_nBins,_ptBins); 
  TF1 *totalmass;

  TH1D* hPt = new TH1D("hPt","",_nBins,_ptBins);
  TH1D* hPtRecoTruth = new TH1D("hPtRecoTruth","",_nBins,_ptBins);
  TH1D* hPtMC = new TH1D("hPtMC","",_nBins,_ptBins);
  TH1D* hPtGen = new TH1D("hPtGen","",_nBins,_ptBins);

  for(int i=0;i<_nBins;i++)
    {
      //TF1* f = fit(nt,ntMC,_ptBins[i],_ptBins[i+1],isMC,isPbPb, totalmass,centmin, centmax, NPpar);
      
      if(isPbPb)
	{
	  hiBinMin = _ptBins[i], hiBinMax = _ptBins[i+1], centMin = hiBinMin/2, centMax = hiBinMax/2;
	  seldata = Form("%s&&hiBin>=%f&&hiBin<=%f",seldataref.Data(),hiBinMin,hiBinMax);
	  selmceff = Form("%s&&hiBin>=%f&&hiBin<=%f",selmceffref.Data(),hiBinMin,hiBinMax);
	  selmcgen = Form("%s&&hiBin>=%f&&hiBin<=%f",selmcgenref.Data(),hiBinMin,hiBinMax);
	  selmc = Form("%s&&hiBin>=%f&&hiBin<=%f",selmcref.Data(),hiBinMin,hiBinMax);
	}

      std::cout<<"using centrality = "<<centMin<<"-"<<centMax<<"%"<<std::endl;;

      TF1* f = fit(nt,ntMC,ptBinsInc[0],ptBinsInc[1],isMC,isPbPb,totalmass,centMin,centMax,npfit);
      hMean->SetBinContent(i+1,f->GetParameter(1));
      hMean->SetBinError(i+1,f->GetParError(1));  
      //printf("centrality bin %.0f-%.0f of p_t bin %.0f-%.0f     Yield: %f     YieldErr: %f     RelYieldErr: %f\n", _ptBins[i]/2, _ptBins[i+1]/2, ptBinsInc[0], ptBinsInc[1], yield, yieldErr, yieldErr/yield);
      yieldErr = yieldErr*_ErrCor;
      hPt->SetBinContent(i+1,yield/(ptBinsInc[1]-ptBinsInc[0]));
      hPt->SetBinError(i+1,yieldErr/(ptBinsInc[1]-ptBinsInc[0]));
    }  

  if(isPbPb)
    {
      seldata = Form("%s&&hiBin>=%f&&hiBin<=%f",seldataref.Data(),centmin*2,centmax*2);
      selmceff = Form("%s&&hiBin>=%f&&hiBin<=%f",selmceffref.Data(),centmin*2,centmax*2);
      selmcgen = Form("%s&&hiBin>=%f&&hiBin<=%f",selmcgenref.Data(),centmin*2,centmax*2);
      selmc = Form("%s&&hiBin>=%f&&hiBin<=%f",selmcref.Data(),centmin*2,centmax*2);
    }

  ntMC->Project("hPtMC","hiBin",TCut(weight)*(TCut(selmceff.Data())&&"(Bgen==23333)"));
  divideBinWidth(hPtMC);
  ntMC->Project("hPtRecoTruth","hiBin",TCut(selmceff.Data())&&"(Bgen==23333)");
  divideBinWidth(hPtRecoTruth);
  ntGen->Project("hPtGen","hiBin",TCut(weightgen)*(TCut(selmcgen.Data())));
  divideBinWidth(hPtGen);

  TCanvas* cPt =  new TCanvas("cPt","",600,600);
  cPt->SetLogy();
  hPt->SetXTitle("Centrality");
  hPt->SetYTitle("Uncorrected dN(D^{0})/dp_{T}");
  hPt->Sumw2();
  hPt->Draw();
  if(isMC==1)
    {
      hPtMC->Draw("same hist");
      TLegend* legPt = myLegend(0.55,0.80,0.90,0.94);
      legPt->AddEntry(hPt,"Signal extraction","pl");
      legPt->AddEntry(hPtMC,"Matched reco","lf");
      legPt->Draw("same");  
    }
  hPtMC->Sumw2();
  TH1D* hEff = (TH1D*)hPtMC->Clone("hEff");
  hEff->SetTitle(";B^{+} p_{T} (GeV/c);Efficiency");
  hEff->Sumw2();
  hEff->Divide(hPtGen);
  
  TCanvas* cEff = new TCanvas("cEff","",600,600);
  hEff->Draw();
  
  TH1D* hPtCor = (TH1D*)hPt->Clone("hPtCor");
  hPtCor->SetTitle(";B^{+} Centrality;Corrected dN(B^{+})/dp_{T}");
  hPtCor->Divide(hEff);
  TCanvas* cPtCor=  new TCanvas("cCorResult","",600,600);
  cPtCor->SetLogy();
  hPtCor->Draw();
  if(isMC==1)
    {
      hPtGen->Draw("same hist");
      TLegend* legPtCor = myLegend(0.55,0.80,0.90,0.94);
      legPtCor->AddEntry(hPtCor,"Corrected signal","pl");
      legPtCor->AddEntry(hPtGen,"Generated B^{+}","lf");
      legPtCor->Draw("same");  
    }

  TH1D* hPtSigma= (TH1D*)hPtCor->Clone("hPtSigma");
  hPtSigma->SetTitle(";B^{+} p_{T} (GeV/c);d#sigma(B^{+})/dp_{T} (pb/GeV)");
  hPtSigma->Scale(1./(2*luminosity*BRchain));
  TCanvas* cPtSigma=  new TCanvas("cPtSigma","",600,600);
  cPtSigma->SetLogy();
  hPtSigma->Draw();

  TString outputf;
  outputf = Form("%s",outputfile.Data());
  
  TFile* outf = new TFile(outputf.Data(),"recreate");
  outf->cd();
  hPt->Write();
  hEff->Write();
  hMean->Write();
  hPtGen->Write();
  hPtMC->Write();
  hPtCor->Write();
  hPtSigma->Write();
  outf->Close();

}

void clean0(TH1D* h)
{
  for (int i=1;i<=h->GetNbinsX();i++)
    {
      if(h->GetBinContent(i)==0) h->SetBinError(i,1);
    }
}

void getNPFnPar(TString npfname, float par[]){
	TFile* npf = new TFile(npfname.Data());
	TF1* f = (TF1*)npf->Get("f1");
	par[0] = f->GetParameter(1);
	par[1] = f->GetParameter(2);
}

   //TF1 *fit(TTree *nt, TTree *ntMC, Double_t ptmin, Double_t ptmax, int isMC, bool isPbPb, TF1* &total, Float_t centmin, Float_t centmax, float NPpar[])
   TF1 *fit(TTree *nt, TTree *ntMC, Double_t ptmin, Double_t ptmax, int isMC, bool isPbPb, TF1* &total, Float_t centMin, Float_t centMax, TString npfit)
{
   //cout<<cut.Data()<<endl;
   static Int_t count=0;
   count++;
   TCanvas* c= new TCanvas(Form("c%d",count),"",600,600);
   TH1D* h = new TH1D(Form("h-%d",count),"",nbinsmasshisto,minhisto,maxhisto);
   TH1D* hMCSignal = new TH1D(Form("hMCSignal-%d",count),"",nbinsmasshisto,minhisto,maxhisto);

   //TString iNP="7.26667e+00*TMath::Gaus(x,5.10472e+00,2.63158e-02)/(sqrt(2*3.14159)*2.63158e-02)+4.99089e+01*TMath::Gaus(x,4.96473e+00,9.56645e-02)/(sqrt(2*3.14159)*9.56645e-02)+3.94417e-01*(3.74282e+01*TMath::Gaus(x,5.34796e+00,3.11510e-02)+1.14713e+01*TMath::Gaus(x,5.42190e+00,1.00544e-01))";
   //TString iNP=Form("TMath::Erf((x-%f)/%f)+1", NPpar[0], NPpar[1]);
   TString iNP = npfit;
   TF1* f = new TF1(Form("f%d",count),"[0]*([7]*TMath::Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*TMath::Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[5]*("+iNP+")");
   f->SetNpx(5000);
   f->SetLineWidth(5);
   
   if(isMC==1) ntMC->Project(Form("h-%d",count),"Bmass",Form("%s*(%s&&Bpt>%f&&Bpt<%f)*(1/%s)",weight.Data(),seldata.Data(),ptmin,ptmax,weightdata.Data()));
   else nt->Project(Form("h-%d",count),"Bmass",Form("(%s&&Bpt>%f&&Bpt<%f)*(1/%s)",seldata.Data(),ptmin,ptmax,weightdata.Data()));
   ntMC->Project(Form("hMCSignal-%d",count),"Bmass",Form("%s*(%s&&Bgen==23333&&Bpt>%f&&Bpt<%f)*(1/%s)",weight.Data(),selmc.Data(),ptmin,ptmax,weightdata.Data()));

   clean0(h);

   f->SetParLimits(4,-1e5,1e5);
   f->SetParLimits(2,0.01,0.05);
   f->SetParLimits(8,0.01,0.05);
   f->SetParLimits(7,0,1);
   f->SetParLimits(5,0,1e4);
   f->SetParLimits(0,0,1e5);

   //Do the signal fit first

   f->SetParameter(0,setparam0);
   f->SetParameter(1,setparam1);
   f->SetParameter(2,setparam2);
   f->SetParameter(8,setparam3);
   f->FixParameter(1,fixparam1);
   f->FixParameter(5,0);
   //f->FixParameter(3,0);
   //f->FixParameter(4,0);
   if(weightdata != "1"){
     int maxb = h->GetMaximumBin();
     double _max = h->GetBinContent(maxb);
     double _maxE = h->GetBinError(maxb);
     _ErrCor = (_maxE/_max)/(1/sqrt(_max));
     f->SetParLimits(0,0,1e5);
     f->SetParLimits(4,-1e5,1e5);
     f->SetParLimits(5,0,1e4);
   }
   h->GetEntries();

   hMCSignal->Fit(Form("f%d",count),"q","",minhisto,maxhisto);
   hMCSignal->Fit(Form("f%d",count),"q","",minhisto,maxhisto);
   //f->ReleaseParameter(1);
   hMCSignal->Fit(Form("f%d",count),"L q","",minhisto,maxhisto);
   hMCSignal->Fit(Form("f%d",count),"L q","",minhisto,maxhisto);
   hMCSignal->Fit(Form("f%d",count),"L q","",minhisto,maxhisto);
   hMCSignal->Fit(Form("f%d",count),"L m","",minhisto,maxhisto);

   //Fix the signal fit and do the background fits

   f->FixParameter(1,f->GetParameter(1));
   f->FixParameter(2,f->GetParameter(2));
   f->FixParameter(7,f->GetParameter(7));
   f->FixParameter(8,f->GetParameter(8));

   //f->ReleaseParameter(3);
   //f->ReleaseParameter(4);
   f->ReleaseParameter(5);
   f->SetParLimits(5,0,1000);

   printf("Fixed para.:\n");
   printf("%f, %f, %f\n", f->GetParameter(2), f->GetParameter(7), f->GetParameter(8));
   h->Fit(Form("f%d",count),"q","",minhisto,maxhisto);
   h->Fit(Form("f%d",count),"q","",minhisto,maxhisto);
   f->ReleaseParameter(1);
   h->Fit(Form("f%d",count),"L q","",minhisto,maxhisto);
   h->Fit(Form("f%d",count),"L q","",minhisto,maxhisto);
   h->Fit(Form("f%d",count),"L q","",minhisto,maxhisto);
   h->Fit(Form("f%d",count),"L m","",minhisto,maxhisto);
   if(weightdata != "1"){
     h->Fit(Form("f%d",count),"q","",minhisto,maxhisto);
     h->Fit(Form("f%d",count),"q","",minhisto,maxhisto);
     h->Fit(Form("f%d",count),"m","",minhisto,maxhisto);
   }

   TF1 *background = new TF1(Form("background%d",count),"[0]+[1]*x");
   background->SetParameter(0,f->GetParameter(3));
   background->SetParameter(1,f->GetParameter(4));
   background->SetLineColor(4);
   background->SetRange(minhisto,maxhisto);
   background->SetLineStyle(2);
   
   TF1 *Bkpi = new TF1(Form("fBkpi%d",count),"[0]*("+iNP+")");
   Bkpi->SetParameter(0,f->GetParameter(5));
   Bkpi->SetLineColor(kGreen+1);
   Bkpi->SetRange(minhisto,maxhisto);
   Bkpi->SetLineStyle(1);
   Bkpi->SetFillStyle(3004);
   Bkpi->SetFillColor(kGreen+1);

   TF1 *mass = new TF1(Form("fmass%d",count),"[0]*([3]*TMath::Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[3])*TMath::Gaus(x,[1],[4])/(sqrt(2*3.14159)*[4]))");
   mass->SetParameters(f->GetParameter(0),f->GetParameter(1),f->GetParameter(2),f->GetParameter(7),f->GetParameter(8));
   mass->SetParError(0,f->GetParError(0));
   mass->SetParError(1,f->GetParError(1));
   mass->SetParError(2,f->GetParError(2));
   mass->SetParError(7,f->GetParError(7));
   mass->SetParError(8,f->GetParError(8));
   mass->SetLineColor(2);

   //h->SetXTitle("m_{#mu#muK} (GeV/c^{2})");
  h->SetXTitle("m_{B} (GeV/c^{2})");

  double yunit = 1000./nbinsmasshisto;
  h->SetYTitle(Form("Events / (%.0f MeV/c^{2})",yunit));//basically, the unit is equivalent to binwidthmass in MeV

  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  h->SetAxisRange(0,h->GetMaximum()*1.4*1.2,"Y");
  h->GetXaxis()->SetTitleOffset(1.3);
  h->GetYaxis()->SetTitleOffset(1.8);
  h->GetXaxis()->SetLabelOffset(0.007);
  h->GetYaxis()->SetLabelOffset(0.007);
  h->GetXaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetTitleSize(0.045);
  h->GetXaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleFont(42);
  h->GetXaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetLabelSize(0.04);
  h->SetMarkerSize(0.8);
  h->SetMarkerStyle(20);
  h->SetStats(0);
  h->Draw("e");
  Bkpi->Draw("same");
  background->Draw("same");   
  mass->SetRange(minhisto,maxhisto);
  mass->Draw("same");
  mass->SetLineStyle(2);
  mass->SetFillStyle(3004);
  mass->SetFillColor(2);
  f->Draw("same");
  c->RedrawAxis();

  yield = mass->Integral(minhisto,maxhisto)/binwidthmass;
  yieldErr = mass->Integral(minhisto,maxhisto)/binwidthmass*mass->GetParError(0)/mass->GetParameter(0);
  printf("centrality bin %.0f-%.0f     yield: %f     yieldErr: %f\n", centMin, centMax, yield, yieldErr);

  /*
  Double_t Signal = mass->Integral(5.19932,5.35932)/binwidthmass; //B+ mass_pdg=5.27932GeV, signal region = pm 0.08GeV
  Double_t Bkg_comb = background->Integral(5.19932,5.35932)/binwidthmass;
  Double_t Bkg_nonprompt = Bkpi->Integral(5.19932,5.35932)/binwidthmass;
  printf("p_t bin %.0f-%.0f S: %f B(comb): %f B(comb+np): %f\n", ptmin, ptmax, Signal, Bkg_comb, Bkg_comb+Bkg_nonprompt);
  printf("p_t bin %.0f-%.0f sig(comb): %f sig(comb+np): %f\n", ptmin, ptmax, Signal/sqrt(Signal+Bkg_comb), Signal/sqrt(Signal+Bkg_comb+Bkg_nonprompt));
  */

  TLegend* leg = new TLegend(0.65,0.58,0.82,0.88,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->AddEntry(h,"Data","pl");
  leg->AddEntry(f,"Fit","l");
  leg->AddEntry(mass,"B^{+} Signal","f");
  leg->AddEntry(background,"Combinatorial","l");
  leg->AddEntry(Bkpi,"B #rightarrow J/#psi X","l");
  leg->Draw("same");

  TLatex* texYield = new TLatex(0.58,0.55,Form("Yield:%.2f#pm%.2f",yield,yieldErr));
  texYield->SetNDC();
  texYield->SetTextAlign(12);
  texYield->SetTextSize(0.035);
  texYield->SetTextFont(42);
  texYield->Draw();

  TLatex* texChi = new TLatex(0.58,0.515,Form("#chi^{2}/nDOF:%.2f/%d=%.2f",f->GetChisquare(),f->GetNDF(),f->GetChisquare()/f->GetNDF()));
  texChi->SetNDC();
  texChi->SetTextAlign(12);
  texChi->SetTextSize(0.035);
  texChi->SetTextFont(42);
  texChi->Draw();
  printf("NDF: %d, chi2: %f, prob: %f\n", f->GetNDF(), f->GetChisquare(), f->GetProb());

  TLatex* texCms = new TLatex(0.18,0.93, "#scale[1.25]{CMS}");
  texCms->SetNDC();
  texCms->SetTextAlign(12);
  texCms->SetTextSize(0.04);
  texCms->SetTextFont(42);
  texCms->Draw();

  TLatex* texCol;
  if(collisionsystem=="pp"||collisionsystem=="PP") texCol= new TLatex(0.96,0.93, Form("%s #sqrt{s_{NN}} = 5.02 TeV","pp"));
  else texCol= new TLatex(0.96,0.93, Form("%s #sqrt{s_{NN}} = 5.02 TeV","PbPb"));
  texCol->SetNDC();
  texCol->SetTextAlign(32);
  texCol->SetTextSize(0.04);
  texCol->SetTextFont(42);
  texCol->Draw();

  TLatex* tex;
  TString texper="%";
  tex = new TLatex(0.22,0.78,Form("Cent. %.0f-%.0f%s",centMin,centMax,texper.Data()));
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();
  
  tex = new TLatex(0.22,0.73,Form("%.0f<p_{T}<%.0f GeV/c",ptBinsInc[0],ptBinsInc[1]));
  tex->SetNDC();
  tex->SetTextColor(1);
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();

  tex = new TLatex(0.22,0.83,"B^{+}+B^{-} |y| < 2.4");
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();
     
  total=f;

  TString _postfix = "%";
  if(weightdata!="1") _postfix = "_EFFCOR";
  if(isPbPb && isMC==0)
    { 
      c->SaveAs(Form("plotFits/AN_Cent_data_PbPb_%.0f-%.0f.png",centMin,centMax));
      c->SaveAs(Form("plotFits/AN_Cent_data_PbPb_%.0f-%.0f.pdf",centMin,centMax));
    }
  else if(isPbPb && isMC==1) 
    {
      c->SaveAs(Form("plotFits/AN_Cent_mc_PbPb_%.0f%.0f.png",centMin,centMax));
      c->SaveAs(Form("plotFits/AN_Cent_mc_PbPb_%.0f%.0f.pdf",centMin,centMax));
    }
  else if(!isPbPb && isMC==0) 
    {
      c->SaveAs("plotFits/Cent_data_pp.png");
      c->SaveAs("plotFits/Cent_data_pp.pdf");
    }
  else 
    {
      c->SaveAs("plotFits/Cent_mc_pp.png");
      c->SaveAs("plotFits/Cent_mc_pp.pdf");
    }
  return mass;
}


int main(int argc, char *argv[])
{
  if(argc==16)
    {
      fitBCent(atoi(argv[1]),argv[2], argv[3], argv[4], argv[5], argv[6], atoi(argv[7]), atof(argv[8]), atoi(argv[9]),argv[10],argv[11],argv[12],atoi(argv[13]),atof(argv[14]),atof(argv[15]));
      return 0;
    }
  else if(argc==14)
    {
      fitBCent(atoi(argv[1]), argv[2], argv[3], argv[4], argv[5], argv[6], atoi(argv[7]), atof(argv[8]), atoi(argv[9]),argv[10],argv[11], argv[12], atoi(argv[13]));
      return 0;
    }
  else
    {
      std::cout << "Wrong number of inputs" << std::endl;
      return 1;
    }
}



double ErrorOnSigma(double width, double errwidth, double smear, double errsmearing){
double squarederroronsigma=(1+smear)*(1+smear)*errwidth*errwidth+width*width*errsmearing*errsmearing;
double erroronsigma=TMath::Sqrt(squarederroronsigma);
return erroronsigma;

}

