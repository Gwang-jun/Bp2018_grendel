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
TString seldata;
TString selmc;
TString selmceff;
TString selmcgen;
TString collisionsystem;
Float_t hiBinMin,hiBinMax,centMin,centMax;
double _ErrCor=1;

int _nBins = nBinsY;
double *_ptBins = ptBinsY;

Double_t yield;
Double_t yieldErr;

void fitBY(int usePbPb=0, TString inputdata="" , TString inputmc="", TString trgselection="",  TString cut="", TString cutmcgen="", int isMC=0, Double_t luminosity=1., int doweight=0, TString collsyst="PbPb", TString outputfile="", TString npfit="0", int doDataCor = 0, Float_t centmin=0., Float_t centmax=90.)
{
  collisionsystem=collsyst;
  if(collisionsystem=="ppInc"||collisionsystem=="PbPbInc"){
    _nBins = nBinsInc;
    _ptBins = ptBinsInc;
  }
  
  hiBinMin = centmin*2;
  hiBinMax = centmax*2;
  centMin = centmin;
  centMax = centmax;
  
  double ErrorOnSigma(double width, double errwidth, double smear, double errsmearing);
  
  if (!(usePbPb==1||usePbPb==0)) std::cout<<"ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!, you are using a non valid isPbPb option"<<std::endl;
  bool isPbPb=(bool)(usePbPb);

  if(!isPbPb)
    {
      seldata = Form("%s&&%s&&Bpt>%f&&Bpt<%f",trgselection.Data(),cut.Data(),ptBinsInc[0],ptBinsInc[1]);
      selmceff = Form("%s&&%s&&Bpt>%f&&Bpt<%f",trgselection.Data(),cut.Data(),ptBinsInc[0],ptBinsInc[1]);
      selmcgen = Form("%s&&Gpt>%f&&Gpt<%f",cutmcgen.Data(),ptBinsInc[0],ptBinsInc[1]);
      selmc = Form("%s&&%s&&Bpt>%f&&Bpt<%f",trgselection.Data(),cut.Data(),ptBinsInc[0],ptBinsInc[1]);
    }
  else
    {
      seldata = Form("%s&&%s&&hiBin>=%f&&hiBin<=%f&&Bpt>%f&&Bpt<%f",trgselection.Data(),cut.Data(),hiBinMin,hiBinMax,ptBinsInc[0],ptBinsInc[1]);
      selmceff = Form("%s&&%s&&hiBin>=%f&&hiBin<=%f&&Bpt>%f&&Bpt<%f",trgselection.Data(),cut.Data(),hiBinMin,hiBinMax,ptBinsInc[0],ptBinsInc[1]);
      selmcgen = Form("%s&&hiBin>=%f&&hiBin<=%f&&Gpt>%f&&Gpt<%f",cutmcgen.Data(),hiBinMin,hiBinMax,ptBinsInc[0],ptBinsInc[1]);
      selmc = Form("%s&&%s&&hiBin>=%f&&hiBin<=%f&&Bpt>%f&&Bpt<%f",trgselection.Data(),cut.Data(),hiBinMin,hiBinMax,ptBinsInc[0],ptBinsInc[1]);
    }

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
  TF1* fit (TTree* nt, TTree* ntMC, double ptmin, double ptmax, int isMC,bool, TF1* &total,Float_t centmin, Float_t centmax, TString npfit);
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
      weightgen="pthatweight*(3.00448277-0.35865276*Gpt+0.01997413*Gpt*Gpt-0.00042585*Gpt*Gpt*Gpt+0.00000315*Gpt*Gpt*Gpt*Gpt)";
      weight="pthatweight*Ncoll*(TMath::Gaus(PVz,0.427450,4.873825)/(sqrt(2*3.14159)*4.873825))/(TMath::Gaus(PVz,0.909938,4.970989)/(sqrt(2*3.14159)*4.970989))*(3.00448277-0.35865276*Bgenpt+0.01997413*Bgenpt*Bgenpt-0.00042585*Bgenpt*Bgenpt*Bgenpt+0.00000315*Bgenpt*Bgenpt*Bgenpt*Bgenpt)";
      //weightgen="pthatweight*((3.506006+0.963473*Gpt+-0.258731*Gpt*Gpt)*TMath::Exp(-0.386065*Gpt)+1.139897)";
      //weight="pthatweight*Ncoll*(TMath::Gaus(PVz,0.427450,4.873825)/(sqrt(2*3.14159)*4.873825))/(TMath::Gaus(PVz,0.909938,4.970989)/(sqrt(2*3.14159)*4.970989))*((3.506006+0.963473*Bgenpt+-0.258731*Bgenpt*Bgenpt)*TMath::Exp(-0.386065*Bgenpt)+1.139897)";
    }  
  std::cout<<"we are using weight="<<weight<<std::endl;
  std::cout<<"we are using weightdata="<<weightdata<<std::endl;
  
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
  ntMC->AddFriend("ntGen");
  ntMC->AddFriend("mvaTree");
  */
  
  TH1D* hPt = new TH1D("hPt","",_nBins,_ptBins);
  TH1D* hPtRecoTruth = new TH1D("hPtRecoTruth","",_nBins,_ptBins);
  TH1D* hPtMC = new TH1D("hPtMC","",_nBins,_ptBins);
  TH1D* hPtGen = new TH1D("hPtGen","",_nBins,_ptBins);

  TH1D* hMean = new TH1D("hMean","",_nBins,_ptBins);                       
  TH1D* hSigmaGaus1 = new TH1D("hSigmaGaus1","",_nBins,_ptBins); 
  TH1D* hSigmaGaus2 = new TH1D("hSigmaGaus2","",_nBins,_ptBins); 
  TF1 *totalmass;

  for(int i=0;i<_nBins;i++)
    {
      //TF1* f = fit(nt,ntMC,_ptBins[i],_ptBins[i+1],isMC,isPbPb, totalmass,centmin, centmax, NPpar);
      TF1* f = fit(nt,ntMC,_ptBins[i],_ptBins[i+1],isMC,isPbPb, totalmass,centmin, centmax, npfit);
      hMean->SetBinContent(i+1,f->GetParameter(1));
      hMean->SetBinError(i+1,f->GetParError(1));  
      yieldErr = yieldErr*_ErrCor;
      hPt->SetBinContent(i+1,yield/(_ptBins[i+1]-_ptBins[i]));
      hPt->SetBinError(i+1,yieldErr/(_ptBins[i+1]-_ptBins[i]));
    }  

  ntMC->Project("hPtMC","TMath::Abs(By)",TCut(weight)*(TCut(selmceff.Data())&&"(Bgen==23333)"));
  divideBinWidth(hPtMC);
  ntMC->Project("hPtRecoTruth","TMath::Abs(By)",TCut(selmceff.Data())&&"(Bgen==23333)");
  divideBinWidth(hPtRecoTruth);
  ntGen->Project("hPtGen","Gy",TCut(weightgen)*(TCut(selmcgen.Data())));
  divideBinWidth(hPtGen);

  TCanvas* cPt =  new TCanvas("cPt","",600,600);
  cPt->SetLogy();
  hPt->SetXTitle("B^{+} y");
  hPt->SetYTitle("Uncorrected dN(B^{+})/dy");
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
  hEff->SetTitle(";B^{+} y;Efficiency");
  hEff->Sumw2();
  hEff->Divide(hPtGen);
  TCanvas* cEff = new TCanvas("cEff","",600,600);
  hEff->Draw();
  
  TH1D* hPtCor = (TH1D*)hPt->Clone("hPtCor");
  hPtCor->SetTitle(";B^{+} y;Corrected dN(B^{+})/dy");
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
  hPtSigma->SetTitle(";B^{+} y;d#sigma(B^{+})/dy (pb)");
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

//TF1 *fit(TTree *nt, TTree *ntMC, Double_t ptmin, Double_t ptmax, int isMC,bool isPbPb,TF1* &total,Float_t centmin, Float_t centmax, float NPpar[])
TF1 *fit(TTree *nt, TTree *ntMC, Double_t ptmin, Double_t ptmax, int isMC,bool isPbPb,TF1* &total,Float_t centmin, Float_t centmax, TString npfit)
{
   static Int_t count=0;
   count++;
   TCanvas* c= new TCanvas(Form("c%d",count),"",600,600);
   TH1D* h = new TH1D(Form("h-%d",count),"",nbinsmasshisto,minhisto,maxhisto);
   TH1D* hMCSignal = new TH1D(Form("hMCSignal-%d",count),"",nbinsmasshisto,minhisto,maxhisto);

   //TString iNP="7.26667e+00*Gaus(x,5.10472e+00,2.63158e-02)/(sqrt(2*3.14159)*2.63158e-02)+4.99089e+01*Gaus(x,4.96473e+00,9.56645e-02)/(sqrt(2*3.14159)*9.56645e-02)+3.94417e-01*(3.74282e+01*Gaus(x,5.34796e+00,3.11510e-02)+1.14713e+01*Gaus(x,5.42190e+00,1.00544e-01))";
   //TString iNP=Form("TMath::Erf((x-%f)/%f)+1", NPpar[0], NPpar[1]);
   TString iNP = npfit;
   TF1* f = new TF1(Form("f%d",count),"[0]*([7]*TMath::Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*TMath::Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[5]*("+iNP+")");

   if(isMC==1) nt->Project(Form("h-%d",count),"Bmass",Form("%s*(%s&&TMath::Abs(By)>%f&&TMath::Abs(By)<%f)*(1/%s)","1",seldata.Data(),ptmin,ptmax,weightdata.Data()));   
   else nt->Project(Form("h-%d",count),"Bmass",Form("(%s&&TMath::Abs(By)>%f&&TMath::Abs(By)<%f)*(1/%s)",seldata.Data(),ptmin,ptmax,weightdata.Data()));   
   ntMC->Project(Form("hMCSignal-%d",count),"Bmass",Form("%s&&Bgen==23333&&TMath::Abs(By)>%f&&TMath::Abs(By)<%f",selmc.Data(),ptmin,ptmax));
   clean0(h);
  
   h->Draw();
   f->SetParLimits(4,-1000,1000);
   f->SetParLimits(2,0.01,0.05);
   f->SetParLimits(8,0.01,0.05);
   f->SetParLimits(7,0,1);
   f->SetParLimits(5,0,1000);

   f->SetParameter(0,setparam0);
   f->SetParameter(1,setparam1);
   f->SetParameter(2,setparam2);
   f->SetParameter(8,setparam3);
   f->FixParameter(1,fixparam1);
   f->FixParameter(5,0);
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
   f->ReleaseParameter(1);
   hMCSignal->Fit(Form("f%d",count),"L q","",minhisto,maxhisto);
   hMCSignal->Fit(Form("f%d",count),"L q","",minhisto,maxhisto);
   hMCSignal->Fit(Form("f%d",count),"L q","",minhisto,maxhisto);
   hMCSignal->Fit(Form("f%d",count),"L m","",minhisto,maxhisto);

   f->FixParameter(1,f->GetParameter(1));
   f->FixParameter(2,f->GetParameter(2));
   f->FixParameter(7,f->GetParameter(7));
   f->FixParameter(8,f->GetParameter(8));
   f->ReleaseParameter(5);
   f->SetParLimits(5,0,1000);
   
   h->Fit(Form("f%d",count),"q","",minhisto,maxhisto);
   h->Fit(Form("f%d",count),"q","",minhisto,maxhisto);
   f->ReleaseParameter(1);
   h->Fit(Form("f%d",count),"L q","",minhisto,maxhisto);
   h->Fit(Form("f%d",count),"L q","",minhisto,maxhisto);
   h->Fit(Form("f%d",count),"L q","",minhisto,maxhisto);
   h->Fit(Form("f%d",count),"L m","",minhisto,maxhisto);

   h->SetMarkerSize(0.8);
   h->SetMarkerStyle(20);

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

   yield = mass->Integral(minhisto,maxhisto)/binwidthmass;
   yieldErr = mass->Integral(minhisto,maxhisto)/binwidthmass*mass->GetParError(0)/mass->GetParameter(0);
   printf("Y bin %.1f-%.1f     yield: %f     yieldErr: %f\n", ptmin, ptmax, yield, yieldErr);

  h->SetXTitle("m_{#mu#muK} (GeV/c^{2})");
  h->SetYTitle("Entries / (20 MeV/c^{2})");
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

   Double_t yield = mass->Integral(minhisto,maxhisto)/binwidthmass;
   Double_t yieldErr = mass->Integral(minhisto,maxhisto)/binwidthmass*mass->GetParError(0)/mass->GetParameter(0);

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

  TLatex* texChi = new TLatex(0.58,0.55, Form("#chi^{2}/nDOF: %.2f/%d = %.2f", f->GetChisquare(), f->GetNDF(), f->GetChisquare()/f->GetNDF()));
  texChi->SetNDC();
  texChi->SetTextAlign(12);
  texChi->SetTextSize(0.03);
  texChi->SetTextFont(42);
  texChi->Draw();

  TLatex* texCms = new TLatex(0.18,0.93, "#scale[1.25]{CMS} Preliminary");
  texCms->SetNDC();
  texCms->SetTextAlign(12);
  texCms->SetTextSize(0.04);
  texCms->SetTextFont(42);
  texCms->Draw();

  TLatex* texCol;
  if(collisionsystem=="pp"||collisionsystem=="PP"||collisionsystem=="ppInc"||collisionsystem=="PbPbInc") texCol= new TLatex(0.96,0.93, Form("%s #sqrt{s_{NN}} = 5.02 TeV","pp"));
  else texCol= new TLatex(0.96,0.93, Form("%s #sqrt{s_{NN}} = 5.02 TeV","PbPb"));
  texCol->SetNDC();
  texCol->SetTextAlign(32);
  texCol->SetTextSize(0.04);
  texCol->SetTextFont(42);
  texCol->Draw();

  TLatex* tex;

  tex = new TLatex(0.22,0.78,Form("%.1f < |y| < %.1f",ptmin,ptmax));
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();
  
  if(centMax>0){
  TString texper="%";
  tex = new TLatex(0.22,0.71,Form("Cent. %.0f-%.0f%s",centMin,centMax,texper.Data()));//0.2612903,0.8425793
  tex->SetNDC();
  tex->SetTextColor(1);
  tex->SetTextFont(42);
  tex->SetTextSize(0.045);
  tex->SetLineWidth(2);
  tex->Draw();
  }

  tex = new TLatex(0.22,0.83,Form("%.0f < p_{T} < %.0f (GeV/c)",ptBinsInc[0], ptBinsInc[1]));
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();
     
   total=f;
   
  //if(!isPbPb) c->SaveAs(Form("plotFits/BMass%s_%d.pdf",collisionsystem.Data(),count));
  //else c->SaveAs(Form("plotFits/BMass%s_%.0f_%.0f_%d.pdf",collisionsystem.Data(),centMin,centMax,count));

  TString _postfix = "";
  if(weightdata!="1") _postfix = "_EFFCOR";
  if(isPbPb && isMC==0) 
      c->SaveAs(Form("plotFits/data_PbPb_Y_%.1f_%.1f%s.pdf",ptmin,ptmax,_postfix.Data()));
  else if(isPbPb && isMC==1) 
      c->SaveAs(Form("plotFits/mc_PbPb_Y_%.1f_%.1f%s.pdf",ptmin,ptmax,_postfix.Data()));
  else if(!isPbPb && isMC==0) 
      c->SaveAs(Form("plotFits/data_pp_Y_%.1f_%.1f%s.pdf",ptmin,ptmax,_postfix.Data()));
  else 
      c->SaveAs(Form("plotFits/mc_pp_Y_%.1f_%.1f%s.pdf",ptmin,ptmax,_postfix.Data()));

  return mass;
}


int main(int argc, char *argv[])
{
  if(argc==16)
    {
      fitBY(atoi(argv[1]),argv[2], argv[3], argv[4], argv[5], argv[6], atoi(argv[7]), atof(argv[8]), atoi(argv[9]),argv[10],argv[11],argv[12],atoi(argv[13]),atof(argv[14]),atof(argv[15]));
      return 0;
    }
  else if(argc==14)
    {
      fitBY(atoi(argv[1]), argv[2], argv[3], argv[4], argv[5], argv[6], atoi(argv[7]), atof(argv[8]), atoi(argv[9]),argv[10],argv[11], argv[12], atoi(argv[13]));
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

