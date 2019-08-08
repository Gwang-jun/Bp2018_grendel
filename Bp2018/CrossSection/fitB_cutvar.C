#include "uti.h"
#include "parameters_cutvar.h"
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

int _nBins = nBins;
double *_ptBins = ptBins;

Double_t yield;
Double_t yieldErr;
Int_t VARNUM;
Int_t stage;

void fitB(int usePbPb=0, TString inputdata="" , TString inputmc="", TString trgselection="1",  TString cut="", TString cutmcgen="", int isMC=0, Double_t luminosity=1., int doweight=0, TString collsyst="PbPb", TString outputfile="", TString npfit="0", int doDataCor = 0, Float_t centmin=0., Float_t centmax=100., int varnum=0)
{
  gStyle->SetOptStat(0);

  VARNUM = varnum;

  cutvarname[0] = "dls3D";
  cutvarname[1] = "costheta";
  cutvarname[2] = "dxysig";
  cutvarname[3] = "trkpt";
  cutvarname[4] = "chisq";
  //cutvarname[5] = "dzsig";

  cutvar[0] = "(BsvpvDistance/BsvpvDisErr)";
  cutvar[1] = "cos(Bdtheta)";
  cutvar[2] = "TMath::Abs(Btrk1Dxy1/Btrk1DxyError1)";
  cutvar[3] = "Btrk1Pt";
  cutvar[4] = "Bchi2cl";
  //cutvar[5] = "TMath::Abs(Btrk1Dz1/Btrk1DzError1)";

  cutspacing = (cutvarmax[varnum]-cutvarmin[varnum])/(ncutvar-1);

  hiBinMin = centmin*2;
  hiBinMax = centmax*2;
  centMin = centmin;
  centMax = centmax;

  double ErrorOnSigma(double width, double errwidth, double smear, double errsmearing);

  if (!(usePbPb==1||usePbPb==0)) std::cout<<"ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!, you are using a non valid isPbPb option"<<std::endl;
  bool isPbPb=(bool)(usePbPb);

  gStyle->SetTextSize(0.05);
  gStyle->SetTextFont(42);
  gStyle->SetPadRightMargin(0.043);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.145);
  gStyle->SetTitleX(.0f);

  void clean0 (TH1D* h);
  void getNPFnPar(TString npfname, float par[]);
  TF1* fit (TTree* nt, TTree* ntMC, double ptmin, double ptmax, int isMC,bool, TF1* &total,Float_t centmin, Float_t centmax, TString npfit);
 
  weightdata="1";
  if(!isPbPb)
    {
      weightgen="pthatweight*(0.599212+-0.020703*Gpt+0.003143*Gpt*Gpt+-0.000034*Gpt*Gpt*Gpt)*((1.055564*TMath::Exp(-0.001720*(PVz+2.375584)*(PVz+2.375584))))";
      weight="pthatweight*(0.599212+-0.020703*Bgenpt+0.003143*Bgenpt*Bgenpt+-0.000034*Bgenpt*Bgenpt*Bgenpt)*(1.055564*TMath::Exp(-0.001720*(PVz+2.375584)*(PVz+2.375584)))";
    }
  else
    {
      weightgen="pthatweight*((2.907795+-0.436572*Gpt+0.006372*Gpt*Gpt)*TMath::Exp(-0.157563*Gpt)+1.01308)";
      weight="pthatweight*Ncoll*(TMath::Gaus(PVz,0.427450,4.873825)/(sqrt(2*3.14159)*4.873825))/(TMath::Gaus(PVz,0.909938,4.970989)/(sqrt(2*3.14159)*4.970989))*((2.907795+-0.436572*Bgenpt+0.006372*Bgenpt*Bgenpt)*TMath::Exp(-0.157563*Bgenpt)+1.01308)";
      weightgen="pthatweight*(0.889175+0.000791*Gpt+0.000015*Gpt*Gpt)";
      weight="pthatweight*Ncoll*(TMath::Gaus(PVz,0.427450,4.873825)/(sqrt(2*3.14159)*4.873825))/(TMath::Gaus(PVz,0.909938,4.970989)/(sqrt(2*3.14159)*4.970989))*(0.889175+0.000791*Bgenpt+0.000015*Bgenpt*Bgenpt)";
      //weightgen="pthatweight*Ncoll*(1.034350*TMath::Exp(-0.000844*(PVz+3.502992)*(PVz+3.502992)))*(0.715021+0.039896*Gpt-0.000834*Gpt*Gpt+0.000006*Gpt*Gpt*Gpt)"; // private MC Gpt
      //weight="pthatweight*Ncoll*(1.034350*TMath::Exp(-0.000844*(PVz+3.502992)*(PVz+3.502992)))*(0.715021+0.039896*Bgenpt-0.000834*Bgenpt*Bgenpt+0.000006*Bgenpt*Bgenpt*Bgenpt)"; // private MC Gpt
    }

  std::cout<<"we are using weight="<<weight<<std::endl;
  std::cout<<"we are using weightdata="<<weightdata<<std::endl;
  std::cout<<"we are using centrality="<<centmin<<"-"<<centmax<<"%"<<std::endl;

  TFile* inf = new TFile(inputdata.Data());
  TFile* infMC = new TFile(inputmc.Data());

  //For 2018 PbPb data, MC
  TTree* nt = (TTree*)inf->Get("Bfinder/ntKp");
  nt->AddFriend("hltanalysis/HltTree");
  nt->AddFriend("hiEvtAnalyzer/HiTree");
  nt->AddFriend("skimanalysis/HltTree");
  TTree* ntGen = (TTree*)infMC->Get("Bfinder/ntGen");
  ntGen->AddFriend("hltanalysis/HltTree");
  ntGen->AddFriend("hiEvtAnalyzer/HiTree");
  ntGen->AddFriend("Bfinder/ntKp"); //call PVz
  ntGen->AddFriend("skimanalysis/HltTree");
  TTree* ntMC = (TTree*)infMC->Get("Bfinder/ntKp");
  ntMC->AddFriend("hltanalysis/HltTree");
  ntMC->AddFriend("hiEvtAnalyzer/HiTree");
  ntMC->AddFriend("Bfinder/ntGen"); //call Bgen
  ntMC->AddFriend("skimanalysis/HltTree");

  TString cutforothers = "";
  for(int i=0; i<nvar; i++)
    {
      if(i!=varnum) cutforothers = Form("%s&&(%s>%f)",cutforothers.Data(),cutvar[i].Data(),cutvarref[i]);
    }

  cut = Form("%s%s",cut.Data(),cutforothers.Data());

  TString cutforref = "";
  cutforref = Form("%s&&(%s>%f)",cut.Data(),cutvar[varnum].Data(),cutvarref[varnum]);
  if(!isPbPb)
    {
      seldata = Form("%s&&%s",trgselection.Data(),cutforref.Data());
      selmc = Form("%s&&%s",trgselection.Data(),cutforref.Data());
    }
  else
    {
      seldata = Form("%s&&%s&&hiBin>=%f&&hiBin<=%f",trgselection.Data(),cutforref.Data(),hiBinMin,hiBinMax);
      selmc = Form("%s&&%s&&hiBin>=%f&&hiBin<=%f",trgselection.Data(),cutforref.Data(),hiBinMin,hiBinMax);
    }

  std::cout<<""<<std::endl;
  std::cout<<"processing Data reference (optimal cut)"<<std::endl;
  std::cout<<""<<std::endl;

  fData_ref = fit(nt,ntMC,_ptBins[0],_ptBins[1],0,isPbPb,totalmass,centmin,centmax,npfit);
  yieldData_ref = yield;
  yieldErrData_ref = yieldErr;

  std::cout<<""<<std::endl;
  std::cout<<"processing MC reference (optimal cut)"<<std::endl;
  std::cout<<""<<std::endl;

  fMC_ref = fit(ntMC,ntMC,_ptBins[0],_ptBins[1],1,isPbPb,totalmass,centmin,centmax,npfit);
  yieldMC_ref = yield;
  yieldErrMC_ref = yieldErr;

  TH1D * cutvarhis = new TH1D(Form("cutvarhis_%s",cutvarname[varnum].Data()),"",ncutvar,cutvarmin[varnum]-0.5*cutspacing,cutvarmax[varnum]+0.5*cutspacing);
  cutvarhis->GetXaxis()->SetTitle(Form("%s",cutvarname[varnum].Data()));
  cutvarhis->GetYaxis()->SetTitle("Double Ratio");
  cutvarhis->GetXaxis()->CenterTitle();
  //cutvarhis->SetTitle(Form("Cut Variation of %s",cutvarname[varnum].Data()));
  cutvarhis->SetMarkerStyle(20);
  cutvarhis->SetMarkerSize(1);
  cutvarhis->GetYaxis()->SetRangeUser(0.,2.);  

  TString cutforvar = "";

  for(int i=0;i<ncutvar;i++)
   {
     stage = i;
     
     cutforvar = Form("%s&&(%s>%f)",cut.Data(),cutvar[varnum].Data(),cutvarmin[varnum]+i*cutspacing);

     if(!isPbPb)
       {
	 seldata = Form("%s&&%s",trgselection.Data(),cutforvar.Data());
	 selmc = Form("%s&&%s",trgselection.Data(),cutforvar.Data());
       }
     else
       {
	 seldata = Form("%s&&%s&&hiBin>=%f&&hiBin<=%f",trgselection.Data(),cutforvar.Data(),hiBinMin,hiBinMax);
	 selmc = Form("%s&&%s&&hiBin>=%f&&hiBin<=%f",trgselection.Data(),cutforvar.Data(),hiBinMin,hiBinMax);
       }

     std::cout<<""<<std::endl;     
     std::cout<<"processing Data with "<<i<<"-th cut"<<std::endl;
     std::cout<<""<<std::endl;     
     
     fData[i] = fit(nt,ntMC,_ptBins[0],_ptBins[1],0,isPbPb,totalmass,centmin,centmax,npfit);
     yieldData[i] = yield;
     yieldErrData[i] = yieldErr;

     std::cout<<""<<std::endl;          
     std::cout<<"processing MC with "<<i<<"-th cut"<<std::endl;
     std::cout<<""<<std::endl;            

     fMC[i] = fit(ntMC,ntMC,_ptBins[0],_ptBins[1],1,isPbPb,totalmass,centmin,centmax,npfit);
     yieldMC[i] = yield;
     yieldErrMC[i] = yieldErr;
     
     RatioData[i] = yieldData[i]/yieldData_ref;
     RatioErrData[i] = RatioData[i]*sqrt((yieldErrData[i]/yieldData[i])*(yieldErrData[i]/yieldData[i])+(yieldErrData_ref/yieldData_ref)*(yieldErrData_ref/yieldData_ref));

     RatioMC[i] = yieldMC[i]/yieldMC_ref;
     RatioErrMC[i] = RatioMC[i]*sqrt((yieldErrMC[i]/yieldMC[i])*(yieldErrMC[i]/yieldMC[i])+(yieldErrMC_ref/yieldMC_ref)*(yieldErrMC_ref/yieldMC_ref));

     DoubleRatio[i] = RatioData[i]/RatioMC[i];
     DoubleRatioErr[i] = DoubleRatio[i]*sqrt((RatioErrData[i]/RatioData[i])*(RatioErrData[i]/RatioData[i])+(RatioErrMC[i]/RatioMC[i])*(RatioErrMC[i]/RatioMC[i]));

     std::cout<<""<<std::endl;
     std::cout<<"Double Ratio at "<<i<<"-th cut = "<<DoubleRatio[i]<<" #pm "<<DoubleRatioErr[i]<<std::endl;
     
     cutvarhis->SetBinContent(i+1,DoubleRatio[i]);
     cutvarhis->SetBinError(i+1,DoubleRatioErr[i]);
   }
  
  double fmax, fmin;
  fmax = cutvarmax[varnum];
  if(varnum==1) fmin = -1.0;
  else fmin = 0.0;
  
  TF1* f = new TF1("f",Form("1.0+[0]*(x-%f)",cutvarref[varnum]),fmin,fmax);
  f->SetParLimits(0,-100,100);  
  cutvarhis->Fit(f,"R");
  double intercept = 1.0-(f->GetParameter(0))*(cutvarref[varnum]);
  double slope = (f->GetParameter(0));
  
  std::cout<<"Linear Fit Function = "<<intercept<<"+"<<slope<<"*"<<cutvarname[varnum]<<std::endl;
  
  TFile* outputroot = new TFile(Form("plotCutVar/cutvariation_%s_%s_pt%.0f-%.0f.root",cutvarname[varnum].Data(),collsyst.Data(),_ptBins[0],_ptBins[1]),"recreate");
  outputroot->cd();
  cutvarhis->Write();
  outputroot->Close();
    
  TCanvas* c0 = new TCanvas("","",600,600);
  c0->cd();
  cutvarhis->Draw("ep");
  //f->Draw("same");
  //c0->RedrawAxis();

  TLine* vline = new TLine(cutvarref[varnum],0.0,cutvarref[varnum],1.0);
  vline->SetLineWidth(1);
  vline->SetLineStyle(2);
  vline->SetLineColor(kGreen);
  vline->Draw("same");

  TLatex* texcms = new TLatex(0.22,0.87,"CMS");
  texcms->SetNDC();
  texcms->SetTextAlign(13);
  texcms->SetTextFont(62);
  texcms->SetTextSize(0.08);
  texcms->SetLineWidth(2);
  texcms->Draw();

  TLatex* texB = new TLatex(0.22,0.73,"B^{+}+B^{-}");
  texB->SetNDC();
  texB->SetTextFont(42);
  texB->SetTextSize(0.07);
  texB->SetLineWidth(2);
  texB->Draw();

  TLatex* texCol;
  if(collisionsystem=="pp"||collisionsystem=="PP"||collisionsystem=="ppInc"||collisionsystem=="PbPbInc") texCol= new TLatex(0.945,0.94, Form("28.0 pb^{-1} (%s 5.02 TeV)","pp"));
  else texCol= new TLatex(0.94,0.93, Form("1.5 nb^{-1} (%s 5.02 TeV)","PbPb"));
  texCol->SetNDC();
  texCol->SetTextAlign(32);
  texCol->SetTextSize(0.06);
  texCol->SetTextFont(42);
  texCol->Draw();

  TLatex* tex;
  tex = new TLatex(0.49,0.845,Form("%.0f<p_{T}<%.0f GeV/c",_ptBins[0],_ptBins[1]));
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.045);
  tex->SetLineWidth(2);
  tex->Draw();

  TString texper="%";
  tex = new TLatex(0.49,0.80,Form("Cent. %.0f-%.0f%s  |y|<2.4",centmin,centmax,texper.Data()));
  tex->SetNDC();
  tex->SetTextColor(1);
  tex->SetTextFont(42);
  tex->SetTextSize(0.045);
  tex->SetLineWidth(2);
  tex->Draw();

  tex = new TLatex(0.46,0.730,Form("Deviation from unity at no cut=%.2f%s",100*(f->GetParameter(0))*(fmin-cutvarref[varnum]),texper.Data()));
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.030);
  tex->SetLineWidth(2);
  //tex->Draw();

  c0->SaveAs(Form("plotCutVar/cutvariation_%s_%s_pt%.0f-%.0f.png",cutvarname[varnum].Data(),collsyst.Data(),_ptBins[0],_ptBins[1]));
  c0->SaveAs(Form("plotCutVar/cutvariation_%s_%s_pt%.0f-%.0f.pdf",cutvarname[varnum].Data(),collsyst.Data(),_ptBins[0],_ptBins[1]));
  
  return;
}

void clean0(TH1D* h)
{
  for (int i=1;i<=h->GetNbinsX();i++)
    {
      if(h->GetBinContent(i)==0) h->SetBinError(i,1);
    }
}

  TF1 *fit(TTree *nt, TTree *ntMC, Double_t ptmin, Double_t ptmax, int isMC, bool isPbPb, TF1* &total, Float_t centmin, Float_t centmax, TString npfit)
{
  static Int_t count=0;
  count++;
  TCanvas* c= new TCanvas(Form("c%d",count),"",600,600);
  TH1D* h = new TH1D(Form("h-%d",count),"",nbinsmasshisto,minhisto,maxhisto);
  TH1D* hMCSignal = new TH1D(Form("hMCSignal-%d",count),"",nbinsmasshisto,minhisto,maxhisto);
  
  TString iNP = npfit;
  TF1* f = new TF1(Form("f%d",count),"[0]*([7]*TMath::Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[7])*TMath::Gaus(x,[1],[8])/(sqrt(2*3.14159)*[8]))+[3]+[4]*x+[5]*("+iNP+")");
  f->SetNpx(5000);
  f->SetLineWidth(5);
  
  if(isMC==1) ntMC->Project(Form("h-%d",count),"Bmass",Form("(%s&&Bpt>%f&&Bpt<%f)*(1/%s)",seldata.Data(),ptmin,ptmax,weightdata.Data()));
  else nt->Project(Form("h-%d",count),"Bmass",Form("(%s&&Bpt>%f&&Bpt<%f)*(1/%s)",seldata.Data(),ptmin,ptmax,weightdata.Data()));
  ntMC->Project(Form("hMCSignal-%d",count),"Bmass",Form("(%s&&Bgen==23333&&Bpt>%f&&Bpt<%f)*(1/%s)",selmc.Data(),ptmin,ptmax,weightdata.Data()));

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
  f->FixParameter(3,0);
  f->FixParameter(4,0);
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
  
  //Fix the signal fit and do the background fits
  
  f->FixParameter(1,f->GetParameter(1));
  f->FixParameter(2,f->GetParameter(2));
  f->FixParameter(7,f->GetParameter(7));
  f->FixParameter(8,f->GetParameter(8));
  
  f->ReleaseParameter(3);
  f->ReleaseParameter(4);
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
  //background->SetLineStyle(2);//PAS
  background->SetLineStyle(7);//paper
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
  Bkpi->SetLineWidth(9);
  
  TF1 *mass = new TF1(Form("fmass%d",count),"[0]*([3]*TMath::Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[3])*TMath::Gaus(x,[1],[4])/(sqrt(2*3.14159)*[4]))");
  mass->SetParameters(f->GetParameter(0),f->GetParameter(1),f->GetParameter(2),f->GetParameter(7),f->GetParameter(8));
  mass->SetParError(0,f->GetParError(0));
  mass->SetParError(1,f->GetParError(1));
  mass->SetParError(2,f->GetParError(2));
  mass->SetParError(7,f->GetParError(7));
  mass->SetParError(8,f->GetParError(8));
  //mass->SetLineColor(2);//PAS
  //mass->SetFillColor(2);//PAS
  mass->SetFillColor(kOrange-3);//paper
  mass->SetLineColor(kOrange-3);//paper
  //mass->SetFillStyle(3004);//PAS
  mass->SetFillStyle(3002);//paper
  mass->SetLineWidth(9);
  //mass->SetLineStyle(2);//PAS
  mass->SetLineStyle(7);//paper
  
  //h->SetXTitle("m_{#mu#muK} (GeV/c^{2})");
  h->SetXTitle("m_{B} (GeV/c^{2})");
  
  double yunit = 1000./nbinsmasshisto;
  h->SetYTitle(Form("Events / (%.0f MeV/c^{2})",yunit));//basically, the unit is equivalent to binwidthmass in MeV
  
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  h->SetAxisRange(0,h->GetMaximum()*1.4*1.2,"Y");
  h->GetXaxis()->SetTitleOffset(0.9);
  h->GetYaxis()->SetTitleOffset(1.3);
  //h->GetXaxis()->SetLabelOffset(0.007);
  //h->GetYaxis()->SetLabelOffset(0.007);
  h->GetXaxis()->SetTitleSize(0.07);
  h->GetYaxis()->SetTitleSize(0.07);
  h->GetXaxis()->SetTitleFont(42);
  h->GetYaxis()->SetTitleFont(42);
  h->GetXaxis()->SetLabelFont(42);
  h->GetYaxis()->SetLabelFont(42);
  h->GetXaxis()->SetLabelSize(0.06);
  h->GetYaxis()->SetLabelSize(0.06);
  h->SetMarkerSize(1.55);
  h->SetMarkerStyle(20);
  h->SetLineColor(1);
  h->SetLineWidth(5);
  h->SetStats(0);
  h->GetXaxis()->SetNdivisions(-50205);
  h->Draw("e");
  Bkpi->SetRange(5.00,5.60);
  Bkpi->Draw("same");
  background->Draw("same");   
  mass->SetRange(5.16,5.40);
  mass->Draw("same");
  f->Draw("same");
  c->RedrawAxis();
  
  yield = mass->Integral(minhisto,maxhisto)/binwidthmass;
  yieldErr = mass->Integral(minhisto,maxhisto)/binwidthmass*mass->GetParError(0)/mass->GetParameter(0);
  printf("p_t bin %.0f-%.0f     yield: %f     yieldErr: %f\n", ptmin, ptmax, yield, yieldErr);
  
  Double_t Signal = mass->Integral(5.19932,5.35932)/binwidthmass; //B+ mass_pdg=5.27932GeV, signal region = pm 0.08GeV
  Double_t Bkg_comb = background->Integral(5.19932,5.35932)/binwidthmass;
  Double_t Bkg_nonprompt = Bkpi->Integral(5.19932,5.35932)/binwidthmass;
  
  //printf("p_t bin %.0f-%.0f S: %f B(comb): %f B(comb+np): %f\n", ptmin, ptmax, Signal, Bkg_comb, Bkg_comb+Bkg_nonprompt);
  //printf("p_t bin %.0f-%.0f sig(comb): %f sig(comb+np): %f\n", ptmin, ptmax, Signal/sqrt(Signal+Bkg_comb), Signal/sqrt(Signal+Bkg_comb+Bkg_nonprompt));

  TLegend *leg = new TLegend(0.55,0.55,0.875,0.775,NULL,"brNDC");//paper
  leg->SetBorderSize(0);
  leg->SetTextSize(0.045);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->AddEntry(h,"Data","pl");
  leg->AddEntry(f,"Fit","l");
  leg->AddEntry(mass,"Signal","f");
  leg->AddEntry(background,"Combinatorial","l");
  leg->AddEntry(Bkpi,"B #rightarrow J/#psi X","f");
  leg->Draw("same");
  
  TLatex* texYield = new TLatex(0.55,0.51,Form("Yield:%.2f", yield));
  texYield->SetNDC();
  texYield->SetTextAlign(12);
  texYield->SetTextSize(0.035);
  texYield->SetTextFont(42);
  texYield->Draw();
  
  TLatex* texChi = new TLatex(0.55,0.475,Form("#chi^{2}/nDOF:%.2f/%d=%.2f", f->GetChisquare(), f->GetNDF(), f->GetChisquare()/f->GetNDF()));
  texChi->SetNDC();
  texChi->SetTextAlign(12);
  texChi->SetTextSize(0.035);
  texChi->SetTextFont(42);
  texChi->Draw();
  printf("NDF: %d, chi2: %f, prob: %f\n", f->GetNDF(), f->GetChisquare(), f->GetProb());
  
  TLatex* texcms = new TLatex(0.22,0.87,"CMS");
  texcms->SetNDC();
  texcms->SetTextAlign(13);
  texcms->SetTextFont(62);
  texcms->SetTextSize(0.08);
  texcms->SetLineWidth(2);
  texcms->Draw();
  
  TLatex* texB = new TLatex(0.22,0.73,"B^{+}+B^{-}");
  texB->SetNDC();
  texB->SetTextFont(42);
  texB->SetTextSize(0.07);
  texB->SetLineWidth(2);
  texB->Draw();
  
  TLatex* texCol;
  if(collisionsystem=="pp"||collisionsystem=="PP"||collisionsystem=="ppInc"||collisionsystem=="PbPbInc") texCol= new TLatex(0.945,0.94, Form("28.0 pb^{-1} (%s 5.02 TeV)","pp"));
  else texCol= new TLatex(0.94,0.93, Form("1.5 nb^{-1} (%s 5.02 TeV)","PbPb"));
  texCol->SetNDC();
  texCol->SetTextAlign(32);
  texCol->SetTextSize(0.06);
  texCol->SetTextFont(42);
  texCol->Draw();
  
  TLatex* tex;
  tex = new TLatex(0.49,0.845,Form("%.0f<p_{T}<%.0f GeV/c",ptmin,ptmax));
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.045);
  tex->SetLineWidth(2);
  tex->Draw();
  
  TString texper="%";
  tex = new TLatex(0.488,0.80,Form("Cent. %.0f-%.0f%s  |y|<2.4",centmin,centmax,texper.Data()));
  tex->SetNDC();
  tex->SetTextColor(1);
  tex->SetTextFont(42);
  tex->SetTextSize(0.045);
  tex->SetLineWidth(2);
  tex->Draw();
  
  total=f;
  
  TF1* t = (TF1*)h->GetFunction(Form("f%d",count))->Clone();
  h->GetFunction(Form("f%d",count))->Delete();
  t->Draw("same");
  h->Draw("e same");
  //h->Write();
  //hMCSignal->Write();
  
  TString _postfix = "";
  if(weightdata!="1") _postfix = "_EFFCOR";
  if(isPbPb && isMC==0) 
    {
      c->SaveAs(Form("plotCutVar/Fits/%s/data_PbPb_pt%.0f-%.0f_%s_%d-thcut.png",cutvarname[VARNUM].Data(),ptmin,ptmax,cutvarname[VARNUM].Data(),stage));
      c->SaveAs(Form("plotCutVar/Fits/%s/data_PbPb_pt%.0f-%.0f_%s_%d-thcut.pdf",cutvarname[VARNUM].Data(),ptmin,ptmax,cutvarname[VARNUM].Data(),stage));
    }
  else if(isPbPb && isMC==1) 
    {
      c->SaveAs(Form("plotCutVar/Fits/%s/mc_PbPb_pt%.0f-%.0f_%s_%d-thcut.png",cutvarname[VARNUM].Data(),ptmin,ptmax,cutvarname[VARNUM].Data(),stage));
      c->SaveAs(Form("plotCutVar/Fits/%s/mc_PbPb_pt%.0f-%.0f_%s_%d-thcut.pdf",cutvarname[VARNUM].Data(),ptmin,ptmax,cutvarname[VARNUM].Data(),stage));
    }
  else if(!isPbPb && isMC==0) 
    {
      c->SaveAs(Form("plotCutVar/Fits/%s/data_pp_pt%.0f-%.0f_%s_%d-thcut.png",cutvarname[VARNUM].Data(),ptmin,ptmax,cutvarname[VARNUM].Data(),stage));
      c->SaveAs(Form("plotCutVar/Fits/%s/data_pp_pt%.0f-%.0f_%s_%d-thcut.pdf",cutvarname[VARNUM].Data(),ptmin,ptmax,cutvarname[VARNUM].Data(),stage));
    }
  else 
    {
      c->SaveAs(Form("plotCutVar/Fits/%s/mc_pp_pt%.0f-%.0f_%s_%d-thcut.png",cutvarname[VARNUM].Data(),ptmin,ptmax,cutvarname[VARNUM].Data(),stage));
      c->SaveAs(Form("plotCutVar/Fits/%s/mc_pp_pt%.0f-%.0f_%s_%d-thcut.pdf",cutvarname[VARNUM].Data(),ptmin,ptmax,cutvarname[VARNUM].Data(),stage));
    }  
  return mass;
}

int main(int argc, char *argv[])
{
  if(argc==17)
    {
      fitB(atoi(argv[1]), argv[2], argv[3], argv[4], argv[5], argv[6], atoi(argv[7]), atof(argv[8]), atoi(argv[9]), argv[10], argv[11], argv[12], atoi(argv[13]), atof(argv[14]), atof(argv[15]), atoi(argv[16]));
      return 0;
    }
  else if(argc==15)
    {
      fitB(atoi(argv[1]), argv[2], argv[3], argv[4], argv[5], argv[6], atoi(argv[7]), atof(argv[8]), atoi(argv[9]), argv[10], argv[11], argv[12], atoi(argv[13]), atoi(argv[14]));
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
