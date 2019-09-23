#include "uti.h"
#include "parameters_cutvar_cent.h"
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

  cutvarname = "BDT";
  cutvar[0] = "BDT_5_7";
  cutvar[1] = "BDT_7_10";
  //cutvar[2] = "BDT_10_15";
  //cutvar[3] = "BDT_15_20";
  //cutvar[4] = "BDT_20_30";
  //cutvar[5] = "BDT_30_50";

  /*
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
  */

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
      weightgen="pthatweight*(3.00448277-0.35865276*Gpt+0.01997413*Gpt*Gpt-0.00042585*Gpt*Gpt*Gpt+0.00000315*Gpt*Gpt*Gpt*Gpt)";
      weight="pthatweight*Ncoll*(TMath::Gaus(PVz,0.427450,4.873825)/(sqrt(2*3.14159)*4.873825))/(TMath::Gaus(PVz,0.909938,4.970989)/(sqrt(2*3\
.14159)*4.970989))*(3.00448277-0.35865276*Bgenpt+0.01997413*Bgenpt*Bgenpt-0.00042585*Bgenpt*Bgenpt*Bgenpt+0.00000315*Bgenpt*Bgenpt*Bgenpt*Bgenpt)";
      //weightgen="pthatweight*((3.506006+0.963473*Gpt+-0.258731*Gpt*Gpt)*TMath::Exp(-0.386065*Gpt)+1.139897)";
      //weight="pthatweight*Ncoll*(TMath::Gaus(PVz,0.427450,4.873825)/(sqrt(2*3.14159)*4.873825))/(TMath::Gaus(PVz,0.909938,4.970989)/(sqrt(2*3.14159)*4.970989))*((3.506006+0.963473*Bgenpt+-0.258731*Bgenpt*Bgenpt)*TMath::Exp(-0.386065*Bgenpt)+1.139897)";
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

  TString cutrefforothers1 = "Bpt>5.0 && Bpt<7.0 && (BsvpvDistance/BsvpvDisErr)>12.0 && cos(Bdtheta)>0.95";
  TString cutrefforothers2 = "Bpt>7.0 && Bpt<10.0";
  TString cutrefforothers3 = "Bpt>10.0 && Bpt<15.0";
  TString cutrefforothers4 = "Bpt>15.0 && Bpt<20.0";
  TString cutrefforothers5 = "Bpt>20.0 && Bpt<30.0";
  TString cutrefforothers6 = "Bpt>30.0 && Bpt<50.0";

  TString cutforothers1 = "";
  TString cutforothers2 = "";
  TString cutforothers3 = "";
  TString cutforothers4 = "";
  TString cutforothers5 = "";
  TString cutforothers6 = "";

  cutforothers1 = Form("%s&&(%s>%f)",cutrefforothers1.Data(),cutvar[0].Data(),cutvarref1[0]);
  cutforothers2 = Form("%s&&(%s>%f)",cutrefforothers2.Data(),cutvar[1].Data(),cutvarref2[0]);
  //cutforothers3 = Form("%s&&(%s>%f)",cutrefforothers3.Data(),cutvar[2].Data(),cutvarref3[0]);
  //cutforothers4 = Form("%s&&(%s>%f)",cutrefforothers4.Data(),cutvar[3].Data(),cutvarref4[0]);
  //cutforothers5 = Form("%s&&(%s>%f)",cutrefforothers5.Data(),cutvar[4].Data(),cutvarref5[0]);
  //cutforothers6 = Form("%s&&(%s>%f)",cutrefforothers6.Data(),cutvar[5].Data(),cutvarref6[0]);

  //TString cutref = Form("%s&&(%s||%s||%s||%s||%s||%s)",cut.Data(),cutforothers1.Data(),cutforothers2.Data(),cutforothers3.Data(),cutforothers4.Data(),cutforothers5.Data(),cutforothers6.Data());
  TString cutref = Form("%s&&(%s||%s)",cut.Data(),cutforothers1.Data(),cutforothers2.Data());

  if(!isPbPb)
    {
      seldata = Form("%s&&%s",trgselection.Data(),cutref.Data());
      selmc = Form("%s&&%s",trgselection.Data(),cutref.Data());
    }
  else
    {
      seldata = Form("%s&&%s&&hiBin>=%f&&hiBin<=%f",trgselection.Data(),cutref.Data(),hiBinMin,hiBinMax);
      selmc = Form("%s&&%s&&hiBin>=%f&&hiBin<=%f",trgselection.Data(),cutref.Data(),hiBinMin,hiBinMax);
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

  Ratio_ref = yieldData_ref/yieldMC_ref;
  RatioErr_ref = Ratio_ref*sqrt((yieldErrData_ref/yieldData_ref)*(yieldErrData_ref/yieldData_ref)+(yieldErrMC_ref/yieldMC_ref)*(yieldErrMC_ref/yieldMC_ref));

  TH1D* cutvarhis = new TH1D(Form("cutvarhis_%s",cutvarname.Data()),"",ncutvar,-nleft*cutspacing,nright*cutspacing);
  cutvarhis->GetXaxis()->SetTitle(Form("%s variation",cutvarname.Data()));
  cutvarhis->GetYaxis()->SetTitle("Double Ratio");
  cutvarhis->GetXaxis()->CenterTitle();
  //cutvarhis->SetTitle(Form("Cut Variation of %s",cutvarname.Data()));
  cutvarhis->SetMarkerStyle(20);
  cutvarhis->SetMarkerSize(1);
  cutvarhis->GetYaxis()->SetRangeUser(0.,2.);  

  TString cutvar = "";

  double maxdeviation = 0.0;

  for(int i=0;i<ncutvar;i++)
   {
     stage = i;

     cutforothers1 = Form("%s&&(BDT_5_7>%f)",cutrefforothers1.Data(),cutvarref1[0]+(-nleft+i)*cutspacing);
     cutforothers2 = Form("%s&&(BDT_7_10>%f)",cutrefforothers2.Data(),cutvarref2[0]+(-nleft+i)*cutspacing);
     //cutforothers3 = Form("%s&&(BDT_10_15>%f)",cutrefforothers3.Data(),cutvarref3[0]+(-nleft+i)*cutspacing);
     //cutforothers4 = Form("%s&&(BDT_15_20>%f)",cutrefforothers4.Data(),cutvarref4[0]+(-nleft+i)*cutspacing);
     //cutforothers5 = Form("%s&&(BDT_20_30>%f)",cutrefforothers5.Data(),cutvarref5[0]+(-nleft+i)*cutspacing);
     //cutforothers6 = Form("%s&&(BDT_30_50>%f)",cutrefforothers6.Data(),cutvarref6[0]+(-nleft+i)*cutspacing);

     //cutvar = Form("%s&&(%s||%s||%s||%s||%s||%s)",cut.Data(),cutforothers1.Data(),cutforothers2.Data(),cutforothers3.Data(),cutforothers4.Data(),cutforothers5.Data(),cutforothers6.Data());
     cutvar = Form("%s&&(%s||%s)",cut.Data(),cutforothers1.Data(),cutforothers2.Data());

     if(!isPbPb)
       {
	 seldata = Form("%s&&%s",trgselection.Data(),cutvar.Data());
	 selmc = Form("%s&&%s",trgselection.Data(),cutvar.Data());
       }
     else
       {
	 seldata = Form("%s&&%s&&hiBin>=%f&&hiBin<=%f",trgselection.Data(),cutvar.Data(),hiBinMin,hiBinMax);
	 selmc = Form("%s&&%s&&hiBin>=%f&&hiBin<=%f",trgselection.Data(),cutvar.Data(),hiBinMin,hiBinMax);
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

     Ratio_var[i] = yieldData[i]/yieldMC[i];
     RatioErr_var[i] = Ratio_var[i]*sqrt((yieldErrData[i]/yieldData[i])*(yieldErrData[i]/yieldData[i])+(yieldErrMC[i]/yieldMC[i])*(yieldErrMC[i]/yieldMC[i]));

     DoubleRatio[i] = Ratio_var[i]/Ratio_ref;
     DoubleRatioErr[i] = sqrt(TMath::Abs(RatioErr_var[i]*RatioErr_var[i]-RatioErr_ref*RatioErr_ref))/Ratio_ref;

     std::cout<<""<<std::endl;
     std::cout<<"Double Ratio at "<<i<<"-th cut = "<<DoubleRatio[i]<<" #pm "<<DoubleRatioErr[i]<<std::endl;
     
     cutvarhis->SetBinContent(i+1,DoubleRatio[i]);
     cutvarhis->SetBinError(i+1,DoubleRatioErr[i]);

     //if(TMath::Abs(DoubleRatio[i]-1.0)>maxdeviation && i<=nleft) maxdeviation = TMath::Abs(DoubleRatio[i]-1.0);
   }

  double  DoubleRatio_nocut = DoubleRatio[3];
  maxdeviation = TMath::Abs(DoubleRatio_nocut-1.0);
  std::cout<<"Maximum deviation from unity = "<<maxdeviation*100<<"%"<<std::endl;

  /*
  double fmax, fmin;
  fmax = cutvarmax[varnum];
  fmin = -1.0;
  //if(varnum==1) fmin = -1.0;
  //else fmin = 0.0;
  
  TF1* f = new TF1("f",Form("1.0+[0]*(x-%f)",cutvarref[varnum]),fmin,cutvarref[varnum]);
  f->SetParLimits(0,-100,100);
  cutvarhis->Fit(f,"R");
  double intercept = 1.0-(f->GetParameter(0))*(cutvarref[varnum]);
  double slope = (f->GetParameter(0));

  std::cout<<"Linear Fit Function = "<<intercept<<"+"<<slope<<"*"<<cutvarname[varnum]<<std::endl;
  */

  TFile* outputroot = new TFile(Form("plotCutVar/Cent/cutvariation_%s_%s_pt%.0f-%.0f_cent%.0f-%.0f.root",cutvarname.Data(),collsyst.Data(),_ptBins[0],_ptBins[1],centmin,centmax),"recreate");
  outputroot->cd();
  cutvarhis->Write();
  outputroot->Close();
    
  TCanvas* c0 = new TCanvas("","",600,600);
  c0->cd();
  cutvarhis->Draw("ep");
  //f->Draw("same");
  //c0->RedrawAxis();

  TLine* vline_ref = new TLine(0.0,0.0,0.0,1.0);
  vline_ref->SetLineWidth(1);
  vline_ref->SetLineStyle(2);
  vline_ref->SetLineColor(kGreen);
  vline_ref->Draw("same");

  TLine* hline_ref = new TLine(-nleft*cutspacing,1.0,0.0,1.0);
  hline_ref->SetLineWidth(1);
  hline_ref->SetLineStyle(2);
  hline_ref->SetLineColor(kGreen);
  hline_ref->Draw("same");

  TLine* vline_nocut = new TLine((-nleft+3)*cutspacing,0.0,(-nleft+3)*cutspacing,DoubleRatio_nocut);
  vline_nocut->SetLineWidth(1);
  vline_nocut->SetLineStyle(2);
  vline_nocut->SetLineColor(kRed);
  vline_nocut->Draw("same");

  TLine* hline_nocut = new TLine(-nleft*cutspacing,DoubleRatio_nocut,(-nleft+3)*cutspacing,DoubleRatio_nocut);
  hline_nocut->SetLineWidth(1);
  hline_nocut->SetLineStyle(2);
  hline_nocut->SetLineColor(kRed);
  hline_nocut->Draw("same");

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

  
  //tex = new TLatex(0.46,0.730,Form("Maximum deviation (looser variation)=%.2f%s",maxdeviation*100,texper.Data()));
  tex = new TLatex(0.46,0.730,Form("Deviation at effective no cut=%.2f%s",maxdeviation*100,texper.Data()));
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.030);
  tex->SetLineWidth(2);
  tex->Draw();
  

  c0->SaveAs(Form("plotCutVar/Cent/cutvariation_%s_%s_pt%.0f-%.0f_cent%.0f-%.0f.png",cutvarname.Data(),collsyst.Data(),_ptBins[0],_ptBins[1],centmin,centmax));
  c0->SaveAs(Form("plotCutVar/Cent/cutvariation_%s_%s_pt%.0f-%.0f_cent%.0f-%.0f.pdf",cutvarname.Data(),collsyst.Data(),_ptBins[0],_ptBins[1],centmin,centmax));
  
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
  f->SetParLimits(1,5.25,5.30);
  
  //Do the signal fit first
  
  f->SetParameter(0,setparam0);
  f->SetParameter(1,setparam1);
  f->SetParameter(2,setparam2);
  f->SetParameter(8,setparam3);
  //f->FixParameter(1,fixparam1);

  f->FixParameter(3,0);
  f->FixParameter(4,0);
  f->FixParameter(5,0);

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
  
  f->ReleaseParameter(3);
  f->ReleaseParameter(4);
  f->ReleaseParameter(5);
  f->SetParLimits(5,0,1e4);
  
  printf("Fixed para.:\n");
  printf("%f, %f, %f\n", f->GetParameter(2), f->GetParameter(7), f->GetParameter(8));
  h->Fit(Form("f%d",count),"q","",minhisto,maxhisto);
  h->Fit(Form("f%d",count),"q","",minhisto,maxhisto);
  //f->ReleaseParameter(1);
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
  Bkpi->SetRange(minhisto,maxhisto);
  Bkpi->Draw("same");
  background->Draw("same");
  mass->SetRange(minhisto,maxhisto);
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
      c->SaveAs(Form("plotCutVar/Fits/%s/Cent/data_PbPb_pt%.0f-%.0f_cent%.0f-%.0f_%s_%d-thcut.png",cutvarname.Data(),ptmin,ptmax,centmin,centmax,cutvarname.Data(),stage));
      c->SaveAs(Form("plotCutVar/Fits/%s/Cent/data_PbPb_pt%.0f-%.0f_cent%.0f-%.0f_%s_%d-thcut.pdf",cutvarname.Data(),ptmin,ptmax,centmin,centmax,cutvarname.Data(),stage));
    }
  else if(isPbPb && isMC==1) 
    {
      c->SaveAs(Form("plotCutVar/Fits/%s/Cent/mc_PbPb_pt%.0f-%.0f_cent%.0f-%.0f_%s_%d-thcut.png",cutvarname.Data(),ptmin,ptmax,centmin,centmax,cutvarname.Data(),stage));
      c->SaveAs(Form("plotCutVar/Fits/%s/Cent/mc_PbPb_pt%.0f-%.0f_cent%.0f-%.0f_%s_%d-thcut.pdf",cutvarname.Data(),ptmin,ptmax,centmin,centmax,cutvarname.Data(),stage));
    }
  else if(!isPbPb && isMC==0) 
    {
      c->SaveAs(Form("plotCutVar/Fits/%s/Cent/data_pp_pt%.0f-%.0f_cent%.0f-%.0f_%s_%d-thcut.png",cutvarname.Data(),ptmin,ptmax,centmin,centmax,cutvarname.Data(),stage));
      c->SaveAs(Form("plotCutVar/Fits/%s/Cent/data_pp_pt%.0f-%.0f_cent%.0f-%.0f_%s_%d-thcut.pdf",cutvarname.Data(),ptmin,ptmax,centmin,centmax,cutvarname.Data(),stage));
    }
  else 
    {
      c->SaveAs(Form("plotCutVar/Fits/%s/Cent/mc_pp_pt%.0f-%.0f_cent%.0f-%.0f_%s_%d-thcut.png",cutvarname.Data(),ptmin,ptmax,centmin,centmax,cutvarname.Data(),stage));
      c->SaveAs(Form("plotCutVar/Fits/%s/Cent/mc_pp_pt%.0f-%.0f_cent%.0f-%.0f_%s_%d-thcut.pdf",cutvarname.Data(),ptmin,ptmax,centmin,centmax,cutvarname.Data(),stage));
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
