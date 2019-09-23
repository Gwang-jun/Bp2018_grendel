#include "uti.h"
#include "parameters.h"

Double_t setparam0=100.;
Double_t setparam1=1.865;
Double_t setparam2=0.03;
Double_t setparam10=0.005;
Double_t setparam8=0.1;
Double_t setparam9=0.1;
Double_t fixparam1=1.865;
Double_t minhisto=1.7;
Double_t maxhisto=2.0;
Double_t nbinsmasshisto=60;
Double_t binwidthmass=(maxhisto-minhisto)/nbinsmasshisto;

Float_t hiBinMin,hiBinMax,centMin,centMax;

TString selmcgenref, selmcgenacceptanceref, cut_recoonlyref, cutref;

int _nBins = nBinsCent;
double *_ptBins = ptBinsCent;

void MCefficiencyCent(int isPbPb=0,TString inputmc="", TString selmcgen="",TString selmcgenacceptance="", TString cut_recoonly="", TString cut="",TString label="",TString outputfile="", int useweight=1,Float_t centmin=0., Float_t centmax=100.)
{ 
  hiBinMin = centmin*2;
  hiBinMax = centmax*2;
  centMin = centmin;
  centMax = centmax;

  selmcgen = selmcgen+Form("&&hiBin>=%f&&hiBin<=%f&&Gpt>%f&&Gpt<%f",hiBinMin,hiBinMax,ptBinsInc[0],ptBinsInc[1]);
  selmcgenacceptance=selmcgenacceptance+Form("&&hiBin>=%f&&hiBin<=%f&&Gpt>%f&&Gpt<%f",hiBinMin,hiBinMax,ptBinsInc[0],ptBinsInc[1]);
  cut_recoonly=cut_recoonly+Form("&&hiBin>=%f&&hiBin<=%f&&Bpt>%f&&Bpt<%f",hiBinMin,hiBinMax,ptBinsInc[0],ptBinsInc[1]);
  cut=cut+Form("&&hiBin>=%f&&hiBin<=%f&&Bpt>%f&&Bpt<%f",hiBinMin,hiBinMax,ptBinsInc[0],ptBinsInc[1]);

  //selmcgenref = selmcgen;
  //selmcgenacceptanceref = selmcgenacceptance;
  //cut_recoonlyref = cut_recoonly;
  //cutref = cut;

  //std::cout<<"selmcgen="<<selmcgen<<std::endl;
  //std::cout<<"selmcgenacceptance="<<selmcgenacceptance<<std::endl;
  //std::cout<<"cut_recoonly"<<cut_recoonly<<std::endl;
  //std::cout<<"cut="<<cut<<std::endl;

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);
  gStyle->SetMarkerStyle(20);

  //For 2018 PbPb MC
  TFile* infMC = new TFile(inputmc.Data());
  TTree* ntMC = (TTree*)infMC->Get("Bfinder/ntKp");
  ntMC->AddFriend("hltanalysis/HltTree");
  ntMC->AddFriend("hiEvtAnalyzer/HiTree");
  ntMC->AddFriend("Bfinder/ntGen");
  ntMC->AddFriend("skimanalysis/HltTree");
  ntMC->AddFriend("BDT");
  TTree* ntGen = (TTree*)infMC->Get("Bfinder/ntGen");
  ntGen->AddFriend("hltanalysis/HltTree");
  ntGen->AddFriend("hiEvtAnalyzer/HiTree");
  ntGen->AddFriend("Bfinder/ntKp");
  ntGen->AddFriend("skimanalysis/HltTree");
  ntGen->AddFriend("BDT");

  // optimal weigths
  TCut weighpthat = "1";
  TCut weightGpt = "1";
  TCut weightBgenpt = "1";
  TCut weightHiBin = "1";
  TCut weightPVz = "1";
  if(useweight==0) { //pp
    weighpthat = "pthatweight";
    weightPVz = "1.055564*TMath::Exp(-0.001720*(PVz+2.375584)*(PVz+2.375584))";
    weightGpt = "0.599212+-0.020703*Gpt+0.003143*Gpt*Gpt+-0.000034*Gpt*Gpt*Gpt";
    weightBgenpt = "0.599212+-0.020703*Bgenpt+0.003143*Bgenpt*Bgenpt+-0.000034*Bgenpt*Bgenpt*Bgenpt";
  }

  if(useweight==1) { //PbPb
    weighpthat = "pthatweight";
    weightHiBin = "Ncoll";
    weightPVz = "(TMath::Gaus(PVz,0.427450,4.873825)/(sqrt(2*3.14159)*4.873825))/(TMath::Gaus(PVz,0.909938,4.970989)/(sqrt(2*3.14159)*4.970989))";
    //weightGpt = "(2.907795+-0.436572*Gpt+0.006372*Gpt*Gpt)*TMath::Exp(-0.157563*Gpt)+1.01308";
    //weightBgenpt = "(2.907795+-0.436572*Bgenpt+0.006372*Bgenpt*Bgenpt)*TMath::Exp(-0.157563*Bgenpt)+1.01308";
    //weightGpt = "(3.506006+0.963473*Gpt+-0.258731*Gpt*Gpt)*TMath::Exp(-0.386065*Gpt)+1.139897";
    //weightBgenpt = "(3.506006+0.963473*Bgenpt+-0.258731*Bgenpt*Bgenpt)*TMath::Exp(-0.386065*Bgenpt)+1.139897";
    //weightGpt = "(1.883180+-0.290677*Gpt+0.000225*Gpt*Gpt)*TMath::Exp(-0.161669*Gpt)+1.171923";
    //weightBgenpt = "(1.883180+-0.290677*Bgenpt+0.000225*Bgenpt*Bgenpt)*TMath::Exp(-0.161669*Bgenpt)+1.171923";
    weightGpt = "(3.00448277-0.35865276*Gpt+0.01997413*Gpt*Gpt-0.00042585*Gpt*Gpt*Gpt+0.00000315*Gpt*Gpt*Gpt*Gpt)";
    weightBgenpt = "(3.00448277-0.35865276*Bgenpt+0.01997413*Bgenpt*Bgenpt-0.00042585*Bgenpt*Bgenpt*Bgenpt+0.00000315*Bgenpt*Bgenpt*Bgenpt*Bgenpt)";
    //weightGpt = "1";
    //weightBgenpt = "1";
    //weightGpt = "(3.76547732-0.48262502*Gpt+0.02740408*Gpt*Gpt-0.00060885*Gpt*Gpt*Gpt+0.00000478*Gpt*Gpt*Gpt*Gpt)";
    //weightBgenpt = "(3.76547732-0.48262502*Bgenpt+0.02740408*Bgenpt*Bgenpt-0.00060885*Bgenpt*Bgenpt*Bgenpt+0.00000478*Bgenpt*Bgenpt*Bgenpt*Bgenpt)";
    //weightGpt = "(2.32807290-0.24845949*Gpt+0.01337001*Gpt*Gpt-0.00026320*Gpt*Gpt*Gpt+0.00000170*Gpt*Gpt*Gpt*Gpt)";
    //weightBgenpt = "(2.32807290-0.24845949*Bgenpt+0.01337001*Bgenpt*Bgenpt-0.00026320*Bgenpt*Bgenpt*Bgenpt+0.00000170*Bgenpt*Bgenpt*Bgenpt*Bgenpt)";
  }

  TH1D* hPtMC = new TH1D("hPtMC","",_nBins,_ptBins);
  TH1D* hPtMCrecoonly = new TH1D("hPtMCrecoonly","",_nBins,_ptBins);
  TH1D* hPtGen = new TH1D("hPtGen","",_nBins,_ptBins);
  TH1D* hPtGenAcc = new TH1D("hPtGenAcc","",_nBins,_ptBins);
  TH1D* hPtGenAccWeighted = new TH1D("hPtGenAccWeighted","",_nBins,_ptBins);
  TH1D* hPthat = new TH1D("hPthat","",100,0,500);
  TH1D* hPthatweight = new TH1D("hPthatweight","",100,0,500);

  ntMC->Project("hPtMC","hiBin",TCut(weighpthat)*TCut(weightBgenpt)*TCut(weightHiBin)*TCut(weightPVz)*(TCut(cut.Data())&&"(Bgen==23333)"));
  ntMC->Project("hPtMCrecoonly","hiBin",TCut(weighpthat)*TCut(weightBgenpt)*TCut(weightHiBin)*TCut(weightPVz)*(TCut(cut_recoonly.Data())&&"(Bgen==23333)"));
  ntGen->Project("hPtGen","hiBin",TCut(weighpthat)*TCut(weightGpt)*(TCut(selmcgen.Data())));
  ntGen->Project("hPtGenAcc","hiBin",TCut(weighpthat)*TCut(weightGpt)*(TCut(selmcgenacceptance.Data())));
  ntGen->Project("hPtGenAccWeighted","hiBin",TCut(weighpthat)*TCut(weightGpt)*TCut(weightHiBin)*TCut(weightPVz)*(TCut(selmcgenacceptance.Data())));

  ////// tag & probe scaling factor
  for(int i = 0; i < _nBins; i++){printf("%.2f, ", hPtMC->GetBinContent(i+1));}printf("\n");
  double sf_pp[2] = {1., 1.};
  //double sf_pbpb[2] = {1.0911, 1.1013};
  double sf_pbpb[1] = {1.0936};

  for(int i = 0; i < _nBins; i++){
    if(label == "pp"){
      hPtMC->SetBinContent(i+1, hPtMC->GetBinContent(i+1)*sf_pp[i]);
      hPtMCrecoonly->SetBinContent(i+1, hPtMCrecoonly->GetBinContent(i+1)*sf_pp[i]);
    }
    if(label == "PbPb"){
      hPtMC->SetBinContent(i+1, hPtMC->GetBinContent(i+1)*sf_pbpb[i]);
      hPtMCrecoonly->SetBinContent(i+1, hPtMCrecoonly->GetBinContent(i+1)*sf_pbpb[i]);
    }
  }
  //for(int i = 0; i < _nBins; i++){printf("%.2f, ", hPtMC->GetBinContent(i+1));}printf("\n");

  divideBinWidth(hPtMC);
  divideBinWidth(hPtMCrecoonly);
  divideBinWidth(hPtGen);
  divideBinWidth(hPtGenAcc);
  divideBinWidth(hPtGenAccWeighted);

  ntMC->Project("hPthat","pthat","1");
  ntMC->Project("hPthatweight","pthat","pthatweight");

  hPtMC->Sumw2();
  hPtGenAcc->Sumw2();
  hPtMCrecoonly->Sumw2();

  //hEffAcc = hPtGenAcc / hPtGen
  //hEffReco = hPtMCrecoonly / hPtGenAcc
  //hEffSelection = hPtMC / hPtMCrecoonly
  //hEff = hPtMC / hPtGen

  //Acceptance
  TH1D* hEffAcc = (TH1D*)hPtGenAcc->Clone("hEffAcc");
  hEffAcc->Sumw2();
  hEffAcc->Divide(hEffAcc,hPtGen,1,1,"b");
  //Selection
  TH1D* hEffSelection = (TH1D*)hPtMC->Clone("hEffSelection");
  hEffSelection->Sumw2();
  hEffSelection->Divide(hEffSelection,hPtGenAccWeighted,1,1,"b");
  //Acc * Eff (one shot)
  TH1D* hEffReco = (TH1D*)hPtMCrecoonly->Clone("hEffReco");
  hEffReco->Sumw2();
  hEffReco->Divide(hEffReco,hPtGen,1,1,"b");
  //Acc * Eff
  TH1D* hEff = (TH1D*)hEffSelection->Clone("hEff");
  hEff->Sumw2();
  //hEff->Divide(hPtMC,hPtGen,1,1,"");
  hEff->Multiply(hEff,hEffAcc,1,1);

  /*
  TFile* filenominalCent1 = new TFile("ptshape/BDT/MCstudiesPbPb_nominal_Cent0-30-90.root");
  TFile* fileplusCent1 = new TFile("ptshape/BDT/MCstudiesPbPb_plus_Cent0-30-90.root");
  TFile* fileminusCent1 = new TFile("ptshape/BDT/MCstudiesPbPb_minus_Cent0-30-90.root");
  TFile* filenominalCent2 = new TFile("ptshape/BDT/MCstudiesPbPb_nominal_Cent0-90.root");
  TFile* fileplusCent2 = new TFile("ptshape/BDT/MCstudiesPbPb_plus_Cent0-90.root");
  TFile* fileminusCent2 = new TFile("ptshape/BDT/MCstudiesPbPb_minus_Cent0-90.root");

  TH1D* hEffnominalCent1 = (TH1D*)filenominalCent1->Get("hEff");
  hEffnominalCent1->GetXaxis()->CenterTitle();
  hEffnominalCent1->GetYaxis()->CenterTitle();
  hEffnominalCent1->GetXaxis()->SetTitle("hiBin");
  hEffnominalCent1->GetYaxis()->SetTitle("#alpha x #epsilon");
  hEffnominalCent1->GetXaxis()->SetTitleOffset(0.9);
  hEffnominalCent1->GetYaxis()->SetTitleOffset(0.95);
  hEffnominalCent1->GetXaxis()->SetTitleSize(0.05);
  hEffnominalCent1->GetYaxis()->SetTitleSize(0.05);
  hEffnominalCent1->GetXaxis()->SetTitleFont(42);
  hEffnominalCent1->GetYaxis()->SetTitleFont(42);
  hEffnominalCent1->GetXaxis()->SetLabelFont(42);
  hEffnominalCent1->GetYaxis()->SetLabelFont(42);
  hEffnominalCent1->GetXaxis()->SetLabelSize(0.035);
  hEffnominalCent1->GetYaxis()->SetLabelSize(0.035);
  hEffnominalCent1->SetLineColor(kRed);
  TH1D* hEffplusCent1 = (TH1D*)fileplusCent1->Get("hEff");
  hEffplusCent1->GetXaxis()->CenterTitle();
  hEffplusCent1->GetYaxis()->CenterTitle();
  hEffplusCent1->GetXaxis()->SetTitle("hiBin");
  hEffplusCent1->GetYaxis()->SetTitle("#alpha x #epsilon");
  hEffplusCent1->GetXaxis()->SetTitleOffset(0.9);
  hEffplusCent1->GetYaxis()->SetTitleOffset(0.95);
  hEffplusCent1->GetXaxis()->SetTitleSize(0.05);
  hEffplusCent1->GetYaxis()->SetTitleSize(0.05);
  hEffplusCent1->GetXaxis()->SetTitleFont(42);
  hEffplusCent1->GetYaxis()->SetTitleFont(42);
  hEffplusCent1->GetXaxis()->SetLabelFont(42);
  hEffplusCent1->GetYaxis()->SetLabelFont(42);
  hEffplusCent1->GetXaxis()->SetLabelSize(0.035);
  hEffplusCent1->GetYaxis()->SetLabelSize(0.035);
  hEffplusCent1->SetLineColor(kBlue);
  TH1D* hEffminusCent1 = (TH1D*)fileminusCent1->Get("hEff");
  hEffminusCent1->GetXaxis()->CenterTitle();
  hEffminusCent1->GetYaxis()->CenterTitle();
  hEffminusCent1->GetXaxis()->SetTitle("hiBin");
  hEffminusCent1->GetYaxis()->SetTitle("#alpha x #epsilon");
  hEffminusCent1->GetXaxis()->SetTitleOffset(0.9);
  hEffminusCent1->GetYaxis()->SetTitleOffset(0.95);
  hEffminusCent1->GetXaxis()->SetTitleSize(0.05);
  hEffminusCent1->GetYaxis()->SetTitleSize(0.05);
  hEffminusCent1->GetXaxis()->SetTitleFont(42);
  hEffminusCent1->GetYaxis()->SetTitleFont(42);
  hEffminusCent1->GetXaxis()->SetLabelFont(42);
  hEffminusCent1->GetYaxis()->SetLabelFont(42);
  hEffminusCent1->GetXaxis()->SetLabelSize(0.035);
  hEffminusCent1->GetYaxis()->SetLabelSize(0.035);
  hEffminusCent1->SetLineColor(kBlue);
  TH1D* hEffnominalCent2 = (TH1D*)filenominalCent2->Get("hEff");
  hEffnominalCent2->GetXaxis()->CenterTitle();
  hEffnominalCent2->GetYaxis()->CenterTitle();
  hEffnominalCent2->GetXaxis()->SetTitle("hiBin");
  hEffnominalCent2->GetYaxis()->SetTitle("#alpha x #epsilon");
  hEffnominalCent2->GetXaxis()->SetTitleOffset(0.9);
  hEffnominalCent2->GetYaxis()->SetTitleOffset(0.95);
  hEffnominalCent2->GetXaxis()->SetTitleSize(0.05);
  hEffnominalCent2->GetYaxis()->SetTitleSize(0.05);
  hEffnominalCent2->GetXaxis()->SetTitleFont(42);
  hEffnominalCent2->GetYaxis()->SetTitleFont(42);
  hEffnominalCent2->GetXaxis()->SetLabelFont(42);
  hEffnominalCent2->GetYaxis()->SetLabelFont(42);
  hEffnominalCent2->GetXaxis()->SetLabelSize(0.035);
  hEffnominalCent2->GetYaxis()->SetLabelSize(0.035);
  hEffnominalCent2->SetLineColor(kRed);
  TH1D* hEffplusCent2 = (TH1D*)fileplusCent2->Get("hEff");
  hEffplusCent2->GetXaxis()->CenterTitle();
  hEffplusCent2->GetYaxis()->CenterTitle();
  hEffplusCent2->GetXaxis()->SetTitle("hiBin");
  hEffplusCent2->GetYaxis()->SetTitle("#alpha x #epsilon");
  hEffplusCent2->GetXaxis()->SetTitleOffset(0.9);
  hEffplusCent2->GetYaxis()->SetTitleOffset(0.95);
  hEffplusCent2->GetXaxis()->SetTitleSize(0.05);
  hEffplusCent2->GetYaxis()->SetTitleSize(0.05);
  hEffplusCent2->GetXaxis()->SetTitleFont(42);
  hEffplusCent2->GetYaxis()->SetTitleFont(42);
  hEffplusCent2->GetXaxis()->SetLabelFont(42);
  hEffplusCent2->GetYaxis()->SetLabelFont(42);
  hEffplusCent2->GetXaxis()->SetLabelSize(0.035);
  hEffplusCent2->GetYaxis()->SetLabelSize(0.035);
  hEffplusCent2->SetLineColor(kBlue);
  TH1D* hEffminusCent2 = (TH1D*)fileminusCent2->Get("hEff");
  hEffminusCent2->GetXaxis()->CenterTitle();
  hEffminusCent2->GetYaxis()->CenterTitle();
  hEffminusCent2->GetXaxis()->SetTitle("hiBin");
  hEffminusCent2->GetYaxis()->SetTitle("#alpha x #epsilon");
  hEffminusCent2->GetXaxis()->SetTitleOffset(0.9);
  hEffminusCent2->GetYaxis()->SetTitleOffset(0.95);
  hEffminusCent2->GetXaxis()->SetTitleSize(0.05);
  hEffminusCent2->GetYaxis()->SetTitleSize(0.05);
  hEffminusCent2->GetXaxis()->SetTitleFont(42);
  hEffminusCent2->GetYaxis()->SetTitleFont(42);
  hEffminusCent2->GetXaxis()->SetLabelFont(42);
  hEffminusCent2->GetYaxis()->SetLabelFont(42);
  hEffminusCent2->GetXaxis()->SetLabelSize(0.035);
  hEffminusCent2->GetYaxis()->SetLabelSize(0.035);
  hEffminusCent2->SetLineColor(kBlue);
  */
  
  /*
  TH1D* ptshapeplusCent1 = (TH1D*)hEffnominalCent1->Clone("ptshapeplusCent1");
  ptshapeplusCent1->Divide(hEffplusCent1);
  ptshapeplusCent1->SetLineColor(kBlue);

  TH1D* ptshapeminusCent1 = (TH1D*)hEffnominalCent1->Clone("ptshapeminusCent1");
  ptshapeminusCent1->Divide(hEffminusCent1);
  ptshapeminusCent1->SetLineColor(kRed);

  TCanvas* c100 = new TCanvas("","",600,600);
  c100->cd();
  ptshapeplusCent1->SetMaximum(1.3);
  ptshapeplusCent1->SetMinimum(0.7);  
  ptshapeplusCent1->GetXaxis()->SetTitle("hiBin");
  ptshapeplusCent1->GetYaxis()->SetTitle("#alpha x #epsilon ratio");
  ptshapeplusCent1->Draw();
  ptshapeminusCent1->Draw("same");

  TLegend *leg100 = new TLegend(0.45,0.70,0.75,0.80,NULL,"brNDC");                                                                               leg100->SetBorderSize(0);
  leg100->SetTextSize(0.04);
  leg100->SetTextFont(42);
  leg100->SetFillStyle(0);
  leg100->AddEntry(ptshapeplusCent1,"Data/MC +1 #sigma","l");
  leg100->AddEntry(ptshapeminusCent1,"Data/MC -1 #sigma","l");
  leg100->Draw("same");                                                                                                                        

  c100->SaveAs("ptshape/ptshape_Cent1_BDT_datadriven.png");
  c100->SaveAs("ptshape/ptshape_Cent1_BDT_datadriven.pdf");

  for(int j=0;j<2;j++)
    {
      printf("Cent bins %.0f-%.0f ptshape uncertainty(plus): %f (percent)\n",_ptBins[j],_ptBins[j+1],100.0*(ptshapeplusCent1->GetBinContent(j+1)-1.0));
      printf("Cent bins %.0f-%.0f ptshape uncertainty(minus): %f (percent)\n",_ptBins[j],_ptBins[j+1],100.0*(ptshapeminusCent1->GetBinContent(j+1)-1.0));
      printf("Cent bins %.0f-%.0f ptshape uncertainty(maximum): %f (percent)\n",_ptBins[j],_ptBins[j+1],TMath::Max(TMath::Abs(100.0*(ptshapeplusCent1->GetBinContent(j+1)-1.0)),TMath::Abs(100.0*(ptshapeminusCent1->GetBinContent(j+1)-1.0))));
    }
  */
  
  /*
  TH1D* ptshapeplusCent2 = (TH1D*)hEffnominalCent2->Clone("ptshapeplusCent2");
  ptshapeplusCent2->Divide(hEffplusCent2);
  ptshapeplusCent2->SetLineColor(kBlue);

  TH1D* ptshapeminusCent2 = (TH1D*)hEffnominalCent2->Clone("ptshapeminusCent2");
  ptshapeminusCent2->Divide(hEffminusCent2);
  ptshapeminusCent2->SetLineColor(kRed);

  TCanvas* c200 = new TCanvas("","",600,600);
  c200->cd();
  ptshapeplusCent2->SetMaximum(1.3);
  ptshapeplusCent2->SetMinimum(0.7);  
  ptshapeplusCent2->GetXaxis()->SetTitle("hiBin");
  ptshapeplusCent2->GetYaxis()->SetTitle("#alpha x #epsilon ratio");
  ptshapeplusCent2->Draw();
  ptshapeminusCent2->Draw("same");

  TLegend *leg200 = new TLegend(0.45,0.70,0.75,0.80,NULL,"brNDC");                                                                               leg200->SetBorderSize(0);
  leg200->SetTextSize(0.04);
  leg200->SetTextFont(42);
  leg200->SetFillStyle(0);
  leg200->AddEntry(ptshapeplusCent2,"Data/MC +1 #sigma","l");
  leg200->AddEntry(ptshapeminusCent2,"Data/MC -1 #sigma","l");
  leg200->Draw("same");                                                                                                                        

  c200->SaveAs("ptshape/ptshape_Cent2_BDT_datadriven.png");
  c200->SaveAs("ptshape/ptshape_Cent2_BDT_datadriven.pdf");

  for(int j=0;j<1;j++)
    {
      printf("Cent bins %.0f-%.0f ptshape uncertainty(plus): %f (percent)\n",_ptBins[j],_ptBins[j+1],100.0*(ptshapeplusCent2->GetBinContent(j+1)-1.0));
      printf("Cent bins %.0f-%.0f ptshape uncertainty(minus): %f (percent)\n",_ptBins[j],_ptBins[j+1],100.0*(ptshapeminusCent2->GetBinContent(j+1)-1.0));
      printf("Cent bins %.0f-%.0f ptshape uncertainty(maximum): %f (percent)\n",_ptBins[j],_ptBins[j+1],TMath::Max(TMath::Abs(100.0*(ptshapeplusCent2->GetBinContent(j+1)-1.0)),TMath::Abs(100.0*(ptshapeminusCent2->GetBinContent(j+1)-1.0))));
    }
  */


  TFile* filenominalCent1 = new TFile("ptshape/BDT/MCstudiesPbPb_nominal_Cent0-30-90.root");
  TFile* filenoweightCent1 = new TFile("ptshape/BDT/MCstudiesPbPb_noweight_Cent0-30-90.root");
  TFile* filenominalCent2 = new TFile("ptshape/BDT/MCstudiesPbPb_nominal_Cent0-90.root");
  TFile* filenoweightCent2 = new TFile("ptshape/BDT/MCstudiesPbPb_noweight_Cent0-90.root");

  TH1D* hEffnominalCent1 = (TH1D*)filenominalCent1->Get("hEff");
  hEffnominalCent1->GetXaxis()->CenterTitle();
  hEffnominalCent1->GetYaxis()->CenterTitle();
  hEffnominalCent1->GetXaxis()->SetTitle("hiBin");
  hEffnominalCent1->GetYaxis()->SetTitle("#alpha x #epsilon");
  hEffnominalCent1->GetXaxis()->SetTitleOffset(0.9);
  hEffnominalCent1->GetYaxis()->SetTitleOffset(0.95);
  hEffnominalCent1->GetXaxis()->SetTitleSize(0.05);
  hEffnominalCent1->GetYaxis()->SetTitleSize(0.05);
  hEffnominalCent1->GetXaxis()->SetTitleFont(42);
  hEffnominalCent1->GetYaxis()->SetTitleFont(42);
  hEffnominalCent1->GetXaxis()->SetLabelFont(42);
  hEffnominalCent1->GetYaxis()->SetLabelFont(42);
  hEffnominalCent1->GetXaxis()->SetLabelSize(0.035);
  hEffnominalCent1->GetYaxis()->SetLabelSize(0.035);
  hEffnominalCent1->SetLineColor(kRed);
  TH1D* hEffnoweightCent1 = (TH1D*)filenoweightCent1->Get("hEff");
  hEffnoweightCent1->GetXaxis()->CenterTitle();
  hEffnoweightCent1->GetYaxis()->CenterTitle();
  hEffnoweightCent1->GetXaxis()->SetTitle("hiBin");
  hEffnoweightCent1->GetYaxis()->SetTitle("#alpha x #epsilon");
  hEffnoweightCent1->GetXaxis()->SetTitleOffset(0.9);
  hEffnoweightCent1->GetYaxis()->SetTitleOffset(0.95);
  hEffnoweightCent1->GetXaxis()->SetTitleSize(0.05);
  hEffnoweightCent1->GetYaxis()->SetTitleSize(0.05);
  hEffnoweightCent1->GetXaxis()->SetTitleFont(42);
  hEffnoweightCent1->GetYaxis()->SetTitleFont(42);
  hEffnoweightCent1->GetXaxis()->SetLabelFont(42);
  hEffnoweightCent1->GetYaxis()->SetLabelFont(42);
  hEffnoweightCent1->GetXaxis()->SetLabelSize(0.035);
  hEffnoweightCent1->GetYaxis()->SetLabelSize(0.035);
  hEffnoweightCent1->SetLineColor(kBlue);
  TH1D* hEffnominalCent2 = (TH1D*)filenominalCent2->Get("hEff");
  hEffnominalCent2->GetXaxis()->CenterTitle();
  hEffnominalCent2->GetYaxis()->CenterTitle();
  hEffnominalCent2->GetXaxis()->SetTitle("hiBin");
  hEffnominalCent2->GetYaxis()->SetTitle("#alpha x #epsilon");
  hEffnominalCent2->GetXaxis()->SetTitleOffset(0.9);
  hEffnominalCent2->GetYaxis()->SetTitleOffset(0.95);
  hEffnominalCent2->GetXaxis()->SetTitleSize(0.05);
  hEffnominalCent2->GetYaxis()->SetTitleSize(0.05);
  hEffnominalCent2->GetXaxis()->SetTitleFont(42);
  hEffnominalCent2->GetYaxis()->SetTitleFont(42);
  hEffnominalCent2->GetXaxis()->SetLabelFont(42);
  hEffnominalCent2->GetYaxis()->SetLabelFont(42);
  hEffnominalCent2->GetXaxis()->SetLabelSize(0.035);
  hEffnominalCent2->GetYaxis()->SetLabelSize(0.035);
  hEffnominalCent2->SetLineColor(kRed);
  TH1D* hEffnoweightCent2 = (TH1D*)filenoweightCent2->Get("hEff");
  hEffnoweightCent2->GetXaxis()->CenterTitle();
  hEffnoweightCent2->GetYaxis()->CenterTitle();
  hEffnoweightCent2->GetXaxis()->SetTitle("hiBin");
  hEffnoweightCent2->GetYaxis()->SetTitle("#alpha x #epsilon");
  hEffnoweightCent2->GetXaxis()->SetTitleOffset(0.9);
  hEffnoweightCent2->GetYaxis()->SetTitleOffset(0.95);
  hEffnoweightCent2->GetXaxis()->SetTitleSize(0.05);
  hEffnoweightCent2->GetYaxis()->SetTitleSize(0.05);
  hEffnoweightCent2->GetXaxis()->SetTitleFont(42);
  hEffnoweightCent2->GetYaxis()->SetTitleFont(42);
  hEffnoweightCent2->GetXaxis()->SetLabelFont(42);
  hEffnoweightCent2->GetYaxis()->SetLabelFont(42);
  hEffnoweightCent2->GetXaxis()->SetLabelSize(0.035);
  hEffnoweightCent2->GetYaxis()->SetLabelSize(0.035);
  hEffnoweightCent2->SetLineColor(kBlue);
  
  /*  
  TH1D* ptshapeCent1 = (TH1D*)hEffnoweightCent1->Clone("ptshapeCent1");
  ptshapeCent1->Divide(hEffnominalCent1);

  TCanvas* c100 = new TCanvas("","",600,600);
  c100->cd();
  ptshapeCent1->SetMaximum(1.5);
  ptshapeCent1->SetMinimum(0.5);  
  ptshapeCent1->GetXaxis()->SetTitle("hiBin");
  ptshapeCent1->GetYaxis()->SetTitle("#alpha x #epsilon nominal/noweight");
  ptshapeCent1->Draw();
  c100->SaveAs("ptshape/ptshape_Cent1_BDT_datadriven.png");
  c100->SaveAs("ptshape/ptshape_Cent1_BDT_datadriven.pdf");

  for(int j=0;j<2;j++)
    {
      printf("Cent bins %.0f-%.0f ptshape uncertainty: %f (percent)\n",_ptBins[j],_ptBins[j+1],100.0*(ptshapeCent1->GetBinContent(j+1)-1.0));
    }
  */
  
  
  TH1D* ptshapeCent2 = (TH1D*)hEffnoweightCent2->Clone("ptshapeCent2");
  ptshapeCent2->Divide(hEffnominalCent2);

  TCanvas* c200 = new TCanvas("","",600,600);
  c200->cd();
  ptshapeCent2->SetMaximum(1.5);
  ptshapeCent2->SetMinimum(0.5);  
  ptshapeCent2->GetXaxis()->SetTitle("hiBin");
  ptshapeCent2->GetYaxis()->SetTitle("#alpha x #epsilon nominal/noweight");
  ptshapeCent2->Draw();
  c200->SaveAs("ptshape/ptshape_Cent2_BDT_datadriven.png");
  c200->SaveAs("ptshape/ptshape_Cent2_BDT_datadriven.pdf");

  for(int j=0;j<1;j++)
    {
      printf("Cent bins %.0f-%.0f ptshape uncertainty: %f (percent)\n",_ptBins[j],_ptBins[j+1],100.0*(ptshapeCent2->GetBinContent(j+1)-1.0));
    }  


  TH2F* hemptyEff=new TH2F("hemptyEff","",50,0.,200.,20,0.,0.2);  
  hemptyEff->GetXaxis()->CenterTitle();
  hemptyEff->GetYaxis()->CenterTitle();
  //hemptyEff->GetYaxis()->SetTitle("acceptance x #epsilon_{reco} x #epsilon_{sel} ");
  hemptyEff->GetYaxis()->SetTitle("#alpha x #epsilon");
  hemptyEff->GetXaxis()->SetTitle("hiBin");
  hemptyEff->GetXaxis()->SetTitleOffset(0.9);
  hemptyEff->GetYaxis()->SetTitleOffset(0.95);
  hemptyEff->GetXaxis()->SetTitleSize(0.05);
  hemptyEff->GetYaxis()->SetTitleSize(0.05);
  hemptyEff->GetXaxis()->SetTitleFont(42);
  hemptyEff->GetYaxis()->SetTitleFont(42);
  hemptyEff->GetXaxis()->SetLabelFont(42);
  hemptyEff->GetYaxis()->SetLabelFont(42);
  hemptyEff->GetXaxis()->SetLabelSize(0.035);
  hemptyEff->GetYaxis()->SetLabelSize(0.035);  
  hemptyEff->SetMaximum(0.5);
  hemptyEff->SetMinimum(0.);
  hemptyEff->Draw();

  TH2F* hemptyEffAcc=(TH2F*)hemptyEff->Clone("hemptyEffAcc");
  TH2F* hemptyEffReco=(TH2F*)hemptyEff->Clone("hemptyEffReco");
  TH2F* hemptyEffSelection=(TH2F*)hemptyEff->Clone("hemptyEffSelection"); 

  TCanvas*canvasEff=new TCanvas("canvasEff","canvasEff",1000.,500);
  canvasEff->Divide(2,1);
  canvasEff->cd(1);
  //gPad->SetLogy();
  hemptyEffAcc->SetYTitle("#alpha");
  hemptyEffAcc->Draw();
  hEffAcc->Draw("same");
  
  canvasEff->cd(2);
  //gPad->SetLogy();
  hemptyEff->Draw();
  hEff->Draw("same");
  canvasEff->SaveAs(Form("plotEffCent/canvasEff_study%s_Cent.png",Form(label.Data())));
  canvasEff->SaveAs(Form("plotEffCent/canvasEff_study%s_Cent.pdf",Form(label.Data())));
  
  TH2F* hemptyPthat=new TH2F("hemptyPthat","",50,0.,200.,10,1e-5,1e9);  
  hemptyPthat->GetXaxis()->CenterTitle();
  hemptyPthat->GetYaxis()->CenterTitle();
  hemptyPthat->GetYaxis()->SetTitle("Entries");
  hemptyPthat->GetXaxis()->SetTitle("pthat");
  hemptyPthat->GetXaxis()->SetTitleOffset(0.9);
  hemptyPthat->GetYaxis()->SetTitleOffset(0.95);
  hemptyPthat->GetXaxis()->SetTitleSize(0.05);
  hemptyPthat->GetYaxis()->SetTitleSize(0.05);
  hemptyPthat->GetXaxis()->SetTitleFont(42);
  hemptyPthat->GetYaxis()->SetTitleFont(42);
  hemptyPthat->GetXaxis()->SetLabelFont(42);
  hemptyPthat->GetYaxis()->SetLabelFont(42);
  hemptyPthat->GetXaxis()->SetLabelSize(0.035);
  hemptyPthat->GetYaxis()->SetLabelSize(0.035);  
  hemptyPthat->SetMaximum(2);
  hemptyPthat->SetMinimum(0.);

  TH2F* hemptySpectra=new TH2F("hemptySpectra","",50,0.,200.,10,1,1e9);  
  hemptySpectra->GetXaxis()->CenterTitle();
  hemptySpectra->GetYaxis()->CenterTitle();
  hemptySpectra->GetYaxis()->SetTitle("Entries");
  hemptySpectra->GetXaxis()->SetTitle("hiBin");
  hemptySpectra->GetXaxis()->SetTitleOffset(0.9);
  hemptySpectra->GetYaxis()->SetTitleOffset(0.95);
  hemptySpectra->GetXaxis()->SetTitleSize(0.05);
  hemptySpectra->GetYaxis()->SetTitleSize(0.05);
  hemptySpectra->GetXaxis()->SetTitleFont(42);
  hemptySpectra->GetYaxis()->SetTitleFont(42);
  hemptySpectra->GetXaxis()->SetLabelFont(42);
  hemptySpectra->GetYaxis()->SetLabelFont(42);
  hemptySpectra->GetXaxis()->SetLabelSize(0.035);
  hemptySpectra->GetYaxis()->SetLabelSize(0.035);  

  TH2F* hemptyPthatWeighted=(TH2F*)hemptyPthat->Clone("hemptyPthatWeighted");
  hemptyPthatWeighted->GetXaxis()->SetTitle("pthat reweighted");
  
  TCanvas*canvasPthat=new TCanvas("canvasPthat","canvasPthat",1000.,500);
  canvasPthat->Divide(2,1);
  canvasPthat->cd(1);
  gPad->SetLogy();
  hemptyPthat->Draw("same");
  hPthat->Draw("same");
  canvasPthat->cd(2);
  gPad->SetLogy();
  hemptyPthatWeighted->Draw();
  hPthatweight->Draw("same");
  canvasPthat->SaveAs(Form("plotEffCent/canvasPthat_%s_Cent.pdf",Form(label.Data())));
  
  TCanvas*canvasSpectra=new TCanvas("canvasSpectra","canvasSpectra",1000.,500);
  canvasSpectra->Divide(2,1);
  canvasSpectra->cd(1);
  gPad->SetLogy();
  hemptySpectra->Draw();
  hPtMC->Draw("same");
  canvasSpectra->cd(2);
  gPad->SetLogy();
  hemptySpectra->Draw();
  hPtGen->Draw("same");
  canvasSpectra->SaveAs(Form("plotEffCent/canvasSpectra_%s_Cent.pdf",Form(label.Data())));

  TCanvas*canvas1D=new TCanvas("canvas1D","",600,600);
  canvas1D->cd();
  gPad->SetLogy();
  hemptySpectra->SetYTitle("Entries of hPtMC");
  hemptySpectra->Draw(); 
  hPtMC->Draw("same");
  canvas1D->SaveAs(Form("plotEffCent/canvas1DhPtMC_%s_Cent.pdf",Form(label.Data())));
  canvas1D->Clear();

  canvas1D=new TCanvas("canvas1D","",600,600);
  canvas1D->cd();
  gPad->SetLogy();
  hemptySpectra->SetYTitle("Entries of hPtMCrecoonly");
  hemptySpectra->Draw(); 
  hPtMCrecoonly->Draw("same");
  canvas1D->SaveAs(Form("plotEffCent/canvas1DhPtMCrecoonly_%s_Cent.pdf",Form(label.Data())));
  canvas1D->Clear();

  canvas1D=new TCanvas("canvas1D","",600,600);
  canvas1D->cd();
  gPad->SetLogy();
  hemptySpectra->SetYTitle("Entries of hPtGen");
  hemptySpectra->Draw(); 
  hPtGen->Draw("same");
  canvas1D->SaveAs(Form("plotEffCent/canvas1DhPtGen_%s_Cent.pdf",Form(label.Data())));
  canvas1D->Clear();

  canvas1D=new TCanvas("canvas1D","",600,600);
  canvas1D->cd();
  gPad->SetLogy();
  hemptySpectra->SetYTitle("Entries of hPtGenAcc");
  hemptySpectra->Draw(); 
  hPtGenAcc->Draw("same");
  canvas1D->SaveAs(Form("plotEffCent/canvas1DhPtGenAcc_%s_Cent.pdf",Form(label.Data())));
  canvas1D->Clear();

  canvas1D=new TCanvas("canvas1D","",600,600);
  canvas1D->cd();
  gPad->SetLogy(0);
  hemptyEff->SetYTitle("hPtMC / hPtGen");
  hemptyEff->Draw(); 
  hEff->Draw("same");
  canvas1D->SaveAs(Form("plotEffCent/canvas1DhEff_%s_Cent.pdf",Form(label.Data())));
  canvas1D->Clear();

  canvas1D=new TCanvas("canvas1D","",600,600);
  canvas1D->cd();
  hemptyEff->SetYTitle("hPtMCrecoonly / hPtGen");
  hemptyEff->Draw(); 
  hEffReco->Draw("same");
  canvas1D->SaveAs(Form("plotEffCent/canvas1DhEffReco_%s_Cent.pdf",Form(label.Data())));
  canvas1D->Clear();

  canvas1D=new TCanvas("canvas1D","",600,600);
  canvas1D->cd();
  hemptyEff->SetYTitle("hPtGenAcc / hPtGen");
  hemptyEff->Draw(); 
  hEffAcc->Draw("same");
  canvas1D->SaveAs(Form("plotEffCent/canvas1DhEffAcc_%s_Cent.pdf",Form(label.Data())));
  canvas1D->Clear();

  canvas1D=new TCanvas("canvas1D","",600,600);
  canvas1D->cd();
  hemptyEff->SetYTitle("hPtMC / hPtGenAcc");
  hemptyEff->Draw(); 
  hEffSelection->Draw("same");
  canvas1D->SaveAs(Form("plotEffCent/canvas1DhEffSelection_%s_Cent.pdf",Form(label.Data())));
  canvas1D->Clear();

  gStyle->SetPalette(55);
  TCanvas* canvas2D=new TCanvas("canvas2D","",600,600);

  TFile *fout=new TFile(outputfile.Data(),"recreate");
  fout->cd();
  hPtGen->Write();
  hEffAcc->Write();
  hEffReco->Write();
  hEffSelection->Write();
  hEff->Write();
  hPtMC->Write();
  fout->Close();  

}

int main(int argc, char *argv[])
{
  if((argc !=12))
  {
    std::cout << "Wrong number of inputs" << std::endl;
    return 1;
  }
  
  if(argc == 12)
    MCefficiencyCent(atoi(argv[1]),argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],atoi(argv[9]),atof(argv[10]),atof(argv[11]));
  return 0;
}
