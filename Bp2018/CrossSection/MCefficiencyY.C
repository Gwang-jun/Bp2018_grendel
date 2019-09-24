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

TString weightfunctiongen = "1";
TString weightfunctionreco = "1";
Float_t hiBinMin,hiBinMax,centMin,centMax;

int _nBins = nBinsYFine;
double *_ptBins = ptBinsYFine;

void MCefficiencyY(int isPbPb=1,TString inputmc="", TString selmcgen="",TString selmcgenacceptance="", TString cut_recoonly="", TString cut="",TString label="",TString outputfile="", int useweight=1, Float_t centmin=0., Float_t centmax=90.)
{ 
  if(label=="ppInc"){
    _nBins = nBinsInc;
    _ptBins = ptBinsInc;
  }
  
  std::string str = label.Data();
  std::size_t found = str.find("Fine");
  //if(label=="ppFine"||label=="PbPbFine"){
  if (found!=std::string::npos){
    _nBins = nBinsFine;
    _ptBins = ptBinsFine;
  }

  hiBinMin = centmin*2;
  hiBinMax = centmax*2;
  centMin = centmin;
  centMax = centmax;
  
  if(isPbPb==1)
    {
      selmcgen = selmcgen+Form("&&hiBin>=%f&&hiBin<=%f&&Gpt>%f&&Gpt<%f",hiBinMin,hiBinMax,ptBinsInc[0],ptBinsInc[1]);
      selmcgenacceptance=selmcgenacceptance+Form("&&hiBin>=%f&&hiBin<=%f&&Gpt>%f&&Gpt<%f",hiBinMin,hiBinMax,ptBinsInc[0],ptBinsInc[1]);
      cut_recoonly=cut_recoonly+Form("&&hiBin>=%f&&hiBin<=%f&&Bpt>%f&&Bpt<%f",hiBinMin,hiBinMax,ptBinsInc[0],ptBinsInc[1]);
      cut=cut+Form("&&hiBin>=%f&&hiBin<=%f&&Bpt>%f&&Bpt<%f",hiBinMin,hiBinMax,ptBinsInc[0],ptBinsInc[1]);
    }

     std::cout<<"selmcgen="<<selmcgen<<std::endl;
     std::cout<<"selmcgenacceptance="<<selmcgenacceptance<<std::endl;
     std::cout<<"cut_recoonly"<<cut_recoonly<<std::endl;
     std::cout<<"cut="<<cut<<std::endl;
     std::cout<<"option="<<useweight<<std::endl;

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

  /*
  //For 2015 PbPb, pp MC
  TFile* infMC = new TFile(inputmc.Data());
  TTree* ntMC = (TTree*)infMC->Get("ntKp");
  ntMC->AddFriend("ntHlt");
  ntMC->AddFriend("ntHi");
  ntMC->AddFriend("ntGen");
  ntMC->AddFriend("ntSkim");
  TTree* ntGen = (TTree*)infMC->Get("ntGen");
  ntGen->AddFriend("ntHlt");
  ntGen->AddFriend("ntHi");
  ntGen->AddFriend("ntKp");
  ntGen->AddFriend("ntSkim");
  */

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
    weightGpt = "(3.00448277-0.35865276*Gpt+0.01997413*Gpt*Gpt-0.00042585*Gpt*Gpt*Gpt+0.00000315*Gpt*Gpt*Gpt*Gpt)";
    weightBgenpt = "(3.00448277-0.35865276*Bgenpt+0.01997413*Bgenpt*Bgenpt-0.00042585*Bgenpt*Bgenpt*Bgenpt+0.00000315*Bgenpt*Bgenpt*Bgenpt*Bgenpt)";
    //weightGpt = "(3.506006+0.963473*Gpt+-0.258731*Gpt*Gpt)*TMath::Exp(-0.386065*Gpt)+1.139897";
    //weightBgenpt = "(3.506006+0.963473*Bgenpt+-0.258731*Bgenpt*Bgenpt)*TMath::Exp(-0.386065*Bgenpt)+1.139897";
  }

  std::cout<<"fit function parameters="<<weightfunctiongen<<std::endl;

  TH1D* hPtMC = new TH1D("hPtMC","",_nBins,_ptBins);
  TH1D* hPtMCrecoonly = new TH1D("hPtMCrecoonly","",_nBins,_ptBins);
  TH1D* hPtGen = new TH1D("hPtGen","",_nBins,_ptBins);
  TH1D* hPtGenAcc = new TH1D("hPtGenAcc","",_nBins,_ptBins);
  TH1D* hPtGenAccWeighted = new TH1D("hPtGenAccWeighted","",_nBins,_ptBins);
  TH1D* hPthat = new TH1D("hPthat","",100,0,500);
  TH1D* hPthatweight = new TH1D("hPthatweight","",100,0,500);

  ntMC->Project("hPtMC","TMath::Abs(By)",TCut(weighpthat)*TCut(weightBgenpt)*TCut(weightHiBin)*TCut(weightPVz)*(TCut(cut.Data())&&"(Bgen==23333)"));
  ntMC->Project("hPtMCrecoonly","TMath::Abs(By)",TCut(weighpthat)*TCut(weightBgenpt)*TCut(weightHiBin)*TCut(weightPVz)*(TCut(cut_recoonly.Data())&&"(Bgen==23333)"));
  ntGen->Project("hPtGen","TMath::Abs(Gy)",TCut(weighpthat)*TCut(weightGpt)*(TCut(selmcgen.Data())));
  ntGen->Project("hPtGenAcc","TMath::Abs(Gy)",TCut(weighpthat)*TCut(weightGpt)*(TCut(selmcgenacceptance.Data())));
  ntGen->Project("hPtGenAccWeighted","TMath::Abs(Gy)",TCut(weighpthat)*TCut(weightGpt)*TCut(weightHiBin)*TCut(weightPVz)*(TCut(selmcgenacceptance.Data())));

  divideBinWidth(hPtMC);
  divideBinWidth(hPtMCrecoonly);
  divideBinWidth(hPtGen);
  divideBinWidth(hPtGenAcc);
  divideBinWidth(hPtGenAccWeighted);

  ntMC->Project("hPthat","pthat","1");
  ntMC->Project("hPthatweight","pthat","pthatweight");

  //Acceptance
  TH1D* hEffAcc = (TH1D*)hPtGenAcc->Clone("hEffAcc");
  hEffAcc->Divide(hEffAcc,hPtGen,1,1,"b");
  //Selection
  TH1D* hEffSelection = (TH1D*)hPtMC->Clone("hEffSelection");
  hEffSelection->Divide(hEffSelection,hPtGenAccWeighted,1,1,"b");
  //Acc * Eff (one shot)
  TH1D* hEffReco = (TH1D*)hPtMCrecoonly->Clone("hEffReco");
  hEffReco->Divide(hEffReco,hPtGen,1,1,"b");
  //Acc * Eff
  TH1D* hEff = (TH1D*)hEffSelection->Clone("hEff");
  hEff->Multiply(hEff,hEffAcc,1,1);

  TH2F* hemptyEff=new TH2F("hemptyEff","",50,_ptBins[0]-1.,_ptBins[_nBins]+1.,10.,0.,1.);  
  hemptyEff->GetXaxis()->CenterTitle();
  hemptyEff->GetYaxis()->CenterTitle();
  //hemptyEff->GetYaxis()->SetTitle("acceptance x #epsilon_{reco} x #epsilon_{sel} ");
  hemptyEff->GetYaxis()->SetTitle("#alpha x #epsilon");
  hemptyEff->GetXaxis()->SetTitle("|y|");
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
  hemptyEff->SetMaximum(1.0);
  hemptyEff->SetMinimum(0.);
  hemptyEff->Draw();

  TH2F* hemptyEffAcc=(TH2F*)hemptyEff->Clone("hemptyEffAcc");
  TH2F* hemptyEffReco=(TH2F*)hemptyEff->Clone("hemptyEffReco");
  TH2F* hemptyEffSelection=(TH2F*)hemptyEff->Clone("hemptyEffSelection");
 
  TCanvas*canvasEff=new TCanvas("canvasEff","canvasEff",1000.,500);
  canvasEff->Divide(2,1);
  canvasEff->cd(1);
  
  hemptyEffAcc->SetYTitle("#alpha");
  hemptyEffAcc->Draw();
  hEffAcc->Draw("same");
  canvasEff->cd(2);
  hemptyEff->Draw();
  hEff->Draw("same");
  canvasEff->SaveAs(Form("plotEffY/canvasEff_study%s.pdf",Form(label.Data())));

  for(int j=0;j<_nBins;j++)
    {
      printf("Y bins %.1f-%.1f hEff: %f hEffErr: %f\n",_ptBins[j],_ptBins[j+1],hEff->GetBinContent(j+1),hEff->GetBinError(j+1));
    }
  
  TH2F* hemptyPthat=new TH2F("hemptyPthat","",50,0.,500.,10,1e-5,1e9);  
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
  hemptyPthat->SetMaximum(1.0);
  hemptyPthat->SetMinimum(0.);


  TH2F* hemptySpectra=new TH2F("hemptySpectra","",50,0.,130.,10,1,1e5);  
  hemptySpectra->GetXaxis()->CenterTitle();
  hemptySpectra->GetYaxis()->CenterTitle();
  hemptySpectra->GetYaxis()->SetTitle("Entries");
  hemptySpectra->GetXaxis()->SetTitle("|y|");
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
  canvasPthat->SaveAs(Form("plotEffY/canvasPthat_%s.pdf",Form(label.Data())));
  
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
  canvasSpectra->SaveAs(Form("plotEffY/canvasSpectra_%s.pdf",Form(label.Data())));

  //### 1D histogram
  //hEff = hPtMC / hPtGen
  //hEffReco = hPtMCrecoonly / hPtGen
  //hEffAcc = hPtGenAcc / hPtGen
  //hEffSelection = hPtMC / hPtMCrecoonly
 
/*
  ntMC->Project("hPtMC","TMath::Abs(By)",TCut(weightfunctionreco)*(TCut(cut.Data())&&"(Bgen==23333)"));
  ntMC->Project("hPtMCrecoonly","TMath::Abs(By)",TCut(weightfunctionreco)*(TCut(cut_recoonly.Data())&&"(Bgen==23333)"));
  ntGen->Project("hPtGen","Gy",TCut(weightfunctiongen)*(TCut(selmcgen.Data())));
  ntGen->Project("hPtGenAcc","Gy",TCut(weightfunctiongen)*(TCut(selmcgenacceptance.Data())));
*/  
  TString text;
  if (isPbPb) { text="1.5 nb^{-1} (5.02 TeV PbPb)";}
  else {text="28.0 pb^{-1} (5.02 TeV pp)";}
  TLatex* texlumi = new TLatex(0.9,0.92,text.Data());
  texlumi->SetNDC();
  texlumi->SetTextAlign(31);
  texlumi->SetTextFont(42);
  texlumi->SetTextSize(0.038);
  texlumi->SetLineWidth(2);

  TString texper="%";
  TLatex* texCent = new TLatex(0.5,0.815,Form("Cent. %.0f - %.0f%s",centMin,centMax,texper.Data()));
  texCent->SetNDC();
  texCent->SetTextFont(42);
  texCent->SetTextSize(0.05);

  TLatex* texY = new TLatex(0.5,0.750,Form("%.0f < p_{T} < %.0f (GeV/c)",ptBinsInc[0], ptBinsInc[1]));
  texY->SetNDC();
  texY->SetTextFont(42);
  texY->SetTextSize(0.05);
  texY->SetLineWidth(2);

  TLatex* texCms = new TLatex(0.08,0.95, "CMS Preliminary");
  texCms->SetNDC();
  texCms->SetTextAlign(13);
  texCms->SetTextSize(0.038);
  texCms->SetTextFont(62);

  TCanvas*canvas1D=new TCanvas("canvas1D","",600,600);
  canvas1D->cd();
  gPad->SetLogy(0);
  //hemptyEff->SetYTitle("hPtMC / hPtGen");
  hemptyEff->SetMaximum(0.6);
  hemptyEff->SetYTitle("#alpha x #epsilon");
  hemptyEff->GetXaxis()->SetTitle("|y|");
  hemptyEff->Draw();
  hEff->SetLineColor(2);
  hEff->SetMarkerColor(2);
  hEff->Draw("same");
  if(isPbPb) texCent->Draw();
  texY->Draw();
  texCms->Draw();
  texlumi->Draw();
  canvas1D->SaveAs(Form("plotEffY/canvas1DhEff_%s.pdf",Form(label.Data())));
  canvas1D->SaveAs(Form("plotEffY/canvas1DhEff_%s.C",Form(label.Data())));
  canvas1D->Clear();

  canvas1D=new TCanvas("canvas1D","",600,600);
  canvas1D->cd();
  gPad->SetLogy();
  hemptySpectra->SetYTitle("Entries of hPtMC");
  hemptySpectra->Draw(); 
  hPtMC->Draw("same");
  canvas1D->SaveAs(Form("plotEffY/canvas1DhPtMC_%s.pdf",Form(label.Data())));
  canvas1D->Clear();

  canvas1D=new TCanvas("canvas1D","",600,600);
  canvas1D->cd();
  gPad->SetLogy();
  hemptySpectra->SetYTitle("Entries of hPtMCrecoonly");
  hemptySpectra->Draw(); 
  hPtMCrecoonly->Draw("same");
  canvas1D->SaveAs(Form("plotEffY/canvas1DhPtMCrecoonly_%s.pdf",Form(label.Data())));
  canvas1D->Clear();

  canvas1D=new TCanvas("canvas1D","",600,600);
  canvas1D->cd();
  gPad->SetLogy();
  hemptySpectra->SetYTitle("Entries of hPtGen");
  hemptySpectra->Draw(); 
  hPtGen->Draw("same");
  canvas1D->SaveAs(Form("plotEffY/canvas1DhPtGen_%s.pdf",Form(label.Data())));
  canvas1D->Clear();

  canvas1D=new TCanvas("canvas1D","",600,600);
  canvas1D->cd();
  gPad->SetLogy();
  hemptySpectra->SetYTitle("Entries of hPtGenAcc");
  hemptySpectra->Draw(); 
  hPtGenAcc->Draw("same");
  canvas1D->SaveAs(Form("plotEffY/canvas1DhPtGenAcc_%s.pdf",Form(label.Data())));
  canvas1D->Clear();

  canvas1D=new TCanvas("canvas1D","",600,600);
  canvas1D->cd();
  hemptyEff->SetYTitle("hPtMCrecoonly / hPtGen");
  hemptyEff->Draw(); 
  hEffReco->Draw("same");
  canvas1D->SaveAs(Form("plotEffY/canvas1DhEffReco_%s.pdf",Form(label.Data())));
  canvas1D->Clear();

  canvas1D=new TCanvas("canvas1D","",600,600);
  canvas1D->cd();
  hemptyEff->SetYTitle("hPtGenAcc / hPtGen");
  hemptyEff->Draw(); 
  hEffAcc->Draw("same");
  canvas1D->SaveAs(Form("plotEffY/canvas1DhEffAcc_%s.pdf",Form(label.Data())));
  canvas1D->Clear();

  canvas1D=new TCanvas("canvas1D","",600,600);
  canvas1D->cd();
  hemptyEff->SetYTitle("hPtMC / hPtGenAcc");
  hemptyEff->Draw(); 
  hEffSelection->Draw("same");
  canvas1D->SaveAs(Form("plotEffY/canvas1DhEffSelection_%s.pdf",Form(label.Data())));
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
    MCefficiencyY(atoi(argv[1]),argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],atoi(argv[9]),atof(argv[10]),atof(argv[11]));
  return 0;
}




