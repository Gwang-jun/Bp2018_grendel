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

int _nBins = nBins;
double *_ptBins = ptBins;

void MCefficiency(int isPbPb=0,TString inputmc="", TString selmcgen="",TString selmcgenacceptance="", TString cut_recoonly="", TString cut="",TString label="PP",TString outputfile="test", int useweight=1,Float_t centmin=0., Float_t centmax=100.)
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
      selmcgen = selmcgen+Form("&&hiBin>=%f&&hiBin<=%f",hiBinMin,hiBinMax);
      selmcgenacceptance=selmcgenacceptance+Form("&&hiBin>=%f&&hiBin<=%f",hiBinMin,hiBinMax);
      cut_recoonly=cut_recoonly+Form("&&hiBin>=%f&&hiBin<=%f",hiBinMin,hiBinMax);
      cut=cut+Form("&&hiBin>=%f&&hiBin<=%f",hiBinMin,hiBinMax);
    }
  
     std::cout<<"selmcgen="<<selmcgen<<std::endl;
     std::cout<<"selmcgenacceptance="<<selmcgenacceptance<<std::endl;
     std::cout<<"cut_recoonly="<<cut_recoonly<<std::endl;
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
  TTree* ntGen = (TTree*)infMC->Get("Bfinder/ntGen");
  ntGen->AddFriend("hltanalysis/HltTree");
  ntGen->AddFriend("hiEvtAnalyzer/HiTree");
  ntGen->AddFriend("Bfinder/ntKp");
  ntGen->AddFriend("skimanalysis/HltTree");
  
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
  
  // optimal weights
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
    weightGpt = "(2.907795+-0.436572*Gpt+0.006372*Gpt*Gpt)*TMath::Exp(-0.157563*Gpt)+1.01308";
    weightBgenpt = "(2.907795+-0.436572*Bgenpt+0.006372*Bgenpt*Bgenpt)*TMath::Exp(-0.157563*Bgenpt)+1.01308";    
    //weightGpt = "0.889175+0.000791*Gpt+0.000015*Gpt*Gpt";
    //weightBgenpt = "0.889175+0.000791*Bgenpt+0.000015*Bgenpt*Bgenpt";
    //weightGpt = "0.094376+0.028350*Gpt+-0.000225*Gpt*Gpt+5.369348/Gpt";
    //weightBgenpt = "0.094376+0.028350*Bgenpt+-0.000225*Bgenpt*Bgenpt+5.369348/Bgenpt";
    //weightPVz = "1.034350*TMath::Exp(-0.000844*(PVz+3.502992)*(PVz+3.502992))"; // private MC
    //weightGpt = "0.715021+0.039896*Gpt-0.000834*Gpt*Gpt+0.000006*Gpt*Gpt*Gpt"; // private MC Gpt
    //weightBgenpt = "0.715021+0.039896*Bgenpt-0.000834*Bgenpt*Bgenpt+0.000006*Bgenpt*Bgenpt*Bgenpt"; // private MC Gpt
    }
  
  TH1D* hPtMC = new TH1D("hPtMC","",_nBins,_ptBins);
  TH1D* hPtMCrecoonly = new TH1D("hPtMCrecoonly","",_nBins,_ptBins);
  TH1D* hPtGen = new TH1D("hPtGen","",_nBins,_ptBins);
  TH1D* hPtGenAcc = new TH1D("hPtGenAcc","",_nBins,_ptBins);
  TH1D* hPtGenAccWeighted = new TH1D("hPtGenAccWeighted","",_nBins,_ptBins);
  TH1D* hPthat = new TH1D("hPthat","",100,0,500);
  TH1D* hPthatweight = new TH1D("hPthatweight","",100,0,500);

  ntMC->Project("hPtMC","Bpt",TCut(weighpthat)*TCut(weightBgenpt)*TCut(weightHiBin)*TCut(weightPVz)*(TCut(cut.Data())&&"(Bgen==23333)"));
  ntMC->Project("hPtMCrecoonly","Bpt",TCut(weighpthat)*TCut(weightBgenpt)*TCut(weightHiBin)*TCut(weightPVz)*(TCut(cut_recoonly.Data())&&"(Bgen==23333)"));
  ntGen->Project("hPtGen","Gpt",TCut(weighpthat)*TCut(weightGpt)*(TCut(selmcgen.Data())));
  ntGen->Project("hPtGenAcc","Gpt",TCut(weighpthat)*TCut(weightGpt)*(TCut(selmcgenacceptance.Data())));
  ntGen->Project("hPtGenAccWeighted","Gpt",TCut(weighpthat)*TCut(weightGpt)*TCut(weightHiBin)*TCut(weightPVz)*(TCut(selmcgenacceptance.Data())));
  
  ////// tag & probe scaling factor
  for(int i = 0; i < _nBins; i++){printf("%.2f, ", hPtMC->GetBinContent(i+1));}printf("\n");
  double sf_pp[7] = {1., 227977.02/207901.56, 632864.53/606712.42, 408808.07/399583.96, 275911.15/272909.48, 85362.85/85846.52, 1.};//SUBJECT TO CHANGE!
  double sf_pbpb[7] = {1., 64266.93/59877.61, 189085.92/187577.36, 156249.18/158435.51, 138502.42/141744.96, 53617.98/55196.42, 1.};//SUBJECT TO CHANGE!
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
  for(int i = 0; i < _nBins; i++){printf("%.2f, ", hPtMC->GetBinContent(i+1));}printf("\n");

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
  /*
  //hEffAcc = hPtGenAcc / hPtGen
  //hEffReco = hPtMCrecoonly / hPtGenAcc
  //hEffSelection = hPtMC / hPtMCrecoonly
  //hEff = hPtMC / hPtGen

  //Acceptance only
  TH1D* hEffAcc = (TH1D*)hPtGenAcc->Clone("hEffAcc");
  hEffAcc->Sumw2();
  hEffAcc->Divide(hEffAcc,hPtGen,1,1,"b");
  //Eff_reco only
  TH1D* hEffReco = (TH1D*)hPtMCrecoonly->Clone("hEffReco");
  hEffReco->Sumw2();
  hEffReco->Divide(hEffReco,hPtGenAcc,1,1,"b");
  hEffReco->Divide(hEffReco,hPtGen,1,1,"b");
  //Eff_selection only
  TH1D* hEffSelection = (TH1D*)hPtMC->Clone("hEffSelection");
  hEffSelection->Sumw2();
  hEffSelection->Divide(hEffSelection,hPtMCrecoonly,1,1,"b");
  //Full efficiency = Acc * Eff_reco * Eff_selection
  TH1D* hEff = (TH1D*)hPtMC->Clone("hEff");
  hEff->Sumw2();
  hEff->Divide(hPtGen);
  */

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

  TFile* filenoweight = new TFile("MCstudiesPbPb_noweight.root");
  TH1D* hEffnoweight = (TH1D*)filenoweight->Get("hEff");
  TH1D* ptshape = (TH1D*)hEff->Clone("ptshape");
  ptshape->Sumw2();
  ptshape->Divide(hEffnoweight);

  TCanvas* c100 = new TCanvas("","",600,600);
  c100->cd();
  ptshape->GetXaxis()->SetTitle("Genpt (GeV/c)");
  ptshape->GetYaxis()->SetTitle("Efficiency pt weight/noweight");
  ptshape->Draw();
  c100->SaveAs("ptshape.png");

  for(int j=0;j<_nBins;j++)
    {
      printf("p_t bins %.0f-%.0f ptshape uncertainty: %.2f (percent)\n",_ptBins[j],_ptBins[j+1],100.0*(ptshape->GetBinContent(j+1)-1.0));
    }

  TH2F* hemptyEff=new TH2F("hemptyEff","",50,_ptBins[0]-5.,_ptBins[_nBins]+5.,10,0.,1.);  
  hemptyEff->GetXaxis()->CenterTitle();
  hemptyEff->GetYaxis()->CenterTitle();
  //hemptyEff->GetYaxis()->SetTitle("acceptance x #epsilon_{reco} x #epsilon_{sel} ");
  hemptyEff->GetYaxis()->SetTitle("#alpha x #epsilon");
  hemptyEff->GetXaxis()->SetTitle("p_{T} (GeV/c)");
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
  canvasEff->SaveAs(Form("plotEff/canvasEff_study%s.png",Form(label.Data())));
  canvasEff->SaveAs(Form("plotEff/canvasEff_study%s.pdf",Form(label.Data())));
  
  for(int j=0;j<_nBins;j++)
    {
      printf("p_t bins %.0f-%.0f hEff: %f hEffErr: %f\n",_ptBins[j],_ptBins[j+1],hEff->GetBinContent(j+1),hEff->GetBinError(j+1));
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
  
  TH2F* hemptySpectra=new TH2F("hemptySpectra","",50,0.,130.,10,1,1e9);  
  hemptySpectra->GetXaxis()->CenterTitle();
  hemptySpectra->GetYaxis()->CenterTitle();
  hemptySpectra->GetYaxis()->SetTitle("Entries");
  hemptySpectra->GetXaxis()->SetTitle("p_{T}");
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
  canvasPthat->SaveAs(Form("plotEff/canvasPthat_%s.pdf",Form(label.Data())));
  
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
  canvasSpectra->SaveAs(Form("plotEff/canvasSpectra_%s.pdf",Form(label.Data())));
  
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
  TLatex* texCent = new TLatex(0.5,0.815,Form("Centrality %.0f-%.0f%s",centMin,centMax,texper.Data()));
  texCent->SetNDC();
  texCent->SetTextFont(42);
  texCent->SetTextSize(0.05);
  
  TLatex* texY = new TLatex(0.5,0.750,"|y| < 2.4");
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
  hemptyEff->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hemptyEff->Draw();
  hEff->SetLineColor(2);
  hEff->SetMarkerColor(2);
  hEff->Draw("same");
  if(isPbPb) texCent->Draw();
  texY->Draw();
  texCms->Draw();
  texlumi->Draw();
  canvas1D->SaveAs(Form("plotEff/canvas1DhEff_%s.pdf",Form(label.Data())));
  canvas1D->SaveAs(Form("plotEff/canvas1DhEff_%s.C",Form(label.Data())));
  canvas1D->Clear();
  
  canvas1D=new TCanvas("canvas1D","",600,600);
  canvas1D->cd();
  gPad->SetLogy();
  hemptySpectra->SetYTitle("Entries of hPtMC");
  hemptySpectra->Draw(); 
  hPtMC->Draw("same");
  canvas1D->SaveAs(Form("plotEff/canvas1DhPtMC_%s.pdf",Form(label.Data())));
  canvas1D->Clear();
  
  canvas1D=new TCanvas("canvas1D","",600,600);
  canvas1D->cd();
  gPad->SetLogy();
  hemptySpectra->SetYTitle("Entries of hPtMCrecoonly");
  hemptySpectra->Draw(); 
  hPtMCrecoonly->Draw("same");
  canvas1D->SaveAs(Form("plotEff/canvas1DhPtMCrecoonly_%s.pdf",Form(label.Data())));
  canvas1D->Clear();
  
  canvas1D=new TCanvas("canvas1D","",600,600);
  canvas1D->cd();
  gPad->SetLogy();
  hemptySpectra->SetYTitle("Entries of hPtGen");
  hemptySpectra->Draw(); 
  hPtGen->Draw("same");
  canvas1D->SaveAs(Form("plotEff/canvas1DhPtGen_%s.pdf",Form(label.Data())));
  canvas1D->Clear();
  
  canvas1D=new TCanvas("canvas1D","",600,600);
  canvas1D->cd();
  gPad->SetLogy();
  hemptySpectra->SetYTitle("Entries of hPtGenAcc");
  hemptySpectra->Draw(); 
  hPtGenAcc->Draw("same");
  canvas1D->SaveAs(Form("plotEff/canvas1DhPtGenAcc_%s.pdf",Form(label.Data())));
  canvas1D->Clear();

  canvas1D=new TCanvas("canvas1D","",600,600);
  canvas1D->cd();
  gPad->SetLogy();
  hemptyEff->SetYTitle("hPtMC / hPtGen");
  hemptyEff->Draw();
  hEff->Draw("same");
  canvas1D->SaveAs(Form("plotEff/canvas1DhEff_%s.pdf",Form(label.Data())));
  canvas1D->Clear();
  
  canvas1D=new TCanvas("canvas1D","",600,600);
  canvas1D->cd();
  hemptyEff->SetYTitle("hPtMCrecoonly / hPtGen");
  hemptyEff->Draw(); 
  hEffReco->Draw("same");
  canvas1D->SaveAs(Form("plotEff/canvas1DhEffReco_%s.pdf",Form(label.Data())));
  canvas1D->Clear();
  
  canvas1D=new TCanvas("canvas1D","",600,600);
  canvas1D->cd();
  hemptyEff->SetYTitle("hPtGenAcc / hPtGen");
  hemptyEff->Draw(); 
  hEffAcc->Draw("same");
  canvas1D->SaveAs(Form("plotEff/canvas1DhEffAcc_%s.pdf",Form(label.Data())));
  canvas1D->Clear();
  
  canvas1D=new TCanvas("canvas1D","",600,600);
  canvas1D->cd();
  hemptyEff->SetYTitle("hPtMC / hPtGenAcc");
  hemptyEff->Draw(); 
  hEffSelection->Draw("same");
  canvas1D->SaveAs(Form("plotEff/canvas1DhEffSelection_%s.pdf",Form(label.Data())));
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
    MCefficiency(atoi(argv[1]),argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],atoi(argv[9]),atof(argv[10]),atof(argv[11]));
  return 0;
}
