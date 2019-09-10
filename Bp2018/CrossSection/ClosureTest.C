#include "uti.h"
#include "TLegendEntry.h"
#include "parameters.h"

void ClosureTest(TString inputfile="",TString label="")
{
  gStyle->SetOptStat(0);

  
  TFile *fyield = new TFile(inputfile.Data());

  TH1D *hPtCor = (TH1D*)fyield->Get("hPtCor");
  TH1D *hPtGen = (TH1D*)fyield->Get("hPtGen");

  TH1D *hClosure = new TH1D("","",nBins,ptBins);

  for(int i=0;i<nBins;i++)
    {
      hClosure->SetBinContent(i+1,hPtCor->GetBinContent(i+1)/hPtGen->GetBinContent(i+1));
      hClosure->SetBinError(i+1,sqrt(TMath::Abs((hPtCor->GetBinError(i+1))*(hPtCor->GetBinError(i+1))-(hPtGen->GetBinError(i+1))*(hPtGen->GetBinError(i+1))))/hPtGen->GetBinContent(i+1));

      printf("Corrected Yield/Generated Yield of pt %.0f-%.0f: %f #pm %f\n",ptBins[i],ptBins[i+1],hClosure->GetBinContent(i+1),hClosure->GetBinError(i+1));
    }

  TCanvas*canvas=new TCanvas("canvas","canvas",550,500);
  canvas->cd();
  TH2F* hemptyClosure=new TH2F("hemptyClosure","",50,ptBins[0]-5,ptBins[nBins]+10,10,0.5,1.5);  
  hemptyClosure->GetXaxis()->CenterTitle();
  hemptyClosure->GetYaxis()->CenterTitle();
  hemptyClosure->GetYaxis()->SetTitle("Corrected Yield/Generated Yield");
  hemptyClosure->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hemptyClosure->GetXaxis()->SetTitleOffset(.9);
  hemptyClosure->GetYaxis()->SetTitleOffset(1.);
  hemptyClosure->GetXaxis()->SetTitleSize(0.05);
  hemptyClosure->GetYaxis()->SetTitleSize(0.05);
  hemptyClosure->GetXaxis()->SetTitleFont(42);
  hemptyClosure->GetYaxis()->SetTitleFont(42);
  hemptyClosure->GetXaxis()->SetLabelFont(42);
  hemptyClosure->GetYaxis()->SetLabelFont(42);
  hemptyClosure->GetXaxis()->SetLabelSize(0.035);
  hemptyClosure->GetYaxis()->SetLabelSize(0.035);  
  hemptyClosure->Draw();
  hClosure->SetLineWidth(2);
  hClosure->SetLineColor(1);
  hClosure->Draw("same");
  canvas->SaveAs(Form("canvasClosure%s.pdf",label.Data()));
  

  /*
  TFile *fyield = new TFile(inputfile.Data());
  TFile *fEff = new TFile("ROOTfiles/MCstudiesPbPb.root");

  TH1D *hPt = (TH1D*)fyield->Get("hPt");
  TH1D *hEff = (TH1D*)fEff->Get("hEffdef");
  TH1D *hPtGen = (TH1D*)fEff->Get("hPtGen");

  hPt->Sumw2();
  hEff->Sumw2();
  hPt->Divide(hEff);

  TH1D *hClosure = new TH1D("","",nBins,ptBins);

  for(int i=0;i<nBins;i++)
    {
      hClosure->SetBinContent(i+1,hPt->GetBinContent(i+1)/hPtGen->GetBinContent(i+1));
      hClosure->SetBinError(i+1,sqrt(TMath::Abs((hPt->GetBinError(i+1))*(hPt->GetBinError(i+1))-(hPtGen->GetBinError(i+1))*(hPtGen->GetBinError(i+1))))/hPtGen->GetBinContent(i+1));

      printf("Corrected Yield/Generated Yield of pt %.0f-%.0f: %f #pm %f\n",ptBins[i],ptBins[i+1],hClosure->GetBinContent(i+1),hClosure->GetBinError(i+1));
    }

  TCanvas*canvas=new TCanvas("canvas","canvas",550,500);
  canvas->cd();
  TH2F* hemptyClosure=new TH2F("hemptyClosure","",50,ptBins[0]-5,ptBins[nBins]+10,10,0.5,1.5);  
  hemptyClosure->GetXaxis()->CenterTitle();
  hemptyClosure->GetYaxis()->CenterTitle();
  hemptyClosure->GetYaxis()->SetTitle("Corrected Yield/Generated Yield");
  hemptyClosure->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hemptyClosure->GetXaxis()->SetTitleOffset(.9);
  hemptyClosure->GetYaxis()->SetTitleOffset(1.);
  hemptyClosure->GetXaxis()->SetTitleSize(0.05);
  hemptyClosure->GetYaxis()->SetTitleSize(0.05);
  hemptyClosure->GetXaxis()->SetTitleFont(42);
  hemptyClosure->GetYaxis()->SetTitleFont(42);
  hemptyClosure->GetXaxis()->SetLabelFont(42);
  hemptyClosure->GetYaxis()->SetLabelFont(42);
  hemptyClosure->GetXaxis()->SetLabelSize(0.035);
  hemptyClosure->GetYaxis()->SetLabelSize(0.035);  
  hemptyClosure->Draw();
  hClosure->SetLineWidth(2);
  hClosure->SetLineColor(1);
  hClosure->Draw("same");
  canvas->SaveAs(Form("canvasClosure%s.pdf",label.Data()));
  */
}

int main(int argc, char *argv[])
{
  if((argc !=3))
  {
    std::cout << "Wrong number of inputs" << std::endl;
    return 1;
  }
  
  if(argc == 3)
    ClosureTest(argv[1],argv[2]);
  return 0;
}
