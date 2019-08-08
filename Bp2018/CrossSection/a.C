#include "uti.h"
#include "parameters.h"
#include "TF1.h"

int main()
{
  TString inputdata = "/raid5/data/gwangjun/crab_Bfinder_20190221_Pythia8_BuToJpsiK_Bpt0p0_1032_NoJSON_pthatweight.root";
  TFile* input = new TFile(inputdata.Data());
  TTree* nt = (TTree*) input->Get("hiEvtAnalyzer/HiTree");
  TH1D* hat = new TH1D("hat","",nBinsFine,ptBinsFine);
  
  nt->Project("hat","pthat", "pthatweight");

  TCanvas* c = new TCanvas("c","",600,600);
  c->cd();
  hat->Draw();
  
  c->SaveAs("test_a.pdf");
  


  return 0;
}
