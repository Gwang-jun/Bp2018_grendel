#include "uti.h"
#include "parameters.h"
#include "TF1.h"

Int_t istotal_xsection = 1;

TString outputfile_total="fonll_total.root";
TString outputfile_pt="fonll_pt.root";

TGraphAsymmErrors *fonll_total()
{
  TString outputf;
  outputf = Form("%s", outputfile_total.Data());
  TFile* outf = new TFile(outputf.Data(),"recreate");
  outf->cd();

  TCanvas *c = new TCanvas("c", "FONLL result", 600, 600);
 
  const Int_t n = nBins;
  Double_t x[n], y[n], exl[n], eyl[n], exh[n], eyh[n];

  const Double_t preset1[n] = {6., 8.5, 12.5, 17.5, 25., 40., 75.};
  const Double_t preset2[n] = {2.2580e+07, 1.6790e+07, 9.4230e+06, 2.7090e+06, 1.3370e+06, 3.3470e+05, 4.1610e+04};
  const Double_t preset3[n] = {1., 1.5, 2.5, 2.5, 5., 10., 25.};
  const Double_t preset4[n] = {6.9050e+06, 4.3930e+06, 2.2631e+06, 5.8100e+05, 2.5260e+05, 5.2610e+04, 5.3000e+03};
  const Double_t preset5[n] = {1., 1.5, 2.5, 2.5, 5., 10., 25.};
  const Double_t preset6[n] = {9.7900e+06, 6.9300e+06, 3.4970e+06, 8.5610e+05, 3.5100e+05, 6.7000e+04, 6.0190e+03};

  for(int i=0;i<n;i++){
    x[i]=preset1[i], y[i]=preset2[i], exl[i]=preset3[i], eyl[i]=preset4[i], exh[i]=preset5[i], eyh[i]=preset6[i];
  }

  TGraphAsymmErrors *graph = new TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);
  graph->SetTitle("FONLL result");
  graph->Draw();
  graph->Write();
  return graph;
}

TGraphAsymmErrors *fonll_pt()
{
  TString outputf;
  outputf = Form("%s",outputfile_pt.Data());
  TFile* outf = new TFile(outputf.Data(),"recreate");
  outf->cd();

  TCanvas *c = new TCanvas("c","FONLL result", 600, 600);

  const Int_t n = nBins;
  Double_t x[n], y[n], exl[n], eyl[n], exh[n], eyh[n];
  
  const Double_t preset1[n] = {6., 8.5, 12.5, 17.5, 25., 40., 75.};
  const Double_t preset2[n] = {1.0331e+07, 5.1920e+06, 1.7012e+06, 5.0562e+05, 1.1445e+05, 12320, 417.09};
  const Double_t preset3[n] = {1., 1.5, 2.5, 2.5, 5., 10., 25.};
  const Double_t preset4[n] = {3.1635e+06, 1.3573e+06, 4.0430e+05, 1.0742e+05, 2.1017e+04, 1842, 49.19};
  const Double_t preset5[n] = {1., 1.5, 2.5, 2.5, 5., 10., 25.};
  const Double_t preset6[n] = {4.5050e+06, 2.1385e+06, 6.2090e+05, 1.5741e+05, 2.8820e+04, 2279, 53.47};

  for(int i=0;i<n;i++){
    x[i]=preset1[i], y[i]=preset2[i], exl[i]=preset3[i], eyl[i]=preset4[i], exh[i]=preset5[i], eyh[i]=preset6[i];
  }
  
  TGraphAsymmErrors *graph = new TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);
  graph->SetTitle("FONLL result");
  graph->Draw();
  graph->Write();
  return graph;
}

int main()
{
  if(!istotal_xsection) fonll_pt();
  else fonll_total();
  return 1;
}
