#include <TCanvas.h>
#include <TPad.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TBox.h>
#include <TCut.h>
#include <TMath.h>
#include <TStyle.h>
#include <TF1.h>
#include <TH1.h>
#include <TFile.h>

double minhisto = 5.0;
double maxhisto = 6.0;

double ptmin = 5.0;
double ptmax = 100.0;
double centMin = 0.0;
double centMax = 90.0;

TString infilepp = "";
TString infilePbPb = "NP_nominal.root";

TString histname;

void fitIndi(){
  void fit(int fitType, bool ispp);
  for(int i = 1; i < 5; i++){
    //fit(i, true);
    fit(i, false);
  }
}

void fit(int fitType, bool ispp)
{
  
  if(fitType==1) histname = "BmassBpPi";
  if(fitType==2) histname = "BmassBpK_tkmatch";
  if(fitType==3) histname = "BmassB0K_tkmatch";
  if(fitType==4) histname = "Bmass_nominal";

  TFile* outf;
  if(ispp) outf = new TFile(Form("plotspp/%s.root", histname.Data()), "recreate");
  else outf = new TFile(Form("plotsPbPb/%s.root", histname.Data()), "recreate");
  
  gStyle->SetOptStat(0);
  TFile* inf;
  if(ispp) inf = new TFile(infilepp.Data());
  else inf = new TFile(infilePbPb.Data());
  
  TH1D* h = (TH1D*)inf->Get(histname.Data());
  h->GetXaxis()->SetRangeUser(minhisto, maxhisto);
  TH1D* hempty = new TH1D("", "", 50, minhisto, maxhisto);
  TCanvas* c = new TCanvas("c", "", 600, 600);
  c->cd();
  
  TF1* f = new TF1(Form("f"),"[0]*TMath::Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2]) + [3]*TMath::Gaus(x,[4],[5])/(sqrt(2*3.14159)*[5]) + [9]+[10]*x ");
  f->SetParLimits(0, 1e-3, 1e3);
  f->SetParLimits(1, 4.9, 5.2);
  f->SetParLimits(2, 0.001, 0.1);
  
  f->SetParLimits(3, 1e-3, 1e3);
  f->SetParLimits(4, 4.9, 5.2);
  f->SetParLimits(5, 0.001, 0.1);
  
  f->SetParLimits(9, 0, 1e2);
  f->SetParLimits(10, -10,  10);
  
  f->SetParameter(0,1e2);
  f->SetParameter(1,5.11);
  f->SetParameter(2,0.05);
  
  f->SetParameter(3,1e2);
  f->SetParameter(4,4.98);
  f->SetParameter(5,0.05);

  f->SetParameter(9,1e3);
  f->SetParameter(10,0);
  
  if(fitType==1) {
    f->SetParLimits(0, 0, 1e4);
    f->SetParLimits(1, 5.2, 5.6);
    f->SetParLimits(2, 0.0001, 0.5);
    f->SetParameter(0,1e2);
    f->SetParameter(1,5.35);
    f->SetParameter(2,0.05);
    
    f->SetParLimits(3, 0, 1e4);
    f->SetParLimits(4, 5.2, 5.6);
    f->SetParLimits(5, 0.0001, 0.5);
    f->SetParameter(3,1e2);
    f->SetParameter(4,5.35);
    f->SetParameter(5,0.05);
  }
  
  //h->Fit(Form("f"),"q","",minhisto,maxhisto);
  //h->Fit(Form("f"),"q","",minhisto,maxhisto);
  //h->Fit(Form("f"),"L q","",minhisto,maxhisto);
  //h->Fit(Form("f"),"L m","",minhisto,maxhisto);
  h->SetMarkerSize(0.8);
  h->SetMarkerStyle(20);
  
  hempty->SetMaximum(h->GetMaximum()*1.4);
  hempty->SetXTitle("m_{#mu#muK} (GeV/c^{2})");
  hempty->SetYTitle("Entries / (20 MeV/c^{2})");
  hempty->GetXaxis()->CenterTitle();
  hempty->GetYaxis()->CenterTitle();
  hempty->SetAxisRange(0,h->GetMaximum()*1.4*1.2,"Y");
  hempty->GetXaxis()->SetTitleOffset(1.0);
  hempty->GetYaxis()->SetTitleOffset(1.1);
  hempty->GetXaxis()->SetLabelOffset(0.007);
  hempty->GetYaxis()->SetLabelOffset(0.007);
  hempty->GetXaxis()->SetTitleSize(0.045);
  hempty->GetYaxis()->SetTitleSize(0.045);
  hempty->GetXaxis()->SetTitleFont(42);
  hempty->GetYaxis()->SetTitleFont(42);
  hempty->GetXaxis()->SetLabelFont(42);
  hempty->GetYaxis()->SetLabelFont(42);
  hempty->GetXaxis()->SetLabelSize(0.04);
  hempty->GetYaxis()->SetLabelSize(0.04);
  hempty->SetMarkerSize(0.01);
  hempty->SetMarkerStyle(20);
  
  hempty->Draw();
  h->Draw("same e");
  //f->Draw("same");
  
  TLegend* leg = new TLegend(0.55,0.51,0.75,0.71,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  if(fitType==1) leg->AddEntry(h,"B^{+} #rightarrow J/#psi + #pi","pl");
  if(fitType==2) leg->AddEntry(h,"B^{+} #rightarrow J/#psi + K/K^{*}/#phi","pl");
  if(fitType==3) leg->AddEntry(h,"B^{0} #rightarrow J/#psi + K/K^{*}/#phi","pl");
  if(fitType==4) leg->AddEntry(h,"B^{+} #rightarrow J/#psi + K^{+}","pl");
  leg->Draw("same");
  
  TLatex* texChi = new TLatex(0.55,0.475,Form("#chi^{2}/nDOF:%.2f/%d=%.2f", f->GetChisquare(), f->GetNDF(), f->GetChisquare()/f->GetNDF()));
  texChi->SetNDC();
  texChi->SetTextAlign(12);
  texChi->SetTextSize(0.035);
  texChi->SetTextFont(42);
  //texChi->Draw();
  
  TLatex* texCms = new TLatex(0.18,0.93, "#scale[1.25]{CMS} Preliminary");
  texCms->SetNDC();
  texCms->SetTextAlign(12);
  texCms->SetTextSize(0.04);
  texCms->SetTextFont(42);
  texCms->Draw();
  
  TLatex* texCol;
  if(ispp) texCol= new TLatex(0.96,0.93, Form("%s #sqrt{s_{NN}} = 5.02 TeV","pp"));
  else texCol= new TLatex(0.96,0.93, Form("%s #sqrt{s_{NN}} = 5.02 TeV","PbPb"));
  texCol->SetNDC();
  texCol->SetTextAlign(32);
  texCol->SetTextSize(0.04);
  texCol->SetTextFont(42);
  texCol->Draw();
  
  TLatex* tex;
  tex = new TLatex(0.22,0.78,Form("%.0f < p_{T} < %.0f GeV/c",ptmin,ptmax));
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();
  
  TString texper="%";
  tex = new TLatex(0.22,0.71,Form("Centrality %.0f-%.0f%s",centMin,centMax,texper.Data()));//0.2612903,0.8425793                          
  tex->SetNDC();
  tex->SetTextColor(1);
  tex->SetTextFont(42);
  tex->SetTextSize(0.045);
  tex->SetLineWidth(2);
  tex->Draw();
  
  tex = new TLatex(0.22,0.83,"|y| < 2.4");
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();
  
  outf->cd();
  h->Write();
  f->Write();
  if(ispp) c->SaveAs(Form("plotspp/%s.pdf", histname.Data()));
  else c->SaveAs(Form("plotsPbPb/%s.pdf", histname.Data()));
}

int main()
{
  fitIndi();
  return 0;
}
