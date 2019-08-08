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
bool ispp = false;

void fitAll(){
  TFile* outf;
  if(ispp) outf = new TFile(Form("plotspp/fitNP_pp.root"), "recreate");
  else outf = new TFile(Form("plotsPbPb/fitNP_PbPb.root"), "recreate");
  
  gStyle->SetOptStat(0);
  TFile* inf;
  if(ispp) inf = new TFile(infilepp.Data());
  else inf = new TFile(infilePbPb.Data());
  
  TH1D* Bnosig = (TH1D*)inf->Get("Bmass_nosig");
  TH1D* cp1 = (TH1D*)inf->Get("BmassBpPi");
  TH1D* cp2 = (TH1D*)inf->Get("BmassBpK_tkmatch");
  TH1D* cp3 = (TH1D*)inf->Get("BmassB0K_tkmatch");
  TH1D* h = (TH1D*)cp1->Clone();
  h->Add(cp2);
  h->Add(cp3);
  h->GetXaxis()->SetRangeUser(minhisto, maxhisto);
  TH1D* hempty = new TH1D("", "", 50, minhisto, maxhisto);
  TCanvas* c = new TCanvas("c", "", 600, 600);
  c->cd();
  
  //TF1* f = new TF1(Form("f"),"[0]*TMath::Erf((x-[1])/[2]) + [0] + [3]*([4]*TMath::Gaus(x,[5],[6])/(sqrt(2*3.14159)*[6])+(1-[4])*TMath::Gaus(x,[5],[7])/(sqrt(2*3.14159)*[7])) + [9]+[10]*x ");
  TF1* f = new TF1(Form("f"),"[0]*TMath::Erf((x-[1])/[2])+[0]+[3]*([4]*TMath::Gaus(x,[5],[6])/(sqrt(2*3.14159)*[6])+(1-[4])*TMath::Gaus(x,[5],[7])/(sqrt(2*3.14159)*[7]))+[8]*TMath::Gaus(x,[9],[10])/(sqrt(2*3.14159)*[10])");
  
  //error fn (left)
  f->SetParLimits(0,1e0,1e3);
  f->SetParLimits(1,4.9,5.3);
  f->SetParLimits(2,-10,-0.05);
  f->SetParameter(0,20);
  f->SetParameter(1,5.14);
  f->SetParameter(2,-1);
  
  //double gaussian peak (right)
  f->SetParLimits(3,0,1e3);
  f->SetParLimits(4,0,1);
  f->SetParLimits(5,5.30,5.35);
  f->SetParLimits(6,0.005,0.07);
  f->SetParLimits(7,0.005,0.3);
  f->SetParameter(3,10);
  f->SetParameter(4,0.5);
  f->SetParameter(5,5.35);
  f->SetParameter(6,0.05);
  f->SetParameter(7,0.3);

  //Gaussian peak (left)
  f->SetParLimits(8,0,1e2);
  f->SetParLimits(9,5.10,5.12);
  f->SetParLimits(10,0.01,0.05);
  f->SetParameter(8,10);
  f->SetParameter(9,5.10);
  f->SetParameter(10,0.05);
  
  /*
  //Combinatorial Background (linear)
  f->SetParLimits(9, 0, 1e5);
  f->SetParLimits(10, -500,  100);
  f->SetParameter(9,1e3);
  f->SetParameter(10,-1);
  */

  h->Fit(Form("f"),"q","",minhisto,maxhisto);
  h->Fit(Form("f"),"q","",minhisto,maxhisto);
  h->Fit(Form("f"),"L q","",minhisto,maxhisto);
  h->Fit(Form("f"),"L q","",minhisto,maxhisto);
  h->Fit(Form("f"),"L q","",minhisto,maxhisto);
  h->Fit(Form("f"),"L m","",minhisto,maxhisto);
  h->SetMarkerSize(0.8);
  h->SetMarkerStyle(20);
  
  printf("%f*TMath::Erf((x-%f)/%f)+%f+%f*(%f*TMath::Gaus(x,%f,%f)/(sqrt(2*3.14159)*%f)+(1-%f)*TMath::Gaus(x,%f,%f)/(sqrt(2*3.14159)*%f))+%f*TMath::Gaus(x,%f,%f)/(sqrt(2*3.14159)*%f) \n", f->GetParameter(0), f->GetParameter(1), f->GetParameter(2), f->GetParameter(0), f->GetParameter(3), f->GetParameter(4), f->GetParameter(5), f->GetParameter(6), f->GetParameter(6), f->GetParameter(4), f->GetParameter(5), f->GetParameter(7), f->GetParameter(7), f->GetParameter(8), f->GetParameter(9), f->GetParameter(10), f->GetParameter(10));

  std::cout<<"chisq = "<<f->GetChisquare()/f->GetNDF()<<std::endl;
  
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
  hempty->SetStats(0);
  
  hempty->SetMaximum(h->GetMaximum()*1.2);
  //hempty->SetMaximum(Bnosig->GetMaximum()*1.2);
  hempty->Draw();
  Bnosig->Draw("same e");
  h->Draw("same e");
  f->Draw("same");
  
  TLegend* leg = new TLegend(0.55,0.52,0.70,0.65,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);  
  leg->AddEntry(h,"Non-prompt events","pl");
  //leg->AddEntry(h,"Bs #rightarrow J/#psi + X","pl");                                                                                   
  //leg->AddEntry(hempty,"B^{+} #rightarrow J/#psi + #pi","pl");
  //leg->AddEntry(hempty,"B^{+} #rightarrow J/#psi + various K","pl");
  //leg->AddEntry(hempty,"B^{0} #rightarrow J/#psi + various K","pl");
  leg->Draw("same");
  
  TLatex* texChi = new TLatex(0.55,0.475,Form("#chi^{2}/nDOF:%.2f/%d=%.2f", f->GetChisquare(), f->GetNDF(), f->GetChisquare()/f->GetNDF()));
  texChi->SetNDC();
  texChi->SetTextAlign(12);
  texChi->SetTextSize(0.035);
  texChi->SetTextFont(42);
  texChi->Draw();
  
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
  if(ispp) c->SaveAs(Form("plotspp/fitNP_pp.png"));
  else c->SaveAs(Form("plotsPbPb/fitNP_PbPb.png"));
}

int main()
{
  fitAll();
  return 0;
}
