#include "uti.h"

const int nBins = 10;
const float FF = 0.340; //Fragmentation Fration of b -> B+, to be multiplied for FONLL
double ptBins[nBins+1] = {5., 10., 13., 17., 24., 30., 40., 50., 60., 70., 100.};
double yrange[nBins+1];
TString inputname_ref_13 = "ppref_real_13.txt";
TString inputname_ref_7 = "ppref_real_7.txt";
TString outputname_ref = "ppref_extrapolation.root";
TString inputmcname = "/mnt/T2_US_MIT/submit-hi2/scratch/gwangjun/crab_Bfinder_20190520_Hydjet_Pythia8_BuToJpsiK_1033p1_pt3tkpt0p7dls2_v2_addSamplePthat_pthatweight.root";

using namespace std;

void ppref_real()
{
  gStyle->SetOptStat(0);

  ifstream getdata(inputname_ref_13.Data());
  // dsigmadpt units in [micr/obarn/GeV] 
  double ptmin[9], ptmax[9], dsigmadpt_real_13[9], e1[9], e2[9], e3[9], dsigmadpt_FONLL_13[9], e4[9], e5[9], dsigmadpt_FONLL_5[9], e6[9], e7[9], e_real_13[9];
  int i;
  
  for(i=1;i<10;i++)
    {
      getdata>>ptmin[i];
      getdata>>ptmax[i];
      getdata>>yrange[i];
      getdata>>dsigmadpt_real_13[i];
      getdata>>e1[i];
      getdata>>e2[i];
      getdata>>e3[i];
      getdata>>dsigmadpt_FONLL_13[i];
      getdata>>e4[i];
      getdata>>e5[i];
      getdata>>dsigmadpt_FONLL_5[i];
      getdata>>e6[i];
      getdata>>e7[i];

      e_real_13[i] = sqrt(e1[i]*e1[i]+e2[i]*e2[i]+e3[i]*e3[i]);
    }

  ifstream getdata(inputname_ref_7.Data());
  // dsigmadpt units in [microbarn/GeV] 
  double dsigmadpt_real_7[5], e8[5], e9[5], dsigmadpt_FONLL_7[5], e10[5], e11[5], e_real_7[5];
  
  for(i=0;i<5;i++)
    {
      getdata>>ptmin[i];
      getdata>>ptmax[i];
      getdata>>yrange[i];
      getdata>>dsigmadpt_real_7[i];
      getdata>>e8[i];
      getdata>>e9[i];
      getdata>>dsigmadpt_FONLL_7[i];
      getdata>>e10[i];
      getdata>>e11[i];
      getdata>>dsigmadpt_FONLL_5[i];
      getdata>>e6[i];
      getdata>>e7[i];

      e_real_7[i] = sqrt(e8[i]*e8[i]+e9[i]*e9[i]);
    }
  
  TH1F* d_r = new TH1F("d_r","",nBins,ptBins);
  TH1F* d_f = new TH1F("d_f","",nBins,ptBins);
  TH1F* d_f5 = new TH1F("d_f5","",nBins,ptBins);

  for(i=0;i<nBins;i++)
    {
      d_r->SetBinContent(i+1,dsigmadpt_real[i]);
      d_f->SetBinContent(i+1,dsigmadpt_FONLL[i]);
      d_f5->SetBinContent(i+1,dsigmadpt_FONLL_5p02[i]);
      pythia_ref->SetBinContent(i+1,dsigmadpt_pythia[i]);
      
      d_r->SetBinError(i+1,e_real[i]);
      d_f->SetBinError(i+1,e_FONLL[i]);
      d_f5->SetBinError(i+1,e_FONLL_5p02[i]);
      pythia_ref->SetBinError(i+1,e_pythia[i]);
    }
  
  d_f->Scale(FF);
  d_f5->Scale(FF);
  
  d_r->Sumw2();
  d_f->Sumw2();
  d_f5->Sumw2();
  pythia_ref->Sumw2();
  
  TH1F* real_ref = (TH1F*) d_r->Clone("real_ref");
  real_ref->Sumw2();
  real_ref->Multiply(d_f5);
  real_ref->Divide(d_f);
  
  real_ref->SetTitle("Bpt distribution comparison");
  real_ref->GetXaxis()->SetTitle("p_{t} (GeV/c)");
  real_ref->GetYaxis()->SetTitle("Normalized Distribution");
  real_ref->GetYaxis()->SetTitleOffset(1.0);
  
  real_ref->Scale(1.0/real_ref->Integral("width"));
  pythia_ref->Scale(1.0/pythia_ref->Integral("width"));
  d_f5->Scale(1.0/d_f5->Integral("width"));

  TFile* inputMC = new TFile(inputmcname.Data());
  TTree* Gen = (TTree*) inputMC->Get("Bfinder/ntGen");
  Gen->AddFriend("hiEvtAnalyzer/HiTree");
  TString weightgen = "pthatweight";
  TString Bselection = "(GisSignal==1)&&(GcollisionId==0)&&((Gpt>5&&Gpt<10&&TMath::Abs(Gy)<2.4)||(Gpt>10&&Gpt<17&&TMath::Abs(Gy)<1.45)||(Gpt>17&&Gpt<100&&TMath::Abs(Gy)<2.1))";
  TH1F* Genpt = new TH1F("Genpt","",nBins,ptBins);
  Gen->Project("Genpt","Gpt",TCut(weightgen)*TCut(Bselection));
  //Genpt->Scale(1.0/Genpt->Integral(),"width");
  divideBinWidth(Genpt);
  Genpt->Scale(1.0/Genpt->Integral("width"));
  
  TH1F* Ratiopt = (TH1F*) real_ref->Clone("Ratiopt");
  Ratiopt->Sumw2();
  Ratiopt->Divide(Genpt);
  
  TF1 *f = new TF1("f1","[0]+[1]*x+[2]*x*x+[3]*x*x*x",5,100);
  f->SetParameters(1,0.1,0.1,0.1);
  f->SetParLimits(0,0,10);
  f->SetParLimits(1,-10,10);
  f->SetParLimits(2,-10,10);
  f->SetParLimits(3,-10,10);
  
  Ratiopt->Fit(f,"R");
  
  printf("Bpt weight function(MC Gpt): %f+%f*Bpt+%f*Bpt*Bpt+%f*Bpt*Bpt*Bpt\n",f->GetParameter(0),f->GetParameter(1),f->GetParameter(2),f->GetParameter(3));
  
  Ratiopt->SetTitle("Real_reference/Genpt");
  Ratiopt->GetXaxis()->SetTitle("p_{t} (GeV/c)");
  Ratiopt->GetYaxis()->SetTitle("Real/Gen");
  Ratiopt->GetYaxis()->SetTitleOffset(1.0);
  
  TCanvas* c1 = new TCanvas("","",600,600);
  c1->cd();
  real_ref->SetLineColor(kRed);
  pythia_ref->SetLineColor(kGreen);
  Genpt->SetLineColor(kBlue);
  d_f5->SetLineColor(kMagenta);
  real_ref->Draw();
  pythia_ref->Draw("same");
  Genpt->Draw("same");
  d_f5->Draw("same");

  TLegend *leg = new TLegend(0.50,0.60,0.80,0.75,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->AddEntry(real_ref,"real_reference","l");
  leg->AddEntry(pythia_ref,"pythia_reference","l");
  leg->AddEntry(Genpt,"Genpt_MC","l");
  leg->AddEntry(d_f5,"FONLL_5.02TeV","l");
  leg->Draw("same");

  //c1->RedrawAxis();
  c1->SetLogy();
  c1->SaveAs("Bpt_comparison.png");

  TCanvas* c2 = new TCanvas("","",600,600);
  c2->cd();
  Ratiopt->Draw("");
  f->Draw("same");
  
  TLatex* fitftn = new TLatex(0.15,0.8,Form("Bpt weight: %f+%f*Bpt+%f*Bpt^2+%f*Bpt^3",f->GetParameter(0),f->GetParameter(1),f->GetParameter(2),f->GetParameter(3)));
  fitftn->SetNDC();
  fitftn->SetTextFont(42);
  fitftn->SetTextSize(0.020);
  fitftn->SetLineWidth(1);
  fitftn->Draw("same");

  c2->SaveAs("Bpt_weight.png");  

  TH1F* Ratio_pythia = (TH1F*) real_ref->Clone("Ratio_pythia");
  Ratio_pythia->Sumw2();
  Ratio_pythia->Divide(pythia_ref);
  
  Ratio_pythia->Fit(f,"R");
  
  printf("Bpt weight function(pythia ref): %f+%f*Bpt+%f*Bpt*Bpt+%f*Bpt*Bpt*Bpt\n",f->GetParameter(0),f->GetParameter(1),f->GetParameter(2),f->GetParameter(3));

  Ratio_pythia->SetTitle("Real_reference/pythia_reference");
  Ratio_pythia->SetLineColor(kBlue);
  
  TCanvas* c3 = new TCanvas("","",600,600);
  c3->cd();
  Ratio_pythia->Draw("");

  TLatex* fitftn2 = new TLatex(0.15,0.8,Form("Bpt weight: %f+%f*Bpt+%f*Bpt^2+%f*Bpt^3",f->GetParameter(0),f->GetParameter(1),f->GetParameter(2),f->GetParameter(3)));
  fitftn2->SetNDC();
  fitftn2->SetTextFont(42);
  fitftn2->SetTextSize(0.020);
  fitftn2->SetLineWidth(1);
  fitftn2->Draw("same");

  c3->SaveAs("Bpt_pythia.png");
  
  TFile* output_ref = new TFile(outputname_ref.Data(),"recreate");
  output_ref->cd();
  real_ref->Write();
  pythia_ref->Write();
  Genpt->Write();
  Ratiopt->Write();
  output_ref->Close();
  
  return;
}

int main()
{
  ppref_real();
  return 0;
}
