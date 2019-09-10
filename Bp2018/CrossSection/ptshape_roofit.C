#include "uti.h"
#include "parameters.h"

const int iteration = 0;
bool dofit = 1;
bool examinevar = 1;

TString inputyieldname_data = "ptshape/BDT/yields_Bp_binned_pt.root";
TString inputyieldname_MC = Form("ptshape/BDT/yields_Bp_mc_binned_pt_iter%d.root",iteration);
TString outputBptweight = Form("ptshape/BDT/Bptweight_roofit_iter%d.root",iteration);

void ptshape_roofit()
{
  gStyle->SetOptStat(0);

  TFile* inputyield_data = new TFile(inputyieldname_data.Data());
  TFile* inputyield_MC = new TFile(inputyieldname_MC.Data());
  TH1D* yield_data = (TH1D*)inputyield_data->Get("hPt");
  yield_data->GetYaxis()->SetTitle("Normalized Raw Yield");
  yield_data->GetXaxis()->SetTitle("p_{t} (GeV)");
  yield_data->GetYaxis()->SetTitleOffset(1.1);
  yield_data->SetLineColor(kRed);
  //yield_data->SetMarkerSize(1);
  TH1D* yield_MC = (TH1D*)inputyield_MC->Get("hPtMC");
  yield_MC->GetYaxis()->SetTitle("Normalized Raw Yield");
  yield_MC->GetXaxis()->SetTitle("p_{t} (GeV)");
  yield_MC->GetYaxis()->SetTitleOffset(1.1);
  yield_MC->SetLineColor(kGreen);
  yield_MC->SetMarkerSize(0);

  double ptwidth[nBins];
  for(int i=0;i<nBins;i++)
    {
      ptwidth[i]=ptBins[i+1]-ptBins[i];
      yield_data->SetBinContent(i+1,yield_data->GetBinContent(i+1)*ptwidth[i]);
      yield_data->SetBinError(i+1,yield_data->GetBinError(i+1)*ptwidth[i]);
      yield_MC->SetBinContent(i+1,yield_MC->GetBinContent(i+1)*ptwidth[i]);
      yield_MC->SetBinError(i+1,yield_MC->GetBinError(i+1)*ptwidth[i]);
    }

  yield_data->Scale(1.0/yield_data->Integral());
  yield_MC->Scale(1.0/yield_MC->Integral());
  
  TCanvas* Cyield = new TCanvas("","",600,600);
  Cyield->cd();
  yield_data->Draw("ep");
  yield_MC->Draw("ep same");
  
  TLegend *legyield = new TLegend(0.45,0.70,0.75,0.80,NULL,"brNDC");
  legyield->SetBorderSize(0);
  legyield->SetTextSize(0.04);
  legyield->SetTextFont(42);
  legyield->SetFillStyle(0);
  legyield->AddEntry(yield_data,"Data Raw Yield","l");
  legyield->AddEntry(yield_MC,"MC Raw Yield","l");
  legyield->Draw("same");

  Cyield->SaveAs(Form("ptshape/ptshape_roofit_yield_iter%d.png",iteration));
  Cyield->SaveAs(Form("ptshape/ptshape_roofit_yield_iter%d.pdf",iteration));

  TH1D* Ratio = (TH1D*) yield_data->Clone("Ratio");
  Ratio->Divide(yield_MC);
  Ratio->GetYaxis()->SetTitle("Data/MC Raw Yield");
  Ratio->SetLineColor(kBlack);

  TF1* f = new TF1("f","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",5,60);
  f->SetParameters(1,0,0,0,0);
  f->SetParLimits(0,-10,10);
  f->SetParLimits(1,-10,10);
  f->SetParLimits(2,-10,10);
  f->SetParLimits(3,-10,10);
  f->SetParLimits(4,-10,10);
  f->SetLineColor(kRed);
  if(dofit) Ratio->Fit(f,"I");

  if(dofit) printf("Bpt weight(iteration %d): (%.8f+%.8f*Gpt+%.8f*Gpt*Gpt+%.8f*Gpt*Gpt*Gpt+%.8f*Gpt*Gpt*Gpt*Gpt)\n",iteration,f->GetParameter(0),f->GetParameter(1),f->GetParameter(2),f->GetParameter(3),f->GetParameter(4));
  if(dofit) printf("Bpt weight(iteration %d): (%.8f+%.8f*Bgenpt+%.8f*Bgenpt*Bgenpt+%.8f*Bgenpt*Bgenpt*Bgenpt+%.8f*Bgenpt*Bgenpt*Bgenpt*Bgenpt)\n",iteration,f->GetParameter(0),f->GetParameter(1),f->GetParameter(2),f->GetParameter(3),f->GetParameter(4));
  
  TH1D* Ratio_plus = new TH1D("Ratio_plus","",nBins,ptBins);
  TH1D* Ratio_minus = new TH1D("Ratio_minus","",nBins,ptBins);
  TF1* f_plus = new TF1("f_plus","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",5,60);
  TF1* f_minus = new TF1("f_minus","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",5,60);

  if(examinevar)
    {
      for(int i=0;i<nBins;i++)
	{
	  Ratio_plus->SetBinContent(i+1,Ratio->GetBinContent(i+1)+Ratio->GetBinError(i+1));
	  Ratio_plus->SetBinError(i+1,Ratio->GetBinError(i+1));
	  Ratio_minus->SetBinContent(i+1,Ratio->GetBinContent(i+1)-Ratio->GetBinError(i+1));
	  Ratio_minus->SetBinError(i+1,Ratio->GetBinError(i+1));
	}
      
      f_plus->SetParameters(1,0,0,0,0);
      f_plus->SetParLimits(0,-10,10);
      f_plus->SetParLimits(1,-10,10);
      f_plus->SetParLimits(2,-10,10);
      f_plus->SetParLimits(3,-10,10);
      f_plus->SetParLimits(4,-10,10);
      f_plus->SetLineColor(kBlue);
      if(dofit) Ratio_plus->Fit(f_plus,"I");
      
      if(dofit) printf("Bpt weight +1 sigma(iteration %d): (%.8f+%.8f*Gpt+%.8f*Gpt*Gpt+%.8f*Gpt*Gpt*Gpt+%.8f*Gpt*Gpt*Gpt*Gpt)\n",iteration,f_plus->GetParameter(0),f_plus->GetParameter(1),f_plus->GetParameter(2),f_plus->GetParameter(3),f_plus->GetParameter(4));  
      if(dofit) printf("Bpt weight +1 sigma(iteration %d): (%.8f+%.8f*Bgenpt+%.8f*Bgenpt*Bgenpt+%.8f*Bgenpt*Bgenpt*Bgenpt+%.8f*Bgenpt*Bgenpt*Bgenpt*Bgenpt)\n",iteration,f_plus->GetParameter(0),f_plus->GetParameter(1),f_plus->GetParameter(2),f_plus->GetParameter(3),f_plus->GetParameter(4));

      f_minus->SetParameters(1,0,0,0,0);
      f_minus->SetParLimits(0,-10,10);
      f_minus->SetParLimits(1,-10,10);
      f_minus->SetParLimits(2,-10,10);
      f_minus->SetParLimits(3,-10,10);
      f_minus->SetParLimits(4,-10,10);
      f_minus->SetLineColor(kGreen);
      if(dofit) Ratio_minus->Fit(f_minus,"I");
      
      if(dofit) printf("Bpt weight -1 sigma(iteration %d): (%.8f+%.8f*Gpt+%.8f*Gpt*Gpt+%.8f*Gpt*Gpt*Gpt+%.8f*Gpt*Gpt*Gpt*Gpt)\n",iteration,f_minus->GetParameter(0),f_minus->GetParameter(1),f_minus->GetParameter(2),f_minus->GetParameter(3),f_minus->GetParameter(4));
      if(dofit) printf("Bpt weight -1 sigma(iteration %d): (%.8f+%.8f*Bgenpt+%.8f*Bgenpt*Bgenpt+%.8f*Bgenpt*Bgenpt*Bgenpt+%.8f*Bgenpt*Bgenpt*Bgenpt*Bgenpt)\n",iteration,f_minus->GetParameter(0),f_minus->GetParameter(1),f_minus->GetParameter(2),f_minus->GetParameter(3),f_minus->GetParameter(4));
    }

  TCanvas* CRatio = new TCanvas("","",600,600);
  CRatio->cd();
  Ratio->Draw("ep");
  f->Draw("same");
  if(examinevar)
    {
      f_plus->Draw("same");
      f_minus->Draw("same");
    }
  
  TLegend *legRatio = new TLegend(0.45,0.70,0.75,0.80,NULL,"brNDC");
  legRatio->SetBorderSize(0);
  legRatio->SetTextSize(0.04);
  legRatio->SetTextFont(42);
  legRatio->SetFillStyle(0);
  legRatio->AddEntry(f,"Data/MC","l");
  if(examinevar)
    {
      legRatio->AddEntry(f_plus,"Data/MC +1 #sigma","l");
      legRatio->AddEntry(f_minus,"Data/MC -1 #sigma","l");
    }
  legRatio->Draw("same");

  CRatio->SaveAs(Form("ptshape/ptshape_roofit_ratio_iter%d.png",iteration));
  CRatio->SaveAs(Form("ptshape/ptshape_roofit_ratio_iter%d.pdf",iteration));

  TFile* output = new TFile(outputBptweight.Data(),"recreate");
  output->cd();
  yield_data->Write();
  yield_MC->Write();
  Ratio->Write();
  if(examinevar)
    {
      Ratio_plus->Write();
      Ratio_minus->Write();
    }
  output->Close();

  return;
}

int main()
{
  ptshape_roofit();
  return 0;
}
