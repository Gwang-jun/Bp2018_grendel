#include "uti.h"
#include <TMultiGraph.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>

const float FF = 0.405; //Fragmentation Fration of b -> B+, to be multiplied for FONLL

const int nBins = 10;
double ptBins[nBins+1] = {5., 10., 13., 17., 24., 30., 40., 50., 60., 70., 100.};
double ptcent[nBins] = {7.5, 11.5, 15., 20.5, 27., 35., 45., 55., 65., 85.};
double ptwidth[nBins] = {2.5, 1.5, 2., 3.5, 3., 5., 5., 5., 5., 15.};

double ptcent_13[9] = {11.5, 15., 20.5, 27., 35., 45., 55., 65., 85.};
double ptwidth_13[9] = {1.5, 2., 3.5, 3., 5., 5., 5., 5., 15.};

double ptcent_7[5] = {7.5, 11.5, 15., 20.5, 27.};
double ptwidth_7[5] = {2.5, 1.5, 2., 3.5, 3.};

const int nBins_new = 7;
double ptBins_new[nBins_new+1] = {5., 7., 10., 15., 20., 30., 50., 100.};
double ptcent_new[nBins_new] = {6., 8.5, 12.5, 17.5, 25., 40., 75.};
double ptwidth_new[nBins_new] = {1., 1.5, 2.5, 2.5, 5., 10., 25.};

const int nBins_2015 = 5;
double ptBins_2015[nBins_2015+1] = {7., 10., 15., 20., 30., 50.};

TString inputname_ref_13 = "ppref_real_13.txt";
TString inputname_ref_7 = "ppref_real_7.txt";
TString outputname_ref = "ppref_extrapolation.root";
TString inputmcname = "/raid5/data/gwangjun/crab_Bfinder_20190624_Hydjet_Pythia8_Official_BuToJpsiK_1033p1_pt3tkpt0p7dls2_allpthat_pthatweight.root";
//TString inputmcname = "/mnt/T2_US_MIT/submit-hi2/scratch/gwangjun/crab_Bfinder_20190624_Hydjet_Pythia8_Official_BuToJpsiK_1033p1_pt3tkpt0p7dls2_allpthat_pthatweight.root";
//TString inputmcname = "/export/d00/scratch/zzshi/CMSSW_7_5_8_patch3/Merge/2018Ana/Samples/FinalAnaSamples/PrivateMC-Data-Official/MC_Bs_PbPb_TMVA_BDT_PbPb.root";


using namespace std;

void ppref_real()
{
  gStyle->SetOptStat(0);

  gStyle->SetTextSize(0.05);
  gStyle->SetTextFont(42);
  gStyle->SetPadRightMargin(0.043);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.145);
  gStyle->SetTitleX(.0f);

  ifstream getdata_13(inputname_ref_13.Data());
  // dsigmadpt units in [picobarn/GeV] 
  double ptmin_13[9], ptmax_13[9], yrange_13[9], dsigmadpt_real_13[9], e1[9], e2[9], e3[9], dsigmadpt_FONLL_13[9], e4[9], e5[9], dsigmadpt_FONLL_5for13[9], e6[9], e7[9], e_real_13[9], ptcent_13[9], ptwidth_13[9];
  int i;
  
  for(i=0;i<9;i++)
    {
      getdata_13>>ptmin_13[i];
      getdata_13>>ptmax_13[i];
      getdata_13>>yrange_13[i];
      getdata_13>>dsigmadpt_real_13[i], dsigmadpt_real_13[i] *= 1e+06;
      getdata_13>>e1[i], e1[i] *= 1e+06;
      getdata_13>>e2[i], e2[i] *= 1e+06;
      getdata_13>>e3[i], e3[i] *= 1e+06;
      getdata_13>>dsigmadpt_FONLL_13[i], dsigmadpt_FONLL_13[i] *= 1e+06;
      getdata_13>>e4[i], e4[i] *= 1e+06;
      getdata_13>>e5[i], e5[i] *= 1e+06;
      getdata_13>>dsigmadpt_FONLL_5for13[i], dsigmadpt_FONLL_5for13[i] *= 1e+06;
      getdata_13>>e6[i], e6[i] *= 1e+06;
      getdata_13>>e7[i], e7[i] *= 1e+06;

      e_real_13[i] = sqrt(e1[i]*e1[i]+e2[i]*e2[i]+e3[i]*e3[i]);
      dsigmadpt_FONLL_13[i] *= FF;
      e4[i] *= FF;
      e5[i] *= FF;
      dsigmadpt_FONLL_5for13[i] *= FF;
      e6[i] *= FF;
      e7[i] *= FF;
    }

  ifstream getdata_7(inputname_ref_7.Data());
  // dsigmadpt units in [picobarn/GeV] 
  double ptmin_7[5], ptmax_7[5], yrange_7[5], dsigmadpt_real_7[5], e8[5], e9[5], dsigmadpt_FONLL_7[5], e10[5], e11[5], dsigmadpt_FONLL_5for7[5], e12[5], e13[5], e_real_7[5], ptcent_7[5], ptwidth_7[5];
  
  for(i=0;i<5;i++)
    {
      getdata_7>>ptmin_7[i];
      getdata_7>>ptmax_7[i];
      getdata_7>>yrange_7[i];
      getdata_7>>dsigmadpt_real_7[i], dsigmadpt_real_7[i] *= 1e+06;
      getdata_7>>e8[i], e8[i] *= 1e+06;
      getdata_7>>e9[i], e9[i] *= 1e+06;
      getdata_7>>dsigmadpt_FONLL_7[i], dsigmadpt_FONLL_7[i] *= 1e+06;
      getdata_7>>e10[i], e10[i] *= 1e+06;
      getdata_7>>e11[i], e11[i] *= 1e+06;
      getdata_7>>dsigmadpt_FONLL_5for7[i], dsigmadpt_FONLL_5for7[i] *= 1e+06;
      getdata_7>>e12[i], e12[i] *= 1e+06;
      getdata_7>>e13[i], e13[i] *= 1e+06;

      e_real_7[i] = sqrt(e8[i]*e8[i]+e9[i]*e9[i]);    
      dsigmadpt_FONLL_7[i] *= FF;
      e10[i] *= FF;
      e11[i] *= FF;
      dsigmadpt_FONLL_5for7[i] *= FF;
      e12[i] *= FF;
      e13[i] *= FF;
    }

  double y_13[9], y_13_l[9], y_13_h[9], y_13_l_d[9], y_13_l_f[9], y_13_h_d[9], y_13_h_f[9];
  double y_7[5], y_7_l[5], y_7_h[5], y_7_l_d[5], y_7_l_f[5], y_7_h_d[5], y_7_h_f[5];
  double y[nBins], y_l[nBins], y_h[nBins], y_l_d[nBins], y_l_f[nBins], y_h_d[nBins], y_h_f[nBins];
  double y_new[nBins_new], y_new_l[nBins_new], y_new_h[nBins_new], y_new_l_d[nBins_new], y_new_l_f[nBins_new], y_new_h_d[nBins_new], y_new_h_f[nBins_new];
  double p, q, r, s;

  for(i=0;i<9;i++)
    {
      p = dsigmadpt_FONLL_5for13[i]/dsigmadpt_FONLL_13[i];
      q = (e7[i]/e5[i]-p)/p;
      r = (e6[i]/e4[i]-p)/p;
      s = e_real_13[i]/dsigmadpt_real_13[i];

      y_13[i] = dsigmadpt_real_13[i]*p;
      y_13_l[i] = y_13[i]*sqrt(s*s+r*r), y_13_l_d[i] = y_13[i]*s, y_13_l_f[i] = y_13[i]*r;
      y_13_h[i] = y_13[i]*sqrt(s*s+q*q), y_13_h_d[i] = y_13[i]*s, y_13_h_f[i] = y_13[i]*q;
    }

  for(i=0;i<5;i++)
    {
      p = dsigmadpt_FONLL_5for7[i]/dsigmadpt_FONLL_7[i];
      q = (e13[i]/e11[i]-p)/p;
      r = (e12[i]/e10[i]-p)/p;
      s = e_real_7[i]/dsigmadpt_real_7[i];

      y_7[i] = dsigmadpt_real_7[i]*p;
      y_7_l[i] = y_7[i]*sqrt(s*s+r*r), y_7_l_d[i] = y_7[i]*s, y_7_l_f[i] = y_7[i]*r;
      y_7_h[i] = y_7[i]*sqrt(s*s+q*q), y_7_h_d[i] = y_7[i]*s, y_7_h_f[i] = y_7[i]*q;
    }

  //printf("\n");
  //printf("x-section and error before rebinning\n");

  y[0] = y_7[0];
  y_l[0] = y_7_l[0], y_l_d[0] = y_7_l_d[0], y_l_f[0] = y_7_l_f[0];
  y_h[0] = y_7_h[0], y_h_d[0] = y_7_h_d[0], y_h_f[0] = y_7_h_f[0];
  //printf("pt %.0f-%.0f: %e (central) %e (low) %e (high)\n",ptBins[0],ptBins[1],y[0],y_l[0],y_h[0]);
  //printf("pt %.0f-%.0f: %e (low_data) %e (low_FONLL) %e (high_data) %e (high_FONLL)\n",ptBins[0],ptBins[1],y_l_d[0],y_l_f[0],y_h_d[0],y_h_f[0]);

  for(i=1;i<nBins;i++)
    {
      y[i] = y_13[i-1];
      y_l[i] = y_13_l[i-1], y_l_d[i] = y_13_l_d[i-1], y_l_f[i] = y_13_l_f[i-1];
      y_h[i] = y_13_h[i-1], y_h_d[i] = y_13_h_d[i-1], y_h_f[i] = y_13_h_f[i-1];
      //printf("pt %.0f-%.0f: %e (central) %e (low) %e (high)\n",ptBins[i],ptBins[i+1],y[i],y_l[i],y_h[i]);
      //printf("pt %.0f-%.0f: %e (low_data) %e (low_FONLL) %e (high_data) %e (high_FONLL)\n",ptBins[i],ptBins[i+1],y_l_d[i],y_l_f[i],y_h_d[i],y_h_f[i]);
    }
  

  TH1D* yhis = new TH1D("yhis","",nBins,ptBins);
  yhis->GetXaxis()->SetTitle("p_{t} (GeV)");
  yhis->GetYaxis()->SetTitle("d#sigma/dp_{t} (pb/GeV)");

  for(i=0;i<nBins;i++)
    {
      yhis->SetBinContent(i+1,y[i]);
      yhis->SetBinError(i+1,0.5*(y_l[i]+y_h[i]));
    }
  
  //TF1* g = new TF1("g","[0]+TMath::Exp([1]+[2]*x+[3]/x+[4]/(x*x))",5,100);
  TF1* g = new TF1("g","TMath::Exp([0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x)",5,100);
  g->SetParLimits(0,-100,100);
  g->SetParLimits(1,-100,100);
  g->SetParLimits(2,-100,100);
  g->SetParLimits(3,-100,100);
  g->SetParLimits(4,-100,100);
  yhis->Fit(g,"IR");

  TCanvas* c0 = new TCanvas("","",600,600);
  c0->cd();
  c0->SetLogy();
  yhis->Draw();
  g->Draw("same");
  c0->SaveAs("pp_continuous.png");
  c0->SaveAs("pp_continuous.pdf");
  
  
  // REBINNING, be careful
  
  double x, dxl, dxh, dxl_d, dxl_f, dxh_d, dxh_f;

  x = g->Integral(5.,7.)/2.;
  dxl = x*(g->Integral(5.,7.)/g->Integral(5.,7.)*y_l[0]/y[0]);
  dxl_d =x*(g->Integral(5.,7.)/g->Integral(5.,7.)*y_l_d[0]/y[0]);
  dxl_f =x*(g->Integral(5.,7.)/g->Integral(5.,7.)*y_l_f[0]/y[0]);
  dxh = x*(g->Integral(5.,7.)/g->Integral(5.,7.)*y_h[0]/y[0]);
  dxh_d =x*(g->Integral(5.,7.)/g->Integral(5.,7.)*y_h_d[0]/y[0]);
  dxh_f =x*(g->Integral(5.,7.)/g->Integral(5.,7.)*y_h_f[0]/y[0]);
  y_new[0] = x, y_new_l[0] = dxl, y_new_h[0] = dxh, y_new_l_d[0] = dxl_d, y_new_h_d[0] = dxh_d, y_new_l_f[0] = dxl_f, y_new_h_f[0] = dxh_f;

  x = g->Integral(7.,10.)/3.;
  dxl = x*(g->Integral(7.,10.)/g->Integral(7.,10.)*y_l[0]/y[0]);
  dxl_d =x*(g->Integral(7.,10.)/g->Integral(7.,10.)*y_l_d[0]/y[0]);
  dxl_f =x*(g->Integral(7.,10.)/g->Integral(7.,10.)*y_l_f[0]/y[0]);
  dxh = x*(g->Integral(7.,10.)/g->Integral(7.,10.)*y_h[0]/y[0]);
  dxh_d =x*(g->Integral(7.,10.)/g->Integral(7.,10.)*y_h_d[0]/y[0]);
  dxh_f =x*(g->Integral(7.,10.)/g->Integral(7.,10.)*y_h_f[0]/y[0]);
  y_new[1] = x, y_new_l[1] = dxl, y_new_h[1] = dxh, y_new_l_d[1] = dxl_d, y_new_h_d[1] = dxh_d, y_new_l_f[1] = dxl_f, y_new_h_f[1] = dxh_f;

  x = g->Integral(10.,15.)/5.;
  dxl = x*(g->Integral(10.,13.)/g->Integral(10.,15.)*y_l[1]/y[1]+g->Integral(13.,15.)/g->Integral(10.,15.)*y_l[2]/y[2]);
  dxl_d = x*(g->Integral(10.,13.)/g->Integral(10.,15.)*y_l_d[1]/y[1]+g->Integral(13.,15.)/g->Integral(10.,15.)*y_l_d[2]/y[2]);
  dxl_f = x*(g->Integral(10.,13.)/g->Integral(10.,15.)*y_l_f[1]/y[1]+g->Integral(13.,15.)/g->Integral(10.,15.)*y_l_f[2]/y[2]);
  dxh = x*(g->Integral(10.,13.)/g->Integral(10.,15.)*y_h[1]/y[1]+g->Integral(13.,15.)/g->Integral(10.,15.)*y_h[2]/y[2]);
  dxh_d = x*(g->Integral(10.,13.)/g->Integral(10.,15.)*y_h_d[1]/y[1]+g->Integral(13.,15.)/g->Integral(10.,15.)*y_h_d[2]/y[2]);
  dxh_f = x*(g->Integral(10.,13.)/g->Integral(10.,15.)*y_h_f[1]/y[1]+g->Integral(13.,15.)/g->Integral(10.,15.)*y_h_f[2]/y[2]);
  y_new[2] = x, y_new_l[2] = dxl, y_new_h[2] = dxh, y_new_l_d[2] = dxl_d, y_new_h_d[2] = dxh_d, y_new_l_f[2] = dxl_f, y_new_h_f[2] = dxh_f;

  x = g->Integral(15.,20.)/5.;
  dxl = x*(g->Integral(15.,17.)/g->Integral(15.,20.)*y_l[2]/y[2]+g->Integral(17.,20.)/g->Integral(15.,20.)*y_l[3]/y[3]);
  dxl_d = x*(g->Integral(15.,17.)/g->Integral(15.,20.)*y_l_d[2]/y[2]+g->Integral(17.,20.)/g->Integral(15.,20.)*y_l_d[3]/y[3]);
  dxl_f = x*(g->Integral(15.,17.)/g->Integral(15.,20.)*y_l_f[2]/y[2]+g->Integral(17.,20.)/g->Integral(15.,20.)*y_l_f[3]/y[3]);
  dxh = x*(g->Integral(15.,17.)/g->Integral(15.,20.)*y_h[2]/y[2]+g->Integral(17.,20.)/g->Integral(15.,20.)*y_h[3]/y[3]);
  dxh_d = x*(g->Integral(15.,17.)/g->Integral(15.,20.)*y_h_d[2]/y[2]+g->Integral(17.,20.)/g->Integral(15.,20.)*y_h_d[3]/y[3]);
  dxh_f = x*(g->Integral(15.,17.)/g->Integral(15.,20.)*y_h_f[2]/y[2]+g->Integral(17.,20.)/g->Integral(15.,20.)*y_h_f[3]/y[3]);
  y_new[3] = x, y_new_l[3] = dxl, y_new_h[3] = dxh, y_new_l_d[3] = dxl_d, y_new_h_d[3] = dxh_d, y_new_l_f[3] = dxl_f, y_new_h_f[3] = dxh_f;  

  x = g->Integral(20.,30.)/10.;
  dxl = x*(g->Integral(20.,24.)/g->Integral(20.,30.)*y_l[3]/y[3]+g->Integral(24.,30.)/g->Integral(20.,30.)*y_l[4]/y[4]);
  dxl_d = x*(g->Integral(20.,24.)/g->Integral(20.,30.)*y_l_d[3]/y[3]+g->Integral(24.,30.)/g->Integral(20.,30.)*y_l_d[4]/y[4]);
  dxl_f = x*(g->Integral(20.,24.)/g->Integral(20.,30.)*y_l_f[3]/y[3]+g->Integral(24.,30.)/g->Integral(20.,30.)*y_l_f[4]/y[4]);
  dxh = x*(g->Integral(20.,24.)/g->Integral(20.,30.)*y_h[3]/y[3]+g->Integral(24.,30.)/g->Integral(20.,30.)*y_h[4]/y[4]);
  dxh_d = x*(g->Integral(20.,24.)/g->Integral(20.,30.)*y_h_d[3]/y[3]+g->Integral(24.,30.)/g->Integral(20.,30.)*y_h_d[4]/y[4]);
  dxh_f = x*(g->Integral(20.,24.)/g->Integral(20.,30.)*y_h_f[3]/y[3]+g->Integral(24.,30.)/g->Integral(20.,30.)*y_h_f[4]/y[4]);
  y_new[4] = x, y_new_l[4] = dxl, y_new_h[4] = dxh, y_new_l_d[4] = dxl_d, y_new_h_d[4] = dxh_d, y_new_l_f[4] = dxl_f, y_new_h_f[4] = dxh_f;  

  x = g->Integral(30.,50.)/20.;
  dxl = x*(g->Integral(30.,40.)/g->Integral(30.,50.)*y_l[5]/y[5]+g->Integral(40.,50.)/g->Integral(30.,50.)*y_l[6]/y[6]);
  dxl_d = x*(g->Integral(30.,40.)/g->Integral(30.,50.)*y_l_d[5]/y[5]+g->Integral(40.,50.)/g->Integral(30.,50.)*y_l_d[6]/y[6]);
  dxl_f = x*(g->Integral(30.,40.)/g->Integral(30.,50.)*y_l_f[5]/y[5]+g->Integral(40.,50.)/g->Integral(30.,50.)*y_l_f[6]/y[6]);
  dxh = x*(g->Integral(30.,40.)/g->Integral(30.,50.)*y_h[5]/y[5]+g->Integral(40.,50.)/g->Integral(30.,50.)*y_h[6]/y[6]);
  dxh_d = x*(g->Integral(30.,40.)/g->Integral(30.,50.)*y_h_d[5]/y[5]+g->Integral(40.,50.)/g->Integral(30.,50.)*y_h_d[6]/y[6]);
  dxh_f = x*(g->Integral(30.,40.)/g->Integral(30.,50.)*y_h_f[5]/y[5]+g->Integral(40.,50.)/g->Integral(30.,50.)*y_h_f[6]/y[6]);
  y_new[5] = x, y_new_l[5] = dxl, y_new_h[5] = dxh, y_new_l_d[5] = dxl_d, y_new_h_d[5] = dxh_d, y_new_l_f[5] = dxl_f, y_new_h_f[5] = dxh_f;  

  x = g->Integral(50.,100.)/50.;
  dxl = x*(g->Integral(50.,60.)/g->Integral(50.,100.)*y_l[7]/y[7]+g->Integral(60.,70.)/g->Integral(50.,100.)*y_l[8]/y[8]+g->Integral(70.,100.)/g->Integral(50.,100.)*y_l[9]/y[9]);
  dxl_d = x*(g->Integral(50.,60.)/g->Integral(50.,100.)*y_l_d[7]/y[7]+g->Integral(60.,70.)/g->Integral(50.,100.)*y_l_d[8]/y[8]+g->Integral(70.,100.)/g->Integral(50.,100.)*y_l_d[9]/y[9]);
  dxl_f = x*(g->Integral(50.,60.)/g->Integral(50.,100.)*y_l_f[7]/y[7]+g->Integral(60.,70.)/g->Integral(50.,100.)*y_l_f[8]/y[8]+g->Integral(70.,100.)/g->Integral(50.,100.)*y_l_f[9]/y[9]);
  dxh = x*(g->Integral(50.,60.)/g->Integral(50.,100.)*y_h[7]/y[7]+g->Integral(60.,70.)/g->Integral(50.,100.)*y_h[8]/y[8]+g->Integral(70.,100.)/g->Integral(50.,100.)*y_h[9]/y[9]);
  dxh_d = x*(g->Integral(50.,60.)/g->Integral(50.,100.)*y_h_d[7]/y[7]+g->Integral(60.,70.)/g->Integral(50.,100.)*y_h_d[8]/y[8]+g->Integral(70.,100.)/g->Integral(50.,100.)*y_h_d[9]/y[9]);
  dxh_f = x*(g->Integral(50.,60.)/g->Integral(50.,100.)*y_h_f[7]/y[7]+g->Integral(60.,70.)/g->Integral(50.,100.)*y_h_f[8]/y[8]+g->Integral(70.,100.)/g->Integral(50.,100.)*y_h_f[9]/y[9]);
  y_new[6] = x, y_new_l[6] = dxl, y_new_h[6] = dxh, y_new_l_d[6] = dxl_d, y_new_h_d[6] = dxh_d, y_new_l_f[6] = dxl_f, y_new_h_f[6] = dxh_f;  

  //printf("\n");
  //printf("x-section and error after rebinning\n");

  for(i=0;i<nBins_new;i++)
    {
      //printf("pt %.0f-%.0f: %e (central) %e (low) %e (high)\n",ptBins_new[i],ptBins_new[i+1],y_new[i],y_new_l[i],y_new_h[i]);
      //printf("pt %.0f-%.0f: %e (low_data) %e (low_FONLL) %e (high_data) %e (high_FONLL)\n",ptBins_new[i],ptBins_new[i+1],y_new_l_d[i],y_new_l_f[i],y_new_h_d[i],y_new_h_f[i]);
    }

  TGraphAsymmErrors* real_ref_13 = new TGraphAsymmErrors(9, ptcent_13, y_13, ptwidth_13, ptwidth_13, y_13_l, y_13_h);
  real_ref_13->SetName("real_ref_13");
  real_ref_13->SetTitle("");
  real_ref_13->GetXaxis()->SetTitle("p_{t} (GeV)");
  real_ref_13->GetYaxis()->SetTitle("d#sigma/dp_{t} (pb/GeV)");
  real_ref_13->GetYaxis()->SetTitleOffset(1.1);
  real_ref_13->GetXaxis()->SetRange(0,110);
  real_ref_13->SetLineColor(kRed);

  TGraphAsymmErrors* real_ref_7 = new TGraphAsymmErrors(5, ptcent_7, y_7, ptwidth_7, ptwidth_7, y_7_l, y_7_h);
  real_ref_7->SetName("real_ref_7");
  real_ref_7->SetTitle("");
  real_ref_7->GetXaxis()->SetTitle("p_{t} (GeV)");
  real_ref_7->GetYaxis()->SetTitle("d#sigma/dp_{t} (pb/GeV)");
  real_ref_7->GetYaxis()->SetTitleOffset(1.1);
  real_ref_7->GetXaxis()->SetRange(0,110);
  real_ref_7->SetLineColor(kGreen);

  TGraphAsymmErrors* real_ref = new TGraphAsymmErrors(nBins, ptcent, y, ptwidth, ptwidth, y_l, y_h);
  real_ref->SetName("real_ref");
  real_ref->SetTitle("");
  real_ref->GetXaxis()->SetTitle("p_{t} (GeV)");
  real_ref->GetYaxis()->SetTitle("d#sigma/dp_{t} (pb/GeV)");
  real_ref->GetYaxis()->SetTitleOffset(1.1);
  real_ref->GetXaxis()->SetRange(0,110);
  real_ref->SetLineColor(kBlack);

  TGraphAsymmErrors* real_ref_new = new TGraphAsymmErrors(nBins_new, ptcent_new, y_new, ptwidth_new, ptwidth_new, y_new_l, y_new_h);
  real_ref_new->SetName("real_ref_new");
  real_ref_new->SetTitle("");
  real_ref_new->GetXaxis()->SetTitle("p_{t} (GeV)");
  real_ref_new->GetYaxis()->SetTitle("d#sigma/dp_{t} (pb/GeV)");
  real_ref_new->GetYaxis()->SetTitleOffset(1.1);
  real_ref_new->GetXaxis()->SetRange(0,110);
  real_ref_new->SetLineColor(kRed);
  real_ref_new->SetFillColor(kOrange);

  TGraphAsymmErrors* real_ref_new_data = new TGraphAsymmErrors(nBins_new, ptcent_new, y_new, ptwidth_new, ptwidth_new, y_new_l_d, y_new_h_d);
  real_ref_new_data->SetName("real_ref_new_data");
  real_ref_new_data->SetTitle("");
  real_ref_new_data->GetXaxis()->SetTitle("p_{t} (GeV)");
  real_ref_new_data->GetYaxis()->SetTitle("d#sigma/dp_{t} (pb/GeV)");
  real_ref_new_data->GetYaxis()->SetTitleOffset(1.1);
  real_ref_new_data->GetXaxis()->SetRange(0,110);
  real_ref_new_data->SetLineColor(kRed);
  real_ref_new_data->SetFillColor(kOrange);

  TGraphAsymmErrors* real_ref_new_FONLL = new TGraphAsymmErrors(nBins_new, ptcent_new, y_new, ptwidth_new, ptwidth_new, y_new_l_f, y_new_h_f);
  real_ref_new_FONLL->SetName("real_ref_new_FONLL");
  real_ref_new_FONLL->SetTitle("");
  real_ref_new_FONLL->GetXaxis()->SetTitle("p_{t} (GeV)");
  real_ref_new_FONLL->GetYaxis()->SetTitle("d#sigma/dp_{t} (pb/GeV)");
  real_ref_new_FONLL->GetYaxis()->SetTitleOffset(1.1);
  real_ref_new_FONLL->GetXaxis()->SetRange(0,110);
  real_ref_new_FONLL->SetLineColor(kRed);
  real_ref_new_FONLL->SetFillColor(kOrange);

  // Bpt reweight start

  TH1D* real_ref_new_sym = new TH1D("","",nBins_new,ptBins_new);
  real_ref_new_sym->SetName("real_ref_new_sym");
  real_ref_new_sym->SetTitle("");
  real_ref_new_sym->GetXaxis()->SetTitle("p_{t} (GeV)");
  real_ref_new_sym->GetYaxis()->SetTitle("d#sigma/dp_{t} (pb/GeV)");
  real_ref_new_sym->GetYaxis()->SetTitleOffset(1.1);
  real_ref_new_sym->SetLineColor(kBlue);

  for(i=0;i<nBins_new;i++)
    {
      real_ref_new_sym->SetBinContent(i+1,y_new[i]);
      real_ref_new_sym->SetBinError(i+1,0.5*(y_new_l[i]+y_new_h[i]));
    }

  //TF1* g2 = new TF1("g2","[0]+TMath::Exp([1]+[2]*x+[3]/x+[4]/(x*x))",5,100);
  TF1* g2 = new TF1("g2","TMath::Exp([0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x)",5,100);
  g2->SetParLimits(0,-100,100);
  g2->SetParLimits(1,-100,100);
  g2->SetParLimits(2,-100,100);
  g2->SetParLimits(3,-100,100);
  g2->SetParLimits(4,-100,100);
  real_ref_new_sym->Fit(g2,"IR");

  TCanvas* cg2 = new TCanvas("","",600,600);
  cg2->cd();
  cg2->SetLogy();
  real_ref_new_sym->Draw();
  g2->Draw("same");
  cg2->SaveAs("pp_continuous2.png");
  cg2->SaveAs("pp_continuous2.pdf");

  printf("\n");
  printf("x-section and error (Symmetric) after rebinning\n");

  for(i=0;i<nBins_new;i++)
    {
      printf("pt %.0f-%.0f: %e (central) %e (error)\n",ptBins_new[i],ptBins_new[i+1],real_ref_new_sym->GetBinContent(i+1),real_ref_new_sym->GetBinError(i+1));
    }

  TH1D* real_ref_new_sym_record = (TH1D*) real_ref_new_sym->Clone("real_ref_new_sym_record");

  real_ref_new_sym->Sumw2();
  real_ref_new_sym->Scale(1.0/real_ref_new_sym->Integral("width"));

  TFile* inputFONLL = new TFile("fonll_Bplus_dsigmadpt.root");
  TH1D* FONLL_sym = (TH1D*) inputFONLL->Get("gaeSigmaBplus_sym");
  FONLL_sym->Sumw2();
  FONLL_sym->Scale(1.0/FONLL_sym->Integral("width"));

  TFile* input2015PbPb = new TFile("ptshape/CrossSectionPbPb_2015PbPb.root");
  TGraphAsymmErrors* Bplus2015PbPb = (TGraphAsymmErrors*) input2015PbPb->Get("gaeSigmaBplus");
  TH1D* Bplus2015PbPb_sym = new TH1D("","",nBins_2015,ptBins_2015);
  Bplus2015PbPb_sym->SetName("Bplus2015PbPb_sym");
  Bplus2015PbPb_sym->SetTitle("");
  Bplus2015PbPb_sym->GetXaxis()->SetTitle("p_{t} (GeV)");
  Bplus2015PbPb_sym->GetYaxis()->SetTitle("d#sigma/dp_{t} (pb/GeV)");
  Bplus2015PbPb_sym->GetYaxis()->SetTitleOffset(1.1);

  for (i=0;i<Bplus2015PbPb->GetN();i++){
    //cout<<Bplus2015PbPb->GetY()[i]<<endl;
    Bplus2015PbPb_sym->SetBinContent(i+1,Bplus2015PbPb->GetY()[i]);
    Bplus2015PbPb_sym->SetBinError(i+1,0.5*((Bplus2015PbPb->GetErrorYhigh(i))+(Bplus2015PbPb->GetErrorYlow(i))));
    //cout<<Bplus2015PbPb_sym->GetBinContent(i+1)<<endl;
  }  
  Bplus2015PbPb_sym->Sumw2();
  Bplus2015PbPb_sym->Scale(1.0/Bplus2015PbPb_sym->Integral("width"));

  TCanvas* c1 = new TCanvas("c1","",600,600);
  c1->cd();
  real_ref_new_sym->Draw("ep");
  c1->SetLogy();                                                                                                                             

  c1->SaveAs("Bptweight1.png");
  c1->SaveAs("Bptweight1.pdf");  

  TFile* inputMC = new TFile(inputmcname.Data());
  TTree* Gen = (TTree*) inputMC->Get("Bfinder/ntGen");
  Gen->AddFriend("hiEvtAnalyzer/HiTree");
  Gen->AddFriend("hltanalysis/HltTree");
  Gen->AddFriend("skimanalysis/HltTree");
  TString weightgen = "pthatweight";
  //TString Bselection = "(GisSignal==1)&&(GcollisionId==0)&&((Gpt>5&&Gpt<10&&TMath::Abs(Gy)<2.4)||(Gpt>10&&Gpt<17&&TMath::Abs(Gy)<1.45)||(Gpt>17&&Gpt<100&&TMath::Abs(Gy)<2.1))";
  TString Bselection = "TMath::Abs(Gy)<2.4 && TMath::Abs(GpdgId)==521 && GisSignal==1 && GcollisionId==0";
  TH1D* Genpt = new TH1D("Genpt","",nBins_new,ptBins_new);
  //TH1D* Genpt = new TH1D("Genpt","",nBins_2015,ptBins_2015);
  Gen->Project("Genpt","Gpt",TCut(weightgen)*TCut(Bselection));
  Genpt->GetXaxis()->SetTitle("p_{t} (GeV)");
  Genpt->GetYaxis()->SetTitle("d#sigma/dp_{t} (pb/GeV)");
  //Genpt->Scale(1.0/Genpt->Integral(),"width");
  divideBinWidth(Genpt);
  TFile* inputMC = new TFile(inputmcname.Data());
  TTree* Gen = (TTree*) inputMC->Get("Bfinder/ntGen");
  Gen->AddFriend("hiEvtAnalyzer/HiTree");
  Gen->AddFriend("hltanalysis/HltTree");
  Gen->AddFriend("skimanalysis/HltTree");
  TString weightgen = "pthatweight";
  //TString Bselection = "(GisSignal==1)&&(GcollisionId==0)&&((Gpt>5&&Gpt<10&&TMath::Abs(Gy)<2.4)||(Gpt>10&&Gpt<17&&TMath::Abs(Gy)<1.45)||(Gpt>17&&Gpt<100&&TMath::Abs(Gy)<2.1))";                                                                                                             
  TString Bselection = "TMath::Abs(Gy)<2.4 && TMath::Abs(GpdgId)==521 && GisSignal==1 && GcollisionId==0";
  TH1D* Genpt = new TH1D("Genpt","",nBins_new,ptBins_new);
  //TH1D* Genpt = new TH1D("Genpt","",nBins_2015,ptBins_2015);                                                                                  
  Gen->Project("Genpt","Gpt",TCut(weightgen)*TCut(Bselection));
  Genpt->GetXaxis()->SetTitle("p_{t} (GeV)");
  Genpt->GetYaxis()->SetTitle("d#sigma/dp_{t} (pb/GeV)");
  //Genpt->Scale(1.0/Genpt->Integral(),"width");                                                                                                
  divideBinWidth(Genpt);
  Genpt->Scale(1.0/Genpt->Integral("width"));
  Genpt->Scale(1.0/Genpt->Integral("width"));

  /*
  TFile* MCoff = new TFile("MCstudiesPbPb.root");
  TH1D* hEffoff = (TH1D*) MCoff->Get("hEff");
  TFile* MCpri = new TFile("MCstudiesPbPb_private.root");
  TH1D* hEffpri = (TH1D*) MCpri->Get("hEff");
  hEffoff->Sumw2();
  hEffpri->Sumw2();
  hEffoff->Divide(hEffpri);

  TCanvas* c100 = new TCanvas("c100","",600,600);
  c100->cd();
  hEffoff->Draw("ep");
  c100->SaveAs("effcheck.png");
  */

  TCanvas* c2 = new TCanvas("c2","",600,600);
  c2->cd();
  Genpt->Draw("ep");
  //real_ref_new_sym->Draw("ep same");
  c2->SetLogy();                                                                                                                             

  c2->SaveAs("Bptweight2.png");
  c2->SaveAs("Bptweight2.pdf");

  //TH1D* Ratiopt = (TH1D*) FONLL_sym->Clone("Ratiopt");
  TH1D* Ratiopt = (TH1D*) real_ref_new_sym->Clone("Ratiopt");
  //TH1D* Ratiopt = (TH1D*) Bplus2015PbPb_sym->Clone("Ratiopt");
  Ratiopt->Sumw2();
  Ratiopt->Divide(Genpt);
  //Ratiopt->Scale(1.0/Ratiopt->Integral("width"));

  TF1 *f = new TF1("f1","([0]+[1]*x+[2]*x*x)*TMath::Exp(-[3]*x)+[4]",5,100);
  //f->SetParameters(3,1,0.1,0.1,1);
  f->SetParLimits(0,-10,10);
  f->SetParLimits(1,-1,1);
  f->SetParLimits(2,-1,1);
  f->SetParLimits(3,-1,1);
  f->SetParLimits(4,0,2);

  Ratiopt->Fit(f,"IR");

  printf("Bpt weight function(MC Gpt): (%f+%f*Bpt+%f*Bpt*Bpt)*TMath::Exp(-%f*Bpt)+%f\n",f->GetParameter(0),f->GetParameter(1),f->GetParameter(2),f->GetParameter(3),f->GetParameter(4));

  Ratiopt->SetTitle("");
  Ratiopt->GetXaxis()->SetTitle("p_{t} (GeV/c)");
  Ratiopt->GetYaxis()->SetTitle("Ref/Gen");
  Ratiopt->GetYaxis()->SetTitleOffset(1.0);  

  TCanvas* c3 = new TCanvas("c3","",600,600);
  c3->cd();
  f->SetLineColor(kRed);
  //Ratiopt->SetAxisRange(0.1,10,"Y"); 
  Ratiopt->Draw();  
  f->Draw("same");
  //c3->SetLogy();
  c3->SaveAs("Bptweight3.png");
  c3->SaveAs("Bptweight3.pdf");

  // Bpt reweight finished

  TFile* pp2015file = new TFile("ROOTfiles/CrossSectionPP.root");
  TH1D* pp2015_data = (TH1D*) pp2015file->Get("hPtSigma");
  pp2015_data->SetName("pp2015_data");
  pp2015_data->SetTitle("");
  pp2015_data->GetXaxis()->SetTitle("p_{t} (GeV)");
  pp2015_data->GetYaxis()->SetTitle("d#sigma/dp_{t} (pb/GeV)");
  pp2015_data->GetYaxis()->SetTitleOffset(1.1);
  pp2015_data->GetXaxis()->SetRange(0,110);
  pp2015_data->SetLineColor(kBlue);
  //pp2015_data->SetFillColor(kOrange);
  pp2015_data->SetLineWidth(1.0);
  pp2015_data->SetMarkerStyle(20);
  pp2015_data->SetMarkerSize(0.1);

  TGraphAsymmErrors* pp2015_syst = (TGraphAsymmErrors*) pp2015file->Get("gaeCrossSyst");
  pp2015_syst->SetName("pp2015_syst");
  pp2015_syst->SetTitle("");
  pp2015_syst->GetXaxis()->SetTitle("p_{t} (GeV)");
  pp2015_syst->GetYaxis()->SetTitle("d#sigma/dp_{t} (pb/GeV)");
  pp2015_syst->GetYaxis()->SetTitleOffset(1.1);
  pp2015_syst->GetXaxis()->SetRange(0,110);
  pp2015_syst->SetLineColor(kBlue);
  pp2015_syst->SetFillColor(kBlack);
  pp2015_syst->SetLineWidth(1.0);
  pp2015_syst->SetMarkerStyle(20);
  pp2015_syst->SetMarkerSize(0.1);

  TCanvas* c4 = new TCanvas("c4","",600,600);
  c4->cd();
  //real_ref_13->Draw("||");
  //real_ref_7->Draw("|| same");
  //real_ref->Draw("|| same");
  //pp2015->Draw("A2");
  //pp2015->Draw("p same");
  real_ref_new_data->Draw("A2");
  real_ref_new_FONLL->Draw("ep same");
  //pp2015_data->Draw("ep same");
  //pp2015_syst->Draw("ep same");
  FONLL_sym->Draw("ep same");

  TLegend *leg = new TLegend(0.45,0.70,0.75,0.80,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  //leg->AddEntry(real_ref_13,"scaled reference from 13TeV","l");
  //leg->AddEntry(real_ref_7,"scaled reference from 7TeV","l");
  //leg->AddEntry(real_ref,"scaled reference combined","l");
  leg->AddEntry(real_ref_new_data,"scaled_data error","f");
  leg->AddEntry(real_ref_new_FONLL,"scaled_FONLL error","l");
  //leg->AddEntry(pp2015_syst,"2015 pp cross-section","l");
  leg->AddEntry(FONLL_sym,"FONLL reference","l");
  leg->Draw("same");

  TLatex* texcms = new TLatex(0.18,0.95,"CMS");
  texcms->SetNDC();
  texcms->SetTextAlign(13);
  texcms->SetTextFont(62);
  texcms->SetTextSize(0.06);
  texcms->SetLineWidth(2);
  texcms->Draw();

  TLatex* texB = new TLatex(0.65,0.84,"B^{+}+B^{-} |y|<2.4");
  texB->SetNDC();  
  texB->SetTextFont(42);
  texB->SetTextSize(0.05);
  texB->SetLineWidth(2);
  texB->Draw();

  TLatex* texCol;
  texCol= new TLatex(0.96,0.94, Form("28.0 pb^{-1} (%s 5.02 TeV)","pp"));
  texCol->SetNDC();
  texCol->SetTextAlign(32);
  texCol->SetTextSize(0.05);
  texCol->SetTextFont(42);
  texCol->Draw();

  c4->RedrawAxis();
  c4->SetLogy();
  c4->SaveAs("pp_reference.png");
  c4->SaveAs("pp_reference.pdf");
  
  TFile* output_ref = new TFile(outputname_ref.Data(),"recreate");
  output_ref->cd();
  real_ref_13->Write();
  real_ref_7->Write();
  real_ref->Write();
  real_ref_new->Write();
  real_ref_new_data->Write();
  real_ref_new_FONLL->Write();
  pp2015_data->Write();
  pp2015_syst->Write();
  real_ref_new_sym_record->Write();
  real_ref_new_sym->Write();
  Genpt->Write();
  Ratiopt->Write();
  FONLL_sym->Write();
  output_ref->Close();


  
  return;
}

int main()
{
  ppref_real();
  return 0;
}
