#include "uti.h"
#include "TF1.h"
#include <TFitResultPtr.h>
#include "ntuple2.h"

TString inputdataname = "/raid5/data/gwangjun/selected_data_ntKp_PbPb_2018_MuonJSON_BDT.root";//Muon JSON
TString inputeffname1, inputeffname2, inputTNPname1, inputTNPname2, outputeffname;
Float_t hiBinMin,hiBinMax,centMin,centMax;

bool useTNP = 1;
bool useFiducial = 1;
bool useToy = 1;

void AverageEff(int usePbPb=0, TString inputdata="" , TString inputmc="", TString trgselection="1",  TString cut="", TString cutmcgen="", int isMC=0, Double_t luminosity=1., int doweight=0, TString collsyst="PbPb", TString outputfile="", TString npfit="0", int doDataCor = 0, Float_t centmin=0., Float_t centmax=90.)
{
  outputeffname = Form("ROOTfiles/MCstudiesPbPbAverage_toy_pt1050_Cent%.0f-%.0f.root",centmin,centmax);

  inputeffname1 = "ROOTfiles/GenStatSyst_Cent0-30.root";
  inputeffname2 = "ROOTfiles/GenStatSyst_Cent30-90.root";

  inputTNPname1 = "ROOTfiles/TNP2D_Bplus_Cent0-30.root";
  inputTNPname2 = "ROOTfiles/TNP2D_Bplus_Cent30-90.root";

  hiBinMin = centmin*2;
  hiBinMax = centmax*2;
  centMin = centmin;
  centMax = centmax;

  if (!(usePbPb==1||usePbPb==0)) std::cout<<"ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!, you are using a non valid isPbPb option"<<std::endl;
  bool isPbPb=(bool)(usePbPb);

  gStyle->SetTextSize(0.05);
  gStyle->SetTextFont(42);
  gStyle->SetPadRightMargin(0.043);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.145);
  gStyle->SetTitleX(.0f);

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  
  TFile* inf = new TFile(inputdataname.Data());
  TTree* nt = (TTree*)inf->Get("ntKp");
  
  setbranchaddress(inf,nt);

  TFile* Eff2Dfile1 = new TFile(inputeffname1.Data());
  TFile* Eff2Dfile2 = new TFile(inputeffname2.Data());
  TH2D* invEff2D1;
  TH2D* invEff2D2;

  int ntoy = 10000;

  TFile* TNP2D1 = new TFile(inputTNPname1.Data());
  TH2D* tnp_scale1 = (TH2D*)TNP2D1->Get("tnp_scale");
  TH2D* tnp_total_d1 = (TH2D*)TNP2D1->Get("tnp_total_d");
  TH2D* tnp_total_u1 = (TH2D*)TNP2D1->Get("tnp_total_u");
  TFile* TNP2D2 = new TFile(inputTNPname2.Data());
  TH2D* tnp_scale2 = (TH2D*)TNP2D2->Get("tnp_scale");
  TH2D* tnp_total_d2 = (TH2D*)TNP2D2->Get("tnp_total_d");
  TH2D* tnp_total_u2 = (TH2D*)TNP2D2->Get("tnp_total_u");

  
  int nevts[ntoy][nBins];
  float sum[ntoy][nBins], sumerr[ntoy][nBins];
  for(int i=0;i<ntoy;i++)
    {
      for(int j=0;j<nBins;j++)
	{
	  nevts[i][j]=0;
	  sum[i][j]=0.0;
	  sumerr[i][j]=0.0;
	}
    }

  for(int i=0;i<nt->GetEntries();i++)
    {
      nt->GetEntry(i);
      if(!(TMath::Abs(Bmass-5.27932)<0.08)) continue;
      if(!(HiBin>=centMin && HiBin<centMax)) continue;
      if(useFiducial) {if(!((Bpt>5 && Bpt<10 && TMath::Abs(By)>1.5) || (Bpt>10))) continue;}// check if already implemented

      for(int j=0;j<ntoy;j++)
	{
	  invEff2D1 = (TH2D*)Eff2Dfile1->Get(Form("EffBptByInvTrial%d",j));
	  invEff2D2 = (TH2D*)Eff2Dfile2->Get(Form("EffBptByInvTrial%d",j));
	  //invEff2D1 = (TH2D*)Eff2Dfile1->Get(Form("AccBptByInvTrial%d",j));
	  //invEff2D2 = (TH2D*)Eff2Dfile2->Get(Form("AccBptByInvTrial%d",j));
	  //invEff2D1 = (TH2D*)Eff2Dfile1->Get(Form("EffonlyBptByInvTrial%d",j));
	  //invEff2D2 = (TH2D*)Eff2Dfile2->Get(Form("EffonlyBptByInvTrial%d",j));

	  for(int k=0;k<nBins;k++)
	    {
	      if(!(Bpt>ptBins[k] && Bpt<ptBins[k+1])) continue;
	      
	      if(HiBin<60)
		{
		  if(invEff2D1->GetBinContent(invEff2D1->GetBin(findBptbin(Bpt),findBybin(By)))!=0)
		    {
		      nevts[j][k]++;
		      
		      if(!useTNP)
			{
			  sum[j][k]+=invEff2D1->GetBinContent(invEff2D1->GetBin(findBptbin(Bpt),findBybin(By)));
			  sumerr[j][k]+=invEff2D1->GetBinError(invEff2D1->GetBin(findBptbin(Bpt),findBybin(By)))*invEff2D1->GetBinError(invEff2D1->GetBin(findBptbin(Bpt),findBybin(By)));
			}
		      
		      if(useTNP)
			{
			  sum[j][k]+=invEff2D1->GetBinContent(invEff2D1->GetBin(findBptbin(Bpt),findBybin(By)))/tnp_scale1->GetBinContent(tnp_scale1->GetBin(findBptbin(Bpt),findBybin(By)));
			  sumerr[j][k]+=(invEff2D1->GetBinError(invEff2D1->GetBin(findBptbin(Bpt),findBybin(By)))/tnp_scale1->GetBinContent(tnp_scale1->GetBin(findBptbin(Bpt),findBybin(By))))*(invEff2D1->GetBinError(invEff2D1->GetBin(findBptbin(Bpt),findBybin(By)))/tnp_scale1->GetBinContent(tnp_scale1->GetBin(findBptbin(Bpt),findBybin(By))));
			}
		    }
		}
	      
	      
	      if(HiBin>=60)
		{
		  if(invEff2D2->GetBinContent(invEff2D2->GetBin(findBptbin(Bpt),findBybin(By)))!=0)
		    {
		      nevts[j][k]++;
		      
		      if(!useTNP)
			{
			  sum[j][k]+=invEff2D2->GetBinContent(invEff2D2->GetBin(findBptbin(Bpt),findBybin(By)));
			  sumerr[j][k]+=invEff2D2->GetBinError(invEff2D2->GetBin(findBptbin(Bpt),findBybin(By)))*invEff2D2->GetBinError(invEff2D2->GetBin(findBptbin(Bpt),findBybin(By)));
		    }
		      
		      if(useTNP)
			{
			  sum[j][k]+=invEff2D2->GetBinContent(invEff2D2->GetBin(findBptbin(Bpt),findBybin(By)))/tnp_scale2->GetBinContent(tnp_scale2->GetBin(findBptbin(Bpt),findBybin(By)));
			  sumerr[j][k]+=(invEff2D2->GetBinError(invEff2D2->GetBin(findBptbin(Bpt),findBybin(By)))/tnp_scale2->GetBinContent(tnp_scale2->GetBin(findBptbin(Bpt),findBybin(By))))*(invEff2D2->GetBinError(invEff2D2->GetBin(findBptbin(Bpt),findBybin(By)))/tnp_scale2->GetBinContent(tnp_scale2->GetBin(findBptbin(Bpt),findBybin(By))));
			}
		    }
		}
	    }	      
	  
	}
    }

  TFile* outf = new TFile(outputeffname.Data(),"recreate");
  outf->cd();
  
  TH1F* invEffavetoy[nBins];
  TCanvas* cinvEfftoy[nBins];
  TLine* invEffnominal[nBins];
  TString texper = "%";

  //Inclusive pt
  //double Nominal[nBins] = {26.54};//pt7-50 Cent0-30%
  //double Nominal[nBins] = {17.76};//pt9-50 Cent30-90%
  //double Nominal[nBins] = {25.38};//pt7-50 Cent0-90%
  double Nominal[nBins] = {17.55};//pt10-50 Cent0-30%
  //double Nominal[nBins] = {16.75};//pt10-50 Cent30-90%
  //double Nominal[nBins] = {17.07};//pt10-50 Cent0-90%

  double Min[nBins] = {0};
  double Max[nBins] = {30};

  //Bsbin
  //double Nominal[nBins] = {115.34, 32.79, 10.15, 5.64};//0-90%
  //double Min[nBins] = {0, 20, 5, 0};
  //double Max[nBins] = {200, 50, 15, 10};

  //Bplusbin
  //double Nominal[nBins] = {494.35, 115.34, 32.79, 10.15, 6.22, 4.33, 4.09, 4.38};//0-90%
  //double Min[nBins] = {0, 0, 20, 5, 0, 0, 0, 0};
  //double Max[nBins] = {1000, 200, 50, 15, 10, 10, 10, 10};

  //double Nominal[nBins] = {563.58, 100.27, 32.79, 10.15, 6.22, 4.33, 4.09, 4.38};//0-90%
  //double Min[nBins] = {0, 0, 20, 5, 0, 0, 0, 0};
  //double Max[nBins] = {1000, 200, 50, 15, 10, 10, 10, 10};


  for(int i=0;i<nBins;i++)
    {
      invEffavetoy[i] = new TH1F(Form("invEffavetoy%d",i+1),"",200,Min[i],Max[i]);
      invEffavetoy[i]->GetXaxis()->SetTitle("<1/(acc x eff)>");
      invEffavetoy[i]->GetYaxis()->SetTitle("Counts");
      invEffavetoy[i]->SetTitle(Form("MC Smeared Distribution for %.0f < Bpt < %.0f",ptBins[i],ptBins[i+1]));

      invEffavetoy[i]->GetXaxis()->CenterTitle();
      invEffavetoy[i]->GetYaxis()->CenterTitle();
      invEffavetoy[i]->SetTitleOffset(1.4);

      //std::cout<<"sum: "<<sum[i]<<" nevts: "<<nevts[i]<<std::endl;
      for(int j=0;j<ntoy;j++)
	{
	  sum[j][i]=sum[j][i]/nevts[j][i];
	  sumerr[j][i]=sqrt(sumerr[j][i])/nevts[j][i];

	  invEffavetoy[i]->Fill(sum[j][i]);
	}
      invEffavetoy[i]->Write();
      
      invEffnominal[i] = new TLine(Nominal[i],0,Nominal[i],invEffavetoy[i]->GetMaximum());
      invEffnominal[i]->SetLineStyle(2);
      invEffnominal[i]->SetLineWidth(2);
      invEffnominal[i]->SetLineColor(2);

      cinvEfftoy[i] = new TCanvas("","",600,600);
      cinvEfftoy[i]->cd();
      invEffavetoy[i]->Draw();
      invEffnominal[i]->Draw();
      cinvEfftoy[i]->SaveAs(Form("plotAverageEff/toyMC/toy_invEff_pt%.0f-%.0f_Cent%.0f-%.0f.png",ptBins[i],ptBins[i+1],centmin,centmax));
      cinvEfftoy[i]->SaveAs(Form("plotAverageEff/toyMC/toy_invEff_pt%.0f-%.0f_Cent%.0f-%.0f.pdf",ptBins[i],ptBins[i+1],centmin,centmax));
    }

  for(int i=0;i<nBins;i++)
    {
      printf("pt %.0f-%.0f Cent %.0f-%.0f%s toy syst = %f%s\n",ptBins[i],ptBins[i+1],centmin,centmax,texper.Data(),100*invEffavetoy[i]->GetRMS()/Nominal[i],texper.Data());
      //printf("pt %.0f-%.0f Cent %.0f-%.0f%s toy syst = %f%s\n",ptBins[i],ptBins[i+1],centmin,centmax,texper.Data(),100*invEffavetoy[i]->GetRMS()/invEffavetoy[i]->GetMean(),texper.Data());
    }

  outf->Close(); 
}


int main(int argc, char *argv[])
{
  if(argc==16)
    {
      AverageEff(atoi(argv[1]), argv[2], argv[3], argv[4], argv[5], argv[6], atoi(argv[7]), atof(argv[8]), atoi(argv[9]), argv[10], argv[11], argv[12], atoi(argv[13]), atof(argv[14]), atof(argv[15]));
      return 0;
    }
  else if(argc==14)
    {
      AverageEff(atoi(argv[1]), argv[2], argv[3], argv[4], argv[5], argv[6], atoi(argv[7]), atof(argv[8]), atoi(argv[9]), argv[10], argv[11], argv[12], atoi(argv[13]));
      return 0;
    }
  else
    {
      std::cout << "Wrong number of inputs" << std::endl;
      return 1;
    }
}
