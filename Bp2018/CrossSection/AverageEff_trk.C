#include "uti.h"
#include "TF1.h"
#include <TFitResultPtr.h>
#include "ntuple2.h"

TString inputdataname = "/raid5/data/gwangjun/selected_data_ntKp_PbPb_2018_MuonJSON_BDT.root";//Muon JSON
Float_t hiBinMin,hiBinMax,centMin,centMax;

bool uselogy = 0;

bool useTNP = 1;
bool useFiducial = 0;
bool compareBpt = 0;
bool compareSplot = 0;

TString trkvarname;

void AverageEff(int usePbPb=0, TString inputdata="" , TString inputmc="", TString trgselection="1",  TString cut="", TString cutmcgen="", int isMC=0, Double_t luminosity=1., int doweight=0, TString collsyst="PbPb", TString outputfile="", TString npfit="0", int doDataCor = 0, Float_t centmin=0., Float_t centmax=90.)
{
  TString plotname = Form("plotAverageEff/Ave/invEff_pt1050_Cent%.0f-%.0f",centmin,centmax);
  TString inputTNPname1 = "ROOTfiles/TNP2D_Bplus_Cent0-30.root";
  TString inputTNPname2 = "ROOTfiles/TNP2D_Bplus_Cent30-90.root";

  TString outputeffname;
  TString nominaleffname;
  TString inputeffname1;
  TString inputeffname2;

  int nevts[nBins];
  float sum[nBins], sumerr[nBins], sum_d[nBins], sumerr_d[nBins], sum_u[nBins], sumerr_u[nBins];

  hiBinMin = centmin*2;
  hiBinMax = centmax*2;
  centMin = centmin;
  centMax = centmax;

  gStyle->SetOptStat(0);

  gStyle->SetTextSize(0.05);
  gStyle->SetTextFont(42);
  gStyle->SetPadRightMargin(0.043);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.145);
  gStyle->SetTitleX(.0f);

  TFile* inf = new TFile(inputdataname.Data());
  TTree* nt = (TTree*)inf->Get("ntKp");
  setbranchaddress(inf,nt);

  //////////////////////////////////////////////////////////////////////////////////////////////
  
  float maxdev[nBins];
  for(int n=0;n<nBins;n++)
    {
      maxdev[n]=0;
    }

  for(int trk=1;trk<8;trk++)
    {
      if(trk==1) trkvarname = "Btrk1Eta";
      if(trk==2) trkvarname = "Btrk1Y";
      if(trk==3) trkvarname = "Btrk1Pt";
      if(trk==4) trkvarname = "Btrk1Dz1";
      if(trk==5) trkvarname = "Btrk1DzError1";
      if(trk==6) trkvarname = "Btrk1Dxy1";
      if(trk==7) trkvarname = "Btrk1DxyError1";

      outputeffname = Form("ROOTfiles/obsolete/MCstudiesPbPbAverage_%s_pt1050_Cent%.0f-%.0f.root",trkvarname.Data(),centmin,centmax);
      nominaleffname = Form("ROOTfiles/MCstudiesPbPbAverage_pt1050_Cent%.0f-%.0f.root",centmin,centmax);
      
      inputeffname1 = Form("ROOTfiles/obsolete/MCstudiesPbPb_%s_Fid2D_Cent0-30.root",trkvarname.Data());
      inputeffname2 = Form("ROOTfiles/obsolete/MCstudiesPbPb_%s_Fid2D_Cent30-90.root",trkvarname.Data());
 
      TFile* Eff2Dfile1 = new TFile(inputeffname1.Data());
      TH2D* hEff2D1 = (TH2D*)Eff2Dfile1->Get("hEff2D");
      TH2D* invEff2D1 = (TH2D*)Eff2Dfile1->Get("invEff2D");
      TFile* TNP2D1 = new TFile(inputTNPname1.Data());
      TH2D* tnp_scale1 = (TH2D*)TNP2D1->Get("tnp_scale");
      TH2D* tnp_total_d1 = (TH2D*)TNP2D1->Get("tnp_total_d");
      TH2D* tnp_total_u1 = (TH2D*)TNP2D1->Get("tnp_total_u");
      
      TFile* Eff2Dfile2 = new TFile(inputeffname2.Data());
      TH2D* hEff2D2 = (TH2D*)Eff2Dfile2->Get("hEff2D");
      TH2D* invEff2D2 = (TH2D*)Eff2Dfile2->Get("invEff2D");
      TFile* TNP2D2 = new TFile(inputTNPname2.Data());
      TH2D* tnp_scale2 = (TH2D*)TNP2D2->Get("tnp_scale");
      TH2D* tnp_total_d2 = (TH2D*)TNP2D2->Get("tnp_total_d");
      TH2D* tnp_total_u2 = (TH2D*)TNP2D2->Get("tnp_total_u");
  
      for(int i=0;i<nBins;i++)
	{
	  nevts[i]=0;
	  sum[i]=0.0;
	  sumerr[i]=0.0;
	  sum_d[i]=0.0;
	  sumerr_d[i]=0.0;
	  sum_u[i]=0.0;
	  sumerr_u[i]=0.0;
	}
      
      for(int i=0;i<nt->GetEntries();i++)
	{
	  nt->GetEntry(i);
	  if(!(TMath::Abs(Bmass-5.27932)<0.08)) continue;
	  if(!(HiBin>=centMin && HiBin<centMax)) continue;
	  if(useFiducial) {if(!((Bpt>5 && Bpt<10 && TMath::Abs(By)>1.5) || (Bpt>10))) continue;}// check if already implemented
	  
	  for(int k=0;k<nBins;k++)
	    {
	      if(!(Bpt>ptBins[k] && Bpt<ptBins[k+1])) continue;
	      
	      if(HiBin<60)
		{
		  if(hEff2D1->GetBinContent(invEff2D1->GetBin(findBptbin(Bpt),findBybin(By)))!=0)
		    {
		      nevts[k]++;
		  
		      if(!useTNP)
			{
			  sum[k]+=invEff2D1->GetBinContent(invEff2D1->GetBin(findBptbin(Bpt),findBybin(By)));
			  sumerr[k]+=invEff2D1->GetBinError(invEff2D1->GetBin(findBptbin(Bpt),findBybin(By)))*invEff2D1->GetBinError(invEff2D1->GetBin(findBptbin(Bpt),findBybin(By)));
			}
		      
		      if(useTNP)
			{
			  sum[k]+=invEff2D1->GetBinContent(invEff2D1->GetBin(findBptbin(Bpt),findBybin(By)))/tnp_scale1->GetBinContent(tnp_scale1->GetBin(findBptbin(Bpt),findBybin(By)));
			  sumerr[k]+=(invEff2D1->GetBinError(invEff2D1->GetBin(findBptbin(Bpt),findBybin(By)))/tnp_scale1->GetBinContent(tnp_scale1->GetBin(findBptbin(Bpt),findBybin(By))))*(invEff2D1->GetBinError(invEff2D1->GetBin(findBptbin(Bpt),findBybin(By)))/tnp_scale1->GetBinContent(tnp_scale1->GetBin(findBptbin(Bpt),findBybin(By))));
			  
			  sum_d[k]+=invEff2D1->GetBinContent(invEff2D1->GetBin(findBptbin(Bpt),findBybin(By)))/(tnp_scale1->GetBinContent(tnp_scale1->GetBin(findBptbin(Bpt),findBybin(By)))*(1.0-tnp_total_d1->GetBinContent(tnp_total_d1->GetBin(findBptbin(Bpt),findBybin(By)))));
			  sumerr_d[k]+=(invEff2D1->GetBinError(invEff2D1->GetBin(findBptbin(Bpt),findBybin(By)))/(tnp_scale1->GetBinContent(tnp_scale1->GetBin(findBptbin(Bpt),findBybin(By)))*(1.0-tnp_total_d1->GetBinContent(tnp_total_d1->GetBin(findBptbin(Bpt),findBybin(By))))))*(invEff2D1->GetBinError(invEff2D1->GetBin(findBptbin(Bpt),findBybin(By)))/(tnp_scale1->GetBinContent(tnp_scale1->GetBin(findBptbin(Bpt),findBybin(By)))*(1.0-tnp_total_d1->GetBinContent(tnp_total_d1->GetBin(findBptbin(Bpt),findBybin(By))))));
			  
			  sum_u[k]+=invEff2D1->GetBinContent(invEff2D1->GetBin(findBptbin(Bpt),findBybin(By)))/(tnp_scale1->GetBinContent(tnp_scale1->GetBin(findBptbin(Bpt),findBybin(By)))*(1.0+tnp_total_u1->GetBinContent(tnp_total_u1->GetBin(findBptbin(Bpt),findBybin(By)))));
			  sumerr_u[k]+=(invEff2D1->GetBinError(invEff2D1->GetBin(findBptbin(Bpt),findBybin(By)))/(tnp_scale1->GetBinContent(tnp_scale1->GetBin(findBptbin(Bpt),findBybin(By)))*(1.0+tnp_total_u1->GetBinContent(tnp_total_u1->GetBin(findBptbin(Bpt),findBybin(By))))))*(invEff2D1->GetBinError(invEff2D1->GetBin(findBptbin(Bpt),findBybin(By)))/(tnp_scale1->GetBinContent(tnp_scale1->GetBin(findBptbin(Bpt),findBybin(By)))*(1.0+tnp_total_u1->GetBinContent(tnp_total_u1->GetBin(findBptbin(Bpt),findBybin(By))))));
			}
		    }
		}
	      
	      
	      if(HiBin>=60)
		{
		  if(hEff2D2->GetBinContent(invEff2D2->GetBin(findBptbin(Bpt),findBybin(By)))!=0)
		    {
		      nevts[k]++;
		      
		      if(!useTNP)
			{
			  sum[k]+=invEff2D2->GetBinContent(invEff2D2->GetBin(findBptbin(Bpt),findBybin(By)));
			  sumerr[k]+=invEff2D2->GetBinError(invEff2D2->GetBin(findBptbin(Bpt),findBybin(By)))*invEff2D2->GetBinError(invEff2D2->GetBin(findBptbin(Bpt),findBybin(By)));
			}
		      
		      if(useTNP)
			{
			  sum[k]+=invEff2D2->GetBinContent(invEff2D2->GetBin(findBptbin(Bpt),findBybin(By)))/tnp_scale2->GetBinContent(tnp_scale2->GetBin(findBptbin(Bpt),findBybin(By)));
			  sumerr[k]+=(invEff2D2->GetBinError(invEff2D2->GetBin(findBptbin(Bpt),findBybin(By)))/tnp_scale2->GetBinContent(tnp_scale2->GetBin(findBptbin(Bpt),findBybin(By))))*(invEff2D2->GetBinError(invEff2D2->GetBin(findBptbin(Bpt),findBybin(By)))/tnp_scale2->GetBinContent(tnp_scale2->GetBin(findBptbin(Bpt),findBybin(By))));
			  
			  sum_d[k]+=invEff2D2->GetBinContent(invEff2D2->GetBin(findBptbin(Bpt),findBybin(By)))/(tnp_scale2->GetBinContent(tnp_scale2->GetBin(findBptbin(Bpt),findBybin(By)))*(1.0-tnp_total_d2->GetBinContent(tnp_total_d2->GetBin(findBptbin(Bpt),findBybin(By)))));
			  sumerr_d[k]+=(invEff2D2->GetBinError(invEff2D2->GetBin(findBptbin(Bpt),findBybin(By)))/(tnp_scale2->GetBinContent(tnp_scale2->GetBin(findBptbin(Bpt),findBybin(By)))*(1.0-tnp_total_d2->GetBinContent(tnp_total_d2->GetBin(findBptbin(Bpt),findBybin(By))))))*(invEff2D2->GetBinError(invEff2D2->GetBin(findBptbin(Bpt),findBybin(By)))/(tnp_scale2->GetBinContent(tnp_scale2->GetBin(findBptbin(Bpt),findBybin(By)))*(1.0-tnp_total_d2->GetBinContent(tnp_total_d2->GetBin(findBptbin(Bpt),findBybin(By))))));
			  
			  sum_u[k]+=invEff2D2->GetBinContent(invEff2D2->GetBin(findBptbin(Bpt),findBybin(By)))/(tnp_scale2->GetBinContent(tnp_scale2->GetBin(findBptbin(Bpt),findBybin(By)))*(1.0+tnp_total_u2->GetBinContent(tnp_total_u2->GetBin(findBptbin(Bpt),findBybin(By)))));
			  sumerr_u[k]+=(invEff2D2->GetBinError(invEff2D2->GetBin(findBptbin(Bpt),findBybin(By)))/(tnp_scale2->GetBinContent(tnp_scale2->GetBin(findBptbin(Bpt),findBybin(By)))*(1.0+tnp_total_u2->GetBinContent(tnp_total_u2->GetBin(findBptbin(Bpt),findBybin(By))))))*(invEff2D2->GetBinError(invEff2D2->GetBin(findBptbin(Bpt),findBybin(By)))/(tnp_scale2->GetBinContent(tnp_scale2->GetBin(findBptbin(Bpt),findBybin(By)))*(1.0+tnp_total_u2->GetBinContent(tnp_total_u2->GetBin(findBptbin(Bpt),findBybin(By))))));
			}
		    }
		}
	      
	  
	    }
	}

      TH1F* invEffave = new TH1F(Form("invEffave%d",trk),Form("invEffave%d",trk),nBins,ptBins);
      TH1F* invEffave_d = new TH1F(Form("invEffave_d%d",trk),Form("invEffave_d%d",trk),nBins,ptBins);
      TH1F* invEffave_u = new TH1F(Form("invEffave_u%d",trk),Form("invEffave_u%d",trk),nBins,ptBins);

      invEffave->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      invEffave->GetYaxis()->SetTitle("<1/(#alpha x #epsilon)>");
      invEffave->GetXaxis()->CenterTitle();
      invEffave->GetYaxis()->CenterTitle();
      invEffave->SetTitleOffset(1.4);
      
      for(int i=0;i<nBins;i++)
	{
	  //std::cout<<"sum: "<<sum[i]<<" nevts: "<<nevts[i]<<std::endl;
	  
	  sum[i]=sum[i]/nevts[i];
	  sumerr[i]=sqrt(sumerr[i])/nevts[i];
	  invEffave->SetBinContent(i+1,sum[i]);
	  invEffave->SetBinError(i+1,sumerr[i]);
	  
	  if(useTNP)
	    {
	      sum_d[i]=sum_d[i]/nevts[i];
	      sumerr_d[i]=sqrt(sumerr_d[i])/nevts[i];
	      sum_u[i]=sum_u[i]/nevts[i];
	      sumerr_u[i]=sqrt(sumerr_u[i])/nevts[i];
	      
	      invEffave_d->SetBinContent(i+1,sum_u[i]);
	      invEffave_d->SetBinError(i+1,sumerr_u[i]);
	      invEffave_u->SetBinContent(i+1,sum_d[i]);
	      invEffave_u->SetBinError(i+1,sumerr_d[i]);
	    }
	  
	  /*
	  std::cout<<"pt "<<ptBins[i]<<"-"<<ptBins[i+1]<<" Cent "<<centMin<<"-"<<centMax<<"%"<<" Nevts: "<<nevts[i]<<", Average 1/(acc*eff) (Cent): "<<sum[i]<<" pm "<<sumerr[i]<<std::endl;
	  
	  if(useTNP)
	    {
	      std::cout<<"pt "<<ptBins[i]<<"-"<<ptBins[i+1]<<" Cent "<<centMin<<"-"<<centMax<<"%"<<" Nevts: "<<nevts[i]<<", Average 1/(acc*eff) (Low): "<<sum_u[i]<<" pm "<<sumerr_u[i]<<", TNP syst (Low): "<<((sum[i]-sum_u[i])/sum[i])*100<<"%"<<std::endl;
	      std::cout<<"pt "<<ptBins[i]<<"-"<<ptBins[i+1]<<" Cent "<<centMin<<"-"<<centMax<<"%"<<" Nevts: "<<nevts[i]<<", Average 1/(acc*eff) (High): "<<sum_d[i]<<" pm "<<sumerr_d[i]<<", TNP syst (High): "<<((sum_d[i]-sum[i])/sum[i])*100<<"%"<<std::endl;
	    }
	  */
	}
      
      
      TFile* nominalfile = new TFile(nominaleffname.Data());
      TH1F* invEffnominal = (TH1F*)nominalfile->Get("invEffave");
      
      invEffnominal->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      invEffnominal->GetYaxis()->SetTitle("<1/(#alpha x #epsilon)>");
      invEffnominal->GetXaxis()->CenterTitle();
      invEffnominal->GetYaxis()->CenterTitle();
      invEffnominal->SetTitleOffset(1.4);
      
      invEffnominal->SetMaximum(30);
      invEffnominal->SetMinimum(0);
      
      TCanvas* cinvcompare = new TCanvas("","",600,600);
      cinvcompare->cd();
      if(uselogy) cinvcompare->SetLogy();

      TLegend *leg = new TLegend(0.70,0.70,0.90,0.80,NULL,"brNDC");
      leg->SetBorderSize(0);
      leg->SetTextSize(0.04);
      leg->SetTextFont(42);
      leg->SetFillStyle(0);
      
      std::cout<<""<<std::endl;
      for(int i=0;i<nBins;i++)
	{
	  std::cout<<"pt "<<ptBins[i]<<"-"<<ptBins[i+1]<<" Cent "<<centMin<<"-"<<centMax<<"% "<<trkvarname<<" systematics: "<<100*(sum[i]-invEffnominal->GetBinContent(i+1))/invEffnominal->GetBinContent(i+1)<<"%"<<std::endl;
	  if(TMath::Abs(100*(sum[i]-invEffnominal->GetBinContent(i+1))/invEffnominal->GetBinContent(i+1))>maxdev[i]) maxdev[i] = 100*(sum[i]-invEffnominal->GetBinContent(i+1))/invEffnominal->GetBinContent(i+1) ;
	}

      invEffnominal->SetLineColor(kBlack);
      invEffave->SetLineColor(kRed);
      invEffnominal->Draw();
      invEffave->Draw("same");

      leg->AddEntry(invEffnominal,"nominal","l");
      leg->AddEntry(invEffave,Form("%s",trkvarname.Data()),"l");
      leg->Draw("same");

      cinvcompare->SaveAs(Form("%s_%s.png",plotname.Data(),trkvarname.Data()));
      cinvcompare->SaveAs(Form("%s_%s.pdf",plotname.Data(),trkvarname.Data()));

      TFile* outf = new TFile(outputeffname.Data(),"recreate");
      outf->cd();
      invEffave->Write();
      if(useTNP)
	{
	  invEffave_d->Write();
	  invEffave_u->Write();
	} 
      outf->Close(); 
    }
  
  for(int i=0;i<nBins;i++)
    {
      std::cout<<"pt "<<ptBins[i]<<"-"<<ptBins[i+1]<<" Cent "<<centMin<<"-"<<centMax<<"%"<<" Track variable systematics(Max deviation): "<<maxdev[i]<<"%"<<std::endl;
    }
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
