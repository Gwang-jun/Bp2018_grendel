#include "uti.h"
#include "parameters.h"
#include "TRandom.h"

TString inputxsecname = "ppref_extrapolation.root";
TString inputmcname = "/raid5/data/gwangjun/crab_Bfinder_20190624_Hydjet_Pythia8_Official_BuToJpsiK_1033p1_pt3tkpt0p7dls2_allpthat_pthatweight_BDT.root";
TString inputeffname = "plottoyMC/MCstudiesPbPb.root";
TString outputname = "toyMC.root";

TString weightGpt = "1";
TString weightBgenpt = "1";
TString weighpthat = "pthatweight";
TString weightHiBin = "Ncoll";
TString weightPVz = "(TMath::Gaus(PVz,0.427450,4.873825)/(sqrt(2*3.14159)*4.873825))/(TMath::Gaus(PVz,0.909938,4.970989)/(sqrt(2*3.14159)*4.970989))";

TString selmcgen = "TMath::Abs(Gy)<2.4 && TMath::Abs(GpdgId)==521 && GisSignal==1 && GcollisionId==0";
TString selmcgenacceptance = "TMath::Abs(Gy)<2.4 && TMath::Abs(GpdgId)==521 && GisSignal==1 && GcollisionId==0 && Gtk1pt>0.7 && TMath::Abs(Gtk1eta)<2.4 && ((TMath::Abs(Gmu1eta)<1.2 && Gmu1pt>3.5) || (TMath::Abs(Gmu1eta)>1.2 && TMath::Abs(Gmu1eta)<2.1 && Gmu1pt>5.47-1.89*TMath::Abs(Gmu1eta)) || (TMath::Abs(Gmu1eta)>2.1 && TMath::Abs(Gmu2eta)<2.4 && Gmu1pt>1.5)) && ((TMath::Abs(Gmu2eta)<1.2 && Gmu2pt>3.5) || (TMath::Abs(Gmu2eta)>1.2 && TMath::Abs(Gmu2eta)<2.1 && Gmu2pt>5.47-1.89*TMath::Abs(Gmu2eta)) || (TMath::Abs(Gmu2eta)>2.1 && TMath::Abs(Gmu2eta)<2.4 && Gmu2pt>1.5))";
TString cut = "pprimaryVertexFilter && phfCoincFilter2Th4 && pclusterCompatibilityFilter && Btrk1Pt>0.9 && Bpt>5.0 && (BsvpvDistance/BsvpvDisErr)>2.0 && Bchi2cl>0.05 && TMath::Abs(Btrk1Eta)<2.4 && TMath::Abs(By)<2.4 && TMath::Abs(PVz)<15 && Bmass>5 && Bmass<6 && TMath::Abs(Bmumumass-3.096900)<0.15 && Bmu1SoftMuID && Bmu2SoftMuID && ((TMath::Abs(Bmu1eta)<1.2 && Bmu1pt>3.5) || (TMath::Abs(Bmu1eta)>1.2 && TMath::Abs(Bmu1eta)<2.1 && Bmu1pt>5.47-1.89*TMath::Abs(Bmu1eta)) || (TMath::Abs(Bmu1eta)>2.1 && TMath::Abs(Bmu1eta)<2.4 && Bmu1pt>1.5)) && ((TMath::Abs(Bmu2eta)<1.2 && Bmu2pt>3.5) || (TMath::Abs(Bmu2eta)>1.2 && TMath::Abs(Bmu2eta)<2.1 && Bmu2pt>5.47-1.89*TMath::Abs(Bmu2eta)) || (TMath::Abs(Bmu2eta)>2.1 && TMath::Abs(Bmu2eta)<2.4 && Bmu2pt>1.5)) && Bmu1isTriggered && Bmu2isTriggered && (Btrk1PixelHit+Btrk1StripHit)>=11 && (Btrk1Chi2ndf/(Btrk1nStripLayer+Btrk1nPixelLayer))<0.18 && TMath::Abs(Btrk1PtErr/Btrk1Pt)<0.1 && ((Bpt>5 && Bpt<7 && (BsvpvDistance/BsvpvDisErr)>12.0 && cos(Bdtheta)>0.95 && BDT_5_7>0.07) || (Bpt>7 && Bpt<10 && BDT_7_10>0.08) || (Bpt>10 && Bpt<15 && BDT_10_15>0.09) || (Bpt>15 && Bpt<20 && BDT_15_20>0.09) || (Bpt>20 && Bpt<30 && BDT_20_30>0.07) || (Bpt>30 && Bpt<50 && BDT_30_50>0.12) || (Bpt>50 && Bpt<100 && BDT_50_100>0.24))";

void toyMC()
{
  gStyle->SetOptStat(0);

  float hiBinMin = 0.;
  float hiBinMax = 180.;
  selmcgen = selmcgen+Form("&&hiBin>=%f&&hiBin<=%f",hiBinMin,hiBinMax);
  selmcgenacceptance = selmcgenacceptance+Form("&&hiBin>=%f&&hiBin<=%f",hiBinMin,hiBinMax);
  cut=cut+Form("&&hiBin>=%f&&hiBin<=%f",hiBinMin,hiBinMax);
  
  TFile* inputxsec = new TFile(inputxsecname.Data());
  TH1D* nominal = (TH1D*) inputxsec->Get("real_ref_new_sym_record");
  
  TFile* inputMC = new TFile(inputmcname.Data());
  TTree* ntMC = (TTree*) inputMC->Get("Bfinder/ntKp");
  ntMC->AddFriend("hltanalysis/HltTree");
  ntMC->AddFriend("hiEvtAnalyzer/HiTree");
  ntMC->AddFriend("Bfinder/ntGen");
  ntMC->AddFriend("skimanalysis/HltTree");
  ntMC->AddFriend("BDT");
  TTree* ntGen = (TTree*) inputMC->Get("Bfinder/ntGen");
  ntGen->AddFriend("Bfinder/ntKp");
  ntGen->AddFriend("hiEvtAnalyzer/HiTree");
  ntGen->AddFriend("hltanalysis/HltTree");
  ntGen->AddFriend("skimanalysis/HltTree");
  ntGen->AddFriend("BDT");
  TH1D* Genpt = new TH1D("Genpt","",nBins,ptBins);
  ntGen->Project("Genpt","Gpt",TCut(weighpthat)*TCut(selmcgen));
  Genpt->GetXaxis()->SetTitle("p_{t} (GeV)");
  Genpt->GetYaxis()->SetTitle("d#sigma/dp_{t} (pb/GeV)");
  Genpt->GetYaxis()->SetTitleOffset(0.9);
  divideBinWidth(Genpt);
  Genpt->Scale(1.0/Genpt->Integral("width"));

  TF1* fGen = new TF1("fGen","TMath::Exp([0]+[1]*x+[2]*x*x+[3]*x*x*x)",5,100);
  fGen->SetParameters(0,-1,0,0);
  fGen->SetParLimits(0,-100,100);
  fGen->SetParLimits(1,-100,100);
  fGen->SetParLimits(2,-100,100);
  fGen->SetParLimits(3,-100,100);
  Genpt->Fit(fGen,"R");
  
  TCanvas* cGen = new TCanvas("cGen","",600,600);                                                                                  
  cGen->cd();                                                                                                                                  
  Genpt->Draw();
  fGen->Draw("same");
  cGen->SetLogy();
  cGen->SaveAs("plottoyMC/toyMC_Gen.png");
  cGen->SaveAs("plottoyMC/toyMC_Gen.pdf");

  TFile* output = new TFile(outputname.Data(),"recreate");

  TFile* inputeff = new TFile(inputeffname.Data());
  TH1D* hEffnominal = (TH1D*) inputeff->Get("hEff");
  
  TH1D* htoy1 = new TH1D("htoy1","",100,hEffnominal->GetBinContent(1)-5*hEffnominal->GetBinError(1),hEffnominal->GetBinContent(1)+5*hEffnominal->GetBinError(1));
  TH1D* htoy2 = new TH1D("htoy2","",100,hEffnominal->GetBinContent(2)-5*hEffnominal->GetBinError(2),hEffnominal->GetBinContent(2)+5*hEffnominal->GetBinError(2));
  TH1D* htoy3 = new TH1D("htoy3","",100,hEffnominal->GetBinContent(3)-5*hEffnominal->GetBinError(3),hEffnominal->GetBinContent(3)+5*hEffnominal->GetBinError(3));
  TH1D* htoy4 = new TH1D("htoy4","",100,hEffnominal->GetBinContent(4)-5*hEffnominal->GetBinError(4),hEffnominal->GetBinContent(4)+5*hEffnominal->GetBinError(4));
  TH1D* htoy5 = new TH1D("htoy5","",100,hEffnominal->GetBinContent(5)-5*hEffnominal->GetBinError(5),hEffnominal->GetBinContent(5)+5*hEffnominal->GetBinError(5));
  TH1D* htoy6 = new TH1D("htoy6","",100,hEffnominal->GetBinContent(6)-5*hEffnominal->GetBinError(6),hEffnominal->GetBinContent(6)+5*hEffnominal->GetBinError(6));
  TH1D* htoy7 = new TH1D("htoy7","",100,hEffnominal->GetBinContent(7)-5*hEffnominal->GetBinError(7),hEffnominal->GetBinContent(7)+5*hEffnominal->GetBinError(7));
  
  for(int i=0;i<1;i++)
    {
      static Int_t count=0;
      count++;
      TH1D* h = new TH1D(Form("h%d",count),"",nBins,ptBins);
      h->GetXaxis()->SetTitle("p_{t} (GeV)");
      h->GetYaxis()->SetTitle("d#sigma/dp_{t} (pb/GeV)");
      h->GetXaxis()->CenterTitle();
      h->GetYaxis()->CenterTitle();
      h->GetYaxis()->SetTitleOffset(0.9);

      for(int j=0;j<nBins;j++)
	{
	  h->SetBinContent(j+1,gRandom->Gaus(nominal->GetBinContent(j+1),nominal->GetBinError(j+1)));
	  h->SetBinError(j+1,nominal->GetBinError(j+1));
	  //std::cout<<h->GetBinContent(j+1)<<std::endl;
	}

      h->Sumw2();
      h->Scale(1.0/h->Integral("width"));

      TF1* f = new TF1(Form("f%d",count),"TMath::Exp([0]+[1]*x+[2]*x*x+[3]*x*x*x)",5,100);
      f->SetParameters(fGen->GetParameter(0),fGen->GetParameter(1),fGen->GetParameter(2),fGen->GetParameter(3));
      f->SetParLimits(0,-100,100);
      f->SetParLimits(1,-100,100);
      f->SetParLimits(2,-100,100);
      f->SetParLimits(3,-100,100);
      h->Fit(f,"R");

      TCanvas* cxfit = new TCanvas(Form("cxfit%d",count),"",600,600);
      cxfit->cd();
      h->Draw();
      f->Draw("same");
      cxfit->SetLogy();
      cxfit->SaveAs(Form("plottoyMC/toyMC_xsecfit_%d.png",count));
      cxfit->SaveAs(Form("plottoyMC/toyMC_xsecfit_%d.pdf",count));

      /*
      h->Divide(Genpt);
      TF1 *f = new TF1(Form("f%d",count),"([0]+[1]*x+[2]*x*x)*TMath::Exp(-[3]*x)+[4]",5,100);      
      f->SetParLimits(0,-100,100);
      f->SetParLimits(1,-10,10);
      f->SetParLimits(2,-10,10);
      f->SetParLimits(3,0,10);
      f->SetParLimits(4,0,2);
      h->Fit(f,"IRQ");

      TCanvas* cratio = new TCanvas(Form("cratio%d",count),"",600,600);
      cratio->cd();
      h->Draw();
      f->Draw("same");
      cratio->SaveAs(Form("plottoyMC/toyMC_ptweight_%d.png",count));

      weightGpt = Form("(%f+%f*Gpt+%f*Gpt*Gpt)*TMath::Exp(-%f*Gpt)+%f",f->GetParameter(0),f->GetParameter(1),f->GetParameter(2),f->GetParameter(3),f->GetParameter(4));
      weightBgenpt = Form("(%f+%f*Bgenpt+%f*Bgenpt*Bgenpt)*TMath::Exp(-%f*Bgenpt)+%f",f->GetParameter(0),f->GetParameter(1),f->GetParameter(2),f->GetParameter(3),f->GetParameter(4));
      */
      
      //weightGpt = Form("(TMath::Exp(%f+%f*Gpt+%f*Gpt*Gpt+%f*Gpt*Gpt*Gpt))/(TMath::Exp(%f+%f*Gpt+%f*Gpt*Gpt+%f*Gpt*Gpt*Gpt))",f->GetParameter(0),f->GetParameter(1),f->GetParameter(2),f->GetParameter(3),fGen->GetParameter(0),fGen->GetParameter(1),fGen->GetParameter(2),fGen->GetParameter(3));
      //weightBgenpt = Form("(TMath::Exp(%f+%f*Bgenpt+%f*Bgenpt*Bgenpt+%f*Bgenpt*Bgenpt*Bgenpt))/(TMath::Exp(%f+%f*Bgenpt+%f*Bgenpt*Bgenpt+%f*Bgenpt*Bgenpt*Bgenpt))",f->GetParameter(0),f->GetParameter(1),f->GetParameter(2),f->GetParameter(3),fGen->GetParameter(0),fGen->GetParameter(1),fGen->GetParameter(2),fGen->GetParameter(3));

      //weightGpt = Form("(TMath::Exp((%f-%f)+(%f-%f)*Gpt+(%f-%f)*Gpt*Gpt+(%f-%f)*Gpt*Gpt*Gpt))",f->GetParameter(0),fGen->GetParameter(0),f->GetParameter(1),fGen->GetParameter(1),f->GetParameter(2),fGen->GetParameter(2),f->GetParameter(3),fGen->GetParameter(3));
      //weightBgenpt = Form("(TMath::Exp((%f-%f)+(%f-%f)*Bgenpt+(%f-%f)*Bgenpt*Bgenpt+(%f-%f)*Bgenpt*Bgenpt*Bgenpt))",f->GetParameter(0),fGen->GetParameter(0),f->GetParameter(1),fGen->GetParameter(1),f->GetParameter(2),fGen->GetParameter(2),f->GetParameter(3),fGen->GetParameter(3));

      weightGpt = "TMath::Exp(-2.981465+0.653199*Gpt-0.043233*Gpt*Gpt+0.001103*Gpt*Gpt*Gpt)";
      weightBgenpt = "TMath::Exp(-2.981465+0.653199*Bgenpt-0.043233*Bgenpt*Bgenpt+0.001103*Bgenpt*Bgenpt*Bgenpt)";
      
      std::cout<<weightGpt<<std::endl;
      
      TH1D* hPtMC = new TH1D(Form("hPtMC%d",count),"",nBins,ptBins);
      TH1D* hPtGen = new TH1D(Form("hPtGen%d",count),"",nBins,ptBins);
      TH1D* hPtGenAcc = new TH1D(Form("hPtGenAcc%d",count),"",nBins,ptBins);
      TH1D* hPtGenAccWeighted = new TH1D(Form("hPtGenAccWeighted%d",count),"",nBins,ptBins);

      ntMC->Project(Form("hPtMC%d",count),"Bpt",TCut(weighpthat)*TCut(weightBgenpt)*TCut(weightHiBin)*TCut(weightPVz)*(TCut(cut.Data())&&"(Bgen==23333)"));
      ntGen->Project(Form("hPtGen%d",count),"Gpt",TCut(weighpthat)*TCut(weightGpt)*(TCut(selmcgen.Data())));
      ntGen->Project(Form("hPtGenAcc%d",count),"Gpt",TCut(weighpthat)*TCut(weightGpt)*(TCut(selmcgenacceptance.Data())));
      ntGen->Project(Form("hPtGenAccWeighted%d",count),"Gpt",TCut(weighpthat)*TCut(weightGpt)*TCut(weightHiBin)*TCut(weightPVz)*(TCut(selmcgenacceptance.Data())));

      output->cd();
      hPtMC->Write();
      hPtGen->Write();
      hPtGenAcc->Write();
      hPtGenAccWeighted->Write();

      /*
      double sf_pbpb[7] = {1.0499, 1.1294, 1.1292, 1.0927, 1.0610, 1.0203, 1.0005};
      for(int j=0;j<nBins;j++)
	{
	  hPtMC->SetBinContent(i+1, hPtMC->GetBinContent(i+1)*sf_pbpb[i]);
	}
      */

      /*
      divideBinWidth(hPtMC);
      divideBinWidth(hPtGen);
      divideBinWidth(hPtGenAcc);
      divideBinWidth(hPtGenAccWeighted);
      
      TH1D* hEfftoy = new TH1D(Form("hEfftoy%d",count),"",nBins,ptBins);
      hEfftoy->Divide(hPtGenAcc,hPtGen,1,1,"");
      hEfftoy->Multiply(hPtMC);
      hEfftoy->Divide(hPtGenAccWeighted);
      
      htoy1->Fill(hEfftoy->GetBinContent(1));
      htoy2->Fill(hEfftoy->GetBinContent(2));
      htoy3->Fill(hEfftoy->GetBinContent(3));
      htoy4->Fill(hEfftoy->GetBinContent(4));
      htoy5->Fill(hEfftoy->GetBinContent(5));
      htoy6->Fill(hEfftoy->GetBinContent(6));
      htoy7->Fill(hEfftoy->GetBinContent(7));
      
      output->cd();
      hPtMC->Write();
      hPtGen->Write();
      hPtGenAcc->Write();
      hPtGenAccWeighted->Write();
      hEfftoy->Write();            
      */
    }
  /*
  htoy1->Write();
  htoy2->Write();
  htoy3->Write();
  htoy4->Write();
  htoy5->Write();
  htoy6->Write();
  htoy7->Write();

  TF1* ftoy1 = new TF1("ftoy1","TMath::Gaus(x,[0],[1])/(sqrt(2*3.14)*[1])",0,1);
  ftoy1->SetParLimits(0,0,1);
  ftoy1->SetParLimits(1,0,1);
  htoy1->Fit(ftoy1,"R");
  std::cout<<"Mean:"<<ftoy1->GetParameter(0)<<" width:"<<ftoy1->GetParameter(1)<<std::endl;
  std::cout<<"Deviation:"<<TMath::Abs(ftoy1->GetParameter(0)-hEffnominal->GetBinContent(1))/hEffnominal->GetBinContent(1)<<std::endl;

  TCanvas* ctoy1 = new TCanvas("ctoy1","",600,600);
  ctoy1->cd();
  htoy1->Draw();
  ftoy1->Draw("same");
  ctoy1->SaveAs("plottoyMC/toyMC_eff_1.png");
  ctoy1->SaveAs("plottoyMC/toyMC_eff_1.pdf");

  TF1* ftoy2 = new TF1("ftoy2","TMath::Gaus(x,[0],[1])/(sqrt(2*3.14)*[1])",0,1);
  ftoy2->SetParLimits(0,0,1);
  ftoy2->SetParLimits(1,0,1);
  htoy2->Fit(ftoy2,"R");
  std::cout<<"Mean:"<<ftoy2->GetParameter(0)<<" width:"<<ftoy2->GetParameter(1)<<std::endl;
  std::cout<<"Deviation:"<<TMath::Abs(ftoy2->GetParameter(0)-hEffnominal->GetBinContent(2))/hEffnominal->GetBinContent(2)<<std::endl;
  
  TCanvas* ctoy2 = new TCanvas("ctoy2","",600,600);
  ctoy2->cd();
  htoy2->Draw();
  ftoy2->Draw("same");
  ctoy2->SaveAs("plottoyMC/toyMC_eff_2.png");
  ctoy2->SaveAs("plottoyMC/toyMC_eff_2.pdf");

  TF1* ftoy3 = new TF1("ftoy3","TMath::Gaus(x,[0],[1])/(sqrt(2*3.14)*[1])",0,1);
  ftoy3->SetParLimits(0,0,1);
  ftoy3->SetParLimits(1,0,1);
  htoy3->Fit(ftoy3,"R");
  std::cout<<"Mean:"<<ftoy3->GetParameter(0)<<" width:"<<ftoy3->GetParameter(1)<<std::endl;
  std::cout<<"Deviation:"<<TMath::Abs(ftoy3->GetParameter(0)-hEffnominal->GetBinContent(3))/hEffnominal->GetBinContent(3)<<std::endl;

  TCanvas* ctoy3 = new TCanvas("ctoy3","",600,600);
  ctoy3->cd();
  htoy3->Draw();
  ftoy3->Draw("same");
  ctoy3->SaveAs("plottoyMC/toyMC_eff_3.png");
  ctoy3->SaveAs("plottoyMC/toyMC_eff_3.pdf");

  TF1* ftoy4 = new TF1("ftoy4","TMath::Gaus(x,[0],[1])/(sqrt(2*3.14)*[1])",0,1);
  ftoy4->SetParLimits(0,0,1);
  ftoy4->SetParLimits(1,0,1);
  htoy4->Fit(ftoy4,"R");
  std::cout<<"Mean:"<<ftoy4->GetParameter(0)<<" width:"<<ftoy4->GetParameter(1)<<std::endl;
  std::cout<<"Deviation:"<<TMath::Abs(ftoy4->GetParameter(0)-hEffnominal->GetBinContent(4))/hEffnominal->GetBinContent(4)<<std::endl;
  
  TCanvas* ctoy4 = new TCanvas("ctoy4","",600,600);
  ctoy4->cd();
  htoy4->Draw();
  ftoy4->Draw("same");
  ctoy4->SaveAs("plottoyMC/toyMC_eff_4.png");
  ctoy4->SaveAs("plottoyMC/toyMC_eff_4.pdf");

  TF1* ftoy5 = new TF1("ftoy5","TMath::Gaus(x,[0],[1])/(sqrt(2*3.14)*[1])",0,1);
  ftoy5->SetParLimits(0,0,1);
  ftoy5->SetParLimits(1,0,1);
  htoy5->Fit(ftoy5,"R");
  std::cout<<"Mean:"<<ftoy5->GetParameter(0)<<" width:"<<ftoy5->GetParameter(1)<<std::endl;
  std::cout<<"Deviation:"<<TMath::Abs(ftoy5->GetParameter(0)-hEffnominal->GetBinContent(5))/hEffnominal->GetBinContent(5)<<std::endl;  

  TCanvas* ctoy5 = new TCanvas("ctoy5","",600,600);
  ctoy5->cd();
  htoy5->Draw();
  ftoy5->Draw("same");
  ctoy5->SaveAs("plottoyMC/toyMC_eff_5.png");
  ctoy5->SaveAs("plottoyMC/toyMC_eff_5.pdf");

  TF1* ftoy6 = new TF1("ftoy6","TMath::Gaus(x,[0],[1])/(sqrt(2*3.14)*[1])",0,1);
  ftoy6->SetParLimits(0,0,1);
  ftoy6->SetParLimits(1,0,1);
  htoy6->Fit(ftoy6,"R");
  std::cout<<"Mean:"<<ftoy6->GetParameter(0)<<" width:"<<ftoy6->GetParameter(1)<<std::endl;
  std::cout<<"Deviation:"<<TMath::Abs(ftoy6->GetParameter(0)-hEffnominal->GetBinContent(6))/hEffnominal->GetBinContent(6)<<std::endl;  

  TCanvas* ctoy6 = new TCanvas("ctoy6","",600,600);
  ctoy6->cd();
  htoy6->Draw();
  ftoy6->Draw("same");
  ctoy6->SaveAs("plottoyMC/toyMC_eff_6.png");
  ctoy6->SaveAs("plottoyMC/toyMC_eff_6.pdf");

  TF1* ftoy7 = new TF1("ftoy7","TMath::Gaus(x,[0],[1])/(sqrt(2*3.14)*[1])",0,1);
  ftoy7->SetParLimits(0,0,1);
  ftoy7->SetParLimits(1,0,1);
  htoy7->Fit(ftoy7,"R");
  std::cout<<"Mean:"<<ftoy7->GetParameter(0)<<" width:"<<ftoy7->GetParameter(1)<<std::endl;
  std::cout<<"Deviation:"<<TMath::Abs(ftoy7->GetParameter(0)-hEffnominal->GetBinContent(7))/hEffnominal->GetBinContent(7)<<std::endl;  

  TCanvas* ctoy7 = new TCanvas("ctoy7","",600,600);
  ctoy7->cd();
  htoy7->Draw();
  ftoy7->Draw("same");
  ctoy7->SaveAs("plottoyMC/toyMC_eff_7.png");
  ctoy7->SaveAs("plottoyMC/toyMC_eff_7.pdf");
  */

  output->Close();  
  
  return;
}

int main()
{
  toyMC();
  return 0;
}
