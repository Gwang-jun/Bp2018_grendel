void canvasSigmaBplusRatiopp()
{
//=========Macro generated from canvas: cSigma/
//=========  (Sun Jul 21 15:17:33 2019) by ROOT version 6.14/04
   TCanvas *cSigma = new TCanvas("cSigma", "",0,0,600,600);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   cSigma->SetHighLightColor(2);
   cSigma->Range(-1.989924,-0.2917772,25.49622,2.212202);
   cSigma->SetFillColor(0);
   cSigma->SetBorderMode(0);
   cSigma->SetBorderSize(2);
   cSigma->SetLogy();
   cSigma->SetLeftMargin(0.1451613);
   cSigma->SetRightMargin(0.05443548);
   cSigma->SetTopMargin(0.08474576);
   cSigma->SetBottomMargin(0);
   cSigma->SetFrameBorderMode(0);
  
// ------------>Primitives in pad: pSigma
   TPad *pSigma = new TPad("pSigma", "",0,0.3,1,1);
   pSigma->Draw();
   pSigma->cd();
   pSigma->Range(-19.04282,2,112.1411,8.555556);
   pSigma->SetFillColor(0);
   pSigma->SetBorderMode(0);
   pSigma->SetBorderSize(2);
   pSigma->SetLogy();
   pSigma->SetLeftMargin(0.1451613);
   pSigma->SetRightMargin(0.05443548);
   pSigma->SetTopMargin(0.08474576);
   pSigma->SetBottomMargin(0);
   pSigma->SetFrameBorderMode(0);
   pSigma->SetFrameBorderMode(0);
   
   TH2F *hemptySigma__1 = new TH2F("hemptySigma__1","",50,0,105,10,100,1e+08);
   hemptySigma__1->SetMinimum(0);
   hemptySigma__1->SetMaximum(2);
   hemptySigma__1->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   hemptySigma__1->SetLineColor(ci);
   hemptySigma__1->SetMarkerStyle(20);
   hemptySigma__1->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   hemptySigma__1->GetXaxis()->CenterTitle(true);
   hemptySigma__1->GetXaxis()->SetLabelFont(42);
   hemptySigma__1->GetXaxis()->SetLabelOffset(0.0015);
   hemptySigma__1->GetXaxis()->SetLabelSize(0.036);
   hemptySigma__1->GetXaxis()->SetTitleSize(0.036);
   hemptySigma__1->GetXaxis()->SetTitleFont(42);
   hemptySigma__1->GetYaxis()->SetTitle("#frac{d#sigma}{dp_{T}} ( pb GeV^{-1}c)");
   hemptySigma__1->GetYaxis()->CenterTitle(true);
   hemptySigma__1->GetYaxis()->SetLabelFont(42);
   hemptySigma__1->GetYaxis()->SetLabelSize(0.06);
   hemptySigma__1->GetYaxis()->SetTitleSize(0.06);
   hemptySigma__1->GetYaxis()->SetTitleFont(42);
   hemptySigma__1->GetZaxis()->SetLabelFont(42);
   hemptySigma__1->GetZaxis()->SetLabelSize(0.035);
   hemptySigma__1->GetZaxis()->SetTitleSize(0.035);
   hemptySigma__1->GetZaxis()->SetTitleFont(42);
   hemptySigma__1->Draw("");
   
   Double_t gaeSigmaBplus_fx3001[7] = {
   6,
   8.5,
   12.5,
   17.5,
   25,
   40,
   75};
   Double_t gaeSigmaBplus_fy3001[7] = {
   4634232,
   2167553,
   685281.8,
   186414.7,
   44056.79,
   5247.477,
   248.2524};
   Double_t gaeSigmaBplus_felx3001[7] = {
   1,
   1.5,
   2.5,
   2.5,
   5,
   10,
   25};
   Double_t gaeSigmaBplus_fely3001[7] = {
   1279894,
   531575.1,
   157767.3,
   38660.19,
   8088.552,
   805.5551,
   30.76349};
   Double_t gaeSigmaBplus_fehx3001[7] = {
   1,
   1.5,
   2.5,
   2.5,
   5,
   10,
   25};
   Double_t gaeSigmaBplus_fehy3001[7] = {
   1908619,
   863802.8,
   246333.8,
   56879.75,
   11149.93,
   1010.418,
   34.60447};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(7,gaeSigmaBplus_fx3001,gaeSigmaBplus_fy3001,gaeSigmaBplus_felx3001,gaeSigmaBplus_fehx3001,gaeSigmaBplus_fely3001,gaeSigmaBplus_fehy3001);
   grae->SetName("gaeSigmaBplus");
   grae->SetTitle(";p_{T}(GeV/c);d#sigma/dp_{T} (B^{+}) (pb GeV-1c)");

   ci = TColor::GetColor("#ffcc00");
   grae->SetFillColor(ci);

   ci = TColor::GetColor("#ffcc00");
   grae->SetLineColor(ci);
   grae->SetLineWidth(3);
   
   TH1F *Graph_gaeSigmaBplus3001 = new TH1F("Graph_gaeSigmaBplus3001","",100,0,109.5);
   Graph_gaeSigmaBplus3001->SetMinimum(195.74);
   Graph_gaeSigmaBplus3001->SetMaximum(7197114);
   Graph_gaeSigmaBplus3001->SetDirectory(0);
   Graph_gaeSigmaBplus3001->SetStats(0);
   Graph_gaeSigmaBplus3001->GetXaxis()->SetTitle("p_{T}(GeV/c)");
   Graph_gaeSigmaBplus3001->GetYaxis()->SetTitle("d#sigma/dp_{T} (B^{+}) (pb GeV-1c)");
   grae->SetHistogram(Graph_gaeSigmaBplus3001);
   
   grae->Draw("5");
   Double_t xAxis1[8] = {5, 7, 10, 15, 20, 30, 50, 100}; 
   
   TH1D *hPtSigma__2 = new TH1D("hPtSigma__2","",7, xAxis1);
   hPtSigma__2->SetLineWidth(2);
   hPtSigma__2->SetMarkerStyle(20);
   hPtSigma__2->SetMarkerSize(1.2);
   hPtSigma__2->GetXaxis()->SetTitle("B^{+} p_{T} (GeV/c)");
   hPtSigma__2->GetXaxis()->SetLabelFont(42);
   hPtSigma__2->GetXaxis()->SetLabelSize(0.035);
   hPtSigma__2->GetXaxis()->SetTitleSize(0.035);
   hPtSigma__2->GetXaxis()->SetTitleFont(42);
   hPtSigma__2->GetYaxis()->SetTitle("Uncorrected dN(B^{+})/dp_{T}");
   hPtSigma__2->GetYaxis()->SetLabelFont(42);
   hPtSigma__2->GetYaxis()->SetLabelSize(0.035);
   hPtSigma__2->GetYaxis()->SetTitleSize(0.035);
   hPtSigma__2->GetYaxis()->SetTitleOffset(0);
   hPtSigma__2->GetYaxis()->SetTitleFont(42);
   hPtSigma__2->GetZaxis()->SetLabelFont(42);
   hPtSigma__2->GetZaxis()->SetLabelSize(0.035);
   hPtSigma__2->GetZaxis()->SetTitleSize(0.035);
   hPtSigma__2->GetZaxis()->SetTitleFont(42);
   hPtSigma__2->Draw("epsame");
   
   Double_t gaeCrossSyst_fx3002[7] = {
   6,
   8.5,
   12.5,
   17.5,
   25,
   40,
   75};
   Double_t gaeCrossSyst_fy3002[7] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t gaeCrossSyst_felx3002[7] = {
   1,
   1.5,
   2.5,
   2.5,
   5,
   10,
   25};
   Double_t gaeCrossSyst_fely3002[7] = {
   -nan,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t gaeCrossSyst_fehx3002[7] = {
   1,
   1.5,
   2.5,
   2.5,
   5,
   10,
   25};
   Double_t gaeCrossSyst_fehy3002[7] = {
   -nan,
   0,
   0,
   0,
   0,
   0,
   0};
   grae = new TGraphAsymmErrors(7,gaeCrossSyst_fx3002,gaeCrossSyst_fy3002,gaeCrossSyst_felx3002,gaeCrossSyst_fehx3002,gaeCrossSyst_fely3002,gaeCrossSyst_fehy3002);
   grae->SetName("gaeCrossSyst");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);
   grae->SetFillStyle(0);
   grae->SetLineWidth(2);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(0.8);
   
   TH1F *Graph_gaeCrossSyst3002 = new TH1F("Graph_gaeCrossSyst3002","Graph",100,0,109.5);
   Graph_gaeCrossSyst3002->SetMinimum(0.0011);
   Graph_gaeCrossSyst3002->SetMaximum(1.1);
   Graph_gaeCrossSyst3002->SetDirectory(0);
   Graph_gaeCrossSyst3002->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_gaeCrossSyst3002->SetLineColor(ci);
   Graph_gaeCrossSyst3002->SetMarkerStyle(20);
   Graph_gaeCrossSyst3002->GetXaxis()->SetLabelFont(42);
   Graph_gaeCrossSyst3002->GetXaxis()->SetLabelSize(0.035);
   Graph_gaeCrossSyst3002->GetXaxis()->SetTitleSize(0.035);
   Graph_gaeCrossSyst3002->GetXaxis()->SetTitleFont(42);
   Graph_gaeCrossSyst3002->GetYaxis()->SetLabelFont(42);
   Graph_gaeCrossSyst3002->GetYaxis()->SetLabelSize(0.035);
   Graph_gaeCrossSyst3002->GetYaxis()->SetTitleSize(0.035);
   Graph_gaeCrossSyst3002->GetYaxis()->SetTitleOffset(0);
   Graph_gaeCrossSyst3002->GetYaxis()->SetTitleFont(42);
   Graph_gaeCrossSyst3002->GetZaxis()->SetLabelFont(42);
   Graph_gaeCrossSyst3002->GetZaxis()->SetLabelSize(0.035);
   Graph_gaeCrossSyst3002->GetZaxis()->SetTitleSize(0.035);
   Graph_gaeCrossSyst3002->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_gaeCrossSyst3002);
   
   grae->Draw("5");
   TLatex *   tex = new TLatex(0.18,0.87,"CMS");
tex->SetNDC();
   tex->SetTextAlign(13);
   tex->SetTextSize(0.08);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.9,0.936,"28.0 pb^{-1} (pp 5.02 TeV)");
tex->SetNDC();
   tex->SetTextAlign(31);
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.94,0.95,"pp #sqrt{s_{NN}} = 5.02 TeV");
tex->SetNDC();
   tex->SetTextAlign(32);
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.53,0.65,"|y| < 2.4");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.77,0.82,"B^{+}");
tex->SetNDC();
   tex->SetTextSize(0.09);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.53,0.59,"Global uncert. 3.4%");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->Draw();
   
   TLegend *leg = new TLegend(0.52,0.7,0.85,0.8,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("hPtSigma","Data","pf");
   entry->SetFillStyle(1001);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("gaeSigmaBplus","FONLL","f");

   ci = TColor::GetColor("#ffcc00");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);

   ci = TColor::GetColor("#ffcc00");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   leg->Draw();
   pSigma->Modified();
   cSigma->cd();
  
// ------------>Primitives in pad: pRatio
   TPad *pRatio = new TPad("pRatio", "",0,0,1,0.3);
   pRatio->Draw();
   pRatio->cd();
   pRatio->Range(-19.04282,-0.4857143,112.1411,1.8);
   pRatio->SetFillColor(0);
   pRatio->SetBorderMode(0);
   pRatio->SetBorderSize(2);
   pRatio->SetLeftMargin(0.1451613);
   pRatio->SetRightMargin(0.05443548);
   pRatio->SetTopMargin(0);
   pRatio->SetBottomMargin(0.3);
   pRatio->SetFrameBorderMode(0);
   pRatio->SetFrameBorderMode(0);
   
   TH2F *hemptyRatio__3 = new TH2F("hemptyRatio__3","",50,0,105,10,0.2,1.8);
   hemptyRatio__3->SetStats(0);

   ci = TColor::GetColor("#000099");
   hemptyRatio__3->SetLineColor(ci);
   hemptyRatio__3->SetMarkerStyle(20);
   hemptyRatio__3->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   hemptyRatio__3->GetXaxis()->CenterTitle(true);
   hemptyRatio__3->GetXaxis()->SetLabelFont(42);
   hemptyRatio__3->GetXaxis()->SetLabelSize(0.12);
   hemptyRatio__3->GetXaxis()->SetTitleSize(0.12);
   hemptyRatio__3->GetXaxis()->SetTickLength(0.09999999);
   hemptyRatio__3->GetXaxis()->SetTitleFont(42);
   hemptyRatio__3->GetYaxis()->SetTitle("Data / FONLL");
   hemptyRatio__3->GetYaxis()->CenterTitle(true);
   hemptyRatio__3->GetYaxis()->SetNdivisions(505);
   hemptyRatio__3->GetYaxis()->SetLabelFont(42);
   hemptyRatio__3->GetYaxis()->SetLabelSize(0.12);
   hemptyRatio__3->GetYaxis()->SetTitleSize(0.12);
   hemptyRatio__3->GetYaxis()->SetTitleOffset(0.5);
   hemptyRatio__3->GetYaxis()->SetTitleFont(42);
   hemptyRatio__3->GetZaxis()->SetLabelFont(42);
   hemptyRatio__3->GetZaxis()->SetLabelSize(0.035);
   hemptyRatio__3->GetZaxis()->SetTitleSize(0.035);
   hemptyRatio__3->GetZaxis()->SetTitleFont(42);
   hemptyRatio__3->Draw("");
   TLine *line = new TLine(0,1,105,1);
   line->SetLineStyle(2);
   line->Draw();
   
   Double_t gaeRatioCrossFONLLunity_fx3003[7] = {
   6,
   8.5,
   12.5,
   17.5,
   25,
   40,
   75};
   Double_t gaeRatioCrossFONLLunity_fy3003[7] = {
   1,
   1,
   1,
   1,
   1,
   1,
   1};
   Double_t gaeRatioCrossFONLLunity_felx3003[7] = {
   1,
   1.5,
   2.5,
   2.5,
   5,
   10,
   25};
   Double_t gaeRatioCrossFONLLunity_fely3003[7] = {
   0.2761826,
   0.245242,
   0.2302225,
   0.2073881,
   0.1835938,
   0.1535128,
   0.1239202};
   Double_t gaeRatioCrossFONLLunity_fehx3003[7] = {
   1,
   1.5,
   2.5,
   2.5,
   5,
   10,
   25};
   Double_t gaeRatioCrossFONLLunity_fehy3003[7] = {
   0.4118522,
   0.3985152,
   0.3594636,
   0.3051248,
   0.2530808,
   0.1925531,
   0.1393923};
   grae = new TGraphAsymmErrors(7,gaeRatioCrossFONLLunity_fx3003,gaeRatioCrossFONLLunity_fy3003,gaeRatioCrossFONLLunity_felx3003,gaeRatioCrossFONLLunity_fehx3003,gaeRatioCrossFONLLunity_fely3003,gaeRatioCrossFONLLunity_fehy3003);
   grae->SetName("gaeRatioCrossFONLLunity");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#ffcc00");
   grae->SetFillColor(ci);

   ci = TColor::GetColor("#ffcc00");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_gaeRatioCrossFONLLunity3003 = new TH1F("Graph_gaeRatioCrossFONLLunity3003","Graph",100,0,109.5);
   Graph_gaeRatioCrossFONLLunity3003->SetMinimum(0.6550139);
   Graph_gaeRatioCrossFONLLunity3003->SetMaximum(1.480656);
   Graph_gaeRatioCrossFONLLunity3003->SetDirectory(0);
   Graph_gaeRatioCrossFONLLunity3003->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_gaeRatioCrossFONLLunity3003->SetLineColor(ci);
   Graph_gaeRatioCrossFONLLunity3003->SetMarkerStyle(20);
   Graph_gaeRatioCrossFONLLunity3003->GetXaxis()->SetLabelFont(42);
   Graph_gaeRatioCrossFONLLunity3003->GetXaxis()->SetLabelSize(0.035);
   Graph_gaeRatioCrossFONLLunity3003->GetXaxis()->SetTitleSize(0.035);
   Graph_gaeRatioCrossFONLLunity3003->GetXaxis()->SetTitleFont(42);
   Graph_gaeRatioCrossFONLLunity3003->GetYaxis()->SetLabelFont(42);
   Graph_gaeRatioCrossFONLLunity3003->GetYaxis()->SetLabelSize(0.035);
   Graph_gaeRatioCrossFONLLunity3003->GetYaxis()->SetTitleSize(0.035);
   Graph_gaeRatioCrossFONLLunity3003->GetYaxis()->SetTitleOffset(0);
   Graph_gaeRatioCrossFONLLunity3003->GetYaxis()->SetTitleFont(42);
   Graph_gaeRatioCrossFONLLunity3003->GetZaxis()->SetLabelFont(42);
   Graph_gaeRatioCrossFONLLunity3003->GetZaxis()->SetLabelSize(0.035);
   Graph_gaeRatioCrossFONLLunity3003->GetZaxis()->SetTitleSize(0.035);
   Graph_gaeRatioCrossFONLLunity3003->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_gaeRatioCrossFONLLunity3003);
   
   grae->Draw("5");
   
   Double_t gaeRatioCrossFONLLstat_fx3004[7] = {
   6,
   8.5,
   12.5,
   17.5,
   25,
   40,
   75};
   Double_t gaeRatioCrossFONLLstat_fy3004[7] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t gaeRatioCrossFONLLstat_felx3004[7] = {
   1,
   1.5,
   2.5,
   2.5,
   5,
   10,
   25};
   Double_t gaeRatioCrossFONLLstat_fely3004[7] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t gaeRatioCrossFONLLstat_fehx3004[7] = {
   1,
   1.5,
   2.5,
   2.5,
   5,
   10,
   25};
   Double_t gaeRatioCrossFONLLstat_fehy3004[7] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   grae = new TGraphAsymmErrors(7,gaeRatioCrossFONLLstat_fx3004,gaeRatioCrossFONLLstat_fy3004,gaeRatioCrossFONLLstat_felx3004,gaeRatioCrossFONLLstat_fehx3004,gaeRatioCrossFONLLstat_fely3004,gaeRatioCrossFONLLstat_fehy3004);
   grae->SetName("gaeRatioCrossFONLLstat");
   grae->SetTitle("Graph");
   grae->SetFillStyle(1000);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(0.8);
   
   TH1F *Graph_gaeRatioCrossFONLLstat3004 = new TH1F("Graph_gaeRatioCrossFONLLstat3004","Graph",100,0,109.5);
   Graph_gaeRatioCrossFONLLstat3004->SetMinimum(0);
   Graph_gaeRatioCrossFONLLstat3004->SetMaximum(1.1);
   Graph_gaeRatioCrossFONLLstat3004->SetDirectory(0);
   Graph_gaeRatioCrossFONLLstat3004->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_gaeRatioCrossFONLLstat3004->SetLineColor(ci);
   Graph_gaeRatioCrossFONLLstat3004->SetMarkerStyle(20);
   Graph_gaeRatioCrossFONLLstat3004->GetXaxis()->SetLabelFont(42);
   Graph_gaeRatioCrossFONLLstat3004->GetXaxis()->SetLabelSize(0.035);
   Graph_gaeRatioCrossFONLLstat3004->GetXaxis()->SetTitleSize(0.035);
   Graph_gaeRatioCrossFONLLstat3004->GetXaxis()->SetTitleFont(42);
   Graph_gaeRatioCrossFONLLstat3004->GetYaxis()->SetLabelFont(42);
   Graph_gaeRatioCrossFONLLstat3004->GetYaxis()->SetLabelSize(0.035);
   Graph_gaeRatioCrossFONLLstat3004->GetYaxis()->SetTitleSize(0.035);
   Graph_gaeRatioCrossFONLLstat3004->GetYaxis()->SetTitleOffset(0);
   Graph_gaeRatioCrossFONLLstat3004->GetYaxis()->SetTitleFont(42);
   Graph_gaeRatioCrossFONLLstat3004->GetZaxis()->SetLabelFont(42);
   Graph_gaeRatioCrossFONLLstat3004->GetZaxis()->SetLabelSize(0.035);
   Graph_gaeRatioCrossFONLLstat3004->GetZaxis()->SetTitleSize(0.035);
   Graph_gaeRatioCrossFONLLstat3004->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_gaeRatioCrossFONLLstat3004);
   
   grae->Draw("ep");
   
   Double_t gaeRatioCrossFONLLsyst_fx3005[7] = {
   6,
   8.5,
   12.5,
   17.5,
   25,
   40,
   75};
   Double_t gaeRatioCrossFONLLsyst_fy3005[7] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t gaeRatioCrossFONLLsyst_felx3005[7] = {
   1,
   1.5,
   2.5,
   2.5,
   5,
   10,
   25};
   Double_t gaeRatioCrossFONLLsyst_fely3005[7] = {
   -nan,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t gaeRatioCrossFONLLsyst_fehx3005[7] = {
   1,
   1.5,
   2.5,
   2.5,
   5,
   10,
   25};
   Double_t gaeRatioCrossFONLLsyst_fehy3005[7] = {
   -nan,
   0,
   0,
   0,
   0,
   0,
   0};
   grae = new TGraphAsymmErrors(7,gaeRatioCrossFONLLsyst_fx3005,gaeRatioCrossFONLLsyst_fy3005,gaeRatioCrossFONLLsyst_felx3005,gaeRatioCrossFONLLsyst_fehx3005,gaeRatioCrossFONLLsyst_fely3005,gaeRatioCrossFONLLsyst_fehy3005);
   grae->SetName("gaeRatioCrossFONLLsyst");
   grae->SetTitle("Graph");
   grae->SetFillColor(5);
   grae->SetFillStyle(0);
   grae->SetLineWidth(2);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_gaeRatioCrossFONLLsyst3005 = new TH1F("Graph_gaeRatioCrossFONLLsyst3005","Graph",100,0,109.5);
   Graph_gaeRatioCrossFONLLsyst3005->SetMinimum(0);
   Graph_gaeRatioCrossFONLLsyst3005->SetMaximum(1.1);
   Graph_gaeRatioCrossFONLLsyst3005->SetDirectory(0);
   Graph_gaeRatioCrossFONLLsyst3005->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_gaeRatioCrossFONLLsyst3005->SetLineColor(ci);
   Graph_gaeRatioCrossFONLLsyst3005->SetMarkerStyle(20);
   Graph_gaeRatioCrossFONLLsyst3005->GetXaxis()->SetLabelFont(42);
   Graph_gaeRatioCrossFONLLsyst3005->GetXaxis()->SetLabelSize(0.035);
   Graph_gaeRatioCrossFONLLsyst3005->GetXaxis()->SetTitleSize(0.035);
   Graph_gaeRatioCrossFONLLsyst3005->GetXaxis()->SetTitleFont(42);
   Graph_gaeRatioCrossFONLLsyst3005->GetYaxis()->SetLabelFont(42);
   Graph_gaeRatioCrossFONLLsyst3005->GetYaxis()->SetLabelSize(0.035);
   Graph_gaeRatioCrossFONLLsyst3005->GetYaxis()->SetTitleSize(0.035);
   Graph_gaeRatioCrossFONLLsyst3005->GetYaxis()->SetTitleOffset(0);
   Graph_gaeRatioCrossFONLLsyst3005->GetYaxis()->SetTitleFont(42);
   Graph_gaeRatioCrossFONLLsyst3005->GetZaxis()->SetLabelFont(42);
   Graph_gaeRatioCrossFONLLsyst3005->GetZaxis()->SetLabelSize(0.035);
   Graph_gaeRatioCrossFONLLsyst3005->GetZaxis()->SetTitleSize(0.035);
   Graph_gaeRatioCrossFONLLsyst3005->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_gaeRatioCrossFONLLsyst3005);
   
   grae->Draw("5");
   pRatio->Modified();
   cSigma->cd();
   cSigma->Modified();
   cSigma->cd();
   cSigma->SetSelected(cSigma);
}
