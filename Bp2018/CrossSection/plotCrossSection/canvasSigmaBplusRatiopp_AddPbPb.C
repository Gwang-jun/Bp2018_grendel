void canvasSigmaBplusRatiopp_AddPbPb()
{
//=========Macro generated from canvas: cSigma/
//=========  (Sun Mar 10 22:37:53 2019) by ROOT version 6.14/04
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
   pSigma->Range(-19.04282,3,112.1411,7.37037);
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
   
   TH2F *hemptySigma__1 = new TH2F("hemptySigma__1","",50,0,105,10,1000,1e+07);
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
   
   Double_t Graph0_fx3001[7] = {
   6,
   8.5,
   12.5,
   17.5,
   25,
   40,
   75};
   Double_t Graph0_fy3001[7] = {
   4275200,
   1935200,
   585640,
   165770,
   36029,
   3715.8,
   120.19};
   Double_t Graph0_felx3001[7] = {
   1,
   1.5,
   2.5,
   2.5,
   5,
   10,
   25};
   Double_t Graph0_fely3001[7] = {
   1236800,
   501900,
   138040,
   34820,
   6525,
   545.5,
   13.78};
   Double_t Graph0_fehx3001[7] = {
   1,
   1.5,
   2.5,
   2.5,
   5,
   10,
   25};
   Double_t Graph0_fehy3001[7] = {
   1848200,
   786900,
   209010,
   50110,
   8761,
   662,
   14.86};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(7,Graph0_fx3001,Graph0_fy3001,Graph0_felx3001,Graph0_fehx3001,Graph0_fely3001,Graph0_fehy3001);
   grae->SetName("Graph0");
   grae->SetTitle("FONLL result with FF=0.405");

   ci = TColor::GetColor("#ffcc00");
   grae->SetFillColor(ci);

   ci = TColor::GetColor("#ffcc00");
   grae->SetLineColor(ci);
   grae->SetLineWidth(3);
   
   TH1F *Graph_Graph03001 = new TH1F("Graph_Graph03001","FONLL result with FF=0.405",100,0,109.5);
   Graph_Graph03001->SetMinimum(95.769);
   Graph_Graph03001->SetMaximum(6735729);
   Graph_Graph03001->SetDirectory(0);
   Graph_Graph03001->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph03001->SetLineColor(ci);
   Graph_Graph03001->SetMarkerStyle(20);
   Graph_Graph03001->GetXaxis()->SetLabelFont(42);
   Graph_Graph03001->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph03001->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph03001->GetXaxis()->SetTitleFont(42);
   Graph_Graph03001->GetYaxis()->SetLabelFont(42);
   Graph_Graph03001->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph03001->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph03001->GetYaxis()->SetTitleOffset(0);
   Graph_Graph03001->GetYaxis()->SetTitleFont(42);
   Graph_Graph03001->GetZaxis()->SetLabelFont(42);
   Graph_Graph03001->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph03001->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph03001->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph03001);
   
   grae->Draw("5");
   Double_t xAxis1[7] = {7, 10, 15, 20, 30, 50, 100}; 
   
   TH1D *hPtSigma__2 = new TH1D("hPtSigma__2","",6, xAxis1);
   hPtSigma__2->SetBinContent(1,2639149);
   hPtSigma__2->SetBinContent(2,753623.6);
   hPtSigma__2->SetBinContent(3,199443.7);
   hPtSigma__2->SetBinContent(4,47072.32);
   hPtSigma__2->SetBinContent(5,5409.234);
   hPtSigma__2->SetBinContent(6,307.8495);
   hPtSigma__2->SetBinError(1,169170.8);
   hPtSigma__2->SetBinError(2,29443.71);
   hPtSigma__2->SetBinError(3,9142.84);
   hPtSigma__2->SetBinError(4,2469.697);
   hPtSigma__2->SetBinError(5,501.2322);
   hPtSigma__2->SetBinError(6,70.55491);
   hPtSigma__2->SetEntries(449.2236);
   hPtSigma__2->SetLineWidth(2);
   hPtSigma__2->SetMarkerStyle(20);
   hPtSigma__2->SetMarkerSize(1.2);
   hPtSigma__2->GetXaxis()->SetTitle("D^{0} p_{T} (GeV/c)");
   hPtSigma__2->GetXaxis()->SetLabelFont(42);
   hPtSigma__2->GetXaxis()->SetLabelSize(0.035);
   hPtSigma__2->GetXaxis()->SetTitleSize(0.035);
   hPtSigma__2->GetXaxis()->SetTitleFont(42);
   hPtSigma__2->GetYaxis()->SetTitle("Uncorrected dN(D^{0})/dp_{T}");
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
   2639149,
   753623.6,
   199443.7,
   47072.32,
   5409.234,
   307.8495,
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
   67773.1,
   15912.01,
   3634.567,
   411.1435,
   23.15373,
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
   67773.1,
   15912.01,
   3634.567,
   411.1435,
   23.15373,
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
   Graph_gaeCrossSyst3002->SetMinimum(2903.064);
   Graph_gaeCrossSyst3002->SetMaximum(2903064);
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
   entry=leg->AddEntry("Graph0","FONLL","f");

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
   
   Double_t gaeCrossSyst_fx3003[7] = {
   6,
   8.5,
   12.5,
   17.5,
   25,
   40,
   75};
   Double_t gaeCrossSyst_fy3003[7] = {
   -184430.3,
   652606.9,
   175290.8,
   64935.93,
   18255.62,
   2999.816,
   150.4674};
   Double_t gaeCrossSyst_felx3003[7] = {
   1,
   1.5,
   2.5,
   2.5,
   5,
   10,
   25};
   Double_t gaeCrossSyst_fely3003[7] = {
   -107505.6,
   133167.2,
   28030.16,
   9760.031,
   2689.442,
   443.5832,
   11.74071};
   Double_t gaeCrossSyst_fehx3003[7] = {
   1,
   1.5,
   2.5,
   2.5,
   5,
   10,
   25};
   Double_t gaeCrossSyst_fehy3003[7] = {
   -107505.6,
   133167.2,
   28030.16,
   9760.031,
   2689.442,
   443.5832,
   11.74071};
   grae = new TGraphAsymmErrors(7,gaeCrossSyst_fx3003,gaeCrossSyst_fy3003,gaeCrossSyst_felx3003,gaeCrossSyst_fehx3003,gaeCrossSyst_fely3003,gaeCrossSyst_fehy3003);
   grae->SetName("gaeCrossSyst");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);
   grae->SetFillStyle(0);
   grae->SetLineColor(2);
   grae->SetLineWidth(2);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(0.8);
   
   TH1F *Graph_Graph_gaeCrossSyst30013003 = new TH1F("Graph_Graph_gaeCrossSyst30013003","Graph",100,0,109.5);
   Graph_Graph_gaeCrossSyst30013003->SetMinimum(882.7946);
   Graph_Graph_gaeCrossSyst30013003->SetMaximum(882794.6);
   Graph_Graph_gaeCrossSyst30013003->SetDirectory(0);
   Graph_Graph_gaeCrossSyst30013003->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_gaeCrossSyst30013003->SetLineColor(ci);
   Graph_Graph_gaeCrossSyst30013003->SetMarkerStyle(20);
   Graph_Graph_gaeCrossSyst30013003->GetXaxis()->SetLabelFont(42);
   Graph_Graph_gaeCrossSyst30013003->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_gaeCrossSyst30013003->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_gaeCrossSyst30013003->GetXaxis()->SetTitleFont(42);
   Graph_Graph_gaeCrossSyst30013003->GetYaxis()->SetLabelFont(42);
   Graph_Graph_gaeCrossSyst30013003->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_gaeCrossSyst30013003->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_gaeCrossSyst30013003->GetYaxis()->SetTitleOffset(0);
   Graph_Graph_gaeCrossSyst30013003->GetYaxis()->SetTitleFont(42);
   Graph_Graph_gaeCrossSyst30013003->GetZaxis()->SetLabelFont(42);
   Graph_Graph_gaeCrossSyst30013003->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_gaeCrossSyst30013003->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_gaeCrossSyst30013003->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_gaeCrossSyst30013003);
   
   grae->Draw("5");
   Double_t xAxis2[8] = {5, 7, 10, 15, 20, 30, 50, 100}; 
   
   TH1D *hPtSigma__3 = new TH1D("hPtSigma__3","",7, xAxis2);
   hPtSigma__3->SetBinContent(1,-184430.3);
   hPtSigma__3->SetBinContent(2,652606.9);
   hPtSigma__3->SetBinContent(3,175290.8);
   hPtSigma__3->SetBinContent(4,64935.93);
   hPtSigma__3->SetBinContent(5,18255.62);
   hPtSigma__3->SetBinContent(6,2999.816);
   hPtSigma__3->SetBinContent(7,150.4674);
   hPtSigma__3->SetBinError(1,1554527);
   hPtSigma__3->SetBinError(2,143900.8);
   hPtSigma__3->SetBinError(3,18864.27);
   hPtSigma__3->SetBinError(4,5378.366);
   hPtSigma__3->SetBinError(5,1436.001);
   hPtSigma__3->SetBinError(6,322.814);
   hPtSigma__3->SetBinError(7,23.08711);
   hPtSigma__3->SetEntries(0.2184981);
   hPtSigma__3->SetLineColor(2);
   hPtSigma__3->SetLineWidth(2);
   hPtSigma__3->SetMarkerColor(2);
   hPtSigma__3->SetMarkerStyle(21);
   hPtSigma__3->SetMarkerSize(1.2);
   hPtSigma__3->GetXaxis()->SetTitle("D^{0} p_{T} (GeV/c)");
   hPtSigma__3->GetXaxis()->SetLabelFont(42);
   hPtSigma__3->GetXaxis()->SetLabelSize(0.035);
   hPtSigma__3->GetXaxis()->SetTitleSize(0.035);
   hPtSigma__3->GetXaxis()->SetTitleFont(42);
   hPtSigma__3->GetYaxis()->SetTitle("Uncorrected dN(D^{0})/dp_{T}");
   hPtSigma__3->GetYaxis()->SetLabelFont(42);
   hPtSigma__3->GetYaxis()->SetLabelSize(0.035);
   hPtSigma__3->GetYaxis()->SetTitleSize(0.035);
   hPtSigma__3->GetYaxis()->SetTitleOffset(0);
   hPtSigma__3->GetYaxis()->SetTitleFont(42);
   hPtSigma__3->GetZaxis()->SetLabelFont(42);
   hPtSigma__3->GetZaxis()->SetLabelSize(0.035);
   hPtSigma__3->GetZaxis()->SetTitleSize(0.035);
   hPtSigma__3->GetZaxis()->SetTitleFont(42);
   hPtSigma__3->Draw("epsame");
   
   leg = new TLegend(0.52,0.7,0.85,0.85,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   entry=leg->AddEntry("hPtSigma","Data pp","pf");
   entry->SetFillStyle(1001);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("hPtSigma","Data PbPb","pf");
   entry->SetFillStyle(1001);
   entry->SetLineColor(2);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(2);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph0","FONLL pp ref.","f");

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
      tex = new TLatex(0.53,0.594,"Global uncert.");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.53,0.55,"pp: 3.8%");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.53,0.5,"PbPb: #plus4.6, #minus5.0%");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->Draw();
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
   
   TH2F *hemptyRatio__4 = new TH2F("hemptyRatio__4","",50,0,105,10,0.2,1.8);
   hemptyRatio__4->SetStats(0);

   ci = TColor::GetColor("#000099");
   hemptyRatio__4->SetLineColor(ci);
   hemptyRatio__4->SetMarkerStyle(20);
   hemptyRatio__4->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   hemptyRatio__4->GetXaxis()->CenterTitle(true);
   hemptyRatio__4->GetXaxis()->SetLabelFont(42);
   hemptyRatio__4->GetXaxis()->SetLabelSize(0.12);
   hemptyRatio__4->GetXaxis()->SetTitleSize(0.12);
   hemptyRatio__4->GetXaxis()->SetTickLength(0.09999999);
   hemptyRatio__4->GetXaxis()->SetTitleFont(42);
   hemptyRatio__4->GetYaxis()->SetTitle("pp / FONLL");
   hemptyRatio__4->GetYaxis()->CenterTitle(true);
   hemptyRatio__4->GetYaxis()->SetNdivisions(505);
   hemptyRatio__4->GetYaxis()->SetLabelFont(42);
   hemptyRatio__4->GetYaxis()->SetLabelSize(0.12);
   hemptyRatio__4->GetYaxis()->SetTitleSize(0.12);
   hemptyRatio__4->GetYaxis()->SetTitleOffset(0.5);
   hemptyRatio__4->GetYaxis()->SetTitleFont(42);
   hemptyRatio__4->GetZaxis()->SetLabelFont(42);
   hemptyRatio__4->GetZaxis()->SetLabelSize(0.035);
   hemptyRatio__4->GetZaxis()->SetTitleSize(0.035);
   hemptyRatio__4->GetZaxis()->SetTitleFont(42);
   hemptyRatio__4->Draw("");
   TLine *line = new TLine(0,1,105,1);
   line->SetLineStyle(2);
   line->Draw();
   
   Double_t gaeRatioCrossFONLLunity_fx3004[7] = {
   6,
   8.5,
   12.5,
   17.5,
   25,
   40,
   75};
   Double_t gaeRatioCrossFONLLunity_fy3004[7] = {
   1,
   1,
   1,
   1,
   1,
   1,
   1};
   Double_t gaeRatioCrossFONLLunity_felx3004[7] = {
   1,
   1.5,
   2.5,
   2.5,
   5,
   10,
   25};
   Double_t gaeRatioCrossFONLLunity_fely3004[7] = {
   0.2892964,
   0.259353,
   0.2357079,
   0.2100501,
   0.1811041,
   0.1468055,
   0.1146518};
   Double_t gaeRatioCrossFONLLunity_fehx3004[7] = {
   1,
   1.5,
   2.5,
   2.5,
   5,
   10,
   25};
   Double_t gaeRatioCrossFONLLunity_fehy3004[7] = {
   0.4323073,
   0.4066246,
   0.3568916,
   0.3022863,
   0.2431652,
   0.1781581,
   0.1236376};
   grae = new TGraphAsymmErrors(7,gaeRatioCrossFONLLunity_fx3004,gaeRatioCrossFONLLunity_fy3004,gaeRatioCrossFONLLunity_felx3004,gaeRatioCrossFONLLunity_fehx3004,gaeRatioCrossFONLLunity_fely3004,gaeRatioCrossFONLLunity_fehy3004);
   grae->SetName("gaeRatioCrossFONLLunity");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#ffcc00");
   grae->SetFillColor(ci);

   ci = TColor::GetColor("#ffcc00");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_gaeRatioCrossFONLLunity3004 = new TH1F("Graph_gaeRatioCrossFONLLunity3004","Graph",100,0,109.5);
   Graph_gaeRatioCrossFONLLunity3004->SetMinimum(0.6385432);
   Graph_gaeRatioCrossFONLLunity3004->SetMaximum(1.504468);
   Graph_gaeRatioCrossFONLLunity3004->SetDirectory(0);
   Graph_gaeRatioCrossFONLLunity3004->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_gaeRatioCrossFONLLunity3004->SetLineColor(ci);
   Graph_gaeRatioCrossFONLLunity3004->SetMarkerStyle(20);
   Graph_gaeRatioCrossFONLLunity3004->GetXaxis()->SetLabelFont(42);
   Graph_gaeRatioCrossFONLLunity3004->GetXaxis()->SetLabelSize(0.035);
   Graph_gaeRatioCrossFONLLunity3004->GetXaxis()->SetTitleSize(0.035);
   Graph_gaeRatioCrossFONLLunity3004->GetXaxis()->SetTitleFont(42);
   Graph_gaeRatioCrossFONLLunity3004->GetYaxis()->SetLabelFont(42);
   Graph_gaeRatioCrossFONLLunity3004->GetYaxis()->SetLabelSize(0.035);
   Graph_gaeRatioCrossFONLLunity3004->GetYaxis()->SetTitleSize(0.035);
   Graph_gaeRatioCrossFONLLunity3004->GetYaxis()->SetTitleOffset(0);
   Graph_gaeRatioCrossFONLLunity3004->GetYaxis()->SetTitleFont(42);
   Graph_gaeRatioCrossFONLLunity3004->GetZaxis()->SetLabelFont(42);
   Graph_gaeRatioCrossFONLLunity3004->GetZaxis()->SetLabelSize(0.035);
   Graph_gaeRatioCrossFONLLunity3004->GetZaxis()->SetTitleSize(0.035);
   Graph_gaeRatioCrossFONLLunity3004->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_gaeRatioCrossFONLLunity3004);
   
   grae->Draw("5");
   
   Double_t gaeRatioCrossFONLLstat_fx3005[7] = {
   6,
   8.5,
   12.5,
   17.5,
   25,
   40,
   75};
   Double_t gaeRatioCrossFONLLstat_fy3005[7] = {
   0.617316,
   0.3894293,
   0.3405568,
   0.2839617,
   0.1501356,
   0.08284879,
   0};
   Double_t gaeRatioCrossFONLLstat_felx3005[7] = {
   1,
   1.5,
   2.5,
   2.5,
   5,
   10,
   25};
   Double_t gaeRatioCrossFONLLstat_fely3005[7] = {
   0.03957026,
   0.01521482,
   0.01561171,
   0.01489834,
   0.01391191,
   0.01898781,
   0};
   Double_t gaeRatioCrossFONLLstat_fehx3005[7] = {
   1,
   1.5,
   2.5,
   2.5,
   5,
   10,
   25};
   Double_t gaeRatioCrossFONLLstat_fehy3005[7] = {
   0.03957026,
   0.01521482,
   0.01561171,
   0.01489834,
   0.01391191,
   0.01898781,
   0};
   grae = new TGraphAsymmErrors(7,gaeRatioCrossFONLLstat_fx3005,gaeRatioCrossFONLLstat_fy3005,gaeRatioCrossFONLLstat_felx3005,gaeRatioCrossFONLLstat_fehx3005,gaeRatioCrossFONLLstat_fely3005,gaeRatioCrossFONLLstat_fehy3005);
   grae->SetName("gaeRatioCrossFONLLstat");
   grae->SetTitle("Graph");
   grae->SetFillStyle(1000);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(0.8);
   
   TH1F *Graph_gaeRatioCrossFONLLstat3005 = new TH1F("Graph_gaeRatioCrossFONLLstat3005","Graph",100,0,109.5);
   Graph_gaeRatioCrossFONLLstat3005->SetMinimum(0);
   Graph_gaeRatioCrossFONLLstat3005->SetMaximum(0.7225749);
   Graph_gaeRatioCrossFONLLstat3005->SetDirectory(0);
   Graph_gaeRatioCrossFONLLstat3005->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_gaeRatioCrossFONLLstat3005->SetLineColor(ci);
   Graph_gaeRatioCrossFONLLstat3005->SetMarkerStyle(20);
   Graph_gaeRatioCrossFONLLstat3005->GetXaxis()->SetLabelFont(42);
   Graph_gaeRatioCrossFONLLstat3005->GetXaxis()->SetLabelSize(0.035);
   Graph_gaeRatioCrossFONLLstat3005->GetXaxis()->SetTitleSize(0.035);
   Graph_gaeRatioCrossFONLLstat3005->GetXaxis()->SetTitleFont(42);
   Graph_gaeRatioCrossFONLLstat3005->GetYaxis()->SetLabelFont(42);
   Graph_gaeRatioCrossFONLLstat3005->GetYaxis()->SetLabelSize(0.035);
   Graph_gaeRatioCrossFONLLstat3005->GetYaxis()->SetTitleSize(0.035);
   Graph_gaeRatioCrossFONLLstat3005->GetYaxis()->SetTitleOffset(0);
   Graph_gaeRatioCrossFONLLstat3005->GetYaxis()->SetTitleFont(42);
   Graph_gaeRatioCrossFONLLstat3005->GetZaxis()->SetLabelFont(42);
   Graph_gaeRatioCrossFONLLstat3005->GetZaxis()->SetLabelSize(0.035);
   Graph_gaeRatioCrossFONLLstat3005->GetZaxis()->SetTitleSize(0.035);
   Graph_gaeRatioCrossFONLLstat3005->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_gaeRatioCrossFONLLstat3005);
   
   grae->Draw("ep");
   
   Double_t gaeRatioCrossFONLLsyst_fx3006[7] = {
   6,
   8.5,
   12.5,
   17.5,
   25,
   40,
   75};
   Double_t gaeRatioCrossFONLLsyst_fy3006[7] = {
   0.617316,
   0.3894293,
   0.3405568,
   0.2839617,
   0.1501356,
   0.08284879,
   0};
   Double_t gaeRatioCrossFONLLsyst_felx3006[7] = {
   1,
   1.5,
   2.5,
   2.5,
   5,
   10,
   25};
   Double_t gaeRatioCrossFONLLsyst_fely3006[7] = {
   -nan,
   0.03502124,
   0.02717029,
   0.02192536,
   0.01141146,
   0.006231156,
   0};
   Double_t gaeRatioCrossFONLLsyst_fehx3006[7] = {
   1,
   1.5,
   2.5,
   2.5,
   5,
   10,
   25};
   Double_t gaeRatioCrossFONLLsyst_fehy3006[7] = {
   -nan,
   0.03502124,
   0.02717029,
   0.02192536,
   0.01141146,
   0.006231156,
   0};
   grae = new TGraphAsymmErrors(7,gaeRatioCrossFONLLsyst_fx3006,gaeRatioCrossFONLLsyst_fy3006,gaeRatioCrossFONLLsyst_felx3006,gaeRatioCrossFONLLsyst_fehx3006,gaeRatioCrossFONLLsyst_fely3006,gaeRatioCrossFONLLsyst_fehy3006);
   grae->SetName("gaeRatioCrossFONLLsyst");
   grae->SetTitle("Graph");
   grae->SetFillColor(5);
   grae->SetFillStyle(0);
   grae->SetLineWidth(2);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_gaeRatioCrossFONLLsyst3006 = new TH1F("Graph_gaeRatioCrossFONLLsyst3006","Graph",100,0,109.5);
   Graph_gaeRatioCrossFONLLsyst3006->SetMinimum(0);
   Graph_gaeRatioCrossFONLLsyst3006->SetMaximum(0.6790476);
   Graph_gaeRatioCrossFONLLsyst3006->SetDirectory(0);
   Graph_gaeRatioCrossFONLLsyst3006->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_gaeRatioCrossFONLLsyst3006->SetLineColor(ci);
   Graph_gaeRatioCrossFONLLsyst3006->SetMarkerStyle(20);
   Graph_gaeRatioCrossFONLLsyst3006->GetXaxis()->SetLabelFont(42);
   Graph_gaeRatioCrossFONLLsyst3006->GetXaxis()->SetLabelSize(0.035);
   Graph_gaeRatioCrossFONLLsyst3006->GetXaxis()->SetTitleSize(0.035);
   Graph_gaeRatioCrossFONLLsyst3006->GetXaxis()->SetTitleFont(42);
   Graph_gaeRatioCrossFONLLsyst3006->GetYaxis()->SetLabelFont(42);
   Graph_gaeRatioCrossFONLLsyst3006->GetYaxis()->SetLabelSize(0.035);
   Graph_gaeRatioCrossFONLLsyst3006->GetYaxis()->SetTitleSize(0.035);
   Graph_gaeRatioCrossFONLLsyst3006->GetYaxis()->SetTitleOffset(0);
   Graph_gaeRatioCrossFONLLsyst3006->GetYaxis()->SetTitleFont(42);
   Graph_gaeRatioCrossFONLLsyst3006->GetZaxis()->SetLabelFont(42);
   Graph_gaeRatioCrossFONLLsyst3006->GetZaxis()->SetLabelSize(0.035);
   Graph_gaeRatioCrossFONLLsyst3006->GetZaxis()->SetTitleSize(0.035);
   Graph_gaeRatioCrossFONLLsyst3006->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_gaeRatioCrossFONLLsyst3006);
   
   grae->Draw("5");
   pRatio->Modified();
   cSigma->cd();
   cSigma->Modified();
   cSigma->cd();
   cSigma->SetSelected(cSigma);
}
