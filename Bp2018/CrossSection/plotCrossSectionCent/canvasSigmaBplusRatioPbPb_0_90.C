void canvasSigmaBplusRatioPbPb_0_90()
{
//=========Macro generated from canvas: cSigma/
//=========  (Tue Oct  8 20:45:47 2019) by ROOT version 6.16/00
   TCanvas *cSigma = new TCanvas("cSigma", "",0,0,600,600);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   cSigma->SetHighLightColor(2);
   cSigma->Range(-60,4.788404,340,7.13947);
   cSigma->SetFillColor(0);
   cSigma->SetBorderMode(0);
   cSigma->SetBorderSize(2);
   cSigma->SetLogy();
   cSigma->SetLeftMargin(0.15);
   cSigma->SetTopMargin(0.05932203);
   cSigma->SetFrameBorderMode(0);
   cSigma->SetFrameBorderMode(0);
   
   TH2F *hemptySigma__1 = new TH2F("hemptySigma__1","",50,0,300,10,100000,1e+07);
   hemptySigma__1->SetMinimum(0);
   hemptySigma__1->SetMaximum(2);
   hemptySigma__1->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   hemptySigma__1->SetLineColor(ci);
   hemptySigma__1->SetMarkerStyle(20);
   hemptySigma__1->GetXaxis()->SetTitle("<N_{part}>");
   hemptySigma__1->GetXaxis()->CenterTitle(true);
   hemptySigma__1->GetXaxis()->SetLabelFont(42);
   hemptySigma__1->GetXaxis()->SetLabelOffset(0.0015);
   hemptySigma__1->GetXaxis()->SetLabelSize(0.036);
   hemptySigma__1->GetXaxis()->SetTitleSize(0.036);
   hemptySigma__1->GetXaxis()->SetTitleOffset(1);
   hemptySigma__1->GetXaxis()->SetTitleFont(42);
   hemptySigma__1->GetYaxis()->SetTitle("#frac{1}{T_{AA}} #frac{dN}{dp_{T}} (mb*GeV^{-1}c)");
   hemptySigma__1->GetYaxis()->CenterTitle(true);
   hemptySigma__1->GetYaxis()->SetLabelFont(42);
   hemptySigma__1->GetYaxis()->SetLabelSize(0.042);
   hemptySigma__1->GetYaxis()->SetTitleSize(0.042);
   hemptySigma__1->GetYaxis()->SetTitleOffset(1.428571);
   hemptySigma__1->GetYaxis()->SetTitleFont(42);
   hemptySigma__1->GetZaxis()->SetLabelFont(42);
   hemptySigma__1->GetZaxis()->SetLabelSize(0.035);
   hemptySigma__1->GetZaxis()->SetTitleSize(0.035);
   hemptySigma__1->GetZaxis()->SetTitleOffset(1);
   hemptySigma__1->GetZaxis()->SetTitleFont(42);
   hemptySigma__1->Draw("");
   
   Double_t hPtSigma_fx3001[1] = {
   126};
   Double_t hPtSigma_fy3001[1] = {
   30930.94};
   Double_t hPtSigma_felx3001[1] = {
   2.5};
   Double_t hPtSigma_fely3001[1] = {
   1592.587};
   Double_t hPtSigma_fehx3001[1] = {
   2.5};
   Double_t hPtSigma_fehy3001[1] = {
   1592.587};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(1,hPtSigma_fx3001,hPtSigma_fy3001,hPtSigma_felx3001,hPtSigma_fehx3001,hPtSigma_fely3001,hPtSigma_fehy3001);
   grae->SetName("hPtSigma");
   grae->SetTitle("Graph");
   grae->SetFillStyle(1000);
   grae->SetLineWidth(2);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(0.84);
   
   TH1F *Graph_hPtSigma3001 = new TH1F("Graph_hPtSigma3001","Graph",100,123,129);
   Graph_hPtSigma3001->SetMinimum(29019.84);
   Graph_hPtSigma3001->SetMaximum(32842.04);
   Graph_hPtSigma3001->SetDirectory(0);
   Graph_hPtSigma3001->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_hPtSigma3001->SetLineColor(ci);
   Graph_hPtSigma3001->SetMarkerStyle(20);
   Graph_hPtSigma3001->GetXaxis()->SetLabelFont(42);
   Graph_hPtSigma3001->GetXaxis()->SetLabelSize(0.035);
   Graph_hPtSigma3001->GetXaxis()->SetTitleSize(0.035);
   Graph_hPtSigma3001->GetXaxis()->SetTitleOffset(1);
   Graph_hPtSigma3001->GetXaxis()->SetTitleFont(42);
   Graph_hPtSigma3001->GetYaxis()->SetLabelFont(42);
   Graph_hPtSigma3001->GetYaxis()->SetLabelSize(0.035);
   Graph_hPtSigma3001->GetYaxis()->SetTitleSize(0.035);
   Graph_hPtSigma3001->GetYaxis()->SetTitleFont(42);
   Graph_hPtSigma3001->GetZaxis()->SetLabelFont(42);
   Graph_hPtSigma3001->GetZaxis()->SetLabelSize(0.035);
   Graph_hPtSigma3001->GetZaxis()->SetTitleSize(0.035);
   Graph_hPtSigma3001->GetZaxis()->SetTitleOffset(1);
   Graph_hPtSigma3001->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_hPtSigma3001);
   
   grae->Draw("ep");
   
   Double_t gaeCrossSyst_fx3002[1] = {
   126};
   Double_t gaeCrossSyst_fy3002[1] = {
   30930.94};
   Double_t gaeCrossSyst_felx3002[1] = {
   2.5};
   Double_t gaeCrossSyst_fely3002[1] = {
   3380.286};
   Double_t gaeCrossSyst_fehx3002[1] = {
   2.5};
   Double_t gaeCrossSyst_fehy3002[1] = {
   3380.286};
   grae = new TGraphAsymmErrors(1,gaeCrossSyst_fx3002,gaeCrossSyst_fy3002,gaeCrossSyst_felx3002,gaeCrossSyst_fehx3002,gaeCrossSyst_fely3002,gaeCrossSyst_fehy3002);
   grae->SetName("gaeCrossSyst");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);
   grae->SetFillStyle(0);
   grae->SetLineWidth(2);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(1.142857);
   
   TH1F *Graph_gaeCrossSyst3002 = new TH1F("Graph_gaeCrossSyst3002","Graph",100,123,129);
   Graph_gaeCrossSyst3002->SetMinimum(26874.6);
   Graph_gaeCrossSyst3002->SetMaximum(34987.28);
   Graph_gaeCrossSyst3002->SetDirectory(0);
   Graph_gaeCrossSyst3002->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_gaeCrossSyst3002->SetLineColor(ci);
   Graph_gaeCrossSyst3002->SetMarkerStyle(20);
   Graph_gaeCrossSyst3002->GetXaxis()->SetLabelFont(42);
   Graph_gaeCrossSyst3002->GetXaxis()->SetLabelSize(0.035);
   Graph_gaeCrossSyst3002->GetXaxis()->SetTitleSize(0.035);
   Graph_gaeCrossSyst3002->GetXaxis()->SetTitleOffset(1);
   Graph_gaeCrossSyst3002->GetXaxis()->SetTitleFont(42);
   Graph_gaeCrossSyst3002->GetYaxis()->SetLabelFont(42);
   Graph_gaeCrossSyst3002->GetYaxis()->SetLabelSize(0.035);
   Graph_gaeCrossSyst3002->GetYaxis()->SetTitleSize(0.035);
   Graph_gaeCrossSyst3002->GetYaxis()->SetTitleFont(42);
   Graph_gaeCrossSyst3002->GetZaxis()->SetLabelFont(42);
   Graph_gaeCrossSyst3002->GetZaxis()->SetLabelSize(0.035);
   Graph_gaeCrossSyst3002->GetZaxis()->SetTitleSize(0.035);
   Graph_gaeCrossSyst3002->GetZaxis()->SetTitleOffset(1);
   Graph_gaeCrossSyst3002->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_gaeCrossSyst3002);
   
   grae->Draw("5");
   TLatex *   tex = new TLatex(0.52,0.916,"CMS");
tex->SetNDC();
   tex->SetTextAlign(13);
   tex->SetTextSize(0.056);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.9,0.9552,"1.5 nb^{-1} (PbPb 5.02 TeV)");
tex->SetNDC();
   tex->SetTextAlign(31);
   tex->SetTextFont(42);
   tex->SetTextSize(0.035);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.53,0.755,"|y| < 2.4");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.035);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.7,0.874,"B^{+}");
tex->SetNDC();
   tex->SetTextSize(0.063);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.53,0.713,"Global uncert. #plus3.6, #minus3.6%");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.035);
   tex->SetLineWidth(2);
   tex->Draw();
   
   TLegend *leg = new TLegend(0.52,0.79,0.85,0.86,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.035);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("hPtSigma","Data","pf");
   entry->SetFillStyle(1000);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(0.84);
   entry->SetTextFont(42);
   leg->Draw();
   cSigma->Modified();
   cSigma->cd();
   cSigma->SetSelected(cSigma);
}
