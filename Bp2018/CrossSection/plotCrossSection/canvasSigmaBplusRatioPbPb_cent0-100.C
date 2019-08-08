void canvasSigmaBplusRatioPbPb_cent0-100()
{
//=========Macro generated from canvas: cSigma/
//=========  (Fri May 10 19:24:24 2019) by ROOT version 6.12/07
   TCanvas *cSigma = new TCanvas("cSigma", "",0,0,600,600);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   cSigma->SetHighLightColor(2);
   cSigma->Range(-19.04282,1.365212,112.1411,8.41841);
   cSigma->SetFillColor(0);
   cSigma->SetBorderMode(0);
   cSigma->SetBorderSize(2);
   cSigma->SetLogy();
   cSigma->SetLeftMargin(0.1451613);
   cSigma->SetRightMargin(0.05443548);
   cSigma->SetTopMargin(0.05932203);
   cSigma->SetFrameBorderMode(0);
   cSigma->SetFrameBorderMode(0);
   
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
   hemptySigma__1->GetYaxis()->SetTitle("#frac{1}{T_{AA}} #frac{dN}{dp_{T}} ( pb GeV^{-1}c)");
   hemptySigma__1->GetYaxis()->CenterTitle(true);
   hemptySigma__1->GetYaxis()->SetLabelFont(42);
   hemptySigma__1->GetYaxis()->SetLabelSize(0.042);
   hemptySigma__1->GetYaxis()->SetTitleSize(0.042);
   hemptySigma__1->GetYaxis()->SetTitleOffset(1.428571);
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
   4405494,
   2052997,
   652850.7,
   179221.3,
   42701.85,
   5133.226,
   244.7777};
   Double_t gaeSigmaBplus_felx3001[7] = {
   1,
   1.5,
   2.5,
   2.5,
   5,
   10,
   25};
   Double_t gaeSigmaBplus_fely3001[7] = {
   1195828,
   502245.6,
   149751.2,
   37027.43,
   7813.787,
   786.1306,
   30.29715};
   Double_t gaeSigmaBplus_fehx3001[7] = {
   1,
   1.5,
   2.5,
   2.5,
   5,
   10,
   25};
   Double_t gaeSigmaBplus_fehy3001[7] = {
   1812953,
   814792.5,
   233301.5,
   54366.09,
   10752.1,
   984.792,
   34.06047};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(7,gaeSigmaBplus_fx3001,gaeSigmaBplus_fy3001,gaeSigmaBplus_felx3001,gaeSigmaBplus_fehx3001,gaeSigmaBplus_fely3001,gaeSigmaBplus_fehy3001);
   grae->SetName("gaeSigmaBplus");
   grae->SetTitle(";p_{T}(GeV/c);d#sigma/dp_{T} (B^{+}) (pb GeV-1c)");

   ci = TColor::GetColor("#ffcc00");
   grae->SetFillColor(ci);

   ci = TColor::GetColor("#ffcc00");
   grae->SetLineColor(ci);
   grae->SetLineWidth(3);
   
   TH1F *Graph_gaeSigmaBplus3001 = new TH1F("Graph_gaeSigmaBplus3001","",100,0,109.5);
   Graph_gaeSigmaBplus3001->SetMinimum(193.0325);
   Graph_gaeSigmaBplus3001->SetMaximum(6840269);
   Graph_gaeSigmaBplus3001->SetDirectory(0);
   Graph_gaeSigmaBplus3001->SetStats(0);
   Graph_gaeSigmaBplus3001->GetXaxis()->SetTitle("p_{T}(GeV/c)");
   Graph_gaeSigmaBplus3001->GetYaxis()->SetTitle("d#sigma/dp_{T} (B^{+}) (pb GeV-1c)");
   grae->SetHistogram(Graph_gaeSigmaBplus3001);
   
   grae->Draw("5");
   Double_t xAxis1[8] = {5, 7, 10, 15, 20, 30, 50, 100}; 
   
   TH1D *hPtSigma__2 = new TH1D("hPtSigma__2","",7, xAxis1);
   hPtSigma__2->SetBinContent(1,1354179);
   hPtSigma__2->SetBinContent(2,496754.7);
   hPtSigma__2->SetBinContent(3,176196.7);
   hPtSigma__2->SetBinContent(4,61092.56);
   hPtSigma__2->SetBinContent(5,17079.29);
   hPtSigma__2->SetBinContent(6,3456.027);
   hPtSigma__2->SetBinContent(7,174.441);
   hPtSigma__2->SetBinError(1,654957.7);
   hPtSigma__2->SetBinError(2,86364.98);
   hPtSigma__2->SetBinError(3,17112.35);
   hPtSigma__2->SetBinError(4,5379.892);
   hPtSigma__2->SetBinError(5,1434.97);
   hPtSigma__2->SetBinError(6,364.6749);
   hPtSigma__2->SetBinError(7,50.78909);
   hPtSigma__2->SetEntries(10.18334);
   hPtSigma__2->SetLineWidth(2);
   hPtSigma__2->SetMarkerStyle(20);
   hPtSigma__2->SetMarkerSize(0.84);
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
   1354179,
   496754.7,
   176196.7,
   61092.56,
   17079.29,
   3456.027,
   174.441};
   Double_t gaeCrossSyst_felx3002[7] = {
   1,
   1.5,
   2.5,
   2.5,
   5,
   10,
   25};
   Double_t gaeCrossSyst_fely3002[7] = {
   538216.9,
   92076.86,
   27614.28,
   9218.037,
   2527.343,
   510.6905,
   13.63121};
   Double_t gaeCrossSyst_fehx3002[7] = {
   1,
   1.5,
   2.5,
   2.5,
   5,
   10,
   25};
   Double_t gaeCrossSyst_fehy3002[7] = {
   538216.9,
   92076.86,
   27614.28,
   9218.037,
   2527.343,
   510.6905,
   13.63121};
   grae = new TGraphAsymmErrors(7,gaeCrossSyst_fx3002,gaeCrossSyst_fy3002,gaeCrossSyst_felx3002,gaeCrossSyst_fehx3002,gaeCrossSyst_fely3002,gaeCrossSyst_fehy3002);
   grae->SetName("gaeCrossSyst");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);
   grae->SetFillStyle(0);
   grae->SetLineWidth(2);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(0.8);
   
   TH1F *Graph_gaeCrossSyst3002 = new TH1F("Graph_gaeCrossSyst3002","Graph",100,0,109.5);
   Graph_gaeCrossSyst3002->SetMinimum(144.7288);
   Graph_gaeCrossSyst3002->SetMaximum(2081620);
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
   TLatex *   tex = new TLatex(0.18,0.909,"CMS");
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
      tex = new TLatex(0.77,0.874,"B^{+}");
tex->SetNDC();
   tex->SetTextSize(0.063);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.53,0.713,"Global uncert. #plus4.6, #minus5.0%");
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
   entry->SetFillStyle(1001);
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
