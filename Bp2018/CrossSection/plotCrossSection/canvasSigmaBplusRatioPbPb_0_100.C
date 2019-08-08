void canvasSigmaBplusRatioPbPb_0_100()
{
//=========Macro generated from canvas: cSigma/
//=========  (Thu Feb  7 19:15:19 2019) by ROOT version 6.14/04
   TCanvas *cSigma = new TCanvas("cSigma", "",0,0,600,600);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   cSigma->SetHighLightColor(2);
   cSigma->Range(-16.6801,1.47101,112.005,7.348675);
   cSigma->SetFillColor(0);
   cSigma->SetBorderMode(0);
   cSigma->SetBorderSize(2);
   cSigma->SetLogy();
   cSigma->SetLeftMargin(0.1451613);
   cSigma->SetRightMargin(0.05443548);
   cSigma->SetTopMargin(0.05932203);
   cSigma->SetFrameBorderMode(0);
   cSigma->SetFrameBorderMode(0);
   
   TH2F *hemptySigma__1 = new TH2F("hemptySigma__1","",50,2,105,10,100,1e+07);
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
   Double_t xAxis1[7] = {7, 10, 15, 20, 30, 50, 100}; 
   
   TH1D *hPtSigma__2 = new TH1D("hPtSigma__2","",6, xAxis1);
   hPtSigma__2->SetBinContent(1,1176867);
   hPtSigma__2->SetBinContent(2,256846.9);
   hPtSigma__2->SetBinContent(3,89175.6);
   hPtSigma__2->SetBinContent(4,23541.95);
   hPtSigma__2->SetBinContent(5,3621.783);
   hPtSigma__2->SetBinContent(6,176.3634);
   hPtSigma__2->SetBinError(1,234479.1);
   hPtSigma__2->SetBinError(2,24888.58);
   hPtSigma__2->SetBinError(3,7272.709);
   hPtSigma__2->SetBinError(4,1833.6);
   hPtSigma__2->SetBinError(5,382.1873);
   hPtSigma__2->SetBinError(6,25.0572);
   hPtSigma__2->SetEntries(43.17955);
   hPtSigma__2->SetLineWidth(2);
   hPtSigma__2->SetMarkerStyle(20);
   hPtSigma__2->SetMarkerSize(0.84);
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
   
   Double_t gaeCrossSyst_fx3001[6] = {
   8.5,
   12.5,
   17.5,
   25,
   40,
   6.95291e-310};
   Double_t gaeCrossSyst_fy3001[6] = {
   1176867,
   256846.9,
   89175.6,
   23541.95,
   3621.783,
   176.3634};
   Double_t gaeCrossSyst_felx3001[6] = {
   1.5,
   2.5,
   2.5,
   5,
   10,
   4.001932e-322};
   Double_t gaeCrossSyst_fely3001[6] = {
   203120.8,
   39147.86,
   13327.83,
   3453.615,
   529.6768,
   -nan};
   Double_t gaeCrossSyst_fehx3001[6] = {
   1.5,
   2.5,
   2.5,
   5,
   10,
   2.420922e-322};
   Double_t gaeCrossSyst_fehy3001[6] = {
   203120.8,
   39147.86,
   13327.83,
   3453.615,
   529.6768,
   -nan};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(6,gaeCrossSyst_fx3001,gaeCrossSyst_fy3001,gaeCrossSyst_felx3001,gaeCrossSyst_fehx3001,gaeCrossSyst_fely3001,gaeCrossSyst_fehy3001);
   grae->SetName("gaeCrossSyst");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);
   grae->SetFillStyle(0);
   grae->SetLineWidth(2);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(0.8);
   
   TH1F *Graph_gaeCrossSyst3001 = new TH1F("Graph_gaeCrossSyst3001","Graph",100,0,55);
   Graph_gaeCrossSyst3001->SetMinimum(158.7271);
   Graph_gaeCrossSyst3001->SetMaximum(1517969);
   Graph_gaeCrossSyst3001->SetDirectory(0);
   Graph_gaeCrossSyst3001->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_gaeCrossSyst3001->SetLineColor(ci);
   Graph_gaeCrossSyst3001->SetMarkerStyle(20);
   Graph_gaeCrossSyst3001->GetXaxis()->SetLabelFont(42);
   Graph_gaeCrossSyst3001->GetXaxis()->SetLabelSize(0.035);
   Graph_gaeCrossSyst3001->GetXaxis()->SetTitleSize(0.035);
   Graph_gaeCrossSyst3001->GetXaxis()->SetTitleFont(42);
   Graph_gaeCrossSyst3001->GetYaxis()->SetLabelFont(42);
   Graph_gaeCrossSyst3001->GetYaxis()->SetLabelSize(0.035);
   Graph_gaeCrossSyst3001->GetYaxis()->SetTitleSize(0.035);
   Graph_gaeCrossSyst3001->GetYaxis()->SetTitleOffset(0);
   Graph_gaeCrossSyst3001->GetYaxis()->SetTitleFont(42);
   Graph_gaeCrossSyst3001->GetZaxis()->SetLabelFont(42);
   Graph_gaeCrossSyst3001->GetZaxis()->SetLabelSize(0.035);
   Graph_gaeCrossSyst3001->GetZaxis()->SetTitleSize(0.035);
   Graph_gaeCrossSyst3001->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_gaeCrossSyst3001);
   
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
