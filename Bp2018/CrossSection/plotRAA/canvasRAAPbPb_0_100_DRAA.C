void canvasRAAPbPb_0_100_DRAA()
{
//=========Macro generated from canvas: canvasRAA/canvasRAA
//=========  (Mon Feb 11 09:30:08 2019) by ROOT version 6.14/04
   TCanvas *canvasRAA = new TCanvas("canvasRAA", "canvasRAA",0,0,600,600);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   canvasRAA->SetHighLightColor(2);
   canvasRAA->Range(-14.95882,-0.2310559,113.8647,1.69441);
   canvasRAA->SetFillColor(0);
   canvasRAA->SetBorderMode(0);
   canvasRAA->SetBorderSize(2);
   canvasRAA->SetLeftMargin(0.12);
   canvasRAA->SetRightMargin(0.03);
   canvasRAA->SetTopMargin(0.075);
   canvasRAA->SetBottomMargin(0.12);
   canvasRAA->SetFrameBorderMode(0);
   canvasRAA->SetFrameBorderMode(0);
   
   TH2F *hemptyEff__1 = new TH2F("hemptyEff__1","",50,0.5,110,10,0,1.55);
   hemptyEff__1->SetMinimum(0);
   hemptyEff__1->SetMaximum(2);
   hemptyEff__1->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   hemptyEff__1->SetLineColor(ci);
   hemptyEff__1->SetMarkerStyle(20);
   hemptyEff__1->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   hemptyEff__1->GetXaxis()->CenterTitle(true);
   hemptyEff__1->GetXaxis()->SetLabelFont(42);
   hemptyEff__1->GetXaxis()->SetTitleSize(0.05);
   hemptyEff__1->GetXaxis()->SetTitleFont(42);
   hemptyEff__1->GetYaxis()->SetTitle("R_{AA}");
   hemptyEff__1->GetYaxis()->CenterTitle(true);
   hemptyEff__1->GetYaxis()->SetLabelFont(42);
   hemptyEff__1->GetYaxis()->SetTitleSize(0.05);
   hemptyEff__1->GetYaxis()->SetTitleOffset(1.1);
   hemptyEff__1->GetYaxis()->SetTitleFont(42);
   hemptyEff__1->GetZaxis()->SetLabelFont(42);
   hemptyEff__1->GetZaxis()->SetLabelSize(0.035);
   hemptyEff__1->GetZaxis()->SetTitleSize(0.035);
   hemptyEff__1->GetZaxis()->SetTitleFont(42);
   hemptyEff__1->Draw("");
   TLine *line = new TLine(0.5,1,110,1);
   line->SetLineStyle(2);
   line->SetLineWidth(2);
   line->Draw();
   TLatex *   tex = new TLatex(0.4,0.595,"Centrality 0-100%");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.2,0.936,"28.0 pb^{-1} (pp 5.02 TeV) + 1.5 nb^{-1} (PbPb 5.02 TeV)");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.16,0.76,"B^{+}");
tex->SetNDC();
   tex->SetTextSize(0.08);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.16,0.9,"CMS");
tex->SetNDC();
   tex->SetTextAlign(13);
   tex->SetTextSize(0.06);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.16,0.72,"|y| < 2.4");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
   
   Double_t gNuclearModification_fx3001[5] = {
   22.5,
   27.5,
   35,
   50,
   80};
   Double_t gNuclearModification_fy3001[5] = {
   0.4279886,
   0.4447959,
   0.5041472,
   0.6796131,
   0.7424828};
   Double_t gNuclearModification_felx3001[5] = {
   2.5,
   2.5,
   5,
   10,
   20};
   Double_t gNuclearModification_fely3001[5] = {
   0.08040356,
   0.08777721,
   0.09566393,
   0.1257512,
   0.1437027};
   Double_t gNuclearModification_fehx3001[5] = {
   2.5,
   2.5,
   5,
   10,
   20};
   Double_t gNuclearModification_fehy3001[5] = {
   0.08040356,
   0.08777721,
   0.09566393,
   0.1257512,
   0.1437027};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(5,gNuclearModification_fx3001,gNuclearModification_fy3001,gNuclearModification_felx3001,gNuclearModification_fehx3001,gNuclearModification_fely3001,gNuclearModification_fehy3001);
   grae->SetName("gNuclearModification");
   grae->SetTitle("Graph");

   ci = 1182;
   color = new TColor(ci, 0.6, 1, 0.6, " ", 0.5);
   grae->SetFillColor(ci);
   grae->SetFillStyle(1000);
   grae->SetMarkerStyle(21);
   
   TH1F *Graph_Graph_gNuclearModification30013001 = new TH1F("Graph_Graph_gNuclearModification30013001","Graph",100,12,108);
   Graph_Graph_gNuclearModification30013001->SetMinimum(0.293725);
   Graph_Graph_gNuclearModification30013001->SetMaximum(0.9400455);
   Graph_Graph_gNuclearModification30013001->SetDirectory(0);
   Graph_Graph_gNuclearModification30013001->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_gNuclearModification30013001->SetLineColor(ci);
   Graph_Graph_gNuclearModification30013001->SetMarkerStyle(20);
   Graph_Graph_gNuclearModification30013001->GetXaxis()->SetLabelFont(42);
   Graph_Graph_gNuclearModification30013001->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_gNuclearModification30013001->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_gNuclearModification30013001->GetXaxis()->SetTitleFont(42);
   Graph_Graph_gNuclearModification30013001->GetYaxis()->SetLabelFont(42);
   Graph_Graph_gNuclearModification30013001->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_gNuclearModification30013001->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_gNuclearModification30013001->GetYaxis()->SetTitleOffset(0);
   Graph_Graph_gNuclearModification30013001->GetYaxis()->SetTitleFont(42);
   Graph_Graph_gNuclearModification30013001->GetZaxis()->SetLabelFont(42);
   Graph_Graph_gNuclearModification30013001->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_gNuclearModification30013001->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_gNuclearModification30013001->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_gNuclearModification30013001);
   
   grae->Draw("5");
   
   Double_t gNuclearModification_fx3002[9] = {
   2.5,
   3.5,
   4.5,
   5.5,
   7,
   9,
   11.25,
   13.75,
   17.5};
   Double_t gNuclearModification_fy3002[9] = {
   0.6315195,
   0.4432341,
   0.3674997,
   0.2758138,
   0.2623725,
   0.2722143,
   0.3291268,
   0.3207651,
   0.3937182};
   Double_t gNuclearModification_felx3002[9] = {
   0.5,
   0.5,
   0.5,
   0.5,
   1,
   1,
   1.25,
   1.25,
   2.5};
   Double_t gNuclearModification_fely3002[9] = {
   0.1385446,
   0.08615785,
   0.06750189,
   0.05013722,
   0.04752139,
   0.04884336,
   0.05908503,
   0.05758429,
   0.07618566};
   Double_t gNuclearModification_fehx3002[9] = {
   0.5,
   0.5,
   0.5,
   0.5,
   1,
   1,
   1.25,
   1.25,
   2.5};
   Double_t gNuclearModification_fehy3002[9] = {
   0.1385446,
   0.08615785,
   0.06750189,
   0.05013722,
   0.04752139,
   0.04884336,
   0.05908503,
   0.05758429,
   0.07618566};
   grae = new TGraphAsymmErrors(9,gNuclearModification_fx3002,gNuclearModification_fy3002,gNuclearModification_felx3002,gNuclearModification_fehx3002,gNuclearModification_fely3002,gNuclearModification_fehy3002);
   grae->SetName("gNuclearModification");
   grae->SetTitle("Graph");

   ci = 1183;
   color = new TColor(ci, 0.6, 1, 0.6, " ", 0.5);
   grae->SetFillColor(ci);
   grae->SetFillStyle(1000);
   grae->SetMarkerStyle(21);
   
   TH1F *Graph_Graph_gNuclearModification30023002 = new TH1F("Graph_Graph_gNuclearModification30023002","Graph",100,0.2,21.8);
   Graph_Graph_gNuclearModification30023002->SetMinimum(0.1593298);
   Graph_Graph_gNuclearModification30023002->SetMaximum(0.8255854);
   Graph_Graph_gNuclearModification30023002->SetDirectory(0);
   Graph_Graph_gNuclearModification30023002->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_gNuclearModification30023002->SetLineColor(ci);
   Graph_Graph_gNuclearModification30023002->SetMarkerStyle(20);
   Graph_Graph_gNuclearModification30023002->GetXaxis()->SetLabelFont(42);
   Graph_Graph_gNuclearModification30023002->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_gNuclearModification30023002->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_gNuclearModification30023002->GetXaxis()->SetTitleFont(42);
   Graph_Graph_gNuclearModification30023002->GetYaxis()->SetLabelFont(42);
   Graph_Graph_gNuclearModification30023002->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_gNuclearModification30023002->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_gNuclearModification30023002->GetYaxis()->SetTitleOffset(0);
   Graph_Graph_gNuclearModification30023002->GetYaxis()->SetTitleFont(42);
   Graph_Graph_gNuclearModification30023002->GetZaxis()->SetLabelFont(42);
   Graph_Graph_gNuclearModification30023002->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_gNuclearModification30023002->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_gNuclearModification30023002->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_gNuclearModification30023002);
   
   grae->Draw("5");
   Double_t xAxis1[6] = {20, 25, 30, 40, 60, 100}; 
   
   TH1D *hNuclearModification2__2 = new TH1D("hNuclearModification2__2","",5, xAxis1);
   hNuclearModification2__2->SetBinContent(1,0.4279886);
   hNuclearModification2__2->SetBinContent(2,0.4447959);
   hNuclearModification2__2->SetBinContent(3,0.5041472);
   hNuclearModification2__2->SetBinContent(4,0.6796131);
   hNuclearModification2__2->SetBinContent(5,0.7424828);
   hNuclearModification2__2->SetBinError(1,0.05025093);
   hNuclearModification2__2->SetBinError(2,0.04192271);
   hNuclearModification2__2->SetBinError(3,0.04691724);
   hNuclearModification2__2->SetBinError(4,0.03941607);
   hNuclearModification2__2->SetBinError(5,0.07611912);
   hNuclearModification2__2->SetEntries(566.4226);
   hNuclearModification2__2->SetDirectory(0);
   hNuclearModification2__2->SetStats(0);
   hNuclearModification2__2->SetLineWidth(3);

   ci = TColor::GetColor("#006600");
   hNuclearModification2__2->SetMarkerColor(ci);
   hNuclearModification2__2->SetMarkerStyle(21);
   hNuclearModification2__2->SetMarkerSize(1.2);
   hNuclearModification2__2->GetXaxis()->SetTitle("D^{0} p_{T} (GeV/c)");
   hNuclearModification2__2->GetXaxis()->SetLabelFont(42);
   hNuclearModification2__2->GetXaxis()->SetLabelSize(0.035);
   hNuclearModification2__2->GetXaxis()->SetTitleSize(0.035);
   hNuclearModification2__2->GetXaxis()->SetTitleFont(42);
   hNuclearModification2__2->GetYaxis()->SetTitle("Uncorrected dN(D^{0})/dp_{T}");
   hNuclearModification2__2->GetYaxis()->SetLabelFont(42);
   hNuclearModification2__2->GetYaxis()->SetLabelSize(0.035);
   hNuclearModification2__2->GetYaxis()->SetTitleSize(0.035);
   hNuclearModification2__2->GetYaxis()->SetTitleOffset(0);
   hNuclearModification2__2->GetYaxis()->SetTitleFont(42);
   hNuclearModification2__2->GetZaxis()->SetLabelFont(42);
   hNuclearModification2__2->GetZaxis()->SetLabelSize(0.035);
   hNuclearModification2__2->GetZaxis()->SetTitleSize(0.035);
   hNuclearModification2__2->GetZaxis()->SetTitleFont(42);
   hNuclearModification2__2->Draw("psame");
   Double_t xAxis2[10] = {2, 3, 4, 5, 6, 8, 10, 12.5, 15, 20}; 
   
   TH1D *hNuclearModification3__3 = new TH1D("hNuclearModification3__3","",9, xAxis2);
   hNuclearModification3__3->SetBinContent(1,0.6315195);
   hNuclearModification3__3->SetBinContent(2,0.4432341);
   hNuclearModification3__3->SetBinContent(3,0.3674997);
   hNuclearModification3__3->SetBinContent(4,0.2758138);
   hNuclearModification3__3->SetBinContent(5,0.2623725);
   hNuclearModification3__3->SetBinContent(6,0.2722143);
   hNuclearModification3__3->SetBinContent(7,0.3291268);
   hNuclearModification3__3->SetBinContent(8,0.3207651);
   hNuclearModification3__3->SetBinContent(9,0.3937182);
   hNuclearModification3__3->SetBinError(1,0.1116572);
   hNuclearModification3__3->SetBinError(2,0.04535789);
   hNuclearModification3__3->SetBinError(3,0.03419174);
   hNuclearModification3__3->SetBinError(4,0.02501397);
   hNuclearModification3__3->SetBinError(5,0.01804332);
   hNuclearModification3__3->SetBinError(6,0.01714575);
   hNuclearModification3__3->SetBinError(7,0.02192858);
   hNuclearModification3__3->SetBinError(8,0.03017596);
   hNuclearModification3__3->SetBinError(9,0.04454357);
   hNuclearModification3__3->SetEntries(534.8555);
   hNuclearModification3__3->SetDirectory(0);
   hNuclearModification3__3->SetStats(0);
   hNuclearModification3__3->SetLineWidth(3);

   ci = TColor::GetColor("#006600");
   hNuclearModification3__3->SetMarkerColor(ci);
   hNuclearModification3__3->SetMarkerStyle(21);
   hNuclearModification3__3->SetMarkerSize(1.2);
   hNuclearModification3__3->GetXaxis()->SetTitle("D^{0} p_{T} (GeV/c)");
   hNuclearModification3__3->GetXaxis()->SetLabelFont(42);
   hNuclearModification3__3->GetXaxis()->SetLabelSize(0.035);
   hNuclearModification3__3->GetXaxis()->SetTitleSize(0.035);
   hNuclearModification3__3->GetXaxis()->SetTitleFont(42);
   hNuclearModification3__3->GetYaxis()->SetTitle("Uncorrected dN(D^{0})/dp_{T}");
   hNuclearModification3__3->GetYaxis()->SetLabelFont(42);
   hNuclearModification3__3->GetYaxis()->SetLabelSize(0.035);
   hNuclearModification3__3->GetYaxis()->SetTitleSize(0.035);
   hNuclearModification3__3->GetYaxis()->SetTitleOffset(0);
   hNuclearModification3__3->GetYaxis()->SetTitleFont(42);
   hNuclearModification3__3->GetZaxis()->SetLabelFont(42);
   hNuclearModification3__3->GetZaxis()->SetLabelSize(0.035);
   hNuclearModification3__3->GetZaxis()->SetTitleSize(0.035);
   hNuclearModification3__3->GetZaxis()->SetTitleFont(42);
   hNuclearModification3__3->Draw("psame");
   TBox *box = new TBox(0.62,0.9586116,0.752,1.045662);

   ci = 1185;
   color = new TColor(ci, 0.6, 1, 0.6, " ", 0.5);
   box->SetFillColor(ci);

   ci = TColor::GetColor("#99ff99");
   box->SetLineColor(ci);
   box->Draw();
   
   TLegend *leg = new TLegend(0.3936242,0.6574695,0.812953,0.9157592,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.04);
   leg->SetLineColor(0);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("","D^{0} |y| < 1.0","pf");

   ci = 1184;
   color = new TColor(ci, 0.6, 1, 0.6, " ", 0.5);
   entry->SetFillColor(ci);
   entry->SetFillStyle(1000);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#006600");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry->SetTextSize(0.038);
   entry=leg->AddEntry("gNuclearModification_UnCor","R_{AA}","pe");
   entry->SetLineStyle(1);

   ci = TColor::GetColor("#0033cc");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry->SetTextSize(0.038);
   entry=leg->AddEntry("gNuclearModification","Syst. uncert.","f");

   ci = 1179;
   color = new TColor(ci, 0, 0.6, 1, " ", 0.5);
   entry->SetFillColor(ci);
   entry->SetFillStyle(1000);
   entry->SetLineStyle(1);
   entry->SetMarkerColor(4);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry->SetTextSize(0.03);
   entry=leg->AddEntry("TBox","Global uncert.","f");
   entry->SetFillColor(16);
   entry->SetFillStyle(1001);
   entry->SetLineColor(16);
   entry->SetLineStyle(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry->SetTextSize(0.03);
   leg->Draw();
   
   Double_t gNuclearModification_fx3003[6] = {
   8.5,
   12.5,
   17.5,
   25,
   40,
   75};
   Double_t gNuclearModification_fy3003[6] = {
   0.4515137,
   0.3450858,
   0.4527235,
   0.5063888,
   0.6779442,
   0.580066};
   Double_t gNuclearModification_felx3003[6] = {
   1.5,
   2.5,
   2.5,
   5,
   10,
   25};
   Double_t gNuclearModification_fely3003[6] = {
   0.0878727,
   0.05936694,
   0.07615833,
   0.08366641,
   0.1114905,
   0.09539504};
   Double_t gNuclearModification_fehx3003[6] = {
   1.5,
   2.5,
   2.5,
   5,
   10,
   25};
   Double_t gNuclearModification_fehy3003[6] = {
   0.0878727,
   0.05936694,
   0.07615833,
   0.08366641,
   0.1114905,
   0.09539504};
   grae = new TGraphAsymmErrors(6,gNuclearModification_fx3003,gNuclearModification_fy3003,gNuclearModification_felx3003,gNuclearModification_fehx3003,gNuclearModification_fely3003,gNuclearModification_fehy3003);
   grae->SetName("gNuclearModification");
   grae->SetTitle("Graph");

   ci = 1179;
   color = new TColor(ci, 0, 0.6, 1, " ", 0.5);
   grae->SetFillColor(ci);
   grae->SetFillStyle(1000);
   grae->SetLineColor(0);
   grae->SetLineWidth(0);

   ci = TColor::GetColor("#0033cc");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(21);
   grae->SetMarkerSize(1.2);
   
   TH1F *Graph_gNuclearModification3003 = new TH1F("Graph_gNuclearModification3003","Graph",100,0,109.3);
   Graph_gNuclearModification3003->SetMinimum(0.2353473);
   Graph_gNuclearModification3003->SetMaximum(0.8398063);
   Graph_gNuclearModification3003->SetDirectory(0);
   Graph_gNuclearModification3003->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_gNuclearModification3003->SetLineColor(ci);
   Graph_gNuclearModification3003->SetMarkerStyle(20);
   Graph_gNuclearModification3003->GetXaxis()->SetLabelFont(42);
   Graph_gNuclearModification3003->GetXaxis()->SetLabelSize(0.035);
   Graph_gNuclearModification3003->GetXaxis()->SetTitleSize(0.035);
   Graph_gNuclearModification3003->GetXaxis()->SetTitleFont(42);
   Graph_gNuclearModification3003->GetYaxis()->SetLabelFont(42);
   Graph_gNuclearModification3003->GetYaxis()->SetLabelSize(0.035);
   Graph_gNuclearModification3003->GetYaxis()->SetTitleSize(0.035);
   Graph_gNuclearModification3003->GetYaxis()->SetTitleOffset(0);
   Graph_gNuclearModification3003->GetYaxis()->SetTitleFont(42);
   Graph_gNuclearModification3003->GetZaxis()->SetLabelFont(42);
   Graph_gNuclearModification3003->GetZaxis()->SetLabelSize(0.035);
   Graph_gNuclearModification3003->GetZaxis()->SetTitleSize(0.035);
   Graph_gNuclearModification3003->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_gNuclearModification3003);
   
   grae->Draw("5");
   Double_t xAxis3[7] = {7, 10, 15, 20, 30, 50, 100}; 
   
   TH1D *hNuclearModification__4 = new TH1D("hNuclearModification__4","",6, xAxis3);
   hNuclearModification__4->SetBinContent(1,0.4515137);
   hNuclearModification__4->SetBinContent(2,0.3450858);
   hNuclearModification__4->SetBinContent(3,0.4527235);
   hNuclearModification__4->SetBinContent(4,0.5063888);
   hNuclearModification__4->SetBinContent(5,0.6779442);
   hNuclearModification__4->SetBinContent(6,0.580066);
   hNuclearModification__4->SetBinError(1,0.09450069);
   hNuclearModification__4->SetBinError(2,0.03605466);
   hNuclearModification__4->SetBinError(3,0.04235486);
   hNuclearModification__4->SetBinError(4,0.0475547);
   hNuclearModification__4->SetBinError(5,0.09520653);
   hNuclearModification__4->SetBinError(6,0.156416);
   hNuclearModification__4->SetEntries(189.9475);

   ci = TColor::GetColor("#0033cc");
   hNuclearModification__4->SetLineColor(ci);
   hNuclearModification__4->SetLineWidth(3);

   ci = TColor::GetColor("#0033cc");
   hNuclearModification__4->SetMarkerColor(ci);
   hNuclearModification__4->SetMarkerStyle(21);
   hNuclearModification__4->SetMarkerSize(1.2);
   hNuclearModification__4->GetXaxis()->SetTitle("D^{0} p_{T} (GeV/c)");
   hNuclearModification__4->GetXaxis()->SetLabelFont(42);
   hNuclearModification__4->GetXaxis()->SetLabelSize(0.035);
   hNuclearModification__4->GetXaxis()->SetTitleSize(0.035);
   hNuclearModification__4->GetXaxis()->SetTitleFont(42);
   hNuclearModification__4->GetYaxis()->SetTitle("Uncorrected dN(D^{0})/dp_{T}");
   hNuclearModification__4->GetYaxis()->SetLabelFont(42);
   hNuclearModification__4->GetYaxis()->SetLabelSize(0.035);
   hNuclearModification__4->GetYaxis()->SetTitleSize(0.035);
   hNuclearModification__4->GetYaxis()->SetTitleOffset(0);
   hNuclearModification__4->GetYaxis()->SetTitleFont(42);
   hNuclearModification__4->GetZaxis()->SetLabelFont(42);
   hNuclearModification__4->GetZaxis()->SetLabelSize(0.035);
   hNuclearModification__4->GetZaxis()->SetTitleSize(0.035);
   hNuclearModification__4->GetZaxis()->SetTitleFont(42);
   hNuclearModification__4->Draw("same");
   box = new TBox(0.62,0.9586116,0.752,1.045662);

   ci = 1185;
   color = new TColor(ci, 0.6, 1, 0.6, " ", 0.5);
   box->SetFillColor(ci);

   ci = TColor::GetColor("#99ff99");
   box->SetLineColor(ci);
   box->Draw();
   Double_t xAxis4[6] = {20, 25, 30, 40, 60, 100}; 
   
   TH1D *hNuclearModification2__5 = new TH1D("hNuclearModification2__5","",5, xAxis4);
   hNuclearModification2__5->SetBinContent(1,0.4279886);
   hNuclearModification2__5->SetBinContent(2,0.4447959);
   hNuclearModification2__5->SetBinContent(3,0.5041472);
   hNuclearModification2__5->SetBinContent(4,0.6796131);
   hNuclearModification2__5->SetBinContent(5,0.7424828);
   hNuclearModification2__5->SetBinError(1,0.05025093);
   hNuclearModification2__5->SetBinError(2,0.04192271);
   hNuclearModification2__5->SetBinError(3,0.04691724);
   hNuclearModification2__5->SetBinError(4,0.03941607);
   hNuclearModification2__5->SetBinError(5,0.07611912);
   hNuclearModification2__5->SetEntries(566.4226);
   hNuclearModification2__5->SetStats(0);
   hNuclearModification2__5->SetLineWidth(3);

   ci = TColor::GetColor("#006600");
   hNuclearModification2__5->SetMarkerColor(ci);
   hNuclearModification2__5->SetMarkerStyle(21);
   hNuclearModification2__5->SetMarkerSize(1.2);
   hNuclearModification2__5->GetXaxis()->SetTitle("D^{0} p_{T} (GeV/c)");
   hNuclearModification2__5->GetXaxis()->SetLabelFont(42);
   hNuclearModification2__5->GetXaxis()->SetLabelSize(0.035);
   hNuclearModification2__5->GetXaxis()->SetTitleSize(0.035);
   hNuclearModification2__5->GetXaxis()->SetTitleFont(42);
   hNuclearModification2__5->GetYaxis()->SetTitle("Uncorrected dN(D^{0})/dp_{T}");
   hNuclearModification2__5->GetYaxis()->SetLabelFont(42);
   hNuclearModification2__5->GetYaxis()->SetLabelSize(0.035);
   hNuclearModification2__5->GetYaxis()->SetTitleSize(0.035);
   hNuclearModification2__5->GetYaxis()->SetTitleOffset(0);
   hNuclearModification2__5->GetYaxis()->SetTitleFont(42);
   hNuclearModification2__5->GetZaxis()->SetLabelFont(42);
   hNuclearModification2__5->GetZaxis()->SetLabelSize(0.035);
   hNuclearModification2__5->GetZaxis()->SetTitleSize(0.035);
   hNuclearModification2__5->GetZaxis()->SetTitleFont(42);
   hNuclearModification2__5->Draw("psame");
   Double_t xAxis5[10] = {2, 3, 4, 5, 6, 8, 10, 12.5, 15, 20}; 
   
   TH1D *hNuclearModification3__6 = new TH1D("hNuclearModification3__6","",9, xAxis5);
   hNuclearModification3__6->SetBinContent(1,0.6315195);
   hNuclearModification3__6->SetBinContent(2,0.4432341);
   hNuclearModification3__6->SetBinContent(3,0.3674997);
   hNuclearModification3__6->SetBinContent(4,0.2758138);
   hNuclearModification3__6->SetBinContent(5,0.2623725);
   hNuclearModification3__6->SetBinContent(6,0.2722143);
   hNuclearModification3__6->SetBinContent(7,0.3291268);
   hNuclearModification3__6->SetBinContent(8,0.3207651);
   hNuclearModification3__6->SetBinContent(9,0.3937182);
   hNuclearModification3__6->SetBinError(1,0.1116572);
   hNuclearModification3__6->SetBinError(2,0.04535789);
   hNuclearModification3__6->SetBinError(3,0.03419174);
   hNuclearModification3__6->SetBinError(4,0.02501397);
   hNuclearModification3__6->SetBinError(5,0.01804332);
   hNuclearModification3__6->SetBinError(6,0.01714575);
   hNuclearModification3__6->SetBinError(7,0.02192858);
   hNuclearModification3__6->SetBinError(8,0.03017596);
   hNuclearModification3__6->SetBinError(9,0.04454357);
   hNuclearModification3__6->SetEntries(534.8555);
   hNuclearModification3__6->SetStats(0);
   hNuclearModification3__6->SetLineWidth(3);

   ci = TColor::GetColor("#006600");
   hNuclearModification3__6->SetMarkerColor(ci);
   hNuclearModification3__6->SetMarkerStyle(21);
   hNuclearModification3__6->SetMarkerSize(1.2);
   hNuclearModification3__6->GetXaxis()->SetTitle("D^{0} p_{T} (GeV/c)");
   hNuclearModification3__6->GetXaxis()->SetLabelFont(42);
   hNuclearModification3__6->GetXaxis()->SetLabelSize(0.035);
   hNuclearModification3__6->GetXaxis()->SetTitleSize(0.035);
   hNuclearModification3__6->GetXaxis()->SetTitleFont(42);
   hNuclearModification3__6->GetYaxis()->SetTitle("Uncorrected dN(D^{0})/dp_{T}");
   hNuclearModification3__6->GetYaxis()->SetLabelFont(42);
   hNuclearModification3__6->GetYaxis()->SetLabelSize(0.035);
   hNuclearModification3__6->GetYaxis()->SetTitleSize(0.035);
   hNuclearModification3__6->GetYaxis()->SetTitleOffset(0);
   hNuclearModification3__6->GetYaxis()->SetTitleFont(42);
   hNuclearModification3__6->GetZaxis()->SetLabelFont(42);
   hNuclearModification3__6->GetZaxis()->SetLabelSize(0.035);
   hNuclearModification3__6->GetZaxis()->SetTitleSize(0.035);
   hNuclearModification3__6->GetZaxis()->SetTitleFont(42);
   hNuclearModification3__6->Draw("psame");
   Double_t xAxis6[7] = {7, 10, 15, 20, 30, 50, 100}; 
   
   TH1D *hNuclearModification__7 = new TH1D("hNuclearModification__7","",6, xAxis6);
   hNuclearModification__7->SetBinContent(1,0.4515137);
   hNuclearModification__7->SetBinContent(2,0.3450858);
   hNuclearModification__7->SetBinContent(3,0.4527235);
   hNuclearModification__7->SetBinContent(4,0.5063888);
   hNuclearModification__7->SetBinContent(5,0.6779442);
   hNuclearModification__7->SetBinContent(6,0.580066);
   hNuclearModification__7->SetBinError(1,0.09450069);
   hNuclearModification__7->SetBinError(2,0.03605466);
   hNuclearModification__7->SetBinError(3,0.04235486);
   hNuclearModification__7->SetBinError(4,0.0475547);
   hNuclearModification__7->SetBinError(5,0.09520653);
   hNuclearModification__7->SetBinError(6,0.156416);
   hNuclearModification__7->SetEntries(189.9475);

   ci = TColor::GetColor("#0033cc");
   hNuclearModification__7->SetLineColor(ci);
   hNuclearModification__7->SetLineWidth(3);

   ci = TColor::GetColor("#0033cc");
   hNuclearModification__7->SetMarkerColor(ci);
   hNuclearModification__7->SetMarkerStyle(21);
   hNuclearModification__7->SetMarkerSize(1.2);
   hNuclearModification__7->GetXaxis()->SetTitle("D^{0} p_{T} (GeV/c)");
   hNuclearModification__7->GetXaxis()->SetLabelFont(42);
   hNuclearModification__7->GetXaxis()->SetLabelSize(0.035);
   hNuclearModification__7->GetXaxis()->SetTitleSize(0.035);
   hNuclearModification__7->GetXaxis()->SetTitleFont(42);
   hNuclearModification__7->GetYaxis()->SetTitle("Uncorrected dN(D^{0})/dp_{T}");
   hNuclearModification__7->GetYaxis()->SetLabelFont(42);
   hNuclearModification__7->GetYaxis()->SetLabelSize(0.035);
   hNuclearModification__7->GetYaxis()->SetTitleSize(0.035);
   hNuclearModification__7->GetYaxis()->SetTitleOffset(0);
   hNuclearModification__7->GetYaxis()->SetTitleFont(42);
   hNuclearModification__7->GetZaxis()->SetLabelFont(42);
   hNuclearModification__7->GetZaxis()->SetLabelSize(0.035);
   hNuclearModification__7->GetZaxis()->SetTitleSize(0.035);
   hNuclearModification__7->GetZaxis()->SetTitleFont(42);
   hNuclearModification__7->Draw("same p");
      tex = new TLatex(0.16,0.76,"B^{+}");
tex->SetNDC();
   tex->SetTextSize(0.08);
   tex->SetLineWidth(2);
   tex->Draw();
   
   TH2F *hemptyEff_copy__8 = new TH2F("hemptyEff_copy__8","",50,0.5,110,10,0,1.55);
   hemptyEff_copy__8->SetMinimum(0);
   hemptyEff_copy__8->SetMaximum(2);
   hemptyEff_copy__8->SetDirectory(0);
   hemptyEff_copy__8->SetStats(0);

   ci = TColor::GetColor("#000099");
   hemptyEff_copy__8->SetLineColor(ci);
   hemptyEff_copy__8->SetMarkerStyle(20);
   hemptyEff_copy__8->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   hemptyEff_copy__8->GetXaxis()->CenterTitle(true);
   hemptyEff_copy__8->GetXaxis()->SetLabelFont(42);
   hemptyEff_copy__8->GetXaxis()->SetTitleSize(0.05);
   hemptyEff_copy__8->GetXaxis()->SetTitleFont(42);
   hemptyEff_copy__8->GetYaxis()->SetTitle("R_{AA}");
   hemptyEff_copy__8->GetYaxis()->CenterTitle(true);
   hemptyEff_copy__8->GetYaxis()->SetLabelFont(42);
   hemptyEff_copy__8->GetYaxis()->SetTitleSize(0.05);
   hemptyEff_copy__8->GetYaxis()->SetTitleOffset(1.1);
   hemptyEff_copy__8->GetYaxis()->SetTitleFont(42);
   hemptyEff_copy__8->GetZaxis()->SetLabelFont(42);
   hemptyEff_copy__8->GetZaxis()->SetLabelSize(0.035);
   hemptyEff_copy__8->GetZaxis()->SetTitleSize(0.035);
   hemptyEff_copy__8->GetZaxis()->SetTitleFont(42);
   hemptyEff_copy__8->Draw("sameaxis");
   canvasRAA->Modified();
   canvasRAA->cd();
   canvasRAA->SetSelected(canvasRAA);
}
