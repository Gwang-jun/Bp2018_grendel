void canvasRAAPbPb_0_100_SepUnc_DRAA_ChHadRAA_NPJpsiRAA()
{
//=========Macro generated from canvas: canvasRAA/canvasRAA
//=========  (Thu Jan 31 12:37:06 2019) by ROOT version 6.14/04
   TCanvas *canvasRAA = new TCanvas("canvasRAA", "canvasRAA",0,0,600,600);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   canvasRAA->SetHighLightColor(2);
   canvasRAA->Range(-0.7357379,-0.2310559,2.886828,1.69441);
   canvasRAA->SetFillColor(0);
   canvasRAA->SetBorderMode(0);
   canvasRAA->SetBorderSize(2);
   canvasRAA->SetLogx();
   canvasRAA->SetLeftMargin(0.12);
   canvasRAA->SetRightMargin(0.03);
   canvasRAA->SetTopMargin(0.075);
   canvasRAA->SetBottomMargin(0.12);
   canvasRAA->SetFrameBorderMode(0);
   canvasRAA->SetFrameBorderMode(0);
   
   TH2F *hemptyEff__1 = new TH2F("hemptyEff__1","",50,0.5,600,10,0,1.55);
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
   TLine *line = new TLine(0.5,1,600,1);
   line->SetLineStyle(2);
   line->SetLineWidth(2);
   line->Draw();
   TLatex *   tex = new TLatex(0.2,0.936,"28.0 pb^{-1} (pp 5.02 TeV) + 351 #mub^{-1} (PbPb 5.02 TeV)");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.038);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.81,0.21,"B^{#pm}");
tex->SetNDC();
   tex->SetTextSize(0.08);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.15,0.9,"CMS");
tex->SetNDC();
   tex->SetTextAlign(13);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.81,0.16,"|y| < 2.4");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.15,0.845,"Supplementary");
tex->SetNDC();
   tex->SetTextAlign(13);
   tex->SetTextFont(52);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
   TBox *box = new TBox(0.7,0.3322697,0.8,0.3887397);

   ci = TColor::GetColor("#ffcc00");
   box->SetFillColor(ci);
   box->Draw();
   box = new TBox(0.8,0.3547214,0.9,0.4153951);

   ci = TColor::GetColor("#ffcc00");
   box->SetFillColor(ci);
   box->Draw();
   box = new TBox(0.9,0.3856854,1,0.4519942);

   ci = TColor::GetColor("#ffcc00");
   box->SetFillColor(ci);
   box->Draw();
   box = new TBox(1,0.4035821,1.1,0.4786815);

   ci = TColor::GetColor("#ffcc00");
   box->SetFillColor(ci);
   box->Draw();
   box = new TBox(1.1,0.4178419,1.2,0.4951344);

   ci = TColor::GetColor("#ffcc00");
   box->SetFillColor(ci);
   box->Draw();
   box = new TBox(1.2,0.4336979,1.4,0.5260356);

   ci = TColor::GetColor("#ffcc00");
   box->SetFillColor(ci);
   box->Draw();
   box = new TBox(1.4,0.4471445,1.6,0.5499409);

   ci = TColor::GetColor("#ffcc00");
   box->SetFillColor(ci);
   box->Draw();
   box = new TBox(1.6,0.454559,1.8,0.5776292);

   ci = TColor::GetColor("#ffcc00");
   box->SetFillColor(ci);
   box->Draw();
   box = new TBox(1.8,0.4545544,2,0.5775976);

   ci = TColor::GetColor("#ffcc00");
   box->SetFillColor(ci);
   box->Draw();
   box = new TBox(2,0.4458492,2.2,0.5856694);

   ci = TColor::GetColor("#ffcc00");
   box->SetFillColor(ci);
   box->Draw();
   box = new TBox(2.2,0.4400381,2.4,0.5780321);

   ci = TColor::GetColor("#ffcc00");
   box->SetFillColor(ci);
   box->Draw();
   box = new TBox(2.4,0.4004238,3.2,0.5432211);

   ci = TColor::GetColor("#ffcc00");
   box->SetFillColor(ci);
   box->Draw();
   box = new TBox(3.2,0.3194642,4,0.443112);

   ci = TColor::GetColor("#ffcc00");
   box->SetFillColor(ci);
   box->Draw();
   box = new TBox(4,0.2555583,4.8,0.3477606);

   ci = TColor::GetColor("#ffcc00");
   box->SetFillColor(ci);
   box->Draw();
   box = new TBox(4.8,0.2236467,5.6,0.290029);

   ci = TColor::GetColor("#ffcc00");
   box->SetFillColor(ci);
   box->Draw();
   box = new TBox(5.6,0.2091373,6.4,0.2633415);

   ci = TColor::GetColor("#ffcc00");
   box->SetFillColor(ci);
   box->Draw();
   box = new TBox(6.4,0.2084873,7.2,0.2542985);

   ci = TColor::GetColor("#ffcc00");
   box->SetFillColor(ci);
   box->Draw();
   box = new TBox(7.2,0.218131,9.6,0.2613273);

   ci = TColor::GetColor("#ffcc00");
   box->SetFillColor(ci);
   box->Draw();
   box = new TBox(9.6,0.254349,12,0.3047338);

   ci = TColor::GetColor("#ffcc00");
   box->SetFillColor(ci);
   box->Draw();
   box = new TBox(12,0.279773,14.4,0.3421784);

   ci = TColor::GetColor("#ffcc00");
   box->SetFillColor(ci);
   box->Draw();
   box = new TBox(14.4,0.3230158,19.2,0.395078);

   ci = TColor::GetColor("#ffcc00");
   box->SetFillColor(ci);
   box->Draw();
   box = new TBox(19.2,0.3691498,24,0.4514966);

   ci = TColor::GetColor("#ffcc00");
   box->SetFillColor(ci);
   box->Draw();
   box = new TBox(24,0.4081148,28.8,0.4990331);

   ci = TColor::GetColor("#ffcc00");
   box->SetFillColor(ci);
   box->Draw();
   box = new TBox(28.8,0.4382444,35.2,0.5359827);

   ci = TColor::GetColor("#ffcc00");
   box->SetFillColor(ci);
   box->Draw();
   box = new TBox(35.2,0.4849575,41.6,0.5931196);

   ci = TColor::GetColor("#ffcc00");
   box->SetFillColor(ci);
   box->Draw();
   box = new TBox(41.6,0.54671,48,0.6660902);

   ci = TColor::GetColor("#ffcc00");
   box->SetFillColor(ci);
   box->Draw();
   box = new TBox(48,0.5934783,60.8,0.7232848);

   ci = TColor::GetColor("#ffcc00");
   box->SetFillColor(ci);
   box->Draw();
   box = new TBox(60.8,0.6229389,73.6,0.7600401);

   ci = TColor::GetColor("#ffcc00");
   box->SetFillColor(ci);
   box->Draw();
   box = new TBox(73.6,0.6765564,86.4,0.8254808);

   ci = TColor::GetColor("#ffcc00");
   box->SetFillColor(ci);
   box->Draw();
   box = new TBox(86.4,0.7227269,103.6,0.8827653);

   ci = TColor::GetColor("#ffcc00");
   box->SetFillColor(ci);
   box->Draw();
   box = new TBox(103.6,0.7255117,120.8,0.8873466);

   ci = TColor::GetColor("#ffcc00");
   box->SetFillColor(ci);
   box->Draw();
   box = new TBox(120.8,0.732652,140,0.8963203);

   ci = TColor::GetColor("#ffcc00");
   box->SetFillColor(ci);
   box->Draw();
   box = new TBox(140,0.7948418,165,0.9722788);

   ci = TColor::GetColor("#ffcc00");
   box->SetFillColor(ci);
   box->Draw();
   box = new TBox(165,0.8158651,250,0.9995696);

   ci = TColor::GetColor("#ffcc00");
   box->SetFillColor(ci);
   box->Draw();
   box = new TBox(250,0.7867229,400,0.9667282);

   ci = TColor::GetColor("#ffcc00");
   box->SetFillColor(ci);
   box->Draw();
   Double_t xAxis1[38] = {0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 3.2, 4, 4.8, 5.6, 6.4, 7.2, 9.6, 12, 14.4, 19.2, 24, 28.8, 35.2, 41.6, 48, 60.8, 73.6, 86.4, 103.6, 120.8, 140, 165, 250, 400}; 
   
   TH1D *RAA_0_100584__2__2__2__2 = new TH1D("RAA_0_100584__2__2__2__2","",37, xAxis1);
   RAA_0_100584__2__2__2__2->SetBinContent(1,-1);
   RAA_0_100584__2__2__2__2->SetBinContent(2,-1);
   RAA_0_100584__2__2__2__2->SetBinContent(3,0.3605047);
   RAA_0_100584__2__2__2__2->SetBinContent(4,0.3850582);
   RAA_0_100584__2__2__2__2->SetBinContent(5,0.4188398);
   RAA_0_100584__2__2__2__2->SetBinContent(6,0.4411318);
   RAA_0_100584__2__2__2__2->SetBinContent(7,0.4564881);
   RAA_0_100584__2__2__2__2->SetBinContent(8,0.4798668);
   RAA_0_100584__2__2__2__2->SetBinContent(9,0.4985427);
   RAA_0_100584__2__2__2__2->SetBinContent(10,0.516094);
   RAA_0_100584__2__2__2__2->SetBinContent(11,0.5160761);
   RAA_0_100584__2__2__2__2->SetBinContent(12,0.5157593);
   RAA_0_100584__2__2__2__2->SetBinContent(13,0.5090351);
   RAA_0_100584__2__2__2__2->SetBinContent(14,0.4718224);
   RAA_0_100584__2__2__2__2->SetBinContent(15,0.381288);
   RAA_0_100584__2__2__2__2->SetBinContent(16,0.3016595);
   RAA_0_100584__2__2__2__2->SetBinContent(17,0.2568378);
   RAA_0_100584__2__2__2__2->SetBinContent(18,0.2362394);
   RAA_0_100584__2__2__2__2->SetBinContent(19,0.2313929);
   RAA_0_100584__2__2__2__2->SetBinContent(20,0.2397292);
   RAA_0_100584__2__2__2__2->SetBinContent(21,0.2795414);
   RAA_0_100584__2__2__2__2->SetBinContent(22,0.3109757);
   RAA_0_100584__2__2__2__2->SetBinContent(23,0.359047);
   RAA_0_100584__2__2__2__2->SetBinContent(24,0.4103232);
   RAA_0_100584__2__2__2__2->SetBinContent(25,0.4535739);
   RAA_0_100584__2__2__2__2->SetBinContent(26,0.4871136);
   RAA_0_100584__2__2__2__2->SetBinContent(27,0.5390386);
   RAA_0_100584__2__2__2__2->SetBinContent(28,0.6064001);
   RAA_0_100584__2__2__2__2->SetBinContent(29,0.6583816);
   RAA_0_100584__2__2__2__2->SetBinContent(30,0.6914895);
   RAA_0_100584__2__2__2__2->SetBinContent(31,0.7510186);
   RAA_0_100584__2__2__2__2->SetBinContent(32,0.8027461);
   RAA_0_100584__2__2__2__2->SetBinContent(33,0.8064291);
   RAA_0_100584__2__2__2__2->SetBinContent(34,0.8144861);
   RAA_0_100584__2__2__2__2->SetBinContent(35,0.8835603);
   RAA_0_100584__2__2__2__2->SetBinContent(36,0.9077175);
   RAA_0_100584__2__2__2__2->SetBinContent(37,0.8767255);
   RAA_0_100584__2__2__2__2->SetBinError(1,1.279903e-05);
   RAA_0_100584__2__2__2__2->SetBinError(2,1.387861e-05);
   RAA_0_100584__2__2__2__2->SetBinError(3,1.439234e-05);
   RAA_0_100584__2__2__2__2->SetBinError(4,1.673993e-05);
   RAA_0_100584__2__2__2__2->SetBinError(5,1.996438e-05);
   RAA_0_100584__2__2__2__2->SetBinError(6,2.296277e-05);
   RAA_0_100584__2__2__2__2->SetBinError(7,2.652961e-05);
   RAA_0_100584__2__2__2__2->SetBinError(8,2.303589e-05);
   RAA_0_100584__2__2__2__2->SetBinError(9,2.956636e-05);
   RAA_0_100584__2__2__2__2->SetBinError(10,3.708897e-05);
   RAA_0_100584__2__2__2__2->SetBinError(11,4.504344e-05);
   RAA_0_100584__2__2__2__2->SetBinError(12,5.385318e-05);
   RAA_0_100584__2__2__2__2->SetBinError(13,6.332594e-05);
   RAA_0_100584__2__2__2__2->SetBinError(14,4.278537e-05);
   RAA_0_100584__2__2__2__2->SetBinError(15,6.266848e-05);
   RAA_0_100584__2__2__2__2->SetBinError(16,8.299302e-05);
   RAA_0_100584__2__2__2__2->SetBinError(17,0.0001099248);
   RAA_0_100584__2__2__2__2->SetBinError(18,0.0001477243);
   RAA_0_100584__2__2__2__2->SetBinError(19,0.0002002684);
   RAA_0_100584__2__2__2__2->SetBinError(20,0.0001962782);
   RAA_0_100584__2__2__2__2->SetBinError(21,0.00043599);
   RAA_0_100584__2__2__2__2->SetBinError(22,0.0007955915);
   RAA_0_100584__2__2__2__2->SetBinError(23,0.001121086);
   RAA_0_100584__2__2__2__2->SetBinError(24,0.002349183);
   RAA_0_100584__2__2__2__2->SetBinError(25,0.003934258);
   RAA_0_100584__2__2__2__2->SetBinError(26,0.005431548);
   RAA_0_100584__2__2__2__2->SetBinError(27,0.007593713);
   RAA_0_100584__2__2__2__2->SetBinError(28,0.008532134);
   RAA_0_100584__2__2__2__2->SetBinError(29,0.006423681);
   RAA_0_100584__2__2__2__2->SetBinError(30,0.009668077);
   RAA_0_100584__2__2__2__2->SetBinError(31,0.009013466);
   RAA_0_100584__2__2__2__2->SetBinError(32,0.01230905);
   RAA_0_100584__2__2__2__2->SetBinError(33,0.01985833);
   RAA_0_100584__2__2__2__2->SetBinError(34,0.02497965);
   RAA_0_100584__2__2__2__2->SetBinError(35,0.04145952);
   RAA_0_100584__2__2__2__2->SetBinError(36,0.0459072);
   RAA_0_100584__2__2__2__2->SetBinError(37,0.152905);
   RAA_0_100584__2__2__2__2->SetMinimum(0);
   RAA_0_100584__2__2__2__2->SetMaximum(1.8);
   RAA_0_100584__2__2__2__2->SetEntries(11078.9);
   RAA_0_100584__2__2__2__2->SetDirectory(0);
   RAA_0_100584__2__2__2__2->SetStats(0);

   ci = TColor::GetColor("#ffcc00");
   RAA_0_100584__2__2__2__2->SetFillColor(ci);

   ci = TColor::GetColor("#000099");
   RAA_0_100584__2__2__2__2->SetLineColor(ci);
   RAA_0_100584__2__2__2__2->SetLineWidth(2);
   RAA_0_100584__2__2__2__2->SetMarkerStyle(20);
   RAA_0_100584__2__2__2__2->SetMarkerSize(1.3);
   RAA_0_100584__2__2__2__2->GetXaxis()->SetTitle("p_{T} (GeV)");
   RAA_0_100584__2__2__2__2->GetXaxis()->SetRange(1,37);
   RAA_0_100584__2__2__2__2->GetXaxis()->CenterTitle(true);
   RAA_0_100584__2__2__2__2->GetXaxis()->SetLabelFont(42);
   RAA_0_100584__2__2__2__2->GetXaxis()->SetLabelOffset(-0.002);
   RAA_0_100584__2__2__2__2->GetXaxis()->SetLabelSize(0.045);
   RAA_0_100584__2__2__2__2->GetXaxis()->SetTitleSize(0.05);
   RAA_0_100584__2__2__2__2->GetXaxis()->SetTitleFont(42);
   RAA_0_100584__2__2__2__2->GetYaxis()->SetTitle("R_{AA}, R_{pA}");
   RAA_0_100584__2__2__2__2->GetYaxis()->CenterTitle(true);
   RAA_0_100584__2__2__2__2->GetYaxis()->SetNdivisions(508);
   RAA_0_100584__2__2__2__2->GetYaxis()->SetLabelFont(42);
   RAA_0_100584__2__2__2__2->GetYaxis()->SetLabelSize(0.045);
   RAA_0_100584__2__2__2__2->GetYaxis()->SetTitleSize(0.05);
   RAA_0_100584__2__2__2__2->GetYaxis()->SetTitleOffset(0);
   RAA_0_100584__2__2__2__2->GetYaxis()->SetTitleFont(42);
   RAA_0_100584__2__2__2__2->GetZaxis()->SetLabelFont(42);
   RAA_0_100584__2__2__2__2->GetZaxis()->SetLabelSize(0.05);
   RAA_0_100584__2__2__2__2->GetZaxis()->SetTitleSize(0.09);
   RAA_0_100584__2__2__2__2->GetZaxis()->SetTitleFont(42);
   RAA_0_100584__2__2__2__2->Draw("same");
   box = new TBox(0.5,0.9080761,0.62,1.091924);

   ci = TColor::GetColor("#ffcc00");
   box->SetFillColor(ci);

   ci = TColor::GetColor("#ffcc00");
   box->SetLineColor(ci);
   box->Draw();
   
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

   ci = 1183;
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

   ci = 1184;
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
   Double_t xAxis2[6] = {20, 25, 30, 40, 60, 100}; 
   
   TH1D *hNuclearModification2__3 = new TH1D("hNuclearModification2__3","",5, xAxis2);
   hNuclearModification2__3->SetBinContent(1,0.4279886);
   hNuclearModification2__3->SetBinContent(2,0.4447959);
   hNuclearModification2__3->SetBinContent(3,0.5041472);
   hNuclearModification2__3->SetBinContent(4,0.6796131);
   hNuclearModification2__3->SetBinContent(5,0.7424828);
   hNuclearModification2__3->SetBinError(1,0.05025093);
   hNuclearModification2__3->SetBinError(2,0.04192271);
   hNuclearModification2__3->SetBinError(3,0.04691724);
   hNuclearModification2__3->SetBinError(4,0.03941607);
   hNuclearModification2__3->SetBinError(5,0.07611912);
   hNuclearModification2__3->SetEntries(566.4226);
   hNuclearModification2__3->SetDirectory(0);
   hNuclearModification2__3->SetStats(0);
   hNuclearModification2__3->SetLineWidth(3);

   ci = TColor::GetColor("#006600");
   hNuclearModification2__3->SetMarkerColor(ci);
   hNuclearModification2__3->SetMarkerStyle(21);
   hNuclearModification2__3->SetMarkerSize(1.2);
   hNuclearModification2__3->GetXaxis()->SetTitle("D^{0} p_{T} (GeV/c)");
   hNuclearModification2__3->GetXaxis()->SetLabelFont(42);
   hNuclearModification2__3->GetXaxis()->SetLabelSize(0.035);
   hNuclearModification2__3->GetXaxis()->SetTitleSize(0.035);
   hNuclearModification2__3->GetXaxis()->SetTitleFont(42);
   hNuclearModification2__3->GetYaxis()->SetTitle("Uncorrected dN(D^{0})/dp_{T}");
   hNuclearModification2__3->GetYaxis()->SetLabelFont(42);
   hNuclearModification2__3->GetYaxis()->SetLabelSize(0.035);
   hNuclearModification2__3->GetYaxis()->SetTitleSize(0.035);
   hNuclearModification2__3->GetYaxis()->SetTitleOffset(0);
   hNuclearModification2__3->GetYaxis()->SetTitleFont(42);
   hNuclearModification2__3->GetZaxis()->SetLabelFont(42);
   hNuclearModification2__3->GetZaxis()->SetLabelSize(0.035);
   hNuclearModification2__3->GetZaxis()->SetTitleSize(0.035);
   hNuclearModification2__3->GetZaxis()->SetTitleFont(42);
   hNuclearModification2__3->Draw("psame");
   Double_t xAxis3[10] = {2, 3, 4, 5, 6, 8, 10, 12.5, 15, 20}; 
   
   TH1D *hNuclearModification3__4 = new TH1D("hNuclearModification3__4","",9, xAxis3);
   hNuclearModification3__4->SetBinContent(1,0.6315195);
   hNuclearModification3__4->SetBinContent(2,0.4432341);
   hNuclearModification3__4->SetBinContent(3,0.3674997);
   hNuclearModification3__4->SetBinContent(4,0.2758138);
   hNuclearModification3__4->SetBinContent(5,0.2623725);
   hNuclearModification3__4->SetBinContent(6,0.2722143);
   hNuclearModification3__4->SetBinContent(7,0.3291268);
   hNuclearModification3__4->SetBinContent(8,0.3207651);
   hNuclearModification3__4->SetBinContent(9,0.3937182);
   hNuclearModification3__4->SetBinError(1,0.1116572);
   hNuclearModification3__4->SetBinError(2,0.04535789);
   hNuclearModification3__4->SetBinError(3,0.03419174);
   hNuclearModification3__4->SetBinError(4,0.02501397);
   hNuclearModification3__4->SetBinError(5,0.01804332);
   hNuclearModification3__4->SetBinError(6,0.01714575);
   hNuclearModification3__4->SetBinError(7,0.02192858);
   hNuclearModification3__4->SetBinError(8,0.03017596);
   hNuclearModification3__4->SetBinError(9,0.04454357);
   hNuclearModification3__4->SetEntries(534.8555);
   hNuclearModification3__4->SetDirectory(0);
   hNuclearModification3__4->SetStats(0);
   hNuclearModification3__4->SetLineWidth(3);

   ci = TColor::GetColor("#006600");
   hNuclearModification3__4->SetMarkerColor(ci);
   hNuclearModification3__4->SetMarkerStyle(21);
   hNuclearModification3__4->SetMarkerSize(1.2);
   hNuclearModification3__4->GetXaxis()->SetTitle("D^{0} p_{T} (GeV/c)");
   hNuclearModification3__4->GetXaxis()->SetLabelFont(42);
   hNuclearModification3__4->GetXaxis()->SetLabelSize(0.035);
   hNuclearModification3__4->GetXaxis()->SetTitleSize(0.035);
   hNuclearModification3__4->GetXaxis()->SetTitleFont(42);
   hNuclearModification3__4->GetYaxis()->SetTitle("Uncorrected dN(D^{0})/dp_{T}");
   hNuclearModification3__4->GetYaxis()->SetLabelFont(42);
   hNuclearModification3__4->GetYaxis()->SetLabelSize(0.035);
   hNuclearModification3__4->GetYaxis()->SetTitleSize(0.035);
   hNuclearModification3__4->GetYaxis()->SetTitleOffset(0);
   hNuclearModification3__4->GetYaxis()->SetTitleFont(42);
   hNuclearModification3__4->GetZaxis()->SetLabelFont(42);
   hNuclearModification3__4->GetZaxis()->SetLabelSize(0.035);
   hNuclearModification3__4->GetZaxis()->SetTitleSize(0.035);
   hNuclearModification3__4->GetZaxis()->SetTitleFont(42);
   hNuclearModification3__4->Draw("psame");
   box = new TBox(0.62,0.9586116,0.752,1.045662);

   ci = 1186;
   color = new TColor(ci, 0.6, 1, 0.6, " ", 0.5);
   box->SetFillColor(ci);

   ci = TColor::GetColor("#99ff99");
   box->SetLineColor(ci);
   box->Draw();
   
   Double_t Graph0_fx1001[6] = {
   7.5,
   9,
   10.25,
   12,
   14.5,
   23};
   Double_t Graph0_fy1001[6] = {
   0.4664189,
   0.4856157,
   0.4065521,
   0.4536073,
   0.3694549,
   0.3521646};
   Double_t Graph0_fex1001[6] = {
   0.5,
   0.5,
   0.5,
   0.5,
   0.5,
   0.5};
   Double_t Graph0_fey1001[6] = {
   0.05427119,
   0.04438349,
   0.03578226,
   0.05860495,
   0.0518269,
   0.06848064};
   TGraphErrors *gre = new TGraphErrors(6,Graph0_fx1001,Graph0_fy1001,Graph0_fex1001,Graph0_fey1001);
   gre->SetName("Graph0");
   gre->SetTitle("Graph");

   ci = 1187;
   color = new TColor(ci, 0.4, 0.4, 0.4, " ", 0.5);
   gre->SetFillColor(ci);
   gre->SetFillStyle(1000);
   gre->SetMarkerStyle(20);
   
   TH1F *Graph_Graph_Graph11001 = new TH1F("Graph_Graph_Graph11001","Graph",100,5.35,25.15);
   Graph_Graph_Graph11001->SetMinimum(0.2590524);
   Graph_Graph_Graph11001->SetMaximum(0.5546307);
   Graph_Graph_Graph11001->SetDirectory(0);
   Graph_Graph_Graph11001->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph11001->SetLineColor(ci);
   Graph_Graph_Graph11001->SetLineStyle(0);
   Graph_Graph_Graph11001->SetMarkerStyle(20);
   Graph_Graph_Graph11001->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph11001->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph_Graph11001->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph_Graph11001->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph_Graph11001->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph_Graph11001->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph11001->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph11001->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph_Graph11001->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph_Graph11001->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph_Graph11001->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph_Graph11001->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph11001->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph11001->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph_Graph11001->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph_Graph11001->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph_Graph11001->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph_Graph11001);
   
   gre->Draw("2");
   
   Double_t Graph1_fx1002[6] = {
   7.5,
   9,
   10.25,
   12,
   14.5,
   23};
   Double_t Graph1_fy1002[6] = {
   0.4664189,
   0.4856157,
   0.4065521,
   0.4536073,
   0.3694549,
   0.3521646};
   Double_t Graph1_fex1002[6] = {
   1,
   0.5,
   0.75,
   1,
   1.5,
   7};
   Double_t Graph1_fey1002[6] = {
   0.02632702,
   0.0353226,
   0.0288783,
   0.03347202,
   0.03358354,
   0.0311077};
   gre = new TGraphErrors(6,Graph1_fx1002,Graph1_fy1002,Graph1_fex1002,Graph1_fey1002);
   gre->SetName("Graph1");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);
   gre->SetFillStyle(1000);

   ci = TColor::GetColor("#333333");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(29);
   gre->SetMarkerSize(2);
   
   TH1F *Graph_Graph_Graph21002 = new TH1F("Graph_Graph_Graph21002","Graph",100,4.15,32.35);
   Graph_Graph_Graph21002->SetMinimum(0.3010687);
   Graph_Graph_Graph21002->SetMaximum(0.5409264);
   Graph_Graph_Graph21002->SetDirectory(0);
   Graph_Graph_Graph21002->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph21002->SetLineColor(ci);
   Graph_Graph_Graph21002->SetLineStyle(0);
   Graph_Graph_Graph21002->SetMarkerStyle(20);
   Graph_Graph_Graph21002->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph21002->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph_Graph21002->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph_Graph21002->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph_Graph21002->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph_Graph21002->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph21002->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph21002->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph_Graph21002->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph_Graph21002->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph_Graph21002->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph_Graph21002->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph21002->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph21002->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph_Graph21002->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph_Graph21002->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph_Graph21002->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph_Graph21002);
   
   gre->Draw("p");
   
   Double_t Graph3_fx1003[3] = {
   3.75,
   5,
   6};
   Double_t Graph3_fy1003[3] = {
   0.7535804,
   0.6891029,
   0.5425812};
   Double_t Graph3_fex1003[3] = {
   0.4,
   0.4,
   0.4};
   Double_t Graph3_fey1003[3] = {
   0.21432,
   0.1218303,
   0.08921681};
   gre = new TGraphErrors(3,Graph3_fx1003,Graph3_fy1003,Graph3_fex1003,Graph3_fey1003);
   gre->SetName("Graph3");
   gre->SetTitle("Graph");

   ci = 1188;
   color = new TColor(ci, 0.4, 0.4, 0.4, " ", 0.5);
   gre->SetFillColor(ci);
   gre->SetFillStyle(1000);
   gre->SetMarkerStyle(20);
   
   TH1F *Graph_Graph_Graph41003 = new TH1F("Graph_Graph_Graph41003","Graph",100,3.045,6.705);
   Graph_Graph_Graph41003->SetMinimum(0.4019107);
   Graph_Graph_Graph41003->SetMaximum(1.019354);
   Graph_Graph_Graph41003->SetDirectory(0);
   Graph_Graph_Graph41003->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph41003->SetLineColor(ci);
   Graph_Graph_Graph41003->SetLineStyle(0);
   Graph_Graph_Graph41003->SetMarkerStyle(20);
   Graph_Graph_Graph41003->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph41003->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph_Graph41003->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph_Graph41003->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph_Graph41003->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph_Graph41003->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph41003->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph41003->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph_Graph41003->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph_Graph41003->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph_Graph41003->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph_Graph41003->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph41003->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph41003->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph_Graph41003->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph_Graph41003->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph_Graph41003->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph_Graph41003);
   
   gre->Draw("2");
   
   Double_t Graph4_fx1004[3] = {
   3.75,
   5,
   6};
   Double_t Graph4_fy1004[3] = {
   0.7535804,
   0.6891029,
   0.5425812};
   Double_t Graph4_fex1004[3] = {
   0.75,
   0.5,
   0.5};
   Double_t Graph4_fey1004[3] = {
   0.1209346,
   0.1059286,
   0.07363414};
   gre = new TGraphErrors(3,Graph4_fx1004,Graph4_fy1004,Graph4_fex1004,Graph4_fey1004);
   gre->SetName("Graph4");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);
   gre->SetFillStyle(1000);

   ci = TColor::GetColor("#333333");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(34);
   gre->SetMarkerSize(1.7);
   
   TH1F *Graph_Graph_Graph51004 = new TH1F("Graph_Graph_Graph51004","Graph",100,2.65,6.85);
   Graph_Graph_Graph51004->SetMinimum(0.4283902);
   Graph_Graph_Graph51004->SetMaximum(0.9150718);
   Graph_Graph_Graph51004->SetDirectory(0);
   Graph_Graph_Graph51004->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph51004->SetLineColor(ci);
   Graph_Graph_Graph51004->SetLineStyle(0);
   Graph_Graph_Graph51004->SetMarkerStyle(20);
   Graph_Graph_Graph51004->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph51004->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph_Graph51004->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph_Graph51004->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph_Graph51004->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph_Graph51004->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph51004->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph51004->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph_Graph51004->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph_Graph51004->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph_Graph51004->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph_Graph51004->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph51004->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph51004->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph_Graph51004->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph_Graph51004->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph_Graph51004->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph_Graph51004);
   
   gre->Draw("p");
   box = new TBox(0.92,0.9249388,1.1,1.075061);

   ci = 1191;
   color = new TColor(ci, 0.4, 0.4, 0.4, " ", 0.5);
   box->SetFillColor(ci);

   ci = TColor::GetColor("#666666");
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
   TLegendEntry *entry=leg->AddEntry("","charged hadrons |y| < 1.0","pf");

   ci = TColor::GetColor("#ffcc00");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1000);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry->SetTextSize(0.038);
   entry=leg->AddEntry("","D^{0} |y| < 1.0","pf");

   ci = 1185;
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
   entry=leg->AddEntry("gNuclearModification_Cor","Correlated syst. uncert.","f");

   ci = 1181;
   color = new TColor(ci, 1, 1, 1, " ", 0.01);
   entry->SetFillColor(ci);
   entry->SetFillStyle(1000);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(4);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry->SetTextSize(0.03);
   entry=leg->AddEntry("gNuclearModification_UnCor","Uncorrelated syst. uncert.","f");

   ci = 1180;
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
   entry=leg->AddEntry("","nonprompt J/#psi 1.6 < |y| < 2.4 (2.76 TeV)","pf");

   ci = 1189;
   color = new TColor(ci, 0.4, 0.4, 0.4, " ", 0.5);
   entry->SetFillColor(ci);
   entry->SetFillStyle(1000);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#333333");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(34);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);
   entry->SetTextSize(0.029);
   entry=leg->AddEntry("","nonprompt J/#psi |y| < 2.4 (2.76 TeV)","pf");

   ci = 1190;
   color = new TColor(ci, 0.4, 0.4, 0.4, " ", 0.5);
   entry->SetFillColor(ci);
   entry->SetFillStyle(1000);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#333333");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(29);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);
   entry->SetTextSize(0.029);
   leg->Draw();
   
   Double_t gNuclearModification_UnCor_fx3003[5] = {
   8.5,
   12.5,
   17.5,
   25,
   40};
   Double_t gNuclearModification_UnCor_fy3003[5] = {
   0.454751,
   0.4312376,
   0.4305237,
   0.5976104,
   0.3416682};
   Double_t gNuclearModification_UnCor_felx3003[5] = {
   1.071429,
   1.785714,
   1.785714,
   3.571429,
   7.142857};
   Double_t gNuclearModification_UnCor_fely3003[5] = {
   0.05328389,
   0.03183323,
   0.02774072,
   0.03354718,
   0.01839324};
   Double_t gNuclearModification_UnCor_fehx3003[5] = {
   1.071429,
   1.785714,
   1.785714,
   3.571429,
   7.142857};
   Double_t gNuclearModification_UnCor_fehy3003[5] = {
   0.05328389,
   0.03183323,
   0.02774072,
   0.03354718,
   0.01839324};
   grae = new TGraphAsymmErrors(5,gNuclearModification_UnCor_fx3003,gNuclearModification_UnCor_fy3003,gNuclearModification_UnCor_felx3003,gNuclearModification_UnCor_fehx3003,gNuclearModification_UnCor_fely3003,gNuclearModification_UnCor_fehy3003);
   grae->SetName("gNuclearModification_UnCor");
   grae->SetTitle("Graph");

   ci = 1180;
   color = new TColor(ci, 0, 0.6, 1, " ", 0.5);
   grae->SetFillColor(ci);
   grae->SetFillStyle(1000);
   grae->SetLineColor(0);
   grae->SetLineWidth(0);

   ci = TColor::GetColor("#0033cc");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(21);
   grae->SetMarkerSize(1.2);
   
   TH1F *Graph_gNuclearModification_UnCor3003 = new TH1F("Graph_gNuclearModification_UnCor3003","Graph",100,3.457143,51.11429);
   Graph_gNuclearModification_UnCor3003->SetMinimum(0.2924867);
   Graph_gNuclearModification_UnCor3003->SetMaximum(0.6619459);
   Graph_gNuclearModification_UnCor3003->SetDirectory(0);
   Graph_gNuclearModification_UnCor3003->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_gNuclearModification_UnCor3003->SetLineColor(ci);
   Graph_gNuclearModification_UnCor3003->SetMarkerStyle(20);
   Graph_gNuclearModification_UnCor3003->GetXaxis()->SetLabelFont(42);
   Graph_gNuclearModification_UnCor3003->GetXaxis()->SetLabelSize(0.035);
   Graph_gNuclearModification_UnCor3003->GetXaxis()->SetTitleSize(0.035);
   Graph_gNuclearModification_UnCor3003->GetXaxis()->SetTitleFont(42);
   Graph_gNuclearModification_UnCor3003->GetYaxis()->SetLabelFont(42);
   Graph_gNuclearModification_UnCor3003->GetYaxis()->SetLabelSize(0.035);
   Graph_gNuclearModification_UnCor3003->GetYaxis()->SetTitleSize(0.035);
   Graph_gNuclearModification_UnCor3003->GetYaxis()->SetTitleOffset(0);
   Graph_gNuclearModification_UnCor3003->GetYaxis()->SetTitleFont(42);
   Graph_gNuclearModification_UnCor3003->GetZaxis()->SetLabelFont(42);
   Graph_gNuclearModification_UnCor3003->GetZaxis()->SetLabelSize(0.035);
   Graph_gNuclearModification_UnCor3003->GetZaxis()->SetTitleSize(0.035);
   Graph_gNuclearModification_UnCor3003->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_gNuclearModification_UnCor3003);
   
   grae->Draw("5");
   
   Double_t gNuclearModification_Cor_fx3004[5] = {
   8.5,
   12.5,
   17.5,
   25,
   40};
   Double_t gNuclearModification_Cor_fy3004[5] = {
   0.454751,
   0.4312376,
   0.4305237,
   0.5976104,
   0.3416682};
   Double_t gNuclearModification_Cor_felx3004[5] = {
   1.5,
   2.5,
   2.5,
   5,
   10};
   Double_t gNuclearModification_Cor_fely3004[5] = {
   0.07066514,
   0.06701132,
   0.06690039,
   0.09286451,
   0.05309286};
   Double_t gNuclearModification_Cor_fehx3004[5] = {
   1.5,
   2.5,
   2.5,
   5,
   10};
   Double_t gNuclearModification_Cor_fehy3004[5] = {
   0.07066514,
   0.06701132,
   0.06690039,
   0.09286451,
   0.05309286};
   grae = new TGraphAsymmErrors(5,gNuclearModification_Cor_fx3004,gNuclearModification_Cor_fy3004,gNuclearModification_Cor_felx3004,gNuclearModification_Cor_fehx3004,gNuclearModification_Cor_fely3004,gNuclearModification_Cor_fehy3004);
   grae->SetName("gNuclearModification_Cor");
   grae->SetTitle("Graph");

   ci = 1181;
   color = new TColor(ci, 1, 1, 1, " ", 0.01);
   grae->SetFillColor(ci);
   grae->SetFillStyle(1000);
   grae->SetLineWidth(3);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_gNuclearModification_Cor3004 = new TH1F("Graph_gNuclearModification_Cor3004","Graph",100,2.7,54.3);
   Graph_gNuclearModification_Cor3004->SetMinimum(0.2483854);
   Graph_gNuclearModification_Cor3004->SetMaximum(0.7306649);
   Graph_gNuclearModification_Cor3004->SetDirectory(0);
   Graph_gNuclearModification_Cor3004->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_gNuclearModification_Cor3004->SetLineColor(ci);
   Graph_gNuclearModification_Cor3004->SetMarkerStyle(20);
   Graph_gNuclearModification_Cor3004->GetXaxis()->SetLabelFont(42);
   Graph_gNuclearModification_Cor3004->GetXaxis()->SetLabelSize(0.035);
   Graph_gNuclearModification_Cor3004->GetXaxis()->SetTitleSize(0.035);
   Graph_gNuclearModification_Cor3004->GetXaxis()->SetTitleFont(42);
   Graph_gNuclearModification_Cor3004->GetYaxis()->SetLabelFont(42);
   Graph_gNuclearModification_Cor3004->GetYaxis()->SetLabelSize(0.035);
   Graph_gNuclearModification_Cor3004->GetYaxis()->SetTitleSize(0.035);
   Graph_gNuclearModification_Cor3004->GetYaxis()->SetTitleOffset(0);
   Graph_gNuclearModification_Cor3004->GetYaxis()->SetTitleFont(42);
   Graph_gNuclearModification_Cor3004->GetZaxis()->SetLabelFont(42);
   Graph_gNuclearModification_Cor3004->GetZaxis()->SetLabelSize(0.035);
   Graph_gNuclearModification_Cor3004->GetZaxis()->SetTitleSize(0.035);
   Graph_gNuclearModification_Cor3004->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_gNuclearModification_Cor3004);
   
   grae->Draw("5");
   Double_t xAxis4[6] = {7, 10, 15, 20, 30, 50}; 
   
   TH1D *hNuclearModification__5 = new TH1D("hNuclearModification__5","",5, xAxis4);
   hNuclearModification__5->SetBinContent(1,0.454751);
   hNuclearModification__5->SetBinContent(2,0.4312376);
   hNuclearModification__5->SetBinContent(3,0.4305237);
   hNuclearModification__5->SetBinContent(4,0.5976104);
   hNuclearModification__5->SetBinContent(5,0.3416682);
   hNuclearModification__5->SetBinError(1,0.1442009);
   hNuclearModification__5->SetBinError(2,0.07163568);
   hNuclearModification__5->SetBinError(3,0.0733829);
   hNuclearModification__5->SetBinError(4,0.08926584);
   hNuclearModification__5->SetBinError(5,0.1097644);
   hNuclearModification__5->SetEntries(99.14019);

   ci = TColor::GetColor("#0033cc");
   hNuclearModification__5->SetLineColor(ci);
   hNuclearModification__5->SetLineWidth(3);

   ci = TColor::GetColor("#0033cc");
   hNuclearModification__5->SetMarkerColor(ci);
   hNuclearModification__5->SetMarkerStyle(21);
   hNuclearModification__5->SetMarkerSize(1.2);
   hNuclearModification__5->GetXaxis()->SetTitle("D^{0} p_{T} (GeV/c)");
   hNuclearModification__5->GetXaxis()->SetLabelFont(42);
   hNuclearModification__5->GetXaxis()->SetLabelSize(0.035);
   hNuclearModification__5->GetXaxis()->SetTitleSize(0.035);
   hNuclearModification__5->GetXaxis()->SetTitleFont(42);
   hNuclearModification__5->GetYaxis()->SetTitle("Uncorrected dN(D^{0})/dp_{T}");
   hNuclearModification__5->GetYaxis()->SetLabelFont(42);
   hNuclearModification__5->GetYaxis()->SetLabelSize(0.035);
   hNuclearModification__5->GetYaxis()->SetTitleSize(0.035);
   hNuclearModification__5->GetYaxis()->SetTitleFont(42);
   hNuclearModification__5->GetZaxis()->SetLabelFont(42);
   hNuclearModification__5->GetZaxis()->SetLabelSize(0.035);
   hNuclearModification__5->GetZaxis()->SetTitleSize(0.035);
   hNuclearModification__5->GetZaxis()->SetTitleFont(42);
   hNuclearModification__5->Draw("same");
   box = new TBox(0.92,0.9249388,1.1,1.075061);

   ci = 1191;
   color = new TColor(ci, 0.4, 0.4, 0.4, " ", 0.5);
   box->SetFillColor(ci);

   ci = TColor::GetColor("#666666");
   box->SetLineColor(ci);
   box->Draw();
   Double_t xAxis5[38] = {0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 3.2, 4, 4.8, 5.6, 6.4, 7.2, 9.6, 12, 14.4, 19.2, 24, 28.8, 35.2, 41.6, 48, 60.8, 73.6, 86.4, 103.6, 120.8, 140, 165, 250, 400}; 
   
   TH1D *RAA_0_100584__2__2__2__6 = new TH1D("RAA_0_100584__2__2__2__6","",37, xAxis5);
   RAA_0_100584__2__2__2__6->SetBinContent(1,-1);
   RAA_0_100584__2__2__2__6->SetBinContent(2,-1);
   RAA_0_100584__2__2__2__6->SetBinContent(3,0.3605047);
   RAA_0_100584__2__2__2__6->SetBinContent(4,0.3850582);
   RAA_0_100584__2__2__2__6->SetBinContent(5,0.4188398);
   RAA_0_100584__2__2__2__6->SetBinContent(6,0.4411318);
   RAA_0_100584__2__2__2__6->SetBinContent(7,0.4564881);
   RAA_0_100584__2__2__2__6->SetBinContent(8,0.4798668);
   RAA_0_100584__2__2__2__6->SetBinContent(9,0.4985427);
   RAA_0_100584__2__2__2__6->SetBinContent(10,0.516094);
   RAA_0_100584__2__2__2__6->SetBinContent(11,0.5160761);
   RAA_0_100584__2__2__2__6->SetBinContent(12,0.5157593);
   RAA_0_100584__2__2__2__6->SetBinContent(13,0.5090351);
   RAA_0_100584__2__2__2__6->SetBinContent(14,0.4718224);
   RAA_0_100584__2__2__2__6->SetBinContent(15,0.381288);
   RAA_0_100584__2__2__2__6->SetBinContent(16,0.3016595);
   RAA_0_100584__2__2__2__6->SetBinContent(17,0.2568378);
   RAA_0_100584__2__2__2__6->SetBinContent(18,0.2362394);
   RAA_0_100584__2__2__2__6->SetBinContent(19,0.2313929);
   RAA_0_100584__2__2__2__6->SetBinContent(20,0.2397292);
   RAA_0_100584__2__2__2__6->SetBinContent(21,0.2795414);
   RAA_0_100584__2__2__2__6->SetBinContent(22,0.3109757);
   RAA_0_100584__2__2__2__6->SetBinContent(23,0.359047);
   RAA_0_100584__2__2__2__6->SetBinContent(24,0.4103232);
   RAA_0_100584__2__2__2__6->SetBinContent(25,0.4535739);
   RAA_0_100584__2__2__2__6->SetBinContent(26,0.4871136);
   RAA_0_100584__2__2__2__6->SetBinContent(27,0.5390386);
   RAA_0_100584__2__2__2__6->SetBinContent(28,0.6064001);
   RAA_0_100584__2__2__2__6->SetBinContent(29,0.6583816);
   RAA_0_100584__2__2__2__6->SetBinContent(30,0.6914895);
   RAA_0_100584__2__2__2__6->SetBinContent(31,0.7510186);
   RAA_0_100584__2__2__2__6->SetBinContent(32,0.8027461);
   RAA_0_100584__2__2__2__6->SetBinContent(33,0.8064291);
   RAA_0_100584__2__2__2__6->SetBinContent(34,0.8144861);
   RAA_0_100584__2__2__2__6->SetBinContent(35,0.8835603);
   RAA_0_100584__2__2__2__6->SetBinContent(36,0.9077175);
   RAA_0_100584__2__2__2__6->SetBinContent(37,0.8767255);
   RAA_0_100584__2__2__2__6->SetBinError(1,1.279903e-05);
   RAA_0_100584__2__2__2__6->SetBinError(2,1.387861e-05);
   RAA_0_100584__2__2__2__6->SetBinError(3,1.439234e-05);
   RAA_0_100584__2__2__2__6->SetBinError(4,1.673993e-05);
   RAA_0_100584__2__2__2__6->SetBinError(5,1.996438e-05);
   RAA_0_100584__2__2__2__6->SetBinError(6,2.296277e-05);
   RAA_0_100584__2__2__2__6->SetBinError(7,2.652961e-05);
   RAA_0_100584__2__2__2__6->SetBinError(8,2.303589e-05);
   RAA_0_100584__2__2__2__6->SetBinError(9,2.956636e-05);
   RAA_0_100584__2__2__2__6->SetBinError(10,3.708897e-05);
   RAA_0_100584__2__2__2__6->SetBinError(11,4.504344e-05);
   RAA_0_100584__2__2__2__6->SetBinError(12,5.385318e-05);
   RAA_0_100584__2__2__2__6->SetBinError(13,6.332594e-05);
   RAA_0_100584__2__2__2__6->SetBinError(14,4.278537e-05);
   RAA_0_100584__2__2__2__6->SetBinError(15,6.266848e-05);
   RAA_0_100584__2__2__2__6->SetBinError(16,8.299302e-05);
   RAA_0_100584__2__2__2__6->SetBinError(17,0.0001099248);
   RAA_0_100584__2__2__2__6->SetBinError(18,0.0001477243);
   RAA_0_100584__2__2__2__6->SetBinError(19,0.0002002684);
   RAA_0_100584__2__2__2__6->SetBinError(20,0.0001962782);
   RAA_0_100584__2__2__2__6->SetBinError(21,0.00043599);
   RAA_0_100584__2__2__2__6->SetBinError(22,0.0007955915);
   RAA_0_100584__2__2__2__6->SetBinError(23,0.001121086);
   RAA_0_100584__2__2__2__6->SetBinError(24,0.002349183);
   RAA_0_100584__2__2__2__6->SetBinError(25,0.003934258);
   RAA_0_100584__2__2__2__6->SetBinError(26,0.005431548);
   RAA_0_100584__2__2__2__6->SetBinError(27,0.007593713);
   RAA_0_100584__2__2__2__6->SetBinError(28,0.008532134);
   RAA_0_100584__2__2__2__6->SetBinError(29,0.006423681);
   RAA_0_100584__2__2__2__6->SetBinError(30,0.009668077);
   RAA_0_100584__2__2__2__6->SetBinError(31,0.009013466);
   RAA_0_100584__2__2__2__6->SetBinError(32,0.01230905);
   RAA_0_100584__2__2__2__6->SetBinError(33,0.01985833);
   RAA_0_100584__2__2__2__6->SetBinError(34,0.02497965);
   RAA_0_100584__2__2__2__6->SetBinError(35,0.04145952);
   RAA_0_100584__2__2__2__6->SetBinError(36,0.0459072);
   RAA_0_100584__2__2__2__6->SetBinError(37,0.152905);
   RAA_0_100584__2__2__2__6->SetMinimum(0);
   RAA_0_100584__2__2__2__6->SetMaximum(1.8);
   RAA_0_100584__2__2__2__6->SetEntries(11078.9);
   RAA_0_100584__2__2__2__6->SetDirectory(0);
   RAA_0_100584__2__2__2__6->SetStats(0);

   ci = TColor::GetColor("#ffcc00");
   RAA_0_100584__2__2__2__6->SetFillColor(ci);

   ci = TColor::GetColor("#000099");
   RAA_0_100584__2__2__2__6->SetLineColor(ci);
   RAA_0_100584__2__2__2__6->SetLineWidth(2);
   RAA_0_100584__2__2__2__6->SetMarkerStyle(20);
   RAA_0_100584__2__2__2__6->SetMarkerSize(1.3);
   RAA_0_100584__2__2__2__6->GetXaxis()->SetTitle("p_{T} (GeV)");
   RAA_0_100584__2__2__2__6->GetXaxis()->SetRange(1,37);
   RAA_0_100584__2__2__2__6->GetXaxis()->CenterTitle(true);
   RAA_0_100584__2__2__2__6->GetXaxis()->SetLabelFont(42);
   RAA_0_100584__2__2__2__6->GetXaxis()->SetLabelOffset(-0.002);
   RAA_0_100584__2__2__2__6->GetXaxis()->SetLabelSize(0.045);
   RAA_0_100584__2__2__2__6->GetXaxis()->SetTitleSize(0.05);
   RAA_0_100584__2__2__2__6->GetXaxis()->SetTitleFont(42);
   RAA_0_100584__2__2__2__6->GetYaxis()->SetTitle("R_{AA}, R_{pA}");
   RAA_0_100584__2__2__2__6->GetYaxis()->CenterTitle(true);
   RAA_0_100584__2__2__2__6->GetYaxis()->SetNdivisions(508);
   RAA_0_100584__2__2__2__6->GetYaxis()->SetLabelFont(42);
   RAA_0_100584__2__2__2__6->GetYaxis()->SetLabelSize(0.045);
   RAA_0_100584__2__2__2__6->GetYaxis()->SetTitleSize(0.05);
   RAA_0_100584__2__2__2__6->GetYaxis()->SetTitleOffset(0);
   RAA_0_100584__2__2__2__6->GetYaxis()->SetTitleFont(42);
   RAA_0_100584__2__2__2__6->GetZaxis()->SetLabelFont(42);
   RAA_0_100584__2__2__2__6->GetZaxis()->SetLabelSize(0.05);
   RAA_0_100584__2__2__2__6->GetZaxis()->SetTitleSize(0.09);
   RAA_0_100584__2__2__2__6->GetZaxis()->SetTitleFont(42);
   RAA_0_100584__2__2__2__6->Draw("same");
   Double_t xAxis6[6] = {20, 25, 30, 40, 60, 100}; 
   
   TH1D *hNuclearModification2__7 = new TH1D("hNuclearModification2__7","",5, xAxis6);
   hNuclearModification2__7->SetBinContent(1,0.4279886);
   hNuclearModification2__7->SetBinContent(2,0.4447959);
   hNuclearModification2__7->SetBinContent(3,0.5041472);
   hNuclearModification2__7->SetBinContent(4,0.6796131);
   hNuclearModification2__7->SetBinContent(5,0.7424828);
   hNuclearModification2__7->SetBinError(1,0.05025093);
   hNuclearModification2__7->SetBinError(2,0.04192271);
   hNuclearModification2__7->SetBinError(3,0.04691724);
   hNuclearModification2__7->SetBinError(4,0.03941607);
   hNuclearModification2__7->SetBinError(5,0.07611912);
   hNuclearModification2__7->SetEntries(566.4226);
   hNuclearModification2__7->SetStats(0);
   hNuclearModification2__7->SetLineWidth(3);

   ci = TColor::GetColor("#006600");
   hNuclearModification2__7->SetMarkerColor(ci);
   hNuclearModification2__7->SetMarkerStyle(21);
   hNuclearModification2__7->SetMarkerSize(1.2);
   hNuclearModification2__7->GetXaxis()->SetTitle("D^{0} p_{T} (GeV/c)");
   hNuclearModification2__7->GetXaxis()->SetLabelFont(42);
   hNuclearModification2__7->GetXaxis()->SetLabelSize(0.035);
   hNuclearModification2__7->GetXaxis()->SetTitleSize(0.035);
   hNuclearModification2__7->GetXaxis()->SetTitleFont(42);
   hNuclearModification2__7->GetYaxis()->SetTitle("Uncorrected dN(D^{0})/dp_{T}");
   hNuclearModification2__7->GetYaxis()->SetLabelFont(42);
   hNuclearModification2__7->GetYaxis()->SetLabelSize(0.035);
   hNuclearModification2__7->GetYaxis()->SetTitleSize(0.035);
   hNuclearModification2__7->GetYaxis()->SetTitleOffset(0);
   hNuclearModification2__7->GetYaxis()->SetTitleFont(42);
   hNuclearModification2__7->GetZaxis()->SetLabelFont(42);
   hNuclearModification2__7->GetZaxis()->SetLabelSize(0.035);
   hNuclearModification2__7->GetZaxis()->SetTitleSize(0.035);
   hNuclearModification2__7->GetZaxis()->SetTitleFont(42);
   hNuclearModification2__7->Draw("psame");
   Double_t xAxis7[10] = {2, 3, 4, 5, 6, 8, 10, 12.5, 15, 20}; 
   
   TH1D *hNuclearModification3__8 = new TH1D("hNuclearModification3__8","",9, xAxis7);
   hNuclearModification3__8->SetBinContent(1,0.6315195);
   hNuclearModification3__8->SetBinContent(2,0.4432341);
   hNuclearModification3__8->SetBinContent(3,0.3674997);
   hNuclearModification3__8->SetBinContent(4,0.2758138);
   hNuclearModification3__8->SetBinContent(5,0.2623725);
   hNuclearModification3__8->SetBinContent(6,0.2722143);
   hNuclearModification3__8->SetBinContent(7,0.3291268);
   hNuclearModification3__8->SetBinContent(8,0.3207651);
   hNuclearModification3__8->SetBinContent(9,0.3937182);
   hNuclearModification3__8->SetBinError(1,0.1116572);
   hNuclearModification3__8->SetBinError(2,0.04535789);
   hNuclearModification3__8->SetBinError(3,0.03419174);
   hNuclearModification3__8->SetBinError(4,0.02501397);
   hNuclearModification3__8->SetBinError(5,0.01804332);
   hNuclearModification3__8->SetBinError(6,0.01714575);
   hNuclearModification3__8->SetBinError(7,0.02192858);
   hNuclearModification3__8->SetBinError(8,0.03017596);
   hNuclearModification3__8->SetBinError(9,0.04454357);
   hNuclearModification3__8->SetEntries(534.8555);
   hNuclearModification3__8->SetStats(0);
   hNuclearModification3__8->SetLineWidth(3);

   ci = TColor::GetColor("#006600");
   hNuclearModification3__8->SetMarkerColor(ci);
   hNuclearModification3__8->SetMarkerStyle(21);
   hNuclearModification3__8->SetMarkerSize(1.2);
   hNuclearModification3__8->GetXaxis()->SetTitle("D^{0} p_{T} (GeV/c)");
   hNuclearModification3__8->GetXaxis()->SetLabelFont(42);
   hNuclearModification3__8->GetXaxis()->SetLabelSize(0.035);
   hNuclearModification3__8->GetXaxis()->SetTitleSize(0.035);
   hNuclearModification3__8->GetXaxis()->SetTitleFont(42);
   hNuclearModification3__8->GetYaxis()->SetTitle("Uncorrected dN(D^{0})/dp_{T}");
   hNuclearModification3__8->GetYaxis()->SetLabelFont(42);
   hNuclearModification3__8->GetYaxis()->SetLabelSize(0.035);
   hNuclearModification3__8->GetYaxis()->SetTitleSize(0.035);
   hNuclearModification3__8->GetYaxis()->SetTitleOffset(0);
   hNuclearModification3__8->GetYaxis()->SetTitleFont(42);
   hNuclearModification3__8->GetZaxis()->SetLabelFont(42);
   hNuclearModification3__8->GetZaxis()->SetLabelSize(0.035);
   hNuclearModification3__8->GetZaxis()->SetTitleSize(0.035);
   hNuclearModification3__8->GetZaxis()->SetTitleFont(42);
   hNuclearModification3__8->Draw("psame");
   
   Double_t Graph0_fx1005[6] = {
   7.5,
   9,
   10.25,
   12,
   14.5,
   23};
   Double_t Graph0_fy1005[6] = {
   0.4664189,
   0.4856157,
   0.4065521,
   0.4536073,
   0.3694549,
   0.3521646};
   Double_t Graph0_fex1005[6] = {
   1,
   0.5,
   0.75,
   1,
   1.5,
   7};
   Double_t Graph0_fey1005[6] = {
   0.02632702,
   0.0353226,
   0.0288783,
   0.03347202,
   0.03358354,
   0.0311077};
   gre = new TGraphErrors(6,Graph0_fx1005,Graph0_fy1005,Graph0_fex1005,Graph0_fey1005);
   gre->SetName("Graph0");
   gre->SetTitle("Graph");
   gre->SetFillStyle(1000);

   ci = TColor::GetColor("#333333");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(29);
   gre->SetMarkerSize(2);
   
   TH1F *Graph_Graph01005 = new TH1F("Graph_Graph01005","Graph",100,4.15,32.35);
   Graph_Graph01005->SetMinimum(0.3010687);
   Graph_Graph01005->SetMaximum(0.5409264);
   Graph_Graph01005->SetDirectory(0);
   Graph_Graph01005->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph01005->SetLineColor(ci);
   Graph_Graph01005->SetLineStyle(0);
   Graph_Graph01005->SetMarkerStyle(20);
   Graph_Graph01005->GetXaxis()->SetLabelFont(42);
   Graph_Graph01005->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph01005->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph01005->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph01005->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph01005->GetXaxis()->SetTitleFont(42);
   Graph_Graph01005->GetYaxis()->SetLabelFont(42);
   Graph_Graph01005->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph01005->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph01005->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph01005->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph01005->GetYaxis()->SetTitleFont(42);
   Graph_Graph01005->GetZaxis()->SetLabelFont(42);
   Graph_Graph01005->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph01005->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph01005->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph01005->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph01005);
   
   gre->Draw("p");
   
   Double_t Graph1_fx1006[3] = {
   3.75,
   5,
   6};
   Double_t Graph1_fy1006[3] = {
   0.7535804,
   0.6891029,
   0.5425812};
   Double_t Graph1_fex1006[3] = {
   0.75,
   0.5,
   0.5};
   Double_t Graph1_fey1006[3] = {
   0.1209346,
   0.1059286,
   0.07363414};
   gre = new TGraphErrors(3,Graph1_fx1006,Graph1_fy1006,Graph1_fex1006,Graph1_fey1006);
   gre->SetName("Graph1");
   gre->SetTitle("Graph");
   gre->SetFillStyle(1000);

   ci = TColor::GetColor("#333333");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(34);
   gre->SetMarkerSize(1.7);
   
   TH1F *Graph_Graph11006 = new TH1F("Graph_Graph11006","Graph",100,2.65,6.85);
   Graph_Graph11006->SetMinimum(0.4283902);
   Graph_Graph11006->SetMaximum(0.9150718);
   Graph_Graph11006->SetDirectory(0);
   Graph_Graph11006->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph11006->SetLineColor(ci);
   Graph_Graph11006->SetLineStyle(0);
   Graph_Graph11006->SetMarkerStyle(20);
   Graph_Graph11006->GetXaxis()->SetLabelFont(42);
   Graph_Graph11006->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph11006->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph11006->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph11006->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph11006->GetXaxis()->SetTitleFont(42);
   Graph_Graph11006->GetYaxis()->SetLabelFont(42);
   Graph_Graph11006->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph11006->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph11006->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph11006->GetYaxis()->SetTitleOffset(1.25);
   Graph_Graph11006->GetYaxis()->SetTitleFont(42);
   Graph_Graph11006->GetZaxis()->SetLabelFont(42);
   Graph_Graph11006->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph11006->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph11006->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph11006->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph11006);
   
   gre->Draw("p");
   Double_t xAxis8[6] = {7, 10, 15, 20, 30, 50}; 
   
   TH1D *hNuclearModification__9 = new TH1D("hNuclearModification__9","",5, xAxis8);
   hNuclearModification__9->SetBinContent(1,0.454751);
   hNuclearModification__9->SetBinContent(2,0.4312376);
   hNuclearModification__9->SetBinContent(3,0.4305237);
   hNuclearModification__9->SetBinContent(4,0.5976104);
   hNuclearModification__9->SetBinContent(5,0.3416682);
   hNuclearModification__9->SetBinError(1,0.1442009);
   hNuclearModification__9->SetBinError(2,0.07163568);
   hNuclearModification__9->SetBinError(3,0.0733829);
   hNuclearModification__9->SetBinError(4,0.08926584);
   hNuclearModification__9->SetBinError(5,0.1097644);
   hNuclearModification__9->SetEntries(99.14019);

   ci = TColor::GetColor("#0033cc");
   hNuclearModification__9->SetLineColor(ci);
   hNuclearModification__9->SetLineWidth(3);

   ci = TColor::GetColor("#0033cc");
   hNuclearModification__9->SetMarkerColor(ci);
   hNuclearModification__9->SetMarkerStyle(21);
   hNuclearModification__9->SetMarkerSize(1.2);
   hNuclearModification__9->GetXaxis()->SetTitle("D^{0} p_{T} (GeV/c)");
   hNuclearModification__9->GetXaxis()->SetLabelFont(42);
   hNuclearModification__9->GetXaxis()->SetLabelSize(0.035);
   hNuclearModification__9->GetXaxis()->SetTitleSize(0.035);
   hNuclearModification__9->GetXaxis()->SetTitleFont(42);
   hNuclearModification__9->GetYaxis()->SetTitle("Uncorrected dN(D^{0})/dp_{T}");
   hNuclearModification__9->GetYaxis()->SetLabelFont(42);
   hNuclearModification__9->GetYaxis()->SetLabelSize(0.035);
   hNuclearModification__9->GetYaxis()->SetTitleSize(0.035);
   hNuclearModification__9->GetYaxis()->SetTitleFont(42);
   hNuclearModification__9->GetZaxis()->SetLabelFont(42);
   hNuclearModification__9->GetZaxis()->SetLabelSize(0.035);
   hNuclearModification__9->GetZaxis()->SetTitleSize(0.035);
   hNuclearModification__9->GetZaxis()->SetTitleFont(42);
   hNuclearModification__9->Draw("same p");
      tex = new TLatex(0.81,0.21,"B^{#pm}");
tex->SetNDC();
   tex->SetTextSize(0.08);
   tex->SetLineWidth(2);
   tex->Draw();
   
   TH2F *hemptyEff_copy__10 = new TH2F("hemptyEff_copy__10","",50,0.5,600,10,0,1.55);
   hemptyEff_copy__10->SetMinimum(0);
   hemptyEff_copy__10->SetMaximum(2);
   hemptyEff_copy__10->SetDirectory(0);
   hemptyEff_copy__10->SetStats(0);

   ci = TColor::GetColor("#000099");
   hemptyEff_copy__10->SetLineColor(ci);
   hemptyEff_copy__10->SetMarkerStyle(20);
   hemptyEff_copy__10->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   hemptyEff_copy__10->GetXaxis()->CenterTitle(true);
   hemptyEff_copy__10->GetXaxis()->SetLabelFont(42);
   hemptyEff_copy__10->GetXaxis()->SetTitleSize(0.05);
   hemptyEff_copy__10->GetXaxis()->SetTitleFont(42);
   hemptyEff_copy__10->GetYaxis()->SetTitle("R_{AA}");
   hemptyEff_copy__10->GetYaxis()->CenterTitle(true);
   hemptyEff_copy__10->GetYaxis()->SetLabelFont(42);
   hemptyEff_copy__10->GetYaxis()->SetTitleSize(0.05);
   hemptyEff_copy__10->GetYaxis()->SetTitleOffset(1.1);
   hemptyEff_copy__10->GetYaxis()->SetTitleFont(42);
   hemptyEff_copy__10->GetZaxis()->SetLabelFont(42);
   hemptyEff_copy__10->GetZaxis()->SetLabelSize(0.035);
   hemptyEff_copy__10->GetZaxis()->SetTitleSize(0.035);
   hemptyEff_copy__10->GetZaxis()->SetTitleFont(42);
   hemptyEff_copy__10->Draw("sameaxis");
   canvasRAA->Modified();
   canvasRAA->cd();
   canvasRAA->SetSelected(canvasRAA);
}
