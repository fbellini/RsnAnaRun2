{
//=========Macro generated from canvas: Canvas_1/Canvas_1
//=========  (Tue Nov 19 16:28:18 2013) by ROOT version5.34/11
   TCanvas *Canvas_1 = new TCanvas("Canvas_1", "Canvas_1",258,73,705,497);
   Canvas_1->Range(-0.08616262,-0.2507614,1.788302,0.2497462);
   Canvas_1->SetFillColor(0);
   Canvas_1->SetBorderMode(0);
   Canvas_1->SetBorderSize(2);
   Canvas_1->SetFrameBorderMode(0);
   Canvas_1->SetFrameBorderMode(0);
   Double_t xAxis1[36] = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3}; 
   
   TH1D *hPionMinusMinBiasSysErrRel = new TH1D("hPionMinusMinBiasSysErrRel","hPionMinusMinBiasSysErrRel",35, xAxis1);
   hPionMinusMinBiasSysErrRel->SetBinContent(5,0.02593945);
   hPionMinusMinBiasSysErrRel->SetBinContent(6,0.01247019);
   hPionMinusMinBiasSysErrRel->SetBinContent(7,0.007665002);
   hPionMinusMinBiasSysErrRel->SetBinContent(8,0.005950744);
   hPionMinusMinBiasSysErrRel->SetBinContent(9,0.00533918);
   hPionMinusMinBiasSysErrRel->SetBinContent(10,0.005121003);
   hPionMinusMinBiasSysErrRel->SetBinContent(11,0.005043168);
   hPionMinusMinBiasSysErrRel->SetBinContent(12,0.0050154);
   hPionMinusMinBiasSysErrRel->SetBinContent(13,0.005005494);
   hPionMinusMinBiasSysErrRel->SetBinContent(14,0.00500196);
   hPionMinusMinBiasSysErrRel->SetBinContent(15,0.005000699);
   hPionMinusMinBiasSysErrRel->SetBinContent(16,0.005000249);
   hPionMinusMinBiasSysErrRel->SetBinContent(17,0.005000089);
   hPionMinusMinBiasSysErrRel->SetBinContent(18,0.005000032);
   hPionMinusMinBiasSysErrRel->SetBinContent(19,0.005000011);
   hPionMinusMinBiasSysErrRel->SetBinContent(20,0.005000004);
   hPionMinusMinBiasSysErrRel->SetBinContent(21,0.005000001);
   hPionMinusMinBiasSysErrRel->SetBinContent(22,0.005);
   hPionMinusMinBiasSysErrRel->SetBinContent(23,0.005);
   hPionMinusMinBiasSysErrRel->SetBinContent(24,0.005);
   hPionMinusMinBiasSysErrRel->SetBinContent(25,0.005);
   hPionMinusMinBiasSysErrRel->SetMinimum(-0.2);
   hPionMinusMinBiasSysErrRel->SetMaximum(0.2);
   hPionMinusMinBiasSysErrRel->SetEntries(21);
   
   TPaveStats *ptstats = new TPaveStats(0.78,0.775,0.98,0.935,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetFillColor(0);
   ptstats->SetTextAlign(12);
   ptstats->SetTextFont(42);
   TText *text = ptstats->AddText("hPionMinusMinBiasSysErrRel");
   text->SetTextSize(0.0368);
   text = ptstats->AddText("Entries = 21     ");
   text = ptstats->AddText("Mean  =  0.636");
   text = ptstats->AddText("RMS   = 0.3729");
   ptstats->SetOptStat(1111);
   ptstats->SetOptFit(0);
   ptstats->Draw();
   hPionMinusMinBiasSysErrRel->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(hPionMinusMinBiasSysErrRel);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#000099");
   hPionMinusMinBiasSysErrRel->SetLineColor(ci);
   hPionMinusMinBiasSysErrRel->SetMarkerColor(2);
   hPionMinusMinBiasSysErrRel->SetMarkerStyle(28);
   hPionMinusMinBiasSysErrRel->GetXaxis()->SetTitle("p_T");
   hPionMinusMinBiasSysErrRel->GetXaxis()->SetRange(3,26);
   hPionMinusMinBiasSysErrRel->GetXaxis()->SetLabelFont(42);
   hPionMinusMinBiasSysErrRel->GetXaxis()->SetLabelSize(0.035);
   hPionMinusMinBiasSysErrRel->GetXaxis()->SetTitleSize(0.035);
   hPionMinusMinBiasSysErrRel->GetXaxis()->SetTitleFont(42);
   hPionMinusMinBiasSysErrRel->GetYaxis()->SetTitle("rel. sys. err.");
   hPionMinusMinBiasSysErrRel->GetYaxis()->SetLabelFont(42);
   hPionMinusMinBiasSysErrRel->GetYaxis()->SetLabelSize(0.035);
   hPionMinusMinBiasSysErrRel->GetYaxis()->SetTitleSize(0.035);
   hPionMinusMinBiasSysErrRel->GetYaxis()->SetTitleFont(42);
   hPionMinusMinBiasSysErrRel->GetZaxis()->SetLabelFont(42);
   hPionMinusMinBiasSysErrRel->GetZaxis()->SetLabelSize(0.035);
   hPionMinusMinBiasSysErrRel->GetZaxis()->SetTitleSize(0.035);
   hPionMinusMinBiasSysErrRel->GetZaxis()->SetTitleFont(42);
   hPionMinusMinBiasSysErrRel->Draw("");
   
   TPaveText *pt = new TPaveText(0.236077,0.9339002,0.763923,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   text = pt->AddText("hPionMinusMinBiasSysErrRel");
   pt->Draw();
   Canvas_1->Modified();
   Canvas_1->cd();
   Canvas_1->SetSelected(Canvas_1);
}
