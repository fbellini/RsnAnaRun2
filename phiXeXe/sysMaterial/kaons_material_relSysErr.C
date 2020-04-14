{
//=========Macro generated from canvas: Canvas_1/Canvas_1
//=========  (Tue Nov 19 16:57:56 2013) by ROOT version5.34/11
   TCanvas *Canvas_1 = new TCanvas("Canvas_1", "Canvas_1",258,73,705,497);
   Canvas_1->Range(-0.005135513,-0.2507614,1.556919,0.2497462);
   Canvas_1->SetFillColor(0);
   Canvas_1->SetBorderMode(0);
   Canvas_1->SetBorderSize(2);
   Canvas_1->SetFrameBorderMode(0);
   Canvas_1->SetFrameBorderMode(0);
   Double_t xAxis2[36] = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3}; 
   
   TH1D *hKaonMinusMinBiasSysErrRel = new TH1D("hKaonMinusMinBiasSysErrRel","hKaonMinusMinBiasSysErrRel",35, xAxis2);
   hKaonMinusMinBiasSysErrRel->SetBinContent(6,0.1106438);
   hKaonMinusMinBiasSysErrRel->SetBinContent(7,0.02137536);
   hKaonMinusMinBiasSysErrRel->SetBinContent(8,0.00753827);
   hKaonMinusMinBiasSysErrRel->SetBinContent(9,0.005393446);
   hKaonMinusMinBiasSysErrRel->SetBinContent(10,0.005060986);
   hKaonMinusMinBiasSysErrRel->SetBinContent(11,0.005009453);
   hKaonMinusMinBiasSysErrRel->SetBinContent(12,0.005001465);
   hKaonMinusMinBiasSysErrRel->SetBinContent(13,0.005000227);
   hKaonMinusMinBiasSysErrRel->SetBinContent(14,0.005000035);
   hKaonMinusMinBiasSysErrRel->SetBinContent(15,0.005000005);
   hKaonMinusMinBiasSysErrRel->SetBinContent(16,0.005000001);
   hKaonMinusMinBiasSysErrRel->SetBinContent(17,0.005);
   hKaonMinusMinBiasSysErrRel->SetBinContent(18,0.005);
   hKaonMinusMinBiasSysErrRel->SetBinContent(19,0.005);
   hKaonMinusMinBiasSysErrRel->SetBinContent(20,0.005);
   hKaonMinusMinBiasSysErrRel->SetBinContent(21,0.005);
   hKaonMinusMinBiasSysErrRel->SetBinContent(22,0.005);
   hKaonMinusMinBiasSysErrRel->SetBinContent(23,0.005);
   hKaonMinusMinBiasSysErrRel->SetMinimum(-0.2);
   hKaonMinusMinBiasSysErrRel->SetMaximum(0.2);
   hKaonMinusMinBiasSysErrRel->SetEntries(18);
   hKaonMinusMinBiasSysErrRel->SetStats(0);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#000099");
   hKaonMinusMinBiasSysErrRel->SetLineColor(ci);
   hKaonMinusMinBiasSysErrRel->SetMarkerColor(2);
   hKaonMinusMinBiasSysErrRel->SetMarkerStyle(28);
   hKaonMinusMinBiasSysErrRel->GetXaxis()->SetTitle("p_T");
   hKaonMinusMinBiasSysErrRel->GetXaxis()->SetRange(4,24);
   hKaonMinusMinBiasSysErrRel->GetXaxis()->SetLabelFont(42);
   hKaonMinusMinBiasSysErrRel->GetXaxis()->SetLabelSize(0.035);
   hKaonMinusMinBiasSysErrRel->GetXaxis()->SetTitleSize(0.035);
   hKaonMinusMinBiasSysErrRel->GetXaxis()->SetTitleFont(42);
   hKaonMinusMinBiasSysErrRel->GetYaxis()->SetTitle("rel. sys. err.");
   hKaonMinusMinBiasSysErrRel->GetYaxis()->SetLabelFont(42);
   hKaonMinusMinBiasSysErrRel->GetYaxis()->SetLabelSize(0.035);
   hKaonMinusMinBiasSysErrRel->GetYaxis()->SetTitleSize(0.035);
   hKaonMinusMinBiasSysErrRel->GetYaxis()->SetTitleFont(42);
   hKaonMinusMinBiasSysErrRel->GetZaxis()->SetLabelFont(42);
   hKaonMinusMinBiasSysErrRel->GetZaxis()->SetLabelSize(0.035);
   hKaonMinusMinBiasSysErrRel->GetZaxis()->SetTitleSize(0.035);
   hKaonMinusMinBiasSysErrRel->GetZaxis()->SetTitleFont(42);
   hKaonMinusMinBiasSysErrRel->Draw("");
   
   TPaveText *pt = new TPaveText(0.2303709,0.9339002,0.7696291,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   TText *text = pt->AddText("hKaonMinusMinBiasSysErrRel");
   pt->Draw();
   Canvas_1->Modified();
   Canvas_1->cd();
   Canvas_1->SetSelected(Canvas_1);
}
