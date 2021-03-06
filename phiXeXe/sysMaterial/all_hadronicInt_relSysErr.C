{
//=========Macro generated from canvas: Canvas_1/Canvas_1
//=========  (Tue Nov 19 17:02:53 2013) by ROOT version5.34/11
   TCanvas *Canvas_1 = new TCanvas("Canvas_1", "Canvas_1",258,73,815,537);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   Canvas_1->Range(-0.1459893,-0.007292754,1.728342,0.05508658);
   Canvas_1->SetFillColor(0);
   Canvas_1->SetBorderMode(0);
   Canvas_1->SetBorderSize(2);
   Canvas_1->SetLeftMargin(0.1312411);
   Canvas_1->SetRightMargin(0.06847361);
   Canvas_1->SetTopMargin(0.06313646);
   Canvas_1->SetBottomMargin(0.1364562);
   Canvas_1->SetFrameBorderMode(0);
   Canvas_1->SetFrameBorderMode(0);
   Double_t xAxis3[36] = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3}; 
   
   TH1D *hPionMinusMinBiasSysErrRel = new TH1D("hPionMinusMinBiasSysErrRel","hPionMinusMinBiasSysErrRel",35, xAxis3);
   hPionMinusMinBiasSysErrRel->SetBinContent(5,0.01708287);
   hPionMinusMinBiasSysErrRel->SetBinContent(6,0.01471419);
   hPionMinusMinBiasSysErrRel->SetBinContent(7,0.01354882);
   hPionMinusMinBiasSysErrRel->SetBinContent(8,0.01241823);
   hPionMinusMinBiasSysErrRel->SetBinContent(9,0.01150448);
   hPionMinusMinBiasSysErrRel->SetBinContent(10,0.01081309);
   hPionMinusMinBiasSysErrRel->SetBinContent(11,0.01066405);
   hPionMinusMinBiasSysErrRel->SetBinContent(12,0.01093415);
   hPionMinusMinBiasSysErrRel->SetBinContent(13,0.01125151);
   hPionMinusMinBiasSysErrRel->SetBinContent(14,0.01143714);
   hPionMinusMinBiasSysErrRel->SetBinContent(15,0.01160452);
   hPionMinusMinBiasSysErrRel->SetBinContent(16,0.01169667);
   hPionMinusMinBiasSysErrRel->SetBinContent(17,0.01203264);
   hPionMinusMinBiasSysErrRel->SetBinContent(18,0.01214947);
   hPionMinusMinBiasSysErrRel->SetBinContent(19,0.012206);
   hPionMinusMinBiasSysErrRel->SetBinContent(20,0.012206);
   hPionMinusMinBiasSysErrRel->SetBinContent(21,0.012206);
   hPionMinusMinBiasSysErrRel->SetBinContent(22,0.012206);
   hPionMinusMinBiasSysErrRel->SetBinContent(23,0.012206);
   hPionMinusMinBiasSysErrRel->SetBinContent(24,0.012206);
   hPionMinusMinBiasSysErrRel->SetBinContent(25,0.012206);
   hPionMinusMinBiasSysErrRel->SetMinimum(0.001219293);
   hPionMinusMinBiasSysErrRel->SetMaximum(0.05114817);
   hPionMinusMinBiasSysErrRel->SetEntries(21);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#000099");
   hPionMinusMinBiasSysErrRel->SetLineColor(ci);
   hPionMinusMinBiasSysErrRel->SetLineWidth(2);
   hPionMinusMinBiasSysErrRel->SetMarkerColor(4);
   hPionMinusMinBiasSysErrRel->GetXaxis()->SetTitle("p_T");
   hPionMinusMinBiasSysErrRel->GetXaxis()->SetRange(3,26);
   hPionMinusMinBiasSysErrRel->GetXaxis()->SetLabelFont(42);
   hPionMinusMinBiasSysErrRel->GetXaxis()->SetLabelSize(0.035);
   hPionMinusMinBiasSysErrRel->GetXaxis()->SetTitleSize(0.035);
   hPionMinusMinBiasSysErrRel->GetXaxis()->SetTitleFont(42);
   hPionMinusMinBiasSysErrRel->GetYaxis()->SetTitle("rel. sys. err.");
   hPionMinusMinBiasSysErrRel->GetYaxis()->SetLabelFont(42);
   hPionMinusMinBiasSysErrRel->GetYaxis()->SetLabelSize(0.035);
   hPionMinusMinBiasSysErrRel->GetYaxis()->SetTitleSize(0.05);
   hPionMinusMinBiasSysErrRel->GetYaxis()->SetTitleFont(42);
   hPionMinusMinBiasSysErrRel->GetZaxis()->SetLabelFont(42);
   hPionMinusMinBiasSysErrRel->GetZaxis()->SetLabelSize(0.035);
   hPionMinusMinBiasSysErrRel->GetZaxis()->SetTitleSize(0.035);
   hPionMinusMinBiasSysErrRel->GetZaxis()->SetTitleFont(42);
   hPionMinusMinBiasSysErrRel->Draw("");
   Double_t xAxis4[36] = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3}; 
   
   TH1D *hPionPlusMinBiasSysErrRel = new TH1D("hPionPlusMinBiasSysErrRel","hPionPlusMinBiasSysErrRel",35, xAxis4);
   hPionPlusMinBiasSysErrRel->SetBinContent(5,0.01868457);
   hPionPlusMinBiasSysErrRel->SetBinContent(6,0.0161216);
   hPionPlusMinBiasSysErrRel->SetBinContent(7,0.01431061);
   hPionPlusMinBiasSysErrRel->SetBinContent(8,0.01278948);
   hPionPlusMinBiasSysErrRel->SetBinContent(9,0.01157455);
   hPionPlusMinBiasSysErrRel->SetBinContent(10,0.01064357);
   hPionPlusMinBiasSysErrRel->SetBinContent(11,0.01038358);
   hPionPlusMinBiasSysErrRel->SetBinContent(12,0.01040321);
   hPionPlusMinBiasSysErrRel->SetBinContent(13,0.01065067);
   hPionPlusMinBiasSysErrRel->SetBinContent(14,0.01073758);
   hPionPlusMinBiasSysErrRel->SetBinContent(15,0.01093908);
   hPionPlusMinBiasSysErrRel->SetBinContent(16,0.0110599);
   hPionPlusMinBiasSysErrRel->SetBinContent(17,0.01122509);
   hPionPlusMinBiasSysErrRel->SetBinContent(18,0.01136315);
   hPionPlusMinBiasSysErrRel->SetBinContent(19,0.01137771);
   hPionPlusMinBiasSysErrRel->SetBinContent(20,0.01137771);
   hPionPlusMinBiasSysErrRel->SetBinContent(21,0.01137771);
   hPionPlusMinBiasSysErrRel->SetBinContent(22,0.01137771);
   hPionPlusMinBiasSysErrRel->SetBinContent(23,0.01137771);
   hPionPlusMinBiasSysErrRel->SetBinContent(24,0.01137771);
   hPionPlusMinBiasSysErrRel->SetBinContent(25,0.01137771);
   hPionPlusMinBiasSysErrRel->SetMinimum(0.001219293);
   hPionPlusMinBiasSysErrRel->SetMaximum(0.05114817);
   hPionPlusMinBiasSysErrRel->SetEntries(21);
   hPionPlusMinBiasSysErrRel->SetStats(0);
   hPionPlusMinBiasSysErrRel->SetLineColor(2);
   hPionPlusMinBiasSysErrRel->SetLineWidth(2);
   hPionPlusMinBiasSysErrRel->SetMarkerColor(2);
   hPionPlusMinBiasSysErrRel->GetXaxis()->SetTitle("p_T");
   hPionPlusMinBiasSysErrRel->GetXaxis()->SetRange(3,26);
   hPionPlusMinBiasSysErrRel->GetXaxis()->SetLabelFont(42);
   hPionPlusMinBiasSysErrRel->GetXaxis()->SetLabelSize(0.035);
   hPionPlusMinBiasSysErrRel->GetXaxis()->SetTitleSize(0.035);
   hPionPlusMinBiasSysErrRel->GetXaxis()->SetTitleFont(42);
   hPionPlusMinBiasSysErrRel->GetYaxis()->SetTitle("rel. sys. err.");
   hPionPlusMinBiasSysErrRel->GetYaxis()->SetLabelFont(42);
   hPionPlusMinBiasSysErrRel->GetYaxis()->SetLabelSize(0.035);
   hPionPlusMinBiasSysErrRel->GetYaxis()->SetTitleSize(0.035);
   hPionPlusMinBiasSysErrRel->GetYaxis()->SetTitleFont(42);
   hPionPlusMinBiasSysErrRel->GetZaxis()->SetLabelFont(42);
   hPionPlusMinBiasSysErrRel->GetZaxis()->SetLabelSize(0.035);
   hPionPlusMinBiasSysErrRel->GetZaxis()->SetTitleSize(0.035);
   hPionPlusMinBiasSysErrRel->GetZaxis()->SetTitleFont(42);
   hPionPlusMinBiasSysErrRel->Draw("same");
   Double_t xAxis5[36] = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3}; 
   
   TH1D *hKaonMinusMinBiasSysErrRel = new TH1D("hKaonMinusMinBiasSysErrRel","hKaonMinusMinBiasSysErrRel",35, xAxis5);
   hKaonMinusMinBiasSysErrRel->SetBinContent(6,0.02101081);
   hKaonMinusMinBiasSysErrRel->SetBinContent(7,0.01966228);
   hKaonMinusMinBiasSysErrRel->SetBinContent(8,0.01836725);
   hKaonMinusMinBiasSysErrRel->SetBinContent(9,0.01654397);
   hKaonMinusMinBiasSysErrRel->SetBinContent(10,0.01492521);
   hKaonMinusMinBiasSysErrRel->SetBinContent(11,0.01329219);
   hKaonMinusMinBiasSysErrRel->SetBinContent(12,0.01197453);
   hKaonMinusMinBiasSysErrRel->SetBinContent(13,0.01136758);
   hKaonMinusMinBiasSysErrRel->SetBinContent(14,0.01149933);
   hKaonMinusMinBiasSysErrRel->SetBinContent(15,0.01209399);
   hKaonMinusMinBiasSysErrRel->SetBinContent(16,0.0129423);
   hKaonMinusMinBiasSysErrRel->SetBinContent(17,0.01319937);
   hKaonMinusMinBiasSysErrRel->SetBinContent(18,0.01338563);
   hKaonMinusMinBiasSysErrRel->SetBinContent(19,0.01305338);
   hKaonMinusMinBiasSysErrRel->SetBinContent(20,0.01305338);
   hKaonMinusMinBiasSysErrRel->SetBinContent(21,0.01305338);
   hKaonMinusMinBiasSysErrRel->SetBinContent(22,0.01305338);
   hKaonMinusMinBiasSysErrRel->SetBinContent(23,0.01305338);
   hKaonMinusMinBiasSysErrRel->SetMinimum(-0.2);
   hKaonMinusMinBiasSysErrRel->SetMaximum(0.2);
   hKaonMinusMinBiasSysErrRel->SetEntries(18);
   hKaonMinusMinBiasSysErrRel->SetStats(0);

   ci = TColor::GetColor("#0099ff");
   hKaonMinusMinBiasSysErrRel->SetLineColor(ci);
   hKaonMinusMinBiasSysErrRel->SetLineWidth(2);
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
   hKaonMinusMinBiasSysErrRel->Draw("same");
   Double_t xAxis6[36] = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3}; 
   
   TH1D *hKaonPlusMinBiasSysErrRel = new TH1D("hKaonPlusMinBiasSysErrRel","hKaonPlusMinBiasSysErrRel",35, xAxis6);
   hKaonPlusMinBiasSysErrRel->SetBinContent(6,0.004342853);
   hKaonPlusMinBiasSysErrRel->SetBinContent(7,0.004994155);
   hKaonPlusMinBiasSysErrRel->SetBinContent(8,0.005399486);
   hKaonPlusMinBiasSysErrRel->SetBinContent(9,0.005463311);
   hKaonPlusMinBiasSysErrRel->SetBinContent(10,0.005426551);
   hKaonPlusMinBiasSysErrRel->SetBinContent(11,0.005190788);
   hKaonPlusMinBiasSysErrRel->SetBinContent(12,0.004973586);
   hKaonPlusMinBiasSysErrRel->SetBinContent(13,0.004767482);
   hKaonPlusMinBiasSysErrRel->SetBinContent(14,0.004511727);
   hKaonPlusMinBiasSysErrRel->SetBinContent(15,0.00452386);
   hKaonPlusMinBiasSysErrRel->SetBinContent(16,0.004648428);
   hKaonPlusMinBiasSysErrRel->SetBinContent(17,0.004730747);
   hKaonPlusMinBiasSysErrRel->SetBinContent(18,0.004991366);
   hKaonPlusMinBiasSysErrRel->SetBinContent(19,0.005138622);
   hKaonPlusMinBiasSysErrRel->SetBinContent(20,0.005138622);
   hKaonPlusMinBiasSysErrRel->SetBinContent(21,0.005138622);
   hKaonPlusMinBiasSysErrRel->SetBinContent(22,0.005138622);
   hKaonPlusMinBiasSysErrRel->SetBinContent(23,0.005138622);
   hKaonPlusMinBiasSysErrRel->SetMinimum(-0.2);
   hKaonPlusMinBiasSysErrRel->SetMaximum(0.2);
   hKaonPlusMinBiasSysErrRel->SetEntries(18);
   hKaonPlusMinBiasSysErrRel->SetStats(0);

   ci = TColor::GetColor("#ff9900");
   hKaonPlusMinBiasSysErrRel->SetLineColor(ci);
   hKaonPlusMinBiasSysErrRel->SetLineWidth(2);

   ci = TColor::GetColor("#ff9900");
   hKaonPlusMinBiasSysErrRel->SetMarkerColor(ci);
   hKaonPlusMinBiasSysErrRel->GetXaxis()->SetTitle("p_T");
   hKaonPlusMinBiasSysErrRel->GetXaxis()->SetRange(4,24);
   hKaonPlusMinBiasSysErrRel->GetXaxis()->SetLabelFont(42);
   hKaonPlusMinBiasSysErrRel->GetXaxis()->SetLabelSize(0.035);
   hKaonPlusMinBiasSysErrRel->GetXaxis()->SetTitleSize(0.035);
   hKaonPlusMinBiasSysErrRel->GetXaxis()->SetTitleFont(42);
   hKaonPlusMinBiasSysErrRel->GetYaxis()->SetTitle("rel. sys. err.");
   hKaonPlusMinBiasSysErrRel->GetYaxis()->SetLabelFont(42);
   hKaonPlusMinBiasSysErrRel->GetYaxis()->SetLabelSize(0.035);
   hKaonPlusMinBiasSysErrRel->GetYaxis()->SetTitleSize(0.035);
   hKaonPlusMinBiasSysErrRel->GetYaxis()->SetTitleFont(42);
   hKaonPlusMinBiasSysErrRel->GetZaxis()->SetLabelFont(42);
   hKaonPlusMinBiasSysErrRel->GetZaxis()->SetLabelSize(0.035);
   hKaonPlusMinBiasSysErrRel->GetZaxis()->SetTitleSize(0.035);
   hKaonPlusMinBiasSysErrRel->GetZaxis()->SetTitleFont(42);
   hKaonPlusMinBiasSysErrRel->Draw("same");
   
   TLegend *leg = new TLegend(0.4771887,0.5536723,0.9235512,0.9227872,NULL,"brNDC");
   leg->SetBorderSize(1);
   leg->SetLineColor(0);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("NULL","Inelastic cross section contribution","h");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("hPionMinusMinBiasSysErrRel","hPionMinusMinBiasSysErrRel","lpf");
   entry->SetFillStyle(1001);

   ci = TColor::GetColor("#000099");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(4);
   entry->SetMarkerStyle(1);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("hPionPlusMinBiasSysErrRel","hPionPlusMinBiasSysErrRel","lpf");
   entry->SetFillStyle(1001);
   entry->SetLineColor(2);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(2);
   entry->SetMarkerStyle(1);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("hKaonMinusMinBiasSysErrRel","hKaonMinusMinBiasSysErrRel","lpf");
   entry->SetFillStyle(1001);

   ci = TColor::GetColor("#0099ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(1);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("hKaonPlusMinBiasSysErrRel","hKaonPlusMinBiasSysErrRel","lpf");
   entry->SetFillStyle(1001);

   ci = TColor::GetColor("#ff9900");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#ff9900");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(1);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   leg->Draw();
   Canvas_1->Modified();
   Canvas_1->cd();
   Canvas_1->SetSelected(Canvas_1);
}
