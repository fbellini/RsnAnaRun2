#include "RebinSpectrum.C"
#include "GetPlotRatio.C"
#include "/Users/fbellini/alice/macros/ResonAnT/phiXeXe/LevyTsallis.h"

//published data PbPb2.76 TeV
TGraphAsymmErrors * getPro2phi_PbPb276_1020(Bool_t sys = 0);
TGraphAsymmErrors * getPro2phi_PbPb276_2040(Bool_t sys = 0);
TGraphAsymmErrors * getPro2phi_PbPb276_6080(Bool_t sys = 0);

//preliminary PbPb5 TeV
TH1F * getPro2phi_PbPb502(Bool_t sys = 0, Int_t centLow = 10, Int_t centUp = 20);

//ratios
void ProtonToPhiXeXe();
TH1F * GetProtonToPhiRatio(Int_t icent = 0, Int_t syst = kFALSE);
TH1F * GetPhiXeXeSpectrum(Int_t cent = 0, Bool_t syst = kFALSE);
TH1F * GetProXeXeSpectrum(Int_t cent = 0, Bool_t syst = kFALSE);

TH1F * getRebinnedHisto(TH1F* h1);
TH1F * SumUncertInQuadrature(TH1F * hstat = 0, TH1F* hsys = 0);


void BaryonToMeson(Bool_t plotPbPb502 = 1)
{
  SetStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetLineWidth(1);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleOffset(0.8,"y");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  TGaxis::SetMaxDigits(4);

  TH1F * hFrame = new TH1F("hFrame", "hFrame; #it{p}_{T} (GeV/#it{c}); Ratio", 50, 0., 5.0);
  hFrame->GetYaxis()->SetRangeUser(0.01, 2.0);
  hFrame->GetYaxis()->SetTitleOffset(1.3);
  hFrame->GetXaxis()->SetNdivisions(509);

  //----------------
  // proton to phi
  //----------------
  TFile * fin = TFile::Open("protonToPhiXeXe.root");
  TH1F * hRatioProtonToPhi[2];
  Float_t scalingFactor = 1./10.;
  hRatioProtonToPhi[0] = (TH1F*)fin->Get("ProtonToPhiXeXe_0_stat");
  hRatioProtonToPhi[0]->Scale(scalingFactor);
  hRatioProtonToPhi[1] = (TH1F*)fin->Get("ProtonToPhiXeXe_0_syst");
  hRatioProtonToPhi[1]->Scale(scalingFactor);
  Beautify(hRatioProtonToPhi[0], kRed+1, 1, 2, 20, 1.3);  
  Beautify(hRatioProtonToPhi[1], kRed+1, 1, 2, 20, 1.3);

  TH1F * hPro2phi_PbPb502_1020_stat = (TH1F*) getPro2phi_PbPb502(0, 10, 20);
  TH1F * hPro2phi_PbPb502_1020_syst = (TH1F*) getPro2phi_PbPb502(1, 10, 20);
  hPro2phi_PbPb502_1020_stat->Scale(scalingFactor);
  hPro2phi_PbPb502_1020_syst->Scale(scalingFactor);
  Beautify(hPro2phi_PbPb502_1020_stat, kRed-5, 1, 2, 24, 1.3);
  Beautify(hPro2phi_PbPb502_1020_syst, kRed-9, 1, 2, 24, 1.3);
  hPro2phi_PbPb502_1020_syst->SetFillStyle(1001);
  hPro2phi_PbPb502_1020_syst->SetFillColorAlpha(kRed-9, 0.5);

  //----------------
  //Lambda to K0s
  //----------------
  TFile * finStrange = TFile::Open("LambdaOverK0S_20180420143637.root");
  TH1F * hRatioLambdaToK0s[2];
  hRatioLambdaToK0s[0] = (TH1F*) ((TH1F*) finStrange->Get("h_LamSumOv2K0S[0]"))->Clone("L2K0s_XeXe544TeV_010_stat");
  hRatioLambdaToK0s[0]->Add((TH1F*) finStrange->Get("h_LamSumOv2K0S[1]"));
  hRatioLambdaToK0s[0]->Scale(0.5);
  hRatioLambdaToK0s[1] = (TH1F*) ((TH1F*) finStrange->Get("h_LamSumOv2K0S_syst[0]"))->Clone("L2K0s_XeXe544TeV_010_syst");
  hRatioLambdaToK0s[1]->Add((TH1F*) finStrange->Get("h_LamSumOv2K0S_syst[1]"));
  hRatioLambdaToK0s[1]->Scale(0.5);
  Beautify(hRatioLambdaToK0s[0], kBlue, 1, 2, 29, 1.7);  
  Beautify(hRatioLambdaToK0s[1], kBlue, 1, 2, 29, 1.7);

  TFile * finStrangePb = TFile::Open("Preliminary_Lambda2K0SRatio_PbPb502TeV.root");
  TDirectoryFile * dir = (TDirectoryFile *) finStrangePb->Get("lCentrality10-20%");
  TH1F * hRatioLambdaToK0s_PbPb502[2];
  hRatioLambdaToK0s_PbPb502[0] = (TH1F*) ((TH1F*) dir->Get("fHistLKStat2"))->Clone("L2K0s_PbPb502TeV_1020_stat");
  hRatioLambdaToK0s_PbPb502[1] = (TH1F*) ((TH1F*) dir->Get("fHistLKSyst2"))->Clone("L2K0s_PbPb502TeV_1020_stat");
  Beautify(hRatioLambdaToK0s_PbPb502[0], kBlue-6, 1, 2, 30, 1.7);  
  Beautify(hRatioLambdaToK0s_PbPb502[1], kBlue-9, 1, 2, 30, 1.7);
  hRatioLambdaToK0s_PbPb502[1]->SetFillStyle(1001);
  hRatioLambdaToK0s_PbPb502[1]->SetFillColorAlpha(kBlue-9, 0.4);
  //---------------
  // proton to pion
  //---------------
  TFile * finProton = TFile::Open("Preliminary_ProtonToPionVsPt_XeXe.root");
  TCanvas * cProton = (TCanvas *) finProton->Get("ProtonToPionTwoEnergies");
  //XeXe
  TH1F * hRatioProtonToPion[2];
  hRatioProtonToPion[0] = (TH1F*) ((TH1F*) cProton->FindObject("hSpectraSummedProton_XeXe_Combined_0.00to5.00_Stat"))->Clone("ProtonToPion_010_stat");
  hRatioProtonToPion[0]->Add((TH1F*) cProton->FindObject("hSpectraSummedProton_XeXe_Combined_5.00to10.00_Stat"));
  hRatioProtonToPion[0]->Scale(0.5);
  hRatioProtonToPion[1] = (TH1F*) ((TH1F*) cProton->FindObject("hSpectraSummedProton_XeXe_Combined_0.00to5.00_Syst"))->Clone("ProtonToPion_010_syst");
  hRatioProtonToPion[1]->Add((TH1F*) cProton->FindObject("hSpectraSummedProton_XeXe_Combined_5.00to10.00_Syst"));
  hRatioProtonToPion[1]->Scale(0.5);
  Beautify(hRatioProtonToPion[0], kBlack, 1, 2, 21, 1.3);  
  Beautify(hRatioProtonToPion[1], kBlack, 1, 2, 21, 1.3);

  TH1F * hRatioProtonToPion_PbPb502[2];
  hRatioProtonToPion_PbPb502[0] = (TH1F*) ((TH1F*) cProton->FindObject("hProtonToPion_1020"))->Clone("ProtonToPion_1020_stat");
  hRatioProtonToPion_PbPb502[1] = (TH1F*) ((TH1F*) cProton->FindObject("hSystProtonToPion_1020"))->Clone("ProtonToPion_1020_syst");
  Beautify(hRatioProtonToPion_PbPb502[0], kGray+2, 1, 2, 25, 1.3);  
  Beautify(hRatioProtonToPion_PbPb502[1], kGray, 1, 2, 25, 1.3);
  hRatioProtonToPion_PbPb502[1]->SetFillStyle(1001);
  hRatioProtonToPion_PbPb502[1]->SetFillColorAlpha(kGray+1, 0.4);
  
  //legend for most central XeXe
  TLegend * leg1 = new TLegend(0.18, 0.57, 0.42, 0.72); //"Xe-Xe #sqrt{#it{s}_{NN}} = 5.44 TeV, 0-10%"
  myLegendSetUp(leg1, 0.035);
  leg1->AddEntry(hRatioProtonToPhi[0], Form("(p+#bar{p})/#phi #times %2.1f", scalingFactor), "p");
  leg1->AddEntry(hRatioProtonToPion[0], "(p+#bar{p})/(#pi^{+}+#pi^{-})", "p");
  leg1->AddEntry(hRatioLambdaToK0s[0], "(#Lambda+#bar{#Lambda})/2K^{0}_{S}", "p");

  //legend for central PbPb
  // TLegend * legPb = new TLegend(0.55, 0.8, 0.95, 0.92,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, 10-20%");
  // myLegendSetUp(legPb, 0.03);
  // legPb->AddEntry(hPro2phi_PbPb502_1020_stat, Form("(p+#bar{p})/#phi #times %2.1f", scalingFactor), "p");
  // legPb->AddEntry(hRatioLambdaToK0s_PbPb502[0], "(#Lambda+#bar{#Lambda})/2K^{0}_{S}", "p");
  // legPb->AddEntry(hRatioProtonToPion_PbPb502[0], "(p+#bar{p})/(#pi^{+}+#pi^{-})", "p");

  TLegend *l1 = new TLegend(0.13,0.75,0.42,0.9);
  myLegendSetUp(l1, 0.035);
  l1->AddEntry((TObject*)0,"#bf{ALICE Preliminary}","");
  l1->AddEntry((TObject*)0,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, 10-20% (open markers)","");
  l1->AddEntry((TObject*)0,"Xe-Xe #sqrt{#it{s}_{NN}} = 5.44 TeV, 0-10% (full markers)","");
  l1->Draw();
  
  TPaveText * pave = new TPaveText(0.18, 0.65, 0.6, 0.7, "NDC");
  pave->SetBorderSize(0);
  pave->SetFillStyle(0);
  pave->SetTextSize(0.035);
  pave->SetTextAlign(12);
  pave->AddText("#bf{ALICE Preliminary}");
  
  TCanvas * c1 = new TCanvas("B2M010","B2M010", 800, 800);
  c1->cd();
  hFrame->Draw();
  if (plotPbPb502){
    hRatioLambdaToK0s_PbPb502[1]->Draw("E2same");
    hRatioLambdaToK0s_PbPb502[0]->Draw("E X0same");
    hRatioProtonToPion_PbPb502[1]->Draw("E2same");
    hRatioProtonToPion_PbPb502[0]->Draw("E X0same");
    hPro2phi_PbPb502_1020_syst->Draw("E2 same");
    hPro2phi_PbPb502_1020_stat->Draw("E X0same");
  }
  // hRatioLambdaToK0s[1]->Draw("E2 same");
  // hRatioLambdaToK0s[0]->Draw("same");
  hRatioProtonToPion[1]->Draw("E2 same");
  hRatioProtonToPion[0]->Draw("E X0 same");
  hRatioProtonToPhi[1]->Draw("E2 same");
  hRatioProtonToPhi[0]->Draw("E X0 same");
  //pave->Draw();
  l1->Draw();
  leg1->Draw();
  //if (plotPbPb502) legPb->Draw();
  c1->Print("BaryonToMeson_XeXe010.pdf");
  c1->Print("BaryonToMeson_XeXe010.eps");

  return;
}

void ProtonToPhiXeXe(){
  
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetLineWidth(1);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleOffset(0.8,"y");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  SetStyle();
  TGaxis::SetMaxDigits(4);

  TFile * fout = new TFile("protonToPhiXeXe.root","recreate");
  TH1F * hRatioProtonToPhi[4][2];
  TH1F * hFrame = new TH1F("hFrame", "hFrame; #it{p}_{T} (GeV/#it{c}); (p+#bar{p}) / #phi", 50, 0., 5.0);
  hFrame->GetYaxis()->SetRangeUser(0.1, 11.9);
  hFrame->GetXaxis()->SetNdivisions(509);

  for (Int_t j=0; j<4; j++) {
    for (Int_t i=0; i<2; i++) {
      TString name;
      if (i>0) name = Form("ProtonToPhiXeXe_%i_syst", j);
      else name = Form("ProtonToPhiXeXe_%i_stat", j);
      hRatioProtonToPhi[j][i] = (TH1F*)  ((TH1F*) GetProtonToPhiRatio(j, i))->Clone(name.Data());
      hRatioProtonToPhi[j][i]->SetTitle(name.Data());
      fout->cd();
      hRatioProtonToPhi[j][i]->Write();
     }
  }

  Beautify(hRatioProtonToPhi[0][0], kRed+1, 1, 2, 20, 1.3);  
  Beautify(hRatioProtonToPhi[0][1], kRed+1, 1, 2, 20, 1.3);
  hRatioProtonToPhi[0][1]->SetFillStyle(0);
  hRatioProtonToPhi[0][1]->SetFillColorAlpha(kRed, 0.3);

  Beautify(hRatioProtonToPhi[1][0], kOrange-3, 1, 2, 20, 1.3);  
  Beautify(hRatioProtonToPhi[1][1], kOrange-3, 1, 2, 20, 1.3);
  hRatioProtonToPhi[1][1]->SetFillStyle(0);
  hRatioProtonToPhi[1][1]->SetFillColorAlpha(kOrange, 0.3);
  
  Beautify(hRatioProtonToPhi[3][0], kBlue+1, 1, 2, 20, 1.3);  
  Beautify(hRatioProtonToPhi[3][1], kBlue+1, 1, 2, 20, 1.3);
  hRatioProtonToPhi[3][1]->SetFillStyle(0);
  hRatioProtonToPhi[3][1]->SetFillColorAlpha(kBlue, 0.2);

  //plot other systems/energies
  TGraphAsymmErrors * gPro2phi_PbPb276_1020_stat = (TGraphAsymmErrors *) getPro2phi_PbPb276_1020(0);
  TGraphAsymmErrors * gPro2phi_PbPb276_1020_syst = (TGraphAsymmErrors *) getPro2phi_PbPb276_1020(1);
  BeautifyGraphAsymmErrors(gPro2phi_PbPb276_1020_stat, kBlack, kBlack, 0, 1, 1, 24, 1.3); 
  BeautifyGraphAsymmErrors(gPro2phi_PbPb276_1020_syst, kBlack, kBlack, 0, 1, 1, 24, 1.3); 

  TGraphAsymmErrors * gPro2phi_PbPb276_2040_stat = (TGraphAsymmErrors *) getPro2phi_PbPb276_2040(0);
  TGraphAsymmErrors * gPro2phi_PbPb276_2040_syst = (TGraphAsymmErrors *) getPro2phi_PbPb276_2040(1);
  BeautifyGraphAsymmErrors(gPro2phi_PbPb276_2040_stat, kBlack, kBlack, 0, 1, 1, 24, 1.3); 
  BeautifyGraphAsymmErrors(gPro2phi_PbPb276_2040_syst, kBlack, kBlack, 0, 1, 1, 24, 1.3); 

  TGraphAsymmErrors * gPro2phi_PbPb276_6080_stat = (TGraphAsymmErrors *) getPro2phi_PbPb276_6080(0);
  TGraphAsymmErrors * gPro2phi_PbPb276_6080_syst = (TGraphAsymmErrors *) getPro2phi_PbPb276_6080(1);
  BeautifyGraphAsymmErrors(gPro2phi_PbPb276_6080_stat, kBlack, kBlack, 0, 1, 1, 24, 1.3); 
  BeautifyGraphAsymmErrors(gPro2phi_PbPb276_6080_syst, kBlack, kBlack, 0, 1, 1, 24, 1.3); 

  TH1F * hPro2phi_PbPb502_1020_stat = (TH1F*) getPro2phi_PbPb502(0, 10, 20);
  TH1F * hPro2phi_PbPb502_1020_syst = (TH1F*) getPro2phi_PbPb502(1, 10, 20);
  Beautify(hPro2phi_PbPb502_1020_stat, kRed-6, 1, 2, 25, 1.3);
  Beautify(hPro2phi_PbPb502_1020_syst, kRed-9, 1, 2, 25, 1.3);
  hPro2phi_PbPb502_1020_syst->SetFillStyle(1001);
  hPro2phi_PbPb502_1020_syst->SetFillColorAlpha(kRed-9, 0.5);
  
  TH1F * hPro2phi_PbPb502_3040_stat = (TH1F*) getPro2phi_PbPb502(0, 30, 40);
  TH1F * hPro2phi_PbPb502_3040_syst = (TH1F*) getPro2phi_PbPb502(1, 30, 40);
  Beautify(hPro2phi_PbPb502_3040_stat, kOrange-1, 1, 2, 25, 1.3);
  Beautify(hPro2phi_PbPb502_3040_syst, kOrange, 1, 2, 25, 1.3);
  hPro2phi_PbPb502_3040_syst->SetFillStyle(1001);
  hPro2phi_PbPb502_3040_syst->SetFillColorAlpha(kOrange, 0.6);
  
  TH1F * hPro2phi_PbPb502_7080_stat = (TH1F*) getPro2phi_PbPb502(0, 70, 80);
  TH1F * hPro2phi_PbPb502_7080_syst = (TH1F*) getPro2phi_PbPb502(1, 70, 80);
  Beautify(hPro2phi_PbPb502_7080_stat, kBlue-6, 1, 2, 25, 1.3);
  Beautify(hPro2phi_PbPb502_7080_syst, kBlue-9, 1, 2, 25, 1.3);
  hPro2phi_PbPb502_7080_syst->SetFillStyle(1001);
  hPro2phi_PbPb502_7080_syst->SetFillColorAlpha(kBlue-9, 0.4);
  //most central XeXe
  TLegend * leg1 = new TLegend(0.22, 0.7, 0.6, 0.9, "#bf{ALICE Preliminary}");
  myLegendSetUp(leg1, 0.04);
  leg1->AddEntry(hRatioProtonToPhi[0][0], "Xe-Xe #sqrt{#it{s}_{NN}} = 5.44 TeV, 0-10%", "p");
  leg1->AddEntry(hPro2phi_PbPb502_1020_stat, "Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, 10-20%", "p");
  leg1->AddEntry(gPro2phi_PbPb276_1020_stat, "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV, 10-20%", "p");
  
  TCanvas * c1 = new TCanvas("p2phi010","p2Phi", 800, 800);
  c1->cd();
  hFrame->Draw();
  gPro2phi_PbPb276_1020_syst->Draw("e2 same");
  gPro2phi_PbPb276_1020_stat->Draw("pz same");
  hPro2phi_PbPb502_1020_syst->Draw("E2 same");
  hPro2phi_PbPb502_1020_stat->Draw("same");
  hRatioProtonToPhi[0][1]->Draw("E2 same");
  hRatioProtonToPhi[0][0]->Draw("X0same");
  leg1->Draw();
  c1->Print("ProtonToPhi_XeXe010.pdf");
  c1->Print("ProtonToPhi_XeXe010.eps");
  
  //10-30 central XeXe
  TLegend * leg2 = new TLegend(0.22, 0.7, 0.6, 0.9, "#bf{ALICE Preliminary}");
  myLegendSetUp(leg2, 0.04);
  leg2->AddEntry(hRatioProtonToPhi[1][0], "Xe-Xe #sqrt{#it{s}_{NN}} = 5.44 TeV, 10-30%", "p");
  leg2->AddEntry(hPro2phi_PbPb502_3040_stat, "Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, 30-40%", "p");
  leg2->AddEntry(gPro2phi_PbPb276_1020_stat, "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV, 20-40%", "p");
  
  TCanvas * c2 = new TCanvas("p2phi1030","p2Phi", 800, 800);
  c2->cd();
  hFrame->Draw();
  gPro2phi_PbPb276_2040_syst->Draw("e2 same");
  gPro2phi_PbPb276_2040_stat->Draw("pz same");
  hPro2phi_PbPb502_3040_syst->Draw("E2 same");
  hPro2phi_PbPb502_3040_stat->Draw("same");
  hRatioProtonToPhi[1][1]->Draw("E2 same");
  hRatioProtonToPhi[1][0]->Draw("X0same");
  leg2->Draw();
  c2->Print("ProtonToPhi_XeXe1030.pdf");
  c2->Print("ProtonToPhi_XeXe1030.eps");
  
  //60-90 XeXe
  TLegend * leg4 = new TLegend(0.22, 0.7, 0.6, 0.9, "#bf{ALICE Preliminary}");
  myLegendSetUp(leg4, 0.04);
  leg4->AddEntry(hRatioProtonToPhi[3][0], "Xe-Xe #sqrt{#it{s}_{NN}} = 5.44 TeV, 60-90%", "p");
  leg4->AddEntry(hPro2phi_PbPb502_7080_stat, "Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, 70-80%", "p");
  leg4->AddEntry(gPro2phi_PbPb276_6080_stat, "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV, 60-80%", "p");
  
  TCanvas * c4 = new TCanvas("p2phi6090","p2Phi", 800, 800);
  c4->cd();
  hFrame->Draw();
  gPro2phi_PbPb276_6080_syst->Draw("e2 same");
  gPro2phi_PbPb276_6080_stat->Draw("pz same");
  hPro2phi_PbPb502_7080_syst->Draw("E2 same");
  hPro2phi_PbPb502_7080_stat->Draw("same");
  hRatioProtonToPhi[3][1]->Draw("E2 same");
  hRatioProtonToPhi[3][0]->Draw("X0same");
  leg4->Draw();
  c4->Print("ProtonToPhi_XeXe6090.pdf");
  c4->Print("ProtonToPhi_XeXe6090.eps");
  
  return;
}


TH1F * GetProtonToPhiRatio(Int_t icent, Int_t syst){
  
  Int_t centrality[] = {0, 10, 30, 60, 90};
  Color_t color[4] = {kRed+2, kSpring+3, kBlue+2, kBlack};
  Int_t centLow = centrality[icent];
  Int_t centHigh = centrality[icent];
  
  //--------------------------
  // PHI
  //--------------------------
  TH1F * phiXeXe[1];
  phiXeXe[0] = (TH1F *) GetPhiXeXeSpectrum(icent, syst);
  
  //--------------------------
  //PROTONS with rebin
  //--------------------------
  TH1F * proXeXe[2];   
  proXeXe[0] = GetProXeXeSpectrum(icent, 0);
  proXeXe[1] = GetProXeXeSpectrum(icent, 1);
  
  TH1F * proXeXeRebin[1];
  proXeXeRebin[0] = (TH1F *) phiXeXe[0]->Clone(Form("p_rebinned_%i", icent));
  proXeXeRebin[0]->Reset("ICES");

  TH1F * proXeXeTot = (TH1F *) SumUncertInQuadrature(proXeXe[0], proXeXe[1]);
  TF1 * myLevyTsallis = LevyTsallis("levy", 0.938);
  RebinSpectrum(proXeXe, proXeXeRebin, myLevyTsallis, 1, 1, proXeXeTot);

  //--------------------------
  //PROTON / PHI
  //--------------------------
  TH1F* ratio = (TH1F*) GetPlotRatio(proXeXeRebin[0], phiXeXe[0]);
  ratio->GetYaxis()->SetRangeUser(0.0, 12.);
  ratio->GetXaxis()->SetRangeUser(0.0, 5.0);
  ratio->GetYaxis()->SetTitle("(p+#bar{p}) / #phi");
  return ratio;

  // TCanvas * cr = 0x0;
  // TPad * pad1 = 0x0;
  // TPad * pad2 = 0x0;
  // cr = new TCanvas("cr","compare",800, 1000);
  // cr->cd();
  // pad1 = new TPad("pad1","This is pad1",0.001,0.5,0.999,0.999);
  // pad1->SetFillColor(0);
  // pad1->SetBorderMode(0);
  // pad1->SetBorderSize(0);
  // pad1->SetMargin(0.15,0.05,0.01,0.02);

  // pad2 = new TPad("pad2","This is pad2",0.001,0.001, 0.999,0.5);
  // pad2->SetFillColor(0);
  // pad2->SetBorderMode(0);
  // pad2->SetBorderSize(0);
  // pad2->SetMargin(0.15,0.05,0.2,0.01);
  // pad1->Draw();
  // pad2->Draw();
  // pad1->cd();
  // pad1->SetLogy();

  // TLegend * leg = new TLegend(0.2, 0.15,0.4,0.30, "Xe-Xe, #sqrt{#it{s}_{NN}} = 5.44 TeV");
  // leg->AddEntry(proXeXe[0], "p", "p");
  // leg->AddEntry(proXeXeRebin[0], "p, rebinned", "p");
  // leg->AddEntry(phiXeXe[0], "#phi", "p");
  // leg->SetFillStyle(0);
  // leg->SetBorderSize(0);
  // leg->SetTextFont(42);
  // leg->SetTextSize(0.05);

  // Beautify(proXeXe[0], kGray, 1, 2, 24, 1.0);
  // Beautify(proXeXe[1], kGray, 1, 2, 24, 1.0);
  // Beautify(proXeXeRebin[0], kRed, 1, 2, 20, 1.0);
  // Beautify(phiXeXe[0], kBlack, 1, 2, 28, 1.0);
  // Beautify(ratio, color[icent], 1, 2, 20, 1.0);
  // ratio->SetFillStyle(0);

  // pad1->cd();
  // proXeXe[0]->GetYaxis()->SetRangeUser(1e-4, 100.);
  // proXeXe[0]->GetXaxis()->SetRangeUser(0, 5.);
  // proXeXe[1]->Draw("E2");
  // proXeXe[0]->Draw("same");
  // proXeXeRebin[0]->Draw("same");
  // phiXeXe[0]->Draw("same");
  // leg->Draw();
  
  // pad2->cd();  
  // ratio->GetYaxis()->SetRangeUser(0.0, 12.);
  // ratio->GetXaxis()->SetRangeUser(0.0, 5.0);
  // ratio->GetYaxis()->SetTitle("(p+#bar{p}) / #phi");
  //ratio->Draw();
 
}


TH1F * SumUncertInQuadrature(TH1F * hstat, TH1F* hsys)
{
  if (!hstat || !hsys) return 0;
  TH1F * hnew = (TH1F *) hstat->Clone();

  for (int i = 1; i<hnew->GetNbinsX()+1; i++){
    Float_t errInQuad =  TMath::Sqrt(TMath::Power(hstat->GetBinError(i), 2.0) + TMath::Power(hsys->GetBinError(i), 2.0));
    hnew->SetBinError(i, errInQuad);
  }
  return hnew;
}

TH1F* GetPhiXeXeSpectrum(Int_t icent, Bool_t syst)
{
  //
  // phi spectra in XeXe
  //
  Int_t cent[] = {0, 10, 30, 60, 90};
  TFile* fin = TFile::Open("Preliminary_Spectra_phi_XeXe544TeV.root");
  if (!fin) {Printf("Invalid file for phi in XeXe."); return 0x0;}
  Printf(":::: File %s", fin->GetName());
  TH1F * spectrum;
  TString hname;
  if (syst) hname = Form("hCorrected_%i%i_syst", icent, icent);
  else hname = Form("hCorrected_%i", icent);
  spectrum = (TH1F*) fin->Get(hname.Data());
  if (!spectrum) Printf("Invalid phi XeXe spectrum");
  return spectrum;
}


TH1F* GetProXeXeSpectrum(Int_t cent, Bool_t syst)
{
  //
  // proton spectrum in XeXe
  //
  TFile *fin = TFile::Open("Preliminary_Spectra_XeXeLHC17n_Combined_Histograms.root");
  if (!fin) {Printf("Invalid file for protons in XeXe."); return 0x0;}
  TList * lin;
  if (syst) lin = (TList *) fin->Get("Summed_Proton_Sys");
  else lin = (TList *) fin->Get("Summed_Proton");
  TH1F * spectrum = 0x0;
  TH1F * dummy = 0x0;
  switch (cent)
    {
    case 0:
      spectrum = (TH1F*) ((TH1F*) lin->FindObject("hSpectraSummedProton_XeXe_Combined_0.00to5.00"))->Clone("hSpectraSummedProton_XeXe_Combined_0");
      dummy = (TH1F*) lin->FindObject("hSpectraSummedProton_XeXe_Combined_5.00to10.00");
      spectrum->Add(dummy, 1.0);
      spectrum->Scale(0.5);
      break;
    case 1:
      spectrum = (TH1F*) ((TH1F*) lin->FindObject("hSpectraSummedProton_XeXe_Combined_10.00to20.00"))->Clone("hSpectraSummedProton_XeXe_Combined_1");
      dummy = (TH1F*) lin->FindObject("hSpectraSummedProton_XeXe_Combined_20.00to30.00");
      spectrum->Add(dummy, 1.0);
      spectrum->Scale(0.5);
      break;
    case 2:
      spectrum = (TH1F*) ((TH1F*) lin->FindObject("hSpectraSummedProton_XeXe_Combined_30.00to40.00"))->Clone("hSpectraSummedProton_XeXe_Combined_2");
      dummy = (TH1F*) lin->FindObject("hSpectraSummedProton_XeXe_Combined_40.00to50.00");
      spectrum->Add(dummy, 1.0);
      dummy = (TH1F*) lin->FindObject("hSpectraSummedProton_XeXe_Combined_50.00to60.00");
      spectrum->Add(dummy, 1.0);
      spectrum->Scale(1./3.);
      break;
    case 3:
      spectrum = (TH1F*) ((TH1F*) lin->FindObject("hSpectraSummedProton_XeXe_Combined_60.00to70.00"))->Clone("hSpectraSummedProton_XeXe_Combined_3");
      dummy = (TH1F*) lin->FindObject("hSpectraSummedProton_XeXe_Combined_70.00to90.00");
      spectrum->Add(dummy, 1.0);
      spectrum->Scale(1./2.);
      break;
    default:
      spectrum = 0x0;
      break;
    }
  
  if (!spectrum) Printf("Invalid phi XeXe spectrum");
  return spectrum;
}
//----------------------------------------------------
TH1F* getRebinnedHisto(TH1F* h1)
{
  if(!h1) return NULL;
  
  Double_t pTbins [] = {0., 0.3, 0.6, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0};
  Int_t npT = sizeof(pTbins) / sizeof(pTbins[0]) - 1;
  TH1F* hResult = new TH1F(Form("%s_rebin", h1->GetName()), h1->GetTitle(), npT, pTbins);
  for (Int_t i = 0; i<npT; i++) {
    Printf(":::: Bin new %i - low edge = %f - up edge = %f", i, pTbins[i], pTbins[i+1]);
    Double_t err1 = 0.0;
    Int_t lowEdgeInt = h1->GetXaxis()->FindBin(pTbins[i]);
    Int_t upEdgeInt = h1->GetXaxis()->FindBin(pTbins[i+1])-1;
    Double_t integral = h1->IntegralAndError(lowEdgeInt, upEdgeInt, err1, "width");
    hResult->SetBinContent(i+1, integral/(pTbins[i+1]-pTbins[i]));
    hResult->SetBinError(i+1, err1/(pTbins[i+1]-pTbins[i]));
    Printf(">> Bin old low edge %i - up edge %i ", lowEdgeInt, upEdgeInt);
  }
  
  return hResult;
}

TGraphAsymmErrors * getPro2phi_PbPb276_1020(Bool_t sys)
{
  //http://hepdata.cedar.ac.uk/view/ins1288320/d28
// Plot: p8573_d28x1y1
  double p8573_d28x1y1_xval[] = { 0.65, 0.9, 1.25, 1.75, 2.25, 2.75, 3.5, 4.5 };
  double p8573_d28x1y1_xerrminus[] = { 0.15000000000000002, 0.09999999999999998, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5 };
  double p8573_d28x1y1_xerrplus[] = { 0.15000000000000002, 0.09999999999999998, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5 };
  double p8573_d28x1y1_yval[] = { 4.683, 4.856, 4.812, 4.558, 4.489, 4.905, 3.953, 2.605 };
  double p8573_d28x1y1_yerrminus[] = { 0.49155365932927403, 0.538004646820081, 0.4109622853742178, 0.6024060092661759, 0.5061076960489733, 0.6004639872631831, 0.6478248220005158, 0.591304490089497 };
  double p8573_d28x1y1_yerrplus[] = { 0.49155365932927403, 0.538004646820081, 0.4109622853742178, 0.6024060092661759, 0.5061076960489733, 0.6004639872631831, 0.6478248220005158, 0.591304490089497 };
  double p8573_d28x1y1_ystatminus[] = { 0.227, 0.18, 0.149, 0.227, 0.257, 0.334, 0.229, 0.204 };
  double p8573_d28x1y1_ystatplus[] = { 0.227, 0.18, 0.149, 0.227, 0.257, 0.334, 0.229, 0.204 };
  int p8573_d28x1y1_numpoints = 8;
  TGraphAsymmErrors * graph =0;
  if (sys) graph = new TGraphAsymmErrors(p8573_d28x1y1_numpoints, p8573_d28x1y1_xval, p8573_d28x1y1_yval, p8573_d28x1y1_xerrminus, p8573_d28x1y1_xerrplus, p8573_d28x1y1_yerrminus, p8573_d28x1y1_yerrplus);
  else graph = new TGraphAsymmErrors(p8573_d28x1y1_numpoints, p8573_d28x1y1_xval, p8573_d28x1y1_yval, p8573_d28x1y1_xerrminus, p8573_d28x1y1_xerrplus, p8573_d28x1y1_ystatminus, p8573_d28x1y1_ystatplus);
  graph->SetName("/HepData/8573/d28x1y1");
  graph->SetTitle("/HepData/8573/d28x1y1");
  
  return graph;

}

TGraphAsymmErrors * getPro2phi_PbPb276_2040(Bool_t sys)
{
  //http://hepdata.cedar.ac.uk/view/ins1288320/all
  //http://hepdata.cedar.ac.uk/view/ins1288320/d29
  // Plot: p8573_d27x1y1
  double p8573_d27x1y1_xval[] = { 0.65, 0.9, 1.25, 1.75, 2.25, 2.75, 3.5, 4.5 };
  double p8573_d27x1y1_xerrminus[] = { 0.15000000000000002, 0.09999999999999998, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5 };
  double p8573_d27x1y1_xerrplus[] = { 0.15000000000000002, 0.09999999999999998, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5 };
  double p8573_d27x1y1_yval[] = { 4.901, 4.648, 5.169, 5.154, 5.094, 5.004, 4.559, 3.518 };
  double p8573_d27x1y1_yerrminus[] = { 0.5886662891656019, 0.429288946980935, 0.4868695923961569, 0.61730057508478, 0.5752842775532806, 0.5764234554561429, 0.727303237996367, 0.7842168067569071 };
  double p8573_d27x1y1_yerrplus[] = { 0.5886662891656019, 0.429288946980935, 0.4868695923961569, 0.61730057508478, 0.5752842775532806, 0.5764234554561429, 0.727303237996367, 0.7842168067569071 };
  double p8573_d27x1y1_ystatminus[] = { 0.252, 0.175, 0.261, 0.288, 0.326, 0.35, 0.299, 0.38 };
  double p8573_d27x1y1_ystatplus[] = { 0.252, 0.175, 0.261, 0.288, 0.326, 0.35, 0.299, 0.38 };
  
  int p8573_d27x1y1_numpoints = 8;
  TGraphAsymmErrors * graph = 0x0;
  if (sys) graph = new TGraphAsymmErrors(p8573_d27x1y1_numpoints, p8573_d27x1y1_xval, p8573_d27x1y1_yval, p8573_d27x1y1_xerrminus, p8573_d27x1y1_xerrplus, p8573_d27x1y1_yerrminus, p8573_d27x1y1_yerrplus);
  else graph = new TGraphAsymmErrors(p8573_d27x1y1_numpoints, p8573_d27x1y1_xval, p8573_d27x1y1_yval, p8573_d27x1y1_xerrminus, p8573_d27x1y1_xerrplus, p8573_d27x1y1_ystatminus, p8573_d27x1y1_ystatplus);
  
  graph->SetName("/HepData/8573/d27x1y1");
  graph->SetTitle("/HepData/8573/d27x1y1");
  return graph;
  
}

TGraphAsymmErrors * getPro2phi_PbPb276_6080(Bool_t sys)
{

  double p8573_d31x1y1_xval[] = { 0.65, 0.9, 1.25, 1.75, 2.25, 2.75, 3.5, 4.5 };
  double p8573_d31x1y1_xerrminus[] = { 0.15000000000000002, 0.09999999999999998, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5 };
  double p8573_d31x1y1_xerrplus[] = { 0.15000000000000002, 0.09999999999999998, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5 };
  double p8573_d31x1y1_yval[] = { 6.203, 6.16, 5.748, 4.556, 4.176, 3.542, 2.799, 1.968 };
  double p8573_d31x1y1_yerrminus[] = { 0.7872312239742527, 0.5413917250937624, 0.5805385430787521, 0.6345210792400833, 0.5313915693723414, 0.5183483384751996, 0.5545637925432926, 0.4651892088172295 };
  double p8573_d31x1y1_yerrplus[] = { 0.7872312239742527, 0.5413917250937624, 0.5805385430787521, 0.6345210792400833, 0.5313915693723414, 0.5183483384751996, 0.5545637925432926, 0.4651892088172295 };
  double p8573_d31x1y1_ystatminus[] = { 0.263, 0.217, 0.167, 0.236, 0.259, 0.227, 0.146, 0.151 };
  double p8573_d31x1y1_ystatplus[] = { 0.263, 0.217, 0.167, 0.236, 0.259, 0.227, 0.146, 0.151 };

  int p8573_d31x1y1_numpoints = 8;
  TGraphAsymmErrors * graph = 0x0;
  if (sys)  graph = new TGraphAsymmErrors(p8573_d31x1y1_numpoints, p8573_d31x1y1_xval, p8573_d31x1y1_yval, p8573_d31x1y1_xerrminus, p8573_d31x1y1_xerrplus, p8573_d31x1y1_yerrminus, p8573_d31x1y1_yerrplus);
  else   graph = new TGraphAsymmErrors(p8573_d31x1y1_numpoints, p8573_d31x1y1_xval, p8573_d31x1y1_yval, p8573_d31x1y1_xerrminus, p8573_d31x1y1_xerrplus, p8573_d31x1y1_ystatminus, p8573_d31x1y1_ystatplus);
  graph->SetName("/HepData/8573/d31x1y1");
  graph->SetTitle("/HepData/8573/d31x1y1");
  return graph;
}

  
TH1F * getPro2phi_PbPb502(Bool_t sys, Int_t centLow, Int_t centUp)
{
  //Preliminary p/phi from Ajay /// 01.05.2018
  TFile * fin = TFile::Open("Preliminary_PhiBinProtonSpectra_PbPb5TeV.root");
  if (!fin) return 0;
  TH1F * histo = 0x0;

  if (sys) histo = (TH1F*) fin->Get(Form("ProtonByPhi_syst_V0M%i_%i", centLow, centUp));
  else histo = (TH1F*) fin->Get(Form("ProtonByPhi_stat_V0M%i_%i", centLow, centUp));
  return histo;

}

