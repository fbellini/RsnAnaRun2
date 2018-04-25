#include "RebinSpectrum.C"
#include "GetPlotRatio.C"

Double_t LevyTsallis_Func(const Double_t *x, const Double_t *p);
TF1 * LevyTsallis(const Char_t *name, Double_t mass, Double_t n = 5., Double_t C = 0.1, Double_t norm = 1.);
TH1F * SumUncertInQuadrature(TH1F * hstat = 0, TH1F* hsys = 0);
void PlotRatioProtonToF0(Bool_t correctByBR = 1);
void PlotRatiof0ToPhi(Bool_t correctByBR = 1);

void PlotRatioToPion(Bool_t correctByBR = 1){
  
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetLineWidth(1);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleOffset(0.8,"y");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  TGaxis::SetMaxDigits(4);
  
  TCanvas * cr = 0x0;
  TPad * pad1 = 0x0;
  TPad * pad2 = 0x0;
  cr = new TCanvas("cr","compare",800, 1000);
  cr->cd();
  pad1 = new TPad("pad1","This is pad1",0.001,0.5,0.999,0.999);
  pad1->SetFillColor(0);
  pad1->SetBorderMode(0);
  pad1->SetBorderSize(0);
  pad1->SetMargin(0.15,0.05,0.01,0.02);

  pad2 = new TPad("pad2","This is pad2",0.001,0.001, 0.999,0.5);
  pad2->SetFillColor(0);
  pad2->SetBorderMode(0);
  pad2->SetBorderSize(0);
  pad2->SetMargin(0.15,0.05,0.2,0.01);
  pad1->Draw();
  pad2->Draw();

  pad1->cd();
  pad1->SetLogy();

  //-------------
  // F0
  //-------------
  TFile* myfile=TFile::Open("prel_finalWsys.root");
  if (!myfile) {Printf("Invalid file."); return;}
  TH1F * f0_Tot[1];
  f0_Tot[0] = (TH1F*) myfile->Get("hCorrectedSpectrum_syst_stat");
  if (!correctByBR) f0_Tot[0]->Scale(0.46);
  
  TH1F * rebinned[1];
  rebinned[0] = (TH1F*) f0_Tot[0]->Clone("dummy");
  rebinned[0]->Reset("ICES");
  
  //-------------
  //PIONS
  //-------------
  ///Users/fbellini/alice/pwglf-piKp5teV/SpectraAnalysisRun2/results/spectra/spectra-pag/Preliminaries/SQM2017/
  TString filePiName = "Spectra_ppLHC15n_Combined_Histograms.root";
  TFile *filePi = TFile::Open(filePiName.Data());
  if (!filePi) {Printf("Invalid file."); return;}
  TH1F *pi_syst = (TH1F*) ((TList*) filePi->Get("Summed_Pion_Sys"))->FindObject("hSpectraSummedPion_pp_Combined_MB");
  TH1F *pi_stat = (TH1F*) ((TList*) filePi->Get("Summed_Pion"))->FindObject("hSpectraSummedPion_pp_Combined_MB");
  
  if (!pi_stat) {Printf("missing plot pion data_stat"); return;}
  if (!pi_syst) {Printf("missing plot pion data_syst"); return;}
  
  TH1F * pi_Tot[1];
  pi_Tot[0] = (TH1F*) SumUncertInQuadrature(pi_stat, pi_syst);
  pi_Tot[0]->Scale(0.5); // (pi+ + pi-) / 2
  
  TF1 * myLevyTsallis = LevyTsallis("levy", 0.139);
  RebinSpectrum(pi_Tot, rebinned, myLevyTsallis, 1, 1,  pi_Tot[0]);
  Beautify(rebinned[0], kRed+1, 1, 2, 24, 1.0);
  Beautify(f0_Tot[0], kBlack, 1, 2, 20, 1.0);
  
  TH1F* ratio = (TH1F*) GetPlotRatio(f0_Tot[0], rebinned[0]);//, 0, "", "f_{0}(980)", "(#pi^{+}+#pi^{-})/2.0", "d^{2}#it{N}_{INEL}/(d#it{y}d#it{p}_{T}) (GeV/#it{c})^{-1}");

  TString text="Uncertainties: #sqrt{stat.^{2} + sys.^{2}}";
  TPaveText * pave = new TPaveText(0.4,0.1,0.4,0.14,"NDC");
  pave->SetBorderSize(0);
  pave->SetFillStyle(0);
  pave->SetTextColor(kBlack);
  pave->SetTextFont(42);
  pave->SetTextSize(0.05);
  pave->InsertText(text.Data());

  TLegend * leg = new TLegend(0.4, 0.6,0.8,0.9, "pp, #sqrt{#it{s}} = 5.02 TeV");
  leg->AddEntry(pi_syst, "(#pi^{+}+#pi^{-})/2 ALICE preliminary", "p");
  leg->AddEntry(rebinned[0], "(#pi^{+}+#pi^{-})/2 rebinned", "p");
  leg->AddEntry(f0_Tot[0], "f_{0}(980), this analysis", "p");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  
  TPaveText * pave2 = new TPaveText(0.2, 0.8,0.5,0.9, "NDC");
  //pave2->AddText("#bf{ALICE preliminary}");
  pave2->AddText("pp, #sqrt{#it{s}} = 5.02 TeV");
  if (!correctByBR) pave2->AddText("no BR correction");
  pave2->SetFillStyle(0);
  pave2->SetBorderSize(0);
  pave2->SetTextFont(42);
  pave2->SetTextSize(0.05);
  pave2->SetTextAlign(12);

  pad1->cd();
  rebinned[0]->GetYaxis()->SetRangeUser(5e-6, 9.);
  rebinned[0]->Draw();
  Beautify(pi_syst, kRed, 1, 2, 25, 0.7);
  pi_syst->Scale(0.5);
  pi_syst->GetYaxis()->SetRangeUser(1e-6, 9.);
  pi_syst->Draw("same");
  f0_Tot[0]->Draw("same");
  pave->Draw();
  leg->Draw();
  
  pad2->cd();  
  ratio->GetYaxis()->SetRangeUser(0.0, 0.32);
  ratio->GetYaxis()->SetTitle("2f_{0}/(#pi^{+}+#pi^{-})");
  pave->DrawClone();

  ratio->Draw();
  pave2->Draw();
  cr->Print("ratioToPionVsPt.pdf");
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


/*****************************************************************/
/* LEVY-TSALLIS */
/*****************************************************************/
Double_t LevyTsallis_Func(const Double_t *x, const Double_t *p)
{
  /* dN/dpt */
  
  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t n = p[1];
  Double_t C = p[2];
  Double_t norm = p[3];

  Double_t part1 = (n - 1.) * (n - 2.);
  Double_t part2 = n * C * (n * C + mass * (n - 2.));
  Double_t part3 = part1 / part2;
  Double_t part4 = 1. + (mt - mass) / n / C;
  Double_t part5 = TMath::Power(part4, -n);
  return pt * norm * part3 * part5;
}

TF1 * LevyTsallis(const Char_t *name, Double_t mass, Double_t n, Double_t C, Double_t norm)
{
  TF1 *fLevyTsallis = new TF1(name, LevyTsallis_Func, 0., 15., 4);
  fLevyTsallis->SetParameters(mass, n, C, norm);
  fLevyTsallis->SetParNames("mass", "n", "C", "norm");
  fLevyTsallis->FixParameter(0, mass);
  fLevyTsallis->SetParLimits(1, 1.e-3, 1.e3);
  fLevyTsallis->SetParLimits(2, 1.e-3, 1.e3);
  fLevyTsallis->SetParLimits(3, 1.e-6, 1.e6);
  return fLevyTsallis;
}



void PlotRatioProtonToF0(Bool_t correctByBR)
{
  
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetLineWidth(1);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleOffset(0.8,"y");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  TGaxis::SetMaxDigits(4);
  
  TCanvas * cr = 0x0;
  TPad * pad1 = 0x0;
  TPad * pad2 = 0x0;
  cr = new TCanvas("cr","compare",800, 1000);
  cr->cd();
  pad1 = new TPad("pad1","This is pad1",0.001,0.5,0.999,0.999);
  pad1->SetFillColor(0);
  pad1->SetBorderMode(0);
  pad1->SetBorderSize(0);
  pad1->SetMargin(0.15,0.05,0.01,0.02);

  pad2 = new TPad("pad2","This is pad2",0.001,0.001, 0.999,0.5);
  pad2->SetFillColor(0);
  pad2->SetBorderMode(0);
  pad2->SetBorderSize(0);
  pad2->SetMargin(0.15,0.05,0.2,0.01);
  pad1->Draw();
  pad2->Draw();

  pad1->cd();
  pad1->SetLogy();

  //-------------
  // F0
  //-------------
  TFile* myfile=TFile::Open("prel_finalWsys.root");
  if (!myfile) {Printf("Invalid file."); return;}
  TH1F * f0_Tot[1];
  f0_Tot[0] = (TH1F*) myfile->Get("hCorrectedSpectrum_syst_stat");
  if (!correctByBR) f0_Tot[0]->Scale(0.46);
  
  TH1F * rebinned[1];
  rebinned[0] = (TH1F*) f0_Tot[0]->Clone("dummy");
  rebinned[0]->Reset("ICES");
  
  //-------------
  //PROTONS
  //-------------
  ///Users/fbellini/alice/pwglf-piKp5teV/SpectraAnalysisRun2/results/spectra/spectra-pag/Preliminaries/SQM2017/
  TString fileProName = "Spectra_ppLHC15n_Combined_Histograms.root";
  TFile *filePro = TFile::Open(fileProName.Data());
  if (!filePro) {Printf("Invalid file."); return;}
  TH1F *pro_syst = (TH1F*) ((TList*) filePro->Get("Summed_Proton_Sys"))->FindObject("hSpectraSummedProton_pp_Combined_MB");
  TH1F *pro_stat = (TH1F*) ((TList*) filePro->Get("Summed_Proton"))->FindObject("hSpectraSummedProton_pp_Combined_MB");
  
  if (!pro_stat) {Printf("missing plot proon data_stat"); return;}
  if (!pro_syst) {Printf("missing plot proon data_syst"); return;}
  
  TH1F * pro_Tot[1];
  pro_Tot[0] = (TH1F*) SumUncertInQuadrature(pro_stat, pro_syst);
  //pro_Tot[0]->Scale(0.5); // (pro+ + pro-) / 2
  
  TF1 * myLevyTsallis = LevyTsallis("levy", 0.938);
  RebinSpectrum(pro_Tot, rebinned, myLevyTsallis, 1, 1,  pro_Tot[0]);
  Beautify(rebinned[0], kBlue+1, 1, 2, 24, 1.0);
  Beautify(f0_Tot[0], kBlack, 1, 2, 20, 1.0);
  
  TH1F* ratio = (TH1F*) GetPlotRatio(rebinned[0],f0_Tot[0]);

  TString text="Uncertainties: #sqrt{stat.^{2} + sys.^{2}}";
  TPaveText * pave = new TPaveText(0.4,0.1,0.4,0.14,"NDC");
  pave->SetBorderSize(0);
  pave->SetFillStyle(0);
  pave->SetTextColor(kBlack);
  pave->SetTextFont(42);
  pave->SetTextSize(0.05);
  pave->InsertText(text.Data());

  TLegend * leg = new TLegend(0.4, 0.6,0.8,0.9, "pp, #sqrt{#it{s}} = 5.02 TeV");
  leg->AddEntry(pro_syst, "(p+#bar{p}) ALICE preliminary", "p");
  leg->AddEntry(rebinned[0], "(p+#bar{p}) rebinned", "p");
  leg->AddEntry(f0_Tot[0], "f_{0}(980), this analysis", "p");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  
  TPaveText * pave2 = new TPaveText(0.2, 0.8,0.5,0.9, "NDC");
  pave2->AddText("pp, #sqrt{#it{s}} = 5.02 TeV");
  if (!correctByBR) pave2->AddText("no BR correction");
  pave2->SetFillStyle(0);
  pave2->SetBorderSize(0);
  pave2->SetTextFont(42);
  pave2->SetTextSize(0.05);
  pave2->SetTextAlign(12);

  pad1->cd();
  rebinned[0]->GetYaxis()->SetRangeUser(5e-6, 9.);
  rebinned[0]->Draw();
  Beautify(pro_syst, kBlue+1, 1, 2, 25, 0.7);
  Beautify(rebinned[0], kBlack, 1, 2, 21, 0.7);
  Beautify(f0_Tot[0], kBlack, 1, 2, 20, 0.7);

  //pro_syst->Scale(0.5);
  pro_syst->GetYaxis()->SetRangeUser(1e-6, 9.);
  pro_syst->Draw("same");
  f0_Tot[0]->Draw("same");
  pave->Draw();
  leg->Draw();
  
  pad2->cd();
  Beautify(ratio, kBlack, 1, 2, 20, 0.7);
  ratio->GetYaxis()->SetRangeUser(0.0, 10.2);
  ratio->GetYaxis()->SetTitle("(p+#bar{p})/f_{0}");
  pave->DrawClone();

  ratio->Draw();
  pave2->Draw();
  cr->Print("ratioProtonToF0VsPt.pdf");

}




void PlotRatiof0ToPhi(Bool_t correctByBR)
{
  
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetLineWidth(1);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleOffset(0.8,"y");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  TGaxis::SetMaxDigits(4);
  
  TCanvas * cr = 0x0;
  TPad * pad1 = 0x0;
  TPad * pad2 = 0x0;
  cr = new TCanvas("cr","compare",800, 1000);
  cr->cd();
  pad1 = new TPad("pad1","This is pad1",0.001,0.5,0.999,0.999);
  pad1->SetFillColor(0);
  pad1->SetBorderMode(0);
  pad1->SetBorderSize(0);
  pad1->SetMargin(0.15,0.05,0.01,0.02);

  pad2 = new TPad("pad2","This is pad2",0.001,0.001, 0.999,0.5);
  pad2->SetFillColor(0);
  pad2->SetBorderMode(0);
  pad2->SetBorderSize(0);
  pad2->SetMargin(0.15,0.05,0.2,0.01);
  pad1->Draw();
  pad2->Draw();

  pad1->cd();
  pad1->SetLogy();

  //-------------
  // F0
  //-------------
  TFile* myfile=TFile::Open("prel_finalWsys.root");
  if (!myfile) {Printf("Invalid file."); return;}
  TH1F * f0_Tot[1];
  f0_Tot[0] = (TH1F*) myfile->Get("hCorrectedSpectrum_syst_stat");
  if (!correctByBR) f0_Tot[0]->Scale(0.46);
  
  TH1F * rebinned[1];
  rebinned[0] = (TH1F*) f0_Tot[0]->Clone("dummy");
  rebinned[0]->Reset("ICES");
  
  //-------------
  //phi
  //-------------
  TString filePhiName = "phi_pp5_prelim.root";
  TFile *filePhi = TFile::Open(filePhiName.Data());
  if (!filePhi) {Printf("Invalid file."); return;}
  TH1F *phi_syst = (TH1F*) filePhi->Get("h2");
  TH1F *phi_stat = (TH1F*) filePhi->Get("h1");
  
  if (!phi_stat) {Printf("missing plot phi data_stat"); return;}
  if (!phi_syst) {Printf("missing plot phi data_syst"); return;}
  
  TH1F * phi_Tot[1];
  phi_Tot[0] = (TH1F*) filePhi->Get("h3");
  
  TF1 * myLevyTsallis = LevyTsallis("levy", 1.01965);
  RebinSpectrum(phi_Tot, rebinned, myLevyTsallis, 1, 1,  phi_Tot[0]);
  Beautify(rebinned[0], kGreen+1, 1, 2, 28, 1.0);
  Beautify(f0_Tot[0], kBlack, 1, 2, 20, 1.0);
  myLevyTsallis->SetLineStyle(2);
  myLevyTsallis->SetLineColor(kBlack);
    
  TH1F* ratio = (TH1F*) GetPlotRatio(f0_Tot[0],rebinned[0]);

  TString text="Uncertainties: #sqrt{stat.^{2} + sys.^{2}}";
  TPaveText * pave = new TPaveText(0.4,0.1,0.4,0.14,"NDC");
  pave->SetBorderSize(0);
  pave->SetFillStyle(0);
  pave->SetTextColor(kBlack);
  pave->SetTextFont(42);
  pave->SetTextSize(0.05);
  pave->InsertText(text.Data());

  TLegend * leg = new TLegend(0.4, 0.6,0.8,0.9, "pp, #sqrt{#it{s}} = 5.02 TeV");
  leg->AddEntry(phi_syst, "#phi, ALICE preliminary", "p");
  leg->AddEntry(rebinned[0], "#phi rebinned", "p");
  leg->AddEntry(f0_Tot[0], "f_{0}/#phi, this analysis", "p");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  
  TPaveText * pave2 = new TPaveText(0.2, 0.8,0.5,0.9, "NDC");
  //pave2->AddText("#bf{ALICE preliminary}");
  pave2->AddText("pp, #sqrt{#it{s}} = 5.02 TeV");
  if (!correctByBR) pave2->AddText("no BR correction");
  pave2->SetFillStyle(0);
  pave2->SetBorderSize(0);
  pave2->SetTextFont(42);
  pave2->SetTextSize(0.05);
  pave2->SetTextAlign(12);

  pad1->cd();
  rebinned[0]->GetYaxis()->SetRangeUser(5e-6, 9.);
  Beautify(phi_syst, kGreen+1, 1, 2, 25, 1.);
  Beautify(rebinned[0], kGreen+3, 1, 2, 28, 1.3);
  Beautify(f0_Tot[0], kBlack, 1, 2, 20, 1.);
  phi_syst->GetYaxis()->SetRangeUser(1e-6, 9.);

  rebinned[0]->Draw();
  phi_syst->Draw("same");
  myLevyTsallis->Draw("same");
  f0_Tot[0]->Draw("same");
  pave->Draw();
  leg->Draw();
  
  pad2->cd();  
  ratio->GetYaxis()->SetRangeUser(0.0, 2.5);
  Beautify(ratio, kBlack, 1, 2, 20, 1.);
  ratio->GetYaxis()->SetTitle("f_{0}/#phi");
  pave->DrawClone();

  ratio->Draw();
  pave2->Draw();
  cr->Print("ratioF0ToPhiVsPt.pdf");

}
