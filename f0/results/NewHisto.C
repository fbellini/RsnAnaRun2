#include "SetStyle.C"
#include "RebinSpectrum.C"

TH1F * SumErrorsInQuadrature(TH1F * hstat, TH1F* hsys);

void NewHisto(
  TString f0Spectrum = "prel_finalWsys.root",
  TString phiSpectrum = "phi_pp5_prelim.root",
  TString protonsSpectrum = "Spectra_ppLHC15n_Combined_Histograms.root",
  TString pionsSpectrum = "Spectra_ppLHC15n_Combined_Histograms.root",
  Double_t phiValue = 0.0257552,
  Double_t phiError = 0.0001,
  Double_t protonsValue = 0.0257552,
  Double_t protonsError = 0.0001,
  Double_t pionsValue = 0.0257552,
  Double_t pionsError = 0.0001
)
{

  TGaxis::SetMaxDigits(3);

  TFile* fout = TFile::Open("newHistos.root", "RECREATE");

  /* ---------------------------------------- f0(980) ---------------------------------------- */
  TFile *f0FSpec = TFile::Open(Form("%s", f0Spectrum.Data()));
  if (!f0FSpec) {Printf("Invalid file."); return;}
  TH1F * f0_Stat_Syst[1];
  f0_Stat_Syst[0] = (TH1F*)f0FSpec->Get("hCorrectedSpectrum_syst_stat");


  /* ---------------------------------------- phi(1020) ---------------------------------------- */
  TFile *phiFSpec = TFile::Open(Form("%s", phiSpectrum.Data()));
  if (!phiFSpec) {Printf("Invalid file."); return;}
  TH1F *hPhi_Stat = (TH1F*) phiFSpec->Get("h1");
  TH1F *hPhi_Syst = (TH1F*) phiFSpec->Get("h2");
  TH1F *hPhi_Stat_Syst = (TH1F*) phiFSpec->Get("h3");

  if (!hPhi_Stat) {Printf("missing plot phi data_stat"); return;}
  if (!hPhi_Syst) {Printf("missing plot phi data_syst"); return;}
  if (!hPhi_Stat_Syst) {Printf("missing plot phi data_syst"); return;}

  Double_t binsPhi [] = {0., 0.4, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.1, 2.4, 2.7, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 9.0, 12.0, 16.0, 20.0};
  Int_t nBinsPhi = sizeof(binsPhi) / sizeof(binsPhi[0]) - 1;
  TH1F* hNewSpectrumStatSystPhi = new TH1F("hNewSpectrumStatSystPhi", "hNewSpectrumStatSystPhi", nBinsPhi, binsPhi);
  for (Int_t ibin=0; ibin<26; ibin++){
    if (ibin==0){
      hNewSpectrumStatSystPhi->SetBinContent(ibin+1, phiValue);
      hNewSpectrumStatSystPhi->SetBinError(ibin+1, phiError);
    }
    else{
      hNewSpectrumStatSystPhi->SetBinContent(ibin+1, hPhi_Stat_Syst->GetBinContent(ibin));
      hNewSpectrumStatSystPhi->SetBinError(ibin+1, hPhi_Stat_Syst->GetBinError(ibin));
    }
  }


  /* ---------------------------------------- protons ---------------------------------------- */
  TFile *protonsFSpec = TFile::Open(Form("%s", protonsSpectrum.Data()));
  if (!protonsFSpec) {Printf("Invalid file."); return;}
  TH1F *pro_syst = (TH1F*) ((TList*) protonsFSpec->Get("Summed_Proton_Sys"))->FindObject("hSpectraSummedProton_pp_Combined_MB");
  TH1F *pro_stat = (TH1F*) ((TList*) protonsFSpec->Get("Summed_Proton"))->FindObject("hSpectraSummedProton_pp_Combined_MB");

  if (!pro_stat) {Printf("missing plot proon data_stat"); return;}
  if (!pro_syst) {Printf("missing plot proon data_syst"); return;}

  TH1F * pro_Tot[1];
  pro_Tot[0] = (TH1F*) SumErrorsInQuadrature(pro_stat, pro_syst);
  pro_Tot[0]->Scale(0.5); // (pi+ + pi-) / 2

  Double_t binsProtons [] = {0., 0.4, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.1, 2.4, 2.7, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 9.0, 12.0, 16.0, 20.0};
  Int_t nBinsProtons = sizeof(binsProtons) / sizeof(binsProtons[0]) - 1;
  TH1F* hNewSpectrumStatSystProtons = new TH1F("hNewSpectrumStatSystProtons", "hNewSpectrumStatSystProtons", nBinsProtons, binsProtons);
  for (Int_t ibin2=0; ibin2<26; ibin2++){
    if (ibin2==0){
      hNewSpectrumStatSystProtons->SetBinContent(ibin2+1, 0.0257552);
      hNewSpectrumStatSystProtons->SetBinError(ibin2+1, 0.0001);
    }
    else{
    hNewSpectrumStatSystProtons->SetBinContent(ibin2+1, pro_Tot[0]->GetBinContent(ibin2));
    hNewSpectrumStatSystProtons->SetBinError(ibin2+1, pro_Tot[0]->GetBinError(ibin2));
    }
  }


  /* ---------------------------------------- pions ---------------------------------------- */
  TFile *pionsFSpec = TFile::Open(Form("%s", pionsSpectrum.Data()));
  if (!pionsFSpec) {Printf("Invalid file."); return;}
  TH1F *pi_syst = (TH1F*) ((TList*) pionsFSpec->Get("Summed_Pion_Sys"))->FindObject("hSpectraSummedPion_pp_Combined_MB");
  TH1F *pi_stat = (TH1F*) ((TList*) pionsFSpec->Get("Summed_Pion"))->FindObject("hSpectraSummedPion_pp_Combined_MB");

  if (!pi_stat) {Printf("missing plot pion data_stat"); return;}
  if (!pi_syst) {Printf("missing plot pion data_syst"); return;}

  TH1F * pi_Tot[1];
  pi_Tot[0] = (TH1F*) SumErrorsInQuadrature(pi_stat, pi_syst);
  pi_Tot[0]->Scale(0.5); // (pi+ + pi-) / 2

  Double_t binsPions [] = {0., 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0};
  Int_t nBinsPions = sizeof(binsPions) / sizeof(binsPions[0]) - 1;
  TH1F* hNewSpectrumStatSystPions = new TH1F("hNewSpectrumStatSystPions", "hNewSpectrumStatSystPions", nBinsPions, binsPions);
  for (Int_t ibin3=0; ibin3<54; ibin3++){
    if (ibin3==0){
      hNewSpectrumStatSystPions->SetBinContent(ibin3+1, 0.0257552);
      hNewSpectrumStatSystPions->SetBinError(ibin3+1, 0.0001);
    }
    else{
    hNewSpectrumStatSystPions->SetBinContent(ibin3, pi_Tot[0]->GetBinContent(ibin3));
    hNewSpectrumStatSystPions->SetBinError(ibin3, pi_Tot[0]->GetBinError(ibin3));
    }
  }

/* ---------------------------------------- new histos ---------------------------------------- */
  /*TH1F * hNewSpectrum[1];
  hNewSpectrum[0] = (TH1F*) f0_Stat_Syst[0]->Clone("dummy");
  hNewSpectrum[0]->Reset("ICES");

  Double_t minBin, maxBin;
  minBin=hNewSpectrum[0]->GetXaxis()->GetBinLowEdge(1);
  maxBin=hNewSpectrum[0]->GetXaxis()->GetBinUpEdge(1);*/




  //printf("%f, %f\n", minBin, maxBin);


  TCanvas* c1 = new TCanvas("c1", "c1", 600, 600);
  c1->cd();
  hPhi_Stat_Syst->Draw();
  hNewSpectrumStatSystPhi->Draw("same");

  TCanvas* c2 = new TCanvas("c2", "c2", 600, 600);
  c2->cd();
  pro_Tot[0]->Draw();
  hNewSpectrumStatSystProtons->Draw("same");

  TCanvas* c3 = new TCanvas("c3", "c3", 600, 600);
  c3->cd();
  pi_Tot[0]->Draw();
  hNewSpectrumStatSystPions->Draw("same");


  fout->cd();
  hNewSpectrumStatSystPhi->Write("hPhi_StatSyst");
  hNewSpectrumStatSystProtons->Write("hProtons_StatSyst");
  hNewSpectrumStatSystPions->Write("hPions_StatSyst");
  fout->Close();
  return;
}


TH1F * SumErrorsInQuadrature(TH1F * hstat, TH1F* hsys)
{
  if (!hstat || !hsys) return 0;
  TH1F * hnew = (TH1F *) hstat->Clone();

  for (int i = 1; i<hnew->GetNbinsX()+1; i++){
    Float_t errInQuad =  TMath::Sqrt(TMath::Power(hstat->GetBinError(i), 2.0) + TMath::Power(hsys->GetBinError(i), 2.0));
    hnew->SetBinError(i, errInQuad);
  }
  return hnew;
}
