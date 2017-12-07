//////////////////////////////////////////////////////////////////////////////////////////////////////
//effciencies ratios - 30.11.2017                                                                   //
//                                                                                                  //
//before using the macro, you need HistoMakeUp.C and SetStyle.C in the same directory. you also     //
//need the following lines in a rootlogon.C                                                         //
//{                                                                                                 //
//gROOT->LoadMacro("HistoMakeUp.C+g");                                                              //
//gROOT->LoadMacro("SetStyle.C+g");                                                                 //
//SetStyle();                                                                                       //
//gROOT->LoadMacro("f0_eff.C+g");                                                                   //
//}                                                                                                 //
//                                                                                                  //
//////////////////////////////////////////////////////////////////////////////////////////////////////

#include "HistoMakeUp.C"
#include "SetStyle.C"
#include "TCanvas.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TSystem.h"
#include "f0Params.h"

TH1D* ratio_ef(TH1D* h1, TH1D* h2);
TH1D* divis(TH1D* h1, TH1D* h2);

void f0_eff()
{
  TGaxis::SetMaxDigits(3);

  TFile* file = TFile::Open("AnalysisResults.root");
  TDirectory* directory1 = (TDirectory*)file->Get("F0_bkg");
  directory1->cd();
  TDirectory* directory2 = (TDirectory*)directory1->Get("MyOutputContainer");
  directory2->cd();

  TCanvas* canvas[10];
  TFile* fout = TFile::Open("efficiencies.root", "RECREATE");

  enum PartSpec { kF0 = 0,
    kOmega = 1,
    kRho = 2,
    kEta = 3,
    kEtaPr = 4,
    kF1
    = 5,
    kF2 = 6,
    kKStar = 7,
    kK0s = 8,
    kPhi = 9,
    kAll = 10 };
  const Char_t* particles[kAll] = { "f0", "omega", "rho", "eta", "etaPr", "f1", "f2", "kStar", "k0s", "phi" };
  enum { kPtBins = 11 };
  const Double_t pT[kPtBins] = { 0.5, 1., 1.5, 2., 2.5, 3., 4., 5., 7., 9., 11. };

  TH2F* hPtGen[kAll] = { 0x0 };
  TH2F* hPtReco[kAll] = { 0x0 };
  TH1D* hPtGen1d[kAll] = { 0x0 };
  TH1D* hPtReco1d[kAll] = { 0x0 };
  TH1D* hEff[kAll] = { 0x0 };

  for (Int_t i = 0; i < kAll; i++) {
    hPtGen[i] = (TH2F*)directory2->FindObject(Form("fHistPtGen_%s", particles[i]));
    hPtReco[i] = (TH2F*)directory2->FindObject(Form("fHistPtReco_%s", particles[i]));
    hPtGen1d[i] = (TH1D*)hPtGen[i]->ProjectionY(Form("hPtGen1d%s", particles[i]), 0., -1.);
    hPtReco1d[i] = (TH1D*)hPtReco[i]->ProjectionY(Form("hPtReco1d%s", particles[i]), 0., -1.);
    hEff[i] = ratio_ef(hPtReco1d[i], hPtGen1d[i]);
    canvas[i] = new TCanvas(Form("c%d", i), Form("%s efficiency", particles[i]), 1500., 500.);
    canvas[i]->Divide(3, 1);
    canvas[i]->cd(1);
    hPtGen1d[i]->Draw("PE");
    canvas[i]->cd(2);
    hPtReco1d[i]->Draw("PE");
    canvas[i]->cd(3);
    hEff[i]->Draw("PE");
  }

  return;
}

TH1D* ratio_ef(TH1D* h1, TH1D* h2)
{
  if (!h1 || !h2)
    return NULL;
  TH1D* hResult = (TH1D*)h1->Clone(Form("%s_eff", h1->GetName()));
  hResult->Reset("ICES");
  //Printf("%s", hResult->GetName());
  const Double_t momT[12] = { 0., 0.5, 1., 1.5, 2., 2.5, 3., 4., 5., 7., 9., 11. };
  Double_t cont1 = 0.0, cont2 = 0.0, error1 = 0.0, error2 = 0.0, error3 = 0.0, ratio_eff = 0.0;
    for(Int_t index=0; index<12; index++){
    h1->RebinX(h1->GetXaxis()->FindBin(momT[index+1]) - h1->GetXaxis()->FindBin(momT[index]));
    h2->RebinX(h2->GetXaxis()->FindBin(momT[index+1]) - h2->GetXaxis()->FindBin(momT[index]));
    for (Int_t j = 1; j < h1->GetNbinsX() + 1; j++) {
    cont1 = h1->GetBinContent(j);
    cont2 = h2->GetBinContent(j);
    error1 = h1->GetBinError(j);
    error2 = h2->GetBinError(j);
    ratio_eff = cont1/cont2;
    error3 = pow((pow(error1, 2) + pow(error2, 2)), 1. / 2);
    hResult->SetBinContent(j, ratio_eff);
    hResult->SetBinError(j, error3);
  }
}

  return hResult;
}

TH1D* divis(TH1D* h1, TH1D* h2)
{
  TH1D* hResult = NULL;
  //returns division of two histograms
  if (!h1 || !h2)
    return NULL;
  Int_t nBin1 = h1->GetXaxis()->GetNbins();
  Int_t nBin2 = h2->GetXaxis()->GetNbins();
  if (nBin1 == nBin2) {
    hResult = (TH1D*)h1->Clone();
    hResult->Divide(h2);
  } else {
    Printf("Histograms do not have the same binning.  nBin1 = %d, nBin2 = %d", nBin1, nBin2);
    return NULL;
  }
  return hResult;
}
