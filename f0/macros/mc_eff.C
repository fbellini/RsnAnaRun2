//////////////////////////////////////////////////////////////////////////////////////////////////////
//MC efficiencies - 31.10.2017                                                                      //
//                                                                                                  //
//before using the macro, you need HistoMakeUp.C and SetStyle.C in the same directory. you also     //
//need the following lines in a rootlogon.C                                                         //
//{                                                                                                 //
//gROOT->LoadMacro("HistoMakeUp.C+g");                                                              //
//gROOT->LoadMacro("SetStyle.C+g");                                                                 //
//SetStyle();                                                                                       //
//gROOT->LoadMacro("mc_eff.C+g");                                                                   //
//}                                                                                                 //
//                                                                                                  //                                                                  //
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
TH1D* div(TH1D* h1, TH1D* h2);


void mc_eff(Int_t rebinVar = 5 /*rebinning factor*/,
                  Double_t chosenPt = -1 /*chosen pT allows to chose the pT bin, if chosenPt = -1 you cover all pT bins
                  list of matches
                  0.5 < pT < 1.0 -> chosenPt=0; 1.0 < pT < 1.5 -> chosenPt=1; 1.5 < pT < 2.0 -> chosenPt=2;
                  2.0 < pT < 2.5 -> chosenPt=3; 2.5 < pT < 3.0 -> chosenPt=4; 3.0 < pT < 4.0 -> chosenPt=5;
                  4.0 < pT < 5.0 -> chosenPt=6; 5.0 < pT < 7.0 -> chosenPt=7; 7.0 < pT < 9.0 -> chosenPt=8;
                  9.0 < pT < 11.0 -> chosenPt=9 */
)
{
  TGaxis::SetMaxDigits(3);

  /* getting histos from root file */
  TFile* file = TFile::Open("AnalysisResultsMC.root");
  TList* list = (TList*)file->Get("RsnOut_f0");
  //list->ls();
  TH2F* hTrueF0 = (TH2F*)list->FindObject("RsnTaskF0_truef0_f0");
  TH2F* hMotherF0 = (TH2F*)list->FindObject("RsnTaskF0_motherf0_f0");


  /* binning: 500 MeV bins up to 3 GeV/c, then 1 GeV bins up to 5 GeV/c, then 2 GeV up to 11 GeV/c */
  Double_t pT[11] = {0.5, 1., 1.5, 2., 2.5, 3., 4., 5., 7., 9., 11.};

  TH1D* hTrueF0px[11];
  TH1D* hMotherF0px[11];
  TH1D* hRatio[11];

  for (Int_t i = 0; i < 10; i++) {
    hTrueF0px[i] = 0x0;
    hMotherF0px[i] = 0x0;
    hRatio[i] = 0x0;
  }

  TFile* fout = TFile::Open("mc_efficiencies.root", "RECREATE");

  TCanvas* canvas[10];

  Int_t nBin1 = hTrueF0->GetXaxis()->GetNbins();
  Int_t nBin2 = hMotherF0->GetXaxis()->GetNbins();

  if (nBin1 == nBin2) {
    for (Int_t ibin = 0; ibin < 10; ibin++) {
      if (chosenPt >= 0 && chosenPt != ibin)
        continue;

      Int_t iMinBinPt = hTrueF0->GetYaxis()->FindBin(pT[ibin]);
      Int_t iMaxBinPt = hTrueF0->GetYaxis()->FindBin(pT[ibin + 1]);
      //printf("min: %d - max:%d; pT min: %5.2f - pT max: %05.2f\n", iMinBinPt, iMaxBinPt, pT[ibin], pT[ibin+1]);

      hTrueF0px[ibin] = (TH1D*)hTrueF0->ProjectionX(Form("hUSP_pT_%2.1f-%2.1f", pT[ibin], pT[ibin + 1]), iMinBinPt, iMaxBinPt);
      hTrueF0px[ibin]->SetTitle(Form("True reconstructed pairs %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin + 1]));
      hTrueF0px[ibin]->RebinX(rebinVar);

      hMotherF0px[ibin] = (TH1D*)hMotherF0->ProjectionX(Form("hLikeMM_pT_%2.1f-%2.1f", pT[ibin], pT[ibin + 1]), iMinBinPt, iMaxBinPt);
      hMotherF0px[ibin]->SetTitle(Form("F0 mothers - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin + 1]));
      hMotherF0px[ibin]->RebinX(rebinVar);

      /* efficiencies */
      hRatio[ibin] = div(hTrueF0px[ibin], hMotherF0px[ibin]);
      hRatio[ibin]->SetTitle(Form("Efficiency - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin + 1]));
      //hRatio[ibin]->GetYaxis()->SetRangeUser(0.,1.);


      /* Histo maquillage */
      HistoMakeUp(hTrueF0px[ibin], kRed, 20, "#it{M}_{#pi#pi}(GeV/#it{c^{2}})", "pairs");
      HistoMakeUp(hMotherF0px[ibin], kBlue, 20, "#it{M}_{#pi#pi}(GeV/#it{c^{2}})", "pairs");
      HistoMakeUp(hRatio[ibin], kBlack, 20, "#it{M}_{#pi#pi}(GeV/#it{c^{2}})", "pairs");

      TLegend* legend1 = new TLegend(0.75, 0.75, 0.9, 0.9);
      legend1->AddEntry(hTrueF0px[ibin], "True", "lpf");
      legend1->AddEntry(hMotherF0px[ibin], "Mothers", "lpf");

      /* drawing and printing output */
      if (chosenPt == -1) {
        canvas[ibin] = new TCanvas(Form("c%d", ibin), Form("canvas %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin + 1]), 1600, 800);
        canvas[ibin]->Divide(2, 1);
        canvas[ibin]->cd(1);
        hTrueF0px[ibin]->Draw("e");
        hMotherF0px[ibin]->Draw("e same");
        legend1->Draw();
        canvas[ibin]->cd(2);
        hRatio[ibin]->Draw("e");
    } else if (chosenPt == ibin) {
      canvas[ibin] = new TCanvas(Form("c%d", ibin), Form("canvas %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin + 1]), 1600, 800);
      canvas[ibin]->Divide(2, 1);
      canvas[ibin]->cd(1);
      hTrueF0px[ibin]->Draw("e");
      hMotherF0px[ibin]->Draw("e same");
      legend1->Draw();
      canvas[ibin]->cd(2);
      hRatio[ibin]->Draw("e");

      /* writing output on file.root*/
      fout->cd();
      hTrueF0px[ibin]->Write("true_reco_F0");
      hMotherF0px[ibin]->Write("F0_mothers");
      hRatio[ibin]->Write("efficiency");
    }
  }
  } else {
    Printf("Histograms do not have the same binning.  nBin1 = %d, nBin2 = %d", nBin1, nBin2);
  }
    //fout->Close();
    return;
}

TH1D* div(TH1D* h1, TH1D* h2)
{
  //returns 2*(geometrical mean) of two histograms
  if (!h1 || !h2)
    return NULL;
  TH1D* hRatio = (TH1D*)h1->Clone();
  hRatio->Reset("ICES");
  Double_t cont1 = 0.0, cont2 = 0.0, error1 = 0.0, error2 = 0.0, error3 = 0.0, ratio = 0.0;
  for (Int_t j = 1; j < h1->GetNbinsX() + 1; j++) {
    cont1 = h1->GetBinContent(j);
    cont2 = h2->GetBinContent(j);
    error1 = h1->GetBinError(j);
    error2 = h2->GetBinError(j);
    ratio = (cont1/cont2);
    Printf("Ratio=%f", ratio);
    error3 = (error1/error2);
    hRatio->SetBinContent(j, ratio);
    hRatio->SetBinError(j, error3);
  }
  return hRatio;
}
