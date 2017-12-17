//////////////////////////////////////////////////////////////////////////////////////////////////////
//bg subtraction - 29.08.2017                                                                       //
//                                                                                                  //
//before using the macro, you need HistoMakeUp.C and SetStyle.C in the same directory. you also     //
//need the following lines in a rootlogon.C                                                         //
//{                                                                                                 //
//gROOT->LoadMacro("HistoMakeUp.C+g");                                                              //
//gROOT->LoadMacro("SetStyle.C+g");                                                                 //
//SetStyle();                                                                                       //
//gROOT->LoadMacro("f0_bkg.C+g");                                                                   //
//}                                                                                                 //
//                                                                                                  //
//AnalysisResults.root --> data                                                                     //
//AnalysisResultsMC.root --> MC                                                                     //
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
TH1D* subtrac(TH1D* h1, TH1D* h2);
TH1D* divis(TH1D* h1, TH1D* h2);
TH1D* summa(TH1D* h1, TH1D* h2);
TH1D* geoM(TH1D* h1, TH1D* h2);
TH1D* normM(TH1D* h1, TH1D* h2, Int_t method);

void f0_bkg(Int_t rebinVar = 2 /*rebinning factor*/,
                  Double_t chosenPt = -1 /*chosenPt allows to chose the pT bin, if chosenPt = -1 you cover all pT bins
                  list of matches
                  0.5 < pT < 1.0 -> chosenPt=0; 1.0 < pT < 1.5 -> chosenPt=1; 1.5 < pT < 2.0 -> chosenPt=2;
                  2.0 < pT < 2.5 -> chosenPt=3; 2.5 < pT < 3.0 -> chosenPt=4; 3.0 < pT < 4.0 -> chosenPt=5;
                  4.0 < pT < 5.0 -> chosenPt=6; 5.0 < pT < 7.0 -> chosenPt=7; 7.0 < pT < 9.0 -> chosenPt=8;
                  9.0 < pT < 11.0 -> chosenPt=9 */,
                  Bool_t isMC = kFALSE
)
{
  TGaxis::SetMaxDigits(3);

  /* getting histos from root file */
  TFile* file = TFile::Open(Form("AnalysisResults%s.root", isMC ? "MC" : ""));
  TList* list = (TList*)file->Get("RsnOut_f0");
  //list->ls();
  TH2F* hUSP = (TH2F*)list->FindObject("RsnTaskF0_f0_UnlikePM"); //same event unlike sign pairs
  TH2F* hLikePP = (TH2F*)list->FindObject("RsnTaskF0_f0_LikePP"); //like sign pairs: (pi+)+(pi+)
  TH2F* hLikeMM = (TH2F*)list->FindObject("RsnTaskF0_f0_LikeMM"); //like sign pairs: (pi-)+(pi-)
  TH2F* hMELPP = (TH2F*)list->FindObject("RsnTaskF0_f0_MixingPP"); //(pi+)+(pi+) from 5 similar events
  TH2F* hMELMM = (TH2F*)list->FindObject("RsnTaskF0_f0_MixingMM"); //(pi-)+(pi-) from 5 similar events
  TH2F* hMEB = (TH2F*)list->FindObject("RsnTaskF0_f0_MixingPM"); //(pi+)+(pi-) from 5 similar events

  /* projecting histos in 2d (losing information on multiplicity)
  TH2D * hUSP = hOrUnlikePM->Project3D("xy");
  TH2D * hLikeMM = hOrLikeMM->Project3D("xy");
  TH2D * hLikePP = hOrLikePP->Project3D("xy");
  TH2D * hMEB = hOrMixingPM->Project3D("xy");
  TH2D * hMELMM = hOrMixingMM->Project3D("xy");
  TH2D * hMELPP = hOrMixingPP->Project3D("xy");*/

  /* binning: 500 MeV bins up to 3 GeV/c, then 1 GeV bins up to 5 GeV/c, then 2 GeV up to 11 GeV/c */
  Double_t pT[11] = {0.5, 1., 1.5, 2., 2.5, 3., 4., 5., 7., 9., 11.};

  Int_t   npT  = sizeof(pT) / sizeof(pT[0]) - 1;
  TAxis *ptBins = new TAxis(npT, pT);

  TH1D* hLikePPpy[10];
  TH1D* hLikeMMpy[10];
  TH1D* hMEBpy[10];
  TH1D* hUSPpy[10];
  TH1D* hLSBRatio[10];
  TH1D* hGeoMean[10];
  TH1D* hSum[10];
  TH1D* hnMEB[10];
  TH1D* hMELRatio[10];
  TH1D* hMELMEBRatio[10];
  TH1D* hUSPminusLSB1[10];
  TH1D* hUSPminusLSB2[10];
  TH1D* hUSPminusMEB[10];
  TH1D* hUSPminusMEL[10];
  TH1D* hMELMMpy[10];
  TH1D* hMELPPpy[10];
  TH1D* hMEL[10];
  TH1D* hnMEL[10];
  TH1D* hLSBnMELRatio[10];
  TH1D* hLSBnMEBRatio[10];

  for (Int_t i = 0; i < 10; i++) {
    hLikePPpy[i] = 0x0;
    hLikeMMpy[i] = 0x0;
    hMEBpy[i] = 0x0;
    hUSPpy[i] = 0x0;
    hLSBRatio[i] = 0x0;
    hGeoMean[i] = 0x0;
    hSum[i] = 0x0;
    hnMEB[i] = 0x0;
    hMELRatio[i] = 0x0;
    hMELMEBRatio[i] = 0x0;
    hUSPminusLSB1[i] = 0x0;
    hUSPminusLSB2[i] = 0x0;
    hUSPminusMEB[i] = 0x0;
    hUSPminusMEL[i] = 0x0;
    hMELMMpy[i] = 0x0;
    hMELPPpy[i] = 0x0;
    hMEL[i] = 0x0;
    hnMEL[i] = 0x0;
    hLSBnMELRatio[i] = 0x0;
    hLSBnMEBRatio[i] = 0x0;
  }

  TFile* fout = TFile::Open("bgSubtraction.root", "RECREATE");

  /* canvas to print histograms for each bin in pT*/
  TCanvas* canvas[10];

  /*  bg subtraction for each bin in pT*/
  for (Int_t ibin = 0; ibin < 10; ibin++) {
    if (chosenPt >= 0 && chosenPt != ibin)
      continue;
    Int_t iMinBinPt = hUSP->GetYaxis()->FindBin(pT[ibin]);
    Int_t iMaxBinPt = hUSP->GetYaxis()->FindBin(pT[ibin + 1]);
    //printf("min: %d - max:%d; pT min: %5.2f - pT max: %05.2f\n", iMinBinPt, iMaxBinPt, pT[ibin], pT[ibin+1]);

    /*  projecting histos in 1d to get invariant mass spectra */
    hUSPpy[ibin] = (TH1D*)hUSP->ProjectionX(Form("hUSP_pT_%2.1f-%2.1f", pT[ibin], pT[ibin + 1]), iMinBinPt, iMaxBinPt);
    hUSPpy[ibin]->SetTitle(Form("USP %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin + 1]));
    hUSPpy[ibin]->RebinX(rebinVar);

    hLikeMMpy[ibin] = (TH1D*)hLikeMM->ProjectionX(Form("hLikeMM_pT_%2.1f-%2.1f", pT[ibin], pT[ibin + 1]), iMinBinPt, iMaxBinPt);
    hLikeMMpy[ibin]->SetTitle(Form("LSB (--) - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin + 1]));
    hLikeMMpy[ibin]->RebinX(rebinVar);

    hLikePPpy[ibin] = (TH1D*)hLikePP->ProjectionX(Form("hLikePP_pT_%2.1f-%2.1f", pT[ibin], pT[ibin + 1]), iMinBinPt, iMaxBinPt);
    hLikePPpy[ibin]->SetTitle(Form("LSB (++) - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin + 1]));
    hLikePPpy[ibin]->RebinX(rebinVar);

    hMEBpy[ibin] = (TH1D*)hMEB->ProjectionX(Form("hMEB_pT_%2.1f-%2.1f", pT[ibin], pT[ibin + 1]), iMinBinPt, iMaxBinPt);
    hMEBpy[ibin]->SetTitle(Form("MEB - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin + 1]));
    hMEBpy[ibin]->RebinX(rebinVar);

    hMELMMpy[ibin] = (TH1D*)hMELMM->ProjectionX(Form("hMELMM_pT_%2.1f-%2.1f", pT[ibin], pT[ibin + 1]), iMinBinPt, iMaxBinPt);
    hMELMMpy[ibin]->SetTitle(Form("MEL (--) - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin + 1]));
    hMELMMpy[ibin]->RebinX(rebinVar);

    hMELPPpy[ibin] = (TH1D*)hMELPP->ProjectionX(Form("hMELPP_pT_%2.1f-%2.1f", pT[ibin], pT[ibin + 1]), iMinBinPt, iMaxBinPt);
    hMELPPpy[ibin]->SetTitle(Form("MEL (++) - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin + 1]));
    hMELPPpy[ibin]->RebinX(rebinVar);

    /* Bg estimation */
    hGeoMean[ibin] = geoM(hLikePPpy[ibin], hLikeMMpy[ibin]); //LSB calculated with geometric mean
    hGeoMean[ibin]->SetTitle(Form("LSB geometric mean - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin + 1]));
    //hGeoMean[ibin]->GetYaxis()->SetRangeUser(0.,0.7e6);

    hSum[ibin] = summa(hLikePPpy[ibin], hLikeMMpy[ibin]); //LSB calculated with arithmetic mean
    hSum[ibin]->SetTitle(Form("LSB sum - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin + 1]));

    hLSBRatio[ibin] = divis(hLikeMMpy[ibin], hLikePPpy[ibin]);
    hLSBRatio[ibin]->SetTitle(Form("Like Sign Ratios - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin + 1]));
    //hLSBRatio[ibin]->GetYaxis()->SetRangeUser(0.8,1.2);

    hnMEB[ibin] = normM(hMEBpy[ibin], hUSPpy[ibin], 2);
    hnMEB[ibin]->SetTitle(Form("nMEB - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin + 1]));

    hMELRatio[ibin] = divis(hMELMMpy[ibin], hMELPPpy[ibin]); // MEL -- / MEL ++
    hMELRatio[ibin]->SetTitle(Form("MEL (--)/ MEL (++) - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin + 1]));

    hMEL[ibin] = geoM(hMELPPpy[ibin], hMELMMpy[ibin]); // MEL calculated with geometric mean
    hMEL[ibin]->SetTitle(Form("MEL geoMean - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin + 1]));

    hnMEL[ibin] = normM(hMEL[ibin], hUSPpy[ibin], 2); // MEL
    hnMEL[ibin]->SetTitle(Form("nMEL - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin + 1]));

    hMELMEBRatio[ibin] = divis(hMEL[ibin], hMEBpy[ibin]); // MEL / MEB
    hMELMEBRatio[ibin]->SetTitle(Form("MEL / MEB - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin + 1]));

    hLSBnMEBRatio[ibin] = divis(hGeoMean[ibin], hnMEB[ibin]); // LSB/nMEB
    hLSBnMEBRatio[ibin]->SetTitle(Form("LSB / nMEB - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin + 1]));

    hLSBnMELRatio[ibin] = divis(hGeoMean[ibin], hnMEL[ibin]); // LSB/nMEL
    hLSBnMELRatio[ibin]->SetTitle(Form("LSB/nMEL - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin + 1]));
    //hLSBnMELRatio[ibin]->GetYaxis()->SetRangeUser(0.8,1.2);

    /* Bg subtraction */
    hUSPminusLSB1[ibin] = subtrac(hUSPpy[ibin], hSum[ibin]);
    hUSPminusLSB1[ibin]->SetTitle(Form("USP-LSB (Sum) - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin + 1]));
    hUSPminusLSB2[ibin] = subtrac(hUSPpy[ibin], hGeoMean[ibin]);
    hUSPminusLSB2[ibin]->SetTitle(Form("USP-LSB (geo mean) - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin + 1]));
    //hUSPminusLSB2[ibin]->GetYaxis()->SetRangeUser(0.,6.4e4);
    hUSPminusMEB[ibin] = subtrac(hUSPpy[ibin], hnMEB[ibin]);
    hUSPminusMEB[ibin]->SetTitle(Form("USP-MEB - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin + 1]));
    hUSPminusMEL[ibin] = subtrac(hUSPpy[ibin], hnMEL[ibin]);
    hUSPminusMEL[ibin]->SetTitle(Form("USP-MEL - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin + 1]));

    /* Histo maquillage */
    HistoMakeUp(hGeoMean[ibin], kMagenta - 2, 20, "#it{M}_{#pi#pi}(GeV/#it{c^{2}})", "pairs");
    HistoMakeUp(hLSBRatio[ibin], kMagenta - 2, 20, "#it{M}_{#pi#pi}(GeV/#it{c^{2}})", "pairs");
    HistoMakeUp(hnMEB[ibin], kGreen - 9, 20, "#it{M}_{#pi#pi}(GeV/#it{c^{2}})", "pairs");
    HistoMakeUp(hnMEL[ibin], kGreen + 2, 20, "#it{M}_{#pi#pi}(GeV/#it{c^{2}})", "pairs");
    HistoMakeUp(hMELRatio[ibin], kGreen + 2, 20, "#it{M}_{#pi#pi}(GeV/#it{c^{2}})", "ratio");
    HistoMakeUp(hMELMEBRatio[ibin], kGreen, 20, "#it{M}_{#pi#pi}(GeV/#it{c^{2}})", "ratio");
    HistoMakeUp(hLSBnMELRatio[ibin], kRed, 20, "#it{M}_{#pi#pi}(GeV/#it{c^{2}})", "ratio");
    HistoMakeUp(hLSBnMEBRatio[ibin], kBlue, 20, "#it{M}_{#pi#pi}(GeV/#it{c^{2}})", "ratio");
    HistoMakeUp(hUSPminusLSB2[ibin], kMagenta - 2, 20, "#it{M}_{#pi#pi}(GeV/#it{c^{2}})", "pairs");
    HistoMakeUp(hUSPminusMEB[ibin], kGreen - 9, 20, "#it{M}_{#pi#pi}(GeV/#it{c^{2}})", "pairs");
    HistoMakeUp(hUSPminusMEL[ibin], kGreen + 2, 20, "#it{M}_{#pi#pi}(GeV/#it{c^{2}})", "pairs");
    HistoMakeUp(hUSPpy[ibin], kBlack, 20, "#it{M}_{#pi#pi}(GeV/#it{c^{2}})", "pairs");

    TLegend* legend1 = new TLegend(0.75, 0.75, 0.9, 0.9);
    legend1->AddEntry(hGeoMean[ibin], "LSB", "lpf");
    legend1->AddEntry(hnMEB[ibin], "nMEB", "lpf");
    legend1->AddEntry(hnMEL[ibin], "nMEL", "lpf");
    legend1->AddEntry(hUSPpy[ibin], "USP", "lpf");

    TLegend* legend2 = new TLegend(0.75, 0.75, 0.9, 0.9);
    legend2->AddEntry(hUSPminusLSB2[ibin], "USP-LSB", "lpf");
    legend2->AddEntry(hUSPminusMEB[ibin], "USP-MEB", "lpf");
    legend2->AddEntry(hUSPminusMEL[ibin], "USP-MEL", "lpf");

    TLegend* legend3 = new TLegend(0.75, 0.75, 0.9, 0.9);
    legend3->AddEntry(hLSBRatio[ibin], "LSB (--)/ LSB (++)", "lpf");
    legend3->AddEntry(hMELRatio[ibin], "MEL (--)/ MEL (++)", "lpf");

    TLegend* legend4 = new TLegend(0.75, 0.75, 0.9, 0.9);
    legend4->AddEntry(hLSBnMELRatio[ibin], "LSB/nMEL", "lpf");
    legend4->AddEntry(hLSBnMEBRatio[ibin], "LSB/nMEB", "lpf");


    /* drawing and printing output*/
    canvas[ibin] = new TCanvas(Form("c%d", ibin), Form("canvas %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin + 1]), 1600, 800);
    canvas[ibin]->Divide(2, 1);
    canvas[ibin]->cd(1);
    hUSPpy[ibin]->Draw("e");
    hGeoMean[ibin]->Draw("e same");
    hnMEB[ibin]->Draw("e same");
    hnMEL[ibin]->Draw("e same");
    legend1->Draw();
    canvas[ibin]->cd(2);
    hUSPminusLSB2[ibin]->Draw("e");
    hUSPminusMEB[ibin]->Draw("e same");
    hUSPminusMEL[ibin]->Draw("e same");
    legend2->Draw();
    canvas[ibin] = new TCanvas(Form("c%d_ratio", ibin), Form("canvas %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin + 1]), 1600, 800);
    canvas[ibin]->Divide(2, 1);
    canvas[ibin]->cd(1);
    hLSBRatio[ibin]->Draw("e");
    hMELRatio[ibin]->Draw("e same");
    legend3->Draw();
    canvas[ibin]->cd(2);
    hLSBnMELRatio[ibin]->Draw("e");
    hLSBnMEBRatio[ibin]->Draw("e same");
    legend4->Draw();

    /* writing output on file.root*/
    fout->cd();
    ptBins->Write(Form("ptBins_%d", ibin));
    hUSPpy[ibin]->Write(Form("USP_%d", ibin));
    hLikeMMpy[ibin]->Write(Form("LikeMM_%d", ibin));
    hLikePPpy[ibin]->Write(Form("LikePP_%d", ibin));
    hMEBpy[ibin]->Write(Form("MEB_%d", ibin));
    hMELMMpy[ibin]->Write(Form("MELmm_%d", ibin));
    hMELPPpy[ibin]->Write(Form("MELpp_%d", ibin));
    hGeoMean[ibin]->Write(Form("LSBGeoMean_%d", ibin));
    hSum[ibin]->Write(Form("LSBSum_%d", ibin));
    hLSBRatio[ibin]->Write(Form("LSBRatio_%d", ibin));
    hnMEB[ibin]->Write(Form("nMEB_%d", ibin));
    hMELRatio[ibin]->Write(Form("MELRatio_%d", ibin));
    hMEL[ibin]->Write(Form("MEL_%d", ibin));
    hnMEL[ibin]->Write(Form("nMEL_%d", ibin));
    hMELMEBRatio[ibin]->Write(Form("MEL/MEB_%d", ibin));
    hLSBnMEBRatio[ibin]->Write(Form("LSB/MEB_%d", ibin));
    hLSBnMELRatio[ibin]->Write(Form("nMELRatio_%d", ibin));
    hUSPminusLSB1[ibin]->Write(Form("USP-LSBSum_%d", ibin));
    hUSPminusLSB2[ibin]->Write(Form("USP-LSBGeoMean_%d", ibin));
    hUSPminusMEB[ibin]->Write(Form("USP-MEB_%d", ibin));
    hUSPminusMEL[ibin]->Write(Form("USP-MEL_%d", ibin));
  }
  //fout->Close();
  // return;
}

TH1D* geoM(TH1D* h1, TH1D* h2)
{
  //returns 2*(geometrical mean) of two histograms
  if (!h1 || !h2)
    return NULL;
  TH1D* hGeoMean = (TH1D*)h1->Clone();
  hGeoMean->Reset("ICES");
  Double_t cont1 = 0.0, cont2 = 0.0, error1 = 0.0, error2 = 0.0, error3 = 0.0, product = 0.0, geoMean = 0.0;
  for (Int_t j = 1; j < h1->GetNbinsX() + 1; j++) {
    cont1 = h1->GetBinContent(j);
    cont2 = h2->GetBinContent(j);
    error1 = h1->GetBinError(j);
    error2 = h2->GetBinError(j);
    geoMean = pow((cont1 * cont2), 1. / 2);
    product = 2 * geoMean;
    error3 = pow(cont2/cont1*(pow(error1, 2)) + cont1/cont2*(pow(error2, 2)), 1. / 2);
    //error3 = pow((pow(error1, 2) + pow(error2, 2)), 1. / 2); //NNNNOOOOOOO!!!!!!
    hGeoMean->SetBinContent(j, product);
    hGeoMean->SetBinError(j, error3);
  }
  return hGeoMean;
}

TH1D* summa(TH1D* h1, TH1D* h2)
{
  //returns sum of two histograms
  if (!h1 || !h2)
    return NULL;
  TH1D* hResult = (TH1D*)h1->Clone();
  hResult->Add(h2);
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

TH1D* normM(TH1D* h1, TH1D* h2, Int_t method /*chosen method for normalizing MEB*/)
{
  //returns normalized histograms
  //if method = 1: MEB normalized with a specific value
  //if method = 2: MEB normalized to the integral of the USP
  if (method == 1) {
    if (!h1)
      return NULL;
    Double_t normVal = 5.;
    TH1D* hResult = (TH1D*)h1->Clone();
    hResult->Scale(1. / normVal);
    return hResult;
  } else if (method == 2) {
    if (!h1 || !h2)
      return NULL;
    TH1D* hResult = (TH1D*)h1->Clone(Form("%s_MEBBB", h1->GetName()));
    //Printf("%s", hResult->GetName());
    const Double_t l[2] = { 1.15, 1.2 };
    const Int_t iMinBin = h1->GetXaxis()->FindBin(l[0]);
    const Int_t iMaxBin = h1->GetXaxis()->FindBin(l[1]);
    const Int_t iMinBin2 = h2->GetXaxis()->FindBin(l[0]);
    const Int_t iMaxBin2 = h2->GetXaxis()->FindBin(l[1]);
    //if (iMinBin != iMinBin2 || iMaxBin != iMaxBin2)
      //::Fatal("asd", "lol");
    Double_t area1 = h1->Integral(iMinBin, iMaxBin);
    Double_t area2 = h2->Integral(iMinBin2, iMaxBin2);
    Double_t yScale = area2 / area1;
    hResult->Scale(yScale);
    //Printf("%f: %f", yScale, area1 * yScale - area2);
    //getchar();
    return hResult;
  }
  return NULL;
}

TH1D* subtrac(TH1D* h1, TH1D* h2)
{
  //returns subtraction of two histograms
  if (!h1)
    return NULL;
  TH1D* hResult = (TH1D*)h1->Clone();
  hResult->Add(h2, -1);
  return hResult;
}
