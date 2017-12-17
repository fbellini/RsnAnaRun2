/////////////////////////////////////////////
//       f0 signal: fit - 21.09.2017       //
/////////////////////////////////////////////

#include "HistoMakeUp.C"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooGlobalFunc.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooRelBW.h"
#include "RooFlatte.h"
#include "SetStyle.C"
#include "TCanvas.h"
#include "TH1.h"
#include "TMath.h"
#include "TRandom.h"
#include "TTree.h"

using namespace RooFit;
using namespace RooStats;

Color_t color[] = { kRed, kGreen, kBlue, kMagenta, kBlack };

void f0_fit(
    TString infilename = "bgSubtraction.root",
    Double_t xMinRange = 0.8,
    Double_t xMaxRange = 1.2,
    Double_t iMinBinPt = 1.0,
    Double_t iMaxBinPt = 1.5,
    Int_t fitMethod = 3
    /*fitMethod allows to chose the function to perform f0 fit
  list of matches
  1 -> Breit Wigner
  2 -> relativistic Breit Wigner
  3 -> Voigtian
  4 -> Flatte' */
)

{
  TGaxis::SetMaxDigits(3);

#ifdef __CINT__
  gROOT->ProcessLineSync(".x RooRelBW.cxx+");
  gROOT->ProcessLineSync(".x RooFlatte.cxx+");
#endif

  //input file
  if (!infilename) {
    Printf("ERROR: invalid file name.");
    return;
  }
  TFile* fin = TFile::Open(infilename.Data());
  if (!fin) {
    Printf("ERROR: invalid input file.");
    return;
  }

  TH1D* hUSPminusLSB = (TH1D*)fin->Get("USP-LSBGeoMean");
  if (!hUSPminusLSB) {
    Printf("Input histogram error.");
    return;
  }

  Double_t histo_integral = hUSPminusLSB->Integral();

  TAxis* ptBins = (TAxis*)fin->Get("ptBins");
  Int_t nPt = ptBins->GetNbins();
  //Double_t binWidth = ptBins->GetBinWidth();
  Int_t dimPt = nPt + 1;
  Double_t pT[dimPt];
  for (Int_t k = 0; k < dimPt; k++) {
    pT[k] = ptBins->GetBinLowEdge(k + 1);
  }

  Double_t fitParams[15];

  // Setup components
  RooRealVar x("x", "x", xMinRange, xMaxRange);
  RooDataHist dh("dh", "dh", x, Import(*hUSPminusLSB));
  RooPlot* frame = x.frame(Title(hUSPminusLSB->GetTitle()));
  dh.plotOn(frame, DataError(RooAbsData::SumW2));
  //dh.statOn(frame,Layout(0.55,0.99,0.8));

  // f0 signal
  RooRealVar mF0("mF0", "mF0", 0.99, 0.97, 1.01); //f0(980) invariant mass = 990 /pm 20 MeV
  RooRealVar sigmaF0("sigmaF0", "sigmaF0", 0.02, 0., 0.1);
  RooRealVar widthF0("widthF0", "widthF0", 0.04, 0.01, 0.1);
  RooBreitWigner sigBW("sigBW", "sigBW", x, mF0, widthF0);
  RooRelBW sigRelBW("relBW", "relBW", x, mF0, widthF0);
  RooVoigtian sigVoig("sigVoig", "sigVoig", x, mF0, widthF0, sigmaF0, kFALSE);

  // residual bkg - exponential(e^alpha) + exponential(e^beta)
  RooRealVar alpha("alpha", "alpha", -9.5, -17.5, -8.5);
  RooExponential bkg("bkg", "Background 1", x, alpha);
  RooRealVar beta("beta", "beta", -0.8, -4., 25.);
  RooExponential bkg2("bkg2", "Background 2", x, beta);
  RooRealVar exp1Frac("alphaFrac", "fraction of exp 1 in background", 0.5, 0., 1.);
  RooAddPdf sumbkg("sumbkg", "Background", RooArgList(bkg, bkg2), exp1Frac);
  RooRealVar bkgFrac("bkgFrac", "fraction of background", 0.5, 0., 1.);

  //model for signal + background for different p.d.f.
  RooRealVar nsig("nsig", "signal fraction", histo_integral * 0.5, 0., histo_integral);
  RooRealVar nbkg("nbkg", "background fraction", histo_integral * 0.5, 0., histo_integral);
  RooAddPdf funcBW("model1", "sig+bgBW", RooArgList(sumbkg, sigBW), RooArgList(nsig, nbkg));
  RooAddPdf funcRelBW("model2", "sig+bgRelBW", RooArgList(sumbkg, sigRelBW), RooArgList(nsig, nbkg));
  RooAddPdf funcVoig("model3", "sig+bgVoig", RooArgList(sumbkg, sigVoig), RooArgList(nsig, nbkg));

  //fit
  switch (fitMethod) {
  case 1: {
    funcBW.fitTo(dh /*, Extended(), Save()*/);
    funcBW.plotOn(frame);
    funcBW.paramOn(frame, Layout(0.4));
    funcBW.plotOn(frame, Components(bkg), LineStyle(kDashed), LineColor(color[0]), LineWidth(2), Range(xMinRange, xMaxRange));
    funcBW.plotOn(frame, Components(bkg2), LineStyle(kDashed), LineColor(color[1]), LineWidth(2), Range(xMinRange, xMaxRange));
    funcBW.plotOn(frame, Components(sigBW), LineStyle(kDashed), LineColor(color[2]), LineWidth(2), Range(xMinRange, xMaxRange));
    RooFitResult* r1 = funcBW.fitTo(dh, Save());
    r1->Print("v");
    r1->correlationMatrix().Print();
    //RooArgSet* params = funcBW.getParameters(x);
    break;
  }
  case 2: {
    funcRelBW.fitTo(dh /*, Extended(), Save()*/);
    funcRelBW.plotOn(frame);
    funcRelBW.paramOn(frame, Layout(0.4));
    funcRelBW.plotOn(frame, Components(bkg), LineStyle(kDashed), LineColor(color[0]), LineWidth(2), Range(xMinRange, xMaxRange));
    funcRelBW.plotOn(frame, Components(bkg2), LineStyle(kDashed), LineColor(color[1]), LineWidth(2), Range(xMinRange, xMaxRange));
    funcRelBW.plotOn(frame, Components(sigRelBW), LineStyle(kDashed), LineColor(color[2]), LineWidth(2), Range(xMinRange, xMaxRange));
    RooFitResult* r2 = funcRelBW.fitTo(dh, Save());
    r2->Print("v");
    r2->correlationMatrix().Print();
    //RooArgSet* params = funcRelBW.getParameters(x);
    break;
  }
  case 3: {
    funcVoig.fitTo(dh /*, Extended(), Save()*/);
    funcVoig.plotOn(frame);
    funcVoig.paramOn(frame, Layout(0.4));
    funcVoig.plotOn(frame, Components(bkg), LineStyle(kDashed), LineColor(color[0]), LineWidth(2), Range(xMinRange, xMaxRange));
    funcVoig.plotOn(frame, Components(bkg2), LineStyle(kDashed), LineColor(color[1]), LineWidth(2), Range(xMinRange, xMaxRange));
    funcVoig.plotOn(frame, Components(sigVoig), LineStyle(kDashed), LineColor(color[2]), LineWidth(2), Range(xMinRange, xMaxRange));
    RooFitResult* r3 = funcVoig.fitTo(dh, Save());
    r3->Print("v");
    r3->correlationMatrix().Print();
    Double_t sigmaVoig = sigmaF0.getVal();
    Double_t sigmaVoigErr = sigmaF0.getError();
    fitParams[13] = sigmaVoig;
    fitParams[14] = sigmaVoigErr;
    //RooArgSet* params = funcVoig.getParameters(x);
    break;
  }
  default: {
    Printf("Invalid method to perform the fit");
    break;
  }
  }

  //output
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);

  TCanvas* c = new TCanvas("c", "f0_fit", 1200, 600);
  c->Divide(2, 1);

  c->cd(1);
  gPad->SetLeftMargin(0.2);
  frame->GetYaxis()->SetTitleOffset(1);
  frame->GetXaxis()->SetTitle("#it{M}_{#pi#pi}(GeV/#it{c^{2}})");
  frame->Draw();

  c->cd(2);
  Double_t signalMass = mF0.getVal();
  Double_t signalMassErr = mF0.getError();
  Double_t signalWidth = widthF0.getVal();
  Double_t signalWidthErr = widthF0.getError();
  Double_t nSignal = nsig.getVal();
  Double_t nSignalErr = nsig.getError();
  Double_t nBkg = nbkg.getVal();
  Double_t nBkgErr = nbkg.getError();
  Double_t chi2 = frame->chiSquare();
  Double_t alphaVal = alpha.getVal();
  Double_t alphaValErr = alpha.getError();
  Double_t betaVal = beta.getVal();
  Double_t betaValErr = beta.getError();

  fitParams[0] = signalMass;
  fitParams[1] = signalMassErr;
  fitParams[2] = signalWidth;
  fitParams[3] = signalWidthErr;
  fitParams[4] = nSignal;
  fitParams[5] = nSignalErr;
  fitParams[6] = nBkg;
  fitParams[7] = nBkgErr;
  fitParams[8] = chi2;
  fitParams[9] = alphaVal;
  fitParams[10] = alphaValErr;
  fitParams[11] = betaVal;
  fitParams[12] = betaValErr;

  TPaveText* pt = new TPaveText(.05, .1, .95, .8);
  pt->SetLabel("Fit Parameters");
  pt->SetBorderSize(1);
  pt->SetFillColor(kWhite);
  pt->SetTextFont(42);
  pt->SetTextColor(kBlack);
  pt->AddText(Form("#it{M}_{#pi#pi} (GeV/#it{c^{2}}) = %6.4f #pm %6.4f", fitParams[0], fitParams[1]));
  pt->AddText(Form("#Gamma (GeV) = %6.4f #pm %6.4f", fitParams[2], fitParams[3]));
  if (fitMethod == 3) {
    pt->AddText(Form("#sigma_{Voigtian} = %6.4f #pm %6.4f", fitParams[13], fitParams[14]));
  }
  pt->AddText(Form("N_{sig} =  %8.0f #pm %6.0f \n", fitParams[4], fitParams[5]));
  pt->AddText(Form("N_{bkg} = %8.0f #pm %6.0f \n", fitParams[6], fitParams[7]));
  pt->AddText(Form("#it{#alpha} = %6.4f #pm %6.4f", fitParams[9], fitParams[10]));
  pt->AddText(Form("#it{#beta} = %6.4f #pm %6.4f", fitParams[11], fitParams[12]));
  pt->AddText(Form("#it{#chi}^{2} = %6.4f", fitParams[8]));
  pt->Draw();

  /*TH2F* hMvsPt = new TH2F("hMVsPt", "f_{0} mass; p_{T} (GeV/#it{c}); M (GeV/#it{c^{2}})", 400, iMinBinPt, iMaxBinPt, 400, xMinRange, xMaxRange);
  TH1F* hWidthVsPt = new TH1F("hWidthVsPt", "f_{0} width; p_{T} (GeV/#it{c}); #Gamma (GeV)", 400, iMinBinPt, iMaxBinPt);
  TH1F* hSigVsPt = new TH1F("hSigVsPt", "f_{0} raw yield; p_{T} (GeV/#it{c}); raw yield, dN/dp_{T}", 400, iMinBinPt, iMaxBinPt);
  TH1F* hSigmaVoigVsPt = new TH1F("hSigmaVoigVsPt", "f_{0} sigma Voigtian; p_{T} (GeV/#it{c}); ", 400, iMinBinPt, iMaxBinPt);

  hMvsPt->SetBinContent(iMaxBinPt,fitParams[0]);
  hMvsPt->SetBinError(iMaxBinPt,fitParams[1]);

  TCanvas* c2 = new TCanvas("c2", "c2", 1200, 1200);
  c2->Divide(2, 2);
  c2->cd(1); hMvsPt->Draw("hist");
  c2->cd(2); hWidthVsPt->Draw();
  c2->cd(3); hSigVsPt->Draw();
  c2->cd(4); hSigmaVoigVsPt->Draw();*/

  return;
}
