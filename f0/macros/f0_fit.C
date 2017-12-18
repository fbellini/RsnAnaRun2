/////////////////////////////////////////////
//       f0 signal: fit - 21.09.2017       //
/////////////////////////////////////////////

#include "HistoMakeUp.C"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFlatte.h"
#include "RooGaussian.h"
#include "RooGlobalFunc.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooRelBW.h"
#include "SetStyle.C"
#include "TCanvas.h"
#include "TH1.h"
#include "TMath.h"
#include "TRandom.h"
#include "TTree.h"

using namespace RooFit;
using namespace RooStats;

RooPlot* fit(TH1D* h1, Double_t xMinRange, Double_t xMaxRange, Int_t fitMethod);
TPaveText* textFitResults(Int_t fitMethod, Float_t x1, Float_t y1, Float_t x2, Float_t y2);

void f0_fit(
    TString infilename = "bgSubtraction.root",
    Double_t xMinRange = 0.8,
    Double_t xMaxRange = 1.2,
    Int_t fitMethod = 2
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

  TH1D* hUSPminusLSB[10];

  for (Int_t i = 0; i < 10; i++) {
    hUSPminusLSB[i] = 0x0;
  }

  Double_t pT[11] = { 0.5, 1., 1.5, 2., 2.5, 3., 4., 5., 7., 9., 11. };
  TCanvas* canvas[10];
  RooPlot* plotframe;
  TPaveText *textFit = 0x0;
  Double_t fitParameters[15];

  for (Int_t ibin = 0; ibin < 10; ibin++) {
    hUSPminusLSB[ibin] = (TH1D*)fin->Get(Form("USP-LSBGeoMean_%d", ibin));
    if (!hUSPminusLSB[ibin]) {
      Printf("Input histogram error.");
      return;
    }


    plotframe = (RooPlot*)fit(hUSPminusLSB[ibin], xMinRange, xMaxRange, fitMethod);
    textFit = (TPaveText*)textFitResults(fitMethod, 0.05, 0.1, 0.95, 0.8);

    canvas[ibin] = new TCanvas(Form("c%d", ibin), Form("canvas %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin + 1]), 1600, 800);
    canvas[ibin]->Divide(2, 1);
    canvas[ibin]->cd(1);
    plotframe->GetXaxis()->SetTitle("#it{M}_{#pi#pi}(GeV/#it{c^{2}})");
    plotframe->Draw("e");
    canvas[ibin]->cd(2);
    textFit->Draw();
  }
  return;
}

RooPlot* fit(TH1D* h1, Double_t xMinRange, Double_t xMaxRange, Int_t fitMethod)
{
  Color_t color[] = { kRed, kGreen+2, kBlue, kMagenta, kBlack };
  Double_t histo_integral = h1->Integral();

  RooRealVar x("x", "x", xMinRange, xMaxRange);
  RooDataHist dh("dh", "dh", x, Import(*h1));
  RooPlot* frame = x.frame(Title(h1->GetTitle()));
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
  RooRealVar alpha("alpha", "alpha", -9.5, -22., 0.);
  RooExponential bkg("bkg", "Background 1", x, alpha);
  RooRealVar beta("beta", "beta", -0.8, -22., 22.);
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

  Double_t fitParameters[15];

  //fit
  switch (fitMethod) {
  case 1: {
    funcBW.fitTo(dh /*, Extended(), Save()*/);
    funcBW.plotOn(frame);
    funcBW.paramOn(frame, Layout(0.4));
    funcBW.plotOn(frame, Components(bkg), LineStyle(kDashed), LineColor(color[0]), LineWidth(3), Range(xMinRange, xMaxRange));
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
    funcRelBW.plotOn(frame, Components(bkg), LineStyle(kDashed), LineColor(color[0]), LineWidth(3), Range(xMinRange, xMaxRange));
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
    funcVoig.plotOn(frame, Components(bkg), LineStyle(kDashed), LineColor(color[0]), LineWidth(3), Range(xMinRange, xMaxRange));
    funcVoig.plotOn(frame, Components(bkg2), LineStyle(kDashed), LineColor(color[1]), LineWidth(2), Range(xMinRange, xMaxRange));
    funcVoig.plotOn(frame, Components(sigVoig), LineStyle(kDashed), LineColor(color[2]), LineWidth(2), Range(xMinRange, xMaxRange));
    RooFitResult* r3 = funcVoig.fitTo(dh, Save());
    r3->Print("v");
    r3->correlationMatrix().Print();
    Double_t sigmaVoig = sigmaF0.getVal();
    Double_t sigmaVoigErr = sigmaF0.getError();
    fitParameters[13] = sigmaVoig;
    fitParameters[14] = sigmaVoigErr;
    //RooArgSet* params = funcVoig.getParameters(x);
    break;
  }
  default: {
    Printf("Invalid method to perform the fit");
    break;
  }
  }

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

  if(fitParameters){
  fitParameters[0] = signalMass;
  fitParameters[1] = signalMassErr;
  fitParameters[2] = signalWidth;
  fitParameters[3] = signalWidthErr;
  fitParameters[4] = nSignal;
  fitParameters[5] = nSignalErr;
  fitParameters[6] = nBkg;
  fitParameters[7] = nBkgErr;
  fitParameters[8] = chi2;
  fitParameters[9] = alphaVal;
  fitParameters[10] = alphaValErr;
  fitParameters[11] = betaVal;
  fitParameters[12] = betaValErr;
  }

  return frame;
}

TPaveText* textFitResults(Int_t fitMethod, Float_t x1, Float_t y1, Float_t x2, Float_t y2)
{
  TPaveText * text = new TPaveText(x1, y1, x2, y2);
  text->SetLabel("Fit Parameters");
  text->SetBorderSize(1);
  text->SetTextAlign(21);
  text->SetFillColor(kWhite);
  text->SetTextFont(42);
  text->SetTextColor(kBlack);
  text->SetTextSize(0.07);
  text->AddText("ciao");
  /*text->AddText(Form("#it{M}_{#pi#pi} (GeV/#it{c^{2}}) = %6.4f #pm %6.4f", fitParameters[0], fitParameters[1]));
  text->AddText(Form("#Gamma (GeV) = %6.4f #pm %6.4f", fitParameters[2], fitParameters[3]));
  if (fitMethod == 3) {
    text->AddText(Form("#sigma_{Voigtian} = %6.4f #pm %6.4f", fitParameters[13], fitParameters[14]));
  }
  text->AddText(Form("N_{sig} =  %8.0f #pm %6.0f \n", fitParameters[4], fitParameters[5]));
  text->AddText(Form("N_{bkg} = %8.0f #pm %6.0f \n", fitParameters[6], fitParameters[7]));
  text->AddText(Form("#it{#alpha} = %6.4f #pm %6.4f", fitParameters[9], fitParameters[10]));
  text->AddText(Form("#it{#beta} = %6.4f #pm %6.4f", fitParameters[11], fitParameters[12]));
  text->AddText(Form("#it{#chi}^{2} = %6.4f", fitParameters[8]));*/
  return text;
}
