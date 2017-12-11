/////////////////////////////////////////////
//       f0 signal: fit - 21.09.2017       //
/////////////////////////////////////////////

#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TMath.h"
#include "TRandom.h"

using namespace RooFit ;

void fit()

{
  #ifdef __CINT__
  gROOT->ProcessLine(".x RooRelBW.cxx+");
  #endif

  TFile *file = TFile::Open("bgSubtraction.root");
  TH1D * hUSPminusLSB = (TH1D*)file->Get("USP-LSBGeoMean");

  // Setup components
  RooRealVar x("x","x",0.6,1.2);
  RooDataHist dh("dh","dh",x,Import(*hUSPminusLSB));
  RooPlot* frame = x.frame(Title("USP-LSB 2.0 < p_{T} < 2.5 GeV/#it{c}"));
  dh.plotOn(frame,DataError(RooAbsData::SumW2));

  // Signal (Voigtian PDFs)
  RooRealVar mRho("mRho","mRho",0.77526); //rho(770) invariant mass = 775.26 /pm 0.25 MeV
  RooRealVar sigmaRho("sigmaRho","sigmaRho",0.0025);
  RooRealVar widthRho( "widthRho", "widthRho", 0.149);
  RooVoigtian sigRho("sigRho","sigRho", x, mRho, widthRho, sigmaRho, kFALSE);
  //RooRelBW sigRho("sigRho","sigRho", x, mRho, widthRho);

  RooRealVar mF0("mF0","mF0",0.99); //f0(980) invariant mass = 990 /pm 20 MeV
  RooRealVar sigmaF0("sigmaF0","sigmaF0",0.02);
  RooRealVar widthF0("widthF0", "widthF0",0.05);
  RooVoigtian sigF0("sigF0", "sigF0", x, mF0, widthF0, sigmaF0, kFALSE);
  //RooRelBW sigF0("sigF0","sigF0", x, mF0, widthF0);

  RooRealVar f("f","f",0.6,1.2) ;
  RooAddPdf sig("sig", "Signal", RooArgList(sigRho, sigF0),f);

  // Background (Maxwell Boltzmann PDF)
  RooGenericPdf bkg("bkg","bkg","(2*sqrt(x/3.14)*pow(1,1.5)*exp(-x))",RooArgSet(x)) ;

  RooRealVar res("res","res", 0.6, 1.2);
  RooAddPdf func("func","func", RooArgList(sig, bkg),res);
  func.fitTo(dh, Extended());
  func.plotOn(frame);
  func.plotOn(frame,Components(bkg),LineStyle(kDashed));

  RooFitResult* r = func.fitTo(dh,Save());
  r->Print("v");

  TCanvas* c = new TCanvas("FitResult","FitResult",1200,600);
  c->Divide(2,1);
  c->cd(1);
  gPad->SetLeftMargin(0.15);
  frame->GetYaxis()->SetTitleOffset(1.4);
  frame->GetXaxis()->SetTitle("#it{M}_{#pi#pi}(GeV/#it{c^{2}})");
  frame->Draw();

/*TPaveText *pt = new TPaveText(.05,.1,.95,.8);
  pt->AddText("#frac{2s}{#pi#alpha^{2}}  ESEMPIO");
  pt->SetLabel("Born equation");
  pt->Draw();*/

return;

}
