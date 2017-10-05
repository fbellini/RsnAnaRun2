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
  TH1D * hUSPminusLSB = (TH1D*)file.Get("USP-LSBGeoMean");

  // Setup components
  RooRealVar x("x","x",0.6,1.2);
  RooDataHist dh("dh","dh",x,Import(*hUSPminusLSB));
  RooPlot* frame = x.frame(Title(hUSPminusLSB->GetTitle()));
  dh.plotOn(frame,DataError(RooAbsData::SumW2));

  // Signal 
  RooRealVar mRho("mRho","mRho",0.77526, 0.77501, 0.77551); //rho(770) invariant mass = 775.26 /pm 0.25 MeV
  RooRealVar widthRho( "widthRho", "widthRho", 0.149);
  RooRelBW sigRho("sigRho","sigRho", x, mRho, widthRho);

  RooRealVar mF0("mF0","mF0",0.97, 1.01); //f0(980) invariant mass = 990 /pm 20 MeV
  RooRealVar widthF0("widthF0", "widthF0",0.05);
  RooRelBW sigF0("sigF0","sigF0", x, mF0, widthF0);

  RooRealVar sig1frac("sig1frac","fraction of component 1 in signal",0.1,0.0001,0.2) ;
  RooAddPdf sig("sig","Signal",RooArgList(sigRho,sigF0),sig1frac) ;


  // Background (pol2 + exponential)
  RooRealVar p0("p0", "p0", 17500., 15000., 20000.);
  RooRealVar p1("p1", "p1", 500., 1., 1000.);
  RooRealVar p2("p2", "p2", 500., 1., 1000.);
  RooPolynomial bkg("bkb","bkg", x, RooArgList(p0, p1, p2));


  //RooRealVar alpha("alpha","alpha",-1) ;
  //RooExponential bkg2("bkg2","Background 2",x,alpha) ;

  RooRealVar MBNorm("MBNorm", "MBNorm", 2.,0., 5.);
  RooGenericPdf bkg2("bkg2","bkg2","(MBNorm*sqrt(x/TMath::pi)*pow(1,1.5)*exp(-x))",RooArgSet(x, MBNorm)) ;

  RooRealVar sum2("sum2", "sum2", 0.6, 1.2);
  RooAddPdf sumbkg("sumbkg","sumbkg",RooArgList(bkg,bkg2),sum2) ;

  RooRealVar res("res","res", 0.6, 1.2);
  RooAddPdf func("func","func", RooArgList(sig, sumbkg),res);
  func.fitTo(dh, Extended());
  func.plotOn(frame);
  func.plotOn(frame,Components(sumbkg),LineStyle(kDashed), LineColor(kRed));
  func.plotOn(frame,Components(bkg2),LineStyle(kDashed), LineColor(kYellow));
  func.plotOn(frame,Components(sigF0),LineStyle(kDashed), LineColor(kBlue));
  func.plotOn(frame,Components(sigRho),LineStyle(kDashed), LineColor(kOrange)) RooRealVar res("res","res", 0.6, 1.2);

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
