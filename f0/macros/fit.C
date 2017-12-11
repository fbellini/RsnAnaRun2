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
  RooRealVar x("x","x",0.8,1.2);
  RooDataHist dh("dh","dh",x,Import(*hUSPminusLSB));
  RooPlot* frame = x.frame(Title(hUSPminusLSB->GetTitle()));
  dh.plotOn(frame,DataError(RooAbsData::SumW2));
  //dh->statOn(frame,Layout(0.55,0.99,0.8));

  // F0 signal
  RooRealVar mF0("mF0","mF0",0.99,0.97,1.01); //f0(980) invariant mass = 990 /pm 20 MeV
  RooRealVar sigmaF0("sigmaF0","sigmaF0",0.02);
  RooRealVar widthF0("widthF0", "widthF0",0.05,0.01,0.1);
  //RooGaussian gauss("gauss","gauss",x,mF0,sigmaF0);
  //RooRelBW relBW("relBW","relBW", x, mF0, widthF0);
  //x.setBins(10000,"cache") ;
  //RooFFTConvPdf sigF0("sigF0","Rel (X) gauss",x,relBW,gauss) ;

  RooVoigtian sigF0("sigF0", "sigF0", x, mF0, widthF0, sigmaF0, kFALSE);
  //RooBreitWigner sigF0("sigF0","sigF0", x, mF0, widthF0);

  // Bkg
  RooRealVar alpha("alpha","alpha",-8.) ;
  RooExponential bkg("bkg","Background 1",x,alpha);

  RooRealVar beta("beta", "beta", -0.5);
  RooExponential bkg2("bkg2","Background 2",x,beta);

  // Sum the background components into a composite background p.d.f.
  RooRealVar sigFrac("sigFrac","fraction of component 1 in background",0.2,0.,1.) ;
  RooAddPdf sumbkg("sumbkg","Signal",RooArgList(bkg,bkg2),sigFrac) ;


  RooRealVar bkgFrac("bkgFrac","fraction of background",0.5,0.,1.) ;
  RooAddPdf  func("model","g1+g2+a",RooArgList(sumbkg,sigF0),bkgFrac) ;


  func.fitTo(dh/*, Extended()*/);
  func.plotOn(frame);
  func.paramOn(frame,Layout(0.4));
  func.plotOn(frame,Components(bkg),LineStyle(kDashed), LineColor(kRed));
  func.plotOn(frame,Components(bkg2),LineStyle(kDashed), LineColor(kYellow));
  func.plotOn(frame,Components(sigF0),LineStyle(kDashed), LineColor(kBlue));
  //  sig.plotOn(frame,Components(sigRho),LineStyle(kDashed), LineColor(kOrange));


  RooFitResult* r = func.fitTo(dh,Save());
  r->Print("v");


  // Draw all frames on a canvas
  TCanvas* c = new TCanvas("FitResult","FitResult",1200,600);
  c->Divide(2,1);
  c->cd(1);
  gPad->SetLeftMargin(0.15);
  frame->GetYaxis()->SetTitleOffset(1.4);
  frame->GetXaxis()->SetTitle("#it{M}_{#pi#pi}(GeV/#it{c^{2}})");
  frame->Draw();

  RooArgSet* params = func.getParameters(x) ;
  RooArgSet* initParams = (RooArgSet*) params->snapshot() ;

  /*TPaveText *pt = new TPaveText(.05,.1,.95,.8);
  pt->AddText(params->printLatex(Sibling(*initParams)));
  pt->SetLabel("Born equation");
  c->cd(2);
  pt->Draw();*/

return;

}
