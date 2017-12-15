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
#include "TParticlePDG.h"
#include "HistoMakeUp.C"
#include "SetStyle.C"


using namespace RooFit;

Color_t color[]={kRed, kGreen, kBlue, kMagenta, kBlack};

Int_t f0PDG = 9010221;
TDatabasePDG *pdg = TDatabasePDG::Instance();
TParticlePDG *f0 = pdg->GetParticle(f0PDG);
//Double_t pdgMass = f0->Mass();
//Double_t pdgWidth = f0->Width();

void f0_fit(

  Int_t method = 1
  /*method allows to chose the function to perform f0 fit
  list of matches
  1 -> Breit Wigner
  2 -> relativistic Breit Wigner
  3 -> Voigtian
  4 -> Flatte' */
)

{

  #ifdef __CINT__
  gROOT->ProcessLineSync(".x RooRelBW.cxx+") ;
#endif

  TGaxis::SetMaxDigits(3);

  //input file
  TFile *file = TFile::Open("bgSubtraction.root");
  if (!file){
    Printf("Input file error.");
    return;
  }
  TH1D * hUSPminusLSB = (TH1D*)file->Get("USP-LSBGeoMean");
  if (!hUSPminusLSB){
    Printf("Input histogram error.");
    return;
  }

  // Setup components
  Double_t xMinRange=0.8;
  Double_t xMaxRange=1.2;
  RooRealVar x("x","x",xMinRange,xMaxRange);
  RooDataHist dh("dh","dh",x,Import(*hUSPminusLSB));
  RooPlot* frame = x.frame(Title(hUSPminusLSB->GetTitle()));
  dh.plotOn(frame,DataError(RooAbsData::SumW2));
  //dh->statOn(frame,Layout(0.55,0.99,0.8));

  // f0 signal
  RooRealVar mF0("mF0","mF0",0.99,0.97,1.01); //f0(980) invariant mass = 990 /pm 20 MeV
  RooRealVar sigmaF0("sigmaF0","sigmaF0",0.02);
  RooRealVar widthF0("widthF0", "widthF0",0.04,0.01,0.1);

  
  RooBreitWigner sigBW("sigF0","sigF0", x, mF0, widthF0);
  RooRelBW sigRelBW("relBW","relBW", x, mF0, widthF0);
  RooVoigtian sigVoig("sigF0", "sigF0", x, mF0, widthF0, sigmaF0, kFALSE);

  
  // residual bkg - exponential(e^alpha) + exponential(e^beta)
  RooRealVar alpha("alpha","alpha",-9.5,-11.0,-8.00);
  RooExponential bkg("bkg","Background 1",x,alpha);
  RooRealVar beta("beta", "beta", -0.8,-11.0,-7.00);
  RooExponential bkg2("bkg2","Background 2",x,beta);
  RooRealVar exp1Frac("alphaFrac","fraction of exp 1 in background",0.5,0.,1.);
  RooAddPdf sumbkg("sumbkg","Background",RooArgList(bkg,bkg2),exp1Frac);
  
  RooRealVar bkgFrac("bkgFrac","fraction of background",0.5,0.,1.) ;
  
  //model for signal + background for different p.d.f.
  RooAddPdf funcBW("model1","sig+bgBW",RooArgList(sumbkg,sigBW),bkgFrac);
  RooAddPdf funcRelBW("model2","sig+bgRelBW",RooArgList(sumbkg,sigRelBW),bkgFrac);
  RooAddPdf funcVoig("model3","sig+bgVoig",RooArgList(sumbkg,sigVoig),bkgFrac);

  //fit
  if(method==1){
    funcBW.fitTo(dh/*, Extended()*/);
    funcBW.plotOn(frame);
    funcBW.paramOn(frame,Layout(0.4));
    funcBW.plotOn(frame,Components(bkg),LineStyle(kDashed), LineColor(kRed));
    funcBW.plotOn(frame,Components(bkg2),LineStyle(kDashed), LineColor(kGreen+2));
    funcBW.plotOn(frame,Components(sigBW),LineStyle(kDashed), LineColor(kBlue));
    //chi2=xframe->chiSquare();
    RooFitResult* r1 = funcBW.fitTo(dh,Save());
    r1->Print("v");
    RooArgSet* params = funcBW.getParameters(x);
  }
  else if(method==2){
    funcRelBW.fitTo(dh/*, Extended()*/);
    funcRelBW.plotOn(frame);
    funcRelBW.paramOn(frame,Layout(0.4));
    funcRelBW.plotOn(frame,Components(bkg),LineStyle(kDashed), LineColor(kRed));
    funcRelBW.plotOn(frame,Components(bkg2),LineStyle(kDashed), LineColor(kGreen+2));
    funcRelBW.plotOn(frame,Components(sigRelBW),LineStyle(kDashed), LineColor(kBlue));
    //chi2=xframe->chiSquare();
    RooFitResult* r2 = funcRelBW.fitTo(dh,Save());
    r2->Print("v");
    RooArgSet* params = funcRelBW.getParameters(x);
  }
  else if(method==3){
    funcVoig.fitTo(dh/*, Extended()*/);
    funcVoig.plotOn(frame);
    funcVoig.paramOn(frame,Layout(0.4));
    funcVoig.plotOn(frame,Components(bkg),LineStyle(kDashed), LineColor(kRed));
    funcVoig.plotOn(frame,Components(bkg2),LineStyle(kDashed), LineColor(kGreen+2));
    funcVoig.plotOn(frame,Components(sigVoig),LineStyle(kDashed), LineColor(kBlue));
    //chi2=xframe->chiSquare();
    RooFitResult* r3 = funcVoig.fitTo(dh,Save());
    r3->Print("v");
    RooArgSet* params = funcVoig.getParameters(x);
  }

  // Draw all frames on a canvas
  TCanvas* c = new TCanvas("FitResult","FitResult",1200,600);
  c->Divide(2,1);
  c->cd(1);
  gPad->SetLeftMargin(0.2);
  frame->GetYaxis()->SetTitleOffset(1.4);
  frame->GetXaxis()->SetTitle("#it{M}_{#pi#pi}(GeV/#it{c^{2}})");
  frame->Draw();


  //RooArgSet* initParams = (RooArgSet*) params->snapshot() ;

return;

}
