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

RooPlot* fit(TH1D* h1, Double_t xMinRange, Double_t xMaxRange, Int_t fitMethod, Double_t* fitParameters, Bool_t useChi2, Int_t* fitStatus, Double_t alphaStart, Double_t alphaMin, Double_t alphaMax, Double_t betaStart, Double_t betaMin, Double_t betaMax);
Bool_t isFitFailed(Int_t* fitStatus);
RooFitResult* fitChi2(RooChi2Var* chi2Var, Int_t* fitStatus, Bool_t minos = 1, Bool_t improve = 0);
//RooFitResult* fitLikelihood(RooNLLVar* nll, Int_t* fitStatus, Bool_t minos = 1, Bool_t improve = 0);
TPaveText* textFitResults(Int_t fitMethod, Float_t x1, Float_t y1, Float_t x2, Float_t y2, Double_t* fitParameters);

void f0_fit(
    //TString inputFile = "bgSubtraction_RsnOut_f0_2sTOF.root",
    TString inputFile = "bgSubtraction_RsnOut_f0_2sTPC_3sTOFveto.root",
    //TString inputFile = "bgSubtraction_RsnOut_f0_3sTPC_3sTOFveto.root",
    //TString inputFile = "bgSubtraction_RsnOut_f0_2sTPC_4sTOFveto.root",
    Bool_t useChi2 = kTRUE,
    Double_t xMinRange = 0.8,
    Double_t xMaxRange = 1.3,
    Int_t fitMethod = 4
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

  gROOT->ForceStyle();
  gStyle->SetOptStat(0);
  gStyle->SetLineWidth(1);
  gStyle->SetMarkerStyle(3); //full triangle down = 22; star = 3; wide circle = 24; wide triangle up = 26
  //gStyle->SetMarkerColor(kRed+1);
  gStyle->SetHistLineWidth(1.5);
  gStyle->SetHistLineColor(kBlue + 1); //BW = kRed+1; relBW = kBlue+1; Voigt = kGreen+2; Flatte' = kMagenta;

  //input file
  if (!inputFile) {
    Printf("Invalid file name.");
    return;
  }
  TFile* fin = TFile::Open(inputFile.Data());
  if (!fin) {
    Printf("Invalid input file.");
    return;
  }

  TCanvas* canvas[9];
  RooPlot* plot = 0x0;
  TPaveText* textFit = 0x0;

  TAxis* bins = (TAxis*)fin->Get("bins");
  Int_t numPt = bins->GetNbins();

  const Int_t dimPt = numPt + 1;
  Double_t pT[dimPt];
  for (Int_t z = 0; z < dimPt; z++) {
    pT[z] = bins->GetBinLowEdge(z + 1);
  }

  TFile* fout = TFile::Open("f0_fit_2sTPC_3sTOFveto_Flatte_mFree_GammaFree.root", "RECREATE");
  //TFile* fout = TFile::Open("f0_fit_3sTPC_3sTOFveto.root", "RECREATE");
  //TFile* fout = TFile::Open("f0_fit_2sTPC_4sTOFveto.root", "RECREATE");
  //TFile* fout = TFile::Open("f0_fit_2sTOF.root", "RECREATE");

  //input histograms
  TH1D* hUSPminusLSB[9];
  for (Int_t i = 0; i < 9; i++) {
    hUSPminusLSB[i] = 0x0;
  }

  //define histograms for fit parameters vs pt
  TH1F* hMvsPt = new TH1F("hMvsPt", "f_{0} mass; p_{T} (GeV/#it{c}); M (GeV/#it{c^{2}})", numPt, pT);
  TH1F* hWidthVsPt = new TH1F("hWidthVsPt", "f_{0} width; p_{T} (GeV/#it{c}); #Gamma (GeV)", numPt, pT);
  TH1F* hSigVsPt = new TH1F("hSigVsPt", "f_{0} raw yield; p_{T} (GeV/#it{c}); raw yield, dN/dp_{T}", numPt, pT);
  TH1F* hSigmaVoigVsPt = new TH1F("hSigmaVoigVsPt", "f_{0} #sigma_{Voigtian}; p_{T} (GeV/#it{c}); #sigma_{Voigtian}", numPt, pT);
  TH1F* hSigOverBkgVsPt = new TH1F("hSigOverBkgVsPt", "f_{0} sig/bkg; p_{T} (GeV/#it{c}); S/B", numPt, pT);
  TH1F* hSignificanceVsPt = new TH1F("hSignificanceVsPt", "f_{0} sig/#sqrt(sig+bkg); p_{T} (GeV/#it{c}); S/#sqrt(S+B)", numPt, pT);
  //TH1F* hSigEvtVsPt = new TH1F("hSigEvtVsPt", "f_{0} per-event yield; p_{T} (GeV/#it{c}); per-event yield, dN/dp_{T}", numPt, pT);

  //Loop on pt bins
  for (Int_t ibin = 0; ibin < 9; ibin++) {

    hUSPminusLSB[ibin] = (TH1D*)fin->Get(Form("USP-LSBGeoMean_%d", ibin));
    if (!hUSPminusLSB[ibin]) {
      Printf("Input histogram error.");
      return;
    }

    Int_t iMinBinPt = pT[ibin];
    Int_t iMaxBinPt = pT[ibin + 1];

    //define fit parameters and status flags
    Double_t bWidth;
    Double_t fitParameters[17];
    for (Int_t ipar = 0; ipar < 17; ipar++) {
      fitParameters[ipar] = -1.0;
    }
    Double_t sigOverBkg = -1.0;
    Double_t significance = -1.0;
    Int_t fitStatus[4] = { -1, -1, -1, -1 };

    //fit parameters
    //                         0     1     2     3      4      5     6      7     8
    Float_t alphaStart[9] = { -8., -12., -12.,  -12.,  -11., -11.,  -7.,  -10.,  -7. };
    Float_t alphaMin[9] = {  -22., -22., -16.,  -13.,  -13., -13., -13.,  -22.,  -9. };
    Float_t alphaMax[9] = {   -7.,  -7.,  -4.,   -4.,   -4.,  -4.,  -4.,   -4.,  -3. };
    Float_t betaStart[9] = {  -4.,  -4.,  -2.,    1.,    1.,   1.,   2.,    3.,   3. };
    Float_t betaMin[9] = {    -5.,  -5.,  -4.,   -2.,   -2.,  -2.,  -2.,    0.,   0. };
    Float_t betaMax[9] = {     0.,   0.,   2.,    2.,    5.,   5.,   5.,    7.,   7. };

    //run fit
    plot = (RooPlot*)fit(hUSPminusLSB[ibin], xMinRange, xMaxRange, fitMethod, fitParameters, kTRUE, fitStatus, alphaStart[ibin], alphaMin[ibin], alphaMax[ibin], betaStart[ibin], betaMin[ibin], betaMax[ibin]);

    //check fit status and save fit output if fit succeeded
    if (!isFitFailed(fitStatus)) {

      textFit = (TPaveText*)textFitResults(fitMethod, 0.05, 0.1, 0.95, 0.8, fitParameters);

      canvas[ibin] = new TCanvas(Form("c%d", ibin), Form("canvas %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin + 1]), 1600, 800);
      canvas[ibin]->Divide(2, 1);
      canvas[ibin]->cd(1);
      plot->GetXaxis()->SetTitle("#it{M}_{#pi#pi}(GeV/#it{c^{2}})");
      plot->Draw("e");
      canvas[ibin]->cd(2);
      textFit->Draw();
      canvas[ibin]->Print(Form("fit_method%d_%2.1f<pT<%2.1f_mFree_GammaFree.png", fitMethod, pT[ibin], pT[ibin + 1]));

      //protection: check that fitParameters exists
      if (fitParameters) {

        //get significance and S/B
        sigOverBkg = fitParameters[4] / fitParameters[6];
        significance = fitParameters[4] / TMath::Sqrt(fitParameters[4] + fitParameters[6]);

        //fill histograms
        hMvsPt->SetBinContent(ibin + 1, fitParameters[0]);
        hMvsPt->SetBinError(ibin + 1, fitParameters[1]);
        hWidthVsPt->SetBinContent(ibin + 1, fitParameters[2]);
        hWidthVsPt->SetBinError(ibin + 1, fitParameters[3]);
        bWidth = bins->GetBinWidth(ibin + 1);
        hSigVsPt->SetBinContent(ibin + 1, fitParameters[4] / bWidth);
        hSigVsPt->SetBinError(ibin + 1, fitParameters[5] / bWidth);
        if (fitMethod == 3) {
          hSigmaVoigVsPt->SetBinContent(ibin + 1, fitParameters[13]);
          hSigmaVoigVsPt->SetBinError(ibin + 1, fitParameters[14]);
        }
        hSigOverBkgVsPt->SetBinContent(ibin + 1, sigOverBkg);
        hSigOverBkgVsPt->SetBinError(ibin + 1, sigOverBkg);
        hSignificanceVsPt->SetBinContent(ibin + 1, significance);
        hSignificanceVsPt->SetBinError(ibin + 1, significance);
        //hSigEvtVsPt->SetBinContent(ibin + 1, fitParameters[4] / bWidth);
        //hSigEvtVsPt->SetBinError(ibin + 1, fitParameters[5] / bWidth);

      } else {
        Printf("Ivalid fitParameters array. Check!!!");
      }
    } else {
      Printf("Fit failed!");
    }

    Int_t n;
    printf("Please enter 1 to continue \n");
    scanf("%d", &n);
    printf("That's all folks! \n");

  } // end loop on pt bins

  Double_t yMin = 0.9;
  Double_t yMax = 1.08;
  Double_t yMin2 = 0.0;
  Double_t yMax2 = 0.11;

  //plot fitted parameters vs pt
  TCanvas* c_histos = new TCanvas("c_histos", "c_histos", 1200, 600);
  c_histos->Divide(3, 2);
  c_histos->cd(1);
  hMvsPt->GetYaxis()->SetRangeUser(yMin, yMax);
  hMvsPt->Draw();
  TLine* lineMass1 = new TLine(0.5, 0.96, 9., 0.96);
  lineMass1->SetLineColor(kGray + 1);
  lineMass1->SetLineStyle(3);
  lineMass1->Draw();
  TLine* lineMass2 = new TLine(0.5, 1.01, 9., 1.01);
  lineMass2->SetLineColor(kGray + 1);
  lineMass2->SetLineStyle(3);
  lineMass2->Draw();
  c_histos->cd(2);
  hWidthVsPt->GetYaxis()->SetRangeUser(yMin2, yMax2);
  hWidthVsPt->Draw();
  TLine* lineWidth1 = new TLine(0.5, 0.01, 9., 0.01);
  lineWidth1->SetLineColor(kGray + 1);
  lineWidth1->SetLineStyle(3);
  lineWidth1->Draw();
  TLine* lineWidth2 = new TLine(0.5, 0.1, 9., 0.1);
  lineWidth2->SetLineColor(kGray + 1);
  lineWidth2->SetLineStyle(3);
  lineWidth2->Draw();
  c_histos->cd(3);
  gPad->SetLogy();
  hSigVsPt->Draw();
  c_histos->cd(4);
  hSigmaVoigVsPt->Draw();
  c_histos->cd(5);
  hSigOverBkgVsPt->Draw("hist");
  c_histos->cd(6);
  hSignificanceVsPt->Draw("hist");
  c_histos->Print("fitParameters.png");

  //save output to file
  fout->cd();
  hMvsPt->Write();
  hWidthVsPt->Write();
  hSigVsPt->Write();
  hSigmaVoigVsPt->Write();
  hSigOverBkgVsPt->Write();
  hSignificanceVsPt->Write();
  bins->Write();
  //fout->Close();

  return;
}

RooPlot* fit(TH1D* h1, Double_t xMinRange, Double_t xMaxRange, Int_t fitMethod, Double_t* fitParameters, Bool_t useChi2, Int_t* fitStatus, Double_t alphaStart, Double_t alphaMin, Double_t alphaMax, Double_t betaStart, Double_t betaMin, Double_t betaMax)
{
  Color_t color[] = { kRed, kGreen + 2, kBlue, kMagenta, kBlack };
  Double_t histo_integral = h1->Integral();
  RooRealVar x("x", "x", xMinRange, xMaxRange);
  RooDataHist dh("dh", "dh", x, Import(*h1));
  RooPlot* frame = x.frame(Title(h1->GetTitle()));
  dh.plotOn(frame, DataError(RooAbsData::SumW2));
  //dh.statOn(frame,Layout(0.55,0.99,0.8));

  // f0 signal
  RooRealVar mF0("mF0", "mF0", 0.99, 0.96, 1.01); //f0(980) invariant mass = 990 /pm 20 MeV
  RooRealVar widthF0("widthF0", "widthF0", 0.05/*, 0.01, 0.1*/);
  RooRealVar sigmaF0("sigmaF0", "sigmaF0", 0.003, 0., 0.01);
  RooRealVar g0("g0", "g0", 199);
  RooRealVar g1("g1", "g1", 199*3);
  RooRealVar m0a("m0a", "m0a", 0.1396);
  RooRealVar m0b("m0b", "m0b", 0.1396);
  RooRealVar m1a("m1a", "m1a", 0.4937);
  RooRealVar m1b("m1b", "m1b", 0.4937);

  RooBreitWigner sigBW("sigBW", "sigBW", x, mF0, widthF0);
  RooRelBW sigRelBW("relBW", "relBW", x, mF0, widthF0);
  RooVoigtian sigVoig("sigVoig", "sigVoig", x, mF0, widthF0, sigmaF0, kFALSE);
  RooFlatte sigFlatte("sigFlatte", "sigFlatte", x, mF0, g0, m0a, m0b, g1, m1a, m1b);


  // residual bkg - exponential(e^alpha) + exponential(e^beta)
  RooRealVar alpha("alpha", "alpha", alphaStart, alphaMin, alphaMax);
  RooExponential bkg("bkg", "Background 1", x, alpha);
  RooRealVar beta("beta", "beta", betaStart, betaMin, betaMax);
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
  RooAddPdf funcFlatte("model4", "sig+bgFlatte", RooArgList(sumbkg, sigFlatte), RooArgList(nsig, nbkg));

  RooFitResult* r1;
  RooFitResult* r2;
  RooFitResult* r3;
  RooFitResult* r4;

  //take fitStatus from outside -- do not redefine with: Int_t fitStatus[4] = { 1, 1, 1, 1 };
  if (!fitStatus) {
    Printf("ERROR: fitStatus array not defined");
    return 0x0;
  }
  Bool_t improve = kFALSE;
  Bool_t minos = kTRUE;

  //fit
  switch (fitMethod) {
  case 1: {
    RooChi2Var* chi2Var = new RooChi2Var("chi2Var", "chi2Var", funcBW, dh, kTRUE);
    RooNLLVar* nll = new RooNLLVar("nll", "nll", funcBW, dh);
    if (useChi2)
      r1 = (RooFitResult*)fitChi2(chi2Var, fitStatus, minos, improve);
    //  else
    //  r1 = (RooFitResult*)fitLikelihood(nll, fitStatus, minos, improve);
    //funcBW.fitTo(dh /*, Extended(), Save()*/);
    //RooNLLVar *nll = new RooNLLVar("nll","-log(L)",funcBW,dh,Extended(kTRUE),Verbose(kTRUE), NumCPU(1));
    //RooChi2Var *chi2 = new RooChi2Var("chi2","chi2",funcBW,dh,Extended(kTRUE),Verbose(kTRUE), NumCPU(1));
    //RooMinuit* m = new RooMinuit(*chi2);
    funcBW.plotOn(frame);
    funcBW.paramOn(frame, Layout(0.4));
    funcBW.plotOn(frame, Components(bkg), LineStyle(kDashed), LineColor(color[0]), LineWidth(4), Range(xMinRange, xMaxRange));
    funcBW.plotOn(frame, Components(bkg2), LineStyle(kDashed), LineColor(color[1]), LineWidth(4), Range(xMinRange, xMaxRange));
    funcBW.plotOn(frame, Components(sigBW), LineStyle(kDashed), LineColor(color[2]), LineWidth(4), Range(xMinRange, xMaxRange));
    /*RooFitResult* */ r1 = funcBW.fitTo(dh, Save());
    r1->Print("v");
    r1->correlationMatrix().Print();
    //RooArgSet* params = funcBW.getParameters(x);
    break;
  }
  case 2: {
    RooChi2Var* chi2Var = new RooChi2Var("chi2Var", "chi2Var", funcRelBW, dh, kTRUE);
    RooNLLVar* nll = new RooNLLVar("nll", "nll", funcRelBW, dh);
    if (useChi2)
      r2 = (RooFitResult*)fitChi2(chi2Var, fitStatus, minos, improve);

    funcRelBW.plotOn(frame);
    funcRelBW.paramOn(frame, Layout(0.4));
    funcRelBW.plotOn(frame, Components(bkg), LineStyle(kDashed), LineColor(color[0]), LineWidth(4), Range(xMinRange, xMaxRange));
    funcRelBW.plotOn(frame, Components(bkg2), LineStyle(kDashed), LineColor(color[1]), LineWidth(4), Range(xMinRange, xMaxRange));
    funcRelBW.plotOn(frame, Components(sigRelBW), LineStyle(kDashed), LineColor(color[2]), LineWidth(4), Range(xMinRange, xMaxRange));
    /*RooFitResult* */ r2 = funcRelBW.fitTo(dh, Save());
    r2->Print("v");
    r2->correlationMatrix().Print();
    //RooArgSet* params = funcRelBW.getParameters(x);
    break;
  }
  case 3: {
    RooChi2Var* chi2Var = new RooChi2Var("chi2Var", "chi2Var", funcVoig, dh, kTRUE);
    RooNLLVar* nll = new RooNLLVar("nll", "nll", funcVoig, dh);
    if (useChi2)
      r3 = (RooFitResult*)fitChi2(chi2Var, fitStatus, minos, improve);

    funcVoig.plotOn(frame);
    funcVoig.paramOn(frame, Layout(0.4));
    funcVoig.plotOn(frame, Components(bkg), LineStyle(kDashed), LineColor(color[0]), LineWidth(4), Range(xMinRange, xMaxRange));
    funcVoig.plotOn(frame, Components(bkg2), LineStyle(kDashed), LineColor(color[1]), LineWidth(4), Range(xMinRange, xMaxRange));
    funcVoig.plotOn(frame, Components(sigVoig), LineStyle(kDashed), LineColor(color[2]), LineWidth(4), Range(xMinRange, xMaxRange));
    /*RooFitResult**/ r3 = funcVoig.fitTo(dh, Save());
    r3->Print("v");
    r3->correlationMatrix().Print();
    Double_t sigmaVoig = sigmaF0.getVal();
    Double_t sigmaVoigErr = sigmaF0.getError();
    fitParameters[13] = sigmaVoig;
    fitParameters[14] = sigmaVoigErr;
    //RooArgSet* params = funcVoig.getParameters(x);
    break;
  }
  case 4: {
    RooChi2Var* chi2Var = new RooChi2Var("chi2Var", "chi2Var", funcFlatte, dh, kTRUE);
    RooNLLVar* nll = new RooNLLVar("nll", "nll", funcFlatte, dh);
    if (useChi2)
      r4 = (RooFitResult*)fitChi2(chi2Var, fitStatus, minos, improve);

    funcFlatte.plotOn(frame);
    funcFlatte.paramOn(frame, Layout(0.4));
    funcFlatte.plotOn(frame, Components(bkg), LineStyle(kDashed), LineColor(color[0]), LineWidth(4), Range(xMinRange, xMaxRange));
    funcFlatte.plotOn(frame, Components(bkg2), LineStyle(kDashed), LineColor(color[1]), LineWidth(4), Range(xMinRange, xMaxRange));
    funcFlatte.plotOn(frame, Components(sigFlatte), LineStyle(kDashed), LineColor(color[2]), LineWidth(4), Range(xMinRange, xMaxRange));
    /*RooFitResult**/ r4 = funcFlatte.fitTo(dh, Save());
    r4->Print("v");
    r4->correlationMatrix().Print();
    //RooArgSet* params = funcVoig.getParameters(x);
    break;
  }
  default: {
    Printf("Invalid method to perform the fit");
    break;
  }
  }

  //chi2=frame->chiSquare("funcBW", "dh", 3);

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
  fitParameters[15] = fitParameters[4] / fitParameters[6];
  fitParameters[16] = fitParameters[4] / TMath::Sqrt(fitParameters[4] + fitParameters[6]);

  ofstream fitOutput;
  fitOutput.open("fitOutput.txt", std::ios_base::app);
  fitOutput << "Fit parameters" << endl;
  fitOutput << "#it{M}_{#pi#pi} (GeV/#it{c^{2}} = "
            << "\t" << fitParameters[0] << "\t #pm\t" << fitParameters[1] << endl;
  fitOutput << "#Gamma (GeV) = "
            << "\t" << fitParameters[2] << "\t #pm\t" << fitParameters[3] << endl;
  if (fitMethod == 3) {
    fitOutput << "#sigma_Voigtian = "
              << "\t" << fitParameters[13] << "\t #pm\t" << fitParameters[14] << endl;
  }
  fitOutput << "N_sig = "
            << "\t" << fitParameters[4] << "\t #pm\t" << fitParameters[5] << endl;
  fitOutput << "N_bkg = "
            << "\t" << fitParameters[6] << "\t #pm\t" << fitParameters[7] << endl;
  fitOutput << "#alpha = "
            << "\t" << fitParameters[9] << "\t #pm\t" << fitParameters[10] << endl;
  fitOutput << "#alpha = "
            << "\t" << fitParameters[11] << "\t #pm\t" << fitParameters[12] << endl;
  fitOutput << "#chi^2 = "
            << "\t" << fitParameters[8] << endl;
  fitOutput << "------------------------------------------------------------------------------------" << endl;
  fitOutput.close();

  return frame;
}

Bool_t isFitFailed(Int_t* fitStatus)
{
  if (!fitStatus) {
    Printf("Invalid array.");
    return kTRUE;
  }
  for (Int_t o = 0; o < 4; o++) {
    if (fitStatus[o] < 0) {
      Printf("Fit status not evaluated.");
      return kTRUE;
    }
    if (fitStatus[o] == 4) {
      Printf("Minuit status = 4. Fit failed");
      return kTRUE;
    } else if (fitStatus[o] == -1)
      Printf("Minuit status unset.");
    else if (fitStatus[o] != 0)
      Printf("Minuit status != 0 or !=4.");
  }
  return 0;
}

RooFitResult* fitChi2(RooChi2Var* chi2Var, Int_t* fitStatus, Bool_t minos = 1, Bool_t improve = 0)
{
  if (!chi2Var)
    return 0x0;
  RooFitResult *temp, *result;
  RooMinuit* model2 = new RooMinuit(*chi2Var);
  model2->setStrategy(2);
  model2->migrad();
  temp = (RooFitResult*)model2->save();
  fitStatus[0] = temp->status();
  model2->hesse();
  temp = (RooFitResult*)model2->save();
  fitStatus[1] = temp->status();

  if (improve) {
    model2->improve();
    temp = (RooFitResult*)model2->save();
    fitStatus[2] = temp->status();
  } else
    fitStatus[2] = 0;

  if (minos) {
    model2->minos();
    temp = (RooFitResult*)model2->save();
    fitStatus[3] = temp->status();
  } else
    fitStatus[3] = 0;

  result = (RooFitResult*)model2->save();
  result->Print("v");
  return result;
}

/*RooFitResult* fitLikelihood(RooNLLVar* nll, Int_t* fitStatus, Bool_t minos = 1, Bool_t improve = 0)
{
  if (!nll)
    return 0x0;
  RooFitResult *temp, *result;
  RooRealVar* offset = new RooRealVar("offset", "offset", -10000000);
  RooAbsReal* L = new RooFormulaVar("L", "L", "L", RooArgSet(*nll, *offset));
  RooMinuit* model1 = new RooMinuit(*L);
  //model1->setPrintLevel(4);
  model1->setStrategy(1);
  model1->setEps(1e-16);
  model1->migrad();
  model1->setStrategy(2);
  model1->migrad();
  temp = (RooFitResult*)model1->save();
  fitStatus[0] = temp->status();
  model1->hesse();
  temp = (RooFitResult*)model1->save();
  fitStatus[1] = temp->status();
  if (improve) {
    model1->improve();
    temp = (RooFitResult*)model1->save();
    fitStatus[2] = temp->status();
  } else
    fitStatus[2] = 0;
  if (minos) {
    model1->minos();
    temp = (RooFitResult*)model1->save();
    fitStatus[3] = temp->status();
  } else
    fitStatus[3] = 0;
  result = (RooFitResult*)model1->save();
  return result;
}*/

TPaveText* textFitResults(Int_t fitMethod, Float_t x1, Float_t y1, Float_t x2, Float_t y2, Double_t* fitParameters)
{
  TPaveText* text = new TPaveText(x1, y1, x2, y2);
  text->SetLabel("Fit Parameters");
  text->SetBorderSize(1);
  text->SetTextAlign(21);
  text->SetFillColor(kWhite);
  text->SetTextFont(42);
  text->SetTextColor(kBlack);
  text->SetTextSize(0.07);
  text->AddText("             ");
  text->AddText(Form("#it{M}_{#pi#pi} (GeV/#it{c^{2}}) = %6.4f #pm %6.4f", fitParameters[0], fitParameters[1]));
  text->AddText(Form("#Gamma (GeV) = %6.4f #pm %6.4f", fitParameters[2], fitParameters[3]));
  if (fitMethod == 3) {
    text->AddText(Form("#sigma_{Voigtian} = %6.4f #pm %6.4f", fitParameters[13], fitParameters[14]));
  }
  text->AddText(Form("N_{sig} =  %8.0f #pm %6.0f \n", fitParameters[4], fitParameters[5]));
  text->AddText(Form("N_{bkg} = %8.0f #pm %6.0f \n", fitParameters[6], fitParameters[7]));
  text->AddText(Form("#it{#alpha} = %6.4f #pm %6.4f", fitParameters[9], fitParameters[10]));
  text->AddText(Form("#it{#beta} = %6.4f #pm %6.4f", fitParameters[11], fitParameters[12]));
  text->AddText(Form("#it{#chi}^{2} = %6.4f", fitParameters[8]));
  return text;
}
