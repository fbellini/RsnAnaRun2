////////////////////////////////////////////////////////////////////
//       f0 signal: fit - 21.09.2017                              //
//                                                                //
//  check after changing the method used to run the fit:          //
//  - output.root name                                            //
//  - markers style, histos style, ecc                            //
//  - depending on the variables, what to print on the TPaveText  //
//  - fit Parameters (evernote - 18_01_29)                        //
//  - summary canvas name                                         //
//  - ...to be continued                                          //
////////////////////////////////////////////////////////////////////

//#include "HistoMakeUp.C"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooGlobalFunc.h"
#include "RooPlot.h"
#include "RooRealVar.h"
//#include "SetStyle.C"
#include "TLegend.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TMath.h"
#include "TRandom.h"
#include "TTree.h"
#include "TLine.h"
#include "TFile.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
//#include "RooComplex.h"
#include "RooNLLVar.h"
//#include "RooFlatte.h"
#include "RooChi2Var.h"
#include "RooRelBW.h"
#include "RooRelBW_f2.h"
#include "TPaveText.h"
#include "RooMinuit.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooVoigtian.h"
#include "RooBreitWigner.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TGaxis.h"
#include <fstream>

using namespace RooFit;

RooPlot* fit(TH1D* h1, Double_t xMinRange, Double_t xMaxRange, Int_t fitMethod, Double_t* fitParameters, Bool_t useChi2, Int_t* fitStatus, Double_t nStart, Double_t nMin, Double_t nMax, Double_t BStart, Double_t BMin, Double_t BMax, Double_t CStart, Double_t CMin, Double_t CMax, Double_t alphaStart, Double_t alphaMin, Double_t alphaMax, Double_t betaStart, Double_t betaMin, Double_t betaMax, Double_t sigmaVoig, Double_t binningPT, Int_t ibin);

Bool_t isFitFailed(Int_t* fitStatus);
RooFitResult* fitChi2(RooChi2Var* chi2Var, Int_t* fitStatus, Bool_t minos, Bool_t improve);
//RooFitResult* fitLikelihood(RooNLLVar* nll, Int_t* fitStatus, Bool_t minos = 1, Bool_t improve = 0);
TPaveText* textFitResults(Int_t fitMethod, Float_t x1, Float_t y1, Float_t x2, Float_t y2, Double_t* fitParameters, Double_t xMinRange, Double_t xMaxRange, Int_t ndf);

void f0_fit_f2_rho(
    Int_t iSelBin=4,
    TString inputFile = "bgSubtraction_RsnOut_f0_2sTPC_3sTOFveto_extended_rebin5_definitivo.root",
    //TString inputFile = "bgSubtraction_RsnOut_f0_3sTPC_3sTOFveto_extended_rebin5_definitivo.root",
    //TString inputFile = "bgSubtraction_RsnOut_f0_2sTPC_4sTOFveto_extended_rebin5_definitivo.root",
    //TString inputFile = "bgSubtraction_RsnOut_f0_2sTPC_extended_rebin5_definitivo.root",
  //  TString inputFile = "bgSubtraction_RsnOut_f0_2sTPC_3sTOFveto_extended_rebin5_definitivo_LHC17p.root",
    //TString inputFile = "bgSubtraction_RsnOut_f0_2sTPC_3sTOFveto_LHC17p.root",
    Bool_t useChi2 = kTRUE,
    Double_t xMinRange = 0.82,
    Double_t xMaxRange = 1.76,
    Int_t ndf = 188,
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

// #ifdef __CINT__
// gROOT->ProcessLine(".x RooRelBW.cxx");
// gROOT->ProcessLine(".x RooComplex.cxx");
// gROOT->ProcessLine(".x RooFlatte.cxx");
// #endif

  TString folderName = Form("fitMethod_%i", fitMethod);
  //gSystem->Exec(Form("mkdir %s", folderName.Data()));
  gROOT->ForceStyle();
  gStyle->SetOptStat(0);
  gStyle->SetLineWidth(1);
  Int_t markerstyle[] = {21, 8, 33, 25};
  gStyle->SetMarkerStyle(markerstyle[fitMethod+1]); //wide square = 21; diamond = 33; full circle = 24; wide triangle up = 25
  //gStyle->SetMarkerColor(kRed+1);
  gStyle->SetMarkerSize(0.5);
  gStyle->SetHistLineWidth(1.5);
  Color_t linecolor[] = { kRed+1, kBlue+1, kGreen+2, kMagenta };
  //gStyle->SetHistLineColor(linecolor[fitMethod]); //BW = kRed+1; relBW = kBlue+1; Voigt = kGreen+2; Flatte' = kMagenta;
  //gStyle->SetMarkerColor(linecolor[fitMethod]);
  gStyle->SetTitleAlign(33);
  gStyle->SetTitleX(.95);
  gStyle->SetTitleY(1.);
  gStyle->SetStatColor(kWhite);
  gStyle->SetStatStyle(1001);
  gStyle->SetStatFontSize(0.05);
  gStyle->SetStatFont(42);

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

  TCanvas* canvas[20];
  RooPlot* plot = 0x0;
  TPaveText* textFit = 0x0;

  TAxis* bins = (TAxis*)fin->Get("bins");
  Int_t numPt = bins->GetNbins();

  const Int_t dimPt = numPt + 1;
  Double_t pT[dimPt];
  for (Int_t z = 0; z < dimPt; z++) {
    pT[z] = bins->GetBinLowEdge(z + 1);
  }

  printf("%d\n", dimPt);

  TFile* fout = TFile::Open("alternative_bg2exp_massFixed_widthFree_LHC15n_080_176.root", "RECREATE");

  //input histograms
  TH1D* hUSPminusLSB[19];
  for (Int_t i = 0; i < 19; i++) {
    hUSPminusLSB[i] = 0x0;
  }

  //define histograms for fit parameters vs pT
  TH1F* hMvsPtF0 = new TH1F("hMvsPtF0", "f_{0} mass; #it{p}_{T} (GeV/#it{c}); M (GeV/#it{c^{2}})", numPt, pT);
  TH1F* hWidthVsPtF0 = new TH1F("hWidthVsPtF0", "f_{0} width; #it{p}_{T} (GeV/#it{c}); #Gamma (GeV)", numPt, pT);
  //TH1F* hWidthVsPtF0 = new TH1F("hWidthVsPtF0", "f_{0} width; #it{p}_{T} (GeV/#it{c}); #Gamma (GeV)", 50, 0.0, 0.1);
  TH1F* hSigVsPtF0 = new TH1F("hSigVsPtF0", "f_{0} raw yield; #it{p}_{T} (GeV/#it{c}); raw yield, dN/d#it{p}_{T}", numPt, pT);

  TH1F* hMvsPtF2 = new TH1F("hMvsPtF2", "f_{2} mass; #it{p}_{T} (GeV/#it{c}); M (GeV/#it{c^{2}})", numPt, pT);
  TH1F* hWidthVsPtF2 = new TH1F("hWidthVsPtF2", "f_{2} width; #it{p}_{T} (GeV/#it{c}); #Gamma (GeV)", numPt, pT);
  TH1F* hSigVsPtF2 = new TH1F("hSigVsPtF2", "f_{2} raw yield; #it{p}_{T} (GeV/#it{c}); raw yield, dN/d#it{p}_{T}", numPt, pT);

  TH1F* hSigmaVoigVsPtF0 = new TH1F("hSigmaVoigVsPtF0", "f_{0} #sigma_{Voigtian}; #it{p}_{T} (GeV/#it{c}); #sigma_{Voigtian}", numPt, pT);


  TH1F* hSigOverBkgVsPtF0 = new TH1F("hSigOverBkgVsPtF0", "f_{0} sig/bkg; #it{p}_{T} (GeV/#it{c}); S/B", numPt, pT);
  TH1F* hSigOverBkgVsPtF2 = new TH1F("hSigOverBkgVsPtF2", "f_{2} sig/bkg; #it{p}_{T} (GeV/#it{c}); S/B", numPt, pT);

  TH1F* hSignificanceVsPtF0 = new TH1F("hSignificanceVsPtF0", "f_{0} sig/#sqrt(sig+bkg); #it{p}_{T} (GeV/#it{c}); S/#sqrt(S+B)", numPt, pT);
  TH1F* hSignificanceVsPtF2 = new TH1F("hSignificanceVsPtF2", "f_{2} sig/#sqrt(sig+bkg); #it{p}_{T} (GeV/#it{c}); S/#sqrt(S+B)", numPt, pT);


  TH1F* hg0VsPt = new TH1F("hg0VsPt", "f_{0} coupling constant to channel 1; #it{p}_{T} (GeV/#it{c}); g_{#pi#pi}", numPt, pT);
  TH1F* hg1VsPt = new TH1F("hg1VsPt", "f_{0} coupling constant to channel 2; #it{p}_{T} (GeV/#it{c}); g_{KK}", numPt, pT);

  TH1F* hChi2OverNDFVsPtF0 = new TH1F("hChi2OverNDFVsPtF0", "#it{#chi}^{2}/n.d.f.; #it{p}_{T} (GeV/#it{c}); #it{#chi}^{2}/n.d.f.", numPt, pT);
  TH1F* hProb = new TH1F("hProb", "fit probability; #it{p}_{T} (GeV/#it{c}); probability", numPt, pT);

  TCanvas* c_frame = new TCanvas("c_frame", "c_frame", 800, 600);

  //Loop on pt bins
  for (Int_t ibin = 0; ibin < 13; ibin++) {
    //if(ibin!=iSelBin) continue;
    hUSPminusLSB[ibin] = (TH1D*)fin->Get(Form("USP-LSBGeoMean_%d", ibin));
    if (!hUSPminusLSB[ibin]) {
      Printf("Input histogram error.");
      return;
    }

    Int_t iMinBinPt = pT[ibin];
    Int_t iMaxBinPt = pT[ibin + 1];

    //define fit parameters and status flags
    Double_t bWidth;
    Double_t fitParameters[30];
    for (Int_t ipar = 0; ipar < 30; ipar++) {
      fitParameters[ipar] = -1.0;
    }
    Double_t sigOverBkgF0 = -1.0;
    Double_t sigOverBkgF2 = -1.0;
    Double_t significanceF0 = -1.0;
    Double_t significanceF2 = -1.0;
    Double_t chi2OverNDF = -9999999.;
    Double_t probability = -9999999.;
    Int_t fitStatus[4] = { -1, -1, -1, -1 };


    Float_t nStart[18] = { 1.,  1.,  1.,   1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1., 1.,  1.,  1./*,  1.,  1.*/};
    Float_t nMin[18] = { 0.,  0.,  0.,  0.,  0.,  1.,  0., 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0., 0.,  0.,  0./*,  0.,  0.*/};
    Float_t nMax[18] = { 2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.,  2./*,  4.,  4.*/};

    Float_t BStart[18] = { 3.,  3.,  3.,  3.,  3.,  3.,  3.,  3.,  3.,  3.,  3.,  3.,  3.,  3.,  3., 3.,  3.,  3./*,  3.,  3.*/};
    Float_t BMin[18] = { 0.,  0.,  0.,  0.,  0.,  1.,  1., 1.,  1.,  1.,  1.,  1.,  0.,  0.,  0., 0.,  0.,  0./*,  0.,  0.*/};
    Float_t BMax[18] = { 5.,  4.,  5.,  4.,   4.,  4.,  4.,  4.,  5.,  5.,  5.,  4.,  5.,  5.,  5.,  5.,  5.,  5./*,  5.,  5.*/};

    Float_t CStart[18] = {2.,  2.,  2.,  2.,  1.,  1.,  1.,  2.,  1.,  1.,  1.,  1.,  1.,  1.,  1., 1.,  1.,  1./*,  1.,  1.*/};
    Float_t CMin[18] = { 0.,  0.,  0.,  0.,  0.,  0.,  0., 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0., 0.,  0.,  0./*,  0.,  0.*/};
    //Float_t CMax[19] = { 10.,  3.,  10.,  10.,  10.,  10.,  10.,  10.,  10.,  10.,  10.,  10.,  10.,  10.,  10.,  10.,  10.,  10.,  10./*,  5.,  5.*/};
    Float_t CMax[18] = { 5.,  5.,  5.,  5.,   3.,  3.,  5.,  5.,  5.,  5.,  4.,  4.,  5.,  5.,  5.,  5.,  5.,  5./*,  5.,  5.*/};


    Float_t alphaStart[18] = { -9.,   -9.,   -5.,   -8.,   -8.,   -8.,  -8.,   -4.,  -4., -9.,   -8.,   -5.,   -8.,   -8.,   -8.,  -8.,   -4.,  -4.};
    Float_t alphaMin[18] = {  -14.,  -14.,  -10.,  -10.,  -10.,   -10., -10.,   -7.,  -8., -14.,  -10.,  -10.,  -10.,  -10.,   -8., -10.,   -7.,  -8.};
    Float_t alphaMax[18] = {   -1.,   -1.,   -3.,   -1.,   -1.,   -2.,  -2.,    0.,   0.,  -4.,   -3.,   -3.,   -1.,   -1.,   -2.,  -2.,    0.,   0.};

    Float_t betaStart[18] = {   0.,   -2.,   -2.,    2.,    2.,    2.,   2.,    2.,   0., 0.,   -2.,   -2.,    2.,    2.,    2.,   2.,    2.,   0.};
    Float_t betaMin[18] = {    -1.,   -4.,   -4.,    0.,    0.,    0.,   0.,    0.,  -1.,  -1.,   -4.,   -4.,    0.,    0.,    0.,   0.,    0.,  -1.};
    Float_t betaMax[18] = {     1.,    4.,    4.,    4.,    4.,    4.,   4.,    4.,   1., 1.,    4.,    4.,    4.,    4.,    4.,   4.,    4.,   1.};

    Float_t sigmaVoig[14] = {0.005, 0.005, 0.005, 0.005, 0.006, 0.006, 0.006, 0.006, 0.007, 0.007, 0.007, 0.008, 0.008, 0.008};

    Double_t binningPT[14] = {0., 0.3, 0.6, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0};

    //run fit
   plot = (RooPlot*)fit(hUSPminusLSB[ibin], xMinRange, xMaxRange, fitMethod, fitParameters, kTRUE, fitStatus, nStart[ibin], nMin[ibin], nMax[ibin], BStart[ibin], BMin[ibin], BMax[ibin], CStart[ibin], CMin[ibin], CMax[ibin], alphaStart[ibin], alphaMin[ibin], alphaMax[ibin], betaStart[ibin], betaMin[ibin], betaMax[ibin], sigmaVoig[ibin], pT[ibin], ibin);

   c_frame->cd();
   hUSPminusLSB[ibin]->Draw();
   hUSPminusLSB[ibin]->GetYaxis()->SetRangeUser(0.,hUSPminusLSB[ibin]->GetMaximum()*1.3);
   plot->Draw("same");

   //check fit status and save fit output if fit succeeded
   if (!isFitFailed(fitStatus)) {

      textFit = (TPaveText*)textFitResults(fitMethod, 0.05, 0.1, 0.95, 0.8, fitParameters, xMinRange, xMaxRange, ndf);

      canvas[ibin] = new TCanvas(Form("c%d", ibin), Form("canvas %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin + 1]), 800, 800);
      //canvas[ibin]->Divide(2, 1);
      canvas[ibin]->cd();
      plot->GetXaxis()->SetTitle("#it{M}_{#pi#pi}(GeV/#it{c^{2}})");
      plot->Draw();
      /*TLegend *leg1 = new TLegend(0.65,0.73,0.86,0.87);
      leg1->SetFillColor(kWhite);
      leg1->SetLineColor(kWhite);
      leg1->AddEntry("dh", "Data", "lpf");
      leg1->AddEntry("model2", "Signal + background", "lpf");
      leg1->AddEntry("sigRelBW", "f0(980)", "lpf");
      leg1->AddEntry("sigf2", "f2(1270)", "lpf");
      leg1->AddEntry("bkg", "exp(a)", "lpf");
      leg1->Draw();*/
    //canvas[ibin]->cd(2);
      //textFit->Draw();
      //canvas[ibin]->Print(Form("fit_method%d_%2.1f<pT<%2.1f_mFree_g0ANDg1Fix.png", fitMethod, pT[ibin], pT[ibin + 1]));

      //protection: check that fitParameters exists
      if (fitParameters) {

        //get significance, S/B and chi2/ndf
        sigOverBkgF0 = fitParameters[4] / fitParameters[6];
        significanceF0 = fitParameters[4] / TMath::Sqrt(fitParameters[4] + fitParameters[6]);
        chi2OverNDF = fitParameters[8]/ndf;
        probability = TMath::Prob(fitParameters[8], ndf);
        //sigOverBkgF2 = fitParameters[4] / fitParameters[6];
        //significanceF2 = fitParameters[4] / TMath::Sqrt(fitParameters[4] + fitParameters[6]);

        //fill histograms
        hMvsPtF0->SetBinContent(ibin + 1, fitParameters[0]);
        hMvsPtF0->SetBinError(ibin + 1, fitParameters[1]);
        hMvsPtF2->SetBinContent(ibin + 1, fitParameters[17]);
        hMvsPtF2->SetBinError(ibin + 1, fitParameters[18]);

        hWidthVsPtF0->SetBinContent(ibin + 1, fitParameters[2]);
        hWidthVsPtF0->SetBinError(ibin + 1, fitParameters[3]);
        //hWidthVsPtF0->Fill(fitParameters[2]);
        hWidthVsPtF2->SetBinContent(ibin + 1, fitParameters[19]);
        hWidthVsPtF2->SetBinError(ibin + 1, fitParameters[20]);

        bWidth = bins->GetBinWidth(ibin + 1);
        hSigVsPtF0->SetBinContent(ibin + 1, fitParameters[4] / bWidth);
        hSigVsPtF0->SetBinError(ibin + 1, fitParameters[5] / bWidth);
        hSigVsPtF2->SetBinContent(ibin + 1, fitParameters[25] / bWidth);
        hSigVsPtF2->SetBinError(ibin + 1, fitParameters[26] / bWidth);

        if (fitMethod == 3) {
          hSigmaVoigVsPtF0->SetBinContent(ibin + 1, fitParameters[13]);
          hSigmaVoigVsPtF0->SetBinError(ibin + 1, fitParameters[14]);
        }

        hSigOverBkgVsPtF0->SetBinContent(ibin + 1, sigOverBkgF0);
        hSigOverBkgVsPtF0->SetBinError(ibin + 1, sigOverBkgF0);
        //hSigOverBkgVsPtF2->SetBinContent(ibin + 1, sigOverBkgF2);
        //hSigOverBkgVsPtF2->SetBinError(ibin + 1, sigOverBkgF2);
        hSignificanceVsPtF0->SetBinContent(ibin + 1, significanceF0);
        hSignificanceVsPtF0->SetBinError(ibin + 1, significanceF0);
        //hSignificanceVsPtF2->SetBinContent(ibin + 1, significanceF2);
        //hSignificanceVsPtF2->SetBinError(ibin + 1, significanceF2);
        hChi2OverNDFVsPtF0->SetBinContent(ibin+1, chi2OverNDF);
        hChi2OverNDFVsPtF0->SetBinError(ibin+1, chi2OverNDF);
        hProb->SetBinContent(ibin+1, probability);
        hProb->SetBinError(ibin+1, probability);


        if (fitMethod == 4){
          hg0VsPt->SetBinContent(ibin + 1, fitParameters[21]);
          hg0VsPt->SetBinError(ibin + 1, fitParameters[22]);
          hg1VsPt->SetBinContent(ibin + 1, fitParameters[23]);
          hg1VsPt->SetBinError(ibin + 1, fitParameters[24]);
        }

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
  Double_t yMin2 = 0.01;
  Double_t yMax2 = 0.11;

  //------------------------------------------------------------------------------------------//

  //plot fitted parameters vs pt
  TCanvas* c_histos1 = new TCanvas("c_histos1", "c_histos1", 1200, 600);
  c_histos1->Divide(4, 2);
  c_histos1->cd(1);
  hMvsPtF0->GetYaxis()->SetRangeUser(yMin, yMax);
  hMvsPtF0->Draw();
  TLine* lineMass1 = new TLine(0.5, 0.96, 9., 0.96);
  lineMass1->SetLineColor(kGray + 1);
  lineMass1->SetLineStyle(3);
  lineMass1->Draw();
  TLine* lineMass2 = new TLine(0.5, 1.01, 9., 1.01);
  lineMass2->SetLineColor(kGray + 1);
  lineMass2->SetLineStyle(3);
  lineMass2->Draw();
  c_histos1->cd(2);
  hWidthVsPtF0->GetYaxis()->SetRangeUser(yMin2, yMax2);
  hWidthVsPtF0->Draw();
  TLine* lineWidth1 = new TLine(0.5, 0.01, 9., 0.01);
  lineWidth1->SetLineColor(kGray + 1);
  lineWidth1->SetLineStyle(3);
  lineWidth1->Draw();
  TLine* lineWidth2 = new TLine(0.5, 0.1, 9., 0.1);
  lineWidth2->SetLineColor(kGray + 1);
  lineWidth2->SetLineStyle(3);
  lineWidth2->Draw();
  c_histos1->cd(3);
  gPad->SetLogy();
  hSigVsPtF0->Draw();
  c_histos1->cd(4);
  hSigmaVoigVsPtF0->Draw();
  c_histos1->cd(5);
  hSigOverBkgVsPtF0->Draw("hist");
  c_histos1->cd(6);
  hSignificanceVsPtF0->Draw("hist");
  c_histos1->cd(7);
  hChi2OverNDFVsPtF0->Draw("hist");
  c_histos1->cd(8);
  hProb->Draw("hist");
  //c_histos1->cd(7);
  //hg0VsPt->Draw();
  //c_histos1->cd(8);
  //hg1VsPt->Draw();
  //c_histos1->Print("fitParameters.png");

  TCanvas* c_histos2 = new TCanvas("c_histos2", "c_histos2", 1200, 600);
  c_histos2->Divide(3, 1);
  c_histos2->cd(1);
  //hMvsPtF2->GetYaxis()->SetRangeUser(yMin, yMax);
  hMvsPtF2->Draw();
  c_histos2->cd(2);
  //hWidthVsPtF2->GetYaxis()->SetRangeUser(yMin2, yMax2);
  hWidthVsPtF2->Draw();
  c_histos2->cd(3);
  gPad->SetLogy();
  hSigVsPtF2->Draw();
  //c_histos2->Print("fitParameters2.png");

  //save output to file
  fout->cd();
  hMvsPtF0->Write();
  hWidthVsPtF0->Write();
  hSigVsPtF0->Write();
  hSigmaVoigVsPtF0->Write();
  hSigOverBkgVsPtF0->Write();
  hSignificanceVsPtF0->Write();
  hChi2OverNDFVsPtF0->Write();
  hProb->Write();
  //hg0VsPtF0->Write();
  //hg1VsPtF0->Write();
  hMvsPtF2->Write();
  hWidthVsPtF2->Write();
  hSigVsPtF2->Write();
  bins->Write();
  //fout->Close();

  //return;
}

RooPlot* fit(TH1D* h1, Double_t xMinRange, Double_t xMaxRange, Int_t fitMethod, Double_t* fitParameters, Bool_t useChi2, Int_t* fitStatus, Double_t nStart, Double_t nMin, Double_t nMax, Double_t BStart, Double_t BMin, Double_t BMax, Double_t CStart, Double_t CMin, Double_t CMax, Double_t alphaStart, Double_t alphaMin, Double_t alphaMax, Double_t betaStart, Double_t betaMin, Double_t betaMax, Double_t sigmaVoig, Double_t binningPT, Int_t ibin)
{
  Color_t color[] = { kBlack, kRed, kGreen + 1, kMagenta, kBlue, kOrange};
  const Int_t minBinInvMass = h1->GetXaxis()->FindBin(xMinRange);
  const Int_t maxBinInvMass = h1->GetXaxis()->FindBin(xMaxRange);
  Double_t histo_integral = h1->Integral(minBinInvMass, maxBinInvMass);
  RooRealVar x("x", "x", xMinRange, xMaxRange);
  RooDataHist dh("dh", "dh", RooArgList(x), Import(*h1));
  RooPlot* frame = x.frame(Title(h1->GetTitle()));
  dh.plotOn(frame, Name("dh"), DataError(RooAbsData::SumW2)); ////Name("dh") allows to get chi2
  //dh.statOn(frame,Layout(0.55,0.99,0.8));

  // f0 signal
  RooRealVar mF0("mF0", "mF0", 0.97/*, 0.96, 1.01*/);//0.98, 0.958, 1.01);
  RooRealVar widthF0("widthF0", "widthF0", 0.04/*, 0.01, 0.1*/);//c3Min, c3Max);
  //RooRealVar sigmaF0("sigmaF0", "sigmaF0", 0.003/*, 0.001, 0.01*/);
  RooRealVar g0("g0", "g0", 0.199);
  RooRealVar g1("g1", "g1", 0.597);
  RooRealVar m0a("m0a", "m0a", 0.1396);
  RooRealVar m0b("m0b", "m0b", 0.1396);
  RooRealVar m1a("m1a", "m1a", 0.4937);
  RooRealVar m1b("m1b", "m1b", 0.4937);
  RooRealVar n("n", "n", nStart, nMin, nMax);
  RooRealVar B("B", "B", BStart, BMin, BMax);
  RooRealVar C("C", "C", CStart, CMin, CMax);

  RooBreitWigner sigBW("sigBW", "sigBW", x, mF0, widthF0);
  RooRelBW sigRelBW("relBW", "relBW", x, mF0, widthF0);
  RooVoigtian sigVoig("sigVoig", "sigVoig", x, mF0, widthF0, sigmaVoig, kFALSE);
  //RooFlatte sigFlatte("sigFlatte", "sigFlatte", x, mF0, g0, m0a, m0b, g1, m1a, m1b);

  // bkg
  RooMB bkg("bkg", "bkg", x, n, B, C);
  //RooRealVar alpha("alpha", "alpha", alphaStart, alphaMin, alphaMax);
  //RooRealVar beta("beta", "beta", betaStart, betaMin, betaMax);
  //RooExponential bkg3("bkg3", "Background 1", x, alpha);
  //RooExponential bkg2("bkg2", "Background 2", x, beta);
  //RooRealVar bkg4("bkg4", "background", 0.5, 0., 1.);
  //RooAddPdf bkg3("bkg3", "bkg3", RooArgList(bkg1, bkg2), bkg4);
  RooRealVar bkgFrac("bkgFrac", "fraction of background", 0.8, 0., 1.);

  /*//////////////////////////////////////////////////////////////////////////////////////////////////////
  //ROBA INUTILE
  RooRealVar alpha("alpha", "alpha", alphaStart, alphaMin, alphaMax);
  RooRealVar beta("beta", "beta", betaStart, betaMin, betaMax);
  //RooExponential bkg2("bkg2", "Background 2", x, beta);
  RooRealVar coef_x1("coef_x1", "coef_x0", c1Start, c1Min, c1Max);
  RooRealVar coef_x2("coef_x2", "coef_x2", c2Start, c2Min, c2Max);
  //RooRealVar coef_x3("coef_x3", "coef_x3", c3Start, c3Min, c3Max);
  RooPolynomial bkgdpoly("bkgdpoly","bkgdpoly", x ,RooArgList(coef_x1, coef_x2/*, coef_x3));
  //RooAddPdf sumbkg("sumbkg", "Background", RooArgList(sigf2, sigRho, bkg, bkg2), bkgFrac);

  //////////////////////////////////////////////////////////////////////////////////////////////////////*/

  RooRealVar mF2("mF2", "mF2", 1.2755/*, 1.2747, 1.2763*/);
  RooRealVar widthF2("widthF2", "widthF2", 0.1867/*, 0.1845, 0.1889*/);
  RooRelBW_f2 sigf2("sigf2", "sigf2", x, mF2, widthF2);

  RooRealVar mRho("mRho", "mRho", 0.7753/*, 0.775, 0.7756*/);
  RooRealVar widthRho("widthRho", "widthRho", 0.1491/*, 0.1483, 0.1499*/);
  RooRelBW_rho sigRho("sigRho", "sigRho", x, mRho, widthRho);

  /*RooRealVar mOmega3("mOmega3", "mOmega3", 0.55, 0.5, 0.6);
  RooRealVar widthOmega3("widthOmega3", "widthOmega3", 0.2, 0.1483, 0.1499);
  RooRelBW_rho sigOmega3("sigOmega3", "sigOmega3", x, mOmega3, widthOmega3);
  */

  RooRealVar meanGaus("meanGaus", "meanGaus", 1.15/*, 1.1, 1.2*/);
  RooRealVar sigmaGaus("sigmaGaus", "sigmaGaus", 0.005, 0.001, 0.1);
  RooGaussian sigNearf2("sigNearf2", "sigNearf2", x, meanGaus, sigmaGaus);

  RooRealVar mKStar("mKStar", "mKStar", 0.9/*, 0.89, 1.*/);
  RooRealVar widthKStar("widthKStar", "widthKStar", 0.05);
  RooBreitWigner sigKStar("sigKStar", "sigKStar", x, mKStar, widthKStar);

  RooRealVar funcFrac("funcFrac", "fraction of functions", 0.5, 0., 1.);
  RooRealVar funcFrac2("funcFrac2", "fraction of functions2", 0.5, 0., 1.);
  RooRealVar funcFrac3("funcFrac3", "fraction of functions3", 0.8, 0., 1.);

  //RooAddPdf sumfunc2("sumfunc2", "sumfunc2", RooArgList(sigf2, sigNearf2), funcFrac2);
  //RooAddPdf sumfunc3("sumfunc3", "sumfunc3", RooArgList(sigf2, sigKStar), funcFrac3);
  RooAddPdf sumfunc("sumfunc", "sumfunc", RooArgList(sigf2, sigRho), funcFrac);
  RooAddPdf sumbkg("sumbkg", "Background", RooArgList(sumfunc, bkg), bkgFrac);


  //model for signal + background for different p.d.f.
  RooRealVar nsig("nsig", "signal fraction", histo_integral * 0.5, 0., histo_integral);
  RooRealVar nbkg("nbkg", "background fraction", histo_integral * 0.5, 0., histo_integral);
  RooAddPdf funcBW("model1", "sig+bgBW", RooArgList(sigBW, sumbkg), RooArgList(nsig, nbkg));
  RooAddPdf funcRelBW("model2", "sig+bgRelBW", RooArgList(sigRelBW, sumbkg), RooArgList(nsig, nbkg));
  RooAddPdf funcVoig("model3", "sig+bgVoig", RooArgList(sigVoig, sumbkg), RooArgList(nsig, nbkg));
  //RooAddPdf funcFlatte("model4", "sig+bgFlatte", RooArgList(sigFlatte, sumbkg), RooArgList(nsig, nbkg));

  RooRealVar nsigF2("nsigF2", "signal fraction F2", histo_integral * 0.5, 0., histo_integral);
  RooRealVar nbkgF2("nBKGF2", "background fraction F2", histo_integral * 0.5, 0., histo_integral);
  RooAddPdf funcf2("funcf2", "funcf2", RooArgList(sigf2, sumbkg), RooArgList(nsigF2, nbkgF2));

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
    funcBW.plotOn(frame, Components(sigKStar), LineStyle(kDashed), LineColor(color[1]), LineWidth(2), Range(xMinRange, xMaxRange));
    funcBW.plotOn(frame, Components(sigRho), LineStyle(kDashed), LineColor(color[2]), LineWidth(2), Range(xMinRange, xMaxRange));
    funcBW.plotOn(frame, Components(sigf2), LineStyle(kDashed), LineColor(color[3]), LineWidth(2), Range(xMinRange, xMaxRange));
    funcBW.plotOn(frame, Components(sigNearf2), LineStyle(kDashed), LineColor(color[4]), LineWidth(2), Range(xMinRange, xMaxRange));
    funcBW.plotOn(frame, Components(bkg), LineStyle(kDashed), LineColor(color[5]), LineWidth(2), Range(xMinRange, xMaxRange));
    funcBW.plotOn(frame, Components(sigBW), Name("funcBW"), LineStyle(kDashed), LineColor(color[0]), LineWidth(3), Range(xMinRange, xMaxRange));

    /*RooFitResult* */ r1 = funcBW.fitTo(dh, Save());
    r1 = funcBW.fitTo(dh, Save());
    r1->Print("v");
    r1->correlationMatrix().Print();
    Double_t chi2 = frame->chiSquare("funcBW", "dh", 3);
    fitParameters[8] = chi2;
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
    funcRelBW.plotOn(frame, Components(sigKStar), LineStyle(kDashed), LineColor(color[1]), LineWidth(2), Range(xMinRange, xMaxRange));
    funcRelBW.plotOn(frame, Components(sigRho), LineStyle(kDashed), LineColor(color[2]), LineWidth(2), Range(xMinRange, xMaxRange));
    funcRelBW.plotOn(frame, Components(sigf2), LineStyle(kDashed), LineColor(color[3]), LineWidth(2), Range(xMinRange, xMaxRange));
    funcRelBW.plotOn(frame, Components(sigNearf2), LineStyle(kDashed), LineColor(color[4]), LineWidth(2), Range(xMinRange, xMaxRange));
    funcRelBW.plotOn(frame, Components(bkg), LineStyle(kDashed), LineColor(color[5]), LineWidth(2), Range(xMinRange, xMaxRange));
    funcRelBW.plotOn(frame, Components(sigRelBW), Name("funcRelBW"), LineStyle(kDashed), LineColor(color[0]), LineWidth(3), Range(xMinRange, xMaxRange));

    /*RooFitResult* */ r2 = funcRelBW.fitTo(dh, Save());
    r2 = funcRelBW.fitTo(dh, Save());
    r2->Print("v");
    r2->correlationMatrix().Print();
    Double_t chi2 = frame->chiSquare("funcRelBW", "dh", 3);
    fitParameters[8] = chi2;
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
    funcVoig.plotOn(frame, Components(sigKStar), LineStyle(kDashed), LineColor(color[1]), LineWidth(2), Range(xMinRange, xMaxRange));
    funcVoig.plotOn(frame, Components(sigRho), LineStyle(kDashed), LineColor(color[2]), LineWidth(2), Range(xMinRange, xMaxRange));
    funcVoig.plotOn(frame, Components(sigf2), LineStyle(kDashed), LineColor(color[3]), LineWidth(2), Range(xMinRange, xMaxRange));
    funcVoig.plotOn(frame, Components(sigNearf2), LineStyle(kDashed), LineColor(color[4]), LineWidth(2), Range(xMinRange, xMaxRange));
    funcVoig.plotOn(frame, Components(bkg), LineStyle(kDashed), LineColor(color[5]), LineWidth(2), Range(xMinRange, xMaxRange));
    funcVoig.plotOn(frame, Components(sigVoig), Name("funcVoig"), LineStyle(kDashed), LineColor(color[0]), LineWidth(3), Range(xMinRange, xMaxRange));


    /*RooFitResult**/ r3 = funcVoig.fitTo(dh, Save());
    r3 = funcVoig.fitTo(dh, Save());
    r3->Print("v");
    r3->correlationMatrix().Print();
    Double_t chi2 = frame->chiSquare("funcVoig", "dh", 3);
    fitParameters[8] = chi2;
    Double_t sigmaVoig = sigmaF0.getVal();
    Double_t sigmaVoigErr = sigmaF0.getError();
    fitParameters[13] = sigmaVoig;
    fitParameters[14] = sigmaVoigErr;

    //RooArgSet* params = funcVoig.getParameters(x);
    break;
  }
  /*case 4: {
    RooChi2Var* chi2Var = new RooChi2Var("chi2Var", "chi2Var", funcFlatte, dh, kTRUE);
    RooNLLVar* nll = new RooNLLVar("nll", "nll", funcFlatte, dh);
    if (useChi2)
      r4 = (RooFitResult*)fitChi2(chi2Var, fitStatus, minos, improve);

    funcFlatte.plotOn(frame);
    funcFlatte.paramOn(frame, Layout(0.4));
    funcFlatte.plotOn(frame, Components(bkg), LineStyle(kDashed), LineColor(color[0]), LineWidth(4), Range(xMinRange, xMaxRange));
    funcFlatte.plotOn(frame, Components(sigf2), LineStyle(kDashed), LineColor(color[3]), LineWidth(4), Range(xMinRange, xMaxRange));
    funcFlatte.plotOn(frame, Components(sigFlatte), Name("funcFlatte"), LineStyle(kDashed), LineColor(color[2]), LineWidth(4), Range(xMinRange, xMaxRange));
  //  RooFitResult*
  r4 = funcFlatte.fitTo(dh, Save());
    r4->Print("v");
    r4->correlationMatrix().Print();
    Double_t chi2 = frame->chiSquare("funcFlatte", "dh", 3);
    fitParameters[8] = chi2;
    Double_t g_pipi = g0.getVal();
    Double_t g_pipiErr = g0.getError();
    Double_t g_KK = g1.getVal();
    Double_t g_KKErr = g1.getError();
    fitParameters[21] = g_pipi;
    fitParameters[22] = g_pipiErr;
    fitParameters[23] = g_KK;
    fitParameters[24] = g_KKErr;
    //RooArgSet* params = funcVoig.getParameters(x);
    break;
  }*/
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
  //Double_t alphaVal = alpha.getVal();
  //Double_t alphaValErr = alpha.getError();
  //Double_t betaVal = beta.getVal();
  //Double_t betaValErr = beta.getError();
  Double_t signalMassF2 = mF2.getVal();
  Double_t signalMassErrF2 = mF2.getError();
  Double_t signalWidthF2 = widthF2.getVal();
  Double_t signalWidthErrF2 = widthF2.getError();
  Double_t signalMassRho = mRho.getVal();
  Double_t signalMassErrRho = mRho.getError();
  Double_t signalWidthRho = widthRho.getVal();
  Double_t signalWidthErrRho = widthRho.getError();
  Double_t signalMassKStar = mKStar.getVal();
  Double_t signalMassErrKStar = mKStar.getError();
  Double_t signalWidthKStar = widthKStar.getVal();
  Double_t signalWidthErrKStar = widthKStar.getError();

  Double_t nSignalF2 = nsigF2.getVal();
  Double_t nSignalErrF2 = nsigF2.getError();
  Double_t nBkgF2 = nbkgF2.getVal();
  Double_t nBkgErrF2 = nbkgF2.getError();


  fitParameters[0] = signalMass;
  fitParameters[1] = signalMassErr;
  fitParameters[2] = signalWidth;
  fitParameters[3] = signalWidthErr;
  fitParameters[4] = nSignal;
  fitParameters[5] = nSignalErr;
  fitParameters[6] = nBkg;
  fitParameters[7] = nBkgErr;
  //fitParameters[9] = alphaVal;
  //fitParameters[10] = alphaValErr;
  //fitParameters[11] = betaVal;
  //fitParameters[12] = betaValErr;
  fitParameters[15] = fitParameters[4] / fitParameters[6];
  fitParameters[16] = fitParameters[4] / TMath::Sqrt(fitParameters[4] + fitParameters[6]);
  fitParameters[17] = signalMassF2;
  fitParameters[18] = signalMassErrF2;
  fitParameters[19] = signalWidthF2;
  fitParameters[20] = signalWidthErrF2;
  fitParameters[25] = nSignalF2;
  fitParameters[26] = nSignalErrF2;
  fitParameters[27] = signalMassRho;
  fitParameters[28] = signalMassErrRho;
  fitParameters[29] = signalMassKStar;
  //fitParameters[30] = signalMassErrKStar;
  //fitParameters[26] = sigf2Err;



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
  if (fitMethod == 4) {
    fitOutput << "#g_{#pi#pi} = "
              << "\t" << fitParameters[21] << "\t #pm\t" << fitParameters[22] << endl;
    fitOutput << "#g_{KK} = "
              << "\t" << fitParameters[23] << "\t #pm\t" << fitParameters[24] << endl;
  }
  fitOutput << "N_sig = "
            << "\t" << fitParameters[4] << "\t #pm\t" << fitParameters[5] << endl;
  fitOutput << "N_bkg = "
            << "\t" << fitParameters[6] << "\t #pm\t" << fitParameters[7] << endl;
  //fitOutput << "#alpha = "
  //          << "\t" << fitParameters[9] << "\t #pm\t" << fitParameters[10] << endl;
  //fitOutput << "#alpha = "
  //          << "\t" << fitParameters[11] << "\t #pm\t" << fitParameters[12] << endl;
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

TPaveText* textFitResults(Int_t fitMethod, Float_t x1, Float_t y1, Float_t x2, Float_t y2, Double_t* fitParameters, Double_t xMinRange, Double_t xMaxRange, Int_t ndf)
{
  TPaveText* text = new TPaveText(x1, y1, x2, y2);
  text->SetLabel("Fit Parameters");
  text->SetBorderSize(1);
  text->SetTextAlign(21);
  text->SetFillColor(kWhite);
  text->SetTextFont(42);
  text->SetTextColor(kBlack);
  text->SetTextSize(0.07);
  text->SetBorderSize(0);
  text->AddText("             ");
  text->AddText(Form("#it{M}_{f_{0}} (GeV/#it{c^{2}}) = %6.4f (fixed)", fitParameters[0]));
  text->AddText(Form("#Gamma_{f_{0}} (GeV/#it{c^{2}}) = %6.4f #pm %6.4f", fitParameters[2], fitParameters[3]));
  if (fitMethod == 3) {
    text->AddText(Form("#sigma_{Voigtian} = %6.4f #pm %6.4f", fitParameters[13], fitParameters[14]));
  }
  /*if (fitMethod == 4) {
    text->AddText(Form("g_{#pi#pi} = %6.4f #pm %6.4f", fitParameters[21], fitParameters[22]));
    text->AddText(Form("g_{KK} = %6.4f #pm %6.4f", fitParameters[23], fitParameters[24]));
  }*/
  text->AddText(Form("N_{sig} =  %8.0f #pm %6.0f \n", fitParameters[4], fitParameters[5]));
  text->AddText(Form("N_{bkg} = %8.0f #pm %6.0f \n", fitParameters[6], fitParameters[7]));
  //text->AddText(Form("#it{#alpha} = %6.4f #pm %6.4f", fitParameters[9], fitParameters[10]));
  //text->AddText(Form("#it{#beta} = %6.4f #pm %6.4f", fitParameters[11], fitParameters[12]));
  text->AddText(Form("#it{M}_{f_{2}} (GeV/#it{c^{2}}) = %6.4f (fixed)", fitParameters[17]/*, fitParameters[18]*/));
  text->AddText(Form("#it{#Gamma}_{f_{2}} (GeV/#it{c^{2}}) = %6.4f (fixed)", fitParameters[19]/*, fitParameters[20]*/));
  text->AddText(Form("#it{M}_{#rho} (GeV/#it{c^{2}}) = %6.4f (fixed)", fitParameters[17]/*, fitParameters[18]*/));

  text->AddText(Form("#it{#chi}^{2}/n.d.f. = %6.4f", fitParameters[8]/ndf));
  text->AddText(Form("fit range: %3.2f - %3.2f", xMinRange, xMaxRange));
  return text;
}
