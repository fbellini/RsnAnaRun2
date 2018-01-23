//////////////////////////////////////////////////////////////////////////////////////////////////////
//effciencies ratios - 30.11.2017                                                                   //
//                                                                                                  //
//before using the macro, you need HistoMakeUp.C and SetStyle.C in the same directory. you also     //
//need the following lines in a rootlogon.C                                                         //
//{                                                                                                 //
//gROOT->LoadMacro("HistoMakeUp.C+g");                                                              //
//gROOT->LoadMacro("SetStyle.C+g");                                                                 //
//SetStyle();                                                                                       //
//gROOT->LoadMacro("f0_eff.C+g");                                                                   //
//}                                                                                                 //
//                                                                                                  //
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
#include "f0Params.h"

void f0_eff()
{
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat("1111");
  gStyle->SetTextFont(42);

  TFile* finAnalysis = TFile::Open("AnalysisResults.root");
  TDirectory* directory1 = (TDirectory*)finAnalysis->Get("F0_bkg");
  directory1->cd();
  TDirectory* directory2 = (TDirectory*)directory1->Get("MyOutputContainer");
  directory2->cd();

  TFile* fout = TFile::Open("efficiencies.root", "RECREATE");

  Double_t pT[] = {0.5, 1., 1.5, 2., 2.5, 3., 4., 5., 7., 9., 11.};
  Int_t   npT  = sizeof(pT) / sizeof(pT[0]) - 1;
  TAxis *pTbins = new TAxis(npT, pT);

  enum PartSpec { kF0 = 0,
   kOmega = 1,
   kRho = 2,
   kEta = 3,
   kEtaPr = 4,
   kF1 = 5,
   kF2 = 6,
   kKStar = 7,
   kK0s = 8,
   kPhi = 9,
   kAll = 10 };

  const Char_t* particles[kAll] = { "f0", "omega", "rho", "eta", "etaPr", "f1", "f2", "kStar", "k0s", "phi" };

  TH2F* hGenMvsPt[kAll] = { 0x0 };
  TH2F* hRecoMvsPt[kAll] = { 0x0 };
  TH1D* hGenPt[kAll] = { 0x0 };
  TH1D* hRecoPt[kAll] = { 0x0 };
  TH1D* hEffVsPt[kAll] = { 0x0 };

  TCanvas* canvas[kAll];

  Double_t all[kAll][11], trueReco[kAll][11], efficiency[kAll][11], errEfficiency[kAll][11];
  Double_t iMinBin1[11], iMinBin2[11], iMaxBin1[11], iMaxBin2[11];

  //Double_t all[kAll][npT+1], trueReco[kAll][npT+1], efficiency[kAll][npT+1], errAll[kAll][npT+1], errTrueReco[kAll][npT+1], errEfficiency[kAll][npT+1];

  for (Int_t i = 0; i < kAll; i++) {
    hGenMvsPt[i] = (TH2F*)directory2->FindObject(Form("fGenMassVsPt_%s", particles[i]));
    hRecoMvsPt[i] = (TH2F*)directory2->FindObject(Form("fRecoMassVsPt_%s", particles[i]));
    hGenPt[i] = (TH1D*)hGenMvsPt[i]->ProjectionY(Form("hGenPt%s", particles[i]));
    hRecoPt[i] = (TH1D*)hRecoMvsPt[i]->ProjectionY(Form("hRecoPt%s", particles[i]));
    hEffVsPt[i] = new TH1D(Form("hEffVsPt_%s", particles[i]), Form("%s efficiency in pp@5.02 TeV; p_{T} (GeV/#it{c}); #epsilon = reco true / generated", particles[i]), npT, pT);


    for(Int_t j=0; j<npT+1; j++){
      iMinBin1[j] = hGenPt[0]->GetXaxis()->FindBin(pT[j]);
      iMaxBin1[j] = hGenPt[0]->GetXaxis()->FindBin(pT[j+1]);
      iMinBin2[j] = hRecoPt[0]->GetXaxis()->FindBin(pT[j]);
      iMaxBin2[j] = hRecoPt[0]->GetXaxis()->FindBin(pT[j+1]);

      if (iMinBin1[j] != iMinBin2[j] || iMaxBin1[j] != iMaxBin2[j])
      ::Fatal("Error in rebinning histos!", "Error in rebinning histos!");

      all[0][j] = hGenPt[i]->Integral(iMinBin1[j], iMaxBin1[j]);
      trueReco[0][j] = hRecoPt[i]->Integral(iMinBin2[j], iMaxBin2[j]);
      efficiency[0][j] = trueReco[0][j]/all[0][j];
      printf("%f\n", all[0][j]);
      //errEfficiency[0][j] = efficiency[i][j]* TMath::Sqrt((errTrueReco[i][j]/trueReco[i][j]) * (errTrueReco[i][j]/trueReco[i][j])+(errAll[i][j]/all[i][j])*(errAll[i][j]/all[i][j]));
      hEffVsPt[0]->SetBinContent(j, efficiency[0][j]);
      //hEffVsPt[0]->SetBinError(j, errEfficiency[0][j]);

    }
    canvas[0] = new TCanvas(Form("c%d", 0), Form("%s efficiency", particles[0]), 1600, 800);
    canvas[0]->cd();
    hEffVsPt[0]->Draw();
    fout->cd();
    hEffVsPt[0]->Write(Form("%s efficiency", particles[0]));

  }//for on i as index
  //fout->Close();
return;

}
