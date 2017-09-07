/////////////////////////////////////////////
//       bg subtraction - 29.08.2017       //
/////////////////////////////////////////////

#include "TLegend.h"

void f0_bkg (){
  
  gROOT->LoadMacro("HistoMakeUp.C");
  
  /* getting histos from root file */
  
  TFile * file = TFile::Open("RsnTask_f0.root");
  TList * list = (TList*)file->Get("RsnOut_f0");
  //list->ls();
  TH3F * hOrLikePP = (TH3F*)list->FindObject("RsnTaskF0_f0_LikePP"); //like sign pairs: (pi+)+(pi+)
  TH3F * hOrLikeMM = (TH3F*)list->FindObject("RsnTaskF0_f0_LikeMM"); //like sign pairs: (pi-)+(pi-)
  TH3F * hOrUnlikePM = (TH3F*)list->FindObject("RsnTaskF0_f0_UnlikePM"); //same event unlike sign pairs 
  TH3F * hOrMixingPM = (TH3F*)list->FindObject("RsnTaskF0_f0_MixingPM"); //bg: (pi+)+(pi-) from 5 similar events


 /* projecting histos in 2d (losing information on multiplicity) */

  TH2D * hLikeMM = hOrLikeMM->Project3D("xy");
  TH2D * hLikePP = hOrLikePP->Project3D("xy");
  TH2D * hMEB = hOrMixingPM->Project3D("xy");
  TH2D * hUSP = hOrUnlikePM->Project3D("xy");
  
  Double_t pT[9] = {0.0, 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4.};
  TH1D * hUSPpy[8]; TH1D * hUSPpy2[8]; TH1D * hLikePPpy[8]; TH1D * hLikeMMpy[8]; TH1D * hMEBpy[8]; TH1D * hUSPClone[8];
  for(i=0; i<8; i++){
  hUSPpy[i]=0x0; hUSPpy2[i]=0x0; hLikePPpy[i]=0x0; hLikeMMpy[i]=0x0; hMEBpy[i]=0x0; hUSPClone[i]=0x0;
  }
  
  TCanvas * c1 = new TCanvas("c1","LSB Subtracted", 1600, 800);
  c1->Divide(4,2);

  TCanvas * c2 = new TCanvas("c2","Compared", 1600, 800);
  c2->Divide(4,2);
  
  for(ibin=0; ibin<8; ibin++){
  
  Int_t iMinBinPt = hUSP->GetXaxis()->FindBin(pT[ibin]);
  Int_t iMaxBinPt = hUSP->GetXaxis()->FindBin(pT[ibin+1]);
  printf("min: %d - max:%d; pT min: %5.2f - pT max: %05.2f\n", iMinBinPt, iMaxBinPt, pT[ibin], pT[ibin+1]);

  hUSPpy[ibin] = (TH1D*) hUSP->ProjectionY(Form("hUSPminusLSB_pT_%2.1f-%2.1f", pT[ibin], pT[ibin+1]), iMinBinPt, iMaxBinPt);
  hUSPpy[ibin]->SetTitle(Form("USP - LSB %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));
  hUSPpy[ibin]->GetXaxis()->SetTitle("#it{M}_{#pi#pi}(GeV/#it{c^{2}})");
  hUSPpy[ibin]->GetYaxis()->SetTitle("pairs");

  hUSPpy2[ibin] = (TH1D*) hUSP->ProjectionY(Form("hUSP2minusMEB_pT_%2.1f-%2.1f", pT[ibin], pT[ibin+1]), iMinBinPt, iMaxBinPt);
  hUSPpy2[ibin]->SetTitle(Form("USP - MEB %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));
  hUSPpy2[ibin]->GetXaxis()->SetTitle("#it{M}_{#pi#pi}(GeV/#it{c^{2}})");
  hUSPpy2[ibin]->GetYaxis()->SetTitle("pairs");

  hUSPClone[ibin] = (TH1D*) hUSPpy2[ibin]->Clone(Form("hUSP_pT_%2.1f-%2.1f", pT[ibin], pT[ibin+1]));
  hUSPClone[ibin]->SetTitle(Form("USP %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));


  hLikePPpy[ibin] = (TH1D*) hLikePP->ProjectionY(Form("hLikePP_pT_%2.1f-%2.1f", pT[ibin], pT[ibin+1]), iMinBinPt, iMaxBinPt);
  hLikePPpy[ibin]->SetTitle(Form("LSB (++) - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));
  hLikePPpy[ibin]->GetXaxis()->SetTitle("#it{M}_{#pi#pi}(GeV/#it{c^{2}})");
  hLikePPpy[ibin]->GetYaxis()->SetTitle("pairs");

  hLikeMMpy[ibin] = (TH1D*) hLikeMM->ProjectionY(Form("hLikeMM_pT_%2.1f-%2.1f", pT[ibin], pT[ibin+1]), iMinBinPt, iMaxBinPt);
  hLikeMMpy[ibin]->SetTitle(Form("LSB (--) - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));
  hLikeMMpy[ibin]->GetXaxis()->SetTitle("#it{M}_{#pi#pi}(GeV/#it{c^{2}})");
  hLikeMMpy[ibin]->GetYaxis()->SetTitle("pairs");

  hMEBpy[ibin] = (TH1D*) hMEB->ProjectionY(Form("hMEB_pT_%2.1f-%2.1f", pT[ibin], pT[ibin+1]), iMinBinPt, iMaxBinPt);
  hMEBpy[ibin]->SetTitle(Form("MEB - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));
  hMEBpy[ibin]->GetXaxis()->SetTitle("#it{M}_{#pi#pi}(GeV/#it{c^{2}})");
  hMEBpy[ibin]->GetYaxis()->SetTitle("pairs");

  hLikePPpy[ibin]->Add(hLikeMMpy[ibin]);
  hUSPpy[ibin]->Add(hLikePPpy[ibin], -1.);
  hUSPpy[ibin]->RebinX(5);

  hMEBpy[ibin]->Scale(1./5.);
  hUSPpy2[ibin]->Add(hMEBpy[ibin], -1);
  hUSPpy2[ibin]->RebinX(5);
 
  HistoMakeUp(hLikePPpy[ibin], kGreen+2, 25);
  //HistoMakeUp(hLikeMMpy[ibin], kOrange, 25);
  HistoMakeUp(hMEBpy[ibin], kRed, 24);
  HistoMakeUp(hUSPClone[ibin], kBlack, 20);

  HistoMakeUp(hUSPpy[ibin], kBlue, 25);
  //HistoMakeUp(hLikeMMpy[ibin], kOrange, 25);
  HistoMakeUp(hUSPpy2[ibin], kRed, 24);
  hUSPpy[ibin]->GetYaxis()->SetRangeUser(-2e5,2e5);
  hUSPpy2[ibin]->GetYaxis()->SetRangeUser(-2e5,2e5);
  
  c1->cd(ibin+1);
  hMEBpy[ibin]->Draw("hist");
  hLikeMMpy[ibin]->Draw("hist same");
  hLikePPpy[ibin]->Draw("hist same");
  hUSPClone[ibin]->Draw("hist same");
  c1->Print("LSB.pdf");

  c2->cd(ibin+1);
  hUSPpy[ibin]->Draw("hist");
  hUSPpy2[ibin]->Draw("hist same");
  c2->Print("Compared (LSB and MEB).pdf");
    
  }
  
/*TFile *fout = TFile::Open("bgSubtraction.root", "RECREATE");
  fout->cd();
  hLSB1->Write();
  hLSB2->Write();
  hLSB3->Write();
  hLSB4->Write();
  hMEB1->Write();
  hMEB2->Write();
  hMEB3->Write();
  hMEB4->Write();
  fout->Close();*/

}
