/////////////////////////////////////////////
//       bg subtraction - 29.08.2017       //
/////////////////////////////////////////////

#include "TLegend.h"
void f0_bkg(Double_t chosenPt = 1.0, //GeV/c --> chosen pT for plotting; if chosenPt = -1 plot all pTs
            Int_t rebinVar = 5) //rebinning factor
{

  gROOT->LoadMacro("HistoMakeUp.C");
  gROOT->LoadMacro("SetStyle.C");

  SetStyle();
  TGaxis::SetMaxDigits(3);

  /* getting histos from root file */
  TFile * file = TFile::Open("RsnTask_f0.root");
  TList * list = (TList*)file->Get("RsnOut_f0");
  //list->ls();
  TH3F * hOrUnlikePM = (TH3F*)list->FindObject("RsnTaskF0_f0_UnlikePM"); //same event unlike sign pairs
  TH3F * hOrLikePP = (TH3F*)list->FindObject("RsnTaskF0_f0_LikePP");     //like sign pairs: (pi+)+(pi+)
  TH3F * hOrLikeMM = (TH3F*)list->FindObject("RsnTaskF0_f0_LikeMM");     //like sign pairs: (pi-)+(pi-)
  TH3F * hOrMixingPP = (TH3F*)list->FindObject("RsnTaskF0_f0_MixingPP"); //(pi+)+(pi+) from 5 similar events
  TH3F * hOrMixingMM = (TH3F*)list->FindObject("RsnTaskF0_f0_MixingMM"); //(pi-)+(pi-) from 5 similar events
  TH3F * hOrMixingPM = (TH3F*)list->FindObject("RsnTaskF0_f0_MixingPM"); //(pi+)+(pi-) from 5 similar events

  /* projecting histos in 2d (losing information on multiplicity) */
  TH2D * hUSP = hOrUnlikePM->Project3D("xy");
  TH2D * hLikeMM = hOrLikeMM->Project3D("xy");
  TH2D * hLikePP = hOrLikePP->Project3D("xy");
  TH2D * hMEB = hOrMixingPM->Project3D("xy");
  TH2D * hMELMM = hOrMixingMM->Project3D("xy");
  TH2D * hMELPP = hOrMixingPP->Project3D("xy");

  /* binning: 200 MeV bins up to 2. GeV/c, then 500 MeV bins up to 5 GeV/c, one big bin 5 to 7 GeV/c */
  Double_t pT[18] = {0., 0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.5, 3., 3.5, 4., 4.5, 5., 7.};

  TH1D * hLikePPpy[18]; TH1D * hLikeMMpy[18]; TH1D * hMEBpy[18]; TH1D * hUSPpy[18]; TH1D * hLSBRatio[18];
  TH1D * hGeoMean[18]; TH1D * hSum[18]; TH1D * hnMEB[18]; TH1D * hnMEB2[18]; TH1D * hMELRatio[18];
  TH1D * hUSPminusLSB1[18]; TH1D * hUSPminusLSB2[18]; TH1D * hUSPminusMEB1[18]; TH1D * hUSPminusMEB2[18];
  TH1D * hMELMMpy[18]; TH1D * hMELPPpy[18]; TH1D * hMEL[18]; TH1D * hnMEL1[18]; TH1D * hnMEL2[18];
  TH1D * hUSPminusMEL1[18]; TH1D * hUSPminusMEL2[18]; TH1D * hMELMEBRatio[18]; TH1D * hnMELLSBRatio[18];

  for(Int_t i=0; i<17; i++){
    hLikePPpy[i]=0x0; hLikeMMpy[i]=0x0; hMEBpy[i]=0x0; hUSPpy[i]=0x0; hLSBRatio[i]=0x0;
    hGeoMean[i]=0x0; hSum[i]=0x0; hnMEB[i]=0x0; hnMEB2[i]=0x0; hMELRatio[i]=0x0;
    hUSPminusLSB1[i]=0x0; hUSPminusLSB2[i]=0x0; hUSPminusMEB1[i]=0x0; hUSPminusMEB2[i]=0x0;
    hMELPPpy[i]=0x0; hMELMMpy[i]=0x0; hMEL[i]=0x0; hnMEL1[i]=0x0; hnMEL2[i]=0x0;
    hUSPminusMEL1[i]=0x0; hUSPminusMEL2[i]=0x0; hMELMEBRatio[i]=0x0; hnMELLSBRatio[i]=0x0;
  }

  TFile *fout = TFile::Open("bgSubtraction.root", "RECREATE");

  /* canvas to print histograms for each bin in pT*/
  TCanvas *canvas[17];

  /*  bg subtraction for each bin in pT*/
  for(Int_t ibin=0; ibin<17; ibin++){

  Int_t iMinBinPt = hUSP->GetXaxis()->FindBin(pT[ibin]);
  Int_t iMaxBinPt = hUSP->GetXaxis()->FindBin(pT[ibin+1]);
  //printf("min: %d - max:%d; pT min: %5.2f - pT max: %05.2f\n", iMinBinPt, iMaxBinPt, pT[ibin], pT[ibin+1]);

  /*  projecting histos in 1d to get invariant mass spectra */
  hUSPpy[ibin] = (TH1D*) hUSP->ProjectionY(Form("hUSP_pT_%2.1f-%2.1f", pT[ibin], pT[ibin+1]), iMinBinPt, iMaxBinPt);
  hUSPpy[ibin]->SetTitle(Form("USP - LSB %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));
  hUSPpy[ibin]->RebinX(rebinVar);

  hLikeMMpy[ibin] = (TH1D*) hLikeMM->ProjectionY(Form("hLikeMM_pT_%2.1f-%2.1f", pT[ibin], pT[ibin+1]), iMinBinPt, iMaxBinPt);
  hLikeMMpy[ibin]->SetTitle(Form("LSB (--) - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));
  hLikeMMpy[ibin]->RebinX(rebinVar);

  hLikePPpy[ibin] = (TH1D*) hLikePP->ProjectionY(Form("hLikePP_pT_%2.1f-%2.1f", pT[ibin], pT[ibin+1]), iMinBinPt, iMaxBinPt);
  hLikePPpy[ibin]->SetTitle(Form("LSB (++) - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));
  hLikePPpy[ibin]->RebinX(rebinVar);

  hMEBpy[ibin] = (TH1D*) hMEB->ProjectionY(Form("hMEB_pT_%2.1f-%2.1f", pT[ibin], pT[ibin+1]), iMinBinPt, iMaxBinPt);
  hMEBpy[ibin]->SetTitle(Form("MEB - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));
  hMEBpy[ibin]->RebinX(rebinVar);

  hMELMMpy[ibin] = (TH1D*) hMELMM->ProjectionY(Form("hMELMM_pT_%2.1f-%2.1f", pT[ibin], pT[ibin+1]), iMinBinPt, iMaxBinPt);
  hMELMMpy[ibin]->SetTitle(Form("MEL (--) - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));
  hMELMMpy[ibin]->RebinX(rebinVar);

  hMELPPpy[ibin] = (TH1D*) hMELPP->ProjectionY(Form("hMELPP_pT_%2.1f-%2.1f", pT[ibin], pT[ibin+1]), iMinBinPt, iMaxBinPt);
  hMELPPpy[ibin]->SetTitle(Form("MEL (++) - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));
  hMELPPpy[ibin]->RebinX(rebinVar);

  /* Bg estimation */
  hGeoMean[ibin]=geoMean(hLikePPpy[ibin], hLikeMMpy[ibin]);   //LSB calculated with geometric mean
  hGeoMean[ibin]->SetTitle(Form("Estimated bg - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));
  //hGeoMean[ibin]->GetYaxis()->SetRangeUser(0.,4e3);
  hSum[ibin]=sum(hLikePPpy[ibin], hLikeMMpy[ibin]);           //LSB calculated with arithmetic mean
  hSum[ibin]->SetTitle(Form("LSB sum - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));

  hLSBRatio[ibin]=division(hLikeMMpy[ibin], hLikePPpy[ibin]);
  hLSBRatio[ibin]->SetTitle(Form("LS (--)/LS (++) - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));

  hnMEB[ibin]=normMEB(hMEBpy[ibin]);                          // MEB (normalized with a specific value)
  hnMEB[ibin]->SetTitle(Form("nMEB (1/5) - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));

  hnMEB2[ibin]=normMEB2(hMEBpy[ibin], hUSPpy[ibin]);          // MEB normalized to the integral of the USP
  hnMEB2[ibin]->SetTitle(Form("nMEB (integral) - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));

  hMELRatio[ibin]=division(hMELMMpy[ibin], hMELPPpy[ibin]);   // MEL -- / MEL ++
  hMELRatio[ibin]->SetTitle(Form("MEL (--)/ MEL (++) - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));

  hMEL[ibin]=geoMean(hMELPPpy[ibin],hMELMMpy[ibin]);          // MEL calculated with geometric mean
  hMEL[ibin]->SetTitle(Form("MEL geoMean - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));

  hnMEL1[ibin]=normMEB(hMEL[ibin]);                           // MEL (normalized with a specific value)
  hnMEL1[ibin]->SetTitle(Form("nMEL (1/5) - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));

  hnMEL2[ibin]=normMEB2(hMEL[ibin], hUSPpy[ibin]);            // MEL normalized to the integral of the USP
  hnMEL2[ibin]->SetTitle(Form("nMEL (integral) - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));

  hMELMEBRatio[ibin]=division(hMEL[ibin], hMEBpy[ibin]);      // MEL / MEB
  hMELMEBRatio[ibin]->SetTitle(Form("MEL / MEB - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));

  hnMELLSBRatio[ibin]=division(hnMEL2[ibin], hGeoMean[ibin]); //nMEL / LSB
  hnMELLSBRatio[ibin]->SetTitle(Form("nMEL / MEB - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));

  /* Bg subtraction */
  hUSPminusLSB1[ibin]=subtraction(hUSPpy[ibin], hSum[ibin]);
  hUSPminusLSB1[ibin]->SetTitle(Form("USP-LSB - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));
  hUSPminusLSB2[ibin]=subtraction(hUSPpy[ibin], hGeoMean[ibin]);
  hUSPminusLSB2[ibin]->SetTitle(Form("USP-bg - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));
  //hUSPminusLSB2[ibin]->GetYaxis()->SetRangeUser(0.,3.5e3);
  hUSPminusMEB1[ibin]=subtraction(hUSPpy[ibin], hnMEB[ibin]);
  hUSPminusMEB1[ibin]->SetTitle(Form("USP-MEB - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));
  hUSPminusMEB2[ibin]=subtraction(hUSPpy[ibin], hnMEB2[ibin]);
  hUSPminusMEB2[ibin]->SetTitle(Form("USP-MEB - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));
  hUSPminusMEL1[ibin]=subtraction(hUSPpy[ibin], hnMEL1[ibin]);
  hUSPminusMEL1[ibin]->SetTitle(Form("USP-MEL - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));
  hUSPminusMEL2[ibin]=subtraction(hUSPpy[ibin], hnMEL2[ibin]);
  hUSPminusMEL2[ibin]->SetTitle(Form("USP-MEL - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));

  /* Histo maquillage */
  HistoMakeUp(hGeoMean[ibin], kMagenta-2, 20, "#it{M}_{#pi#pi}(GeV/#it{c^{2}})", "pairs");
  HistoMakeUp(hLSBRatio[ibin], kAzure+7, 20, "#it{M}_{#pi#pi}(GeV/#it{c^{2}})", "pairs");
  HistoMakeUp(hnMEB2[ibin], kGreen-9, 20, "#it{M}_{#pi#pi}(GeV/#it{c^{2}})", "pairs");
  HistoMakeUp(hnMEL2[ibin], kGreen+2, 20, "#it{M}_{#pi#pi}(GeV/#it{c^{2}})", "pairs");
  HistoMakeUp(hMELRatio[ibin], kMagenta, 20, "#it{M}_{#pi#pi}(GeV/#it{c^{2}})", "ratio");
  HistoMakeUp(hMELMEBRatio[ibin], kGreen+2, 20, "#it{M}_{#pi#pi}(GeV/#it{c^{2}})", "ratio");
  HistoMakeUp(hnMELLSBRatio[ibin], kRed, 20, "#it{M}_{#pi#pi}(GeV/#it{c^{2}})", "ratio");
  HistoMakeUp(hUSPminusLSB2[ibin], kRed, 20, "#it{M}_{#pi#pi}(GeV/#it{c^{2}})", "pairs");
  HistoMakeUp(hUSPminusMEB2[ibin], kBlue+1, 20, "#it{M}_{#pi#pi}(GeV/#it{c^{2}})", "pairs");
  HistoMakeUp(hUSPminusMEL2[ibin], kGreen+2, 20, "#it{M}_{#pi#pi}(GeV/#it{c^{2}})", "pairs");
  HistoMakeUp(hUSPpy[ibin], kBlack, 20, "#it{M}_{#pi#pi}(GeV/#it{c^{2}})", "pairs");

  TLegend *legend1 = new TLegend(0.75,0.75,0.9,0.9);
  legend1->AddEntry(hGeoMean[ibin],"LSB","lpf");
  legend1->AddEntry(hnMEB2[ibin],"MEB","lpf");
  legend1->AddEntry(hnMEL2[ibin],"MEL","lpf");
  legend1->AddEntry(hUSPpy[ibin],"USP","lpf");

  TLegend *legend2 = new TLegend(0.75,0.75,0.9,0.9);
  legend2->AddEntry(hUSPminusLSB2[ibin],"LSB", "lpf");
  legend2->AddEntry(hUSPminusMEB2[ibin],"MEB","lpf");
  legend2->AddEntry(hUSPminusMEL2[ibin],"MEL","lpf");

  /* drawing and printing output */
  canvas[ibin]= new TCanvas(Form("c%d", ibin), Form("canvas %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]), 1600, 800);
  canvas[ibin]->Divide(3,2);
  canvas[ibin]->cd(1); hGeoMean[ibin]->Draw("hist");
  canvas[ibin]->cd(1); hnMEB2[ibin]->Draw("hist same");
  canvas[ibin]->cd(1); hnMEL2[ibin]->Draw("hist same");
  canvas[ibin]->cd(1); hUSPpy[ibin]->Draw("hist same");
  canvas[ibin]->cd(1); legend1->Draw();
  canvas[ibin]->cd(2); hUSPminusLSB2[ibin]->Draw("hist");
  canvas[ibin]->cd(2); hUSPminusMEB2[ibin]->Draw("hist same");
  canvas[ibin]->cd(2); hUSPminusMEL2[ibin]->Draw("hist same");
  canvas[ibin]->cd(2); legend2->Draw();
  canvas[ibin]->cd(3); hLSBRatio[ibin]->Draw("hist");
  canvas[ibin]->cd(4); hMELRatio[ibin]->Draw("hist");
  canvas[ibin]->cd(5); hMELMEBRatio[ibin]->Draw("hist");
  canvas[ibin]->cd(6); hnMELLSBRatio[ibin]->Draw("hist");


  /* writing output on file.root */
  fout->cd();
  hUSPpy[ibin]->Write();
  hLikeMMpy[ibin]->Write();
  hLikePPpy[ibin]->Write();
  hMEBpy[ibin]->Write();
  hMELMMpy[ibin]->Write();
  hMELPPpy[ibin]->Write();
  hGeoMean[ibin]->Write();
  hSum[ibin]->Write();
  hLSBRatio[ibin]->Write();
  hnMEB[ibin]->Write();
  hnMEB2[ibin]->Write();
  hMELRatio[ibin]->Write();
  hMEL[ibin]->Write();
  hnMEL1[ibin]->Write();
  hnMEL2[ibin]->Write();
  hMELMEBRatio[ibin]->Write();
  hnMELLSBRatio[ibin]->Write();
  hUSPminusLSB1[ibin]->Write();
  hUSPminusLSB2[ibin]->Write();
  hUSPminusMEB1[ibin]->Write();
  hUSPminusMEB2[ibin]->Write();
  hUSPminusMEL1[ibin]->Write();
  hUSPminusMEL2[ibin]->Write();
}
//fout->Close();
return;
}


TH1D * geoMean(TH1D *h1, TH1D *h2)
{
  //returs geometrical mean of two histograms
  if(!h1 || !h2) return;
  TH1D * hGeoMean = (TH1D*) h1->Clone();
  hGeoMean->Reset("ICES");
  Double_t cont1=0.0, cont2=0.0, product=0.0, geoMean=0.0, product2=0.0;
  for(Int_t j=1; j<h1->GetNbinsX()+1; j++){
    cont1 = h1->GetBinContent(j);
    cont2 = h2->GetBinContent(j);
    product = cont1*cont2;
    geoMean = pow(product, 1./2);
    product2 = 2*geoMean;
    hGeoMean->SetBinContent(j, product2);
    //error calculation to be implemented
    }
  return hGeoMean;
}

TH1D * sum(TH1D *h1, TH1D *h2)
{
  if(!h1 || !h2) return;
  TH1D * hResult = (TH1D*) h1->Clone();
  hResult->Add(h2);
  return hResult;
}

TH1D * division(TH1D *h1, TH1D *h2)
{
  if(!h1 || !h2) return;
  // check and raise an alarm if histograms do not have the same binning --to be imlemented
  TH1D * hResult = (TH1D*) h1->Clone();
  hResult->Divide(h2);
  return hResult;
}

TH1D * normMEB(TH1D *h1){
  if(!h1) return;
  TH1D * hResult = (TH1D*) h1->Clone();
  hResult->Scale(1./5.);
  return hResult;
}

TH1D * normMEB2(TH1D *h1, TH1D *h2){
  if(!h1 || !h2) return;
  TH1D * hResult = (TH1D*) h1->Clone();
  TH1D * h3 = (TH1D*) h1->Clone();
  TH1D * h4 = (TH1D*) h2->Clone();
  Int_t iMinBin = h3->GetXaxis()->FindBin(1.15);
  Int_t iMaxBin = h3->GetXaxis()->FindBin(1.20);
  Int_t iMinBin2 = h4->GetXaxis()->FindBin(1.15);
  Int_t iMaxBin2 = h4->GetXaxis()->FindBin(1.20);
  Double_t area1 = h3->Integral(iMinBin, iMaxBin);
  Double_t area2 = h4->Integral(iMinBin2, iMaxBin2);
  Double_t yScale = area2/area1;
  hResult->Scale(yScale);
  return hResult;
}

TH1D * subtraction(TH1D *h1, TH1D *h2){
  if(!h1) return;
  TH1D * hResult = (TH1D*) h1->Clone();
  hResult->Add(h2, -1);
  return hResult;
}
