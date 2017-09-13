/////////////////////////////////////////////
//       bg subtraction - 29.08.2017       //
/////////////////////////////////////////////

#include "TLegend.h"

TH1D * geoMean(TH1D *h1, TH1D *h2){
  if(!h1) return;
  TH1D * hGeoMean = (TH1D*) h1->Clone();
  Double_t cont1=0, cont2=0, product=0, geoMean=0, product2=0;
  hGeoMean->Reset("ICES");
  Double_t nBinRange=h1->GetXaxis()->FindBin(1.2) - h1->GetXaxis()->FindBin(0.6); // 0.6 < M_inv < 1.2 GeV/c2
  for(Int_t j=0; j<nBinRange; j++){
    cont1 = h1->GetBinContent(j);
    cont2 = h2->GetBinContent(j);
    product = cont1*cont2;
    geoMean = pow(product, 1./2);
    product2 = 2*geoMean;
    hGeoMean->SetBinContent(j, product2);
    }
  return hGeoMean;
}

TH1D * sum(TH1D *h1, TH1D *h2){
  if(!h1) return;
  TH1D * hResult = (TH1D*) h1->Clone();
  hResult->Add(h2);
  return hResult;
}

TH1D * division(TH1D *h1, TH1D *h2){
  if(!h1) return;
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
  if(!h1) return;
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

void f0_bkg (){
  
  gROOT->LoadMacro("HistoMakeUp.C");
  gROOT->LoadMacro("SetStyle.C");

  SetStyle();
  TGaxis::SetMaxDigits(3);

  Int_t rebinVar=5; //rebinning factor
  
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
  
  /* binning: 200 MeV bins up to 2. GeV/c, then 500 MeV bins up to 5 GeV/c, one big bin 5 to 7 GeV/c */
  Double_t pT[18] = {0., 0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.5, 3., 3.5, 4., 4.5, 5., 7.};
  
  TH1D * hLikePPpy[18]; TH1D * hLikeMMpy[18]; TH1D * hMEBpy[18]; TH1D * hUSPpy[18]; TH1D * hLSBRatio[18];
  TH1D * hGeoMean[18]; TH1D * hSum[18]; TH1D * hNormMEB[18]; TH1D * hNormMEB2[18];
  TH1D * hUSPminusLSB1[18]; TH1D * hUSPminusLSB2[18]; TH1D * hUSPminusMEB1[18]; TH1D * hUSPminusMEB2[18];

  for(Int_t i=0; i<17; i++){
  hLikePPpy[i]=0x0; hLikeMMpy[i]=0x0; hMEBpy[i]=0x0; hUSPpy[i]=0x0; hLSBRatio[i]=0x0;
  hGeoMean[i]=0x0; hSum[i]=0x0; hNormMEB[i]=0x0; hNormMEB2[i]=0x0;
  hUSPminusLSB1[i]=0x0; hUSPminusLSB2[i]=0x0; hUSPminusMEB1[i]=0x0; hUSPminusMEB2[i]=0x0;
  }
  
  TFile *fout = TFile::Open("bgSubtraction.root", "RECREATE");

  TCanvas * c1 = new TCanvas("c1","EstimatedBg1", 1600, 800);
  c1->Divide(5,2);

  TCanvas * c2 = new TCanvas("c2","EstimatedBg2", 1600, 800);
  c2->Divide(4,2);

  TCanvas * c3 = new TCanvas("c3","USP-Bg1", 1600, 800);
  c3->Divide(5,2);

  TCanvas * c4 = new TCanvas("c4","USP-Bg2", 1600, 800);
  c4->Divide(4,2);

  TCanvas * c5 = new TCanvas("c5","LSB_Ratio1", 1600, 800);
  c5->Divide(5,2);

  TCanvas * c6 = new TCanvas("c6","LSB_Ratio2", 1600, 800);
  c6->Divide(4,2);

  TCanvas * c7 = new TCanvas("c7","LSBComparison1", 1600, 800);
  c7->Divide(5,2);

  TCanvas * c8 = new TCanvas("c8","LSBComparison2", 1600, 800);
  c8->Divide(4,2);

  TCanvas * c9 = new TCanvas("c9","MEBComparison1", 1600, 800);
  c9->Divide(5,2);

  TCanvas * c10 = new TCanvas("c10","MEBComparison2", 1600, 800);
  c10->Divide(4,2);
  
  /*  bg subtraction for each bin in pT*/
  for(Int_t ibin=0; ibin<17; ibin++){
  
  Int_t iMinBinPt = hUSP->GetXaxis()->FindBin(pT[ibin]);
  Int_t iMaxBinPt = hUSP->GetXaxis()->FindBin(pT[ibin+1]);
  //printf("min: %d - max:%d; pT min: %5.2f - pT max: %05.2f\n", iMinBinPt, iMaxBinPt, pT[ibin], pT[ibin+1]);

  /*  projecting histos in 1d to get invariant mass spectra */
  hLikePPpy[ibin] = (TH1D*) hLikePP->ProjectionY(Form("hLikePP_pT_%2.1f-%2.1f", pT[ibin], pT[ibin+1]), iMinBinPt, iMaxBinPt);
  hLikePPpy[ibin]->SetTitle(Form("LSB (++) - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));
  hLikePPpy[ibin]->RebinX(rebinVar);
    
  hLikeMMpy[ibin] = (TH1D*) hLikeMM->ProjectionY(Form("hLikeMM_pT_%2.1f-%2.1f", pT[ibin], pT[ibin+1]), iMinBinPt, iMaxBinPt);
  hLikeMMpy[ibin]->SetTitle(Form("LSB (--) - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));
  hLikeMMpy[ibin]->RebinX(rebinVar);
    
  hMEBpy[ibin] = (TH1D*) hMEB->ProjectionY(Form("hMEB_pT_%2.1f-%2.1f", pT[ibin], pT[ibin+1]), iMinBinPt, iMaxBinPt);
  hMEBpy[ibin]->SetTitle(Form("MEB - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));
  hMEBpy[ibin]->RebinX(rebinVar);
    
  hUSPpy[ibin] = (TH1D*) hUSP->ProjectionY(Form("hUSP_pT_%2.1f-%2.1f", pT[ibin], pT[ibin+1]), iMinBinPt, iMaxBinPt);
  hUSPpy[ibin]->SetTitle(Form("USP - LSB %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));
  hUSPpy[ibin]->RebinX(rebinVar);
    
  /* Bg estimation */
  hGeoMean[ibin]=geoMean(hLikePPpy[ibin], hLikeMMpy[ibin]);   //LSB calculated with geometric mean
  hGeoMean[ibin]->SetTitle(Form("LSB calculated with geometric mean - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));
  hSum[ibin]=sum(hLikePPpy[ibin], hLikeMMpy[ibin]);           //LSB calculated with arithmetic mean
  hSum[ibin]->SetTitle(Form("LSB calculated adding LS (++) and LS (--) - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));
  hLSBRatio[ibin]=division(hLikeMMpy[ibin], hLikePPpy[ibin]);
  hLSBRatio[ibin]->SetTitle(Form("LS (--)/LS (++) - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));
  hNormMEB[ibin]=normMEB(hMEBpy[ibin]);                       // MEB (normalized with a specific value)
  hNormMEB[ibin]->SetTitle(Form("MEB normalized with 1/5 - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));
  hNormMEB2[ibin]=normMEB2(hMEBpy[ibin], hUSPpy[ibin]);       // MEB normalized to the integral of the USP
  hNormMEB2[ibin]->SetTitle(Form("MEB normalized to the integral of the USP - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));

  /* Bg subtraction */
  hUSPminusLSB1[ibin]=subtraction(hUSPpy[ibin], hSum[ibin]);
  hUSPminusLSB1[ibin]->SetTitle(Form("USP-LSB - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));
  hUSPminusLSB2[ibin]=subtraction(hUSPpy[ibin], hGeoMean[ibin]);
  hUSPminusLSB2[ibin]->SetTitle(Form("USP-LSB - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));
  hUSPminusMEB1[ibin]=subtraction(hUSPpy[ibin], hNormMEB[ibin]);
  hUSPminusMEB1[ibin]->SetTitle(Form("USP-MEB - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));
  hUSPminusMEB2[ibin]=subtraction(hUSPpy[ibin], hNormMEB2[ibin]);
  hUSPminusMEB2[ibin]->SetTitle(Form("USP-MEB - %2.1f < p_{T} < %2.1f GeV/#it{c}", pT[ibin], pT[ibin+1]));

  /* Histo maquillage */
  HistoMakeUp(hGeoMean[ibin], kMagenta-2, 20, "#it{M}_{#pi#pi}(GeV/#it{c^{2}})", "pairs");
  HistoMakeUp(hSum[ibin], kMagenta-9, 20, "#it{M}_{#pi#pi}(GeV/#it{c^{2}})", "pairs");
  HistoMakeUp(hLSBRatio[ibin], kAzure+7, 20, "#it{M}_{#pi#pi}(GeV/#it{c^{2}})", "pairs");
  HistoMakeUp(hNormMEB[ibin], kGreen+2, 20, "#it{M}_{#pi#pi}(GeV/#it{c^{2}})", "pairs");
  HistoMakeUp(hNormMEB2[ibin], kGreen-9, 20, "#it{M}_{#pi#pi}(GeV/#it{c^{2}})", "pairs");
  HistoMakeUp(hUSPminusLSB1[ibin], kRed-7, 20, "#it{M}_{#pi#pi}(GeV/#it{c^{2}})", "pairs");
  HistoMakeUp(hUSPminusLSB2[ibin], kRed, 20, "#it{M}_{#pi#pi}(GeV/#it{c^{2}})", "pairs");
  HistoMakeUp(hUSPminusMEB1[ibin], kBlue-7, 20, "#it{M}_{#pi#pi}(GeV/#it{c^{2}})", "pairs");
  HistoMakeUp(hUSPminusMEB2[ibin], kBlue+1, 20, "#it{M}_{#pi#pi}(GeV/#it{c^{2}})", "pairs");
  HistoMakeUp(hUSPpy[ibin], kBlack, 20, "#it{M}_{#pi#pi}(GeV/#it{c^{2}})", "pairs");
  
  //hUSPpy[ibin]->GetYaxis()->SetRangeUser(-2e5,2e5);
  //hUSPpy2[ibin]->GetYaxis()->SetRangeUser(-2e5,2e5);
  
  /* drawing and printing output */
  if(ibin<10){
  c1->cd(ibin+1);
  hGeoMean[ibin]->Draw("hist");
  hNormMEB[ibin]->Draw("hist same");
  hUSPpy[ibin]->Draw("hist same");
  TLegend *legend1 = new TLegend(0.75,0.75,0.9,0.9);
  legend1->AddEntry(hNormMEB[ibin],"MEB", "lpf");
  legend1->AddEntry(hGeoMean[ibin],"LSB","lpf");
  legend1->AddEntry(hUSPpy[ibin],"USP","lpf");
  legend1->Draw();

  c3->cd(ibin+1);
  hUSPminusLSB1[ibin]->Draw("hist");
  hUSPminusMEB2[ibin]->Draw("hist same");
  TLegend *legend3 = new TLegend(0.75,0.75,0.9,0.9);
  legend3->AddEntry(hUSPminusLSB1[ibin],"LSB", "lpf");
  legend3->AddEntry(hUSPminusMEB2[ibin],"MEB","lpf");
  legend3->Draw();

  c5->cd(ibin+1);
  hLSBRatio[ibin]->Draw("hist");

  c7->cd(ibin+1);
  hUSPminusLSB1[ibin]->Draw("hist");
  hUSPminusLSB2[ibin]->Draw("hist same");
  TLegend *legend7 = new TLegend(0.75,0.75,0.9,0.9);
  legend7->AddEntry(hUSPminusLSB1[ibin],"LSB 1", "lpf");
  legend7->AddEntry(hUSPminusLSB2[ibin],"LSB 2","lpf");
  legend7->Draw();

  c9->cd(ibin+1);
  hUSPminusMEB1[ibin]->Draw("hist");
  hUSPminusMEB2[ibin]->Draw("hist same");
  TLegend *legend8 = new TLegend(0.75,0.75,0.9,0.9);
  legend8->AddEntry(hUSPminusMEB1[ibin],"MEB 1", "lpf");
  legend8->AddEntry(hUSPminusMEB2[ibin],"MEB 2","lpf");
  legend8->Draw();
  }
  else{
  c2->cd(ibin-9);
  hGeoMean[ibin]->Draw("hist");
  hNormMEB[ibin]->Draw("hist same");
  hUSPpy[ibin]->Draw("hist same");
  TLegend *legend2 = new TLegend(0.75,0.75,0.9,0.9);
  legend2->AddEntry(hNormMEB[ibin],"MEB", "lpf");
  legend2->AddEntry(hGeoMean[ibin],"LSB","lpf");
  legend2->AddEntry(hUSPpy[ibin],"USP","lpf");
  legend2->Draw();

  c4->cd(ibin-9);
  hUSPminusLSB1[ibin]->Draw("hist");
  hUSPminusMEB2[ibin]->Draw("hist same");
  TLegend *legend4 = new TLegend(0.75,0.75,0.9,0.9);
  legend4->AddEntry(hUSPminusLSB1[ibin],"LSB", "lpf");
  legend4->AddEntry(hUSPminusMEB2[ibin],"MEB","lpf");
  legend4->Draw();

  c6->cd(ibin-9);
  hLSBRatio[ibin]->Draw("hist");

  c8->cd(ibin-9);
  hUSPminusLSB1[ibin]->Draw("hist");
  hUSPminusLSB2[ibin]->Draw("hist same");
  TLegend *legend8 = new TLegend(0.75,0.75,0.9,0.9);
  legend8->AddEntry(hUSPminusLSB1[ibin],"LSB 1", "lpf");
  legend8->AddEntry(hUSPminusLSB2[ibin],"LSB 2","lpf");
  legend8->Draw();

  c10->cd(ibin-9);
  hUSPminusMEB1[ibin]->Draw("hist");
  hUSPminusMEB2[ibin]->Draw("hist same");
  TLegend *legend10 = new TLegend(0.75,0.75,0.9,0.9);
  legend10->AddEntry(hUSPminusMEB1[ibin],"MEB 1", "lpf");
  legend10->AddEntry(hUSPminusMEB2[ibin],"MEB 2","lpf");
  legend10->Draw();
}
  c1->Print("EstimatedBg1.pdf");
  c2->Print("EstimatedBg2.pdf");
  c3->Print("USPSubtracted1.pdf");
  c4->Print("USPSubtracted2.pdf");
  c5->Print("RatioLSB1.pdf");
  c6->Print("RatioLSB2.pdf");
  c7->Print("LSBComparison1.pdf");
  c8->Print("LSBComparison2.pdf");
  c9->Print("MEBComparison1.pdf");
  c10->Print("MEBComparison2.pdf");
    
  /* writing output on file.root */
  fout->cd();
  hLikePPpy[ibin]->Write();
  hLikeMMpy[ibin]->Write();
  hMEBpy[ibin]->Write();
  hUSPpy[ibin]->Write();
  hGeoMean[ibin]->Write();
  hSum[ibin]->Write();
  hLSBRatio[ibin]->Write();
  hNormMEB[ibin]->Write();
  hNormMEB2[ibin]->Write();
  hUSPminusLSB1[ibin]->Write();
  hUSPminusLSB2[ibin]->Write();
  hUSPminusMEB1[ibin]->Write();
  hUSPminusMEB2[ibin]->Write();  
}
//fout->Close();
return;
}
