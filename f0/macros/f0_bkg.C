/////////////////////////////////////////////
//       bg subtraction - 29.08.2017       //
/////////////////////////////////////////////

#include "TLegend.h"

void f0_bkg (){
  
  /* getting histos from root file */
  
  TFile * file = TFile::Open("RsnTask_f0.root");
  TList * list = (TList*)file->Get("RsnOut_f0");
  //list->ls();
  TH3F * hOrLikePP = (TH3F*)list->FindObject("RsnTaskF0_f0_LikePP"); //like sign pairs: (pi+)+(pi+)
  TH3F * hOrLikeMM = (TH3F*)list->FindObject("RsnTaskF0_f0_LikeMM"); //like sign pairs: (pi-)+(pi-)
  TH3F * hOrUnlikePM = (TH3F*)list->FindObject("RsnTaskF0_f0_UnlikePM"); //same event unlike sign pairs 
  TH3F * hOrMixingPM = (TH3F*)list->FindObject("RsnTaskF0_f0_MixingPM"); //bg: (pi+)+(pi-) from 5 similar events


  /* cloning useful histos */

  TH3F * hSumLikePPMM = hOrLikeMM->Clone("hSumLikePPMM");
  TH3F * hUnlikePM1 = hOrUnlikePM->Clone("hUnlikePM1");
  TH3F * hUnlikePM2 = hOrUnlikePM->Clone("hUnlikePM2");
  
  
  /* like-sign pairs combination */
  
  hSumLikePPMM->Add(hOrLikePP);
  hSumLikePPMM->Scale(1./2.);
  hUnlikePM1->Add(hSumLikePPMM, -1);
  
  TH2D * hProjXYUnlikePM = hUnlikePM1->Project3D("xy");
  //Int_t minBinPt = hProjXYUnlikePM->GetXaxis()->FindBin(2);
  //Int_t maxBinPt = hProjXYUnlikePM->GetXaxis()->FindBin(3.5);
  //printf("min: %d - max:%d\n", minBinPt, maxBinPt);
  TH1  * hProjYUnlikePM1 = hProjXYUnlikePM->ProjectionY("_LSB1y", 6, 11);  //0.5 GeV < pT < 1 GeV
  TH1  * hProjYUnlikePM2 = hProjXYUnlikePM->ProjectionY("_LSB2y", 15, 20); //1.5 GeV < pT < 2 GeV
  TH1  * hProjYUnlikePM3 = hProjXYUnlikePM->ProjectionY("_LSB3y", 20, 35); //2 GeV < pT < 3.5 GeV
  TH1  * hProjYUnlikePM4 = hProjXYUnlikePM->ProjectionY("_LSB4y", 35, 40); //3.5 GeV < pT < 4 GeV
  // Int_t nBinRange=hProjYUnlikePM->GetXaxis()->FindBin(1.2) - hProjYUnlikePM->GetXaxis()->FindBin(0.6);
  // printf("Range 0.6 - 1.2 GeV -> total number of bins = %d\n", nBinRange);
  TH1 * hLSB1  = hProjYUnlikePM1->RebinX(5,"hLSB1");
  TH1 * hLSB2  = hProjYUnlikePM2->RebinX(5,"hLSB2");
  TH1 * hLSB3  = hProjYUnlikePM3->RebinX(5,"hLSB3");
  TH1 * hLSB4  = hProjYUnlikePM4->RebinX(5,"hLSB4");
  
  TH1 * hLSB1b = hLSB1->Clone("hLSB1b");
  hLSB1b->SetTitle("Estimated background - 0.5 < p_{T} < 1 (GeV/#it{c})");
  hLSB1b->GetXaxis()->SetTitle("M_{inv}(GeV/#it{c^{2}})");
  hLSB1b->GetYaxis()->SetTitle("pairs");
  TH1 * hLSB2b = hLSB2->Clone("hLSB2b");
  hLSB2b->SetTitle("Estimated background - 1.5 < p_{T} < 2 (GeV/#it{c})");
  hLSB2b->GetXaxis()->SetTitle("M_{inv}(GeV/#it{c^{2}})");
  hLSB2b->GetYaxis()->SetTitle("pairs");
  TH1 * hLSB3b = hLSB3->Clone("hLSB3b");
  hLSB3b->SetTitle("Estimated background - 2 < p_{T} < 3.5 (GeV/#it{c})");
  hLSB3b->GetXaxis()->SetTitle("M_{inv}(GeV/#it{c^{2}})");
  hLSB3b->GetYaxis()->SetTitle("pairs");
  TH1 * hLSB4b = hLSB4->Clone("hLSB4b");
  hLSB4b->SetTitle("Estimated background - 3.5 < p_{T} < 4 (GeV/#it{c})");
  hLSB4b->GetXaxis()->SetTitle("M_{inv}(GeV/#it{c^{2}})");
  hLSB4b->GetYaxis()->SetTitle("pairs");
  
  
  /* event mixing */
  
  hOrMixingPM->Scale(1./5.);
  hUnlikePM2->Add(hOrMixingPM, -1);
  TH2D * hProjXYMEB = hUnlikePM2->Project3D("xy");
  TH1 * hProjYMEB1 = hProjXYMEB->ProjectionY("_MEB1y", 6, 11); //0.5 GeV < pT < 1 GeV
  TH1 * hProjYMEB2 = hProjXYMEB->ProjectionY("_MEB2y", 15, 20);//1.5 GeV < pT < 2 GeV
  TH1 * hProjYMEB3 = hProjXYMEB->ProjectionY("_MEB3y", 20, 35);//2 GeV < pT < 3.5 GeV
  TH1 * hProjYMEB4 = hProjXYMEB->ProjectionY("_MEB4y", 35, 40);//3.5 GeV < pT < 4 GeV
  TH1 * hMEB1  = hProjYMEB1->RebinX(5,"hMEB1");
  TH1 * hMEB2  = hProjYMEB2->RebinX(5,"hMEB2");
  TH1 * hMEB3  = hProjYMEB3->RebinX(5,"hMEB3");
  TH1 * hMEB4  = hProjYMEB4->RebinX(5,"hMEB4");
  
  TH1 * hMEB1b = hMEB1->Clone("hMEB1b");
  TH1 * hMEB2b = hMEB2->Clone("hLSB2b");
  TH1 * hMEB3b = hMEB3->Clone("hLSB3b");
  TH1 * hMEB4b = hMEB4->Clone("hLSB4b");
  
  
  /* write output */
  
  TCanvas * c1 = new TCanvas("c1","LSB Subtracted", 800, 800);
  c1->Divide(2,2);
  c1->cd(1);
  hLSB1->SetTitle("LSB Subtracted - 0.5 < p_{T} < 1 (GeV/#it{c})");
  hLSB1->GetXaxis()->SetTitle("M_{inv}(GeV/#it{c^{2}})");
  hLSB1->GetYaxis()->SetTitle("pairs");
  hLSB1->Draw();
  c1->cd(2);
  hLSB2->SetTitle("LSB Subtracted - 1.5 < p_{T} < 2 (GeV/#it{c})");
  hLSB2->GetXaxis()->SetTitle("M_{inv}(GeV/#it{c^{2}})");
  hLSB2->GetYaxis()->SetTitle("pairs");
  hLSB2->Draw();
  c1->cd(3);
  hLSB3->SetTitle("LSB Subtracted - 2 < p_{T} < 3.5 (GeV/#it{c})");
  hLSB3->GetXaxis()->SetTitle("M_{inv}(GeV/#it{c^{2}})");
  hLSB3->GetYaxis()->SetTitle("pairs");
  hLSB3->Draw();
  c1->cd(4);
  hLSB4->SetTitle("LSB Subtracted - 3.5 < p_{T} < 4 (GeV/#it{c})");
  hLSB4->GetXaxis()->SetTitle("M_{inv}(GeV/#it{c^{2}})");
  hLSB4->GetYaxis()->SetTitle("pairs");
  hLSB4->Draw();

  TCanvas * c2 = new TCanvas("c2","MEB Subtracted", 800, 800);
  c2->Divide(2,2);
  c2->cd(1);
  hMEB1->SetTitle("MEB Subtracted - 0.5 < p_{T} < 1 (GeV/#it{c})");
  hMEB1->GetXaxis()->SetTitle("M_{inv}(GeV/#it{c^{2}})");
  hMEB1->GetYaxis()->SetTitle("pairs");
  hMEB1->Draw();
  c2->cd(2);
  hMEB2->SetTitle("MEB Subtracted - 1.5 < p_{T} < 2 (GeV/#it{c})");
  hMEB2->GetXaxis()->SetTitle("M_{inv}(GeV/#it{c^{2}})");
  hMEB2->GetYaxis()->SetTitle("pairs");
  hMEB2->Draw();
  c2->cd(3);
  hMEB3->SetTitle("MEB Subtracted - 2 < p_{T} < 3.5 (GeV/#it{c})");
  hMEB3->GetXaxis()->SetTitle("M_{inv}(GeV/#it{c^{2}})");
  hMEB3->GetYaxis()->SetTitle("pairs");
  hMEB3->Draw();
  c2->cd(4);
  hMEB4->SetTitle("MEB Subtracted - 3.5 < p_{T} < 4 (GeV/#it{c})");
  hMEB4->GetXaxis()->SetTitle("M_{inv}(GeV/#it{c^{2}})");
  hMEB4->GetYaxis()->SetTitle("pairs");
  hMEB4->Draw();

  TCanvas * c3 = new TCanvas("c3","Compared (LSB and MEB)", 800, 800);
  c3->Divide(2,2);
  c3->cd(1);
  gStyle->SetOptStat(0);
  hLSB1b->SetLineColor(kRed);
  hLSB1b->SetLineWidth(1);
  hLSB1b->SetMarkerColor(kRed);
  hLSB1b->SetFillColor(kRed);
  hLSB1b->Draw();
  hMEB1b->SetLineColor(kBlue);
  hMEB1b->SetLineWidth(1);
  hMEB1b->SetMarkerColor(kBlue);
  hMEB1b->SetFillColor(kBlue);
  hMEB1b->Draw("SAME");
  TLegend *legend = new TLegend(0.75,0.75,0.9,0.9);
  legend->AddEntry(hLSB1b,"LS bkg", "lpf");
  legend->AddEntry(hMEB1b,"ME bkg","lpf");
  legend->Draw();
  c3->cd(2);
  gStyle->SetOptStat(0);
  hLSB2b->SetLineColor(kRed);
  hLSB2b->SetLineWidth(1);
  hLSB2b->SetMarkerColor(kRed);
  hLSB2b->SetFillColor(kRed);
  hLSB2b->Draw();
  hMEB2b->SetLineColor(kBlue);
  hMEB2b->SetLineWidth(1);
  hMEB2b->SetMarkerColor(kBlue);
  hMEB2b->SetFillColor(kBlue);
  hMEB2b->Draw("SAME");
  legend->Draw();
  c3->cd(3);
  gStyle->SetOptStat(0);
  hLSB3b->SetLineColor(kRed);
  hLSB3b->SetLineWidth(1);
  hLSB3b->SetMarkerColor(kRed);
  hLSB3b->SetFillColor(kRed);
  hLSB3b->Draw();
  hMEB3b->SetLineColor(kBlue);
  hMEB3b->SetLineWidth(1);
  hMEB3b->SetMarkerColor(kBlue);
  hMEB3b->SetFillColor(kBlue);
  hMEB3b->Draw("SAME");
  legend->Draw();
  c3->cd(4);
  gStyle->SetOptStat(0);
  hLSB4b->SetLineColor(kRed);
  hLSB4b->SetLineWidth(1);
  hLSB4b->SetMarkerColor(kRed);
  hLSB4b->SetFillColor(kRed);
  hLSB4b->Draw();
  hMEB4b->SetLineColor(kBlue);
  hMEB4b->SetLineWidth(1);
  hMEB4b->SetMarkerColor(kBlue);
  hMEB4b->SetFillColor(kBlue);
  hMEB4b->Draw("SAME");
  legend->Draw();
  
  c1->Print("LSB.pdf");
  c2->Print("MEB.pdf");
  c3->Print("Compared (LSB and MEB).pdf");

  TFile *fout = TFile::Open("bgSubtraction.root", "RECREATE");
  fout->cd();
  hLSB1->Write();
  hLSB2->Write();
  hLSB3->Write();
  hLSB4->Write();
  hMEB1->Write();
  hMEB2->Write();
  hMEB3->Write();
  hMEB4->Write();
  fout->Close();

}
