/////////////////////////////////////////////
//       bg subtraction - 29.08.2017       //
/////////////////////////////////////////////


void f0_bkg (){
  
  Int_t nBinRange;

  /* getting histos from root file */
  TFile * file = TFile::Open("RsnTask_f0.root");
  TList * list = (TList*)file->Get("RsnOut_f0");
  //list->ls();
  TH3F * hOrLikePP = (TH3F*)list->FindObject("RsnTaskF0_f0_LikePP"); //like sign pairs: (pi+)+(pi+)
  TH3F * hOrLikeMM = (TH3F*)list->FindObject("RsnTaskF0_f0_LikeMM"); //like sign pairs: (pi-)+(pi-)
  TH3F * hOrUnlikePM = (TH3F*)list->FindObject("RsnTaskF0_f0_UnlikePM"); //same event unlike sign pairs 
  TH3F * hOrMixingPP = (TH3F*)list->FindObject("RsnTaskF0_f0_MixingPP"); //bg: (pi+)+(pi+) from 5 similar events
  TH3F * hOrMixingMM = (TH3F*)list->FindObject("RsnTaskF0_f0_MixingMM"); //bg: (pi-)+(pi-) from 5 similar events
  TH3F * hOrMixingPM = (TH3F*)list->FindObject("RsnTaskF0_f0_MixingPM"); //bg: (pi+)+(pi-) from 5 similar events


  /* like-sign pairs combination */
  
  TH3F * hLikePP = hOrLikePP->Clone("hLikePP");
  TH3F * hLikeMM = hOrLikeMM->Clone("hLikeMM");
  TH3F * hUnlikePM = hOrUnlikePM->Clone("hUnlikePM");
  
  hLikeMM->Add(hLikePP);
  TH3F * hSumLikeMMPP = hLikeMM->Clone("hSumLikeMMPP");
  hSumLikeMMPP->Scale(1/2);
  hUnlikePM->Add(hSumLikeMMPP, -1);
  
  TH2D * hProjXYUnlikePM = hUnlikePM->Project3D("xy");
  TH1  * hProjYUnlikePM = hProjXYUnlikePM->ProjectionY();
  nBinRange=hProjYUnlikePM->GetXaxis()->FindBin(1.2) - hProjYUnlikePM->GetXaxis()->FindBin(0.6);
  printf("Range 0.6 - 1.2 GeV -> total number of bins = %d\n", nBinRange);
  TH1 * hLSB  = hProjYUnlikePM->RebinX(5,"hLSB");
  hLSB->SetTitle("LSB Subtracted");
  hLSB->GetXaxis()->SetTitle("M_{inv}(GeV/#it{c^{2}})");
  hLSB->GetYaxis()->SetTitle("pairs");
  hLSB->SetMarkerStyle(kFullCircle);


  /* event mixing */
  
  TH3F * hMixingPM = hOrMixingPM->Clone("hMixingPM");
  TH3F * hUnlikePM2 = hOrUnlikePM->Clone("hUnlikePM2");
  
  hMixingPM->Scale(1/5);
  hUnlikePM2->Add(hMixingPM, -1);
  TH2D * hProjXYMEB = hUnlikePM2->Project3D("xy");
  TH1 * hProjYMEB = hProjXYMEB->ProjectionY();
  TH1 * hMEB  = hProjYMEB->RebinX(5,"hMEB");
  
  hMEB->SetTitle("MEB Subtracted");
  hMEB->GetXaxis()->SetTitle("M_{inv}(GeV/#it{c^{2}})");
  hMEB->GetXaxis()->SetRangeUser(0.6,1.2);
  hMEB->GetYaxis()->SetTitle("pairs");
  hMEB->SetMarkerStyle(kOpenCircle);


  /* write output */
  TCanvas * c1 = new TCanvas("c1","LSB Subtracted", 1200, 800);
  hLSB->Draw();
  TCanvas * c2 = new TCanvas("c2","MEB Subtracted", 1200, 800);
  hMEB->Draw();
  TCanvas * c3 = new TCanvas("c3","Compared (LSB and MEB)", 1200, 800);
  hLSB->SetLineColor(kRed);
  hLSB->SetLineWidth(2);
  hLSB->SetMarkerColor(kRed);
  hLSB->SetFillColor(kRed);
  hLSB->Draw();
  hMEB->SetLineColor(kBlue);
  hMEB->SetLineWidth(2);
  hMEB->SetMarkerColor(kBlue);
  hMEB->SetFillColor(kBlue);
  hMEB->Draw("SAME");
  hOrUnlikePM->SetLineColor(kBlack);
  hOrUnlikePM->SetLineWidth(2);
  hOrUnlikePM->SetMarkerColor(kBlack);
  hOrUnlikePM->SetFillColor(kBlack);
  hOrUnlikePM->Draw("SAME");
  gPad->BuildLegend();
  
  c1->Print("LSB.pdf");
  c2->Print("MEB.pdf");
  c3->Print("Compared (LSB and MEB).pdf");

  TFile *fout = TFile::Open("bgSubtraction.root", "RECREATE");
  fout->cd();
  hLSB->Write();
  hMEB->Write();
  fout->Close();

}
