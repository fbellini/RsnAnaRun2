/////////////////////////////////////////////
//       bg subtraction - 29.08.2017       //
/////////////////////////////////////////////


void f0_bkg (){

  /* getting histos from root file */
  TFile * file = TFile::Open("RsnTask_f0test.root");
  TList * list = (TList*)file->Get("RsnOut_f0");
  //list->ls();


  /* like-sign pairs combination */
  TH3F * hOrLikePP = (TH3F*)list->FindObject("RsnTaskF0_f0_LikePP");
  TH3F * hOrLikeMM = (TH3F*)list->FindObject("RsnTaskF0_f0_LikeMM");
  TH3F * hOrUnlikePM = (TH3F*)list->FindObject("RsnTaskF0_f0_UnlikePM");
  TH3F * hLikePP = hOrLikePP->Clone("hLikePP");
  TH3F * hLikeMM = hOrLikeMM->Clone("hLikeMM");
  TH3F * hUnlikePM = hOrUnlikePM->Clone("hUnlikePM");
  hLikeMM->Add(hLikePP);
  TH3F * hSumLikeMMPP = hLikeMM->Clone("hSumLikeMMPP");
  hUnlikePM->Add(hSumLikeMMPP, -1);
  TH2D * hProjXYUnlikePM = hUnlikePM->Project3D("xy");
  TH1 * hLastProjYLikePM = hProjXYUnlikePM->ProjectionY();
  TH1 * hLastProjYRebLikePM  = hLastProjYLikePM->RebinX(5,"hLastProjYRebLikePM");
  hLastProjYRebLikePM->SetTitle("LSB Subtracted");
  hLastProjYRebLikePM->GetXaxis()->SetTitle("M_{inv}(GeV/#it{c^{2}})");
  hLastProjYRebLikePM->GetYaxis()->SetTitle("pairs");
  hLastProjYRebLikePM->SetMarkerStyle(kFullCircle);


  /* event mixing */
  TH3F * hOrMixingPP = (TH3F*)list->FindObject("RsnTaskF0_f0_MixingPP");
  TH3F * hOrMixingMM = (TH3F*)list->FindObject("RsnTaskF0_f0_MixingMM");
  TH3F * hOrMixingPM = (TH3F*)list->FindObject("RsnTaskF0_f0_MixingPM");
  TH3F * hMixingPP = hOrMixingPP->Clone("hMixingPP");
  TH3F * hMixingMM = hOrMixingMM->Clone("hMixingMM");
  TH3F * hMixingPM = hOrMixingPM->Clone("hMixingPM");
  hMixingMM->Add(hMixingPP);
  TH3F * hSumMixingMMPP = hMixingMM->Clone("hSumMixingMMPP");
  hMixingPM->Add(hSumMixingMMPP, -1);
  TH2D * hProjXYMixingPM = hMixingPM->Project3D("xy");
  TH1 * hLastProjYMixingPM = hProjXYMixingPM->ProjectionY();
  TH1 * hLastProjYRebMixingPM  = hLastProjYMixingPM->RebinX(5,"hLastProjYRebMixingPM");
  hLastProjYRebMixingPM->SetTitle("MEB Subtracted");
  hLastProjYRebMixingPM->GetXaxis()->SetTitle("M_{inv}(GeV/#it{c^{2}})");
  hLastProjYRebMixingPM->GetXaxis()->SetRangeUser(0.6,1.2);
  hLastProjYRebMixingPM->GetYaxis()->SetTitle("pairs");
  hLastProjYRebMixingPM->SetMarkerStyle(kOpenCircle);


  /* write output */
  TCanvas * c1 = new TCanvas("c1","LSB Subtracted", 1200, 800);
  hLastProjYRebLikePM->Draw();
  TCanvas * c2 = new TCanvas("c2","MEB Subtracted", 1200, 800);
  hLastProjYRebMixingPM->Draw();
  TCanvas * c3 = new TCanvas("c3","Compared (LSB and MEB)", 1200, 800);
  hLastProjYRebLikePM->SetLineColor(kRed);
  hLastProjYRebLikePM->SetLineWidth(2);
  hLastProjYRebLikePM->SetMarkerColor(kRed);
  hLastProjYRebLikePM->SetFillColor(kRed);
  hLastProjYRebLikePM->Draw();
  hLastProjYRebMixingPM->SetLineColor(kBlue);
  hLastProjYRebMixingPM->SetLineWidth(2);
  hLastProjYRebMixingPM->SetMarkerColor(kBlue);
  hLastProjYRebMixingPM->SetFillColor(kBlue);
  hLastProjYRebMixingPM->Draw("SAME");
  gPad->BuildLegend();
  
  c1->Print("LSB.pdf");
  c2->Print("MEB.pdf");
  c3->Print("Compared (LSB and MEB).pdf");

  TFile *fout = TFile::Open("bgSubtraction.root", "RECREATE");
  fout->cd();
  hLastProjYRebLikePM->Write();
  hLastProjYRebMixingPM->Write();
  fout->Close();

}
