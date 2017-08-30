/////////////////////////////////////////////
//       bg subtraction - 29.08.2017       //
/////////////////////////////////////////////


void f0_bkg (){

  /* getting histos from root file */
  TFile * file = TFile::Open("RsnTask_f0test.root");
  TList * list = (TList*)file->Get("RsnOut_f0");
  list->ls();


  /* like-sign pairs combination */
  TH3F * hLikePP = (TH3F*)list->FindObject("RsnTaskF0_f0_LikePP");
  TH3F * hLikeMM = (TH3F*)list->FindObject("RsnTaskF0_f0_LikeMM");
  TH3F * hUnlikePM = (TH3F*)list->FindObject("RsnTaskF0_f0_UnlikePM");
  hLikeMM->Add(hLikePP);
  hUnlikePM->Add(hLikeMM, -1);

  /* event mixing */
  TH3F * hMixingPP = (TH3F*)list->FindObject("RsnTaskF0_f0_MixingPP");
  TH3F * hMixingMM = (TH3F*)list->FindObject("RsnTaskF0_f0_MixingMM");
  TH3F * hMixingPM = (TH3F*)list->FindObject("RsnTaskF0_f0_MixingPM");
  hMixingMM->Add(hMixingPP);
  hMixingPM->Add(hMixingMM, -1);

  /* write output */
  TCanvas * c1 = new TCanvas("c1","c1", 1200, 800);
  c1->Divide(3,2);
  c1->cd(1); hLikeMM->Draw("colz");
  c1->cd(2); hLikePP->Draw("colz");
  c1->cd(3); hUnlikePM->Draw("colz");
  c1->cd(4); hMixingMM->Draw("colz");
  c1->cd(5); hMixingPP->Draw("colz");
  c1->cd(6); hMixingPM->Draw("colz");
  c1->Print("bgSubtraction.pdf");

  TFile *fout = TFile::Open("bgSubtraction.root", "RECREATE");
  fout->cd();
  hLikeMM->Write();
  hLikePP->Write();
  hUnlikePM->Write();
  hMixingMM->Write();
  hMixingPP->Write();
  hMixingPM->Write();
  fout->Close();

}
