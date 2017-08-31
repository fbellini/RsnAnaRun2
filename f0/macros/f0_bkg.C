/////////////////////////////////////////////
//       bg subtraction - 29.08.2017       //
/////////////////////////////////////////////


void f0_bkg (){

  /* getting histos from root file */
  TFile * file = TFile::Open("RsnTask_f0test.root");
  TList * list = (TList*)file->Get("RsnOut_f0");
  //list->ls();


  /* like-sign pairs combination */
  TH3F * hLikePP = (TH3F*)list->FindObject("RsnTaskF0_f0_LikePP");
  TH3F * hLikeMM = (TH3F*)list->FindObject("RsnTaskF0_f0_LikeMM");
  TH3F * hUnlikePM = (TH3F*)list->FindObject("RsnTaskF0_f0_UnlikePM");
  hLikeMM->Add(hLikePP);
  TH1 * hProjXYLikePPMM = hLikeMM->Project3D("xy");
  hUnlikePM->Add(hLikeMM, -1);
  TH1 * hProjXYUnlikePM = hUnlikePM->Project3D("xy");


  /* event mixing */
  TH3F * hMixingPP = (TH3F*)list->FindObject("RsnTaskF0_f0_MixingPP");
  TH3F * hMixingMM = (TH3F*)list->FindObject("RsnTaskF0_f0_MixingMM");
  TH3F * hMixingPM = (TH3F*)list->FindObject("RsnTaskF0_f0_MixingPM");
  hMixingMM->Add(hMixingPP);
  TH1 * hProjXYMixingPPMM = hMixingMM->Project3D("xy");
  hMixingPM->Add(hMixingMM, -1);
  TH1 * hProjXYMixingPM = hMixingPM->Project3D("xy");


  /* write output */
  TCanvas * c1 = new TCanvas("c1","c1", 1200, 800);
  c1->Divide(3,2);
  c1->cd(1); hProjXYLikePPMM->Draw("colz");
  //c1->cd(2); hLikePP->Draw("colz");
  c1->cd(3); hProjXYUnlikePM->Draw("colz");
  c1->cd(4); hProjXYMixingPPMM->Draw("colz");
  //c1->cd(5); hMixingPP->Draw("colz");
  c1->cd(6); hProjXYMixingPM->Draw("colz");
  c1->Print("bgSubtraction.pdf");

  TFile *fout = TFile::Open("bgSubtraction.root", "RECREATE");
  fout->cd();
  hProjXYLikePPMM->Write();
  hProjXYUnlikePM->Write();
  hProjXYMixingPPMM->Write();
  hProjXYMixingPM->Write();
  fout->Close();

}
