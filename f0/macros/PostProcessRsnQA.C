Int_t PostProcessRsnQA(TString fname = "AnalysisResults.root",
		       TString production = "LHC16d6a",
		       Color_t color = kBlue,
		       Bool_t save2file = 0)
{
  
  TFile * fin = TFile::Open(fname.Data());
  if (!fin) return 1;

  TList * lrsn = (TList*) fin->Get("RsnQA_phi");
  if (!lrsn) return 2;

  TH2D * hGenEtaPt = (TH2D*) lrsn->FindObject("taskRsnQA_phi_Gen_EtaPt");
  if (!hGenEtaPt) return 3;
  
  TH1D * hGenPt = hGenEtaPt->ProjectionY("hGenPt");  //bins 6 to 15 to select eta < 0.5
  if (!hGenPt) return 4;
      
  TH2D * hRecEtaPt = (TH2D*) lrsn->FindObject("taskRsnQA_phi_Trues_EtaPt");
  if (!hRecEtaPt) return 5;
  
  TH1D * hRecPt = hRecEtaPt->ProjectionY("hRecPt"); //bins 6 to 15 to select eta < 0.5
  if (!hRecEtaPt) return 6;

   TH2D * hGenMPt = (TH2D*) lrsn->FindObject("taskRsnQA_phi_Gen_MPt");
  if (!hGenMPt) return 7;

  TH2D * hRecMPt = (TH2D*) lrsn->FindObject("taskRsnQA_phi_Trues_MPt");
  if (!hRecMPt) return 8;

  /*compute efficiencies */
  TH2D * recEffEtaPt = hRecEtaPt->Clone("EffEtaPt");
  recEffEtaPt->Divide(hGenEtaPt);
  recEffEtaPt->SetMinimum(0.);
  recEffEtaPt->SetMaximum(1.);
  recEffEtaPt->SetTitle(production.Data());
  recEffEtaPt->GetXaxis()->SetTitle("#eta");
  recEffEtaPt->GetYaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
   
  TH2D * recEffMPt = hRecMPt->Clone("EffMPt");
  recEffMPt->Divide(hGenMPt);
  recEffMPt->SetMinimum(0.);
  recEffMPt->SetMaximum(1.);
  recEffMPt->SetTitle(production.Data());
  recEffMPt->GetXaxis()->SetTitle("M (GeV/#it{c}^{2})");
  recEffMPt->GetYaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  
  TH1D * recEffPt = hRecPt->Clone("EffPt");
  recEffPt->Divide(hGenPt);
  recEffPt->GetYaxis()->SetRangeUser(0., 1.);
  recEffPt->SetTitle(production.Data());
  recEffPt->SetMarkerStyle(20);
  recEffPt->SetMarkerSize(1.);
  recEffPt->SetMarkerColor(color);
  recEffPt->SetLineColor(color);  
  recEffPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  recEffPt->GetYaxis()->SetTitle("Acceptance #times efficiency");
      
  TCanvas * cEffPt = new TCanvas("EffPt", "EffPt", 1200, 500);
  cEffPt->Divide(3,1);
  cEffPt->cd(1); recEffMPt->Draw("colz");
  cEffPt->cd(2); recEffEtaPt->Draw("colz");
  cEffPt->cd(3); gPad->SetLogx(); recEffPt->Draw();

  cEffPt->Print(Form("RsnQA_Phi_Eff_%s.pdf", production.Data()));

  if (save2file) {
    TFile * fout = new TFile(Form("RsnQA_Phi_Eff_%s.root", production.Data()), "recreate");
    fout->cd();
    recEffMPt->Write();
    recEffEtaPt->Write();
    recEffPt->Write();
    fout->Close();
  }
	       
  
  return 0;
}
