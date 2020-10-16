void FixSpectraRange()
{
    TFile fin("FINAL_Spectra_phi_XeXe544TeV.root", "read");
    TH1D * hin[5];
    TH1D * hinsy[5];

    TFile fout("Reranged_FINAL_Spectra_phi_XeXe544TeV.root", "recreate");
    for (int j=0;j<5;j++){ 
      hin[j] = (TH1D*) fin.Get(Form("hCorrected_%i", j));   
      hinsy[j] = (TH1D*) fin.Get(Form("hCorrected_%i%i_syst", j, j));  

      if (j==4){
          TH1D * h4 = (TH1D*) hin[j]->Clone(); h4->Reset("ICES");
          TH1D * h4sy = (TH1D*) hinsy[j]->Clone(); h4sy->Reset("ICES");
            for (int i=1; i<h4->GetXaxis()->GetNbins();i++){
                if (h4->GetXaxis()->GetBinLowEdge(i)<0.699999 || h4->GetXaxis()->GetBinUpEdge(i)>7.00001) continue;
                h4->SetBinContent(i, hin[j]->GetBinContent(i));
                h4->SetBinError(i, hin[j]->GetBinError(i));
                h4sy->SetBinContent(i, hinsy[j]->GetBinContent(i));
                h4sy->SetBinError(i, hinsy[j]->GetBinError(i));
            }
        fout.cd();
        h4->Write();
        h4sy->Write();
      } else {
        fout.cd();
        hin[j]->Write();
        hinsy[j]->Write();
      } 
    }
    return;
}