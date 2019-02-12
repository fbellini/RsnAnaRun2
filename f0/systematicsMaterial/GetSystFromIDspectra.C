enum ETypeUnc { kMaterial= 0, 
		kHadrInt};

Double_t getDaughterPeakPtSymmetricDecay(Double_t motherPt = 1.0);
void GetSystFromIDspectra(Int_t type = ETypeUnc::kMaterial , Bool_t usePeakPt = 0)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TString filename;
  if (type == kMaterial) filename = "COMMON_sysErrRelMaterial.root";
  if (type == kHadrInt) filename = "COMMON_sysErrRelInleasticXSec.root";
  //
  //get pt-dependent systematic uncert due to material budget estimate, for K* analysis
  //
  // define binning in pT    
  Double_t pt[] = {0.0, 0.3, 0.6, 1., 1.5, 2., 2.5, 3., 3.5, 4., 5., 7., 8.};
  Int_t   npt  = sizeof(pt) / sizeof(pt[0]) - 1;
  TAxis *ptbins = new TAxis(npt, pt);
  TFile * systMaterialFile = TFile::Open(filename.Data());
  TH1D * hMaterialPi[2] = {0x0, 0x0}; // pi-, pi+
  TH1D* hTotalSys = new TH1D("hSystVsPtPercentageOfCentral","hSystVsPtPercentageOfCentral; p_{T} (GeV/c); rel. syst. uncert. (%)", npt, pt);
  
  hMaterialPi[1] = (TH1D*) systMaterialFile->Get("hPionMinusMinBiasSysErrRel");
  hMaterialPi[0] = (TH1D*) systMaterialFile->Get("hPionPlusMinBiasSysErrRel");
  
  Int_t NbinsPiP = hMaterialPi[0]->GetNbinsX();
  Int_t NbinsPiM = hMaterialPi[1]->GetNbinsX();
  Printf("Nbins_PiM = %i", NbinsPiP);
  Printf("Nbins_PiP = %i", NbinsPiM);
  
  //loop
  Double_t peakPtPiplus, peakPtPiminus;
  for (Int_t ibin = 0; ibin < npt; ibin++) {

    peakPtPiplus = getDaughterPeakPtSymmetricDecay(pt[ibin]);//hPeakPtKa[0]->GetBinContent(ibin);
    peakPtPiminus = getDaughterPeakPtSymmetricDecay(pt[ibin]);//hPeakPtPi[0]->GetBinContent(ibin);
    Printf("mother pt = %3.2f, pi+ = %3.2f, pi- = %3.2f",pt[ibin], peakPtPiplus, peakPtPiminus);
    //retrieve bin in syst plot corresponding to most probable value of pt of each daughter
    Int_t Mat_PiplusBin = hMaterialPi[0]->GetXaxis()->FindBin(peakPtPiplus);      
    Int_t Mat_PiminusBin = hMaterialPi[1]->GetXaxis()->FindBin(peakPtPiminus);
    //retrieve syst for pt of each daughter
    Float_t systMat_Piplus = hMaterialPi[0]->GetBinContent(Mat_PiplusBin);
    Float_t systMat_Piminus = hMaterialPi[1]->GetBinContent(Mat_PiminusBin);
    hTotalSys->SetBinContent(ibin, systMat_Piplus+systMat_Piminus);
    hTotalSys->SetBinError(ibin, 0.0);
  }
  
  hTotalSys->SetLineColor(kRed); 

  TCanvas * c1 = new TCanvas("c1","c1",800, 600);
  c1->cd();
  hTotalSys->Draw("h");
  hTotalSys->Scale(100.0);
  hTotalSys->GetYaxis()->SetRangeUser(0.0, 20.0);

  hTotalSys->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hTotalSys->GetYaxis()->SetTitle("rel. syst. uncert. (%)");

  TFile * fout = new TFile(Form("syst%s.root",((type == 0)?"Material":"HadrInt")),"recreate");		  
  fout->cd();
  c1->Write();
  hTotalSys->Write();
  return;
}


Double_t getDaughterPeakPtSymmetricDecay(Double_t motherPt)
{
  Double_t totalcm2 = motherPt*motherPt + 0.990*0.990;
  return TMath::Sqrt((totalcm2 - 2*0.139*0.139) / 2.0);
}
