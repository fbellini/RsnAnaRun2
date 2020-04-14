enum ETypeUnc { kMaterial= 0, 
		kHadrInt = 1,
		kITSTPCmatch = 2,
		kTrackCuts};

TH1F * getITSTPCMatchingUncertXeXe();
void getTrackCutsUncertaintyFromChargedXeXe(Int_t ncentbins = 4);
TString GetTypeName(Int_t type = 0);

Double_t getDaughterPeakPtSymmetricDecay(Double_t motherPt = 1.0);
Double_t getDaughterPeakPtPhaseSpace(Int_t motherPtBin = 0, Int_t daughId = 1);

void GetSystFromIDspectra(Int_t type = ETypeUnc::kMaterial , Bool_t usePhaseSpace = 0, Int_t icent = -1)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TString filename;
  if (type == kMaterial) filename = "/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/sysMaterial/COMMON_sysErrRelMaterial.root";
  if (type == kHadrInt) filename = "/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/sysMaterial/COMMON_sysErrRelInleasticXSec.root";
  if (type == kITSTPCmatch) filename = "/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/sysMaterial/COMMON_sysITSTPCmatchingXeXe.root";
  if (type == kTrackCuts) filename = "COMMON_sysTrackCutsXeXe.root";
  //
  //get pt-dependent systematic uncert due to material budget estimate, for phi analysis
  //
  // define binning in pT    
  Double_t pt[] = {0.0, 0.3, 0.5, 0.7, 0.9, 1.10, 1.30, 1.50, 2.00, 3.00, 4.00, 5.0, 7.0, 10.0, 11.};
  Int_t   npt  = sizeof(pt) / sizeof(pt[0]) - 1;
  TAxis *ptbins = new TAxis(npt, pt);
  TFile * systMaterialFile = TFile::Open(filename.Data());
  TH1D * hMaterialKa[2] = {0x0, 0x0}; // ka-, ka+
  TH1D* hTotalSys = new TH1D("hSystVsPtPercentageOfCentral","hSystVsPtPercentageOfCentral; p_{T} (GeV/c); rel. syst. uncert. (%)", npt, pt);
  
  if (icent<0) {
    hMaterialKa[1] = (TH1D*) systMaterialFile->Get(Form("hKaonMinusMinBiasSysErrRel"));
    hMaterialKa[0] = (TH1D*) systMaterialFile->Get(Form("hKaonPlusMinBiasSysErrRel"));
  } else {
    hMaterialKa[1] = (TH1D*) systMaterialFile->Get(Form("hKaonMinusMinBiasSysErrRel_%i", icent));
    hMaterialKa[0] = (TH1D*) systMaterialFile->Get(Form("hKaonPlusMinBiasSysErrRel_%i", icent));
  }
  if (!hMaterialKa[0]) {Printf("We have a problem"); return;}
  Int_t NbinsKaP = hMaterialKa[0]->GetNbinsX();
  Int_t NbinsKaM = hMaterialKa[1]->GetNbinsX();
  Printf("Nbins_KaM = %i", NbinsKaP);
  Printf("Nbins_KaP = %i", NbinsKaM);
  
  //loop
  Double_t peakPtKaplus, peakPtKaminus;
  for (Int_t ibin = 0; ibin < npt; ibin++) {
    if (usePhaseSpace) {
      peakPtKaplus = getDaughterPeakPtPhaseSpace(ibin, 1);
      peakPtKaminus = getDaughterPeakPtPhaseSpace(ibin, 2);
    } else {
      peakPtKaplus = getDaughterPeakPtSymmetricDecay(pt[ibin]);//hPeakPtKa[0]->GetBinContent(ibin);
      peakPtKaminus = getDaughterPeakPtSymmetricDecay(pt[ibin]);//hPeakPtKa[0]->GetBinContent(ibin);
    }
    Printf("mother pt = %3.2f, ka+ = %3.2f, ka- = %3.2f",pt[ibin], peakPtKaplus, peakPtKaminus);
    //retrieve bin in syst plot corresponding to most probable value of pt of each daughter
    Int_t Mat_KaplusBin = hMaterialKa[0]->GetXaxis()->FindBin(peakPtKaplus);      
    Int_t Mat_KaminusBin = hMaterialKa[1]->GetXaxis()->FindBin(peakPtKaminus);
    //retrieve syst for pt of each daughter
    Float_t systMat_Kaplus = hMaterialKa[0]->GetBinContent(Mat_KaplusBin);
    Float_t systMat_Kaminus = hMaterialKa[1]->GetBinContent(Mat_KaminusBin);
    hTotalSys->SetBinContent(ibin, systMat_Kaplus+systMat_Kaminus);
    hTotalSys->SetBinError(ibin, 0.0);
    Printf(">>> sys uncert %f", systMat_Kaplus+systMat_Kaminus);
  }
  
  hTotalSys->SetLineColor(kRed); 

  TCanvas * c1 = new TCanvas("c1","c1",800, 600);
  c1->cd();
  hTotalSys->Draw("h");
  hTotalSys->Scale(100.0);
  hTotalSys->GetYaxis()->SetRangeUser(0.0, 20.0);

  hTotalSys->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hTotalSys->GetYaxis()->SetTitle("rel. syst. uncert. (%)");

  TString foutname = ((icent<0)?Form("syst%s.root", GetTypeName(type).Data()) : Form("syst%s%i.root", GetTypeName(type).Data(), icent));
  TFile * fout = new TFile(foutname.Data(),"recreate");		  
  fout->cd();
  c1->Write();
  hTotalSys->Write();
  return;
}

TString GetTypeName(Int_t type)
{
  switch (type) {
  case 0:  
    return "Material";
  case 1:
    return "HadrInt";
  case 2:
    return "ITSTPCmatch";
  case 3:  
    return "TrackCuts";
  }
  return "";
}

Double_t getDaughterPeakPtSymmetricDecay(Double_t motherPt)
{
  Double_t totalcm2 = motherPt*motherPt + 1.01995*1.01995;
  return TMath::Sqrt((totalcm2 - 2*0.4937*0.4937) / 2.0);
}

Double_t getDaughterPeakPtPhaseSpace(Int_t motherPtBin, Int_t daughId)
{
  TFile * fin = TFile::Open("phi_PhaseSpace.root");
  TString histname = Form("phi_1std_pt%i", motherPtBin);
  if (daughId==2) histname = Form("phi_2nsd_pt%i", motherPtBin);
  TH1D * hin = (TH1D*) fin->Get(histname.Data());
  if (!hin) return 0.0;
  
  return hin->GetMean();
}

TH1F* getITSTPCMatchingUncertXeXe()
{
  //from DPG: https://twiki.cern.ch/twiki/bin/viewauth/ALICE/AliDPGtoolsTrackSystematicUncertaintyBookkeping
  TH1F * hist = new TH1F("hITSTPCmatchingSys", "hITSTPCmatchingSys", 140, 0., 14.); 
  for (int j=0; j<hist->GetNbinsX()+1; j++){
    if (hist->GetXaxis()->GetBinUpEdge(j)<1.0) hist->SetBinContent(j, 0.032);
    else if (hist->GetXaxis()->GetBinUpEdge(j)<2.0) hist->SetBinContent(j, 0.043);
    else if (hist->GetXaxis()->GetBinUpEdge(j)<3.0) hist->SetBinContent(j, 0.055);
    else if (hist->GetXaxis()->GetBinUpEdge(j)<4.0) hist->SetBinContent(j, 0.050);
    else if (hist->GetXaxis()->GetBinUpEdge(j)<5.0) hist->SetBinContent(j, 0.047);
    else if (hist->GetXaxis()->GetBinUpEdge(j)<6.0) hist->SetBinContent(j, 0.035);
    else if (hist->GetXaxis()->GetBinUpEdge(j)<7.0) hist->SetBinContent(j, 0.037);
    else if (hist->GetXaxis()->GetBinUpEdge(j)<8.0) hist->SetBinContent(j, 0.040);
    else if (hist->GetXaxis()->GetBinUpEdge(j)<9.0) hist->SetBinContent(j, 0.042);
    else if (hist->GetXaxis()->GetBinUpEdge(j)<10.0) hist->SetBinContent(j, 0.033);
    else if (hist->GetXaxis()->GetBinUpEdge(j)<11.0) hist->SetBinContent(j, 0.052);
    else if (hist->GetXaxis()->GetBinUpEdge(j)<12.0) hist->SetBinContent(j, 0.025);
    else if (hist->GetXaxis()->GetBinUpEdge(j)<14.2) hist->SetBinContent(j, 0.040);
  }
  TFile * fout = new TFile("COMMON_sysITSTPCmatchingXeXe.root", "recreate");
  fout->cd();
  hist->Write("hKaonMinusMinBiasSysErrRel");
  hist->Write("hKaonPlusMinBiasSysErrRel");
  hist->Draw();
  return hist;
}

void getTrackCutsUncertaintyFromChargedXeXe(Int_t ncentbins)
{
  //this method is designed for centrality bins as 0-10, 10-30, 30-60 and 60-90%
  TFile * fin[10];   TH1D * hin[10];

  fin[0] = TFile::Open("/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/sysMaterial/systContributions-XeXe-5TeV_c05.root");
  hin[0] = (TH1D *) fin[0]->Get("uncTrkCuts");
  fin[1] = TFile::Open("/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/sysMaterial/systContributions-XeXe-5TeV_c510.root");
  hin[1] = (TH1D *) fin[1]->Get("uncTrkCuts");
  for (int ic = 1; ic<9; ic++){
    fin[ic+1] = TFile::Open(Form("/Users/fbellini/alice/resonances/RsnAnaRun2/phiXeXe/sysMaterial/systContributions-XeXe-5TeV_c%i%i.root", 10*ic, 10*(ic+1)));
    hin[ic+1] = (TH1D *) fin[ic+1]->Get("uncTrkCuts");
  }
  Color_t color[4] = {kRed+1, kOrange, kSpring+5, kBlue+1};
  Float_t weight[10] = {1167.*0.5, 939.*0.5, 706., 478., 315., 198., 118., 64.7, 32.0, 13.3};
  TFile * fout = new TFile("COMMON_sysTrackCutsXeXe.root", "recreate");
   
  for (int ic = 0; ic<ncentbins; ic++){  
    TH1D * hist = 0x0;
    if (ncentbins==3) {
      if (ic==0){
	hin[0]->Scale(weight[0]);
	hin[1]->Scale(weight[1]);
	hin[2]->Scale(weight[2]);
	hin[3]->Scale(weight[3]);
	hist = (TH1D*)hin[0]->Clone("hSys");
	hist->Add(hin[1]);
	hist->Add(hin[2]);
	hist->Add(hin[3]);
	hist->Scale(1./(weight[0]+weight[1]+weight[2]+weight[3]));
      } else if (ic==1){
	hin[4]->Scale(weight[4]);
	hin[5]->Scale(weight[5]);
	hin[6]->Scale(weight[6]);
	hist = (TH1D*)hin[4]->Clone("hSys");
	hist->Add(hin[5]);
	hist->Add(hin[6]);
	hist->Scale(1./(weight[4]+weight[5]+weight[6]));
      }	else if (ic==2) {
	hin[7]->Scale(weight[7]);
	hin[8]->Scale(weight[8]);
	hin[9]->Scale(weight[9]);
	hist = (TH1D*)hin[7]->Clone("hSys");
	hist->Add(hin[8]);
	hist->Add(hin[9]);
	hist->Scale(1./(weight[7]+weight[8]+weight[9]));
      }
    } else {
      if (ic==0){
	hin[0]->Scale(weight[0]);
	hin[1]->Scale(weight[1]);
	hist = (TH1D*)hin[0]->Clone("hSys");
	hist->Add(hin[1]);
	hist->Add(hin[2]);
	hist->Scale(1./(weight[0]+weight[1]));
      } else if (ic==1){
	hin[2]->Scale(weight[2]);
	hin[3]->Scale(weight[3]);
	hist = (TH1D*)hin[0]->Clone("hSys");
	hist->Add(hin[2]);
	hist->Add(hin[3]);
	hist->Scale(1./(weight[2]+weight[3]));
      } else if (ic==2) {
	hin[4]->Scale(weight[4]);
	hin[5]->Scale(weight[5]);
	hin[6]->Scale(weight[6]);
	hist = (TH1D*)hin[4]->Clone("hSys");
	hist->Add(hin[5]);
	hist->Add(hin[6]);
	hist->Scale(1./(weight[4]+weight[5]+weight[6]));
      }	else if (ic==3) {
	hin[7]->Scale(weight[7]);
	hin[8]->Scale(weight[8]);
	hin[9]->Scale(weight[9]);
	hist = (TH1D*)hin[7]->Clone("hSys");
	hist->Add(hin[8]);
	hist->Add(hin[9]);
	hist->Scale(1./(weight[7]+weight[8]+weight[9]));
      }
    }      
    hist->Scale(0.01);
    hist->SetLineColor(color[ic]);
    fout->cd();
    hist->Write(Form("hKaonMinusMinBiasSysErrRel_%i", ic));
    hist->Write(Form("hKaonPlusMinBiasSysErrRel_%i", ic));
  }
  return;
}
