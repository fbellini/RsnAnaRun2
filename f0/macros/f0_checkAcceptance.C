void f0_checkAcceptance(Int_t ptbin = 0,
			Int_t etabin = 0,
			TString nameData = "MCout_acceptanceCheck.root",
			TString listName = "RsnOut_f0_MC_2sTPC_3sTOFveto")
{

    // initial setup
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  gStyle->SetTextFont(42);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetLineWidth(1);
  gStyle->SetLabelSize(0.06,"xyz");
  gStyle->SetLabelOffset(0.005,"yx");
  gStyle->SetTitleSize(0.07,"xyz");
  gStyle->SetTitleOffset(1.1,"y");
  gStyle->SetTitleOffset(1.,"x");
  gStyle->SetEndErrorSize(0); //sets in #of pixels the lenght of the tick at the end of the error bar
  gStyle->SetTitleAlign(33);
  gStyle->SetTitleX(.70);
  gStyle->SetTitleY(.96);
  TGaxis::SetMaxDigits(3);


   // open input file  
  TFile * fileData = TFile::Open(nameData.Data());
  if (!(fileData && fileData->IsOpen())) { 
    Printf("ERROR: cannot open file");
    return; 
  }
    
  TList * listData = (TList*)fileData->Get(listName.Data());
  if (! listData) return;
  
  TH2F * hTrueRes =  (TH2F*) listData->FindObject("RsnTaskF0_Mres_f0");
  if (!hTrueRes) return;

  TH1F * hTrueF0px[14];
  
  TF1 * fg = new TF1("fg","gaus",-0.02,.02); // fit range +- 2 sigma
  TLine l;
  TObjArray arr;
  
  Double_t pT[15] = {0.0, 0.3, 0.6, 1., 1.5, 2., 2.5, 3., 3.5, 4., 5., 6., 7., 8., 10.};
  Double_t eta[6] = {-0.5,-0.3,-0.1, 0.1, 0.3, 0.5};
  
  TH1F *hM = new TH1F("mean", "mean", 11, pT);
  hM->SetMarkerStyle(20);
  hM->SetMarkerSize(1.);
  
  hM->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hM->GetYaxis()->SetTitle("mean(#it{m}_{rec}-#it{m}_{PDG}) (GeV/#it{c}^{2})");
  hM->GetYaxis()->SetRangeUser(-0.010, 0.010);

  
  TH1F *hS = new TH1F("sigma", "sigma", 11, pT);
  hS->SetMarkerStyle(20);
  hS->SetMarkerSize(.5);
  hS->SetMarkerColor(kRed);
  hS->SetLineColor(kRed);
  hS->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hS->GetYaxis()->SetTitle("#sigma(#it{m}_{rec}-#it{m}_{PDG}) (GeV/#it{c}^{2})");
  hS->GetYaxis()->SetRangeUser(0., 0.010);
  
  for (Int_t ibin = 0; ibin < 14; ibin++) {
    Int_t iMinBinPt = hTrueRes->GetYaxis()->FindBin(pT[ibin]);
    Int_t iMaxBinPt = hTrueRes->GetYaxis()->FindBin(pT[ibin + 1]);
    hTrueF0px[ibin] = (TH1F*)hTrueRes->ProjectionX(Form("hRes_%2.1f-%2.1f", pT[ibin], pT[ibin + 1]), iMinBinPt, iMaxBinPt);
    hTrueF0px[ibin]->Fit(fg, "QS");
    hM->SetBinContent(ibin, fg->GetParameter(1));
    hM->SetBinError(ibin, fg->GetParError(1));
    hS->SetBinContent(ibin, fg->GetParameter(2));
    hS->SetBinError(ibin, fg->GetParError(2));
    Printf("%4.3f, ", hS->GetBinContent(ibin));
  }
  
  TCanvas * mean = new TCanvas("mean", "mean", 900, 800);
  mean->cd();
  hM->Draw("same"); 

  TCanvas * sigma = new TCanvas("sigma", "sigma", 900, 800);
  //  hTrueRes->Draw("colz");
  // hM->Draw("same");
  sigma->cd();
  hS->Draw("same"); 

  TFile* resfile = new TFile("F0massResolution.root", "recreate");
  resfile->cd();
  hM->Write();
  hS->Write();

  TCanvas * acc = new TCanvas("acc", "acc", 900, 600);
  acc->Divide(2,2);

  TH3F * hTrueYIn =  (TH3F*) listData->FindObject("RsnTaskF0_trueRap_f0");
  if (!hTrueYIn) return;

  
  TH3F * hMotherY =  (TH3F*) listData->FindObject("RsnTaskF0_motherRap_f0");
  if (!hMotherY) return;
  acc->cd(2);
  hMotherY->Draw("lego");


  TH2F * hAccTrue[14];
  TH2F * hAccMom[14];
  TH2F * hEtaPt[12];
  TH2F * hEtaPtMom[12];

  for (Int_t ibin = 0; ibin < 14; ibin++) {
    Int_t iMinBinPt = hTrueYIn->GetXaxis()->FindBin(pT[ibin]);
    Int_t iMaxBinPt = hTrueYIn->GetXaxis()->FindBin(pT[ibin + 1]);
    hTrueYIn->GetXaxis()->SetRange(pT[ibin], pT[ibin]+1);
    hAccTrue[ibin] = (TH2F*)hTrueYIn->Project3D("zy");
    hAccTrue[ibin]->SetName(Form("hAcc_%2.1f-%2.1f", pT[ibin], pT[ibin + 1]));
    hAccTrue[ibin]->SetTitle(Form("True f_{0}, %2.1f < #it{p}_{T} (GeV/#it{c}) < %2.1f", pT[ibin], pT[ibin + 1]));
    hAccTrue[ibin]->GetYaxis()->SetTitle("#eta");
    hAccTrue[ibin]->GetXaxis()->SetTitle("#it{y}");

    hMotherY->GetXaxis()->SetRange(pT[ibin], pT[ibin]+1);
    hAccMom[ibin] = (TH2F*)hMotherY->Project3D("zy");
    hAccMom[ibin]->SetName(Form("hAccMom_%2.1f-%2.1f", pT[ibin], pT[ibin + 1]));
    hAccMom[ibin]->SetTitle(Form("Generated f_{0}, %2.1f < #it{p}_{T} (GeV/#it{c}) < %2.1f", pT[ibin], pT[ibin + 1]));
    hAccMom[ibin]->GetYaxis()->SetTitle("#eta");
    hAccMom[ibin]->GetXaxis()->SetTitle("#it{y}");
    hAccMom[ibin]->GetXaxis()->SetRangeUser(-0.6, 0.6);
  }
  
  hTrueYIn->GetXaxis()->SetRangeUser(0.01, 15.);
  hMotherY->GetXaxis()->SetRangeUser(0.01, 15.);
  
  for (Int_t ibin = 0; ibin < 5; ibin++) {
    Int_t iMinBinEta = hTrueYIn->GetYaxis()->FindBin(eta[ibin]);
    Int_t iMaxBinEta = hTrueYIn->GetYaxis()->FindBin(eta[ibin + 1]);
    hEtaPt[ibin] = (TH2F*)hTrueYIn->Project3D("zx");
    hEtaPt[ibin]->SetName(Form("hEtaPt_%2.1f-%2.1f", pT[ibin], pT[ibin + 1]));
    hEtaPt[ibin]->SetTitle(Form("True f_{0}, %2.1f < #it{y} < %2.1f", eta[ibin], eta[ibin + 1]));
    hEtaPt[ibin]->GetYaxis()->SetTitle("#eta");
    hEtaPt[ibin]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");

    hEtaPtMom[ibin] = (TH2F*)hMotherY->Project3D("zx");
    hEtaPtMom[ibin]->SetName(Form("hEtaPtMom_%2.1f-%2.1f", pT[ibin], pT[ibin + 1]));
    hEtaPtMom[ibin]->SetTitle(Form("Generated f_{0}, %2.1f < #it{y} < %2.1f", eta[ibin], eta[ibin + 1]));
    hEtaPtMom[ibin]->GetYaxis()->SetTitle("#eta");
    hEtaPtMom[ibin]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  }

  acc->cd(1);
  hAccTrue[ptbin]->Draw("colz");
  acc->cd(2);
  hAccMom[ptbin]->Draw("colz");
  acc->cd(3); gPad->SetLogx();
  hEtaPt[etabin]->Draw("colz");
  acc->cd(4); gPad->SetLogx();
  hEtaPtMom[etabin]->Draw("colz");

}
