void syst_contributions(TString date="27feb18")
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  //set input name
  TString corrFile = "f0_corrSpec_2sTPC_3sTOFveto.root"; //CHANGE-ME
  TString hCorrYieldName = Form("hCorrectedSpectrum"); //CHANGE-ME

  //pp analysis
  Double_t pt[] = {0.5, 1., 1.5, 2., 2.5, 3., 4., 5., 7., 9.};  //CHANGE-ME
  Int_t   npt  = sizeof(pt) / sizeof(pt[0]) - 1;
 //Int_t   ncent  = sizeof(cent) / sizeof(cent[0]) - 1;
  TString centLabel= "INEL";

  //cosmetics
  Color_t color = kBlack;
  Int_t marker = 20;

  //create axis to reproduce the binning
  TAxis *ptbins = new TAxis(npt, pt);
//  TAxis *centbins = new TAxis(ncent, cent);

  //statistical
  TH1F * statunc = new TH1F("statunc","Statistical uncertainty", npt, pt);
  statunc->SetLineWidth(3);
  statunc->SetLineColor(kGreen+2);
  statunc->SetMarkerColor(kGreen+2);
  statunc->SetLineStyle(1);
  statunc->SetMarkerStyle(0);

  //total sys
  TH1F * sum2 = new TH1F("sum2","Total sys. uncert.", npt, pt);
  sum2->SetLineWidth(3);
  sum2->SetLineColor(kRed);
  sum2->SetMarkerColor(kRed);
  sum2->SetLineStyle(1);
  sum2->SetMarkerStyle(0);

  //pt dependent
  TH1F * material = new TH1F("material", "Material Budget", npt, pt);
  TH1F * hadrint = new TH1F("hadrint", "Hadronic inter.", npt, pt);
  TH1F * function = new TH1F("functioni", "Res. bg. fit function", npt, pt);
  TH1F * fit = new TH1F(Form("fit_%i","Fit settings", npt, pt);
  TH1F * pid = new TH1F(Form("pid_%i","PID", npt, pt);

  //pt independent
  TH1F * tracking = new TH1F("tracking","Global tracking",npt, pt);
  tracking->SetLineWidth(2);
  tracking->SetLineColor(kOrange);
  tracking->SetMarkerColor(kOrange);
  tracking->SetLineStyle(9);
  tracking->SetMarkerStyle(0);

  TH1F * trackcuts = new TH1F("trackcuts","Track cuts",npt, pt);
  trackcuts->SetLineWidth(2);
  trackcuts->SetLineColor(kAzure+10);
  trackcuts->SetMarkerColor(kAzure+10);
  trackcuts->SetLineStyle(2);
  trackcuts->SetMarkerStyle(0);

  TFile * fMaterial = TFile::Open(Form("systMaterial.root"));
  TH1F * dummyMT = (TH1F*) fMaterial->Get(Form("hSystVsPtPercentageOfCentral"));
  material = (TH1F*) dummyMT->Clone("material");
  material->SetTitle("Material budget");
  material->SetLineWidth(4);
  material->SetLineColor(kPink+2);
  material->SetMarkerColor(kPink+2);
  material->SetLineStyle(3);
  material->SetMarkerStyle(0);

  TFile * fHadrSyst = TFile::Open(Form("systHadrInt.root"));
  TH1F * dummyHI = (TH1F*) fHadrSyst->Get(Form("hSystVsPtPercentageOfCentral"));
  hadrint = (TH1F*) dummyHI->Clone("hadrint");
  hadrint->SetTitle("Hadronic int.");
  hadrint->SetLineWidth(2);
  hadrint->SetLineColor(kOrange+1);
  hadrint->SetMarkerColor(kOrange+1);
  hadrint->SetLineStyle(1);
  hadrint->SetMarkerStyle(0);


  TString fcnfile = "f0_systematics_functions.root";
  sysHistoName = Form("hFunctions"); //CHANGE ME
  TFile * fFuncSyst = TFile::Open(fcnfile.Data());
  TH1F * dummyF = (TH1F*) fFuncSyst->Get(sysHistoName.Data());
  Printf("Using histo %s / %s",fFuncSyst->GetName(),dummyF->GetName());
  function = (TH1F*) dummyF->Clone("function");
  function->SetTitle("Res. bg. fit function");
  function->SetLineWidth(2);
  function->SetLineColor(kBlue+2);
  function->SetMarkerColor(kBlue+2);
  function->SetLineStyle(1);
  function->SetMarkerStyle(0);

  TString PIDfile = "f0_systematics_PID.root";
  sysHistoName = Form("hPID"); //CHANGE ME
  TFile * fPidSyst = TFile::Open(PIDfile.Data());
  TH1F * dummyP = (TH1F*) fPidSyst->Get(sysHistoName.Data());
  Printf("Using histo %s / %s",fPidSyst->GetName(),dummyP->GetName());
  pid = (TH1F*) dummyP->Clone("PID");
  pid->SetTitle("PID");
  pid->SetLineWidth(3);
  pid->SetLineColor(kBlack);
  pid->SetMarkerColor(kBlack);
  pid->SetLineStyle(2);
  pid->SetMarkerStyle(0);

  //ITS-TPC matching efficiency sys uncertainty in pp 5 TeV LHC15n
  //https://twiki.cern.ch/twiki/bin/view/ALICE/AliDPGtoolsTrackSystematicUncertaintyBookkeping
  //https://twiki.cern.ch/twiki/bin/view/ALICE/AliDPGtoolsTrackSystematicUncertainty
  //Summary of recommended values of the systematic error for the single track:
  //1% for 0.5<pT<2GeV/c, 2% for 2<pT<7, 1% for 7<pT<10, 2% for for 10<pT<15
  //When considering a resonance in two-body decay, multiply the uncert. on the single track times 2
  //as the two daughters are considered independently tracked thus independently affected by the systemtics
  Double_t tracking_1trkBelow2GeV = 0.01;
  Double_t tracking_1trk2to7GeV = 0.02;
  Double_t tracking_1trkAbove7GeV = 0.02;

  //track cuts systematics uncertainty inherited from the analysis of pi,K,p and pt-independent
  //When considering a resonance in two-body decay, multiply the uncert. on the single track times 2
  //as the two daughters are considered independently tracked thus independently affected by the systemtics
  Double_t trackCuts_rsn = 0.025;

  for (Int_t ii = 0;ii<npt;ii++){
    Int_t ibin = ii+1;
    if (pt[ibin] < 4.0) {
      tracking->SetBinContent(ibin, tracking_1trkBelow2GeV*2.0*100.0);
    } else if (pt[ibin] < 14.0){
      tracking->SetBinContent(ibin, tracking_1trk2to7GeV*2.0*100.0);
    } else
      tracking->SetBinContent(ibin, tracking_1trkAbove7GeV*2.0*100.0);

    trackcuts->SetBinContent(ibin, trackCuts_rsn*2.0*100.0);
  }

  //sum correlated contributions in quadrature
  Double_t syst_ptFullyCorr_sum2 = trackCuts_rsn*trackCuts_rsn;

  //estimate uncertainty per each pt bin
  for (Int_t ii = 0;ii<npt;ii++){
    Int_t ibin = ii+1;

    Double_t tracking_sys = tracking->GetBinContent(ibin)/100.;
    Double_t material_syst = dummyMT->GetBinContent(ibin)/100.;
    Double_t hadrint_syst = dummyHI->GetBinContent(ibin)/100.;
    Double_t range_syst = dummyR->GetBinContent(ibin)/100.;
    Double_t func_syst = dummyF->GetBinContent(ibin)/100.;
    Double_t pid_syst = dummyP->GetBinContent(ibin)/100.;

    //tracking
    Double_t ptuncorr2 = tracking_sys*tracking_sys;

   //track cuts
    ptuncorr2+=(trackCuts_rsn*trackCuts_rsn);

    //fit range
    if (range_syst>0.0) {
      ptuncorr2+=range_syst*range_syst;
    }

    //function res bg
    if (func_syst>0.0){
      ptuncorr2+=func_syst*func_syst;
   }

    //PID
    if (pid_syst>0.0) {
      ptuncorr2+=pid_syst*pid_syst;
    }

   //material budget
    if (material_syst>0.0) {
      ptuncorr2+=material_syst*material_syst;
    }

    //hadronic interaction cross section
    if (hadrint_syst>0.0) {
      ptuncorr2+=hadrint_syst*hadrint_syst;
    }

    Double_t totsyst = TMath::Sqrt(syst_ptFullyCorr_sum2 + ptuncorr2);
    sum2->SetBinContent(ibin, totsyst*100);
    Printf("bin %i tot. perc. %6.4f", ibin, totsyst*100.);

    Double_t totsystUncorr = TMath::Sqrt(ptuncorr2);
    //Printf("totsystUncorr = %f", totsystUncorr);
    sum2_uncorr->SetBinContent(ibin, totsystUncorr*100.);
  }

  //assign syst err to data
  TFile * fdata = TFile::Open(Form("%s", corrFile.Data()));
  TH1F * data = (TH1F*) fdata->Get(Form("%s",hCorrYieldName.Data()));
  TH1F * data_Wsyst = (TH1F*) data->Clone(Form("%s_syst",hCorrYieldName.Data()));
  TH1F * data_Wsyst_uncorr = (TH1F*) data->Clone(Form("%s_syst_uncorr",hCorrYieldName.Data()));
  TH1F * data_Wsyst_Wstat = (TH1F*) data->Clone(Form("%s_syst_stat",hCorrYieldName.Data()));

  for (Int_t ii = 1;ii<npt+1;ii++){
    Int_t ibin = ii;
    Double_t yd = data->GetBinContent(ibin);
    Double_t yd_stat = data->GetBinError(ibin);
    Double_t perc = sum2->GetBinContent(ibin);
    Double_t yd_syst = yd*perc/100.;
    statunc->SetBinContent(ibin, (yd>0? (yd_stat*100.0/yd) : 0.0));
    data_Wsyst->SetBinContent(ibin,yd);
    data_Wsyst->SetBinError(ibin, yd_syst);
    data_Wsyst_Wstat->SetBinContent(ibin,yd);
    data_Wsyst_Wstat->SetBinError(ibin,TMath::Sqrt(yd_syst*yd_syst+yd_stat*yd_stat));
    Printf("bin %i     yield = %e     syst err = %e    stat err = %e", ibin, yd, yd_syst, yd_stat);

    //uncorrelated wrt pt
    Double_t perc_uncorr = sum2_uncorr->GetBinContent(ibin);
    Double_t yd_syst_uncorr = yd*perc_uncorr/100.;
    data_Wsyst_uncorr->SetBinContent(ibin,yd);
    data_Wsyst_uncorr->SetBinError(ibin, yd_syst_uncorr);
  }

  //systematic uncertainty plot
  TCanvas *cs=new TCanvas("cs","Systematic uncertainty vs #it{p}_{T}", 750,600);
  cs->cd();
  sum2->SetTitle("Total syst. uncert.");
  sum2->GetYaxis()->SetTitle("relative uncert. (%)");
  sum2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  sum2->GetYaxis()->SetRangeUser(0.01, 40);
  sum2->GetXaxis()->SetRangeUser(0.0, 6.9);

  sum2->Draw();
  tracking->Draw("same");
  range->Draw("same");
  trackcuts->Draw("same");
  function->Draw("same");
  material->Draw("same");
  pid->Draw("same");
  hadrint->Draw("same");

  TLegend * autolegry = (TLegend*)gPad->BuildLegend(0.25,0.65,0.88,0.88);
  autolegry->SetFillColor(kWhite);
  autolegry->SetLineColor(kWhite);
  autolegry->SetTextFont(42);
  autolegry->SetNColumns(2);
  data->SetMarkerColor(color);
  data->SetLineColor(color);
  data->SetMarkerStyle(marker);
  data->SetMarkerSize(0.7);
  data->SetLineWidth(2);
  data->SetTitle(Form("%s (stat. uncert.)",centLabel.Data()));
  data_Wsyst->SetMarkerColor(color[ipid]);
  data_Wsyst->SetFillStyle(0); //Color(color2[icent]);
  data_Wsyst->SetLineWidth(1); //Color(color2[icent]);
  data_Wsyst->SetLineColor(color[ipid]);
  data_Wsyst->SetMarkerStyle(0);
  data_Wsyst->SetOption("E2");
  data_Wsyst->SetTitle(Form("%s (syst. uncert.)",centLabel.Data()));

  TString imagefilename = Form("summaryAllSystUncert") ; //CHANGE ME if you want
  cs->SaveAs(Form("%s.png", imagefilename.Data()));
  cs->SaveAs(Form("%s.eps", imagefilename.Data()));

  //make-up
  data_Wsyst_uncorr->SetMarkerColor(color);
  data_Wsyst_uncorr->SetFillStyle(0); //Color(color2[icent]);
  data_Wsyst_uncorr->SetLineWidth(1); //Color(color2[icent]);
  data_Wsyst_uncorr->SetLineColor(color);
  data_Wsyst_uncorr->SetMarkerStyle(0);
  data_Wsyst_uncorr->SetOption("E2");
  data_Wsyst_uncorr->SetTitle(Form("%s (syst. uncert., #it{p}_{T}-uncorr.)",centLabel.Data()));

  data_Wsyst_Wstat->SetMarkerColor(color);
  data_Wsyst_Wstat->SetFillStyle(0); //Color(color2[icents]);
  data_Wsyst_Wstat->SetLineWidth(1); //Color(color2[icent]);
  data_Wsyst_Wstat->SetLineColor(color);
  data_Wsyst_Wstat->SetMarkerStyle(0);
  data_Wsyst_Wstat->SetTitle(Form("%s (#sqrt{syst^{2}+stat^{2}})",centLabel.Data()));

  TCanvas *cunc=new TCanvas("cunc","Summary of uncertainty vs p_{t}", 750,600);
  cunc->cd();
  sum2->Draw();
  sum2_uncorr->Draw("HIST same");
  statunc->Draw("HIST same");
  TLegend * autolegry2 = (TLegend*)gPad->BuildLegend(0.25,0.65,0.88,0.88);
  autolegry2->SetFillColor(kWhite);
  autolegry2->SetLineColor(kWhite);
  autolegry2->SetTextFont(42);
  autolegry2->SetNColumns(1);
  autolegry2->Draw();

  //save to out file
  TString outfilename = Form("finalWsys.root"); //CHANGE ME if you want
  TFile * fout = new TFile(outfilename.Data(),"recreate");
  fout->cd();
  material->Write();
  tracking->Write();
  trackcuts->Write();
  pid->Write();
  range->Write();
  function->Write();
  hadrint->Write();
  statunc->Write();
  data->Write();
  sum2->Write();
  sum2_uncorr->Write();
  data_Wsyst->Write();
  data_Wsyst_Wstat->Write();
  cunc->Write();
  cs->Write();
  fout->Close();

  return;

}
