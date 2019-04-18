#include "/Users/fbellini/alice/macros/ReadInputFromFile.C"
#include "/Users/fbellini/alice/macros/SetStyle.C"
#include "/Users/fbellini/alice/macros/Beautify.C"
void GetRatioSysStatSq_XeXe(Double_t &x, Double_t &y, Int_t phiCent, TH1F * kaStat, TH1F * kaSys, TGraphErrors * phiStat, TGraphErrors * phiSys);
void GetRatioSysStatSq_PbPb(Double_t &x, Double_t &y, Int_t phiCent, TH1F * kaStat, TH1F * kaSys, TGraphErrors * phiStat, TGraphErrors * phiSys);
void GetRatio_XeXe(Double_t &x, Double_t &y, Int_t phiCent, TH1F * kaStat, TGraphErrors * phiStat);
void GetRatio_PbPb(Double_t &x, Double_t &y, Int_t phiCent, TH1F * kaStat, TGraphErrors * phiStat);
TGraphErrors * GetPhi2Ka_XeXe544(Bool_t sys = 0);
TGraphErrors * GetPhi2Ka_PbPb502(Bool_t sys = 0);
TGraphErrors * GetRatio_pp13(Int_t sys = 0);
TGraphAsymmErrors *  GetRatio_PbPb276(Bool_t sys = 0);
TGraphAsymmErrors *  GetRatio_pPb502(Bool_t sys = 0);

void Ratio2Kaons(Bool_t plotThermal = 0)
{
  SetStyle();

  //GetRatio Xe-Xe 5.44 TeV preliminary
  TGraphErrors * gPhi2KaXeXe =  (TGraphErrors *) GetPhi2Ka_XeXe544(0);
  TGraphErrors * gPhi2KaXeXe_syst = (TGraphErrors *) GetPhi2Ka_XeXe544(1);
  BeautifyGraph(gPhi2KaXeXe, kRed, kRed, 0, 1, 1, 33, 1.7); 
  BeautifyGraph(gPhi2KaXeXe_syst, kRed, kRed, 0, 1, 1, 33, 1.7); 
 
  //GetRatio Pb-Pb 5.02 TeV preliminary
  TGraphErrors * gPhi2KaPbPb = (TGraphErrors *) GetPhi2Ka_PbPb502(0);
  TGraphErrors * gPhi2KaPbPb_syst = (TGraphErrors *) GetPhi2Ka_PbPb502(1);
  BeautifyGraph(gPhi2KaPbPb, kBlack, kBlack, 0, 1, 1, 21, 1.3); 
  BeautifyGraph(gPhi2KaPbPb_syst, kBlack, kBlack, 0, 1, 1, 21, 1.3); 
  
  //Get ratio PbPb 2.76 TeV published
  TGraphAsymmErrors *gPhi2Ka_PbPb276TeV_stat = (TGraphAsymmErrors*) GetRatio_PbPb276(0);
  TGraphAsymmErrors *gPhi2Ka_PbPb276TeV_syst = (TGraphAsymmErrors*) GetRatio_PbPb276(1);
  BeautifyGraphAsymmErrors(gPhi2Ka_PbPb276TeV_stat, kGray+1, kGray+1, 0, 1, 1, 25, 1.3); 
  BeautifyGraphAsymmErrors(gPhi2Ka_PbPb276TeV_syst, kGray+1, kGray+1, 0, 1, 1, 25, 1.3); 

  //GetRatio pPb 5.02 TeV published
  TGraphAsymmErrors * gPhi2Ka_pPb502TeV_stat = (TGraphAsymmErrors *) GetRatio_pPb502(0);
  TGraphAsymmErrors * gPhi2Ka_pPb502TeV_syst = (TGraphAsymmErrors *) GetRatio_pPb502(1);
  gPhi2Ka_pPb502TeV_stat->SetName("prel_phi2ka_pPb502TeV_stat");
  gPhi2Ka_pPb502TeV_syst->SetName("prel_phi2ka_pPb502TeV_syst");
  BeautifyGraphAsymmErrors(gPhi2Ka_pPb502TeV_stat, kBlue+1, kBlue+1, 0, 1, 1, 28, 1.5); 
  BeautifyGraphAsymmErrors(gPhi2Ka_pPb502TeV_syst, kBlue+1, kBlue+1, 0, 1, 1, 28, 1.5); 

  //GetRatio in pp 13 TeV preliminary
  TGraphAsymmErrors * gPhi2Ka_pp13_stat = (TGraphAsymmErrors *) GetRatio_pp13(0);
  TGraphAsymmErrors * gPhi2Ka_pp13_syst = (TGraphAsymmErrors *) GetRatio_pp13(1);
  TGraphAsymmErrors * gPhi2Ka_pp13_uncor = (TGraphAsymmErrors *) GetRatio_pp13(2);
  BeautifyGraphAsymmErrors(gPhi2Ka_pp13_stat, kMagenta+1, kMagenta+1, 0, 1, 1, 20, 1.3); 
  BeautifyGraphAsymmErrors(gPhi2Ka_pp13_syst, kMagenta+1, kMagenta+1, 0, 1, 1, 20, 1.3); 
  BeautifyGraphAsymmErrors(gPhi2Ka_pp13_uncor, kMagenta+1, kMagenta+1, 0, 1, 1, 20, 1.3); 
  gPhi2Ka_pp13_uncor->SetFillStyle(1001);
  
  TH1F* hf = new TH1F("frame","", 3000, 1.5, 3.E3);
  hf->SetMinimum(1.e-3);
  hf->SetMaximum(4.e-2);
  hf->SetXTitle("#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta} | < 0.5}");
  hf->SetLineColor(1);
  hf->GetXaxis()->SetTitleSize(0.06);
  hf->GetXaxis()->SetLabelSize(0.045);
  hf->GetXaxis()->SetTitleOffset(1.2);
  hf->SetNdivisions(509,"x");
  hf->GetXaxis()->SetMoreLogLabels(1);
  hf->SetYTitle("2#phi/(K^{+}+K^{-})");
  hf->GetYaxis()->SetTitleOffset(1.2);
  hf->GetYaxis()->SetTitleSize(0.06);
  hf->GetYaxis()->SetLabelSize(0.045);
  hf->GetYaxis()->SetTitleOffset(1.2);
  hf->SetNdivisions(507,"y");
  hf->GetYaxis()->SetRangeUser(0.0, 0.25);
  hf->GetXaxis()->SetRangeUser(1.9, 3000.0);

  gPhi2KaXeXe->SetName("gPhi2Ka_XeXe544_stat");
  gPhi2KaXeXe_syst->SetName("gPhi2Ka_XeXe544_syst");
  gPhi2KaXeXe->SetTitle("Xe-Xe #sqrt{#it{s}_{NN}} = 5.44 TeV");
  gPhi2KaXeXe_syst->SetTitle("Xe-Xe #sqrt{#it{s}_{NN}} = 5.44 TeV");
  gPhi2KaPbPb->SetName("gPhi2Ka_PbPb502_stat");
  gPhi2KaPbPb_syst->SetName("gPhi2Ka_PbPb502_syst");
  gPhi2KaPbPb->SetTitle("Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  gPhi2KaPbPb_syst->SetTitle("Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  gPhi2Ka_PbPb276TeV_stat->SetName("gPhi2Ka_PbPb276_stat");
  gPhi2Ka_PbPb276TeV_syst->SetName("gPhi2Ka_PbPb276_syst");
  gPhi2Ka_PbPb276TeV_stat->SetTitle("Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV");
  gPhi2Ka_PbPb276TeV_syst->SetTitle("Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV");

  TPaveText *titletextunc = new TPaveText(0.2,0.83,0.5,0.93,"brNDC");
  titletextunc->SetBorderSize(0);
  titletextunc->SetFillColor(0);
  titletextunc->SetFillStyle(0);
  titletextunc->SetTextAlign(12);
  titletextunc->SetTextSize(0.04);
  titletextunc->SetTextFont(42);
  titletextunc->InsertText("V0 multiplicity classes");
  titletextunc->InsertText("Uncertainties: stat.(bar), syst.(box)");  

  // TLegend * leg2 = new TLegend(0.2, 0.35, 0.9, 0.45, "ALICE, V0 multiplicity classes");
  // myLegendSetUp(leg2, 0.04);
  // leg2->SetNColumns(2);
  // leg2->AddEntry(gPhi2Ka_PbPb276TeV_stat, "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV", "pf"); // [PRC 91 (2015) 024609]
  // leg2->AddEntry(gPhi2Ka_pPb502TeV_stat, "p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV", "pf"); //  [EPJC 76 (2016) 245]

  // TLegend * leg = new TLegend(0.2, 0.2, 0.9, 0.35, "#bf{Preliminary}");
  // myLegendSetUp(leg, 0.04);
  // leg->SetNColumns(2);
  // leg->AddEntry(gPhi2KaXeXe, "Xe-Xe #sqrt{#it{s}_{NN}} = 5.44 TeV", "pf");
  // leg->AddEntry(gPhi2KaPbPb, "Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV", "pf");
  // leg->AddEntry(gPhi2Ka_pp13_stat, "pp #sqrt{#it{s}} = 13 TeV", "pf");

  //Thermal model, A.Andronic et al.
  //https://arxiv.org/abs/1710.09425
  Float_t thermalModelPredGSI = 0.12;
  TLine * thermalGSI = new TLine(700., thermalModelPredGSI, 2500., thermalModelPredGSI);
  thermalGSI->SetLineStyle(1);
  thermalGSI->SetLineColor(kSpring-5);
  thermalGSI->SetLineWidth(5);

  TLegend * legModel = new TLegend(0.2, 0.7, 0.45, 0.8);
  myLegendSetUp(legModel, 0.04);
  legModel->AddEntry(thermalGSI, "GSI-Heidelberg", "l");
  legModel->AddEntry(thermalGSI, "#it{T}_{ch} = 156 MeV", "");
  
  TLegend * leg = new TLegend(0.2, 0.18, 0.45, 0.38, "#bf{ALICE Preliminary}");
  myLegendSetUp(leg, 0.04);
  leg->AddEntry(gPhi2KaXeXe, "Xe-Xe #sqrt{#it{s}_{NN}} = 5.44 TeV", "p");
  leg->AddEntry(gPhi2KaPbPb, "Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV", "p");
  //leg->AddEntry(gPhi2Ka_pp13_stat, "pp #sqrt{#it{s}} = 13 TeV", "p");
   
  TLegend * leg2 = new TLegend(0.55, 0.23, 0.8, 0.38, "ALICE");
  myLegendSetUp(leg2, 0.04);
  leg2->AddEntry(gPhi2Ka_PbPb276TeV_stat, "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV", "p"); // []
  //leg2->AddEntry(gPhi2Ka_PbPb276TeV_stat, "PRC 91 (2015) 024609", "");
  leg2->AddEntry(gPhi2Ka_pPb502TeV_stat, "p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV", "p"); //  [EPJC 76 (2016) 245]
  //leg2->AddEntry(gPhi2Ka_pPb502TeV_stat, "EPJ C76 (2016) 245","");
  
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.17);
  gStyle->SetPadRightMargin(0.05);
  TCanvas * c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();
  gPad->SetLogx();
  hf->Draw();
  // gPhi2Ka_PbPb276TeV_syst->Draw("e2 same");
  // gPhi2Ka_PbPb276TeV_stat->Draw("pz same");
  // gPhi2Ka_pPb502TeV_syst->Draw("e2 same");
  // gPhi2Ka_pPb502TeV_stat->Draw("pz same");
  //gPhi2Ka_pp13_uncor->Draw("e2 same");
  // gPhi2Ka_pp13_syst->Draw("e2 same");
  // gPhi2Ka_pp13_stat->Draw("pz same");
  gPhi2KaPbPb_syst->Draw("e2 same");
  gPhi2KaPbPb->Draw("pz same");
  gPhi2KaXeXe_syst->Draw("e2 same");
  gPhi2KaXeXe->Draw("pz same");
  leg->Draw();
  //leg2->Draw();
  titletextunc->Draw();
  if (plotThermal) {
   thermalGSI->Draw();
   legModel->Draw();
  }
  
  TFile * fout = new TFile("Preliminary_PhiToKaonRatio.root","recreate");
  fout->cd();
  gPhi2Ka_PbPb276TeV_syst->Write();
  gPhi2Ka_PbPb276TeV_stat->Write();
  gPhi2KaPbPb_syst->Write();
  gPhi2KaPbPb->Write();
  gPhi2KaXeXe_syst->Write();
  gPhi2KaXeXe->Write();
  gPhi2Ka_pPb502TeV_stat->Write();
  gPhi2Ka_pPb502TeV_syst->Write();
  gPhi2Ka_pp13_stat->Write();
  gPhi2Ka_pp13_uncor->Write();
  gPhi2Ka_pp13_syst->Write();
  c1->Write();
  c1->Print(Form("phiToKaon_vsdNdeta%s.pdf", (plotThermal? "_wThermal" : "")));
  c1->Print(Form("phiToKaon_vsdNdeta%s.eps", (plotThermal? "_wThermal" : "")));
  fout->Close();
  return;
}

void GetRatioSysStatSq_PbPb(Double_t &x, Double_t &y, Int_t phiCent, TH1F * kaStat, TH1F * kaSys, TGraphErrors * phiStat, TGraphErrors * phiSys)
{
  if (!kaStat  || !kaSys || !phiStat || !phiSys) return 0.0;
  Printf(":::: PhiCent = %i", phiCent);
  Double_t kaYmerged = 0.0;
  Double_t kaYmergedErr = 0.0;
  Double_t phiY = 0.0;
  Double_t phiYerr = 0.0;

  //centrality sums for kaon yields
  if (phiCent == 0){ //phi 0-10%, kaons 0-5%, 5-10%
    kaYmerged = kaStat->GetBinContent(1)*0.5 + kaStat->GetBinContent(2)*0.5;
    kaYmergedErr  = TMath::Sqrt(TMath::Power(kaStat->GetBinError(1), 2.0) + TMath::Power(kaSys->GetBinError(1), 2.0))*0.5;
    kaYmergedErr  += TMath::Sqrt(TMath::Power(kaStat->GetBinError(2), 2.0) + TMath::Power(kaSys->GetBinError(2), 2.0))*0.5;
  } else {
    kaYmerged = kaStat->GetBinContent(3+phiCent);    
    kaYmergedErr = TMath::Sqrt(TMath::Power(kaStat->GetBinError(3+phiCent), 2.0) + TMath::Power(kaSys->GetBinError(3+phiCent), 2.0));  
  }
  
  //Get Phi values
  Printf("::::: Centrality PbPb: %i, kaon = %f +/- %f", phiCent, kaYmerged, kaYmergedErr);
  phiY = phiStat->GetY()[phiCent];
  phiYerr = TMath::Sqrt(TMath::Power(phiStat->GetErrorY(phiCent), 2.0) + TMath::Power(phiSys->GetErrorY(phiCent), 2.0));
  Printf("::::: Centrality PbPb : %i, phi = %f +/- %f", phiCent, phiY, phiYerr);

  if (!x || !y) return;
  x = phiY / kaYmerged;
  y = TMath::Sqrt(TMath::Power(phiYerr/phiY, 2.0) + TMath::Power(kaYmergedErr/kaYmerged, 2.0)) * x;
  Printf("::::: Centrality PbPb: %i, ratio = %f +/- %f", phiCent, x, y);
  return;
}

void GetRatioSysStatSq_XeXe(Double_t &x, Double_t &y, Int_t phiCent, TH1F * kaStat, TH1F * kaSys, TGraphErrors * phiStat, TGraphErrors * phiSys)
{
  if (!kaStat  || !kaSys || !phiStat || !phiSys) return 0.0;
  Printf(":::: PhiCent = %i", phiCent);
  Double_t kaYmerged = 0.0;
  Double_t kaYmergedErr = 0.0;
  Double_t phiY = 0.0;
  Double_t phiYerr = 0.0;

  //centrality sums for kaon yields
  if (phiCent == 0){
    
    kaYmerged = kaStat->GetBinContent(1)*5.0 + kaStat->GetBinContent(2)*5.0 + kaStat->GetBinContent(3)*10.0 + kaStat->GetBinContent(4)*10.0;
    kaYmerged /= 30.0;
    kaYmergedErr  = TMath::Sqrt(TMath::Power(kaStat->GetBinError(1), 2.0) + TMath::Power(kaSys->GetBinError(1), 2.0))*5.0;
    kaYmergedErr += TMath::Sqrt(TMath::Power(kaStat->GetBinError(2), 2.0) + TMath::Power(kaSys->GetBinError(2), 2.0))*5.0;
    kaYmergedErr += TMath::Sqrt(TMath::Power(kaStat->GetBinError(3), 2.0) + TMath::Power(kaSys->GetBinError(3), 2.0))*5.0;
    kaYmergedErr += TMath::Sqrt(TMath::Power(kaStat->GetBinError(4), 2.0) + TMath::Power(kaSys->GetBinError(4), 2.0))*5.0;
    kaYmergedErr/=30.0;
    
  } else if (phiCent == 1) {
    kaYmerged = kaStat->GetBinContent(5) + kaStat->GetBinContent(6) + kaStat->GetBinContent(7);
    kaYmerged/= 3.0;
    kaYmergedErr  = TMath::Sqrt(TMath::Power(kaStat->GetBinError(5), 2.0) + TMath::Power(kaSys->GetBinError(5), 2.0))*5.0;
    kaYmergedErr += TMath::Sqrt(TMath::Power(kaStat->GetBinError(6), 2.0) + TMath::Power(kaSys->GetBinError(6), 2.0))*5.0;
    kaYmergedErr += TMath::Sqrt(TMath::Power(kaStat->GetBinError(7), 2.0) + TMath::Power(kaSys->GetBinError(7), 2.0))*5.0;
    kaYmergedErr/=30.0;
  } else if (phiCent == 2) {
    kaYmerged = kaStat->GetBinContent(8) + kaStat->GetBinContent(9) + kaStat->GetBinContent(10);
    kaYmerged/= 3.0;
    kaYmergedErr  = TMath::Sqrt(TMath::Power(kaStat->GetBinError(8), 2.0) + TMath::Power(kaSys->GetBinError(8), 2.0))*5.0;
    kaYmergedErr += TMath::Sqrt(TMath::Power(kaStat->GetBinError(9), 2.0) + TMath::Power(kaSys->GetBinError(9), 2.0))*5.0;
    kaYmergedErr += TMath::Sqrt(TMath::Power(kaStat->GetBinError(10), 2.0) + TMath::Power(kaSys->GetBinError(10), 2.0))*5.0;
    kaYmergedErr/=30.0;
  } else {
    kaYmerged = 0.0;
    kaYmergedErr = 0.0;
  }
  
  //Get Phi values
  Printf("::::: Centrality XeXe: %i, kaon = %f +/- %f", phiCent, kaYmerged, kaYmergedErr);
  phiY = phiStat->GetY()[phiCent];
  phiYerr = TMath::Sqrt(TMath::Power(phiStat->GetErrorY(phiCent), 2.0) + TMath::Power(phiSys->GetErrorY(phiCent), 2.0));
  Printf("::::: Centrality XeXe: %i, phi = %f +/- %f", phiCent, phiY, phiYerr);

  if (!x || !y) return;
  x = phiY / kaYmerged;
  y = TMath::Sqrt(TMath::Power(phiYerr/phiY, 2.0) + TMath::Power(kaYmergedErr/kaYmerged, 2.0)) * x;
  Printf("::::: Centrality XeXe: %i, ratio = %f +/- %f", phiCent, x, y);
  return;
}

void GetRatio_XeXe(Double_t &x, Double_t &y, Int_t phiCent, TH1F * kaStat, TGraphErrors * phiStat)
{
  if (!kaStat || !phiStat) return 0.0;
  Printf(":::: PhiCent = %i", phiCent);
  Double_t kaYmerged = 0.0;
  Double_t kaYmergedErr = 0.0;
  Double_t phiY = 0.0;
  Double_t phiYerr = 0.0;

  //centrality sums for kaon yields
  if (phiCent == 0){
    kaYmerged    = (kaStat->GetBinContent(1) + kaStat->GetBinContent(2)) * 0.5;
    kaYmergedErr = (kaStat->GetBinError(1) + kaStat->GetBinError(2)) * 0.5;
  } else if (phiCent == 1) {
    kaYmerged    = (kaStat->GetBinContent(3) + kaStat->GetBinContent(4)) * 0.5;
    kaYmergedErr = (kaStat->GetBinError(3) + kaStat->GetBinError(4)) * 0.5;
  } else if (phiCent == 2) {
    kaYmerged    = (kaStat->GetBinContent(5) + kaStat->GetBinContent(6) + kaStat->GetBinContent(7)) / 3.0;
    kaYmergedErr = (kaStat->GetBinError(5) + kaStat->GetBinError(6) + kaStat->GetBinError(7)) / 3.0;
  } else if (phiCent == 3) {
    kaYmerged    = (kaStat->GetBinContent(8)*0.5 + kaStat->GetBinContent(9)) / 1.5;
    kaYmergedErr = (kaStat->GetBinError(8)*0.5 + kaStat->GetBinError(9)) / 1.5;
  }
  
  //Get Phi values
  Printf("::::: Centrality XeXe: %i, kaon = %f +/- %f", phiCent, kaYmerged, kaYmergedErr);
  phiY = phiStat->GetY()[phiCent];
  phiYerr = phiStat->GetEY()[phiCent]; 
  Printf("::::: Centrality XeXe: %i, phi = %f +/- %f", phiCent, phiY, phiYerr);

  if (!x || !y) return;
  x = 2* phiY / kaYmerged;
  y = 2* TMath::Sqrt(TMath::Power(phiYerr/phiY, 2.0) + TMath::Power(kaYmergedErr/kaYmerged, 2.0)) * x;
  Printf("::::: Centrality XeXe: %i, ratio = %f +/- %f", phiCent, x, y);
  return;
}

void GetRatio_PbPb(Double_t &x, Double_t &y, Int_t phiCent, TH1F * kaStat, TGraphErrors * phiStat)
{
  if (!kaStat  || !phiStat) return 0.0;
  Printf(":::: PhiCent = %i", phiCent);
  Double_t kaYmerged = 0.0;
  Double_t kaYmergedErr = 0.0;
  Double_t phiY = 0.0;
  Double_t phiYerr = 0.0;

  //centrality sums for kaon yields
  if (phiCent == 0){ //phi 0-10%, kaons 0-5%, 5-10%
    kaYmerged = kaStat->GetBinContent(1)*0.5 + kaStat->GetBinContent(2)*0.5;
    kaYmergedErr = kaStat->GetBinError(1)*0.5 + kaStat->GetBinError(2)*0.5;
  } else {
    kaYmerged = kaStat->GetBinContent(2+phiCent);    
    kaYmergedErr = kaStat->GetBinError(2+phiCent);  
  }
  
  //Get Phi values
  Printf("::::: Centrality PbPb: %i, kaon = %f +/- %f", phiCent, kaYmerged, kaYmergedErr);
  phiY = phiStat->GetY()[phiCent];
  phiYerr = phiStat->GetEY()[phiCent]; 
  Printf("::::: Centrality PbPb : %i, phi = %f +/- %f", phiCent, phiY, phiYerr);

  if (!x || !y) return;
  x = 2* phiY / kaYmerged;
  y = 2* TMath::Sqrt(TMath::Power(phiYerr/phiY, 2.0) + TMath::Power(kaYmergedErr/kaYmerged, 2.0)) * x;
  Printf("::::: Centrality PbPb: %i, ratio = %f +/- %f", phiCent, x, y);
  return;

}

TGraphAsymmErrors *  GetRatio_PbPb276(Bool_t sys)
{
  // PRC 91, 024609 (2015) (DOI:10.1103/PhysRevC.91.024609) 
  //http://hepdata.cedar.ac.uk/view/ins1288320/d36
  TGraphAsymmErrors * p8573_d36x1y1 = 0x0;
  int p8573_d36x1y1_numpoints = 10;
  
  double p8573_d36x1y1_xval[] = {1601., 1294., 966., 649., 426., 261., 149., 76., 35., 17.52};
  double p8573_d36x1y1_xerrminus[] = {60., 49., 37., 23., 15., 9., 6., 4., 2., 1.84};
  double p8573_d36x1y1_xerrplus[] = {60., 49., 37., 23., 15., 9., 6., 4., 2., 1.84};
  
  double p8573_d36x1y1_yval[] = { 0.009383, 0.009676, 0.009957, 0.01136, 0.01069, 0.0108, 0.01046, 0.009676, 0.008851, 0.007349 };
  double p8573_d36x1y1_yerrminus[] = { 9.199739126736149E-4, 8.787104187387334E-4, 7.671779454598523E-4, 8.640023148117139E-4, 7.21110255092798E-4, 6.992138442565337E-4, 6.767569726275453E-4, 6.122826144845205E-4, 7.77618158224202E-4, 6.785226599016425E-4 };
  double p8573_d36x1y1_yerrplus[] = { 9.199739126736149E-4, 8.787104187387334E-4, 7.671779454598523E-4, 8.640023148117139E-4, 7.21110255092798E-4, 6.992138442565337E-4, 6.767569726275453E-4, 6.122826144845205E-4, 7.77618158224202E-4, 6.785226599016425E-4 };
  double p8573_d36x1y1_ystatminus[] = { 3.16E-4, 3.26E-4, 2.39E-4, 2.4E-4, 2.4E-4, 2.0E-4, 2.2E-4, 2.67E-4, 3.01E-4,  2.93E-4 };
  double p8573_d36x1y1_ystatplus[] = { 3.16E-4, 3.26E-4, 2.39E-4, 2.4E-4, 2.4E-4, 2.0E-4, 2.2E-4, 2.67E-4, 3.01E-4, 2.93E-4 };

  
  if (sys) p8573_d36x1y1 = new TGraphAsymmErrors(p8573_d36x1y1_numpoints, p8573_d36x1y1_xval, p8573_d36x1y1_yval, p8573_d36x1y1_xerrminus, p8573_d36x1y1_xerrplus, p8573_d36x1y1_yerrminus, p8573_d36x1y1_yerrplus);
  else p8573_d36x1y1 = new TGraphAsymmErrors(p8573_d36x1y1_numpoints, p8573_d36x1y1_xval, p8573_d36x1y1_yval, p8573_d36x1y1_xerrminus, p8573_d36x1y1_xerrplus, p8573_d36x1y1_ystatminus, p8573_d36x1y1_ystatplus);
  p8573_d36x1y1->SetName("/HepData/8573/d36x1y1");
  p8573_d36x1y1->SetTitle("/HepData/8573/d36x1y1");
  return p8573_d36x1y1;
  
}

TGraphAsymmErrors *  GetRatio_pPb502(Bool_t sys)
{
  //ALICE 	Eur. Phys. J. C 76 (2016) 245
  //http://hepdata.cedar.ac.uk/view/ins1418181/d30
  //2phi/ka
  double p9065_d30x1y1_xval[] = { 45.0, 36.2, 30.5, 23.2, 16.1, 9.8, 4.16 };
  double p9065_d30x1y1_xerrminus[] = { 1.0, 0.8, 0.7, 0.5, 0.4, 0.2, 0.09 };
  double p9065_d30x1y1_xerrplus[] = { 1.0, 0.8, 0.7, 0.5, 0.4, 0.2, 0.09 };
  
  double p9065_d30x1y1_yval[] = { 0.0185, 0.01743, 0.01741, 0.01704, 0.01612, 0.01474, 0.01426 };
  double p9065_d30x1y1_yerrminus[] = { 0.0016835082417380675, 0.0014335968749965941, 0.0013756816492197603, 0.0013773525329413671, 0.0013773525329413671, 0.0012300812981262661, 0.0016333401360402555 };
  double p9065_d30x1y1_yerrplus[] = { 0.0016835082417380675, 0.0014335968749965941, 0.0013756816492197603, 0.0013773525329413671, 0.0013773525329413671, 0.0012300812981262661, 0.0016333401360402555 };
  double p9065_d30x1y1_ystatminus[] = { 1.8E-4, 1.8E-4, 1.4E-4, 1.1E-4, 1.1E-4, 1.5E-4, 2.1E-4 };
  double p9065_d30x1y1_ystatplus[] = { 1.8E-4, 1.8E-4, 1.4E-4, 1.1E-4, 1.1E-4, 1.5E-4, 2.1E-4 };

  int p9065_d30x1y1_numpoints = 7;

  //divide by 2 to go to phi/ka
  //Eur.Phys.J. C76 (2016) no.5, 245 
  for (int i = 0; i<7; i++){
    p9065_d30x1y1_yval[i]/=2.0;
    p9065_d30x1y1_yerrminus[i]/=2.0;
    p9065_d30x1y1_yerrplus[i]/=2.0;
    p9065_d30x1y1_ystatminus[i]/=2.0;
    p9065_d30x1y1_ystatplus[i]/=2.0;
  }
  
  TGraphAsymmErrors * graph = 0;
  if (sys) graph = new TGraphAsymmErrors(p9065_d30x1y1_numpoints, p9065_d30x1y1_xval, p9065_d30x1y1_yval, p9065_d30x1y1_xerrminus, p9065_d30x1y1_xerrplus, p9065_d30x1y1_yerrminus, p9065_d30x1y1_yerrplus);
  else graph = new TGraphAsymmErrors(p9065_d30x1y1_numpoints, p9065_d30x1y1_xval, p9065_d30x1y1_yval, p9065_d30x1y1_xerrminus, p9065_d30x1y1_xerrplus, p9065_d30x1y1_ystatminus, p9065_d30x1y1_ystatplus);
  graph->SetName("/HepData/9065/d30x1y1");
  graph->SetTitle("/HepData/9065/d30x1y1");

  return graph;
}

TGraphErrors *  GetRatio_pp13(Int_t sys)
{
  //Preliminary QM 2018
  // phi AN
  // kaon AN
  TFile * fin = TFile::Open("./SummaryPlotsFromAjay/IdentifiedParticlesRatioPlots_AllPart.root");
  TGraphErrors * graph_stat = (TGraphErrors*) fin->Get("gr[5][3][0]");
  TGraphErrors * graph_syst = (TGraphErrors*) fin->Get("gr[5][3][1]");
  TGraphErrors * graph_uncor = (TGraphErrors*) fin->Get("gr[5][3][2]");

  graph_stat->SetName("prel_phi2ka_pp13TeV_stat");
  graph_syst->SetName("prel_phi2ka_pp13TeV_syst");
  graph_uncor->SetName("prel_phi2ka_pp13TeV_uncor");

  if (sys) return graph_syst;
  return graph_stat;
}

TGraphErrors * GetPhi2Ka_XeXe544(Bool_t sys)
{
  TString KaonXeFileName = "Preliminary_XeXe_YieldAndMeanPtSumKwSys.root";  
  TString XeFileName = ReadInputFromFile("phiXeXe_dNdy_preliminary.dat", "graph", 1, kRed, "XeXe", 5.44);
  
  //Get phi 
  TFile * XeFile = TFile::Open(XeFileName.Data());
  if (!XeFile) return 0;
  TGraphErrors *gPhi_dNdy_XeXe5TeV_stat = (TGraphErrors*) XeFile->Get("stat");
  TGraphErrors *gPhi_dNdy_XeXe5TeV_syst = (TGraphErrors*) XeFile->Get("sys");

  //Get ka+ + ka-
  TFile * KaonXeFile = TFile::Open(KaonXeFileName.Data());
  if (!KaonXeFile) return 0;
  TH1F *gKa_dNdy_XeXe5TeV_stat = (TH1F*) KaonXeFile->Get("hStdYieldSummedKaon");
  TH1F *gKa_dNdy_XeXe5TeV_syst = (TH1F*) KaonXeFile->Get("hStdYieldSysSummedKaon");
  
  TGraphErrors * gPhi2KaXeXe = new TGraphErrors(4);
  TGraphErrors * gPhi2KaXeXe_syst = new TGraphErrors(4);
  
  for (int cc=0; cc<4; cc++){
    Double_t y, ey;
    //stat
    GetRatio_XeXe(y, ey, cc, gKa_dNdy_XeXe5TeV_stat, gPhi_dNdy_XeXe5TeV_stat);
    gPhi2KaXeXe->SetPoint(cc, gPhi_dNdy_XeXe5TeV_stat->GetX()[cc], y);
    gPhi2KaXeXe->SetPointError(cc, gPhi_dNdy_XeXe5TeV_stat->GetEX()[cc], ey);
    //syst
    GetRatio_XeXe(y, ey, cc, gKa_dNdy_XeXe5TeV_syst, gPhi_dNdy_XeXe5TeV_syst);
    gPhi2KaXeXe_syst->SetPoint(cc, gPhi_dNdy_XeXe5TeV_syst->GetX()[cc], y);
    gPhi2KaXeXe_syst->SetPointError(cc, gPhi_dNdy_XeXe5TeV_syst->GetEX()[cc], ey);
  }
  
  if (sys) return gPhi2KaXeXe_syst;
  return gPhi2KaXeXe;
}


TGraphErrors * GetPhi2Ka_PbPb502(Bool_t sys)
{
  TString KaonPbFileName = "Preliminary_YieldAndMeanPtSumKawSys_PbPb5TeV.root";  
  TString PbFileName = "Preliminary_PhidNdYMeanPtV0M_PbPb5TeV.root";

  //Get phi 
  TFile * PbFile = TFile::Open(PbFileName.Data());
  if (!PbFile) return 0; 
  TGraphErrors *gPhi_dNdy_PbPb5TeV_stat = (TGraphErrors*) PbFile->Get("gYield_stat_phi");
  TGraphErrors *gPhi_dNdy_PbPb5TeV_syst = (TGraphErrors*) PbFile->Get("gYield_syst_phi");

  //Get ka+ + ka-
  TFile * KaonPbFile = TFile::Open(KaonPbFileName.Data());
  if (!KaonPbFile) return 0;
  TH1F *gKa_dNdy_PbPb5TeV_stat = (TH1F*) KaonPbFile->Get("hStdYieldSummedKaon");
  TH1F *gKa_dNdy_PbPb5TeV_syst = (TH1F*) KaonPbFile->Get("hStdYieldSysSummedKaon");

  TGraphErrors * gPhi2KaPbPb = new TGraphErrors(8);
  TGraphErrors * gPhi2KaPbPb_syst = new TGraphErrors(8);

  for (int cc=0; cc<8; cc++){
    Double_t y, ey;
    //stat
    GetRatio_PbPb(y, ey, cc, gKa_dNdy_PbPb5TeV_stat, gPhi_dNdy_PbPb5TeV_stat);
    gPhi2KaPbPb->SetPoint(cc, gPhi_dNdy_PbPb5TeV_stat->GetX()[cc], y);
    gPhi2KaPbPb->SetPointError(cc, gPhi_dNdy_PbPb5TeV_stat->GetEX()[cc], ey);
    //syst
    GetRatio_PbPb(y, ey, cc, gKa_dNdy_PbPb5TeV_syst, gPhi_dNdy_PbPb5TeV_syst);
    gPhi2KaPbPb_syst->SetPoint(cc, gPhi_dNdy_PbPb5TeV_stat->GetX()[cc], y);
    gPhi2KaPbPb_syst->SetPointError(cc, gPhi_dNdy_PbPb5TeV_stat->GetEX()[cc], ey);
  }

  if (sys) return gPhi2KaPbPb_syst;
  return gPhi2KaPbPb;
}


