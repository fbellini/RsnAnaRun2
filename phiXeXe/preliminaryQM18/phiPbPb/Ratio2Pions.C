#include "/Users/fbellini/alice/macros/ReadInputFromFile.C"
#include "/Users/fbellini/alice/macros/SetStyle.C"
#include "/Users/fbellini/alice/macros/Beautify.C"
void GetRatioSysStatSq_XeXe(Double_t &x, Double_t &y, Int_t phiCent, TH1F * piStat, TH1F * piSys, TGraphErrors * phiStat, TGraphErrors * phiSys);
void GetRatioSysStatSq_PbPb(Double_t &x, Double_t &y, Int_t phiCent, TH1F * piStat, TH1F * piSys, TGraphErrors * phiStat, TGraphErrors * phiSys);
void GetRatio_XeXe(Double_t &x, Double_t &y, Int_t phiCent, TH1F * piStat, TGraphErrors * phiStat);
void GetRatio_PbPb(Double_t &x, Double_t &y, Int_t phiCent, TH1F * piStat, TGraphErrors * phiStat);
TGraphAsymmErrors *  GetRatio_PbPb276(Bool_t sys = 0);

void Ratio2Pions()
{
  SetStyle();
 
  TString PionXeFileName = "20180502_YieldAndMeanPtSumPiwSys.root";  
  TString PionPbFileName = "Preliminary_YieldAndMeanPtSumPiwSys_PbPb5TeV.root";  

  TString XeFileName = ReadInputFromFile("phiXeXe_dNdy_preliminary.dat", "graph", 1, kRed, "XeXe", 5.44);
  TString PbFileName = "Preliminary_PhidNdYMeanPtV0M_PbPb5TeV.root";

  //Get phi 
  TFile * XeFile = TFile::Open(XeFileName.Data());
  if (!XeFile) return;
  TGraphErrors *gPhi_dNdy_XeXe5TeV_stat = (TGraphErrors*) XeFile->Get("stat");
  TGraphErrors *gPhi_dNdy_XeXe5TeV_syst = (TGraphErrors*) XeFile->Get("sys");

  TFile * PbFile = TFile::Open(PbFileName.Data());
  if (!PbFile) return;
  TGraphErrors *gPhi_dNdy_PbPb5TeV_stat = (TGraphErrors*) PbFile->Get("gYield_stat_phi");
  TGraphErrors *gPhi_dNdy_PbPb5TeV_syst = (TGraphErrors*) PbFile->Get("gYield_syst_phi");


  //Get pi+ + pi-
  TFile * PionXeFile = TFile::Open(PionXeFileName.Data());
  if (!PionXeFile) return;
  TH1F *gPi_dNdy_XeXe5TeV_stat = (TH1F*) PionXeFile->Get("hStdYieldSummedPion");
  TH1F *gPi_dNdy_XeXe5TeV_syst = (TH1F*) PionXeFile->Get("hStdYieldSysSummedPion");
  
  TFile * PionPbFile = TFile::Open(PionPbFileName.Data());
  if (!PionPbFile) return;
  TH1F *gPi_dNdy_PbPb5TeV_stat = (TH1F*) PionPbFile->Get("hStdYieldSummedPion");
  TH1F *gPi_dNdy_PbPb5TeV_syst = (TH1F*) PionPbFile->Get("hStdYieldSysSummedPion");

  TGraphErrors * gPhi2PiXeXe = new TGraphErrors(4);
  TGraphErrors * gPhi2PiXeXe_syst = new TGraphErrors(4);
  TGraphErrors * gPhi2PiPbPb = new TGraphErrors(8);
  TGraphErrors * gPhi2PiPbPb_syst = new TGraphErrors(8);

  for (int cc=0; cc<8; cc++){
    Double_t y, ey;
    if (cc<4){
      //GetRatioSysStatSq_XeXe(y, ey, cc, gPi_dNdy_XeXe5TeV_stat, gPi_dNdy_XeXe5TeV_syst, gPhi_dNdy_XeXe5TeV_stat, gPhi_dNdy_XeXe5TeV_syst);
      //stat
      GetRatio_XeXe(y, ey, cc, gPi_dNdy_XeXe5TeV_stat, gPhi_dNdy_XeXe5TeV_stat);
      gPhi2PiXeXe->SetPoint(cc+1, gPhi_dNdy_XeXe5TeV_stat->GetX()[cc], y);
      gPhi2PiXeXe->SetPointError(cc+1, gPhi_dNdy_XeXe5TeV_stat->GetEX()[cc], ey);
      //syst
      GetRatio_XeXe(y, ey, cc, gPi_dNdy_XeXe5TeV_syst, gPhi_dNdy_XeXe5TeV_syst);
      gPhi2PiXeXe_syst->SetPoint(cc+1, gPhi_dNdy_XeXe5TeV_syst->GetX()[cc], y);
      gPhi2PiXeXe_syst->SetPointError(cc+1, gPhi_dNdy_XeXe5TeV_syst->GetEX()[cc], ey);
    }
    
    //GetRatioSysStatSq_PbPb(y, ey, cc, gPi_dNdy_PbPb5TeV_stat, gPi_dNdy_PbPb5TeV_syst, gPhi_dNdy_PbPb5TeV_stat, gPhi_dNdy_PbPb5TeV_syst);
    //stat
    GetRatio_PbPb(y, ey, cc, gPi_dNdy_PbPb5TeV_stat, gPhi_dNdy_PbPb5TeV_stat);
    gPhi2PiPbPb->SetPoint(cc+1, gPhi_dNdy_PbPb5TeV_stat->GetX()[cc], y);
    gPhi2PiPbPb->SetPointError(cc+1, gPhi_dNdy_PbPb5TeV_stat->GetEX()[cc], ey);
    //syst
    GetRatio_PbPb(y, ey, cc, gPi_dNdy_PbPb5TeV_syst, gPhi_dNdy_PbPb5TeV_syst);
    gPhi2PiPbPb_syst->SetPoint(cc+1, gPhi_dNdy_PbPb5TeV_stat->GetX()[cc], y);
    gPhi2PiPbPb_syst->SetPointError(cc+1, gPhi_dNdy_PbPb5TeV_stat->GetEX()[cc], ey);
  }

  //Get ratio PbPb 276
  TGraphAsymmErrors *gPhi2Pi_PbPb276TeV_stat = (TGraphAsymmErrors*) GetRatio_PbPb276(0);
  TGraphAsymmErrors *gPhi2Pi_PbPb276TeV_syst = (TGraphAsymmErrors*) GetRatio_PbPb276(1);

  BeautifyGraph(gPhi2PiXeXe, kRed, kRed, 0, 1, 1, 20, 1.3); 
  BeautifyGraph(gPhi2PiXeXe_syst, kRed, kRed, 0, 1, 1, 20, 1.3); 
  BeautifyGraph(gPhi2PiPbPb, kBlue+1, kBlue+1, 0, 1, 1, 21, 1.3); 
  BeautifyGraph(gPhi2PiPbPb_syst, kBlue+1, kBlue+1, 0, 1, 1, 21, 1.3); 
  BeautifyGraphAsymmErrors(gPhi2Pi_PbPb276TeV_stat, kGray+1, kGray+1, 0, 1, 1, 24, 1.3); 
  BeautifyGraphAsymmErrors(gPhi2Pi_PbPb276TeV_syst, kGray+1, kGray+1, 0, 1, 1, 24, 1.3); 
  
  TH1F* hf = new TH1F("frame","", 300, 0.0, 3.E3);
  hf->SetMinimum(1.e-3);
  hf->SetMaximum(4.e-2);
  hf->SetXTitle("#LTd#it{N}_{ch}/d#eta#GT");
  hf->SetLineColor(1);
  hf->GetXaxis()->SetTitleSize(0.06);
  hf->GetXaxis()->SetLabelSize(0.045);
  hf->GetXaxis()->SetTitleOffset(1.);
  hf->SetNdivisions(509,"x");
  hf->SetYTitle("#phi/(#pi^{+}+#pi^{-})");
  hf->GetYaxis()->SetTitleOffset(1.2);
  hf->GetYaxis()->SetTitleSize(0.06);
  hf->GetYaxis()->SetLabelSize(0.045);
  hf->GetYaxis()->SetTitleOffset(1.5);
  hf->SetNdivisions(507,"y");
  hf->GetYaxis()->SetRangeUser(3.e-3, 0.018);
  hf->GetXaxis()->SetRangeUser(10, 3000.0);

  gPhi2PiXeXe->SetName("gPhi2Pi_XeXe544_stat");
  gPhi2PiXeXe_syst->SetName("gPhi2Pi_XeXe544_syst");
  gPhi2PiXeXe->SetTitle("Xe-Xe #sqrt{#it{s}_{NN}} = 5.44 TeV");
  gPhi2PiXeXe_syst->SetTitle("Xe-Xe #sqrt{#it{s}_{NN}} = 5.44 TeV");
  gPhi2PiPbPb->SetName("gPhi2Pi_PbPb502_stat");
  gPhi2PiPbPb_syst->SetName("gPhi2Pi_PbPb502_syst");
  gPhi2PiPbPb->SetTitle("Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  gPhi2PiPbPb_syst->SetTitle("Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  gPhi2Pi_PbPb276TeV_stat->SetName("gPhi2Pi_PbPb276_stat");
  gPhi2Pi_PbPb276TeV_syst->SetName("gPhi2Pi_PbPb276_syst");
  gPhi2Pi_PbPb276TeV_stat->SetTitle("Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV");
  gPhi2Pi_PbPb276TeV_syst->SetTitle("Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV");
  
  TLegend * leg = new TLegend(0.22, 0.65, 0.6, 0.9, "#bf{ALICE Preliminary}");
  myLegendSetUp(leg, 0.04);
  leg->AddEntry(gPhi2PiXeXe, "Xe-Xe #sqrt{#it{s}_{NN}} = 5.44 TeV", "p");
  leg->AddEntry(gPhi2PiPbPb, "Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV", "p");
  leg->AddEntry(gPhi2Pi_PbPb276TeV_stat, "Pb-Pb #sqrt{#it{s}_{NN}} = 2.76 TeV", "p");
  
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetPadRightMargin(0.05);
  TCanvas * c1 = new TCanvas("c1", "c1", 600, 600);
  c1->cd();
  gPad->SetLogx();
  hf->Draw();
  gPhi2Pi_PbPb276TeV_syst->Draw("e2 same");
  gPhi2Pi_PbPb276TeV_stat->Draw("pz same");
  gPhi2PiPbPb_syst->Draw("e2 same");
  gPhi2PiPbPb->Draw("pz same");
  gPhi2PiXeXe_syst->Draw("e2 same");
  gPhi2PiXeXe->Draw("pz same");
  leg->Draw();

  TFile * fout = new TFile("PhiToPionRatio.root");
  fout->cd();
  gPhi2Pi_PbPb276TeV_syst->Write();
  gPhi2Pi_PbPb276TeV_stat->Write();
  gPhi2PiPbPb_syst->Write();
  gPhi2PiPbPb->Write();
  gPhi2PiXeXe_syst->Write();
  gPhi2PiXeXe->Write();
  c1->Write();
  c1->Print("comparePhi2Pi_PbXe.pdf");
  c1->Print("comparePhi2Pi_PbXe.eps");
  fout->Close();
  return;
}

void GetRatioSysStatSq_PbPb(Double_t &x, Double_t &y, Int_t phiCent, TH1F * piStat, TH1F * piSys, TGraphErrors * phiStat, TGraphErrors * phiSys)
{
  if (!piStat  || !piSys || !phiStat || !phiSys) return 0.0;
  Printf(":::: PhiCent = %i", phiCent);
  Double_t piYmerged = 0.0;
  Double_t piYmergedErr = 0.0;
  Double_t phiY = 0.0;
  Double_t phiYerr = 0.0;

  //centrality sums for pion yields
  if (phiCent == 0){ //phi 0-10%, pions 0-5%, 5-10%
    piYmerged = piStat->GetBinContent(1)*0.5 + piStat->GetBinContent(2)*0.5;
    piYmergedErr  = TMath::Sqrt(TMath::Power(piStat->GetBinError(1), 2.0) + TMath::Power(piSys->GetBinError(1), 2.0))*0.5;
    piYmergedErr  += TMath::Sqrt(TMath::Power(piStat->GetBinError(2), 2.0) + TMath::Power(piSys->GetBinError(2), 2.0))*0.5;
  } else {
    piYmerged = piStat->GetBinContent(3+phiCent);    
    piYmergedErr = TMath::Sqrt(TMath::Power(piStat->GetBinError(3+phiCent), 2.0) + TMath::Power(piSys->GetBinError(3+phiCent), 2.0));  
  }
  
  //Get Phi values
  Printf("::::: Centrality PbPb: %i, pion = %f +/- %f", phiCent, piYmerged, piYmergedErr);
  phiY = phiStat->GetY()[phiCent];
  phiYerr = TMath::Sqrt(TMath::Power(phiStat->GetErrorY(phiCent), 2.0) + TMath::Power(phiSys->GetErrorY(phiCent), 2.0));
  Printf("::::: Centrality PbPb : %i, phi = %f +/- %f", phiCent, phiY, phiYerr);

  if (!x || !y) return;
  x = phiY / piYmerged;
  y = TMath::Sqrt(TMath::Power(phiYerr/phiY, 2.0) + TMath::Power(piYmergedErr/piYmerged, 2.0)) * x;
  Printf("::::: Centrality PbPb: %i, ratio = %f +/- %f", phiCent, x, y);
  return;
}

void GetRatioSysStatSq_XeXe(Double_t &x, Double_t &y, Int_t phiCent, TH1F * piStat, TH1F * piSys, TGraphErrors * phiStat, TGraphErrors * phiSys)
{
  if (!piStat  || !piSys || !phiStat || !phiSys) return 0.0;
  Printf(":::: PhiCent = %i", phiCent);
  Double_t piYmerged = 0.0;
  Double_t piYmergedErr = 0.0;
  Double_t phiY = 0.0;
  Double_t phiYerr = 0.0;

  //centrality sums for pion yields
  if (phiCent == 0){
    
    piYmerged = piStat->GetBinContent(1)*5.0 + piStat->GetBinContent(2)*5.0 + piStat->GetBinContent(3)*10.0 + piStat->GetBinContent(4)*10.0;
    piYmerged /= 30.0;
    piYmergedErr  = TMath::Sqrt(TMath::Power(piStat->GetBinError(1), 2.0) + TMath::Power(piSys->GetBinError(1), 2.0))*5.0;
    piYmergedErr += TMath::Sqrt(TMath::Power(piStat->GetBinError(2), 2.0) + TMath::Power(piSys->GetBinError(2), 2.0))*5.0;
    piYmergedErr += TMath::Sqrt(TMath::Power(piStat->GetBinError(3), 2.0) + TMath::Power(piSys->GetBinError(3), 2.0))*5.0;
    piYmergedErr += TMath::Sqrt(TMath::Power(piStat->GetBinError(4), 2.0) + TMath::Power(piSys->GetBinError(4), 2.0))*5.0;
    piYmergedErr/=30.0;
    
  } else if (phiCent == 1) {
    piYmerged = piStat->GetBinContent(5) + piStat->GetBinContent(6) + piStat->GetBinContent(7);
    piYmerged/= 3.0;
    piYmergedErr  = TMath::Sqrt(TMath::Power(piStat->GetBinError(5), 2.0) + TMath::Power(piSys->GetBinError(5), 2.0))*5.0;
    piYmergedErr += TMath::Sqrt(TMath::Power(piStat->GetBinError(6), 2.0) + TMath::Power(piSys->GetBinError(6), 2.0))*5.0;
    piYmergedErr += TMath::Sqrt(TMath::Power(piStat->GetBinError(7), 2.0) + TMath::Power(piSys->GetBinError(7), 2.0))*5.0;
    piYmergedErr/=30.0;
  } else if (phiCent == 2) {
    piYmerged = piStat->GetBinContent(8) + piStat->GetBinContent(9) + piStat->GetBinContent(10);
    piYmerged/= 3.0;
    piYmergedErr  = TMath::Sqrt(TMath::Power(piStat->GetBinError(8), 2.0) + TMath::Power(piSys->GetBinError(8), 2.0))*5.0;
    piYmergedErr += TMath::Sqrt(TMath::Power(piStat->GetBinError(9), 2.0) + TMath::Power(piSys->GetBinError(9), 2.0))*5.0;
    piYmergedErr += TMath::Sqrt(TMath::Power(piStat->GetBinError(10), 2.0) + TMath::Power(piSys->GetBinError(10), 2.0))*5.0;
    piYmergedErr/=30.0;
  } else {
    piYmerged = 0.0;
    piYmergedErr = 0.0;
  }
  
  //Get Phi values
  Printf("::::: Centrality XeXe: %i, pion = %f +/- %f", phiCent, piYmerged, piYmergedErr);
  phiY = phiStat->GetY()[phiCent];
  phiYerr = TMath::Sqrt(TMath::Power(phiStat->GetErrorY(phiCent), 2.0) + TMath::Power(phiSys->GetErrorY(phiCent), 2.0));
  Printf("::::: Centrality XeXe: %i, phi = %f +/- %f", phiCent, phiY, phiYerr);

  if (!x || !y) return;
  x = phiY / piYmerged;
  y = TMath::Sqrt(TMath::Power(phiYerr/phiY, 2.0) + TMath::Power(piYmergedErr/piYmerged, 2.0)) * x;
  Printf("::::: Centrality XeXe: %i, ratio = %f +/- %f", phiCent, x, y);
  return;
}

void GetRatio_XeXe(Double_t &x, Double_t &y, Int_t phiCent, TH1F * piStat, TGraphErrors * phiStat)
{
  if (!piStat || !phiStat) return 0.0;
  Printf(":::: PhiCent = %i", phiCent);
  Double_t piYmerged = 0.0;
  Double_t piYmergedErr = 0.0;
  Double_t phiY = 0.0;
  Double_t phiYerr = 0.0;

  //centrality sums for pion yields
  if (phiCent == 0){
    piYmerged    = (piStat->GetBinContent(1) + piStat->GetBinContent(2)) * 0.5;
    piYmergedErr = (piStat->GetBinError(1) + piStat->GetBinError(2)) * 0.5;
  } else if (phiCent == 1) {
    piYmerged    = (piStat->GetBinContent(3) + piStat->GetBinContent(4)) * 0.5;
    piYmergedErr = (piStat->GetBinError(3) + piStat->GetBinError(4)) * 0.5;
  } else if (phiCent == 2) {
    piYmerged    = (piStat->GetBinContent(5) + piStat->GetBinContent(6) + piStat->GetBinContent(7)) / 3.0;
    piYmergedErr = (piStat->GetBinError(5) + piStat->GetBinError(6) + piStat->GetBinError(7)) / 3.0;
  } else if (phiCent == 3) {
    piYmerged    = (piStat->GetBinContent(8) + piStat->GetBinContent(9) + piStat->GetBinContent(10)) / 3.0;
    piYmergedErr = (piStat->GetBinError(8) + piStat->GetBinError(9) + piStat->GetBinError(10)) / 3.0;
  } else {
    piYmerged = 0.0;
    piYmergedErr = 0.0;
  }
  
  //Get Phi values
  //Printf("::::: Centrality XeXe: %i, pion = %f +/- %f", phiCent, piYmerged, piYmergedErr);
  phiY = phiStat->GetY()[phiCent];
  phiYerr = phiStat->GetEY()[phiCent]; 
  //Printf("::::: Centrality XeXe: %i, phi = %f +/- %f", phiCent, phiY, phiYerr);

  if (!x || !y) return;
  x = phiY / piYmerged;
  y = TMath::Sqrt(TMath::Power(phiYerr/phiY, 2.0) + TMath::Power(piYmergedErr/piYmerged, 2.0)) * x;
  //Printf("::::: Centrality XeXe: %i, ratio = %f +/- %f", phiCent, x, y);
  return;
}

void GetRatio_PbPb(Double_t &x, Double_t &y, Int_t phiCent, TH1F * piStat, TGraphErrors * phiStat)
{
  if (!piStat  || !phiStat) return 0.0;
  Printf(":::: PhiCent = %i", phiCent);
  Double_t piYmerged = 0.0;
  Double_t piYmergedErr = 0.0;
  Double_t phiY = 0.0;
  Double_t phiYerr = 0.0;

  //centrality sums for pion yields
  if (phiCent == 0){ //phi 0-10%, pions 0-5%, 5-10%
    piYmerged = piStat->GetBinContent(1)*0.5 + piStat->GetBinContent(2)*0.5;
    piYmergedErr = piStat->GetBinError(1)*0.5 + piStat->GetBinError(2)*0.5;
  } else {
    piYmerged = piStat->GetBinContent(2+phiCent);    
    piYmergedErr = piStat->GetBinError(2+phiCent);  
  }
  
  //Get Phi values
  Printf("::::: Centrality PbPb: %i, pion = %f +/- %f", phiCent, piYmerged, piYmergedErr);
  phiY = phiStat->GetY()[phiCent];
  phiYerr = phiStat->GetEY()[phiCent]; 
  Printf("::::: Centrality PbPb : %i, phi = %f +/- %f", phiCent, phiY, phiYerr);

  if (!x || !y) return;
  x = phiY / piYmerged;
  y = TMath::Sqrt(TMath::Power(phiYerr/phiY, 2.0) + TMath::Power(piYmergedErr/piYmerged, 2.0)) * x;
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
