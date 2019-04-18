#include "/Users/fbellini/alice/macros/ReadInputFromFile.C"
#include "/Users/fbellini/alice/macros/SetStyle.C"
#include "/Users/fbellini/alice/macros/Beautify.C"

void compareXe2Pb(Bool_t loglog = 0, TString date = "20180430")
{
  SetStyle();
  TString PbFileName = "Preliminary_PhidNdYMeanPtV0M_PbPb5TeV.root";  
  TString XeFileName = ReadInputFromFile("phiXeXe_dNdy_preliminary.dat", "graph", 1, kRed, "XeXe", 5.44);

  TFile * XeFile = TFile::Open(XeFileName.Data());
  if (!XeFile) return;

  TGraphErrors *gPhi_dNdy_XeXe5TeV_stat = (TGraphErrors*) XeFile->Get("stat");
  TGraphErrors *gPhi_dNdy_XeXe5TeV_syst = (TGraphErrors*) XeFile->Get("sys");

  BeautifyGraph(gPhi_dNdy_XeXe5TeV_stat, kRed, kRed, 0, 1, 1, 21, 1.); 
  BeautifyGraph(gPhi_dNdy_XeXe5TeV_syst, kRed, kRed, 0, 1, 1, 21, 1.);
  
  TFile * PbFile = TFile::Open(PbFileName.Data());
  if (!PbFile) return;

  TGraphErrors *gPhi_dNdy_PbPb5TeV_stat = (TGraphErrors*) PbFile->Get("gYield_stat_phi");
  TGraphErrors *gPhi_dNdy_PbPb5TeV_syst = (TGraphErrors*) PbFile->Get("gYield_syst_phi");

  BeautifyGraph(gPhi_dNdy_PbPb5TeV_stat, kGray+1, kGray+1, 0, 1, 1, 20, 1.); 
  BeautifyGraph(gPhi_dNdy_PbPb5TeV_syst, kGray+1, kGray+1, 0, 1, 1, 20, 1.); 

  TH1F* hf = new TH1F("frame","", 2.E3, 0., 2.0E3);
  hf->SetMinimum(0.1);
  hf->SetMaximum(17.);
  hf->SetXTitle("#LTd#it{N}_{ch}/d#eta#GT");
  hf->SetLineColor(1);
  hf->GetXaxis()->SetTitleSize(0.06);
  hf->GetXaxis()->SetLabelSize(0.045);
  hf->GetXaxis()->SetTitleOffset(1.);
  if (loglog) {
    hf->SetNdivisions(510,"x");
    hf->SetNdivisions(509,"y");
  } else {
    hf->SetNdivisions(408,"x");
    hf->SetNdivisions(407,"y");
  }
  hf->SetYTitle("d#it{N}/d#it{y}");
  hf->GetYaxis()->SetTitleOffset(1.2);
  hf->GetYaxis()->SetTitleSize(0.06);
  hf->GetYaxis()->SetLabelSize(0.045);
  hf->GetYaxis()->SetTitleOffset(1.);
  hf->GetYaxis()->SetRangeUser(0.1, 20.0);
  hf->GetXaxis()->SetRangeUser(10, 2000.0);

  TLegend * leg = new TLegend(0.22, 0.7, 0.52, 0.90, "#bf{ALICE Preliminary}");
  myLegendSetUp(leg, 0.05);
  leg->AddEntry(gPhi_dNdy_PbPb5TeV_stat, "Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV", "p");
  leg->AddEntry(gPhi_dNdy_XeXe5TeV_stat, "Xe-Xe #sqrt{#it{s}_{NN}} = 5.44 TeV", "p");
   
  TCanvas * c1 = new TCanvas("c1", "c1", 600, 600);
  c1->cd();
  if (loglog) {
    gPad->SetLogx();
    gPad->SetLogy();
  }
  hf->Draw();
  gPhi_dNdy_PbPb5TeV_syst->Draw("E2 same");
  gPhi_dNdy_PbPb5TeV_stat->Draw("pz same");
  gPhi_dNdy_XeXe5TeV_syst->Draw("E2 same");
  gPhi_dNdy_XeXe5TeV_stat->Draw("pz same");
  leg->Draw();
  c1->Print(Form("compareXe2Pb_phidNdy_%s.pdf", date.Data()));
  c1->Print(Form("compareXe2Pb_phidNdy_%s.eps", date.Data()));
  return;
}

void compareXe2PbMeanPt(Bool_t loglog = 0, TString date = "20180430")
{
  SetStyle();
  TString PbFileName = "Preliminary_PhidNdYMeanPtV0M_PbPb5TeV.root";  
  TString XeFileName = ReadInputFromFile("phiXeXe_meanpt_preliminary.dat", "graph", 1, kRed, "XeXe", 5.44);

  TFile * XeFile = TFile::Open(XeFileName.Data());
  if (!XeFile) return;

  TGraphErrors *gPhi_meanPt_XeXe5TeV_stat = (TGraphErrors*) XeFile->Get("stat");
  TGraphErrors *gPhi_meanPt_XeXe5TeV_syst = (TGraphErrors*) XeFile->Get("sys");

  BeautifyGraph(gPhi_meanPt_XeXe5TeV_stat, kRed, kRed, 0, 1, 1, 21, 1.); 
  BeautifyGraph(gPhi_meanPt_XeXe5TeV_syst, kRed, kRed, 0, 1, 1, 21, 1.);
  
  TFile * PbFile = TFile::Open(PbFileName.Data());
  if (!PbFile) return;

  TGraphErrors *gPhi_meanPt_PbPb5TeV_stat = (TGraphErrors*) PbFile->Get("gMeanPt_stat_phi");
  TGraphErrors *gPhi_meanPt_PbPb5TeV_syst = (TGraphErrors*) PbFile->Get("gMeanPt_syst_phi");

  BeautifyGraph(gPhi_meanPt_PbPb5TeV_stat, kGray+1, kGray+1, 0, 1, 1, 20, 1.3); 
  BeautifyGraph(gPhi_meanPt_PbPb5TeV_syst, kGray+1, kGray+1, 0, 1, 1, 20, 1.3); 

  TH1F* hf = new TH1F("frame","", 2.E3, 0., 2.0E3);
  hf->SetMinimum(0.5);
  hf->SetMaximum(1.5);
  hf->SetXTitle("#LTd#it{N}_{ch}/d#eta#GT");
  hf->SetLineColor(1);
  hf->GetXaxis()->SetTitleSize(0.06);
  hf->GetXaxis()->SetLabelSize(0.045);
  hf->GetXaxis()->SetTitleOffset(1.);
  hf->SetYTitle("#LT#it{p}_{T}#GT (GeV/#it{c})");
  hf->GetYaxis()->SetTitleOffset(1.2);
  hf->GetYaxis()->SetTitleSize(0.06);
  hf->GetYaxis()->SetLabelSize(0.045);
  hf->GetYaxis()->SetTitleOffset(1.);
  hf->GetYaxis()->SetRangeUser(0.97, 1.78);
  hf->GetXaxis()->SetRangeUser(10, 3000.0);
  if (loglog) {
    hf->SetNdivisions(510,"x");
    hf->SetNdivisions(509,"y");
  } else {
    hf->SetNdivisions(408,"x");
    hf->SetNdivisions(306,"y");
  }
  
  TLegend * leg = new TLegend(0.22, 0.75, 0.52, 0.90, "#bf{ALICE Preliminary}");
  myLegendSetUp(leg, 0.05);
  leg->AddEntry(gPhi_meanPt_PbPb5TeV_stat, "Pb-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV", "p");
  leg->AddEntry(gPhi_meanPt_XeXe5TeV_stat, "Xe-Xe, #sqrt{#it{s}_{NN}} = 5.44 TeV", "p");
  
  TCanvas * c1 = new TCanvas("c1", "c1", 600, 600);
  c1->cd();
  //gPad->SetLogx();
  hf->Draw();
  gPhi_meanPt_PbPb5TeV_syst->Draw("E2 same");
  gPhi_meanPt_PbPb5TeV_stat->Draw("pz same");
  gPhi_meanPt_XeXe5TeV_syst->Draw("E2 same");
  gPhi_meanPt_XeXe5TeV_stat->Draw("pz same");
  leg->Draw();
  c1->Print(Form("compareXe2Pb_phiMeanPt_%s.pdf", date.Data()));
  c1->Print(Form("compareXe2Pb_phiMeanPt_%s.eps", date.Data()));

  TFile * foutXe = new TFile("Preliminary_PhiMeanPt_XeXe544TeV.root","recreate");
  foutXe->cd();
  gPhi_meanPt_XeXe5TeV_syst->Write();
  gPhi_meanPt_XeXe5TeV_stat->Write();
  foutXe->Close();
  return;
}


