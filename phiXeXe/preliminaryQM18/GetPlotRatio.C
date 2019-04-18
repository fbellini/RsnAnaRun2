#include "/Users/fbellini/alice/macros/SetStyle.C"
#include "/Users/fbellini/alice/macros/Beautify.C"
TH1F * GetPlotRatio(TH1F* num = 0x0,
		    TH1F* den = 0x0);

TH1D * GetPlotRatio(TH1D* num = 0x0,
		    TH1D* den = 0x0,
		    Bool_t display = 1,
		    TString savepngName = "LSB2MEB.eps",
		    TString numName = "LSB",
		    TString denName = "MEB",
		    TString titleY = "d^{2}#it{N}_{raw}/(d#it{y}d#it{p}_{T})",
		    Int_t rebinFactorNum = 0,
		    Int_t rebinFactorDen = 0,
		    Float_t showMin=0.0,
		    Float_t showMax=-1.e10,
		    TString errorType="sum2")
{

  /* returns ratio of two plots passed as arguments and displays if requested */
  if (!num) {
    Printf("invalid numerator passed as argument");
    return 0;
  }
  if (!den) {
    Printf("invalid numerator passed as argument");
    return 0;
  }

  //  SetStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleOffset(0.8,"y");

  Beautify(num, num->GetLineColor(), 1, 2, 20, 1.0);
  Beautify(den, num->GetLineColor(), 1, 2, 25, 1.0);
 
  TCanvas * cr = 0x0;
  TPad * pad1 = 0x0;
  TPad * pad2 = 0x0;
  if (display){
    cr = new TCanvas("cr","compare",800, 800);
    cr->cd();
    pad1 = new TPad("pad1","This is pad1",0.001,0.5,0.999,0.999);
    pad1->SetFillColor(0);
    pad1->SetBorderMode(0);
    pad1->SetBorderSize(0);
    pad1->SetMargin(0.15,0.05,0.01,0.02);

    pad2 = new TPad("pad2","This is pad2",0.001,0.001, 0.999,0.5);
    pad2->SetFillColor(0);
    pad2->SetBorderMode(0);
    pad2->SetBorderSize(0);
    pad2->SetMargin(0.15,0.05,0.2,0.01);
    pad1->Draw();
    pad2->Draw();

    pad1->cd();
    pad1->SetLogy();
    if (numName.CompareTo("",TString::kExact)==0) numName = num->GetName();
    if (numName.CompareTo("",TString::kExact)==0) denName = den->GetName();
    num->SetTitle(numName.Data());

    den->SetTitle(denName.Data());
    num->SetFillStyle(0);
    den->SetFillStyle(0);
    den->GetYaxis()->SetTitle(titleY.Data());
  
    if (showMax>showMin) den->GetXaxis()->SetRangeUser(showMin,showMax);
    den->Draw("");
    num->Draw("same");
    den->GetYaxis()->SetTitle(titleY.Data());

    TLegend * leg = (TLegend*) pad1->BuildLegend(0.7,0.75,0.93,0.93);
    myLegendSetUp(leg, 0.06);
    leg->SetEntryOption("pl");
  }

  //get ratio
  // num->Sumw2();
  // den->Sumw2();
  if (rebinFactorNum>0) num->Rebin(rebinFactorNum);
  if (rebinFactorDen>0) den->Rebin(rebinFactorDen);
  num->SetTitle(numName.Data());
  den->SetTitle(denName.Data());

  TH1D * ratio = (TH1D*) num->Clone("LSB");
  ratio->SetTitle(Form("%s/%s",num->GetTitle(), den->GetTitle()));
  ratio->SetMarkerStyle(20);
  if (errorType.Contains("num")) {
    ratio->Reset();
    for (Int_t j=1; j<num->GetNbinsX()+1;j++){
      if (den->GetBinContent(j)>0) {
	ratio->SetBinContent(j, num->GetBinContent(j)/den->GetBinContent(j));
	ratio->SetBinError(j, num->GetBinError(j)/num->GetBinContent(j));
      } else {
	ratio->SetBinContent(j,0);
      	ratio->SetBinError(j,0);
      }
    }
  } else {
    ratio->Divide(num, den, 1., 1., errorType.Data());
  }
  
  if (display) {
    Float_t xmax = ratio->GetXaxis()->GetBinUpEdge(ratio->GetNbinsX());
    TLine * at1 = new TLine(ratio->GetXaxis()->GetBinLowEdge(1), 1., xmax, 1.);
    at1->SetLineStyle(7);
    at1->SetLineWidth(1);

    TLine * at120 = new TLine(ratio->GetXaxis()->GetBinLowEdge(1), 1.2, xmax, 1.2);
    at120->SetLineStyle(2);
    at120->SetLineWidth(1);
    
    TLine * at080 = new TLine(ratio->GetXaxis()->GetBinLowEdge(1), 0.8, xmax, 0.8);
    at080->SetLineStyle(2);
    at080->SetLineWidth(1);

    ratio->GetYaxis()->SetTitle(Form("%s / %s",num->GetTitle(), den->GetTitle()));
    ratio->GetXaxis()->SetTitle("#it{p}_{T}(GeV/#it{c})");
    if (showMax>showMin) ratio->GetXaxis()->SetRangeUser(showMin,showMax);
    pad2->cd();
    ratio->Draw("");
    TPaveText *pave = new TPaveText(0.6, 0.25, 0.8, 0.32, "NDC");
    if (errorType.Contains("num"))
      pave->InsertText(Form("Relative uncert. of numerator."));
    else
      pave->InsertText(Form("%s uncert.", errorType.Data()));
    pave->SetBorderSize(0);
    pave->SetTextFont(42);
    pave->SetFillStyle(0);
    pave->Draw();
    at1->Draw();
    at120->Draw();
    at080->Draw();
    // TLegend * leg2 = (TLegend*) pad2->BuildLegend(0.6,0.88, 0.88,0.94);
    // leg2->SetEntryOption("p");
    // myLegendSetUp(leg2, 0.06);
  }
  
  if (savepngName.Contains(".png") || savepngName.Contains(".eps")) cr->Print(savepngName.Data());
  return ratio;
}
TH1F * GetPlotRatio(TH1F* num,
		    TH1F* den)
{
  Bool_t display = 0;
  TString savepngName = "";
  TString numName = "";
  TString denName = "";
  TString titleY = "d^{2}#it{N}/(d#it{y}d#it{p}_{T})";
  Int_t rebinFactorNum = 0;
  Int_t rebinFactorDen = 0;
  Float_t showMin=0.0;
  Float_t showMax=-1.e10;
  TString errorType="sum2";
  /* returns ratio of two plots passed as arguments and displays if requested */
  if (!num) {
    Printf("invalid numerator passed as argument");
    return 0;
  }
  if (!den) {
    Printf("invalid numerator passed as argument");
    return 0;
  }

  //  SetStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleOffset(0.8,"y");

  Beautify(num, num->GetLineColor(), 1, 2, 20, 1.0);
  Beautify(den, num->GetLineColor(), 1, 2, 25, 1.0);
 
  TCanvas * cr = 0x0;
  TPad * pad1 = 0x0;
  TPad * pad2 = 0x0;
  if (display){
    cr = new TCanvas("cr","compare",800, 800);
    cr->cd();
    pad1 = new TPad("pad1","This is pad1",0.001,0.5,0.999,0.999);
    pad1->SetFillColor(0);
    pad1->SetBorderMode(0);
    pad1->SetBorderSize(0);
    pad1->SetMargin(0.15,0.05,0.01,0.02);

    pad2 = new TPad("pad2","This is pad2",0.001,0.001, 0.999,0.5);
    pad2->SetFillColor(0);
    pad2->SetBorderMode(0);
    pad2->SetBorderSize(0);
    pad2->SetMargin(0.15,0.05,0.2,0.01);
    pad1->Draw();
    pad2->Draw();

    pad1->cd();
    pad1->SetLogy();
    if (numName.CompareTo("",TString::kExact)==0) numName = num->GetName();
    if (numName.CompareTo("",TString::kExact)==0) denName = den->GetName();
    num->SetTitle(numName.Data());

    den->SetTitle(denName.Data());
    num->SetFillStyle(0);
    den->SetFillStyle(0);
    den->GetYaxis()->SetTitle(titleY.Data());
  
    if (showMax>showMin) den->GetXaxis()->SetRangeUser(showMin,showMax);
    den->Draw("");
    num->Draw("same");
    den->GetYaxis()->SetTitle(titleY.Data());

    TLegend * leg = (TLegend*) pad1->BuildLegend(0.7,0.75,0.93,0.93);
    myLegendSetUp(leg, 0.06);
    leg->SetEntryOption("pl");
  }

  //get ratio
  // num->Sumw2();
  // den->Sumw2();
  if (rebinFactorNum>0) num->Rebin(rebinFactorNum);
  if (rebinFactorDen>0) den->Rebin(rebinFactorDen);
  num->SetTitle(numName.Data());
  den->SetTitle(denName.Data());

  TH1F * ratio = (TH1F*) num->Clone("LSB");
  ratio->SetTitle(Form("%s/%s",num->GetTitle(), den->GetTitle()));
  ratio->SetMarkerStyle(20);
  if (errorType.Contains("num")) {
    ratio->Reset();
    for (Int_t j=1; j<num->GetNbinsX()+1;j++){
      if (den->GetBinContent(j)>0) {
	ratio->SetBinContent(j, num->GetBinContent(j)/den->GetBinContent(j));
	ratio->SetBinError(j, num->GetBinError(j)/num->GetBinContent(j));
      } else {
	ratio->SetBinContent(j,0);
      	ratio->SetBinError(j,0);
      }
    }
  } else {
    ratio->Divide(num, den, 1., 1., errorType.Data());
  }
  
  if (display) {
    Float_t xmax = ratio->GetXaxis()->GetBinUpEdge(ratio->GetNbinsX());
    TLine * at1 = new TLine(ratio->GetXaxis()->GetBinLowEdge(1), 1., xmax, 1.);
    at1->SetLineStyle(7);
    at1->SetLineWidth(1);

    TLine * at120 = new TLine(ratio->GetXaxis()->GetBinLowEdge(1), 1.2, xmax, 1.2);
    at120->SetLineStyle(2);
    at120->SetLineWidth(1);
    
    TLine * at080 = new TLine(ratio->GetXaxis()->GetBinLowEdge(1), 0.8, xmax, 0.8);
    at080->SetLineStyle(2);
    at080->SetLineWidth(1);

    ratio->GetYaxis()->SetTitle(Form("%s / %s",num->GetTitle(), den->GetTitle()));
    ratio->GetXaxis()->SetTitle("#it{p}_{T}(GeV/#it{c})");
    if (showMax>showMin) ratio->GetXaxis()->SetRangeUser(showMin,showMax);
    pad2->cd();
    ratio->Draw("");
    TPaveText *pave = new TPaveText(0.6, 0.25, 0.8, 0.32, "NDC");
    if (errorType.Contains("num"))
      pave->InsertText(Form("Relative uncert. of numerator."));
    else
      pave->InsertText(Form("%s uncert.", errorType.Data()));
    pave->SetBorderSize(0);
    pave->SetTextFont(42);
    pave->SetFillStyle(0);
    pave->Draw();
    at1->Draw();
    at120->Draw();
    at080->Draw();
    // TLegend * leg2 = (TLegend*) pad2->BuildLegend(0.6,0.88, 0.88,0.94);
    // leg2->SetEntryOption("p");
    // myLegendSetUp(leg2, 0.06);
  }
  
  if (savepngName.Contains(".png") || savepngName.Contains(".eps")) cr->Print(savepngName.Data());
  return ratio;
}
