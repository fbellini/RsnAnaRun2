#if !defined (__CINT__) || (defined(__MAKECINT__))
#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TList.h"
#include "TROOT.h"
#include "TColor.h"
#include "TSystem.h"
#include "TStyle.h"
#endif

//---------------------------------------------
void BeautifyLine(TF1* obj = 0x0,
              Color_t color = kBlack,
              Int_t Line_Style = 1,
              Int_t Line_Width = 2)
{
  if (!obj) return;
  obj->SetLineColor(color);
  obj->SetLineStyle(Line_Style);
  obj->SetLineWidth(Line_Width);
  return;
}

//----------------------------------------------
void BeautifyGraph(TGraphErrors* obj = 0x0,
              Color_t color = kBlack,
              Color_t Fill_Color = kGray+1,
              Int_t Fill_Style = 1001,
              Int_t Line_Style = 1,
              Int_t Line_Width = 2,
              Int_t Marker_Style = 1,
              Float_t Marker_Size = 1.3)
{
  if (!obj) return;
  obj->SetLineColor(color);
  obj->SetFillColor(Fill_Color);
  obj->SetFillStyle(Fill_Style);
  obj->SetFillColorAlpha(Fill_Color, 0.4);
  obj->SetLineStyle(Line_Style);
  obj->SetLineWidth(Line_Width);
  obj->SetMarkerColor(color);
  obj->SetMarkerStyle(Marker_Style);
  obj->SetMarkerSize(Marker_Size);
  return;
}

//----------------------------------------------
void BeautifyGraphAsymmErrors(TGraphAsymmErrors* obj = 0x0,
              Color_t color = kBlack,
              Color_t Fill_Color = kGray+1,
              Int_t Fill_Style = 1001,
              Int_t Line_Style = 1,
              Int_t Line_Width = 2,
              Int_t Marker_Style = 1,
              Float_t Marker_Size = 1.3,
              Float_t Opacity = 0.4)
{
  if (!obj) return;
  obj->SetLineColor(color);
  obj->SetFillStyle(Fill_Style);
  obj->SetFillColorAlpha(Fill_Color, Opacity);
  obj->SetLineStyle(Line_Style);
  obj->SetLineWidth(Line_Width);
  obj->SetMarkerColor(color);
  obj->SetMarkerStyle(Marker_Style);
  obj->SetMarkerSize(Marker_Size);
  return;
}

//---------------------------------------------
void HistoMakeUp(TH1D * obj = 0x0,
		 Color_t color = kBlack,
		 Int_t Marker_Style = 1,
		 TString titleX = "#it{p}_{T} (GeV/#it{c})",
		 TString titleY = "")
{

  if (!obj) return;
  obj->SetLineColor(color);
  obj->SetLineStyle(1);
  obj->SetLineWidth(2);
  obj->SetMarkerColor(color);
  obj->SetMarkerStyle(Marker_Style);
  obj->SetMarkerSize(1.2);
  if (!titleX.IsNull()) obj->GetXaxis()->SetTitle(titleX.Data());
  if (!titleY.IsNull()) obj->GetYaxis()->SetTitle(titleY.Data());
  return;

}
//----------------------------------------------
void Beautify(TH1* obj = 0x0,
              Color_t color = kBlack,
              Int_t Line_Style = 1,
              Int_t Line_Width = 2,
              Int_t Marker_Style = 1,
              Float_t Marker_Size = 1.3)
{
  if (!obj) return;
  obj->SetLineColor(color);
  obj->SetLineStyle(Line_Style);
  obj->SetFillStyle(0);
  obj->SetLineWidth(Line_Width);
  obj->SetMarkerColor(color);
  obj->SetMarkerStyle(Marker_Style);
  obj->SetMarkerSize(Marker_Size);
  return;
}

//----------------------------------------------
void myLegendSetUp(TLegend *currentLegend=0,float currentTextSize=0.07){
  currentLegend->SetTextFont(42);
  currentLegend->SetBorderSize(0);
  currentLegend->SetFillStyle(0);
  currentLegend->SetFillColor(0);
  currentLegend->SetMargin(0.25);
  currentLegend->SetTextSize(currentTextSize);
  currentLegend->SetEntrySeparation(0.5);
  return;
}

//----------------------------------------------
void myPadSetUp(TPad *currentPad, float currentLeft=0.11, float currentTop=0.04, float currentRight=0.04, float currentBottom=0.15){
  currentPad->SetLeftMargin(currentLeft);
  currentPad->SetTopMargin(currentTop);
  currentPad->SetRightMargin(currentRight);
  currentPad->SetBottomMargin(currentBottom);
  return;
}

//----------------------------------------------
void myGraphSetUp(TGraphErrors *currentGraph=0, Float_t currentMarkerSize = 1.4,
		  int currentMarkerStyle=20, Color_t currentMarkerColor=kBlue,
		  int currentLineStyle=1, int currentLineColor=kBlue, Float_t Opacity = 0.4){
  currentGraph->SetMarkerSize(currentMarkerSize);
  currentGraph->SetMarkerStyle(currentMarkerStyle);
  currentGraph->SetMarkerColor(currentMarkerColor);
  currentGraph->SetLineStyle(currentLineStyle);
  currentGraph->SetLineColor(currentLineColor);
  currentGraph->SetFillStyle(1001);
  currentGraph->SetFillColorAlpha(currentLineColor, Opacity);

  return;
}

//----------------------------------------------
void myOptions(Int_t lStat=0){
  // Set gStyle
  int font = 42;
  // From plain
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(10);
  gStyle->SetCanvasColor(10);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleBorderSize(1);
  gStyle->SetStatColor(10);
  gStyle->SetStatBorderSize(1);
  gStyle->SetLegendBorderSize(1);
  //
  gStyle->SetDrawBorder(0);
  gStyle->SetTextFont(font);
  gStyle->SetStatFont(font);
  gStyle->SetStatFontSize(0.05);
  gStyle->SetStatX(0.97);
  gStyle->SetStatY(0.98);
  gStyle->SetStatH(0.03);
  gStyle->SetStatW(0.3);
  gStyle->SetTickLength(0.02,"y");
  gStyle->SetEndErrorSize(3);
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetLabelFont(font,"xyz");
  gStyle->SetLabelOffset(0.01,"xyz");
  gStyle->SetTitleFont(font,"xyz");
  gStyle->SetTitleOffset(1.0,"xyz");
  gStyle->SetLineWidth(2);
  gStyle->SetTitleSize(0.06,"xyz");
  gStyle->SetMarkerSize(1);
  gStyle->SetPalette(1,0);
  if (lStat){
    gStyle->SetOptTitle(1);
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);
    }
  else {
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
  }
}
