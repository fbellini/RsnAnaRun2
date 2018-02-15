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

void BeautifyLine(TF1* obj = 0x0, Color_t color = kBlack, Int_t Line_Style = 1, Int_t Line_Width = 2);
void BeautifyGraph(TGraphAsymmErrors* obj = 0x0, Color_t color = kBlack, Color_t Fill_Color = kGray+1, Int_t Fill_Style = 1001, Int_t Line_Style = 1, Int_t Line_Width = 2, Int_t Marker_Style = 1, Float_t Marker_Size = 1.3);
void BeautifyGraph(TGraphErrors* obj = 0x0, Color_t color = kBlack, Color_t Fill_Color = kGray+1, Int_t Fill_Style = 1001, Int_t Line_Style = 1, Int_t Line_Width = 2, Int_t Marker_Style = 1, Float_t Marker_Size = 1.3);
void Beautify(TH1D* obj = 0x0, Color_t color = kBlack, Int_t Line_Style = 1, Int_t Line_Width = 2, Int_t Marker_Style = 1, Float_t Marker_Size = 1.3);
void Beautify(TH1F* obj = 0x0, Color_t color = kBlack, Int_t Line_Style = 1, Int_t Line_Width = 2, Int_t Marker_Style = 1, Float_t Marker_Size = 1.3);
void myLegendSetUp(TLegend *currentLegend=0,float currentTextSize=0.07);
void myPadSetUp(TPad *currentPad, float currentLeft = 0.11, float currentTop = 0.04, float currentRight = 0.04, float currentBottom = 0.15);
void myGraphSetUp(TGraphErrors *currentGraph = 0, Float_t currentMarkerSize = 1.4, int currentMarkerStyle = 20, Color_t currentMarkerColor = kBlue, int currentLineStyle = 1, int currentLineColor = kBlue);
void myOptions(Int_t lStat = 0);
void SetLabels(TH1F * h = 0, Float_t labelSize = 0.06, Float_t titleSize = 0.07, Float_t titlexoff = 1.1, Float_t titleyoff = 1.1);
void SetLabels(TH1D * h = 0, Float_t labelSize = 0.06, Float_t titleSize = 0.07, Float_t titlexoff = 1.1, Float_t titleyoff = 1.1);
void SetPadProperties(TPad * pad = 0, Float_t leftMargin = 0.2, Float_t rightMargin = 0.02, Float_t lowMargin = 0.18, Float_t upMargin = 0.005);


//---------------------------------------------
void BeautifyLine(TF1* obj, Color_t color, Int_t Line_Style, Int_t Line_Width)
{
  if (!obj) return;
  obj->SetLineColor(color);
  obj->SetLineStyle(Line_Style);
  obj->SetLineWidth(Line_Width);
  return;
}

//----------------------------------------------
void BeautifyGraph(TGraphAsymmErrors* obj, Color_t color, Color_t Fill_Color, Int_t Fill_Style, Int_t Line_Style, Int_t Line_Width, Int_t Marker_Style, Float_t Marker_Size)
{
  if (!obj) return;
  obj->SetLineColor(color);
  obj->SetFillColor(Fill_Color);
  obj->SetFillStyle(Fill_Style);
  obj->SetFillColorAlpha(Fill_Color, 0.6);
  obj->SetLineStyle(Line_Style);
  obj->SetLineWidth(Line_Width);
  obj->SetMarkerColor(color);
  obj->SetMarkerStyle(Marker_Style);
  obj->SetMarkerSize(Marker_Size);
  return;
}

//----------------------------------------------
void BeautifyGraph(TGraphErrors* obj, Color_t color, Color_t Fill_Color, Int_t Fill_Style, Int_t Line_Style, Int_t Line_Width, Int_t Marker_Style, Float_t Marker_Size)
{
  if (!obj) return;
  obj->SetLineColor(color);
  obj->SetFillColor(Fill_Color);
  obj->SetFillStyle(Fill_Style);
  obj->SetFillColorAlpha(Fill_Color, 0.6);
  obj->SetLineStyle(Line_Style);
  obj->SetLineWidth(Line_Width);
  obj->SetMarkerColor(color);
  obj->SetMarkerStyle(Marker_Style);
  obj->SetMarkerSize(Marker_Size);
  return;
}

//----------------------------------------------
void Beautify(TH1F* obj, Color_t color, Int_t Line_Style, Int_t Line_Width, Int_t Marker_Style, Float_t Marker_Size)
{
  if (!obj) return;
  obj->SetLineColor(color);
  obj->SetLineStyle(Line_Style);
  obj->SetLineWidth(Line_Width);
  obj->SetMarkerColor(color);
  obj->SetMarkerStyle(Marker_Style);
  obj->SetMarkerSize(Marker_Size);
  return;
}
//----------------------------------------------
void Beautify(TH1D* obj, Color_t color, Int_t Line_Style, Int_t Line_Width, Int_t Marker_Style, Float_t Marker_Size)
{
  if (!obj) return;
  obj->SetLineColor(color);
  obj->SetLineStyle(Line_Style);
  obj->SetLineWidth(Line_Width);
  obj->SetMarkerColor(color);
  obj->SetMarkerStyle(Marker_Style);
  obj->SetMarkerSize(Marker_Size);
  return;
}

//----------------------------------------------
void myLegendSetUp(TLegend *currentLegend,float currentTextSize)
{
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
void myPadSetUp(TPad *currentPad, float currentLeft, float currentTop, float currentRight, float currentBottom)
{
  currentPad->SetLeftMargin(currentLeft);
  currentPad->SetTopMargin(currentTop);
  currentPad->SetRightMargin(currentRight);
  currentPad->SetBottomMargin(currentBottom);
  return;
}

//----------------------------------------------
void myGraphSetUp(TGraphErrors *currentGraph, Float_t currentMarkerSize, int currentMarkerStyle, Color_t currentMarkerColor, int currentLineStyle, int currentLineColor)
{
  currentGraph->SetMarkerSize(currentMarkerSize);
  currentGraph->SetMarkerStyle(currentMarkerStyle);
  currentGraph->SetMarkerColor(currentMarkerColor);
  currentGraph->SetLineStyle(currentLineStyle);
  currentGraph->SetLineColor(currentLineColor);
  currentGraph->SetFillStyle(1001);
	currentGraph->SetFillColorAlpha(currentLineColor, 0.4);

  return;
}

//----------------------------------------------
void myOptions(Int_t lStat)
{
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

void SetPadProperties(TPad * pad = 0, Float_t leftMargin = 0.2, Float_t rightMargin = 0.02, Float_t lowMargin = 0.18, Float_t upMargin = 0.005)
{
  pad->SetFillColor(kWhite);
  pad->SetFrameFillColor(kWhite);
  pad->SetBorderSize(0);
  pad->SetLeftMargin(leftMargin);
  pad->SetRightMargin(rightMargin);
  pad->SetBottomMargin(lowMargin);
  pad->SetTopMargin(upMargin);
  return;
}

void SetLabels(TH1F * h, Float_t labelSize, Float_t titleSize, Float_t titlexoff, Float_t titleyoff)
{
  if (!h) return;
  h->GetXaxis()->SetLabelSize(labelSize);
  h->GetYaxis()->SetLabelSize(labelSize);
  h->GetXaxis()->SetTitleSize(titleSize);
  h->GetYaxis()->SetTitleSize(titleSize);
  h->GetYaxis()->SetTitleOffset(titleyoff);
  h->GetXaxis()->SetTitleOffset(titlexoff);
  return;
}

void SetLabels(TH1D * h, Float_t labelSize, Float_t titleSize, Float_t titlexoff, Float_t titleyoff)
{
  if (!h) return;
  h->GetXaxis()->SetLabelSize(labelSize);
  h->GetYaxis()->SetLabelSize(labelSize);
  h->GetXaxis()->SetTitleSize(titleSize);
  h->GetYaxis()->SetTitleSize(titleSize);
  h->GetYaxis()->SetTitleOffset(titleyoff);
  h->GetXaxis()->SetTitleOffset(titlexoff);
  return;
}
