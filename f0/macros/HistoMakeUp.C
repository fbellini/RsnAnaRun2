#include "TAxis.h"
#include "TH1.h"
#ifndef HistoMakeUp_C
#define HistoMakeUp_C

void HistoMakeUp(TH1 * histo, 
		 Color_t color = kBlack,
		 Int_t marker = 20,
		 TString titleX = "", 
		 TString titleY = "")
{
  if (!histo) return;
  TAxis * xaxis = (TAxis*)histo->GetXaxis();
  if (!titleX.IsNull()) xaxis->SetTitle(titleX.Data());
  xaxis->SetTitleOffset(1.1);
  xaxis->SetTitleSize(0.05);
  xaxis->SetLabelSize(0.05);
  TAxis * yaxis = (TAxis*)histo->GetYaxis();
  if (!titleX.IsNull()) yaxis->SetTitle(titleY.Data());
  yaxis->SetTitleOffset(1.2);
  yaxis->SetTitleSize(0.05);
  yaxis->SetLabelSize(0.05);
  histo->SetLineColor(color);
  histo->SetLineWidth(2);
  histo->SetMarkerStyle(marker);
  histo->SetMarkerColor(color);
  histo->SetFillColor(color);
  histo->SetFillStyle(0);
  histo->SetMarkerSize(1.);
  return;
}
#endif
