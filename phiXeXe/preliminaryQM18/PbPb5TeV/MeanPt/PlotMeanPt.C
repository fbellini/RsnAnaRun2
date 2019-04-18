#if !defined (__CINT__) || defined (__CLING__)
#include "/Users/fbellini/alice/macros/SetStyle.C"
#endif


TCanvas *DrawFrame(TString o ="c");
Double_t RatioError(Double_t X,Double_t ErX,Double_t Y,Double_t ErY);
TLine *DrawLine(Double_t xmin = 0,Double_t ymin = 1, Double_t xmax = 1, Double_t ymax = 1,Int_t lStyle = 1, Int_t lColor = 1, Int_t lWidth = 1);
TLatex *DrawText(Double_t x = 0, Double_t y = 0,Int_t tColor = 2,TString name = "label");
TLegend *DrawLegend(Double_t x1,Double_t y1,Double_t x2,Double_t y2);


void PlotMeanPt()
{
  //Pb-Pb 2.76
  //
  gStyle->SetErrorX(0);
  
  //
  //Average dNch/dEta Value
  /*
    0-2.5       2035   52    398    2
    2.5-5       1850   55    372    3
    5-7.5       1666   48    346    4
    7.5-10      1505   44    320    4
    10-20       1180   31    263    4
    20-30        786   20    188    3
    30-40        512   15    131    2
    40-50        318   12     86.3  1.7
    50-60        183    8     53.6 1.2
    60-70         96.3  5.8   30.4 0.8
    70-80         44.9  3.4   15.6 0.5
  */
  /*
    0-5         1942.5   53.5
    5-10        1585.5   46.04
    10-20       1180     31    
    20-30        786     20    
    30-40        512     15    
    40-50        318     12    
    50-60        183      8    
    60-70         96.3    5.8  
    70-80         44.9    3.4
  
  Double_t AvgdNchBydEta[] = {1942.5,1585.5,1180,786,512,318,183,96.3,44.9};
  Double_t AvgdNchBydEtaEr[] = {53.5,46.0,31,20,15,12,8,5.8,3.4};
  */
  TFile *fileks = TFile::Open("FinalKStardNdYMeanPtV0MNew.root");
  TGraphErrors * grKStat = (TGraphErrors*)fileks->Get("gMeanPt_stat_KStar");
  TGraphErrors * grKSyst = (TGraphErrors*)fileks->Get("gMeanPt_syst_KStar");

  TFile *filephi = TFile::Open("FinalPhidNdYMeanPtV0MNew.root");
  TGraphErrors * grPhiStat = (TGraphErrors*)filephi->Get("gMeanPt_stat_phi");
  TGraphErrors * grPhiSyst = (TGraphErrors*)filephi->Get("gMeanPt_syst_phi");
  for(Int_t i = 0; i<grPhiStat->GetN();i++)grPhiStat->SetPointError(i,0,grPhiStat->GetErrorY(i));
  for(Int_t i = 0; i<grKStat->GetN();i++)grKStat->SetPointError(i,0,grKStat->GetErrorY(i));
  
  TFile *f5TeV = TFile::Open("YieldAndMeanPtPbPb5TeV.root");
  f5TeV->ls();
  TGraphErrors  *grPiStat5TeV  = (TGraphErrors*)f5TeV->Get("PionSumByTwoMeanPtStat");
  TGraphErrors  *grPiSyat5TeV    = (TGraphErrors*)f5TeV->Get("PionSumByTwoMeanPtSyst");
  TGraphErrors  *grKaStat5TeV  = (TGraphErrors*)f5TeV->Get("KaonSumByTwoMeanPtStat");
  TGraphErrors  *grKaSyat5TeV    = (TGraphErrors*)f5TeV->Get("KaonSumByTwoMeanPtSyst");
  TGraphErrors  *grPrStat5TeV  = (TGraphErrors*)f5TeV->Get("ProtonSumByTwoMeanPtStat");
  TGraphErrors  *grPrSyat5TeV    = (TGraphErrors*)f5TeV->Get("ProtonSumByTwoMeanPtSyst");

  TCanvas *c1 = new TCanvas("c1","c",10,10,600,600);
  c1->SetBottomMargin(0.15);
  c1->SetLeftMargin(0.16);
  c1->SetTopMargin(0.02);
  c1->SetRightMargin(0.02);
  c1->SetTicks(1,1);
  c1->cd();
  
  grPiStat5TeV->SetMaximum(1.55);
  grPiStat5TeV->SetMinimum(0.01);
  grPiStat5TeV->GetXaxis()->SetLimits(0,2100);
  grPiStat5TeV->GetXaxis()->SetTitle("#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#eta|<0.5}");
  grPiStat5TeV->GetXaxis()->SetTitleSize(0.06);
  grPiStat5TeV->GetXaxis()->SetTitleOffset(1.1);
  grPiStat5TeV->GetXaxis()->SetTitleFont(42);
  //grPiStat5TeV->GetXaxis()->CenterTitle(true);
  grPiStat5TeV->GetXaxis()->SetLabelSize(0.05);
  grPiStat5TeV->GetXaxis()->SetLabelFont(42);
  grPiStat5TeV->GetXaxis()->SetNdivisions(505);
  grPiStat5TeV->GetXaxis()->SetTickLength(0.04);
  grPiStat5TeV->GetYaxis()->SetTitle("#LT#it{p}_{T}#GT (GeV/#it{c})");
  grPiStat5TeV->GetYaxis()->SetTitleFont(42);
  //grPiStat5TeV->GetYaxis()->CenterTitle(true);
  grPiStat5TeV->GetYaxis()->SetTitleSize(0.06);
  grPiStat5TeV->GetYaxis()->SetTitleOffset(1.15);
  grPiStat5TeV->GetYaxis()->SetLabelSize(0.05);
  grPiStat5TeV->GetYaxis()->SetLabelFont(42);
  grPiStat5TeV->GetYaxis()->SetNdivisions(508);
  grPiStat5TeV->GetYaxis()->SetTickLength(0.04);
  
  grPiStat5TeV->SetMarkerColor(2);
  grPiStat5TeV->SetMarkerStyle(33);
  grPiStat5TeV->SetMarkerSize(1.7);
  grPiStat5TeV->SetLineColor(grPiStat5TeV->GetMarkerColor());
  grPiSyat5TeV->SetMarkerStyle(29);
  grPiStat5TeV->Draw("APZ");
  grPiSyat5TeV->SetMarkerColor(grPiStat5TeV->GetMarkerColor());
  grPiSyat5TeV->SetMarkerStyle(grPiStat5TeV->GetMarkerStyle());
  grPiSyat5TeV->SetMarkerSize(grPiStat5TeV->GetMarkerSize());
  grPiSyat5TeV->SetFillStyle(0);
  grPiSyat5TeV->SetLineColor(grPiStat5TeV->GetMarkerColor());
  grPiSyat5TeV->Draw("5P");

  grKaStat5TeV->SetMarkerColor(kGreen+3);
  grKaStat5TeV->SetMarkerStyle(22);
  grKaStat5TeV->SetMarkerSize(1.7);
  grKaStat5TeV->Draw("PZ");

  grKaSyat5TeV->SetMarkerColor(grKaStat5TeV->GetMarkerColor());
  grKaSyat5TeV->SetMarkerStyle(grKaStat5TeV->GetMarkerStyle());
  grKaSyat5TeV->SetMarkerSize(grKaStat5TeV->GetMarkerSize());
  grKaSyat5TeV->SetFillStyle(0);
  grKaSyat5TeV->SetLineColor(grKaStat5TeV->GetMarkerColor());
  grKaSyat5TeV->Draw("5P");

  grPrStat5TeV->SetMarkerColor(kBlue);
  grPrStat5TeV->SetMarkerStyle(29);
  grPrStat5TeV->SetMarkerSize(1.7);
  grPrStat5TeV->SetLineColor(grPrStat5TeV->GetMarkerColor());
  grPrStat5TeV->Draw("PZ");

  grPrSyat5TeV->SetMarkerColor(grPrStat5TeV->GetMarkerColor());
  grPrSyat5TeV->SetMarkerStyle(grPrStat5TeV->GetMarkerStyle());
  grPrSyat5TeV->SetMarkerSize(grPrStat5TeV->GetMarkerSize());
  grPrSyat5TeV->SetLineColor(grPrStat5TeV->GetMarkerColor());
  grPrSyat5TeV->SetFillStyle(0);
  grPrSyat5TeV->Draw("5P");

  

  grKStat->SetMarkerStyle(21);
  grKStat->SetMarkerColor(6);
  grKStat->SetLineColor(6);
  grKStat->Draw("PZ");
  
  grKSyst->SetFillStyle(0);
  grKSyst->SetMarkerStyle(21);
  grKSyst->SetMarkerColor(6);
  grKSyst->SetLineColor(6);
  grKSyst->Draw("5P");
  

  grPhiStat->SetMarkerColor(1);
  grPhiStat->SetLineColor(1);
  grPhiStat->Draw("PZ");

  grPhiSyst->SetFillStyle(0);
  //grPhiSyst->SetMarkerStyle(21);
  grPhiSyst->SetMarkerColor(1);
  grPhiSyst->SetLineColor(1);
  grPhiSyst->Draw("5P");


  TLegend *l1 = DrawLegend(0.2,0.23,0.32,0.33);
  l1->SetTextSize(0.04);
  l1->SetTextFont(42);
  l1->AddEntry((TObject*)0,"#bf{ALICE Preliminary}","");
  l1->AddEntry((TObject*)0,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","");
  l1->Draw();
  
  TLegend *l01 = DrawLegend(0.2,0.18,0.3,0.2);
  l01->SetTextSize(0.04);
  l01->SetTextFont(42);
  l01->AddEntry((TObject*)0,"#scale[1.0]{Uncertainties: stat. (bars), sys. (boxes)}", "");
  l01->Draw();

  TLegend *l = DrawLegend(0.22,0.35,0.9,0.38);
  //l->SetHeader("Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  l->SetTextSize(0.05);
  l->SetNColumns(5);
  l->SetColumnSeparation(0.3);
  l->SetMargin(0.4);
  l->AddEntry(grPiStat5TeV,"#pi^{+}+#pi^{-}", "p");
  l->AddEntry(grKaStat5TeV,"K^{+}+K^{-}", "p");
  l->AddEntry(grKStat,"K*^{0}", "p");
  l->AddEntry(grPrStat5TeV,"p+#bar{p}", "p");
  l->AddEntry(grPhiStat,"#phi", "p");  
  l->Draw();
  

  c1->SaveAs("MeanPtComparison.gif");
  c1->Print("MeanPtComparison.eps");
  c1->Print("MeanPtComparison_pic.root");
  c1->cd();
}
TGraphErrors * PlotGraph(Int_t NdataPoint, Int_t MarkerColor = 1, Int_t MarkerStyle = 20, Double_t *X=0, Double_t *ErX=0,Double_t *Y=0, Double_t *ErY=0)
{
  TGraphErrors *gr = new TGraphErrors(NdataPoint, X ,Y , ErX, ErY);
  gr->SetTitle("");
  gr->SetMarkerStyle(MarkerStyle);
  gr->SetMarkerColor(MarkerColor);
  gr->SetMarkerSize(1.5);
  gr->GetYaxis()->SetTitle("#LT#it{p}_{T}#GT ( GeV/#{c})");
  gr->GetXaxis()->SetTitleSize(0.07);
  gr->GetXaxis()->SetTitleOffset(1.05);
  gr->GetXaxis()->SetTitleFont(42);
  //gr->GetXaxis()->CenterTitle(true);
  gr->GetXaxis()->SetLabelSize(0.05);
  gr->GetXaxis()->SetLabelFont(42);
  gr->GetXaxis()->SetNdivisions(505);
  gr->GetXaxis()->SetTitle("#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#eta|<0.5}");
  gr->GetYaxis()->SetTitleFont(42);
  //gr->GetYaxis()->CenterTitle(true);
  gr->GetYaxis()->SetTitleSize(0.07);
  gr->GetYaxis()->SetTitleOffset(1.2);
  gr->GetYaxis()->SetLabelSize(0.05);
  gr->GetYaxis()->SetLabelFont(42);
  gr->GetYaxis()->SetNdivisions(505);
  return gr;
}

//-----------------------------------------------------------
TCanvas *DrawFrame(TString o)
{
  TCanvas *c = new TCanvas(o.Data(),"Yield",10,10,600,600);
  c->SetLeftMargin(0.2);
  c->SetRightMargin(0.03);
  c->SetTopMargin(0.03);
  c->SetBottomMargin(0.18);
  return c;
}
//--------------------------------------------------------------
Double_t RatioError(Double_t X,Double_t ErX,Double_t Y,Double_t ErY)
{
  return (X/Y)*TMath::Sqrt((ErX/X)*(ErX/X) + (ErY/Y)*(ErY/Y));
}
//==========================================
TLine *DrawLine(Double_t xmin,Double_t ymin, Double_t xmax, Double_t ymax,Int_t lStyle, Int_t lColor, Int_t lWidth)
{
  TLine *line = new TLine(xmin,ymin,xmax,ymax);
  line->SetLineStyle(lStyle);
  line->SetLineColor(lColor);
  line->SetLineWidth(lWidth);
  //line->Draw();
  return line;
}
//==========================================
TLatex *DrawText(Double_t x, Double_t y,Int_t tColor,TString name)
{
  TLatex* tex = new TLatex(x,y,name.Data());
  tex->SetTextSize(0.04);
  tex->SetTextColor(tColor);
  tex->SetTextFont(42);
  //tex->Draw();
  return tex;
}
//=================================
TLegend *DrawLegend(Double_t x1,Double_t y1,Double_t x2,Double_t y2)
{
  TLegend *legend = new TLegend(x1,y1,x2,y2);
  legend->SetTextFont(42);
  legend->SetTextSize(0.04);
  legend->SetLineColor(0);
  legend->SetShadowColor(0);
  //legend->AddEntry(gr1,"(0 - 100) %","p");
  //legend->AddEntry(func1,"p_{0}[ 1 + 2 v_{2}^{obs} cos(2(#Phi - #Psi))]","l");
  return legend;
}
