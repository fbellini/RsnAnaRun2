#if !defined (__CINT__) || defined (__CLING__)
#include "/Users/fbellini/alice/macros/SetStyle.C"
#include "/Users/fbellini/alice/macros/Beautify.C"
#endif


TCanvas *DrawFrame(TString o ="c");
Double_t RatioError(Double_t X,Double_t ErX,Double_t Y,Double_t ErY);
TLine *DrawLine(Double_t xmin = 0,Double_t ymin = 1, Double_t xmax = 1, Double_t ymax = 1,Int_t lStyle = 1, Int_t lColor = 1, Int_t lWidth = 1);
TLatex *DrawText(Double_t x = 0, Double_t y = 0,Int_t tColor = 2,TString name = "label");
TLegend *DrawLegend(Double_t x1,Double_t y1,Double_t x2,Double_t y2);
TGraphErrors * MeanPt_Xe(TString part = "Pi", Bool_t sys = 0);
TGraph * GetBayesPredictionsMeanPt_PbPb502TeV(TString particle = "Pion");

void MeanPt(Bool_t plotPbPb = 1, Bool_t plotBayes = 1, Bool_t plotKstar = 0)
{
  SetStyle();
  gStyle->SetErrorX(0);
  gStyle->SetOptStat(0);

  
  // Color_t color[5][2] = {kRed-8, kRed+2,
  // 			 kSpring+2, kGreen+3,
  // 			 kBlue-8,  kBlue+1,
  // 			 kMagenta-8, kMagenta+1,
  // 			 kGray, kBlack};

  // Color_t color[5][2] = {kBlack, kRed+1,
  // 			 kBlack, kGreen+2,
  // 			 kBlack, kBlue+1,
  // 			 kBlack, kAzure+1,
  // 			 kBlack, kOrange+1};
  
  Color_t color[5][2] = {kRed-3, kRed+1,
			 kGreen-3, kGreen+2,
			 kBlue-3, kBlue+1,
			 kMagenta-3, kMagenta+2,
			 kOrange-3, kOrange+1};

  Int_t marker[5][2] = { 27, 33,
  			 28, 34,
  			 30, 29,
  			 25, 21,
  			 24, 20};
  
  //--------------------
  //Xe-Xe 5.44
  //--------------------
  //
  //K* from Sourav, 10.05.2018 by email
  //TFile *fileks_Xe = TFile::Open("FinalKStardNdYMeanPtV0MNew.root");
  Double_t MeanPT[6]={1.39199,1.29832,1.20329,1.04989,1.28681,1.15611};
  Double_t MeanPT_Err_Stat[6]={0.0855817,0.0730252,0.0680304,0.0626365,0.0615943,0.0565969};
  Double_t MeanPT_Err_Syst[6] = {0.0463404,0.0488587,0.0478421,0.055248,0.0471142,0.051664};
  Double_t MNch[6] = {748.3,258.4,92,23,212,37.12};
  Double_t MNchErr[6] = {12.84,4,2.0,1.0,3.66,1.136};
  
  TGraphErrors * gKs_meanPT_XeXe5TeV_stat = new TGraphErrors(6, MNch, MeanPT, MNchErr, MeanPT_Err_Stat);
  TGraphErrors * gKs_meanPT_XeXe5TeV_syst = new TGraphErrors(6, MNch, MeanPT, MNchErr, MeanPT_Err_Syst);
  
  TFile *filephi_Xe = TFile::Open("Preliminary_PhiMeanPt_XeXe544TeV.root");
  TGraphErrors * gPhi_meanPT_XeXe5TeV_stat = (TGraphErrors*)filephi_Xe->Get("stat");
  TGraphErrors * gPhi_meanPT_XeXe5TeV_syst = (TGraphErrors*)filephi_Xe->Get("sys");
  
  //TFile *fPiXeXe544 = TFile::Open("20180502_YieldAndMeanPtSumPiwSys.root");
  TGraphErrors  *grPiStatXeXe544  = (TGraphErrors*) MeanPt_Xe("Pi", 0);
  TGraphErrors  *grPiSystXeXe544  = (TGraphErrors*) MeanPt_Xe("Pi", 1);
 
  //TFile *fProXeXe544 = TFile::Open("20180502_YieldAndMeanPtSumPwSys.root");
  TGraphErrors  *grProStatXeXe544  = (TGraphErrors*) MeanPt_Xe("P", 0);
  TGraphErrors  *grProSystXeXe544  = (TGraphErrors*) MeanPt_Xe("P", 1);

  //TFile *fKaXeXe544 = TFile::Open("20180502_YieldAndMeanPtSumKwSys.root");
  TGraphErrors  *grKaStatXeXe544  = (TGraphErrors*) MeanPt_Xe("K", 0);
  TGraphErrors  *grKaSystXeXe544  = (TGraphErrors*) MeanPt_Xe("K", 1);

  
  //--------------------
  //Pb-Pb 2.76
  //--------------------
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
  TFile *fileks = TFile::Open("PbPb5TeV/MeanPt/FinalKStardNdYMeanPtV0MNew.root");
  TGraphErrors * grKStat = (TGraphErrors*)fileks->Get("gMeanPt_stat_KStar");
  TGraphErrors * grKSyst = (TGraphErrors*)fileks->Get("gMeanPt_syst_KStar");

  TFile *filephi = TFile::Open("PbPb5TeV/MeanPt/FinalPhidNdYMeanPtV0MNew.root");
  TGraphErrors * grPhiStat = (TGraphErrors*)filephi->Get("gMeanPt_stat_phi");
  TGraphErrors * grPhiSyst = (TGraphErrors*)filephi->Get("gMeanPt_syst_phi");
  for(Int_t i = 0; i<grPhiStat->GetN();i++)grPhiStat->SetPointError(i,0,grPhiStat->GetErrorY(i));
  for(Int_t i = 0; i<grKStat->GetN();i++)grKStat->SetPointError(i,0,grKStat->GetErrorY(i));
  
  TFile *f5TeV = TFile::Open("PbPb5TeV/MeanPt/YieldAndMeanPtPbPb5TeV.root");
  f5TeV->ls();
  TGraphErrors  *grPiStat5TeV  = (TGraphErrors*)f5TeV->Get("PionSumByTwoMeanPtStat");
  TGraphErrors  *grPiSyat5TeV    = (TGraphErrors*)f5TeV->Get("PionSumByTwoMeanPtSyst");
  TGraphErrors  *grKaStat5TeV  = (TGraphErrors*)f5TeV->Get("KaonSumByTwoMeanPtStat");
  TGraphErrors  *grKaSyat5TeV    = (TGraphErrors*)f5TeV->Get("KaonSumByTwoMeanPtSyst");
  TGraphErrors  *grPrStat5TeV  = (TGraphErrors*)f5TeV->Get("ProtonSumByTwoMeanPtStat");
  TGraphErrors  *grPrSyat5TeV    = (TGraphErrors*)f5TeV->Get("ProtonSumByTwoMeanPtSyst");

  // Bayesian analysis
  TGraph * gBayesPion = (TGraph *) GetBayesPredictionsMeanPt_PbPb502TeV("Pion");
  gBayesPion->SetLineColor(color[0][0]);
  gBayesPion->SetLineWidth(3);
  TGraph * gBayesKaon = (TGraph *) GetBayesPredictionsMeanPt_PbPb502TeV("Kaon");
  gBayesKaon->SetLineColor(color[1][0]);
  gBayesKaon->SetLineWidth(3);
  TGraph * gBayesProton = (TGraph *) GetBayesPredictionsMeanPt_PbPb502TeV("Proton");
  gBayesProton->SetLineColor(color[2][0]);
  gBayesProton->SetLineWidth(3);

  TH1F* hf = new TH1F("frame","", 2.E3, 1., 3.0E3);
  hf->SetMinimum(0.01);
  hf->SetMaximum(1.55);
  hf->SetXTitle("#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta} | < 0.5}");
  hf->SetLineColor(1);
  hf->GetXaxis()->SetTitleSize(0.06);
  hf->GetXaxis()->SetLabelSize(0.045);
  hf->GetXaxis()->SetTitleOffset(1.1);
  hf->SetYTitle("#LT#it{p}_{T}#GT (GeV/#it{c})");
  hf->GetYaxis()->SetTitleOffset(1.2);
  hf->GetYaxis()->SetTitleSize(0.06);
  hf->GetYaxis()->SetLabelSize(0.045);
  hf->GetYaxis()->SetTitleOffset(1.);
  hf->SetNdivisions(409,"x");
  hf->SetNdivisions(407,"y");
  if (plotPbPb) hf->GetYaxis()->SetRangeUser(0.01, 2.25);
  else hf->GetYaxis()->SetRangeUser(0.97, 1.78);
  hf->GetXaxis()->SetRangeUser(13., 3.e3);
  
  TCanvas *c1 = new TCanvas("c1","c",10,10, 700,800);
  c1->SetBottomMargin(0.15);
  c1->SetLeftMargin(0.16);
  c1->SetTopMargin(0.02);
  c1->SetRightMargin(0.02);
  c1->SetTicks(1,1);
  c1->cd();
  gPad->SetLogx();
  hf->Draw();
  
  if (plotPbPb) {
    BeautifyGraph(grPiStat5TeV, color[0][0], color[0][0], 0, 1, 1, marker[0][0], 2.3); 
    BeautifyGraph(grPiSyat5TeV, color[0][0], color[0][0], 0, 1, 1, marker[0][0], 2.3); 
    grPiStat5TeV->Draw("PZ same");
    grPiSyat5TeV->Draw("5P same");

    BeautifyGraph(grKaStat5TeV, color[1][0], color[1][0], 0, 1, 1, marker[1][0], 2.3); 
    BeautifyGraph(grKaSyat5TeV, color[1][0], color[1][0], 0, 1, 1, marker[1][0], 2.3); 
    grKaStat5TeV->Draw("PZ same");
    grKaSyat5TeV->Draw("5P same");

    BeautifyGraph(grPrStat5TeV, color[2][0], color[2][0], 0, 1, 1, marker[2][0], 2.3); 
    BeautifyGraph(grPrSyat5TeV, color[2][0], color[2][0], 0, 1, 1, marker[2][0], 2.3); 
    grPrStat5TeV->Draw("PZ same");
    grPrSyat5TeV->Draw("5P same");

    if (!plotBayes) {
      BeautifyGraph(grKStat, color[3][0], color[3][0], 0, 1, 1, marker[3][0], 1.5); 
      BeautifyGraph(grKSyst, color[3][0], color[3][0], 0, 1, 1, marker[3][0], 1.5); 
      if (plotKstar) {
	grKStat->Draw("PZ same");
	grKSyst->Draw("5P same");
      }
      BeautifyGraph(grPhiStat, color[4][0], color[4][0], 0, 1, 1, marker[4][0], 1.7); 
      BeautifyGraph(grPhiSyst, color[4][0], color[4][0], 0, 1, 1, marker[4][0], 1.7);
      grPhiStat->Draw("PZ same");
      grPhiSyst->Draw("5P same");
    }
  }

  if (plotBayes){
    gBayesPion->Draw("l same");
    gBayesKaon->Draw("l same");
    gBayesProton->Draw("l same");
  }
  
  BeautifyGraph(grPiStatXeXe544, color[0][1], color[0][1], 0, 1, 1, marker[0][1], 2.3); 
  BeautifyGraph(grPiSystXeXe544, color[0][1], color[0][1], 0, 1, 1, marker[0][1], 2.3); 
  grPiSystXeXe544->Draw("E2 same");
  grPiStatXeXe544->Draw("PZ same");

  BeautifyGraph(grKaStatXeXe544, color[1][1], color[1][1], 0, 1, 1, marker[1][1], 2.3); 
  BeautifyGraph(grKaSystXeXe544, color[1][1], color[1][1], 0, 1, 1, marker[1][1], 2.3); 
  grKaSystXeXe544->Draw("E2 same");
  grKaStatXeXe544->Draw("PZ same");

  BeautifyGraph(grProStatXeXe544, color[2][1], color[1][1], 0, 1, 1, marker[2][1], 2.3); 
  BeautifyGraph(grProSystXeXe544, color[2][1], color[1][1], 0, 1, 1, marker[2][1], 2.3); 
  grProSystXeXe544->Draw("E2 same");
  grProStatXeXe544->Draw("PZ same");

  if (!plotBayes) {
    BeautifyGraph(gKs_meanPT_XeXe5TeV_stat, color[3][1], color[3][1], 0, 1, 1, marker[3][1], 1.5); 
    BeautifyGraph(gKs_meanPT_XeXe5TeV_syst, color[3][1], color[3][1], 0, 1, 1, marker[3][1], 1.5);
    if (plotKstar) {
      gKs_meanPT_XeXe5TeV_syst->Draw("e2 same");
      gKs_meanPT_XeXe5TeV_stat->Draw("pz same");
    }
    BeautifyGraph(gPhi_meanPT_XeXe5TeV_stat, color[4][1], color[4][1], 0, 1, 1, marker[4][1], 1.7); 
    BeautifyGraph(gPhi_meanPT_XeXe5TeV_syst, color[4][1], color[4][1], 0, 1, 1, marker[4][1], 1.7);
    gPhi_meanPT_XeXe5TeV_syst->Draw("e2 same");
    gPhi_meanPT_XeXe5TeV_stat->Draw("pz same");
  }

  TLegend *l1 = DrawLegend(0.16,0.8,0.42,0.95);
  l1->SetTextSize(0.035);
  l1->SetTextFont(42);
  l1->SetBorderSize(0);
  l1->SetFillStyle(0);
  l1->AddEntry((TObject*)0,"#bf{ALICE Preliminary}","");
  l1->AddEntry((TObject*)0,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV (open markers)","");
  l1->AddEntry((TObject*)0,"Xe-Xe #sqrt{#it{s}_{NN}} = 5.44 TeV (full markers)","");
  l1->Draw();

  TGraph * dummy = new TGraph(1); dummy->SetLineWidth(3);

  TLegend *lduke = DrawLegend(0.2,0.18,0.5,0.28);
  lduke->SetHeader("Bayesian analysis (Duke)");
  lduke->SetTextSize(0.035);
  lduke->SetTextFont(42);
  lduke->SetBorderSize(0);
  lduke->SetFillStyle(0);
  lduke->AddEntry(dummy, "Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","l");
  if (plotBayes) lduke->Draw();

  
  TLegend *l01 = DrawLegend(0.2,0.19,0.3,0.2);
  l01->SetTextSize(0.035);
  l01->SetTextFont(42);
  l01->AddEntry((TObject*)0,"#scale[1.0]{Uncertainties: stat. (bars), sys. (boxes)}", "");
  //l01->Draw();

  TLegend *l = DrawLegend(0.22,0.72,0.9,0.8);
  //l->SetHeader("Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  l->SetTextSize(0.04);
  l->SetNColumns(5);
  l->SetFillStyle(0);
  l->SetColumnSeparation(0.3);
  l->SetMargin(0.4);
  l->AddEntry(grPiStatXeXe544,"#pi^{+}+#pi^{-}", "p");
  l->AddEntry(grKaStatXeXe544,"K^{+}+K^{-}", "p");
  if (plotKstar && !plotBayes) l->AddEntry(gKs_meanPT_XeXe5TeV_stat,"K*^{0}", "p");
  l->AddEntry(grProStatXeXe544,"p+#bar{p}", "p");
  if (!plotBayes) l->AddEntry(gPhi_meanPT_XeXe5TeV_stat,"#phi", "p");
  l->Draw();

  TString nameimg = "IdMeanPt_CompareXe2Pb";
  if (plotBayes) nameimg.Append("_wBayes");
  c1->Print(Form("%s.eps", nameimg.Data()));
  c1->Print(Form("%s.pdf", nameimg.Data()));
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

TGraphErrors * MeanPt_Xe(TString part, Bool_t sys)
{
  TString XeFileName = Form("Preliminary_XeXe_YieldAndMeanPtSum%swSys.root", part.Data());
  TFile * XeFile = TFile::Open(XeFileName.Data());
  if (!XeFile) return 0;
  TString statName = "hStdPtSummed";
  TString sysName = "hStdPtSysSummed";
  if (part.Contains("K")) {
    statName.Append("Kaon");
    sysName.Append("Kaon");
  } else if (part.Contains("Pi")) {
    statName.Append("Pion");
    sysName.Append("Pion");
  } else {
    statName.Append("Proton");
    sysName.Append("Proton");
  }

  TH1F *stat = (TH1F*) XeFile->Get(statName.Data());
  TH1F *syst = (TH1F*) XeFile->Get(sysName.Data());

  Float_t dNdeta[9] = {1167., 939., 706., 478., 315., 198., 118., 64.7, (32.0 + 13.3) / 2.0};
  Float_t dNdetaErr[9] = {26., 24., 17., 11., 8., 5., 3., 2., (1.3 + 0.9) / 2.0};
  Float_t meanPt[9];
  Float_t meanPtStat[9];
  Float_t meanPtSys[9];

  for (int j=0;j<9;j++){
    meanPt[j] = stat->GetBinContent(j+1);
    meanPtStat[j] = stat->GetBinError(j+1);
    meanPtSys[j] = syst->GetBinError(j+1);
  }
  
  TGraphErrors * graph;
  if (sys) graph = new TGraphErrors(9, dNdeta, meanPt, dNdetaErr, meanPtSys);
  else  graph = new TGraphErrors(9, dNdeta, meanPt, dNdetaErr, meanPtStat);

  return graph;
}
  
TGraph * GetBayesPredictionsMeanPt_PbPb502TeV(TString particle)
{
  /****************************************
   // bayesian analysis, Duke group      
   //Data from
   //****************************************/
  TFile * fBayes = TFile::Open(Form("Bayesian/PbPb5020/BayesianPredictionMeanPt%s.root", particle.Data()));
  TGraph * gBayes = (TGraph*) fBayes->Get("Graph");
  return gBayes;
}
