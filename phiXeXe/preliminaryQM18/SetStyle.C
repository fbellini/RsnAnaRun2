void SetStyle(UInt_t statOpt = 1111, Bool_t logz = 1)
{
  //
  //set style
  //
  gStyle->Reset("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(8,0);
  gStyle->SetPalette(kBird); //kDeepSea, kBlueGreenYellow, kRainBow, kBird =default
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.17);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetHistLineColor(kBlack);
  gStyle->SetHistLineWidth(0);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(1);

  gStyle->SetLabelOffset(0.01,"yz");
  gStyle->SetTitleOffset(1.25,"yz");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetTitleOffset(1.2,"x");

  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.06,"xyz");
  gStyle->SetLabelSize(0.05,"xyz");

  gStyle->SetTitleFillColor(kWhite);
  //gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  gStyle->SetEndErrorSize(0); //sets in #of pixels the lenght of the tick at the end of the error bar
  //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y"); 

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  //  gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);
  gStyle->SetOptStat(statOpt);

}
