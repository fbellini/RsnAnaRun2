// fignum:
// 1 = spectra
// 2 = to-pion ratios
// 3 = scaling plot

void LoadPlotsStyle(int fignum){
  
  // general
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);

  // pad margins
  gStyle->SetPadTopMargin(0.015);
  gStyle->SetPadBottomMargin(0.08);
  if(fignum==3) {
    gStyle->SetPadTopMargin(0.015/0.5475);
    gStyle->SetPadBottomMargin(0.08/0.5475);
  }
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadLeftMargin(0.18);

  // figure-dependent settings
  if(fignum==1) {
    gStyle->SetOptLogy();
  }
  else{
    //TGaxis::SetMaxDigits(2);
  }
  
  // frame ticks
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetTickLength(0.02,"X");
  if(fignum==3) gStyle->SetTickLength(0.02/0.5475,"X");
  gStyle->SetTickLength(0.03,"Y");

  // x-axis style
  gStyle->SetTitleFont(42,"X");
  gStyle->SetTitleSize(0.055,"X");
  if(fignum==3) gStyle->SetTitleSize(0.055/0.5475/1.55,"X");
  gStyle->SetTitleOffset(0.65,"X");
  if(fignum==3) gStyle->SetTitleOffset(0.65*1.55,"X");
  gStyle->SetLabelFont(42,"X");
  gStyle->SetLabelSize(0.050,"X");
  if(fignum==3) gStyle->SetLabelSize(0.050/0.5475/1.55,"X");
  gStyle->SetLabelOffset(-0.012,"X");

  // y-axis style
  gStyle->SetTitleFont(42,"Y");
  gStyle->SetTitleSize(0.055,"Y");
  if(fignum==3) gStyle->SetTitleSize(0.055/0.5475/1.55,"Y");
  gStyle->SetTitleOffset(1.45,"Y");
  if(fignum==3) gStyle->SetTitleOffset(1.45*0.9,"Y");
  gStyle->SetLabelFont(42,"Y");
  gStyle->SetLabelSize(0.050,"Y");
  if(fignum==3) gStyle->SetLabelSize(0.050/0.5475/1.55,"Y");
  gStyle->SetLabelOffset(0.01,"Y");

}
