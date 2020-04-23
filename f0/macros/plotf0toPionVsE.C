#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include <iostream>
#include "/Users/fbellini/alice/macros/cosmetics/SetStyle.C"
#include "/Users/fbellini/alice/macros/cosmetics/MakeUp.C"

void plotf0toPionVsE(Bool_t plotLogy = 1)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  //Set global style
  SetStyle();
  gStyle->SetOptStat(0);
  //TGaxis::SetMaxDigits(2);
  
  Double_t eee[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
  
  TCanvas *c1 = new TCanvas("c1","c1",900,700);
  //c1->SetLogx();
  c1->SetBottomMargin(0.15);
  c1->SetLeftMargin(0.18);
  c1->SetRightMargin(0.07);
  c1->SetTopMargin(0.07);
  if (plotLogy) gPad->SetLogy();
  gPad->SetLogx();
  TH1F *frame = new TH1F("frame","frame",1300,9,13000);
  frame->SetMinimum(0.00008);
  if (plotLogy) frame->SetMaximum(.2);
  else frame->SetMaximum(.1);
  frame->GetYaxis()->SetTitleOffset(1.45);
  frame->GetXaxis()->SetTitleOffset(1.25);
  frame->GetXaxis()->SetLabelSize(0.045*1.2);
  frame->GetYaxis()->SetLabelSize(0.045*1.2);
  frame->GetXaxis()->SetTitleSize(0.045*1.2);
  frame->GetYaxis()->SetTitleSize(0.045*1.2);
  frame->SetTitle("; #sqrt{#it{s}_{NN}} (GeV); Ratio");
  frame->Draw();
  
  //pp
  //SPS pp 27.5 GeV: Aguilar-Benitez M et al (1991) Z. Phys. C 50 405,
  //http://cds.cern.ch/record/216837/files/cer-000129214.pdf
  //F0(980) sigma_T = (0.74 ± 0.26) mb
  //π+ sigma_T = (134.4 ± 3.4) mb
  Double_t ratio_pp_275 = 0.74/134.4;
  Double_t ratio_pp_275_err = TMath::Sqrt(TMath::Power(0.26/0.74, 2.0) + TMath::Power(3.4/134.4, 2.0)) * ratio_pp_275;
  Double_t x_pp_275[1] = {27.5};
  Double_t y_pp_275[1] = {ratio_pp_275};  //f0/π+
  Double_t ey_pp_275[1] = {ratio_pp_275_err};
  TGraphErrors *gr_pp_275 = new TGraphErrors(1, x_pp_275, y_pp_275, eee, ey_pp_275);
  gr_pp_275->SetMarkerStyle(24);
  gr_pp_275->SetMarkerSize(1.7);
  gr_pp_275->Draw("P");

  //ee 29 GeV
  //f0 = 0.06 ±0.03 Abachi, S., et aI., Phys. Rev. Lett. 57:1 990 (1986) 
  //π+ = 10.3 ± 0.4  Derrick, M., et aI., Phys. Rev. D35:2639 (1987) 
  //https://www-annualreviews-org.ezproxy.cern.ch/doi/pdf/10.1146/annurev.ns.38.120188.001431
  Double_t ratio_pp_29 = 0.06/10.3;
  Double_t ratio_pp_29_err = TMath::Sqrt(TMath::Power(0.03/0.06, 2.0) + TMath::Power(0.4/10.3, 2.0)) * ratio_pp_29;
  Double_t x_ee_29[1] = {29};
  Double_t y_ee_29[1] = {ratio_pp_29};
  Double_t ey_ee_29[1] = {ratio_pp_29_err};
  TGraphErrors *gr_ee_29 = new TGraphErrors(1, x_ee_29, y_ee_29, eee, ey_ee_29);
  gr_ee_29->SetMarkerStyle(25);
  gr_ee_29->SetMarkerSize(1.7);
  gr_ee_29->Draw("P");

  //pp
  //ISR pp 52.5 GeV: Drijard D et al (1981) Z. Phys. C 9 293,
  //https://link.springer.com/article/10.1007/BF01548764
  //Sigma(pp —> f +x) / Sigma(pp —> π- + x) = 0.025 ± 0.005
  Double_t x_pp_525[1] = {52.5};
  Double_t y_pp_525[1] = {0.025};
  Double_t ey_pp_525[1] = {0.005};
  TGraphErrors *gr_pp_525 = new TGraphErrors(1, x_pp_525, y_pp_525, eee, ey_pp_525);
  gr_pp_525->SetMarkerStyle(24);
  gr_pp_525->SetMarkerSize(1.7);
  gr_pp_525->Draw("P");


  //ee √s = 91 TeV
  //344 P.V. Chliapnikov, Physics Letters B 462 1999 341–353
  //https://ac.els-cdn.com/S0370269399009107/1-s2.0-S0370269399009107-main.pdf?_tid=2b212de3-7a15-4e72-9ee1-cd5fb6e8f10a&acdnat=1523958001_fa3ae38d9c067fa1a25503e0fcbfc3fd
  Double_t ratio_ee_91 = 0.147/9.53;
  Double_t ratio_ee_91_err = TMath::Sqrt(TMath::Power(0.011/0.147, 2.0) + TMath::Power(0.12/8.53, 2.0)) * ratio_ee_91;
  Double_t x_ee_91[1] = {91};
  Double_t y_ee_91[1] = {ratio_ee_91};
  Double_t ey_ee_91[1] = {ratio_ee_91_err};
  TGraphErrors *gr_ee_91 = new TGraphErrors(1, x_ee_91, y_ee_91, eee, ey_ee_91);
  gr_ee_91->SetMarkerStyle(25);
  gr_ee_91->SetMarkerSize(1.7);
  gr_ee_91->Draw("P");

  //pp
  //http://iopscience.iop.org/article/10.1088/0954-3899/30/1/069/pdf P.Fachini
  //STAR preliminary 2003 - unpublished
  Double_t x_pp_200[1] = {200};
  Double_t y_pp_200[1] = {0.042};
  Double_t ey_pp_200[1] = {0.022};
  TGraphErrors *gr_pp_200 = new TGraphErrors(1, x_pp_200, y_pp_200, eee, ey_pp_200);
  gr_pp_200->SetMarkerStyle(34);
  gr_pp_200->SetMarkerSize(1.7);
  gr_pp_200->Draw("P");

  /*
  //AuAu
  Double_t x_AuAu_200[1] = {200};
  Double_t y_AuAu_200[1] = {0.169};
  Double_t ey_AuAu_200[1] = {0.037};
  TGraphErrors *gr_AuAu_200 = new TGraphErrors(1, x_AuAu_200, y_AuAu_200, eee, ey_AuAu_200);
  gr_AuAu_200->SetMarkerStyle(29);
  gr_AuAu_200->SetMarkerSize(1.7);
  gr_AuAu_200->Draw("P");
  */

  //pp
  //ALICE preliminary candidate  //π++π- = 4.1375 ± 0.000504 ± 0.167941 (preliminary values, from paper draft)
  //0.03780 ± 0.0002 (stat.) ± 0.0065 (sys.)
  Double_t ratio_pp_5200 = (0.46/0.5)*2.0*0.03780/4.137502;
  Double_t err_f0_tot = TMath::Sqrt(TMath::Power(0.0002, 2.0) + TMath::Power(0.0065, 2.0));
  Double_t err_pi_tot = TMath::Sqrt(TMath::Power(0.000504, 2.0) + TMath::Power(0.167941, 2.0));
  Double_t ratio_pp_5200_err = TMath::Sqrt(TMath::Power(err_f0_tot/0.03780, 2.0) + TMath::Power( err_pi_tot/4.137502, 2.0)) * ratio_pp_5200;
  Double_t x_pp_5200[1] = {5200};
  Double_t y_pp_5200[1] = {ratio_pp_5200};
  Double_t ey1_pp_5200[1] = {ratio_pp_5200_err};
  TGraphErrors *gr_pp_5200 = new TGraphErrors(1, x_pp_5200, y_pp_5200, eee, ey1_pp_5200);
  gr_pp_5200->SetLineColor(2);
  gr_pp_5200->SetMarkerColor(kRed);
  gr_pp_5200->SetMarkerStyle(20);
  gr_pp_5200->SetMarkerSize(1.7);
  gr_pp_5200->Draw("P");
  Printf("f0/pi = %6.4f +/- %6.4f", ratio_pp_5200, ratio_pp_5200_err);


  //Becattini
  //Pon yield Eur.Phys.J.C56:493-510,2008 // https://arxiv.org/abs/0805.0964
  //f0 sent by priv. comm. to P. Abreu et al, PLB 449 (1999) 364–382 is slightly different, not used
  TLine * shmBecattini_ee_91 = new TLine(70, 0.0779/8.48, 130, 0.0772/8.48);
  shmBecattini_ee_91 ->SetLineColor(kAzure+7);
  shmBecattini_ee_91 ->SetLineWidth(5);
  shmBecattini_ee_91 ->SetLineStyle(3);
  shmBecattini_ee_91->Draw();
  
  //Becattini, Heinz https://link.springer.com/content/pdf/10.1007%2Fs002880050551.pdf
  TLine * shmBecattini_pp_7000 = new TLine(4500, 0.0094, 8000, 0.0094);
  shmBecattini_pp_7000 ->SetLineColor(kBlue);
  shmBecattini_pp_7000 ->SetLineWidth(5);
  shmBecattini_pp_7000->Draw();


  TLegend * legmodel = new TLegend(0.2,0.83,0.6,0.9);
  legmodel->AddEntry(shmBecattini_ee_91,"Becattini et al., e^{+}e^{-} [Eur.Phys.J. C56 (2008) 493]", "l");
  legmodel->AddEntry(shmBecattini_pp_7000,"Becattini, pp [priv. comm.] ", "l"); //[Z.Phys. C76 (1997) 269]
  legmodel->SetFillColor(kWhite);
  legmodel->SetBorderSize(0);
  legmodel->SetTextSize(0.032);
  legmodel->SetFillStyle(0);
  legmodel->Draw();

    //Create legend:0.435,0.82,0.95,0.95
  TLegend * leg = new TLegend(0.2,0.18,0.6,0.45);
  leg->AddEntry(gr_pp_5200, "2f_{0}/(#pi^{+}+#pi^{-}), pp [this analysis]","p");
  leg->AddEntry(gr_pp_200, "f_{0}/#pi, pp [STAR preliminary (2003), unpublished]","p");
  leg->AddEntry(gr_pp_525, "f_{0}/#pi^{-}, pp [Z.Phys. C9 (1981) 293]","p");
  leg->AddEntry(gr_pp_275, "f_{0}/#pi^{+}, pp [Z.Phys. C50 (1991) 405]","p");
  leg->AddEntry(gr_ee_91, "f_{0}/#pi^{+}, e^{+}e^{-} [Phys. Lett. B 462 (1999) 341]","p");
  leg->AddEntry(gr_ee_29, "f_{0}/#pi^{+}, e^{+}e^{-} [Ann.Rev.Nucl.Part.Sci. 38 (1988) 279]","p");

  // leg->AddEntry(gr_pp_684, "#rho/#pi (pp)","p");
  // leg->AddEntry(gr_ee_1045, "#rho/#pi^{-} (e^{+}e^{-})","p");
  // leg->AddEntry(gr_pip_196, "#rho/#pi^{-} (#pi^{-}+p)","p");
  // leg->AddEntry(gr_kp_782, "#rho/#pi (K^{+}p)","p");
  // leg->AddEntry(gr_pp_2760, "#rho/#pi (pp, ALICE)","p");
  leg->SetNColumns(1);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.032);
  leg->SetFillStyle(0);
  leg->Draw();

  TString text="Uncertainties: #sqrt{stat.^{2} + sys.^{2}}";
  TPaveText * pave = new TPaveText(0.4,0.2,0.85,0.26,"brNDC");
  pave->SetBorderSize(0);
  pave->SetFillColor(kWhite);
  pave->SetTextColor(kBlack);
  pave->SetTextFont(42);
  pave->SetTextSize(0.0);
  pave->InsertText(text.Data());
  //pave->Draw();

  TString nameImg("f0toPionVsE");
  if (plotLogy) nameImg.Append(".log");
  c1->SaveAs((nameImg.Append(".pdf")));
  return;
}
