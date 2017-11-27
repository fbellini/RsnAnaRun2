/****************************************************************
Author: fbellinicern.ch
First version: 27.01.2017
****************************************************************
This macro reads data from an external simple text file, that has to be formatted as
cent-min cent-max value stat sys uncorr

According to requested format it can display the data in several different formats:
- "array"
- "graph" (TGraph in root file)
- "histo"
- "latex" (table) --> to be implemented
****************************************************************/

#include "Riostream.h"

void ReaddNdetaFromFile(TString file = "input.dat",
		       TString outputOpt = "array",
		       Int_t icolor = 1,
		       TString system = "PbPb",
		       Float_t energy = 2.76)
{
  Color_t colE[4] = {kBlack, kBlue+1, kGreen+1, kRed+1};
  Color_t colF[4] = {kGray+1, kAzure+1, kSpring-5, kRed-7};
  Int_t  markE[4] = {24, 20, 21, 34};

  const Int_t Np = 20; 
  //centrlaity_low, centrality_up, y, stat_y, sys_y, uncorr_y
  Double_t cl[Np], cu[Np], ccentral[Np], chalfwidth[Np];
  Double_t y[Np], uy[Np], sy[Np], np[Np], unp[Np];
      
  TString foutname = "";
  TFile * f = 0x0;

  //Read data from an external text file
  //The file has to be formatted as
  //cent-min cent-max value stat sys uncorr
  ifstream in;
  in.open(Form("%s",file.Data()));
  Int_t nlines = 0;
  TNtuple *ntuple = new TNtuple("ntuple","data from file","cl:cu:y:uy:sy:np:unp");
  
  while (nlines<Np) {
    in >> cl[nlines] >> cu[nlines] >> y[nlines] >> uy[nlines] >> sy[nlines] >> np[nlines] >> unp[nlines];
    chalfwidth[nlines] = (cu[nlines]-cl[nlines])/2;
    ccentral[nlines] = cl[nlines]+chalfwidth[nlines];
    
    if (!in.good()) break;
    Printf("cent_low=%8f, cent_up=%8f, val=%8f, stat=%8f, sys=%8f, npart=%8f, err=%8f",cl[nlines],cu[nlines],y[nlines],uy[nlines],sy[nlines],np[nlines], unp[nlines]);
    ntuple->Fill(cl[nlines],cu[nlines],y[nlines],uy[nlines],sy[nlines],np[nlines], unp[nlines]);
    nlines++;
  }
 
  Printf("Reading file %s: found %d points\n", file.Data(), nlines); 
  
  if (outputOpt.Contains("array")) {
    printf("Double_t x[%i] = {", nlines);
    for (Int_t j = 0; j<nlines; j++){ printf("%8.5f, ", cl[j]);}
    printf("}; \nDouble_t y[%i] = {", nlines);
    for (Int_t j = 0; j<nlines; j++){ printf("%8.5f, ", y[j]);}
    printf("}; \nDouble_t uy[%i] = {");
    for (Int_t j = 0; j<nlines; j++){ printf("%8.5f, ", uy[j]);}
    printf("}; \nDouble_t sy[%i] = {", nlines);
    for (Int_t j = 0; j<nlines; j++){ printf("%8.5f, ", sy[j]);}
    printf("}; \nDouble_t np[%i] = {", nlines);
    for (Int_t j = 0; j<nlines; j++){ printf("%8.5f, ", np[j]);}
    printf("}; \nDouble_t unp[%i] = {", nlines);
    for (Int_t j = 0; j<nlines; j++){ printf("%8.5f, ", unp[j]);}
    printf("};\n");
  } 
  
  if (outputOpt.Contains("graph")) {
    foutname = Form("%s_%s%3.2fTeV_graph.root", file.Data(), system.Data(), energy);
    f = new TFile(foutname.Data(),"RECREATE"); 

    //create TGraphs
    TGraphErrors *gstat = new TGraphErrors(nlines, ccentral, y, chalfwidth, uy);
    gstat->SetName("stat");
    gstat->SetTitle("value vs cent., statistical unc");
    gstat->SetMarkerStyle(markE[icolor]);
    gstat->SetMarkerColor(colE[icolor]);
    gstat->SetMarkerSize(1.4);
    gstat->SetFillStyle(0);
    gstat->SetFillColorAlpha(colF[icolor], 0.4);
    gstat->SetLineWidth(1);
    gstat->SetLineColor(colE[icolor]);

    TGraphErrors *gsys = new TGraphErrors(nlines, ccentral, y, chalfwidth, sy);
    gsys->SetName("sys");
    gsys->SetTitle("value vs cent., systematics unc");
    gsys->SetMarkerStyle(markE[icolor]);
    gsys->SetMarkerColor(colE[icolor]);
    gsys->SetMarkerSize(1.4);
    gsys->SetFillStyle(0);
    gsys->SetFillColorAlpha(colF[icolor], 0.4);
    gsys->SetLineWidth(1);
    gsys->SetLineColor(colE[icolor]);

    //Draw
    //Create frame
    TH1F* hf = new TH1F("frame","", 20, 0., 100.0.);
    hf->SetXTitle("centrality (%)");
    hf->SetLineColor(1);
    hf->GetXaxis()->SetTitleSize(0.06);
    hf->GetXaxis()->SetLabelSize(0.045);
    hf->GetXaxis()->SetTitleOffset(1.);
    hf->SetNdivisions(510,"x");
    hf->SetYTitle("value");
    hf->GetYaxis()->SetTitleOffset(1.2);
    hf->GetYaxis()->SetTitleSize(0.06);
    hf->GetYaxis()->SetLabelSize(0.045);
    hf->GetYaxis()->SetTitleOffset(1.);
    hf->SetNdivisions(509,"y");
      
    //display
    gStyle->SetOptStat(0);
    TCanvas* c43=new TCanvas("c3","",12,12,800,600);
    c43->SetFillColor(0);
    c43->cd();
    hf->Draw();
    gsys->Draw("lfsame");
    guncorr->Draw("e2same");
    gstat->Draw("pzsame");

    Int_t nslots = 1;
    Float_t legslot = 0.07;
    TLegend *l1=new TLegend(0.19, 0.88-nslots*legslot, 0.48, 0.88);
    l1->SetBorderSize(0);
    l1->AddEntry(gstat,Form("%s %3.2f TeV", system.Data(), energy),"p");
    l1->Draw();
    f->cd();
    gstat->Write();
    gsys->Write();      
  }
    
  if (outputOpt.Contains("histo")) {
    foutname = Form("%s_%s%3.2fTeV_histo.root", file.Data(), system.Data(), energy);
    f = new TFile(foutname.Data(),"RECREATE"); 

    const int nbins = nlines+1;
    Double_t cbins[nbins];
    for (Int_t ii=0; ii<nbins; ii++) {if (ii<nlines) cbins[ii] = cl[ii]; else cbins[ii] = cu[ii-1];} 
    
    //create histos
    TH1D *hstat = new TH1D("stat", "stat", nlines, cbins);
    hstat->SetMarkerStyle(markE[icolor]);
    hstat->SetMarkerColor(colE[icolor]);
    hstat->SetMarkerSize(1.4);
    hstat->SetFillStyle(0);
    hstat->SetFillColorAlpha(colF[icolor], 0.4);
    hstat->SetLineWidth(1);
    hstat->SetLineColor(colE[icolor]);

    TH1D *hsys = new TH1D("sys", "sys", nlines, cbins);
    hsys->SetMarkerStyle(markE[icolor]);
    hsys->SetMarkerColor(colE[icolor]);
    hsys->SetMarkerSize(1.4);
    hsys->SetFillStyle(0);
    hsys->SetFillColorAlpha(colF[icolor], 0.4);
    hsys->SetLineWidth(1);
    hsys->SetLineColor(colE[icolor]);
	
    TH1D *huncorr = new TH1D("uncorr", "uncorr", nlines, cbins);
    huncorr->SetMarkerStyle(markE[icolor]);
    huncorr->SetMarkerColor(colE[icolor]);
    huncorr->SetMarkerSize(1.4);
    huncorr->SetFillStyle(3001);
    huncorr->SetFillColorAlpha(colF[icolor], 0.4);
    huncorr->SetLineWidth(1);
    huncorr->SetLineColor(colE[icolor]);
	
    for (int i = 0; i<nlines; i++) {
      int ibin = hstat->GetXaxis()->FindBin(ccentral[i]);
      hstat->SetBinContent(ibin, y[i]);
      hstat->SetBinError(ibin, uy[i]);
      hsys->SetBinContent(ibin, sy[i]);
      hsys->SetBinError(ibin, sy[i]);
    }
	
    //Draw
    gStyle->SetOptStat(0);
    TCanvas* c44=new TCanvas("c44","",12,12,800,600);
    c44->SetFillColor(0);
    c44->cd();
    hsys->Draw("E2");
    hstat->Draw("same");
	
    Int_t nslots = 1;
    Float_t legslot = 0.07;
    TLegend *l1=new TLegend(0.19, 0.88-nslots*legslot, 0.48, 0.88);
    l1->SetBorderSize(0);
    l1->AddEntry(hstat,Form("%s %3.2f TeV", system.Data(), energy),"p");
    l1->Draw();
    f->cd();
    hstat->Write();
    hsys->Write();      
    huncorr->Write();      
  }
      
  //close file
  in.close();
  return;
      

}

 
