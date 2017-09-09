///////////////////////////////////////////////////////////////////////////////////////////////
//       exercise: check how Root treats the rebinning of multi-dimensional histograms       //
///////////////////////////////////////////////////////////////////////////////////////////////


#include "TH2.h"
#include "TCanvas.h"

void exercise(){

Int_t a=0;
TH2D * histo = new TH2D("histo", "histo", 10, 0., 10., 10, 0., 10.);
TH2D * histo2 = new TH2D("histo2", "histo2", 10, 0., 10., 10, 0., 10.);
for(Int_t i=0; i<10; i++){
  a=10*i+1;
  for(Int_t j=0; j<10; j++){
    histo->Fill(i, j, a);
    histo2->Fill(i, j, a);
    a++;
  }
}

histo2->RebinX(2);
//histo2->RebinY(2);

TCanvas* c1 = new TCanvas("c1", "c1", 800, 400);
gStyle->SetOptStat(111111);
c1->Divide(2,1);
c1->cd(1);
histo->Draw("lego");
c1->cd(2);
histo2->Draw("lego");

}
