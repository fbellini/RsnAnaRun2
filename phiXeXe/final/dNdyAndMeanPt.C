void dNdyAndMeanPt()
{
  //get dN/dy and mean pt with syst sude to extrapolation
  const Short_t nfunc = 5;
  const Short_t ncent = 5;
  TFile * file[nfunc];
  file[0] = TFile::Open("FITSPECTRUM_bgbw_0.5-10.0_16sep20.root");
  file[1] = TFile::Open("FITSPECTRUM_boltz_0.5-10.0_16sep20.root");
  file[2] = TFile::Open("FITSPECTRUM_bose_0.5-10.0_16sep20.root");
  file[3] = TFile::Open("FITSPECTRUM_levy_0.5-10.0_16sep20.root");
  file[4] = TFile::Open("FITSPECTRUM_mtexp_0.5-10.0_16sep20.root");

  TH1F * fitResults[ncent][nfunc];

  for (int ic = 0; ic<ncent; ic++){

    Printf("--------------------- CENTRALITY %i", ic);
    Float_t yieldMax = 0.0;
    Float_t yieldMin = 1e10;
    Float_t meanptMax = 0.0;
    Float_t meanptMin = 1e10;
    Float_t yield[nfunc];
    Float_t meanpt[nfunc];

    for (int ifu = 0; ifu<nfunc; ifu++){
      fitResults[ic][ifu] = (TH1F*) file[ifu]->Get(Form("fitResults%i",ic));
      if (!fitResults[ic][ifu]) Printf("Invalid histo cent %i func %i", ic, ifu);
      yield[ifu] = fitResults[ic][ifu]->GetBinContent(1);
      meanpt[ifu] = fitResults[ic][ifu]->GetBinContent(5);
    }
    
    for (int ifu = 0; ifu < nfunc; ifu++){
      if (yield[ifu] > yieldMax) yieldMax = yield[ifu];
      else if (yield[ifu] < yieldMin) yieldMin = yield[ifu];

      if (meanpt[ifu] > meanptMax) meanptMax = meanpt[ifu];
      else if (meanpt[ifu] < meanptMin) meanptMin = meanpt[ifu]; 
    }

    Float_t yieldExtrapErr = 0.5*(yieldMax-yieldMin);
    Float_t meanptExtrapErr = 0.5*(meanptMax-meanptMin);

    Float_t totSys = TMath::Sqrt( TMath::Power(fitResults[ic][0]->GetBinContent(3), 2.0) + TMath::Power(yieldExtrapErr, 2.0) );
    Float_t totSysMpt = TMath::Sqrt( TMath::Power(fitResults[ic][0]->GetBinContent(7), 2.0) + TMath::Power(meanptExtrapErr, 2.0) );
    Printf("dN/dy = %8.4f +- %8.4f (stat) +- %8.4f (sys) +- %8.4f (extrap, %3.1f%%) -- tot sys = %8.4f", 
    	   fitResults[ic][0]->GetBinContent(1), fitResults[ic][0]->GetBinContent(2), fitResults[ic][0]->GetBinContent(3), yieldExtrapErr, yieldExtrapErr*100./fitResults[ic][0]->GetBinContent(1), totSys);
    
    Printf(" <pT> = %8.4f +- %8.4f (stat) +- %8.4f (sys) +- %8.4f (extrap, %3.1f%%) -- tot sys = %8.4f", fitResults[ic][0]->GetBinContent(5), fitResults[ic][0]->GetBinContent(6), fitResults[ic][0]->GetBinContent(7), meanptExtrapErr, meanptExtrapErr*100./fitResults[ic][0]->GetBinContent(5), totSysMpt);
  }
  return;

}
