void CheckKine(TString fname = "galice.root",
	       Int_t motherPDG = 321, //321=kaon
	       Int_t firstDaughterPDG = 13, //13=muon
	       Int_t secondDaughterPDG = -1) 
{
  const Long_t f0PDG = 9010221;
  const Int_t pionPDG = 211; //charged pion
  const Int_t kaonPDG = 321; //charged kaon

  /* check if AliEn connection is needed */
  if (fname.BeginsWith("alien"))
    TGrid::Connect("alien://");
  
  /* load libpythia */
  gSystem->Load("libpythia6.so");
  
  AliRunLoader *rl = AliRunLoader::Open(fname.Data());
  if (!rl) {
    printf("cannot open the run loader \n");
    return;
  }
  rl->LoadgAlice();
  if (!gAlice) {
    printf("cannot find galice object \n");
  }
  rl->LoadKinematics();

  /* calculation tools */
  TVector3 muonV, kaonV;
  TLorentzVector pion1V, pion2V, momV;
	
  Double_t pt2, OpeningAngle, OpeningAngle2, OpeningAngle3;

  /* output */
  TFile *fout = TFile::Open("OutKine.root", "RECREATE");  
  ofstream outKine;
  outKine.open("outKine.txt");
  TH1F *hqT = new TH1F("hqT", "", 1000, 0., 1.);
  TH1F *histoPDGCode = new TH1F("histoPDGCode", "histoPDGCode; PDG code; entries", 10001, -1.0, 10000.0);
  TH1F *histof0pt = new TH1F("histof0pt", "generated f_{0}(980) p_{T}; p_{T, gen} (GeV/#it{c}); entries", 100, 0., 10.0);
  TH1F *histof0pt2 = new TH1F("histof0pt2", "generated f_{0}(980) p_{T}; p_{T, gen} (GeV/#it{c}); entries", 100, 0., 10.0);	
  TH1F *histof0y = new TH1F("histof0y", "generated f_{0}(980) rapidity; #it{y}_{gen}; entries", 100, -1.0, 1.0);
  TH1F *histof0mother = new TH1F("histof0mother", "label of f0 mother; label of f0 mother; entries", 102, -2.0, 100);	
  TH2F *histoOpeningAngle = new TH2F("histoOpeningAngle","Daughters opening angle (#alpha); #alpha; p_{T, gen} (GeV/#it{c})", 400, 0., 4.0, 1000, 0., 10.0);
  TH2F *histoOpeningAngle2 = new TH2F("histoOpeningAngle2","Daughters opening angle (#alpha); #alpha; #it{y}_{gen}", 400, 0., 4.0, 200, -1., 1.0);
  TH2F *histoOpeningAngle3 = new TH2F("histoOpeningAngle3","Daughters opening angle (#alpha); #alpha; #it{#phi}_{gen}", 400, 0., 4.0, 700, 0., 7.0);
  TH2F *histof0pipi = new TH2F("histof0pipi","f_{0}(980) #rightarrow #pi^{+}#pi^{-}; p_{T, 1} (GeV/#it{c}); p_{T, 2} (GeV/#it{c})", 750, 0., 7.5, 750, 0., 7.5);
  TH2F *histof0kk = new TH2F("histof0kk","f_{0}(980) #rightarrow K^{+}K^{-}; p_{T, 1} (GeV/#it{c}); p_{T, 2} (GeV/#it{c})", 750, 0., 7.5, 750, 0., 7.5);

  gStyle->SetOptStat(111111); //bin overflow	
  outKine<<"Table of PDG Codes"<<endl;
 
  /* loop over events */
  for (Int_t iev = 0; iev < rl->GetNumberOfEvents(); iev++) {
    outKine<<"Event #: "<<iev<<endl;
    outKine<<"PDG Code"<<"\t\t"<<"Particle label"<<"\t\t"<<"Mother PDG Code"<<"\t\t"<<"Mother label"<<"\t\t"<<"First Daughter"<<"\t\t"<<"Last daughter"<<endl;
    /* get event */
    rl->GetEvent(iev);

    /* loop over particles */
    AliStack *stack = rl->Stack();
    for (Int_t ipart = 0; ipart < stack->GetNtrack(); ipart++) {
      
      /* get particle */
      TParticle *particle = stack->Particle(ipart);
      if (!particle) continue;

      Float_t firstDaugPt = 0.0; 
      Float_t secondDaugPt = 0.0;

      //Get mother
      Int_t momL = particle->GetMother(0);
      if (momL <= 0) continue;
      TParticle *mom = stack->Particle(momL);

      //Get first daughter
      Int_t firstdaugL = particle->GetFirstDaughter();
      if (firstdaugL <= 0) continue;
      TParticle *firstdaug = stack->Particle(firstdaugL);

      //Get second daughter
      Int_t seconddaugL = particle->GetDaughter(1);
      if (seconddaugL <= 0) continue;
      TParticle *seconddaug = stack->Particle(seconddaugL);

      //Get last daughter
      Int_t lastdaugL = particle->GetLastDaughter();
      if (lastdaugL <= 0) continue;
      TParticle *lastdaug = stack->Particle(lastdaugL);

      if (firstdaug) firstDaugPt = firstdaug->Pt();
      if (seconddaug) secondDaugPt = seconddaug->Pt();

      histoPDGCode->Fill(TMath::Abs(particle->GetPdgCode())); /* histogram of the PDG codes */

      /*selecting f0 with PDG code*/
      if (mom->GetPdgCode()==f0PDG){		  
      
	outKine<<particle->GetPdgCode()<<"\t\t\t"<<particle->GetName()<<"\t\t\t"<<mom->GetPdgCode()<<"\t\t\t"<<mom->GetName()<<"\t\t\t"<<firstdaug->GetPdgCode()<<"\t\t\t"<<lastdaug->GetPdgCode()<<endl;
	
	/*histograms for pt and y of f0(980)*/      
	histof0pt->Fill(mom->Pt());
	histof0y->Fill(mom->Y());
	histof0mother->Fill(mom->GetMother(0));

	//Get f0 first daughter
	Int_t f0firstdaugL = mom->GetFirstDaughter();
	if (f0firstdaugL <= 0) continue;
	TParticle *f0firstdaug = stack->Particle(f0firstdaugL);

	//Get f0 second daughter
	Int_t f0seconddaugL = mom->GetDaughter(1);
	if (f0seconddaugL <= 0) continue;
	TParticle *f0seconddaug = stack->Particle(f0seconddaugL);

	//Get f0 last daughter
	Int_t f0lastdaugL = mom->GetLastDaughter();
	if (f0lastdaugL <= 0) continue;
	TParticle *f0lastdaug = stack->Particle(f0lastdaugL);

	if (!f0firstdaug || !f0seconddaug) continue;
	Float_t f0firstDaugPt = f0firstdaug->Pt();
	Float_t f0secondDaugPt = f0seconddaug->Pt();

	
	/* 2D histograms for pt of first and second daughter of the f0(980)—> ππ */
	if (TMath::Abs(f0firstdaug->GetPdgCode()) == pionPDG ||
	    TMath::Abs(f0seconddaug->GetPdgCode()) == pionPDG) {
	  histof0pipi->Fill(f0firstDaugPt, f0secondDaugPt);
          OpeningAngle = pion1V.Angle(pion2V.Vect()); //Angle between products of decay
	}	  
	  
	/* 2D histograms for pt of first and second daughter of the f0(980)—> KK   */
	if(TMath::Abs(f0firstdaug->GetPdgCode()) == kaonPDG &&
	   TMath::Abs(f0seconddaug->GetPdgCode()) == kaonPDG) {
	  histof0kk->Fill(f0firstDaugPt, f0secondDaugPt);
	}
      	
	/* Use the energy-momentum quadrivector to check the method TParticle::Pt() */
	mom->Momentum(momV);
	f0firstdaug->Momentum(pion1V);
	f0seconddaug->Momentum(pion2V);
	pt2=TMath::Sqrt((momV.Px())**2+(momV.Py())**2);
	histof0pt2->Fill(pt2);

	/* information on the daughters opening angle vs pT, eta and phi of the generated mother */
	histoOpeningAngle->Fill(OpeningAngle, pt2);
	histoOpeningAngle2->Fill(OpeningAngle, mom->Y());
	histoOpeningAngle3->Fill(OpeningAngle, mom->Phi());
	
	
      } //end checks on f0 PDG 

	      
      /* K --> µ + nu: find requested mother-daughter relationship */
      if (TMath::Abs(particle->GetPdgCode()) != firstDaughterPDG) continue;
      Int_t motherL = particle->GetMother(0);
      if (motherL <= 0) continue;
      if (motherL >= stack->GetNtrack()) {
        printf(">>> WEIRD: mother label is larger than stack size: %d \n", motherL);
        continue;
      }
      TParticle *mother = stack->Particle(motherL);
      if (!mother) continue;
      if (TMath::Abs(mother->GetPdgCode()) != motherPDG) continue;

      /* compute daughter qT */
      muonV.SetXYZ(particle->Px(), particle->Py(), particle->Pz());
      kaonV.SetXYZ(mother->Px(), mother->Py(), mother->Pz());
      Double_t qT = muonV.Perp(kaonV);
      hqT->Fill(qT);	  	  

    } /* loop over particles */

  } /* loop over events */

  /* write output */
  TCanvas * c1 = new TCanvas("c1","c1", 1200, 800);
  c1->Divide(3,2);
  c1->cd(1); histoPDGCode->Draw(); gPad->SetLogy();
  c1->cd(2); histof0pt->Draw("hist");
  c1->cd(3); histof0y->Draw("hist");
  c1->cd(4); histof0mother->Draw("hist");		
  c1->cd(5); histof0pipi->Draw("colz");  gPad->SetLogy();  gPad->SetLogx();
  c1->cd(6); histof0kk->Draw("colz");  gPad->SetLogy();  gPad->SetLogx();
	
  TCanvas * c2 = new TCanvas("c2","c2", 1200, 800);
  c2->Divide(3,1);
  c2->cd(1); histoOpeningAngle->Draw("colz");
  c2->cd(2); histoOpeningAngle2->Draw("colz");
  c2->cd(3); histoOpeningAngle3->Draw("colz");
	
  c1->Print("checkKineOut.pdf");
  c2->Print("checkKineOut2.pdf");	
	
  fout->cd();
  hqT->Write();
  histoPDGCode->Write();
  histof0pt->Write();
  histof0pt2->Write();
  histof0y->Write();
  histof0mother->Write();
  histof0pipi->Write();
  histoOpeningAngle->Write();
  histoOpeningAngle2->Write();
  histoOpeningAngle3->Write();
  histof0kk->Write();
  fout->Close();

 
  return;
}
