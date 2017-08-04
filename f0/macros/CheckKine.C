void CheckKine(TString fname = "galice.root",
	       Int_t motherPDG = 321, //321=kaon
	       Int_t firstDaughterPDG = 13, //13=muon
	       Int_t secondDaughterPDG = -1,
	       Bool_t verbose = 0) 
{
  const Long_t f0PDG =  113;// 9010221;
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
	
  Double_t Pt2, OpeningAngle, OpeningAngle2, OpeningAngle3;

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
  TH2F *histoOpeningAngle = new TH2F("histoOpeningAngle","Daughters opening angle (#alpha); p_{T, gen} (GeV/#it{c}); #alpha (rad)", 1000, 0., 10.0, 400, 0., 4.0);
  TH2F *histoOpeningAngle2 = new TH2F("histoOpeningAngle2","Daughters opening angle (#alpha); #it{y}_{gen}; #alpha (rad)", 200, -1., 1.0, 400, 0., 4.0);
  TH2F *histoOpeningAngle3 = new TH2F("histoOpeningAngle3","Daughters opening angle (#alpha); #it{#phi}_{gen}; #alpha (rad)", 700, 0., 7.0, 400, 0., 4.0);
  TH2F *histof0pipi = new TH2F("histof0pipi","f_{0}(980) #rightarrow #pi^{+}#pi^{-}; p_{T, 1} (GeV/#it{c}); p_{T, 2} (GeV/#it{c})", 1000, 0., 10., 1000, 0., 10.);
  TH2F *histof0kk = new TH2F("histof0kk","f_{0}(980) #rightarrow K^{+}K^{-}; p_{T, 1} (GeV/#it{c}); p_{T, 2} (GeV/#it{c})", 750, 0., 7.5, 750, 0., 7.5);

  gStyle->SetOptStat(111111); //bin overflow	
  outKine<<"Table of PDG Codes"<<endl;
  
  /* CHECK ESDs for daughters of a decayed f0 */
  /* open esds */
  TFile *fesd = TFile::Open(Form("AliESDs.root"));
  if (!fesd || !fesd->IsOpen()) continue;
  TTree *esdT = (TTree *)fesd->Get("esdTree");
  
  AliESDEvent *esd = new AliESDEvent();
  esd->ReadFromTree(esdT);
  for (Int_t iev = 0; iev < esdT->GetEntries(); iev++) {
    esdT->GetEvent(iev);
    if (verbose) Printf("########################### Event %i - Ntracks = %i", iev, esd->GetNumberOfTracks());

    for (Int_t it = 0; it < esd->GetNumberOfTracks(); it++) {
      /* get track */
      AliESDtrack *track = (AliESDtrack*) esd->GetTrack(it);
      if (!track) continue;
      Int_t label = track->GetLabel();
      Float_t recPt = track->Pt();
	
      //fill here histos for pt, y and phi for all tracks
	
      //look for the pions from f0 decay from kinematics
      // if (track->GetLabel() == pione nello stack) {
      // }
      //fill here histos for pt, y and phi for the tracks of pions from f0
    }   
  }


  /* **************** */
  /* CHECK KINEMATICS */
  /* **************** */
  /* loop over events */
  for (Int_t iev = 0; iev < rl->GetNumberOfEvents(); iev++) {
    outKine<<"Event #: "<<iev<<endl;
    outKine<<"PDG Code"<<"\t\t"<<"Particle label"<<"\t\t"<<"Mother PDG Code"<<"\t\t"<<"Mother label"<<"\t\t"<<"First Daughter"<<"\t\t"<<"Last daughter"<<endl;
    /* get event */
    rl->GetEvent(iev);
    if (verbose) Printf("########################### Event %i", iev);
    
    /* loop over particles */
    AliStack *stack = rl->Stack();

    if (verbose) Printf("KINEMATICS: N tracks = %i   -- N primary = %i    -- N transported = %i",
			stack->GetNtrack(), stack->GetNprimary(), stack->GetNtransported());

    /* look at the stack searching for f0 among primary particles, then 
       searching for the daughters and filling histograms */
    for (Int_t ipart = 0; ipart < stack->GetNprimary(); ipart++) {
      /* get particle */
      TParticle *particle = stack->Particle(ipart);
      if (!particle) continue;
      if (particle->GetPdgCode() == f0PDG) {
	if (verbose) {
	  Printf("f0 found in stack entry %i ::::: Unique ID = %i - Status = %i", ipart, particle->GetUniqueID(), particle->GetStatusCode());
	  //Add check on the daughters
	  //...
	}
      }
    }

    /* look at the full stack for the daughters of an f0 and filling histograms */    
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
      Int_t firstDaugL = particle->GetFirstDaughter();
      if (firstDaugL <= 0) continue;
      TParticle *firstDaug = stack->Particle(firstDaugL);

      //Get second daughter
      Int_t secondDaugL = particle->GetDaughter(1);
      if (secondDaugL <= 0) continue;
      TParticle *secondDaug = stack->Particle(secondDaugL);

      //Get last daughter
      Int_t lastDaugL = particle->GetLastDaughter();
      if (lastDaugL <= 0) continue;
      TParticle *lastDaug = stack->Particle(lastDaugL);

      if (firstDaug) firstDaugPt = firstDaug->Pt();
      if (secondDaug) secondDaugPt = secondDaug->Pt();

      histoPDGCode->Fill(TMath::Abs(particle->GetPdgCode())); /* histogram of the PDG codes */

      /*selecting f0 with PDG code*/
      if (mom->GetPdgCode()==f0PDG) {		  
      
	outKine<<particle->GetPdgCode()<<"\t\t\t"<<particle->GetName()<<"\t\t\t"<<mom->GetPdgCode()<<"\t\t\t"<<mom->GetName()<<"\t\t\t"<<firstDaug->GetPdgCode()<<"\t\t\t"<<lastDaug->GetPdgCode()<<endl;
	
	/*histograms for pt and y of f0(980)*/      
	histof0pt->Fill(mom->Pt());
	histof0y->Fill(mom->Y());
	histof0mother->Fill(mom->GetMother(0));

	//Get f0 first daughter
	Int_t f0firstDaugL = mom->GetFirstDaughter();
	if (f0firstDaugL <= 0) continue;
	TParticle *f0firstDaug = stack->Particle(f0firstDaugL);

	//Get f0 second daughter
	Int_t f0secondDaugL = mom->GetDaughter(1);
	if (f0secondDaugL <= 0) continue;
	TParticle *f0secondDaug = stack->Particle(f0secondDaugL);

	//Get f0 last daughter
	Int_t f0lastDaugL = mom->GetLastDaughter();
	if (f0lastDaugL <= 0) continue;
	TParticle *f0lastDaug = stack->Particle(f0lastDaugL);

	if (!f0firstDaug || !f0secondDaug) continue;
	Float_t f0firstDaugPt = f0firstDaug->Pt();
	Float_t f0secondDaugPt = f0secondDaug->Pt();

	
	/* 2D histograms for pt of first and second daughter of the f0(980)—> ππ */
	if (TMath::Abs(f0firstDaug->GetPdgCode()) == pionPDG ||
	    TMath::Abs(f0secondDaug->GetPdgCode()) == pionPDG) {
	  histof0pipi->Fill(f0firstDaugPt, f0secondDaugPt);
	  OpeningAngle = pion1V.Angle(pion2V.Vect()); //Angle between products of decay
	}	  
	  
	/* 2D histograms for pt of first and second daughter of the f0(980)—> KK   */
	if(TMath::Abs(f0firstDaug->GetPdgCode()) == kaonPDG &&
	   TMath::Abs(f0secondDaug->GetPdgCode()) == kaonPDG) {
	  histof0kk->Fill(f0firstDaugPt, f0secondDaugPt);
	}
      	
	/* Use the energy-momentum quadrivector to check the method TParticle::Pt() */
	mom->Momentum(momV);
	f0firstDaug->Momentum(pion1V);
	f0secondDaug->Momentum(pion2V);
	Pt2=TMath::Sqrt((momV.Px())**2+(momV.Py())**2);
	histof0pt2->Fill(Pt2);

	/* pT, eta and phi of the generated mother vs daughters opening angle */
	histoOpeningAngle->Fill(Pt2, OpeningAngle);
	histoOpeningAngle2->Fill(mom->Y(), OpeningAngle);
	histoOpeningAngle3->Fill(mom->Phi(), OpeningAngle);	
	
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
