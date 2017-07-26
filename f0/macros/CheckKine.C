void CheckKine(TString fname = "galice.root",
	       Int_t motherPDG = 321, //321=kaon
		   Int_t motherf0 = 9010221, //9010221=f0
	       Int_t daughterPDG = 13, //13=muon
		   Int_t daughterpi = 211) //211=pion
{

  /* check if AliEn connection is needed */
  if (fname.BeginsWith("alien"))
    TGrid::Connect("alien");
  
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
  TVector3 muonV, kaonV, pion1V, pion2V, momV;

  /* output */
  TFile *fout = TFile::Open("OutKine.root", "RECREATE");  
  ofstream outKine;
  outKine.open("outKine.txt");
  TH1 *hqT = new TH1F("hqT", "", 1000, 0., 1.);
  TH1F *histoPDGCode = new TH1F("histoPDGCode", "histoPDGCode", 300, -1000.0, 3000.0);
  TH1F *histof0pt = new TH1F("histof0pt", "histof0pt", 150, 0., 15.0);
  TH1F *histof0y = new TH1F("histof0y", "histof0y", 100, -1.0, 1.0);
  TH2F *histof0pipi = new TH2F("histof0pipi","histof0pipi", 150, 0., 0.2, 150, 0., 0.2);
  TH2F *histof0kk = new TH2F("histof0kk","histof0kk", 150, 0., 0.2, 150, 0., 0.2);
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
	  Int_t momL = particle->GetMother(0);
	  if (momL <= 0) continue;
	  TParticle *mom = stack->Particle(momL);
	  Int_t firstdaugL = particle->GetFirstDaughter();
	  if (firstdaugL <= 0) continue;
	  TParticle *firstdaug = stack->Particle(firstdaugL);
	  Int_t seconddaugL = particle->GetDaughter(1);
	  if (seconddaugL <= 0) continue;
	  TParticle *seconddaug = stack->Particle(seconddaugL);
	  Int_t lastdaugL = particle->GetLastDaughter();
	  if (lastdaugL <= 0) continue;
	  TParticle *lastdaug = stack->Particle(lastdaugL);
	  
	  outKine<<particle->GetPdgCode()<<"\t\t\t"<<particle->GetName()<<"\t\t\t"<<mom->GetPdgCode()<<"\t\t\t"<<mom->GetName()<<"\t\t\t"<<firstdaug->GetPdgCode()<<"\t\t\t"<<lastdaug->GetPdgCode()<<endl;
      
	  histoPDGCode->Fill(particle->GetPdgCode()); /* histogram of the PDG codes */
    
	  /*histograms for pt and y of f0(980), by selecting f0 with PDG code*/
	  if (mom->GetPdgCode() == motherf0){
		  
		  histof0pt->Fill(mom->Pt());
		  histof0y->Fill(mom->Y());
      	
	  }
   
	 
	  
	  
	 /* 2D histograms for pt of first and second daughter of the f0(980)—> ππ */
	  if(TMath::Abs(firstdaug->GetPdgCode()) != daughterpi && TMath::Abs(seconddaug->GetPdgCode()) != daughterpi) continue;
      if (TMath::Abs(mom->GetPdgCode()) != motherf0) continue;
	 
	  Double_t x = firstdaug->Pt();
	  Double_t y = seconddaug->Pt();
	  

      histof0pipi->Fill(x,y);
	  
	  
	  
	 /* 2D histograms for pt of first and second daughter of the f0(980)—> KK   */
	  if(TMath::Abs(firstdaug->GetPdgCode()) != motherPDG && TMath::Abs(seconddaug->GetPdgCode()) != motherPDG) continue;
      if (TMath::Abs(mom->GetPdgCode()) != motherf0) continue;
	  Double_t a = firstdaug->Pt();
	  Double_t b = seconddaug->Pt();
	  histof0kk->Fill(a,b);   
	  
	

	      
      /* find requested mother-daughter relationship */
      if (TMath::Abs(particle->GetPdgCode()) != daughterPDG) continue;
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
  fout->cd();
  hqT->Write();
  histoPDGCode->Write();
  histof0pt->Write();
  histof0y->Write();
  histof0pipi->Write();
  histof0kk->Write();
  fout->Close();

  TCanvas * c1 = new TCanvas("c1","c1", 1200, 800);
  c1->Divide(3,2);
  c1->cd(1); histoPDGCode->Draw(); 
  c1->cd(2); histof0pt->Draw("hist");
  c1->cd(3); histof0y->Draw("hist");
  c1->cd(4); histof0pipi->Draw("colz");
  c1->cd(5); histof0kk->Draw("colz");

  c1->Print("checkKineOut.pdf");
  return;
}
