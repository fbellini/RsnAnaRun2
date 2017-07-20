void CheckKine(TString fname = "galice.root",
	       Int_t motherPDG = 321,
	       Int_t daughterPDG = 13)
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
  TVector3 muonV, kaonV;

  /* output */
  TFile *fout = TFile::Open("OutKine.root", "RECREATE");
  TH1 *hqT = new TH1F("hqT", "", 1000, 0., 1.);

  /* loop over events */
  for (Int_t iev = 0; iev < rl->GetNumberOfEvents(); iev++) {

    /* get event */
    rl->GetEvent(iev);

    /* loop over particles */
    AliStack *stack = rl->Stack();
    for (Int_t ipart = 0; ipart < stack->GetNtrack(); ipart++) {
      
      /* get particle */
      TParticle *particle = stack->Particle(ipart);
      if (!particle) continue;
      
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
  fout->Close();

return;
}
