/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/


/////////////////////////////////////////////////////////////////////////
//                                                                     //
//      Task for f0 background (MC efficiencies and template)          //
//      analysis in pp collisions                                      //
//      27.11.2017                                                     //
//	                                                                   //
/////////////////////////////////////////////////////////////////////////


#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskF0Bkg.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliMCParticle.h"
#include "AliESDpid.h"
#include "AliESDtrackCuts.h"
#include "AliGenEventHeader.h"
#include "AliVHeader.h"
#include "AliAnalysisCuts.h"
#include "AliAnalysisFilter.h"
#include "AliAnalysisVertexingHF.h"
#include "AliPIDResponse.h"
#include "AliTPCPIDResponse.h"
#include "AliTOFPIDResponse.h"
#include "AliESDpid.h"
#include "AliTOFPIDParams.h"


class AliAnalysisTaskF0Bkg;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskF0Bkg) // classimp: necessary for root

//f0PDG, omegaPDG, rhoPDG, etaPDG, etaPrimePDG, f1PDG, f2PDG, kstarPDG, k0s, phiPDG, pionPDG};
const ULong_t AliAnalysisTaskF0Bkg::fPdgArray[]={9010221, 223, 113, 221, 331, 20223, 225, 313, 310, 333, 211};

const Char_t AliAnalysisTaskF0Bkg::fParticleName[][6]={"f0","omega","rho","eta","etaPr","f1","f2","kStar","k0s", "phi", "pion"};


AliAnalysisTaskF0Bkg::AliAnalysisTaskF0Bkg() : AliAnalysisTaskSE(),
fESD(0), fOutputList(0), fNEvents(0), fMCEvent(0), fMCStack(0), fTrackFilter(0x0), fPIDResponse(0), fHistPID1tpc(0), fHistPID2tpc(0), fHistPID1tof(0), fHistPID2tof(0), fNoSelTracksEtaPt(0), fAcceptedGeomAccEtaPt(0), fAcceptedQualityEtaPt(0), fAcceptedTOFMatchEtaPt(0), fAcceptedTracksEtaPt(0)
{
  // default constructor, don't allocate memory here!
  // this is used by root for IO purposes, it needs to remain empty

  for(Int_t i=0; i<10; i++){
    fGenMassVsPt[i]=0x0;
    fRecoMassVsPt[i]=0x0;
    fGenYVsPt[i]=0x0;
    fRecoYVsPt[i]=0x0;
    fGenEtaVsPt[i]=0x0;
    fRecoEtaVsPt[i]=0x0;
  }
}
//_____________________________________________________________________________
AliAnalysisTaskF0Bkg::AliAnalysisTaskF0Bkg(const char* name) : AliAnalysisTaskSE(name),
fESD(0), fOutputList(0), fNEvents(0), fMCEvent(0), fMCStack(0), fTrackFilter(0x0), fPIDResponse(0), fHistPID1tpc(0), fHistPID2tpc(0), fHistPID1tof(0), fHistPID2tof(0), fNoSelTracksEtaPt(0), fAcceptedGeomAccEtaPt(0), fAcceptedQualityEtaPt(0), fAcceptedTOFMatchEtaPt(0), fAcceptedTracksEtaPt(0)
{
  // constructor

  for(Int_t i=0; i<10; i++){
    fGenMassVsPt[i]=0x0;
    fRecoMassVsPt[i]=0x0;
    fGenYVsPt[i]=0x0;
    fRecoYVsPt[i]=0x0;
    fGenEtaVsPt[i]=0x0;
    fRecoEtaVsPt[i]=0x0;
  }

  DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
  // this chain is created by the analysis manager, so no need to worry about it,
  // it does its work automatically
  DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms
  // you can add more output objects by calling DefineOutput(2, classname::Class())
  // if you add more output objects, make sure to call PostData for all of them, and to
  // make changes to your AddTask macro!
}
//_____________________________________________________________________________
AliAnalysisTaskF0Bkg::~AliAnalysisTaskF0Bkg()
{
  // destructor
  if (fPIDResponse) delete fPIDResponse;
  if (fESD) delete fESD;
  if (fMCEvent) delete fMCEvent;
  if (fMCStack) delete fMCStack;
  if (fTrackFilter) delete fTrackFilter;
  /*if (fHlistPID){
   delete fHlistPID;
   fHlistPID = 0;
 }*/
 if (fNEvents){
  delete fNEvents;
  fNEvents = 0;
 }
 if (fHistPID1tpc){
  delete fHistPID1tpc;
  fHistPID1tpc = 0;
 }
 if (fHistPID2tpc){
  delete fHistPID2tpc;
  fHistPID2tpc = 0;
 }
 if (fHistPID1tof){
  delete fHistPID1tpc;
  fHistPID1tpc = 0;
 }
 if (fHistPID2tof){
  delete fHistPID2tpc;
  fHistPID2tpc = 0;
 }
 if (fNoSelTracksEtaPt){
  delete fNoSelTracksEtaPt;
  fNoSelTracksEtaPt = 0;
 }
 if (fAcceptedGeomAccEtaPt){
  delete fAcceptedGeomAccEtaPt;
  fAcceptedGeomAccEtaPt = 0;
 }
 if (fAcceptedQualityEtaPt){
  delete fAcceptedQualityEtaPt;
  fAcceptedQualityEtaPt = 0;
 }
 if (fAcceptedTOFMatchEtaPt){
  delete fAcceptedTOFMatchEtaPt;
  fAcceptedTOFMatchEtaPt = 0;
 }
 if (fAcceptedTracksEtaPt){
  delete fAcceptedTracksEtaPt;
  fAcceptedTracksEtaPt = 0;
 }
 if(fOutputList) {
    delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskF0Bkg::UserCreateOutputObjects()
{
  // create output objects

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  if (!man)  AliFatal("Analysis manager needed");
  AliInputEventHandler *inputHandler=dynamic_cast<AliInputEventHandler*>(man->GetInputEventHandler());
  if (!inputHandler) AliFatal("Input handler needed");

  fPIDResponse = inputHandler->GetPIDResponse();
  if (!fPIDResponse) AliError("PID Response object was not created");

  fOutputList = new TList();          // this is a list which will contain all of your histograms
  // at the end of the analysis, the contents of this list are written
  // to the output file
  fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
  // if requested (dont worry about this now)

  fNEvents = new TH1I("fNEvents", "Event selection monitoring", 10, 0, 10);
  fNEvents->GetXaxis()->SetBinLabel(1, "accepted");
  fNEvents->GetXaxis()->SetBinLabel(2, "input");
  fNEvents->GetXaxis()->SetBinLabel(3, "no contributors");
  fNEvents->GetXaxis()->SetBinLabel(4, "no trk vtx");
  fNEvents->GetXaxis()->SetBinLabel(5, "no SPD vtx");
  fNEvents->GetXaxis()->SetBinLabel(6, "bad SPD resol");
  fNEvents->GetXaxis()->SetBinLabel(7, "no proximity");
  fNEvents->GetXaxis()->SetBinLabel(8, "z_{vtx} out range");
  fOutputList->Add(fNEvents);

  fHistPID1tpc = new TH2F("fHistPID1tpc", "N sigma vs pT; #it{p}_{T} (GeV/#it{c}); #it{N #sigma}", 220, 0., 11., 800, -4., 4.);
  fOutputList->Add(fHistPID1tpc);
  fHistPID2tpc = new TH2F("fHistPID2tpc", "N sigma vs pT; #it{p}_{T} (GeV/#it{c}); #it{N #sigma}", 220, 0., 11., 800, -4., 4.);
  fOutputList->Add(fHistPID2tpc);
  fHistPID1tof = new TH2F("fHistPID1tof", "N sigma vs pT; #it{p}_{T} (GeV/#it{c}); #it{N #sigma}", 220, 0., 11., 800, -4., 4.);
  fOutputList->Add(fHistPID1tof);
  fHistPID2tof = new TH2F("fHistPID2tof", "N sigma vs pT; #it{p}_{T} (GeV/#it{c}); #it{N #sigma}", 220, 0., 11., 800, -4., 4.);
  fOutputList->Add(fHistPID2tof);

  fNoSelTracksEtaPt = new TH2F("fNoSelTracksEtaPt", "#eta vs. pT - no selected tracks; #eta; #it{p}_{T} (GeV/#it{c})", 600, -3., 3., 220, 0., 11.);
  fOutputList->Add(fNoSelTracksEtaPt);

  fAcceptedGeomAccEtaPt = new TH2F("fAcceptedGeomAccEtaPt", "#eta vs. pT - accepted tracks (geom acceptance); #eta; #it{p}_{T} (GeV/#it{c})", 600, -3., 3., 220, 0., 11.);
  fOutputList->Add(fAcceptedGeomAccEtaPt);

  fAcceptedQualityEtaPt = new TH2F("fAcceptedQualityEtaPt", "#eta vs. pT - accepted tracks (geom acceptance + 2011 std cuts); #eta; #it{p}_{T} (GeV/#it{c})", 600, -3., 3., 220, 0., 11.);
  fOutputList->Add(fAcceptedQualityEtaPt);

  fAcceptedTOFMatchEtaPt = new TH2F("fAcceptedTOFMatchEtaPt", "#eta vs. pT - accepted tracks (geom acceptance + 2011 std cuts + TOF match); #eta; #it{p}_{T} (GeV/#it{c})", 600, -3., 3., 220, 0., 11.);
  fOutputList->Add(fAcceptedTOFMatchEtaPt);

  fAcceptedTracksEtaPt = new TH2F("fAcceptedTracksEtaPt", "#eta vs. pT - accepted tracks (geom acceptance + 2011 std cuts + TOF match + PID); #eta; #it{p}_{T} (GeV/#it{c})", 600, -3., 3., 220, 0., 11.);
  fOutputList->Add(fAcceptedTracksEtaPt);

  for(Int_t j=0; j<10; j++){
    fGenMassVsPt[j]= new TH2F(Form("fGenMassVsPt_%s",fParticleName[j]), Form("generated %s; M_{#pi#pi} (GeV/#it{c}^{2}); #it{p}_{T} (GeV/#it{c})", fParticleName[j]), 1200, 0.3, 1.5, 220, 0., 11.);
    fRecoMassVsPt[j] = new TH2F(Form("fRecoMassVsPt_%s",fParticleName[j]), Form("reconstructed %s; M_{#pi#pi} (GeV/#it{c}^{2}); #it{p}_{T} (GeV/#it{c})", fParticleName[j]), 1200, 0.3, 1.5, 220, 0., 11.);
    fGenYVsPt[j] = new TH2F(Form("fGenYVsPt_%s",fParticleName[j]), Form("generated %s rapidity; #it{p}_{T} (GeV/#it{c}); #it{y}", fParticleName[j]), 220, 0., 11., 14, -0.7, 0.7);
    fRecoYVsPt[j] = new TH2F(Form("fRecoYVsPt_%s",fParticleName[j]), Form("reconstructed %s rapidity; #it{p}_{T} (GeV/#it{c}); #it{y}", fParticleName[j]), 220, 0., 11., 14, -0.7, 0.7);
    fGenEtaVsPt[j] = new TH2F(Form("fGenEtaVsPt_%s",fParticleName[j]), Form("generated %s pseudorapidity; #it{p}_{T} (GeV/#it{c}); #it{#eta}", fParticleName[j]), 220, 0., 11., 20, -1., 1.);
    fRecoEtaVsPt[j] = new TH2F(Form("fRecoEtaVsPt_%s",fParticleName[j]), Form("reconstructed %s pseudorapidity; #it{p}_{T} (GeV/#it{c}); #it{#eta}", fParticleName[j]), 220, 0., 11., 20, -1., 1.);

    fOutputList->Add(fGenMassVsPt[j]);          // don't forget to add it to the list! the list will be written to file, so if you want
    fOutputList->Add(fRecoMassVsPt[j]);         // your histogram in the output file, add it to the list!
    fOutputList->Add(fGenYVsPt[j]);
    fOutputList->Add(fRecoYVsPt[j]);
    fOutputList->Add(fGenEtaVsPt[j]);
    fOutputList->Add(fRecoEtaVsPt[j]);
  }
  PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the
  // fOutputList object. the manager will in the end take care of writing your output to file
  // so it needs to know what's in the output
  //PostData(2, fHlistPID);
}
//_____________________________________________________________________________
void AliAnalysisTaskF0Bkg::UserExec(Option_t *)
{
  // user exec
  // this function is called once for each event
  // the manager will take care of reading the events from file, and with the static function InputEvent() you
  // have access to the current event.
  // once you return from the UserExec function, the manager will retrieve the next event from the chain

  //printf("------------------------------------------------------------\n");

  AliVEvent *event = InputEvent();
  if(!event) return;

  Bool_t isESD = kFALSE;
  if (dynamic_cast<AliESDEvent*>(event)) isESD = kTRUE;
  else if (!dynamic_cast<AliAODEvent*>(event))
  AliFatal("I don't find the AOD event nor the ESD one, aborting.");

  //apply vertex selection and fill histogram with number of accepted events
  if (!SelectVertex2015pp(event, kTRUE, kFALSE, kTRUE, kTRUE)) return;
  fNEvents->Fill(0);

  ///////////// ------ generated f0s ------ /////////////
  AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());

  if (!eventHandler) {
    AliError("Could not retrieve MC event handler");
    return;
  }

  fMCEvent = (AliMCEvent*) eventHandler->MCEvent();
  if (!fMCEvent) {
    AliError("Could not retrieve MC event");
    return;
  }

  fMCStack = fMCEvent->Stack();
  if (!fMCStack) {
    AliError("Could not retrieve MC Stack");
    return;
  }

  TParticle* trackMC = 0x0;
  TLorentzVector v1, v2, vMC, vSum;
  const Float_t piMass = 0.13957;  // GeV/c^2

  for (Int_t iStack = 0; iStack < fMCStack->GetNtrack(); iStack++) {

    trackMC = (TParticle*) fMCStack->Particle(iStack);
    ULong_t pdgTest=trackMC->GetPdgCode();
    if (TMath::Abs((trackMC->Eta()))>0.8) continue;
    if (TMath::Abs((trackMC->Pt()))<0.15) continue;

    for(Int_t index=0; index<10; index++){
      if (pdgTest==fPdgArray[index]){
        vMC.SetXYZM(trackMC->Px(), trackMC->Py(), trackMC->Pz(), piMass);
        if (fabs(vMC.Rapidity())>0.5) continue; //rapidity cut
        fGenMassVsPt[index]->Fill(trackMC->GetCalcMass(), trackMC->Pt());
        fGenYVsPt[index]->Fill(trackMC->Pt(), trackMC->Y());
        fGenEtaVsPt[index]->Fill(trackMC->Pt(), trackMC->Eta());
      } //endif
    } //pdg code loop
  } //stack loop

  ///////////// ------ reconstructed f0s ------ /////////////
  if(!fPIDResponse) Printf("::::: ERROR: no PID available.");

  Int_t iTracks = event->GetNumberOfTracks();           // see how many tracks there are in the event
  for(Int_t k1 = 0; k1 < iTracks; k1++) {                 // loop over all these tracks

    AliESDtrack* track1 = static_cast<AliESDtrack*>(event->GetTrack(k1));
    if(!track1) continue;

    fNoSelTracksEtaPt->Fill(track1->Eta(), track1->Pt());

    if (TMath::Abs((track1->Eta()))>0.8) continue;
    if (TMath::Abs((track1->Pt()))<0.15) continue;

    fAcceptedGeomAccEtaPt->Fill(track1->Eta(), track1->Pt());

    if (!fTrackFilter->IsSelected(track1)) continue;

    fAcceptedQualityEtaPt->Fill(track1->Eta(), track1->Pt());

    Int_t daug1Label = track1->GetLabel();
    if (daug1Label < 0) continue;
    TParticle*  daughter1 = fMCStack->Particle(daug1Label);
    Int_t mother1Label = daughter1->GetFirstMother();
    TParticle* mother1 = fMCStack->Particle(mother1Label);
    Long_t mother1PDG = mother1->GetPdgCode();
    Int_t daug1PDG = daughter1->GetPdgCode();
    //printf("%ld\n", daug1PDG);
    if (TMath::Abs(daug1PDG) != 211) continue;

    for(Int_t k2 = k1+1; k2 < iTracks; k2++){
      if (k2!=k1){
        AliESDtrack* track2 = static_cast<AliESDtrack*>(event->GetTrack(k2));
        if(!track2)continue;

        if (TMath::Abs((track2->Eta()))>0.8) continue;
        if (TMath::Abs((track2->Pt()))<0.15) continue;
        if (!fTrackFilter->IsSelected(track2)) continue;

        Int_t daug2Label = track2->GetLabel();
        if (daug2Label < 0) continue;
        if(daug1Label == daug2Label) continue;
        TParticle*  daughter2 = fMCStack->Particle(daug2Label);
        Int_t mother2Label = daughter2->GetFirstMother();
        TParticle* mother2 = fMCStack->Particle(mother2Label);
        Long_t mother2PDG = mother2->GetPdgCode();
        Int_t daug2PDG = daughter2->GetPdgCode();

        if(mother1Label == mother2Label){

          Double_t nSigmaTPCkPion1=-9999., nSigmaTPCkPion2=-9999., nSigmaTOFkPion1=-9999., nSigmaTOFkPion2=-9999., nSigmaTPCAcceptedkPion1=-9999., nSigmaTPCAcceptedkPion2=-9999., nSigmaTOFAcceptedkPion1=-9999., nSigmaTOFAcceptedkPion2=-9999.;

          Int_t PIDMethod = 4;

          /* PIDMethod = 1: 2sTPC_3sTOFveto; PIDMethod = 2: 3sTPC_3sTOFveto; PIDMethod = 3: 2sTPC_4sTOFveto; PIDMethod = 4: 2sTPC; PIDMethod = 5: 2sTOF */

          nSigmaTPCkPion1 = fPIDResponse->NumberOfSigmasTPC(track1, AliPID::kPion);
          nSigmaTOFkPion1 = fPIDResponse->NumberOfSigmasTOF(track1, AliPID::kPion);
          Bool_t TOFmatch1 = IsTOFMatched(track1);
          nSigmaTPCkPion2 = fPIDResponse->NumberOfSigmasTPC(track2, AliPID::kPion);
          nSigmaTOFkPion2 = fPIDResponse->NumberOfSigmasTOF(track2, AliPID::kPion);
          Bool_t TOFmatch2 = IsTOFMatched(track2);

          if(TOFmatch1){
            fAcceptedTOFMatchEtaPt->Fill(track1->Eta(), track1->Pt());
          }

          switch (PIDMethod) {
            case 1 /*2sTPC_3sTOFveto*/ : {
              if ( ((TMath::Abs(nSigmaTPCkPion1) < 2.) & (!TOFmatch1)) | ((TMath::Abs(nSigmaTPCkPion1) < 2.) & (TMath::Abs(nSigmaTOFkPion1) < 3.) & TOFmatch1) ){
                 v1.SetXYZM(track1->Px(), track1->Py(), track1->Pz(), piMass);
                 nSigmaTPCAcceptedkPion1 = fPIDResponse->NumberOfSigmasTPC(track1, AliPID::kPion);
                 nSigmaTOFAcceptedkPion1 = fPIDResponse->NumberOfSigmasTOF(track1, AliPID::kPion);
               }
              if ( ((TMath::Abs(nSigmaTPCkPion2) < 2.) & (!TOFmatch2)) | ((TMath::Abs(nSigmaTPCkPion2) < 2.) & (TMath::Abs(nSigmaTOFkPion2) < 3.) & TOFmatch2) ){
                 v2.SetXYZM(track2->Px(), track2->Py(), track2->Pz(), piMass);
                 nSigmaTPCAcceptedkPion2 = fPIDResponse->NumberOfSigmasTPC(track2, AliPID::kPion);
                 nSigmaTOFAcceptedkPion2 = fPIDResponse->NumberOfSigmasTOF(track2, AliPID::kPion);
               }
              break;
            }
            case 2 /*3sTPC_3sTOFveto*/ : {
              if ( ((TMath::Abs(nSigmaTPCkPion1) < 3.) & (!TOFmatch1)) | ((TMath::Abs(nSigmaTPCkPion1) < 3.) & (TMath::Abs(nSigmaTOFkPion1) < 3.) & TOFmatch1) ){
                v1.SetXYZM(track1->Px(), track1->Py(), track1->Pz(), piMass);
                nSigmaTPCAcceptedkPion1 = fPIDResponse->NumberOfSigmasTPC(track1, AliPID::kPion);
                nSigmaTOFAcceptedkPion1 = fPIDResponse->NumberOfSigmasTOF(track1, AliPID::kPion);
              }
              if ( ((TMath::Abs(nSigmaTPCkPion2) < 3.) & (!TOFmatch2)) | ((TMath::Abs(nSigmaTPCkPion2) < 3.) & (TMath::Abs(nSigmaTOFkPion2) < 3.) & TOFmatch2) ){
                v2.SetXYZM(track2->Px(), track2->Py(), track2->Pz(), piMass);
                nSigmaTPCAcceptedkPion2 = fPIDResponse->NumberOfSigmasTPC(track2, AliPID::kPion);
                nSigmaTOFAcceptedkPion2 = fPIDResponse->NumberOfSigmasTOF(track2, AliPID::kPion);
              }
              break;
            }
            case 3 /*2sTPC_4sTOFveto*/ : {
              if ( ((TMath::Abs(nSigmaTPCkPion1) < 2.) & (!TOFmatch1)) | ((TMath::Abs(nSigmaTPCkPion1) < 2.) & (TMath::Abs(nSigmaTOFkPion1) < 4.) & TOFmatch1) ){
                v1.SetXYZM(track1->Px(), track1->Py(), track1->Pz(), piMass);
                nSigmaTPCAcceptedkPion1 = fPIDResponse->NumberOfSigmasTPC(track1, AliPID::kPion);
                nSigmaTOFAcceptedkPion1 = fPIDResponse->NumberOfSigmasTOF(track1, AliPID::kPion);
              }
              if ( ((TMath::Abs(nSigmaTPCkPion2) < 2.) & (!TOFmatch2)) | ((TMath::Abs(nSigmaTPCkPion2) < 2.) & (TMath::Abs(nSigmaTOFkPion2) < 4.) & TOFmatch2) ){
                v2.SetXYZM(track2->Px(), track2->Py(), track2->Pz(), piMass);
                nSigmaTPCAcceptedkPion2 = fPIDResponse->NumberOfSigmasTPC(track2, AliPID::kPion);
                nSigmaTOFAcceptedkPion2 = fPIDResponse->NumberOfSigmasTOF(track2, AliPID::kPion);
              }
              break;
            }
            case 4 /*2sTPC*/ : {
              if ((TMath::Abs(nSigmaTPCkPion1) < 2.)) {
                v1.SetXYZM(track1->Px(), track1->Py(), track1->Pz(), piMass);
                nSigmaTPCAcceptedkPion1 = fPIDResponse->NumberOfSigmasTPC(track1, AliPID::kPion);
                }
              if ((TMath::Abs(nSigmaTPCkPion2) < 2.)) {
                v2.SetXYZM(track2->Px(), track2->Py(), track2->Pz(), piMass);
                nSigmaTPCAcceptedkPion2 = fPIDResponse->NumberOfSigmasTPC(track2, AliPID::kPion);
                }
              break;
            }
            case 5 /*2sTOF*/ : {
              if ((TMath::Abs(nSigmaTOFkPion1) < 2.) & TOFmatch1) {
                v1.SetXYZM(track1->Px(), track1->Py(), track1->Pz(), piMass);
                nSigmaTOFAcceptedkPion1 = fPIDResponse->NumberOfSigmasTOF(track1, AliPID::kPion);
                }
              if ((TMath::Abs(nSigmaTOFkPion2) < 2.) & TOFmatch2) {
                v2.SetXYZM(track2->Px(), track2->Py(), track2->Pz(), piMass);
                nSigmaTOFAcceptedkPion2 = fPIDResponse->NumberOfSigmasTOF(track2, AliPID::kPion);
                }
              break;
            }
            default: {
              Printf("Invalid method to perform the PID");
              break;
            }
          }

          vSum = v1+v2;

          fHistPID1tpc->Fill(v1.Pt(), nSigmaTPCAcceptedkPion1);
          fHistPID2tpc->Fill(v2.Pt(), nSigmaTPCAcceptedkPion2);
          fHistPID1tof->Fill(v1.Pt(), nSigmaTOFAcceptedkPion1);
          fHistPID2tof->Fill(v2.Pt(), nSigmaTOFAcceptedkPion2);
          fAcceptedTracksEtaPt->Fill(v1.PseudoRapidity(), v1.Pt());

          if ((TMath::Abs(vSum.Rapidity()))<0.5){ //rapidity cut
            for (Int_t iReco=0; iReco<10; iReco++){
              if(TMath::Abs(mother1->GetPdgCode()) == fPdgArray[iReco]){
                fRecoMassVsPt[iReco]->Fill(vSum.Mag(), vSum.Pt());
                fRecoYVsPt[iReco]->Fill(vSum.Pt(), vSum.Rapidity());
                fRecoEtaVsPt[iReco]->Fill(vSum.Pt(), vSum.PseudoRapidity());
              } //histos
            } //iReco
          } //rapidity cut
        } //if mom1=mom2
      } //i2Tracks
    }//if (k2!=k1)
  } //i1Tracks

  // continue until all the tracks are processed

  PostData(1, fOutputList);
} //esd track loop

//_____________________________________________________________________________
void AliAnalysisTaskF0Bkg::Terminate(Option_t *)
{
  // terminate
  // called at the END of the analysis (when all events are processed)
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskF0Bkg::SelectVertex2015pp(AliVEvent *event,
						Bool_t checkSPDres, //enable check on vtx resolution
						Bool_t requireSPDandTrk, //ask for both trk and SPD vertex
						Bool_t checkProximity, //apply cut on relative position of spd and trk verteces
						Bool_t enaMonitor) //enable monitoring and counters
{

  if (!event) return kFALSE;
  AliESDEvent *esd = dynamic_cast<AliESDEvent*>(event);
  if (!esd) return kFALSE;
  if (enaMonitor) fNEvents->Fill(1);

  const AliESDVertex * trkVertex = esd->GetPrimaryVertexTracks();
  const AliESDVertex * spdVertex = esd->GetPrimaryVertexSPD();
  Bool_t hasSPD = spdVertex->GetStatus();
  Bool_t hasTrk = trkVertex->GetStatus();

  //Note that AliVertex::GetStatus checks that N_contributors is > 0
  //reject events if both are explicitly requested and none is available
  if (requireSPDandTrk && !(hasSPD && hasTrk)) {
    if (enaMonitor) fNEvents->Fill(2);
    return kFALSE;
  }

  //reject events if none between the SPD or track verteces are available
  //if no trk vertex, try to fall back to SPD vertex;
  if (!hasTrk) {
    if (enaMonitor) fNEvents->Fill(3);
    if (!hasSPD) {
      if (enaMonitor) fNEvents->Fill(4);
      return kFALSE;
    }
    //on demand check the spd vertex resolution and reject if not satisfied
    if (checkSPDres && !IsGoodSPDvertexRes(spdVertex)) {
      if (enaMonitor) fNEvents->Fill(5);
      return kFALSE;
    }
  } else {
    if (hasSPD) {
      //if enabled check the spd vertex resolution and reject if not satisfied
      //if enabled, check the proximity between the spd vertex and trak vertex, and reject if not satisfied
      if (checkSPDres && !IsGoodSPDvertexRes(spdVertex)) {
	if (enaMonitor) fNEvents->Fill(5);
	return kFALSE;
      }
      if ((checkProximity && TMath::Abs(spdVertex->GetZ() - trkVertex->GetZ())>0.5)) {
	if (enaMonitor) fNEvents->Fill(6);
	return kFALSE;
      }
    }
  }

  //Cut on the vertex z position
  const AliESDVertex * vertex = esd->GetPrimaryVertex();
  if (TMath::Abs(vertex->GetZ())>10) {
    if (enaMonitor) fNEvents->Fill(7);
    return kFALSE;
  }
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskF0Bkg::IsGoodSPDvertexRes(const AliESDVertex * spdVertex)
{
  if (!spdVertex) return kFALSE;
  if (spdVertex->IsFromVertexerZ() && !(spdVertex->GetDispersion()<0.04 && spdVertex->GetZRes()<0.25)) return kFALSE;
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskF0Bkg::IsTOFMatched(AliESDtrack *track) {
  //This function checks whether a track has or has not a prolongation in the TOF. It returns true if the track has a matching hit in the TOF.
  bool hasTOFout  = track->GetStatus() & AliESDtrack::kTOFout;
  bool hasTOFtime = track->GetStatus() & AliESDtrack::kTIME;
  bool hasTPCout = track->GetStatus() & AliESDtrack::kTPCout;
  return (hasTOFout && hasTOFtime && hasTPCout);
}
