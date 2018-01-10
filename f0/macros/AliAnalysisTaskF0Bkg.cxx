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


////////////////////////////////////////////////////////////////////////
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
fESD(0), fOutputList(0), fNEvents(0), fMCEvent(0), fMCStack(0), fTrackFilter(0x0), fPID(0)
{
  // default constructor, don't allocate memory here!
  // this is used by root for IO purposes, it needs to remain empty

  for(Int_t i=0; i<10; i++){
    fHistPtGen[i]=0x0;
    fHistPtReco[i]=0x0;
    fHistYGen[i]=0x0;
    fHistYReco[i]=0x0;
    fHistEtaGen[i]=0x0;
    fHistEtaReco[i]=0x0;
  }
}
//_____________________________________________________________________________
AliAnalysisTaskF0Bkg::AliAnalysisTaskF0Bkg(const char* name) : AliAnalysisTaskSE(name),
fESD(0), fOutputList(0), fNEvents(0), fMCEvent(0), fMCStack(0), fTrackFilter(0x0), fPID(0)
{
  // constructor

  for(Int_t i=0; i<10; i++){
    fHistPtGen[i]=0x0;
    fHistPtReco[i]=0x0;
    fHistYGen[i]=0x0;
    fHistYReco[i]=0x0;
    fHistEtaGen[i]=0x0;
    fHistEtaReco[i]=0x0;
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
  if (fPID) delete fPID;
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

  fPID = inputHandler->GetPIDResponse();
  if (!fPID) AliError("PID Response object was not created");


  fOutputList = new TList();          // this is a list which will contain all of your histograms
  // at the end of the analysis, the contents of this list are written
  // to the output file
  fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
  // if requested (dont worry about this now)
  //fHlistPID = new TList();
  //fHlistPID->SetOwner(kTRUE);

  //  fPdgArray[10]={f0PDG, omegaPDG, rhoPDG, etaPDG, etaPrimePDG, f1PDG, f2PDG, kstarPDG, k0sPDG, phiPDG};

  //  fParticleName[10]={"f0","omega","rho","eta","etaPr","f1","f2","kStar", "k0s", "phi"};
  // fNEvents = new TH1I("fNEvents", "fNEvents", 4, 0, 4);
  // fNEvents->GetXaxis()->SetBinLabel(1, "accepted");
  // fNEvents->GetXaxis()->SetBinLabel(2, "processed");
  // fNEvents->GetXaxis()->SetBinLabel(3, "with vtx");
  // fNEvents->GetXaxis()->SetBinLabel(4, "|z_{vtx}|<10cm");
  // fOutputList->Add(fNEvents);

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


  for(Int_t j=0; j<10; j++){
    fHistPtGen[j]= new TH2F(Form("fHistPtGen_%s",fParticleName[j]), Form("generated %s; M_{#pi#pi} (GeV/#it{c}^{2}); #it{p}_{T} (GeV/#it{c})", fParticleName[j]), 1000, 0.3, 1.3, 220, 0., 11.);
    fHistPtReco[j] = new TH2F(Form("fHistPtReco_%s",fParticleName[j]), Form("reconstructed %s; M_{#pi#pi} (GeV/#it{c}^{2}); #it{p}_{T} (GeV/#it{c})", fParticleName[j]), 1000, 0.3, 1.3, 220, 0., 11.);
    fHistYGen[j] = new TH2F(Form("fHistYGen_%s",fParticleName[j]), Form("generated %s rapidity; #it{p}_{T} (GeV/#it{c}); #it{y}", fParticleName[j]), 220, 0., 11., 14, -0.7, 0.7);
    fHistYReco[j] = new TH2F(Form("fHistYReco_%s",fParticleName[j]), Form("reconstructed %s rapidity; #it{p}_{T} (GeV/#it{c}); #it{y}", fParticleName[j]), 220, 0., 11., 14, -0.7, 0.7);
    fHistEtaGen[j] = new TH2F(Form("fHistEtaGen_%s",fParticleName[j]), Form("generated %s pseudorapidity; #it{p}_{T} (GeV/#it{c}); #it{#eta}", fParticleName[j]), 220, 0., 11., 20, -1., 1.);
    fHistEtaReco[j] = new TH2F(Form("fHistEtaReco_%s",fParticleName[j]), Form("reconstructed %s pseudorapidity; #it{p}_{T} (GeV/#it{c}); #it{#eta}", fParticleName[j]), 220, 0., 11., 20, -1., 1.);


    fOutputList->Add(fHistPtGen[j]);          // don't forget to add it to the list! the list will be written to file, so if you want
    fOutputList->Add(fHistPtReco[j]);         // your histogram in the output file, add it to the list!
    fOutputList->Add(fHistYGen[j]);
    fOutputList->Add(fHistYReco[j]);
    fOutputList->Add(fHistEtaGen[j]);
    fOutputList->Add(fHistEtaReco[j]);


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
  
  /*AliPIDResponse::EStartTimeType_t startTimeMethodDefault = AliPIDResponse::kBest_T0;
  if (fESDpid->GetTPCPIDParams()) {  // during reconstruction OADB not yet available
    startTimeMethodDefault = ((AliTOFPIDParams *)fESDpid->GetTPCPIDParams())->GetStartTimeMethod();
  }*/

  //apply vertex selection and fill histogram with number of accepted events
  if (!SelectVertex2015pp(event, kTRUE, kFALSE, kTRUE, kTRUE)) return;
  fNEvents->Fill(0);

  //  fNEvents->Fill(1);
  // const AliVVertex* vertex = event->GetPrimaryVertex();
  // if (!vertex) {
  //   Printf("No primary vertex found!");
  //   fNEvents->Fill(2);
  //   return;
  // }
  // //fNEvents->Fill(2);

  // const AliVVertex* vertexSPD = event->GetPrimaryVertexSPD();
  // if (!vertexSPD) {
  //   Printf("No primary vertexSPD found!");
  //   fNEvents->Fill(3);
  //   return;
  // }

  // //Note that AliVertex::GetStatus checks that N_contributors is > 0
  // //reject events if both are explicitly requested and none is available
  // Bool_t hasSPD = vertexSPD->GetStatus();
  // Bool_t hasVtk = vertex->GetStatus();
  // if (requireSPDandTrk && !(hasSPD && hasTrk)) return kFALSE;
  

  // //vertex proximity cut
  // if (TMath::Abs(vertex->GetZ()-vertexSPD->GetZ()) <= 0.5){
  //   fNEvents->Fill(4);
  // } else {
  //   return;
  // }
  
  // if (!(TMath::Abs(vertex->GetZ())<10.)){
  //   fNEvents->Fill(5);
  //   return;
  // }
  //fNEvents->Fill(3);
  
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
  printf ("ciao\n");
  for (Int_t iStack = 0; iStack < fMCStack->GetNtrack(); iStack++) {

    trackMC = (TParticle*) fMCStack->Particle(iStack);
    ULong_t pdgTest=trackMC->GetPdgCode();

    for(Int_t index=0; index<10; index++){

      if (pdgTest==fPdgArray[index]){
        fHistPtGen[index]->Fill(trackMC->GetCalcMass(), trackMC->Pt());
        fHistYGen[index]->Fill(trackMC->Pt(), trackMC->Y());
        fHistEtaGen[index]->Fill(trackMC->Pt(), trackMC->Eta());
      } //endif
    } //pdg code loop
  } //stack loop

  ///////////// ------ reconstructed f0s ------ /////////////
  if(!fPID) Printf("::::: ERROR: no PID available.");

  TLorentzVector v1, v2;
  TLorentzVector vSum;
  const Float_t piMass = 0.13957;  // GeV/c^2

  Int_t iTracks = event->GetNumberOfTracks();           // see how many tracks there are in the event
  for(Int_t k1 = 0; k1 < iTracks; k1++) {                 // loop over all these tracks

    AliESDtrack* track1 = static_cast<AliESDtrack*>(event->GetTrack(k1));
    if(!track1) continue;
    if (!fTrackFilter->IsSelected(track1)) continue;
    if (TMath::Abs(fPID->NumberOfSigmasTPC(track1, AliPID::kPion) > 2.)) continue;
    v1.SetXYZM(track1->Px(), track1->Py(), track1->Pz(), piMass);
    Int_t daug1Label = track1->GetLabel();
    if (daug1Label < 0) continue;
    TParticle*  daughter1 = fMCStack->Particle(daug1Label);
    Int_t mother1Label = daughter1->GetFirstMother();
    TParticle* mother1 = fMCStack->Particle(mother1Label);
    Long_t mother1PDG = mother1->GetPdgCode();
    Int_t daug1PDG = daughter1->GetPdgCode();
    //printf("%ld\n", daug1PDG);
    if (TMath::Abs(daug1PDG) != 211) continue;
    
    //printf("******* k: %d, dau1label: %d, mother1label: %d, mother1PDG: %ld, daug1PDG: %d\n", k1, daug1Label, mother1Label, mother1PDG, daug1PDG);
    
    for(Int_t k2 = k1+1; k2 < iTracks; k2++){
      if (k2!=k1){
        AliESDtrack* track2 = static_cast<AliESDtrack*>(event->GetTrack(k2));
        if(!track2)continue;
        if (!fTrackFilter->IsSelected(track2)) continue;
        if (TMath::Abs(fPID->NumberOfSigmasTPC(track2, AliPID::kPion) > 2.)) continue;

        v2.SetXYZM(track2->Px(),track2->Py(), track2->Pz(), piMass);

	//apply rapidity cut to pair
	vSum = v1+v2;
        if (fabs(vSum.Rapidity())>0.5) continue; //rapidity cut
        Int_t daug2Label = track2->GetLabel();
        if (daug2Label < 0) continue;
        if(daug1Label == daug2Label) continue;
        TParticle*  daughter2 = fMCStack->Particle(daug2Label);
        Int_t mother2Label = daughter2->GetFirstMother();
        TParticle* mother2 = fMCStack->Particle(mother2Label);
        Long_t mother2PDG = mother2->GetPdgCode();
        Int_t daug2PDG = daughter2->GetPdgCode();

        //printf("Particella 2: dau2label: %d, mother2label: %d, mother2PDG: %ld, daug2PDG: %d\n", daug2Label, mother2Label, mother2PDG, daug2PDG);

        if(mother1Label == mother2Label){
          for (Int_t iReco=0; iReco<10; iReco++){
	    //compute invariant mass
            if(TMath::Abs(mother1->GetPdgCode()) == fPdgArray[iReco]){
              fHistPtReco[iReco]->Fill(vSum.Mag(), PairPt(track1,track2));
              fHistYReco[iReco]->Fill(PairPt(track1,track2), PairY(track1, track2));
              fHistEtaReco[iReco]->Fill(PairPt(track1, track2), PairEta(track1, track2));
            } //histos
          } //iReco
        } //if mom1=mom2
      } //i2Tracks
    }//if (k2!=k1)
  } //i1Tracks
  
  // continue until all the tracks are processed
  
  PostData(1, fOutputList);
} //esd track loop

//_____________________________________________________________________________
Float_t AliAnalysisTaskF0Bkg::PairPt(AliESDtrack* track1, AliESDtrack* track2){

  if ((!track1)||(!track2)){
    return 0x0;
  }
  Float_t py1, px1;
  Float_t py2, px2;
  px1 = track1->Px();
  py1 = track1->Py();
  px2 = track2->Px();
  py2 = track2->Py();
  Float_t pairPt=TMath::Sqrt((px1+px2)*(px1+px2)+(py1+py2)*(py1+py2));
  return pairPt;
}
//_____________________________________________________________________________
Float_t AliAnalysisTaskF0Bkg::PairEta(AliESDtrack* track1, AliESDtrack* track2){

  if ((!track1)||(!track2)){
    return 0x0;
  }
  Float_t theta1;
  Float_t theta2;

  theta1 = track1->Theta();
  theta2 = track2->Theta();
  Float_t pairEta=-1.*(TMath::Log((theta1+theta2)/2));
  return pairEta;
}
//_____________________________________________________________________________
Float_t AliAnalysisTaskF0Bkg::PairY(AliESDtrack* track1, AliESDtrack* track2){

  if ((!track1)||(!track2)){
    return 0x0;
  }
  Float_t e1, pz1;
  Float_t e2, pz2;

  e1 = track1->M();
  pz1 = track1->Pz();
  e2 = track2->M();
  pz2 = track2->Pz();

  Float_t PairY = 0.5*TMath::Log(((e1+e2)+(pz1+pz2))/((e1+e2)-(pz1+pz2)));
  return PairY;
}
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
