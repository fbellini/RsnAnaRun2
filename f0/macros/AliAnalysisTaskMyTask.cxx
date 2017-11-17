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

/* AliAnaysisTaskMyTask
 *
 * empty task which can serve as a starting point for building an analysis
 * as an example, one histogram is filled
 */

#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliAnalysisTaskMyTask.h"
#include <iostream>

class AliAnalysisTaskMyTask;
using namespace std;
ClassImp(AliAnalysisTaskMyTask)

AliAnalysisTaskMyTask::AliAnalysisTaskMyTask() : AliAnalysisTaskSE(),
  fPrimaryVertex(0),
  fOutputList(0),
  fHistF0Gen(0),
  fHistF0Reco(0),
  fHistF0daughtersPhaseSpace(0)
{
}
//_____________________________________________________________________________
AliAnalysisTaskMyTask::AliAnalysisTaskMyTask(const char* name) :
  AliAnalysisTaskSE(name),
  fPrimaryVertex(0),
  fOutputList(0),
  fHistF0Gen(0),
  fHistF0Reco(0),
  fHistF0daughtersPhaseSpace(0)
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskMyTask::~AliAnalysisTaskMyTask()
{

  delete fPrimaryVertex;
  if(fOutputList) {
    delete fOutputList;
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskMyTask::UserCreateOutputObjects()
{

  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);
  fHistF0Gen = new TH2F("fHistF0Gen","f0 generated; M_{#pi#pi} (GeV/#it{c}^{2}); #it{p}_{T} (GeV/#it{c})", 500, 0.7, 1.2, 200, 0., 10.);
  fHistF0Reco = new TH2F("fHistF0Reco","f0 reconstructed; M_{#pi#pi} (GeV/#it{c}^{2}); #it{p}_{T} (GeV/#it{c})", 500, 0.7, 1.2, 200, 0., 10.);
  //add here similar plots for reco and gen f0 filled with rapidity distributions in the range (-0.7, 0.7) vs pT (same range as above)
  //add here similar plots for reco and gen f0 filled with pseudorapidity distributions in the range (-1., 1) vs pT (same range as above)
  fHistF0daughtersPhaseSpace = new TH2F("fHistF0Reco","f0 generated; #it{p}_{T,1} (GeV/#it{c}); #it{p}_{T,2} (GeV/#it{c})", 200, 0.0, 10.0, 200, 0., 10.);
  //add plot with pT of 1st daughter vs pT of 2nd daughter for reconstructed f0 (1 plot)
  //add plot with opening angle of the two daughters
    
  fOutputList->Add(fHistF0Gen);
  fOutputList->Add(fHistF0Reco);

  PostData(1, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskMyTask::UserExec(Option_t *)
{
  const Long_t f0PDG = 9010221;
  const Long_t f2PDG = 225;
  const Int_t pionPDG = 211; //charged pion

  AliVEvent *event = InputEvent();
  if(!event) return; // if the pointer to the event is empty (getting it failed) skip this event

  //check here if the event is ESD or AOD and set an appropriate flag
  Bool_t isESD = kFALSE;
  if (dynamic_cast<AliESDEvent*>(event))
    isESD = kTRUE;
  else if (!dynamic_cast<AliAODEvent*>(event))
    AliFatal("I don't find the AOD event nor the ESD one, aborting.");


  // Retrieve primary vertex and check its existence
  const AliVVertex* vertex = event->GetPrimaryVertex();
  if (!vertex) {
    Printf("No primary vertex found!");
    return;
  }

  //Check vertex position along the beam axis
  if (TMath::Abs(vertex->GetZ())>10) return;


  //Add here additional vertex cuts if needed
  //Bool_t IncompleteDAQ = event->IsIncompleteDAQ();
  if(fPrimaryVertex) {
    //&& IncompleteDAQ){
    return;
  }

  //loop over ESD tracks
  //  for (Int_t iev = 0; iev < event->GetNumberOfTracks(); iev++) {

    //add here link between ESD track and the kinematics
    
   //loop on the MC generated (kinematics)
   AliStack *stack = event->Stack();
    if (!stack) AliFatal("AliStack not found.");
    
    for (Int_t ipart = 0; ipart < stack->GetNtrack(); ipart++) {
      
      /* get particle */
      TParticle *particle = stack->Particle(ipart);
      if (!particle) continue;
      
      Int_t momL = particle->GetMother(0);
      TParticle *mom = stack->Particle(momL);
      Long_t momPdg = mom->GetPdgCode();

      switch (momPdg) {
      case f0PDG :

	fHistF0Gen->Fill(mom->Pt());
	
	FillDaughterPhaseSpaceHisto(stack, particle);
	Printf("::: f0(980) found. Filling template histograms");
	break;
	
      case f2PDG :
	Printf("::: f2(1270) found. Filling template histograms");
	break;
	
      default:
	break;
      }
      
    }
    //}

  PostData(1, fOutputList);                           // stream the results the analysis of this event to
  // the output manager which will take care of writing
  // it to a file
}
//_____________________________________________________________________________
void AliAnalysisTaskMyTask::Terminate(Option_t *)
{
  // terminate
  // called at the END of the analysis (when all events are processed)
}

//_____________________________________________________________________________
void AliAnalysisTaskMyTask::FillDaughterPhaseSpaceHisto(AliStack * stack, TParticle * particle);
{
  //gets daughters 
  
  if (!stack) return;
  if (!particle) return;
  
  Int_t firstDaugL = particle->GetFirstDaughter();
  if (firstDaugL <= 0) continue;
  TParticle *firstDaug = stack->Particle(firstDaugL);
  Int_t secondDaugL = particle->GetDaughter(1);
  if (secondDaugL <= 0) continue;
  TParticle *secondDaug = stack->Particle(secondDaugL);
  Int_t lastDaugL = particle->GetLastDaughter();
  if (lastDaugL <= 0) continue;
  TParticle *lastDaug = stack->Particle(lastDaugL);
  
  Float_t firstDaugPt = 0.0;
  Float_t secondDaugPt = 0.0;
  
  if (firstDaug) firstDaugPt = firstDaug->Pt();
  if (secondDaug) secondDaugPt = secondDaug->Pt();

  if (fHistF0daughtersPhaseSpace)
    fHistF0daughtersPhaseSpace->Fill(firstDaugPt, secondDaugPt);
  
  return;
}
