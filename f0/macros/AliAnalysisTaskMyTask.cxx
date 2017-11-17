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
fPrimaryVertex(0), fOutputList(0), fHistF0Gen(0), fHistF0Reco(0)
{
}
//_____________________________________________________________________________
AliAnalysisTaskMyTask::AliAnalysisTaskMyTask(const char* name) : AliAnalysisTaskSE(name),
fPrimaryVertex(0), fOutputList(0), fHistF0Gen(0), fHistF0Reco(0)
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
    fHistF0Gen = new TH2F("fHistF0Gen","fHistF0Gen", 1000, 0., 10., 1000, 0., 10.);
    fHistF0Reco = new TH2F("fHistF0Reco","fHistF0Reco", 1000, 0., 10., 1000, 0., 10.);
    fOutputList->Add(fHistF0Gen);
    fOutputList->Add(fHistF0Reco);

    PostData(1, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskMyTask::UserExec(Option_t *)
{
    const Long_t f0PDG = 9010221;
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

    Bool_t IncompleteDAQ = event->IsIncompleteDAQ();
    if(fPrimaryVertex && IncompleteDAQ){
        return;
    }

    for (Int_t iev = 0; iev < event->GetNumberOfEvents(); iev++) {

    AliStack *stack = event->Stack();
    for (Int_t ipart = 0; ipart < stack->GetNtrack(); ipart++) {

    /* get particle */
    TParticle *particle = stack->Particle(ipart);
    if (!particle) continue;
    Float_t secondDaugPt = 0.0;
    Int_t momL = particle->GetMother(0);
    if (momL <= 0) continue;
    TParticle *mom = stack->Particle(momL);
    Int_t firstDaugL = particle->GetFirstDaughter();
    if (firstDaugL <= 0) continue;
    TParticle *firstDaug = stack->Particle(firstDaugL);
    Int_t secondDaugL = particle->GetDaughter(1);
    if (secondDaugL <= 0) continue;
    TParticle *secondDaug = stack->Particle(secondDaugL);
    Int_t lastDaugL = particle->GetLastDaughter();
    if (lastDaugL <= 0) continue;
    TParticle *lastDaug = stack->Particle(lastDaugL);
    if (firstDaug) firstDaugPt = firstDaug->Pt();
    if (secondDaug) secondDaugPt = secondDaug->Pt();

    if (mom->GetPdgCode()==f0PDG){
	     fHistF0Reco->Fill(mom->Pt());
     }
}
}

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
