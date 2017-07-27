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
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskMyTask.h"
#include <iostream>

class AliAnalysisTaskMyTask;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskMyTask) // classimp: necessary for root

AliAnalysisTaskMyTask::AliAnalysisTaskMyTask() : AliAnalysisTaskSE(), 
  fPrimaryVertex(0), fOutputList(0), fHistPt(0), fNEvents(0), fHistEtaPhi(0), fTPCSignal(0), fTOFSignal(0)
{
  // default constructor, don't allocate memory here!
  // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskMyTask::AliAnalysisTaskMyTask(const char* name) : AliAnalysisTaskSE(name),
								 fPrimaryVertex(0), fOutputList(0), fHistPt(0), fNEvents(0), fHistEtaPhi(0), fTPCSignal(0), fTOFSignal(0)
{
  // constructor
  DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
  // this chain is created by the analysis manager, so no need to worry about it, 
  // it does its work automatically
  DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 
  // you can add more output objects by calling DefineOutput(2, classname::Class())
  // if you add more output objects, make sure to call PostData for all of them, and to
  // make changes to your AddTask macro!
}
//_____________________________________________________________________________
AliAnalysisTaskMyTask::~AliAnalysisTaskMyTask()
{
  // destructor
  delete fPrimaryVertex;
  if(fOutputList) {
    delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskMyTask::UserCreateOutputObjects()
{
  // create output objects
  //
  // this function is called ONCE at the start of your analysis (RUNTIME)
  // here you ceate the histograms that you want to use 
  //
  // the histograms are in this case added to a tlist, this list is in the end saved
  // to an output file
  //
  fOutputList = new TList();          // this is a list which will contain all of your histograms
  // at the end of the analysis, the contents of this list are written
  // to the output file
  fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
  // if requested (dont worry about this now)

  // example of a histogram
  fHistPt = new TH1F("fHistPt", "fHistPt", 100, 0, 10);       // create your histogra
  fNEvents = new TH1I("fNEvents", "fNEvents", 4, 0, 4);
  fNEvents->GetXaxis()->SetBinLabel(1, "accepted");
  fNEvents->GetXaxis()->SetBinLabel(2, "processed");
  fNEvents->GetXaxis()->SetBinLabel(3, "with vtx");
  fNEvents->GetXaxis()->SetBinLabel(4, "|z_{vtx}|<10cm");
  fHistEtaPhi = new TH2F("fHistEtaPhi", "fHistEtaPhi", 400, -2.0, 2.0, 700, 0.0, 7.0);
  fHistEtaPhi->GetXaxis()->SetTitle("Eta");
  fHistEtaPhi->GetYaxis()->SetTitle("Phi");
  fTPCSignal = new TH2F("fPTPCSignal", "fPTPCSignal", 300, 0.0, 30.0, 10, -5.0, 5.0);
  fTPCSignal->GetXaxis()->SetTitle("Pions");
  fTPCSignal->GetYaxis()->SetTitle("TPC Signal");
  fTOFSignal = new TH2F("fTOFSignal", "fTOFSignal", 500, 0.0, 50.0, 100, -5.0, -5.0);
  fTOFSignal->GetXaxis()->SetTitle("Pions");
  fTOFSignal->GetYaxis()->SetTitle("TOF Signal");
    
    
  fOutputList->Add(fHistPt);          // don't forget to add it to the list! the list will be written to file, so if you want            // your histogram in the output file, add it to the list!
    
    
  fOutputList->Add(fNEvents);
  fOutputList->Add(fHistEtaPhi);
  fOutputList->Add(fTPCSignal);
  fOutputList->Add(fTOFSignal);

    
  PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the
  // fOutputList object. the manager will in the end take care of writing your output to file
  // so it needs to know what's in the output
}
//_____________________________________________________________________________
void AliAnalysisTaskMyTask::UserExec(Option_t *)
{
  // user exec
  // this function is called once for each event
  // the manager will take care of reading the events from file, and with the static function InputEvent() you 
  // have access to the current event. 
  // once you return from the UserExec function, the manager will retrieve the next event from the chain

  AliVEvent *event = InputEvent();
  if(!event) return; // if the pointer to the event is empty (getting it failed) skip this event
  fNEvents->Fill(1);

  //check here if the event is ESD or AOD and set an appropriate flag
  Bool_t isAOD = kFALSE;
  if (dynamic_cast<AliAODEvent*>(event))
    isAOD = kTRUE;
  else if (!dynamic_cast<AliESDEvent*>(event))
    AliFatal("I don't find the AOD event nor the ESD one, aborting.");
  
  // Retrieve primary vertex and check its existence
  const AliVVertex* vertex = event->GetPrimaryVertex();
  if (!vertex) {
    Printf("No primary vertex found!");
    return;
  }
  fNEvents->Fill(2);
  
  //Check vertex position along the beam axis
  if (TMath::Abs(vertex->GetZ())>10) return;
  fNEvents->Fill(3);

  //Add here additional vertex cuts if needed
	
   Bool_t IncompleteDAQ = event->IsIncompleteDAQ();
   if(fPrimaryVertex && IncompleteDAQ){
   return;
    }

  //fill histogram with number of accepted events
  fNEvents->Fill(0);
    
  // loop over the tracks in an event 
  // and extract some information from them which we'll store in a histogram
  Int_t iTracks(event->GetNumberOfTracks());           // see how many tracks there are in the event

  for(Int_t i(0); i < iTracks; i++) {                 // loop ove rall these tracks
    AliAODTrack* track = static_cast<AliAODTrack*>(event->GetTrack(i));         // get a track (type AliAODTrack) from the event
    if(!track) continue;                            // if we failed, skip this track
    if(!track->TestFilterBit(1)) continue;

    fHistPt->Fill(track->Pt());                     // plot the pt value of the track in a histogram
        
    //cout << track->Phi() << endl;
    //cout << track->Eta() << endl;
        
    fHistEtaPhi->Fill(track->Eta(),track->Phi());
    fTPCSignal->Fill(track->P(), track->GetTPCsignal());
    fTOFSignal->Fill(track->P(), track->GetTOFsignal());

        
  }                                                   // continue until all the tracks are processed
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
