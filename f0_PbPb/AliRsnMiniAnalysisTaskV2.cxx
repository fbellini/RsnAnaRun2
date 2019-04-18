//
// Analysis task for 'mini' sub-package
// Contains all definitions needed for running an analysis:
// -- global event cut
// -- a list of track cuts (any number)
// -- definitions of output histograms
// -- values to be computed.
// Each one must be defined using the "CREATE" methods, which
// add directly a new element in the task collections, and don't
// need an external object to be passed to the task itself.
//
// Author: F. Bellini (fbellini@cern.ch)

#include <Riostream.h>

#include <TH1.h>
#include <TList.h>
#include <TTree.h>
#include <TStopwatch.h>
#include "TRandom.h"

#include "AliLog.h"
#include "AliEventplane.h"
#include "AliMultiplicity.h"
#include "AliTriggerAnalysis.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisUtils.h"

#include "AliESDtrackCuts.h"
#include "AliESDUtils.h"
#include "AliEventCuts.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"

#include "AliMultSelection.h"

#include "AliAnalysisTaskFlowVectorCorrections.h"
#include "AliQnCorrectionsManager.h"
#include "AliQnCorrectionsQnVector.h"

#include "AliRsnCutSet.h"
#include "AliRsnMiniPair.h"
#include "AliRsnMiniEvent.h"
#include "AliRsnMiniParticle.h"
#include "AliRsnMiniResonanceFinder.h"
#include "AliRsnMiniAnalysisTaskV2.h"
//#include "AliSpherocityUtils.h"

ClassImp(AliRsnMiniAnalysisTaskV2)

//__________________________________________________________________________________________________
AliRsnMiniAnalysisTaskV2::AliRsnMiniAnalysisTaskV2() :
   AliAnalysisTaskSE(),
   fUseMC(kFALSE),
   fEvNum(0),
   fTriggerMask(0),
   fSkipTriggerMask(0),
   fMultiEstimator("QUALITY"),
   fUseAOD049CentralityPatch(kFALSE),
   fUseCentralityPatchPbPb2011(0),
   fFlowQnVectorMgr(0),
   fFlowQnVectorSubDet("VZEROA"),
   fFlowQnVectorExpStep("latest"),
   fContinuousMix(kTRUE),
   fNMix(0),
   fMaxDiffMult(10),
   fMaxDiffVz(1.0),
   fMaxDiffAngle(1E20),
   fOutput(0x0),
   fHistograms("AliRsnMiniOutput", 0),
   fValues("AliRsnMiniValue", 0),
   fHEventStat(0x0),
   fHAEventsVsMulti(0x0),
   fHAEventsVsTracklets(0x0),
   fHAEventVzCent(0x0),
   fHAEventSpherocityCent(0x0),
   fHAEventPlane(0x0),
   fTrackFilter(0x0),
   fEventCuts(),
   fTrackCuts(),
   fRsnEvent(0x0),
   fEvBuffer(0x0),
   fESDtrackCuts(0x0),
   fMiniEvent(0x0),
   fBigOutput(kFALSE),
   fMixPrintRefresh(-1),
   fCheckDecay(kTRUE),
   fMaxNDaughters(-1),
   fCheckP(kFALSE),
   fCheckFeedDown(kFALSE),
   fOriginDselection(kFALSE),
   fKeepDfromB(kFALSE),
   fKeepDfromBOnly(kFALSE),
   fRejectIfNoQuark(kFALSE),
   fMotherAcceptanceCutMinPt(0.0),
   fMotherAcceptanceCutMaxEta(0.9),
   fKeepMotherInAcceptance(kFALSE),
   fRsnTreeInFile(kFALSE),
   fComputeSpherocity(kFALSE),
   fSpherocity(-10),
   fResonanceFinders(0)
{
//
// Dummy constructor ALWAYS needed for I/O.
//
}

//__________________________________________________________________________________________________
AliRsnMiniAnalysisTaskV2::AliRsnMiniAnalysisTaskV2(const char *name, Bool_t useMC,Bool_t saveRsnTreeInFile) :
   AliAnalysisTaskSE(name),
   fUseMC(useMC),
   fEvNum(0),
   fTriggerMask(AliVEvent::kINT7),
   fSkipTriggerMask(0),
   fMultiEstimator("ALIMULTSELECTION_V0M"),
   fUseAOD049CentralityPatch(kFALSE),
   fUseCentralityPatchPbPb2011(0),
   fFlowQnVectorMgr(0),
   fFlowQnVectorSubDet("VZEROA"),
   fFlowQnVectorExpStep("latest"),
   fContinuousMix(kTRUE),
   fNMix(0),
   fMaxDiffMult(10),
   fMaxDiffVz(1.0),
   fMaxDiffAngle(1E20),
   fOutput(0x0),
   fHistograms("AliRsnMiniOutput", 0),
   fValues("AliRsnMiniValue", 0),
   fHEventStat(0x0),
   fHAEventsVsMulti(0x0),
   fHAEventsVsTracklets(0x0),
   fHAEventVzCent(0x0),
   fHAEventSpherocityCent(0x0),
   fHAEventPlane(0x0),
   fEventCuts(0x0),
   fTrackCuts(0x0),
   fRsnEvent(0x0),
   fEvBuffer(0x0),
   fESDtrackCuts(0x0),
   fTrackFilter(0x0),
   fMiniEvent(0x0),
   fBigOutput(kFALSE),
   fMixPrintRefresh(-1),
   fCheckDecay(kTRUE),
   fMaxNDaughters(-1),
   fCheckP(kFALSE),
   fCheckFeedDown(kFALSE),
   fOriginDselection(kFALSE),
   fKeepDfromB(kFALSE),
   fKeepDfromBOnly(kFALSE),
   fRejectIfNoQuark(kFALSE),
   fMotherAcceptanceCutMinPt(0.0),
   fMotherAcceptanceCutMaxEta(0.9),
   fKeepMotherInAcceptance(kFALSE),
   fRsnTreeInFile(saveRsnTreeInFile),
   fComputeSpherocity(kFALSE),
   fSpherocity(-10),
   fResonanceFinders(0)
{
//
// Default constructor.
// Define input and output slots here (never in the dummy constructor)
// Input slot #0 works with a TChain - it is connected to the default input container
// Output slot #1 writes into a TH1 container
//
   DefineOutput(1, TList::Class());
   if (fRsnTreeInFile) DefineOutput(2, TTree::Class());
}

//__________________________________________________________________________________________________
AliRsnMiniAnalysisTaskV2::AliRsnMiniAnalysisTaskV2(const AliRsnMiniAnalysisTaskV2 &copy) :
   AliAnalysisTaskSE(copy),
   fUseMC(copy.fUseMC),
   fEvNum(0),
   fTriggerMask(copy.fTriggerMask),
   fSkipTriggerMask(copy.fSkipTriggerMask),
   fMultiEstimator(copy.fMultiEstimator),
   fUseAOD049CentralityPatch(copy.fUseAOD049CentralityPatch),
   fUseCentralityPatchPbPb2011(copy.fUseCentralityPatchPbPb2011),
   fFlowQnVectorMgr(copy.fFlowQnVectorMgr),
   fFlowQnVectorSubDet(copy.fFlowQnVectorSubDet),
   fFlowQnVectorExpStep(copy.fFlowQnVectorExpStep),
   fContinuousMix(copy.fContinuousMix),
   fNMix(copy.fNMix),
   fMaxDiffMult(copy.fMaxDiffMult),
   fMaxDiffVz(copy.fMaxDiffVz),
   fMaxDiffAngle(copy.fMaxDiffAngle),
   fOutput(0x0),
   fHistograms(copy.fHistograms),
   fValues(copy.fValues),
   fHEventStat(0x0),
   fHAEventsVsMulti(0x0),
   fHAEventsVsTracklets(0x0),
   fHAEventVzCent(0x0),
   fHAEventSpherocityCent(0x0),
   fHAEventPlane(0x0),
   fTrackFilter(copy.fTrackFilter),
   fEventCuts(copy.fEventCuts),
   fTrackCuts(copy.fTrackCuts),
   fRsnEvent(),
   fEvBuffer(0x0),
   fESDtrackCuts(copy.fESDtrackCuts),
   fMiniEvent(0x0),
   fBigOutput(copy.fBigOutput),
   fMixPrintRefresh(copy.fMixPrintRefresh),
   fCheckDecay(copy.fCheckDecay),
   fMaxNDaughters(copy.fMaxNDaughters),
   fCheckP(copy.fCheckP),
   fCheckFeedDown(copy.fCheckFeedDown),
   fOriginDselection(copy.fOriginDselection),
   fKeepDfromB(copy.fOriginDselection),
   fKeepDfromBOnly(copy.fKeepDfromBOnly),
   fRejectIfNoQuark(copy.fRejectIfNoQuark),
   fMotherAcceptanceCutMinPt(copy.fMotherAcceptanceCutMinPt),
   fMotherAcceptanceCutMaxEta(copy.fMotherAcceptanceCutMaxEta),
   fKeepMotherInAcceptance(copy.fKeepMotherInAcceptance),
   fRsnTreeInFile(copy.fRsnTreeInFile),
   fComputeSpherocity(copy.fComputeSpherocity),
   fSpherocity(copy.fSpherocity),
   fResonanceFinders(copy.fResonanceFinders)
{
//
// Copy constructor.
// Implemented as requested by C++ standards.
// Can be used in PROOF and by plugins.
//
}

//__________________________________________________________________________________________________
AliRsnMiniAnalysisTaskV2 &AliRsnMiniAnalysisTaskV2::operator=(const AliRsnMiniAnalysisTaskV2 &copy)
{
//
// Assignment operator.
// Implemented as requested by C++ standards.
// Can be used in PROOF and by plugins.
//

   AliAnalysisTaskSE::operator=(copy);
   if (this == &copy)
      return *this;
   fUseMC = copy.fUseMC;
   fEvNum = copy.fEvNum;
   fTriggerMask = copy.fTriggerMask;
   fSkipTriggerMask = copy.fSkipTriggerMask;
   fMultiEstimator = copy.fMultiEstimator;
   fUseAOD049CentralityPatch = copy.fUseAOD049CentralityPatch;
   fUseCentralityPatchPbPb2011 = copy.fUseCentralityPatchPbPb2011;
   fFlowQnVectorMgr = copy.fFlowQnVectorMgr;
   fFlowQnVectorSubDet = copy.fFlowQnVectorSubDet;
   fFlowQnVectorExpStep = copy.fFlowQnVectorExpStep;
   fContinuousMix = copy.fContinuousMix;
   fNMix = copy.fNMix;
   fMaxDiffMult = copy.fMaxDiffMult;
   fMaxDiffVz = copy.fMaxDiffVz;
   fMaxDiffAngle = copy.fMaxDiffAngle;
   fHistograms = copy.fHistograms;
   fValues = copy.fValues;
   fHEventStat = copy.fHEventStat;
   fHAEventsVsMulti = copy.fHAEventsVsMulti;
   fHAEventsVsTracklets = copy.fHAEventsVsTracklets;
   fHAEventVzCent = copy.fHAEventVzCent;
   fHAEventSpherocityCent = copy.fHAEventSpherocityCent;
   fHAEventPlane = copy.fHAEventPlane;
   fTrackFilter = copy.fTrackFilter;
   fEventCuts = copy.fEventCuts;
   fTrackCuts = copy.fTrackCuts;
   fESDtrackCuts = copy.fESDtrackCuts;
   fBigOutput = copy.fBigOutput;
   fMixPrintRefresh = copy.fMixPrintRefresh;
   fCheckDecay = copy.fCheckDecay;
   fMaxNDaughters = copy.fMaxNDaughters;
   fCheckP = copy.fCheckP;
   fCheckFeedDown = copy.fCheckFeedDown;
   fOriginDselection = copy.fOriginDselection;
   fKeepDfromB = copy.fOriginDselection;
   fKeepDfromBOnly = copy.fKeepDfromBOnly;
   fRejectIfNoQuark = copy.fRejectIfNoQuark;
   fMotherAcceptanceCutMinPt = copy.fMotherAcceptanceCutMinPt;
   fMotherAcceptanceCutMaxEta = copy.fMotherAcceptanceCutMaxEta;
   fKeepMotherInAcceptance = copy.fKeepMotherInAcceptance;
   fRsnTreeInFile = copy.fRsnTreeInFile;
   fComputeSpherocity = copy.fComputeSpherocity;
   fSpherocity = copy.fSpherocity;
   fResonanceFinders = copy.fResonanceFinders;

   return (*this);
}

//__________________________________________________________________________________________________
AliRsnMiniAnalysisTaskV2::~AliRsnMiniAnalysisTaskV2()
{
//
// Destructor.
// Clean-up the output list, but not the histograms that are put inside
// (the list is owner and will clean-up these histograms). Protect in PROOF case.
//
   if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
      delete fOutput;
      delete fEvBuffer;
   }
}

//__________________________________________________________________________________________________
Int_t AliRsnMiniAnalysisTaskV2::AddTrackCuts(AliRsnCutSet *cuts)
{
//
// Add a new cut set for a new criterion for track selection.
// A user can add as many as he wants, and each one corresponds
// to one of the available bits in the AliRsnMiniParticle mask.
// The only check is the following: if a cut set with the same name
// as the argument is there, this is not added.
// Return value is the array position of this set.
//

   TObject *obj = fTrackCuts.FindObject(cuts->GetName());
   Int_t v = 0;

   if (obj) {
      AliInfo(Form("A cut set named '%s' already exists", cuts->GetName()));
      v = fTrackCuts.IndexOf(obj);
   } else {
      fTrackCuts.AddLast(cuts);
      v = fTrackCuts.IndexOf(cuts);
   }

   for (Int_t i=0; i<fResonanceFinders.GetEntries(); i++){
      AliRsnMiniResonanceFinder* f = (AliRsnMiniResonanceFinder*) fResonanceFinders[i];
      if(f) f->IncrementResonanceCutID();
   }
    
   return v;
}

//__________________________________________________________________________________________________
void AliRsnMiniAnalysisTaskV2::UserCreateOutputObjects()
{
//
// Initialization of outputs.
// This is called once per worker node.
//

   // reset counter
   fEvNum = -1;

   // create list and set it as owner of its content (MANDATORY)
   if (fBigOutput) OpenFile(1);
   fOutput = new TList();
   fOutput->SetOwner();

   //initialise and config event selection via AliEventCuts
   fAliEventCuts.AddQAplotsToList(fOutput);
   fAliEventCuts.UseMultSelectionEventSelection(fUseMultSelectionEvSelection);
   AliInfo(Form("Use AliMultSelection Event Selection: %i", fUseMultSelectionEvSelection));

   // initialize ESD quality cuts
   if (fESDtrackCuts) delete fESDtrackCuts;
   fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010();

   //initialize quality trackcuts for spherocity
   if(!fTrackFilter){	
     fTrackFilter = new AliAnalysisFilter("trackFilter2015");
     SetTrackCuts(fTrackFilter);
   }

   // message
   AliInfo(Form("Multiplicity estimator: %s)", fMultiEstimator.Data()));

   // initialize event statistics counter rsn custom
   fHEventStat = new TH1F("hEventStat", "Event statistics", 7, 0.0, 7.0);
   fHEventStat->GetXaxis()->SetBinLabel(1, "Inspected");
   fHEventStat->GetXaxis()->SetBinLabel(2, "AliEventCuts ACCEPTED");
   fHEventStat->GetXaxis()->SetBinLabel(3, "AliEventCuts REJECTED");
   fHEventStat->GetXaxis()->SetBinLabel(4, "RsnCuts ACCEPTED");
   fHEventStat->GetXaxis()->SetBinLabel(5, "RsnCuts REJECTED");
   fHEventStat->GetXaxis()->SetBinLabel(6, "Accepted");
   fHEventStat->GetXaxis()->SetBinLabel(7, "Rejected");
   fOutput->Add(fHEventStat);

   Double_t multMin = -1.0;
   Double_t multMax = 3000.;
   fMultiEstimator.ToUpper();
   if (fMultiEstimator.Contains("ALIMULTSELECTION" || "ALICENTRALITY")) multMax = 100.0;
   Int_t nMultBinsX = multMax-multMin;

   if (!fHAEventsVsMulti) 
      fHAEventsVsMulti = new TH1F("hAEventsVsMulti", "Accepted events vs Multiplicity", multMin, multMax, nMultBinsX);
   fOutput->Add(fHAEventsVsMulti);
   
   if (!fHAEventsVsTracklets)
      fHAEventsVsTracklets = new TH1F("hAEventsVsTracklets", "Accepted events vs Tracklet Number", multMin, multMax, nMultBinsX);
   fOutput->Add(fHAEventsVsTracklets);

   if(fHAEventVzCent) fOutput->Add(fHAEventVzCent);
   if(fHAEventSpherocityCent) fOutput->Add(fHAEventSpherocityCent);
   if(fHAEventPlane) fOutput->Add(fHAEventPlane);

   AliAnalysisTaskFlowVectorCorrections *flowQnVectorTask = dynamic_cast<AliAnalysisTaskFlowVectorCorrections *>(AliAnalysisManager::GetAnalysisManager()->GetTask("FlowQnVectorCorrections"));
   if (flowQnVectorTask) {
     AliInfo("Using Flow Qn vector corrections framework task ...");
     fFlowQnVectorMgr = flowQnVectorTask->GetAliQnCorrectionsManager();
   }

   TIter next(&fTrackCuts);
   AliRsnCutSet *cs;
   while ((cs = (AliRsnCutSet *) next())) {
      cs->Init(fOutput);
   }

   // create temporary tree for filtered events
   if (fMiniEvent) SafeDelete(fMiniEvent);
   if (fRsnTreeInFile) OpenFile(2);
   fEvBuffer = new TTree("EventBuffer", "Temporary buffer for mini events");
   fMiniEvent = new AliRsnMiniEvent();
   fEvBuffer->Branch("events", "AliRsnMiniEvent", &fMiniEvent);
   
   // create one histogram per each stored definition (event histograms)
   Int_t i, ndef = fHistograms.GetEntries();
   AliRsnMiniOutput *def = 0x0;
   for (i = 0; i < ndef; i++) {
      def = (AliRsnMiniOutput *)fHistograms[i];
      if (!def) continue;
      if (!def->Init(GetName(), fOutput)) {
         AliError(Form("Def '%s': failed initialization", def->GetName()));
         continue;
      }
   }

   // post data for ALL output slots >0 here, to get at least an empty histogram
   PostData(1, fOutput);
   if (fRsnTreeInFile) PostData(2, fEvBuffer);
}

//__________________________________________________________________________________________________
void AliRsnMiniAnalysisTaskV2::UserExec(Option_t *)
{
//
// Computation loop.
// In this case, it checks if the event is acceptable, and eventually
// creates the corresponding mini-event and stores it in the buffer.
// The real histogram filling is done at the end, in "FinishTaskOutput".
//
   // increment event counter
   fEvNum++;
   
   // check current event against event cuts
   if (!CheckCurrentEvent()) return;

   // setup PID response
   AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
   AliInputEventHandler *inputHandler = (AliInputEventHandler *)man->GetInputEventHandler();
   fRsnEvent.SetPIDResponse(inputHandler->GetPIDResponse());

   // fill a mini-event from current
   // and skip this event if no tracks were accepted
   FillMiniEvent();

   for (Int_t i=0; i<fResonanceFinders.GetEntries(); i++){
      AliRsnMiniResonanceFinder* f = (AliRsnMiniResonanceFinder*) fResonanceFinders[i];
      if(f) f->RunResonanceFinder(fMiniEvent);
   }

   // fill MC based histograms on mothers,
   // which do need the original event
   if (fUseMC) {
      if (fRsnEvent.IsESD() && fMCEvent)
         FillTrueMotherESD(fMiniEvent);
      else if (fRsnEvent.IsAOD() && fRsnEvent.GetAODList())
         FillTrueMotherAOD(fMiniEvent);
   }

   // if the event is not empty, store it
   if (fMiniEvent->IsEmpty()) {
      AliDebugClass(2, Form("Rejecting empty event #%d", fEvNum));
   } else {
      Int_t id = fEvBuffer->GetEntries();
      AliDebugClass(2, Form("Adding event #%d with ID = %d", fEvNum, id));
      fMiniEvent->ID() = id;
      fEvBuffer->Fill();
   }

   // post data for computed stuff
   PostData(1, fOutput);
   if (fRsnTreeInFile) PostData(2, fEvBuffer);
}

//__________________________________________________________________________________________________
void AliRsnMiniAnalysisTaskV2::FinishTaskOutput()
{
//
// This function is called at the end of the loop on available events,
// and then the buffer will be full with all the corresponding mini-events,
// each one containing all tracks selected by each of the available track cuts.
// Here a loop is done on each of these events, and both single-event and mixing are computed
//

   // security code: reassign the buffer to the mini-event cursor
   fEvBuffer->SetBranchAddress("events", &fMiniEvent);
   TStopwatch timer;
   // prepare variables
   Int_t ievt, nEvents = (Int_t)fEvBuffer->GetEntries();
   Int_t idef, nDefs   = fHistograms.GetEntries();
   Int_t imix, iloop, ifill;
   AliRsnMiniOutput *def = 0x0;
   AliRsnMiniOutput::EComputation compType;

   Int_t printNum = fMixPrintRefresh;
   if (printNum < 0) {
      if (nEvents>1e5) printNum=nEvents/100;
      else if (nEvents>1e4) printNum=nEvents/10;
      else printNum = 0;
   }

   // loop on events, and for each one fill all outputs
   // using the appropriate procedure depending on its type
   // only mother-related histograms are filled in UserExec,
   // since they require direct access to MC event
   timer.Start();
   for (ievt = 0; ievt < nEvents; ievt++) {
      // get next entry
      fEvBuffer->GetEntry(ievt);
      if (printNum&&(ievt%printNum==0)) {
         AliInfo(Form("[%s] Std.Event %d/%d",GetName(), ievt,nEvents));
         timer.Stop(); timer.Print(); fflush(stdout); timer.Start(kFALSE);
      }
      // fill
      for (idef = 0; idef < nDefs; idef++) {
         def = (AliRsnMiniOutput *)fHistograms[idef];
         if (!def) continue;
         compType = def->GetComputation();
         // execute computation in the appropriate way
         switch (compType) {
            case AliRsnMiniOutput::kEventOnly:
               //AliDebugClass(1, Form("Event %d, def '%s': event-value histogram filling", ievt, def->GetName()));
               ifill = 1;
               def->FillEvent(fMiniEvent, &fValues);
               break;
            case AliRsnMiniOutput::kTruePair:
               //AliDebugClass(1, Form("Event %d, def '%s': true-pair histogram filling", ievt, def->GetName()));
               ifill = def->FillPair(fMiniEvent, fMiniEvent, &fValues);
               break;
            case AliRsnMiniOutput::kTrackPair:
               //AliDebugClass(1, Form("Event %d, def '%s': pair-value histogram filling", ievt, def->GetName()));
               ifill = def->FillPair(fMiniEvent, fMiniEvent, &fValues);
               break;
            case AliRsnMiniOutput::kTrackPairRotated1:
               //AliDebugClass(1, Form("Event %d, def '%s': rotated (1) background histogram filling", ievt, def->GetName()));
               ifill = def->FillPair(fMiniEvent, fMiniEvent, &fValues);
               break;
            case AliRsnMiniOutput::kTrackPairRotated2:
               //AliDebugClass(1, Form("Event %d, def '%s': rotated (2) background histogram filling", ievt, def->GetName()));
               ifill = def->FillPair(fMiniEvent, fMiniEvent, &fValues);
               break;
            default:
               // other kinds are processed elsewhere
               ifill = 0;
               AliDebugClass(2, Form("Computation = %d", (Int_t)compType));
         }
         // message
         AliDebugClass(1, Form("Event %6d: def = '%15s' -- fills = %5d", ievt, def->GetName(), ifill));
      }
   }

   // if no mixing is required, stop here and post the output
   if (fNMix < 1) {
      AliDebugClass(2, "Stopping here, since no mixing is required");
      PostData(1, fOutput);
      return;
   }

   // initialize mixing counter
   Int_t    nmatched[nEvents];
   TString *smatched = new TString[nEvents];
   for (ievt = 0; ievt < nEvents; ievt++) {
      smatched[ievt] = "|";
      nmatched[ievt] = 0;
   }


   AliInfo(Form("[%s] Std.Event %d/%d",GetName(), nEvents,nEvents));
   timer.Stop(); timer.Print(); timer.Start(); fflush(stdout);

   // search for good matchings
   for (ievt = 0; ievt < nEvents; ievt++) {
      if (printNum&&(ievt%printNum==0)) {
         AliInfo(Form("[%s] EventMixing searching %d/%d",GetName(),ievt,nEvents));
         timer.Stop(); timer.Print(); timer.Start(kFALSE); fflush(stdout);
      }
      if (nmatched[ievt] >= fNMix) continue;
      fEvBuffer->GetEntry(ievt);
      AliRsnMiniEvent evMain(*fMiniEvent);
      for (iloop = 1; iloop < nEvents; iloop++) {
         imix = ievt + iloop;
         if (imix >= nEvents) imix -= nEvents;
         if (imix == ievt) continue;
         // text next entry
         fEvBuffer->GetEntry(imix);
         // skip if events are not matched
         if (!EventsMatch(&evMain, fMiniEvent)) continue;
         // check that the array of good matches for mixed does not already contain main event
         if (smatched[imix].Contains(Form("|%d|", ievt))) continue;
         // check that the found good events has not enough matches already
         if (nmatched[imix] >= fNMix) continue;
         // add new mixing candidate
         smatched[ievt].Append(Form("%d|", imix));
         nmatched[ievt]++;
         nmatched[imix]++;
         if (nmatched[ievt] >= fNMix) break;
      }
      AliDebugClass(1, Form("Matches for event %5d = %d [%s] (missing are declared above)", evMain.ID(), nmatched[ievt], smatched[ievt].Data()));
   }

   AliInfo(Form("[%s] EventMixing searching %d/%d",GetName(),nEvents,nEvents));
   timer.Stop(); timer.Print(); fflush(stdout); timer.Start();

   // perform mixing
   TObjArray *list = 0x0;
   TObjString *os = 0x0;
   for (ievt = 0; ievt < nEvents; ievt++) {
      if (printNum&&(ievt%printNum==0)) {
         AliInfo(Form("[%s] EventMixing %d/%d",GetName(),ievt,nEvents));
         timer.Stop(); timer.Print(); timer.Start(kFALSE); fflush(stdout);
      }
      ifill = 0;
      fEvBuffer->GetEntry(ievt);
      AliRsnMiniEvent evMain(*fMiniEvent);
      list = smatched[ievt].Tokenize("|");
      TObjArrayIter next(list);
      while ( (os = (TObjString *)next()) ) {
         imix = os->GetString().Atoi();
         fEvBuffer->GetEntry(imix);
         for (idef = 0; idef < nDefs; idef++) {
            def = (AliRsnMiniOutput *)fHistograms[idef];
            if (!def) continue;
            if (!def->IsTrackPairMix()) continue;
            ifill += def->FillPair(&evMain, fMiniEvent, &fValues, kTRUE);
            if (!def->IsSymmetric()) {
               AliDebugClass(2, "Reflecting non symmetric pair");
               ifill += def->FillPair(fMiniEvent, &evMain, &fValues, kFALSE);
            }
         }
      }
      delete list;
   }

   delete [] smatched;

   AliInfo(Form("[%s] EventMixing %d/%d",GetName(),nEvents,nEvents));
   timer.Stop(); timer.Print(); fflush(stdout);

   // post computed data
   PostData(1, fOutput);
   if (fRsnTreeInFile) PostData(2, fEvBuffer);
}

//__________________________________________________________________________________________________
void AliRsnMiniAnalysisTaskV2::Terminate(Option_t *)
{
//
// Draw result to screen, or perform fitting, normalizations
// Called once at the end of the query
//

   fOutput = dynamic_cast<TList *>(GetOutputData(1));
   if (!fOutput) {
      AliError("Could not retrieve TList fOutput");
      return;
   }
}

//__________________________________________________________________________________________________
Bool_t AliRsnMiniAnalysisTaskV2::CheckCurrentEvent()
{
//
// This method checks if current event is OK for analysis.
// In case it is, the pointers of the local AliRsnEvent data member
// will point to it, in order to allow cut checking, otherwise the
// function exits with a failure message.
// 

   // string to sum messages
   TString msg("");

   // fill counter of inspected events
   fHEventStat->Fill(0.1);
   Bool_t isSelected = kFALSE;

   if (fAliEventCuts.AcceptEvent(fInputEvent)) {
      isSelected = kTRUE;  
      fHEventStat->Fill(1.1);
      msg += " -- AliEventCuts ACCEPTED";
   } else {
      msg += " -- AliEventCuts REJECTED";
      fHEventStat->Fill(2.1);
   }

   // if rsn event cuts are defined, they are checked here
   // final decision on the event depends on this
   if (fEventCuts) {
      if (fEventCuts->IsSelected(&fRsnEvent)) {
         isSelected = kTRUE;
         msg += " -- AliRsnCuts event ACCEPTED";
         fHEventStat->Fill(3.1);
      } else {
         isSelected = kFALSE;
         msg += " -- AliRsnCuts event REJECTED";
   	   fHEventStat->Fill(4.1);
      } 
   } else {
      isSelected = kTRUE;
      msg += " -- AliRsnCuts event NONE";
   }
     
   if (!isSelected) {
      fHEventStat->Fill(6.1);
      return msg;
   } 
   // if the above exit point is not taken, the event is accepted
   fHEventStat->Fill(5.1);
   // set reference for rsn event to input
   fRsnEvent.SetRef(fInputEvent);
   // add MC if requested and available
   if (fUseMC && fMCEvent) fRsnEvent.SetRefMC(fMCEvent);
   else fRsnEvent.SetRefMC(0x0);


   Double_t multi = ComputeMultiplicity();
   Printf(":::DEBUG AliRsnMiniAnalysisTaskV2::CheckCurrentEvent() -- multi = %5.2f ", multi);
   Double_t tracklets = ComputeTracklets();
   fSpherocity = (fComputeSpherocity) ? ComputeSpherocity() : -10;
   if (fHAEventsVsMulti) fHAEventsVsMulti->Fill(multi);
   if (fHAEventsVsTracklets) fHAEventsVsTracklets->Fill(tracklets);
   if (fHAEventVzCent) fHAEventVzCent->Fill(multi,fInputEvent->GetPrimaryVertex()->GetZ());
   if (fHAEventSpherocityCent && fComputeSpherocity) fHAEventSpherocityCent->Fill(multi,fSpherocity);
   if (fHAEventPlane) fHAEventPlane->Fill(multi,ComputeAngle());     
   return msg;
}

//__________________________________________________________________________________________________
void AliRsnMiniAnalysisTaskV2::FillMiniEvent()
{
//
// Refresh cursor mini-event data member to fill with current event.
// Returns the total number of tracks selected.
//
   fMiniEvent->Clear();
   // assign event-related values
   fMiniEvent->SetRef(fRsnEvent.GetRef());
   fMiniEvent->SetRefMC(fRsnEvent.GetRefMC());
   fMiniEvent->Vz() = fInputEvent->GetPrimaryVertex()->GetZ();
   fMiniEvent->Spherocity() = fSpherocity;
   fMiniEvent->Angle() = ComputeAngle();
   fMiniEvent->Mult() = ComputeMultiplicity();
   fMiniEvent->Tracklets() = ComputeTracklets();

   Printf(":::DEBUG AliRsnMiniAnalysisTaskV2::FillMiniEvent() -- multi = %5.2f ", fMiniEvent->Mult());

   
   if (fFlowQnVectorMgr) {
      TList *qnlist = fFlowQnVectorMgr->GetQnVectorList();
      if (qnlist) {
         fMiniEvent->SetQnVector(GetQnVectorFromList(qnlist, fFlowQnVectorSubDet.Data(), fFlowQnVectorExpStep.Data()));
      }
  }
   // loop on daughters and assign track-related values
   Int_t ic, ncuts = fTrackCuts.GetEntries();
   Int_t ip, npart = fRsnEvent.GetAbsoluteSum();
   Int_t npos = 0, nneg = 0, nneu = 0;
   AliRsnDaughter cursor;
  //  AliRsnMiniParticle miniParticle;
   AliRsnMiniParticle *miniParticlePtr;
   for (ip = 0; ip < npart; ip++) {
      // point cursor to next particle
      fRsnEvent.SetDaughter(cursor, ip);
      miniParticlePtr = fMiniEvent->AddParticle();
      miniParticlePtr->CopyDaughter(&cursor);
      miniParticlePtr->Index() = ip;
      
      // copy momentum and MC info if present
      // miniParticle.CopyDaughter(&cursor);
      // miniParticle.Index() = ip;
      // switch on the bits corresponding to passed cuts
      for (ic = 0; ic < ncuts; ic++) {
         AliRsnCutSet *cuts = (AliRsnCutSet *)fTrackCuts[ic];
        //  if (cuts->IsSelected(&cursor)) miniParticle.SetCutBit(ic);
         if (cuts->IsSelected(&cursor)) miniParticlePtr->SetCutBit(ic);
        }
        // continue;
       
        // if a track passes at least one track cut, it is added to the pool
      // if (miniParticle.CutBits()) {
        if (miniParticlePtr->CutBits()) {
          // fMiniEvent->AddParticle(miniParticle);
        // if (miniParticle.Charge() == '+') npos++;
        //  else if (miniParticle.Charge() == '-') nneg++;
        //  else nneu++;
        if (miniParticlePtr->Charge() == '+') npos++;
        else if (miniParticlePtr->Charge() == '-') nneg++;
        else nneu++;
     } else {
        TClonesArray &arr = fMiniEvent->Particles();
        // Printf("B %d",arr.GetEntries());
        arr.RemoveAt(arr.GetEntries()-1);
        // Printf("A %d",arr.GetEntries());
      }
   }
   // get number of accepted tracks
   AliDebugClass(1, Form("Event %6d: total = %5d, accepted = %4d (pos %4d, neg %4d, neu %4d)", fEvNum, npart, (Int_t)fMiniEvent->Particles().GetEntriesFast(), npos, nneg, nneu));
}

//__________________________________________________________________________________________________
Double_t AliRsnMiniAnalysisTaskV2::ComputeAngle()
{
//
// Get the plane angle
//

   AliEventplane *plane = 0x0;

   if (fInputEvent->InheritsFrom(AliESDEvent::Class()))
      plane = fInputEvent->GetEventplane();
   else if (fInputEvent->InheritsFrom(AliAODEvent::Class())) {
      AliAODEvent *aodEvent = (AliAODEvent *)fInputEvent;
      plane = ((AliVAODHeader*)aodEvent->GetHeader())->GetEventplaneP();
   }

   if (plane)
      return plane->GetEventplane("Q");
   else {
      AliWarning("No event plane defined");
      return 1E20;
   }
}

//__________________________________________________________________________________________________
Double_t AliRsnMiniAnalysisTaskV2::ComputeMultiplicity()
{
//
// Computes event centrality/multiplicity according to the crtiteria 
// defined in AliMultSelection framework. 
//Optionally, compute centrality with old frameworks for backward compatibility
fMultiEstimator.ToUpper();
if (fMultiEstimator.IsNull()) {
   AliWarning("Multiplicity/centrality estimator not defined - setting to -1.");
   return -1.0;
}

//Use AliMultSelection + estimator e.g. V0M or other
if (fMultiEstimator.Contains("ALIMULTSELECTION")) {
	TObjArray* a = (TObjArray*) fMultiEstimator.Tokenize("_");
	if (!a) {
      AliWarning("Problem with AliRsnMiniAnalysisTaskV2::fMultiEstimator string");
      return -2.0;
	}
   TObjString* o = (TObjString*) a->At(1);
   if (!o) {
      AliWarning("Problem with AliRsnMiniAnalysisTaskV2::fMultiEstimator string");
      return -3.0;
   }
   TString s = o->GetString();
   AliMultSelection *MultSelection = (AliMultSelection*) fInputEvent -> FindListObject("MultSelection");
	if (!MultSelection) {
      AliWarning("Could not find MultSelection object");
      return -4.0;
   }
   if (s.EqualTo("TEST")) {
	   MultSelection->PrintInfo();
	   return -10.;
	}
   return MultSelection->GetMultiplicityPercentile(s.Data());
} else 

//all charged tracks, no quality cuts
if (!fMultiEstimator.CompareTo("TRACKS"))
 return fInputEvent->GetNumberOfTracks();

//all charged tracks, passing standard quality cuts
if (!fMultiEstimator.CompareTo("QUALITY")){
   Double_t count = 0.;
   Int_t iTrack, ntracksLoop = fInputEvent->GetNumberOfTracks();
   for(Int_t iTrack = 0; iTrack < ntracksLoop; iTrack++){
      AliVTrack   *track = (AliVTrack *)fInputEvent->GetTrack(iTrack);
      AliAODTrack *aodt  = dynamic_cast<AliAODTrack *>(track);
      AliESDtrack *esdt  = dynamic_cast<AliESDtrack *>(track);
      if (aodt) if (!aodt->TestFilterBit(5)) continue;
      if (esdt) if (!fTrackFilter->IsSelected(esdt)) continue;
      count++;
   } 
 return count;
}

//AliCentrality framework, used in Run I analyses
if (!fMultiEstimator.Contains("ALICENTRALITY")) {
   if ((!fUseMC) && (fUseCentralityPatchPbPb2011))
      return ApplyCentralityPatchPbPb2011();
   
   if ((!fUseMC) && (!fInputEvent->InheritsFrom(AliESDEvent::Class())) && (fUseAOD049CentralityPatch)) 
      return ApplyCentralityPatchAOD049();
   
   AliCentrality *centrality = fInputEvent->GetCentrality();
   if (!centrality) {
      AliError("Cannot compute centrality!");
      return -5.0;
   }
   TObjArray* a = (TObjArray*) fMultiEstimator.Tokenize("_");
	if (!a) {
      AliWarning("Problem with AliRsnMiniAnalysisTaskV2::fMultiEstimator string");
      return -6.0;
	}
   TObjString* o = (TObjString*) a->At(1);
   if (!o) {
      AliWarning("Problem with AliRsnMiniAnalysisTaskV2::fMultiEstimator string");
      return -7.0;
   }
   TString estim = o->GetString();
   return centrality->GetCentralityPercentile(estim.Data());   
}

//Using tracklets
if (!fMultiEstimator.CompareTo("TRACKLETS")) {
   return ComputeTracklets();
}

AliError(Form("String '%s' does not define a possible multiplicity/centrality computation", fMultiEstimator.Data()));
return -10.0;
}

//__________________________________________________________________________________________________
Double_t AliRsnMiniAnalysisTaskV2::ComputeTracklets()
{
//
// Get number of tracklets in event
//
  Double_t nTr = 100;
  Double_t count = 0.0;

  if (fInputEvent->InheritsFrom(AliESDEvent::Class())){
    AliESDEvent *esdEvent = (AliESDEvent *)fInputEvent;
    const AliMultiplicity *spdmult = esdEvent->GetMultiplicity();
    nTr = 1.0*spdmult->GetNumberOfTracklets();
    for(Int_t iTr=0; iTr<nTr; iTr++){
      Double_t theta=spdmult->GetTheta(iTr);
      Double_t eta=-TMath::Log(TMath::Tan(theta/2.));
      if(TMath::Abs(eta)<1.0) count++;
    }
  }
  else if (fInputEvent->InheritsFrom(AliAODEvent::Class())) {
    AliAODEvent *aodEvent = (AliAODEvent *)fInputEvent;
    AliAODTracklets *spdmult = aodEvent->GetTracklets();
    nTr = 1.0*spdmult->GetNumberOfTracklets();
    for(Int_t iTr=0; iTr<nTr; iTr++){
      Double_t theta=spdmult->GetTheta(iTr);
      Double_t eta=-TMath::Log(TMath::Tan(theta/2.));
      if(TMath::Abs(eta)<1.0) count++;
    }
  }
  return count;
}

//__________________________________________________________________________________________________
Double_t AliRsnMiniAnalysisTaskV2::ComputeSpherocity()
{
//
// Get spherocity
//
  AliVEvent * evTypeS = InputEvent();
  Int_t ntracksLoop = evTypeS->GetNumberOfTracks();
  Int_t GoodTracks = 0;
  Float_t spherocity = -10.0;
  Float_t pFull = 0;
  Float_t Spherocity = 2;
  Float_t pt[10000],phi[1000];
  
  //computing total pt
  Float_t sumapt = 0;
  for(Int_t i1 = 0; i1 < ntracksLoop; ++i1){
    AliVTrack   *track = (AliVTrack *)evTypeS->GetTrack(i1);
    AliAODTrack *aodt  = dynamic_cast<AliAODTrack *>(track);
    AliESDtrack *esdt  = dynamic_cast<AliESDtrack *>(track);
    if (aodt) if (!aodt->TestFilterBit(5)) continue;
    if (esdt) if (!fTrackFilter->IsSelected(esdt)) continue;
    if (track->Pt() < 0.15) continue;
    if(TMath::Abs(track->Eta()) > 0.8) continue;
    pt[i1] = track->Pt();
    sumapt += pt[i1];
    GoodTracks++;
  }

  //Getting thrust
  for(Int_t i = 0; i < 360/0.1; ++i){
	Float_t numerador = 0;
	Float_t phiparam  = 0;
	Float_t nx = 0;
	Float_t ny = 0;
	phiparam=( (TMath::Pi()) * i * 0.1 ) / 180; // parametrization of the angle
	nx = TMath::Cos(phiparam);            // x component of an unitary vector n
	ny = TMath::Sin(phiparam);            // y component of an unitary vector n
	for(Int_t i1 = 0; i1 < ntracksLoop; ++i1){
	  AliVTrack   *track = (AliVTrack *)evTypeS->GetTrack(i1);
	  AliAODTrack *aodt  = dynamic_cast<AliAODTrack *>(track);
	  AliESDtrack *esdt  = dynamic_cast<AliESDtrack *>(track);
	  if (aodt) if (!aodt->TestFilterBit(5)) continue;
	  if (esdt) if (!fTrackFilter->IsSelected(esdt)) continue;
	  if (track->Pt() < 0.15) continue;
	  if(TMath::Abs(track->Eta()) > 0.8) continue;
	  pt[i1] = track->Pt();
	  phi[i1] = track->Phi();
	  Float_t pxA = pt[i1] * TMath::Cos( phi[i1] );
	  Float_t pyA = pt[i1] * TMath::Sin( phi[i1] );
	  numerador += TMath::Abs( ny * pxA - nx * pyA );//product between p  proyection in XY plane and the unitary vector
	}
	pFull=TMath::Power( (numerador / sumapt),2 );
	if(pFull < Spherocity)//maximization of pFull
	  {
	    Spherocity = pFull;
	  }
  }
  spherocity=((Spherocity)*TMath::Pi()*TMath::Pi())/4.0;
  if (GoodTracks > 2) return spherocity;
  else return -1.0;
}

//__________________________________________________________________________________________________
void AliRsnMiniAnalysisTaskV2::FillTrueMotherESD(AliRsnMiniEvent *miniEvent)
{
//
// Fills the histograms with generated mother (ESD version)
//

   Bool_t okMatch;
   Int_t id, ndef = fHistograms.GetEntries();
   Int_t ip, label1, label2, npart = fMCEvent->GetNumberOfTracks();
   static AliRsnMiniPair miniPair;
   AliMCParticle *daughter1, *daughter2;
   TLorentzVector p1, p2;
   AliRsnMiniOutput *def = 0x0;

   for (Int_t i=0; i<fResonanceFinders.GetEntries(); i++){
      AliRsnMiniResonanceFinder* f = (AliRsnMiniResonanceFinder*) fResonanceFinders[i];
      if(f) f->FillMother(fMCEvent, miniEvent);
   }

   for (id = 0; id < ndef; id++) {
      def = (AliRsnMiniOutput *)fHistograms[id];
      if (!def) continue;
      if (!def->IsMother() && !def->IsMotherInAcc()) continue;
      for (ip = 0; ip < npart; ip++) {
         AliMCParticle *part = (AliMCParticle *)fMCEvent->GetTrack(ip);
	 
         //get mother pdg code
         if (!AliRsnDaughter::IsEquivalentPDGCode(part->Particle()->GetPdgCode() , def->GetMotherPDG())) continue;
         // check that daughters match expected species
         if (part->Particle()->GetNDaughters() < 2) continue;
	 if (fMaxNDaughters > 0 && part->Particle()->GetNDaughters() > fMaxNDaughters) continue;
         label1 = part->Particle()->GetDaughter(0);
         label2 = part->Particle()->GetDaughter(1);
         daughter1 = (AliMCParticle *)fMCEvent->GetTrack(label1);
         daughter2 = (AliMCParticle *)fMCEvent->GetTrack(label2);
         okMatch = kFALSE;
         if (AliRsnDaughter::IsEquivalentPDGCode(TMath::Abs(daughter1->Particle()->GetPdgCode()) , def->GetPDG(0))
	     && AliRsnDaughter::IsEquivalentPDGCode(TMath::Abs(daughter2->Particle()->GetPdgCode()) , def->GetPDG(1))) {
            okMatch = kTRUE;
            p1.SetXYZM(daughter1->Px(), daughter1->Py(), daughter1->Pz(), def->GetMass(0));
            p2.SetXYZM(daughter2->Px(), daughter2->Py(), daughter2->Pz(), def->GetMass(1));
         } else if (AliRsnDaughter::IsEquivalentPDGCode(TMath::Abs(daughter1->Particle()->GetPdgCode()) , def->GetPDG(1))
		    && AliRsnDaughter::IsEquivalentPDGCode(TMath::Abs(daughter2->Particle()->GetPdgCode()) , def->GetPDG(0))) {
            okMatch = kTRUE;
            p1.SetXYZM(daughter1->Px(), daughter1->Py(), daughter1->Pz(), def->GetMass(1));
            p2.SetXYZM(daughter2->Px(), daughter2->Py(), daughter2->Pz(), def->GetMass(0));
         }
         if (fCheckDecay && !okMatch) continue;
	 if(fCheckP && (TMath::Abs(part->Px()-(daughter1->Px()+daughter2->Px()))/(TMath::Abs(part->Px())+1.e-13)) > 0.00001 &&
     				(TMath::Abs(part->Py()-(daughter1->Py()+daughter2->Py()))/(TMath::Abs(part->Py())+1.e-13)) > 0.00001 &&
     				(TMath::Abs(part->Pz()-(daughter1->Pz()+daughter2->Pz()))/(TMath::Abs(part->Pz())+1.e-13)) > 0.00001 ) continue;
	 if(fCheckFeedDown){
		Int_t pdgGranma = 0;
		Int_t mother = 0;
		mother = part->GetMother();
		Int_t istep = 0;
		Int_t abspdgGranma =0;
		Bool_t isFromB=kFALSE;
		Bool_t isQuarkFound=kFALSE;
		while (mother >=0 ){
			istep++;
			AliDebug(2,Form("mother at step %d = %d", istep, mother));
			AliMCParticle* mcGranma = dynamic_cast<AliMCParticle*>(fMCEvent->GetTrack(mother));
			if (mcGranma){
				pdgGranma = mcGranma->PdgCode();
				AliDebug(2,Form("Pdg mother at step %d = %d", istep, pdgGranma));
				abspdgGranma = TMath::Abs(pdgGranma);
				if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)){
				  isFromB=kTRUE;
				}
				if(abspdgGranma==4 || abspdgGranma==5) isQuarkFound=kTRUE;
				mother = mcGranma->GetMother();
			}else{
				AliError("Failed casting the mother particle!");
				break;
			}
		}
		if(fRejectIfNoQuark && !isQuarkFound) pdgGranma = -99999;
		if(isFromB){
		  if (!fKeepDfromB) pdgGranma = -9999; //skip particle if come from a B meson.
		}
		else{
		  if (fKeepDfromBOnly) pdgGranma = -999;
		  
		if (pdgGranma == -99999){
			AliDebug(2,"This particle does not have a quark in his genealogy\n");
			continue;
		}
		if (pdgGranma == -9999){
			AliDebug(2,"This particle come from a B decay channel but according to the settings of the task, we keep only the prompt charm particles\n");
			continue;
		}
		
		if (pdgGranma == -999){
			AliDebug(2,"This particle come from a prompt charm particles but according to the settings of the task, we want only the ones coming from B\n");
			continue;
		}

	      }
	 }
         // assign momenta to computation object
         miniPair.Sum(0) = miniPair.Sum(1) = (p1 + p2);
         miniPair.FillRef(def->GetMotherMass());
	 miniPair.P1(1) = p1;
	 miniPair.P2(1) = p2;

         // do computations and fill output
         def->FillMother(&miniPair, miniEvent, &fValues);
	 if(fKeepMotherInAcceptance){
	      if(daughter1->Pt()<fMotherAcceptanceCutMinPt || daughter2->Pt()<fMotherAcceptanceCutMinPt || TMath::Abs(daughter1->Eta())>fMotherAcceptanceCutMaxEta ||  TMath::Abs(daughter2->Eta())>fMotherAcceptanceCutMaxEta) continue;
	      def->FillMotherInAcceptance(&miniPair, miniEvent, &fValues);
	 }
	 
	 
      }
   }
}

//__________________________________________________________________________________________________
void AliRsnMiniAnalysisTaskV2::FillTrueMotherAOD(AliRsnMiniEvent *miniEvent)
{
//
// Fills the histograms with generated mother (AOD version)
//

   Bool_t okMatch;
   TClonesArray *list = fRsnEvent.GetAODList();
   Int_t id, ndef = fHistograms.GetEntries();
   Int_t ip, label1, label2, npart = list->GetEntries();
   static AliRsnMiniPair miniPair;
   AliAODMCParticle *daughter1, *daughter2;
   TLorentzVector p1, p2;
   AliRsnMiniOutput *def = 0x0;

   for (Int_t i=0; i<fResonanceFinders.GetEntries(); i++){
      AliRsnMiniResonanceFinder* f = (AliRsnMiniResonanceFinder*) fResonanceFinders[i];
      if(f) f->FillMother(list, miniEvent);
   }

   for (id = 0; id < ndef; id++) {
      def = (AliRsnMiniOutput *)fHistograms[id];
      if (!def) continue;
      if (!def->IsMother() && !def->IsMotherInAcc()) continue;
      for (ip = 0; ip < npart; ip++) {
         AliAODMCParticle *part = (AliAODMCParticle *)list->At(ip);
	 
         if (!AliRsnDaughter::IsEquivalentPDGCode(part->GetPdgCode() , def->GetMotherPDG())) continue;
         // check that daughters match expected species
         if (part->GetNDaughters() < 2) continue;
	 if (fMaxNDaughters > 0 && part->GetNDaughters() > fMaxNDaughters) continue;
         label1 = part->GetDaughterLabel(0);
         label2 = part->GetDaughterLabel(1);
         daughter1 = (AliAODMCParticle *)list->At(label1);
         daughter2 = (AliAODMCParticle *)list->At(label2);
         okMatch = kFALSE;
         if (AliRsnDaughter::IsEquivalentPDGCode(TMath::Abs(daughter1->GetPdgCode()) , def->GetPDG(0))
	     && AliRsnDaughter::IsEquivalentPDGCode(TMath::Abs(daughter2->GetPdgCode()) , def->GetPDG(1))) {
            okMatch = kTRUE;
            p1.SetXYZM(daughter1->Px(), daughter1->Py(), daughter1->Pz(), def->GetMass(0));
            p2.SetXYZM(daughter2->Px(), daughter2->Py(), daughter2->Pz(), def->GetMass(1));
         } else if (AliRsnDaughter::IsEquivalentPDGCode(TMath::Abs(daughter1->GetPdgCode()) , def->GetPDG(1))
		    && AliRsnDaughter::IsEquivalentPDGCode(TMath::Abs(daughter2->GetPdgCode()) , def->GetPDG(0))) {
            okMatch = kTRUE;
            p1.SetXYZM(daughter1->Px(), daughter1->Py(), daughter1->Pz(), def->GetMass(1));
            p2.SetXYZM(daughter2->Px(), daughter2->Py(), daughter2->Pz(), def->GetMass(0));
         }
         if (fCheckDecay && !okMatch) continue;
	 if(fCheckP && (TMath::Abs(part->Px()-(daughter1->Px()+daughter2->Px()))/(TMath::Abs(part->Px())+1.e-13)) > 0.00001 &&
     				(TMath::Abs(part->Py()-(daughter1->Py()+daughter2->Py()))/(TMath::Abs(part->Py())+1.e-13)) > 0.00001 &&
     				(TMath::Abs(part->Pz()-(daughter1->Pz()+daughter2->Pz()))/(TMath::Abs(part->Pz())+1.e-13)) > 0.00001 ) continue;
	 if(fCheckFeedDown){
		Int_t pdgGranma = 0;
		Int_t mother = 0;
		mother = part->GetMother();
		Int_t istep = 0;
		Int_t abspdgGranma =0;
		Bool_t isFromB=kFALSE;
		Bool_t isQuarkFound=kFALSE;
		while (mother >=0 ){
			istep++;
			AliDebug(2,Form("mother at step %d = %d", istep, mother));
			AliAODMCParticle* mcGranma = dynamic_cast<AliAODMCParticle*>(list->At(mother));
			if (mcGranma){
				pdgGranma = mcGranma->GetPdgCode();
				AliDebug(2,Form("Pdg mother at step %d = %d", istep, pdgGranma));
				abspdgGranma = TMath::Abs(pdgGranma);
				if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)){
				  isFromB=kTRUE;
				}
				if(abspdgGranma==4 || abspdgGranma==5) isQuarkFound=kTRUE;
				mother = mcGranma->GetMother();
			}else{
				AliError("Failed casting the mother particle!");
				break;
			}
		}
		if(fRejectIfNoQuark && !isQuarkFound) pdgGranma = -99999;
		if(isFromB){
		  if (!fKeepDfromB) pdgGranma = -9999; //skip particle if come from a B meson.
		}
		else{
		  if (fKeepDfromBOnly) pdgGranma = -999;
		  }
		  
		if (pdgGranma == -99999){
			AliDebug(2,"This particle does not have a quark in his genealogy\n");
			continue;
		}
		if (pdgGranma == -9999){
			AliDebug(2,"This particle come from a B decay channel but according to the settings of the task, we keep only the prompt charm particles\n");
			continue;
		}
		
		if (pdgGranma == -999){
			AliDebug(2,"This particle come from a prompt charm particles but according to the settings of the task, we want only the ones coming from B\n");
			continue;
		}
	 }
	 // assign momenta to computation object
         miniPair.Sum(0) = miniPair.Sum(1) = (p1 + p2);
         miniPair.FillRef(def->GetMotherMass());
	 miniPair.P1(1) = p1;
	 miniPair.P2(1) = p2;

         // do computations
         def->FillMother(&miniPair, miniEvent, &fValues);
	 if(fKeepMotherInAcceptance){
	      if(daughter1->Pt()<fMotherAcceptanceCutMinPt || daughter2->Pt()<fMotherAcceptanceCutMinPt || TMath::Abs(daughter1->Eta())>fMotherAcceptanceCutMaxEta ||  TMath::Abs(daughter2->Eta())>fMotherAcceptanceCutMaxEta) continue;
	      def->FillMotherInAcceptance(&miniPair, miniEvent, &fValues);
	 }
      }
   }
}
//___________________________________________________________
void AliRsnMiniAnalysisTaskV2::SetDselection(UShort_t originDselection)
{
	// setting the way the D0 will be selected
	// 0 --> only from c quarks
	// 1 --> only from b quarks
	// 2 --> from both c quarks and b quarks
		
	fOriginDselection = originDselection;
	
	if (fOriginDselection == 0) {
		fKeepDfromB = kFALSE;
		fKeepDfromBOnly = kFALSE;
	}
	
	if (fOriginDselection == 1) {
		fKeepDfromB = kTRUE;
		fKeepDfromBOnly = kTRUE;
	}
	
	if (fOriginDselection == 2) {
		fKeepDfromB = kTRUE;
		fKeepDfromBOnly = kFALSE;
	}
	
	return;
}
//__________________________________________________________________________________________________
Bool_t AliRsnMiniAnalysisTaskV2::EventsMatch(AliRsnMiniEvent *event1, AliRsnMiniEvent *event2)
{
//
// Check if two events are compatible.
// If the mixing is continuous, this is true if differences in vz, mult and angle are smaller than
// the specified values.
// If the mixing is binned, this is true if the events are in the same bin.
//

   if (!event1 || !event2) return kFALSE;
   Int_t ivz1, ivz2, imult1, imult2, iangle1, iangle2;
   Double_t dv, dm, da;

   if (fContinuousMix) {
      dv = TMath::Abs(event1->Vz()    - event2->Vz()   );
      dm = TMath::Abs(event1->Mult()  - event2->Mult() );
      da = TMath::Abs(event1->Angle() - event2->Angle());
      if (dv > fMaxDiffVz) {
         //AliDebugClass(2, Form("Events #%4d and #%4d don't match due to a too large diff in Vz = %f", event1->ID(), event2->ID(), dv));
         return kFALSE;
      }
      if (dm > fMaxDiffMult ) {
         //AliDebugClass(2, Form("Events #%4d and #%4d don't match due to a too large diff in Mult = %f", event1->ID(), event2->ID(), dm));
         return kFALSE;
      }
      if (da > fMaxDiffAngle) {
         //AliDebugClass(2, Form("Events #%4d and #%4d don't match due to a too large diff in Angle = %f", event1->ID(), event2->ID(), da));
         return kFALSE;
      }
      return kTRUE;
   } else {
      ivz1 = (Int_t)(event1->Vz() / fMaxDiffVz);
      ivz2 = (Int_t)(event2->Vz() / fMaxDiffVz);
      imult1 = (Int_t)(event1->Mult() / fMaxDiffMult);
      imult2 = (Int_t)(event2->Mult() / fMaxDiffMult);
      iangle1 = (Int_t)(event1->Angle() / fMaxDiffAngle);
      iangle2 = (Int_t)(event2->Angle() / fMaxDiffAngle);
      if (ivz1 != ivz2) return kFALSE;
      if (imult1 != imult2) return kFALSE;
      if (iangle1 != iangle2) return kFALSE;
      return kTRUE;
   }
}

//---------------------------------------------------------------------
Double_t AliRsnMiniAnalysisTaskV2::ApplyCentralityPatchPbPb2011(){
  //This part rejects randomly events such that the centrality gets flat for LHC11h Pb-Pb data
  //for 0-5% and 10-20% centrality bin
  
  if (fMultiEstimator!="V0M") {
    AliWarning("Wrong value (not centrality from V0).");
    return -999.0;
  }
  
  AliCentrality *centrality = fInputEvent->GetCentrality();
  if (!centrality) {
    AliWarning("Cannot get centrality from AOD event.");
    return -999.0;
  }
  
  Double_t cent = (Float_t)(centrality->GetCentralityPercentile("V0M"));
  Double_t rnd_hc = -1., testf = 0.0, ff = 0, N1 = -1., N2 = -1.;

  if(fUseCentralityPatchPbPb2011==510){
    N1 = 1.9404e+06;
    N2 = 1.56435e+06; //N2 is the reference
    ff = 5.04167e+06 - 1.49885e+06*cent + 2.35998e+05*cent*cent -1.22873e+04*cent*cent*cent;
  } else {
    if(fUseCentralityPatchPbPb2011==1020){
      N2 = 2.0e+05; //N2 is the reference
      N1 = 3.7e+05;
      ff = -1.73979e+06 - 3.05316e+06*cent + 1.05517e+06*cent*cent - 133205*cent*cent*cent + 8187.45*cent*cent*cent*cent - 247.875*cent*cent*cent*cent*cent + 2.9676*cent*cent*cent*cent*cent*cent;
    } else {
      AliError(Form("Patch for the requested centrality (%i) is not available", fUseCentralityPatchPbPb2011));
      return -999.0;
    }
  }
  testf = ( N2 + (N1-ff) ) / N1;
  rnd_hc = gRandom->Rndm();
  //AliDebugClass(1, Form("Flat Centrality %d", fUseCentralityPatchPbPb2011));
  if (rnd_hc < 0 || rnd_hc > 1 )
    {
      AliWarning("Wrong Random number generated");
      return -999.0;
    }
  
  if (rnd_hc < testf)
    return cent;
  else
    return -999.0;
}
//---------------------------------------------------------------------
Double_t AliRsnMiniAnalysisTaskV2::ApplyCentralityPatchAOD049()
{
   //
   //Apply centrality patch for AOD049 outliers
   //
   if (fInputEvent->InheritsFrom(AliESDEvent::Class())) {
      AliWarning("Requested patch for AOD049 for ESD. ");
      return -999.0;
   }

   if (fMultiEstimator!="V0M") {
      AliWarning("Requested patch forAOD049 for wrong value (not centrality from V0).");
      return -999.0;
   }

   AliCentrality *centrality = fInputEvent->GetCentrality();
   if (!centrality) {
      AliWarning("Cannot get centrality from AOD event.");
      return -999.0;
   }

   Float_t cent = (Float_t)(centrality->GetCentralityPercentile("V0M"));
   if(cent>=0.0) {
      Float_t v0 = 0.0;
      AliAODEvent *aodEvent = (AliAODEvent *)fInputEvent;
      AliAODVZERO *aodV0 = (AliAODVZERO *) aodEvent->GetVZEROData();
      v0+=aodV0->GetMTotV0A();
      v0+=aodV0->GetMTotV0C();
      if ( (cent==0) && (v0<19500) ) {
         AliDebug(3, Form("Filtering issue in centrality -> cent = %5.2f",cent));
         return -999.0;
      }
      Float_t tkl = (Float_t)(aodEvent->GetTracklets()->GetNumberOfTracklets());
      Float_t val = 1.30552 +  0.147931 * v0;

      Float_t tklSigma[101] = {176.644, 156.401, 153.789, 153.015, 142.476, 137.951, 136.127, 129.852, 127.436, 124.86,
                               120.788, 115.611, 113.172, 110.496, 109.127, 104.421, 102.479, 99.9766, 97.5152, 94.0654,
                               92.4602, 89.3364, 87.1342, 83.3497, 82.6216, 81.1084, 78.0793, 76.1234, 72.9434, 72.1334,
                               68.0056, 68.2755, 66.0376, 62.9666, 62.4274, 59.65, 58.3776, 56.6361, 54.5184, 53.4224,
                               51.932, 50.8922, 48.2848, 47.912, 46.5717, 43.4114, 43.2083, 41.3065, 40.1863, 38.5255,
                               37.2851, 37.5396, 34.4949, 33.8366, 31.8043, 31.7412, 30.8392, 30.0274, 28.8793, 27.6398,
                               26.6488, 25.0183, 25.1489, 24.4185, 22.9107, 21.2002, 21.6977, 20.1242, 20.4963, 19.0235,
                               19.298, 17.4103, 16.868, 15.2939, 15.2939, 16.0295, 14.186, 14.186, 15.2173, 12.9504, 12.9504,
                               12.9504, 15.264, 12.3674, 12.3674, 12.3674, 12.3674, 12.3674, 18.3811, 13.7544, 13.7544,
                               13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544
                              };

      if ( TMath::Abs(tkl-val) > 6.*tklSigma[(Int_t)cent] )  {
         AliDebug(3, Form("Outlier event in centrality -> cent = %5.2f",cent));
         return -999.0;
      }
   } else {
      //force it to be -999. whatever the negative value was
      cent = -999.;
   }
   return cent;
}

//----------------------------------------------------------------------------------
void AliRsnMiniAnalysisTaskV2::SetEventQAHist(TString type,TH1 *histo)
{
   if(!histo) {
      AliWarning(Form("event QA histogram pointer not defined for slot %s",type.Data()));
      return;
   }

   type.ToLower();
   TString multitype(histo->GetYaxis()->GetTitle());
   multitype.ToUpper();
   
   if(!type.CompareTo("eventsvsmulti")) fHAEventsVsMulti = (TH1F*) histo;
   else if(!type.CompareTo("eventsvstracklets")) fHAEventsVsTracklets = (TH1F*) histo;
   else if(!type.CompareTo("vz")) fHAEventVzCent = (TH2F*) histo;
   else if(!type.CompareTo("spherocitycent")){
      fHAEventSpherocityCent = (TH2F*) histo;
      fComputeSpherocity = kTRUE;
   }
   else if(!type.CompareTo("multicent")) {
      if(multitype.CompareTo("QUALITY") && multitype.CompareTo("TRACKS") && multitype.CompareTo("TRACKLETS")) {
         AliWarning(Form("multiplicity vs. centrality histogram y-axis %s unknown, setting to TRACKS",multitype.Data()));
         histo->GetYaxis()->SetTitle("TRACKS");
      }
   }
   else if(!type.CompareTo("refmulti")){
     if ( multitype.CompareTo("GLOBAL") && multitype.CompareTo("TRACKLETS") ) {
       AliWarning(Form("Reference multiplicity vs. centrality histogram y-axis %s unknown, setting to GLOBAL",multitype.Data()));
       histo->GetYaxis()->SetTitle("GLOBAL");
     }
   }
   else if(!type.CompareTo("eventplane")) fHAEventPlane = (TH2F*) histo;
   else AliWarning(Form("event QA histogram slot %s undefined",type.Data()));

   return;
}

//----------------------------------------------------------------------------------
Int_t AliRsnMiniAnalysisTaskV2::CreateValue(AliRsnMiniValue::EType type, Bool_t useMC)
{
//
// Create a new value in the task,
// and returns its ID, which is needed for setting up histograms.
// If that value was already initialized, returns its ID and does not recreate it.
//

   Int_t valID = ValueID(type, useMC);
   if (valID >= 0 && valID < fValues.GetEntries()) {
      AliInfo(Form("Value '%s' is already created in slot #%d", AliRsnMiniValue::ValueName(type, useMC), valID));
   } else {
      valID = fValues.GetEntries();
      AliInfo(Form("Creating value '%s' in slot #%d", AliRsnMiniValue::ValueName(type, useMC), valID));
      new (fValues[valID]) AliRsnMiniValue(type, useMC);
   }

   return valID;
}

//----------------------------------------------------------------------------------
Int_t AliRsnMiniAnalysisTaskV2::ValueID(AliRsnMiniValue::EType type, Bool_t useMC)
{
//
// Searches if a value computation is initialized
//

   const char *name = AliRsnMiniValue::ValueName(type, useMC);
   TObject *obj = fValues.FindObject(name);
   if (obj)
      return fValues.IndexOf(obj);
   else
      return -1;
}

//----------------------------------------------------------------------------------
AliRsnMiniOutput *AliRsnMiniAnalysisTaskV2::CreateOutput(const char *name, AliRsnMiniOutput::EOutputType type, AliRsnMiniOutput::EComputation src)
{
//
// Create a new histogram definition in the task,
// which is then returned to the user for its configuration
//

   Int_t n = fHistograms.GetEntries();
   AliRsnMiniOutput *newDef = new (fHistograms[n]) AliRsnMiniOutput(name, type, src);

   return newDef;
}

//----------------------------------------------------------------------------------
AliRsnMiniOutput *AliRsnMiniAnalysisTaskV2::CreateOutput(const char *name, const char *outType, const char *compType)
{
//
// Create a new histogram definition in the task,
// which is then returned to the user for its configuration
//

   Int_t n = fHistograms.GetEntries();
   AliRsnMiniOutput *newDef = new (fHistograms[n]) AliRsnMiniOutput(name, outType, compType);

   return newDef;
}

AliQnCorrectionsQnVector *AliRsnMiniAnalysisTaskV2::GetQnVectorFromList(
  const TList *list, const char *subdetector, const char *expectedstep) const {

  AliQnCorrectionsQnVector *theQnVector = NULL;

  // TList *pQvecList = dynamic_cast<TList *>(list->FindObject(subdetector));
  TList *pQvecList = (TList*)list->FindObject(subdetector);
  if (pQvecList != NULL) {
    /* the detector is present */
    if (TString(expectedstep).EqualTo("latest"))
      theQnVector = (AliQnCorrectionsQnVector *)pQvecList->First();
    else
      theQnVector =
        (AliQnCorrectionsQnVector *)pQvecList->FindObject(expectedstep);
  }
  if (theQnVector != NULL) {
    /* check the Qn vector quality */
    if (!(theQnVector->IsGoodQuality()) || !(theQnVector->GetN() != 0))
      /* not good quality, discarded */
      theQnVector = NULL;
  }
  return theQnVector;
}

//----------------------------------------------------------------------------------
Int_t AliRsnMiniAnalysisTaskV2::AddResonanceFinder(AliRsnMiniResonanceFinder* f)
{
   //
   // Add a new AliRsnMiniResonanceFinder object.
   // A user can add as many as he wants, and each one corresponds
   // to one of the available bits in the AliRsnMiniParticle mask.
   // The only check is the following: if a ResonanceFinder set with the same name
   // as the argument is there, this is not added.
   // Return value is the cut ID for the ResonanceFinder f.
   //

   TObject *obj = fResonanceFinders.FindObject(f->GetName());
   Int_t v = 0;

   if (obj) {
      AliInfo(Form("A ResonanceFinder named '%s' already exists", f->GetName()));
      v = fResonanceFinders.IndexOf(obj) + GetNumberOfTrackCuts();
   } else {
      fResonanceFinders.AddLast(f);
      v = fResonanceFinders.IndexOf(f) + GetNumberOfTrackCuts();
      f->SetResonanceCutID(v);
   }

   return v;
}

void AliRsnMiniAnalysisTaskV2::SetTrackCuts(AliAnalysisFilter* fTrackFilter)
{
	AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts();
	//TPC Only
	esdTrackCuts->SetMinNClustersTPC(50);
	esdTrackCuts->SetMaxChi2PerClusterTPC(4);
	esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
	esdTrackCuts->SetMaxDCAToVertexZ(3.2);
	esdTrackCuts->SetMaxDCAToVertexXY(2.4);
	esdTrackCuts->SetDCAToVertex2D(kTRUE);
	
	esdTrackCuts->SetRequireTPCRefit(kTRUE);// TPC Refit
	esdTrackCuts->SetRequireITSRefit(kTRUE);// ITS Refit
	fTrackFilter->AddCuts(esdTrackCuts);

}
