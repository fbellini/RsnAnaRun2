/* Author: fbellini@cern.ch 
   Date of creation: 12.11.2017
*/
#if !defined (CINT) || defined (CLING)

#include "AliMultiInputEventHandler.h"
#include "AliMCEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"
#include "AliTaskCDBconnect.h"
#include "AliPhysicsSelectionTask.h"
#include "AliAnalysisDataContainer.h"
#include "AliMultSelectionTask.h"
#include "AliRsnMiniAnalysisTask.h"

#endif

TString AnalysisSetup(const char *options,
		      const char *outputFileName,
		      Int_t      ppass        = 1,        //reco pass needed for PIDresponse tune
		      Bool_t     enaMultSel   = kTRUE,    //enable multiplicity axis
		      Bool_t     enableMon    = kTRUE)    //enable single-track cuts monitoring
{
  // prepare output
  TString out("");
   
  // === EXAMINE OPTIONS ==========================================================================
  TString opt(options);
  opt.ToUpper();
  
  Bool_t isMC      = opt.Contains("MC") || (!opt.Contains("DATA"));
  Bool_t isPP      = opt.Contains("PP") || (!opt.Contains("PBPB"));
  Bool_t isESD     = opt.Contains("ESD");
   
  // === CREATE ANALYSIS MANAGER ==================================================================
  AliAnalysisManager *mgr = new AliAnalysisManager("RsnAnalysisManager");
  mgr->SetCommonFileName(outputFileName);

  // === INPUT / OUTPUT HANDLER CONFIGURATION =====================================================
  if (isESD) {
    out = "esdTree";
    ::Info("AnalysisSetup", "Creating ESD handler");
    AliESDInputHandler *esdHandler = new AliESDInputHandler();
    mgr->SetInputEventHandler(esdHandler);
    esdHandler->SetNeedField(); //needed to load magnetic field configuration
    
    if (isMC) {
      ::Info("AnalysisSetup", "Creating MC handler");
      AliMCEventHandler *mcHandler = new AliMCEventHandler();
      mcHandler->SetReadTR(false);
      mgr->SetMCtruthEventHandler(mcHandler);
    }
  } else {
    out = "aodTree";
    ::Info("AnalysisSetup", "Creating AOD handler");
    AliAODInputHandler *aodHandler = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodHandler);
  }

#if !defined (CINT) || defined (CLING)

  gInterpreter->LoadMacro("$ALICE_PHYSICS/PWGPP/PilotTrain/AddTaskCDBconnect.C+");
  AliTaskCDBconnect *taskCDB = AddTaskCDBconnect("raw://");
  if (!taskCDB) return;
  
  gInterpreter->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C+");
  AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(isMC);
									      
  AliAnalysisDataContainer *cstatsout = (AliAnalysisDataContainer*)mgr->GetOutputs()->FindObject("cstatsout");
  cstatsout->SetFileName("EventStatPhysSel.root");
  
  gInterpreter->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C+");
  AliMultSelectionTask* taskMult = 0X0;
  if (enaMultSel) {
    taskMult = reinterpret_cast<AliMultSelectionTask*>(gInterpreter->ExecuteMacro(AddTaskMultSelection())); //kTRUE, "B") to run in calibration mode
    if (isMC) taskMult->SetUseDefaultMCCalib(isMC);
  }
  
  AliPIDResponse *pidtask = reinterpret_cast<AliPIDResponse*>(gInterpreter->ExecuteMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C(isMC)"));
  AliAnalysisTask * pidRespTask = reinterpret_cast<AliAnalysisTask*>(gInterpreter->ExecuteMacro(AddTaskPIDResponse(isMC)));  //1 = MC

  gInterpreter->LoadMacro("$(ALICE_ROOT)/ANALYSIS/macros/AddTaskPIDqa.C");
  AliAnalysisTask* pidQAtask = 0x0;
  if (enableMon) pidQAtask = reinterpret_cast<AliAnalysisTask*>(gInterpreter->ExecuteMacro(AddTaskPIDqa()));
  
  gInterpreter->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/qa/AddTaskRsnQA.C");
  AliRsnMiniAnalysisTask * rsnQaTask = reinterpret_cast<AliRsnMiniAnalysisTask*>(gInterpreter->ExecuteMacro(AddTaskRsnQA(isMC, kFALSE, "V0M", AliVEvent::kINT7, "phi", 1, 0)));
   
  // gInterpreter->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddTaskPhiXeXe.C");
  gInterpreter->LoadMacro("$HOME/alice/resonances/RsnAnaRun2/phiXeXe/AddTaskPhiXeXe.C");
  AliRsnMiniAnalysisTask * rsnAnaTask = reinterpret_cast<AliRsnMiniAnalysisTask*>(gInterpreter->ExecuteMacro(AddTaskPhiXeXe(isMC)));
  
#else

  // === CDB connection =============================================================
  // Check that the correct trigger class is selected!!! USER DEFINED!!! 
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/PilotTrain/AddTaskCDBconnect.C");
  AliTaskCDBconnect *taskCDB = AddTaskCDBconnect("raw://");
  if (!taskCDB) return;
  
  // === PHYSICS SELECTION =============================================================
  // Check that the correct trigger class is selected!!! USER DEFINED!!!
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(isMC);
 
  //add also info on physics selection statistics in separate file
  AliAnalysisDataContainer *cstatsout = (AliAnalysisDataContainer*)mgr->GetOutputs()->FindObject("cstatsout");
  cstatsout->SetFileName("EventStatPhysSel.root");
  
  // === CENTRALITY/MULTIPLICITY ==============================================================
  //Use MultSelectionTask to provide centrlaity estimate
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  AliMultSelectionTask* taskMult = 0X0;
  if (enaMultSel) {
    taskMult = (AliMultSelectionTask *) AddTaskMultSelection(); //kTRUE, "B") to run in calibration mode
    if (isMC) taskMult->SetUseDefaultMCCalib(isMC);
  }
  
  // === PID RESPONSE =============================================================================
  // Set PIDresponse task   
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTask * pidRespTask = AddTaskPIDResponse(isMC,   //1 = MC
						     isMC,   //autoMCesd
						     isMC,   //tuneOnData
						     ppass,  //recopass
						     kFALSE, //cachePID, default
						     "",     //detResponse, default
						     kTRUE,  //useTPCEtaCorrection,/*Please use default value! Otherwise splines can be off*/
						     kFALSE, //useTPCMultiplicityCorrection --> default was kTRUE, but not avail for LHC13b2_efix_p1 
						     ppass); //recodatapass, default=-1
   
  // === PID QA ===================================================================================
  // Enable PIDqa - standard configuration  
  gROOT->LoadMacro("$(ALICE_ROOT)/ANALYSIS/macros/AddTaskPIDqa.C ");
  AliAnalysisTask* pidQAtask = 0x0;
  if (enableMon) AddTaskPIDqa();
   
  // === RSN TASKS ==============================================================================
  // Add resonance tasks
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/qa/AddTaskRsnQA.C ");
  AliRsnMiniAnalysisTask * rsnQaTask = AddTaskRsnQA(isMC, kFALSE, "V0M", AliVEvent::kINT7, "phi", 1, 0);
   
  // gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddTaskPhiXeXe.C");
  gROOT->LoadMacro("$HOME/alice/resonances/RsnAnaRun2/phiXeXe/AddTaskPhiXeXe.C");
  AliRsnMiniAnalysisTask * rsnAnaTask = AddTaskPhiXeXe(isMC);

#endif

  ::Info("<<<<<<<<<<<<<<<<< TRAIN CONFIGURATION >>>>>>>>>>>>>>>>");
  ::Info("AnalysisSetup", "Common file name: %s", outputFileName);
  if (physSelTask) ::Info("AnalysisSetup", "Added physics selection");
  if (cstatsout)   ::Info("AnalysisSetup", "Saving Physics Selection output in file EventStatPhysSel.root");
  if (taskMult)    ::Info("AnalysisSetup", "Add centrality/multiplicity computation tasks");
  if (pidRespTask) ::Info("AnalysisSetup", "Added task for PID response");
  if (pidQAtask)   ::Info("AnalysisSetup", "Added task for PID QA");
  if (rsnQaTask)   ::Info("AnalysisSetup", "Added task for Rsn QA"); 
  if (rsnAnaTask)  ::Info("AnalysisSetup", "Added task for Rsn Analysis");
  ::Info("AnalysisSetup", "Setup successful");
  return out;
}

  
  

void LoadLibs(Bool_t useTender = kFALSE)
{
  //
  // === LOAD LIBRARIES ===========================================================================
  // uncomment only if not already done in the steering macro, eg. RunGrid.C
  // load analysis libraries
  
  gSystem->AddIncludePath("-I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/ANALYSISalice");
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libOADB");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libEventMixing.so");
  gSystem->Load("libPWGLFresonances.so");
  //Pythia libs
  // gSystem->Load("libpythia6_4_28.so");
  // gSystem->Load("libpythia6.so");
  gSystem->Load("libAliPythia6.so");
  
  /*
    load development RSN library
    uncomment only for par files usage (not recommended)
    if (!AliAnalysisAlien::SetupPar("PWGLFresonances.par")) return "";
  */
  
  // tender-related libraries
  if (useTender) {
    ::Info("AnalysisSetup", "Loading tender libraries");
    gSystem->Load("libTENDER.so");
    gSystem->Load("libTENDERSupplies.so");
  }
  return;
}
   
