#if !defined (__CINT__) || defined (__CLING__)
#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisTaskF0Bkg.h"
#include "AliAnalysisTaskPIDResponse.h"
#include "AliPIDResponse.h"
#endif

void runAnalysisF0Bkg (Bool_t local=kFALSE, Bool_t gridTest=kTRUE)
{
    // local --> set if you want to run the analysis locally (kTRUE), or on grid (kFALSE)
    // gridTest --> if you run on grid, specify test mode (kTRUE) or full grid model (kFALSE)

    Int_t NtestFiles=3; //number of files you want to use for your pdgTest

    // since we will compile a class, tell root where to look for headers

#if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->ProcessLine(".include $ROOTSYS/include");
    gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
#else
    gROOT->ProcessLine(".include $ROOTSYS/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
#endif

    gSystem->Load("liblhapdf");
    gSystem->Load("libEGPythia6");
    gSystem->Load("libpythia6");
    gSystem->Load("libAliPythia6");

    // create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskExample");
    AliESDInputHandler *esdH = new AliESDInputHandler();
    mgr->SetInputEventHandler(esdH);
    AliMCEventHandler *mcHandler  = new AliMCEventHandler();
    mcHandler->SetReadTR(false);
    mgr->SetMCtruthEventHandler(mcHandler);

    // compile the class and load the add task macro
    // here we have to differentiate between using the just-in-time compiler
    // from root6, or the interpreter of root5


#if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->LoadMacro("${ALICE_PHYSICS}/OADB/macros/AddTaskCDBconnect.C");
    AliTaskCDBconnect *taskCDB = reinterpret_cast<AliTaskCDBconnect*>(gInterpreter->ExecuteMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCDBconnect.C"));

    gInterpreter->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physSelTask = reinterpret_cast<AliPhysicsSelectionTask*>(gInterpreter->ExecuteMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C(kTRUE)"));
    physSelTask->SelectCollisionCandidates(AliVEvent::kINT7);

    gInterpreter->LoadMacro("${ALICE_ROOT}/ANALYSIS/macros/AddTaskPIDResponse.C");
    AliAnalysisTaskPIDResponse* pidTask = reinterpret_cast<AliAnalysisTaskPIDResponse*>(gInterpreter->ExecuteMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"));

    gInterpreter->LoadMacro("AliAnalysisTaskF0Bkg.cxx++g");
    AliAnalysisTaskF0Bkg *task = reinterpret_cast<AliAnalysisTaskF0Bkg*>(gInterpreter->ExecuteMacro("AddTaskF0Bkg.C"));

  #else
      gROOT->LoadMacro("AliAnalysisTaskF0Bkg.cxx++g");
      gROOT->LoadMacro("AddTaskF0Bkg.C");
      AliAnalysisTaskF0Bkg *task = AddTaskF0Bkg();
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
      AliAnalysisTaskPIDResponse *pidTask = AddTaskPIDResponse(isMC);
      gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
      AliPhysicsSelectionTask *physSelTask = SelectCollisionCandidates(AliVEvent::kINT7);
      gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCDBconnect.C");
      AliTaskCDBconnect *taskCDB = AddTaskCDBconnect();
  #endif


    if(!mgr->InitAnalysis()) return;
    mgr->SetDebugLevel(2);
    mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 25);

    if(local) {
        // if you want to run locally, we need to define some input
        TChain* chain = new TChain("esdTree");
        // add a few files to the chain (change this so that your local files are added)
        chain->Add("AliESDs.root");
        // start the analysis locally, reading the events from the tchain
        mgr->StartAnalysis("local", chain);
    } else {
        // if we want to run on grid, we create and configure the plugin
        AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
        // also specify the include (header) paths on grid
        alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
        // make sure your source files get copied to grid
        alienHandler->SetAdditionalLibs("liblhapdf.so libEGPythia6.so libpythia6_4_25.so libAliPythia6.so AliAnalysisTaskF0Bkg.cxx AliAnalysisTaskF0Bkg.h AddTaskF0Bkg.C");
        alienHandler->SetAnalysisSource("AliAnalysisTaskF0Bkg.cxx");
        // select the aliphysics version. all other packages
        // are LOADED AUTOMATICALLY!
        alienHandler->SetAliPhysicsVersion("vAN-20171126-1");
        // set the Alien API version
        alienHandler->SetAPIVersion("V1.1x");
        // select the input data
        alienHandler->SetGridDataDir("/alice/sim/2017/LHC17d8");
        alienHandler->SetDataPattern("*/*AliESDs.root"); // if AOD --> */pass3/AOD/*AliAOD.root
        // MC has no prefix, data has prefix 000
        alienHandler->SetRunPrefix("");
        // runnumber
        /*Int_t run[25]={244628, 244627, 244626, 244619, 244618, 244617, 244542, 244540, 244531, 244484,
        244483, 244482, 244481, 244480, 244456, 244453, 244421, 244416, 244377, 244364,
        244359, 244355, 244351, 244343, 244340};

        for (Int_t iRun=0; iRun<25; iRun++){
        alienHandler->AddRunNumber(run[iRun]);
      }*/
        // to add only one run list -> comment lines 90-96 and use line 98
        alienHandler->AddRunNumber(244628);

        // number of files per subjob
        alienHandler->SetSplitMaxInputFileNumber(50);
        alienHandler->SetExecutable("myTaskF0bkg.sh");
        // specify how many seconds your job may take
        alienHandler->SetTTL(18000);
        alienHandler->SetJDLName("myTaskF0bkg.jdl");

        alienHandler->SetOutputToRunNo(kTRUE);
        alienHandler->SetKeepLogs(kTRUE);
        // merging: run with kTRUE to merge on grid
        // after re-running the jobs in SetRunMode("terminate")
        // (see below) mode, set SetMergeViaJDL(kFALSE)
        // to collect final results
        alienHandler->SetMaxMergeStages(2);
        alienHandler->SetMergeViaJDL(kTRUE);

        // define the output folders
        alienHandler->SetGridWorkingDir("F0bkg");
        alienHandler->SetGridOutputDir("20171127");

        // connect the alien plugin to the manager
        mgr->SetGridHandler(alienHandler);
        if(gridTest) {
            // speficy on how many files you want to run
            alienHandler->SetNtestFiles(NtestFiles);
            // and launch the analysis
            alienHandler->SetRunMode("test");
            mgr->StartAnalysis("grid");
        } else {
            // else launch the full grid analysis
            alienHandler->SetRunMode("full");
            mgr->StartAnalysis("grid");
        }
    }
}
