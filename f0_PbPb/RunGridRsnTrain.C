/*********************************************************************************************************
  Macro to launch analysis plugin for analysis
  Created by fbellini@cern.ch, 22/03/2018

Usage: 
- Modify the main settings (see description below) 
  IMPORTANT: the plugin needs you to configure the global variables,
  please DON'T SKIP THIS PASSAGE!!!
- Launch with AliRoot
  $ aliroot -l -b -q -x RunGridRsnTrain.C"(<mode>, <Grid User>, <suffix>, <ROOT>, <ALIROOT>, <ALIPHYSICS>)"

  <mode> can be "test", "full" or "terminate" - MANDATORY
  <Grid User> is your grid username - MANDATORY
  <suffix> is an additional string which identifies the names of the executables and jdl - OPTIONAL
           if it is empty (e.g., "") all files created by the plugin will have the same name 
  <ROOT>, <ALIROOT>, <ALIPHYSICS> must be compatible versions of the packages (check on Monalisa)
**********************************************************************************************************/

#ifdef __CLING__
#include "AliAODInputHandler.h"
// Tell  ROOT where to find AliRoot headers
R__ADD_INCLUDE_PATH($ALICE_ROOT)
#include <ANALYSIS/macros/train/AddESDHandler.C>
#include <ANALYSIS/macros/AddTaskPIDResponse.C>
#include <ANALYSIS/macros/AddTaskPIDqa.C>
// Tell ROOT where to find AliPhysics headers
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include <PWGPP/PilotTrain/AddTaskCDBconnect.C>
#include <OADB/macros/AddTaskPhysicsSelection.C>
#include <OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C>
#include <PWGLF/RESONANCES/macros/mini/AddTaskF0.C>
#include "AddTaskF0V2.C"
R__LOAD_LIBRARY(./AliRsnMiniAnalysisTaskV2_cxx.so)
#endif

TChain * CreateESDChain(TString esdpath=".", Int_t ifirst=-1, Int_t ilast=-1);

void     RunGridRsnTrain(TString pluginmode="test", Short_t ntest = 1, TString suffix="f0PbTest", Bool_t isLocalMerge = 1,
			                   TString dataset = "LHC15n", Bool_t isMC = 0, Int_t aodN = 0,
                  			 TString username="fbellini", TString aliPhysicsVer = "vAN-20190408_ROOT6-1");


//------------------------------------------------------------------------------------
// main
//------------------------------------------------------------------------------------
void RunGridRsnTrain(TString pluginmode, Short_t ntest, TString suffix, Bool_t isLocalMerge,
		     TString dataset, Bool_t isMC, Int_t aodN,
		     TString username, TString aliPhysicsVer)
{

  //------------------------------------------------------------------------------------
  //--------------------------------- LIBRARIES ---------------------------------------
  //-----------------------------------------------------------------------------------    
  printf("%s\n",gSystem->GetIncludePath());
  // Pythia libraries
  gSystem->Load("liblhapdf");
  gSystem->Load("libEGPythia6");
  gSystem->Load("libpythia6");
  gSystem->Load("libAliPythia6");
  //gSystem->Load("AliRsnMiniAnalysisTaskV2_cxx.so");
  //gSystem->Load("libPWGLFresonances.so");
  //gInterpreter->LoadMacro("AliRsnMiniAnalysisTaskV2.cxx++g");

  //------------------------------------------------------------------------------------
  //--------------------------------- DATASETS, RUNLIST --------------------------------
  //------------------------------------------------------------------------------------     
  Int_t yearMC;
  Int_t yearData;
  TString sim = "";
  TString data = "";
  TString passSuffix = "";
  Int_t ppass = 0;
  Int_t * runList;
  Int_t runNmin = 0;
  Int_t runNmax = 2;
  
  //trigger mask - to be specified for each period
  UInt_t triggerMask = AliVEvent::kINT7;

  // ESD/AOD analysis
  Bool_t isESD = ((aodN<0)? 1 : 0);
  
  //pp 5 TeV 2015 sample
  Int_t LHC15n[27] = {244340, 244343, 244351, 244355, 244359, 244364, 244377, 244411, 244416, 244418,
		      244421, 244453, 244456, 244480, 244481, 244482, 244483, 244484, 244531, 244540,
		      244542, 244617, 244618, 244619, 244626, 244627, 244628};
  if (dataset.Contains("15n")) {
    runList = LHC15n;
    data = "LHC15n";
    ppass = 3;
    sim =  "LHC17d8";
    passSuffix = "";
    yearData = 2015;
    yearMC = 2017;
    triggerMask = AliVEvent::kINT7;
  }

//pp 5 TeV 2017 samples
//run list for CENT_wSDD production: https://twiki.cern.ch/twiki/bin/viewauth/ALICE/AliDPGRunList17p
// remove 282030 for FAST
  Int_t LHC17p[42] = {282343, 282342, 282341, 282340, 282314, 282313, 282312, 282309, 282307, 282306, 
                      282305, 282304, 282303, 282302, 282247, 282230, 282229, 282227, 282224, 282206, 
                      282189, 282147, 282146, 282127, 282126, 282125, 282123, 282122, 282120, 282119, 
                      282118, 282099, 282098, 282078, 282051, 282050, 282031, 282030, 282025, 282021, 
                      282016, 282008};
  
  if (dataset.Contains("17p")) {
    runList = LHC17p;
    data = "LHC17p";
    ppass = 1;
    sim =  "LHC18c8b";
    passSuffix = "";
    yearData = 2017;
    yearMC = 2018;
    triggerMask = AliVEvent::kINT7;
  }

//run list for CENT_wSDD production: https://twiki.cern.ch/twiki/bin/view/ALICE/AliDPGRunList17q 
// run list is the same for all subcycles, so also FAST and CENT_woSDD
  Int_t LHC17q[3] = {282367, 282366, 282365};
  
  if (dataset.Contains("17q")) {
    runList = LHC17q;
    data = "LHC17q";
    ppass = 1;
    sim =  "LHC18c8b";
    yearData = 2017;
    yearMC = 2018;
    triggerMask = AliVEvent::kINT7;
  }

  //Xe-Xe 5.44 TeV 2017 sample
  Int_t LHC17n[2] = {280234, 280235};
      
  if (dataset.Contains("17n")) {
    runList = LHC17n;
    data = "LHC17n";
    ppass = 1;
    sim =  "LHC17j7";//LHC17j7_ZDCfix, LHC17j7_ZDCfix_extra
    yearData = 2017;
    yearMC = 2017;
    triggerMask = AliVEvent::kINT7;
  }
  
//Pb-Pb 5.02 TeV 2018 sample
  Int_t LHC18q[126] = {296623, 296622, 296621, 296619, 296618, 296616, 296615, 296594, 296553, 296552, 
                       296551, 296550, 296549, 296548, 296547, 296516, 296512, 296511, 296510, 296509,
                       296472, 296433, 296424, 296423, 296420, 296419, 296415, 296414, 296383, 296381, 
                       296380, 296379, 296378, 296377, 296376, 296375, 296312, 296309, 296304, 296303, 
                       296280, 296279, 296273, 296270, 296269, 296247, 296246, 296244, 296243, 296242, 
                       296241, 296240, 296198, 296197, 296196, 296195, 296194, 296192, 296191, 296143, 
                       296142, 296135, 296134, 296133, 296132, 296123, 296074, 296066, 296065, 296063, 
                       296062, 296060, 296016, 295942, 295941, 295937, 295936, 295913, 295910, 295909, 
                       295861, 295860, 295859, 295856, 295855, 295854, 295853, 295831, 295829, 295826, 
                       295825, 295822, 295819, 295818, 295816, 295791, 295788, 295786, 295763, 295762, 
                       295759, 295758, 295755, 295754, 295725, 295723, 295721, 295719, 295718, 295717, 
                       295714, 295712, 295676, 295675, 295673, 295668, 295667, 295666, 295615, 295612, 
                       295611, 295610, 295589, 295588, 295586, 295585};
      
  if (dataset.Contains("18q")) {
    runList = LHC18q;
    data = "LHC18q";
    ppass = 1;
    sim = "LHC18l8b2"; //GP central 0-10%   "LHC18l8b2_cent2mb" //GP 0-10% , "LHC18l8c2_mb2cent"//GP 30-50% 
    yearData = 2018;
    yearMC = 2018;
    triggerMask = (AliVEvent::kINT7|AliVEvent::kCentral|AliVEvent::kSemiCentral);
  }

  //------------------------------------------------------------------------------------
  //--------------------------------- DATASET and PATHS -------------------------------
  //------------------------------------------------------------------------------------ 
  /* Define data path */
  TString myGridDataDir="/alice";
  TString myDataPattern="";
  
  /* Define output directory and names of the files generated by the plugin.
     Omit the extension of the files, it will be automatically set by the SetupIO function */
  TString myOutDir  = Form("f0%s_%s", (isMC? sim.Data() : data.Data()), suffix.Data());
  TString myWorkDir = "Resonances";
  TString myJDLname = "jobF0Pb";
  TString myExecutableName = "RunF0Pb";
  TString myMacroName = "RunF0Pb";
  
  //Setup I/O paths and file names
  if (isMC) myGridDataDir.Append(Form("/sim/%i/%s",yearMC, sim.Data()));
  else myGridDataDir.Append(Form("/data/%i/%s",yearData, data.Data()));
  
  if (isESD) {
    if (isMC) myDataPattern="*AliESDs.root";
    else myDataPattern=Form("/pass%i/*AliESDs.root",ppass);
  } else {
    if (aodN>0) myDataPattern=Form("AOD%03i/*AliAOD.root",aodN);
    else myDataPattern="AOD/*AliAOD.root";
    if (!isMC) myDataPattern.Prepend(Form("*/pass%i/",ppass)); 
  }
  
  myJDLname.Append(suffix.Data()); myJDLname.Append(".jdl");
  myExecutableName.Append(suffix.Data()); myExecutableName.Append(".sh");
  myMacroName.Append(suffix.Data()); myMacroName.Append(".C");
  
  if (myOutDir.IsNull()) myOutDir = Form("%s_pass%i_%s", (isMC ? sim.Data() : data.Data()), ppass, suffix.Data());
  
  Printf("=========================================  Setup I/O:");
  Printf("myGridDataDir = %s", myGridDataDir.Data());
  Printf("myDataPattern = %s", myDataPattern.Data());
  Printf("myWorkDir = %s", myWorkDir.Data());
  Printf("myOutDir = %s", myOutDir.Data());
  Printf("myJDLname = %s", myJDLname.Data());
  Printf("myExecutableName = %s", myExecutableName.Data());
  Printf("myMacroName = %s", myMacroName.Data());
  Printf("=======================================================");

  //------------------------------------------------------------------------------------
  //--------------------------------- ANALYSIS MODE ------------------------------------
  //------------------------------------------------------------------------------------ 
  /* Set the usage of the alien plugin to kTRUE to run the analysis on grid or in test mode,
     set it to kFALSE to run test locally.
     If you want to run merging via plugin, disable the local merge when running in "full" mode or "terminate" mode (multiple merging stages are possible).
     Enable the local merge at the last step of the "terminate" to download locally the final merged output.
     Otherwise, local merge is always possible. */
  TString analysisMode  = "grid"; //in alternative: "local"
  Bool_t useAlienPlugin = (pluginmode.IsNull()? kFALSE : kTRUE);
  
  if(analysisMode=="grid") {
    if (!useAlienPlugin) analysisMode = "local";
    TGrid::Connect("alien://");
  }
  
  //needed for running locally
  //The file added to the chain can also be a file on the grid, e.g. alien://<full path on grid>/AliESDs.root
  Long64_t nentries=100000;
  Long64_t firstentry=0; //needed to read input when running locally
  TChain *chain = 0x0;
  AliAnalysisAlien *plugin = 0x0; 
  
  //-----------------------------------------------------------------------------------
  //--------------------------------- ALIEN PLUGIN ------------------------------------
  //----------------------------------------------------------------------------------- 

  if(useAlienPlugin) {
    
    plugin = new AliAnalysisAlien();
    plugin->SetRunMode(pluginmode.Data());
    plugin->SetUser(username.Data());
    plugin->SetAPIVersion("V1.1x");
    plugin->SetAliPhysicsVersion(aliPhysicsVer.Data()); 
    plugin->SetGridWorkingDir(myWorkDir.Data()); 
    plugin->SetGridOutputDir(myOutDir.Data()); 
    plugin->SetCheckCopy(kTRUE);
    plugin->SetGridDataDir(myGridDataDir.Data());
    plugin->SetDataPattern(myDataPattern.Data());
    if (!isMC) plugin->SetRunPrefix("000");
    
    for (Int_t irun=runNmin;irun<runNmax;irun++){
      plugin->AddRunNumber((Int_t )runList[irun]);
    }
    
    //change this to se the number of files to be used in "test" mode
    plugin->SetNtestFiles(ntest); 
    
    // change this not to have the output run-by-run 
    plugin->SetNrunsPerMaster(1);
    plugin->SetOutputToRunNo(1); 
    
    // feel free to modify the list of include paths and libraries according to your need!
    plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include -I$ALICE_PHYSICS/OADB/ -I$ALICE_PHYSICS/OADB/COMMON -I$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY -g");
    //Order matters!!!!
#if 1   
    plugin->SetAdditionalLibs("liblhapdf.so libEGPythia6.so libpythia6.so libAliPythia6.so AddTaskF0V2.C ConfigF0V2.C"); 
    plugin->SetAnalysisSource("AliRsnMiniAnalysisTaskV2.cxx");
#else
    plugin->SetAdditionalLibs("liblhapdf.so libEGPythia6.so libpythia6.so libAliPythia6.so");
#endif 
    // these have to be customized 
    plugin->SetAnalysisMacro(myMacroName.Data()); 
    plugin->SetExecutable(myExecutableName.Data()); 
    plugin->SetJDLName(myJDLname);
  
    // more job option 
    plugin->SetDefaultOutputs(kTRUE);
    plugin->SetExecutableCommand("root -b -q");
    plugin->SetMaxInitFailed(15);
    plugin->SetMasterResubmitThreshold(10);
    plugin->SetTTL(50000);
    plugin->SetInputFormat("xml-single");
    if (isESD) plugin->SetSplitMaxInputFileNumber(20);
    else plugin->SetSplitMaxInputFileNumber(10);
    plugin->SetMaxMergeFiles(25);
    plugin->SetSplitMode("se");
  
    // Set properly the global flac to enable/disable local or plugin merging  
    plugin->SetMergeViaJDL(!isLocalMerge);
  
    // Force a single merging stage or define multiple stages     
    // plugin->SetOneStageMerging(kTRUE);                          
    plugin->SetMaxMergeStages(2); // adapt n to your expected number of files              
    
  } else {

    //run locally over a file chain
    if (isESD){
      chain=new TChain("esdTree");
      chain->Add("AliESDs.root");  
    } else {
      chain=new TChain("aodTree");
      chain->Add("AliAOD.root"); 
    }
    if (chain) Printf("Input chain has %d entries\n",(Int_t)chain->GetEntries());

  }
  
  //------------------------------------------------------------------------------------
  //--------------------------------- ANALYSIS SETTINGS --------------------------------
  //------------------------------------------------------------------------------------    
  /* Settings for the train  */
  Bool_t enaMultSel = kTRUE;
  Bool_t runMonOnly = kFALSE;
  Bool_t enableMon = kTRUE;
  Bool_t cutPU = kFALSE; //enable pile-up cut in physics selection
  /* settings for event mixing */
  Int_t nmix = 5; 

  //------------------------------------------------------------------------------------
  //--------------------------------- CONFIGURE TRAIN ----------------------------------
  //------------------------------------------------------------------------------------ 
  AliAnalysisManager *mgr = new AliAnalysisManager("RsnAnalysisManager");
  mgr->SetCommonFileName("RsnOut.root");
  TString out = "";
  if (isESD) {
    out = "esdTree";
    ::Info("RunGridRsnTrain", "Creating ESD handler");
    AliESDInputHandler *esdHandler = new AliESDInputHandler();
    mgr->SetInputEventHandler(esdHandler);
    esdHandler->SetNeedField(); //needed to load magnetic field configuration
    
    if (isMC) {
      ::Info("RunGridRsnTrain", "Creating MC handler");
      AliMCEventHandler *mcHandler = new AliMCEventHandler();
      mcHandler->SetReadTR(kFALSE);
      mgr->SetMCtruthEventHandler(mcHandler);
    }
  } else {
    out = "aodTree";
    ::Info("RunGridRsnTrain", "Creating AOD handler");
    AliAODInputHandler *aodHandler = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodHandler);
  }
  
  // compile the class and load the add task macro
  AliTaskCDBconnect * taskCDB = (AliTaskCDBconnect*) AddTaskCDBconnect("raw://");
  if (!taskCDB) return;
  
  //Physics selection
  AliPhysicsSelectionTask * physSelTask = (AliPhysicsSelectionTask*) AddTaskPhysicsSelection(isMC, cutPU);
  physSelTask->SelectCollisionCandidates(triggerMask);
  AliAnalysisDataContainer *cstatsout = (AliAnalysisDataContainer*)mgr->GetOutputs()->FindObject("cstatsout");
  cstatsout->SetFileName("EventStatPhysSel.root");

  //Multiplicity selection
  AliMultSelectionTask * taskMult = 0x0;
  if (enaMultSel) {
    taskMult = (AliMultSelectionTask*) AddTaskMultSelection(); 
    taskMult->SetSelectedTriggerClass((AliVEvent::EOfflineTriggerTypes)triggerMask);
    if (isMC) taskMult->SetUseDefaultMCCalib(isMC);
  }
  
  //PID response
  AliAnalysisTaskPIDResponse * pidTask = (AliAnalysisTaskPIDResponse*) AddTaskPIDResponse(isMC, isMC, isMC);
  //AliAnalysisTask * pidQAtask = 0x0;
  //if (enableMon) pidQAtask = (AliAnalysisTask*) AddTaskPIDqa();

  //Resonance task
  AliRsnMiniAnalysisTask * rsnAnaTask1 = (AliRsnMiniAnalysisTask*) AddTaskF0("f0_tpc2s_tofveto", 
                                                                             isMC, 
                                                                             AliPIDResponse::kPBPB, 
                                                                             triggerMask, 
                                                                             enaMultSel,
                                                                             eventCutSet::kEvtDefault,
                                                                             pairYCutSet::kPairDefault,
                                                                             5,
                                                                             AliRsnCutSetDaughterParticle::kTPCpidTOFveto3s,
                                                                             3.0,
                                                                             0.3,
                                                                             1.6,
                                                                             260,
                                                                             enableMon);
  //local test
  AliRsnMiniAnalysisTaskV2 * rsnAnaTaskV2  = 0x0;
  // (AliRsnMiniAnalysisTaskV2*) AddTaskF0V2("f0_devel", 
  //                                                                            isMC, 
  //                                                                            AliPIDResponse::kPBPB, 
  //                                                                            triggerMask, 
  //                                                                            enaMultSel,
  //                                                                            eventCutSet::kEvtDefault,
  //                                                                            pairYCutSet::kPairDefault,
  //                                                                            5,
  //                                                                            AliRsnCutSetDaughterParticle::kTPCpidTOFveto3s,
  //                                                                            2.0,
  //                                                                            0.3,
  //                                                                            1.6,
  //                                                                            650,
  //                                                                            enableMon);

  ::Info("AnalysisSetup", "<<<<<<<<<<<<<<<<< TRAIN CONFIGURATION >>>>>>>>>>>>>>>>");
  if (physSelTask) ::Info("AnalysisSetup", "Added physics selection");
  if (taskMult)    ::Info("AnalysisSetup", "Add centrality/multiplicity computation tasks");
  if (pidTask)     ::Info("AnalysisSetup", "Added task for PID response");
  //if (pidQAtask)   ::Info("AnalysisSetup", "Added task for PID QA");
  if (rsnAnaTask1)  ::Info("AnalysisSetup", "Added 1 task for Rsn Analysis");
  if (rsnAnaTaskV2)  ::Info("AnalysisSetup", "Added 2 task for Rsn Analysis - devel");
  ::Info("AnalysisSetup", "Setup successful"); 
  
  if (out.Length() < 1) return;
  
  // add plugin to analysis manager
  if (useAlienPlugin)
    mgr->SetGridHandler(plugin);
  
  // init analysis
  if(!mgr->InitAnalysis()) {
    printf("Error: Analysis manager failed to be initialized. Nothing done!\n");
    return;
  }
  mgr->PrintStatus();
  
  //start analysis
  mgr->StartAnalysis(analysisMode.Data(),chain,nentries,firstentry);  
  return;
}
 
 

/***************************************************************************************************/
TChain *CreateESDChain(TString esdpath, Int_t ifirst, Int_t ilast)
{
  TChain *chainESD = new TChain("esdTree");

  if(ifirst<0) {
    chainESD->Add("AliESDs.root");
  } else {
    for(Int_t i=ifirst; i<=ilast; i++) {
      TString esdfile=esdpath; esdfile+=i; esdfile.Append("/AliESDs.root");
      chainESD->Add(esdfile.Data());
    }
  }
  
  return chainESD;
}
