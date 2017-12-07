
#if !defined (__CINT__) || defined (__CLING__)
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskF0Bkg.h"
#include "AliESDtrackCuts.h"
#include "AliGenEventHeader.h"
#include "AliVHeader.h"
#include "AliPID.h"
#include "AliAnalysisCuts.h"
#include "AliAnalysisFilter.h"
#include <TString.h>
#include <TList.h>
#endif



AliAnalysisTaskF0Bkg* AddTaskF0Bkg(TString name = "name")
{
    // get the manager via the static access member. since it's static, you don't need
    // an instance of the class to call the function

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
      if (!mgr) {
       ::Error("AddTask", "No analysis manager to connect to.");
        return NULL;
      }

      UInt_t runN = mgr->GetRunFromPath();
      Printf(":::::::  PROCESSING RUN  NUMBER = %i", runN);

    // get the input event handler, again via a static method.
    // this handler is part of the managing system and feeds events
    // to your task
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTask", "This task requires an input event handler");
        return 0x0;
    }

    //TString* type = mgr->GetInputEventHandler()->GetDataType(); //"ESD" or "AOD"
    AliAnalysisFilter * fTrackFilter =  new AliAnalysisFilter("trackF");
    AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts", "myTrackCuts");
    esdTrackCuts->GetStandardITSTPCTrackCuts2011(kTRUE, 1);
    fTrackFilter->AddCuts(esdTrackCuts);

    // by default, a file is open for writing. here, we get the filename
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":F0_bkg";      // create a subfolder in the file
    // now we create an instance of your task
    AliAnalysisTaskF0Bkg* task = new AliAnalysisTaskF0Bkg(name.Data());
    if(!task) return 0x0;

    task->SelectCollisionCandidates(AliVEvent::kINT7);
    task->SetTrackFilter(fTrackFilter);
    // add your task to the manager
    mgr->AddTask(task);
    //mgr->AddTask(pidTask);
    // your task needs input: here we connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    // same for the output
    mgr->ConnectOutput(task,1,mgr->CreateContainer("MyOutputContainer", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));

    // in the end, this macro returns a pointer to your task. this will be convenient later on
    // when you will run your analysis in an analysis train on grid
    return task;
}
