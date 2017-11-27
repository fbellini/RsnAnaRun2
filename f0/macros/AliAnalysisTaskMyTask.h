/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskMyTask_H
#define AliAnalysisTaskMyTask_H

#include "AliAnalysisTaskSE.h"

class AliAODEvent;
class TList;
class TH1F;

class AliAnalysisTaskMyTask : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisTaskMyTask();
                                AliAnalysisTaskMyTask(const char *name);
        virtual                 ~AliAnalysisTaskMyTask();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

 private:
	void                    FillDaughterPhaseSpaceHisto(AliStack * stack, TParticle * mother);
	AliVVertex*             fPrimaryVertex;//!<! Primary vertex pointer

        TList*                  fOutputList;    //! output list
        TH2F*                   fHistF0Reco;
        TH2F*                   fHistF0Gen;
	TH2F*                   fHistF0daughtersPhaseSpace;



        AliAnalysisTaskMyTask(const AliAnalysisTaskMyTask&); // not implemented
        AliAnalysisTaskMyTask& operator=(const AliAnalysisTaskMyTask&); // not implemented

        ClassDef(AliAnalysisTaskMyTask, 1);
};

#endif