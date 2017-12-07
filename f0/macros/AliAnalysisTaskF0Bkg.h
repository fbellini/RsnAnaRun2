/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskF0Bkg_H
#define AliAnalysisTaskF0Bkg_H

#include "AliAnalysisTaskSE.h"
#include "AliPIDResponse.h"
#include "AliPID.h"

class AliAnalysisFilter;
//class AliAnalysisTaskPIDResponse;
class AliAnalysisTaskF0Bkg : public AliAnalysisTaskSE

{
    public:
                                AliAnalysisTaskF0Bkg();
                                AliAnalysisTaskF0Bkg(const char *name);
        virtual                 ~AliAnalysisTaskF0Bkg();
        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);
        virtual void            SetTrackFilter(AliAnalysisFilter* trackF) {fTrackFilter = trackF;}
        Float_t                 PairPt(AliESDtrack* track1, AliESDtrack* track2);
        Float_t                 PairEta(AliESDtrack* track1, AliESDtrack* track2);
        Float_t                 PairY(AliESDtrack* track1, AliESDtrack* track2);


    private:
        AliESDEvent*            fESD;
        AliMCEvent*             fMCEvent;
        AliStack*               fMCStack;
        AliAnalysisFilter*      fTrackFilter;
        //AliESDpid*              fESDpid;
        AliPIDResponse*         fPID;
        TList*                  fOutputList;
        TH1I*                   fNEvents;
        TH2F*                   fHistPtGen[10];
        TH2F*                   fHistPtReco[10];
        TH2F*                   fHistYGen[10];
        TH2F*                   fHistYReco[10];
        TH2F*                   fHistEtaGen[10];
        TH2F*                   fHistEtaReco[10];
        const static Char_t     fParticleName[][6];
        const static ULong_t    fPdgArray[];


        AliAnalysisTaskF0Bkg(const AliAnalysisTaskF0Bkg&); // not implemented
        AliAnalysisTaskF0Bkg& operator=(const AliAnalysisTaskF0Bkg&); // not implemented

        ClassDef(AliAnalysisTaskF0Bkg, 1);
};

#endif
