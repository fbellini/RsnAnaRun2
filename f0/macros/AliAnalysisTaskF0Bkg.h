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
	Bool_t                  IsGoodSPDvertexRes(const AliESDVertex * spdVertex);
	Bool_t                  SelectVertex2015pp(AliVEvent *event, Bool_t checkSPDres = kTRUE, Bool_t requireSPDandTrk = kFALSE, Bool_t checkProximity = kTRUE, Bool_t enaMonitor = kTRUE);
        Bool_t                  IsTOFMatched(AliESDtrack *track);

    private:
        AliESDEvent*            fESD;
        AliMCEvent*             fMCEvent;
        AliStack*               fMCStack;
        AliAnalysisFilter*      fTrackFilter;
        //AliESDpid*              fESDpid;
        AliPIDResponse*         fPIDResponse;
        TList*                  fOutputList;
        TH1I*                   fNEvents;
        TH2F*                   fHistPID1tpc;
        TH2F*                   fHistPID2tpc;
        TH2F*                   fHistPID1tof;
        TH2F*                   fHistPID2tof;
        TH2F*                   fHistoNoSelTracksEtaPt1;
        TH2F*                   fHistoAcceptedTracksEtaPt1;
        TH2F*                   fHistoAcceptedTracksEtaPt2;
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
