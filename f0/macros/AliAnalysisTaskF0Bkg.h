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
        TH2F*                   fHistBeforePID1tpc;
        TH2F*                   fHistBeforePID2tpc;
        TH2F*                   fHistBeforePID1tof;
        TH2F*                   fHistBeforePID2tof;
        TH2F*                   fHistPID1tpc;
        TH2F*                   fHistPID2tpc;
        TH2F*                   fHistPID1tof;
        TH2F*                   fHistPID2tof;
        TH2F*                   fNoSelTracksEtaPt;
        TH2F*                   fAcceptedGeomAccEtaPt;
        TH2F*                   fAcceptedQualityEtaPt;
        TH2F*                   fAcceptedTOFMatchEtaPt;
        TH2F*                   fAcceptedTracksEtaPt;
        TH2F*                   fAcceptedPIDTOFMatchEtaPt;
        TH2F*                   fAcceptedPIDnoTOFMatchEtaPt;
        TH2F*                   fGenMassVsPt[9];
        TH2F*                   fRecoMassVsPt[9];
        TH2F*                   fGenYVsPt[9];
        TH2F*                   fRecoYVsPt[9];
        TH2F*                   fGenEtaVsPt[9];
        TH2F*                   fRecoEtaVsPt[9];
        const static Char_t     fParticleName[][6];
        const static ULong_t    fPdgArray[];

        AliAnalysisTaskF0Bkg(const AliAnalysisTaskF0Bkg&); // not implemented
        AliAnalysisTaskF0Bkg& operator=(const AliAnalysisTaskF0Bkg&); // not implemented

        ClassDef(AliAnalysisTaskF0Bkg, 1);
};

#endif
