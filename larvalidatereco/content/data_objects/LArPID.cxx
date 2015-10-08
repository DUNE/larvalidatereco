#include "LArPID.h"

ClassImp(LArPID)

  void LArPID::Clear()
  {
    EventEventRun =        kUnassigned;
    EventID       =        kUnassigned;
    EventSubRun   =        kUnassigned;

    IsStoppingTrue=        kUnassigned;
    IsStoppingReco=        kUnassigned;
    EigenValues             .clear();
    EigenVectors            .clear();
    Covariance              .clear();
    AnglePrincipalTrueTrack .clear();
    PCAHitsSpacePoints      .clear();
    StdDevDistFromPA        .clear();
    AvgedEdxAmpStart        .clear();
    AvgedEdxAmpEnd          .clear();
    AvgedEdxAreaStart       .clear();
    AvgedEdxAreaEnd         .clear();
    EvalRatio               .clear();
    ChargeRatioCoreHalo     .clear();
    Conicalness             .clear();
    Concentration           .clear();
    ChargeLongRatio         .clear();
    ChargeLongRatioHalf     .clear();
    ChargeEndRatio          .clear();
    ChargeEndRatio10        .clear();

    NTracks=kUnassigned;

    EventEField   = (float)kUnassigned;
    EventEventRun =        kUnassigned;
    EventID       =        kUnassigned;
    EventSubRun   =        kUnassigned;
    EventT0       = (float)kUnassigned;

    NTrajMC           = kUnassigned;
    TrajNHitsMC       .clear();
    TrajChargeMC      .clear();
    TrajEnd4MomMC     .clear();
    TrajEnd4PosMC     .clear();
    TrajIDMC          .clear();
    TrajMassMC        .clear();
    TrajNameMC        .clear();
    TrajPDGMC         .clear();
    TrajParentIDMC    .clear();
    TrajStart4MomMC   .clear();
    TrajStart4PosMC   .clear();
    TrajLengthMC      .clear();

    NTracks           = kUnassigned;
    TrackCharge       .clear();
    TrackEnd4Pos      .clear();
    TrackID           .clear();
    TrackStart4Pos    .clear();
    TrackNHits        .clear();
    TrackPitch        .clear();
    TrackLength       .clear();

    NHits             = kUnassigned;
    HitChannel        .clear();
    HitCharge         .clear();
    HitClusterID      .clear();
    HitPeakT          .clear();
    HitPlane          .clear();
    HitTrackID        .clear();
    HitWire           .clear();
    HitTPC            .clear();
    Hit3Pos           .clear();
    HitIsMatched      .clear();
    //HitdEdxAmp        .clear();
    //HitdEdxArea       .clear();
  }

  LArPID::LArPID():kUnassigned(0xDEADBEEF){this->Clear();}
  LArPID::~LArPID(){}

