#include "LArPID.h"

ClassImp(LArPID)

  void LArPID::Clear()
  {
    EventEventRun =        kUnassigned;
    EventID       =        kUnassigned;
    EventSubRun   =        kUnassigned;

    IsStoppingTrue=        kUnassigned;
    IsStoppingReco=        kUnassigned;
    EigenValues              .clear();
    EigenVectors             .clear();
    Covariance               .clear();
    AnglePATrueTrack         .clear();
    PCAHitsSpacePoints       .clear();
    StdDevDistFromFitLine    .clear();
    AvgedEdxAmpStart         .clear();
    AvgedEdxAmpEnd           .clear();
    AvgedEdxAmpLongRatio     .clear();
    AvgedEdxAmpEndRatio      .clear();
    AvgedEdxAmpEndRatio10    .clear();
    AvgedEdxAmpStartDist     .clear();
    AvgedEdxAmpEndDist       .clear();
    AvgedEdxAmpLongRatioDist .clear();
    AvgedEdxAmpEndRatioDist  .clear();
    AvgedEdxAreaStart        .clear();
    AvgedEdxAreaEnd          .clear();
    AvgedEdxAreaLongRatio    .clear();
    AvgedEdxAreaEndRatio     .clear();
    AvgedEdxAreaEndRatio10   .clear();
    AvgedEdxAreaStartDist    .clear();
    AvgedEdxAreaEndDist      .clear();
    AvgedEdxAreaLongRatioDist.clear();
    AvgedEdxAreaEndRatioDist .clear();
    EvalRatio                .clear();
    ChargeRatioCoreHalo      .clear();
    Conicalness              .clear();
    Concentration            .clear();
    ChargeLongRatio          .clear();
    ChargeLongRatioHalf      .clear();
    ChargeEndRatio           .clear();
    ChargeEndRatio10         .clear();

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

