#include "LArAnalysis.h"

ClassImp(LArAnalysis)

  void LArAnalysis::Clear()
  {
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

    NVtxMC            = kUnassigned;
    Vtx4PosMC         .clear();
    VtxNPrimTrajMC    .clear();
    VtxPrimTrajIDsMC  .clear();
    VtxReactionCodeMC .clear();
    VtxTargetPDGMC    .clear();
    VtxNuMomMC        .clear();
    VtxNuPDGMC        .clear();

    NTracks           = kUnassigned;
    TrackCharge       .clear();
    TrackEnd4Pos      .clear();
    TrackID           .clear();
    TrackStart4Pos    .clear();
    TrackNHits        .clear();
    TrackPitch        .clear();
    TrackLength       .clear();

    NClusters         = kUnassigned;
    ClusterCharge     .clear();
    ClusterID         .clear();
    ClusterNHits      .clear();
    ClusterTrackID    .clear();

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
  }

  LArAnalysis::LArAnalysis():kUnassigned(0xDEADBEEF){this->Clear();}
  LArAnalysis::~LArAnalysis(){}

