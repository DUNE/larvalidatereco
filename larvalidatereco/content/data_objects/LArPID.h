
#ifndef __LAR_PID__
#define __LAR_PID__
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TVector3.h"
#include "TObject.h"
#include <TLorentzVector.h>
#include <TVector3.h>

///
/// The LArAnalysis class definintion.
/// This object will constitute a single TBranch.
///
class LArPID : public TObject{

///
/// Abbreviations
///

  typedef TLorentzVector TLVec;
  typedef TVector3       T3Vec;

public:

  const int               kUnassigned;           ///< unassigned parameters initialised to this

  ///
  /// Event Bookeeping Information
  ///
  float                   EventEField;           ///< electric field strength
  int                     EventEventRun;         ///< run number
  int                     EventID;               ///< event ID
  int                     EventSubRun;           ///< sub run number
  float                   EventT0;               ///< t_0 of the event

  ///
  /// MC Truth Trajectories
  ///
  int                          NTrajMC;               ///< number of true trajectories in the event
  std::vector<int>             TrajNHitsMC;
  std::vector<int>             TrajChargeMC;          ///< true charge             [NTrajMC]
  std::vector<TLVec>           TrajEnd4MomMC;         ///< true momentum at end    [NTrajMC][px,py,pz,E]
  std::vector<TLVec>           TrajEnd4PosMC;         ///< true position at end    [NTrajMC][x,y,z,t]
  std::vector<int>             TrajIDMC;              ///< true ID                 [NTrajMC]
  std::vector<double>          TrajMassMC;            ///< true particle mass      [NTrajMC]
  std::vector<std::string>     TrajNameMC;            ///< true particle name      [NTrajMC]
  std::vector<int>             TrajPDGMC;             ///< true particle pdg code  [NTrajMC]
  std::vector<int>             TrajParentIDMC;        ///< true particle parent ID [NTrajMC]
  std::vector<TLVec>           TrajStart4MomMC;       ///< true momentum at start  [NTrajMC][px,py,pz,E]
  std::vector<TLVec>           TrajStart4PosMC;       ///< true position at start  [NTrajMC][x,y,z,t]
  std::vector<double>          TrajLengthMC;

  ///
  /// Reconstructed Tracks
  ///
  int                          NTracks;               ///< number of reconstructed tracks
  std::vector<int>             TrackCharge;           ///< reconstructed charge            [NTracks]
  std::vector<TLVec>           TrackEnd4Pos;          ///< reconstructed position at end   [NTracks][x,y,z,t]
  std::vector<int>             TrackID;               ///< reconstructed ID                [NTracks]
  std::vector<TLVec>           TrackStart4Pos;        ///< reconstructed position at start [NTracks][x,y,z,t]
  std::vector<int>             TrackNHits;            ///< number of hits along this track [NTracks]
  std::vector<float>           TrackPitch;            ///< reconstructed track pitch       [NTracks]
  std::vector<double>          TrackLength;

  ///
  /// Reconstructed Hits
  ///
  int                          NHits;                 ///< number of reconstructed hits
  std::vector<int>             HitChannel;            ///< reconstructed hit channel         [NHits]
  std::vector<float>           HitCharge;             ///< reconstructed hit charge          [NHits]
  std::vector<int>             HitClusterID;          ///< ID of cluster this hit belongs to [NHits]
  std::vector<float>           HitPeakT;              ///< reconstructed hit peak time       [NHits]
  std::vector<int>             HitPlane;              ///< reconstructed hit plane           [NHits]
  std::vector<int>             HitTrackID;            ///< ID of track this hit belongs to   [NHits]
  std::vector<int>             HitWire;               ///< reconstructed hit wire            [NHits]
  std::vector<int>             HitTPC;                ///< reconstructed hit TPC             [NHits]
  std::vector<T3Vec>           Hit3Pos;               ///< reconstructed 3D-matched position [NHits]
  std::vector<bool>            HitIsMatched;          ///< whether hit has a valid 3D position [NHits]
  //std::vector<double>          HitdEdxAmp;            ///< dE/dx using hit amplitude         [NHits]
  //std::vector<double>          HitdEdxArea;           ///< dE/dx using hit area              [NHits]
  
  std::vector<TVectorD> EigenValues;
  std::vector<TMatrixD> EigenVectors;
  std::vector<TMatrixD> Covariance;
  std::vector<TVectorD> MeanValues;
  std::vector<double>   AnglePrincipalTrueTrack;
  std::vector<std::vector<TVector3>> PCAHitsSpacePoints;
  std::vector<double>   AvgedEdxAmpStart;
  std::vector<double>   AvgedEdxAmpEnd;
  std::vector<double>   AvgedEdxAreaStart;
  std::vector<double>   AvgedEdxAreaEnd;
  std::vector<double>   EvalRatio;
  std::vector<double>   ChargeRatioCoreHalo;
  std::vector<double>   Conicalness;
  std::vector<double>   Concentration;

  ///
  /// Reset parameter values and vector containers
  ///

  virtual void Clear();

  ///
  /// Default constructor/destructor
  ///
  LArPID();
  virtual ~LArPID();
  ClassDef(LArPID,1);
};

#endif
