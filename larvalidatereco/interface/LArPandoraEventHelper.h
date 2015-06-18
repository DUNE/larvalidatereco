/**
 *  @file   larvalidatereco/interface/LArPandoraEventHelper.h
 *
 *  @brief  Helper to load in recon quantities from Pandora reconstruction
 *
 */

#ifndef LARPANDORAEVENTHELPER_H
#define LARPANDORAEVENTHELPER_H 1

// Framework Includes
#include "EventHelper.h"
#include "LArPandoraInterface/LArPandoraCollector.h"

namespace lar_valrec
{

/**
 *  @brief  LArPandoraEventHelper class
 */
  class LArPandoraEventHelper : public EventHelper
{
public:
  LArPandoraEventHelper(const fhicl::ParameterSet& pset);
  ~LArPandoraEventHelper();

  //Getters for all items pulled out of event
  virtual const MCParticleVector& GetMCParticles() const;
  virtual const MCTruthVector& GetMCVertices() const;
  virtual const HitVector& GetHits() const;
  virtual const SpacePointVector& GetSpacePoints() const;
  virtual const MCTruthToMCParticles& GetTruthToParticles() const;
  virtual const MCParticlesToMCTruth& GetParticlesToTruth() const;
  virtual const MCParticlesToHits& GetMCParticleToHitAssociations() const;
  virtual const HitsToMCParticles& GetHitToMCParticleAssociations() const;
  virtual const PFParticleVector& GetPFParticles() const;
  virtual const PFParticlesToHits& GetPFParticleToHitAssociations() const;
  virtual const HitsToPFParticles& GetHitToPFParticleAssociations() const;
  virtual const ClusterVector& GetClusters() const;
  virtual const TrackVector& GetTracks() const;
  virtual const ClustersToHits& GetClustersToHits() const;
  virtual const HitsToClusters& GetHitsToClusters() const;
  virtual const TracksToHits& GetTracksToHits() const;
  virtual const HitsToTracks& GetHitsToTracks() const;
  virtual const HitsToSpacePoints& GetHitsToSpacePoints() const;
  virtual const SpacePointsToHits& GetSpacePointsToHits() const;

 private: 

    ///Collect objects from event
    ///Called by base class public LoadEvent method.
    virtual void LoadEvent_(const art::Event& evt);

    //Input data products
    HitVector fHits;
    MCParticleVector fMCParticles;
    MCTruthVector fMCVertices;
    PFParticleVector fPFParticles;
    ClusterVector fClusters;
    TrackVector fTracks;
    SpacePointVector fSpacePoints;

    //Relationships
    MCTruthToMCParticles fTruthToParticles;
    MCParticlesToMCTruth fParticlesToTruth;
    MCParticlesToHits fMCParticlesToHits;
    HitsToMCParticles fHitsToMCParticles;
    PFParticlesToHits fPFParticlesToHits;
    HitsToPFParticles fHitsToPFParticles;
    ClustersToHits fClustersToHits;
    HitsToClusters fHitsToClusters;
    HitsToTracks fHitsToTracks;
    TracksToHits fTracksToHits;
    SpacePointsToHits fSpacePointsToHits;
    HitsToSpacePoints fHitsToSpacePoints;

    std::string fSimModuleName;
    std::string fRecoModuleName;
    std::string fHitModuleName;
    std::string fStitcherModuleName;
};

  inline const MCParticleVector& LArPandoraEventHelper::GetMCParticles() const{
    return fMCParticles;
  }

  inline const MCTruthVector& LArPandoraEventHelper::GetMCVertices() const{
    return fMCVertices;
  }

  inline const MCTruthToMCParticles& LArPandoraEventHelper::GetTruthToParticles() const{
    return fTruthToParticles;
  }

  inline const MCParticlesToMCTruth& LArPandoraEventHelper::GetParticlesToTruth() const{
    return fParticlesToTruth;
  }

  inline const HitVector& LArPandoraEventHelper::GetHits() const{
    return fHits;
  }

  inline const PFParticleVector& LArPandoraEventHelper::GetPFParticles() const{
    return fPFParticles;
  }

  inline const ClusterVector& LArPandoraEventHelper::GetClusters() const{
    return fClusters;
  }

  inline const TrackVector& LArPandoraEventHelper::GetTracks() const{
    return fTracks;
  }

  inline const MCParticlesToHits& LArPandoraEventHelper::GetMCParticleToHitAssociations() const{
    return fMCParticlesToHits;
  }

  inline const HitsToMCParticles& LArPandoraEventHelper::GetHitToMCParticleAssociations() const{
    return fHitsToMCParticles;
  }

  inline const PFParticlesToHits& LArPandoraEventHelper::GetPFParticleToHitAssociations() const{
    return fPFParticlesToHits;
  }

  inline const HitsToPFParticles& LArPandoraEventHelper::GetHitToPFParticleAssociations() const{
    return fHitsToPFParticles;
  }

  inline const ClustersToHits& LArPandoraEventHelper::GetClustersToHits() const{
    return fClustersToHits;
  }

  inline const HitsToClusters& LArPandoraEventHelper::GetHitsToClusters() const{
    return fHitsToClusters;
  }

  inline const TracksToHits& LArPandoraEventHelper::GetTracksToHits() const{
    return fTracksToHits;
  }

  inline const HitsToTracks& LArPandoraEventHelper::GetHitsToTracks() const{
    return fHitsToTracks;
  }

  inline const SpacePointVector& LArPandoraEventHelper::GetSpacePoints() const{
    return fSpacePoints;
  }

  inline const SpacePointsToHits& LArPandoraEventHelper::GetSpacePointsToHits() const{
    return fSpacePointsToHits;
  }

  inline const HitsToSpacePoints& LArPandoraEventHelper::GetHitsToSpacePoints() const{
    return fHitsToSpacePoints;
  }


} // namespace lar_valrec

#endif // #ifndef LARPANDORAEVENTHELPER_H
