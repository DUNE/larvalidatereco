/**
 *  @file   larvalidatereco/interface/EventHelper.h
 *
 *  @brief  Base for classes pulling out reco objects from different reconstruction tools
 *
 */

#ifndef EVENTHELPER_H
#define EVENTHELPER_H 1

#include "PandoraTypedefs.h"
#include "LArValidateRecon.h"
#include "larvalidatereco/framework/VarHelper.h"

namespace lar_valrec
{

/**
 *  @brief  EventHelper class
 */
class EventHelper
{
public:
  EventHelper(const fhicl::ParameterSet& pset){}

  ///Get event metadata and call LoadEvent_ method implemented by derived class
  virtual void LoadEvent(const art::Event& evt);
  
  //abstract methods to get reconstructed quantities from event.
  //Implemented in derived classes
  virtual const MCParticleVector& GetMCParticles() const=0;
  virtual const MCTruthVector& GetMCVertices() const=0;
  virtual const MCTruthToMCParticles& GetTruthToParticles() const=0;
  virtual const MCParticlesToMCTruth& GetParticlesToTruth() const=0;
  virtual const HitVector& GetHits() const=0;
  virtual const SpacePointVector& GetSpacePoints() const=0;
  virtual const MCParticlesToHits& GetMCParticleToHitAssociations() const=0;
  virtual const HitsToMCParticles& GetHitToMCParticleAssociations() const=0;
  virtual const PFParticleVector& GetPFParticles() const=0;
  virtual const PFParticlesToHits& GetPFParticleToHitAssociations() const=0;
  virtual const HitsToPFParticles& GetHitToPFParticleAssociations() const=0;
  virtual const ClusterVector& GetClusters() const=0;
  virtual const TrackVector& GetTracks() const=0;
  virtual const ClustersToHits& GetClustersToHits() const=0;
  virtual const HitsToClusters& GetHitsToClusters() const=0;
  virtual const TracksToHits& GetTracksToHits() const=0;
  virtual const HitsToTracks& GetHitsToTracks() const=0;
  virtual const HitsToSpacePoints& GetHitsToSpacePoints() const=0;
  virtual const SpacePointsToHits& GetSpacePointsToHits() const=0;

  //Getters for event metadata
  virtual int GetRun() const{return fRun;}
  virtual int GetSubRun() const{return fSubRun;}
  virtual int GetEventID() const{return fEventID;}
  ///Magnitude of E field only for now
  virtual double GetEField() const{return fEField;}
  virtual int GetEventT0() const{return fEventT0;}
  
  virtual ~EventHelper(){}

private:

  //Abstract method; must be overriden in derived classes by
  //code which actually pulls data from event
  virtual void LoadEvent_(const art::Event& evt)=0;

  int     fRun;         ///< run number
  int     fSubRun;      ///< sub run number
  int     fEventID;     ///< event ID
  double  fEField;      ///< electric field strength
  int     fEventT0;     ///< t_0 of the event

};

} // namespace lar_valrec

#endif // #ifndef EVENTHELPER_H