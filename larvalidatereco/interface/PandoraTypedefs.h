#include "LArPandoraInterface/LArPandoraCollector.h"

namespace lar_valrec
{
  //Import shorthands for art containers used by Pandora
  //These are defined in LArPandoraCollector.h
using lar_pandora::WireVector;
using lar_pandora::HitVector;
using lar_pandora::SpacePointVector;
using lar_pandora::ClusterVector;
using lar_pandora::SeedVector;
using lar_pandora::VertexVector;
using lar_pandora::TrackVector;
using lar_pandora::ShowerVector;
using lar_pandora::PFParticleVector;
using lar_pandora::MCTruthVector;
using lar_pandora::MCParticleVector;
using lar_pandora::SimChannelVector;
using lar_pandora::TrackIDEVector;

using lar_pandora::PFParticlesToTracks;
using lar_pandora::PFParticlesToShowers;
using lar_pandora::PFParticlesToClusters;
using lar_pandora::PFParticlesToSeeds;
using lar_pandora::PFParticlesToVertices;
using lar_pandora::PFParticlesToSpacePoints;
using lar_pandora::PFParticlesToHits;
using lar_pandora::ClustersToHits;
using lar_pandora::SpacePointsToHits;
using lar_pandora::MCTruthToMCParticles;
using lar_pandora::MCParticlesToMCTruth;
using lar_pandora::MCParticlesToHits;
using lar_pandora::MCParticlesToPFParticles;
using lar_pandora::HitsToPFParticles;
using lar_pandora::HitsToMCParticles;
using lar_pandora::HitsToTrackIDEs;

using lar_pandora::PFParticleMap;
using lar_pandora::ClusterMap;
using lar_pandora::SpacePointMap;
using lar_pandora::HitMap;
using lar_pandora::MCParticleMap;

//Some new relationships
typedef std::map< art::Ptr<recob::Hit>, art::Ptr<recob::Track> >              HitsToTracks;
typedef std::map< art::Ptr<recob::Track>, std::vector<art::Ptr<recob::Hit> > >  TracksToHits;
typedef std::map< art::Ptr<recob::Hit>,   art::Ptr<recob::Cluster> >            HitsToClusters;
}
