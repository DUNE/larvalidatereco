#include "LArPandoraEventHelper.h"
#include "LArPandoraInterface/LArPandoraCollector.h"
#include "LArPandoraInterface/LArPandoraHelper.h"

namespace lar_valrec{

  LArPandoraEventHelper::LArPandoraEventHelper(const fhicl::ParameterSet& pset):
    EventHelper(pset),
    fSimModuleName(pset.get<std::string>("SimModuleName")),
    fRecoModuleName(pset.get<std::string>("RecoModuleName")),
    fHitModuleName(pset.get<std::string>("HitModuleName")),
    fStitcherModuleName(pset.get<std::string>("StitcherModuleName"))
{
  }

  void LArPandoraEventHelper::LoadEvent_(const art::Event& evt){

    fHits.clear();
    fMCParticles.clear();
    fMCVertices.clear();
    fPFParticles.clear();
    fClusters.clear();
    fTracks.clear();

    //Relationships
    fMCParticlesToHits.clear();
    fHitsToMCParticles.clear();
    fTruthToParticles.clear();
    fParticlesToTruth.clear();
    fPFParticlesToHits.clear();
    fHitsToPFParticles.clear();
    fClustersToHits.clear();
    fHitsToClusters.clear();
    fHitsToTracks.clear();
    fTracksToHits.clear();
    
    //Collect simulation and reconstruction objects
    lar_pandora::LArPandoraCollector::CollectMCParticles(evt, fSimModuleName, fMCParticles);
    lar_pandora::LArPandoraCollector::CollectMCParticles(evt, fSimModuleName, fTruthToParticles,
							 fParticlesToTruth);
    lar_pandora::LArPandoraCollector::CollectHits(evt, fHitModuleName, fHits);
    lar_pandora::LArPandoraCollector::CollectClusters(evt, fRecoModuleName, fClusters, fClustersToHits);
    lar_pandora::LArPandoraCollector::CollectTracks(evt,fStitcherModuleName, fTracks, fTracksToHits);

    //Map hits onto MC particles and pandora PFParticles.
    lar_pandora::LArPandoraCollector::BuildMCParticleHitMaps(evt,fSimModuleName,fHits,fMCParticlesToHits, 
							     fHitsToMCParticles,lar_pandora::LArPandoraCollector::kIgnoreDaughters);
    
    //Build MCTruth vertex vector from vertex->particle map
    for(auto vIter=fTruthToParticles.begin();vIter!=fTruthToParticles.end();++vIter){
      fMCVertices.push_back(vIter->first);
    }

    //Reverse clusters->hits mapping
    fHitsToClusters.clear();
    for(auto clIter=fClustersToHits.begin();clIter!=fClustersToHits.end();++clIter){
      for(auto hitIter=clIter->second.begin();hitIter!=clIter->second.end();++hitIter){
	fHitsToClusters[*hitIter]=clIter->first;
      }
    }

    //Reverse tracks->hits mapping
    fHitsToTracks.clear();
    for(auto trIter=fTracksToHits.begin();trIter!=fTracksToHits.end();++trIter){
      for(auto hitIter=trIter->second.begin();hitIter!=trIter->second.end();++hitIter){
	fHitsToTracks[*hitIter]=trIter->first;
      }
    }    
  }

  LArPandoraEventHelper::~LArPandoraEventHelper(){
  }

}//namespace lar_valrec
