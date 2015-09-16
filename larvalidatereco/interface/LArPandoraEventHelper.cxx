#include "LArPandoraEventHelper.h"
#include "LArPandoraInterface/LArPandoraCollector.h"
#include "LArPandoraInterface/LArPandoraHelper.h"
#include "art/Framework/Core/FindOneP.h"

namespace lar_valrec{

  LArPandoraEventHelper::LArPandoraEventHelper(const fhicl::ParameterSet& pset):
    EventHelper(pset),
    fSimModuleName(pset.get<std::string>("SimModuleName")),
    fRecoModuleName(pset.get<std::string>("RecoModuleName")),
    fHitModuleName(pset.get<std::string>("HitModuleName")),
    fStitcherModuleName(pset.get<std::string>("StitcherModuleName")),
    fCaloModuleName(pset.get<std::string>("CaloModuleName"))
{
  }

  void LArPandoraEventHelper::CollectCalo(const art::Event &evt, const std::string label, CaloVector &caloVector, 
					  TracksToCalo &tracksToCalo)
{
  art::Handle< std::vector<anab::Calorimetry> > caloObjs;
  evt.getByLabel(label, caloObjs);

  if (!caloObjs.isValid())
    return;

  art::FindOneP<recob::Track> caloTrackAssns(caloObjs, evt, label);
  for (unsigned int i = 0; i < caloObjs->size(); ++i)
    {
      const art::Ptr<anab::Calorimetry> calo(caloObjs, i);
      caloVector.push_back(calo);
      const art::Ptr<recob::Track> track = caloTrackAssns.at(i);
      tracksToCalo[track] = calo;
    }
	 }	 
  

  void LArPandoraEventHelper::LoadEvent_(const art::Event& evt){

    fHits.clear();
    fMCParticles.clear();
    fMCVertices.clear();
    fClusters.clear();
    fTracks.clear();
    fSpacePoints.clear();
    fShowers.clear();
    fCalo.clear();
    //Relationships
    fMCParticlesToHits.clear();
    fHitsToMCParticles.clear();
    fTruthToParticles.clear();
    fParticlesToTruth.clear();
    fClustersToHits.clear();
    fHitsToClusters.clear();
    fHitsToTracks.clear();
    fTracksToHits.clear();
    fShowersToHits.clear();
    fSpacePointsToHits.clear();
    fHitsToSpacePoints.clear();
    fTracksToCalo.clear();

    //Collect simulation and reconstruction objects
    CollectCalo(evt,fCaloModuleName,fCalo,fTracksToCalo);

    for(auto iter=fTracksToCalo.begin();iter!=fTracksToCalo.end();++iter){
      std::cout<<"Calo for track "<<iter->first->ID()<<std::endl;
      art::Ptr<anab::Calorimetry> calo=iter->second;

      std::cout<<"Range "<<calo->Range()<<std::endl;
      std::cout<<"KE "<<calo->KineticEnergy()<<std::endl;
    }

    lar_pandora::LArPandoraCollector::CollectMCParticles(evt, fSimModuleName, fMCParticles);
    lar_pandora::LArPandoraCollector::CollectMCParticles(evt, fSimModuleName, fTruthToParticles,
							 fParticlesToTruth);
    lar_pandora::LArPandoraCollector::CollectHits(evt, fHitModuleName, fHits);
    lar_pandora::LArPandoraCollector::CollectClusters(evt, fRecoModuleName, fClusters, fClustersToHits);
    lar_pandora::LArPandoraCollector::CollectTracks(evt,fStitcherModuleName, fTracks, fTracksToHits);
    lar_pandora::LArPandoraCollector::CollectSpacePoints(evt,fRecoModuleName,fSpacePoints, 
							 fSpacePointsToHits, fHitsToSpacePoints);
    lar_pandora::LArPandoraCollector::CollectShowers(evt,fRecoModuleName,fShowers,fShowersToHits);

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

    //Reverse showers->hits mapping
    fHitsToShowers.clear();
    for(auto shIter=fShowersToHits.begin();shIter!=fShowersToHits.end();++shIter){
      for(auto hitIter=shIter->second.begin();hitIter!=shIter->second.end();++hitIter){
	fHitsToShowers[*hitIter]=shIter->first;
      }
    }    
  }

  LArPandoraEventHelper::~LArPandoraEventHelper(){
  }

}//namespace lar_valrec
