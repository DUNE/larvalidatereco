/**
 *  @file   larvalidaterecon/interface/LArValidateReco_module.cc
 *
 *  @brief  Provide proper module filename for LArValidateRecon class
 *
 */

// Framework Includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "PandoraTypedefs.h"
#include "LArPandoraInterface/LArPandoraHelper.h"
#include "Utilities/AssociationUtil.h"
// Local includes

// std includes
#include <string>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_valrec
{

/**
 *  @brief  LArValRec class
 */
class TracksFromPFOs : public art::EDProducer
{
public: 
  TracksFromPFOs(fhicl::ParameterSet const &pset);
  
  virtual ~TracksFromPFOs();
  
  virtual void beginJob();
  
  virtual void produce(art::Event &evt);

private:

  std::string fStitcherModuleName;
  std::string fRecoModuleName;
};

DEFINE_ART_MODULE(TracksFromPFOs)

} // namespace lar_valrec

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows


namespace lar_valrec {

  TracksFromPFOs::TracksFromPFOs(fhicl::ParameterSet const &pset) :
    fStitcherModuleName(pset.get<std::string>("StitcherModuleName")),
    fRecoModuleName(pset.get<std::string>("RecoModuleName"))
  {
    produces< std::vector<recob::Track> >(); 
    produces< art::Assns<recob::Track, recob::Hit> >();
  }

//------------------------------------------------------------------------------------------------------------------------------------------

TracksFromPFOs::~TracksFromPFOs()
{
}

void TracksFromPFOs::beginJob()
{
}

void TracksFromPFOs::produce(art::Event& evt)
{
  std::unique_ptr< std::vector<recob::Track> > outputTracks( new std::vector<recob::Track> );
  std::unique_ptr< art::Assns<recob::Track, recob::Hit> > outputTracksToHits( new art::Assns<recob::Track, recob::Hit> );
  PFParticlesToHits pfParticlesToHits;
  HitsToPFParticles hitsToPFParticles;
  SpacePointVector spacePointVector;
  SpacePointsToHits spacePointsToHits;
  HitsToSpacePoints hitsToSpacePoints;

  lar_pandora::LArPandoraCollector::BuildPFParticleHitMaps(evt, fStitcherModuleName, fRecoModuleName,
							   pfParticlesToHits, hitsToPFParticles,lar_pandora::LArPandoraCollector::kIgnoreDaughters);
  lar_pandora::LArPandoraCollector::CollectSpacePoints(evt,fRecoModuleName,spacePointVector,spacePointsToHits,hitsToSpacePoints);   

  int trackCounter=0;
  std::cout<<"Input PFOs: "<<pfParticlesToHits.size()<<std::endl;
  for(auto pfoMapIter=pfParticlesToHits.begin();pfoMapIter!=pfParticlesToHits.end();++pfoMapIter){
    try{
      HitVector& pfoHits=pfoMapIter->second;
      SpacePointVector pfoSpacePoints;
      std::cout<<"Size of hit vector: "<<pfoHits.size()<<std::endl;
      for(auto hitIter=pfoHits.begin();hitIter!=pfoHits.end();++hitIter){
	if(hitsToSpacePoints[*hitIter].isNonnull()){
	  pfoSpacePoints.push_back(hitsToSpacePoints[*hitIter]);
	}
      }
      std::cout<<"Size of sp vector: "<<pfoSpacePoints.size()<<std::endl;
      recob::Track newTrack(lar_pandora::LArPandoraHelper::BuildTrack(trackCounter++, pfoSpacePoints,false));    
      outputTracks->push_back(newTrack);
      util::CreateAssn(*this, evt, *(outputTracks.get()), pfoHits, *(outputTracksToHits.get()));
    }
    catch (cet::exception &e)
      {
	std::cout<<"Exception thrown!"<<std::endl;
      }
    
  }

  evt.put(std::move(outputTracks));
  evt.put(std::move(outputTracksToHits));
  
}

} // namespace lar_valrec
