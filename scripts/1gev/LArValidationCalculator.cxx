#include "TDatabasePDG.h"

#include "larvalidatereco/framework/VarHelper.h"
#include "larvalidatereco/interface/EventHelper.h"
#include "LArValidationCalculator.h"

#include "TGraph.h"
#include "TFile.h"

namespace lar_valrec{

  class LArValidation_register{
  public:
    LArValidation_register(){
        VarHelper::RegisterOutput<LArValidation,LArValidationCalculator>("Validation");
    }
  };
  static LArValidation_register __LArValidationRegister;


void LArValidationCalculator::Calculate(TObject* tObjectPtr,const VarHelper& varHelper,const EventHelper& evHelper){
  LArValidation* outputPtr=dynamic_cast<LArValidation*>(tObjectPtr);
  if(!outputPtr){
    exit(1);
  }
  outputPtr->Clear();

  const TrackVector& tracks=evHelper.GetTracks();
  const TracksToHits& tracksToHits=evHelper.GetTracksToHits();
  const MCParticlesToHits& particlesToHits=evHelper.GetMCParticleToHitAssociations();
  const HitsToMCParticles& hitsToParticles=evHelper.GetHitToMCParticleAssociations();

  outputPtr->NTracks=0;
  for(auto track=tracks.begin();track!=tracks.end();++track){
    if(!tracksToHits.count(*track)) continue;
    ++outputPtr->NTracks;
    outputPtr->TrackID.push_back((*track)->ID());
    std::map<art::Ptr<simb::MCParticle>,double> particleCharges;
    std::map<int,bool> trajTPCs;
    std::map<int,bool> trackTPCs;
    double trackCharge=0;

    for(auto hit=tracksToHits.at(*track).begin();hit!=tracksToHits.at(*track).end();++hit){
      if(hitsToParticles.count(*hit)&&!hitsToParticles.at(*hit).isNull()){
	particleCharges[hitsToParticles.at(*hit)]+=(*hit)->Integral();
	trackTPCs[(*hit)->WireID().TPC]=true;
      }
      trackCharge+=(*hit)->Integral();
    }

    const art::Ptr<simb::MCParticle>* bestParticle=0;
    double maxCharge=0.;

    for(auto particle=particleCharges.begin();particle!=particleCharges.end();++particle){
      if(!bestParticle||particle->second>maxCharge){
	bestParticle=&(particle->first);
	maxCharge=particle->second;
      }
    }
    if(bestParticle){
      outputPtr->TrackTrajID.push_back((*bestParticle)->TrackId());
      
      double trajCharge=0.;
      for(auto hIter=particlesToHits.at(*bestParticle).begin();hIter!=particlesToHits.at(*bestParticle).end();++hIter){
	trajCharge+=(*hIter)->Integral();
	trajTPCs[(*hIter)->WireID().TPC]=true;
      }

      //Stitching has worked if the track contains all the TPCs included in the trajectory
      bool stitched=-1;
      if(trajTPCs.size()!=1){
	stitched=1;
	for(auto tpcIter=trajTPCs.begin();tpcIter!=trajTPCs.end();++tpcIter){
	  if(!trackTPCs.count(tpcIter->first)){
	    stitched=0;
	  }
	}
      }
      outputPtr->TrackStitched.push_back(stitched);
      outputPtr->TrackCompleteness.push_back(maxCharge/trajCharge);
      outputPtr->TrackPurity.push_back(maxCharge/trackCharge);
    }
    else{
      outputPtr->TrackStitched.push_back(-1);
      outputPtr->TrackTrajID.push_back(-1);
      outputPtr->TrackCompleteness.push_back(-1.);
      outputPtr->TrackPurity.push_back(-1);
    }
  }
}

}
