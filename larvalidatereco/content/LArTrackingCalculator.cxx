#include "TDatabasePDG.h"

#include "Geometry/TPCGeo.h"
#include "Geometry/Geometry.h"
#include "larvalidatereco/framework/VarHelper.h"
#include "larvalidatereco/interface/EventHelper.h"
#include "LArTrackingCalculator.h"

#include "TGraph.h"
#include "TFile.h"

namespace lar_valrec{

  class LArTracking_register{
  public:
    LArTracking_register(){
        VarHelper::RegisterOutput<LArTracking,LArTrackingCalculator>("Tracking");
    }
  };
  static LArTracking_register __LArTrackingRegister;


void LArTrackingCalculator::Calculate(TObject* tObjectPtr,const VarHelper& varHelper,const EventHelper& evHelper){
  LArTracking* outputPtr=dynamic_cast<LArTracking*>(tObjectPtr);
  if(!outputPtr){
    exit(1);
  }
  outputPtr->Clear();

  outputPtr->NTrajMC=0;
  const MCParticleVector& mcParticles=evHelper.GetMCParticles();
  for(auto particle=mcParticles.begin();particle!=mcParticles.end();++particle){
    ++(outputPtr->NTrajMC);
    
    outputPtr->TrajIDMC.push_back((*particle)->TrackId());

    std::vector<TVector3> trackPoints;

    for(unsigned int point=0;point!=(*particle)->NumberTrajectoryPoints();++point){
      trackPoints.push_back((*particle)->Position(point).Vect());
    }
    outputPtr->TrajPointsMC.push_back(trackPoints);
  }

  const TrackVector& tracks=evHelper.GetTracks();

  outputPtr->NTracks=0;
  for(auto track=tracks.begin();track!=tracks.end();++track){
    ++(outputPtr->NTracks);
      outputPtr->TrackID.push_back((*track)->ID());

      std::vector<TVector3> trackPoints;
      
      for(unsigned int point=0;point!=(*track)->NumberTrajectoryPoints();++point){
	trackPoints.push_back((*track)->LocationAtPoint(point));
      }
      outputPtr->TrackPoints.push_back(trackPoints);
  }
 
}
}
