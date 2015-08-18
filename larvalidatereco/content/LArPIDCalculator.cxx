#include "TDatabasePDG.h"

#include "Geometry/TPCGeo.h"
#include "Geometry/Geometry.h"
#include "larvalidatereco/framework/VarHelper.h"
#include "larvalidatereco/interface/EventHelper.h"
#include "LArPIDCalculator.h"

#include "TGraph.h"
#include "TFile.h"

namespace lar_valrec{

  class LArPID_register{
  public:
    LArPID_register(){
        VarHelper::RegisterOutput<LArPID,LArPIDCalculator>("PID");
    }
  };
  static LArPID_register __LArPIDRegister;

  void LArPIDCalculator::Calculate(TObject* tObjectPtr,const VarHelper& varHelper,const EventHelper& evHelper){
  LArPID* outputPtr=dynamic_cast<LArPID*>(tObjectPtr);
  if(!outputPtr){
    exit(1);
  }
  outputPtr->Clear();
  /*
  const TrackVector& tracks=evHelper.GetTracks();
  const TracksToHits& tracksToHits=evHelper.GetTracksToHits();
  const MCParticlesToHits& particlesToHits=evHelper.GetMCParticleToHitAssociations();
  const HitsToMCParticles& hitsToParticles=evHelper.GetHitToMCParticleAssociations();

  FillClusterPCA(outputPtr,evHelper);

  static void LArPIDCalculator::FillClusterPCA(LArPID* outputPtr,const EventHelper& evHelper){
    const ClusterVector& clusters=evHelper.GetClusters();

    for(auto cl=clusters.cbegin();cl!=clusters.cend();++cl){
      
  */
  }
}
