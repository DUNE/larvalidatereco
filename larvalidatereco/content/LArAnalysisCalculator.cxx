#include "TDatabasePDG.h"

#include "Geometry/TPCGeo.h"
#include "Geometry/Geometry.h"
#include "larvalidatereco/framework/VarHelper.h"
#include "larvalidatereco/interface/EventHelper.h"
#include "LArAnalysisCalculator.h"

#include "TGraph.h"
#include "TFile.h"

namespace lar_valrec{

  class LArAnalysis_register{
  public:
    LArAnalysis_register(){
        VarHelper::RegisterOutput<LArAnalysis,LArAnalysisCalculator>("Analysis");
    }
  };
  static LArAnalysis_register __LArAnalysisRegister;

  bool LArAnalysisCalculator::IsInActiveRegion(const TVector3& position){
    
    art::ServiceHandle<geo::Geometry> geomService;
    geo::Geometry::TPC_iterator tpcIter;

    for(;tpcIter++;){
      const geo::TPCGeo* tpcGeom=tpcIter.get();
      if(!tpcGeom) break;
      TVector3 localPos=tpcGeom->WorldToLocal(position);
      const TGeoVolume* vol=tpcGeom->ActiveVolume();
      Double_t xyz[3]={0,0,0};
      localPos.GetXYZ(xyz);
      if(vol->Contains(xyz)) return true;
    }
    return false;
  }
      

  double LArAnalysisCalculator::GetTrackLength(std::vector<TVector3>& points,int trackIndex){

  TVector3 entryPoint(0.,0.,0.);
  TVector3 exitPoint(0.,0.,0.);
  auto firstValidPoint=points.end();
  auto firstInvalidPoint=points.end();
  bool hasEntered=false;
  bool isInDetector=false;

  const TVector3* prevPoint=&*points.begin();

  //Identify entry and exit points of track wrt active detector
  for(auto pointIter=points.begin();pointIter!=points.end();++pointIter){

    //Current point is inside the detector
    if(IsInActiveRegion(*pointIter)){
      isInDetector=true;

      //As far as we know, this track will now be in the detector till its end
      firstInvalidPoint=points.end();      

      //If this is the first point inside the detector, make a note of it, and
      //find the entry point by tracking back towards the previous point
      if(!hasEntered){
	hasEntered=true;
	firstValidPoint=pointIter;

	entryPoint=*prevPoint;
	for(;
	    !IsInActiveRegion(entryPoint)&&(*pointIter-*prevPoint).Angle(*pointIter-entryPoint)<0.1;
	    //entryPoint+=(*pointIter-entryPoint).Unit()*0.1){}
	    entryPoint+=(*pointIter-*prevPoint).Unit()*0.1){}
      }
    }

    //Current point is outside the detector
    else{
      //If the previous point was in the detector, then make a note of it, and
      //find the exit point by tracking back towards the last point
      //NB, this point may be updated later if the particle reenters the detector
      if(isInDetector){
	isInDetector=false;
	firstInvalidPoint=pointIter;      
	exitPoint=*prevPoint;
	for(;
	    IsInActiveRegion(exitPoint);
	    exitPoint+=(*pointIter-*prevPoint).Unit()*0.1){}
      }
    }
    //We need the previous point to allow us to track back to entry/exit points
    prevPoint=&*pointIter;
  }

  double trackLength=0.;
  prevPoint=&entryPoint;
  for(auto pointIter=firstValidPoint;pointIter!=firstInvalidPoint;++pointIter){
    trackLength+=(*pointIter-*prevPoint).Mag();
    prevPoint=&*pointIter;
  }
  if(!isInDetector)trackLength+=(exitPoint-*prevPoint).Mag();
  return trackLength;
  }
  /*auto isInDetector = [] (const TVector3& p) { 
  static const double y_min=-80;
  static const double y_max=110;
  static const double x_min=-30;
  static const double x_max=220;
  static const double z_min=0;
  static const double z_max=155;
  return p.X()>x_min&&p.X()<x_max&&
  p.Y()>y_min&&p.Y()<y_max&&
  p.Z()>z_min&&p.Z()<z_max;};*/
  /*
  TVector3 lastPoint=*(points.begin());
  for(auto pointIter=points.begin();pointIter!=points.end();++pointIter){
    TVector3 point=*pointIter;
    if(IsInActiveRegion(point)){
      if(!hasEntered){
	hasEntered=true;
	for(;
	    !IsInActiveRegion(lastPoint);
	    lastPoint+=(point-lastPoint).Unit()*0.1){}
	gr->SetPoint(iPoint++,lastPoint.Y(),lastPoint.Z());
      }
      trackLength+=(lastPoint-point).Mag();
      gr->SetPoint(iPoint++,point.Y(),point.Z());
    }
    else{
      if(hasEntered){
	for(;
	    IsInActiveRegion(lastPoint);
	    lastPoint+=(point-lastPoint).Unit()*0.1){
	  trackLength+=0.1;
	}
	gr->SetPoint(iPoint++,lastPoint.Y(),lastPoint.Z());
	break;
      }
    }
    lastPoint=point;
  }
  */


void LArAnalysisCalculator::Calculate(TObject* tObjectPtr,const VarHelper& varHelper,const EventHelper& evHelper){
  LArAnalysis* outputPtr=dynamic_cast<LArAnalysis*>(tObjectPtr);
  if(!outputPtr){
    exit(1);
  }
  outputPtr->Clear();

  
  this->FillEventMetadata(outputPtr,evHelper);
  this->FillEventMCTraj(outputPtr,evHelper);
  this->FillEventMCVertices(outputPtr,evHelper);
  this->FillEventRecoTracks(outputPtr,evHelper);
  this->FillEventRecoClusters(outputPtr,evHelper);
  this->FillEventHits(outputPtr,evHelper);
}

void LArAnalysisCalculator::FillEventMetadata(LArAnalysis* outputPtr,const EventHelper& evHelper){
  
  outputPtr->EventEField=evHelper.GetEField();
  outputPtr->EventEventRun=evHelper.GetRun();
  outputPtr->EventID=evHelper.GetEventID();
  outputPtr->EventSubRun=evHelper.GetSubRun();
  outputPtr->EventT0=evHelper.GetEventT0();
}

void LArAnalysisCalculator::FillEventMCTraj(LArAnalysis* outputPtr,const EventHelper& evHelper){

  //Can't get CMake to find root's libEG so commenting this out for now.
  //static TDatabasePDG pdgDB;
  outputPtr->NTrajMC=0;
  const MCParticleVector& mcParticles=evHelper.GetMCParticles();
  const MCParticlesToHits& mcParticlesToHits=evHelper.GetMCParticleToHitAssociations();
  for(auto i=mcParticlesToHits.begin();i!=mcParticlesToHits.end();++i){
  }
  for(auto particle=mcParticles.begin();particle!=mcParticles.end();++particle){
    ++(outputPtr->NTrajMC);

    //Static particle properties
    int pdgCode=(*particle)->PdgCode();
    outputPtr->TrajPDGMC.push_back(pdgCode);
    if(mcParticlesToHits.count(*particle)){
	outputPtr->TrajNHitsMC.push_back(mcParticlesToHits.at(*particle).size());
    }
    else{
      outputPtr->TrajNHitsMC.push_back(0);
    }
    //These need to be implemented using TDatabasePDG or similar
    //    outputPtr->TrajChargeMC.push_back(pdgDB.GetParticle(pdgCode)->Charge());
    //outputPtr->TrajMassMC.push_back(pdgDB.GetParticle(pdgCode)->Mass());
    //outputPtr->TrajNameMC.push_back(pdgDB.GetParticle(pdgCode)->GetName());

    //Particle and parent indices
    outputPtr->TrajIDMC.push_back((*particle)->TrackId());
    outputPtr->TrajParentIDMC.push_back((*particle)->Mother());

    std::vector<TVector3> trackPoints;

    for(unsigned int point=0;point!=(*particle)->NumberTrajectoryPoints();++point){
      trackPoints.push_back((*particle)->Position(point).Vect());
    }

    //Start and end momenta and positions
    outputPtr->TrajStart4MomMC.push_back((*particle)->Momentum());
    outputPtr->TrajStart4PosMC.push_back((*particle)->Position());
    outputPtr->TrajEnd4MomMC.push_back((*particle)->EndMomentum());
    outputPtr->TrajEnd4PosMC.push_back((*particle)->EndPosition());
    outputPtr->TrajLengthMC.push_back(GetTrackLength(trackPoints,outputPtr->NTrajMC==1?outputPtr->EventID:-1));
  }
}

void LArAnalysisCalculator::FillEventMCVertices(LArAnalysis* outputPtr,const EventHelper& evHelper){

  const MCTruthVector& vertices=evHelper.GetMCVertices();

  for(auto vIter=vertices.begin();vIter!=vertices.end();++vIter){
    ++outputPtr->NVtxMC;
    if((*vIter)->NeutrinoSet()){
	outputPtr->Vtx4PosMC.push_back((*vIter)->GetNeutrino().Nu().EndPosition());
	outputPtr->VtxNuMomMC.push_back((*vIter)->GetNeutrino().Nu().Momentum());
	outputPtr->VtxNuPDGMC.push_back((*vIter)->GetNeutrino().Nu().PdgCode());
	outputPtr->VtxTargetPDGMC.push_back((*vIter)->GetNeutrino().Target());
	outputPtr->VtxReactionCodeMC.push_back((*vIter)->GetNeutrino().Mode());
    }
    else{
      outputPtr->Vtx4PosMC.push_back(TLorentzVector(-1,-1,-1,-1));
      outputPtr->VtxNuMomMC.push_back(TLorentzVector(-1,-1,-1,-1));
      outputPtr->VtxNuPDGMC.push_back(-1);
      outputPtr->VtxTargetPDGMC.push_back(-1);
      outputPtr->VtxReactionCodeMC.push_back(-1);
    }
    outputPtr->VtxNPrimTrajMC.push_back((*vIter)->NParticles());
    outputPtr->VtxPrimTrajIDsMC.emplace_back();

    //For some reason the track ID always seems to return -1 (for particle gun anyway).
    //PDG code is correct so not sure why this is...
    for(int i=0;i<(*vIter)->NParticles();++i){
      outputPtr->VtxPrimTrajIDsMC.back().push_back((*vIter)->GetParticle(i).TrackId());
    }
  }
}

void LArAnalysisCalculator::FillEventRecoTracks(LArAnalysis* outputPtr,const EventHelper& evHelper){

  const TrackVector& tracks=evHelper.GetTracks();
  const TracksToHits& tracksToHits=evHelper.GetTracksToHits();

  outputPtr->NTracks=0;
  for(auto track=tracks.begin();track!=tracks.end();++track){
    ++(outputPtr->NTracks);
      outputPtr->TrackID.push_back((*track)->ID());
      outputPtr->TrackStart4Pos.push_back(TLorentzVector((*track)->Vertex(),std::numeric_limits<double>::max()));
      outputPtr->TrackEnd4Pos.push_back(TLorentzVector((*track)->End(),std::numeric_limits<double>::max()));
      outputPtr->TrackNHits.push_back(tracksToHits.at(*track).size());
      //Bits we still need to implement
      outputPtr->TrackPitch.push_back(std::numeric_limits<float>::max());
      outputPtr->TrackCharge.push_back(std::numeric_limits<int>::max());

      std::vector<TVector3> trackPoints;

      for(unsigned int point=0;point!=(*track)->NumberTrajectoryPoints();++point){
	trackPoints.push_back((*track)->LocationAtPoint(point));
      }

      outputPtr->TrackLength.push_back(GetTrackLength(trackPoints,outputPtr->NTracks==1?1000000+evHelper.GetEventID():-1));
  }
}

void LArAnalysisCalculator::FillEventRecoClusters(LArAnalysis* outputPtr,const EventHelper& evHelper){

  const ClusterVector& clusters=evHelper.GetClusters();
  outputPtr->NClusters=0;
  for(auto cluster=clusters.begin();cluster!=clusters.end();++cluster){
    ++(outputPtr->NClusters);
    
    outputPtr->ClusterCharge.push_back((*cluster)->Integral());
    outputPtr->ClusterID.push_back((*cluster)->ID());
    outputPtr->ClusterNHits.push_back((*cluster)->NHits());
    outputPtr->ClusterTrackID.push_back(std::numeric_limits<int>::max());    
  }
}

void LArAnalysisCalculator::FillEventHits(LArAnalysis* outputPtr,const EventHelper& evHelper){
  
    const HitVector& hits=evHelper.GetHits();
    const HitsToClusters& hitsToClusters=evHelper.GetHitsToClusters();
    const HitsToTracks& hitsToTracks=evHelper.GetHitsToTracks();
    const HitsToSpacePoints& hitsToSpacePoints=evHelper.GetHitsToSpacePoints();
    outputPtr->NHits=0;
    for(auto hit=hits.begin();hit!=hits.end();++hit){
      ++(outputPtr->NHits);

      outputPtr->HitChannel.push_back((*hit)->Channel());  
      outputPtr->HitCharge.push_back((*hit)->Integral());
      if(hitsToClusters.count(*hit)){
	outputPtr->HitClusterID.push_back(hitsToClusters.at(*hit)->ID());
      }
      else{
	outputPtr->HitClusterID.push_back(-1);
      }
      
      if(hitsToTracks.count(*hit)){
	outputPtr->HitTrackID.push_back(hitsToTracks.at(*hit)->ID());  
      }
      else{
	outputPtr->HitTrackID.push_back(-1);
      }

      outputPtr->HitPeakT.push_back((*hit)->PeakTime());
      outputPtr->HitPlane.push_back((*hit)->WireID().Plane);
      outputPtr->HitWire.push_back((*hit)->WireID().Wire);      
      outputPtr->HitTPC.push_back((*hit)->WireID().TPC);
      if(hitsToSpacePoints.count(*hit)&&hitsToSpacePoints.at(*hit).isNonnull()){
	outputPtr->Hit3Pos.push_back(TVector3(hitsToSpacePoints.at(*hit)->XYZ()));
	outputPtr->HitIsMatched.push_back(true);
      }
      else{
	outputPtr->Hit3Pos.push_back(TVector3(0,0,0));
	outputPtr->HitIsMatched.push_back(false);
      }
    }
  }
}

