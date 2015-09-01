#include "TDatabasePDG.h"
#include "TPrincipal.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/Geometry.h"
#include "larvalidatereco/framework/VarHelper.h"
#include "larvalidatereco/interface/EventHelper.h"
#include "LArPIDCalculator.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TFile.h"
#include "TMatrix.h"
#include "TVectorD.h"

namespace lar_valrec{

  class LArPID_register{
  public:
    LArPID_register(){
        VarHelper::RegisterOutput<LArPID,LArPIDCalculator>("PID");
    }
  };
  static LArPID_register __LArPIDRegister;


  void LArPIDPCA(TObject* tObjectPtr,const VarHelper& varHelper,const EventHelper& evHelper)
  {
    LArPID* outputPtr=dynamic_cast<LArPID*>(tObjectPtr);
    if(!outputPtr){
      exit(1);
    }
    outputPtr->Clear();

    // Define the TPrincipal
    TPrincipal* principal = new TPrincipal(3,"ND"); 
    // Define variables to hold the eigenvalues and eigenvectors
     const TVectorD* eigenval = new TVectorD();
     const TMatrixD* eigenvec = new TMatrixD();
     const TMatrixD* covar = new TMatrixD();
     const Double_t *hits_spacepoints = new Double_t[3];
     Double_t* pca_hits_spacepoints = new Double_t[3];
     TVector3* pca_hits_spacepoints_output = new TVector3();
     std::vector<TVector3>* pca_hsp_cont = new std::vector<TVector3>;
     // LOOP OVER TRACKS
     //const TrackVector& TrackVector=evHelper.GetTracks();
     const TracksToHits& TracksToHits=evHelper.GetTracksToHits();
     const HitsToSpacePoints& hitsToSpacePoints=evHelper.GetHitsToSpacePoints();

     // LOOP OVER TRACK TO HITS
     for (auto tracks = (TracksToHits.begin()); tracks!=(TracksToHits.end()); tracks++)
       {
		// FILL SPACEPOINTS ARRAY WITH SPACEPOINTS BY LOOPING OVER ALL MATCHED HITS
		for ( auto hits = (tracks->second).begin(); hits!=(tracks->second).end(); hits++ )
		  {
		    // CHECK FOR EXISTENCE OF KEY TO PREVENT MAP OVERRUN
		    if(hitsToSpacePoints.count(*hits)) 
		      {
			principal->AddRow( hitsToSpacePoints.at(*hits)->XYZ() );
		      }
		  }
	// PERFORM PCA
	principal->MakePrincipals();
	// GET EIGENVALUES AND EIGENVECTORS
	eigenval=principal->GetEigenValues();
	eigenvec=principal->GetEigenVectors();
	covar = principal->GetCovarianceMatrix();
	// WRITE EIGENVECTORS AND EIGENVALUES TO NTUPLE
	outputPtr->EigenValues.push_back(*eigenval);
	outputPtr->EigenVectors.push_back(*eigenvec);
	outputPtr->Covariance.push_back(*covar);
	// *** ITERATE OVER THE HITS AGAIN TO CALCULATE THE NEW SPACEPOINTS UNDER PCA TRANSFORMATION ***
                // FILL SPACEPOINTS ARRAY WITH SPACEPOINTS BY LOOPING OVER ALL MATCHED HITS     
                for ( auto hits = (tracks->second).begin(); hits!=(tracks->second).end(); hits++ )
                  {
                    // CHECK FOR EXISTENCE OF KEY TO PREVENT MAP OVERRUN                        
                    if(hitsToSpacePoints.count(*hits)) {
		      // CALCULATE PCA ON SPACEPOINTS AND WRITE TO NTUPLE
		      hits_spacepoints = hitsToSpacePoints.at(*hits)->XYZ();
		      principal->X2P(hits_spacepoints,pca_hits_spacepoints);
		      pca_hits_spacepoints_output->SetXYZ(pca_hits_spacepoints[0], pca_hits_spacepoints[1], pca_hits_spacepoints[2] );
		      pca_hsp_cont->push_back(*pca_hits_spacepoints_output);	
                    }
		  }
	outputPtr->PCAHitsSpacePoints.push_back(*pca_hsp_cont);
	// Clear the data from the TPrincipal
	principal->Clear();
		  
       }
  }

  void LArPIDCalculator::Calculate(TObject* tObjectPtr,const VarHelper& varHelper,const EventHelper& evHelper){
  LArPID* outputPtr=dynamic_cast<LArPID*>(tObjectPtr);
  if(!outputPtr){
    exit(1);
  }
  outputPtr->Clear();
  
  // Call PCA Method
  LArPIDPCA(tObjectPtr, varHelper, evHelper);

  this->FillEventMetadata(outputPtr,evHelper);
  this->FillEventMCTraj(outputPtr,evHelper);
  this->FillEventdEdx(outputPtr,evHelper);

  }

  bool LArPIDCalculator::IsInActiveRegion(const TVector3& position){
    
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

  double LArPIDCalculator::GetTrackLength(std::vector<TVector3>& points,int trackIndex){

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

  void LArPIDCalculator::FillEventMetadata(LArPID* outputPtr,const EventHelper& evHelper){
  
    outputPtr->EventEField=evHelper.GetEField();
    outputPtr->EventEventRun=evHelper.GetRun();
    outputPtr->EventID=evHelper.GetEventID();
    outputPtr->EventSubRun=evHelper.GetSubRun();
    outputPtr->EventT0=evHelper.GetEventT0();
  }

  void LArPIDCalculator::FillEventMCTraj(LArPID* outputPtr,const EventHelper& evHelper){

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
  void LArPIDCalculator::FillEventdEdx(LArPID* outputPtr,const EventHelper& evHelper){

    //CalorimetryAlg has art::ServiceHandle<geo::Geometry> geom as a private attribute.
    //Please see http://nusoft.fnal.gov/larsoft/doxsvn/html/classcalo_1_1CalorimetryAlg.html.
    //Looks as though you are meant to use it but tried calling it from public method and this did not compile.
    //Therefore instantiate a Geometry object in order to get wire and plane pitches.
    art::ServiceHandle<geo::Geometry> geom;

    calo::CalorimetryAlg CaloAlg(evHelper.GetParameterSet());

    const TrackVector& tracks=evHelper.GetTracks();
    const TracksToHits& tracksToHits=evHelper.GetTracksToHits();
    const HitsToSpacePoints& hitsToSpacePoints=evHelper.GetHitsToSpacePoints();
    const HitsToTracks& hitsToTracks=evHelper.GetHitsToTracks();

    unsigned int planeLastHit = 0;

    outputPtr->NHits=0;
    outputPtr->NTracks=0;

    for(auto track=tracks.begin();track!=tracks.end();++track){
      ++(outputPtr->NTracks);
      outputPtr->TrackID.push_back((*track)->ID());
      outputPtr->TrackStart4Pos.push_back(TLorentzVector((*track)->Vertex(),std::numeric_limits<double>::max()));
      outputPtr->TrackEnd4Pos.push_back(TLorentzVector((*track)->End(),std::numeric_limits<double>::max()));
      outputPtr->TrackNHits.push_back(tracksToHits.at(*track).size());

      std::vector<TVector3> trackPoints;

      for(unsigned int point=0;point!=(*track)->NumberTrajectoryPoints();++point){
	trackPoints.push_back((*track)->LocationAtPoint(point));
      }

      outputPtr->TrackLength.push_back(GetTrackLength(trackPoints,outputPtr->NTracks==1?1000000+evHelper.GetEventID():-1));

      for(auto hit=tracksToHits.at(*track).begin();hit!=tracksToHits.at(*track).end();++hit)
        {
	  ++(outputPtr->NHits);

	  outputPtr->HitChannel.push_back((*hit)->Channel());
	  outputPtr->HitCharge.push_back((*hit)->Integral());

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
	  
          //geom->WirePitch() and geom->PlanePitch() are defined at http://nusoft.fnal.gov/larsoft/doxsvn/html.1.7.1/classgeo_1_1Geometry.html
          //Use WirePitch() if next hit is in same plane, or PlanePitch() if it is in a different plane.
          //First 2 arguments in WirePitch() are 0 and 1 to give pitch between 2 adjacent wires, third argument is wire plane, fourth argument is wire TPC.
	  //Wire pitches agree with those given in DocDB 7550.
          //First 2 arguments in PlanePitch() are the planes, third argument is the TPC. 
	  if((*hit)->WireID().Plane == planeLastHit)
	    {
            outputPtr->HitdEdxAmp.push_back(CaloAlg.dEdx_AMP(*hit, geom->WirePitch(0,1,(*hit)->WireID().Plane, (*hit)->WireID().TPC), evHelper.GetEventT0()));
	    outputPtr->HitdEdxArea.push_back(CaloAlg.dEdx_AREA(*hit, geom->WirePitch(0,1,(*hit)->WireID().Plane, (*hit)->WireID().TPC), evHelper.GetEventT0()));
	    }
          else if((*hit)->WireID().Plane < planeLastHit)
	    {
	      outputPtr->HitdEdxAmp.push_back(CaloAlg.dEdx_AMP(*hit, geom->PlanePitch((*hit)->WireID().Plane, planeLastHit, (*hit)->WireID().TPC), evHelper.GetEventT0()));
	      outputPtr->HitdEdxArea.push_back(CaloAlg.dEdx_AREA(*hit, geom->PlanePitch((*hit)->WireID().Plane, planeLastHit, (*hit)->WireID().TPC), evHelper.GetEventT0()));
	    }
          else if((*hit)->WireID().Plane > planeLastHit)
	    {
	      outputPtr->HitdEdxAmp.push_back(CaloAlg.dEdx_AMP(*hit, geom->PlanePitch(planeLastHit, (*hit)->WireID().Plane, (*hit)->WireID().TPC), evHelper.GetEventT0()));
	      outputPtr->HitdEdxArea.push_back(CaloAlg.dEdx_AREA(*hit, geom->PlanePitch(planeLastHit, (*hit)->WireID().Plane, (*hit)->WireID().TPC), evHelper.GetEventT0()));
            }
	  if(hitsToSpacePoints.count(*hit)&&hitsToSpacePoints.at(*hit).isNonnull())
            {
	    outputPtr->Hit3Pos.push_back(TVector3(hitsToSpacePoints.at(*hit)->XYZ()));
	    outputPtr->HitIsMatched.push_back(true);
	    }
	  else
            {
	    outputPtr->Hit3Pos.push_back(TVector3(0,0,0));
	    outputPtr->HitIsMatched.push_back(false);
	    }

	  planeLastHit = (*hit)->WireID().Plane;

        }//end of for(auto hit=tracksToHits.at(*track).begin();hit!=tracksToHits.at(*track).end();++hit){

      }//end of for(auto track=tracks.begin();track!=tracks.end();++track){

  }

}
