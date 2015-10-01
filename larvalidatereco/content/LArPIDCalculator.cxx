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


  void LArPIDCalculator::FillEventPID(TObject* tObjectPtr,const VarHelper& varHelper,const EventHelper& evHelper)
  {
    LArPID* outputPtr=dynamic_cast<LArPID*>(tObjectPtr);
    if(!outputPtr){
      exit(1);
    }

    const double MoliereRadius = 10.1;
    const double MoliereRadiusFraction = 0.5;
    const double trackFraction = 0.2;
    const double pca_sp_check_threshold = 0.01;
    const int min_pca_spacepoints = 25;
    const double max_dEdx = 50.0;
    const double min_trackpitch = 0.1;

    // Define the TPrincipal
    TPrincipal* principal = new TPrincipal(3,"D"); 
    // Define variables to hold the eigenvalues and eigenvectors
     const TVectorD* eigenval = new TVectorD();
     const TMatrixD* eigenvec = new TMatrixD();
     const TMatrixD* covar = new TMatrixD();
     const TVectorD* meanval = new TVectorD();
     const Double_t *hits_spacepoints = new Double_t[3];
     Double_t* pca_spacepoints = new Double_t[3];
     TVector3* pca_spacepoints_output = new TVector3();
     double *pca_spacepoints_check = new double[3];
     std::vector<TVector3>* pca_hsp_cont = new std::vector<TVector3>;
     // LOOP OVER TRACKS
     //const TrackVector& TrackVector=evHelper.GetTracks();
     const TracksToHits& TracksToHits=evHelper.GetTracksToHits();
     const HitsToSpacePoints& hitsToSpacePoints=evHelper.GetHitsToSpacePoints();
     const MCParticleVector& mcParticles=evHelper.GetMCParticles();
     const TracksToCalo& tracksToCalo=evHelper.GetTracksToCalo();

     calo::CalorimetryAlg CaloAlg(evHelper.GetParameterSet());

     std::vector<double> pca_spacepoints_0;
     std::vector<double>::iterator pcaIt;

     //Loop over hits from first reconstructed track
     if(TracksToHits.size() >= 1)
       {
         auto tracks = TracksToHits.begin();

	   // FILL SPACEPOINTS ARRAY WITH SPACEPOINTS BY LOOPING OVER ALL MATCHED HITS
	   for ( auto hit = (tracks->second).begin(); hit!=(tracks->second).end(); hit++ )
	      {
	      // CHECK FOR EXISTENCE OF KEY TO PREVENT MAP OVERRUN
	      if(hitsToSpacePoints.count(*hit)) 
	        {
		  TVector3 xyz=hitsToSpacePoints.at(*hit)->XYZ();
		  principal->AddRow( hitsToSpacePoints.at(*hit)->XYZ() );
		}
	      }
	// PERFORM PCA
	principal->MakePrincipals();
	// GET EIGENVALUES AND EIGENVECTORS
	eigenval=principal->GetEigenValues();
	eigenvec=principal->GetEigenVectors();
	covar = principal->GetCovarianceMatrix();
        meanval=principal->GetMeanValues();

	const double *eval = new double[3];
        eval = eigenval->GetMatrixArray();
	const double *evec = new double[9];
        evec = eigenvec->GetMatrixArray();
	const double *data_means = new double[3];
        data_means = meanval->GetMatrixArray();

	//Calculate angle between principal component and true direction of particle (validation check)
	auto particle=mcParticles.begin();
        double momNorm = sqrt((*particle)->Momentum().X() * (*particle)->Momentum().X()
			      + (*particle)->Momentum().Y() * (*particle)->Momentum().Y()
			      + (*particle)->Momentum().Z() * (*particle)->Momentum().Z());
	//Angle is arccos of scalar product between principal component and (normalised) true direction of particle
	double angleWithTrueParticle = acos((evec[0] * (*particle)->Momentum().X()
                                             + evec[3] * (*particle)->Momentum().Y()
                                             + evec[6] * (*particle)->Momentum().Z()) / momNorm);

	std::cout<<"Angle between principal component and true direction of particle = "<<angleWithTrueParticle<<" rad"<<std::endl;

	// WRITE TO NTUPLE
	outputPtr->EigenValues.push_back(*eigenval);
	outputPtr->EigenVectors.push_back(*eigenvec);
	outputPtr->Covariance.push_back(*covar);
	outputPtr->MeanValues.push_back(*meanval);
        outputPtr->AnglePrincipalTrueTrack.push_back(angleWithTrueParticle);
	outputPtr->EvalRatio.push_back(sqrt(eval[1] * eval[1] + eval[2] * eval[2]) / eval[0]);

	double chargeCore = 0.0;
	double chargeHalo = 0.0;

	double chargeCon = 0.0;

	pca_spacepoints_0.clear();

	// *** ITERATE OVER THE HITS AGAIN TO CALCULATE THE NEW SPACEPOINTS UNDER PCA TRANSFORMATION ***
                // FILL SPACEPOINTS ARRAY WITH SPACEPOINTS BY LOOPING OVER ALL MATCHED HITS     
                for ( auto hit = (tracks->second).begin(); hit!=(tracks->second).end(); hit++ )
                  {
                    // CHECK FOR EXISTENCE OF KEY TO PREVENT MAP OVERRUN                        
                    if(hitsToSpacePoints.count(*hit)) {
		      // CALCULATE PCA ON SPACEPOINTS AND WRITE TO NTUPLE
		      hits_spacepoints = hitsToSpacePoints.at(*hit)->XYZ();
		      principal->X2P(hits_spacepoints,pca_spacepoints);
		      pca_spacepoints_output->SetXYZ(pca_spacepoints[0], pca_spacepoints[1], pca_spacepoints[2] );
		      pca_hsp_cont->push_back(*pca_spacepoints_output);

		      pca_spacepoints_0.push_back(pca_spacepoints[0]);

		      //Subtract means from hit spacepoints and multiply by transpose of matrix of eigenvectors 
		      //This should agree with the PCA spacepoints calculated by TPrincipal (validation check)
		      pca_spacepoints_check[0] = evec[0] * (hits_spacepoints[0] - data_means[0]) 
			                       + evec[3] * (hits_spacepoints[1] - data_means[1])
                                               + evec[6] * (hits_spacepoints[2] - data_means[2]);
		      pca_spacepoints_check[1] = evec[1] * (hits_spacepoints[0] - data_means[0])
			                       + evec[4] * (hits_spacepoints[1] - data_means[1])
                                               + evec[7] * (hits_spacepoints[2] - data_means[2]);
		      pca_spacepoints_check[2] = evec[2] * (hits_spacepoints[0] - data_means[0])
			                       + evec[5] * (hits_spacepoints[1] - data_means[1])
                                               + evec[8] * (hits_spacepoints[2] - data_means[2]);
                      
		      if(fabs(pca_spacepoints_check[0] - pca_spacepoints[0]) > pca_sp_check_threshold
			 || fabs(pca_spacepoints_check[1] - pca_spacepoints[1]) > pca_sp_check_threshold	
			 || fabs(pca_spacepoints_check[2] - pca_spacepoints[2]) > pca_sp_check_threshold)
			{
			  std::cerr<<"LArPIDCalculator::LArPIDPCA() PCA spacepoint validation check failed ...... exiting."<<std::endl;
			  exit(1);
			}

		      if(sqrt(pca_spacepoints[1] * pca_spacepoints[1] + pca_spacepoints[2] * pca_spacepoints[2]) < MoliereRadiusFraction * MoliereRadius)
			chargeCore += (*hit)->Integral();
		      else
			chargeHalo += (*hit)->Integral();

                      chargeCon += (*hit)->Integral() / sqrt(pca_spacepoints[1] * pca_spacepoints[1] + pca_spacepoints[2] * pca_spacepoints[2]);
                    }
		  }

	outputPtr->PCAHitsSpacePoints.push_back(*pca_hsp_cont);
	outputPtr->ChargeRatioCoreHalo.push_back(chargeHalo / chargeCore);
	outputPtr->Concentration.push_back(chargeCon);

	double trackPitchC = 1.0;

	if(tracksToCalo.size() >= 1)
          {
  	  auto caloIt=tracksToCalo.begin();
	  art::Ptr<anab::Calorimetry> calo=caloIt->second;
	  trackPitchC = calo->TrkPitchC();
	  }

	double track_pca_start = 0.0;
        double track_pca_end = 0.0;

	if(pca_spacepoints_0.size() >= min_pca_spacepoints)
	  { 
	  sort(pca_spacepoints_0.begin(), pca_spacepoints_0.end());

	  pcaIt = pca_spacepoints_0.begin();
          track_pca_start = (*pcaIt);
          pcaIt = pca_spacepoints_0.end()-1;
	  track_pca_end = (*pcaIt);
	  }

	double dEdxAmpStart = 0.0;
	double dEdxAmpEnd = 0.0;
	double dEdxAreaStart = 0.0;
	double dEdxAreaEnd = 0.0;

	double stdDevStart = 0.0;
        double stdDevEnd = 0.0;

	int nhits_dEdx_amp_start = 0;
        int nhits_dEdx_amp_end = 0;
        int nhits_dEdx_area_start = 0;
        int nhits_dEdx_area_end = 0;
	int nhits_con_start = 0;
        int nhits_con_end = 0;

	//Loop over hits again to calculate average dE/dx and conicalness 
	for ( auto hit = (tracks->second).begin(); hit!=(tracks->second).end(); hit++ )
	  {
	    if(hitsToSpacePoints.count(*hit)) {

	      hits_spacepoints = hitsToSpacePoints.at(*hit)->XYZ();
	      principal->X2P(hits_spacepoints,pca_spacepoints);

	      if(pca_spacepoints_0.size() >= min_pca_spacepoints && pca_spacepoints[0] - track_pca_start < trackFraction * (track_pca_end - track_pca_start))
		{
		if(tracksToCalo.size() >= 1 && trackPitchC > min_trackpitch && CaloAlg.dEdx_AMP(*hit, trackPitchC, evHelper.GetEventT0()) < max_dEdx)
		  {
		  nhits_dEdx_amp_start++;
		  dEdxAmpStart += CaloAlg.dEdx_AMP(*hit, trackPitchC, evHelper.GetEventT0());
		  }
		if(tracksToCalo.size() >= 1 && trackPitchC > min_trackpitch && CaloAlg.dEdx_AREA(*hit, trackPitchC, evHelper.GetEventT0()) < max_dEdx)
		  {
		  nhits_dEdx_area_start++;
		  dEdxAreaStart += CaloAlg.dEdx_AREA(*hit, trackPitchC, evHelper.GetEventT0());
		  }
		nhits_con_start++;
  		stdDevStart += pca_spacepoints[1] * pca_spacepoints[1] + pca_spacepoints[2] * pca_spacepoints[2];
         	}
	      if(pca_spacepoints_0.size() >= min_pca_spacepoints && track_pca_end - pca_spacepoints[0] < trackFraction * (track_pca_end - track_pca_start))
		{
		if(tracksToCalo.size() >= 1 && trackPitchC > min_trackpitch && CaloAlg.dEdx_AMP(*hit, trackPitchC, evHelper.GetEventT0()) < max_dEdx)
		  {
		  nhits_dEdx_amp_end++;
		  dEdxAmpEnd += CaloAlg.dEdx_AMP(*hit, trackPitchC, evHelper.GetEventT0());
		  }
                if(tracksToCalo.size() >= 1 && trackPitchC > min_trackpitch && CaloAlg.dEdx_AREA(*hit, trackPitchC, evHelper.GetEventT0()) < max_dEdx)
		  {
		  nhits_dEdx_area_end++;
                  dEdxAreaEnd += CaloAlg.dEdx_AREA(*hit, trackPitchC, evHelper.GetEventT0());
		  }
		nhits_con_end++;
		stdDevEnd += pca_spacepoints[1] * pca_spacepoints[1] + pca_spacepoints[2] * pca_spacepoints[2];
		}
	    }
	  }

	if(tracksToCalo.size() >= 1 && pca_spacepoints_0.size() >= min_pca_spacepoints && trackPitchC > min_trackpitch && nhits_dEdx_amp_start >= 1)
	  outputPtr->AvgedEdxAmpStart.push_back(dEdxAmpStart / nhits_dEdx_amp_start);
	else
	  outputPtr->AvgedEdxAmpStart.push_back(-999.9);

	if(tracksToCalo.size() >= 1 && pca_spacepoints_0.size() >= min_pca_spacepoints && trackPitchC > min_trackpitch && nhits_dEdx_area_start >= 1)
          outputPtr->AvgedEdxAreaStart.push_back(dEdxAreaStart / nhits_dEdx_area_start);
        else
          outputPtr->AvgedEdxAreaStart.push_back(-999.9);

	if(tracksToCalo.size() >= 1 && pca_spacepoints_0.size() >= min_pca_spacepoints && trackPitchC > min_trackpitch && nhits_dEdx_amp_end >= 1)
          outputPtr->AvgedEdxAmpEnd.push_back(dEdxAmpEnd / nhits_dEdx_amp_end);
        else
          outputPtr->AvgedEdxAmpEnd.push_back(-999.9);

        if(tracksToCalo.size() >= 1 && pca_spacepoints_0.size() >= min_pca_spacepoints && trackPitchC > min_trackpitch && nhits_dEdx_area_end >= 1)
          outputPtr->AvgedEdxAreaEnd.push_back(dEdxAreaEnd / nhits_dEdx_area_end);
        else
          outputPtr->AvgedEdxAreaEnd.push_back(-999.9);

	if(pca_spacepoints_0.size() >= min_pca_spacepoints && nhits_con_start >= 2)
	  stdDevStart = sqrt(stdDevStart / (nhits_con_start - 1));

	if(pca_spacepoints_0.size() >= min_pca_spacepoints && nhits_con_end >= 2)
	  stdDevEnd = sqrt(stdDevEnd / (nhits_con_end - 1));

	if(pca_spacepoints_0.size() < min_pca_spacepoints || nhits_con_start <= 1 || nhits_con_end <= 1)
          outputPtr->Conicalness.push_back(-999.9);
	else
	  outputPtr->Conicalness.push_back(stdDevStart / stdDevEnd);

	// Clear the data from the TPrincipal
	principal->Clear();
		  
       }//end of if(TracksToHits.size() >= 1) 

  }

  void LArPIDCalculator::Calculate(TObject* tObjectPtr,const VarHelper& varHelper,const EventHelper& evHelper){
  LArPID* outputPtr=dynamic_cast<LArPID*>(tObjectPtr);
  if(!outputPtr){
    exit(1);
  }
  outputPtr->Clear();
  
  this->FillEventPID(tObjectPtr, varHelper, evHelper);
  this->FillEventMetadata(outputPtr,evHelper);
  this->FillEventMCTraj(outputPtr,evHelper);
  this->FillEventTracksAndHits(outputPtr,evHelper);

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

  void LArPIDCalculator::FillEventTracksAndHits(LArPID* outputPtr,const EventHelper& evHelper){

    const TrackVector& tracks=evHelper.GetTracks();
    const TracksToHits& tracksToHits=evHelper.GetTracksToHits();
    const HitsToSpacePoints& hitsToSpacePoints=evHelper.GetHitsToSpacePoints();
    const HitsToTracks& hitsToTracks=evHelper.GetHitsToTracks();
    const TracksToCalo& tracksToCalo=evHelper.GetTracksToCalo();

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

        }//end of for(auto hit=tracksToHits.at(*track).begin();hit!=tracksToHits.at(*track).end();++hit){

      }//end of for(auto track=tracks.begin();track!=tracks.end();++track){


    //std::cout<<"Tracks / tracksToCalo size    "<<tracks.size()<<"    "<<tracksToCalo.size()<<std::endl;
 
    //for(auto iter=tracksToCalo.begin();iter!=tracksToCalo.end();++iter){
    // art::Ptr<anab::Calorimetry> calo=iter->second;

      /*
      std::vector<double> vdEdx = calo->dEdx();
      std::cout<<std::endl;
      std::cout<<"LArPIDCalculator dE/dx"<<std::endl;
      std::cout<<std::endl;
      std::cout<<"vdedx size =    "<<vdEdx.size()<<"    "<<(double)vdEdx.size() / nhits<<std::endl;

      for(auto dEdxIt= vdEdx.begin(); dEdxIt!=vdEdx.end(); ++dEdxIt)
	std::cout<<"    "<<calo->PlaneID()<<"    "<<(*dEdxIt);

      std::cout<<std::endl;
      */
      /*
      std::vector<double> vResRange = calo->ResidualRange();
      std::cout<<std::endl;
      std::cout<<"vResRange size =    "<<vResRange.size()<<std::endl;
      std::cout<<std::endl;

      for(auto resRangeIt= vResRange.begin(); resRangeIt!=vResRange.end(); ++resRangeIt)
	std::cout<<"    "<<calo->PlaneID()<<"    "<<(*resRangeIt);

      std::cout<<std::endl;
 

      std::vector<TVector3> vXYZ = calo->XYZ();
      std::cout<<std::endl;
      std::cout<<"vXYZ size =    "<<vXYZ.size()<<std::endl;
      std::cout<<std::endl;

      for(auto XYZIt= vXYZ.begin(); XYZIt!=vXYZ.end(); ++XYZIt)
	std::cout<<"    "<<calo->PlaneID()<<"    "<<(*XYZIt)[0]<<"    "<<(*XYZIt)[1]<<"    "<<(*XYZIt)[2]<<std::endl;
      */
      /*
      std::cout<<std::endl;
      std::cout<<"Track pitch C = "<<calo->TrkPitchC()<<std::endl;
      std::cout<<std::endl;
      
      std::vector<double> vTrkPitch = calo->TrkPitchVec();
      std::cout<<std::endl;
      std::cout<<"vTrkPitch size =    "<<vTrkPitch.size()<<std::endl;
      std::cout<<std::endl;
      */
      /*
      for(auto trkPitchIt= vTrkPitch.begin(); trkPitchIt!=vTrkPitch.end(); ++trkPitchIt)
	std::cout<<"    "<<calo->PlaneID()<<"    "<<(*trkPitchIt);
      */
      //std::cout<<std::endl;
   
    //}
    

  }//end of void LArPIDCalculator::FillEventdEdx()

}//end of namespace lar_valrec
