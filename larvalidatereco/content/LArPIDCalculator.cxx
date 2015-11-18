
#include "TDatabasePDG.h"
#include "TPrincipal.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/Geometry.h"
#include "larvalidatereco/framework/VarHelper.h"
#include "larvalidatereco/interface/EventHelper.h"
#include "LArPIDCalculator.h"
//#include "TVector3.h"
#include "TGraph.h"
//#include "TGraph2D.h"
#include "TFile.h"
#include "TF1.h"
#include "TF2.h"
#include "TMatrix.h"
#include "TVectorD.h"
#include "TVirtualFitter.h"

namespace lar_valrec{

  //Active volume of 35 ton detector taken from slide 8 in talk by Tristan Blackburn on 19 August 2015, slides are at
  //https://indico.fnal.gov/getFile.py/access?contribId=4&resId=0&materialId=slides&confId=10284

  const double activeVolMinX = -35.18;
  const double activeVolMaxX = 222.46;
  const double activeVolMinY = -84.22;
  const double activeVolMaxY = 115.09;
  const double activeVolMinZ = -2.04;
  const double activeVolMaxZ = 156.78;

  //Allow a few cm since last hit will not be exactly at edge of active volume 
  const double fiducialDist = 4.0;

  //Unit normal vectors to wires in each plane in the 35 ton detector
  //These are taken from the diagram in slide 27 of talk by Jeff Hartnell at UCL on 15 May 2015 
  //Slides are at http://indico.hep.manchester.ac.uk/getFile.py/access?contribId=5&resId=0&materialId=slides&confId=4755
  const TVector3 normToWiresU(0.0, 0.7071, -0.7071);
  const TVector3 normToWiresV(0.0, 0.7071, 0.7071);
  const TVector3 normToWiresW(0.0, 0.0, 1.0);

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
    const double MoliereRadiusFraction = 0.2;
    const double dEdxTrackFraction = 0.2;
    const double chargeTrackFraction = 0.2;
    const double coreHaloFraction = 0.25;
    const double pca_sp_check_threshold = 0.01;
    const int min_spacepoints = 30;
    const double max_dEdx = 50.0;
    const double min_chargecore = 10.0;    
    const double min_dist_from_fitline = 0.01;
    const double dEdxTrackDist = 5.0;

    trackFitMade = false;

    // Define the TPrincipal
    TPrincipal* principal = new TPrincipal(3,"D"); 
    // Define variables to hold the eigenvalues and eigenvectors
     const TVectorD* eigenval = new TVectorD();
     const TMatrixD* eigenvec = new TMatrixD();
     const TMatrixD* covar = new TMatrixD();
     const TVectorD* meanval = new TVectorD();
     const double *hits_spacepoints = new double[3];
     double* pca_spacepoints = new double[3];
     TVector3* pca_spacepoints_output = new TVector3();
     double *pca_spacepoints_check = new double[3];
     std::vector<TVector3>* pca_hsp_cont = new std::vector<TVector3>;

     // LOOP OVER TRACKS
     const TrackVector& tracks=evHelper.GetTracks();
     const TracksToHits& tracksToHits=evHelper.GetTracksToHits();
     const HitsToSpacePoints& hitsToSpacePoints=evHelper.GetHitsToSpacePoints();
     const MCParticleVector& mcParticles=evHelper.GetMCParticles();

     art::ServiceHandle<geo::Geometry> geom;

     calo::CalorimetryAlg CaloAlg(evHelper.GetParameterSet());

     unsigned int maxHits = 0;
     unsigned int trackIndex = 0;
     unsigned int trackMaxHits = 0;

     TVector3 recoTrackStart;
     TVector3 recoTrackEnd;

     for(auto track=tracks.begin();track!=tracks.end();++track)
       {
	 if(tracksToHits.at(*track).size() > maxHits)
	   {
	   maxHits = tracksToHits.at(*track).size();
	   trackMaxHits = trackIndex;
	   }
	 trackIndex++;
       }

     if(tracksToHits.size() >= trackMaxHits + 1)
       {          
       auto track=tracks.begin();
       for(unsigned int iterate=0; iterate<trackMaxHits; iterate++)
         ++track;

       //Track reconstruction sometimes reverses direction of track
       if((*track)->End().Z() > (*track)->Vertex().Z())
         {
	   recoTrackStart.SetX((*track)->Vertex().X());
           recoTrackStart.SetY((*track)->Vertex().Y());
	   recoTrackStart.SetZ((*track)->Vertex().Z());
	   recoTrackEnd.SetX((*track)->End().X());
           recoTrackEnd.SetY((*track)->End().Y());
	   recoTrackEnd.SetZ((*track)->End().Z());
         }
       else
         {
	   recoTrackStart.SetX((*track)->End().X());
           recoTrackStart.SetY((*track)->End().Y());
           recoTrackStart.SetZ((*track)->End().Z());
	   recoTrackEnd.SetX((*track)->Vertex().X());
	   recoTrackEnd.SetY((*track)->Vertex().Y());
	   recoTrackEnd.SetZ((*track)->Vertex().Z());
         }

       }//end of if(tracksToHits.size() >= trackMaxHits + 1)

     int nhits = 0;
     
     //Loop over hits from track with largest number of hits
     if(tracksToHits.size() >= trackMaxHits + 1)
       {
	 //Iterate to track with largest number of hits
         auto trackHits = tracksToHits.begin();
	 for(unsigned int iterate=0; iterate<trackMaxHits; iterate++)
	   ++trackHits;

	   // FILL SPACEPOINTS ARRAY WITH SPACEPOINTS BY LOOPING OVER ALL MATCHED HITS
	   for ( auto hit = (trackHits->second).begin(); hit!=(trackHits->second).end(); hit++ )
	      {
	      // CHECK FOR EXISTENCE OF KEY TO PREVENT MAP OVERRUN
	      if(hitsToSpacePoints.count(*hit)) 
	        {
		  TVector3 xyz=hitsToSpacePoints.at(*hit)->XYZ();
		  principal->AddRow( hitsToSpacePoints.at(*hit)->XYZ() );
		  nhits++;
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

	std::cout<<"Eigenvalues    "<<eval[0]<<"    "<<eval[1]<<"    "<<eval[2]<<std::endl;
	std::cout<<std::endl;
	std::cout<<"Eigenvector 1   "<<evec[0]<<"    "<<evec[3]<<"    "<<evec[6]<<std::endl;
	std::cout<<"Eigenvector 2   "<<evec[1]<<"    "<<evec[4]<<"    "<<evec[7]<<std::endl;
        std::cout<<"Eigenvector 3   "<<evec[2]<<"    "<<evec[5]<<"    "<<evec[8]<<std::endl;
	std::cout<<std::endl;

	const double *data_means = new double[3];
        data_means = meanval->GetMatrixArray();

	//Calculate angle between principal component and true direction of particle (validation check)
	auto particle=mcParticles.begin();
        double startMomNorm = sqrt((*particle)->Momentum().X() * (*particle)->Momentum().X()
			      + (*particle)->Momentum().Y() * (*particle)->Momentum().Y()
			      + (*particle)->Momentum().Z() * (*particle)->Momentum().Z());
	//Angle is arccos of scalar product between principal component and (normalised) true direction of particle
	double angleWithTrueParticle = acos((evec[0] * (*particle)->Momentum().X()
                                             + evec[3] * (*particle)->Momentum().Y()
                                             + evec[6] * (*particle)->Momentum().Z()) / startMomNorm);

	std::cout<<"Angle between principal component and true start momentum of particle = "<<angleWithTrueParticle<<" rad"<<std::endl;

        if((*particle)->EndPosition().X() > (activeVolMinX + fiducialDist) && (*particle)->EndPosition().X() < (activeVolMaxX - fiducialDist)
	   && (*particle)->EndPosition().Y() > (activeVolMinY + fiducialDist) && (*particle)->EndPosition().Y() < (activeVolMaxY - fiducialDist)
	   && (*particle)->EndPosition().Z() > (activeVolMinZ + fiducialDist) && (*particle)->EndPosition().Z() < (activeVolMaxZ - fiducialDist))
	  outputPtr->IsStoppingTrue = true;
	else
	  outputPtr->IsStoppingTrue = false;

	// WRITE TO NTUPLE
	outputPtr->EigenValues.push_back(*eigenval);
	outputPtr->EigenVectors.push_back(*eigenvec);
	outputPtr->Covariance.push_back(*covar);
	outputPtr->MeanValues.push_back(*meanval);
        outputPtr->AnglePATrueTrack.push_back(angleWithTrueParticle);
	outputPtr->EvalRatio.push_back(sqrt(eval[1] * eval[1] + eval[2] * eval[2]) / eval[0]);

	TH1D *hHitsSpacepoints = new TH1D("hHitsSpacepoints", "", 3 * nhits, 0.0, 3.0 * nhits);

	int spacept = 0;

	// *** ITERATE OVER THE HITS AGAIN TO CALCULATE THE NEW SPACEPOINTS UNDER PCA TRANSFORMATION ***
                // FILL SPACEPOINTS ARRAY WITH SPACEPOINTS BY LOOPING OVER ALL MATCHED HITS     
                for ( auto hit = (trackHits->second).begin(); hit!=(trackHits->second).end(); hit++ )
                  {
                    // CHECK FOR EXISTENCE OF KEY TO PREVENT MAP OVERRUN                        
                    if(hitsToSpacePoints.count(*hit)) {
		      // CALCULATE PCA ON SPACEPOINTS AND WRITE TO NTUPLE
		      hits_spacepoints = hitsToSpacePoints.at(*hit)->XYZ();
		      principal->X2P(hits_spacepoints,pca_spacepoints);

		      pca_spacepoints_output->SetXYZ(pca_spacepoints[0], pca_spacepoints[1], pca_spacepoints[2] );
		      pca_hsp_cont->push_back(*pca_spacepoints_output);

		      for(int i=0; i<3; i++)
		        hHitsSpacepoints->SetBinContent(3*spacept+i+1, hits_spacepoints[i]);

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
			  std::cerr<<"LArPIDCalculator::FillEventPID() PCA spacepoint validation check failed ...... exiting."<<std::endl;
			  exit(1);
			}

		      spacept++;
                    }
		  }

	outputPtr->PCAHitsSpacePoints.push_back(*pca_hsp_cont);

	//Fit line to track by minimising sum of squared residuals between hits and fitted line.
	//It would be more logical to put the hit spacepoints into a TGraph2D and pass that to the fit method.
	//However I tried this and it gave memory issues in the function CalcSumSqResidual().
	//For this reason, I put the spacepoints into a TH1D using SetBinContent() and pass that to the fit method.
	int status = this->FitTrack(hHitsSpacepoints);

	std::cout<<"Fit status = "<<status<<std::endl;

	trackFitMade = true;
	
	//Calculate scalar product of unit fitted track vector with unit norm to wire planes
	double scalarProdTrackWires[3] = {0};
        scalarProdTrackWires[0] = fittedTrackVector.Unit().Dot(normToWiresU);
        scalarProdTrackWires[1] = fittedTrackVector.Unit().Dot(normToWiresV);
        scalarProdTrackWires[2] = fittedTrackVector.Unit().Dot(normToWiresW);

	//Calculate nearest points to reconstructed start and end of track on line fitted to track
	TVector3 nearestPointStart = this->CalcNearestPointOnLine(recoTrackStart);
        TVector3 nearestPointEnd = this->CalcNearestPointOnLine(recoTrackEnd);

	double yzPitch, xComponent;
        double pitch3D = 0.45;

	double dEdxAmpStart = 0.0;
	double dEdxAmpEnd = 0.0;
        double dEdxAmpPenultimate = 0.0;
        double dEdxAmpEnd10 = 0.0;
        double dEdxAmpPenultimate10 = 0.0;
	double dEdxAmpStartDist = 0.0;
        double dEdxAmpEndDist = 0.0;
        double dEdxAmpPenultimateDist = 0.0;
	double dEdxAreaStart = 0.0;
	double dEdxAreaEnd = 0.0;
	double dEdxAreaEnd10 = 0.0;
        double dEdxAreaPenultimate = 0.0;
	double dEdxAreaPenultimate10 = 0.0;
	double dEdxAreaStartDist = 0.0;
        double dEdxAreaEndDist = 0.0;
        double dEdxAreaPenultimateDist = 0.0;

	double stdDevDist = 0.0;

	double stdDevStart = 0.0;
        double stdDevEnd = 0.0;

	int nhits_dEdx_amp_start = 0;
        int nhits_dEdx_amp_end = 0;
	int nhits_dEdx_amp_end10 = 0;
        int nhits_dEdx_amp_penultimate = 0;
        int nhits_dEdx_amp_penultimate10 = 0;
	int nhits_dEdx_amp_start_dist = 0;
        int nhits_dEdx_amp_end_dist = 0;
        int nhits_dEdx_amp_penultimate_dist = 0;
        int nhits_dEdx_area_start = 0;
        int nhits_dEdx_area_end = 0;
        int nhits_dEdx_area_end10 = 0;
        int nhits_dEdx_area_penultimate = 0;
	int nhits_dEdx_area_penultimate10 = 0;
	int nhits_dEdx_area_start_dist = 0;
        int nhits_dEdx_area_end_dist = 0;
        int nhits_dEdx_area_penultimate_dist = 0;

	int nhits_con_start = 0;
        int nhits_con_end = 0;

	double chargeCore = 0.0;
        double chargeHalo = 0.0;

	double chargeCoreStart = 0.0;
        double chargeCoreEnd = 0.0;
        double chargeHaloStart = 0.0;
        double chargeHaloEnd = 0.0;

	double chargeCon = 0.0;

        double chargeStart = 0.0;
        double chargeEnd = 0.0;
	double chargeFirstHalf = 0.0;
	double chargeSecondHalf = 0.0;
        double chargePenultimate = 0.0;
	double chargeEnd10 = 0.0;
        double chargePenultimate10 = 0.0;

	//Loop over hits again to calculate average dE/dx and shape variables
	for ( auto hit = (trackHits->second).begin(); hit!=(trackHits->second).end(); hit++ )
	  {
	    if(hitsToSpacePoints.count(*hit)) {

	      hits_spacepoints = hitsToSpacePoints.at(*hit)->XYZ();
	      principal->X2P(hits_spacepoints,pca_spacepoints);

	      TVector3 spacepoint(hits_spacepoints[0], hits_spacepoints[1], hits_spacepoints[2]);

	      double distSqFromTrackFit = this->CalcDistSqPointLine(spacepoint, fittedTrackPoint, fittedTrackVector);
	      double resRange = this->CalcResRange(spacepoint, nearestPointEnd);
	      double resRangeFrac = this->CalcResRangeFraction(spacepoint, nearestPointStart, nearestPointEnd);

	      //geom->WirePitch() is defined at http://nusoft.fnal.gov/larsoft/doxsvn/html.1.7.1/classgeo_1_1Geometry.html
	      //First 2 arguments are 0 and 1 to give true pitch between 2 adjacent wires, third argument is wire plane, fourth argument is wire TPC.
	      //True wire pitches agree with those given in DocDB 7550.
	      yzPitch = geom->WirePitch(0,1,(*hit)->WireID().Plane, (*hit)->WireID().TPC) / fabs(scalarProdTrackWires[(*hit)->WireID().Plane]);

	      xComponent = yzPitch * fittedTrackVector[0] / sqrt(fittedTrackVector[1] * fittedTrackVector[1] + fittedTrackVector[2] * fittedTrackVector[2]);
              pitch3D = sqrt(xComponent * xComponent + yzPitch * yzPitch);

	      if(nhits >= min_spacepoints && resRangeFrac > (1.0 - dEdxTrackFraction))
		{
		  if(CaloAlg.dEdx_AMP(*hit, pitch3D, evHelper.GetEventT0()) < max_dEdx)
		  {
		    nhits_dEdx_amp_start++;
		    dEdxAmpStart += CaloAlg.dEdx_AMP(*hit, pitch3D, evHelper.GetEventT0());
		  }
		  if(CaloAlg.dEdx_AREA(*hit, pitch3D, evHelper.GetEventT0()) < max_dEdx)
                  {
                    nhits_dEdx_area_start++;
                    dEdxAreaStart += CaloAlg.dEdx_AREA(*hit, pitch3D, evHelper.GetEventT0());
                  }
		nhits_con_start++;
  		stdDevStart += distSqFromTrackFit;
         	}
	      if(nhits >= min_spacepoints && resRangeFrac < dEdxTrackFraction)
		{
		  if(CaloAlg.dEdx_AMP(*hit, pitch3D, evHelper.GetEventT0()) < max_dEdx)
                  {
                    nhits_dEdx_amp_end++;
                    dEdxAmpEnd += CaloAlg.dEdx_AMP(*hit, pitch3D, evHelper.GetEventT0());
                  }
		  if(CaloAlg.dEdx_AREA(*hit, pitch3D, evHelper.GetEventT0()) < max_dEdx)
                  {
                    nhits_dEdx_area_end++;
                    dEdxAreaEnd += CaloAlg.dEdx_AREA(*hit, pitch3D, evHelper.GetEventT0());
                  }
		nhits_con_end++;
		stdDevEnd += distSqFromTrackFit;
		}
	      if(nhits >= min_spacepoints && resRangeFrac < 0.1)
                {
                  if(CaloAlg.dEdx_AMP(*hit, pitch3D, evHelper.GetEventT0()) < max_dEdx)
		    {
		      nhits_dEdx_amp_end10++;
		      dEdxAmpEnd10 += CaloAlg.dEdx_AMP(*hit, pitch3D, evHelper.GetEventT0());
		    }
                  if(CaloAlg.dEdx_AREA(*hit, pitch3D, evHelper.GetEventT0()) < max_dEdx)
		    {
		      nhits_dEdx_area_end10++;
		      dEdxAreaEnd10 += CaloAlg.dEdx_AREA(*hit, pitch3D, evHelper.GetEventT0());
		    }
                }
	      if(nhits >= min_spacepoints && resRangeFrac > dEdxTrackFraction && resRangeFrac < 2.0 * dEdxTrackFraction)
                {
                  if(CaloAlg.dEdx_AMP(*hit, pitch3D, evHelper.GetEventT0()) < max_dEdx)
		    {
		      nhits_dEdx_amp_penultimate++;
		      dEdxAmpPenultimate += CaloAlg.dEdx_AMP(*hit, pitch3D, evHelper.GetEventT0());
		    }
                  if(CaloAlg.dEdx_AREA(*hit, pitch3D, evHelper.GetEventT0()) < max_dEdx)
		    {
		      nhits_dEdx_area_penultimate++;
		      dEdxAreaPenultimate += CaloAlg.dEdx_AREA(*hit, pitch3D, evHelper.GetEventT0());
		    }
                }
	      if(nhits >= min_spacepoints && resRangeFrac > 0.1 && resRangeFrac < 0.2)
                {
                  if(CaloAlg.dEdx_AMP(*hit, pitch3D, evHelper.GetEventT0()) < max_dEdx)
                    {
                      nhits_dEdx_amp_penultimate10++;
                      dEdxAmpPenultimate10 += CaloAlg.dEdx_AMP(*hit, pitch3D, evHelper.GetEventT0());
                    }
                  if(CaloAlg.dEdx_AREA(*hit, pitch3D, evHelper.GetEventT0()) < max_dEdx)
                    {
                      nhits_dEdx_area_penultimate10++;
                      dEdxAreaPenultimate10 += CaloAlg.dEdx_AREA(*hit, pitch3D, evHelper.GetEventT0());
                    }
                }


	      if((nearestPointEnd - nearestPointStart).Mag() > 2.0 * dEdxTrackDist && (nearestPointEnd - nearestPointStart).Mag() - resRange < dEdxTrackDist)
                {
                  if(CaloAlg.dEdx_AMP(*hit, pitch3D, evHelper.GetEventT0()) < max_dEdx)
                    {
                      nhits_dEdx_amp_start_dist++;
                      dEdxAmpStartDist += CaloAlg.dEdx_AMP(*hit, pitch3D, evHelper.GetEventT0());
                    }
                  if(CaloAlg.dEdx_AREA(*hit, pitch3D, evHelper.GetEventT0()) < max_dEdx)
                    {
                      nhits_dEdx_area_start_dist++;
                      dEdxAreaStartDist += CaloAlg.dEdx_AREA(*hit, pitch3D, evHelper.GetEventT0());
                    }
                }
	      if((nearestPointEnd - nearestPointStart).Mag() > 2.0 * dEdxTrackDist && resRange < dEdxTrackDist)
                {
		  if(CaloAlg.dEdx_AMP(*hit, pitch3D, evHelper.GetEventT0()) < max_dEdx)
                    {
                      nhits_dEdx_amp_end_dist++;
                      dEdxAmpEndDist += CaloAlg.dEdx_AMP(*hit, pitch3D, evHelper.GetEventT0());
                    }
                  if(CaloAlg.dEdx_AREA(*hit, pitch3D, evHelper.GetEventT0()) < max_dEdx)
                    {
                      nhits_dEdx_area_end_dist++;
                      dEdxAreaEndDist += CaloAlg.dEdx_AREA(*hit, pitch3D, evHelper.GetEventT0());
                    }
		}
	      if((nearestPointEnd - nearestPointStart).Mag() > 2.0 * dEdxTrackDist && resRange > dEdxTrackDist && resRange < 2.0 * dEdxTrackDist)
                {
                  if(CaloAlg.dEdx_AMP(*hit, pitch3D, evHelper.GetEventT0()) < max_dEdx)
                    {
                      nhits_dEdx_amp_penultimate_dist++;
                      dEdxAmpPenultimateDist += CaloAlg.dEdx_AMP(*hit, pitch3D, evHelper.GetEventT0());
                    }
                  if(CaloAlg.dEdx_AREA(*hit, pitch3D, evHelper.GetEventT0()) < max_dEdx)
                    {
                      nhits_dEdx_area_penultimate_dist++;
		      dEdxAreaPenultimateDist += CaloAlg.dEdx_AREA(*hit, pitch3D, evHelper.GetEventT0());
                    }
                }


	      stdDevDist += distSqFromTrackFit;

	      if(sqrt(distSqFromTrackFit) < MoliereRadiusFraction * MoliereRadius)
		chargeCore += (*hit)->Integral();
	      else
		chargeHalo += (*hit)->Integral();

	      //Avoid division by zero by requiring a minimum distance from fit line
	      chargeCon += (*hit)->Integral() / TMath::Max(min_dist_from_fitline, sqrt(distSqFromTrackFit));

	      if(nhits >= min_spacepoints)
                {
		  if(resRangeFrac > (1.0 - coreHaloFraction))
		    {
		    if(sqrt(distSqFromTrackFit) < MoliereRadiusFraction * MoliereRadius)
		      chargeCoreStart += (*hit)->Integral();
		    else
		      chargeHaloStart += (*hit)->Integral();
		    }

		  if(resRangeFrac < coreHaloFraction)
		    {
		    if(sqrt(distSqFromTrackFit) < MoliereRadiusFraction * MoliereRadius)
		      chargeCoreEnd += (*hit)->Integral();
		    else
		      chargeHaloEnd += (*hit)->Integral();
		    }

		  if(resRangeFrac > (1.0 - chargeTrackFraction))
		    chargeStart += (*hit)->Integral();
		  if(resRangeFrac < chargeTrackFraction)
		    chargeEnd += (*hit)->Integral();
                  if(resRangeFrac > 0.5)
		    chargeFirstHalf += (*hit)->Integral();
		  if(resRangeFrac < 0.5)
                    chargeSecondHalf += (*hit)->Integral();

		  if(resRangeFrac > chargeTrackFraction && resRangeFrac < 2.0 * chargeTrackFraction)
                    chargePenultimate += (*hit)->Integral();

		  if(resRangeFrac < 0.1)
                    chargeEnd10 += (*hit)->Integral();
                  if(resRangeFrac > 0.1 && resRangeFrac < 0.2)
                    chargePenultimate10 += (*hit)->Integral();

		}

	    }
	  }

	if(recoTrackEnd.X() > (activeVolMinX + fiducialDist) && recoTrackEnd.X() < (activeVolMaxX - fiducialDist)
	   && recoTrackEnd.Y() > (activeVolMinY + fiducialDist) && recoTrackEnd.Y() < (activeVolMaxY - fiducialDist)
	   && recoTrackEnd.Z() > (activeVolMinZ + fiducialDist) && recoTrackEnd.Z() < (activeVolMaxZ - fiducialDist))
	  outputPtr->IsStoppingReco = true;
	else
	  outputPtr->IsStoppingReco = false;



	if(nhits >= min_spacepoints && nhits_dEdx_amp_start >= 1)
	  outputPtr->AvgedEdxAmpStart.push_back(dEdxAmpStart / nhits_dEdx_amp_start);
	else
	  outputPtr->AvgedEdxAmpStart.push_back(0.0);

	if(nhits >= min_spacepoints && nhits_dEdx_amp_end >= 1)
          outputPtr->AvgedEdxAmpEnd.push_back(dEdxAmpEnd / nhits_dEdx_amp_end);
        else
          outputPtr->AvgedEdxAmpEnd.push_back(0.0);

	if(nhits >= min_spacepoints && nhits_dEdx_amp_start >= 1 && nhits_dEdx_amp_end >= 1)
          outputPtr->AvgedEdxAmpLongRatio.push_back((dEdxAmpEnd / nhits_dEdx_amp_end) / (dEdxAmpStart / nhits_dEdx_amp_start));
        else
          outputPtr->AvgedEdxAmpLongRatio.push_back(1.0);

        if(nhits >= min_spacepoints && nhits_dEdx_amp_end >= 1 && nhits_dEdx_amp_penultimate >= 1)
          outputPtr->AvgedEdxAmpEndRatio.push_back((dEdxAmpEnd / nhits_dEdx_amp_end) / (dEdxAmpPenultimate / nhits_dEdx_amp_penultimate));
        else
          outputPtr->AvgedEdxAmpEndRatio.push_back(1.0);

	if(nhits >= min_spacepoints && nhits_dEdx_amp_end10 >= 1 && nhits_dEdx_amp_penultimate10 >= 1)
          outputPtr->AvgedEdxAmpEndRatio10.push_back((dEdxAmpEnd10 / nhits_dEdx_amp_end10) / (dEdxAmpPenultimate10 / nhits_dEdx_amp_penultimate10));
        else
          outputPtr->AvgedEdxAmpEndRatio10.push_back(1.0);

	if((nearestPointEnd - nearestPointStart).Mag() > 2.0 * dEdxTrackDist && nhits_dEdx_amp_start_dist >= 1)
          outputPtr->AvgedEdxAmpStartDist.push_back(dEdxAmpStartDist / nhits_dEdx_amp_start_dist);
	else
          outputPtr->AvgedEdxAmpStartDist.push_back(0.0);

	if((nearestPointEnd - nearestPointStart).Mag() > 2.0 * dEdxTrackDist && nhits_dEdx_amp_end_dist >= 1)
          outputPtr->AvgedEdxAmpEndDist.push_back(dEdxAmpEndDist / nhits_dEdx_amp_end_dist);
	else
          outputPtr->AvgedEdxAmpEndDist.push_back(0.0);

	if((nearestPointEnd - nearestPointStart).Mag() > 2.0 * dEdxTrackDist && nhits_dEdx_amp_start_dist >= 1 && nhits_dEdx_amp_end_dist >= 1)
          outputPtr->AvgedEdxAmpLongRatioDist.push_back((dEdxAmpEndDist / nhits_dEdx_amp_end_dist) / (dEdxAmpStartDist / nhits_dEdx_amp_start_dist));
        else
	  outputPtr->AvgedEdxAmpLongRatioDist.push_back(1.0);

        if((nearestPointEnd - nearestPointStart).Mag() > 2.0 * dEdxTrackDist && nhits_dEdx_amp_end_dist >= 1 && nhits_dEdx_amp_penultimate_dist >= 1)
          outputPtr->AvgedEdxAmpEndRatioDist.push_back((dEdxAmpEndDist / nhits_dEdx_amp_end_dist) / (dEdxAmpPenultimateDist / nhits_dEdx_amp_penultimate_dist));
        else
          outputPtr->AvgedEdxAmpEndRatioDist.push_back(1.0);



	if(nhits >= min_spacepoints && nhits_dEdx_area_start >= 1)
          outputPtr->AvgedEdxAreaStart.push_back(dEdxAreaStart / nhits_dEdx_area_start);
        else
          outputPtr->AvgedEdxAreaStart.push_back(0.0);

        if(nhits >= min_spacepoints && nhits_dEdx_area_end >= 1)
          outputPtr->AvgedEdxAreaEnd.push_back(dEdxAreaEnd / nhits_dEdx_area_end);
        else
          outputPtr->AvgedEdxAreaEnd.push_back(0.0);

	if(nhits >= min_spacepoints && nhits_dEdx_area_start >= 1 && nhits_dEdx_area_end >= 1)
          outputPtr->AvgedEdxAreaLongRatio.push_back((dEdxAreaEnd / nhits_dEdx_area_end) / (dEdxAreaStart / nhits_dEdx_area_start));
        else
          outputPtr->AvgedEdxAreaLongRatio.push_back(1.0);

	if(nhits >= min_spacepoints && nhits_dEdx_area_end >= 1 && nhits_dEdx_area_penultimate >= 1)
          outputPtr->AvgedEdxAreaEndRatio.push_back((dEdxAreaEnd / nhits_dEdx_area_end) / (dEdxAreaPenultimate / nhits_dEdx_area_penultimate));
        else
          outputPtr->AvgedEdxAreaEndRatio.push_back(1.0);

        if(nhits >= min_spacepoints && nhits_dEdx_area_end10 >= 1 && nhits_dEdx_area_penultimate10 >= 1)
          outputPtr->AvgedEdxAreaEndRatio10.push_back((dEdxAreaEnd10 / nhits_dEdx_area_end10) / (dEdxAreaPenultimate10 / nhits_dEdx_area_penultimate10));
        else
          outputPtr->AvgedEdxAreaEndRatio10.push_back(1.0);

	if((nearestPointEnd - nearestPointStart).Mag() > 2.0 * dEdxTrackDist && nhits_dEdx_area_start_dist >= 1)
          outputPtr->AvgedEdxAreaStartDist.push_back(dEdxAreaStartDist / nhits_dEdx_area_start_dist);
        else
          outputPtr->AvgedEdxAreaStartDist.push_back(0.0);

        if((nearestPointEnd - nearestPointStart).Mag() > 2.0 * dEdxTrackDist && nhits_dEdx_area_end_dist >= 1)
          outputPtr->AvgedEdxAreaEndDist.push_back(dEdxAreaEndDist / nhits_dEdx_area_end_dist);
        else
          outputPtr->AvgedEdxAreaEndDist.push_back(0.0);

        if((nearestPointEnd - nearestPointStart).Mag() > 2.0 * dEdxTrackDist && nhits_dEdx_area_start_dist >= 1 && nhits_dEdx_area_end_dist >= 1)
          outputPtr->AvgedEdxAreaLongRatioDist.push_back((dEdxAreaEndDist / nhits_dEdx_area_end_dist) / (dEdxAreaStartDist / nhits_dEdx_area_start_dist));
	else
          outputPtr->AvgedEdxAreaLongRatioDist.push_back(1.0);

        if((nearestPointEnd - nearestPointStart).Mag() > 2.0 * dEdxTrackDist && nhits_dEdx_area_end_dist >= 1 && nhits_dEdx_area_penultimate_dist >= 1)
          outputPtr->AvgedEdxAreaEndRatioDist.push_back((dEdxAreaEndDist / nhits_dEdx_area_end_dist) / (dEdxAreaPenultimateDist / nhits_dEdx_area_penultimate_dist));
        else
          outputPtr->AvgedEdxAreaEndRatioDist.push_back(1.0);



	if(nhits >= 2)
          outputPtr->StdDevDistFromFitLine.push_back(sqrt(stdDevDist / (nhits-1)));
        else
          outputPtr->StdDevDistFromFitLine.push_back(0.0);

	if(nhits >= min_spacepoints && nhits_con_start >= 2)
	  stdDevStart = sqrt(stdDevStart / (nhits_con_start - 1));

	if(nhits >= min_spacepoints && nhits_con_end >= 2)
	  stdDevEnd = sqrt(stdDevEnd / (nhits_con_end - 1));

	if(nhits < min_spacepoints || nhits_con_start <= 1 || nhits_con_end <= 1)
          outputPtr->Conicalness.push_back(1.0);
	else
	  outputPtr->Conicalness.push_back(stdDevEnd / stdDevStart);

	if(chargeCore < min_chargecore)
          outputPtr->ChargeRatioCoreHalo.push_back(1.0);
        else
          outputPtr->ChargeRatioCoreHalo.push_back(chargeHalo / chargeCore);

        outputPtr->Concentration.push_back(chargeCon);

	if(nhits < min_spacepoints)
	  {
	  outputPtr->ChargeLongRatio.push_back(1.0);
	  outputPtr->ChargeLongRatioHalf.push_back(1.0);
	  outputPtr->ChargeEndRatio.push_back(1.0);
	  outputPtr->ChargeEndRatio10.push_back(1.0);
	  }
	else
	  {
	  if(chargeStart == 0.0)
	    outputPtr->ChargeLongRatio.push_back(1.0);
	  else
	    outputPtr->ChargeLongRatio.push_back(chargeEnd / chargeStart);
	  if(chargeFirstHalf == 0.0)
	    outputPtr->ChargeLongRatioHalf.push_back(1.0);
	  else
            outputPtr->ChargeLongRatioHalf.push_back(chargeSecondHalf / chargeFirstHalf);
	  if(chargePenultimate == 0.0)
	    outputPtr->ChargeEndRatio.push_back(1.0);
	  else
            outputPtr->ChargeEndRatio.push_back(chargeEnd / chargePenultimate);
	  if(chargePenultimate10 == 0.0)
	    outputPtr->ChargeEndRatio10.push_back(1.0);
	  else
            outputPtr->ChargeEndRatio10.push_back(chargeEnd10 / chargePenultimate10);
	  }

	delete hHitsSpacepoints;

       }//end of if(tracksToHits.size() >= 1) 

     // Clear the data from the TPrincipal
     principal->Clear();

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

 double LArPIDCalculator::CalcDistSqPointLine(const TVector3& point, const TVector3& linePoint, const TVector3& lineDirection)
  {
    //Calculate the square of the shortest distance between a point and a line.
    //This distance is given by |AC x AB| / |AB|, where A is a point on the line, AB is the direction of the line, C is the point 
    //and "x" is the cross product. This is method 2 in http://www.qc.edu.hk/math/Advanced%20Level/Point_to_line.htm
    
    TVector3 unitDirection = lineDirection.Unit();
    double distanceSq = ((point - linePoint).Cross(unitDirection)).Mag2();
    
    return distanceSq;    
  }

 TVector3 LArPIDCalculator::CalcNearestPointOnLine(const TVector3& point)
  {

    if(!trackFitMade)
      { 
	std::cerr<<"LArPIDCalculator::CalcNearestPointOnLine() is called before the fit of the track is made. Please make this call after the track is fitted."<<std::endl;
        exit(1);
      }

    //Calculate the nearest point on a line from a given point in 3D space
    //This is adapted from method 3 in http://www.qc.edu.hk/math/Advanced%20Level/Point_to_line.htm 
    //The parametric equation of the line is A + t*B where A corresponds to fittedTrackPoint, B to fittedTrackVector and t is a parameter
    //Let the given point be C and the nearest point on the line be D; then CD must be perpendicular to B
    //D = A + x*B
    //CD = D - C = A + x*B - C
    //Solving B . (A + x*B - C) = 0 for x gives:  
    double x = (-fittedTrackVector.X() * fittedTrackPoint.X() + fittedTrackVector.X() * point.X() - fittedTrackVector.Y() * fittedTrackPoint.Y()
		+ fittedTrackVector.Y() * point.Y() - fittedTrackVector.Z() * fittedTrackPoint.Z() + fittedTrackVector.Z() * point.Z())
      / (fittedTrackVector.X() * fittedTrackVector.X() + fittedTrackVector.Y() * fittedTrackVector.Y() + fittedTrackVector.Z() * fittedTrackVector.Z());

    TVector3 nearestPoint = fittedTrackPoint + x * fittedTrackVector;

    return nearestPoint;
  }

  double LArPIDCalculator::CalcResRange(const TVector3& point, const TVector3& trackEnd)
  {
    if(!trackFitMade)
      {
	std::cerr<<"LArPIDCalculator::CalcResRangeFraction() is called before the fit of the track is made. Please make this call after the track is fitted."<<std::endl;
        exit(1);
      }

    //Calculate residual range as a distance; this range is calculated along the line fitted to the track.
    TVector3 nearestPoint = this->CalcNearestPointOnLine(point);
    double resRange = (trackEnd - nearestPoint).Mag();

    return resRange;
  }

  double LArPIDCalculator::CalcResRangeFraction(const TVector3& point, const TVector3& trackStart, const TVector3& trackEnd)
  {
    if(!trackFitMade)
      {
	std::cerr<<"LArPIDCalculator::CalcResRangeFraction() is called before the fit of the track is made. Please make this call after the track is fitted."<<std::endl;
        exit(1);
      }

    //Calculate residual range as a fraction of track length; this range is calculated along the line fitted to the track.
    TVector3 nearestPoint = this->CalcNearestPointOnLine(point);
    double resRangeFraction = (trackEnd - nearestPoint).Mag() / (trackEnd - trackStart).Mag();

    return resRangeFraction;
  }

  int LArPIDCalculator::FitTrack(TH1D *hHitsSpacepoints)
  {

    const int nparams = 6;
    
    TVirtualFitter::SetDefaultFitter("Minuit");
    TVirtualFitter * fitter = TVirtualFitter::Fitter(0, nparams);
    assert(fitter);

    fitter->SetObjectFit(hHitsSpacepoints);

    //Set print level
    Double_t arglist[5];
    arglist[0] = 0;
    fitter->ExecuteCommand("SET PRINT",arglist,1);

    fitter->Clear();

    //Set fit parameters
    //Arguments are parameter index, parameter name, start value, step, minimum allowed value, maximum allowed value     
    fitter->SetParameter(0, "TrackPointX", 100.0, 1.0, activeVolMinX, activeVolMaxX);
    fitter->SetParameter(1, "TrackPointY", 10.0, 1.0, activeVolMinY, activeVolMaxY);    
    fitter->SetParameter(2, "TrackPointZ", 75.0, 1.0, activeVolMinZ, activeVolMaxZ);
    fitter->SetParameter(3, "TrackVectorX", 0.0, 0.2, -100.0, 100.0);
    fitter->SetParameter(4, "TrackVectorY", 0.0, 0.2, -100.0, 100.0);
    //Restrict fitted z direction to positive values. This is a cheating way of resolving the ambiguity in the fitted direction
    //(there are 2 possible best-fit values of direction which are exactly opposite to each other).
    //fitter->SetParameter(5, "TrackVectorZ", 10.0, 0.2, -100.0, 100.0);
    fitter->SetParameter(5, "TrackVectorZ", 10.0, 0.2, 0.0, 100.0);

    //Set minimization function 
    fitter->SetFCN(CalcSumSqResidual);

    // MINUIT's minimization step
    arglist[0] = 5000; // number of FCN calls
    arglist[1] = 0.01; // tolerance
    int migrad_status = fitter->ExecuteCommand("MIGRAD",arglist, 2);

    std::cout<<"MIGRAD fit status = "<<migrad_status<<std::endl;

    if(migrad_status != 0)
      {
      std::cout<<"MIGRAD fit failed, calling SIMPLEX...."<<std::endl;
      arglist[1] = 0.1;  //increase tolerance
      int simplex_status = fitter->ExecuteCommand("SIMPLEX",arglist, 2);
      std::cout<<"SIMPLEX fit status = "<<simplex_status<<std::endl;
      }

    TVector3 fitDirection(fitter->GetParameter(3), fitter->GetParameter(4), fitter->GetParameter(5));
    std::cout<<"Best-fit direction    "<<fitDirection.Unit().X()<<"    "<<fitDirection.Unit().Y()<<"    "<<fitDirection.Unit().Z()<<std::endl;

    fittedTrackPoint.SetX(fitter->GetParameter(0));
    fittedTrackPoint.SetY(fitter->GetParameter(1));
    fittedTrackPoint.SetZ(fitter->GetParameter(2));
    fittedTrackVector.SetX(fitter->GetParameter(3));
    fittedTrackVector.SetY(fitter->GetParameter(4));
    fittedTrackVector.SetZ(fitter->GetParameter(5));

    return migrad_status;

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
    //const TracksToCalo& tracksToCalo=evHelper.GetTracksToCalo();

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

  }//end of void LArPIDCalculator::FillEventTracksAndHits()

  void CalcSumSqResidual(Int_t &npar, Double_t *, Double_t &f, Double_t *par, Int_t iflag)
  {
    /*
    This function is required to have 5 arguments as above. These arguments are
    1. Number of parameters allowed to vary in the fit
    2. Array of gradients or first derivatives of parameters - input or output or both ?
    3. Calculated function value - in this case it is the sum of squares of residuals between track hits and fit to them
    4. Array of fit parameters - includes both fixed and variable parameters
    5. Indicates what is to be calculated

    This information is taken from http://root.cern.ch/root/html/TMinuit.html#TMinuit:Eval

    It is a bit vague about the meaning of argument 5
    */

    TVirtualFitter * fitter = TVirtualFitter::GetFitter();
    TH1D * hit_spacepoints = (TH1D*) fitter->GetObjectFit();
    if(!hit_spacepoints) 
      {
	std::cout << " *** No hit spacepoints! ***" << std::endl;
	exit(1);
      }
    
    double trackPointX = 0.0;
    double trackPointY = 0.0;
    double trackPointZ = 0.0;
    double trackVectorX = 0.0;
    double trackVectorY = 0.0;
    double trackVectorZ = 0.0;

    trackPointX = par[0];
    trackPointY = par[1];
    trackPointZ = par[2];
    trackVectorX = par[3];
    trackVectorY = par[4];
    trackVectorZ = par[5];

    const TVector3 trackPoint(trackPointX, trackPointY, trackPointZ);
    const TVector3 trackVector(trackVectorX, trackVectorY, trackVectorZ);

    int npoints = hit_spacepoints->GetNbinsX()/3;
    const int nspacepoints = npoints;

    lar_valrec::LArPIDCalculator LArPIDCalc;

    double sumSqResidual = 0.0;

    for(int point=0; point<nspacepoints; point++)
      {
	double X = hit_spacepoints->GetBinContent(3*point+1);
        double Y = hit_spacepoints->GetBinContent(3*point+2);
	double Z = hit_spacepoints->GetBinContent(3*point+3);

        const TVector3 spacepoint(X, Y, Z);

        sumSqResidual += LArPIDCalc.CalcDistSqPointLine(spacepoint, trackPoint, trackVector);
      }

    f = sumSqResidual;

  }

}//end of namespace lar_valrec
