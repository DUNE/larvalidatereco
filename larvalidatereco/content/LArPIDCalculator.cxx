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
  }

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
