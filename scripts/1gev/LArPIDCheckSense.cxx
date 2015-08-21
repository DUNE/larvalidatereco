// Check for the sense of the tracks output by the PCA, ensuring they are aligned along the positive Z axis in the lab frame.

#include <string>
#include "TSystem.h"
#include "TCanvas.h"
#include "TObject.h"
#include "../data_objects/LArPID.h"
#include "iostream"
#include "TTree.h"
#include "TFile.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TH2D.h"
#include "TAxis.h"
#include <vector>
#include "TLorentzVector.h"
#include </data/t2k/phumhh/lbne/srcs/larvalidatereco/larvalidatereco/content/data_objects/LArAnalysis.h>
#include </data/t2k/phumhh/lbne/srcs/larvalidatereco/larvalidatereco/content/data_objects/LArValidation.h>

using namespace std;



// Fill vector of XYZ spacepoints from ntuple
std::vector<std::vector<TVector3>> FillSpacePointsSense(const char* filename)
{
  TFile *ntuple_PID = new TFile(filename);
  ntuple_PID->cd();
  if ( !ntuple_PID ) { std::cout << std::endl << "Unable to open muon ntuple file" << std::endl; }
  else{ std::cout << std::endl <<  "ntuple open: " << filename << std::endl; }
  TDirectoryFile *valrecfolder = (TDirectoryFile*)ntuple_PID->Get("valrec");
  TTree *ValRecTree = (TTree*)valrecfolder->Get("valrec");
  LArAnalysis* ana=0;
  ValRecTree->SetBranchAddress("Analysis",&ana);
  if ( ValRecTree && valrecfolder ) { std::cout << "data retrieved" << std::endl; }
  Long64_t nentries = ValRecTree->GetBranch("Hit3Pos")->GetEntries();
  std::cout << std::endl << "Entries read: " << nentries << std::endl;
  std::vector<std::vector<TVector3>>* spacepoint_con = new std::vector<std::vector<TVector3>>;
  std::vector<TVector3>* hit_con = new std::vector<TVector3>;

  // Fill vector of original hit spacepoints
  std::cout << "calculating" << std::endl;
  for(int i=0;i<nentries;++i)
    {
      ValRecTree->GetBranch("Hit3Pos")->GetEntry(i);
      if ( ana->Hit3Pos.size() )
        {
          if ( i%500==0) { std::cout << "Event: " << i << std::endl;}
          for (auto j = ana->Hit3Pos.begin(); j != ana->Hit3Pos.end(); j++ )
            {
                  //std::cout << "Event: " << i << " | Track: " << track << " | Hit: " << hit << endl;
                  hit_con->push_back(*j);
                  //cout << "X: " << k->X() << " Y: " << k->Y() << " Z: " << k->Z() << endl 
            }

        }
      spacepoint_con->push_back(*hit_con);
      hit_con->clear();
    }

  ntuple_PID->Close();
  return *spacepoint_con;
}

 // Fill track objects from ntuple
std::vector<std::vector<TVector3>> FillTracksSense(const char* filename)
{
  TFile *ntuple_PID = new TFile(filename);
  ntuple_PID->cd();
  if ( !ntuple_PID ) { std::cout << std::endl << "Unable to open muon ntuple file" << std::endl; }
  else{ std::cout << std::endl <<  "ntuple open: " << filename << std::endl; }
  TDirectoryFile *valrecfolder = (TDirectoryFile*)ntuple_PID->Get("valrec");
  TTree *ValRecTree = (TTree*)valrecfolder->Get("valrec");
  LArPID* pid=0;
  LArAnalysis* ana=0;
  ValRecTree->SetBranchAddress("PID",&pid);
  ValRecTree->SetBranchAddress("Analysis",&ana);
  if ( ValRecTree && valrecfolder ) { std::cout << "data retrieved" << std::endl; }
  Long64_t nentries = ValRecTree->GetBranch("PCAHitsSpacePoints")->GetEntries();
  std::cout << std::endl << "Entries read: " << nentries << std::endl;
  std::vector<std::vector<TVector3>>* track_con = new std::vector<std::vector<TVector3>>;
  std::vector<TVector3>* hit_con = new std::vector<TVector3>;

  // Fill vector of pca transformed tracks
  std::cout << "calculating" << std::endl;
  for(int i=0;i<nentries;++i)
    {
      ValRecTree->GetBranch("PCAHitsSpacePoints")->GetEntry(i);
      if ( pid->PCAHitsSpacePoints.size() && ana->NHits )
	{
	  if ( i%500==0) { std::cout << "Event: " << i << std::endl;}
	 
	  for (auto j = pid->PCAHitsSpacePoints.begin(); j != pid->PCAHitsSpacePoints.end(); j++ )
	    {
	      for ( auto k= j->begin(); k!=j->end(); k++ )
		{

		  //std::cout << "Event: " << i << " | Track: " << track << " | Hit: " << hit << endl;
	      hit_con->push_back(*k);
		  //cout << "X: " << k->X() << " Y: " << k->Y() << " Z: " << k->Z() << endl;
		}
	    }

	}
      track_con->push_back(*hit_con);
      hit_con->clear();
	   
    }

  ntuple_PID->Close();
  return *track_con;
}

std::vector<TMatrixD> FillEigenVectorsSense(const char* filename)
{
  TFile *ntuple_PID = new TFile(filename);
  ntuple_PID->cd();
  if ( !ntuple_PID ) { std::cout << std::endl << "Unable to open muon ntuple file" << std::endl; }
  else{ std::cout << std::endl <<  "ntuple open: " << filename << std::endl; }
  TDirectoryFile *valrecfolder = (TDirectoryFile*)ntuple_PID->Get("valrec");
  TTree *ValRecTree = (TTree*)valrecfolder->Get("valrec");
  LArPID* pid=0;
  LArAnalysis* ana=0;
  ValRecTree->SetBranchAddress("PID",&pid);
  ValRecTree->SetBranchAddress("Analysis",&ana);
  if ( ValRecTree && valrecfolder ) { std::cout << "data retrieved" << std::endl; }
  Long64_t nentries = ValRecTree->GetBranch("EigenVectors")->GetEntries();
  std::cout << std::endl << "Entries read: " << nentries << std::endl;
  std::vector<TMatrixD>* track_con = new std::vector<TMatrixD>;

  // Fill vector of eigenvector matrices
  std::cout << "calculating" << std::endl;
  for(int i=0;i<nentries;i++)
    {
      ValRecTree->GetBranch("EigenVectors")->GetEntry(i);
      if ( pid->EigenVectors.size() )
	{
	  for (auto j = pid->EigenVectors.begin(); j!=pid->EigenVectors.end(); j++ )
	    {
	      if ( i%500==0) { std::cout << "Event: " << i << std::endl;}
	      track_con->push_back(*j);    	  
	    }
	}
    }
  ntuple_PID->Close();
  return *track_con;
}




void LArPIDCheckSense()

{
  // Fill vectors 
  std::vector<std::vector<TVector3>> trackcon = FillTracksSense("/data/t2k/phumhh/lbne/ntuple_1gev_muon.root");
  std::vector<std::vector<TVector3>> spacepointcon = FillSpacePointsSense("/data/t2k/phumhh/lbne/ntuple_1gev_muon.root");
  std::vector<TMatrixD> evcon = FillEigenVectorsSense("/data/t2k/phumhh/lbne/ntuple_1gev_muon.root");
  Int_t tracknum= 0;
  Double_t prev = 0; Double_t errnum = 0; Double_t hitnum = 0;
  // iterate over tracks for pca values
  for (auto i = trackcon.begin(); i!=trackcon.end(); i++ )
    {

      for (auto j = i->begin(); j!=i->end(); j++)
	{
	  if ( j->X() < prev && (TMath::Abs(prev-j->X()) < 0.2) ) { errnum++; }

	  hitnum++;
	}
      tracknum++;
      prev = 0;
    }
  cout << "\tTracks read: " << tracknum << endl;
  cout << "\tErrors: " << (errnum/hitnum)*100 << "%" << endl;
  
  for (auto j = evcon.begin(); j!= evcon.end(); j++ ) 
    {
      j->Print();
    }

}
