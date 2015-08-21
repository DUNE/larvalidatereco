// Comparison between the transformed PCA spacepoints and the initial spacepoints, and an initial attempt to reconstruct the PCA via the covariance matrix from the TPrincipal

#include <string>
#include "TSystem.h"
#include "TCanvas.h"
#include "TObject.h"
#include "data_objects/LArPID.h"
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
std::vector<std::vector<TVector3>> FillSpacePointsPCAreco(const char* filename)
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
std::vector<std::vector<TVector3>> FillTracksPCAreco(const char* filename)
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

std::vector<TMatrixD> FillEigenVectorsPCAreco(const char* filename)
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



TVector3 FindTrackEnd(std::vector<TVector3> track)

{
  TVector3* max = new TVector3();
  for ( auto hit_pntr = track.begin() ; hit_pntr!=track.end(); hit_pntr++ )
    {

      if ( hit_pntr->X() > max->X() ) { *max = *hit_pntr; }
    }

  return *max;
}
// Calculate the minimum track coordinate along the principal axis
TVector3 FindTrackStart(std::vector<TVector3> track)

{
  TVector3 *min = new TVector3();
  for ( auto hit_pntr = track.begin() ; hit_pntr!=track.end(); hit_pntr++ )
    {

      if ( hit_pntr->X() < min->X() ) { *min = *hit_pntr; }
    }

  return *min;
}



void LArPIDCheckPCAreco()

{
  // Fill vectors 
  std::vector<std::vector<TVector3>> trackcon = FillTracksPCAreco("/data/t2k/phumhh/lbne/ntuple_1gev_muon.root");
  std::vector<std::vector<TVector3>> spacepointcon = FillSpacePointsPCAreco("/data/t2k/phumhh/lbne/ntuple_1gev_muon.root");
  std::vector<TMatrixD> evcon = FillEigenVectorsPCAreco("/data/t2k/phumhh/lbne/ntuple_1gev_muon.root");
 
  // Load in first values of each
  TVector3 orig = spacepointcon[0][0];
  TVector3 pca = trackcon[0][0];
  TMatrixD trans = evcon[0];
  // Declare variables for first and last trackhit positions
  TVector3 *pcatrackstart = new TVector3();
  TVector3 *pcatrackend = new TVector3();
  TVector3 *sptrackstart = new TVector3();
  TVector3 *sptrackend = new TVector3();
  TVector3 *pca_diff_holding = new TVector3();
  TVector3 *sp_diff_holding = new TVector3();
  // Declare vectors to hold the enpoint coordiantes
  std::vector<TVector3>* pcamin = new std::vector<TVector3>;
  std::vector<TVector3>* pcamax = new std::vector<TVector3>;
  std::vector<TVector3>* spmin = new std::vector<TVector3>;
  std::vector<TVector3>* spmax = new std::vector<TVector3>;
  std::vector<TVector3>* pcadiff = new std::vector<TVector3>;
  std::vector<TVector3>* spdiff = new std::vector<TVector3>;
  
  // iterate over tracks for pca values
  for (auto i = trackcon.begin(); i!=trackcon.end(); i++ )
    {	  
	  // Read out track vector in terms of pca and the position of its endpoints
	  *pcatrackstart = FindTrackStart(*i);
	  *pcatrackend = FindTrackEnd(*i);
	  pcamin->push_back(*pcatrackstart);
	  pcamax->push_back(*pcatrackend);
	  *pca_diff_holding = *pcatrackend-*pcatrackstart;
	  // Normalise the path difference
	  Double_t max = pca_diff_holding->Mag();
	  pca_diff_holding->SetX(pca_diff_holding->X() / max);
	  pca_diff_holding->SetY(pca_diff_holding->Y() / max);
	  pca_diff_holding->SetZ(pca_diff_holding->Z() / max);
	  pcadiff->push_back(*pca_diff_holding);
    }
  cout << "Length of PCA Track vector: " << pcadiff->size() << endl;;
  // iterate over tracks for spacepoint values
  for (auto j = spacepointcon.begin(); j!=spacepointcon.end(); j++)
    {
     
      // Read out track vector in terms of spacepoints and the positon of its endpoints
      *sptrackstart = FindTrackStart(*j);
      *sptrackend = FindTrackEnd(*j);
      spmin->push_back(*sptrackstart);
      spmax->push_back(*sptrackend);
      *sp_diff_holding = *sptrackend-*sptrackstart;
      // Normalise the path difference
      Double_t max = sp_diff_holding->Mag();
      sp_diff_holding->SetX(sp_diff_holding->X() / max);
      sp_diff_holding->SetY(sp_diff_holding->Y() / max);
      sp_diff_holding->SetZ(sp_diff_holding->Z() / max);
      spdiff->push_back(*sp_diff_holding);
    }
  cout << "Length of Spacepoints Track vector: " << spdiff->size() << endl;

  // CHECK READOUT
  for (unsigned int m = 0; m!=spdiff->size(); m++ ) { if(m%20==0){cout<<"SP |  " << spdiff->at(m).Z() << endl ;cout<<"PCA |  " <<  pcadiff->at(m).Z() << endl;} }


  std::vector<Int_t> x_match; std::vector<Int_t> y_match; std::vector<Int_t> z_match;
  // iterate over the endpoints, find their differences, and match axes
  for (auto k = spdiff->begin(); k!=spdiff->end(); k++)
    {
      for (auto l = pcadiff->begin(); l!=pcadiff->end(); l++)
	{
	  // If spacepoint difference matches pca difference
	  if ( k->Z() == l->X() ) { x_match.push_back(1); }
	  if ( k->Z() == l->Y() ) { y_match.push_back(1); }
	  if ( k->Z() == l->Z() ) { z_match.push_back(1); }
	}	  
    }
  cout << endl << "Tracks matched Z-X: " << x_match.size() << " Tracks matched Z-Y: " << y_match.size() << " Tracks matched Z-Z: " << z_match.size() << endl;




  /*
  TMatrixD* pcamat = new TMatrixD(3,1);
  Double_t* pca_cont = new Double_t[3];
  pca_cont[0] = pca.X(); pca_cont[1] = pca.Y(); pca_cont[2] = pca.Z();
  pcamat->Allocate(0,2,0,0);
  pcamat->Use(0,2,0,0,pca_cont);

  // Multiply the original vector by the rotation matrix to obtain the new coordinates
  TMatrixD* end = new TMatrixD(1,3);
  TMatrixD* transpose = new TMatrixD(3,3);
  transpose->Transpose(trans);
  end->Mult(*transpose,*origmat);
  end->Print();
  pca.Print();
  Double_t* newvector = new Double_t[3];
  newvector = end->GetMatrixArray();
  TVector3 finalvec(newvector[0],newvector[1],newvector[2]);
  Double_t angle = finalvec.Angle(orig);

  // Stored pca direction cosines
  pca_cos_x = 180*(std::acos(pca.X() / pca.Mag()))/3.1415926;
  pca_cos_y = 180*(std::acos(pca.Y() / pca.Mag()))/3.1415926;
  pca_cos_z = 180*(std::acos(pca.Z() / pca.Mag()))/3.1415926;
  */
}
