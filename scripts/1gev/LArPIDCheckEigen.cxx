// An initial check for the validity of the TPrincipal PCA output, printing the eigenvalues of any badly transformed tracks by the PCA

#include <string>
#include "TSystem.h"
#include "TCanvas.h"
#include "TObject.h"
#include "../data_objects/LArPID.h"
#include "iostream"
#include "TTree.h"
#include "TFile.h"
#include "TVectorD.h"
#include "TH2D.h"
#include "TAxis.h"
#include <vector>
#include "TLorentzVector.h"
#include </data/t2k/phumhh/lbne/srcs/larvalidatereco/larvalidatereco/content/data_objects/LArAnalysis.h>
#include </data/t2k/phumhh/lbne/srcs/larvalidatereco/larvalidatereco/content/data_objects/LArValidation.h>

using namespace std;


void LArPIDCheckEigen() 
{

// Muon Calculation
TFile *ntuple_muon_PID = new TFile("/data/t2k/phumhh/lbne/ntuple_1gev_muon.root");
if ( !ntuple_muon_PID ) { cout << endl << "Unable to open muon ntuple file" << endl; }
 else{ cout << endl <<  "muon ntuple open" << endl; }
TDirectoryFile *valrecfolder_muon = (TDirectoryFile*)ntuple_muon_PID->Get("valrec");
TTree *ValRecTree_muon = (TTree*)valrecfolder_muon->Get("valrec");

LArPID* pid_muon=0;
LArAnalysis* ana_muon=0;
ValRecTree_muon->SetBranchAddress("PID",&pid_muon);
ValRecTree_muon->SetBranchAddress("Analysis",&ana_muon);
 if ( ValRecTree_muon && valrecfolder_muon ) { cout << "data retrieved" << endl; }
Long64_t nentries_muon = ValRecTree_muon->GetEntries();
 std::cout << std::endl << "read entries: " << nentries_muon << std::endl;
 Double_t av_eig_muon=0;
 cout << "calculating" << endl;
for(int i=0;i<nentries_muon;++i)
  {
    if (i%500==0) {cout << "completed entry: " << i << endl;}
    ValRecTree_muon->GetEntry(i);
    if ( pid_muon->EigenValues.size() && ana_muon->NHits )
      {
	av_eig_muon+=pid_muon->EigenValues[0][0];
      }
	
  }
 ntuple_muon_PID->Close();

 // electron Calculation
 TFile *ntuple_electron_PID = new TFile("/data/t2k/phumhh/lbne/ntuple_1gev_electron_tracks.root");
 if ( !ntuple_electron_PID ) { cout << endl << "Unable to open electron ntuple file" << endl; }
 else{ cout << endl <<  "electron ntuple open" << endl; }
 TDirectoryFile *valrecfolder_electron = (TDirectoryFile*)ntuple_electron_PID->Get("valrec");
 TTree *ValRecTree_electron = (TTree*)valrecfolder_electron->Get("valrec");

 LArPID* pid_electron=0;
 LArAnalysis* ana_electron=0;
 ValRecTree_electron->SetBranchAddress("PID",&pid_electron);
 ValRecTree_electron->SetBranchAddress("Analysis",&ana_electron);
 if ( ValRecTree_electron && valrecfolder_electron ) { cout << "data retrieved" << endl; }
 Long64_t nentries_electron = ValRecTree_electron->GetEntries();
 Double_t av_eig_electron=0;
 cout << "calculating" << endl;
 for(int i=0;i<nentries_electron;++i)
   {
     if (i%500==0) {cout << "completed entry: " << i << endl;}
     ValRecTree_electron->GetEntry(i);
     if ( pid_electron->EigenValues.size() && ana_electron->NHits )
       {
	 av_eig_electron+=pid_electron->EigenValues[0][0];
       }
   }
 ntuple_electron_PID->Close();


 // proton Calculation
 TFile *ntuple_proton_PID = new TFile("/data/t2k/phumhh/lbne/ntuple_1gev_proton.root");
 if ( !ntuple_proton_PID ) { cout << endl << "Unable to open proton ntuple file" << endl; }
 else{ cout << endl <<  "proton ntuple open" << endl; }
 TDirectoryFile *valrecfolder_proton = (TDirectoryFile*)ntuple_proton_PID->Get("valrec");
 TTree *ValRecTree_proton = (TTree*)valrecfolder_proton->Get("valrec");

 LArPID* pid_proton=0;
 LArAnalysis* ana_proton=0;
 ValRecTree_proton->SetBranchAddress("PID",&pid_proton);
 ValRecTree_proton->SetBranchAddress("Analysis",&ana_proton);
 if ( ValRecTree_proton && valrecfolder_proton ) { cout << "data retrieved" << endl; }
 Long64_t nentries_proton = ValRecTree_proton->GetEntries();
 Double_t av_eig_proton=0;
 cout << "calculating" << endl;
 for(int i=0;i<nentries_proton;++i)
   {
     if (i%500==0) {cout << "completed entry: " << i << endl;}
     ValRecTree_proton->GetEntry(i);
     if ( pid_proton->EigenValues.size() && ana_proton->NHits )
       {
	 av_eig_proton+=pid_muon->EigenValues[0][0];
       }
   }
 ntuple_proton_PID->Close();

 // piplus Calculation
 TFile *ntuple_piplus_PID = new TFile("/data/t2k/phumhh/lbne/ntuple_1gev_piplus.root");
 if ( !ntuple_piplus_PID ) { cout << endl << "Unable to open piplus ntuple file" << endl; }
 else{ cout << endl <<  "piplus ntuple open" << endl; }
 TDirectoryFile *valrecfolder_piplus = (TDirectoryFile*)ntuple_piplus_PID->Get("valrec");
 TTree *ValRecTree_piplus = (TTree*)valrecfolder_piplus->Get("valrec");

 LArPID* pid_piplus=0;
 LArAnalysis* ana_piplus=0;
 ValRecTree_piplus->SetBranchAddress("PID",&pid_piplus);
 ValRecTree_piplus->SetBranchAddress("Analysis",&ana_piplus);
 if ( ValRecTree_piplus && valrecfolder_piplus ) { cout << "data retrieved" << endl; }
 Long64_t nentries_piplus = ValRecTree_piplus->GetEntries();
 Double_t av_eigen_piplus=0;
 cout << "calculating" << endl;
 for(int i=0;i<nentries_piplus;++i)
   {
     if (i%500==0) {cout << "completed entry: " << i << endl;}
     ValRecTree_piplus->GetEntry(i);
     if ( pid_piplus->EigenValues.size() && ana_piplus->NHits )
       {
	 av_eigen_piplus+=pid_muon->EigenValues[0][0];
       }
   }
 std::cout << std::endl << "Average Muon Principal Eigenvalue: " << av_eig_muon/nentries_muon << std::endl;
 std::cout << std::endl << "Average Electron Principal Eigenvalue: " << av_eig_electron/nentries_electron << std::endl;
 std::cout << std::endl << "Average Proton Principal Eigenvalue: " << av_eig_proton/nentries_proton << std::endl;
 std::cout << std::endl << "Average Pi Plus Principal Eigenvalue: " << av_eigen_piplus/nentries_piplus << std::endl << std::endl;
 ntuple_piplus_PID->Close();

   
   
}
