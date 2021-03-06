// Initial check of the validity of the TPrincipal PCA output, plotting a 2D histogram of the first two eigenvalues.

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

void LArPIDCheck3() 
{

  TFile* HistFile = new TFile("pid_media/EigenPlot.root","RECREATE");
  TH2D* MuonHist = new TH2D("MuonHist","Muon Sample",10,0,1,10,0,1);
  TH2D* ElectronHist = new TH2D("ElectronHist","Electron Sample",10,0,1,10,0,1);
  MuonHist->GetXaxis()->SetTitle("X");
  MuonHist->GetYaxis()->SetTitle("Y");
  ElectronHist->GetXaxis()->SetTitle("X");
  ElectronHist->GetYaxis()->SetTitle("Y");

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

 Double_t* Eig_Met_Con_co_muon = new Double_t[3];

 cout << "calculating" << endl;
for(int i=0;i<nentries_muon;++i)
  {
    if (i%500==0) {cout << "completed entry: " << i << endl;}
    ValRecTree_muon->GetEntry(i);
    if ( pid_muon->EigenValues.size() && ana_muon->NHits )
      {
	for (auto j = pid_muon->EigenValues.begin(); j!=pid_muon->EigenValues.end(); j++)
	  {
	    
	    Eig_Met_Con_co_muon = j->GetMatrixArray();
	    //Eig_Met_Con_muon->push_back(*Eig_Met_Con_co_muon);
	    //if (Eig_Met_Con_muon->at(vector_count) > muonhist_max) { muonhist_max = Eig_Met_Con_muon->at(vector_count); }
	    MuonHist->Fill(Eig_Met_Con_co_muon[0], Eig_Met_Con_co_muon[1]);
	  
	  }
	  }
  }
 ntuple_muon_PID->Close();


 HistFile->cd();
 MuonHist->Write();

   
   
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

 Double_t* Eig_Met_Con_co_electron = new Double_t[3];

 cout << "calculating" << endl;
 for(int i=0;i<nentries_electron;++i)
   {
     if (i%500==0) {cout << "completed entry: " << i << endl;}
     ValRecTree_electron->GetEntry(i);
     if ( pid_electron->EigenValues.size() && ana_electron->NHits )
       {
	 for (auto j = pid_electron->EigenValues.begin(); j!=pid_electron->EigenValues.end(); j++)
	   {

	     Eig_Met_Con_co_electron = j->GetMatrixArray();
	     //Eig_Met_Con_electron->push_back(*Eig_Met_Con_co_electron);
	     //if (Eig_Met_Con_electron->at(vector_count) > electronhist_max) { electronhist_max = Eig_Met_Con_electron->at(vector_count); }
	     ElectronHist->Fill(Eig_Met_Con_co_electron[0], Eig_Met_Con_co_electron[1]);
	 
	   }
       }
   }
 ntuple_electron_PID->Close();


 HistFile->cd();
 ElectronHist->Write();


 MuonHist->SetDirectory(0); ElectronHist->SetDirectory(0);
 HistFile->Close();
 TCanvas* mycan = new TCanvas("mycan","Eigenvalues",1600,100,1200,1200);
 mycan->Divide(1,2);
 mycan->cd(1);
 MuonHist->Draw("COLZ");
 mycan->cd(2);
 ElectronHist->Draw("COLZ");
 mycan->cd(); mycan->Update();

}
