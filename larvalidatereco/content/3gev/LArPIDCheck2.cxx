// An initial check for the validity of the TPrincipal PCA output, using the square sum of the first two eigenvalues as a metric, plotted against momentum. 


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

Double_t EigenMetMag( Double_t Eigen_1, Double_t Eigen_2 ) 
{ 
  Double_t eigenmag = sqrt((Eigen_1*Eigen_1) + (Eigen_2*Eigen_2));
  if ( (eigenmag < 1) && (eigenmag > 0) ) 
    {return eigenmag; }
  else {return 0; }
}


void LArPIDCheck2() 
{

  TFile* HistFile = new TFile("pid_files/MomVEigenHist2.root","RECREATE");
 TH2D* MuonHist = new TH2D("MuonHist","Muon Sample",10,0,2,10,0,1);
 TH2D* ElectronHist = new TH2D("ElectronHist","Electron Sample",10,0,2,10,0,1);
 TH2D* ProtonHist = new TH2D("ProtonHist","Proton Sample",10,0,2,10,0,1);
 TH2D* PiPlusHist = new TH2D("PiPlusHist","PiPlus Sample",10,0,2,10,0,1);

 MuonHist->GetXaxis()->SetTitle("Momentum");
 MuonHist->GetYaxis()->SetTitle("EigenMetric");
 ElectronHist->GetXaxis()->SetTitle("Momentum");
 ElectronHist->GetYaxis()->SetTitle("EigenMetric");
 ProtonHist->GetXaxis()->SetTitle("Momentum");
 ProtonHist->GetYaxis()->SetTitle("EigenMetric");
 PiPlusHist->GetXaxis()->SetTitle("Momentum");
 PiPlusHist->GetYaxis()->SetTitle("EigenMetric");

// Muon Calculation
TFile *ntuple_muon_PID = new TFile("/data/t2k/phumhh/lbne/ntuple_3gev_muon.root");
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
Double_t* Eig_Met_Con_muon = new Double_t[nentries_muon];
 Double_t* Mom_Con_muon = new Double_t[nentries_muon];

 cout << "calculating" << endl;
for(int i=0;i<nentries_muon;++i)
  {
    if (i%500==0) {cout << "completed entry: " << i << endl;}
    ValRecTree_muon->GetEntry(i);
    if ( pid_muon->EigenValues.size() && ana_muon->NHits )
      {
	Eig_Met_Con_muon[i] = EigenMetMag(pid_muon->EigenValues[0][0],pid_muon->EigenValues[0][1]);
	Mom_Con_muon[i] = ana_muon->TrajStart4MomMC[0].Mag();
	
	MuonHist->Fill(Mom_Con_muon[i],Eig_Met_Con_muon[i]);

	  }
  }
 ntuple_muon_PID->Close();

// electron Calculation
 TFile *ntuple_electron_PID = new TFile("/data/t2k/phumhh/lbne/ntuple_3gev_electron_tracks.root");
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
 Double_t* Eig_Met_Con_electron = new Double_t[nentries_electron];
 Double_t* Mom_Con_electron = new Double_t[nentries_electron];

 cout << "calculating" << endl;
 for(int j=0;j<nentries_electron;++j)
   {
     if (j%500==0) {cout << "completed entry: " << j << endl; }
     ValRecTree_electron->GetEntry(j);
     if ( pid_electron->EigenValues.size() && ana_electron->NHits )
       {
	 Eig_Met_Con_electron[j] = EigenMetMag(pid_electron->EigenValues[0][0],pid_electron->EigenValues[0][1]);
	 Mom_Con_electron[j] = ana_electron->TrajStart4MomMC[0].Mag();
	     ElectronHist->Fill(Mom_Con_electron[j],Eig_Met_Con_electron[j]);

       }

   }

 ntuple_electron_PID->Close();


 // proton Calculation
 TFile *ntuple_proton_PID = new TFile("/data/t2k/phumhh/lbne/ntuple_3gev_proton.root");
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
 Double_t* Eig_Met_Con_proton = new Double_t[nentries_proton];
 Double_t* Mom_Con_proton = new Double_t[nentries_proton];

 cout << "calculating" << endl;
 for(int i=0;i<nentries_proton;++i)
   {
     if (i%500==0) {cout << "completed entry: " << i << endl;}
     ValRecTree_proton->GetEntry(i);
     if ( pid_proton->EigenValues.size() && ana_proton->NHits )
       {
	 Eig_Met_Con_proton[i] = EigenMetMag(pid_proton->EigenValues[0][0],pid_proton->EigenValues[0][1]);
	 Mom_Con_proton[i] = ana_proton->TrajStart4MomMC[0].Mag();

	 ProtonHist->Fill(Mom_Con_proton[i],Eig_Met_Con_proton[i]);

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
 Double_t* Eig_Met_Con_piplus = new Double_t[nentries_piplus];
 Double_t* Mom_Con_piplus = new Double_t[nentries_piplus];

 cout << "calculating" << endl;
 for(int i=0;i<nentries_piplus;++i)
   {
     if (i%500==0) {cout << "completed entry: " << i << endl;}
     ValRecTree_piplus->GetEntry(i);
     if ( pid_piplus->EigenValues.size() && ana_piplus->NHits )
       {
	 Eig_Met_Con_piplus[i] = EigenMetMag(pid_piplus->EigenValues[0][0],pid_piplus->EigenValues[0][1]);
	 Mom_Con_piplus[i] = ana_piplus->TrajStart4MomMC[0].Mag();

	 PiPlusHist->Fill(Mom_Con_piplus[i],Eig_Met_Con_piplus[i]);

       }
   }
 ntuple_piplus_PID->Close();

 HistFile->cd();
 MuonHist->Write();
 ElectronHist->Write();
 ProtonHist->Write();
 PiPlusHist->Write();

 HistFile->Close();
   
   
}
