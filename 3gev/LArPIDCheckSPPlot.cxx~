// Plotting PCA transformed spacepoints to check shower shape
// Shows electrons more diffuse than muons, indicating the PCA does retain shower shape

#include <string>
#include "TSystem.h"
#include "TCanvas.h"
#include "TObject.h"
#include "../data_objects/LArPID.h"
#include "iostream"
#include "TTree.h"
#include "TFile.h"
#include "TVectorD.h"
#include "TH3D.h"
#include "TAxis.h"
#include <vector>
#include "TLorentzVector.h"
#include </data/t2k/phumhh/lbne/srcs/larvalidatereco/larvalidatereco/content/data_objects/LArAnalysis.h>
#include </data/t2k/phumhh/lbne/srcs/larvalidatereco/larvalidatereco/content/data_objects/LArValidation.h>

using namespace std;

void LArPIDCheckSPPlot() 
{

  TFile* HistFile = new TFile("pid_files/SpacePointsPlot.root","RECREATE");

  TH3D* PCAHist = new TH3D("PCAHist","PCA Hits",10,-5,5,10,-3,3,10,-3,3);
  TH3D* SPHist = new TH3D("SPHist","Real Hits",100,50,150,100,-100,100,100,0,200);
 PCAHist->GetXaxis()->SetTitle("X");
 PCAHist->GetYaxis()->SetTitle("Y");
 PCAHist->GetZaxis()->SetTitle("Z");
 SPHist->GetXaxis()->SetTitle("X");
 SPHist->GetYaxis()->SetTitle("Y");
 SPHist->GetZaxis()->SetTitle("Z");
// Fill PCA Plot
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

 Double_t* spc_con_co_muon = new Double_t[3];

 cout << "calculating pca spcpnts" << endl;
for(int i=0;i<nentries_muon;++i)
  {
    if (i%500==0) {cout << "completed entry: " << i << endl;}
    ValRecTree_muon->GetEntry(i);
    if ( pid_muon->PCAHitsSpacePoints.size() && ana_muon->NHits )
      {
	for (auto j = pid_muon->PCAHitsSpacePoints.begin(); j!=pid_muon->PCAHitsSpacePoints.end(); j++)
	  {
	    for (auto k = j->begin(); k!=j->end(); k++ )
	      {
		
		
		k->GetXYZ(spc_con_co_muon);
		
		PCAHist->Fill(spc_con_co_muon[2], spc_con_co_muon[1], spc_con_co_muon[0]);
	      }
	  }
      }
  }

     Long64_t nentries = ValRecTree_muon->GetBranch("Hit3Pos")->GetEntries();
     Double_t * spc_con_co_muon_real = new Double_t[3];
     // Fill vector of original hit spacepoints
     cout << "calculating real spcpnts" << endl;
     for(int i=0;i<nentries;++i)
       {
	 ValRecTree_muon->GetEntry(i);
	 if ( ana_muon->Hit3Pos.size() )
	   {
	     if ( i%500==0) { std::cout << "completed entry: " << i << endl;}
	     for (auto j = ana_muon->Hit3Pos.begin(); j != ana_muon->Hit3Pos.end(); j++ )
	       {
		 j->GetXYZ(spc_con_co_muon_real);
		 
		 SPHist->Fill(spc_con_co_muon_real[0],spc_con_co_muon_real[1],spc_con_co_muon_real[2]);
	       }

	   }
       }

 ntuple_muon_PID->Close();


 HistFile->cd();
 PCAHist->Write();
 SPHist->Write();
 PCAHist->SetDirectory(0); SPHist->SetDirectory(0);
 HistFile->Close();


}
