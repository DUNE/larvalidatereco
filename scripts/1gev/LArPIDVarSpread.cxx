// Design variable using PCA eigenvalues to give the width/length ratio of the shower     
// This file generates a plot and ntuple tree (in /pid_files) of signal (electrons) vs. background (protons & piplus) 
#include <string>
#include "TSystem.h"
#include "TCanvas.h"
#include "TObject.h"
#include "../data_objects/LArPID.h"
#include "iostream"
#include "TTree.h"
#include "TFile.h"
#include "TVectorD.h"
#include "TH1D.h"
#include "TAxis.h"
#include "TLegend.h"
#include <vector>
#include "TLorentzVector.h"
#include </data/t2k/phumhh/lbne/srcs/larvalidatereco/larvalidatereco/content/data_objects/LArAnalysis.h>
#include </data/t2k/phumhh/lbne/srcs/larvalidatereco/larvalidatereco/content/data_objects/LArValidation.h>

using namespace std;


void LArPIDVarSpread() 
{

  TFile* HistFile = new TFile("pid_files/VarSpread.root", "RECREATE");
  TH1D* hist_bg_spread = new TH1D("hist_bg_spread","Background Signal",50,0,1.4);
  TH1D* hist_signal_spread = new TH1D("hist_signal_spread","Electron Signal",50,0,1.4);
 hist_bg_spread->GetXaxis()->SetTitle("Width / Length");
 hist_bg_spread->GetYaxis()->SetTitle("");
 hist_bg_spread->SetStats(0);
 Double_t spread_bg = 0; 
 Double_t spread_signal = 0;

 std::vector<Double_t>* spread_hold_bg = new std::vector<Double_t>;

 // Proton Calculation
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
 Double_t* eigen_hold_proton = new Double_t[3];

 cout << "calculating" << endl;
 for(int i=0;i<nentries_proton;++i)
   {
     if (i%500==0) {cout << "completed entry: " << i << endl;}
     ValRecTree_proton->GetEntry(i);
     if ( pid_proton->EigenValues.size() && ana_proton->NHits )
       {
         for (auto j = pid_proton->EigenValues.begin(); j!=pid_proton->EigenValues.end(); j++)
           {
             eigen_hold_proton = j->GetMatrixArray();
             spread_hold_bg->push_back(std::sqrt(TMath::Power(eigen_hold_proton[1],2) + TMath::Power(eigen_hold_proton[2],2))/eigen_hold_proton[0]);
	     
             hist_bg_spread->Fill(std::sqrt(TMath::Power(eigen_hold_proton[1],2) + TMath::Power(eigen_hold_proton[2],2))/eigen_hold_proton[0]);
           }
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
 Double_t* eigen_hold_piplus = new Double_t[3];

 cout << "calculating" << endl;
 for(int i=0;i<nentries_piplus;++i)
   {
     if (i%500==0) {cout << "completed entry: " << i << endl;}
     ValRecTree_piplus->GetEntry(i);
     if ( pid_piplus->EigenValues.size() && ana_piplus->NHits )
       {
         for (auto j = pid_piplus->EigenValues.begin(); j!=pid_piplus->EigenValues.end(); j++)
           {
             eigen_hold_piplus = j->GetMatrixArray();
             spread_hold_bg->push_back(std::sqrt(TMath::Power(eigen_hold_piplus[1],2) + TMath::Power(eigen_hold_piplus[2],2))/eigen_hold_piplus[0]);
	     
             hist_bg_spread->Fill(std::sqrt(TMath::Power(eigen_hold_piplus[1],2) + TMath::Power(eigen_hold_piplus[2],2))/eigen_hold_piplus[0]);
           }
       }
   }
 ntuple_piplus_PID->Close();




 hist_signal_spread->GetXaxis()->SetTitle("Spread Variable");
 hist_signal_spread->GetYaxis()->SetTitle("Arbitrary Unit");
 hist_signal_spread->SetStats(0);

 std::vector<Double_t>* spread_hold_electron = new std::vector<Double_t>;

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
 Double_t* eigen_hold_electron = new Double_t[3];

 cout << "calculating" << endl;
 for(int i=0;i<nentries_electron;++i)
   {
     if (i%500==0) {cout << "completed entry: " << i << endl;}
     ValRecTree_electron->GetEntry(i);
     if ( pid_electron->EigenValues.size() && ana_electron->NHits )
       {
	 for (auto j = pid_electron->EigenValues.begin(); j!=pid_electron->EigenValues.end(); j++)
	   {
	     eigen_hold_electron = j->GetMatrixArray();
	     spread_hold_electron->push_back(std::sqrt(TMath::Power(eigen_hold_electron[1],2) + TMath::Power(eigen_hold_electron[2],2))/eigen_hold_electron[0]);
	     
	     hist_signal_spread->Fill(std::sqrt(TMath::Power(eigen_hold_electron[1],2) + TMath::Power(eigen_hold_electron[2],2))/eigen_hold_electron[0]);
	   }
       }
   }  
 ntuple_electron_PID->Close();

  TTree* VarTree = new TTree("SpreadTree","SpreadTree");
  VarTree->Branch("spread_bg", &spread_bg, "spread_bg/D");
  VarTree->Branch("spread_signal", &spread_signal, "spread_signal/D");
  Int_t nentries_fill;
  if (spread_hold_bg->size() > spread_hold_electron->size() ) { nentries_fill = spread_hold_electron->size(); }
  else { nentries_fill = spread_hold_bg->size(); }
  for (int i=0; i<nentries_fill; i++)
    {
      if (spread_hold_bg->at(i) && spread_hold_electron->at(i) )
	{
	  spread_bg = spread_hold_bg->at(i);
	  spread_signal = spread_hold_electron->at(i);
	  VarTree->Fill();
	}
    }
  HistFile->cd();
  VarTree->Write();
  


 Double_t norm = hist_bg_spread->GetEntries() + hist_signal_spread->GetEntries();
 hist_bg_spread->Scale(1/norm);
 hist_signal_spread->Scale(1/norm);
 hist_bg_spread->Write();
 hist_bg_spread->SetDirectory(0);
 hist_signal_spread->SetMaximum(0.25);
 hist_signal_spread->Write();
 hist_signal_spread->SetDirectory(0);
 
 hist_signal_spread->SetFillColor(kRed);
 hist_bg_spread->SetFillColor(kBlue);
 

 TCanvas *mycan = new TCanvas("mycan","Spread Variable", 1600,100,1200,1200);
 TLegend* myleg = new TLegend(0.75,0.75,0.95,0.95);
 myleg->AddEntry(hist_signal_spread);
 myleg->AddEntry(hist_bg_spread);
 
 hist_signal_spread->Draw();
 hist_bg_spread->Draw("SAME"); 
 myleg->Draw();
 HistFile->Close();
}
