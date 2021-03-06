// This file generates a separate plot and ntuple tree (in /pid_files) of electrons vs. protons vs. piplus)    

// Design variable to compare the core of the track shower to its halo
// As LArSOFT uses cm as its distance unit, we can define some deviation from the pca principal axis
// as a threshold boundary between the 'core' and 'halo' of the shower / track object
// This was defined as half the Moliere radius in LAr: 9.61cm
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
#include "TLegend.h"
#include "TAxis.h"
#include <vector>
#include "TLorentzVector.h"
#include </data/t2k/phumhh/lbne/srcs/larvalidatereco/larvalidatereco/content/data_objects/LArAnalysis.h>
#include </data/t2k/phumhh/lbne/srcs/larvalidatereco/larvalidatereco/content/data_objects/LArValidation.h>

using namespace std;
#define MOLIERE_RAD 9.61/2

double HaloNum(std::vector<TVector3> track)
{
  Double_t halo = 0;
  for ( auto i = track.begin(); i!=track.end(); i++)
    {
      if ( TMath::Abs(i->Y()) > MOLIERE_RAD || TMath::Abs(i->Z()) > MOLIERE_RAD  )  {  halo++; }
    }
  return halo;
}

double CoreNum(std::vector<TVector3> track)
{
  Double_t core = 0;
  for ( auto i = track.begin(); i!=track.end(); i++ )
    {
      if ( TMath::Abs(i->Y()) < MOLIERE_RAD || TMath::Abs(i->Z()) < MOLIERE_RAD ) { core++; }
    }
  return core;
}


void LArPIDVarCoreHaloAll() 
{

  TFile* HistFile = new TFile("pid_files/VarCoreHaloAll.root","RECREATE");
  
  TH1D* hist_proton_corehalo = new TH1D("hist_proton_corehalo","Proton",50,0,0.05);
  TH1D* hist_piplus_corehalo = new TH1D("hist_piplus_corehalo", "Piplus",50,0,0.05);
  TH1D* hist_signal_corehalo = new TH1D("hist_signal_corehalo","Signal",50,0,0.05);
  Double_t corehalo_proton, corehalo_piplus, corehalo_signal;
 hist_proton_corehalo->GetXaxis()->SetTitle("Halo / Core Hits");
 hist_proton_corehalo->GetYaxis()->SetTitle("Arbitrary Units");
 hist_piplus_corehalo->GetXaxis()->SetTitle("Halo / Core Hits");
 hist_piplus_corehalo->GetYaxis()->SetTitle("Arbitrary Units");
 hist_signal_corehalo->GetXaxis()->SetTitle("CoreHalo Variable");
 hist_signal_corehalo->GetYaxis()->SetTitle("Arbitrary Units");
 hist_proton_corehalo->SetStats(0);
 hist_signal_corehalo->SetStats(0);
 hist_piplus_corehalo->SetStats(0);
 // Electron Calculation
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
 Double_t* halo_electron = new Double_t;
 Double_t* core_electron = new Double_t;
 
 std::vector<Double_t>* corehalo_hold_electron = new std::vector<Double_t>;

 cout << "calculating" << endl;
 for(int i=0;i<nentries_electron;++i)
   {
     if (i%500==0) {cout << "completed entry: " << i << endl;}
     ValRecTree_electron->GetEntry(i);
     if ( pid_electron->PCAHitsSpacePoints.size() && ana_electron->NHits )
       {
	 for (auto j = pid_electron->PCAHitsSpacePoints.begin(); j!=pid_electron->PCAHitsSpacePoints.end(); j++)
	   {
	     *halo_electron = HaloNum(*j);
	     *core_electron = CoreNum(*j);
	     if ( *halo_electron && *core_electron) 
	       { 
		 hist_signal_corehalo->Fill((*halo_electron) / (*core_electron)); 
		 corehalo_signal = ( (*halo_electron) / (*core_electron) );
		 corehalo_hold_electron->push_back(*halo_electron / *core_electron);
	       }
	   }
       }
   }
 ntuple_electron_PID->Close();
 
 HistFile->cd();
 hist_signal_corehalo->Write();

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

 TFile *ntuple_piplus_PID = new TFile("/data/t2k/phumhh/lbne/ntuple_1gev_piplus.root");
 if ( !ntuple_piplus_PID ) { cout << endl << "Unable to open piplus ntuple file" << endl; }
 else{ cout << endl <<  "piplus ntuple open" << endl; }
 TDirectoryFile *valrecfolder_piplus = (TDirectoryFile*)ntuple_piplus_PID->Get("valrec");
 TTree *ValRecTree_piplus = (TTree*)valrecfolder_piplus->Get("valrec");

 Long64_t nentries_piplus = ValRecTree_piplus->GetEntries();

 Double_t* halo_proton = new Double_t;
 Double_t* core_proton = new Double_t;

 std::vector<Double_t>* corehalo_hold_proton = new std::vector<Double_t>;
 cout << "calculating" << endl;
 for(int i=0;i<nentries_proton;++i)
   {
     if (i%500==0) {cout << "completed entry: " << i << endl;}
     ValRecTree_proton->GetEntry(i);
     if ( pid_proton->PCAHitsSpacePoints.size() && ana_proton->NHits )
       {
	 for (auto j = pid_proton->PCAHitsSpacePoints.begin(); j!=pid_proton->PCAHitsSpacePoints.end(); j++)
	   {
	     *halo_proton = HaloNum(*j);
	     *core_proton = CoreNum(*j);
	     if ( *halo_proton && *core_proton) 
	       { 
		 hist_proton_corehalo->Fill((*halo_proton) / (*core_proton)); 
		 corehalo_hold_proton->push_back( (*halo_proton) / (*core_proton) );
		 
	       }
	   }
       }
   }
 ntuple_proton_PID->Close();

 HistFile->cd();


 // PiPlus Calculation

 LArPID* pid_piplus=0;
 LArAnalysis* ana_piplus=0;
 ValRecTree_piplus->SetBranchAddress("PID",&pid_piplus);
 ValRecTree_piplus->SetBranchAddress("Analysis",&ana_piplus);
 if ( ValRecTree_piplus && valrecfolder_piplus ) { cout << "data retrieved" << endl; }
 Double_t* halo_piplus = new Double_t;
 Double_t* core_piplus = new Double_t;

 std::vector<Double_t>* corehalo_hold_piplus = new std::vector<Double_t>;
 cout << "calculating" << endl;
 for(int i=0;i<nentries_piplus;++i)
   {
     if (i%500==0) {cout << "completed entry: " << i << endl;}
     ValRecTree_piplus->GetEntry(i);
     if ( pid_piplus->PCAHitsSpacePoints.size() && ana_piplus->NHits )
       {
	 for (auto j = pid_piplus->PCAHitsSpacePoints.begin(); j!=pid_piplus->PCAHitsSpacePoints.end(); j++)
	   {
	     *halo_piplus = HaloNum(*j);
	     *core_piplus = CoreNum(*j);
	     if ( *halo_piplus && *core_piplus) 
	       { 
		 hist_piplus_corehalo->Fill((*halo_piplus) / (*core_piplus)); 
		 corehalo_hold_piplus->push_back( (*halo_piplus) / (*core_piplus) );
		 
	       }    
	   }
       }
   }
 ntuple_piplus_PID->Close();
 
 TTree* VarTree = new TTree("CoreHaloTree","CoreHaloTree");

  VarTree->Branch("corehalo_proton", &corehalo_proton, "corehalo_proton/D");
  VarTree->Branch("corehalo_piplus", &corehalo_piplus, "corehalo_piplus/D");
  VarTree->Branch("corehalo_signal", &corehalo_signal, "corehalo_signal/D");
  Int_t nentries_fill;
  if (corehalo_hold_proton->size()> corehalo_hold_electron->size() && corehalo_hold_piplus->size()>corehalo_hold_electron->size()) { nentries_fill = corehalo_hold_electron->size(); }
  else { nentries_fill = corehalo_hold_proton->size(); }

 for (int i = 0 ; i<nentries_fill; i++)
   {
     if (corehalo_hold_proton->at(i) && corehalo_hold_piplus->at(i) && corehalo_hold_electron->at(i) )
       { 
	 corehalo_proton = corehalo_hold_proton->at(i);
	 corehalo_piplus = corehalo_hold_piplus->at(i);
	 corehalo_signal = corehalo_hold_electron->at(i);
	 VarTree->Fill();
       }
   }

 HistFile->cd();
 VarTree->Write();
 hist_proton_corehalo->Write();
 hist_piplus_corehalo->Write();

 hist_proton_corehalo->SetDirectory(0); hist_signal_corehalo->SetDirectory(0); hist_piplus_corehalo->SetDirectory(0);
 hist_proton_corehalo->SetFillColor(kBlue);
 hist_proton_corehalo->SetLineColor(kBlue);
 hist_proton_corehalo->SetFillStyle(3001);
 hist_signal_corehalo->SetFillColor(kRed);
 hist_signal_corehalo->SetLineColor(kRed);
 hist_signal_corehalo->SetFillStyle(3004);
 hist_piplus_corehalo->SetFillColor(kGreen);
 hist_piplus_corehalo->SetFillStyle(3005);
 hist_piplus_corehalo->SetLineColor(kGreen);
 HistFile->Close();
 TCanvas* mycan = new TCanvas("mycan","Core to Halo Variable",1600,100,1200,1200);
 TLegend* myleg = new TLegend(0.75,0.75,0.95,0.95);
 myleg->AddEntry(hist_proton_corehalo);
 myleg->AddEntry(hist_signal_corehalo);
 myleg->AddEntry(hist_piplus_corehalo);
 hist_signal_corehalo->Draw();
 myleg->Draw("SAME");
 hist_piplus_corehalo->Draw("SAME");
 hist_proton_corehalo->Draw("SAME");
 mycan->Update();
}
