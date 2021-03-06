// Shower-shape PID variable to measure 'conicalness' of the tracks. Takes the sd of hits from the principal axis in the first 20% of the track and the last 20%, and takes their ratio.                                                        
// This file generates a separate plot and ntuple tree (in /pid_files) of electrons vs. protons vs. piplus) 

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

// Fill track objects from ntuple
std::vector<std::vector<TVector3>> FillCon(const char* filename)
{
TFile *ntuple_PID = new TFile(filename);
 ntuple_PID->cd();
if ( !ntuple_PID ) { cout << endl << "Unable to open muon ntuple file" << endl; }
 else{ cout << endl <<  "ntuple open: " << filename << endl; }
TDirectoryFile *valrecfolder = (TDirectoryFile*)ntuple_PID->Get("valrec");
TTree *ValRecTree = (TTree*)valrecfolder->Get("valrec");

LArPID* pid=0;
LArAnalysis* ana=0;
ValRecTree->SetBranchAddress("PID",&pid);
ValRecTree->SetBranchAddress("Analysis",&ana);
 if ( ValRecTree && valrecfolder ) { cout << "data retrieved" << endl; }
 Long64_t nentries = ValRecTree->GetBranch("PCAHitsSpacePoints")->GetEntries();
 cout << endl << "Entries read: " << nentries << endl;
 std::vector<std::vector<TVector3>>* track_con = new std::vector<std::vector<TVector3>>;
 std::vector<TVector3>* hit_con = new std::vector<TVector3>;

 // Fill vector of pca transformed tracks
 cout << "calculating" << endl;
for(int i=0;i<nentries;++i)
  {
    ValRecTree->GetBranch("PCAHitsSpacePoints")->GetEntry(i);
    if ( pid->PCAHitsSpacePoints.size() )	  
      {
	if ( i%500==0) { std::cout << "Event: " << i << endl;}
		for (auto j = pid->PCAHitsSpacePoints.begin(); j != pid->PCAHitsSpacePoints.end(); j++ )
		  { 
		    for (auto k = j->begin(); k!=j->end(); k++) 
		      {
			hit_con->push_back(*k); 
		      }
		    track_con->push_back(*hit_con);
		    hit_con->clear();
		  }
	  }
  }
 ntuple_PID->Close();
 return *track_con;
}

// Calculate the maximum track coordiante along the principal axis
double FindTrackMax(std::vector<TVector3> track)

{
  Double_t max = track.begin()->X();
  for ( auto hit_pntr = track.begin() ; hit_pntr!=track.end(); hit_pntr++ )
    {
     
      if ( hit_pntr->X() > max ) { max = hit_pntr->X(); }
    }
  return max;
}
// Calculate the minimum track coordinate along the principal axis
double FindTrackMin(std::vector<TVector3> track)

{
  Double_t min = track.begin()->X();
  for ( auto hit_pntr = track.begin() ; hit_pntr!=track.end(); hit_pntr++ )
    {

      if ( hit_pntr->X() < min ) { min = hit_pntr->X(); }
    }
  return min;
}
// Calculate the sd dev of hits from the principal axis along the first 20% of the track
double CalcSdStart(Double_t trackmax, Double_t trackmin,std::vector<TVector3> track)
{
  Double_t mag[track.size()];
  Double_t tot=0;
  Int_t hitnum=0;
  Double_t mean;
  Double_t sd = 0;
  Double_t* tracklen = new Double_t;
  if ( trackmin<0 && trackmax>0 ) {*tracklen = trackmax+TMath::Abs(trackmin);} 
  else if ( trackmin>0 && trackmax>0 ) {*tracklen = trackmax-trackmin;}
  else if (trackmin<0 && trackmax<0 ) { *tracklen = TMath::Abs(TMath::Abs(trackmax)-TMath::Abs(trackmin));}
  Double_t startmax = trackmin + (*tracklen)/5;
  for (auto i = track.begin(); i!=track.end(); i++) 
    {
      if ( i->X() < startmax )
	{
	  mag[hitnum]= std::sqrt( (i->Y()*i->Y()) + (i->Z())*(i->Z()) ); 
	  hitnum++;
	}
    }
  for (unsigned int j = 0; j!=track.size(); j++) {tot+=mag[j];}
  mean = tot/track.size();
  for (unsigned int k = 0; k!=track.size(); k++) {  sd += (mag[k]-mean)*(mag[k]-mean); }
  return std::sqrt(sd/(track.size()-1));
}
double CalcSdEnd(Double_t trackmax, Double_t trackmin,std::vector<TVector3> track)
{
  Double_t mag[track.size()];
  Double_t tot=0;
  Int_t hitnum=0;
  Double_t mean;
  Double_t sd = 0;
  Double_t* tracklen = new Double_t; 
  if ( trackmin<0 && trackmax>0 ) {*tracklen = trackmax + TMath::Abs(trackmin);}
  else if (trackmin>0 && trackmax>0) { *tracklen = trackmax - trackmin;}
  else if (trackmin<0 && trackmax < 0) { *tracklen = TMath::Abs(TMath::Abs(trackmax)-TMath::Abs(trackmin));}
Double_t endmin = trackmax - (*tracklen)/5;
  for (auto i = track.begin(); i!=track.end(); i++)
    {
      if ( i->X() > endmin )
        {
          mag[hitnum]= std::sqrt( (i->Y()*i->Y()) + (i->Z())*(i->Z())  );
          hitnum++;
        }
    }
  for (unsigned int j = 0; j!=track.size(); j++) {tot+=mag[j];};
  mean = tot/track.size();
  for (unsigned int k = 0; k!=track.size(); k++) { sd+= (mag[k]-mean)*(mag[k]-mean); }
  return std::sqrt(sd/(track.size()-1));
}


void LArPIDVarConAll() 
{
  Int_t nbins = 40;
  TFile* histfile = new TFile("pid_files/VarConAll.root","RECREATE");
  TH1D *hist_proton_con = new TH1D("hist_proton_con","Proton Signal - Conicalness", nbins,0.5,1);
  TH1D* hist_piplus_con = new TH1D("hist_piplus_con","Piplus Signal - Conicalness", nbins,0.5,1);
  // ITERATE OVER TRACKS AND SORT HITS INTO ASCENDING ORDER OF PRINCIPAL AXES
  std::vector<Double_t>* con_hold_proton = new std::vector<Double_t>;
  std::vector<Double_t>* con_hold_piplus = new std::vector<Double_t>;
  // Proton Calculation
  Double_t trackmax_proton, trackmin_proton,  sdstart_proton, sdend_proton;
  Double_t con_proton, con_piplus, con_signal;

  // ITERATE OVER TRACKS AND SORT HITS INTO ASCENDING ORDER OF PRINCIPAL AXES
  std::vector<std::vector<TVector3>> tracks_proton = FillCon("/data/t2k/phumhh/lbne/ntuple_1gev_proton.root");
  for ( auto track_count = tracks_proton.begin(); track_count!=tracks_proton.end(); track_count++ )
    {
      trackmax_proton = FindTrackMax(*track_count);
      trackmin_proton = FindTrackMin(*track_count);
      sdstart_proton = CalcSdStart(trackmax_proton, trackmin_proton, *track_count);
      sdend_proton = CalcSdEnd(trackmax_proton, trackmin_proton, *track_count);
      if ( sdend_proton >= sdstart_proton) { 
	hist_proton_con->Fill(sdstart_proton/sdend_proton); 
	con_hold_proton->push_back(sdstart_proton/sdend_proton);
      }
    }


  // Piplus Calculation
  Double_t trackmax_piplus, trackmin_piplus, sdstart_piplus, sdend_piplus;
  // ITERATE OVER TRACKS AND SORT HITS INTO ASCENDING ORDER OF PRINCIPAL AXES
  std::vector<std::vector<TVector3>> tracks_piplus = FillCon("/data/t2k/phumhh/lbne/ntuple_1gev_piplus.root");
  for ( auto track_count = tracks_piplus.begin(); track_count!=tracks_piplus.end(); track_count++ )
    {
      trackmax_piplus = FindTrackMax(*track_count);
      trackmin_piplus = FindTrackMin(*track_count);
      sdstart_piplus = CalcSdStart(trackmax_piplus, trackmin_piplus, *track_count);
      sdend_piplus = CalcSdEnd(trackmax_piplus, trackmin_piplus, *track_count);
      if ( sdend_piplus >= sdstart_piplus ) 
	{
	  hist_piplus_con->Fill(sdstart_piplus/sdend_piplus);
	  con_hold_piplus->push_back(sdstart_piplus/sdend_piplus);
	}
    }
  
  std::vector<Double_t>* con_hold_electron = new std::vector<Double_t>;
 
  // Electron Calculation
  Double_t trackmax_electron, trackmin_electron, sdstart_electron, sdend_electron;
  TH1D *hist_signal_con = new TH1D("hist_signal_con","Electron Signal - Conicalness", nbins,0.5,1);
  // ITERATE OVER TRACKS AND SORT HITS INTO ASCENDING ORDER OF PRINCIPAL AXES
  std::vector<std::vector<TVector3>> tracks_electron = FillCon("/data/t2k/phumhh/lbne/ntuple_1gev_electron_tracks.root");
  for ( auto track_count = tracks_electron.begin(); track_count!=tracks_electron.end(); track_count++ )
    {
      trackmax_electron = FindTrackMax(*track_count);
      trackmin_electron = FindTrackMin(*track_count);
      sdstart_electron = CalcSdStart(trackmax_electron, trackmin_electron, *track_count);
      sdend_electron = CalcSdEnd(trackmax_electron, trackmin_electron,*track_count);
      if (sdend_electron >= sdstart_electron ) 
	{ 
	  hist_signal_con->Fill(sdstart_electron/sdend_electron); 
	  con_hold_electron->push_back(sdstart_electron/sdend_electron);
	}
    }

  TTree* VarTree = new TTree("ConTree", "ConTree");
  VarTree->Branch("con_proton", &con_proton, "con_proton/D");
  VarTree->Branch("con_signal", &con_signal, "con_signal/D");
  VarTree->Branch("con_piplus",&con_piplus,"con_piplus/D");
  Int_t nentries_fill;
  if (con_hold_proton->size() > con_hold_electron->size() && con_hold_piplus->size()>con_hold_electron->size() ) { nentries_fill = con_hold_electron->size();}
  else { nentries_fill = con_hold_proton->size(); }
  
  for (int i = 0; i<nentries_fill; i++ ) 
    {
      con_proton = con_hold_proton->at(i);
      con_piplus = con_hold_piplus->at(i);
      con_signal = con_hold_electron->at(i);
      VarTree->Fill();
    }
  histfile->cd();
  VarTree->Write();

  hist_proton_con->SetDirectory(0);
  hist_piplus_con->SetDirectory(0);
  hist_signal_con->SetDirectory(0);
  TCanvas* mycan = new TCanvas("mycan","Track 'Conicalness' Variable",1600,100,1200,1200);
  //mycan->Divide(2,1);
  TLegend* myleg = new TLegend(0.5,0.65,0.75,0.85);
  myleg->AddEntry(hist_signal_con);
  myleg->AddEntry(hist_proton_con);
  myleg->AddEntry(hist_piplus_con);

  //  mycan->cd(2);
  hist_signal_con->SetStats(0);
  hist_signal_con->GetXaxis()->CenterTitle();
  hist_signal_con->GetXaxis()->SetTitle("Conicalness Variable"); 
  hist_signal_con->GetYaxis()->SetTitle("Arbitrary Units");
  hist_signal_con->GetYaxis()->CenterTitle();
  hist_signal_con->SetFillColor(kBlue);
  hist_signal_con->SetFillStyle(3001);
  hist_signal_con->SetLineColor(kBlue);
  hist_signal_con->Write();

  //mycan->cd(1);
  hist_proton_con->SetStats(0);
  hist_proton_con->GetXaxis()->SetTitle("Conicalness Variable");
  hist_proton_con->GetXaxis()->CenterTitle();
  hist_proton_con->GetYaxis()->SetTitle("Arbitrary Units");
  hist_proton_con->GetYaxis()->CenterTitle();
  hist_proton_con->SetFillColor(kRed);
  hist_proton_con->SetFillStyle(3004);
  hist_proton_con->SetLineColor(kRed);

  hist_piplus_con->SetStats(0);
  hist_piplus_con->GetXaxis()->SetTitle("Conicalness Variable");
  hist_piplus_con->GetXaxis()->CenterTitle();
  hist_piplus_con->GetYaxis()->SetTitle("Arbitrary Units");
  hist_piplus_con->GetYaxis()->CenterTitle();
  hist_piplus_con->SetFillColor(kGreen);
  hist_piplus_con->SetFillStyle(3005);
  hist_piplus_con->SetLineColor(kGreen);

  hist_signal_con->Draw();
  hist_piplus_con->Draw("SAME");  
  hist_proton_con->Draw("SAME");
  myleg->Draw("SAME");
  hist_proton_con->Write();
  Double_t max = hist_signal_con->GetMaximum();
  hist_signal_con->SetMaximum(max);
  
  mycan->Update();

  histfile->Close();

      
}

