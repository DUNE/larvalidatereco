// PID variable to compute the gradient of the SD of distance of hits from the principal axis as a function of the longitudinal position along the track.                                                      
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
#include "TMatrixD.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TAxis.h"
#include <vector>
#include "TLorentzVector.h"
#include </data/t2k/phumhh/lbne/srcs/larvalidatereco/larvalidatereco/content/data_objects/LArAnalysis.h>
#include </data/t2k/phumhh/lbne/srcs/larvalidatereco/larvalidatereco/content/data_objects/LArValidation.h>

using namespace std;

#define SLICE_DENS 10;

 // Fill track objects from ntuple
std::vector<std::vector<TVector3>> FillTracksGrad(const char* filename)
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
		  hit_con->push_back(*k);
		}
	    }

	}
      track_con->push_back(*hit_con);
      hit_con->clear();
	   
    }

  ntuple_PID->Close();
  return *track_con;
}


TVector3 FindTrackEndGrad(std::vector<TVector3> track)

{
  TVector3* max = new TVector3();
  for ( auto hit_pntr = track.begin() ; hit_pntr!=track.end(); hit_pntr++ )
    {

      if ( hit_pntr->X() > max->X() ) { *max = *hit_pntr; }
    }

  return *max;
}
// Calculate the minimum track coordinate along the principal axis
TVector3 FindTrackStartGrad(std::vector<TVector3> track)

{
  TVector3 *min = new TVector3();
  for ( auto hit_pntr = track.begin() ; hit_pntr!=track.end(); hit_pntr++ )
    {

      if ( hit_pntr->X() < min->X() ) { *min = *hit_pntr; }
    }

  return *min;
}



void LArPIDVarGrad()

{

  TFile* GradFile = new TFile("pid_files/VarGrad.root","RECREATE");
  TH1D* hist_bg_grad = new TH1D("hist_bg_grad","Background",20,0,20);
  TH1D* hist_signal_grad = new TH1D("hist_signal_grad","Signal",20,0,20);

  Double_t grad_bg, grad_signal;




  // Declare variables for first and last hit positions per track
  TVector3 *trackstart = new TVector3();
  TVector3 *trackend = new TVector3();
  // Declare vectors to hold the enpoint coordiantes
  Double_t sd_grad_proton;
  Double_t sd_grad_piplus;
  Double_t sd_grad_electron;
  //  GradTree->Branch("sd_grad_proton",&sd_grad_proton,"sd_grad_proton/D");
  //GradTree->Branch("sd_grad_piplus",&sd_grad_piplus,"sd_grad_piplus/D");
  //GradTree->Branch("sd_grad_electron",&sd_grad_electron,"sd_grad_electron/D");
  Int_t slice_num;
  Double_t range, slice_width, slice_lb, slice_up, slice_sd, slice_mean;
    Double_t hit_tot = 0; Double_t slide_sq_tot = 0; Double_t slice_sd_prev = 0;
  std::vector<Double_t>* hit_rad = new std::vector<Double_t>;

  // Fill vectors 
  std::vector<std::vector<TVector3>> trackcon_proton = FillTracksGrad("/data/t2k/phumhh/lbne/ntuple_1gev_proton.root");
  std::vector<std::vector<TVector3>> trackcon_piplus = FillTracksGrad("/data/t2k/phumhh/lbne/ntuple_1gev_piplus.root");
  std::vector<std::vector<TVector3>> trackcon_electron = FillTracksGrad("/data/t2k/phumhh/lbne/ntuple_1gev_electron_tracks.root");


  std::vector<Double_t>* grad_hold_bg = new std::vector<Double_t>;

  // iterate over tracks for pca values

  // proton calculation
  for (auto i = trackcon_proton.begin(); i!=trackcon_proton.end(); i++ )
    {	  
	  // Read out track vector in terms of pca and the position of its endpoints
	  *trackstart = FindTrackStartGrad(*i);
	  *trackend = FindTrackEndGrad(*i);
	  range = TMath::Abs(trackend->X()-trackstart->X());
	  slice_num= (Int_t)range*SLICE_DENS;
	  slice_width = range/slice_num;
  for (Double_t slice=0; slice!= slice_num; slice++)
    {
      slice_lb = trackstart->X() + (slice*slice_width);
      slice_up = trackstart->X() + ((slice+1)*slice_width);
      // ITERATE OVER HITS
      for (auto j = i->begin(); j!=i->end(); j++)
	{
	  if ( j->X() < slice_up && j->X() > slice_lb ) 
	    {
	      hit_rad->push_back(std::sqrt(TMath::Power(j->Y(),2)+TMath::Power(j->Z(),2)));
	    }

	}
      for (unsigned int k = 0; k!=hit_rad->size(); k++)
	{
	  hit_tot+=hit_rad->at(k);
	}
      slice_mean = hit_tot/hit_rad->size();
	for (unsigned int l = 0; l!=hit_rad->size(); l++)
	  {
	    slide_sq_tot+=(TMath::Power((hit_rad->at(l)-slice_mean),2));
	  }
  slice_sd = std::sqrt(slide_sq_tot/hit_rad->size());
  sd_grad_proton = (slice_sd - slice_sd_prev)/slice_width;

  // ADD SD TO NTUPLE HERE
  hist_bg_grad->Fill(sd_grad_proton);
  grad_hold_bg->push_back(sd_grad_proton);
  slice_sd_prev = slice_sd;
  hit_rad->clear();
  hit_tot = 0; slice_sd_prev = 0; slide_sq_tot = 0;
    }
 
    }

  // piplus calculation
  for (auto i = trackcon_piplus.begin(); i!=trackcon_piplus.end(); i++ )
    {
      // Read out track vector in terms of pca and the position of its endpoints
      *trackstart = FindTrackStartGrad(*i);
      *trackend = FindTrackEndGrad(*i);
      range = TMath::Abs(trackend->X()-trackstart->X());
      slice_num= (Int_t)range*SLICE_DENS;
      slice_width = range/slice_num;
      for (Double_t slice=0; slice!= slice_num; slice++)
	{
	  slice_lb = trackstart->X() + (slice*slice_width);
	  slice_up = trackstart->X() + ((slice+1)*slice_width);
	  // ITERATE OVER HITS
	  for (auto j = i->begin(); j!=i->end(); j++)
	    {
	      if ( j->X() < slice_up && j->X() > slice_lb )
		{
		  hit_rad->push_back(std::sqrt(TMath::Power(j->Y(),2)+TMath::Power(j->Z(),2)));
		}

	    }
	  for (unsigned int k = 0; k!=hit_rad->size(); k++)
	    {
	      hit_tot+=hit_rad->at(k);
	    }
	  slice_mean = hit_tot/hit_rad->size();
	  for (unsigned int l = 0; l!=hit_rad->size(); l++)
	    {
	      slide_sq_tot+=(TMath::Power((hit_rad->at(l)-slice_mean),2));
	    }
	  slice_sd = std::sqrt(slide_sq_tot/hit_rad->size());
	  sd_grad_piplus = (slice_sd - slice_sd_prev)/slice_width;

	  // ADD SD TO NTUPLE HERE
	  hist_bg_grad->Fill(sd_grad_piplus);
	  grad_hold_bg->push_back(sd_grad_piplus);
	  slice_sd_prev = slice_sd;
	  hit_rad->clear();
	  hit_tot = 0; slice_sd_prev = 0; slide_sq_tot = 0;
	}

    }

  std::vector<Double_t>* grad_hold_electron = new std::vector<Double_t>;

  // electron calculation
  for (auto i = trackcon_electron.begin(); i!=trackcon_electron.end(); i++ )
    {
      // Read out track vector in terms of pca and the position of its endpoints
      *trackstart = FindTrackStartGrad(*i);
      *trackend = FindTrackEndGrad(*i);
      range = TMath::Abs(trackend->X()-trackstart->X());
      slice_num= (Int_t)range*SLICE_DENS;
      slice_width = range/slice_num;
      for (Double_t slice=0; slice!= slice_num; slice++)
	{
	  slice_lb = trackstart->X() + (slice*slice_width);
	  slice_up = trackstart->X() + ((slice+1)*slice_width);
	  // ITERATE OVER HITS
	  for (auto j = i->begin(); j!=i->end(); j++)
	    {
	      if ( j->X() < slice_up && j->X() > slice_lb )
		{
		  hit_rad->push_back(std::sqrt(TMath::Power(j->Y(),2)+TMath::Power(j->Z(),2)));
		}

	    }
	  for (unsigned int k = 0; k!=hit_rad->size(); k++)
	    {
	      hit_tot+=hit_rad->at(k);
	    }
	  slice_mean = hit_tot/hit_rad->size();
	  for (unsigned int l = 0; l!=hit_rad->size(); l++)
	    {
	      slide_sq_tot+=(TMath::Power((hit_rad->at(l)-slice_mean),2));
	    }
	  slice_sd = std::sqrt(slide_sq_tot/hit_rad->size());
	  sd_grad_electron = (slice_sd - slice_sd_prev)/slice_width;

	  // ADD SD TO NTUPLE HERE
	  hist_signal_grad->Fill(sd_grad_electron);
	  grad_hold_electron->push_back(sd_grad_electron);
	  slice_sd_prev = slice_sd;
	  hit_rad->clear();
	  hit_tot = 0; slice_sd_prev = 0; slide_sq_tot = 0;
	}

    }
  TTree* VarTree = new TTree("GradTree","GradTree");
  VarTree->Branch("grad_bg", &grad_bg, "grad_bg/D");
  VarTree->Branch("grad_signal", &grad_signal, "grad_signal/D");
  Int_t nentries_fill;
  if(grad_hold_bg->size()>grad_hold_electron->size()) { nentries_fill =  grad_hold_electron->size(); }
  else { nentries_fill = grad_hold_bg->size(); }

  for (int i = 0; i<nentries_fill; i++ )
    {
      if (grad_hold_bg->at(i) && grad_hold_electron->at(i) ) 
	{
	  grad_signal = grad_hold_electron->at(i);
	  grad_bg = grad_hold_bg->at(i);
	  VarTree->Fill();
	}
    }
  GradFile->cd();
  VarTree->Write();
  

  TCanvas* mycan = new TCanvas("mycan","Grad Variable",1600,100,1200,1200);
  TLegend* myleg = new TLegend(0.75,0.75,0.95,0.95);
  myleg->AddEntry(hist_signal_grad);
  myleg->AddEntry(hist_bg_grad);
  
  hist_bg_grad->SetFillColor(kBlue);
  hist_signal_grad->SetFillColor(kRed);
  hist_signal_grad->SetMaximum(1);
  hist_signal_grad->GetXaxis()->SetTitle("Gradient Variable");

  if (hist_bg_grad->Integral()!=0) { hist_bg_grad->Scale(1/hist_bg_grad->Integral()); }
  if (hist_signal_grad->Integral()!=0) {hist_signal_grad->Scale(1/hist_signal_grad->Integral()); }
  hist_bg_grad->SetDirectory(0);
  hist_signal_grad->SetDirectory(0);
  hist_signal_grad->Draw();
  hist_bg_grad->Draw("SAME");
  hist_bg_grad->Write();
  hist_signal_grad->Write();
  GradFile->Close();
    
}
