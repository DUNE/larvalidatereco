#include "stdlib.h"
#include "iostream"
#include "string"
#include "TMultiLayerPerceptron.h"
#include "TMLPAnalyzer.h"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TBranch.h"

void LArPIDMultiVariate(Int_t ntrain = 800)
{
  gSystem->Load("libMLP");

  // Load variable files
  // Load variable trees
  // Set variable tree addresses
  std::vector<Double_t>* con_arr_sig =    new std::vector<Double_t>;
  std::vector<Double_t>* con_arr_bg =     new std::vector<Double_t>;
  std::vector<Double_t>* spread_arr_sig = new std::vector<Double_t>;
  std::vector<Double_t>* spread_arr_bg =  new std::vector<Double_t>;
  std::vector<Double_t>* grad_arr_sig =   new std::vector<Double_t>;
  std::vector<Double_t>* grad_arr_bg =    new std::vector<Double_t>;

  Double_t con_cont_sig, con_cont_bg, grad_cont_sig, grad_cont_bg, spread_cont_sig, spread_cont_bg;
  Double_t con_load, spread_load, grad_load;
  Int_t type;

  // Load con
  TFile* con = new TFile("pid_files/VarCon.root");
  TTree* con_tree = (TTree*)con->Get("ConTree");
  con->cd();
  con_tree->SetBranchAddress("con_bg",&con_cont_sig);
  con_tree->SetBranchAddress("con_signal",&con_cont_bg);
  for (int i = 0 ; i < con_tree->GetEntries(); i++ )
    {
      con_tree->GetEntry(i);
      con_arr_sig->push_back(con_cont_sig);

    }
  for ( int j = 0; j < con_tree->GetEntries(); j++ ) 
    {
      con_tree->GetEntry(j);
      con_arr_bg->push_back(con_cont_bg);
    }
  con->Close();
  // Load Spread
  TFile* spread = new TFile("pid_files/VarSpread.root");
  TTree* spread_tree = (TTree*)spread->Get("SpreadTree");
  spread->cd();
  spread_tree->SetBranchAddress("spread_bg",&spread_cont_sig);
  spread_tree->SetBranchAddress("spread_signal",&spread_cont_bg);
  for ( int k = 0 ; k < spread_tree->GetEntries(); k++ )
    {
      spread_tree->GetEntry(k);
      spread_arr_sig->push_back(spread_cont_sig);
    }
  type=0;
  for ( int l = 0; l < spread_tree->GetEntries(); l++ )
    {
      spread_tree->GetEntry(l);
      spread_arr_bg->push_back(spread_cont_bg);
    }
  spread->Close();

  // Load Grad
  TFile* grad = new TFile("pid_files/VarGrad.root");
  TTree* grad_tree = (TTree*)grad->Get("GradTree");
  grad->cd();
  grad_tree->SetBranchAddress("grad_bg",&grad_cont_sig);
  grad_tree->SetBranchAddress("grad_signal",&grad_cont_bg);
  for ( int m = 0 ; m< grad_tree->GetEntries(); m++ )
    {
      grad_tree->GetEntry(m);
      grad_arr_sig->push_back(grad_cont_sig);
    }

  for ( int n = 0; n < grad_tree->GetEntries(); n++ )
    {
      grad_tree->GetEntry(n);
      grad_arr_bg->push_back(grad_cont_bg);
    }
  grad->Close();

  // Create new combined tree
  TTree* nn_tree = new TTree("nn_tree","Input Tree");
  // Create combined tree branches
  nn_tree->Branch("con",&con_load,"con/D");
  nn_tree->Branch("spread",&spread_load,"spread/D");
  nn_tree->Branch("grad",&grad_load,"grad/D");
  nn_tree->Branch("type",&type,"type/I");


  // Now Load All into Input Tree
  for (unsigned int i = 0; i < con_arr_sig->size(); i++ ) 
    {
      type=1;
      con_load = con_arr_sig->at(i);
      nn_tree->Fill();
    }   
  for (unsigned int j = 0; j< con_arr_bg->size(); j++)
    {
      type=0;
      con_load = con_arr_bg->at(j);
      nn_tree->Fill();
    }
  for (unsigned int k = 0; k< grad_arr_sig->size(); k++)
    {
      type=1;
      grad_load = grad_arr_sig->at(k);
      nn_tree->Fill();
      }
      for (unsigned int l = 0 ; l< grad_arr_bg->size(); l++)
      {
      type=0;
      grad_load = grad_arr_bg->at(l);
      nn_tree->Fill();
    }
  for (unsigned int m=0; m< spread_arr_sig->size(); m++)
    {
      type=1;
      spread_load = spread_arr_sig->at(m);
      nn_tree->Fill();
      }
      for (unsigned int n=0; n<spread_arr_bg->size(); n++)
      {
      type=0;
      spread_load = spread_arr_bg->at(n);
      nn_tree->Fill();
    }


  std::cout << "TOTAL ENTRIES: " << nn_tree->GetEntries() << std::endl << "Half: " << nn_tree->GetEntries()/2 << std::endl;

  TMultiLayerPerceptron * mlp = 
    new TMultiLayerPerceptron("@con,@spread,@grad:3:type",nn_tree,"Entry$%2","(Entry$+1)%2");
  mlp->Train(ntrain, "text,graph,update=10");
  mlp->Export( "NNtest","C++" );

  TCanvas* mlpa_canvas = new TCanvas("mlpa_canvas","Network Analysis");
  mlpa_canvas->Divide(2,2);
  TMLPAnalyzer ana(mlp);
  
  ana.GatherInformations();
  ana.CheckNetwork();
  // DRAWS THE EFFECT OF EACH VARIABLE OVER EACH CYCLE
  mlpa_canvas->cd(1);
  ana.DrawDInputs();
  // DRAWS THE OUTPUT OF THE SIGNAL AND BACKGROUND DISCRIMINATION
  mlpa_canvas->cd(2);
  mlp->Draw();
  // DRAWS THE CONTRIBUTIONS OF EACH VARIABLE TO THE NETWORK
  mlpa_canvas->cd(3);
  ana.DrawNetwork(0,"type==1","type==0");

  //mlpa_canvas->cd(4);
  


}
