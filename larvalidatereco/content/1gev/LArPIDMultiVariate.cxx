#include "stdlib.h"
#include "string"
#include "TMultiLayerPerceptron.h"
#include "TMLPAnalyzer.h"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TSystem.h"

void LArPIDMultiVariate(Int_t ntrain = 100)
{
  gSystem->Load("libMLP");

  // Load variable files
  TFile* con = new TFile("pid_files/VarCon.root");
  TFile* grad = new TFile("pid_files/VarGrad.root");
  TFile* spread = new TFile("pid_files/VarSpread.root");
  // Load variable trees
  TTree* con_tree = (TTree*)con->Get("ConTree");
  TTree* grad_tree = (TTree*)grad->Get("GradTree");
  TTree* spread_tree = (TTree*)spread->Get("SpreadTree");
  // Create new combined tree
  TTree* nn_tree = new TTree("nn_tree","Input Tree");
  Double_t con_cont, grad_cont, spread_cont;
  Int_t type;
  // Set variable tree addresses





;
  // Create combined tree branches
  nn_tree->Branch("con",&con_cont,"con/D");
  nn_tree->Branch("spread",&spread_cont,"spread/D");
  nn_tree->Branch("grad",&grad_cont,"grad/D");
  nn_tree->Branch("type",&type,"type/I");


  Int_t i,j,k,l,m,n;
  // Load con
  con->cd();
  type=1;
  con_tree->SetBranchAddress("con_sig",&con_cont);
  //  con_tree->SetBranchStatus("con_bg",0);
  for ( i = 0 ; i < con_tree->GetEntries(); i++ )
    {
      con_cont = con_tree->GetEntry(i);
      if (con_tree->GetEntry(i)){nn_tree->Fill();}
    }
  type=0;
  con_tree->SetBranchAddress("con_bg",&con_cont);
  //con_tree->SetBranchStatus("con_bg",1);
  //con_tree->SetBranchStatus("con_sig",0);
  for ( j = 0; j < con_tree->GetEntries(); j++ ) 
    {
      con_cont = con_tree->GetEntry(j);
      if(con_tree->GetEntry(j)){nn_tree->Fill();}
    }
  spread->cd();
  type=1;
  spread_tree->SetBranchAddress("spread_sig",&spread_cont)
  // Load Spread
  //spread_tree->SetBranchStatus("spread_bg",0);
  for ( i = 0 ; i < spread_tree->GetEntries(); i++ )
    {
      spread_cont = spread_tree->GetEntry(i);
      if(spread_tree->GetEntry(i)){nn_tree->Fill();}
    }
  type=0;
  spread_tree->SetBranchAddress("spread_bg",&spread_cont);
  //spread_tree->SetBranchStatus("spread_bg",1);
  //spread_tree->SetBranchStatus("spread_sig",0);
  for ( i = 0; i < spread_tree->GetEntries(); i++ )
    {
      if(spread_tree->GetEntry(i)){spread_cont = spread_tree->GetEntry(i);}
      nn_tree->Fill();
    }
  grad->cd();
  type=1;
  grad_tree->SetBranchAddress("grad_sig",&grad_cont);
  // Load Grad
  //grad_tree->SetBranchStatus("grad_bg",0);
  for ( i = 0 ; i < grad_tree->GetEntries(); i++ )
    {
      if(grad_tree->GetEntry(i)){grad_cont = grad_tree->GetEntry(i);}
      nn_tree->Fill();
    }
  type=0;
  grad_tree->SetBranchAddress("grad_bg",&grad_cont);
  //grad_tree->SetBranchStatus("grad_bg",1);
  //grad_tree->SetBranchStatus("grad_sig",0);
  for ( i = 0; i < grad_tree->GetEntries(); i++ )
    {
      if(grad_tree->GetEntry(i)){grad_cont = grad_tree->GetEntry(i);}
      nn_tree->Fill();
    }
  
  TMultiLayerPerceptron * mlp = 
    new TMultiLayerPerceptron("@con,@spread,@grad:3:3:type","spread",nn_tree,"Entry$%2","(Entry$+1)%2");
  mlp->Train(ntrain, "text,graph,update=10");
  //  mlp->Export( WHAT GOES HERE );

  TCanvas* mlpa_canvas = new TCanvas("mlpa_canvas","Network Analysis");
  mlpa_canvas->Divide(2,2);
  TMLPAnalyzer ana(mlp);
  
  ana.GatherInformations();
  ana.CheckNetwork();
  
  mlpa_canvas->cd(1);
  ana.DrawDInputs();
  
  mlpa_canvas->cd(2);
  mlp->Draw();
  
  mlpa_canvas->cd(3);
  ana.DrawNetwork(0,"type==1","type==0");

  //mlpa_canvas->cd(4);
  


}
