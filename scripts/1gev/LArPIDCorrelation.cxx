// A plot of the four proposed 'shower-shape' PID variables correlation

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
#include "TH2D.h"
#include "TLegend.h"
#include "TAxis.h"
#include <vector>
#include "TLorentzVector.h"
#include "TRandom.h"
#include </data/t2k/phumhh/lbne/srcs/larvalidatereco/larvalidatereco/content/data_objects/LArAnalysis.h>
#include </data/t2k/phumhh/lbne/srcs/larvalidatereco/larvalidatereco/content/data_objects/LArValidation.h>

using namespace std;

void LArPIDCorrelation()

{

  TH2D* grad_spread = new TH2D("grad_spread", "Grad v. Spread", 50,0,10,50,0,0.6);
  grad_spread->SetStats(0);
  grad_spread->GetXaxis()->SetTitle("Grad Variable");
  grad_spread->GetYaxis()->SetTitle("Spread Variable");
  grad_spread->GetYaxis()->SetTitleOffset(1.4);
  TH2D* grad_con = new TH2D("grad_con", "Grad v. Con", 50,0,10,50,0.6,1);
  grad_con->SetStats(0);
  grad_con->GetXaxis()->SetTitle("Grad Variable");
  grad_con->GetYaxis()->SetTitle("Conicalness Variable");
  grad_con->GetYaxis()->SetTitleOffset(1.4);
  TH2D* grad_corehalo = new TH2D("grad_corehalo", "Grad v. CoreHalo", 50,0,10,50,0,0.02);
  grad_corehalo->SetStats(0);
  grad_corehalo->GetXaxis()->SetTitle("Grad Variable");
  grad_corehalo->GetYaxis()->SetTitle("CoreHalo Variable");
  grad_corehalo->GetYaxis()->SetTitleOffset(1.6);
  TH2D* spread_con = new TH2D("spread_con", "Spread v. Con", 50,0,0.6,50,0.6,1);
  spread_con->SetStats(0);
  spread_con->GetXaxis()->SetTitle("Spread Variable");
  spread_con->GetYaxis()->SetTitle("Conicalness Variable");
  spread_con->GetYaxis()->SetTitleOffset(1.4);
  TH2D* spread_corehalo = new TH2D("spread_corehalo", "Spread v. CoreHalo",50,0,0.6,50,0,0.02);
  spread_corehalo->SetStats(0);
  spread_corehalo->GetXaxis()->SetTitle("Spread Variable");
  spread_corehalo->GetYaxis()->SetTitle("CoreHalo Variable");
  spread_corehalo->GetYaxis()->SetTitleOffset(1.6);
  TH2D* con_corehalo = new TH2D("con_corehalo", "Con v. CoreHalo",50,0.6,1,50,0,0.02);
  con_corehalo->SetStats(0);
  con_corehalo->GetXaxis()->SetTitle("Conicalness Variable");
  con_corehalo->GetYaxis()->SetTitle("CoreHalo Variable");
  con_corehalo->GetYaxis()->SetTitleOffset(1.6);
  
  Double_t grad_h;
  Double_t spread_h; 
  Double_t corehalo_h; 
  Double_t con_h;  


  TFile* grad_f = new TFile("pid_files/VarGrad.root");
  TFile* spread_f = new TFile("pid_files/VarSpread.root");
  TTree* grad_t = (TTree*)grad_f->Get("GradTree");
  TTree* spread_t = (TTree*)spread_f->Get("SpreadTree");
  grad_t->SetBranchAddress("grad_signal",&grad_h);
  spread_t->SetBranchAddress("spread_signal",&spread_h);
  Int_t nentries = 0;
  if ( grad_t->GetEntries() < spread_t->GetEntries()) {nentries=grad_t->GetEntries();}
  else {nentries=spread_t->GetEntries();}
  for (int i=0; i<nentries; i++)
    {
      spread_t->GetEntry(i);
      grad_t->GetEntry(i);
      grad_spread->Fill(grad_h,spread_h);
    }
  grad_spread->SetDirectory(0);
  grad_f->Close();
  spread_f->Close();


  grad_f = new TFile("pid_files/VarGrad.root");
  TFile* con_f = new TFile("pid_files/VarCon.root");
  TTree* con_t = (TTree*)con_f->Get("ConTree");
  grad_t = (TTree*)grad_f->Get("GradTree");
  con_t->SetBranchAddress("con_signal",&con_h);
  grad_t->SetBranchAddress("grad_signal",&grad_h);

  if ( grad_t->GetEntries() < con_t->GetEntries()) {nentries=grad_t->GetEntries();}
  else {nentries=con_t->GetEntries();}

for (int j=0; j<nentries; j++ )
    {
      grad_t->GetEntry(j);
      con_t->GetEntry(j);
      grad_con->Fill(grad_h,con_h);
    }
  grad_con->SetDirectory(0);

  grad_f = new TFile("pid_files/VarGrad.root");
  TFile * corehalo_f = new TFile("pid_files/VarCoreHalo.root");
  grad_t = (TTree*)grad_f->Get("GradTree");
  TTree* corehalo_t = (TTree*)corehalo_f->Get("CoreHaloTree");
  grad_t->SetBranchAddress("grad_signal",&grad_h);
  corehalo_t->SetBranchAddress("corehalo_signal",&corehalo_h);
  if ( grad_t->GetEntries() < corehalo_t->GetEntries()) {nentries=grad_t->GetEntries();}
  else {nentries=corehalo_t->GetEntries();}

  for (int k=0; k<nentries; k++ )
    {
      grad_t->GetEntry(k);
      corehalo_t->GetEntry(k);
      grad_corehalo->Fill(grad_h,corehalo_h);
    }
  grad_corehalo->SetDirectory(0);
  grad_f->Close();
  corehalo_f->Close();

  con_f = new TFile("pid_files/VarCon.root");
  con_t = (TTree*)con_f->Get("ConTree");
  spread_f = new TFile("pid_files/VarSpread.root");
  spread_t = (TTree*)spread_f->Get("SpreadTree");
  con_t->SetBranchAddress("con_signal",&con_h);
  spread_t->SetBranchAddress("spread_signal",&spread_h);
  if ( con_t->GetEntries() < spread_t->GetEntries()) {nentries=con_t->GetEntries();}
  else {nentries=spread_t->GetEntries();}

  for (int l=0; l<nentries; l++)
    {
      spread_t->GetEntry(l);
      con_t->GetEntry(l);
      spread_con->Fill(spread_h,con_h);
    }
  spread_con->SetDirectory(0);
  spread_f->Close();
  con_f->Close();

  spread_f = new TFile("pid_files/VarSpread.root");
  spread_t = (TTree*)spread_f->Get("SpreadTree");
  corehalo_f = new TFile("pid_files/VarCoreHalo.root");
  corehalo_t = (TTree*)corehalo_f->Get("CoreHaloTree");
  spread_t->SetBranchAddress("spread_signal",&spread_h);
  corehalo_t->SetBranchAddress("corehalo_signal",&corehalo_h);
  if ( corehalo_t->GetEntries() < spread_t->GetEntries()) {nentries=corehalo_t->GetEntries();}
  else {nentries=spread_t->GetEntries();}

  for (int m=0; m<nentries; m++)
    {
      spread_t->GetEntry(m);
      corehalo_t->GetEntry(m);
      spread_corehalo->Fill(spread_h,corehalo_h);
    }
  spread_corehalo->SetDirectory(0);
  spread_f->Close();
  corehalo_f->Close();

  corehalo_f = new TFile("pid_files/VarCoreHalo.root");
  corehalo_t = (TTree*)corehalo_f->Get("CoreHaloTree");
  con_f = new TFile("pid_files/VarCon.root");
  con_t = (TTree*)con_f->Get("ConTree");
  con_t->SetBranchAddress("con_signal",&con_h);
  corehalo_t->SetBranchAddress("corehalo_signal",&corehalo_h);
  if ( con_t->GetEntries() < corehalo_t->GetEntries()) {nentries=con_t->GetEntries();}
  else {nentries=corehalo_t->GetEntries();}
  for (int n=0; n<nentries; n++ )
    { 
      con_t->GetEntry(n);
      corehalo_t->GetEntry(n);
      con_corehalo->Fill(con_h,corehalo_h);
    }
  con_corehalo->SetDirectory(0);
  con_f->Close();
  corehalo_f->Close();


  TCanvas* mycan = new TCanvas("mycan", "Variable Correlations", 1600,100,1500,1200);
  mycan->Divide(3,2);
  mycan->cd(1);
  grad_spread->Draw("COLZ");
  mycan->cd(2);
  grad_con->Draw("COLZ");
  mycan->cd(3);
  grad_corehalo->Draw("COLZ");
  mycan->cd(4);
  spread_con->Draw("COLZ");
  mycan->cd(5);
  spread_corehalo->Draw("COLZ");
  mycan->cd(6);
  con_corehalo->Draw("COLZ");
  mycan->cd();
  mycan->Update();
  


  TFile* corr_f = new TFile("pid_files/VarCorrelate.root", "RECREATE");
  corr_f->cd();
  grad_spread->Write();
  grad_con->Write();
  grad_corehalo->Write();
  spread_con->Write();
  spread_corehalo->Write();
  con_corehalo->Write();
  corr_f->Close();
  
  con_f->Close();
  corehalo_f->Close();
  spread_f->Close();
  grad_f->Close();
}
