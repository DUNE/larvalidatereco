#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TAxis.h"

void LArPIDCheck2Draw()

{
  /*
  TFile HistFile("pid_files/MomVEigenHist2.root");
  TH2D* MuonHist = HistFile.Get("MuonHist");
  TH2D* ElectronHist = HistFile.Get("ElectronHist");
  TH2D* ProtonHist = HistFile.Get("ProtonHist");
  TH2D* PiPlusHist = HistFile.Get("PiPlusHist");
  MuonHist->SetDirectory(0);
  MuonHist->SetStats(0);
  ElectronHist->SetDirectory(0);
  ElectronHist->SetStats(0);
  ProtonHist->SetDirectory(0);
  ProtonHist->SetStats(0);
  PiPlusHist->SetDirectory(0);
  PiPlusHist->SetStats(0);
  TCanvas* mycan = new TCanvas("mycan","MomVEigen Histograms",1600,200,1200,1200);
  mycan->Divide(2,2);
  mycan->cd(1);
  MuonHist->Draw("COLZ");
  mycan->cd(2);
  ElectronHist->Draw("COLZ");
  mycan->cd(3);
  ProtonHist->Draw("COLZ");
  mycan->cd(4);
  PiPlusHist->Draw("COLZ");
  */
}
