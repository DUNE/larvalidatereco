
// This file generates an optimized cut of the data stored in /pid_files/VarCoreHalo.root  
#include "TCut.h"
#include "TH1D.h"
#include "TAxis.h"
#include "string"
#include "stdio.h"
#include "iostream"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TCanvas.h"
#include "TLegend.h"

#define CUT_THRESH 6
#define STEP_SIZE 0.4
#define BIN_NUM 50
#define MAX 0.04
void LArPIDCutCoreHalo()
{

  TFile* f = new TFile("pid_files/VarCoreHalo.root");
  TTree* t = (TTree*)f->Get("CoreHaloTree");
  TH1D* bg = new TH1D("bg","CoreHalo Variable", BIN_NUM,0,MAX);
  TH1D* sig = new TH1D("sig","Signal",BIN_NUM,0,MAX);
  bg->GetXaxis()->SetTitle("CoreHalo");
  bg->GetYaxis()->SetTitle("Arbitrary Units");
  bg->SetStats(0);
  bg->SetMaximum(300);
  bg->SetFillColor(kBlue);
  sig->SetFillColor(kRed);
  sig->SetFillStyle(3004);
  sig->SetLineColor(kRed);
  
  Double_t bg_h, sig_h;
  t->SetBranchAddress("corehalo_bg",&bg_h);
  t->SetBranchAddress("corehalo_signal",&sig_h);

  Int_t entries = t->GetEntries();
  for (int i = 0; i<entries; i++)
    {
      t->GetEntry(i);
      bg->Fill(bg_h);
      sig->Fill(sig_h);
    }

  //Int_t hist_tot = bg->Integral();
  //  sig->Scale(1/hist_tot);

  // Compute optimal cut by taking the difference between signal and bg bin height, then dividing by the number of entries in the proposed cut.

  Double_t diff[BIN_NUM];
  Double_t cut_num[BIN_NUM]; 
  Double_t diff_err[BIN_NUM];
  Double_t sig_tot, bg_tot;
  Double_t sig_cont, bg_cont;
  for ( int l = 25; l>1; l--)
	{
	  for (int m=l; m<25; m++)
	    {
	      sig_cont = sig->GetBinContent(m);
	      bg_cont = bg->GetBinContent(m);
	      if (sig_cont > bg_cont )
		{
		  sig_tot+=sig_cont;
		  bg_tot+=bg_cont;
		}
	    }
	  //sig_tot = sig->Integral(l,BIN_NUM);
	  //bg_tot = bg->Integral(l,BIN_NUM);
	  if ( sig_tot > bg_tot ) { diff[l] = (sig_tot-bg_tot)/(sig->Integral(l,BIN_NUM)+bg->Integral(l,BIN_NUM)); }
	  else { diff[l] = 0; 
	  }
	  if (bg->GetBinContent(l) < sig->GetBinContent(l)) {
            diff_err[l] =
              std::sqrt(
			TMath::Power
			(
			 (
			  std::sqrt(
				    sig->GetBinContent(l) + bg->GetBinContent(l))
			  )/
			 (sig->GetBinContent(l)-bg->GetBinContent(l))
			 ,2)
              +
			TMath::Power
			(
			 (
			  std::sqrt(
				    sig->GetBinContent(l) + bg->GetBinContent(l))
			  )/
			 (sig->GetBinContent(l)+bg->GetBinContent(l))
			 ,2)
                        )/5;
	  } else { diff_err[l] = 0; }

	  cut_num[l] = l*0.0008;
	}

  // Produce a graph of the cuts, and find the maximum enry in the graph, and use this as the optimal cut. 

  std::cout << std::endl << "Complete loop" << std::endl;
  TGraphErrors* optimal = new TGraphErrors(BIN_NUM,cut_num,diff,NULL,diff_err);
  optimal->GetXaxis()->SetTitle("# Cut");
  optimal->GetYaxis()->SetTitle("Sig - Back");
  optimal->GetXaxis()->CenterTitle();
  optimal->GetYaxis()->CenterTitle();
  optimal->SetTitle("Cut Optimization");

  optimal->SetPointError(5,0,0.3);
  optimal->SetPointError(15,0,0.3);
  optimal->SetPointError(23,0,0);

  Double_t cut_optimal;
  Double_t opt = TMath::MaxElement(optimal->GetN(),optimal->GetY());
  for (int m = 0; m<optimal->GetN(); m++ ) { if ( diff[m] == opt ) { cut_optimal = cut_num[m]; } }
  
  std::cout << std::endl << "Optimal Cut: " << cut_optimal << std::endl << std::endl;
    
  // Fill histograms with the cut in place

  TH1D* bg_cut = new TH1D("bg_cut","Post-Cut CoreHalo Variable",BIN_NUM,0,MAX);
  TH1D* sig_cut = new TH1D("sig_cut","Signal CoreHalo",BIN_NUM,0,MAX);
  bg_cut->GetXaxis()->SetTitle("CoreHalo");
  bg_cut->GetYaxis()->SetTitle("Arbitrary Units");
  bg_cut->SetStats(0);
  bg_cut->SetMaximum(300);
  for (int k = 0; k<entries; k++ )
    {
      t->GetEntry(k);
      if ( bg_h>cut_optimal) { bg_cut->Fill(bg_h); }
      if ( sig_h>cut_optimal) { sig_cut->Fill(sig_h); }
    }
  // Compute purity and efficiency of plots   
  Double_t pur, eff, eff_err, pur_err;
  pur = (sig_cut->GetEntries())/(sig_cut->GetEntries()+bg_cut->GetEntries());
  eff = (sig_cut->GetEntries())/(sig->GetEntries());
  eff_err = std::sqrt((eff*(1-eff))/(bg->GetEntries()+sig->GetEntries()));
  pur_err = std::sqrt((pur*(1-pur))/(bg->GetEntries()+sig->GetEntries()));
  std::cout << std::endl << "Purity: " << pur << " | Err: " << pur_err << std::endl << "Efficiency: " << eff*100 << " | err: "   << eff_err*100 << std::endl;

  TCanvas *mycan = new TCanvas("mycan","Grad Cut", 1600,100,1500,500);
  TLegend *leg = new TLegend(0.75,0.75,0.95,0.95);
  leg->AddEntry(bg, "Background");
  leg->AddEntry(sig, "Signal");
  TLegend *leg_cut = new TLegend(0.75,0.75,0.95,0.95);
  leg_cut->AddEntry(bg_cut, "Background");
  leg_cut->AddEntry(sig_cut, "Signal");
  TPaveText* cut_param = new TPaveText(0.4,0.4,0.6,0.6);
  cut_param->AddText("Purity:");
  cut_param->AddText("Efficiency");

  mycan->Divide(2,1);
  mycan->cd(1);
  bg->Draw();
  sig->Draw("SAME");
  leg->Draw("SAME");
  
  
  mycan->cd(2);

  bg_cut->SetFillColor(kBlue);
  sig_cut->SetFillColor(kRed);
  sig_cut->SetFillStyle(3004);
  sig_cut->SetLineColor(kRed);
  bg_cut->Draw();
  sig_cut->Draw("SAME");
  leg_cut->Draw("SAME");
  cut_param->Draw("SAME");
  mycan->cd();
  mycan->Update();
  TCanvas* r_can = new TCanvas("r_can", "Cut Optimization", 700,500);
  r_can->cd();
  optimal->Draw("AP*");


}
