// USAGE: Fitter("../histograms/experimental/hist_ditau_mass_topocut.root", "savename")
// Fitter("../histograms/experimental/hist_ditau_mass.root", "savename")


#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TLegend.h>
#include <TStyle.h>

#include <cmath>
#include <vector>
#include <string>
#include <utility>
#include <iostream>


void Fitter(const char *rootfile, const char *outname, int nbins) {
  gStyle->SetOptStat(0);
  //gStyle->SetOptFit(2222);
  gROOT->SetBatch(kTRUE);


  TCanvas *c1 = new TCanvas("c3", "", 1280, 960);
  TCanvas *c2 = new TCanvas("c4", "", 1280, 960);

  // Load in the TTree root file.
  TFile *file_in = TFile::Open(rootfile);

  // Read in the background TH1's with fine binning.
  THStack *bg_stacks = new THStack("bg","");
  TH1D *h1 = (TH1D*)file_in->Get("ttZ");
  TH1D *h2 = (TH1D*)file_in->Get("ttW");
  TH1D *backgrounds = (TH1D*)file_in->Get("ttZ");
  backgrounds->Add(h2);




  // Load signal
  TH1D *s40 = (TH1D*)file_in->Get("ttZp-40GeV");
  TH1D *s50 = (TH1D*)file_in->Get("ttZp-50GeV");
  TH1D *s70 = (TH1D*)file_in->Get("ttZp-70GeV");
  TH1D *s100 = (TH1D*)file_in->Get("ttZp-100GeV");
  TH1D *s120 = (TH1D*)file_in->Get("ttZp-120GeV");
  TH1D *s150 = (TH1D*)file_in->Get("ttZp-150GeV");


  // Rebin
  TH1F *backgrounds_flat = (TH1F*)backgrounds->Rebin(nbins, "Flat");   ///// NOTE: combining backgrounds here.
  TH1F *h1_flat = (TH1F*)h1->Rebin(nbins, "Flat");
  TH1F *h2_flat = (TH1F*)h2->Rebin(nbins, "Flat");
  TH1F *s40_flat = (TH1F*)s40->Rebin(nbins, "Flat");
  TH1F *s50_flat = (TH1F*)s50->Rebin(nbins, "Flat");
  TH1F *s70_flat = (TH1F*)s70->Rebin(nbins, "Flat");
  TH1F *s100_flat = (TH1F*)s100->Rebin(nbins, "Flat");
  TH1F *s120_flat = (TH1F*)s120->Rebin(nbins, "Flat");
  TH1F *s150_flat = (TH1F*)s150->Rebin(nbins, "Flat");


  // Perform fits.

  // Possibilities: Gaussian, Beta distribution, or Landau distribution.
  TF1 *f40 = new TF1("f40","[0]*ROOT::Math::crystalball_function(x,[1],[2],[3],[4])",15,200);
  TF1 *f50 = new TF1("f50","[0]*ROOT::Math::crystalball_function(x,[1],[2],[3],[4])",15,200);
  TF1 *f70 = new TF1("f70","[0]*ROOT::Math::crystalball_function(x,[1],[2],[3],[4])",15,200);
  TF1 *f100 = new TF1("f100","[0]*ROOT::Math::crystalball_function(x,[1],[2],[3],[4])",15,200);
  TF1 *f120 = new TF1("f120","[0]*ROOT::Math::crystalball_function(x,[1],[2],[3],[4])",15,200);
  TF1 *f150 = new TF1("f150","[0]*ROOT::Math::crystalball_function(x,[1],[2],[3],[4])",15,200);
  TF1 *fTTV = new TF1("fTTV","[0]*ROOT::Math::crystalball_function(x,[1],[2],[3],[4])",15,200);
  TF1 *fTTW = new TF1("fTTW","[0]*ROOT::Math::crystalball_function(x,[1],[2],[3],[4])",15,200);
  TF1 *fBG = new TF1("fBG","[0]*ROOT::Math::crystalball_function(x,[1],[2],[3],[4])",15,200);

  fBG->SetParNames("A","#alpha","N","#sigma","mean");
  fTTV->SetParNames("A","#alpha","N","#sigma","mean");
  fTTW->SetParNames("A","#alpha","N","#sigma","mean");
  f40->SetParNames("A","#alpha","N","#sigma","mean");
  f50->SetParNames("A","#alpha","N","#sigma","mean");
  f70->SetParNames("A","#alpha","N","#sigma","mean");
  f100->SetParNames("A","#alpha","N","#sigma","mean");
  f120->SetParNames("A","#alpha","N","#sigma","mean");
  f150->SetParNames("A","#alpha","N","#sigma","mean");
  fBG->SetParameters(20.,-1.2,2.,30.0,58.);
  fTTV->SetParameters(20.,-1.2,2.,30.0,58.);
  fTTW->SetParameters(4.,-0.1,2.,50.0,100.);
  f40->SetParameters(14.8,-.9,.8,4.0,26.);
  f50->SetParameters(12.,-1.3,.9,8.6,34.);
  f70->SetParameters(8.,-1.2,.9,12.0,44.);
  f100->SetParameters(15.,-.7,1.9,15.0,58.);
  f120->SetParameters(10.,-0.8,2.9,20.0,70.);
  f150->SetParameters(8.,-1.5,1.5,28.0,95.);
  float fixed_n = 1;
  fBG->SetParameter(2,fixed_n);
  fBG->SetParLimits(2,fixed_n,fixed_n);
  fTTV->SetParameter(2,fixed_n);
  fTTV->SetParLimits(2,fixed_n,fixed_n);
  fTTW->SetParameter(2,fixed_n);
  fTTW->SetParLimits(2,fixed_n,fixed_n);
  f40->SetParameter(2,fixed_n);
  f40->SetParLimits(2,fixed_n,fixed_n);
  f50->SetParameter(2,fixed_n);
  f50->SetParLimits(2,fixed_n,fixed_n);
  f70->SetParameter(2,fixed_n);
  f70->SetParLimits(2,fixed_n,fixed_n);
  f100->SetParameter(2,fixed_n);
  f100->SetParLimits(2,fixed_n,fixed_n);
  f120->SetParameter(2,fixed_n);
  f120->SetParLimits(2,fixed_n,fixed_n);
  f150->SetParameter(2,fixed_n);
  f150->SetParLimits(2,fixed_n,fixed_n);


  backgrounds_flat->Fit(fBG, "RF");
  h1_flat->Fit(fTTV, "RF");
  h2_flat->Fit(fTTW, "RF");
  s40_flat->Fit(f40,"RF");
  s50_flat->Fit(f50,"RF");
  s70_flat->Fit(f70,"RF");
  s100_flat->Fit(f100,"RF");
  s120_flat->Fit(f120,"RF");
  s150_flat->Fit(f150,"RF");


  printf("MZp=40 GeV: X^2/ndf = %f --- p = %f --- Integral = %f \n", f40->GetChisquare() / f40->GetNDF(),
         f40->GetProb(), s40_flat->Integral());
  printf("MZp=50 GeV: X^2/ndf = %f --- p = %f --- Integral = %f \n", f50->GetChisquare() / f50->GetNDF(),
         f50->GetProb(), s50_flat->Integral());
  printf("MZp=70 GeV: X^2/ndf = %f --- p = %f --- Integral = %f \n", f70->GetChisquare() / f70->GetNDF(),
         f70->GetProb(), s70_flat->Integral());
  printf("MZp=100 GeV: X^2/ndf = %f --- p = %f --- Integral = %f \n", f100->GetChisquare() / f100->GetNDF(),
         f100->GetProb(), s100_flat->Integral());
  printf("MZp=120 GeV: X^2/ndf = %f --- p = %f --- Integral = %f \n", f120->GetChisquare() / f120->GetNDF(),
         f120->GetProb(), s120_flat->Integral());
  printf("MZp=150 GeV: X^2/ndf = %f --- p = %f --- Integral = %f \n", f150->GetChisquare() / f150->GetNDF(),
         f150->GetProb(), s150_flat->Integral());
  printf("ttV: X^2/ndf = %f --- p = %f --- Integral = %f \n", fTTV->GetChisquare() / fTTV->GetNDF(),
         fTTV->GetProb(), h1_flat->Integral());
  printf("ttW: X^2/ndf = %f --- p = %f --- Integral = %f \n", fTTW->GetChisquare() / fTTW->GetNDF(),
         fTTW->GetProb(), h2_flat->Integral());
  printf("BG: X^2/ndf = %f --- p = %f --- Integral = %f \n", fBG->GetChisquare() / fBG->GetNDF(),
         fBG->GetProb(), backgrounds_flat->Integral());


  printf("alpha = [%f, %f, %f, %f, %f, %f]\n", f40->GetParameter(1), f50->GetParameter(1), f70->GetParameter(1),
                                               f100->GetParameter(1), f120->GetParameter(1), f150->GetParameter(1));
  printf("sigma = [%f, %f, %f, %f, %f, %f]\n", f40->GetParameter(3), f50->GetParameter(3), f70->GetParameter(3),
                                               f100->GetParameter(3), f120->GetParameter(3), f150->GetParameter(3));
  printf("mu = [%f, %f, %f, %f, %f, %f]\n", f40->GetParameter(4),
                                            f50->GetParameter(4),
                                            f70->GetParameter(4),
                                            f100->GetParameter(4),
                                            f120->GetParameter(4),
                                            f150->GetParameter(4));


  // Draw histograms.

  h1_flat->SetTitle("");
  h1_flat->SetLineColor(1);
  h1_flat->GetYaxis()->SetRange(0., 100.);


  // Flat dist's
  backgrounds_flat->SetFillColorAlpha(57,1);
  h1_flat->SetFillColorAlpha(42,0.6);
  h2_flat->SetFillColorAlpha(38,0.6);
  s40_flat->SetLineWidth(3);
  s50_flat->SetLineWidth(3);
  s70_flat->SetLineWidth(3);
  s100_flat->SetLineWidth(3);
  s120_flat->SetLineWidth(3);
  s150_flat->SetLineWidth(3);

  fBG->SetLineWidth(4);
  fTTV->SetLineWidth(4);
  fTTW->SetLineWidth(4);
  f40->SetLineWidth(4);
  f50->SetLineWidth(4);
  f70->SetLineWidth(4);
  f100->SetLineWidth(4);
  f120->SetLineWidth(4);
  f150->SetLineWidth(4);

  s40_flat->SetLineColor(kAzure-3);
  s50_flat->SetLineColor(3);
  s70_flat->SetLineColor(kOrange-3);
  s100_flat->SetLineColor(870);
  s120_flat->SetLineColor(6);
  s150_flat->SetLineColor(kRed-7);

  fBG->SetLineColor(13);
  fTTV->SetLineColor(42);
  fTTW->SetLineColor(38);
  f40->SetLineColor(kAzure-3);
  f50->SetLineColor(3);
  f70->SetLineColor(kOrange-3);
  f100->SetLineColor(870);
  f120->SetLineColor(6);
  f150->SetLineColor(kRed-7);

  // Draw the two backgrounds stacked together (but the fit is on the sum):
  bg_stacks->Add(h1_flat);
  bg_stacks->Add(h2_flat);


  TLegend *legend = new TLegend(0.56, 0.60, 0.78, 0.88);
  legend->SetBorderSize(0.);
  legend->SetTextSize(0.04); // % of pad size

  legend->AddEntry(h1_flat, "ttZ/H/#gamma*", "f");
  legend->AddEntry(h2_flat, "ttW", "f");
  legend->AddEntry(f40, "ttZ' (40 GeV, g_{b}=0.1)", "l");
  //legend->AddEntry(s50_flat, "ttZ' (50 GeV)", "l");
  legend->AddEntry(f70, "ttZ' (70 GeV, g_{b}=0.1)", "l");
  //legend->AddEntry(s100_flat, "ttZ' (100 GeV)", "l");
  //legend->AddEntry(s120_flat, "ttZ' (120 GeV)", "l");
  legend->AddEntry(f150, "ttZ' (150 GeV, g_{b}=0.2)", "l");



  c1->cd();
  //h1_flat->Draw("HIST E");
  //h1_flat->SetMaximum(20.);
  //h2_flat->Draw("SAME HIST E");
  bg_stacks->Draw("HIST E");
  bg_stacks->SetMaximum(35.);
  bg_stacks->GetXaxis()->SetRangeUser(20.,200.);
  bg_stacks->GetYaxis()->SetTitle("Events @3000 fb^{-1}");
  bg_stacks->GetXaxis()->SetTitle("M(#tau_{h}^{#pm} , l^{#mp})");
  //backgrounds_flat->Draw("HIST SAME");  ////// hide the added background th1f

  s40_flat->Draw("SAME HIST E");
  //s50_flat->Draw("SAME HIST E");
  s70_flat->Draw("SAME HIST E");
  //s100_flat->Draw("SAME HIST E");
  //s120_flat->Draw("SAME HIST E");
  s150_flat->Draw("SAME HIST E");
  //fBG->Draw("SAME C");
  //fTTW->Draw("SAME C");
  f40->Draw("SAME C");
  //f50->Draw("SAME C");
  f70->Draw("SAME C");
  //f100->Draw("SAME C");
  //f120->Draw("SAME C");
  f150->Draw("SAME C");
  legend->Draw();

  // Save out.
  string out_dir = "out/";
  string png = ".png";
  string pdf = ".pdf";
  string png_string = out_dir + outname +  png;
  string pdf_string = out_dir + outname +  pdf;
  c1->Print(png_string.c_str());
  c1->Print(pdf_string.c_str());


/*
    w->factory("CBShape::ttV(x,muV[60.5],sigmaV[16.3],alphaV[-0.99],nV[1.6])");
    w->factory("CBShape::ttW(x,muW[96.5],sigmaW[60.0],alphaW[-0.826],nW[1.75])");
*/


  // Testing
  TF1 *f1 = new TF1("f1","12.2*ROOT::Math::crystalball_function(x,-0.99,1.6,16.3,60.5)",0,200);
  TF1 *f2 = new TF1("f2","2.8*ROOT::Math::crystalball_function(x,-0.826,1.75,60.0,96.5)",0,200);
  c2->cd();
  f1->Draw("C");
  f2->Draw("SAME C");
  c2->Print("out/test.png");




}
