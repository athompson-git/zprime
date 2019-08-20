#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TEfficiency.h>
#include <TCanvas.h>
#include <TMath.h>
#include "lester_mt2_bisect.h"

#include <cmath>
#include <vector>
#include <string>
#include <utility>
#include <iostream>


// Set constants.
double kTauMass = 1.77682; // GeV
double kTopMass = 172.44;  // GeV
double kWMass = 80.385;  // GeV



unsigned int factorial(unsigned int n) {
  if (n == 0) return 1;
  return n * factorial(n - 1);
}


int sgn(Double_t number) {
  if (number < 0) return -1;
  return 1;
}


// Get the transverse mass.
double Mt(TLorentzVector *pt_miss, TLorentzVector *k) {
  return TMath::Sqrt(2 * (k->Pt()) * (pt_miss->Pt()) * (1 - TMath::Cos(k->DeltaPhi(*pt_miss))));
}


// Total transverse mass.
double GetTotalMT(TLorentzVector *p1, TLorentzVector *p2,
                  TLorentzVector *met) {
  double scalar = p1->Pt() + p2->Pt() + met->Pt();
  return sqrt(pow(scalar, 2) - ((*p1) + (*p2) + (*met)).M2());
}


// Returns the 4-vector magnitude p_zeta^miss - 0.85 * p_zeta^vis
double GetDZeta(TLorentzVector *vis1, TLorentzVector *vis2,
                TLorentzVector *pt_miss) {
  TVector3 pt_vis1(vis1->Px(), vis1->Py(), 0.);
  TVector3 pt_vis2(vis2->Px(), vis2->Py(), 0.);
  TVector3 pt_miss_v3(pt_miss->Px(), pt_miss->Py(), 0.);

  TVector3 zeta = (pt_vis2.Mag())*pt_vis1 + (pt_vis1.Mag())*pt_vis2;
  TVector3 pt_vis = pt_vis1 + pt_vis2;

  // Take the projections.
  Double_t pt_vis_zeta = pt_vis.Dot(zeta) * (1/zeta.Mag());
  Double_t pt_miss_zeta = pt_miss_v3.Dot(zeta) * (1/zeta.Mag());

  return pt_miss_zeta - 0.85 * pt_vis_zeta;
}



// Returns attempted reconstructed mass hypothesis.
// top = 1, w = 2
double MassHypothesis(TLorentzVector *tau, TLorentzVector *ell,
                      TLorentzVector *jet1, TLorentzVector *jet2,
                      int hyp) {
  double mass, hyp_mass;

  switch (hyp) {
    case 1: hyp_mass = kTopMass;
      break;
    case 2: hyp_mass = kWMass;
      break;
  }

  // Construct permutations.
  TLorentzVector p11 = *tau + *jet1;
  TLorentzVector p22 = *ell + *jet2;
  TLorentzVector p12 = *tau + *jet2;
  TLorentzVector p21 = *ell + *jet1;

  double m11 = p11.M();
  double m22 = p22.M();
  double m12 = p12.M();
  double m21 = p21.M();

  if ((m11 - m22) < (m12 - m21)) {
    mass = max(m11,m22);
  } else {
    mass = max(m12,m21);
  }

  return mass;
}


// MT2 calc.
double mt2(TLorentzVector *vis1, TLorentzVector *vis2, TLorentzVector *miss) {
  const double mVisA = 0.0; // Mass of visible object on side A.
  const double pxA = vis1->Px(); // x momentum of visible object on side A.
  const double pyA = vis1->Py(); // y momentum of visible object on side A.

  const double mVisB = 0.0; // Mass of visible object on side B.
  const double pxB = vis2->Px(); // x momentum of visible object on side B.
  const double pyB = vis2->Py(); // y momentum of visible object on side B.

  const double pxMiss = miss->Px(); // x component of missing transverse momentum.
  const double pyMiss = miss->Py(); // y component of missing transverse momentum.

  const double chiA = 0.0; // Hypothesised mass of invisible on side A.
  const double chiB = 0.0; // Hypothesised mass of invisible on side B.
  //printf("%f, %f, %f, %f, %f, %f, %f, %f", mVisA,pxA,pyA,mVisB,pxB,pyB,pxMiss,pyMiss);
  const double desiredPrecisionOnMt2 = 0.0; // Algo aims for machine precision.

  asymm_mt2_lester_bisect lester;
  return lester.get_mT2(mVisA, pxA, pyA, mVisB, pxB, pyB,
                        pxMiss, pyMiss, chiA, chiB,
                        desiredPrecisionOnMt2);
}


class Parton {
public:
  TLorentzVector p;
  Int_t q = 0;
  bool btag = false;
  bool tautag = false;
  bool muon = false;
  bool electron = false;

};


// Custom event class.
class Event : public TObject {
public:

  Parton tau;
  Parton l1;
  Parton b1;
  Parton met;

  ClassDef(Event,1)
};

ClassImp(Event)


////////////////////////////////////////////////////////////////////////////////


// Main macro.
void Analyzer(const char *sample_desc, const char *target_tree, int nbins, bool skip_histograms = false) {
  TFile *file_in = TFile::Open("/fdata/hepx/store/user/thompson/zprime_ditau/analysis/ntuples/analysis_tree.root");

  // Book histograms.
  TH1F *th1f_lepton_pt = new TH1F("lepton_pt", "", nbins, 0., 500.);
  TH1F *th1f_lepton_eta = new TH1F("lepton_eta", "", nbins, -2.5, 2.5);
  TH1F *th1f_lepton_phi = new TH1F("lepton_phi", "", nbins, 0.0, 3.2);
  TH1F *th1f_b_pt = new TH1F("b_pt", "", nbins, 0., 500.);
  TH1F *th1f_b2_pt = new TH1F("b2_pt", "", nbins, 0., 500.);
  TH1F *th1f_b_eta = new TH1F("b_eta", "", nbins, -2.5, 2.5);
  TH1F *th1f_b2_eta = new TH1F("b2_eta", "", nbins, -2.5, 2.5);
  TH1F *th1f_b_phi = new TH1F("b_phi", "", nbins, 0.0, 3.2);
  TH1F *th1f_b2_phi = new TH1F("b2_phi", "", nbins, 0.0, 3.2);
  TH1F *th1f_tau_h_pt = new TH1F("tauh_pt", "", nbins, 0., 500.);
  TH1F *th1f_tau_h_eta = new TH1F("tauh_eta", "", nbins, -2.5, 2.5);
  TH1F *th1f_tau_h_phi = new TH1F("tauh_phi", "", nbins, 0.0, 3.2);
  TH1F *th1f_deltaR_ditau = new TH1F("dR_ditau", "", nbins, 0.0, 6.0);
  TH1F *th1f_ditau_mass = new TH1F("ditau_mass", "", nbins, 0.0, 200.);
  th1f_ditau_mass->GetYaxis()->SetTitle("Events");
  TH1F *th1f_ditau_mass_topocut = new TH1F("ditau_mass_topocut", "", nbins, 0.0, 200.);
  th1f_ditau_mass_topocut->GetYaxis()->SetTitle("Events");
  TH1F *th1f_ditau_mass_fitter = new TH1F("ditau_mass_fitter", "", 200, 0.0, 200.);
  TH1F *th1f_ditau_mass_fitter_cut = new TH1F("ditau_mass_fitter_cut", "", 200, 0.0, 200.);
  TH1F *th1f_ditau_pt = new TH1F("ditau_pt", "", nbins, 0.0, 500.);
  TH1F *th1f_met = new TH1F("met", "", nbins, 0.0, 350.);
  TH1F *th1f_dzeta = new TH1F("dzeta", "", nbins, -300., 300.);
  TH1F *th1f_mt = new TH1F("mt", "", nbins, 0., 200.);
  TH1F *th1f_DPhi_btau = new TH1F("DPhi_btau", "", nbins, -3.5, 3.5);
  TH1F *th1f_DPhi_b2tau = new TH1F("DPhi_b2tau", "", nbins, -3.5, 3.5);
  TH1F *th1f_DPhi_bell = new TH1F("DPhi_bell", "", nbins, -3.5, 3.5);
  TH1F *th1f_DPhi_b2ell = new TH1F("DPhi_b2ell", "", nbins, -3.5, 3.5);
  TH1F *th1f_DPhi_b2b = new TH1F("DPhi_b2ell", "", nbins, -3.5, 3.5);
  TH1F *th1f_DPhi_elltau = new TH1F("DPhi_elltau", "", nbins, -3.5, 3.5);
  TH1F *th1f_DPhi_mettau = new TH1F("DPhi_mettau", "", nbins, -3.5, 3.5);
  TH1F *th1f_DPhi_metb = new TH1F("DPhi_metb", "", nbins, -3.5, 3.5);
  TH1F *th1f_DPhi_metb2 = new TH1F("DPhi_metb2", "", nbins, -3.5, 3.5);
  TH1F *th1f_DPhi_metell = new TH1F("DPhi_metell", "", nbins, -3.5, 3.5);
  TH1F *th1f_cosDPhi_btau = new TH1F("cosDPhi_btau", "", nbins, -1.2, 1.2);
  TH1F *th1f_cosDPhi_b2tau = new TH1F("cosDPhi_b2tau", "", nbins, -1.2, 1.2);
  TH1F *th1f_cosDPhi_bell = new TH1F("cosDPhi_bell", "", nbins, -1.2, 1.2);
  TH1F *th1f_cosDPhi_b2ell = new TH1F("cosDPhi_b2ell", "", nbins, -1.2, 1.2);
  TH1F *th1f_cosDPhi_b2b = new TH1F("cosDPhi_b2ell", "", nbins, -1.2, 1.2);
  TH1F *th1f_cosDPhi_elltau = new TH1F("cosDPhi_elltau", "", nbins, -1.2, 1.2);
  TH1F *th1f_cosDPhi_mettau = new TH1F("cosDPhi_mettau", "", nbins, -1.2, 1.2);
  TH1F *th1f_cosDPhi_metb = new TH1F("cosDPhi_metb", "", nbins, -1.2, 1.2);
  TH1F *th1f_cosDPhi_metb2 = new TH1F("cosDPhi_metb2", "", nbins, -1.2, 1.2);
  TH1F *th1f_cosDPhi_metell = new TH1F("cosDPhi_metell", "", nbins, -1.2, 1.2);
  TH1F *th1f_top_mass = new TH1F("top_mass", "", nbins, 0., 500.);
  TH1F *th1f_w_mass = new TH1F("w_mass", "", nbins, 0., 500.);
  TH1F *th1f_deltaR_b_tau = new TH1F("dR_b_tau", "", nbins, 0.0, 6.0);
  TH1F *th1f_deltaR_b_b2 = new TH1F("dR_b_b2", "", nbins, 0.0, 6.0);
  TH1F *th1f_deltaR_b2_tau = new TH1F("dR_b2_tau", "", nbins, 0.0, 6.0);
  TH1F *th1f_deltaR_b_lepton = new TH1F("dR_b_lep", "", nbins, 0.0, 6.0);
  TH1F *th1f_deltaR_b2_lepton = new TH1F("dR_b2_lep", "", nbins, 0.0, 6.0);
  TH1F *th1f_total_mt = new TH1F("total_mt", "", nbins, 0.0, 550.0);
  TH1F *th1f_total_mt_topocut = new TH1F("total_mt", "", nbins, 0.0, 550.0);
  th1f_total_mt->GetYaxis()->SetTitle("Events / 16 GeV");
  TH1F *th1f_equilibrant = new TH1F("vector_syst", "", nbins, 0., 500.);
  TH1F *th1f_MT2 = new TH1F("mt2", "", nbins, 0., 150.);


  // Other 2D histograms.
  TH2F *th2f_dzeta_M = new TH2F("dzeta_M", "", nbins, 0., 300., nbins, -300., 300.);
  TH2F *th2f_dR_M = new TH2F("dr_vs_M", "", nbins, 0., 300., nbins, 0., 6.);
  // Multiplicity histograms.
  TH1F *th1f_electron_multiplicity = new TH1F("el_mult", "", 5, -0.5, 4.5);
  TH1F *th1f_muon_multiplicity = new TH1F("mu_mult", "", 5, -0.5, 4.5);
  //TH1F *n_btags = new TH1F("bjets", "", 4, -.5, 3.5);


  // TTree infrastructure
  TTree *t2 = (TTree*)file_in->Get(target_tree);

  TBranch *b_tau = t2->GetBranch("TauBranch");
  TBranch *b_lep1 = t2->GetBranch("Lep1Branch");
  TBranch *b_lep2 = t2->GetBranch("Lep2Branch");
  TBranch *b_btag = t2->GetBranch("BTagBranch");
  TBranch *b_jet = t2->GetBranch("JetBranch");
  TBranch *b_met = t2->GetBranch("METBranch");
  TBranch *b_wgt = t2->GetBranch("Weight");
  //TBranch *b_nb = t2->GetBranch("Nb");

  Float_t tau_arr[4];
  Float_t lep1_arr[4];
  Float_t lep2_arr[4];
  Float_t btag_arr[4];
  Float_t jet_arr[4];
  Float_t met_arr[4];
  Float_t reweight;
  int nb;

  b_tau->SetAddress(&tau_arr);
  b_lep1->SetAddress(&lep1_arr);
  b_lep2->SetAddress(&lep2_arr);
  b_btag->SetAddress(&btag_arr);
  b_jet->SetAddress(&jet_arr);
  b_met->SetAddress(&met_arr);
  b_wgt->SetAddress(&reweight);
  //b_nb->SetAddress(&nb);


  // Declare cutflow variable efficiencies.
  float delta_r_eff = 0;
  float cdphi_metell_eff = 0;
  float cdphi_mettau_eff = 0;



  // EVENT LOOP.
  Long64_t nentries = t2->GetEntries();
  for (Long64_t i=0; i<nentries; i++) {
    if (i % 50 == 0) printf("On event %lld / %lld \n", i, nentries);


    b_tau->GetEntry(i);
    b_lep1->GetEntry(i);
    b_lep2->GetEntry(i);
    b_btag->GetEntry(i);
    b_jet->GetEntry(i);
    b_met->GetEntry(i);
    b_wgt->GetEntry(i);
    //b_nb->GetEntry(i);

    // Make TLorentzVectors
    TLorentzVector tau_h_p4, lepton_p4, lepton2_p4, met_p4, b_p4, b2_p4;
    tau_h_p4.SetPtEtaPhiE(tau_arr[0], tau_arr[1], tau_arr[2], tau_arr[3]);
    lepton_p4.SetPtEtaPhiE(lep1_arr[0], lep1_arr[1], lep1_arr[2], lep1_arr[3]);
    lepton2_p4.SetPtEtaPhiE(lep2_arr[0], lep2_arr[1], lep2_arr[2], lep2_arr[3]);
    b_p4.SetPtEtaPhiE(btag_arr[0], btag_arr[1], btag_arr[2], btag_arr[3]);
    b2_p4.SetPtEtaPhiE(jet_arr[0], jet_arr[1], jet_arr[2], jet_arr[3]);
    met_p4.SetPtEtaPhiE(met_arr[0], met_arr[1], met_arr[2], met_arr[3]);

    TLorentzVector ditau = tau_h_p4 + lepton_p4;
    TLorentzVector boosted_syst = ditau + met_p4;
    TLorentzVector recoil_syst = b_p4 + b2_p4;
    TLorentzVector balanced_syst = boosted_syst + recoil_syst;

    // KINEMATICS //////////////////////////////////////////////////////////////

    // Topological Cuts.
    if (cos(lepton_p4.DeltaPhi(met_p4)) > 0.0) {
      cdphi_metell_eff += 1;
      if (tau_h_p4.DeltaR(lepton_p4) < 2.) {
        delta_r_eff += 1;
      }
    }

    // High-level.
    th1f_MT2->Fill(mt2(&tau_h_p4, &lepton_p4, &met_p4));
    th1f_top_mass->Fill(MassHypothesis(&tau_h_p4, &lepton_p4, &b_p4,
                                       &b2_p4, 1), reweight*3000000);
    th1f_w_mass->Fill(MassHypothesis(&tau_h_p4, &lepton_p4, &b_p4,
                                     &b2_p4, 2));
    th1f_dzeta->Fill(GetDZeta(&tau_h_p4, &lepton_p4, &met_p4));
    th1f_mt->Fill(Mt(&met_p4, &lepton_p4));
    th2f_dR_M->Fill(ditau.M(), tau_h_p4.DeltaR(lepton_p4));
    th2f_dzeta_M->Fill(ditau.M(), GetDZeta(&tau_h_p4, &lepton_p4, &met_p4));
    th1f_equilibrant->Fill(balanced_syst.Pt());

    // topology
    th1f_cosDPhi_elltau->Fill(cos(lepton_p4.DeltaPhi(tau_h_p4)));
    th1f_cosDPhi_btau->Fill(cos(b_p4.DeltaPhi(tau_h_p4)));
    th1f_cosDPhi_bell->Fill(cos(lepton_p4.DeltaPhi(b_p4)));
    th1f_cosDPhi_b2ell->Fill(cos(b2_p4.DeltaPhi(lepton_p4)));
    th1f_cosDPhi_b2tau->Fill(cos(b2_p4.DeltaPhi(tau_h_p4)));
    th1f_cosDPhi_b2b->Fill(cos(b_p4.DeltaPhi(b2_p4)));
    th1f_cosDPhi_mettau->Fill(cos(met_p4.DeltaPhi(tau_h_p4)));
    th1f_cosDPhi_metb->Fill(cos(met_p4.DeltaPhi(b_p4)));
    th1f_cosDPhi_metb2->Fill(cos(met_p4.DeltaPhi(b2_p4)));
    th1f_cosDPhi_metell->Fill(cos(met_p4.DeltaPhi(lepton_p4)), reweight*3000000);
    th1f_DPhi_elltau->Fill(lepton_p4.DeltaPhi(tau_h_p4));
    th1f_DPhi_btau->Fill(b_p4.DeltaPhi(tau_h_p4));
    th1f_DPhi_bell->Fill(lepton_p4.DeltaPhi(b_p4));
    th1f_DPhi_b2ell->Fill(b2_p4.DeltaPhi(lepton_p4));
    th1f_DPhi_b2b->Fill(b_p4.DeltaPhi(b2_p4));
    th1f_DPhi_b2tau->Fill(b2_p4.DeltaPhi(tau_h_p4));
    th1f_DPhi_mettau->Fill(met_p4.DeltaPhi(tau_h_p4));
    th1f_DPhi_metb->Fill(met_p4.DeltaPhi(b_p4));
    th1f_DPhi_metb2->Fill(met_p4.DeltaPhi(b2_p4));
    th1f_DPhi_metell->Fill(met_p4.DeltaPhi(lepton_p4));

    // pT plots.
    th1f_b_pt->Fill(b_p4.Pt());
    th1f_b2_pt->Fill(b2_p4.Pt());
    th1f_tau_h_pt->Fill(tau_h_p4.Pt());
    th1f_lepton_pt->Fill(lepton_p4.Pt());
    th1f_met->Fill(met_p4.Pt());

    // Eta plots.
    th1f_b_eta->Fill(b_p4.Eta());
    th1f_b2_eta->Fill(b_p4.Eta());
    th1f_tau_h_eta->Fill(tau_h_p4.Eta());
    th1f_lepton_eta->Fill(lepton_p4.Eta());

    // Phi plots.
    th1f_b_phi->Fill(b_p4.Phi());
    th1f_b2_phi->Fill(b2_p4.Phi());
    th1f_tau_h_phi->Fill(tau_h_p4.Phi());
    th1f_lepton_phi->Fill(lepton_p4.Phi());

    // Delta R plots.
    th1f_deltaR_ditau->Fill(tau_h_p4.DeltaR(lepton_p4), reweight*3000000);
    th1f_deltaR_b_tau->Fill(tau_h_p4.DeltaR(b_p4));
    th1f_deltaR_b_b2->Fill(b_p4.DeltaR(b2_p4));
    th1f_deltaR_b2_tau->Fill(tau_h_p4.DeltaR(b2_p4));
    th1f_deltaR_b_lepton->Fill(lepton_p4.DeltaR(b_p4));
    th1f_deltaR_b2_lepton->Fill(lepton_p4.DeltaR(b2_p4));


    // Invariant masses.
    TLorentzVector total_mt = tau_h_p4 + lepton_p4 + met_p4;
    th1f_ditau_mass->Fill(ditau.M(), reweight*3000000);
    th1f_ditau_mass_fitter->Fill(ditau.M(), reweight*3000000);
    th1f_total_mt->Fill(total_mt.Mt(), reweight*3000000);
    if (cos(lepton_p4.DeltaPhi(met_p4)) > 0.0) {
      th1f_ditau_mass_topocut->Fill(ditau.M(), reweight*3000000);
      th1f_ditau_mass_fitter_cut->Fill(ditau.M(), reweight*3000000);
      th1f_total_mt_topocut->Fill(total_mt.Mt(), reweight*3000000);
    }
    th1f_ditau_pt->Fill(ditau.Pt());

    //n_btags->Fill(nb);


  } // End event loop.

  printf("delta_r_eff = %f \n", delta_r_eff / nentries);
  printf("cdphi_met_ell_eff = %f \n", cdphi_metell_eff / nentries);
  printf("cdphi_met_tau_eff = %f \n", cdphi_mettau_eff / nentries);


  gSystem->cd("/fdata/hepx/store/user/thompson/zprime_ditau/analysis/histograms/semileptonic");

  th1f_lepton_pt->Scale(1/th1f_lepton_pt->Integral());
  TFile *f1 = new TFile("hist_lepton_pt.root", "UPDATE");
  th1f_lepton_pt->Write(sample_desc, TObject::kOverwrite);

  th1f_lepton_eta->Scale(1/th1f_lepton_eta->Integral());
  TFile *f2 = new TFile("hist_lepton_eta.root", "UPDATE");
  th1f_lepton_eta->Write(sample_desc, TObject::kOverwrite);

  th1f_lepton_phi->Scale(1/th1f_lepton_phi->Integral());
  TFile *f3 = new TFile("hist_lepton_phi.root", "UPDATE");
  th1f_lepton_phi->Write(sample_desc, TObject::kOverwrite);

  th1f_b_pt->Scale(1/th1f_b_pt->Integral());
  TFile *f4 = new TFile("hist_b_pt.root", "UPDATE");
  th1f_b_pt->Write(sample_desc, TObject::kOverwrite);

  th1f_b_eta->Scale(1/th1f_b_eta->Integral());
  TFile *f5 = new TFile("hist_b_eta.root", "UPDATE");
  th1f_b_eta->Write(sample_desc, TObject::kOverwrite);

  th1f_b_phi->Scale(1/th1f_b_phi->Integral());
  TFile *f6 = new TFile("hist_b_phi.root", "UPDATE");
  th1f_b_phi->Write(sample_desc, TObject::kOverwrite);

  th1f_tau_h_pt->Scale(1/th1f_tau_h_pt->Integral());
  TFile *f16 = new TFile("hist_tau_h_pt.root", "UPDATE");
  th1f_tau_h_pt->Write(sample_desc, TObject::kOverwrite);

  th1f_tau_h_eta->Scale(1/th1f_tau_h_eta->Integral());
  TFile *f17 = new TFile("hist_tau_h_eta.root", "UPDATE");
  th1f_tau_h_eta->Write(sample_desc, TObject::kOverwrite);

  th1f_tau_h_phi->Scale(1/th1f_tau_h_phi->Integral());
  TFile *f18 = new TFile("hist_tau_h_phi.root", "UPDATE");
  th1f_tau_h_phi->Write(sample_desc, TObject::kOverwrite);

  th1f_deltaR_ditau->Scale(1/th1f_deltaR_ditau->Integral());
  TFile *f19 = new TFile("hist_deltaR_ditau.root", "UPDATE");
  th1f_deltaR_ditau->Write(sample_desc, TObject::kOverwrite);

  th1f_met->Scale(1/th1f_met->Integral());
  TFile *f24 = new TFile("hist_met.root", "UPDATE");
  th1f_met->Write(sample_desc, TObject::kOverwrite);

  th1f_ditau_pt->Scale(1/th1f_ditau_pt->Integral());
  TFile *f29 = new TFile("hist_ditau_pt.root", "UPDATE");
  th1f_ditau_pt->Write(sample_desc, TObject::kOverwrite);

  th1f_dzeta->Scale(1/th1f_dzeta->Integral());
  TFile *f30 = new TFile("hist_dzeta.root", "UPDATE");
  th1f_dzeta->Write(sample_desc, TObject::kOverwrite);

  th1f_mt->Scale(1/th1f_mt->Integral());
  TFile *f31 = new TFile("hist_mt.root", "UPDATE");
  th1f_mt->Write(sample_desc, TObject::kOverwrite);

  th1f_cosDPhi_btau->Scale(1/th1f_cosDPhi_btau->Integral());
  TFile *f34 = new TFile("hist_cosDPhi_btau.root", "UPDATE");
  th1f_cosDPhi_btau->Write(sample_desc, TObject::kOverwrite);

  th1f_cosDPhi_elltau->Scale(1/th1f_cosDPhi_elltau->Integral());
  TFile *f35 = new TFile("hist_cosDPhi_elltau.root", "UPDATE");
  th1f_cosDPhi_elltau->Write(sample_desc, TObject::kOverwrite);

  TFile *f36 = new TFile("hist_ditau_mass.root", "UPDATE");
  th1f_ditau_mass->Write(sample_desc, TObject::kOverwrite);

  TFile *f37 = new TFile("hist_ditau_mass_topocut.root", "UPDATE");
  th1f_ditau_mass_topocut->Write(sample_desc, TObject::kOverwrite);

  TFile *f39 = new TFile("hist_top_mass.root", "UPDATE");
  th1f_top_mass->Write(sample_desc, TObject::kOverwrite);

  th1f_w_mass->Scale(1/th1f_w_mass->Integral());
  TFile *f40 = new TFile("hist_w_mass.root", "UPDATE");
  th1f_w_mass->Write(sample_desc, TObject::kOverwrite);

  th1f_MT2->Scale(1/th1f_MT2->Integral());
  TFile *f41 = new TFile("hist_MT2.root", "UPDATE");
  th1f_MT2->Write(sample_desc, TObject::kOverwrite);

  th1f_muon_multiplicity->Scale(1/th1f_muon_multiplicity->Integral());
  TFile *f42 = new TFile("hist_mu_mult.root", "UPDATE");
  th1f_muon_multiplicity->Write(sample_desc, TObject::kOverwrite);

  th1f_electron_multiplicity->Scale(1/th1f_electron_multiplicity->Integral());
  TFile *f43 = new TFile("hist_e_mult.root", "UPDATE");
  th1f_electron_multiplicity->Write(sample_desc, TObject::kOverwrite);

  th1f_deltaR_b_tau->Scale(1/th1f_deltaR_b_tau->Integral());
  TFile *f44 = new TFile("hist_dR_b_tau.root", "UPDATE");
  th1f_deltaR_b_tau->Write(sample_desc, TObject::kOverwrite);
  th1f_deltaR_b_b2->Scale(1/th1f_deltaR_b_b2->Integral());
  TFile *f45 = new TFile("hist_dR_b_b2.root", "UPDATE");
  th1f_deltaR_b_b2->Write(sample_desc, TObject::kOverwrite);
  th1f_deltaR_b2_tau->Scale(1/th1f_deltaR_b2_tau->Integral());
  TFile *f46 = new TFile("hist_dR_b2_tau.root", "UPDATE");
  th1f_deltaR_b2_tau->Write(sample_desc, TObject::kOverwrite);
  th1f_deltaR_b_lepton->Scale(1/th1f_deltaR_b_lepton->Integral());
  TFile *f47 = new TFile("hist_dR_b_lepton.root", "UPDATE");
  th1f_deltaR_b_lepton->Write(sample_desc, TObject::kOverwrite);
  th1f_deltaR_b2_lepton->Scale(1/th1f_deltaR_b2_lepton->Integral());
  TFile *f48 = new TFile("hist_dR_b2_lepton.root", "UPDATE");
  th1f_deltaR_b2_lepton->Write(sample_desc, TObject::kOverwrite);

  th1f_b2_pt->Scale(1/th1f_b2_pt->Integral());
  TFile *f49 = new TFile("hist_b2_pt.root", "UPDATE");
  th1f_b2_pt->Write(sample_desc, TObject::kOverwrite);

  th1f_b2_eta->Scale(1/th1f_b2_eta->Integral());
  TFile *f50 = new TFile("hist_b2_eta.root", "UPDATE");
  th1f_b2_eta->Write(sample_desc, TObject::kOverwrite);

  th1f_b2_phi->Scale(1/th1f_b2_phi->Integral());
  TFile *f51 = new TFile("hist_b2_phi.root", "UPDATE");
  th1f_b2_phi->Write(sample_desc, TObject::kOverwrite);


  th1f_cosDPhi_b2tau->Scale(1/th1f_cosDPhi_b2tau->Integral());
  TFile *f52 = new TFile("hist_cosDPhi_b2tau.root", "UPDATE");
  th1f_cosDPhi_b2tau->Write(sample_desc, TObject::kOverwrite);

  th1f_cosDPhi_bell->Scale(1/th1f_cosDPhi_bell->Integral());
  TFile *f53 = new TFile("hist_cosDPhi_bell.root", "UPDATE");
  th1f_cosDPhi_bell->Write(sample_desc, TObject::kOverwrite);

  th1f_cosDPhi_b2ell->Scale(1/th1f_cosDPhi_b2ell->Integral());
  TFile *f54 = new TFile("hist_cosDPhi_b2ell.root", "UPDATE");
  th1f_cosDPhi_b2ell->Write(sample_desc, TObject::kOverwrite);

  th1f_cosDPhi_mettau->Scale(1/th1f_cosDPhi_mettau->Integral());
  TFile *f55 = new TFile("hist_cosDPhi_mettau.root", "UPDATE");
  th1f_cosDPhi_mettau->Write(sample_desc, TObject::kOverwrite);

  th1f_cosDPhi_metb->Scale(1/th1f_cosDPhi_metb->Integral());
  TFile *f56 = new TFile("hist_cosDPhi_metb.root", "UPDATE");
  th1f_cosDPhi_metb->Write(sample_desc, TObject::kOverwrite);

  th1f_cosDPhi_metb2->Scale(1/th1f_cosDPhi_metb2->Integral());
  TFile *f57 = new TFile("hist_cosDPhi_metb2.root", "UPDATE");
  th1f_cosDPhi_metb2->Write(sample_desc, TObject::kOverwrite);

  th1f_cosDPhi_metell->Scale(1/th1f_cosDPhi_metell->Integral());
  TFile *f58 = new TFile("hist_cosDPhi_metell.root", "UPDATE");
  th1f_cosDPhi_metell->Write(sample_desc, TObject::kOverwrite);

  th1f_cosDPhi_b2b->Scale(1/th1f_cosDPhi_b2b->Integral());
  TFile *f59 = new TFile("hist_cosDPhi_b2b.root", "UPDATE");
  th1f_cosDPhi_b2b->Write(sample_desc, TObject::kOverwrite);

// dPhi
  th1f_DPhi_b2tau->Scale(1/th1f_DPhi_b2tau->Integral());
  TFile *f60 = new TFile("hist_DPhi_b2tau.root", "UPDATE");
  th1f_DPhi_b2tau->Write(sample_desc, TObject::kOverwrite);

  th1f_DPhi_bell->Scale(1/th1f_DPhi_bell->Integral());
  TFile *f61 = new TFile("hist_DPhi_bell.root", "UPDATE");
  th1f_DPhi_bell->Write(sample_desc, TObject::kOverwrite);

  th1f_DPhi_b2ell->Scale(1/th1f_DPhi_b2ell->Integral());
  TFile *f62 = new TFile("hist_DPhi_b2ell.root", "UPDATE");
  th1f_DPhi_b2ell->Write(sample_desc, TObject::kOverwrite);

  th1f_DPhi_mettau->Scale(1/th1f_DPhi_mettau->Integral());
  TFile *f63 = new TFile("hist_DPhi_mettau.root", "UPDATE");
  th1f_DPhi_mettau->Write(sample_desc, TObject::kOverwrite);

  th1f_DPhi_metb->Scale(1/th1f_DPhi_metb->Integral());
  TFile *f64 = new TFile("hist_DPhi_metb.root", "UPDATE");
  th1f_DPhi_metb->Write(sample_desc, TObject::kOverwrite);

  th1f_DPhi_metb2->Scale(1/th1f_DPhi_metb2->Integral());
  TFile *f65 = new TFile("hist_DPhi_metb2.root", "UPDATE");
  th1f_DPhi_metb2->Write(sample_desc, TObject::kOverwrite);

  th1f_DPhi_metell->Scale(1/th1f_DPhi_metell->Integral());
  TFile *f66 = new TFile("hist_DPhi_metell.root", "UPDATE");
  th1f_DPhi_metell->Write(sample_desc, TObject::kOverwrite);

  th1f_DPhi_btau->Scale(1/th1f_DPhi_btau->Integral());
  TFile *f67 = new TFile("hist_DPhi_btau.root", "UPDATE");
  th1f_DPhi_btau->Write(sample_desc, TObject::kOverwrite);

  th1f_DPhi_elltau->Scale(1/th1f_DPhi_elltau->Integral());
  TFile *f68 = new TFile("hist_DPhi_elltau.root", "UPDATE");
  th1f_DPhi_elltau->Write(sample_desc, TObject::kOverwrite);

  th1f_DPhi_b2b->Scale(1/th1f_DPhi_b2b->Integral());
  TFile *f69 = new TFile("hist_DPhi_b2b.root", "UPDATE");
  th1f_DPhi_b2b->Write(sample_desc, TObject::kOverwrite);

  TFile *f70 = new TFile("hist_total_mt.root", "UPDATE");
  th1f_total_mt->Write(sample_desc, TObject::kOverwrite);

  th1f_equilibrant->Scale(1/th1f_equilibrant->Integral());
  TFile *f72 = new TFile("hist_equilibrant.root", "UPDATE");
  th1f_equilibrant->Write(sample_desc, TObject::kOverwrite);

  TFile *f74 = new TFile("hist_total_mt_topocut.root", "UPDATE");
  th1f_total_mt_topocut->Write(sample_desc, TObject::kOverwrite);


  // Save TH2Fs

  th2f_dR_M->Scale(1/th2f_dR_M->Integral());
  TFile *f2D_5 = new TFile("th2f_dR_M.root", "UPDATE");
  th2f_dR_M->Write(sample_desc, TObject::kOverwrite);

  th2f_dzeta_M->Scale(1/th2f_dzeta_M->Integral());
  TFile *f2D_6 = new TFile("th2f_dzeta_M.root", "UPDATE");
  th2f_dzeta_M->Write(sample_desc, TObject::kOverwrite);


  // Save special TH1F mass spectra for Fitter.C
  gSystem->cd("/fdata/hepx/store/user/thompson/zprime_ditau/analysis/histograms/experimental");
  TFile *s1 = new TFile("hist_ditau_mass.root", "UPDATE");
  th1f_ditau_mass_fitter->Write(sample_desc, TObject::kOverwrite);

  TFile *s2 = new TFile("hist_ditau_mass_topocut.root", "UPDATE");
  th1f_ditau_mass_fitter_cut->Write(sample_desc, TObject::kOverwrite);

  gSystem->cd("/fdata/hepx/store/user/thompson/zprime_ditau/analysis/src");


}  // End macro.


