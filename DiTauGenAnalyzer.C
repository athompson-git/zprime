#ifdef __CLING__
#include "/Users/adrianthompson/MG5_aMC_v2_6_1/ExRootAnalysis/ExRootAnalysis/ExRootTreeReader.h"
#include "/Users/adrianthompson/MG5_aMC_v2_6_1/ExRootAnalysis/ExRootAnalysis/ExRootResult.h"
#include "/Users/adrianthompson/MG5_aMC_v2_6_1/ExRootAnalysis/ExRootAnalysis/ExRootClasses.h"
#include "lester_mt2_bisect.h"
#include <utility> // std::pair, std::make_pair
#else
class ExRootTreeReader;
class ExRootResult;
#endif

// Usage:
// gSystem->Load("~/MG5_aMC_v2_6_1/ExRootAnalysis/libExRootAnalysis");
// .x DiTauGenAnalyzer.C("zp_ditau_sample01/Events/350GeV_80kEvents/350GeV_80k_unweighted_events.root", "Z' (350 GeV)", 40, true);

void DiTauGenAnalyzer(const char *file_name, const char *sample_desc,
                      int nbins, bool apply_cuts = false) {

  TChain chain("LHEF");
  chain.Add(file_name);

  // Create object of class ExRootTreeReader.
  ExRootTreeReader *tree_reader = new ExRootTreeReader(&chain);
  Long64_t number_of_entries = tree_reader->GetEntries();

  // Get pointers to branches used in this analysis.
  TClonesArray *branch_particle = tree_reader->UseBranch("Particle");

  // Book histograms.
  TH1F *hist_pair_mass = new TH1F("pair_mass", "Max{M(tau,j)}", nbins, 0., 500.);
  TH1F *hist_MET = new TH1F("met", "Normalized Missing ET", nbins, 0., 2.);
  TH1F *hist_MT2 = new TH1F("tmass", "MT2 (Ditau + MET)", nbins, 0., 150.);
  TH1F *hist_unboosted_MT2 = new TH1F("unboost_mt2", "", nbins, 0., 150.);
  TH1F *hist_HT_LT = new TH1F("HT_LT", "HT - LT", nbins, -500., 500.);
  TH1F *hist_ditau_mass = new TH1F("ditau_mass", "M(tau+,tau-)", nbins, 0., 500.);
  TH1F *hist_pt_b = new TH1F("pt_btag", "BTag Pt", nbins, 0., 300.);
  TH1F *hist_pt_sub_b = new TH1F("pt_jet", "Secondary Jet Pt", nbins, 0., 300.);
  TH1F *hist_pt_tau_p = new TH1F("pt_tau_p", "Tau+ Pt", nbins, 0., 300.);
  TH1F *hist_pt_tau_m = new TH1F("pt_tau_m", "Tau- Pt", nbins, 0., 300.);
  TH1F *hist_pt_ditau = new TH1F("pt_ditau", "Ditau Pt", nbins, 0., 300.);
  TH1F *hist_pt_nu = new TH1F("pt_nu", "#nu_{#tau} Pt", nbins, 0., 300.);
  TH1F *hist_pt_antinu = new TH1F("pt_anti_nu", "anti - #nu_{#tau} Pt", nbins, 0., 300.);
  TH1F *hist_e_b = new TH1F("e_btag", "BTag Jet E", nbins, 0., 300.);
  TH1F *hist_e_sub_b = new TH1F("e_jet", "Secondary Jet E", nbins, 0., 300.);
  TH1F *hist_e_tau_p = new TH1F("e_tau_p", "E(#tau^{+}", nbins, 0., 300.);
  TH1F *hist_e_tau_m = new TH1F("e_tau_m", "E(#tau^{-})", nbins, 0., 300.);
  TH1F *hist_denis = new TH1F("hist_denis", "", nbins, -.3, .3);
  TH1F *hist_topology = new TH1F("hist_unboosted_topo", "", nbins, -3.14, 3.14);

  // Main Event loop.
  for (Int_t entry = 0; entry < number_of_entries; ++entry) {
    tree_reader->ReadEntry(entry);

    // Declare physics objects.
    TRootLHEFParticle *particle;
    TLorentzVector tlv_tau_plus;
    TLorentzVector tlv_nu;
    TLorentzVector tlv_tau_minus;
    TLorentzVector tlv_anti_nu;
    TLorentzVector tlv_b;
    TLorentzVector tlv_sub_b;

    Int_t particle_size = branch_particle->GetEntries();

    bool found_tau_plus = false;
    bool found_tau_minus = false;
    bool found_b = false;
    bool found_sub_b = false;
    bool found_nu = false;
    bool found_anti_nu = false;
    Int_t leading_b_id = -1;

    // Loop over to find OS Tau's.
    for (Int_t i = 0; i < particle_size; i++) {
      particle = (TRootLHEFParticle*) branch_particle->At(i);
      printf("%d th entry, particle ID = %f \n", i, particle->PT);
      // ID Tau (-).
      if (particle->PID == 15 && particle->PT > tlv_tau_minus.Pt()) {
        tlv_tau_minus.SetPtEtaPhiM(particle->PT, particle->Eta, particle->Phi,
                                   particle->M);
        printf("tau minus PID = %d, PT = %f \n", particle->PID, particle->PT);
        found_tau_minus = true;
      }

      // ID anti-Tau-neutrino, if it exists.
      if (particle->PID == -16 && particle->PT > tlv_anti_nu.Pt()) {
        tlv_anti_nu.SetPtEtaPhiM(particle->PT, particle->Eta, particle->Phi,
                                 particle->M);
        found_anti_nu = true;
      }

      // ID Tau (+).
      if (particle->PID == -15 && particle->PT > tlv_tau_plus.Pt()) {
        tlv_tau_plus.SetPtEtaPhiM(particle->PT, particle->Eta, particle->Phi,
                                  particle->M);
        found_tau_plus = true;
      }

      // ID Tau-neutrino, if it exists.
      if (particle->PID == 16 && particle->PT > tlv_nu.Pt()) {
        tlv_nu.SetPtEtaPhiM(particle->PT, particle->Eta, particle->Phi,
                            particle->M);
        found_nu = true;
      }

      // ID b quark.
      if (abs(particle->PID) == 5 && particle->PT > tlv_b.Pt()) {
        tlv_b.SetPtEtaPhiM(particle->PT, particle->Eta, particle->Phi,
                           particle->M);
        found_b = true;
        leading_b_id = i;
      }
    }

    for (Int_t i = 0; i < particle_size; i++) {
      particle = (TRootLHEFParticle*) branch_particle->At(i);
      if (abs(particle->PID) == 5 && leading_b_id != i) {
        tlv_sub_b.SetPtEtaPhiM(particle->PT, particle->Eta, particle->Phi,
                               particle->M);
        found_sub_b = true;
      }
    }

    // Make selection cuts.
    if (!found_b || !found_tau_plus || !found_tau_minus) continue;
    if (apply_cuts) {
      if (tlv_b.Pt() < 30.) continue;
      if (tlv_sub_b.Pt() < 30.) continue;
      if (tlv_tau_plus.Pt() < 70.) continue;
      if (tlv_tau_minus.Pt() < 70.) continue;
    }

    // Build composite physics objects.
    TLorentzVector tau_p_b = tlv_tau_plus + tlv_b;
    TLorentzVector tau_m_j = tlv_tau_minus + tlv_sub_b;
    TLorentzVector tau_p_j = tlv_tau_plus + tlv_sub_b;
    TLorentzVector tau_m_b = tlv_tau_minus + tlv_b;
    TLorentzVector ditau = tlv_tau_plus + tlv_tau_minus;
    TLorentzVector tlv_met = tlv_nu + tlv_anti_nu;
    std::pair <TLorentzVector, TLorentzVector> b_tau_pair; // Matched pair.

    TLorentzVector jet_recoil = tlv_b + tlv_sub_b;
    TLorentzVector unboost_tau_plus = tlv_tau_plus + jet_recoil;
    TLorentzVector unboost_tau_minus = tlv_tau_minus + jet_recoil;
    TLorentzVector unboost_met = tlv_met + jet_recoil;

    // Fill histograms.

    // Find the right tau-jet pairing by taking the permutation of tau-jet
    // pairs (where jet = b or non-b) with the smallest mass difference, then
    // picking the highest pair mass of that permutation.
    Double_t choice_pair_mass;
    if (abs(tau_p_b.M() - tau_m_j.M()) < abs(tau_p_j.M() - tau_m_b.M())) {
      b_tau_pair = std::make_pair(tlv_b, tlv_tau_plus);
      choice_pair_mass = std::max(tau_p_b.M(), tau_m_j.M());
    } else {
      b_tau_pair = std::make_pair(tlv_b, tlv_tau_minus);
      choice_pair_mass = std::max(tau_p_j.M(), tau_m_b.M());
    }
    hist_pair_mass->Fill(choice_pair_mass);

    // MET
    hist_MET->Fill(tlv_met.Pt());

    // HT - LT
    Double_t H_T = tlv_b.Pt() + tlv_sub_b.Pt();
    Double_t L_T = tlv_tau_plus.Pt() + tlv_tau_minus.Pt();
    hist_HT_LT->Fill(H_T - L_T);

    // Ditau pair mass.
    hist_ditau_mass->Fill(ditau.M());

    // Calculate MT2 (http://www.hep.phy.cam.ac.uk/~lester/mt2/).
    Double_t mVisA = tlv_tau_plus.M(); // Mass of visible object on side A.
    Double_t pxA = tlv_tau_plus.Px(); // x momentum of visible object on side A.
    Double_t pyA = tlv_tau_plus.Py(); // y momentum of visible object on side A.

    Double_t mVisB = tlv_tau_minus.M(); // Mass of visible object on side B.
    Double_t pxB = tlv_tau_minus.Px(); // x momentum of visible object on side B.
    Double_t pyB =tlv_tau_minus.Py(); // y momentum of visible object on side B.

    Double_t pxMiss = tlv_met.Px(); // x component of missing transverse momentu$
    Double_t pyMiss = tlv_met.Py(); // y component of missing transverse momentu$

    Double_t chiA = 0.0; // Hypothesised mass of invisible on side A.
    Double_t chiB = 0.0; // Hypothesised mass of invisible on side B.

    Double_t desiredPrecisionOnMt2 = 0; // Algo aims for machine precision.

    Double_t MT2 = asymm_mt2_lester_bisect::get_mT2(
                   mVisA, pxA, pyA,
                   mVisB, pxB, pyB,
                   pxMiss, pyMiss,
                   chiA, chiB,
                   desiredPrecisionOnMt2);
    hist_MT2->Fill(MT2);

    // Calculate UNBOOSTED MT2.
    Double_t unboost_pxA = unboost_tau_plus.Px();
    Double_t unboost_pyA = unboost_tau_plus.Py();
    Double_t unboost_pxB = unboost_tau_minus.Px();
    Double_t unboost_pyB = unboost_tau_minus.Py();
    Double_t unboost_pxMiss = unboost_met.Px();
    Double_t unboost_pyMiss = unboost_met.Py();

    Double_t unboost_MT2 = asymm_mt2_lester_bisect::get_mT2(
                           mVisA, unboost_pxA, unboost_pyA,
                           mVisB, unboost_pxB, unboost_pyB,
                           unboost_pxMiss, unboost_pyMiss, chiA,
                           chiB, desiredPrecisionOnMt2);
    hist_unboosted_MT2->Fill(unboost_MT2);

    // Basic kinematic variables.

    hist_pt_tau_p->Fill(tlv_tau_plus.Pt());
    hist_pt_tau_m->Fill(tlv_tau_minus.Pt());
    hist_pt_nu->Fill(tlv_nu.Pt());
    hist_pt_antinu->Fill(tlv_anti_nu.Pt());
    hist_pt_b->Fill(tlv_b.Pt());
    hist_pt_sub_b->Fill(tlv_sub_b.Pt());


  } // End event loop.

  // Normalize and save histograms.
  hist_pair_mass->Scale(1/hist_pair_mass->Integral());
  TFile *f1 = new TFile("ghist_taujet_mass.root", "UPDATE");
  hist_pair_mass->Write(sample_desc, TObject::kOverwrite);

  hist_HT_LT->Scale(1/hist_HT_LT->Integral());
  TFile *f2 = new TFile("ghist_ht_lt.root", "UPDATE");
  hist_HT_LT->Write(sample_desc, TObject::kOverwrite);

  hist_MET->Scale(1/hist_MET->Integral());
  TFile *f3 = new TFile("ghist_met.root", "UPDATE");
  hist_MET->Write(sample_desc, TObject::kOverwrite);

  hist_ditau_mass->Scale(1/hist_ditau_mass->Integral());
  TFile *f4 = new TFile("ghist_ditau_mass.root", "UPDATE");
  hist_ditau_mass->Write(sample_desc, TObject::kOverwrite);

  hist_MT2->Scale(1/hist_MT2->Integral());
  TFile *f5 = new TFile("ghist_mt2.root", "UPDATE");
  hist_MT2->Write(sample_desc, TObject::kOverwrite);

  hist_unboosted_MT2->Scale(1/hist_unboosted_MT2->Integral());
  TFile *f6 = new TFile("ghist_unboost_mt2.root", "UPDATE");
  hist_unboosted_MT2->Write(sample_desc, TObject::kOverwrite);

  hist_pt_b->Scale(1/hist_pt_b->Integral());
  TFile *f7 = new TFile("ghist_pt_b.root", "UPDATE");
  hist_pt_b->Write(sample_desc, TObject::kOverwrite);

  hist_pt_sub_b->Scale(1/hist_pt_sub_b->Integral());
  TFile *f8 = new TFile("ghist_pt_sub_b.root", "UPDATE");
  hist_pt_sub_b->Write(sample_desc, TObject::kOverwrite);

  hist_pt_tau_p->Scale(1/hist_pt_tau_p->Integral());
  TFile *f9 = new TFile("ghist_pt_tau_p.root", "UPDATE");
  hist_pt_tau_p->Write(sample_desc, TObject::kOverwrite);

  hist_pt_tau_m->Scale(1/hist_pt_tau_m->Integral());
  TFile *f10 = new TFile("ghist_pt_tau_m.root", "UPDATE");
  hist_pt_tau_m->Write(sample_desc, TObject::kOverwrite);

  hist_pt_nu->Scale(1/hist_pt_nu->Integral());
  TFile *f11 = new TFile("ghist_pt_nu.root", "UPDATE");
  hist_pt_nu->Write(sample_desc, TObject::kOverwrite);

  hist_pt_b->Scale(1/hist_pt_antinu->Integral());
  TFile *f12 = new TFile("ghist_pt_antinu.root", "UPDATE");
  hist_pt_antinu->Write(sample_desc, TObject::kOverwrite);

} // End void.
