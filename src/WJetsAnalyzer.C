// An analyzer to produce key histograms to distinguish Z'(ditau) from
// W + Jets background processes.
// arXiv:1707.07016v1
// USAGE:
// gSystem->Load("<path_to_delphes>/Delphes/libDelphes.so");
// .x WJetsAnalyzer.C("input_file.root", "Z' 500 GeV", nbins);

#ifdef __CLING__
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#else
class ExRootTreeReader;
class ExRootResult;
#endif

void WJetsAnalyzer(const char *file_name, const char *sample_desc, int nbins) {
  gSystem->Load("libDelphes.so");

  // Create a chain of root trees.
  TChain chain("Delphes");
  chain.Add(file_name);

  // Create object of class ExRootTreeReader.
  ExRootTreeReader *tree_reader = new ExRootTreeReader(&chain);
  Long64_t number_of_entries = tree_reader->GetEntries();

  // Get pointers to branches used in this analysis.
  TClonesArray *branch_jet = tree_reader->UseBranch("Jet");
  TClonesArray *branch_met = tree_reader->UseBranch("MissingET");

  // Book histograms.
  TH1F *hist_dzeta_85 = new TH1F("dzeta_85", "D#zeta, #alpha = 0.85", nbins,
                                 -300., 200.);
  hist_dzeta_85->GetXaxis()->SetTitle("DZeta [GeV]");
  hist_dzeta_85->GetYaxis()->SetTitle("a.u.");

  TH1F *hist_dzeta_50 = new TH1F("dzeta_50", "D#zeta, #alpha = 0.50", nbins,
                                 -300., 200.);
  hist_dzeta_50->GetXaxis()->SetTitle("DZeta [GeV]");
  hist_dzeta_50->GetYaxis()->SetTitle("a.u.");

  TH1F *hist_dzeta_15 = new TH1F("dzeta_15", "D#zeta, #alpha = 0.15", nbins,
                                 -300., 200.);
  hist_dzeta_15->GetXaxis()->SetTitle("DZeta [GeV]");
  hist_dzeta_15->GetYaxis()->SetTitle("a.u.");

  Int_t accepted_events = 0;

  // Main event loop.
  for(Int_t entry = 0; entry < number_of_entries; ++entry) {

    tree_reader->ReadEntry(entry);

    // Declare physics objects.
    Jet *jet;
    MissingET *met;
    TLorentzVector btag_jet;
    TLorentzVector second_jet;
    TLorentzVector tau_plus;
    TLorentzVector tau_minus;

    // Declare running indices.
    Int_t jet_size = branch_jet->GetEntries();
    Int_t met_size = branch_met->GetEntries();

    // Skip events if they do not have leading/subleading jets and
    // both OS Tau's.
    bool found_btag = false;
    bool found_sub_jet = false;
    bool found_tau_plus = false;
    bool found_tau_minus = false;

    if (jet_size > 3) printf("jet size = %d \n", jet_size);
    if (jet_size < 4) continue;

    // Loop over jets and find the highest pt BTag, the second highest pt, non-TauTag jet,
    // and two opposite sign TauTag jets.
    Int_t leading_btag_id = -1;
    for (Int_t ii = 0; ii < jet_size; ii++) {
      jet = (Jet*) branch_jet->At(ii);
      // Find the leading b-tagged jet.
      if (jet->BTag == 1 && jet->PT > btag_jet.Pt()) {
        btag_jet.SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, jet->Mass);
        leading_btag_id = ii;
        found_btag = true;
      }
      // Find the leading OS Tau pair.
      if (jet->TauTag == 1) {
        if (jet->Charge == 1 && jet->PT > tau_plus.Pt()) {
          tau_plus.SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, jet->Mass);
          found_tau_plus = true;
        }
        if (jet->Charge == -1 && jet->PT > tau_minus.Pt()) {
          tau_minus.SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, jet->Mass);
          found_tau_minus = true;
        }
      }
    }

    if (found_tau_plus && found_tau_minus) printf("found both");
    // Loop to find the other highest-pt jet with ID different than the leading BTag.
    for (Int_t ii = 0; ii < jet_size; ii++) {
      jet = (Jet*) branch_jet->At(ii);
      if (jet->TauTag == 0 && jet->PT > second_jet.Pt() && ii != leading_btag_id) {
        second_jet.SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, jet->Mass);
        found_sub_jet = true;
      }
    }

    // Make selection cuts.
    if (!found_tau_plus || !found_tau_minus) continue;
    if (!found_btag) continue;
    if (btag_jet.Pt() < 30.) continue;
    if (second_jet.Pt() < 30.) continue;
    if (tau_plus.Pt() < 70.) continue;
    if (tau_minus.Pt() < 70.) continue;
    accepted_events++;

    // Now that we have ID'd our 4 particles, calculate kinematics.
    // p = +, m = -
    // b = btag jet, j = other leading jet (btag or non-btag)
    TLorentzVector tau_p_b = tau_plus + btag_jet;
    TLorentzVector tau_m_j = tau_minus + second_jet;
    TLorentzVector tau_p_j = tau_plus + second_jet;
    TLorentzVector tau_m_b = tau_minus + btag_jet;
    TLorentzVector ditau = tau_plus + tau_minus;
    TLorentzVector met_p4, transverse_tau_plus, transverse_tau_minus;

    met = (MissingET*) branch_met->At(0);
    met_p4.SetPtEtaPhiE(met->MET, 0, met->Phi, met->MET);
    transverse_tau_plus.SetPtEtaPhiE(tau_plus.Pt(), 0., tau_plus.Phi(),
                                      tau_plus.Et());
    transverse_tau_minus.SetPtEtaPhiE(tau_minus.Pt(), 0., tau_minus.Phi(),
                                       tau_minus.Et());

    // Calculate D-Zeta to distinguish W+jets backgrounds.
    // (these are actually 2-component vectors in the Eta=0 plane,
    // but using TVector3 for convenience.)
    TVector3 p_vis_tau_plus = transverse_tau_plus.Vect();
    TVector3 p_vis_tau_minus = transverse_tau_minus.Vect();
    TVector3 zeta = (p_vis_tau_minus.Mag() * p_vis_tau_plus
                     + p_vis_tau_plus.Mag() * p_vis_tau_minus);  // Bisector.
    TVector3 zeta_hat = zeta * (1 / zeta.Mag());  // Unit bisector.
    Double_t p_vis = p_vis_tau_plus.Dot(zeta_hat) + p_vis_tau_minus.Dot(zeta_hat);
    Double_t p_miss = (met_p4.Vect()).Dot(zeta_hat);

    hist_dzeta_85->Fill(p_miss - 0.85 * p_vis); // alpha = 0.85
    hist_dzeta_50->Fill(p_miss - 0.50 * p_vis); // alpha = 0.50
    hist_dzeta_15->Fill(p_miss - 0.15 * p_vis); // alpha = 0.15

  } // End event loop.

  printf("Out of 40k events %d were accepted. \n", accepted_events);

  // Draw histograms and save them in a .root format.
  TCanvas *c1 = new TCanvas();
  hist_dzeta_85->Scale(1/hist_dzeta_85->Integral());
  hist_dzeta_85->Draw("HIST");

  TCanvas *c2 = new TCanvas();
  hist_dzeta_50->Scale(1/hist_dzeta_50->Integral());
  hist_dzeta_50->Draw("HIST");

  TCanvas *c3 = new TCanvas();
  hist_dzeta_15->Scale(1/hist_dzeta_15->Integral());
  hist_dzeta_15->Draw("HIST");

  TFile *f1 = new TFile("hist_dzeta_85.root", "UPDATE");
  hist_dzeta_85->Write(sample_desc, TObject::kOverwrite);

  TFile *f2 = new TFile("hist_dzeta_50.root", "UPDATE");
  hist_dzeta_50->Write(sample_desc, TObject::kOverwrite);

  TFile *f3 = new TFile("hist_dzeta_15.root", "UPDATE");
  hist_dzeta_15->Write(sample_desc, TObject::kOverwrite);

}

