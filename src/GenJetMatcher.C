// An analyzer to produce key selection criteria plots for Z' -> ditau processes.
// arXiv:1707.07016v1
// USAGE:
// gSystem->Load("<path_to_delphes>/Delphes/libDelphes.so");
// .x DiTauAnalyzer.C("input_file.root", "Z' 500 GeV", nbins);

#ifdef __CLING__
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif

#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TEfficiency.h>
#include <TCanvas.h>

#include "lester_mt2_bisect.h"
#include <vector>
#include <string>
#include <utility>
#include <iostream>




// List of accepted Pythia status codes.
bool GoodStatus(Int_t x) {
  Int_t codes[] = {23, 24};
  return std::find(std::begin(codes), std::end(codes), x) != std::end(codes);
}




// Main macro.
void GenJetMatcher(string sample_name, int nbins, bool verbose = false) {
  gSystem->Load("libDelphes.so");

  // Create a chain of root trees.
  TChain chain("Delphes");

  // 200 GeV Z' samples.
  if (sample_name == "200") {
    chain.Add("/Users/adrianthompson/physics/zprime/zp_ditau_sample01/Events/200GeV_50k_gtau05_run01_NarrowWidth/tag_1_delphes_events.root");
    chain.Add("/Users/adrianthompson/physics/zprime/zp_ditau_sample01/Events/200GeV_50k_gtau05_run02_NarrowWidth/tag_1_delphes_events.root");
    chain.Add("/Users/adrianthompson/physics/zprime/zp_ditau_sample01/Events/200GeV_50k_gtau05_run03_NarrowWidth/tag_1_delphes_events.root");
    chain.Add("/Users/adrianthompson/physics/zprime/zp_ditau_sample01/Events/200GeV_50k_gtau05_run04_NarrowWidth/tag_1_delphes_events.root");
  }

  // 350 GeV Z' samples.
  if (sample_name == "350") {
    chain.Add("/Users/adrianthompson/physics/zprime/zp_ditau_sample01/Events/350GeV_50k_gtau05_run01_NarrowWidth/tag_1_delphes_events.root");
    chain.Add("/Users/adrianthompson/physics/zprime/zp_ditau_sample01/Events/350GeV_50k_gtau05_run02_NarrowWidth/tag_1_delphes_events.root");
    chain.Add("/Users/adrianthompson/physics/zprime/zp_ditau_sample01/Events/350GeV_50k_gtau05_run03_NarrowWidth/tag_1_delphes_events.root");
  }

  // 500 GeV Z' samples.
  if (sample_name == "500") {
    chain.Add("/Users/adrianthompson/physics/zprime/zp_ditau_sample01/Events/500GeV_50k_gtau05_run01_NarrowWidth/tag_3_delphes_events.root");
    chain.Add("/Users/adrianthompson/physics/zprime/zp_ditau_sample01/Events/500GeV_50k_gtau05_run02_NarrowWidth/tag_1_delphes_events.root");
  }

  // ttbar samples
  if (sample_name == "ttbar") {
    chain.Add("/Users/adrianthompson/physics/zprime/ttbar_ditau_02/Events/160kEvents_run/tag_1_delphes_events.root");
    chain.Add("/Users/adrianthompson/physics/zprime/ttbar_ditau_02/Events/160kEvents_run02_xqcut30_qcut60_pdf306000/tag_1_delphes_events.root");
  }



  // Create object of class ExRootTreeReader.
  ExRootTreeReader *tree_reader = new ExRootTreeReader(&chain);
  Long64_t number_of_entries = tree_reader->GetEntries();
  Int_t accepted_events = 0;

  // Get pointers to branches used in this analysis.
  TClonesArray *branch_jet = tree_reader->UseBranch("Jet");
  TClonesArray *branch_met = tree_reader->UseBranch("MissingET");
  TClonesArray *branch_event = tree_reader->UseBranch("Event");
  TClonesArray *branch_particle = tree_reader->UseBranch("Particle");


  // Define histograms.
  // 1) 2d profile of GEN quark flavor vs. Delphes Tag type between matched pairs
  TH1F *th1f_pt_correction = new TH1F("pt_correction", "", nbins, -40., 40.);

  // Main event loop. #number_of_entries
  for(Int_t entry = 0; entry < number_of_entries; ++entry) {
    if (entry % 50000 == 0) printf("On event %d / %lld \n", entry, number_of_entries);
    if (verbose && entry == 500) break;

    tree_reader->ReadEntry(entry);
    HepMCEvent *event = (HepMCEvent*) branch_event->At(0);

    vector<Jet*> delphes_jets;
    vector<GenParticle*> gen_quarks;
    vector<Jet*> matched_jets;
    vector<GenParticle*> matched_partons;

    // Collect Delphes jets into a vector.
    for (unsigned j = 0; j < branch_jet->GetEntries(); ++j) {
      Jet *jet = (Jet*) branch_jet->At(j);
      if (jet->TauTag) continue;
      delphes_jets.push_back(jet);
    }

    // Collect GEN quarks into a vector.
    for (unsigned q = 0; q < branch_particle->GetEntries(); ++q) {
      GenParticle *quark = (GenParticle*) branch_particle->At(q);
      if (GoodStatus(quark->Status) && (fabs(quark->PID) < 7 || quark->PID == 21) ) {
        gen_quarks.push_back(quark);
      }
    }


    // Build the sorted list of partons.
    for (int j = 0; j<delphes_jets.size(); ++j) {
      Jet *j_j = (Jet*) branch_jet->At(j);
      Double_t best_deltaR = 999.0;
      GenParticle *matched_parton;

      for (int i = 0; i < gen_quarks.size(); ++i) {
        GenParticle *p_i = (GenParticle*) branch_particle->At(i);
        if ((p_i->P4()).DeltaR(j_j->P4()) < best_deltaR) {
          best_deltaR = (p_i->P4()).DeltaR(j_j->P4());
          matched_parton = p_i;
        }
      }

      if (best_deltaR < 0.2) {
        // Push back a list of pairs of ID's.
        matched_jets.push_back(j_j);
        matched_partons.push_back(matched_parton);
      }
    }

    if (matched_jets.size() == 0) continue;
    if (verbose) {
      printf("Event # %d \n", entry);
      cout << "-------------------------" << endl;
      cout << "     - Non-Tau Jets -      - Best-match Parton -         - DeltaR -" << endl;
      for (int i = 0; i < matched_jets.size(); ++i) {
        Jet *j = matched_jets[i];
        GenParticle *p = matched_partons[i];
        printf("#%d | BTag=%d pT=%f |  PID=%d Status=%d pT=%f | %f \n", i, j->BTag,
               j->PT, p->PID, p->Status, p->PT, (p->P4()).DeltaR(j->P4()));
      }

      cout << "-------------------------" << endl;
    }
    // End of selection. ////////////////////////////////////////


    // Sort jet/parton pairs by pT, delta(R), delta(pT), and fill histos.
    Jet *jet1 = matched_jets[0];
    GenParticle *parton1 = matched_partons[0];
    TLorentzVector jet1_tlv = jet1->P4();
    TLorentzVector parton1_tlv = parton1->P4();

    th1f_pt_correction->Fill(parton1_tlv.Pt() - jet1_tlv.Pt());


  } // End event loop.

  printf("%d / %lld accepted \n", accepted_events, number_of_entries);
  th1f_pt_correction->DrawNormalized();


}  // End macro.
