// Skeletal version of CutFlow.C. Made for validation on ttbar and Z' samples.
// USAGE on Brazos:
// gSystem->Load("<path/to/your/Delphes/libDelphes");
// .x cutflow.C("/fdata/scratch/mexanick/ttbar");
//
// On my pc:
// .x cutflow_MT2.C("/Users/adrianthompson/physics/zprime/samples/ttbar");


#ifdef __CLING__
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TStyle.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TEfficiency.h>

#include "lester_mt2_bisect.h"
#include <vector>
#include <string>
#include <iostream>
#include <cmath>


double mt2(TLorentzVector *vis1, TLorentzVector *vis2, TLorentzVector *miss) {
  double mVisA = 0.1056; //vis1.M(); // Mass of visible object on side A.
  double pxA = vis1->Px(); // x momentum of visible object on side A.
  double pyA = vis1->Py(); // y momentum of visible object on side A.

  double mVisB = 0.1056; //vis2.M(); // Mass of visible object on side B.
  double pxB = vis2->Px(); // x momentum of visible object on side B.
  double pyB = vis2->Py(); // y momentum of visible object on side B.

  double pxMiss = miss->Px(); // x component of missing transverse momentum.
  double pyMiss = miss->Py(); // y component of missing transverse momentum.

  double chiA = 0.0; // Hypothesised mass of invisible on side A.
  double chiB = 0.0; // Hypothesised mass of invisible on side B.

  double desiredPrecisionOnMt2 = 0.0; // Algo aims for machine precision.

  asymm_mt2_lester_bisect lester;
  return lester.get_mT2(mVisA, pxA, pyA, mVisB, pxB, pyB,
                        pxMiss, pyMiss, chiA, chiB,
                        desiredPrecisionOnMt2);
}


float correlation(float i, float j, float ij) {
  if (i == j) ij = 1.000000;
  return ((ij - i) * j) / std::sqrt(i * j * (1 - i) * (1 - j));
}



void cutflow_MT2(const char *sample_directory = 0) {
  gSystem->Load("libDelphes.so");

  std::string ttbar_samples;
  ttbar_samples += sample_directory;
  ttbar_samples += "/tt_*.root";

  std::cout << "taking samples from " << ttbar_samples.c_str() << std::endl;

  TChain chain_("Delphes");
  chain_.Add(ttbar_samples.c_str());
//  chain_.Add("/Users/adrianthompson/physics/zprime/samples/zp_500GeV/*");

  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain_);
  Long64_t numberOfEntries = treeReader->GetEntries();

  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchMPT = treeReader->UseBranch("MissingET");

  // Book TProfiles and TEfficiencies.
  TProfile *cutflow = new TProfile("cutflow",
                                   "(2b & 2OS muons) & max(SBM)>170, E_{T}^{miss} / M(#mu^{+}, #mu^{-}), HTLT<0, UnbMT2 < 10;cut number;percent passing",
                                   5,-.5,4.5);

  // Make an array of TProfiles.
  TProfile e1 = TProfile("e1", "max{M(#mu, b)}", 4, -.5, 3.5);
  TProfile e2 = TProfile("e2", "E_{T}^{miss} / M(#mu^{+}, #mu^{-})", 4, -.5, 3.5);
  TProfile e3 = TProfile("e3", "H_{T} - L_{T}", 4, -.5, 3.5);
  TProfile e4 = TProfile("e3", "Unboosted MT2", 4, -.5, 3.5);
  TProfile *efficiency_matrix [4] = {&e1, &e2, &e3, &e4};

  // Begin event loop.
  for (Int_t entry = 0; entry < numberOfEntries; ++entry) {
    if (entry % 10000 == 0) printf("On event %d / %lld \n", entry, numberOfEntries);

    treeReader->ReadEntry(entry);

    std::vector<int> bPosition (0,0);
    int nonbPosition = -1;
    float nonbPT = 0.;


    // If the event contains at least one jet, push back the vector of b's and
    // record the branch position of the highest-pT non-b jet.
    for (int j = 0; j < branchJet->GetEntries(); ++j) {
      Jet *jet = (Jet*) branchJet->At(j);
      if (jet->PT < 30.) continue;  // Implicitly require pT > 30 GeV jets.

      if (jet->BTag) {
        bPosition.push_back(j);
      } else if (nonbPT < jet->PT) {
        nonbPosition = j;
        nonbPT = jet->PT;
      }
    }

    bool DiBottom = false;
    bool BottomJet = false;
    if (bPosition.size() >= 2) DiBottom = true;
    if (bPosition.size() == 1 && nonbPT > 0.) BottomJet = true;

    Muon *mu1, *mu2;
    bool OSMuons = false;

    // If the event contains at least 2 muons, check if they are OS.
    if (branchMuon->GetEntries() > 1) {
      mu1 = (Muon*) branchMuon->At(0);
      mu2 = (Muon*) branchMuon->At(1);

      if((mu1->Charge) != (mu2->Charge)) OSMuons = true;
    }

    // Skip to the next event if preselection requirements are unmet.
    cutflow->Fill(0., (OSMuons && (DiBottom || BottomJet)));
    if (!OSMuons) continue;
    if (!(DiBottom || BottomJet)) continue;

    // Jet assignments.
    Jet *b1, *b2;
    if(DiBottom){
      b1 = (Jet*) branchJet->At(bPosition[0]);
      b2 = (Jet*) branchJet->At(bPosition[1]);
    }
    else if(BottomJet){
      b1 = (Jet*) branchJet->At(bPosition[0]);
      b2 = (Jet*) branchJet->At(nonbPosition);
    }

    // Initialize cuts.
    bool PassMissingET = false;
    bool PassHTLT = false;
    bool PassMuJetMass = false;
    bool PassUnboostedMT2 = false;

    // Cut (1): Calculate max{M(mu,j)}.
    float Pm1_1 = ((mu1->P4()) + (b1->P4())).M();
    float Pm1_2 = ((mu2->P4()) + (b2->P4())).M();
    float Pm2_1 = ((mu2->P4()) + (b1->P4())).M();
    float Pm2_2 = ((mu1->P4()) + (b2->P4())).M();
    float MassPair1 = 0.;
    float MassPair2 = 0.;

    if (fabs(Pm1_1 - Pm1_2) < fabs(Pm2_1 - Pm2_2)) {
      MassPair1 = Pm1_1;
      MassPair2 = Pm1_2;
    } else {
      MassPair1 = Pm2_1;
      MassPair2 = Pm2_2;
    }

    if ((MassPair1 > 170.) || (MassPair2 > 170.)) PassMuJetMass = true;

    // Cut (2): Calculate Dimuon mass and normalized Missing ET.
    float Mmm = ((mu1->P4()) + (mu2->P4())).M();
    MissingET *MPT = (MissingET*) branchMPT->At(0);
    if ((5. * MPT->MET) < Mmm) PassMissingET = true;

    // Cut (3): Calculate HT - LT.
    float HTLT = (b1->PT) + (b2->PT) - (mu1->PT) - (mu2->PT);
    if (HTLT < 0) PassHTLT = true;



    TLorentzVector mu1_p4 = mu1->P4();
    TLorentzVector mu2_p4 = mu2->P4();
    TLorentzVector met_p4;
    met_p4.SetPtEtaPhiE(MPT->MET, 0., MPT->Phi, MPT->MET);
    TLorentzVector b1_p4 = b1->P4();
    TLorentzVector b2_p4 = b2->P4();

    // (4) Calculate MT2 (http://www.hep.phy.cam.ac.uk/~lester/mt2/).
    double MT2 = mt2(&mu1_p4, &mu2_p4, &met_p4);

    // UNBOOSTED system.
    // First calculate recoil in the transverse plane.
    TLorentzVector jet_recoil = b1_p4 + b2_p4;
    TLorentzVector unboost_mu1 = mu1_p4 + jet_recoil;
    TLorentzVector unboost_mu2 = mu2_p4 + jet_recoil;
    TLorentzVector unboost_met = met_p4 + jet_recoil;

    // (5) Calculate UNBOOSTED MT2.
    double unboosted_MT2 = mt2(&unboost_mu1, &unboost_mu2, &unboost_met);
    if (unboosted_MT2 < 10.) PassUnboostedMT2 = true;

    bool cut_vector [4] = {PassMuJetMass, PassMissingET, PassHTLT, PassUnboostedMT2};

    // Fill TProfiles.
    cutflow->Fill(1., PassMuJetMass);
    if (PassMuJetMass) {
      cutflow->Fill(2., PassMissingET);
      if (PassMissingET) {
        cutflow->Fill(3., PassHTLT);
        if (PassHTLT) {
          cutflow->Fill(4., PassUnboostedMT2);
        }
      }
    }


    // Fill correlation matrix.
    // matrix elements = (TProfile #, bin #)
    for (int i = 0; i < 4; ++i) {
      // Fill the diagonal elements.
      efficiency_matrix[i]->Fill(i, cut_vector[i]);
      for (int j = 0; j < 4; ++j) {
        // Fill the off-diagonal elements.
        if (j == i) continue;
        if (cut_vector[i]) {
          // Fill the correlation element (i,j).
          efficiency_matrix[i]->Fill(j, cut_vector[j]);
        }
      }
    }

  } // End event loop.

  // Print out efficiency information.
  for (int i=0; i < cutflow->GetNbinsX(); ++i) {
    printf("%d: %f +/- %f \n" ,i,
           cutflow->GetBinContent(i+1),
           cutflow->GetBinError(i+1));
  }


  for (unsigned k = 0; k < 4; ++k) {
    printf("(%f +/- %f) (%f +/- %f) (%f +/- %f) (%f +/- %f)\n",
           efficiency_matrix[k]->GetBinContent(1),
           efficiency_matrix[k]->GetBinError(1),
           efficiency_matrix[k]->GetBinContent(2),
           efficiency_matrix[k]->GetBinError(2),
           efficiency_matrix[k]->GetBinContent(3),
           efficiency_matrix[k]->GetBinError(3),
           efficiency_matrix[k]->GetBinContent(4),
           efficiency_matrix[k]->GetBinError(4));
  }

  // Print correlation matrix.
  float rho [4][4];

  for (int x = 0; x < 4; ++x) {
    for (int y = 0; y < 4; ++y) {
      rho[x][y] = efficiency_matrix[x]->GetBinContent(y+1);
    }
  }

  for (int l = 0; l < 4; ++l) {
    printf("%f %f %f %f \n", correlation(rho[l][l], rho[0][0], rho[0][l]),
           correlation(rho[l][l], rho[1][1], rho[1][l]),
           correlation(rho[l][l], rho[2][2], rho[2][l]),
           correlation(rho[l][l], rho[3][3], rho[3][l]));
  }

}
