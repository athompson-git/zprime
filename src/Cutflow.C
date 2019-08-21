#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
R__LOAD_LIBRARY(libExRootAnalysis)
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
#include <TMath.h>

#include <cmath>
#include <vector>
#include <string>
#include <utility>
#include <iostream>

// Main macro.
void Cutflow(string run_name) {
  gSystem->Load("libDelphes.so");
  gSystem->Load("libExRootAnalysis");



  // Create a chain of root trees.
  TChain chain("Delphes");

  // WZ(tautau)
  if (run_name == "wztautau") {
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_lo_madspin_z-tautau/Events/run_01_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_lo_madspin_z-tautau/Events/run_02_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/samples/wz/wz_runs_01-09.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/samples/wz/wz_runs_10-19.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/samples/wz/wz_runs_20-29.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/samples/wz/wz_runs_30-39.root");
  }

  if (run_name == "wz_trimmed") {
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/samples/wz/vanilla/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/samples/wz/vanilla/tag_2_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/samples/wz/vanilla/tag_3_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/samples/wz/vanilla/tag_4_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/samples/wz/vanilla/tag_5_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/samples/wz/vanilla/tag_6_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/samples/wz/vanilla/tag_7_delphes_events.root");
  }

  // ttZp samples.
  if (run_name == "signal_50") {
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_02_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_03_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_04_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_05_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_06_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_07_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_08_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_09_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_46_decayed_1/tag_2_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_47_decayed_1/tag_2_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_48_decayed_1/tag_2_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_49_decayed_1/tag_2_delphes_events.root");
  }

  if (run_name == "signal_40") {
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_lo/Events/run_02_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_14_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_15_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_16_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_17_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_42_decayed_1/tag_2_delphes_events.root"); // Guess!
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_43_decayed_1/tag_2_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_44_decayed_1/tag_2_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_45_decayed_1/tag_2_delphes_events.root");
  }

  if (run_name == "signal_70") {
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_10_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_11_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_12_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_13_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_32_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_33_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_lo/Events/run_03_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_50_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_51_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_52_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_53_decayed_1/tag_1_delphes_events.root");
  }

  if (run_name == "signal_100") {
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_20_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_21_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_22_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_23_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_34_decayed_1/tag_2_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_35_decayed_1/tag_2_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_36_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_37_decayed_1/tag_2_delphes_events.root");
  }

  if (run_name == "signal_120") {
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_24_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_25_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_26_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_27_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_38_decayed_1/tag_2_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_39_decayed_1/tag_2_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_40_decayed_1/tag_2_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_41_decayed_1/tag_2_delphes_events.root");
  }

  if (run_name == "signal_150") {
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_28_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_29_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_30_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin-ditau_redux/Events/run_31_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin_ditau_v3/Events/run_02_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin_ditau_v3/Events/run_03_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin_ditau_v3/Events/run_04_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttzp_madspin_ditau_v3/Events/run_05_decayed_1/tag_1_delphes_events.root");
  }


  // ttH samples
  if (run_name == "ttH") {
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttH/Events/run_02/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttH/Events/run_03/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttH/Events/run_04/tag_1_delphes_events.root");
  }

  // ttZ samples
  if (run_name == "ttZ") {
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttZ/Events/run_01/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttZ/Events/run_02/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttZ/Events/run_03/tag_1_delphes_events.root");
  }
  if (run_name == "ttz-off-shell") {
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tttata_lo/Events/run_03_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tttata_lo/Events/run_04_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tttata_lo/Events/run_05_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tttata_lo/Events/run_06_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tttata_v2/Events/run_02/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tttata_v2/Events/run_03/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tttata_v2/Events/run_04/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tttata_v2/Events/run_05/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tttata_v2/Events/run_06/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tttata_v2/Events/run_07/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tttata_v2/Events/run_08/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tttata_v2/Events/run_09/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tttata_v2/Events/run_10/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tttata_v2/Events/run_11/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tttata_v2/Events/run_12/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tttata_v2/Events/run_13/tag_1_delphes_events.root");
  }

  if (run_name == "ttW") {
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_01/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_02/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_03/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_04/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_06/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_07/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_08/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_09/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_11/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_12/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_13/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_14/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_16/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_17/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_18/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_19/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_21/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_22/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_23/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_24/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_26/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_27/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_28/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_29/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_30/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_31/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_32/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_33/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_34/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_35/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_36/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_37/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_38/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_39/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_40/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_41/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_42/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_43/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_44/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_45/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_46/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_47/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_48/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_49/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_50/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_51/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_52/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_53/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_54/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_55/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_56/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_57/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_58/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_59/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_60/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_61/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_62/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_63/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_64/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_65/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_66/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_67/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_68/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_69/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_70/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_71/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_72/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_73/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_74/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_75/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_76/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_77/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_78/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_79/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_80/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_81/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_82/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_83/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_84/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_85/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_86/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_87/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_88/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_89/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_90/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_91/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_92/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_93/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_94/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_95/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_96/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_97/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_98/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_99/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_100/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_101/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_102/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_103/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_104/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_105/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_106/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_107/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_108/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_109/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_110/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_111/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_112/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_113/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_114/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_115/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_116/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_117/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_118/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_119/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_120/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_121/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_122/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_123/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_124/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_125/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_126/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ttw_012jets_5f_lo/Events/run_127/tag_1_delphes_events.root");
  }

  if (run_name == "ttbar") {
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tt0123jets_5f_madspin/Events/run_01_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tt0123jets_5f_madspin/Events/run_02_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tt0123jets_5f_madspin/Events/run_03_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tt0123jets_5f_madspin/Events/run_05_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tt0123jets_5f_madspin_redux/Events/run_01_decayed_1/tag_2_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tt0123jets_5f_madspin_redux/Events/run_02_decayed_1/tag_2_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tt0123jets_5f_madspin_redux/Events/run_03_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tt0123jets_5f_madspin_redux/Events/run_04_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tt0123jets_5f_madspin_redux/Events/run_05_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tt0123jets_5f_madspin_redux/Events/run_06_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tt0123jets_5f_madspin_redux/Events/run_07_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tt0123jets_5f_madspin_redux/Events/run_08_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tt0123jets_5f_madspin_redux/Events/run_09_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tt0123jets_5f_madspin_redux/Events/run_10_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tt0123jets_5f_madspin_redux/Events/run_11_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tt0123jets_5f_madspin_redux/Events/run_12_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tt0123jets_5f_madspin_v3/Events/run_04_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tt0123jets_5f_madspin_v3/Events/run_05_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tt0123jets_5f_madspin_v3/Events/run_06_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tt0123jets_5f_madspin_v3/Events/run_07_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tt0123jets_5f_madspin_v3/Events/run_08_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/tt0123jets_5f_madspin_v3/Events/run_09_decayed_1/tag_1_delphes_events.root");
  }

  if (run_name == "zz" || run_name == "EWK") {
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/zz_012jets_5f_lo/Events/run_01/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/zz_012jets_5f_lo/Events/run_02/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/zz_012jets_5f_lo/Events/run_03/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/zz_012jets_5f_lo/Events/run_04/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/zz_012jets_5f_lo/Events/run_05/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/zz_012jets_5f_lo/Events/run_06/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/zz_012jets_5f_lo/Events/run_07/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/zz_012jets_5f_lo/Events/run_08/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/zz_012jets_5f_lo/Events/run_09/tag_1_delphes_events.root");
  }

  if (run_name == "wz" || run_name == "EWK" || run_name == "wz_trimmed") {
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_05/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_06/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_07/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_08/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_09/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_10/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_11/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_12/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_13/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_14/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_15/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_16/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_17/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_18/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_19/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_20/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_21/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_22/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_23/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_24/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_25/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_26/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_27/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_28/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_29/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_30/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_31/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_32/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_33/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_34/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_35/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_36/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_37/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_38/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_39/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_40/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_41/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_42/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_43/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_44/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_45/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_46/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_47/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_48/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_49/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_50/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_51/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_52/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_53/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/wz_012jets_5f_madspin_lo/Events/run_54/tag_1_delphes_events.root");
  }

  if (run_name == "ww" || run_name == "EWK") {
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_01_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_02_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_03_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_04_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_05_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_06_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_07_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_08_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_09_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_10_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_11_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_12_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_13_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_14_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_15_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_16_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_17_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_18_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_19_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_20_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_21_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_22_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_23_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_24_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_25_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_26_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_27_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_28_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_29_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_30_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_31_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_32_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_33_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_34_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_35_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_36_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_37_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_38_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_39_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_40_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_41_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_42_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_43_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_44_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_45_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_46_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_47_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_48_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_49_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_50_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_51_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_52_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_53_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_54_decayed_1/tag_1_delphes_events.root");
    chain.Add("/fdata/hepx/store/user/thompson/zprime_ditau/ww_012jets_5f_madspin_lo/Events/run_55_decayed_1/tag_1_delphes_events.root");
  }


  // Set the remaining branch addresses to zero.
  //chain.Print();
  chain.SetBranchStatus("*", 0);
  chain.SetBranchStatus("Jet*", 1);
  chain.SetBranchStatus("MissingET*", 1);
  chain.SetBranchStatus("Event*", 1);
  chain.SetBranchStatus("Muon*", 1);
  chain.SetBranchStatus("Electron*", 1);

  // Create object of class ExRootTreeReader.
  ExRootTreeReader *tree_reader = new ExRootTreeReader(&chain);
  Long64_t number_of_entries = tree_reader->GetEntries();

  // Get pointers to branches used in this analysis.
  TClonesArray *branch_jet = tree_reader->UseBranch("Jet");
  TClonesArray *branch_met = tree_reader->UseBranch("MissingET");
  TClonesArray *branch_event = tree_reader->UseBranch("Event");
  TClonesArray *branch_electron = tree_reader->UseBranch("Electron");
  TClonesArray *branch_muon = tree_reader->UseBranch("Muon");


  // Write out a trimmed TTree with only the objects in events that pass selection.
  TTree *out_tree = new TTree("out_tree","DataTree");

  Event *evt = new Event();
  Float_t tau_arr[4];
  Float_t btag_arr[4];
  Float_t jet_arr[4];
  Float_t met_arr[4];
  Float_t lep1_arr[4];
  Float_t lep2_arr[4];
  Float_t wgt;
  //Int_t nb;
  out_tree->Branch("TauBranch", &tau_arr, "tau_arr[4]/F");
  out_tree->Branch("Lep1Branch", &lep1_arr, "lep1_arr[4]/F");
  out_tree->Branch("Lep2Branch", &lep2_arr, "lep2_arr[4]/F");
  out_tree->Branch("BTagBranch", &btag_arr, "btag_arr[4]/F");
  out_tree->Branch("JetBranch", &jet_arr, "jet_arr[4]/F");
  out_tree->Branch("METBranch", &met_arr, "met_arr[4]/F");
  out_tree->Branch("Weight", &wgt, "weight/F");
  //out_tree->Branch("Nb", &nb, "Nb/I");

  // EVENT LOOP.
  Int_t accepted_events = 0;
  Int_t accepted_events_before_ss = 0;
  for (Int_t entry = 0; entry < number_of_entries; ++entry) {
    if (entry % 10000 == 0) printf("On event %d / %lld \n", entry, number_of_entries);
    //if (accepted_events == 2000) break;
    tree_reader->ReadEntry(entry);
    HepMCEvent *event = (HepMCEvent*) branch_event->At(0);


    Double_t nentries = number_of_entries;
    Double_t weight = event->Weight;
    Double_t root_file_event_size = chain.GetTree()->GetEntries();
    Double_t reweight = weight * (root_file_event_size / nentries);

    // Apply MET cut right away.
    MissingET *ETMiss = (MissingET *) branch_met->At(0);
    if (ETMiss->MET < 30.) continue;


    vector<int> bottom_jets;
    vector<int> tau_jets;
    vector<int> light_jets;
    vector<int> e_candidates;
    vector<int> mu_candidates;


    bool os = false;
    bool found_jets = false;
    bool found_btag = false;
    bool found_e = false;
    bool found_mu = false;
    Double_t n_mu = 0;
    Double_t n_e = 0;


    // PRESELECTION: 1 hadronic tau, 1 oppositely charged lepton (e/mu), 1 btag
    // JET LOOP
    for (unsigned j = 0; j < branch_jet->GetEntries(); ++j) {
      Jet *jet = (Jet*) branch_jet->At(j);
      if (fabs(jet->Eta) > 2.4) continue;
      if (jet->BTag && !(jet->TauTag)) {
        if (jet->PT < 20.) continue;
        bottom_jets.push_back(j);
      } else if (!(jet->BTag) && !(jet->TauTag)) {
        if (jet->PT < 30.) continue;
        light_jets.push_back(j);
      } else if (jet->TauTag && !(jet->BTag)) {
        if (jet->PT < 30.) continue;
        tau_jets.push_back(j);
      }
    }

    if (tau_jets.size() != 1) continue;
    if (bottom_jets.size() < 1) continue;
    if ((light_jets.size() + bottom_jets.size()) < 2) continue; // mult req

    // Electron and Muon loops.
    for (unsigned i = 0; i < branch_electron->GetEntries(); ++i) {
      Electron *e = (Electron*) branch_electron->At(i);
      if (e->PT < 26. || fabs(e->Eta) > 2.1) continue;
      e_candidates.push_back(i);
    }
    for (unsigned i = 0; i < branch_muon->GetEntries(); ++i) {
      Muon *m = (Muon*) branch_muon->At(i);
      if (m->PT < 23. || fabs(m->Eta) > 2.4) continue;
      mu_candidates.push_back(i);
    }
    if ((e_candidates.size() + mu_candidates.size()) != 2) continue;
    n_e = e_candidates.size();
    n_mu = mu_candidates.size();


    // BTAG SELECTION
    Jet *btag, *btag_2;
    btag = (Jet *) branch_jet->At(bottom_jets[0]);
    if (bottom_jets.size() > 1) {
      btag_2 = (Jet *) branch_jet->At(bottom_jets[1]);
    }
    if (bottom_jets.size() == 1 && light_jets.size() > 0) {
      btag_2 = (Jet *) branch_jet->At(light_jets[0]);
    }

    // HADRONIC TAU SELECTION.
    Jet *tau_h;
    tau_h = (Jet *) branch_jet->At(tau_jets[0]);

    // LEPTON SELECTION
    Electron *electron;
    Muon *muon;
    Double_t best_dr_el = 999.;
    Double_t best_dr_mu = 999.;
    if (n_e > 0) {
      for (Int_t e_i = 0; e_i < n_e; ++e_i) {
        Electron *this_e = (Electron *) branch_electron->At(e_candidates[e_i]);
        Double_t ell_tau_dr = (this_e->P4()).DeltaR(tau_h->P4());
        if (this_e->Charge != tau_h->Charge && ell_tau_dr < best_dr_el) {
          found_e = true;
          best_dr_el = ell_tau_dr;
          electron = (Electron *) branch_electron->At(e_i);
        }
      }
    }
    if (n_mu > 0) {
      for (Int_t mu_i = 0; mu_i < n_mu; ++mu_i) {
        Muon *this_mu = (Muon *) branch_muon->At(mu_candidates[mu_i]);
        Double_t mu_tau_dr = (this_mu->P4()).DeltaR(tau_h->P4());
        if (this_mu->Charge != tau_h->Charge && mu_tau_dr < best_dr_mu) {
          found_mu = true;
          best_dr_mu = mu_tau_dr;
          muon = (Muon *) branch_muon->At(mu_i);
        }
      }
    }
    if (!(found_mu || found_e)) continue;  // preselection continue
    // Decide which to use.
    bool use_el = !found_mu;
    bool use_mu = !found_e;
    if (found_mu && found_e) {
      if (cos(electron->Phi - ETMiss->Phi) > cos(muon->Phi - ETMiss->Phi)) {
        use_el = true;
      } else {
        use_mu = true;
      }
    }

    int lepton_charge = 0;
    TLorentzVector *lepton_p4 = new TLorentzVector();
    TLorentzVector *lepton2_p4 = new TLorentzVector();

    if (use_el) {
      lepton_p4 = new TLorentzVector(electron->P4());
      lepton_charge = electron->Charge;
    } else {
      lepton_charge = muon->Charge;
      lepton_p4 = new TLorentzVector(muon->P4());
    }

    // SS dilepton requirement and second lepton identification.
    int ss_lep = 0;
    for (int jj = 0; jj < e_candidates.size(); ++jj) {
      Electron *e = (Electron *) branch_electron->At(e_candidates[jj]);
      if (e->Charge == lepton_charge) {
        ss_lep++;
        if (ss_lep == 2) {
          lepton2_p4 = new TLorentzVector(e->P4());
        }
        break;
      }
    }
    for (int ii = 0; ii < mu_candidates.size(); ++ii) {
      Muon *mu = (Muon *) branch_muon->At(mu_candidates[ii]);
      if (mu->Charge == lepton_charge) {
        ss_lep++;
        if (ss_lep == 2) {
          lepton2_p4 = new TLorentzVector(mu->P4());
        }
        break;
      }
    }
    accepted_events_before_ss++;
    if (ss_lep < 2) continue;

    // END OF PRESELECTION

    // Fill TTree with all particles.
    TLorentzVector *tau_p4 = new TLorentzVector(tau_h->P4());
    TLorentzVector *b_p4 = new TLorentzVector(btag->P4());
    TLorentzVector *b2_p4 = new TLorentzVector(btag_2->P4());
    TLorentzVector *met_p4 = new TLorentzVector();
    met_p4->SetPtEtaPhiE(ETMiss->MET, 0, ETMiss->Phi, ETMiss->MET);

    accepted_events++;

    wgt = reweight;
    //nb = bottom_jets.size();

    tau_arr[0] = tau_p4->Pt();
    tau_arr[1] = tau_p4->Eta();
    tau_arr[2] = tau_p4->Phi();
    tau_arr[3] = tau_p4->E();

    btag_arr[0] = b_p4->Pt();
    btag_arr[1] = b_p4->Eta();
    btag_arr[2] = b_p4->Phi();
    btag_arr[3] = b_p4->E();

    jet_arr[0] = b2_p4->Pt();
    jet_arr[1] = b2_p4->Eta();
    jet_arr[2] = b2_p4->Phi();
    jet_arr[3] = b2_p4->E();

    met_arr[0] = met_p4->Pt();
    met_arr[1] = met_p4->Eta();
    met_arr[2] = met_p4->Phi();
    met_arr[3] = met_p4->E();

    lep1_arr[0] = lepton_p4->Pt();
    lep1_arr[1] = lepton_p4->Eta();
    lep1_arr[2] = lepton_p4->Phi();
    lep1_arr[3] = lepton_p4->E();

    lep2_arr[0] = lepton2_p4->Pt();
    lep2_arr[1] = lepton2_p4->Eta();
    lep2_arr[2] = lepton2_p4->Phi();
    lep2_arr[3] = lepton2_p4->E();

    out_tree->Fill();


    delete tau_p4;
    delete b_p4;
    delete b2_p4;
    delete met_p4;
    delete lepton_p4;
    delete lepton2_p4;


  } // End event loop.

  printf("%d / %lld accepted \n", accepted_events, number_of_entries);
  //printf("%d / %lld accepted before SS lepton req\n", accepted_events_before_ss, number_of_entries);

  // Write NTuples to file.
  TFile *f = new TFile("../ntuples/analysis_tree.root", "UPDATE");
  out_tree->Write(run_name.c_str(), TObject::kOverwrite);


}  // End macro.


