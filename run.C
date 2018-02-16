//#ifdef __CLING__
//#include "classes/DelphesClasses.h"
//#include "external/ExRootAnalysis/ExRootTreeReader.h"
//#endif

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
#include "DiTauAnalyzer.C"

void run() {
  gSystem->Load("/Users/adrianthompson/MG5_aMC_v2_6_1/Delphes/libDelphes");
  gROOT->LoadMacro("DiTauAnalyzer.C");
  DiTauAnalyzer("zp_ditau_sample01/Events/200GeV_80kEvents/tag_1_delphes_events.root", "Z' (200 GeV)", 20);
//  if (sample == 0) {
   // gROOT->LoadMacro("DiTauAnalyzer.C(\"zp_ditau_sample01/Events/200GeV_80kEvents/tag_1_delphes_events.root\", \"Z\' (200 GeV)\", 20)");
//  }
/*  if (sample == 1) {
    gROOT->LoadMacro("DiTauAnalyzer.C(\"zp_ditau_sample01/Events/350GeV_80kEvents/tag_1_delphes_events.root\", \"Z\' (350 GeV)\", 20);");
  }
  if (sample == 2) {
    gROOT->LoadMacro("DiTauAnalyzer.C(\"sample_40k_zp_ditau.root\", \"Z\' (500 GeV)\", 20);");
  }
  if (sample == 3) {
    gROOT->LoadMacro("DiTauAnalyzer.C(\"zp_ditau_sample01/Events/1000GeV_80kEvents/tag_1_delphes_events.root\", \"Z\' (1 TeV)\", 20);");
  }
  if (sample == 4) {
    gROOT->LoadMacro("DiTauAnalyzer.C(\"ttbar_ditau_02/Events/160kEvents_run/tag_1_delphes_events.root\", \"tt~ ditau\", 20);");
  }
  if (sample == 5) {
    gROOT->LoadMacro("DiTauAnalyzer.C(\"wjets_sample01/Events/run_05/tag_1_delphes_events.root\", \"W + Jets\", 20);");
  }*/
}
