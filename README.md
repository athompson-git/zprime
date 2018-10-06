# zprime

#### Paper:
[Bottom-quark fusion processes at the LHC for probing Z' models and B-meson decay anomalies](https://arxiv.org/pdf/1707.07016.pdf)

#### Macros and Tools:

* DiTauAnalyzer.C	: Delphes TTree macro with some custom functions and an MT2 calculator
* DiTauGenAnalyzer.C	: LHEF TTree macro with an MT2 calculator
* GenJetMatcher.C	: Matches Pythia8 `GenParticle` class Truth-level objects with Delphes physics objects.
* NPlotOverlay.C	: Takes a `.root` file with histograms inside, loops over the histos and plots them on the same graph
* TopoPlotter.C	: Specialized plotting tool for making (*pTcos(dPhi)*, *pTsin(dPhi)*) topological plots 
* WJetsAnalyzer.C	: Specialized TTree macro for analyzing W+jets and calculating DZeta
* cutflow.C	: simplified TTree analyzer with no plotting; just cut efficiency calculations
* cutflow_MT2.C	: same as cutflow.C but with an MT2 calculator
* lester_mt2_bisect.h	: necessary header file for MT2 calculation
* overlap.h	: custom measure of the overlap between two histograms
