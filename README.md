# Installing MadGraph + Pythia8 + Delphes Software Suite
## MadGraph Installation (Brazos Cluster)
Login to brazos via
```
ssh <user>@login.brazos.tamu.edu
```

You may install MadGraph in your home directory (`/home/<user>/`). This is the default working directory you end up at after SSH’ing into the cluster.
Next, go to to <https://launchpad.net/mg5amcnlo>. On the right hand side, right click on the link for “latest version” and select copy link address. This should be a link such as 
```
<https://launchpad.net/mg5amcnlo/2.0/2.6.x/+download/MG5_aMC_v2.6.3.tar.gz>
```
In your terminal on Brazos, call
```
wget https://launchpad.net/mg5amcnlo/2.0/2.6.x/+download/MG5_aMC_v2.6.3.tar.gz
```
```
tar xvzf MG5_aMC_v2.6.3.tar.gz
```
At this point there is no need to make/compile MadGraph; its binaries are already compiled.
The MadGraph5 executable is located in `~/MG5_aMC_v2_6_3/bin/mg5_aMC`
Copy this path and edit your .bashrc with your favorite editor (vim, nano etc.);
```
vim ~/.bashrc
```
Somewhere towards the bottom of the file add
```
alias mg5='/home/thompson/MG5_aMC_v2_6_3/bin/mg5_aMC’
```
making sure to change the `vX_Y_Z` to the version you grabbed where appropriate.
Now you can run madgraph in any directory with the command `mg5`.





# Top-production and Ditau-decay Zprime

#### Previous Paper:
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
