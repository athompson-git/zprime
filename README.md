# Installing MadGraph + Pythia8 + Delphes Software Suite
## MadGraph Installation (Brazos Cluster)
Login to brazos via
```
$ ssh <user>@login.brazos.tamu.edu
```

First you need to setup ROOT and a working GCC version.
```
source /home/hepxadmin/root-6/root-6.12.06-install/bin/thisroot.sh
module load gcc/5.2.0
```


You may install MadGraph in your home directory (`/home/<user>/`). This is the default working directory you end up at after SSH’ing into the cluster.
Next, go to to <https://launchpad.net/mg5amcnlo>. On the right hand side, right click on the link for “latest version” and select copy link address. This should be a link such as 
```
<https://launchpad.net/mg5amcnlo/2.0/2.6.x/+download/MG5_aMC_v2.6.3.tar.gz>
```
In your terminal on Brazos, call
```
$ wget https://launchpad.net/mg5amcnlo/2.0/2.6.x/+download/MG5_aMC_v2.6.3.tar.gz
```
```
$ tar xvzf MG5_aMC_v2.6.3.tar.gz
```
At this point there is no need to make/compile MadGraph; its binaries are already compiled.
The MadGraph5 executable is located in `~/MG5_aMC_v2_6_3/bin/mg5_aMC`
Copy this path and edit your .bashrc with your favorite editor (vim, nano etc.);
```
$ vim ~/.bashrc
```
Somewhere towards the bottom of the file add
```
$ alias mg5='~/MG5_aMC_v2_6_3/bin/mg5_aMC’
```
making sure to change the `vX_Y_Z` to the version you grabbed where appropriate.
Now you can run madgraph in any directory with the command `mg5`. Test it out!

```
$ mg5
MG5_aMC> generate p p > t t~
MG5_aMC> add process p p > t t~ j
MG5_aMC> output ttbar_test
MG5_aMC> exit
```
This should create a folder called `ttbar_test` that contains the feynman diagrams for the process you specify, the relevant parameter cards (`ttbar_test/Cards/`), and a folder to hold any data files from your simulation runs (`ttbar_test/Events`). You can delete this folder if you like.

Next, download the model and put it in the models folder in your madgraph installation directory.
```
$ cd ~/MG5_aMC_v2_6_3/models
$ git clone https://github.com/athompson-tamu/zprime_UFO.git

```

Now you can generate samples with the Zprime model!
```
$ mg5
MG5_aMC> import model zprime_UFO
MG5_aMC> generate p p > zp t t~, zp > ta+ ta-
MG5_aMC> output ttzp_ditau
MG5_aMC> exit
```

You should also install ExRootAnalysis, Pythia8 and DELPHES (and their dependencies) at some point. You can install them interactively through MG5:
```
$ mg5
MG5_aMC> install pythia8
MG5_aMC> install Delphes
MG5_aMC> install ExRootAnalysis
```
Running each of these install commands will no doubt try to install their various dependencies (hepmc, boost, etc.). Go ahead and agree to these installations. This whole set of installs typically takes 20 minutes to an hour at most.

Next you can try simulating a ttZp sample with pythia and Delphes.
```
$ cd ttzp_ditau
$ ./bin/madevent
```
Answer the interactive questions with the numberpad/keypad. Make sure Pythia8 and DELPHES are available and turned on. Also make sure ExRootAnalysis is turned on. Then hit ENTER to continue through the questions until the simulation begins. The beginning of the simulation should look like this:
```
INFO: Update the dependent parameter of the param_card.dat
Generating 10000 events with run name run_01
survey  run_01
INFO: compile directory
```

Then at some point MG5 should finish generating GEN-level events and pass things to Pythia8:
```
  === Results Summary for run: run_01 tag: tag_1 ===

     Cross-section :   1274 +- 2.201 pb
     Nb of events :  10000

store_events
INFO: Storing parton level results
INFO: End Parton
reweight -from_cards
decay_events -from_cards
INFO: Running Pythia8 [arXiv:1410.3012]
Splitting .lhe event file for PY8 parallelization...
Submitting Pythia8 jobs...
```
and finally, Pythia8 will pass the parton shower to Delphes for detector modelling.

Your finished job should be saved as a ROOT file in
```ttzp_ditau/Events/run_01/tag_1_delphes_events.root```



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
