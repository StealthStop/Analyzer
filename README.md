# Analyzer

A repository that is comprised of user-customized "Analyzer" modules for looping over events and producing histograms for whichever quantities desired.
Supporting plotting scripts and tools are provided.

## Setup code

To have easy access to TensorFlow and UpRoot, we need to work in a CMSSW_11_2_0_pre5 release:
```
export SCRAM_ARCH=slc7_amd64_gcc820
cmsrel CMSSW_11_2_0_pre5
cd CMSSW_11_2_0_pre5/src/
cmsenv
```

Then, check out the latest tagged version of the top tagger repository. 

```
git clone git@github.com:susy2015/TopTagger.git
cd TopTagger/TopTagger/test
./configure 
make -j4
```

Now also check out our repository if not done already:
```
cd $CMSSW_BASE/src
git clone git@github.com:StealthStop/Framework.git
git clone -b Stealth git@github.com:susy2015/TopTaggerTools.git
git clone git@github.com:susy2015/NTupleReader.git
git clone git@github.com:StealthStop/Analyzer.git
cd Analyzer/Analyzer/test
source setup.sh #.csh if in tcsh
./configure
make -j4
```

We set up separate top tagger cfg files for each year, because different b tagger working points (WPs) are used.
Now we switch to DeepFlavor for b-tagging and top tagger also implements the medium DeepFlavor working point for each year.
So, we have two sets of releases, implementing DeepFlavor WPs.
The first set is the normal version, where the merged and resolved top tagger working points are passed to the tagger at run time.

A second set of configs is available.
Note that in the second set, we set up both DeepResolved and DeepAK8 WPs as 0.00 to get any top candidates in between WP 0 and 1 for calculating the SFs.

In the Framework/Framework/include/RunTopTagger.h, the resolved and merged working points are explicity applied when counting the number of tops.
```
cmsenv
getTaggerCfg.sh -t StealthStop_DeepFlavorWp0.2598_DeepResolvedwp0.95_DeepAK8wp0.937_2016preVFPUL -f TopTaggerCfg_2016preVFP.cfg -o
getTaggerCfg.sh -t StealthStop_DeepFlavorWp0.2489_DeepResolvedwp0.95_DeepAK8wp0.937_2016postVFPUL -f TopTaggerCfg_2016postVFP.cfg -o
getTaggerCfg.sh -t StealthStop_DeepFlavorWp0.3040_DeepResolvedwp0.95_DeepAK8wp0.895_2017UL -f TopTaggerCfg_2017.cfg -o
getTaggerCfg.sh -t StealthStop_DeepFlavorWp0.2783_DeepResolvedwp0.95_DeepAK8wp0.895_2018UL -f TopTaggerCfg_2018.cfg -o
```

We have two set of Double DisCo NN, one for RPV model and another one SYY model.
To get all them, run these comments below. Note that any relese with patch number 1 (e.g. v3.0.1) contains optimized bin edges whereas patch number 0 (e.g. v3.0.0) has non-optimized bin edges.

```
cmsenv
getDeepESMCfg.sh -t DoubleDisCo_Reg_0l_Run2_RPV_v3.4.1_MassExclusion -o -m DoubleDisCo_Reg.cfg -M DoubleDisCo_Reg_NonIsoMuon.cfg -f Keras_Tensorflow -F Keras_Tensorflow_NonIsoMuon -s DoubleDisCo_Reg_0l_RPV_Run2_MassExclusion
getDeepESMCfg.sh -t DoubleDisCo_Reg_0l_Run2_SYY_v3.4.1_MassExclusion -o -m DoubleDisCo_Reg.cfg -M DoubleDisCo_Reg_NonIsoMuon.cfg -f Keras_Tensorflow -F Keras_Tensorflow_NonIsoMuon -s DoubleDisCo_Reg_0l_SYY_Run2_MassExclusion
getDeepESMCfg.sh -t DoubleDisCo_Reg_1l_Run2_RPV_v3.4.1_MassExclusion -o -m DoubleDisCo_Reg.cfg -M DoubleDisCo_Reg_NonIsoMuon.cfg -f Keras_Tensorflow -F Keras_Tensorflow_NonIsoMuon -s DoubleDisCo_Reg_1l_RPV_Run2_MassExclusion
getDeepESMCfg.sh -t DoubleDisCo_Reg_1l_Run2_SYY_v3.4.1_MassExclusion -o -m DoubleDisCo_Reg.cfg -M DoubleDisCo_Reg_NonIsoMuon.cfg -f Keras_Tensorflow -F Keras_Tensorflow_NonIsoMuon -s DoubleDisCo_Reg_1l_SYY_Run2_MassExclusion
getDeepESMCfg.sh -t DoubleDisCo_Reg_2l_Run2_RPV_v3.4.1_MassExclusion -o -m DoubleDisCo_Reg.cfg -M DoubleDisCo_Reg_NonIsoMuon.cfg -f Keras_Tensorflow -F Keras_Tensorflow_NonIsoMuon -s DoubleDisCo_Reg_2l_RPV_Run2_MassExclusion
getDeepESMCfg.sh -t DoubleDisCo_Reg_2l_Run2_SYY_v3.5.1_MassExclusion -o -m DoubleDisCo_Reg.cfg -M DoubleDisCo_Reg_NonIsoMuon.cfg -f Keras_Tensorflow -F Keras_Tensorflow_NonIsoMuon -s DoubleDisCo_Reg_2l_SYY_Run2_MassExclusion
getDeepESMCfg.sh -t DoubleDisCo_Reg_0l_Run2_RPV_v3.4.1_MaxSig -o -m DoubleDisCo_Reg.cfg -M DoubleDisCo_Reg_NonIsoMuon.cfg -f Keras_Tensorflow -F Keras_Tensorflow_NonIsoMuon -s DoubleDisCo_Reg_0l_RPV_Run2_MaxSig
getDeepESMCfg.sh -t DoubleDisCo_Reg_0l_Run2_SYY_v3.4.1_MaxSig -o -m DoubleDisCo_Reg.cfg -M DoubleDisCo_Reg_NonIsoMuon.cfg -f Keras_Tensorflow -F Keras_Tensorflow_NonIsoMuon -s DoubleDisCo_Reg_0l_SYY_Run2_MaxSig
getDeepESMCfg.sh -t DoubleDisCo_Reg_1l_Run2_RPV_v3.4.1_MaxSig -o -m DoubleDisCo_Reg.cfg -M DoubleDisCo_Reg_NonIsoMuon.cfg -f Keras_Tensorflow -F Keras_Tensorflow_NonIsoMuon -s DoubleDisCo_Reg_1l_RPV_Run2_MaxSig
getDeepESMCfg.sh -t DoubleDisCo_Reg_1l_Run2_SYY_v3.4.1_MaxSig -o -m DoubleDisCo_Reg.cfg -M DoubleDisCo_Reg_NonIsoMuon.cfg -f Keras_Tensorflow -F Keras_Tensorflow_NonIsoMuon -s DoubleDisCo_Reg_1l_SYY_Run2_MaxSig
getDeepESMCfg.sh -t DoubleDisCo_Reg_2l_Run2_RPV_v3.5.1_MaxSig -o -m DoubleDisCo_Reg.cfg -M DoubleDisCo_Reg_NonIsoMuon.cfg -f Keras_Tensorflow -F Keras_Tensorflow_NonIsoMuon -s DoubleDisCo_Reg_2l_RPV_Run2_MaxSig
getDeepESMCfg.sh -t DoubleDisCo_Reg_2l_Run2_SYY_v3.5.1_MaxSig -o -m DoubleDisCo_Reg.cfg -M DoubleDisCo_Reg_NonIsoMuon.cfg -f Keras_Tensorflow -F Keras_Tensorflow_NonIsoMuon -s DoubleDisCo_Reg_2l_SYY_Run2_MaxSig
```

## Running an Analyzer Locally

Analyzer modules are run via the main "controlling" script `MyAnalysis.C`, which handles several different arguments.

```
Options:
   -c : runOnCondor  An internally specified argument which is used to signify when running on a condor cluster node
   -v : isQuiet      When specified, print logging information while running (custom to a specific analyzer) 
   -s : fastMode     Run the analyzer in fast mode, where the module pipeline can be terminated early
   -A : analyzer     The name of the analyzer to be executed
   -H : histFile     Name of the output ROOT file to contain any histograms created by the analyzer
   -D : dataSets     Comma-separated list of data set names to run over
   -N : nFiles       Number of files (per data set) to process
   -M : startFile    Which file in the data set filelist to start processing at
   -E : maxEvts      Absolute maximum number of events to process for each data set
```

An example of running MyAnalysis interactively is
```
cd $CMSSW_BASE/src/Analyzer/Analyzer/test/
./MyAnalysis -A AnalyzeDoubleDisCo -H myoutputfile.root -D 2017_TTToSemiLeptonic -E 1001 -M 2 -s 
```

## Condor submission

The `condor` subdirectory contains some scripts to help submit jobs via condor on the cmslpc cluster. 

The requirements for condor submission are: 
```
 1. A shell script to run on the worker node. This script should set up the working area, copy any needed files, call `MyAnalysis.C` with the right options, and make sure the output gets copied to the user's EOS area.
     - The example included here is [run_Analyzer_condor.sh](Analyzer/test/condor/run_Analyzer_condor.sh)
 2. One or more tarballs to unpack on the worker node, these usually contain a slimmed down CMSSW area, and the `MyAnalysis` executable with any needed libraries
 3. A so-called jdl file that contains the condor setup and specifies the jobs to be submitted
     - The last two items are produced by a python script called [condorSubmit.py](Analyzer/test/condor/condorSubmit.py). 
```

An example call to the condor submission script would be:

```
python condorSubmit.py --analyze AnalyzeDoubleDisCo --output DisCoAnaOutput -d "2016preVFP_TT,2016preVFP_QCD" -n 20
```
where possible arguments are

```
Usage: condorSubmit.py [options]

Options:
  -h, --help         show this help message and exit
  -n NUMFILE         number of files per job
  -d DATASETS        List of datasets, comma separated
  -l                 List all datacollections
  -L                 List all datacollections and sub collections
  -c                 Do not submit jobs.  Only create condor_submit.txt.
  -s                 Run Analyzer in fast mode
  -u USEROVERRIDE    Override username with something else
  --output=OUTPATH   Name of directory where output of each condor job goes
  --analyze=ANALYZE  AnalyzeBackground, AnalyzeEventSelection, Analyze0Lep,
                     Analyze1Lep, MakeNJetDists
```

With the `-n` option one can specify how many files to run over per job.
The `--analyze` option lets the user pick which analyzer to use. 
MyAnalysis incorporates dedicator code modules to keep track of datasets, their cross sections, and their names. 
To see a list of available datasets, one can call the submission script with the `-l` or `-L` options.
Pass the list of datasets desired to run over to the script with the option `-d`. 
Before submitting jobs, make sure to have called `voms-proxy-init`. 

In the event that jobs fail and do not send their output to EOS, a script is provided that can resubmit these missing jobs.
All that is required is the original job folder in `condor` as well as the corresponding output folder in EOS.
The cleanup script compares the total possible jobs (via `.log` files) to the output `.root` files to identify jobs that did not complete correctly.
Given the example call to `condorSubmit.py` up above, a corresponding call to the cleanup script would be:

```
python condorSubmit.py --analyze AnalyzeDoubleDisCo --jobdir DisCoAnaOutput
```
with available options given as

```
Usage: cleanupSubmit.py [options]

Options:
  -h, --help         show this help message and exit
  -c                 Do not submit jobs.  Only create condor_submit.txt.
  -s                 Run Analyzer in fast mode
  -u USEROVERRIDE    Override username with something else
  --jobdir=JOBDIR    Name of directory where output of each condor job goes
  --analyze=ANALYZE  AnalyzeBackground, AnalyzeEventSelection, Analyze0Lep,
                     Analyze1Lep, MakeNJetDists
```

## Plotting from TTrees

A script that wraps around the `TTree->Draw()` concept is provided in the form of `miniTupleDrawer.py`.
This gives easy abilities to plot from TTrees produced from an analyzer derived from the `MiniTupleMaker` class.
Current analyzers that produce simple TTrees are `MakeMiniTree`, `MakeNNVariables`, and `MakeQCDValTree`.
The TTree drawer script requires a "sidecar" auxiliary file that specifies a dictionary of histogram names mapped to a subdictionary of options.
An example aux file is `miniTupleDrawer_aux.py`.
```
usage: %miniTupleDrawer [options] [-h] --inputDir INPUTDIR
                                  [--outputDir OUTPUTDIR] [--tree TREE]
                                  [--year YEAR] [--options OPTIONS]

optional arguments:
  -h, --help            show this help message and exit
  --inputDir INPUTDIR   Path to ntuples
  --outputDir OUTPUTDIR
                        path to output
  --tree TREE           TTree name to draw
  --year YEAR           which year
  --options OPTIONS     options file
```
An example call to this script would be:
```
python Plotters/General/miniTupleDrawer.py --options miniTupleDrawer_aux \
                                           --inputDir ~/path/to/minituples/ \
                                           --outputDir subdir/structure/in/condor/folder \
                                           --tree PreSelection \
                                           --year Run2UL
```
Output ROOT files with the drawn histograms are contained in the `outputDir` folder subdirectory structure and placed automatically in the `condor` folder.
This placement of the output makes it intuitive to then use the other plotting tools (`stackPlotter` mentioned below) for making final, pretty plots.

## Making stack plots

A generic plotter is available for making stack plots with or without a data/MC ratio panel.
Information on which backgrounds and signals to plot in addition to visualization should be specified
in the `stackPlotter_aux.py` file. There too, one specifies which histograms to extract from the ROOT
file and plot.

```
usage: usage: %stackPlotter [options] [-h] [--noRatio] [--approved]
                                      [--printNEvents] [--normMC2Data]
                                      [--normalize] [--printSign] --inpath
                                      INPATH --outpath OUTPATH --year YEAR
                                      [--options OPTIONS]

optional arguments:
  -h, --help         show this help message and exit
  --noRatio          No ratio plot
  --approved         Plot is approved
  --printNEvents     Show number of events
  --normMC2Data      Normalize MC to data
  --normalize        Normalize all to unity
  --printSign        Print simple significance
  --inpath INPATH    Path to root files
  --outpath OUTPATH  Where to put plots
  --year YEAR        which year
  --options OPTIONS  options file
usage: usage: %stackPlotter [options] [-h] [--noRatio] [--approved]
                                      [--printNEvents] [--normMC]
                                      [--printSign] --inpath INPATH --outpath
                                      OUTPATH --year YEAR
```

An example call to the stack plotter could be:

```
python stackPlotter.py --year 2016 --inpath ./condor/2016_DisCo_0L_1L_hadd/ --outpath plot_histos --normMC2Data
```

## Make Yields Tables

A tool is provided to load information from ROOT files output by the `AnalyzeDoubleDisCo` analyzer and make a table of event yields (and fractional yields) for different backgrounds and signal processes.
The user can choose which channel as well as toggle between QCD CR or not.
The tables are output in standalone `.tex` format that can be simply input into another `.tex` document.

```
usage: usage: %tableYields [options] [-h] --channel CHANNEL --inputDir
                                     INPUTDIR --outputDir OUTPUTDIR [--QCDCR]
                                     [--year YEAR]

optional arguments:
  -h, --help            show this help message and exit
  --channel CHANNEL     which channel to process
  --inputDir INPUTDIR   directory for input ROOT
  --outputDir OUTPUTDIR
                        where to put tex table files
  --QCDCR               do for QCD CR
  --year YEAR           which year to process
```

An example call to the tool would be:

```
python Tools/makeYieldsTables.py --inputDir ~/nobackup/outputsPath/DataVsMC_Run2 --outputDir MyOutput --channel 1l --year Run2
```


## Generating Filelists, Checking Ntuples, and Checking Event Numbers

### Produce the Filelist and sampleSet.cfg

The main script for generating file lists and the corresponding sample set is `Tools/makeFilelist.py`.

```
usage: makeFilelist.py [-h] [--prod PROD] [--tag TAG] [--skim]

optional arguments:
  -h, --help   show this help message and exit
  --prod PROD  Unique tag for output
  --tag TAG    Path to PU file
  --skim       Make for skim
```

The `--prod` argument refers to the versioning of the folder on EOS that contains ntuple ROOT files for all MC and data samples.
The latest version is `V20` and can be found at `/eos/uscms/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/`.

Running the script without any arguments will generate a folder `filelists_Kevin_UL_v2` and 
`sampleSet_UL_v2.cfg`, where the default `--tag="UL_v2"` has been used.

The filelist folder contains a text file for each MC sample (found in the directory of ntuple files), which lists the full paths to all the ntuple ROOT files
for the corresponding sample.
This folder is intended to be placed in the `StealthStop` area of the `lpcsusystealth` group space on EOS.

The `sampleSet_UL_v2.cfg` file contains a mapping between "friendly" sample names and the corresponding text file listing all ntuple ROOT files for the sample.
This config file should be placed in the `cfg` area of `Framework`.
In order to pick up this new `cfg` in the `Analyzer` area automatically, the appropriate lines in the `getSamplesCfg.sh` script in `Framework/scripts` need to be modified and `source setup.sh` rerun.

Additionally, a new `sampleCollection_UL_v2.cfg` needs to be constructed (easiest by hand), which creates groups of samples that are to be referenced when running analyzers.

Note, when running `makeFilelist.py` it is most effective to have an up-to-date `TreeMaker` to reference in the script.
This allows population of each sample line with total event numbers, cross sections, k factors.

When specifying the option `--skim`, the script is configured to look for a `Skims/{2016preVFP,2016postVFP,2017,2018}` folders in the `lpcsusystealth/StealthStop` area.
These folders contain ROOT files from the `MakeMiniTree` analyzer which contain a smaller size event TTree.
The filelist script handles these differences when looking for skimmed ROOT files, but returns the same sort of outputs as described above.

### Checking Positive and Negative Events

With a new `sampleSets.cfg` symlink in `Analyzer/test` pointing to `sampleSets_UL_v2.cfg` in `Framework`, the number of events in each sample can be measured.
This is useful for verifying exactly how many events are present for a given sample in order to guarantee accurate calculation of an event weight.

`nEvt.py` jobs can be submitted with `nEvtsCondorSubmit.py`, which will read in `sampleSets.cfg` and spawn a job for each sample and loop through all its files.
An output text file is generated in the end and returned to the user which reads total positive and negative events counts for the sample.
Available arguments for the script are provided below.

```
Usage: nEvtsCondorSubmit.py [options]


Options:
  -h, --help            show this help message and exit
  --noSubmit            Do not submit jobs. Only create condor_submit.txt
  --outPath=OUTPATH     Name of directory where output of each condor job goes
  --sampleSets=SAMPLESETS
                        Sample sets config file
  --wildcard=WILDCARD   Wildcard expression for picking only some sample sets
```

An example call to the script would be
```
python nEvtsCondorSubmit.py --outPath nEvtsOutput --sampleSet sampleSets --wildcard "*2016pre*"
```

A helper script `checkNevents.py` is available to compare the numbers reported in the `nEvt.py` job output and the original `sampleSets.cfg`, 
whose numbers were sourced directly from the `TreeMaker` repository.
Discrepancies are printed to screen for investigation.
Additionally, a new `sampleSets_new.cfg` is written with the numbers measured by `nEvt.py` inserted into the original `sampleSets.cfg`.

```
usage: checkNevents.py [-h] [--sampleSet SAMPLESET] [--nEvtsDir NEVTDIR]

optional arguments:
  -h, --help            show this help message and exit
  --sampleSet SAMPLESET
                        Path to sample set file
  --nEvtsDir NEVTDIR     Directory to nEvt output
```

## Splitting Signal Files

The UL signal samples are produced where all mass points from 300 to 1400 can appear in a single ROOT file.
To restore the behavior from the legacy analysis, where a given ROOT file only contains events for a single mass point
some code infrastructure is provided to disentangle the mass points.

Jobs can be submitted to condor to run over ROOT files in a user-specified directory on EOS.
The output ROOT files are sent to a user-specified area in their own EOS area.

```
usage: submitSplitSignal.py [-h] --eosPath EOSPATH [--outPath OUTPATH]
                            [--era ERA] [--model MODEL]
                            [--ttreePath TTREEPATH] [--noSubmit]

optional arguments:
  -h, --help            show this help message and exit
  --eosPath EOSPATH     Path to files on EOS
  --outPath OUTPATH     Output path for jobs
  --era ERA             Era for signal samples
  --model MODEL         Signal model to split
  --ttreePath TTREEPATH
                        TTree name to read
  --noSubmit            Do not submit to condor
```

Companion `runSplitSignal.py` and `runSplitSignal.sh` are provided to process each file.
A full suite of calls to this script to split our signals would be:

```
python submitSplitSignal.py --eosPath /eos/uscms/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/ --outPath SusyRA2Analysis2015/Run2ProductionV20/ --era Summer20UL16 --model RPV
python submitSplitSignal.py --eosPath /eos/uscms/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/ --outPath SusyRA2Analysis2015/Run2ProductionV20/ --era Summer20UL16APV --model RPV
python submitSplitSignal.py --eosPath /eos/uscms/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/ --outPath SusyRA2Analysis2015/Run2ProductionV20/ --era Summer20UL17 --model RPV
python submitSplitSignal.py --eosPath /eos/uscms/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/ --outPath SusyRA2Analysis2015/Run2ProductionV20/ --era Summer20UL18 --model RPV

python submitSplitSignal.py --eosPath /eos/uscms/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/ --outPath SusyRA2Analysis2015/Run2ProductionV20/ --era Summer20UL16 --model StealthSYY
python submitSplitSignal.py --eosPath /eos/uscms/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/ --outPath SusyRA2Analysis2015/Run2ProductionV20/ --era Summer20UL16APV --model StealthSYY
python submitSplitSignal.py --eosPath /eos/uscms/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/ --outPath SusyRA2Analysis2015/Run2ProductionV20/ --era Summer20UL17 --model StealthSYY
python submitSplitSignal.py --eosPath /eos/uscms/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/ --outPath SusyRA2Analysis2015/Run2ProductionV20/ --era Summer20UL18 --model StealthSYY

python submitSplitSignal.py --eosPath /eos/uscms/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/ --outPath SusyRA2Analysis2015/Run2ProductionV20/ --era Summer20UL16 --model StealthSHH
python submitSplitSignal.py --eosPath /eos/uscms/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/ --outPath SusyRA2Analysis2015/Run2ProductionV20/ --era Summer20UL16APV --model StealthSHH
python submitSplitSignal.py --eosPath /eos/uscms/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/ --outPath SusyRA2Analysis2015/Run2ProductionV20/ --era Summer20UL17 --model StealthSHH
python submitSplitSignal.py --eosPath /eos/uscms/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/ --outPath SusyRA2Analysis2015/Run2ProductionV20/ --era Summer20UL18 --model StealthSHH
```

Here the split signals are automatically returned to the `lpcsusyhad` EOS space in the `SusyRA2Analysis2015/Run2ProductionV20/` subdirectory.
