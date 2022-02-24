# Analyzer
Makes input histograms for our combine setup

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
git clone -b Run2_UL git@github.com:StealthStop/Framework.git
git clone -b Stealth git@github.com:susy2015/TopTaggerTools.git
git clone git@github.com:susy2015/NTupleReader.git
git clone -b Run2_UL git@github.com:StealthStop/Analyzer.git
cd Analyzer/Analyzer/test
source setup.sh #.csh if in tcsh
./configure
make -j4
```

We set up the top tagger cfg files for per year, because per year has different b-tagger working points (WPs).
Last step is to get the cfg and model files for the top tagger, deepESM, and mass regression.
```
cmsenv
getTaggerCfg.sh -t StealthStop_DeepCSV_DeepResolved_DeepAK8_wp0.98_2016_v1 -f TopTaggerCfg_2016.cfg -o
getTaggerCfg.sh -t StealthStop_DeepCSV_DeepResolved_DeepAK8_wp0.98_2017_v1 -f TopTaggerCfg_2017.cfg -o
getTaggerCfg.sh -t StealthStop_DeepCSV_DeepResolved_DeepAK8_wp0.98_2018_v1 -f TopTaggerCfg_2018.cfg -o
getDeepESMCfg.sh -t Keras_Tensorflow_2016_v1.2 -o -s 2016
getDeepESMCfg.sh -t Keras_Tensorflow_2017_v1.2 -o -s 2017
getDeepESMCfg.sh -t Keras_Tensorflow_2018pre_v1.2 -o -s 2018pre
getDeepESMCfg.sh -t Keras_Tensorflow_2018post_v1.2 -o -s 2018post
getDeepESMCfg.sh -t DoubleDisCo_Reg_0l_RPV_2016_v3.0 -o -m DoubleDisCo_Reg.cfg -M DoubleDisCo_Reg_NonIsoMuon.cfg -f Keras_Tensorflow -F Keras_Tensorflow_NonIsoMuon -s DoubleDisCo_Reg_0l_RPV_2016
getDeepESMCfg.sh -t DoubleDisCo_Reg_1l_RPV_2016_v4.0 -o -m DoubleDisCo_Reg.cfg -M DoubleDisCo_Reg_NonIsoMuon.cfg -f Keras_Tensorflow -F Keras_Tensorflow_NonIsoMuon -s DoubleDisCo_Reg_1l_RPV_2016
```

Example of running MyAnalysis interactively
```
cd $CMSSW_BASE/src/Analyzer/Analyzer/test/
./MyAnalysis -A AnalyzeTest -H myoutputfile.root -D 2018post_TTToSemiLeptonic -E 1001 
```


## Condor submission

The condor directory contains some scripts to help submit jobs via condor on the cmslpc cluster. 
The requirements for condor submission are: 
 - A script to run on the worker node. This script should set up the area, copy any needed files, call your executable with the right options, and make sure the output gets copied to where you want. The example included here is [run_Analyzer_condor.tcsh](Analyzer/test/condor/run_Analyzer_condor.tcsh)
 - One or more tarballs to unpack on the worker node, these usually contain a slimmed down CMSSW area, and your executable with any needed libraries
 - A so-called jdl file that contains the condor setup and specifies the jobs to be submitted
The last two items are produced by a python script called [condorSubmit.py](Analyzer/test/condor/condorSubmit.py). 

```
[condor]$ python condorSubmit.py -h
Usage: condorSubmit.py [options]


Options:
  -h, --help         show this help message and exit
  -n NUMFILE         number of files per job
  -d DATASETS        List of datasets, comma separated
  -l                 List all datacollections
  -L                 List all datacollections and sub collections
  -c                 Do not submit jobs.  Only create condor_submit.txt.
  --output=OUTPATH   Name of directory where output of each condor job goes
  --analyze=ANALYZE  AnalyzeBackground (b), AnalyzeEventSelection (s),
                     Analyze0Lep (z), Analyze1Lep (o), MakeNJetDists (n)
```
As you can see from the above help menu, there are a number of options. 
With the `-n` option you can specify how many files to run over per job. The `--analyze` option lets you pick which analyzer to use. 
The MyAnalysis program has been updated to have these same switches. 
MyAnalysis now also uses the samples code to keep track of datasets, their cross sections, 6nd their names. 
To see a list of available datasets, you can call the submission script with the `-l` or `-L` options. Pass the list of datasets you want to run over to the script with the option `-d`. 
Before submitting jobs, make sure to have called `voms-proxy-init`. 

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

## Generating Filelists, Checking Ntuples, and Checking Event Numbers

### Produce the Filelist and sampleSet.cfg

The main script for generating file lists and the corresponding sample set is `makefilelist.py`.

```
usage: makefilelist.py [-h] [--prod PROD] [--tag TAG]

optional arguments:
  -h, --help   show this help message and exit
  --prod PROD  Unique tag for output
  --tag TAG    Path to PU file
```

The `prod` argument refers to the versioning of the folder on EOS that contains ntuple ROOT files for all MC and data samples.
The latest version is `V20` and can be found at `/eos/uscms/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/`.

Running the script without any arguments will generate a folder `filelists_Kevin_UL_v1` and 
`sampleSet_UL_v1.cfg`, where the default `--tag="UL_v1"` has been used.

The filelist folder contains a text file for each MC sample (found in the directory of ntuple files), which lists the full paths to all the ntuple ROOT files
for the corresponding sample.
This folder is intended to be placed in the `StealthStop` area of the `lpcsusyhad` group space on EOS.

The `sampleSet_UL_v1.cfg` file contains a mapping between "friendly" sample names and the corresponding text file listing all ntuple ROOT files for the sample.
This config file should be placed in the `cfg` area of `Framework`.
In order to pick up this new `cfg` in the `Analyzer` area automatically, the appropriate lines in the `getSamplesCfg.sh` script in `Framework/scripts` need to be modified and `source setup.sh` rerun.

Additionally, a new `sampleCollection_UL_v1.cfg` needs to be constructed (easiest by hand), which creates groups of samples that are to be referenced when running analyzers.

Note, when running `makefilelist.py` it is most effective to have an up-to-date `TreeMaker` to reference in the script.
This allows population of each sample line with total event numbers, cross sections, k factors.
These additional pieces of information are not technically necessary, but maintained in case of needing to use them.

### Checking Positive and Negative Events

With a new `sampleSet.cfg` symlink in `Analyzer/test` pointing to `sampleSet_UL_v1.cfg` in `Framework`, the number of events in each sample can be measured.
This is useful for verifying that all reported MINIAOD files in a sample were run on successfully by `TreeMaker` and hence that the event weight reflects the correct number of events.

`nEvt.py` jobs can be submitted with `nEvtsCondorSubmit.py`, which will read in `sampleSet.cfg` and spawn a job for each sample and loop through all its files.
An output text file is generated in the end and returned to the user which reads total positive and negative events counts for the sample.

```
Usage: nEvtsCondorSubmit.py [options]

Options:
  -h, --help         show this help message and exit
  -c                 Do not submit jobs. Only create condor_submit.txt
  --output=OUTPATH   Name of directory where output of each condor job goes
  -s SAMPLESETSFILE  Sample sets config file
```

A helper script `checkNevents.py` is available to compare the numbers reported in the `nEvt.py` job output and the original `sampleSet.cfg`, 
whose numbers were sourced directly from the `TreeMaker` repository.
Discrepancies are printed to screen for investigation.

```
usage: checkNevents.py [-h] [--sampleSet SAMPLESET] [--nEvtDir NEVTDIR]

optional arguments:
  -h, --help            show this help message and exit
  --sampleSet SAMPLESET
                        Path to sample set file
  --nEvtDir NEVTDIR     Directory to nEvt output
```

## Making inputs for the fit

Running the condor jobs to produce the input histograms for the fit.

```
cd $CMSSW_BASE/src/Analyzer/Analyzer/test/condor
python condorSubmit.py --analyze MakeNJetDists -d 2016_Data_SingleElectron,2016_Data_SingleMuon,2016_TT,2016_TT_erdOn,2016_TT_hdampUp,2016_TT_hdampDown,2016_TT_underlyingEvtUp,2016_TT_underlyingEvtDown,2016_TT_fsrUp,2016_TT_fsrDown,2016_TT_isrUp,2016_TT_isrDown,2016_WJets,2016_DYJetsToLL_M-50,2016_QCD,2016_ST,2016_Diboson,2016_Triboson,2016_TTX,2016_AllSignal -n 15 --output CondorOutput_2016_v1.2
python condorSubmit.py --analyze MakeNJetDists -d 2017_Data_SingleElectron,2017_Data_SingleMuon,2017_TT,2017_TT_erdOn,2017_TT_hdampUp,2017_TT_hdampDown,2017_TT_underlyingEvtUp,2017_TT_underlyingEvtDown,2017_WJets,2017_DYJetsToLL_M-50,2017_QCD,2017_ST,2017_Diboson,2017_Triboson,2017_TTX,2017_AllSignal                                                             -n 15 --output CondorOutput_2017_v1.2
python condorSubmit.py --analyze MakeNJetDists -d 2018pre_Data_SingleElectron,2018pre_Data_SingleMuon,2018pre_TT,2018pre_TT_erdOn,2018pre_TT_hdampUp,2018pre_TT_hdampDown,2018pre_TT_underlyingEvtUp,2018pre_TT_underlyingEvtDown,2018pre_WJets,2018pre_DYJetsToLL_M-50,2018pre_QCD,2018pre_ST,2018pre_Diboson,2018pre_Triboson,2018pre_TTX,2018pre_AllSignal                            -n 15 --output CondorOutput_2018pre_v1.2
python condorSubmit.py --analyze MakeNJetDists -d 2018post_Data_SingleElectron,2018post_Data_SingleMuon,2018post_TT,2018post_TT_erdOn,2018post_TT_hdampUp,2018post_TT_hdampDown,2018post_TT_underlyingEvtUp,2018post_TT_underlyingEvtDown,2018post_WJets,2018post_DYJetsToLL_M-50,2018post_QCD,2018post_ST,2018post_Diboson,2018post_Triboson,2018post_TTX,2018post_AllSignal                 -n 15 --output CondorOutput_2018post_v1.2
```

Now hadd the outputs when the jobs are done.

```
cd $CMSSW_BASE/src/Analyzer/Analyzer/test/condor
python hadder.py -d  2016_Data_SingleElectron,2016_Data_SingleMuon,2016_TT,2016_TT_erdOn,2016_TT_hdampUp,2016_TT_hdampDown,2016_TT_underlyingEvtUp,2016_TT_underlyingEvtDown,2016_TT_fsrUp,2016_TT_fsrDown,2016_TT_isrUp,2016_TT_isrDown,2016_WJets,2016_DYJetsToLL_M-50,2016_QCD,2016_ST,2016_Diboson,2016_Triboson,2016_TTX,2016_AllSignal  -H MakeNJetsDists_2016_v1.2     -p CondorOutput_2016_v1.2/output-files     -y 2016     --haddOther --haddData
python hadder.py -d  2017_Data_SingleElectron,2017_Data_SingleMuon,2017_TT,2017_TT_erdOn,2017_TT_hdampUp,2017_TT_hdampDown,2017_TT_underlyingEvtUp,2017_TT_underlyingEvtDown,2017_WJets,2017_DYJetsToLL_M-50,2017_QCD,2017_ST,2017_Diboson,2017_Triboson,2017_TTX,2017_AllSignal                                                              -H MakeNJetsDists_2017_v1.2     -p CondorOutput_2017_v1.2/output-files     -y 2017     --haddOther --haddData
python hadder.py -d  2018pre_Data_SingleElectron,2018pre_Data_SingleMuon,2018pre_TT,2018pre_TT_erdOn,2018pre_TT_hdampUp,2018pre_TT_hdampDown,2018pre_TT_underlyingEvtUp,2018pre_TT_underlyingEvtDown,2018pre_WJets,2018pre_DYJetsToLL_M-50,2018pre_QCD,2018pre_ST,2018pre_Diboson,2018pre_Triboson,2018pre_TTX,2018pre_AllSignal                             -H MakeNJetsDists_2018pre_v1.2  -p CondorOutput_2018pre_v1.2/output-files  -y 2018pre  --haddOther --haddData
python hadder.py -d  2018post_Data_SingleElectron,2018post_Data_SingleMuon,2018post_TT,2018post_TT_erdOn,2018post_TT_hdampUp,2018post_TT_hdampDown,2018post_TT_underlyingEvtUp,2018post_TT_underlyingEvtDown,2018post_WJets,2018post_DYJetsToLL_M-50,2018post_QCD,2018post_ST,2018post_Diboson,2018post_Triboson,2018post_TTX,2018post_AllSignal                  -H MakeNJetsDists_2018post_v1.2 -p CondorOutput_2018post_v1.2/output-files -y 2018post --haddOther --haddData
```

If there are missing jobs you should see a message in red after each sample is hadded.
After hadding now put all the samples together into one file for the fit input.

```
cd $CMSSW_BASE/src/Analyzer/Analyzer/test/
python write_fit_input.py -d condor/MakeNJetsDists_2016_v1.2     -H FitInput/Keras_2016_v1.2     -y 2016
python write_fit_input.py -d condor/MakeNJetsDists_2017_v1.2     -H FitInput/Keras_2017_v1.2     -y 2017
python write_fit_input.py -d condor/MakeNJetsDists_2018pre_v1.2  -H FitInput/Keras_2018pre_v1.2  -y 2018pre
python write_fit_input.py -d condor/MakeNJetsDists_2018post_v1.2 -H FitInput/Keras_2018post_v1.2 -y 2018post
```

## Deriving ttbar shape systematics

Before running `njets_systs_comp.py` one must run `run_fits4ttbar_systematics.sh` in the HiggsAnalysis-CombinedLimit repository. This will stash away results in directories called TAG_YEAR_SYSTEMATIC

Once these have been generated, one can run `njets_systs_comp.py` with three arguments:

```
--fittag TAG
--year YEAR
--fitdir Base directory containing TAG_YEAR_SYSTEMATIC results
```

An example workflow to run would be something like:

```
cd $HOME/../../CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/

./run_fits4ttbar_systematics.sh

# Let this run for several hours...
# After it is finished successfully

cd $HOME/../../CMSSW_9_3_3/src/Analyzer/Analyzer/test/

# For example running on 2018pre
python njets_systs_comp.py --year 2018pre --fitdir $HOME/../../src/HiggsAnalysis/CombinedLimit/ --fittag Approval_StatErrPlusFullDev_12JetFix --outputdir MYDIR
```
