# README for QCD Estimate Derivation

## Last updated by Kelvin Mei on December 30th, 2020

The purpose of the scripts in this folder is to generate the qcd estimate for the analysis (1 lepton) through the qcd control region. It takes the qcd MC in both the control region and signal region and calculates a transfer factor that is then applied to the data in the control region, with all the other background subtracted as predicted by the MC (i.e. TT and Other).

To run these scripts, you will first need:

* The mini tree that has totalEventWeight (all SFs applied) for the QCD and TT MC that pass the signal region baseline (also TT+X and Other if you want to do extra comparison studies)
* The mini tree that has totalEventWeightNIM (just CR SFs applied) for the QCD MC, TT MC, TT+X MC, and Other that pass the control region baseline
* The mini tree that has the data (Single Muon data set) in the QCD control region that pass the control region baseline
* The njets\_for\_Aron.root file just for QCD passing the signal region baseline (normally output as QCD.root) that has the default QCD MC njets shape and all the njets shapes with the systematics applied. This needs to be copied to this directory and be named 201\*_QCD_MC.root

The order of the scripts are as follows:

1. `python makeQcdCrComparisonHistos.py`
2. `python makeDataQcdOnly.py`
3. `python makeUpdatedQcdMakeNJetsRootFile.py`
