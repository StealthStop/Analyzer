# README for Data driven systematic from the QCD Control Region

## Last updated by Kelvin Mei on January 12th, 2021

The purpose of the scripts in this folder is to create a systematic from looking at data in the QCD control region (CR) and the Powheg TT MC in the signal region (SR) for the SUS-19-004 analysis. 

The general idea is that the trends in high njets in the Snn distribution of the NN trainings are similar for QCD in the CR, data in the CR and Powheg TT MC in the SR, and so it is reasonable to expect similar trends in the TT dominated SR.

Thus, an additional multiplicative weight can be derived from data in the CR and from MC in both the SR and CR The ratio of the Snn shapes between the data in the CR and the TT MC in the SR can be the backbone of a new systematic, giving an indirect measurement of how the NN impacts the data over the MC differently.

To run these scripts, you will first need:

* The mini tree that has the totalEventWeight (all SFs applied) for the TT MC that pass the SR baseline
* The mini trees that have the data (Single Muon data set) and the QCD MC that pass the CR baseline
* The njets\_for\_Aron.root file that has the final TT njets shape (this can be made from the mini trees if you do not have it).
* (Optional) The mini tree that has the QCD MC that pass the SR baseline and TT MC for events that pass the CR baseline. The code currently is set up so that it can run with these files, and it produces plots that compare the Snn distributions between these samples and the other MC samples/CR data, but it is not strictly necessary.

The order of the scripts are as follows:

1. `python makeSnnHistograms.py --signal`
2. `python makeSnnHistograms.py --signal --ttbar`
3. `python makeSnnCompPlots.py`
4. `python makeRatioOfRatios.py`
5. `python throwToysForSystematic.py`
6. `python makeQCDCRSystematic.py`
