#! /bin/env/python
import copy

selections = {"pass0Lbaseline" : "passBaseline0l_Good",
}

histosInfo = {"h_NMrgTops_withTopTagSF" : {"weight" : "TotalWeight_0l",                      "selection" : "${SELECTION}", "variable" : "ntops_1jet",            "xbins" : 6, "xmin" : -0.5, "xmax" : 5.5}, 
              "h_NMrgTops_noTopTagSF"   : {"weight" : "TotalWeight_0l/topTaggerScaleFactor", "selection" : "${SELECTION}", "variable" : "ntops_1jet",            "xbins" : 6, "xmin" : -0.5, "xmax" : 5.5}, 
              "h_NResTops_withTopTagSF" : {"weight" : "TotalWeight_0l",                      "selection" : "${SELECTION}", "variable" : "ntops_3jet",            "xbins" : 6, "xmin" : -0.5, "xmax" : 5.5}, 
              "h_NResTops_noTopTagSF"   : {"weight" : "TotalWeight_0l/topTaggerScaleFactor", "selection" : "${SELECTION}", "variable" : "ntops_3jet",            "xbins" : 6, "xmin" : -0.5, "xmax" : 5.5}, 
              "h_NTops_withTopTagSF"    : {"weight" : "TotalWeight_0l",                      "selection" : "${SELECTION}", "variable" : "ntops_1jet+ntops_3jet", "xbins" : 6, "xmin" : -0.5, "xmax" : 5.5}, 
              "h_NTops_noTopTagSF"      : {"weight" : "TotalWeight_0l/topTaggerScaleFactor", "selection" : "${SELECTION}", "variable" : "ntops_1jet+ntops_3jet", "xbins" : 6, "xmin" : -0.5, "xmax" : 5.5}, 
}

histograms = {}

for aStr, aExp in selections.items():
    for histoName, histoOps in histosInfo.items():
        hopsCopy = copy.copy(histoOps)
        hopsCopy["selection"] = aExp
        histograms["%s_%s"%(histoName,aStr)] = hopsCopy

processes = [
    "TT",
    "QCD",
    "TTX",
    "BG_OTHER",
    "Data"
]
