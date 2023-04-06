#! /bin/env/python
import copy

njetsSelections = {"Njets7incl" : "NGoodJets_pt30>=7&&passBaseline1l_Good",
}

histosInfo = {"h_DoubleDisCo_RPV_disc1_disc2_1l" : {"weight" : "TotalWeight_1l", "selection" : "${SELECTION}", "variable" : "DoubleDisCo_disc2_1l_RPV:DoubleDisCo_disc1_1l_RPV", "xbins" : 100, "xmin" : 0, "xmax" : 1, "ybins" : 100, "ymin" : 0, "ymax" : 1}, 
}

histograms = {}

for njetStr, njetExp in njetsSelections.items():
    for histoName, histoOps in histosInfo.items():
        hopsCopy = copy.copy(histoOps)
        hopsCopy["selection"] = njetExp
        histograms["%s_%s"%(histoName,njetStr)] = hopsCopy

processes = [
    "TT",
    "QCD",
    "TTX",
    "BG_OTHER",
    "Data"
]

for mass in [550]:
    processes.append("RPV_2t6j_mStop-%d"%(mass))
    processes.append("StealthSYY_2t6j_mStop-%d"%(mass))
