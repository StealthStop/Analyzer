#! /bin/env/python
import copy

edgesLim = {
    "RPV" : {
        "0l" : (0.74, 0.80),
        "1l" : (0.80, 0.72)
    },
    "StealthSYY" : {
        "0l" : (0.54, 0.56),
        "1l" : (0.68, 0.82)
    }
}

edgesSig = {
    "RPV" : {
        "0l" : (0.52, 0.54),
        "1l" : (0.84, 0.42)
    },
    "StealthSYY" : {
        "0l" : (0.76, 0.70),
        "1l" : (0.44, 0.42)
    }
}

channels = ["0l", "1l"]
models   = ["RPV", "StealthSYY"]
regops   = [("A", ">", ">"), ("B", "<", ">"), ("C", ">", "<"), ("D", "<", "<")]

regions = {}

for regop in regops:
    for channel in channels:
        for model in models:

            chanmod = "%s_%s"%(channel, model)

            disc1 = edgesSig[model][channel][0]
            disc2 = edgesSig[model][channel][1]

            regions["%s_%s"%(chanmod, regop[0])] = "DoubleDisCo_disc2_%s%s%f&&DoubleDisCo_disc1_%s%s%f"%(chanmod.replace("Stealth", ""), regop[2], disc2, chanmod.replace("Stealth", ""), regop[1], disc1)

nbjets = {"Nbjets2incl" : "NGoodBJets_pt30>=2",
}

njets = {}

njetVals = range(7, 14)
for njet in njetVals:
    njets["Njets%d"%(njet)]     = "NGoodJets_pt30==%d"%(njet)
    njets["Njets%dincl"%(njet)] = "NGoodJets_pt30>=%d"%(njet)

histos = {#"h_DoubleDisCo_RPV_disc1_disc2_0l" : {"weight" : "TotalWeight_QCDCR", "selection" : "${SELECTION}", "variable" : "DoubleDisCo_disc2_NonIsoMuon_0l_RPV:DoubleDisCo_disc1_NonIsoMuon_0l_RPV", "xbins" : 100, "xmin" : 0, "xmax" : 1, "ybins" : 100, "ymin" : 0, "ymax" : 1  }, 
          #"h_DoubleDisCo_SYY_disc1_disc2_0l" : {"weight" : "TotalWeight_QCDCR", "selection" : "${SELECTION}", "variable" : "DoubleDisCo_disc2_NonIsoMuon_0l_SYY:DoubleDisCo_disc1_NonIsoMuon_0l_SYY", "xbins" : 100, "xmin" : 0, "xmax" : 1, "ybins" : 100, "ymin" : 0, "ymax" : 1  }, 
          #"h_DoubleDisCo_RPV_disc1_disc2_1l" : {"weight" : "TotalWeight_QCDCR", "selection" : "${SELECTION}", "variable" : "DoubleDisCo_disc2_NonIsoMuon_1l_RPV:DoubleDisCo_disc1_NonIsoMuon_1l_RPV", "xbins" : 100, "xmin" : 0, "xmax" : 1, "ybins" : 100, "ymin" : 0, "ymax" : 1  },
          #"h_DoubleDisCo_SYY_disc1_disc2_1l" : {"weight" : "TotalWeight_QCDCR", "selection" : "${SELECTION}", "variable" : "DoubleDisCo_disc2_NonIsoMuon_1l_SYY:DoubleDisCo_disc1_NonIsoMuon_1l_SYY", "xbins" : 100, "xmin" : 0, "xmax" : 1, "ybins" : 100, "ymin" : 0, "ymax" : 1  },
          #"h_Mbb"                            : {"weight" : "TotalWeight_QCDCR", "selection" : "${SELECTION}", "variable" : "Mbb",                                               "xbins" : 100, "xmin" : 0, "xmax" : 500                                       }, 
          "h_dRbjets"                        : {"weight" : "TotalWeight_QCDCR", "selection" : "${SELECTION}", "variable" : "dR_bjets",                                          "xbins" : 50,  "xmin" : 0, "xmax" : 5                                         },
}

histograms = {}

for njetStr, njetExp in njets.items():
    for nbjetStr, nbjetExp in nbjets.items():
        for regionStr, regionExp in regions.items():
            cutStr = "&&".join([njetExp, nbjetExp, regionExp])
            for histoName, histoOps in histos.items():
                hopsCopy = copy.copy(histoOps)
                hopsCopy["selection"] = cutStr
                histograms["%s_%s_%s_%s"%(histoName,regionStr,njetStr,nbjetStr)] = hopsCopy

processes = [
    "TT",
    "QCD",
    "Non_QCD",
    "Data_SingleMuon"
]

for mass in [400]:
    processes.append("RPV_2t6j_mStop-%d"%(mass))
    processes.append("StealthSYY_2t6j_mStop-%d"%(mass))
