#! /bin/env/python
import copy

regions = {"0l_RPV_A"        : "DoubleDisCo_disc2_0l_RPV>0.6&&DoubleDisCo_disc1_0l_RPV>0.6",
           "0l_StealthSYY_A" : "DoubleDisCo_disc2_0l_SYY>0.6&&DoubleDisCo_disc1_0l_SYY>0.6",
           "1l_RPV_A"        : "DoubleDisCo_disc2_1l_RPV>0.6&&DoubleDisCo_disc1_1l_RPV>0.6",
           "1l_StealthSYY_A" : "DoubleDisCo_disc2_1l_SYY>0.6&&DoubleDisCo_disc1_1l_SYY>0.6",

           "0l_RPV_B"        : "DoubleDisCo_disc2_0l_RPV>0.6&&DoubleDisCo_disc1_0l_RPV<0.6",
           "0l_StealthSYY_B" : "DoubleDisCo_disc2_0l_SYY>0.6&&DoubleDisCo_disc1_0l_SYY<0.6",
           "1l_RPV_B"        : "DoubleDisCo_disc2_1l_RPV>0.6&&DoubleDisCo_disc1_1l_RPV<0.6",
           "1l_StealthSYY_B" : "DoubleDisCo_disc2_1l_SYY>0.6&&DoubleDisCo_disc1_1l_SYY<0.6",

           "0l_RPV_C"        : "DoubleDisCo_disc2_0l_RPV<0.6&&DoubleDisCo_disc1_0l_RPV>0.6",
           "0l_StealthSYY_C" : "DoubleDisCo_disc2_0l_SYY<0.6&&DoubleDisCo_disc1_0l_SYY>0.6",
           "1l_RPV_C"        : "DoubleDisCo_disc2_1l_RPV<0.6&&DoubleDisCo_disc1_1l_RPV>0.6",
           "1l_StealthSYY_C" : "DoubleDisCo_disc2_1l_SYY<0.6&&DoubleDisCo_disc1_1l_SYY>0.6",

           "0l_RPV_D"        : "DoubleDisCo_disc2_0l_RPV<0.6&&DoubleDisCo_disc1_0l_RPV<0.6",
           "0l_StealthSYY_D" : "DoubleDisCo_disc2_0l_SYY<0.6&&DoubleDisCo_disc1_0l_SYY<0.6",
           "1l_RPV_D"        : "DoubleDisCo_disc2_1l_RPV<0.6&&DoubleDisCo_disc1_1l_RPV<0.6",
           "1l_StealthSYY_D" : "DoubleDisCo_disc2_1l_SYY<0.6&&DoubleDisCo_disc1_1l_SYY<0.6",
}

nbjets = {"Nbjets2incl" : "NGoodBJets_pt30>=2",
}

njets = {"Njets7"      : "NGoodJets_pt30==7",
         "Njets7incl"  : "NGoodJets_pt30>=7",
         "Njets8"      : "NGoodJets_pt30==8",
         "Njets8incl"  : "NGoodJets_pt30>=8",
         "Njets9"      : "NGoodJets_pt30==9",
         "Njets9incl"  : "NGoodJets_pt30>=9",
         "Njets10"     : "NGoodJets_pt30==10",
         "Njets10incl" : "NGoodJets_pt30>=10",
         "Njets11"     : "NGoodJets_pt30==11",
         "Njets11incl" : "NGoodJets_pt30>=11",
         "Njets12"     : "NGoodJets_pt30==12",
         "Njets12incl" : "NGoodJets_pt30>=12",
         "Njets13"     : "NGoodJets_pt30==13",
         "Njets13incl" : "NGoodJets_pt30>=13",
}

histos = {"h_DoubleDisCo_RPV_disc1_disc2_0l" : {"selection" : "${SELECTION}", "variable" : "DoubleDisCo_disc2_0l_RPV:DoubleDisCo_disc1_0l_RPV", "xbins" : 100, "xmin" : 0, "xmax" : 1, "ybins" : 100, "ymin" : 0, "ymax" : 1  }, 
          "h_DoubleDisCo_SYY_disc1_disc2_0l" : {"selection" : "${SELECTION}", "variable" : "DoubleDisCo_disc2_0l_SYY:DoubleDisCo_disc1_0l_SYY", "xbins" : 100, "xmin" : 0, "xmax" : 1, "ybins" : 100, "ymin" : 0, "ymax" : 1  }, 
          "h_DoubleDisCo_RPV_disc1_disc2_1l" : {"selection" : "${SELECTION}", "variable" : "DoubleDisCo_disc2_1l_RPV:DoubleDisCo_disc1_1l_RPV", "xbins" : 100, "xmin" : 0, "xmax" : 1, "ybins" : 100, "ymin" : 0, "ymax" : 1  },
          "h_DoubleDisCo_SYY_disc1_disc2_1l" : {"selection" : "${SELECTION}", "variable" : "DoubleDisCo_disc2_1l_SYY:DoubleDisCo_disc1_1l_SYY", "xbins" : 100, "xmin" : 0, "xmax" : 1, "ybins" : 100, "ymin" : 0, "ymax" : 1  },
          "h_Mbb"                            : {"selection" : "${SELECTION}", "variable" : "Mbb",                                               "xbins" : 100, "xmin" : 0, "xmax" : 500                                       }, 
          "h_dRbjets"                        : {"selection" : "${SELECTION}", "variable" : "dR_bjets",                                          "xbins" : 50,  "xmin" : 0, "xmax" : 5                                         },
          "h_dRbjets_Mbb"                    : {"selection" : "${SELECTION}", "variable" : "Mbb:dR_bjets",                                      "xbins" : 100, "xmin" : 0, "xmax" : 5, "ybins" : 100, "ymin" : 0, "ymax" : 500}, 
}

histograms = {}

for njetStr, njetExp in njets.items():
    for nbjetStr, nbjetExp in nbjets.items():
        for regionStr, regionExp in regions.items():
            cutStr = "&&".join([njetExp, nbjetExp, regionExp])
            for histoName, histoOps in histos.items():
                hopsCopy = copy.copy(histoOps)
                hopsCopy["selection"] = "${WEIGHT}*(%s)"%("".join(cutStr))
                histograms["%s_%s_%s_%s"%(histoName,regionStr,njetStr,nbjetStr)] = hopsCopy

processes = [
    "TT",
    "QCD",
    "Non_QCD",
    "Data_SingleMuon"
]

for mass in [550]:
    processes.append("RPV_2t6j_mStop-%d"%(mass))
    processes.append("StealthSYY_2t6j_mStop-%d"%(mass))
