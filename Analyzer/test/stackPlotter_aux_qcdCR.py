#! /bin/env/python

import ROOT

from collections import OrderedDict

# ----------------
# 1l QCD CR for 0l
# ----------------
qcdCR1_list1 = ["ge7NonIsoMuonJet", "7NonIsoMuonJet", "8NonIsoMuonJet", "9NonIsoMuonJet", "10NonIsoMuonJet", "ge11NonIsoMuonJet"]
qcdCR1_list2 = [
            "qcdCR_0l_HT300_1NonIsoMuon_?",
]

qcdCR1_selections = []

for cr1_list2 in qcdCR1_list2:

    for cr1_list1 in qcdCR1_list1:

         qcdCR1_selections.append( cr1_list2.replace("?", cr1_list1) )

# -----------
# 0l Baseline
# -----------
temp_list1 = ["ge7j", "7j", "8j", "9j", "10j", "ge11j"]

temp_list2 = [
        "0l_HT500_0NonIsoMuon_?_ge2t",
        #"0l_HT500_0NonIsoMuon_?_ge2t_ge1b",
        #"0l_HT500_0NonIsoMuon_?_ge2t_ge2b",
        #"0l_HT500_0NonIsoMuon_?_ge2t_ge1dRbjets",
]

base_selections = []

for list2 in temp_list2:

    for list1 in temp_list1:

        base_selections.append( list2.replace("?", list1) )

# ---------
# variables
# ---------
histograms = {
    "h_DoubleDisCo_disc1_?"   : {"logY" : False, "Y" : {"title" : "A.U.",   "min" : 2e-3}, "X" : {"title" : "Disc 1",                "rebin" : 20, "min" :  0,    "max" :    1}},
    "h_DoubleDisCo_disc2_?"   : {"logY" : False, "Y" : {"title" : "A.U.",   "min" : 2e-3}, "X" : {"title" : "Disc 2",                "rebin" : 20, "min" :  0,    "max" :    1}},
    "h_DoubleDisCo_massReg_?" : {"logY" : False, "Y" : {"title" : "A.U.",   "min" : 2e-3}, "X" : {"title" : "Regression Mass [GeV]", "rebin" : 5,  "min" :  0,    "max" : 1000}},
}

backgrounds = {
    "QCD" : {"name" : "QCD multijet", "color" : 30,  "lstyle" : 1, "mstyle" : 8, "lsize" : 0, "msize" : 0},
}

signals = OrderedDict({
    #"RPV_2t6j_mStop-350" : {"name" : "RPV m_{ #tilde{t}} = 350 GeV", "color" : 2, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    #"RPV_2t6j_mStop-550" : {"name" : "RPV m_{ #tilde{t}} = 550 GeV", "color" : 7, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    #"RPV_2t6j_mStop-850" : {"name" : "RPV m_{ #tilde{t}} = 850 GeV", "color" : 4, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
})

data = {
    #"Data" : {"name" : "Data", "color" : ROOT.kBlack, "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 1.3}
    #"pseudoDataS" : {"name" : "pseudoDataS", "color" : ROOT.kBlack, "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 1.3}
}
