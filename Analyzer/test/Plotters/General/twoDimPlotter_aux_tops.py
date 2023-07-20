#! /bin/env/python

import ROOT

from collections import OrderedDict

selections = [
    "",
]

histograms = {
    "h_nRtops_vs_nMtops_0l_ABCD" : {"logX" : False, "logY" : False, "logZ" : False, "Y" : {"title" : "Number of Merged Top Tags", "min" : -0.5, "max" : 5.5}, "X" : {"title" : "Number of Resolved Top Tags", "min" : -0.5, "max" : 5.5}},
}

backgrounds = {
    "TT"       : {"name" : "t#bar{t} + jets", "option" : "COLZ E TEXT"},
    "QCD"      : {"name" : "QCD multijet",    "option" : "COLZ E TEXT"},
    #"TTX"      : {"name" : "t#bar{t} + X",    "option" : "COLZ"},
    #"BG_OTHER" : {"name" : "Other",           "option" : "COLZ"}
}

signals = OrderedDict({
    "StealthSYY_2t6j_mStop-350"  : {"name" : "Stealth SYY m_{ #tilde{t}} = 350 GeV",  "option" : "TEXT E COLZ"},
    "StealthSYY_2t6j_mStop-550"  : {"name" : "Stealth SYY m_{ #tilde{t}} = 550 GeV",  "option" : "TEXT E COLZ"},
    "StealthSYY_2t6j_mStop-850"  : {"name" : "Stealth SYY m_{ #tilde{t}} = 850 GeV",  "option" : "TEXT E COLZ"},
    "StealthSYY_2t6j_mStop-1150" : {"name" : "Stealth SYY m_{ #tilde{t}} = 1150 GeV", "option" : "TEXT E COLZ"}
})

data = {
}
