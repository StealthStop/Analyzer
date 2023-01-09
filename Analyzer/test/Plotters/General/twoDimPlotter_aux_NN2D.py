#! /bin/env/python

import ROOT

from collections import OrderedDict

selections = [
    "",
]

histograms = {
    "h_DoubleDisCo_SYY_disc1_disc2_0l_Njets@_ABCD" : {"orders" : map(str, list(range(8, 13)))+["13incl"], "logX" : False, "logY" : False, "logZ" : False, "Y" : {"rebin" : 5, "title" : "Neural Network Discriminant 2", "min" : 0.0, "max" : 1.0}, "X" : {"rebin" : 5, "title" : "Neural Network Discriminant 1", "min" : 0.0, "max" : 1.0}},
    "h_DoubleDisCo_SYY_disc1_disc2_1l_Njets@_ABCD" : {"orders" : map(str, list(range(7, 12)))+["12incl"], "logX" : False, "logY" : False, "logZ" : False, "Y" : {"rebin" : 5, "title" : "Neural Network Discriminant 2", "min" : 0.0, "max" : 1.0}, "X" : {"rebin" : 5, "title" : "Neural Network Discriminant 1", "min" : 0.0, "max" : 1.0}},
    "h_DoubleDisCo_SYY_disc1_disc2_2l_Njets@_ABCD" : {"orders" : map(str, list(range(6, 11)))+["11incl"], "logX" : False, "logY" : False, "logZ" : False, "Y" : {"rebin" : 5, "title" : "Neural Network Discriminant 2", "min" : 0.0, "max" : 1.0}, "X" : {"rebin" : 5, "title" : "Neural Network Discriminant 1", "min" : 0.0, "max" : 1.0}},
    "h_DoubleDisCo_RPV_disc1_disc2_0l_Njets@_ABCD" : {"orders" : map(str, list(range(8, 13)))+["13incl"], "logX" : False, "logY" : False, "logZ" : False, "Y" : {"rebin" : 5, "title" : "Neural Network Discriminant 2", "min" : 0.0, "max" : 1.0}, "X" : {"rebin" : 5, "title" : "Neural Network Discriminant 1", "min" : 0.0, "max" : 1.0}},
    "h_DoubleDisCo_RPV_disc1_disc2_1l_Njets@_ABCD" : {"orders" : map(str, list(range(7, 12)))+["12incl"], "logX" : False, "logY" : False, "logZ" : False, "Y" : {"rebin" : 5, "title" : "Neural Network Discriminant 2", "min" : 0.0, "max" : 1.0}, "X" : {"rebin" : 5, "title" : "Neural Network Discriminant 1", "min" : 0.0, "max" : 1.0}},
    "h_DoubleDisCo_RPV_disc1_disc2_2l_Njets@_ABCD" : {"orders" : map(str, list(range(6, 11)))+["11incl"], "logX" : False, "logY" : False, "logZ" : False, "Y" : {"rebin" : 5, "title" : "Neural Network Discriminant 2", "min" : 0.0, "max" : 1.0}, "X" : {"rebin" : 5, "title" : "Neural Network Discriminant 1", "min" : 0.0, "max" : 1.0}},
}

backgrounds = {
    "TT"       : {"name" : "t#bar{t} + jets", "option" : "COLZ"},
    #"QCD"      : {"name" : "QCD multijet",    "option" : "COLZ"},
    #"TTX"      : {"name" : "t#bar{t} + X",    "option" : "COLZ"},
    #"BG_OTHER" : {"name" : "Other",           "option" : "COLZ"}
}

signals = OrderedDict({
    "RPV_2t6j_mStop-350" : {"name" : "RPV m_{ #tilde{t}} = 350 GeV",  "option" : "COLZ"},
    "RPV_2t6j_mStop-550" : {"name" : "RPV m_{ #tilde{t}} = 550 GeV",  "option" : "COLZ"},
    "RPV_2t6j_mStop-850" : {"name" : "RPV m_{ #tilde{t}} = 850 GeV",  "option" : "COLZ"},
})

data = {
}
