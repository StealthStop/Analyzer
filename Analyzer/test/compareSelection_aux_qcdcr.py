#! /bin/env/python

import ROOT

from collections import OrderedDict

selections = [
    ("_0l_ABCD",30),
    ("_1l_ABCD",40),
    ("_0l_QCDCR_ABCD",2),
    ("_1l_QCDCR_ABCD",7),
]

histograms = {
    "h_Njets?"                                : {"logY" : True,                                 "Y" : {"title" : "A.U.", "min" : 2e-3}, "X" : {"title" : "N_{jets}",              "rebin" : 1, "min" :  4, "max" :   18}},
    "h_njets_12incl?"                         : {"logY" : True,                                 "Y" : {"title" : "A.U.", "min" : 2e-3}, "X" : {"title" : "N_{jets} ABCD",         "rebin" : 1, "min" :  0, "max" :   23}},
    "h_njets_13incl?"                         : {"logY" : True,                                 "Y" : {"title" : "A.U.", "min" : 2e-3}, "X" : {"title" : "N_{jets} ABCD",         "rebin" : 1, "min" :  0, "max" :   23}},
}

samples = OrderedDict({
    #"TT"       : {"name" : "t#bar{t} + jets", "color" : 10,  "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "QCD_skim"      : {"name" : "QCD multijet",    "color" : 0,  "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    # "TTX"      : {"name" : "t#bar{t} + X",    "color" : 38,  "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    # "BG_OTHER" : {"name" : "Other",           "color" : 41,  "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    # "Data" : {"name" : "Data", "color" : ROOT.kBlack, "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 1.3},
})
