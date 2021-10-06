#! /bin/env/python

import ROOT

selections = ["_0l",
              "_1l",
]

histograms = {
    "h_DoubleDisCo_massReg?_Njets@"  : {"logY" : False,  "orders" : list(xrange(6,13)), "Y" : {"title" : "A.U.", "min" : 2e-3}, "X" : {"title" : "Regression Mass [GeV]",   "rebin" : 2, "min" : 0,  "max" : 1500}},
}

backgrounds = {
    "TT"       : {"name" : "t#bar{t} + jets", "color" : 40,  "lstyle" : 1, "mstyle" : 8, "lsize" : 0, "msize" : 0},
}

signals = {
    "RPV_2t6j_mStop-350"        : {"name" : "RPV m_{ #tilde{t}} = 350 GeV",  "color" : 2,              "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "RPV_2t6j_mStop-550"        : {"name" : "RPV m_{ #tilde{t}} = 550 GeV",  "color" : 4,              "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "RPV_2t6j_mStop-850"        : {"name" : "RPV m_{ #tilde{t}} = 850 GeV",  "color" : 6,              "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "RPV_2t6j_mStop-1250"       : {"name" : "RPV m_{ #tilde{t}} = 1250 GeV", "color" : ROOT.kCyan+1,   "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
}

data = {
}
