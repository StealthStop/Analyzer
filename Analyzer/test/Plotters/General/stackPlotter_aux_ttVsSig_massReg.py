#! /bin/env/python

import ROOT

selections = ["_0l",
              "_1l",
]

histograms = {
    "h_DoubleDisCo_massReg?_Njets@_ABCD"  : {"logY" : False,  "orders" : list(xrange(7,12)), "Y" : {"title" : "A.U.", "min" : 2e-3}, "X" : {"title" : "Regression Mass [GeV]",   "rebin" : 2, "min" : 0,  "max" : 1500}},
    "h_DoubleDisCo_massReg?_Njets11incl_ABCD"  : {"logY" : False, "Y" : {"title" : "A.U.", "min" : 2e-3}, "X" : {"title" : "Regression Mass [GeV]",   "rebin" : 2, "min" : 0,  "max" : 1500}},
    "h_DoubleDisCo_massReg?_Njets12incl_ABCD"  : {"logY" : False, "Y" : {"title" : "A.U.", "min" : 2e-3}, "X" : {"title" : "Regression Mass [GeV]",   "rebin" : 2, "min" : 0,  "max" : 1500}},
    "h_DoubleDisCo_disc1?_Njets@"  : {"logY" : False,  "orders" : list(xrange(7,12)), "Y" : {"title" : "A.U.", "min" : 2e-3}, "X" : {"title" : "Neural Network Disc. 1",   "rebin" : 2, "min" : 0,  "max" : 1}},
    "h_DoubleDisCo_disc2?_Njets@"  : {"logY" : False,  "orders" : list(xrange(7,12)), "Y" : {"title" : "A.U.", "min" : 2e-3}, "X" : {"title" : "Neural Network Disc. 2",   "rebin" : 2, "min" : 0,  "max" : 1}},
    "h_DoubleDisCo_disc1?_Njets11incl"  : {"logY" : False,  "Y" : {"title" : "A.U.", "min" : 2e-3}, "X" : {"title" : "Neural Network Disc. 1",   "rebin" : 2, "min" : 0,  "max" : 1}},
    "h_DoubleDisCo_disc2?_Njets11incl"  : {"logY" : False,  "Y" : {"title" : "A.U.", "min" : 2e-3}, "X" : {"title" : "Neural Network Disc. 2",   "rebin" : 2, "min" : 0,  "max" : 1}},
    "h_DoubleDisCo_disc1?_Njets12incl"  : {"logY" : False,  "Y" : {"title" : "A.U.", "min" : 2e-3}, "X" : {"title" : "Neural Network Disc. 1",   "rebin" : 2, "min" : 0,  "max" : 1}},
    "h_DoubleDisCo_disc2?_Njets12incl"  : {"logY" : False,  "Y" : {"title" : "A.U.", "min" : 2e-3}, "X" : {"title" : "Neural Network Disc. 2",   "rebin" : 2, "min" : 0,  "max" : 1}},

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
