#! /bin/env/python

import ROOT

from collections import OrderedDict

selections = ["0l", "1l", "2l"]

njets = [str(njet) for njet in range(6,13)] + ["11incl", "12incl", "13incl"]

histograms = {

    # Regression Mass
    "h_DoubleDisCo_%s_massReg?_Njets@_ABCD"%(model) : {"logY" : False,  "orders" : njets, "Y" : {"title" : "A.U.", "min" : 2e-3}, "X" : {"title" : "Regression Mass [GeV]",  "rebin" : 2, "min" : 0,  "max" : 1500}},
    # NN Discs.
    "h_DoubleDisCo_%s_disc1?_Njets@_ABCD"%(model)   : {"logY" : False,  "orders" : njets, "Y" : {"title" : "A.U.", "min" : 2e-3}, "X" : {"title" : "Neural Network Disc. 1", "rebin" : 2, "min" : 0,  "max" : 1   }},
    "h_DoubleDisCo_%s_disc2?_Njets@_ABCD"%(model)   : {"logY" : False,  "orders" : njets, "Y" : {"title" : "A.U.", "min" : 2e-3}, "X" : {"title" : "Neural Network Disc. 2", "rebin" : 2, "min" : 0,  "max" : 1   }},
}

backgrounds = {
    "TT" : {"name" : "${YEAR} t#bar{t} + jets", "loption" : "F", "color" : 40,   "lstyle" : 1, "mstyle" : 8, "lsize" : 6, "msize" : 0},
}

signals = OrderedDict({
})

data = {
    "TT" : {"name" : "Run2UL t#bar{t} + jets", "color" : ROOT.kBlack, "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 1.3}
}
