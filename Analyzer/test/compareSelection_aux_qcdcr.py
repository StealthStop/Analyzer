#! /bin/env/python

import ROOT

from collections import OrderedDict

selections = [
    #("_ABCD",       "QCD_{MC}^{SR}", 1, 8, ROOT.TColor.GetColor("#524585")),
    #("_QCDCR_ABCD", "QCD_{MC}^{CR}", 1, 8, ROOT.TColor.GetColor("#c9a8c4")),
    ("_ABCD",       "QCD_{MC}^{SR}",  2, 8, ROOT.TColor.GetColor("#88258c")),
    ("_QCDCR_ABCD", "QCD_{MC}^{CR}",  2, 21,ROOT.TColor.GetColor("#5cb4e8")),
]

histograms = {
    "h_Njets_1l?"                 : {"logY" : True,  "Y" : {"title" : "A.U.", "min" : 2e-3}, "X" : {"title" : "N_{ jets}",         "rebin" : 1, "min" :  5, "max" :   17}},
    "h_DoubleDisCo_SYY_disc1_1l?" : {"logY" : False, "Y" : {"title" : "A.U.", "min" : 2e-3}, "X" : {"title" : "NN Discriminant 1", "rebin" : 4, "min" :  0, "max" :    1}},
    "h_DoubleDisCo_SYY_disc2_1l?" : {"logY" : False, "Y" : {"title" : "A.U.", "min" : 2e-3}, "X" : {"title" : "NN Discriminant 2", "rebin" : 4, "min" :  0, "max" :    1}},
    "h_Njets_0l?"                 : {"logY" : True,  "Y" : {"title" : "A.U.", "min" : 2e-3}, "X" : {"title" : "N_{ jets}",         "rebin" : 1, "min" :  5, "max" :   17}},
    "h_DoubleDisCo_SYY_disc1_0l?" : {"logY" : False, "Y" : {"title" : "A.U.", "min" : 2e-3}, "X" : {"title" : "NN Discriminant 1", "rebin" : 4, "min" :  0, "max" :    1}},
    "h_DoubleDisCo_SYY_disc2_0l?" : {"logY" : False, "Y" : {"title" : "A.U.", "min" : 2e-3}, "X" : {"title" : "NN Discriminant 2", "rebin" : 4, "min" :  0, "max" :    1}},

}

samples = OrderedDict({
    #"TT"       : {"name" : "t#bar{t} + jets", "color" : 10,  "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "QCD"      : {"name" : "QCD multijet",    "color" : 0,  "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 3},
    # "TTX"      : {"name" : "t#bar{t} + X",    "color" : 38,  "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    # "BG_OTHER" : {"name" : "Other",           "color" : 41,  "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    # "Data" : {"name" : "Data", "color" : ROOT.kBlack, "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 1.3},
})
