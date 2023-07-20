#! /bin/env/python

import ROOT

from collections import OrderedDict

selections = ["0l", "1l", "2l"]

njets = [str(njet) for njet in range(6,13)] + ["11incl", "12incl", "13incl"]

histograms = {

    "h_DoubleDisCo_RPV_disc1?_Njets@_ABCD"   : {"logY" : False, "orders" : njets, "Y" : {"title" : "A.U.",   "min" : 2e-3}, "X" : {"title" : "Discriminant 1",        "rebin" : 10, "min" :  0, "max" :    1}},
    "h_DoubleDisCo_RPV_disc2?_Njets@_ABCD"   : {"logY" : False, "orders" : njets, "Y" : {"title" : "A.U.",   "min" : 2e-3}, "X" : {"title" : "Discriminant 2",        "rebin" : 10, "min" :  0, "max" :    1}},
    "h_DoubleDisCo_RPV_massReg?"             : {"logY" : False, "orders" : njets, "Y" : {"title" : "A.U.",   "min" : 2e-3}, "X" : {"title" : "Regression Mass [GeV]", "rebin" :  2, "min" :  0, "max" : 1500}},

    "h_DoubleDisCo_SYY_disc1?_Njets@_ABCD"   : {"logY" : False, "orders" : njets, "Y" : {"title" : "A.U.",   "min" : 2e-3}, "X" : {"title" : "Discriminant 1",        "rebin" : 10, "min" :  0, "max" :    1}},
    "h_DoubleDisCo_SYY_disc2?_Njets@_ABCD"   : {"logY" : False, "orders" : njets, "Y" : {"title" : "A.U.",   "min" : 2e-3}, "X" : {"title" : "Discriminant 2",        "rebin" : 10, "min" :  0, "max" :    1}},
    "h_DoubleDisCo_SYY_massReg?_Njets@_ABCD" : {"logY" : False, "orders" : njets, "Y" : {"title" : "A.U.",   "min" : 2e-3}, "X" : {"title" : "Regression Mass [GeV]", "rebin" :  2, "min" :  0, "max" : 1500}},
}

backgrounds = {
    "TT"              : {"name" : "t#bar{t} + jets",                   "loption" : "F", "color" : 40,   "lstyle" : 1, "mstyle" : 8, "lsize" : 6, "msize" : 0},

    #"QCD"             : {"name" : "QCD multijet",                      "loption" : "L", "color" : 30,   "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},

    #"WJets"           : {"name" : "W + jets",                          "loption" : "L", "color" : ROOT.TColor.GetColor("#fb8072"), "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},
    #"DYJetsToLL_M-50" : {"name" : "DY + jets",                         "loption" : "L", "color" : ROOT.TColor.GetColor("#80b1d3"), "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},
    #"Diboson"         : {"name" : "Diboson",                           "loption" : "L", "color" : ROOT.TColor.GetColor("#fdb462"), "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},
    #"Triboson"        : {"name" : "Triboson",                          "loption" : "L", "color" : ROOT.TColor.GetColor("#ffcc33"), "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},

    "TTZ"             : {"name" : "t#bar{t}Z",                         "loption" : "L",  "color" : ROOT.TColor.GetColor("#fb8072"), "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},
    "ttHJets"         : {"name" : "t#bar{t}H + jets",                  "loption" : "L",  "color" : ROOT.TColor.GetColor("#80b1d3"), "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},
    "TTWJets"         : {"name" : "t#bar{t}W + jets",                  "loption" : "L",  "color" : ROOT.TColor.GetColor("#fdb462"), "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},
    "TTTT"            : {"name" : "t#bar{t}t#bar{t}",                  "loption" : "L",  "color" : ROOT.TColor.GetColor("#b3de69"), "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},
    "TTWW"            : {"name" : "t#bar{t} + WW",                     "loption" : "L",  "color" : ROOT.TColor.GetColor("#ffcc33"), "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},
    "TTWZ"            : {"name" : "t#bar{t} + WZ",                     "loption" : "L",  "color" : ROOT.TColor.GetColor("#fccde5"), "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},
    "TTTW"            : {"name" : "t#bar{t} + tW",                     "loption" : "L",  "color" : ROOT.TColor.GetColor("#8dd3c7"), "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},
    "TTZH"            : {"name" : "t#bar{t} + ZH",                     "loption" : "L",  "color" : 28,                              "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},
    "TTWH"            : {"name" : "t#bar{t} + WH",                     "loption" : "L",  "color" : 22,                              "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},
    "TTZZ"            : {"name" : "t#bar{t} + ZZ",                     "loption" : "L",  "color" : 13,                              "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},
    "TTHH"            : {"name" : "t#bar{t} + HH",                     "loption" : "L",  "color" : 9,                               "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},
    "TTTJ"            : {"name" : "t#bar{t} + t + jets",               "loption" : "L",  "color" : 49,                              "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},
}

signals = OrderedDict({
})

data = {
}
