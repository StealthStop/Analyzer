#! /bin/env/python

import ROOT

from collections import OrderedDict


selections = ["pass0Lbaseline"] 

histograms = {

    # jets
    "h_NMrgTops_withTopTagSF_?" : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.0002}, "X" : {"title" : "N_{ tops} (merged)",   "rebin" : 1, "min" : -0.5, "max" : 5.5}},
    "h_NMrgTops_noTopTagSF_?"   : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.0002}, "X" : {"title" : "N_{ tops} (merged)",   "rebin" : 1, "min" : -0.5, "max" : 5.5}},
    "h_NResTops_withTopTagSF_?" : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.0002}, "X" : {"title" : "N_{ tops} (resolved)", "rebin" : 1, "min" : -0.5, "max" : 5.5}},
    "h_NResTops_noTopTagSF_?"   : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.0002}, "X" : {"title" : "N_{ tops} (resolved)", "rebin" : 1, "min" : -0.5, "max" : 5.5}},
    "h_NTops_withTopTagSF_?"    : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.0002}, "X" : {"title" : "N_{ tops}",            "rebin" : 1, "min" : -0.5, "max" : 5.5}},
    "h_NTops_noTopTagSF_?"      : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.0002}, "X" : {"title" : "N_{ tops}",            "rebin" : 1, "min" : -0.5, "max" : 5.5}},
}

data = {
    "Data" : {"name" : "Data", "color" : ROOT.kBlack, "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 1.3}
} 

backgrounds = {
    "TT"       : {"name" : "t#bar{t} + jets", "color" : 40,  "lstyle" : 1, "mstyle" : 8, "lsize" : 0, "msize" : 0},
    "QCD"      : {"name" : "QCD multijet",    "color" : 30,  "lstyle" : 1, "mstyle" : 8, "lsize" : 0, "msize" : 0},
    "TTX"      : {"name" : "t#bar{t} + X",    "color" : 38,  "lstyle" : 1, "mstyle" : 8, "lsize" : 0, "msize" : 0},
    "BG_OTHER" : {"name" : "Other",           "color" : 41,  "lstyle" : 1, "mstyle" : 8, "lsize" : 0, "msize" : 0},
}

signals = OrderedDict([])
