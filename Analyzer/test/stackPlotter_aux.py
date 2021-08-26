#! /bin/env/python

import ROOT

histograms = {
    "Jet_cm_pt_@_0l" : {"logY" : True, "orders" : list(xrange(1,13)), "Y" : {"title" : "Weighted Events"}, "X" : {"title" : "Jet @ p_{T} [GeV]", "rebin" : 5, "min" : 0,  "max" : 1500}},
}

backgrounds = {
    "TT"       : {"name" : "t#bar{t} + jets", "color" : 38,  "lstyle" : 1, "mstyle" : 8, "lsize" : 0, "msize" : 0},
    "QCD"      : {"name" : "QCD multijet",    "color" : 30,  "lstyle" : 1, "mstyle" : 8, "lsize" : 0, "msize" : 0},
    "TTX"      : {"name" : "t#bar{t} + X",    "color" : 45,  "lstyle" : 1, "mstyle" : 8, "lsize" : 0, "msize" : 0},
    "BG_OTHER" : {"name" : "Other",           "color" : 40,  "lstyle" : 1, "mstyle" : 8, "lsize" : 0, "msize" : 0}
}

signals = {
    "RPV_2t6j_mStop-350" : {"name" : "RPV m_{ #tilde{t}} = 350 GeV", "color" : 2, "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "RPV_2t6j_mStop-550" : {"name" : "RPV m_{ #tilde{t}} = 550 GeV", "color" : 7, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "RPV_2t6j_mStop-850" : {"name" : "RPV m_{ #tilde{t}} = 850 GeV", "color" : 4, "lstyle" : 3, "mstyle" : 8, "lsize" : 3, "msize" : 0}
}

data = {
    "Data" : {"name" : "Data", "color" : ROOT.kBlack, "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 1.3}
}
