#! /bin/env/python

import ROOT

from collections import OrderedDict

controlRegions = {
    "h_njets_12incl_1l_QCDCR_ABCD"  : {"logY" : True, "orders" : list(xrange(1,2)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "N_{Jets}",    "rebin" : 1,  "min" : 0,  "max" :    24}},
}

signalRegions = {
    "h_njets_12incl_0l_ABCD"  : {"logY" : True, "orders" : list(xrange(1,2)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "N_{Jets}",    "rebin" : 1,  "min" : 0,  "max" :    24}},
    "h_njets_12incl_1l_ABCD"  : {"logY" : True, "orders" : list(xrange(1,2)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "N_{Jets}",    "rebin" : 1,  "min" : 0,  "max" :    24}},
    #"h_njets_12incl_2l_ABCD"  : {"logY" : True, "orders" : list(xrange(1,2)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "N_{Jets}",    "rebin" : 1,  "min" : 0,  "max" :    24}},
}

backgrounds = {
    "TT"       : {"name" : "t#bar{t} + jets", "color" : 40,  "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "QCD"      : {"name" : "QCD multijet",    "color" : 30,  "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "TTX"      : {"name" : "t#bar{t} + X",    "color" : 38,  "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "BG_OTHER" : {"name" : "Other",           "color" : 41,  "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
}

signals = OrderedDict({
    #"RPV_350"        : {"name" : "RPV m_{ #tilde{t}} = 350 GeV",               "color" : 2, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    #"SYY_550"        : {"name" : "SYY m_{ #tilde{t}} = 550 GeV",               "color" : 7, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    #"RPV_850"        : {"name" : "RPV m_{ #tilde{t}} = 850 GeV",               "color" : 4, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "RPV_2t6j_mStop-350"        : {"name" : "RPV m_{ #tilde{t}} = 350 GeV",                 "color" : 2, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "RPV_2t6j_mStop-850"       : {"name" : "RPV m_{ #tilde{t}} = 850 GeV",                  "color" : 7, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "StealthSYY_2t6j_mStop-550" : {"name" : "Stealth SY#bar{Y} m_{ #tilde{t}} = 550 GeV",   "color" : 4, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0}
})

data = {
    "Data" : {"name" : "Data", "color" : ROOT.kBlack, "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 1.3}
}
