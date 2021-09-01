#! /bin/env/python

import ROOT

histograms = {
    "Jet_cm_pt_@_?"           : {"logY" : True,  "channels" : ["0l", "1l"], "orders" : list(xrange(1,8)), "Y" : {"title" : "Weighted Events"}, "X" : {"title" : "Jet @ p_{T} [GeV]",    "rebin" : 2, "min" : 0,  "max" : 1500}},
    "Jet_cm_eta_@_?"          : {"logY" : False, "channels" : ["0l", "1l"], "orders" : list(xrange(1,8)), "Y" : {"title" : "Weighted Events"}, "X" : {"title" : "Jet @ #eta",           "rebin" : 1, "min" : -4, "max" : 4}},
    "Jet_cm_phi_@_?"          : {"logY" : False, "channels" : ["0l", "1l"], "orders" : list(xrange(1,8)), "Y" : {"title" : "Weighted Events"}, "X" : {"title" : "Jet @ #phi",           "rebin" : 1, "min" : -6, "max" : 6}},
    "Jet_cm_m_@_?"            : {"logY" : True,  "channels" : ["0l", "1l"], "orders" : list(xrange(1,8)), "Y" : {"title" : "Weighted Events"}, "X" : {"title" : "Jet @ mass [GeV]",     "rebin" : 1, "min" : 0,  "max" : 1500}},
    #"Jet_cm_flavb_@_?"        : {"logY" : False, "channels" : ["0l", "1l"], "orders" : list(xrange(1,8)), "Y" : {"title" : "Weighted Events"}, "X" : {"title" : "Jet @ DeepFlavor b",   "rebin" : 1, "min" : 0,  "max" : 1}},
    #"Jet_cm_flavc_@_?"        : {"logY" : False, "channels" : ["0l", "1l"], "orders" : list(xrange(1,8)), "Y" : {"title" : "Weighted Events"}, "X" : {"title" : "Jet @ DeepFlavor c",   "rebin" : 1, "min" : 0,  "max" : 1}},
    #"Jet_cm_flavuds_@_?"      : {"logY" : False, "channels" : ["0l", "1l"], "orders" : list(xrange(1,8)), "Y" : {"title" : "Weighted Events"}, "X" : {"title" : "Jet @ DeepFlavor uds", "rebin" : 1, "min" : 0,  "max" : 1}},
    #"Jet_cm_flavq_@_?"        : {"logY" : True, "channels" : ["0l", "1l"], "orders" : list(xrange(1,8)), "Y" : {"title" : "Weighted Events"}, "X" : {"title" : "Jet @ DeepFlavor q", "rebin" : 1, "min" : 0,  "max" : 1}},
    #"Jet_cm_flavg_@_?"        : {"logY" : False, "channels" : ["0l", "1l"], "orders" : list(xrange(1,8)), "Y" : {"title" : "Weighted Events"}, "X" : {"title" : "Jet @ DeepFlavor g",   "rebin" : 1, "min" : 0,  "max" : 1}},
    "h_ht_1l"                 : {"logY" : True,  "channels" : ["0l", "1l"], "orders" : list(xrange(1,8)), "Y" : {"title" : "Weighted Events"}, "X" : {"title" : "H_{T} [GeV]",          "rebin" : 5, "min" : 0,  "max" : 3000}},
    "h_ht_0l"                 : {"logY" : True,  "channels" : ["0l", "1l"], "orders" : list(xrange(1,8)), "Y" : {"title" : "Weighted Events"}, "X" : {"title" : "H_{T} [GeV]",          "rebin" : 5, "min" : 0,  "max" : 3000}},
    "fwm@_top6_?"             : {"logY" : False, "channels" : ["0l", "1l"], "orders" : list(xrange(2,6)), "Y" : {"title" : "Weighted Events"}, "X" : {"title" : "Fox-Wolfram Moment @", "rebin" : 1, "min" : 0,  "max" : 1}},
    "jmt_ev@_top6_?"          : {"logY" : False, "channels" : ["0l", "1l"], "orders" : list(xrange(0,3)), "Y" : {"title" : "Weighted Events"}, "X" : {"title" : "Jet p-E Eigenvalue @", "rebin" : 1, "min" : 0,  "max" : 1}},
    "fwm@_top6_?"             : {"logY" : False, "channels" : ["0l", "1l"], "orders" : list(xrange(2,6)), "Y" : {"title" : "Weighted Events"}, "X" : {"title" : "Fox-Wolfram Moment @", "rebin" : 1, "min" : 0,  "max" : 1}},
    "jmt_ev@_top6_?"          : {"logY" : False, "channels" : ["0l", "1l"], "orders" : list(xrange(0,3)), "Y" : {"title" : "Weighted Events"}, "X" : {"title" : "Jet p-E Eigenvalue @", "rebin" : 1, "min" : 0,  "max" : 1}},
    "Stop@_mass_cm_OldSeed_?" : {"logY" : True,  "channels" : ["0l", "1l"], "orders" : list(xrange(1,3)), "Y" : {"title" : "Weighted Events"}, "X" : {"title" : "Stop @ p_{T} [GeV]",   "rebin" : 1, "min" : 0,  "max" : 1500}},
    "Stop@_pt_cm_OldSeed_?"   : {"logY" : True,  "channels" : ["0l", "1l"], "orders" : list(xrange(1,3)), "Y" : {"title" : "Weighted Events"}, "X" : {"title" : "Stop @ #eta",          "rebin" : 1, "min" : -6,  "max" : 6}},
    "Stop@_eta_cm_OldSeed_?"  : {"logY" : False, "channels" : ["0l", "1l"], "orders" : list(xrange(1,3)), "Y" : {"title" : "Weighted Events"}, "X" : {"title" : "Stop @ #phi",          "rebin" : 1, "min" : -4,  "max" : 4}},
    "Stop@_phi_cm_OldSeed_?"  : {"logY" : False, "channels" : ["0l", "1l"], "orders" : list(xrange(1,3)), "Y" : {"title" : "Weighted Events"}, "X" : {"title" : "Stop @ mass [GeV]",    "rebin" : 1, "min" : 0,  "max" : 1500}},
}

backgrounds = {
    "TT"       : {"name" : "t#bar{t} + jets", "color" : 40,  "lstyle" : 1, "mstyle" : 8, "lsize" : 0, "msize" : 0},
    "QCD"      : {"name" : "QCD multijet",    "color" : 30,  "lstyle" : 1, "mstyle" : 8, "lsize" : 0, "msize" : 0},
    "TTX"      : {"name" : "t#bar{t} + X",    "color" : 38,  "lstyle" : 1, "mstyle" : 8, "lsize" : 0, "msize" : 0},
    "BG_OTHER" : {"name" : "Other",           "color" : 41,  "lstyle" : 1, "mstyle" : 8, "lsize" : 0, "msize" : 0}
}

signals = {
    "RPV_2t6j_mStop-350"        : {"name" : "RPV m_{ #tilde{t}} = 350 GeV",               "color" : 2, "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "RPV_2t6j_mStop-550"        : {"name" : "RPV m_{ #tilde{t}} = 550 GeV",               "color" : 7, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "StealthSYY_2t6j_mStop-850" : {"name" : "Stealth SY#bar{Y} m_{ #tilde{t}} = 850 GeV", "color" : 4, "lstyle" : 3, "mstyle" : 8, "lsize" : 3, "msize" : 0}
}

data = {
    "pseudoDataS" : {"name" : "pseudoDataS", "color" : ROOT.kBlack, "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 1.3}
}
