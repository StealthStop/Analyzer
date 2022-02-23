#! /bin/env/python

import ROOT

selections = ["_1l"]

histograms = {
    "Jet_cm_pt_@?"           : {"logY" : True,  "orders" : list(xrange(1,7)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ p_{T} [GeV]",    "rebin" : 2, "min" :  0,  "max" : 1500}},
    "Jet_cm_eta_@?"          : {"logY" : False, "orders" : list(xrange(1,7)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ #eta",           "rebin" : 1, "min" : -4,  "max" :    4}},
    "Jet_cm_phi_@?"          : {"logY" : False, "orders" : list(xrange(1,7)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ #phi",           "rebin" : 1, "min" : -6,  "max" :    6}},
    "Jet_cm_m_@?"            : {"logY" : False, "orders" : list(xrange(1,7)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ mass [GeV]",     "rebin" : 2, "min" :  0,  "max" :  150}},
    "Jet_cm_flavb_@?"        : {"logY" : True,  "orders" : list(xrange(1,7)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ DeepFlavor b",   "rebin" : 1, "min" :  0,  "max" :    1}},
    "Jet_cm_flavc_@?"        : {"logY" : True,  "orders" : list(xrange(1,7)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ DeepFlavor c",   "rebin" : 1, "min" :  0,  "max" :    1}},
    "Jet_cm_flavuds_@?"      : {"logY" : True,  "orders" : list(xrange(1,7)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ DeepFlavor uds", "rebin" : 1, "min" :  0,  "max" :    1}},
    "Jet_cm_flavq_@?"        : {"logY" : True,  "orders" : list(xrange(1,7)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ DeepFlavor q",   "rebin" : 1, "min" :  0,  "max" :    1}},
    "Jet_cm_flavg_@?"        : {"logY" : True,  "orders" : list(xrange(1,7)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ DeepFlavor g",   "rebin" : 1, "min" :  0,  "max" :    1}},
    "h_ht?"                  : {"logY" : True,                                "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "H_{T} [GeV]",          "rebin" : 5, "min" :  0,  "max" : 3500}},
    "fwm@_top6?"             : {"logY" : False, "orders" : list(xrange(2,6)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Fox-Wolfram Moment @", "rebin" : 1, "min" :  0,  "max" :    1}},
    "jmt_ev@_top6?"          : {"logY" : False, "orders" : list(xrange(0,3)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet p-E Eigenvalue @", "rebin" : 1, "min" :  0,  "max" :    1}},
    "Stop@_pt_cm_OldSeed?"   : {"logY" : True,  "orders" : list(xrange(1,3)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Stop @ p_{T} [GeV]",   "rebin" : 1, "min" :  0,  "max" : 1500}},
    "Stop@_mass_cm_OldSeed?" : {"logY" : False, "orders" : list(xrange(1,3)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Stop @ mass [GeV]",    "rebin" : 1, "min" :  0,  "max" : 1500}},
    "Stop@_eta_cm_OldSeed?"  : {"logY" : False, "orders" : list(xrange(1,3)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Stop @ #eta",          "rebin" : 1, "min" : -4,  "max" :    4}},
    "Stop@_phi_cm_OldSeed?"  : {"logY" : False, "orders" : list(xrange(1,3)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Stop @ #phi",          "rebin" : 1, "min" : -6,  "max" :    6}},
}

backgrounds = {
    "TT"       : {"name" : "t#bar{t} + jets", "color" : ROOT.TColor.GetColor("#9999FF"),  "fstyle" : 1001, "lcolor" : ROOT.TColor.GetColor("#010199"), "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0, "option" : "FHIST", "loption" : "F"},
}

signals = {
    "RPV_2t6j_mStop-350"        : {"name" : "RPV m_{ #tilde{t}} = 350 GeV",  "color" : ROOT.TColor.GetColor("#D325D3"), "fstyle" : 3004,  "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0, "option" : "FHIST", "loption" : "F"},
    "RPV_2t6j_mStop-850"        : {"name" : "RPV m_{ #tilde{t}} = 850 GeV",  "color" : ROOT.TColor.GetColor("#FFA851"), "fstyle" : 3005,  "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0, "option" : "FHIST", "loption" : "F"},
}

data = {
}
