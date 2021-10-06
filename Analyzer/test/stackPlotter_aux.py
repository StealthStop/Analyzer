#! /bin/env/python

import ROOT

selections = ["_0l",
              "_1l",
              "_1l_QCDCR"
]

histograms = {
    "Jet_cm_pt_@?"                  : {"logY" : True,  "orders" : list(xrange(1,8)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "Jet @ p_{T} [GeV]",     "rebin" : 2, "min" :  0, "max" : 1500}},
    "Jet_cm_eta_@?"                 : {"logY" : False, "orders" : list(xrange(1,8)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "Jet @ #eta",            "rebin" : 1, "min" : -4, "max" :    4}},
    "Jet_cm_phi_@?"                 : {"logY" : False, "orders" : list(xrange(1,8)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "Jet @ #phi",            "rebin" : 1, "min" : -6, "max" :    6}},
    "Jet_cm_m_@?"                   : {"logY" : False, "orders" : list(xrange(1,8)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "Jet @ mass [GeV]",      "rebin" : 2, "min" :  0, "max" :  150}},
    "Jet_cm_E_@?"                   : {"logY" : False, "orders" : list(xrange(1,8)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "Jet @ energy [GeV]",    "rebin" : 2, "min" :  0, "max" : 1500}},
    "Jet_cm_flavb_@?"               : {"logY" : True,  "orders" : list(xrange(1,8)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "Jet @ DeepFlavor b",    "rebin" : 1, "min" :  0, "max" :    1}},
    "Jet_cm_flavc_@?"               : {"logY" : True,  "orders" : list(xrange(1,8)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "Jet @ DeepFlavor c",    "rebin" : 1, "min" :  0, "max" :    1}},
    "Jet_cm_flavuds_@?"             : {"logY" : True,  "orders" : list(xrange(1,8)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "Jet @ DeepFlavor uds",  "rebin" : 1, "min" :  0, "max" :    1}},
    "Jet_cm_flavq_@?"               : {"logY" : True,  "orders" : list(xrange(1,8)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "Jet @ DeepFlavor q",    "rebin" : 1, "min" :  0, "max" :    1}},
    "Jet_cm_flavg_@?"               : {"logY" : True,  "orders" : list(xrange(1,8)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "Jet @ DeepFlavor g",    "rebin" : 1, "min" :  0, "max" :    1}},
    "h_ht?"                         : {"logY" : True,  "orders" : list(xrange(1,8)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "H_{T} [GeV]",           "rebin" : 5, "min" :  0, "max" : 3500}},
    "fwm@_top6?"                    : {"logY" : False, "orders" : list(xrange(2,6)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "Fox-Wolfram Moment @",  "rebin" : 1, "min" :  0, "max" :    1}},
    "jmt_ev@_top6?"                 : {"logY" : False, "orders" : list(xrange(0,3)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "Jet p-E Eigenvalue @",  "rebin" : 1, "min" :  0, "max" :    1}},
    "fwm@_top6?"                    : {"logY" : False, "orders" : list(xrange(2,6)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "Fox-Wolfram Moment @",  "rebin" : 1, "min" :  0, "max" :    1}},
    "jmt_ev@_top6?"                 : {"logY" : False, "orders" : list(xrange(0,3)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "Jet p-E Eigenvalue @",  "rebin" : 1, "min" :  0, "max" :    1}},
    "Stop@_pt_cm_OldSeed?"          : {"logY" : True,  "orders" : list(xrange(1,3)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "Stop @ p_{T} [GeV]",    "rebin" : 1, "min" :  0, "max" : 1500}},
    "Stop@_mass_cm_OldSeed?"        : {"logY" : False, "orders" : list(xrange(1,3)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "Stop @ mass [GeV]",     "rebin" : 1, "min" :  0, "max" : 1500}},
    "Stop@_eta_cm_OldSeed?"         : {"logY" : False, "orders" : list(xrange(1,3)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "Stop @ #eta",           "rebin" : 1, "min" : -4, "max" :    4}},
    "Stop@_phi_cm_OldSeed?"         : {"logY" : False, "orders" : list(xrange(1,3)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "Stop @ #phi",           "rebin" : 1, "min" : -6, "max" :    6}},
    "h_njets?"                      : {"logY" : True,                                 "Y" : {"title" : "Weighted Events", "min" : 2e-3}, "X" : {"title" : "N_{jets}",              "rebin" : 1, "min" :  6, "max" :   18}},
    "h_DoubleDisCo_massReg?_Njets@" : {"logY" : False, "orders" : list(xrange(7,13)), "Y" : {"title" : "A.U.",            "min" : 2e-3}, "X" : {"title" : "Regression Mass [GeV]", "rebin" : 2, "min" :  0, "max" : 1500}},
    "h_njets?"                      : {"logY" : True,                                 "Y" : {"title" : "A.U.",            "min" : 2e-3}, "X" : {"title" : "N_{jets}", "rebin" : 1, "min" :  -0.5, "max" : 20.5}},

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
    "Data" : {"name" : "Data", "color" : ROOT.kBlack, "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 1.3}
    #"pseudoDataS" : {"name" : "pseudoDataS", "color" : ROOT.kBlack, "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 1.3}
}
