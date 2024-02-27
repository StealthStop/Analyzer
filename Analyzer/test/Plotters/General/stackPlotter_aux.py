#! /bin/env/python

import ROOT

from collections import OrderedDict

selections = [
        #"_0l",
        "_1l",
        #"_1l_QCDCR",
        
        # >= 6 jets
        "0l_0NonIsoMuon_HT500_ge6j",
        "0l_0NonIsoMuon_HT500_ge6j_ge1t",
        "0l_0NonIsoMuon_HT500_ge6j_ge1t_ge1b",
        "0l_0NonIsoMuon_HT500_ge6j_ge1t_ge1b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge6j_ge1tR_ge1b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge6j_ge1tM_ge1b_ge1dRbjets",

        "0l_0NonIsoMuon_HT500_ge6j_ge1t_ge2b",
        "0l_0NonIsoMuon_HT500_ge6j_ge1t_ge2b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge6j_ge1tR_ge2b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge6j_ge1tM_ge2b_ge1dRbjets",

        "0l_0NonIsoMuon_HT500_ge6j_ge2t",
        "0l_0NonIsoMuon_HT500_ge6j_ge2t_ge1b",
        "0l_0NonIsoMuon_HT500_ge6j_ge2t_ge1b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge6j_ge2tR_ge1b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge6j_ge2tM_ge1b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge6j_ge2tRM_ge1b_ge1dRbjets",

        "0l_0NonIsoMuon_HT500_ge6j_ge2t_ge2b",
        "0l_0NonIsoMuon_HT500_ge6j_ge2t_ge2b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge6j_ge2tR_ge2b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge6j_ge2tM_ge2b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge6j_ge2tRM_ge2b_ge1dRbjets",

        # >= 7 jets
        "0l_0NonIsoMuon_HT500_ge7j",
        "0l_0NonIsoMuon_HT500_ge7j_ge1t",
        "0l_0NonIsoMuon_HT500_ge7j_ge1t_ge1b",
        "0l_0NonIsoMuon_HT500_ge7j_ge1t_ge1b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge7j_ge1tR_ge1b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge7j_ge1tM_ge1b_ge1dRbjets",

        "0l_0NonIsoMuon_HT500_ge7j_ge1t_ge2b",
        "0l_0NonIsoMuon_HT500_ge7j_ge1t_ge2b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge7j_ge1tR_ge2b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge7j_ge1tM_ge2b_ge1dRbjets",

        "0l_0NonIsoMuon_HT500_ge7j_ge2t",
        "0l_0NonIsoMuon_HT500_ge7j_ge2t_ge1b",
        "0l_0NonIsoMuon_HT500_ge7j_ge2t_ge1b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge7j_ge2tR_ge1b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge7j_ge2tM_ge1b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge7j_ge2tRM_ge1b_ge1dRbjets",

        "0l_0NonIsoMuon_HT500_ge7j_ge2t_ge2b",
        "0l_0NonIsoMuon_HT500_ge7j_ge2t_ge2b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge7j_ge2tR_ge2b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge7j_ge2tM_ge2b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge7j_ge2tRM_ge2b_ge1dRbjets",

        # >= 8 jets
        "0l_0NonIsoMuon_HT500_ge8j",
        "0l_0NonIsoMuon_HT500_ge8j_ge1t",
        "0l_0NonIsoMuon_HT500_ge8j_ge1t_ge1b",
        "0l_0NonIsoMuon_HT500_ge8j_ge1t_ge1b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge8j_ge1tR_ge1b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge8j_ge1tM_ge1b_ge1dRbjets",

        "0l_0NonIsoMuon_HT500_ge8j_ge1t_ge2b",
        "0l_0NonIsoMuon_HT500_ge8j_ge1t_ge2b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge8j_ge1tR_ge2b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge8j_ge1tM_ge2b_ge1dRbjets",

        "0l_0NonIsoMuon_HT500_ge8j_ge2t",
        "0l_0NonIsoMuon_HT500_ge8j_ge2t_ge1b",
        "0l_0NonIsoMuon_HT500_ge8j_ge2t_ge1b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge8j_ge2tR_ge1b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge8j_ge2tM_ge1b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge8j_ge2tRM_ge1b_ge1dRbjets",

        "0l_0NonIsoMuon_HT500_ge8j_ge2t_ge2b",
        "0l_0NonIsoMuon_HT500_ge8j_ge2t_ge2b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge8j_ge2tR_ge2b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge8j_ge2tM_ge2b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge8j_ge2tRM_ge2b_ge1dRbjets",

        # >= 9 jets
        "0l_0NonIsoMuon_HT500_ge9j",
        "0l_0NonIsoMuon_HT500_ge9j_ge1t",
        "0l_0NonIsoMuon_HT500_ge9j_ge1t_ge1b",
        "0l_0NonIsoMuon_HT500_ge9j_ge1t_ge1b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge9j_ge1tR_ge1b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge9j_ge1tM_ge1b_ge1dRbjets",

        "0l_0NonIsoMuon_HT500_ge9j_ge1t_ge2b",
        "0l_0NonIsoMuon_HT500_ge9j_ge1t_ge2b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge9j_ge1tR_ge2b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge9j_ge1tM_ge2b_ge1dRbjets",

        "0l_0NonIsoMuon_HT500_ge9j_ge2t",
        "0l_0NonIsoMuon_HT500_ge9j_ge2t_ge1b",
        "0l_0NonIsoMuon_HT500_ge9j_ge2t_ge1b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge9j_ge2tR_ge1b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge9j_ge2tM_ge1b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge9j_ge2tRM_ge1b_ge1dRbjets",

        "0l_0NonIsoMuon_HT500_ge9j_ge2t_ge2b",
        "0l_0NonIsoMuon_HT500_ge9j_ge2t_ge2b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge9j_ge2tR_ge2b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge9j_ge2tM_ge2b_ge1dRbjets",
        "0l_0NonIsoMuon_HT500_ge9j_ge2tRM_ge2b_ge1dRbjets",
]

histograms = {
    #"Jet_cm_pt_@?"                  : {"logY" : True,  "orders" : list(xrange(1,8)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "Jet @ p_{T} [GeV]",     "rebin" : 2, "min" :  0, "max" : 1500}},
    #"Jet_cm_eta_@?"                 : {"logY" : False, "orders" : list(xrange(1,8)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "Jet @ #eta",            "rebin" : 1, "min" : -4, "max" :    4}},
    #"Jet_cm_phi_@?"                 : {"logY" : False, "orders" : list(xrange(1,8)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "Jet @ #phi",            "rebin" : 1, "min" : -6, "max" :    6}},
    #"Jet_cm_m_@?"                   : {"logY" : False, "orders" : list(xrange(1,8)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "Jet @ mass [GeV]",      "rebin" : 2, "min" :  0, "max" :  150}},
    #"Jet_cm_E_@?"                   : {"logY" : False, "orders" : list(xrange(1,8)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "Jet @ energy [GeV]",    "rebin" : 2, "min" :  0, "max" : 1500}},
    #"Jet_cm_flavb_@?"               : {"logY" : True,  "orders" : list(xrange(1,8)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "Jet @ DeepFlavor b",    "rebin" : 1, "min" :  0, "max" :    1}},
    #"Jet_cm_flavc_@?"               : {"logY" : True,  "orders" : list(xrange(1,8)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "Jet @ DeepFlavor c",    "rebin" : 1, "min" :  0, "max" :    1}},
    #"Jet_cm_flavuds_@?"             : {"logY" : True,  "orders" : list(xrange(1,8)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "Jet @ DeepFlavor uds",  "rebin" : 1, "min" :  0, "max" :    1}},
    #"Jet_cm_flavq_@?"               : {"logY" : True,  "orders" : list(xrange(1,8)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "Jet @ DeepFlavor q",    "rebin" : 1, "min" :  0, "max" :    1}},
    #"Jet_cm_flavg_@?"               : {"logY" : True,  "orders" : list(xrange(1,8)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "Jet @ DeepFlavor g",    "rebin" : 1, "min" :  0, "max" :    1}},
    #"h_ht?"                         : {"logY" : True,  "orders" : list(xrange(1,8)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "H_{T} [GeV]",           "rebin" : 5, "min" :  0, "max" : 3500}},
    #"fwm@_top6?"                    : {"logY" : False, "orders" : list(xrange(2,6)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "Fox-Wolfram Moment @",  "rebin" : 1, "min" :  0, "max" :    1}},
    #"jmt_ev@_top6?"                 : {"logY" : False, "orders" : list(xrange(0,3)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "Jet p-E Eigenvalue @",  "rebin" : 1, "min" :  0, "max" :    1}},
    #"fwm@_top6?"                    : {"logY" : False, "orders" : list(xrange(2,6)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "Fox-Wolfram Moment @",  "rebin" : 1, "min" :  0, "max" :    1}},
    #"jmt_ev@_top6?"                 : {"logY" : False, "orders" : list(xrange(0,3)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "Jet p-E Eigenvalue @",  "rebin" : 1, "min" :  0, "max" :    1}},
    #"Stop@_pt_cm_OldSeed?"          : {"logY" : True,  "orders" : list(xrange(1,3)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "Stop @ p_{T} [GeV]",    "rebin" : 1, "min" :  0, "max" : 1500}},
    #"Stop@_mass_cm_OldSeed?"        : {"logY" : False, "orders" : list(xrange(1,3)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "Stop @ mass [GeV]",     "rebin" : 1, "min" :  0, "max" : 1500}},
    #"Stop@_eta_cm_OldSeed?"         : {"logY" : False, "orders" : list(xrange(1,3)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "Stop @ #eta",           "rebin" : 1, "min" : -4, "max" :    4}},
    #"Stop@_phi_cm_OldSeed?"         : {"logY" : False, "orders" : list(xrange(1,3)),  "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "Stop @ #phi",           "rebin" : 1, "min" : -6, "max" :    6}},
    #"h_njets?"                      : {"logY" : True,                                 "Y" : {"title" : "Weighted Events", "min" : 2e-3}, "X" : {"title" : "N_{jets}",              "rebin" : 1, "min" :  6, "max" :   18}},
    #"h_DoubleDisCo_massReg?_Njets@" : {"logY" : False, "orders" : list(xrange(7,13)), "Y" : {"title" : "A.U.",            "min" : 2e-3}, "X" : {"title" : "Regression Mass [GeV]", "rebin" : 2, "min" :  0, "max" : 1500}},
    "h_ntops?"                      : {"logY" : True,                                 "Y" : {"title" : "Events",         "min" : 2e-3}, "X" : {"title" : "N_{tops}", "rebin" : 1, "min" :  -0.5, "max" : 10.5}},
    "h_njets?"                      : {"logY" : True,                                 "Y" : {"title" : "Events",         "min" : 2e-3}, "X" : {"title" : "N_{jets}", "rebin" : 1, "min" :  -0.5, "max" : 20.5}},
}

backgrounds = {
    "TT"       : {"name" : "t#bar{t} + jets", "color" : 40,  "lstyle" : 1, "mstyle" : 8, "lsize" : 0, "msize" : 0},
    "QCD"      : {"name" : "QCD multijet",    "color" : 30,  "lstyle" : 1, "mstyle" : 8, "lsize" : 0, "msize" : 0},
    "TTX"      : {"name" : "t#bar{t} + X",    "color" : 38,  "lstyle" : 1, "mstyle" : 8, "lsize" : 0, "msize" : 0},
    "BG_OTHER" : {"name" : "Other",           "color" : 41,  "lstyle" : 1, "mstyle" : 8, "lsize" : 0, "msize" : 0}
}

signals = OrderedDict({
    "RPV_2t6j_mStop-350"        : {"name" : "RPV m_{ #tilde{t}} = 350 GeV",               "color" : 2, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "RPV_2t6j_mStop-550"        : {"name" : "RPV m_{ #tilde{t}} = 550 GeV",               "color" : 7, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "RPV_2t6j_mStop-850"        : {"name" : "RPV m_{ #tilde{t}} = 850 GeV",               "color" : 4, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    #"RPV_2t6j_mStop-1050"       : {"name" : "RPV m_{ #tilde{t}} = 1050 GeV",              "color" : 5, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    #"StealthSYY_2t6j_mStop-900" : {"name" : "Stealth SY#bar{Y} m_{ #tilde{t}} = 900 GeV", "color" : 6, "lstyle" : 3, "mstyle" : 8, "lsize" : 3, "msize" : 0}
})

data = {
    #"Data" : {"name" : "Data", "color" : ROOT.kBlack, "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 1.3}
    #"pseudoDataS" : {"name" : "pseudoDataS", "color" : ROOT.kBlack, "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 1.3}
}
