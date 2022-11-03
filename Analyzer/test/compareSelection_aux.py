#! /bin/env/python

import ROOT

from collections import OrderedDict

selections = [
    ("_0l_ABCD",30),
    ("_1l_ABCD",40),
    ("_2l_ABCD",38),
    ("_QCDCR_ABCD",41),
    ("_QCDCR_1b_ABCD",2),
    ("_QCDCR_1t_ABCD",7),
    ("_QCDCR_1b_1t_ABCD",4),
    ("_QCDCR_2b_ABCD",ROOT.kBlack),
]

histograms = {
    "Jet_cm_ptrHT_@?"        : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ p_{T}/H_{T} [GeV]",    "rebin" : 5, "min" :  0,  "max" : 1500}},
    "Jet_cm_pt_@?"           : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ p_{T} [GeV]",          "rebin" : 2, "min" :  0,  "max" : 1500}},
    "Jet_cm_eta_@?"          : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ #eta",                 "rebin" : 1, "min" : -4,  "max" :    4}},
    "Jet_cm_phi_@?"          : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ #phi",                 "rebin" : 1, "min" : -6,  "max" :    6}},
    "Jet_cm_m_@?"            : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ mass [GeV]",           "rebin" : 2, "min" :  0,  "max" :  150}},
    "Jet_cm_E_@?"            : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ Energy [GeV]",         "rebin" : 5, "min" :  0,  "max" : 1500}},
    "Jet_cm_flavb_@?"        : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ DeepFlavor b",         "rebin" : 5, "min" :  0,  "max" :    1}},
    #"Jet_cm_flavc_@?"        : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ DeepFlavor c",   "rebin" : 1, "min" :  0,  "max" :    1}},
    #"Jet_cm_flavuds_@?"      : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ DeepFlavor uds", "rebin" : 1, "min" :  0,  "max" :    1}},
    #"Jet_cm_flavq_@?"        : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ DeepFlavor q",   "rebin" : 1, "min" :  0,  "max" :    1}},
    #"Jet_cm_flavg_@?"        : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ DeepFlavor g",   "rebin" : 1, "min" :  0,  "max" :    1}},
    #"Jet_cm_CSVb_@?"        : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ DeepCSV b",   "rebin" : 1, "min" :  0,  "max" :    1}},
    #"Jet_cm_CSVc_@?"        : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ DeepCSV c",   "rebin" : 1, "min" :  0,  "max" :    1}},
    #"Jet_cm_CSVudsg_@?"      : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ DeepCSV udsg", "rebin" : 1, "min" :  0,  "max" :    1}},
    #"HT_trigger_pt30?"                  : {"logY" : True,                                "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "H_{T} [GeV]",          "rebin" : 5, "min" :  0,  "max" : 3500}},
    "fwm@_top6?"             : {"logY" : True, "orders" : list(xrange(2,6)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Fox-Wolfram Moment @", "rebin" : 1, "min" :  0,  "max" :    1}},
    "jmt_ev@_top6?"          : {"logY" : True, "orders" : list(xrange(0,3)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet p-E Eigenvalue @", "rebin" : 1, "min" :  0,  "max" :    1}},
    "Stop@_pt_cm_OldSeed?"   : {"logY" : True,  "orders" : list(xrange(1,3)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Stop @ p_{T} [GeV]",   "rebin" : 5, "min" :  0,  "max" : 1500}},
    #"Stop@_ptrHT_cm_OldSeed?"   : {"logY" : True,  "orders" : list(xrange(1,3)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Stop @ p_{T} [GeV]",   "rebin" : 1, "min" :  0,  "max" : 1}},
    "Stop@_mass_cm_OldSeed?" : {"logY" : True, "orders" : list(xrange(1,3)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Stop @ mass [GeV]",    "rebin" : 5, "min" :  0,  "max" : 1500}},
    "Stop@_eta_cm_OldSeed?"  : {"logY" : True, "orders" : list(xrange(1,3)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Stop @ #eta",          "rebin" : 2, "min" : -4,  "max" :    4}},
    "Stop@_phi_cm_OldSeed?"  : {"logY" : True, "orders" : list(xrange(1,3)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Stop @ #phi",          "rebin" : 2, "min" : -6,  "max" :    6}},
    #"Stop@_pt_cm_TopSeed?"   : {"logY" : True,  "orders" : list(xrange(1,3)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Stop @ p_{T} [GeV]",   "rebin" : 1, "min" :  0,  "max" : 1500}},
    #"Stop@_ptrHT_cm_TopSeed?"   : {"logY" : True,  "orders" : list(xrange(1,3)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Stop @ p_{T} [GeV]",   "rebin" : 1, "min" :  0,  "max" : 1}},
    #"Stop@_mass_cm_TopSeed?" : {"logY" : False, "orders" : list(xrange(1,3)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Stop @ mass [GeV]",    "rebin" : 1, "min" :  0,  "max" : 1500}},
    #"Stop@_eta_cm_TopSeed?"  : {"logY" : False, "orders" : list(xrange(1,3)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Stop @ #eta",          "rebin" : 1, "min" : -4,  "max" :    4}},
    #"Stop@_phi_cm_TopSeed?"  : {"logY" : False, "orders" : list(xrange(1,3)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Stop @ #phi",          "rebin" : 1, "min" : -6,  "max" :    6}},
    #"h_nb?"  : {"logY" : False, "orders" : list(xrange(1,2)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "N B_{Jets}",          "rebin" : 1, "min" : 0,  "max" :    6}},
    #"h_dRbjet?"  : {"logY" : False, "orders" : list(xrange(1,2)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "dR B_{Jets}",          "rebin" : 1, "min" : 0,  "max" :    6}},
    #"h_DoubleDisCo_disc1?"        : {"logY" : True, "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Disc1",   "rebin" : 5, "min" :  0,  "max" :    1}},
    #"h_DoubleDisCo_disc2?"        : {"logY" : True, "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Disc2",   "rebin" : 5, "min" :  0,  "max" :    1}},
    "h_ntops?"  : {"logY" : True, "orders" : list(xrange(1,2)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "N Top Tags",          "rebin" : 1, "min" : 0,  "max" :    6}},
    #"h_mMiniIso?"  : {"logY" : False, "orders" : list(xrange(1,2)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "N Top Tags",          "rebin" : 1, "min" : 0.01,  "max" :    1.0}},
    "h_ht?"                                   : {"logY" : True,                                 "Y" : {"title" : "A.U.", "min" : 0.2},  "X" : {"title" : "H_{T} [GeV]",           "rebin" : 10, "min" :  0, "max" : 4000}},
    "h_njets?"                                : {"logY" : True,                                 "Y" : {"title" : "A.U.", "min" : 2e-3}, "X" : {"title" : "N_{jets}",              "rebin" : 1, "min" :  6, "max" :   18}},
}

samples = OrderedDict({
    #"RPV_350"        : {"name" : "RPV m_{ #tilde{t}} = 350 GeV",               "color" : 2, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    #"SYY_550"        : {"name" : "SYY m_{ #tilde{t}} = 550 GeV",               "color" : 7, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    #"RPV_850"        : {"name" : "RPV m_{ #tilde{t}} = 850 GeV",               "color" : 4, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    # "RPV_2t6j_mStop-350"        : {"name" : "RPV m_{ #tilde{t}} = 350 GeV",                 "color" : 2, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    # "RPV_2t6j_mStop-850"       : {"name" : "RPV m_{ #tilde{t}} = 850 GeV",                  "color" : 7, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    # "StealthSYY_2t6j_mStop-550" : {"name" : "Stealth SY#bar{Y} m_{ #tilde{t}} = 550 GeV",   "color" : 4, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    # "TT"       : {"name" : "t#bar{t} + jets", "color" : 40,  "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "QCD"      : {"name" : "QCD multijet",    "color" : 30,  "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    # "TTX"      : {"name" : "t#bar{t} + X",    "color" : 38,  "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    # "BG_OTHER" : {"name" : "Other",           "color" : 41,  "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    # "Data" : {"name" : "Data", "color" : ROOT.kBlack, "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 1.3},
    #"pseudoDataS" : {"name" : "pseudoDataS", "color" : ROOT.kBlack, "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 1.3},
})
