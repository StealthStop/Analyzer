#! /bin/env/python

import ROOT

from collections import OrderedDict


selections = ["0l_blind_ABCD", "1l_blind_ABCD", "2l_blind_ABCD", "0l_QCDCR_ABCD", "1l_QCDCR_ABCD", "2l_QCDCR_ABCD"] 

histograms = {

    # jets
    "h_Njets_?"    : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "N_{ jets}",        "rebin" : 1, "min" : -0.5, "max" : 20.5}},

    # bjets
    "h_Nbjets_?"    : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "N_{ b jets}",       "rebin" : 1, "min" : -0.5, "max" : 20.5}},

    # tops
    "h_Ntops_?"     : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "N_{ tops}",                         "rebin" : 1,  "min" : -0.5,  "max" :  6.5}},
    "h_Top1_Mass_?" : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Leading p_{T} Top mass [GeV]",     "rebin" : 1,  "min" :  50,   "max" :  350}},
    "h_Top1_Pt_?"   : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Leading p_{T} Top p_{T} [GeV]",    "rebin" : 20, "min" :  0,    "max" : 1800}},
    "h_Top1_Eta_?"  : {"logY" : False, "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Leading p_{T} Top #eta",           "rebin" : 1,  "min" : -6,    "max" :    6}},
    "h_Top1_Phi_?"  : {"logY" : False, "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Leading p_{T} Top #phi",           "rebin" : 1,  "min" : -4,    "max" :    4}},
    "h_Top2_Mass_?" : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Subleading p_{T} Top mass [GeV]",  "rebin" : 1,  "min" :  50,   "max" :  350}},
    "h_Top2_Pt_?"   : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Subleading p_{T} Top p_{T} [GeV]", "rebin" : 20, "min" :  0,    "max" : 1800}},
    "h_Top2_Eta_?"  : {"logY" : False, "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Subleading p_{T} Top #eta",        "rebin" : 1,  "min" : -6,    "max" :    6}},
    "h_Top2_Phi_?"  : {"logY" : False, "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Subleading p_{T} Top #phi",        "rebin" : 1,  "min" : -4,    "max" :    4}},

    # others
    "h_dRbjets_?"  : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "#DeltaR_{BJets}", "rebin" : 5,  "min" : 0, "max" :    6}},
    "h_HT_?"       : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "H_{T} [GeV]",     "rebin" : 10, "min" : 0, "max" : 3500}},
    "h_Mbl_?"      : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "M_{b,l} [GeV]",   "rebin" : 3,  "min" : 0, "max" :  360}}, 
    "h_Mll_?"      : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "M_{l,l} [GeV]",   "rebin" : 3,  "min" : 0, "max" :  360}},
    "h_Jet1_Flavb_?" : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Jet 1 Flavb",    "rebin" : 1, "min" :  0,    "max" : 1}},
    "h_Jet2_Flavb_?" : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Jet 2 Flavb",    "rebin" : 1, "min" :  0,    "max" : 1}},
    "h_Jet3_Flavb_?" : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Jet 3 Flavb",    "rebin" : 1, "min" :  0,    "max" : 1}},
    "h_Jet4_Flavb_?" : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Jet 4 Flavb",    "rebin" : 1, "min" :  0,    "max" : 1}},
    "h_Jet5_Flavb_?" : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Jet 5 Flavb",    "rebin" : 1, "min" :  0,    "max" : 1}},
    "h_Jet6_Flavb_?" : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Jet 6 Flavb",    "rebin" : 1, "min" :  0,    "max" : 1}},


    "h_njets_10incl_SYY_2l_QCDCR_ABCD"    : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "N_{jets} in each A,B,C,D region",        "rebin" : 1, "min" : -0.5, "max" : 23.5}},
    "h_njets_11incl_SYY_1l_QCDCR_ABCD"    : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "N_{jets} in each A,B,C,D region",        "rebin" : 1, "min" : -0.5, "max" : 23.5}},
    "h_njets_12incl_SYY_0l_QCDCR_ABCD"    : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "N_{jets} in each A,B,C,D region",        "rebin" : 1, "min" : -0.5, "max" : 23.5}},

    "h_njets_10incl_RPV_2l_QCDCR_ABCD"    : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "N_{jets} in each A,B,C,D region",        "rebin" : 1, "min" : -0.5, "max" : 23.5}},
    "h_njets_11incl_RPV_1l_QCDCR_ABCD"    : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "N_{jets} in each A,B,C,D region",        "rebin" : 1, "min" : -0.5, "max" : 23.5}},
    "h_njets_12incl_RPV_0l_QCDCR_ABCD"    : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "N_{jets} in each A,B,C,D region",        "rebin" : 1, "min" : -0.5, "max" : 23.5}},

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

signals = OrderedDict([
    ("RPV_2t6j_mStop-300",        {"name" : "RPV m_{ #tilde{t}} = 300 GeV",               "color" : 2, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0}),
    ("RPV_2t6j_mStop-800",        {"name" : "RPV m_{ #tilde{t}} = 800 GeV",               "color" : 7, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0}),
    #("StealthSYY_2t6j_mStop-300", {"name" : "Stealth SY#bar{Y} m_{ #tilde{t}} = 300 GeV", "color" : 4, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0}),
    #("StealthSYY_2t6j_mStop-800", {"name" : "Stealth SY#bar{Y} m_{ #tilde{t}} = 800 GeV", "color" : 6, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0})
])
