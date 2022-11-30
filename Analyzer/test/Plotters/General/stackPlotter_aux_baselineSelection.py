#! /bin/env/python

import ROOT

from collections import OrderedDict


selections = ["0l_blind_ABCD", "1l_blind_ABCD", "2l_blind_ABCD", "0l_QCDCR_ABCD", "1l_QCDCR_ABCD"] 

histograms = {

    # jets
    "h_Njets_?"    : {"logY" : False, "Y" : {"title" : "Events", "min" : 0.2}, "X" : {"title" : "N_{Jets}",        "rebin" : 1, "min" : -0.5, "max" : 20.5}},
    #"h_jetsMass_?" : {"logY" : True,  "Y" : {"title" : "Events", "min" : 0.2}, "X" : {"title" : "Jet mass [GeV]",  "rebin" : 5, "min" :  0,   "max" :  500}},
    #"h_jetsPt_?"   : {"logY" : True,  "Y" : {"title" : "Events", "min" : 0.2}, "X" : {"title" : "Jet p_{T} [GeV]", "rebin" : 5, "min" :  0,   "max" : 2000}},
    #"h_jetsEta_?"  : {"logY" : False, "Y" : {"title" : "Events", "min" : 0.2}, "X" : {"title" : "Jet #eta",        "rebin" : 1, "min" : -6,   "max" :    6}},
    #"h_jetsPhi_?"  : {"logY" : False, "Y" : {"title" : "Events", "min" : 0.2}, "X" : {"title" : "Jet #phi",        "rebin" : 1, "min" : -4,   "max" :    4}},

    # bjets
    "h_Nbjets_?"    : {"logY" : False, "Y" : {"title" : "Events", "min" : 0.2}, "X" : {"title" : "N B_{Jets}",       "rebin" : 1, "min" : -0.5, "max" : 20.5}},
    #"h_bjetsMass_?" : {"logY" : True,  "Y" : {"title" : "Events", "min" : 0.2}, "X" : {"title" : "BJet mass [GeV]",  "rebin" : 5, "min" :  0,   "max" :  500}},
    #"h_bjetsPt_?"   : {"logY" : True,  "Y" : {"title" : "Events", "min" : 0.2}, "X" : {"title" : "BJet p_{T} [GeV]", "rebin" : 5, "min" :  0,   "max" : 2000}},
    #"h_bjetsEta_?"  : {"logY" : False, "Y" : {"title" : "Events", "min" : 0.2}, "X" : {"title" : "BJet #eta",        "rebin" : 1, "min" : -6,   "max" :    6}},
    #"h_bjetsPhi_?"  : {"logY" : False, "Y" : {"title" : "Events", "min" : 0.2}, "X" : {"title" : "BJet #phi",        "rebin" : 1, "min" : -4,   "max" :    4}},

    # tops
    "h_Ntops_?"     : {"logY" : False, "Y" : {"title" : "Events", "min" : 0.2}, "X" : {"title" : "N_{Tops}",        "rebin" : 1,  "min" : -0.5,  "max" : 6.5 }},
    "h_topsMass_?" : {"logY" : True,  "Y" : {"title" : "Events", "min" : 0.2}, "X" : {"title" : "Top mass [GeV]",  "rebin" : 10, "min" :  50,   "max" :  350}},
    "h_topsPt_?"   : {"logY" : True,  "Y" : {"title" : "Events", "min" : 0.2}, "X" : {"title" : "Top p_{T} [GeV]", "rebin" : 20, "min" :  0,    "max" : 1500}},
    "h_topsEta_?"  : {"logY" : False, "Y" : {"title" : "Events", "min" : 0.2}, "X" : {"title" : "Top #eta",        "rebin" : 1,  "min" : -6,    "max" :    6}},
    "h_topsPhi_?"  : {"logY" : False, "Y" : {"title" : "Events", "min" : 0.2}, "X" : {"title" : "Top #phi",        "rebin" : 1,  "min" : -4,    "max" :    4}},

    # others
    "h_dRbjets_?"    : {"logY" : False, "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "#Delta R B_{Jets}", "rebin" : 1,  "min" : 0,    "max" :    6}},
    "h_HT_?"       : {"logY" : True,  "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "H_{T} [GeV]",       "rebin" : 20, "min" : 0,    "max" : 3500}},
    "h_Mbl_?"      : {"logY" : True,  "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "M_{b,l} [GeV]",     "rebin" : 10, "min" : 0,    "max" : 360 }}, 
    "h_Mll_?"      : {"logY" : True,  "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "M_{l,l} [GeV]",     "rebin" : 10, "min" : 0,    "max" : 360 }},

}

data = {
    "Data_skim" : {"name" : "Data", "color" : ROOT.kBlack, "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 1.3}
} 

backgrounds = {
    "TT_skim"       : {"name" : "t#bar{t} + jets", "color" : 40,  "lstyle" : 1, "mstyle" : 8, "lsize" : 0, "msize" : 0},
    "QCD_skim"      : {"name" : "QCD multijet",    "color" : 30,  "lstyle" : 1, "mstyle" : 8, "lsize" : 0, "msize" : 0},
    "TTX_skim"      : {"name" : "t#bar{t} + X",    "color" : 38,  "lstyle" : 1, "mstyle" : 8, "lsize" : 0, "msize" : 0},
    "BG_OTHER_skim" : {"name" : "Other",           "color" : 41,  "lstyle" : 1, "mstyle" : 8, "lsize" : 0, "msize" : 0},
    #"DYJetsToLL_M-50" : {"name" : "DY + jets",       "color" : ROOT.TColor.GetColor("#80b1d3"), "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},
    #"Diboson"         : {"name" : "Diboson",         "color" : ROOT.TColor.GetColor("#fdb462"), "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},
    #"Triboson"        : {"name" : "Triboson",        "color" : ROOT.TColor.GetColor("#ffcc33"), "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},   
    #"WJets"           : {"name" : "W + jets",        "color" : ROOT.TColor.GetColor("#fb8072"), "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},
}

signals = OrderedDict({
    "RPV_2t6j_mStop-350"        : {"name" : "RPV m_{ #tilde{t}} = 350 GeV",               "color" : 2, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "RPV_2t6j_mStop-550"        : {"name" : "RPV m_{ #tilde{t}} = 550 GeV",               "color" : 7, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "RPV_2t6j_mStop-850"        : {"name" : "RPV m_{ #tilde{t}} = 850 GeV",               "color" : 4, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "RPV_2t6j_mStop-1050"       : {"name" : "RPV m_{ #tilde{t}} = 1050 GeV",              "color" : 5, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "StealthSYY_2t6j_mStop-900" : {"name" : "Stealth SY#bar{Y} m_{ #tilde{t}} = 900 GeV", "color" : 6, "lstyle" : 3, "mstyle" : 8, "lsize" : 3, "msize" : 0}
})

