#! /bin/env/python

import ROOT

from collections import OrderedDict

selections = [
    "_1l",
    "_1l_QCDCR",
    #"_0l_HT500_0NonIsoMuon_ge7j_ge2t_ge1dRbjets",
    #"_0l_HT500_0NonIsoMuon_7j_ge2t_ge1dRbjets",
    #"_0l_HT500_0NonIsoMuon_8j_ge2t_ge1dRbjets",
    #"_0l_HT500_0NonIsoMuon_9j_ge2t_ge1dRbjets",
    #"_0l_HT500_0NonIsoMuon_10j_ge2t_ge1dRbjets",
    #"_0l_HT500_0NonIsoMuon_11j_ge2t_ge1dRbjets",
    #"_0l_HT500_0NonIsoMuon_ge12j_ge2t_ge1dRbjets",
    #"_0l_HT500_ge7j_ge2t_ge1dRbjets",
    #"_0l_HT500_7j_ge2t_ge1dRbjets",
    #"_0l_HT500_8j_ge2t_ge1dRbjets",
    #"_0l_HT500_9j_ge2t_ge1dRbjets",
    #"_0l_HT500_10j_ge2t_ge1dRbjets",
    #"_0l_HT500_11j_ge2t_ge1dRbjets",
    #"_0l_HT500_ge12j_ge2t_ge1dRbjets",
    #"_0l_HT500_ge7j_ge2b_ge2t_ge1dRbjets",
    #"_0l_HT500_7j_ge2b_ge2t_ge1dRbjets",
    #"_0l_HT500_8j_ge2b_ge2t_ge1dRbjets",
    #"_0l_HT500_9j_ge2b_ge2t_ge1dRbjets",
    #"_0l_HT500_10j_ge2b_ge2t_ge1dRbjets",
    #"_0l_HT500_11j_ge2b_ge2t_ge1dRbjets",
    #"_0l_HT500_ge12j_ge2b_ge2t_ge1dRbjets",
]

histograms = {
    #"h_cm_pt_jetRank_1l_QCDCR" : {"logX" : False, "logY" : False, "logZ" : False, "Y" : {"title" : "Number of Resolved Tops", "min" : -0.5, "max" : 5.5}, "X" : {"title" : "Number of Merged Tops", "min" :  -0.5, "max" : 5.5}},
    "h_cm_pt_jetRank?"    : {"logX" : False, "logY" : False, "logZ" : False, "Y" : {"title" : "Jet Rank by pT", "min" : 0.0, "max" : 10.0}, "X" : {"title" : "Jet pT [GeV]", "min" :  0.0, "max" : 1500.0, "rebin" : 5}},
    "h_cm_ptrHT_jetRank?" : {"logX" : False, "logY" : False, "logZ" : False, "Y" : {"title" : "Jet Rank by pT", "min" : 0.0, "max" : 10.0}, "X" : {"title" : "Jet pT / HT",  "min" :  0.0, "max" :    1.0, "rebin" : 5}},
}

backgrounds = {
    "AllBG"    : {"name" : "AllBg", "option" : "COLZ"},
    #"TT"       : {"name" : "t#bar{t} + jets", "option" : "COLZ"},
    #"QCD"      : {"name" : "QCD multijet",    "option" : "COLZ"},
    #"TTX"      : {"name" : "t#bar{t} + X",    "option" : "COLZ"},
    #"BG_OTHER" : {"name" : "Other",           "option" : "COLZ"}
}

signals = OrderedDict({
    #"RPV_2t6j_mStop-300"  : {"name" : "RPV m_{ #tilde{t}} = 300 GeV",  "option" : "COLZ E TEXT"},
    #"RPV_2t6j_mStop-350"  : {"name" : "RPV m_{ #tilde{t}} = 350 GeV",  "option" : "COLZ E TEXT"},
    #"RPV_2t6j_mStop-400"  : {"name" : "RPV m_{ #tilde{t}} = 400 GeV",  "option" : "COLZ E TEXT"},
    #"RPV_2t6j_mStop-450"  : {"name" : "RPV m_{ #tilde{t}} = 450 GeV",  "option" : "COLZ E TEXT"},
    #"RPV_2t6j_mStop-500"  : {"name" : "RPV m_{ #tilde{t}} = 500 GeV",  "option" : "COLZ E TEXT"},
    #"RPV_2t6j_mStop-550"  : {"name" : "RPV m_{ #tilde{t}} = 550 GeV",  "option" : "COLZ E TEXT"},
    #"RPV_2t6j_mStop-600"  : {"name" : "RPV m_{ #tilde{t}} = 600 GeV",  "option" : "COLZ E TEXT"},
    #"RPV_2t6j_mStop-650"  : {"name" : "RPV m_{ #tilde{t}} = 650 GeV",  "option" : "COLZ E TEXT"},
    #"RPV_2t6j_mStop-700"  : {"name" : "RPV m_{ #tilde{t}} = 700 GeV",  "option" : "COLZ E TEXT"},
    #"RPV_2t6j_mStop-750"  : {"name" : "RPV m_{ #tilde{t}} = 750 GeV",  "option" : "COLZ E TEXT"},
    #"RPV_2t6j_mStop-800"  : {"name" : "RPV m_{ #tilde{t}} = 800 GeV",  "option" : "COLZ E TEXT"},
    #"RPV_2t6j_mStop-850"  : {"name" : "RPV m_{ #tilde{t}} = 850 GeV",  "option" : "COLZ E TEXT"},
    #"RPV_2t6j_mStop-900"  : {"name" : "RPV m_{ #tilde{t}} = 900 GeV",  "option" : "COLZ E TEXT"},
    #"RPV_2t6j_mStop-950"  : {"name" : "RPV m_{ #tilde{t}} = 950 GeV",  "option" : "COLZ E TEXT"},
    #"RPV_2t6j_mStop-1000" : {"name" : "RPV m_{ #tilde{t}} = 1000 GeV",  "option" : "COLZ E TEXT"},
    #"RPV_2t6j_mStop-1050" : {"name" : "RPV m_{ #tilde{t}} = 1050 GeV",  "option" : "COLZ E TEXT"},
    #"RPV_2t6j_mStop-1100" : {"name" : "RPV m_{ #tilde{t}} = 1100 GeV",  "option" : "COLZ E TEXT"},
    #"RPV_2t6j_mStop-1150" : {"name" : "RPV m_{ #tilde{t}} = 1150 GeV",  "option" : "COLZ E TEXT"},
    #"RPV_2t6j_mStop-1200" : {"name" : "RPV m_{ #tilde{t}} = 1200 GeV",  "option" : "COLZ E TEXT"},
    #"RPV_2t6j_mStop-1250" : {"name" : "RPV m_{ #tilde{t}} = 1250 GeV",  "option" : "COLZ E TEXT"},
    #"RPV_2t6j_mStop-1300" : {"name" : "RPV m_{ #tilde{t}} = 1300 GeV",  "option" : "COLZ E TEXT"},
    #"RPV_2t6j_mStop-1350" : {"name" : "RPV m_{ #tilde{t}} = 1350 GeV",  "option" : "COLZ E TEXT"},
    #"RPV_2t6j_mStop-1400" : {"name" : "RPV m_{ #tilde{t}} = 1400 GeV", "option" : "COLZ E TEXT"},
})

data = {
    "Data" : {"name" : "Data", "color" : ROOT.kBlack, "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 1.3, "option" : "COLZ"}
    #"pseudoDataS" : {"name" : "pseudoDataS", "color" : ROOT.kBlack, "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 1.3}
}
