#! /bin/env/python

import ROOT

from collections import OrderedDict

selections = [
    "",
]

histograms = {
    "h_DoubleDisCo_disc1_disc2_1l_Njets7_ABCD" :      {"logX" : False, "logY" : False, "logZ" : False, "Y" : {"title" : "Neural Network Discriminant 2", "min" : 0.0, "max" : 1.0}, "X" : {"title" : "Neural Network Discriminant 1",  "min" :  0.0, "max" :    1.0}},
    "h_DoubleDisCo_disc1_disc2_1l_Njets8_ABCD" :      {"logX" : False, "logY" : False, "logZ" : False, "Y" : {"title" : "Neural Network Discriminant 2", "min" : 0.0, "max" : 1.0}, "X" : {"title" : "Neural Network Discriminant 1",  "min" :  0.0, "max" :    1.0}},
    "h_DoubleDisCo_disc1_disc2_1l_Njets9_ABCD" :      {"logX" : False, "logY" : False, "logZ" : False, "Y" : {"title" : "Neural Network Discriminant 2", "min" : 0.0, "max" : 1.0}, "X" : {"title" : "Neural Network Discriminant 1",  "min" :  0.0, "max" :    1.0}},
    "h_DoubleDisCo_disc1_disc2_1l_Njets10_ABCD" :     {"logX" : False, "logY" : False, "logZ" : False, "Y" : {"title" : "Neural Network Discriminant 2", "min" : 0.0, "max" : 1.0}, "X" : {"title" : "Neural Network Discriminant 1",  "min" :  0.0, "max" :    1.0}},
    "h_DoubleDisCo_disc1_disc2_1l_Njets11_ABCD" :     {"logX" : False, "logY" : False, "logZ" : False, "Y" : {"title" : "Neural Network Discriminant 2", "min" : 0.0, "max" : 1.0}, "X" : {"title" : "Neural Network Discriminant 1",  "min" :  0.0, "max" :    1.0}},
    "h_DoubleDisCo_disc1_disc2_1l_Njets12incl_ABCD" : {"logX" : False, "logY" : False, "logZ" : False, "Y" : {"title" : "Neural Network Discriminant 2", "min" : 0.0, "max" : 1.0}, "X" : {"title" : "Neural Network Discriminant 1",  "min" :  0.0, "max" :    1.0}},
}

backgrounds = {
    "TT"       : {"name" : "t#bar{t} + jets", "option" : "COLZ"},
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
    "RPV_2t6j_mStop-550"  : {"name" : "RPV m_{ #tilde{t}} = 550 GeV",  "option" : "COLZ"},
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
}
