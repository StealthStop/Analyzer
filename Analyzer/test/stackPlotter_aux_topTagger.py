#! /bin/env/python

import ROOT

from collections import OrderedDict

selections = [
        "",
        "_1mu",
        "_1mu_7j",
        "_1mu_7j_0el",
        "_1mu_7j_0el_1b",
        "_1mu_7j_0el_1b_jMETdPhi",
        "_1mu_7j_0el_1b_jMETdPhi_1bl",
        "_1mu_7j_0el_1b_jMETdPhi_1bl_mubdR_mubM",
        "_1mu_7j_0el_1b_jMETdPhi_1bl_mubdR_mubM_muMETphi",
        "_1mu_7j_0el_1b_jMETdPhi_1bl_mubdR_mubM_muMETphi_muMETmT",
        "_1mu_7j_0el_1b_jMETdPhi_1bl_mubdR_mubM_muMETphi_muMETmT_HT500",
        "_1mu_7j_0el_1b_jMETdPhi_1bl_mubdR_mubM_muMETphi_muMETmT_HT500_MET50",
        "_TTbarCR",
        "_QCDCR"
]

histograms = {
    "h_jet1METdPhi?_@_Njet7incl"  : {"logY" : False, "orders" : ["num", "den"], "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "#Delta#phi(jet 1, MET)",          "rebin" :  1, "min" :  -4.0, "max" :  4.0}},
    "h_jet2METdPhi?_@_Njet7incl"  : {"logY" : False, "orders" : ["num", "den"], "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "#Delta#phi(jet 2, MET)",          "rebin" :  1, "min" :  -4.0, "max" :  4.0}},
    "h_jet3METdPhi?_@_Njet7incl"  : {"logY" : False, "orders" : ["num", "den"], "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "#Delta#phi(jet 3, MET)",          "rebin" :  1, "min" :  -4.0, "max" :  4.0}},
    "h_met?_@_Njet7incl"          : {"logY" : True,  "orders" : ["num", "den"], "Y" : {"title" : "Weighted Events", "min" : 2e-4},  "X" : {"title" : "E_{T}^{miss} [GeV]",             "rebin" :  6, "min" :     0, "max" : 1500}},
    "h_ht?_@_Njet7incl"           : {"logY" : True,  "orders" : ["num", "den"], "Y" : {"title" : "Weighted Events", "min" : 2e-4},  "X" : {"title" : "H_{T} [GeV]",                    "rebin" : 15, "min" :     0, "max" : 5000}},
    "h_muonBjetDR?_@_Njet7incl"   : {"logY" : False, "orders" : ["num", "den"], "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "min #DeltaR(muon, loose b jet)",  "rebin" :  3, "min" :     0, "max" :    5}},
    "h_muonBjetMass?_@_Njet7incl" : {"logY" : False, "orders" : ["num", "den"], "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "m(muon, loose b jet) [GeV]",      "rebin" : 4, "min" :     0, "max" : 500}},
    "h_muonMETdPhi?_@_Njet7incl"  : {"logY" : False, "orders" : ["num", "den"], "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "#Delta#phi(muon, E_{T}^{miss})",  "rebin" :  1, "min" :  -4.0, "max" :  4.0}},
    "h_muonMETmT?_@_Njet7incl"    : {"logY" : False, "orders" : ["num", "den"], "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "m_{T}(muon, E_{T}^{miss}) [GeV]", "rebin" :  2, "min" :     0, "max" : 200}},
    "h_nJets?_@_Njet7incl"        : {"logY" : True,  "orders" : ["num", "den"], "Y" : {"title" : "Events",          "min" : 2e-4}, "X" : {"title" : "N_{jets}",                        "rebin" :  1, "min" :  -0.5, "max" : 20.5}},
    "h_nBJets?_@_Njet7incl"       : {"logY" : True,  "orders" : ["num", "den"], "Y" : {"title" : "Events",          "min" : 2e-4}, "X" : {"title" : "N_{b jets}",                      "rebin" :  1, "min" :  -0.5, "max" : 10.5}},
    "h_nElectrons?_@_Njet7incl"   : {"logY" : True,  "orders" : ["num", "den"], "Y" : {"title" : "Events",          "min" : 2e-4}, "X" : {"title" : "N_{electrons}",                   "rebin" :  1, "min" :  -0.5, "max" : 20.5}},
    "h_nMuons?_@_Njet7incl"       : {"logY" : True,  "orders" : ["num", "den"], "Y" : {"title" : "Events",          "min" : 2e-4}, "X" : {"title" : "N_{muons}",                       "rebin" :  1, "min" :  -0.5, "max" : 20.5}},
    "h_topPt?_num_Njet@"        : {"logY" : True,  "orders" : ["7", "8", "9", "10", "11", "12incl", "7incl"], "Y" : {"title" : "Weighted Events", "min" : 2e-4},  "X" : {"title" : "Best Candidate top p_{T} [GeV]",           "rebin" : 12, "min" :     0, "max" : 2000}},
    "h_topPt?_den_Njet@"        : {"logY" : True,  "orders" : ["7", "8", "9", "10", "11", "12incl", "7incl"], "Y" : {"title" : "Weighted Events", "min" : 2e-4},  "X" : {"title" : "Best Candidate top p_{T} [GeV]",           "rebin" : 12, "min" :     0, "max" : 2000}},


}

backgrounds = {
    "TT"       : {"name" : "t#bar{t} + jets", "color" : 40,  "lstyle" : 1, "mstyle" : 8, "lsize" : 0, "msize" : 0},
    "QCD"      : {"name" : "QCD multijet",    "color" : 30,  "lstyle" : 1, "mstyle" : 8, "lsize" : 0, "msize" : 0},
    "TTX"      : {"name" : "t#bar{t} + X",    "color" : 38,  "lstyle" : 1, "mstyle" : 8, "lsize" : 0, "msize" : 0},
    "BG_OTHER" : {"name" : "Other",           "color" : 41,  "lstyle" : 1, "mstyle" : 8, "lsize" : 0, "msize" : 0}
}

signals = OrderedDict({
    #"RPV_2t6j_mStop-350"        : {"name" : "RPV m_{ #tilde{t}} = 350 GeV",               "color" : 2, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    #"RPV_2t6j_mStop-550"        : {"name" : "RPV m_{ #tilde{t}} = 550 GeV",               "color" : 7, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    #"RPV_2t6j_mStop-850"        : {"name" : "RPV m_{ #tilde{t}} = 850 GeV",               "color" : 4, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
})

data = {
    #"Data" : {"name" : "Data", "color" : ROOT.kBlack, "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 1.3}
    #"pseudoDataS" : {"name" : "pseudoDataS", "color" : ROOT.kBlack, "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 1.3}
}
