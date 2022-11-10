#! /bin/env/python

import ROOT

from collections import OrderedDict

selections = [
        "",
        "_1mu",
        "_1mu_ge6j",
        "_1mu_ge6j_0el",
        "_1mu_ge6j_0el_1b",
        "_1mu_ge6j_0el_1b_jMETdPhi",
        "_1mu_ge6j_0el_1b_jMETdPhi_1bl",
        "_1mu_ge6j_0el_1b_jMETdPhi_1bl_mubdR_mubM",
        "_1mu_ge6j_0el_1b_jMETdPhi_1bl_mubdR_mubM_muMETphi",
        "_1mu_ge6j_0el_1b_jMETdPhi_1bl_mubdR_mubM_muMETphi_muMETmT",
        "_1mu_ge6j_0el_1b_jMETdPhi_1bl_mubdR_mubM_muMETphi_muMETmT_HT200_pt30",
        "_1mu_ge6j_0el_1b_jMETdPhi_1bl_mubdR_mubM_muMETphi_muMETmT_HT200_pt30_MET50",
        "_TTbarCR",
        #"_QCDCR"
]

histograms = {
    "jet1METdPhi?"  : {"logY" : False, "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "#Delta#phi(jet 1, MET)",          "rebin" :  2, "min" :  -4.0, "max" :  4.0}},
    "jet2METdPhi?"  : {"logY" : False, "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "#Delta#phi(jet 2, MET)",          "rebin" :  2, "min" :  -4.0, "max" :  4.0}},
    "jet3METdPhi?"  : {"logY" : False, "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "#Delta#phi(jet 3, MET)",          "rebin" :  2, "min" :  -4.0, "max" :  4.0}},
    "met?"          : {"logY" : True,  "Y" : {"title" : "Weighted Events", "min" : 2e-3}, "X" : {"title" : "E_{T}^{miss} [GeV]",              "rebin" :  6, "min" :     0, "max" : 1500}},
    "ht?"           : {"logY" : True,  "Y" : {"title" : "Weighted Events", "min" : 2e-3}, "X" : {"title" : "H_{T} [GeV]",                     "rebin" : 15, "min" :     0, "max" : 5000}},
    "muonBjetDR?"   : {"logY" : False, "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "min #DeltaR(muon, loose b jet)",  "rebin" :  3, "min" :     0, "max" :    5}},
    "muonBjetMass?" : {"logY" : False, "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "m(muon, loose b jet) [GeV]",      "rebin" :  6, "min" :     0, "max" : 500}},
    "muonMETdPhi?"  : {"logY" : False, "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "#Delta#phi(muon, E_{T}^{miss})",  "rebin" :  2, "min" :  -4.0, "max" :  4.0}},
    "muonMETmT?"    : {"logY" : False, "Y" : {"title" : "Weighted Events", "min" : 0.2},  "X" : {"title" : "m_{T}(muon, E_{T}^{miss}) [GeV]", "rebin" :  3, "min" :     0, "max" : 200}},
    "nJets?"        : {"logY" : True,  "Y" : {"title" : "Weighted Events", "min" : 2e-3}, "X" : {"title" : "N_{jets}",                        "rebin" :  1, "min" :  -0.5, "max" : 20.5}},
    "nBJets?"       : {"logY" : True,  "Y" : {"title" : "Weighted Events", "min" : 2e-3}, "X" : {"title" : "N_{b jets}",                      "rebin" :  1, "min" :  -0.5, "max" : 10.5}},
    "nElectrons?"   : {"logY" : True,  "Y" : {"title" : "Weighted Events", "min" : 2e-3}, "X" : {"title" : "N_{electrons}",                   "rebin" :  1, "min" :  -0.5, "max" : 20.5}},
    "nMuons?"       : {"logY" : True,  "Y" : {"title" : "Weighted Events", "min" : 2e-3}, "X" : {"title" : "N_{muons}",                       "rebin" :  1, "min" :  -0.5, "max" : 20.5}},
    "topPt?"        : {"logY" : True,  "Y" : {"title" : "Weighted Events", "min" : 2e-3}, "X" : {"title" : "Best Candidate top p_{T} [GeV]",  "rebin" : 12, "min" :     0, "max" : 2000}},
    "topPt?"        : {"logY" : True,  "Y" : {"title" : "Weighted Events", "min" : 2e-3}, "X" : {"title" : "Best Candidate top p_{T} [GeV]",  "rebin" : 12, "min" :     0, "max" : 2000}},


}

backgrounds = {
    "TT"             : {"name" : "t#bar{t} + jets", "color" : 40,  "lstyle" : 1, "mstyle" : 8, "lsize" : 0, "msize" : 0},
    "QCD"            : {"name" : "QCD multijet",    "color" : 30,  "lstyle" : 1, "mstyle" : 8, "lsize" : 0, "msize" : 0},
    "WJets"          : {"name" : "W + jets",        "color" : ROOT.TColor.GetColor("#fb8072"),   "lstyle" : 1, "mstyle" : 8, "lsize" : 0, "msize" : 0},
    "DYJetsToLL_M-50": {"name" : "DY + jets",       "color" : ROOT.TColor.GetColor("#80b1d3"),   "lstyle" : 1, "mstyle" : 8, "lsize" : 0, "msize" : 0},
    "Diboson"        : {"name" : "Diboson",         "color" : ROOT.TColor.GetColor("#fdb462"),   "lstyle" : 1, "mstyle" : 8, "lsize" : 0, "msize" : 0},
    "Triboson"       : {"name" : "Triboson",        "color" : ROOT.TColor.GetColor("#ffcc33"),   "lstyle" : 1, "mstyle" : 8, "lsize" : 0, "msize" : 0},
    "TTX"            : {"name" : "t#bar{t} + X",    "color" : 38,  "lstyle" : 1, "mstyle" : 8, "lsize" : 0, "msize" : 0},
    #"BG_OTHER"      : {"name" : "Other",           "color" : 41,  "lstyle" : 1, "mstyle" : 8, "lsize" : 0, "msize" : 0}
}

signals = OrderedDict({
    #"RPV_2t6j_mStop-350"        : {"name" : "RPV m_{ #tilde{t}} = 350 GeV",               "color" : 2, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    #"RPV_2t6j_mStop-550"        : {"name" : "RPV m_{ #tilde{t}} = 550 GeV",               "color" : 7, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    #"RPV_2t6j_mStop-850"        : {"name" : "RPV m_{ #tilde{t}} = 850 GeV",               "color" : 4, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
})

data = {
    "Data_SingleMuon" : {"name" : "Data", "color" : ROOT.kBlack, "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 1.3}
    #"pseudoDataS" : {"name" : "pseudoDataS", "color" : ROOT.kBlack, "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 1.3}
}
