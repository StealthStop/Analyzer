#! /bin/env/python

import ROOT

from collections import OrderedDict


selections = ["1l"]

histograms = {

    "h_DoubleDisCo_disc1_?_Njets@_ABCD"   : {"logY" : False, "orders" : [str(njet) for njet in range(7,12)] + ["12incl"], "Y" : {"title" : "A.U.",   "min" : 2e-3}, "X" : {"title" : "Disc 1",                "rebin" : 10, "min" :  0,    "max" :    1}},
    "h_DoubleDisCo_disc2_?_Njets@_ABCD"   : {"logY" : False, "orders" : [str(njet) for njet in range(7,12)] + ["12incl"], "Y" : {"title" : "A.U.",   "min" : 2e-3}, "X" : {"title" : "Disc 2",                "rebin" : 10, "min" :  0,    "max" :    1}},
    #"h_DoubleDisCo_massReg_?" : {"logY" : False, "Y" : {"title" : "A.U.",   "min" : 2e-3}, "X" : {"title" : "Regression Mass [GeV]", "rebin" : 2, "min" :  0,    "max" : 1500}},

    #"h_njets_?"         : {"logY" : True,  "Y" : {"title" : "Events", "min" : 1e-1}, "X" : {"title" : "N_{jets}",        "rebin" : 1, "min" : 5.5, "max" : 20.5}},
    #"h_jetsMass_?"      : {"logY" : True,  "Y" : {"title" : "Events", "min" : 0.2},  "X" : {"title" : "Jet mass [GeV]",  "rebin" : 5, "min" :  0,   "max" :  500}},
    #"h_jetsEta_?"       : {"logY" : False, "Y" : {"title" : "Events", "min" : 0.2},  "X" : {"title" : "Jet #eta",        "rebin" : 1, "min" : -6,   "max" :    6}},
    #"h_jetsPhi_?"       : {"logY" : False, "Y" : {"title" : "Events", "min" : 0.2},  "X" : {"title" : "Jet #phi",        "rebin" : 1, "min" : -4,   "max" :    4}},
    #"h_jetsPt_?"        : {"logY" : True,  "Y" : {"title" : "Events", "min" : 0.2},  "X" : {"title" : "Jet p_{T} [GeV]", "rebin" : 5, "min" :  0,   "max" : 2000}},

    #"h_nbjets_?"        : {"logY" : True,  "Y" : {"title" : "Events", "min" : 2e-3}, "X" : {"title" : "N_{bjets}",        "rebin" : 1, "min" : -0.5, "max" : 20.5}},
    #"h_bjetsMass_?"     : {"logY" : True,  "Y" : {"title" : "Events", "min" : 0.2},  "X" : {"title" : "BJet mass [GeV]",  "rebin" : 5, "min" :  0,   "max" :  500}},
    #"h_bjetsEta_?"      : {"logY" : False, "Y" : {"title" : "Events", "min" : 0.2},  "X" : {"title" : "BJet #eta",        "rebin" : 1, "min" : -6,   "max" :    6}},
    #"h_bjetsPhi_?"      : {"logY" : False, "Y" : {"title" : "Events", "min" : 0.2},  "X" : {"title" : "BJet #phi",        "rebin" : 1, "min" : -4,   "max" :    4}},
    #"h_bjetsPt_?"       : {"logY" : True,  "Y" : {"title" : "Events", "min" : 0.2},  "X" : {"title" : "BJet p_{T} [GeV]", "rebin" : 5, "min" :  0,   "max" : 2000}},

    #"h_ntops_?"         : {"logY" : True,  "Y" : {"title" : "Events", "min" : 1e-1}, "X" : {"title" : "N_{tops}",        "rebin" : 1,  "min" :  -0.5, "max" : 10.5}},
    #"h_nRtops_?"        : {"logY" : True,  "Y" : {"title" : "Events", "min" : 2e-3}, "X" : {"title" : "N_{Res. tops}",   "rebin" : 1,  "min" :  -0.5, "max" : 10.5}},
    #"h_nMtops_?"        : {"logY" : True,  "Y" : {"title" : "Events", "min" : 2e-3}, "X" : {"title" : "N_{Mer. tops}",   "rebin" : 1,  "min" :  -0.5, "max" : 10.5}},
    #"h_topsMass_?"      : {"logY" : True,  "Y" : {"title" : "Events", "min" : 0.2},  "X" : {"title" : "Top mass [GeV]",  "rebin" : 5,  "min" : 50,    "max" :  350}},
    #"h_topsEta_?"       : {"logY" : False, "Y" : {"title" : "Events", "min" : 0.2},  "X" : {"title" : "Top #eta",        "rebin" : 1,  "min" : -6,    "max" :    6}},
    #"h_topsPhi_?"       : {"logY" : False, "Y" : {"title" : "Events", "min" : 0.2},  "X" : {"title" : "Top #phi",        "rebin" : 1,  "min" : -4,    "max" :    4}},
    #"h_topsPt_?"        : {"logY" : True,  "Y" : {"title" : "Events", "min" : 0.2},  "X" : {"title" : "Top p_{T} [GeV]", "rebin" : 10, "min" :  0,    "max" : 1500}},
    #
    #"h_ht_?"            : {"logY" : True,  "Y" : {"title" : "Events", "min" : 2e-3}, "X" : {"title" : "H_{T} [GeV]",                "rebin" : 1, "min" :  0, "max" : 3000}},  
    #"h_met_?"           : {"logY" : True,  "Y" : {"title" : "Events", "min" : 2e-3}, "X" : {"title" : "MET [GeV]",                  "rebin" : 1, "min" :  0, "max" : 2000}},  
    #"h_dR_bjets_?"      : {"logY" : False, "Y" : {"title" : "Events", "min" : 2e-3}, "X" : {"title" : "#DeltaR_{bjets} [GeV]",      "rebin" : 1, "min" :  0, "max" :   10}},  
    #"h_dR_top1_top2_?"  : {"logY" : False, "Y" : {"title" : "Events", "min" : 2e-3}, "X" : {"title" : "#DeltaR_{top1,top2} [GeV]",  "rebin" : 1, "min" :  0, "max" :   10}},  
    #"h_dR_tops_bjets_?" : {"logY" : False, "Y" : {"title" : "Events", "min" : 2e-3}, "X" : {"title" : "#DeltaR_{tops,bjets} [GeV]", "rebin" : 1, "min" :  0, "max" :   10}},  
}

backgrounds = {
    "TT"              : {"name" : "t#bar{t} + jets",                   "loption" : "F", "color" : 40,   "lstyle" : 1, "mstyle" : 8, "lsize" : 6, "msize" : 0},
    #"QCD"             : {"name" : "QCD multijet",                      "loption" : "L", "color" : 30,   "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},

    #"WJets"           : {"name" : "W + jets",                          "loption" : "L", "color" : ROOT.TColor.GetColor("#fb8072"),   "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},
    #"DYJetsToLL_M-50" : {"name" : "DY + jets",                         "loption" : "L", "color" : ROOT.TColor.GetColor("#80b1d3"),   "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},
    #"Diboson"         : {"name" : "Diboson",                           "loption" : "L", "color" : ROOT.TColor.GetColor("#fdb462"),   "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},
    #"Triboson"        : {"name" : "Triboson",                          "loption" : "L", "color" : ROOT.TColor.GetColor("#ffcc33"),   "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},

    "TTZ"             : {"name" : "t#bar{t}Z",                         "loption" : "L",  "color" : ROOT.TColor.GetColor("#fb8072"),   "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},
    "ttHJets"         : {"name" : "t#bar{t}H + jets",                  "loption" : "L",  "color" : ROOT.TColor.GetColor("#80b1d3"),   "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},
    "TTWJets"         : {"name" : "t#bar{t}W + jets",                  "loption" : "L",  "color" : ROOT.TColor.GetColor("#fdb462"),   "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},
    "TTTT"            : {"name" : "t#bar{t}t#bar{t}",                  "loption" : "L",  "color" : ROOT.TColor.GetColor("#b3de69"),   "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},
    "TTWW"            : {"name" : "t#bar{t} + WW",                     "loption" : "L",  "color" : ROOT.TColor.GetColor("#ffcc33"),   "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},
    "TTWZ"            : {"name" : "t#bar{t} + WZ",                     "loption" : "L",  "color" : ROOT.TColor.GetColor("#fccde5"),   "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},
    "TTTW"            : {"name" : "t#bar{t} + tW",                     "loption" : "L",  "color" : ROOT.TColor.GetColor("#8dd3c7"),   "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},
    "TTZH"            : {"name" : "t#bar{t} + ZH",                     "loption" : "L",  "color" : 28,   "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},
    "TTWH"            : {"name" : "t#bar{t} + WH",                     "loption" : "L",  "color" : 22,   "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},
    "TTZZ"            : {"name" : "t#bar{t} + ZZ",                     "loption" : "L",  "color" : 13,   "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},
    "TTHH"            : {"name" : "t#bar{t} + HH",                     "loption" : "L",  "color" : 9,   "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},
    "TTTJ"            : {"name" : "t#bar{t} + t + jets",               "loption" : "L",  "color" : 49,    "lstyle" : 1, "mstyle" : 8, "lsize" : 4, "msize" : 0},
}

signals = OrderedDict({
})

data = {
}
