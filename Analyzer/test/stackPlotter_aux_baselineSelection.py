#! /bin/env/python

import ROOT

from collections import OrderedDict


#selections = ["0l_HT500_ge6j_ge2t_ge2b_ge1dRbjets"]

temp_list1 = ["", "_0NonIsoMuon"] 

temp_list2 = ["ge7j", "7j", "8j", "9j", "10j", "ge11j"]

temp_list3 = [
        #"0l_HT500@_?_ge2t",                   
        #"0l_HT500@_?_ge2t_ge1b",             
        #"0l_HT500@_?_ge2t_ge2b", 
        "0l_HT500@_?_ge2t_ge1dRbjets",        
]

selections = []

for list1 in temp_list1:

    for list2 in temp_list2:

        for list3 in temp_list3:

            selections.append( list3.replace("@", list1).replace("?", list2) )

histograms = {

    #"h_DoubleDisCo_disc1_?"   : {"logY" : False, "Y" : {"title" : "A.U.",   "min" : 2e-3}, "X" : {"title" : "Disc 1",                "rebin" : 1, "min" :  0,    "max" :    1}},
    #"h_DoubleDisCo_disc2_?"   : {"logY" : False, "Y" : {"title" : "A.U.",   "min" : 2e-3}, "X" : {"title" : "Disc 2",                "rebin" : 1, "min" :  0,    "max" :    1}},
    #"h_DoubleDisCo_massReg_?" : {"logY" : False, "Y" : {"title" : "A.U.",   "min" : 2e-3}, "X" : {"title" : "Regression Mass [GeV]", "rebin" : 2, "min" :  0,    "max" : 1500}},

    "h_njets_?"         : {"logY" : True,  "Y" : {"title" : "Events", "min" : 2e-3}, "X" : {"title" : "N_{jets}",        "rebin" : 1, "min" : -0.5, "max" : 20.5}},
    "h_jetsMass_?"      : {"logY" : True,  "Y" : {"title" : "Events", "min" : 0.2},  "X" : {"title" : "Jet mass [GeV]",  "rebin" : 5, "min" :  0,   "max" :  500}},
    "h_jetsEta_?"       : {"logY" : False, "Y" : {"title" : "Events", "min" : 0.2},  "X" : {"title" : "Jet #eta",        "rebin" : 1, "min" : -6,   "max" :    6}},
    "h_jetsPhi_?"       : {"logY" : False, "Y" : {"title" : "Events", "min" : 0.2},  "X" : {"title" : "Jet #phi",        "rebin" : 1, "min" : -4,   "max" :    4}},
    "h_jetsPt_?"        : {"logY" : True,  "Y" : {"title" : "Events", "min" : 0.2},  "X" : {"title" : "Jet p_{T} [GeV]", "rebin" : 5, "min" :  0,   "max" : 2000}},

    "h_nbjets_?"        : {"logY" : True,  "Y" : {"title" : "Events", "min" : 2e-3}, "X" : {"title" : "N_{bjets}",        "rebin" : 1, "min" : -0.5, "max" : 20.5}},
    #"h_bjetsMass_?"     : {"logY" : True,  "Y" : {"title" : "Events", "min" : 0.2},  "X" : {"title" : "BJet mass [GeV]",  "rebin" : 5, "min" :  0,   "max" :  500}},
    #"h_bjetsEta_?"      : {"logY" : False, "Y" : {"title" : "Events", "min" : 0.2},  "X" : {"title" : "BJet #eta",        "rebin" : 1, "min" : -6,   "max" :    6}},
    #"h_bjetsPhi_?"      : {"logY" : False, "Y" : {"title" : "Events", "min" : 0.2},  "X" : {"title" : "BJet #phi",        "rebin" : 1, "min" : -4,   "max" :    4}},
    #"h_bjetsPt_?"       : {"logY" : True,  "Y" : {"title" : "Events", "min" : 0.2},  "X" : {"title" : "BJet p_{T} [GeV]", "rebin" : 5, "min" :  0,   "max" : 2000}},

    "h_ntops_?"         : {"logY" : True,  "Y" : {"title" : "Events", "min" : 2e-3}, "X" : {"title" : "N_{tops}",        "rebin" : 1, "min" :  -0.5, "max" : 10.5}},
    #"h_nRtops_?"        : {"logY" : True,  "Y" : {"title" : "Events", "min" : 2e-3}, "X" : {"title" : "N_{Res. tops}",   "rebin" : 1, "min" :  -0.5, "max" : 10.5}},
    #"h_nMtops_?"        : {"logY" : True,  "Y" : {"title" : "Events", "min" : 2e-3}, "X" : {"title" : "N_{Mer. tops}",   "rebin" : 1, "min" :  -0.5, "max" : 10.5}},
    #"h_topsMass_?"      : {"logY" : True,  "Y" : {"title" : "Events", "min" : 0.2},  "X" : {"title" : "Top mass [GeV]",  "rebin" : 5, "min" :  0,    "max" :  500}},
    #"h_topsEta_?"       : {"logY" : False, "Y" : {"title" : "Events", "min" : 0.2},  "X" : {"title" : "Top #eta",        "rebin" : 1, "min" : -6,    "max" :    6}},
    #"h_topsPhi_?"       : {"logY" : False, "Y" : {"title" : "Events", "min" : 0.2},  "X" : {"title" : "Top #phi",        "rebin" : 1, "min" : -4,    "max" :    4}},
    #"h_topsPt_?"        : {"logY" : True,  "Y" : {"title" : "Events", "min" : 0.2},  "X" : {"title" : "Top p_{T} [GeV]", "rebin" : 5, "min" :  0,    "max" : 2000}},
    #
    #"h_ht_?"            : {"logY" : True,  "Y" : {"title" : "Events", "min" : 2e-3}, "X" : {"title" : "H_{T} [GeV]",                "rebin" : 1, "min" :  0, "max" : 3000}},  
    #"h_met_?"           : {"logY" : True,  "Y" : {"title" : "Events", "min" : 2e-3}, "X" : {"title" : "MET [GeV]",                  "rebin" : 1, "min" :  0, "max" : 2000}},  
    #"h_dR_bjets_?"      : {"logY" : False, "Y" : {"title" : "Events", "min" : 2e-3}, "X" : {"title" : "#DeltaR_{bjets} [GeV]",      "rebin" : 1, "min" :  0, "max" :   10}},  
    #"h_dR_top1_top2_?"  : {"logY" : False, "Y" : {"title" : "Events", "min" : 2e-3}, "X" : {"title" : "#DeltaR_{top1,top2} [GeV]",  "rebin" : 1, "min" :  0, "max" :   10}},  
    #"h_dR_tops_bjets_?" : {"logY" : False, "Y" : {"title" : "Events", "min" : 2e-3}, "X" : {"title" : "#DeltaR_{tops,bjets} [GeV]", "rebin" : 1, "min" :  0, "max" :   10}},  
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
