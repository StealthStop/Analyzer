#! /bin/env/python

import ROOT

from collections import OrderedDict

selections    = ["0l_ABCD", "1l_ABCD", "2l_ABCD"]
combinedNjets = ["6thJet",  "7thJet",  "8thJet" ]
variables     = []
for selection in selections:
    for combinedNjet in combinedNjets:
        variables.append( selection.replace("?", combinedNjet) )


histograms = {

    "h_Jet@_cm_PtrHT_?"          : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ p_{T}/H_{T} [GeV]", "rebin" : 2, "min" :  0,    "max" : 1500}},
    #"h_Jet@_cm_Pt_?"             : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ p_{T} [GeV]",       "rebin" : 5, "min" :  0,    "max" : 1500}},
    "h_Jet@_cm_Eta_?"            : {"logY" : False, "orders" : list(xrange(1,8)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ #eta",              "rebin" : 1, "min" : -4,    "max" :    4}},
    "h_Jet@_cm_Phi_?"            : {"logY" : False, "orders" : list(xrange(1,8)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ #phi",              "rebin" : 1, "min" : -6,    "max" :    6}},
    "h_Jet@_cm_Mass_?"           : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ mass [GeV]",        "rebin" : 5, "min" :  0,    "max" :  150}},
    "h_Jet@_cm_Energy_?"         : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ Energy [GeV]",      "rebin" : 5, "min" :  0,    "max" : 1500}},
    "h_Jet@_cm_Flavb_?"          : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ DeepFlavor b",      "rebin" : 1, "min" :  0,    "max" :    1}},
    "h_Jet@_cm_Flavc_?"          : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ DeepFlavor c",      "rebin" : 1, "min" :  0,    "max" :    1}},
    "h_Jet@_cm_Flavuds_?"        : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ DeepFlavor uds",    "rebin" : 1, "min" :  0,    "max" :    1}},
    "h_Jet@_cm_Flavq_?"          : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ DeepFlavor q",      "rebin" : 1, "min" :  0,    "max" :    1}},
    "h_Jet@_cm_Flavg_?"          : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ DeepFlavor g",      "rebin" : 1, "min" :  0,    "max" :    1}},
    
    "h_FWM@_top6_?"              : {"logY" : False, "orders" : list(xrange(2,6)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Fox-Wolfram Moment @",    "rebin" : 1, "min" :  0,    "max" :    1}},
    "h_JMT_ev@_top6_?"           : {"logY" : False, "orders" : list(xrange(0,3)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet p-E Eigenvalue @",    "rebin" : 1, "min" :  0,    "max" :    1}},
   
    #"h_Stop@_PtrHT_cm_OldSeed_?" : {"logY" : True,  "orders" : list(xrange(1,3)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Stop @ p_{T} [GeV]",      "rebin" : 1, "min" :  0,    "max" :    1}}, 
    #"h_Stop@_Pt_cm_OldSeed_?"    : {"logY" : True,  "orders" : list(xrange(1,3)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Stop @ p_{T} [GeV]",      "rebin" : 1, "min" :  0,    "max" : 1500}},
    "h_Stop@_Mass_cm_OldSeed_?"  : {"logY" : False, "orders" : list(xrange(1,3)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Stop @ mass [GeV]",       "rebin" : 1, "min" :  0,    "max" : 1500}},
    "h_Stop@_Eta_cm_OldSeed_?"   : {"logY" : False, "orders" : list(xrange(1,3)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Stop @ #eta",             "rebin" : 1, "min" : -4,    "max" :    4}},
    "h_Stop@_Phi_cm_OldSeed_?"   : {"logY" : False, "orders" : list(xrange(1,3)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Stop @ #phi",             "rebin" : 1, "min" : -6,    "max" :    6}},
   
    #"h_Stop@_PtrHT_cm_TopSeed_?" : {"logY" : True,  "orders" : list(xrange(1,3)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Stop @ p_{T} [GeV]",      "rebin" : 1, "min" :  0,    "max" :    1}}, 
    #"h_Stop@_Pt_cm_TopSeed_?"    : {"logY" : True,  "orders" : list(xrange(1,3)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Stop @ p_{T} [GeV]",      "rebin" : 1, "min" :  0,    "max" : 1500}},
    "h_Stop@_Mass_cm_TopSeed_?"  : {"logY" : False, "orders" : list(xrange(1,3)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Stop @ mass [GeV]",       "rebin" : 1, "min" :  0,    "max" : 1500}},
    "h_Stop@_Eta_cm_TopSeed_?"   : {"logY" : False, "orders" : list(xrange(1,3)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Stop @ #eta",             "rebin" : 1, "min" : -4,    "max" :    4}},
    "h_Stop@_Phi_cm_TopSeed_?"   : {"logY" : False, "orders" : list(xrange(1,3)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Stop @ #phi",             "rebin" : 1, "min" : -6,    "max" :    6}},
   
    #"h_combined?ToLastJet_cm_PtrHT_?" : {"logY" : True,  "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ p_{T}/H_{T} [GeV]", "rebin" : 2, "min" :  0,    "max" : 1500}},
    #"h_combined?ToLastJet_cm_Pt_?"    : {"logY" : True,  "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ p_{T} [GeV]",       "rebin" : 5, "min" :  0,    "max" : 1500}},
    #"h_combined?ToLastJet_cm_Eta_?"   : {"logY" : False, "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ #eta",              "rebin" : 1, "min" : -4,    "max" :    4}},
    #"h_combined?ToLastJet_cm_Phi_?"   : {"logY" : False, "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ #phi",              "rebin" : 1, "min" : -6,    "max" :    6}},
    #"h_combined?ToLastJet_cm_Mass_?"  : {"logY" : True,  "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ mass [GeV]",        "rebin" : 5, "min" :  0,    "max" :  150}},
    #"h_combined?ToLastJet_cm_Energy_?": {"logY" : True,  "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ Energy [GeV]",      "rebin" : 5, "min" :  0,    "max" : 1500}},

 
}

data = {
    "Data_skim" : {"name" : "Data", "color" : ROOT.kBlack, "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 1.3}
}
 

backgrounds = {
    "TT_skim"       : {"name" : "t#bar{t} + jets", "color" : 40,  "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "QCD_skim"      : {"name" : "QCD multijet",    "color" : 30,  "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "TTX_skim"      : {"name" : "t#bar{t} + X",    "color" : 38,  "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "BG_OTHER_skim" : {"name" : "Other",           "color" : 41,  "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0}
}

signals = OrderedDict({
    "RPV_2t6j_mStop-350"        : {"name" : "RPV m_{ #tilde{t}} = 350 GeV",                 "color" : 2, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "RPV_2t6j_mStop-550"        : {"name" : "RPV m_{ #tilde{t}} = 550 GeV",                 "color" : 7, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "RPV_2t6j_mStop-850"        : {"name" : "RPV m_{ #tilde{t}} = 850 GeV",                 "color" : 4, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "RPV_2t6j_mStop-1050"       : {"name" : "RPV m_{ #tilde{t}} = 1050 GeV",                "color" : 5, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "StealthSYY_2t6j_mStop-900" : {"name" : "Stealth SY#bar{Y} m_{ #tilde{t}} = 900 GeV",   "color" : 6, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0}
})

