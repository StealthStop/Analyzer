#! /bin/env/python

import ROOT

from collections import OrderedDict

selections = ["0l_blind_ABCD", "1l_blind_ABCD", "2l_blind_ABCD", "0l_ABCD", "1l_ABCD", "2l_ABCD", "0l_QCDCR_ABCD", "1l_QCDCR_ABCD"]
comboJets  = ["8",             "7",             "6",             "8",       "7",       "6",       "8",             "7"            ]

histograms = {
    "h_Jet@_PtrHT_cm_?"                  : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Jet @ p_{T}/H_{T}",                   "rebin" : 4, "min" :  0,    "max" :    1}},
    "h_Jet@_Pt_cm_?"                     : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Jet @ p_{T} [GeV]",                   "rebin" : 8, "min" :  0,    "max" : 1500}},
    "h_Jet@_Eta_cm_?"                    : {"logY" : False, "orders" : list(xrange(1,8)), "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Jet @ #eta",                          "rebin" : 1, "min" : -4,    "max" :    4}},
    "h_Jet@_Phi_cm_?"                    : {"logY" : False, "orders" : list(xrange(1,8)), "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Jet @ #phi",                          "rebin" : 1, "min" : -6,    "max" :    6}},
    "h_Jet@_Mass_cm_?"                   : {"logY" : False, "orders" : list(xrange(1,8)), "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Jet @ mass [GeV]",                    "rebin" : 2, "min" :  0,    "max" :  150}},
    "h_Jet@_Energy_cm_?"                 : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Jet @ Energy [GeV]",                  "rebin" : 8, "min" :  0,    "max" : 1500}},
    "h_Jet@_Flavb_cm_?"                  : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Jet @ DeepFlavor b",                  "rebin" : 1, "min" :  0,    "max" :    1}},
    "h_Jet@_Flavc_cm_?"                  : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Jet @ DeepFlavor c",                  "rebin" : 1, "min" :  0,    "max" :    1}},
    "h_Jet@_Flavuds_cm_?"                : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Jet @ DeepFlavor uds",                "rebin" : 1, "min" :  0,    "max" :    1}},
    "h_Jet@_Flavg_cm_?"                  : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Jet @ DeepFlavor g",                  "rebin" : 1, "min" :  0,    "max" :    1}},
    
    "h_FWM@_top6_?"                      : {"logY" : False, "orders" : list(xrange(2,6)), "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Fox-Wolfram Moment @",                "rebin" : 1, "min" :  0,    "max" :    1}},
    "h_JMT_ev@_top6_?"                   : {"logY" : False, "orders" : list(xrange(0,3)), "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Jet p-E Eigenvalue @",                "rebin" : 1, "min" :  0,    "max" :    1}},
   
    "h_Stop@_PtrHT_cm_OldSeed_?"         : {"logY" : True,  "orders" : list(xrange(1,3)), "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Stop @ p_{T}/H_{T}",                  "rebin" : 4, "min" :  0,    "max" :    1}}, 
    "h_Stop@_Pt_cm_OldSeed_?"            : {"logY" : True,  "orders" : list(xrange(1,3)), "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Stop @ p_{T} [GeV]",                  "rebin" : 8, "min" :  0,    "max" : 1500}},
    "h_Stop@_Mass_cm_OldSeed_?"          : {"logY" : False, "orders" : list(xrange(1,3)), "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Stop @ mass [GeV]",                   "rebin" : 2, "min" :  0,    "max" : 1500}},
    "h_Stop@_Eta_cm_OldSeed_?"           : {"logY" : False, "orders" : list(xrange(1,3)), "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Stop @ #eta",                         "rebin" : 1, "min" : -4,    "max" :    4}},
    "h_Stop@_Phi_cm_OldSeed_?"           : {"logY" : False, "orders" : list(xrange(1,3)), "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Stop @ #phi",                         "rebin" : 1, "min" : -6,    "max" :    6}},
   
    "h_Stop@_PtrHT_cm_TopSeed_?"         : {"logY" : True,  "orders" : list(xrange(1,3)), "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Stop @ p_{T}/H_{T}",                  "rebin" : 4, "min" :  0,    "max" :    1}}, 
    "h_Stop@_Pt_cm_TopSeed_?"            : {"logY" : True,  "orders" : list(xrange(1,3)), "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Stop @ p_{T} [GeV]",                  "rebin" : 8, "min" :  0,    "max" : 1500}},
    "h_Stop@_Mass_cm_TopSeed_?"          : {"logY" : False, "orders" : list(xrange(1,3)), "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Stop @ mass [GeV]",                   "rebin" : 2, "min" :  0,    "max" : 1500}},
    "h_Stop@_Eta_cm_TopSeed_?"           : {"logY" : False, "orders" : list(xrange(1,3)), "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Stop @ #eta",                         "rebin" : 1, "min" : -4,    "max" :    4}},
    "h_Stop@_Phi_cm_TopSeed_?"           : {"logY" : False, "orders" : list(xrange(1,3)), "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Stop @ #phi",                         "rebin" : 1, "min" : -6,    "max" :    6}},
}

for selection in selections:
    comboJet = comboJets[selections.index(selection)]

    histograms["h_combined%sthJet_PtrHT_cm_%s"%(comboJet, selection)]   = {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Combined %s to Last Jet p_{T}/H_{T}"%(comboJet),  "rebin" : 4, "min" :  0,    "max" :    1}},
    histograms["h_combined%sthJet_Pt_cm_%s"%(comboJet, selection)]      = {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Combined %s to Last Jet p_{T} [GeV]"%(comboJet),  "rebin" : 8, "min" :  0,    "max" : 1500}},
    histograms["h_combined%sthJet_Eta_cm_%s"%(comboJet, selection)]     = {"logY" : False, "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Combined %s to Last Jet #eta"%(comboJet),         "rebin" : 1, "min" : -4,    "max" :    4}},
    histograms["h_combined%sthJet_Phi_cm_%s"%(comboJet, selection)]     = {"logY" : False, "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Combined %s to Last Jet #phi"%(comboJet),         "rebin" : 1, "min" : -6,    "max" :    6}},
    histograms["h_combined%sthJet_Mass_cm_%s"%(comboJet, selection)]    = {"logY" : False, "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Combined %s to Last Jet mass [GeV]"%(comboJet),   "rebin" : 2, "min" :  0,    "max" :  150}},
    histograms["h_combined%sthJet_Energy_cm_%s"%(comboJet, selection)]  = {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Combined %s to Last Jet Energy [GeV]"%(comboJet), "rebin" : 8, "min" :  0,    "max" : 1500}},


data = {
    "Data" : {"name" : "Data", "color" : ROOT.kBlack, "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 1.3}
}
 

backgrounds = {
    "TT"       : {"name" : "t#bar{t} + jets", "color" : 40,  "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "QCD"      : {"name" : "QCD multijet",    "color" : 30,  "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "TTX"      : {"name" : "t#bar{t} + X",    "color" : 38,  "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "BG_OTHER" : {"name" : "Other",           "color" : 41,  "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0}
}

signals = OrderedDict({
    "RPV_2t6j_mStop-350"         : {"name" : "RPV m_{ #tilde{t}} = 350 GeV",                 "color" : 2, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "StealthSYY_2t6j_mStop-550"  : {"name" : "Stealth SY#bar{Y} m_{ #tilde{t}} = 550 GeV",   "color" : 7, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "StealthSYY_2t6j_mStop-850"  : {"name" : "Stealth SY#bar{Y} m_{ #tilde{t}} = 850 GeV",   "color" : 4, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "StealthSYY_2t6j_mStop-1250" : {"name" : "Stealth SY#bar{Y} m_{ #tilde{t}} = 1250 GeV",  "color" : 6, "lstyle" : 2, "mstyle" : 8, "lsize" : 3, "msize" : 0}
})

