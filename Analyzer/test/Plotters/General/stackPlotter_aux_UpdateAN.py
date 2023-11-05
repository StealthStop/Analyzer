#! /bin/env/python

import ROOT

from collections import OrderedDict


selections = ["0l_blind_ABCD", "1l_blind_ABCD", "2l_blind_ABCD", "0l_QCDCR_ABCD", "1l_QCDCR_ABCD", "2l_QCDCR_ABCD"] 
comboJets  = ["8",             "7",             "6",             "8",       "7",       "6",       "8",             "7"            ]

histograms = {

    # jets
    "h_Njets_?"    : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "N_{Jets}",        "rebin" : 1, "min" : -0.5, "max" : 20.5}},
    "h_njets_10incl_RPV_?"    : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "N_{Jets}",        "rebin" : 1, "min" : -0.5, "max" : 19.5}},
    "h_njets_11incl_RPV_?"    : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "N_{Jets}",        "rebin" : 1, "min" : -0.5, "max" : 19.5}},
    "h_njets_12incl_RPV_?"    : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "N_{Jets}",        "rebin" : 1, "min" : -0.5, "max" : 19.5}},
    "h_njets_10incl_SYY_?"    : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "N_{Jets}",        "rebin" : 1, "min" : -0.5, "max" : 20.5}},
    "h_njets_11incl_SYY_?"    : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "N_{Jets}",        "rebin" : 1, "min" : -0.5, "max" : 20.5}},
    "h_njets_12incl_SYY_?"    : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "N_{Jets}",        "rebin" : 1, "min" : -0.5, "max" : 20.5}},

    # bjets
    "h_Nbjets_?"    : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "N_{BJets}",       "rebin" : 1, "min" : -0.5, "max" : 20.5}},

    # tops
    "h_Ntops_?"     : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "N_{Tops}",                         "rebin" : 1,  "min" : -0.5,  "max" :  6.5}},
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

    # Jet
    "h_Jet@_PtrHT_cm_?"                  : {"logY" : False, "orders" : list(xrange(1,8)), "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Jet @ p_{T}/H_{T}",                   "rebin" : 4, "min" :  0,    "max" :    1}},
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
   
    "h_Stop@_PtrHT_cm_OldSeed_?"         : {"logY" : False, "orders" : list(xrange(1,3)), "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Stop @ p_{T}/H_{T}",                  "rebin" : 4, "min" :  0,    "max" :    1}}, 
    "h_Stop@_Pt_cm_OldSeed_?"            : {"logY" : True,  "orders" : list(xrange(1,3)), "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Stop @ p_{T} [GeV]",                  "rebin" : 8, "min" :  0,    "max" : 1500}},
    "h_Stop@_Mass_cm_OldSeed_?"          : {"logY" : False, "orders" : list(xrange(1,3)), "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Stop @ mass [GeV]",                   "rebin" : 2, "min" :  0,    "max" : 1500}},
    "h_Stop@_Eta_cm_OldSeed_?"           : {"logY" : False, "orders" : list(xrange(1,3)), "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Stop @ #eta",                         "rebin" : 1, "min" : -4,    "max" :    4}},
    "h_Stop@_Phi_cm_OldSeed_?"           : {"logY" : False, "orders" : list(xrange(1,3)), "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Stop @ #phi",                         "rebin" : 1, "min" : -6,    "max" :    6}},
   
    "h_Stop@_PtrHT_cm_TopSeed_?"         : {"logY" : False, "orders" : list(xrange(1,3)), "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Stop @ p_{T}/H_{T}",                  "rebin" : 4, "min" :  0,    "max" :    1}}, 
    "h_Stop@_Pt_cm_TopSeed_?"            : {"logY" : True,  "orders" : list(xrange(1,3)), "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Stop @ p_{T} [GeV]",                  "rebin" : 8, "min" :  0,    "max" : 1500}},
    "h_Stop@_Mass_cm_TopSeed_?"          : {"logY" : False, "orders" : list(xrange(1,3)), "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Stop @ mass [GeV]",                   "rebin" : 2, "min" :  0,    "max" : 1500}},
    "h_Stop@_Eta_cm_TopSeed_?"           : {"logY" : False, "orders" : list(xrange(1,3)), "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Stop @ #eta",                         "rebin" : 1, "min" : -4,    "max" :    4}},
    "h_Stop@_Phi_cm_TopSeed_?"           : {"logY" : False, "orders" : list(xrange(1,3)), "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Stop @ #phi",                         "rebin" : 1, "min" : -6,    "max" :    6}},
    
    "GoodLeptons_pt_@"      : {"logY" : False, "range" : [1,2],              "X" : {"title" : "Lepton @ p_{T} [GeV]",              "bins" : 80, "min" : 0, "max" : 360}},
    "GoodLeptons_phi_@"     : {"logY" : False, "range" : [1,2],              "X" : {"title" : "Lepton @ #phi",                     "bins" : 64, "min" : -4, "max" : 4}},
    "GoodLeptons_m_@"       : {"logY" : True,  "range" : [1,2],              "X" : {"title" : "Lepton @ Mass [GeV]",               "bins" : 72, "min" : 0, "max" : 0.12}},
    "GoodLeptons_eta_@"     : {"logY" : False, "range" : [1,2],              "X" : {"title" : "Lepton @ #eta",                     "bins" : 80, "min" : -6, "max" : 6}},
    "GoodLeptons_miniIso_@" : {"logY" : True, "range" : [1,2],              "X" : {"title" : "Lepton @ miniIso",                  "bins" : 80, "min" : 0, "max" : 0.3}},
    "Jet_pt_@"              : {"logY" : False, "range" : list(xrange(1,13)), "X" : {"title" : "Leading Jet p_{T} [GeV]",           "bins" : 72, "min" : 0, "max" : 1500}},
    "Jet_phi_@"             : {"logY" : False, "range" : list(xrange(1,13)), "X" : {"title" : "Jet #phi",                          "bins" : 64, "min" : -4, "max" : 4}},
    "Jet_m_@"               : {"logY" : False, "range" : list(xrange(1,13)), "X" : {"title" : "Jet Mass [GeV]",                    "bins" : 80, "min" : 0, "max" : 120}},
    "Jet_eta_@"             : {"logY" : False, "range" : list(xrange(1,13)), "X" : {"title" : "Jet #eta",                          "bins" : 80, "min" : -6, "max" : 6}},
    "Jet_ptD_@"             : {"logY" : False, "range" : list(xrange(1,13)), "X" : {"title" : "Jet p_{T,D}",                       "bins" : 80, "min" : 0, "max" : 1}}, 
    "Jet_axismajor_@"       : {"logY" : False, "range" : list(xrange(1,13)), "X" : {"title" : "Jet Axismajor",                     "bins" : 80, "min" : 0, "max" : 0.4}}, 
    "Jet_axisminor_@"       : {"logY" : False, "range" : list(xrange(1,13)), "X" : {"title" : "Jet Axisminor",                     "bins" : 80, "min" : 0, "max" : 0.4}}, 
    "Jet_flavg_@"           : {"logY" : False, "range" : list(xrange(1,13)), "X" : {"title" : "Jet DeepFlavour g",                 "bins" : 80, "min" : 0.0, "max" : 1.0}}, 
    "Jet_flavb_@"           : {"logY" : False, "range" : list(xrange(1,13)), "X" : {"title" : "Jet DeepFlavour b",                 "bins" : 80, "min" : 0.0, "max" : 1.0}}, 
    "Jet_flavc_@"           : {"logY" : False, "range" : list(xrange(1,13)), "X" : {"title" : "Jet DeepFlavour c",                 "bins" : 80, "min" : 0.0, "max" : 1.0}}, 
    "Jet_flavuds_@"         : {"logY" : False, "range" : list(xrange(1,13)), "X" : {"title" : "Jet DeepFlavour uds",               "bins" : 80, "min" : 0.0, "max" : 1.0}}, 
    "Jet_cEF_@"             : {"logY" : True,  "range" : list(xrange(1,13)), "X" : {"title" : "Jet Charged EM Fraction",           "bins" : 80, "min" : 0.0, "max" : 1.0}}, 
    "Jet_nEF_@"             : {"logY" : True,  "range" : list(xrange(1,13)), "X" : {"title" : "Jet Neutral EM Fraction",           "bins" : 80, "min" : 0.0, "max" : 1.0}}, 
    "Jet_cHF_@"             : {"logY" : True,  "range" : list(xrange(1,13)), "X" : {"title" : "Jet Charged Had. Fraction",         "bins" : 80, "min" : 0.0, "max" : 1.0}}, 
    "Jet_nHF_@"             : {"logY" : True,  "range" : list(xrange(1,13)), "X" : {"title" : "Jet Neutral Had. Fraction",         "bins" : 80, "min" : 0.0, "max" : 1.0}}, 
    "fwm@_top6"             : {"logY" : False, "range" : list(xrange(2,11)), "X" : {"title" : "Fox-Wolfram Moment @",              "bins" : 80, "min" : 0, "max" : 1}}, 
    "jmt_ev@_top6"          : {"logY" : False, "range" : list(xrange(0,3)),  "X" : {"title" : "Jet E-#vec{p} Tensor Eigenvalue @", "bins" : 80, "min" : 0, "max" : 1}},
    "lvMET_cm_pt"           : {"logY" : False, "range" : [-1],               "X" : {"title" : "|#vec{E}_{T}^{miss}| [GeV]",        "bins" : 80, "min" : 0, "max" : 500}},
    "lvMET_cm_phi"          : {"logY" : False, "range" : [-1],               "X" : {"title" : "#vec{E}_{T}^{miss} #phi",           "bins" : 64, "min" : -4, "max" : 4}},
    "lvMET_cm_m"            : {"logY" : False, "range" : [-1],               "X" : {"title" : "#vec{E}_{T}^{miss} Mass [GeV]",     "bins" : 80, "min" : 0, "max" : 1}},
    "lvMET_cm_eta"          : {"logY" : False, "range" : [-1],               "X" : {"title" : "#vec{E}_{T}^{miss} #eta",           "bins" : 80, "min" : -6, "max" : 6}},
    "Stop@_mass_cm_OldSeed" : {"logY" : False, "range" : [1,2], "X" : {"title" : "Stop@ Mass (p_{T}-Ranked, OldSeed) [GeV]",  "bins" : 72, "min" : 0, "max" : 1500}},
    "Stop@_eta_cm_OldSeed"  : {"logY" : False, "range" : [1,2], "X" : {"title" : "Stop@ #eta (p_{T}-Ranked, OldSeed)",        "bins" : 80, "min" : -6, "max" : 6}},
    "Stop@_phi_cm_OldSeed"  : {"logY" : False, "range" : [1,2], "X" : {"title" : "Stop@ #phi (p_{T}-Ranked, OldSeed)",        "bins" : 64, "min" : -4, "max" : 4}},
    "Stop@_pt_cm_OldSeed"   : {"logY" : False, "range" : [1,2], "X" : {"title" : "Stop@ p_{T} (p_{T}-Ranked, OldSeed) [GeV]", "bins" : 72, "min" : 0, "max" : 1500}},
}

for selection in selections:
    comboJet = comboJets[selections.index(selection)]

    histograms["h_combined%sthJet_PtrHT_cm_%s"%(comboJet, selection)]   = {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Combined %s to Last Jet p_{T}/H_{T}"%(comboJet),  "rebin" : 4, "min" :  0,    "max" :    1}}
    histograms["h_combined%sthJet_Pt_cm_%s"%(comboJet, selection)]      = {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Combined %s to Last Jet p_{T} [GeV]"%(comboJet),  "rebin" : 8, "min" :  0,    "max" : 1500}}
    histograms["h_combined%sthJet_Eta_cm_%s"%(comboJet, selection)]     = {"logY" : False, "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Combined %s to Last Jet #eta"%(comboJet),         "rebin" : 1, "min" : -4,    "max" :    4}}
    histograms["h_combined%sthJet_Phi_cm_%s"%(comboJet, selection)]     = {"logY" : False, "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Combined %s to Last Jet #phi"%(comboJet),         "rebin" : 1, "min" : -6,    "max" :    6}}
    histograms["h_combined%sthJet_Mass_cm_%s"%(comboJet, selection)]    = {"logY" : False, "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Combined %s to Last Jet mass [GeV]"%(comboJet),   "rebin" : 2, "min" :  0,    "max" :  150}}
    histograms["h_combined%sthJet_Energy_cm_%s"%(comboJet, selection)]  = {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "Combined %s to Last Jet Energy [GeV]"%(comboJet), "rebin" : 8, "min" :  0,    "max" : 1500}}

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
