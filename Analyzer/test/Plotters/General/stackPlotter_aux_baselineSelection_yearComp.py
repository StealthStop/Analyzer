#! /bin/env/python

import ROOT

from collections import OrderedDict


selections = ["0l_ABCD", "1l_ABCD", "2l_ABCD", "0l_QCDCR_ABCD", "1l_QCDCR_ABCD"] 

histograms = {

    # jets
    "h_Njets_?"    : {"logY" : True,  "Y" : {"title" : "Number of Events", "min" : 0.2}, "X" : {"title" : "N_{Jets}",        "rebin" : 1, "min" : -0.5, "max" : 20.5}},

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

}

data = {
    "TT" : {"name" : "Run2UL t#bar{t} + jets", "color" : ROOT.kBlack, "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 1.3}
} 

backgrounds = {
    "TT" : {"name" : "${YEAR} t#bar{t} + jets", "color" : 40,  "lstyle" : 1, "mstyle" : 8, "lsize" : 0, "msize" : 0},
}

signals = OrderedDict({
})

