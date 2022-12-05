#! /bin/env/python

import ROOT

from collections import OrderedDict

selections = [
    ("_0l_ABCD",30),
    ("_1l_ABCD",40),
    ("_2l_ABCD",38),
    ("_0l_QCDCR_ABCD",2),
    ("_1l_QCDCR_ABCD",7),
]

histograms = {
     "Jet_cm_ptrHT_@?"        : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ p_{T}/H_{T} [GeV]",    "rebin" : 5, "min" :  0,  "max" : 1500}},
     "Jet_cm_pt_@?"           : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ p_{T} [GeV]",          "rebin" : 2, "min" :  0,  "max" : 1500}},
     "Jet_cm_eta_@?"          : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ #eta",                 "rebin" : 1, "min" : -4,  "max" :    4}},
     "Jet_cm_phi_@?"          : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ #phi",                 "rebin" : 1, "min" : -6,  "max" :    6}},
     "Jet_cm_m_@?"            : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ mass [GeV]",           "rebin" : 2, "min" :  0,  "max" :  150}},
     "Jet_cm_E_@?"            : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ Energy [GeV]",         "rebin" : 5, "min" :  0,  "max" : 1500}},
     "Jet_cm_flavb_@?"        : {"logY" : True,  "orders" : list(xrange(1,8)), "Y" : {"title" : "A.U.", "min" : 0.2}, "X" : {"title" : "Jet @ DeepFlavor b",         "rebin" : 5, "min" :  0,  "max" :    1}},
    "h_Njets?"                                : {"logY" : True,                                 "Y" : {"title" : "A.U.", "min" : 2e-3}, "X" : {"title" : "N_{jets}",              "rebin" : 1, "min" :  4, "max" :   18}},
}

samples = OrderedDict({
    "TT"       : {"name" : "t#bar{t} + jets", "color" : 10,  "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
})
