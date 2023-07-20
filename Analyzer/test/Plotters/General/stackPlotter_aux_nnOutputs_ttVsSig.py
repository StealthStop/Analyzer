#! /bin/env/python

import ROOT

model    = "SYY"
modelStr = "StealthSYY"

selections = ["_0l",
              "_1l",
              "_2l",
]

histograms = {
    # Regression Mass
    "h_DoubleDisCo_%s_massReg?_Njets@_ABCD"%(model)      : {"logY" : False,  "orders" : list(xrange(6,13)), "Y" : {"title" : "A.U.", "min" : 2e-3, "max" : 0.15}, "X" : {"title" : "Regression Mass [GeV]",  "rebin" : 2, "min" : 0,  "max" : 1500}},
    "h_DoubleDisCo_%s_massReg?_Njets@_ABCD"%(model)      : {"logY" : False,  "orders" : list(xrange(6,13)), "Y" : {"title" : "A.U.", "min" : 2e-3, "max" : 0.15}, "X" : {"title" : "Regression Mass [GeV]",  "rebin" : 2, "min" : 0,  "max" : 1500}},
    "h_DoubleDisCo_%s_massReg?_Njets@_ABCD"%(model)      : {"logY" : False,  "orders" : list(xrange(6,13)), "Y" : {"title" : "A.U.", "min" : 2e-3, "max" : 0.15}, "X" : {"title" : "Regression Mass [GeV]",  "rebin" : 2, "min" : 0,  "max" : 1500}},
    "h_DoubleDisCo_%s_massReg?_Njets11incl_ABCD"%(model) : {"logY" : False,                                 "Y" : {"title" : "A.U.", "min" : 2e-3, "max" : 0.15}, "X" : {"title" : "Regression Mass [GeV]",  "rebin" : 2, "min" : 0,  "max" : 1500}},
    "h_DoubleDisCo_%s_massReg?_Njets12incl_ABCD"%(model) : {"logY" : False,                                 "Y" : {"title" : "A.U.", "min" : 2e-3, "max" : 0.15}, "X" : {"title" : "Regression Mass [GeV]",  "rebin" : 2, "min" : 0,  "max" : 1500}},
    "h_DoubleDisCo_%s_massReg?_Njets13incl_ABCD"%(model) : {"logY" : False,                                 "Y" : {"title" : "A.U.", "min" : 2e-3, "max" : 0.15}, "X" : {"title" : "Regression Mass [GeV]",  "rebin" : 2, "min" : 0,  "max" : 1500}},
    # NN Discs.
    "h_DoubleDisCo_%s_disc1?_Njets@_ABCD"%(model)        : {"logY" : False,  "orders" : list(xrange(6,13)), "Y" : {"title" : "A.U.", "min" : 2e-3, "max" : 0.40}, "X" : {"title" : "Neural Network Disc. 1", "rebin" : 2, "min" : 0,  "max" : 1   }},
    "h_DoubleDisCo_%s_disc2?_Njets@_ABCD"%(model)        : {"logY" : False,  "orders" : list(xrange(6,13)), "Y" : {"title" : "A.U.", "min" : 2e-3, "max" : 0.40}, "X" : {"title" : "Neural Network Disc. 2", "rebin" : 2, "min" : 0,  "max" : 1   }},
    "h_DoubleDisCo_%s_disc1?_Njets11incl_ABCD"%(model)   : {"logY" : False,                                 "Y" : {"title" : "A.U.", "min" : 2e-3, "max" : 0.40}, "X" : {"title" : "Neural Network Disc. 1", "rebin" : 2, "min" : 0,  "max" : 1   }},
    "h_DoubleDisCo_%s_disc2?_Njets11incl_ABCD"%(model)   : {"logY" : False,                                 "Y" : {"title" : "A.U.", "min" : 2e-3, "max" : 0.40}, "X" : {"title" : "Neural Network Disc. 2", "rebin" : 2, "min" : 0,  "max" : 1   }},
    "h_DoubleDisCo_%s_disc1?_Njets12incl_ABCD"%(model)   : {"logY" : False,                                 "Y" : {"title" : "A.U.", "min" : 2e-3, "max" : 0.40}, "X" : {"title" : "Neural Network Disc. 1", "rebin" : 2, "min" : 0,  "max" : 1   }},
    "h_DoubleDisCo_%s_disc2?_Njets12incl_ABCD"%(model)   : {"logY" : False,                                 "Y" : {"title" : "A.U.", "min" : 2e-3, "max" : 0.40}, "X" : {"title" : "Neural Network Disc. 2", "rebin" : 2, "min" : 0,  "max" : 1   }},
    "h_DoubleDisCo_%s_disc1?_Njets13incl_ABCD"%(model)   : {"logY" : False,                                 "Y" : {"title" : "A.U.", "min" : 2e-3, "max" : 0.40}, "X" : {"title" : "Neural Network Disc. 1", "rebin" : 2, "min" : 0,  "max" : 1   }},
    "h_DoubleDisCo_%s_disc2?_Njets13incl_ABCD"%(model)   : {"logY" : False,                                 "Y" : {"title" : "A.U.", "min" : 2e-3, "max" : 0.40}, "X" : {"title" : "Neural Network Disc. 2", "rebin" : 2, "min" : 0,  "max" : 1   }},
}

backgrounds = {
    "TT" : {"name" : "t#bar{t} + jets", "color" : 40,  "lstyle" : 1, "mstyle" : 8, "lsize" : 0, "msize" : 0},
}

signals = {
    "%s_2t6j_mStop-350"%(modelStr)  : {"name" : "m_{ #tilde{t}} = 350 GeV",  "color" : 2,            "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "%s_2t6j_mStop-550"%(modelStr)  : {"name" : "m_{ #tilde{t}} = 550 GeV",  "color" : 4,            "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "%s_2t6j_mStop-850"%(modelStr)  : {"name" : "m_{ #tilde{t}} = 850 GeV",  "color" : 6,            "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    "%s_2t6j_mStop-1150"%(modelStr) : {"name" : "m_{ #tilde{t}} = 1150 GeV", "color" : ROOT.kCyan+1, "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
}

data = {
}
