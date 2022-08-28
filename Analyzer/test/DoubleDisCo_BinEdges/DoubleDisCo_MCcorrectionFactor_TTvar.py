import numpy as np

from collections               import defaultdict
from common_Regions            import *
from common_Plotter            import *
from common_Aggregator         import *
from common_HiggsCombineInputs import *

class MCcorrectionFactor_TTvar():

    def __init__(self, year, channel, model, mass, translator):

        self.year       = year
        self.channel    = channel
        self.sig        = model
        self.mass       = mass
        self.translator = translator

        self.ttVarsModel = ["TT_fsrUp",            
                            "TT_fsrDown",          
                            "TT_isrUp",            
                            "TT_isrDown",          
                            "TT_hdampUp",          
                            "TT_hdampDown",        
                            "TT_underlyingEvtUp",
                            "TT_underlyingEvtDown",
                            "TT_erdOn"
        ]

        self.ttVarsDetect = ["TT_JECup",            
                             "TT_JECdown",          
                             "TT_JERup",            
                             "TT_JERdown"          
        ]

        self.ttVars = self.ttVarsModel + self.ttVarsDetect

        # make colors for each TT variance
        self.colors = { "TT"                   : "#525252",
                        "None"                 : "#525252",
                        "TT_erdOn"             : "#A6CEE3",
                        "TT_fsrUp"             : "#B2DF8A",
                        "TT_fsrDown"           : "#B2DF8A",
                        "TT_isrUp"             : "#FB9A99",
                        "TT_isrDown"           : "#FB9A99",
                        "TT_hdampUp"           : "#CAB2D6",
                        "TT_hdampDown"         : "#CAB2D6",
                        "TT_JECup"             : "#1F78B4",
                        "TT_JECdown"           : "#1F78B4",
                        "TT_JERup"             : "#33A02C",
                        "TT_JERdown"           : "#33A02C",
                        "TT_underlyingEvtUp"   : "#FDBF6F",
                        "TT_underlyingEvtDown" : "#FDBF6F",
        }

        self.correctionLabels = { "TT"                   : "Nominal",
                                  "TT_erdOn"             : "Color Reconn.",
                                  "TT_fsrUp"             : "FSR Up",
                                  "TT_fsrDown"           : "FSR Down",
                                  "TT_isrUp"             : "ISR Up",
                                  "TT_isrDown"           : "ISR Down",
                                  "TT_hdampUp"           : "ME-PS Up",
                                  "TT_hdampDown"         : "ME-PS Down",
                                  "TT_JECup"             : "JEC Up",
                                  "TT_JECdown"           : "JEC Down",
                                  "TT_JERup"             : "JER Up",
                                  "TT_JERdown"           : "JER Down",
                                  "TT_underlyingEvtUp"   : "Tune Up",
                                  "TT_underlyingEvtDown" : "Tune Down",
        }

        self.closureLabels = { "TT"                   : "w/ Nominal Corr.",
                               "None"                 : "w/o Corr.",
                               "TT_erdOn"             : "w/ Color Reconn. Corr.",
                               "TT_fsrUp"             : "w/ FSR Up Corr.",
                               "TT_fsrDown"           : "w/ FSR Down Corr.",
                               "TT_isrUp"             : "w/ ISR Up Corr.",
                               "TT_isrDown"           : "w/ ISR Down Corr.",
                               "TT_hdampUp"           : "w/ ME-PS Up Corr.",
                               "TT_hdampDown"         : "w/ ME-PS Down Corr.",
                               "TT_JECup"             : "w/ JEC Up Corr.",
                               "TT_JECdown"           : "w/ JEC Down Corr.",
                               "TT_JERup"             : "w/ JER Up Corr.",
                               "TT_JERdown"           : "w/ JER Down Corr.",
                               "TT_underlyingEvtUp"   : "w/ Tune Up Corr.",
                               "TT_underlyingEvtDown" : "w/ Tune Down Corr.",
        }

        self.valColors = { "Val_BD" : "#DDBB87",
                           "Val_CD" : "#429c93", 
                           "Val_D"  : "#990099",
        }

        self.regionGridWidth = 0.05

        # ------------------------------------------------
        # make the lists for rightBoundary and topBoundary
        # ------------------------------------------------
        self.list_boundaries = {"Val_BD" : np.arange(0.40, 1.05, self.regionGridWidth),
                                "Val_CD" : np.arange(0.40, 1.05, self.regionGridWidth),
                                "Val_D"  : np.arange(0.60, 1.05, self.regionGridWidth)
        }


    def run(self, disc1edge=None, disc2edge=None, fastMode=False, **kwargs):

        tablesPath     = kwargs["tablesPath"]["TT"]
        plotter        = kwargs["plotter"]
        regions        = kwargs["regions"]
        njets          = kwargs["njets"]
        samples        = kwargs["samples"] + self.ttVars
        histName       = kwargs["histName"]
        files          = kwargs["files"]

        # Will hold the pair of final edge values for the full ABCD region
        abcdFinalEdges = None

        # -------------------------------------------------------
        # Make an Aggregator object to hold everything convenient
        # -------------------------------------------------------
        theAggy = Aggregator(samples, njets, regions, self.list_boundaries["Val_BD"])

        # ---------------
        # loop over njets
        # --------------- 
        for njet in njets: 

            hist_lists = {}

            for sample in samples:

                hist_lists[sample] = files[sample].Get(histName + njet)

            minEdge  = hist_lists["TT"].GetXaxis().GetBinLowEdge(1) 
            maxEdge  = hist_lists["TT"].GetXaxis().GetBinLowEdge(hist_lists["TT"].GetNbinsX()+1)
            binWidth = hist_lists["TT"].GetXaxis().GetBinWidth(1)

            for ttVar in self.ttVars:

                # ------------------------------------------------------------------
                # Loop through the regions and make the set of plots for each
                # Make sure ABCD goes first so that the val regions can use its info
                # ------------------------------------------------------------------
                for region in regions:

                    theEdgesClass = None

                    # -----------------
                    # make ABCD regions
                    # -----------------
                    if region == "ABCD":

                        theEdgesClass = All_Regions(hist_lists, Sig=self.sig, ttVar=ttVar, disc1Edge=disc1edge, disc2Edge=disc2edge, fastMode=fastMode)
                        theAggy.aggregate(theEdgesClass, region = region, njet = njet)

                        abcdFinalEdges = theEdgesClass.getFinal("edges", "TT")

                    # ----------------------------
                    # make Sub-Division BD - Val I
                    # ----------------------------
                    elif region == "Val_BD":

                        for r in self.list_boundaries[region]:
                       
                            disc1_edge = ((float(abcdFinalEdges[0]) - 0.2) / (1.0 - 0.4)) * (float(r) - 0.4) + 0.2  
                            theEdgesClass = All_Regions(hist_lists, Sig=self.sig, ttVar=ttVar, disc2Edge=abcdFinalEdges[1], rightBoundary=float(r), disc1Edge=float(disc1_edge), fastMode=fastMode)
                            theAggy.aggregate(theEdgesClass, region = region, njet = njet, boundary = r)

                    # -----------------------------
                    # make Sub-Division CD - Val II
                    # -----------------------------
                    elif region == "Val_CD":

                        for t in self.list_boundaries[region]:

                            disc2_edge = ((float(abcdFinalEdges[1]) - 0.2) / (1.0 - 0.4)) * (float(t) - 0.4) + 0.2
                            theEdgesClass = All_Regions(hist_lists, Sig=self.sig, ttVar=ttVar, disc1Edge=abcdFinalEdges[0], topBoundary=float(t), disc2Edge=float(disc2_edge), fastMode=fastMode)
                            theAggy.aggregate(theEdgesClass, region = region, njet = njet, boundary = t)

                    # -----------------------------
                    # make Sub-Division D - Val III
                    # -----------------------------
                    elif region == "Val_D":

                        for d in self.list_boundaries[region]:            

                            disc1 = (float(abcdFinalEdges[0]) - (float(abcdFinalEdges[0]) / 2.0)) / (1.0 - float(abcdFinalEdges[0])) * (float(d) - float(abcdFinalEdges[0])) + (float(abcdFinalEdges[0]) / 2.0)
                            disc2 = (float(abcdFinalEdges[1]) - (float(abcdFinalEdges[1]) / 2.0)) / (1.0 - float(abcdFinalEdges[1])) * (float(d) - float(abcdFinalEdges[1])) + (float(abcdFinalEdges[1]) / 2.0)
                            theEdgesClass = All_Regions(hist_lists, Sig=self.sig, ttVar=ttVar, rightBoundary=d, topBoundary=d, disc1Edge=float(disc1), disc2Edge=float(disc2), fastMode=fastMode)
                            theAggy.aggregate(theEdgesClass, region = region, njet = njet, boundary = d)

        # -----------------------------------
        # Plot bkg subtracted data-MC closure
        # -----------------------------------
        for region in regions:
            if "Val_" not in region:
                continue

            for b in self.list_boundaries[region]:
                eventsTT = {}; eventsData = {}; edgesPerNjets = {}
                for njet in njets:

                    nEventsTT   = {subregion : theAggy.get("nEvents%s"%(subregion), region = region, njet = njet, boundary = b, sample = "TT"  ) for subregion in ["A", "B", "C", "D"]}
                    nEventsData = {subregion : theAggy.get("nEvents%s"%(subregion), region = region, njet = njet, boundary = b, sample = "Data") for subregion in ["A", "B", "C", "D"]}

                    eventsTT.setdefault(njet,      {}).setdefault(region, nEventsTT)
                    eventsData.setdefault(njet,    {}).setdefault(region, nEventsData)
                    edgesPerNjets.setdefault(njet, {}).setdefault(region, theAggy.get("finalEdges", region = region, njet = njet, boundary = b                 ))

        # ----------------------------------------
        # Plot non-closure as function of boundary
        # ----------------------------------------
        # initilaize a dictionary to include all variables for sys to put in Higgs Combine
        TT_TTvar_sys = { "MCcorr_TT"        : {},
                         "MCcorr_TTvar"     : {}, 
                         "MCcorr_Ratio_MC"  : {},
                         #"MCcoorRatio_Data" : {},
        }

        # -------------------
        # loop over the njets
        # -------------------
        for njet in njets:

            # initialize a sub-dictionaries for sys
            for key, value in TT_TTvar_sys.items():
                if njet not in value:
                    value[njet] = {}        
                
            # ---------------------
            # loop over the regions
            # -------------------
            for region in regions:
                if "Val_" not in region:
                    continue

                # initialize the dictionaries
                closureCorrPerBoundaryTT                        = {}
                closurePerBoundaryTTinData                      = {}
                MC_TT_corrected_dataClosure_PerBoundaryTTinData = {}
                MCcorrRatio_MC_BoundaryTT                       = {}

                closureCorrPerBoundaryTT["TT"]                          = theAggy.getPerBoundary(variable = "closureCorr",              sample = "TT",       region = region, njet = njet)            
                MC_TT_corrected_dataClosure_PerBoundaryTTinData["TT"]   = theAggy.getPerBoundary(variable = "MC_corrected_dataClosure", sample = "TTinData", region = region, njet = njet)
                MC_TT_corrected_dataClosure_PerBoundaryTTinData["None"] = theAggy.getPerBoundary(variable = "Closure",                  sample = "TTinData", region = region, njet = njet)

                # -------------------
                # loop over the ttVar
                # -------------------
                for ttVar in self.ttVars:
                
                    closureCorrPerBoundaryTT[ttVar]                        = theAggy.getPerBoundary(variable = "closureCorr",                    sample = ttVar,                 region = region, njet = njet)  
                    MC_TT_corrected_dataClosure_PerBoundaryTTinData[ttVar] = theAggy.getPerBoundary(variable = "MC_ttVar_corrected_dataClosure", sample = "TTinData_%s"%(ttVar), region = region, njet = njet)
                    MCcorrRatio_MC_BoundaryTT[ttVar]                       = theAggy.getPerBoundary(variable = "MCcorrRatio_MC",                 sample = ttVar,                 region = region, njet = njet)

                    # fill the TT_TTvar_sys dictionary
                    TT_TTvar_sys["MCcorr_TT"][njet]["TT"]        = closureCorrPerBoundaryTT["TT"][1.00]
                    TT_TTvar_sys["MCcorr_TTvar"][njet][ttVar]    = closureCorrPerBoundaryTT[ttVar][1.00]
                    TT_TTvar_sys["MCcorr_Ratio_MC"][njet][ttVar] = MCcorrRatio_MC_BoundaryTT[ttVar][1.00]

                # set y axis for higher njets bins
                yMin = None; yMax = None
                if int(njet.replace("incl", "")) > 8:
                    yMin = 0.7; yMax = 1.7

                else:
                    yMin = 0.65; yMax = 1.5

                # -------------------------------
                # Make plots for all TT variances
                # -------------------------------
                # TT
                plotter["TT"].plot_VarVsBoundary_MCcorrectionFactor_TTvar(closureCorrPerBoundaryTT,          self.ttVars + ["TT"], self.correctionLabels, self.regionGridWidth/2.0, yMin, yMax, 1.0, region,  "Closure Correction [TT]",   "TT_ClosureCorrection_PerBoundary",   njet, self.colors, self.valColors[region])

                # TTinData = Data - NonTT
                plotter["Data"].plot_VarVsBoundary_MCcorrectionFactor_TTvar(MC_TT_corrected_dataClosure_PerBoundaryTTinData, self.ttVars + ["TT", "None"], self.closureLabels, self.regionGridWidth/2.0, yMin, yMax,  1.0, region, "MC Corrected Data", "TTinData_MC_corrected_dataClosure_PerBoundary", njet, self.colors, self.valColors[region])

                # ----------------------------------------
                # Make plots for all TT modeling variances
                # ----------------------------------------
                # TT
                plotter["TT"].plot_VarVsBoundary_MCcorrectionFactor_TTvar(closureCorrPerBoundaryTT,          self.ttVarsModel + ["TT"], self.correctionLabels, self.regionGridWidth/2.0, yMin, yMax, 1.0, region,  "Closure Correction [TT]",   "TT_ClosureCorrection_PerBoundary_ModelVars",   njet, self.colors, self.valColors[region])

                # TTinData = Data - NonTT
                plotter["Data"].plot_VarVsBoundary_MCcorrectionFactor_TTvar(MC_TT_corrected_dataClosure_PerBoundaryTTinData, self.ttVarsModel + ["TT", "None"], self.closureLabels, self.regionGridWidth/2.0, yMin, yMax,  1.0, region, "MC Corrected Data", "TTinData_MC_corrected_dataClosure_PerBoundary_ModelVars", njet, self.colors, self.valColors[region])

                # -------------------------------------
                # Make plots for all detector variances
                # -------------------------------------
                # TT
                plotter["TT"].plot_VarVsBoundary_MCcorrectionFactor_TTvar(closureCorrPerBoundaryTT,          self.ttVarsDetect + ["TT"], self.correctionLabels, self.regionGridWidth/2.0, yMin, yMax, 1.0, region, "Closure Correction [TT]",   "TT_ClosureCorrection_PerBoundary_DetectorVars",   njet, self.colors, self.valColors[region])

                # TTinData = Data - NonTT
                plotter["Data"].plot_VarVsBoundary_MCcorrectionFactor_TTvar(MC_TT_corrected_dataClosure_PerBoundaryTTinData, self.ttVarsDetect + ["TT", "None"], self.closureLabels, self.regionGridWidth/2.0, yMin, yMax,  1.0, region, "MC Corrected Data", "TTinData_MC_corrected_dataClosure_PerBoundary_DetectorVars", njet, self.colors, self.valColors[region])

        # -------------------------------------------------------------
        # put each variable of TT_TTvar_sys dictionary to the root file
        # -------------------------------------------------------------
        higgsCombine = HiggsCombineInputs(self.year, njets, self.channel, samples, regions)        
        higgsCombine.put_HiggsCombineInputs_toRootFiles(TT_TTvar_sys)
        higgsCombine.close_HiggsCombineInputs_RootFiles()
