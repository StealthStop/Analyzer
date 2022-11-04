import numpy as np

from collections               import defaultdict
from common_Regions            import *
from common_Plotter            import *
from common_Aggregator         import *
from common_TableWriter        import *
from common_HiggsCombineInputs import *


class MCcorrectionFactor_TTvar():

    def __init__(self, year, channel, model, mass, translator, edges):

        self.year        = year
        self.channel     = channel
        self.sig         = model
        self.mass        = mass
        self.translator  = translator
        self.edges       = edges

        self.ttVarsModel = [
                            "TT_fsrUp",            
                            "TT_fsrDown",          
                            "TT_isrUp",            
                            "TT_isrDown",          
                            "TT_hdampUP",          
                            "TT_hdampDOWN",        
                            "TT_TuneCP5up",
                            "TT_TuneCP5down",
                            "TT_erdON"
        ]

        self.ttVarsDetect = [
                             "TT_JECup",            
                             "TT_JECdown",          
                             "TT_JERup",            
                             "TT_JERdown"          
        ]

        self.ttVars = self.ttVarsModel + self.ttVarsDetect

        # make colors for each TT variance
        self.colors = { "TT"             : "#525252",
                        "None"           : "#525252",
                        "TT_erdON"       : "#A6CEE3",
                        "TT_fsrUp"       : "#B2DF8A",
                        "TT_fsrDown"     : "#B2DF8A",
                        "TT_isrUp"       : "#FB9A99",
                        "TT_isrDown"     : "#FB9A99",
                        "TT_hdampUP"     : "#CAB2D6",
                        "TT_hdampDOWN"   : "#CAB2D6",
                        "TT_JECup"       : "#1F78B4",
                        "TT_JECdown"     : "#1F78B4",
                        "TT_JERup"       : "#33A02C",
                        "TT_JERdown"     : "#33A02C",
                        "TT_TuneCP5up"   : "#FDBF6F",
                        "TT_TuneCP5down" : "#FDBF6F",
        }

        self.correctionLabels = { "TT"             : "Nominal",
                                  "TT_erdON"       : "Color Reconn.",
                                  "TT_fsrUp"       : "FSR Up",
                                  "TT_fsrDown"     : "FSR Down",
                                  "TT_isrUp"       : "ISR Up",
                                  "TT_isrDown"     : "ISR Down",
                                  "TT_hdampUP"     : "ME-PS Up",
                                  "TT_hdampDOWN"   : "ME-PS Down",
                                  "TT_JECup"       : "JEC Up",
                                  "TT_JECdown"     : "JEC Down",
                                  "TT_JERup"       : "JER Up",
                                  "TT_JERdown"     : "JER Down",
                                  "TT_TuneCP5up"   : "Tune Up",
                                  "TT_TuneCP5down" : "Tune Down",
        }

        self.closureLabels = { "TT"             : "w/ Nominal Corr.",
                               "None"           : "w/o Corr.",
                               "TT_erdON"       : "w/ Color Reconn. Corr.",
                               "TT_fsrUp"       : "w/ FSR Up Corr.",
                               "TT_fsrDown"     : "w/ FSR Down Corr.",
                               "TT_isrUp"       : "w/ ISR Up Corr.",
                               "TT_isrDown"     : "w/ ISR Down Corr.",
                               "TT_hdampUP"     : "w/ ME-PS Up Corr.",
                               "TT_hdampDOWN"   : "w/ ME-PS Down Corr.",
                               "TT_JECup"       : "w/ JEC Up Corr.",
                               "TT_JECdown"     : "w/ JEC Down Corr.",
                               "TT_JERup"       : "w/ JER Up Corr.",
                               "TT_JERdown"     : "w/ JER Down Corr.",
                               "TT_TuneCP5up"   : "w/ Tune Up Corr.",
                               "TT_TuneCP5down" : "w/ Tune Down Corr.",
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
        # create tex file
        # ---------------
        maxCorrData_ttSyst = maximumCorrectedData_ttSyst(tablesPath, self.channel, self.year, "TT_TTvar_Syst", self.sig)

        # ---------------
        # loop over njets
        # --------------- 
        for njet in njets: 

            hist_lists = {}

            for sample in samples:

                # get the fsr/isr higtograms from TT root file
                ttvarStr = ""

                if sample == "TT_fsrDown":
                    ttvarStr = "_fsrDown"

                elif sample == "TT_fsrUp":
                    ttvarStr = "_fsrUp"

                elif sample == "TT_isrDown":
                    ttvarStr = "_isrDown"

                elif sample == "TT_isrUp":
                    ttvarStr = "_isrUp"

                hist_lists[sample] = files[sample].Get(histName.replace("${NJET}", njet) + ttvarStr)


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
                    edgesPerNjets.setdefault(njet, {}).setdefault(region, theAggy.get("finalEdges", region = region, njet = njet, boundary = b))

        # ----------------------------------------
        # Plot non-closure as function of boundary
        # ----------------------------------------
        # initilaize a dictionary to include all variables for sys to put in Higgs Combine
        TT_TTvar_sys = { "MCcorr_TT"       : {},
                         "MCcorr_TTvar"    : {}, 
                         "MCcorr_Ratio_MC" : {},
        }

        # initialitze a dictionary to get the average and maximum values of MC corrected data 
        # for all variances and all val regions to put Higgs combine
        averageCorrectedDataValue_forAllRegions  = {}
        maximum_CorrectedDataValue_forAllRegions = {}
        best_correctedDataValue                  = {}
        best_correctedDataValue_Diff             = {}
        calculatedSys                            = {}

        for region in regions:
            averageCorrectedDataValue_forAllRegions[region]  = {}
            maximum_CorrectedDataValue_forAllRegions[region] = {}
            best_correctedDataValue[region]                  = {}
            calculatedSys[region]                            = {}

        averageCorrectedDataValue_forAllRegions["All"]  = {}
        maximum_CorrectedDataValue_forAllRegions["All"] = {}
        best_correctedDataValue["All"]                  = {}
        best_correctedDataValue_Diff["All"]             = {}
        calculatedSys["All"]                            = {}

        # -------------------
        # loop over the njets
        # -------------------
        for njet in njets:

            # initialize a sub-dictionaries for sys
            for key, value in TT_TTvar_sys.items():
                if njet not in value:
                    value[njet] = {}        
            
            # initialize variables to get the average and maximum values of MC corrected data 
            # for all variances and all val regions to put Higgs combine
            summedCorrectedDataValues_forAllRegions       = {}
            summedCorrectedDataValues_forAllRegions_count = {}         
            best_correctedDataValue                       = {}
            best_correctedDataValue_Diff                  = {}

            for region in regions:
                summedCorrectedDataValues_forAllRegions[region]       = 0.0
                summedCorrectedDataValues_forAllRegions_count[region] = 0.0     

            summedCorrectedDataValues_forAllRegions["All"]       = 0.0
            summedCorrectedDataValues_forAllRegions_count["All"] = 0.0
            best_correctedDataValue["All"]                       = None      
            best_correctedDataValue_Diff["All"]                  = -999.0  

            for region in regions:
                if "Val_" not in region:
                    continue

                MC_TT_corrected_dataClosure_PerBoundaryTTinData = {}
                MC_TT_corrected_dataClosure_PerBoundaryTTinData["TT"] = theAggy.getPerBoundary(variable = "MC_corrected_dataClosure", sample = "TTinData", region = region, njet = njet)

                for ttVar in self.ttVars:
                    MC_TT_corrected_dataClosure_PerBoundaryTTinData[ttVar] = theAggy.getPerBoundary(variable = "MC_ttVar_corrected_dataClosure", sample = "TTinData_%s"%(ttVar), region = region, njet = njet)

                # ----------------------------------------------------------------------------
                # fill a dictionary to get the average and maximum values of MC corrected data 
                # for all variances and all val regions to put Higgs combine
                # ----------------------------------------------------------------------------
                for ttProcess, boundaryDictionary in MC_TT_corrected_dataClosure_PerBoundaryTTinData.items():

                    if ("TT_" in ttProcess) or ("None" in ttProcess): continue

                    for boundaryValue, MC_corrected_dataClosure in boundaryDictionary.items(): 

                        # get the average values of MC corrected data
                        if MC_corrected_dataClosure[0] != -999.0:

                            summedCorrectedDataValues_forAllRegions[region]       += (1.0 / MC_corrected_dataClosure[0])
                            summedCorrectedDataValues_forAllRegions_count[region] += 1.0

                            summedCorrectedDataValues_forAllRegions["All"]       += (1.0 / MC_corrected_dataClosure[0])
                            summedCorrectedDataValues_forAllRegions_count["All"] += 1.0

                            # get the maximum values of MC corrected data
                            correctedDataValue_Diff = abs(MC_corrected_dataClosure[0] - 1.0)

                            if (correctedDataValue_Diff > best_correctedDataValue_Diff["All"]):

                                best_correctedDataValue_Diff["All"] = correctedDataValue_Diff
                                best_correctedDataValue["All"]      = MC_corrected_dataClosure[0]

                if summedCorrectedDataValues_forAllRegions_count[region] != 0.0:
                    averageCorrectedDataValue_forAllRegions[region][njet] = (summedCorrectedDataValues_forAllRegions[region] / summedCorrectedDataValues_forAllRegions_count[region])

                else:
                    averageCorrectedDataValue_forAllRegions[region][njet] = 1.0
           
            # -------------------------------------------------------------------------------
            # put the average MC corrected data value to the sub dictionary for all njet bins
            # ------------------------------------------------------------------------------- 
            if summedCorrectedDataValues_forAllRegions_count["All"] != 0.0:
                averageCorrectedDataValue_forAllRegions["All"][njet] = (summedCorrectedDataValues_forAllRegions["All"] / summedCorrectedDataValues_forAllRegions_count["All"])

            else:
                averageCorrectedDataValue_forAllRegions["All"][njet] = 1.0

            # -------------------------------------------------------------------------------
            # put the maximum MC corrected data value to the sub dictionary for all njet bins
            # ------------------------------------------------------------------------------- 
            maximum_CorrectedDataValue_forAllRegions["All"][njet] = best_correctedDataValue["All"]

            # ---------------------
            # loop over the regions
            # ---------------------
            for region in regions:
                if "Val_" not in region:
                    continue

                # initialize the dictionaries
                closureCorrPerBoundaryTT                        = {}
                closurePerBoundaryTTinData                      = {}
                MC_TT_corrected_dataClosure_PerBoundaryTTinData = {}
                MCcorrRatio_MC_BoundaryTT                       = {}
                MCcorrRatio_MC_Unc_BoundaryTT                   = {}

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
                    MCcorrRatio_MC_Unc_BoundaryTT[ttVar]                   = theAggy.getPerBoundary(variable = "MCcorrRatio_MC_Unc",             sample = ttVar,                 region = region, njet = njet)

                    # fill the TT_TTvar_sys dictionary
                    TT_TTvar_sys["MCcorr_TT"][njet]["TT"]        = closureCorrPerBoundaryTT["TT"][1.00]
                    TT_TTvar_sys["MCcorr_TTvar"][njet][ttVar]    = closureCorrPerBoundaryTT[ttVar][1.00]
                    TT_TTvar_sys["MCcorr_Ratio_MC"][njet][ttVar] = MCcorrRatio_MC_BoundaryTT[ttVar][1.00]

                # set y axis for higher njets bins
                yMin = None; yMax = None
                if int(njet.replace("incl", "")) > 8:
                    yMin = 0.7; yMax = 1.7

                else:
                    yMin = 0.7; yMax = 1.3

                # -------------------------------
                # Make plots for all TT variances
                # -------------------------------
                # TT
                plotter["TT"].plot_VarVsBoundary_MCcorrectionFactor_TTvar(closureCorrPerBoundaryTT,      self.ttVars + ["TT"], self.correctionLabels, self.regionGridWidth/2.0, yMin, yMax, 1.0, region,  "Closure Correction [TT]",                    "TT_ClosureCorrection_PerBoundary",   njet, self.colors, self.valColors[region])
                plotter["TT"].plot_VarVsBoundary_MCcorrectionFactor_TTvar(MCcorrRatio_MC_Unc_BoundaryTT, self.ttVars, self.correctionLabels, self.regionGridWidth/2.0, yMin, yMax, 1.0, region,  "Unc. on Closure Correction Ratio[TTvar/TT]", "MCcorrRatio_MC_Unc_PerBoundary",   njet, self.colors, self.valColors[region])

                # TTinData = Data - NonTT
                plotter["Data"].plot_VarVsBoundary_MCcorrectionFactor_TTvar(MC_TT_corrected_dataClosure_PerBoundaryTTinData, self.ttVars + ["TT", "None"], self.closureLabels, self.regionGridWidth/2.0, yMin, yMax,  1.0, region, "MC Corrected Data", "TTinData_MC_corrected_dataClosure_PerBoundary", njet, self.colors, self.valColors[region])

                # ----------------------------------------
                # Make plots for all TT modeling variances
                # ----------------------------------------
                # TT
                plotter["TT"].plot_VarVsBoundary_MCcorrectionFactor_TTvar(closureCorrPerBoundaryTT, self.ttVarsModel + ["TT"], self.correctionLabels, self.regionGridWidth/2.0, yMin, yMax, 1.0, region,  "Closure Correction [TT]",   "TT_ClosureCorrection_PerBoundary_ModelVars",   njet, self.colors, self.valColors[region])

                # TTinData = Data - NonTT
                plotter["Data"].plot_VarVsBoundary_MCcorrectionFactor_TTvar(MC_TT_corrected_dataClosure_PerBoundaryTTinData, self.ttVarsModel + ["TT", "None"], self.closureLabels, self.regionGridWidth/2.0, yMin, yMax,  1.0, region, "MC Corrected Data", "TTinData_MC_corrected_dataClosure_PerBoundary_ModelVars", njet, self.colors, self.valColors[region])

                # -------------------------------------
                # Make plots for all detector variances
                # -------------------------------------
                # TT
                plotter["TT"].plot_VarVsBoundary_MCcorrectionFactor_TTvar(closureCorrPerBoundaryTT, self.ttVarsDetect + ["TT"], self.correctionLabels, self.regionGridWidth/2.0, yMin, yMax, 1.0, region, "Closure Correction [TT]",   "TT_ClosureCorrection_PerBoundary_DetectorVars",   njet, self.colors, self.valColors[region])

                # TTinData = Data - NonTT
                plotter["Data"].plot_VarVsBoundary_MCcorrectionFactor_TTvar(MC_TT_corrected_dataClosure_PerBoundaryTTinData, self.ttVarsDetect + ["TT", "None"], self.closureLabels, self.regionGridWidth/2.0, yMin, yMax,  1.0, region, "MC Corrected Data", "TTinData_MC_corrected_dataClosure_PerBoundary_DetectorVars", njet, self.colors, self.valColors[region])


                # ----------------------------------------------
                # Make corrected MC corrected data closure plots
                # ----------------------------------------------
                #plotter["Data"].plot_VarVsBoundary_MCcorrectionFactor_TTvar(MC_TT_corrected_dataClosure_PerBoundaryTTinData, self.ttVars + ["TT", "None"], self.closureLabels, self.regionGridWidth/2.0, yMin, yMax,  1.0, region, "MC Corrected Data", "TTinData_MC_corrected_dataClosure_PerBoundary_ScaledbyTotalAve", njet, self.colors, self.valColors[region], scaleFactor=averageCorrectedDataValue_forAllRegions["All"][njet])
                #plotter["Data"].plot_VarVsBoundary_MCcorrectionFactor_TTvar(MC_TT_corrected_dataClosure_PerBoundaryTTinData, self.ttVarsModel + ["TT", "None"], self.closureLabels, self.regionGridWidth/2.0, yMin, yMax,  1.0, region, "MC Corrected Data", "TTinData_MC_corrected_dataClosure_PerBoundary_ScaledbyTotalAve_ModelVars", njet, self.colors, self.valColors[region], scaleFactor=averageCorrectedDataValue_forAllRegions["All"][njet])
                #plotter["Data"].plot_VarVsBoundary_MCcorrectionFactor_TTvar(MC_TT_corrected_dataClosure_PerBoundaryTTinData, self.ttVarsDetect + ["TT", "None"], self.closureLabels, self.regionGridWidth/2.0, yMin, yMax,  1.0, region, "MC Corrected Data", "TTinData_MC_corrected_dataClosure_PerBoundary_ScaledbyTotalAve_DetectorVars", njet, self.colors, self.valColors[region], scaleFactor=averageCorrectedDataValue_forAllRegions["All"][njet])

                #plotter["Data"].plot_VarVsBoundary_MCcorrectionFactor_TTvar(MC_TT_corrected_dataClosure_PerBoundaryTTinData, self.ttVars + ["TT", "None"], self.closureLabels, self.regionGridWidth/2.0, yMin, yMax,  1.0, region, "MC Corrected Data", "TTinData_MC_corrected_dataClosure_PerBoundary_ScaledbyRegionAve", njet, self.colors, self.valColors[region], scaleFactor=averageCorrectedDataValue_forAllRegions[region][njet])
                #plotter["Data"].plot_VarVsBoundary_MCcorrectionFactor_TTvar(MC_TT_corrected_dataClosure_PerBoundaryTTinData, self.ttVarsModel + ["TT", "None"], self.closureLabels, self.regionGridWidth/2.0, yMin, yMax,  1.0, region, "MC Corrected Data", "TTinData_MC_corrected_dataClosure_PerBoundary_ScaledbyRegionAve_ModelVars", njet, self.colors, self.valColors[region], scaleFactor=averageCorrectedDataValue_forAllRegions[region][njet])
                #plotter["Data"].plot_VarVsBoundary_MCcorrectionFactor_TTvar(MC_TT_corrected_dataClosure_PerBoundaryTTinData, self.ttVarsDetect + ["TT", "None"], self.closureLabels, self.regionGridWidth/2.0, yMin, yMax,  1.0, region, "MC Corrected Data", "TTinData_MC_corrected_dataClosure_PerBoundary_ScaledbyRegionAve_DetectorVars", njet, self.colors, self.valColors[region], scaleFactor=averageCorrectedDataValue_forAllRegions[region][njet])

        # -------------------------------------------------------------
        # calculate the sys. via the maximum value of MC corrected data
        # sys = (1.0 / maximum value of MC corrected data)
        # -------------------------------------------------------------
        for key_njet, maximum_correctedDataValue in maximum_CorrectedDataValue_forAllRegions["All"].items():

            if maximum_correctedDataValue != -999 and maximum_correctedDataValue != None:
                calculatedSys["All"][key_njet] =  (1.0 / maximum_CorrectedDataValue_forAllRegions["All"][key_njet])
            else:
                absoluteMax = -999.0
                someNjets = ["7", "8", "9"]
                for someNjet in someNjets:
                    if maximum_CorrectedDataValue_forAllRegions["All"][someNjet] == None: continue 
                    if abs(1- maximum_CorrectedDataValue_forAllRegions["All"][someNjet]) > absoluteMax:
                        absoluteMax = maximum_CorrectedDataValue_forAllRegions["All"][someNjet]

                calculatedSys["All"][key_njet] =  (1.0 / absoluteMax)
            
            # ------------------
            # write the tex file
            # ------------------
            #print "njet       : ", key_njet
            #print "MCcorr_TT  : ", TT_TTvar_sys["MCcorr_TT"][key_njet]["TT"][0] 
            #print "maxCorrData: ", maximum_CorrectedDataValue_forAllRegions["All"][key_njet]
            #print "ttSyst     : ", calculatedSys["All"][key_njet]

            if maximum_correctedDataValue != -999 and maximum_correctedDataValue != None:
                maxCorrData_ttSyst.writeLine(njet=key_njet, maxCorrData=maximum_CorrectedDataValue_forAllRegions["All"][key_njet], ttSyst=calculatedSys["All"][key_njet])

            else:
                maxCorrData_ttSyst.writeLine(njet=njet, maxCorrData=0.000, ttSyst=calculatedSys["All"][key_njet])


        # -------------------------------------------------------------
        # put each variable of TT_TTvar_sys dictionary to the root file
        # -------------------------------------------------------------
        higgsCombine = HiggsCombineInputs(self.year, njets, self.channel, samples, regions, self.edges)        
        higgsCombine.put_MCcorrFactors_toRootFiles(TT_TTvar_sys)
        
        for region in regions:

            if "Val" not in region: continue

            higgsCombine.put_averageCorrectedDataValue_toRootFiles(averageCorrectedDataValue_forAllRegions[region], region)
        
        higgsCombine.put_averageCorrectedDataValue_toRootFiles(averageCorrectedDataValue_forAllRegions["All"], "All")
        higgsCombine.put_maximumCorrectedDataValue_toRootFiles(calculatedSys["All"], "All")

        higgsCombine.close_HiggsCombineInputs_RootFiles()

        # ------------------
        # close the tex file
        # ------------------
        maxCorrData_ttSyst.writeClose()


