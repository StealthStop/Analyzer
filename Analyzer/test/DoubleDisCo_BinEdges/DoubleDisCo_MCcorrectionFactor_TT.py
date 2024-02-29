import numpy as np

from collections               import defaultdict
from common_Regions            import *
from common_Plotter            import *
from common_Aggregator         import *
from common_HiggsCombineInputs import *

class MCcorrectionFactor_TT():

    def __init__(self, year, channel, model, mass, ttVar, translator):
        
        self.year       = year
        self.channel    = channel
        self.sig        = model
        self.mass       = mass
        self.ttVar      = ttVar
        self.translator = translator

        # make colors for each validation regions
        self.valColors = { "Val_BD" : "#DDBB87",
                           "Val_CD" : "#429c93", 
                           "Val_D"  : "#990099",
        }

        self.regionGridWidth = 0.05

    # -------------
    # main function
    # ------------- 
    def run(self, disc1edge=None, disc2edge=None, fastMode=False, **kwargs):

        tablesPath        = kwargs["tablesPath"]["TT"]
        plotter           = kwargs["plotter"]
        regions           = kwargs["regions"]
        njets             = kwargs["njets"]
        samples           = kwargs["samples"]
        histName          = kwargs["histName"]
        files             = kwargs["files"]
        plotVars2D        = kwargs["plotVars2D"]
        plotVarVsBoundary = kwargs["plotVarVsBoundary"]

        # ------------------------------------------------
        # make the lists for rightBoundary and topBoundary
        # ------------------------------------------------
        self.list_boundaries = {"Val_BD" : np.arange(0.40, 1.05, self.regionGridWidth),
                                "Val_CD" : np.arange(0.40, 1.05, self.regionGridWidth),
                                "Val_D"  : np.arange(0.60, 1.05, self.regionGridWidth)
        }

        # Will hold the pair of final edge values for the full ABCD region
        abcdFinalEdges = None

        # -------------------------------------------------------
        # Make an Aggregator object to hold everything convenient
        # -------------------------------------------------------
        theAggy = Aggregator(samples, njets, regions, self.list_boundaries["Val_BD"])

        #QCDCRInfo = {}
        # -----------------
        # Setup QCDCR stuff
        # -----------------
        #for njet in njets:
        #    
        #    hist_lists = {}

        #    for sample in samples:

        #        # get the fsr/isr, jec/jer higtograms from TT root file
        #        ttvarStr = ""

        #        if sample == "TT_fsrDown":
        #            ttvarStr = "_fsrDown"

        #        elif sample == "TT_fsrUp":
        #            ttvarStr = "_fsrUp"

        #        elif sample == "TT_isrDown":
        #            ttvarStr = "_isrDown"

        #        elif sample == "TT_isrUp":
        #            ttvarStr = "_isrUp"

        #        elif sample == "TT_JECdown":
        #            ttvarStr = "_JECdown"

        #        elif sample == "TT_JECup":
        #            ttvarStr = "_JECup"

        #        elif sample == "TT_JERdown":
        #            ttvarStr = "_JERdown"

        #        elif sample == "TT_JERup":
        #            ttvarStr = "_JERup"

        #        hist_lists[sample] = files[sample].Get(histName.replace("${NJET}", njet) + ttvarStr)
        #        if ttvarStr == "":
        #            histNameQCDCR = histName.split("Njets")[0] + "QCDCR_Njets" + histName.split("Njets")[1]
        #            hist_lists[sample+"_QCDCR"] = files[sample].Get(histNameQCDCR.replace("${NJET}", njet))


        #    theEdgesClass = All_Regions(hist_lists, Sig=self.sig, ttVar="", disc1Edge=disc1edge, disc2Edge=disc2edge, fastMode=fastMode, justEvents=True)

        #    QCDCRInfo[njet] = theEdgesClass.getQCDCRValues()

        #    del theEdgesClass

        # ---------------
        # loop over njets
        # --------------- 
        for njet in njets: 

            hist_lists = {}

            for sample in samples:

                # get the fsr/isr, jec/jer higtograms from TT root file
                ttvarStr = ""

                if sample == "TT_fsrDown":
                    ttvarStr = "_fsrDown"

                elif sample == "TT_fsrUp":
                    ttvarStr = "_fsrUp"

                elif sample == "TT_isrDown":
                    ttvarStr = "_isrDown"

                elif sample == "TT_isrUp":
                    ttvarStr = "_isrUp"
            
                elif sample == "TT_JECdown":
                    ttvarStr = "_JECdown"

                elif sample == "TT_JECup":
                    ttvarStr = "_JECup"

                elif sample == "TT_JERdown":
                    ttvarStr = "_JERdown"
            
                elif sample == "TT_JERup":
                    ttvarStr = "_JERup"

                hist_lists[sample] = files[sample].Get(histName.replace("${NJET}", njet) + ttvarStr)
                #if ttvarStr == "":
                #    histNameQCDCR = histName.split("Njets")[0] + "QCDCR_Njets" + histName.split("Njets")[1]
                #    hist_lists[sample + "_QCDCR"] = files[sample].Get(histNameQCDCR.replace("${NJET}", njet))


            minEdge  = hist_lists["TT"].GetXaxis().GetBinLowEdge(1) 
            maxEdge  = hist_lists["TT"].GetXaxis().GetBinLowEdge(hist_lists["TT"].GetNbinsX()+1)
            binWidth = hist_lists["TT"].GetXaxis().GetBinWidth(1)

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

                    theEdgesClass = All_Regions(hist_lists, Sig=self.sig, ttVar=self.ttVar, disc1Edge=disc1edge, disc2Edge=disc2edge, fastMode=fastMode)#, QCDCRInfo=QCDCRInfo)
                    theAggy.aggregate(theEdgesClass, region = region, njet = njet)

                    abcdFinalEdges = theEdgesClass.getFinal("edges", "TT")

                # ----------------------------
                # make Sub-Division BD - Val I
                # ----------------------------
                elif region == "Val_BD":

                    for r in self.list_boundaries[region]:
                   
                        disc1_edge = ((float(abcdFinalEdges[0]) - 0.2) / (1.0 - 0.4)) * (float(r) - 0.4) + 0.2  
                        theEdgesClass = All_Regions(hist_lists, Sig=self.sig, ttVar=self.ttVar, disc2Edge=abcdFinalEdges[1], rightBoundary=float(r), disc1Edge=float(disc1_edge), fastMode=fastMode)#, QCDCRInfo=QCDCRInfo)
                        theAggy.aggregate(theEdgesClass, region = region, njet = njet, boundary = r)

                # -----------------------------
                # make Sub-Division CD - Val II
                # -----------------------------
                elif region == "Val_CD":

                    for t in self.list_boundaries[region]:

                        disc2_edge = ((float(abcdFinalEdges[1]) - 0.2) / (1.0 - 0.4)) * (float(t) - 0.4) + 0.2
                        theEdgesClass = All_Regions(hist_lists, Sig=self.sig, ttVar=self.ttVar, disc1Edge=abcdFinalEdges[0], topBoundary=float(t), disc2Edge=float(disc2_edge), fastMode=fastMode)#, QCDCRInfo=QCDCRInfo)
                        theAggy.aggregate(theEdgesClass, region = region, njet = njet, boundary = t)

                # -----------------------------
                # make Sub-Division D - Val III
                # -----------------------------
                elif region == "Val_D":

                    for d in self.list_boundaries[region]:            

                        disc1_edge = ((float(abcdFinalEdges[0]) - 0.3) / (1.0 - 0.6)) * (float(d) - 0.6) + 0.3
                        disc2_edge = ((float(abcdFinalEdges[1]) - 0.3) / (1.0 - 0.6)) * (float(d) - 0.6) + 0.3

                        theEdgesClass = All_Regions(hist_lists, Sig=self.sig, ttVar=self.ttVar, rightBoundary=float(d), topBoundary=float(d), disc1Edge=float(disc1_edge), disc2Edge=float(disc2_edge), fastMode=fastMode)#, QCDCRInfo=QCDCRInfo)
                        theAggy.aggregate(theEdgesClass, region = region, njet = njet, boundary = d)

            # ---------------------------
            # Plot variable vs Disc1Disc2
            # ---------------------------
            for region in regions:

                if "Val_" not in region:
                    continue

                for b in self.list_boundaries[region]:

                    boundary = b

                    kwgs = {"region" : region, "njet" : njet, "boundary" : boundary}
  
                    edges                 = np.array(theAggy.get("edges", **kwgs), dtype=float)
                    finalEdges            = theAggy.get("finalEdges",                        **kwgs)
                    closures              = theAggy.get("nonClosures",  sample = "TT",       **kwgs)
                    closureCorr_TT        = theAggy.get("closureCorrs", sample = "TT",       **kwgs)
                    closureCorr_TTinData  = theAggy.get("closureCorrs", sample = "TTinData", **kwgs)

                    if plotVars2D:
                        plotter["TT"].plot_Var_vsDisc1Disc2(closures[:,0],               edges, float(finalEdges[0]), float(finalEdges[1]), minEdge, maxEdge, binWidth, 10e-10, 20.0, 0.0,  0.5, njet, name=region+"_%.2f"%(boundary), variable="NonClosure"             )
                        plotter["TT"].plot_Var_vsDisc1Disc2(closures[:,1],               edges, float(finalEdges[0]), float(finalEdges[1]), minEdge, maxEdge, binWidth, 10e-10, 20.0, 0.0,  0.5, njet, name=region+"_%.2f"%(boundary), variable="NonClosureUnc"          )
                        plotter["TT"].plot_Var_vsDisc1Disc2(closureCorr_TT[:,0],         edges, float(finalEdges[0]), float(finalEdges[1]), minEdge, maxEdge, binWidth, 10e-10, 20.0, 0.7,  1.3, njet, name=region+"_%.2f"%(boundary), variable="closureCorr_TT"         )
                        plotter["TT"].plot_Var_vsDisc1Disc2(closureCorr_TT[:,1],         edges, float(finalEdges[0]), float(finalEdges[1]), minEdge, maxEdge, binWidth, 10e-10, 20.0, 0.7,  1.3, njet, name=region+"_%.2f"%(boundary), variable="closureCorrUnc_TT"      )
                        plotter["Data"].plot_Var_vsDisc1Disc2(closureCorr_TTinData[:,0], edges, float(finalEdges[0]), float(finalEdges[1]), minEdge, maxEdge, binWidth, 10e-10, 20.0, 0.7,  1.3, njet, name=region+"_%.2f"%(boundary), variable="closureCorr_TTinData"   )
                        plotter["Data"].plot_Var_vsDisc1Disc2(closureCorr_TTinData[:,1], edges, float(finalEdges[0]), float(finalEdges[1]), minEdge, maxEdge, binWidth, 10e-10, 20.0, 0.7,  1.3, njet, name=region+"_%.2f"%(boundary), variable="closureCorrUnc_TTinData")
 
        # -----------------------------------
        # Plot bkg subtracted data-MC closure
        # -----------------------------------
        for region in regions:
            if "Val_" not in region:
                continue

            for b in self.list_boundaries[region]:
                boundary = b

                eventsTT = {}; eventsData = {}; edgesPerNjets = {}
                for njet in njets:

                    nEventsTT   = {subregion : theAggy.get("nEvents%s"%(subregion), region = region, njet = njet, boundary = boundary, sample = "TT"  ) for subregion in ["A", "B", "C", "D"]}
                    nEventsData = {subregion : theAggy.get("nEvents%s"%(subregion), region = region, njet = njet, boundary = boundary, sample = "Data") for subregion in ["A", "B", "C", "D"]}

                    eventsTT.setdefault(njet,      {}).setdefault(region, nEventsTT)
                    eventsData.setdefault(njet,    {}).setdefault(region, nEventsData)
                    edgesPerNjets.setdefault(njet, {}).setdefault(region, theAggy.get("finalEdges", region = region, njet = njet, boundary = boundary             ))

                plotter["Data"].make_allClosures(edgesPerNjets, eventsTT, None, None, eventsData, njets, name = region, closureTag = "b_%.2f"%(boundary), bkgTag = "forTT")

        # ----------------------------------------
        # Plot non-closure as function of boundary
        # ----------------------------------------
        for njet in njets:

            # TT
            weighted_nTTeventsA_PerBoundaryTT = {}; nonClosurePerBoundaryTT    = {}; 
            sigFractionA_PerBoundaryTT        = {}; sigFractionB_PerBoundaryTT = {}; sigFractionC_PerBoundaryTT = {}; sigFractionD_PerBoundaryTT = {}
            closureCorrPerBoundaryTT          = {}
            closureCorrPerBoundaryTTvar       = {}
            closurePerBoundaryTTvar           = {}

            # TTinData = Data - NonTT
            weighted_nTTEventsA_PerBoundaryTTinData = {}; nonClosurePerBoundaryTTinData = {}; closurePerBoundaryTTinData = {}; CorrectedDataClosure_PerBoundaryTTinData = {} 
            ttVar_CorrectedDataClosure_PerBoundaryTTinData = {} 

            for region in regions:
                if "Val_" not in region:
                    continue

                # TT
                weighted_nTTeventsA_PerBoundaryTT[region] = theAggy.getPerBoundary(variable = "nEventsA",    sample = "TT",       region = region, njet = njet)
                nonClosurePerBoundaryTT[region]           = theAggy.getPerBoundary(variable = "nonClosure",  sample = "TT",       region = region, njet = njet)
                sigFractionA_PerBoundaryTT[region]        = theAggy.getPerBoundary(variable = "sigFractionA",                     region = region, njet = njet)
                sigFractionB_PerBoundaryTT[region]        = theAggy.getPerBoundary(variable = "sigFractionB",                     region = region, njet = njet)
                sigFractionC_PerBoundaryTT[region]        = theAggy.getPerBoundary(variable = "sigFractionC",                     region = region, njet = njet)
                sigFractionD_PerBoundaryTT[region]        = theAggy.getPerBoundary(variable = "sigFractionD",                     region = region, njet = njet)
                closureCorrPerBoundaryTT[region]          = theAggy.getPerBoundary(variable = "closureCorr", sample = "TT",       region = region, njet = njet)            
                closureCorrPerBoundaryTTvar[region]       = theAggy.getPerBoundary(variable = "closureCorr", sample = self.ttVar, region = region, njet = njet)            
                closurePerBoundaryTTvar[region]           = theAggy.getPerBoundary(variable = "Closure",     sample = self.ttVar, region = region, njet = njet)

                closureCorrRatiosPerBoundary = {}
                closureCorrDiffsPerBoundary = {}
                for key_region in closureCorrPerBoundaryTT.keys():
                    closureCorrRatiosPerBoundary[key_region] = {}
                    closureCorrDiffsPerBoundary[key_region] = {}
                    for key_bin in closureCorrPerBoundaryTT[key_region].keys():
                        temp_val = closureCorrPerBoundaryTTvar[key_region][key_bin][0] / closureCorrPerBoundaryTT[key_region][key_bin][0]
                        temp_diff = closureCorrPerBoundaryTTvar[key_region][key_bin][0] - closureCorrPerBoundaryTT[key_region][key_bin][0]
                        temp_err = temp_val * math.sqrt((closureCorrPerBoundaryTTvar[key_region][key_bin][1]/closureCorrPerBoundaryTTvar[key_region][key_bin][0])**2 + (closureCorrPerBoundaryTT[key_region][key_bin][1]/closureCorrPerBoundaryTT[key_region][key_bin][0])**2)
                        closureCorrRatiosPerBoundary[key_region][key_bin] = (temp_val, temp_err)
                        closureCorrDiffsPerBoundary[key_region][key_bin] = (temp_diff, temp_err)

                correctedTTVarClosurePerBoundary = {}
                for key_region in closurePerBoundaryTTvar.keys():
                    correctedTTVarClosurePerBoundary[key_region] = {}
                    for key_bin in closurePerBoundaryTTvar[key_region].keys():
                        temp_val = closurePerBoundaryTTvar[key_region][key_bin][0] * closureCorrPerBoundaryTT[key_region][key_bin][0]
                        temp_err = temp_val * math.sqrt((closurePerBoundaryTTvar[key_region][key_bin][1] / closurePerBoundaryTTvar[key_region][key_bin][0])**2 + (closureCorrPerBoundaryTT[key_region][key_bin][1] / closureCorrPerBoundaryTT[key_region][key_bin][0])**2)
                        correctedTTVarClosurePerBoundary[key_region][key_bin] = (temp_val, temp_err)

                # TTinData = Data - NonTT
                weighted_nTTEventsA_PerBoundaryTTinData[region]        = theAggy.getPerBoundary(variable = "nEventsA",                   sample = "TTinData", region = region, njet = njet)
                nonClosurePerBoundaryTTinData[region]                  = theAggy.getPerBoundary(variable = "nonClosure",                 sample = "TTinData", region = region, njet = njet)
                closurePerBoundaryTTinData[region]                     = theAggy.getPerBoundary(variable = "Closure",                    sample = "TTinData", region = region, njet = njet)#, manualOverride=correctedTTVarClosurePerBoundary)
                CorrectedDataClosure_PerBoundaryTTinData[region]       = theAggy.getPerBoundary(variable = "CorrectedDataClosure",       sample = "TTinData", region = region, njet = njet)#, manualOverride=correctedTTVarClosurePerBoundary)
                ttVar_CorrectedDataClosure_PerBoundaryTTinData[region] = theAggy.getPerBoundary(variable = "ttVar_CorrectedDataClosure", sample = "TTinData", region = region, njet = njet)

            # set y axis for higher njets bins
            yMin = None; yMax = None
            if njet > 8:
                yMin = 0.5; yMax = 1.5

            else:
                yMin = 0.7; yMax = 1.3

            if plotVarVsBoundary:
                # TT
                #plotter["TT"].plot_VarVsBoundary(weighted_nTTeventsA_PerBoundaryTT, self.regionGridWidth/2.0, None, None, None, "Weighted Events in A [MC]", "TT_weighted_nTTeventsA_PerBoundary", njet, self.valColors)
                #plotter["TT"].plot_VarVsBoundary(nonClosurePerBoundaryTT,           self.regionGridWidth/2.0, 0.0,  0.3,  None, "Non-Closure [MC]",          "TT_NonClosure_PerBoundary",          njet, self.valColors)
                plotter["TT"].plot_VarVsBoundary(sigFractionA_PerBoundaryTT,        self.regionGridWidth/2.0,  0.0,  1.0,  1.0, "SigFrac\'A\'",              "TT_SigFracA_PerBoundary",            njet, self.valColors) 
                plotter["TT"].plot_VarVsBoundary(sigFractionB_PerBoundaryTT,        self.regionGridWidth/2.0, None, None, None, "SigFrac\'B\'",              "TT_SigFracB_PerBoundary",            njet, self.valColors)
                plotter["TT"].plot_VarVsBoundary(sigFractionC_PerBoundaryTT,        self.regionGridWidth/2.0, None, None, None, "SigFrac\'C\'",              "TT_SigFracC_PerBoundary",            njet, self.valColors)
                plotter["TT"].plot_VarVsBoundary(sigFractionD_PerBoundaryTT,        self.regionGridWidth/2.0, None, None, None, "SigFrac\'D\'",              "TT_SigFracD_PerBoundary",            njet, self.valColors)
                plotter["TT"].plot_VarVsBoundary(closureCorrPerBoundaryTT,          self.regionGridWidth/2.0, yMin, yMax, 1.0,   "Closure Correction [TT]",  "TT_ClosureCorrection_PerBoundary",   njet, self.valColors)
                plotter[self.ttVar].plot_VarVsBoundary(closureCorrPerBoundaryTTvar, self.regionGridWidth/2.0, yMin, yMax, 1.0,   "Closure Correction [%s]"%(self.ttVar), "%s_ClosureCorrection_PerBoundary"%(self.ttVar),   njet, self.valColors)

                plotter[self.ttVar].plot_VarVsBoundary(correctedTTVarClosurePerBoundary, self.regionGridWidth/2.0, yMin, yMax, 1.0,   "Residual Non-Closure (%s)"%(self.ttVar), "%s_CorrectedDataClosure"%(self.ttVar),   njet, self.valColors)
                plotter[self.ttVar].plot_VarVsBoundary(closureCorrRatiosPerBoundary, self.regionGridWidth/2.0, yMin, yMax, 1.0,   "Closure Correction %s / Closure Correciton Nominal"%(self.ttVar), "%s_ClosureCorrection_Ratio"%(self.ttVar),   njet, self.valColors)
                plotter[self.ttVar].plot_VarVsBoundary(closureCorrDiffsPerBoundary, self.regionGridWidth/2.0, None, None, None,   "Closure Correction %s - Closure Correciton Nominal"%(self.ttVar), "%s_ClosureCorrection_Diff"%(self.ttVar),   njet, self.valColors)
                    
                # TTinData = Data - NonTT
                #plotter["Data"].plot_VarVsBoundary(weighted_nTTEventsA_PerBoundaryTTinData,            self.regionGridWidth/2.0, None, None, None, "Weighted Events in A [Data]",  "TTinData_weighted_nTTeventsA_PerBoundary",     njet, self.valColors)
                #plotter["Data"].plot_VarVsBoundary(nonClosurePerBoundaryTTinData,                      self.regionGridWidth/2.0, 0.0, 0.3,   None, "Non-Closure [Data]",           "TTinData_NonClosure_PerBoundary",              njet, self.valColors)
                plotter["Data"].plot_VarVsBoundary(closurePerBoundaryTTinData,                     self.regionGridWidth/2.0, yMin, yMax,  1.0, "Data Closure",           "TTinData_DataClosure_PerBoundary",                     njet, self.valColors)
                plotter["Data"].plot_VarVsBoundary(CorrectedDataClosure_PerBoundaryTTinData,       self.regionGridWidth/2.0, yMin, yMax,  1.0, "Residual Non-Closure", "TTinData_CorrectedDataClosure_PerBoundary",            njet, self.valColors)
                #plotter["Data"].plot_VarVsBoundary(ttVar_CorrectedDataClosure_PerBoundaryTTinData, self.regionGridWidth/2.0, yMin, yMax,  1.0, "Residual Non-Closure", "TTinData_MC_%s_corrected_dataClosure_PerBoundary"%(self.ttVar), njet, self.valColors)
