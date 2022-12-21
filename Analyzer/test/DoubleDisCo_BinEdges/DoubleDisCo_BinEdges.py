from common_Regions            import *
from common_Plotter            import *
from common_TableWriter        import *

class BinEdges():

    def __init__(self, year, channel, model, mass, ttVar, translator):
        
        self.year       = year
        self.channel    = channel
        self.sig        = model
        self.mass       = mass
        self.ttVar      = ttVar
        self.translator = translator

    def run(self, disc1edge=None, disc2edge=None, fastMode=False, **kwargs):

        tablesPath     = kwargs["tablesPath"]["TT"]
        plotter        = kwargs["plotter"]
        regions        = kwargs["regions"]
        njets          = kwargs["njets"]
        samples        = kwargs["samples"]
        histName       = kwargs["histName"]
        files          = kwargs["files"]

        # ------------------
        # make all tex files
        # ------------------
        # get the signal fracs for each region
        sigFracsTable_AllRegions = SignalFractionsAllRegionsTable(tablesPath, self.channel, self.year, "Sig_fracs_BinEdges_AllRegions", self.sig, self.mass)

        # get the fracs for each ABCD and Validation region
        abcdFracsTable = ABCDfracsTable(tablesPath, self.channel, self.year, "Sig_TT_BinEdges_fixed_ABCD", self.sig)
        valFracsTable  = ValFracsTable(tablesPath,  self.channel, self.year, "Sig_BinEdges_fixed_Val",     self.sig)

        EventsPerNjets = {}
        for sample in samples:
            # Hold onto per Njet things so we can plot altogether after the Njets loop
            EventsPerNjets[sample] = {njet : None for njet in njets}

        # hold on edges per njet
        edgesPerNjets = {njet : None for njet in njets}

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
            
    
            minEdge  = hist_lists["TT"].GetXaxis().GetBinLowEdge(1) 
            maxEdge  = hist_lists["TT"].GetXaxis().GetBinLowEdge(hist_lists["TT"].GetNbinsX()+1)
            binWidth = hist_lists["TT"].GetXaxis().GetBinWidth(1)

            # initialize the dictionaries of quantities and variables with any combination of bin edges 
            # initialize the dictionaries of quantities and variables with the final choice of bin edges
            allRegionsFinalEdges  = {}
            allRegionsSigFracs_TT = {}; allRegionsFinalSigFracs_TT = {}
            allRegionsTTFracs     = {}; allRegionsFinalTTFracs     = {}
            allRegionsEvents      = {}; allRegionsFinalEvents      = {} 

            # loop over for initialize the big dictionaries to get all regions' events
            for hist_key in hist_lists.keys():

                allRegionsEvents[hist_key]      = {}
                allRegionsFinalEvents[hist_key] = {}

            # loop over for populating the dictionaries for all regions
            for region in regions:

                allRegionsSigFracs_TT[region]      = {"A" : {}, "B" : {}, "C" : {}, "D" : {}}  
                allRegionsFinalSigFracs_TT[region] = {"A" : {}, "B" : {}, "C" : {}, "D" : {}}
                allRegionsTTFracs[region]          = {"A" : {}, "B" : {}, "C" : {}, "D" : {}}
                allRegionsFinalTTFracs[region]     = {"A" : {}, "B" : {}, "C" : {}, "D" : {}}

                # loop over for TT, NonTT, Data to get all regions' events
                for hist_key in hist_lists.keys():

                    allRegionsEvents[hist_key][region]      = {"A" : {}, "B" : {}, "C" : {}, "D" : {}}
                    allRegionsFinalEvents[hist_key][region] = {"A" : {}, "B" : {}, "C" : {}, "D" : {}}

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

                    theEdgesClass = All_Regions(hist_lists, Sig=self.sig, ttVar=self.ttVar, disc1Edge=disc1edge, disc2Edge=disc2edge, fastMode=fastMode)

                # --------------------------
                #        *    ||   |                       
                #    vB  * vA ||   |    A                
                #  ______*____||___|________             
                #        *    ||   |                     
                #    vD  * vC ||   |    C                
                #        *                  
                # --------------------------             
                elif region == "Val_BD":

                   theEdgesClass = All_Regions(hist_lists, Sig=self.sig, ttVar=self.ttVar, disc2Edge=allRegionsFinalEdges["ABCD"][1], rightBoundary=float(0.4), disc1Edge=float(0.2), fastMode=fastMode)

                # -------------------
                #          B  |  A   
                #       ______|______
                #             |
                #       ------------- (0.4)
                #             |
                #         hB  |  hA 
                #       *************
                #             |
                #         hD  |  hC
                # --------------------
                elif region == "Val_CD":

                    theEdgesClass = All_Regions(hist_lists, Sig=self.sig, ttVar=self.ttVar, disc1Edge=allRegionsFinalEdges["ABCD"][0], topBoundary=float(0.4), disc2Edge=float(0.2), fastMode=fastMode)

                # ------------------------------------------------
                # make the validation regions as sub-division of D 
                #               |
                #        B      |   A
                #               |
                #     __________|_______
                #         *     |
                #      dB * dA  | 
                #     ********* |   C
                #      dD * dC  |
                #         *     |
                # ------------------------------------------------
                elif region == "Val_D":
                    
                    theEdgesClass = All_Regions(hist_lists, Sig=self.sig, ttVar=self.ttVar, rightBoundary=allRegionsFinalEdges["ABCD"][0], topBoundary=allRegionsFinalEdges["ABCD"][1], ABCDdisc1=float(allRegionsFinalEdges["ABCD"][0])/2.0, ABCDdisc2=float(allRegionsFinalEdges["ABCD"][1])/2.0, fastMode=fastMode)

                # ----------------------------
                # Optimization based on TT !!!
                # ----------------------------
                # get all final bin edges
                allRegionsFinalEdges[region] = theEdgesClass.getFinal("edges", "TT")

                # quantities and variables with any combination of bin edges
                significances                     = theEdgesClass.get("significance",                     None, None, "TT")
                significances_includingNonClosure = theEdgesClass.get("significance_includingNonClosure", None, None, "TT")
                significances_nonSimplified       = theEdgesClass.get("significance_nonSimplified",       None, None, "TT")

                allRegionsSigFracs_TT[region]     = {"A" : theEdgesClass.get("sigFractionA", None, None, self.sig),
                                                     "B" : theEdgesClass.get("sigFractionB", None, None, self.sig),
                                                     "C" : theEdgesClass.get("sigFractionC", None, None, self.sig),
                                                     "D" : theEdgesClass.get("sigFractionD", None, None, self.sig)}
                allRegionsTTFracs[region]         = {"A" : theEdgesClass.get("ttFractionA",  None, None, "TT"),
                                                     "B" : theEdgesClass.get("ttFractionB",  None, None, "TT"),
                                                     "C" : theEdgesClass.get("ttFractionC",  None, None, "TT"),
                                                     "D" : theEdgesClass.get("ttFractionD",  None, None, "TT")}

                # quantities and variables with the final choice of bin edges
                finalSignificance                  = theEdgesClass.getFinal("significance", "TT")
                allRegionsFinalSigFracs_TT[region] = {"A" : theEdgesClass.getFinal("sigFractionA", self.sig),
                                                      "B" : theEdgesClass.getFinal("sigFractionB", self.sig),
                                                      "C" : theEdgesClass.getFinal("sigFractionC", self.sig),
                                                      "D" : theEdgesClass.getFinal("sigFractionD", self.sig)}
                allRegionsFinalTTFracs[region]     = {"A" : theEdgesClass.getFinal("ttFractionA",  "TT"),
                                                      "B" : theEdgesClass.getFinal("ttFractionB",  "TT"),
                                                      "C" : theEdgesClass.getFinal("ttFractionC",  "TT"),
                                                      "D" : theEdgesClass.getFinal("ttFractionD",  "TT")}
                
                # -----------------------------------------------
                # loop over for getting plots for TT, NonTT, Data
                # -----------------------------------------------
                for hist_key in hist_lists.keys():

                    edges                                   = np.array(theEdgesClass.get("edges"), dtype=float)
                    allRegionsEvents[hist_key][region]      = {"A" : theEdgesClass.get("nEventsA", None, None, hist_key)}
                    allRegionsFinalEvents[hist_key][region] = {"A" : theEdgesClass.getFinal("nEventsA", hist_key),
                                                               "B" : theEdgesClass.getFinal("nEventsB", hist_key),
                                                               "C" : theEdgesClass.getFinal("nEventsC", hist_key),
                                                               "D" : theEdgesClass.getFinal("nEventsD", hist_key)}
                    if hist_key != self.sig:
                        nonClosures      = theEdgesClass.get("nonClosure",      None, None, hist_key) # vars with any combination of bin edges
                        pull             = theEdgesClass.get("pull",            None, None, hist_key)
                        finalNonClosure  = theEdgesClass.getFinal("nonClosure",             hist_key) # vars with the final choice of bin edges
                        finalPull        = theEdgesClass.getFinal("pull",                   hist_key)
                   
                    # ---------------------------  
                    # plot variable vs disc as 1D
                    # --------------------------- 
                    if kwargs["plotVars1D"]:
                        for disc in [1, 2]:

                            if hist_key == self.sig:
                                plotter["TT"].plot_VarVsDisc(allRegionsEvents[hist_key][region]["A"], edges, binWidth/2.0, -1.0, "Weighted Signal Events", "wSigEvts", disc, njet, name = region)

                            elif hist_key != self.ttVar:
                                plotter[hist_key].plot_VarVsDisc(allRegionsEvents[hist_key][region]["A"], edges, binWidth/2.0, -1.0, "Weighted Background Events", "wBkgEvts", disc, njet, name = region)
                                plotter[hist_key].plot_VarVsDisc(nonClosures,                             edges, binWidth/2.0, 1.0,  "Closure",                    "Closure",  disc, njet, name = region)

                            if hist_key == "TT":
                                plotter[hist_key].plot_VarVsDisc(significances, edges, binWidth/2.0, 5.0, "%s Significance"%(region), "Significance", disc, njet, name = region)

                                for subregion in self.translator["ABCD"].keys():
                                    plotter[hist_key].plot_VarVsDisc(allRegionsSigFracs_TT[region][subregion], edges, binWidth/2.0, 0.8, "Signal Contamination %s"%(self.translator[region][subregion]), "SigFracs%s"%(self.translator[region][subregion]), disc, njet, name = region)

                    # ---------------------------
                    # plot variable vs Disc1Disc2
                    # ---------------------------
                    if kwargs["plotVars2D"]:
                        absMin = 0.0
                        if hist_key == "TT":
                            # significance including non-closure
                            plotter[hist_key].plot_Var_vsDisc1Disc2(significances_includingNonClosure[:,0], edges, float(allRegionsFinalEdges[region][0]), float(allRegionsFinalEdges[region][1]), minEdge, maxEdge, binWidth, absMin, 20.0, 0.0,  5.0, njet, name=region, variable="Sign_includingNonClosure"        )
                            #plotter[hist_key].plot_Var_vsDisc1Disc2(significances_includingNonClosure[:,1], edges, float(allRegionsFinalEdges[region][0]), float(allRegionsFinalEdges[region][1]), minEdge, maxEdge, binWidth, absMin, 20.0, 0.0,  5.0, njet, name=region, variable="SignUnc_includingNonClosure"     )
                            # sigFracs A, B, C, D
                            if region == "ABCD":
                                #plotter[hist_key].plot_Var_vsDisc1Disc2(allRegionsSigFracs_TT[region]["A"][:,0],   edges, float(allRegionsFinalEdges[region][0]), float(allRegionsFinalEdges[region][1]), minEdge, maxEdge, binWidth, absMin, 20.0, 0.0,  0.8, njet, name=region, variable="SigFrac%s"%(self.translator[region]["A"]))
                                plotter[hist_key].plot_Var_vsDisc1Disc2(allRegionsSigFracs_TT[region]["B"][:,0],   edges, float(allRegionsFinalEdges[region][0]), float(allRegionsFinalEdges[region][1]), minEdge, maxEdge, binWidth, absMin, 20.0, 0.0,  0.8, njet, name=region, variable="SigFrac%s"%(self.translator[region]["B"]))
                                plotter[hist_key].plot_Var_vsDisc1Disc2(allRegionsSigFracs_TT[region]["C"][:,0],   edges, float(allRegionsFinalEdges[region][0]), float(allRegionsFinalEdges[region][1]), minEdge, maxEdge, binWidth, absMin, 20.0, 0.0,  0.8, njet, name=region, variable="SigFrac%s"%(self.translator[region]["C"]))
                                plotter[hist_key].plot_Var_vsDisc1Disc2(allRegionsSigFracs_TT[region]["D"][:,0],   edges, float(allRegionsFinalEdges[region][0]), float(allRegionsFinalEdges[region][1]), minEdge, maxEdge, binWidth, absMin, 20.0, 0.0,  0.8, njet, name=region, variable="SigFrac%s"%(self.translator[region]["D"]))
                            # inverse significance
                            #plotter[hist_key].plot_inverseSignificance_vsNonClosure(finalSignificance, finalNonClosure, significances, nonClosures, edges, allRegionsFinalEdges[region], njet, name=region)

                        if hist_key != self.sig:
                            plotter[hist_key].plot_Var_vsDisc1Disc2(nonClosures[:,0], edges, float(allRegionsFinalEdges[region][0]), float(allRegionsFinalEdges[region][1]), minEdge, maxEdge, binWidth, absMin, 20.0, 0.0,  0.5, njet, name=region, variable="NonClosure"   )
                            #plotter[hist_key].plot_Var_vsDisc1Disc2(nonClosures[:,1], edges, float(allRegionsFinalEdges[region][0]), float(allRegionsFinalEdges[region][1]), minEdge, maxEdge, binWidth, absMin, 20.0, 0.0,  0.5, njet, name=region, variable="NonClosureUnc")
                            plotter[hist_key].plot_Var_vsDisc1Disc2(pull[:,0],        edges, float(allRegionsFinalEdges[region][0]), float(allRegionsFinalEdges[region][1]), minEdge, maxEdge, binWidth, -20.0,  20.0, -5.0, 5.0, njet, name=region, variable="Pull"         )

                    EventsPerNjets[hist_key][njet] = allRegionsFinalEvents[hist_key]
                
                edgesPerNjets[njet] = allRegionsFinalEdges

            # ------------------------------------
            # plot Disc1s vs Disc2s with all edges
            # ------------------------------------
            if kwargs["plotDisc1VsDisc2"]:
                plotter["TT"].plot_Disc1VsDisc2(hist_lists[self.sig],  allRegionsFinalEdges, njet, tag = self.sig, name = "Val_BD_CD")
                plotter["TT"].plot_Disc1VsDisc2(hist_lists["TT"],      allRegionsFinalEdges, njet, tag = "TT",     name = "Val_BD_CD")

            # get the nEvents for each ABCD region
            sigFracsTable_AllRegions.writeLine(region="ABCD", njet=njet, finalSigFracs=allRegionsFinalSigFracs_TT["ABCD"], nEvents_AC=allRegionsFinalEvents["TT"]["ABCD"]["A"][0]+allRegionsFinalEvents["TT"]["ABCD"]["C"][0], nEvents_AB=allRegionsFinalEvents["TT"]["ABCD"]["A"][0]+allRegionsFinalEvents["TT"]["ABCD"]["B"][0])

            # get the nEvents for each B'D'EF region 
            sigFracsTable_AllRegions.writeLine(region="Val_BD", njet=njet, finalSigFracs=allRegionsFinalSigFracs_TT["Val_BD"], nEvents_AC=allRegionsFinalEvents["TT"]["Val_BD"]["A"][0]+allRegionsFinalEvents["TT"]["Val_BD"]["C"][0], nEvents_AB=allRegionsFinalEvents["TT"]["Val_BD"]["A"][0]+allRegionsFinalEvents["TT"]["Val_BD"]["B"][0])

            # get the nEvents for each C'D'GH region 
            sigFracsTable_AllRegions.writeLine(region="Val_CD", njet=njet, finalSigFracs=allRegionsFinalSigFracs_TT["Val_CD"], nEvents_AC=allRegionsFinalEvents["TT"]["Val_CD"]["A"][0]+allRegionsFinalEvents["TT"]["Val_CD"]["C"][0], nEvents_AB=allRegionsFinalEvents["TT"]["Val_CD"]["A"][0]+allRegionsFinalEvents["TT"]["Val_CD"]["B"][0])

            sigFracsTable_AllRegions.writeLine(region="Val_D", njet=njet, finalSigFracs=allRegionsFinalSigFracs_TT["Val_D"], nEvents_AC=allRegionsFinalEvents["TT"]["Val_D"]["A"][0]+allRegionsFinalEvents["TT"]["Val_D"]["C"][0], nEvents_AB=allRegionsFinalEvents["TT"]["Val_D"]["A"][0]+allRegionsFinalEvents["TT"]["Val_D"]["B"][0]) 

            abcdFracsTable.writeLine(njet=njet, finalTTfracs=allRegionsFinalTTFracs["ABCD"], finalSigFracs=allRegionsFinalSigFracs_TT["ABCD"])

            valFracsTable.writeLine(njet=njet, finalSigFracs=allRegionsFinalSigFracs_TT)

        # ----------------------
        # make all closure plots
        # ----------------------
        for region in regions:

            # usual closure
            plotter["TT"].make_allClosures(edgesPerNjets,    EventsPerNjets["TT"], None,                    None,                None,                   njets, name = region, closureTag = "b",  bkgTag = "TT")
            #plotter["TT"].make_allClosures(edgesPerNjets,    EventsPerNjets["TT"], None,                    EventsPerNjets[self.sig], None,              njets, name = region, closureTag = "sb", bkgTag = "TT")
            #plotter["NonTT"].make_allClosures(edgesPerNjets, None,                 EventsPerNjets["NonTT"], None,                None,                   njets, name = region, closureTag = "b",  bkgTag = "NonTT")
            #plotter["Data"].make_allClosures(edgesPerNjets,  None,                 None,                    None,                EventsPerNjets["Data"], njets, name = region, closureTag = "b",  bkgTag = "InitialData")
            #plotter[self.ttVar].make_ttVariances_allClosures(edgesPerNjets, EventsPerNjets["TT"], EventsPerNjets[self.ttVar], njets, name = region, closureTag = "b", varTag = self.ttVar)
            plotter["Data"].make_allClosures(edgesPerNjets, EventsPerNjets["TT"], EventsPerNjets["NonTT"], None, EventsPerNjets["Data"], njets, name = region, closureTag = "b", bkgTag = "TTinData")

            # data-MC closure
            #if region != "ABCD":
            #    # data-MC closure & pull
            #    plotter["Data"].make_allClosures(edgesPerNjets, EventsPerNjets["TT"], EventsPerNjets["NonTT"], None, EventsPerNjets["Data"], njets, name = region, closureTag = "b", bkgTag = "TTinData")
            #    # MC corrected data closure 
            #    plotter["Data"].make_allClosures(edgesPerNjets, EventsPerNjets[self.ttVar], EventsPerNjets["NonTT"], None, EventsPerNjets["Data"], njets, name = region, closureTag = "b", bkgTag = self.ttVar)

        # ------------------
        # make all tex files
        # ------------------
        sigFracsTable_AllRegions.writeClose()

        abcdFracsTable.writeClose()

        valFracsTable.writeClose()
