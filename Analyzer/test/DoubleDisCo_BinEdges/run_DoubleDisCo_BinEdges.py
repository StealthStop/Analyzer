import ROOT
import os
import argparse

from DoubleDisCo_Regions     import *
from DoubleDisCo_Plotter     import *
from DoubleDisCo_Variances   import *
from DoubleDisCo_TableWriter import *

def main():

    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # get the plots with non-optimized ABCD edges:
    #   -- python run_DoubleDisCo_BinEdges.py --year 2016 --tt TT --nontt NonTT --ttVar TT_erdOn --data Data --sig RPV --mass 550 --channel 0l
    #   -- python run_DoubleDisCo_BinEdges.py --year 2016 --tt TT --nontt NonTT --ttVar TT_erdOn --data Data --sig RPV --mass 550 --channel 1l
    # get the plots with fixed ABCD and Validation edges:
    #   -- python run_DoubleDisCo_BinEdges.py --year 2016 --tt TT --nontt NonTT --ttVar TT_erdOn --data Data --sig RPV --mass 550 --channel 0l --fixedDisc1edge 0.6 --fixedDisc2edge 0.6
    #   -- python run_DoubleDisCo_BinEdges.py --year 2016 --tt TT --nontt NonTT --ttVar TT_erdOn --data Data --sig RPV --mass 550 --channel 1l --fixedDisc1edge 0.6 --fixedDisc2edge 0.6
    # get the plots with optimized ABCD edges:
    #   -- first run the "run_DoubleDisCo_Optimized_BinEdges.py" to print the optimized ABCD edges and use them as fixed edges on the commend line
    #   -- python run_DoubleDisCo_BinEdges.py --year 2016 --tt TT --nontt NonTT --ttVar TT_erdOn --data Data --sig RPV --mass 550 --channel 0l --fixedDisc1edge 0.69 --fixedDisc2edge 0.74
    #   -- python run_DoubleDisCo_BinEdges.py --year 2016 --tt TT --nontt NonTT --ttVar TT_erdOn --data Data --sig RPV --mass 550 --channel 1l --fixedDisc1edge 0.59 --fixedDisc2edge 0.70
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
    usage  = "usage: %prog [options]"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("--year",             dest="year",             help="which year",                        required=True)
    parser.add_argument("--path",             dest="path",             help="Input dir with histos",             default="/uscms_data/d3/jhiltb/PO_Boxes/shared/2016_DisCo_0L_Cand1_1L") # Both 0L & 1l with OldSeed
    parser.add_argument("--tt",               dest="tt",               help="TT",                                required=True)
    parser.add_argument("--nontt",            dest="nontt",            help="NonTT",                             required=True)
    parser.add_argument("--ttVar",            dest="ttVar",            help="TT Variances",                      required=True)
    parser.add_argument("--sig",              dest="sig",              help="RPV, SYY",                          default="RPV")
    parser.add_argument("--mass",             dest="mass",             help="signal mass",                       default="550")
    parser.add_argument("--data",             dest="data",             help="JetHT, SingleLepton",               required=True)
    parser.add_argument("--channel",          dest="channel",          help="0l or 1l",                          required=True)
    parser.add_argument("--fixedDisc1edge",   dest="fixedDisc1edge",   help="fixed d1 edge",                     default=None,  type=float)
    parser.add_argument("--fixedDisc2edge",   dest="fixedDisc2edge",   help="fixed d2 edge",                     default=None,  type=float)
    parser.add_argument("--fastMode",         dest="fastMode",         help="Fast mode, don't scan all choices", default=False, action="store_true") 
    parser.add_argument("--njets",            dest="njets",            help="which njet bins to run on",         nargs="+",     default=["7", "8", "9", "10", "11", "12incl"], type=str) 
    parser.add_argument("--plotVars1D",       dest="plotVars1D",       help="Plot 1D var vs disc (slices)",      default=False, action="store_true") 
    parser.add_argument("--plotVars2D",       dest="plotVars2D",       help="Plot var vs disc1 and disc2 (2D)",  default=False, action="store_true") 
    parser.add_argument("--plotDisc1VsDisc2", dest="plotDisc1VsDisc2", help="Plot disc1 and disc2 (2D)",         default=False, action="store_true") 
    parser.add_argument("--makeEventsTables", dest="makeEventsTables", help="Make nEvents LaTeX tables",         default=False, action="store_true") 
    args = parser.parse_args()

    Sig     = "%s_%s"%(args.sig, args.mass)
    ttVar   = "%s"%(args.ttVar)

    # Names of samples/processes/data whose 2D disc1 vs disc2 histos will be analyzed
    samples = [args.tt, args.nontt, ttVar, Sig, args.data]

    # Make the output directories if they do not already exist
    plotsPath = {}; tablesPath = {}

    for sample in samples:

        # make directories to save plots and tables        
        if sample == Sig: continue

        if args.fixedDisc1edge != None or args.fixedDisc2edge != None:
            plotsPath[sample]  = "plots_BinEdges_%s_%s_%s/%s_%s/%s/"%(args.fixedDisc1edge, args.fixedDisc2edge, sample, args.sig, args.mass, args.channel)
            
            if sample == "TT":
                tablesPath[sample] = "tables_BinEdges_%s_%s_%s/%s"%(args.fixedDisc1edge, args.fixedDisc2edge, sample, args.channel)

        # Make all paths for saving plots and tables
        # if paths do not already exist
        if not os.path.exists(plotsPath[sample]):
            os.makedirs(plotsPath[sample])

        if not os.path.exists(tablesPath["TT"]):
            os.makedirs(tablesPath["TT"])

    # get SYY histogram
    modelDecay = "2t6j"
    if ("SHH" in args.sig):
        modelDecay = "2t4b"

    # root files
    files = {
        "TT"                   : ROOT.TFile.Open(args.path + "/" + args.year + "_TT.root"),
        "TT_erdOn"             : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_erdOn.root"),
        "TT_fsrDown"           : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_fsrDown.root"),
        "TT_fsrUp"             : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_fsrUp.root"),
        "TT_isrDown"           : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_isrDown.root"),
        "TT_isrUp"             : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_isrUp.root"),
        "TT_hdampDown"         : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_hdampDown.root"),
        "TT_hdampUp"           : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_hdampUp.root"),
        "TT_JECdown"           : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_JECdown.root"),
        "TT_JECup"             : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_JECup.root"),
        "TT_JERdown"           : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_JERdown.root"),
        "TT_JERup"             : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_JERup.root"),
        "TT_underlyingEvtDown" : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_underlyingEvtDown.root"),
        "TT_underlyingEvtUp"   : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_underlyingEvtUp.root"), 
        "TT_UL"                : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_UL.root"), 
        "NonTT"                : ROOT.TFile.Open(args.path + "/" + args.year + "_Non_TT.root"),
        #"NonTT_withoutQCD"     : ROOT.TFile.Open(args.path + "/" + args.year + "_NonTT_withoutQCD.root"),
        "Data"                 : ROOT.TFile.Open(args.path + "/" + args.year + "_Data.root"),
        Sig                    : ROOT.TFile.Open(args.path + "/" + args.year + "_%s_%s_mStop-%s.root"%(args.sig,modelDecay,args.mass)),
    }

    # get the 2D histograms 
    histNames = "h_DoubleDisCo_disc1_disc2_%s_Njets"%(args.channel)

    # Njets bins to process
    njets = args.njets

    # ------------------
    # make all tex files
    # ------------------
    # get the signal fracs for each region
    sigFracsTable_AllRegions = SignalFractionsAllRegionsTable(tablesPath["TT"], args.channel, args.year, "Sig_fracs_fixed_AllRegions", args.sig, args.mass)

    # get the fracs for each ABCD region
    abcdFracsTable = ABCDfracsTable(tablesPath["TT"], args.channel, args.year, "Sig_TT_fracs_fixed_ABCD", args.sig, args.mass)

    valFracsTable  = ValFracsTable(tablesPath["TT"], args.channel, args.year, "Sig_fracs_fixed_Val", args.sig, args.mass)

    # hold on edges per njet
    edgesPerNjets = {njet : None for njet in njets}

    # make regionis list for adding all edges to DoubleDisCo cfg file
    regions = ["ABCD",
               "Val_BD",
               "Val_CD",
               "Val_D", 
    ]

    # initialize the dictionaries of any regions
    translator = {"ABCD"   : {"A" : "A",  "B" : "B",  "C" : "C",  "D" : "D" },
                  "Val_BD" : {"A" : "b",  "B" : "E",  "C" : "d",  "D" : "F" },
                  "Val_CD" : {"A" : "c",  "B" : "di", "C" : "G",  "D" : "H" },
                  "Val_D"  : {"A" : "dA", "B" : "dB", "C" : "dC", "D" : "dD"},
    }

    # --------------------------------
    # Common calculations and plotters
    # --------------------------------
    plotter = {}; plotter_ttVar = {}; EventsPerNjets = {}

    for sample in samples:

        # Hold onto per Njet things so we can plot altogether after the Njets loop
        EventsPerNjets[sample] = {njet : None for njet in njets}

        if sample == Sig: continue

        if sample == ttVar:
            plotter[sample]       = Common_Calculations_Plotters(plotsPath[sample], args.year, args.sig, args.mass, args.channel)
            plotter_ttVar[sample] = ttVariances_Plotters(plotsPath[sample], args.year, args.channel)

        else:
            plotter[sample] = Common_Calculations_Plotters(plotsPath[sample], args.year, args.sig, args.mass, args.channel)

    # ---------------
    # loop over njets
    # --------------- 
    for njet in njets: 

        hist_lists = {}

        for sample in samples:

            hist_lists[sample] = files[sample].Get(histNames + njet)

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

                theEdgesClass = All_Regions(hist_lists, Sig=Sig, ttVar=ttVar, fixedDisc1Edge=args.fixedDisc1edge, fixedDisc2Edge=args.fixedDisc2edge, fastMode=args.fastMode)

            # --------------------------
            #        *    ||   |                       
            #    vB  * vA ||   |    A                
            #  ______*____||___|________             
            #        *    ||   |                     
            #    vD  * vC ||   |    C                
            #        *                  
            # --------------------------             
            elif region == "Val_BD":

               theEdgesClass = All_Regions(hist_lists, Sig=Sig, ttVar=ttVar, fixedDisc2Edge=allRegionsFinalEdges["ABCD"][1], rightBoundary=float(0.4), fixedDisc1Edge=float(0.2), fastMode=args.fastMode)

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

                theEdgesClass = All_Regions(hist_lists, Sig=Sig, ttVar=ttVar, fixedDisc1Edge=allRegionsFinalEdges["ABCD"][0], topBoundary=float(0.4), fixedDisc2Edge=float(0.2), fastMode=args.fastMode)

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
                
                theEdgesClass = All_Regions(hist_lists, Sig=Sig, ttVar=ttVar, rightBoundary=allRegionsFinalEdges["ABCD"][0], topBoundary=allRegionsFinalEdges["ABCD"][1], ABCDdisc1=float(allRegionsFinalEdges["ABCD"][0])/2.0, ABCDdisc2=float(allRegionsFinalEdges["ABCD"][1])/2.0, fastMode=args.fastMode)

            # ----------------------------
            # Optimization based on TT !!!
            # ----------------------------
            # get all final bin edges
            allRegionsFinalEdges[region] = theEdgesClass.getFinal("edges", "TT")

            # quantities and variables with any combination of bin edges
            significances                     = theEdgesClass.get("significance",                     None, None, "TT")
            significances_includingNonClosure = theEdgesClass.get("significance_includingNonClosure", None, None, "TT")
            significances_nonSimplified       = theEdgesClass.get("significance_nonSimplified",       None, None, "TT")

            allRegionsSigFracs_TT[region]     = {"A" : theEdgesClass.get("sigFractionA", None, None, Sig),
                                                 "B" : theEdgesClass.get("sigFractionB", None, None, Sig),
                                                 "C" : theEdgesClass.get("sigFractionC", None, None, Sig),
                                                 "D" : theEdgesClass.get("sigFractionD", None, None, Sig)}
            allRegionsTTFracs[region]         = {"A" : theEdgesClass.get("ttFractionA", None, None, "TT"),
                                                 "B" : theEdgesClass.get("ttFractionB", None, None, "TT"),
                                                 "C" : theEdgesClass.get("ttFractionC", None, None, "TT"),
                                                 "D" : theEdgesClass.get("ttFractionD", None, None, "TT")}

            # quantities and variables with the final choice of bin edges
            finalSignificance               = theEdgesClass.getFinal("significance", "TT")
            allRegionsFinalSigFracs_TT[region] = {"A" : theEdgesClass.getFinal("sigFractionA", Sig),
                                                  "B" : theEdgesClass.getFinal("sigFractionB", Sig),
                                                  "C" : theEdgesClass.getFinal("sigFractionC", Sig),
                                                  "D" : theEdgesClass.getFinal("sigFractionD", Sig)}
            allRegionsFinalTTFracs[region]     = {"A" : theEdgesClass.getFinal("ttFractionA", "TT"),
                                                  "B" : theEdgesClass.getFinal("ttFractionB", "TT"),
                                                  "C" : theEdgesClass.getFinal("ttFractionC", "TT"),
                                                  "D" : theEdgesClass.getFinal("ttFractionD", "TT")}
            
            # -----------------------------------------------
            # loop over for getting plots for TT, NonTT, Data
            # -----------------------------------------------
            for hist_key in hist_lists.keys():

                edges                                = np.array(theEdgesClass.get("edges"), dtype=float)
                allRegionsEvents[hist_key][region]      = {"A" : theEdgesClass.get("nEventsA", None, None, hist_key)}
                allRegionsFinalEvents[hist_key][region] = {"A" : theEdgesClass.getFinal("nEventsA", hist_key),
                                                           "B" : theEdgesClass.getFinal("nEventsB", hist_key),
                                                           "C" : theEdgesClass.getFinal("nEventsC", hist_key),
                                                           "D" : theEdgesClass.getFinal("nEventsD", hist_key)}
                if hist_key != Sig:
                    nonClosures      = theEdgesClass.get("nonClosure",      None, None, hist_key) # vars with any combination of bin edges
                    pull             = theEdgesClass.get("pull",            None, None, hist_key)
                    finalNonClosure  = theEdgesClass.getFinal("nonClosure",             hist_key) # vars with the final choice of bin edges
                    finalPull        = theEdgesClass.getFinal("pull",                   hist_key)
               
                # ---------------------------  
                # plot variable vs disc as 1D
                # --------------------------- 
                if args.plotVars1D:
                    for disc in [1, 2]:

                        if hist_key == Sig:
                            plotter["TT"].plot_VarVsDisc(allRegionsEvents[hist_key][region]["A"], edges, binWidth/2.0, -1.0, "Weighted Signal Events", "wSigEvts", disc, njet, name = region)

                        elif hist_key != ttVar:
                            plotter[hist_key].plot_VarVsDisc(allRegionsEvents[hist_key][region]["A"], edges, binWidth/2.0, -1.0, "Weighted Background Events", "wBkgEvts", disc, njet, name = region)
                            plotter[hist_key].plot_VarVsDisc(nonClosures,                          edges, binWidth/2.0, 1.0,  "Closure",                    "Closure",  disc, njet, name = region)

                        if hist_key == "TT":
                            plotter[hist_key].plot_VarVsDisc(significances, edges, binWidth/2.0, 5.0, "%s Significance"%(region), "Significance", disc, njet, name = region)

                            for subregion in translator["ABCD"].keys():
                                plotter[hist_key].plot_VarVsDisc(allRegionsSigFracs_TT[region][subregion], edges, binWidth/2.0, 0.8, "Signal Contamination %s"%(translator[region][subregion]), "SigFracs%s"%(translator[region][subregion]), disc, njet, name = region)

                # ---------------------------
                # plot variable vs Disc1Disc2
                # ---------------------------
                if args.plotVars2D:
                    if hist_key == "TT":
                        plotter[hist_key].plot_Var_vsDisc1Disc2(significances[:,0],                     edges, float(allRegionsFinalEdges[region][0]), float(allRegionsFinalEdges[region][1]), minEdge, maxEdge, binWidth, 10e-10, 20.0, 0.0,  5.0, njet, name=region, variable="Sign"                            )
                        plotter[hist_key].plot_Var_vsDisc1Disc2(significances[:,1],                     edges, float(allRegionsFinalEdges[region][0]), float(allRegionsFinalEdges[region][1]), minEdge, maxEdge, binWidth, 10e-10, 20.0, 0.0,  5.0, njet, name=region, variable="SignUnc"                         )
                        # significance including non-closure
                        plotter[hist_key].plot_Var_vsDisc1Disc2(significances_includingNonClosure[:,0], edges, float(allRegionsFinalEdges[region][0]), float(allRegionsFinalEdges[region][1]), minEdge, maxEdge, binWidth, 10e-10, 20.0, 0.0,  5.0, njet, name=region, variable="Sign_includingNonClosure"        )
                        plotter[hist_key].plot_Var_vsDisc1Disc2(significances_includingNonClosure[:,1], edges, float(allRegionsFinalEdges[region][0]), float(allRegionsFinalEdges[region][1]), minEdge, maxEdge, binWidth, 10e-10, 20.0, 0.0,  5.0, njet, name=region, variable="SignUnc_includingNonClosure"     )
                        # significance, non-simplified version
                        #plotter[hist_key].plot_Var_vsDisc1Disc2(significances_nonSimplified[:,0],       edges, float(allRegionsFinalEdges[region][0]), float(allRegionsFinalEdges[region][1]), minEdge, maxEdge, binWidth, 10e-10, 20.0, 0.0,  5.0, njet, name=region, variable="Sign_nonSimplified"              )
                        #plotter[hist_key].plot_Var_vsDisc1Disc2(allRegionsSigFracs_TT[key]["A"][:,0],   edges, float(allRegionsFinalEdges[region][0]), float(allRegionsFinalEdges[region][1]), minEdge, maxEdge, binWidth, 10e-10, 20.0, 0.0,  0.8, njet, name=region, variable="SigFrac%s"%(translator[region]["A"]))
                        #plotter[hist_key].plot_Var_vsDisc1Disc2(allRegionsSigFracs_TT[key]["B"][:,0],   edges, float(allRegionsFinalEdges[region][0]), float(allRegionsFinalEdges[region][1]), minEdge, maxEdge, binWidth, 10e-10, 20.0, 0.0,  0.8, njet, name=region, variable="SigFrac%s"%(translator[region]["B"]))
                        #plotter[hist_key].plot_Var_vsDisc1Disc2(allRegionsSigFracs_TT[key]["C"][:,0],   edges, float(allRegionsFinalEdges[region][0]), float(allRegionsFinalEdges[region][1]), minEdge, maxEdge, binWidth, 10e-10, 20.0, 0.0,  0.8, njet, name=region, variable="SigFrac%s"%(translator[region]["C"]))
                        #plotter[hist_key].plot_Var_vsDisc1Disc2(allRegionsSigFracs_TT[key]["D"][:,0],   edges, float(allRegionsFinalEdges[region][0]), float(allRegionsFinalEdges[region][1]), minEdge, maxEdge, binWidth, 10e-10, 20.0, 0.0,  0.8, njet, name=region, variable="SigFrac%s"%(translator[region]["D"]))
                        #plotter[hist_key].plot_inverseSignificance_vsNonClosure(finalSignificance, finalNonClosure, significances, nonClosures, edges, allRegionsFinalEdges[region], njet, name=region)

                    if hist_key != Sig:
                        plotter[hist_key].plot_Var_vsDisc1Disc2(nonClosures[:,0], edges, float(allRegionsFinalEdges[region][0]), float(allRegionsFinalEdges[region][1]), minEdge, maxEdge, binWidth, 10e-10, 20.0, 0.0,  0.5, njet, name=region, variable="NonClosure"   )
                        plotter[hist_key].plot_Var_vsDisc1Disc2(nonClosures[:,1], edges, float(allRegionsFinalEdges[region][0]), float(allRegionsFinalEdges[region][1]), minEdge, maxEdge, binWidth, 10e-10, 20.0, 0.0,  0.5, njet, name=region, variable="NonClosureUnc")
                        plotter[hist_key].plot_Var_vsDisc1Disc2(pull[:,0],        edges, float(allRegionsFinalEdges[region][0]), float(allRegionsFinalEdges[region][1]), minEdge, maxEdge, binWidth, -20.0,  20.0, -5.0, 5.0, njet, name=region, variable="Pull"         )

                EventsPerNjets[hist_key][njet] = allRegionsFinalEvents[hist_key]
            
            edgesPerNjets[njet] = allRegionsFinalEdges

        # ------------------------------------
        # plot Disc1s vs Disc2s with all edges
        # ------------------------------------
        if args.plotDisc1VsDisc2:
            plotter["TT"].plot_Disc1VsDisc2(hist_lists[Sig],  allRegionsFinalEdges, njet, tag = "RPV550", name = "Val_BD_CD", col1="yellow", col2="lime" )
            plotter["TT"].plot_Disc1VsDisc2(hist_lists["TT"], allRegionsFinalEdges, njet, tag = "TT",     name = "Val_BD_CD", col1="yellow", col2="lime" )
            plotter["TT"].plot_Disc1VsDisc2(hist_lists[Sig],  allRegionsFinalEdges, njet, tag = "RPV550", name = "Val_D",     col1="white",  col2="white")
            plotter["TT"].plot_Disc1VsDisc2(hist_lists["TT"], allRegionsFinalEdges, njet, tag = "TT",     name = "Val_D",     col1="white",  col2="white")

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
        #plotter["TT"].make_allClosures(edgesPerNjets,    EventsPerNjets["TT"], None,                    None,                None,                   njets, name = region, closureTag = "b",  bkgTag = "TT")
        #plotter["TT"].make_allClosures(edgesPerNjets,    EventsPerNjets["TT"], None,                    EventsPerNjets[Sig], None,                   njets, name = region, closureTag = "sb", bkgTag = "TT")
        #plotter["NonTT"].make_allClosures(edgesPerNjets, None,                 EventsPerNjets["NonTT"], None,                None,                   njets, name = region, closureTag = "b",  bkgTag = "NonTT")
        #plotter["TT"].make_allClosures(edgesPerNjets,    EventsPerNjets["TT"], EventsPerNjets["NonTT"], None,                EventsPerNjets["Data"], njets, name = region, closureTag = "b",  bkgTag = "AllBkg")
        #plotter["Data"].make_allClosures(edgesPerNjets,  None,                 None,                    None,                EventsPerNjets["Data"], njets, name = region, closureTag = "b",  bkgTag = "InitialData")
        #plotter_ttVar[ttVar].make_ttVariances_allClosures(edgesPerNjets, EventsPerNjets["TT"], EventsPerNjets[ttVar], njets, name = region, closureTag = "b", varTag = ttVar)

        # data-MC closure
        if region != "ABCD":
            # data-MC closure & pull
            # MC corrected data closure 
            plotter["Data"].make_allClosures(edgesPerNjets, EventsPerNjets["TT"], EventsPerNjets["NonTT"], None, EventsPerNjets["Data"], njets, name = region, closureTag = "b", bkgTag = "TTinData")
            plotter["Data"].make_allClosures(edgesPerNjets, EventsPerNjets[ttVar], EventsPerNjets["NonTT"], None, EventsPerNjets["Data"], njets, name = region, closureTag = "b", bkgTag = ttVar)

    # ------------------
    # make all tex files
    # ------------------
    sigFracsTable_AllRegions.writeClose()

    abcdFracsTable.writeClose()

    valFracsTable.writeClose()

if __name__ == '__main__':
    main()
