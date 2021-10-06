import ROOT
import os
import argparse

from DoubleDisCo_Regions     import *
from DoubleDisCo_Plotter     import *
from DoubleDisCo_TableWriter import *

def main():

    # -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # command to run this script
    #   -- python run_DoubleDisCo_BinEdges.py --year 2016 --tt TT --nontt NonTT --data Data --sig RPV --mass 550 --channel 0l --metric New  --fixedDisc1edge 0.6 --fixedDisc2edge 0.6
    #   -- python run_DoubleDisCo_BinEdges.py --year 2016 --tt TT --nontt NonTT --data Data --sig RPV --mass 550 --channel 1l --metric New  --fixedDisc1edge 0.6 --fixedDisc2edge 0.6
    #   -- python run_DoubleDisCo_BinEdges.py --year 2016 --tt TT --nontt NonTT --data Data --sig RPV --mass 550 --channel 0l --metric New  
    #   -- python run_DoubleDisCo_BinEdges.py --year 2016 --tt TT --nontt NonTT --data Data --sig RPV --mass 550 --channel 1l --metric New  
    # --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
    usage  = "usage: %prog [options]"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("--year",           dest="year",           help="which year",            required=True)
    parser.add_argument("--path",           dest="path",           help="Input dir with histos", default="/uscms/home/jhiltb/nobackup/PO_Boxes/DoubleDisCo_Reg_0L_loose_1L_RPV_2016_20210920_Output/")
    parser.add_argument("--tt",             dest="tt",             help="TT",                    required=True)
    parser.add_argument("--nontt",          dest="nontt",          help="NonTT",                 required=True)
    parser.add_argument("--sig",            dest="sig",            help="RPV, SYY",              default="RPV")
    parser.add_argument("--mass",           dest="mass",           help="signal mass",           default="550")
    parser.add_argument("--data",           dest="data",           help="JetHT, SingleLepton",   required=True)
    parser.add_argument("--channel",        dest="channel",        help="0l or 1l",              required=True)
    parser.add_argument("--metric",         dest="metric",         help="NN,New",                required=True, type=str)
    parser.add_argument("--fixedDisc1edge", dest="fixedDisc1edge", help="fixed d1 edge",         default=None, type=float)
    parser.add_argument("--fixedDisc2edge", dest="fixedDisc2edge", help="fixed d2 edge",         default=None, type=float)
    args = parser.parse_args()

    Sig = "%s_%s"%(args.sig, args.mass)
     
    samples = [args.tt, args.nontt, Sig, args.data]

    # Make the output directories if they do not already exist
    plotsPath = {}; tablesPath = {}

    for sample in samples:

        # make directories to save plots and tables        
        if sample == Sig: continue

        if args.fixedDisc1edge != None or args.fixedDisc2edge != None:
            plotsPath[sample]  = "plots_fixedEdges_%s_%s/%s_%s/%s/"%(args.fixedDisc1edge, sample, args.sig, args.mass, args.channel)
            
            if sample == "TT":
                tablesPath[sample] = "tables_fixedEdges_%s_%s/%s"%(args.fixedDisc1edge, sample, args.channel)

        else:
            plotsPath[sample]  = "plots_%s/%s_%s/%s/"%(sample, args.sig, args.mass, args.channel)

            if sample == "TT":
                tablesPath[sample] = "tables_%s/%s"%(sample, args.channel)

        # 
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
        "TT"               : ROOT.TFile.Open(args.path + "/" + args.year + "_TT.root"),
        "NonTT"            : ROOT.TFile.Open(args.path + "/" + args.year + "_Non_TT.root"),
        "NonTT_withoutQCD" : ROOT.TFile.Open(args.path + "/" + args.year + "_NonTT_withoutQCD.root"),
        "Data"             : ROOT.TFile.Open(args.path + "/" + args.year + "_Data.root"),
        Sig                : ROOT.TFile.Open(args.path + "/" + args.year + "_%s_%s_mStop-%s.root"%(args.sig,modelDecay,args.mass)),
    }

    # get the 2D histograms 
    njets = None
    histNames = "h_DoubleDisCo_disc1_disc2_%s_Njets"%(args.channel)
    
    # for 0-lepton 
    if args.channel == "0l":
        njets = ["6", "7", "8", "9", "10", "11", "12"]

    # for 1-lepton
    else:
        njets = ["7", "8", "9", "10", "11"]
    
    # ------------------
    # make all tex files
    # ------------------
    tag = ""
    if args.fixedDisc1edge != None or args.fixedDisc2edge != None:
        tag = "fixed"
    else:
        tag = "final"

    # put all latest bin edges to tex file
    binEdgesTable = BinEdgesTable(tablesPath["TT"], args.channel, args.year, "All_%sBinEdges"%(tag), args.sig, args.mass)

    # get the nEvents for each ABCD region
    abcdEventsTable = ABCDeventsTable(tablesPath["TT"], args.channel, args.year, "nEvents_%s_ABCD"%(tag), args.sig, args.mass)

    # get the nEvents for each B'D'EF region 
    bdefEventsTable       = BDEFeventsTable(tablesPath["TT"], args.channel, args.year, "nEvents_%s_Val_bdEF"%(tag), args.sig, args.mass)
    fixedBDEF_EventsTable = BDEFeventsTable(tablesPath["TT"], args.channel, args.year, "nEvents_%s_Val_fixedBDEF"%(tag), args.sig, args.mass)

    # get the nEvents for each C'D'GH region
    cdghEventsTable       = CDGHeventsTable(tablesPath["TT"], args.channel, args.year, "nEvents_%s_Val_cdiGH"%(tag), args.sig, args.mass)
    fixedCDGH_EventsTable = CDGHeventsTable(tablesPath["TT"], args.channel, args.year, "nEvents_%s_Val_fixedCDGH"%(tag), args.sig, args.mass)

    # get the table for Validation region sub-division D
    subdBinEdgesTable = SubDivDeventsTable(tablesPath["TT"], args.channel, args.year, "%s_FinalBinEdges_Val_subDivD"%(tag), args.sig, args.mass)

    # get the table for bkg and sig+bkg events in ABCD regions
    bkgEventsTable = BkgTotEvents(tablesPath["TT"], args.channel, args.year, "TT_Fracs_%s_ABCD"%(tag), args.sig, args.mass)

    # hold on edges per njet
    edgesPerNjets = {njet : None for njet in njets}

    # make regionis list for adding all edges to DoubleDisCo cfg file
    regions = {"ABCD"          : "ABCD",
               "Val_bdEF"      : "bdEF",
               "Val_fixedBDEF" : "fixedBDEF",
               "Val_cdiGH"     : "cdGH",
               "Val_fixedCDGH" : "fixedCDGH",
               "Val_subDivD"   : "subDivD",
    }

    # initialize the dictionaries of any regions
    translator = {"ABCD"          : {"A" : "A",  "B" : "B",  "C" : "C",  "D" : "D" },
                  "Val_bdEF"      : {"A" : "b",  "B" : "E",  "C" : "d",  "D" : "F" },
                  "Val_cdiGH"     : {"A" : "c",  "B" : "di", "C" : "G",  "D" : "H" },
                  "Val_subDivD"   : {"A" : "dA", "B" : "dB", "C" : "dC", "D" : "dD"},
                  "Val_fixedBDEF" : {"A" : "vA", "B" : "vB", "C" : "vC", "D" : "vD"},
                  "Val_fixedCDGH" : {"A" : "hA", "B" : "hB", "C" : "hC", "D" : "hD"},
    }

    # --------------------------------
    # Common calculations and plotters
    # --------------------------------
    plotter = {}; EventsPerNjets = {}

    for sample in samples:

        # Hold onto per Njet things so we can plot altogether after the Njets loop
        EventsPerNjets[sample] = {njet : None for njet in njets}

        if sample == Sig: continue

        plotter[sample] = Common_Calculations_Plotters(plotsPath[sample], args.metric, args.year, args.sig, args.mass, args.channel)

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
        for key, region in regions.items():

            allRegionsSigFracs_TT[key]      = {"A" : {}, "B" : {}, "C" : {}, "D" : {}}  
            allRegionsFinalSigFracs_TT[key] = {"A" : {}, "B" : {}, "C" : {}, "D" : {}}
            allRegionsTTFracs[key]          = {"A" : {}, "B" : {}, "C" : {}, "D" : {}}
            allRegionsFinalTTFracs[key]     = {"A" : {}, "B" : {}, "C" : {}, "D" : {}}

            # loop over for TT, NonTT, Data to get all regions' events
            for hist_key in hist_lists.keys():

                allRegionsEvents[hist_key][key]      = {"A" : {}, "B" : {}, "C" : {}, "D" : {}}
                allRegionsFinalEvents[hist_key][key] = {"A" : {}, "B" : {}, "C" : {}, "D" : {}}

        # ------------------------------------------------------------------
        # Loop through the regions and make the set of plots for each
        # Make sure ABCD goes first so that the val regions can use its info
        # ------------------------------------------------------------------
        for key, region in regions.items():

            theEdgesClass = None

            # -----------------
            # make ABCD regions
            # -----------------
            if key == "ABCD":

                theEdgesClass = ABCDedges(hist_lists, Sig=Sig, fixedDisc1Edge=args.fixedDisc1edge, fixedDisc2Edge=args.fixedDisc2edge, metric=args.metric)

            # ---------------------------------------------------------------------------
            # make the validation regions in BD : bdEF
            #   -- define validation disc1 edge between 0 and the final bin edge of disc1
            #   -- final disc2 edge is constant
            #          *                              
            #      E   *    B'  |    A                
            #   _______*________|________             
            #          *        |                     
            #      F   *    D'  |    C                
            #          *                     
            # ---------------------------------------------------------------------------
            elif key == "Val_bdEF":

                theEdgesClass = bdEFedges(hist_lists, Sig=Sig, fixedDisc2Edge=allRegionsFinalEdges["ABCD"][1], rightBoundary=allRegionsFinalEdges["ABCD"][0], nBkgA=allRegionsFinalEvents["TT"]["ABCD"]["A"][0], nBkgC=allRegionsFinalEvents["TT"]["ABCD"]["C"][0], metric=args.metric)

            # ---------------------------------------------------------------------------
            # make the validation regions in CD : cdiGH
            #   -- define validation disc2 edge between 0 and the final bin edge of disc2
            #   -- final disc1 edge is constant
            #
            #          B  |  A   
            #       ______|______
            #             |
            #          D' |  C' 
            #       *************
            #             |
            #          H  |  G
            # --------------------------------------------------------------------------- 
            elif key == "Val_cdiGH":

                theEdgesClass = cdGHedges(hist_lists, Sig=Sig, fixedDisc1Edge=allRegionsFinalEdges["ABCD"][0], topBoundary=allRegionsFinalEdges["ABCD"][1], nBkgA=allRegionsFinalEvents["TT"]["ABCD"]["A"][0], nBkgB=allRegionsFinalEvents["TT"]["ABCD"]["B"][0], metric=args.metric)

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
            elif key == "Val_subDivD":
                
                theEdgesClass = subDivDedges(hist_lists, Sig=Sig, rightBoundary=allRegionsFinalEdges["ABCD"][0], topBoundary=allRegionsFinalEdges["ABCD"][1], ABCDdisc1=float(allRegionsFinalEdges["ABCD"][0]), ABCDdisc2=float(allRegionsFinalEdges["ABCD"][1]), metric=args.metric)

            # --------------------------
            #        *    ||   |                       
            #    vB  * vA ||   |    A                
            #  ______*____||___|________             
            #        *    ||   |                     
            #    vD  * vC ||   |    C                
            #        *                  
            # --------------------------             
            elif key == "Val_fixedBDEF":

               theEdgesClass = bdEFedges(hist_lists, Sig=Sig, fixedDisc2Edge=allRegionsFinalEdges["ABCD"][1], rightBoundary=float(0.4), fixedDisc1Edge=float(0.2), metrc=args.metric)

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
            elif key == "Val_fixedCDGH":

                theEdgesClass = cdGHedges(hist_lists, Sig=Sig, fixedDisc1Edge=allRegionsFinalEdges["ABCD"][0], topBoundary=float(0.4), fixedDisc2Edge=float(0.2), metric=args.metric)

            # ----------------------------
            # Optimization based on TT !!!
            # ----------------------------
            # get all final bin edges
            allRegionsFinalEdges[key] = theEdgesClass.getFinal("edges", "TT")

            # quantities and variables with any combination of bin edges
            significances              = theEdgesClass.get("significance", None, None, "TT")
            allRegionsSigFracs_TT[key] = {"A" : theEdgesClass.get("sigFractionA", None, None, Sig),
                                          "B" : theEdgesClass.get("sigFractionB", None, None, Sig),
                                          "C" : theEdgesClass.get("sigFractionC", None, None, Sig),
                                          "D" : theEdgesClass.get("sigFractionD", None, None, Sig)}
            allRegionsTTFracs[key]     = {"A" : theEdgesClass.get("ttFractionA", None, None, "TT"),
                                          "B" : theEdgesClass.get("ttFractionB", None, None, "TT"),
                                          "C" : theEdgesClass.get("ttFractionC", None, None, "TT"),
                                          "D" : theEdgesClass.get("ttFractionD", None, None, "TT")}

            # quantities and variables with the final choice of bin edges
            finalSignificance               = theEdgesClass.getFinal("significance", "TT")
            allRegionsFinalSigFracs_TT[key] = {"A" : theEdgesClass.getFinal("sigFractionA", Sig),
                                               "B" : theEdgesClass.getFinal("sigFractionB", Sig),
                                               "C" : theEdgesClass.getFinal("sigFractionC", Sig),
                                               "D" : theEdgesClass.getFinal("sigFractionD", Sig)}
            allRegionsFinalTTFracs[key]     = {"A" : theEdgesClass.getFinal("ttFractionA", "TT"),
                                               "B" : theEdgesClass.getFinal("ttFractionB", "TT"),
                                               "C" : theEdgesClass.getFinal("ttFractionC", "TT"),
                                               "D" : theEdgesClass.getFinal("ttFractionD", "TT")}
            
            # -----------------------------------------------
            # loop over for getting plots for TT, NonTT, Data
            # -----------------------------------------------
            for hist_key in hist_lists.keys():

                edges                                = np.array(theEdgesClass.get("edges"), dtype=float)
                allRegionsEvents[hist_key][key]      = {"A" : theEdgesClass.get("nEventsA", None, None, hist_key)}
                allRegionsFinalEvents[hist_key][key] = {"A" : theEdgesClass.getFinal("nEventsA", hist_key),
                                                        "B" : theEdgesClass.getFinal("nEventsB", hist_key),
                                                        "C" : theEdgesClass.getFinal("nEventsC", hist_key),
                                                        "D" : theEdgesClass.getFinal("nEventsD", hist_key)}
                if hist_key != Sig:
                    closureErrs      = theEdgesClass.get("closureError",      None, None, hist_key) # vars with any combination of bin edges
                    pull             = theEdgesClass.get("pull",              None, None, hist_key)
                    finalClosureErr  = theEdgesClass.getFinal("closureError", hist_key            ) # vars with the final choice of bin edges
                    finalPull        = theEdgesClass.getFinal("pull",         hist_key            )
                
                # plot variable vs disc as 1D 
                for disc in [1, 2]:

                    if hist_key == Sig:
                        plotter["TT"].plot_VarVsDisc(allRegionsEvents[hist_key][key]["A"], edges, binWidth/2.0, -1.0, "Weighted Signal Events",     "wSigEvts", disc, njet, name = key)

                    else:
                        plotter[hist_key].plot_VarVsDisc(allRegionsEvents[hist_key][key]["A"], edges, binWidth/2.0, -1.0, "Weighted Background Events", "wBkgEvts", disc, njet, name = key)
                        plotter[hist_key].plot_VarVsDisc(closureErrs,                          edges, binWidth/2.0, 1.0,  "%s Closure"%(key),           "Closure",  disc, njet, name = key)

                    if hist_key == "TT":
                        plotter[hist_key].plot_VarVsDisc(significances, edges, binWidth/2.0, 5.0, "%s Significance"%(key), "Significance", disc, njet, name = key)

                        for subregion in translator["ABCD"].keys():
                            plotter[hist_key].plot_VarVsDisc(allRegionsSigFracs_TT[key][subregion], edges, binWidth/2.0, 0.8, "Signal Contamination %s"%(translator[key][subregion]), "SigFracs%s"%(translator[key][subregion]), disc, njet, name = key)

                # plot variable vs Disc1Disc2
                if hist_key == "TT":
                    plotter[hist_key].plot_Var_vsDisc1Disc2(significances[:,0],                   edges, float(allRegionsFinalEdges[key][0]), float(allRegionsFinalEdges[key][1]), minEdge, maxEdge, binWidth, 10e-10, 20.0, 0.0,  5.0, njet, name=key, variable="Sign"                            )
                    plotter[hist_key].plot_Var_vsDisc1Disc2(allRegionsSigFracs_TT[key]["A"][:,0], edges, float(allRegionsFinalEdges[key][0]), float(allRegionsFinalEdges[key][1]), minEdge, maxEdge, binWidth, 10e-10, 20.0, 0.0,  0.8, njet, name=key, variable="SigFrac%s"%(translator[key]["A"]))
                    plotter[hist_key].plot_Var_vsDisc1Disc2(allRegionsSigFracs_TT[key]["B"][:,0], edges, float(allRegionsFinalEdges[key][0]), float(allRegionsFinalEdges[key][1]), minEdge, maxEdge, binWidth, 10e-10, 20.0, 0.0,  0.8, njet, name=key, variable="SigFrac%s"%(translator[key]["B"]))
                    plotter[hist_key].plot_Var_vsDisc1Disc2(allRegionsSigFracs_TT[key]["C"][:,0], edges, float(allRegionsFinalEdges[key][0]), float(allRegionsFinalEdges[key][1]), minEdge, maxEdge, binWidth, 10e-10, 20.0, 0.0,  0.8, njet, name=key, variable="SigFrac%s"%(translator[key]["C"]))
                    plotter[hist_key].plot_Var_vsDisc1Disc2(allRegionsSigFracs_TT[key]["D"][:,0], edges, float(allRegionsFinalEdges[key][0]), float(allRegionsFinalEdges[key][1]), minEdge, maxEdge, binWidth, 10e-10, 20.0, 0.0,  0.8, njet, name=key, variable="SigFrac%s"%(translator[key]["D"]))
                    plotter[hist_key].plot_inverseSignificance_vsClosureErr(finalSignificance, finalClosureErr, significances, closureErrs, edges, allRegionsFinalEdges[key], njet, name=key)

                if hist_key != Sig:
                    plotter[hist_key].plot_Var_vsDisc1Disc2(closureErrs[:,0], edges, float(allRegionsFinalEdges[key][0]), float(allRegionsFinalEdges[key][1]), minEdge, maxEdge, binWidth, 10e-10, 20.0, 0.0,  0.3, njet, name=key, variable="ClosureErr"   )
                    plotter[hist_key].plot_Var_vsDisc1Disc2(closureErrs[:,1], edges, float(allRegionsFinalEdges[key][0]), float(allRegionsFinalEdges[key][1]), minEdge, maxEdge, binWidth, 10e-10, 20.0, 0.0,  0.3, njet, name=key, variable="ClosureErrUnc")
                    plotter[hist_key].plot_Var_vsDisc1Disc2(pull[:,0],        edges, float(allRegionsFinalEdges[key][0]), float(allRegionsFinalEdges[key][1]), minEdge, maxEdge, binWidth, -20.0,  20.0, -5.0, 5.0, njet, name=key, variable="Pull"         )
                    plotter[hist_key].plot_Var_vsDisc1Disc2(pull[:,1],        edges, float(allRegionsFinalEdges[key][0]), float(allRegionsFinalEdges[key][1]), minEdge, maxEdge, binWidth, -20.0,  20.0, -5.0, 5.0, njet, name=key, variable="PullUnc"      )


                EventsPerNjets[hist_key][njet] = allRegionsFinalEvents[hist_key]
            
            edgesPerNjets[njet] = allRegionsFinalEdges

        # ------------------------------------
        # plot Disc1s vs Disc2s with all edges
        # ------------------------------------
        plotter["TT"].plot_Disc1VsDisc2(hist_lists[Sig],  allRegionsFinalEdges, njet, tag = "RPV550", name = "Val_BD_CD",   col1="yellow", col2="lime" )
        plotter["TT"].plot_Disc1VsDisc2(hist_lists["TT"], allRegionsFinalEdges, njet, tag = "TT",     name = "Val_BD_CD",   col1="yellow", col2="lime" )
        plotter["TT"].plot_Disc1VsDisc2(hist_lists[Sig],  allRegionsFinalEdges, njet, tag = "RPV550", name = "Val_SubDivD", col1="white",  col2="white")
        plotter["TT"].plot_Disc1VsDisc2(hist_lists["TT"], allRegionsFinalEdges, njet, tag = "TT",     name = "Val_SubDivD", col1="white",  col2="white")

        # put all latest bin edges to txt file
        binEdgesTable.writeLine(njet=njet, finalDiscs=allRegionsFinalEdges)

        # get the nEvents for each ABCD region
        abcdEventsTable.writeLine(njet=njet, finalSigFracs=allRegionsFinalSigFracs_TT["ABCD"], nEvents_AC=allRegionsFinalEvents["TT"]["ABCD"]["A"][0]+allRegionsFinalEvents["TT"]["ABCD"]["C"][0], nEvents_AB=allRegionsFinalEvents["TT"]["ABCD"]["A"][0]+allRegionsFinalEvents["TT"]["ABCD"]["B"][0])

        # get the nEvents for each B'D'EF region 
        bdefEventsTable.writeLine(njet=njet, finalSigFracs=allRegionsFinalSigFracs_TT["Val_bdEF"], nEvents_AC=allRegionsFinalEvents["TT"]["Val_bdEF"]["A"][0]+allRegionsFinalEvents["TT"]["Val_bdEF"]["C"][0], nEvents_AB=allRegionsFinalEvents["TT"]["Val_bdEF"]["A"][0]+allRegionsFinalEvents["TT"]["Val_bdEF"]["B"][0])
        fixedBDEF_EventsTable.writeLine(njet=njet, finalSigFracs=allRegionsFinalSigFracs_TT["Val_fixedBDEF"], nEvents_AC=allRegionsFinalEvents["TT"]["Val_fixedBDEF"]["A"][0]+allRegionsFinalEvents["TT"]["Val_fixedBDEF"]["C"][0], nEvents_AB=allRegionsFinalEvents["TT"]["Val_fixedBDEF"]["A"][0]+allRegionsFinalEvents["TT"]["Val_fixedBDEF"]["B"][0])

        # get the nEvents for each C'D'GH region 
        cdghEventsTable.writeLine(njet=njet, finalSigFracs=allRegionsFinalSigFracs_TT["Val_cdiGH"], nEvents_AC=allRegionsFinalEvents["TT"]["Val_cdiGH"]["A"][0]+allRegionsFinalEvents["TT"]["Val_cdiGH"]["C"][0], nEvents_AB=allRegionsFinalEvents["TT"]["Val_cdiGH"]["A"][0]+allRegionsFinalEvents["TT"]["Val_cdiGH"]["B"][0])
        fixedCDGH_EventsTable.writeLine(njet=njet, finalSigFracs=allRegionsFinalSigFracs_TT["Val_fixedCDGH"], nEvents_AC=allRegionsFinalEvents["TT"]["Val_fixedCDGH"]["A"][0]+allRegionsFinalEvents["TT"]["Val_fixedCDGH"]["C"][0], nEvents_AB=allRegionsFinalEvents["TT"]["Val_fixedCDGH"]["A"][0]+allRegionsFinalEvents["TT"]["Val_fixedCDGH"]["B"][0])

        # get the table for Validation region sub-division D
        subdBinEdgesTable.writeLine(njet=njet, finalDiscs=allRegionsFinalEdges["Val_subDivD"], finalSigFracs=allRegionsFinalSigFracs_TT["Val_subDivD"], final_nTot_Sig=allRegionsFinalEvents[Sig]["Val_subDivD"], final_nTot_Bkg=allRegionsFinalEvents["TT"]["Val_subDivD"]) 

        # get the table for bkg and sig+bkg events in ABCD regions
        bkgEventsTable.writeLine(njet=njet, finalBkgFracs=allRegionsFinalTTFracs["ABCD"])

    # ----------------------
    # make all closure plots
    # ----------------------
    for key, region in regions.items():

        # usual closure
        plotter["TT"].make_allClosures(edgesPerNjets,    EventsPerNjets["TT"], None,                    None,                None,                   njets, name = key, closureTag = "b",  bkgTag = "TT")
        plotter["TT"].make_allClosures(edgesPerNjets,    EventsPerNjets["TT"], None,                    EventsPerNjets[Sig], None,                   njets, name = key, closureTag = "sb", bkgTag = "TT")
        plotter["NonTT"].make_allClosures(edgesPerNjets, None,                 EventsPerNjets["NonTT"], None,                None,                   njets, name = key, closureTag = "b",  bkgTag = "NonTT")
        plotter["TT"].make_allClosures(edgesPerNjets,    EventsPerNjets["TT"], EventsPerNjets["NonTT"], None,                EventsPerNjets["Data"], njets, name = key, closureTag = "b",  bkgTag = "AllBkg")

        # data-MC closure
        if key != "ABCD":
            plotter["Data"].make_allClosures(edgesPerNjets, EventsPerNjets["TT"], None, None, EventsPerNjets["Data"], njets, name = key, closureTag = "b", bkgTag = "forTT")

    # -------------------------------------
    # add all edges to DoubleDisCo cfg file
    # -------------------------------------
    addEdges = addEdges_toDoubleDisco(args.year, args.sig, args.mass, args.channel, regions)
    addEdges.addEdges_toDoubleDiscoCfg(edgesPerNjets, njets)

    # ------------------
    # make all tex files
    # ------------------
    # put all latest bin edges to txt file
    binEdgesTable.writeClose()
    
    # get the nEvents for each ABCD region
    abcdEventsTable.writeClose()

    # get the nEvents for each B'D'EF region 
    bdefEventsTable.writeClose()
    fixedBDEF_EventsTable.writeClose()

    # get the nEvents for each C'D'GH region 
    cdghEventsTable.writeClose()
    fixedCDGH_EventsTable.writeClose()    

    # get the table for Validation region sub-division D
    subdBinEdgesTable.writeClose()    

    # get the table for bkg and sig+bkg events in ABCD regions
    bkgEventsTable.writeClose()

if __name__ == '__main__':
    main()
