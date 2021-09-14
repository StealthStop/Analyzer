import ROOT
import os
import argparse

from DoubleDisCo_Regions     import *
from DoubleDisCo_Plotter     import *
from DoubleDisCo_TableWriter import *

def main():

    # ----------------------------------------------------------------------------------------------------------------------------------------------
    # command to run this script
    #   -- python run_DoubleDisCo_BinEdges.py --year 2016 --model RPV --mass 550 --channel 0l --metric NN  --fixedDisc1edge 0.6 --fixedDisc2edge 0.6
    #   -- python run_DoubleDisCo_BinEdges.py --year 2016 --model RPV --mass 550 --channel 1l --metric NN  --fixedDisc1edge 0.6 --fixedDisc2edge 0.6
    #   -- python run_DoubleDisCo_BinEdges.py --year 2016 --model RPV --mass 550 --channel 0l --metric New --fixedDisc1edge 0.6 --fixedDisc2edge 0.6
    #   -- python run_DoubleDisCo_BinEdges.py --year 2016 --model RPV --mass 550 --channel 1l --metric New --fixedDisc1edge 0.6 --fixedDisc2edge 0.6
    #   -- python run_DoubleDisCo_BinEdges.py --year 2016 --model RPV --mass 550 --channel 0l --metric NN  
    #   -- python run_DoubleDisCo_BinEdges.py --year 2016 --model RPV --mass 550 --channel 1l --metric NN  
    #   -- python run_DoubleDisCo_BinEdges.py --year 2016 --model RPV --mass 550 --channel 0l --metric New 
    #   -- python run_DoubleDisCo_BinEdges.py --year 2016 --model RPV --mass 550 --channel 1l --metric New 
    # ----------------------------------------------------------------------------------------------------------------------------------------------  
    usage  = "usage: %prog [options]"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("--year",           dest="year",           help="which year",            required=True)
    parser.add_argument("--path",           dest="path",           help="Input dir with histos", default="/uscms/home/jhiltb/nobackup/PO_Boxes/DoubleDisCo_Reg_0L_loose_1L_RPV_2016_20210818_Output/")
    parser.add_argument("--model",          dest="model",          help="signal model",          default="RPV")
    parser.add_argument("--mass",           dest="mass",           help="signal mass",           default="550")
    parser.add_argument("--channel",        dest="channel",        help="0l or 1l",              required=True)
    parser.add_argument("--metric",         dest="metric",         help="NN,New",                required=True, type=str)
    parser.add_argument("--fixedDisc1edge", dest="fixedDisc1edge", help="fixed d1 edge",         default=None, type=float)
    parser.add_argument("--fixedDisc2edge", dest="fixedDisc2edge", help="fixed d2 edge",         default=None, type=float)
    args = parser.parse_args()

    # Make the output directories if they do not already exist
    plotsPath  = ""
    tablesPath = ""
    
    if args.fixedDisc1edge != None or args.fixedDisc2edge != None:
        plotsPath  = "plots_fixedEdges_%s/%s_%s/%s/"%(args.fixedDisc1edge, args.model, args.mass, args.channel)
        tablesPath = "tables_%s/%s"%(args.fixedDisc1edge, args.channel)
    else:
        plotsPath  = "plots/%s_%s/%s/"%(args.model, args.mass, args.channel)
        tablesPath = "tables/%s"%(args.channel)

    if not os.path.exists(tablesPath):
        os.makedirs(tablesPath)

    if not os.path.exists(plotsPath):
        os.makedirs(plotsPath)

    modelDecay = "2t6j"
    if ("SHH" in args.model):
        modelDecay = "2t4b"

    files = {
        "TT"                          : ROOT.TFile.Open(args.path + "/" + args.year + "_TT.root"),
        "QCD"                         : ROOT.TFile.Open(args.path + "/" + args.year + "_QCD.root"),  
        "%s%s"%(args.model,args.mass) : ROOT.TFile.Open(args.path + "/" + args.year + "_%s_%s_mStop-%s.root"%(args.model,modelDecay,args.mass)),
    }

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
    # -----------------
    tag = ""
    if args.fixedDisc1edge != None or args.fixedDisc2edge != None:
        tag = "fixed"
    else:
        tag = "final"
    # put all latest bin edges to tex file
    binEdgesTable = BinEdgesTable(tablesPath, args.channel, args.year, "All_%sBinEdges"%(tag), args.model, args.mass)

    # get the nEvents for each ABCD region
    abcdEventsTable = ABCDeventsTable(tablesPath, args.channel, args.year, "nEvents_%s_ABCD"%(tag), args.model, args.mass)

    # get the nEvents for each B'D'EF region 
    bdefEventsTable       = BDEFeventsTable(tablesPath, args.channel, args.year, "nEvents_%s_Val_bdEF"%(tag), args.model, args.mass)
    fixedBDEF_EventsTable = BDEFeventsTable(tablesPath, args.channel, args.year, "nEvents_%s_Val_fixedBDEF"%(tag), args.model, args.mass)

    # get the nEvents for each C'D'GH region
    cdghEventsTable       = CDGHeventsTable(tablesPath, args.channel, args.year, "nEvents_%s_Val_cdiGH"%(tag), args.model, args.mass)
    fixedCDGH_EventsTable = CDGHeventsTable(tablesPath, args.channel, args.year, "nEvents_%s_Val_fixedCDGH"%(tag), args.model, args.mass)

    # get the table for Validation region sub-division D
    subdBinEdgesTable = SubDivDeventsTable(tablesPath, args.channel, args.year, "%s_FinalBinEdges_Val_subDivD"%(tag), args.model, args.mass)

    # --------------------------------
    # Common calculations and plotters
    # --------------------------------
    plotter = Common_Calculations_Plotters(plotsPath, args.metric, args.year, args.model, args.mass, args.channel)

    # Hold onto per Njet things so we can plot altogether after the Njets loop
    bkgEventsPerNjets = {njet : None for njet in njets}
    sigEventsPerNjets = {njet : None for njet in njets}
    edgesPerNjets     = {njet : None for njet in njets}

    # make regionis list for adding all edges to DoubleDisCo cfg file
    regions = {"ABCD"          : "ABCD", 
               "Val_bdEF"      : "bdEF",
               "Val_fixedBDEF" : "fixedBDEF", 
               "Val_cdiGH"     : "cdGH", 
               "Val_fixedCDGH" : "fixedCDGH",
               "Val_subDivD"   : "subDivD", 
    }

    # ---------------
    # loop over njets
    # --------------- 
    for njet in njets: 

        histBkg = files["TT"].Get(histNames + njet)
        histSig = files["%s%s"%(args.model,args.mass)].Get(histNames + njet)

        minEdge  = histSig.GetXaxis().GetBinLowEdge(1) 
        maxEdge  = histSig.GetXaxis().GetBinLowEdge(histBkg.GetNbinsX()+1)
        binWidth = histSig.GetXaxis().GetBinWidth(1)

        # initialize the dictionaries of quantities and variables with any combination of bin edges
        allRegionsEdges = {}
        allRegionsSigEvents = {"ABCD"          : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                               "Val_bdEF"      : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                               "Val_cdiGH"     : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                               "Val_subDivD"   : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                               "Val_fixedBDEF" : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                               "Val_fixedCDGH" : {"A" : {}, "B" : {}, "C" : {}, "D" : {}}, 
        }

        allRegionsBkgEvents = {"ABCD"          : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                               "Val_bdEF"      : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                               "Val_cdiGH"     : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                               "Val_subDivD"   : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                               "Val_fixedBDEF" : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                               "Val_fixedCDGH" : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
        }
        
        allRegionsSigFracs = {"ABCD"          : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                              "Val_bdEF"      : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                              "Val_cdiGH"     : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                              "Val_subDivD"   : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                              "Val_fixedBDEF" : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                              "Val_fixedCDGH" : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
        }

        # initialize the dictionaries of quantities and variables with the final choice of bin edges
        allRegionsFinalSigEvents = {"ABCD"          : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                                    "Val_bdEF"      : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                                    "Val_cdiGH"     : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                                    "Val_subDivD"   : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                                    "Val_fixedBDEF" : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                                    "Val_fixedCDGH" : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
        }

        allRegionsFinalBkgEvents = {"ABCD"          : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                                    "Val_bdEF"      : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                                    "Val_cdiGH"     : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                                    "Val_subDivD"   : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                                    "Val_fixedBDEF" : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                                    "Val_fixedCDGH" : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
        }

        allRegionsFinalSigFracs = {"ABCD"          : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                                   "Val_bdEF"      : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                                   "Val_cdiGH"     : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                                   "Val_subDivD"   : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                                   "Val_fixedBDEF" : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                                   "Val_fixedCDGH" : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
        }

        # initialize the dictionaries of any regions
        translator = {"ABCD"          : {"A" : "A",  "B" : "B",  "C" : "C",  "D" : "D" },
                      "Val_bdEF"      : {"A" : "b",  "B" : "E",  "C" : "d",  "D" : "F" },
                      "Val_cdiGH"     : {"A" : "c",  "B" : "di", "C" : "G",  "D" : "H" },
                      "Val_subDivD"   : {"A" : "dA", "B" : "dB", "C" : "dC", "D" : "dD"},
                      "Val_fixedBDEF" : {"A" : "vA", "B" : "vB", "C" : "vC", "D" : "vD"},
                      "Val_fixedCDGH" : {"A" : "hA", "B" : "hB", "C" : "hC", "D" : "hD"},
        }

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

                theEdgesClass = ABCDedges(histBkg, histSig, fixedDisc1Edge=args.fixedDisc1edge, fixedDisc2Edge=args.fixedDisc2edge, metric=args.metric)

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

                theEdgesClass = bdEFedges(histBkg=histBkg, histSig=histSig, fixedDisc2Edge=allRegionsEdges["ABCD"][1], rightBoundary=allRegionsEdges["ABCD"][0], nBkgA=allRegionsFinalBkgEvents["ABCD"]["A"][0], nBkgC=allRegionsFinalBkgEvents["ABCD"]["C"][0], metric=args.metric)

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

                theEdgesClass = cdGHedges(histBkg=histBkg, histSig=histSig, fixedDisc1Edge=allRegionsEdges["ABCD"][0], topBoundary=allRegionsEdges["ABCD"][1], nBkgA=allRegionsFinalBkgEvents["ABCD"]["A"][0], nBkgB=allRegionsFinalBkgEvents["ABCD"]["B"][0], metric=args.metric)

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
                
                theEdgesClass = subDivDedges(histBkg, histSig, rightBoundary=allRegionsEdges["ABCD"][0], topBoundary=allRegionsEdges["ABCD"][1], ABCDdisc1=float(allRegionsEdges["ABCD"][0]), ABCDdisc2=float(allRegionsEdges["ABCD"][1]), metric=args.metric)

            # --------------------------
            #        *    ||   |                       
            #    vB  * vA ||   |    A                
            #  ______*____||___|________             
            #        *    ||   |                     
            #    vD  * vC ||   |    C                
            #        *                  
            # --------------------------             
            elif key == "Val_fixedBDEF":

               theEdgesClass = bdEFedges(histBkg=histBkg, histSig=histSig, fixedDisc2Edge=allRegionsEdges["ABCD"][1], rightBoundary=float(0.4), fixedDisc1Edge=float(0.2), metrc=args.metric)

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

                theEdgesClass = cdGHedges(histBkg=histBkg, histSig=histSig, fixedDisc1Edge=allRegionsEdges["ABCD"][0], topBoundary=float(0.4), fixedDisc2Edge=float(0.2), metric=args.metric)

            # quantities and variables with any combination of bin edges
            edges                    = np.array(theEdgesClass.get("edges"), dtype=float)
            significances            = theEdgesClass.get("significance")
            closureErrs              = theEdgesClass.get("closureError")
            allRegionsEdges[key]     = theEdgesClass.getFinal("edges")
            allRegionsSigEvents[key] = {"A" : theEdgesClass.get("nSigEventsA")}
            allRegionsBkgEvents[key] = {"A" : theEdgesClass.get("nBkgEventsA")}
            allRegionsSigFracs[key]  = {"A" : theEdgesClass.get("sigFractionA"),
                                        "B" : theEdgesClass.get("sigFractionB"),
                                        "C" : theEdgesClass.get("sigFractionC"),
                                        "D" : theEdgesClass.get("sigFractionD")}

            # quantities and variables with the final choice of bin edges
            finalSignificance             = theEdgesClass.getFinal("significance")
            finalClosureErr               = theEdgesClass.getFinal("closureError")
            allRegionsFinalSigEvents[key] = {"A" : theEdgesClass.getFinal("nSigEventsA"),
                                             "B" : theEdgesClass.getFinal("nSigEventsB"),
                                             "C" : theEdgesClass.getFinal("nSigEventsC"),
                                             "D" : theEdgesClass.getFinal("nSigEventsD")}
            allRegionsFinalBkgEvents[key] = {"A" : theEdgesClass.getFinal("nBkgEventsA"),
                                             "B" : theEdgesClass.getFinal("nBkgEventsB"),
                                             "C" : theEdgesClass.getFinal("nBkgEventsC"), 
                                             "D" : theEdgesClass.getFinal("nBkgEventsD")}
            allRegionsFinalSigFracs[key]  = {"A" : theEdgesClass.getFinal("sigFractionA"),
                                             "B" : theEdgesClass.getFinal("sigFractionB"),
                                             "C" : theEdgesClass.getFinal("sigFractionC"),
                                             "D" : theEdgesClass.getFinal("sigFractionD")}    

            # plot variable vs disc as 1D 
            for disc in [1, 2]:
                plotter.plot_VarVsDisc(closureErrs,   edges, binWidth/2.0, 1.0, "%s Closure"%(key),      "Closure",      disc, njet, name = key)
                plotter.plot_VarVsDisc(significances, edges, binWidth/2.0, 5.0, "%s Significance"%(key), "Significance", disc, njet, name = key)
                plotter.plot_VarVsDisc(allRegionsSigEvents[key]["A"], edges, binWidth/2.0, -1.0, "Weighted Signal Events",     "wSigEvts", disc, njet, name = key)
                plotter.plot_VarVsDisc(allRegionsBkgEvents[key]["A"], edges, binWidth/2.0, -1.0, "Weighted Background Events", "wBkgEvts", disc, njet, name = key)

                for subregion in translator["ABCD"].keys():
                    plotter.plot_VarVsDisc(allRegionsSigFracs[key][subregion], edges, binWidth/2.0, 0.8, "Signal Contamination %s"%(translator[key][subregion]), "SigFracs%s"%(translator[key][subregion]), disc, njet, name = key)

            # plot 2Ds
            plotter.plot_Var_vsDisc1Disc2(significances[:,0],                edges, float(allRegionsEdges[key][0]), float(allRegionsEdges[key][1]), minEdge, maxEdge, binWidth, 20.0, 5.0, njet, name=key, variable="Sign"                               )
            plotter.plot_Var_vsDisc1Disc2(closureErrs[:,0],                  edges, float(allRegionsEdges[key][0]), float(allRegionsEdges[key][1]), minEdge, maxEdge, binWidth, 20.0, 0.3, njet, name=key, variable="ClosureErr"                         )
            plotter.plot_Var_vsDisc1Disc2(closureErrs[:,1],                  edges, float(allRegionsEdges[key][0]), float(allRegionsEdges[key][1]), minEdge, maxEdge, binWidth, 20.0, 0.3, njet, name=key, variable="ClosureErrUnc"                      )
            plotter.plot_Var_vsDisc1Disc2(allRegionsSigFracs[key]["A"][:,0], edges, float(allRegionsEdges[key][0]), float(allRegionsEdges[key][1]), minEdge, maxEdge, binWidth, 20.0, 0.8, njet, name=key, variable="SigFrac%s"%(translator[key]["A"]))
            plotter.plot_Var_vsDisc1Disc2(allRegionsSigFracs[key]["B"][:,0], edges, float(allRegionsEdges[key][0]), float(allRegionsEdges[key][1]), minEdge, maxEdge, binWidth, 20.0, 0.8, njet, name=key, variable="SigFrac%s"%(translator[key]["B"]))
            plotter.plot_Var_vsDisc1Disc2(allRegionsSigFracs[key]["C"][:,0], edges, float(allRegionsEdges[key][0]), float(allRegionsEdges[key][1]), minEdge, maxEdge, binWidth, 20.0, 0.8, njet, name=key, variable="SigFrac%s"%(translator[key]["C"]))
            plotter.plot_Var_vsDisc1Disc2(allRegionsSigFracs[key]["D"][:,0], edges, float(allRegionsEdges[key][0]), float(allRegionsEdges[key][1]), minEdge, maxEdge, binWidth, 20.0, 0.8, njet, name=key, variable="SigFrac%s"%(translator[key]["D"]))
            plotter.plot_inverseSignificance_vsClosureErr(finalSignificance, finalClosureErr, significances, closureErrs, edges, allRegionsEdges[key], njet, name=key)

        bkgEventsPerNjets[njet] = allRegionsFinalBkgEvents
        sigEventsPerNjets[njet] = allRegionsFinalSigEvents
        edgesPerNjets[njet]     = allRegionsEdges

        # ------------------------------------
        # plot Disc1s vs Disc2s with all edges
        # ------------------------------------
        plotter.plot_Disc1VsDisc2(histSig, allRegionsEdges, njet, tag = "RPV550", name = "Val_BD_CD",   col1="yellow", col2="lime")
        plotter.plot_Disc1VsDisc2(histBkg, allRegionsEdges, njet, tag = "TT",     name = "Val_BD_CD",   col1="yellow", col2="lime")
        plotter.plot_Disc1VsDisc2(histSig, allRegionsEdges, njet, tag = "RPV550", name = "Val_SubDivD", col1="white",  col2="white")
        plotter.plot_Disc1VsDisc2(histBkg, allRegionsEdges, njet, tag = "TT",     name = "Val_SubDivD", col1="white",  col2="white")

        # put all latest bin edges to txt file
        binEdgesTable.writeLine(njet=njet, finalDiscs=allRegionsEdges)

        # get the nEvents for each ABCD region
        abcdEventsTable.writeLine(njet=njet, finalSigFracs=allRegionsFinalSigFracs["ABCD"], nEvents_AC=allRegionsFinalBkgEvents["ABCD"]["A"][0]+allRegionsFinalBkgEvents["ABCD"]["C"][0], nEvents_AB=allRegionsFinalBkgEvents["ABCD"]["A"][0]+allRegionsFinalBkgEvents["ABCD"]["B"][0])

        # get the nEvents for each B'D'EF region 
        bdefEventsTable.writeLine(njet=njet, finalSigFracs=allRegionsFinalSigFracs["Val_bdEF"], nEvents_AC=allRegionsFinalBkgEvents["Val_bdEF"]["A"][0]+allRegionsFinalBkgEvents["Val_bdEF"]["C"][0], nEvents_AB=allRegionsFinalBkgEvents["Val_bdEF"]["A"][0]+allRegionsFinalBkgEvents["Val_bdEF"]["B"][0])
        fixedBDEF_EventsTable.writeLine(njet=njet, finalSigFracs=allRegionsFinalSigFracs["Val_fixedBDEF"], nEvents_AC=allRegionsFinalBkgEvents["Val_fixedBDEF"]["A"][0]+allRegionsFinalBkgEvents["Val_fixedBDEF"]["C"][0], nEvents_AB=allRegionsFinalBkgEvents["Val_fixedBDEF"]["A"][0]+allRegionsFinalBkgEvents["Val_fixedBDEF"]["B"][0])

        # get the nEvents for each C'D'GH region 
        cdghEventsTable.writeLine(njet=njet, finalSigFracs=allRegionsFinalSigFracs["Val_cdiGH"], nEvents_AC=allRegionsFinalBkgEvents["Val_cdiGH"]["A"][0]+allRegionsFinalBkgEvents["Val_cdiGH"]["C"][0], nEvents_AB=allRegionsFinalBkgEvents["Val_cdiGH"]["A"][0]+allRegionsFinalBkgEvents["Val_cdiGH"]["B"][0])
        fixedCDGH_EventsTable.writeLine(njet=njet, finalSigFracs=allRegionsFinalSigFracs["Val_fixedCDGH"], nEvents_AC=allRegionsFinalBkgEvents["Val_fixedCDGH"]["A"][0]+allRegionsFinalBkgEvents["Val_fixedCDGH"]["C"][0], nEvents_AB=allRegionsFinalBkgEvents["Val_fixedCDGH"]["A"][0]+allRegionsFinalBkgEvents["Val_fixedCDGH"]["B"][0])

        # get the table for Validation region sub-division D
        subdBinEdgesTable.writeLine(njet=njet, finalDiscs=allRegionsEdges["Val_subDivD"], finalSigFracs=allRegionsFinalSigFracs["Val_subDivD"], final_nTot_Sig=allRegionsFinalSigEvents["Val_subDivD"], final_nTot_Bkg=allRegionsFinalBkgEvents["Val_subDivD"]) 

    # ----------------------
    # make all closure plots
    # ----------------------
    for key, region in regions.items():

        plotter.make_allClosures(edgesPerNjets, bkgEventsPerNjets, None, njets, name = key, closureTag = "b")
        plotter.make_allClosures(edgesPerNjets, bkgEventsPerNjets, sigEventsPerNjets, njets, name = key, closureTag = "sb") 

    # -------------------------------------
    # add all edges to DoubleDisCo cfg file
    # -------------------------------------
    addEdges = addEdges_toDoubleDisco(args.year, args.model, args.mass, args.channel, regions)
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

if __name__ == '__main__':
    main()
