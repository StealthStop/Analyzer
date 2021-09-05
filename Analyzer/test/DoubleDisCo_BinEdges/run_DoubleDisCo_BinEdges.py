import ROOT
import os
import argparse

from DoubleDisCo_Regions     import *
from DoubleDisCo_Plotter     import *
from DoubleDisCo_TableWriter import *

def main():

    # ----------------------------------------------------------------------------------------------------------------------
    # command to run this script
    #   -- python run_DoubleDisCo_BinEdges.py --year 2016 --model RPV --mass 550 --channel {0l, 1l} --tag Fixed --fixedDisc1edge 0.6 --fixedDisc2edge 0.6
    #   -- python run_DoubleDisCo_BinEdges.py --year 2016 --model RPV --mass 550 --channel {0l, 1l} --tag New
    # -----------------------------------------------------------------------------------------------------------------------  
    usage  = "usage: %prog [options]"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("--year",           dest="year",           help="which year",            required=True)
    parser.add_argument("--path",           dest="path",           help="Input dir with histos", default="/uscms/home/jhiltb/nobackup/PO_Boxes/DoubleDisCo_Reg_0L_loose_1L_RPV_2016_20210818_Output/")
    parser.add_argument("--model",          dest="model",          help="signal model",          default="RPV")
    parser.add_argument("--mass",           dest="mass",           help="signal mass",           default="550")
    parser.add_argument("--channel",        dest="channel",        help="0l, 1l",                required=True)
    parser.add_argument("--fixedDisc1edge", dest="fixedDisc1edge", help="fixed d1 edge",         default=None, type=float)
    parser.add_argument("--fixedDisc2edge", dest="fixedDisc2edge", help="fixed d2 edge",         default=None, type=float)
    parser.add_argument("--tag",            dest="tag",            help="tag for naming",        required=True, type=str)
    args = parser.parse_args()

    # Make the output directories if they do not already exist
    tablesPath = "tables/%s"%(args.channel)
    plotsPath  = "plots_%s/%s_%s/%s/"%(args.tag, args.model, args.mass, args.channel)
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
    # put all latest bin edges to tex file
    binEdgesTable = BinEdgesTable(args.channel, args.year, "All_%sBinEdges"%(args.tag), args.model, args.mass)

    # get the nEvents for each ABCD region
    abcdEventsTable = ABCDeventsTable(args.channel, args.year, "nEvents_%s_ABCD"%(args.tag), args.model, args.mass)

    # get the nEvents for each B'D'EF region 
    bdefEventsTable = BDEFeventsTable(args.channel, args.year, "nEvents_%s_Val_bdEF"%(args.tag), args.model, args.mass)

    # get the nEvents for each C'D'GH region
    cdghEventsTable = CDGHeventsTable(args.channel, args.year, "nEvents_%s_Val_cdiGH"%(args.tag), args.model, args.mass)

    # get the table for Validation region sub-division D
    subdBinEdgesTable = SubDivDeventsTable(args.channel, args.year, "%s_FinalBinEdges_Val_subDivD"%(args.tag), args.model, args.mass)

    # --------------------------------
    # Common calculations and plotters
    # --------------------------------
    plotter = Common_Calculations_Plotters(args.tag, args.year, args.model, args.mass, args.channel)

    # Hold onto per Njet things so we can plot altogether after the Njets loop
    bkgEventsPerNjets = {njet : None for njet in njets}
    sigEventsPerNjets = {njet : None for njet in njets}
    edgesPerNjets     = {njet : None for njet in njets}

    # ---------------
    # loop over njets
    # --------------- 
    for njet in njets: 

        histBkg = files["TT"].Get(histNames + njet)
        histSig = files["%s%s"%(args.model,args.mass)].Get(histNames + njet)

        minEdge  = histSig.GetXaxis().GetBinLowEdge(1) 
        maxEdge  = histSig.GetXaxis().GetBinLowEdge(histBkg.GetNbinsX()+1)
        binWidth = histSig.GetXaxis().GetBinWidth(1)

        allRegionsEdges = {}
        allRegionsSigEvents = {"ABCD"         : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                               "Val_bdEF"     : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                               "Val_cdiGH"    : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                               "Val_SubDiv_D" : {"A" : {}, "B" : {}, "C" : {}, "D" : {}}
        }
        allRegionsBkgEvents = {"ABCD"         : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                               "Val_bdEF"     : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                               "Val_cdiGH"    : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                               "Val_SubDiv_D" : {"A" : {}, "B" : {}, "C" : {}, "D" : {}}
        }
        allRegionsFinalSigEvents = {"ABCD"         : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                                    "Val_bdEF"     : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                                    "Val_cdiGH"    : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                                    "Val_SubDiv_D" : {"A" : {}, "B" : {}, "C" : {}, "D" : {}}
        }
        allRegionsSigFracs = {"ABCD"         : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                              "Val_bdEF"     : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                              "Val_cdiGH"    : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                              "Val_SubDiv_D" : {"A" : {}, "B" : {}, "C" : {}, "D" : {}}
        }
        allRegionsFinalSigFracs = {"ABCD"         : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                                   "Val_bdEF"     : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                                   "Val_cdiGH"    : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                                   "Val_SubDiv_D" : {"A" : {}, "B" : {}, "C" : {}, "D" : {}}
        }
        allRegionsFinalBkgEvents = {"ABCD"         : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                                    "Val_bdEF"     : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                                    "Val_cdiGH"    : {"A" : {}, "B" : {}, "C" : {}, "D" : {}},
                                    "Val_SubDiv_D" : {"A" : {}, "B" : {}, "C" : {}, "D" : {}}
        }


        translator = {"ABCD"         : {"A" : "A",  "B" : "B",  "C" : "C",  "D" : "D"},
                      "Val_bdEF"     : {"A" : "b",  "B" : "E",  "C" : "d",  "D" : "F"},
                      "Val_cdiGH"    : {"A" : "c",  "B" : "di", "C" : "G",  "D" : "H"},
                      "Val_SubDiv_D" : {"A" : "dA", "B" : "dB", "C" : "dC", "D" : "dD"}}

        # Loop through the regions and make the set of plots for each
        # Make sure ABCD goes first so that the val regions can use its info
        for region in ["ABCD", "Val_bdEF", "Val_cdiGH", "Val_SubDiv_D"]:

            theEdgesClass = None

            if region == "ABCD":

                # ---------------
                # ABCD Bin Edges
                # ---------------
                theEdgesClass = ABCDedges(histBkg, histSig, fixedDisc1Edge=args.fixedDisc1edge, fixedDisc2Edge=args.fixedDisc2edge)

            elif region == "Val_bdEF":

                # ------------------
                # bdEF Regions
                # ------------------
                #   Validation Regions in BD:             
                #          *                              
                #      E   *    B'  |    A                
                #   _______*________|________             
                #          *        |                     
                #      F   *    D'  |    C                
                #          *                              
                # Here, the ABCD disc1 edge separating A,C from B',D' is the right boundary for the bdEF validation region
                # In addition, the ABCD disc2 edge is fixed with respect to the bdEF validation region
                theEdgesClass = bdEFedges(histBkg=histBkg, histSig=histSig, fixedDisc2Edge=allRegionsEdges["ABCD"][1], rightBoundary=allRegionsEdges["ABCD"][0], nBkgA=allRegionsFinalBkgEvents["ABCD"]["A"][0], nBkgC=allRegionsFinalBkgEvents["ABCD"]["C"][0])

            elif region == "Val_cdiGH":

                # ------------------
                # cdGH Regions
                # ------------------
                #   Validation Regions in CD: 
                #
                #          B  |  A   
                #       ______|______
                #             |
                #          D' |  C' 
                #       *************
                #             |
                #          H  |  G
                # Here, the ABCD disc2 edge separating A,B from D',C' is the top boundary for the cdiGH validation region
                # In addition, the ABCD disc1 edge is fixed with respect to the cdiGH validation region
                theEdgesClass = cdGHedges(histBkg=histBkg, histSig=histSig, fixedDisc1Edge=allRegionsEdges["ABCD"][0], topBoundary=allRegionsEdges["ABCD"][1], nBkgA=allRegionsFinalBkgEvents["ABCD"]["A"][0], nBkgB=allRegionsFinalBkgEvents["ABCD"]["B"][0])

            elif region == "Val_SubDiv_D":

                # -----------------------------------------------------
                # Validation Regions as SubDivisions of each D
                # -----------------------------------------------------
                #            |
                #        B   |   A
                #            |
                #     _______|_______
                #        *   |
                #      dB*dA | 
                #     *******|   C
                #      dD*dC |
                #        *   |
                theEdgesClass = subDivDedges(histBkg, histSig, rightBoundary=allRegionsEdges["ABCD"][0], topBoundary=allRegionsEdges["ABCD"][1], ABCDdisc1=float(allRegionsEdges["ABCD"][0]), ABCDdisc2=float(allRegionsEdges["ABCD"][1]))

            edges         = np.array(theEdgesClass.get("edges"), dtype=float)
            significances = theEdgesClass.get("significance")
            closureErrs   = theEdgesClass.get("closureError")

            allRegionsSigEvents[region] = {"A" : theEdgesClass.get("nSigEventsA")}
            allRegionsBkgEvents[region] = {"A" : theEdgesClass.get("nBkgEventsA")}

            finalSignificance = theEdgesClass.getFinal("significance")
            finalClosureErr   = theEdgesClass.getFinal("closureError")

            allRegionsEdges[region]    = theEdgesClass.getFinal("edges") 
            allRegionsSigFracs[region] = {"A" : theEdgesClass.get("sigFractionA"),
                                          "B" : theEdgesClass.get("sigFractionB"),
                                          "C" : theEdgesClass.get("sigFractionC"),
                                          "D" : theEdgesClass.get("sigFractionD")}
            allRegionsFinalSigEvents[region] = {"A" : theEdgesClass.getFinal("nSigEventsA"),
                                                "B" : theEdgesClass.getFinal("nSigEventsB"),
                                                "C" : theEdgesClass.getFinal("nSigEventsC"),
                                                "D" : theEdgesClass.getFinal("nSigEventsD")}
            allRegionsFinalSigFracs[region]  = {"A" : theEdgesClass.getFinal("sigFractionA"),
                                                "B" : theEdgesClass.getFinal("sigFractionB"),
                                                "C" : theEdgesClass.getFinal("sigFractionC"),
                                                "D" : theEdgesClass.getFinal("sigFractionD")}
            allRegionsFinalBkgEvents[region] = {"A" : theEdgesClass.getFinal("nBkgEventsA"),
                                                "B" : theEdgesClass.getFinal("nBkgEventsB"),
                                                "C" : theEdgesClass.getFinal("nBkgEventsC"), 
                                                "D" : theEdgesClass.getFinal("nBkgEventsD")}
    
            # plot variable vs disc as 1D 
            for disc in [1, 2]:
                plotter.plot_VarVsDisc(closureErrs,    edges, binWidth/2.0, 1.0,  "ABCD Closure", "Closure",      disc, njet, region = region)
                plotter.plot_VarVsDisc(significances,  edges, binWidth/2.0, 5.0,  "Significance", "Significance", disc, njet, region = region)

                plotter.plot_VarVsDisc(allRegionsSigEvents[region]["A"], edges, binWidth/2.0, -1.0, "Weighted Signal Events",     "wSigEvts", disc, njet, region = region)
                plotter.plot_VarVsDisc(allRegionsBkgEvents[region]["A"], edges, binWidth/2.0, -1.0, "Weighted Background Events", "wBkgEvts", disc, njet, region = region)

                for subregion in translator["ABCD"].keys():
                    plotter.plot_VarVsDisc(allRegionsSigFracs[region][subregion], edges, binWidth/2.0, 0.8, "Signal Contamination %s"%(translator[region][subregion]), "SigFracs%s"%(translator[region][subregion]), disc, njet, region = region)

            # plot 2Ds
            plotter.plot_Var_vsDisc1Disc2(significances[:,0],                   edges, float(allRegionsEdges[region][0]), float(allRegionsEdges[region][1]), minEdge, maxEdge, binWidth, 20.0, 5.0, njet, name="Significance", region = region)
            plotter.plot_Var_vsDisc1Disc2(closureErrs[:,0],                     edges, float(allRegionsEdges[region][0]), float(allRegionsEdges[region][1]), minEdge, maxEdge, binWidth, 20.0, 0.3, njet, name="ClosureErr",   region = region)
            plotter.plot_Var_vsDisc1Disc2(allRegionsSigFracs[region]["A"][:,0], edges, float(allRegionsEdges[region][0]), float(allRegionsEdges[region][1]), minEdge, maxEdge, binWidth, 20.0, 0.8, njet, name="SigFrac%s"%(translator[region]["A"]),     region = region)

            plotter.plot_inverseSignificance_vsClosureErr(finalSignificance, finalClosureErr, significances, closureErrs, edges, allRegionsEdges[region], njet, name=region)

        bkgEventsPerNjets[njet] = allRegionsFinalBkgEvents
        sigEventsPerNjets[njet] = allRegionsFinalSigEvents
        edgesPerNjets[njet]     = allRegionsEdges

        # ------------------------------------
        # plot Disc1s vs Disc2s with all edges
        # ------------------------------------
        plotter.plot_Disc1VsDisc2(histSig, allRegionsEdges, njet, tag = "RPV550", name = "Val_BD_CD",   col1="yellow", col2="lime")
        plotter.plot_Disc1VsDisc2(histBkg, allRegionsEdges, njet, tag = "TT",     name = "Val_BD_CD",   col1="yellow", col2="lime")
        plotter.plot_Disc1VsDisc2(histSig, allRegionsEdges, njet, tag = "RPV550", name = "Val_SubDiv_D", col1="white", col2="white")
        plotter.plot_Disc1VsDisc2(histBkg, allRegionsEdges, njet, tag = "TT",     name = "Val_SubDiv_D", col1="white", col2="white")

        # put all latest bin edges to txt file
        binEdgesTable.writeLine(njet=njet, finalDiscs=allRegionsEdges)

        # get the nEvents for each ABCD region
        abcdEventsTable.writeLine(njet=njet, finalSigFracs=allRegionsFinalSigFracs["ABCD"],      nEvents_AC=allRegionsFinalBkgEvents["ABCD"]["A"][0]+allRegionsFinalBkgEvents["ABCD"]["C"][0], nEvents_AB=allRegionsFinalBkgEvents["ABCD"]["A"][0]+allRegionsFinalBkgEvents["ABCD"]["B"][0])

        # get the nEvents for each B'D'EF region 
        bdefEventsTable.writeLine(njet=njet, finalSigFracs=allRegionsFinalSigFracs["Val_bdEF"],  nEvents_AC=allRegionsFinalBkgEvents["Val_bdEF"]["A"][0]+allRegionsFinalBkgEvents["Val_bdEF"]["C"][0], nEvents_AB=allRegionsFinalBkgEvents["Val_bdEF"]["A"][0]+allRegionsFinalBkgEvents["Val_bdEF"]["B"][0])

        # get the nEvents for each C'D'GH region 
        cdghEventsTable.writeLine(njet=njet, finalSigFracs=allRegionsFinalSigFracs["Val_cdiGH"], nEvents_AC=allRegionsFinalBkgEvents["Val_cdiGH"]["A"][0]+allRegionsFinalBkgEvents["Val_cdiGH"]["C"][0], nEvents_AB=allRegionsFinalBkgEvents["Val_cdiGH"]["A"][0]+allRegionsFinalBkgEvents["Val_cdiGH"]["B"][0])

        # get the table for Validation region sub-division D
        subdBinEdgesTable.writeLine(njet=njet, finalDiscs=allRegionsEdges["Val_SubDiv_D"], finalSigFracs=allRegionsFinalSigFracs["Val_SubDiv_D"], final_nTot_Sig=allRegionsFinalSigEvents["Val_SubDiv_D"], final_nTot_Bkg=allRegionsFinalBkgEvents["Val_SubDiv_D"]) 

    # ----------------------
    # make all closure plots
    # ----------------------
    plotter.make_allClosures(edgesPerNjets, bkgEventsPerNjets, sigEventsPerNjets, njets, region = "ABCD",         tag = "sb")
    plotter.make_allClosures(edgesPerNjets, bkgEventsPerNjets, sigEventsPerNjets, njets, region = "Val_bdEF",     tag = "sb")    
    plotter.make_allClosures(edgesPerNjets, bkgEventsPerNjets, sigEventsPerNjets, njets, region = "Val_cdiGH",    tag = "sb")
    plotter.make_allClosures(edgesPerNjets, bkgEventsPerNjets, sigEventsPerNjets, njets, region = "Val_SubDiv_D", tag = "sb")

    plotter.make_allClosures(edgesPerNjets, bkgEventsPerNjets, None, njets, region = "ABCD",         tag = "b")
    plotter.make_allClosures(edgesPerNjets, bkgEventsPerNjets, None, njets, region = "Val_bdEF",     tag = "b")    
    plotter.make_allClosures(edgesPerNjets, bkgEventsPerNjets, None, njets, region = "Val_cdiGH",    tag = "b")
    plotter.make_allClosures(edgesPerNjets, bkgEventsPerNjets, None, njets, region = "Val_SubDiv_D", tag = "b")

    # ------------------
    # make all tex files
    # ------------------
    # put all latest bin edges to txt file
    binEdgesTable.writeClose()
    
    # get the nEvents for each ABCD region
    abcdEventsTable.writeClose()

    # get the nEvents for each B'D'EF region 
    bdefEventsTable.writeClose()

    # get the nEvents for each C'D'GH region 
    cdghEventsTable.writeClose()

    # get the table for Validation region sub-division D
    subdBinEdgesTable.writeClose()    

if __name__ == '__main__':
    main()
