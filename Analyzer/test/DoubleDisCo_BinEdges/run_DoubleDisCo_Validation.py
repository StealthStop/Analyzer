import os
import ROOT
import argparse

from DoubleDisCo_BinEdges                 import *
from DoubleDisCo_Optimized_BinEdges       import *
from DoubleDisCo_MCcorrectionFactor_TT    import *
from DoubleDisCo_MCcorrectionFactor_TTvar import *

def main():
    # Running for Optimized_BinEdges
    # --------------------------------------------------------------------------------------------------------------
    # command to run optimized bin edges:
    # to make a tex file in the "tables_Optimized_BinEdges_TT/", including top 5 edge pair
    # to print the optimized ABCD edges and use them as fixed edges on the command line
    #   -- python run_DoubleDisCo_Validation.py --run Optimized_BinEdges --year 2016 --channel 0l --numEdgeChoices 5
    #   -- python run_DoubleDisCo_Validation.py --run Optimized_BinEdges --year 2016 --channel 1l --numEdgeChoices 5
    # --------------------------------------------------------------------------------------------------------------

    # Running for BinEdges
    # --------------------------------------------------------------------------------------------------------------------------------
    # get the plots with optimized ABCD edges:
    #   -- python run_DoubleDisCo_Validation.py --run BinEdges --year 2016 --channel 0l --disc1edge 0.71 --disc2edge 0.65 --plotVars2D
    #   -- python run_DoubleDisCo_Validation.py --run BinEdges --year 2016 --channel 1l --disc1edge 0.57 --disc2edge 0.65 --plotVars2D
    # --------------------------------------------------------------------------------------------------------------------------------

    # Running for MCcorrectionFactor_TT
    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------
    # command to run this script
    #   -- python run_DoubleDisCo_Validation.py --run MCcorrectionFactor_TT --year 2016 --channel 0l --disc1edge 0.71 --disc2edge 0.65 --plotVarVsBoundary --fastMode
    #   -- python run_DoubleDisCo_Validation.py --run MCcorrectionFactor_TT --year 2016 --channel 1l --disc1edge 0.57 --disc2edge 0.65 --plotVarVsBoundary --fastMode
    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------  

    # -------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # command to run this script
    #   -- python run_DoubleDisCo_Validation.py --run MCcorrectionFactor_TTvar --year 2016 --channel 0l  --disc1edge 0.71 --disc2edge 0.65 --plotVarVsBoundary --fastMode  
    #   -- python run_DoubleDisCo_Validation.py --run MCcorrectionFactor_TTvar --year 2016 --channel 1l  --disc1edge 0.57 --disc2edge 0.65 --plotVarVsBoundary --fastMode
    # -------------------------------------------------------------------------------------------------------------------------------------------------------------------  

    usage  = "usage: %prog [options]"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("--run",               dest="run",               help="which code to run",                            required=True                                                         )
    parser.add_argument("--year",              dest="year",              help="which year",                                   required=True                                                         )
    #parser.add_argument("--path",              dest="path",              help="Input dir with histos",                        default="/uscms_data/d3/jhiltb/PO_Boxes/shared/2016_DisCo_0L_Cand1_1L") # with OldSeed - Old Ntuples
    parser.add_argument("--path",              dest="path",              help="Input dir with histos",                        default="/uscms/home/bcrossma/nobackup/analysis/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/hadd_DoubleDisCo_MassExclusion_Min3_Run2UL") # Run2UL 
    parser.add_argument("--outpath",           dest="outpath",           help="output dir where group dirs",                  default="Run2UL_MassExclusion_RPV"                                    )
    parser.add_argument("--tt",                dest="tt",                help="name of TT sample",                            default="TT"                                                          )
    parser.add_argument("--nontt",             dest="nontt",             help="name of NonTT sample",                         default="NonTT"                                                       )
    parser.add_argument("--ttVar",             dest="ttVar",             help="TT MCcorrectionFactor_TTvar (default no var)", default="TT"                                                          )
    parser.add_argument("--sig",               dest="sig",               help="signal model RPV, SYY",                        default="RPV"                                                         )
    parser.add_argument("--mass",              dest="mass",              help="signal mass",                                  default="550"                                                         )
    parser.add_argument("--data",              dest="data",              help="name of Data sample",                          default="Data"                                                        )
    parser.add_argument("--channel",           dest="channel",           help="0l or 1l",                                     required=True                                                         )
    parser.add_argument("--cmslabel",          dest="cmslabel",          help="CMS label: Preliminary/Work in Progress",      default="Work in Progress"                                            )
    parser.add_argument("--disc1edge",         dest="disc1edge",         help="fixed d1 edge",                                default=None,  type=float                                             )
    parser.add_argument("--disc2edge",         dest="disc2edge",         help="fixed d2 edge",                                default=None,  type=float                                             )
    parser.add_argument("--fastMode",          dest="fastMode",          help="Fast mode, don't scan all choices",            default=False, action="store_true"                                    ) 
    parser.add_argument("--njets",             dest="njets",             help="which njet bins to run on",     nargs="+",     default=["7", "8", "9", "10", "11", "12incl"], type=str               ) 
    parser.add_argument("--plotVars1D",        dest="plotVars1D",        help="Plot 1D var vs disc (slices)",                 default=False, action="store_true"                                    ) 
    parser.add_argument("--plotVars2D",        dest="plotVars2D",        help="Plot var vs disc1 and disc2 (2D)",             default=False, action="store_true"                                    ) 
    parser.add_argument("--plotDisc1VsDisc2",  dest="plotDisc1VsDisc2",  help="Plot disc1 and disc2 (2D)",                    default=False, action="store_true"                                    ) 
    parser.add_argument("--plotVarVsBoundary", dest="plotVarVsBoundary", help="Plot var vs boundary",                         default=False, action="store_true"                                    ) 
    parser.add_argument("--numEdgeChoices",    dest="numEdgeChoices",    help="number of edge choices to print",              default=4,     type=int                                               )
    parser.add_argument("--makeDiscoCfg",      dest="makeDiscoCfg",      help="make cfg file including opt. ABCD edges",      default=True,  action="store_true"                                    )
    args = parser.parse_args()

    # -------------------------------------------
    # Construct signal name from process and mass
    # -------------------------------------------
    Sig = "%s_%s"%(args.sig, args.mass)

    # -------------------------------------------------------------------------------
    # Names of samples/processes/data whose 2D disc1 vs disc2 histos will be analyzed
    # -------------------------------------------------------------------------------
    samples = [args.tt, "QCD", "TTX", "BG_OTHER", args.ttVar, Sig, args.data]
    #samples = [args.tt, args.nontt, args.ttVar, Sig, args.data]

    # --------------------------------------------------------
    # Make the output directories if they do not already exist
    # --------------------------------------------------------
    plotsPath = {}; tablesPath = {}

    # --------------------------------
    # Common calculations and plotters
    # --------------------------------
    plotter = {}

    for sample in samples:

        # make directories to save plots and tables        
        if sample == Sig: continue

        if sample in ["TTX", "QCD", "BG_OTHER"]:
            sample = "NonTT"

        if args.disc1edge != None or args.disc2edge != None:
            plotsPath[sample]  = "%s/%s_plots_%s_%s_%s_%s/%s_%s/%s/"%(args.outpath, args.year, args.run, args.disc1edge, args.disc2edge, sample, args.sig, args.mass, args.channel)
            
            # Make all paths for saving plots and tables
            # if paths do not already exist
            if not os.path.exists(plotsPath[sample]):
                os.makedirs(plotsPath[sample])

            plotter[sample] = Common_Calculations_Plotters(plotsPath[sample], args.year, args.sig, args.mass, args.channel, args.cmslabel)

    edgeStub = ""
    if args.disc1edge != None and args.disc2edge != None:
        edgeStub += "_" + str(args.disc1edge) + "_" + str(args.disc2edge)

    # ------------------------------------------------------
    # Create a path to save LaTeX tables (only in TT folder)
    # ------------------------------------------------------
    tablesPath["TT"] = "%s/%s_tables_%s%s_TT/%s"%(args.outpath, args.year, args.run, edgeStub, args.channel)
    if not os.path.exists(tablesPath["TT"]):
        os.makedirs(tablesPath["TT"])

    # -----------------
    # get SYY histogram
    # -----------------
    modelLabel = None
    if (args.sig != "RPV"):
        modelLabel = "Stealth"
    else:
        modelLabel = ""


    modelDecay = "2t6j"
    if ("SHH" in args.sig):
        modelDecay = "2t4b"

    # ----------
    # root files
    # ----------
    files = {
        "TT"             : ROOT.TFile.Open(args.path + "/" + args.year + "_TT.root"            ), # TT.root file includes also fsrUp/Down, isrUp/Down JECup/down, JERup/down histograms anymore
        "TT_fsrDown"     : ROOT.TFile.Open(args.path + "/" + args.year + "_TT.root"            ),
        "TT_fsrUp"       : ROOT.TFile.Open(args.path + "/" + args.year + "_TT.root"            ),
        "TT_isrDown"     : ROOT.TFile.Open(args.path + "/" + args.year + "_TT.root"            ),
        "TT_isrUp"       : ROOT.TFile.Open(args.path + "/" + args.year + "_TT.root"            ),
        "TT_erdON"       : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_erdON.root"      ),
        "TT_hdampDOWN"   : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_hdampDOWN.root"  ),
        "TT_hdampUP"     : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_hdampUP.root"    ),
        "TT_TuneCP5down" : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_TuneCP5down.root"),
        "TT_TuneCP5up"   : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_TuneCP5up.root"  ),
        "TT_JECdown"     : ROOT.TFile.Open(args.path + "/" + args.year + "_TT.root"            ),
        "TT_JECup"       : ROOT.TFile.Open(args.path + "/" + args.year + "_TT.root"            ),
        "TT_JERdown"     : ROOT.TFile.Open(args.path + "/" + args.year + "_TT.root"            ),
        "TT_JERup"       : ROOT.TFile.Open(args.path + "/" + args.year + "_TT.root"            ),
        #"NonTT"          : ROOT.TFile.Open(args.path + "/" + args.year + "_Non_TT.root"        ),
        "TTX"            : ROOT.TFile.Open(args.path + "/" + args.year + "_TTX.root"           ),
        "BG_OTHER"       : ROOT.TFile.Open(args.path + "/" + args.year + "_BG_OTHER.root"      ),
        "QCD"            : ROOT.TFile.Open(args.path + "/" + args.year + "_QCD.root"           ),
        "Data"           : ROOT.TFile.Open(args.path + "/" + args.year + "_Data.root"          ),
        Sig              : ROOT.TFile.Open(args.path + "/" + args.year + "_%s%s_%s_mStop-%s.root"%(modelLabel, args.sig, modelDecay, args.mass)),
    }

    # ---------------------
    # get the 2D histograms
    # --------------------- 
    histName = "h_DoubleDisCo_%s_disc1_disc2_%s_Njets${NJET}_ABCD"%(args.sig,args.channel)


    # ---------------------------------------------------------------
    # make regionis list for adding all edges to DoubleDisCo cfg file
    # ---------------------------------------------------------------
    regions = ["ABCD",
               "Val_BD",
               "Val_CD",
               "Val_D", 
    ]

    # ------------------------------------------
    # initialize the dictionaries of any regions
    # ------------------------------------------
    translator = {"ABCD"   : {"A" : "A",  "B" : "B",  "C" : "C",  "D" : "D" },
                  "Val_BD" : {"A" : "b",  "B" : "E",  "C" : "d",  "D" : "F" },
                  "Val_CD" : {"A" : "c",  "B" : "di", "C" : "G",  "D" : "H" },
                  "Val_D"  : {"A" : "dA", "B" : "dB", "C" : "dC", "D" : "dD"},
    }

    theApp = None
    if  args.run == "Optimized_BinEdges":
        theApp = Optimized_BinEdges(args.year, args.channel, Sig, args.mass, args.ttVar, translator, args.makeDiscoCfg)
        theApp.run(disc1edge=None, disc2edge=None, fastMode=False, samples=samples, files=files, histName=histName, njets=args.njets, regions=regions, tablesPath=tablesPath, numEdgeChoices=args.numEdgeChoices)

    elif args.run == "BinEdges":
        theApp = BinEdges(args.year, args.channel, Sig, args.mass, args.ttVar, translator)
        theApp.run(args.disc1edge, args.disc2edge, args.fastMode, plotVars1D=args.plotVars1D, plotVars2D=args.plotVars2D, plotDisc1VsDisc2=args.plotDisc1VsDisc2, samples=samples, files=files, histName=histName, njets=args.njets, regions=regions, plotter=plotter, tablesPath=tablesPath)

    elif args.run == "MCcorrectionFactor_TT":
        theApp = MCcorrectionFactor_TT(args.year, args.channel, Sig, args.mass, args.ttVar, translator)
        theApp.run(args.disc1edge, args.disc2edge, args.fastMode, plotVars2D=args.plotVars2D, plotVarVsBoundary=args.plotVarVsBoundary, samples=samples, files=files, histName=histName, njets=args.njets, regions=regions, plotter=plotter, tablesPath=tablesPath)

    elif args.run == "MCcorrectionFactor_TTvar":
        theApp = MCcorrectionFactor_TTvar(args.year, args.channel, Sig, args.mass, translator, edges=edgeStub, outpath=args.outpath)
        theApp.run(args.disc1edge, args.disc2edge, args.fastMode, samples=samples, files=files, histName=histName, njets=args.njets, regions=regions, plotter=plotter, tablesPath=tablesPath, edges=edgeStub)

def BryansHack(files, channel, Sig, mass, histName, regions, translator, disc1, disc2):

    if "0l" in channel:
        njets = ["8", "9", "10", "11", "12", "13incl"]
    if "1l" in channel:
        njets = ["7", "8", "9", "10", "11", "12incl"]
    if "2l" in channel:
        njets = ["6", "7", "8", "9", "10", "11incl"]

    hack = True

    samples = ["TT", "QCD", "TTX", "BG_OTHER", "TT", Sig, "Data"]
    tablesPath = {"TT": "test"}
    plotter = None

    theApp = MCcorrectionFactor_TTvar("Run2UL", channel, Sig, mass, translator, edges="", hack=hack)
    output = theApp.run(disc1, disc2, True, samples=samples, files=files, histName=histName, njets=njets, regions=regions, plotter=plotter, tablesPath=tablesPath, edges="")
   
    return output

if __name__ == "__main__":
    main() 
