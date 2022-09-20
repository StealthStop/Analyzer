import os
import ROOT
import argparse

from DoubleDisCo_BinEdges                 import *
from DoubleDisCo_Optimized_BinEdges       import *
from DoubleDisCo_MCcorrectionFactor_TT    import *
from DoubleDisCo_MCcorrectionFactor_TTvar import *

# Running for Optimized_BinEdges
# -------------------------------------------------------------------------------------------
# command to run optimized bin edges:
#   -- python run_DoubleDisCo_Validation.py --run Optimized_BinEdges --year 2016 --channel 0l
#   -- python run_DoubleDisCo_Validation.py --run Optimized_BinEdges --year 2016 --channel 1l
# ------------------------------------------------------------------------------------------- 

# Running for MCcorrectionFactor_TT
# -----------------------------------------------------------------------------------------------------------------------------------------
# command to run this script
#   -- python run_DoubleDisCo_Validation.py --run MCcorrectionFactor_TT --year 2016 --channel 0l --disc1edge 0.6 --disc2edge 0.6 --fastMode
#   -- python run_DoubleDisCo_Validation.py --run MCcorrectionFactor_TT --year 2016 --channel 1l --disc1edge 0.6 --disc2edge 0.6 --fastMode
# -----------------------------------------------------------------------------------------------------------------------------------------  

# Running for MCcorrectionFactor_TTvar
# ------------------------------------------------------------------------------------------------------------------------------------
# command to run this script
#   -- python run_DoubleDisCo_Validation.py --run MCcorrectionFactor_TTvar --year 2016 --channel 0l  --disc1edge 0.69 --disc2edge 0.74
#   -- python run_DoubleDisCo_Validation.py --run MCcorrectionFactor_TTvar --year 2016 --channel 1l  --disc1edge 0.59 --disc2edge 0.70
# ------------------------------------------------------------------------------------------------------------------------------------  

# Running for BinEdges
# --------------------------------------------------------------------------------------------------------------------------------------------
# get the plots with fixed ABCD and Validation edges:
#   -- python run_DoubleDisCo_Validation.py --run BinEdges --year 2016 --channel 0l --disc1edge 0.6 --disc2edge 0.6
#   -- python run_DoubleDisCo_Validation.py --run BinEdges --year 2016 --channel 1l --disc1edge 0.6 --disc2edge 0.6

# get the plots with optimized ABCD edges:
#   -- first run the "run_DoubleDisCo_Optimized_BinEdges.py" to print the optimized ABCD edges and use them as fixed edges on the command line
#   -- python run_DoubleDisCo_Validation.py --run BinEdges --year 2016 --channel 0l --disc1edge 0.69 --disc2edge 0.74
#   -- python run_DoubleDisCo_Validation.py --run BinEdges --year 2016 --channel 1l --disc1edge 0.59 --disc2edge 0.70
# --------------------------------------------------------------------------------------------------------------------------------------------  

usage  = "usage: %prog [options]"
parser = argparse.ArgumentParser(usage)
parser.add_argument("--run",               dest="run",               help="which code to run",                                  required=True)
parser.add_argument("--year",              dest="year",              help="which year",                                         required=True)
parser.add_argument("--path",              dest="path",              help="Input dir with histos",             default="/uscms_data/d3/jhiltb/PO_Boxes/shared/2016_DisCo_0L_Cand1_1L") # Both 0L & 1l with OldSeed
parser.add_argument("--tt",                dest="tt",                help="name of TT sample",                 default="TT"                  )
parser.add_argument("--nontt",             dest="nontt",             help="name of NonTT sample",              default="NonTT"               )
parser.add_argument("--ttVar",             dest="ttVar",             help="TT MCcorrectionFactor_TTvar (default no var)",     default="TT"                  )
parser.add_argument("--sig",               dest="sig",               help="signal model RPV, SYY",             default="RPV"                 )
parser.add_argument("--mass",              dest="mass",              help="signal mass",                       default="550"                 )
parser.add_argument("--data",              dest="data",              help="name of Data sample",               default="Data"                )
parser.add_argument("--channel",           dest="channel",           help="0l or 1l",                                           required=True)
parser.add_argument("--disc1edge",         dest="disc1edge",         help="fixed d1 edge",                     default=None,    type=float)
parser.add_argument("--disc2edge",         dest="disc2edge",         help="fixed d2 edge",                     default=None,    type=float)
parser.add_argument("--fastMode",          dest="fastMode",          help="Fast mode, don't scan all choices", default=False,   action="store_true") 
parser.add_argument("--njets",             dest="njets",             help="which njet bins to run on",         nargs="+",       default=["7", "8", "9", "10", "11", "12incl"], type=str) 
parser.add_argument("--plotVars1D",        dest="plotVars1D",        help="Plot 1D var vs disc (slices)",      default=False,   action="store_true") 
parser.add_argument("--plotVars2D",        dest="plotVars2D",        help="Plot var vs disc1 and disc2 (2D)",  default=False,   action="store_true") 
parser.add_argument("--plotDisc1VsDisc2",  dest="plotDisc1VsDisc2",  help="Plot disc1 and disc2 (2D)",         default=False,   action="store_true") 
parser.add_argument("--plotVarVsBoundary", dest="plotVarVsBoundary", help="Plot var vs boundary",              default=False,   action="store_true") 
parser.add_argument("--numEdgeChoices",    dest="numEdgeChoices",    help="number of edge choices to print",   default=4,       type=int)
args = parser.parse_args()

# Construct signal name from process and mass
Sig     = "%s_%s"%(args.sig, args.mass)

# Names of samples/processes/data whose 2D disc1 vs disc2 histos will be analyzed
samples = [args.tt, args.nontt, args.ttVar, Sig, args.data]

# Make the output directories if they do not already exist
plotsPath = {}; tablesPath = {}
# --------------------------------
# Common calculations and plotters
# --------------------------------
plotter = {}

for sample in samples:

    # make directories to save plots and tables        
    if sample == Sig: continue

    if args.disc1edge != None or args.disc2edge != None:
        plotsPath[sample]  = "plots_%s_%s_%s_%s/%s_%s/%s/"%(args.run, args.disc1edge, args.disc2edge, sample, args.sig, args.mass, args.channel)
        
        # Make all paths for saving plots and tables
        # if paths do not already exist
        if not os.path.exists(plotsPath[sample]):
            os.makedirs(plotsPath[sample])

        plotter[sample] = Common_Calculations_Plotters(plotsPath[sample], args.year, args.sig, args.mass, args.channel)

edgeStub = ""
if args.disc1edge != None and args.disc2edge != None:
    edgeStub += "_" + str(args.disc1edge) + "_" + str(args.disc2edge)

# Create a path to save LaTeX tables (only in TT folder)
tablesPath["TT"] = "tables_%s%s_TT/%s"%(args.run, edgeStub, args.channel)
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
histName = "h_DoubleDisCo_disc1_disc2_%s_Njets"%(args.channel)

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

theApp = None
if  args.run == "Optimized_BinEdges":
    theApp = Optimized_BinEdges(args.year, args.channel, Sig, args.mass, args.ttVar, translator)
    theApp.run(disc1edge=None, disc2edge=None, fastMode=False, samples=samples, files=files, histName=histName, njets=args.njets, regions=regions, tablesPath=tablesPath, numEdgeChoices=args.numEdgeChoices)

elif args.run == "MCcorrectionFactor_TT":
    theApp = MCcorrectionFactor_TT(args.year, args.channel, Sig, args.mass, args.ttVar, translator)
    theApp.run(args.disc1edge, args.disc2edge, args.fastMode, plotVars2D=args.plotVars2D, plotVarVsBoundary=args.plotVarVsBoundary, samples=samples, files=files, histName=histName, njets=args.njets, regions=regions, plotter=plotter, tablesPath=tablesPath)

elif args.run == "MCcorrectionFactor_TTvar":
    theApp = MCcorrectionFactor_TTvar(args.year, args.channel, Sig, args.mass, translator)
    theApp.run(args.disc1edge, args.disc2edge, args.fastMode, samples=samples, files=files, histName=histName, njets=args.njets, regions=regions, plotter=plotter, tablesPath=tablesPath)

elif args.run == "BinEdges":
    theApp = BinEdges(args.year, args.channel, Sig, args.mass, args.ttVar, translator)
    theApp.run(args.disc1edge, args.disc2edge, args.fastMode, plotVars1D=args.plotVars1D, plotVars2D=args.plotVars2D, plotDisc1VsDisc2=args.plotDisc1VsDisc2, samples=samples, files=files, histName=histName, njets=args.njets, regions=regions, plotter=plotter, tablesPath=tablesPath)
