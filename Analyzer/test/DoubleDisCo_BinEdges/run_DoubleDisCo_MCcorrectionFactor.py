import ROOT
import os
import argparse
import numpy as np
from collections import defaultdict

from DoubleDisCo_Regions     import *
from DoubleDisCo_Plotter     import *
from DoubleDisCo_Variances   import *
from DoubleDisCo_TableWriter import *

# The Aggregator class is fed Regions objects, extracts useful information from them
# and stores everything in a single dictionary, allowing information to be retreived
# in bulk
class Aggregator:

    # Initialize a master dictionary e.g. "data", to hold lots of useful things 
    def __init__(self, samples, njets, regions, boundaries):

        self.data       = {}
        self.samples    = samples
        self.njets      = njets
        self.boundaries = boundaries
        self.regions    = regions
        self.subregions = ["A", "B", "C", "D"]

    # Construct a specifically formatted string that is used
    # for accessing information from "data"
    def makeKey(self, variable = None, **kwargs):

        theKey = "" 

        if "region"   in kwargs: theKey += kwargs["region"]   + "_"
        if "njet"     in kwargs: theKey += kwargs["njet"]     + "_"
        if "boundary" in kwargs: theKey += "%.2f"%(kwargs["boundary"]) + "_"
        if "sample"   in kwargs: theKey += kwargs["sample"]   + "_"

        if "variable" != None:   theKey += variable           + "_"
  
        return theKey[:-1]
        
    # Main method to eat an instance of the regions class
    # and grab all the things we want from it for safe keeping
    def aggregate(self, regionObj, **kwargs):

        # get all final bin edges
        self.data[self.makeKey(variable = "edges",      **kwargs)] = regionObj.get("edges",      None, None, "TT")
        self.data[self.makeKey(variable = "finalEdges", **kwargs)] = regionObj.getFinal("edges",             "TT")

        Sig = None
        for sample in self.samples:
        
            if Sig != None:
                break

            for s in ["RPV", "SYY", "SHH"]:
                if s in sample:
                    Sig = sample
                    break

        # quantities and variables with any combination of bin edges
        for subregion in self.subregions:

            self.data[self.makeKey(variable = "sigFractions", **kwargs)] = {subregion : regionObj.get("sigFraction%s"%(subregion), None, None, Sig)  for subregion in self.subregions}
            self.data[self.makeKey(variable = "ttFractions",  **kwargs)] = {subregion : regionObj.get("ttFraction%s"%(subregion),  None, None, "TT") for subregion in self.subregions}

            # quantities and variables with the final choice of bin edges
            self.data[self.makeKey(variable = "sigFraction", **kwargs)] = {subregion : regionObj.getFinal("sigFractionA", Sig) for subregion in self.subregions}
            self.data[self.makeKey(variable = "ttFraction",  **kwargs)] = {subregion : regionObj.getFinal("ttFractionA", "TT") for subregion in self.subregions}

        # -----------------------------------------------
        # loop over for getting plots for TT, NonTT, Data
        # -----------------------------------------------
        for sample in self.samples:

            self.data[self.makeKey(variable = "nEvents",  sample = sample, **kwargs)] = {subregion : regionObj.getFinal("nEvents%s"%(subregion), sample) for subregion in self.subregions}
            self.data[self.makeKey(variable = "nEventsA", sample = sample, **kwargs)] = regionObj.get("nEventsA", None, None, sample)

            if not any(s in sample for s in ["RPV", "SYY", "SHH"]):
                self.data[self.makeKey(variable = "nonClosures", sample = sample, **kwargs)] = regionObj.get("closureError",      None, None, sample) # vars with any combination of bin edges
                self.data[self.makeKey(variable = "pulls",        sample = sample, **kwargs)] = regionObj.get("pull",              None, None, sample)
                self.data[self.makeKey(variable = "nonClosure",   sample = sample, **kwargs)] = regionObj.getFinal("closureError",             sample) # vars with the final choice of bin edges
                self.data[self.makeKey(variable = "pull",         sample = sample, **kwargs)] = regionObj.getFinal("pull",                     sample)

                self.data[self.makeKey(variable = "mcCorrection", sample = sample, **kwargs)] = regionObj.getFinal("mcCorrection",             sample) # vars with the final choice of bin edges

        self.data[self.makeKey(variable = "significances", **kwargs)] = regionObj.get("significance",      None, None, "TT")
        self.data[self.makeKey(variable = "significance",  **kwargs)] = regionObj.getFinal("significance",             "TT")

    def getPerNjets(self, variable, **kwargs):

        payload = {}
        for njet in self.njets:

            masterKey = self.makeKey(variable = variable, njet = njet, **kwargs)
            payload[njet] = self.data[masterKey]
              
        return payload

    def getPerBoundary(self, variable, **kwargs):
        
        payload = {}
        for boundary in self.boundaries:

            masterKey = self.makeKey(variable = variable, boundary = boundary, **kwargs)

            if masterKey not in self.data:
                print("Skipping key \"%s\""%(masterKey))
                continue

            payload[boundary] = self.data[masterKey]

        return payload

    def get(self, variable, **kwargs):
       
        masterKey = self.makeKey(variable = variable, **kwargs)

        return self.data[masterKey]

def main():

    # -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # command to run this script
    #   -- python run_DoubleDisCo_BinEdges.py --year 2016 --tt TT --nontt NonTT --ttVar TT_erdOn --data Data --sig RPV --mass 550 --channel 0l --metric New  --fixedDisc1edge 0.6 --fixedDisc2edge 0.6
    #   -- python run_DoubleDisCo_BinEdges.py --year 2016 --tt TT --nontt NonTT --ttVar TT_erdOn --data Data --sig RPV --mass 550 --channel 1l --metric New  --fixedDisc1edge 0.6 --fixedDisc2edge 0.6
    #   -- python run_DoubleDisCo_BinEdges.py --year 2016 --tt TT --nontt NonTT --ttVar TT_erdOn --data Data --sig RPV --mass 550 --channel 0l --metric New  
    #   -- python run_DoubleDisCo_BinEdges.py --year 2016 --tt TT --nontt NonTT --ttVar TT_erdOn --data Data --sig RPV --mass 550 --channel 1l --metric New  
    # --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
    usage  = "usage: %prog [options]"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("--year",           dest="year",           help="which year",            required=True)
    parser.add_argument("--path",           dest="path",           help="Input dir with histos", default="/uscms_data/d3/jhiltb/PO_Boxes/shared/2016_DisCo_0L_Cand1_1L") # Both 0L & 1l with OldSeed
    #parser.add_argument("--path",           dest="path",           help="Input dir with histos", default="/uscms_data/d3/jhiltb/PO_Boxes/shared/2016_DisCo_0L_Cand1_TopSeed_1L") # 0L with TopSeed
    parser.add_argument("--tt",             dest="tt",             help="TT",                    required=True)
    parser.add_argument("--nontt",          dest="nontt",          help="NonTT",                 required=True)
    parser.add_argument("--ttVar",          dest="ttVar",          help="TT Variances",          required=True)
    parser.add_argument("--sig",            dest="sig",            help="RPV, SYY",              default="RPV")
    parser.add_argument("--mass",           dest="mass",           help="signal mass",           default="550")
    parser.add_argument("--data",           dest="data",           help="JetHT, SingleLepton",   required=True)
    parser.add_argument("--channel",        dest="channel",        help="0l or 1l",              required=True)
    parser.add_argument("--metric",         dest="metric",         help="NN,New",                required=None, type=str)
    parser.add_argument("--fixedDisc1edge", dest="fixedDisc1edge", help="fixed d1 edge",         default=None, type=float)
    parser.add_argument("--fixedDisc2edge", dest="fixedDisc2edge", help="fixed d2 edge",         default=None, type=float)
    args = parser.parse_args()

    Sig     = "%s_%s"%(args.sig, args.mass)
    ttVar   = "%s"%(args.ttVar)
    samples = [args.tt, args.nontt, ttVar, Sig, args.data]

    # Make the output directories if they do not already exist
    plotsPath = {}; tablesPath = {}

    for sample in samples:

        # make directories to save plots and tables        
        if sample == Sig: continue

        if args.fixedDisc1edge != None or args.fixedDisc2edge != None:
            plotsPath[sample]  = "plots_MCcorrFactor_%s_%s_%s/%s_%s/%s/"%(args.fixedDisc1edge, args.fixedDisc2edge, sample, args.sig, args.mass, args.channel)
            
            if sample == "TT":
                tablesPath[sample] = "tables_MCcorrFactor_%s_%s/%s"%(args.fixedDisc1edge, sample, args.channel)


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
    njets = None
    histNames = "h_DoubleDisCo_disc1_disc2_%s_Njets"%(args.channel)
    njets = ["7"]#, "8"]#, "9", "10", "11", "12"]

    # ------------------
    # make all tex files
    # ------------------
    tag = ""
    if args.fixedDisc1edge != None or args.fixedDisc2edge != None:
        tag = "fixed"
    else:
        tag = "final"

    # make regionis list for adding all edges to DoubleDisCo cfg file
    regions = {"ABCD"          : "ABCD",
               "Val_fixedBDEF" : "fixedBDEF",
               "Val_fixedCDGH" : "fixedCDGH",
               "Val_subDivD"   : "subDivD",
    }

    # initialize the dictionaries of any regions
    translator = {"ABCD"          : {"A" : "A",  "B" : "B",  "C" : "C",  "D" : "D" },
                  "Val_fixedBDEF" : {"A" : "vA", "B" : "vB", "C" : "vC", "D" : "vD"},
                  "Val_fixedCDGH" : {"A" : "hA", "B" : "hB", "C" : "hC", "D" : "hD"},
                  "Val_subDivD"   : {"A" : "dA", "B" : "dB", "C" : "dC", "D" : "dD"},
    }

    # Will hold the pair of final edge values for the full ABCD region
    abcdFinalEdges = None

    # --------------------------------
    # Common calculations and plotters
    # --------------------------------
    plotter = {}; plotter_ttVar = {}

    for sample in samples:

        if sample == Sig: continue

        if sample == ttVar:
            plotter[sample]       = Common_Calculations_Plotters(plotsPath[sample], args.metric, args.year, args.sig, args.mass, args.channel)
            plotter_ttVar[sample] = ttVariances_Plotters(plotsPath[sample], args.metric, args.year, args.channel)

        else:
            plotter[sample] = Common_Calculations_Plotters(plotsPath[sample], args.metric, args.year, args.sig, args.mass, args.channel)

    regionGridWidth = 0.05

    # make the lists for rightBoundary and topBoundary
    list_boundaries = {"Val_fixedBDEF" : np.arange(0.4, 1.05, regionGridWidth),
                       "Val_fixedCDGH" : np.arange(0.4, 1.05, regionGridWidth),
                       "Val_subDivD"   : np.arange(0.6, 1.05, regionGridWidth)
    }

    # Make an Aggregator object to hold everything convenient
    theAggy = Aggregator(samples, njets, regions.keys(), list_boundaries["Val_fixedBDEF"])

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

                theEdgesClass = ABCDedges(hist_lists, Sig=Sig, ttVar=ttVar, fixedDisc1Edge=args.fixedDisc1edge, fixedDisc2Edge=args.fixedDisc2edge, metric=args.metric)
                theAggy.aggregate(theEdgesClass, region = key, njet = njet)

                abcdFinalEdges = theEdgesClass.getFinal("edges", "TT")

            # ----------------------------
            # make Sub-Division BD - Val I
            # ----------------------------
            elif key == "Val_fixedBDEF":

                for r in list_boundaries[key]:
                
                    theEdgesClass = bdEFedges(hist_lists, Sig=Sig, ttVar=ttVar, fixedDisc2Edge=abcdFinalEdges[1], rightBoundary=float(r), fixedDisc1Edge=float(r * 0.5), metrc=args.metric)
                    theAggy.aggregate(theEdgesClass, region = key, njet = njet, boundary = r)

            # -----------------------------
            # make Sub-Division CD - Val II
            # -----------------------------
            elif key == "Val_fixedCDGH":

                for t in list_boundaries[key]:
                
                    theEdgesClass = cdGHedges(hist_lists, Sig=Sig, ttVar=ttVar, fixedDisc1Edge=abcdFinalEdges[0], topBoundary=float(t), fixedDisc2Edge=float(t * 0.5), metric=args.metric)
                    theAggy.aggregate(theEdgesClass, region = key, njet = njet, boundary = t)

            # -----------------------------
            # make Sub-Division D - Val III
            # -----------------------------
            elif key == "Val_subDivD":

                for d in list_boundaries[key]:            

                    theEdgesClass = subDivDedges(hist_lists, Sig=Sig, ttVar=ttVar, rightBoundary=d, topBoundary=d, ABCDdisc1=float(d), ABCDdisc2=float(d), metric=args.metric)
                    theAggy.aggregate(theEdgesClass, region = key, njet = njet, boundary = d)


        # Now go back over the regions and make simple 2D plots
        for key, region in regions.items():

            if "Val_" not in key:
                continue

            for b in list_boundaries[key]:
                kwargs = {"region" : key, "njet" : njet, "boundary" : b}
  
                edges         = np.array(theAggy.get("edges", **kwargs), dtype=float)
                finalEdges    = theAggy.get("finalEdges",     **kwargs)
                closures      = theAggy.get("nonClosures",    sample = "TT", **kwargs)

                # ---------------------------
                # plot variable vs Disc1Disc2
                # ---------------------------
                plotter["TT"].plot_Var_vsDisc1Disc2(closures[:,0], edges, float(finalEdges[0]), float(finalEdges[1]), minEdge, maxEdge, binWidth, 10e-10, 20.0, 0.0,  0.5, njet, name=key+"_%.2f"%(b), variable="ClosureErr"   )
                plotter["TT"].plot_Var_vsDisc1Disc2(closures[:,1], edges, float(finalEdges[0]), float(finalEdges[1]), minEdge, maxEdge, binWidth, 10e-10, 20.0, 0.0,  0.5, njet, name=key+"_%.2f"%(b), variable="ClosureErrUnc")
    
    # ------------------------------
    # bkg subtracted data-MC closure
    # ------------------------------
    for key, region in regions.items():
        if "Val_" not in key:
            continue

        for b in list_boundaries[key]:
            eventsTT = {}; eventsData = {}; edgesPerNjets = {}
            for njet in njets:

                eventsTT.setdefault(njet, {}).setdefault(key, theAggy.get("nEvents", region = key, njet = njet, boundary = b, sample = "TT"))
                eventsData.setdefault(njet, {}).setdefault(key, theAggy.get("nEvents", region = key, njet = njet, boundary = b, sample = "Data"))
                edgesPerNjets.setdefault(njet, {}).setdefault(key, theAggy.get("finalEdges", region = key, njet = njet, boundary = b))

            plotter["Data"].make_allClosures(edgesPerNjets, eventsTT, None, None, eventsData, njets, name = key, closureTag = "b_%.2f"%(b), bkgTag = "forTT")

    # ----------------------------------------
    # Plot non-closure as function of boundary
    # ----------------------------------------
    for njet in njets:

        nonClosurePerBoundaryTT   = {}
        mcCorrectionPerBoundaryTT = {}

        nonClosurePerBoundaryData   = {}
        mcCorrectionPerBoundaryData = {}

        for key, region in regions.items():
            if "Val_" not in key:
                continue

            nonClosurePerBoundaryTT[key]   = theAggy.getPerBoundary(variable = "nonClosure",   sample = "TT", region = key, njet = njet)
            mcCorrectionPerBoundaryTT[key] = theAggy.getPerBoundary(variable = "mcCorrection", sample = "TT", region = key, njet = njet)

            nonClosurePerBoundaryData[key]   = theAggy.getPerBoundary(variable = "nonClosure",   sample = "Data", region = key, njet = njet)
            mcCorrectionPerBoundaryData[key] = theAggy.getPerBoundary(variable = "mcCorrection", sample = "Data", region = key, njet = njet)

        # Add function to plotter
        plotter["TT"].plot_VarVsBoundary(nonClosurePerBoundaryTT,   regionGridWidth/2.0, 0.0, 0.3, "Non-closure",   "NonClosureExt", njet)
        plotter["TT"].plot_VarVsBoundary(mcCorrectionPerBoundaryTT, regionGridWidth/2.0, 0.7, 1.3, "MC Correction", "MCcorrectionExt", njet)

        if njet == "7":
            plotter["Data"].plot_VarVsBoundary(nonClosurePerBoundaryData,   regionGridWidth/2.0, 0.0, 0.3, "Non-closure",   "NonClosureExt", njet)
            plotter["Data"].plot_VarVsBoundary(mcCorrectionPerBoundaryData, regionGridWidth/2.0, 0.7, 1.3, "MC Correction", "MCcorrectionExt", njet)

if __name__ == '__main__':
    main()
