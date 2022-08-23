import ROOT
import os
import argparse
import numpy as np

from collections import defaultdict
from DoubleDisCo_Regions import *
from DoubleDisCo_Plotter import *

# The Aggregator class is fed Regions objects, extracts useful information from them
# and stores everything in a single dictionary, allowing information to be retreived
# in bulk
class Aggregator:

    # -------------------------------------------------------------------------
    # Initialize a master dictionary e.g. "data", to hold lots of useful things
    # ------------------------------------------------------------------------- 
    def __init__(self, samples, njets, regions, boundaries):

        self.data       = {}
        self.samples    = samples + ["TTinData"]
        self.njets      = njets
        self.boundaries = boundaries
        self.regions    = regions
        self.subregions = ["A", "B", "C", "D"]

    # ------------------------------------------------------
    # Construct a specifically formatted string that is used
    # for accessing information from "data"
    # ------------------------------------------------------
    def makeKey(self, variable = None, **kwargs):

        theKey = "" 

        if "region"   in kwargs: theKey += kwargs["region"]   + "_"
        if "njet"     in kwargs: theKey += kwargs["njet"]     + "_"
        if "boundary" in kwargs: theKey += "%.2f"%(kwargs["boundary"]) + "_"
        if "sample"   in kwargs: theKey += kwargs["sample"]   + "_"

        if "variable" != None:   theKey += variable           + "_"
  
        return theKey[:-1]

    # --------------------------------------------------------    
    # Main method to eat an instance of the regions class
    # and grab all the things we want from it for safe keeping
    # --------------------------------------------------------
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

            self.data[self.makeKey(variable = "sigFractions%s"%(subregion), **kwargs)] = regionObj.get("sigFraction%s"%(subregion), None, None, Sig)
            self.data[self.makeKey(variable = "ttFractions%s"%(subregion),  **kwargs)] = regionObj.get("ttFraction%s"%(subregion),  None, None, "TT")

            # quantities and variables with the final choice of bin edges
            self.data[self.makeKey(variable = "sigFraction%s"%(subregion), **kwargs)] = regionObj.getFinal("sigFraction%s"%(subregion), Sig)
            self.data[self.makeKey(variable = "ttFraction%s"%(subregion),  **kwargs)] = regionObj.getFinal("ttFraction%s"%(subregion), "TT")

        # loop over for getting plots for TT, NonTT, Data
        for sample in self.samples:

            for subregion in self.subregions:
                self.data[self.makeKey(variable = "nEvents%s"%(subregion),  sample = sample, **kwargs)] = regionObj.getFinal("nEvents%s"%(subregion), sample)

            if not any(s in sample for s in ["RPV", "SYY", "SHH"]):
                self.data[self.makeKey(variable = "nonClosures",  sample = sample, **kwargs)] = regionObj.get("nonClosure",       None, None, sample) # vars with any combination of bin edges
                self.data[self.makeKey(variable = "pulls",        sample = sample, **kwargs)] = regionObj.get("pull",             None, None, sample)
                self.data[self.makeKey(variable = "closureCorrs", sample = sample, **kwargs)] = regionObj.get("closureCorr",      None, None, sample)
                self.data[self.makeKey(variable = "Closure",      sample = sample, **kwargs)] = regionObj.getFinal("Closure",                 sample)
                self.data[self.makeKey(variable = "nonClosure",   sample = sample, **kwargs)] = regionObj.getFinal("nonClosure",              sample) # vars with the final choice of bin edges
                self.data[self.makeKey(variable = "pull",         sample = sample, **kwargs)] = regionObj.getFinal("pull",                    sample)
                self.data[self.makeKey(variable = "closureCorr",  sample = sample, **kwargs)] = regionObj.getFinal("closureCorr",             sample) # vars with the final choice of bin edges

                if sample == "TTinData":
                    self.data[self.makeKey(variable = "closureCorrTTinDataVsTT",        sample = sample, **kwargs)] = regionObj.getFinal("closureCorrTTinDataVsTT",        sample)
                    self.data[self.makeKey(variable = "MC_corrected_dataClosure",       sample = sample, **kwargs)] = regionObj.getFinal("MC_corrected_dataClosure",       sample)
                    self.data[self.makeKey(variable = "closureCorrTTinDataVsTTvar",     sample = sample, **kwargs)] = regionObj.getFinal("closureCorrTTinDataVsTTvar",     sample)
                    self.data[self.makeKey(variable = "MC_ttVar_corrected_dataClosure", sample = sample, **kwargs)] = regionObj.getFinal("MC_ttVar_corrected_dataClosure", sample)

        self.data[self.makeKey(variable = "significances", **kwargs)] = regionObj.get("significance",      None, None, "TT")
        self.data[self.makeKey(variable = "significance",  **kwargs)] = regionObj.getFinal("significance",             "TT")

    # --------------------------------
    # Get variables for each Njets bin
    # --------------------------------
    def getPerNjets(self, variable, **kwargs):

        payload = {}
        for njet in self.njets:

            masterKey = self.makeKey(variable = variable, njet = njet, **kwargs)
            payload[njet] = self.data[masterKey]
              
        return payload

    # -------------------------------
    # Get variables for each boundary
    # -------------------------------
    def getPerBoundary(self, variable, **kwargs):
        
        payload = {}
        newKwargs = kwargs.copy()
        sample = newKwargs.pop("sample", None)

        for boundary in self.boundaries:

            masterKey = self.makeKey(variable = variable, boundary = boundary, **kwargs)
            chefKey   = self.makeKey(variable = "sigFractionA", boundary = boundary, **newKwargs)

            if masterKey not in self.data:
                print("Skipping key \"%s\""%(masterKey))
                continue

            # this statement for data and data/MC closure correction
            sigFracA = self.data[chefKey][0]
            if (sigFracA >= 0.05 and sample == "TTinData"):
                payload[boundary] = (-999.0, 0.0)

            else: 
                payload[boundary] = self.data[masterKey]

        return payload 

    # -------------
    # Get variables
    # -------------
    def get(self, variable, **kwargs):
       
        masterKey = self.makeKey(variable = variable, **kwargs)

        return self.data[masterKey]

def main():

    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # command to run this script
    #   -- python run_DoubleDisCo_MCcorrectionFactor.py --year 2016 --tt TT --nontt NonTT --ttVar TT_erdOn --data Data --sig RPV --mass 550 --channel 0l --fixedDisc1edge 0.69 --fixedDisc2edge 0.74 --fastMode   
    #   -- python run_DoubleDisCo_MCcorrectionFactor.py --year 2016 --tt TT --nontt NonTT --ttVar TT_erdOn --data Data --sig RPV --mass 550 --channel 1l --fixedDisc1edge 0.59 --fixedDisc2edge 0.70 --fastMode
    # ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------  
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
    parser.add_argument("--fixedDisc1edge", dest="fixedDisc1edge", help="fixed d1 edge",         default=None, type=float)
    parser.add_argument("--fixedDisc2edge", dest="fixedDisc2edge", help="fixed d2 edge",         default=None, type=float)
    parser.add_argument("--fastMode",       dest="fastMode",       help="Fast mode, don't scan all choices ", default=False,  action="store_true") 
    args = parser.parse_args()

    Sig     = "%s_%s"%(args.sig, args.mass)
    ttVar   = "%s"%(args.ttVar)
    samples = [args.tt, args.nontt, ttVar, Sig, args.data]

    # Make the output directories if they do not already exist
    plotsPath = {}

    for sample in samples:

        # make directories to save plots and tables        
        if sample == Sig: continue

        if args.fixedDisc1edge != None or args.fixedDisc2edge != None:
            plotsPath[sample]  = "plots_MCcorrFactor_%s_%s_%s/%s_%s/%s/"%(args.fixedDisc1edge, args.fixedDisc2edge, sample, args.sig, args.mass, args.channel)
            
        # 
        if not os.path.exists(plotsPath[sample]):
            os.makedirs(plotsPath[sample])

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
    njets = ["7", "8", "9", "10", "11", "12incl"]

    # make regionis list for adding all edges to DoubleDisCo cfg file
    regions = ["ABCD",
               "Val_BD",
               "Val_CD",
               "Val_D", 
    ]

    # initialize the dictionaries of any regions
    translator = {"ABCD"   : {"A" : "A",  "B" : "B",  "C" : "C",  "D" : "D" },
                  "Val_BD" : {"A" : "vA", "B" : "vB", "C" : "vC", "D" : "vD"},
                  "Val_CD" : {"A" : "hA", "B" : "hB", "C" : "hC", "D" : "hD"},
                  "Val_D"  : {"A" : "dA", "B" : "dB", "C" : "dC", "D" : "dD"},
    }

    # Will hold the pair of final edge values for the full ABCD region
    abcdFinalEdges = None

    # make colors for each validation regions
    color = { "Val_BD" : "#DDBB87",
              "Val_CD" : "#429c93", 
              "Val_D"  : "#990099",
    }

    # --------------------------------
    # Common calculations and plotters
    # --------------------------------
    plotter = {}; plotter_ttVar = {}

    for sample in samples:

        if sample == Sig: continue

        if sample == ttVar:
            plotter[sample] = Common_Calculations_Plotters(plotsPath[sample], args.year, args.sig, args.mass, args.channel)
        else:
            plotter[sample] = Common_Calculations_Plotters(plotsPath[sample], args.year, args.sig, args.mass, args.channel)

    regionGridWidth = 0.05

    # ------------------------------------------------
    # make the lists for rightBoundary and topBoundary
    # ------------------------------------------------
    list_boundaries = {"Val_BD" : np.arange(0.40, 1.05, regionGridWidth),
                       "Val_CD" : np.arange(0.40, 1.05, regionGridWidth),
                       "Val_D"  : np.arange(0.60, 1.05, regionGridWidth)
    }

    # -------------------------------------------------------
    # Make an Aggregator object to hold everything convenient
    # -------------------------------------------------------
    theAggy = Aggregator(samples, njets, regions, list_boundaries["Val_BD"])

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
        for region in regions:

            theEdgesClass = None

            # -----------------
            # make ABCD regions
            # -----------------
            if region == "ABCD":

                theEdgesClass = All_Regions(hist_lists, Sig=Sig, ttVar=ttVar, fixedDisc1Edge=args.fixedDisc1edge, fixedDisc2Edge=args.fixedDisc2edge, fastMode=args.fastMode)
                theAggy.aggregate(theEdgesClass, region = region, njet = njet)

                abcdFinalEdges = theEdgesClass.getFinal("edges", "TT")

            # ----------------------------
            # make Sub-Division BD - Val I
            # ----------------------------
            elif region == "Val_BD":

                for r in list_boundaries[region]:
               
                    disc1_edge = ((float(abcdFinalEdges[0]) - 0.2) / (1.0 - 0.4)) * (float(r) - 0.4) + 0.2  
                    theEdgesClass = All_Regions(hist_lists, Sig=Sig, ttVar=ttVar, fixedDisc2Edge=abcdFinalEdges[1], rightBoundary=float(r), fixedDisc1Edge=float(disc1_edge), fastMode=args.fastMode)
                    theAggy.aggregate(theEdgesClass, region = region, njet = njet, boundary = r)

            # -----------------------------
            # make Sub-Division CD - Val II
            # -----------------------------
            elif region == "Val_CD":

                for t in list_boundaries[region]:

                    disc2_edge = ((float(abcdFinalEdges[1]) - 0.2) / (1.0 - 0.4)) * (float(t) - 0.4) + 0.2
                    theEdgesClass = All_Regions(hist_lists, Sig=Sig, ttVar=ttVar, fixedDisc1Edge=abcdFinalEdges[0], topBoundary=float(t), fixedDisc2Edge=float(disc2_edge), fastMode=args.fastMode)
                    theAggy.aggregate(theEdgesClass, region = region, njet = njet, boundary = t)

            # -----------------------------
            # make Sub-Division D - Val III
            # -----------------------------
            elif region == "Val_D":

                for d in list_boundaries[region]:            

                    disc1 = (float(abcdFinalEdges[0]) - (float(abcdFinalEdges[0]) / 2.0)) / (1.0 - float(abcdFinalEdges[0])) * (float(d) - float(abcdFinalEdges[0])) + (float(abcdFinalEdges[0]) / 2.0)
                    disc2 = (float(abcdFinalEdges[1]) - (float(abcdFinalEdges[1]) / 2.0)) / (1.0 - float(abcdFinalEdges[1])) * (float(d) - float(abcdFinalEdges[1])) + (float(abcdFinalEdges[1]) / 2.0)
                    theEdgesClass = All_Regions(hist_lists, Sig=Sig, ttVar=ttVar, rightBoundary=d, topBoundary=d, fixedDisc1Edge=float(disc1), fixedDisc2Edge=float(disc2), fastMode=args.fastMode)
                    theAggy.aggregate(theEdgesClass, region = region, njet = njet, boundary = d)


        # ---------------------------
        # Plot variable vs Disc1Disc2
        # ---------------------------
        for region in regions:

            if "Val_" not in region:
                continue

            for b in list_boundaries[region]:
                kwargs = {"region" : region, "njet" : njet, "boundary" : b}
  
                edges                 = np.array(theAggy.get("edges", **kwargs), dtype=float)
                finalEdges            = theAggy.get("finalEdges",                        **kwargs)
                closures              = theAggy.get("nonClosures",  sample = "TT",       **kwargs)
                closureCorr_TT        = theAggy.get("closureCorrs", sample = "TT",       **kwargs)
                closureCorr_TTinData  = theAggy.get("closureCorrs", sample = "TTinData", **kwargs)

                #plotter["TT"].plot_Var_vsDisc1Disc2(closures[:,0],               edges, float(finalEdges[0]), float(finalEdges[1]), minEdge, maxEdge, binWidth, 10e-10, 20.0, 0.0,  0.5, njet, name=region+"_%.2f"%(b), variable="NonClosure"             )
                #plotter["TT"].plot_Var_vsDisc1Disc2(closures[:,1],               edges, float(finalEdges[0]), float(finalEdges[1]), minEdge, maxEdge, binWidth, 10e-10, 20.0, 0.0,  0.5, njet, name=region+"_%.2f"%(b), variable="NonClosureUnc"          )
                #plotter["TT"].plot_Var_vsDisc1Disc2(closureCorr_TT[:,0],         edges, float(finalEdges[0]), float(finalEdges[1]), minEdge, maxEdge, binWidth, 10e-10, 20.0, 0.7,  1.3, njet, name=region+"_%.2f"%(b), variable="closureCorr_TT"         )
                #plotter["TT"].plot_Var_vsDisc1Disc2(closureCorr_TT[:,1],         edges, float(finalEdges[0]), float(finalEdges[1]), minEdge, maxEdge, binWidth, 10e-10, 20.0, 0.7,  1.3, njet, name=region+"_%.2f"%(b), variable="closureCorrUnc_TT"      )
                #plotter["Data"].plot_Var_vsDisc1Disc2(closureCorr_TTinData[:,0], edges, float(finalEdges[0]), float(finalEdges[1]), minEdge, maxEdge, binWidth, 10e-10, 20.0, 0.7,  1.3, njet, name=region+"_%.2f"%(b), variable="closureCorr_TTinData"   )
                #plotter["Data"].plot_Var_vsDisc1Disc2(closureCorr_TTinData[:,1], edges, float(finalEdges[0]), float(finalEdges[1]), minEdge, maxEdge, binWidth, 10e-10, 20.0, 0.7,  1.3, njet, name=region+"_%.2f"%(b), variable="closureCorrUnc_TTinData")
 
    # -----------------------------------
    # Plot bkg subtracted data-MC closure
    # -----------------------------------
    for region in regions:
        if "Val_" not in region:
            continue

        for b in list_boundaries[region]:
            eventsTT = {}; eventsData = {}; edgesPerNjets = {}
            for njet in njets:

                nEventsTT   = {subregion : theAggy.get("nEvents%s"%(subregion),    region = region, njet = njet, boundary = b, sample = "TT"  ) for subregion in ["A", "B", "C", "D"]}
                nEventsData = {subregion : theAggy.get("nEvents%s"%(subregion),    region = region, njet = njet, boundary = b, sample = "Data") for subregion in ["A", "B", "C", "D"]}

                eventsTT.setdefault(njet,      {}).setdefault(region, nEventsTT)
                eventsData.setdefault(njet,    {}).setdefault(region, nEventsData)
                edgesPerNjets.setdefault(njet, {}).setdefault(region, theAggy.get("finalEdges", region = region, njet = njet, boundary = b                 ))

            plotter["Data"].make_allClosures(edgesPerNjets, eventsTT, None, None, eventsData, njets, name = region, closureTag = "b_%.2f"%(b), bkgTag = "forTT")

    # ----------------------------------------
    # Plot non-closure as function of boundary
    # ----------------------------------------
    for njet in njets:

        # TT
        weighted_nTTeventsA_PerBoundaryTT = {}; nonClosurePerBoundaryTT    = {}; 
        sigFractionA_PerBoundaryTT        = {}; sigFractionB_PerBoundaryTT = {}; sigFractionC_PerBoundaryTT = {}; sigFractionD_PerBoundaryTT = {}
        closureCorrPerBoundaryTT          = {}
        closureCorrPerBoundaryTTvar       = {}

        # TTinData = Data - NonTT
        weighted_nTTEventsA_PerBoundaryTTinData = {}; nonClosurePerBoundaryTTinData = {}; closurePerBoundaryTTinData = {}; MC_corrected_dataClosure_PerBoundaryTTinData = {} 
        MC_ttVar_corrected_dataClosure_PerBoundaryTTinData = {} 


        for region in regions:
            if "Val_" not in region:
                continue

            # TT
            weighted_nTTeventsA_PerBoundaryTT[region] = theAggy.getPerBoundary(variable = "nEventsA",    sample = "TT",  region = region, njet = njet)
            nonClosurePerBoundaryTT[region]           = theAggy.getPerBoundary(variable = "nonClosure",  sample = "TT",  region = region, njet = njet)
            sigFractionA_PerBoundaryTT[region]        = theAggy.getPerBoundary(variable = "sigFractionA",                region = region, njet = njet)
            sigFractionB_PerBoundaryTT[region]        = theAggy.getPerBoundary(variable = "sigFractionB",                region = region, njet = njet)
            sigFractionC_PerBoundaryTT[region]        = theAggy.getPerBoundary(variable = "sigFractionC",                region = region, njet = njet)
            sigFractionD_PerBoundaryTT[region]        = theAggy.getPerBoundary(variable = "sigFractionD",                region = region, njet = njet)
            closureCorrPerBoundaryTT[region]          = theAggy.getPerBoundary(variable = "closureCorr", sample = "TT",  region = region, njet = njet)            
            closureCorrPerBoundaryTTvar[region]       = theAggy.getPerBoundary(variable = "closureCorr", sample = ttVar, region = region, njet = njet)            

            # TTinData = Data - NonTT
            weighted_nTTEventsA_PerBoundaryTTinData[region]            = theAggy.getPerBoundary(variable = "nEventsA",                       sample = "TTinData", region = region, njet = njet)
            nonClosurePerBoundaryTTinData[region]                      = theAggy.getPerBoundary(variable = "nonClosure",                     sample = "TTinData", region = region, njet = njet)
            closurePerBoundaryTTinData[region]                         = theAggy.getPerBoundary(variable = "Closure",                        sample = "TTinData", region = region, njet = njet)
            MC_corrected_dataClosure_PerBoundaryTTinData[region]       = theAggy.getPerBoundary(variable = "MC_corrected_dataClosure",       sample = "TTinData", region = region, njet = njet)
            MC_ttVar_corrected_dataClosure_PerBoundaryTTinData[region] = theAggy.getPerBoundary(variable = "MC_ttVar_corrected_dataClosure", sample = "TTinData", region = region, njet = njet)

        # set y axis for higher njets bins
        yMin = None; yMax = None
        if njet > 8:
            yMin = 0.5; yMax = 1.5

        else:
            yMin = 0.7; yMax = 1.3

        # TT
        #plotter["TT"].plot_VarVsBoundary(weighted_nTTeventsA_PerBoundaryTT, regionGridWidth/2.0, None, None, None, "Weighted Events in A [MC]", "TT_weighted_nTTeventsA_PerBoundary", njet, color)
        #plotter["TT"].plot_VarVsBoundary(nonClosurePerBoundaryTT,           regionGridWidth/2.0, 0.0,  0.3,  None, "Non-Closure [MC]",          "TT_NonClosure_PerBoundary",          njet, color)
        #plotter["TT"].plot_VarVsBoundary(sigFractionA_PerBoundaryTT,        regionGridWidth/2.0, None, None, None, "SigFrac\'A\'",              "TT_SigFracA_PerBoundary",            njet, color) 
        #plotter["TT"].plot_VarVsBoundary(sigFractionB_PerBoundaryTT,        regionGridWidth/2.0, None, None, None, "SigFrac\'B\'",              "TT_SigFracB_PerBoundary",            njet, color)
        #plotter["TT"].plot_VarVsBoundary(sigFractionC_PerBoundaryTT,        regionGridWidth/2.0, None, None, None, "SigFrac\'C\'",              "TT_SigFracC_PerBoundary",            njet, color)
        #plotter["TT"].plot_VarVsBoundary(sigFractionD_PerBoundaryTT,        regionGridWidth/2.0, None, None, None, "SigFrac\'D\'",              "TT_SigFracD_PerBoundary",            njet, color)
        #plotter["TT"].plot_VarVsBoundary(closureCorrPerBoundaryTT,          regionGridWidth/2.0, yMin, yMax, 1.0,   "Closure Correction [TT]",   "TT_ClosureCorrection_PerBoundary",   njet, color)
        #plotter[ttVar].plot_VarVsBoundary(closureCorrPerBoundaryTTvar,      regionGridWidth/2.0, yMin, yMax, 1.0,   "Closure Correction [%s]"%(ttVar),   "%s_ClosureCorrection_PerBoundary"%(ttVar),   njet, color)

        # TTinData = Data - NonTT
        #plotter["Data"].plot_VarVsBoundary(weighted_nTTEventsA_PerBoundaryTTinData,      regionGridWidth/2.0, None, None, None, "Weighted Events in A [Data]",  "TTinData_weighted_nTTeventsA_PerBoundary",     njet, color)
        #plotter["Data"].plot_VarVsBoundary(nonClosurePerBoundaryTTinData,                regionGridWidth/2.0, 0.0, 0.3,   None, "Non-Closure [Data]",           "TTinData_NonClosure_PerBoundary",              njet, color)
        #plotter["Data"].plot_VarVsBoundary(closurePerBoundaryTTinData,                         regionGridWidth/2.0, yMin, yMax,  1.0, "Data Closure",      "TTinData_DataClosure_PerBoundary",                         njet, color)
        #plotter["Data"].plot_VarVsBoundary(MC_corrected_dataClosure_PerBoundaryTTinData,       regionGridWidth/2.0, yMin, yMax,  1.0, "MC Corrected Data", "TTinData_MC_corrected_dataClosure_PerBoundary",            njet, color)
        #plotter["Data"].plot_VarVsBoundary(MC_ttVar_corrected_dataClosure_PerBoundaryTTinData, regionGridWidth/2.0, yMin, yMax,  1.0, "MC Corrected Data", "TTinData_MC_%s_corrected_dataClosure_PerBoundary"%(ttVar), njet, color)

if __name__ == '__main__':
    main()
