import ROOT
import math
import ctypes
import numpy as np

from DoubleDisCo_Helpers    import *
from DoubleDisCo_Quantities import *

# -------------------------------------------------------------------------
# Base class for performing integrals and determining bin edges
# of all possible choices in an arbitrary region of the full "ABCD" plane
# The arbitrary region is defined by the boundary variables and fixed edges
# Derived classes for a particular region are provided below
# -------------------------------------------------------------------------

# histBkg        : ROOT TH2 corresponding to background process
# histSig        : ROOT TH2 corresponding to signal process
# fixedDisc1Edge : fix the disc 1 edge with provided value while let disc 2 float
# fixedDisc2Edge : fix the disc 2 edge with provided value while let disc 1 float
# leftBoundary   : value corresponding to left edge defining the "ABCD" region
# rightBoundary  : value corresponding to right edge defining the "ABCD" region
# topBoundary    : value corresponding to the top edge defining the "ABCD" region
# bottomBoundary : value corresponding to the bottom edge defining the "ABCD" region
class Regions:

    def __init__(self, histBkg=None, histSig=None, fixedDisc1Edge=None, fixedDisc2Edge=None, leftBoundary=None, rightBoundary=None, topBoundary=None, bottomBoundary=None, **kwargs):

        self.histBkg = histBkg
        self.histSig = histSig

        self.fixedDisc1Edge = fixedDisc1Edge
        self.fixedDisc2Edge = fixedDisc2Edge
        self.leftBoundary   = leftBoundary
        self.rightBoundary  = rightBoundary
        self.topBoundary    = topBoundary
        self.bottomBoundary = bottomBoundary

        self.info = Quantities()

        self.extraArgs = kwargs

        # First calculate events counts for all possible choices of bin edges
        self.count_Events_inBinEdges()

        # Then determine the "final" choice of bin edges for the region
        # based on the optimization_metric function
        # Likewise, store useful quantities per choice of bin edges
        self.get_BinEdges()

    # -----------------------------------------------------
    # Optimization metric function to be overridden by derived class
    # -----------------------------------------------------
    def optimization_metric(self, **kwargs):
        return 0.0

    # -------------------------------------------------------
    # get signal and background histograms' counts
    #   -- count both signal and background events separately 
    #   -- in each A, B, C, D regions
    #   -- put them to the dictionaries
    # -------------------------------------------------------
    def count_Events_inBinEdges(self):
    
        # Define the absolute bounds that our integral takes place in
        # if boundaries are not defined, then assume the edge of the histogram
        firstXBin = 1 if self.leftBoundary   == None else self.histBkg.GetXaxis().FindBin(float(self.leftBoundary))
        firstYBin = 1 if self.bottomBoundary == None else self.histBkg.GetYaxis().FindBin(float(self.bottomBoundary))
 
        lastXBin = self.histBkg.GetNbinsX()+1 if self.rightBoundary == None else self.histBkg.GetXaxis().FindBin(float(self.rightBoundary))
        lastYBin = self.histBkg.GetNbinsY()+1 if self.topBoundary   == None else self.histBkg.GetYaxis().FindBin(float(self.topBoundary))

        nXBins = range(firstXBin+1, lastXBin)
        nYBins = range(firstYBin+1, lastYBin)

        # loop over the x bins i.e. choice of disc 1 as an edge
        for xBin in nXBins:

            # Store disc 1 edge as string with three digits of accuracy for now
            xLowBinEdge = self.histBkg.GetXaxis().GetBinLowEdge(xBin)
            xBinKey = "%.3f"%(xLowBinEdge)
    
            # loop over the y bins
            for yBin in nYBins:
    
                # Store disc 2 edge as string with three digits of accuracy for now
                yLowBinEdge = self.histBkg.GetYaxis().GetBinLowEdge(yBin)
                yBinKey     = "%.3f"%(yLowBinEdge)

                # count signal and background events and errors in bin edges
                nSigEventsErr_A = ROOT.Double(0.0); nSigEventsErr_B = ROOT.Double(0.0); nSigEventsErr_C = ROOT.Double(0.0); nSigEventsErr_D = ROOT.Double(0.0)
                nBkgEventsErr_A = ROOT.Double(0.0); nBkgEventsErr_B = ROOT.Double(0.0); nBkgEventsErr_C = ROOT.Double(0.0); nBkgEventsErr_D = ROOT.Double(0.0)

                # last      | 
                #        B  |  A
                # yBin _____|_____
                #           |
                #        D  |  C
                #           |
                #    
                # first   xBin   last
                nSigEvents_A = math.ceil(self.histSig.IntegralAndError(xBin,      lastXBin, yBin,      lastYBin, nSigEventsErr_A))
                nSigEvents_B = math.ceil(self.histSig.IntegralAndError(firstXBin, xBin-1,   yBin,      lastYBin, nSigEventsErr_B))
                nSigEvents_C = math.ceil(self.histSig.IntegralAndError(xBin,      lastXBin, firstYBin, yBin-1,   nSigEventsErr_C))
                nSigEvents_D = math.ceil(self.histSig.IntegralAndError(firstXBin, xBin-1,   firstYBin, yBin-1,   nSigEventsErr_D))

                nBkgEvents_A = math.ceil(self.histBkg.IntegralAndError(xBin,      lastXBin, yBin,      lastYBin, nBkgEventsErr_A))
                nBkgEvents_B = math.ceil(self.histBkg.IntegralAndError(firstXBin, xBin-1,   yBin,      lastYBin, nBkgEventsErr_B))
                nBkgEvents_C = math.ceil(self.histBkg.IntegralAndError(xBin,      lastXBin, firstYBin, yBin-1,   nBkgEventsErr_C))
                nBkgEvents_D = math.ceil(self.histBkg.IntegralAndError(firstXBin, xBin-1,   firstYBin, yBin-1,   nBkgEventsErr_D))

                self.info.add("nSigEventsA", xBinKey, yBinKey, (nSigEvents_A, nSigEvents_A**0.5))
                self.info.add("nSigEventsB", xBinKey, yBinKey, (nSigEvents_B, nSigEvents_B**0.5))
                self.info.add("nSigEventsC", xBinKey, yBinKey, (nSigEvents_C, nSigEvents_C**0.5))
                self.info.add("nSigEventsD", xBinKey, yBinKey, (nSigEvents_D, nSigEvents_D**0.5))

                self.info.add("nBkgEventsA", xBinKey, yBinKey, (nBkgEvents_A, nBkgEvents_A**0.5))
                self.info.add("nBkgEventsB", xBinKey, yBinKey, (nBkgEvents_B, nBkgEvents_B**0.5))
                self.info.add("nBkgEventsC", xBinKey, yBinKey, (nBkgEvents_C, nBkgEvents_C**0.5))
                self.info.add("nBkgEventsD", xBinKey, yBinKey, (nBkgEvents_D, nBkgEvents_D**0.5))

    # -------------------------------------------------------------------------------
    # Region by region calculation of signal fraction, closure Err, significance, etc
    # look at (SIG + BKG) events in each region, how much signal is in each region
    #   -- SigFragA = Nsig / (Nbkg + Nsig) 
    #   -- SigFragB = Nsig / (Nbkg + Nsig) 
    #   -- SigFragC = Nsig / (Nbkg + Nsig) 
    #   -- SigFragD = Nsig / (Nbkg + Nsig) 
    # ----------------------------------------------------------------------------
    def get_BinEdges(self, bkgNormUnc = 0.3, minBkgEvents = 5, minSigEvents = 5):
      
        optMetric = 999.0

        # loop over the disc1 and disc2 to get any possible combination of them
        for disc1Key, disc2Key in self.info.get("edges"):

            # number of signal and background events in aech A, B, C, D region
            nSigEvents_A, nSigEventsErr_A = self.info.get("nSigEventsA", disc1Key, disc2Key) 
            nSigEvents_B, nSigEventsErr_B = self.info.get("nSigEventsB", disc1Key, disc2Key) 
            nSigEvents_C, nSigEventsErr_C = self.info.get("nSigEventsC", disc1Key, disc2Key) 
            nSigEvents_D, nSigEventsErr_D = self.info.get("nSigEventsD", disc1Key, disc2Key) 

            nBkgEvents_A, nBkgEventsErr_A = self.info.get("nBkgEventsA", disc1Key, disc2Key)
            nBkgEvents_B, nBkgEventsErr_B = self.info.get("nBkgEventsB", disc1Key, disc2Key)
            nBkgEvents_C, nBkgEventsErr_C = self.info.get("nBkgEventsC", disc1Key, disc2Key)
            nBkgEvents_D, nBkgEventsErr_D = self.info.get("nBkgEventsD", disc1Key, disc2Key)

            # Region by region signal fraction calculation
            # get some plots based on region by region signal fraction
            nTot_SigBkg_A = nSigEvents_A + nBkgEvents_A
            nTot_SigBkg_B = nSigEvents_B + nBkgEvents_B
            nTot_SigBkg_C = nSigEvents_C + nBkgEvents_C
            nTot_SigBkg_D = nSigEvents_D + nBkgEvents_D

            sigFracsA    = -1.0; sigFracsB    = -1.0; sigFracsC    = -1.0; sigFracsD    = -1.0
            sigFracsErrA = -1.0; sigFracsErrB = -1.0; sigFracsErrC = -1.0; sigFracsErrD = -1.0
    
            if nTot_SigBkg_A > 0.0: 
                sigFracsA    = nSigEvents_A / nTot_SigBkg_A
                sigFracsErrA = ((nBkgEvents_A * nSigEventsErr_A**0.5 / (nTot_SigBkg_A)**2.0)**2.0 + \
                                (nSigEvents_A * nBkgEventsErr_A**0.5 / (nTot_SigBkg_A)**2.0)**2.0)**0.5
            
            if nTot_SigBkg_B > 0.0: 
                sigFracsB    = nSigEvents_B / nTot_SigBkg_B
                sigFracsErrB = ((nBkgEvents_B * nSigEventsErr_B**0.5 / (nTot_SigBkg_B)**2.0)**2.0 + \
                                (nSigEvents_B * nBkgEventsErr_B**0.5 / (nTot_SigBkg_B)**2.0)**2.0)**0.5
            
            if nTot_SigBkg_C > 0.0: 
                sigFracsC    = nSigEvents_C / nTot_SigBkg_C
                sigFracsErrC = ((nBkgEvents_C * nSigEventsErr_C**0.5 / (nTot_SigBkg_C)**2.0)**2.0 + \
                                (nSigEvents_C * nBkgEventsErr_C**0.5 / (nTot_SigBkg_C)**2.0)**2.0)**0.5
    
            if nTot_SigBkg_D > 0.0: 
                sigFracsD    = nSigEvents_D / nTot_SigBkg_D
                sigFracsErrD = ((nBkgEvents_D * nSigEventsErr_D**0.5 / (nTot_SigBkg_D)**2.0)**2.0 + \
                                (nSigEvents_D * nBkgEventsErr_D**0.5 / (nTot_SigBkg_D)**2.0)**2.0)**0.5
    
            # significance and closure error for optimization of bin edges
            significance = 0.0; significanceUnc = 0.0; closureErr = -999.0; closureErrUnc = -999.0; tempOptMetric = 999.0
    
            closureErr, closureErrUnc = cal_ClosureError(nBkgEvents_A, nBkgEvents_B, nBkgEvents_C, nBkgEvents_D, nBkgEventsErr_A, nBkgEventsErr_B, nBkgEventsErr_C, nBkgEventsErr_D)

            significance += cal_Significance(nSigEvents_A, nBkgEvents_A)**2.0
            significance = significance**0.5
   
            # get the significanceUncs to plot variable vs disc as 1D 
            if nBkgEvents_A > 0.0:
                significanceUnc = (   ( nSigEventsErr_A / (nBkgEvents_A + (bkgNormUnc * nBkgEvents_A)**2.0 + (closureErr * nBkgEvents_A)**2.0)**0.5 )**2.0
                                  + ( ( nSigEvents_A * nBkgEventsErr_A * (2.0 * nBkgEvents_A * closureErr**2.0 + 2.0 * bkgNormUnc**2.0 * nBkgEvents_A + 1) ) / ( nBkgEvents_A + (bkgNormUnc * nBkgEvents_A)**2.0 + (closureErr * nBkgEvents_A)**2.0 )**1.5 )**2.0
                                  + ( ( nBkgEvents_A**2.0 * closureErr * nSigEvents_A * closureErrUnc) / ( nBkgEvents_A * ( nBkgEvents_A * (closureErr**2.0 + bkgNormUnc**2.0) + 1) )**1.5 )**2.0 )**0.5
 
            # Store significance and closure error for current choice of bin edges
            self.info.add("significance", disc1Key, disc2Key, (significance, significanceUnc))
            self.info.add("closureError", disc1Key, disc2Key, (closureErr,   closureErrUnc))

            # Region by region fraction and fraction error
            self.info.add("sigFractionA", disc1Key, disc2Key, (sigFracsA, sigFracsErrA))
            self.info.add("sigFractionB", disc1Key, disc2Key, (sigFracsB, sigFracsErrB))
            self.info.add("sigFractionC", disc1Key, disc2Key, (sigFracsC, sigFracsErrC))
            self.info.add("sigFractionD", disc1Key, disc2Key, (sigFracsD, sigFracsErrD))
   
            # If in region with minimum number of bkg events, compute the metric value
            if (nBkgEvents_A > minBkgEvents):
                tempOptMetric = self.optimization_metric(significance=significance, closureError=closureErr, 
                                                         nBkgA=nBkgEvents_A, nBkgB=nBkgEvents_B, nBkgC=nBkgEvents_C, nBkgD=nBkgEvents_D,
                                                         disc1=float(disc1Key), disc2=float(disc2Key))

            # After calculating interesting quantities above for current choice of bin edges,
            # we check here if this choice is the new best choice
            # If we care about fixed edges check if at fixed edge before checking optimization
            # If our region contains a fixed vertical edge (disc 1), then skip all other choices of x
            # that do not correspond to that edge
            if self.fixedDisc1Edge != None and abs(float(self.fixedDisc1Edge) - float(disc1Key)) > 0.01: continue

            # If our region contains a fixed horizontal edge (disc 2), then skip all other choices of y
            # that do not correspond to that edge
            if self.fixedDisc2Edge != None and abs(float(self.fixedDisc2Edge) - float(disc2Key)) > 0.01: continue

            # Determine based on the metric value if the current
            # choice of bin edges is better and if so, save them
            if tempOptMetric < optMetric:
                self.info.finalEdges = (disc1Key, disc2Key)
                optMetric = tempOptMetric

    # A couple of simple helper functions to retrive a quantity
    # i.e. closureErr, significance, etc.
    def get(self, name):
        return self.info.get(name)

    def getFinal(self, name):
        return self.info.getFinal(name)
    
# Derived from the Regions class, this class is for 
# doing things for the full ABCD region.
# The optimization function takes into account closure err
# and significance.
class ABCDedges(Regions):

    def optimization_metric(self, **kwargs):

        optimizationMetric  = (5 * kwargs["closureError"])**2 + (1.0 / kwargs["significance"])**2
        return optimizationMetric

# Derived from the Regions class, this class is for            *               
# doing things for the bdEF validation region.             E   *   B'  |   A   
# The optimization function tries to match bkg A+C to   _______*_______|_______
# to bkg B'+D'.                                                *       |       
#                                                          F   *   D'  |   C   
#                                                              *                
# The extraArgs contains info from ABCD region
class bdEFedges(Regions):

    def optimization_metric(self, **kwargs):

        optimizationMetric  = abs(1.0 - (kwargs["nBkgA"]+kwargs["nBkgC"])/(self.extraArgs["nBkgA"]+self.extraArgs["nBkgC"]))
        return optimizationMetric

# Derived from the Regions class, this class is for      B  |  A 
# doing things for the cdiEF validation region.        _____|_____
# The optimization function tries to match bkg A+B to       |    
# to bkg C'+D'.                                          D' |  C'
#                                                      ***********
#                                                        H  |  G 
#                                                           |   
# The extraArgs contains info from ABCD region
class cdGHedges(Regions):

    def optimization_metric(self, **kwargs):

        optimizationMetric  = abs(1.0 - (kwargs["nBkgA"]+kwargs["nBkgB"])/(self.extraArgs["nBkgA"]+self.extraArgs["nBkgB"]))
        return optimizationMetric

# Derived from the Regions class, this class is for  
# doing things for the SubDiv D validation region.      
# The optimization function simply checks if
# current edges are at half the value of ABCD edges
class subDivDedges(Regions):

    def optimization_metric(self, **kwargs):

        optimizationMetric  = 999.0

        if abs(kwargs["disc1"] - self.extraArgs["ABCDdisc1"]/2.0) < 0.01 and abs(kwargs["disc2"] - self.extraArgs["ABCDdisc2"]/2.0) < 0.01:
            optimizationMetric = 1.0

        return optimizationMetric
