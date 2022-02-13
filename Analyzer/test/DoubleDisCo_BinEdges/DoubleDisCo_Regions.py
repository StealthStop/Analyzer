import ROOT
import os
import math
import ctypes
import numpy as np


# ----------------------------------------------------------------------------------
# histBkg        : ROOT TH2 corresponding to background process
# histSig        : ROOT TH2 corresponding to signal process
# fixedDisc1Edge : fix the disc 1 edge with provided value 
# fixedDisc2Edge : fix the disc 2 edge with provided value 
# leftBoundary   : value corresponding to the left edge defining the "ABCD" region
# rightBoundary  : value corresponding to the right edge defining the "ABCD" region
# topBoundary    : value corresponding to the top edge defining the "ABCD" region
# bottomBoundary : value corresponding to the bottom edge defining the "ABCD" region
# ----------------------------------------------------------------------------------
class All_Regions:

    def __init__(self, hist=None, Sig=None, ttVar=None, fixedDisc1Edge=None, fixedDisc2Edge=None, leftBoundary=None, rightBoundary=None, topBoundary=None, bottomBoundary=None, metric=None, **kwargs):

        self.hist           = hist
        self.Sig            = Sig
        self.ttVar          = ttVar
        self.fixedDisc1Edge = fixedDisc1Edge
        self.fixedDisc2Edge = fixedDisc2Edge
        self.leftBoundary   = leftBoundary
        self.rightBoundary  = rightBoundary
        self.topBoundary    = topBoundary
        self.bottomBoundary = bottomBoundary
        self.metric         = metric

        self.extraArgs = kwargs
        self.finalEdges = (fixedDisc1Edge, fixedDisc2Edge)
        self.quantities = {}

        for key in hist.keys():
            self.quantities[key] = {}

        # First calculate events counts for all possible choices of bin edges
        self.count_Events_inBinEdges()

        # Then determine the "final" choice of bin edges for the region
        # based on the optimization_metric function
        self.get_nEvents_Quantities()

    # -------------------------------------
    # Significance calculation with only TT
    # -------------------------------------
    def cal_Significance(self, nSigEvents, nTTEvents, sys=0.3):
        if (nTTEvents == 0.0):
            return 0.0
    
        significance = nSigEvents / ( nTTEvents + (sys * nTTEvents)**2.0 )**0.5
        return significance
    
    # -----------------------
    # Non-Closure calculation
    # -----------------------
    def cal_NonClosure(self, nEvents_A, nEvents_B, nEvents_C, nEvents_D, nEventsErr_A, nEventsErr_B, nEventsErr_C, nEventsErr_D):
    
        if nEvents_A == 0.0 or nEvents_D == 0.0:
            return -999.0, -999.0
        
        nonClosure = abs(1.0 - ( (nEvents_B * nEvents_C) / (nEvents_A * nEvents_D) ) )
        
        nonClosureUnc = ( ( ( nEvents_C * nEventsErr_B ) / ( nEvents_A * nEvents_D) )**2.0 
                        + ( ( nEvents_B * nEventsErr_C ) / ( nEvents_A * nEvents_D) )**2.0 
                        + ( ( nEvents_B * nEvents_C * nEventsErr_A ) / ( nEvents_A**2.0 * nEvents_D ) )**2.0 
                        + ( ( nEvents_B * nEvents_C * nEventsErr_D ) / ( nEvents_A * nEvents_D**2.0 ) )**2.0 )**0.5
    
        return nonClosure, nonClosureUnc

    def cal_McCorrection(self, nEvents_A, nEvents_B, nEvents_C, nEvents_D, nEventsErr_A, nEventsErr_B, nEventsErr_C, nEventsErr_D):
    
        if nEvents_B == 0.0 or nEvents_C == 0.0:
            return -999.0, -999.0
        
        mcCorr = (nEvents_A * nEvents_D) / (nEvents_B * nEvents_C)
        
        mcCorrUnc = ( ( ( nEvents_A * nEventsErr_D ) / ( nEvents_B * nEvents_C) )**2.0 
                        + ( ( nEvents_D * nEventsErr_A ) / ( nEvents_B * nEvents_C) )**2.0 
                        + ( ( nEvents_D * nEvents_A * nEventsErr_A ) / ( nEvents_B**2.0 * nEvents_C ) )**2.0 
                        + ( ( nEvents_D * nEvents_A * nEventsErr_C ) / ( nEvents_B * nEvents_C**2.0 ) )**2.0 )**0.5
    
        return mcCorr, mcCorrUnc

    # ------------------------
    # Closure Pull calculation
    # ------------------------
    def cal_Pull(self, nEvents_A, nEvents_B, nEvents_C, nEvents_D, nEventsErr_A, nEventsErr_B, nEventsErr_C, nEventsErr_D):

        if nEvents_D == 0.0 or nEventsErr_A == 0.0:
            return -999.0, -999.0

        pull    = ( ( ((nEvents_B * nEvents_C) / nEvents_D) - nEvents_A ) / nEventsErr_A )
        pullUnc = 1.0

        return pull, pullUnc
    
    # --------------------------------------------------------------
    # Optimization metric function to be overridden by derived class
    # --------------------------------------------------------------
    def optimization_metric(self, **kwargs):
        return 0.0

    # -----------------------------------------------------
    # get signal and background histograms' counts
    #   count both signal and background events separately:
    #   -- in each A, B, C, D regions
    #   -- B'D'EF, C'D'GH validation regions
    #   -- Sub-division D validation regions
    # -----------------------------------------------------
    def count_Events_inBinEdges(self):
    
        # Define the absolute bounds that our integral takes place in
        # If boundaries are not defined, then assume the edge of the histogram
        firstXBin = 1 if self.leftBoundary   == None else self.hist["TT"].GetXaxis().FindBin(float(self.leftBoundary))
        firstYBin = 1 if self.bottomBoundary == None else self.hist["TT"].GetYaxis().FindBin(float(self.bottomBoundary))
 
        lastXBin = self.hist["TT"].GetNbinsX()+1 if self.rightBoundary == None else self.hist["TT"].GetXaxis().FindBin(float(self.rightBoundary)) + 1
        lastYBin = self.hist["TT"].GetNbinsY()+1 if self.topBoundary   == None else self.hist["TT"].GetYaxis().FindBin(float(self.topBoundary))   + 1

        nXBins = range(firstXBin+1, lastXBin)
        nYBins = range(firstYBin+1, lastYBin)

        # loop over the x bins i.e. choice of disc 1 as an edge
        for xBin in nXBins:

            # Store disc 1 edge as string with three digits of accuracy for now
            xLowBinEdge = self.hist["TT"].GetXaxis().GetBinLowEdge(xBin)
            xBinKey     = "%.3f"%(xLowBinEdge)
    
            # loop over the y bins
            for yBin in nYBins:
    
                # Store disc 2 edge as string with three digits of accuracy for now
                yLowBinEdge = self.hist["TT"].GetYaxis().GetBinLowEdge(yBin)
                yBinKey     = "%.3f"%(yLowBinEdge)

                # count signal and background events and errors in bin edges
                for key, h1 in self.hist.items():
                    
                    nEventsErr_A = ROOT.Double(0.0); nEventsErr_B = ROOT.Double(0.0); nEventsErr_C = ROOT.Double(0.0); nEventsErr_D = ROOT.Double(0.0)

                    # last      | 
                    #        B  |  A
                    # yBin _____|_____
                    #           |
                    #        D  |  C
                    #           |
                    #    
                    # first   xBin   last
                    nEvents_A = math.ceil(h1.IntegralAndError(xBin,      lastXBin, yBin,      lastYBin, nEventsErr_A))
                    nEvents_B = math.ceil(h1.IntegralAndError(firstXBin, xBin-1,   yBin,      lastYBin, nEventsErr_B))
                    nEvents_C = math.ceil(h1.IntegralAndError(xBin,      lastXBin, firstYBin, yBin-1,   nEventsErr_C))
                    nEvents_D = math.ceil(h1.IntegralAndError(firstXBin, xBin-1,   firstYBin, yBin-1,   nEventsErr_D))

                    self.add("nEventsA", xBinKey, yBinKey, (nEvents_A, nEvents_A**0.5), key)
                    self.add("nEventsB", xBinKey, yBinKey, (nEvents_B, nEvents_B**0.5), key)
                    self.add("nEventsC", xBinKey, yBinKey, (nEvents_C, nEvents_C**0.5), key)
                    self.add("nEventsD", xBinKey, yBinKey, (nEvents_D, nEvents_D**0.5), key)

    # -------------------------------------------------------------------------------
    # Region by region calculation of signal fraction, closure Err, significance, etc
    # look at (SIG + BKG) events in each region, how much signal is in each region
    #   -- SigFragA = Nsig / (Nbkg + Nsig) 
    #   -- SigFragB = Nsig / (Nbkg + Nsig) 
    #   -- SigFragC = Nsig / (Nbkg + Nsig) 
    #   -- SigFragD = Nsig / (Nbkg + Nsig)
    # Total (A+B+C+D) signal fraction calculation
    # look at (SIG) events in each region, how many events are signal
    #   NtotalSig = NA + NB + NC + ND
    #   -- SigTotFragA = NA / Ntotal
    #   -- SigTotFragB = NB / Ntotal
    #   -- SigTotFragC = NC / Ntotal
    #   -- SigTotFragD = ND / Ntotal  
    # ----------------------------------------------------------------------------
    def get_nEvents_Quantities(self, bkgNormUnc = 0.3, minBkgEvents = 1, minSigEvents = 5):

        # loop over the disc1 and disc2 to get any possible combination of them
        for disc1Key, disc2Key in self.get("edges",None,None,"TT"):

            # number of signal and background events in aech A, B, C, D region
            nTTEvents_A,    nTTEventsErr_A    = self.get("nEventsA", disc1Key, disc2Key, "TT")
            nTTEvents_B,    nTTEventsErr_B    = self.get("nEventsB", disc1Key, disc2Key, "TT")
            nTTEvents_C,    nTTEventsErr_C    = self.get("nEventsC", disc1Key, disc2Key, "TT")
            nTTEvents_D,    nTTEventsErr_D    = self.get("nEventsD", disc1Key, disc2Key, "TT")

            nNonTTEvents_A, nNonTTEventsErr_A = self.get("nEventsA", disc1Key, disc2Key, "NonTT")
            nNonTTEvents_B, nNonTTEventsErr_B = self.get("nEventsB", disc1Key, disc2Key, "NonTT")
            nNonTTEvents_C, nNonTTEventsErr_C = self.get("nEventsC", disc1Key, disc2Key, "NonTT")
            nNonTTEvents_D, nNonTTEventsErr_D = self.get("nEventsD", disc1Key, disc2Key, "NonTT")

            nTTvarEvents_A, nTTvarEventsErr_A = self.get("nEventsA", disc1Key, disc2Key, self.ttVar)
            nTTvarEvents_B, nTTvarEventsErr_B = self.get("nEventsB", disc1Key, disc2Key, self.ttVar)
            nTTvarEvents_C, nTTvarEventsErr_C = self.get("nEventsC", disc1Key, disc2Key, self.ttVar)
            nTTvarEvents_D, nTTvarEventsErr_D = self.get("nEventsD", disc1Key, disc2Key, self.ttVar)

            nSigEvents_A,   nSigEventsErr_A   = self.get("nEventsA", disc1Key, disc2Key, self.Sig) 
            nSigEvents_B,   nSigEventsErr_B   = self.get("nEventsB", disc1Key, disc2Key, self.Sig) 
            nSigEvents_C,   nSigEventsErr_C   = self.get("nEventsC", disc1Key, disc2Key, self.Sig) 
            nSigEvents_D,   nSigEventsErr_D   = self.get("nEventsD", disc1Key, disc2Key, self.Sig) 

            nDataEvents_A,  nDataEventsErr_A  = self.get("nEventsA", disc1Key, disc2Key, "Data")
            nDataEvents_B,  nDataEventsErr_B  = self.get("nEventsB", disc1Key, disc2Key, "Data")
            nDataEvents_C,  nDataEventsErr_C  = self.get("nEventsC", disc1Key, disc2Key, "Data")
            nDataEvents_D,  nDataEventsErr_D  = self.get("nEventsD", disc1Key, disc2Key, "Data")

            # Region by region signal and bkg fractions for 1D-2D plots / for only TT !!!
            nTot_SigTT_A = nSigEvents_A + nTTEvents_A
            nTot_SigTT_B = nSigEvents_B + nTTEvents_B
            nTot_SigTT_C = nSigEvents_C + nTTEvents_C
            nTot_SigTT_D = nSigEvents_D + nTTEvents_D

            sigFracsA    = -1.0; sigFracsB    = -1.0; sigFracsC    = -1.0; sigFracsD    = -1.0
            sigFracsErrA = -1.0; sigFracsErrB = -1.0; sigFracsErrC = -1.0; sigFracsErrD = -1.0
            bkgFracsA    = -1.0; bkgFracsB    = -1.0; bkgFracsC    = -1.0; bkgFracsD    = -1.0
            bkgFracsErrA = -1.0; bkgFracsErrB = -1.0; bkgFracsErrC = -1.0; bkgFracsErrD = -1.0
 
            if nTot_SigTT_A > 0.0: 
                sigFracsA    = nSigEvents_A / nTot_SigTT_A
                sigFracsErrA = ((nTTEvents_A * nSigEventsErr_A**0.5 / (nTot_SigTT_A)**2.0)**2.0 + \
                                (nSigEvents_A * nTTEventsErr_A**0.5 / (nTot_SigTT_A)**2.0)**2.0 )**0.5
                bkgFracsA    = 1 - sigFracsA
                bkgFracsErrA = ((nTTEvents_A * nSigEventsErr_A**0.5 / (nTot_SigTT_A)**2.0)**2.0 + \
                                (nSigEvents_A * nTTEventsErr_A**0.5 / (nTot_SigTT_A)**2.0)**2.0)**0.5
 
            if nTot_SigTT_B > 0.0: 
                sigFracsB    = nSigEvents_B / nTot_SigTT_B
                sigFracsErrB = ((nTTEvents_B * nSigEventsErr_B**0.5 / (nTot_SigTT_B)**2.0)**2.0 + \
                                (nSigEvents_B * nTTEventsErr_B**0.5 / (nTot_SigTT_B)**2.0)**2.0)**0.5
                bkgFracsB    = 1 - sigFracsB
                bkgFracsErrB = ((nTTEvents_B * nSigEventsErr_B**0.5 / (nTot_SigTT_B)**2.0)**2.0 + \
                                (nSigEvents_B * nTTEventsErr_B**0.5 / (nTot_SigTT_B)**2.0)**2.0)**0.5           

            if nTot_SigTT_C > 0.0: 
                sigFracsC    = nSigEvents_C / nTot_SigTT_C
                sigFracsErrC = ((nTTEvents_C * nSigEventsErr_C**0.5 / (nTot_SigTT_C)**2.0)**2.0 + \
                                (nSigEvents_C * nTTEventsErr_C**0.5 / (nTot_SigTT_C)**2.0)**2.0)**0.5
                bkgFracsC    = 1 - sigFracsC 
                bkgFracsErrC = ((nTTEvents_C * nSigEventsErr_C**0.5 / (nTot_SigTT_C)**2.0)**2.0 + \
                                (nSigEvents_C * nTTEventsErr_C**0.5 / (nTot_SigTT_C)**2.0)**2.0)**0.5
    
            if nTot_SigTT_D > 0.0: 
                sigFracsD    = nSigEvents_D / nTot_SigTT_D
                sigFracsErrD = ((nTTEvents_D * nSigEventsErr_D**0.5 / (nTot_SigTT_D)**2.0)**2.0 + \
                                (nSigEvents_D * nTTEventsErr_D**0.5 / (nTot_SigTT_D)**2.0)**2.0)**0.5
                bkgFracsD    = 1 - sigFracsD
                bkgFracsErrD = ((nTTEvents_D * nSigEventsErr_D**0.5 / (nTot_SigTT_D)**2.0)**2.0 + \
                                (nSigEvents_D * nTTEventsErr_D**0.5 / (nTot_SigTT_D)**2.0)**2.0)**0.5   
         
            self.add("sigFractionA", disc1Key, disc2Key, (sigFracsA, sigFracsErrA), self.Sig)
            self.add("sigFractionB", disc1Key, disc2Key, (sigFracsB, sigFracsErrB), self.Sig)
            self.add("sigFractionC", disc1Key, disc2Key, (sigFracsC, sigFracsErrC), self.Sig)
            self.add("sigFractionD", disc1Key, disc2Key, (sigFracsD, sigFracsErrD), self.Sig)   

            self.add("bkgFractionA", disc1Key, disc2Key, (bkgFracsA, bkgFracsErrA), "TT")
            self.add("bkgFractionB", disc1Key, disc2Key, (bkgFracsB, bkgFracsErrB), "TT")
            self.add("bkgFractionC", disc1Key, disc2Key, (bkgFracsC, bkgFracsErrC), "TT")
            self.add("bkgFractionD", disc1Key, disc2Key, (bkgFracsD, bkgFracsErrD), "TT")
    
            # TT fraction by looking at (TT + NonTT) events in each A,B,C,D region
            # to understand how NonTT events dominated in each region
            nTot_TT_NonTT_A = nTTEvents_A + nNonTTEvents_A
            nTot_TT_NonTT_B = nTTEvents_B + nNonTTEvents_B
            nTot_TT_NonTT_C = nTTEvents_C + nNonTTEvents_C
            nTot_TT_NonTT_D = nTTEvents_D + nNonTTEvents_D

            ttFracsA    = -1.0; ttFracsB    = -1.0; ttFracsC    = -1.0; ttFracsD    = -1.0
            ttFracsErrA = -1.0; ttFracsErrB = -1.0; ttFracsErrC = -1.0; ttFracsErrD = -1.0

            if nTot_TT_NonTT_A > 0.0:
                ttFracsA    = nTTEvents_A / nTot_TT_NonTT_A
                ttFracsErrA = ( (nTTEvents_A * nNonTTEventsErr_A**0.5 / (nTot_TT_NonTT_A)**2.0 )**2.0 + \
                               (nNonTTEvents_A * nTTEventsErr_A**0.5 / (nTot_TT_NonTT_A)**2.0)**2.0 )**0.5

            if nTot_TT_NonTT_B > 0.0:
                ttFracsB    = nTTEvents_B / nTot_TT_NonTT_B
                ttFracsErrB = ( (nTTEvents_B * nNonTTEventsErr_B**0.5 / (nTot_TT_NonTT_B)**2.0 )**2.0 + \
                               (nNonTTEvents_B * nTTEventsErr_B**0.5 / (nTot_TT_NonTT_B)**2.0)**2.0 )**0.5

            if nTot_TT_NonTT_C > 0.0:
                ttFracsC    = nTTEvents_C / nTot_TT_NonTT_C
                ttFracsErrC = ( (nTTEvents_C * nNonTTEventsErr_C**0.5 / (nTot_TT_NonTT_C)**2.0 )**2.0 + \
                               (nNonTTEvents_C * nTTEventsErr_C**0.5 / (nTot_TT_NonTT_C)**2.0)**2.0 )**0.5

            if nTot_TT_NonTT_D > 0.0:
                ttFracsD    = nTTEvents_D / nTot_TT_NonTT_D
                ttFracsErrD = ( (nTTEvents_D * nNonTTEventsErr_D**0.5 / (nTot_TT_NonTT_D)**2.0 )**2.0 + \
                               (nNonTTEvents_D * nTTEventsErr_D**0.5 / (nTot_TT_NonTT_D)**2.0)**2.0 )**0.5

            self.add("ttFractionA", disc1Key, disc2Key, (ttFracsA, ttFracsErrA), "TT")
            self.add("ttFractionB", disc1Key, disc2Key, (ttFracsB, ttFracsErrB), "TT")
            self.add("ttFractionC", disc1Key, disc2Key, (ttFracsC, ttFracsErrC), "TT")
            self.add("ttFractionD", disc1Key, disc2Key, (ttFracsD, ttFracsErrD), "TT")

            # closure for optimization metric for only TT
            # closure error and pull for 2D plots / for TT, NonTT, Data
            closureErr_TT    = -999.0; closureErrUnc_TT    = -999.0; pull_TT    = -999.0; pullUnc_TT    = -999.0 
            closureErr_NonTT = -999.0; closureErrUnc_NonTT = -999.0; pull_NonTT = -999.0; pullUnc_NonTT = -999.0 
            closureErr_TTvar = -999.0; closureErrUnc_TTvar = -999.0; pull_TTvar = -999.0; pullUnc_TTvar = -999.0
            closureErr_Data  = -999.0; closureErrUnc_Data  = -999.0; pull_Data  = -999.0; pullUnc_Data  = -999.0
            mcCorrection_TT    = -999.0; mcCorrectionUnc_TT    = -999.0
            mcCorrection_NonTT = -999.0; mcCorrectionUnc_NonTT = -999.0
            mcCorrection_TTvar = -999.0; mcCorrectionUnc_TTvar = -999.0
            mcCorrection_Data  = -999.0; mcCorrectionUnc_Data  = -999.0
    
            nonClosure_TT,    nonClosureUnc_TT    = self.cal_NonClosure(nTTEvents_A,    nTTEvents_B,    nTTEvents_C,    nTTEvents_D,    nTTEventsErr_A,    nTTEventsErr_B,    nTTEventsErr_C,    nTTEventsErr_D   )
            nonClosure_NonTT, nonClosureUnc_NonTT = self.cal_NonClosure(nNonTTEvents_A, nNonTTEvents_B, nNonTTEvents_C, nNonTTEvents_D, nNonTTEventsErr_A, nNonTTEventsErr_B, nNonTTEventsErr_C, nNonTTEventsErr_D)
            nonClosure_TTvar, nonClosureUnc_TTvar = self.cal_NonClosure(nTTvarEvents_A, nTTvarEvents_B, nTTvarEvents_C, nTTvarEvents_D, nTTvarEventsErr_A, nTTvarEventsErr_B, nTTvarEventsErr_C, nTTvarEventsErr_D)
            nonClosure_Data,  nonClosureUnc_Data  = self.cal_NonClosure(nDataEvents_A,  nDataEvents_B,  nDataEvents_C,  nDataEvents_D,  nDataEventsErr_A,  nDataEventsErr_B,  nDataEventsErr_C,  nDataEventsErr_D )

            mcCorrection_TT,    mcCorrectionUnc_TT    = self.cal_McCorrection(nTTEvents_A,    nTTEvents_B,    nTTEvents_C,    nTTEvents_D,    nTTEventsErr_A,    nTTEventsErr_B,    nTTEventsErr_C,    nTTEventsErr_D   )
            mcCorrection_Data,  mcCorrectionUnc_Data  = self.cal_McCorrection(nDataEvents_A,  nDataEvents_B,  nDataEvents_C,  nDataEvents_D,  nDataEventsErr_A,  nDataEventsErr_B,  nDataEventsErr_C,  nDataEventsErr_D )
            mcCorrection_NonTT, mcCorrectionUnc_NonTT = self.cal_McCorrection(nNonTTEvents_A, nNonTTEvents_B, nNonTTEvents_C, nNonTTEvents_D, nNonTTEventsErr_A, nNonTTEventsErr_B, nNonTTEventsErr_C, nNonTTEventsErr_D   )
            mcCorrection_TTvar, mcCorrectionUnc_TTvar = self.cal_McCorrection(nTTvarEvents_A, nTTvarEvents_B, nTTvarEvents_C, nTTvarEvents_D, nTTvarEventsErr_A, nTTvarEventsErr_B, nTTvarEventsErr_C, nTTvarEventsErr_D )

            pull_TT,    pullUnc_TT    = self.cal_Pull(nTTEvents_A,    nTTEvents_B,    nTTEvents_C,    nTTEvents_D,    nTTEventsErr_A,    nTTEventsErr_B,    nTTEventsErr_C,    nTTEventsErr_D   )
            pull_NonTT, pullUnc_NonTT = self.cal_Pull(nNonTTEvents_A, nNonTTEvents_B, nNonTTEvents_C, nNonTTEvents_D, nNonTTEventsErr_A, nNonTTEventsErr_B, nNonTTEventsErr_C, nNonTTEventsErr_D)
            pull_TTvar, pullUnc_TTvar = self.cal_Pull(nTTvarEvents_A, nTTvarEvents_B, nTTvarEvents_C, nTTvarEvents_D, nTTvarEventsErr_A, nTTvarEventsErr_B, nTTvarEventsErr_C, nTTvarEventsErr_D)
            pull_Data,  pullUnc_Data  = self.cal_Pull(nDataEvents_A,  nDataEvents_B,  nDataEvents_C,  nDataEvents_D,  nDataEventsErr_A,  nDataEventsErr_B,  nDataEventsErr_C,  nDataEventsErr_D )

            self.add("nonClosure", disc1Key, disc2Key, (nonClosure_TT,    nonClosureUnc_TT) ,   "TT"      )
            self.add("nonClosure", disc1Key, disc2Key, (nonClosure_NonTT, nonClosureUnc_NonTT), "NonTT"   )
            self.add("nonClosure", disc1Key, disc2Key, (nonClosure_TTvar, nonClosureUnc_TTvar), self.ttVar)
            self.add("nonClosure", disc1Key, disc2Key, (nonClosure_Data,  nonClosureUnc_Data ), "Data"    )
            self.add("pull",         disc1Key, disc2Key, (pull_TT,    pullUnc_TT),    "TT"      )
            self.add("pull",         disc1Key, disc2Key, (pull_NonTT, pullUnc_NonTT), "NonTT"   )
            self.add("pull",         disc1Key, disc2Key, (pull_TTvar, pullUnc_TTvar), self.ttVar)
            self.add("pull",         disc1Key, disc2Key, (pull_Data,  pullUnc_Data),  "Data"    )
       
            self.add("mcCorrection", disc1Key, disc2Key, (mcCorrection_TT,    mcCorrectionUnc_TT) ,    "TT"      )
            self.add("mcCorrection", disc1Key, disc2Key, (mcCorrection_Data,  mcCorrectionUnc_Data ),  "Data"    )
            self.add("mcCorrection", disc1Key, disc2Key, (mcCorrection_NonTT, mcCorrectionUnc_NonTT) , "NonTT"   )
            self.add("mcCorrection", disc1Key, disc2Key, (mcCorrection_TTvar, mcCorrectionUnc_TTvar ), self.ttVar)

            # significance for optimization metric for only TT !!! 
            # significance, significanceUnc for 2D plots
            significance_TT = 0.0; significanceUnc_TT = 0.0; tempOptMetric = 999.0

            significance_TT += self.cal_Significance(nSigEvents_A, nTTEvents_A)**2.0
            significance_TT = significance_TT**0.5

            if nTTEvents_A > 0.0:
                significanceUnc_TT = (   ( nSigEventsErr_A / (nTTEvents_A + (bkgNormUnc * nTTEvents_A)**2.0 + (nonClosure_TT * nTTEvents_A)**2.0)**0.5 )**2.0
                                     + ( ( nSigEvents_A * nTTEventsErr_A * (2.0 * nTTEvents_A * nonClosure_TT**2.0 + 2.0 * bkgNormUnc**2.0 * nTTEvents_A + 1) ) / ( nTTEvents_A + (bkgNormUnc * nTTEvents_A)**2.0 + (nonClosure_TT * nTTEvents_A)**2.0 )**1.5 )**2.0
                                     + ( ( nTTEvents_A**2.0 * nonClosure_TT * nSigEvents_A * nonClosureUnc_TT) / ( nTTEvents_A * ( nTTEvents_A * (nonClosure_TT**2.0 + bkgNormUnc**2.0) + 1) )**1.5 )**2.0 )**0.5

            self.add("significance", disc1Key, disc2Key, (significance_TT, significanceUnc_TT), "TT") 

            # use this statement for cdGH regions if  fixed disc1 edge (vertivcal edge)
            if self.fixedDisc1Edge != None and abs(float(self.fixedDisc1Edge) - float(disc1Key)) > 0.01: continue

            # use this statement for bdEF regions if fixed disc2 edge (horizontal edge)
            if self.fixedDisc2Edge != None and abs(float(self.fixedDisc2Edge) - float(disc2Key)) > 0.01: continue

            self.finalEdges = (disc1Key, disc2Key)

    # ----------------------------------
    # store quantities to make any plots
    # with any combination of bin edges
    # i.e. significance, closure, etc.
    # ----------------------------------
    def add(self, var, disc1, disc2, val, hist=""):

        if hist not in self.quantities:
            self.quantities[hist] = {}

        if disc1 not in self.quantities[hist]:
            self.quantities[hist][disc1] = {}

        if disc2 not in self.quantities[hist][disc1]:
            self.quantities[hist][disc1][disc2] = {}

        self.quantities[hist][disc1][disc2][var] = val

    # -------------------------------------------------
    # get specific variables to make any plots
    # with any combination of bin edges
    # i.e. nSigEventsA, nBkgEventsA, sigFractionA, etc.
    # -------------------------------------------------
    def get(self, var, disc1=None, disc2=None, hist="TT"):

        if var == "edges":
            if disc1 != None and disc2 != None:
                return self.finalEdges
            else:
                payload = []                
                for disc1, disc2s in self.quantities[hist].items():
                    for disc2, _ in disc2s.items():
                        payload.append((disc1,disc2))
                return payload

        if disc1 == None and disc2 == None:
            payload = []
            for disc1, disc2s in self.quantities[hist].items():
                for disc2, q in disc2s.items():
                    payload.append(q[var])
            return np.array(payload)

        if disc1 != None and disc1 not in self.quantities[hist]:
            raise LookupError("Cannot find the key \"%s\" in the dictionary !"%(disc1))

        if disc2 != None and disc2 not in self.quantities[hist][disc1]:
            raise LookupError("Cannot find the key \"%s\" in the dictionary !"%(disc2))

        return self.quantities[hist][disc1][disc2][var]
        
    # -------------------------------------------------------------
    # store quatities and get spesific variables to make the tables 
    # with the final choice of bin edges
    # -------------------------------------------------------------
    def getFinal(self, name, hist=""):
        return self.get(name, self.finalEdges[0], self.finalEdges[1], hist)


# -------------------------------------------
# Calculate optimization metric of bin edges
#   This function use the command line option
#   -- NN optimization metric
#   -- New optimization metric
# ------------------------------------------- 
class ABCDedges(All_Regions):

    def optimization_metric(self, **kwargs):

        if (kwargs["nBkgA"] > kwargs["minBkgEvents"] and kwargs["significance"] != 0.0):
        
            # NN optimization metric
            if self.metric == "NN":
                optimizationMetric  = (kwargs["nonClosure"])**2 + (1.0 / kwargs["significance"])**2
                   
            # New optimization metric
            else: 
                optimizationMetric  = (5 * kwargs["nonClosure"])**2 + (1.0 / kwargs["significance"])**2
        
            return optimizationMetric

        else:

            return 999.0

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
# Calculate the metric of Validation Regions in BD: bdEF
# ---------------------------------------------------------------------------
class bdEFedges(All_Regions):

    def optimization_metric(self, **kwargs):

        optimizationMetric = None              

        if "nBkgA" in self.extraArgs:
            optimizationMetric  = abs(1.0 - (kwargs["nBkgA"]+kwargs["nBkgC"])/(self.extraArgs["nBkgA"]+self.extraArgs["nBkgC"]))
        
        else:
            optimizationMetric = 1.0

        return optimizationMetric


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
# Calculate the metric of Validation Regions in CD: cdiGH
# --------------------------------------------------------------------------- 
class cdGHedges(All_Regions):

    def optimization_metric(self, **kwargs):

        optimizationMetric = None

        if "nBkgA" in self.extraArgs:
            optimizationMetric  = abs(1.0 - (kwargs["nBkgA"]+kwargs["nBkgB"])/(self.extraArgs["nBkgA"]+self.extraArgs["nBkgB"]))
        
        else:
            optimizationMetric = 1.0

        return optimizationMetric


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
# Calculate the metric of sub-division of D 
# ------------------------------------------------
class subDivDedges(All_Regions):

    def optimization_metric(self, **kwargs):

        optimizationMetric  = 999.0

        if abs(kwargs["disc1"] - self.extraArgs["ABCDdisc1"]/2.0) < 0.01 and abs(kwargs["disc2"] - self.extraArgs["ABCDdisc2"]/2.0) < 0.01:
            optimizationMetric = 1.0

        return optimizationMetric

# -----------------------------------------
# add the all edges to DoubleDisCo Cfg file
# -----------------------------------------
class addEdges_toDoubleDisco():

    def __init__(self, year, model, mass, channel, regions):
        self.year      = year
        self.model     = model
        self.mass      = mass
        self.channel   = channel
        self.regions   = regions

    def addEdges_toDoubleDiscoCfg(self, edgesPerNjets=None, Njets=None):

        cfgVer = ""
        
        if self.channel == "0l":
            cfgVer = "v2.0"
        else:
            cfgVer = "v4.0"

        outputDir = "DeepESMCfg_DoubleDisCo_Reg_%s_%s_2016_%s/"%(self.channel, self.model, str(cfgVer))
  
        if not os.path.exists(outputDir):
            os.makedirs(outputDir)
 
        f = open("../DeepESMCfg_DoubleDisCo_Reg_%s_%s_2016_%s/DoubleDisCo_Reg.cfg"%(self.channel, self.model, str(cfgVer)), "r")
      
        g = open("%s/DoubleDisCo_Reg.cfg"%(outputDir), "w") 
        
        # get the all information from DoubleDisCo cfg file
        lines = f.readlines()
        f.close()
        
        for line in lines:

            if "}" in line:
                continue

            g.write(line)

        # add the all ABCD and Validation edges to DoubleDisCo cfg file
        for key, region in self.regions.items():

            g.write("\n")
            g.write(" # %s Bin Edges\n"%(region))                        

            i = 0
            for njet in Njets:

                if (njet < 10):
                    g.write("   binEdges_%s[%d] = %s \n" %(region, i, edgesPerNjets[njet][key][0]))
                    i += 1
                    g.write("   binEdges_%s[%d] = %s \n" %(region, i, edgesPerNjets[njet][key][1]))
                else:
                    g.write("   binEdges_%s[%d] = %s \n" %(region, i, edgesPerNjets[njet][key][0]))
                    i += 1
                    g.write("   binEdges_%s[%d] = %s \n" %(region, i, edgesPerNjets[njet][key][1]))
                i += 1

        g.write("} \n")
        g.close()

                    
