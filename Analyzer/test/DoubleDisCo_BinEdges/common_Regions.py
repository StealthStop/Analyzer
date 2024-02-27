import ROOT
import os
import math
import ctypes
import numpy as np

ROOT.TH1.AddDirectory(False)

# ----------------------------------------------------------------------------------
# histBkg        : ROOT TH2 corresponding to background process
# histSig        : ROOT TH2 corresponding to signal process
# disc1Edge : fix the disc 1 edge with provided value 
# disc2Edge : fix the disc 2 edge with provided value 
# leftBoundary   : value corresponding to the left edge defining the "ABCD" region
# rightBoundary  : value corresponding to the right edge defining the "ABCD" region
# topBoundary    : value corresponding to the top edge defining the "ABCD" region
# bottomBoundary : value corresponding to the bottom edge defining the "ABCD" region
# ----------------------------------------------------------------------------------
class All_Regions:

    def __init__(self, hist=None, Sig=None, ttVar=None, disc1Edge=None, disc2Edge=None, leftBoundary=None, rightBoundary=None, topBoundary=None, bottomBoundary=None, fastMode=False, step=None, justEvents=False, QCDCRInfo=None, binStart=None, binEnd=None, **kwargs):

        self.hist           = hist
        self.Sig            = Sig
        self.ttVar          = ttVar
        self.disc1Edge      = round(float(disc1Edge), 2) if disc1Edge != None else None
        self.disc2Edge      = round(float(disc2Edge), 2) if disc2Edge != None else None
        self.leftBoundary   = round(float(leftBoundary), 2)   if leftBoundary   != None else None
        self.rightBoundary  = round(float(rightBoundary), 2)  if rightBoundary  != None else None
        self.topBoundary    = round(float(topBoundary), 2)    if topBoundary    != None else None
        self.bottomBoundary = round(float(bottomBoundary), 2) if bottomBoundary != None else None
        self.step           = step
        self.binStart       = binStart
        self.binEnd         = binEnd
        self.fastMode       = fastMode
        self.QCDCRInfo      = QCDCRInfo

        self.extraArgs = kwargs

        disc1EdgeStr = None
        disc2EdgeStr = None
        if disc1Edge != None:
            disc1EdgeStr = "%0.3f"%(float(disc1Edge))
        if disc2Edge != None:
            disc2EdgeStr = "%0.3f"%(float(disc2Edge))

        self.finalEdges = (disc1EdgeStr, disc2EdgeStr)
        self.quantities = {}

        for key in hist.keys():
            self.quantities[key] = {}

        # First calculate events counts for all possible choices of bin edges
        self.count_Events_inBinEdges()

        dilep = True
        # Need to create the QCD regions using the QCDCR prediction
        #if QCDCRInfo is not None and "QCD" in hist.keys():
        #    self.make_QCD_Regions(QCDCRInfo)
        #    dilep = False
        #if QCDCRInfo is not None:
        self.make_NonTT_Regions()

        # Then determine the "final" choice of bin edges for the region
        if not justEvents:
            self.get_nEvents_Quantities()

    def get_ttVar_Name(self):
        return self.ttVar

    # -------------------------------------
    # Significance calculation with only TT
    # -------------------------------------
    def cal_Significance(self, nSigEvents, nTTEvents, sys=0.3):
        
        if nSigEvents == 0.0:
            if nTTEvents <= 0.0:
                return -999.0
            return 0.0
        elif nTTEvents <= 0.0:
            return -999.0
    
        significance = nSigEvents / ( nTTEvents + (sys * nTTEvents)**2.0 )**0.5
        
        return significance
   
    # -------------------------------------------------
    # Significance calculation with only TT
    #   -- added non-closure to optimize the ABCD edges
    # -------------------------------------------------
    def cal_Significance_includingNonClosure(self, nSigEvents, nTTEvents, nonClosure=None, sys=0.3):
        
        if nSigEvents == 0.0:
            if nTTEvents <= 0.0:
                return -999.0
            return 0.0
        elif nTTEvents <= 0.0:
            return -999.0

        significance = nSigEvents / ( nTTEvents + (sys * nTTEvents)**2.0 + (nonClosure * nTTEvents)**2.0 )**0.5
        
        return significance
 
    # -------------------------------------
    # Significance calculation with only TT
    #   -- non simplified version
    # -------------------------------------
    def cal_Significance_nonSimplified(self, nSigEvents, nTTEvents, sys=0.3):
        
        if nSigEvents == 0.0:
            if nTTEvents <= 0.0:
                return -999.0
            return 0.0
        elif nTTEvents <= 0.0:
            return -999.0

        b = nTTEvents;     s = nSigEvents
        n = (s+b);     sigma = (b * sys)**2.0
       
        significance = (2 * ( 
                              ( n *  math.log( n * (b + sigma) / (b**2.0 + (n * sigma)) ) ) - 
                              ( (b**2.0 / sigma) * math.log( 1 + (sigma * s / (b * (b + sigma)) ) ) )  
                       ) )**0.5

        return significance

    # -------------------
    # Closure calculation
    # -------------------
    def cal_Closure(self, nEvents_A, nEvents_B, nEvents_C, nEvents_D, nEventsErr_A, nEventsErr_B, nEventsErr_C, nEventsErr_D):
    
        if nEvents_A <= 0.0 or nEvents_D <= 0.0:
            return -999.0, -999.0
        
        Closure = (nEvents_B * nEvents_C) / (nEvents_A * nEvents_D) 
        
        ClosureUnc = ( ( ( nEvents_C * nEventsErr_B ) / ( nEvents_A * nEvents_D) )**2.0 
                     + ( ( nEvents_B * nEventsErr_C ) / ( nEvents_A * nEvents_D) )**2.0 
                     + ( ( nEvents_B * nEvents_C * nEventsErr_A ) / ( nEvents_A**2.0 * nEvents_D ) )**2.0 
                     + ( ( nEvents_B * nEvents_C * nEventsErr_D ) / ( nEvents_A * nEvents_D**2.0 ) )**2.0 )**0.5
    
        return Closure, ClosureUnc

    # -----------------------
    # Non-Closure calculation
    # -----------------------
    def cal_NonClosure(self, nEvents_A, nEvents_B, nEvents_C, nEvents_D, nEventsErr_A, nEventsErr_B, nEventsErr_C, nEventsErr_D):
    
        if nEvents_A <= 0.0 or nEvents_D <= 0.0:
            return -999.0, -999.0
        
        nonClosure = abs(1.0 - ( (nEvents_B * nEvents_C) / (nEvents_A * nEvents_D) ) )
        
        nonClosureUnc = ( ( ( nEvents_C * nEventsErr_B ) / ( nEvents_A * nEvents_D) )**2.0 
                        + ( ( nEvents_B * nEventsErr_C ) / ( nEvents_A * nEvents_D) )**2.0 
                        + ( ( nEvents_B * nEvents_C * nEventsErr_A ) / ( nEvents_A**2.0 * nEvents_D ) )**2.0 
                        + ( ( nEvents_B * nEvents_C * nEventsErr_D ) / ( nEvents_A * nEvents_D**2.0 ) )**2.0 )**0.5
    
        return nonClosure, nonClosureUnc

    # ------------------------------
    # Closure Correction calculation
    # ------------------------------
    def cal_ClosureCorr(self, nEvents_A, nEvents_B, nEvents_C, nEvents_D, nEventsErr_A, nEventsErr_B, nEventsErr_C, nEventsErr_D):
 
        if nEvents_B <= 0.0 or nEvents_C <= 0.0:
            return -999.0, -999.0
        
        closureCorr = (nEvents_A * nEvents_D) / (nEvents_B * nEvents_C)
        
        closureCorrUnc = ( ( ( nEvents_A * nEventsErr_D ) / ( nEvents_B * nEvents_C) )**2.0 
                        + ( ( nEvents_D * nEventsErr_A ) / ( nEvents_B * nEvents_C) )**2.0 
                        + ( ( nEvents_D * nEvents_A * nEventsErr_A ) / ( nEvents_B**2.0 * nEvents_C ) )**2.0 
                        + ( ( nEvents_D * nEvents_A * nEventsErr_C ) / ( nEvents_B * nEvents_C**2.0 ) )**2.0 )**0.5
   
        return closureCorr, closureCorrUnc

    # ------------------------
    # Closure Pull calculation
    # ------------------------
    def cal_Pull(self, nEvents_A, nEvents_B, nEvents_C, nEvents_D, nEventsErr_A, nEventsErr_B, nEventsErr_C, nEventsErr_D):

        if nEvents_D <= 0.0 or nEventsErr_A <= 0.0:
            return -999.0, -999.0

        pull    = ( ( ((nEvents_B * nEvents_C) / nEvents_D) - nEvents_A ) / nEventsErr_A )
        pullUnc = 1.0

        return pull, pullUnc
    
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

        binStart = 1 if self.binStart == None else self.hist["TT"].GetXaxis().FindBin(float(self.binStart))
        binEnd = 1 if self.binEnd == None else self.hist["TT"].GetXaxis().FindBin(float(self.binEnd))

        totalXbins = self.hist["TT"].GetNbinsX()
        totalYbins = self.hist["TT"].GetNbinsY()

        lastXBin = totalXbins+1 if self.rightBoundary == None else self.hist["TT"].GetXaxis().FindBin(float(self.rightBoundary)) + 1
        lastYBin = totalYbins+1 if self.topBoundary   == None else self.hist["TT"].GetYaxis().FindBin(float(self.topBoundary))   + 1

        if lastXBin > totalXbins: lastXBin = totalXbins + 1
        if lastYBin > totalYbins: lastYBin = totalYbins + 1

        if self.step is not None:
            nXBins = range(binStart, binEnd, int(round(totalXbins * self.step)))
            nYBins = range(binStart, binEnd, int(round(totalYbins * self.step)))
        else:
            nXBins = range(firstXBin + 1, lastXBin)
            nYBins = range(firstYBin + 1, lastYBin)

        if not self.fastMode:
            for key, h1 in self.hist.items():
                for xBin in range(1, totalXbins+1):
                    xLowBinEdge = self.hist["TT"].GetXaxis().GetBinLowEdge(xBin)
                    xBinKey     = "%.3f"%(xLowBinEdge)

                    self.add("nEventsA", xBinKey, "1.00", (0.0, 0.0), key)
                    self.add("nEventsB", xBinKey, "1.00", (0.0, 0.0), key)
                    self.add("nEventsC", xBinKey, "1.00", (0.0, 0.0), key)
                    self.add("nEventsD", xBinKey, "1.00", (0.0, 0.0), key)

                    for yBin in range(1, totalYbins+1):
                        yLowBinEdge = self.hist["TT"].GetYaxis().GetBinLowEdge(yBin)
                        yBinKey     = "%.3f"%(yLowBinEdge)

                        self.add("nEventsA", xBinKey, yBinKey, (0.0, 0.0), key)
                        self.add("nEventsB", xBinKey, yBinKey, (0.0, 0.0), key)
                        self.add("nEventsC", xBinKey, yBinKey, (0.0, 0.0), key)
                        self.add("nEventsD", xBinKey, yBinKey, (0.0, 0.0), key)

                        self.add("nEventsA", "1.00", yBinKey, (0.0, 0.0), key)
                        self.add("nEventsB", "1.00", yBinKey, (0.0, 0.0), key)
                        self.add("nEventsC", "1.00", yBinKey, (0.0, 0.0), key)
                        self.add("nEventsD", "1.00", yBinKey, (0.0, 0.0), key)

                self.add("nEventsA", "1.00", "1.00", (0.0, 0.0), key)
                self.add("nEventsB", "1.00", "1.00", (0.0, 0.0), key)
                self.add("nEventsC", "1.00", "1.00", (0.0, 0.0), key)
                self.add("nEventsD", "1.00", "1.00", (0.0, 0.0), key)

        # count signal and background events and errors in bin edges
        for key, h1 in self.hist.items():

            #if "QCD" in key and "QCDCR" not in key and self.QCDCRInfo is not None:
            #    continue
            # loop over the x bins i.e. choice of disc 1 as an edge
            for xBin in nXBins:

                # Store disc 1 edge as string with three digits of accuracy for now
                xLowBinEdge = self.hist["TT"].GetXaxis().GetBinLowEdge(xBin)
                xBinKey     = "%.3f"%(xLowBinEdge)

                # Only care about actual choice of bin edges
                if self.fastMode and self.disc1Edge != None and (abs(self.disc1Edge - xLowBinEdge) >= 10e-3):
                    continue
                
                # For each choice of xBin (vertical divider in ABCD plane),
                # initialize counts for the four regions
                startOfScan = True
                nEvents_A = 0.0
                nEvents_B = 0.0
                nEvents_C = 0.0
                nEvents_D = 0.0

                # loop over the y bins
                for yBin in nYBins:

                    # Store disc 2 edge as string with three digits of accuracy for now
                    yLowBinEdge = self.hist["TT"].GetYaxis().GetBinLowEdge(yBin)
                    yBinKey     = "%.3f"%(yLowBinEdge)

                    # Only care about actual choice of bin edges
                    if self.fastMode and self.disc2Edge != None and (abs(self.disc2Edge - yLowBinEdge) >= 10e-3):
                        continue


                    nEventsErr_A = ctypes.c_double(0.0); nEventsErr_B = ctypes.c_double(0.0); nEventsErr_C = ctypes.c_double(0.0); nEventsErr_D = ctypes.c_double(0.0); nEventsErr_Tot = ctypes.c_double(0.0)

                    # last      | 
                    #        B  |  A
                    # yBin _____|_____
                    #           |
                    #        D  |  C
                    #           |
                    #    
                    # first   xBin   last
                    #if startOfScan:
                    nEvents_A = h1.IntegralAndError(xBin,      lastXBin, yBin,      lastYBin, nEventsErr_A)
                    nEvents_B = h1.IntegralAndError(firstXBin, xBin-1,   yBin,      lastYBin, nEventsErr_B)
                    nEvents_C = h1.IntegralAndError(xBin,      lastXBin, firstYBin, yBin-1,   nEventsErr_C)
                    nEvents_D = h1.IntegralAndError(firstXBin, xBin-1,   firstYBin, yBin-1,   nEventsErr_D)
                    nEvents_Tot = h1.IntegralAndError(firstXBin, lastXBin, firstYBin, lastYBin, nEventsErr_Tot)

                    #startOfScan = False

                    #else:
                    #    incrementBD = h1.IntegralAndError(firstXBin, xBin-1,   yBin-1, yBin-1, nEventsErr_A)
                    #    incrementAC = h1.IntegralAndError(xBin,      lastXBin, yBin-1, yBin-1, nEventsErr_B)
                   
                    #    nEvents_A -= incrementAC
                    #    nEvents_B -= incrementBD
                    #    nEvents_C += incrementAC
                    #    nEvents_D += incrementBD

                    self.add("nEventsA", xBinKey, yBinKey, (nEvents_A, nEventsErr_A.value), key)
                    self.add("nEventsB", xBinKey, yBinKey, (nEvents_B, nEventsErr_B.value), key)
                    self.add("nEventsC", xBinKey, yBinKey, (nEvents_C, nEventsErr_C.value), key)
                    self.add("nEventsD", xBinKey, yBinKey, (nEvents_D, nEventsErr_D.value), key)

                    # Needed to get the average weight for MC stats
                    entries = h1.GetEntries()

                    if entries > 0:
                        weight = h1.GetSumOfWeights() / entries
                    else:
                        weight = 1.0
                    
                    self.add("weight", xBinKey, yBinKey, (weight), key)
                    if nEvents_Tot > 0.0:
                        self.add("nEntriesA", xBinKey, yBinKey, (round(entries * (nEvents_A / nEvents_Tot))), key)
                        self.add("nEntriesB", xBinKey, yBinKey, (round(entries * (nEvents_B / nEvents_Tot))), key)
                        self.add("nEntriesC", xBinKey, yBinKey, (round(entries * (nEvents_C / nEvents_Tot))), key)
                        self.add("nEntriesD", xBinKey, yBinKey, (round(entries * (nEvents_D / nEvents_Tot))), key)
                    else:
                        self.add("nEntriesA", xBinKey, yBinKey, (0.0), key)
                        self.add("nEntriesB", xBinKey, yBinKey, (0.0), key)
                        self.add("nEntriesC", xBinKey, yBinKey, (0.0), key)
                        self.add("nEntriesD", xBinKey, yBinKey, (0.0), key)

        for h in self.hist.values():
            ROOT.SetOwnership(h, True)
            del h


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

            if type(self.ttVar) == str:
                nTTvarEvents_A, nTTvarEventsErr_A = self.get("nEventsA", disc1Key, disc2Key, self.ttVar)
                nTTvarEvents_B, nTTvarEventsErr_B = self.get("nEventsB", disc1Key, disc2Key, self.ttVar)
                nTTvarEvents_C, nTTvarEventsErr_C = self.get("nEventsC", disc1Key, disc2Key, self.ttVar)
                nTTvarEvents_D, nTTvarEventsErr_D = self.get("nEventsD", disc1Key, disc2Key, self.ttVar)
            else:
                nTTvarEvents_A = {}; nTTvarEvents_B = {}; nTTvarEvents_C = {}; nTTvarEvents_D = {}; 
                nTTvarEventsErr_A = {}; nTTvarEventsErr_B = {}; nTTvarEventsErr_C = {}; nTTvarEventsErr_D = {}; 

                for key in self.ttVar.keys():
                    nTTvarEvents_A[key], nTTvarEventsErr_A[key] = self.get("nEventsA", disc1Key, disc2Key, key)
                    nTTvarEvents_B[key], nTTvarEventsErr_B[key] = self.get("nEventsB", disc1Key, disc2Key, key)
                    nTTvarEvents_C[key], nTTvarEventsErr_C[key] = self.get("nEventsC", disc1Key, disc2Key, key)
                    nTTvarEvents_D[key], nTTvarEventsErr_D[key] = self.get("nEventsD", disc1Key, disc2Key, key)

            nSigEvents_A,   nSigEventsErr_A   = self.get("nEventsA", disc1Key, disc2Key, self.Sig) 
            nSigEvents_B,   nSigEventsErr_B   = self.get("nEventsB", disc1Key, disc2Key, self.Sig) 
            nSigEvents_C,   nSigEventsErr_C   = self.get("nEventsC", disc1Key, disc2Key, self.Sig) 
            nSigEvents_D,   nSigEventsErr_D   = self.get("nEventsD", disc1Key, disc2Key, self.Sig) 

            nDataEvents_A,  nDataEventsErr_A  = self.get("nEventsA", disc1Key, disc2Key, "Data")
            nDataEvents_B,  nDataEventsErr_B  = self.get("nEventsB", disc1Key, disc2Key, "Data")
            nDataEvents_C,  nDataEventsErr_C  = self.get("nEventsC", disc1Key, disc2Key, "Data")
            nDataEvents_D,  nDataEventsErr_D  = self.get("nEventsD", disc1Key, disc2Key, "Data")

            # Compute background subtracted data to give the TT in data for each of the four ABCD bins
            # To be used for computing the data closure correction
            nTTinDataEvents_A = nDataEvents_A - nNonTTEvents_A
            nTTinDataEvents_B = nDataEvents_B - nNonTTEvents_B
            nTTinDataEvents_C = nDataEvents_C - nNonTTEvents_C
            nTTinDataEvents_D = nDataEvents_D - nNonTTEvents_D

            nTTinDataEventsErr_A = (nDataEventsErr_A**2.0 + nNonTTEventsErr_A**2.0)**0.5
            nTTinDataEventsErr_B = (nDataEventsErr_B**2.0 + nNonTTEventsErr_B**2.0)**0.5
            nTTinDataEventsErr_C = (nDataEventsErr_C**2.0 + nNonTTEventsErr_C**2.0)**0.5
            nTTinDataEventsErr_D = (nDataEventsErr_D**2.0 + nNonTTEventsErr_D**2.0)**0.5

            self.add("nEventsA", disc1Key, disc2Key, (nTTinDataEvents_A, nTTinDataEventsErr_A), "TTinData")
            self.add("nEventsB", disc1Key, disc2Key, (nTTinDataEvents_B, nTTinDataEventsErr_B), "TTinData")
            self.add("nEventsC", disc1Key, disc2Key, (nTTinDataEvents_C, nTTinDataEventsErr_C), "TTinData")
            self.add("nEventsD", disc1Key, disc2Key, (nTTinDataEvents_D, nTTinDataEventsErr_D), "TTinData")

            # Region by region signal and bkg fractions for 1D-2D plots / for only TT !!!
            nTot_SigTT_A = nSigEvents_A + nTTEvents_A
            nTot_SigTT_B = nSigEvents_B + nTTEvents_B
            nTot_SigTT_C = nSigEvents_C + nTTEvents_C
            nTot_SigTT_D = nSigEvents_D + nTTEvents_D

            sigFracsA    = -999.0; sigFracsB    = -999.0; sigFracsC    = -999.0; sigFracsD    = -999.0
            sigFracsErrA = -999.0; sigFracsErrB = -999.0; sigFracsErrC = -999.0; sigFracsErrD = -999.0
            bkgFracsA    = -999.0; bkgFracsB    = -999.0; bkgFracsC    = -999.0; bkgFracsD    = -999.0
            bkgFracsErrA = -999.0; bkgFracsErrB = -999.0; bkgFracsErrC = -999.0; bkgFracsErrD = -999.0
 
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
            Closure_TT               = -999.0; ClosureUnc_TT       = -999.0
            Closure_NonTT            = -999.0; ClosureUnc_NonTT    = -999.0
            Closure_TTvar            = -999.0; ClosureUnc_TTvar    = -999.0
            Closure_Data             = -999.0; ClosureUnc_Data     = -999.0
            Closure_TTinData         = -999.0; ClosureUnc_TTinData = -999.0

            nonClosure_TT            = -999.0; nonClosureUnc_TT       = -999.0
            nonClosure_NonTT         = -999.0; nonClosureUnc_NonTT    = -999.0
            nonClosure_TTvar         = -999.0; nonClosureUnc_TTvar    = -999.0
            nonClosure_Data          = -999.0; nonClosureUnc_Data     = -999.0
            nonClosure_TTinData      = -999.0; nonClosureUnc_TTinData = -999.0

            pull_TT                  = -999.0; pullUnc_TT       = -999.0   
            pull_NonTT               = -999.0; pullUnc_NonTT    = -999.0
            pull_TTvar               = -999.0; pullUnc_TTvar    = -999.0
            pull_Data                = -999.0; pullUnc_Data     = -999.0
            pull_TTinData            = -999.0; pullUnc_TTinData = -999.0

            closureCorr_TT           = -999.0; closureCorrUnc_TT           = -999.0
            closureCorr_NonTT        = -999.0; closureCorrUnc_NonTT        = -999.0
            closureCorr_TTvar        = -999.0; closureCorrUnc_TTvar        = -999.0
            closureCorr_Data         = -999.0; closureCorrUnc_Data         = -999.0
            closureCorr_TTinData     = -999.0; closureCorrUnc_TTinData     = -999.0
            closureCorr_TTinDataVsTT = -999.0; closureCorrUnc_TTinDataVsTT = -999.0

            #
            Closure_TT,       ClosureUnc_TT       = self.cal_Closure(nTTEvents_A,    nTTEvents_B,    nTTEvents_C,    nTTEvents_D,    nTTEventsErr_A,    nTTEventsErr_B,    nTTEventsErr_C,    nTTEventsErr_D   )
            Closure_NonTT,    ClosureUnc_NonTT    = self.cal_Closure(nNonTTEvents_A, nNonTTEvents_B, nNonTTEvents_C, nNonTTEvents_D, nNonTTEventsErr_A, nNonTTEventsErr_B, nNonTTEventsErr_C, nNonTTEventsErr_D)
            if type(self.ttVar) == str:
                Closure_TTvar,    ClosureUnc_TTvar    = self.cal_Closure(nTTvarEvents_A, nTTvarEvents_B, nTTvarEvents_C, nTTvarEvents_D, nTTvarEventsErr_A, nTTvarEventsErr_B, nTTvarEventsErr_C, nTTvarEventsErr_D)
            else:
                Closure_TTvar = {}; ClosureUnc_TTvar = {}
                for key in self.ttVar.keys():
                    Closure_TTvar[key], ClosureUnc_TTvar[key] = self.cal_Closure(nTTvarEvents_A[key], nTTvarEvents_B[key], nTTvarEvents_C[key], nTTvarEvents_D[key], nTTvarEventsErr_A[key], nTTvarEventsErr_B[key], nTTvarEventsErr_C[key], nTTvarEventsErr_D[key])

            Closure_Data,     ClosureUnc_Data     = self.cal_Closure(nDataEvents_A,  nDataEvents_B,  nDataEvents_C,  nDataEvents_D,  nDataEventsErr_A,  nDataEventsErr_B,  nDataEventsErr_C,  nDataEventsErr_D )

            Closure_TTinData, ClosureUnc_TTinData = self.cal_Closure(nTTinDataEvents_A,  nTTinDataEvents_B,  nTTinDataEvents_C,  nTTinDataEvents_D,  nTTinDataEventsErr_A,  nTTinDataEventsErr_B,  nTTinDataEventsErr_C,  nTTinDataEventsErr_D )

            nonClosure_TT,     nonClosureUnc_TT     = self.cal_NonClosure(nTTEvents_A,    nTTEvents_B,    nTTEvents_C,    nTTEvents_D,    nTTEventsErr_A,    nTTEventsErr_B,    nTTEventsErr_C,    nTTEventsErr_D   )
            nonClosure_NonTT,  nonClosureUnc_NonTT  = self.cal_NonClosure(nNonTTEvents_A, nNonTTEvents_B, nNonTTEvents_C, nNonTTEvents_D, nNonTTEventsErr_A, nNonTTEventsErr_B, nNonTTEventsErr_C, nNonTTEventsErr_D)
            if type(self.ttVar) == str:
                nonClosure_TTvar,  nonClosureUnc_TTvar  = self.cal_NonClosure(nTTvarEvents_A, nTTvarEvents_B, nTTvarEvents_C, nTTvarEvents_D, nTTvarEventsErr_A, nTTvarEventsErr_B, nTTvarEventsErr_C, nTTvarEventsErr_D)
            else:
                nonClosure_TTvar = {};  nonClosureUnc_TTvar = {}
                for key in self.ttVar.keys():
                    nonClosure_TTvar[key],  nonClosureUnc_TTvar[key]  = self.cal_NonClosure(nTTvarEvents_A[key], nTTvarEvents_B[key], nTTvarEvents_C[key], nTTvarEvents_D[key], nTTvarEventsErr_A[key], nTTvarEventsErr_B[key], nTTvarEventsErr_C[key], nTTvarEventsErr_D[key])
                    
            nonClosure_Data,   nonClosureUnc_Data   = self.cal_NonClosure(nDataEvents_A,  nDataEvents_B,  nDataEvents_C,  nDataEvents_D,  nDataEventsErr_A,  nDataEventsErr_B,  nDataEventsErr_C,  nDataEventsErr_D )
            nonClosure_TTinData, nonClosureUnc_TTinData = self.cal_NonClosure(nTTinDataEvents_A,  nTTinDataEvents_B,  nTTinDataEvents_C,  nTTinDataEvents_D,  nTTinDataEventsErr_A,  nTTinDataEventsErr_B,  nTTinDataEventsErr_C,  nTTinDataEventsErr_D )
            
            pull_TT,    pullUnc_TT    = self.cal_Pull(nTTEvents_A,    nTTEvents_B,    nTTEvents_C,    nTTEvents_D,    nTTEventsErr_A,    nTTEventsErr_B,    nTTEventsErr_C,    nTTEventsErr_D   )
            pull_NonTT, pullUnc_NonTT = self.cal_Pull(nNonTTEvents_A, nNonTTEvents_B, nNonTTEvents_C, nNonTTEvents_D, nNonTTEventsErr_A, nNonTTEventsErr_B, nNonTTEventsErr_C, nNonTTEventsErr_D)
            if type(self.ttVar) == str:
                pull_TTvar, pullUnc_TTvar = self.cal_Pull(nTTvarEvents_A, nTTvarEvents_B, nTTvarEvents_C, nTTvarEvents_D, nTTvarEventsErr_A, nTTvarEventsErr_B, nTTvarEventsErr_C, nTTvarEventsErr_D)
            else:
                pull_TTvar = {}; pullUnc_TTvar = {}
                for key in self.ttVar.keys():
                    pull_TTvar[key], pullUnc_TTvar[key] = self.cal_Pull(nTTvarEvents_A[key], nTTvarEvents_B[key], nTTvarEvents_C[key], nTTvarEvents_D[key], nTTvarEventsErr_A[key], nTTvarEventsErr_B[key], nTTvarEventsErr_C[key], nTTvarEventsErr_D[key])
            pull_Data,  pullUnc_Data  = self.cal_Pull(nDataEvents_A,  nDataEvents_B,  nDataEvents_C,  nDataEvents_D,  nDataEventsErr_A,  nDataEventsErr_B,  nDataEventsErr_C,  nDataEventsErr_D )
            pull_TTinData,  pullUnc_TTinData  = self.cal_Pull(nTTinDataEvents_A,  nTTinDataEvents_B,  nTTinDataEvents_C,  nTTinDataEvents_D,  nTTinDataEventsErr_A,  nTTinDataEventsErr_B,  nTTinDataEventsErr_C,  nTTinDataEventsErr_D )

            closureCorr_TT,    closureCorrUnc_TT    = self.cal_ClosureCorr(nTTEvents_A,    nTTEvents_B,    nTTEvents_C,    nTTEvents_D,    nTTEventsErr_A,    nTTEventsErr_B,    nTTEventsErr_C,    nTTEventsErr_D   )
            closureCorr_Data,  closureCorrUnc_Data  = self.cal_ClosureCorr(nDataEvents_A,  nDataEvents_B,  nDataEvents_C,  nDataEvents_D,  nDataEventsErr_A,  nDataEventsErr_B,  nDataEventsErr_C,  nDataEventsErr_D )
            closureCorr_NonTT, closureCorrUnc_NonTT = self.cal_ClosureCorr(nNonTTEvents_A, nNonTTEvents_B, nNonTTEvents_C, nNonTTEvents_D, nNonTTEventsErr_A, nNonTTEventsErr_B, nNonTTEventsErr_C, nNonTTEventsErr_D   )
            if type(self.ttVar) == str:
                closureCorr_TTvar, closureCorrUnc_TTvar = self.cal_ClosureCorr(nTTvarEvents_A, nTTvarEvents_B, nTTvarEvents_C, nTTvarEvents_D, nTTvarEventsErr_A, nTTvarEventsErr_B, nTTvarEventsErr_C, nTTvarEventsErr_D )
            else:
                closureCorr_TTvar = {}; closureCorrUnc_TTvar = {} 
                for key in self.ttVar.keys():
                    closureCorr_TTvar[key], closureCorrUnc_TTvar[key] = self.cal_ClosureCorr(nTTvarEvents_A[key], nTTvarEvents_B[key], nTTvarEvents_C[key], nTTvarEvents_D[key], nTTvarEventsErr_A[key], nTTvarEventsErr_B[key], nTTvarEventsErr_C[key], nTTvarEventsErr_D[key] )
                
            closureCorr_TTinData, closureCorrUnc_TTinData = self.cal_ClosureCorr(nTTinDataEvents_A, nTTinDataEvents_B, nTTinDataEvents_C, nTTinDataEvents_D, nTTinDataEventsErr_A, nTTinDataEventsErr_B, nTTinDataEventsErr_C, nTTinDataEventsErr_D )

            #
            self.add("Closure",    disc1Key, disc2Key, (Closure_TT,       ClosureUnc_TT) ,       "TT"      )
            self.add("Closure",    disc1Key, disc2Key, (Closure_NonTT,    ClosureUnc_NonTT),     "NonTT"   )
            if type(self.ttVar) == str:
                self.add("Closure",    disc1Key, disc2Key, (Closure_TTvar,    ClosureUnc_TTvar),     self.ttVar)
            else:
                for key in self.ttVar.keys():
                    self.add("Closure",    disc1Key, disc2Key, (Closure_TTvar[key],    ClosureUnc_TTvar[key]),     key)
            self.add("Closure",    disc1Key, disc2Key, (Closure_Data,     ClosureUnc_Data ),     "Data"    )
            self.add("Closure",    disc1Key, disc2Key, (Closure_TTinData, ClosureUnc_TTinData ), "TTinData")
            
            self.add("nonClosure", disc1Key, disc2Key, (nonClosure_TT,       nonClosureUnc_TT) ,       "TT"      )
            self.add("nonClosure", disc1Key, disc2Key, (nonClosure_NonTT,    nonClosureUnc_NonTT),     "NonTT"   )
            if type(self.ttVar) == str:
                self.add("nonClosure", disc1Key, disc2Key, (nonClosure_TTvar,    nonClosureUnc_TTvar),     self.ttVar)
            else:
                for key in self.ttVar.keys():
                    self.add("nonClosure", disc1Key, disc2Key, (nonClosure_TTvar[key],    nonClosureUnc_TTvar[key]),     key)
            self.add("nonClosure", disc1Key, disc2Key, (nonClosure_Data,     nonClosureUnc_Data ),     "Data"    )
            self.add("nonClosure", disc1Key, disc2Key, (nonClosure_TTinData, nonClosureUnc_TTinData ), "TTinData")

            self.add("pull",       disc1Key, disc2Key, (pull_TT,       pullUnc_TT),       "TT"      )
            self.add("pull",       disc1Key, disc2Key, (pull_NonTT,    pullUnc_NonTT),    "NonTT"   )
            if type(self.ttVar) == str:
                self.add("pull",       disc1Key, disc2Key, (pull_TTvar,    pullUnc_TTvar),    self.ttVar)
            else:
                for key in self.ttVar.keys():
                    self.add("pull", disc1Key, disc2Key, (nonClosure_TTvar[key],    nonClosureUnc_TTvar[key]),     key)
            self.add("pull",       disc1Key, disc2Key, (pull_Data,     pullUnc_Data),     "Data"    )
            self.add("pull",       disc1Key, disc2Key, (pull_TTinData, pullUnc_TTinData), "TTinData")

            # MC closure correction factor     
            self.add("closureCorr", disc1Key, disc2Key, (closureCorr_TT,       closureCorrUnc_TT) ,       "TT"      )
            self.add("closureCorr", disc1Key, disc2Key, (closureCorr_NonTT,    closureCorrUnc_NonTT) ,    "NonTT"   )
            if type(self.ttVar) == str:
                self.add("closureCorr", disc1Key, disc2Key, (closureCorr_TTvar,    closureCorrUnc_TTvar ),    self.ttVar)
            else:
                for key in self.ttVar.keys():
                    self.add("closureCorr", disc1Key, disc2Key, (nonClosure_TTvar[key],    nonClosureUnc_TTvar[key]),     key)
            self.add("closureCorr", disc1Key, disc2Key, (closureCorr_Data,     closureCorrUnc_Data ),     "Data"    )
            self.add("closureCorr", disc1Key, disc2Key, (closureCorr_TTinData, closureCorrUnc_TTinData ), "TTinData")

            # ------------------------------------
            # MC correction factor ratio: TT/TTvar
            #   -- in MC level
            # ------------------------------------
            if type(self.ttVar) == str:
                temp_closureCorr_TTvar = closureCorr_TTvar
                temp_closureCorr_TT    = closureCorr_TT

                if (closureCorr_TTvar == 0.0): 
                    temp_closureCorr_TTvar = 10e-10
                    
                if (closureCorr_TT == 0.0):
                    temp_closureCorr_TT = 10e-10

                MCcorrRatio_MC            = (closureCorr_TTvar / temp_closureCorr_TT)
                closureCorrUnc_TT_withSys = abs(1.0 - MCcorrRatio_MC)

                MCcorrRatio_MC_Unc         = MCcorrRatio_MC * ( (closureCorrUnc_TTvar / temp_closureCorr_TTvar)**2.0 + (closureCorrUnc_TT / temp_closureCorr_TT)**2.0 )**0.5
                MCcorrRatio_MC_Unc_withSys = MCcorrRatio_MC * ( (closureCorrUnc_TTvar / temp_closureCorr_TTvar)**2.0 + ((closureCorrUnc_TT**2.0 + closureCorrUnc_TT_withSys**2.0)**0.5 / temp_closureCorr_TT)**2.0 )**0.5

            else:
                MCcorrRatio_MC = {}; MCcorrRatio_MC_Unc = {}; MCcorrRatio_MC_Unc_withSys = {}
                for key in self.ttVar.keys():
                    temp_closureCorr_TTvar = closureCorr_TTvar[key]
                    temp_closureCorr_TT    = closureCorr_TT

                    if (closureCorr_TTvar[key] == 0.0): 
                        temp_closureCorr_TTvar = 10e-10
                        
                    if (closureCorr_TT == 0.0):
                        temp_closureCorr_TT = 10e-10

                    MCcorrRatio_MC[key]            = (closureCorr_TTvar[key] / temp_closureCorr_TT)
                    closureCorrUnc_TT_withSys = abs(1.0 - MCcorrRatio_MC[key])

                    MCcorrRatio_MC_Unc[key]         = MCcorrRatio_MC[key] * ( (closureCorrUnc_TTvar[key] / temp_closureCorr_TTvar)**2.0 + (closureCorrUnc_TT / temp_closureCorr_TT)**2.0 )**0.5
                    MCcorrRatio_MC_Unc_withSys[key] = MCcorrRatio_MC[key] * ( (closureCorrUnc_TTvar[key] / temp_closureCorr_TTvar)**2.0 + ((closureCorrUnc_TT**2.0 + closureCorrUnc_TT_withSys**2.0)**0.5 / temp_closureCorr_TT)**2.0 )**0.5
                

            if type(self.ttVar) == str:
                self.add("MCcorrRatio_MC", disc1Key, disc2Key, (MCcorrRatio_MC, MCcorrRatio_MC_Unc, MCcorrRatio_MC_Unc_withSys), self.ttVar) 
            else:
                for key in self.ttVar.keys():
                    self.add("MCcorrRatio_MC", disc1Key, disc2Key, (MCcorrRatio_MC, MCcorrRatio_MC_Unc, MCcorrRatio_MC_Unc_withSys), key) 

            # ---------------------------------------------------------
            # MC corrected Data Closure for TT
            #   -- using MC correction factor to calculate Data Closure
            # ---------------------------------------------------------
            CorrectedDataClosure             = (Closure_TTinData * closureCorr_TT)
            CorrectedDataClosure_Unc         = math.sqrt((Closure_TTinData * closureCorrUnc_TT)**2.0 + (ClosureUnc_TTinData * closureCorr_TT)**2.0)
            CorrectedDataClosure_Unc_withSys = math.sqrt((Closure_TTinData * closureCorrUnc_TT_withSys)**2.0 + (ClosureUnc_TTinData * closureCorr_TT)**2.0)
            self.add("CorrectedDataClosure", disc1Key, disc2Key, (CorrectedDataClosure, CorrectedDataClosure_Unc, CorrectedDataClosure_Unc_withSys), "TTinData")

            # ----------------------------------------------------
            # MC corrected Data Closure for TTvar
            # using MC correction factor to calculate Data Closure
            # ----------------------------------------------------
            if type(self.ttVar) == str:
                ttVar_CorrectedDataClosure     = (Closure_TTinData * closureCorr_TTvar)
                ttVar_CorrectedDataClosure_Unc = math.sqrt((Closure_TTinData * closureCorrUnc_TTvar)**2.0 + (ClosureUnc_TTinData * closureCorr_TTvar)**2.0)
                self.add("ttVar_CorrectedDataClosure", disc1Key, disc2Key, (ttVar_CorrectedDataClosure, ttVar_CorrectedDataClosure_Unc), "TTinData")
            else:
                for key in self.ttVar.keys():
                    ttVar_CorrectedDataClosure     = (Closure_TTinData * closureCorr_TTvar[key])
                    ttVar_CorrectedDataClosure_Unc = math.sqrt((Closure_TTinData * closureCorrUnc_TTvar[key])**2.0 + (ClosureUnc_TTinData * closureCorr_TTvar[key])**2.0)
                    self.add("ttVar_CorrectedDataClosure", disc1Key, disc2Key, (ttVar_CorrectedDataClosure, ttVar_CorrectedDataClosure_Unc), "TTinData"+key)
                    

            # ----------------------------------------------------
            # significance for optimization metric for only TT !!! 
            # significance, significanceUnc for 2D plots
            # ----------------------------------------------------
            significance_TT = 0.0; significanceUnc_TT = 0.0; tempOptMetric = 999.0

            significance_TT = self.cal_Significance(nSigEvents_A, nTTEvents_A)**2.0
            significance_TT = significance_TT**0.5

            if nTTEvents_A > 0.0:
                D = ( nTTEvents_A + (bkgNormUnc * nTTEvents_A)**2.0 )**0.5
                significanceUnc_TT = ( (nSigEventsErr_A / D )**2.0 + ( (nSigEvents_A * nTTEventsErr_A) * (1 + (2.0 * nTTEvents_A * bkgNormUnc**2.) / 2 * D**3.0) ) )**0.5
 
            self.add("significance", disc1Key, disc2Key, (significance_TT, significanceUnc_TT), "TT") 

            # significance for optimizing ABCD edges
            # this one including also non-closure
            significance_includingNonClosure = 0.0; significanceUnc_includingNonClosure = 0.0
            significance_includingNonClosure = self.cal_Significance_includingNonClosure(nSigEvents_A, nTTEvents_A, nonClosure_TT)**2.0
            significance_includingNonClosure = significance_includingNonClosure**0.5

            if nTTEvents_A > 0.0:
                significanceUnc_includingNonClosure = ( ( nSigEventsErr_A / (nTTEvents_A + (bkgNormUnc * nTTEvents_A)**2.0 + (nonClosure_TT * nTTEvents_A)**2.0)**0.5 )**2.0
                                                    + ( ( nSigEvents_A * nTTEventsErr_A * (2.0 * nTTEvents_A * nonClosure_TT**2.0 + 2.0 * bkgNormUnc**2.0 * nTTEvents_A + 1) ) / ( nTTEvents_A + (bkgNormUnc * nTTEvents_A)**2.0 + (nonClosure_TT * nTTEvents_A)**2.0 )**1.5 )**2.0
                                                    + ( ( nTTEvents_A**2.0 * nonClosure_TT * nSigEvents_A * nonClosureUnc_TT) / ( nTTEvents_A * ( nTTEvents_A * (nonClosure_TT**2.0 + bkgNormUnc**2.0) + 1) )**1.5 )**2.0 )**0.5

            self.add("significance_includingNonClosure", disc1Key, disc2Key, (significance_includingNonClosure,significanceUnc_includingNonClosure), "TT")

            # this one non-simplified version
            significance_nonSimplified = 0.0; significanceUnc_nonSimplified = 0.0
            significance_nonSimplified = self.cal_Significance_nonSimplified(nSigEvents_A, nTTEvents_A)**2.0
            significance_nonSimplified = significance_nonSimplified**0.5
            self.add("significance_nonSimplified", disc1Key, disc2Key, (significance_nonSimplified,significanceUnc_nonSimplified), "TT")


            # use this statement for cdGH regions if  fixed disc1 edge (vertivcal edge)
            if self.disc1Edge != None and abs(self.disc1Edge - float(disc1Key)) > 0.01: continue

            # use this statement for bdEF regions if fixed disc2 edge (horizontal edge)
            if self.disc2Edge != None and abs(self.disc2Edge - float(disc2Key)) > 0.01: continue

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

    def getQCDCRValues(self):
        
        # Surrogate function for calculating the QCDCR transfer factor and propagating to QCD region event yields

        QCDCR_dict = {"nQCDEvents_CR": [], "nQCDEvents_SR": []}

        for disc1Key, disc2Key in self.get("edges",None,None,"TT"):

            # Getting the TF per bin edge
            nQCDEvents_A_CR,    nQCDEventsErr_A_CR    = self.get("nEventsA", disc1Key, disc2Key, "QCD_QCDCR")
            nQCDEvents_B_CR,    nQCDEventsErr_B_CR    = self.get("nEventsB", disc1Key, disc2Key, "QCD_QCDCR")
            nQCDEvents_C_CR,    nQCDEventsErr_C_CR    = self.get("nEventsC", disc1Key, disc2Key, "QCD_QCDCR")
            nQCDEvents_D_CR,    nQCDEventsErr_D_CR    = self.get("nEventsD", disc1Key, disc2Key, "QCD_QCDCR")

            QCDCR_dict["nQCDEvents_CR"] = [nQCDEvents_A_CR, nQCDEvents_B_CR, nQCDEvents_C_CR, nQCDEvents_D_CR]

            nQCDEvents_A_SR,    nQCDEventsErr_A_SR    = self.get("nEventsA", disc1Key, disc2Key, "QCD")
            nQCDEvents_B_SR,    nQCDEventsErr_B_SR    = self.get("nEventsB", disc1Key, disc2Key, "QCD")
            nQCDEvents_C_SR,    nQCDEventsErr_C_SR    = self.get("nEventsC", disc1Key, disc2Key, "QCD")
            nQCDEvents_D_SR,    nQCDEventsErr_D_SR    = self.get("nEventsD", disc1Key, disc2Key, "QCD")

            QCDCR_dict["nQCDEvents_SR"] = [nQCDEvents_A_SR, nQCDEvents_B_SR, nQCDEvents_C_SR, nQCDEvents_D_SR]

        return QCDCR_dict
    
    def make_QCD_Regions(self, QCDCRInfo):
        
        # Distill down the QCD info into final QCDCR SR prediction

        total_CR = [0.0, 0.0, 0.0, 0.0]
        total_SR = [0.0, 0.0, 0.0, 0.0]

        TFs = []

        if self.disc1Edge is not None and self.disc2Edge is not None:

            for njet in QCDCRInfo.keys():

                temp_CR = QCDCRInfo[njet]["nQCDEvents_CR"]
                temp_SR = QCDCRInfo[njet]["nQCDEvents_SR"]
           
                for i in range(4):
                    
                    total_CR[i] += temp_CR[i]
                    total_SR[i] += temp_SR[i]
                     
            for i in range(len(total_CR)):
                
                TFs.append(total_SR[i]/total_CR[i])

            for disc1Key, disc2Key in self.get("edges",None,None,"TT"):

                # Getting number of events in the QCDCR in data - other MC backgrounds
                nDataEvents_A,    nDataEventsErr_A    = self.get("nEventsA", disc1Key, disc2Key, "Data_QCDCR")
                nDataEvents_B,    nDataEventsErr_B    = self.get("nEventsB", disc1Key, disc2Key, "Data_QCDCR")
                nDataEvents_C,    nDataEventsErr_C    = self.get("nEventsC", disc1Key, disc2Key, "Data_QCDCR")
                nDataEvents_D,    nDataEventsErr_D    = self.get("nEventsD", disc1Key, disc2Key, "Data_QCDCR")
            
                nTTXEvents_A,    nTTXEventsErr_A    = self.get("nEventsA", disc1Key, disc2Key, "TTX_QCDCR")
                nTTXEvents_B,    nTTXEventsErr_B    = self.get("nEventsB", disc1Key, disc2Key, "TTX_QCDCR")
                nTTXEvents_C,    nTTXEventsErr_C    = self.get("nEventsC", disc1Key, disc2Key, "TTX_QCDCR")
                nTTXEvents_D,    nTTXEventsErr_D    = self.get("nEventsD", disc1Key, disc2Key, "TTX_QCDCR")
            
                nOtherEvents_A,    nOtherEventsErr_A    = self.get("nEventsA", disc1Key, disc2Key, "BG_OTHER_QCDCR")
                nOtherEvents_B,    nOtherEventsErr_B    = self.get("nEventsB", disc1Key, disc2Key, "BG_OTHER_QCDCR")
                nOtherEvents_C,    nOtherEventsErr_C    = self.get("nEventsC", disc1Key, disc2Key, "BG_OTHER_QCDCR")
                nOtherEvents_D,    nOtherEventsErr_D    = self.get("nEventsD", disc1Key, disc2Key, "BG_OTHER_QCDCR")
            
                nTTEvents_A,    nTTEventsErr_A    = self.get("nEventsA", disc1Key, disc2Key, "TT_QCDCR")
                nTTEvents_B,    nTTEventsErr_B    = self.get("nEventsB", disc1Key, disc2Key, "TT_QCDCR")
                nTTEvents_C,    nTTEventsErr_C    = self.get("nEventsC", disc1Key, disc2Key, "TT_QCDCR")
                nTTEvents_D,    nTTEventsErr_D    = self.get("nEventsD", disc1Key, disc2Key, "TT_QCDCR")
           
                nQCDinData_A = nDataEvents_A - nTTXEvents_A - nOtherEvents_A - nTTEvents_A
                nQCDinData_B = nDataEvents_B - nTTXEvents_B - nOtherEvents_B - nTTEvents_B
                nQCDinData_C = nDataEvents_C - nTTXEvents_C - nOtherEvents_C - nTTEvents_C
                nQCDinData_D = nDataEvents_D - nTTXEvents_D - nOtherEvents_D - nTTEvents_D

                nQCDinDataErr_A = math.sqrt(nDataEventsErr_A ** 2 + nTTXEventsErr_A ** 2 + nOtherEventsErr_A ** 2 + nTTEventsErr_A ** 2) * TFs[0]
                nQCDinDataErr_B = math.sqrt(nDataEventsErr_B ** 2 + nTTXEventsErr_B ** 2 + nOtherEventsErr_B ** 2 + nTTEventsErr_B ** 2) * TFs[1]
                nQCDinDataErr_C = math.sqrt(nDataEventsErr_C ** 2 + nTTXEventsErr_C ** 2 + nOtherEventsErr_C ** 2 + nTTEventsErr_C ** 2) * TFs[2]
                nQCDinDataErr_D = math.sqrt(nDataEventsErr_D ** 2 + nTTXEventsErr_D ** 2 + nOtherEventsErr_D ** 2 + nTTEventsErr_D ** 2) * TFs[3]

                nQCDFinal_A = nQCDinData_A * TFs[0]
                nQCDFinal_B = nQCDinData_B * TFs[1]
                nQCDFinal_C = nQCDinData_C * TFs[2]
                nQCDFinal_D = nQCDinData_D * TFs[3]
     
                self.add("nEventsA", disc1Key, disc2Key, (nQCDFinal_A, nQCDinDataErr_A), "QCDFinal")
                self.add("nEventsB", disc1Key, disc2Key, (nQCDFinal_B, nQCDinDataErr_B), "QCDFinal")
                self.add("nEventsC", disc1Key, disc2Key, (nQCDFinal_C, nQCDinDataErr_C), "QCDFinal")
                self.add("nEventsD", disc1Key, disc2Key, (nQCDFinal_D, nQCDinDataErr_D), "QCDFinal")

                self.add("nEventsA", disc1Key, disc2Key, (nQCDinData_A), "QCDCR")
                self.add("nEventsB", disc1Key, disc2Key, (nQCDinData_B), "QCDCR")
                self.add("nEventsC", disc1Key, disc2Key, (nQCDinData_C), "QCDCR")
                self.add("nEventsD", disc1Key, disc2Key, (nQCDinData_D), "QCDCR")

                self.add("QCDTF_A", disc1Key, disc2Key, (TFs[0]), "QCD")
                self.add("QCDTF_B", disc1Key, disc2Key, (TFs[1]), "QCD")
                self.add("QCDTF_C", disc1Key, disc2Key, (TFs[2]), "QCD")
                self.add("QCDTF_D", disc1Key, disc2Key, (TFs[3]), "QCD")
        else:


            for disc1Key, disc2Key in self.get("edges",None,None,"TT"):

                TFs = []
                total_CR = [0.0, 0.0, 0.0, 0.0]
                total_SR = [0.0, 0.0, 0.0, 0.0]

                for njet in QCDCRInfo["({},{})".format(disc1Key, disc2Key)].keys():

                    temp_CR = QCDCRInfo["({},{})".format(disc1Key, disc2Key)][njet]["nQCDEvents_CR"]
                    temp_SR = QCDCRInfo["({},{})".format(disc1Key, disc2Key)][njet]["nQCDEvents_SR"]
               
                    for i in range(4):
                        
                        total_CR[i] += temp_CR[i]
                        total_SR[i] += temp_SR[i]
                         
                for i in range(len(total_CR)):
                
                    TFs.append(total_SR[i]/total_CR[i])

                # Getting number of events in the QCDCR in data - other MC backgrounds
                nDataEvents_A,    nDataEventsErr_A    = self.get("nEventsA", disc1Key, disc2Key, "Data_QCDCR")
                nDataEvents_B,    nDataEventsErr_B    = self.get("nEventsB", disc1Key, disc2Key, "Data_QCDCR")
                nDataEvents_C,    nDataEventsErr_C    = self.get("nEventsC", disc1Key, disc2Key, "Data_QCDCR")
                nDataEvents_D,    nDataEventsErr_D    = self.get("nEventsD", disc1Key, disc2Key, "Data_QCDCR")
            
                nTTXEvents_A,    nTTXEventsErr_A    = self.get("nEventsA", disc1Key, disc2Key, "TTX_QCDCR")
                nTTXEvents_B,    nTTXEventsErr_B    = self.get("nEventsB", disc1Key, disc2Key, "TTX_QCDCR")
                nTTXEvents_C,    nTTXEventsErr_C    = self.get("nEventsC", disc1Key, disc2Key, "TTX_QCDCR")
                nTTXEvents_D,    nTTXEventsErr_D    = self.get("nEventsD", disc1Key, disc2Key, "TTX_QCDCR")
            
                nOtherEvents_A,    nOtherEventsErr_A    = self.get("nEventsA", disc1Key, disc2Key, "BG_OTHER_QCDCR")
                nOtherEvents_B,    nOtherEventsErr_B    = self.get("nEventsB", disc1Key, disc2Key, "BG_OTHER_QCDCR")
                nOtherEvents_C,    nOtherEventsErr_C    = self.get("nEventsC", disc1Key, disc2Key, "BG_OTHER_QCDCR")
                nOtherEvents_D,    nOtherEventsErr_D    = self.get("nEventsD", disc1Key, disc2Key, "BG_OTHER_QCDCR")
            
                nTTEvents_A,    nTTEventsErr_A    = self.get("nEventsA", disc1Key, disc2Key, "TT_QCDCR")
                nTTEvents_B,    nTTEventsErr_B    = self.get("nEventsB", disc1Key, disc2Key, "TT_QCDCR")
                nTTEvents_C,    nTTEventsErr_C    = self.get("nEventsC", disc1Key, disc2Key, "TT_QCDCR")
                nTTEvents_D,    nTTEventsErr_D    = self.get("nEventsD", disc1Key, disc2Key, "TT_QCDCR")
           
                nQCDinData_A = nDataEvents_A - nTTXEvents_A - nOtherEvents_A - nTTEvents_A
                nQCDinData_B = nDataEvents_B - nTTXEvents_B - nOtherEvents_B - nTTEvents_B
                nQCDinData_C = nDataEvents_C - nTTXEvents_C - nOtherEvents_C - nTTEvents_C
                nQCDinData_D = nDataEvents_D - nTTXEvents_D - nOtherEvents_D - nTTEvents_D

                nQCDinDataErr_A = math.sqrt(nDataEventsErr_A ** 2 + nTTXEventsErr_A ** 2 + nOtherEventsErr_A ** 2 + nTTEventsErr_A ** 2) * TFs[0]
                nQCDinDataErr_B = math.sqrt(nDataEventsErr_B ** 2 + nTTXEventsErr_B ** 2 + nOtherEventsErr_B ** 2 + nTTEventsErr_B ** 2) * TFs[1]
                nQCDinDataErr_C = math.sqrt(nDataEventsErr_C ** 2 + nTTXEventsErr_C ** 2 + nOtherEventsErr_C ** 2 + nTTEventsErr_C ** 2) * TFs[2]
                nQCDinDataErr_D = math.sqrt(nDataEventsErr_D ** 2 + nTTXEventsErr_D ** 2 + nOtherEventsErr_D ** 2 + nTTEventsErr_D ** 2) * TFs[3]

                nQCDFinal_A = nQCDinData_A * TFs[0]
                nQCDFinal_B = nQCDinData_B * TFs[1]
                nQCDFinal_C = nQCDinData_C * TFs[2]
                nQCDFinal_D = nQCDinData_D * TFs[3]
     
                self.add("nEventsA", disc1Key, disc2Key, (nQCDFinal_A, nQCDinDataErr_A), "QCDFinal")
                self.add("nEventsB", disc1Key, disc2Key, (nQCDFinal_B, nQCDinDataErr_B), "QCDFinal")
                self.add("nEventsC", disc1Key, disc2Key, (nQCDFinal_C, nQCDinDataErr_C), "QCDFinal")
                self.add("nEventsD", disc1Key, disc2Key, (nQCDFinal_D, nQCDinDataErr_D), "QCDFinal")

                self.add("QCDTF_A", disc1Key, disc2Key, (TFs[0]), "QCD")
                self.add("QCDTF_B", disc1Key, disc2Key, (TFs[1]), "QCD")
                self.add("QCDTF_C", disc1Key, disc2Key, (TFs[2]), "QCD")
                self.add("QCDTF_D", disc1Key, disc2Key, (TFs[3]), "QCD")

    def make_NonTT_Regions(self):

        for disc1Key, disc2Key in self.get("edges",None,None,"TT"):

            nTTXEvents_A,    nTTXEventsErr_A    = self.get("nEventsA", disc1Key, disc2Key, "TTX")
            nTTXEvents_B,    nTTXEventsErr_B    = self.get("nEventsB", disc1Key, disc2Key, "TTX")
            nTTXEvents_C,    nTTXEventsErr_C    = self.get("nEventsC", disc1Key, disc2Key, "TTX")
            nTTXEvents_D,    nTTXEventsErr_D    = self.get("nEventsD", disc1Key, disc2Key, "TTX")
        
            nOtherEvents_A,    nOtherEventsErr_A    = self.get("nEventsA", disc1Key, disc2Key, "BG_OTHER")
            nOtherEvents_B,    nOtherEventsErr_B    = self.get("nEventsB", disc1Key, disc2Key, "BG_OTHER")
            nOtherEvents_C,    nOtherEventsErr_C    = self.get("nEventsC", disc1Key, disc2Key, "BG_OTHER")
            nOtherEvents_D,    nOtherEventsErr_D    = self.get("nEventsD", disc1Key, disc2Key, "BG_OTHER")
        
            nQCDEvents_A,    nQCDEventsErr_A    = self.get("nEventsA", disc1Key, disc2Key, "QCD")
            nQCDEvents_B,    nQCDEventsErr_B    = self.get("nEventsB", disc1Key, disc2Key, "QCD")
            nQCDEvents_C,    nQCDEventsErr_C    = self.get("nEventsC", disc1Key, disc2Key, "QCD")
            nQCDEvents_D,    nQCDEventsErr_D    = self.get("nEventsD", disc1Key, disc2Key, "QCD")

            nNonTTEvents_A = nTTXEvents_A + nOtherEvents_A + nQCDEvents_A
            nNonTTEvents_B = nTTXEvents_B + nOtherEvents_B + nQCDEvents_B
            nNonTTEvents_C = nTTXEvents_C + nOtherEvents_C + nQCDEvents_C
            nNonTTEvents_D = nTTXEvents_D + nOtherEvents_D + nQCDEvents_D

            nNonTTEventsErr_A = math.sqrt(nTTXEventsErr_A ** 2 + nOtherEventsErr_A ** 2 + nQCDEventsErr_A ** 2)
            nNonTTEventsErr_B = math.sqrt(nTTXEventsErr_B ** 2 + nOtherEventsErr_B ** 2 + nQCDEventsErr_B ** 2)
            nNonTTEventsErr_C = math.sqrt(nTTXEventsErr_C ** 2 + nOtherEventsErr_C ** 2 + nQCDEventsErr_C ** 2)
            nNonTTEventsErr_D = math.sqrt(nTTXEventsErr_D ** 2 + nOtherEventsErr_D ** 2 + nQCDEventsErr_D ** 2)

            self.add("nEventsA", disc1Key, disc2Key, (nNonTTEvents_A, nNonTTEventsErr_A), "NonTT")
            self.add("nEventsB", disc1Key, disc2Key, (nNonTTEvents_B, nNonTTEventsErr_B), "NonTT")
            self.add("nEventsC", disc1Key, disc2Key, (nNonTTEvents_C, nNonTTEventsErr_C), "NonTT")
            self.add("nEventsD", disc1Key, disc2Key, (nNonTTEvents_D, nNonTTEventsErr_D), "NonTT")
       

# -----------------------------------------
# add the all edges to DoubleDisCo Cfg file
# -----------------------------------------
class ConfigWriter():

    def __init__(self, model, channel, regions):
        self.model     = model.partition("_")[0]
        self.channel   = channel
        self.regions   = regions

    def addEdges_toDoubleDiscoCfg(self, isNonIso=False, bestEdges=None):
    
        nonIsoStr = ""
        if isNonIso:
            nonIsoStr = "_NonIsoMuon"
    
        configName = "Keras_Tensorflow%s_DoubleDisCo_Reg_%s_%s_Run2.cfg"%(nonIsoStr, self.channel, self.model)
        realPath   = os.path.realpath("../" + configName).rpartition("/")[0]
    
        f = open("%s/DoubleDisCo_Reg%s.cfg"%(realPath,     nonIsoStr), "r")
        g = open("%s/DoubleDisCo_Reg%s_Opt.cfg"%(realPath, nonIsoStr), "w") 
    
        # get the all information from DoubleDisCo cfg file
        lines = f.readlines()
        f.close()
    
        for iLine in range(0, len(lines)):
    
            line       = lines[iLine]
    
            if "binEdges_ABCD" not in line:
                g.write(line)
            else:
                futureLine = lines[iLine+1]
    
                edgeCand1 = line.rpartition(" = ")[-1].rstrip()
                edgeCand2 = futureLine.rpartition(" = ")[-1].rstrip()
    
                if (edgeCand1 != "1.00" and edgeCand1 != "0.00") and \
                   (edgeCand2 != "1.00" and edgeCand2 != "0.00"):
                    g.write(line.replace(edgeCand1, bestEdges[0]))
                    g.write(futureLine.replace(edgeCand2, bestEdges[1]))
                elif (edgeCand1 != "1.00" and edgeCand1 != "0.00") and \
                     (edgeCand2 == "1.00" or edgeCand2 == "0.00"):
                    continue
                else:
                    g.write(line)
        
        g.close()
