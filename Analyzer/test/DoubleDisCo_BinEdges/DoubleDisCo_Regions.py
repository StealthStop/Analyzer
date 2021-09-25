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

    def __init__(self, histBkg=None, histSig=None, fixedDisc1Edge=None, fixedDisc2Edge=None, leftBoundary=None, rightBoundary=None, topBoundary=None, bottomBoundary=None, metric=None, **kwargs):

        self.quantities = {}
        self.finalEdges = ()

        self.histBkg        = histBkg
        self.histSig        = histSig
        self.fixedDisc1Edge = fixedDisc1Edge
        self.fixedDisc2Edge = fixedDisc2Edge
        self.leftBoundary   = leftBoundary
        self.rightBoundary  = rightBoundary
        self.topBoundary    = topBoundary
        self.bottomBoundary = bottomBoundary
        self.metric         = metric

        self.extraArgs = kwargs

        # First calculate events counts for all possible choices of bin edges
        self.count_Events_inBinEdges()

        # Then determine the "final" choice of bin edges for the region
        # based on the optimization_metric function
        self.get_BinEdges()

    # ------------------------
    # Significance calculation
    # ------------------------
    def cal_Significance(self, nSigEvents, nBkgEvents, sys=0.3):
        if (nBkgEvents == 0.0):
            return 0.0
    
        significance = nSigEvents / ( nBkgEvents + (sys * nBkgEvents)**2.0 )**0.5
        return significance
    
    # -------------------------
    # Closure error calculation
    # -------------------------
    def cal_ClosureError(self, nEvents_A, nEvents_B, nEvents_C, nEvents_D, nEventsErr_A, nEventsErr_B, nEventsErr_C, nEventsErr_D):
    
        if nEvents_A == 0.0 or nEvents_D == 0.0:
            return -999.0, -999.0
        
        closureError = abs(1.0 - ( (nEvents_B * nEvents_C) / (nEvents_A * nEvents_D) ) )
        
        closureErrUnc = ( ( ( nEvents_C * nEventsErr_B ) / ( nEvents_A * nEvents_D) )**2.0 
                        + ( ( nEvents_B * nEventsErr_C ) / ( nEvents_A * nEvents_D) )**2.0 
                        + ( ( nEvents_B * nEvents_C * nEventsErr_A ) / ( nEvents_A**2.0 * nEvents_D ) )**2.0 
                        + ( ( nEvents_B * nEvents_C * nEventsErr_D ) / ( nEvents_A * nEvents_D**2.0 ) )**2.0 )**0.5
    
        return closureError, closureErrUnc

    # ----------------
    # Pull calculation
    # ----------------
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
            xBinKey     = "%.3f"%(xLowBinEdge)
    
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

                self.add("nSigEventsA", xBinKey, yBinKey, (nSigEvents_A, nSigEvents_A**0.5))
                self.add("nSigEventsB", xBinKey, yBinKey, (nSigEvents_B, nSigEvents_B**0.5))
                self.add("nSigEventsC", xBinKey, yBinKey, (nSigEvents_C, nSigEvents_C**0.5))
                self.add("nSigEventsD", xBinKey, yBinKey, (nSigEvents_D, nSigEvents_D**0.5))

                self.add("nBkgEventsA", xBinKey, yBinKey, (nBkgEvents_A, nBkgEvents_A**0.5))
                self.add("nBkgEventsB", xBinKey, yBinKey, (nBkgEvents_B, nBkgEvents_B**0.5))
                self.add("nBkgEventsC", xBinKey, yBinKey, (nBkgEvents_C, nBkgEvents_C**0.5))
                self.add("nBkgEventsD", xBinKey, yBinKey, (nBkgEvents_D, nBkgEvents_D**0.5))

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
    def get_BinEdges(self, bkgNormUnc = 0.3, minBkgEvents = 1, minSigEvents = 5):
      
        optMetric = 999.0

        # loop over the disc1 and disc2 to get any possible combination of them
        for disc1Key, disc2Key in self.get("edges"):

            # number of signal and background events in aech A, B, C, D region
            nSigEvents_A, nSigEventsErr_A = self.get("nSigEventsA", disc1Key, disc2Key) 
            nSigEvents_B, nSigEventsErr_B = self.get("nSigEventsB", disc1Key, disc2Key) 
            nSigEvents_C, nSigEventsErr_C = self.get("nSigEventsC", disc1Key, disc2Key) 
            nSigEvents_D, nSigEventsErr_D = self.get("nSigEventsD", disc1Key, disc2Key) 

            nBkgEvents_A, nBkgEventsErr_A = self.get("nBkgEventsA", disc1Key, disc2Key)
            nBkgEvents_B, nBkgEventsErr_B = self.get("nBkgEventsB", disc1Key, disc2Key)
            nBkgEvents_C, nBkgEventsErr_C = self.get("nBkgEventsC", disc1Key, disc2Key)
            nBkgEvents_D, nBkgEventsErr_D = self.get("nBkgEventsD", disc1Key, disc2Key)

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
    
            # significance, closure error and pull for optimization of bin edges
            significance = 0.0; significanceUnc = 0.0; closureErr = -999.0; closureErrUnc = -999.0; pull = -999.0; pullUnc = -999.0; tempOptMetric = 999.0
    
            closureErr, closureErrUnc = self.cal_ClosureError(nBkgEvents_A, nBkgEvents_B, nBkgEvents_C, nBkgEvents_D, nBkgEventsErr_A, nBkgEventsErr_B, nBkgEventsErr_C, nBkgEventsErr_D)
            pull, pullUnc             = self.cal_Pull(nBkgEvents_A, nBkgEvents_B, nBkgEvents_C, nBkgEvents_D, nBkgEventsErr_A, nBkgEventsErr_B, nBkgEventsErr_C, nBkgEventsErr_D)

            significance += self.cal_Significance(nSigEvents_A, nBkgEvents_A)**2.0
            significance = significance**0.5
   
            # get the significanceUncs to plot variable vs disc as 1D 
            if nBkgEvents_A > 0.0:
                significanceUnc = (   ( nSigEventsErr_A / (nBkgEvents_A + (bkgNormUnc * nBkgEvents_A)**2.0 + (closureErr * nBkgEvents_A)**2.0)**0.5 )**2.0
                                  + ( ( nSigEvents_A * nBkgEventsErr_A * (2.0 * nBkgEvents_A * closureErr**2.0 + 2.0 * bkgNormUnc**2.0 * nBkgEvents_A + 1) ) / ( nBkgEvents_A + (bkgNormUnc * nBkgEvents_A)**2.0 + (closureErr * nBkgEvents_A)**2.0 )**1.5 )**2.0
                                  + ( ( nBkgEvents_A**2.0 * closureErr * nSigEvents_A * closureErrUnc) / ( nBkgEvents_A * ( nBkgEvents_A * (closureErr**2.0 + bkgNormUnc**2.0) + 1) )**1.5 )**2.0 )**0.5
 
            # Store significance, closure error and pull
            self.add("significance", disc1Key, disc2Key, (significance, significanceUnc))
            self.add("closureError", disc1Key, disc2Key, (closureErr,   closureErrUnc  ))
            self.add("pull",         disc1Key, disc2Key, (pull,         pullUnc        ))

            # Region by region fraction and fraction error
            self.add("sigFractionA", disc1Key, disc2Key, (sigFracsA, sigFracsErrA))
            self.add("sigFractionB", disc1Key, disc2Key, (sigFracsB, sigFracsErrB))
            self.add("sigFractionC", disc1Key, disc2Key, (sigFracsC, sigFracsErrC))
            self.add("sigFractionD", disc1Key, disc2Key, (sigFracsD, sigFracsErrD))
   
            # If in region with minimum number of bkg events, compute the metric value
            tempOptMetric = self.optimization_metric(significance=significance, closureError=closureErr, minBkgEvents=minBkgEvents,  
                                                     nBkgA=nBkgEvents_A, nBkgB=nBkgEvents_B, nBkgC=nBkgEvents_C, nBkgD=nBkgEvents_D,
                                                     disc1=float(disc1Key), disc2=float(disc2Key))
            
            # use this statement for cdGH regions if  fixed disc1 edge (vertivcal edge)
            if self.fixedDisc1Edge != None and abs(float(self.fixedDisc1Edge) - float(disc1Key)) > 0.01: continue

            # use this statement for bdEF regions if fixed disc2 edge (horizontal edge)
            if self.fixedDisc2Edge != None and abs(float(self.fixedDisc2Edge) - float(disc2Key)) > 0.01: continue

            # Determine based on the metric value if the current
            # choice of bin edges is better and if so, save them
            if tempOptMetric < optMetric:
                self.finalEdges = (disc1Key, disc2Key)
                optMetric = tempOptMetric

    # ----------------------------------
    # store quantities to make any plots
    # with any combination of bin edges
    # i.e. significance, closure, etc.
    # ----------------------------------
    def add(self, name, disc1, disc2, val):

        if disc1 not in self.quantities:
            self.quantities[disc1] = {}

        if disc2 not in self.quantities[disc1]:
            self.quantities[disc1][disc2] = {}

        self.quantities[disc1][disc2][name] = val

    # -------------------------------------------------
    # get specific variables to make any plots
    # with any combination of bin edges
    # i.e. nSigEventsA, nBkgEventsA, sigFractionA, etc.
    # -------------------------------------------------
    def get(self, name, disc1=None, disc2=None):

        if name == "edges":
            if disc1 != None and disc2 != None:
                return self.finalEdges
            else:
                payload = []
                for disc1, disc2s in self.quantities.items():
                    for disc2, _ in disc2s.items():
                        payload.append((disc1,disc2))
                return payload

        if disc1 == None and disc2 == None:
            payload = []
            for disc1, disc2s in self.quantities.items():
                for disc2, q in disc2s.items():
                    payload.append(q[name])
                    
            return np.array(payload)

        if disc1 != None and disc1 not in self.quantities:
            raise LookupError("Cannot find the key \"%s\" in the dictionary !"%(disc1))

        if disc2 != None and disc2 not in self.quantities[disc1]:
            raise LookupError("Cannot find the key \"%s\" in the dictionary !"%(disc2))

        return self.quantities[disc1][disc2][name]
        
    # -------------------------------------------------------------
    # store quatities and get spesific variables to make the tables 
    # with the final choice of bin edges
    # -------------------------------------------------------------
    def getFinal(self, name):
        return self.get(name, self.finalEdges[0], self.finalEdges[1])


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
                optimizationMetric  = (kwargs["closureError"])**2 + (1.0 / kwargs["significance"])**2
                   
            # New optimization metric
            else: 
                optimizationMetric  = (5 * kwargs["closureError"])**2 + (1.0 / kwargs["significance"])**2
        
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

                    
