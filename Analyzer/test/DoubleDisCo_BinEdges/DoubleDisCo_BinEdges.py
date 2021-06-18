import ROOT
import os
import sys
import math
import argparse
import numpy as np
import array as arr

import matplotlib as mpl
mpl.use('Agg')

import matplotlib.lines as ml
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
#import mplhep as hep
#plt.style.use([hep.style.ROOT,hep.style.CMS]) # For now ROOT defaults to CMS
#plt.style.use({'legend.frameon':False,'legend.fontsize':16,'legend.edgecolor':'black'})

from matplotlib.colors import LogNorm



# ------------------------
# Significance calculation
# ------------------------
def cal_Significance(nSigEvents, nBkgEvents, sys=0.3):
    if (nBkgEvents == 0.0):
        return 0

    significance = ( nSigEvents / ( nBkgEvents + (sys * nBkgEvents)**2.0 )**0.5 )**2.0
    return significance

# -------------------------
# Closure error calculation
# -------------------------
def cal_ClosureError(nBkgEvents_A, nBkgEvents_B, nBkgEvents_C, nBkgEvents_D):
    closureError = abs(1.0 - ( (nBkgEvents_B * nBkgEvents_C) / (nBkgEvents_A * nBkgEvents_D) ) )
    return closureError

# --------------------------------------------
# Optimization metric of bin edges calculation
# --------------------------------------------
def cal_OptMetricOfBinEdges(inverseSignificance, closureError):
    #inverseSignificance = (1.0 / significance)
    optimizationMetric  = (closureError)**2 + (inverseSignificance)**2
    return optimizationMetric

# -------------------------------------------------------
# get signal and background histograms' counts
#   -- count both signal and background events separately 
#   -- in each A, B, C, D regions
#   -- put them to the dictionaries
# -------------------------------------------------------
def count_Events_inBinEdges(histBkg, histSig):

    lastXBin = histBkg.GetNbinsX()
    lastYBin = histBkg.GetNbinsY()
    nXBins   = range(2, lastXBin+1)
    nYBins   = range(2, lastYBin+1)

    nTotSigCount_ABCD = { 
        "nSigEvents_A" : {},    "nSigEvents_B" : {},    "nSigEvents_C" : {},    "nSigEvents_D" : {}, 
        "nSigEventsErr_A" : {}, "nSigEventsErr_B" : {}, "nSigEventsErr_C" : {}, "nSigEventsErr_D" : {}
    }
    
    nTotBkgCount_ABCD = { 
        "nBkgEvents_A" : {},    "nBkgEvents_B" : {},    "nBkgEvents_C" : {},    "nBkgEvents_D" : {}, 
        "nBkgEventsErr_A" : {}, "nBkgEventsErr_B" : {}, "nBkgEventsErr_C" : {}, "nBkgEventsErr_D" : {} 
    }

    # loop over the x bins
    for xBin in nXBins:

        xLowBinEdge = histBkg.GetXaxis().GetBinLowEdge(xBin)
        xBinKey     = "%.2f"%xLowBinEdge

        if xBinKey not in nTotSigCount_ABCD["nSigEventsErr_A"]:
            nTotSigCount_ABCD["nSigEvents_A"][xBinKey] = {}
            nTotSigCount_ABCD["nSigEvents_B"][xBinKey] = {}
            nTotSigCount_ABCD["nSigEvents_C"][xBinKey] = {}
            nTotSigCount_ABCD["nSigEvents_D"][xBinKey] = {}
            nTotSigCount_ABCD["nSigEventsErr_A"][xBinKey] = {}
            nTotSigCount_ABCD["nSigEventsErr_B"][xBinKey] = {}
            nTotSigCount_ABCD["nSigEventsErr_C"][xBinKey] = {}
            nTotSigCount_ABCD["nSigEventsErr_D"][xBinKey] = {}


        if xBinKey not in nTotBkgCount_ABCD["nBkgEventsErr_A"]:
            nTotBkgCount_ABCD["nBkgEvents_A"][xBinKey] = {}
            nTotBkgCount_ABCD["nBkgEvents_B"][xBinKey] = {}
            nTotBkgCount_ABCD["nBkgEvents_C"][xBinKey] = {}
            nTotBkgCount_ABCD["nBkgEvents_D"][xBinKey] = {}
            nTotBkgCount_ABCD["nBkgEventsErr_A"][xBinKey] = {}
            nTotBkgCount_ABCD["nBkgEventsErr_B"][xBinKey] = {}
            nTotBkgCount_ABCD["nBkgEventsErr_C"][xBinKey] = {}
            nTotBkgCount_ABCD["nBkgEventsErr_D"][xBinKey] = {}

        # loop over the y bins
        for yBin in nYBins:

            yLowBinEdge = histBkg.GetYaxis().GetBinLowEdge(yBin)
            yBinKey     = "%.2f"%yLowBinEdge

            if yBinKey not in nTotSigCount_ABCD["nSigEventsErr_A"]:
                nTotSigCount_ABCD["nSigEvents_A"][xBinKey][yBinKey] = 0.0
                nTotSigCount_ABCD["nSigEvents_B"][xBinKey][yBinKey] = 0.0
                nTotSigCount_ABCD["nSigEvents_C"][xBinKey][yBinKey] = 0.0
                nTotSigCount_ABCD["nSigEvents_D"][xBinKey][yBinKey] = 0.0
                nTotSigCount_ABCD["nSigEventsErr_A"][xBinKey][yBinKey] = 0.0
                nTotSigCount_ABCD["nSigEventsErr_B"][xBinKey][yBinKey] = 0.0
                nTotSigCount_ABCD["nSigEventsErr_C"][xBinKey][yBinKey] = 0.0
                nTotSigCount_ABCD["nSigEventsErr_D"][xBinKey][yBinKey] = 0.0

            if yBinKey not in nTotBkgCount_ABCD["nBkgEventsErr_A"]:
                nTotBkgCount_ABCD["nBkgEvents_A"][xBinKey][yBinKey] = 0.0
                nTotBkgCount_ABCD["nBkgEvents_B"][xBinKey][yBinKey] = 0.0
                nTotBkgCount_ABCD["nBkgEvents_C"][xBinKey][yBinKey] = 0.0
                nTotBkgCount_ABCD["nBkgEvents_D"][xBinKey][yBinKey] = 0.0
                nTotBkgCount_ABCD["nBkgEventsErr_A"][xBinKey][yBinKey] = 0.0
                nTotBkgCount_ABCD["nBkgEventsErr_B"][xBinKey][yBinKey] = 0.0
                nTotBkgCount_ABCD["nBkgEventsErr_C"][xBinKey][yBinKey] = 0.0
                nTotBkgCount_ABCD["nBkgEventsErr_D"][xBinKey][yBinKey] = 0.0

            # count signal and background events and errors in bin edges
            nSigEventsErr_A = ROOT.Double(0.0); nSigEventsErr_B = ROOT.Double(0.0); nSigEventsErr_C = ROOT.Double(0.0); nSigEventsErr_D = ROOT.Double(0.0)
            nBkgEventsErr_A = ROOT.Double(0.0); nBkgEventsErr_B = ROOT.Double(0.0); nBkgEventsErr_C = ROOT.Double(0.0); nBkgEventsErr_D = ROOT.Double(0.0)
            nSigEvents_A = ( histSig.IntegralAndError(xBin, lastXBin, yBin, lastYBin, nSigEventsErr_A) )
            nSigEvents_B = ( histSig.IntegralAndError(1, xBin-1, yBin, lastYBin, nSigEventsErr_B) )
            nSigEvents_C = ( histSig.IntegralAndError(xBin, lastXBin, 1, yBin-1, nSigEventsErr_C) )
            nSigEvents_D = ( histSig.IntegralAndError(1, xBin-1, 1, yBin-1, nSigEventsErr_D) )
            nBkgEvents_A = ( histBkg.IntegralAndError(xBin, lastXBin, yBin, lastYBin, nBkgEventsErr_A) )
            nBkgEvents_B = ( histBkg.IntegralAndError(1, xBin-1, yBin, lastYBin, nBkgEventsErr_B) )
            nBkgEvents_C = ( histBkg.IntegralAndError(xBin, lastXBin, 1, yBin-1, nBkgEventsErr_C) )
            nBkgEvents_D = ( histBkg.IntegralAndError(1, xBin-1, 1, yBin-1, nBkgEventsErr_D) )
            
            nTotSigCount_ABCD["nSigEvents_A"][xBinKey][yBinKey] = nSigEvents_A
            nTotSigCount_ABCD["nSigEvents_B"][xBinKey][yBinKey] = nSigEvents_B
            nTotSigCount_ABCD["nSigEvents_C"][xBinKey][yBinKey] = nSigEvents_C
            nTotSigCount_ABCD["nSigEvents_D"][xBinKey][yBinKey] = nSigEvents_D
            nTotSigCount_ABCD["nSigEventsErr_A"][xBinKey][yBinKey] = nSigEventsErr_A
            nTotSigCount_ABCD["nSigEventsErr_B"][xBinKey][yBinKey] = nSigEventsErr_B
            nTotSigCount_ABCD["nSigEventsErr_C"][xBinKey][yBinKey] = nSigEventsErr_C
            nTotSigCount_ABCD["nSigEventsErr_D"][xBinKey][yBinKey] = nSigEventsErr_D

            nTotBkgCount_ABCD["nBkgEvents_A"][xBinKey][yBinKey] = nBkgEvents_A
            nTotBkgCount_ABCD["nBkgEvents_B"][xBinKey][yBinKey] = nBkgEvents_B
            nTotBkgCount_ABCD["nBkgEvents_C"][xBinKey][yBinKey] = nBkgEvents_C
            nTotBkgCount_ABCD["nBkgEvents_D"][xBinKey][yBinKey] = nBkgEvents_D
            nTotBkgCount_ABCD["nBkgEventsErr_A"][xBinKey][yBinKey] = nBkgEventsErr_A
            nTotBkgCount_ABCD["nBkgEventsErr_B"][xBinKey][yBinKey] = nBkgEventsErr_B
            nTotBkgCount_ABCD["nBkgEventsErr_C"][xBinKey][yBinKey] = nBkgEventsErr_C
            nTotBkgCount_ABCD["nBkgEventsErr_D"][xBinKey][yBinKey] = nBkgEventsErr_D

    return nTotSigCount_ABCD, nTotBkgCount_ABCD

# -------------------------------------------------------------------------
# Region by region signal fraction calculation
# look at (SIG + BKG) events in each region, how much signal is in each bin
#   -- in region A : Nsig / (Nbkg + Nsig) 
#   -- in region B : Nsig / (Nbkg + Nsig) 
#   -- in region C : Nsig / (Nbkg + Nsig) 
#   -- in region D : Nsig / (Nbkg + Nsig) 
#
# Total (A+B+C+D) signal fraction calculation
# look at (SIG) events in each region, how many events are signal
#   Ntotal = NA + NB + NC + ND
#   -- SigFragA = NA / Ntotal
#   -- SigFragB = NB / Ntotal
#   -- SigFragC = NC / Ntotal
#   -- SigFragD = ND / Ntotal 
# -------------------------------------------------------------------------
def calc_Sig_SigBkg_Fractions(nTotSigCount_ABCD, nTotBkgCount_ABCD, minBkgFrac = 0.01, minSigFrac = 0.1):
  
    significance = 0.0; finalDisc1Key = -1.0; finalDisc2Key = -1.0; closureErr   = 0.0; optMetric  = 999.0
    inverseSignificance = []; closureErrsList = []; disc1KeyOut  = []; disc2KeyOut  = []
    sigFracsA           = []; sigFracsB       = []; sigFracsC    = []; sigFracsD    = []
    sigTotFracsA        = []; sigTotFracsB    = []; sigTotFracsC = []; sigTotFracsD = []
    bkgTotFracsA        = []; bkgTotFracsB    = []; bkgTotFracsC = []; bkgTotFracsD = []

    # loop over the disc1 and disc2 to get any possible combination of them
    for disc1Key, disc2s in nTotBkgCount_ABCD["nBkgEvents_A"].items():
        
        for disc2Key, nEvents in disc2s.items():

            # number of signal and background events in aech A, B, C, D region
            nSigEvents_A = nTotSigCount_ABCD["nSigEvents_A"][disc1Key][disc1Key] 
            nSigEvents_B = nTotSigCount_ABCD["nSigEvents_B"][disc1Key][disc1Key]
            nSigEvents_C = nTotSigCount_ABCD["nSigEvents_C"][disc1Key][disc1Key]
            nSigEvents_D = nTotSigCount_ABCD["nSigEvents_D"][disc1Key][disc1Key]
            
            nBkgEvents_A = nTotBkgCount_ABCD["nBkgEvents_A"][disc1Key][disc1Key]
            nBkgEvents_B = nTotBkgCount_ABCD["nBkgEvents_B"][disc1Key][disc1Key]
            nBkgEvents_C = nTotBkgCount_ABCD["nBkgEvents_C"][disc1Key][disc1Key]
            nBkgEvents_D = nTotBkgCount_ABCD["nBkgEvents_D"][disc1Key][disc1Key]

            # Region by region signal fraction calculation
            # get some plots based on region by region signal fraction
            nTot_SigBkg_A = nSigEvents_A + nBkgEvents_A
            nTot_SigBkg_B = nSigEvents_B + nBkgEvents_B
            nTot_SigBkg_C = nSigEvents_C + nBkgEvents_C
            nTot_SigBkg_D = nSigEvents_D + nBkgEvents_D
            tempSigFracA = -1.0; tempSigFracB = -1.0; tempSigFracC = -1.0; tempSigFracD = -1.0

            if nTot_SigBkg_A > 0.0: tempSigFracsA = nSigEvents_A / nTot_SigBkg_A
            if nTot_SigBkg_B > 0.0: tempSigFracsB = nSigEvents_B / nTot_SigBkg_B
            if nTot_SigBkg_C > 0.0: tempSigFracsC = nSigEvents_C / nTot_SigBkg_C
            if nTot_SigBkg_D > 0.0: tempSigFracsD = nSigEvents_D / nTot_SigBkg_D

            sigFracsA.append(float(tempSigFracA))
            sigFracsB.append(float(tempSigFracB))
            sigFracsC.append(float(tempSigFracC))
            sigFracsD.append(float(tempSigFracD))

            # Total signal (and background) fractions in aech A, B, C, D region
            # get the latest bin edges based on total signal fraction
            nTot_Sig_ABCD = nSigEvents_A + nSigEvents_B + nSigEvents_C + nSigEvents_D
            nTot_Bkg_ABCD = nBkgEvents_A + nBkgEvents_B + nBkgEvents_C + nBkgEvents_D
            tempSigTotFracsA = -1.0; tempSigTotFracsB = -1.0; tempSigTotFracsC = -1.0; tempSigTotFracsD = -1.0
            tempBkgTotFracsA = -1.0; tempBkgTotFracsB = -1.0; tempBkgTotFracsC = -1.0; tempBkgTotFracsD = -1.0
            
            tempSigTotFracsA = nSigEvents_A / nTot_Sig_ABCD 
            tempSigTotFracsB = nSigEvents_B / nTot_Sig_ABCD
            tempSigTotFracsC = nSigEvents_C / nTot_Sig_ABCD 
            tempSigTotFracsD = nSigEvents_D / nTot_Sig_ABCD
            
            tempBkgTotFracsA = nBkgEvents_A / nTot_Bkg_ABCD   
            tempBkgTotFracsB = nBkgEvents_B / nTot_Bkg_ABCD  
            tempBkgTotFracsC = nBkgEvents_C / nTot_Bkg_ABCD 
            tempBkgTotFracsD = nBkgEvents_D / nTot_Bkg_ABCD

            # significance and closure error for optimization of bin edges
            tempSignificance = 0.0; tempClosureErr = -999.0; tempOptMetric = 999.0

            if nBkgEvents_A > 0.0:
                tempSignificance += ( nSigEvents_A / ( nBkgEvents_A + (0.3 * nBkgEvents_A)**2.0 )**0.5 )**2.0 
            if nBkgEvents_B > 0.0:
                tempSignificance += ( nSigEvents_B / ( nBkgEvents_B + (0.3 * nBkgEvents_B)**2.0 )**0.5 )**2.0
            if nBkgEvents_C > 0.0:
                tempSignificance += ( nSigEvents_C / ( nBkgEvents_C + (0.3 * nBkgEvents_C)**2.0 )**0.5 )**2.0
            if nBkgEvents_D > 0.0:
                tempSignificance += ( nSigEvents_D / ( nBkgEvents_D + (0.3 * nBkgEvents_D)**2.0 )**0.5 )**2.0
            
            if nBkgEvents_A > 0.0 and nBkgEvents_D > 0.0: 
                tempClosureErr = cal_ClosureError(nBkgEvents_A, nBkgEvents_B, nBkgEvents_C, nBkgEvents_D)

            tempSignificance = tempSignificance**0.5

            if tempSignificance > 0.0 and tempClosureErr > 0.0:
                inverseSignificance.append(1.0 / tempSignificance) 
                closureErrsList.append(abs(tempClosureErr))
                disc1KeyOut.append(float(disc1Key))
                disc2KeyOut.append(float(disc2Key)) 
        
                sigTotFracsA.append(float(tempSigTotFracsA))
                sigTotFracsB.append(float(tempSigTotFracsB))
                sigTotFracsC.append(float(tempSigTotFracsC))
                sigTotFracsD.append(float(tempSigTotFracsD))

                bkgTotFracsA.append(float(tempBkgTotFracsA))
                bkgTotFracsB.append(float(tempBkgTotFracsB))
                bkgTotFracsC.append(float(tempBkgTotFracsC))
                bkgTotFracsD.append(float(tempBkgTotFracsD))

            if (tempBkgTotFracsA > minBkgFrac) and (tempBkgTotFracsB > minBkgFrac) and (tempBkgTotFracsC > minBkgFrac) and (tempBkgTotFracsD > minBkgFrac):                
                tempOptMetric = cal_OptMetricOfBinEdges(tempSignificance, tempClosureErr)

            if tempOptMetric < optMetric:
                finalDisc1Key = disc1Key 
                finalDisc2Key = disc2Key
                significance  = tempSignificance
                closureErr    = tempClosureErr
                optMetric     = tempOptMetric

    # print out the disc1 and disc2 edges
    #print "tempSignificance", tempSignificance
    #print "disc1 (x bin) low bin edges: ", finalDisc1Key
    #print "disc2 (y bin) low bin edges: ", finalDisc2Key

    return finalDisc1Key, finalDisc2Key, significance, closureErr, inverseSignificance, closureErrsList, disc1KeyOut, disc2KeyOut, sigFracsA, sigFracsB, sigFracsC, sigFracsD, sigTotFracsA, sigTotFracsB, sigTotFracsC, sigTotFracsD, bkgTotFracsA, bkgTotFracsB, bkgTotFracsC, bkgTotFracsD 


# --------------------------------
# calculate everything for closure
# --------------------------------
def cal_simpleClosureABCD(nBkgEvents_A, nBkgEvents_B, nBkgEvents_C, nBkgEvents_D, nBkgEventsErr_A, nBkgEventsErr_B, nBkgEventsErr_C, nBkgEventsErr_D):

    # Define A: > c1, > c2        C    |    A    
    # Define B: > c1, < c2   __________|__________        
    # Define C: < c1, > c2             |        
    # Define D: < c1, < c2        D    |    B   

    numerator   = nBkgEvents_B * nBkgEvents_C 
    denominator = nBkgEvents_A * nBkgEvents_D

    nPredBkgEvents_A = -1.0; nPredBkgUncEvents_A = 0.0

    if nBkgEvents_D > 0.0:
        nPredBkgEvents_A    = ( numerator / nBkgEvents_D )
        nPredBkgUncEvents_A = ( ( (nBkgEvents_C * nBkgEventsErr_B) / nBkgEvents_D )**2.0 
                              + ( (nBkgEvents_B * nBkgEventsErr_C) / nBkgEvents_B )**2.0 
                              + ( (nBkgEvents_B * nBkgEvents_C * nBkgEventsErr_D) / nBkgEvents_D**2.0 )**2.0 )**0.5

    if denominator > 0.0:
        closureErr = ( ( (nBkgEvents_B * nBkgEventsErr_C) / denominator )**2.0 
                     + ( (nBkgEvents_C * nBkgEventsErr_B) / denominator )**2.0 
                     + ( (numerator * nBkgEventsErr_A) / (denominator * nBkgEvents_A) )**2.0 
                     + ( (numerator * nBkgEventsErr_D) / (denominator * nBkgEvents_D) )**2.0 )**0.5
        closure = numerator / denominator
    
    else:
        closureErr = -999.0
        closure    = -999.0

    return closure, closureErr, nPredBkgEvents_A, nPredBkgUncEvents_A 

# ---------------------------------------------------------
# plot closure
#   -- chi2 is way to measure the closure
#   -- totalChi2 = [ (A_predicted - A_actual) / delta_A ]^2
# ---------------------------------------------------------
def plot_ClosureNjets(bkg, bkgUnc, bkgPred, bkgPredUnc, Njets, year, channel):

    binCenters = []; xErr    = []; abcdPull = []; abcdError = []
    unc        = []; predUnc = []; obs      = []; pred      = []

    totalChi2 = 0.0; ndof = 0

    print "bkg Unc: ", bkgUnc

    for i in range(0, len(Njets)):
        
        if bkgUnc[i] != 0.0:

            binCenters.append(Njets[i])
            unc.append(bkgUnc[i])
            predUnc.append(bkgPredUnc[i])
            obs.append(bkg[i])
            pred.append(bkgPred[i])
            xErr.append(0.5)
            pull              = ( bkgPred[i] - bkg[i] ) / bkgUnc[i] 
            closureError      = 1.0 - ( bkgPred[i] / bkg[i] )
            abcdPull.append(pull)
            abcdError.append(closureError)
            totalChi2         += pull**2.0
            ndof              += 1    
        
    fig = plt.figure()
    gs  = fig.add_gridspec(8, 1)
    ax  = fig.add_subplot(111)
    ax1 = fig.add_subplot(gs[0:6])
    ax2 = fig.add_subplot(gs[6:8], sharex=ax1)    
    
    fig.subplots_adjust(hspace=0)

    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)

    x1.set_yscale("log")
    ax1.text(0.05, 0.1, "$\chi^2$ / ndof = %3.2f"%(totalChi2/float(ndof)), horizontalalignment="left", verticalalignment="center", transform=ax1.transAxes, fontsize=18)
       
    ax1.errorbar(binCenters, pred, yerr=predUnc, label="Observed",  xerr=xErr, fmt='', color="red",   lw=0, elinewidth=2, marker="o", markerfacecolor="red")
    ax1.errorbar(binCenters, obs,  yerr=unc,     label="Predicted", xerr=xErr, fmt='', color="black", lw=0, elinewidth=2, marker="o", markerfacecolor="black")

    lowerNjets  = Njets[0]
    higherNjets = Njets[-1]

    ax1.set_xlim([lowerNjets-0.5,higherNjets+0.5])
        
    plt.xticks(Njets)

    ax2.errorbar(binCenters, abcdError, yerr=None, xerr=xErr, fmt='', color="blue",  lw=0, elinewidth=2, marker="o", markerfacecolor="blue")
    ax2.axhline(y=0.0, color="black", linestyle="dashed", lw=1)
    ax2.grid(color="black", which="both", axis="y")

    #hep.cms.label(data=True, paper=False, year=year, ax=ax1)
    
    ax2.set_xlabel('Number of jets')
    ax2.set_ylabel('1 - Pred./Obs.', fontsize="small")
    ax1.set_ylabel('Unweighted Event Counts')
    ax1.legend(loc='best')
    ax2.set_ylim([-1.6, 1.6])
    
    fig.savefig("plots/%s/Njets_Region_A_PredVsActual.pdf" %(channel))

    plt.close(fig)

    return totalChi2, ndof


# --------------------------------
# plot Disc1VsDisc2 - significance
# --------------------------------
#def plot_Disc1vsDisc2_Significance(inverseSignificance, closureErr, disc1Edges, disc2Edges, c1, c2, minEdge, maxEdge, edgeWidth, Njets = -1):

    




# ------------------------------
# plot Disc1VsDisc2 - comparison
# ------------------------------
#def plot_Disc1vsDisc2_Comparison(finalSignificance, finalClosureErr, inverseSignificance, closureErr, edges, disc1Edge, disc2Edge, Njets = -1):



# ------------------------------
# plot Disc1VsDisc2 - percentage
# ------------------------------
#def plot_Disc1vsDisc2_Percentage(sigFracsA, sigFracsB, sigFracsC, sigFracsD,
#                                 sigTotFracsA, sigTotFracsB, sigTotFracsC, sigTotFracsD, 
#                                 bkgTotFracsA, bkgTotFracsB, bkgTotFracsC, bkgTotFracsD, 
#                                 disc1Edges, disc2Edges, c1, c2, minEdge, maxEdge, edgeWidth, Njets = -1): 










def main():

    # -----------------------------------------------------------------------------------
    # command to run this script
    #   -- python DoubleDisCo_BinEdges.py --year 2016 --model RPV --mass 550 --channel 0l
    #   -- python DoubleDisCo_BinEdges.py --year 2016 --model RPV --mass 550 --channel 1l
    # -----------------------------------------------------------------------------------  
    usage  = "usage: %prog [options]"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("--year",    dest="year",    help="which year",            required=True)
    parser.add_argument("--path",    dest="path",    help="Input dir with histos", default="/uscms_data/d3/jhiltb/PO_Boxes/Bryan/2016_DisCo_0l_1l_Inputs/")
    parser.add_argument("--model",   dest="model",   help="signal model",          default="RPV")
    parser.add_argument("--mass",    dest="mass",    help="signal mass",           default="550")
    parser.add_argument("--channel", dest="channel", help="0l, 1l",                required=True)
    args = parser.parse_args()

    modelDecay = "2t6j"
    if ("SHH" in args.model):
        modelDecay = "2t4b"

    files = {
        "TT"                          : ROOT.TFile.Open(args.path + "/" + args.year + "_TT.root"),
        "QCD"                         : ROOT.TFile.Open(args.path + "/" + args.year + "_QCD.root"),  
        "%s%s"%(args.model,args.mass) : ROOT.TFile.Open(args.path + "/" + args.year + "_%s_%s_mStop-%s.root"%(args.model,modelDecay,args.mass)),
    }

    #os.system("rm BinEdges_%s_%s_%s.txt" %(args.model, args.mass, args.channel))

    # for 0-lepton 
    if args.channel == "0l":
        histNames = "h_DoubleDisCo_disc1_disc2_0l"
        njets = [
            "_Njets7",
            "_Njets8",
            "_Njets9",
            "_Njets10",
            "_Njets11",
            "_Njets12",
        ]

    # for 1-lepton
    else:
        histNames = "h_DoubleDisCo_disc1_disc2_1l"
        njets = [
            "_Njets7",
            "_Njets8",
            "_Njets9",
            "_Njets10",
            "_Njets11",
        ]

    # initialize all variables
    sigNjets       = { "A" : [], "B" : [], "C" : [], "D" : [] }
    sigNjetsErr    = { "A" : [], "B" : [], "C" : [], "D" : [] }
    bkgNjets       = { "A" : [], "B" : [], "C" : [], "D" : [] }
    bkgNjetsErr    = { "A" : [], "B" : [], "C" : [], "D" : [] }
    bkgNjetsPred_A = { "value" : [], "error" : []}
    Njets = None

    # loop over njets 
    for njet in njets: 

        # ------------------------
        # get the latest bin edges
        # ------------------------
        histBkg = files["TT"].Get(histNames + njet)
        histSig = files["%s%s"%(args.model,args.mass)].Get(histNames + njet)

        nTotSigCount_ABCD, nTotBkgCount_ABCD = count_Events_inBinEdges(histBkg, histSig)
        
        finalDisc1Key, finalDisc2Key, significance, closureErr, inverseSignificance, closureErrsList, disc1KeyOut, disc2KeyOut, sigFracsA, sigFracsB, sigFracsC, sigFracsD, sigTotFracsA, sigTotFracsB, sigTotFracsC, sigTotFracsD, bkgTotFracsA, bkgTotFracsB, bkgTotFracsC, bkgTotFracsD = calc_Sig_SigBkg_Fractions(nTotSigCount_ABCD, nTotBkgCount_ABCD, minBkgFrac = 0.01, minSigFrac = 0.1)

        # -----------------------------
        # calculate simple closure ABCD
        # -----------------------------
        if finalDisc1Key == -1.0 or finalDisc2Key == -1.0:
           
            sigNjets["A"].append(0.0); sigNjetsErr["A"].append(0.0) 
            sigNjets["B"].append(0.0); sigNjetsErr["B"].append(0.0)
            sigNjets["C"].append(0.0); sigNjetsErr["C"].append(0.0)
            sigNjets["D"].append(0.0); sigNjetsErr["D"].append(0.0)
 
            bkgNjets["A"].append(0.0); bkgNjetsErr["A"].append(0.0) 
            bkgNjets["B"].append(0.0); bkgNjetsErr["B"].append(0.0)
            bkgNjets["C"].append(0.0); bkgNjetsErr["C"].append(0.0)
            bkgNjets["D"].append(0.0); bkgNjetsErr["D"].append(0.0)

            bkgNjetsAPred["value"].append(0.0); bkgNjetsAPred["error"].append(0.0)
        
        else:

            closure, closureUnc, pred_A, predUnc_A = cal_simpleClosureABCD( nTotBkgCount_ABCD["nBkgEvents_A"][finalDisc1Key][finalDisc2Key], 
                                                                            nTotBkgCount_ABCD["nBkgEvents_B"][finalDisc1Key][finalDisc2Key], 
                                                                            nTotBkgCount_ABCD["nBkgEvents_C"][finalDisc1Key][finalDisc2Key], 
                                                                            nTotBkgCount_ABCD["nBkgEvents_D"][finalDisc1Key][finalDisc2Key], 
                                                                            nTotBkgCount_ABCD["nBkgEventsErr_A"][finalDisc1Key][finalDisc2Key], 
                                                                            nTotBkgCount_ABCD["nBkgEventsErr_B"][finalDisc1Key][finalDisc2Key], 
                                                                            nTotBkgCount_ABCD["nBkgEventsErr_C"][finalDisc1Key][finalDisc2Key], 
                                                                            nTotBkgCount_ABCD["nBkgEventsErr_D"][finalDisc1Key][finalDisc2Key] )

            # ---------------------
            # make the closure plot
            # ---------------------
            sigNjets["A"].append(nTotSigCount_ABCD["nSigEvents_A"][finalDisc1Key][finalDisc2Key])
            sigNjets["B"].append(nTotSigCount_ABCD["nSigEvents_B"][finalDisc1Key][finalDisc2Key])
            sigNjets["C"].append(nTotSigCount_ABCD["nSigEvents_C"][finalDisc1Key][finalDisc2Key])
            sigNjets["D"].append(nTotSigCount_ABCD["nSigEvents_D"][finalDisc1Key][finalDisc2Key])
            sigNjetsErr["A"].append(nTotSigCount_ABCD["nSigEventsErr_A"][finalDisc1Key][finalDisc2Key])
            sigNjetsErr["B"].append(nTotSigCount_ABCD["nSigEventsErr_B"][finalDisc1Key][finalDisc2Key])
            sigNjetsErr["C"].append(nTotSigCount_ABCD["nSigEventsErr_C"][finalDisc1Key][finalDisc2Key])
            sigNjetsErr["D"].append(nTotSigCount_ABCD["nSigEventsErr_D"][finalDisc1Key][finalDisc2Key])

            bkgNjets["A"].append(nTotBkgCount_ABCD["nBkgEvents_A"][finalDisc1Key][finalDisc2Key]) 
            bkgNjets["B"].append(nTotBkgCount_ABCD["nBkgEvents_B"][finalDisc1Key][finalDisc2Key]) 
            bkgNjets["C"].append(nTotBkgCount_ABCD["nBkgEvents_C"][finalDisc1Key][finalDisc2Key]) 
            bkgNjets["D"].append(nTotBkgCount_ABCD["nBkgEvents_D"][finalDisc1Key][finalDisc2Key])
            bkgNjetsErr["A"].append(nTotBkgCount_ABCD["nBkgEventsErr_A"][finalDisc1Key][finalDisc2Key]) 
            bkgNjetsErr["B"].append(nTotBkgCount_ABCD["nBkgEventsErr_B"][finalDisc1Key][finalDisc2Key]) 
            bkgNjetsErr["C"].append(nTotBkgCount_ABCD["nBkgEventsErr_C"][finalDisc1Key][finalDisc2Key]) 
            bkgNjetsErr["D"].append(nTotBkgCount_ABCD["nBkgEventsErr_D"][finalDisc1Key][finalDisc2Key])  

            bkgNjetsPred_A["value"].append(pred_A); bkgNjetsPred_A["error"].append(predUnc_A)

    if args.channel == "0l":
        Njets = [7, 8, 9, 10, 11, 12]
    else:
        Njets = [7, 8, 9, 10, 11]
    totalChi2, ndof = plot_ClosureNjets(bkgNjets["A"], bkgNjetsErr["A"], bkgNjetsPred_A["value"], bkgNjetsPred_A["error"], Njets, args.year, args.channel)


 
if __name__ == '__main__':
    main()
