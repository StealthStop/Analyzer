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
#mpl.__version__
#print mpl.__version__
# mpl version is 2.2.5

from matplotlib.colors import LogNorm
from ROOT import TFile, gROOT, gStyle, TLatex


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
def cal_OptMetricOfBinEdges(significance, closureError):
    #inverseSignificance = (1.0 / significance)
    optimizationMetric  = (closureError)**2 + (1.0 / significance)**2
    return optimizationMetric

# -------------------------------------------------------
# get signal and background histograms' counts
#   -- count both signal and background events separately 
#   -- in each A, B, C, D regions
#   -- put them to the dictionaries
# -------------------------------------------------------
def count_Events_inBinEdges(histBkg, histSig):

    lastXBin = histBkg.GetNbinsX();  lastYBin = histBkg.GetNbinsY()
    nXBins   = range(2, lastXBin+1); nYBins   = range(2, lastYBin+1)

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

        xLowBinEdge = histBkg.GetXaxis().GetBinCenter(xBin)
        xBinKey     = "%.2f"%xLowBinEdge

        if xBinKey not in nTotSigCount_ABCD["nSigEventsErr_A"]:
            nTotSigCount_ABCD["nSigEvents_A"][xBinKey] = {}; nTotSigCount_ABCD["nSigEventsErr_A"][xBinKey] = {} 
            nTotSigCount_ABCD["nSigEvents_B"][xBinKey] = {}; nTotSigCount_ABCD["nSigEventsErr_B"][xBinKey] = {}
            nTotSigCount_ABCD["nSigEvents_C"][xBinKey] = {}; nTotSigCount_ABCD["nSigEventsErr_C"][xBinKey] = {}
            nTotSigCount_ABCD["nSigEvents_D"][xBinKey] = {}; nTotSigCount_ABCD["nSigEventsErr_D"][xBinKey] = {}

        if xBinKey not in nTotBkgCount_ABCD["nBkgEventsErr_A"]:
            nTotBkgCount_ABCD["nBkgEvents_A"][xBinKey] = {}; nTotBkgCount_ABCD["nBkgEventsErr_A"][xBinKey] = {} 
            nTotBkgCount_ABCD["nBkgEvents_B"][xBinKey] = {}; nTotBkgCount_ABCD["nBkgEventsErr_B"][xBinKey] = {}
            nTotBkgCount_ABCD["nBkgEvents_C"][xBinKey] = {}; nTotBkgCount_ABCD["nBkgEventsErr_C"][xBinKey] = {}
            nTotBkgCount_ABCD["nBkgEvents_D"][xBinKey] = {}; nTotBkgCount_ABCD["nBkgEventsErr_D"][xBinKey] = {}

        # loop over the y bins
        for yBin in nYBins:

            yLowBinEdge = histBkg.GetYaxis().GetBinCenter(yBin)
            yBinKey     = "%.2f"%yLowBinEdge

            if yBinKey not in nTotSigCount_ABCD["nSigEventsErr_A"]:
                nTotSigCount_ABCD["nSigEvents_A"][xBinKey][yBinKey] = 0.0; nTotSigCount_ABCD["nSigEventsErr_A"][xBinKey][yBinKey] = 0.0  
                nTotSigCount_ABCD["nSigEvents_B"][xBinKey][yBinKey] = 0.0; nTotSigCount_ABCD["nSigEventsErr_B"][xBinKey][yBinKey] = 0.0
                nTotSigCount_ABCD["nSigEvents_C"][xBinKey][yBinKey] = 0.0; nTotSigCount_ABCD["nSigEventsErr_C"][xBinKey][yBinKey] = 0.0
                nTotSigCount_ABCD["nSigEvents_D"][xBinKey][yBinKey] = 0.0; nTotSigCount_ABCD["nSigEventsErr_D"][xBinKey][yBinKey] = 0.0

            if yBinKey not in nTotBkgCount_ABCD["nBkgEventsErr_A"]:
                nTotBkgCount_ABCD["nBkgEvents_A"][xBinKey][yBinKey] = 0.0; nTotBkgCount_ABCD["nBkgEventsErr_A"][xBinKey][yBinKey] = 0.0 
                nTotBkgCount_ABCD["nBkgEvents_B"][xBinKey][yBinKey] = 0.0; nTotBkgCount_ABCD["nBkgEventsErr_B"][xBinKey][yBinKey] = 0.0 
                nTotBkgCount_ABCD["nBkgEvents_C"][xBinKey][yBinKey] = 0.0; nTotBkgCount_ABCD["nBkgEventsErr_C"][xBinKey][yBinKey] = 0.0
                nTotBkgCount_ABCD["nBkgEvents_D"][xBinKey][yBinKey] = 0.0; nTotBkgCount_ABCD["nBkgEventsErr_D"][xBinKey][yBinKey] = 0.0

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
            
            nTotSigCount_ABCD["nSigEvents_A"][xBinKey][yBinKey] = nSigEvents_A; nTotSigCount_ABCD["nSigEventsErr_A"][xBinKey][yBinKey] = nSigEventsErr_A 
            nTotSigCount_ABCD["nSigEvents_B"][xBinKey][yBinKey] = nSigEvents_B; nTotSigCount_ABCD["nSigEventsErr_B"][xBinKey][yBinKey] = nSigEventsErr_B
            nTotSigCount_ABCD["nSigEvents_C"][xBinKey][yBinKey] = nSigEvents_C; nTotSigCount_ABCD["nSigEventsErr_C"][xBinKey][yBinKey] = nSigEventsErr_C
            nTotSigCount_ABCD["nSigEvents_D"][xBinKey][yBinKey] = nSigEvents_D; nTotSigCount_ABCD["nSigEventsErr_D"][xBinKey][yBinKey] = nSigEventsErr_D

            nTotBkgCount_ABCD["nBkgEvents_A"][xBinKey][yBinKey] = nBkgEvents_A; nTotBkgCount_ABCD["nBkgEventsErr_A"][xBinKey][yBinKey] = nBkgEventsErr_A 
            nTotBkgCount_ABCD["nBkgEvents_B"][xBinKey][yBinKey] = nBkgEvents_B; nTotBkgCount_ABCD["nBkgEventsErr_B"][xBinKey][yBinKey] = nBkgEventsErr_B
            nTotBkgCount_ABCD["nBkgEvents_C"][xBinKey][yBinKey] = nBkgEvents_C; nTotBkgCount_ABCD["nBkgEventsErr_C"][xBinKey][yBinKey] = nBkgEventsErr_C
            nTotBkgCount_ABCD["nBkgEvents_D"][xBinKey][yBinKey] = nBkgEvents_D; nTotBkgCount_ABCD["nBkgEventsErr_D"][xBinKey][yBinKey] = nBkgEventsErr_D

    return nTotSigCount_ABCD, nTotBkgCount_ABCD

# -------------------------------------------------------------------------
# Region by region signal fraction calculation
# look at (SIG + BKG) events in each region, how much signal is in each bin
#   -- in region A : Nsig / (Nbkg + Nsig) 
#   -- in region B : Nsig / (Nbkg + Nsig) 
#   -- in region C : Nsig / (Nbkg + Nsig) 
#   -- in region D : Nsig / (Nbkg + Nsig) 
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
            nSigEvents_A = nTotSigCount_ABCD["nSigEvents_A"][disc1Key][disc2Key]; nBkgEvents_A = nTotBkgCount_ABCD["nBkgEvents_A"][disc1Key][disc2Key] 
            nSigEvents_B = nTotSigCount_ABCD["nSigEvents_B"][disc1Key][disc2Key]; nBkgEvents_B = nTotBkgCount_ABCD["nBkgEvents_B"][disc1Key][disc2Key]
            nSigEvents_C = nTotSigCount_ABCD["nSigEvents_C"][disc1Key][disc2Key]; nBkgEvents_C = nTotBkgCount_ABCD["nBkgEvents_C"][disc1Key][disc2Key]
            nSigEvents_D = nTotSigCount_ABCD["nSigEvents_D"][disc1Key][disc2Key]; nBkgEvents_D = nTotBkgCount_ABCD["nBkgEvents_D"][disc1Key][disc2Key]

            # Region by region signal fraction calculation
            # get some plots based on region by region signal fraction
            nTot_SigBkg_A = nSigEvents_A + nBkgEvents_A
            nTot_SigBkg_B = nSigEvents_B + nBkgEvents_B
            nTot_SigBkg_C = nSigEvents_C + nBkgEvents_C
            nTot_SigBkg_D = nSigEvents_D + nBkgEvents_D
            
            tempSigFracsA = -1.0; tempSigFracsB = -1.0; tempSigFracsC = -1.0; tempSigFracsD = -1.0

            if nTot_SigBkg_A > 0.0: 
                tempSigFracsA = nSigEvents_A / nTot_SigBkg_A
            if nTot_SigBkg_B > 0.0: 
                tempSigFracsB = nSigEvents_B / nTot_SigBkg_B
            if nTot_SigBkg_C > 0.0: 
                tempSigFracsC = nSigEvents_C / nTot_SigBkg_C
            if nTot_SigBkg_D > 0.0: 
                tempSigFracsD = nSigEvents_D / nTot_SigBkg_D

            #sigFracsA.append(float(tempSigFracsA))
            #sigFracsB.append(float(tempSigFracsB))
            #sigFracsC.append(float(tempSigFracsC))
            #sigFracsD.append(float(tempSigFracsD))

            # Total signal (and background) fractions in aech A, B, C, D region
            # get the latest bin edges based on total signal fraction
            nTot_Sig_ABCD = nSigEvents_A + nSigEvents_B + nSigEvents_C + nSigEvents_D
            nTot_Bkg_ABCD = nBkgEvents_A + nBkgEvents_B + nBkgEvents_C + nBkgEvents_D
            
            tempSigTotFracsA = -1.0; tempSigTotFracsB = -1.0; tempSigTotFracsC = -1.0; tempSigTotFracsD = -1.0
            tempBkgTotFracsA = -1.0; tempBkgTotFracsB = -1.0; tempBkgTotFracsC = -1.0; tempBkgTotFracsD = -1.0
            
            tempSigTotFracsA = nSigEvents_A / nTot_Sig_ABCD; tempBkgTotFracsA = nBkgEvents_A / nTot_Bkg_ABCD 
            tempSigTotFracsB = nSigEvents_B / nTot_Sig_ABCD; tempBkgTotFracsB = nBkgEvents_B / nTot_Bkg_ABCD
            tempSigTotFracsC = nSigEvents_C / nTot_Sig_ABCD; tempBkgTotFracsC = nBkgEvents_C / nTot_Bkg_ABCD
            tempSigTotFracsD = nSigEvents_D / nTot_Sig_ABCD; tempBkgTotFracsD = nBkgEvents_D / nTot_Bkg_ABCD

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
        
                # Region by region fraction
                sigFracsA.append(float(tempSigFracsA))
                sigFracsB.append(float(tempSigFracsB))
                sigFracsC.append(float(tempSigFracsC))
                sigFracsD.append(float(tempSigFracsD))

                # Total fraction
                sigTotFracsA.append(float(tempSigTotFracsA)); bkgTotFracsA.append(float(tempBkgTotFracsA)) 
                sigTotFracsB.append(float(tempSigTotFracsB)); bkgTotFracsB.append(float(tempBkgTotFracsB))
                sigTotFracsC.append(float(tempSigTotFracsC)); bkgTotFracsC.append(float(tempBkgTotFracsC))
                sigTotFracsD.append(float(tempSigTotFracsD)); bkgTotFracsD.append(float(tempBkgTotFracsD))

            if (tempBkgTotFracsA > minBkgFrac) and (tempBkgTotFracsB > minBkgFrac) and (tempBkgTotFracsC > minBkgFrac) and (tempBkgTotFracsD > minBkgFrac):                
                tempOptMetric = cal_OptMetricOfBinEdges(tempSignificance, tempClosureErr)

            if tempOptMetric < optMetric:
                finalDisc1Key = disc1Key 
                finalDisc2Key = disc2Key
                significance  = tempSignificance
                closureErr    = tempClosureErr
                optMetric     = tempOptMetric

    # print out the disc1 and disc2 edges
    #print "disc1 (x bin) low bin edges: ", finalDisc1Key
    #print "disc2 (y bin) low bin edges: ", finalDisc2Key

    return finalDisc1Key, finalDisc2Key, significance, closureErr, inverseSignificance, closureErrsList, disc1KeyOut, disc2KeyOut, sigFracsA, sigFracsB, sigFracsC, sigFracsD, sigTotFracsA, sigTotFracsB, sigTotFracsC, sigTotFracsD, bkgTotFracsA, bkgTotFracsB, bkgTotFracsC, bkgTotFracsD 

# --------------------------------------------------------------
# plot significance and closure error as a function of bin edges 
# --------------------------------------------------------------
def plot_SignificanceClosure_BinEdges(nBins, inverseSignificance, closureErrsList, disc1Edges, disc2Edges, c1, c2, minEdge, maxEdge, binWidth, year, model, mass, channel, Njets = -1):

    #nBins = int( (1.0 + binWidth) / binWidth )

    # significance as a function of bin edges
    fig = plt.figure()
    plt.hist2d(disc1Edges, disc2Edges, bins=[nBins, nBins], range=[[0.0, 1.0], [0.0, 1.0]], cmap=plt.cm.jet, weights=np.reciprocal(inverseSignificance), cmin=10e-10, cmax=10.0)
    plt.colorbar()
    ax = plt.gca()
    ax.set_xlabel("Disc. 1 Bin Edge")
    ax.set_ylabel("Disc. 2 Bin Edge")
    ax.text(0.12, 1.05, "CMS",                transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
    ax.text(0.33, 1.04, "Preliminary",        transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
    ax.text(0.99, 1.04, "%s (13 TeV)"%(year), transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right')
    l1 = ml.Line2D([c1, c1], [0.0, 1.0], color="black", linewidth=2, linestyle="dashed")
    l2 = ml.Line2D([0.0, 1.0], [c2, c2], color="black", linewidth=2, linestyle="dashed")
    ax.add_line(l1)
    ax.add_line(l2)

    if Njets == -1: 
        fig.savefig("plots/%s_%s/%s/%s_Sign_vs_Disc1Disc2_%s.pdf"%(model, mass, channel, year, channel), dpi=fig.dpi)
    else:           
        fig.savefig("plots/%s_%s/%s/%s_Sign_vs_Disc1Disc2_Njets%s_%s.pdf"%(model, mass, channel, year, Njets, channel), dpi=fig.dpi)
    plt.close(fig)

    # closure error as a function of bin edges
    fig = plt.figure()
    plt.hist2d(disc1Edges, disc2Edges, bins=[nBins, nBins], range=[[0.0, 1.0], [0.0, 1.0]], cmap=plt.cm.jet, weights=closureErrsList, cmin=10e-10, cmax=2.5)
    plt.colorbar()
    ax = plt.gca()
    ax.set_xlabel("Disc. 1 Bin Edge")
    ax.set_ylabel("Disc. 2 Bin Edge")
    ax.text(0.12, 1.05, "CMS",                transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
    ax.text(0.33, 1.04, "Preliminary",        transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
    ax.text(0.99, 1.04, "%s (13 TeV)"%(year), transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right')
    l1 = ml.Line2D([c1, c1], [0.0, 1.0], color="black", linewidth=2, linestyle="dashed")
    l2 = ml.Line2D([0.0, 1.0], [c2, c2], color="black", linewidth=2, linestyle="dashed")
    ax.add_line(l1)
    ax.add_line(l2)

    if Njets == -1: 
        fig.savefig("plots/%s_%s/%s/%s_CloseErr_vs_Disc1Disc2_%s.pdf"%(model, mass, channel, year, channel), dpi=fig.dpi)
    else:           
        fig.savefig("plots/%s_%s/%s/%s_CloseErr_vs_Disc1Disc2_Njets%s_%s.pdf"%(model, mass, channel, year, Njets, channel), dpi=fig.dpi)
    plt.close(fig)

# --------------------------------------
# plot inverseSignificance vs ClosureErr 
# --------------------------------------
def plot_inverseSignificance_vsClosureErr(significance, closureErr, inverseSignificance, closureErrsList, edges, disc1Edge, disc2Edge, year, model, mass, channel, Njets = -1):

    fig = plt.figure(figsize=(5,5))
    ax = plt.gca()
    ax.text(0.12, 1.05, "CMS",                transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
    ax.text(0.33, 1.04, "Preliminary",        transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
    ax.text(0.99, 1.04, "%s (13 TeV)"%(year), transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right')
    plt.scatter(inverseSignificance, closureErrsList, color='xkcd:black', marker="o", label="1 - Pred./Obs. vs 1 / Significance")

    if significance != 0.0: 
        plt.scatter([1.0 / significance], [closureErr], color='xkcd:red', marker="o", label="Chosen Solution")
    plt.xlabel('1 / Significance')
    plt.xlim(left=0)
    plt.ylabel('|1 - Pred./Obs.|')
    plt.legend(loc='best', frameon=False)
    plt.ylim(bottom=0)
    plt.gca().invert_yaxis()
    plt.text(0.40, 0.85, r"$%.2f < \bf{Disc.\;1\;Edge}$ = %s < %.2f"%(edges[0],disc1Edge,edges[-1]), transform=ax.transAxes, fontsize=8)
    plt.text(0.40, 0.80, r"$%.2f < \bf{Disc.\;2\;Edge}$ = %s < %.2f"%(edges[0],disc2Edge,edges[-1]), transform=ax.transAxes, fontsize=8)

    if Njets == -1:
        fig.savefig("plots/%s_%s/%s/%s_InvSign_vs_CloseErr_%s.pdf"%(model, mass, channel, year, channel), dpi=fig.dpi)        
    else:
        fig.savefig("plots/%s_%s/%s/%s_InvSign_vs_CloseErr_Njets%s_%s.pdf"%(model, mass, channel, year, Njets, channel), dpi=fig.dpi)        

    plt.close(fig)

# ---------------------------------------------
# plot Frac vs Edges
#   -- SigFrac vs Disc1Disc2 in each A, B, C, D
#   -- TotSigFrac vs Disc1Disc2
#   -- BkgFrac vs Disc1Disc2 in each A, B, C, D
#   -- TotBkgFrac vs Disc1Disc2
# ---------------------------------------------
def plot_SigBkgFrac_vsEdges(nBins, sigFracsA, sigFracsB, sigFracsC, sigFracsD, sigTotFracsA, sigTotFracsB, sigTotFracsC, sigTotFracsD, bkgTotFracsA, bkgTotFracsB, bkgTotFracsC, bkgTotFracsD, disc1Edges, disc2Edges, c1, c2, minEdge, maxEdge, binWidth, year, model, mass, channel, Njets = -1):

    #nBins = int( (1.0 + binWidth) / binWidth )

    # SigFracA vs Disc1Disc2
    fig = plt.figure() 
    plt.hist2d(disc1Edges, disc2Edges, bins=[nBins, nBins], range=[[0.0, 1.0], [0.0, 1.0]], cmap=plt.cm.jet, weights=sigFracsA, cmin = 0.00001, cmax = 1.0)
    plt.colorbar()
    ax = plt.gca()
    ax.set_xlabel("Disc. 1 Bin Edge")
    ax.set_ylabel("Disc. 2 Bin Edge")
    ax.text(0.12, 1.05, "CMS",                transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
    ax.text(0.33, 1.04, "Preliminary",        transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
    ax.text(0.99, 1.04, "%s (13 TeV)"%(year), transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right')
    l1 = ml.Line2D([c1, c1], [0.0, 1.0], color="black", linewidth=2, linestyle="dashed"); l2 = ml.Line2D([0.0, 1.0], [c2, c2], color="black", linewidth=2, linestyle="dashed")
    ax.add_line(l1); ax.add_line(l2)

    if Njets == -1: 
        fig.savefig("plots/%s_%s/%s/%s_SigFracA_vs_Disc1Disc2_%s.pdf"%(model, mass, channel, year, channel), dpi=fig.dpi)
    else:           
        fig.savefig("plots/%s_%s/%s/%s_SigFracA_vs_Disc1Disc2_Njets%s_%s.pdf"%(model, mass, channel, year, Njets, channel), dpi=fig.dpi)
    plt.close(fig) 

    # SigFracB vs Disc1Disc2
    fig = plt.figure() 
    plt.hist2d(disc1Edges, disc2Edges, bins=[nBins, nBins], range=[[0.0, 1.0], [0.0, 1.0]], cmap=plt.cm.jet, weights=sigFracsB, cmin = 0.00001, cmax = 1.0)
    plt.colorbar()
    ax = plt.gca()
    ax.set_xlabel("Disc. 1 Bin Edge")
    ax.set_ylabel("Disc. 2 Bin Edge")
    ax.text(0.12, 1.05, "CMS",                transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
    ax.text(0.33, 1.04, "Preliminary",        transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
    ax.text(0.99, 1.04, "%s (13 TeV)"%(year), transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right')
    l1 = ml.Line2D([c1, c1], [0.0, 1.0], color="black", linewidth=2, linestyle="dashed"); l2 = ml.Line2D([0.0, 1.0], [c2, c2], color="black", linewidth=2, linestyle="dashed")
    ax.add_line(l1); ax.add_line(l2)
    
    if Njets == -1: 
        fig.savefig("plots/%s_%s/%s/%s_SigFracB_vs_Disc1Disc2_%s.pdf"%(model, mass, channel, year, channel), dpi=fig.dpi)
    else:           
        fig.savefig("plots/%s_%s/%s/%s_SigFracB_vs_Disc1Disc2_Njets%s_%s.pdf"%(model, mass, channel, year, Njets, channel), dpi=fig.dpi)
    plt.close(fig) 

    # SigFracC vs Disc1Disc2
    fig = plt.figure() 
    plt.hist2d(disc1Edges, disc2Edges, bins=[nBins, nBins], range=[[0.0, 1.0], [0.0, 1.0]], cmap=plt.cm.jet, weights=sigFracsC, cmin = 0.00001, cmax = 1.0)
    plt.colorbar()
    ax = plt.gca()
    ax.set_xlabel("Disc. 1 Bin Edge")
    ax.set_ylabel("Disc. 2 Bin Edge")
    ax.text(0.12, 1.05, "CMS",                transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
    ax.text(0.33, 1.04, "Preliminary",        transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
    ax.text(0.99, 1.04, "%s (13 TeV)"%(year), transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right')
    l1 = ml.Line2D([c1, c1], [0.0, 1.0], color="black", linewidth=2, linestyle="dashed"); l2 = ml.Line2D([0.0, 1.0], [c2, c2], color="black", linewidth=2, linestyle="dashed")
    ax.add_line(l1); ax.add_line(l2)
    
    if Njets == -1: 
        fig.savefig("plots/%s_%s/%s/%s_SigFracC_vs_Disc1Disc2_%s.pdf"%(model, mass, channel, year, channel), dpi=fig.dpi)
    else:           
        fig.savefig("plots/%s_%s/%s/%s_SigFracC_vs_Disc1Disc2_Njets%s_%s.pdf"%(model, mass, channel, year, Njets, channel), dpi=fig.dpi)
    plt.close(fig) 
  
    # SigFracD vs Disc1Disc2
    fig = plt.figure() 
    plt.hist2d(disc1Edges, disc2Edges, bins=[nBins, nBins], range=[[0.0, 1.0], [0.0, 1.0]], cmap=plt.cm.jet, weights=sigFracsD, cmin = 0.00001, cmax = 1.0)
    plt.colorbar()
    ax = plt.gca()
    ax.set_xlabel("Disc. 1 Bin Edge")
    ax.set_ylabel("Disc. 2 Bin Edge")
    ax.text(0.12, 1.05, "CMS",                transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
    ax.text(0.33, 1.04, "Preliminary",        transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
    ax.text(0.99, 1.04, "%s (13 TeV)"%(year), transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right')
    l1 = ml.Line2D([c1, c1], [0.0, 1.0], color="black", linewidth=2, linestyle="dashed"); l2 = ml.Line2D([0.0, 1.0], [c2, c2], color="black", linewidth=2, linestyle="dashed")
    ax.add_line(l1); ax.add_line(l2)
    
    if Njets == -1: 
        fig.savefig("plots/%s_%s/%s/%s_SigFracD_vs_Disc1Disc2_%s.pdf"%(model, mass, channel, year, channel), dpi=fig.dpi)
    else:           
        fig.savefig("plots/%s_%s/%s/%s_SigFracD_vs_Disc1Disc2_Njets%s_%s.pdf"%(model, mass, channel, year, Njets, channel), dpi=fig.dpi)
    plt.close(fig) 
    
    ############################

    # sigTotFracsA vs Disc1Disc2
    fig = plt.figure()
    plt.hist2d(disc1Edges, disc2Edges, bins=[nBins, nBins], range=[[0.0, 1.0], [0.0, 1.0]], cmap=plt.cm.jet, weights=sigTotFracsA, cmin = 0.00001, cmax = 1.0)
    plt.colorbar()
    ax = plt.gca()
    ax.set_xlabel("Disc. 1 Bin Edge")
    ax.set_ylabel("Disc. 2 Bin Edge")
    ax.text(0.12, 1.05, "CMS",                transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
    ax.text(0.33, 1.04, "Preliminary",        transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
    ax.text(0.99, 1.04, "%s (13 TeV)"%(year), transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right')
    l1 = ml.Line2D([c1, c1], [0.0, 1.0], color="black", linewidth=2, linestyle="dashed"); l2 = ml.Line2D([0.0, 1.0], [c2, c2], color="black", linewidth=2, linestyle="dashed")
    ax.add_line(l1); ax.add_line(l2)

    if Njets == -1:
        fig.savefig("plots/%s_%s/%s/%s_SigTotFracA_vs_Disc1Disc2_%s.pdf"%(model, mass, channel, year, channel), dpi=fig.dpi)
    else:
        fig.savefig("plots/%s_%s/%s/%s_SigTotFracA_vs_Disc1Disc2_Njets%s_%s.pdf"%(model, mass, channel, year, Njets, channel), dpi=fig.dpi)
    plt.close(fig)
 
    # sigTotFracsB vs Disc1Disc2 
    fig = plt.figure()
    plt.hist2d(disc1Edges, disc2Edges, bins=[nBins, nBins], range=[[0.0, 1.0], [0.0, 1.0]], cmap=plt.cm.jet, weights=sigTotFracsB, cmin = 0.00001, cmax = 1.0)
    plt.colorbar()
    ax = plt.gca()
    ax.set_xlabel("Disc. 1 Bin Edge")
    ax.set_ylabel("Disc. 2 Bin Edge")
    ax.text(0.12, 1.05, "CMS",                transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
    ax.text(0.33, 1.04, "Preliminary",        transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
    ax.text(0.99, 1.04, "%s (13 TeV)"%(year), transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right')
    l1 = ml.Line2D([c1, c1], [0.0, 1.0], color="black", linewidth=2, linestyle="dashed"); l2 = ml.Line2D([0.0, 1.0], [c2, c2], color="black", linewidth=2, linestyle="dashed")
    ax.add_line(l1); ax.add_line(l2)

    if Njets == -1:
        fig.savefig("plots/%s_%s/%s/%s_SigTotFracB_vs_Disc1Disc2_%s.pdf"%(model, mass, channel, year, channel), dpi=fig.dpi)
    else:
        fig.savefig("plots/%s_%s/%s/%s_SigTotFracB_vs_Disc1Disc2_Njets%s_%s.pdf"%(model, mass, channel, year, Njets, channel), dpi=fig.dpi)
    plt.close(fig)

    # sigTotFracsC vs Disc1Disc2
    fig = plt.figure()
    plt.hist2d(disc1Edges, disc2Edges, bins=[nBins, nBins], range=[[0.0, 1.0], [0.0, 1.0]], cmap=plt.cm.jet, weights=sigTotFracsC, cmin = 0.00001, cmax = 1.0)
    plt.colorbar()
    ax = plt.gca()
    ax.set_xlabel("Disc. 1 Bin Edge")
    ax.set_ylabel("Disc. 2 Bin Edge")
    ax.text(0.12, 1.05, "CMS",                transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
    ax.text(0.33, 1.04, "Preliminary",        transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
    ax.text(0.99, 1.04, "%s (13 TeV)"%(year), transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right')
    l1 = ml.Line2D([c1, c1], [0.0, 1.0], color="black", linewidth=2, linestyle="dashed"); l2 = ml.Line2D([0.0, 1.0], [c2, c2], color="black", linewidth=2, linestyle="dashed")
    ax.add_line(l1); ax.add_line(l2)

    if Njets == -1:
        fig.savefig("plots/%s_%s/%s/%s_SigTotFracC_vs_Disc1Disc2_%s.pdf"%(model, mass, channel, year, channel), dpi=fig.dpi)
    else:
        fig.savefig("plots/%s_%s/%s/%s_SigTotFracC_vs_Disc1Disc2_Njets%s_%s.pdf"%(model, mass, channel, year, Njets, channel), dpi=fig.dpi)
    plt.close(fig)
 
    # sigTotFracsD vs Disc1Disc2 
    fig = plt.figure()
    plt.hist2d(disc1Edges, disc2Edges, bins=[nBins, nBins], range=[[0.0, 1.0], [0.0, 1.0]], cmap=plt.cm.jet, weights=sigTotFracsD, cmin = 0.00001, cmax = 1.0)
    plt.colorbar()
    ax = plt.gca()
    ax.set_xlabel("Disc. 1 Bin Edge")
    ax.set_ylabel("Disc. 2 Bin Edge")
    ax.text(0.12, 1.05, "CMS",                transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
    ax.text(0.33, 1.04, "Preliminary",        transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
    ax.text(0.99, 1.04, "%s (13 TeV)"%(year), transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right')
    l1 = ml.Line2D([c1, c1], [0.0, 1.0], color="black", linewidth=2, linestyle="dashed"); l2 = ml.Line2D([0.0, 1.0], [c2, c2], color="black", linewidth=2, linestyle="dashed")
    ax.add_line(l1); ax.add_line(l2)

    if Njets == -1:
        fig.savefig("plots/%s_%s/%s/%s_SigTotFracD_vs_Disc1Disc2_%s.pdf"%(model, mass, channel, year, channel), dpi=fig.dpi)
    else:
        fig.savefig("plots/%s_%s/%s/%s_SigTotFracD_vs_Disc1Disc2_Njets%s_%s.pdf"%(model, mass, channel, year, Njets, channel), dpi=fig.dpi)
    plt.close(fig)  

    # bkgTotFracsA vs Disc1Disc2
    fig = plt.figure()
    plt.hist2d(disc1Edges, disc2Edges, bins=[nBins, nBins], range=[[0.0, 1.0], [0.0, 1.0]], cmap=plt.cm.jet, weights=bkgTotFracsA, cmin = 0.00001, cmax = 1.0)
    plt.colorbar()
    ax = plt.gca()
    ax.set_xlabel("Disc. 1 Bin Edge")
    ax.set_ylabel("Disc. 2 Bin Edge")
    ax.text(0.12, 1.05, "CMS",                transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
    ax.text(0.33, 1.04, "Preliminary",        transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
    ax.text(0.99, 1.04, "%s (13 TeV)"%(year), transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right')
    l1 = ml.Line2D([c1, c1], [0.0, 1.0], color="black", linewidth=2, linestyle="dashed"); l2 = ml.Line2D([0.0, 1.0], [c2, c2], color="black", linewidth=2, linestyle="dashed")
    ax.add_line(l1); ax.add_line(l2)

    if Njets == -1:
        fig.savefig("plots/%s_%s/%s/%s_BkgTotFracA_vs_Disc1Disc2_%s.pdf"%(model, mass, channel, year, channel), dpi=fig.dpi)
    else:
        fig.savefig("plots/%s_%s/%s/%s_BkgTotFracA_vs_Disc1Disc2_Njets%s_%s.pdf"%(model, mass, channel, year, Njets, channel), dpi=fig.dpi)
    plt.close(fig)

    # bkgTotFracsB vs Disc1Disc2 
    fig = plt.figure()
    plt.hist2d(disc1Edges, disc2Edges, bins=[nBins, nBins], range=[[0.0, 1.0], [0.0, 1.0]], cmap=plt.cm.jet, weights=bkgTotFracsB, cmin = 0.00001, cmax = 1.0)
    plt.colorbar()
    ax = plt.gca()
    ax.set_xlabel("Disc. 1 Bin Edge")
    ax.set_ylabel("Disc. 2 Bin Edge")
    ax.text(0.12, 1.05, "CMS",                transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
    ax.text(0.33, 1.04, "Preliminary",        transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
    ax.text(0.99, 1.04, "%s (13 TeV)"%(year), transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right')
    l1 = ml.Line2D([c1, c1], [0.0, 1.0], color="black", linewidth=2, linestyle="dashed"); l2 = ml.Line2D([0.0, 1.0], [c2, c2], color="black", linewidth=2, linestyle="dashed")
    ax.add_line(l1); ax.add_line(l2)

    if Njets == -1:
        fig.savefig("plots/%s_%s/%s/%s_BkgTotFracB_vs_Disc1Disc2_%s.pdf"%(model, mass, channel, year, channel), dpi=fig.dpi)
    else:
        fig.savefig("plots/%s_%s/%s/%s_BkgTotFracB_vs_Disc1Disc2_Njets%s_%s.pdf"%(model, mass, channel, year, Njets, channel), dpi=fig.dpi)
    plt.close(fig)

    # bkgTotFracsC vs Disc1Disc2
    fig = plt.figure()
    plt.hist2d(disc1Edges, disc2Edges, bins=[nBins, nBins], range=[[0.0, 1.0], [0.0, 1.0]], cmap=plt.cm.jet, weights=bkgTotFracsC, cmin = 0.00001, cmax = 1.0)
    plt.colorbar()
    ax = plt.gca()
    ax.set_xlabel("Disc. 1 Bin Edge")
    ax.set_ylabel("Disc. 2 Bin Edge")
    ax.text(0.12, 1.05, "CMS",                transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
    ax.text(0.33, 1.04, "Preliminary",        transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
    ax.text(0.99, 1.04, "%s (13 TeV)"%(year), transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right')
    l1 = ml.Line2D([c1, c1], [0.0, 1.0], color="black", linewidth=2, linestyle="dashed"); l2 = ml.Line2D([0.0, 1.0], [c2, c2], color="black", linewidth=2, linestyle="dashed")
    ax.add_line(l1); ax.add_line(l2)

    if Njets == -1:
        fig.savefig("plots/%s_%s/%s/%s_BkgTotFracC_vs_Disc1Disc2_%s.pdf"%(model, mass, year, channel), dpi=fig.dpi)
    else:
        fig.savefig("plots/%s_%s/%s/%s_BkgTotFracC_vs_Disc1Disc2_Njets%s_%s.pdf"%(model, mass, channel, year, Njets, channel), dpi=fig.dpi)
    plt.close(fig)

    # bkgTotFracsD vs Disc1Disc2 
    fig = plt.figure()
    plt.hist2d(disc1Edges, disc2Edges, bins=[nBins, nBins], range=[[0.0, 1.0], [0.0, 1.0]],cmap=plt.cm.jet, weights=bkgTotFracsD, cmin = 0.00001, cmax = 1.0)
    plt.colorbar()
    ax = plt.gca()
    ax.set_xlabel("Disc. 1 Bin Edge")
    ax.set_ylabel("Disc. 2 Bin Edge")
    ax.text(0.12, 1.05, "CMS",                transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
    ax.text(0.33, 1.04, "Preliminary",        transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
    ax.text(0.99, 1.04, "%s (13 TeV)"%(year), transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right')
    l1 = ml.Line2D([c1, c1], [0.0, 1.0], color="black", linewidth=2, linestyle="dashed"); l2 = ml.Line2D([0.0, 1.0], [c2, c2], color="black", linewidth=2, linestyle="dashed")
    ax.add_line(l1); ax.add_line(l2)

    if Njets == -1:
        fig.savefig("plots/%s_%s/%s/%s_BkgTotFracD_vs_Disc1Disc2_%s.pdf"%(model, mass, year, channel), dpi=fig.dpi)
    else:
        fig.savefig("plots/%s_%s/%s/%s_BkgTotFracD_vs_Disc1Disc2_Njets%s_%s.pdf"%(model, mass, channel, year, Njets, channel), dpi=fig.dpi)
    plt.close(fig)

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
def plot_ClosureNjets(bkg, bkgUnc, bkgPred, bkgPredUnc, Njets, year, model, mass, channel):

    binCenters = []; xErr    = []; abcdPull = []; abcdError = []
    unc        = []; predUnc = []; obs      = []; pred      = []

    totalChi2 = 0.0; ndof = 0

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
        
    fig = plt.figure(figsize=(5,5))
    ax  = fig.add_subplot(111)
    fig.subplots_adjust(hspace=0)

    lowerNjets  = Njets[0]
    higherNjets = Njets[-1]

    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)

    ax1 = fig.add_subplot(3,1,(1,2))
    ax1.set_yscale("log")
    ax1.set_xlim([lowerNjets-0.5,higherNjets+0.5])
    ax1.text(0.05, 0.1,  "$\chi^2$ / ndof = %3.2f"%(totalChi2/float(ndof)), horizontalalignment="left", verticalalignment="center", transform=ax1.transAxes, fontsize=10)
    ax1.text(0.12, 1.05, "CMS",                transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right') 
    ax1.text(0.33, 1.04, "Preliminary",        transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
    ax1.text(0.99, 1.04, "%s (13 TeV)"%(year), transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right') 
    ax1.set_ylabel('Unweighted Event Counts')
    ax1.errorbar(binCenters, pred, yerr=predUnc, label="Observed",  xerr=xErr, fmt='', color="red",   lw=0, elinewidth=2, marker="o", markerfacecolor="red",   markersize=4.0)
    ax1.errorbar(binCenters, obs,  yerr=unc,     label="Predicted", xerr=xErr, fmt='', color="black", lw=0, elinewidth=2, marker="o", markerfacecolor="black", markersize=4.0)

    ax2 = fig.add_subplot(3,1,3)
    ax2.set_xlim([lowerNjets-0.5,higherNjets+0.5])
    ax2.errorbar(binCenters, abcdError, yerr=None, xerr=xErr, fmt='', color="blue",  lw=0, elinewidth=2, marker="o", markerfacecolor="blue", markersize=4.0)
    ax2.axhline(y=0.0, color="black", linestyle="dashed", lw=1)
    ax2.grid(axis="y", color="black", linestyle="dashed", which="both")
    ax2.set_xlabel('Number of jets')
    ax2.set_ylabel('1 - Pred./Obs.', fontsize="small")
    ax2.set_ylim([-1.6, 1.6])
 
    plt.xticks(Njets)

    ax1.legend(loc='best', frameon=False)
 
    fig.savefig("plots/%s_%s/%s/%s_Njets_Region_A_PredVsActual_%s.pdf" %(model, mass, channel, year, channel))

    plt.close(fig)

    return totalChi2, ndof


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

    os.system("rm BinEdges_%s_%s_%s.txt" %(args.model, args.mass, args.channel))

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

        # put the latest bin edges to txt file
        d = open("BinEdges_%s_%s_%s.txt" %(args.model, args.mass, args.channel), "a")
        d.write("%s bin edges: \n" %(njet))
        d.write("x bin edges: %s \n" %(finalDisc1Key))
        d.write("y bin edges: %s \n" %(finalDisc2Key))
        d.write("\n")

        # --------------------------------------------------------------------
        # make significance and closure error as a function of bin edges plots
        # --------------------------------------------------------------------
        minEdge  = histSig.GetXaxis().GetBinLowEdge(1) 
        maxEdge  = histSig.GetXaxis().GetBinUpEdge(histBkg.GetNbinsX())
        binWidth = histSig.GetXaxis().GetBinWidth(1)
        nBins    = histSig.GetNbinsX()
        plot_SignificanceClosure_BinEdges(nBins, inverseSignificance, closureErrsList, disc1KeyOut, disc2KeyOut, float(finalDisc1Key), float(finalDisc2Key), minEdge, maxEdge, binWidth, args.year, args.model, args.mass, args.channel, njet)

        # --------------------------------------
        # plot inverseSignificance vs ClosureErr 
        # --------------------------------------
        edges = np.arange(minEdge, maxEdge, binWidth)
        plot_inverseSignificance_vsClosureErr(significance, closureErr, inverseSignificance, closureErrsList, edges, finalDisc1Key, finalDisc2Key, args.year, args.model, args.mass, args.channel, njet)
        
        # ---------------------------------------------
        # plot Frac vs Edges
        #   -- SigFrac vs Disc1Disc2 in each A, B, C, D
        #   -- TotSigFrac vs Disc1Disc2
        #   -- BkgFrac vs Disc1Disc2 in each A, B, C, D
        #   -- TotBkgFrac vs Disc1Disc2
        # ---------------------------------------------
        plot_SigBkgFrac_vsEdges(nBins, sigFracsA, sigFracsB, sigFracsC, sigFracsD, sigTotFracsA, sigTotFracsB, sigTotFracsC, sigTotFracsD, bkgTotFracsA, bkgTotFracsB, bkgTotFracsC, bkgTotFracsD, disc1KeyOut, disc2KeyOut, float(finalDisc1Key), float(finalDisc2Key), minEdge, maxEdge, binWidth, args.year, args.model, args.mass, args.channel, njet)

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

            closure, closureUnc, pred_A, predUnc_A = cal_simpleClosureABCD( nTotBkgCount_ABCD["nBkgEvents_A"][finalDisc1Key][finalDisc2Key],    nTotBkgCount_ABCD["nBkgEvents_B"][finalDisc1Key][finalDisc2Key], 
                                                                            nTotBkgCount_ABCD["nBkgEvents_C"][finalDisc1Key][finalDisc2Key],    nTotBkgCount_ABCD["nBkgEvents_D"][finalDisc1Key][finalDisc2Key], 
                                                                            nTotBkgCount_ABCD["nBkgEventsErr_A"][finalDisc1Key][finalDisc2Key], nTotBkgCount_ABCD["nBkgEventsErr_B"][finalDisc1Key][finalDisc2Key], 
                                                                            nTotBkgCount_ABCD["nBkgEventsErr_C"][finalDisc1Key][finalDisc2Key], nTotBkgCount_ABCD["nBkgEventsErr_D"][finalDisc1Key][finalDisc2Key] )

            # ---------------------
            # make the closure plot
            # ---------------------
            sigNjets["A"].append(nTotSigCount_ABCD["nSigEvents_A"][finalDisc1Key][finalDisc2Key]); sigNjetsErr["A"].append(nTotSigCount_ABCD["nSigEventsErr_A"][finalDisc1Key][finalDisc2Key]) 
            sigNjets["B"].append(nTotSigCount_ABCD["nSigEvents_B"][finalDisc1Key][finalDisc2Key]); sigNjetsErr["B"].append(nTotSigCount_ABCD["nSigEventsErr_B"][finalDisc1Key][finalDisc2Key])
            sigNjets["C"].append(nTotSigCount_ABCD["nSigEvents_C"][finalDisc1Key][finalDisc2Key]); sigNjetsErr["C"].append(nTotSigCount_ABCD["nSigEventsErr_C"][finalDisc1Key][finalDisc2Key])
            sigNjets["D"].append(nTotSigCount_ABCD["nSigEvents_D"][finalDisc1Key][finalDisc2Key]); sigNjetsErr["D"].append(nTotSigCount_ABCD["nSigEventsErr_D"][finalDisc1Key][finalDisc2Key])

            bkgNjets["A"].append(nTotBkgCount_ABCD["nBkgEvents_A"][finalDisc1Key][finalDisc2Key]); bkgNjetsErr["A"].append(nTotBkgCount_ABCD["nBkgEventsErr_A"][finalDisc1Key][finalDisc2Key]) 
            bkgNjets["B"].append(nTotBkgCount_ABCD["nBkgEvents_B"][finalDisc1Key][finalDisc2Key]); bkgNjetsErr["B"].append(nTotBkgCount_ABCD["nBkgEventsErr_B"][finalDisc1Key][finalDisc2Key])
            bkgNjets["C"].append(nTotBkgCount_ABCD["nBkgEvents_C"][finalDisc1Key][finalDisc2Key]); bkgNjetsErr["C"].append(nTotBkgCount_ABCD["nBkgEventsErr_C"][finalDisc1Key][finalDisc2Key])
            bkgNjets["D"].append(nTotBkgCount_ABCD["nBkgEvents_D"][finalDisc1Key][finalDisc2Key]); bkgNjetsErr["D"].append(nTotBkgCount_ABCD["nBkgEventsErr_D"][finalDisc1Key][finalDisc2Key])

            bkgNjetsPred_A["value"].append(pred_A); bkgNjetsPred_A["error"].append(predUnc_A)

    if args.channel == "0l":
        Njets = [7, 8, 9, 10, 11, 12]
    else:
        Njets = [7, 8, 9, 10, 11]
    
    totalChi2, ndof = plot_ClosureNjets(bkgNjets["A"], bkgNjetsErr["A"], bkgNjetsPred_A["value"], bkgNjetsPred_A["error"], Njets, args.year, args.model, args.mass, args.channel)


    d.close()
 
if __name__ == '__main__':
    main()
