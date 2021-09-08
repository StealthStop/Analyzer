import math
import numpy as np

import matplotlib as mpl
mpl.use('Agg')

import matplotlib.lines as ml
import matplotlib.pyplot as plt

from matplotlib.colors import LogNorm

from DoubleDisCo_Helpers import *

# Class that makes all the fancy plots
# Information needed to make the plots
# is not owned by this class, just passed in on the call
class Common_Calculations_Plotters:

    def __init__(self, tag, year, model, mass, channel):
        self.tag     = tag
        self.year    = year
        self.model   = model
        self.mass    = mass
        self.channel = channel

        self.outputDir = "plots_%s/%s_%s/%s/"%(self.tag, self.model, self.mass, self.channel)

    # -----------------
    # plot all closures
    # -----------------
    def plot_ClosureNjets(self, bkg, bkgPred, Njets, region = '', tag = ''):
        
        x         = []; xUnc         = []
        obs       = []; obsUnc       = [] 
        pred      = []; predUnc      = []
        abcdError = []; abcdErrorUnc = []
        
        totalChi2 = 0.0; ndof = 0
    
        for i in range(0, len(Njets)):
    
            if bkg[i][1] != 0.0:
    
                x.append(float(Njets[i]));  xUnc.append(0.5)
                pred.append(bkgPred[i][0]); predUnc.append(bkgPred[i][1])
                obs.append(bkg[i][0]);      obsUnc.append(bkg[i][1])

                pull         = (bkgPred[i][0] - bkg[i][0]) / bkg[i][1]
                closureError    = 1.0 - bkgPred[i][0] / bkg[i][0]
                closureErrorUnc = ((bkgPred[i][1] / bkg[i][0])**2.0 + (bkg[i][1] * bkgPred[i][0] / bkg[i][0]**2.0)**2.0)**0.5

                abcdError.append(closureError); abcdErrorUnc.append(closureErrorUnc)
                totalChi2    += pull ** 2.0
                ndof         += 1
    
        fig = plt.figure(figsize=(5, 5))
        ax  = fig.add_subplot(111)
        fig.subplots_adjust(hspace=0, left=0.15, right=0.95, top=0.94)
    
        lowerNjets  = float(Njets[0])
        higherNjets = float(Njets[(-1)])
    
        ax.spines['top'].set_color('none')
        ax.spines['bottom'].set_color('none')
        ax.spines['left'].set_color('none')
        ax.spines['right'].set_color('none')
        ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
    
        ax1 = fig.add_subplot(3, 1, (1, 2))
        fig.subplots_adjust(left=0.15, right=0.95)
        ax1.set_yscale('log')
        ax1.set_xlim([lowerNjets - 0.5, higherNjets + 0.5])
        ax1.text(0.05, 0.1, '$\\chi^2$ / ndof = %3.2f' % (totalChi2 / float(ndof)), horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes, fontsize=10)
        ax1.text(0.05, 0.25,  '%s Metric'%(self.tag),     horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes, fontsize=14, fontweight='bold')
        ax1.text(0.16, 1.065, 'CMS',                      transform=ax.transAxes, fontsize=20, fontweight='bold',   va='top', ha='right')
        ax1.text(0.50, 1.055, 'Preliminary',              transform=ax.transAxes, fontsize=16, fontstyle='italic',  va='top', ha='right')
        ax1.text(1.0,  1.055, '%s (13 TeV)' %(self.year), transform=ax.transAxes, fontsize=14, fontweight='normal', va='top', ha='right')
        ax1.set_ylabel('Unweighted Event Counts')
        ax1.errorbar(x, pred, yerr=predUnc, label='Predicted', xerr=xUnc, fmt='', capsize=0, color='red',   lw=0, elinewidth=2, marker='o', markeredgecolor='red',   markerfacecolor='red',   markersize=5.0)
        ax1.errorbar(x, obs,  yerr=obsUnc,  label='Observed',  xerr=xUnc, fmt='', capsize=0, color='black', lw=0, elinewidth=2, marker='o', markeredgecolor='black', markerfacecolor='black', markersize=5.0)
    
        ax2 = fig.add_subplot(3, 1, 3)
        ax2.set_xlim([lowerNjets - 0.5, higherNjets + 0.5])
        ax2.errorbar(x, abcdError, yerr=abcdErrorUnc, xerr=xUnc, fmt='', capsize=0, color='blue', lw=0, elinewidth=2, marker='o', markeredgecolor='blue', markerfacecolor='blue', markersize=5.0)
        ax2.axhline(y=0.0, color='black', linestyle='dashed', lw=1)
        ax2.grid(axis='y', color='black', linestyle='dashed', which='both')
        ax2.set_xlabel('Number of jets')
        ax2.set_ylabel('1 - Pred./Obs.')
        ax2.set_ylim([-1.4, 1.4])
    
        plt.xticks([int(Njet) for Njet in Njets])
    
        ax1.legend(loc='best', numpoints=1, frameon=False)
    
        fig.savefig('%s/%s_Njets_Region_A_PredVsActual_%s_%s_%s.pdf' % (self.outputDir, self.year, region, self.channel, tag))
    
        plt.close(fig)
    
    # ----------------------
    # make all closure plots
    # ----------------------    
    def make_allClosures(self, edgesPerNjets=None, bkgEventsPerNjets=None, sigEventsPerNjets=None, Njets=None, region = "", tag = ""):
    
        evtsNjets     = []
        evtsNjetsPred = []
    
        for Njet in Njets:
    
            if edgesPerNjets[Njet][region][0] == -1.0 or edgesPerNjets[Njet][region][1] == -1.0:
    
                evtsNjets.append((0.0, 0.0))
                evtsNjetsPred.append((0.0, 0.0))
    
            else:
    
                evtsA = bkgEventsPerNjets[Njet][region]["A"][0]; evtsAunc = bkgEventsPerNjets[Njet][region]["A"][1]**2.0
                evtsB = bkgEventsPerNjets[Njet][region]["B"][0]; evtsBunc = bkgEventsPerNjets[Njet][region]["B"][1]**2.0
                evtsC = bkgEventsPerNjets[Njet][region]["C"][0]; evtsCunc = bkgEventsPerNjets[Njet][region]["C"][1]**2.0
                evtsD = bkgEventsPerNjets[Njet][region]["D"][0]; evtsDunc = bkgEventsPerNjets[Njet][region]["D"][1]**2.0
    
                if sigEventsPerNjets != None:
                    evtsA += sigEventsPerNjets[Njet][region]["A"][0]; evtsAunc += sigEventsPerNjets[Njet][region]["A"][1]**2.0
                    evtsB += sigEventsPerNjets[Njet][region]["B"][0]; evtsBunc += sigEventsPerNjets[Njet][region]["B"][1]**2.0
                    evtsC += sigEventsPerNjets[Njet][region]["C"][0]; evtsCunc += sigEventsPerNjets[Njet][region]["C"][1]**2.0
                    evtsD += sigEventsPerNjets[Njet][region]["D"][0]; evtsDunc += sigEventsPerNjets[Njet][region]["D"][1]**2.0
    
                pred_A, predUnc_A = cal_SimpleABCD(evtsA, evtsB, evtsC, evtsD, evtsAunc**0.5, evtsBunc**0.5, evtsCunc**0.5, evtsDunc**0.5)
    
                evtsNjets.append((evtsA, evtsAunc**0.5))
                evtsNjetsPred.append((pred_A, predUnc_A))
    
        self.plot_ClosureNjets(np.array(evtsNjets), np.array(evtsNjetsPred), Njets, region, tag)

    # ----------------------------------------------------------------
    # plot whichever variable as a function of the choice of bin edges
    # ----------------------------------------------------------------
    def plot_Var_vsDisc1Disc2(self, var, edges, c1, c2, minEdge, maxEdge, edgeWidth, cmax, vmax, Njets = -1, name = "", region = "", tag = ""):

        nBins = math.ceil((1.0 + edgeWidth)/edgeWidth)

        fig = plt.figure() 
        fig.subplots_adjust(top = 0.94)
        plt.hist2d(edges[:,0], edges[:,1], bins=[nBins, nBins], range=[[-edgeWidth/2.0, 1+edgeWidth/2.0], [-edgeWidth/2.0, 1+edgeWidth/2.0]], cmap=plt.cm.jet, weights=var, cmin=10e-10, cmax=cmax, vmin = 0.0, vmax = vmax)
        plt.colorbar()
        ax = plt.gca()
        ax.set_ylabel("Disc. 2 Bin Edge"); ax.set_xlabel("Disc. 1 Bin Edge")
        ax.text(0.16, 1.065, 'CMS',                     transform=ax.transAxes, fontsize=20, fontweight='bold',   va='top', ha='right')
        ax.text(0.50, 1.055, 'Preliminary',             transform=ax.transAxes, fontsize=16, fontstyle='italic',  va='top', ha='right')
        ax.text(1.0,  1.055, '%s (13 TeV)' % self.year, transform=ax.transAxes, fontsize=14, fontweight='normal', va='top', ha='right')

        l1 = ml.Line2D([c1, c1], [0.0, 1.0], color="black", linewidth=2, linestyle="dashed"); l2 = ml.Line2D([0.0, 1.0], [c2, c2], color="black", linewidth=2, linestyle="dashed")
        ax.add_line(l1); ax.add_line(l2)
        fig.tight_layout()

        fig.savefig(self.outputDir+"/%s_%s_vs_Disc1Disc2_Njets%s_%s_%s.pdf"%(self.year, name, Njets, region, tag), dpi=fig.dpi)

        plt.close(fig)

    # --------------------------------------
    # plot inverseSignificance vs ClosureErr
    # -------------------------------------- 
    def plot_inverseSignificance_vsClosureErr(self, finalSignificance, finalClosureErr, significances, closureErrs, edges, finalDiscEdges, Njets = -1, name=""):

        fig = plt.figure(figsize=(5, 5))
        ax = plt.gca()
        ax.text(0.12, 1.05, 'CMS',                     transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
        ax.text(0.33, 1.04, 'Preliminary',             transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
        ax.text(0.99, 1.04, '%s (13 TeV)' % self.year, transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right')
        plt.scatter(np.reciprocal(significances), closureErrs, color='grey', marker='o', label='1 - Pred./Obs. vs 1 / Significance')

        if finalSignificance != 0.0:
            plt.scatter([1.0 / finalSignificance[0]], [finalClosureErr[0]], color='red', marker='o', label='Chosen Solution')
        plt.xlabel('1 / Significance')
        plt.xlim(left=0)
        plt.ylabel('|1 - Pred./Obs.|')
        plt.legend(loc='best', numpoints=1, frameon=False)
        plt.ylim(bottom=0)
        plt.gca().invert_yaxis()
        plt.text(0.4, 0.85, '$%.2f < \\bf{Disc.\\;1\\;Edge}$ = %s < %.2f' % (edges[0][0], finalDiscEdges[0], edges[-1][0]), transform=ax.transAxes, fontsize=8)
        plt.text(0.4, 0.8, '$%.2f < \\bf{Disc.\\;2\\;Edge}$ = %s < %.2f' % (edges[0][1],  finalDiscEdges[1], edges[-1][1]), transform=ax.transAxes, fontsize=8)

        fig.savefig('%s/%s_InvSign_vs_ClosureErr_Njets%s_%s.pdf' % (self.outputDir, self.year, Njets, name), dpi=fig.dpi)

        plt.close(fig)

    # -------------------------------------------------
    # plot Disc1 vs Disc2 in each Region with all edges
    # -------------------------------------------------
    def plot_Disc1VsDisc2(self, hist, allRegionsEdges, Njets = -1, tag = "", name="", col1="", col2=""):

        nBins   = hist.GetXaxis().GetNbins()
        nEvents = []; disc1Edges = []; disc2Edges = []

        for x in range(0, hist.GetXaxis().GetNbins()+1):
        
            for y in range(0, hist.GetYaxis().GetNbins()+1):
        
                nEvents.append(hist.GetBinContent(x,y))
                disc1Edges.append(hist.GetXaxis().GetBinCenter(x))
                disc2Edges.append(hist.GetYaxis().GetBinCenter(y))

        fig = plt.figure()
        plt.hist2d(disc1Edges, disc2Edges, bins=[nBins, nBins], range=[[0.0, 1.0], [0.0, 1.0]], cmap=plt.cm.jet, weights=nEvents, cmin=10e-10)
        plt.colorbar()
        ax = plt.gca()
        ax.set_xlabel('Disc. 1')
        ax.set_ylabel('Disc. 2')
        ax.text(0.12, 1.05, 'CMS',                     transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
        ax.text(0.33, 1.04, 'Preliminary',             transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
        ax.text(0.99, 1.04, '%s (13 TeV)' % self.year, transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right')
        l1 = ml.Line2D([float(allRegionsEdges["ABCD"][0]), float(allRegionsEdges["ABCD"][0])], [0.0, 1.0],   color="darkviolet", linewidth=4, linestyle='solid')
        l2 = ml.Line2D([0.0, 1.0], [float(allRegionsEdges["ABCD"][1]), float(allRegionsEdges["ABCD"][1])],   color="darkviolet", linewidth=4, linestyle='solid')
        l3 = ml.Line2D([float(allRegionsEdges["Val_bdEF"][0]), float(allRegionsEdges["Val_bdEF"][0])], [0.0, 1.0],   color=col1, linewidth=4, linestyle='solid')
        l4 = ml.Line2D([0.0, 1.0], [float(allRegionsEdges["Val_cdiGH"][1]), float(allRegionsEdges["Val_cdiGH"][1])], color=col2, linewidth=4, linestyle='solid')
        ax.add_line(l1)
        ax.add_line(l2)
        ax.add_line(l3)
        ax.add_line(l4)

        fig.tight_layout()
        fig.savefig('%s/%s_Disc1VsDisc2_%s_Njets%s_%s.pdf' % (self.outputDir, self.year, tag, Njets, name), dpi=fig.dpi)

        plt.close(fig)

    # ----------------------------------------------------------------------
    # plot variable vs disc as 1D
    #   -- Closure, Significance, weightedEventCounts vs disc1, disc2 slices
    # ----------------------------------------------------------------------
    def plot_VarVsDisc(self, var, edges, edgeWidth, ylim = -1.0, ylabel = "", name = "", disc = -1, Njets = -1, channel="", region="", tag=""):

        x25 = []; x50 = []; x75 = []; xDiag = []
        y25 = []; y50 = []; y75 = []; yDiag = []

        y25unc = []; y50unc = []; y75unc = []; yDiagUnc = []

        for i in range(0, len(var[:,0])):

            if abs(edges[i][disc-1] - 0.25) < 0.01:
               x25.append(edges[i][2-disc])
               y25.append(var[i][0])
               y25unc.append(var[i][1])
           
            elif abs(edges[i][disc-1] - 0.50) < 0.01:
               x50.append(edges[i][2-disc])
               y50.append(var[i][0])
               y50unc.append(var[i][1])
           
            elif abs(edges[i][disc-1] - 0.75) < 0.01: 
               x75.append(edges[i][2-disc])
               y75.append(var[i][0])
               y75unc.append(var[i][1])

            if edges[i][0] == edges[i][1]: 
               xDiag.append(edges[i][2-disc])
               yDiag.append(var[i][0])
               yDiagUnc.append(var[i][1])

        fig = plt.figure(figsize=(5, 5))
        ax = plt.gca()
        
        xWidths25   = [edgeWidth for i in range(0, len(x25))]
        xWidths50   = [edgeWidth for i in range(0, len(x50))]
        xWidths75   = [edgeWidth for i in range(0, len(x75))]
        xWidthsDiag = [edgeWidth for i in range(0, len(xDiag))]

        if x25 != []:   ax.errorbar(x25,   y25,   yerr=y25unc,   label="Disc. %d = 0.25"    %(disc),        xerr=xWidths25,   fmt='', capsize=0, color="red",    lw=0, elinewidth=2, marker="o", markeredgecolor="red",    markerfacecolor="red"   )
        if x50 != []:   ax.errorbar(x50,   y50,   yerr=y50unc,   label="Disc. %d = 0.50"    %(disc),        xerr=xWidths50,   fmt='', capsize=0, color="blue",   lw=0, elinewidth=2, marker="o", markeredgecolor="blue",   markerfacecolor="blue"  )
        if x75 != []:   ax.errorbar(x75,   y75,   yerr=y75unc,   label="Disc. %d = 0.75"    %(disc),        xerr=xWidths75,   fmt='', capsize=0, color="green",  lw=0, elinewidth=2, marker="o", markeredgecolor="green",  markerfacecolor="green" )
        if xDiag != []: ax.errorbar(xDiag, yDiag, yerr=yDiagUnc, label="Disc. %d = Disc. %d"%(disc,3-disc), xerr=xWidthsDiag, fmt='', capsize=0, color="purple", lw=0, elinewidth=2, marker="o", markeredgecolor="purple", markerfacecolor="purple")

        if ylim != -1.0:
             ax.set_ylim((0.0, ylim))

        ax.set_ylabel(ylabel); ax.set_xlabel("Disc. %d Value"%(3-disc))
        plt.legend(loc='best', numpoints=1)

        ax.text(0.12, 1.05, 'CMS',                     transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
        ax.text(0.33, 1.04, 'Preliminary',             transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
        ax.text(0.99, 1.04, '%s (13 TeV)' % self.year, transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right') 

        fig.tight_layout()
        fig.savefig('%s/%s_%s_Slices_Disc%d_Njets%s_%s_%s_%s.pdf' % (self.outputDir, self.year, name, disc, Njets, region, channel, tag), dpi=fig.dpi)

        plt.close(fig)
