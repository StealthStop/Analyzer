import ROOT
import os
import sys
import math
import ctypes
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


class Common_Calculations_Plotters():

    def __init__(self, year, model, mass, channel, metric, Njets=-1):
        self.Njets   = Njets
        self.year    = year
        self.model   = model
        self.mass    = mass
        self.channel = channel
        self.metric  = metric

    def __del__(self):
        del self.Njets
        del self.year
        del self.model
        del self.mass
        del self.channel
        del self.metric 

    # ----------------------
    # calculate all closures
    # ----------------------
    def cal_simpleClosure_ABCD(self, nBkgEvents_A, nBkgEvents_B, nBkgEvents_C, nBkgEvents_D, nBkgEventsErr_A, nBkgEventsErr_B, nBkgEventsErr_C, nBkgEventsErr_D):
       
        #      B    |    A    
        # __________|__________        
        #           |        
        #      D    |    C  
 
        numerator   = nBkgEvents_B * nBkgEvents_C
        denominator = nBkgEvents_A * nBkgEvents_D

        nPredBkgEvents_A = -1.0; nPredBkgUncEvents_A = 0.0

        if nBkgEvents_D > 0.0:
            nPredBkgEvents_A    = ( numerator / nBkgEvents_D )
            nPredBkgUncEvents_A = ( ( (nBkgEvents_C * nBkgEventsErr_B) / nBkgEvents_D )**2.0
                                  + ( (nBkgEvents_B * nBkgEventsErr_C) / nBkgEvents_D )**2.0
                                  + ( (nBkgEvents_B * nBkgEvents_C * nBkgEventsErr_D) / nBkgEvents_D**2.0 )**2.0 )**0.5

        return nPredBkgEvents_A, nPredBkgUncEvents_A

    # -----------------
    # plot all closures
    # -----------------
    def plot_ClosureNjets(self, bkg, bkgUnc, bkgPred, bkgPredUnc, Njets, name=''):
        
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
                pull         = (bkgPred[i] - bkg[i]) / bkgUnc[i]
                closureError = 1.0 - bkgPred[i] / bkg[i]
                abcdPull.append(pull)
                abcdError.append(closureError)
                totalChi2    += pull ** 2.0
                ndof         += 1

        fig = plt.figure(figsize=(5, 5))
        ax  = fig.add_subplot(111)
        fig.subplots_adjust(hspace=0)

        lowerNjets  = Njets[0]
        higherNjets = Njets[(-1)]

        ax.spines['top'].set_color('none')
        ax.spines['bottom'].set_color('none')
        ax.spines['left'].set_color('none')
        ax.spines['right'].set_color('none')
        ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)

        ax1 = fig.add_subplot(3, 1, (1, 2))
        ax1.set_yscale('log')
        ax1.set_xlim([lowerNjets - 0.5, higherNjets + 0.5])
        ax1.text(0.05, 0.1, '$\\chi^2$ / ndof = %3.2f' % (totalChi2 / float(ndof)), horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes, fontsize=10)
        ax1.text(0.05, 0.25, '%s Metric'% (self.metric), horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes, fontsize=14, fontweight='bold')
        ax1.text(0.12, 1.05, 'CMS',                     transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
        ax1.text(0.33, 1.04, 'Preliminary',             transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
        ax1.text(0.99, 1.04, '%s (13 TeV)' % self.year, transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right')
        ax1.set_ylabel('Unweighted Event Counts')
        ax1.errorbar(binCenters, pred, yerr=predUnc, label='Predicted', xerr=xErr, fmt='', color='red', lw=0, elinewidth=2, marker='o', markerfacecolor='red', markersize=4.0)
        ax1.errorbar(binCenters, obs, yerr=unc, label='Observed', xerr=xErr, fmt='', color='black', lw=0, elinewidth=2, marker='o', markerfacecolor='black', markersize=4.0)

        ax2 = fig.add_subplot(3, 1, 3)
        ax2.set_xlim([lowerNjets - 0.5, higherNjets + 0.5])
        ax2.errorbar(binCenters, abcdError, yerr=None, xerr=xErr, fmt='', color='blue', lw=0, elinewidth=2, marker='o', markerfacecolor='blue', markersize=4.0)
        ax2.axhline(y=0.0, color='black', linestyle='dashed', lw=1)
        ax2.grid(axis='y', color='black', linestyle='dashed', which='both')
        ax2.set_xlabel('Number of jets')
        ax2.set_ylabel('1 - Pred./Obs.', fontsize='small')
        ax2.set_ylim([-1.6, 1.6])

        plt.xticks(Njets)

        ax1.legend(loc='best', frameon=False)

        fig.savefig('plots/%s_%s/%s/%s_Njets_Region_A_PredVsActual_%s_%s_%s.pdf' % (self.model, self.mass, self.channel, self.year, name, self.channel, self.metric))

        plt.close(fig)

        return totalChi2, ndof
   
    # ----------------------
    # make all closure plots
    # ----------------------    
    def make_allClosures(self, FinalBinEdges, list_nTotBkgCount_ABCD, name = ""):

        bkgNjets       = {'A': [], 'B': [], 'C': [], 'D': [], 'A2': [], 'B2': [], 'C2': [], 'D2': []}
        bkgNjetsErr    = {'A': [], 'B': [], 'C': [], 'D': [], 'A2': [], 'B2': [], 'C2': [], 'D2': []}
        bkgNjetsPred_A = {'value': [], 'error': [], 'value_val': [], 'error_val': []}
        Njets          = None

        for i in range(0, len(FinalBinEdges)):

            nTotBkgCount_ABCD = list_nTotBkgCount_ABCD[i]

            if FinalBinEdges[i][0] == -1.0 or FinalBinEdges[i][1] == -1.0:

                bkgNjets['A'].append(0.0)
                bkgNjetsErr['A'].append(0.0)
                bkgNjetsPred_A['value'].append(0.0)
                bkgNjetsPred_A['error'].append(0.0)

            else:
                pred_A, predUnc_A = self.cal_simpleClosure_ABCD( nTotBkgCount_ABCD['nBkgEvents_A'   ][FinalBinEdges[i][0]][FinalBinEdges[i][1]], nTotBkgCount_ABCD['nBkgEvents_B'   ][FinalBinEdges[i][0]][FinalBinEdges[i][1]],
                                                                 nTotBkgCount_ABCD['nBkgEvents_C'   ][FinalBinEdges[i][0]][FinalBinEdges[i][1]], nTotBkgCount_ABCD['nBkgEvents_D'   ][FinalBinEdges[i][0]][FinalBinEdges[i][1]],
                                                                 nTotBkgCount_ABCD['nBkgEventsErr_A'][FinalBinEdges[i][0]][FinalBinEdges[i][1]], nTotBkgCount_ABCD['nBkgEventsErr_B'][FinalBinEdges[i][0]][FinalBinEdges[i][1]],
                                                                 nTotBkgCount_ABCD['nBkgEventsErr_C'][FinalBinEdges[i][0]][FinalBinEdges[i][1]], nTotBkgCount_ABCD['nBkgEventsErr_D'][FinalBinEdges[i][0]][FinalBinEdges[i][1]])

                bkgNjets['A'].append(nTotBkgCount_ABCD['nBkgEvents_A'][FinalBinEdges[i][0]][FinalBinEdges[i][1]])
                bkgNjetsErr['A'].append(nTotBkgCount_ABCD['nBkgEventsErr_A'][FinalBinEdges[i][0]][FinalBinEdges[i][1]])
                bkgNjetsPred_A['value'].append(pred_A)
                bkgNjetsPred_A['error'].append(predUnc_A)

        if self.channel == '0l':
            Njets = [
             7, 8, 9, 10, 11, 12]
        else:
            Njets = [
             7, 8, 9, 10, 11]

        self.plot_ClosureNjets(bkgNjets['A'], bkgNjetsErr['A'], bkgNjetsPred_A['value'], bkgNjetsPred_A['error'], Njets, name)       

    # --------------------------------------------
    # plot significance as a function of bin edges 
    # --------------------------------------------
    def plot_Significance_vsDisc1Disc2(self, nBins, inverseSignificance, disc1Edges, disc2Edges, c1, c2, minEdge, maxEdge, binWidth, Njets = -1, name=""):

        # significance as a function of bin edges
        fig = plt.figure()
        plt.hist2d(disc1Edges, disc2Edges, bins=[nBins, nBins], range=[[0.0, 1.0], [0.0, 1.0]], cmap=plt.cm.jet, weights=np.reciprocal(inverseSignificance), cmin=1e-09, cmax=10.0)
        plt.colorbar()
        ax = plt.gca()
        ax.set_xlabel('Disc. 1 Bin Edge')
        ax.set_ylabel('Disc. 2 Bin Edge')
        ax.text(0.12, 1.05, 'CMS',                     transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
        ax.text(0.33, 1.04, 'Preliminary',             transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
        ax.text(0.99, 1.04, '%s (13 TeV)' % self.year, transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right')
        l1 = ml.Line2D([c1, c1], [0.0, 1.0], color='black', linewidth=2, linestyle='dashed')
        l2 = ml.Line2D([0.0, 1.0], [c2, c2], color='black', linewidth=2, linestyle='dashed')
        ax.add_line(l1)
        ax.add_line(l2)

        if Njets == -1:
            fig.savefig('plots/%s_%s/%s/%s_Sign_vs_Disc1Disc2_%s_%s_%s.pdf' % (self.model, self.mass, self.channel, self.year, name, self.channel, self.metric), dpi=fig.dpi)
        else:
            fig.savefig('plots/%s_%s/%s/%s_Sign_vs_Disc1Disc2%s_%s_%s_%s.pdf' % (self.model, self.mass, self.channel, self.year, Njets, name, self.channel, self.metric), dpi=fig.dpi)

        plt.close(fig)

    # -----------------------------------------------------------------
    # plot closure error and closure error unc. a function of bin edges 
    # -----------------------------------------------------------------
    def plot_ClosureError_vsDisc1Disc2(self, nBins, closureErrsList, closureErrUncList, disc1Edges, disc2Edges, c1, c2, minEdge, maxEdge, binWidth, Njets = -1, name=""):

        # closure error as a function of bin edges
        fig = plt.figure()
        plt.hist2d(disc1Edges, disc2Edges, bins=[nBins, nBins], range=[[0.0, 1.0], [0.0, 1.0]], cmap=plt.cm.jet, weights=closureErrsList, cmin=1e-09, cmax=2.5)
        plt.colorbar()
        ax = plt.gca()
        ax.set_xlabel('Disc. 1 Bin Edge')
        ax.set_ylabel('Disc. 2 Bin Edge')
        ax.text(0.12, 1.05, 'CMS',                     transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
        ax.text(0.33, 1.04, 'Preliminary',             transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
        ax.text(0.99, 1.04, '%s (13 TeV)' % self.year, transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right')
        l1 = ml.Line2D([c1, c1], [0.0, 1.0], color='black', linewidth=2, linestyle='dashed')
        l2 = ml.Line2D([0.0, 1.0], [c2, c2], color='black', linewidth=2, linestyle='dashed')
        ax.add_line(l1)
        ax.add_line(l2)

        if Njets == -1:
            fig.savefig('plots/%s_%s/%s/%s_ClosureErr_vs_Disc1Disc2_%s_%s_%s.pdf' % (self.model, self.mass, self.channel, self.year, name, self.channel, self.metric), dpi=fig.dpi)
        else:
            fig.savefig('plots/%s_%s/%s/%s_ClosureErr_vs_Disc1Disc2%s_%s_%s_%s.pdf' % (self.model, self.mass, self.channel, self.year, self.Njets, name, self.channel, self.metric), dpi=fig.dpi)

        plt.close(fig)

        # closure error uncertainty as a function of bin edges
        fig = plt.figure()
        plt.hist2d(disc1Edges, disc2Edges, bins=[nBins, nBins], range=[[0.0, 1.0], [0.0, 1.0]], cmap=plt.cm.jet, weights=closureErrUncList, cmin=1e-09, cmax=2.5)
        plt.colorbar()
        ax = plt.gca()
        ax.set_xlabel('Disc. 1 Bin Edge')
        ax.set_ylabel('Disc. 2 Bin Edge')
        ax.text(0.12, 1.05, 'CMS',                     transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
        ax.text(0.33, 1.04, 'Preliminary',             transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
        ax.text(0.99, 1.04, '%s (13 TeV)' % self.year, transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right')
        l1 = ml.Line2D([c1, c1], [0.0, 1.0], color='black', linewidth=2, linestyle='dashed')
        l2 = ml.Line2D([0.0, 1.0], [c2, c2], color='black', linewidth=2, linestyle='dashed')
        ax.add_line(l1)
        ax.add_line(l2)

        if Njets == -1:
            fig.savefig('plots/%s_%s/%s/%s_ClosureErrUnc_vs_Disc1Disc2_%s_%s_%s.pdf' % (self.model, self.mass, self.channel, self.year, name, self.channel, self.metric), dpi=fig.dpi)
        else:
            fig.savefig('plots/%s_%s/%s/%s_ClosureErrUnc_vs_Disc1Disc2%s_%s_%s_%s.pdf' % (self.model, self.mass, self.channel, self.year, Njets, name, self.channel, self.metric), dpi=fig.dpi)

        plt.close(fig)

    # --------------------------------------
    # plot inverseSignificance vs ClosureErr
    # -------------------------------------- 
    def plot_inverseSignificance_vsClosureErr(self, significance, closureErr, inverseSignificance, closureErrsList, edges, disc1Edge, disc2Edge, Njets = -1, name=""):

        fig = plt.figure(figsize=(5, 5))
        ax = plt.gca()
        ax.text(0.12, 1.05, 'CMS',                     transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
        ax.text(0.33, 1.04, 'Preliminary',             transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
        ax.text(0.99, 1.04, '%s (13 TeV)' % self.year, transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right')
        plt.scatter(inverseSignificance, closureErrsList, color='xkcd:grey', marker='o', label='1 - Pred./Obs. vs 1 / Significance')

        if significance != 0.0:
            plt.scatter([1.0 / significance], [closureErr], color='xkcd:red', marker='o', label='Chosen Solution')
        plt.xlabel('1 / Significance')
        plt.xlim(left=0)
        plt.ylabel('|1 - Pred./Obs.|')
        plt.legend(loc='best', frameon=False)
        plt.ylim(bottom=0)
        plt.gca().invert_yaxis()
        plt.text(0.4, 0.85, '$%.2f < \\bf{Disc.\\;1\\;Edge}$ = %s < %.2f' % (edges[0], disc1Edge, edges[(-1)]), transform=ax.transAxes, fontsize=8)
        plt.text(0.4, 0.8, '$%.2f < \\bf{Disc.\\;2\\;Edge}$ = %s < %.2f' % (edges[0], disc2Edge, edges[(-1)]), transform=ax.transAxes, fontsize=8)

        if Njets == -1:
            fig.savefig('plots/%s_%s/%s/%s_InvSign_vs_ClosureErr_%s_%s_%s.pdf' % (self.model, self.mass, self.channel, self.year, name, self.channel, self.metric), dpi=fig.dpi)
        else:
            fig.savefig('plots/%s_%s/%s/%s_InvSign_vs_ClosureErr%s_%s_%s_%s.pdf' % (self.model, self.mass, self.channel, self.year, Njets, name, self.channel, self.metric), dpi=fig.dpi)

        plt.close(fig)

    # ----------------------------------
    # plot Disc1 vs Disc2 in each Region
    # ----------------------------------
    def plot_Disc1VsDisc2(self, hist, disc1Edge_ABCD, disc2Edge_ABCD, disc1Edge_bdEF, disc2Edge_cdiGH, tag = ""):

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
        l1 = ml.Line2D([disc1Edge_ABCD, disc1Edge_ABCD], [0.0, 1.0],   color="darkviolet", linewidth=4, linestyle='solid')
        l2 = ml.Line2D([0.0, 1.0], [disc2Edge_ABCD, disc2Edge_ABCD],   color="darkviolet", linewidth=4, linestyle='solid')
        l3 = ml.Line2D([disc1Edge_bdEF, disc1Edge_bdEF], [0.0, 1.0],   color="yellow", linewidth=4, linestyle='solid')
        l4 = ml.Line2D([0.0, 1.0], [disc2Edge_cdiGH, disc2Edge_cdiGH], color="lime",   linewidth=4, linestyle='solid')
        ax.add_line(l1)
        ax.add_line(l2)
        ax.add_line(l3)
        ax.add_line(l4)

        if self.Njets == -1:
            fig.savefig('plots/%s_%s/%s/%s_Disc1VsDisc2_%s_%s.pdf' % (self.model, self.mass, self.channel, self.year, tag, self.channel), dpi=fig.dpi)
        else:
            fig.savefig('plots/%s_%s/%s/%s_Disc1VsDisc2_%s%s_%s.pdf' % (self.model, self.mass, self.channel, self.year, tag, self.Njets, self.channel), dpi=fig.dpi)
        plt.close(fig)

    # ---------------------------------------------
    # plot SigFrac vs Disc1Disc2 in each A, B, C, D
    # ---------------------------------------------
    def plot_SigFrac_vsDisc1Disc2(self, nBins, sigFracsA, sigFracsB, sigFracsC, sigFracsD, disc1Edges, disc2Edges, c1, c2, minEdge, maxEdge, binWidth, name=''):

        # SigFracA vs Disc1Disc2
        fig = plt.figure()
        plt.hist2d(disc1Edges, disc2Edges, bins=[nBins, nBins], range=[[0.0, 1.0], [0.0, 1.0]], cmap=plt.cm.jet, weights=sigFracsA, cmin=1e-05, cmax=1.0)
        plt.colorbar()
        ax = plt.gca()
        ax.set_xlabel('Disc. 1 Bin Edge')
        ax.set_ylabel('Disc. 2 Bin Edge')
        ax.text(0.12, 1.05, 'CMS',                     transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
        ax.text(0.33, 1.04, 'Preliminary',             transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
        ax.text(0.99, 1.04, '%s (13 TeV)' % self.year, transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right')
        l1 = ml.Line2D([c1, c1], [0.0, 1.0], color='black', linewidth=2, linestyle='dashed')
        l2 = ml.Line2D([0.0, 1.0], [c2, c2], color='black', linewidth=2, linestyle='dashed')
        ax.add_line(l1)
        ax.add_line(l2)
        if self.Njets == -1:
            fig.savefig('plots/%s_%s/%s/%s_SigFracA_vs_Disc1Disc2_%s%s.pdf' % (self.model, self.mass, self.channel, self.year, self.channel, name), dpi=fig.dpi)
        else:
            fig.savefig('plots/%s_%s/%s/%s_SigFracA_vs_Disc1Disc2%s_%s%s.pdf' % (self.model, self.mass, self.channel, self.year, self.Njets, self.channel, name), dpi=fig.dpi)
        plt.close(fig)

        # SigFracB vs Disc1Disc2
        fig = plt.figure()
        plt.hist2d(disc1Edges, disc2Edges, bins=[nBins, nBins], range=[[0.0, 1.0], [0.0, 1.0]], cmap=plt.cm.jet, weights=sigFracsB, cmin=1e-05, cmax=1.0)
        plt.colorbar()
        ax = plt.gca()
        ax.set_xlabel('Disc. 1 Bin Edge')
        ax.set_ylabel('Disc. 2 Bin Edge')
        ax.text(0.12, 1.05, 'CMS',                     transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
        ax.text(0.33, 1.04, 'Preliminary',             transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
        ax.text(0.99, 1.04, '%s (13 TeV)' % self.year, transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right')
        l1 = ml.Line2D([c1, c1], [0.0, 1.0], color='black', linewidth=2, linestyle='dashed')
        l2 = ml.Line2D([0.0, 1.0], [c2, c2], color='black', linewidth=2, linestyle='dashed')
        ax.add_line(l1)
        ax.add_line(l2)

        if self.Njets == -1:
            fig.savefig('plots/%s_%s/%s/%s_SigFracB_vs_Disc1Disc2_%s%s.pdf' % (self.model, self.mass, self.channel, self.year, self.channel, name), dpi=fig.dpi)
        else:
            fig.savefig('plots/%s_%s/%s/%s_SigFracB_vs_Disc1Disc2%s_%s%s.pdf' % (self.model, self.mass, self.channel, self.year, self.Njets, self.channel, name), dpi=fig.dpi)
        plt.close(fig)
 
        # SigFracC vs Disc1Disc2
        fig = plt.figure()
        plt.hist2d(disc1Edges, disc2Edges, bins=[nBins, nBins], range=[[0.0, 1.0], [0.0, 1.0]], cmap=plt.cm.jet, weights=sigFracsC, cmin=1e-05, cmax=1.0)
        plt.colorbar()
        ax = plt.gca()
        ax.set_xlabel('Disc. 1 Bin Edge')
        ax.set_ylabel('Disc. 2 Bin Edge')
        ax.text(0.12, 1.05, 'CMS',                     transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
        ax.text(0.33, 1.04, 'Preliminary',             transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
        ax.text(0.99, 1.04, '%s (13 TeV)' % self.year, transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right')
        l1 = ml.Line2D([c1, c1], [0.0, 1.0], color='black', linewidth=2, linestyle='dashed')
        l2 = ml.Line2D([0.0, 1.0], [c2, c2], color='black', linewidth=2, linestyle='dashed')
        ax.add_line(l1)
        ax.add_line(l2)

        if self.Njets == -1:
            fig.savefig('plots/%s_%s/%s/%s_SigFracC_vs_Disc1Disc2_%s%s.pdf' % (self.model, self.mass, self.channel, self.year, self.channel, name), dpi=fig.dpi)
        else:
            fig.savefig('plots/%s_%s/%s/%s_SigFracC_vs_Disc1Disc2%s_%s%s.pdf' % (self.model, self.mass, self.channel, self.year, self.Njets, self.channel, name), dpi=fig.dpi)
        plt.close(fig)

        # SigFracD vs Disc1Disc2
        fig = plt.figure()
        plt.hist2d(disc1Edges, disc2Edges, bins=[nBins, nBins], range=[[0.0, 1.0], [0.0, 1.0]], cmap=plt.cm.jet, weights=sigFracsD, cmin=1e-05, cmax=1.0)
        plt.colorbar()
        ax = plt.gca()
        ax.set_xlabel('Disc. 1 Bin Edge')
        ax.set_ylabel('Disc. 2 Bin Edge')
        ax.text(0.12, 1.05, 'CMS',                     transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
        ax.text(0.33, 1.04, 'Preliminary',             transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
        ax.text(0.99, 1.04, '%s (13 TeV)' % self.year, transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right')
        l1 = ml.Line2D([c1, c1], [0.0, 1.0], color='black', linewidth=2, linestyle='dashed')
        l2 = ml.Line2D([0.0, 1.0], [c2, c2], color='black', linewidth=2, linestyle='dashed')
        ax.add_line(l1)
        ax.add_line(l2)

        if self.Njets == -1:
            fig.savefig('plots/%s_%s/%s/%s_SigFracD_vs_Disc1Disc2_%s%s.pdf' % (self.model, self.mass, self.channel, self.year, self.channel, name), dpi=fig.dpi)
        else:
            fig.savefig('plots/%s_%s/%s/%s_SigFracD_vs_Disc1Disc2%s_%s%s.pdf' % (self.model, self.mass, self.channel, self.year, self.Njets, self.channel, name), dpi=fig.dpi)
        plt.close(fig)

    # -------------------------------------------------
    # plot SigTotFracs vs Disc1Disc2 in each A, B, C, D
    # ------------------------------------------------- 
    def plot_SigTotFrac_vsDisc1Disc2(self, nBins, sigTotFracsA, sigTotFracsB, sigTotFracsC, sigTotFracsD, disc1Edges, disc2Edges, c1, c2, minEdge, maxEdge, binWidth):
       
        # SigTotFracA vs Disc1Disc2
        fig = plt.figure()
        plt.hist2d(disc1Edges, disc2Edges, bins=[nBins, nBins], range=[[0.0, 1.0], [0.0, 1.0]], cmap=plt.cm.jet, weights=sigTotFracsA, cmin=1e-05, cmax=1.0)
        plt.colorbar()
        ax = plt.gca()
        ax.set_xlabel('Disc. 1 Bin Edge')
        ax.set_ylabel('Disc. 2 Bin Edge')
        ax.text(0.12, 1.05, 'CMS',                transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
        ax.text(0.33, 1.04, 'Preliminary',        transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
        ax.text(0.99, 1.04, '%s (13 TeV)' % year, transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right')
        l1 = ml.Line2D([c1, c1], [0.0, 1.0], color='black', linewidth=2, linestyle='dashed')
        l2 = ml.Line2D([0.0, 1.0], [c2, c2], color='black', linewidth=2, linestyle='dashed')
        ax.add_line(l1)
        ax.add_line(l2)
        
        if Njets == -1:
            fig.savefig('plots/%s_%s/%s/%s_SigTotFracA_vs_Disc1Disc2_%s.pdf' % (self.model, self.mass, self.channel, self.year, self.channel), dpi=fig.dpi)
        else:
            fig.savefig('plots/%s_%s/%s/%s_SigTotFracA_vs_Disc1Disc2%s_%s.pdf' % (self.model, self.mass, self.channel, self.year, self.Njets, self.channel), dpi=fig.dpi)
        plt.close(fig)
        
        # SigTotFracB vs Disc1Disc2
        fig = plt.figure()
        plt.hist2d(disc1Edges, disc2Edges, bins=[nBins, nBins], range=[[0.0, 1.0], [0.0, 1.0]], cmap=plt.cm.jet, weights=sigTotFracsB, cmin=1e-05, cmax=1.0)
        plt.colorbar()
        ax = plt.gca()
        ax.set_xlabel('Disc. 1 Bin Edge')
        ax.set_ylabel('Disc. 2 Bin Edge')
        ax.text(0.12, 1.05, 'CMS',                transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
        ax.text(0.33, 1.04, 'Preliminary',        transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
        ax.text(0.99, 1.04, '%s (13 TeV)' % year, transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right')
        l1 = ml.Line2D([c1, c1], [0.0, 1.0], color='black', linewidth=2, linestyle='dashed')
        l2 = ml.Line2D([0.0, 1.0], [c2, c2], color='black', linewidth=2, linestyle='dashed')
        ax.add_line(l1)
        ax.add_line(l2)

        if Njets == -1:
            fig.savefig('plots/%s_%s/%s/%s_SigTotFracB_vs_Disc1Disc2_%s.pdf' % (self.model, self.mass, self.channel, self.year, self.channel), dpi=fig.dpi)
        else:
            fig.savefig('plots/%s_%s/%s/%s_SigTotFracB_vs_Disc1Disc2%s_%s.pdf' % (self.model, self.mass, self.channel, self.year, self.Njets, self.channel), dpi=fig.dpi)
        plt.close(fig)

        # SigTotFracC vs Disc1Disc2
        fig = plt.figure()
        plt.hist2d(disc1Edges, disc2Edges, bins=[nBins, nBins], range=[[0.0, 1.0], [0.0, 1.0]], cmap=plt.cm.jet, weights=sigTotFracsC, cmin=1e-05, cmax=1.0)
        plt.colorbar()
        ax = plt.gca()
        ax.set_xlabel('Disc. 1 Bin Edge')
        ax.set_ylabel('Disc. 2 Bin Edge')
        ax.text(0.12, 1.05, 'CMS',                transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
        ax.text(0.33, 1.04, 'Preliminary',        transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
        ax.text(0.99, 1.04, '%s (13 TeV)' % year, transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right')
        l1 = ml.Line2D([c1, c1], [0.0, 1.0], color='black', linewidth=2, linestyle='dashed')
        l2 = ml.Line2D([0.0, 1.0], [c2, c2], color='black', linewidth=2, linestyle='dashed')
        ax.add_line(l1)
        ax.add_line(l2)

        if Njets == -1:
            fig.savefig('plots/%s_%s/%s/%s_SigTotFracC_vs_Disc1Disc2_%s.pdf' % (self.model, self.mass, self.channel, self.year, self.channel), dpi=fig.dpi)
        else:
            fig.savefig('plots/%s_%s/%s/%s_SigTotFracC_vs_Disc1Disc2%s_%s.pdf' % (self.model, self.mass, self.channel, self.year, self.Njets, self.channel), dpi=fig.dpi)
        plt.close(fig)

        # SigTotFracD vs Disc1Disc2
        fig = plt.figure()
        plt.hist2d(disc1Edges, disc2Edges, bins=[nBins, nBins], range=[[0.0, 1.0], [0.0, 1.0]], cmap=plt.cm.jet, weights=sigTotFracsD, cmin=1e-05, cmax=1.0)
        plt.colorbar()
        ax = plt.gca()
        ax.set_xlabel('Disc. 1 Bin Edge')
        ax.set_ylabel('Disc. 2 Bin Edge')
        ax.text(0.12, 1.05, 'CMS',                transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
        ax.text(0.33, 1.04, 'Preliminary',        transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
        ax.text(0.99, 1.04, '%s (13 TeV)' % year, transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right')
        l1 = ml.Line2D([c1, c1], [0.0, 1.0], color='black', linewidth=2, linestyle='dashed')
        l2 = ml.Line2D([0.0, 1.0], [c2, c2], color='black', linewidth=2, linestyle='dashed')
        ax.add_line(l1)
        ax.add_line(l2)

        if Njets == -1:
            fig.savefig('plots/%s_%s/%s/%s_SigTotFracD_vs_Disc1Disc2_%s.pdf' % (self.model, self.mass, self.channel, self.year, self.channel), dpi=fig.dpi)
        else:
            fig.savefig('plots/%s_%s/%s/%s_SigTotFracD_vs_Disc1Disc2%s_%s.pdf' % (self.model, self.mass, self.channel, self.year, self.Njets, self.channel), dpi=fig.dpi)
        plt.close(fig)

    # -------------------------------------------------
    # plot BkgTotFracs vs Disc1Disc2 in each A, B, C, D
    # -------------------------------------------------
    def plot_BkgTotFrac_vsDisc1Disc2(self, nBins, bkgTotFracsA, bkgTotFracsB, bkgTotFracsC, bkgTotFracsD, disc1Edges, disc2Edges, c1, c2, minEdge, maxEdge, binWidth):
        
        # BkgTotFracA vs Disc1Disc2
        fig = plt.figure()
        plt.hist2d(disc1Edges, disc2Edges, bins=[nBins, nBins], range=[[0.0, 1.0], [0.0, 1.0]], cmap=plt.cm.jet, weights=bkgTotFracsA, cmin=1e-05, cmax=1.0)
        plt.colorbar()
        ax = plt.gca()
        ax.set_xlabel('Disc. 1 Bin Edge')
        ax.set_ylabel('Disc. 2 Bin Edge')
        ax.text(0.12, 1.05, 'CMS',                transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
        ax.text(0.33, 1.04, 'Preliminary',        transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
        ax.text(0.99, 1.04, '%s (13 TeV)' % year, transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right')
        l1 = ml.Line2D([c1, c1], [0.0, 1.0], color='black', linewidth=2, linestyle='dashed')
        l2 = ml.Line2D([0.0, 1.0], [c2, c2], color='black', linewidth=2, linestyle='dashed')
        ax.add_line(l1)
        ax.add_line(l2)

        if Njets == -1:
            fig.savefig('plots/%s_%s/%s/%s_BkgTotFracA_vs_Disc1Disc2_%s.pdf' % (self.model, self.mass, self.channel, self.year, self.channel), dpi=fig.dpi)
        else:
            fig.savefig('plots/%s_%s/%s/%s_BkgTotFracA_vs_Disc1Disc2%s_%s.pdf' % (self.model, self.mass, self.channel, self.year, self.Njets, self.channel), dpi=fig.dpi)
        plt.close(fig)

        # BkgTotFracB vs Disc1Disc2
        fig = plt.figure()
        plt.hist2d(disc1Edges, disc2Edges, bins=[nBins, nBins], range=[[0.0, 1.0], [0.0, 1.0]], cmap=plt.cm.jet, weights=bkgTotFracsB, cmin=1e-05, cmax=1.0)
        plt.colorbar()
        ax = plt.gca()
        ax.set_xlabel('Disc. 1 Bin Edge')
        ax.set_ylabel('Disc. 2 Bin Edge')
        ax.text(0.12, 1.05, 'CMS',                transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
        ax.text(0.33, 1.04, 'Preliminary',        transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
        ax.text(0.99, 1.04, '%s (13 TeV)' % year, transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right')
        l1 = ml.Line2D([c1, c1], [0.0, 1.0], color='black', linewidth=2, linestyle='dashed')
        l2 = ml.Line2D([0.0, 1.0], [c2, c2], color='black', linewidth=2, linestyle='dashed')
        ax.add_line(l1)
        ax.add_line(l2)

        if Njets == -1:
            fig.savefig('plots/%s_%s/%s/%s_BkgTotFracB_vs_Disc1Disc2_%s.pdf' % (self.model, self.mass, self.channel, self.year, self.channel), dpi=fig.dpi)
        else:
            fig.savefig('plots/%s_%s/%s/%s_BkgTotFracB_vs_Disc1Disc2%s_%s.pdf' % (self.model, self.mass, self.channel, self.year, self.Njets, self.channel), dpi=fig.dpi)
        plt.close(fig)

        # BkgTotFracC vs Disc1Disc2
        fig = plt.figure()
        plt.hist2d(disc1Edges, disc2Edges, bins=[nBins, nBins], range=[[0.0, 1.0], [0.0, 1.0]], cmap=plt.cm.jet, weights=bkgTotFracsC, cmin=1e-05, cmax=1.0)
        plt.colorbar()
        ax = plt.gca()
        ax.set_xlabel('Disc. 1 Bin Edge')
        ax.set_ylabel('Disc. 2 Bin Edge')
        ax.text(0.12, 1.05, 'CMS',                transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
        ax.text(0.33, 1.04, 'Preliminary',        transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
        ax.text(0.99, 1.04, '%s (13 TeV)' % year, transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right')
        l1 = ml.Line2D([c1, c1], [0.0, 1.0], color='black', linewidth=2, linestyle='dashed')
        l2 = ml.Line2D([0.0, 1.0], [c2, c2], color='black', linewidth=2, linestyle='dashed')
        ax.add_line(l1)
        ax.add_line(l2)

        if Njets == -1:
            fig.savefig('plots/%s_%s/%s/%s_BkgTotFracC_vs_Disc1Disc2_%s.pdf' % (self.model, self.mass, self.year, self.channel), dpi=fig.dpi)
        else:
            fig.savefig('plots/%s_%s/%s/%s_BkgTotFracC_vs_Disc1Disc2%s_%s.pdf' % (self.model, self.mass, self.channel, self.year, self.Njets, self.channel), dpi=fig.dpi)
        plt.close(fig)

        # BkgTotFracD vs Disc1Disc2
        fig = plt.figure()
        plt.hist2d(disc1Edges, disc2Edges, bins=[nBins, nBins], range=[[0.0, 1.0], [0.0, 1.0]], cmap=plt.cm.jet, weights=bkgTotFracsD, cmin=1e-05, cmax=1.0)
        plt.colorbar()
        ax = plt.gca()
        ax.set_xlabel('Disc. 1 Bin Edge')
        ax.set_ylabel('Disc. 2 Bin Edge')
        ax.text(0.12, 1.05, 'CMS',                transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
        ax.text(0.33, 1.04, 'Preliminary',        transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
        ax.text(0.99, 1.04, '%s (13 TeV)' % year, transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right')
        l1 = ml.Line2D([c1, c1], [0.0, 1.0], color='black', linewidth=2, linestyle='dashed')
        l2 = ml.Line2D([0.0, 1.0], [c2, c2], color='black', linewidth=2, linestyle='dashed')
        ax.add_line(l1)
        ax.add_line(l2)
        if Njets == -1:
            fig.savefig('plots/%s_%s/%s/%s_BkgTotFracD_vs_Disc1Disc2_%s.pdf' % (self.model, self.mass, self.year, self.channel), dpi=fig.dpi)
        else:
            fig.savefig('plots/%s_%s/%s/%s_BkgTotFracD_vs_Disc1Disc2%s_%s.pdf' % (self.model, self.mass, self.channel, self.year, self.Njets, self.channel), dpi=fig.dpi)
        plt.close(fig)

    # ------------------------------------------------
    # plot variable vs disc as 1D for ABCD regions
    #   -- Closure vs disc1, disc2 slices
    #   -- Closure vs disc2, disc1 slices
    #   -- Significance vs disc1, disc2 slices
    #   -- Significance vs disc2, disc1 slices
    #   -- weightedEventCounts vs disc1, disc2 slices
    #   -- weightedEventCounts vs disc2, disc1 slices
    # -----------------------------------------------
    def plot_VarVsDisc(self, var, varUncs, d1edges, d2edges, edgeWidth, ylim = -1.0, ylabel = "", tag = "", disc = -1, Njets = -1, name=""):

        x25 = []; x50 = []; x75 = []; xDiag = []
        y25 = []; y50 = []; y75 = []; yDiag = []

        y25unc = []; y50unc = []; y75unc = []; yDiagUnc = []

        edges = (d1edges, d2edges)

        for i in range(0, len(var)):

            if abs(edges[disc-1][i] - 0.25) < 0.01:
               x25.append(edges[2-disc][i])
               y25.append(var[i])
               y25unc.append(varUncs[i])
           
            elif abs(edges[disc-1][i] - 0.50) < 0.01:
               x50.append(edges[2-disc][i])
               y50.append(var[i])
               y50unc.append(varUncs[i])
           
            elif abs(edges[disc-1][i] - 0.75) < 0.01: 
               x75.append(edges[2-disc][i])
               y75.append(var[i])
               y75unc.append(varUncs[i])

            if edges[0][i] == edges[1][i]: 
               xDiag.append(edges[2-disc][i])
               yDiag.append(var[i])
               yDiagUnc.append(varUncs[i])

        fig = plt.figure(figsize=(5, 5))
        ax = plt.gca()
        
        xWidths25   = [edgeWidth for i in range(0, len(x25))]
        xWidths50   = [edgeWidth for i in range(0, len(x50))]
        xWidths75   = [edgeWidth for i in range(0, len(x75))]
        xWidthsDiag = [edgeWidth for i in range(0, len(xDiag))]

        ax.errorbar(x25,   y25,   yerr=y25unc,   label="Disc. %d = 0.25"    %(disc),        xerr=xWidths25,   fmt='', color="red",    lw=0, elinewidth=2, marker="o", markerfacecolor="red"   )
        ax.errorbar(x50,   y50,   yerr=y50unc,   label="Disc. %d = 0.50"    %(disc),        xerr=xWidths50,   fmt='', color="blue",   lw=0, elinewidth=2, marker="o", markerfacecolor="blue"  )
        ax.errorbar(x75,   y75,   yerr=y75unc,   label="Disc. %d = 0.75"    %(disc),        xerr=xWidths75,   fmt='', color="green",  lw=0, elinewidth=2, marker="o", markerfacecolor="green" )
        ax.errorbar(xDiag, yDiag, yerr=yDiagUnc, label="Disc. %d = Disc. %d"%(disc,3-disc), xerr=xWidthsDiag, fmt='', color="purple", lw=0, elinewidth=2, marker="o", markerfacecolor="purple")

        if ylim != -1.0:
             ax.set_ylim((0.0, ylim))

        ax.set_ylabel(ylabel); ax.set_xlabel("Disc. %d Value"%(3-disc))
        plt.legend(loc='best')

        ax.text(0.12, 1.05, 'CMS',                     transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
        ax.text(0.33, 1.04, 'Preliminary',             transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
        ax.text(0.99, 1.04, '%s (13 TeV)' % self.year, transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right') 

        fig.tight_layout()

        if Njets == -1: 
            fig.savefig('plots/%s_%s/%s/%s_%s_Slices_Disc%d_%s_%s.pdf' % (self.model, self.mass, self.year, tag, disc, self.channel, self.metric), dpi=fig.dpi)
        else:        
            fig.savefig('plots/%s_%s/%s/%s_%s_Slices_Disc%d%s_%s_%s_%s.pdf' % (self.model, self.mass, self.channel, self.year, tag, disc, Njets, name, self.channel, self.metric), dpi=fig.dpi)   

        plt.close(fig)

    # --------------------------------------------------
    # plot variable vs disc as 1D for Validation regions
    #   -- Closure vs disc1, disc2 slices
    #   -- Closure vs disc2, disc1 slices
    #   -- weightedEventCounts vs disc1, disc2 slices
    #   -- weightedEventCounts vs disc2, disc1 slices
    # -------------------------------------------------
    def plot_VarVsDisc_Val(self, var, varUncs, valDiscEdges, finalBinEdge, edgeWidth, ylim = -1.0, ylabel = "", tag = "", disc = -1, Njets = -1, name="", col=""):

        fig = plt.figure(figsize=(5, 5))
        ax = plt.gca()

        xWidths = [edgeWidth for i in range(0, len(valDiscEdges))]

        ax.errorbar(valDiscEdges, var,  yerr=varUncs, label="ABCD Disc. %d = %.2f"%(disc, finalBinEdge), xerr=xWidths, fmt='', color=col, lw=0, elinewidth=2, marker="o", markerfacecolor=col)

        if ylim != -1.0:
             ax.set_ylim((0.0, ylim))

        ax.set_ylabel(ylabel); ax.set_xlabel("Disc. %d Value "%(3-disc))
        plt.legend(loc='best')

        ax.text(0.12, 1.05, 'CMS',                     transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
        ax.text(0.33, 1.04, 'Preliminary',             transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
        ax.text(0.99, 1.04, '%s (13 TeV)' % self.year, transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right')

        fig.tight_layout()

        if Njets == -1:
            fig.savefig('plots/%s_%s/%s/%s_%s_Slices_Disc%d_%s_%s.pdf' % (self.model, self.mass, self.year, tag, disc, self.channel, self.metric), dpi=fig.dpi)
        else:
            fig.savefig('plots/%s_%s/%s/%s_%s_Slices_Disc%d%s_%s_%s_%s.pdf' % (self.model, self.mass, self.channel, self.year, tag, disc, self.Njets, name, self.channel, self.metric), dpi=fig.dpi)

        plt.close(fig)



class FinalBinEdges:

    def __init__(self, year, model, mass, channel, metric, Njets = -1):
        self.Njets   = Njets
        self.year    = year
        self.model   = model
        self.mass    = mass
        self.channel = channel
        self.metric  = metric

    def __del__(self):
        del self.Njets
        del self.year
        del self.model
        del self.mass
        del self.channel
        del self.metric

    # ------------------------
    # Significance calculation
    # ------------------------
    def cal_Significance(self, nSigEvents, nBkgEvents, sys=0.3):
        if (nBkgEvents == 0.0):
            return 0
    
        significance = ( nSigEvents / ( nBkgEvents + (sys * nBkgEvents)**2.0 )**0.5 )**2.0
        return significance
    
    # -------------------------
    # Closure error calculation
    # -------------------------
    def cal_ClosureError(self, nBkgEvents_A, nBkgEvents_B, nBkgEvents_C, nBkgEvents_D, nBkgEventsErr_A, nBkgEventsErr_B, nBkgEventsErr_C, nBkgEventsErr_D):
        
        closureError = abs(1.0 - ( (nBkgEvents_B * nBkgEvents_C) / (nBkgEvents_A * nBkgEvents_D) ) )
        
        closureErrUnc = ( ( ( nBkgEvents_C * nBkgEventsErr_B ) / ( nBkgEvents_A * nBkgEvents_D) )**2.0 
                        + ( ( nBkgEvents_B * nBkgEventsErr_C ) / ( nBkgEvents_A * nBkgEvents_D) )**2.0 
                        + ( ( nBkgEvents_B * nBkgEvents_C * nBkgEventsErr_A ) / ( nBkgEvents_A**2.0 * nBkgEvents_D ) )**2.0 
                        + ( ( nBkgEvents_B * nBkgEvents_C * nBkgEventsErr_D ) / ( nBkgEvents_A * nBkgEvents_D**2.0 ) )**2.0 )**0.5

        return closureError, closureErrUnc
   
    # -------------------------------------------
    # Calculate optimization metric of bin edges
    #   This function use the command line option
    #   -- NN optimization metric
    #   -- New optimization metric
    # ------------------------------------------- 
    def cal_OptMetric_ofBinEdges(self, significance, closureError):

        optimizationMetric = None        

        # NN optimization metric
        if self.metric == "NN":
            optimizationMetric  = (closureError)**2 + (1.0 / significance)**2
        
        # New optimization metric
        else: 
            optimizationMetric  = (5 * closureError)**2 + (1.0 / significance)**2

        return optimizationMetric

    # -------------------------------------------------------
    # get signal and background histograms' counts
    #   -- count both signal and background events separately 
    #   -- in each A, B, C, D regions
    #   -- put them to the dictionaries
    # -------------------------------------------------------
    def count_Events_inBinEdges(self, histBkg, histSig):
    
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
            xBinKey     = "%.3f"%xLowBinEdge
    
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
                yBinKey     = "%.3f"%yLowBinEdge
    
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
                nSigEvents_B = ( histSig.IntegralAndError(1,    xBin-1,   yBin, lastYBin, nSigEventsErr_B) )
                nSigEvents_C = ( histSig.IntegralAndError(xBin, lastXBin, 1,    yBin-1,   nSigEventsErr_C) )
                nSigEvents_D = ( histSig.IntegralAndError(1,    xBin-1,   1,    yBin-1,   nSigEventsErr_D) )
                nBkgEvents_A = ( histBkg.IntegralAndError(xBin, lastXBin, yBin, lastYBin, nBkgEventsErr_A) )
                nBkgEvents_B = ( histBkg.IntegralAndError(1,    xBin-1,   yBin, lastYBin, nBkgEventsErr_B) )
                nBkgEvents_C = ( histBkg.IntegralAndError(xBin, lastXBin, 1,    yBin-1,   nBkgEventsErr_C) )
                nBkgEvents_D = ( histBkg.IntegralAndError(1,    xBin-1,   1,    yBin-1,   nBkgEventsErr_D) )
                
                nTotSigCount_ABCD["nSigEvents_A"][xBinKey][yBinKey] = nSigEvents_A; nTotSigCount_ABCD["nSigEventsErr_A"][xBinKey][yBinKey] = nSigEventsErr_A 
                nTotSigCount_ABCD["nSigEvents_B"][xBinKey][yBinKey] = nSigEvents_B; nTotSigCount_ABCD["nSigEventsErr_B"][xBinKey][yBinKey] = nSigEventsErr_B
                nTotSigCount_ABCD["nSigEvents_C"][xBinKey][yBinKey] = nSigEvents_C; nTotSigCount_ABCD["nSigEventsErr_C"][xBinKey][yBinKey] = nSigEventsErr_C
                nTotSigCount_ABCD["nSigEvents_D"][xBinKey][yBinKey] = nSigEvents_D; nTotSigCount_ABCD["nSigEventsErr_D"][xBinKey][yBinKey] = nSigEventsErr_D
    
                nTotBkgCount_ABCD["nBkgEvents_A"][xBinKey][yBinKey] = nBkgEvents_A; nTotBkgCount_ABCD["nBkgEventsErr_A"][xBinKey][yBinKey] = nBkgEventsErr_A 
                nTotBkgCount_ABCD["nBkgEvents_B"][xBinKey][yBinKey] = nBkgEvents_B; nTotBkgCount_ABCD["nBkgEventsErr_B"][xBinKey][yBinKey] = nBkgEventsErr_B
                nTotBkgCount_ABCD["nBkgEvents_C"][xBinKey][yBinKey] = nBkgEvents_C; nTotBkgCount_ABCD["nBkgEventsErr_C"][xBinKey][yBinKey] = nBkgEventsErr_C
                nTotBkgCount_ABCD["nBkgEvents_D"][xBinKey][yBinKey] = nBkgEvents_D; nTotBkgCount_ABCD["nBkgEventsErr_D"][xBinKey][yBinKey] = nBkgEventsErr_D
    
        return nTotSigCount_ABCD, nTotBkgCount_ABCD

    # ----------------------------------------------------------------------------
    # Region by region signal fraction calculation
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
    def get_FinalBinEdges(self, nTotSigCount_ABCD, nTotBkgCount_ABCD, bkgNormUnc = 0.3, minBkgFrac = 0.01, minSigFrac = 0.1):
      
        significance  = 0.0; finalDisc1Key = -1.0; finalDisc2Key = -1.0; closureErr    = 0.0; optMetric  = 999.0
        finalSigFracA = 0.0; finalSigFracB = 0.0;  finalSigFracC = 0.0;  finalSigFracD = 0.0
        nEvents_AB    = 0.0; nEvents_AC    = 0.0 
        final_nBkgEvents_A = 0.0; final_nBkgEvents_B = 0.0; final_nBkgEvents_C = 0.0; final_nBkgEvents_D = 0.0;
        
        inverseSignificance = []; closureErrsList = []; closureErrUncList = []; disc1KeyOut  = []; disc2KeyOut  = []
        sigFracsA           = []; sigFracsB       = []; sigFracsC    = []; sigFracsD    = []
        sigTotFracsA        = []; sigTotFracsB    = []; sigTotFracsC = []; sigTotFracsD = []
        bkgTotFracsA        = []; bkgTotFracsB    = []; bkgTotFracsC = []; bkgTotFracsD = []
  
        sigUncs = [] 
        weighted_Sig_A   = []; weighted_Bkg_A   = []; weighted_SigUnc_A   = []; weighted_BkgUnc_A   = []

        # loop over the disc1 and disc2 to get any possible combination of them
        for disc1Key, disc2s in nTotBkgCount_ABCD["nBkgEvents_A"].items():
            
            for disc2Key, nEvents in disc2s.items():
    
                # number of signal and background events in aech A, B, C, D region
                nSigEvents_A = nTotSigCount_ABCD["nSigEvents_A"][disc1Key][disc2Key]; nBkgEvents_A = nTotBkgCount_ABCD["nBkgEvents_A"][disc1Key][disc2Key] 
                nSigEvents_B = nTotSigCount_ABCD["nSigEvents_B"][disc1Key][disc2Key]; nBkgEvents_B = nTotBkgCount_ABCD["nBkgEvents_B"][disc1Key][disc2Key]
                nSigEvents_C = nTotSigCount_ABCD["nSigEvents_C"][disc1Key][disc2Key]; nBkgEvents_C = nTotBkgCount_ABCD["nBkgEvents_C"][disc1Key][disc2Key]
                nSigEvents_D = nTotSigCount_ABCD["nSigEvents_D"][disc1Key][disc2Key]; nBkgEvents_D = nTotBkgCount_ABCD["nBkgEvents_D"][disc1Key][disc2Key]
            
                nSigEventsErr_A = nTotSigCount_ABCD["nSigEventsErr_A"][disc1Key][disc2Key]; nBkgEventsErr_A = nTotBkgCount_ABCD["nBkgEventsErr_A"][disc1Key][disc2Key]
                nSigEventsErr_B = nTotSigCount_ABCD["nSigEventsErr_B"][disc1Key][disc2Key]; nBkgEventsErr_B = nTotBkgCount_ABCD["nBkgEventsErr_B"][disc1Key][disc2Key]
                nSigEventsErr_C = nTotSigCount_ABCD["nSigEventsErr_C"][disc1Key][disc2Key]; nBkgEventsErr_C = nTotBkgCount_ABCD["nBkgEventsErr_C"][disc1Key][disc2Key]
                nSigEventsErr_D = nTotSigCount_ABCD["nSigEventsErr_D"][disc1Key][disc2Key]; nBkgEventsErr_D = nTotBkgCount_ABCD["nBkgEventsErr_D"][disc1Key][disc2Key]   
 
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
                tempSignificance = 0.0; tempClosureErr = -999.0; tempClosureErrUnc = -999.0; tempOptMetric = 999.0; tempSigUnc = 0.0
    
                if nBkgEvents_A > 0.0:
                    tempSignificance += ( nSigEvents_A / ( nBkgEvents_A + (0.3 * nBkgEvents_A)**2.0 )**0.5 )**2.0 
                if nBkgEvents_B > 0.0:
                    tempSignificance += ( nSigEvents_B / ( nBkgEvents_B + (0.3 * nBkgEvents_B)**2.0 )**0.5 )**2.0
                if nBkgEvents_C > 0.0:
                    tempSignificance += ( nSigEvents_C / ( nBkgEvents_C + (0.3 * nBkgEvents_C)**2.0 )**0.5 )**2.0
                if nBkgEvents_D > 0.0:
                    tempSignificance += ( nSigEvents_D / ( nBkgEvents_D + (0.3 * nBkgEvents_D)**2.0 )**0.5 )**2.0
                
                if nBkgEvents_A > 0.0 and nBkgEvents_D > 0.0: 
                    tempClosureErr, tempClosureErrUnc = self.cal_ClosureError(nBkgEvents_A, nBkgEvents_B, nBkgEvents_C, nBkgEvents_D, nBkgEventsErr_A, nBkgEventsErr_B, nBkgEventsErr_C, nBkgEventsErr_D)
    
                tempSignificance = tempSignificance**0.5
   
                # get the sigUncs to plot variable vs disc as 1D 
                if nBkgEvents_A > 0.0:
                    tempSigUnc += (   ( nSigEventsErr_A / (nBkgEvents_A + (bkgNormUnc * nBkgEvents_A)**2.0 + (tempClosureErr * nBkgEvents_A)**2.0) **0.5 )**2.0
                                  + ( ( nSigEvents_A * nBkgEventsErr_A * (2.0 * nBkgEvents_A * tempClosureErr**2.0 + 2.0 * bkgNormUnc**2.0 * nBkgEvents_A + 1) ) / ( nBkgEvents_A + (bkgNormUnc * nBkgEvents_A)**2.0 + (tempClosureErr * nBkgEvents_A)**2.0 )**1.5 )**2.0
                                  + ( ( nBkgEvents_A**2.0 * tempClosureErr * nSigEvents_A * tempClosureErrUnc) / ( nBkgEvents_A * ( nBkgEvents_A * (tempClosureErr**2.0 + bkgNormUnc**2.0) + 1) )**1.5 )**2.0 )**0.5
 
                if tempSignificance > 0.0 and tempClosureErr > 0.0:
                    inverseSignificance.append(1.0 / tempSignificance) 
                    closureErrsList.append(abs(tempClosureErr))
                    closureErrUncList.append(tempClosureErrUnc)
                    sigUncs.append(tempSigUnc)
                    disc1KeyOut.append(float(disc1Key))
                    disc2KeyOut.append(float(disc2Key)) 

                    # get the weighted signal and backgoround events to plot variable vs disc as 1D 
                    weighted_Sig_A.append(nSigEvents_A)
                    weighted_Bkg_A.append(nBkgEvents_A)   
                    weighted_SigUnc_A.append(nSigEventsErr_A)    
                    weighted_BkgUnc_A.append(nBkgEventsErr_A)  

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
                    tempOptMetric = self.cal_OptMetric_ofBinEdges(tempSignificance, tempClosureErr)
    
                if tempOptMetric < optMetric:
                    finalDisc1Key      = disc1Key 
                    finalDisc2Key      = disc2Key
                    significance       = tempSignificance
                    closureErr         = tempClosureErr
                    optMetric          = tempOptMetric
                    finalSigFracA      = tempSigFracsA
                    finalSigFracB      = tempSigFracsB
                    finalSigFracC      = tempSigFracsC
                    finalSigFracD      = tempSigFracsD
                    nEvents_AB         = nBkgEvents_A + nBkgEvents_B
                    nEvents_AC         = nBkgEvents_A + nBkgEvents_C
                    final_nBkgEvents_A = nBkgEvents_A
                    final_nBkgEvents_B = nBkgEvents_B
                    final_nBkgEvents_C = nBkgEvents_C
                    final_nBkgEvents_D = nBkgEvents_D 

        return finalDisc1Key, finalDisc2Key, significance, closureErr, inverseSignificance, closureErrsList, closureErrUncList, sigUncs, disc1KeyOut, disc2KeyOut, weighted_Sig_A, weighted_Bkg_A, weighted_SigUnc_A, weighted_BkgUnc_A, sigFracsA, sigFracsB, sigFracsC, sigFracsD, sigTotFracsA, sigTotFracsB, sigTotFracsC, sigTotFracsD, bkgTotFracsA, bkgTotFracsB, bkgTotFracsC, bkgTotFracsD, finalSigFracA, finalSigFracB, finalSigFracC, finalSigFracD, nEvents_AB, nEvents_AC, final_nBkgEvents_A, final_nBkgEvents_B, final_nBkgEvents_C, final_nBkgEvents_D 



class ValidationRegions(FinalBinEdges):

    def __init__(self, year, model, mass, channel, Njets = -1):
        self.Njets   = Njets
        self.year    = year
        self.model   = model
        self.mass    = mass
        self.channel = channel

    def __del__(self):
        del self.Njets
        del self.year
        del self.model
        del self.mass
        del self.channel

    # --------------------------------------
    # Calculate metric of validation regions
    # 1. Validation Regions in BD: bdEF
    # 2. Validation Regions in CD: cdiGH
    # --------------------------------------
    def cal_MetricOfValidationRegion_bdEF(self, nBkgEvents_A, nBkgEvents_C, nBkgEvents_b, nBkgEvents_d):

        #   Validation Regions in BD:             
        #          *                              
        #      E   *    B'  |    A                
        #   _______*________|________             
        #          *        |                     
        #      F   *    D'  |    C                
        #          *                              
        
        validationMetric = 999.0

        if (nBkgEvents_b + nBkgEvents_d) > 0.0:
            validationMetric = abs(1.0 - ( (nBkgEvents_A + nBkgEvents_C) / (nBkgEvents_b + nBkgEvents_d) ) )

        return validationMetric

    def cal_MetricOfValidationRegion_cdiGH(self, nBkgEvents_A, nBkgEvents_B, nBkgEvents_c, nBkgEvents_di):

        #   Validation Regions in CD: 
        #
        #          B  |    A   
        #       ______|________
        #             |
        #          D' |    C' 
        #       ***************
        #             |
        #          H  |    G

        validationMetric = 999.0
    
        if (nBkgEvents_c + nBkgEvents_di) > 0.0:
            validationMetric = abs(1.0 - ( (nBkgEvents_A + nBkgEvents_B) / (nBkgEvents_c + nBkgEvents_di) ) )

        return validationMetric

    # --------------------------------------------------
    # get the number of events in each Validation Region
    #   there are two validation regions:
    #   -- Validation region in BD : bdEF
    #   -- Validation region in CD : cdGH
    # --------------------------------------------------
    def count_Events_inValidationRegions(self, histBkg, histSig, finalDisc1Edge, finalDisc2Edge):
    
        lastXBin   = histBkg.GetNbinsX()
        nXBins     = range(2, lastXBin+1) 
        xBin_final = histBkg.GetXaxis().FindBin(float(finalDisc1Edge))

        lastYBin   = histBkg.GetNbinsY()
        nYBins     = range(2, lastYBin+1)
        yBin_final = histBkg.GetYaxis().FindBin(float(finalDisc2Edge))
        
        # ------------------------------
        # Validation region in BD : bdEF
        # ------------------------------
        nTotSigCount_bdEF = {
            "nSigEvents_A" : {},    "nSigEvents_B" : {},    "nSigEvents_C" : {},    "nSigEvents_D" : {}, 
            "nSigEventsErr_A" : {}, "nSigEventsErr_B" : {}, "nSigEventsErr_C" : {}, "nSigEventsErr_D" : {}
        }
    
        nTotBkgCount_bdEF = {
            "nBkgEvents_A" : {},    "nBkgEvents_B" : {},    "nBkgEvents_C" : {},    "nBkgEvents_D" : {}, 
            "nBkgEventsErr_A" : {}, "nBkgEventsErr_B" : {}, "nBkgEventsErr_C" : {}, "nBkgEventsErr_D" : {} 
        }
    
        for xBin in nXBins:
    
            xLowBinEdge = histBkg.GetXaxis().GetBinCenter(xBin)
    
            if (xLowBinEdge > float(finalDisc1Edge)): break
    
            xBinKey = "%.3f"%xLowBinEdge
    
            if xBinKey not in nTotSigCount_bdEF["nSigEventsErr_A"]:
                nTotSigCount_bdEF["nSigEvents_A"][xBinKey] = {}; nTotSigCount_bdEF["nSigEventsErr_A"][xBinKey] = {}
                nTotSigCount_bdEF["nSigEvents_B"][xBinKey] = {}; nTotSigCount_bdEF["nSigEventsErr_B"][xBinKey] = {}
                nTotSigCount_bdEF["nSigEvents_C"][xBinKey] = {}; nTotSigCount_bdEF["nSigEventsErr_C"][xBinKey] = {}
                nTotSigCount_bdEF["nSigEvents_D"][xBinKey] = {}; nTotSigCount_bdEF["nSigEventsErr_D"][xBinKey] = {}
    
            if xBinKey not in nTotBkgCount_bdEF["nBkgEventsErr_A"]:
                nTotBkgCount_bdEF["nBkgEvents_A"][xBinKey] = {}; nTotBkgCount_bdEF["nBkgEventsErr_A"][xBinKey] = {}
                nTotBkgCount_bdEF["nBkgEvents_B"][xBinKey] = {}; nTotBkgCount_bdEF["nBkgEventsErr_B"][xBinKey] = {}
                nTotBkgCount_bdEF["nBkgEvents_C"][xBinKey] = {}; nTotBkgCount_bdEF["nBkgEventsErr_C"][xBinKey] = {}
                nTotBkgCount_bdEF["nBkgEvents_D"][xBinKey] = {}; nTotBkgCount_bdEF["nBkgEventsErr_D"][xBinKey] = {}
    
            # count signal and  background events and errors in validation regions (b, d, E, F)
            nSigEventsErr_b = ROOT.Double(0.0); nSigEventsErr_d = ROOT.Double(0.0); nSigEventsErr_E = ROOT.Double(0.0); nSigEventsErr_F = ROOT.Double(0.0)
            nBkgEventsErr_b = ROOT.Double(0.0); nBkgEventsErr_d = ROOT.Double(0.0); nBkgEventsErr_E = ROOT.Double(0.0); nBkgEventsErr_F = ROOT.Double(0.0)
    
            nSigEvents_b = ( histSig.IntegralAndError(xBin, xBin_final, yBin_final, lastYBin,     nSigEventsErr_b) )
            nSigEvents_E = ( histSig.IntegralAndError(1,    xBin-1,     yBin_final, lastYBin,     nSigEventsErr_E) )
            nSigEvents_d = ( histSig.IntegralAndError(xBin, xBin_final, 1,          yBin_final-1, nSigEventsErr_d) )
            nSigEvents_F = ( histSig.IntegralAndError(1,    xBin-1,     1,          yBin_final-1, nSigEventsErr_F) )
            nBkgEvents_b = ( histBkg.IntegralAndError(xBin, xBin_final, yBin_final, lastYBin,     nBkgEventsErr_b) )
            nBkgEvents_E = ( histBkg.IntegralAndError(1,    xBin-1,     yBin_final, lastYBin,     nBkgEventsErr_E) )
            nBkgEvents_d = ( histBkg.IntegralAndError(xBin, xBin_final, 1,          yBin_final-1, nBkgEventsErr_d) )
            nBkgEvents_F = ( histBkg.IntegralAndError(1,    xBin-1,     1,          yBin_final-1, nBkgEventsErr_F) )
    
            nTotSigCount_bdEF["nSigEvents_A"][xBinKey][finalDisc2Edge] = nSigEvents_b; nTotSigCount_bdEF["nSigEventsErr_A"][xBinKey][finalDisc2Edge] = nSigEventsErr_b
            nTotSigCount_bdEF["nSigEvents_B"][xBinKey][finalDisc2Edge] = nSigEvents_d; nTotSigCount_bdEF["nSigEventsErr_B"][xBinKey][finalDisc2Edge] = nSigEventsErr_d
            nTotSigCount_bdEF["nSigEvents_C"][xBinKey][finalDisc2Edge] = nSigEvents_E; nTotSigCount_bdEF["nSigEventsErr_C"][xBinKey][finalDisc2Edge] = nSigEventsErr_E
            nTotSigCount_bdEF["nSigEvents_D"][xBinKey][finalDisc2Edge] = nSigEvents_F; nTotSigCount_bdEF["nSigEventsErr_D"][xBinKey][finalDisc2Edge] = nSigEventsErr_F
    
            nTotBkgCount_bdEF["nBkgEvents_A"][xBinKey][finalDisc2Edge] = nBkgEvents_b; nTotBkgCount_bdEF["nBkgEventsErr_A"][xBinKey][finalDisc2Edge] = nBkgEventsErr_b
            nTotBkgCount_bdEF["nBkgEvents_B"][xBinKey][finalDisc2Edge] = nBkgEvents_d; nTotBkgCount_bdEF["nBkgEventsErr_B"][xBinKey][finalDisc2Edge] = nBkgEventsErr_d
            nTotBkgCount_bdEF["nBkgEvents_C"][xBinKey][finalDisc2Edge] = nBkgEvents_E; nTotBkgCount_bdEF["nBkgEventsErr_C"][xBinKey][finalDisc2Edge] = nBkgEventsErr_E
            nTotBkgCount_bdEF["nBkgEvents_D"][xBinKey][finalDisc2Edge] = nBkgEvents_F; nTotBkgCount_bdEF["nBkgEventsErr_D"][xBinKey][finalDisc2Edge] = nBkgEventsErr_F

        # -------------------------------
        # Validation region in CD : cdiGH
        # -------------------------------
        nTotSigCount_cdiGH = {
            "nSigEvents_A" : {},    "nSigEvents_B" : {},    "nSigEvents_C" : {},    "nSigEvents_D" : {}, 
            "nSigEventsErr_A" : {}, "nSigEventsErr_B" : {}, "nSigEventsErr_C" : {}, "nSigEventsErr_D" : {}
        }

        nTotBkgCount_cdiGH = {
            "nBkgEvents_A" : {},    "nBkgEvents_B" : {},    "nBkgEvents_C" : {},    "nBkgEvents_D" : {}, 
            "nBkgEventsErr_A" : {}, "nBkgEventsErr_B" : {}, "nBkgEventsErr_C" : {}, "nBkgEventsErr_D" : {} 
        }

        for yBin in nYBins:

            yLowBinEdge = histBkg.GetYaxis().GetBinCenter(yBin)

            if (yLowBinEdge > float(finalDisc2Edge)): break

            yBinKey = "%.3f"%yLowBinEdge

            if finalDisc1Edge not in nTotSigCount_cdiGH["nSigEventsErr_A"]:
                nTotSigCount_cdiGH["nSigEvents_A"][finalDisc1Edge] = {}; nTotSigCount_cdiGH["nSigEventsErr_A"][finalDisc1Edge] = {}
                nTotSigCount_cdiGH["nSigEvents_B"][finalDisc1Edge] = {}; nTotSigCount_cdiGH["nSigEventsErr_B"][finalDisc1Edge] = {}
                nTotSigCount_cdiGH["nSigEvents_C"][finalDisc1Edge] = {}; nTotSigCount_cdiGH["nSigEventsErr_C"][finalDisc1Edge] = {}
                nTotSigCount_cdiGH["nSigEvents_D"][finalDisc1Edge] = {}; nTotSigCount_cdiGH["nSigEventsErr_D"][finalDisc1Edge] = {}

            if finalDisc1Edge not in nTotBkgCount_cdiGH["nBkgEventsErr_A"]:
                nTotBkgCount_cdiGH["nBkgEvents_A"][finalDisc1Edge] = {}; nTotBkgCount_cdiGH["nBkgEventsErr_A"][finalDisc1Edge] = {}
                nTotBkgCount_cdiGH["nBkgEvents_B"][finalDisc1Edge] = {}; nTotBkgCount_cdiGH["nBkgEventsErr_B"][finalDisc1Edge] = {}
                nTotBkgCount_cdiGH["nBkgEvents_C"][finalDisc1Edge] = {}; nTotBkgCount_cdiGH["nBkgEventsErr_C"][finalDisc1Edge] = {}
                nTotBkgCount_cdiGH["nBkgEvents_D"][finalDisc1Edge] = {}; nTotBkgCount_cdiGH["nBkgEventsErr_D"][finalDisc1Edge] = {}

            # count signal and  background events and errors in validation regions (c, di, G, H)
            nSigEventsErr_c = ROOT.Double(0.0); nSigEventsErr_di = ROOT.Double(0.0); nSigEventsErr_G = ROOT.Double(0.0); nSigEventsErr_H = ROOT.Double(0.0)
            nBkgEventsErr_c = ROOT.Double(0.0); nBkgEventsErr_di = ROOT.Double(0.0); nBkgEventsErr_G = ROOT.Double(0.0); nBkgEventsErr_H = ROOT.Double(0.0)

            nSigEvents_c  = ( histSig.IntegralAndError(xBin_final, lastXBin,     yBin, yBin_final, nSigEventsErr_c ) )
            nSigEvents_di = ( histSig.IntegralAndError(1,          xBin_final-1, yBin, yBin_final, nSigEventsErr_di) )
            nSigEvents_G  = ( histSig.IntegralAndError(xBin_final, lastXBin,     1,    yBin-1,     nSigEventsErr_G ) )
            nSigEvents_H  = ( histSig.IntegralAndError(1,          xBin_final-1, 1,    yBin-1,     nSigEventsErr_H ) )
            nBkgEvents_c  = ( histBkg.IntegralAndError(xBin_final, lastXBin,     yBin, yBin_final, nBkgEventsErr_c ) )
            nBkgEvents_di = ( histBkg.IntegralAndError(1,          xBin_final-1, yBin, yBin_final, nBkgEventsErr_di) )
            nBkgEvents_G  = ( histBkg.IntegralAndError(xBin_final, lastXBin,     1,    yBin-1,     nBkgEventsErr_G ) )
            nBkgEvents_H  = ( histBkg.IntegralAndError(1,          xBin_final-1, 1,    yBin-1,     nBkgEventsErr_H ) )

            nTotSigCount_cdiGH["nSigEvents_A"][finalDisc1Edge][yBinKey] = nSigEvents_c ; nTotSigCount_cdiGH["nSigEventsErr_A"][finalDisc1Edge][yBinKey] = nSigEventsErr_c
            nTotSigCount_cdiGH["nSigEvents_B"][finalDisc1Edge][yBinKey] = nSigEvents_di; nTotSigCount_cdiGH["nSigEventsErr_B"][finalDisc1Edge][yBinKey] = nSigEventsErr_di
            nTotSigCount_cdiGH["nSigEvents_C"][finalDisc1Edge][yBinKey] = nSigEvents_G ; nTotSigCount_cdiGH["nSigEventsErr_C"][finalDisc1Edge][yBinKey] = nSigEventsErr_G
            nTotSigCount_cdiGH["nSigEvents_D"][finalDisc1Edge][yBinKey] = nSigEvents_H ; nTotSigCount_cdiGH["nSigEventsErr_D"][finalDisc1Edge][yBinKey] = nSigEventsErr_H

            nTotBkgCount_cdiGH["nBkgEvents_A"][finalDisc1Edge][yBinKey] = nBkgEvents_c ; nTotBkgCount_cdiGH["nBkgEventsErr_A"][finalDisc1Edge][yBinKey] = nBkgEventsErr_c
            nTotBkgCount_cdiGH["nBkgEvents_B"][finalDisc1Edge][yBinKey] = nBkgEvents_di; nTotBkgCount_cdiGH["nBkgEventsErr_B"][finalDisc1Edge][yBinKey] = nBkgEventsErr_di
            nTotBkgCount_cdiGH["nBkgEvents_C"][finalDisc1Edge][yBinKey] = nBkgEvents_G ; nTotBkgCount_cdiGH["nBkgEventsErr_C"][finalDisc1Edge][yBinKey] = nBkgEventsErr_G
            nTotBkgCount_cdiGH["nBkgEvents_D"][finalDisc1Edge][yBinKey] = nBkgEvents_H ; nTotBkgCount_cdiGH["nBkgEventsErr_D"][finalDisc1Edge][yBinKey] = nBkgEventsErr_H
        
        return nTotSigCount_bdEF, nTotBkgCount_bdEF, nTotSigCount_cdiGH, nTotBkgCount_cdiGH

    # ---------------------------------------------------------------------------
    # make the validation region in BD : bdEF
    #   -- define validation disc1 edge between 0 and the final bin edge of disc1
    #   -- final disc2 edge is constant
    # ---------------------------------------------------------------------------
    def make_ValidationRegionEdges_bdEF(self, nTotSigCount_bdEF, nTotBkgCount_bdEF, final_nBkgEvents_A, final_nBkgEvents_C):
    
        disc1KeyOut_bdEF = [];  disc2KeyOut_bdEF = [];  finalDisc1Key_bdEF = -1.0; finalDisc2Key_bdEF = -1.0; validationMetric_bdEF  = 999.0
        finalSigFracb    = 0.0; finalSigFracd    = 0.0; finalSigFracE      = 0.0;  finalSigFracF      = 0.0;  nEvents_bd             = 0.0; nEvents_EF = 0.0
        sigFracsb        = [];  sigFracsd        = [];  sigFracsE          = [];   sigFracsF          = []
        
        closureErrsList_bdEF = []; closureErrUncList_bdEF = []
        weighted_Sig_b = []; weighted_Bkg_b = []; weighted_SigUnc_b = []; weighted_BkgUnc_b = []

        for disc1CutKey, disc2s in nTotBkgCount_bdEF["nBkgEvents_A"].items():
    
            for disc2Key, nEvents in disc2s.items():
    
                # number of signal and background events in aech b, d, E, F region
                nSigEvents_b = nTotSigCount_bdEF["nSigEvents_A"][disc1CutKey][disc2Key]; nBkgEvents_b = nTotBkgCount_bdEF["nBkgEvents_A"][disc1CutKey][disc2Key]
                nSigEvents_d = nTotSigCount_bdEF["nSigEvents_B"][disc1CutKey][disc2Key]; nBkgEvents_d = nTotBkgCount_bdEF["nBkgEvents_B"][disc1CutKey][disc2Key]
                nSigEvents_E = nTotSigCount_bdEF["nSigEvents_C"][disc1CutKey][disc2Key]; nBkgEvents_E = nTotBkgCount_bdEF["nBkgEvents_C"][disc1CutKey][disc2Key]
                nSigEvents_F = nTotSigCount_bdEF["nSigEvents_D"][disc1CutKey][disc2Key]; nBkgEvents_F = nTotBkgCount_bdEF["nBkgEvents_D"][disc1CutKey][disc2Key]
   
                nSigEventsErr_b = nTotSigCount_bdEF["nSigEventsErr_A"][disc1CutKey][disc2Key]; nBkgEventsErr_b = nTotBkgCount_bdEF["nBkgEventsErr_A"][disc1CutKey][disc2Key]
                nSigEventsErr_d = nTotSigCount_bdEF["nSigEventsErr_B"][disc1CutKey][disc2Key]; nBkgEventsErr_d = nTotBkgCount_bdEF["nBkgEventsErr_B"][disc1CutKey][disc2Key]
                nSigEventsErr_E = nTotSigCount_bdEF["nSigEventsErr_C"][disc1CutKey][disc2Key]; nBkgEventsErr_E = nTotBkgCount_bdEF["nBkgEventsErr_C"][disc1CutKey][disc2Key]
                nSigEventsErr_F = nTotSigCount_bdEF["nSigEventsErr_D"][disc1CutKey][disc2Key]; nBkgEventsErr_F = nTotBkgCount_bdEF["nBkgEventsErr_D"][disc1CutKey][disc2Key]
 
                # Signal fractions in each b, d, E, F region
                nTot_SigBkg_b = nSigEvents_b + nBkgEvents_b 
                nTot_SigBkg_d = nSigEvents_d + nBkgEvents_d
                nTot_SigBkg_E = nSigEvents_E + nBkgEvents_E
                nTot_SigBkg_F = nSigEvents_F + nBkgEvents_F
    
                tempSigFracsb = -1.0; tempSigFracsd = -1.0; tempSigFracsE = -1.0; tempSigFracsF = -1.0
    
                if nTot_SigBkg_b > 0.0: 
                    tempSigFracsb = nSigEvents_b / nTot_SigBkg_b 
                if nTot_SigBkg_d > 0.0: 
                    tempSigFracsd = nSigEvents_d / nTot_SigBkg_d
                if nTot_SigBkg_E > 0.0: 
                    tempSigFracsE = nSigEvents_E / nTot_SigBkg_E
                if nTot_SigBkg_F > 0.0: 
                    tempSigFracsF = nSigEvents_F / nTot_SigBkg_F

                # get the closure error and closure error unc. to plot variable vs disc as 1D 
                tempClosureErr = -999.0; tempClosureErrUnc = -999.0;
                if nBkgEvents_b > 0.0 and nBkgEvents_F > 0.0: 
                    tempClosureErr, tempClosureErrUnc = self.cal_ClosureError(nBkgEvents_b, nBkgEvents_d, nBkgEvents_E, nBkgEvents_F, nBkgEventsErr_b, nBkgEventsErr_d, nBkgEventsErr_E, nBkgEventsErr_F)

                if tempClosureErr > 0.0:
                    closureErrsList_bdEF.append(abs(tempClosureErr))
                    closureErrUncList_bdEF.append(tempClosureErrUnc)
                    disc1KeyOut_bdEF.append(float(disc1CutKey))
                    disc2KeyOut_bdEF.append(float(disc2Key))

                    # get the weighted signal and backgoround events to plot variable vs disc as 1D 
                    weighted_Sig_b.append(nSigEvents_b)
                    weighted_Bkg_b.append(nBkgEvents_b)   
                    weighted_SigUnc_b.append(nBkgEventsErr_b)    
                    weighted_BkgUnc_b.append(nBkgEventsErr_b) 

                    sigFracsb.append(float(tempSigFracsb))
                    sigFracsd.append(float(tempSigFracsd))
                    sigFracsE.append(float(tempSigFracsE))
                    sigFracsF.append(float(tempSigFracsF))
    
                # calculate the validation metric
                temp_ValidationMetric = 999.0
                temp_ValidationMetric = self.cal_MetricOfValidationRegion_bdEF(final_nBkgEvents_A, final_nBkgEvents_C, nBkgEvents_b, nBkgEvents_d)
    
                if temp_ValidationMetric < validationMetric_bdEF:
                    finalDisc1Key_bdEF    = disc1CutKey
                    finalDisc2Key_bdEF    = disc2Key
                    validationMetric_bdEF = temp_ValidationMetric
                    finalSigFracb         = tempSigFracsb
                    finalSigFracd         = tempSigFracsd
                    finalSigFracE         = tempSigFracsE
                    finalSigFracF         = tempSigFracsF
                    nEvents_bd            = nBkgEvents_b + nBkgEvents_d
                    nEvents_EF            = nBkgEvents_E + nBkgEvents_F

        return finalDisc1Key_bdEF, finalDisc2Key_bdEF, disc1KeyOut_bdEF, disc2KeyOut_bdEF, closureErrsList_bdEF, closureErrUncList_bdEF, weighted_Sig_b, weighted_Bkg_b, weighted_SigUnc_b, weighted_BkgUnc_b, sigFracsb, sigFracsd, sigFracsE, sigFracsF, finalSigFracb, finalSigFracd, finalSigFracE, finalSigFracF, nEvents_bd, nEvents_EF

    # ---------------------------------------------------------------------------
    # make the validation region in CD : cdiGH
    #   -- define validation disc2 edge between 0 and the final bin edge of disc2
    #   -- final disc1 edge is constant
    # --------------------------------------------------------------------------- 
    def make_ValidationRegionEdges_cdiGH(self, nTotSigCount_cdiGH, nTotBkgCount_cdiGH, final_nBkgEvents_A, final_nBkgEvents_B):    

        disc1KeyOut_cdiGH = [];  disc2KeyOut_cdiGH = [];  finalDisc1Key_cdiGH = -1.0; finalDisc2Key_cdiGH = -1.0; validationMetric_cdiGH  = 999.0
        finalSigFracc     = 0.0; finalSigFracdi    = 0.0; finalSigFracG       = 0.0;  finalSigFracH       = 0.0;  nEvents_cdi             = 0.0; nEvents_GH = 0.0
        sigFracsc         = [];  sigFracsdi        = [];  sigFracsG           = [];   sigFracsH           = []

        closureErrsList_cdiGH = []; closureErrUncList_cdiGH = []
        weighted_Sig_c = []; weighted_Bkg_c = []; weighted_SigUnc_c = []; weighted_BkgUnc_c = []

        for disc1Key, disc2s in nTotBkgCount_cdiGH["nBkgEvents_A"].items():

            for disc2CutKey, nEvents in disc2s.items():

                # number of signal and background events in aech c, di, G, H region
                nSigEvents_c  = nTotSigCount_cdiGH["nSigEvents_A"][disc1Key][disc2CutKey]; nBkgEvents_c  = nTotBkgCount_cdiGH["nBkgEvents_A"][disc1Key][disc2CutKey] 
                nSigEvents_di = nTotSigCount_cdiGH["nSigEvents_B"][disc1Key][disc2CutKey]; nBkgEvents_di = nTotBkgCount_cdiGH["nBkgEvents_B"][disc1Key][disc2CutKey]
                nSigEvents_G  = nTotSigCount_cdiGH["nSigEvents_C"][disc1Key][disc2CutKey]; nBkgEvents_G  = nTotBkgCount_cdiGH["nBkgEvents_C"][disc1Key][disc2CutKey]
                nSigEvents_H  = nTotSigCount_cdiGH["nSigEvents_D"][disc1Key][disc2CutKey]; nBkgEvents_H  = nTotBkgCount_cdiGH["nBkgEvents_D"][disc1Key][disc2CutKey]

                nSigEventsErr_c  = nTotSigCount_cdiGH["nSigEventsErr_A"][disc1Key][disc2CutKey]; nBkgEventsErr_c  = nTotBkgCount_cdiGH["nBkgEventsErr_A"][disc1Key][disc2CutKey]
                nSigEventsErr_di = nTotSigCount_cdiGH["nSigEventsErr_B"][disc1Key][disc2CutKey]; nBkgEventsErr_di = nTotBkgCount_cdiGH["nBkgEventsErr_B"][disc1Key][disc2CutKey]
                nSigEventsErr_G  = nTotSigCount_cdiGH["nSigEventsErr_C"][disc1Key][disc2CutKey]; nBkgEventsErr_G  = nTotBkgCount_cdiGH["nBkgEventsErr_C"][disc1Key][disc2CutKey]
                nSigEventsErr_H  = nTotSigCount_cdiGH["nSigEventsErr_D"][disc1Key][disc2CutKey]; nBkgEventsErr_H  = nTotBkgCount_cdiGH["nBkgEventsErr_D"][disc1Key][disc2CutKey]

                # Signal fractions in each c, di, G, H region
                nTot_SigBkg_c  = nSigEvents_c  + nBkgEvents_c
                nTot_SigBkg_di = nSigEvents_di + nBkgEvents_di
                nTot_SigBkg_G  = nSigEvents_G  + nBkgEvents_G
                nTot_SigBkg_H  = nSigEvents_H  + nBkgEvents_H

                tempSigFracsc = -1.0; tempSigFracsdi = -1.0; tempSigFracsG = -1.0; tempSigFracsH = -1.0

                if nTot_SigBkg_c > 0.0:
                    tempSigFracsc = nSigEvents_c  / nTot_SigBkg_c
                if nTot_SigBkg_di > 0.0:
                    tempSigFracsdi = nSigEvents_di / nTot_SigBkg_di
                if nTot_SigBkg_G > 0.0:
                    tempSigFracsG = nSigEvents_G  / nTot_SigBkg_G
                if nTot_SigBkg_H > 0.0:
                    tempSigFracsH = nSigEvents_H  / nTot_SigBkg_H

                # get the closure error and closure error unc. to plot variable vs disc as 1D 
                tempClosureErr = -999.0; tempClosureErrUnc = -999.0;
                if nBkgEvents_c > 0.0 and nBkgEvents_H > 0.0:
                    tempClosureErr, tempClosureErrUnc = self.cal_ClosureError(nBkgEvents_c, nBkgEvents_di, nBkgEvents_G, nBkgEvents_H, nBkgEventsErr_c, nBkgEventsErr_di, nBkgEventsErr_G, nBkgEventsErr_H)

                if tempClosureErr > 0.0:
                    closureErrsList_cdiGH.append(abs(tempClosureErr))
                    closureErrUncList_cdiGH.append(tempClosureErrUnc)
                    disc1KeyOut_cdiGH.append(float(disc1Key))
                    disc2KeyOut_cdiGH.append(float(disc2CutKey))

                    # get the weighted signal and backgoround events to plot variable vs disc as 1D 
                    weighted_Sig_c.append(nSigEvents_c)
                    weighted_Bkg_c.append(nBkgEvents_c)
                    weighted_SigUnc_c.append(nBkgEventsErr_c)
                    weighted_BkgUnc_c.append(nBkgEventsErr_c)

                    sigFracsc.append(float(tempSigFracsc))
                    sigFracsdi.append(float(tempSigFracsdi))
                    sigFracsG.append(float(tempSigFracsG))
                    sigFracsH.append(float(tempSigFracsH))

                # calculate the validation metric
                temp_ValidationMetric = 999.0
                temp_ValidationMetric = self.cal_MetricOfValidationRegion_cdiGH(final_nBkgEvents_A, final_nBkgEvents_B, nBkgEvents_c, nBkgEvents_di)

                if temp_ValidationMetric < validationMetric_cdiGH:
                    finalDisc1Key_cdiGH    = disc1Key
                    finalDisc2Key_cdiGH    = disc2CutKey
                    validationMetric_cdiGH = temp_ValidationMetric
                    finalSigFracc          = tempSigFracsc
                    finalSigFracdi         = tempSigFracsdi
                    finalSigFracG          = tempSigFracsG
                    finalSigFracH          = tempSigFracsH
                    nEvents_cdi            = nBkgEvents_c + nBkgEvents_di
                    nEvents_GH             = nBkgEvents_G + nBkgEvents_H

        return finalDisc1Key_cdiGH, finalDisc2Key_cdiGH, disc1KeyOut_cdiGH, disc2KeyOut_cdiGH, closureErrsList_cdiGH, closureErrUncList_cdiGH, weighted_Sig_c, weighted_Bkg_c, weighted_SigUnc_c, weighted_BkgUnc_c, sigFracsc, sigFracsdi, sigFracsG, sigFracsH, finalSigFracc, finalSigFracdi, finalSigFracG, finalSigFracH, nEvents_cdi, nEvents_GH

