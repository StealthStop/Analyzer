import ROOT
import os
import math
import argparse
import numpy as np

import matplotlib as mpl
mpl.use('Agg')

import matplotlib.lines as ml
import matplotlib.pyplot as plt

from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

from DoubleDisCo_Regions import *
from DoubleDisCo_Plotter import *


class ttVariances_Plotters:

    def __init__(self, outputDir, metric, year, channel):
        self.outputDir = outputDir
        self.metric    = metric 
        self.year      = year
        self.channel   = channel


    # ----------------------
    # calculate all closures
    # ----------------------
    def cal_simpleClosure(self, nEvents_A, nEvents_B, nEvents_C, nEvents_D, nEventsErr_A, nEventsErr_B, nEventsErr_C, nEventsErr_D):

        if nEvents_D == 0.0:
            return -999.0, -999.0

        nPred_A    = (nEvents_B * nEvents_C) / nEvents_D
        nPred_Aunc = ((nEvents_C * nEventsErr_B / nEvents_D)**2.0 + (nEventsErr_C * nEvents_B / nEvents_D)**2.0 + (nEvents_C * nEvents_B * nEventsErr_D / nEvents_D**2.0)**2.0)**0.5

        return nPred_A, nPred_Aunc

    # -----------------
    # plot all closures
    # -----------------
    def plot_ttVar_ClosureNjets(self, bkgObs, bkgPred, bkgObs_Var, bkgPred_Var, Njets, name = '', closureTag = '', varTag = ''):
        
        x             = []; xUnc             = []
        obs           = []; obsUnc           = [] 
        pred          = []; predUnc          = []
        obs_var       = []; obsUnc_var       = []
        pred_var      = []; predUnc_var      = []
        abcdError     = []; abcdErrorUnc     = []
        abcdError_var = []; abcdErrorUnc_var = []
        
        Ratio1 = []; Ratio1Unc = []; Ratio2 = []; Ratio2Unc = []
    
        totalChi2 = 0.0; ndof = 0
    
        # ------------------------------------
        # calculate pull, ratio and ratio unc.
        # ------------------------------------
        def get_ratioUnc(numerator, denominator):
    
            pull     = ( numerator[0] - denominator[0] ) / ( (numerator[1])**2 + (denominator[1])**2 )**0.5
            ratio    = ( numerator[0] / denominator[0] )
            ratioUnc = ( (numerator[1] / denominator[0])**2.0 + ( denominator[1] * numerator[0] / denominator[0]**2.0 )**2.0 )**0.5
    
            return pull, ratio, ratioUnc
    
        # -------------------
        # loop over njet bins
        # -------------------
        for i in range(0, len(Njets)):
       
            if bkgObs[i][1] != 0.0:
      
                #
                x.append(float(Njets[i]))
                xUnc.append(0.5) 
                
                # closure for tt
                pred.append(bkgPred[i][0])
                predUnc.append(bkgPred[i][1])
                obs.append(bkgObs[i][0])
                obsUnc.append(bkgObs[i][1])
    
                # closure for tt variances
                pred_var.append(bkgPred_Var[i][0])
                predUnc_var.append(bkgPred_Var[i][1])
                obs_var.append(bkgObs_Var[i][0])
                obsUnc_var.append(bkgObs_Var[i][1])
    
                # non-closure for tt
                closureError    = 1.0 - ( bkgPred[i][0] / bkgObs[i][0] )
                closureErrorUnc = ((bkgPred[i][1] / bkgObs[i][0])**2.0 + (bkgObs[i][1] * bkgPred[i][0] / bkgObs[i][0]**2.0)**2.0)**0.5
                abcdError.append(closureError)
                abcdErrorUnc.append(closureErrorUnc)
    
                # non-closure for tt variances
                closureError_var    = 1.0 - ( bkgPred_Var[i][0] / bkgObs_Var[i][0] )
                closureErrorUnc_var = ((bkgPred_Var[i][1] / bkgObs_Var[i][0])**2.0 + (bkgObs_Var[i][1] * bkgPred_Var[i][0] / bkgObs_Var[i][0]**2.0)**2.0)**0.5
                abcdError_var.append(closureError_var)
                abcdErrorUnc_var.append(closureErrorUnc_var)
    
                # MC correction for tt
                _, ratio1, ratio1Unc = get_ratioUnc(bkgObs[i], bkgPred[i]) 
                MC_correction = ratio1
                Ratio1.append(ratio1)
                Ratio1Unc.append(ratio1Unc)
    
                # MC correction for tt variances
                _, ratio2, ratio2Unc = get_ratioUnc(bkgObs_Var[i], bkgPred_Var[i])
                MC_correction_var = ratio2
                Ratio2.append(ratio2)
                Ratio2Unc.append(ratio2Unc)
    
                #
                pull, _, _ = get_ratioUnc(bkgPred[i], bkgObs[i])
                totalChi2 += pull ** 2.0
                ndof      += 1  
    
        # ------------
        # plot closure
        # ------------ 
        fig = plt.figure(figsize=(15,15))
        ax  = fig.add_subplot(111)
        fig.subplots_adjust(hspace=0.4, left=0.15, right=0.95, top=0.94)       
 
        lowerNjets  = float(Njets[0])
        higherNjets = float(Njets[(-1)])
        
        ax.spines['top'].set_color('none')
        ax.spines['bottom'].set_color('none')
        ax.spines['left'].set_color('none')
        ax.spines['right'].set_color('none')
        ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
    
        # closure 
        ax1 = fig.add_subplot(14, 1, (1,6))
        fig.subplots_adjust(hspace=0.8, left=0.15, right=0.95)
        ax1.set_xlim([lowerNjets - 0.5, higherNjets + 0.5])
        #ax1.text(0.05, 0.1, '$\\chi^2$ / ndof = %3.2f' % (totalChi2 / float(ndof)), horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes, fontsize=10)
        ax1.text(0.16, 1.065, 'CMS',                      transform=ax.transAxes, fontsize=40, fontweight='bold',   va='top', ha='right')
        ax1.text(0.50, 1.055, 'Preliminary',              transform=ax.transAxes, fontsize=36, fontstyle='italic',  va='top', ha='right')
        ax1.text(1.0,  1.055, '%s (13 TeV)' %(self.year), transform=ax.transAxes, fontsize=34, fontweight='normal', va='top', ha='right')
        ax1.errorbar(x, pred,     yerr=predUnc,     label='Predicted TT',          xerr=xUnc, fmt='', capsize=0, color='red',     lw=0, elinewidth=2, marker='o', markeredgecolor='red',     markerfacecolor='red',     markersize=5.0)
        ax1.errorbar(x, obs,      yerr=obsUnc,      label='Observed TT',           xerr=xUnc, fmt='', capsize=0, color='black',   lw=0, elinewidth=2, marker='o', markeredgecolor='black',   markerfacecolor='black',   markersize=5.0)
        ax1.errorbar(x, pred_var, yerr=predUnc_var, label='Predicted %s'%(varTag), xerr=xUnc, fmt='', capsize=0, color='cyan',    lw=0, elinewidth=2, marker='o', markeredgecolor='cyan',    markerfacecolor='cyan',    markersize=5.0)
        ax1.errorbar(x, obs_var,  yerr=obsUnc_var,  label='Observed %s'%(varTag),  xerr=xUnc, fmt='', capsize=0, color='dimgrey', lw=0, elinewidth=2, marker='o', markeredgecolor='dimgrey', markerfacecolor='dimgrey', markersize=5.0)
        ax1.set_xticklabels([])
        ax1.set_ylabel('Unweighted Event Counts', fontsize=28)
        ax1.set_yscale('log')
        ax1.set_ylim([5.0, 2e4])
 
        # non-closure
        ax2 = fig.add_subplot(14, 1, (7,10))
        fig.subplots_adjust(hspace=0.8, left=0.15, right=0.95)
        ax2.set_xlim([lowerNjets - 0.5, higherNjets + 0.5])
        ax2.errorbar(x, abcdError,     yerr=abcdErrorUnc,     label='Non-Closure TT',          xerr=xUnc, fmt='', capsize=0, color='cornflowerblue', lw=0, elinewidth=3, marker='o', markeredgecolor='cornflowerblue', markerfacecolor='cornflowerblue', markersize=6.0)
        ax2.errorbar(x, abcdError_var, yerr=abcdErrorUnc_var, label='Non-Closure %s'%(varTag), xerr=xUnc, fmt='', capsize=0, color='blue',           lw=0, elinewidth=3, marker='o', markeredgecolor='blue',           markerfacecolor='blue',           markersize=6.0)
        ax2.axhline(y=0.0, color='black', linestyle='dashed', lw=1.5)
        ax2.grid(axis='y', color='black', linestyle='dashed', which='both')
        ax2.set_xticklabels([])
        ax2.set_ylabel('1 - Pred./Obs.', fontsize=20)
        ax2.set_ylim([-0.5, 0.5])   
    
        # MC-correction
        ax3 = fig.add_subplot(14, 1, (11,14))
        fig.subplots_adjust(hspace=0.8, left=0.15, right=0.95)
        ax3.set_xlim([lowerNjets - 0.5, higherNjets + 0.5])
        ax3.errorbar(x, Ratio1, yerr=Ratio1Unc, xerr=xUnc, label='MC corr. TT',          fmt='', capsize=0, color='mediumseagreen', lw=0, elinewidth=3, marker='o', markeredgecolor='mediumseagreen', markerfacecolor='mediumseagreen', markersize=6.0)
        ax3.errorbar(x, Ratio2, yerr=Ratio2Unc, xerr=xUnc, label='MC corr. %s'%(varTag), fmt='', capsize=0, color='darkgreen',      lw=0, elinewidth=3, marker='o', markeredgecolor='darkgreen',      markerfacecolor='darkgreen',      markersize=6.0)
        ax3.axhline(y=1.0, color='black', linestyle='dashed', lw=1.5)
        ax3.grid(axis='y', color='black', linestyle='dashed', which='both')
        ax3.set_ylabel('(Obs./Pred.)', fontsize=20)
        ax3.set_ylim([0.5, 1.5])
    
        plt.xticks([int(Njet) for Njet in Njets], fontsize=26)
        plt.xlabel('Number of jets', fontsize=28)    
   
        ax1.legend(loc='upper right', numpoints=1, frameon=False, fontsize=20)
        ax2.legend(loc='upper right', numpoints=1, frameon=False, fontsize=20)
        ax3.legend(loc='upper right', numpoints=1, frameon=False, fontsize=20) 
        
        fig.savefig('%s/%s_Njets_Region_A_PredVsActual_%s_%s_%s_%s_%s.pdf' % (self.outputDir, self.year, name, self.channel, self.metric, closureTag, varTag))
       
        plt.close(fig)

    # ----------------------
    # make all closure plots
    # ----------------------
    def make_ttVariances_allClosures(self, edgesPerNjets=None, bkgEventsPerNjets=None, bkgVarEventsPerNjets=None, Njets=None, name="", closureTag="", varTag=""):

        evtsNjets     = []; evtsNjetsPred     = []
        evtsNjets_Var = []; evtsNjetsPred_Var = []    

        for Njet in Njets:
    
            if edgesPerNjets[Njet][name][0] == -1.0 or edgesPerNjets[Njet][name][1] == -1.0:
    
                evtsNjets.append((0.0, 0.0))
                evtsNjetsPred.append((0.0, 0.0))
                evtsNjets_Var.append((0.0, 0.0))
                evtsNjetsPred_Var.append((0.0, 0.0))
    
            else:

                # for TT
                evtsA = bkgEventsPerNjets[Njet][name]["A"][0]; evtsAunc = bkgEventsPerNjets[Njet][name]["A"][1]**2.0
                evtsB = bkgEventsPerNjets[Njet][name]["B"][0]; evtsBunc = bkgEventsPerNjets[Njet][name]["B"][1]**2.0
                evtsC = bkgEventsPerNjets[Njet][name]["C"][0]; evtsCunc = bkgEventsPerNjets[Njet][name]["C"][1]**2.0
                evtsD = bkgEventsPerNjets[Njet][name]["D"][0]; evtsDunc = bkgEventsPerNjets[Njet][name]["D"][1]**2.0

                pred_A, predUnc_A = self.cal_simpleClosure(evtsA, evtsB, evtsC, evtsD, evtsAunc**0.5, evtsBunc**0.5, evtsCunc**0.5, evtsDunc**0.5)

                evtsNjets.append((evtsA, evtsAunc**0.5))
                evtsNjetsPred.append((pred_A, predUnc_A))
    
                # for TT Variances
                evtsA_var = bkgVarEventsPerNjets[Njet][name]["A"][0]; evtsAunc_var = bkgVarEventsPerNjets[Njet][name]["A"][1]**2.0
                evtsB_var = bkgVarEventsPerNjets[Njet][name]["B"][0]; evtsBunc_var = bkgVarEventsPerNjets[Njet][name]["B"][1]**2.0
                evtsC_var = bkgVarEventsPerNjets[Njet][name]["C"][0]; evtsCunc_var = bkgVarEventsPerNjets[Njet][name]["C"][1]**2.0
                evtsD_var = bkgVarEventsPerNjets[Njet][name]["D"][0]; evtsDunc_var = bkgVarEventsPerNjets[Njet][name]["D"][1]**2.0

                predVar_A, predUncVar_A = self.cal_simpleClosure(evtsA_var, evtsB_var, evtsC_var, evtsD_var, evtsAunc_var**0.5, evtsBunc_var**0.5, evtsCunc_var**0.5, evtsDunc_var**0.5)

                evtsNjets_Var.append((evtsA_var, evtsAunc_var**0.5))
                evtsNjetsPred_Var.append((predVar_A, predUncVar_A))

        self.plot_ttVar_ClosureNjets(np.array(evtsNjets), np.array(evtsNjetsPred), np.array(evtsNjets_Var), np.array(evtsNjetsPred_Var), Njets, name, closureTag, varTag)

