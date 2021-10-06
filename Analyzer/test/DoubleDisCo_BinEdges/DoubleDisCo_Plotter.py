import math
import numpy as np

import matplotlib as mpl
mpl.use('Agg')

import matplotlib.lines as ml
import matplotlib.pyplot as plt

from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

class Common_Calculations_Plotters:

    def __init__(self, outputDir, metric, year, model, mass, channel):
        self.metric  = metric 
        self.year    = year
        self.model   = model
        self.mass    = mass
        self.channel = channel

        self.outputDir = outputDir

    # ----------------------
    # calculate all closures
    # ----------------------
    def cal_simpleClosure_ABCD(self, nEvents_A, nEvents_B, nEvents_C, nEvents_D, nEventsErr_A, nEventsErr_B, nEventsErr_C, nEventsErr_D):

        if nEvents_D == 0.0:
            return -999.0, -999.0

        nPred_A    = (nEvents_B * nEvents_C) / nEvents_D
        nPred_Aunc = ((nEvents_C * nEventsErr_B / nEvents_D)**2.0 + (nEventsErr_C * nEvents_B / nEvents_D)**2.0 + (nEvents_C * nEvents_B * nEventsErr_D / nEvents_D**2.0)**2.0)**0.5

        return nPred_A, nPred_Aunc

    # -----------------
    # plot all closures
    # -----------------
    def plot_ClosureNjets(self, bkgObs, bkgPred, Njets, name = '', closureTag = '', bkgTag = '', isBlind=False):
        
        x         = []; xUnc         = []
        obs       = []; obsUnc       = [] 
        pred      = []; predUnc      = []
        abcdError = []; abcdErrorUnc = []
        abcdPull  = []; abcdPullUnc  = []
        pullDenom = []; zeros        = []
        
        totalChi2 = 0.0; ndof = 0

        # ------------------------------------------   
        # make a band with (pull denominator) / pred
        # ------------------------------------------  
        def makeErrorBoxes(x, zeros, xUnc, pullDenom):

            # list for all the error patches
            errorBoxes = []

            # loop over data points / create box from errors at each point
            for xc, yc, xe, ye in zip(x, zeros, xUnc.T, pullDenom.T):
                rect = Rectangle( (xc-xe[0], yc-ye[0]), xe.sum(), ye.sum() )
                errorBoxes.append(rect)

            # create patch collection with specified colour/alpha
            pc = PatchCollection(errorBoxes, facecolor='mediumslateblue', alpha=0.5, edgecolor='none')

            return pc
 
        # -------------------
        # loop over the njets
        # -------------------
        for i in range(0, len(Njets)):
    
            if bkgObs[i][1] != 0.0:
    
                x.append(float(Njets[i]))  
                xUnc.append(0.5)

                pull            = (bkgPred[i][0] - bkgObs[i][0]) / ( bkgPred[i][1]**2 + bkgObs[i][1]**2 )**0.5
                pullUnc         = 1.0
                closureError    = 1.0 - ( bkgPred[i][0] / bkgObs[i][0] )
                closureErrorUnc = ((bkgPred[i][1] / bkgObs[i][0])**2.0 + (bkgObs[i][1] * bkgPred[i][0] / bkgObs[i][0]**2.0)**2.0)**0.5
                pullDenominator = ( bkgPred[i][1]**2 + bkgObs[i][1]**2 )**0.5 / bkgPred[i][0]
                zero            = 0.0

                # for the data closure / blind the 9-10-11 njet bins
                if isBlind and i >= 2:
                    pred.append(0.0)
                    predUnc.append(0.0)
                    obs.append(0.0)
                    obsUnc.append(0.0)
                    abcdPull.append(999)
                    abcdPullUnc.append(0)
                    abcdError.append(999)
                    abcdErrorUnc.append(0)
                    pullDenom.append(0)
                    zeros.append(999)

                else:
                    pred.append(bkgPred[i][0])
                    predUnc.append(bkgPred[i][1])
                    obs.append(bkgObs[i][0])
                    obsUnc.append(bkgObs[i][1])
                    abcdPull.append(pull)
                    abcdPullUnc.append(pullUnc)
                    abcdError.append(closureError)
                    abcdErrorUnc.append(closureErrorUnc)
                    pullDenom.append(pullDenominator)
                    zeros.append(zero)
                    totalChi2 += pull ** 2.0
                    ndof      += 1                    

        # ------------------
        # plot usual closure
        # ------------------
        fig = plt.figure(figsize=(6, 6))
        ax  = fig.add_subplot(111)
        fig.subplots_adjust(hspace=0, left=0.15, right=0.85, top=0.94)
    
        lowerNjets  = float(Njets[0])
        higherNjets = float(Njets[(-1)])
    
        ax.spines['top'].set_color('none')
        ax.spines['bottom'].set_color('none')
        ax.spines['left'].set_color('none')
        ax.spines['right'].set_color('none')
        ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
    
        # closure plot
        ax1 = fig.add_subplot(4, 1, (1, 2))
        fig.subplots_adjust(left=0.15, right=0.95)
        ax1.set_yscale('log')
        ax1.set_xlim([lowerNjets - 0.5, higherNjets + 0.5])
        ax1.text(0.05, 0.1, '$\\chi^2$ / ndof = %3.2f' % (totalChi2 / float(ndof)), horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes, fontsize=10)
        #ax1.text(0.05, 0.25,  '%s Metric'%(self.metric),  horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes, fontsize=14, fontweight='bold')
        ax1.text(0.16, 1.065, 'CMS',                      transform=ax.transAxes, fontsize=20, fontweight='bold',   va='top', ha='right')
        ax1.text(0.50, 1.055, 'Preliminary',              transform=ax.transAxes, fontsize=16, fontstyle='italic',  va='top', ha='right')
        ax1.text(1.0,  1.055, '%s (13 TeV)' %(self.year), transform=ax.transAxes, fontsize=14, fontweight='normal', va='top', ha='right')
        ax1.set_ylabel('Unweighted Event Counts')
        ax1.errorbar(x, pred, yerr=predUnc, label='Predicted', xerr=xUnc, fmt='', capsize=0, color='red',   lw=0, elinewidth=2, marker='o', markeredgecolor='red',   markerfacecolor='red',   markersize=5.0)
        ax1.errorbar(x, obs,  yerr=obsUnc,  label='Observed',  xerr=xUnc, fmt='', capsize=0, color='black', lw=0, elinewidth=2, marker='o', markeredgecolor='black', markerfacecolor='black', markersize=5.0)

        # simple ratio plot    
        ax2 = fig.add_subplot(4, 1, 3)
        ax2.set_xlim([lowerNjets - 0.5, higherNjets + 0.5])
        ax2.errorbar(x, abcdError, yerr=abcdErrorUnc, xerr=xUnc, fmt='', capsize=0, color='blue', lw=0, elinewidth=2, marker='o', markeredgecolor='blue', markerfacecolor='blue', markersize=5.0)
        ax2.axhline(y=0.0, color='black', linestyle='dashed', lw=1.5)
        ax2.grid(axis='y', color='black', linestyle='dashed', which='both')
        ax2.set_ylabel('1 - Pred./Obs.')
        ax2.set_ylim([-1.4, 1.4])
   
        # pull plot
        ax3 = fig.add_subplot(4, 1, 4)
        ax3.set_xlim([lowerNjets - 0.5, higherNjets + 0.5])
        ax3.errorbar(x, abcdPull, yerr=abcdPullUnc, xerr=xUnc, fmt='', capsize=0, color='purple', lw=0, elinewidth=2, marker='o', markeredgecolor='purple', markerfacecolor='purple', markersize=5.0)
        ax3.axhline(y=0.0, color='black', linestyle='dashed', lw=1.5)
        ax3.grid(axis='y', color='black', linestyle='dashed', which='both')
        ax3.set_xlabel('Number of jets')
        ax3.set_ylabel('Pull') # (Pred - Obs) / ObsUnc') 
        ax3.set_ylim([-5.4, 5.4])
 
        # ( obsUnc^2 + predUnc^2 )^0.5 / pred
        pc = makeErrorBoxes(np.array(x), np.array(zeros), np.array([xUnc, xUnc]), np.array([pullDenom, pullDenom]))
        ax4 = ax3.twinx()
        ax4.set_ylabel('pullDen / pred', color='mediumslateblue')
        ax4.add_collection(pc)
        ax4.set_ylim([-0.2, 0.2])        
        #ax4.set_yscale('log')
        plt.show()

        plt.xticks([int(Njet) for Njet in Njets])
    
        ax1.legend(loc='best', numpoints=1, frameon=False)
    
        fig.savefig('%s/%s_Njets_Region_A_PredVsActual_%s_%s_%s_%s_%s.pdf' % (self.outputDir, self.year, name, self.channel, self.metric, closureTag, bkgTag))
   
        plt.close(fig)


    # ----------------------------------------
    # plot all closures for data-MC comparison
    #   -- MC is TT
    #   -- MC is TT + NonTT
    # ----------------------------------------
    def plot_dataVsMC_ClosureNjets(self, bkgPred, bkgObs, dataPred, dataObs, Njets, name = '', closureTag = '', bkgTag = ''):

        x         = []; xUnc         = []
        bkg_obs   = []; bkg_obsUnc   = [] 
        bkg_pred  = []; bkg_predUnc  = []
        data_obs  = []; data_obsUnc  = [] 
        data_pred = []; data_predUnc = []
        PullUnc   = [] 
        Pull1     = []; Pull2     = []; Pull3     = []; Pull4     = []; Pull5     = []; Pull6     = []
        Ratio1    = []; Ratio2    = []; Ratio3    = []; Ratio4    = []; Ratio5    = []; Ratio6    = []
        Ratio1Unc = []; Ratio2Unc = []; Ratio3Unc = []; Ratio4Unc = []; Ratio5Unc = []; Ratio6Unc = [] 
        Ratio     = []; RatioUnc  = []; DataPred_corrected = []

        # ------------------------------------
        # calculate pull, ratio and ratio unc.
        # ------------------------------------
        def cal_ratioUnc(numerator, denominator):

            pull     = ( numerator[0] - denominator[0] ) / ( (numerator[1])**2 + (denominator[1])**2 )**0.5
            ratio    = ( numerator[0] / denominator[0] )
            ratioUnc = ( (numerator[1] / denominator[0])**2.0 + ( denominator[1] * numerator[0] / denominator[0]**2.0 )**2.0 )**0.5

            return pull, ratio, ratioUnc

        # -------------------
        # loop over njet bins
        # -------------------
        for i in range(0, len(Njets)):

            if bkgObs[i][1] != 0.0:

                x.append(float(Njets[i]))
                xUnc.append(0.5)

                # bkg
                bkg_pred.append(bkgPred[i][0])
                bkg_predUnc.append(bkgPred[i][1])
                bkg_obs.append(bkgObs[i][0])
                bkg_obsUnc.append(bkgObs[i][1])

                # ------------------------------ 
                # pull  = (Obs - Pred) / PredUnc
                # ratio = (Obs / Pred) 
                # ratioUnc = ((bkgObs[i][1] / bkgPred[i][0])**2.0 + (bkgPred[i][1] * bkgObs[i][0] / bkgPred[i][0]**2.0)**2.0)**0.5
                # ------------------------------
                # bkgObs / bkgPred 
                pullUnc   = 1.0
                pull1, ratio1, ratio1Unc = cal_ratioUnc(bkgObs[i], bkgPred[i]) 
            
                # dataObs / dataPred 
                pull2, ratio2, ratio2Unc = cal_ratioUnc(dataObs[i], dataPred[i]) 

                # dataPred / bkgPred
                pull3, ratio3, ratio3Unc = cal_ratioUnc(dataPred[i], bkgPred[i]) 

                # dataObs / bkgObs
                pull4, ratio4, ratio4Unc = cal_ratioUnc(dataObs[i], bkgObs[i])                

                # dataPred / bkgObs
                pull5, ratio5, ratio5Unc = cal_ratioUnc(dataPred[i], bkgObs[i]) 

                # dataObs / bkgPred
                pull6, ratio6, ratio6Unc = cal_ratioUnc(dataObs[i], bkgPred[i])

                # MC correction
                MC_correction_val  = ratio1
                dataPred_corrected = MC_correction_val * dataPred[i][0]
                _, ratio, ratioUnc = cal_ratioUnc( dataObs[i], (dataPred_corrected, dataPred[i][1]) )
               
                PullUnc.append(pullUnc)
                Pull1.append(pull1)
                Ratio1.append(ratio1)
                Ratio1Unc.append(ratio1Unc)
 
                # -------------------------------------------- 
                # append the information of 7th, 8th Njet bins
                # -------------------------------------------- 
                if i <= 1:
                    # data
                    data_pred.append(dataPred[i][0])
                    data_predUnc.append(dataPred[i][1])
                    data_obs.append(dataObs[i][0])
                    data_obsUnc.append(dataObs[i][1])
    
                    Pull2.append(pull2)
                    Ratio2.append(ratio2)
                    Ratio2Unc.append(ratio2Unc)

                    Pull3.append(pull3)
                    Ratio3.append(ratio3)
                    Ratio3Unc.append(ratio3Unc)
            
                    Pull4.append(pull4)
                    Ratio4.append(ratio4)
                    Ratio4Unc.append(ratio4Unc)

                    Pull5.append(pull5)
                    Ratio5.append(ratio5)
                    Ratio5Unc.append(ratio5Unc)

                    Pull6.append(pull6)
                    Ratio6.append(ratio6)
                    Ratio6Unc.append(ratio6Unc)

                    DataPred_corrected.append(dataPred_corrected)
                    Ratio.append(ratio)
                    RatioUnc.append(ratioUnc)

                else:
                    data_pred.append(999999999999)
                    data_predUnc.append(0)
                    data_obs.append(999999999999)
                    data_obsUnc.append(0)

                    Pull2.append(-999)
                    Ratio2.append(-999)
                    Ratio2Unc.append(0)

                    Pull3.append(-999)
                    Ratio3.append(-999)
                    Ratio3Unc.append(0)

                    Pull4.append(-999)
                    Ratio4.append(-999)
                    Ratio4Unc.append(0)

                    Pull5.append(-999)
                    Ratio5.append(-999)
                    Ratio5Unc.append(0)

                    Pull6.append(-999)
                    Ratio6.append(-999)
                    Ratio6Unc.append(0)

                    DataPred_corrected.append(-999)
                    Ratio.append(-999)
                    RatioUnc.append(0)

            # --------------------------
            # plot Chris's fancy closure
            # --------------------------
            fig = plt.figure(figsize=(15, 15))
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
            #ax1.text(0.05, 0.25,  '%s Metric'%(self.metric),  horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes, fontsize=14, fontweight='bold')
            ax1.text(0.16, 1.065, 'CMS',                      transform=ax.transAxes, fontsize=40, fontweight='bold',   va='top', ha='right')
            ax1.text(0.50, 1.055, 'Preliminary',              transform=ax.transAxes, fontsize=36, fontstyle='italic',  va='top', ha='right')
            ax1.text(1.0,  1.055, '%s (13 TeV)' %(self.year), transform=ax.transAxes, fontsize=34, fontweight='normal', va='top', ha='right')
            ax1.errorbar(x, bkg_pred,  yerr=bkg_predUnc,  label='Predicted MC',   xerr=xUnc, fmt='', capsize=0, color='red',       lw=0, elinewidth=3, marker='o', markeredgecolor='red',       markerfacecolor='red',       markersize=6.0)
            ax1.errorbar(x, bkg_obs,   yerr=bkg_obsUnc,   label='Observed MC',    xerr=xUnc, fmt='', capsize=0, color='black',     lw=0, elinewidth=3, marker='o', markeredgecolor='black',     markerfacecolor='black',     markersize=6.0)
            ax1.errorbar(x, data_pred, yerr=data_predUnc, label='Predicted Data', xerr=xUnc, fmt='', capsize=0, color='palegreen', lw=0, elinewidth=3, marker='o', markeredgecolor='palegreen', markerfacecolor='palegreen', markersize=6.0)
            ax1.errorbar(x, data_obs,  yerr=data_obsUnc,  label='Observed Data',  xerr=xUnc, fmt='', capsize=0, color='green',     lw=0, elinewidth=3, marker='o', markeredgecolor='green',     markerfacecolor='green',     markersize=6.0)
            ax1.set_xticklabels([])
            ax1.set_ylabel('Unweighted Event Counts', fontsize=28)
            ax1.set_yscale('log')
            ax1.set_ylim([5.0, 2e4])

            # pull
            ax2 = fig.add_subplot(14, 1, (7,10))
            fig.subplots_adjust(hspace=0.8, left=0.15, right=0.95)
            ax2.set_xlim([lowerNjets - 0.5, higherNjets + 0.5])
            ax2.errorbar(x, Pull3, yerr=pullUnc, xerr=xUnc, label='PredData vs PredMC', fmt='', capsize=0, color='rosybrown',       lw=0, elinewidth=3, marker='o', markeredgecolor='rosybrown',       markerfacecolor='rosybrown',       markersize=6.0)
            ax2.errorbar(x, Pull4, yerr=pullUnc, xerr=xUnc, label='ObsData vs ObsMC',   fmt='', capsize=0, color='lightgreen',      lw=0, elinewidth=3, marker='o', markeredgecolor='lightgreen',      markerfacecolor='lightgreen',      markersize=6.0)
            ax2.errorbar(x, Pull5, yerr=pullUnc, xerr=xUnc, label='PredData vs ObsMC',  fmt='', capsize=0, color='lightskyblue',    lw=0, elinewidth=3, marker='o', markeredgecolor='lightskyblue',    markerfacecolor='lightskyblue',    markersize=6.0)
            ax2.errorbar(x, Pull6, yerr=pullUnc, xerr=xUnc, label='ObsData vs PredMC',  fmt='', capsize=0, color='mediumslateblue', lw=0, elinewidth=3, marker='o', markeredgecolor='mediumslateblue', markerfacecolor='mediumslateblue', markersize=6.0)
            ax2.axhline(y=0.0, color='black', linestyle='dashed', lw=1.5)
            ax2.grid(axis='y', color='black', linestyle='dashed', which='both')
            ax2.set_xticklabels([])
            ax2.set_ylabel('Pull', fontsize=20)
            ax2.set_ylim([-9.4, 9.4])

            # MC correction
            ax3 = fig.add_subplot(14, 1, (11,14))
            fig.subplots_adjust(hspace=0.8, left=0.15, right=0.95)
            ax3.set_xlim([lowerNjets - 0.5, higherNjets + 0.5])
            ax3.errorbar(x, Ratio2, yerr=Ratio2Unc, xerr=xUnc, label='ObsData / PredData',                      fmt='', capsize=0, color='limegreen',      lw=0, elinewidth=3, marker='o', markeredgecolor='limegreen',      markerfacecolor='limegreen',      markersize=6.0)
            ax3.errorbar(x, Ratio1, yerr=Ratio1Unc, xerr=xUnc, label='MC corr. = ObsMC / PredMC',               fmt='', capsize=0, color='crimson',        lw=0, elinewidth=3, marker='o', markeredgecolor='crimson',        markerfacecolor='crimson',        markersize=6.0)
            ax3.errorbar(x, Ratio,  yerr=RatioUnc,  xerr=xUnc, label='Ratio = ObsData / (MC corr. * PredData)', fmt='', capsize=0, color='cornflowerblue', lw=0, elinewidth=3, marker='o', markeredgecolor='cornflowerblue', markerfacecolor='cornflowerblue', markersize=6.0)
            ax3.axhline(y=1.0, color='black', linestyle='dashed', lw=1.5)
            ax3.grid(axis='y', color='black', linestyle='dashed', which='both')
            #ax3.set_xlabel('Number of jets', fontsize=28)
            ax3.set_ylabel('MC Correction', fontsize=20)
            ax3.set_ylim([0.7, 1.3])

            plt.xticks([int(Njet) for Njet in Njets], fontsize=26)
            plt.xlabel('Number of jets', fontsize=28)    

            ax1.legend(loc='upper right', numpoints=1, frameon=False, fontsize=20)
            ax2.legend(loc='upper right', numpoints=1, frameon=False, fontsize=20)
            ax3.legend(loc='upper right', numpoints=1, frameon=False, fontsize=20)
    
            fig.savefig('%s/%s_Njets_Region_A_PredVsActual_dataVsMC_%s_%s_%s_%s_%s.pdf' % (self.outputDir, self.year, name, self.channel, self.metric, closureTag, bkgTag))
   
            plt.close(fig)

    # ----------------------
    # make all closure plots
    # ----------------------      
    def make_allClosures(self, edgesPerNjets=None, TT_EventsPerNjets=None, NonTT_EventsPerNjets=None, sig_EventsPerNjets=None, data_EventsPerNjets=None, Njets=None, name="", closureTag="", bkgTag=""):
    
        TT_EventsA     = 0.0; TT_EventsB         = 0.0; TT_EventsC    = 0.0; TT_EventsD    = 0.0;
        TT_EventsUncA  = 0.0; TT_EventsUncB      = 0.0; TT_EventsUncC = 0.0; TT_EventsUncD = 0.0;
        TT_Pred_A      = 0.0; TT_PredUnc_A       = 0.0
        TT_EventsNjets = [];  TT_EventsNjetsPred = []
        
        NonTT_EventsA     = 0.0; NonTT_EventsB         = 0.0; NonTT_EventsC    = 0.0; NonTT_EventsD    = 0.0;
        NonTT_EventsUncA  = 0.0; NonTT_EventsUncB      = 0.0; NonTT_EventsUncC = 0.0; NonTT_EventsUncD = 0.0;
        NonTT_Pred_A      = 0.0; NonTT_PredUnc_A       = 0.0
        NonTT_EventsNjets = [];  NonTT_EventsNjetsPred = []
        
        data_EventsA     = 0.0; data_EventsB         = 0.0; data_EventsC    = 0.0; data_EventsD    = 0.0;
        data_EventsUncA  = 0.0; data_EventsUncB      = 0.0; data_EventsUncC = 0.0; data_EventsUncD = 0.0;
        data_Pred_A      = 0.0; data_PredUnc_A       = 0.0
        data_EventsNjets = [];  data_EventsNjetsPred = []

        TT_inData_EventsB     = 0.0; TT_inData_EventsC        = 0.0; TT_inData_EventsD    = 0.0
        TT_inData_EventsUncB  = 0.0; TT_inData_EventsUncC     = 0.0; TT_inData_EventsUncD = 0.0
        TT_inData_Pred_A      = 0.0; TT_inData_PredUnc_A      = 0.0
        TT_NonTT_EventsA_pred = 0.0; TT_NonTT_EventsUncA_pred = 0.0; TT_NonTT_EventsNjetsPred = [] 

        # --------------------------------------
        # get all TT, NonTT, Signal, Data events
        # --------------------------------------
        for Njet in Njets:
    
            if edgesPerNjets[Njet][name][0] == -1.0 or edgesPerNjets[Njet][name][1] == -1.0:
                TT_EventsNjets.append((0.0, 0.0));    TT_EventsNjetsPred.append((0.0, 0.0))
                NonTT_EventsNjets.append((0.0, 0.0)); NonTT_EventsNjetsPred.append((0.0, 0.0))
                data_EventsNjets.append((0.0, 0.0));  data_EventsNjetsPred.append((0.0, 0.0))   
                TT_NonTT_EventsNjetsPred((0.0, 0.0))
 
            else:
                # ------------------------------------------------------
                # get TT and signal injected TT events for usual closure
                # ------------------------------------------------------
                if TT_EventsPerNjets != None:
                    TT_EventsA = TT_EventsPerNjets[Njet][name]["A"][0]; TT_EventsUncA = TT_EventsPerNjets[Njet][name]["A"][1]**2.0
                    TT_EventsB = TT_EventsPerNjets[Njet][name]["B"][0]; TT_EventsUncB = TT_EventsPerNjets[Njet][name]["B"][1]**2.0
                    TT_EventsC = TT_EventsPerNjets[Njet][name]["C"][0]; TT_EventsUncC = TT_EventsPerNjets[Njet][name]["C"][1]**2.0
                    TT_EventsD = TT_EventsPerNjets[Njet][name]["D"][0]; TT_EventsUncD = TT_EventsPerNjets[Njet][name]["D"][1]**2.0

                    if sig_EventsPerNjets != None:
                        TT_EventsA += sig_EventsPerNjets[Njet][name]["A"][0]; TT_EventsUncA += sig_EventsPerNjets[Njet][name]["A"][1]**2.0
                        TT_EventsB += sig_EventsPerNjets[Njet][name]["B"][0]; TT_EventsUncB += sig_EventsPerNjets[Njet][name]["B"][1]**2.0
                        TT_EventsC += sig_EventsPerNjets[Njet][name]["C"][0]; TT_EventsUncC += sig_EventsPerNjets[Njet][name]["C"][1]**2.0
                        TT_EventsD += sig_EventsPerNjets[Njet][name]["D"][0]; TT_EventsUncD += sig_EventsPerNjets[Njet][name]["D"][1]**2.0

                    TT_Pred_A, TT_PredUnc_A = self.cal_simpleClosure_ABCD(TT_EventsA, TT_EventsB, TT_EventsC, TT_EventsD, TT_EventsUncA**0.5, TT_EventsUncB**0.5, TT_EventsUncC**0.5, TT_EventsUncD**0.5)
                    TT_EventsNjets.append((TT_EventsA, TT_EventsUncA**0.5))
                    TT_EventsNjetsPred.append((TT_Pred_A, TT_PredUnc_A))

                # --------------------------------------
                # get the NonTT events for usual closure
                # --------------------------------------
                if  NonTT_EventsPerNjets != None:
                    NonTT_EventsA = NonTT_EventsPerNjets[Njet][name]["A"][0]; NonTT_EventsUncA = NonTT_EventsPerNjets[Njet][name]["A"][1]**2.0
                    NonTT_EventsB = NonTT_EventsPerNjets[Njet][name]["B"][0]; NonTT_EventsUncB = NonTT_EventsPerNjets[Njet][name]["B"][1]**2.0
                    NonTT_EventsC = NonTT_EventsPerNjets[Njet][name]["C"][0]; NonTT_EventsUncC = NonTT_EventsPerNjets[Njet][name]["C"][1]**2.0
                    NonTT_EventsD = NonTT_EventsPerNjets[Njet][name]["D"][0]; NonTT_EventsUncD = NonTT_EventsPerNjets[Njet][name]["D"][1]**2.0                

                    NonTT_Pred_A, NonTT_PredUnc_A = self.cal_simpleClosure_ABCD(NonTT_EventsA, NonTT_EventsB, NonTT_EventsC, NonTT_EventsD, NonTT_EventsUncA**0.5, NonTT_EventsUncB**0.5, NonTT_EventsUncC**0.5, NonTT_EventsUncD**0.5)
                    NonTT_EventsNjets.append((NonTT_EventsA, NonTT_EventsUncA**0.5))
                    NonTT_EventsNjetsPred.append((NonTT_Pred_A, NonTT_PredUnc_A))

                # ---------------------------------------------
                # get data events for usual and data-MC closure
                # ---------------------------------------------
                if data_EventsPerNjets != None:
                    data_EventsA = data_EventsPerNjets[Njet][name]["A"][0]; data_EventsUncA = data_EventsPerNjets[Njet][name]["A"][1]**2.0
                    data_EventsB = data_EventsPerNjets[Njet][name]["B"][0]; data_EventsUncB = data_EventsPerNjets[Njet][name]["B"][1]**2.0
                    data_EventsC = data_EventsPerNjets[Njet][name]["C"][0]; data_EventsUncC = data_EventsPerNjets[Njet][name]["C"][1]**2.0
                    data_EventsD = data_EventsPerNjets[Njet][name]["D"][0]; data_EventsUncD = data_EventsPerNjets[Njet][name]["D"][1]**2.0

                    data_Pred_A, data_PredUnc_A = self.cal_simpleClosure_ABCD(data_EventsA, data_EventsB, data_EventsC, data_EventsD, data_EventsUncA**0.5, data_EventsUncB**0.5, data_EventsUncC**0.5, data_EventsUncD**0.5)
                    data_EventsNjets.append((data_EventsA, data_EventsUncA**0.5))
                    data_EventsNjetsPred.append((data_Pred_A, data_PredUnc_A))                

                # ---------------------------------------------------------
                # get the TT events in data, by subtracting NonTT from data
                #   and get (TT + NonTT) as all bkgPred
                #   -- bkgObs  = data_EventsNjets
                #   -- bkgPred = TT_NonTT_EventsNjetsPred  
                # ---------------------------------------------------------
                if TT_EventsPerNjets != None and NonTT_EventsPerNjets != None and data_EventsPerNjets != None:
                    TT_inData_EventsB = data_EventsB - NonTT_EventsB; TT_inData_EventsUncB = data_EventsUncB + NonTT_EventsUncB
                    TT_inData_EventsC = data_EventsC - NonTT_EventsC; TT_inData_EventsUncC = data_EventsUncC + NonTT_EventsUncC
                    TT_inData_EventsD = data_EventsD - NonTT_EventsD; TT_inData_EventsUncD = data_EventsUncD + NonTT_EventsUncD
                    
                    TT_inData_Pred_A, TT_inData_PredUnc_A = self.cal_simpleClosure_ABCD(TT_EventsA, TT_inData_EventsB, TT_inData_EventsC, TT_inData_EventsD, TT_EventsUncA**0.5, TT_inData_EventsUncB**0.5, TT_inData_EventsUncC**0.5, TT_inData_EventsUncD**0.5)
                    TT_NonTT_EventsA_pred    = TT_inData_Pred_A    + NonTT_EventsA
                    TT_NonTT_EventsUncA_pred = TT_inData_PredUnc_A + NonTT_EventsUncA
                    TT_NonTT_EventsNjetsPred.append((TT_NonTT_EventsA_pred, TT_NonTT_EventsUncA_pred**0.5))

        # ------------------
        # make closure plots
        # ------------------ 
        # usual closure
        if TT_EventsPerNjets != None and NonTT_EventsPerNjets == None and data_EventsPerNjets == None: 
            self.plot_ClosureNjets(np.array(TT_EventsNjets), np.array(TT_EventsNjetsPred), Njets, name, closureTag, bkgTag)
        
        if NonTT_EventsPerNjets != None and TT_EventsPerNjets == None and  data_EventsPerNjets == None:
            self.plot_ClosureNjets(np.array(NonTT_EventsNjets), np.array(NonTT_EventsNjetsPred), Njets, name, closureTag, bkgTag)

        if TT_EventsPerNjets != None and NonTT_EventsPerNjets != None and data_EventsPerNjets != None:
            self.plot_ClosureNjets(np.array(data_EventsNjets), np.array(TT_NonTT_EventsNjetsPred), Njets, name, closureTag, bkgTag, isBlind=True) 

        # data-MC closure
        if TT_EventsPerNjets != None and data_EventsPerNjets != None and NonTT_EventsPerNjets == None: 
            self.plot_dataVsMC_ClosureNjets(np.array(TT_EventsNjets), np.array(TT_EventsNjetsPred), np.array(data_EventsNjets), np.array(data_EventsNjetsPred), Njets, name, closureTag, bkgTag) # for TT


    # ----------------------------------------------------------------
    # plot whichever variable as a function of the choice of bin edges
    # ----------------------------------------------------------------
    def plot_Var_vsDisc1Disc2(self, var, edges, c1, c2, minEdge, maxEdge, binWidth, cmin, cmax, vmin, vmax, Njets = -1, name = "", variable = ""):

        nBins = math.ceil((1.0 + binWidth)/binWidth)

        fig = plt.figure() 
        fig.subplots_adjust(top = 0.94)
        plt.hist2d(edges[:,0], edges[:,1], bins=[nBins, nBins], range=[[-binWidth/2.0, 1+binWidth/2.0], [-binWidth/2.0, 1+binWidth/2.0]], cmap=plt.cm.jet, weights=var, cmin=cmin, cmax=cmax, vmin=vmin, vmax=vmax)
        plt.colorbar()
        ax = plt.gca()
        ax.set_xlabel("Disc. 1 Bin Edge")
        ax.set_ylabel("Disc. 2 Bin Edge")
        ax.text(0.16, 1.065, 'CMS',                     transform=ax.transAxes, fontsize=20, fontweight='bold',   va='top', ha='right')
        ax.text(0.50, 1.055, 'Preliminary',             transform=ax.transAxes, fontsize=16, fontstyle='italic',  va='top', ha='right')
        ax.text(1.0,  1.055, '%s (13 TeV)' % self.year, transform=ax.transAxes, fontsize=14, fontweight='normal', va='top', ha='right')

        l1 = ml.Line2D([c1, c1], [0.0, 1.0], color="black", linewidth=2, linestyle="dashed"); l2 = ml.Line2D([0.0, 1.0], [c2, c2], color="black", linewidth=2, linestyle="dashed")
        ax.add_line(l1); ax.add_line(l2)
        fig.tight_layout()

        fig.savefig(self.outputDir+"/%s_%s_vs_Disc1Disc2_Njets%s_%s_%s_%s.pdf"%(self.year, variable, Njets, name, self.channel, self.metric), dpi=fig.dpi)

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
        ax.set_xlim(0.0, 20.0)
        ax.set_ylim(0.0, 0.5)
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

        fig.savefig('%s/%s_InvSign_vs_ClosureErr_Njets%s_%s_%s_%s.pdf' % (self.outputDir, self.year, Njets, name, self.channel, self.metric), dpi=fig.dpi)

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
        
        l3 = None; l4 = None
        if name == "Val_SubDivD":
            l3 = ml.Line2D([float(allRegionsEdges["Val_subDivD"][0]), float(allRegionsEdges["Val_subDivD"][0])], [0.0, 1.0],   color=col1, linewidth=4, linestyle='solid')
            l4 = ml.Line2D([0.0, 1.0], [float(allRegionsEdges["Val_subDivD"][1]), float(allRegionsEdges["Val_subDivD"][1])], color=col2, linewidth=4, linestyle='solid')
        else:
            l3 = ml.Line2D([float(allRegionsEdges["Val_bdEF"][0]), float(allRegionsEdges["Val_bdEF"][0])], [0.0, 1.0],   color=col1, linewidth=4, linestyle='solid')
            l4 = ml.Line2D([0.0, 1.0], [float(allRegionsEdges["Val_cdiGH"][1]), float(allRegionsEdges["Val_cdiGH"][1])], color=col2, linewidth=4, linestyle='solid')
           
        ax.add_line(l1)
        ax.add_line(l2)
        ax.add_line(l3)
        ax.add_line(l4)

        fig.tight_layout()
        fig.savefig('%s/%s_Disc1VsDisc2_%s_Njets%s_%s_%s_%s.pdf' % (self.outputDir, self.year, tag, Njets, name, self.channel, self.metric), dpi=fig.dpi)

        plt.close(fig)

    # ----------------------------------------------------------------------
    # plot variable vs disc as 1D
    #   -- Closure, Significance, weightedEventCounts vs disc1, disc2 slices
    # ----------------------------------------------------------------------
    def plot_VarVsDisc(self, var, edges, edgeWidth, ylim = -1.0, ylabel = "", tag = "", disc = -1, Njets = -1, name=""):

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

        ax.set_ylabel(ylabel)
        ax.set_xlabel("Disc. %d Value"%(3-disc))
        plt.legend(loc='best', numpoints=1)

        ax.text(0.12, 1.05, 'CMS',                     transform=ax.transAxes, fontsize=14, fontweight='bold',   va='top', ha='right')
        ax.text(0.33, 1.04, 'Preliminary',             transform=ax.transAxes, fontsize=10, fontstyle='italic',  va='top', ha='right')
        ax.text(0.99, 1.04, '%s (13 TeV)' % self.year, transform=ax.transAxes, fontsize=10, fontweight='normal', va='top', ha='right') 

        fig.tight_layout()
        fig.savefig('%s/%s_%s_Slices_Disc%d_Njets%s_%s_%s_%s.pdf' % (self.outputDir, self.year, tag, disc, Njets, name, self.channel, self.metric), dpi=fig.dpi)

        plt.close(fig)