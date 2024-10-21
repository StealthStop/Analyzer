import math
import numpy as np
from sklearn.linear_model import LinearRegression

import matplotlib as mpl
mpl.use('Agg')

import matplotlib.lines as ml
import matplotlib.pyplot as plt

from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import matplotlib.ticker as ticker

from collections import OrderedDict
import mplhep as hep
plt.style.use([hep.style.ROOT,hep.style.CMS]) # For now ROOT defaults to CMS
plt.style.use({'legend.frameon':False,'legend.fontsize':16,'legend.edgecolor':'black'})

class Common_Calculations_Plotters:

    def __init__(self, outputDir, year, model, mass, channel, cmsLabel):
        self.year            = year
        self.model           = model
        self.mass            = mass
        self.channel         = channel
        self.cmsLabel        = cmsLabel
        self.cmap            = plt.cm.jet
        self.cmap.set_under("w", 1)

        self.outputDir = outputDir

    def addCMSlabel(self, ax, fontsize=None):

        if fontsize == None:
            ax.text(0.0,  1.003, 'CMS',                     transform=ax.transAxes, fontsize=16, fontweight='bold',   va='bottom', ha='left')
            ax.text(0.15, 1.010, '%s'%(self.cmsLabel),      transform=ax.transAxes, fontsize=11, fontstyle='italic',  va='bottom', ha='left')
            ax.text(1.0,  1.010, '%s (13 TeV)'%(self.year), transform=ax.transAxes, fontsize=11, fontweight='normal', va='bottom', ha='right')

        else:
            ax.text(0.0,  1.003, 'CMS',                     transform=ax.transAxes, fontsize=16*fontsize, fontweight='bold',   va='bottom', ha='left')
            #ax.text(0.08, 1.010, '%s'%(self.cmsLabel),      transform=ax.transAxes, fontsize=11*fontsize, fontstyle='italic',  va='bottom', ha='left')
            ax.text(1.0,  1.010, 'Run 2 (13 TeV)', transform=ax.transAxes, fontsize=11*fontsize, fontweight='normal', va='bottom', ha='right')
    
        return ax

    # ----------------------
    # calculate all closures
    # ----------------------
    def cal_simpleClosure_ABCD(self, nEvents_A, nEvents_B, nEvents_C, nEvents_D, nEventsErr_A, nEventsErr_B, nEventsErr_C, nEventsErr_D):

        if nEvents_D == 0.0:
            return -999.0, 0.0

        nPred_A    = (nEvents_B * nEvents_C) / nEvents_D
        nPred_Aunc = ((nEvents_C * nEventsErr_B / nEvents_D)**2.0 + (nEventsErr_C * nEvents_B / nEvents_D)**2.0 + (nEvents_C * nEvents_B * nEventsErr_D / nEvents_D**2.0)**2.0)**0.5

        return nPred_A, nPred_Aunc

    # -----------------
    # plot all closures
    # -----------------
    def plot_ClosureNjets(self, bkgObs, bkgPred, Njets, name = '', closureTag = '', bkgTag = '', valColor=None, isBlind=False):
        
        x         = []; xUnc         = []
        obs       = []; obsUnc       = [] 
        pred      = []; predUnc      = []
        abcdError = []; abcdErrorUnc = []
        abcdPull  = []; abcdPullUnc  = []
        pullDenom = []; zeros        = []
        closureCorrections = []; closureCorrectionsUnc = []
        
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
    
                x.append(float((Njets[i]).replace("incl", "")))  
                xUnc.append(0.5)

                pull            = (bkgPred[i][0] - bkgObs[i][0]) / ( bkgPred[i][1]**2 + bkgObs[i][1]**2 )**0.5
                pullUnc         = 1.0
                closureError    = 1.0 - ( bkgPred[i][0] / bkgObs[i][0] )
                closureErrorUnc = ((bkgPred[i][1] / bkgObs[i][0])**2.0 + (bkgObs[i][1] * bkgPred[i][0] / bkgObs[i][0]**2.0)**2.0)**0.5
                closureCorrection = ( bkgPred[i][0] / bkgObs[i][0] )
                closureCorrectionUnc = closureErrorUnc
                pullDenominator = ( bkgPred[i][1]**2 + bkgObs[i][1]**2 )**0.5 / bkgPred[i][0]
                zero            = 0.0

                # for the data closure / blind the 9-10-11 njet bins
                if isBlind and i > 1:
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
                    closureCorrections.append(999)
                    closureCorrectionsUnc.append(999)

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
                    closureCorrections.append(closureCorrection)
                    closureCorrectionsUnc.append(closureCorrectionUnc)
                    totalChi2 += pull ** 2.0
                    ndof      += 1                    

        closureCorrections = [float(1./cc) for cc in closureCorrections]

        print(closureCorrections)
        print(closureCorrectionsUnc)

        # ------------------
        # plot usual closure
        # ------------------
        fig = plt.figure(figsize=(5, 5))
        ax  = fig.add_subplot(111)
        fig.subplots_adjust(hspace=0, left=0.15, right=0.95, top=0.95)
    
        lowerNjets  = float(Njets[0])
        higherNjets = float( (Njets[(-1)]).replace("incl","") )
    
        ax.spines['top'].set_color('none')
        ax.spines['bottom'].set_color('none')
        ax.spines['left'].set_color('none')
        ax.spines['right'].set_color('none')
        ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
    
        # unweighted event counts
        ax1 = fig.add_subplot(4, 1, (1, 3))  
        fig.subplots_adjust(left=0.15, right=0.95)
        ax1.set_yscale('log')
        ax1.set_xlim([lowerNjets - 0.5, higherNjets + 0.5])
        
        ax1 = self.addCMSlabel(ax1)

        # put model, channel labels
        md = ""
        if self.model == "SYY":
            md = "Stealth SYY"
        else:
            md = self.model

        ch = ""
        if self.channel == "0l":
            ch = "0L"
        elif self.channel == "1l":
            ch = "1L"
        elif self.channel == "2l":
            ch = "2L"

        textLabel = '\n'.join(( "%s"%(name), "%s"%(md), "%s"%(ch) ))
        ax1.text(0.66, 0.80, textLabel, transform=ax.transAxes, color="black",  fontsize=9,  fontweight='normal', va='center', ha='left' )
        #ax1.text(0.5, 0.95, name,      transform=ax.transAxes, color=valColor, fontsize=11, fontweight='bold',   va='center', ha='center') 

        ax1.set_ylabel('Num. Events', fontsize=11)
        ax1.errorbar(x, pred, yerr=predUnc, label=r'$N_{A,Pred.}$', xerr=xUnc, fmt='', capsize=0, color='red',   lw=0, elinewidth=2, marker='o', markeredgecolor='red',   markerfacecolor='red',   markersize=5.0)
        ax1.errorbar(x, obs,  yerr=obsUnc,  label=r'$N_{A}$',  xerr=xUnc, fmt='', capsize=0, color='black', lw=0, elinewidth=2, marker='o', markeredgecolor='black', markerfacecolor='black', markersize=5.0)

        # closure   
        #ax2 = fig.add_subplot(4, 1, 3)
        #ax2.set_xlim([lowerNjets - 0.5, higherNjets + 0.5])
        #ax2.errorbar(x, abcdError, yerr=abcdErrorUnc, xerr=xUnc, fmt='', capsize=0, color='blue', lw=0, elinewidth=2, marker='o', markeredgecolor='blue', markerfacecolor='blue', markersize=5.0)
        #ax2.axhline(y=0.0, color='black', linestyle='dashed', lw=1.5)
        #ax2.grid(axis='y', color='black', linestyle='dashed', which='both')
        #ax2.set_ylabel('Non-Closure', fontsize=12)
        #ax2.set_ylim([-0.59, 0.59])   

        # Closure Correction Value
        ax3 = fig.add_subplot(4, 1, 4)
        ax3.set_xlim([lowerNjets - 0.5, higherNjets + 0.5])
        ax3.errorbar(x, closureCorrections, yerr=closureCorrectionsUnc, xerr=xUnc, fmt='', capsize=0, color='purple', lw=0, elinewidth=2, marker='o', markeredgecolor='purple', markerfacecolor='purple', markersize=5.0)
        ax3.axhline(y=1.0, color='black', linestyle='dashed', lw=1.5)
        ax3.grid(axis='y', color='black', linestyle='dashed', which='both')
        ax3.set_xlabel('Number of jets', fontsize=12)
        ax3.set_ylabel(r'$\kappa$', fontsize=12) 
        ax3.set_ylim([-0.2, 2.2])

        # pull 
        #ax3 = fig.add_subplot(4, 1, 4)
        #ax3.set_xlim([lowerNjets - 0.5, higherNjets + 0.5])
        #ax3.errorbar(x, abcdPull, yerr=abcdPullUnc, xerr=xUnc, fmt='', capsize=0, color='purple', lw=0, elinewidth=2, marker='o', markeredgecolor='purple', markerfacecolor='purple', markersize=5.0)
        #ax3.axhline(y=0.0, color='black', linestyle='dashed', lw=1.5)
        #ax3.grid(axis='y', color='black', linestyle='dashed', which='both')
        #ax3.set_xlabel('Number of jets', fontsize=12)
        #ax3.set_ylabel('Pull',           fontsize=12) 
        #ax3.set_ylim([-5.9, 5.9])
 
        #plt.show()
        plt.xticks([int(Njet.replace("incl","")) for Njet in Njets])
    
        ax1.legend(loc='best', numpoints=1, frameon=False)
    
        fig.savefig('%s/%s_Njets_Region_A_PredVsActual_%s_%s_%s_%s.pdf' % (self.outputDir, self.year, name, closureTag, bkgTag, self.channel))
   
        plt.close(fig)

    # ----------------------
    # make all closure plots
    # ----------------------      
    def make_allClosures(self, edgesPerNjets=None, TT_EventsPerNjets=None, NonTT_EventsPerNjets=None, sig_EventsPerNjets=None, data_EventsPerNjets=None, Njets=None, name="", closureTag="", bkgTag="", valColor=None):
    
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

                # --------------------------------------
                # get data events for usual data closure
                # --------------------------------------
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
                # (predicted Data = predicted TT in data + observed NonTT)
                # ---------------------------------------------------------
                if TT_EventsPerNjets != None and NonTT_EventsPerNjets != None and data_EventsPerNjets != None:
                    TT_inData_EventsB = data_EventsB - NonTT_EventsB; TT_inData_EventsUncB = (data_EventsUncB)**2.0 + (NonTT_EventsUncB)**2.0
                    TT_inData_EventsC = data_EventsC - NonTT_EventsC; TT_inData_EventsUncC = (data_EventsUncC)**2.0 + (NonTT_EventsUncC)**2.0
                    TT_inData_EventsD = data_EventsD - NonTT_EventsD; TT_inData_EventsUncD = (data_EventsUncD)**2.0 + (NonTT_EventsUncD)**2.0
                    
                    TT_inData_Pred_A, TT_inData_PredUnc_A = self.cal_simpleClosure_ABCD(TT_EventsA, TT_inData_EventsB, TT_inData_EventsC, TT_inData_EventsD, TT_EventsUncA**0.5, TT_inData_EventsUncB**0.5, TT_inData_EventsUncC**0.5, TT_inData_EventsUncD**0.5)
                    TT_NonTT_EventsA_pred    = TT_inData_Pred_A    + NonTT_EventsA
                    TT_NonTT_EventsUncA_pred = TT_inData_PredUnc_A + NonTT_EventsUncA
                    TT_NonTT_EventsNjetsPred.append((TT_NonTT_EventsA_pred, TT_NonTT_EventsUncA_pred**0.5))

        # ------------------
        # make closure plots
        # ------------------ 
        # usual closure
        if TT_EventsPerNjets != None and NonTT_EventsPerNjets == None and data_EventsPerNjets == None: 
            self.plot_ClosureNjets(np.array(TT_EventsNjets), np.array(TT_EventsNjetsPred), Njets, name, closureTag, bkgTag, valColor)
        
        # MC correction factor related closure
        if TT_EventsPerNjets != None and data_EventsPerNjets != None and NonTT_EventsPerNjets != None: 
            self.plot_MCcorr_ClosureNjets(np.array(TT_EventsNjets), np.array(TT_EventsNjetsPred), np.array(data_EventsNjets), np.array(TT_NonTT_EventsNjetsPred), Njets, name, closureTag, bkgTag, valColor)



    # ----------------------------------------------------------------
    # plot whichever variable as a function of the choice of bin edges
    #   -- call it inside DoubleDisCo_BinEdges module
    #   -- make significance, non-closure, pull, sigFrac plots
    # ----------------------------------------------------------------
    def plot_Var_vsDisc1Disc2(self, var, edges, c1, c2, minEdge, maxEdge, binWidth, cmin, cmax, vmin, vmax, Njets = -1, name = "", variable = ""):

        nBins = math.ceil((1.0 + binWidth)/binWidth)

        fig = plt.figure(figsize=(6, 5)) 
        plt.hist2d(edges[:,0], edges[:,1], bins=[nBins, nBins], range=[[-binWidth/2.0, 1+binWidth/2.0], [-binWidth/2.0, 1+binWidth/2.0]], cmap=self.cmap, weights=var, cmin=cmin, cmax=cmax, vmin=vmin, vmax=vmax)
        plt.colorbar()
        ax = plt.gca()
        ax.set_xlabel("Disc. 1 Bin Edge", fontsize=14)
        ax.set_ylabel("Disc. 2 Bin Edge", fontsize=14)

        ax = self.addCMSlabel(ax)

        # put model, channel, njet labels
        md = ""
        if self.model == "SYY":
            md = "Stealth SYY"
        else:
            md = self.model

        ch = ""
        if self.channel == "0l":
            ch = "0L"
        elif self.channel == "1l":
            ch = "1L"
        elif self.channel == "2l":
            ch = "2L"

        nj = ""
        if "incl" in Njets:
            nj = "$N_{jets} \geq$ %s"%(Njets.replace("incl",""))
        else:
            nj = "$N_{jets}$ = %s"%(Njets)

        vr = ""
        if variable == "Sign_includingNonClosure":
            vr = "Significance"
        elif variable == "NonClosure":
            vr = "Non-Closure"
        else:
            vr = "%s"%(variable)

        textLabel = '\n'.join(( "%s in %s"%(vr,name), "%s"%(md), "%s"%(ch), nj ))

        if name == "ABCD":
            ax.text(0.05, 0.10, textLabel, transform=ax.transAxes, color="cadetblue", fontsize=9, fontweight='bold',  va='center', ha='left')
        if name == "Val_BD":
            ax.text(0.95, 0.10, textLabel, transform=ax.transAxes, color="cadetblue", fontsize=9, fontweight='bold',  va='center', ha='right')
        if name == "Val_CD":
            ax.text(0.05, 0.90, textLabel, transform=ax.transAxes, color="cadetblue", fontsize=9, fontweight='bold',  va='center', ha='left')
        if name == "Val_D":
            ax.text(0.95, 0.90, textLabel, transform=ax.transAxes, color="cadetblue", fontsize=9, fontweight='bold',  va='center', ha='right')

        #
        l1 = ml.Line2D([c1, c1], [0.0, 1.0], color="cadetblue", linewidth=2, linestyle="dashed"); l2 = ml.Line2D([0.0, 1.0], [c2, c2], color="cadetblue", linewidth=2, linestyle="dashed")
        ax.add_line(l1); 
        ax.add_line(l2)

        fig.tight_layout()
        fig.subplots_adjust(top=0.95, right=0.99)

        fig.savefig(self.outputDir+"/%s_%s_vs_Disc1Disc2_Njets%s_%s_%s.pdf"%(self.year, variable, Njets, name, self.channel), dpi=fig.dpi)

        plt.close(fig)

    # -------------------------------------------------
    # plot Disc1 vs Disc2 in each Region with all edges
    # -------------------------------------------------
    def plot_Disc1VsDisc2(self, hist, allRegionsEdges, Njets = -1, tag = "", name=""):

        nBins   = hist.GetXaxis().GetNbins()
        nEvents = []; disc1Edges = []; disc2Edges = []

        for x in range(0, hist.GetXaxis().GetNbins()+1):
        
            for y in range(0, hist.GetYaxis().GetNbins()+1):
        
                nEvents.append(hist.GetBinContent(x,y))
                disc1Edges.append(hist.GetXaxis().GetBinCenter(x))
                disc2Edges.append(hist.GetYaxis().GetBinCenter(y))

        fig = plt.figure(figsize=(6, 5))
        plt.hist2d(disc1Edges, disc2Edges, bins=[nBins, nBins], range=[[0.0, 1.0], [0.0, 1.0]], cmap=plt.cm.jet, weights=nEvents, cmin=10e-10)
        plt.colorbar()
        ax = plt.gca()
        ax.set_xlabel('Disc. 1', fontsize=14)
        ax.set_ylabel('Disc. 2', fontsize=14)

        ax = self.addCMSlabel(ax)

        fig.tight_layout()
        fig.subplots_adjust(top=0.95, right=0.99)

        fig.savefig('%s/%s_Disc1VsDisc2_%s_Njets%s_%s_%s.pdf' % (self.outputDir, self.year, tag, Njets, name, self.channel), dpi=fig.dpi)

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

        ax.set_xlabel("Disc. %d Value"%(3-disc), fontsize=14)
        ax.set_ylabel(ylabel, fontsize=14)
        plt.legend(loc='best', numpoints=1)

        ax = self.addCMSlabel(ax)

        fig.tight_layout()
        fig.subplots_adjust(top=0.95, right=0.95)

        fig.savefig('%s/%s_%s_Slices_Disc%d_Njets%s_%s_%s.pdf' % (self.outputDir, self.year, tag, disc, Njets, name, self.channel), dpi=fig.dpi)

        plt.close(fig)

    # -----------------------------------------------
    # plot variable as a function of per boundary
    #   -- as part of closure correction factor study
    # -----------------------------------------------
    def plot_VarVsBoundary(self, var, xWidth, yMin = None, yMax = None, lineY = None, ylabel = "", tag = "", Njets = -1, color=None):

        regions = var.keys()

        nice_region_labels = {
            "Val_BD": "Val. I",
            "Val_CD": "Val. II",
            "Val_D": "Val. III",
        }

        x    = {region : [] for region in regions}
        y    = {region : [] for region in regions}
        yUnc = {region : [] for region in regions} 
        xUnc = {region : [] for region in regions}

        fig = plt.figure(figsize=(5, 5))
        ax = plt.gca()

        regions.reverse()

        for region in regions:
            for boundary, value in var[region].items():
                x[region].append(boundary)
                y[region].append(value[0])
                yUnc[region].append(value[1])
                xUnc[region].append(xWidth)

            ax.errorbar(x[region], y[region], yerr=yUnc[region], label="%s"%(nice_region_labels[region]), xerr=xUnc[region], fmt='', capsize=0, color=color[region], lw=0, elinewidth=2, marker="o", markeredgecolor=color[region], markerfacecolor=color[region])

            if yMin != None and yMax != None:
                ax.set_ylim((yMin, yMax))

        ax.set_xlabel("Boundary Value", fontsize=14)
        ax.set_ylabel(ylabel, fontsize=14)
        plt.legend(loc='upper right', numpoints=1)

        ax = self.addCMSlabel(ax)

        # put model, channel, njet labels
        md = ""
        if self.model == "SYY":
            md = "Stealth SYY"
        else:
            md = self.model

        ch = ""
        if self.channel == "0l":
            ch = "0L"
        elif self.channel == "1l":
            ch = "1L"
        elif self.channel == "2l":
            ch = "2L"

        nj = ""
        if "incl" in Njets:
            nj = "$N_{jets} \geq$ %s"%(Njets.replace("incl",""))
        else:
            nj = "$N_{jets}$ = %s"%(Njets)

        textLabel = '\n'.join(( "%s"%(md), "%s"%(ch), nj ))
        ax.text(0.03, 0.9, textLabel, transform=ax.transAxes, color="black", fontsize=11, fontweight='normal',  va='center', ha='left')

        #
        if lineY != None:
            l1 = ml.Line2D([0.0, 1.05], [lineY, lineY], color="black", linewidth=2, linestyle="dashed")
            ax.add_line(l1); 

        fig.tight_layout()
        fig.subplots_adjust(top=0.95, right=0.95)
        fig.savefig('%s/%s_%s_Njets%s_%s.pdf' % (self.outputDir, self.year, tag, Njets, self.channel), dpi=fig.dpi)

        plt.close(fig)

    # --------------------------------------------------------------------
    # Same plotting function as above but with a ratio between Data and MC
    # --------------------------------------------------------------------
    def plot_VarVsBoundaryRatio(self, var1, var2, xWidth, yMin = None, yMax = None, lineY = None, ylabel = "", tag = "", Njets = -1, color=None, extraSys=False):

        regions = var1.keys()

        label_map = {
            "Val_BD": "Val I",
            "Val_CD": "Val II",
            "Val_D": "Val III",
            }

        x    = {region : [] for region in regions}
        y    = {region : [] for region in regions}
        yUnc = {region : [] for region in regions} 
        xUnc = {region : [] for region in regions}

        yMin = -0.15
        yMax = 0.15

        fig, ax = plt.subplots(2,1, figsize=(12,10), gridspec_kw={"height_ratios": [3,1]})
        #ax = plt.gca()

        for idx, region in enumerate(regions):
            for boundary, value in var1[region].items():
                x[region].append(boundary + 0.05 * (2.*(idx-1.5) + 0.5)/8.)
                y[region].append(value[0])
                yUnc[region].append(value[1])
                xUnc[region].append(xWidth)

                yMin = min(yMin, min(y[region])*1.4)
                yMax = max(yMax, max(y[region])*1.4)

                if abs(yMin) > yMax:
                    yMax = -yMin
                else:
                    yMin = -yMax   

                #ax[0].errorbar(x[region], y[region], yerr=yUnc[region], label="%s"%(region), xerr=xUnc[region], fmt='', capsize=0, color=color[region], lw=0, elinewidth=2, marker="o", markeredgecolor=color[region], markerfacecolor=color[region])
            if not extraSys:
                ax[0].errorbar(x[region], y[region], yerr=yUnc[region], label="%s (Sim.)"%(label_map[region]), fmt='', color=color[region], lw=0, elinewidth=2, markersize=10, marker=".", markeredgecolor=color[region], markerfacecolor=color[region])

            if yMin != None and yMax != None:
                ax[0].set_ylim((yMin, yMax))
                ax[0].set_xlim((0.375, 1.025))
                #ax[0].xaxis.set_ticks(np.arange(0.4, 1.05, 0.05))
                #ax[0].xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))

                #ax[0].tick_params(axis='x', labelsize=16)
                ax[0].tick_params(axis='y', labelsize=16)
                ax[1].set_xlim((0.375, 1.025))
                ax[1].tick_params(axis='x', labelsize=16)
                ax[1].tick_params(axis='y', labelsize=16)

        print(x)
        print(y)

        x2    = {region : [] for region in regions}
        y2    = {region : [] for region in regions}
        yUnc2 = {region : [] for region in regions} 
        xUnc2 = {region : [] for region in regions}

        ratio = {region : [] for region in regions}
        ratioUnc = {region : [] for region in regions}

        temp_sys = {region : [] for region in regions}

        yMinRatio = -0.15
        yMaxRatio = 0.15

        for idx, region in enumerate(regions):
            for i, (boundary, value) in enumerate(var2[region].items()):
                x2[region].append(boundary + 0.05 *(2. * (idx-1.5) + 1.5)/8.)
                y2[region].append(value[0])
                yUnc2[region].append(value[1])
                xUnc2[region].append(xWidth)

                if y[region][-1] != -999.0 and y2[region][-1] != -999.0:
                    ratio[region].append(y2[region][i] - y[region][i])
                    ratioUnc[region].append(math.sqrt(yUnc2[region][i]**2 + yUnc[region][i]**2))
                    temp_sys[region].append(abs(1-y2[region][i]) / abs(1-y[region][i]))
                else:
                    ratio[region].append(0.0)
                    ratioUnc[region].append(0.0)
                    temp_sys[region].append(1.0)
             
                yMinRatio = min(yMinRatio, min(ratio[region])*1.4)
                yMaxRatio = max(yMaxRatio, max(ratio[region])*1.4)

                if abs(yMinRatio) > yMaxRatio:
                    yMaxRatio = -yMinRatio
                else:
                    yMinRatio = -yMaxRatio   
                

                #ax[0].errorbar(x[region], y[region], yerr=yUnc[region], label="%s"%(region), xerr=xUnc[region], fmt='', capsize=0, color=color[region], lw=0, elinewidth=2, marker="^", markeredgecolor="black", markerfacecolor=color[region])
            ax[0].errorbar(x2[region], y2[region], yerr=yUnc2[region], label="%s (Data)"%(label_map[region]), fmt='', color=color[region], lw=0, markersize=10, elinewidth=2, marker="^", markeredgecolor=color[region], markerfacecolor=color[region])

            if yMin != None and yMax != None:
                ax[0].set_ylim((yMin, yMax))
                #ax[0].xaxis.set_ticks(np.arange(0.4, 1.05, 0.05))
                #ax[0].xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
                ax[0].tick_params(
                    axis='x',          # changes apply to the x-axis
                    which='both',      # both major and minor ticks are affected
                    bottom=False,      # ticks along the bottom edge are off
                    top=False,         # ticks along the top edge are off
                    labelbottom=False)
                ax[1].set_ylim((yMinRatio, yMaxRatio))
                ax[1].xaxis.set_ticks(np.arange(0.4, 1.05, 0.05))
                ax[1].xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))

    
            if extraSys:
                dataSys = 0.0
                for idx, region in enumerate(regions):
                    dataSys = max(max(abs(1-max(temp_sys[region])), abs(1-min(temp_sys[region]))), dataSys)

                    if dataSys == 999.0:
                        dataSys = 0.0

                #for idx, region in enumerate(regions):
                #    yUnc[region] = [math.sqrt(val**2 + dataSys**2) for val in yUnc[region]]
                #    ratioUnc[region] = [math.sqrt(val**2 + dataSys**2) for val in ratioUnc[region]]

                    ax[0].errorbar(x[region], y[region], yerr=yUnc[region], label="%s (Sim.)"%(label_map[region]), fmt='', color=color[region], lw=0, elinewidth=2, markersize=10, marker=".", markeredgecolor=color[region], markerfacecolor=color[region])

            for idx, region in enumerate(regions):
                ratio[region] = [x if x != 0.0 else -999.0 for x in ratio[region]]

        for xval in np.arange(0.375,1.075,0.05):
            ax[0].axvline(x=xval, linestyle="dotted")
            ax[1].axvline(x=xval, linestyle="dotted")

        ax[1].set_xlabel("Boundary Value", fontsize=20)
        ax[1].set_ylabel("Data - Sim.", fontsize=20)
        ax[0].set_ylabel("Non-Closure", fontsize=20)
        ax[0].legend(loc='lower right', numpoints=1, ncol=2, fontsize=20)

        ax[0] = self.addCMSlabel(ax[0], fontsize = 1.4)

        for region in regions:
            ax[1].errorbar(x2[region], ratio[region], yerr=ratioUnc[region], label="%s"%(region), fmt='', color=color[region], lw=0, markersize=10, elinewidth=2, marker="^", markeredgecolor=color[region], markerfacecolor=color[region])
            if extraSys:
                #ax[1].plot(x2[region], ratio[region])
                x_fill = np.arange(0.1, 1.5, 0.05)
                ax[1].fill_between(x_fill, [0-dataSys for x in x_fill], [0+dataSys for x in x_fill], alpha=0.2, color="grey")

        # put model, channel, njet labels
        md = ""
        if self.model == "SYY":
            md = "Stealth SYY"
        else:
            md = self.model

        ch = ""
        if self.channel == "0l":
            ch = "0L"
        elif self.channel == "1l":
            ch = "1L"
        elif self.channel == "2l":
            ch = "2L"

        nj = ""
        if "incl" in Njets:
            nj = "$N_{jets} \geq$ %s"%(Njets.replace("incl",""))
        else:
            nj = "$N_{jets}$ = %s"%(Njets)

        textLabel = '\n'.join(( "%s"%(md), "%s"%(ch), nj ))
        ax[0].text(0.03, 0.9, textLabel, transform=ax[0].transAxes, color="black", fontsize=20, fontweight='normal',  va='center', ha='left')

        #
        if lineY != None:
            l1 = ml.Line2D([0.0, 1.10], [lineY, lineY], color="black", linewidth=2, linestyle="dashed")
            l2 = ml.Line2D([0.0, 1.10], [lineY, lineY], color="black", linewidth=2, linestyle="dashed")
            ax[0].add_line(l1); 
            ax[1].add_line(l2); 

        fig.tight_layout()
        fig.subplots_adjust(top=0.95, right=0.95)
        fig.savefig('%s/%s_%s_Njets%s_%s_Ratio.pdf' % (self.outputDir, self.year, tag, Njets, self.channel), dpi=fig.dpi)

        plt.close(fig)

    # -------------------------------------------------------------
    # plot variable as a function of per boundary for all variances
    #   -- as part of closure correction factor study
    # -------------------------------------------------------------
    def plot_VarVsBoundary_MCcorrectionFactor_TTvar(self, var, ttVars, labels, xWidth, yMin = None, yMax = None, lineY = None, region = "", ylabel = "", tag = "", Njets = -1, colors=None, valColor=None, scaleFactor=1.0):

        x    = {ttVar : [] for ttVar in ttVars}
        y    = {ttVar : [] for ttVar in ttVars}
        yUnc = {ttVar : [] for ttVar in ttVars} 
        xUnc = {ttVar : [] for ttVar in ttVars}

        fig = plt.figure(figsize=(5, 5))
        ax = plt.gca()

        for ttVar in ttVars:
            for boundary, value in var[ttVar].items():
                x[ttVar].append(boundary)
                y[ttVar].append(value[0]*scaleFactor)
                yUnc[ttVar].append(value[1])
                xUnc[ttVar].append(xWidth)

            marker = "o"
            if "Up" in ttVar or "UP" in ttVar or "up" in ttVar:
                marker = "^"
            elif "Down" in ttVar or "DOWN" in ttVar or "down" in ttVar:
                marker = "v"

            if ttVar == "TT":
                ax.errorbar(x[ttVar], y[ttVar], yerr=yUnc[ttVar], xerr=xUnc[ttVar], label=labels[ttVar], fmt='', capsize=0, color=colors[ttVar], lw=0, elinewidth=2, marker=marker, markersize=6.0, markeredgecolor=colors[ttVar], markerfacecolor=colors[ttVar])
            elif ttVar == "None":
                ax.errorbar(x[ttVar], y[ttVar], yerr=yUnc[ttVar], label=labels[ttVar], fmt='', capsize=0, color=colors[ttVar], lw=0, elinewidth=2, marker=marker, markersize=6.0, markeredgecolor=colors[ttVar], markerfacecolor="white")
            else:
                ax.errorbar(x[ttVar], y[ttVar], label=labels[ttVar], fmt='', capsize=0, color=colors[ttVar], lw=0, elinewidth=2, marker=marker, markersize=5.0, markeredgecolor=colors[ttVar], markerfacecolor=colors[ttVar])

            if yMin != None and yMax != None:
                ax.set_ylim((yMin, yMax))

        ax.set_xlabel("Boundary Value", fontsize=14)
        ax.set_ylabel(ylabel, fontsize=14)

        iamLegend = plt.legend(ncol=2, loc='upper right', numpoints=1, frameon=False, fontsize=7, markerscale=0.8)
        ax.text(0.03, 0.96, region, transform=ax.transAxes, color=valColor, fontsize=12, fontweight='bold',   va='center', ha='left')

        ax = self.addCMSlabel(ax)

        # put model, channel, njet labels
        md = ""
        if self.model == "SYY":
            md = "Stealth SYY"
        else:
            md = self.model
            
        ch = ""
        if self.channel == "0l":
            ch = "0L"
        elif self.channel == "1l":
            ch = "1L"
        elif self.channel == "2l":
            ch = "2L"
        
        nj = ""
        if "incl" in Njets:
            nj = "$N_{jets} \geq$ %s"%(Njets.replace("incl",""))
        else:
            nj = "$N_{jets}$ = %s"%(Njets)
    
        textLabel = '\n'.join(( "%s"%(md), "%s"%(ch), nj ))
        ax.text(0.03, 0.85, textLabel, transform=ax.transAxes, color=valColor, fontsize=9, fontweight='normal',  va='center', ha='left')

        #
        if lineY != None:
            l1 = ml.Line2D([0.0, 1.05], [lineY, lineY], color="black", linewidth=2, linestyle="dashed")
            ax.add_line(l1); 

        fig.tight_layout()
        fig.subplots_adjust(top=0.95, right=0.95)

        fig.savefig('%s/%s_%s_Njets%s_%s_%s_Variances.pdf' % (self.outputDir, self.year, tag, Njets, region, self.channel), dpi=fig.dpi)

        plt.close(fig)


    # -------------------------------------------------------
    # plot closure per njets for MC closure correction factor
    #   -- MC correction factor
    #   -- Data closure (Data with non-TT subtruction)
    #   -- MC corrected Data Closure
    # -------------------------------------------------------
    def plot_MCcorr_ClosureNjets(self, bkgPred, bkgObs, dataPred, dataObs, Njets, name = '', closureTag = '', bkgTag = '', valColor=None):

        x         = []; xUnc         = []
        bkg_obs   = []; bkg_obsUnc   = [] 
        bkg_pred  = []; bkg_predUnc  = []
        data_obs  = []; data_obsUnc  = [] 
        data_pred = []; data_predUnc = []

        MC_closureCorrection     = []; MC_closureCorrection_Unc     = []
        MC_corrected_PredData    = []; MC_corrected_PredData_Unc    = []
        DataClosure              = []; DataClosure_Unc              = [] 
        MC_corrected_DataClosure = []; MC_corrected_DataClosure_Unc = []


        # ------------------------------
        # calculate ratio and ratio unc.
        # ------------------------------
        def cal_ratioUnc(numerator, denominator):

            ratio    = ( numerator[0] / denominator[0] )
            ratioUnc = ( (numerator[1] / denominator[0])**2.0 + ( denominator[1] * numerator[0] / denominator[0]**2.0 )**2.0 )**0.5

            return ratio, ratioUnc

        # -------------------
        # loop over njet bins
        # -------------------
        for i in range(0, len(Njets)):

            if bkgObs[i][1] != 0.0:

                x.append(float((Njets[i]).replace("incl","")))
                xUnc.append(0.5)

                # bkg
                bkg_pred.append(bkgPred[i][0])
                bkg_predUnc.append(bkgPred[i][1])
                bkg_obs.append(bkgObs[i][0])
                bkg_obsUnc.append(bkgObs[i][1])

                # MC closure correction 
                mc_closureCorrection, mc_closureCorrection_Unc = cal_ratioUnc(bkgObs[i], bkgPred[i]) 
                MC_closureCorrection.append(mc_closureCorrection)
                MC_closureCorrection_Unc.append(mc_closureCorrection_Unc)

                # MC corrected Predicted Data
                mc_corrected_PredData        = ( mc_closureCorrection * dataPred[i][0] )
                mc_corrected_PredData_Unc    = math.sqrt( (mc_closureCorrection * dataPred[i][1])**2.0 + (mc_closureCorrection_Unc * dataPred[i][0])**2.0  )
                MC_corrected_PredData.append(mc_corrected_PredData)
                MC_corrected_PredData_Unc.append(mc_corrected_PredData_Unc)

                # -------------------------------------------- 
                # append the information of 7th, 8th Njet bins
                # -------------------------------------------- 
                if i <= 1: # for 7th, 8th Njet bins
                #if i < 1:   # for 7th Njet bin 
                    # data 
                    data_pred.append(dataPred[i][0])
                    data_predUnc.append(dataPred[i][1])
                    data_obs.append(dataObs[i][0])
                    data_obsUnc.append(dataObs[i][1])

                    # Data Closure
                    dataClosure, dataClosure_Unc = cal_ratioUnc(dataPred[i], dataObs[i])
                    DataClosure.append(dataClosure)
                    DataClosure_Unc.append(dataClosure_Unc)

                    # MC corrected Data closure
                    mc_corrected_DataClosure, mc_corrected_DataClosure_Unc = cal_ratioUnc((mc_corrected_PredData,mc_corrected_PredData_Unc), dataObs[i])
                    MC_corrected_DataClosure.append(mc_corrected_DataClosure)
                    MC_corrected_DataClosure_Unc.append(mc_corrected_DataClosure_Unc)

                else:
                    data_pred.append(999999999999)
                    data_predUnc.append(0)
                    data_obs.append(999999999999)
                    data_obsUnc.append(0)
                    
                    DataClosure.append(-999)
                    DataClosure_Unc.append(0)
        
                    MC_corrected_DataClosure.append(-999)
                    MC_corrected_DataClosure_Unc.append(0) 

                # -----------------------
                # plotting all them above
                # -----------------------
                fig = plt.figure(figsize=(5, 5))
                ax  = fig.add_subplot(111)
                fig.subplots_adjust(hspace=0.0, left=0.15, right=0.95, top=0.95)

                lowerNjets  = float(Njets[0])
                higherNjets = float( (Njets[(-1)]).replace("incl","") )
    
                ax.spines['top'].set_color('none')
                ax.spines['bottom'].set_color('none')
                ax.spines['left'].set_color('none')
                ax.spines['right'].set_color('none')
                ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
        
                # unweighted event counts
                ax1 = fig.add_subplot(10, 1, (1,6))
                fig.subplots_adjust(left=0.15, right=0.95)
                ax1.set_xlim([lowerNjets - 0.5, higherNjets + 0.5])

                ax1 = self.addCMSlabel(ax1)

                # put model, channel labels
                md = ""
                if self.model == "SYY":
                    md = "Stealth SYY"
                else:
                    md = self.model

                ch = ""
                if self.channel == "0l":
                    ch = "0L"
                elif self.channel == "1l":
                    ch = "1L"
                elif self.channel == "2l":
                    ch = "2L"

                textLabel = '\n'.join(( "%s"%(name), "%s"%(md), "%s"%(ch) ))
                ax1.text(0.66, 0.70, textLabel, transform=ax.transAxes, color="black",  fontsize=9,  fontweight='normal', va='center', ha='left' )
                #ax1.text(0.5, 0.95, name,      transform=ax.transAxes, color=valColor, fontsize=11, fontweight='bold',   va='center', ha='center') 

                ax1.errorbar(x, bkg_pred,  yerr=bkg_predUnc,  label=r'Predicted $t\bar{t}+jets$', xerr=xUnc, fmt='', capsize=0, color='red',       lw=0, elinewidth=2, marker='o', markeredgecolor='red',       markerfacecolor='red',       markersize=5.0)
                ax1.errorbar(x, bkg_obs,   yerr=bkg_obsUnc,   label=r'Observed $t\bar{t}+jets$',  xerr=xUnc, fmt='', capsize=0, color='black',     lw=0, elinewidth=2, marker='o', markeredgecolor='black',     markerfacecolor='black',     markersize=5.0)
                #ax1.errorbar(x, data_pred, yerr=data_predUnc, label='Predicted Data',             xerr=xUnc, fmt='', capsize=0, color='palegreen', lw=0, elinewidth=2, marker='o', markeredgecolor='palegreen', markerfacecolor='palegreen', markersize=5.0)
                #ax1.errorbar(x, data_obs,  yerr=data_obsUnc,  label='Observed Data',              xerr=xUnc, fmt='', capsize=0, color='green',     lw=0, elinewidth=2, marker='o', markeredgecolor='green',     markerfacecolor='green',     markersize=5.0)
                ax1.set_xticklabels([])
                ax1.set_ylabel('Num. Events', fontsize=11)
                ax1.set_yscale('log')
                ax1.set_ylim([5.0, 2e4])
            
                # MC closure correction & Data Closure & MC corrected Data Closure
                ax2 = fig.add_subplot(10, 1, (7,10))
                fig.subplots_adjust(left=0.15, right=0.95)
                ax2.set_xlim([lowerNjets - 0.5, higherNjets + 0.5]) 
                ax2.errorbar(x, MC_closureCorrection,     yerr=MC_closureCorrection_Unc,     xerr=xUnc, label='Closure Correction',     fmt='', capsize=0, color='limegreen',      lw=0, elinewidth=2, marker='o', markeredgecolor='limegreen',      markerfacecolor='limegreen',      markersize=5.0)
                ax2.errorbar(x, DataClosure,              yerr=DataClosure_Unc,              xerr=xUnc, label='Data Closure',           fmt='', capsize=0, color='crimson',        lw=0, elinewidth=2, marker='o', markeredgecolor='crimson',        markerfacecolor='crimson',        markersize=5.0)
                #ax2.errorbar(x, MC_corrected_DataClosure, yerr=MC_corrected_DataClosure_Unc, xerr=xUnc, label='Corrected Data Closure', fmt='', capsize=0, color='cornflowerblue', lw=0, elinewidth=2, marker='o', markeredgecolor='cornflowerblue', markerfacecolor='cornflowerblue', markersize=5.0)
                ax2.axhline(y=1.0, color='black', linestyle='dashed', lw=1.5)
                ax2.grid(axis='y', color='black', linestyle='dashed', which='both')
                ax2.set_ylabel('Closure Correction [MC]', fontsize=11)
                ax2.set_ylim([0.7, 1.3])

                plt.xticks([int(Njet.replace("incl","")) for Njet in Njets], fontsize=10)
                plt.xlabel('Number of jets', fontsize=12)    

                ax1.legend(loc='upper right', numpoints=1, frameon=False)
                ax2.legend(loc='upper right', numpoints=1, frameon=False, fontsize=7)
    
                fig.savefig('%s/%s_Njets_Region_A_PredVsActual_ClosureCorr_dataVsMC_%s_%s_%s_%s.pdf' % (self.outputDir, self.year, name, closureTag, bkgTag, self.channel))
   
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


    # ---------------------------------------------------------------------------------------
    # plot correction, data-based systematic, FSR systematic, and signal contamination effect
    # ---------------------------------------------------------------------------------------
    #def plot_ClosureAll(self, closeCorr, dataSys, fsrSys, sigContam, dataStatUnc, Njets, name = '', closureTag = '', bkgTag = '', valColor=None, isBlind=False):
    def plot_ClosureAll(self, closeCorr, dataSys, fsrSys, dataStatUnc, Njets, name = '', closureTag = '', bkgTag = '', valColor=None, isBlind=False):
        
        print("closeCorr: ", closeCorr)
        print("dataSys: ", dataSys)
        print("fsrSys: ", fsrSys)
        print("Njets: ", Njets)

        plt_closeCorr = []; plt_closeCorr_UncY = [];
        plt_dataSys = [];   plt_dataSys_UncY = [];  
        plt_dataSys_OnlyUnc = []; plt_dataSys_OnlyUncY = [];
        plt_fsrSys = [];    plt_fsrSys_UncY = [];   
        #plt_sigContam = []; plt_sigContam_UncY = [];
        plt_totalSysUnc = []
        plt_totalUnc = []

        firstNJ = Njets[0]

        skipRest = False
        for nj in Njets:
            plt_closeCorr.append(closeCorr[nj][0])
            plt_closeCorr_UncY.append(closeCorr[nj][1])
            if dataSys[nj][0] == 1.0 or skipRest:
                skipRest = True
                plt_dataSys_OnlyUnc.append(1.0)
                plt_dataSys_OnlyUncY.append(dataSys[nj][1])
            else:
                plt_dataSys.append(dataSys[nj][0])
                plt_dataSys_UncY.append(dataSys[nj][1])
            plt_fsrSys.append(fsrSys[nj][0])
            plt_fsrSys_UncY.append(fsrSys[nj][1])
            #plt_sigContam.append(sigContam[nj][0])
            #plt_sigContam_UncY.append(sigContam[nj][1])
            plt_totalSysUnc.append(math.sqrt((1-fsrSys[nj][0])**2 + (closeCorr[nj][1])**2 + (1-dataSys[firstNJ][0])**2))
            plt_totalUnc.append(math.sqrt((1-fsrSys[nj][0])**2 + (closeCorr[nj][1])**2 + dataStatUnc[nj]**2 + (1-dataSys[firstNJ][0])**2))


        xe = [0.5] * 5
        x_noe = [0.0] * 5
        y_noe = [0.0] * 5
            
        # ------------------------------------------   
        # make a band for signal contam
        # ------------------------------------------  
        def makeErrorBoxes(x, zeros, xUnc, yUnc, hatch=""):

            # list for all the error patches
            errorBoxes = []

            # loop over data points / create box from errors at each point
            for xc, yc, xe, ye in zip(x, zeros, xUnc, yUnc):
                rect = Rectangle( (xc-xe, yc-ye), 2*xe, 2*ye )
                errorBoxes.append(rect)

            # create patch collection with specified colour/alpha
            pc = PatchCollection(errorBoxes, facecolor='gray', alpha=0.5, edgecolor='none')

            return pc
 
        # ------------------
        # plot usual closure
        # ------------------
        fig = plt.figure(figsize=(5, 5))
        ax  = fig.add_subplot(111)
        fig.subplots_adjust(hspace=0, left=0.15, right=0.95, top=0.95)
    
        lowerNjets  = float(Njets[0])
        higherNjets = float( (Njets[(-1)]).replace("incl","") )
   
        x = range(int(lowerNjets), int(higherNjets) + 1)
        x_ds = []
        x_dsu = []

        skipRest = False
        for nj in Njets:
            if dataSys[nj][0] == 1.0 or skipRest:
                skipRest = True
                x_dsu.append([i for i in x if str(i) in nj][0])
            else:
                x_ds.append([i for i in x if str(i) in nj][0])
 
        #ax.spines['top'].set_color('none')
        #ax.spines['bottom'].set_color('none')
        #ax.spines['left'].set_color('none')
        #ax.spines['right'].set_color('none')
        #ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
    
        # unweighted event counts
        #ax1 = fig.add_subplot(4, 1, (1, 2))  
        #fig.subplots_adjust(left=0.15, right=0.95)
        #ax1.set_xlim([lowerNjets - 0.5, higherNjets + 0.5])
        
        ax = self.addCMSlabel(ax)

        # put model, channel labels
        md = ""
        if self.model == "SYY":
            md = "Stealth SYY"
        else:
            md = self.model

        ch = ""
        if self.channel == "0l":
            ch = "0L"
        elif self.channel == "1l":
            ch = "1L"
        elif self.channel == "2l":
            ch = "2L"

        #if ch == "1L" and md == "Stealth SYY" and plt_sigContam[-1] > 1.0:
        #    plt_sigContam[-1] = 0.9
            

        textLabel = '\n'.join(( "%s"%(name), "%s"%(md), "%s"%(ch) ))
        ax.text(0.02, 0.90, textLabel, transform=ax.transAxes, color="black",  fontsize=9,  fontweight='normal', va='center', ha='left' )
        #ax.text(0.5, 0.95, name,      transform=ax.transAxes, color=valColor, fontsize=11, fontweight='bold',   va='center', ha='center') 

        ax.set_ylabel('Closure Effect', fontsize=11)
        ax.set_xlabel('$N_{Jets}$', fontsize=11)
        ax.errorbar([i-0.25 for i in x], plt_closeCorr, yerr=plt_closeCorr_UncY, label=r'Closure Correction', xerr=x_noe, fmt='', capsize=0, color='red',   lw=0, elinewidth=1, marker='o', markeredgecolor='red',   markerfacecolor='red',   markersize=8.0)
        ax.errorbar([x_ds[0]], [plt_dataSys[0]], yerr=[plt_dataSys_UncY[0]], label=r'Data Residual Non-Closure (Systematic)', xerr=[x_noe[0]], fmt='', capsize=0, color='blue',   lw=0, elinewidth=1, marker='o', markeredgecolor='blue',   markerfacecolor='blue',   markersize=8.0)
        ax.errorbar(x_ds[1:len(plt_dataSys)], plt_dataSys[1:len(plt_dataSys)], yerr=plt_dataSys_UncY[1:len(plt_dataSys)], label=r'Data Residual Non-Closure', xerr=x_noe[1:len(plt_dataSys)], fmt='', capsize=0, color='blue',   lw=0, elinewidth=1, marker='^', markeredgecolor='blue',   markerfacecolor='blue',   markersize=8.0)
        ax.errorbar(x_dsu, plt_dataSys_OnlyUnc, yerr=plt_dataSys_OnlyUncY, label=r'Blinded Data Residual Non-Closure Unc.', xerr=x_noe[len(plt_dataSys):], fmt='', capsize=0, color='blue',   lw=0, elinewidth=1, marker='^', markeredgecolor='blue',   markerfacecolor='none',   markersize=8.0)
        ax.errorbar([i+0.25 for i in x], plt_fsrSys, yerr=plt_fsrSys_UncY, label=r'FSR Systematic', xerr=x_noe, fmt='', capsize=0, color='green',   lw=0, elinewidth=1, marker='o', markeredgecolor='green',   markerfacecolor='green',   markersize=8.0)
        #eb = ax.errorbar(x, plt_sigContam, label=r'Effect of Signal Contam.', xerr=xe, fmt='', capsize=0, color='black',   lw=0, elinewidth=1, marker='', markersize=0.0)
        #eb[-1][0].set_linestyle('--')

        pc = makeErrorBoxes(x, [1.0]*5, xe, plt_totalSysUnc)
        ax.add_collection(pc)

        pc = makeErrorBoxes(x, [1.0]*5, xe, plt_totalUnc)
        ax.add_collection(pc)

        sys = mpl.patches.Patch(color='gray', label='Syst. Unc')
        statSys = mpl.patches.Patch(color='gray', label='Stat. + Syst. Unc', alpha=0.5)
        
        handles, labels = plt.gca().get_legend_handles_labels()
        handles.append(sys)
        handles.append(statSys)
        
        #plt.show()
        plt.xticks([int(Njet.replace("incl","")) for Njet in Njets])
        plt.grid(axis = 'y', linestyle = '-')
   
        #box = ax.get_position()
        #ax.set_position([box.x0, box.y0, box.width * 0.9, box.height*0.85])

        lgd = ax.legend(handles=handles, loc='center left', numpoints=1, frameon=False, bbox_to_anchor=(1, 0.5))
   
        #plt.subplots_adjust(top=0.88)
 
        fig.savefig('%s/%s_AllClosure_Njets_%s.pdf' % (self.outputDir, self.year, self.channel), bbox_inches='tight', dpi=300)
   
        plt.ylim(0.5, 2.0)

        fig.savefig('%s/%s_AllClosure_Njets_%s_zoom.pdf' % (self.outputDir, self.year, self.channel), bbox_inches='tight', dpi=300)
 
        plt.close(fig)

        return ([i-0.25 for i in x], plt_closeCorr, plt_closeCorr_UncY), (x_ds, plt_dataSys, plt_dataSys_UncY), (x_dsu, plt_dataSys_OnlyUnc, plt_dataSys_OnlyUncY), ([i+0.25 for i in x], plt_fsrSys, plt_fsrSys_UncY), (x, [1.0]*5, plt_totalSysUnc), (x, [1.0]*5, plt_totalUnc)
