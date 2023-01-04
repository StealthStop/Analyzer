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

    def __init__(self, outputDir, year, model, mass, channel, cmsLabel):
        self.year            = year
        self.model           = model
        self.mass            = mass
        self.channel         = channel
        self.cmsLabel        = cmsLabel
        self.cmap            = plt.cm.jet
        self.cmap.set_under("w", 1)

        self.outputDir = outputDir

    def addCMSlabel(self, ax):

        ax.text(0.0,  1.003, 'CMS',                     transform=ax.transAxes, fontsize=16, fontweight='bold',   va='bottom', ha='left')
        ax.text(0.15, 1.010, '%s'%(self.cmsLabel),      transform=ax.transAxes, fontsize=11, fontstyle='italic',  va='bottom', ha='left')
        ax.text(1.0,  1.010, '%s (13 TeV)'%(self.year), transform=ax.transAxes, fontsize=11, fontweight='normal', va='bottom', ha='right')
    
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
        ax1 = fig.add_subplot(4, 1, (1, 2))  
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
            ch = "Fully-Hadronic"
        elif self.channel == "1l":
            ch = "Semi-Leptonic"
        elif self.channel == "2l":
            ch = "Fully-Leptonic"

        textLabel = '\n'.join(( "%s"%(name), "%s"%(md), "%s"%(ch) ))
        ax1.text(0.66, 0.80, textLabel, transform=ax.transAxes, color="black",  fontsize=9,  fontweight='normal', va='center', ha='left' )
        #ax1.text(0.5, 0.95, name,      transform=ax.transAxes, color=valColor, fontsize=11, fontweight='bold',   va='center', ha='center') 

        ax1.set_ylabel('Unweighted Event Counts', fontsize=11)
        ax1.errorbar(x, pred, yerr=predUnc, label=r'Predicted $t\bar{t}+jets$', xerr=xUnc, fmt='', capsize=0, color='red',   lw=0, elinewidth=2, marker='o', markeredgecolor='red',   markerfacecolor='red',   markersize=5.0)
        ax1.errorbar(x, obs,  yerr=obsUnc,  label=r'Observed $t\bar{t}+jets$',  xerr=xUnc, fmt='', capsize=0, color='black', lw=0, elinewidth=2, marker='o', markeredgecolor='black', markerfacecolor='black', markersize=5.0)

        # closure   
        ax2 = fig.add_subplot(4, 1, 3)
        ax2.set_xlim([lowerNjets - 0.5, higherNjets + 0.5])
        ax2.errorbar(x, abcdError, yerr=abcdErrorUnc, xerr=xUnc, fmt='', capsize=0, color='blue', lw=0, elinewidth=2, marker='o', markeredgecolor='blue', markerfacecolor='blue', markersize=5.0)
        ax2.axhline(y=0.0, color='black', linestyle='dashed', lw=1.5)
        ax2.grid(axis='y', color='black', linestyle='dashed', which='both')
        ax2.set_ylabel('Non-Closure', fontsize=12)
        ax2.set_ylim([-0.59, 0.59])   

        # pull 
        ax3 = fig.add_subplot(4, 1, 4)
        ax3.set_xlim([lowerNjets - 0.5, higherNjets + 0.5])
        ax3.errorbar(x, abcdPull, yerr=abcdPullUnc, xerr=xUnc, fmt='', capsize=0, color='purple', lw=0, elinewidth=2, marker='o', markeredgecolor='purple', markerfacecolor='purple', markersize=5.0)
        ax3.axhline(y=0.0, color='black', linestyle='dashed', lw=1.5)
        ax3.grid(axis='y', color='black', linestyle='dashed', which='both')
        ax3.set_xlabel('Number of jets', fontsize=12)
        ax3.set_ylabel('Pull',           fontsize=12) 
        ax3.set_ylim([-5.9, 5.9])
 
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
            ch = "Fully-Hadronic"
        elif self.channel == "1l":
            ch = "Semi-Leptonic"
        elif self.channel == "2l":
            ch = "Fully-Leptonic"

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
            ax.text(0.05, 0.10, textLabel, transform=ax.transAxes, color="maroon", fontsize=9, fontweight='bold',  va='center', ha='left')
        if name == "Val_BD":
            ax.text(0.95, 0.10, textLabel, transform=ax.transAxes, color="maroon", fontsize=9, fontweight='bold',  va='center', ha='right')
        if name == "Val_CD":
            ax.text(0.05, 0.90, textLabel, transform=ax.transAxes, color="maroon", fontsize=9, fontweight='bold',  va='center', ha='left')
        if name == "Val_D":
            ax.text(0.95, 0.90, textLabel, transform=ax.transAxes, color="maroon", fontsize=9, fontweight='bold',  va='center', ha='right')

        #
        l1 = ml.Line2D([c1, c1], [0.0, 1.0], color="maroon", linewidth=2, linestyle="dashed"); l2 = ml.Line2D([0.0, 1.0], [c2, c2], color="maroon", linewidth=2, linestyle="dashed")
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

        x    = {region : [] for region in regions}
        y    = {region : [] for region in regions}
        yUnc = {region : [] for region in regions} 
        xUnc = {region : [] for region in regions}

        fig = plt.figure(figsize=(5, 5))
        ax = plt.gca()

        for region in regions:
            for boundary, value in var[region].items():
                x[region].append(boundary)
                y[region].append(value[0])
                yUnc[region].append(value[1])
                xUnc[region].append(xWidth)

            ax.errorbar(x[region], y[region], yerr=yUnc[region], label="%s"%(region), xerr=xUnc[region], fmt='', capsize=0, color=color[region], lw=0, elinewidth=2, marker="o", markeredgecolor=color[region], markerfacecolor=color[region])

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
            ch = "Fully-Hadronic"
        elif self.channel == "1l":
            ch = "Semi-Leptonic"
        elif self.channel == "2l":
            ch = "Fully-Leptonic"

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
            ch = "Fully-Hadronic"
        elif self.channel == "1l":
            ch = "Semi-Leptonic"
        elif self.channel == "2l":
            ch = "Fully-Leptonic"
        
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
                    ch = "Fully-Hadronic"
                elif self.channel == "1l":
                    ch = "Semi-Leptonic"
                elif self.channel == "2l":
                    ch = "Fully-Leptonic"

                textLabel = '\n'.join(( "%s"%(name), "%s"%(md), "%s"%(ch) ))
                ax1.text(0.66, 0.70, textLabel, transform=ax.transAxes, color="black",  fontsize=9,  fontweight='normal', va='center', ha='left' )
                #ax1.text(0.5, 0.95, name,      transform=ax.transAxes, color=valColor, fontsize=11, fontweight='bold',   va='center', ha='center') 

                ax1.errorbar(x, bkg_pred,  yerr=bkg_predUnc,  label=r'Predicted $t\bar{t}+jets$', xerr=xUnc, fmt='', capsize=0, color='red',       lw=0, elinewidth=2, marker='o', markeredgecolor='red',       markerfacecolor='red',       markersize=5.0)
                ax1.errorbar(x, bkg_obs,   yerr=bkg_obsUnc,   label=r'Observed $t\bar{t}+jets$',  xerr=xUnc, fmt='', capsize=0, color='black',     lw=0, elinewidth=2, marker='o', markeredgecolor='black',     markerfacecolor='black',     markersize=5.0)
                ax1.errorbar(x, data_pred, yerr=data_predUnc, label='Predicted Data',             xerr=xUnc, fmt='', capsize=0, color='palegreen', lw=0, elinewidth=2, marker='o', markeredgecolor='palegreen', markerfacecolor='palegreen', markersize=5.0)
                ax1.errorbar(x, data_obs,  yerr=data_obsUnc,  label='Observed Data',              xerr=xUnc, fmt='', capsize=0, color='green',     lw=0, elinewidth=2, marker='o', markeredgecolor='green',     markerfacecolor='green',     markersize=5.0)
                ax1.set_xticklabels([])
                ax1.set_ylabel('Unweighted Event Counts', fontsize=11)
                ax1.set_yscale('log')
                ax1.set_ylim([5.0, 2e4])
            
                # MC closure correction & Data Closure & MC corrected Data Closure
                ax2 = fig.add_subplot(10, 1, (7,10))
                fig.subplots_adjust(left=0.15, right=0.95)
                ax2.set_xlim([lowerNjets - 0.5, higherNjets + 0.5]) 
                ax2.errorbar(x, MC_closureCorrection,     yerr=MC_closureCorrection_Unc,     xerr=xUnc, label='Closure Correction',     fmt='', capsize=0, color='limegreen',      lw=0, elinewidth=2, marker='o', markeredgecolor='limegreen',      markerfacecolor='limegreen',      markersize=5.0)
                ax2.errorbar(x, DataClosure,              yerr=DataClosure_Unc,              xerr=xUnc, label='Data Closure',           fmt='', capsize=0, color='crimson',        lw=0, elinewidth=2, marker='o', markeredgecolor='crimson',        markerfacecolor='crimson',        markersize=5.0)
                ax2.errorbar(x, MC_corrected_DataClosure, yerr=MC_corrected_DataClosure_Unc, xerr=xUnc, label='Corrected Data Closure', fmt='', capsize=0, color='cornflowerblue', lw=0, elinewidth=2, marker='o', markeredgecolor='cornflowerblue', markerfacecolor='cornflowerblue', markersize=5.0)
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
