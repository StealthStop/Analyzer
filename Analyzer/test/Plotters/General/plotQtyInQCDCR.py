import os
import glob
import argparse
import numpy as np
import multiprocessing as mp

import matplotlib as mpl
mpl.use("Agg")

import matplotlib.pyplot as plt

# Looks very similar to Arial/Helvetica
plt.rc("font", family = "Nimbus Sans")

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH2.SetDefaultSumw2()

def combineHistos(file, histoStub, process, pattern, njets, regions, channel, training, histos):

    tfile = ROOT.TFile(file, "READ")

    # Extract all base histograms from file, which are TH1 for different A, B, C, D regions
    # and all as a function of selection on dR(b1,b2)
    tempHistos = {}
    for key in tfile.GetListOfKeys():
        obj  = key.ReadObj()
        name = obj.GetName()
        if obj.InheritsFrom(ROOT.TH1.Class()):
            obj.SetDirectory(0)
            tempHistos[name] = obj
    
    # From the base histograms, where each bin is a total number of events,
    # Create equivalent histograms for A*D, B*C, closure, and closure correction
    for njet in njets:
        for region in regions:
    
            histoName = "h_%s_%s_%s_%s_Njets%s_Nbjets2incl"%(histoStub,channel,training,region,njet)
       
            tempHisto = tempHistos[histoName].Clone("%s_%s_%s"%(process,njet,region))
            tempHisto.Reset()

            content = 0.0; errorSqr = 0.0
            for iBin in range(1, tempHisto.GetNbinsX()+1):

                content  += tempHistos[histoName].GetBinContent(iBin)
                errorSqr += tempHistos[histoName].GetBinError(iBin)**2.0

                tempHisto.SetBinContent(iBin, content)
                tempHisto.SetBinError(iBin, errorSqr**0.5)
   
            histos["%s_%s_%s"%(process,njet,region)] = tempHisto

        histoBC = histos["%s_%s_B"%(process,njet)].Clone("%s_%s_%s_BC"%(process,njet,region)); histoBC.Reset()
        histoAD = histos["%s_%s_D"%(process,njet)].Clone("%s_%s_%s_AD"%(process,njet,region)); histoAD.Reset()
    
        histoBC.Multiply(histos["%s_%s_B"%(process,njet)], histos["%s_%s_C"%(process,njet)])
        histoAD.Multiply(histos["%s_%s_A"%(process,njet)], histos["%s_%s_D"%(process,njet)])

        histoCorr    = histoAD.Clone("%s_%s_Corr"%(process,njet));    histoCorr.Reset()
        histoClosure = histoAD.Clone("%s_%s_Closure"%(process,njet)); histoClosure.Reset()
        
        histoCorr.Divide(histoAD, histoBC)
        histoClosure.Divide(histoBC, histoAD)

        histos["%s_%s_Corr"%(process,njet)]    = histoCorr
        histos["%s_%s_Closure"%(process,njet)] = histoClosure

    tfile.Close()

def hist2Array(histo):
    arr = np.array([[histo.GetBinContent(i), histo.GetBinError(i)] for i in range(1, histo.GetNbinsX()+2)])
    return arr

def addCMSlabel(ax, scale = 1.0, **kwargs):

    # Put CMS by the preferred standard of within the frame
    ax.text(0.03, 0.970, "CMS",              transform = ax.transAxes, fontsize = 18*scale, fontweight = "bold",   va = "top", ha = "left")
    ax.text(0.03, 0.900, "Work in Progress", transform = ax.transAxes, fontsize = 11*scale, fontstyle  = "italic", va = "top", ha = "left")

    if "year" in kwargs:
        ax.text(1.00, 1.007, "%s (13 TeV)"%(kwargs["year"]), transform = ax.transAxes, fontsize = 12*scale, fontweight = "normal", va = "bottom", ha = "right")

    return ax

def axisOptions(ax, yMin, yMax, xMin, xMax, yLabel, doCMS = False, doGrid = False):
    ax.set_ylabel(yLabel, fontsize = 18)

    ax.set_ylim([yMin, yMax])
    ax.set_xlim([xMin, xMax])

    if doCMS:
        ax = addCMSlabel(ax, scale = 1.35, year = "Run2UL")

    if doGrid:
        ax.tick_params(axis = "x", labelsize = 16)
    else:
        ax.tick_params(axis = "x", labelsize = 0)
       
    ax.tick_params(axis = "y", labelsize = 16)

    ax.xaxis.set_tick_params(which = "both", direction = "in", bottom = True, top   = True)
    ax.yaxis.set_tick_params(which = "both", direction = "in", left   = True, right = True)

    ax.minorticks_on()

    ax.set_ylabel(yLabel, fontsize = 18)

    if doGrid:
        ax.set_axisbelow(True)
        ax.grid(axis = "both", color = "0.90")

    return ax

# Make line histos with partially transparent error band
def plot2Panel(ax, zippedTuple, xEdges):
    handles = []; labels = []
    for histo, label, color in zippedTuple:
        histoFill  = ax.fill_between(xEdges, np.subtract(histo[:, 0], histo[:, 1]), np.sum(histo, axis = 1), step = "post", lw = 0, facecolor = color, edgecolor = (1,1,1,0), alpha = 0.3)
        histoStep, = ax.step(xEdges, histo[:,0], where = "post", color = color, lw = 2)

        handles.append((histoFill, histoStep))
        labels.append(label)

    return handles, labels

def plotTwoPanel(outputName, topPanelTuple, bottomPanelTuple, xEdges, mcBasedSyst = None, dataBasedSyst = None, xMin = None, xMax = None, yMin1 = None, yMax1 = None, yMin2 = None, yMax2 = None, yLabel1 = None, yLabel2 = None):

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize = (6, 10))
    fig.subplots_adjust(left = 0.18, top = 0.97, right = 0.97, bottom = 0.06, hspace = 0)

    md = ""
    if "SYY" in outputName:
        md = "Stealth SYY"
    else:
        md = "RPV"

    ch = ""
    if "0l" in outputName:
        ch = "Fully-Hadronic"
    elif "1l" in outputName:
        ch = "Semi-Leptonic"
    elif "2l" in outputName:
        ch = "Fully-Leptonic"

    aux = ""
    op = " = "
    if "incl" in outputName:
        op = "\geq"
    aux = r"$N_{\rm{jets}} %s %s$"%(op, outputName.split("Njets")[-1].split("_")[0].replace("incl", ""))

    textLabel = "\n".join(( "%s"%(md), "%s"%(ch), "%s"%(aux) ))
    ax1.text(0.03, 0.82, textLabel, transform = ax1.transAxes, color = "black", fontsize = 14, fontweight = "normal", va = "top", ha = "left")

    if "dR" in outputName:
        ax2.set_xlabel(r"$\Delta R(b_1, b_2)<x$", fontsize = 18)
    else:
        ax2.set_xlabel(r"$M(b_1, b_2)<x$", fontsize = 18)

    ax1 = axisOptions(ax1, yMin1, yMax1, xMin, xMax, yLabel1, doCMS = True)
    ax2 = axisOptions(ax2, yMin2, yMax2, xMin, xMax, yLabel2, doGrid = True)

    if "FracsEvts" in outputName:
        plt.axhline(y = 0.1, color = "midnightblue", linewidth = 2, linestyle = "dotted")
    else:
        plt.axhline(y = 1.0, color = "midnightblue", linewidth = 2, linestyle = "dotted")

    handles = []; labels = []
    if "None" not in mcBasedSyst.__class__.__name__:
        mcBasedFill = ax2.fill_between(xEdges, np.subtract(mcBasedSyst[:, 0], mcBasedSyst[:, 1]), np.sum(mcBasedSyst, axis = 1), step = "post", hatch = "\\\\\\", lw = 0, edgecolor = "gray", alpha = 0.5, facecolor = "none")

        handles.append((mcBasedFill))
        labels.append("Correction Stat. Syst.")
    if "None" not in dataBasedSyst.__class__.__name__:
        dataBasedFill = ax2.fill_between(xEdges, np.subtract(dataBasedSyst[:, 0], dataBasedSyst[:, 1]), np.sum(dataBasedSyst, axis = 1), step = "post", hatch = "///", lw = 0, edgecolor = "gray", alpha = 0.5, facecolor = "none")

        handles.append((dataBasedFill))
        labels.append("Data-Based Syst.")

    handles1, labels1 = plot2Panel(ax1, topPanelTuple,    xEdges)
    handles2, labels2 = plot2Panel(ax2, bottomPanelTuple, xEdges)

    ax1.legend(handles1,           labels1,          loc = "upper right", numpoints = 1, frameon = False, fontsize = 15)
    ax2.legend(handles2 + handles, labels2 + labels, loc = "upper right", numpoints = 1, frameon = False, fontsize = 15)

    fig.savefig(outputName)

    plt.close(fig)

if __name__ == "__main__":
    usage = "usage: %plotQtyInQCDCR [options]"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("--path",     dest = "path",     help = "Path to ntuples",  required = True)
    parser.add_argument("--output",   dest = "output",   help = "Path to outputs",  required = True)
    parser.add_argument("--pattern",  dest = "pattern",  help = "string pattern",   required = True)
    parser.add_argument("--channel",  dest = "channel",  help = "channel",          required = True)
    parser.add_argument("--training", dest = "training", help = "training",         required = True)
    parser.add_argument("--sigMass",  dest = "sigMass",  help = "sigMass",          default  = "400")
    parser.add_argument("--histo",    dest = "histo",    help = "stub histo name",  required = True)
    parser.add_argument("--disc1",    dest = "disc1",    help = "disc1 edge",       default  = "0.6")
    parser.add_argument("--disc2",    dest = "disc2",    help = "disc2 edge",       default  = "0.6")
    
    args = parser.parse_args()
    
    histo     = args.histo
    path      = args.path
    disc1edge = args.disc1
    disc2edge = args.disc2
    pattern   = args.pattern
    output    = args.output
    channel   = args.channel
    training  = args.training
    sigMass   = args.sigMass
    
    if not os.path.exists(output):
        os.makedirs(output)
    
    cmsswDir = os.getenv("CMSSW_BASE")
    
    inputPath = os.path.abspath(path)
    
    if channel == "1l":
        njets = [str(njet) for njet in range(7, 13)]
    else:
        njets = [str(njet) for njet in range(8, 14)]
    njets[-1] += "incl"
    
    regions = ["A", "B", "C", "D"]
    
    procs = ["QCD", "Non_QCD", "Data_SingleMuon", "%s_2t6j_mStop-%s"%(training, sigMass)]
    
    # Extract and create needed histograms for weighted events in parallel
    manager = mp.Manager()
    histos = manager.dict()
    
    for proc in procs:
        histos[proc] = manager.dict()
    
    pool = mp.Pool(processes = min(4, len(procs)))
    
    for proc in procs:
        pool.apply_async(combineHistos, (inputPath + "/Run2UL_%s.root"%(proc), histo, proc, pattern, njets, regions, channel, training, histos))
    
    pool.close()
    pool.join()
    
    # Extract correction statistical uncertainty systematic and residual, post-correction systematic uncertainty
    systRootFileName = "%s/src/Analyzer/Analyzer/test/DoubleDisCo_BinEdges/Run2UL_TT_TTvar_Syst_%s_%s_%s_%s_%s.root"%(cmsswDir, training[-3:], sigMass, channel, disc1edge, disc2edge)
    
    f = ROOT.TFile.Open(systRootFileName, "READ")
    mcBasedHist   = f.Get("Run2UL_MCcorr_TTvar_TT")
    dataBasedHist = f.Get("Run2UL_maximum_MCcorrectedData_Syst_All")
    
    mcBasedSyst   = {}
    dataBasedSyst = {}
    
    for iBin, njet in zip(list(range(1, mcBasedHist.GetNbinsX()+1)), njets):
        mcBasedSystVal      = mcBasedHist.GetBinError(iBin)
        dataBasedSystVal    = abs(1.0-dataBasedHist.GetBinContent(iBin))
        mcBasedSyst[njet]   = mcBasedSystVal
        dataBasedSyst[njet] = dataBasedSystVal 
        
    # From the base histos extracted above (individually for each process)
    # Create additional histograms when combining specific processes together
    for njet in njets:
        for region in regions:
    
            histoQCDinData = histos["Data_SingleMuon_%s_%s"%(njet,region)].Clone("QCDinData_%s_%s"%(njet,region))
            histoQCDinData.Add(histos["Non_QCD_%s_%s"%(njet,region)], -1.0)
    
            histos["QCDinData_%s_%s"%(njet,region)] = histoQCDinData
    
            histoBkgSig = histos["QCD_%s_%s"%(njet,region)].Clone("BkgSig_%s_%s"%(njet,region))
            histoBkgSig.Add(histos["Non_QCD_%s_%s"%(njet,region)])
            histoBkgSig.Add(histos["%s_2t6j_mStop-%s_%s_%s"%(training,sigMass,njet,region)], 0.21)
    
            histos["BkgSig_%s_%s"%(njet,region)] = histoBkgSig
    
            histoBkg = histos["QCD_%s_%s"%(njet,region)].Clone("Bkg_%s_%s"%(njet,region))
            histoBkg.Add(histos["Non_QCD_%s_%s"%(njet,region)])
    
            histos["Bkg_%s_%s"%(njet,region)] = histoBkg
    
            histoQCDSig = histos["QCD_%s_%s"%(njet,region)].Clone("QCDSig_%s_%s"%(njet,region))
            histoQCDSig.Add(histos["%s_2t6j_mStop-%s_%s_%s"%(training,sigMass,njet,region)], 0.21)
    
            histos["QCDSig_%s_%s"%(njet,region)] = histoQCDSig
    
        for proc in ["QCDinData", "BkgSig", "Bkg", "QCDSig"]:
    
            histoBC = histos["%s_%s_B"%(proc,njet)].Clone("%s_%s_%s_BC"%(proc,njet,region)); histoBC.Reset()
            histoAD = histos["%s_%s_D"%(proc,njet)].Clone("%s_%s_%s_AD"%(proc,njet,region)); histoAD.Reset()
    
            histoBC.Multiply(histos["%s_%s_B"%(proc,njet)], histos["%s_%s_C"%(proc,njet)])
            histoAD.Multiply(histos["%s_%s_A"%(proc,njet)], histos["%s_%s_D"%(proc,njet)])
    
            histoCorr    = histoAD.Clone("%s_%s_Corr"%(proc,njet));    histoCorr.Reset()
            histoClosure = histoAD.Clone("%s_%s_Closure"%(proc,njet)); histoClosure.Reset()
            
            histoCorr.Divide(histoAD, histoBC)
            histoClosure.Divide(histoBC, histoAD)
        
            histos["%s_%s_Corr"%(proc,njet)]    = histoCorr
            histos["%s_%s_Closure"%(proc,njet)] = histoClosure
           
        for region in regions:
    
            histoQCDfrac = histos["QCD_%s_%s"%(njet,region)].Clone("QCD_%s_%s_frac"%(njet,region)); histoQCDfrac.Reset()
            histoQCDfrac.Divide(histos["QCD_%s_%s"%(njet,region)], histos["Bkg_%s_%s"%(njet,region)])
            histos["QCD_%s_%s_frac"%(njet,region)] = histoQCDfrac
    
            histoSigFrac = histos["%s_2t6j_mStop-%s_%s_%s"%(training,sigMass,njet,region)].Clone("Sig_%s_%s_frac"%(njet,region)); histoSigFrac.Reset()
            histoSigFrac.Divide(histos["%s_2t6j_mStop-%s_%s_%s"%(training,sigMass,njet,region)], histos["BkgSig_%s_%s"%(njet,region)])
            histos["%s_2t6j_mStop-%s_%s_%s_frac"%(training,sigMass,njet,region)] = histoSigFrac
    
        for proc in procs + ["QCDinData", "BkgSig", "Bkg", "QCDSig"]:
    
            histoCorrClosure = histos["%s_%s_Closure"%(proc,njet)].Clone("%s_%s_CorrClosure"%(proc,njet)); histoCorrClosure.Reset()
            histoCorrClosure.Multiply(histos["%s_%s_Closure"%(proc,njet)], histos["QCD_%s_Corr"%(njet)])
    
            histos["%s_%s_CorrClosure"%(proc,njet)] = histoCorrClosure
    
        # Any histogram will do for getting a set of bin edges
        tempHisto = histos.values()[0]
        lowEdges = [tempHisto.GetXaxis().GetBinLowEdge(1)]
        for iBin in range(2, tempHisto.GetNbinsX()+1):
            lowEdges.append(tempHisto.GetXaxis().GetBinLowEdge(iBin))
        lowEdges.append(lowEdges[-1]+tempHisto.GetXaxis().GetBinWidth(1))
        lowEdges = np.array(lowEdges)
    
        mcBasedSystVals   = np.array([[1.0, mcBasedSyst[njet]]   for edge in lowEdges])
        dataBasedSystVals = np.array([[1.0, dataBasedSyst[njet]] for edge in lowEdges])
    
        yScale1 = 1.5
        yScale2 = 1.3

        # Plot weighted events and weighted event fractions
        qcdValsAll = []; qcdInDataValsAll = []; sigFracsAll = []
        for region in regions:
            qcdVals       = hist2Array(histos["QCD_%s_%s"%(njet,region)])
            qcdInDataVals = hist2Array(histos["QCDinData_%s_%s"%(njet,region)])
            qcdValsAll.append(qcdVals)
            qcdInDataValsAll.append(qcdInDataVals)
    
            bkgVals     = hist2Array(histos["Bkg_%s_%s"%(njet,region)])
            sigVals     = hist2Array(histos["%s_2t6j_mStop-%s_%s_%s"%(training,sigMass,njet,region)])
    
            qcdFracVals = hist2Array(histos["QCD_%s_%s_frac"%(njet,region)])
            sigFracVals = hist2Array(histos["%s_2t6j_mStop-%s_%s_%s_frac"%(training,sigMass,njet,region)])
            sigFracsAll.append(sigFracVals)
    
            topTuple    = zip([qcdVals, bkgVals, sigVals], ["QCD (MC)", "Total Bkg (MC)", "Signal (MC)"], ["#85C2A3", "gray", "darkturquoise"])
            bottomTuple = zip([qcdFracVals, sigFracVals],  ["QCD Purity", "Signal Contamination"],        ["#85C2A3", "darkturquoise"])
            outputName  = "%s/%s_%s_%s_Njets%s_%s_FracsEvts.pdf"%(output, histo, training, channel, njet, region)
            plotTwoPanel(outputName, topTuple, bottomTuple, lowEdges, xMin = 0.95, xMax = 5.0, yMin1 = 0, yMax1 = yScale1*np.amax(np.sum(bkgVals, axis = 1)), yMin2 = 0, yMax2 = 1.19, yLabel1 = "Weighted Events", yLabel2 = "Weighted Event Fractions")

        # Not too useful to look below dR(b_1, b_2) < 1.0, due to limited number of events
        dRfilter = lowEdges>1.0
    
        qcdColors = ["#6666CD", "#04A000", "#CD6600", "#990099"]
    
        # Plot signal contamination and closure post correction (with signal contamination present)
        qcdSigCorrClosure = hist2Array(histos["QCDSig_%s_CorrClosure"%(njet)])
        yRange2           = abs(1.0-np.amax(np.sum(qcdSigCorrClosure[dRfilter], axis = 1)))
        topTuple          = zip(sigFracsAll,         ["Signal Cont. in A", "Signal Cont. in B", "Signal Cont. in C", "Signal Cont. in D"], qcdColors)
        bottomTuple       = zip([qcdSigCorrClosure], ["Closure in MC (QCD CR, post correction)"], ["turquoise"])
        outputName        = "%s/%s_%s_%s_Njets%s_FracCorrClosure.pdf"%(output, histo, training, channel, njet)
        plotTwoPanel(outputName, topTuple, bottomTuple, lowEdges, xMin = 0.95, xMax = 5.0, yMin1 = -0.05, yMax1 = yScale1*max([np.amax(np.sum(sigFracsAll[i][dRfilter], axis = 1)) for i in range(0, len(sigFracsAll))]), yMin2 = 0.0, yMax2 = 2.19, yLabel1 = "Signal Contamination", yLabel2 = "Closure in MC (QCD CR, post-correction)")
    
        # Plot weighted events and the MC correction factor
        qcdCorr     = hist2Array(histos["QCD_%s_Corr"%(njet)])
        topTuple    = zip(qcdValsAll, ["QCD (MC) in A", "QCD (MC) in B", "QCD (MC) in C", "QCD (MC) in D"], qcdColors)
        bottomTuple = zip([qcdCorr],  ["Correction (QCD MC)"], ["lightsteelblue"])
        outputName  = "%s/%s_%s_%s_Njets%s_Correction.pdf"%(output, histo, training, channel, njet)
        plotTwoPanel(outputName, topTuple, bottomTuple, lowEdges, xMin = 0.95, xMax = 5.0, yMin1 = 0, yMax1 = yScale1*np.amax(np.sum(qcdValsAll[-1], axis = 1)), yMin2 = 0.0, yMax2 = 3.20, yLabel1 = "Weighted Events", yLabel2 = "Correction (QCD MC)")
    
        # Plot number of weighted QCD events and closure in data pre-correction
        qcdInDataClosure     = hist2Array(histos["QCDinData_%s_Closure"%(njet)])
        qcdInDataCorrClosure = hist2Array(histos["QCDinData_%s_CorrClosure"%(njet)])
        yRange2              = max(np.amax(np.sum(qcdInDataClosure[dRfilter], axis = 1)), np.amax(np.sum(qcdInDataCorrClosure[dRfilter], axis = 1)))
        qcdInDataLabels      = ["QCD (Data) in A",     "QCD (Data) in B", "QCD (Data) in C",      "QCD (Data) in D"]
        topTuple             = zip(qcdInDataValsAll,   qcdInDataLabels,                           qcdColors)
        bottomTuple          = zip([qcdInDataClosure], ["Closure in Data (QCD), pre-correction"], ["turquoise"])
        outputName           = "%s/%s_%s_%s_Njets%s_ClosureInData.pdf"%(output, histo, training, channel, njet)
        plotTwoPanel(outputName, topTuple, bottomTuple, lowEdges, mcBasedSystVals, dataBasedSystVals, xMin = 0.95, xMax = 5.0, yMin1 = 0, yMax1 = yScale1*np.amax(np.sum(qcdInDataValsAll[-1], axis = 1)), yMin2 = 0.0, yMax2 = min(3.19,yRange2*yScale2), yLabel1 = "Weighted Events", yLabel2 = "Closure in Data (QCD), pre-correction")
    
        # Plot number of weighted QCD events and closure in data post-correction
        topTuple             = zip(qcdInDataValsAll,       qcdInDataLabels,                            qcdColors)
        bottomTuple          = zip([qcdInDataCorrClosure], ["Closure in Data (QCD), post-correction"], ["magenta"])
        outputName           = "%s/%s_%s_%s_Njets%s_CorrClosureInData.pdf"%(output, histo, training, channel, njet)
        plotTwoPanel(outputName, topTuple, bottomTuple, lowEdges, mcBasedSystVals, dataBasedSystVals, xMin = 0.95, xMax = 5.0, yMin1 = 0, yMax1 = yScale1*np.amax(np.sum(qcdInDataValsAll[-1], axis = 1)), yMin2 = 0.0, yMax2 = min(3.19,yRange2*yScale2), yLabel1 = "Weighted Events", yLabel2 = "Closure in Data (QCD), post-correction")
