#! /usr/bin/env python

import ROOT, random, os, argparse, string, copy, math, ctypes
ROOT.PyConfig.IgnoreCommandLineOptions = True

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat("")
ROOT.gStyle.SetPaintTextFormat("3.2f")
ROOT.gStyle.SetFrameLineWidth(2)
ROOT.gStyle.SetEndErrorSize(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH2.SetDefaultSumw2()

# This class owns a histogram and can take a ROOT TH1 or
# a path to extract the histogram from
#    histo      : ROOT TH1 object to hold onto
#    rootFile   : input path to ROOT file containing histos
#    upperSplit : fraction of canvas taken by upper pad
#    lowerSplit : fraction of canvas taken by lower pad
#    histoName  : name of histo to be extracted from file
#    hinfo      : dictionary of params for customizing histo
#    pinfo      : additional params for customizing histo
#    aux        : auxiliary multiplier for offsetting y-axis label
class Histogram:

    def __init__(self, histo=None, rootFile=None, upperSplit=None, lowerSplit=None, histoName=None, hinfo=None, pinfo=None, aux = 1.0):

        self.xDim = {}
        self.yDim = {}

        self.info = copy.deepcopy(pinfo)
        self.info.update(hinfo)

        self.histoName = histoName
        self.filePath  = rootFile

        # Custom tuned values for these dimensions/sizes
        # Are scaled accordingly if margins change
        self.xDim["title"]  = 0.050; self.yDim["title"]  = 0.050
        self.xDim["label"]  = 0.039; self.yDim["label"]  = 0.039
        self.xDim["offset"] = 1.1;   self.yDim["offset"] = 1.7 

        self.histogram = histo

        self.setupHisto(upperSplit, lowerSplit, aux)
        
    # Take input hist as a background hist and assume self is signal
    # Calculate simple S / sqrt(B) significance in this scenario
    def Significance(self, otherHist):

        significance = 0.0
        for xBin in range(1,otherHist.GetNbinsX()+1):

            S = self.histogram.GetBinContent(xBin)
            B = otherHist.GetBinContent(xBin)

            if B > 1.0 and S > 1.0:
                significance += (S / math.sqrt(B + (0.3*B)**2.0))**2.0

        return significance**0.5

    def Integral(self):
        return self.histogram.Integral()
        
    def IntegralError(self):
        error = ctypes.c_double(-999)
        integral = self.histogram.IntegralAndError(0, self.histogram.GetNbinsX(), error)
        return error.value

    def Clone(self, name):
        return self.histogram.Clone(name)

    def Scale(self, fraction):
        return self.histogram.Scale(fraction)

    def IsGood(self):
        return self.histogram != -1
    
    def Write(self):
        self.histogram.Write()

    def Draw(self, canvas, showNevents, firstDraw, nLegend, legend):

        if self.histogram != -1:
            lname = self.info["name"]
            if showNevents:
                lname += " (%.1f)"%(self.histogram.Integral())

            if not firstDraw:
                self.histogram.Draw(self.info["option"])
                firstDraw = True
            else:
                self.histogram.Draw("%s SAME"%(self.info["option"]))
    
            legend.AddEntry(self.histogram, lname, self.info["loption"])
            nLegend += 1     

        return nLegend, firstDraw

    # Simply get the raw histogram from the input ROOT file
    def getHisto(self):

        if self.histogram != None:
            return 0

        f = None
        try:
            f = ROOT.TFile.Open(self.filePath, "READ")
        except Exception as e:
            print("WARNING: Could not open file \"%s\" !"%(self.filePath), e)
            return -1
        
        if f == None:
            print("WARNING: Skipping file \"%s\" for histo \"%s\" !"%(self.filePath, self.histoName))
            return -1

        histo = f.Get(self.histoName)

        if histo == None:
            print("WARNING: Skipping histo \"%s\" from \"%s\""%(self.histoName, self.filePath))
            return -1

        histo.SetDirectory(0)

        f.Close()

        self.histogram = histo
        return 0

    def setupHisto(self, upperSplit, lowerSplit, aux):

        code = self.getHisto()

        if code != -1:
            self.histogram.GetXaxis().SetTitleSize(self.xDim["title"] / upperSplit); self.histogram.GetYaxis().SetTitleSize(self.yDim["title"] / upperSplit)
            self.histogram.GetXaxis().SetLabelSize(self.xDim["label"] / upperSplit); self.histogram.GetYaxis().SetLabelSize(self.yDim["label"] / upperSplit)
            self.histogram.GetXaxis().SetTitleOffset(self.xDim["offset"]);           self.histogram.GetYaxis().SetTitleOffset(self.yDim["offset"] * upperSplit * aux)
            self.histogram.GetXaxis().SetTitle(self.info["X"]["title"]);             self.histogram.GetYaxis().SetTitle(self.info["Y"]["title"])

            if "lcolor" in self.info:
                self.histogram.SetLineColor(self.info["lcolor"])
                self.histogram.SetMarkerColor(self.info["lcolor"])
            else:
                self.histogram.SetLineColor(self.info["color"])
                self.histogram.SetMarkerColor(self.info["color"])

            self.histogram.SetMarkerSize(self.info["msize"]); self.histogram.SetMarkerStyle(self.info["mstyle"])
            self.histogram.SetLineWidth(self.info["lsize"]);  self.histogram.SetLineStyle(self.info["lstyle"])
    
            if "loption" in self.info and "L" not in self.info["loption"]:
                self.histogram.SetFillColorAlpha(self.info["color"], 1.0)

            if "fstyle" in self.info:
                self.histogram.SetFillStyle(self.info["fstyle"])

            self.histogram.RebinX(self.info["X"]["rebin"])
            self.histogram.GetXaxis().SetRangeUser(self.info["X"]["min"], self.info["X"]["max"])

            self.histogram.SetTitle("")
        else:
            self.histogram = -1

# The StackPlotter class oversees the creation of all stack plots
#     approved     : are these plots approved, then no "Preliminary"
#     noRatio      : do not make a ratio plot in bottom panel
#     printNEvents : show the number of events in legend
#     printInfo    : show significance and cut label
#     year         : corresponding year for the plots/inputs
#     outpath      : where to put the plots, path is created if missing
#     inpath       : where the input ROOT files are located
#     normMC2Data  : normalize MC to the data
#     normalize    : normalize all processes to unity area
#     histograms   : dictionary containing config info for desired histos
#     selections   : list of cut strings appearing in the names of plots
#     backgrounds  : dictionary containing config info for desired backgrounds
#     signals      : dictionary containing config info for desired signals
#     data         : dictionary containing config info for data
class StackPlotter:

    def __init__(self, approved, noRatio, printNEvents, printInfo, year, outpath, inpath, normMC2Data, normalize, histograms, selections, backgrounds, signals, data):

        self.histograms  = histograms
        self.selections  = selections
        self.backgrounds = backgrounds
        self.signals     = signals
        self.data        = data

        self.year        = year
        self.inpath      = inpath
        self.outpath     = outpath
        
        self.njetsTableDictionary = {}

        if not os.path.exists(self.outpath):
            os.makedirs(self.outpath)

        self.approved     = approved
        self.noRatio      = noRatio
        self.normMC2Data  = normMC2Data
        self.normalize    = normalize
        self.printNEvents = printNEvents
        self.printInfo    = printInfo

        # Customized numbers that are scaled
        # to work with or without a ratio plot
        self.TopMargin    = 0.06
        self.BottomMargin = 0.12
        self.RightMargin  = 0.04
        self.LeftMargin   = 0.16

    # Draw the significance and cut labels on the plot
    def add_Significance_CutLabel(self, canvas, significance_550, selections):

        mark = ROOT.TLatex()
        mark.SetNDC(True)

        mark.SetTextAlign(11)
        mark.SetTextFont(52); mark.SetTextSize(0.025)
        mark.DrawLatex(0.2, 0.92, "Significance = %.3f"%(significance_550))
        mark.DrawLatex(0.45, 0.92, "%s"%(selections))

    # Create a canvas and determine if it should be split for a ratio plot
    # Margins are scaled on-the-fly so that distances are the same in either
    # scenario.
    def makeCanvas(self, doLogY):

        randStr = ''.join([random.choice(string.ascii_letters + string.digits) for n in xrange(8)])
       
        canvas = ROOT.TCanvas(randStr, randStr, 900, 900)

        # Split the canvas 70 / 30 by default if doing ratio
        # scale parameter keeps text sizes in ratio panel the
        # same as in the upper panel
        split           = 0.3
        self.upperSplit = 1.0
        self.lowerSplit = 1.0
        self.scale      = 1.0

        if not self.noRatio:
            self.upperSplit = 1.0-split
            self.lowerSplit = split
            self.scale = self.upperSplit / self.lowerSplit

            canvas.Divide(1,2)

            canvas.cd(1)
            ROOT.gPad.SetPad(0.0, split, 1.0, 1.0)
            ROOT.gPad.SetTopMargin(self.TopMargin / self.upperSplit)
            ROOT.gPad.SetBottomMargin(0)
            ROOT.gPad.SetLeftMargin(self.LeftMargin)
            ROOT.gPad.SetRightMargin(self.RightMargin)
            if doLogY:
                ROOT.gPad.SetLogy()
   
            canvas.cd(2)
            ROOT.gPad.SetPad(0.0, 0.0, 1.0, split)
            ROOT.gPad.SetTopMargin(0)
            ROOT.gPad.SetBottomMargin(self.BottomMargin / self.lowerSplit)
            ROOT.gPad.SetLeftMargin(self.LeftMargin)
            ROOT.gPad.SetRightMargin(self.RightMargin)

        else:
            ROOT.gPad.SetTopMargin(self.TopMargin)
            ROOT.gPad.SetBottomMargin(self.BottomMargin)
            ROOT.gPad.SetLeftMargin(self.LeftMargin)
            ROOT.gPad.SetRightMargin(self.RightMargin)
            if doLogY:
                ROOT.gPad.SetLogy()

        return canvas

    # Some customized sizes and distances but scaled
    # when going between ratio plot and pure stack with no ratio
    def makeLegends(self, nBkgs, nSigs, doLogY, theMin, theMax):

        textSize = 0.025 / self.upperSplit
        space    = 0.015

        bkgXmin = 0.70 - self.RightMargin
        if self.printNEvents:
            bkgXmin = 0.50 - self.RightMargin
        if tooManyBkgds:
            bkgXmin = self.LeftMargin + 0.05
            
        bkgYmax = 1.0-(self.TopMargin/self.upperSplit)-0.01
        bkgXmax = 1.0-self.RightMargin-0.01
        bkgYmin = bkgYmax-nBkgs*(textSize+space)
        
        if self.printInfo:
            bkgYmax -= 0.025

        bkgYFrac = (1.0-self.TopMargin-bkgYmin) / (1.0 - self.TopMargin - self.BottomMargin)

        bkgLegend = ROOT.TLegend(bkgXmin, bkgYmin, bkgXmax, bkgYmax)
        bkgLegend.SetBorderSize(0)
        bkgLegend.SetTextSize(textSize)

        sigXmin = self.LeftMargin+0.03
        sigYmax = bkgYmax
        sigXmax = bkgXmin
        sigYmin = bkgYmax-nSigs*(textSize+space) 

        if self.printInfo:
            bkgYmax -= 0.025

        sigYFrac = (1.0-self.TopMargin-sigYmin) / (1.0 - self.TopMargin - self.BottomMargin)

        sigLegend = ROOT.TLegend(sigXmin, sigYmin, sigXmax, sigYmax)
        sigLegend.SetBorderSize(0)
        sigLegend.SetMargin(0.10)
        sigLegend.SetTextSize(textSize)

        yMax = 1.0; factor = 1.0; power = 1.0
        if doLogY and theMax != 0.0 and theMin != 0.0:
            power = math.log10(theMax / theMin) * 5.0

        theFrac = bkgYFrac if bkgYFrac > sigYFrac else sigYFrac

        if self.printInfo:
            yMax = (theMax-theMin) * (1.1 - theFrac)**(-power) * factor
        else:                             
            yMax = (theMax-theMin) * (1.0 - theFrac)**(-power) * factor

        return bkgLegend, sigLegend, yMax

    def addCMSlogo(self, canvas):

        canvas.cd()

        mark = ROOT.TLatex()
        mark.SetNDC(True)

        mark.SetTextAlign(11)
        mark.SetTextSize(0.055)
        mark.SetTextFont(61)
        mark.DrawLatex(self.LeftMargin, 1 - (self.TopMargin - 0.015), "CMS")

        mark.SetTextFont(52)
        mark.SetTextSize(0.040)

        simStr = ""
        if self.data == {}:
            simStr = "Simulation "

        if self.approved:
            mark.DrawLatex(self.LeftMargin + 0.12, 1 - (self.TopMargin - 0.017), "%sSupplementary"%(simStr))
        else:
            mark.DrawLatex(self.LeftMargin + 0.12, 1 - (self.TopMargin - 0.017), "%sPreliminary"%(simStr))

        mark.SetTextFont(42)
        mark.SetTextAlign(31)
        if "Run 2" in self.year:
            mark.DrawLatex(1 - self.RightMargin, 1 - (self.TopMargin - 0.017), "Run 2 (13 TeV)")
        else:
            mark.DrawLatex(1 - self.RightMargin, 1 - (self.TopMargin - 0.017), "%s (13 TeV)"%(self.year))

    # Main function to compose the full stack plot with or without a ratio panel
    def makePlots(self):

        # Top loop begins going over each histo-to-be-plotted
        for hname, hinfo in self.histograms.items():

            orders = [-1]; selections = ["None"]
            if "orders" in hinfo:
                orders = hinfo["orders"]

            if "selections" in hinfo:
                selections = hinfo["selections"]

            newName = hname
            for order in orders:
                for selection in self.selections:
         
                    newName = hname.replace("@", "%d"%(order)).replace("?", "%s"%(selection))
                    newInfo = copy.deepcopy(hinfo)
                    newInfo["X"]["title"] = hinfo["X"]["title"].replace("@", "%d"%(order))

                    canvas               = self.makeCanvas(newInfo["logY"])
                    if self.noRatio:
                        canvas.cd()
                    else:
                        canvas.cd(1)

                    firstDraw = False
                    bstack = ROOT.THStack("stack%s"%(newName), "stack%s"%(newName))
                    dummy = None
                    totalMC = None
                    bhistos = {}
    
                    nBkgLegend = 0
                    nSigLegend = 0
                    theSignificance_350 = 0.0
                    theSignificance_550 = 0.0
                    theSignificance_850 = 0.0

                    dataScale = 0.0
                    mcScale = 0.0
                    theMax = 0.0

                    # Preemptively get data counts to be used for normalzing the histograms later
                    for dname, dinfo in self.data.items():

                        rootFile = "%s/%s_%s.root"%(self.inpath, self.year, dname)

                        Hobj = Histogram(None, rootFile, self.upperSplit, self.lowerSplit, newName, newInfo, dinfo)
                        if Hobj.IsGood():
                            dataScale = Hobj.Integral()

                            if self.normalize:
                                Hobj.Scale(1.0 / dataScale)

                            tempMax = Hobj.histogram.GetMaximum()
                            if tempMax > theMax:
                                theMax = tempMax


                    # Preemptively get MC counts to be used for normalizing the histograms later
                    totalWegBkg = 0.0; wegTT    = 0.0; wegQCD = 0.0; ttFrac = 0.0; qcdFrac = 0.0                     
                    unWegTT     = 0.0; unWegQCD = 0.0

                    for bname, binfo in self.backgrounds.items(): 

                        rootFile = "%s/%s_%s.root"%(self.inpath, self.year, bname)
    
                        Hobj = Histogram(None, rootFile, self.upperSplit, self.lowerSplit, newName, newInfo, binfo)
                        
                        if Hobj.IsGood(): mcScale += Hobj.Integral()

                        # get fractions and unweighted events for TT and QCD
                        totalWegBkg += Hobj.Integral()
            
                        if bname == "TT":
                            wegTT   = Hobj.Integral()
                            unWegTT = Hobj.histogram.GetEntries()                        

                        if bname == "QCD":
                            wegQCD   = Hobj.Integral()
                            unWegQCD = Hobj.histogram.GetEntries()

                    if totalWegBkg != 0.0:
                        ttFrac  = wegTT / totalWegBkg
                        qcdFrac = wegQCD / totalWegBkg
          
                    if "njets" in hname:
                        self.njetsTableDictionary[selection] = {}
                        self.njetsTableDictionary[selection]["TT Frac"]               = ttFrac
                        self.njetsTableDictionary[selection]["QCD Frac"]              = qcdFrac
                        self.njetsTableDictionary[selection]["TT Unweighted Events"]  = unWegTT
                        self.njetsTableDictionary[selection]["QCD Unweighted Events"] = unWegQCD

                    # Preemptively loop over signal to determine maximums
                    wegRPV350   = 0.0; wegRPV550   = 0.0; wegRPV850   = 0.0; rpv350Frac = 0.0; rpv550Frac = 0.0; rpv850Frac = 0.0;
                    unWegRPV350 = 0.0; unWegRPV550 = 0.0; unWegRPV850 = 0.0;

                    for sname, sinfo in self.signals.items():

                        rootFile = "%s/%s_%s.root"%(self.inpath, self.year, sname)

                        Hobj = Histogram(None, rootFile, self.upperSplit, self.lowerSplit, newName, newInfo, sinfo)
                        if Hobj.IsGood():
                            sigScale = Hobj.Integral()

                            if self.normalize:
                                Hobj.Scale(1.0 / sigScale)

                            tempMax = Hobj.histogram.GetMaximum()
                            if tempMax > theMax:
                                theMax = tempMax

                            # get fractions and unweighted events for RPV 350, 550, 850
                            if "350" in sname:
                                wegRPV350   = Hobj.Integral()
                                unWegRPV350 = Hobj.histogram.GetEntries()

                            if "550" in sname:
                                wegRPV550   = Hobj.Integral()
                                unWegRPV550 = Hobj.histogram.GetEntries()

                            if "850" in sname:
                                wegRPV850   = Hobj.Integral()
                                unWegRPV850 = Hobj.histogram.GetEntries()

                    if wegRPV350 + totalWegBkg != 0.0: rpv350Frac = wegRPV350 / (wegRPV350 + totalWegBkg) 
                    if wegRPV550 + totalWegBkg != 0.0: rpv550Frac = wegRPV550 / (wegRPV550 + totalWegBkg)
                    if wegRPV850 + totalWegBkg != 0.0: rpv850Frac = wegRPV850 / (wegRPV850 + totalWegBkg)
                   
                    if "njets" in hname: 
                        self.njetsTableDictionary[selection]["RPV_2t6j_mStop-350 Frac"]              = rpv350Frac
                        self.njetsTableDictionary[selection]["RPV_2t6j_mStop-550 Frac"]              = rpv550Frac
                        self.njetsTableDictionary[selection]["RPV_2t6j_mStop-850 Frac"]              = rpv850Frac
                        self.njetsTableDictionary[selection]["RPV_2t6j_mStop-350 Unweighted Events"] = unWegRPV350
                        self.njetsTableDictionary[selection]["RPV_2t6j_mStop-550 Unweighted Events"] = unWegRPV550
                        self.njetsTableDictionary[selection]["RPV_2t6j_mStop-850 Unweighted Events"] = unWegRPV850

                    # Loop over each background and get their respective histo, scale if necessary
                    option = "HIST"; loption = "F"
                    for bname, binfo in self.backgrounds.items(): 

                        rootFile = "%s/%s_%s.root"%(self.inpath, self.year, bname)
    
                        if "option"  not in binfo: binfo["option"]  = option
                        if "loption" not in binfo: binfo["loption"] = loption

                        Hobj = Histogram(None, rootFile, self.upperSplit, self.lowerSplit, newName, newInfo, binfo)

                        if Hobj.IsGood(): 
                            if mcScale != 0.0:
                                if self.normMC2Data:
                                    Hobj.Scale(dataScale / mcScale)
                                elif self.normalize:
                                    Hobj.Scale(1.0 / mcScale)
                            if not firstDraw:
                                totalMC = Hobj.Clone("totalMC%s"%(newName))
                                firstDraw = True
                            else:
                                totalMC.Add(Hobj.histogram)

                            bhistos[Hobj.Integral()] = (binfo["name"], Hobj.histogram, binfo["option"], binfo["loption"])
                            dummy = Hobj.Clone("dummy%s"%(hname)); dummy.Reset("ICESM")

                    # Add each background histo to the stack in order based on number of entries
                    for count, h in sorted(bhistos.items(), key=lambda x: x[0], reverse=False): 
                        bstack.Add(h[1], h[2])
                        nBkgLegend += 1

                    tempMax = bstack.GetMaximum()
                    if tempMax > theMax:
                        theMax = tempMax

                    theMin = newInfo["Y"]["min"] if "min" in newInfo["Y"] else 0.0

                    if self.normalize and mcScale != 0.0:
                        theMin /= mcScale

                    # Here we get the bkgd and sig legends as well as a tuned maximum for the canvas to avoid overlap
                    bkgLegend, sigLegend, yMax = self.makeLegends(len(self.backgrounds), len(self.signals), newInfo["logY"], theMin, theMax)

                    for count, h in sorted(bhistos.items(), key=lambda x: x[0], reverse=True):
                        lname = h[0]
                        if self.printNEvents:
                            lname += " (%.1f)"%(h[1].Integral())

                        bkgLegend.AddEntry(h[1], lname, h[3])

                    dummy.SetMaximum(yMax)
                    dummy.SetMinimum(theMin)
                    dummy.Draw()
                    bstack.Draw("SAME")

                    # Loop over each signal and get their respective histo
                    option = "HIST"; loption = "L"
                    for sname, sinfo in self.signals.items(): 

                        rootFile = "%s/%s_%s.root"%(self.inpath, self.year, sname)

                        if "option"  not in sinfo: sinfo["option"]  = option 
                        if "loption" not in sinfo: sinfo["loption"] = loption

                        Hobj = Histogram(None, rootFile, self.upperSplit, self.lowerSplit, newName, newInfo, sinfo)

                        if "350" in sname or len(self.signals) == 1:
                            theSignificance_350 = Hobj.Significance(totalMC)
                       
                        elif "550" in sname or len(self.signals) == 1:
                            theSignificance_550 = Hobj.Significance(totalMC)
                
                        elif "850" in sname or len(self.signals) == 1:
                            theSignificance_850 = Hobj.Significance(totalMC)
 
                        scale = Hobj.Integral()
                        if self.normalize and scale != 0.0:
                            Hobj.Scale(1.0 / scale)

                        nSigLegend, firstDraw = Hobj.Draw(canvas, self.printNEvents, firstDraw, nSigLegend, sigLegend)

                        if "njets" in hname:
                            self.njetsTableDictionary[selection]["significance_350"] = theSignificance_350
                            self.njetsTableDictionary[selection]["significance_550"] = theSignificance_550
                            self.njetsTableDictionary[selection]["significance_850"] = theSignificance_850

                    # Loop over the data and get their respective histo
                    option = "E0P"; loption = "ELP"
                    for dname, dinfo in self.data.items():

                        rootFile = "%s/%s_%s.root"%(self.inpath, self.year, dname)

                        if "option"  not in dinfo: dinfo["option"]  = option 
                        if "loption" not in dinfo: dinfo["loption"] = loption

                        Hobj = Histogram(None, rootFile, self.upperSplit, self.lowerSplit, newName, newInfo, dinfo)

                        scale = Hobj.Integral()
                        if self.normalize and scale != 0.0:
                            Hobj.Scale(1.0 / scale)

                        nBkgLegend, firstDraw = Hobj.Draw(canvas, self.printNEvents, firstDraw, nBkgLegend, bkgLegend)

                        if Hobj.IsGood():

                            ratio = Hobj.Clone("ratio")
                            ratio.Divide(totalMC)

                            for xbin in range(1, ratio.GetNbinsX()+1):
                                if totalMC.GetBinContent(xbin) == 0.0 or Hobj.histogram.GetBinContent(xbin) == 0.0:
                                    ratio.SetBinContent(xbin, -999.0)


                    if nBkgLegend != 0: bkgLegend.Draw("SAME")
                    if nSigLegend != 0: sigLegend.Draw("SAME")

                    self.addCMSlogo(canvas)

                    # put the siginificance and cut label on the canvas
                    if self.printInfo:
                        self.add_Significance_CutLabel(canvas, theSignificance_350, selection)

                    # Here we go into bottom panel if drawing the ratio
                    if not self.noRatio:

                        rinfo = {"name" : "ratio", "color" : ROOT.kBlack,  "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 1 / self.upperSplit}
                        rnewInfo = copy.deepcopy(newInfo)
                        rnewInfo["Y"]["title"] = "#frac{Data}{Pred.}"
                        rnewInfo["X"]["rebin"] = 1

                        canvas.cd(2)
                        ROOT.gPad.SetGridy()
                        ratio = Histogram(ratio, None, self.upperSplit/self.scale, self.lowerSplit/self.scale, None, rnewInfo, rinfo, 0.8).histogram

                        ratio.SetMinimum(0.0); ratio.SetMaximum(2.2)
                        ratio.GetYaxis().SetNdivisions(5, 5, 0)

                        ratio.Draw("E0P")

                    canvas.SaveAs("%s/%s_%s.png"%(self.outpath, self.year, newName))


    # function to put the raw nEvents significance and njet bin by njet bin bkg & sig fractions for baseline selection
    def make_njetsTable(self):        

        f = open("%s/%s_nRawEvents_significance_bkgSigFracs_njets_stackPlots_0l_.tex" %(self.outpath, self.year), "w")
        f.write("\\resizebox{\linewidth}{!}{%")
        f.write("\n")
        f.write("    \def\\arraystretch{0.6}")
        f.write("\n")
        f.write("    \\begin{tabular}{| l | c | c | c | c | c | c | c | c | c | c | c | c | c |}")
        f.write("\n")
        f.write("        \hline")
        f.write("\n")
        f.write("        & \multicolumn{2}{c|}{\\textcolor{ttjetscol}{TT}} & \multicolumn{2}{c|}{\\textcolor{qcdcol}{QCD}} & \multicolumn{2}{c|}{\\textcolor{rpvcol}{RPV350}} & \multicolumn{2}{c|}{\\textcolor{rpvcol}{RPV550}} & \multicolumn{2}{c|}{\\textcolor{rpvcol}{RPV850}} & \multicolumn{3}{c|}{Significance} \\\\")
        f.write("\n")
        f.write("        \hline")
        f.write("\n")
        f.write("        Cuts & \\textcolor{ttjetscol}{Raw Events} & \\textcolor{ttjetscol}{Frac} & \\textcolor{qcdcol}{Raw Events} & \\textcolor{qcdcol}{Frac} & \\textcolor{rpvcol}{Raw Events}    & \\textcolor{rpvcol}{Frac} & \\textcolor{rpvcol}{Raw Events} & \\textcolor{rpvcol}{Frac} & \\textcolor{rpvcol}{Raw Events} & \\textcolor{rpvcol}{Frac} & RPV350 & RPV550 & RPV850 \\\\")
        f.write("\n")
        f.write("        \hline")
        f.write("\n")
 
        # loops for putting the table
        for cut in self.selections:  
 
            f.write("    \scriptsize %s & \scriptsize \\textcolor{ttjetscol}{%d} & \scriptsize \\textcolor{ttjetscol}{%.3f} & \scriptsize \\textcolor{qcdcol}{%d} & \scriptsize \\textcolor{qcdcol}{%.3f} & \scriptsize \\textcolor{rpvcol}{%d} & \scriptsize \\textcolor{rpvcol}{%.3f} & \scriptsize \\textcolor{rpvcol}{%d} & \scriptsize \\textcolor{rpvcol}{%.3f} & \scriptsize \\textcolor{rpvcol}{%d} & \scriptsize \\textcolor{rpvcol}{%.3f} & \scriptsize %.3f & \scriptsize %.3f & \scriptsize %.3f \\\\"
                %(cut.replace("_", "\_"), self.njetsTableDictionary[cut]["TT Unweighted Events"],  self.njetsTableDictionary[cut]["TT Frac"], 
                       self.njetsTableDictionary[cut]["QCD Unweighted Events"], self.njetsTableDictionary[cut]["QCD Frac"],
                       self.njetsTableDictionary[cut]["RPV_2t6j_mStop-350 Unweighted Events"], self.njetsTableDictionary[cut]["RPV_2t6j_mStop-350 Frac"], 
                       self.njetsTableDictionary[cut]["RPV_2t6j_mStop-550 Unweighted Events"], self.njetsTableDictionary[cut]["RPV_2t6j_mStop-550 Frac"],
                       self.njetsTableDictionary[cut]["RPV_2t6j_mStop-850 Unweighted Events"], self.njetsTableDictionary[cut]["RPV_2t6j_mStop-850 Frac"],
                       self.njetsTableDictionary[cut]["significance_350"], self.njetsTableDictionary[cut]["significance_550"], self.njetsTableDictionary[cut]["significance_850"]))
            f.write("\n")
            f.write("        \hline")
            f.write("\n")

        f.write("    \end{tabular}")
        f.write("\n")
        f.write("}")
        f.write("\n")
        f.close()                             


if __name__ == "__main__":

    usage = "usage: %stackPlotter [options]"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("--noRatio",      dest="noRatio",      help="No ratio plot",               default=False,  action="store_true") 
    parser.add_argument("--approved",     dest="approved",     help="Plot is approved",            default=False,  action="store_true") 
    parser.add_argument("--printNEvents", dest="printNEvents", help="Show number of events",       default=False,  action="store_true") 
    parser.add_argument("--normMC2Data",  dest="normMC2Data",  help="Normalize MC to data",        default=False,  action="store_true") 
    parser.add_argument("--normalize",    dest="normalize",    help="Normalize all to unity",      default=False,  action="store_true") 
    parser.add_argument("--printInfo",    dest="printInfo",    help="Print significance and cuts", default=False,  action="store_true")
    parser.add_argument("--makeTable",    dest="makeTable",    help="Make the table of fracs",     default=False,  action="store_true")
    parser.add_argument("--inpath",       dest="inpath",       help="Path to root files",          default="NULL", required=True)
    parser.add_argument("--outpath",      dest="outpath",      help="Where to put plots",          default="NULL", required=True)
    parser.add_argument("--year",         dest="year",         help="which year",                                  required=True)
    parser.add_argument("--options",      dest="options",      help="options file",                default="stackPlotter_aux", type=str)
    args = parser.parse_args()

    # The auxiliary file contains many "hardcoded" items
    # describing which histograms to get and how to draw
    # them. These things are changed often by the user
    # and thus are kept in separate sidecar file.
    importedGoods = __import__(args.options)

    # Names of histograms, rebinning, titles, ranges, etc.
    histograms  = importedGoods.histograms

    # Background/signal/data categories to plot, and plotting options
    selections  = importedGoods.selections
    backgrounds = importedGoods.backgrounds
    signals     = importedGoods.signals
    data        = importedGoods.data

    plotter = StackPlotter(args.approved, args.noRatio, args.printNEvents, args.printInfo, args.year, args.outpath, args.inpath, args.normMC2Data, args.normalize, histograms, selections, backgrounds, signals, data)
    plotter.makePlots()

    if args.makeTable:
        plotter.make_njetsTable()
        plotter.makeTable()
