#! /usr/bin/env python

import ROOT, random, os, argparse, string, copy, math
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

        self.info = hinfo
        self.info.update(pinfo)

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

            if B > 5.0 and S > 5.0:
                significance += (S / math.sqrt(B + (0.3*B)**2.0))**2.0

        return significance**0.5

    def Integral(self):
        return self.histogram.Integral()
        
    def Clone(self, name):
        return self.histogram.Clone(name)

    def Scale(self, fraction):
        return self.histogram.Scale(fraction)

    def IsGood(self):
        return self.histogram != -1

    def Draw(self, canvas, showNevents, firstDraw, nLegend, legend, histOpt, legOpt):

        if self.histogram != -1:
            lname = self.info["name"]
            if showNevents:
                lname += " (%.1f)"%(self.histogram.Integral())

            if not firstDraw:
                self.histogram.Draw(histOpt)
                firstDraw = True
            else:
                self.histogram.Draw("%s SAME"%(histOpt))
    
            legend.AddEntry(self.histogram, lname, legOpt)
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

            self.histogram.SetMarkerColor(self.info["color"]); self.histogram.SetMarkerSize(self.info["msize"]); self.histogram.SetMarkerStyle(self.info["mstyle"])
            self.histogram.SetLineColor(self.info["color"]);   self.histogram.SetLineWidth(self.info["lsize"]);  self.histogram.SetLineStyle(self.info["lstyle"])
    
            # Only fill self.histogram with color when no line is present
            if self.info["lsize"] == 0:
                self.histogram.SetFillColor(self.info["color"])

            self.histogram.RebinX(self.info["X"]["rebin"])
            self.histogram.GetXaxis().SetRangeUser(self.info["X"]["min"], self.info["X"]["max"])

            self.histogram.SetTitle("")
        else:
            self.histogram = -1

# The StackPlotter class oversees the creation of all stack plots
#     approved    : are these plots approved, then no "Preliminary"
#     noRatio           : do not make a ratio plot in bottom panel
#     printNEvents      : show the number of events in legend
#     printSignificance : show a simple S / sqrt(B) significance
#     year              : corresponding year for the plots/inputs
#     outpath           : where to put the plots, path is created if missing
#     inpath            : where the input ROOT files are located
#     normMC2Data       : normalize MC to the data
#     normalize         : normalize all processes to unity area
#     histograms        : dictionary containing config info for desired histos
#     selections        : list of cut strings appearing in the names of plots
#     backgrounds       : dictionary containing config info for desired backgrounds
#     signals           : dictionary containing config info for desired signals
#     data              : dictionary containing config info for data
class StackPlotter:

    def __init__(self, approved, noRatio, printNEvents, printSignificance, year, outpath, inpath, normMC2Data, normalize, histograms, selections, backgrounds, signals, data):

        self.histograms  = histograms
        self.selections  = selections
        self.backgrounds = backgrounds
        self.signals     = signals
        self.data        = data

        self.year        = year
        self.inpath      = inpath
        self.outpath     = outpath

        if not os.path.exists(self.outpath):
            os.makedirs(self.outpath)

        self.approved          = approved
        self.noRatio           = noRatio
        self.normMC2Data       = normMC2Data
        self.normalize         = normalize
        self.printNEvents      = printNEvents
        self.printSignificance = printSignificance

        # Customized numbers that are scaled
        # to work with or without a ratio plot
        self.TopMargin    = 0.06
        self.BottomMargin = 0.12
        self.RightMargin  = 0.04
        self.LeftMargin   = 0.16

    # Draw the significance value on the plot
    def addSignificance(self, canvas, sign, legBottom):

        mark = ROOT.TLatex()
        mark.SetNDC(True)

        mark.SetTextAlign(11)
        mark.SetTextFont(52); mark.SetTextSize(0.025)
        mark.DrawLatex(self.LeftMargin + 0.05, legBottom / (1.0-self.TopMargin) - 0.01, "S / #sqrt{B} = %.2f"%(sign))

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

        bkgXmin = 0.60 if self.printNEvents else 0.75
        bkgYmax = 1.0-(self.TopMargin/self.upperSplit)-0.01
        bkgXmax = 1.0-self.RightMargin-0.01
        bkgYmin = bkgYmax-nBkgs*(textSize+space)

        bkgYFrac = (1.0-self.TopMargin-bkgYmin) / (1.0 - self.TopMargin - self.BottomMargin)

        bkgLegend = ROOT.TLegend(bkgXmin, bkgYmin, bkgXmax, bkgYmax)
        bkgLegend.SetBorderSize(0)
        bkgLegend.SetTextSize(textSize)

        sigXmin = self.LeftMargin+0.03
        sigYmax = bkgYmax
        sigXmax = bkgXmin
        sigYmin = bkgYmax-nSigs*(textSize+space) 

        sigYFrac = (1.0-self.TopMargin-sigYmin) / (1.0 - self.TopMargin - self.BottomMargin)

        sigLegend = ROOT.TLegend(sigXmin, sigYmin, sigXmax, sigYmax)
        sigLegend.SetBorderSize(0)
        sigLegend.SetMargin(0.10)
        sigLegend.SetTextSize(textSize)

        yMax = 1.0; factor = 1.05; power = 1.0
        if doLogY and theMax != 0.0 and theMin != 0.0:
            power = math.log10(theMax / theMin) * 3.0

        theFrac = bkgYFrac if bkgYFrac > sigYFrac else sigYFrac

        if self.printSignificance:
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
                    theSignificance = 0.0

                    dataScale = 0.0
                    mcScale = 0.0
                    # Preemptively get data counts
                    for dname, dinfo in self.data.items():

                        rootFile = "%s/%s_%s.root"%(self.inpath, self.year, dname)

                        Hobj = Histogram(None, rootFile, self.upperSplit, self.lowerSplit, newName, newInfo, dinfo)
                        if Hobj.IsGood(): dataScale = Hobj.Integral()

                    # Preemptively get MC counts
                    for bname, binfo in self.backgrounds.items(): 

                        rootFile = "%s/%s_%s.root"%(self.inpath, self.year, bname)
    
                        Hobj = Histogram(None, rootFile, self.upperSplit, self.lowerSplit, newName, newInfo, binfo)
                        
                        if Hobj.IsGood(): mcScale += Hobj.Integral()

                    # Loop over each background and get their respective histo 
                    for bname, binfo in self.backgrounds.items(): 

                        rootFile = "%s/%s_%s.root"%(self.inpath, self.year, bname)
    
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

                            bhistos[Hobj.Integral()] = (binfo["name"], Hobj.histogram)
                            dummy = Hobj.Clone("dummy%s"%(hname)); dummy.Reset("ICESM")
                       
                    for count, h in sorted(bhistos.items(), key=lambda x: x[0], reverse=False): 
                        bstack.Add(h[1], "HIST")
                        nBkgLegend += 1

                    theMax = bstack.GetMaximum()

                    bkgLegend, sigLegend, yMax = self.makeLegends(len(self.backgrounds), len(self.signals), newInfo["logY"], newInfo["Y"]["min"], theMax)

                    for count, h in sorted(bhistos.items(), key=lambda x: x[0], reverse=True):
                        lname = h[0]
                        if self.printNEvents:
                            lname += " (%.1f)"%(h[1].Integral())

                        bkgLegend.AddEntry(h[1], lname, "F")

                    dummy.SetMaximum(yMax); dummy.SetMinimum(hinfo["Y"]["min"]) if "min" in hinfo["Y"] else dummy.SetMinimum(0.2)
                    dummy.Draw()
                    bstack.Draw("SAME")

                    # Loop over each signal and get their respective histo
                    for sname, sinfo in self.signals.items(): 

                        rootFile = "%s/%s_%s.root"%(self.inpath, self.year, sname)

                        Hobj = Histogram(None, rootFile, self.upperSplit, self.lowerSplit, newName, newInfo, sinfo)

                        if "550" in sname or len(self.signals) == 1:
                            theSignificance = Hobj.Significance(totalMC)
                        
                        scale = Hobj.Integral()
                        if self.normalize and scale != 0.0:
                            Hobj.Scale(1.0 / scale)

                        nSigLegend, firstDraw = Hobj.Draw(canvas, self.printNEvents, firstDraw, nSigLegend, sigLegend, "HIST", "L")

                    # Loop over the data and get their respective histo
                    for dname, dinfo in self.data.items():

                        rootFile = "%s/%s_%s.root"%(self.inpath, self.year, dname)

                        Hobj = Histogram(None, rootFile, self.upperSplit, self.lowerSplit, newName, newInfo, dinfo)

                        scale = Hobj.Integral()
                        if self.normalize and scale != 0.0:
                            Hobj.Scale(1.0 / scale)

                        nBkgLegend, firstDraw = Hobj.Draw(canvas, self.printNEvents, firstDraw, nBkgLegend, bkgLegend, "E0P", "ELP")

                        if Hobj.IsGood():

                            ratio = Hobj.Clone("ratio")
                            ratio.Divide(totalMC)

                            for xbin in range(1, ratio.GetNbinsX()+1):
                                if totalMC.GetBinContent(xbin) == 0.0 or Hobj.histogram.GetBinContent(xbin) == 0.0:
                                    ratio.SetBinContent(xbin, -999.0)

                    if nBkgLegend != 0: bkgLegend.Draw("SAME")
                    if nSigLegend != 0: sigLegend.Draw("SAME")

                    self.addCMSlogo(canvas)
                    if self.printSignificance:
                        self.addSignificance(canvas, theSignificance, sigLegend.GetY1())

                    # Here we go into bottom panel if drawing the ratio
                    if not self.noRatio:

                        rinfo = {"name" : "ratio", "color" : ROOT.kBlack,  "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 1 / self.upperSplit}
                        rnewInfo = copy.deepcopy(newInfo)
                        rnewInfo["Y"]["title"] = "#frac{Data}{Pred.}"
                        rnewInfo["X"]["rebin"] = 1

                        canvas.cd(2)
                        ROOT.gPad.SetGridy()
                        ratio = Histogram(ratio, None, self.upperSplit/self.scale, self.lowerSplit/self.scale, None, rnewInfo, rinfo, 0.8).histogram

                        ratio.SetMinimum(0.4); ratio.SetMaximum(1.6)
                        ratio.GetYaxis().SetNdivisions(5, 5, 0)

                        ratio.Draw("E0P")

                    canvas.SaveAs("%s/%s_%s.pdf"%(self.outpath, self.year, newName))
                    canvas.SaveAs("%s/%s_%s.png"%(self.outpath, self.year, newName))

if __name__ == "__main__":

    usage = "usage: %stackPlotter [options]"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("--noRatio",      dest="noRatio",      help="No ratio plot",             default=False,  action="store_true") 
    parser.add_argument("--approved",     dest="approved",     help="Plot is approved",          default=False,  action="store_true") 
    parser.add_argument("--printNEvents", dest="printNEvents", help="Show number of events",     default=False,  action="store_true") 
    parser.add_argument("--normMC2Data",  dest="normMC2Data",  help="Normalize MC to data",      default=False,  action="store_true") 
    parser.add_argument("--normalize",    dest="normalize",    help="Normalize all to unity",    default=False,  action="store_true") 
    parser.add_argument("--printSign",    dest="printSign",    help="Print simple significance", default=False,  action="store_true") 
    parser.add_argument("--inpath",       dest="inpath",       help="Path to root files",        default="NULL", required=True)
    parser.add_argument("--outpath",      dest="outpath",      help="Where to put plots",        default="NULL", required=True)
    parser.add_argument("--year",         dest="year",         help="which year",                                required=True)
    parser.add_argument("--options",      dest="options",      help="options file",              default="stackPlotter_aux", type=str)
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

    plotter = StackPlotter(args.approved, args.noRatio, args.printNEvents, args.printSign, args.year, args.outpath, args.inpath, args.normMC2Data, args.normalize, histograms, selections, backgrounds, signals, data)
    plotter.makePlots()
