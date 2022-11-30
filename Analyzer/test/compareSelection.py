#! /usr/bin/env python

import ROOT, random, os, argparse, string, copy, math, ctypes
from stackPlotter import Histogram
ROOT.PyConfig.IgnoreCommandLineOptions = True

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat("")
ROOT.gStyle.SetPaintTextFormat("3.2f")
ROOT.gStyle.SetFrameLineWidth(2)
ROOT.gStyle.SetEndErrorSize(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH2.SetDefaultSumw2()

class comparePlotter:

    def __init__(self, year, outpath, inpath, histograms, selections, samples):

        self.histograms  = histograms
        self.selections  = selections
        self.samples     = samples

        self.year        = year
        self.inpath      = inpath
        self.outpath     = outpath
        self.noRatio     = True
        
        if not os.path.exists(self.outpath):
            os.makedirs(self.outpath)

        # Customized numbers that are scaled
        # to work with or without a ratio plot
        self.TopMargin    = 0.06
        self.BottomMargin = 0.12
        self.RightMargin  = 0.04
        self.LeftMargin   = 0.16

    # Create a canvas and determine if it should be split for a ratio plot
    # Margins are scaled on-the-fly so that distances are the same in either
    # scenario.
    def makeCanvas(self, doLogY):

        randStr = ''.join([random.choice(string.ascii_letters + string.digits) for n in range(8)])
       
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

        textSize = 0.030 / self.upperSplit
        nColumns = 1
        maxBkgLegendFrac = 0.25 * (1.0 - self.TopMargin - self.BottomMargin)

        tooManyBkgds = nBkgs * textSize * 1.2 > maxBkgLegendFrac
        if tooManyBkgds:
            nColumns = 1

        space    = 0.015

        bkgXmin = 0.70 - self.RightMargin
        if tooManyBkgds:
            bkgXmin = self.LeftMargin + 0.05
            
        bkgYmax = 1.0-(self.TopMargin/self.upperSplit)-0.01
        bkgXmax = 1.0-self.RightMargin-0.01
        bkgYmin = bkgYmax-(float(nBkgs)/float(nColumns+5))*(1.2*textSize+space)
        bkgYFrac = (1.0-self.TopMargin-bkgYmin) / (1.0 - self.TopMargin - self.BottomMargin)

        bkgLegend = ROOT.TLegend(bkgXmin, bkgYmin, bkgXmax, bkgYmax)
        bkgLegend.SetBorderSize(0)
        bkgLegend.SetTextSize(textSize)
        if tooManyBkgds:
            bkgLegend.SetNColumns(nColumns)

        sigXmin = self.LeftMargin+0.03
        sigYmax = bkgYmax
        sigXmax = bkgXmin
        sigYmin = bkgYmax-nSigs*(textSize+space) 
        sigYFrac = (1.0-self.TopMargin-sigYmin) / (1.0 - self.TopMargin - self.BottomMargin)

        sigLegend = ROOT.TLegend(sigXmin, sigYmin, sigXmax, sigYmax)
        sigLegend.SetBorderSize(0)
        sigLegend.SetMargin(0.10)
        sigLegend.SetTextSize(textSize)

        yMax = 1.0; factor = 1.0; power = 1.0
        if doLogY and theMax != 0.0 and theMin != 0.0:
            power = math.log10(theMax / theMin) * 3.0

        theFrac = bkgYFrac if bkgYFrac > sigYFrac else sigYFrac                           
        #yMax = (theMax-theMin) * (0.95 - theFrac)**(-power) * factor

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

        # if self.approved:
        #     mark.DrawLatex(self.LeftMargin + 0.12, 1 - (self.TopMargin - 0.017), "")
        # elif (self.wip):
        #     mark.DrawLatex(self.LeftMargin + 0.12, 1 - (self.TopMargin - 0.017), "Work in Progress")
        # else:
        mark.DrawLatex(self.LeftMargin + 0.12, 1 - (self.TopMargin - 0.017), "Preliminary")

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
            orders = [""]
            if "orders" in hinfo:
                orders = hinfo["orders"]

            for order in orders:
                order = str(order)
                newInfo = copy.deepcopy(hinfo)
                newInfo["X"]["title"] = hinfo["X"]["title"].replace("@", order)

                canvas = self.makeCanvas(newInfo["logY"])
                if self.noRatio:
                    canvas.cd()
                else:
                    canvas.cd(1)

                firstDraw = False
                dummy = None

                nSigLegend = 0
                theMin = 2e-6
                theMax = 0.0

                for sname, sinfo in self.samples.items():
                    for selection, color in self.selections:
                        newName = hname.replace("@", order).replace("?", "%s"%(selection))
                        #newName = hname.replace("?", "%s"%(selection))

                        rootFile = "%s/%s_%s.root"%(self.inpath, self.year, sname)
                        Hobj = Histogram(None, rootFile, self.upperSplit, self.lowerSplit, newName, newInfo, sinfo)
                        dummy = Hobj.Clone("dummy%s"%(hname)); dummy.Reset("ICESM")
                        if Hobj.IsGood():
                            scale = Hobj.Integral()
                            if scale != 0.0: 
                                Hobj.Scale(1.0 / scale)
                                tempMax = Hobj.histogram.GetMaximum()/scale
                            if tempMax > theMax:
                                theMax = tempMax

                # Here we get the bkgd and sig legends as well as a tuned maximum for the canvas to avoid overlap
                sigLegend, _, yMax = self.makeLegends(len(self.samples)*len(selection), len(self.samples)*len(selection), newInfo["logY"], theMin, theMax)

                if newInfo["logY"]:
                    yMax=5000
                else:
                    yMax=1
                dummy.SetMaximum(yMax)
                dummy.SetMinimum(theMin)
                dummy.Draw()

                # Loop over each signal and get their respective histo
                option = "HIST"; loption = "L"
                for sname, sinfo in self.samples.items():
                    for selection, color in self.selections:
                        newName = hname.replace("@", order).replace("?", "%s"%(selection))
                        #newName = hname.replace("?", "%s"%(selection))
                        newInfo["color"] = sinfo["color"] + color
                        
                        rootFile = "%s/%s_%s.root"%(self.inpath, self.year, sname)
                        if "option"  not in sinfo: sinfo["option"]  = option 
                        if "loption" not in sinfo: sinfo["loption"] = loption
                        Hobj = Histogram(None, rootFile, self.upperSplit, self.lowerSplit, newName, newInfo, sinfo) 
                        scale = Hobj.Integral()
                        if scale != 0.0: 
                            Hobj.Scale(1.0 / scale)
                        Hobj.histogram.Draw("same")
                        #sigLegend.AddEntry(Hobj.histogram, selection[1:].replace("_ABCD",""), "l")
                        sigLegend.AddEntry(Hobj.histogram, selection[1:].replace("_ABCD","")+"_{}".format(sname), "l")
                        #nSigLegend, firstDraw = Hobj.Draw(canvas, False, firstDraw, nSigLegend, sigLegend)

                if Hobj.IsGood():
                    ratio = Hobj.Clone("ratio")
                    #ratio.Divide(totalMC)

                sigLegend.Draw("SAME")

                self.addCMSlogo(canvas)

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

                canvas.SaveAs("%s/%s_%s.png"%(self.outpath, self.year, hname.replace("?","").replace("@", order)))

if __name__ == "__main__":

    usage = "usage: %stackPlotter [options]"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("--year",         dest="year",         help="which year",                  default="Run2UL"               )
    parser.add_argument("--outpath",      dest="outpath",      help="Where to put plots",          default="output_QCDCRPrediction")
    parser.add_argument("--inpath",       dest="inpath",       help="Path to root files",          default="condor/Hadded_2016preVFP_DoubleDisco_22-10-24")
    parser.add_argument("--options",      dest="options",      help="options file",                default="compareSelection_aux", type=str)
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
    samples = importedGoods.samples

    plotter = comparePlotter(args.year, args.outpath, args.inpath, histograms, selections, samples)
    plotter.makePlots()
