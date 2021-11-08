#! /usr/bin/env python

import ROOT, random, os, argparse, string, copy, math
ROOT.PyConfig.IgnoreCommandLineOptions = True

from stackPlotter import Histogram

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat("")
ROOT.gStyle.SetPaintTextFormat("3.2f")
ROOT.gStyle.SetFrameLineWidth(2)
ROOT.gStyle.SetEndErrorSize(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH2.SetDefaultSumw2()

class Superimpose:

    def __init__(self, year, inpath, outpath, backgrounds, signals, data, histograms, qcdCR1_selections, base_selections, printNEvents, printInfo, normalize, normMC2Data):

        self.year              = year
        self.inpath            = inpath
        self.outpath           = outpath
        self.backgrounds       = backgrounds
        self.signals           = signals
        self.data              = data
        self.histograms        = histograms
        self.qcdCR1_selections = qcdCR1_selections
        self.base_selections   = base_selections
        self.printNEvents      = printNEvents
        self.printInfo         = printInfo
        self.normalize         = normalize
        self.normMC2Data       = normMC2Data

        if not os.path.exists(self.outpath):
            os.makedirs(self.outpath)

        # Customized numbers that are scaled
        # to work with or without a ratio plot
        self.TopMargin    = 0.06
        self.BottomMargin = 0.12
        self.RightMargin  = 0.04
        self.LeftMargin   = 0.16

    # -----------
    # make canvas
    # -----------
    def makeCanvas(self, doLogY):

        randStr = ''.join([random.choice(string.ascii_letters + string.digits) for n in xrange(8)])
       
        canvas = ROOT.TCanvas(randStr, randStr, 900, 900)

        ROOT.gPad.SetTopMargin(self.TopMargin)
        ROOT.gPad.SetBottomMargin(self.BottomMargin)
        ROOT.gPad.SetLeftMargin(self.LeftMargin)
        ROOT.gPad.SetRightMargin(self.RightMargin)
        if doLogY:
            ROOT.gPad.SetLogy()

        return canvas
    
    # ---------------------    
    # add cms logo to plots
    # --------------------- 
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

        mark.DrawLatex(self.LeftMargin + 0.12, 1 - (self.TopMargin - 0.017), "%sPreliminary"%(simStr))
        mark.SetTextFont(42)
        mark.SetTextAlign(31)
        if "Run 2" in self.year:
            mark.DrawLatex(1 - self.RightMargin, 1 - (self.TopMargin - 0.017), "Run 2 (13 TeV)")
        else:
            mark.DrawLatex(1 - self.RightMargin, 1 - (self.TopMargin - 0.017), "%s (13 TeV)"%(self.year))


    # ---------------
    # make the legend
    # ---------------
    def makeLegends(self, nBkgs):

        textSize = 0.035 
        space    = 0.015

        bkgXmin = 0.5
        bkgYmax = 1.0-(self.TopMargin)-0.01
        bkgXmax = 1.0-self.RightMargin-0.01
        bkgYmin = bkgYmax-nBkgs*(textSize+space)
        bkgYFrac = (1.0-self.TopMargin-bkgYmin) / (1.0 - self.TopMargin - self.BottomMargin)

        legend = ROOT.TLegend(bkgXmin, bkgYmin, bkgXmax, bkgYmax)
        legend.SetBorderSize(0)
        #legend.SetMargin(0.10)
        legend.SetTextSize(textSize)

        return legend 


    # -----------------------------------------
    # superimpose qcdCR and baseline histograms
    # ----------------------------------------- 
    def superimpose_baseline_qcdCR(self):

        for hname, hinfo in self.histograms.items():

            orders            = [-1]
            qcdCR1_selections = []; hist_qcdCR = []
            base_selections   = []; hist_base  = []

            if "orders" in hinfo:
                orders = hinfo["orders"]

            if "qcdCR1_selections" in hinfo:
                selections = hinfo["qcdCR1_selections"]

            if "base_selections" in hinfo:
                selections = hinfo["base_selections"]

            for order in orders:

                # put the njet by njet qcd histograms to a list
                for qcdCR1 in range(0, len(self.qcdCR1_selections)):
                
                    hist_qcdCR.append( hname.replace("@", "%d"%(order)).replace("?", "%s"%(self.qcdCR1_selections[qcdCR1])) )

                # put the njet by njet baseline histograms to a list 
                for base in range(0, len(self.base_selections)):
         
                    hist_base.append( hname.replace("@", "%d"%(order)).replace("?", "%s"%(self.base_selections[base])) )
                    
                    # get the histogram information to set up axes labels 
                    newInfo = copy.deepcopy(hinfo)
                    newInfo["X"]["title"] = hinfo["X"]["title"].replace("@", "%d"%(order))

                # superimpose qcdCR and baseline histograms njet by njet
                for i in range(0, len(hist_qcdCR)):

                    hist1 = hist_qcdCR[i] 
                    hist2 = hist_base[i]

                    # get canvas and legend
                    canvas = self.makeCanvas(newInfo["logY"])
                    legend = self.makeLegends(2)

                    # open QCD root file and get histograms
                    filename = self.inpath + "/" + self.year + "_QCD.root"
                    f = ROOT.TFile.Open(filename, "READ")
                    qcdHist  = f.Get("%s"%(hist1))
                    baseHist = f.Get("%s"%(hist2))

                    # normalize the histogram
                    scale_qcdCR = qcdHist.Integral()
                    scale_base  = baseHist.Integral()

                    if self.normalize and scale_qcdCR != 0.0 and scale_base != 0.0:
                        qcdHist.Scale(1.0 / scale_qcdCR)
                        baseHist.Scale(1.0 / scale_base)

                    hist1_label = "_".join(hist1.split("_")[3:])
                    hist2_label = "_".join(hist2.split("_")[3:])

                    # get qcd & baseline histograms
                    qcdHist.SetTitle("")
                    qcdHist.Rebin(newInfo["X"]["rebin"])
                    qcdHist.GetXaxis().SetRangeUser(newInfo["X"]["min"], newInfo["X"]["max"])
                    qcdHist.GetXaxis().SetTitle(newInfo["X"]["title"])
                    qcdHist.GetYaxis().SetTitle("A.U.")
                    qcdHist.SetMaximum(0.75)
                    qcdHist.SetLineColor(30)                        
                    qcdHist.SetLineWidth(4)
                    baseHist.Rebin(newInfo["X"]["rebin"])
                    baseHist.SetLineColor(38)                       
                    baseHist.SetLineWidth(4) 
                    qcdHist.Draw("E HIST")                       
                    baseHist.Draw("SAME E HIST")
                    #legend.AddEntry(qcdHist,  "%s, %.2f"%(hist1_label,scale_qcdCR), "l")
                    #legend.AddEntry(baseHist, "%s, %.2f"%(hist2_label, scale_base), "l")
                    legend.AddEntry(qcdHist,  "0l QCD CR, %.2f"%(scale_qcdCR),  "l")
                    legend.AddEntry(baseHist, "0l Baseline, %.2f"%(scale_base), "l")
                    legend.Draw("SAME")
                    self.addCMSlogo(canvas)
                    canvas.SaveAs("%s/%s_%s.pdf"%(self.outpath, self.year, hist1))


if __name__ == "__main__":

    usage = "usage: %superimpose_baseline_qcdCR [options]"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("--noRatio",      dest="noRatio",      help="No ratio plot",               default=False,  action="store_true") 
    parser.add_argument("--approved",     dest="approved",     help="Plot is approved",            default=False,  action="store_true") 
    parser.add_argument("--printNEvents", dest="printNEvents", help="Show number of events",       default=False,  action="store_true") 
    parser.add_argument("--normMC2Data",  dest="normMC2Data",  help="Normalize MC to data",        default=False,  action="store_true") 
    parser.add_argument("--normalize",    dest="normalize",    help="Normalize all to unity",      default=False,  action="store_true") 
    parser.add_argument("--printInfo",    dest="printInfo",    help="Print significance and cuts", default=False,  action="store_true")
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
    qcdCR1_selections = importedGoods.qcdCR1_selections
    base_selections   = importedGoods.base_selections 
    backgrounds       = importedGoods.backgrounds
    signals           = importedGoods.signals
    data              = importedGoods.data

    superimposer = Superimpose(args.year, args.inpath, args.outpath, backgrounds, signals, data, histograms, qcdCR1_selections, base_selections, args.printNEvents, args.printInfo, args.normalize, args.normMC2Data)
    superimposer.superimpose_baseline_qcdCR()




