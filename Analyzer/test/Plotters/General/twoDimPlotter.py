#! /usr/bin/env python

import ROOT, random, os, argparse, string, copy, math
ROOT.PyConfig.IgnoreCommandLineOptions = True

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat("")
ROOT.gStyle.SetPaintTextFormat("3.3f")
ROOT.gStyle.SetFrameLineWidth(2)
ROOT.gStyle.SetEndErrorSize(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH2.SetDefaultSumw2()

# This class owns a histogram and can take a ROOT TH1 or
# a path to extract the histogram from
#    histo      : ROOT TH1 object to hold onto
#    rootFile   : input path to ROOT file containing histos
#    histoName  : name of histo to be extracted from file
#    hinfo      : dictionary of params for customizing histo
#    pinfo      : additional params for customizing histo
#    aux        : auxiliary multiplier for offsetting y-axis label
class Histogram:

    def __init__(self, histo=None, rootFile=None, histoName=None, hinfo=None, pinfo=None, aux = 1.0):

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
        self.xDim["offset"] = 1.1;   self.yDim["offset"] = 1.1 

        self.histogram = histo

        self.setupHisto(aux)

    def Integral(self):
        return self.histogram.Integral()
        
    def Clone(self, name):
        return self.histogram.Clone(name)

    def Scale(self, fraction):
        return self.histogram.Scale(fraction)

    def IsGood(self):
        return self.histogram != -1

    def Draw(self, canvas):

        if self.histogram != -1:
            lname = self.info["name"]

            self.histogram.Draw("%s SAME"%(self.info["option"]))
            
    def Write(self, name):
        n = "{}_{}".format(name,self.histoName)
        self.histogram.SetName(n)
        self.histogram.SetTitle(n)
        self.histogram.Write()

    def Divide(self, hDen):
        self.histogram.Divide(hDen.histogram)
    
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

    def setupHisto(self, aux):

        code = self.getHisto()

        if code != -1:
            self.histogram.GetXaxis().SetTitleSize(self.xDim["title"]);    self.histogram.GetYaxis().SetTitleSize(self.yDim["title"])
            self.histogram.GetXaxis().SetLabelSize(self.xDim["label"]);    self.histogram.GetYaxis().SetLabelSize(self.yDim["label"])
            self.histogram.GetXaxis().SetTitleOffset(self.xDim["offset"]); self.histogram.GetYaxis().SetTitleOffset(self.yDim["offset"] * aux)
            self.histogram.GetXaxis().SetTitle(self.info["X"]["title"]);   self.histogram.GetYaxis().SetTitle(self.info["Y"]["title"])

            if "rebin" in self.info["X"]:
                self.histogram.RebinX(self.info["X"]["rebin"])
            if "rebin" in self.info["Y"]:
                self.histogram.RebinY(self.info["Y"]["rebin"])

            self.histogram.GetXaxis().SetRangeUser(self.info["X"]["min"], self.info["X"]["max"])
            self.histogram.GetYaxis().SetRangeUser(self.info["Y"]["min"], self.info["Y"]["max"])

            self.histogram.SetMarkerSize(1.7)

            self.histogram.SetTitle("")

            #self.histogram.GetXaxis().SetNdivisions(-4, "I")
            #self.histogram.GetYaxis().SetNdivisions(-4, "I")

        else:
            self.histogram = -1

# The StackPlotter class oversees the creation of all stack plots
#     approved     : are these plots approved, then no "Preliminary"
#     wip          : a work in progress, put to plot heading
#     year         : corresponding year for the plots/inputs
#     outpath      : where to put the plots, path is created if missing
#     inpath       : where the input ROOT files are located
#     normalize    : normalize all processes to unity area
#     histograms   : dictionary containing config info for desired histos
#     selections   : list of cut strings appearing in the names of plots
#     backgrounds  : dictionary containing config info for desired backgrounds
#     signals      : dictionary containing config info for desired signals
#     data         : dictionary containing config info for data
class TwoDimPlotter:

    def __init__(self, approved, wip, year, outpath, inpath, normalize, histograms, selections, backgrounds, signals, data, outrootname):

        self.histograms  = histograms
        self.selections  = selections
        self.backgrounds = backgrounds
        self.signals     = signals
        self.data        = data

        self.year        = year
        self.inpath      = inpath
        self.outpath     = outpath
        self.outrootname = outrootname
        self.outrootfile = ROOT.TFile("{}/{}.root".format(self.outpath, self.outrootname),"RECREATE") if self.outrootname else None
        
        if not os.path.exists(self.outpath):
            os.makedirs(self.outpath)

        self.approved     = approved
        self.wip          = wip
        self.normalize    = normalize

        # Customized numbers that are scaled
        # to work with or without a ratio plot
        self.TopMargin    = 0.06
        self.BottomMargin = 0.12
        self.RightMargin  = 0.12
        self.LeftMargin   = 0.12

    def __del__(self):
        if self.outrootfile:
            self.outrootfile.Close()

    def WriteHisto(self, h, name):
        if self.outrootfile:
            self.outrootfile.cd()
            h.Write(name)

    # Create a canvas and determine if it should be split for a ratio plot
    # Margins are scaled on-the-fly so that distances are the same in either
    # scenario.
    def makeCanvas(self, doLogX, doLogY, doLogZ):

        randStr = ''.join([random.choice(string.ascii_letters + string.digits) for n in xrange(8)])
       
        canvas = ROOT.TCanvas(randStr, randStr, 900, 900)

        ROOT.gPad.SetTopMargin(self.TopMargin)
        ROOT.gPad.SetBottomMargin(self.BottomMargin)
        ROOT.gPad.SetLeftMargin(self.LeftMargin)
        ROOT.gPad.SetRightMargin(self.RightMargin)
        if doLogX:
            ROOT.gPad.SetLogx()
        if doLogY:
            ROOT.gPad.SetLogy()
        if doLogZ:
            ROOT.gPad.SetLogz()

        return canvas

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

        if self.approved:
            mark.DrawLatex(self.LeftMargin + 0.12, 1 - (self.TopMargin - 0.017), "Supplementary")
        elif self.wip:
            mark.DrawLatex(self.LeftMargin + 0.12, 1 - (self.TopMargin - 0.017), "Work in Progress")
        else:
            mark.DrawLatex(self.LeftMargin + 0.12, 1 - (self.TopMargin - 0.017), "Preliminary")

        mark.SetTextFont(42)
        mark.SetTextAlign(31)
        if "Run 2" in self.year:
            mark.DrawLatex(1 - self.RightMargin, 1 - (self.TopMargin - 0.017), "Run 2 (13 TeV)")
        else:
            mark.DrawLatex(1 - self.RightMargin, 1 - (self.TopMargin - 0.017), "%s (13 TeV)"%(self.year))

    def addExtraInfo(self, canvas, packedInfo):

        njetStr = ""
        if "Njets" in packedInfo:
            njets = packedInfo.partition("Njets")[-1].partition("_")[0]
    
            op = "="
            if "incl" in njets:
                op = "#geq"

            njetStr = "N_{ jets} %s %s"%(op, njets.replace("incl", ""))

        modelStr = ""
        if "RPV" in packedInfo:
            modelStr = "RPV"
        elif "SYY" in packedInfo:
            modelStr = "Stealth SYY"

        channelStr = ""
        if "0l" in packedInfo:
            channelStr = "Fully-Hadronic"
        elif "1l" in packedInfo:
            channelStr = "Semi-Leptonic"
        elif "2l" in packedInfo:
            channelStr = "Fully-Leptonic"

        canvas.cd()

        text = ROOT.TLatex()
        text.SetNDC(True)

        text.SetTextAlign(13)
        text.SetTextSize(0.04)
        text.SetTextFont(42)
        text.SetTextColor(ROOT.TColor.GetColor("#7C99D1"))

        offset = 0.0
        if channelStr != "":
            text.DrawLatex(self.LeftMargin + 0.04, 1 - (self.TopMargin + 0.02 + offset), channelStr)
            offset += 0.05
        if modelStr != "":
            text.DrawLatex(self.LeftMargin + 0.04, 1 - (self.TopMargin + 0.02 + offset), modelStr)
            offset += 0.05
        if njetStr != "":
            text.DrawLatex(self.LeftMargin + 0.04, 1 - (self.TopMargin + 0.02 + offset), njetStr)

    def preprocess(self, container, hname, hinfo, theMax):

        # Preemptively get data counts to be used for normalzing the histograms later
        for name, info in container.items():

            rootFile = "%s/%s_%s.root"%(self.inpath, self.year, name)

            Hobj = Histogram(None, rootFile, hname, hinfo, info)
            if Hobj.IsGood():
                if self.normalize:
                    Hobj.Scale(1.0 / Hobj.Integral())

                tempMax = Hobj.histogram.GetMaximum()
                if tempMax > theMax:
                    theMax = tempMax

        return theMax

    def plotLoop(self, container, hname, hinfo, theMax):
        # Loop over each background and get their respective histo, scale if necessary
        for name, info in container.items(): 

            rootFile = "%s/%s_%s.root"%(self.inpath, self.year, name)
    
            Hobj = Histogram(None, rootFile, hname, hinfo, info)

            if Hobj.IsGood(): 
                if self.normalize:
                    Hobj.Scale(1.0 / Hobj.Integral())

            canvas = self.makeCanvas(hinfo["logX"], hinfo["logY"], hinfo["logZ"])
            canvas.cd()

            Hobj.histogram.SetContour(255)
            Hobj.Draw(canvas)
            self.addCMSlogo(canvas)
            self.addExtraInfo(canvas, hname)
            saveName = "%s/%s_%s_%s"%(self.outpath, self.year, name, hname)
            canvas.SaveAs(saveName+".pdf")
            canvas.SaveAs(saveName+".gif")

            self.WriteHisto(Hobj, name)

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
         
                    newName = hname.replace("@", str(order)).replace("?", selection)
                    newInfo = copy.deepcopy(hinfo)
                    newInfo["X"]["title"] = hinfo["X"]["title"].replace("@", str(order))

                    theMax = 0.0
                    theMax = self.preprocess(self.data,        newName, newInfo, theMax)
                    theMax = self.preprocess(self.backgrounds, newName, newInfo, theMax)
                    theMax = self.preprocess(self.signals,     newName, newInfo, theMax)

                    self.plotLoop(self.data,        newName, newInfo, theMax)
                    self.plotLoop(self.backgrounds, newName, newInfo, theMax)
                    self.plotLoop(self.signals,     newName, newInfo, theMax)

if __name__ == "__main__":

    usage = "usage: %twoDimPlotter [options]"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("--approved",     dest="approved",     help="Plot is approved",         default=False,  action="store_true") 
    parser.add_argument("--wip",          dest="wip",          help="Plot is work in progress", default=False,  action="store_true") 
    parser.add_argument("--normalize",    dest="normalize",    help="Normalize all to unity",   default=False,  action="store_true") 
    parser.add_argument("--inpath",       dest="inpath",       help="Path to root files",       default="NULL", required=True)
    parser.add_argument("--outpath",      dest="outpath",      help="Where to put plots",       default="NULL", required=True)
    parser.add_argument("--year",         dest="year",         help="which year",                               required=True)
    parser.add_argument("--options",      dest="options",      help="options file",             default="twoDimPlotter_aux", type=str)
    parser.add_argument("--outrootname",  dest="outrootname",  help="outrootname",              default=None)
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

    plotter = TwoDimPlotter(args.approved, args.wip, args.year, args.outpath, args.inpath, args.normalize, histograms, selections, backgrounds, signals, data, args.outrootname)
    plotter.makePlots()
