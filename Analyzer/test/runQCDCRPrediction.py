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
                self.histogram.SetFillColor(self.info["color"])

            if "fstyle" in self.info:
                self.histogram.SetFillStyle(self.info["fstyle"])

            self.histogram.RebinX(self.info["X"]["rebin"])
            self.histogram.GetXaxis().SetRangeUser(self.info["X"]["min"], self.info["X"]["max"])

            self.histogram.SetTitle("")
        else:
            self.histogram = -1

class CRProducer:

    def __init__(self, approved, noRatio, printNEvents, printInfo, year, outpath, inpath, normMC2Data, normalize, controlRegions, signalRegions, backgrounds, signals, data):

        self.controlRegions  = controlRegions
        self.signalRegions   = signalRegions
        self.backgrounds     = backgrounds
        self.signals         = signals
        self.data            = data

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

        self.mainBG = "QCD"
        self.dataQCDOnly = None
        self.transforFactorsHisto = {}

    # Main function to compose the full stack plot with or without a ratio panel
    def getCRData(self):
        # Top loop begins going over each histo-to-be-plotted
        for hname, hinfo in self.controlRegions.items():
            newName = hname
            newInfo = copy.deepcopy(hinfo)
            newInfo["X"]["title"] = hinfo["X"]["title"]
            totalMC = None
            firstBG = False

            # Loop over each background and get their respective histo, scale if necessary
            for bname, binfo in self.backgrounds.items():
                if bname == self.mainBG: 
                    continue
                rootFile = "%s/%s_%s.root"%(self.inpath, self.year, bname)
                Hobj = Histogram(None, rootFile, 1.0, 1.0, newName, newInfo, binfo)
                if not firstBG:
                    totalMC = Hobj.Clone("totalMC%s"%(newName))
                    firstBG = True
                else:
                    totalMC.Add(Hobj.histogram)
            
            # Loop over the data and get their respective histo
            dataHist = None
            for dname, dinfo in self.data.items():
                rootFile = "%s/%s_%s.root"%(self.inpath, self.year, dname)
                dataHist = Histogram(None, rootFile, 1.0, 1.0, newName, newInfo, dinfo)

            self.dataQCDOnly = dataHist.Clone("{}_Data_only_{}_{}".format(self.year, self.mainBG, hname.replace("h_njets_12incl_","").replace("_ABCD","")))
            self.dataQCDOnly.Add(totalMC, -1)
            print("Data QCD only", self.dataQCDOnly.Integral())

    def getTransferFactors(self):
        rootFile = "%s/%s_%s.root"%(self.inpath, self.year, self.mainBG)

        den = {}
        for hname, hinfo in self.controlRegions.items():
            newName = hname
            newInfo = copy.deepcopy(hinfo)
            newInfo["X"]["title"] = hinfo["X"]["title"]
            bname, binfo = self.mainBG, self.backgrounds[self.mainBG]
            Hobj = Histogram(None, rootFile, 1.0, 1.0, newName, newInfo, binfo)
            nEvents = Hobj.Integral()
            error = Hobj.IntegralError()
            print(hname, nEvents, "+/-", error)
            den[hname.replace("h_njets_12incl_", "").replace("_ABCD","")] = (nEvents, error, Hobj)

        num = {}
        for hname, hinfo in self.signalRegions.items():
            newName = hname
            newInfo = copy.deepcopy(hinfo)
            newInfo["X"]["title"] = hinfo["X"]["title"]
            bname, binfo = self.mainBG, self.backgrounds[self.mainBG]
            Hobj = Histogram(None, rootFile, 1.0, 1.0, newName, newInfo, binfo)
            nEvents = Hobj.Integral()
            error = Hobj.IntegralError()
            print(hname, nEvents, "+/-", error)
            num[hname.replace("h_njets_12incl_", "").replace("_ABCD","")] = (nEvents, error, Hobj)

        transferFactors = {}
        for nameNum, n in num.items():
            for nameDen, d in den.items():
                r = n[0]/d[0]
                transferFactors["{}_TF_{}Over{}".format(self.year,nameNum,nameDen)] = (r, r*math.sqrt((n[1]/n[0])**2 + (d[1]/d[0])**2), n[2], d[2])

        for name, tf in transferFactors.items():
            h = ROOT.TH1D(name, name, 1, 0, 1)
            h.SetBinContent(1, tf[0])
            h.SetBinError(  1, tf[1])
            self.transforFactorsHisto[name] = (h, tf[2])

    def write(self):
        outfile = ROOT.TFile.Open("{}/QCD_Prediction.root".format(self.outpath),"RECREATE")
        self.dataQCDOnly.Write()

        for name, h in self.transforFactorsHisto.items():
            h[0].Write()
            h[1].Write()

        outfile.Close()

    def makeOutputPlots(self):
        canvas = ROOT.TCanvas("", "", 900, 900)
        ROOT.gPad.SetLogy()

        for name, h in self.transforFactorsHisto.items():
            hCRPred = self.dataQCDOnly.Clone("")
            hCRPred.Scale(h[0].GetBinContent(1))
            hCRPred.Draw()
            h[1].histogram.Draw("same")
            canvas.SaveAs("%s/%s_%s.png"%(self.outpath, self.year, name))


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
    parser.add_argument("--inpath",       dest="inpath",       help="Path to root files",          default="condor/Run2UL_DisCo_outputs_0L_1L_04.10.2022/")
    parser.add_argument("--outpath",      dest="outpath",      help="Where to put plots",          default="output_QCDCRPrediction")
    parser.add_argument("--year",         dest="year",         help="which year",                  default="Run2UL"               )
    parser.add_argument("--options",      dest="options",      help="options file",                default="runQCDCRPrediction_aux", type=str)
    args = parser.parse_args()

    # The auxiliary file contains many "hardcoded" items
    # describing which histograms to get and how to draw
    # them. These things are changed often by the user
    # and thus are kept in separate sidecar file.
    importedGoods = __import__(args.options)

    # Names of histograms, rebinning, titles, ranges, etc.
    controlRegions = importedGoods.controlRegions
    signalRegions  = importedGoods.signalRegions

    # Background/signal/data categories to plot, and plotting options
    backgrounds = importedGoods.backgrounds
    signals     = importedGoods.signals
    data        = importedGoods.data

    crProducer = CRProducer(args.approved, args.noRatio, args.printNEvents, args.printInfo, args.year, args.outpath, args.inpath, args.normMC2Data, args.normalize, controlRegions, signalRegions, backgrounds, signals, data)
    print("-----Make small MC subtracted Data histograms from the CR-----")
    crProducer.getCRData()
    print("-----Calcualte all transfer factors-----")
    crProducer.getTransferFactors()
    print("-----Save all output histograms-----")
    crProducer.write()
    crProducer.makeOutputPlots()

    if args.makeTable:
        plotter.make_njetsTable()
        plotter.makeTable()
