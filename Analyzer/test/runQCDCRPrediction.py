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

class ControlRegionProducer:

    def __init__(self, year, outpath, inpath, controlRegions, signalRegions, backgrounds, data):

        self.year        = year
        self.outpath     = outpath
        self.inpath      = inpath        
        self.controlRegions  = controlRegions
        self.signalRegions   = signalRegions
        self.backgrounds     = backgrounds
        self.data            = data
        self.approved        = None
        
        if not os.path.exists(self.outpath):
            os.makedirs(self.outpath)

        self.mainBG = "QCD"
        self.dataQCDOnly = None
        self.transforFactorsHisto = {}

        self.TopMargin    = 0.06
        self.BottomMargin = 0.12
        self.RightMargin  = 0.04
        self.LeftMargin   = 0.16

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

    def getCRData(self):
        print("-----Make small MC subtracted Data histograms from the CR-----")
        for hname, hinfo in self.controlRegions.items():
            newName = hname
            newInfo = copy.deepcopy(hinfo)
            newInfo["X"]["title"] = hinfo["X"]["title"]
            totalMC = None
            firstBG = False
            mainBGHisto = None

            for bname, binfo in self.backgrounds.items():
                rootFile = "%s/%s_%s.root"%(self.inpath, self.year, bname)
                Hobj = Histogram(None, rootFile, 1.0, 1.0, newName, newInfo, binfo)
                if bname == self.mainBG: 
                    mainBGHisto = Hobj
                    continue
                if not firstBG:
                    totalMC = Hobj.Clone("totalMC%s"%(newName))
                    firstBG = True
                else:
                    totalMC.Add(Hobj.histogram)
            
            dataHist = None
            for dname, dinfo in self.data.items():
                rootFile = "%s/%s_%s.root"%(self.inpath, self.year, dname)
                dataHist = Histogram(None, rootFile, 1.0, 1.0, newName, newInfo, dinfo)

            self.dataQCDOnly = dataHist.Clone("{}_Data_only_{}_{}".format(self.year, self.mainBG, hname.replace("h_njets_12incl_","").replace("_ABCD","")))
            self.dataQCDOnly.Add(totalMC, -1)

            #Correct QCD normalization per ABCD bin
            #print("Data QCD only", self.dataQCDOnly.Integral())
            #print("MC   QCD only", mainBGHisto.Integral())
            #CRCorrection = self.dataQCDOnly.Integral()/mainBGHisto.Integral()
            #print("CR correction", CRCorrection)
            #self.dataQCDOnly.Scale(1/CRCorrection)
            #self.dataQCDOnly.Print("all")

    def getABCDCounts(self, h):
        cA, eA = h.Integral( 0, 6), ctypes.c_double(-999)
        cB, eB = h.Integral( 7,12), ctypes.c_double(-999)
        cC, eC = h.Integral(13,18), ctypes.c_double(-999)
        cD, eD = h.Integral(19,25), ctypes.c_double(-999)
        h.IntegralAndError( 0, 6, eA)
        h.IntegralAndError( 7,12, eB)
        h.IntegralAndError(13,18, eC)
        h.IntegralAndError(19,25, eD)

        return ([cA,cB,cC,cD],[eA.value,eB.value,eC.value,eD.value])

    def calculateTF(self, n, d, nE, dE):
        r = n/d
        e = r*math.sqrt( (nE/n)**2 + (dE/d)**2 )
        return (r, e)

    def getTransferFactors(self):
        print("-----Calcualte all transfer factors-----")
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
            cABCD, eABCD = self.getABCDCounts(Hobj.Clone("temp"))
            print(hname, nEvents, "+/-", error)
            den[hname.replace("h_njets_12incl_", "").replace("_ABCD","")] = (nEvents, error, Hobj, cABCD, eABCD)

        num = {}
        for hname, hinfo in self.signalRegions.items():
            newName = hname
            newInfo = copy.deepcopy(hinfo)
            newInfo["X"]["title"] = hinfo["X"]["title"]
            bname, binfo = self.mainBG, self.backgrounds[self.mainBG]
            Hobj = Histogram(None, rootFile, 1.0, 1.0, newName, newInfo, binfo)
            nEvents = Hobj.Integral()
            error = Hobj.IntegralError()
            cABCD, eABCD = self.getABCDCounts(Hobj.Clone("temp"))
            print(hname, nEvents, "+/-", error)
            num[hname.replace("h_njets_12incl_", "").replace("_ABCD","")] = (nEvents, error, Hobj, cABCD, eABCD)

        transferFactors = {}
        for nameNum, n in num.items():
            for nameDen, d in den.items():
                tfABCD, tfErrABCD = [],[]
                for i in range(len(n[3])):
                    tf, tfError = self.calculateTF(n[3][i], d[3][i], n[4][i], d[4][i])
                    tfABCD.append(tf)
                    tfErrABCD.append(tfError)
                tf, tfError = self.calculateTF(n[0], d[0], n[1], d[1])
                transferFactors["{}_TF_{}Over{}".format(self.year,nameNum,nameDen)] = (tf, tfError, n[2], d[2], tfABCD, tfErrABCD)

        for name, tf in transferFactors.items():
            h = ROOT.TH1D(name, name, 1, 0, 1)
            h.SetBinContent(1, tf[0])
            h.SetBinError(  1, tf[1])
            print(name, "TF:", round(tf[0],4), "+/-", round(tf[1],4))
            hABCD = ROOT.TH1D(name+"ABCD", name+"ABCD", 4, 0, 4)
            for i in range(len(tf[4])):
                print(tf[4][i], tf[5][i])
                hABCD.SetBinContent(i+1, tf[4][i])
                hABCD.SetBinError(  i+1, tf[5][i])
            self.transforFactorsHisto[name] = (h, tf[2], hABCD)

    def getTFPerBin(self, tfABCDHisto, i, repeat=6):
        tf, etf = 0.0, 0.0
        if   i <=   repeat: tf, etf = tfABCDHisto.GetBinContent(1), tfABCDHisto.GetBinError(1) 
        elif i <= 2*repeat: tf, etf = tfABCDHisto.GetBinContent(2), tfABCDHisto.GetBinError(2) 
        elif i <= 3*repeat: tf, etf = tfABCDHisto.GetBinContent(3), tfABCDHisto.GetBinError(3) 
        elif i <= 4*repeat: tf, etf = tfABCDHisto.GetBinContent(4), tfABCDHisto.GetBinError(4) 
        return tf, etf

    def write(self, crName):
        print("-----Save all output histograms-----")
        outfile = ROOT.TFile.Open("{}/{}_{}_Prediction.root".format(self.outpath, self.year, crName),"RECREATE")
        self.dataQCDOnly.Write()

        for name, h in self.transforFactorsHisto.items():
            h[0].Write()
            h[1].Write()
            h[2].Write()

        outfile.Close()

    def makeOutputPlots(self):
        for name, h in self.transforFactorsHisto.items():
            tfHisto = h[0]
            qcdSRMC = h[1]
            tfABCDHisto = h[2]
            hCRPred     = self.dataQCDOnly.Clone(""    )
            hCRPredUp   = self.dataQCDOnly.Clone("Up"  )
            hCRPredDown = self.dataQCDOnly.Clone("Down")
            hCRPred    .Scale(tfHisto.GetBinContent(1))
            hCRPredUp  .Scale(tfHisto.GetBinContent(1)+tfHisto.GetBinError(1))
            hCRPredDown.Scale(tfHisto.GetBinContent(1)-tfHisto.GetBinError(1))

            for iBin in range(1,hCRPred.GetNbinsX()+1):
                sigmaStat = hCRPred.GetBinError(iBin)/hCRPred.GetBinContent(iBin)
                sigmaTF   = (hCRPredUp.GetBinContent(iBin) - hCRPred.GetBinContent(iBin))/hCRPred.GetBinContent(iBin)
                sigmaTot = math.sqrt(sigmaStat**2 + sigmaTF**2)
                errorTot = sigmaTot*hCRPred.GetBinContent(iBin)
                hCRPred.SetBinError(iBin, errorTot)

            hCRPredABCD   = self.dataQCDOnly.Clone("ABCD")
            hCRPredABCD.SetMarkerColor(ROOT.kBlue)
            hCRPredABCD.SetLineColor(ROOT.kBlue)
            for iBin in range(1,hCRPredABCD.GetNbinsX()+1):
                tf, etf   = self.getTFPerBin(tfABCDHisto, iBin)
                val       = tf*hCRPredABCD.GetBinContent(iBin)
                valUp     = (tf+etf)*hCRPredABCD.GetBinContent(iBin)
                sigmaStat = tf*hCRPredABCD.GetBinError(iBin) / val
                sigmaTF   = (valUp - val) /val
                sigmaTot  = math.sqrt(sigmaStat**2 + sigmaTF**2)
                errorTot  = sigmaTot*val
                statError = tf*hCRPredABCD.GetBinError(iBin)
                hCRPredABCD.SetBinContent(iBin, val)
                hCRPredABCD.SetBinError(iBin, errorTot)

            canvas = ROOT.TCanvas("", "", 900, 900)
            split           = 0.3
            self.upperSplit = 1.0
            self.lowerSplit = 1.0
            self.scale      = 1.0
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
            ROOT.gPad.SetLogy()
   
            canvas.cd(2)
            ROOT.gPad.SetPad(0.0, 0.0, 1.0, split)
            ROOT.gPad.SetTopMargin(0)
            ROOT.gPad.SetBottomMargin(self.BottomMargin / self.lowerSplit)
            ROOT.gPad.SetLeftMargin(self.LeftMargin)
            ROOT.gPad.SetRightMargin(self.RightMargin)

            #cA = qcdSRMC.Integral(0,6)  /hCRPred.Integral(0,6)
            #cB = qcdSRMC.Integral(7,12) /hCRPred.Integral(7,12)
            #cC = qcdSRMC.Integral(13,18)/hCRPred.Integral(13,18)
            #cD = qcdSRMC.Integral(19,25)/hCRPred.Integral(19,25)
            #for i in range(1,self.dataQCDOnly.GetNbinsX()+1):
            #    s = cA if i <= 6 else (cB if 7<=i and i<=12 else ( cC if 13<=i and i<=18 else cD))
            #    x = hCRPred.GetBinContent(i)
            #    hCRPred.SetBinContent(i, s*x)

            canvas.cd(1)
            max = qcdSRMC.histogram.GetMaximum()
            qcdSRMC.histogram.SetMaximum(max*25.0)
            qcdSRMC.histogram.SetMinimum(2e-2)
            qcdSRMC.histogram.Draw()
            #hCRPredUp.Draw("same")
            #hCRPredDown.Draw("same")
            hCRPred.Draw("same")
            hCRPredABCD.Draw("same")
            print(hCRPred.Integral(), qcdSRMC.histogram.Integral())

            leg = ROOT.TLegend(0.2, 0.75, 0.5, 0.85)
            leg.SetBorderSize(0)
            leg.SetTextSize(0.05)
            leg.AddEntry(hCRPred,     "Scaled QCD CR Data", "p")
            leg.AddEntry(qcdSRMC.histogram,    "QCD SR MC",    "l")
            leg.AddEntry(hCRPredABCD, "ABCD Scaled QCD CR Data", "p")
            leg.Draw()

            canvas.cd(2)
            ROOT.gPad.SetGridy()
            ratio = hCRPredABCD.Clone("ratio")
            ratio.Divide(qcdSRMC.histogram)
            ratio2 = hCRPred.Clone("ratio")
            ratio2.Divide(qcdSRMC.histogram)
            ratio.SetMaximum(2.2)
            ratio.SetMinimum(0.0)
            ratio.GetYaxis().SetNdivisions(5, 5, 0)
            ratio.GetYaxis().SetTitle("QCD Pred/QCD MC")
            ratio.GetYaxis().SetTitleSize(0.1)
            ratio.GetXaxis().SetTitleSize(0.15)
            ratio.GetYaxis().SetTitleOffset(0.8)
            ratio.GetYaxis().SetLabelSize(0.1)
            ratio.GetXaxis().SetLabelSize(0.1)
            ratio.Draw()
            ratio2.Draw("same")

            self.addCMSlogo(canvas)

            canvas.SaveAs("%s/%s_%s.png"%(self.outpath, self.year, name))

if __name__ == "__main__":

    usage = "usage: %stackPlotter [options]"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("--year",         dest="year",         help="which year",                  default="Run2UL"               )
    parser.add_argument("--outpath",      dest="outpath",      help="Where to put plots",          default="output_QCDCRPrediction")
    parser.add_argument("--inpath",       dest="inpath",       help="Path to root files",          default="condor/Hadded_2016preVFP_DoubleDisco_22-10-24")
    args = parser.parse_args()

    controlRegions = [
        {"h_njets_12incl_QCDCR_ABCD"        : {"logY" : True, "orders" : list(xrange(1,2)), "Y" : {"title" : "Events / bin", "min" : 0.2}, "X" : {"title" : "N_{Jets} (ABCD bins)",    "rebin" : 1,  "min" : 0,  "max" :    24}},},
        {"h_njets_12incl_QCDCR_1b_ABCD"     : {"logY" : True, "orders" : list(xrange(1,2)), "Y" : {"title" : "Events / bin", "min" : 0.2}, "X" : {"title" : "N_{Jets} (ABCD bins)",    "rebin" : 1,  "min" : 0,  "max" :    24}},},
        {"h_njets_12incl_QCDCR_1t_ABCD"     : {"logY" : True, "orders" : list(xrange(1,2)), "Y" : {"title" : "Events / bin", "min" : 0.2}, "X" : {"title" : "N_{Jets} (ABCD bins)",    "rebin" : 1,  "min" : 0,  "max" :    24}},},
        {"h_njets_12incl_QCDCR_1b_1t_ABCD"  : {"logY" : True, "orders" : list(xrange(1,2)), "Y" : {"title" : "Events / bin", "min" : 0.2}, "X" : {"title" : "N_{Jets} (ABCD bins)",    "rebin" : 1,  "min" : 0,  "max" :    24}},},
        {"h_njets_12incl_QCDCR_2b_ABCD"     : {"logY" : True, "orders" : list(xrange(1,2)), "Y" : {"title" : "Events / bin", "min" : 0.2}, "X" : {"title" : "N_{Jets} (ABCD bins)",    "rebin" : 1,  "min" : 0,  "max" :    24}},},
    ]

    signalRegions = {
        "h_njets_12incl_0l_ABCD"  : {"logY" : True, "orders" : list(xrange(1,2)), "Y" : {"title" : "Events / bin", "min" : 0.2}, "X" : {"title" : "N_{Jets}",    "rebin" : 1,  "min" : 0,  "max" :    24}},
        #"h_njets_12incl_1l_ABCD"  : {"logY" : True, "orders" : list(xrange(1,2)), "Y" : {"title" : "Events / bin", "min" : 0.2}, "X" : {"title" : "N_{Jets}",    "rebin" : 1,  "min" : 0,  "max" :    24}},
        #"h_njets_12incl_2l_ABCD"  : {"logY" : True, "orders" : list(xrange(1,2)), "Y" : {"title" : "Events / bin", "min" : 0.2}, "X" : {"title" : "N_{Jets}",    "rebin" : 1,  "min" : 0,  "max" :    24}},
    }

    backgrounds = {
        "TT"       : {"name" : "t#bar{t} + jets", "color" : 40,  "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
        "QCD"      : {"name" : "QCD multijet",    "color" : ROOT.kRed,  "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
        "TTX"      : {"name" : "t#bar{t} + X",    "color" : 38,  "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
        "BG_OTHER" : {"name" : "Other",           "color" : 41,  "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    }

    data = {
        "Data" : {"name" : "Data", "color" : ROOT.kBlack, "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 1.3}
    }

    for cr in controlRegions:
        crName = cr.keys()[0].replace("h_njets_12incl_","").replace("_ABCD","")
        crProducer = ControlRegionProducer(args.year, args.outpath, args.inpath, cr, signalRegions, backgrounds, data)
        crProducer.getCRData()
        crProducer.getTransferFactors()
        crProducer.write(crName)
        crProducer.makeOutputPlots()
        print("-"*150+"\n")
