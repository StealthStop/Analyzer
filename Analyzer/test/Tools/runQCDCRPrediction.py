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

    def __init__(self, year, outpath, inpath, controlRegions, signalRegions, backgrounds, data, mainBG, inclBin, channel=None, model=None, sysName=""):
        self.year            = year
        self.outpath         = outpath
        self.inpath          = inpath        
        self.controlRegions  = controlRegions
        self.signalRegions   = signalRegions
        self.backgrounds     = backgrounds
        self.data            = data
        self.approved        = None
        self.channel         = channel
        self.model           = model
        self.sysName         = sysName
        
        if not os.path.exists(self.outpath):
            os.makedirs(self.outpath)

        self.mainBG = mainBG
        self.inclBin = inclBin
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
            mark.DrawLatex(self.LeftMargin + 0.12, 1 - (self.TopMargin - 0.017), "%sWork in Progress"%(simStr))

        mark.SetTextFont(42)
        mark.SetTextAlign(31)
        if "Run 2" in self.year:
            mark.DrawLatex(1 - self.RightMargin, 1 - (self.TopMargin - 0.017), "Run 2 (13 TeV)")
        else:
            mark.DrawLatex(1 - self.RightMargin, 1 - (self.TopMargin - 0.017), "%s (13 TeV)"%(self.year))

    def addExtraInfo(self, canvas, packedInfo):

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

        canvas.cd(1)

        text = ROOT.TLatex()
        text.SetNDC(True)

        text.SetTextAlign(31)
        text.SetTextSize(0.045)
        text.SetTextFont(42)
        text.SetTextColor(ROOT.TColor.GetColor("#7C99D1"))

        offset = 0.0
        if modelStr != "":
            text.DrawLatex(1 - self.RightMargin - 0.03, 1 - (self.TopMargin + 0.09 + offset), modelStr)
            offset += 0.05
        if channelStr != "":
            text.DrawLatex(1 - self.RightMargin - 0.03, 1 - (self.TopMargin + 0.09 + offset), channelStr)
            offset += 0.05
    
    def getCRData(self):
        print("-----Make small MC subtracted Data histograms from the CR-----")
        for hname, hinfo in self.controlRegions.items():
            newName = hname + self.sysName
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
                dataHist = Histogram(None, rootFile, 1.0, 1.0, hname, newInfo, dinfo)

            self.dataQCDOnly = dataHist.Clone("{}_Data_only_{}_{}{}".format(self.year, self.mainBG, hname.replace("h_njets_%s_"%(self.inclBin),"").replace("_ABCD",""),self.sysName))
            print("{}_Data_only_{}_{}{}".format(self.year, self.mainBG, hname.replace("h_njets_%s_"%(self.inclBin),"").replace("_ABCD",""),self.sysName))
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
        if n == 0 or d == 0:
            e = r
        else:
            e = r*math.sqrt( (nE/n)**2 + (dE/d)**2 )
        return (r, e)

    def getTransferFactors(self):
        print("-----Calcualte all transfer factors-----")
        rootFile = "%s/%s_%s.root"%(self.inpath, self.year, self.mainBG)

        qcdCRMChisto = None
        qcdSRMChisto = None
        den = {}
        for hname, hinfo in self.controlRegions.items():
            newName = hname + self.sysName
            newInfo = copy.deepcopy(hinfo)
            newInfo["X"]["title"] = hinfo["X"]["title"]
            bname, binfo = self.mainBG, self.backgrounds[self.mainBG]
            Hobj = Histogram(None, rootFile, 1.0, 1.0, newName, newInfo, binfo)
            nEvents = Hobj.Integral()
            error = Hobj.IntegralError()
            cABCD, eABCD = self.getABCDCounts(Hobj.Clone("temp"))
            print(newName, nEvents, "+/-", error)
            den[newName.replace("h_njets_%s_"%(self.inclBin), "").replace("_ABCD","")] = (nEvents, error, Hobj, cABCD, eABCD)

            qcdCRMChisto = Hobj.histogram

        num = {}
        for hname, hinfo in self.signalRegions.items():
            newName = hname + self.sysName
            newInfo = copy.deepcopy(hinfo)
            newInfo["X"]["title"] = hinfo["X"]["title"]
            bname, binfo = self.mainBG, self.backgrounds[self.mainBG]
            Hobj = Histogram(None, rootFile, 1.0, 1.0, newName, newInfo, binfo)
            nEvents = Hobj.Integral()
            error = Hobj.IntegralError()
            cABCD, eABCD = self.getABCDCounts(Hobj.Clone("temp"))
            print(newName, nEvents, "+/-", error)
            num[newName.replace("h_njets_%s_"%(self.inclBin), "").replace("_ABCD","")] = (nEvents, error, Hobj, cABCD, eABCD)

            qcdSRMChisto = Hobj.histogram

        # Only make QCD systematic when getting nominal histograms and not already varied ones
        tfShapeHisto = self.makeQCDshapeCorr(qcdCRMChisto, qcdSRMChisto)

        transferFactors = {}
        for nameNum, n in num.items():
            for nameDen, d in den.items():
                tfABCD, tfErrABCD, tfShape = [],[],[]
                for i in range(len(n[3])):
                    tf, tfError = self.calculateTF(n[3][i], d[3][i], n[4][i], d[4][i])
                    tfABCD.append(tf)
                    tfErrABCD.append(tfError)
                for ibin in range(1, tfShapeHisto.GetNbinsX()+1):
                    tfShape.append(tfShapeHisto.GetBinContent(ibin))
                tf, tfError = self.calculateTF(n[0], d[0], n[1], d[1])
                transferFactors["{}_TF_{}Over{}".format(self.year,nameNum,nameDen)] = (tf, tfError, n[2], d[2], tfABCD, tfErrABCD, tfShape)

        for name, tf in transferFactors.items():
            h = ROOT.TH1D(name, name, 1, 0, 1)
            h.SetBinContent(1, tf[0])
            h.SetBinError(  1, tf[1])
            print(name, "TF:", round(tf[0],4), "+/-", round(tf[1],4))
            hABCD = ROOT.TH1D(name+"ABCD", name+"ABCD", 4, 0, 4)
            tfPrint = ""

            shapeBins = tfShapeHisto.GetNbinsX()
            nbins = shapeBins*4
            hABCDshape = ROOT.TH1D(name+"ABCD_perNjets", name+"ABCD_perNjets", nbins, 0, nbins)

            for i in range(len(tf[4])):
                print(tf[4][i], tf[5][i])
                tfPrint += " {} \pm {} &".format(round(tf[4][i]*100, 1), round(tf[5][i]*100, 1))
                hABCD.SetBinContent(i+1, tf[4][i])
                hABCD.SetBinError(  i+1, tf[5][i])

                for j in range(1, tfShapeHisto.GetNbinsX()+1):
                    hABCDshape.SetBinContent(i*shapeBins + j, tfShapeHisto.GetBinContent(j))
                    hABCDshape.SetBinError(i*shapeBins + j,   tfShapeHisto.GetBinError(j))

            print(tfPrint)

            self.transforFactorsHisto[name] = (h, tf[2], hABCD, tf[3], hABCDshape)

    def getTFPerBin(self, tfABCDHisto, i, repeat=6):
        tf, etf = 0.0, 0.0
        if   i <=   repeat: tf, etf = tfABCDHisto.GetBinContent(1), tfABCDHisto.GetBinError(1) 
        elif i <= 2*repeat: tf, etf = tfABCDHisto.GetBinContent(2), tfABCDHisto.GetBinError(2) 
        elif i <= 3*repeat: tf, etf = tfABCDHisto.GetBinContent(3), tfABCDHisto.GetBinError(3) 
        elif i <= 4*repeat: tf, etf = tfABCDHisto.GetBinContent(4), tfABCDHisto.GetBinError(4) 
        return tf, etf

    def write(self, crName):
        print("-----Save all output histograms-----")
        outfileName = "{}/{}_{}_Prediction.root".format(self.outpath, self.year, crName)
        if "10_10" in crName:
            outfile = ROOT.TFile.Open(outfileName,"RECREATE")
        else:
            outfile = ROOT.TFile.Open(outfileName,"UPDATE")
        self.dataQCDOnly.Write()

        for name, h in self.transforFactorsHisto.items():
            h[0].Write()
            h[1].Write()
            h[2].Write()
            h[3].Write()
            h[4].Write()

        outfile.Close()
        return outfileName

    def write_QCD_Data(self):

        outfile = ROOT.TFile.Open(self.outpath + "/Run2UL_QCD_Data.root", "UPDATE")
        histName = "h_njets_{}_{}_{}_ABCD{}".format(self.inclBin, self.model, self.channel,self.sysName)
        self.QCDPred = ROOT.TH1D(histName, histName, 24, -0.5, 23.5)

        for bin in range(1,self.QCDPred.GetNbinsX()):
            tf, etf = self.getTFPerBin(self.transforFactorsHisto["{}_TF_{}_{}{}Over{}_{}_QCDCR{}".format(self.year, self.model, self.channel, self.sysName, self.model, self.channel, self.sysName)][2], bin)
        
            self.QCDPred.SetBinContent(bin, self.dataQCDOnly.GetBinContent(bin)*tf)
            self.QCDPred.SetBinError(bin, math.sqrt(self.dataQCDOnly.GetBinError(bin)**2+etf**2))

        self.QCDPred.Write()

        outfile.Close()

    def makeQCDshapeCorr(self, crHisto, srHisto):

        def sumBins(histo, name):
            histoCollapsedName = histo.GetName() + "_" + name + "_collapsed"
            nbins = histo.GetNbinsX() / 4
            histoCollapsed = ROOT.TH1F(histoCollapsedName, histoCollapsedName, histo.GetNbinsX()/4, 0, histo.GetNbinsX()/4)

            for ibin in range(1, nbins+1):

                content = 0.0
                error   = 0.0
                for jbin in range(4):
                    content += histo.GetBinContent(ibin + jbin * nbins)
                    error   += histo.GetBinError(ibin + jbin * nbins)**2.0

                histoCollapsed.SetBinContent(ibin, content)
                histoCollapsed.SetBinError(ibin, error**0.5)

            return histoCollapsed

        def getFitHisto(histo, hName, name):

            nbins = histo.GetNbinsX()

            func = "expo(0)"
            if self.channel == "0l":
                func = "[0] + expo(1)"

            fit = ROOT.TF1(name, func, 0.5, nbins-0.5)

            if self.channel == "0l":
                fit.SetParameter(0, -50.0)
                fit.SetParameter(1, 10.0)
                fit.SetParameter(2, -1.0)
            else:
                fit.SetParameter(0, 1.0)
                fit.SetParameter(1, -1.0)
               
            histo.Fit(name, "WLR")

            graphErrors = ROOT.TGraphErrors(nbins)
            for i in range(0, nbins):
                graphErrors.SetPoint(i+1, i+0.5, 0)
            ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(graphErrors, 0.68)

            fitHisto = histo.Clone(hName + "_fit"); fitHisto.Reset()
            nbins = fitHisto.GetNbinsX()
            for ibin in range(1, nbins+1):
                binX    = fitHisto.GetBinCenter(ibin)
                fitVal  = fit.Eval(binX)
                fitErr  = graphErrors.GetErrorY(ibin-1)
                fitHisto.SetBinContent(ibin, fitVal)
                fitHisto.SetBinError(ibin, fitErr)

            N = 100*nbins
            graphErrorsAux = ROOT.TGraphErrors(N)
            actualNpoints = 0
            for i in range(0, N):
                x = float(nbins)*(float(i)/float(N))
                if x < 0.49 or x > 4.5:
                    continue
                actualNpoints += 1
                graphErrorsAux.SetPoint(i, x, 0)
            ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(graphErrorsAux, 0.68)

            returnGraph = ROOT.TGraph(2*actualNpoints)

            j = 0
            for i in range(0, N):
                x1 = float(nbins)*(float(i)/float(N))
                if x1 < 0.49 or x1 > 4.5:
                    continue
                fitVal    = fit.Eval(x1)
                fitErr    = graphErrorsAux.GetErrorY(i)

                x2 = float(nbins)*(float(N-i-1)/float(N))
                if x2 < 0.49 or x2 > 4.5:
                    continue
                fitValAux = fit.Eval(x2)
                fitErrAux = graphErrorsAux.GetErrorY(N-i-1)

                returnGraph.SetPoint(j,               x1, fitVal + fitErr)
                returnGraph.SetPoint(actualNpoints+j, x2, fitValAux - fitErrAux)

                j += 1

            return fitHisto, fit, returnGraph

        crHistoCollapsed = sumBins(crHisto, "total")
        srHistoCollapsed = sumBins(srHisto, "total")

        crFitHisto, crFit, crGraph = getFitHisto(crHistoCollapsed, crHisto.GetName(), "crFit")
        srFitHisto, srFit, srGraph = getFitHisto(srHistoCollapsed, srHisto.GetName(), "srFit")

        qcdFitCorr = crFitHisto.Clone(crHisto.GetName() + "_fit_qcdShapeTF"); qcdFitCorr.Reset()
        qcdFitCorr.Divide(srFitHisto, crFitHisto)

        name = "QCD_Njets_Shape_Corr_%s_%s%s"%(self.model, self.channel, self.sysName)
        canvas = ROOT.TCanvas(name, name, 800, 800)
        canvas.Divide(1,2)

        split      = 0.3
        upperSplit = 1.0
        lowerSplit = 1.0
        scale      = 1.0
        RightMargin  = 0.1
        LeftMargin   = 0.1
        TopMargin    = 0.1
        BottomMargin = 0.1

        self.addCMSlogo(canvas)

        upperSplit = 1.0-split
        lowerSplit = split
        scale = upperSplit / lowerSplit

        canvas.cd(1)
        ROOT.gPad.SetPad(0.0, split, 1.0, 1.0)
        ROOT.gPad.SetTopMargin(self.TopMargin / upperSplit)
        ROOT.gPad.SetBottomMargin(0)
        ROOT.gPad.SetLeftMargin(self.LeftMargin)
        ROOT.gPad.SetRightMargin(self.RightMargin)
        ROOT.gPad.SetLogy()

        titleSize = 0.07
        labelSize = 0.06

        minimum = srHistoCollapsed.GetBinContent(srHistoCollapsed.GetNbinsX())
        maximum = crHistoCollapsed.GetBinContent(1)

        crHistoCollapsed.SetTitle("")
        crHistoCollapsed.GetYaxis().SetTitle("Number of Events")
        crHistoCollapsed.GetXaxis().SetTitle("N_{jets}")
        crHistoCollapsed.GetYaxis().SetTitleSize(titleSize)
        crHistoCollapsed.GetYaxis().SetLabelSize(labelSize)
        crHistoCollapsed.GetYaxis().SetTitleOffset(0.9)
        crHistoCollapsed.GetYaxis().SetRangeUser(minimum/10.0, maximum*10.0)
        crHistoCollapsed.SetLineWidth(4)
        crHistoCollapsed.SetLineColor(30)
    
        srHistoCollapsed.SetLineWidth(4)
        srHistoCollapsed.SetLineColor(ROOT.kGreen+3)

        crHistoCollapsed.Draw("EHIST")
        srHistoCollapsed.Draw("EHIST SAME")

        iamLegend = ROOT.TLegend(0.50, 0.70, 0.70, 0.9)
        iamLegend.AddEntry(crHistoCollapsed, "QCD^{MC}_{CR}", "EL")
        iamLegend.AddEntry(srHistoCollapsed, "QCD^{MC}_{SR}", "EL")
        iamLegend.SetTextSize(0.06)

        crFit.SetLineWidth(4)
        crFit.SetLineColor(30)
        crFit.SetLineStyle(7)
        srFit.SetLineWidth(4)
        srFit.SetLineColor(ROOT.kGreen+3)
        srFit.SetLineStyle(3)

        crGraph.SetFillColorAlpha(30, 0.5)
        srGraph.SetFillColorAlpha(ROOT.kGreen+3, 0.3)

        crFit.Draw("LSAME")
        srFit.Draw("LSAME")

        crGraph.Draw("F SAME")
        srGraph.Draw("F SAME")

        iamLegend.Draw("SAME")

        self.addExtraInfo(canvas, self.channel)

        canvas.cd(2)
        ROOT.gPad.SetPad(0.0, 0.0, 1.0, split)
        ROOT.gPad.SetTopMargin(0)
        ROOT.gPad.SetBottomMargin(self.BottomMargin / lowerSplit)
        ROOT.gPad.SetLeftMargin(self.LeftMargin)
        ROOT.gPad.SetRightMargin(self.RightMargin)
        ROOT.gPad.SetGridy()

        line = ROOT.TLine(0, 1, qcdFitCorr.GetXaxis().GetXmax(), 1)
        line.SetLineColor(ROOT.kBlack)
        line.SetLineStyle(2)
        line.SetLineWidth(2)

        qcdFitCorr.GetYaxis().SetTitle("TF(SR/CR)")
        qcdFitCorr.GetXaxis().SetTitle("N_{ jets}")
        qcdFitCorr.GetYaxis().SetTitleSize(scale*titleSize)
        qcdFitCorr.GetYaxis().SetLabelSize(scale*labelSize)
        qcdFitCorr.GetXaxis().SetTitleSize(scale*titleSize)
        qcdFitCorr.GetYaxis().SetTitleOffset(1.1 / scale)
        qcdFitCorr.GetXaxis().SetTitleOffset(1.1)
        qcdFitCorr.GetYaxis().SetNdivisions(5, 5, 0)

        startNjet = 6
        endNjet = 10
        if self.channel == "1l":
            startNjet = 7
            endNjet = 11
        elif self.channel == "0l":
            startNjet = 8
            endNjet = 12

        qcdFitCorr.GetXaxis().SetLabelSize(scale*labelSize*1.6)
        qcdFitCorr.GetXaxis().SetLabelOffset(0.05 / scale)

        for ibin in range(1, qcdFitCorr.GetNbinsX()+1):

            label = ""
            Njet = startNjet + ibin - 1
            if Njet == endNjet:
                label = "#geq %d"%(Njet)
            else:
                label = str(Njet)

            qcdFitCorr.GetXaxis().SetBinLabel(ibin, label)

        upperBound = 0.85
        lowerBound = -0.05
        if self.channel == "2l":
            upperBound = 0.012
            lowerBound = -0.0015
        elif self.channel == "1l":
            upperBound = 0.14
            lowerBound = -0.01
        qcdFitCorr.GetYaxis().SetRangeUser(lowerBound, upperBound)
        qcdFitCorr.SetLineColor(ROOT.kBlue)
        qcdFitCorr.SetTitle("")
        qcdFitCorr.SetLineWidth(3)

        qcdFitCorrErr = qcdFitCorr.Clone(qcdFitCorr.GetName() + "_justError")
        qcdFitCorrErr.SetFillColorAlpha(ROOT.kBlue, 0.3)

        qcdFitCorr.Draw("HIST")
        qcdFitCorrErr.Draw("E2SAME")
        line.Draw("SAME")

        canvas.SaveAs(self.outpath + "/QCD_Shape_TF_%s_%s.pdf"%(self.model, self.channel))

        return qcdFitCorr

    def makeOutputPlots(self, reducedPlot = True, systHisto = None):
        for name, h in self.transforFactorsHisto.items():
            #print(name)
            #print(h)
            tfHisto = h[0]
            qcdSRMC = h[1]
            tfABCDHisto = h[2]
            qcdCRMC = h[3]
            tfABCDshapeHisto = h[4]
            hCRPred     = self.dataQCDOnly.Clone(""    )
            hCRPredUp   = self.dataQCDOnly.Clone("Up"  )
            hCRPredDown = self.dataQCDOnly.Clone("Down")
            hCRPred    .Scale(tfHisto.GetBinContent(1))
            hCRPredUp  .Scale(tfHisto.GetBinContent(1)+tfHisto.GetBinError(1))
            hCRPredDown.Scale(tfHisto.GetBinContent(1)-tfHisto.GetBinError(1))

            for iBin in range(1,hCRPred.GetNbinsX()+1):
                #print(iBin, hCRPred.GetBinError(iBin), hCRPred.GetBinContent(iBin))
                sigmaStat = 0.0
                sigmaTF = 0.0
                if hCRPred.GetBinContent(iBin)>0.0:
                    sigmaStat = hCRPred.GetBinError(iBin)/hCRPred.GetBinContent(iBin)
                    sigmaTF   = (hCRPredUp.GetBinContent(iBin) - hCRPred.GetBinContent(iBin))/hCRPred.GetBinContent(iBin)
                sigmaTot = math.sqrt(sigmaStat**2 + sigmaTF**2)
                errorTot = sigmaTot*hCRPred.GetBinContent(iBin)
                hCRPred.SetBinError(iBin, errorTot)

            hCRPredABCD   = self.dataQCDOnly.Clone("ABCD")
            hCRPredABCD.SetMarkerColor(ROOT.TColor.GetColor("#006d2c"))
            hCRPredABCD.SetMarkerStyle(8)
            hCRPredABCD.SetMarkerSize(2)
            hCRPredABCD.SetLineColor(ROOT.TColor.GetColor("#006d2c"))
            #hCRPredABCD.Print("ALL")
            for iBin in range(1,hCRPredABCD.GetNbinsX()+1):

                # NOTE: Using single TF here, not per-ABCD TFs
                #tf, etf    = self.getTFPerBin(tfHisto, 1)
                njetShift  = tfABCDshapeHisto.GetNbinsX()+1
                tfShape    = tfABCDshapeHisto.GetBinContent(iBin % njetShift)
                tfShapeErr = tfABCDshapeHisto.GetBinError(iBin % njetShift)

                val       = tfShape*hCRPredABCD.GetBinContent(iBin)
                valUp     = (tfShape+tfShapeErr)*hCRPredABCD.GetBinContent(iBin)
                sigmaStat = 0.0
                sigmaTF = 0.0
                if val > 0.0:
                    sigmaStat = tfShape*hCRPredABCD.GetBinError(iBin) / val
                    sigmaTF   = (valUp - val) /val
                sigmaTot  = math.sqrt(sigmaStat**2 + sigmaTF**2)
                errorTot  = sigmaTot*val
                statError = tfShape*hCRPredABCD.GetBinError(iBin)
                hCRPredABCD.SetBinContent(iBin, val)
                hCRPredABCD.SetBinError(iBin, errorTot)

                if systHisto != None:
                    systHisto.SetBinContent(iBin, 1.0)
                    systHisto.SetBinError(iBin, (systHisto.GetBinError(iBin)**2.0 + (tfShapeErr/tfShape)**2.0)**0.5)

                qcdCRMC.histogram.SetBinContent(iBin, tfShape*qcdCRMC.histogram.GetBinContent(iBin))

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

            canvas.cd(1)
            max = qcdSRMC.histogram.GetMaximum()
            qcdSRMC.histogram.SetMaximum(max*50.0)
            qcdSRMC.histogram.SetMinimum(2e-2)
            qcdSRMC.histogram.SetMarkerStyle(21)
            qcdSRMC.histogram.SetMarkerSize(2)
            qcdSRMC.histogram.GetYaxis().SetTitleOffset(1.2)
            qcdSRMC.histogram.GetXaxis().SetRangeUser(-0.5, 19.5)
            qcdSRMC.histogram.Draw()
            #hCRPredUp.Draw("same")
            #hCRPredDown.Draw("same")
            #hCRPred.Draw("same")
            hCRPredABCD.Draw("same")
            #print(hCRPred.Integral(), qcdSRMC.histogram.Integral())

            if not reducedPlot:
                qcdCRMC.histogram.SetLineColor(ROOT.kBlack)
                qcdCRMC.histogram.Draw("same hist")

            # get model, channel labels
            self.addExtraInfo(canvas,name)         

            leg = ROOT.TLegend(0.2, 0.65, 0.5, 0.9)
            leg.SetBorderSize(0)
            leg.SetTextSize(0.05)
            #leg.AddEntry(hCRPred,     "Scaled QCD CR Data", "p")
            leg.AddEntry(qcdSRMC.histogram, "QCD_{MC}^{SR}",                        "lp")
            if not reducedPlot:
                leg.AddEntry(qcdCRMC.histogram, "QCD_{MC}^{CR} (Scaled by TF)",         "lp")
            leg.AddEntry(hCRPredABCD,       "QCD_{Data}^{Pred.} (Scaled by TF)",    "lp")
            leg.Draw()

            canvas.cd(2)
            ROOT.gPad.SetGridy()
            ratio = hCRPredABCD.Clone(hCRPredABCD.GetName() + "_ratio")
            ratio.Reset()
            ratio.Divide(qcdSRMC.histogram, hCRPredABCD)

            ratio2 = hCRPred.Clone(hCRPred.GetName() + "_ratio")
            ratio2.Reset()
            ratio2.Divide(qcdSRMC.histogram, hCRPred)

            ratio3 = qcdSRMC.histogram.Clone(qcdSRMC.histogram.GetName() + "_ratio")
            ratio3.Reset()
            ratio3.Divide(qcdSRMC.histogram, qcdCRMC.histogram)

            ratio.SetMaximum(2.2)
            ratio.SetMinimum(0.0)
            ratio.SetMarkerColor(ROOT.kBlack)
            #ratio.SetMarkerSize(1)
            ratio.SetLineColor(ROOT.kBlack)
            ratio.GetYaxis().SetNdivisions(5, 5, 0)
            ratio.GetYaxis().SetTitle("QCD_{MC}^{SR} / QCD_{Data}^{CR}")
            ratio.GetYaxis().SetTitleSize(0.1)
            ratio.GetXaxis().SetTitleSize(0.12)
            ratio.GetYaxis().SetTitleOffset(0.6)
            ratio.GetYaxis().SetLabelSize(0.1)
            ratio.GetXaxis().SetLabelSize(0.12)
            ratio.GetXaxis().SetRangeUser(-0.5, 19.5)
            ratio.Draw()

            if systHisto != None:
                systHisto.SetFillColor(ROOT.kGray+2)
                systHisto.SetFillStyle(3354)
                systHisto.Draw("E2 SAME")
            ratio.Draw("SAME")
            if not reducedPlot:
                #ratio2.Draw("same")
                ratio3.Draw("same hist P")

            self.addCMSlogo(canvas)

            canvas.SaveAs("%s/%s_%s.png"%(self.outpath, self.year, name))
            canvas.SaveAs("%s/%s_%s.pdf"%(self.outpath, self.year, name))

if __name__ == "__main__":

    usage = "usage: %stackPlotter [options]"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("--year",         dest="year",         help="which year",                  default="Run2UL"               )
    parser.add_argument("--outpath",      dest="outpath",      help="Where to put plots",          default="output_QCDCRPrediction")
    parser.add_argument("--inpath",       dest="inpath",       help="Path to root files",          default="/uscms/home/jhiltb/nobackup/PO_Boxes/shared/DisCo_outputs_0l_1l_2l_optimizedBinEdges_forSYY_10.12.2022")
    parser.add_argument("--channel",      dest="channel",      help="0l, 1l, 2l",                  default="1l")
    parser.add_argument("--model",        dest="model",        help="Which NN model (RPV, SYY)",   default="RPV")
    args = parser.parse_args()

    if   args.channel == "2l":
        inclBin = "10incl"
    elif args.channel == "1l":
        inclBin = "11incl"
    elif args.channel == "0l":
        inclBin = "12incl"

    controlRegions = [
        {"h_njets_%s_%s_%s_QCDCR_ABCD"%(inclBin, args.model, args.channel) : {"logY" : True, "orders" : list(xrange(1,2)), "Y" : {"title" : "Events / bin", "min" : 0.2}, "X" : {"title" : "N_{Jets} in each A,B,C,D region", "rebin" : 1,  "min" : 0,  "max" : 24}},},
    ]

    signalRegions = {
        "h_njets_%s_%s_%s_ABCD"%(inclBin, args.model, args.channel) : {"logY" : True, "orders" : list(xrange(1,2)), "Y" : {"title" : "Events / bin", "min" : 0.2}, "X" : {"title" : "N_{Jets}", "rebin" : 1, "min" : 0, "max" : 24}},
    }

    backgrounds = {
        "TT"       : {"name" : "t#bar{t} + jets", "color" : 40,                              "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
        "QCD"      : {"name" : "QCD multijet",    "color" : ROOT.TColor.GetColor("#85c2a3"), "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
        "TTX"      : {"name" : "t#bar{t} + X",    "color" : 38,                              "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
        "BG_OTHER" : {"name" : "Other",           "color" : 41,                              "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 0},
    }

    data = {
        "Data" : {"name" : "Data", "color" : ROOT.kBlack, "lstyle" : 1, "mstyle" : 8, "lsize" : 3, "msize" : 1.3}
    }

    mainBG = "QCD" if "QCD" in backgrounds else "QCD_skim"

    sys = ["", "_JECup", "_JECdown", "_JERup", "_JERdown", "_fsrUp", "_fsrDown", "_isrUp", "_isrDown", "_pdfUp", "_pdfDown", "_prfUp", "_prfDown", "_puUp", "_puDown", "_sclUp", "_sclDown"]
    #if args.channel == "0l":
    #    sys += ["_jetUp", "_jetDown"]
    #else:
    #    sys += ["_lepUp", "_lepDown"]

    colors = [ROOT.kBlack, ROOT.kCyan+1, ROOT.kCyan+2, ROOT.kBlue-7, ROOT.kBlue-5, ROOT.kPink-9, ROOT.kPink-8, ROOT.kOrange+4, ROOT.kOrange+5, ROOT.kBlue-7, ROOT.kBlue-5, ROOT.kMagenta+3, ROOT.kMagenta-1, ROOT.kRed, ROOT.kRed+2, ROOT.kGreen-6, ROOT.kGreen+2]

    names = ["nominal", "JEC Up", "JEC Down", "JER Up", "JER Down", "FSR Up", "FSR Down", "ISR Up", "ISR Down", "PDF Up", "PDF Down", "Prf. Up", "Prf. Down", "PU Up", "PU Down", "Scl. Up", "Scl. Down"]

    outfileNames = []
    crProducers = {}
    for sysName in sys:
        for cr in controlRegions:
            print(sysName)
            crName = cr.keys()[0].replace("h_njets_%s_"%(inclBin),"").replace("_ABCD",sysName)
            crProducer = ControlRegionProducer(args.year, args.outpath, args.inpath, cr, signalRegions, backgrounds, data, mainBG, inclBin, args.channel, args.model, sysName)
            crProducer.getCRData()
            crProducer.getTransferFactors()
            if sysName != "": crProducer.makeOutputPlots()
            outfileNames.append(crProducer.write(crName))
            crProducer.write_QCD_Data()
            crProducers[sysName] = crProducer
            print("\t\t\n________________\t\t\n")

    def makeSystPlot(producers, model, channel, outpath):

        canvas = ROOT.TCanvas("", "", 900, 900)
        canvas.cd()
        canvas.SetTopMargin(0.06)
        canvas.SetBottomMargin(0.12)
        canvas.SetRightMargin(0.04)
        canvas.SetLeftMargin(0.16)

        nominal = producers[""].dataQCDOnly.Clone(producers[""].dataQCDOnly.GetName() + "_temp")
        nominal.GetXaxis().SetTitle("N_{ jets} each A,B,C,D region")
        nominal.GetYaxis().SetTitle("A.U.")
        nominal.GetYaxis().SetRangeUser(0.8, 1.25)
        nominal.GetXaxis().SetRangeUser(-0.5, 19.5)
        nominal.SetLineWidth(0)
        nominal.SetMarkerSize(0)
        nominal.Draw("HIST")

        iamLegend = ROOT.TLegend(0.16, 0.7, 0.96, 0.96)
        iamLegend.SetNColumns(3)
        
        garbage = {}
        for sysName in sys:
            if sysName == "": continue

            variation = producers[sysName].dataQCDOnly.Clone(producers[sysName].dataQCDOnly.GetName() + "_temp")

            systematic = variation.Clone(variation.GetName() + "_systematic")
            systematic.Reset()

            systematic.Divide(variation, nominal)

            systematic.SetLineWidth(3)
            systematic.SetMarkerSize(0)
            systematic.SetMarkerStyle(20)
            systematic.SetLineColor(colors[sys.index(sysName)])
            iamLegend.AddEntry(systematic, names[sys.index(sysName)], "L")
            systematic.Draw("SAME HIST ][")

            garbage[sysName] = systematic

        iamLegend.Draw("SAME")

        canvas.SaveAs("%s/Transfered_QCD_Variations_%s_%s.pdf"%(outpath, model, channel))

        quadSyst = [0.0 for n in range(producers[""].dataQCDOnly.GetNbinsX())]
        for sysName in garbage.keys():
            if "up" not in sys and "Up" not in sysName: continue

            sysName.replace("Up", "").replace("up", "")

            upHisto   = garbage[sysName]
            downHisto = garbage[sysName.replace("up", "down").replace("Up", "Down")]

            for ibin in range(1, upHisto.GetNbinsX()+1):
                isyst = max(abs(1.0-upHisto.GetBinContent(ibin)), abs(1.0-downHisto.GetBinContent(ibin)))

                quadSyst[ibin-1] += isyst**2.0

        totalSystHisto = ROOT.TH1F("totalSystHisto", "totalSystHisto", len(quadSyst), 0, len(quadSyst))

        for ival in range(len(quadSyst)):
            totalSystHisto.SetBinContent(ival+1, 1.0)
            totalSystHisto.SetBinError(ival+1, quadSyst[ival]**0.5)

        return totalSystHisto

    totalSystHisto = makeSystPlot(crProducers, args.model, args.channel, args.outpath)

    crProducers[""].makeOutputPlots(reducedPlot = True, systHisto = totalSystHisto)

    ##########
    # Finally hadd all the output systematic variations together to have one root file for combine data card
    finalHadd = "hadd -f {}/Total_{} {}".format(args.outpath, outfileNames[0].replace(args.outpath+"/",""), " ".join(outfileNames))
    print(finalHadd)
    os.system(finalHadd)
