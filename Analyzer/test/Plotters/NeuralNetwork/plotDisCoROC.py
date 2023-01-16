import os, ROOT, argparse

ROOT.TH1.SetDefaultSumw2()
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat("")
ROOT.gStyle.SetPaintTextFormat("3.2f")
ROOT.gStyle.SetFrameLineWidth(2)

class RocAndRoll:

    def __init__(self, approved, inputDir, outputDir, year, channel, background, signal):
    
        self.approved   = approved
        self.inputDir   = inputDir
        self.outputDir  = outputDir
        self.year       = year
        self.channel    = channel
        self.background = background

        self.TopMargin    = 0.08
        self.LeftMargin   = 0.13
        self.BottomMargin = 0.12
        self.RightMargin  = 0.04

        self.signal = signal
        if signal[0] == "S":
            self.signal = "Stealth" + signal

        self.signalShort = signal

        minNjet = 7; maxNjet = 12
        if channel == "0l":
            minNjet = 8
            maxNjet = 13
        elif channel == "2l":
            minNjet = 6
            maxNjet = 11

        self.Njets = [str(Njet) if Njet < maxNjet else str(Njet)+"incl" for Njet in range(minNjet, maxNjet+1)]
        self.masses = ["350", "450", "550", "650", "850", "1150"]
        self.massColors = {"350" : ROOT.kRed,       "450" : ROOT.kGreen+2,  "550"  : ROOT.kBlack,
                           "650" : ROOT.kMagenta+2, "850" : ROOT.kOrange+1, "1150" : ROOT.kBlue
        }

        colors = [ROOT.kRed,
                  ROOT.kGreen+2,
                  ROOT.kBlack,
                  ROOT.kMagenta+2,
                  ROOT.kOrange+1,
                  ROOT.kBlue
        ]

        self.njetsColors = {}
        for i in range(0, len(self.Njets)):
            self.njetsColors[self.Njets[i]] = colors[i]

        self.getHistograms()

        self.ROCs = {1 : {Njet : {} for Njet in self.Njets},
                     2 : {Njet : {} for Njet in self.Njets}
        }

        self.AUCs = {1 : {Njet : {} for Njet in self.Njets},
                     2 : {Njet : {} for Njet in self.Njets}
        }

        self.getROCs(1)
        self.getROCs(2)

    def printSignal(self, signal):

        if "SYY" in signal:
            return "Stealth SY#bar{Y}"
        else:
            return signal

    def printLabels(self, c, disc, text):
        c.cd()
        label = ROOT.TLatex()
        label.SetNDC(True)
        label.SetTextAlign(13)
        label.SetTextFont(42)
        label.SetTextSize(0.040)
        label.SetTextColor(ROOT.TColor.GetColor("#7C99D1"))

        channelStr = ""
        if "0l" in self.channel:
            channelStr = "Fully-Hadronic"
        elif "1l" in self.channel:
            channelStr = "Semi-Leptonic"
        elif "2l" in self.channel:
            channelStr = "Fully-Leptonic"

        modelStr = ""
        if "RPV" in self.signal:
            modelStr = "RPV"
        else:
            modelStr = "Stealth SYY"

        offset = 0.35
        xpos = 0.70
        label.DrawLatex(xpos, offset, "Disc. %d"%(disc))
        offset -= 0.05

        if modelStr != "":
            label.DrawLatex(xpos, offset, modelStr)
            offset -= 0.05
        if channelStr != "":
            label.DrawLatex(xpos, offset, channelStr)
            offset -= 0.05

        label.DrawLatex(xpos, offset, text)

        label.Draw("SAME")

    def drawCMS(self, c):
        c.cd()

        mark = ROOT.TLatex()
        mark.SetNDC(True)

        mark.SetTextAlign(11)
        mark.SetTextSize(0.055)
        mark.SetTextFont(61)
        mark.DrawLatex(self.LeftMargin, 1 - (self.TopMargin - 0.015), "CMS")
        mark.SetTextFont(52)
        mark.SetTextSize(0.040)
        if self.approved:
            mark.DrawLatex(self.LeftMargin + 0.12, 1 - (self.TopMargin - 0.016), "Simulation Supplementary")
        else:
            mark.DrawLatex(self.LeftMargin + 0.12, 1 - (self.TopMargin - 0.016), "Work in Progress")

        mark.SetTextFont(42)
        mark.SetTextAlign(31)
        mark.DrawLatex(1 - self.RightMargin, 1 - (self.TopMargin - 0.015), "%s (13 TeV)"%(self.year))

    def getHistograms(self):
    
        self.histograms = {Njet : {} for Njet in self.Njets} 
    
        f = ROOT.TFile.Open(self.inputDir + "/%s_%s.root"%(self.year, self.background))
    
        for Njet in self.Njets:
    
            self.histograms[Njet][self.background] = f.Get("h_DoubleDisCo_%s_disc1_disc2_%s_Njets%s_ABCD"%(self.signalShort,self.channel,Njet))
            self.histograms[Njet][self.background].SetDirectory(0)
    
        f.Close()
    
        for mass in self.masses:
            g = ROOT.TFile.Open(self.inputDir + "/%s_%s_2t6j_mStop-%s.root"%(self.year, self.signal, mass))
           
            for Njet in self.Njets:
                self.histograms[Njet][self.signal+mass] = g.Get("h_DoubleDisCo_%s_disc1_disc2_%s_Njets%s_ABCD"%(self.signalShort,self.channel,Njet)) 
                self.histograms[Njet][self.signal+mass].SetDirectory(0)
    
    def getROCs(self, disc):
    
        for Njet in self.Njets: 
            for mass in self.masses:

                bkg1D = None; sig1D = None; nBins = None
                minX = None; maxX = None
                if disc == 1:
                    bkg1D = self.histograms[Njet][self.background].ProjectionX("proj_%s_%s_%d_%s"%(Njet,self.background,disc,self.channel))
                    sig1D = self.histograms[Njet][self.signal+mass].ProjectionX("proj_%s_%s_%s_%d_%s"%(Njet,self.signal,mass,disc,self.channel))
                    nBins = bkg1D.GetNbinsX()
                else:
                    bkg1D = self.histograms[Njet][self.background].ProjectionY("proj_%s_%s_%d_%s"%(Njet,self.background,disc,self.channel))
                    sig1D = self.histograms[Njet][self.signal+mass].ProjectionY("proj_%s_%s_%s_%d_%s"%(Njet,self.signal,mass,disc,self.channel))
                    nBins = bkg1D.GetNbinsX()
   
                minX = bkg1D.GetXaxis().GetBinLowEdge(1)
                maxX = bkg1D.GetXaxis().GetBinUpEdge(nBins+1)

                bkg1D.Scale(1.0/bkg1D.Integral())
                sig1D.Scale(1.0/sig1D.Integral())

                newBins = 10 * nBins 
    
                name = "%s_%s_%s_disc%d_%s"%(self.signal,mass,Njet,disc,self.channel)
                roc = ROOT.TH1D(name, ";Background Efficiency;Signal Efficiency", newBins, minX, maxX)
                auc = 0.0
    
                fpSum = 1.0; tpSum = 1.0; currentBin = roc.GetNbinsX()
                fnSum = 0.0; tnSum = 0.0
                for iBin in range(1, bkg1D.GetNbinsX()+1):
                    fpSum -= bkg1D.GetBinContent(iBin)
                    tpSum -= sig1D.GetBinContent(iBin)
                    fnSum += sig1D.GetBinContent(iBin)
                    tnSum += bkg1D.GetBinContent(iBin)

                    xBin = roc.GetXaxis().FindBin(fpSum / (fpSum + tnSum))

                    for jBin in range(xBin, currentBin):
                        roc.SetBinContent(jBin, tpSum / (tpSum + fnSum))
                        auc += tpSum / (tpSum + fnSum) * (1.0/float(newBins))

                    currentBin = xBin
    
                roc.GetYaxis().SetLabelSize(0.040); roc.GetYaxis().SetTitleSize(0.050); roc.GetYaxis().SetTitleOffset(1.1)
                roc.GetXaxis().SetLabelSize(0.040); roc.GetXaxis().SetTitleSize(0.050); roc.GetXaxis().SetTitleOffset(1.0)
                roc.GetZaxis().SetLabelSize(0.040); roc.GetZaxis().SetTitleSize(0.050)
    
                roc.GetYaxis().SetRangeUser(0.0, 1.2)
                roc.GetXaxis().SetRangeUser(0.0, 1.0)
    
                roc.SetLineWidth(3)
                roc.SetLineStyle(0)
                roc.SetMarkerSize(0)
                roc.SetMarkerStyle(4)

                self.AUCs[disc][Njet][self.signal+mass] = auc
                self.ROCs[disc][Njet][self.signal+mass] = roc
    
    def drawPerNjetROC(self, disc):
    
        for Njet in self.Njets:
            
            c = ROOT.TCanvas("c%s_%s"%(Njet,disc), "c%s_%s"%(Njet,disc), 1000, 1000)
    
            c.cd()

            ROOT.gPad.SetTopMargin(self.TopMargin)
            ROOT.gPad.SetLeftMargin(self.LeftMargin)
            ROOT.gPad.SetBottomMargin(self.BottomMargin)
            ROOT.gPad.SetRightMargin(self.RightMargin)

            iamLegend = ROOT.TLegend(self.LeftMargin, (1.0-self.TopMargin+0.2*self.BottomMargin) / 1.2, 1.0-self.RightMargin, 1.0-self.TopMargin)
            iamLegend.SetNColumns(3)
            iamLegend.SetTextSize(0.028)

            jamLegend = ROOT.TLegend(0.65, self.BottomMargin+0.03, 1.0-self.RightMargin-0.03, 0.5)
            jamLegend.SetNColumns(1)
            jamLegend.SetTextSize(0.039)

            for mass in self.masses:

                self.ROCs[disc][Njet][self.signal+mass].SetLineColor(self.massColors[mass])
                self.ROCs[disc][Njet][self.signal+mass].SetMarkerColor(self.massColors[mass])
    
                iamLegend.AddEntry(self.ROCs[disc][Njet][self.signal+mass], "m_{ #tilde{t}} = %s GeV"%(mass)) 
                jamLegend.AddEntry(self.ROCs[disc][Njet][self.signal+mass], "AUC = %.3f"%(self.AUCs[disc][Njet][self.signal+mass]))

                self.ROCs[disc][Njet][self.signal+mass].Draw("SAME ][")

            l = ROOT.TLine(0.0, 0.0, 1.0, 1.0)
            l.SetLineWidth(2)
            l.SetLineColor(ROOT.kBlack)
            l.SetLineStyle(2)
            l.Draw("SAME")
    
            self.drawCMS(c)

            njetLabel = None
            if "incl" in Njet:
                njetLabel = "N_{ jets} #geq %s"%(Njet.replace("incl",""))
            else:
                njetLabel = "N_{ jets} = %s"%(Njet)

            self.printLabels(c, disc, njetLabel)

            iamLegend.Draw("SAME")
            #jamLegend.Draw("SAME")
    
            c.SaveAs(self.outputDir + "/ROC_%s_vs_TT_Njets%s_Disc%d_%s.pdf"%(self.signal,Njet,disc,self.channel))

            del c
    
    def drawPerMassROC(self, disc):
    
        for mass in self.masses:
        
            c = ROOT.TCanvas("c%s_%s"%(mass,disc), "c%s_%s"%(mass,disc), 1000, 1000)

            ROOT.gPad.SetTopMargin(self.TopMargin)
            ROOT.gPad.SetLeftMargin(self.LeftMargin)
            ROOT.gPad.SetBottomMargin(self.BottomMargin)
            ROOT.gPad.SetRightMargin(self.RightMargin)

            iamLegend = ROOT.TLegend(self.LeftMargin, (1.0-self.TopMargin+0.2*self.BottomMargin) / 1.2, 1.0-self.RightMargin, 1.0-self.TopMargin)
            iamLegend.SetNColumns(3)
            iamLegend.SetTextSize(0.034)

            jamLegend = ROOT.TLegend(0.65, self.BottomMargin+0.03, 1.0-self.RightMargin-0.03, 0.5)
            jamLegend.SetNColumns(1)
            jamLegend.SetTextSize(0.039)

            for Njet in self.Njets:
            
                self.ROCs[disc][Njet][self.signal+mass].SetLineColor(self.njetsColors[Njet]) 
                self.ROCs[disc][Njet][self.signal+mass].SetMarkerColor(self.njetsColors[Njet])

                if "incl" in Njet:
                    iamLegend.AddEntry(self.ROCs[disc][Njet][self.signal+mass], "N_{ jets} #geq %s"%(Njet.replace("incl","")))
                else:
                    iamLegend.AddEntry(self.ROCs[disc][Njet][self.signal+mass], "N_{ jets} = %s"%(Njet))

                jamLegend.AddEntry(self.ROCs[disc][Njet][self.signal+mass], "AUC = %.3f"%(self.AUCs[disc][Njet][self.signal+mass]))

                self.ROCs[disc][Njet][self.signal+mass].Draw("SAME ][")

            l = ROOT.TLine(0.0, 0.0, 1.0, 1.0)
            l.SetLineWidth(2)
            l.SetLineColor(ROOT.kBlack)
            l.SetLineStyle(2)
            l.Draw("SAME")

            self.drawCMS(c)
            self.printLabels(c, disc, "m_{ #tilde{t}} = %s GeV"%(mass))

            iamLegend.Draw("SAME")
            #jamLegend.Draw("SAME")

            c.SaveAs(self.outputDir + "/ROC_%s%s_vs_TT_Disc%d_%s.pdf"%(self.signal,mass,disc,self.channel))
    
            del c

if __name__ == "__main__":
    usage = "usage: %prog [options]"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("--approved",  dest="approved",  help="Plot is approved",     action="store_true", default=False) 
    parser.add_argument("--year",      dest="year",      help="which year",                     required=True)
    parser.add_argument("--inputDir",  dest="inputDir",  help="input directory",      type=str, required=True)
    parser.add_argument("--outputDir", dest="outputDir", help="output directory",     type=str, required=True)
    parser.add_argument("--model",     dest="model",     help="Which model to use",   type=str, required=True) 
    parser.add_argument("--channel",   dest="channel",   help="Which channel to use", type=str, required=True) 
    args = parser.parse_args()
    
    # Create the output directory if it already does not exist
    if not os.path.exists(args.outputDir):
        os.makedirs(args.outputDir)
    
    rocPlotter = RocAndRoll(args.approved, args.inputDir, args.outputDir, args.year, args.channel, "TT", args.model)

    rocPlotter.drawPerNjetROC(1)
    rocPlotter.drawPerMassROC(1)

    rocPlotter.drawPerNjetROC(2)
    rocPlotter.drawPerMassROC(2)
