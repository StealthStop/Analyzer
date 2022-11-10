import ROOT, os, sys, argparse

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat("")
ROOT.gStyle.SetLineWidth(2)
ROOT.gStyle.SetFrameLineWidth(1)
ROOT.gStyle.SetEndErrorSize(0)
ROOT.gStyle.SetPaintTextFormat("3.3f")
ROOT.TH1.SetDefaultSumw2()
ROOT.TH2.SetDefaultSumw2()

usage = "usage: %prog [options]"
parser = argparse.ArgumentParser(usage)
parser.add_argument("--tag"     , dest="tag"     , help="Unique tag for output"      , type=str, required=True)
parser.add_argument("--year"    , dest="year"    , help="Which year of run II"       , type=str, required=True)
parser.add_argument("--inputDir", dest="inputDir", help="Dir in base path with files", type=str, required=True)
parser.add_argument("--cr"      , dest="cr"      , help="Which control region"       , type=str, required=True)

arg = parser.parse_args()

def setMinimumErrors(dataMcRatio):

    xbins = dataMcRatio.GetXaxis().GetNbins()+1
    ybins = dataMcRatio.GetYaxis().GetNbins()+1

    for xbin in xrange(xbins):
        for ybin in xrange(ybins):

            content = dataMcRatio.GetBinContent(xbin,ybin)
            error = dataMcRatio.GetBinError(xbin,ybin)

            if error < 1e-10:
                
                newerror = content*0.03
                dataMcRatio.SetBinError(xbin,ybin,newerror)

def make1DRatioPlot(dataNum, dataDen, mcNum, mcDen, cr, goodName, outputFile):

    theName = dataNum.GetName().replace("num","ratio")

    XTitle = ""
    if "Pt" in dataNum.GetName(): XTitle = "Top Candidate p_{T} [GeV]"
    else: XTitle = "N_{jets}"

    XMin = 0;    XMax = 1; RatioXMin = 0; RatioXMax = 1
    YMin = 0.30; YMax = 1; RatioYMin = 0; RatioYMax = 0.30
    PadFactor = (YMax-YMin) / (RatioYMax-RatioYMin)

    XTitleSize = 0.05;  XLabelSize = 0.05;  XTitleOffset = 4.0 
    YTitleSize = 0.05;  YLabelSize = 0.05;  YTitleOffset = 1.2

    dataRatio = dataNum.Clone(); mcRatio = mcNum.Clone(); dataMcRatio = dataNum.Clone()
    dataRatio.Divide(dataNum, dataDen, 1, 1, "B"); mcRatio.Divide(mcNum, mcDen, 1, 1, "B")
    dataMcRatio.Divide(dataRatio, mcRatio)
    dataMcRatio.SetTitle("")

    if "ratioR" in theName:
        dataRatio.GetXaxis().SetRangeUser(0,1000); mcRatio.GetXaxis().SetRangeUser(0,1000); dataMcRatio.GetXaxis().SetRangeUser(0,1000)
    else:
        dataRatio.GetXaxis().SetRangeUser(400,1000); mcRatio.GetXaxis().SetRangeUser(400,1000); dataMcRatio.GetXaxis().SetRangeUser(400,1000)

    if "nJets" in theName:
        dataRatio.GetXaxis().SetRangeUser(6.5,10.5); mcRatio.GetXaxis().SetRangeUser(6.5,10.5); dataMcRatio.GetXaxis().SetRangeUser(6.5,10.5)
        dataMcRatio.GetXaxis().SetBinLabel(1, "7")
        dataMcRatio.GetXaxis().SetBinLabel(2, "8")
        dataMcRatio.GetXaxis().SetBinLabel(3, "#geq 9")

    dataRatio.SetLineColor(ROOT.kBlack); mcRatio.SetLineColor(ROOT.kRed)
    dataRatio.SetMarkerColor(ROOT.kBlack); mcRatio.SetMarkerColor(ROOT.kRed)
    dataRatio.SetLineWidth(3); mcRatio.SetLineWidth(3)
    dataRatio.SetMarkerSize(3.0); mcRatio.SetMarkerSize(3.0)
    dataRatio.SetMarkerStyle(20); mcRatio.SetMarkerStyle(20)
    if cr == "TTbar":
        dataRatio.GetYaxis().SetRangeUser(-0.05,1.1)
    else:
        dataRatio.GetYaxis().SetRangeUser(-0.02,0.2)
           
    dataRatio.SetTitle(""); mcRatio.SetTitle("")

    if cr == "TTbar":
        dataRatio.GetYaxis().SetTitle("Top Tagger Efficiency")
    else:
        dataRatio.GetYaxis().SetTitle("Top Tagger Mistag Rate")
       
    dataRatio.GetYaxis().SetTitleSize(YTitleSize); dataRatio.GetXaxis().SetTitleSize(XTitleSize)
    dataRatio.GetYaxis().SetLabelSize(YLabelSize); dataRatio.GetXaxis().SetLabelSize(0)
    dataRatio.GetYaxis().SetTitleOffset(YTitleOffset); dataRatio.GetXaxis().SetTitleOffset(XTitleOffset)

    aCanvas = ROOT.TCanvas("c_%s"%(theName), "c_%s"%(theName), 1400, 1400)
    ROOT.gPad.Clear()
    aCanvas.Divide(1,2)

    aCanvas.cd(1); ROOT.gPad.SetPad(XMin, YMin, XMax, YMax)
    ROOT.gPad.SetGridy()
    ROOT.gPad.SetGridx()

    TopMargin    = 0.1
    LeftMargin   = 0.12
    BottomMargin = 0.01
    RightMargin  = 0.05

    ROOT.gPad.SetTopMargin(TopMargin)
    ROOT.gPad.SetLeftMargin(LeftMargin)
    ROOT.gPad.SetBottomMargin(BottomMargin)
    ROOT.gPad.SetRightMargin(RightMargin)
    ROOT.gPad.SetTicks()

    dataRatio.Draw("E1 P"); mcRatio.Draw("E1 P SAME")

    iamLegend = None
    if cr == "TTbar":
        iamLegend = ROOT.TLegend(0.58, 0.70, 0.9, 0.85, "", "trNDC")
    else:
        iamLegend = ROOT.TLegend(0.58, 0.70, 0.9, 0.85, "", "trNDC")

    iamLegend.SetTextSize(0.035)
    if "TTbar" in cr:
        iamLegend.AddEntry(dataRatio, "SingleMuon Data", "E1 P")
        iamLegend.AddEntry(mcRatio, "TT MC", "E1 P")
    elif "QCD" in cr:
        iamLegend.AddEntry(dataRatio, "JetHT Data", "E1 P")
        iamLegend.AddEntry(mcRatio, "QCD multijet MC", "E1 P")

    iamLegend.Draw("SAME")

    mark = ROOT.TLatex()
    mark.SetNDC(True)

    mark.SetTextAlign(11)
    mark.SetTextSize(0.09)
    mark.SetTextFont(61)
    mark.DrawLatex(LeftMargin, 1 - (TopMargin - 0.015), "CMS")

    mark.SetTextFont(52)
    mark.SetTextSize(0.060)

    mark.DrawLatex(LeftMargin + 0.14, 1 - (TopMargin - 0.017), "work in progress")

    mark.SetTextSize(0.050)
    mark.SetTextFont(42)
    mark.SetTextAlign(31)
    mark.DrawLatex(1 - RightMargin, 1 - (TopMargin - 0.017), "%s (13 TeV)"%(year))

    aCanvas.cd(2); ROOT.gPad.SetPad(RatioXMin, RatioYMin, RatioXMax, RatioYMax)
    ROOT.gPad.SetTopMargin(0.05)
    ROOT.gPad.SetBottomMargin(0.35)
    ROOT.gPad.SetRightMargin(RightMargin)
    ROOT.gPad.SetLeftMargin(LeftMargin)
    ROOT.gPad.SetTicks()
    ROOT.gPad.SetGridy()

    if "topPt" in goodName and "Njets" in goodName:
        if "QCD" in goodName:
            goodName = goodName.replace("Binned", "MisTagSF_vs").replace("_QCDCR", "")
        else:
            goodName = goodName.replace("Binned", "TagRateSF_vs").replace("_TTbarCR", "")

        goodName = goodName.replace("numR", "Resolved").replace("numM", "Merged")
        outputFile.cd()
        dataMcRatio.SetName(goodName)
        dataMcRatio.Write(goodName)

    dataMcRatio.SetMinimum(0.4)
    dataMcRatio.SetMaximum(1.6)
    dataMcRatio.GetYaxis().SetNdivisions(-304)
    dataMcRatio.SetTitle("")
    dataMcRatio.SetLineWidth(3)
    dataMcRatio.SetMarkerSize(3.0)
    dataMcRatio.SetMarkerStyle(20)
    dataMcRatio.SetMarkerColor(dataMcRatio.GetLineColor())
    dataMcRatio.GetYaxis().SetTitle("Data / MC")
    dataMcRatio.GetXaxis().SetTitle(XTitle)
    dataMcRatio.GetXaxis().SetTitleSize(PadFactor*XTitleSize); dataMcRatio.GetYaxis().SetTitleSize(0.8*PadFactor*YTitleSize)
    if "nJets" in theName:
        dataMcRatio.GetXaxis().SetLabelSize(1.5*PadFactor*XLabelSize); dataMcRatio.GetYaxis().SetLabelSize(PadFactor*YLabelSize)
    else:
        dataMcRatio.GetXaxis().SetLabelSize(PadFactor*XLabelSize); dataMcRatio.GetYaxis().SetLabelSize(PadFactor*YLabelSize)

    dataMcRatio.GetYaxis().SetTitleOffset(1.2*YTitleOffset/PadFactor); dataMcRatio.GetXaxis().SetTitleOffset(0.65*XTitleOffset/PadFactor)

    dataMcRatio.Draw("E1 P")

    aCanvas.SaveAs("Studies/TopTaggerSF/%s/%s_%s.pdf"%(tag,year,theName))

def make2DRatioPlot(dataNum, dataDen, mcNum, mcDen, aName, outputFile):

    XTitle = ""
    if "Pt" in dataNum.GetName(): XTitle = "Top Candidate p_{T} [GeV]"

    XMin = 0;    XMax = 1; RatioXMin = 0; RatioXMax = 1
    YMin = 0.20; YMax = 1; RatioYMin = 0; RatioYMax = 0.20
    PadFactor = (YMax-YMin) / (RatioYMax-RatioYMin)

    XTitleSize = 0.045;  XLabelSize = 0.036;  XTitleOffset = 1.0 
    YTitleSize = 0.045;  YLabelSize = 0.036;  YTitleOffset = 1.3 

    dataRatio = dataNum.Clone(); mcRatio = mcNum.Clone(); dataMcRatio = dataNum.Clone()
    dataRatio.Divide(dataNum, dataDen, 1, 1, "B"); mcRatio.Divide(mcNum, mcDen, 1, 1, "B")
    dataMcRatio.Divide(dataRatio, mcRatio)
    dataMcRatio.SetTitle("")

    dataRatio.GetZaxis().SetRangeUser(0.4,1.6)
    dataRatio.SetTitle(""); mcRatio.SetTitle("")
    dataRatio.GetYaxis().SetTitle("N_{jets}"); dataRatio.GetXaxis().SetTitle(XTitle)
    dataRatio.GetYaxis().SetTitleSize(YTitleSize); dataRatio.GetXaxis().SetTitleSize(XTitleSize)
    dataRatio.GetYaxis().SetLabelSize(YLabelSize); dataRatio.GetXaxis().SetLabelSize(XLabelSize)
    dataRatio.GetYaxis().SetTitleOffset(YTitleOffset); dataRatio.GetXaxis().SetTitleOffset(XTitleOffset)
    dataRatio.SetContour(255)

    mcRatio.GetZaxis().SetRangeUser(0.40,1.6)
    mcRatio.SetTitle(""); mcRatio.SetTitle("")
    mcRatio.GetYaxis().SetTitle("N_{jets}"); mcRatio.GetXaxis().SetTitle(XTitle)
    mcRatio.GetYaxis().SetTitleSize(YTitleSize); mcRatio.GetXaxis().SetTitleSize(XTitleSize)
    mcRatio.GetYaxis().SetLabelSize(YLabelSize); mcRatio.GetXaxis().SetLabelSize(XLabelSize)
    mcRatio.GetYaxis().SetTitleOffset(YTitleOffset); mcRatio.GetXaxis().SetTitleOffset(XTitleOffset)
    mcRatio.SetContour(255)

    dataMcRatio.GetZaxis().SetRangeUser(0.4,1.6)
    dataMcRatio.SetTitle(""); dataMcRatio.SetTitle("")
    dataMcRatio.GetYaxis().SetTitle("N_{jets}"); dataMcRatio.GetXaxis().SetTitle(XTitle)
    dataMcRatio.GetYaxis().SetTitleSize(YTitleSize); dataMcRatio.GetXaxis().SetTitleSize(XTitleSize)
    dataMcRatio.GetYaxis().SetLabelSize(1.5*YLabelSize); dataMcRatio.GetXaxis().SetLabelSize(XLabelSize)
    dataMcRatio.GetYaxis().SetTitleOffset(YTitleOffset); dataMcRatio.GetXaxis().SetTitleOffset(XTitleOffset)
    dataMcRatio.SetContour(255)
    dataMcRatio.GetYaxis().SetNdivisions(-4)

    dataMcRatio.GetYaxis().SetBinLabel(1, "7")
    dataMcRatio.GetYaxis().SetBinLabel(2, "8")
    dataMcRatio.GetYaxis().SetBinLabel(3, "9")
    dataMcRatio.GetYaxis().SetBinLabel(4, "#geq 10")

    setMinimumErrors(dataMcRatio)

    theName = dataNum.GetName().replace("num","ratio")

    aCanvas = ROOT.TCanvas("c_data_%s"%(theName), "c_data_%s"%(theName), 1400, 1400)
    ROOT.gPad.Clear()
    aCanvas.cd()

    TopMargin    = 0.07
    LeftMargin   = 0.12
    BottomMargin = 0.11
    RightMargin  = 0.12

    ROOT.gPad.SetTopMargin(TopMargin)
    ROOT.gPad.SetLeftMargin(LeftMargin)
    ROOT.gPad.SetBottomMargin(BottomMargin)
    ROOT.gPad.SetRightMargin(RightMargin)

    ROOT.gPad.SetTicks()

    dataRatio.Draw("COLZ E TEXT")
    mark = ROOT.TLatex()
    mark.SetNDC(True)

    mark.SetTextAlign(11)
    mark.SetTextSize(0.065)
    mark.SetTextFont(61)
    mark.DrawLatex(LeftMargin, 1 - (TopMargin - 0.015), "CMS")

    mark.SetTextFont(52)
    mark.SetTextSize(0.040)

    mark.DrawLatex(LeftMargin + 0.22, 1 - (TopMargin - 0.017), "work in progress")

    mark.SetTextFont(42)
    mark.SetTextAlign(31)
    mark.DrawLatex(1 - RightMargin, 1 - (TopMargin - 0.017), "%s (13 TeV)"%(year))

    aCanvas.SaveAs("Studies/TopTaggerSF/%s/%s_%s.pdf"%(tag,year,"data_"+theName))

    aCanvas = ROOT.TCanvas("c_mc_%s"%(theName), "c_mc_%s"%(theName), 1500, 1200)
    ROOT.gPad.Clear()
    aCanvas.cd()
    ROOT.gPad.SetTopMargin(TopMargin)
    ROOT.gPad.SetLeftMargin(LeftMargin)
    ROOT.gPad.SetBottomMargin(BottomMargin)
    ROOT.gPad.SetRightMargin(RightMargin)
    ROOT.gPad.SetTicks()

    mcRatio.Draw("COLZ E TEXT")
    mark = ROOT.TLatex()
    mark.SetNDC(True)

    mark.SetTextAlign(11)
    mark.SetTextSize(0.050)
    mark.SetTextFont(61)
    mark.DrawLatex(LeftMargin, 1 - (TopMargin - 0.015), "CMS")

    mark.SetTextFont(52)
    mark.SetTextSize(0.040)

    mark.DrawLatex(LeftMargin + 0.14, 1 - (TopMargin - 0.017), "work in Progress")

    mark.SetTextFont(42)
    mark.SetTextAlign(31)
    mark.DrawLatex(1 - RightMargin, 1 - (TopMargin - 0.017), "%s (13 TeV)"%(year))

    aCanvas.SaveAs("Studies/TopTaggerSF/%s/%s_%s.pdf"%(tag,year,"mc_"+theName))

    aCanvas = ROOT.TCanvas("c_datamc_%s"%(theName), "c_datamc_%s"%(theName), 1400, 1400)
    ROOT.gPad.Clear()
    aCanvas.cd()
    ROOT.gPad.SetTopMargin(TopMargin)
    ROOT.gPad.SetLeftMargin(LeftMargin)
    ROOT.gPad.SetBottomMargin(BottomMargin)
    ROOT.gPad.SetRightMargin(RightMargin)
    ROOT.gPad.SetTicks()

    dataMcRatio.Draw("COLZ E TEXT")
    mark = ROOT.TLatex()
    mark.SetNDC(True)

    mark.SetTextAlign(11)
    mark.SetTextSize(0.059)
    mark.SetTextFont(61)
    mark.DrawLatex(LeftMargin, 1 - (TopMargin - 0.015), "CMS")

    mark.SetTextFont(52)
    mark.SetTextSize(0.039)

    mark.DrawLatex(LeftMargin + 0.13, 1 - (TopMargin - 0.017), "work in progress")

    mark.SetTextSize(0.035)
    mark.SetTextFont(42)
    mark.SetTextAlign(31)
    mark.DrawLatex(1 - RightMargin, 1 - (TopMargin - 0.017), "%s (13 TeV)"%(year))

    aCanvas.SaveAs("Studies/TopTaggerSF/%s/%s_%s.pdf"%(tag, year, "datamc_"+theName))

if __name__ == '__main__':

    varTags      = [ "topPt", "nJets" ]
    topTags      = [ "R", "M" ]
    nJetCutTags  = [ "7", "8", "9incl", "" ]

    inputDir = arg.inputDir
    tag      = arg.tag
    year     = arg.year
    cr       = arg.cr

    fullPath = inputDir

    # For pure top, we include semileptonic ttbar, single top, and tt + W/Z/H
    # Non-pure top is everything else
    mcFilePureTop    = ROOT.TFile.Open(fullPath + "/" + year + "_PureTop.root")
    mcFileNotPureTop = ROOT.TFile.Open(fullPath + "/" + year + "_NotPureTop.root")
    dataFile         = 0
    if   cr == "TTbar":
        dataFile  = ROOT.TFile.Open(fullPath + "/" + year + "_Data_SingleMuon.root")
    elif cr == "QCD":
        dataFile = ROOT.TFile.Open(fullPath + "/" + year + "_Data_JetHT.root")

    outputPath = "./Studies/TopTaggerSF/%s/%s"%(tag,arg.year)
    if not os.path.exists(outputPath):
        os.makedirs(outputPath)

    outputFileName = "%s_TopTaggerSF.root"%(year)
    outputFile = ROOT.TFile.Open("%s/%s"%(outputPath,outputFileName), "UPDATE")

    for topTag in topTags:
        denomTag2 = "Binned_topPt_vs_nJets" + "_" + cr + "CR_den%s"%(topTag)
        numerTag2 = "Binned_topPt_vs_nJets" + "_" + cr + "CR_num%s"%(topTag)

        goodName = "TopTagEff" + "_" + year + "_num"

        dataFile.cd()
        h2DataNum = dataFile.Get(numerTag2)
        h2DataDen = dataFile.Get(denomTag2)

        mcFilePureTop.cd()
        h2McTopNum = mcFilePureTop.Get(numerTag2)
        h2McTopDen = mcFilePureTop.Get(denomTag2)

        mcFileNotPureTop.cd()
        h2McNotPureTopNum = mcFileNotPureTop.Get(numerTag2)
        h2McNotPureTopDen = mcFileNotPureTop.Get(denomTag2)

        mcNum = 0
        mcDen = 0
        dataNum = 0
        dataDen = 0
        if cr == "TTbar":
            mcNum = h2McTopNum
            mcDen = h2McTopDen
            dataNum = h2DataNum
            dataNum.Add(h2McNotPureTopNum, -1.0)
            
            dataDen = h2DataDen
            dataDen.Add(h2McNotPureTopDen, -1.0)
        elif cr == "QCD":
            mcNum = h2McNotPureTopNum
            mcDen = h2McNotPureTopDen
            dataNum = h2DataNum
            dataNum.Add(h2McTopNum, -1.0)
            
            dataDen = h2DataDen
            dataDen.Add(h2McTopDen, -1.0)

        make2DRatioPlot(dataNum, dataDen, mcNum, mcDen, goodName, outputFile)

        for varTag in varTags:
            for nJetCutTag in nJetCutTags:
                nJetStr = ""
                if nJetCutTag != "":
                    nJetStr = "_Njets" + nJetCutTag

                denomTag1 = "Binned_%s"%(varTag) + "_" + cr + "CR_den%s"%(topTag) + nJetStr
                numerTag1 = "Binned_%s"%(varTag) + "_" + cr + "CR_num%s"%(topTag) + nJetStr

                dataFile.cd()
                hDataNum = dataFile.Get(numerTag1)
                hDataDen = dataFile.Get(denomTag1)

                mcFilePureTop.cd()
                hMcTopNum = mcFilePureTop.Get(numerTag1)
                hMcTopDen = mcFilePureTop.Get(denomTag1)

                mcFileNotPureTop.cd()
                hMcNotPureTopNum = mcFileNotPureTop.Get(numerTag1)
                hMcNotPureTopDen = mcFileNotPureTop.Get(denomTag1)

                mcNum = 0
                mcDen = 0
                dataNum = 0
                dataDen = 0
                if cr == "TTbar":
                    mcNum = hMcTopNum
                    mcDen = hMcTopDen
                    dataNum = hDataNum
                    dataNum.Add(hMcNotPureTopNum, -1.0)
                    
                    dataDen = hDataDen
                    dataDen.Add(hMcNotPureTopDen, -1.0)
                elif cr == "QCD":
                    mcNum = hMcNotPureTopNum
                    mcDen = hMcNotPureTopDen
                    dataNum = hDataNum
                    dataNum.Add(hMcTopNum, -1.0)
                    
                    dataDen = hDataDen
                    dataDen.Add(hMcTopDen, -1.0)

                make1DRatioPlot(dataNum, dataDen, mcNum, mcDen, cr, numerTag1, outputFile) 
