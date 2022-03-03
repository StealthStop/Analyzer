import ROOT, os, sys, argparse

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat("")
ROOT.gStyle.SetLineWidth(2)
ROOT.gStyle.SetFrameLineWidth(1)
ROOT.gStyle.SetPaintTextFormat("3.3f")
ROOT.TH1.SetDefaultSumw2()
ROOT.TH2.SetDefaultSumw2()

usage = "usage: %prog [options]"
parser = argparse.ArgumentParser(usage)
parser.add_argument("--tag"      , dest="tag"      , help="Unique tag for output"         , type=str, required=True)
parser.add_argument("--year"     , dest="year"     , help="Which year of run II"          , type=str, required=True)
parser.add_argument("--basePath" , dest="basePath" , help="Base path to input"            , type=str, required=True)
parser.add_argument("--inputDir" , dest="inputDir" , help="Dir in base path with files"   , type=str, required=True)
parser.add_argument("--dataset"  , dest="dataset"  , help="The muon or electron dataset"  , type=str, required=True)

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

def make1DRatioPlot(dataNum, dataDen, mcNum, mcDen):

    theName = dataNum.GetName().replace("num","ratio")[2:]

    XTitle = ""
    if "Pt" in dataNum.GetName(): XTitle = "p_{T} [GeV]"
    else: XTitle = "#eta"

    XMin = 0;    XMax = 1; RatioXMin = 0; RatioXMax = 1
    YMin = 0.30; YMax = 1; RatioYMin = 0; RatioYMax = 0.30
    PadFactor = (YMax-YMin) / (RatioYMax-RatioYMin)

    XTitleSize = 0.045;  XLabelSize = 0.036;  XTitleOffset = 4.0 
    YTitleSize = 0.045;  YLabelSize = 0.036;  YTitleOffset = 0.9

    dataRatio = dataNum.Clone(); mcRatio = mcNum.Clone(); dataMcRatio = dataNum.Clone()
    dataRatio.Divide(dataNum, dataDen, 1, 1, "B"); mcRatio.Divide(mcNum, mcDen, 1, 1, "B")
    dataMcRatio.Divide(dataRatio, mcRatio)
    dataMcRatio.SetTitle("")

    dataRatio.GetXaxis().SetRangeUser(30,150); mcRatio.GetXaxis().SetRangeUser(30,150); dataMcRatio.GetXaxis().SetRangeUser(30,150)

    dataRatio.SetLineColor(ROOT.kBlack); mcRatio.SetLineColor(ROOT.kRed)
    dataRatio.SetMarkerColor(ROOT.kBlack); mcRatio.SetMarkerColor(ROOT.kRed)
    dataRatio.SetLineWidth(3); mcRatio.SetLineWidth(3)
    dataRatio.SetMarkerSize(1.7); mcRatio.SetMarkerSize(1.7)
    dataRatio.SetMarkerStyle(20); mcRatio.SetMarkerStyle(20)
    dataRatio.GetYaxis().SetRangeUser(0.31,1.05)
    dataRatio.SetTitle(""); mcRatio.SetTitle("")
    dataRatio.GetYaxis().SetTitle("L1 + HLT Efficiency")
    dataRatio.GetYaxis().SetTitleSize(YTitleSize); dataRatio.GetXaxis().SetTitleSize(XTitleSize)
    dataRatio.GetYaxis().SetLabelSize(YLabelSize); dataRatio.GetXaxis().SetLabelSize(0)
    dataRatio.GetYaxis().SetTitleOffset(YTitleOffset); dataRatio.GetXaxis().SetTitleOffset(XTitleOffset)

    aCanvas = ROOT.TCanvas("c_%s"%(theName), "c_%s"%(theName), 1200, 1200)
    ROOT.gPad.Clear()
    aCanvas.Divide(1,2)

    aCanvas.cd(1); ROOT.gPad.SetPad(XMin, YMin, XMax, YMax)
    ROOT.gPad.SetGridy()
    ROOT.gPad.SetTopMargin(0.01)
    ROOT.gPad.SetLeftMargin(0.09)
    ROOT.gPad.SetBottomMargin(0.01)
    ROOT.gPad.SetRightMargin(0.025)
    ROOT.gPad.SetTicks()

    dataRatio.Draw("E1 P"); mcRatio.Draw("E1 P SAME")

    iamLegend = ROOT.TLegend(0.58, 0.20, 0.94, 0.42, "", "trNDC")
    iamLegend.SetTextSize(0.035)
    if dataset == "electron":
        iamLegend.AddEntry(dataRatio, "Data_SingleElectron", "E1 P")
    elif dataset == "muon":
        iamLegend.AddEntry(dataRatio, "Data_SingleMuon", "E1 P")
    iamLegend.AddEntry(mcRatio, "MC (TT)", "E1 P")

    iamLegend.Draw("SAME")

    aCanvas.cd(2); ROOT.gPad.SetPad(RatioXMin, RatioYMin, RatioXMax, RatioYMax)
    ROOT.gPad.SetTopMargin(0.05)
    ROOT.gPad.SetBottomMargin(0.26)
    ROOT.gPad.SetRightMargin(0.025)
    ROOT.gPad.SetLeftMargin(0.09)
    ROOT.gPad.SetTicks()
    ROOT.gPad.SetGridy()

    dataMcRatio.SetMinimum(0.8)
    dataMcRatio.SetMaximum(1.2)
    dataMcRatio.GetYaxis().SetNdivisions(-304)
    dataMcRatio.SetTitle("")
    dataMcRatio.SetLineWidth(3)
    dataMcRatio.SetMarkerSize(1.7)
    dataMcRatio.SetMarkerStyle(20)
    dataMcRatio.SetMarkerColor(dataMcRatio.GetLineColor())
    dataMcRatio.GetYaxis().SetTitle("Data / MC")
    dataMcRatio.GetXaxis().SetTitle(XTitle)
    dataMcRatio.GetXaxis().SetTitleSize(PadFactor*XTitleSize); dataMcRatio.GetYaxis().SetTitleSize(0.8*PadFactor*YTitleSize)
    dataMcRatio.GetXaxis().SetLabelSize(PadFactor*XLabelSize); dataMcRatio.GetYaxis().SetLabelSize(PadFactor*YLabelSize)
    dataMcRatio.GetYaxis().SetTitleOffset(1.3*YTitleOffset/PadFactor); dataMcRatio.GetXaxis().SetTitleOffset(0.65*XTitleOffset/PadFactor)

    dataMcRatio.Draw("E1 P")

    aCanvas.SaveAs("plots/%s/%s/1D/%s.pdf"%(tag,year,theName))

def make2DRatioPlot(dataNum, dataDen, mcNum, mcDen, aName, outputFile):

    XTitle = ""
    if "Pt" in dataNum.GetName(): XTitle = "p_{T} [GeV]"
    else: XTitle = "#eta"

    XMin = 0;    XMax = 1; RatioXMin = 0; RatioXMax = 1
    YMin = 0.20; YMax = 1; RatioYMin = 0; RatioYMax = 0.20
    PadFactor = (YMax-YMin) / (RatioYMax-RatioYMin)

    XTitleSize = 0.045;  XLabelSize = 0.036;  XTitleOffset = 1.0 
    YTitleSize = 0.045;  YLabelSize = 0.036;  YTitleOffset = 0.9

    dataRatio = dataNum.Clone(); mcRatio = mcNum.Clone(); dataMcRatio = dataNum.Clone()
    dataRatio.Divide(dataNum, dataDen, 1, 1, "B"); mcRatio.Divide(mcNum, mcDen, 1, 1, "B")
    dataMcRatio.Divide(dataRatio, mcRatio)
    dataMcRatio.SetTitle("")

    dataRatio.GetZaxis().SetRangeUser(0.75,1.1)
    dataRatio.SetTitle(""); mcRatio.SetTitle("")
    dataRatio.GetYaxis().SetTitle("#eta"); dataRatio.GetXaxis().SetTitle("p_{T} [GeV]")
    dataRatio.GetYaxis().SetTitleSize(YTitleSize); dataRatio.GetXaxis().SetTitleSize(XTitleSize)
    dataRatio.GetYaxis().SetLabelSize(YLabelSize); dataRatio.GetXaxis().SetLabelSize(XLabelSize)
    dataRatio.GetYaxis().SetTitleOffset(YTitleOffset); dataRatio.GetXaxis().SetTitleOffset(XTitleOffset)
    dataRatio.SetContour(255)

    mcRatio.GetZaxis().SetRangeUser(0.75,1.1)
    mcRatio.SetTitle(""); mcRatio.SetTitle("")
    mcRatio.GetYaxis().SetTitle("#eta"); mcRatio.GetXaxis().SetTitle("p_{T} [GeV]")
    mcRatio.GetYaxis().SetTitleSize(YTitleSize); mcRatio.GetXaxis().SetTitleSize(XTitleSize)
    mcRatio.GetYaxis().SetLabelSize(YLabelSize); mcRatio.GetXaxis().SetLabelSize(XLabelSize)
    mcRatio.GetYaxis().SetTitleOffset(YTitleOffset); mcRatio.GetXaxis().SetTitleOffset(XTitleOffset)
    mcRatio.SetContour(255)

    dataMcRatio.GetZaxis().SetRangeUser(0.75,1.1)
    dataMcRatio.SetTitle(""); dataMcRatio.SetTitle("")
    dataMcRatio.GetYaxis().SetTitle("#eta"); dataMcRatio.GetXaxis().SetTitle("p_{T} [GeV]")
    dataMcRatio.GetYaxis().SetTitleSize(YTitleSize); dataMcRatio.GetXaxis().SetTitleSize(XTitleSize)
    dataMcRatio.GetYaxis().SetLabelSize(YLabelSize); dataMcRatio.GetXaxis().SetLabelSize(XLabelSize)
    dataMcRatio.GetYaxis().SetTitleOffset(YTitleOffset); dataMcRatio.GetXaxis().SetTitleOffset(XTitleOffset)
    dataMcRatio.SetContour(255)

    setMinimumErrors(dataMcRatio)

    if "5jCut" in aName and "trig" in aName:
        outputFile.cd()
        dataMcRatio.SetName(aName)
        dataMcRatio.Write(aName)

    theName = dataNum.GetName().replace("num","ratio")[3:]

    aCanvas = ROOT.TCanvas("c_data_%s"%(theName), "c_data_%s"%(theName), 1500, 1200)
    ROOT.gPad.Clear()
    aCanvas.cd()
    ROOT.gPad.SetTopMargin(0.02)
    ROOT.gPad.SetLeftMargin(0.09)
    ROOT.gPad.SetBottomMargin(0.11)
    ROOT.gPad.SetRightMargin(0.11)

    ROOT.gPad.SetTicks()

    dataRatio.Draw("COLZ E TEXT")
    aCanvas.SaveAs("plots/%s/%s/2D/%s.pdf"%(tag,year,"data_"+theName))

    aCanvas = ROOT.TCanvas("c_mc_%s"%(theName), "c_mc_%s"%(theName), 1500, 1200)
    ROOT.gPad.Clear()
    aCanvas.cd()
    ROOT.gPad.SetTopMargin(0.02)
    ROOT.gPad.SetLeftMargin(0.09)
    ROOT.gPad.SetBottomMargin(0.11)
    ROOT.gPad.SetRightMargin(0.11)
    ROOT.gPad.SetTicks()

    mcRatio.Draw("COLZ E TEXT")
    aCanvas.SaveAs("plots/%s/%s/2D/%s.pdf"%(tag,year,"mc_"+theName))

    aCanvas = ROOT.TCanvas("c_datamc_%s"%(theName), "c_datamc_%s"%(theName), 1500, 1200)
    ROOT.gPad.Clear()
    aCanvas.cd()
    ROOT.gPad.SetTopMargin(0.02)
    ROOT.gPad.SetLeftMargin(0.09)
    ROOT.gPad.SetBottomMargin(0.11)
    ROOT.gPad.SetRightMargin(0.11)
    ROOT.gPad.SetTicks()

    dataMcRatio.Draw("COLZ E TEXT")
    aCanvas.SaveAs("plots/%s/%s/2D/%s.pdf"%(tag, year, "datamc_"+theName))

if __name__ == '__main__':

    effTags      = [ "den", "num" ] 
    lepTags      = [ "el", "mu" ]
    binTags      = [ "Eta", "Pt" ]
    ptTags       = [ "pt40" ]
    trigTags     = [ "trig" ] 
    nJetCutTags  = [ "5", "6" ]

    basePath = arg.basePath
    inputDir = arg.inputDir
    tag = arg.tag
    year = arg.year
    dataset = arg.dataset

    fullPath = basePath + "/" + inputDir

    mcFile = ROOT.TFile.Open(fullPath + "/" + year + "_TT.root")

    dataFile = 0
    if dataset == "electron": dataFile = ROOT.TFile.Open(fullPath + "/" + year + "_Data_SingleElectron.root")
    elif dataset == "muon": dataFile = ROOT.TFile.Open(fullPath + "/" + year + "_Data_SingleMuon.root")


    outputPath = "./plots/%s/%s"%(tag,arg.year)
    outputFileName = "%s_TrigEff.root"%(year)

    if not os.path.exists(outputPath+"/1D"): os.makedirs(outputPath+"/1D")
    if not os.path.exists(outputPath+"/2D"): os.makedirs(outputPath+"/2D")

    outputFile = ROOT.TFile.Open("%s/%s"%(outputPath,outputFileName), "UPDATE")

    for lepTag in lepTags:
        if (lepTag == "el" and dataset == "electron") or (lepTag == "mu" and dataset == "muon"): continue
        for ptTag in ptTags:
            for trigTag in trigTags:
                for nJetCutTag in nJetCutTags:
                    denomTag2 = "h2_trig_den_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"jCut_wLepPtLepEtaBin"
                    numerTag2 = "h2_trig_num_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"jCut_wLepPtLepEtaBin"

                    goodName = "TrigEff" + "_" + year + "_num_" + lepTag + "_" + ptTag + "_" + trigTag + "_" + nJetCutTag + "jCut_htCut_DeepCSV"

                    dataFile.cd()
                    h2DataNum = dataFile.Get(numerTag2); h2DataDen = dataFile.Get(denomTag2)

                    mcFile.cd()
                    h2McNum = mcFile.Get(numerTag2); h2McDen = mcFile.Get(denomTag2)

                    if h2DataNum == None or h2DataDen == None or h2McNum == None or h2McDen == None: continue

                    make2DRatioPlot(h2DataNum, h2DataDen, h2McNum, h2McDen, goodName, outputFile)

                    for binTag in binTags:

                        denomTag1 = "h_trig_den_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"jCut_wLep%sBin"%(binTag)
                        numerTag1 = "h_trig_num_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"jCut_wLep%sBin"%(binTag)

                        dataFile.cd()
                        hDataNum = dataFile.Get(numerTag1); hDataDen = dataFile.Get(denomTag1)
                        mcFile.cd()
                        hMcNum = mcFile.Get(numerTag1); hMcDen = mcFile.Get(denomTag1)

                        if hDataNum == None or hDataDen == None or hMcNum == None or hMcDen == None: continue

                        make1DRatioPlot(hDataNum, hDataDen, hMcNum, hMcDen) 
