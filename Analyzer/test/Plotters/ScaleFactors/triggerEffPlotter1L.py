import argparse, os, sys, math
from array import array
from ctypes import c_int

# --------------------------------------------
# --------------------------------------------
usage = "usage: %prog [options]"
parser = argparse.ArgumentParser(usage)
parser.add_argument("--basePath", dest="basePath", help="Base path to input"          , type=str, required=True)
parser.add_argument("--inputDir", dest="inputDir", help="Dir in base path with files" , type=str, required=True)
parser.add_argument("--dataset" , dest="dataset" , help="The muon or electron dataset", type=str, required=True)
parser.add_argument("--year"    , dest="year"    , help="Which year of run II"        , type=str, required=True)
parser.add_argument("--tag"     , dest="tag"     , help="Unique tag for output"       , type=str, required=True)
arg = parser.parse_args()

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat("")
ROOT.gStyle.SetLineWidth(2)
ROOT.gStyle.SetFrameLineWidth(1)
#ROOT.gStyle.SetEndErrorSize(0)
ROOT.gStyle.SetPaintTextFormat("3.3f")
ROOT.TH1.SetDefaultSumw2()
ROOT.TH2.SetDefaultSumw2()

# -----------------
# set minimum error
# -----------------
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

# ----------------
# Devin's function
# ----------------
def updateErrors(hist,den):

    xbins = hist.GetXaxis().GetNbins()+1
    ybins = hist.GetYaxis().GetNbins()+1

    for xbin in range(1,xbins):
        for ybin in range(1,ybins):

            content = hist.GetBinContent(xbin,ybin)
            error   = hist.GetBinError(xbin,ybin)

            # If eff = 1 (errors undefined), set conservative, symmetric errors equal to Clopper-Pearson lower interval distance
            if content == 1.0 and error < 1e-10:
                error = 1.0 - (ROOT.Math.normal_cdf(-1,1,0))**(1.0 / den.GetBinContent(xbin,ybin))
                hist.SetBinError(xbin,ybin,error)

            if error < 1e-10:
                print('WARNING: Error in bin ({},{}) is < 1e-10. Setting error to 3% of bin content.'.format(xbin,ybin))
                error = content*0.03
                hist.SetBinError(xbin,ybin,error)

# ------------
# Add CMS logo
# ------------
def addCMSlogo(canvas, year, TopMargin, LeftMargin, RightMargin, SF=1.0):

    canvas.cd()
    mark = ROOT.TLatex()
    mark.SetNDC(True)
    mark.SetTextAlign(11)
    mark.SetTextSize(0.048)
    mark.SetTextFont(61)
    mark.DrawLatex(LeftMargin, (1 - (TopMargin - 0.034)*SF), "CMS")
    mark.SetTextFont(52)
    mark.SetTextSize(0.038)
    mark.DrawLatex(LeftMargin + 0.12, (1 - (TopMargin - 0.034)*SF), "Work in Progress")
    mark.SetTextAlign(31)
    mark.DrawLatex(1 - RightMargin, (1 - (TopMargin - 0.034)*SF), "%s (13 TeV)"%(year))

# -------------------------------------------
# Make Efficiency plots
# they are 1D plots
#   -- denominator = preselections
#   -- numerator   = preselections + triggers
#   -- efficiency  = numerator / denominator
# -------------------------------------------
def make_Efficiency_Plots(dataNum, dataDen, mcNum, mcDen):

    XTitle = ""
    if "Pt" in dataNum.GetName(): XTitle = "p_{T} [GeV]"
    else: XTitle = "#eta"

    XMin = 0;    XMax = 1; RatioXMin = 0; RatioXMax = 1
    YMin = 0.30; YMax = 1; RatioYMin = 0; RatioYMax = 0.30
    PadFactor  = (YMax-YMin) / (RatioYMax-RatioYMin)
    XTitleSize = 0.05;  XLabelSize = 0.05;  XTitleOffset = 4.0 
    YTitleSize = 0.05;  YLabelSize = 0.05;  YTitleOffset = 0.9

    # get data efficiency
    #dataRatio = dataNum.Clone()
    dataNum.GetYaxis().SetTitle("L1+HLT Efficiency")
    dataNum.AddBinContent(dataNum.GetNbinsX(), dataNum.GetBinContent(dataNum.GetNbinsX()+1))
    dataDen.GetYaxis().SetTitle("L1+HLT Efficiency")
    dataDen.AddBinContent(dataDen.GetNbinsX(), dataDen.GetBinContent(dataDen.GetNbinsX()+1))
    dataRatio = ROOT.TEfficiency(dataNum, dataDen)
    #dataRatio.SetTitle("")
    #dataRatio.GetYaxis().SetTitle("L1+HLT Efficiency")
    #dataRatio.GetXaxis().SetLabelSize(0)
    #dataRatio.GetYaxis().SetTitleSize(0.15)
    #dataRatio.GetYaxis().SetLabelSize(0.15)
    dataRatio.SetLineWidth(3)
    dataRatio.SetLineColor(1)
    dataRatio.SetMarkerColor(1)
    dataRatio.SetMarkerSize(2.0)
    dataRatio.SetMarkerStyle(20)

    ##dataRatio.Divide(dataNum, dataDen, 1, 1, "B")
    #dataRatio.GetYaxis().SetRangeUser(0.4, 1.1)
    # get mc efficiency
    #mcRatio = mcNum.Clone()
    mcNum.GetYaxis().SetTitle("L1+HLT Efficiency")
    mcNum.AddBinContent(mcNum.GetNbinsX(), mcNum.GetBinContent(mcNum.GetNbinsX()+1))
    mcDen.GetYaxis().SetTitle("L1+HLT Efficiency")
    mcDen.AddBinContent(mcDen.GetNbinsX(), mcDen.GetBinContent(mcDen.GetNbinsX()+1))
    mcRatio = ROOT.TEfficiency(mcNum, mcDen)
    mcRatio.SetLineWidth(3)
    mcRatio.SetLineColor(2)
    mcRatio.SetMarkerColor(2)
    mcRatio.SetMarkerSize(2.0)
    mcRatio.SetMarkerStyle(20)
    ##mcRatio.Divide(mcNum, mcDen, 1, 1, "B")
    #mcRatio.GetYaxis().SetRangeUser(0.4, 1.1)
    # get errors
    #updateErrors(dataRatio,dataDen)
    #updateErrors(mcRatio,mcDen)
    # get scale factor
    dataRatioGraph = dataRatio.CreateGraph()
    mcRatioGraph = mcRatio.CreateGraph()
    dataMcRatio = dataRatioGraph.Clone()
    for i in range(1,dataMcRatio.GetN()+1):
        x_vals = dataRatioGraph.GetX()
        try:
            temp_err_high = math.sqrt((dataRatioGraph.GetErrorYhigh(i-1) / dataRatioGraph.Eval(x_vals[i-1]))**2 + (mcRatioGraph.GetErrorYhigh(i-1) / mcRatioGraph.Eval(x_vals[i-1]))**2)
            temp_err_low = math.sqrt((dataRatioGraph.GetErrorYlow(i-1) / dataRatioGraph.Eval(x_vals[i-1]))**2 + (mcRatioGraph.GetErrorYlow(i-1) / mcRatioGraph.Eval(x_vals[i-1]))**2)
        except:
            temp_err_high = 1
            temp_err_low = 1

        temp_err_x_low = dataRatioGraph.GetErrorXlow(i-1)
        temp_err_x_high = dataRatioGraph.GetErrorXhigh(i-1)

        temp_y_val = dataRatioGraph.Eval(x_vals[i-1]) / mcRatioGraph.Eval(x_vals[i-1])
        #print(temp_y_val, x_vals[i-1], temp_err_high, temp_err_low, temp_err_x_high, temp_err_x_low)

        if math.isnan(temp_err_high):
            temp_err_high = 0
       
        if "nonisomuon" in dataset and "eta" not in XTitle:
            print("hi")
 
        dataMcRatio.SetPoint(i-1, x_vals[i-1], temp_y_val)
        dataMcRatio.SetPointEYhigh(i-1, temp_err_high)
        dataMcRatio.SetPointEYlow(i-1, temp_err_low)
        dataMcRatio.SetPointEXhigh(i-1, temp_err_x_high)
        dataMcRatio.SetPointEXlow(i-1, temp_err_x_low)
        
    #dataMcRatio.GetXaxis().SetLimits(x_vals[0]-dataMcRatio.GetErrorXlow(0), x_vals[len(x_vals)-1]+dataMcRatio.GetErrorXlow(len(x_vals)-1))
    if "eta" in XTitle:
        dataMcRatio.GetXaxis().SetLimits(-2.9,2.9)
    else:
        dataMcRatio.GetXaxis().SetLimits(30,320)

    #dataMcRatio_top = dataNum.Clone()
    #dataMcRatio_top.Divide(dataDen.Clone())
    #dataMcRatio_bottom = mcNum.Clone()
    #dataMcRatio_bottom.Divide(mcDen.Clone())
    #dataMcRatio_top.Divide(dataMcRatio_bottom.Clone())
    #dataMcRatio = dataMcRatio_top.Clone()
    dataMcRatio.SetMinimum(0.8)
    dataMcRatio.SetMaximum(1.2)
    dataMcRatio.GetYaxis().SetNdivisions(-304)
    dataMcRatio.SetTitle("")
    dataMcRatio.SetLineWidth(3)
    dataMcRatio.SetMarkerSize(0.5)
    dataMcRatio.SetMarkerStyle(20)
    dataMcRatio.SetMarkerColor(dataMcRatio.GetLineColor())
    dataMcRatio.GetXaxis().SetTitleSize(PadFactor*XTitleSize)
    dataMcRatio.GetXaxis().SetLabelSize(PadFactor*XLabelSize)
    dataMcRatio.GetXaxis().SetTitleOffset(0.65*XTitleOffset/PadFactor)
    dataMcRatio.GetXaxis().SetTitle(XTitle)
    dataMcRatio.GetYaxis().SetTitle("Data / MC")
    dataMcRatio.GetYaxis().SetTitleSize(0.8*PadFactor*YTitleSize)
    dataMcRatio.GetYaxis().SetLabelSize(PadFactor*YLabelSize)
    dataMcRatio.GetYaxis().SetTitleOffset(1.4*YTitleOffset/PadFactor)
 
    # fill Eficiency and SF to canvas
    theName = dataNum.GetName().replace("h_num", year)
    aCanvas = ROOT.TCanvas("c_%s"%(theName), "c_%s"%(theName), 1400, 1400)
    ROOT.gPad.Clear()
    aCanvas.Divide(1,2)
    aCanvas.cd(1)
    ROOT.gPad.SetPad(XMin, YMin, XMax, YMax)
    ROOT.gPad.SetGridy()
    ROOT.gPad.SetGridx()
    TopMargin    = 0.1
    LeftMargin   = 0.17
    BottomMargin = 0.01
    RightMargin  = 0.05
    ROOT.gPad.SetTopMargin(TopMargin)
    ROOT.gPad.SetLeftMargin(LeftMargin)
    ROOT.gPad.SetBottomMargin(BottomMargin)
    ROOT.gPad.SetRightMargin(RightMargin)
    ROOT.gPad.SetTicks()

    mcRatio.Draw("AP")
    dataRatio.Draw("P SAME")
    aCanvas.Update()

    # Formating TEfficiency via TGraphAsymmErrors
    mcGraph = mcRatio.GetPaintedGraph()
    mcGraph.SetTitle("")
    mcGraph.GetYaxis().SetTitle("L1+HLT Efficiency")
    mcGraph.GetXaxis().SetLabelSize(0)
    mcGraph.GetXaxis().SetRangeUser(30,320)
    mcGraph.GetYaxis().SetTitleSize(0.05)
    mcGraph.GetYaxis().SetLabelSize(0.05)
    mcGraph.GetYaxis().SetRangeUser(0.45,1.1)
    #mcGraph.GetXaxis().SetRange(1,mcGraph.GetNbinsX()+1)

    dataGraph = dataRatio.GetPaintedGraph()
    dataGraph.SetTitle("")
    dataGraph.GetYaxis().SetTitle("L1+HLT Efficiency")
    dataGraph.GetXaxis().SetLabelSize(0)
    dataGraph.GetXaxis().SetRangeUser(30,320)
    dataGraph.GetYaxis().SetTitleSize(0.05)
    dataGraph.GetYaxis().SetLabelSize(0.05)
    dataGraph.GetYaxis().SetRangeUser(0.45,1.1)
    #dataGraph.GetXaxis().SetRange(1,dataGraph.GetNbinsX()+1)

    aCanvas.Modified()
    aCanvas.Update()

    iamLegend = ROOT.TLegend(0.63, 0.05, 0.99, 0.27, "", "trNDC")
    iamLegend.SetTextSize(0.035)
    iamLegend.SetBorderSize(0)
    iamLegend.SetFillStyle(0)
    if dataset == "electron":
        iamLegend.AddEntry(dataRatio, "Data_SingleElectron", "E1 P")
    elif dataset == "muon":
        iamLegend.AddEntry(dataRatio, "Data_SingleMuon", "E1 P")
    iamLegend.AddEntry(mcRatio, "TT", "E1 P")
    iamLegend.Draw("SAME")
    addCMSlogo(aCanvas, year, TopMargin=0.098, LeftMargin=0.17, RightMargin=0.05, SF=1.0)
    
    aCanvas.cd(2)
    ROOT.gPad.SetPad(RatioXMin, RatioYMin, RatioXMax, RatioYMax)
    ROOT.gPad.SetTopMargin(0.05)
    ROOT.gPad.SetBottomMargin(0.35)
    ROOT.gPad.SetRightMargin(RightMargin)
    ROOT.gPad.SetLeftMargin(LeftMargin)
    ROOT.gPad.SetTicks()
    ROOT.gPad.SetGridy()
    dataMcRatio.Draw("APE")
    dataMcRatio.GetXaxis().SetRangeUser(30,320)

    aCanvas.Modified()
    aCanvas.Update()
    aCanvas.SaveAs("triggerEfficiencySF_plots/%s/%s_TriggerEfficiency.pdf"%(tag,theName))

# --------------------------
# Make 2D Scale Factor plots
# --------------------------
def make_ScaleFactor_Plots(dataNum, dataDen, mcNum, mcDen, aName, outputFile):

    XTitle = ""
    if "Pt" in dataNum.GetName(): XTitle = "p_{T} [GeV]"
    else: XTitle = "#eta"

    XMin = 0;    XMax = 1; RatioXMin = 0; RatioXMax = 1
    YMin = 0.20; YMax = 1; RatioYMin = 0; RatioYMax = 0.20
    PadFactor  = (YMax-YMin) / (RatioYMax-RatioYMin)
    XTitleSize = 0.045;  XLabelSize = 0.036;  XTitleOffset = 1.0 
    YTitleSize = 0.045;  YLabelSize = 0.036;  YTitleOffset = 0.9

    # get data and mc efficiency 
    dataRatio = dataNum.Clone()
    dataRatio.Divide(dataNum, dataDen, 1, 1, "B")
    mcRatio = mcNum.Clone()
    mcRatio.Divide(mcNum, mcDen, 1, 1, "B")
    # get errors
    updateErrors(dataRatio,dataDen)
    updateErrors(mcRatio,mcDen)
    # get scale factor 
    dataMcRatio = dataNum.Clone()
    dataMcRatio.Divide(dataRatio, mcRatio)

    dataRatio.GetZaxis().SetRangeUser(0.7, 1.1)
    dataRatio.SetTitle(""); mcRatio.SetTitle("")
    dataRatio.GetYaxis().SetTitle("#eta"); dataRatio.GetXaxis().SetTitle("p_{T} [GeV]")
    dataRatio.GetYaxis().SetTitleSize(YTitleSize); dataRatio.GetXaxis().SetTitleSize(XTitleSize)
    dataRatio.GetYaxis().SetLabelSize(YLabelSize); dataRatio.GetXaxis().SetLabelSize(XLabelSize)
    dataRatio.GetYaxis().SetTitleOffset(YTitleOffset); dataRatio.GetXaxis().SetTitleOffset(XTitleOffset)
    dataRatio.SetContour(255)

    mcRatio.GetZaxis().SetRangeUser(0.7, 1.1)
    mcRatio.SetTitle(""); mcRatio.SetTitle("")
    mcRatio.GetYaxis().SetTitle("#eta"); mcRatio.GetXaxis().SetTitle("p_{T} [GeV]")
    mcRatio.GetYaxis().SetTitleSize(YTitleSize); mcRatio.GetXaxis().SetTitleSize(XTitleSize)
    mcRatio.GetYaxis().SetLabelSize(YLabelSize); mcRatio.GetXaxis().SetLabelSize(XLabelSize)
    mcRatio.GetYaxis().SetTitleOffset(YTitleOffset); mcRatio.GetXaxis().SetTitleOffset(XTitleOffset)
    mcRatio.SetContour(255)

    dataMcRatio.GetZaxis().SetRangeUser(0.7, 1.1)
    dataMcRatio.SetTitle(""); dataMcRatio.SetTitle("")
    dataMcRatio.GetYaxis().SetTitle("#eta"); dataMcRatio.GetXaxis().SetTitle("p_{T} [GeV]")
    dataMcRatio.GetYaxis().SetTitleSize(YTitleSize); dataMcRatio.GetXaxis().SetTitleSize(XTitleSize)
    dataMcRatio.GetYaxis().SetLabelSize(YLabelSize); dataMcRatio.GetXaxis().SetLabelSize(XLabelSize)
    dataMcRatio.GetYaxis().SetTitleOffset(YTitleOffset); dataMcRatio.GetXaxis().SetTitleOffset(XTitleOffset)
    dataMcRatio.SetContour(255)

    #setMinimumErrors(dataMcRatio)

    if "ge5jetCut" in aName and "trig" in aName:
        outputFile.cd()
        dataMcRatio.SetName(aName+"_TriggerSF")
        dataMcRatio.Write(aName+"_TriggerSF")

    theName = dataNum.GetName().replace("h2_num","")
    TopMargin    = 0.098
    LeftMargin   = 0.09
    BottomMargin = 0.11
    RightMargin  = 0.11

    #aCanvas = ROOT.TCanvas("c_data_%s"%(theName), "c_data_%s"%(theName), 1400, 1400)
    #ROOT.gPad.Clear()
    #aCanvas.cd()
    #ROOT.gPad.SetTopMargin(TopMargin)
    #ROOT.gPad.SetLeftMargin(LeftMargin)
    #ROOT.gPad.SetBottomMargin(BottomMargin)
    #ROOT.gPad.SetRightMargin(RightMargin)
    #ROOT.gPad.SetTicks()
    #dataRatio.Draw("COLZ E TEXT")
    #addCMSlogo(aCanvas, year, TopMargin=0.12, LeftMargin=0.09, RightMargin=0.11, SF=1.0)
    #aCanvas.SaveAs("triggerEfficiencySF_plots/%s/%s_data_%s.pdf"%(tag,year,theName))

    #aCanvas = ROOT.TCanvas("c_mc_%s"%(theName), "c_mc_%s"%(theName), 1500, 1200)
    #ROOT.gPad.Clear()
    #aCanvas.cd()
    #ROOT.gPad.SetTopMargin(TopMargin)
    #ROOT.gPad.SetLeftMargin(LeftMargin)
    #ROOT.gPad.SetBottomMargin(BottomMargin)
    #ROOT.gPad.SetRightMargin(RightMargin)
    #ROOT.gPad.SetTicks()
    #mcRatio.Draw("COLZ E TEXT")
    #addCMSlogo(aCanvas, year, TopMargin=0.098, LeftMargin=0.17, RightMargin=0.05, SF=1.0)
    #aCanvas.SaveAs("triggerEfficiencySF_plots/%s/%s_mc_%s.pdf"%(tag,year,theName))

    aCanvas = ROOT.TCanvas("c_datamc_%s"%(theName), "c_datamc_%s"%(theName), 1400, 1400)
    ROOT.gPad.Clear()
    aCanvas.cd()
    ROOT.gPad.SetTopMargin(TopMargin)
    ROOT.gPad.SetLeftMargin(LeftMargin)
    ROOT.gPad.SetBottomMargin(BottomMargin)
    ROOT.gPad.SetRightMargin(RightMargin)
    ROOT.gPad.SetTicks()
    dataMcRatio.Draw("COLZ E TEXT")
    addCMSlogo(aCanvas, year, TopMargin=0.12, LeftMargin=0.09, RightMargin=0.11, SF=1.0)
    if not os.path.isdir("./triggerEfficiencySF_plots/%s"%(tag)):
        os.makedirs("./triggerEfficiencySF_plots/%s"%(tag))
    aCanvas.SaveAs("triggerEfficiencySF_plots/%s/%s%s_TriggerSF.pdf"%(tag,year,theName))
    print("made pdf: {}".format("triggerEfficiencySF_plots/%s/%s%s_TriggerSF.pdf"%(tag,year,theName)))

def GetOverflow(hist):

    if "Eta" in hist.GetName():
        return hist

    bin_edges = []

    for i in range(1,hist.GetNbinsX()+1):
        bin_edges.append(hist.GetBinLowEdge(i))

    bin_edges.append(bin_edges[-1] + 40)
    
    bin_content = []
    bin_errors = []
    for i in range(1, hist.GetNbinsX()+2):
        bin_content.append(hist.GetBinContent(i))
        bin_errors.append(hist.GetBinError(i))

    print(bin_edges)

    return {'content': bin_content, 'error': bin_errors, "edges": bin_edges}


if __name__ == '__main__':

    effTags      = [ "den", "num" ] 
    #lepTags      = ["nonIsoMu"]
    lepTags      = [ "el", "mu"]
    binTags      = [ "Eta", "Pt"  ]
    ptTags       = [ "pt40"       ]
    trigTags     = [ "trig"       ] 
    nJetCutTags  = [ "ge1jetCut", "ge5jetCut", "le3jetCut", "le5jetCut" ]

    basePath = arg.basePath
    inputDir = arg.inputDir
    tag      = arg.tag
    year     = arg.year
    dataset  = arg.dataset

    fullPath = basePath + "/" + inputDir
    if dataset == "nonisomuon":
        mcFile   = ROOT.TFile.Open(fullPath + "/" + year + "_QCD.root")
    else:
        mcFile   = ROOT.TFile.Open(fullPath + "/" + year + "_TT.root")
    dataFile = 0
    if dataset == "electron": 
        dataFile = ROOT.TFile.Open(fullPath + "/" + year + "_Data_SingleElectron.root")
    elif dataset == "muon": 
        dataFile = ROOT.TFile.Open(fullPath + "/" + year + "_Data_SingleMuon.root")
    elif dataset == "nonisomuon": 
        dataFile = ROOT.TFile.Open(fullPath + "/" + year + "_Data_SingleElectron.root")

    outputPath     = "./triggerEfficiencySF_plots/triggerEffSFS_1l/"
    outputFileName = "%s_Leptonic_Triggers_SF.root"%(year)
    outputFile     = ROOT.TFile.Open("%s/%s"%(outputPath,outputFileName), "WRITE")

    
    for lepTag in lepTags:
    
        if (lepTag == "el" and dataset == "electron") or (lepTag == "mu" and dataset == "muon") or (lepTag == "nonisomuon" and dataset == "muon"): continue
    
        for ptTag in ptTags:
    
            for trigTag in trigTags:
    
                for nJetCutTag in nJetCutTags:
   
                    # make scale factor plots 
                    denomTag2 = "h2_den_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"_wLepPtLepEtaBin"
                    numerTag2 = "h2_num_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"_wLepPtLepEtaBin"
                    goodName  = year+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"_wLepPtLepEtaBin"

                    dataFile.cd()
                    h2DataNum = dataFile.Get(numerTag2)
                    h2DataDen = dataFile.Get(denomTag2)
                    mcFile.cd()
                    h2McNum = mcFile.Get(numerTag2)
                    h2McDen = mcFile.Get(denomTag2)

                    if h2DataNum == None or h2DataDen == None or h2McNum == None or h2McDen == None: continue

                    make_ScaleFactor_Plots(h2DataNum, h2DataDen, h2McNum, h2McDen, goodName, outputFile)

                    for binTag in binTags:

                        # make efficiency plots
                        denomTag1 = "h_den_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"_wLep%sBin"%(binTag)
                        numerTag1 = "h_num_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"_wLep%sBin"%(binTag)

                        dataFile.cd()
                        hDataNum = dataFile.Get(numerTag1)
                        hDataDen = dataFile.Get(denomTag1)
                        mcFile.cd()
                        hMcNum = mcFile.Get(numerTag1)
                        hMcDen = mcFile.Get(denomTag1)
                        #hDataNum_dict = GetOverflow(hDataNum)
                        #hDataDen_dict = GetOverflow(hDataDen)
                        #hMcNum_dict = GetOverflow(hMcNum)
                        #hMcDen_dict = GetOverflow(hMcDen)

                        #if not "Eta" in hDataNum.GetTitle():

                        #    hDataNum = ROOT.TH1D("", "", 10, array('d', hDataNum_dict["edges"]))
                        #    #hDataNum = ROOT.TH1D("{}".format(hDataNum.GetTitle()), "{}".format(hDataNum.GetTitle()), 10, array('d', hDataNum_dict["edges"]))
                        #    for i in range(hDataNum.GetNbinsX()+1):
                        #        hDataNum.SetBinContent(i, hDataNum_dict['content'][i-1])
                        #        hDataNum.SetBinError(i, hDataNum_dict['error'][i-1])

                        #    hDataDen = ROOT.TH1D("", "", 10, array('d', hDataDen_dict["edges"]))
                        #    #hDataDen = ROOT.TH1D("{}".format(hDataDen.GetTitle()), "{}".format(hDataDen.GetTitle()), c_int(hDataDen.GetNbinsX()+1), array('d', hDataDen_dict["edges"]))
                        #    for i in range(hDataDen.GetNbinsX()+1):
                        #        hDataDen.SetBinContent(i, hDataDen_dict['content'][i-1])
                        #        hDataDen.SetBinError(i, hDataDen_dict['error'][i-1])

                        #    hMcNum = ROOT.TH1D("", "", 10, array('d', hDataDen_dict["edges"]))
                        #    #hMcNum = ROOT.TH1D("{}".format(hMcNum.GetTitle()), "{}".format(hMcNum.GetTitle()), c_int(hMcNum.GetNbinsX()+1), array('d', hMcNum_dict["edges"]))
                        #    for i in range(hMcNum.GetNbinsX()+1):
                        #        hMcNum.SetBinContent(i, hMcNum_dict['content'][i-1])
                        #        hMcNum.SetBinError(i, hMcNum_dict['error'][i-1])
                        #        
                        #    hMcDen = ROOT.TH1D("", "", 10, array('d', hDataDen_dict["edges"]))
                        #    #hMcDen = ROOT.TH1D("{}".format(hMcDen.GetTitle()), "{}".format(hMcDen.GetTitle()), c_int(hMcDen.GetNbinsX()+1), array('d', hMcDen_dict["edges"]))
                        #    for i in range(hMcDen.GetNbinsX()+1):
                        #        hMcDen.SetBinContent(i, hMcDen_dict['content'][i-1])
                        #        hMcDen.SetBinError(i, hMcDen_dict['error'][i-1])

                        if hDataNum == None or hDataDen == None or hMcNum == None or hMcDen == None: continue

                        make_Efficiency_Plots(hDataNum, hDataDen, hMcNum, hMcDen) 
