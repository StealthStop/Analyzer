#! /bin/env/python

# This script makes Fig. 5 of the SUS-19-004 paper
# Inputs are the output from fit diagnostics in Higgs Combine

import ROOT, os, math, array, copy, argparse

ROOT.TH1.AddDirectory(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH2.SetDefaultSumw2()
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat("")
ROOT.gStyle.SetPaintTextFormat("3.2f")
ROOT.gStyle.SetFrameLineWidth(2)
ROOT.gStyle.SetEndErrorSize(0)

usage = "usage: %prog [options]"
parser = argparse.ArgumentParser(usage)
parser.add_argument("--approved", dest="approved", help="Plot is approved",           action="store_true", default=False) 
parser.add_argument("--sbfit",    dest="sbfit",    help="Results for s+b fit",        action="store_true", default=False)
parser.add_argument("--inputDir", dest="inputDir", help="Input dir with fit results", type=str,            required=True)
args = parser.parse_args()

fitStr1 = "fit_b"; fitStr2 = "_bonly"
if args.sbfit:
    fitStr1 = "fit_s"
    fitStr2 = ""

inputDir = args.inputDir

outpath = "./PlotsForPaper/"
if not os.path.exists(outpath): os.makedirs(outpath)

yearFile = {"2016" : {}, "2017" : {}, "2018pre" : {}, "2018post" : {}}

masses = ["350", "400", "550", "850"]

for year in yearFile.keys():
    for mass in masses:
        theFile = "%s/Fit_Data_%s/output-files/RPV_%s_%s/fitDiagnostics%sRPV%s.root"%(inputDir,year,mass,year,year,mass)
        yearFile[year][mass] = ROOT.TFile.Open(theFile, "READ")

        theFile = "%s/Fit_Data_%s/output-files/SYY_%s_%s/fitDiagnostics%sSYY%s.root"%(inputDir,year,mass,year,year,mass)
        yearFile[year]["S"+mass] = ROOT.TFile.Open(theFile, "READ")

mvaBins = ["D1", "D2", "D3", "D4"]
aliases = {"TT"     : "t#bar{t}    ",
           "QCD"    : "QCD multijet",
           "TTX"    : "t#bar{t} + X   ",
           "OTHER"  : "Other   ",
           "SIG"    : "Fit Signal",
           "SIG1"   : "RPV m_{ #tilde{t}} = 350 GeV",
           "SIG2"   : "RPV m_{ #tilde{t}} = 850 GeV",
           "SIG3"   : "Stealth SYY m_{ #tilde{t}} = 350 GeV",
           "SIG4"   : "Stealth SYY m_{ #tilde{t}} = 850 GeV"
     
}

procDictionary = {
                   "TTX"   : {"graph" : 0, "h" : 0, "lstyle" : 0, "lsize" : 0, "msize" : 0, "color" : ROOT.kOrange+2},
                   "QCD"   : {"graph" : 0, "h" : 0, "lstyle" : 0, "lsize" : 0, "msize" : 0, "color" : ROOT.kGreen+1},
                   "OTHER" : {"graph" : 0, "h" : 0, "lstyle" : 0, "lsize" : 0, "msize" : 0, "color" : ROOT.kMagenta+2},
                   "TT"    : {"graph" : 0, "h" : 0, "lstyle" : 0, "lsize" : 0, "msize" : 0, "color" : ROOT.kBlue-6},
                   "DATA"  : {"graph" : 0, "h" : 0, "lstyle" : 0, "lsize" : 5, "msize" : 5, "color" : ROOT.kBlack},
                   "SYST"  : {"graph" : 0, "h" : 0, "lstyle" : 0, "lsize" : 0, "msize" : 0, "color" : ROOT.kBlack},
                   "SIG"   : {"graph" : 0, "h" : 0, "lstyle" : 0, "lsize" : 0, "msize" : 0, "color" : ROOT.kGray+1},
                   "SIG1"  : {"graph" : 0, "h" : 0, "lstyle" : 2, "lsize" : 4, "msize" : 0, "color" : ROOT.kMagenta},
                   "SIG2"  : {"graph" : 0, "h" : 0, "lstyle" : 3, "lsize" : 4, "msize" : 0, "color" : ROOT.kRed},
                   "SIG3"  : {"graph" : 0, "h" : 0, "lstyle" : 1, "lsize" : 4, "msize" : 0, "color" : ROOT.kCyan+1},
                   "SIG4"  : {"graph" : 0, "h" : 0, "lstyle" : 4, "lsize" : 4, "msize" : 0, "color" : ROOT.kBlue}
}

mvaDictionary = {"D1" : copy.deepcopy(procDictionary),
                 "D2" : copy.deepcopy(procDictionary),
                 "D3" : copy.deepcopy(procDictionary),
                 "D4" : copy.deepcopy(procDictionary)
}

theDictionary = {"2016"     : copy.deepcopy(mvaDictionary),
                 "2017"     : copy.deepcopy(mvaDictionary),
                 "2018pre"  : copy.deepcopy(mvaDictionary),
                 "2018post" : copy.deepcopy(mvaDictionary)
}

perMVADictionary = {"D1" : copy.deepcopy(procDictionary),
                    "D2" : copy.deepcopy(procDictionary),
                    "D3" : copy.deepcopy(procDictionary),
                    "D4" : copy.deepcopy(procDictionary)
}

run2Dictionary = copy.deepcopy(procDictionary)

for year, mvaDict in theDictionary.iteritems():

    for mva, procDict in mvaDict.iteritems(): 

        for process, d in procDict.iteritems():

            if   process == "SIG1":
                path = "shapes_prefit/%s/total_signal"%(mva)

                d["h"] = yearFile[year]["350"].Get(path); d["h"].SetDirectory(0); d["h"].SetName("%s_%s_RPV350"%(year,mva))

            elif process == "SIG2":
                path = "shapes_prefit/%s/total_signal"%(mva)

                d["h"] = yearFile[year]["850"].Get(path); d["h"].SetDirectory(0); d["h"].SetName("%s_%s_RPV850"%(year,mva))

            elif process == "SIG3":
                 path = "shapes_prefit/%s/total_signal"%(mva)

                 d["h"] = yearFile[year]["S350"].Get(path); d["h"].SetDirectory(0); d["h"].SetName("%s_%s_SYY350"%(year,mva))
       
            elif process == "SIG4":
                 path = "shapes_prefit/%s/total_signal"%(mva)

                 d["h"] = yearFile[year]["S850"].Get(path); d["h"].SetDirectory(0); d["h"].SetName("%s_%s_SYY850"%(year,mva))

            elif process == "SYST":
    
                rooPlot = yearFile[year]["350"].Get("%s_CMS_th1x_%s"%(mva,fitStr1))
                rooCurveErr = rooPlot.getCurve("pdf_bin%s%s_Norm[CMS_th1x]_errorband"%(mva,fitStr2))
                rooCurveVal = rooPlot.getCurve("pdf_bin%s%s_Norm[CMS_th1x]"%(mva,fitStr2))
                
                xband = rooCurveErr.GetX(); yband = rooCurveErr.GetY(); npband = rooCurveErr.GetN()
                x = rooCurveVal.GetX(); y = rooCurveVal.GetY(); np = rooCurveVal.GetN()
                yErrUp = array.array('d', [0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
                yErrDown = array.array('d', [0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
                xErr = array.array('d', [0.5, 0.5, 0.5, 0.5, 0.5, 0.5])
                hx = array.array('d', [0.5, 1.5, 2.5, 3.5, 4.5, 5.5])
                hy = array.array('d', [0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
                
                bins = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5] 
                for i in bins:
                    idx = bins.index(i)
                    out = False
                    for aX in xrange(0,npband/2):
                        if out: break
                        for bX in xrange(0, np):
                            if x[bX] == i and xband[aX] == i:
                                yErrUp[idx] = yband[aX] - y[bX]
                                hy[idx] = y[bX]
                                out = True
                                break
                
                for i in bins:
                    idx = bins.index(i)
                    out = False
                    for aX in xrange(npband/2,npband):
                        if out: break
                        for bX in xrange(0,np):
                            if x[bX] == i and xband[aX] == i:
                                yErrDown[idx] = y[bX]-yband[aX]
                                out = True
                                break
                       
                d["graph"] = ROOT.TGraphAsymmErrors(6, hx, hy, xErr, xErr, yErrDown, yErrUp)

            elif process == "DATA":
    
                rooPlot = yearFile[year]["350"].Get("%s_CMS_th1x_fit_b"%(mva))
                rooHist = rooPlot.getHist("h_%s"%(mva))
    
                x = rooHist.GetX()
                y = rooHist.GetY()
                xErr = array.array('d', [0.5, 0.5, 0.5, 0.5, 0.5, 0.5])
                yErr = array.array('d', [math.sqrt(i) for i in y])
    
                d["graph"] = ROOT.TGraphAsymmErrors(6, x, y, xErr, xErr, yErr, yErr)
    
                d["h"] = ROOT.TH1F("%s_%s_data"%(year,mva), "%s_%s_data"%(year,mva), 6, 0, 6)
                d["h"].SetDirectory(0)
                x = ROOT.Double(0.0); y = ROOT.Double(0.0)
                for iP in xrange(0,d["graph"].GetN()):
                    
                    theBin = d["h"].GetXaxis().FindBin(iP)
    
                    d["graph"].GetPoint(iP, x, y)
                    
                    d["h"].Fill(x,y)
                    d["h"].SetBinError(theBin, math.sqrt(y))

            else:
                path = "shapes_%s/%s/%s"%(fitStr1,mva,process)

                d["h"] = yearFile[year]["400"].Get(path); d["h"].SetDirectory(0); d["h"].SetName("%s_%s_%s"%(year,mva,process)); d["h"].Sumw2()
    

for year, d in yearFile.iteritems():
    for mass, f in d.iteritems():
        f.Close()

for mva in mvaBins:
    for process in ["TTX", "OTHER", "QCD", "TT", "DATA", "SYST", "SIG", "SIG1", "SIG2", "SIG3", "SIG4"]:
        for year in ["2016", "2017", "2018pre", "2018post"]:

            histo = "h"
            if process == "SYST": histo = "graph"

            if perMVADictionary[mva][process][histo] == 0: perMVADictionary[mva][process][histo] = theDictionary[year][mva][process][histo].Clone("perMVA_%s_%s_clone"%(mva,process))
        
            else:
                if process == "SYST":

                    xOld = ROOT.Double(0.0); yOld = ROOT.Double(0.0)
                    xNew = ROOT.Double(0.0); yNew = ROOT.Double(0.0)
                    for i in xrange(0, theDictionary[year][mva][process][histo].GetN()):

                        perMVADictionary[mva][process][histo].GetPoint(i,xOld,yOld)
                        theDictionary[year][mva][process][histo].GetPoint(i,xNew,yNew)

                        oldUpErr = perMVADictionary[mva][process][histo].GetErrorYhigh(i)
                        oldDownErr = perMVADictionary[mva][process][histo].GetErrorYlow(i)

                        newUpErr = theDictionary[year][mva][process][histo].GetErrorYhigh(i)
                        newDownErr = theDictionary[year][mva][process][histo].GetErrorYlow(i)

                        perMVADictionary[mva][process][histo].SetPoint(i,xNew,yNew+yOld)

                        perMVADictionary[mva][process][histo].SetPointEYhigh(i,(oldUpErr**2.0+newUpErr**2.0)**0.5)
                        perMVADictionary[mva][process][histo].SetPointEYlow(i,(oldDownErr**2.0+newDownErr**2.0)**0.5)
    
                else: perMVADictionary[mva][process][histo].Add(theDictionary[year][mva][process][histo])

            if run2Dictionary[process][histo] == 0: run2Dictionary[process][histo] = theDictionary[year][mva][process][histo].Clone("run2_%s_%s_clone"%(mva,process))
            else:
                if process == "SYST":

                    xOld = ROOT.Double(0.0); yOld = ROOT.Double(0.0)
                    xNew = ROOT.Double(0.0); yNew = ROOT.Double(0.0)
                    for i in xrange(0, theDictionary[year][mva][process][histo].GetN()):
                    
                        run2Dictionary[process][histo].GetPoint(i,xOld,yOld)
                        theDictionary[year][mva][process][histo].GetPoint(i,xNew,yNew)

                        oldUpErr = run2Dictionary[process][histo].GetErrorYhigh(i)
                        oldDownErr = run2Dictionary[process][histo].GetErrorYlow(i)

                        newUpErr = theDictionary[year][mva][process][histo].GetErrorYhigh(i)
                        newDownErr = theDictionary[year][mva][process][histo].GetErrorYlow(i)

                        run2Dictionary[process][histo].SetPoint(i,xNew,yNew+yOld)
                        run2Dictionary[process][histo].SetPointEYhigh(i,(oldUpErr**2.0+newUpErr**2.0)**0.5)
                        run2Dictionary[process][histo].SetPointEYlow(i,(oldDownErr**2.0+newDownErr**2.0)**0.5)
    
                else: run2Dictionary[process][histo].Add(theDictionary[year][mva][process][histo])

            if (process != "DATA" and process != "SYST" and "SIG" not in process) or process == "SIG": run2Dictionary[process][histo].SetFillColor(theDictionary[year][mva][process]["color"])
            if process == "SYST":
                run2Dictionary[process][histo].SetFillStyle(3004)
                run2Dictionary[process][histo].SetFillColor(ROOT.kBlack)

            run2Dictionary[process][histo].SetMarkerColor(theDictionary[year][mva][process]["color"])
            run2Dictionary[process][histo].SetLineColor(theDictionary[year][mva][process]["color"])
    
            run2Dictionary[process][histo].SetMarkerSize(theDictionary[year][mva][process]["msize"])
            run2Dictionary[process][histo].SetMarkerStyle(20)
            run2Dictionary[process][histo].SetLineWidth(theDictionary[year][mva][process]["lsize"])
            if "SIG" in process: run2Dictionary[process][histo].SetLineStyle(theDictionary[year][mva][process]["lstyle"])

            if (process != "DATA" and process != "SYST" and "SIG" not in process) or process == "SIG": perMVADictionary[mva][process][histo].SetFillColor(theDictionary[year][mva][process]["color"])
            if process == "SYST":
                perMVADictionary[mva][process][histo].SetFillStyle(3004)
                perMVADictionary[mva][process][histo].SetFillColor(ROOT.kBlack)

            perMVADictionary[mva][process][histo].SetMarkerColor(theDictionary[year][mva][process]["color"])
            perMVADictionary[mva][process][histo].SetLineColor(theDictionary[year][mva][process]["color"])
    
            perMVADictionary[mva][process][histo].SetMarkerSize(theDictionary[year][mva][process]["msize"])
            perMVADictionary[mva][process][histo].SetMarkerStyle(20)
            perMVADictionary[mva][process][histo].SetLineWidth(theDictionary[year][mva][process]["lsize"])
            if "SIG" in process: perMVADictionary[mva][process][histo].SetLineStyle(theDictionary[year][mva][process]["lstyle"])

for mva in mvaBins:

    c = ROOT.TCanvas("%s_c"%(mva), "%s_c"%(mva), 2400, 2400)
    totalBkgd = ROOT.TH1F("%s_totalBkgd"%(mva), "%s_totalBkgd"%(mva), 6, 0, 6)
    totalRatio = ROOT.TH1F("%s_totalRatio"%(mva), "%s_totalRatio"%(mva), 6, 0, 6)
    theStack = ROOT.THStack("%s_bkgd_stack"%(mva), "%s_bkgd_stack"%(mva))

    systRatio = perMVADictionary[mva]["SYST"]["graph"].Clone("%s_ratioSyst_clone"%(mva))

    for process in ["SIG", "TTX", "OTHER", "QCD", "TT"]:
        totalBkgd.Add(perMVADictionary[mva][process]["h"])
        theStack.Add(perMVADictionary[mva][process]["h"])

    mcy = ROOT.Double(0.0)
    datax = ROOT.Double(0.0); datay = ROOT.Double(0.0)
    yErrUp = ROOT.Double(0.0); yErrDown = ROOT.Double(0.0)
    for i in xrange(0, systRatio.GetN()):

        mcy = totalBkgd.GetBinContent(i+1)

        perMVADictionary[mva]["SYST"]["graph"].GetPoint(i, datax, datay)
        yErrUp = perMVADictionary[mva]["SYST"]["graph"].GetErrorYhigh(i)
        yErrDown = perMVADictionary[mva]["SYST"]["graph"].GetErrorYlow(i)

        systRatio.SetPoint(i, datax, 0.0)
        systRatio.SetPointEYhigh(i, yErrUp/datay**0.5)
        systRatio.SetPointEYlow(i, yErrDown/datay**0.5)
      
    for xbin in xrange(1, totalBkgd.GetNbinsX()+1):
        totalRatio.SetBinContent(xbin, (perMVADictionary[mva]["DATA"]["h"].GetBinContent(xbin)-totalBkgd.GetBinContent(xbin))/perMVADictionary[mva]["DATA"]["h"].GetBinContent(xbin)**0.5)
        totalRatio.SetBinError(xbin, 1.0)
    
    XMin = 0;    XMax = 1; RatioXMin = 0; RatioXMax = 1 
    YMin = 0.30; YMax = 1; RatioYMin = 0; RatioYMax = 0.30
    PadFactor = (YMax-YMin) / (RatioYMax-RatioYMin)
    c.Divide(1,2); c.cd(1); ROOT.gPad.SetLogy(); ROOT.gPad.SetLogz(); ROOT.gPad.SetPad(XMin, YMin, XMax, YMax)
    
    ROOT.gPad.SetTopMargin(0.10)
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetBottomMargin(0.00)
    ROOT.gPad.SetRightMargin(0.04)
    
    iamLegend = ROOT.TLegend(ROOT.gPad.GetLeftMargin() + 0.02, 0.70 - ROOT.gPad.GetTopMargin(), ROOT.gPad.GetLeftMargin() + 0.47, 1 - ROOT.gPad.GetTopMargin() - 0.02)
    
    iamLegend.SetTextSize(0.045)
    iamLegend.SetNColumns(2)
    iamLegend.SetBorderSize(2)
    
    sigLegend = ROOT.TLegend(ROOT.gPad.GetLeftMargin() + 0.50, 0.70 - ROOT.gPad.GetTopMargin(), 1. - ROOT.gPad.GetRightMargin() - 0.05, 1 - ROOT.gPad.GetTopMargin() - 0.02)
    sigLegend.SetTextSize(0.045)
    sigLegend.SetNColumns(1)
    sigLegend.SetBorderSize(2)

    for process in ["TTX", "OTHER", "QCD", "TT"]: iamLegend.AddEntry(perMVADictionary[mva][process]["h"], aliases[process], "F")
    for process in ["SIG1", "SIG2", "SIG3", "SIG4"]: sigLegend.AddEntry(perMVADictionary[mva][process]["h"], aliases[process], "L")

    iamLegend.AddEntry(perMVADictionary[mva]["DATA"]["h"], "Data", "EP")

    dummy1 = ROOT.TH1F("dummy1_%s"%(mva), "dummy1_%s"%(mva), 6, 0, 6)
    dummy1.SetTitle("")
    dummy1.GetYaxis().SetTitle("Events / bin")
    dummy1.GetYaxis().SetTitleSize(0.07)
    dummy1.GetXaxis().SetTitleSize(0.07)
    dummy1.GetYaxis().SetLabelSize(0.06)
    dummy1.GetXaxis().SetLabelSize(0.06)
    dummy1.GetYaxis().SetTitleOffset(1.1)
    dummy1.GetXaxis().SetTitleOffset(2.1)
    dummy1.SetMaximum(1e8)
    dummy1.SetMinimum(0.5)

    dummy1.Draw()
    theStack.Draw("SAME")
    perMVADictionary[mva]["DATA"]["h"].Draw("E0P SAME")
    perMVADictionary[mva]["SYST"]["graph"].Draw("2SAME")
    perMVADictionary[mva]["SIG1"]["h"].Draw("SAME")
    perMVADictionary[mva]["SIG2"]["h"].Draw("SAME")
    perMVADictionary[mva]["SIG3"]["h"].Draw("SAME")
    perMVADictionary[mva]["SIG4"]["h"].Draw("SAME")

    iamLegend.Draw("SAME")
    sigLegend.Draw("SAME")

    mark = ROOT.TLatex()
    mark.SetNDC(True)

    mark.SetTextAlign(11);
    mark.SetTextSize(0.065);
    mark.SetTextFont(61);
    mark.DrawLatex(ROOT.gPad.GetLeftMargin(), 1 - (ROOT.gPad.GetTopMargin() - 0.02), "CMS")
    mark.SetTextFont(52);
    if not args.approved: mark.DrawLatex(ROOT.gPad.GetLeftMargin() + 0.103, 1 - (ROOT.gPad.GetTopMargin() - 0.02), "Preliminary")

    mark.SetTextFont(42)
    mark.SetTextAlign(31)
    mark.DrawLatex(1 - ROOT.gPad.GetRightMargin(), 1 - (ROOT.gPad.GetTopMargin() - 0.02), "137.2 fb^{-1} (13 TeV)")
   
    mark2 = ROOT.TLatex()
    mark2.SetNDC(True)
    mark.SetTextAlign(11)
    mark.SetTextSize(0.065)
    mark.SetTextFont(62)
    mark.DrawLatex(ROOT.gPad.GetLeftMargin() + 0.60, 0.65 - ROOT.gPad.GetTopMargin(), "S_{NN,%s}"%(mva[-1]))

    c.cd(2)
    
    ROOT.gPad.SetGridy()
    ROOT.gPad.SetTopMargin(0.00)
    ROOT.gPad.SetBottomMargin(0.35)
    ROOT.gPad.SetRightMargin(0.04)
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetPad(RatioXMin, RatioYMin, RatioXMax, RatioYMax)
    
    for i in xrange(1, totalRatio.GetNbinsX()+1):
        if i != totalRatio.GetNbinsX(): totalRatio.GetXaxis().SetBinLabel(i, str(i+6))
        else:                           totalRatio.GetXaxis().SetBinLabel(i, "#geq %s"%(str(i+6)))

    totalRatio.GetYaxis().SetRangeUser(-1.09, 1.09)
    totalRatio.SetTitle("")
    totalRatio.GetXaxis().SetTitle("Number of jets")
    totalRatio.GetYaxis().SetNdivisions(5, 5, 0)
    totalRatio.GetYaxis().SetTitle("Data / Pred.")
    totalRatio.GetYaxis().SetTitleSize(dummy1.GetYaxis().GetTitleSize()*PadFactor)
    totalRatio.GetXaxis().SetTitleSize(dummy1.GetXaxis().GetTitleSize()*PadFactor)
    totalRatio.GetYaxis().SetLabelSize(dummy1.GetYaxis().GetLabelSize()*PadFactor)
    totalRatio.GetXaxis().SetLabelSize(1.5*dummy1.GetXaxis().GetLabelSize()*PadFactor)
    totalRatio.GetYaxis().SetTitleOffset(0.8*dummy1.GetYaxis().GetTitleOffset()/PadFactor)
    totalRatio.GetXaxis().SetTitleOffset(dummy1.GetXaxis().GetTitleOffset()/PadFactor)
    totalRatio.GetXaxis().SetLabelOffset(0.02)
    totalRatio.SetMarkerStyle(20); totalRatio.SetMarkerSize(5); totalRatio.SetMarkerColor(ROOT.kBlack);
    totalRatio.SetLineWidth(5); totalRatio.SetLineColor(ROOT.kBlack)
    
    totalRatio.Draw("E0P")
    systRatio.Draw("2SAME")
    
    if args.approved:
        c.SaveAs("%s/%s_njets_DataVsMC_%s.pdf"%(outpath,mva,fitStr1))
    else:
        c.SaveAs("%s/%s_njets_DataVsMC_%s_prelim.pdf"%(outpath,mva,fitStr1))

cAll = ROOT.TCanvas("c", "c", 2400, 2400)
totalBkgdAll = ROOT.TH1F("totalBkgdAll", "totalBkgdAll", 6, 0, 6)
totalDataAll = ROOT.TH1F("totalDataAll", "totalDataAll", 6, 0, 6)
totalRatioAll = ROOT.TH1F("totalRatioAll", "totalRatioAll", 6, 0, 6)
theStackAll = ROOT.THStack("theStackAll", "theStackAll")

systRatioAll = run2Dictionary["SYST"]["graph"].Clone("%s_ratioSyst_clone"%(mva))

for process in ["SIG", "TTX", "OTHER", "QCD", "TT"]:
    totalBkgdAll.Add(run2Dictionary[process]["h"])
    theStackAll.Add(run2Dictionary[process]["h"])

for process in ["SIG1", "SIG2", "SIG3", "SIG4"]:
    sigLegend.AddEntry(run2Dictionary[process]["h"], aliases[process], "L")

mcy = ROOT.Double(0.0)
datax = ROOT.Double(0.0); datay = ROOT.Double(0.0)
yErrUp = ROOT.Double(0.0); yErrDown = ROOT.Double(0.0)
for i in xrange(0, systRatioAll.GetN()):

    mcy = totalBkgdAll.GetBinContent(i+1)

    run2Dictionary["SYST"]["graph"].GetPoint(i, datax, datay)
    yErrUp = run2Dictionary["SYST"]["graph"].GetErrorYhigh(i)
    yErrDown = run2Dictionary["SYST"]["graph"].GetErrorYlow(i)

    systRatioAll.SetPoint(i, datax, 1.0)
    systRatioAll.SetPointEYhigh(i, datay/mcy - datay/(mcy+yErrUp))
    systRatioAll.SetPointEYlow(i, datay/(mcy-yErrDown) - datay/mcy)

totalRatioAll.Divide(run2Dictionary["DATA"]["h"],totalBkgdAll)    

XMin = 0;    XMax = 1; RatioXMin = 0; RatioXMax = 1 
YMin = 0.30; YMax = 1; RatioYMin = 0; RatioYMax = 0.30
PadFactor = (YMax-YMin) / (RatioYMax-RatioYMin)
cAll.Divide(1,2); cAll.cd(1); ROOT.gPad.SetLogy(); ROOT.gPad.SetLogz(); ROOT.gPad.SetPad(XMin, YMin, XMax, YMax)

ROOT.gPad.SetTopMargin(0.10)
ROOT.gPad.SetLeftMargin(0.15)
ROOT.gPad.SetBottomMargin(0.00)
ROOT.gPad.SetRightMargin(0.04)

iamLegend = ROOT.TLegend(ROOT.gPad.GetLeftMargin() + 0.02, 0.85 - ROOT.gPad.GetTopMargin(), ROOT.gPad.GetLeftMargin() + 0.8, 1 - ROOT.gPad.GetTopMargin() - 0.02)

iamLegend.SetTextSize(0.040)
iamLegend.SetNColumns(5)
iamLegend.SetBorderSize(2)

sigLegend = ROOT.TLegend(ROOT.gPad.GetLeftMargin() + 0.02, 0.75 - ROOT.gPad.GetTopMargin(), 1. - ROOT.gPad.GetRightMargin() - 0.02, 0.85 - ROOT.gPad.GetTopMargin())
sigLegend.SetTextSize(0.040)
sigLegend.SetNColumns(2)
sigLegend.SetBorderSize(2)

for process in ["TT", "TTX", "OTHER", "QCD"]: iamLegend.AddEntry(run2Dictionary[process]["h"], aliases[process], "F") 
for process in ["SIG1", "SIG3", "SIG2", "SIG4"]: sigLegend.AddEntry(run2Dictionary[process]["h"], aliases[process], "L")

iamLegend.AddEntry(run2Dictionary["DATA"]["h"], "Data", "EP")

dummy1 = ROOT.TH1F("dummy1", "dummy1", 6, 0, 6)
dummy1.SetTitle("")
dummy1.GetYaxis().SetTitle("Events / bin")
dummy1.GetYaxis().SetTitleSize(0.07)
dummy1.GetXaxis().SetTitleSize(0.07)
dummy1.GetYaxis().SetLabelSize(0.06)
dummy1.GetXaxis().SetLabelSize(0.06)
dummy1.GetYaxis().SetTitleOffset(1.1)
dummy1.GetXaxis().SetTitleOffset(2.1)
dummy1.SetMaximum(2e8)
dummy1.SetMinimum(5)

dummy1.Draw()
theStackAll.Draw("SAME")
run2Dictionary["DATA"]["h"].Draw("E0P SAME")
run2Dictionary["SYST"]["graph"].Draw("2SAME")
run2Dictionary["SIG1"]["h"].Draw("SAME")
run2Dictionary["SIG2"]["h"].Draw("SAME")
run2Dictionary["SIG3"]["h"].Draw("SAME")
run2Dictionary["SIG4"]["h"].Draw("SAME")

iamLegend.Draw("SAME")
sigLegend.Draw("SAME")

mark = ROOT.TLatex()
mark.SetNDC(True)

mark.SetTextAlign(11);
mark.SetTextSize(0.065);
mark.SetTextFont(61);
mark.DrawLatex(ROOT.gPad.GetLeftMargin(), 1 - (ROOT.gPad.GetTopMargin() - 0.02), "CMS")
mark.SetTextFont(52);
if not args.approved: mark.DrawLatex(ROOT.gPad.GetLeftMargin() + 0.103, 1 - (ROOT.gPad.GetTopMargin() - 0.02), "Preliminary")

mark.SetTextFont(42)
mark.SetTextAlign(31)
mark.DrawLatex(1 - ROOT.gPad.GetRightMargin(), 1 - (ROOT.gPad.GetTopMargin() - 0.02), "137 fb^{-1} (13 TeV)")

cAll.cd(2)

ROOT.gPad.SetGridy()
ROOT.gPad.SetTopMargin(0.00)
ROOT.gPad.SetBottomMargin(0.35)
ROOT.gPad.SetRightMargin(0.04)
ROOT.gPad.SetLeftMargin(0.15)
ROOT.gPad.SetPad(RatioXMin, RatioYMin, RatioXMax, RatioYMax)

for i in xrange(1, totalRatioAll.GetNbinsX()+1):
    if i != totalRatioAll.GetNbinsX(): totalRatioAll.GetXaxis().SetBinLabel(i, str(i+6))
    else:                           totalRatioAll.GetXaxis().SetBinLabel(i, "#geq %s"%(str(i+6)))

totalRatioAll.GetYaxis().SetRangeUser(0.91,1.09)
totalRatioAll.SetTitle("")
totalRatioAll.GetXaxis().SetTitle("Number of jets")
totalRatioAll.GetYaxis().SetNdivisions(5, 5, 0)
totalRatioAll.GetYaxis().SetTitle("Data / Pred.")
totalRatioAll.GetYaxis().SetTitleSize(dummy1.GetYaxis().GetTitleSize()*PadFactor)
totalRatioAll.GetXaxis().SetTitleSize(dummy1.GetXaxis().GetTitleSize()*PadFactor)
totalRatioAll.GetYaxis().SetLabelSize(dummy1.GetYaxis().GetLabelSize()*PadFactor)
totalRatioAll.GetXaxis().SetLabelSize(1.5*dummy1.GetXaxis().GetLabelSize()*PadFactor)
totalRatioAll.GetYaxis().SetTitleOffset(1.0*dummy1.GetYaxis().GetTitleOffset()/PadFactor)
totalRatioAll.GetXaxis().SetTitleOffset(dummy1.GetXaxis().GetTitleOffset()/PadFactor)
totalRatioAll.GetXaxis().SetLabelOffset(0.02)
totalRatioAll.SetMarkerStyle(20); totalRatioAll.SetMarkerSize(5); totalRatioAll.SetMarkerColor(ROOT.kBlack);
totalRatioAll.SetLineWidth(5); totalRatioAll.SetLineColor(ROOT.kBlack)

totalRatioAll.Draw("E0P")
systRatioAll.Draw("2SAME")

if args.approved:
    cAll.SaveAs("%s/njets_DataVsMC_%s.pdf"%(outpath,fitStr1))
else:
    cAll.SaveAs("%s/njets_DataVsMC_%s_prelim.pdf"%(outpath,fitStr1))
