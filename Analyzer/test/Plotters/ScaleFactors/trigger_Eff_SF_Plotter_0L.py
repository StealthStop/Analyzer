import ROOT
import math
import sys
import argparse
import numpy as np
import array as arr
from ROOT import TFile, gROOT, gStyle

ROOT.gStyle.SetPaintTextFormat("3.3f")

debug = True

def print_db(input):
    if (debug):
        print input


# -----------------
# set minimum error
# ----------------- 
def setMinimumErrors(dataMcRatio):

    xbins = dataMcRatio.GetXaxis().GetNbins()+1
    ybins = dataMcRatio.GetYaxis().GetNbins()+1

    for xbin in xrange(xbins):
        for ybin in xrange(ybins):

            content = dataMcRatio.GetBinContent(xbin,ybin)
            error   = dataMcRatio.GetBinError(xbin,ybin)

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
def make_Efficiency_Plots(dataNum, dataDen, mcNum, mcDen, year, dataset, mc, var, varKey, bjet_cut):

    XMin = 0;    XMax = 1; RatioXMin = 0; RatioXMax = 1
    YMin = 0.30; YMax = 1; RatioYMin = 0; RatioYMax = 0.30
    PadFactor  = (YMax-YMin) / (RatioYMax-RatioYMin)
    XTitleSize = 0.06;  XLabelSize = 0.05;  XTitleOffset = 1.0
    YTitleSize = 0.06;  YLabelSize = 0.05;  YTitleOffset = 0.6
    # get data efficiency
    Efficiency_data = dataNum.Clone()
    Efficiency_data.SetTitle("")
    Efficiency_data.GetYaxis().SetTitle("L1+HLT Efficiency")
    Efficiency_data.GetXaxis().SetLabelSize(0)
    Efficiency_data.GetYaxis().SetTitleSize(0.05)
    Efficiency_data.GetYaxis().SetLabelSize(0.05)
    Efficiency_data.SetLineWidth(3)
    Efficiency_data.SetLineColor(1)
    Efficiency_data.SetMarkerColor(1)
    Efficiency_data.SetMarkerSize(3.0)
    Efficiency_data.SetMarkerStyle(20)
    Efficiency_data.Divide(dataNum, dataDen, 1, 1, "B")
    Efficiency_data.GetYaxis().SetRangeUser(0.4, 1.1)
    # get mc efficiency
    Efficiency_mc = mcNum.Clone()
    Efficiency_mc.SetLineWidth(3)
    Efficiency_mc.SetLineColor(2)
    Efficiency_mc.SetMarkerColor(2)
    Efficiency_mc.SetMarkerSize(3.0)
    Efficiency_mc.SetMarkerStyle(20)
    Efficiency_mc.Divide(mcNum, mcDen, 1, 1, "B")
    Efficiency_mc.GetYaxis().SetRangeUser(0.4, 1.1)
    # get errors
    updateErrors(Efficiency_data, dataDen)
    updateErrors(Efficiency_mc, mcDen)
    # get scale factor
    ScaleFactor = dataNum.Clone()
    ScaleFactor.Divide(Efficiency_data, Efficiency_mc)
    ScaleFactor.SetMinimum(0.7) # 0.8 
    ScaleFactor.SetMaximum(1.1) # 1.2
    ScaleFactor.GetYaxis().SetNdivisions(-304) 
    ScaleFactor.SetTitle("")
    ScaleFactor.SetLineWidth(3)
    ScaleFactor.SetMarkerSize(3.0)
    ScaleFactor.SetMarkerStyle(20)
    ScaleFactor.SetMarkerColor(ScaleFactor.GetLineColor())
    ScaleFactor.GetXaxis().SetTitle(varKey)
    ScaleFactor.GetXaxis().SetTitleSize(PadFactor*XTitleSize)
    ScaleFactor.GetXaxis().SetLabelSize(PadFactor*XLabelSize)   
    ScaleFactor.GetXaxis().SetTitleOffset(1.9*XTitleOffset/PadFactor)
    ScaleFactor.GetYaxis().SetTitle("Data / MC")
    ScaleFactor.GetYaxis().SetTitleSize(0.8*PadFactor*YTitleSize)
    ScaleFactor.GetYaxis().SetLabelSize(PadFactor*YLabelSize)
    ScaleFactor.GetYaxis().SetTitleOffset(1.6*YTitleOffset/PadFactor)
    
    # fill Eficiency and SF to canvas
    canvas = ROOT.TCanvas("c_%s"%(var), "c_%s"%(var), 1400, 1400)
    legend = ROOT.TLegend(0.63, 0.05, 0.99, 0.27, "", "trNDC")
    canvas.Divide(1,2)
    canvas.cd(1)
    ROOT.gPad.SetPad(XMin, YMin, XMax, YMax)
    ROOT.gPad.SetGrid()
    ROOT.gPad.SetTopMargin(0.1)
    ROOT.gPad.SetBottomMargin(0.017)
    ROOT.gPad.SetLeftMargin(0.17)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetTicks()
    Efficiency_data.Draw("E1 P SAME")
    Efficiency_mc.Draw("E1 P SAME")
    legend.SetTextSize(0.035)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(Efficiency_data, "Data" + dataset, "E1 P")
    legend.AddEntry(Efficiency_mc, mc, "E1 P")
    legend.Draw("SAME")
    addCMSlogo(canvas, year, TopMargin=0.098, LeftMargin=0.17, RightMargin=0.05, SF=1.0)
    canvas.cd(2)
    ROOT.gPad.SetPad(RatioXMin, RatioYMin, RatioXMax, RatioYMax)
    ROOT.gPad.SetTopMargin(0.05)
    ROOT.gPad.SetBottomMargin(0.35)
    ROOT.gPad.SetLeftMargin(0.17)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetTicks()
    ROOT.gPad.SetGridy()
    ScaleFactor.Draw("E1 P")
    canvas.SaveAs("triggerEfficiencySF_plots/triggerEffSFS_0l/" + year + "_jet_trig_{}_".format(bjet_cut) + var + "_TriggerEfficiency" + ".pdf")

# --------------------------
# Make 2D Scale Factor plots
# --------------------------
def make_ScaleFactor_Plots(dataNum, dataDen, mcNum, mcDen, year, outputFile, xTitle, yTitle):

    XMin = 0;    XMax = 1; RatioXMin = 0; RatioXMax = 1
    YMin = 0.20; YMax = 1; RatioYMin = 0; RatioYMax = 0.20
    PadFactor  = (YMax-YMin) / (RatioYMax-RatioYMin)
    XTitleSize = 0.045;  XLabelSize = 0.036;  XTitleOffset = 1.0 
    YTitleSize = 0.045;  YLabelSize = 0.036;  YTitleOffset = 0.9

    # get data and mc efficiency 
    Efficiency_data = dataNum.Clone() 
    Efficiency_data.Divide(dataNum, dataDen, 1, 1, "B")
    Efficiency_mc = mcNum.Clone() 
    Efficiency_mc.Divide(mcNum, mcDen, 1, 1, "B")
    # get errors
    updateErrors(Efficiency_data, dataDen)
    updateErrors(Efficiency_mc, mcDen)  
    # get scale factor 
    ScaleFactor = dataNum.Clone()
    ScaleFactor.Divide(Efficiency_data, Efficiency_mc)
    ScaleFactor.SetTitle("")
    ScaleFactor.GetZaxis().SetRangeUser(0.65, 1.05) 
    ScaleFactor.GetXaxis().SetTitle(xTitle)
    ScaleFactor.GetYaxis().SetTitle(yTitle)
    ScaleFactor.GetYaxis().SetTitleSize(YTitleSize)
    ScaleFactor.GetXaxis().SetTitleSize(XTitleSize)
    ScaleFactor.GetYaxis().SetLabelSize(YLabelSize)
    ScaleFactor.GetXaxis().SetLabelSize(XLabelSize)
    ScaleFactor.GetYaxis().SetTitleOffset(YTitleOffset)
    ScaleFactor.GetXaxis().SetTitleOffset(XTitleOffset)
    #ScaleFactor.GetXaxis().SetRangeUser(500,1500)
    #ScaleFactor.GetYaxis().SetRangeUser(45,120) # (40,120)
    ScaleFactor.SetContour(255)
 
    # fill SF to canvas
    theName = dataNum.GetName().replace("h2_num", year) 
    aCanvas = ROOT.TCanvas("c_ScaleFactor_%s"%(theName), "c_ScaleFactor_%s"%(theName), 1400, 1400)
    ROOT.gPad.Clear()
    aCanvas.cd()
    ROOT.gPad.SetTopMargin(0.098)
    ROOT.gPad.SetBottomMargin(0.11)
    ROOT.gPad.SetLeftMargin(0.09)
    ROOT.gPad.SetRightMargin(0.11)
    ROOT.gPad.SetTicks()
    ROOT.gPad.SetLogx()
    ScaleFactor.Draw("COLZ E TEXT")
    addCMSlogo(aCanvas, year, TopMargin=0.12, LeftMargin=0.09, RightMargin=0.11, SF=1.0)
    aCanvas.SaveAs("triggerEfficiencySF_plots/triggerEffSFS_0l/%s.pdf"%(theName+"_TriggerSF"))

    # write SF histograms to the root file
    if outputFile is not None:
        outputFile.cd()
        ScaleFactor.SetName(theName+"_TriggerSF")
        ScaleFactor.Write(theName+"_TriggerSF") 


# -------------
# Main function
# -------------
def main():

    # -------------------------------------------------------------------
    # commandline options
    #   -- python triggerRefAN_Efficiency_SF_0l.py --model RPV --mass 550
    # -------------------------------------------------------------------
    #usage  = "usage: %prog [options]"
    #parser = argparse.ArgumentParser(usage)
    #parser.add_argument("--inpath",  dest="inpath",  help="Path to root files", default="NULL", required=True      )
    #parser.add_argument("--outpath", dest="outpath", help="Where to put plots", default="NULL", required=True      )
    #parser.add_argument("--year",    dest="year",    help="which year",                         required=True      )
    #parser.add_argument("--channel", dest="channel", help="which channel",                      required=True      )
    #parser.add_argument("--wip",     dest="wip",     help="Work in Progress",   default=False,  action="store_true")
    #args = parser.parse_args()


    gROOT.SetBatch(True)
    gStyle.SetOptStat(0)
    #gStyle.SetOptStat(1)    
 
    # -------------------------------
    # root years & paths & histograms
    # ------------------------------- 
    path_jetTriggers    = "/uscms/home/bcrossma/nobackup/analysis/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/Thesis_AN_FullStatusTalk/hadd_Run2UL_HadronicTriggerEfficiencySF/"

    years = [
        "2016preVFP" ,
        "2016postVFP" ,
        "2017" ,
        "2018" ,
    ]
    
    varList_jetTrigEff = [
        "wJetHtBin",
        "w6thJetPtBin",
    ]

    varKeys_jetTrig = {
        "wJetHtBin"    : "H_{T}",
        "w6thJetPtBin" : "6^{th} Jet p_{T}",
    }

    varList_jetTrigSF = [
        "1bjetCut_wJetHt6thJetPtBin",
        "2bjetCut_wJetHt6thJetPtBin",
        "3bjetCut_wJetHt6thJetPtBin",
        "ge1bjetCut_wJetHt6thJetPtBin",
        "ge2bjetCut_wJetHt6thJetPtBin",
        "ge3bjetCut_wJetHt6thJetPtBin",
        "ge4bjetCut_wJetHt6thJetPtBin",
    ]

    # --------------------------
    # loop over & get histograms
    # --------------------------
    for year in years:
        print_db("Processing year " + year)

        # --------------------------------------
        # jet triggers root files and histograms
        # --------------------------------------
        filename1 = path_jetTriggers + year + "_" + "TT" + ".root"
        f1 = ROOT.TFile.Open(filename1, "READ")
        filename2 = path_jetTriggers + year + "_" + "Data_SingleMuon" + ".root"
        f2 = ROOT.TFile.Open(filename2, "READ")
        f3 = ROOT.TFile.Open("triggerEfficiencySF_plots/triggerEffSFS_0l/%s_Hadronic_Triggers_SF.root"%(year), "RECREATE")

        # --------------------------------
        # efficiency plots & scale factors 
        # efficiency = num / den
        # SF = Data Eff / MC Eff
        # --------------------------------
        for var in varList_jetTrigEff:
            ROOT.TH1.AddDirectory(0)                                      
            h_SingleMuon_Denominator = f2.Get("h_den_jet_trig_1bjetCut_%s"%(var))
            h_SingleMuon_Numerator   = f2.Get("h_num_jet_trig_1bjetCut_%s"%(var))
            h_TT_Denominator         = f1.Get("h_den_jet_trig_1bjetCut_%s"%(var))
            h_TT_Numerator           = f1.Get("h_num_jet_trig_1bjetCut_%s"%(var))
            make_Efficiency_Plots(h_SingleMuon_Numerator, h_SingleMuon_Denominator, h_TT_Numerator, h_TT_Denominator, year, "_SingleMuon", "TT", var, varKeys_jetTrig[var], "1bjetCut")

            h_SingleMuon_Denominator = f2.Get("h_den_jet_trig_ge2bjetCut_%s"%(var))
            h_SingleMuon_Numerator   = f2.Get("h_num_jet_trig_ge2bjetCut_%s"%(var))
            h_TT_Denominator         = f1.Get("h_den_jet_trig_ge2bjetCut_%s"%(var))
            h_TT_Numerator           = f1.Get("h_num_jet_trig_ge2bjetCut_%s"%(var))
            make_Efficiency_Plots(h_SingleMuon_Numerator, h_SingleMuon_Denominator, h_TT_Numerator, h_TT_Denominator, year, "_SingleMuon", "TT", var, varKeys_jetTrig[var], "ge2bjetCut")

        for var in varList_jetTrigSF:
            ROOT.TH1.AddDirectory(0)
            h_SingleMuon_Denominator = f2.Get("h2_den_jet_trig_%s"%(var))
            h_SingleMuon_Numerator   = f2.Get("h2_num_jet_trig_%s"%(var))
            h_TT_Denominator         = f1.Get("h2_den_jet_trig_%s"%(var))
            h_TT_Numerator           = f1.Get("h2_num_jet_trig_%s"%(var))
            make_ScaleFactor_Plots(h_SingleMuon_Numerator, h_SingleMuon_Denominator, h_TT_Numerator, h_TT_Denominator, year, f3, "H_{T} [GeV]", "6^{th} Jet p_{T} [GeV]")

    f1.Close()
    f2.Close()
    f3.Close()
 
if __name__ == "__main__":
    main()

