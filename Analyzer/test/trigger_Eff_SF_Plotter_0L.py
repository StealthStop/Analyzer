import ROOT
import math
import sys
import argparse
import numpy as np
import array as arr
from ROOT import TFile, gROOT, gStyle

ROOT.gStyle.SetPaintTextFormat("3.2f")

debug = True

def print_db(input):
    if (debug):
        print input


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


def makeEfficiency(dataNum, dataDen, mcNum, mcDen, year, dataset, mc, var, varKey):

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
    Efficiency_data.Divide(dataNum, dataDen, 1, 1, "B")
    Efficiency_data.SetLineWidth(3)
    Efficiency_data.SetLineColor(1)
    Efficiency_data.SetMarkerColor(1)
    Efficiency_data.SetMarkerSize(3.0)
    Efficiency_data.SetMarkerStyle(20)
    Efficiency_data.GetYaxis().SetRangeUser(0.65, 1.1)
    # get mc efficiency
    Efficiency_mc = mcNum.Clone()
    Efficiency_mc.Divide(mcNum, mcDen, 1, 1, "B")
    Efficiency_mc.SetLineWidth(3)
    Efficiency_mc.SetLineColor(2)
    Efficiency_mc.SetMarkerColor(2)
    Efficiency_mc.SetMarkerSize(3.0)
    Efficiency_mc.SetMarkerStyle(20)
    # calculate scale factor
    ScaleFactor = dataNum.Clone()
    ScaleFactor.Divide(Efficiency_data, Efficiency_mc)
    # fill eff. and sf to canvas
    canvas = ROOT.TCanvas("c_%s"%(var), "c_%s"%(var), 1400, 1400)
    legend = ROOT.TLegend(0.58, 0.20, 0.94, 0.42, "", "trNDC")
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
    ScaleFactor = dataNum.Clone()
    ScaleFactor.Divide(Efficiency_data,Efficiency_mc)
    ScaleFactor.SetMinimum(0.8) 
    ScaleFactor.SetMaximum(1.2) 
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
    ScaleFactor.Draw("E1 P")
    canvas.SaveAs("triggerEfficiencySF_plots/triggerEffSFs_0l/" + year + dataset + "_" + mc + "_" + var + "_TriggerEfficiency" + ".pdf")


def make2DScaleFactor(dataNum, dataDen, mcNum, mcDen, year, dataset, outputFile):

    XMin = 0;    XMax = 1; RatioXMin = 0; RatioXMax = 1
    YMin = 0.20; YMax = 1; RatioYMin = 0; RatioYMax = 0.20
    PadFactor  = (YMax-YMin) / (RatioYMax-RatioYMin)
    XTitleSize = 0.045;  XLabelSize = 0.036;  XTitleOffset = 1.0 
    YTitleSize = 0.045;  YLabelSize = 0.036;  YTitleOffset = 0.9

    Efficiency_data = dataNum.Clone() 
    Efficiency_mc   = mcNum.Clone() 
    ScaleFactor     = dataNum.Clone()
    Efficiency_data.Divide(dataNum, dataDen, 1, 1, "B")
    Efficiency_mc.Divide(mcNum, mcDen, 1, 1, "B")
    ScaleFactor.Divide(Efficiency_data, Efficiency_mc)
    ScaleFactor.SetTitle("")

    #Efficiency_data.GetZaxis().SetRangeUser(0.75,1.1)
    #Efficiency_data.SetTitle("")
    #Efficiency_data.GetYaxis().SetTitle(6thJetPt [GeV)
    #Efficiency_data.GetXaxis().SetTitle("HT [GeV]")
    #Efficiency_data.GetYaxis().SetTitleSize(YTitleSize)
    #Efficiency_data.GetXaxis().SetTitleSize(XTitleSize)
    #Efficiency_data.GetYaxis().SetLabelSize(YLabelSize)
    #Efficiency_data.GetXaxis().SetLabelSize(XLabelSize)
    #Efficiency_data.GetYaxis().SetTitleOffset(YTitleOffset)
    #Efficiency_data.GetXaxis().SetTitleOffset(XTitleOffset)
    #Efficiency_data.SetContour(255)

    #Efficiency_mc.GetZaxis().SetRangeUser(0.75,1.1)
    #Efficiency_mc.SetTitle("")
    #Efficiency_mc.GetYaxis().SetTitle(6thJetPt [GeV)
    #Efficiency_mc.GetXaxis().SetTitle("HT [GeV]")
    #Efficiency_mc.GetYaxis().SetTitleSize(YTitleSize)
    #Efficiency_mc.GetXaxis().SetTitleSize(XTitleSize)
    #Efficiency_mc.GetYaxis().SetLabelSize(YLabelSize)
    #Efficiency_mc.GetXaxis().SetLabelSize(XLabelSize)
    #Efficiency_mc.GetYaxis().SetTitleOffset(YTitleOffset)
    #Efficiency_mc.GetXaxis().SetTitleOffset(XTitleOffset)
    #Efficiency_mc.SetContour(255)

    ScaleFactor.GetZaxis().SetRangeUser(0.65, 1.05) # 0.75, 1.1 / 0.5, 1.2
    #ScaleFactor.GetZaxis().SetRangeUser(0.0, 1.025) #0.65, 1.025
    ScaleFactor.SetTitle("")
    ScaleFactor.GetYaxis().SetTitle("6^{th} Jet p_{T} [GeV]")
    ScaleFactor.GetXaxis().SetTitle("H_{T} [GeV]")
    ScaleFactor.GetYaxis().SetTitleSize(YTitleSize)
    ScaleFactor.GetXaxis().SetTitleSize(XTitleSize)
    ScaleFactor.GetYaxis().SetLabelSize(YLabelSize)
    ScaleFactor.GetXaxis().SetLabelSize(XLabelSize)
    ScaleFactor.GetYaxis().SetTitleOffset(YTitleOffset)
    ScaleFactor.GetXaxis().SetTitleOffset(XTitleOffset)
    #ScaleFactor.GetXaxis().SetRangeUser(500,1500)
    #ScaleFactor.GetYaxis().SetRangeUser(45,120) # (40,120)
    ScaleFactor.SetContour(255)

    #setMinimumErrors(ScaleFactor)

    theName = dataNum.GetName().replace("numerator", year) + dataset

    #aCanvas = ROOT.TCanvas("c_data_%s"%(theName), "c_data_%s"%(theName), 1500, 1200)
    #ROOT.gPad.Clear()
    #aCanvas.cd()
    #ROOT.gPad.SetTopMargin(0.02)
    #ROOT.gPad.SetLeftMargin(0.09)
    #ROOT.gPad.SetBottomMargin(0.11)
    #ROOT.gPad.SetRightMargin(0.11)
    #ROOT.gPad.SetTicks()
    #Efficiency_data.Draw("COLZ E TEXT")
    #aCanvas.SaveAs("plots/plots_refAN_18.03.2021/%s.pdf"%("2016_data_"+theName))

    #aCanvas = ROOT.TCanvas("c_mc_%s"%(theName), "c_mc_%s"%(theName), 1500, 1200)
    #ROOT.gPad.Clear()
    #aCanvas.cd()
    #ROOT.gPad.SetTopMargin(0.02)
    #ROOT.gPad.SetLeftMargin(0.09)
    #ROOT.gPad.SetBottomMargin(0.11)
    #ROOT.gPad.SetRightMargin(0.11)
    #ROOT.gPad.SetTicks()
    #Efficiency_mc.Draw("COLZ E TEXT")
    #aCanvas.SaveAs("plots/plots_refAN_18.03.2021/%s.pdf"%("2016_mc_"+theName))

    aCanvas = ROOT.TCanvas("c_ScaleFactor_%s"%(theName), "c_ScaleFactor_%s"%(theName), 1400, 1400)
    ROOT.gPad.Clear()
    aCanvas.cd()
    ROOT.gPad.SetTopMargin(0.098)
    ROOT.gPad.SetBottomMargin(0.11)
    ROOT.gPad.SetLeftMargin(0.09)
    ROOT.gPad.SetRightMargin(0.11)
    ROOT.gPad.SetTicks()
    ScaleFactor.Draw("COLZ E TEXT")
    addCMSlogo(aCanvas, year, TopMargin=0.12, LeftMargin=0.09, RightMargin=0.11, SF=1.0)
    aCanvas.SaveAs("triggerEfficiencySF_plots/triggerEffSFS_0l/%s.pdf"%(theName+"_TriggerSF"))

    # make the zoom version of the SF plots
    #ScaleFactor.GetXaxis().SetRangeUser(500, 1000)
    #ScaleFactor.GetYaxis().SetRangeUser(45, 90)
    #ScaleFactor.GetZaxis().SetRangeUser(0.65, 1.05)
    #aCanvas.SaveAs("test/latest/%s.pdf"%(theName+"_TriggerSF"+"_zoom"))

    # write SF histograms to the root file
    if outputFile is not None:
        outputFile.cd()
        ScaleFactor.SetName(theName)
        ScaleFactor.Write(theName) 

def main():
    gROOT.SetBatch(True)
    gStyle.SetOptStat(0)
    #gStyle.SetOptStat(1)    
 
    # ---------------------------------
    # root path & years & histograms
    # --------------------------------- 
    path = "/uscms_data/d3/semrat/SUSY/CMSSW_11_2_0_pre5/src/Analyzer/Analyzer/test/condor/Thesis_AN_2022/2_triggerSFs/hadd_year_HadronicTriggerEfficiencySF_12.10.2022/"

    years = [
        "2016preVFP" ,
        "2016postVFP" ,
        "2017" ,
        "2018" ,
    ]

    varList1D = [
        "HT",
        "6thJetPt",
        #"NJet",
        #"NBJet",
    ]

    varKeys = {
        "HT"        : "H_{T}",
        "6thJetPt"  : "6^{th} Jet p_{T}",
        #"NJet"      : "N_{jet}",
        #"NBJet"     : "N_{bjet}",
    }

    varList2D = [
        #"ge2bjetCut_pt45_HTvs6thJetPt",
        "2bjetCut_pt45_HTvs6thJetPt",
        "3bjetCut_pt45_HTvs6thJetPt",
        "ge4bjetCut_pt45_HTvs6thJetPt",

        #"ge6jetCut_pt45_HTvs6thJetPt",
        #"6jetCut_pt45_HTvs6thJetPt",
        #"7jetCut_pt45_HTvs6thJetPt",
        #"8jetCut_pt45_HTvs6thJetPt",
        #"9jetCut_pt45_HTvs6thJetPt",
        #"ge10jetCut_pt45_HTvs6thJetPt",
    ]
   
    # ------------------------------------------------------------------------------
    # commandline options
    #   -- python triggerRefAN_Efficiency_SF_0l.py --model RPV --mass 550 
    # ------------------------------------------------------------------------------
    usage  = "usage: %prog [options]"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("--model", dest="model", help="signal model", default="RPV")
    parser.add_argument("--mass",  dest="mass",  help="signal mass",  default="550")
    args = parser.parse_args()

    # -----------------------------
    # loop over & get histograms 
    # -----------------------------
    for year in years:
        print_db("Processing year " + year)

        filename1 = path.replace("year",year) + year + "_" + args.model + "_2t6j_mStop-" + args.mass + ".root"  
        print_db("Opening file " + filename1)
        f1 = ROOT.TFile.Open(filename1, "READ")

        filename2 = path.replace("year",year) + year + "_" + "TT" + ".root"
        print_db("Opening file " + filename2)
        f2 = ROOT.TFile.Open(filename2, "READ")

        filename3 = path.replace("year",year) + year + "_" + "Data_SingleMuon" + ".root"
        print_db("Opening file " + filename3)
        f3 = ROOT.TFile.Open(filename3, "READ")

        f4 = ROOT.TFile.Open("triggerEfficiencySF_plots/triggerEffSFS_0l/%s_Hadronic_Triggers_SF.root"%(year), "RECREATE")


        for var in varList1D:
            ROOT.TH1.AddDirectory(0)                                      

            h_SingleMuon_Denominator = f3.Get("h_denominator_CombHadIsoMu_noTrig_ge2bjetCut_pt45_%s"%(var))
            h_SingleMuon_Numerator   = f3.Get("h_numerator_CombHadIsoMu_trig_ge2bjetCut_pt45_%s"%(var))
            h_RPV550_Denominator     = f1.Get("h_denominator_CombHadIsoMu_noTrig_ge2bjetCut_pt45_%s"%(var))
            h_RPV550_Numerator       = f1.Get("h_numerator_CombHadIsoMu_trig_ge2bjetCut_pt45_%s"%(var))
            h_TT_Denominator         = f2.Get("h_denominator_CombHadIsoMu_noTrig_ge2bjetCut_pt45_%s"%(var))
            h_TT_Numerator           = f2.Get("h_numerator_CombHadIsoMu_trig_ge2bjetCut_pt45_%s"%(var))
            makeEfficiency(h_SingleMuon_Numerator, h_SingleMuon_Denominator, h_RPV550_Numerator, h_RPV550_Denominator, year, "_SingleMuon", args.model + "_" + args.mass, var, varKeys[var])            
            makeEfficiency(h_SingleMuon_Numerator, h_SingleMuon_Denominator, h_TT_Numerator, h_TT_Denominator, year, "_SingleMuon", "TT", var, varKeys[var])

        for var in varList2D:
            ROOT.TH1.AddDirectory(0)

            h_SingleMuon_Denominator = f3.Get("h_denominator_CombHadIsoMu_noTrig_%s"%(var))
            h_SingleMuon_Numerator   = f3.Get("h_numerator_CombHadIsoMu_trig_%s"%(var))
            h_RPV550_Denominator     = f1.Get("h_denominator_CombHadIsoMu_noTrig_%s"%(var))
            h_RPV550_Numerator       = f1.Get("h_numerator_CombHadIsoMu_trig_%s"%(var))
            h_TT_Denominator         = f2.Get("h_denominator_CombHadIsoMu_noTrig_%s"%(var))
            h_TT_Numerator           = f2.Get("h_numerator_CombHadIsoMu_trig_%s"%(var))
            make2DScaleFactor(h_SingleMuon_Numerator, h_SingleMuon_Denominator, h_RPV550_Numerator, h_RPV550_Denominator, year, "_SingleMuon" + "_" + args.model + "_" + args.mass, None)
            make2DScaleFactor(h_SingleMuon_Numerator, h_SingleMuon_Denominator, h_TT_Numerator, h_TT_Denominator, year, "_SingleMuon_TT", f4)

    f1.Close()
    f2.Close()
    f3.Close()
    f4.Close()
 
if __name__ == "__main__":
    main()

