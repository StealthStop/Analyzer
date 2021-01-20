#! /bin/env/python

# Makes ttbar vs signal plots of NN inputs for 2016 and 2017 for the legacy 1L analysis
# These plots can be found in Fig. 7 and 8 of the in the analysis supplementary plots twiki
# https://twiki.cern.ch/twiki/bin/view/CMS/SUS19004SupplementaryMaterial

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
parser.add_argument("--approved", dest="approved", help="Plot is approved",      action="store_true", default=False) 
parser.add_argument("--year",     dest="year",     help="which year",            required=True)
parser.add_argument("--inputDir", dest="inputDir", help="Input dir with histos", required=True)
args = parser.parse_args()

inputDir = args.inputDir

outpath = "./PlotsForSupplementary/"
if not os.path.exists(outpath): os.makedirs(outpath)

year = args.year

ttf = ROOT.TFile.Open(inputDir + "/%s/%s_TT.root"%(year,year))
sig1f = ROOT.TFile.Open(inputDir + "/%s/%s_RPV_2t6j_mStop-350.root"%(year,year))
sig2f = ROOT.TFile.Open(inputDir + "/%s/%s_StealthSYY_2t6j_mStop-750.root"%(year,year))

histDict = {"Jet_cm_pt_1_1l_ge7j_ge1b" : {"X" : {"title" : "Leading Jet p_{T} [GeV]", "rebin" : 2}},
           "Jet_cm_eta_1_1l_ge7j_ge1b" : {"X" : {"title" : "Leading Jet #eta", "rebin" : 1}},
           "Jet_cm_phi_1_1l_ge7j_ge1b" : {"X" : {"title" : "Leading Jet #phi", "rebin" : 1}},
           "Jet_cm_m_1_1l_ge7j_ge1b"   : {"X" : {"title" : "Leading Jet mass [GeV]", "rebin" : 10}},
           "jmt_ev0_top6_1l_ge7j_ge1b" : {"X" : {"title" : "JMT0", "rebin" : 1}},
           "jmt_ev1_top6_1l_ge7j_ge1b" : {"X" : {"title" : "JMT1", "rebin" : 1}},
           "fwm2_top6_1l_ge7j_ge1b"    : {"X" : {"title" : "FWM2", "rebin" : 1}},
           "fwm4_top6_1l_ge7j_ge1b"    : {"X" : {"title" : "FWM4", "rebin" : 1}}
}

for name, options in histDict.iteritems():

    c = ROOT.TCanvas("%s_c"%(name), "%s_c"%(name), 2400, 2400)
    tth = ttf.Get(name)
    sig1h = sig1f.Get(name)
    sig2h = sig2f.Get(name)

    tth.GetXaxis().SetTitle(options["X"]["title"])

    tth.Rebin(options["X"]["rebin"])
    sig1h.Rebin(options["X"]["rebin"])
    sig2h.Rebin(options["X"]["rebin"])

    sig1col = "#D325D3"
    sig2col = "#FFA851"
    ttcol   = "#9999FF"
    ttcol2  = "#010199"

    tth.Scale(1.0/tth.Integral())
    sig1h.Scale(1.0/sig1h.Integral())
    sig2h.Scale(1.0/sig2h.Integral())

    maxBinTT = tth.GetMaximumBin();   theMaxTT = tth.GetBinContent(maxBinTT)
    maxBinSG = sig2h.GetMaximumBin(); theMaxSG = sig2h.GetBinContent(maxBinSG)

    if theMaxTT > theMaxSG: tth.SetMaximum(theMaxTT*1.5)
    else:                   tth.SetMaximum(theMaxSG*1.5)

    for xbin in xrange(1, tth.GetNbinsX()+1): tth.SetBinError(xbin, 0.00001)
    for xbin in xrange(1, sig1h.GetNbinsX()+1): sig1h.SetBinError(xbin, 0.00001)
    for xbin in xrange(1, sig2h.GetNbinsX()+1): sig2h.SetBinError(xbin, 0.00001)

    tth.SetFillColorAlpha(ROOT.TColor.GetColor(ttcol), 1.0)
    sig1h.SetFillColor(ROOT.TColor.GetColor(sig1col))
    sig2h.SetFillColor(ROOT.TColor.GetColor(sig2col))

    tth.SetLineColor(ROOT.TColor.GetColor(ttcol2))
    sig1h.SetLineColor(ROOT.TColor.GetColor(sig1col))
    sig2h.SetLineColor(ROOT.TColor.GetColor(sig2col))

    tth.SetMarkerColor(ROOT.TColor.GetColor(ttcol))
    sig1h.SetMarkerColor(ROOT.TColor.GetColor(sig1col))
    sig2h.SetMarkerColor(ROOT.TColor.GetColor(sig2col))

    lw = 2
    tth.SetLineWidth(lw)
    sig1h.SetLineWidth(lw)
    sig2h.SetLineWidth(lw)

    ms = 0
    tth.SetMarkerSize(ms)
    sig1h.SetMarkerSize(ms)
    sig2h.SetMarkerSize(ms)

    sig1h.SetFillStyle(3004)
    sig2h.SetFillStyle(3005)

    ROOT.gPad.SetLogz()
    
    ROOT.gPad.SetTopMargin(0.10)
    ROOT.gPad.SetLeftMargin(0.14)
    ROOT.gPad.SetBottomMargin(0.12)
    ROOT.gPad.SetRightMargin(0.04)
    
    iamLegend = ROOT.TLegend(ROOT.gPad.GetLeftMargin() + 0.05, 0.70, 1 - ROOT.gPad.GetRightMargin() - 0.05, 0.8)
    
    iamLegend.SetTextSize(0.03)
    iamLegend.SetNColumns(3)
    iamLegend.SetBorderSize(2)
    
    iamLegend.AddEntry(sig1h, "RPV m_{#tilde{t}} = 350 GeV", "F")
    iamLegend.AddEntry(sig2h, "SYY m_{#tilde{t}} = 750 GeV", "F")
    iamLegend.AddEntry(tth,   "t#bar{t}", "F")

    tth.SetTitle("")
    tth.GetYaxis().SetTitle("A.U.")
    tth.GetYaxis().SetTitleSize(0.05)
    tth.GetXaxis().SetTitleSize(0.05)
    tth.GetYaxis().SetLabelSize(0.04)
    tth.GetXaxis().SetLabelSize(0.04)
    tth.GetYaxis().SetTitleOffset(1.4)
    tth.GetXaxis().SetTitleOffset(1.1)

    tth.Draw("HIST")
    sig1h.Draw("HIST SAME")
    sig2h.Draw("HIST SAME")
    iamLegend.Draw("SAME")

    mark = ROOT.TLatex()
    mark.SetNDC(True)

    mark.SetTextAlign(11);
    mark.SetTextSize(0.045);
    mark.SetTextFont(61);
    mark.DrawLatex(ROOT.gPad.GetLeftMargin() + 0.03, 1 - (ROOT.gPad.GetTopMargin() + 0.05), "CMS")
    mark.SetTextFont(52);
    if args.approved:
        mark.DrawLatex(ROOT.gPad.GetLeftMargin() + 0.13, 1 - (ROOT.gPad.GetTopMargin() + 0.05), "Simulation Supplementary")
    else:
        mark.DrawLatex(ROOT.gPad.GetLeftMargin() + 0.13, 1 - (ROOT.gPad.GetTopMargin() + 0.05), "Simulation Preliminary")

    mark.SetTextAlign(11)
    mark.SetTextSize(0.025)
    mark.SetTextFont(42)
    mark.DrawLatex(ROOT.gPad.GetLeftMargin() + 0.63, 1 - (ROOT.gPad.GetTopMargin() + 0.05), "arXiv:XXXX.XXXX")

    mark.SetTextSize(0.045);
    mark.SetTextFont(42)
    mark.SetTextAlign(31)
    if year == "2016":
        mark.DrawLatex(1 - ROOT.gPad.GetRightMargin(), 1 - (ROOT.gPad.GetTopMargin() - 0.015), "2016 (13 TeV)")
    else:
        mark.DrawLatex(1 - ROOT.gPad.GetRightMargin(), 1 - (ROOT.gPad.GetTopMargin() - 0.015), "2017 (13 TeV)")
  
    if args.approved:
        c.SaveAs("%s/%s_%s.pdf"%(outpath,year,name.replace("_1l_ge7j_ge1b", "")))
        c.SaveAs("%s/%s_%s.png"%(outpath,year,name.replace("_1l_ge7j_ge1b", "")))
    else:
        c.SaveAs("%s/%s_%s_prelim.pdf"%(outpath,year,name.replace("_1l_ge7j_ge1b", "")))
        c.SaveAs("%s/%s_%s_prelim.png"%(outpath,year,name.replace("_1l_ge7j_ge1b", "")))
