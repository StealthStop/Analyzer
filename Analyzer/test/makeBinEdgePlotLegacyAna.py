#!/bin/python
import ROOT
import math
import argparse

lb = ROOT.TColor.GetColor("#CDE4F1")
db = ROOT.TColor.GetColor("#98C8E3")
lg = ROOT.TColor.GetColor("#95C495")
dg = ROOT.TColor.GetColor("#499948")

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat("")

binNames = ["7", "8", "9", "10", "11", "#geq 12"]

usage = "usage: %prog [options]"
parser = argparse.ArgumentParser(usage)
parser.add_argument("--approved", dest="approved", help="Plot is approved", action="store_true", default=False) 

args = parser.parse_args()

def makeBinFill(e2, e1, name, color):

    g = ROOT.TGraphErrors(6); g.SetName(name)

    for edge in xrange(0, len(e1)):

        g.SetPoint(edge, float(edge+0.5), e2[edge] - (e2[edge]-e1[edge])/2.0)
        g.SetPointError(edge, 0.5, (e2[edge]-e1[edge])/2.0) 

    g.SetFillColorAlpha(color,1.0)
    g.SetLineWidth(0)
    g.SetMarkerSize(0)

    return g

def makeBinEdges(e, name):

    h = ROOT.TH1F(name, name, 6, 0, 6); h.SetDirectory(0)
    h.GetYaxis().SetRangeUser(0,1)

    for edge in xrange(1, 7):
        h.SetBinContent(edge, e[edge-1])
        h.SetBinError(edge, 0.0)

    h.GetYaxis().SetTitle("S_{NN}")
    h.GetXaxis().SetTitle("Number of Jets")

    for bin in xrange(1,h.GetNbinsX()+1): h.GetXaxis().SetBinLabel(bin, binNames[bin-1])
    
    h.GetXaxis().SetTitleSize(0.065)
    h.GetYaxis().SetTitleSize(0.065)

    h.GetXaxis().SetLabelSize(0.085)
    h.GetYaxis().SetLabelSize(0.055)

    h.GetYaxis().SetLabelOffset(0.015)

    h.GetXaxis().SetTitleOffset(0.8)
    h.GetYaxis().SetTitleOffset(1.05)

    h.SetTitle("")

    h.SetLineWidth(2)
    h.SetMarkerSize(0)

    h.SetLineColor(ROOT.kBlue+3)
    h.SetMarkerColor(ROOT.kBlue+3)

    return h

def main() :

    edges_0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    edges_1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

    edges16_1 = [0.349, 0.358, 0.371, 0.386, 0.410, 0.410]
    edges16_2 = [0.680, 0.740, 0.787, 0.816, 0.848, 0.853]
    edges16_3 = [0.835, 0.878, 0.901, 0.917, 0.936, 0.936]

    edges17_1 = [0.346, 0.363, 0.376, 0.384, 0.396, 0.396]
    edges17_2 = [0.714, 0.752, 0.780, 0.794, 0.801, 0.801]
    edges17_3 = [0.833, 0.864, 0.885, 0.895, 0.901, 0.901]

    edges18_1 = [0.345, 0.363, 0.374, 0.386, 0.386, 0.386]
    edges18_2 = [0.713, 0.752, 0.779, 0.805, 0.805, 0.805]
    edges18_3 = [0.832, 0.863, 0.883, 0.900, 0.911, 0.911]

    edges19_1 = [0.364, 0.386, 0.405, 0.411, 0.411, 0.411]
    edges19_2 = [0.729, 0.771, 0.798, 0.826, 0.826, 0.826]
    edges19_3 = [0.842, 0.874, 0.895, 0.914, 0.914, 0.914]

    h16_1 = makeBinEdges(edges16_1, "edges16_1")
    h16_2 = makeBinEdges(edges16_2, "edges16_2")
    h16_3 = makeBinEdges(edges16_3, "edges16_3")

    h17_1 = makeBinEdges(edges17_1, "edges17_1")
    h17_2 = makeBinEdges(edges17_2, "edges17_2")
    h17_3 = makeBinEdges(edges17_3, "edges17_3")

    h18_1 = makeBinEdges(edges18_1, "edges18_1")
    h18_2 = makeBinEdges(edges18_2, "edges18_2")
    h18_3 = makeBinEdges(edges18_3, "edges18_3")

    h19_1 = makeBinEdges(edges19_1, "edges19_1")
    h19_2 = makeBinEdges(edges19_2, "edges19_2")
    h19_3 = makeBinEdges(edges19_3, "edges19_3")


    f16_1 = makeBinFill(edges16_1, edges_0,   "fill16_1", db)
    f16_2 = makeBinFill(edges16_2, edges16_1, "fill16_2", lb)
    f16_3 = makeBinFill(edges16_3, edges16_2, "fill16_3", lg)
    f16_4 = makeBinFill(edges_1, edges16_3,   "fill16_4", dg)

    f17_1 = makeBinFill(edges17_1, edges_0,   "fill17_1", db)
    f17_2 = makeBinFill(edges17_2, edges17_1, "fill17_2", lb)
    f17_3 = makeBinFill(edges17_3, edges17_2, "fill17_3", lg)
    f17_4 = makeBinFill(edges_1, edges17_3,   "fill17_4", dg)

    f18_1 = makeBinFill(edges18_1, edges_0,   "fill18_1", db)
    f18_2 = makeBinFill(edges18_2, edges18_1, "fill18_2", lb)
    f18_3 = makeBinFill(edges18_3, edges18_2, "fill18_3", lg)
    f18_4 = makeBinFill(edges_1, edges18_3,   "fill18_4", dg)

    f19_1 = makeBinFill(edges19_1, edges_0,   "fill19_1", db)
    f19_2 = makeBinFill(edges19_2, edges19_1, "fill19_2", lb)
    f19_3 = makeBinFill(edges19_3, edges19_2, "fill19_3", lg)
    f19_4 = makeBinFill(edges_1, edges19_3,   "fill19_4", dg)

    line1 = ROOT.TLine(1.0,0.0,1.0,1.0); line1.SetLineWidth(1); line1.SetLineColor(ROOT.kBlack); line1.SetLineStyle(3)
    line2 = ROOT.TLine(2.0,0.0,2.0,1.0); line2.SetLineWidth(1); line2.SetLineColor(ROOT.kBlack); line2.SetLineStyle(3)
    line3 = ROOT.TLine(3.0,0.0,3.0,1.0); line3.SetLineWidth(1); line3.SetLineColor(ROOT.kBlack); line3.SetLineStyle(3)
    line4 = ROOT.TLine(4.0,0.0,4.0,1.0); line4.SetLineWidth(1); line4.SetLineColor(ROOT.kBlack); line4.SetLineStyle(3)
    line5 = ROOT.TLine(5.0,0.0,5.0,1.0); line5.SetLineWidth(1); line5.SetLineColor(ROOT.kBlack); line5.SetLineStyle(3)
    line6 = ROOT.TLine(6.0,0.0,6.0,1.0); line6.SetLineWidth(1); line6.SetLineColor(ROOT.kBlack); line6.SetLineStyle(3)

    c1 = ROOT.TCanvas( "c1", "c1", 3200, 800 )

    tm = 0.1; bm = 0.12; lm = 0.14; rm = 0.02

    l1 = 0.277967; l2 = 0.239052; l3 = 0.239052; l4 = 0.24393

    c1.Divide(4,1)
    c1.cd(1)

    ROOT.gPad.SetPad(0.0, 0.0, l1, 1.0)
    ROOT.gPad.SetGridx()
    ROOT.gPad.SetLeftMargin(lm)
    ROOT.gPad.SetRightMargin(0.0)
    ROOT.gPad.SetTopMargin(tm)
    ROOT.gPad.SetBottomMargin(bm)


    h16_1.Draw("HIST")
    f16_1.Draw("2SAME")
    f16_2.Draw("2SAME")
    f16_3.Draw("2SAME")
    f16_4.Draw("2SAME")
    h16_1.Draw("HIST SAME")
    h16_2.Draw("HIST SAME")
    h16_3.Draw("HIST SAME")

    line1.Draw("SAME")
    line2.Draw("SAME")
    line3.Draw("SAME")
    line4.Draw("SAME")
    line5.Draw("SAME")
    line6.Draw("SAME")

    mark = ROOT.TLatex();
    mark.SetNDC(True);

    mark.SetTextAlign(11)
    mark.SetTextSize(0.065)
    mark.SetTextFont(61)
    mark.DrawLatex(ROOT.gPad.GetLeftMargin(), 1 - (ROOT.gPad.GetTopMargin() - 0.020), "CMS")
    mark.SetTextAlign(31)
    mark.SetTextSize(0.060)
    mark.DrawLatex(1.0 - ROOT.gPad.GetRightMargin() - 0.01, 1 - (ROOT.gPad.GetTopMargin() - 0.020), "2016")
    mark.SetTextFont(52)
    mark.SetTextAlign(11)
    if args.approved:
        mark.DrawLatex(ROOT.gPad.GetLeftMargin() + 0.135, 1 - (ROOT.gPad.GetTopMargin() - 0.020), "Simulation Supplementary");
    else:
        mark.DrawLatex(ROOT.gPad.GetLeftMargin() + 0.135, 1 - (ROOT.gPad.GetTopMargin() - 0.020), "Simulation Preliminary");


    c1.cd(2)
    ROOT.gPad.SetPad(l1, 0.0, l1+l2, 1.0)
    ROOT.gPad.SetLeftMargin(0.0)
    ROOT.gPad.SetRightMargin(0.0)
    ROOT.gPad.SetTopMargin(tm)
    ROOT.gPad.SetBottomMargin(bm)

    h17_1.Draw("HIST")
    f17_1.Draw("2SAME")
    f17_2.Draw("2SAME")
    f17_3.Draw("2SAME")
    f17_4.Draw("2SAME")
    h17_1.Draw("HIST SAME")
    h17_2.Draw("HIST SAME")
    h17_3.Draw("HIST SAME")

    mark.SetTextAlign(31)
    mark.SetTextSize(0.055)
    mark.SetTextFont(42)
    mark.DrawLatex(1.0 - ROOT.gPad.GetLeftMargin() - 0.34, 1 - (ROOT.gPad.GetTopMargin() - 0.020), "arXiv:XXXX.XXXX")

    mark.SetTextAlign(31)
    mark.SetTextSize(0.065)
    mark.SetTextFont(61)
    mark.DrawLatex(1.0 - ROOT.gPad.GetRightMargin() - 0.01, 1 - (ROOT.gPad.GetTopMargin() - 0.020), "2017")

    line1.Draw("SAME")
    line2.Draw("SAME")
    line3.Draw("SAME")
    line4.Draw("SAME")
    line5.Draw("SAME")
    line6.Draw("SAME")

    c1.cd(3)
    ROOT.gPad.SetPad(l1+l2, 0.0, l1+l2+l3, 1.0)
    ROOT.gPad.SetLeftMargin(0.0)
    ROOT.gPad.SetRightMargin(0.0)
    ROOT.gPad.SetTopMargin(tm)
    ROOT.gPad.SetBottomMargin(bm)

    h18_1.Draw("HIST")
    f18_1.Draw("2SAME")
    f18_2.Draw("2SAME")
    f18_3.Draw("2SAME")
    f18_4.Draw("2SAME")
    h18_1.Draw("HIST SAME")
    h18_2.Draw("HIST SAME")
    h18_3.Draw("HIST SAME")

    mark.SetTextAlign(31)
    mark.SetTextSize(0.065)
    mark.SetTextFont(61)

    mark.DrawLatex(1.0 - ROOT.gPad.GetRightMargin() - 0.01, 1 - (ROOT.gPad.GetTopMargin() - 0.020), "2018A")

    line1.Draw("SAME")
    line2.Draw("SAME")
    line3.Draw("SAME")
    line4.Draw("SAME")
    line5.Draw("SAME")
    line6.Draw("SAME")

    c1.cd(4)
    ROOT.gPad.SetPad(1.0-l4, 0.0, 1.0, 1.0)
    ROOT.gPad.SetLeftMargin(0.0)
    ROOT.gPad.SetRightMargin(rm)
    ROOT.gPad.SetTopMargin(tm)
    ROOT.gPad.SetBottomMargin(bm)

    h19_1.Draw("HIST")
    f19_1.Draw("2SAME")
    f19_2.Draw("2SAME")
    f19_3.Draw("2SAME")
    f19_4.Draw("2SAME")
    h19_1.Draw("HIST SAME")
    h19_2.Draw("HIST SAME")
    h19_3.Draw("HIST SAME")

    mark.SetTextAlign(31)
    mark.SetTextSize(0.065)
    mark.SetTextFont(61)

    mark.DrawLatex(1.0 - ROOT.gPad.GetRightMargin() - 0.01, 1 - (ROOT.gPad.GetTopMargin() - 0.020), "2018B")

    line1.Draw("SAME")
    line2.Draw("SAME")
    line3.Draw("SAME")
    line4.Draw("SAME")
    line5.Draw("SAME")

    if args.approved:
        c1.SaveAs("BinEdgesPlot.pdf")
    else:
        c1.SaveAs("BinEdgesPlot_prelim.pdf")
        c1.SaveAs("BinEdgesPlot_prelim.png")

if __name__ == "__main__" :
    main()
