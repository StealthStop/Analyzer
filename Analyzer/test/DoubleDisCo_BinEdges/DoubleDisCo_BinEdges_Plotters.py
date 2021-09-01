import ROOT
import os
import sys
import math
import ctypes
import argparse
import numpy as np
import array as arr

import matplotlib as mpl
mpl.use('Agg')

import matplotlib.lines as ml
import matplotlib.pyplot as plt
#mpl.__version__
#print mpl.__version__
# mpl version is 2.2.5

from matplotlib.colors import LogNorm
from DoubleDisCo_BinEdges_Classes import Common_Calculations_Plotters, FinalBinEdges, ValidationRegions, ValidationRegions_SubDivisionsOfABCD 


def main():

    # ----------------------------------------------------------------------------------------------------------------------
    # command to run this script
    #   -- python DoubleDisCo_BinEdges_Plotters.py --year 2016 --model RPV --mass 550 --channel 0l --metric NN --edges final
    #   -- python DoubleDisCo_BinEdges_Plotters.py --year 2016 --model RPV --mass 550 --channel 0l --metric NN --edges fixed
    #   -- python DoubleDisCo_BinEdges_Plotters.py --year 2016 --model RPV --mass 550 --channel 0l --metric New --edges final
    #   -- python DoubleDisCo_BinEdges_Plotters.py --year 2016 --model RPV --mass 550 --channel 0l --metric New --edges fixed
    #   -- python DoubleDisCo_BinEdges_Plotters.py --year 2016 --model RPV --mass 550 --channel 1l --metric NN --edges final
    #   -- python DoubleDisCo_BinEdges_Plotters.py --year 2016 --model RPV --mass 550 --channel 1l --metric NN --edges fixed
    #   -- python DoubleDisCo_BinEdges_Plotters.py --year 2016 --model RPV --mass 550 --channel 1l --metric New --edges final
    #   -- python DoubleDisCo_BinEdges_Plotters.py --year 2016 --model RPV --mass 550 --channel 1l --metric New --edges fixed
    # -----------------------------------------------------------------------------------------------------------------------  
    usage  = "usage: %prog [options]"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("--year",    dest="year",    help="which year",            required=True)
    parser.add_argument("--path",    dest="path",    help="Input dir with histos", default="/uscms/home/jhiltb/nobackup/PO_Boxes/DoubleDisCo_Reg_0L_loose_1L_RPV_2016_20210818_Output/")
    parser.add_argument("--model",   dest="model",   help="signal model",          default="RPV")
    parser.add_argument("--mass",    dest="mass",    help="signal mass",           default="550")
    parser.add_argument("--channel", dest="channel", help="0l, 1l",                required=True)
    parser.add_argument("--metric",  dest="metric",  help="NN, New",               required=True)
    parser.add_argument("--edges",   dest="edges",   help="final, fixed",          required=True, type=str)
    args = parser.parse_args()

    modelDecay = "2t6j"
    if ("SHH" in args.model):
        modelDecay = "2t4b"

    files = {
        "TT"                          : ROOT.TFile.Open(args.path + "/" + args.year + "_TT.root"),
        "QCD"                         : ROOT.TFile.Open(args.path + "/" + args.year + "_QCD.root"),  
        "%s%s"%(args.model,args.mass) : ROOT.TFile.Open(args.path + "/" + args.year + "_%s_%s_mStop-%s.root"%(args.model,modelDecay,args.mass)),
    }

    # for 0-lepton 
    if args.channel == "0l":
        histNames = "h_DoubleDisCo_disc1_disc2_0l_Njets"
        njets = ["6", "7", "8", "9", "10", "11", "12"]

    # for 1-lepton
    else:
        histNames = "h_DoubleDisCo_disc1_disc2_1l_Njets"
        njets = ["7", "8", "9", "10", "11"]

    # initialize the things
    finalBinEdges       = []; final_nTotBkgCount_ABCD  = [] 
    bdEF_FinalEdges     = []; final_nTotBkgCount_bdEF  = []
    cdiGH_FinalEdges    = []; final_nTotBkgCount_cdiGH = []
    final_subDivD_Edges = []; final_nTotBkgCount       = []


    # ------------------
    # make all tex files
    # -----------------
    # put all latest bin edges to tex file
    h = open("tables/%s/%s_All_%sBinEdges_%s_%s_%s.tex" %(args.channel, args.year, args.edges, args.model, args.mass, args.channel), "w")
    h.write("\\resizebox{\linewidth}{!}{%")
    h.write("\n")
    h.write("\\begin{tabular}{| c | c | c | c | c | c | c |}")
    h.write("\n")
    h.write("\hline")
    h.write("\n")
    h.write("\\textcolor{ttjetscol}{NJets} &  \multicolumn{2}{c|} {\\textcolor{click}{ABCD}} & \multicolumn{2}{c|} {\\textcolor{rpvcol}{B'D'EF}} & \multicolumn{2}{c|} {\\textcolor{disc}{C'D'GH}} \\\\")
    h.write("\n")
    h.write("\hline")
    h.write("\n")
    h.write("& \scriptsize \\textcolor{click}{disc1 edges} & \scriptsize \\textcolor{click}{disc2 edges} & \scriptsize \\textcolor{rpvcol}{disc1 edges} & \scriptsize \\textcolor{click}{disc2 edges} & \scriptsize \\textcolor{click}{disc1 edges} & \scriptsize \\textcolor{disc}{disc2 edges} \\\\")
    h.write("\n")
    h.write("\hline")
    h.write("\n")

    # get the nEvents for each ABCD region
    nEvents_ABCD = open("tables/%s/%s_nEvents_%s_ABCD_%s_%s_%s.tex" %(args.channel, args.year, args.edges, args.model, args.mass, args.channel), "w")
    nEvents_ABCD.write("\\resizebox{\linewidth}{!}{%")
    nEvents_ABCD.write("\n")
    nEvents_ABCD.write("\\begin{tabular}{| c | c | c | c | c | c | c |}")
    nEvents_ABCD.write("\n")
    nEvents_ABCD.write("\hline")
    nEvents_ABCD.write("\n")
    nEvents_ABCD.write("\scriptsize \\textcolor{ttjetscol}{NJets} & \scriptsize \\textcolor{click}{sigFracA} & \scriptsize \\textcolor{click}{sigFracB} & \scriptsize \\textcolor{click}{sigFracC} & \scriptsize \\textcolor{click}{sigFracD} & \scriptsize \\textcolor{click}{nBkgEvents(A+C)} & \scriptsize \\textcolor{click}{nBkgEvents(A+B)} \\\\")
    nEvents_ABCD.write("\n")
    nEvents_ABCD.write("\hline")
    nEvents_ABCD.write("\n")

    # get the nEvents for each B'D'EF region 
    nEvents_bdEF = open("tables/%s/%s_nEvents_%s_Val_bdEF_%s_%s_%s.tex" %(args.channel, args.year, args.edges, args.model, args.mass, args.channel), "w")
    nEvents_bdEF.write("\\resizebox{\linewidth}{!}{%")
    nEvents_bdEF.write("\n")
    nEvents_bdEF.write("\\begin{tabular}{| c | c | c | c | c | c | c |}")
    nEvents_bdEF.write("\n")
    nEvents_bdEF.write("\hline")
    nEvents_bdEF.write("\n")
    nEvents_bdEF.write("\scriptsize \\textcolor{ttjetscol}{NJets} & \scriptsize \\textcolor{rpvcol}{sigFracB'} & \scriptsize \\textcolor{rpvcol}{sigFracD'} & \scriptsize \\textcolor{rpvcol}{sigFracE} & \scriptsize \\textcolor{rpvcol}{sigFracF} & \scriptsize \\textcolor{rpvcol}{nBkgEvents(B'+D')} & \scriptsize \\textcolor{rpvcol}{nBkgEvents(E+F)} \\\\")
    nEvents_bdEF.write("\n")
    nEvents_bdEF.write("\hline")
    nEvents_bdEF.write("\n")

    # get the nEvents for each C'D'GH region
    nEvents_cdiGH = open("tables/%s/%s_nEvents_%s_Val_cdiGH_%s_%s_%s.tex" %(args.channel, args.year, args.edges, args.model, args.mass, args.channel), "w")
    nEvents_cdiGH.write("\\resizebox{\linewidth}{!}{%")
    nEvents_cdiGH.write("\n")
    nEvents_cdiGH.write("\\begin{tabular}{| c | c | c | c | c | c | c |}")
    nEvents_cdiGH.write("\n")
    nEvents_cdiGH.write("\hline")
    nEvents_cdiGH.write("\n")
    nEvents_cdiGH.write("\scriptsize \\textcolor{ttjetscol}{NJets} & \scriptsize \\textcolor{disc}{sigFracC'} & \scriptsize \\textcolor{disc}{sigFracD'} & \scriptsize \\textcolor{disc}{sigFracG} & \scriptsize \\textcolor{disc}{sigFracH} & \scriptsize \\textcolor{disc}{nBkgEvents(C'+D')} & \scriptsize \\textcolor{disc}{nBkgEvents(G+H)} \\\\")
    nEvents_cdiGH.write("\n")
    nEvents_cdiGH.write("\hline")
    nEvents_cdiGH.write("\n") 

    # get the table for Validation region sub-division D
    d = open("tables/%s/%s_%s_FinalBinEdges_%s_%s_Val_subDivD_%s_%s.tex" %(args.channel, args.year, args.edges, args.model, args.mass, args.channel, args.metric), "w")
    d.write("\\resizebox{\linewidth}{!}{%")
    d.write("\n")
    d.write("\def\\arraystretch{0.6}")
    d.write("\n")
    d.write("\\begin{tabular}{| c | c | c | c | c | c | c | c | c | c | c |}")
    d.write("\n")
    d.write("\hline")
    d.write("\n")
    d.write("\\textcolor{ttjetscol}{NJets} &  \multicolumn{2}{c|} {\\textcolor{click}{Edges}} & \multicolumn{4}{c|} {\\textcolor{massReg}{SigFracs}} & \multicolumn{4}{c|} {\\textcolor{disc}{nEvents(Sig+Bkg)}} \\\\")
    d.write("\n")
    d.write("\hline")
    d.write("\n")
    d.write("& \scriptsize \\textcolor{click}{disc1} & \scriptsize \\textcolor{click}{disc2} & \scriptsize \\textcolor{massReg}{sigFrac\_dA} & \scriptsize \\textcolor{massReg}{sigFrac\_dB} & \scriptsize \\textcolor{massReg}{SigFrac\_dC} & \scriptsize \\textcolor{massReg}{SigFrac\_dD} & \scriptsize \\textcolor{disc}{in region dA} & \scriptsize \\textcolor{disc}{in region dB} & \scriptsize \\textcolor{disc}{in region dC} & \scriptsize \\textcolor{disc}{in region dD} \\\\")
    d.write("\n")
    d.write("\hline")
    d.write("\n")

    # ---------------
    # loop over njets
    # --------------- 
    for njet in njets: 

        histBkg = files["TT"].Get(histNames + njet)
        histSig = files["%s%s"%(args.model,args.mass)].Get(histNames + njet)

        minEdge  = histSig.GetXaxis().GetBinLowEdge(1) 
        maxEdge  = histSig.GetXaxis().GetBinUpEdge(histBkg.GetNbinsX())
        binWidth = histSig.GetXaxis().GetBinWidth(1)
        nBins    = histSig.GetNbinsX()
        edges    = np.arange(minEdge, maxEdge, binWidth)

        # --------------------------------
        # Common calculations and plotters
        # --------------------------------
        plotter = Common_Calculations_Plotters(args.year, args.model, args.mass, args.channel, args.metric, args.edges, njet)

        # ---------------
        # Final Bin Edges
        # ---------------
        binEdges = FinalBinEdges(args.year, args.model, args.mass, args.channel, args.metric, args.edges, njet)
        nTotSigCount_ABCD, nTotBkgCount_ABCD = binEdges.count_Events_inBinEdges(histBkg, histSig)
        finalDisc1Key, finalDisc2Key, significance, closureErr, inverseSignificance, closureErrsList, closureErrUncList, sigUncs, disc1KeyOut, disc2KeyOut, weighted_Sig_A, weighted_Bkg_A, weighted_SigUnc_A, weighted_BkgUnc_A, sigFracsA, sigFracsB, sigFracsC, sigFracsD, sigFracsErrA, sigFracsErrB, sigFracsErrC, sigFracsErrD, sigTotFracsA, sigTotFracsB, sigTotFracsC, sigTotFracsD, bkgTotFracsA, bkgTotFracsB, bkgTotFracsC, bkgTotFracsD, finalSigFracA, finalSigFracB, finalSigFracC, finalSigFracD, nEvents_AB, nEvents_AC, final_nBkgEvents_A, final_nBkgEvents_B, final_nBkgEvents_C, final_nBkgEvents_D = binEdges.get_FinalBinEdges(nTotSigCount_ABCD, nTotBkgCount_ABCD, minBkgFrac = 0.01, minSigFrac = 0.1)   
        
        # plot variable vs disc as 1D 
        plotter.plot_VarVsDisc(closureErrsList, closureErrUncList, disc1KeyOut, disc2KeyOut, binWidth/2.0, 1.0, "ABCD Closure", "Closure", 1, njet, name = "ABCD")
        plotter.plot_VarVsDisc(closureErrsList, closureErrUncList, disc1KeyOut, disc2KeyOut, binWidth/2.0, 1.0, "ABCD Closure", "Closure", 2, njet, name = "ABCD")
        plotter.plot_VarVsDisc(np.reciprocal(inverseSignificance), sigUncs, disc1KeyOut, disc2KeyOut, binWidth/2.0, 5.0, "Significance", "Significance", 1, njet, name = "ABCD")
        plotter.plot_VarVsDisc(np.reciprocal(inverseSignificance), sigUncs, disc1KeyOut, disc2KeyOut, binWidth/2.0, 5.0, "Significance", "Significance", 2, njet, name = "ABCD")
        plotter.plot_VarVsDisc(weighted_Sig_A, weighted_SigUnc_A, disc1KeyOut, disc2KeyOut, binWidth/2.0, -1.0, "Weighted Signal Events", "wSigEvts", 1, njet, name = "ABCD")
        plotter.plot_VarVsDisc(weighted_Sig_A, weighted_SigUnc_A, disc1KeyOut, disc2KeyOut, binWidth/2.0, -1.0, "Weighted Signal Events", "wSigEvts", 2, njet, name = "ABCD")
        plotter.plot_VarVsDisc(weighted_Bkg_A, weighted_BkgUnc_A, disc1KeyOut, disc2KeyOut, binWidth/2.0, -1.0, "Weighted Background Events", "wBkgEvts", 1, njet, name = "ABCD")
        plotter.plot_VarVsDisc(weighted_Bkg_A, weighted_BkgUnc_A, disc1KeyOut, disc2KeyOut, binWidth/2.0, -1.0, "Weighted Background Events", "wBkgEvts", 2, njet, name = "ABCD")
        plotter.plot_VarVsDisc(sigFracsA, sigFracsErrA, disc1KeyOut, disc2KeyOut, binWidth/2.0, 0.8, "SigFracsA", "SigFracsA", 1, njet, name = "ABCD")
        plotter.plot_VarVsDisc(sigFracsB, sigFracsErrB, disc1KeyOut, disc2KeyOut, binWidth/2.0, 0.8, "SigFracsB", "SigFracsB", 1, njet, name = "ABCD")
        plotter.plot_VarVsDisc(sigFracsC, sigFracsErrC, disc1KeyOut, disc2KeyOut, binWidth/2.0, 0.8, "SigFracsC", "SigFracsC", 1, njet, name = "ABCD")
        plotter.plot_VarVsDisc(sigFracsD, sigFracsErrD, disc1KeyOut, disc2KeyOut, binWidth/2.0, 0.8, "SigFracsD", "SigFracsD", 1, njet, name = "ABCD") 
        plotter.plot_VarVsDisc(sigFracsA, sigFracsErrA, disc1KeyOut, disc2KeyOut, binWidth/2.0, 0.8, "SigFracsA", "SigFracsA", 2, njet, name = "ABCD")
        plotter.plot_VarVsDisc(sigFracsB, sigFracsErrB, disc1KeyOut, disc2KeyOut, binWidth/2.0, 0.8, "SigFracsB", "SigFracsB", 2, njet, name = "ABCD")
        plotter.plot_VarVsDisc(sigFracsC, sigFracsErrC, disc1KeyOut, disc2KeyOut, binWidth/2.0, 0.8, "SigFracsC", "SigFracsC", 2, njet, name = "ABCD")
        plotter.plot_VarVsDisc(sigFracsD, sigFracsErrD, disc1KeyOut, disc2KeyOut, binWidth/2.0, 0.8, "SigFracsD", "SigFracsD", 2, njet, name = "ABCD")

        # plot 2Ds
        plotter.plot_Significance_vsDisc1Disc2(nBins, inverseSignificance, disc1KeyOut, disc2KeyOut, float(finalDisc1Key), float(finalDisc2Key), minEdge, maxEdge, binWidth, njet, name="ABCD")
        plotter.plot_ClosureError_vsDisc1Disc2(nBins, closureErrsList, closureErrUncList, disc1KeyOut, disc2KeyOut, float(finalDisc1Key), float(finalDisc2Key), minEdge, maxEdge, binWidth, njet, name="ABCD")
        plotter.plot_inverseSignificance_vsClosureErr(significance, closureErr, inverseSignificance, closureErrsList, edges, float(finalDisc1Key), float(finalDisc2Key), njet, name="ABCD")
        plotter.plot_SigFrac_vsDisc1Disc2(nBins, sigFracsA, sigFracsB, sigFracsC, sigFracsD, disc1KeyOut, disc2KeyOut, float(finalDisc1Key), float(finalDisc2Key), minEdge, maxEdge, binWidth, name = "ABCD")
        plotter.plot_SigTotFrac_vsDisc1Disc2(nBins, sigTotFracsA, sigTotFracsB, sigTotFracsC, sigTotFracsD, disc1KeyOut, disc2KeyOut, float(finalDisc1Key), float(finalDisc2Key), minEdge, maxEdge, binWidth, name="ABCD")
        plotter.plot_BkgTotFrac_vsDisc1Disc2(nBins, bkgTotFracsA, bkgTotFracsB, bkgTotFracsC, bkgTotFracsD, disc1KeyOut, disc2KeyOut, float(finalDisc1Key), float(finalDisc2Key), minEdge, maxEdge, binWidth, name="ABCD")

        # ------------------
        # Validation Regions
        # ------------------
        validationRegions = ValidationRegions(args.year, args.model, args.mass, args.channel, args.metric, args.edges, njet)
        nTotSigCount_bdEF, nTotBkgCount_bdEF, nTotSigCount_cdiGH, nTotBkgCount_cdiGH = validationRegions.count_Events_inValidationRegions(histBkg, histSig, finalDisc1Key, finalDisc2Key)
        finalDisc1Key_bdEF, finalDisc2Key_bdEF, disc1KeyOut_bdEF, disc2KeyOut_bdEF, closureErrsList_bdEF, closureErrUncList_bdEF, weighted_Sig_b, weighted_Bkg_b, weighted_SigUnc_b, weighted_BkgUnc_b, sigFracsb, sigFracsd, sigFracsE, sigFracsF, sigFracsErrb, sigFracsErrd, sigFracsErrE, sigFracsErrF, finalSigFracb, finalSigFracd, finalSigFracE, finalSigFracF, nEvents_bd, nEvents_EF = validationRegions.make_ValidationRegionEdges_bdEF(nTotSigCount_bdEF, nTotBkgCount_bdEF, final_nBkgEvents_A, final_nBkgEvents_C)
        finalDisc1Key_cdiGH, finalDisc2Key_cdiGH, disc1KeyOut_cdiGH, disc2KeyOut_cdiGH, closureErrsList_cdiGH, closureErrUncList_cdiGH, weighted_Sig_c, weighted_Bkg_c, weighted_SigUnc_c, weighted_BkgUnc_c, sigFracsc, sigFracsdi, sigFracsG, sigFracsH, sigFracsErrc, sigFracsErrdi, sigFracsErrG, sigFracsErrH, finalSigFracc, finalSigFracdi, finalSigFracG, finalSigFracH, nEvents_cdi, nEvents_GH = validationRegions.make_ValidationRegionEdges_cdiGH(nTotSigCount_cdiGH, nTotBkgCount_cdiGH, final_nBkgEvents_A, final_nBkgEvents_B)

        # get the lists to make all closure plots 
        finalBinEdges.append((finalDisc1Key, finalDisc2Key))
        bdEF_FinalEdges.append((finalDisc1Key_bdEF, finalDisc2Key_bdEF))
        cdiGH_FinalEdges.append((finalDisc1Key_cdiGH, finalDisc2Key_cdiGH))
        final_nTotBkgCount_ABCD.append(nTotBkgCount_ABCD)
        final_nTotBkgCount_bdEF.append(nTotBkgCount_bdEF)
        final_nTotBkgCount_cdiGH.append(nTotBkgCount_cdiGH)

        # plot variable vs disc as 1D
        plotter.plot_VarVsDisc_Val(closureErrsList_bdEF, closureErrUncList_bdEF, disc1KeyOut_bdEF, float(finalDisc2Key), binWidth/2.0, 1.0, "ABCD Closure",     "Closure",  2, njet, name = "Val_bdEF", col = "#DDBB87")
        plotter.plot_VarVsDisc_Val(weighted_Sig_b, weighted_SigUnc_b, disc1KeyOut_bdEF, float(finalDisc2Key), binWidth/2.0, -1.0, "Weighted Signal Events",     "wSigEvts", 2, njet, name = "Val_bdEF", col = "#DDBB87")
        plotter.plot_VarVsDisc_Val(weighted_Bkg_b, weighted_BkgUnc_b, disc1KeyOut_bdEF, float(finalDisc2Key), binWidth/2.0, -1.0, "Weighted Background Events", "wBkgEvts", 2, njet, name = "Val_bdEF", col = "#DDBB87")
        plotter.plot_VarVsDisc_Val(sigFracsb, sigFracsErrb, disc1KeyOut_bdEF, float(finalDisc2Key), binWidth/2.0, 0.8, "SigFracsb", "SigFracsb", 2, njet, name = "Val_bdEF", col = "#DDBB87")
        plotter.plot_VarVsDisc_Val(sigFracsd, sigFracsErrd, disc1KeyOut_bdEF, float(finalDisc2Key), binWidth/2.0, 0.8, "SigFracsd", "SigFracsd", 2, njet, name = "Val_bdEF", col = "#DDBB87")
        plotter.plot_VarVsDisc_Val(sigFracsE, sigFracsErrE, disc1KeyOut_bdEF, float(finalDisc2Key), binWidth/2.0, 0.8, "SigFracsE", "SigFracsE", 2, njet, name = "Val_bdEF", col = "#DDBB87")
        plotter.plot_VarVsDisc_Val(sigFracsF, sigFracsErrF, disc1KeyOut_bdEF, float(finalDisc2Key), binWidth/2.0, 0.8, "SigFracsF", "SigFracsF", 2, njet, name = "Val_bdEF", col = "#DDBB87")
        plotter.plot_VarVsDisc_Val(closureErrsList_cdiGH, closureErrUncList_cdiGH, disc2KeyOut_cdiGH, float(finalDisc1Key), binWidth/2.0, 1.0, "ABCD Closure",   "Closure",  1, njet, name = "Val_cdiGH", col = "#429c93")
        plotter.plot_VarVsDisc_Val(weighted_Sig_c, weighted_SigUnc_c, disc2KeyOut_cdiGH, float(finalDisc1Key), binWidth/2.0, -1.0, "Weighted Signal Events",     "wSigEvts", 1, njet, name = "Val_cdiGH", col = "#429c93")
        plotter.plot_VarVsDisc_Val(weighted_Bkg_c, weighted_BkgUnc_c, disc2KeyOut_cdiGH, float(finalDisc1Key), binWidth/2.0, -1.0, "Weighted Background Events", "wBkgEvts", 1, njet, name = "Val_cdiGH", col = "#429c93")
        plotter.plot_VarVsDisc_Val(sigFracsc,  sigFracsErrc,  disc2KeyOut_cdiGH, float(finalDisc1Key), binWidth/2.0, 0.8, "SigFracsc",  "SigFracsc",  1, njet, name = "Val_cdiGH", col = "#429c93")
        plotter.plot_VarVsDisc_Val(sigFracsdi, sigFracsErrdi, disc2KeyOut_cdiGH, float(finalDisc1Key), binWidth/2.0, 0.8, "SigFracsdi", "SigFracsdi", 1, njet, name = "Val_cdiGH", col = "#429c93")
        plotter.plot_VarVsDisc_Val(sigFracsG,  sigFracsErrG,  disc2KeyOut_cdiGH, float(finalDisc1Key), binWidth/2.0, 0.8, "SigFracsG",  "SigFracsG",  1, njet, name = "Val_cdiGH", col = "#429c93")
        plotter.plot_VarVsDisc_Val(sigFracsH,  sigFracsErrH,  disc2KeyOut_cdiGH, float(finalDisc1Key), binWidth/2.0, 0.8, "SigFracsH",  "SigFracsH",  1, njet, name = "Val_cdiGH", col = "#429c93")

        # plot 2Ds
        plotter.plot_ClosureError_vsDisc1Disc2(nBins, closureErrsList_bdEF, closureErrUncList_bdEF, disc1KeyOut_bdEF, disc2KeyOut_bdEF, float(finalDisc1Key_bdEF), float(finalDisc2Key_bdEF), minEdge, maxEdge, binWidth, njet, name="Val_bdEF")
        plotter.plot_ClosureError_vsDisc1Disc2(nBins, closureErrsList_cdiGH, closureErrUncList_cdiGH, disc1KeyOut_cdiGH, disc2KeyOut_cdiGH, float(finalDisc1Key_cdiGH), float(finalDisc2Key_cdiGH), minEdge, maxEdge, binWidth, njet, name="Val_cdiGH")
        #plotter.plot_SigFrac_vsDisc1Disc2(nBins, sigFracsb, sigFracsd, sigFracsE, sigFracsF, disc1KeyOut_bdEF, disc2KeyOut_bdEF, float(finalDisc1Key_bdEF), float(finalDisc2Key_bdEF), minEdge, maxEdge, binWidth, name = "_Val_bdEF")
        #plotter.plot_SigFrac_vsDisc1Disc2(nBins, sigFracsc, sigFracsdi, sigFracsG, sigFracsH, disc1KeyOut_cdiGH, disc2KeyOut_cdiGH, float(finalDisc1Key_cdiGH), float(finalDisc2Key_cdiGH), minEdge, maxEdge, binWidth, name = "_Val_cdiGH")

        # put all latest bin edges to txt file
        h.write("    \\tiny \\textcolor{ttjetscol}{%s}  & \\tiny \\textcolor{click}{%s} & \\tiny \\textcolor{click}{%s} & \\tiny \\textcolor{rpvcol}{%s} & \\tiny \\textcolor{click}{%s} & \\tiny \\textcolor{click}{%s} & \\tiny \\textcolor{disc}{%s} \\\\" %(njet, finalDisc1Key, finalDisc2Key, finalDisc1Key_bdEF, finalDisc2Key_bdEF, finalDisc1Key_cdiGH, finalDisc2Key_cdiGH))
        h.write("\n")
        h.write("\hline")
        h.write("\n")

        # get the nEvents for each ABCD region
        nEvents_ABCD.write("\\tiny \\textcolor{ttjetscol}{%s} & \\tiny \\textcolor{click}{%.2f} & \\tiny \\textcolor{click}{%.2f} & \\tiny \\textcolor{click}{%.2f} & \\tiny \\textcolor{click}{%.2f} & \\tiny \\textcolor{click}{%.2f} & \\tiny \\textcolor{click}{%.2f} \\\\" %(njet, finalSigFracA, finalSigFracB, finalSigFracC, finalSigFracD, nEvents_AC, nEvents_AB))
        nEvents_ABCD.write("\n")
        nEvents_ABCD.write("\hline")
        nEvents_ABCD.write("\n")
     
        # get the nEvents for each B'D'EF region 
        nEvents_bdEF.write("\\tiny \\textcolor{ttjetscol}{%s} & \\tiny \\textcolor{rpvcol}{%.2f} & \\tiny \\textcolor{rpvcol}{%.2f} & \\tiny \\textcolor{rpvcol}{%.2f} & \\tiny \\textcolor{rpvcol}{%.2f} & \\tiny \\textcolor{rpvcol}{%.2f} & \\tiny \\textcolor{rpvcol}{%.2f} \\\\" %(njet, finalSigFracb, finalSigFracd, finalSigFracE, finalSigFracF, nEvents_bd, nEvents_EF))
        nEvents_bdEF.write("\n")
        nEvents_bdEF.write("\hline")
        nEvents_bdEF.write("\n")

        # get the nEvents for each C'D'GH region 
        nEvents_cdiGH.write("\\tiny \\textcolor{ttjetscol}{%s} & \\tiny \\textcolor{disc}{%.2f} & \\tiny \\textcolor{disc}{%.2f} & \\tiny \\textcolor{disc}{%.2f} & \\tiny \\textcolor{disc}{%.2f} & \\tiny \\textcolor{disc}{%.2f} & \\tiny \\textcolor{disc}{%.2f} \\\\" %(njet, finalSigFracc, finalSigFracdi, finalSigFracG, finalSigFracH, nEvents_cdi, nEvents_GH))
        nEvents_cdiGH.write("\n")
        nEvents_cdiGH.write("\hline")
        nEvents_cdiGH.write("\n")

        # -----------------------------------------------------
        # Validation Regions as SubDivisions of each A, B, C, D
        # -----------------------------------------------------
        val_SubDivisionsOfABCD = ValidationRegions_SubDivisionsOfABCD(args.year, args.model, args.mass, args.channel, args.metric, args.edges, njet)
        nTotSigCount, nTotBkgCount = val_SubDivisionsOfABCD.count_Events_inSubDivisionsOfABCD(histBkg, histSig, float(finalDisc1Key), float(finalDisc2Key))
        finalDisc1Key_subDivD, finalDisc2Key_subDivD, finalSigFrac_dA, finalSigFrac_dB, finalSigFrac_dC, finalSigFrac_dD, final_nTot_SigBkg_dA, final_nTot_SigBkg_dB, final_nTot_SigBkg_dC, final_nTot_SigBkg_dD, closureErrsList_subDivD, closureErrUncList_subDivD, disc1KeyOut_subDivD, disc2KeyOut_subDivD = val_SubDivisionsOfABCD.get_ValidationEdges_SubDivisionsOfABCD(nTotSigCount, nTotBkgCount, float(finalDisc1Key), float(finalDisc2Key)) 

        # get the lists to make all closure plots 
        final_subDivD_Edges.append((finalDisc1Key_subDivD, finalDisc2Key_subDivD))
        final_nTotBkgCount.append(nTotBkgCount)

        # get the table for Validation region sub-division D
        d.write("\\tiny \\textcolor{ttjetscol}{%s} &  \\tiny \\textcolor{click}{%s} & \\tiny \\textcolor{click}{%s} & \\tiny \\textcolor{massReg}{%.3f} & \\tiny \\textcolor{massReg}{%.3f} & \\tiny \\textcolor{massReg}{%.3f} & \\tiny \\textcolor{massReg}{%.3f} & \\tiny \\textcolor{disc}{%.2f} & \\tiny \\textcolor{disc}{%.2f} & \\tiny \\textcolor{disc}{%.2f} & \\tiny \\textcolor{disc}{%.2f} \\\\" %(njet, finalDisc1Key_subDivD, finalDisc2Key_subDivD, finalSigFrac_dA, finalSigFrac_dB, finalSigFrac_dC, finalSigFrac_dD, final_nTot_SigBkg_dA, final_nTot_SigBkg_dB, final_nTot_SigBkg_dC, final_nTot_SigBkg_dD) )
        d.write("\n")
        d.write("\hline")
        d.write("\n")

        # plot variable vs disc as 1D 
        plotter.plot_VarVsDisc(closureErrsList_subDivD, closureErrUncList_subDivD, disc1KeyOut_subDivD, disc2KeyOut_subDivD, binWidth/2.0, 1.0, "Sub-Division D Closure", "Closure", 1, njet, name = "Val_SubDivD")
        plotter.plot_VarVsDisc(closureErrsList_subDivD, closureErrUncList_subDivD, disc1KeyOut_subDivD, disc2KeyOut_subDivD, binWidth/2.0, 1.0, "Sub-Division D Closure", "Closure", 2, njet, name = "Val_SubDivD")

        # ------------------------------------
        # plot Disc1s vs Disc2s with all edges
        # ------------------------------------
        if args.metric == "New":
            plotter.plot_Disc1VsDisc2(histSig, float(finalDisc1Key), float(finalDisc2Key), float(finalDisc1Key_bdEF), float(finalDisc2Key_cdiGH), tag = "RPV550", name = "Val_BD_CD", col1="yellow", col2="lime")
            plotter.plot_Disc1VsDisc2(histBkg, float(finalDisc1Key), float(finalDisc2Key), float(finalDisc1Key_bdEF), float(finalDisc2Key_cdiGH), tag = "TT",     name = "Val_BD_CD", col1="yellow", col2="lime")
            plotter.plot_Disc1VsDisc2(histSig, float(finalDisc1Key), float(finalDisc2Key), float(finalDisc1Key_subDivD), float(finalDisc2Key_subDivD), tag = "RPV550", name = "Val_SubDivD", col1="white", col2="white")
            plotter.plot_Disc1VsDisc2(histBkg, float(finalDisc1Key), float(finalDisc2Key), float(finalDisc1Key_subDivD), float(finalDisc2Key_subDivD), tag = "TT",     name = "Val_SubDivD", col1="white", col2="white")

    # ----------------------
    # make all closure plots
    # ----------------------
    plotter.make_allClosures(finalBinEdges,       final_nTotBkgCount_ABCD,  name = "ABCD"       )
    plotter.make_allClosures(bdEF_FinalEdges,     final_nTotBkgCount_bdEF,  name = "Val_bdEF"   )    
    plotter.make_allClosures(cdiGH_FinalEdges,    final_nTotBkgCount_cdiGH, name = "Val_cdiGH"  )
    plotter.make_allClosures(final_subDivD_Edges, final_nTotBkgCount,       name = "Val_subDivD")

    # ------------------
    # make all tex files
    # ------------------
    # put all latest bin edges to txt file
    h.write("\end{tabular}")
    h.write("\n")
    h.write("}")
    h.write("\n")
    h.write("\n")
    h.close()
    
    # get the nEvents for each ABCD region
    nEvents_ABCD.write("\end{tabular}")
    nEvents_ABCD.write("\n")
    nEvents_ABCD.write("}")
    nEvents_ABCD.write("\n")
    nEvents_ABCD.write("\n")
    nEvents_ABCD.close()

    # get the nEvents for each B'D'EF region 
    nEvents_bdEF.write("\end{tabular}")
    nEvents_bdEF.write("\n")
    nEvents_bdEF.write("}")
    nEvents_bdEF.write("\n")
    nEvents_bdEF.write("\n")
    nEvents_bdEF.close()

    # get the nEvents for each C'D'GH region 
    nEvents_cdiGH.write("\end{tabular}")
    nEvents_cdiGH.write("\n")
    nEvents_cdiGH.write("}")
    nEvents_cdiGH.write("\n")
    nEvents_cdiGH.write("\n")
    nEvents_cdiGH.close()

    # get the table for Validation region sub-division D
    d.write("\end{tabular}")
    d.write("\n")
    d.write("}")
    d.write("\n")
    d.write("\n")
    d.close()


if __name__ == '__main__':
    main()
