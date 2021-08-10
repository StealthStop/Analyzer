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
from DoubleDisCo_BinEdges_Classes import Common_Calculations_Plotters,FinalBinEdges,ValidationRegions 


def main():

    # ---------------------------------------------------------------------------------------------------------
    # command to run this script
    #   -- python DoubleDisCo_BinEdges_Plotters.py --year 2016 --model RPV --mass 550 --channel 0l --metric NN
    #   -- python DoubleDisCo_BinEdges_Plotters.py --year 2016 --model RPV --mass 550 --channel 0l --metric New
    #   -- python DoubleDisCo_BinEdges_Plotters.py --year 2016 --model RPV --mass 550 --channel 1l --metric NN
    #   -- python DoubleDisCo_BinEdges_Plotters.py --year 2016 --model RPV --mass 550 --channel 1l --metric New
    # ---------------------------------------------------------------------------------------------------------  
    usage  = "usage: %prog [options]"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("--year",    dest="year",    help="which year",            required=True)
    parser.add_argument("--path",    dest="path",    help="Input dir with histos", default="/uscms_data/d3/jhiltb/PO_Boxes/Bryan/2016_DisCo_0l_1l_Inputs/")
    parser.add_argument("--model",   dest="model",   help="signal model",          default="RPV")
    parser.add_argument("--mass",    dest="mass",    help="signal mass",           default="550")
    parser.add_argument("--channel", dest="channel", help="0l, 1l",                required=True)
    parser.add_argument("--metric",  dest="metric", help="NN, New",                required=True)
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
        histNames = "h_DoubleDisCo_disc1_disc2_0l"
        njets = [
            "_Njets7",
            "_Njets8",
            "_Njets9",
            "_Njets10",
            "_Njets11",
            "_Njets12",
        ]

    # for 1-lepton
    else:
        histNames = "h_DoubleDisCo_disc1_disc2_1l"
        njets = [
            "_Njets7",
            "_Njets8",
            "_Njets9",
            "_Njets10",
            "_Njets11",
        ]

    # initialize the things
    finalBinEdges           = []; bdEF_FinalEdges         = []; cdiGH_FinalEdges         = []
    final_nTotBkgCount_ABCD = []; final_nTotBkgCount_bdEF = []; final_nTotBkgCount_cdiGH = []

    # ------------------
    # make all tex files
    # -----------------
    # put all latest bin edges to tex file
    h = open("All_LatestBinEdges_%s_%s_%s.tex" %(args.model, args.mass, args.channel), "w")
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
    nEvents_ABCD = open("nEvents_ABCD_%s_%s_%s.tex" %(args.model, args.mass, args.channel), "w")
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
    nEvents_bdEF = open("nEvents_Val_bdEF_%s_%s_%s.tex" %(args.model, args.mass, args.channel), "w")
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
    nEvents_cdiGH = open("nEvents_Val_cdiGH_%s_%s_%s.tex" %(args.model, args.mass, args.channel), "w")
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

        # --------------------------------
        # Common calculations and plotters
        # --------------------------------
        plotter = Common_Calculations_Plotters(args.year, args.model, args.mass, args.channel, args.metric, njet)

        # ---------------
        # Final Bin Edges
        # ---------------
        binEdges = FinalBinEdges(args.year, args.model, args.mass, args.channel, args.metric, njet)
        nTotSigCount_ABCD, nTotBkgCount_ABCD = binEdges.count_Events_inBinEdges(histBkg, histSig)
        finalDisc1Key, finalDisc2Key, significance, closureErr, inverseSignificance, closureErrsList, closureErrUncList, disc1KeyOut, disc2KeyOut, sigFracsA, sigFracsB, sigFracsC, sigFracsD, sigTotFracsA, sigTotFracsB, sigTotFracsC, sigTotFracsD, bkgTotFracsA, bkgTotFracsB, bkgTotFracsC, bkgTotFracsD, finalSigFracA, finalSigFracB, finalSigFracC, finalSigFracD, nEvents_AB, nEvents_AC = binEdges.get_FinalBinEdges(nTotSigCount_ABCD, nTotBkgCount_ABCD, minBkgFrac = 0.01, minSigFrac = 0.1)   
       
        # plot variable vs disc as 1D 
        plotter.plot_VarVsDisc(closureErrsList, closureErrUncList, disc1KeyOut, disc2KeyOut, binWidth/2.0, 1.0, "ABCD Closure", "Closure", 1, njet)
        plotter.plot_VarVsDisc(closureErrsList, closureErrUncList, disc1KeyOut, disc2KeyOut, binWidth/2.0, 1.0, "ABCD Closure", "Closure", 2, njet)


        # ------------------
        # Validation Regions
        # ------------------
        validationRegions = ValidationRegions(args.year, args.model, args.mass, args.channel, njet)
        nTotSigCount_bdEF, nTotBkgCount_bdEF, nTotSigCount_cdiGH, nTotBkgCount_cdiGH = validationRegions.count_Events_inValidationRegions(histBkg, histSig, finalDisc1Key, finalDisc2Key)
        finalDisc1Key_bdEF, finalDisc2Key_bdEF, disc1KeyOut_bdEF, disc2KeyOut_bdEF, sigFracsb, sigFracsd, sigFracsE, sigFracsF, finalSigFracb, finalSigFracd, finalSigFracE, finalSigFracF, nEvents_bd, nEvents_EF = validationRegions.make_ValidationRegionEdges_bdEF(nTotSigCount_bdEF, nTotBkgCount_bdEF, minBkgFrac = 0.01)
        finalDisc1Key_cdiGH, finalDisc2Key_cdiGH, disc1KeyOut_cdiGH, disc2KeyOut_cdiGH, sigFracsc, sigFracsdi, sigFracsG, sigFracsH, finalSigFracc, finalSigFracdi, finalSigFracG, finalSigFracH, nEvents_cdi, nEvents_GH = validationRegions.make_ValidationRegionEdges_cdiGH(nTotSigCount_cdiGH, nTotBkgCount_cdiGH, minBkgFrac = 0.01)

        # plot the signal fractions
        plotter.plot_SigFrac_vsDisc1Disc2(nBins, sigFracsb, sigFracsd, sigFracsE, sigFracsF, disc1KeyOut_bdEF, disc2KeyOut_bdEF, float(finalDisc1Key_bdEF), float(finalDisc2Key_bdEF), minEdge, maxEdge, binWidth, name = "_Val_bdEF")
        plotter.plot_SigFrac_vsDisc1Disc2(nBins, sigFracsc, sigFracsdi, sigFracsG, sigFracsH, disc1KeyOut_cdiGH, disc2KeyOut_cdiGH, float(finalDisc1Key_cdiGH), float(finalDisc2Key_cdiGH), minEdge, maxEdge, binWidth, name = "_Val_cdiGH")    

        # get the lists to make all closure plots 
        finalBinEdges.append((finalDisc1Key, finalDisc2Key))
        bdEF_FinalEdges.append((finalDisc1Key_bdEF, finalDisc2Key_bdEF))
        cdiGH_FinalEdges.append((finalDisc1Key_cdiGH, finalDisc2Key_cdiGH))

        final_nTotBkgCount_ABCD.append(nTotBkgCount_ABCD)
        final_nTotBkgCount_bdEF.append(nTotBkgCount_bdEF)
        final_nTotBkgCount_cdiGH.append(nTotBkgCount_cdiGH)

        # ------------------
        # make all tex files
        # ------------------
        # put all latest bin edges to txt file
        h.write("    \\tiny \\textcolor{ttjetscol}{%s}  & \\tiny \\textcolor{click}{%s} & \\tiny \\textcolor{click}{%s} & \\tiny \\textcolor{rpvcol}{%s} & \\tiny \\textcolor{click}{%s} & \\tiny \\textcolor{click}{%s} & \\tiny \\textcolor{disc}{%s} \\\\" %(njet, finalDisc1Key, finalDisc2Key, finalDisc1Key_bdEF, finalDisc2Key_bdEF, finalDisc1Key_cdiGH, finalDisc2Key_cdiGH))
        h.write("\n")
        h.write("\hline")
        h.write("\n")

        # get the nEvents for each ABCD region
        nEvents_ABCD.write("\\tiny \\textcolor{ttjetscol}{%s} & \\tiny \\textcolor{click}{%.2f} & \scriptsize \\textcolor{click}{%.2f} & \\tiny \\textcolor{click}{%.2f} & \\tiny \\textcolor{click}{%.2f} & \\tiny \\textcolor{click}{%.2f} & \\tiny \\textcolor{click}{%.2f} \\\\" %(njet, finalSigFracA, finalSigFracB, finalSigFracC, finalSigFracD, nEvents_AC, nEvents_AB))
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

    # ----------------------
    # make all closure plots
    # ----------------------
    plotter.make_allClosures(finalBinEdges, final_nTotBkgCount_ABCD, name = "")
    plotter.make_allClosures(bdEF_FinalEdges, final_nTotBkgCount_bdEF, name = "_Val_bdEF")    
    plotter.make_allClosures(cdiGH_FinalEdges, final_nTotBkgCount_cdiGH, name = "_Val_cdiGH")

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

if __name__ == '__main__':
    main()
