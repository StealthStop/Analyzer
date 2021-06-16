import ROOT
import sys
import math
import argparse
import numpy as np
import array as arr


def calSignificance(nSigEvents, nBkgEvents, sys=0.2):
    if (nBkgEvents == 0.0):
        return -9999

    significance = nSigEvents / math.sqrt( nBkgEvents + ( sys * nBkgEvents)**2.0 )
    return significance


def calClosure(nEvents_A, nEvents_B, nEvents_C, nEvents_D):
    if (nEvents_A == 0.0 or nEvents_D == 0.0):
        return 9999

    A_pred  = (nEvents_B * nEvents_C) / nEvents_D
    delta_A = (nEvents_A)**0.5
    chi2    = ((A_pred - nEvents_A) / delta_A)**2.0
    return chi2


def getBinEdges(histBkg, histSig):
    lastXBin       = histBkg.GetNbinsX()
    lastYBin       = histBkg.GetNbinsY()
    nXBins         = range(2, lastXBin+1)
    nYBins         = range(2, lastYBin+1)
    significance_A = -9999
    closure_A      = 9999 
    final_xBin     = None    
    final_yBin     = None

    for xBin in nXBins:

        for yBin in nYBins:
            
            nSigEvents_A = ( histSig.Integral(xBin, lastXBin, yBin, lastYBin) )
            nSigEvents_B = ( histSig.Integral(1, xBin-1, yBin, lastYBin) )
            nSigEvents_C = ( histSig.Integral(xBin, lastXBin, 1, yBin-1) )
            nSigEvents_D = ( histSig.Integral(1, xBin-1, 1, yBin-1) )
            
            nBkgEvents_A = ( histBkg.Integral(xBin, lastXBin, yBin, lastYBin) )
            nBkgEvents_B = ( histBkg.Integral(1, xBin-1, yBin, lastYBin) )
            nBkgEvents_C = ( histBkg.Integral(xBin, lastXBin, 1, yBin-1) )
            nBkgEvents_D = ( histBkg.Integral(1, xBin-1, 1, yBin-1) )

            temp_significance_A = calSignificance(nSigEvents_A, nBkgEvents_A)
            temp_closure_A      = calClosure(nBkgEvents_A, nBkgEvents_B, nBkgEvents_C, nBkgEvents_D) 
       
            if ( temp_significance_A > significance_A and temp_closure_A < closure_A ):
                significance_A = temp_significance_A
                closure_A      = temp_closure_A
                final_xBin     = xBin
                final_yBin     = yBin
   
    xBinEdges = histBkg.GetXaxis().GetBinLowEdge(final_xBin)
    yBinEdges = histBkg.GetYaxis().GetBinLowEdge(final_yBin)  
        
    print "x bin edges: ", xBinEdges
    print "y bin edges: ", yBinEdges


def main():    
    usage  = "usage: %prog [options]"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("--year",  dest="year",  help="which year",            required=True)
    parser.add_argument("--path",  dest="path",  help="Input dir with histos", default="./2016_DisCo")
    parser.add_argument("--model", dest="model", help="signal model",          default="RPV")
    parser.add_argument("--mass",  dest="mass",  help="signal mass",           default="550")
    args = parser.parse_args()

    modelDecay = "2t6j"
    if ("SHH" in args.model):
        modelDecay = "2t4b"

    files = {
        "TT"                          : ROOT.TFile.Open(args.path + "/" + args.year + "_TT.root"),
        "QCD"                         : ROOT.TFile.Open(args.path + "/" + args.year + "_QCD.root"),  
        "%s%s"%(args.model,args.mass) : ROOT.TFile.Open(args.path + "/" + args.year + "_%s_%s_mStop-%s.root"%(args.model,modelDecay,args.mass)),
    }
   
    histNames = "h_DoubleDisCo_disc1_disc2_1l_HT300_ge7j_ge1b_Mbl" 

    njets = [
        "_Njets7",
        "_Njets8",
        "_Njets9",
        "_Njets10",
        "_Njets11",
    ]

    for njet in njets: 
        histBkg = files["TT"].Get(histNames + njet)
        histSig = files["%s%s"%(args.model,args.mass)].Get(histNames + njet)
        getBinEdges(histBkg,histSig)

        # put the xbin and ybin edges to txt file
        d = open("BinEdges_%s_%s_%s.txt" %(args.model, args.mass, args.channel), "a")
        d.write("x bin edges: %.3f \n" %(xBinEdges))
        d.write("y bin edges: %.3f \n" %(yBinEdges))

    d.close() 
 
if __name__ == '__main__':
    main()
