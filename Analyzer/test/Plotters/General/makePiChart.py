#!/bin/python
import os
import ROOT
import argparse

import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams['font.size'] = 18.0

import matplotlib.pyplot as plt

ROOT.gROOT.SetBatch(True)

def main(inpath, outpath, year, channel) :

    if not os.path.exists(outpath):
        os.makedirs(outpath)

    histoName = "h_Njets_%s_ABCD"%(channel)

    ttFile  = ROOT.TFile.Open(inpath + "/" + year + "_TT.root", "READ")
    qcdFile = ROOT.TFile.Open(inpath + "/" + year + "_QCD.root", "READ")
    ttxFile = ROOT.TFile.Open(inpath + "/" + year + "_TTX.root", "READ")
    othFile = ROOT.TFile.Open(inpath + "/" + year + "_BG_OTHER.root", "READ")

    ttHist  = ttFile.Get(histoName)
    qcdHist = qcdFile.Get(histoName) 
    ttxHist = ttxFile.Get(histoName) 
    othHist = othFile.Get(histoName) 

    ttCounts  = ttHist.Integral()
    qcdCounts = qcdHist.Integral()
    ttxCounts = ttxHist.Integral()
    othCounts = othHist.Integral()

    vals   = [ttCounts, ttxCounts, qcdCounts, othCounts]
    colors = ["#9995AB", "#7d99d1", "#85c2a3", "#d4cf87"]

    labels = [r"$t\bar{t}$ + jets", r"$t\bar{t} + X$", "QCD multijet", "Other"]

    fig, ax = plt.subplots()
    fig.subplots_adjust(top=1, bottom=0, right=0.83, left=0.07, hspace=0, wspace=0)
    ax.pie(vals, colors=colors, labels=labels, autopct='%1.1f%%', shadow=False, startangle=40, textprops={"fontsize" : 16})
    ax.axis('equal')    
    ax.margins(0.2, 0)

    fig.savefig("%s/%sMcComposition_%s.pdf"%(outpath, year, channel), dpi=fig.dpi)

if __name__ == "__main__" :

    usage = "usage: %stackPlotter [options]"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("--inpath",  dest="inpath",  help="Path to root files", default="NULL", required=True)
    parser.add_argument("--outpath", dest="outpath", help="Where to put plots", default="NULL", required=True)
    parser.add_argument("--year",    dest="year",    help="which year",                         required=True)
    parser.add_argument("--channel", dest="channel", help="which channel",                      required=True)
    args = parser.parse_args()

    main(args.inpath, args.outpath, args.year, args.channel)
