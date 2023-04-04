#! /bin/env/python

import os
import ROOT
import argparse
import multiprocessing as mp

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat("")
ROOT.gStyle.SetPaintTextFormat("3.2f")
ROOT.gStyle.SetFrameLineWidth(2)
ROOT.gStyle.SetEndErrorSize(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH2.SetDefaultSumw2()

def makeNDhisto(year, proc, histName, infile, outfile, treeName, histOps):

    selection = histOps["selection"]
    if "Data" in proc:
        selection = selection.replace("${WEIGHT}", "Weight")
    else:
        selection = selection.replace("${WEIGHT}", "TotalWeight_QCDCR")

    infile.cd()
    tree = infile.Get(treeName)

    is2D = False
    if ":" in histOps["variable"]:
        is2D = True

    outfile.cd()

    htemp = None
    if not is2D:
        temph = ROOT.TH1F(histName, "", histOps["xbins"], histOps["xmin"], histOps["xmax"])
    elif h == None:
        temph = ROOT.TH2F(histName, "", histOps["xbins"], histOps["xmin"], histOps["xmax"], histOps["ybins"], histOps["ymin"], histOps["ymax"])

    tree.Draw("%s>>%s"%(histOps["variable"], histName), selection)
       
    temph = ROOT.gDirectory.Get(histName)
    temph.Sumw2()

    temph.Write(histName, ROOT.TObject.kOverwrite)

def processFile(outputDir, inputDir, year, proc, histograms, treeName):
    outfile = ROOT.TFile.Open("%s/%s_%s.root"%(outputDir, year, proc), "UPDATE")
    infile  = ROOT.TFile.Open("%s/%s_%s.root"%(inputDir, year, proc), "READ")
    for histName, histOps in histograms.items():
        if ("SYY" in proc and "SYY" not in histName) or \
           ("RPV" in proc and "RPV" not in histName):
               continue 
        makeNDhisto(year, proc, histName, infile, outfile, treeName, histOps)

    outfile.Close()

usage = "usage: %prog [options]"
parser = argparse.ArgumentParser(usage)
parser.add_argument("--path",     dest="path",     help="Path to ntuples",   default="NULL",        required=True)
parser.add_argument("--tag",      dest="tag",      help="Tag for output",    default="MyTag",                    )
parser.add_argument("--tree",     dest="tree",     help="TTree name to use", default="PreSelection"              )
parser.add_argument("--year",     dest="year",     help="which year",                               required=True)
parser.add_argument("--options",  dest="options",  help="options file",      default="miniTupleDrawer_aux", type=str)
args = parser.parse_args()

# The auxiliary file contains many "hardcoded" items
# describing which histograms to get and how to draw
# them. These things are changed often by the user
# and thus are kept in separate sidecar file.
importedGoods = __import__(args.options)

# Names of histograms, rebinning
histograms = importedGoods.histograms

# Background/signal/data categories to draw
processes  = importedGoods.processes

tag      = args.tag
treeName = args.tree
year     = args.year
inputDir = args.path

base = os.getenv("CMSSW_BASE")

outputDir = "%s/src/Analyzer/Analyzer/test/%s/"%(base,tag)
if not os.path.exists(outputDir):
    os.makedirs(outputDir)

manager = mp.Manager()
pool = mp.Pool(processes=len(processes))

for proc in processes:
    pool.apply_async(processFile, args=(outputDir, inputDir, year, proc, histograms, treeName))

pool.close()
pool.join()
