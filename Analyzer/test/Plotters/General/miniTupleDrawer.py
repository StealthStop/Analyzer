#! /bin/env/python

import os
import re
import argparse
import multiprocessing as mp

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH2.SetDefaultSumw2()

# Routine that is called for each individual histogram that is to be 
# drawn from the input tree. All information about what to draw, selections,
# and weights is contained in the histOps dictionary
def makeNDhisto(year, proc, histName, histOps, outfile, tree, isData):

    # To efficiently TTree->Draw(), we will only "activate"
    # necessary branches. So first, disable all branches
    tree.SetBranchStatus("*", 0)

    selection = histOps["selection"]
    variable  = histOps["variable"]
    weight    = histOps["weight"]
    
    # Make one big string to extract all relevant branch names from
    concatStr = selection + "," + variable
    if not isData:
        concatStr = concatStr + "," + weight

    # The idea is to turn the concatStr into a comma separated list of branches
    # Parenthesis can simply be removed, operators are simply replace with a comma
    # After all replacements, the string is split on the comma and filtered for empty strings
    expressions = re.findall(r'(\w+::\w+)', concatStr)
    for exp in expressions:
        concatStr = concatStr.replace(exp, "")
    concatStr = re.sub("[()]", "", concatStr)
    for replaceStr in ["&&", "||", "==", "<=", ">=", ">", "<", "*", ".", "/", "+", "-", ":"]:
        concatStr = concatStr.replace(replaceStr, ",")
    branches = list(filter(bool, concatStr.split(",")))

    # Here, the branches list will be names of branches and strings of digits
    # The digits are residual cut expressions like NGoodJets_pt30>=7 ==> "NGoodJets_pt30", "7"
    # So if a supposed branch name can be turned into an int, then it is not a legit branch name
    for branch in branches:
        try:
            int(branch)
            continue
        except:
            tree.SetBranchStatus(branch, 1)

    is2D = False
    if re.search(r'^[^:]*:[^:]*$', variable):
        is2D = True

    outfile.cd()

    htemp = None
    if not is2D:
        temph = ROOT.TH1F(histName, "", histOps["xbins"], histOps["xmin"], histOps["xmax"])
    else:
        temph = ROOT.TH2F(histName, "", histOps["xbins"], histOps["xmin"], histOps["xmax"], histOps["ybins"], histOps["ymin"], histOps["ymax"])

    # For MC, we multiply the selection string by our chosen weight in order
    # to fill the histogram with an event's corresponding weight
    drawExpression = "%s>>%s"%(variable, histName)
    if isData:
        tree.Draw(drawExpression, selection)
    else:
        tree.Draw(drawExpression, "(%s)*(%s)"%(weight,selection))
      
    temph = ROOT.gDirectory.Get(histName)
    temph.Sumw2()

    temph.Write(histName, ROOT.TObject.kOverwrite)

# Main function that a given pool process runs, the input TTree is opened
# and the list of requested histograms are drawn to the output ROOT file
def processFile(outputDir, inputDir, year, proc, histograms, treeName):

    inFileName = "%s/%s_%s.root"%(inputDir, year, proc)
    infile = ROOT.TFile.Open(inFileName,  "READ"); infile.cd()
    if infile == None:
        print("Could not open input ROOT file \"%s\""%(inFileName))
        return

    tree = infile.Get(treeName)

    if tree == None:
        print("Could not get tree \"%s\" from input ROOT file \"%s\""%(treeName, inFileName))
        return

    outfile = ROOT.TFile.Open("%s/%s_%s.root"%(outputDir, year, proc), "RECREATE")

    isData = "Data" in proc
    for histName, histOps in histograms.items():
        makeNDhisto(year, proc, histName, histOps, outfile, tree, isData)

    outfile.Close()

if __name__ == "__main__":
    usage = "%miniTupleDrawer [options]"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("--inputDir",  dest="inputDir",  help="Path to ntuples",    required=True                )
    parser.add_argument("--outputDir", dest="outputDir", help="path to output",     default="MyTag",             )
    parser.add_argument("--tree",      dest="tree",      help="TTree name to draw", default="AnaSkim"            )
    parser.add_argument("--year",      dest="year",      help="which year",         default="Run2UL"             )
    parser.add_argument("--options",   dest="options",   help="options file",       default="miniTupleDrawer_aux")
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
    
    inputDir  = args.inputDir
    outputDir = args.outputDir
    treeName  = args.tree
    year      = args.year
    
    base = os.getenv("CMSSW_BASE")
    
    # The draw histograms and their host ROOT files are kept in the output
    # folder in the user's condor folder. This then makes running a plotter
    # on the output exactly like running on histogram output from an analyzer
    outputDir = "%s/src/Analyzer/Analyzer/test/condor/%s/"%(base,outputDir)
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    
    # For speed, histogramming for each specified physics process, e.g. TT, QCD
    # is run in a separate pool process. This is limited to 4 at a time to avoid abuse
    manager = mp.Manager()
    pool = mp.Pool(processes=min(4, len(processes)))
    
    # The processFile function is attached to each process
    for proc in processes:
        pool.apply_async(processFile, args=(outputDir, inputDir, year, proc, histograms, treeName))
    
    pool.close()
    pool.join()
