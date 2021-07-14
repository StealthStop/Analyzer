# All this script does is take the histograms produced from makeQcdCrComparisonHistos.py and 
# takes the data and subtract out the TT and Other MC component to look at what the "pure QCD" data looks like

## Last edited by Kelvin Mei on December 30th, 2020

#!/bin/python
import ROOT
import copy
import os
import array

ROOT.gROOT.SetBatch( True )

def main() :
    
    #Set the defaultsumw2
    ROOT.TH1.SetDefaultSumw2()
    ROOT.TH2.SetDefaultSumw2()
    
    #Make histograms for each year, type of bin edges, and snn bin - these should correspond to the names in makeQcdCrComparisonHistos.py
    yearList        = [ "2016" , "2017", "2018pre", "2018post" ]
    binEdgeList     = [ "CRBE", "SRBE" ]
    mvaBinList      = [ "D1", "D2", "D3", "D4", "All" ]

    #Other variables from makeQcdCRComparisonHistos that is specific to this file (only want to do the data subtraction for CR and with muons)
    region          = "CR"
    lepton          = "mu"

    inputDir                    = "njetRootHistos" ###NOTE: This is the default output directory for makeQcdCrComparisonHistos.py
    if not os.path.exists(inputDir):
        print "Input directory",inputDir,"does not exist."

    inputFileDataName           = inputDir+"/Data.root"
    inputFileMcSubtractName     = inputDir+"/Other.root"
    inputFileMcSubtract2Name    = inputDir+"/TT.root"

    outputFile = ROOT.TFile( inputDir+"/Data_QcdOnly.root", "RECREATE")

    for year in yearList :
        for binEdge in binEdgeList:
            for mvaBin in mvaBinList :
                dataQcdOnlyHist     = ROOT.TH1D()
                inputFileData       = ROOT.TFile.Open( inputFileDataName )
                dataAllHist         = inputFileData.Get( "h_njets_"+region+"_"+year+"_"+binEdge+"_"+lepton+"_"+mvaBin )
                
                dataQcdOnlyHist     = copy.deepcopy( dataAllHist.Clone() )
                inputFileData.Close()
                
                inputFileMC         = ROOT.TFile.Open( inputFileMcSubtractName )
                mcSubtractNameHist  = inputFileMC.Get( "h_njets_"+region+"_"+year+"_"+binEdge+"_"+lepton+"_"+mvaBin ) 
                dataQcdOnlyHist.Add( mcSubtractNameHist, -1 )
                inputFileMC.Close()
                
                inputFileMC2        = ROOT.TFile.Open( inputFileMcSubtract2Name )
                mcSubtract2NameHist = inputFileMC2.Get( "h_njets_"+region+"_"+year+"_"+binEdge+"_"+lepton+"_"+mvaBin ) 
                dataQcdOnlyHist.Add( mcSubtract2NameHist, -1 )
                inputFileMC2.Close()
                
                outputFile.cd()
                dataQcdOnlyHist.Write()
    
    outputFile.Close()

if __name__ == '__main__' :
    main()
