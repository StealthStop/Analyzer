# This script takes the output of MakeMiniTrees and generates the first step to the QCD control region systematic

# Input: Mini trees for both the control region MC, the control region data, and signal region MC
# Output: Histograms that have the Snn score distributions broken down by Njets for each MC/data sample

# Last edited: Kelvin Mei on December 30th, 2020

import math
import ROOT
import array
import random
import copy
import os

ROOT.gROOT.SetBatch(True)

import binEdges as bE

from optparse import OptionParser

parser = OptionParser()

#The signal option uses the bin edges derived in the TT bar MC that will also be used in the signal region. This is the option you will use if you want to make the QCD control region systematic
parser.add_option( '--signal', action = 'store_true',
                    dest = 'signal',
                    default = False, help = 'Use bin edges derived in signal TT region')

#The data option is used when taking the bin edges derived using binEdgeAnalyzer.py on the Single Muon data set. This was used for studies where you derive the bin edges on QCD and see how the fit function we used works.
parser.add_option( '--data', action = 'store_true',
                    dest = 'data',
                    default = False, help = 'Use bin edges derived with Single Muon data in QCD control region')

#Use the signal region baseline mini trees to make TT Snn histograms
parser.add_option( '--ttbar', action = 'store_true',
                    dest = 'ttbar',
                    default = False, help = 'Run over the signal region baseline TT mini trees to make the proper histograms')


(options, args) = parser.parse_args()

def main() :

    #Define variables needed for the analyzer here.
    ROOT.TH1.SetDefaultSumw2()

    #This is the location of the miniTrees
    inputDir            = '../haddfiles/MakeMiniTrees_NIM/'

    if options.ttbar :
        inputDir        = '../haddfiles/MakeMiniTrees/'

    if not os.path.exists(inputDir):
        print "Directory to miniTrees do not exist"

    outputDir           = 'snnPlots'

    if not os.path.exists(outputDir):
        print "Output directory does not exist, so making output directory", outputDir
        os.makedirs(outputDir)
    
    #Utilize this for the year option and to swap between different year files
    yearList            = [ "2016", "2017", "2018pre", "2018post" ]
    
    fileArray           = [ "QCD", "Data_SingleMuon" ]
    njetArray           = [ 7, 8, 9, 10, 11 ] 
    
    if options.ttbar :
        fileArray        = [ "TT" ]
    
    #Choose the correct version of the bin edges (can add more options for more different bin edges)
    if options.data and options.signal :
        print "Cannot have two options for bin edges"
        return 0

    if not options.data and not options.signal :
        print "Must choose a set of bin edges"
        return 0

    regionTag           = ""

    if options.data :
        regionTag       = "CRBE"
    if options.signal :
        regionTag       = "SRBE" 
    
    for yearTag in yearList :
        
        #Create a dictionary to hold onto all the histograms that will be used to write to an output file
        histogramDict                   = {}
        
        #Define bin edges here
        binEdges            = bE.binEdgeDictReturn( regionTag, yearTag )
        
        for mcFile in fileArray :
            #Make all the input Snn histograms
            inputFileName   = inputDir+yearTag+"_"+mcFile+".root"
            inputFile       = ROOT.TFile.Open( inputFileName )
            inputTree       = ROOT.TTree()
            inputTree       = inputFile.Get( "myMiniTree" )
               
            histogramDict[mcFile]             = {}
            histogramDict[mcFile]["All"]      = ROOT.TH1D( "h_snn_all", "h_snn_all", 7, -0.2, 1.2 )
                
            for njet in njetArray :
                histogramDict[mcFile][njet]   = ROOT.TH1D( "h_snn_"+str(njet)+"j", "h_snn_"+str(njet)+"j", 7, -0.2, 1.2 )
                
            for entry in xrange( 0, inputTree.GetEntries() ) :
                inputTree.GetEntry( entry )

                if options.ttbar :
                    deepESM_val                 = inputTree.deepESM_val
                    N_jets                      = inputTree.NGoodJets_pt30
                    eventWeight                 = inputTree.Lumi*inputTree.totalEventWeight
                else:
                    deepESM_val                 = inputTree.deepESM_valNonIsoMuon
                    N_jets                      = inputTree.NNonIsoMuonJets_pt30
                    if "Data" in mcFile :
                        eventWeight             = 1.0
                    else :
                        eventWeight             = inputTree.Lumi*inputTree.totalEventWeightNIM
    
                if options.signal :
                    if N_jets < 7 :
                        continue
                        
                    histogramDict[mcFile]["All"].Fill( deepESM_val, eventWeight)
            
                    if N_jets < 12 :
                        for njet in njetArray :
                            if N_jets == njet :
                                histogramDict[mcFile][njet].Fill( deepESM_val, eventWeight )
                                break
            
                    else:
                        histogramDict[mcFile][11].Fill( deepESM_val, eventWeight )
                    
            #Write the Snn histograms to output files
            outputFile  = ROOT.TFile( outputDir+"/"+yearTag+"_"+mcFile+"_SnnPerNJets.root", "RECREATE" )
            histogramDict[mcFile]["All"].Write()
            for njet in njetArray :
                histogramDict[mcFile][njet].Write()
            outputFile.Close()
            inputFile.Close()
        
if __name__ == '__main__' :
    main()
