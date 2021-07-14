# This script takes the mini Trees with the signal selection with the Njets requirement
# lowered to greater than and equal to 5 and creates files with the histograms needed
# to make the ht nominal scale factor

#!/bin/python

import math
import ROOT
import array
import random
import copy
import os

ROOT.gROOT.SetBatch(True)

from optparse import OptionParser

parser = OptionParser()

parser.add_option( '--data', action = 'store_true',
                    dest = 'data',
                    default = False, help = 'Use 2017 data file instead of 2016' )

(options, args) = parser.parse_args()

yearList        = [ '2016' ]#, '2017', '2018pre', '2018post' ]
njetList        = [ '5j', '6j', '7j', '8j' ]
fileList        = [ 'AllBG', 'Data_SingleLepton' ]

def main() :
        
    inputDir            = '../haddfiles/MakeMiniTrees/'

    if not os.path.exists( inputDir ) :
        print( "Directory with the input mini trees do not exist" )
        return 0

    outputDir           = 'mainHtSystematic'

    if not os.path.exists( outputDir ) :
        print( "Output directory does not exist. Creating directory, ",outputDir )
        os.makedirs( outputDir )
        
    for year in yearList :
        
        for fileType in fileList :
    
            # Make all the input histograms
            histogramDict       = {}
            
            inputFileName   = inputDir+year+"_"+fileType+".root"
            inputFile       = ROOT.TFile.Open( inputFileName )
            inputTree       = inputFile.Get("myMiniTree")
           
            for njet in njetList:
                tempHisto   = ROOT.TH1D( "h_ht_1l_"+njet+"_ge1b", "h_ht_1l_"+njet+"_ge1b", 300, 0.0, 3000 )
                histogramDict[njet] = copy.deepcopy( tempHisto )
            
            for entry in xrange( 0, inputTree.GetEntries() ) :
                
                inputTree.GetEntry( entry )

                weight      = 1.0
                if fileType == 'AllBG':
                    weight  = inputTree.totalEventWeight/inputTree.htDerivedweight*inputTree.Lumi
                for njet in njetList :
                    if str(inputTree.NGoodJets_pt30)+'j' == njet :
                        histogramDict[njet].Fill( inputTree.HT_trigger_pt30, weight )
                        break
    
            inputFile.Close()
        
            outputFileName          = year+"_htStudy_"+fileType+".root"
            outputFile              = ROOT.TFile( outputDir+"/"+outputFileName, "RECREATE" )
    
            for njet in njetList :
                histogramDict[njet].Write()
            outputFile.Close()

            del histogramDict

if __name__ == '__main__' :
    main()
