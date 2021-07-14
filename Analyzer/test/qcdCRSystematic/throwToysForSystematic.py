# This is the script that will throw n (default n = 100 ) number of toys using
# an additional multiplicative weight derived from using the data from the CR 
# and the Powheg TT MC in the SR.

# Input: ROOT file output from makeRatioOfRatios.py
# Output: ROOT file output with all 100 toys and some diagnostic plots.

# Last updated January 12th, 2021 by Kelvin Mei

#!/bin/python
import ROOT
import copy
import os.path
import array
import numpy as np

ROOT.gROOT.SetBatch( True )

import binEdges as bE

def main():

    # Make sure this is set so errors are properly propagated
    ROOT.TH1.SetDefaultSumw2()
    ROOT.TH2.SetDefaultSumw2()

    # Define output directory for ROOT files with all toys
    outputDir           = 'toyPlots'

    if not os.path.exists(outputDir):
        print("Output directory does not exist, so making output directory", outputDir)
        os.makedirs(outputDir)

    # Define some useful variables here
    yearList                            = [ "2016" ] # Each year takes a bit of time to run, so I would recommend running each individual year separately.
    njetArray                           = [ 7, 8, 9, 10, 11, 12 ] 
    snnCutArray                         = []
    snnBinArray                         = [ "D1", "D2", "D3", "D4" ]
    seedNumber                          = 689300203 # This was a random number that was chosen from Google for reproducibility. You can change this number whenever
    
    # This is the location of the signal miniTrees (where the ttbar signal miniTree is)
    inputDir                            = "../haddfiles/MakeMiniTrees/"

    for year in yearList :
        print("Generating systematic for the year:",year)
        inputFile                       = ROOT.TFile.Open( inputDir+year+"_TT.root" )
        inputTree                       = inputFile.Get( "myMiniTree" )
    
        # This corresponds to the binning of the Snn histograms
        # CHECK THIS TO MAKE SURE THIS IS RIGHT IF THE SNN HISTOGRAM BINNING CHANGES (can also streamline)
        #   i.e. if there are 5 bins for the Snn histograms (which is how we do it), each bin is 0.2 of the Snn score, and there is an underflow and overflow bin.

        for i in range( 7 ):
            snnCutArray.append( float(i-1)*.20 )
       
        # Define number of toys that you want to throw here. The output file will have all
        #       of the D1-D4 histograms for each toy
        # For 2016, 100 toys takes about 30 - 60 minutes running locally on a 2015 Macbook
        nToys                               = 100
    
        # Get all the information for the ratio of ratios and place them in a dictionary for easy access. The SF dictionary will hold the central value and the Err dictionary will hold the error bar
        ratioOfRatiosSF                     = {}
        ratioOfRatiosErr                    = {}
        inputFileName                       = "ratioOfRatioPlots/"+year+"_ratioOfRatios.root"
        inputRatioOfRatiosFile              = ROOT.TFile.Open( inputFileName )
            
        for njet in njetArray :
            if njet > 11 : # Decision was made to do the snn plots only up to N_{jets} = 11 for statistics reasons
                continue
            ratioOfRatioHisto               = inputRatioOfRatiosFile.Get( "h_ratio_DataCR_over_TTSR_"+str(njet) )
            ratioOfRatiosSF[njet]           = {}
            ratioOfRatiosErr[njet]          = {}
            
            for binIter in range( ratioOfRatioHisto.GetNbinsX() ) :
               ratioOfRatiosSF[njet][snnCutArray[binIter]] = ratioOfRatioHisto.GetBinContent(binIter)
               ratioOfRatiosErr[njet][snnCutArray[binIter]] = ratioOfRatioHisto.GetBinError(binIter)

        inputRatioOfRatiosFile.Close()
    
        # Define two dictionaries to hold plots to validate procedure and to calculate the systematic. 
        # The 1D histogram dictionary will hold the mean and standard deviation values (just for sanity checks) of this new multiplicative factor
        # The 2D histogram dictionary will hold the plot showing the extra chosen multiplicative factor.
        
        histoDict                           = {}
        histo2dDict                         = {}
    
        # Make histograms and sanity plots here
        for njet in njetArray :
        
            # Histogram filled with all 100 toys with the Snn Score
            histoDict[njet]                 = ROOT.TH1D( "h_snn_100toys_"+str(njet)+"j", "h_snn_100toys_"+str(njet)+"j", 20, 0, 1.0 )
            
            histo2dDict[njet]               = ROOT.TH2D( "h_snn_evtWght_"+str(njet)+"j", "h_snn_evtWght_"+str(njet)+"j", 20, 0, 1.0, 200, 0.0, 2.0 )
        
        histoDict["std"]                    = {}
   
        # Make Njets plots keyed by Snn bin (i.e. njets distribution for Snn bin D1)
        for snn in snnBinArray :
            histoDict[snn]                  = ROOT.TH1D( "h_"+snn+"_qcdCR_meanValue", "h_"+snn+"_qcdCR_meanValue", 6, 7, 13 )
            histoDict["std"][snn]           = ROOT.TH1D( "h_"+snn+"_qcdCR_stdDev", "h_"+snn+"_qcdCR_stdDev", 6, 7, 13 )
        
        #Define bin edges based on the function. "SRBE" means derived from TT signal region.
        binEdges                            = bE.binEdgeDictReturn( "SRBE", year )
    
        #Define some more dictionaries for the second attempt to calculate the systematic
        #   totalEventArray holds all the final number of events in each njet bin and Snn bin (24 bins total)
        #   toyHistoDict holds all the D1-D4 njet histograms for each toy trial
    
        totalEventArray                         = {}
        toyHistoDict                            = {}
        for i in snnBinArray :
            totalEventArray[i]                  = {}
            for j in njetArray :
                totalEventArray[i][j]           = []
                for k in xrange( 0, nToys ) :
                    totalEventArray[i][j].append(0.0)
    
        for i in xrange( 0, nToys ) :
            toyHistoDict[i]                     = {}
            for snn in snnBinArray :
                toyHistoDict[i][snn]            = ROOT.TH1D( "h_"+snn+"_toy"+str(i), "h_"+snn+"_toy"+str(i), 6, 7, 13 )
    
        # Since these processes take a long time, it is good to print the number of entries
        print( "Number of entries in the TT file that will be looped over: "+str(inputTree.GetEntries())+" for "+str(nToys)+" toys" )
        
        for entry in xrange( 0, inputTree.GetEntries() ): 
            inputTree.GetEntry( entry )
            if entry % 10000 == 0:
                print( "Processing event number: ",entry)
    
            #Out of habit, initialize variables (and also to clear the previous value)
            scaleFactorMean                 = 0.0
            scaleFactorSigma                = 1.0
            snnBin                          = ""

            # Sanity safety check
            if inputTree.NGoodJets_pt30 < 7 :
                continue
   
            # Normal case code here - need a special case for events with 12 or more jets
            if inputTree.NGoodJets_pt30 < 12 : 

                # Match njet with the njets in the event
                for njet in njetArray :    
                    if inputTree.NGoodJets_pt30 == njet: 
                        for binEdgeIter in range( len( snnBinArray ) ) :
                            if inputTree.deepESM_val < binEdges[njet][binEdgeIter] :
                                snnBin          = snnBinArray[binEdgeIter]
                                break
                    
                        #Determine what the new event based weight is 
                        for snnBinIter in xrange( 0, len( snnCutArray ) ) :
                            if inputTree.deepESM_val < snnCutArray[snnBinIter] :
                                scaleFactorMean     = ratioOfRatiosSF[njet][snnCutArray[snnBinIter]]
                                scaleFactorSigma    = ratioOfRatiosErr[njet][snnCutArray[snnBinIter]]
                                break
    
                        for i in xrange( 0, nToys ) :
                            rNG                             = ROOT.TRandom3( seedNumber )
                            seedNumber                      += 1
                            newEventWeightSF                = rNG.Gaus( scaleFactorMean, scaleFactorSigma )
                            histo2dDict[njet].Fill( inputTree.deepESM_val, newEventWeightSF )
                            histoDict[njet].Fill( inputTree.deepESM_val, newEventWeightSF )
                            histoDict[snnBin].Fill( inputTree.NGoodJets_pt30, inputTree.Lumi*inputTree.totalEventWeight*newEventWeightSF )
                            totalEventArray[snnBin][inputTree.NGoodJets_pt30][i] += inputTree.Lumi*inputTree.totalEventWeight*newEventWeightSF
                            toyHistoDict[i][snnBin].Fill( inputTree.NGoodJets_pt30, inputTree.Lumi*inputTree.totalEventWeight*newEventWeightSF )
                            
            # For all events with more than 11 jets
            else :
                
                # Determine which Snn bin this event falls into - use njets >= 12 bin edges
                for binEdgeIter in range( len( snnBinArray ) ) :
                    if inputTree.deepESM_val < binEdges[12][binEdgeIter] :
                        snnBin          = snnBinArray[binEdgeIter]
                        break

                # Determine what the new event based weight is using njets >= 11 ratio of ratios histogram
                for snnBinIter in xrange( 0, len( snnCutArray ) ) :
                    if inputTree.deepESM_val < snnCutArray[snnBinIter] :
                        scaleFactorMean                 = ratioOfRatiosSF[11][snnCutArray[snnBinIter]]
                        scaleFactorSigma                = ratioOfRatiosErr[11][snnCutArray[snnBinIter]]
                        break
    
                for i in xrange( 0, nToys ) :
                    rNG                                 = ROOT.TRandom3( seedNumber )
                    seedNumber                          += 1
                    newEventWeightSF                    = rNG.Gaus( scaleFactorMean, scaleFactorSigma )
                    histo2dDict[12].Fill( inputTree.deepESM_val, newEventWeightSF )
                    histoDict[12].Fill( inputTree.deepESM_val, newEventWeightSF )
                    histoDict[snnBin].Fill( 12, inputTree.Lumi*inputTree.totalEventWeight*newEventWeightSF )
                    totalEventArray[snnBin][12][i] += inputTree.Lumi*inputTree.totalEventWeight*newEventWeightSF
                    toyHistoDict[i][snnBin].Fill( 12, inputTree.Lumi*inputTree.totalEventWeight*newEventWeightSF )
    
    
        # This is to make plots of the RMS of the 100 toys to determine the final RMS value of the 100 toys
        rmsHistoDict                        = {}
        for snn in snnBinArray:
            rmsHistoDict[snn]               = {}
            for njet in njetArray:
                y = np.array( totalEventArray[snn][njet] )
    
                print snn, histoDict[snn].GetBinLowEdge( njet-6 ), np.mean(y)
                rmsHistoDict[snn][njet]     = ROOT.TH1D( "h_"+snn+"_rms_"+str(njet)+"j", "h_"+snn+"_rms_"+str(njet)+"j", 100, -10.0, 10.0 )
                for iteration in totalEventArray[snn][njet] :
                    rmsHistoDict[snn][njet].Fill( iteration - np.mean(y) )
                histoDict["std"][snn].SetBinContent( njet-6, np.std(y) )
      
        # Write all all histograms to a ROOT file as well as all the diagnostic histograms
        outputFileName = outputDir+"/"+year+"_evtWghtSyst_allToys.root"
        outputFile  = ROOT.TFile.Open( outputFileName, "RECREATE" )
        
        for njet in njetArray :
            histoDict[njet].Write()
            histo2dDict[njet].Write()
        
        for snn in snnBinArray :
            histoDict[snn].Write()
            histoDict["std"][snn].Write()
        
            for njet in njetArray:
                rmsHistoDict[snn][njet].Write() 
        
        for i in xrange(0, nToys):
            for snn in snnBinArray :
                toyHistoDict[i][snn].Write()
    
        outputFile.Close()

if __name__ == '__main__' :
    main()
