#This script takes the mini tree outputs for the MC and data in both the signal and control region and makes histograms of njets in different categories for comparison.
#Comparisons include:
#   - Using CR derived bin edges (from SingleMuon dataset) and using SR derived bin edges (from TT MC)
#   - Looking at events with just muons in the signal region, just electrons in the signal region, just muons in the control region, or both leptons in the signal region
#   - Between the different time eras: 2016, 2017, 2018 pre-HEM and 2018 post-HEM

## Last edited by Kelvin Mei - Dec 30th, 2020

#!/bin/python
import ROOT
import copy
import os
import sys
import array

ROOT.gROOT.SetBatch( True )

#Location of binEdges.py
#sys.path.insert(1, '/Users/Kmei/SUS-19-004/src/') ###NOTE:Location of your binEdges.py file

import binEdges as bE

def main() :

    #Set the DefaultSumw2 for good error propagation
    ROOT.TH1.SetDefaultSumw2()
    ROOT.TH2.SetDefaultSumw2()
    
    # MC/Data Input Files - make sure you have [2016,2017,2018pre,2018post]_Other.root where Other is the histogram add (ROOT command hadd) of all the other backgrounds (you can break it down further if you want, but for this file, there really isn't a need. One could theoretically add in TT too, but there were some TT CR MC studies that were conducted).
    ### NOTE: For this script to work, copy 201*_Data_SingleMuon.root to 201*_Data.root for the CR case (and h-add 201*Data_SingleMuon/Electron.root to 201*_Data.root for SR)
    categoryList      = [ "QCD", "TT", "Other", "Data" ]
    
    # Output directory check and production
    outputDir       = "njetRootHistos"
    
    if not os.path.exists(outputDir):
        print "Output directory does not exist - making directory",outputDir
        os.makedirs(outputDir)
    
    #Make histograms for each year, MVA bin, region, and type of lepton
    yearList        = [ "2016", "2017", "2018pre", "2018post" ]     # year of data - also prefix of file names
    regionList      = [ "CR", "SR" ]                                # pass CR baseline or pass SR baseline
    binEdgeList     = [ "CRBE", "SRBE" ]                            # using control region derived bin edges (CRBE) or signal region bin edges (SRBE)
    leptonList      = [ "el", "mu", "both" ]                        # which lepton category (for comparisons in the SR)
    mvaBinList      = [ "D1", "D2", "D3", "D4", "All" ]             # which MVA bin (for comparison plots)

    for category in categoryList:                                   #A separate file will be made for each MC breakdown

        # Define dictionary with all the histograms and all the histograms
        njetHistoDict   = {}
    
        for region in regionList :
            njetHistoDict[region] = {}
            
            for year in yearList :
                njetHistoDict[region][year] = {}
    
                for binEdge in binEdgeList:
    
                    # Notable exception: we do not look at signal region using control region derived bin edges, but one may wonder what that looks like (i.e. if QCD and TT have similar enough shapes?)
                    if binEdge == "CRBE" and region == "SR" :
                        continue

                    njetHistoDict[region][year][binEdge] = {}
                
                    for lepton in leptonList :
                        
                        # Single Muon data set - Ignore anything with electrons in the control region (of course)
                        if lepton == "el" and region == "CR" :
                            continue

                        if lepton == "both" and region == "CR" :
                            continue

                        else :
                            njetHistoDict[region][year][binEdge][lepton] = {}
        
                            for mvaBin in mvaBinList :
                                njetHistoDict[region][year][binEdge][lepton][mvaBin] = ROOT.TH1D( "h_njets_"+region+"_"+year+"_"+binEdge+"_"+lepton+"_"+mvaBin, "h_njets_"+region+"_"+year+"_"+binEdge+"_"+lepton+"_"+mvaBin, 6, 7, 13 )
    
        #Get the input file mini trees

        for region in regionList :
            inputDir    = ""
            if region == "CR" :
                inputDir = "../haddfiles/MakeMiniTrees_NIM/"            ###NOTE: Change this location to the location of your Mini Tress that pass the CR baseline
            if region == "SR" :
                inputDir    = "../haddfiles/MakeMiniTrees/"        ###NOTE: Change this location to the location of your Mini Trees that pass the SR baseline (mostly for QCD and TT for secondary comparisons)

            if not os.path.exists(inputDir):
                print "Error: you do not have the correct path/directory to the miniTrees for the ",region

            for year in yearList :
            
                for binEdge in binEdgeList :
                    
                    if binEdge == "CRBE" and region == "SR" :
                        continue
                    
                    print "Initializing and making histograms for",category,"passing the",region,"baseline in the year",year,"using bin edges derived from the",binEdge[:2]
                    
                    binEdges        = bE.binEdgeDictReturn( binEdge, year )
                    inputFile       = ROOT.TFile.Open( inputDir+year+"_"+category+".root")
                    myTree          = inputFile.Get( "myMiniTree" )     ###NOTE: Default mini tree name - if you changed the mini tree name, this has to change too

                    for entry in xrange( 0, myTree.GetEntries() ) : 
                        myTree.GetEntry(entry)
                        eventWeight = 1.0
                        
                        if region == "CR" :

                            if myTree.NNonIsoMuonJets_pt30 < 7: #Legacy requirement where the mini trees had events with less than 7 jets to see trends in the bulk
                                continue

                            if category != "Data" :
                                eventweight = myTree.totalEventWeightNIM*myTree.Lumi ###NOTE: Implementation of the CR totalEventWeightNIM did not have lumi factored in
                            else : 
                                eventweight = 1.0
                        
                            njetHistoDict[region][year][binEdge]["mu"]["All"].Fill( myTree.NNonIsoMuonJets_pt30, eventweight ) 
                            
                            if myTree.NNonIsoMuonJets_pt30 < 13 : #For all events where the njets = the bin number
                            
                                for itBinEdge in xrange( 0, len( binEdges[ myTree.NNonIsoMuonJets_pt30 ] ) ) :
                                    if myTree.deepESM_valNonIsoMuon < binEdges[ myTree.NNonIsoMuonJets_pt30 ][ itBinEdge ] :
                                        njetHistoDict[region][year][binEdge]["mu"][mvaBinList[itBinEdge]].Fill( myTree.NNonIsoMuonJets_pt30, eventweight )
                                        break

                            else :                          #Loop to make the last bin fully inclusive
                                for itBinEdge in xrange( 0, len( binEdges[ 12 ] ) ):
                                    if myTree.deepESM_valNonIsoMuon < binEdges[ 12 ][ itBinEdge ] :
                                        njetHistoDict[region][year][binEdge]["mu"][mvaBinList[itBinEdge]].Fill( 12, eventweight )
                                        break
                        
                        elif region == "SR" :
                            
                            if myTree.NGoodJets_pt30 < 7: #Legacy requirement where the mini trees had events with less than 7 jets to see trends in the bulk
                                continue
                            
                            if category != "Data" :
                                eventweight = myTree.totalEventWeight*myTree.Lumi ###NOTE: Implementation of the SR totalEventWeight did not have lumi factored in
                            else :
                                eventweight = 1.0

                            njetHistoDict[region][year][binEdge]["both"]["All"].Fill( myTree.NGoodJets_pt30, eventweight )
                            
                            if myTree.NGoodMuons == 1:
                                njetHistoDict[region][year][binEdge]["mu"]["All"].Fill( myTree.NGoodJets_pt30, eventweight )
                            
                            if myTree.NGoodElectrons == 1:
                                njetHistoDict[region][year][binEdge]["el"]["All"].Fill( myTree.NGoodJets_pt30, eventweight )
                                
                            if myTree.NGoodJets_pt30 < 13 : #For all events where the njets = the bin number
                                for itBinEdge in xrange( 0, len( binEdges[ myTree.NGoodJets_pt30 ] ) ) :
                                    if myTree.deepESM_val < binEdges[ myTree.NGoodJets_pt30 ][ itBinEdge ] :
                                        njetHistoDict[region][year][binEdge]["both"][mvaBinList[itBinEdge]].Fill( myTree.NGoodJets_pt30, eventweight )
                                        if myTree.NGoodMuons == 1:
                                            njetHistoDict[region][year][binEdge]["mu"][mvaBinList[itBinEdge]].Fill( myTree.NGoodJets_pt30, eventweight )
                                        if myTree.NGoodElectrons == 1:
                                            njetHistoDict[region][year][binEdge]["el"][mvaBinList[itBinEdge]].Fill( myTree.NGoodJets_pt30, eventweight )
                                        break

                            else :                          #Lopp to make the last bin fully inclusive
                                for itBinEdge in xrange( 0, len( binEdges[ 12 ] ) ) :
                                    if myTree.deepESM_val < binEdges[ 12 ][ itBinEdge ] :
                                        njetHistoDict[region][year][binEdge]["both"][mvaBinList[itBinEdge]].Fill( 12, eventweight )
                                        if myTree.NGoodMuons == 1:
                                            njetHistoDict[region][year][binEdge]["mu"][mvaBinList[itBinEdge]].Fill( 12, eventweight )
                                        if myTree.NGoodElectrons == 1:
                                            njetHistoDict[region][year][binEdge]["el"][mvaBinList[itBinEdge]].Fill( 12, eventweight )
                                        break
                        
        
        outputFile = ROOT.TFile( outputDir+"/"+category+".root", "RECREATE")
        for region in regionList :
            for year in yearList :
                for binEdge in binEdgeList:
                    if binEdge == "CRBE" and region == "SR" :
                        continue
                    for lepton in leptonList :
                        if lepton == "el" and region == "CR" :
                            continue
                        if lepton == "both" and region == "CR" :
                            continue
                        else :
                            for mvaBin in mvaBinList :
                                njetHistoDict[region][year][binEdge][lepton][mvaBin].Write()
        
        outputFile.Close()

if __name__ == '__main__' :
    main( )
