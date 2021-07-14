#!/bin/python
import ROOT
import copy
import os.path
import array
import binEdges as bE

ROOT.gROOT.SetBatch( True )

from optparse import OptionParser

parser = OptionParser()
parser.add_option('--log', action='store_true',
                      dest='log',
                      default = False, help='Set log Y')

parser.add_option('--test', action='store_true',
                      dest='test',
                      default = False, help='Flag for debug mode' )
parser.add_option('--blind', action='store_true',
                      dest='blind',
                      default = False, help='Look only at NJets <= 10')
parser.add_option('--snn', action='store_false',
                      dest='snn',
                      default = True, help='Default true, use this flag to turn off making histograms divided by bins in snn')


(options, args) = parser.parse_args()

def main() :

    # It takes a long time to run each year, so it may be better to run each year individually
    # yearList       = [ "2016", "2017", "2018pre", "2018post" ] 
    yearList        = [ "2018pre", "2018post" ]

    # Define the input directory for the mini trees in the SR
    inputDir        = "../haddfiles/MakeMiniTrees/"
    if not os.path.exists(inputDir):
        print("Directory with all the input mini trees do not exist")

    outputDir       = "htExtraRootFiles/"

    if not os.path.exists(outputDir) :
        print("Output directory does not exist, so making output directory",outputDir)
        os.makedirs(outputDir)

    # Set DefaultSumw2()
    ROOT.TH1.SetDefaultSumw2()

    # Bkg Input files
    fileListBkgd    = [ "TT","QCD","TTX", "Other" ]
    
    # Data Input files
    fileListData    = [ "Data" ]
    
    snnBinList      = [ "_D1", "_D2", "_D3", "_D4" ]
    njetList        = [ 7, 8, 9, 10, 11, 12 ]
    htWghtVarList   = [ "noHTWght", "htWght", "htWghtFlat2000", "htWghtNJet7" ]
    # Define the four types of event weights:
    #   1. Applying no HT weight
    #   2. Applying the full HT weight
    #   3. Applying the modified HT weight with the weight kept constant after HT = 2000 GeV
    #   4. Applying the modified HT weight with the njets dependence (using only njet = 7 weights)
   
    for yearTag in yearList: 
        
        # Get the bin edges from the designed binEdges.py file
        binEdges        = bE.binEdgeDictReturn( "SRBE", yearTag )
    
        # A lot of histograms are needed to be mad. The categories are named as such:
        #   - the first variable listed is the variable on the x-axis
        #   - snnBin means that there will be a histogram for each Snn bin (as opposed to all four snn bins in one histogram without this moniker)
        #   - perNjets means that there will be a histogram for each njets
        # For each category, there will be an additional four histograms, each with the different HT weights as mentioned above

        for filename in fileListBkgd:
            histoDict                                   = {}
            histoDict["HT"]                             = {}
            histoDict["HT_snnBin"]                      = {}
            
            histoDict["njets"]                          = {}
            histoDict["njets_snnBin"]                   = {}
            
            histoDict["HT_perNjets"]                    = {}
            histoDict["HT_perNjets_snnBin"]             = {}
    
            histoDict["snn"]                            = {}
    
            for i in range( len ( htWghtVarList ) ):
                h_HT_All            = ROOT.TH1F( "h_HT_"+htWghtVarList[i]+"_Allj", "h_HT_"+htWghtVarList[i]+"_Allj", 60, 0.0, 3000)

                h_njets_All         = ROOT.TH1F( "h_njets_"+htWghtVarList[i]+"_Allj", "h_njets_"+htWghtVarList[i]+"_Allj", 6, 7, 13)

                h_snn_All           = ROOT.TH1F( "h_snn_"+htWghtVarList[i]+"_Allj", "h_snn_"+htWghtVarList[i]+"_Allj", 20, 0.0, 1.0 )

                histoDict["HT"][ htWghtVarList[i] ] = copy.deepcopy( h_HT_All ) 
                histoDict["njets"][ htWghtVarList[i] ] = copy.deepcopy( h_njets_All ) 
                histoDict["snn"][ htWghtVarList[i] ] = copy.deepcopy( h_snn_All )
    
                histoDict["HT_snnBin"][ htWghtVarList[i] ] = {}
                histoDict["njets_snnBin"][ htWghtVarList[i] ] = {}
                
                for snnbin in snnBinList :
    
                    h_HT_snn_All    = ROOT.TH1F( "h_HT_"+htWghtVarList[i]+"_Allj"+snnbin, "h_HT_"+htWghtVarList[i]+"_Allj"+snnbin, 60, 0.0, 3000)       
                    h_njets_snn_All = ROOT.TH1F( "h_njets_"+htWghtVarList[i]+"_Allj"+snnbin, "h_njets_"+htWghtVarList[i]+"_Allj"+snnbin, 6, 7, 13)       
                
                    histoDict["HT_snnBin"][ htWghtVarList[i] ][ snnbin ] = copy.deepcopy( h_HT_snn_All ) 
                    histoDict["njets_snnBin"][ htWghtVarList[i] ][ snnbin ] = copy.deepcopy( h_njets_snn_All ) 
            
                histoDict["HT_perNjets"][ htWghtVarList[i] ] = {}
                histoDict["HT_perNjets_snnBin"][ htWghtVarList[i] ] = {}
    
                for j in xrange(0, len(njetList)):
                    
                    h_HT             = ROOT.TH1F( "h_HT_"+htWghtVarList[i]+"_"+str(njetList[j])+"j", "h_HT_"+str(njetList[j])+"j", 60, 0.0, 3000)       
                    histoDict["HT_perNjets"][ htWghtVarList[i] ][ njetList[j] ] = copy.deepcopy( h_HT )
    
                    histoDict["HT_perNjets_snnBin"][ htWghtVarList[i] ][ njetList[j] ] = {}
                    
                    for snnbin in snnBinList :
                        h_HT_snn     = ROOT.TH1F( "h_HT_"+htWghtVarList[i]+"_"+str(njetList[j])+"j"+snnbin, "h_HT_"+str(njetList[j])+"j"+snnbin, 60, 0.0, 3000)       
                    
                        histoDict["HT_perNjets_snnBin"][ htWghtVarList[i] ][ njetList[j] ][ snnbin ] = copy.deepcopy( h_HT_snn )
                    
            inputFile            = ROOT.TFile.Open(inputDir+yearTag+"_"+filename+".root")
            myTree               = inputFile.Get("myMiniTree")
    
            for entry in xrange( 0, myTree.GetEntries() ) : 
                myTree.GetEntry(entry)
                
                noHTweight          = myTree.totalEventWeight * myTree.Lumi / myTree.htDerivedweight
                HTweight            = myTree.totalEventWeight * myTree.Lumi
                HTweightFlat2000    = noHTweight * myTree.htDerivedweightFlat2000
                HTweightNJet7       = noHTweight * myTree.htDerivedweightNJet7
    
                htWeights           = [ noHTweight, HTweight, HTweightFlat2000, HTweightNJet7 ]
                
                for i in range( len( htWghtVarList ) ) :
                    
                    histoDict["HT"][ htWghtVarList[i] ].Fill( myTree.HT_trigger_pt30, htWeights[i] ) 
                    histoDict["njets"][ htWghtVarList[i] ].Fill( myTree.NGoodJets_pt30, htWeights[i] ) 
                    if myTree.NGoodJets_pt30 == 7 :
                        histoDict["snn"][ htWghtVarList[i] ].Fill( myTree.deepESM_val, htWeights[i] )
    
                    if options.snn :
                        if myTree.NGoodJets_pt30 > 12 :
                            for itBinEdge in xrange( 0, len(binEdges[12]) ) :
                                if myTree.deepESM_val < binEdges[12][itBinEdge] :
                                    histoDict["HT_snnBin"][ htWghtVarList[i] ][ snnBinList[itBinEdge] ].Fill( myTree.HT_trigger_pt30, htWeights[i] )
                                    histoDict["njets_snnBin"][ htWghtVarList[i] ][ snnBinList[itBinEdge] ].Fill( 12, htWeights[i] )
                                    break
                        else :
                            for itBinEdge in xrange( 0, len(binEdges[myTree.NGoodJets_pt30]) ) :
                                if myTree.deepESM_val < binEdges[myTree.NGoodJets_pt30][itBinEdge] :
                                    histoDict["HT_snnBin"][ htWghtVarList[i] ][ snnBinList[itBinEdge] ].Fill( myTree.HT_trigger_pt30, htWeights[i] )
                                    histoDict["njets_snnBin"][ htWghtVarList[i] ][ snnBinList[itBinEdge] ].Fill( myTree.NGoodJets_pt30, htWeights[i] )
                                    break
    
                    for j in range( len( njetList ) ):
                        if njetList[j] != 12: 
                            if myTree.NGoodJets_pt30 != njetList[j]:
                                continue
                        else :
                            if myTree.NGoodJets_pt30 <= 12 :
                                continue
                        
                        histoDict["HT_perNjets"][ htWghtVarList[i] ][ njetList[j] ].Fill( myTree.HT_trigger_pt30, htWeights[i] )
                        if options.snn:
                            for itBinEdge in xrange( 0, len(binEdges[njetList[j]]) ) :
                                if myTree.deepESM_val < binEdges[njetList[j]][itBinEdge] :
                                    histoDict["HT_perNjets_snnBin"][ htWghtVarList[i] ][ njetList[j] ][ snnBinList[itBinEdge] ].Fill( myTree.HT_trigger_pt30, htWeights[i] )
                                    break
            
            myOutputFile = ROOT.TFile( outputDir+yearTag+"_HTsystematic_"+filename+".root", "RECREATE")
    
            for htWghtTag in htWghtVarList:
                histoDict["HT"][ htWghtTag ].Write()
                histoDict["njets"][ htWghtTag ].Write()
                histoDict["snn"][ htWghtTag ].Write()
    
                for snn in snnBinList:
                    histoDict["HT_snnBin"][ htWghtTag ][ snn ].Write()
                    histoDict["njets_snnBin"][ htWghtTag ][ snn ].Write()
    
                for njet in njetList:
                    histoDict["HT_perNjets"][ htWghtTag ][ njet ].Write()
                    for snn in snnBinList:
                        histoDict["HT_perNjets_snnBin"][ htWghtTag ][ njet ][ snn ].Write()
                
            myOutputFile.Close()
        
                            
        for filename in fileListData:
            histoDict["HT"]           = {}
            histoDict["HT_snnBin"]       = {}
            
            histoDict["njets"]        = {}
            histoDict["njets_snnBin"]    = {}
            
            histoDict["HT_perNjets"]               = {}
            histoDict["HT_perNjets_snnBin"]           = {}
    
            histoDict["snn"]          = {}
    
            for i in range( len ( htWghtVarList ) ):
                h_HT_All            = ROOT.TH1F( "h_HT_"+htWghtVarList[i]+"_Allj", "h_HT_"+htWghtVarList[i]+"_Allj", 60, 0.0, 3000)      
                h_njets_All         = ROOT.TH1F( "h_njets_"+htWghtVarList[i]+"_Allj", "h_njets_"+htWghtVarList[i]+"_Allj", 6, 7, 13)      
                h_snn_All           = ROOT.TH1F( "h_snn_"+htWghtVarList[i]+"_Allj", "h_snn_"+htWghtVarList[i]+"_Allj", 20, 0.0, 1.0 )
    
                histoDict["HT"][ htWghtVarList[i] ] = copy.deepcopy( h_HT_All ) 
                histoDict["njets"][ htWghtVarList[i] ] = copy.deepcopy( h_njets_All ) 
                histoDict["snn"][ htWghtVarList[i] ] = copy.deepcopy( h_snn_All )
    
                histoDict["HT_snnBin"][ htWghtVarList[i] ] = {}
                histoDict["njets_snnBin"][ htWghtVarList[i] ] = {}
                
                for snnbin in snnBinList :
    
                    h_HT_snn_All    = ROOT.TH1F( "h_HT_"+htWghtVarList[i]+"_Allj"+snnbin, "h_HT_"+htWghtVarList[i]+"_Allj"+snnbin, 60, 0.0, 3000)       
                    h_njets_snn_All = ROOT.TH1F( "h_njets_"+htWghtVarList[i]+"_Allj"+snnbin, "h_njets_"+htWghtVarList[i]+"_Allj"+snnbin, 6, 7, 13)       
                
                    histoDict["HT_snnBin"][ htWghtVarList[i] ][ snnbin ] = copy.deepcopy( h_HT_snn_All ) 
                    histoDict["njets_snnBin"][ htWghtVarList[i] ][ snnbin ] = copy.deepcopy( h_njets_snn_All ) 
            
                histoDict["HT_perNjets"][ htWghtVarList[i] ] = {}
                histoDict["HT_perNjets_snnBin"][ htWghtVarList[i] ] = {}
    
                for j in xrange(0, len(njetList)):
                    
                    h_HT             = ROOT.TH1F( "h_HT_"+htWghtVarList[i]+"_"+str(njetList[j])+"j", "h_HT_"+str(njetList[j])+"j", 60, 0.0, 3000)       
                    histoDict["HT_perNjets"][ htWghtVarList[i] ][ njetList[j] ] = copy.deepcopy( h_HT )
    
                    histoDict["HT_perNjets_snnBin"][ htWghtVarList[i] ][ njetList[j] ] = {}
                    
                    for snnbin in snnBinList :
                        h_HT_snn     = ROOT.TH1F( "h_HT_"+htWghtVarList[i]+"_"+str(njetList[j])+"j"+snnbin, "h_HT_"+str(njetList[j])+"j"+snnbin, 60, 0.0, 3000)       
                    
                        histoDict["HT_perNjets_snnBin"][ htWghtVarList[i] ][ njetList[j] ][ snnbin ] = copy.deepcopy( h_HT_snn )
                    
            file                 = ROOT.TFile.Open(inputDir+yearTag+"_"+filename+".root")
            myTree               = file.Get("myMiniTree")
    
            for entry in xrange( 0, myTree.GetEntries() ) :        
                myTree.GetEntry(entry)
                
                noHTweight          = 1.0
                HTweight            = 1.0
                HTweightFlat2000    = 1.0
                HTweightNJet7       = 1.0
    
                htWeights           = [ noHTweight, HTweight, HTweightFlat2000, HTweightNJet7 ]
    
                for i in range( len( htWghtVarList ) ) :
                    
                    histoDict["HT"][ htWghtVarList[i] ].Fill( myTree.HT_trigger_pt30, htWeights[i] ) 
                    histoDict["njets"][ htWghtVarList[i] ].Fill( myTree.NGoodJets_pt30, htWeights[i] ) 
                    if myTree.NGoodJets_pt30 == 7 :
                        histoDict["snn"][ htWghtVarList[i] ].Fill( myTree.deepESM_val, htWeights[i] ) 
                    
                    if options.snn :
                        if myTree.NGoodJets_pt30 > 12 :
                            for itBinEdge in xrange( 0, len(binEdges[12]) ) :
                                if myTree.deepESM_val < binEdges[12][itBinEdge] :
                                    histoDict["HT_snnBin"][ htWghtVarList[i] ][ snnBinList[itBinEdge] ].Fill( myTree.HT_trigger_pt30, htWeights[i] )
                                    histoDict["njets_snnBin"][ htWghtVarList[i] ][ snnBinList[itBinEdge] ].Fill( myTree.HT_trigger_pt30, htWeights[i] )
                                    break
                        else :
                            for itBinEdge in xrange( 0, len(binEdges[myTree.NGoodJets_pt30]) ) :
                                if myTree.deepESM_val < binEdges[myTree.NGoodJets_pt30][itBinEdge] :
                                    histoDict["HT_snnBin"][ htWghtVarList[i] ][ snnBinList[itBinEdge] ].Fill( myTree.HT_trigger_pt30, htWeights[i] )
                                    histoDict["njets_snnBin"][ htWghtVarList[i] ][ snnBinList[itBinEdge] ].Fill( myTree.NGoodJets_pt30, htWeights[i] )
                                    break
    
                    for j in range( len( njetList ) ):
                        if njetList[j] != 12: 
                            if myTree.NGoodJets_pt30 != njetList[j]:
                                continue
                        else :
                            if myTree.NGoodJets_pt30 <= 12 :
                                continue
                        
                        histoDict["HT_perNjets"][ htWghtVarList[i] ][ njetList[j] ].Fill( myTree.HT_trigger_pt30, htWeights[i] )
                        if options.snn:
                            for itBinEdge in xrange( 0, len(binEdges[njetList[j]]) ) :
                                if myTree.deepESM_val < binEdges[njetList[j]][itBinEdge] :
                                    histoDict["HT_perNjets_snnBin"][ htWghtVarList[i] ][ njetList[j] ][ snnBinList[itBinEdge] ].Fill( myTree.HT_trigger_pt30, htWeights[i] )
                                    break
            
            myOutputFile = ROOT.TFile( outputDir+yearTag+"_HTsystematic_"+filename+".root", "RECREATE")
    
            for htWghtTag in htWghtVarList:
    
                histoDict["HT"][ htWghtTag ].Write()
                histoDict["njets"][ htWghtTag ].Write()
                histoDict["snn"][ htWghtTag ].Write()
    
                for snn in snnBinList:
                    histoDict["HT_snnBin"][ htWghtTag ][ snn ].Write()
                    histoDict["njets_snnBin"][ htWghtTag ][ snn ].Write()
    
                for njet in njetList:
                    histoDict["HT_perNjets"][ htWghtTag ][ njet ].Write()
                    for snn in snnBinList:
                        histoDict["HT_perNjets_snnBin"][ htWghtTag ][ njet ][ snn ].Write()
                
            myOutputFile.Close()

if __name__ == '__main__' :
    main()
