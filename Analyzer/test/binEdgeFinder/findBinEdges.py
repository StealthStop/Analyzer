import os
import ROOT
import math
import array
import copy
import sys
import numpy as np

ROOT.gROOT.SetBatch( True )

from optparse import OptionParser

parser = OptionParser()

parser.add_option( '--maxTtFrac', type = 'int', action = 'store',
                    dest = 'maxTtFrac',
                    default = 1000, help = 'Determine how finely grained you want to divide your SNN score (100 = .01)' )

parser.add_option( '--mainBkgdSyst', type = 'float', action = 'store',
                    dest = 'mainBkgdSyst',
                    default = 0.2, help = 'Systematic size for the main background - ttbar (set to 0.2 by default)' )
parser.add_option( '--subBkgdSyst', type = 'float', action = 'store',
                    dest = 'subBkgdSyst',
                    default = 0.2, help = 'Systematic size for the sub-dominant background - qcd (set to 0.2 by default)' )

parser.add_option( '--signalMStop', type = 'string', action = 'store',
                    dest = 'signalMStop',
                    default = '550', help = 'Mass of the scalar top quark for the signal model to optimize by' )

parser.add_option( '--year', type = 'string', action = 'store',
                    dest = 'year',
                    default = '2016', help = 'Find bin edge for which year' )

parser.add_option( '--debug', action = 'store_true',
                    dest = 'debug',
                    default = 'False', help = 'Output print statements to track progress' )

( options, args )       =  parser.parse_args()

njetBins                = [ 7, 8, 9, 10, 11, 12 ] #where 12 is >= 12 jets

def main() :
    
    #Define location of miniTrees as input (treeDir) and location to output all plots (outputDir)
    treeDir             = "../../haddfiles/MakeMiniTrees/"
    outputDir           = "binEdgePlots/"+options.year+"/"

    if not os.path.exists( outputDir ) :
        print("Output directory does not exist - making directory",outputDir)
        os.makedirs( outputDir )

    #Create the dictionary for the main background, sub-dominant background, and the signal
    ttbarDict           = makeNJetArrays( treeDir, options.year+"_TT" )
    qcdDict             = makeNJetArrays( treeDir, options.year+"_QCD" )
    signalDict          = makeNJetArrays( treeDir, options.year+"_RPV_2t6j_mStop-"+options.signalMStop )
    
    #Create the ttbar efficiencies array such that they sum up to 100% for the four bins (using the successive bin method)
    fracArray           = []
    maxBins             = 4

    for nBins in xrange(2, maxBins+1) :
        bkgdFracArray       = createBkgdFrac( options.maxTtFrac, nBins, fracArray )
    
        #Create dummy variables
        maxSignificance     = 0.0
        finalBinEdges       = {}
        finalBkgdFrac       = [ 0.0, 0.0, 0.0, 0.0 ]
       
        debugValue          = 0.01
        
        #bkgdFracArray = [[ 0.560, 0.365, 0.051, 0.024 ]] #Bin edge background fractions determined for 80 fb-1 using 2016 simulation
    
        #Find the bin edges corresponding to the background fraction
        for bkgdFrac in bkgdFracArray :
            
            if options.debug :
                if debugValue != bkgdFrac[0] :
                    debugValue = bkgdFrac[0]
                    print bkgdFrac[0]
            binEdgesDict    = findBinEdges( ttbarDict, bkgdFrac )
            tempSig         = calcSignificance( binEdgesDict, bkgdFrac, ttbarDict, qcdDict, signalDict )
        
            if tempSig > maxSignificance :
                maxSignificance = tempSig
                finalBinEdges   = binEdgesDict
                finalBkgdFrac   = bkgdFrac

        fracArray           = finalBkgdFrac[1:]

    #Write out the bin edges as for DeepESM.cfg file
    writeBinEdgesCfg( outputDir, finalBinEdges )

def calcSignificance( binEdgesDict, bkgdFrac, mainBkgdDict, subBkgdDict, signalDict ) :

    #print "Calculating significance with bin edges for bkgd fraction",bkgdFrac
    
    #Create event yield per Njet bin per Snn bin dictionaries
    mainBkgdEventYield          = {}
    subBkgdEventYield           = {}
    signalEventYield            = {}
    
    for itJet in njetBins:
        mainBkgdEventYield[itJet]   = {}
        subBkgdEventYield[itJet]    = {}
        signalEventYield[itJet]     = {}
        
        mainBkgdArray               = zip( mainBkgdDict["snn"][itJet], mainBkgdDict["evtWght"][itJet], np.cumsum(mainBkgdDict["evtWght"][itJet]) )
        subBkgdArray                = zip( subBkgdDict["snn"][itJet], subBkgdDict["evtWght"][itJet], np.cumsum(subBkgdDict["evtWght"][itJet]) )
        signalArray                 = zip( signalDict["snn"][itJet], signalDict["evtWght"][itJet], np.cumsum(signalDict["evtWght"][itJet]) )
        
        itSnnBin                = 1
        for entry in xrange( 0, len(mainBkgdArray) ):
            if mainBkgdArray[entry][0] > binEdgesDict[itJet][itSnnBin] :
                if entry != 0 :
                    mainBkgdEventYield[itJet][itSnnBin] = mainBkgdArray[entry-1][2]
                else :
                    mainBkgdEventYield[itJet][itSnnBin] = mainBkgdArray[0][2] 
                itSnnBin += 1
                continue
        while itSnnBin != len(bkgdFrac):
            mainBkgdEventYield[itJet][itSnnBin]  = mainBkgdArray[len(mainBkgdArray)-1][2] 
            itSnnBin += 1
        mainBkgdEventYield[itJet][len(bkgdFrac)] = mainBkgdArray[len(mainBkgdArray)-1][2]

        itSnnBin                = 1
        for entry in xrange( 0, len(subBkgdArray) ):
            if subBkgdArray[entry][0] > binEdgesDict[itJet][itSnnBin] :
                if entry != 0 :
                    subBkgdEventYield[itJet][itSnnBin] = subBkgdArray[entry-1][2]
                else :
                    subBkgdEventYield[itJet][itSnnBin] = subBkgdArray[0][2] 
                itSnnBin += 1
                continue
        while itSnnBin != len(bkgdFrac):
            subBkgdEventYield[itJet][itSnnBin]  = subBkgdArray[len(subBkgdArray)-1][2]
            itSnnBin += 1
            
        subBkgdEventYield[itJet][len(bkgdFrac)] = subBkgdArray[len(subBkgdArray)-1][2]

        itSnnBin                = 1
        for entry in xrange( 0, len(signalArray) ):
            if signalArray[entry][0] > binEdgesDict[itJet][itSnnBin] :
                signalEventYield[itJet][itSnnBin] = signalArray[entry-1][2]
                if entry != 0 :
                    signalEventYield[itJet][itSnnBin] = signalArray[entry-1][2]
                else :
                    signalEventYield[itJet][itSnnBin] = signalArray[0][2] 
                itSnnBin += 1
                continue
        while itSnnBin != len(bkgdFrac):
            signalEventYield[itJet][itSnnBin]  = signalArray[len(signalArray)-1][2]
            itSnnBin += 1

        signalEventYield[itJet][len(bkgdFrac)] = signalArray[len(signalArray)-1][2]

        for i in xrange( 1, len( mainBkgdEventYield[itJet] )+1 ):
            for j in xrange( i+1, len( mainBkgdEventYield[itJet] )+1 ) :
                mainBkgdEventYield[itJet][j]  -= mainBkgdEventYield[itJet][i]
                subBkgdEventYield[itJet][j]   -= subBkgdEventYield[itJet][i]
                signalEventYield[itJet][j]    -= signalEventYield[itJet][i]

    significance                = {}
    totalSignificance           = 0.0
    for itJet in njetBins :
        significance[itJet]     = {}

#   Commented code is for calculating overall significance, but original calculation was done with each successive bin having the most significance.
#        for itSnnBin in xrange( 1, len(binEdgesDict[itJet])+1 ) :
#            significance[itJet][itSnnBin]        = calcPerBinSignificance( mainBkgdEventYield[itJet][itSnnBin], subBkgdEventYield[itJet][itSnnBin], signalEventYield[itJet][itSnnBin] )
#            totalSignificance                   += significance[itJet][itSnnBin]*significance[itJet][itSnnBin]
        
        significance[itJet]                 = calcPerBinSignificance( mainBkgdEventYield[itJet][2], subBkgdEventYield[itJet][2], signalEventYield[itJet][2] )
        totalSignificance                  += significance[itJet]

    return math.sqrt(totalSignificance)

def calcPerBinSignificance( mainBkgdEY, subBkgdEY, signalEY ):
    
    numerator                   = signalEY
    denominator                 = mainBkgdEY + subBkgdEY + options.mainBkgdSyst*options.mainBkgdSyst*mainBkgdEY*mainBkgdEY + options.subBkgdSyst*options.subBkgdSyst*subBkgdEY*subBkgdEY 
    if denominator <= 0 :
        print "You have a bin with no background!"
        return 0
    else :
        return numerator/math.sqrt(denominator)


def findBinEdges( mainBkgdDict, bkgdFrac ) :

    #print "Calculating bin edges for background fraction...",bkgdFrac 
    
    tempBinEdgesDict    = {}
    cumBkgdFrac         = np.cumsum( bkgdFrac )
    
    for itJet in njetBins :
        tempBinEdgesDict[itJet] = {}
        mainBkgdArray           = zip( mainBkgdDict["snn"][itJet], mainBkgdDict["evtWght"][itJet], np.cumsum(mainBkgdDict["evtWght"][itJet])/float(sum(mainBkgdDict["evtWght"][itJet])) )

        snnBin                  = 1
        for binIter in xrange(0, len(mainBkgdArray) ):
            if mainBkgdArray[binIter][2] > cumBkgdFrac[snnBin-1] :
                tempBinEdgesDict[itJet][snnBin] = float(math.floor(mainBkgdArray[binIter-1][0]*1000))/float(1000.0)
                snnBin          += 1

        tempBinEdgesDict[itJet][len(cumBkgdFrac)] = 1.0

    return tempBinEdgesDict

def makeNJetArrays( treeDir, fileName ) :
    
    #Get arrays from selected miniTree 
    print "Making njet arrays for this file: "+treeDir+fileName+".root"
    myFile              = ROOT.TFile( treeDir+fileName+".root" )
    myTree              = myFile.Get( "myMiniTree" )

    outDict             = {}
    outDict["snn"]      = {} #Output a dictionary with all the SNN scores
    outDict["evtWght"]  = {} #Output a dictionary with all the event weights
    
    for itJet in njetBins :
        outDict["snn"][itJet]       = [] #For each nJet, create an array
        outDict["evtWght"][itJet]   = []

    #lumi                            = myTree.Lumi
    lumi                            = 80.0
    for entry in xrange( 0, myTree.GetEntries() ):
        myTree.GetEntry( entry ) 
        
        if myTree.NGoodJets_pt30 <= 6 : #Depending on which tree you use, you may need to remove events with less than 6 jets
            continue
        
        if myTree.NGoodJets_pt30 <= 12 : #Important when making the njets = 12 bin inclusive
            outDict["snn"][myTree.NGoodJets_pt30].append( myTree.deepESM_val )
            outDict["evtWght"][myTree.NGoodJets_pt30].append( myTree.totalEventWeight * lumi )
        else :
            outDict["snn"][12].append( myTree.deepESM_val )
            outDict["evtWght"][12].append( myTree.totalEventWeight * lumi )

    myFile.Close()

    for itJet in njetBins :
        outDict["snn"][itJet], outDict["evtWght"][itJet] = zip(*sorted(zip(outDict["snn"][itJet], outDict["evtWght"][itJet]))) #Sort both arrays by the Snn score
    
    return outDict

def createBkgdFrac( maxTtFrac, nBins, oldBkgdFrac ):
    if nBins > 2 and len(oldBkgdFrac) != nBins-2 :
        print "There's a bug in your code!"
        raise RuntimeError

    bkgdFracArray       = []
    bkgdFracAccounted   = sum( oldBkgdFrac )

    for i in xrange( 1, maxTtFrac ):
        if i + int(bkgdFracAccounted*maxTtFrac) >= maxTtFrac :
            continue
        
        tempFracArray = [ float(i)/float(maxTtFrac), float(maxTtFrac - i - int(bkgdFracAccounted*maxTtFrac))/float(maxTtFrac) ]
        for element in oldBkgdFrac :
            tempFracArray.append( element )
        bkgdFracArray.append( tempFracArray )
    
    print "Bkgd frac array has",len(bkgdFracArray),"entries"

    return bkgdFracArray

def writeBinEdgesCfg( outputDir, binEdgesDict ) :
    
    count = 0 
    for itJet in njetBins:
        print '    binEdges['+str(count)+'] = 0.00000000'
        count += 1
        for itCut in range( 1, len( binEdgesDict[itJet] ) ):
            print '    binEdges['+str(count)+'] = '+str(binEdgesDict[itJet][itCut])
            count += 1
        print '    binEdges['+str(count)+'] = 1.00000000'
        count += 1

if __name__ == '__main__' :
    main()
