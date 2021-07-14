import ROOT
import math
import array
import copy

ROOT.gROOT.SetBatch( True )

# The following code determines the optimal bin edges for the Stealth Stop Analysis
# For our signal region, the mainBkgd is TT, the subBkgd is QCD, and it is optimized for RPV 550
# For our control region, the mainBkgd is Data_SingleMuon

# If you are deriving bin edges for the first time, you should use the function: matchBkgdEffWithMaxSig
# If you are using already pre-determined bin edges ( from a different bin study ), you should use the function: usePreDeterminedCuts

#Allow options ot be passed into the function
from optparse import OptionParser

parser = OptionParser()

parser.add_option( '--debug', action = 'store_true', 
                    dest = 'debug',
                    default = False, help = 'Print Out All Debug Print Statements' )

parser.add_option( '--plotDesm', action = 'store_false',
                    dest = 'plotDesm',
                    default = True, help = 'Plot simple Deep ESM distributions for main and sub bkgd' )

parser.add_option( '--njets', action = 'store_false',
                    dest = 'njets',
                    default = True, help = 'Plot simple NJets histograms' )

parser.add_option( '--graphs', action = 'store_false',
                    dest = 'graphs',
                    default = True, help = 'Plot all the TGraphs of signifiance and number of events' )

parser.add_option( '--testCuts', action = 'store_false',
                    dest = 'testCuts',
                    default = True, help = 'Used calculated bin edges from one mass point on another signal mass points to see the significance' )

parser.add_option( '--testCutsTable', action = 'store_false',
                    dest = 'testCutsTable',
                    default = True, help = 'Used calculated bin edges from one mass point on another signal mass point and print out all LaTeX formatted event yields' )

parser.add_option( '--testCutsGraphs', action = 'store_false',
                    dest = 'testCutsGraphs',
                    default = True, help = 'Used calculated bin edges from one mass point on another signal mass point and print out all the TGraphs of significance' )

parser.add_option( '--year', action = 'store',
                    dest = 'year',
                    default = '2016', help = 'Run over the files from this year (2016, 2017, 2018)' )

parser.add_option( '--control', action = 'store_true',
                    dest = 'control',
                    default = False, help = 'Calculate the bin edges in the control region (non-isolated muon)' )

parser.add_option( '--noQcd', action = 'store_false',
                    dest = 'noQcd',
                    default = True, help = 'Calculate significance wihtout QCD' )

( options, args )       =  parser.parse_args()

def main() :
   
    # change the version tag if you are doing it over many different cfg files

    versionTag          = ""
    if options.year == "2016" :
        versionTag      = "Keras_V1.2.5"
    if options.year == "2017" :
        versionTag      = "Keras_V3.0.2"
    
    # treeDir is the location of your minitrees
    # outputDir is where you want your plots to be output

    treeDir             = "../../haddFiles/SignalMiniTrees/"
    outputDir           = "../../output/SignalBinStudy/PU/"
    outputDir           = outputDir+options.year+"/"+versionTag

    # define how fine a binning you want in your plots
    bkgdEffCuts         = [ i/1000.0 + 0.01 for i in range(700) ]

    # decide on which mass points you want to produce cuts that maximize signficance for (massPtCuts)
    # decide on how many mass points to run over and calculate the significance for (massPts)

    massPtCuts          = [ "550" ]
    massPts             = [ "350" ]

    # decide on number of bins
    nBins               = 4

    # decide on systematics size ( #FIXME: make this more realistic? )
    mainBkgdSys         = 0.2 # TT in our signal region
    subBkgdSys          = 0.2 # QCD in our signal region
   
    # njet binning
    njetbins        = [ 7, 8, 9, 10, 11, 12, 13, 14 ]

    # decide on what graphs you want to plot
    simpSignGraphs          = []        # plots of simple significance
    signGraphs              = []        # plots of significance
    nSigGraphs              = []        # plots of number of sig events
    njetHistos              = []        # plots of nJets
    sigNJetHistos           = []        # plots of signal Njet Histos
    fullNjetHistos          = []        # plots of nJets
    fullDESMHistos          = []        # histogams of DESM distribution
    
    finalEffs               = []        # final background efficiencies
    finalCuts               = {}
    njetHistosCuts          = []        # plots of nJets made in the cuts

    for njetbin in njetbins :
        finalCuts[njetbin] = []
    
    # get arrays for dataset where we require constant background fraction ( TT for signal, Data_SingleMuon for control region)
    mainBkgdDesmDict        = {}
    mainBkgdWghtDict        = {}
    mainBkgdNjetArray       = []

    mainBkgdName            = "TT"
    if options.control:
        mainBkgdName        = "Data_SingleMuon"
    mainBkgdDesmDict, mainBkgdWghtDict, mainBkgdNjetArray = makeNJetArrays( treeDir, options.year+"_"+mainBkgdName, njetbins )

    # get arrays for dataset that we are interested in ( QCD for signal, QCD for control region )
    #   for signal region, the QCD arrays are necessary to calculate the significance
    #   for control region, the QCD arrays are useful to look at the deepESM distribution in both MC and data and compare ( also run with --noQcd tag if deriving new bin edges )

    subBkgdDesmDict         = {}
    subBkgdWghtDict         = {}
    subBkgdNjetArray        = []
    
    subBkgdName             = "QCD"
    
    subBkgdDesmDict, subBkgdWghtDict, subBkgdNjetArray = makeNJetArrays( treeDir, options.year+"_"+subBkgdName, njetbins )

    # make diagnostic plots of the desm distribution and njet distribution for the bkgds

    if( options.plotDesm ) :
        plotSimpleDESMHists( mainBkgdDesmDict, mainBkgdWghtDict, subBkgdDesmDict, subBkgdWghtDict, njetbins, mainBkgdName, subBkgdName, outputDir )
        
    # for the mass point cut, lets find the correct cut for each efficiency
    for mass in massPtCuts :
        signalDesmDict      = {}
        signalWghtDict      = {}
        signalNjetArray     = []

        signalName          = "RPV_2t6j_mStop-"+mass

        signalDesmDict, signalWghtDict, signalNjetArray = makeNJetArrays( treeDir, options.year+"_"+signalName, njetbins ) 
        lastDesmCut         = {}

        # determine cuts, initialize by choosing the whole data set ( anything less than desm value of 1.0 )
        for njetbin in njetbins :
            lastDesmCut[njetbin] = 1.0
        
        # loop based on how many bins you want ( this really is number of cuts, so it is one less than number of bins )
        for nBinIter in range( nBins - 1 ) :

            # Create arrays for the bins excluding all those higher than the bin cut
            tempSignalDict              = {}
            tempSignalWghtDict          = {}
            tempSubBkgdDict             = {}
            tempSubBkgdWghtDict         = {}
            tempMainBkgdDict            = {}
            tempMainBkgdWghtDict        = {}

            tempSignalDict, tempSignalWghtDict, tempSubBkgdDict, tempSubBkgdWghtDict, tempMainBkgdDict, tempMainBkgdWghtDict = createTempArrays( lastDesmCut, signalDesmDict, signalWghtDict, subBkgdDesmDict, subBkgdWghtDict, mainBkgdDesmDict, mainBkgdWghtDict, njetbins )

            # calculate the total number of events in each bin (changes every cut you make)
            totBkgdEvts                 = {}
            tempCut                     = {}

            for njetbin in njetbins :
                totBkgdEvts[njetbin]    = sum( tempMainBkgdWghtDict[njetbin] )
                tempCut[njetbin]        = lastDesmCut[njetbin]
            
                print str(njetbin),":",totBkgdEvts[njetbin]
            
            # if there are no events in some of the bins, then ignore them in the significance calculation
            # this can also be useful if you have events in low stats njet bins that really impact significance

            ignore13Bin         = ( totBkgdEvts[13] == 0 )
            ignore14Bin         = ( totBkgdEvts[14] == 0 )

            # these six dictionaries will have per njet arrays so that TGraphs can be made of the variables with respects to bkgd eff / bkgd frac
            bkgdEffDict             = {}
            simpSignDict            = {}
            significanceDict        = {}
            nSignalEvtsDict         = {}
            nMainBkgdEvtsDict       = {}
            nSubBkgdEvtsDict        = {}

            # these six variables will temporarily hold the output of the calculate significance function
            desmCuts                = {}
            tempBkgdEff             = {}

            for njetbin in njetbins :
                
                bkgdEffDict[njetbin]        = array.array('d') # in anticipation for making some TGraphs
                simpSignDict[njetbin]       = array.array('d')
                significanceDict[njetbin]   = array.array('d') 
                nSignalEvtsDict[njetbin]    = array.array('d')
                nMainBkgdEvtsDict[njetbin]  = array.array('d')
                nSubBkgdEvtsDict[njetbin]   = array.array('d')

                tempBkgdEff[njetbin] = 0.0
                desmCuts[njetbin] = []

                if njetbin == 13 and ignore13Bin : 
                    desmCuts[13] = -1.0
                    
                    simpSignDict[njetbin].append( 0.0 )
                    significanceDict[njetbin].append( 0.0 )
                    nSignalEvtsDict[njetbin].append( 0.0 )
                    nMainBkgdEvtsDict[njetbin].append( 0.0 )
                    nSubBkgdEvtsDict[njetbin].append( 0.0 )
                    
                    continue
                
                if njetbin == 14 and ignore14Bin :
                    desmCuts[14] = -1.0
                    
                    simpSignDict[njetbin].append( 0.0 )
                    significanceDict[njetbin].append( 0.0 )
                    nSignalEvtsDict[njetbin].append( 0.0 )
                    nMainBkgdEvtsDict[njetbin].append( 0.0 )
                    nSubBkgdEvtsDict[njetbin].append( 0.0 )
                    
                    continue

                for bkgdEffCut in bkgdEffCuts:

                    bkgdEffDict[njetbin].append( float(bkgdEffCut) )
                    
                    while( tempBkgdEff[njetbin] < bkgdEffCut ) :
                        tempCut[njetbin] = tempCut[njetbin] - .001
                        tempBkgdEff[njetbin] = calcWghtEvts( tempCut[njetbin], tempMainBkgdWghtDict[njetbin], tempMainBkgdDict[njetbin] ) / totBkgdEvts[njetbin]
                    desmCuts[njetbin].append( tempCut[njetbin] )

                    tempSimpSignificance, tempSignificance, tempNSignalEvts, tempNMainBkgdEvts, tempNSubBkgdEvts = calcPerBinSignificance( desmCuts[njetbin], tempSignalWghtDict[njetbin], tempSignalDict[njetbin], tempMainBkgdWghtDict[njetbin], tempMainBkgdDict[njetbin], tempSubBkgdWghtDict[njetbin], tempSubBkgdDict[njetbin], mainBkgdSys, subBkgdSys )
                    
                    simpSignDict[njetbin].append( float(tempSimpSignificance) )
                    significanceDict[njetbin].append( float(tempSignificance) )
                    nSignalEvtsDict[njetbin].append( float(tempNSignalEvts) )
                    nMainBkgdEvtsDict[njetbin].append( float(tempNMainBkgdEvts) )
                    nSubBkgdEvtsDict[njetbin].append( float(tempNSubBkgdEvts) )

            simpSignDict["All"]         = calcTotSignificance( simpSignDict, njetbins )
            significanceDict["All"]     = calcTotSignificance( significanceDict, njetbins )
            nSignalEvtsDict["All"]      = sumOverAllJets( nSignalEvtsDict, njetbins )
            nMainBkgdEvtsDict["All"]    = sumOverAllJets( nMainBkgdEvtsDict, njetbins )
            nSubBkgdEvtsDict["All"]     = sumOverAllJets( nSubBkgdEvtsDict, njetbins )

            # here is the meat of the analyzer:
            #    matchBkgdEffwithMaxSig finds the maximum significance and gives the desm cut value per njetbin and the overall bgkd eff as a function of percentage of remaining events
            #    usePreDeterminedCuts finds the desm cut values that correspond to an overall bkgd efficiency

            signArrayIndex, bkgdEffPt, maxSignificance = matchBkgdEffWithMaxSig( significanceDict["All"], bkgdEffCuts )
            #signArrayIndex, bkgdEffPt, maxSignificance = usePreDeterminedCuts( significanceDict["All"], bkgdEffCuts, nBinIter )
            # update the lastest cut for recursion
            
            finalEffs.append( bkgdEffPt )
            
            for njetbin in njetbins:
                lastDesmCut[njetbin]  = desmCuts[njetbin][signArrayIndex]
                finalCuts[njetbin].append( lastDesmCut[njetbin] )

           
            # fill a histogram with the optimal cut and njet distribution and print for each MVA bin
            if( options.njets ) :
                histoName           = 'h_njets_RPV-m'+mass+'_bin'+str( nBins - nBinIter )+'_bkgdCut'+str( bkgdEffPt )
                if( options.debug ) :
                    print 'Making histogram: '+histoName 
                
                # if this is the last iteration, have the bulk histogram (bin 1) also
                if nBinIter == nBins - 2 :
                    histoNameBin1   = 'h_njets_RPV-m'+mass+'_bin'+str( nBins - nBinIter - 1 )
                    if( options.debug ) :
                        print 'Making histogram: '+histoNameBin1 
                    nJetsHistoBin1  = ROOT.TH1F( histoNameBin1, histoNameBin1, 8, 7, 15 )

                nJetsHisto          = ROOT.TH1F( histoName, histoName, 8, 7, 15 )
    
                # fill in histograms
                for njetbin in njetbins :
                    for iDesm in range( len( tempMainBkgdDict[njetbin] ) ) :
                        if tempMainBkgdDict[njetbin][iDesm] > lastDesmCut[njetbin] :
                            nJetsHisto.Fill( njetbin, tempMainBkgdWghtDict[njetbin][iDesm] )
                        else :
                            if nBinIter == nBins - 2 :
                                nJetsHistoBin1.Fill( njetbin, tempMainBkgdWghtDict[njetbin][iDesm] )
                    if not options.noQcd:
                        for iDesm in range( len( tempSubBkgdDict[njetbin] ) ) :
                            if tempSubBkgdDict[njetbin][iDesm] > lastDesmCut[njetbin] :
                                nJetsHisto.Fill( njetbin, tempSubBkgdWghtDict[njetbin][iDesm] )
                            else :
                                if nBinIter == nBins - 2 :
                                    nJetsHistoBin1.Fill( njetbin, tempSubBkgdWghtDict[njetbin][iDesm] )
    
                if options.debug :
                    print "Total number of weighted events in bin "+str( nBins - nBinIter )+": "+nJetsHisto.Integral()
                plotHistSimp( nJetsHisto, outputDir+'/'+histoName+'.pdf' )
                nJetsHisto.Scale( 1.0/nJetsHisto.Integral() )
                plotHistSimp( nJetsHisto, outputDir+'/'+histoName+'_norm.pdf' )
                njetHistosCuts.append( copy.deepcopy( nJetsHisto ) )
                
                if nBinIter == nBins-2 :
                    if options.debug :
                         print "Total number of weighted events in bin "+str( nBins - nBinIter - 1 )+": "+nJetsHisto.Integral()
                    plotHistSimp( nJetsHistoBin1,outputDir+'/'+histoNameBin1+'.pdf' )
                    nJetsHistoBin1.Scale( 1.0/nJetsHistoBin1.Integral() )
                    plotHistSimp( nJetsHistoBin1, outputDir+'/'+histoName+'_norm.pdf' )
                    njetHistosCuts.append( copy.deepcopy( nJetsHistoBin1 ) )

            # make plots of different variables vs bkgd efficiency/bgkd fraction 
            if( options.graphs ) :
                for njetbin in njetbins :
                    tempSignificanceGraph    = ROOT.TGraph( len(bkgdEffCuts), bkgdEffDict[njetbin], significanceDict[njetbin] )
                    tempSignalEvtsGraph     = ROOT.TGraph( len(bkgdEffCuts), bkgdEffDict[njetbin], nSignalEvtsDict[njetbin] )
                    tempMainBkgdEvtsGraph     = ROOT.TGraph( len(bkgdEffCuts), bkgdEffDict[njetbin], nMainBkgdEvtsDict[njetbin] )
                    tempSubBkgdEvtsGraph     = ROOT.TGraph( len(bkgdEffCuts), bkgdEffDict[njetbin], nSubBkgdEvtsDict[njetbin] )
                    
                    plotGraphSimp( tempSignificanceGraph, mass, str( nBins - nBinIter ), str( njetbin ), outputDir, "Significance" )
                    plotGraphSimp( tempSignalEvtsGraph, mass, str( nBins - nBinIter ), str( njetbin ), outputDir, signalName )
                    plotGraphSimp( tempMainBkgdEvtsGraph, mass, str( nBins - nBinIter ), str( njetbin ), outputDir, mainBkgdName )
                    plotGraphSimp( tempSubBkgdEvtsGraph, mass, str( nBins - nBinIter ), str( njetbin ), outputDir, subBkgdName )
    
    # once all the njets distributions per bin is made, overlay all the MVA njet distributions together
    if( options.njets ) :
        plotJetBinShapes( nBins, massPtCuts, njetHistosCuts, outputDir ) 

    # print out the output that is needed for the cfg
    printDeepEventShapeCfg( finalEffs, finalCuts, njetbins )

    if options.testCuts :
        for mass in massPts :
            signalDesmDict      = {}
            signalWghtDict      = {}
            signalNjetArray     = []
    
            signalName          = "RPV_2t6j_mStop-"+mass
    
            signalDesmDict, signalWghtDict, signalNjetArray = makeNJetArrays( treeDir, options.year+"_"+signalName, njetbins ) 
            lastDesmCut         = {}
        
            # determine cuts, initialize by choosing the whole data set ( anything less than desm value of 1.0 )
            for njetbin in njetbins :
                lastDesmCut[njetbin] = 1.0
            
            # loop based on how many bins you want ( this really is number of cuts, so it is one less than number of bins )
            for nBinIter in range( nBins - 1 ) :
    
                # Create arrays for the bins excluding all those higher than the bin cut
                tempSignalDict              = {}
                tempSignalWghtDict          = {}
                tempSubBkgdDict             = {}
                tempSubBkgdWghtDict         = {}
                tempMainBkgdDict            = {}
                tempMainBkgdWghtDict        = {}
    
                tempSignalDict, tempSignalWghtDict, tempSubBkgdDict, tempSubBkgdWghtDict, tempMainBkgdDict, tempMainBkgdWghtDict = createTempArrays( lastDesmCut, signalDesmDict, signalWghtDict, subBkgdDesmDict, subBkgdWghtDict, mainBkgdDesmDict, mainBkgdWghtDict, njetbins )
    
                # calculate the total number of events in each bin (changes every cut you make)
                totBkgdEvts                 = {}
                tempCut                     = {}
    
                for njetbin in njetbins :
                    totBkgdEvts[njetbin]    = sum( tempMainBkgdWghtDict[njetbin] )
                    tempCut[njetbin]        = lastDesmCut[njetbin]
                
                    print str(njetbin),":",totBkgdEvts[njetbin]
                
                # if there are no events in some of the bins, then ignore them in the significance calculation
                # this can also be useful if you have events in low stats njet bins that really impact significance
    
                ignore13Bin         = ( totBkgdEvts[13] == 0 )
                ignore14Bin         = ( totBkgdEvts[14] == 0 )
    
                # these six dictionaries will have per njet arrays so that TGraphs can be made of the variables with respects to bkgd eff / bkgd frac
                bkgdEffDict             = {}
                simpSignDict            = {}
                significanceDict        = {}
                nSignalEvtsDict         = {}
                nMainBkgdEvtsDict       = {}
                nSubBkgdEvtsDict        = {}
    
                # these six variables will temporarily hold the output of the calculate significance function
                desmCuts                = {}
                tempBkgdEff             = {}
    
                for njetbin in njetbins :
                    
                    bkgdEffDict[njetbin]        = array.array('d') # in anticipation for making some TGraphs
                    simpSignDict[njetbin]       = array.array('d')
                    significanceDict[njetbin]   = array.array('d') 
                    nSignalEvtsDict[njetbin]    = array.array('d')
                    nMainBkgdEvtsDict[njetbin]  = array.array('d')
                    nSubBkgdEvtsDict[njetbin]   = array.array('d')
    
                    tempBkgdEff[njetbin] = 0.0
                    desmCuts[njetbin] = []
    
                    if njetbin == 13 and ignore13Bin : 
                        desmCuts[13] = -1.0
                        
                        simpSignDict[njetbin].append( 0.0 )
                        significanceDict[njetbin].append( 0.0 )
                        nSignalEvtsDict[njetbin].append( 0.0 )
                        nMainBkgdEvtsDict[njetbin].append( 0.0 )
                        nSubBkgdEvtsDict[njetbin].append( 0.0 )
                        
                        continue
                    
                    if njetbin == 14 and ignore14Bin :
                        desmCuts[14] = -1.0
                        
                        simpSignDict[njetbin].append( 0.0 )
                        significanceDict[njetbin].append( 0.0 )
                        nSignalEvtsDict[njetbin].append( 0.0 )
                        nMainBkgdEvtsDict[njetbin].append( 0.0 )
                        nSubBkgdEvtsDict[njetbin].append( 0.0 )
                        
                        continue
    
                    for bkgdEffCut in bkgdEffCuts:
    
                        bkgdEffDict[njetbin].append( float(bkgdEffCut) )
                        
                        while( tempBkgdEff[njetbin] < bkgdEffCut ) :
                            tempCut[njetbin] = tempCut[njetbin] - .001
                            tempBkgdEff[njetbin] = calcWghtEvts( tempCut[njetbin], tempMainBkgdWghtDict[njetbin], tempMainBkgdDict[njetbin] ) / totBkgdEvts[njetbin]
                        desmCuts[njetbin].append( tempCut[njetbin] )
    
                        tempSimpSignificance, tempSignificance, tempNSignalEvts, tempNMainBkgdEvts, tempNSubBkgdEvts = calcPerBinSignificance( desmCuts[njetbin], tempSignalWghtDict[njetbin], tempSignalDict[njetbin], tempMainBkgdWghtDict[njetbin], tempMainBkgdDict[njetbin], tempSubBkgdWghtDict[njetbin], tempSubBkgdDict[njetbin], mainBkgdSys, subBkgdSys )
                        
                        simpSignDict[njetbin].append( float(tempSimpSignificance) )
                        significanceDict[njetbin].append( float(tempSignificance) )
                        nSignalEvtsDict[njetbin].append( float(tempNSignalEvts) )
                        nMainBkgdEvtsDict[njetbin].append( float(tempNMainBkgdEvts) )
                        nSubBkgdEvtsDict[njetbin].append( float(tempNSubBkgdEvts) )
                        
                        if( bkgdEffCut == finalEffs[nBinIter] ) :
                            if options.testCutsTable:
                                if njetbin == 7:
                                    print "Event Yield and Significance per NJets for "+signalName+" in Bin "+str(nBins-nBinIter)
                                    print( "%i,%.2f,%.2f,%.2f,%.4f,%.4f" % ( njetbin, nSignalEvtsDict[njetbin][-1], nMainBkgdEvtsDict[njetbin][-1], nSubBkgdEvtsDict[njetbin][-1], simpSignDict[njetbin][-1], significanceDict[njetbin] [-1]) ) 
                                else :
                                    print( "%i,%.2f,%.2f,%.2f,%.4f,%.4f" % ( njetbin, nSignalEvtsDict[njetbin][-1], nMainBkgdEvtsDict[njetbin][-1], nSubBkgdEvtsDict[njetbin][-1], simpSignDict[njetbin][-1], significanceDict[njetbin] [-1]) ) 
                        
                for njetbin in njetbins:
                    lastDesmCut[njetbin]  = finalCuts[njetbin][nBinIter]
            
            
                if( options.testCutsGraphs ) :
                    for njetbin in njetbins :
                        tempSignifianceGraph    = ROOT.TGraph( len(bkgdEffCuts), bkgdEffDict[njetbin], significanceDict[njetbin] )
                        tempSignalEvtsGraph     = ROOT.TGraph( len(bkgdEffCuts), bkgdEffDict[njetbin], nSignalEvtsDict[njetbin] )
                        tempMainBkgdEvtsGraph     = ROOT.TGraph( len(bkgdEffCuts), bkgdEffDict[njetbin], nMainBkgdEvtsDict[njetbin] )
                        tempSubBkgdEvtsGraph     = ROOT.TGraph( len(bkgdEffCuts), bkgdEffDict[njetbin], nSubBkgdEvtsDict[njetbin] )
                        
                        plotGraphSimp( tempSignificanceGraph, mass, str( nBins - nBinIter ), str( njetbin ), outputDir, "Significance" )
                        plotGraphSimp( tempSignalEvtsGraph, mass, str( nBins - nBinIter ), str( njetbin ), outputDir, signalName )
                        plotGraphSimp( tempMainBkgdEvtsGraph, mass, str( nBins - nBinIter ), str( njetbin ), outputDir, mainBkgdName )
                        plotGraphSimp( tempSubBkgdEvtsGraph, mass, str( nBins - nBinIter ), str( njetbin ), outputDir, subBkgdName )

    # plot significances between bins - FIXME: need to update
    #plotBinShapes( nBins, massPts, len( bkgdEffCuts ), signGraphs, 'significance', outputDir )

def plotJetBinShapes( nBins, massPts, myHistoList, outputDir ) :

    # create a file in case these njet distributions have to be manipulated and used later
    myNewOutFile =  ROOT.TFile( outputDir+"/"+options.year+"_binStudy.root", "RECREATE" )

    for binIter in range( nBins ) :
        cs2     = ROOT.TCanvas( "cs2", "cs2", 700, 700 )
        pad2    = ROOT.TPad( "pad2", "pad2", 0, 0, 1.0, 1.0 )
        pad2.SetLeftMargin( 0.15 )
        pad2.SetLogy()
        pad2.Draw()
        leg2    = ROOT.TLegend( 0.8, 0.7, 0.95, 0.90 )
        pad2.cd()
        
        for massIter in range( len( massPts ) ) :
            histo = myHistoList[ (nBins)*massIter + binIter ]
            histo.SetLineColor( decideColor( massPts[massIter] ) )
            histo.Scale( 1.0/histo.Integral() )
            histo.SetStats( 0 )
            histo.SetMaximum( 1.1 )
            histo.SetMinimum( 0.00001 )
            leg2.AddEntry( histo, decideLegendTag( massPts[massIter] ), 'l' )
            
            if massIter == 0 :
                print massIter, binIter, (nBins)*massIter
                histo.GetXaxis().SetTitle( 'Number of Jets' )
                histo.GetYaxis().SetTitle( 'Normalized Number of Events' )
                histo.SetTitle( 'NJets Distribution for Bin '+str( nBins-binIter ) )
                histo.Draw()
                histo.Write()
            else :
                print massIter, binIter, (nBins)*massIter
                histo.Draw( 'SAME' )
                histo.Write()
        leg2.Draw()
        cs2.SaveAs(outputDir+'/njets_bin'+str(nBins-binIter)+'.pdf')
        if 'cs2' in locals() :
            del cs2, pad2, leg2

    myNewOutFile.Close()

    for massIter in range( len( massPts ) ) :

        cs2     = ROOT.TCanvas( "cs2", "cs2", 700, 700 )
        pad2    = ROOT.TPad( "pad2", "pad2", 0, 0, 1.0, 1.0 )
        pad2.SetLeftMargin( 0.15 )
        pad2.SetLogy()
        pad2.Draw()
        leg2    = ROOT.TLegend( 0.8, 0.75, 0.95, 0.90 )
        pad2.cd()
        
        for binIter in range( nBins ) :
            histo = myHistoList[ (nBins)*massIter + binIter ]
            histo.SetLineColor( decideColorBin( nBins-binIter ) )
            histo.Scale( 1.0/histo.Integral() )
            histo.SetStats( 0 )
            histo.SetMaximum( 1.1 )
            histo.SetMinimum( 0.00001 )
            leg2.AddEntry( histo, decideLegendTagBinHisto( nBins-binIter ), 'l ')

            if binIter == 0 :
                print massIter, binIter, (nBins)*massIter
                histo.GetXaxis().SetTitle( 'Number of Jets' )
                histo.GetYaxis().SetTitle( 'Normalized Number of Events' )
                #histo.SetTitle( 'NJets Distribution with Mass '+massPts[massIter] )
                histo.Draw()
            else :
                print massIter, binIter, (nBins)*massIter
                histo.Draw( 'SAME' )

        leg2.Draw()
        cs2.SaveAs( outputDir+'/njets_m'+massPts[massIter]+'.pdf' )
        if 'cs2' in locals() :
            del cs2, pad2, leg2

def decideColorBin( filename ) :
    if filename == 1 :
        return ROOT.kRed
    elif filename == 2 :
        return ROOT.kOrange
    elif filename == 3 :
        return ROOT.kGreen
    elif filename == 4 :
        return ROOT.kCyan
    elif filename == 5 :
        return ROOT.kBlue

def decideColor( filename ) :
    if "350" in filename :
        return ROOT.kRed
    elif "450" in filename :
        return ROOT.kOrange
    elif "550" in filename :
        return ROOT.kGreen
    elif "650" in filename :
        return ROOT.kCyan
    elif "750" in filename :
        return ROOT.kBlue
    elif "850" in filename :
        return ROOT.kViolet
    elif "all" in filename :
        return ROOT.kBlack

def decideLegendTagBinHisto( filename ) :
    return "Bin "+str( filename )

def decideLegendTagBinEff( filename ) :
    return "Bin "+str( filename )+" -> Bin "+str( filename - 1 )

def decideLegendTag( filename ):
    if "all" in filename :
        return "All RPV"
    else :
        return "RPV m"+str( filename )+" GeV"

def plotBinShapes( nBins, massPts, graphLength, myGraphList, nameAppend, outputDir) :
    for binIter in range( nBins - 1 ) :

        cs1     = ROOT.TCanvas( "cs1", "cs1", 700, 700 )
        pad1    = ROOT.TPad( "pad1", "pad1", 0, 0, 1.0, 1.0 )
        pad1.SetLeftMargin( 0.15 )
        pad1.Draw()
        leg1    = ROOT.TLegend( 0.8, 0.65, 0.95, 0.90 )
        pad1.cd()

        graphMax = 0.0
        graphMin = 1.0

        for massIter in range( len( massPts ) ) :
            graph = myGraphList[ (nBins-1)*massIter + binIter ]
            if( ROOT.TMath.MaxElement( graphLength, graph.GetY() ) > graphMax ) :
                graphMax = ROOT.TMath.MaxElement( graphLength, graph.GetY() )
            if( ROOT.TMath.MinElement( graphLength, graph.GetY() ) < graphMin ) :
                graphMin = ROOT.TMath.MinElement( graphLength, graph.GetY() )
            graph.SetLineColor( decideColor( massPts[massIter] ) )
            leg1.AddEntry( graph, decideLegendTag( massPts[massIter] ), 'l' )

        for massIter in range( len( massPts ) ) :
            graph = myGraphList[ (nBins-1)*massIter + binIter ]
            if massIter == 0 :
                graph.SetMaximum( 1.1*graphMax )
                graph.SetMinimum( 0.9*graphMin )
                graph.Draw()
            else :
                graph.Draw( 'SAME' )

        leg1.Draw()
        cs1.SaveAs(outputDir+'/'+nameAppend+'_bin'+str(nBins-binIter)+'.pdf')
        if 'cs1' in locals() :
            del cs1, pad1, leg1
        
    for massIter in range( len( massPts ) ) :
        cs1     = ROOT.TCanvas( "cs1", "cs1", 700, 700 )
        pad1    = ROOT.TPad( "pad1", "pad1", 0, 0, 1.0, 1.0 )
        pad1.SetLeftMargin( 0.15 )
        pad1.Draw()
        leg1    = ROOT.TLegend( 0.8, 0.7, 0.98, 0.90 )
        pad1.cd()
        
        graphMax = 0.0
        graphMin = 1.0

        for binIter in range( nBins - 1 ) :
            graph = myGraphList[ (nBins-1)*massIter + binIter ]
            graph.SetLineColor( decideColorBin( nBins-binIter ) )
            if( ROOT.TMath.MaxElement( graphLength, graph.GetY() ) > graphMax ) :
                graphMax = ROOT.TMath.MaxElement( graphLength, graph.GetY() )
            if( ROOT.TMath.MinElement( graphLength, graph.GetY() ) < graphMin ) :
                graphMin = ROOT.TMath.MinElement( graphLength, graph.GetY() )
            leg1.AddEntry( graph, decideLegendTagBinEff( nBins-binIter ), 'l ')

        for binIter in range( nBins - 1 ) :
            graph = myGraphList[ (nBins-1)*massIter + binIter ]
            if binIter == 0 :
                graph.SetMaximum( 1.1*graphMax )
                graph.SetMinimum( 0.9*graphMin )
                graph.Draw()
            else :
                graph.Draw( 'SAME' )
        leg1.Draw()
        cs1.SaveAs( outputDir+'/'+nameAppend+'_'+massPts[massIter]+'.pdf' )
        if 'cs1' in locals() :
            del cs1, pad1, leg1


def usePreDeterminedCuts( signArray, bkgdEffCuts, binNumber ) :
    #myArray = [ .024, .065, .389 ] # 2016 Optimized Cuts - 80.0 fb-1
    myArray = [ .024, .052, .395 ] # 2016 Optimized Cuts - 35.5 fb-1
    myCut = myArray[binNumber]
    for myIter in range( len( bkgdEffCuts ) ) :
        if abs(bkgdEffCuts[myIter] - myCut) < 0.0001 :
            print bkgdEffCuts[myIter], signArray[myIter]
            return myIter, bkgdEffCuts[myIter], signArray[myIter]
        
def matchBkgdEffWithMaxSig( signArray, bkgdEffCuts ) :
    for binIter in range( len( bkgdEffCuts ) ) :
        if signArray[binIter] - max( signArray ) > -0.00001 :
            print bkgdEffCuts[binIter], max(signArray)
            return binIter, bkgdEffCuts[binIter], max(signArray)


# this function calculates the total weighted events in a list 
def calcWghtEvts( cut, wght, List ) :
    nEvts   = 0.0
    for iEvt in range( len( List ) ) :
        if List[iEvt] > cut :
            nEvts += wght[iEvt]
    return nEvts

# this function takes the important information from the miniTrees and puts them into dictionaries/arrays
def makeNJetArrays( treeDir, fileName, njetbins ) :
    
    # get arrays from selected miniTree 
    print "Making njet arrays for this file: "+treeDir+"/"+fileName+".root"
    myFile                  = ROOT.TFile( treeDir+"/"+fileName+".root" )
    myTree                  = myFile.Get( "myMiniTree" )
    
    desmDict                = {}
    wghtDict                = {}
    njetArray              = []
    
    for njetbin in njetbins :
        desmDict[njetbin]   = []
        wghtDict[njetbin]   = []

    desmDict["All"]         = []
    wghtDict["All"]         = []

    for entry in xrange( 0, myTree.GetEntries() ) :
        myTree.GetEntry( entry )
       
        # determine the weights. If the main bkgd is the data from the control region, then the weight is 1.0. Otherwise it is the weight with all SFs applied ( bTagSF,leptonSF,htSF,prefiringSF(2017),puSF )
        if "Data" in fileName:
            wghtDict["All"].append( 1.0 )
        else:
            wghtDict["All"].append( myTree.totalEventWeight*myTree.Lumi )

        desmDict["All"].append( myTree.deepESM_val )
        njetArray.append( myTree.NGoodJets_pt30 )

        for njetbin in njetbins :
            if myTree.NGoodJets_pt30 == njetbin :
                desmDict[njetbin].append( myTree.deepESM_val )
                if "Data" in fileName :
                    wghtDict[njetbin].append( 1.0 )
                else :
                    wghtDict[njetbin].append( myTree.totalEventWeight*myTree.Lumi )
                break
    
    # do some simple printm statements
    print 'Number of weighted '+fileName+' events: ',sum(wghtDict["All"])
    print 'Number of unweighted '+fileName+' events: ',len(wghtDict["All"])

    if( options.debug ) :
        print 'Closing '+fileName+' File'

    #Close the file
    myFile.Close()

    return desmDict, wghtDict, njetArray
        
# the following function just plots the deep esm distribution for the main bkgd and the sub bkgd - useful for debugging

def plotSimpleDESMHists( mainBkgdDesmDict, mainBkgdWghtDict, subBkgdDesmDict, subBkgdWghtDict, njetbins, mainBkgdName, subBkgdName, outputDir ) :
    
    for njetbin in njetbins :
        mainBkgdHistName    = "h_deepESM_"+mainBkgdName+"_"+str(njetbin)+"j"
        subBkgdHistName     = "h_deepESM_"+subBkgdName+"_"+str(njetbin)+"j"

        mainBkgdHist        = ROOT.TH1F( mainBkgdHistName, mainBkgdHistName, 15, 0.0, 1.0 )
        subBkgdHist         = ROOT.TH1F( subBkgdHistName, subBkgdHistName, 15, 0.0, 1.0 )

        for iDesm in range( len( mainBkgdDesmDict[njetbin] ) ):
            mainBkgdHist.Fill( mainBkgdDesmDict[njetbin][iDesm], mainBkgdWghtDict[njetbin][iDesm] )
        for iDesm in range( len( subBkgdDesmDict[njetbin] ) ):
            subBkgdHist.Fill( subBkgdDesmDict[njetbin][iDesm], subBkgdWghtDict[njetbin][iDesm] )

        plotHistSimp( mainBkgdHist, outputDir+'/'+mainBkgdHistName+'.pdf' )
        plotHistSimp( subBkgdHist, outputDir+'/'+subBkgdHistName+'.pdf' )

# this function creates the arrays that remove the deep esm values / events that are already in a bin 
def createTempArrays( lastDesmCut, signalDesmDict, signalWghtDict, subBkgdDesmDict, subBkgdWghtDict, mainBkgdDesmDict, mainBkgdWghtDict, njetbins ) :

    tempSignalDict          = {}
    tempSignalWghtDict      = {}
    tempSubBkgdDict         = {}
    tempSubBkgdWghtDict     = {}
    tempMainBkgdDict        = {}
    tempMainBkgdWghtDict    = {}

    for njetbin in njetbins :
        
        tempSignalDict[njetbin]         = []
        tempSignalWghtDict[njetbin]     = []
        tempSubBkgdDict[njetbin]        = []
        tempSubBkgdWghtDict[njetbin]    = []
        tempMainBkgdDict[njetbin]       = []
        tempMainBkgdWghtDict[njetbin]   = []

        for iEvt in range( len( signalDesmDict[njetbin] ) ) :
            if signalDesmDict[njetbin][iEvt] <  lastDesmCut[njetbin] :
                tempSignalDict[njetbin].append( signalDesmDict[njetbin][iEvt] )
                tempSignalWghtDict[njetbin].append( signalWghtDict[njetbin][iEvt] )
    
        for iEvt in range( len( subBkgdDesmDict[njetbin] ) ) :
            if subBkgdDesmDict[njetbin][iEvt] <  lastDesmCut[njetbin] :
                tempSubBkgdDict[njetbin].append( subBkgdDesmDict[njetbin][iEvt] )
                tempSubBkgdWghtDict[njetbin].append( subBkgdWghtDict[njetbin][iEvt] )
        
        for iEvt in range( len( mainBkgdDesmDict[njetbin] ) ) :
            if mainBkgdDesmDict[njetbin][iEvt] <  lastDesmCut[njetbin] :
                tempMainBkgdDict[njetbin].append( mainBkgdDesmDict[njetbin][iEvt] )
                tempMainBkgdWghtDict[njetbin].append( mainBkgdWghtDict[njetbin][iEvt] )

    return tempSignalDict, tempSignalWghtDict, tempSubBkgdDict, tempSubBkgdWghtDict, tempMainBkgdDict, tempMainBkgdWghtDict

# this function calculates the significance and number of weighted signal, main bkgd, and sub bkgd events (may be used to make plots )
def calcPerBinSignificance( cut, sigWghtDict, sigDict, mainBkgdWghtDict, mainBkgdDict, subBkgdWghtDict, subBkgdDict, mainBkgdSys, subBkgdSys ) :
        
    nSignalEvts    = 0.0
    nMainBkgdEvts  = 0.0
    nSubBkgdEvts   = 0.0

    for iEvt in range( len( sigDict ) ) :
        if sigDict[iEvt] >  cut[-1] :
            nSignalEvts += sigWghtDict[iEvt]
    for iEvt in range( len( mainBkgdDict ) ) :
        if mainBkgdDict[iEvt] >  cut[-1] :
            nMainBkgdEvts += mainBkgdWghtDict[iEvt]
    for iEvt in range( len( subBkgdDict ) ) :
        if subBkgdDict[iEvt] >  cut[-1] :
            nSubBkgdEvts += subBkgdWghtDict[iEvt]

    if options.noQcd :
        perBinSoverRtB  = nSignalEvts / math.sqrt( nMainBkgdEvts )
        perBinSignificance        = nSignalEvts / math.sqrt( nMainBkgdEvts + mainBkgdSys*mainBkgdSys*nMainBkgdEvts*nMainBkgdEvts )
    else :
        perBinSoverRtB  = nSignalEvts / math.sqrt( nMainBkgdEvts + nSubBkgdEvts )
        perBinSignificance        = nSignalEvts / math.sqrt( nSubBkgdEvts + subBkgdSys*subBkgdSys*nSubBkgdEvts*nSubBkgdEvts + nSubBkgdEvts + subBkgdSys*subBkgdSys*nSubBkgdEvts*nSubBkgdEvts )

    return perBinSoverRtB, perBinSignificance, nSignalEvts, nMainBkgdEvts, nSubBkgdEvts

# calculate total significance
def calcTotSignificance( significanceDict, njetbins ) :
    
    myTotalSignificanceArray  =  array.array( 'd' )

    for i in xrange( 0, len( significanceDict[7] ) ) :
        tempSignificanceSquared     = 0.0
        for njetbin in njetbins :
            tempSignificanceSquared += significanceDict[njetbin][i]*significanceDict[njetbin][i]
        myTotalSignificanceArray.append( math.sqrt( tempSignificanceSquared ) )

    return myTotalSignificanceArray 

# calculate sum total number of events :
def sumOverAllJets( evtsDict, njetbins ) :
    myTotalSumArray           = array.array( 'd' )

    for i in xrange( 0, len( evtsDict[7] ) ) :
        tempSum               = 0.0
        for njetbin in njetbins :
            tempSum += evtsDict[njetbin][i]
        myTotalSumArray.append( tempSum )

    return myTotalSumArray

# as its name suggest, it simply plots a histogram
def plotHistSimp( histo, outName ) :
    cTemp           = ROOT.TCanvas( "cTemp", "cTemp", 0, 0, 700, 700 )
    pTemp           = ROOT.TPad( "pTemp", "pTemp", 0, 0, 1.0, 1.0 )
    pTemp.SetLeftMargin( 0.15 )
    pTemp.SetLogy()
    pTemp.Draw()
    pTemp.cd()
    histo.Draw()
    cTemp.Update()
    cTemp.SaveAs( outName )
    del cTemp, pTemp

# simple graph plotting script for debugging - can elaborate if needed 
def plotGraphSimp( graph, mass, binNum, njetBin, outputDir, yaxisName ) :
    cTemp           = ROOT.TCanvas( "cTemp", "cTemp", 0, 0, 700, 700 )
    pTemp           = ROOT.TPad( "pTemp", "pTemp", 0, 0, 1.0, 1.0 )
    
    if yaxisName == "Signifiance" :
        graph.SetTitle( 'Significance for Mass Pt '+mass+' Bin '+binNum )
        graph.GetYaxis().SetTitle( 'Significance' )
    else :
        graph.SetTitle( 'Number of '+yaxisName+' Events for Mass Pt '+mass+' Bin ' +binNum )
        graph.GetYaxis().SetTitle( 'Number of '+yaxisName+' Events' )
    graph.SetLineWidth( 2 )

    pTemp.SetLeftMargin( 0.15 )
    pTemp.Draw()
    pTemp.cd()
    graph.Draw()
    cTemp.Update()
    cTemp.SaveAs( outputDir+'/'+yaxisName+'_m'+mass+'_bin'+binNum+'_njet'+njetBin+'.pdf' )
    del cTemp, pTemp

# print out the bin edges for the config file
def printDeepEventShapeCfg( finalEffs, finalCutsDict, njetbins ) :
    
    count = 0
    for njetbin in njetbins:
        print '    binEdges['+str(count)+'] = 0.00000000'
        count += 1
        for itCut in range( len( finalCutsDict[njetbin] ) ):
            print '    binEdges['+str(count)+'] = '+str(finalCutsDict[njetbin][len(finalCutsDict[njetbin])-itCut-1])
            count += 1
        print '    binEdges['+str(count)+'] = 1.00000000'
        count += 1

if __name__ == '__main__' :
    main()
