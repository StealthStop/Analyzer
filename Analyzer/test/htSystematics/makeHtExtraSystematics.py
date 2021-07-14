#!/bin/python
import ROOT
import copy
import os.path
import array
import math

ROOT.gROOT.SetBatch( True )

from optparse import OptionParser

parser  = OptionParser()

parser.add_option('--log', action='store_false',
                    dest='log',
                    default = True, help = 'Make all plots log on the y scale' )

parser.add_option('--norm', action='store_false',
                    dest='norm',
                    default = True, help = 'Normalize all distributions to unit area before doing any divisions' )

parser.add_option('--makePlots', action='store_true',
                    dest='makePlotsBool',
                    default = False, help = 'Make pdf versions comparing the N_{jets} distributions for all the H_{T} reweightings ' )

(options,args) = parser.parse_args()

def main():

    yearList            = [ "2016", "2017", "2018pre", "2018post" ]
    yearList            = [ "2016" ]
    
    snnBinArray         = [ "" , "_D1", "_D2", "_D3", "_D4" ]
    wghtArray           = [ "noHTWght", "htWght", "htWghtFlat2000", "htWghtNJet7" ]
    bkgdArray           = [ "TTX", "Other", "QCD", "TT" ]
    histoNameArray      = []

    for year in yearList :
        
        outputDir       = "htExtraSystRootFiles"
        
        if not os.path.exists( outputDir) :
            print "Output directory does not exist - making directory",outputDir
            os.makedirs( outputDir )
        
        if not os.path.exists( outputDir+"/"+year) :
            os.makedirs( outputDir+"/"+year )
        
        histoDict           = {}
        
        for snnBin in snnBinArray:
            
            histoDict[snnBin]               = {}
            histoDict[snnBin]["All"]        = {}
            
            for bkgd in bkgdArray :
                inputFileName       = "htRootFiles/"+year+"_HTSystematic_"+bkgd+".root"
                inputFile           = ROOT.TFile.Open( inputFileName )
                histoDict[snnBin][bkgd] = {}
                for wght in wghtArray : 
                    tempHisto       = inputFile.Get( "h_njets_"+wght+"_Allj"+snnBin )
                    histoDict[snnBin][bkgd][wght] = copy.deepcopy( tempHisto )
            inputFile.Close()
    
            for bkgd in bkgdArray:
                for wght in wghtArray:
                    if bkgd == bkgdArray[0]:
                        histoDict[snnBin]["All"][wght] = histoDict[snnBin][bkgd][wght].Clone( "h_njets_AllMC_Allj"+snnBin )
                    else:
                        histoDict[snnBin]["All"][wght].Add( histoDict[snnBin][bkgd][wght] )
    
        # To plot the N_{jets} shapes for different H_{T} reweightings - uncomment out this line 
        # plotMCDiff( outputDir, year, histoDict, wghtArray, bkgdArray, snnBinArray )

        # To make the extra H_{T} systematics ROOT file
        #       - The bool on the right is for whether you want to also make the corresponding plots
        #       - The "htWght" flag is used as the denominator in the comparison - can change to no HT weight if you want systematics to compare to the totally unweighted H_{T} distributions.

        makeExtraHtSystematics( outputDir, year, histoDict, wghtArray, bkgdArray, histoNameArray, snnBinArray, "htWght" )

    
def makeExtraHtSystematics( outputDir, year, originalHistoDict, wghtArray, bkgdArray, histoNameArray, snnBinArray, standName ) :
    
    outputHistos            = {}
   
    for snnBin in snnBinArray :
        for wght in wghtArray :
            if wght == standName :
                continue
            c1, p1, p2       = makeCanvasAndTwoPads()
            p1.cd()
            
            l1                          = ROOT.TLegend( 0.41, 0.70, 0.85, 0.85 )
            l1.SetBorderSize(0)
            l1.SetFillStyle( 0 )
            l1.SetLineWidth( 3 )
            l1.SetTextFont( 42 )
            l1.SetTextSize( 0.03 )
   
            histoDict = copy.deepcopy( originalHistoDict )
            
            h_njets = histoDict[snnBin]["TT"][ wght ]
            h_njets.SetStats(0)
            h_njets.SetLineColor(ROOT.kRed)
            
            h_stand = histoDict[snnBin]["TT"][ standName ]
            h_stand.SetStats(0)
            h_stand.SetLineColor(ROOT.kBlack)

            if options.norm :
                h_njets.Scale( 1.0 / h_njets.Integral() )
                h_stand.Scale( 1.0 / h_stand.Integral() )
          
            l1.AddEntry(h_njets, getLegendEntry( wght ),"l")
            l1.AddEntry(h_stand, getLegendEntry( standName ),"l")
          
            dummyHist                   = makeDummyHistogram()
            dummyHist.GetYaxis().SetRangeUser( math.pow( 10.0, (math.floor( math.log10( h_njets.GetMaximum()))-3)+.01), math.pow( 10.0, (math.floor( math.log10( h_njets.GetMaximum()))+1)))
            dummyHist.GetXaxis().SetLabelSize( 0.0 )
            
            if snnBin == "" :
                dummyHist.SetTitle( "N_{jets} for All S_{NN} Bins" )
            else:
                dummyHist.SetTitle( "N_{jets} in S_{NN} Bin "+snnBin[1:3] )
            

            dummyHist.Draw()
            h_stand.Draw("SAME E")
            h_njets.Draw("SAME E")

            l1.Draw()
           
            h_ratio = ROOT.TH1D()
            dummyRatioHist              = makeDummyRatioHistogram()

            if wght == "htWght" or wght == "noHTWght":
                h_ratio = h_njets.Clone()
            if wght == "htWghtFlat2000" :
                h_ratio = h_njets.Clone( snnBin[1:]+"_httail" )
            if wght == "htWghtNJet7" :
                h_ratio = h_njets.Clone( snnBin[1:]+"_htnjet" )

            h_ratio.SetLineColor( ROOT.kBlack )
            h_ratio.Divide( h_ratio, h_stand )
            h_ratio.SetLineWidth( 1 )
            h_ratio.SetStats( ROOT.kFALSE )
            p2.cd()
            ROOT.gPad.SetGridy(1)
            ROOT.gPad.SetTicks( 1,1) 
            dummyRatioHist.Draw("")
            h_ratio.Draw("SAME")
            c1.Update()

            if standName == "htWght" and wght != "htWght" and wght != "noHTWght" and snnBin != "" :
                outputHistos[wght+"_"+snnBin] = copy.deepcopy( h_ratio )
           
            if options.makePlotsBool :
                outputPdfFileName               = "ratio_"+year+snnBin+"_"+wght+"_njets"
                if options.log:
                    outputPdfFileName           += "_log"
                if options.norm:
                    outputPdfFileName           += "_norm"
                c1.SaveAs(outputDir+"/"+year+"/"+outputPdfFileName+".pdf")
            c1.Clear()

            del c1, p1, p2, dummyHist, dummyRatioHist

        outputFile    = ROOT.TFile( year+"_htRatioSyst.root", "RECREATE" )
        for key in outputHistos:
            outputHistos[key].Write()
        outputFile.Close()

def plotMCDiff( outputDir, year, originalHistoDict, wghtArray, bkgdArray, snnBinArray ) :

    for snnBin in snnBinArray :

        c1, p1                      = makeCanvasAndPads()
        p1.SetLogy()
        p1.Draw()

        l1                          = ROOT.TLegend( 0.41, 0.70, 0.85, 0.85 )
        l1.SetBorderSize(0)
        l1.SetFillStyle( 0 )
        l1.SetLineWidth( 3 )
        l1.SetTextFont( 42 )
        l1.SetTextSize( 0.03 )

        dummyHist                   = makeDummyHistogram()
        dummyHist.Draw()

        histoDict = copy.deepcopy( originalHistoDict )
        
        for wght in wghtArray :            
            h_njets = histoDict[ snnBin ]["All"][ wght ]
            h_njets.SetStats(0)
            
            if options.norm :
                h_njets.Scale( 1.0 / h_njets.Integral() )
            
            if wght == "noHTWght" :
                h_njets.SetLineColor( ROOT.kBlack )
            elif wght == "htWght" :
                h_njets.SetLineColor( ROOT.kBlue - 3 )
            elif wght == "htWghtFlat2000" :
                h_njets.SetLineColor( ROOT.kGreen )
            elif wght == "htWghtNJet7" :
                h_njets.SetLineColor( ROOT.kRed )
            
            l1.AddEntry(h_njets, getLegendEntry( wght ), "l")
            dummyHist.GetYaxis().SetRangeUser( math.pow( 10.0, (math.floor( math.log10( h_njets.GetMaximum()))-3)), math.pow( 10.0, (math.floor( math.log10( h_njets.GetMaximum()))+1)))

            h_njets.Draw("SAME E")
            
        l1.Draw()

        outputPdfFileName               = year+snnBin+"_njets"
        if options.log:
            outputPdfFileName           += "_log"
        if options.norm:
            outputPdfFileName           += "_norm"
        c1.SaveAs(outputDir+"/"+year+"/"+outputPdfFileName+".pdf")

        del c1, p1, dummyHist

def getLegendEntry( wght ) :

    if wght == "noHTWght" :
        return "No H_{T} Weighting"
    if wght == "htWght" :
        return "Nominal H_{T} Weighting"
    if wght == "htWghtFlat2000" :
        return "H_{T} Weight ( const. at > 2000 GeV)"
    if wght == "htWghtNJet7" :
        return "H_{T} Weight ( no N_{jets} depend.)"
        
def makeDummyRatioHistogram() :

    dummyRatioHist                   = ROOT.TH1D( "dummyRatioHist", "dummyRatioHist", 6, 7, 13 )
    
    dummyRatioHist.GetYaxis().SetRangeUser( 0.3, 1.7 )
    
    dummyRatioHist.SetYTitle( "Ratio" )
    dummyRatioHist.SetTitleSize( 0.15, "y" )
    dummyRatioHist.SetTitleOffset( 0.5, "y" )
    dummyRatioHist.SetLabelOffset( 0.025, "x" )
    dummyRatioHist.SetLabelOffset( 0.025, "y" )
    dummyRatioHist.SetLabelSize( 0.175, "x" )
    dummyRatioHist.SetLabelSize( 0.125, "y" )
    dummyRatioHist.SetNdivisions( 404, "y" )
    dummyRatioHist.SetStats( 0 )
    dummyRatioHist.SetTitle( "" )

    return dummyRatioHist

def makeDummyHistogram() :

    dummyHist                   = ROOT.TH1D( "dummyHist", "dummyHist", 6, 7, 13 )
    dummyHist.GetYaxis().SetRangeUser( 1.0, 3.0e4 )
                
    dummyHist.SetTitle( "N_{jets} For Each H_{T} Scale Factor" )

    dummyHist.SetXTitle( "N_{jets}" )
    dummyHist.GetXaxis().SetLabelSize( 0.045 )
    dummyHist.GetXaxis().SetLabelFont( 42 )
    dummyHist.GetXaxis().SetTitleSize( 0.055 )
    dummyHist.GetXaxis().SetTitleFont( 42 )

    dummyHist.SetYTitle( "Weighted Events" )
    dummyHist.GetYaxis().SetLabelSize( 0.045 )
    dummyHist.GetYaxis().SetLabelFont( 42 )
    dummyHist.GetYaxis().SetTitleSize( 0.055 )
    dummyHist.GetYaxis().SetTitleFont( 42 )
    dummyHist.GetYaxis().SetTitleOffset( 1.15 )

    dummyHist.SetStats( 0 )

    return dummyHist

def makeCanvasAndPads() :
    c1                       = ROOT.TCanvas( "c1", "c1", 0, 0, 400, 400 )
    c1.cd()
    p1                       = ROOT.TPad( "p1", "p1", 0, 0, 1.0, 1.0 )
    if( options.log ):
        p1.SetLogy()
    p1.Draw()
    p1.cd()
    
    p1.SetFrameFillStyle( 1001 )
    p1.SetTicks()
    p1.SetFillColor( 0 ) 
    
    p1.SetLeftMargin( 0.15 )
    p1.SetTopMargin( 0.1 )
    p1.SetBottomMargin( 0.16 )
    p1.SetRightMargin( 0.10 )
        
    return c1, p1

def makeCanvasAndTwoPads() :

    c1              = ROOT.TCanvas( "c1", "c1", 0, 0, 400, 480 )
    p1              = ROOT.TPad( "p1", "p1", 0, 0.30, 1.0, 1.0 )
    p1.SetLeftMargin( 0.15 )
    p1.SetBottomMargin( 0 )
    if( options.log ):
        p1.SetLogy()
    p1.Draw()
    p2              = ROOT.TPad( "p2", "p2", 0, 0, 1.0, 0.3 )
    p2.SetTopMargin( 0 )
    p2.SetLeftMargin( 0.15 )
    p2.SetBottomMargin( 0.3 )
    p2.Draw()
    c1.Update()

    return c1, p1, p2

if __name__ == '__main__' :
    main()
