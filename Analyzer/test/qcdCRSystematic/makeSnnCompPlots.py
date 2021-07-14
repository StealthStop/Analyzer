# This script plots all the MVA distribution histograms per njet and then creates a file 
# with the ratios calculated as (MVA for all events with 7 and greater jets) / (MVA for njets = 7,8,9,10,11) 

# Input: ROOT file output from makeSnnHistograms.py (Snn distributions per N_jets)
# Output: ROOT file wiht histograms with the ratio plots (N_jets = All / N_jets = j)

# Last edited by Kelvin Mei on January 4th, 2021

#!/bin/python
import ROOT
import copy
import os.path
import array

ROOT.gROOT.SetBatch( True )

from optparse import OptionParser

parser  = OptionParser()

parser.add_option('--log', action='store_true',
                    dest='log',
                    default = False, help = 'Make all plots log on the y scale' )

parser.add_option('--pdf', action='store_true',
                    dest='pdf',
                    default = False, help = 'Print out PDF file versions of all the plots' )

(options,args) = parser.parse_args()

def main():

    yearList            = [ "2016", "2017", "2018pre", "2018post" ]

    inputDir            = "snnPlots" #default output directory for makeSnnHistograms.py

    if not os.path.exists(inputDir):
        print( "Input directory does not exist. Please run makeSnnHistograms.py first")

    outputDir           = "snnRatioPlots"
    
    if not os.path.exists( outputDir ):
        print( "Output directory does not exist, so making output directory", outputDir )
        os.makedirs( outputDir )

    if( options.pdf and not os.path.exists( outputDir+"/pdfPlots" ) ):
        print( "Create pdf plots directory", outputDir+"/pdfPlots" )
        os.makedirs( outputDir+"/pdfPlots" )
        

    fileNameArray       = [ "TT", "QCD", "Data_SingleMuon" ]
    njetArray           = [ "7", "8", "9", "10", "11" ]

    for yearTag in yearList :
        histoDict           = {}
        
        #Create dictionary with all the histograms to plot later
        for fileName in fileNameArray:

            histoDict[fileName]             = {}
            inputFileName                   = inputDir+"/"+yearTag+"_"+fileName+"_SnnPerNJets.root"
            myInputFile                     = ROOT.TFile.Open( inputFileName )
            myHisto                         = myInputFile.Get( "h_snn_all" )
            histoDict[fileName]["All"]      = copy.deepcopy( myHisto )
    
            for njet in njetArray :
                myHisto                     = myInputFile.Get( "h_snn_"+njet+"j" )
                histoDict[fileName][njet]   = copy.deepcopy( myHisto )
        
        #This function plots everything and writes it to an output file
        plotRatioPlots( yearTag, outputDir, histoDict, fileNameArray, njetArray, "All" )

def plotRatioPlots( yearTag, outputDir, histoDict, fileNameArray, njetArray, numeratorName ) :

    outputFileName = outputDir+"/"+yearTag+"_snnRatios.root"
    outputFile     = ROOT.TFile( outputFileName, "RECREATE" )
        
    for fileName in fileNameArray:
        for njet in njetArray :
            if( options.pdf ):
                c1, p1, p2 = createCanvasAndPad()
                p1.cd()
               
                leg1 = createLegend()
                   
            h_numerator = histoDict[fileName][numeratorName]
            h_numerator.SetStats( 0 )
            h_numerator.SetLineWidth( 2 )
            h_numerator.SetLineColor( ROOT.kRed )
            
            if h_numerator.Integral() != 0 :
                h_numerator.Scale( 1.0 / h_numerator.Integral() )
            else :
                print("Numerator histogram has no entries and cannot be normalized")
            
            h_denominator = histoDict[fileName][njet]
            h_denominator.SetStats( 0 )
            h_denominator.SetLineWidth( 2 )
            h_denominator.SetLineColor( ROOT.kBlue )
                
            if h_denominator.Integral() != 0 :
                h_denominator.Scale( 1.0 / h_denominator.Integral() )
            else :
                print("Denominator histogram has no entries and cannot be normalized")
            
            h_numerator.SetTitle( "S_{NN} Shape for "+fileName+" for N_{Jets} = "+njet )
            
            if( options.pdf ):
                leg1.AddEntry(h_numerator, "All N_{j}","l")
                leg1.AddEntry(h_denominator, "N_{j} = "+njet,"l")
               
                h_numerator.Draw("HIST")
                h_denominator.Draw("SAME E")
    
                leg1.Draw()
    
            h_ratio = h_numerator.Clone( "h_ratio_"+fileName+"_"+njet )
            h_ratio.SetLineColor(ROOT.kBlue)
            h_ratio.SetTitle("")
            h_ratio.GetXaxis().SetTickLength( 0.1 )
            h_ratio.GetYaxis().SetNdivisions( 12, 5, 0 )

            h_ratio.GetXaxis().SetTitle("S_{NN} Score")
            h_ratio.GetYaxis().SetTitle("Ratio: All N_{Jets} / N_{Jets} = "+njet ) 
            h_ratio.GetXaxis().SetTitleSize(0.075*.8)
            h_ratio.GetYaxis().SetTitleSize(0.075*.8)
            h_ratio.GetYaxis().SetTitleOffset( 0.9 )
            h_ratio.GetXaxis().SetLabelSize(0.08*.8)
            h_ratio.GetYaxis().SetLabelSize(0.08*.8)
            h_ratio.Divide( h_ratio, h_denominator ) 
            h_ratio.SetMinimum(0.0)
            h_ratio.SetMaximum(1.55)
            h_ratio.SetStats(ROOT.kFALSE)
            
            if( options.pdf ) :        
                p2.cd()
            
                ROOT.gPad.SetGridy(1)
                ROOT.gPad.SetTicks( 1,1) 
                h_ratio.Draw("PE")
                
                line = ROOT.TLine( -0.1, 1.0, 1.1, 1.0 )
                line.SetLineWidth( 2 )
                line.SetLineColor( ROOT.kRed )
                line.Draw("SAME")

                c1.Update()
                c1.SaveAs(outputDir+"/pdfPlots/"+yearTag+"_"+fileName+"_"+njet+"_snnRatio.pdf")
                c1.Clear()

            outputFile.cd()
            h_ratio.Write()
            
    outputFile.Close()

def createLegend() :

        leg1            = ROOT.TLegend( 0.62, 0.70, 0.89, 0.88 )

        leg1.SetFillStyle( 0 )
        leg1.SetBorderSize( 0 )
        leg1.SetLineWidth( 1 )
        leg1.SetNColumns( 1 )
        leg1.SetTextFont( 42 )
        leg1.SetTextSize( 0.055 )

        return leg1

def createCanvasAndPad() :
        c1              = ROOT.TCanvas( "c1", "c1", 0, 0, 700, 700 )
        p1              = ROOT.TPad( "p1", "p1", 0, 0.41, 1.0, 1.0 )
        p1.SetLeftMargin( 0.15 )
        p1.SetBottomMargin( 0 )
        if( options.log ):
            p1.SetLogy()
        p1.Draw()
        
        p2              = ROOT.TPad( "p2", "p2", 0, 0, 1.0, 0.4 )
        p2.SetTopMargin( 0 )
        p2.SetLeftMargin( 0.15 )
        p2.SetBottomMargin( 0.3 )
        p2.Draw()
        
        c1.Update()

        return c1, p1, p2


if __name__ == '__main__' :
    main()
