# This script takes the root files that make the per Njet ratios from makeSnnCompPlots.py
# and makes the ratio of ratio between Data CR and TT in the SR. QCD CR is added to plots
# just for comparison.

# The plots here do not use the Clopper Pearson normalization - just the normal area under
# the curve normalization.

# Input: ROOT file output from makeSnnCompPlots.py (ratio plots)
# Output: ROOT file with the ratio of ratio plots (i.e. data CR / ttbar MC SR) used to
#           to create the new QCD CR systematic

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

(options,args) = parser.parse_args()

def main():

    ROOT.TH1.SetDefaultSumw2()
        
    #Check output pdf directory
    outputDir = "ratioOfRatioPlots"
        
    if not os.path.exists(outputDir) :
        print "Output directory does not exist, so making output directory", outputDir
        os.makedirs(outputDir)
    
    if not os.path.exists(outputDir+"/pdfFiles") :
        print "Output directory for PDF plots does not exist, so making output directory", outputDir+"/pdfFiles"
        os.makedirs(outputDir+"/pdfFiles")

    yearList                = [ "2016", "2017", "2018pre", "2018post" ]
    njetList                = [ "7", "8", "9", "10", "11" ]

    histoDict               = {}
    histoDict["TT"]         = {}
    histoDict["QCD"]        = {}
    histoDict["Data"]       = {}

    for year in yearList :
        #Name the input file after what you hadded all the root files in mvaRatioPlots.
        inputFileName                   = "snnRatioPlots/"+year+"_snnRatios.root"
        inputFile                       = ROOT.TFile.Open( inputFileName )
            
        for njet in njetList :
            histoDict["TT"][njet]       = copy.deepcopy( inputFile.Get( "h_ratio_TT_"+njet ) )
            histoDict["QCD"][njet]      = copy.deepcopy( inputFile.Get( "h_ratio_QCD_"+njet ) )
            histoDict["Data"][njet]     = copy.deepcopy( inputFile.Get( "h_ratio_Data_SingleMuon_"+njet ) )
            dataHistoOverTtHisto        = histoDict["Data"][njet].Clone( "h_ratio_DataCR_over_TTSR_"+njet )
            dataHistoOverTtHisto.Divide( dataHistoOverTtHisto, histoDict["TT"][njet], 1, 1, "" )
            histoDict[njet]             = copy.deepcopy( dataHistoOverTtHisto )
        
        inputFile.Close()
    
        #Write all histograms to output file
        outputFileName = outputDir+"/"+year+"_ratioOfRatios.root"
        outputFile  = ROOT.TFile.Open( outputFileName, "RECREATE" )
        for njet in njetList :
            histoDict[njet].Write()
        outputFile.Close()
       
        # Basic plots of the njet ratios - not the ones used in the final paper. 
        # You do not need to plot these if you have the plots for Figure 3.
        plotTrio( year, outputDir, histoDict, njetList )

def plotTrio( year, outputDir, tempHistoDict, njetList ) :
    
    for njet in njetList :
        c1, p1          = createCanvasAndPad()
        p1.cd()
        leg1            = createLegend() 
        
        histoDict       = copy.deepcopy( tempHistoDict )
           
        h1              = histoDict["TT"][njet]
        h1.SetStats(0)
        h1.SetLineWidth(2)
        h1.SetLineColor(ROOT.kGreen+2)

        h2              = histoDict["QCD"][njet]
        h2.SetStats(0)
        h2.SetLineWidth(2)
        h2.SetLineColor(ROOT.kOrange+7)
        
        h3              = histoDict["Data"][njet]
        h3.SetStats(0)
        h3.SetLineWidth(2)
        h3.SetLineColor(ROOT.kBlue)
            
        h2.SetTitle( "N_{Jets} Ratio Shapes for N_{Jets} = "+njet )
        h2.GetXaxis().SetTitleSize( 0.04 )
        h2.GetYaxis().SetTitleSize( 0.04 )
        h2.GetXaxis().SetLabelSize( 0.04 )
        h2.GetYaxis().SetLabelSize( 0.04 )
        h2.GetYaxis().SetTitle("")
        h2.GetYaxis().SetTitleOffset( 1.2 )
        h2.SetMaximum(2.0)
        h2.SetMinimum(0.0)

        leg1.AddEntry(h1, "TT SR","l")
        leg1.AddEntry(h2, "QCD CR","l")
        leg1.AddEntry(h3, "Data CR","l")
            
        h2.Draw("HIST E")
        h1.Draw("HIST E SAME" )
        h3.Draw("HIST E SAME" )

        leg1.Draw()
        c1.SaveAs("ratioOfRatioPlots/pdfFiles/plotTrio_"+year+"_"+njet+"j_ratioOfRatioComp.pdf")
        c1.Clear()

def createCanvasAndPad() :
    c1              = ROOT.TCanvas( "c1", "c1", 0, 0, 1000, 1000 )
    p1              = ROOT.TPad( "p1", "p1", 0, 0, 1.0, 1.0 )
    p1.SetLeftMargin( 0.15 )
    if( options.log ):
        p1.SetLogy()
    p1.Draw()
    c1.Update()
    return c1, p1

def createLegend() :
    leg1            = ROOT.TLegend( 0.30, 0.78, 0.89, 0.88 )

    leg1.SetFillStyle( 0 )
    leg1.SetBorderSize( 0 )
    leg1.SetLineWidth( 1 )
    leg1.SetNColumns( 2 )
    leg1.SetTextFont( 42 )
    leg1.SetTextSize( 0.035 )

    return leg1

def decideColor( njet ) :
    if njet == "7" :
        #return ROOT.kRed
        return ROOT.kBlue
    if njet == "8" :
        return ROOT.kRed
    if njet == "9" :
        #return ROOT.kRed
        return ROOT.kGreen+2
    if njet == "10" :
        return ROOT.kMagenta
    if njet == "11" :
        return ROOT.kOrange+2

if __name__ == '__main__' :
    main()
