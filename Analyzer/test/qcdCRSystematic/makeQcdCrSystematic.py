# This script creates the final distributions used to generate the fit
#   parameters for the CR data-driven systematic

# Input: ROOT files containing all of the toys from throwToysForSystematic.py and
#       the njets_for_Aron.root file for each year
# Output: Histograms per year that serve as the baseline for the CR data-driven 
#   systematic.

# Last edited by Kelvin on January 12th, 2021

#!/bin/python
import ROOT
import copy
import os.path
import array
import math
import numpy as np

ROOT.gROOT.SetBatch( True )

from optparse import OptionParser

parser = OptionParser()

# Choose whether to add the full deviation from 1 to the error

parser.add_option( '--noFullDev', action = 'store_true',
                    dest = 'noFullDev',
                    default = False, help = 'Default to adding full deviation from 1 to the stat error')

(options, args) = parser.parse_args()

def main():

    ROOT.TH1.SetDefaultSumw2()
    ROOT.TH2.SetDefaultSumw2()

    inputDir                                = "toyPlots"
    if not os.path.exists(inputDir):
        print "Directory where all the toy njets histograms does not exist"

    # Define location of the njets_for_Aron.root files
    #   For all intents and purposes, I copied them all over with the year added as a prefix
    #   i.e. 2016_njets_for_Aron.root
    ttFileDir                               = "."
    if not os.path.exists(ttFileDir):
        print "Directory where the njets_for_Aron.root files are does not exist"


    yearList                                = [ "2016" ]#, "2017", "2018pre", "2018post" ]

    snnBinArray                             = [ "D1", "D2", "D3", "D4" ]
    njetNameArray                           = [ "7", "8", "9", "10", "11", "12" ]
    nToys                                   = 100

    for year in yearList :
        #Get the ttbar nominal shape systematic from the njets_for_Aron.root file by taking the actual njets shape and doing the division
        ttFileName                          = ttFileDir+"/"+year+"_njets_for_Aron.root"
        ttFile                              = ROOT.TFile.Open( ttFileName )
    
        ttHistoAll                          = ROOT.TH1D()
        ttHisto                             = {}
    
        #Make output directories for plots
        if not os.path.exists("toyPdfs"):
            print "Making output directory toyPdfs"
            os.makedirs("toyPdfs")
        
        if not os.path.exists("njetPdfs"):
            print "Making output directory njetPdfs"
            os.makedirs("njetPdfs")
    
        for snn in snnBinArray:
            tempHisto                       = ttFile.Get( snn+"_TT_h_njets_pt30_1l" )
            ttHisto[snn]                    = copy.deepcopy( tempHisto.Clone( snn+"_TT" ) )
    
            ttHisto[snn].Scale( 1.0 / ttHisto[snn].Integral() )
            
            if snn == "D1" :
                ttHistoAll                  = copy.deepcopy( tempHisto.Clone( "All_TT" ) )
            else :
                ttHistoAll.Add( tempHisto )
    
        ttHistoAll.Scale( 1.0 / ttHistoAll.Integral() )
    
        for snn in snnBinArray :
            ttHisto[snn].Divide( ttHisto[snn], ttHistoAll, 1, 1, "" )
        
        ttFile.Close()
    
        #Calculate the error based on all the toys - input file is the output of throwToysForSystematic.py
        toyFileName                         = inputDir+"/"+year+"_evtWghtSyst_AllToys.root"
        toyFile                             = ROOT.TFile.Open( toyFileName )
    
        toyHisto                            = {} #Holds all the histograms for each toy
    
        valHisto                            = {}
        errHisto                            = {}
    
        for snn in snnBinArray :
            valHisto[snn]                   = {}
            errHisto[snn]                   = {}
    
            for njet in njetNameArray :
                valHisto[snn][njet]         = []
                errHisto[snn][njet]         = []
    
    
        #Everything had to be rebinned because we were using the shifted histograms from njets_for_Aron.root (0 to 6) file and it would not divide properly with the not shifted histograms from throwToysForSystematic.py (7 to 12) - hence the seemingly redundant histograms.
        for i in range( nToys ) :
            toyHisto[i]                     = {}
            toyHisto[i]["All"]              = ROOT.TH1D( "h_All_shift_toy"+str(i), "h_All_shift_toy"+str(i), 6, 0, 6 )
            for snn in snnBinArray :
                toyHisto[i][snn]            = ROOT.TH1D( "h_"+snn+"_shift_toy"+str(i), "h_"+snn+"_shift_toy"+str(i), 6, 0, 6 )
    
            for snn in snnBinArray :
                tempHisto                   = toyFile.Get( "h_"+snn+"_toy"+str(i) )
                for binIter in range( 0, tempHisto.GetNbinsX() + 1 ) :
                    toyHisto[i][snn].SetBinContent( binIter, tempHisto.GetBinContent( binIter ) )
                    toyHisto[i][snn].SetBinError( binIter, tempHisto.GetBinError( binIter ) )
    
            for snn in snnBinArray :
                toyHisto[i]["All"].Add( toyHisto[i][snn] )
                toyHisto[i][snn].Scale( 1.0 / toyHisto[i][snn].Integral() )
    
            toyHisto[i]["All"].Scale( 1.0 / toyHisto[i]["All"].Integral() )
    
            for snn in snnBinArray :
                #Divide each njets shape for each Snn bin by the njets shape for all the Snn bins
                toyHisto[i][snn].Divide( toyHisto[i][snn], toyHisto[i]["All"], 1, 1, "" )
                #Divide the resulting ratio with the normalized shape from ttbar nominal (like all systematics)
                toyHisto[i][snn].Divide( toyHisto[i][snn], ttHisto[snn], 1, 1, "" )
    
        
        for snn in snnBinArray :
            for njet in njetNameArray :
                for i in range( nToys ) :
                    valHisto[snn][njet].append( toyHisto[i][snn].GetBinContent( int(njet) - 6 ) )
                    errHisto[snn][njet].append( toyHisto[i][snn].GetBinError( int(njet) - 6 ) )
       
        #The next few lines of code were just added to visualize the spread of the value for number of events per Snn bin and per njet bin (may not be useful to see anymore)
        valHistoArray                       = {}
        errHistoArray                       = {}
    
        for snn in snnBinArray :
            valHistoArray[snn]              = {}
            errHistoArray[snn]              = {}
            for njet in njetNameArray :
                valHistoArray[snn][njet]    = ROOT.TH1D( "h_"+snn+"_"+njet+"j", "h_"+snn+"_"+njet+"j", 100, np.mean( valHisto[snn][njet] ) - 3*np.std( valHisto[snn]["12"] ), np.mean( valHisto[snn][njet] ) + 3*np.std( valHisto[snn]["12"] ) )
                errHistoArray[snn][njet]    = ROOT.TH1D( "h_"+snn+"_"+njet+"j_err", "h_"+snn+"_"+njet+"j_err", 300, np.mean( errHisto[snn][njet] ) - .5, np.mean( errHisto[snn][njet] ) + .5 )
    
        for snn in snnBinArray :
            for njet in njetNameArray :
                for i in range( len(valHisto[snn][njet]) ) :
                    valHistoArray[snn][njet].Fill( valHisto[snn][njet][i] )
                    errHistoArray[snn][njet].Fill( errHisto[snn][njet][i] )
    
        outputFileName                  = year+"_QcdCrSyst_StatPlusFullDev.root"
        outputFile                      = ROOT.TFile( outputFileName, "RECREATE" ) 
        for snn in snnBinArray:
            qcdCRHisto                  = ROOT.TH1D( snn+"_qcdCR", snn+"_qcdCR", 6, 0, 6 )
            qcdCRErrHisto               = ROOT.TH1D( snn+"_qcdCRErr", snn+"_qcdCRErr", 6, 0, 6 )
            
            for binIter in xrange( 1, 7 ):
                errValue                = np.std( valHisto[snn][str(binIter+6)] ) #Base error is the standard deviation for each Snn bin each njet bin for the 100 toys
                
                if options.noFullDev :
                    errValueAdded           = 0.0 #No full deviation from 1
                else :
                    errValueAdded           = ( np.mean( valHisto[snn][str(binIter+6)] ) - 1 ) #full deviation from one - will be squared in any meaningful calculation below
    
                #print snn,binIter+6,np.mean( errHisto[snn][str(binIter+6)] ), np.std( valHisto[snn][str(binIter+6)] ) #Just used to compare the size of the average error in each toy with the standard deviation of the toy values
                qcdCRHisto.SetBinContent( binIter, np.mean( valHisto[snn][str(binIter+6)] ) )
                qcdCRHisto.SetBinError( binIter, math.sqrt( errValue**2 + errValueAdded**2 ) )
                qcdCRErrHisto.SetBinContent( binIter, 1 + math.sqrt( errValue**2 + errValueAdded**2 )/ qcdCRHisto.GetBinContent( binIter ) )
                
            qcdCRHisto.Write()
            qcdCRErrHisto.Write()
    
            #Plot error histograms
            c1              = ROOT.TCanvas( "c1", "c1", 0, 0, 1000, 1000 )
            p1              = ROOT.TPad( "p1", "p1", 0, 0, 1.0, 1.0 )
            p1.SetLeftMargin( 0.15 )
            p1.SetLogy(False)
            p1.Draw()
            p1.SetGridy(1)
            c1.Update()
        
            p1.cd()
            leg1            = ROOT.TLegend( 0.30, 0.78, 0.89, 0.88 )
            leg1.SetFillStyle( 0 )
            leg1.SetBorderSize( 0 )
            leg1.SetLineWidth( 1 )
            leg1.SetNColumns( 2 )
            leg1.SetTextFont( 42 )
            leg1.SetTextSize( 0.035 )
               
            qcdCRHisto.SetStats(0)
            qcdCRHisto.SetLineWidth(2)
            qcdCRHisto.SetLineColor(decideColorSnn(snn))
                
            qcdCRHisto.SetTitle( "Data-driven CR Systematic Values "+snn+" in "+year )
            qcdCRHisto.GetXaxis().SetTitleSize( 0.04 )
            qcdCRHisto.GetYaxis().SetTitleSize( 0.04 )
            qcdCRHisto.GetXaxis().SetLabelSize( 0.04 )
            qcdCRHisto.GetYaxis().SetLabelSize( 0.04 )
            qcdCRHisto.GetYaxis().SetTitle("")
            qcdCRHisto.GetYaxis().SetTitleOffset( 1.2 )
            qcdCRHisto.GetYaxis().SetRangeUser( 0.5, 1.5 )
    
            qcdCRHisto.Draw("HIST E")
            c1.SaveAs("njetPdfs/"+year+"_ratio_"+snn+".pdf")
           
            for snn in snnBinArray :
                for njet in njetNameArray :
                    c2              = ROOT.TCanvas( "c2", "c2", 0, 0, 1000, 1000 )
                    p2              = ROOT.TPad( "p2", "p2", 0, 0, 1.0, 1.0 )
                    p2.SetLeftMargin( 0.15 )
                    p2.SetLogy(False)
                    p2.Draw()
                    p2.SetGridy(1)
                    c2.Update()
                
                    p2.cd()
                    leg1            = ROOT.TLegend( 0.30, 0.78, 0.89, 0.88 )
                    leg1.SetFillStyle( 0 )
                    leg1.SetBorderSize( 0 )
                    leg1.SetLineWidth( 1 )
                    leg1.SetNColumns( 2 )
                    leg1.SetTextFont( 42 )
                    leg1.SetTextSize( 0.035 )
            
                    tempHisto       = valHistoArray[snn][njet]
                    tempHisto.SetLineWidth(2)
                    tempHisto.SetLineColor(decideColorSnn(snn))
                        
                    tempHisto.SetTitle( "Data-driven CR Systematic Value for 100 Toys in "+snn+" for "+njet+"j in "+year )
                    tempHisto.GetXaxis().SetTitleSize( 0.04 )
                    tempHisto.GetYaxis().SetTitleSize( 0.04 )
                    tempHisto.GetXaxis().SetLabelSize( 0.04 )
                    tempHisto.GetYaxis().SetLabelSize( 0.04 )
                    tempHisto.GetYaxis().SetTitle("")
                    tempHisto.GetYaxis().SetTitleOffset( 1.2 )
            
                    tempHisto.Draw("HIST E")
                    c2.SaveAs("toyPdfs/"+year+"_toyR_"+snn+"_"+njet+"j.pdf")
                    
                    tempHisto2       = errHistoArray[snn][njet]
                    tempHisto2.SetLineWidth(2)
                    tempHisto2.SetLineColor(decideColorSnn(snn))
                        
                    tempHisto2.SetTitle( "Data-driven CR Systematic Error Bar Values for 100 Toys in "+snn+" for "+njet+"j in "+year )
                    tempHisto2.GetXaxis().SetTitleSize( 0.04 )
                    tempHisto2.GetYaxis().SetTitleSize( 0.04 )
                    tempHisto2.GetXaxis().SetLabelSize( 0.04 )
                    tempHisto2.GetYaxis().SetLabelSize( 0.04 )
                    tempHisto2.GetYaxis().SetTitle("")
                    tempHisto2.GetYaxis().SetTitleOffset( 1.2 )
            
                    tempHisto2.Draw("HIST E")
                    #c2.SaveAs("toyPdfs/"+year+"_toyRprime_"+snn+"_"+njet+"j.pdf")

def decideColorSnn( snn ) :
    if snn == "D1" :
        return ROOT.kBlue
    if snn == "D2" :
        return ROOT.kRed
    if snn == "D3" :
        return ROOT.kGreen+2
    if snn == "D4" :
        return ROOT.kMagenta

if __name__ == '__main__' :
    main()
