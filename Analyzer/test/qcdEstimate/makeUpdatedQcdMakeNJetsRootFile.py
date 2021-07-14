# This is the script that takes the QCD.root file produced by the MakeNJets analyzer and creates the new QCD.root updated file 
# that uses the data in the control region as the actual estimate for QCD. If there are more QCD based systematics, those need
# to be added to the "tagList"

# The original QCD_MC files (from MakeNJets analyzer) are copied to this directory as 201*_QCD_MC.root
# This also outputs the percent error on the transfer factor (as well as the transfer factor components) to the screen.

# Last edited by Kelvin Mei on December 30th, 2020

#!/bin/python
import os
import math
import ROOT
import array
import random
import copy

ROOT.gROOT.SetBatch(True)

from optparse import OptionParser

parser = OptionParser()

(options, args) = parser.parse_args()

def main() :
    yearList            = [ "2016", "2017", "2018pre", "2018post" ]

    for yearTag in yearList :
        print "                    \n\n",yearTag,"\n\n"
        inputDir            = "njetRootHistos"
        if not os.path.exists(inputDir):
            print 'Input directory',inputDir,' does not exist'
       
        mvaBinList          = [ "", "_D1", "_D2", "_D3", "_D4" ]
        shiftList           = [ "njetsShifted" ]
        tagList             = [ "",
                                "_btgUp", "_btgDown",
                                "_lepUp", "_lepDown",
                                "_isrUp", "_isrDown",
                                "_fsrUp", "_fsrDown",
                                "_isr2Up", "_isr2Down",
                                "_fsr2Up", "_fsr2Down",
                                "_pdfUp", "_pdfDown",
                                "_htUp", "_htDown",
                                "_puUp", "_puDown",
                                "_sclUp", "_sclDown",
                                "_JECUp", "_JECDown",
                                "_JERUp", "_JERDown"
                            ]
    
    
        #Derive the Data/MC Scale Factor - first look at the Data in the control region
        dataCrFileName              = inputDir+"/Data_QcdOnly.root"
        dataCrFile                  = ROOT.TFile.Open( dataCrFileName )
        
        qcdCrDataShapeHistoDict     = {}
        qcdCrDataIntegralValues     = {}
    
        qcdCrDataShapeHistoDict["_D1"] = copy.deepcopy( dataCrFile.Get("h_njets_CR_"+yearTag+"_SRBE_mu_D1") )
        qcdCrDataIntegralValues["_D1"] = qcdCrDataShapeHistoDict["_D1"].Integral()
    
        qcdCrDataShapeHistoDict["_D2"] = copy.deepcopy( dataCrFile.Get("h_njets_CR_"+yearTag+"_SRBE_mu_D2") )
        qcdCrDataIntegralValues["_D2"] = qcdCrDataShapeHistoDict["_D2"].Integral()
    
        qcdCrDataShapeHistoDict["_D3"] = copy.deepcopy( dataCrFile.Get("h_njets_CR_"+yearTag+"_SRBE_mu_D3") )
        qcdCrDataIntegralValues["_D3"] = qcdCrDataShapeHistoDict["_D3"].Integral()
    
        qcdCrDataShapeHistoDict["_D4"] = copy.deepcopy( dataCrFile.Get("h_njets_CR_"+yearTag+"_SRBE_mu_D4") )
        qcdCrDataIntegralValues["_D4"] = qcdCrDataShapeHistoDict["_D4"].Integral()
    
        qcdCrDataShapeHistoDict["All"] = copy.deepcopy( qcdCrDataShapeHistoDict["_D1"] )
        qcdCrDataShapeHistoDict["All"].Add( qcdCrDataShapeHistoDict["_D2"] )
        qcdCrDataShapeHistoDict["All"].Add( qcdCrDataShapeHistoDict["_D3"] )
        qcdCrDataShapeHistoDict["All"].Add( qcdCrDataShapeHistoDict["_D4"] )
        qcdCrDataIntegralValues["All"] = qcdCrDataShapeHistoDict["All"].Integral()
    
        dataCrFile.Close()
        
        mcCrFileName                = inputDir+"/QCD.root"
        mcCrFile                    = ROOT.TFile.Open( mcCrFileName )
        
        qcdCrMCShapeHistoDict       = {}
        qcdCrMCIntegralValues       = {}
    
        qcdCrMCShapeHistoDict["_D1"] = copy.deepcopy( mcCrFile.Get("h_njets_CR_"+yearTag+"_SRBE_mu_D1") )
        qcdCrMCIntegralValues["_D1"] = qcdCrMCShapeHistoDict["_D1"].Integral()
    
        qcdCrMCShapeHistoDict["_D2"] = copy.deepcopy( mcCrFile.Get("h_njets_CR_"+yearTag+"_SRBE_mu_D2") )
        qcdCrMCIntegralValues["_D2"] = qcdCrMCShapeHistoDict["_D2"].Integral()
    
        qcdCrMCShapeHistoDict["_D3"] = copy.deepcopy( mcCrFile.Get("h_njets_CR_"+yearTag+"_SRBE_mu_D3") )
        qcdCrMCIntegralValues["_D3"] = qcdCrMCShapeHistoDict["_D3"].Integral()
    
        qcdCrMCShapeHistoDict["_D4"] = copy.deepcopy( mcCrFile.Get("h_njets_CR_"+yearTag+"_SRBE_mu_D4") )
        qcdCrMCIntegralValues["_D4"] = qcdCrMCShapeHistoDict["_D4"].Integral()
    
        qcdCrMCShapeHistoDict["All"] = copy.deepcopy( qcdCrMCShapeHistoDict["_D1"] )
        qcdCrMCShapeHistoDict["All"].Add( qcdCrMCShapeHistoDict["_D2"] )
        qcdCrMCShapeHistoDict["All"].Add( qcdCrMCShapeHistoDict["_D3"] )
        qcdCrMCShapeHistoDict["All"].Add( qcdCrMCShapeHistoDict["_D4"] )
        qcdCrMCIntegralValues["All"] = qcdCrMCShapeHistoDict["All"].Integral()

        print yearTag,"QCD CR MC Integral: ",qcdCrMCIntegralValues["All"]
    
        mcCrFile.Close()
        
        #Now take the histograms from the MakeNJets file output by the analyzer
        makeNjetsForAronFileName    = yearTag+"_QCD_MC.root"
        makeNjetsForAronFile        = ROOT.TFile.Open( makeNjetsForAronFileName )
    
        qcdSrAllNjetsHisto          = makeNjetsForAronFile.Get( "h_njetsShifted_pt30_1l" )
        qcdSrAllNjetsIntegralValue  = qcdSrAllNjetsHisto.Integral()

        print yearTag,"QCD SR MC Integral: ",qcdSrAllNjetsHisto.Integral()
    
        #The overall scale factor is all the QCD events in the signal region divided by all the QCD events in the control region
        srCrSf                      = qcdSrAllNjetsIntegralValue / qcdCrMCIntegralValues["All"]
    
        srIntegralValueErrSq        = 0.0
        crIntegralValueErrSq        = 0.0
    
        for binIter in xrange( 1, qcdCrMCShapeHistoDict["All"].GetNbinsX() + 1) :
            print "QCD CR MC Bin Content(",binIter,"): ",qcdCrMCShapeHistoDict["All"].GetBinError(binIter),";\tQCD SR MC BinContent(",binIter,"): ",qcdSrAllNjetsHisto.GetBinError(binIter)
            srIntegralValueErrSq    += qcdSrAllNjetsHisto.GetBinError(binIter)**2
            crIntegralValueErrSq    += qcdCrMCShapeHistoDict["All"].GetBinError(binIter)**2
            
            print "Numerator (SR): ",srIntegralValueErrSq, math.sqrt(qcdSrAllNjetsHisto.GetBinError(binIter)**2)
            print "Denominator (CR): ",crIntegralValueErrSq, math.sqrt(qcdCrMCShapeHistoDict["All"].GetBinError(binIter)**2)
        
        print "CR Error Squared: ",math.sqrt(crIntegralValueErrSq),"\nSR Error Squared : ",math.sqrt(srIntegralValueErrSq),"\nPercent Error Total: ",math.sqrt(crIntegralValueErrSq/(qcdCrMCIntegralValues["All"]**2)+srIntegralValueErrSq/(qcdSrAllNjetsIntegralValue**2))
        print qcdCrMCIntegralValues["All"], qcdSrAllNjetsIntegralValue
    
        newQcdHistos                = []
    
        for iBin in mvaBinList :
            for iShift in shiftList :
                
                oldQcdHisto         = makeNjetsForAronFile.Get( "h_"+iShift+"_pt30_1l"+iBin )
                newQcdHisto         = ROOT.TH1D()
    
                if iBin == "" :
                    newQcdHisto     = qcdCrDataShapeHistoDict["All"].Clone( "h_"+iShift+"_pt30_1l"+iBin )
                else:
                    newQcdHisto     = qcdCrDataShapeHistoDict[iBin].Clone( "h_"+iShift+"_pt30_1l"+iBin )
                newQcdHisto.Scale( srCrSf )
                
                for iTag in tagList :
                    oldQcdSysHisto  = makeNjetsForAronFile.Get( "h_"+iShift+"_pt30_1l"+iBin+iTag )
                    newQcdSysHisto  = oldQcdSysHisto.Clone( "h_"+iShift+"_pt30_1l"+iBin+iTag )
                    newQcdSysHisto.SetTitle( "h_"+iShift+"_pt30_1l"+iBin+iTag )
    
                    for iNjBin in xrange( 1, oldQcdSysHisto.GetNbinsX()+2 ) :
                        perBinDiff  = oldQcdSysHisto.GetBinContent(iNjBin) - oldQcdHisto.GetBinContent(iNjBin)
                        newQcdSysHisto.SetBinContent( iNjBin, newQcdHisto.GetBinContent( iNjBin ) + perBinDiff )
                        newQcdSysHisto.SetBinError( iNjBin, newQcdHisto.GetBinError( iNjBin ) ) 
                    
                    newQcdHistos.append( copy.deepcopy( newQcdSysHisto ) )
        
        makeNjetsForAronFile.Close()
    
        myOutputFileName = yearTag+"_QCD_Updated.root"
        myOutputFile = ROOT.TFile.Open( myOutputFileName, "RECREATE" )
        for histo in newQcdHistos :
            histo.Write()
        myOutputFile.Close()

if __name__ == '__main__' :
    main()
