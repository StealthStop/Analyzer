# This scripts creates the HT scale factor used in the SUS-19-004 analysis.
# It also produces all the AN plots with the different HT scale factors

import ROOT
import math
import os
import numpy as np
import DataSetInfo as info

ROOT.gROOT.SetBatch( True )
ROOT.TH1.SetDefaultSumw2()

njetList        = [ 5, 6, 7, 8 ]
yearList        = [ '2016', '2017', '2018pre', '2018post' ]

def main():
    
    normCoeff   =   {}
    expoCoeff   =   {}
    
    inputDir    = "mainHtSystematic/" 
    if not os.path.exists( inputDir ) :
        print( "Directory with the input histograms do not exist" )
        return 0
    
    outputDir    = "mainHtSystematic/" #In this case, the input directory is also the output directory
    if not os.path.exists( outputDir ) :
        print( "Output directory does not exist. Creating directory: ",outputDir )
        os.makedirs( outputDir )

    for year in yearList : 
        
        for njet in njetList :
            
            #data        = info.DataSetInfo(basedir=inputDir, fileName=year+"_htStudy_Data_SingleLepton.root",   label="Data")
            #bg          = info.DataSetInfo(basedir=inputDir, fileName=year+"_htStudy_AllBG.root",  label="Background")
            
            data        = info.DataSetInfo(basedir=inputDir, fileName=year+"_Data.root",   label="Data")
            bg          = info.DataSetInfo(basedir=inputDir, fileName=year+"_AllBG.root",  label="Background")
    
            histoName   = 'h_ht_1l_'+str(njet)+'j_ge1b'
    
            dataHist    = data.getHisto( histoName )
            bgHist      = bg.getHisto( histoName )
            
            # The original histograms were binned too finely so had to rebin
            dataHist.Rebin(10)
            bgHist.Rebin(10)

            # Histograms from Analyze1Lep are too finely binned by double
            # dataHist.Rebin(2)
            # bgHist.Rebin(2)
            
            dataHist.GetXaxis().SetRange(1,30)
            bgHist.GetXaxis().SetRange(1,30)
            
            makeBasicDataBgPlots( outputDir, year, histoName, dataHist, bgHist )
    
            # Create the ratio histogram
            c1, p1      = createCanvasAndPad() 
            ratioHist   = dataHist.Clone( "Ratio" )
            ratioHist.Divide( bgHist )
            ratioHist   = formatRatioHist( ratioHist, histoName )
            ratioHist.Draw("EP")
            
            line        = ROOT.TF1("1" ,"1" ,-2000,20000)
            line.SetLineColor( ROOT.kBlack )
            line.Draw( "SAME" )
            
            fit         = ROOT.TF1("fit", "[0]*exp([1]*x)", 300.0, 3000.0, 2)
            fit.SetParLimits(0,  0.0, 5.0)
            fit.SetParLimits(1,  -5.0, 5.0)
            fit.SetLineColor(ROOT.kRed)
            ratioHist.Fit(fit, "", "", 300.0, 3000.0)
            fit.Draw( "SAME" )
    
            normCoeff[njet]  = fit.GetParameter(0) 
            expoCoeff[njet]  = 1000*fit.GetParameter(1)
    
            c1.SaveAs( outputDir+year+"_"+histoName+"_ratio.pdf" )
            del c1 
    
        ####################################################
        #Fit results of all the fits
        ####################################################
        
        # Make np arrays since TGraph can only use np.arrays
        nJets               = np.array( [ 5, 6, 7 ] )
        norm                = np.array( [ normCoeff[i] for i in nJets ] )
        expo                = np.array( [ expoCoeff[i] for i in nJets ] )
    
        s = ROOT.TGraph(nJets.size, nJets.astype(np.double), norm.astype(np.double))
        y = ROOT.TGraph(nJets.size, nJets.astype(np.double), expo.astype(np.double))
       
        paramDict = makeLinearFitPlots( s, y, outputDir, year )
    
        ####################################################
        #Plot Scale Factor per nJet
        ####################################################
        c1, p1              = createCanvasAndPad()
        leg                 = createLegend( 0.6, 0.65, 0.9, 0.88 )
    
        h = ROOT.TH1F("H_{T} Scale Factor", "H_{T} Scale Factor", 1,0,3000)
        h.GetXaxis().SetTitle("H_{T} [GeV]")
        h.GetYaxis().SetTitle("Scale Factor")
        h.SetMaximum(1.6)
        h.Draw()
    
        dummyList = []
        for nJet in [(ROOT.kCyan+1, 7) , (ROOT.kCyan+2, 8), (ROOT.kCyan+3, 9), (ROOT.kCyan+4, 10), (ROOT.kRed, 11), (ROOT.kRed+1, 12)]:       
            
            par     = [ paramDict['Norm']['p0']*nJet[1]+paramDict['Norm']['p1'], paramDict['Exp']['p0']*nJet[1]+paramDict['Exp']['p1'] ] 
            func    = '{}*exp( {}*(x/1000) )'.format( par[0], par[1] )
            sf = ROOT.TF1("sf"+str(nJet[1]), func, 0.0, 3000)
    
            sf.SetLineColor(nJet[0])
            sf.DrawCopy("L SAME")
    
            hd = ROOT.TH1F("HT Scale Factor"+str(nJet[1]), "HT Scale Factor"+str(nJet[1]), 1, 0, 3000)
            hd.Draw("P SAME")
            hd.SetLineColor(nJet[0])
            hd.SetLineWidth(3)
            leg.AddEntry(hd, "N_{jets} = "+str(nJet[1]), "l")
            dummyList.append(hd)
    
        leg.Draw()
        c1.SaveAs(outputDir+year+"_sf.pdf")
    
        del c1
    
        ####################################################
        #Plot 8 Jet bin extrapolation comparison with direct fit
        ####################################################
        c1, p1          = createCanvasAndPad()
        leg             = createLegend( 0.2, 0.7, 0.85, 0.88 )
    
        h = ROOT.TH1F("H_{T} Scale Factor: N_{jets} = 8", "H_{T} Scale Factor: N_{jets} = 8", 1,0,3000)
        h.GetXaxis().SetTitle("H_{T} [GeV]")
        h.GetYaxis().SetTitle("Scale Factor")
        h.SetMaximum( 2.0 )
        h.Draw()
        
        nJet    = (ROOT.kCyan+2, "8")
        par     = [ paramDict['Norm']['p0']*8+paramDict['Norm']['p1'], paramDict['Exp']['p0']*8+paramDict['Exp']['p1'] ] 
        func    = '{}*exp( {}*(x/1000) )'.format( par[0], par[1] )
        sf = ROOT.TF1("a", func, 0.0, 3000)
        sf.SetLineColor(ROOT.kRed)
        leg.AddEntry(sf, "Extrapolation: N_{jets} = 8", "l")
        sf.DrawCopy("L SAME")
        
        func2   = '{}*exp( {}*(x/1000) )'.format( normCoeff[8], expoCoeff[8] )
        sffit = ROOT.TF1("b", func2, 0.0, 3000)
        sffit.SetLineColor(ROOT.kBlack)
        leg.AddEntry(sffit, "Fit: N_{jets} = 8", "l")
        sffit.DrawCopy("L SAME")
        
        leg.Draw()
        c1.SaveAs(outputDir+year+"_sf_8Jetbin.pdf")
        del c1
    
        ####################################################
        #Plot Scale Factor Up and Down variation
        ####################################################
        c1, p1          = createCanvasAndPad()
        leg             = createLegend( 0.5, 0.8, 0.85, 0.88 )
        
        h = ROOT.TH1F("HT Scale Factor Variation", "H_{T} Scale Factor Variation", 1,0,3000)
        h.GetXaxis().SetTitle("H_{T} [GeV]")
        h.GetYaxis().SetTitle("A.U.")
        h.SetMinimum(0.0)
        h.SetMaximum(2.0)
        h.Draw()
        
        ratio = ROOT.TF1("asf", "b/a", 0.0, 3000)
        ratio.DrawCopy("L SAME")
        leg.AddEntry(ratio, "Variation Up", "l")
        
        ratio2 = ROOT.TF1("yep","a/b", 0.0, 3000)
        ratio2.SetLineColor(ROOT.kBlack)
        ratio2.DrawCopy("L SAME")
        leg.AddEntry(ratio2, "Variation Down", "l")
    
        leg.Draw()
        c1.SaveAs(outputDir+year+"_sf_8Jetbin_Up_Down.pdf")
        del c1

def makeLinearFitPlots( s, y, outputDir, year ) :
    plotConstraints = [ ("Norm Term",s, -2.0, 2.0, -2, 5, 3, -1.0, 5.0),
                        ("Exp Term",y, -2.0, 2.0, -5, 5, -1, -10.0, 4.0) ]

    paramDict       = {}

    for l in plotConstraints:
        print( l[0] )
        c1 = ROOT.TCanvas( "c1", "c1", 0, 0, 800, 800)
        ROOT.gPad.SetLeftMargin(0.12)
        ROOT.gPad.SetRightMargin(0.15)
        ROOT.gPad.SetTopMargin(0.08)
        ROOT.gPad.SetBottomMargin(0.12)
        ROOT.gPad.SetTicks(1,1)
        ROOT.TH1.SetDefaultSumw2()
        ROOT.gStyle.SetOptFit();
        
        l[1].SetMaximum(l[3])
        l[1].SetMinimum(l[2])
        l[1].SetLineColor(ROOT.kBlack)
        l[1].SetMarkerColor(ROOT.kBlack)
        l[1].SetMarkerStyle(21)
        l[1].SetTitle(l[0]+" Fits")
        l[1].GetXaxis().SetTitle("NJet")
        l[1].GetYaxis().SetTitle(l[0])
        l[1].Draw("AP")
        
        fit = ROOT.TF1("fit", "[0]*x + [1]", 0.0, 10.0, 3)
        fit.SetParLimits(0,  l[4], l[5])
        fit.SetParLimits(1,  l[7], l[8])
        fit.SetParameter(0,  l[6])
        fit.SetParameter(1,  0.0)        
        fit.SetLineColor(ROOT.kRed)
        l[1].Fit(fit, "M", "", 0.0, 10.0)

        paramDict[l[0][:-5]] = {}
        param = fit.GetParameters()
        paramDict[l[0][:-5]]["p0"] = param[0]
        paramDict[l[0][:-5]]["p1"] = param[1]
        fit.Draw("SAME")
        
        c1.SaveAs(outputDir+year+"_"+l[0][:-5]+"_fits.pdf")
        del c1           

    return paramDict

def makeBasicDataBgPlots( outputDir, year, histoName, dataHist, bgHist ) :
    c1, p1      = createCanvasAndPad() 
    leg         = createLegend()

    bgHist.SetLineColor( ROOT.kRed )
    leg.AddEntry( dataHist, "Data", 'l' )
    leg.AddEntry( bgHist, "All MC", 'l' )
    bgHist.Draw("")
    dataHist.Draw("SAME")
    leg.Draw()
    c1.SaveAs( outputDir+year+"_"+histoName+".pdf")
        
    p1.SetLogy()
    bgHist.DrawNormalized("")
    dataHist.DrawNormalized("SAME")
    leg.Draw()
    c1.SaveAs( outputDir+year+"_"+histoName+"_log_norm.pdf")
    del c1, p1, leg

def formatRatioHist( ratioHist, histoName ) :

    ratioHist.SetMaximum( 2.0 )
    ratioHist.SetMinimum( 0.0 )
    ratioHist.SetLineColor( ROOT.kBlack )
    ratioHist.SetMarkerColor( ROOT.kBlack )
    ratioHist.SetMarkerStyle( 21 )
    ratioHist.SetTitle( histoName+" Ratio" )
    ratioHist.GetYaxis().SetTitle(" Data/Simulation ")
    ratioHist.GetXaxis().SetTitle( "H_{T} [GeV]" )
    ratioHist.SetStats(0)
    return ratioHist

def createLegend( x0=0.6, y0=0.75, x1=0.8, y1=0.88) :
    leg = ROOT.TLegend( x0, y0, x1, y1)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetLineWidth(1)
    leg.SetNColumns(1)
    leg.SetTextFont(42)
    return leg

def createCanvasAndPad() :
    c1          = ROOT.TCanvas( "c", "c", 0, 0, 800, 800 )
    c1.cd() 
    p1          = ROOT.TPad( "p", "p", 0, 0, 1.0, 1.0 )
    p1.Draw()
    p1.cd()

    p1.SetLeftMargin(0.12)
    p1.SetRightMargin(0.08)
    p1.SetTopMargin(0.08)
    p1.SetBottomMargin(0.12)
    p1.SetTicks(1,1)
    
    ROOT.TH1.SetDefaultSumw2()
    ROOT.gStyle.SetOptStat(0)

    return c1, p1

if __name__ == '__main__':
    main()
