import ROOT
import math
import DataSetInfo as info

def main():
    path = "condor/Analyze1Lep_Kerasv1.2.0_MCTrigger_bTag_leptonWeight_ht300/"
    #path = "condor/Analyze1Lep_Kerasv1.2.0_MCTrigger_bTag_leptonWeight/"
    #path = "condor/Analyze1Lep_Kerasv1.2.0_noMCTrigger_PU_bTag_leptonWeight/"
    #histoName = "h_mbl_1l_ge7j_ge1b_noMbl"
    histoName = "h_njets_mbl_1l_ge7j_ge1b_noMbl"
    rebinVal = 5

    bgData = {
        "TT"              : info.DataSetInfo(basedir=path, fileName="TT.root",              sys=0.1, label="Background"),
        #"TTJets"          : info.DataSetInfo(basedir=path, fileName="TTJets.root",          sys=0.3, label="Background"),
        "QCD"             : info.DataSetInfo(basedir=path, fileName="QCD.root",             sys=0.5, label="Background"),
        #"DYJetsToLL_M-50" : info.DataSetInfo(basedir=path, fileName="DYJetsToLL_M-50.root", sys=1.0, label="Background"),
        #"Rare"            : info.DataSetInfo(basedir=path, fileName="Rare.root",            sys=1.0, label="Background"),
        #"Diboson"         : info.DataSetInfo(basedir=path, fileName="Diboson.root",         sys=1.0, label="Background"),
        #"ST"              : info.DataSetInfo(basedir=path, fileName="ST.root",              sys=1.0, label="Background"),
        #"WJetsToLNu"      : info.DataSetInfo(basedir=path, fileName="WJetsToLNu.root",      sys=1.0, label="Background"),
    }

    sgData = {
        "rpv_stop_350" : info.DataSetInfo(basedir=path, fileName="rpv_stop_350.root", sys=-1.0, label="Signal"),
        "rpv_stop_450" : info.DataSetInfo(basedir=path, fileName="rpv_stop_450.root", sys=-1.0, label="Signal"),
        "rpv_stop_550" : info.DataSetInfo(basedir=path, fileName="rpv_stop_550.root", sys=-1.0, label="Signal"),
        "rpv_stop_650" : info.DataSetInfo(basedir=path, fileName="rpv_stop_650.root", sys=-1.0, label="Signal"),
        "rpv_stop_750" : info.DataSetInfo(basedir=path, fileName="rpv_stop_750.root", sys=-1.0, label="Signal"),
        "rpv_stop_850" : info.DataSetInfo(basedir=path, fileName="rpv_stop_850.root", sys=-1.0, label="Signal"),

        "stealth_stop_350_SYY" : info.DataSetInfo(basedir=path, fileName="stealth_stop_350_SYY.root", sys=-1.0, label="Signal"),
        "stealth_stop_450_SYY" : info.DataSetInfo(basedir=path, fileName="stealth_stop_450_SYY.root", sys=-1.0, label="Signal"),
        "stealth_stop_550_SYY" : info.DataSetInfo(basedir=path, fileName="stealth_stop_550_SYY.root", sys=-1.0, label="Signal"),
        "stealth_stop_650_SYY" : info.DataSetInfo(basedir=path, fileName="stealth_stop_650_SYY.root", sys=-1.0, label="Signal"),
        "stealth_stop_750_SYY" : info.DataSetInfo(basedir=path, fileName="stealth_stop_750_SYY.root", sys=-1.0, label="Signal"),
        "stealth_stop_850_SYY" : info.DataSetInfo(basedir=path, fileName="stealth_stop_850_SYY.root", sys=-1.0, label="Signal"),
    }

    # Loop over all background and get their histograms
    bgHistos = {}
    for key in bgData:
        h = bgData[key].getHisto(histoName)
        h.RebinY(rebinVal)
        bgHistos[key] = {"hist" : h}

    # Loop over all signal models and calculate sig. for each one 
    for key in sgData:
        hSg = sgData[key].getHisto(histoName)
        hSg.RebinY(rebinVal)
        nBins = hSg.GetNbinsY()
        low = hSg.GetYaxis().GetXmin()
        high = hSg.GetYaxis().GetXmax()
        binToGeV = (high - low)/nBins
        sigDic = {}
            
        # Find all possible MLB cut values
        c1 = ROOT.TCanvas( "c", "c", 0, 0, 800, 800)
        ROOT.gStyle.SetStatY(0.5)
        ROOT.gStyle.SetStatX(0.85)
        ROOT.gStyle.SetPalette(ROOT.kRainBow)
        ROOT.gPad.SetLeftMargin(0.12)
        ROOT.gPad.SetRightMargin(0.15)
        ROOT.gPad.SetTopMargin(0.08)
        ROOT.gPad.SetBottomMargin(0.12)
        ROOT.gPad.SetTicks(1,1)
        sig2DHisto = ROOT.TH2D(key+"hist",key+"hist", nBins, low, high, nBins, low, high)        
        for lowBin in range(1, nBins+1):       
            for highBin in range(nBins+1):
                if(lowBin < highBin):
                    sigTot2 = 0.0
                    for nJetBin in range(6, hSg.GetNbinsX()+1):            
                        nSg = hSg.Integral(nJetBin, nJetBin, lowBin, highBin)
                        
                        # Loop over all background histos and calculate sigma
                        sigma2 = 0.0
                        for k in bgHistos:
                            n = bgHistos[k]["hist"].Integral(nJetBin, nJetBin, lowBin, highBin)
                            sys = ((bgData[k].sys)*n)**2
                            sigma2 += n + sys
            
                        # Calculate and fill sig.
                        sig = 0.0
                        if(sigma2 != 0): sig = nSg/math.sqrt(sigma2)
                        sigTot2 += sig**2

                    sig2DHisto.SetBinContent(lowBin, highBin, math.sqrt(sigTot2))

        # Make the plot nice
        sig2DHisto.SetTitle(key+" Significance")
        sig2DHisto.GetXaxis().SetTitle("Lower Cut [GeV]")
        sig2DHisto.GetYaxis().SetTitle("Higher Cut [GeV]")
        sig2DHisto.Draw("colz")
        c1.SaveAs(key+"_sig.png")
        del c1                

if __name__ == '__main__':
    main()
