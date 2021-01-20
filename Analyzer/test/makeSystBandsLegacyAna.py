# This script reads in histograms from the analysis framework (specifically for systematically varied ttbar samples)
# and computes a "total" systematic band for the full background stack and the data/MC ratio.
# The output files here are intended to be used by plot_1l_LegacyAna.C

import sys, os, ROOT, argparse, math
from collections import defaultdict

ROOT.TH1.SetDefaultSumw2()
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat("")
ROOT.gStyle.SetPaintTextFormat("3.2f")
ROOT.gStyle.SetFrameLineWidth(2)
ROOT.gStyle.SetEndErrorSize(0)

def doOptions(histo, histoName):

    is1D = "TH1" in histo.ClassName()

    for axis, options in OPTIONSMAP[histoName].iteritems():

        if axis == "X":
            if "rebin" in options:
                if is1D: histo.Rebin(options["rebin"])
                else: histo.RebinX(options["rebin"])
            if "min" in options and "max" in options: histo.GetXaxis().SetRangeUser(options["min"],options["max"])
            if "title" in options: histo.GetXaxis().SetTitle(options["title"])
        if axis == "Y":
            if "rebin" in options:
                if is1D: histo.Rebin(options["rebin"])
                else: histo.RebinY(options["rebin"])
            if "min" in options and "max" in options: histo.GetYaxis().SetRangeUser(options["min"],options["max"])
            if "title" in options: histo.GetYaxis().SetTitle(options["title"])
        if axis == "Z":
            if "min" in options and "max" in options: histo.GetZaxis().SetRangeUser(options["min"],options["max"])

def getSystematic(nominal, variation, normalize = False):

    nom = nominal.Clone(nominal.GetName()+"_nominalClone"); var = variation.Clone(variation.GetName()+"_variationClone")

    if normalize: nom.Scale(1.0/nom.Integral()); var.Scale(1.0/var.Integral())

    ratio = var.Clone(var.GetName()+"_ratio"); ratioSymm = var.Clone(var.GetName()+"+ratioSymm")

    ratio.Divide(nom)

    for xbin in xrange(1, ratio.GetNbinsX()+1):
        
        ratioValAbs = abs(1 - ratio.GetBinContent(xbin))
        ratio.SetBinContent(xbin, ratioValAbs)

    return ratio

def getFullSystematic(*argv):

    nom = argv[0].Clone(argv[0].GetName()+"_fullSyst")
    syst = nom.Clone(argv[0].GetName()+"_Syst")

    for xbin in xrange(1, nom.GetNbinsX()+1):

        iVariation = 0
        for histo in argv[1:]: iVariation += histo.GetBinContent(xbin)**2
        
        iVariation = iVariation**0.5

        current = nom.GetBinContent(xbin)

        errBar = current*iVariation

        syst.SetBinContent(xbin, current)
        syst.SetBinError(xbin, errBar)

    return syst

def getDataMCSyst(data, mc, normalize = False):

    dataMcSyst = mc.Clone("mcSyst")

    scale = 1.0
    if normalize: scale = data.Integral()/mc.Integral()
    
    for xbin in xrange(1, mc.GetNbinsX()+1):
   
        datanom = data.GetBinContent(xbin)
        mcnom  = scale*mc.GetBinContent(xbin)
        mcVar  = scale*(mcnom + mc.GetBinError(xbin))

        dataMCerr = 0.0
        if (mcnom == 0 and mcVar != 0) or (mcnom != 0 and mcVar == 0):
            dataMCerr = 999.0
        elif mcnom != 0 and mcVar != 0:
            dataMC = datanom / mcnom
            dataMCvar = datanom / mcVar
            dataMCerr = abs(dataMC - dataMCvar)

        dataMcSyst.SetBinContent(xbin, 1)
        dataMcSyst.SetBinError(xbin, dataMCerr)

    return dataMcSyst

def constructSyst(mc, normUnc = False):
    mcSyst = mc.Clone(mc.GetName() + "mcstat")

    for xbin in xrange(1, mc.GetNbinsX()+1):
    
        count = mc.GetBinContent(xbin)
        
        normErr = 0.0
        if normUnc: normErr = 0.3*count

        err = mc.GetBinError(xbin)

        mcSyst.SetBinContent(xbin, count+err+normErr)

    return mcSyst

def getMCsystBand(normalize, normUnc, dataNom, ttNom, qcdNom, *argv):

    mcSum = 0.0
    mcSum += ttNom.Integral(); mcSum += qcdNom.Integral()
    for mcVar in argv: mcSum += mcVar.Integral()

    systBand         = ttNom.Clone("systBand")
    systBandNoMCstat = ttNom.Clone("systBandNoMCstat")

    scale = 1.0
    if normalize: scale = dataNom.Integral() / mcSum

    for xbin in xrange(1, systBand.GetNbinsX()+1):

        mcErr = 0.0; mc = 0.0; normErr = 0.0; mcErrNoStat = 0.0

        data = dataNom.GetBinContent(xbin)

        tt    = ttNom.GetBinContent(xbin)
        mc    += scale*tt
        ttErr = ttNom.GetBinError(xbin)
        mcErr += ttErr**2.0; mcErrNoStat += ttErr**2.0

        qcd    = qcdNom.GetBinContent(xbin)
        mc     += scale*qcd
        qcdErr = qcdNom.GetBinError(xbin)
        mcErr  += qcdErr**2.0

        # Only TTX and OTH are left
        for mcSyst in argv:
            normErr = 0.0
            if normUnc: normErr = 0.3*mcSyst.GetBinContent(xbin)
            mcErr += (mcSyst.GetBinError(xbin)+normErr)**2.0
            mcErrNoStat += normErr**2.0
            mc += scale*mcSyst.GetBinContent(xbin)
 
        mcErr = math.sqrt(mcErr); mcErrNoStat = math.sqrt(mcErrNoStat)

        systBand.SetBinContent(xbin, mc);  systBandNoMCstat.SetBinContent(xbin, mc)
        systBand.SetBinError(xbin, mcErr); systBandNoMCstat.SetBinError(xbin, mcErrNoStat)

    return systBand, systBandNoMCstat

def fixNJetsHisto(histo):

    nbinsX = histo.GetNbinsX()
    bin12jet = histo.FindBin(12.2)
    for bin in xrange(bin12jet+1,nbinsX+1):
    
        content = histo.GetBinContent(bin)
        histo.Fill(12.2, content)

        histo.SetBinContent(bin, 0.0)
        histo.SetBinError(bin, 0.0)

if __name__ == '__main__':

    usage = "usage: %prog [options]"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("--inputDir", dest="inputDir", help="Path to input"         , default="NULL", type=str) 
    parser.add_argument("--year"    , dest="year"    , help="Which year?"           , default="NULL", type=str) 
    parser.add_argument("--tag"     , dest="tag"     , help="Unique tag for output" , default="NULL", type=str) 
    
    arg = parser.parse_args()
    
    OPTIONSMAP = {"h_deepESM" : {"X" : {"rebin" : 10}},
                  "h_njets" : {"X" : {}},
                  "h_nb" : {"X" : {}},
                  "h_mbl" : {"X" : {"rebin" : 10}},
                  "h_lPt" : {"X" : {"rebin" : 2, "min" : 0, "max" : 1000}},
                  "h_ht" : {"X" : {"rebin" : 10}},
                  "h_fixedGridRhoFastjetAll" : {"X" : {"rebin" : 6}},
                  "fwm2_top6" : {"X" : {"rebin" : 1}},
                  "fwm3_top6" : {"X" : {"rebin" : 1}},
                  "fwm4_top6" : {"X" : {"rebin" : 1}},
                  "jmt_ev0_top6" : {"X" : {"rebin" : 1}},
                  "jmt_ev1_top6" : {"X" : {"rebin" : 1}},
                  "jmt_ev2_top6" : {"X" : {"rebin" : 1}}
    }

    for i in xrange(1, 8):
        OPTIONSMAP["Jet_cm_pt_%d"%(i)]  = {"X" : {"rebin" : 3}}
        OPTIONSMAP["Jet_cm_eta_%d"%(i)] = {"X" : {"rebin" : 1}}
        OPTIONSMAP["Jet_cm_phi_%d"%(i)] = {"X" : {"rebin" : 1}}
        OPTIONSMAP["Jet_cm_m_%d"%(i)]   = {"X" : {"rebin" : 5}}

    if arg.inputDir == "NULL" or arg.year == "NULL": quit()

    year = arg.year
    tag = arg.tag
    INPUTLOC = arg.inputDir + "/%s"%(year)
       
    BACKGROUNDS = ["Data", "TT", "QCD", "TTX", "Other"]
    BACKGROUNDS += ["TT_noHT", "TT_mpTScaled"]
    BACKGROUNDS += ["TT_isrUp", "TT_isrDown", "TT_fsrUp", "TT_fsrDown"]
    BACKGROUNDS += ["TT_erdOn", "TT_hdampUp", "TT_hdampDown", "TT_underlyingEvtUp", "TT_underlyingEvtDown"]
    BACKGROUNDS += ["TT_JECup", "TT_JECdown", "TT_JERup", "TT_JERdown"]

    # Open up all the backgrounds in ROOT
    BKGDFILES = {}; BKGDHISTOS = {}
    for bkgd in BACKGROUNDS: BKGDFILES[bkgd] = ROOT.TFile.Open(INPUTLOC + "/%s_%s.root"%(year, bkgd))

    selections = ["_1l_ge7j_ge1b", "_1l_HT300_ge7j_ge1b_Mbl"]

    # Setup a dictionary of histograms for each process
    for histName in OPTIONSMAP.keys():
        for bkgd, tfile in BKGDFILES.iteritems():
            for selection in selections:

                if "h_" != histName[0:2] and selections.index(selection) > 0 or \
                   "h_" == histName[0:2] and selections.index(selection) == 0: continue

                fullName = histName + selection

                if "noHT" in tfile.GetName():
                    if "h_" == histName[0:2]: fullName += "_noHTWeight"
                    else: fullName += "_NoHTweight"

                if "fsrUp"   in tfile.GetName() and "2016" not in tfile.GetName(): fullName += "_fsrUp"
                if "fsrDown" in tfile.GetName() and "2016" not in tfile.GetName(): fullName += "_fsrDown"
                if "isrUp"   in tfile.GetName() and "2016" not in tfile.GetName(): fullName += "_isrUp"
                if "isrDown" in tfile.GetName() and "2016" not in tfile.GetName(): fullName += "_isrDown"

                histo = tfile.Get(fullName); histo.SetDirectory(0)

                histo.SetName(bkgd + "_" + fullName)

                doOptions(histo, histName)

                histo.SetDirectory(0)

                BKGDHISTOS.setdefault(histName + selection, {}).setdefault(bkgd, histo) 

    f = ROOT.TFile.Open("%s/%s_MC_Syst_%s.root"%(INPUTLOC,year,tag),       "RECREATE")
    g = ROOT.TFile.Open("%s/%s_MC_Ratio_Syst_%s.root"%(INPUTLOC,year,tag), "RECREATE")

    normUnc = False
    if "wNormUnc" in tag: normUnc = True

    for histName, processDict in BKGDHISTOS.iteritems():

        data                = processDict["Data"]
        ttNominal           = processDict["TT"]
        ttnoHT              = processDict["TT_noHT"]
        ttmpTScaled         = processDict["TT_mpTScaled"]

        tterdOn             = processDict["TT_erdOn"]
        tthdampUp           = processDict["TT_hdampUp"]
        tthdampDown         = processDict["TT_hdampDown"]
        ttunderlyingEvtUp   = processDict["TT_underlyingEvtUp"]
        ttunderlyingEvtDown = processDict["TT_underlyingEvtDown"]

        ttisrUp             = processDict["TT_isrUp"]
        ttisrDown           = processDict["TT_isrDown"]
        ttfsrUp             = processDict["TT_fsrUp"]
        ttfsrDown           = processDict["TT_fsrDown"]

        ttJECup             = processDict["TT_JECup"]
        ttJECdown           = processDict["TT_JECdown"]
        ttJERup             = processDict["TT_JERup"]
        ttJERdown           = processDict["TT_JERdown"]

        qcd                 = processDict["QCD"]
        ttx                 = processDict["TTX"]
        other               = processDict["Other"]

        if "njets" in histName:

            fixNJetsHisto(data)
            fixNJetsHisto(ttNominal)
            fixNJetsHisto(ttnoHT)
            fixNJetsHisto(ttmpTScaled)
            fixNJetsHisto(qcd)
            fixNJetsHisto(ttx)
            fixNJetsHisto(other)

            fixNJetsHisto(tterdOn)           
            fixNJetsHisto(tthdampUp)         
            fixNJetsHisto(tthdampDown)       
            fixNJetsHisto(ttunderlyingEvtUp) 
            fixNJetsHisto(ttunderlyingEvtDown)
                               
            fixNJetsHisto(ttisrUp) 
            fixNJetsHisto(ttisrDown)
            fixNJetsHisto(ttfsrUp)
            fixNJetsHisto(ttfsrDown)
                               
            fixNJetsHisto(ttJECup)
            fixNJetsHisto(ttJECdown)
            fixNJetsHisto(ttJERup)
            fixNJetsHisto(ttJERdown)

        mpTScaledSyst         = getSystematic(ttNominal, ttmpTScaled,         normalize=True)
        noHTSyst              = getSystematic(ttNominal, ttnoHT,              normalize=True)

        erdOnSyst             = getSystematic(ttNominal, tterdOn,             normalize=True)
        hdampUpSyst           = getSystematic(ttNominal, tthdampUp,           normalize=True)
        hdampDownSyst         = getSystematic(ttNominal, tthdampDown,         normalize=True)
        underlyingEvtUpSyst   = getSystematic(ttNominal, ttunderlyingEvtUp,   normalize=True)
        underlyingEvtDownSyst = getSystematic(ttNominal, ttunderlyingEvtDown, normalize=True)

        isrUpSyst             = getSystematic(ttNominal, ttisrUp,             normalize=True)
        isrDownSyst           = getSystematic(ttNominal, ttisrDown,           normalize=True)
        fsrUpSyst             = getSystematic(ttNominal, ttfsrUp,             normalize=True)
        fsrDownSyst           = getSystematic(ttNominal, ttfsrDown,           normalize=True)

        JECupSyst             = getSystematic(ttNominal, ttJECup,             normalize=True)
        JECdownSyst           = getSystematic(ttNominal, ttJECdown,           normalize=True)
        JERupSyst             = getSystematic(ttNominal, ttJERup,             normalize=True)
        JERdownSyst           = getSystematic(ttNominal, ttJERdown,           normalize=True)

        ttFullVariation    = getFullSystematic(ttNominal, noHTSyst, mpTScaledSyst, erdOnSyst, hdampUpSyst, hdampDownSyst, underlyingEvtUpSyst, underlyingEvtDownSyst, isrUpSyst, isrDownSyst, fsrUpSyst, fsrDownSyst, JECupSyst, JECdownSyst, JERupSyst, JERdownSyst)

        mcSystBand, mcSystBandNoMCstat = getMCsystBand(True, normUnc, data, ttFullVariation, qcd, ttx, other)

        dataMCsyst         = getDataMCSyst(data, mcSystBandNoMCstat, normalize=False)

        f.cd()
        mcSystBand.SetName(histName)
        mcSystBand.Write(histName)

        g.cd()
        dataMCsyst.SetName(histName)
        dataMCsyst.Write(histName)

    f.Close()
    g.Close()
