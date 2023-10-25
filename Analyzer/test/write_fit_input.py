import ROOT
ROOT.gStyle.SetOptStat(0)
from math import sqrt
import random
import DataSetInfo as info
import optparse
import os

class WriteNJetPlots:
    def __init__(self):
        self.colors = [ROOT.kCyan+1, ROOT.kMagenta+1, ROOT.kYellow+1,
                       ROOT.kRed+1, ROOT.kGreen+1, ROOT.kBlue+1,
                       ROOT.kBlack]
        self.signalcolors = [ROOT.kViolet-6, ROOT.kTeal-6, ROOT.kPink-6, ROOT.kAzure-6]

    def calcUnc(self,a,unc_a,b,unc_b):
        return  a/b * sqrt( (unc_a/a)**2+(unc_b/b)**2)

    def ratioFunc(self,x,n,m,a0,a1,a2):
        num=(a1-a2)**(x-m)
        den=(a0-a2)**(x-n)
        val=(a2 + (abs(num/den))**(1/(n-m)))
        return val

    def readHistos(self, data, basenameIn, basenameOut, bin, sys):
        histos = []
        for key, dsi in data.iteritems():
            h = None 
            try:
                h = dsi.getHisto(bin+"_pseudodata"+sys+"_"+basenameOut)
            except:
                print "    Didn't find that histo either, the code will now fail"
                pass
            h.Write()
            histos.append(h)
        return histos

    def writeHistos(self, data, basenameIn, basenameOut, bin, sys):
        histos = []
        for key, dsi in data.iteritems():
            h = None
            try:
                h = dsi.getHisto(basenameIn+"_"+bin+sys)
            except:
                try:
                    print "    Using \""+basenameIn+"_"+bin+"\" instead"
                    h = dsi.getHisto(basenameIn+"_"+bin)
                except:
                    print "    Didn't find that histo either, the code will now fail"
                pass
            h.SetName(bin+"_"+dsi.label+"_"+basenameOut+sys)
            h.SetTitle(bin+"_"+dsi.label+"_"+basenameOut+sys)
            h.Write()
            histos.append(h)
        return histos

    def writeHistosSetBins(self, data, suffix, basenameIn, basenameOut, bin, sys, binDic):
        for key, dsi in data.iteritems():
            h = dsi.getHisto(basenameIn+"_"+bin+sys)
            h.SetName(bin+"_"+dsi.label+suffix+"_"+basenameOut+sys)
            h.SetTitle(bin+"_"+dsi.label+suffix+"_"+basenameOut+sys)
            for keyBin, valueList in binDic.iteritems():
                if keyBin == bin:
                    i=0
                    for b in valueList:
                        i+=1
                        if b >= 0.0:
                            h.SetBinContent(i, b)
            h.Write()

    def writeStatHistos(self, data, basenameIn, basenameOut, bin, sys):
        for key, dsi in data.iteritems():
            h = dsi.getHisto(basenameIn+"_"+bin+sys)
            for i in range(1, h.GetNbinsX()+1):
                name = bin+"_"+dsi.label+"_mcStatBin"+str(i)
                hStatUp = h.Clone(name+"Up")
                hStatUp.SetTitle(name+"Up")
                hStatDown = h.Clone(name+"Down")
                hStatDown.SetTitle(name+"Down")
                content = h.GetBinContent(i)
                error = h.GetBinError(i)
                up = content+error if content+error > 0.0 else 0.0
                down = content-error if content-error > 0.0 else 0.0
                hStatUp.SetBinContent(i, up)
                hStatDown.SetBinContent(i, down)
                hStatUp.Write()
                hStatDown.Write()

    def makePseudoData(self, histos, signalhistos, sgData, basename, bin, sys): #make some pseudo_data
        name = bin+"_pseudodata"+sys+"_"+basename
        mynewh = ROOT.TH1D(name, name, histos[0].GetNbinsX(), histos[0].GetBinLowEdge(1), histos[0].GetBinLowEdge(1)+histos[0].GetNbinsX())

        #make some pseudo_data with signal
        pseudodataS_histos = []
        for key, dsi in sgData.iteritems():
            sig = dsi.label
            name = bin+"_pseudodataS"+sys+"_"+sig+"_"+basename
            pseudodataS_histos.append( ROOT.TH1D(name, name, histos[0].GetNbinsX(), histos[0].GetBinLowEdge(1), histos[0].GetBinLowEdge(1)+histos[0].GetNbinsX()) )

        for bin in range(histos[0].GetNbinsX()):
            sumdata = 0
            #sumunc = 0
            for h in histos:
                sumdata += h.GetBinContent(bin+1)
                #sumunc += h.GetBinError(bin+1)
            #sumdata = sumdata + sumunc*2*(random.random()-0.5)
            mynewh.SetBinContent(bin+1, int(round(sumdata)))

            for i, sigh in enumerate(signalhistos):
                sumsig = sigh.GetBinContent(bin+1)
                #sumsigunc = sigh.GetBinError(bin+1)
                sumdatasigh = sumdata + sumsig #+ sumsigunc*2*(random.random()-0.5)
                pseudodataS_histos[i].SetBinContent(bin+1, int(round(sumdatasigh)))

        #mynewh.Write()
        for h in pseudodataS_histos:
            h.Write()

    def makePseudoData_Func(self, histos, suffix, basename, bin, sys, a0, a1, a2):
        name = bin+"_pseudodataFunc"+suffix+sys+"_"+basename
        mynewh = ROOT.TH1D(name, name, histos[0].GetNbinsX(), histos[0].GetBinLowEdge(1), histos[0].GetBinLowEdge(1)+histos[0].GetNbinsX())

        for bin in range(histos[0].GetNbinsX()):
            sumdata = 0
            N7 = 0
            for h in histos:
                if "TT" in h.GetTitle():
                    N7 = h.GetBinContent(1)
                    continue
                sumdata += h.GetBinContent(bin+1)

            nTTbar = N7
            for i in range(1,bin+1):
                nTTbar *= self.ratioFunc(x=float(i-1),n=2.0,m=0.0,a0=a0,a1=a1,a2=a2)
            sumdata += nTTbar
            mynewh.SetBinContent(bin+1, int(round(sumdata)))

        mynewh.Write()

class MakeDataCard:
    def __init__(self, outFile, bgData, otData, sgData, basename, mvaBin):
        self.nMVABins = len(mvaBin)
        self.outFile = outFile
        self.bgData = bgData
        self.otData = otData
        self.sgData = sgData
        self.basename = basename
        self.mvaBin = mvaBin
        bgMVAHistos = self.getMVAHists(bgData, mvaBin, basename)
        otMVAHistos = self.getMVAHists(otData, mvaBin, basename)
        sgMVAHistos = self.getMVAHists(sgData, mvaBin, basename)
        self.dataSets = [(sgMVAHistos, sgData), (bgMVAHistos, bgData), (otMVAHistos, otData)]

    def getMVAHists(self, data, mvaBin, basename):
        histos = []
        for cut in self.mvaBin:
            h = data.getHisto(self.basename+"_"+cut)
            histos.append(h)
        return histos

    def writeDataCard(self, name):
        f = open(name, 'w')

        f.write( "Signal Region 1 Datacard -- signal category" )
        f.write( "\n" )
        f.write( "imax * number of bins\n" )
        f.write( "jmax * number of processes minus 1\n" )
        f.write( "kmax * number of nuisance parameters\n" )
        f.write( "\n" )
        f.write( "-------------------------------------------------------------------------------------------------------------------------------------------\n" )
        f.write( "\n" )
        for i in range(self.nMVABins):
            f.write( "shapes data_obs    sigD{0}    {1} wspace:data_obs_D{0}\n".format(i+1, self.outFile) )
            f.write( "shapes bkg_tt      sigD{0}    {1} wspace:background_tt_D{0}\n".format(i+1, self.outFile) )
            f.write( "shapes bkg_other   sigD{0}    {1} wspace:background_other_D{0}\n".format(i+1, self.outFile) )
            f.write( "shapes signal      sigD{0}    {1} wspace:signal_D{0}\n".format(i+1, self.outFile) )
            f.write( "\n" )
        f.write( "\n" )
        f.write( "-------------------------------------------------------------------------------------------------------------------------------------------\n" )
        bins        = "bin         "
        observation = "observation "
        for i in range(self.nMVABins):
            bins        += " sigD{0} ".format(i+1)
            observation += "  -1   "
        f.write( bins+"\n" )
        f.write( observation+"\n" )
        f.write( "-------------------------------------------------------------------------------------------------------------------------------------------\n" )
        f.write( "# background rate must be taken from _norm param x 1\n" )
        bin      = "bin              "
        process1 = "process          "
        process2 = "process          "
        rate     = "rate             "
        for i in range(self.nMVABins):
            for d in self.dataSets:
                r = 1
                if(d[1].rate):
                    r = d[0][i].Integral()
                bin      += "{0: <16} ".format("sigD"+str(i+1))
                process1 += "{0: <16} ".format(d[1].processName)
                process2 += "{0: <16} ".format(d[1].process)
                rate     += "{0: <16} ".format(r)
        f.write( bin+"\n" )
        f.write( process1+"\n" )
        f.write( process2+"\n" )
        f.write( rate+"\n" )
        f.write( "-------------------------------------------------------------------------------------------------------------------------------------------\n" )
        f.write( "# Normal uncertainties in the signal region\n" )
        lumiSys = "lumi_13TeV      lnN  "
        for i in range(self.nMVABins):
            for d in self.dataSets:
                lumiSys += "{0: <5} ".format(d[1].lumiSys)
        f.write( lumiSys+"\n" )
        f.write( "-------------------------------------------------------------------------------------------------------------------------------------------\n" )
        for i in range(self.nMVABins):
            f.write( "np_tt_D{0} param 0.0 1.0\n".format(i+1) )

if __name__ == "__main__":
    parser = optparse.OptionParser("usage: %prog [options]\n")
    parser.add_option('-d',         dest='basedir',  type='string', default='condor/JECPlots',     help="Path to input rootfiles")
    parser.add_option('--rootFile', dest='rootFile', type='string', default='njets_for_Aron.root', help="Output root file")
    parser.add_option('--dataCard', dest='dataCard', type='string', default='dataCard.txt',        help="Output data card file")
    parser.add_option('-H',         dest='outDir',   type='string', default='',                    help="Can pass in the output directory name")
    parser.add_option('-y',         dest='year',     type='string', default='2016',                help="Can pass in the run year")
    parser.add_option("--sB",       dest='scaleB',   type='float',  default=1.0,                   help="Can scale to a different lumi for background")
    parser.add_option("--sS",       dest='scaleS',   type='float',  default=1.0,                   help="Can scale to a different lumi for signal")
    options, args = parser.parse_args()

    # Where the root files are stored
    basedir = options.basedir + "/"
    outDir = options.outDir
    if os.path.exists(outDir):
        print "Failed: Output directory %s already exits" % ('"'+outDir+'"')
        exit(0)
    else:
        os.makedirs(outDir)
    jettypes = ["pt30_1l"]
    mvaBin = [ "D1", "D2","D3","D4"]
    systypes = ["", "_JECUp", "_JECDown", "_JERUp", "_JERDown", "_btgUp", "_btgDown", "_lepUp", "_lepDown",
                "_isrUp", "_isrDown", "_fsrUp", "_fsrDown", "_isr2Up", "_isr2Down", "_fsr2Up", "_fsr2Down",
                "_pdfUp", "_pdfDown", "_htUp", "_htDown", "_puUp", "_puDown", "_sclUp", "_sclDown", "_prfUp", "_prfDown",
                "_mpTScaled", "_noHT", "_mpTScalednoHT"]
    outputfile = ROOT.TFile.Open(outDir + "/" + options.rootFile,"RECREATE")
    outputDataCard = options.dataCard

    # I hadd my MakeN files into TT.root, TTX, and OTHER, then take QCD from control region
    Data = {
        #"data" : info.DataSetInfo(basedir=basedir, fileName=options.year+"_Data.root", label="data", processName="data", process="1", rate=False, lumiSys="-", scale=-1.0)
    }

    bgData = {
        "pseudodata" : info.DataSetInfo(basedir=basedir, fileName="njets_for_Aron.root",      label="pseudodata",    processName="pseudodata",    process="1", rate=False, lumiSys="-", scale=options.scaleB),
    }

    sgData = {}
    stopMasses = [300, 500, 700, 900, 1100, 1300, 1500]
    ctaus = ["0p01", "0p1", "1", "10", "100", "1000"]
    for stopMass in stopMasses:
        for lspMass in [100, stopMass-225]:
            for ctau in ctaus:

                label = "SYY_%d_%s_%d"%(lspMass, ctau, stopMass)
                sgData[label] = info.DataSetInfo(basedir=basedir, fileName=options.year+"_StealthSYY_2t6j_mStop-%d_mSo-%d_ctau-%s.root"%(stopMass, lspMass, ctau), label=label, processName="signal", process="0", rate=True, lumiSys="1.05", scale=options.scaleS)

                label = "SHH_%d_%s_%d"%(lspMass, ctau, stopMass)
                sgData[label] = info.DataSetInfo(basedir=basedir, fileName=options.year+"_StealthSHH_2t4b_mStop-%d_mSo-%d_ctau-%s.root"%(stopMass, lspMass, ctau), label=label, processName="signal", process="0", rate=True, lumiSys="1.05", scale=options.scaleS)

    if options.year == "2016":
        binDicData = {
            "D1" : [-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0], 
            "D2" : [-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0],
            "D3" : [-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0],
            "D4" : [-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0],
        }

        TTBar_ISR_FSR_SYS_2016 = {
            "TT_isrUp"             : info.DataSetInfo(basedir=basedir, fileName="2016_TT_isrUp.root",             label="TT_isrUp",             processName="bg", process="0", rate=False, lumiSys="-", scale=-1.0),
            "TT_isrDown"           : info.DataSetInfo(basedir=basedir, fileName="2016_TT_isrDown.root",           label="TT_isrDown",           processName="bg", process="0", rate=False, lumiSys="-", scale=-1.0),
            "TT_fsrUp"             : info.DataSetInfo(basedir=basedir, fileName="2016_TT_fsrUp.root",             label="TT_fsrUp",             processName="bg", process="0", rate=False, lumiSys="-", scale=-1.0),
            "TT_fsrDown"           : info.DataSetInfo(basedir=basedir, fileName="2016_TT_fsrDown.root",           label="TT_fsrDown",           processName="bg", process="0", rate=False, lumiSys="-", scale=-1.0),
            "TT_erdOn"             : info.DataSetInfo(basedir=basedir, fileName="2016_TT_erdOn.root",             label="TT_erdOn",             processName="bg", process="0", rate=False, lumiSys="-", scale=-1.0),
            "TT_hdampUp"           : info.DataSetInfo(basedir=basedir, fileName="2016_TT_hdampUp.root",           label="TT_hdampUp",           processName="bg", process="0", rate=False, lumiSys="-", scale=-1.0),
            "TT_hdampDown"         : info.DataSetInfo(basedir=basedir, fileName="2016_TT_hdampDown.root",         label="TT_hdampDown",         processName="bg", process="0", rate=False, lumiSys="-", scale=-1.0),
            "TT_underlyingEvtUp"   : info.DataSetInfo(basedir=basedir, fileName="2016_TT_underlyingEvtUp.root",   label="TT_underlyingEvtUp",   processName="bg", process="0", rate=False, lumiSys="-", scale=-1.0),
            "TT_underlyingEvtDown" : info.DataSetInfo(basedir=basedir, fileName="2016_TT_underlyingEvtDown.root", label="TT_underlyingEvtDown", processName="bg", process="0", rate=False, lumiSys="-", scale=-1.0),
        }

    elif options.year== "2017":
        #sgData.update({
        #    "RPVCP5_350" : info.DataSetInfo(basedir=basedir, fileName=options.year+"_RPV_2t6j_mStop-350_CP5.root", label="RPVCP5_350", processName="signal", process="0", rate=True, lumiSys="1.05", scale=options.scaleS),
        #    "RPVCP5_550" : info.DataSetInfo(basedir=basedir, fileName=options.year+"_RPV_2t6j_mStop-550_CP5.root", label="RPVCP5_550", processName="signal", process="0", rate=True, lumiSys="1.05", scale=options.scaleS),
        #    "RPVCP5_850" : info.DataSetInfo(basedir=basedir, fileName=options.year+"_RPV_2t6j_mStop-850_CP5.root", label="RPVCP5_850", processName="signal", process="0", rate=True, lumiSys="1.05", scale=options.scaleS),
        #    })

        binDicData = {
            "D1" : [-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,2.0], 
            "D2" : [-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,2.0],
            "D3" : [-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,0.0],
            "D4" : [-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,0.0],
        }

        #TTBar_SYS_2017 = {
        #    "TT_erdOn"             : info.DataSetInfo(basedir=basedir, fileName="2017_TT_erdOn.root",             label="TT_erdOn",             processName="bg", process="0", rate=False, lumiSys="-", scale=-1.0),
        #    "TT_hdampUp"           : info.DataSetInfo(basedir=basedir, fileName="2017_TT_hdampUp.root",           label="TT_hdampUp",           processName="bg", process="0", rate=False, lumiSys="-", scale=-1.0),
        #    "TT_hdampDown"         : info.DataSetInfo(basedir=basedir, fileName="2017_TT_hdampDown.root",         label="TT_hdampDown",         processName="bg", process="0", rate=False, lumiSys="-", scale=-1.0),
        #    "TT_underlyingEvtUp"   : info.DataSetInfo(basedir=basedir, fileName="2017_TT_underlyingEvtUp.root",   label="TT_underlyingEvtUp",   processName="bg", process="0", rate=False, lumiSys="-", scale=-1.0),
        #    "TT_underlyingEvtDown" : info.DataSetInfo(basedir=basedir, fileName="2017_TT_underlyingEvtDown.root", label="TT_underlyingEvtDown", processName="bg", process="0", rate=False, lumiSys="-", scale=-1.0),
        #}

    elif options.year== "2018pre":
        binDicData = {
            "D1" : [-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,2.0], 
            "D2" : [-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,2.0],
            "D3" : [-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,0.0],
            "D4" : [-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,0.0],
        }

        TTBar_SYS_2018pre = {
            "TT_erdOn"             : info.DataSetInfo(basedir=basedir, fileName="2018pre_TT_erdOn.root",             label="TT_erdOn",             processName="bg", process="0", rate=False, lumiSys="-", scale=-1.0),
            "TT_hdampUp"           : info.DataSetInfo(basedir=basedir, fileName="2018pre_TT_hdampUp.root",           label="TT_hdampUp",           processName="bg", process="0", rate=False, lumiSys="-", scale=-1.0),
            "TT_hdampDown"         : info.DataSetInfo(basedir=basedir, fileName="2018pre_TT_hdampDown.root",         label="TT_hdampDown",         processName="bg", process="0", rate=False, lumiSys="-", scale=-1.0),
            "TT_underlyingEvtUp"   : info.DataSetInfo(basedir=basedir, fileName="2018pre_TT_underlyingEvtUp.root",   label="TT_underlyingEvtUp",   processName="bg", process="0", rate=False, lumiSys="-", scale=-1.0),
            "TT_underlyingEvtDown" : info.DataSetInfo(basedir=basedir, fileName="2018pre_TT_underlyingEvtDown.root", label="TT_underlyingEvtDown", processName="bg", process="0", rate=False, lumiSys="-", scale=-1.0),
        }

    elif options.year== "2018post":
        binDicData = {
            "D1" : [-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,2.0], 
            "D2" : [-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,2.0],
            "D3" : [-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,0.0],
            "D4" : [-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,0.0],
        }

        TTBar_SYS_2018post = {
            "TT_erdOn"             : info.DataSetInfo(basedir=basedir, fileName="2018post_TT_erdOn.root",             label="TT_erdOn",             processName="bg", process="0", rate=False, lumiSys="-", scale=-1.0),
            "TT_hdampUp"           : info.DataSetInfo(basedir=basedir, fileName="2018post_TT_hdampUp.root",           label="TT_hdampUp",           processName="bg", process="0", rate=False, lumiSys="-", scale=-1.0),
            "TT_hdampDown"         : info.DataSetInfo(basedir=basedir, fileName="2018post_TT_hdampDown.root",         label="TT_hdampDown",         processName="bg", process="0", rate=False, lumiSys="-", scale=-1.0),
            "TT_underlyingEvtUp"   : info.DataSetInfo(basedir=basedir, fileName="2018post_TT_underlyingEvtUp.root",   label="TT_underlyingEvtUp",   processName="bg", process="0", rate=False, lumiSys="-", scale=-1.0),
            "TT_underlyingEvtDown" : info.DataSetInfo(basedir=basedir, fileName="2018post_TT_underlyingEvtDown.root", label="TT_underlyingEvtDown", processName="bg", process="0", rate=False, lumiSys="-", scale=-1.0),
        }

    dicNoD4 = {
        "D1" : [   -1.0,  -1.0,  -1.0, -1.0, -1.0, -1.0, -1.0, -1.0], 
        "D2" : [   -1.0,  -1.0,  -1.0, -1.0, -1.0, -1.0, -1.0, -1.0],
        "D3" : [   -1.0,  -1.0,  -1.0, -1.0, -1.0, -1.0, -1.0, -1.0],
        "D4" : [ 1364.0, 398.0, 109.0, 28.0,  7.0,  2.0,  1.0,  0.0],
    }
    dicNoD3D4 = {
        "D1" : [   -1.0,  -1.0,  -1.0, -1.0, -1.0, -1.0, -1.0, -1.0], 
        "D2" : [   -1.0,  -1.0,  -1.0, -1.0, -1.0, -1.0, -1.0, -1.0],
        "D3" : [ 3250.0, 939.0, 259.0, 63.0, 18.0,  4.0,  1.0,  0.0],
        "D4" : [ 1364.0, 398.0, 109.0, 28.0,  7.0,  2.0,  1.0,  0.0],
    }
    dicNoD1D2D3D4 = {
        "D1" : [29495.0, 8374.0, 2216.0, 563.0, 143.0, 35.0, 9.0, 2.0],
        "D2" : [18308.0, 5295.0, 1454.0, 392.0, 107.0, 24.0, 6.0, 2.0],
        "D3" : [ 3250.0,  939.0,  259.0,  63.0,  18.0,  4.0, 1.0, 0.0],
        "D4" : [ 1364.0,  398.0,  109.0,  28.0,   7.0,  2.0, 1.0, 0.0],
    }    

    #Write all histos to outputfile
    outputfile.cd()
    wp = WriteNJetPlots()
    for jettype in jettypes:
        for bin in mvaBin:
            for sys in systypes:
                basenameIn  = "h_njetsShifted_" + jettype
                basenameOut = "h_njets_" + jettype
                signalhistos = wp.writeHistos(sgData, basenameIn, basenameOut, bin, sys)
                if sys in ["", "_JECUp", "_JECDown", "_JERUp", "_JERDown"]:
                    histos = wp.readHistos(bgData, basenameIn, basenameOut, bin, sys)
                    wp.makePseudoData(histos, signalhistos, sgData, basenameOut, bin, sys)
                    #wp.makePseudoData_Func(histos, "28_24_236", basenameOut, bin, sys, a0=0.28, a1=0.24, a2=0.236)
                    #wp.makePseudoData_Func(histos, "28_24_18",  basenameOut, bin, sys, a0=0.28, a1=0.24, a2=0.18)
                    #wp.makePseudoData_Func(histos, "28_24_-20", basenameOut, bin, sys, a0=0.28, a1=0.24, a2=-0.20)
                    if sys == "":
                    #    wp.writeHistos(Data, basenameIn, basenameOut, bin, sys)
                    #    wp.writeHistosSetBins(Data, "SetBin", basenameIn, basenameOut, bin, sys, binDicData)
                    #    if options.year == "2016":
                    #        wp.writeHistos(TTBar_ISR_FSR_SYS_2016, basenameIn, basenameOut, bin, sys)                            
                    #        wp.writeHistosSetBins(Data, "SetBinNoD4", basenameIn, basenameOut, bin, sys, dicNoD4)
                    #        wp.writeHistosSetBins(Data, "SetBinNoD3D4", basenameIn, basenameOut, bin, sys, dicNoD3D4)
                    #        wp.writeHistosSetBins(Data, "SetBinNoD1D2D3D4", basenameIn, basenameOut, bin, sys, dicNoD1D2D3D4)                            
                    #    if options.year == "2017":
                    #        wp.writeHistos(TTBar_SYS_2017, basenameIn, basenameOut, bin, sys)
                    #    if options.year == "2018pre":
                    #        wp.writeHistos(TTBar_SYS_2018pre, basenameIn, basenameOut, bin, sys)
                    #    if options.year == "2018post":
                    #        wp.writeHistos(TTBar_SYS_2018post, basenameIn, basenameOut, bin, sys)
                    #    wp.writeStatHistos({"OTHER" : bgData["OTHER"]}, basenameIn, basenameOut, bin, sys)
                    #    wp.writeStatHistos({"QCD"   : bgData["QCD"]},   basenameIn, basenameOut, bin, sys)
                    #    wp.writeStatHistos({"TTX"   : bgData["TTX"]}, basenameIn, basenameOut, bin, sys)
                        wp.writeStatHistos(sgData, basenameIn, basenameOut, bin, sys)

    #Close outfile
    outputfile.Close()

    #Old----we no longer put number of events in data card
    ##Make data card
    #for key in sgData:
    #        md = MakeDataCard(outFile="MVA_ws.root", bgData=bgData["TT"], otData=bgData["OTHER"], sgData=sgData[key], basename="h_njets_pt30_1l", mvaBin=mvaBin)
    #        md.writeDataCard(outDir+"/"+key+"_"+outputDataCard)
