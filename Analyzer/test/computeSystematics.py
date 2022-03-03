import ROOT, plot, math, sys, argparse, os

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2()

# This class accesses information from fit results and njets for Aron in order to compute
# single and double ratio per-MVA bin shapes
class SystematicComputer:

    # Things to be passed to the constructor are:
    # rawName: string used for extracting raw njets histos
    # fitName: usually the same as rawName, though overridden when applicable, used for extracting from fit
    # color: color for plotting purposes
    # special: naming structure is a little different for 2016 isr/fsr
    # symmetrize: should double ratios be symmetrized
    # njetsPath: path to njets_for_Aron
    # fitPath: path to fit results for systematic
    def __init__(self, process, rawName, fitName, color, special, symmetrize, njetsPath, fitPath):

        self.process = process

        # The special case of nom is when the member name should be blank
        # "nom" is still used for extracting fit results
        self.rawName = rawName 
        if rawName != "": self.rawName = "_" + self.rawName

        # If nothing is passed for fitname, assume it is same as raw name
        self.fitName = fitName
        if self.fitName == "": self.fitName = self.rawName[1:]

        self.fitPathShared = fitPath.replace("_shared", "_%s_shared"%(self.fitName))
        self.fitPathSep    = fitPath.replace("_shared", "_%s_sep"%(self.fitName))

        self.useFits = fitPath != ""

        self.color = color                    # Color to be used for plotting purposes
        self.base = "h_njets_pt30_1l"         # Piece of histogram name in njets for Aron
        self.mvas = ["D1", "D2", "D3", "D4"]  # MVA bins
        self.special = special                # When naming convention of raw histos is flipped
        self.symmetrize = symmetrize          # Symmetrize the double ratio is set to true
        self.symmDoubleRatioHistosFit = {}    # Map MVA bin to symmetrized double ratio using fit results
        self.symmDoubleRatioHistosRaw = {}    # Map MVA bin to symmetrized double ratio using raw results
        self.njetHistosRaw = {}               # Map MVA bin to raw njets histos
        self.njetSumHistoRaw = 0              # Sum histogram of njets in all four MVA bins
        self.singleRatioHistosRaw = {}        # Map MVA bin to single ratio of raw njets
        self.doubleRatioHistosRaw = {}        # Map MVA bin to double ratio of raw njets

        self.njetHistosFitSep = {}            # Map MVA bin to separate fit of njets 
        self.njetHistosFitShared = {}         # Map MVA bin to shared fit of njets 
        self.njetSumHistoFit = 0              # Sum histogram of njets in all four MVA bins
        self.singleRatioHistosFit = {}        # Map MVA bin to single ratio of fit njets
        self.doubleRatioHistosFit = {}        # Map MVA bin to double ratio of fit njets
        self.tripleRatioHistosFit = {}        # Map MVA bin to triple ratio of fit njets

        # Path to njets_for_Aron file for getting raw systematically varied njets shapes
        fnjets = ROOT.TFile.Open(njetsPath, "READ")

        # Get each of the four raw njets histograms
        for mva in self.mvas: 

            # Special is when the njets histogram name is shuffled
            histoName = ""
            if self.special: histoName = mva + "_" + self.process + self.rawName + "_h_njets_pt30_1l"
            else:            histoName = mva + "_" + self.process + "_" + self.base + self.rawName

            self.njetHistosRaw[mva] = fnjets.Get(histoName); self.njetHistosRaw[mva].SetDirectory(0)

            # A slightly hacky way of making signal-as-a-systematic---add ttbar shape to signal shape
            if "TT" not in self.process: self.njetHistosRaw[mva].Add(fnjets.Get(histoName.replace(self.process,"TT")))

        # Sum the njets histograms over all MVA bins to make a "total" histogram
        for mva in self.mvas:
            if mva == "D1": self.njetSumHistoRaw = self.njetHistosRaw[mva].Clone("Sum" + self.rawName)
            else:           self.njetSumHistoRaw.Add(self.njetHistosRaw[mva])

        # Normalize all of the histos, the per-MVA and the total njets
        self.njetSumHistoRaw.Scale(1.0 / self.njetSumHistoRaw.Integral())
        for mva in self.mvas: self.njetHistosRaw[mva].Scale(1.0 / self.njetHistosRaw[mva].Integral())
        
        # Compare njets shape in given mva bin to total shape and compute ratio
        for mva in self.mvas: 

            ratio = self.njetHistosRaw[mva].Clone(mva + "_ratio" + self.rawName); ratio.SetDirectory(0)
            ratio.Divide(self.njetSumHistoRaw)
            self.singleRatioHistosRaw[mva] = ratio

        # Now get the equivalent histos as done above and compute ratios using the fitted version of things
        # If we are told do so, i.e. if a path to fits was passed
        if self.useFits:
            fshared = ROOT.TFile.Open(self.fitPathShared, "READ")
            for mva in self.mvas:
                h = fshared.Get("shapes_fit_b/%s/TT"%(mva)); h.SetDirectory(0)
                self.njetHistosFitShared[mva] = h
            fshared.Close()
            
            fsep = ROOT.TFile.Open(self.fitPathSep, "READ")
            for mva in self.mvas:
                h = fsep.Get("shapes_fit_b/%s/TT"%(mva)); h.SetDirectory(0)
                self.njetHistosFitSep[mva] = h
            fsep.Close()

            # Make the single ratio using these fit njets shapes
            for mva in self.mvas:
                ratio = self.njetHistosFitSep[mva].Clone(mva + "_ratio%s_fit"%(self.rawName)); ratio.SetDirectory(0)
                ratio.Divide(self.njetHistosFitShared[mva])
                self.singleRatioHistosFit[mva] = ratio

        fnjets.Close()

    # Return bool designating if double ratios should be symmetrized
    def doSymmetrize(self): return self.symmetrize

    # Return bool designating if fit results are anticipated
    def doFits(self): return self.useFits
        
    # Return the color to be used for plotting
    def getColor(self): return self.color

    # Return the raw name to be used for plotting
    def getRawName(self): return self.rawName

    # Return the fit name to be used for plotting
    def getFitName(self): return self.fitName

    # A simple getter to return the dictionary of single ratio histograms (raw)
    def getSingleRatioHistosRaw(self): return self.singleRatioHistosRaw

    # A simple getter to return the dictionary of double ratio histograms (raw)
    def getDoubleRatioHistosRaw(self): return self.doubleRatioHistosRaw

    # A simple getter to return the dictionary of single ratio histograms (fit)
    def getSingleRatioHistosFit(self): return self.singleRatioHistosFit

    # A simple getter to return the dictionary of double ratio histograms (fit)
    def getDoubleRatioHistosFit(self): return self.doubleRatioHistosFit

    # A simple getter to return a single ratio histogram (raw)
    def getSingleRatioHistoRaw(self, mva): return self.singleRatioHistosRaw[mva]

    # A simple getter to return a double ratio histogram (raw)
    def getDoubleRatioHistoRaw(self, mva): return self.doubleRatioHistosRaw[mva]

    # A simple getter to return a single ratio histogram (fit)
    def getSingleRatioHistoFit(self, mva): return self.singleRatioHistosFit[mva]

    # A simple getter to return a double ratio histogram (fit)
    def getDoubleRatioHistoFit(self, mva): return self.doubleRatioHistosFit[mva]

    # A simple getter to return a double ratio histogram (fit)
    def getTripleRatioHistoFit(self, mva): return self.tripleRatioHistosFit[mva]

    # Pass in a systematic computer object for another systematic and use it to compute the double ratio between the two
    # Here the only real use case is to pass in the sysetmatic computer object for the nominal systematic
    def computeDoubleRatios(self, nomComputer):

        if self.useFits:
            # Compute double ratios for fitted things first
            nomSingleRatiosFit = nomComputer.getSingleRatioHistosFit()
            for mva, mvaSingleRatio in nomSingleRatiosFit.iteritems():
                newName = mva + "_ratio" + self.rawName + "_div_fit"
                self.doubleRatioHistosFit[mva] = self.singleRatioHistosFit[mva].Clone(newName)
                self.doubleRatioHistosFit[mva].Divide(mvaSingleRatio)

        # Can always compute double ratios for raw things
        nomSingleRatiosRaw = nomComputer.getSingleRatioHistosRaw()
        for mva, mvaSingleRatio in nomSingleRatiosRaw.iteritems():
            newName = mva + "_ratio_div" + self.rawName
            self.doubleRatioHistosRaw[mva] = self.singleRatioHistosRaw[mva].Clone(newName)
            self.doubleRatioHistosRaw[mva].Divide(mvaSingleRatio)

    # This is really a special case for one systematic and should only be used for such
    # The use case is for dividing the mpTScalednoHT double ratio by the noHT double ratio
    def computeTripleRatios(self, noHTcomputer):

        # Compute triple ratios for fitted things
        noHTdoubleRatiosFit = noHTcomputer.getDoubleRatioHistosFit()
        for mva, mvaDoubleRatio in noHTdoubleRatiosFit.iteritems():
            newName = mva + "_ratio" + self.rawName + "_div_fit"
            self.tripleRatioHistosFit[mva] = self.doubleRatioHistosFit[mva].Clone(newName)
            self.tripleRatioHistosFit[mva].Divide(mvaDoubleRatio)

    # Symmetrize the double ratio---only for certain systematics, so check symmetrize member variable
    # Use case is to call this on the up variation systematic computer object while passing in down variation
    def symmetrizeDoubleRatios(self, downSystComputer):

        # Start with symmetrizing the raw double ratios
        downDoubleRatiosRaw = downSystComputer.getDoubleRatioHistosRaw()

        njetBins = self.doubleRatioHistosRaw["D1"].GetNbinsX()
        
        for mva, doubleRatioRaw in self.doubleRatioHistosRaw.iteritems():

            newName = mva + self.rawName.replace("Up","")
            symmDoubleRatio = ROOT.TH1D(newName, newName, njetBins, 0, njetBins); symmDoubleRatio.SetDirectory(0)

            for bin in xrange(njetBins):

                upval   = doubleRatioRaw.GetBinContent(bin+1)
                downval = downDoubleRatiosRaw[mva].GetBinContent(bin+1)
                if upval > 1.0 and downval < 1.0:
                    if downval > 0: downval = 1.0 / downval
                    else:           downval = upval
                elif upval < 1.0 and downval > 1.0: 
                    downval = 1 / downval

                average = (upval + downval) / 2.0
                if average > 0: symmDoubleRatio.SetBinContent(bin+1, average)
                else:           symmDoubleRatio.SetBinContent(bin+1, symmDoubleRatio.GetBinContent(bin))

            self.symmDoubleRatioHistosRaw[mva] = symmDoubleRatio 
        
        # If fitting shapes to get ratios and double ratios, go through symmetrization procedure for them as well
        if self.useFits:
            # For using the fit info, the symmetrization is using maximum between Up/Down variations, but preserving "sign"
            downDoubleRatiosFit = downSystComputer.getDoubleRatioHistosFit()

            njetBins = self.doubleRatioHistosFit["D1"].GetNbinsX()
            
            for mva, doubleRatioFit in self.doubleRatioHistosFit.iteritems():

                newName = mva + self.fitName.replace("Up","")
                symmDoubleRatio = ROOT.TH1D(newName, newName, njetBins, 0, njetBins); symmDoubleRatio.SetDirectory(0)

                for bin in xrange(njetBins):

                    upval = doubleRatioFit.GetBinContent(bin+1); newupval = upval
                    if upval < 1.0: newupval = 1.0 / upval

                    downval = downDoubleRatiosFit[mva].GetBinContent(bin+1); newdownval = downval
                    if downval < 1.0: newdownval = 1.0 / downval

                    if newdownval > newupval:
                        if doubleRatioFit.GetBinContent(bin+1) > 1: symmDoubleRatio.SetBinContent(bin+1, newdownval)
                        else:                                       symmDoubleRatio.SetBinContent(bin+1, 1.0 / newdownval)
                    else:
                        symmDoubleRatio.SetBinContent(bin+1, upval)

                self.symmDoubleRatioHistosFit[mva] = symmDoubleRatio 

            # Finally, for 2016 FSR, adjust ratio closer to "reduced" value per Top systematics recommendations 
            # Historical note: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopSystematics#Parton_shower_ISR_and_FSR_govern 
            if self.special and "FSR" in self.fitName:
                for mva, symmDoubleRatioFit in self.symmDoubleRatioHistosFit.iteritems():

                    for bin in xrange(symmDoubleRatioFit.GetNbinsX()):
                        oldval = symmDoubleRatioFit.GetBinContent(bin+1)
                        diffval = abs(oldval-1) / math.sqrt(2.0)
                        newval = 1 + diffval if (oldval > 1.0) else 1 - diffval
                        symmDoubleRatioFit.SetBinContent(bin+1, newval)

    # Write relevant double ratio histograms to outputfile that is passed in here
    # For sysetmatics that have been designated to be symmetrized, write those double ratios
    # Otherwise, write the double ratios using the fit information
    def writeRatioHistos(self, outputFile):

        outputFile.cd()
        if self.symmetrize: 

            if self.useFits:
                outName = self.fitName.replace("Down", "").replace("Up", "")
                for mva, symmDoubleRatioHisto in self.symmDoubleRatioHistosFit.iteritems(): symmDoubleRatioHisto.Write(mva + "_" + outName)
            else:
                outName = self.rawName.replace("Down", "").replace("Up", "")
                for mva, symmDoubleRatioHisto in self.symmDoubleRatioHistosRaw.iteritems(): symmDoubleRatioHisto.Write(mva + outName)
   
        else:
            for mva, doubleRatioHisto in self.doubleRatioHistosFit.iteritems(): doubleRatioHisto.Write(mva + "_" + self.fitName)        

if __name__ == '__main__':

    usage = "usage: %prog [options]"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("--year",      dest="year",      help="Which period of Run2",       default="NULL", required=True, type=str) 
    parser.add_argument("--fitdir",    dest="fitdir",    help="Directory containing fits",  default="NULL", required=True, type=str) 
    parser.add_argument("--outputdir", dest="outputdir", help="Output directory",           default="NULL", required=True, type=str) 
    parser.add_argument("--fittag",    dest="fittag",    help="Unique tag for fit results", default="NULL", required=True, type=str) 
    arg = parser.parse_args()
    
    year = arg.year
    
    isrFsrSpecial = "2016" in year
    fitdir = arg.fitdir
    fittag = arg.fittag + "_" + year

    outPath = "./%s/%s/"%(arg.outputdir, year)
    if not os.path.exists(outPath): os.makedirs(outPath)

    USER = os.getenv("USER")
    
    njetsPath = "/uscms_data/d3/%s/susy/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/Keras_%s_v1.2/njets_for_Aron.root"%(USER, year)
    
    # When passed to SystematicComputer class, this path will be modified to include systematic name
    fitPathBasic = "%s/%s_shared/fitDiagnostics%sRPV550.root"%(fitdir, fittag, year)
    
    systComputers = {}
    systComputers["nom"]               = SystematicComputer("TT",      "",              "nom",     ROOT.kBlack,                     False, False,         njetsPath, fitPathBasic)

    systComputers["noHT"]              = SystematicComputer("TT",      "noHT",          "",        ROOT.kGreen-6,                   False, False,         njetsPath, fitPathBasic)
    systComputers["mpTScalednoHT"]     = SystematicComputer("TT",      "mpTScalednoHT", "",        ROOT.kMagenta+3,                 False, False,         njetsPath, fitPathBasic)

    systComputers["erdOn"]             = SystematicComputer("TT",      "erdOn",             "",    ROOT.kBlue,                      True, False,         njetsPath, fitPathBasic)
    systComputers["hdampDown"]         = SystematicComputer("TT",      "hdampDown",         "",    ROOT.kRed,                       True, False,         njetsPath, fitPathBasic)
    systComputers["hdampUp"]           = SystematicComputer("TT",      "hdampUp",           "",    ROOT.kGreen+2,                   True, False,         njetsPath, fitPathBasic)
    systComputers["underlyingEvtDown"] = SystematicComputer("TT",      "underlyingEvtDown", "",    ROOT.kViolet,                    True, False,         njetsPath, fitPathBasic)
    systComputers["underlyingEvtUp"]   = SystematicComputer("TT",      "underlyingEvtUp",   "",    ROOT.kOrange+7,                  True, False,         njetsPath, fitPathBasic)

    # Example how to make signal-as-a-systematic
    systComputers["RPV_350"]           = SystematicComputer("RPV_350", "",              "",        ROOT.TColor.GetColor("#D8D8D8"), False, False,         njetsPath, "")
    systComputers["RPV_550"]           = SystematicComputer("RPV_550", "",              "",        ROOT.TColor.GetColor("#969696"), False, False,         njetsPath, "")
    systComputers["RPV_850"]           = SystematicComputer("RPV_850", "",              "",        ROOT.kBlack,                     False, False,         njetsPath, "")

    systComputers["btgUp"]            = SystematicComputer("TT",     "btgUp",         "",        ROOT.kRed,                       False, True,         njetsPath, "")
    systComputers["btgDown"]          = SystematicComputer("TT",     "btgDown",       "",        ROOT.kRed+2,                     False, True,         njetsPath, "")
    systComputers["lepUp"]            = SystematicComputer("TT",     "lepUp",         "",        ROOT.kGreen-6,                   False, True,         njetsPath, "")
    systComputers["lepDown"]          = SystematicComputer("TT",     "lepDown",       "",        ROOT.kGreen+2,                   False, True,         njetsPath, "")
    systComputers["pdfUp"]            = SystematicComputer("TT",     "pdfUp",         "",        ROOT.kBlue-7,                    False, True,         njetsPath, "")
    systComputers["pdfDown"]          = SystematicComputer("TT",     "pdfDown",       "",        ROOT.kBlue-5,                    False, True,         njetsPath, "")
    systComputers["sclUp"]            = SystematicComputer("TT",     "sclUp",         "",        ROOT.kGreen-6,                   False, True,         njetsPath, "")
    systComputers["sclDown"]          = SystematicComputer("TT",     "sclDown",       "",        ROOT.kGreen+2,                   False, True,         njetsPath, "")
    systComputers["htUp"]             = SystematicComputer("TT",     "htUp",          "",        ROOT.kBlue-7,                    False, True,         njetsPath, "")
    systComputers["htDown"]           = SystematicComputer("TT",     "htDown",        "",        ROOT.kBlue-5,                    False, True,         njetsPath, "")
    systComputers["puUp"]             = SystematicComputer("TT",     "puUp",          "",        ROOT.kRed,                       False, True,         njetsPath, "")
    systComputers["puDown"]           = SystematicComputer("TT",     "puDown",        "",        ROOT.kRed+2,                     False, True,         njetsPath, "")

    systComputers["JECUp"]            = SystematicComputer("TT",     "JECUp",         "",        ROOT.kCyan+1,                    False, False,         njetsPath, fitPathBasic)
    systComputers["JECDown"]          = SystematicComputer("TT",     "JECDown",       "",        ROOT.kCyan+2,                    False, False,         njetsPath, fitPathBasic)
    systComputers["JERUp"]            = SystematicComputer("TT",     "JERUp",         "",        ROOT.kBlue-7,                    False, False,         njetsPath, fitPathBasic)
    systComputers["JERDown"]          = SystematicComputer("TT",     "JERDown",       "",        ROOT.kBlue-5,                    False, False,         njetsPath, fitPathBasic)
    
    if year == "2016":
        systComputers["isrUp"]            = SystematicComputer("TT",     "isrUp",        "ISRUp",   ROOT.kRed,                       isrFsrSpecial, True,  njetsPath, fitPathBasic)
        systComputers["isrDown"]          = SystematicComputer("TT",     "isrDown",      "ISRDown", ROOT.kRed+2,                     isrFsrSpecial, True,  njetsPath, fitPathBasic)
        systComputers["fsrUp"]            = SystematicComputer("TT",     "fsrUp",        "FSRUp",   ROOT.kCyan+1,                    isrFsrSpecial, True,  njetsPath, fitPathBasic)
        systComputers["fsrDown"]          = SystematicComputer("TT",     "fsrDown",      "FSRDown", ROOT.kCyan+2,                    isrFsrSpecial, True,  njetsPath, fitPathBasic)
    else:
        systComputers["isr2Up"]           = SystematicComputer("TT",     "isr2Up",        "ISRUp",   ROOT.kRed,                      False,         True,  njetsPath, fitPathBasic)
        systComputers["isr2Down"]         = SystematicComputer("TT",     "isr2Down",      "ISRDown", ROOT.kRed+2,                    False,         True,  njetsPath, fitPathBasic)
        systComputers["fsr2Up"]           = SystematicComputer("TT",     "fsr2Up",        "FSRUp",   ROOT.kCyan+1,                   False,         True,  njetsPath, fitPathBasic)
        systComputers["fsr2Down"]         = SystematicComputer("TT",     "fsr2Down",      "FSRDown", ROOT.kCyan+2,                   False,         True,  njetsPath, fitPathBasic)
   
    outputFile = ROOT.TFile.Open("%s/ttbar_systematics_%s.root"%(outPath, year), "RECREATE")

    # First loop through and calculate all the double ratios
    for syst in systComputers.keys():

        # Skip nom systematic, will be used by others
        if syst == "nom": continue
    
        # For any other systematic compute the double ratios using nom
        systComputers[syst].computeDoubleRatios(systComputers["nom"])

    for syst in systComputers.keys():

        # Special case of computing triple ratio
        if syst == "mpTScalednoHT": systComputers[syst].computeTripleRatios(systComputers["noHT"])

        useFits = systComputers[syst].doFits()
        doSymm = systComputers[syst].doSymmetrize()

        for mva in ["D1", "D2", "D3", "D4"]:
    
            # Again, skip nom systematic
            if syst == "nom": continue
    
            # For Up/Down systematics skip the Down one as the Up will be used to do everything
            if "Down" in syst: continue

            # For systematic that has Up/Down _and_ should be symmetrized, do that now
            if "Up" in syst and doSymm: systComputers[syst].symmetrizeDoubleRatios(systComputers[syst.replace("Up", "Down")])

            # Init vectors for histos, their corresponding colors and names
            rawSingleRatios = []; rawDoubleRatios = []
            fitSingleRatios = []; fitDoubleRatios = []
            fitRawDouble = []; fitRawNames = []; fitRawColors = [ROOT.kBlack]
            singleNames = []; doubleNames = []; colors = []
            
            # Get the histo colors to use
            colors.append(systComputers["nom"].getColor())
            colors.append(systComputers[syst].getColor()); fitRawColors.append(systComputers[syst].getColor())
            
            singleNames.append("Nominal"); fitRawNames.append(syst + " (Fit)")
            singleNames.append(syst);      fitRawNames.append(syst + " (Raw)")
    
            doubleNames.append(syst)
            
            # Get the raw single ratios
            rawSingleRatios.append(systComputers["nom"].getSingleRatioHistoRaw(mva))
            rawSingleRatios.append(systComputers[syst].getSingleRatioHistoRaw(mva))
    
            # Get the raw double ratios
            rawDoubleRatios.append(systComputers[syst].getDoubleRatioHistoRaw(mva))

            if useFits:
                # Get the fit single ratios
                fitSingleRatios.append(systComputers["nom"].getSingleRatioHistoFit(mva))
                fitSingleRatios.append(systComputers[syst].getSingleRatioHistoFit(mva))
    
                # Get the fit double ratios
                fitDoubleRatios.append(systComputers[syst].getDoubleRatioHistoFit(mva))

                fitRawDouble.append(systComputers[syst].getDoubleRatioHistoFit(mva))
                fitRawDouble.append(systComputers[syst].getDoubleRatioHistoRaw(mva))
   
            # If this systematic is a "Down" variation, then there is an Up, so bring the corresponding information for Up for plotting
            if "Up" in syst: 

                systDown = syst.replace("Up", "Down")
    
                colors.append(systComputers[systDown].getColor())
    
                singleNames.append(systDown)
                doubleNames.append(systDown)
    
                rawSingleRatios.append(systComputers[systDown].getSingleRatioHistoRaw(mva))
                rawDoubleRatios.append(systComputers[systDown].getDoubleRatioHistoRaw(mva))
    
                if useFits:
                    fitSingleRatios.append(systComputers[systDown].getSingleRatioHistoFit(mva))
                    fitDoubleRatios.append(systComputers[systDown].getDoubleRatioHistoFit(mva))

            systFitClean = systComputers[syst].getFitName().replace("Up","")
            systRawClean = systComputers[syst].getRawName().replace("Up","")
   
            # Now start making relevant plots of everything
            # Start with plotting the raw single ratios for Nominal and then syst(Up/Down)
            plot.makeplot(rawSingleRatios, singleNames, "N_{j}-7", mva+systRawClean+"_single_ratio_raw", plotdir=outPath, linear=True, legendColumns=1, 
                         ymin=0.8, ymax=1.3, ylabel="",
                         colors=colors, norm=False, drawstyle="hist")
            
            # For plotting raw double ratios, exclude nominal as it will just be a flat line at 1. 
            plot.makeplot(rawDoubleRatios, doubleNames, "N_{j}-7", mva+systRawClean+"_double_ratio_raw", plotdir=outPath, linear=True, legendColumns=1, 
                         ymin=0.8, ymax=1.3, ylabel="", dropzeroes=False,
                         colors=colors[1:], norm=False, drawstyle="hist")
            
            # Plot single and double ratios from fit shapes
            if useFits: 
                plot.makeplot(fitSingleRatios, singleNames, "N_{j}-7", mva+"_"+systFitClean+"_single_ratio_fit", plotdir=outPath, linear=True, legendColumns=1, 
                             ymin=0.8, ymax=1.3, ylabel="",
                             colors=colors, norm=False, drawstyle="hist")
         
                plot.makeplot(fitDoubleRatios, doubleNames, "N_{j}-7", mva+"_"+systFitClean+"_double_ratio_fit", plotdir=outPath, linear=True, legendColumns=1, 
                             ymin=0.8, ymax=1.3, ylabel="", dropzeroes=False,
                             colors=colors[1:], norm=False, drawstyle="hist")
            
                # Make a special plot of putting the double ratio from raw shapes and fit shapes together to check how the fit is doing 
                plot.makeplot(fitRawDouble, fitRawNames, "N_{j}-7", mva+"_"+systFitClean+"_double_ratio_fit_raw", plotdir=outPath, linear=True, legendColumns=1, 
                             ymin=0.8, ymax=1.8, ylabel="", dropzeroes=False,
                             colors=fitRawColors, norm=False, drawstyle="e hist")

            # Example of plotting systematic with signal-as-systematic
            if syst == "mpTScalednoHT":
                mpTScalednoHT = systComputers["mpTScalednoHT"].getDoubleRatioHistoFit(mva)
                noHT = systComputers["noHT"].getDoubleRatioHistoFit(mva)
                mpTScaled = systComputers["mpTScalednoHT"].getTripleRatioHistoFit(mva)

                RPV350 = systComputers["RPV_350"].getDoubleRatioHistoRaw(mva); RPV350color = systComputers["RPV_350"].getColor()
                RPV550 = systComputers["RPV_550"].getDoubleRatioHistoRaw(mva); RPV550color = systComputers["RPV_550"].getColor()
                RPV850 = systComputers["RPV_850"].getDoubleRatioHistoRaw(mva); RPV850color = systComputers["RPV_850"].getColor()

                plot.makeplot([mpTScalednoHT, noHT, mpTScaled], ["mpTScalednoHT", "noHT", "mpTScaled"], "N_{j}-7", mva+"_mpTScaled_overlay_div", plotdir="./", linear=True, legendColumns=1, 
                             ymin=0.8, ymax=1.3, ylabel="", dropzeroes=False,
                             colors=[ROOT.kViolet, ROOT.kRed, ROOT.kBlue], norm=False, drawstyle="hist")

                plot.makeplot([mpTScaled, RPV350, RPV550, RPV850], ["mpTScaled", "RPV 350", "RPV 550", "RPV 850"], "N_{j}-7", mva+"_mpTScaled_signal_div", plotdir="./", linear=True, legendColumns=1, 
                             ymin=0.8, ymax=1.3, ylabel="", dropzeroes=False,
                             colors=[ROOT.kBlue, RPV350color, RPV550color, RPV850color], norm=False, drawstyle="hist")

        systComputers[syst].writeRatioHistos(outputFile) 

    outputFile.Close()
