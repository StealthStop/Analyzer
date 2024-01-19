# The HiggsCombineInputs class 
# it makes the inputs (root files, txt files, etc.) 
# to the Higgs combine
# we can add both sys and nuisance parameters to this class
import ROOT

class HiggsCombineInputs:

    # --------------------
    # Initialize the class 
    # --------------------
    def __init__(self, year, njets, signal, channel, samples, regions, edges, outpath):

        self.year     = year
        self.njets    = njets
        self.signal   = signal
        self.channel  = channel
        self.samples  = samples 
        self.regions  = regions
        self.edges    = edges
        self.outpath  = outpath 

        # ----------------
        # open a root file
        # ----------------
        fileName = outpath + "/" + year + "_" + "TT_TTvar_Syst" + "_" +  signal + "_" + channel + edges + ".root" 
        self.f   = ROOT.TFile.Open(fileName, "RECREATE")        

    # -------------------------------------------------------
    # make a root file to include:
    #   -- MC correction values for TT in signal region 
    #   -- MC correction values for TTvar in signal region
    #   -- MC correction factor ratio: TT/TTvar (in MC level)
    #   -- average value of MC corrected data 
    # -------------------------------------------------------
    def make_HiggsCombineInputs_RootFiles(self, TH1, varName): 

        self.f.WriteObject(TH1, "%s"%(varName))

    # ------------------------------------------------------
    # put the MC correction factors and ratio to a root file 
    # ------------------------------------------------------
    def put_MCcorrFactors_toRootFiles(self, dictionary):

        # loop over the dictionary
        for var, value in dictionary.items():
            
            ttVars = {}

            # loop over the sub-dictionary to get ttVar
            for key, subDictionary in value.items():
            
                ttVars = subDictionary.keys()
                break

            # loop over ttVar
            for ttVar in ttVars:

                hist = ROOT.TH1F(self.year + "_" + var + "_" + ttVar, var + "_" + ttVar, len(self.njets), 0, len(self.njets))

                # loop over njets
                for njet in self.njets:                        

                    Njet = int(njet.replace("incl", ""))

                    binIndex = None
                    if self.channel == "0l":
                        binIndex = Njet - len(self.njets) - 2
                    elif self.channel == "1l":
                        binIndex = Njet - len(self.njets) - 1
                    elif self.channel == "2l":
                        binIndex = Njet - len(self.njets) 

                    if var == "nEventsA":
                        hist.SetBinContent(binIndex, 1.0)

                        content = round(value[njet][ttVar][0]) 
                        if content == 0.0 and value[njet][ttVar][0] != 0.0:
                            content = 1.0
                        hist.SetBinError(binIndex,   1.0 / content**0.5)
                    elif "nEventsApred" in var:
                        hist.SetBinContent(binIndex, 1.0)

                        content = value[njet][ttVar][0]
                        error   = value[njet][ttVar][1]
                        hist.SetBinError(binIndex,   error / content)
                    else:
                        hist.SetBinContent(binIndex, value[njet][ttVar][0])
                        hist.SetBinError(binIndex,   value[njet][ttVar][1])

                # put the variables to root file
                self.make_HiggsCombineInputs_RootFiles(hist, (self.year + "_" + var + "_" + ttVar))
                        
    # ----------------------------------------------------------------
    # put the average value of MC corrected data to the same root file 
    #   -- in between for all variances and all val regions
    # ----------------------------------------------------------------
    def put_averageCorrectedDataValue_toRootFiles(self, dictionary, region):

        hist = ROOT.TH1F("average_MCcorrectedData_Syst", "average_MCcorrectedData_Syst", len(self.njets), 0, len(self.njets))

        # loop over the njets
        for njet in self.njets:

            Njet = int(njet.replace("incl", ""))

            if self.channel == "0l":
                hist.SetBinContent((Njet - len(self.njets) - 2), dictionary[njet]) 

            if self.channel == "1l":
                hist.SetBinContent((Njet - len(self.njets) - 1), dictionary[njet])

            if self.channel == "2l":
                hist.SetBinContent((Njet - len(self.njets)), dictionary[njet])

        # put the average value of MC corrected data to the dictionary
        self.make_HiggsCombineInputs_RootFiles(hist, ("%s_average_MCcorrectedData_Syst_%s"%(self.year,region)))


    # ----------------------------------------------------------------
    # put the maximum value of MC corrected data to the same root file 
    #   -- in between for all variances and all val regions
    #   -- this value is used to calculate another type sys.
    # ----------------------------------------------------------------
    def put_maximumCorrectedDataValue_toRootFiles(self, dictionary, region):

        hist = ROOT.TH1F("maximum_MCcorrectedData_Syst", "maximum_MCcorrectedData_Syst", len(self.njets), 0, len(self.njets))

        # loop over the njets
        for njet in self.njets:

            Njet = int(njet.replace("incl", ""))

            if self.channel == "0l":
                hist.SetBinContent((Njet - len(self.njets) - 2), dictionary[njet][0])
                hist.SetBinError((Njet - len(self.njets) - 2), dictionary[njet][1])

            if self.channel == "1l":
                hist.SetBinContent((Njet - len(self.njets) - 1), dictionary[njet][0])
                hist.SetBinError((Njet - len(self.njets) - 1), dictionary[njet][1])

            if self.channel == "2l":
                hist.SetBinContent((Njet - len(self.njets)), dictionary[njet][0])
                hist.SetBinError((Njet - len(self.njets)), dictionary[njet][1])

        # put the mximum value of MC corrected data to the dicionary
        self.make_HiggsCombineInputs_RootFiles(hist, ("%s_maximum_MCcorrectedData_Syst_%s"%(self.year,region)))


    # -------------------
    # close the root file
    # -------------------
    def close_HiggsCombineInputs_RootFiles(self):

        self.f.Close()

