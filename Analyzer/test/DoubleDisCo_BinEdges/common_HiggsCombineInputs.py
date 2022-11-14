# The HiggsCombineInputs class 
# it makes the inputs (root files, txt files, etc.) 
# to the Higgs combine
# we can add both sys and nuisance parameters to this class
import ROOT

class HiggsCombineInputs:

    # --------------------
    # Initialize the class 
    # --------------------
    def __init__(self, year, njets, channel, signal, samples, regions, edges):

        self.year     = year
        self.njets    = njets
        self.channel  = channel
        self.signal   = signal
        self.samples  = samples 
        self.regions  = regions
        self.edges    = edges

        # ----------------
        # open a root file
        # ----------------
        fileName = year + "_" + "TT_TTvar_Syst" + "_" + channel + "_" + signal + edges + ".root" 
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

                hist = ROOT.TH1F(var + "_" + ttVar, var + "_" + ttVar, len(self.njets), 0, len(self.njets))

                # loop over njets
                for njet in self.njets:                        

                    Njet = int(njet.replace("incl", ""))

                    hist.SetBinContent((Njet - len(self.njets)), value[njet][ttVar][0])
                    hist.SetBinError((Njet - len(self.njets)), value[njet][ttVar][1])

                # put the variables to root file
                self.make_HiggsCombineInputs_RootFiles(hist, (var + "_" + ttVar))
                        
    # ----------------------------------------------------------------
    # put the average value of MC corrected data to the same root file 
    #   -- in between for all variances and all val regions
    # ----------------------------------------------------------------
    def put_averageCorrectedDataValue_toRootFiles(self, dictionary, region):

        hist = ROOT.TH1F("average_MCcorrectedData_Syst", "average_MCcorrectedData_Syst", len(self.njets), 0, len(self.njets))

        # loop over the njets
        for njet in self.njets:

            Njet = int(njet.replace("incl", ""))

            hist.SetBinContent((Njet - len(self.njets)), dictionary[njet]) 

        # put the average value of MC corrected data to the dictionary
        self.make_HiggsCombineInputs_RootFiles(hist, ("average_MCcorrectedData_Syst_%s"%(region)))


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

            hist.SetBinContent((Njet - len(self.njets)), dictionary[njet])

        # put the mximum value of MC corrected data to the dicionary
        self.make_HiggsCombineInputs_RootFiles(hist, ("maximum_MCcorrectedData_Syst_%s"%(region)))


    # -------------------
    # close the root file
    # -------------------
    def close_HiggsCombineInputs_RootFiles(self):

        self.f.Close()







    



