# The HiggsCombineInputs class 
# it makes the inputs (root files, txt files, etc.) 
# to the Higgs combine
# we can add both sys and nuisance parameters to this class
import ROOT

class HiggsCombineInputs:

    # --------------------
    # Initialize the class 
    # --------------------
    def __init__(self, year, njets, channel, samples, regions):

        self.year       = year
        self.njets      = njets
        self.channel    = channel
        self.samples    = samples 
        self.regions    = regions

        # ----------------
        # open a root file
        # ----------------
        fileName = year + "_" + "TT_TTvar_sys" + "_" + channel + ".root" 
        self.f   = ROOT.TFile.Open(fileName, "RECREATE")        

    # -------------------------------------------------------
    # make a root file to include:
    #   -- MC correction values for TT in signal region 
    #   -- MC correction values for TTvar in signal region
    #   -- MC correction factor ratio: TT/TTvar (in MC level)
    # -------------------------------------------------------
    def make_HiggsCombineInputs_RootFiles(self, TH1, varName): 

        self.f.WriteObject(TH1, "%s"%(varName))

    # ----------------------------------
    # put the variables to the root file
    # ----------------------------------
    def put_HiggsCombineInputs_toRootFiles(self, dictionary):

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

                # put the variables to root file
                self.make_HiggsCombineInputs_RootFiles(hist, (var + "_" + ttVar))
                        
    # ----------------------------------------------------------
    # add the average MC corrected data values for all njet bins 
    # to the same root file like for MC corrections included 
    # ----------------------------------------------------------
    def put_averageCorrectedDataValue_toRootFiles(self, dictionary):

        hist = ROOT.TH1F("average_MCcorrectedData_Value", "average_MCcorrectedData_Value", len(self.njets), 0, len(self.njets))

        # loop over the njets
        for njet in self.njets:

            Njet = int(njet.replace("incl", ""))

            hist.SetBinContent((Njet - len(self.njets)), dictionary[njet]) 

        # put the average value of MC corrected data to the dictionary
        self.make_HiggsCombineInputs_RootFiles(hist, ("average_MCcorrectedData_Value"))

    # ---------------
    # close root file
    # ---------------
    def close_HiggsCombineInputs_RootFiles(self):

        self.f.Close()







    



