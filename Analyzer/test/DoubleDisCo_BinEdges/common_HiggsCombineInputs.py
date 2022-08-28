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
        fileName = path.replace("year",year) + year + "_" + "TT_TTvar_sys" + ".root"  
        self.f   = ROOT.TFile.Open(fileNAme, "RECREATE")        

    # -------------------------------------------------------
    # make a root file to include:
    #   -- MC correction values for TT in signal region 
    #   -- MC correction values for TTvar in signal region
    #   -- MC correction factor ratio: TT/TTvar (in MC level)
    # and close root file 
    # -------------------------------------------------------
    def make_HiggsCombineInputs_RootFiles(self, TH1, varName) 

        self.f.WriteObject(TH1, "%s"%(varName))

    def close_HiggsCombineInputs_RootFiles(self)

        self.f.Close()

