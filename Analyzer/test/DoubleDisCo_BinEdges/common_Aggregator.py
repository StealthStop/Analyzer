# The Aggregator class is fed Regions objects, extracts useful information from them
# and stores everything in a single dictionary, allowing information to be retreived
# in bulk
class Aggregator:

    # -------------------------------------------------------------------------
    # Initialize a master dictionary e.g. "data", to hold lots of useful things
    # ------------------------------------------------------------------------- 
    def __init__(self, samples, njets, regions, boundaries):

        self.data       = {}
        self.samples    = samples + ["TTinData"]
        self.njets      = njets
        self.boundaries = boundaries
        self.regions    = regions
        self.subregions = ["A", "B", "C", "D"]

    # ------------------------------------------------------
    # Construct a specifically formatted string that is used
    # for accessing information from "data"
    # ------------------------------------------------------
    def makeKey(self, variable = None, **kwargs):

        theKey = "" 

        if "region"   in kwargs: theKey += kwargs["region"]   + "_"
        if "njet"     in kwargs: theKey += kwargs["njet"]     + "_"
        if "boundary" in kwargs: theKey += "%.2f"%(kwargs["boundary"]) + "_"
        if "sample"   in kwargs: theKey += kwargs["sample"]   + "_"

        if "variable" != None:   theKey += variable           + "_"
  
        return theKey[:-1]

    # --------------------------------------------------------    
    # Main method to eat an instance of the regions class
    # and grab all the things we want from it for safe keeping
    # --------------------------------------------------------
    def aggregate(self, regionObj, **kwargs):

        # get all final bin edges
        self.data[self.makeKey(variable = "edges",      **kwargs)] = regionObj.get("edges",      None, None, "TT")
        self.data[self.makeKey(variable = "finalEdges", **kwargs)] = regionObj.getFinal("edges",             "TT")

        Sig = None
        for sample in self.samples:
        
            if Sig != None:
                break

            for s in ["RPV", "SYY", "SHH"]:
                if s in sample:
                    Sig = sample
                    break

        # quantities and variables with any combination of bin edges
        for subregion in self.subregions:

            self.data[self.makeKey(variable = "sigFractions%s"%(subregion), **kwargs)] = regionObj.get("sigFraction%s"%(subregion), None, None, Sig)
            self.data[self.makeKey(variable = "ttFractions%s"%(subregion),  **kwargs)] = regionObj.get("ttFraction%s"%(subregion),  None, None, "TT")

            # quantities and variables with the final choice of bin edges
            self.data[self.makeKey(variable = "sigFraction%s"%(subregion), **kwargs)] = regionObj.getFinal("sigFraction%s"%(subregion), Sig)
            self.data[self.makeKey(variable = "ttFraction%s"%(subregion),  **kwargs)] = regionObj.getFinal("ttFraction%s"%(subregion), "TT")

        # loop over for getting plots for TT, NonTT, Data
        for sample in self.samples:

            # Skip all ttVars not corresponding to the one currently used by regionObj
            if "TT_" in sample and sample not in regionObj.get_ttVar_Name():
                continue

            for subregion in self.subregions:
                self.data[self.makeKey(variable = "nEvents%s"%(subregion),  sample = sample, **kwargs)] = regionObj.getFinal("nEvents%s"%(subregion), sample)

            if not any(s in sample for s in ["RPV", "SYY", "SHH"]):
                self.data[self.makeKey(variable = "nonClosures",    sample = sample, **kwargs)] = regionObj.get("nonClosure",       None, None, sample) # vars with any combination of bin edges
                self.data[self.makeKey(variable = "pulls",          sample = sample, **kwargs)] = regionObj.get("pull",             None, None, sample)
                self.data[self.makeKey(variable = "closureCorrs",   sample = sample, **kwargs)] = regionObj.get("closureCorr",      None, None, sample)
                self.data[self.makeKey(variable = "Closure",        sample = sample, **kwargs)] = regionObj.getFinal("Closure",                 sample)
                self.data[self.makeKey(variable = "nonClosure",     sample = sample, **kwargs)] = regionObj.getFinal("nonClosure",              sample) # vars with the final choice of bin edges
                self.data[self.makeKey(variable = "pull",           sample = sample, **kwargs)] = regionObj.getFinal("pull",                    sample)
                self.data[self.makeKey(variable = "closureCorr",    sample = sample, **kwargs)] = regionObj.getFinal("closureCorr",             sample) # vars with the final choice of bin edges
                if sample == "TTinData":
                    self.data[self.makeKey(variable = "CorrectedDataClosure",       sample = sample,                                **kwargs)] = regionObj.getFinal("CorrectedDataClosure",       sample)
                    self.data[self.makeKey(variable = "ttVar_CorrectedDataClosure", sample = sample+"_"+regionObj.get_ttVar_Name(), **kwargs)] = regionObj.getFinal("ttVar_CorrectedDataClosure", sample)

            if "TT_" in sample:
                self.data[self.makeKey(variable = "MCcorrRatio_MC", sample = sample, **kwargs)] = regionObj.getFinal("MCcorrRatio_MC",          sample)

        self.data[self.makeKey(variable = "significances", **kwargs)] = regionObj.get("significance",      None, None, "TT")
        self.data[self.makeKey(variable = "significance",  **kwargs)] = regionObj.getFinal("significance",             "TT")

    # --------------------------------
    # Get variables for each Njets bin
    # --------------------------------
    def getPerNjets(self, variable, **kwargs):

        payload = {}
        for njet in self.njets:

            masterKey = self.makeKey(variable = variable, njet = njet, **kwargs)
            payload[njet] = self.data[masterKey]
              
        return payload

    # -------------------------------
    # Get variables for each boundary
    # -------------------------------
    def getPerBoundary(self, variable, **kwargs):
        
        payload = {}
        newKwargs = kwargs.copy()
        sample = newKwargs.pop("sample", "None")

        for boundary in self.boundaries:

            masterKey = self.makeKey(variable = variable, boundary = boundary, **kwargs)
            chefKey   = self.makeKey(variable = "sigFractionA", boundary = boundary, **newKwargs)

            if masterKey not in self.data:
                #print("Skipping key \"%s\""%(masterKey))
                continue

            # this statement for data and data/MC closure correction
            sigFracA = self.data[chefKey][0]
            if (sigFracA >= 0.05 and "TTinData" in sample):
                payload[round(boundary, 2)] = (-999.0, 0.0)

            else: 
                payload[round(boundary, 2)] = self.data[masterKey]

        return payload 

    # -------------
    # Get variables
    # -------------
    def get(self, variable, **kwargs):
       
        masterKey = self.makeKey(variable = variable, **kwargs)

        return self.data[masterKey]
