from common_Regions     import *
from common_TableWriter import *

# ---------------------------------------------------
# get a table to see how the significance differences 
# between non-optimized and optimized ABCD bin edges
#   -- non-closure, pull, sigFracs and bin edges
# ---------------------------------------------------
class Optimized_BinEdges():

    def __init__(self, year, channel, model, mass, ttVar, translator):
        
        self.year       = year
        self.channel    = channel
        self.sig        = model
        self.mass       = mass
        self.ttVar      = ttVar
        self.translator = translator

    # ---------------------------------------------------
    # get list of all possible optimized ABCD edges
    #   -- put them in a list
    #   -- choice the same edges for all njets bins later
    # ---------------------------------------------------
    def get_Optimized_ABCDedges(self, all_ABCDEdges, significance, nonClosure, nonClosure_Pull, sigFracB, sigFracC, sigFracD, sigFracB_Unc, sigFracC_Unc, sigFracD_Unc, sigFracCut):
    
        njets               = nonClosure.keys()
        max_significance    = 0.0
        i_bestChoice        = -999
        best_choices_params = {}
    
        for i_list in range(0, len(all_ABCDEdges)):
    
            if ( float(all_ABCDEdges[i_list][0]) < 0.5  or  float(all_ABCDEdges[i_list][1]) < 0.5 ): continue
            if ( float(all_ABCDEdges[i_list][0]) > 0.95 or  float(all_ABCDEdges[i_list][1]) > 0.95 ): continue
    
            total_significance = 0.0; pass_nonClosure = True; pass_sigFrac = True; 

            for njet in njets:
    
                # get significance for each njet bin and add them quadrature 
                total_significance += (significance[njet][i_list]**2.0)
 
                # check the non-closure and pull to get the same ABCD edges for all njet bins
                if not (nonClosure[njet][i_list] < 0.3 or abs(nonClosure_Pull[njet][i_list]) < 2): 
                    pass_nonClosure = False
    
                # based on the sigFracs table, 30% for 0l, 40% for 1l
                if not (sigFracB[njet][i_list] < sigFracCut and sigFracC[njet][i_list] < sigFracCut and sigFracD[njet][i_list] < sigFracCut): 
                    pass_sigFrac = False
    
            # Add manually the non-optimized edges and associated quantities
            if all_ABCDEdges[i_list][0] == "0.600" and all_ABCDEdges[i_list][1] == "0.600":

                best_choices_params[999.0] = {}
            
                best_choices_params[999.0]["Significance"] = math.sqrt(total_significance)
                
                for njet in njets:
    
                    best_choices_params[999.0]["nonClosure_njet%s"%njet]      = nonClosure[njet][i_list]
                    best_choices_params[999.0]["nonCllosurePull_njet%s"%njet] = nonClosure_Pull[njet][i_list]
                    best_choices_params[999.0]["sigFracB_njet%s"%njet]        = sigFracB[njet][i_list]
                    best_choices_params[999.0]["sigFracC_njet%s"%njet]        = sigFracC[njet][i_list]  
                    best_choices_params[999.0]["sigFracD_njet%s"%njet]        = sigFracD[njet][i_list]
                    
                best_choices_params[999.0]["ABCDedges"] = all_ABCDEdges[i_list]
           
            # get the current best choice of ABCD edges     
            if (pass_nonClosure and pass_sigFrac):
    
                if (total_significance > max_significance):

                    max_significance = total_significance
                    i_bestChoice     = i_list

                    #print "maximum significance : ", max_significance
                    #print "disc1, disc2         : ", all_ABCDEdges[i_bestChoice] 
 
                # make the table to compare the significances between 
                if math.sqrt(total_significance) > 1.0:
    
                    best_choices_params[total_significance] = {}
            
                    best_choices_params[total_significance]["Significance"] = math.sqrt(total_significance)
                    
                    for njet in njets:
    
                        best_choices_params[total_significance]["nonClosure_njet%s"%njet]      = nonClosure[njet][i_list]
                        best_choices_params[total_significance]["nonCllosurePull_njet%s"%njet] = nonClosure_Pull[njet][i_list]
                        best_choices_params[total_significance]["sigFracB_njet%s"%njet]        = sigFracB[njet][i_list]
                        best_choices_params[total_significance]["sigFracC_njet%s"%njet]        = sigFracC[njet][i_list]  
                        best_choices_params[total_significance]["sigFracD_njet%s"%njet]        = sigFracD[njet][i_list]
                        
                    best_choices_params[total_significance]["ABCDedges"] = all_ABCDEdges[i_list]
    
        return max_significance, all_ABCDEdges[i_bestChoice], best_choices_params

    # --------------
    # run the module
    # -------------- 
    def run(self, disc1edge=None, disc2edge=None, fastMode=False, **kwargs):

        tablesPath     = kwargs["tablesPath"]["TT"]
        regions        = kwargs["regions"]
        njets          = kwargs["njets"]
        samples        = kwargs["samples"]
        histName       = kwargs["histName"]
        files          = kwargs["files"]
        numEdgeChoices = kwargs["numEdgeChoices"]

        optimized_ABCDedges = Optimized_ABCDedges(tablesPath, self.channel, self.year, "ABCDedges", self.sig)

        # ---------------------------------------------------------
        # to optimize the ABCD edges
        # initialize the Njet dictionaries which include quantities 
        # ---------------------------------------------------------
        all_ABCDEdges = {}; significance = {}; significance_includingNonClosure = {}; significance_nonSimplified  = {}; nonClosure = {}; nonClosure_Pull = {}
        sigFracB      = {}; sigFracB_Unc = {}; sigFracC   = {}; sigFracC_Unc    = {}; sigFracD = {}; sigFracD_Unc = {}

        for region in regions:
            all_ABCDEdges[region] = {}; significance[region] = {}; significance_includingNonClosure[region] = {}; significance_nonSimplified[region] = {}; nonClosure[region] = {}; nonClosure_Pull[region] = {}
            sigFracB[region]      = {}; sigFracB_Unc[region] = {}; sigFracC[region]                         = {}; sigFracC_Unc[region]               = {}; sigFracD[region]   = {}; sigFracD_Unc[region]    = {}
 
        # ---------------
        # loop over njets
        # --------------- 
        for njet in njets: 

            hist_lists = {}

            for sample in samples:

                # get the fsr/isr higtograms from TT root file
                ttvarStr = ""

                if sample == "TT_fsrDown":
                    ttvarStr = "_fsrDown"

                elif sample == "TT_fsrUp":
                    ttvarStr = "_fsrUp"

                elif sample == "TT_isrDown":
                    ttvarStr = "_isrDown"

                elif sample == "TT_isrUp":
                    ttvarStr = "_isrUp"

                hist_lists[sample] = files[sample].Get(histName.replace("${NJET}", njet) + ttvarStr)


            minEdge  = hist_lists["TT"].GetXaxis().GetBinLowEdge(1) 
            maxEdge  = hist_lists["TT"].GetXaxis().GetBinLowEdge(hist_lists["TT"].GetNbinsX()+1)
            binWidth = hist_lists["TT"].GetXaxis().GetBinWidth(1)

            # ---------------------
            # loop over the regions
            # ---------------------
            for region in regions:
            
                # -----------------
                # make ABCD regions
                # -----------------
                theEdgesClass = None
                if region == "ABCD":

                    theEdgesClass = All_Regions(hist_lists, Sig=self.sig, ttVar=self.ttVar, disc1Edge=disc1edge, disc2Edge=disc2edge)    

                    # ----------------------------------------------------
                    # fill all the dictionaries to optimize the ABCD edges
                    # ----------------------------------------------------
                    all_ABCDEdges[region][njet]                    = theEdgesClass.get("edges",                                     None, None, "TT") 
                    significance[region][njet]                     = np.array(theEdgesClass.get("significance",                     None, None, "TT"))[:,0]
                    significance_includingNonClosure[region][njet] = np.array(theEdgesClass.get("significance_includingNonClosure", None, None, "TT"))[:,0] # significance including non-closure
                    significance_nonSimplified[region][njet]       = np.array(theEdgesClass.get("significance_nonSimplified",       None, None, "TT"))[:,0] # significance, non-simplified version 
                    nonClosure[region][njet]                       = np.array(theEdgesClass.get("nonClosure",                       None, None, "TT"))[:,0] 
                    nonClosure_Pull[region][njet]                  = np.array(theEdgesClass.get("pull",                             None, None, "TT"))[:,0]
                    sigFracB[region][njet]                         = np.array(theEdgesClass.get("sigFractionB",                     None, None, self.sig ))[:,0]
                    sigFracB_Unc[region][njet]                     = np.array(theEdgesClass.get("sigFractionB",                     None, None, self.sig ))[:,1]
                    sigFracC[region][njet]                         = np.array(theEdgesClass.get("sigFractionC",                     None, None, self.sig ))[:,0]
                    sigFracC_Unc[region][njet]                     = np.array(theEdgesClass.get("sigFractionC",                     None, None, self.sig ))[:,1]
                    sigFracD[region][njet]                         = np.array(theEdgesClass.get("sigFractionD",                     None, None, self.sig ))[:,0]
                    sigFracD_Unc[region][njet]                     = np.array(theEdgesClass.get("sigFractionD",                     None, None, self.sig ))[:,1]
  
        sigFracsCut = 0.4 
        if self.channel == "0l": 
            sigFracsCut = 0.3

        # ------------------------------------------------------------
        # optimized ABCD edges with significance including non-closure
        # ------------------------------------------------------------
        significance, opt_ABCDEdges_1, best_choices_parameters = self.get_Optimized_ABCDedges(all_ABCDEdges["ABCD"]["7"], significance_includingNonClosure["ABCD"], nonClosure["ABCD"], nonClosure_Pull["ABCD"], sigFracB["ABCD"], sigFracC["ABCD"], sigFracD["ABCD"], sigFracB_Unc["ABCD"], sigFracC_Unc["ABCD"], sigFracD_Unc["ABCD"], sigFracsCut)

        optimized_ABCDedges.writeLine(Njets = njets, Best_Parameters = best_choices_parameters, numEdgeChoices = numEdgeChoices)

        # ----------------
        # close the tables
        # ----------------
        optimized_ABCDedges.writeClose()

        # -------------------------------------
        # add all edges to DoubleDisCo cfg file
        # -------------------------------------
        #if args.updateDisCoCfg:
        #    addEdges = addEdges_toDoubleDisco(args.year, args.sig, args.mass, args.channel, regions)
        #    addEdges.addEdges_toDoubleDiscoCfg(edgesPerNjets, njets)
