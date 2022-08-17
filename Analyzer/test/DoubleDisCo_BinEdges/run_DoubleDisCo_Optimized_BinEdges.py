import ROOT
import os
import argparse

from DoubleDisCo_Regions     import *
from DoubleDisCo_Plotter     import *
from DoubleDisCo_TableWriter import *

# ---------------------------------------------------
# get list of all possible optimized ABCD edges
#   -- put them in a list
#   -- choice the same edges for all njets bins later
# ---------------------------------------------------
def get_Optimized_ABCDedges(all_ABCDEdges, significance, nonClosure, nonClosure_Pull, sigFracB, sigFracC, sigFracD, sigFracB_Unc, sigFracC_Unc, sigFracD_Unc, sigFracCut):

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

            # make the table to compare the significances between 
            if math.sqrt(total_significance) > 3.5:

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

# ---------------------------------------------------
# get a table to see how the significance differences 
# between non-optimized and optimized ABCD bin edges
#   -- non-closure, pull, sigFracs and bin edges
# ---------------------------------------------------


def main():

    # --------------------------------------------------------------------------------------------------------------------------------------------------
    # command to run this script:
    #   -- python run_DoubleDisCo_Optimized_BinEdges.py --year 2016 --tt TT --nontt NonTT --ttVar TT_erdOn --data Data --sig RPV --mass 550 --channel 0l
    #   -- python run_DoubleDisCo_Optimized_BinEdges.py --year 2016 --tt TT --nontt NonTT --ttVar TT_erdOn --data Data --sig RPV --mass 550 --channel 1l
    # -------------------------------------------------------------------------------------------------------------------------------------------------- 
    usage  = "usage: %prog [options]"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("--year",           dest="year",           help="which year",            required=True)
    parser.add_argument("--path",           dest="path",           help="Input dir with histos", default="/uscms_data/d3/jhiltb/PO_Boxes/shared/2016_DisCo_0L_Cand1_1L") # Both 0L & 1l with OldSeed
    parser.add_argument("--tt",             dest="tt",             help="TT",                    required=True)
    parser.add_argument("--nontt",          dest="nontt",          help="NonTT",                 required=True)
    parser.add_argument("--ttVar",          dest="ttVar",          help="TT Variances",          required=True)
    parser.add_argument("--sig",            dest="sig",            help="RPV, SYY",              default="RPV")
    parser.add_argument("--mass",           dest="mass",           help="signal mass",           default="550")
    parser.add_argument("--data",           dest="data",           help="JetHT, SingleLepton",   required=True)
    parser.add_argument("--channel",        dest="channel",        help="0l or 1l",              required=True)
    parser.add_argument("--metric",         dest="metric",         help="NN,New",                required=None, type=str)
    parser.add_argument("--fixedDisc1edge", dest="fixedDisc1edge", help="fixed d1 edge",         default=None, type=float)
    parser.add_argument("--fixedDisc2edge", dest="fixedDisc2edge", help="fixed d2 edge",         default=None, type=float)
    parser.add_argument("--numChoices",     dest="numChoices",     help="number of edge choices",default=4, type=int)
    args = parser.parse_args()

    Sig     = "%s_%s"%(args.sig, args.mass)
    ttVar   = "%s"%(args.ttVar)
    samples = [args.tt, args.nontt, ttVar, Sig, args.data]

    # get SYY histogram
    modelDecay = "2t6j"
    if ("SHH" in args.sig):
        modelDecay = "2t4b"

    # root files
    files = {
        "TT"                   : ROOT.TFile.Open(args.path + "/" + args.year + "_TT.root"),
        "TT_erdOn"             : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_erdOn.root"),
        "TT_fsrDown"           : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_fsrDown.root"),
        "TT_fsrUp"             : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_fsrUp.root"),
        "TT_isrDown"           : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_isrDown.root"),
        "TT_isrUp"             : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_isrUp.root"),
        "TT_hdampDown"         : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_hdampDown.root"),
        "TT_hdampUp"           : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_hdampUp.root"),
        "TT_JECdown"           : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_JECdown.root"),
        "TT_JECup"             : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_JECup.root"),
        "TT_JERdown"           : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_JERdown.root"),
        "TT_JERup"             : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_JERup.root"),
        "TT_underlyingEvtDown" : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_underlyingEvtDown.root"),
        "TT_underlyingEvtUp"   : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_underlyingEvtUp.root"), 
        "TT_UL"                : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_UL.root"), 
        "NonTT"                : ROOT.TFile.Open(args.path + "/" + args.year + "_Non_TT.root"),
        #"NonTT_withoutQCD"     : ROOT.TFile.Open(args.path + "/" + args.year + "_NonTT_withoutQCD.root"),
        "Data"                 : ROOT.TFile.Open(args.path + "/" + args.year + "_Data.root"),
        Sig                    : ROOT.TFile.Open(args.path + "/" + args.year + "_%s_%s_mStop-%s.root"%(args.sig,modelDecay,args.mass)),
    }

    # get the 2D histograms 
    njets = None
    histNames = "h_DoubleDisCo_disc1_disc2_%s_Njets"%(args.channel)
    njets = ["7", "8", "9", "10", "11", "12incl"]

    # make regionis list for adding all edges to DoubleDisCo cfg file
    regions = {"ABCD"        : "ABCD",
               "Val_bdEF"    : "bdEF",
               "Val_cdiGH"   : "cdGH",
               "Val_subDivD" : "subDivD",
    }

    # initialize the dictionaries of any regions
    translator = {"ABCD"          : {"A" : "A",  "B" : "B",  "C" : "C",  "D" : "D" },
                  "Val_bdEF"      : {"A" : "b",  "B" : "E",  "C" : "d",  "D" : "F" },
                  "Val_cdiGH"     : {"A" : "c",  "B" : "di", "C" : "G",  "D" : "H" },
                  "Val_subDivD"   : {"A" : "dA", "B" : "dB", "C" : "dC", "D" : "dD"},
    }

    # initialize table path
    thePath = "table_optimizing_ABCDedges"
    
    if not os.path.exists(thePath):
        os.makedirs(thePath)

    optimized_ABCDedges    = Optimized_ABCDedges(thePath, args.channel, args.year, "ABCDedges",    args.sig, args.mass)

    # ---------------------------------------------------------
    # to optimize the ABCD edges
    # initialize the Njet dictionaries which include quantities 
    # ---------------------------------------------------------
    all_ABCDEdges = {}; significance = {}; significance_includingNonClosure = {}; significance_nonSimplified  = {}; nonClosure = {}; nonClosure_Pull = {}
    sigFracB      = {}; sigFracB_Unc = {}; sigFracC   = {}; sigFracC_Unc    = {}; sigFracD = {}; sigFracD_Unc = {}

    for key, region in regions.items():
        all_ABCDEdges[key] = {}; significance[key] = {}; significance_includingNonClosure[key] = {}; significance_nonSimplified[key] = {}; nonClosure[key] = {}; nonClosure_Pull[key] = {}
        sigFracB[key]      = {}; sigFracB_Unc[key] = {}; sigFracC[key]                         = {}; sigFracC_Unc[key]               = {}; sigFracD[key]   = {}; sigFracD_Unc[key]    = {}
 
    # ---------------
    # loop over njets
    # --------------- 
    for njet in njets: 

        hist_lists = {}
        for sample in samples:
            hist_lists[sample] = files[sample].Get(histNames + njet)

        minEdge  = hist_lists["TT"].GetXaxis().GetBinLowEdge(1) 
        maxEdge  = hist_lists["TT"].GetXaxis().GetBinLowEdge(hist_lists["TT"].GetNbinsX()+1)
        binWidth = hist_lists["TT"].GetXaxis().GetBinWidth(1)

        # ---------------------
        # loop over the regions
        # ---------------------
        for key, region in regions.items():
        
            # -----------------
            # make ABCD regions
            # -----------------
            theEdgesClass = None
            if key == "ABCD":

                theEdgesClass = All_Regions(hist_lists, Sig=Sig, ttVar=ttVar, fixedDisc1Edge=args.fixedDisc1edge, fixedDisc2Edge=args.fixedDisc2edge, metric=args.metric)    

                # ----------------------------------------------------
                # fill all the dictionaries to optimize the ABCD edges
                # ----------------------------------------------------
                all_ABCDEdges[key][njet]                    = theEdgesClass.get("edges",                                     None, None, "TT") 
                significance[key][njet]                     = np.array(theEdgesClass.get("significance",                     None, None, "TT"))[:,0]
                significance_includingNonClosure[key][njet] = np.array(theEdgesClass.get("significance_includingNonClosure", None, None, "TT"))[:,0] # significance including non-closure
                significance_nonSimplified[key][njet]       = np.array(theEdgesClass.get("significance_nonSimplified",       None, None, "TT"))[:,0] # significance, non-simplified version 
                nonClosure[key][njet]                       = np.array(theEdgesClass.get("nonClosure",                       None, None, "TT"))[:,0] 
                nonClosure_Pull[key][njet]                  = np.array(theEdgesClass.get("pull",                             None, None, "TT"))[:,0]
                sigFracB[key][njet]                         = np.array(theEdgesClass.get("sigFractionB",                     None, None, Sig ))[:,0]
                sigFracB_Unc[key][njet]                     = np.array(theEdgesClass.get("sigFractionB",                     None, None, Sig ))[:,1]
                sigFracC[key][njet]                         = np.array(theEdgesClass.get("sigFractionC",                     None, None, Sig ))[:,0]
                sigFracC_Unc[key][njet]                     = np.array(theEdgesClass.get("sigFractionC",                     None, None, Sig ))[:,1]
                sigFracD[key][njet]                         = np.array(theEdgesClass.get("sigFractionD",                     None, None, Sig ))[:,0]
                sigFracD_Unc[key][njet]                     = np.array(theEdgesClass.get("sigFractionD",                     None, None, Sig ))[:,1]
  
    sigFracsCut = 0.4 
    if args.channel == "0l": 
        sigFracsCut = 0.3

    # --------------------------------------------------
    # optimized ABCD edges with different significances:
    #   -- basic version
    #   -- non-simplified version
    # --------------------------------------------------
    #opt_ABCDEdges = get_Optimized_ABCDedges(all_ABCDEdges["ABCD"]["7"], significance["ABCD"], nonClosure["ABCD"], nonClosure_Pull["ABCD"], sigFracB["ABCD"], sigFracC["ABCD"], sigFracD["ABCD"], sigFracB_Unc["ABCD"], sigFracC_Unc["ABCD"], sigFracD_Unc["ABCD"], sigFracsCut) 
    #opt_ABCDEdges_2 = get_Optimized_ABCDedges(all_ABCDEdges["ABCD"]["7"], significance_nonSimplified["ABCD"], nonClosure["ABCD"], nonClosure_Pull["ABCD"], sigFracB["ABCD"], sigFracC["ABCD"], sigFracD["ABCD"], sigFracB_Unc["ABCD"], sigFracC_Unc["ABCD"], sigFracD_Unc["ABCD"], sigFracsCut)    
    
    #print "ABCD edges with Significance                      : ", opt_ABCDEdges
    #print "ABCD edges with Significance nonSimplified        : ", opt_ABCDEdges_2

    # ------------------------------------------------------------
    # optimized ABCD edges with significance including non-closure
    # ------------------------------------------------------------
    significance, opt_ABCDEdges_1, best_choices_parameters = get_Optimized_ABCDedges(all_ABCDEdges["ABCD"]["7"], significance_includingNonClosure["ABCD"], nonClosure["ABCD"], nonClosure_Pull["ABCD"], sigFracB["ABCD"], sigFracC["ABCD"], sigFracD["ABCD"], sigFracB_Unc["ABCD"], sigFracC_Unc["ABCD"], sigFracD_Unc["ABCD"], sigFracsCut)

    optimized_ABCDedges.writeLine(Njets = njets, Best_Parameters = best_choices_parameters, numEdgeChoices = args.numChoices)

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


if __name__ == '__main__':
    main()
