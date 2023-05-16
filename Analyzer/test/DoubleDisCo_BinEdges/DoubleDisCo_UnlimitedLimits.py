import ROOT
import os
import sys
import ctypes
import numpy as np

from argparse import ArgumentParser
from common_Regions import All_Regions
from run_DoubleDisCo_Validation import BryansHack
from os import sys, path
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
from runQCDCRPrediction import get_QCDCR


#ROOT.TH1.AddDirectory(False)

# Get necessary histograms from files
def load_hists(files, sig, ch, njet, var_dict):

    hists = {}

    for (sample, f) in files.items():
        
        histName = "h_DoubleDisCo_{}_disc1_disc2_{}_Njets{}_ABCD"
        if ch == "2l" and sample == "QCD": continue
        hists[sample] = f.Get(histName.format(sig, ch, njet))

        if sig in sample or "TTX" in sample or "BG_OTHER" in sample:

            varHistName = "h_DoubleDisCo_{}_disc1_disc2_{}_Njets{}_ABCD_{}"
            for var in var_dict[ch]:
                if var == "": continue
                hists[sample + var] = f.Get(varHistName.format(sig, ch, njet, var))

        #if "TT_" not in sample:
        QCDCRHistName = "h_DoubleDisCo_{}_disc1_disc2_{}_QCDCR_Njets{}_ABCD"
        hists[sample + "_QCDCR"] = f.Get(QCDCRHistName.format(sig, ch, njet))

    return hists

# Instantiate the all regions object and let it find event counts
# for each bin boundary we care about
def get_edge_info(files, sig, njets, var_dict, step, min, max, QCDCRInfo):

    all_ABCDEdges = {}

    print("Getting edge info")

    for ch in njets.keys():
        
        print("Making Regions object for {}".format(ch))

        all_ABCDEdges[ch] = {}
        
        for nj in njets[ch]:

            njet = str(nj)
            if nj + 1 not in njets[ch]:
                njet += "incl" 

            hists = load_hists(files, sig.split("_")[0], ch, njet, var_dict)
            
            region = All_Regions(hists, Sig=sig, step=step, ttVar="TT", QCDCRInfo=QCDCRInfo[ch], leftBoundary=None, rightBoundary=None, topBoundary=None, bottomBoundary=None, fastMode=True, binStart=min, binEnd=max)

            all_ABCDEdges[ch][nj] = region

    return all_ABCDEdges

# Make a histogram for each of the ABCD bin edges with specified granularity and starting/finishing points
def make_ABCD_hists(all_ABCDEdges, files, sig, step, start, finish, njets, var_dict, ttSys, inpath, outpath, year):

    signal = sig.split("_")[0]

    if not os.path.isdir("inputsAll"):
        os.makedirs("inputsAll")

    print("Making hists")

    disc_range = range(0, int(round((finish - start) / step)))

    for process in files.keys():

        if "TT" in process and process is not "TTX" and process is not "NonTT":
            procFile = ROOT.TFile.Open("inputsAll/Run2UL_TT.root".format(process), "UPDATE")
        else:
            procFile = ROOT.TFile.Open("inputsAll/Run2UL_{}.root".format(process), "Update")
        procFile.cd()

        for d1 in disc_range:
            
            for d2 in disc_range:

                disc1 = d1 * step + start
                disc2 = d2 * step + start


                for ch in njets.keys():
                  
                    if process == "QCD" and ch == "2l": continue
 
                    for var in var_dict[ch]:

                        var_str = var

                        if var is not "" and "TTX" not in process and "BG_OTHER" not in process and signal not in process:
                            continue
                        elif var is "":
                            var_str = var
                        else:
                            var_str = "_" + var

                        tempHist = ROOT.TH1D("h_njets_{}incl_{}_{}_ABCD{}_{}_{}".format(njets[ch][-1],signal,ch,var_str,int(100*disc1),int(100*disc2)),"h_njets_{}incl_{}_{}_ABCD{}_{}_{}".format(njets[ch][-1],signal,ch,var_str,int(100*disc1),int(100*disc2)), 24, -0.5, 23.5)
                        if "QCD" in process:
                            tempWeight = ROOT.TH1D("h_njets_{}incl_{}_{}_ABCD{}_{}_{}_weight".format(njets[ch][-1],signal,ch,var_str,int(100*disc1),int(100*disc2)),"h_njets_{}incl_{}_{}_ABCD{}_{}_{}".format(njets[ch][-1],signal,ch,var_str,int(100*disc1),int(100*disc2)), 8, 0, 8)
                        else:
                            tempWeight = ROOT.TH1D("h_njets_{}incl_{}_{}_ABCD{}_{}_{}_weight".format(njets[ch][-1],signal,ch,var_str,int(100*disc1),int(100*disc2)),"h_njets_{}incl_{}_{}_ABCD{}_{}_{}".format(njets[ch][-1],signal,ch,var_str,int(100*disc1),int(100*disc2)), 5, 0, 5)
                        first_jet = 0

                        if "TT" in process and process is not "TTX" or process is not "NonTT":
                            tempClosureSys = ROOT.TH1D("Run2UL_maximum_MCcorrectedData_Syst_All_{}_{}{}_{}_{}".format(signal,ch,var_str,int(100*disc1),int(100*disc2)), "Run2UL_maximum_MCcorrectedData_Syst_All_{}_{}{}_{}_{}".format(signal,ch,var_str,int(100*disc1),int(100*disc2)), 6, 0, 6)
                            tempMCcorr = ROOT.TH1D("Run2UL_MCcorr_TT_TT{}_{}_{}_{}_{}".format(var,signal,ch,int(100*disc1),int(100*disc2)), "Run2UL_MCcorr_TT_TT{}_{}_{}_{}_{}".format(var,signal,ch,int(100*disc1),int(100*disc2)), 6, 0, 6)

                        tempQCDCRHist = ROOT.TH1D("h_njets_{}incl_{}_{}_QCDCR_ABCD{}_{}_{}".format(njets[ch][-1],signal,ch,var_str,int(100*disc1),int(100*disc2)),"h_njets_{}incl_{}_{}_QCDCR_ABCD{}_{}_{}".format(njets[ch][-1],signal,ch,var_str,int(100*disc1),int(100*disc2)), 24, -0.5, 23.5)

                        for nj in njets[ch]: 

                            if first_jet == 0:
                                first_jet = nj

                            bin_val_A, bin_error_A = all_ABCDEdges[ch][nj].get("nEventsA", "{:.3f}".format(disc1), "{:.3f}".format(disc2), process+var)
                            bin_val_B, bin_error_B = all_ABCDEdges[ch][nj].get("nEventsB", "{:.3f}".format(disc1), "{:.3f}".format(disc2), process+var)
                            bin_val_C, bin_error_C = all_ABCDEdges[ch][nj].get("nEventsC", "{:.3f}".format(disc1), "{:.3f}".format(disc2), process+var)
                            bin_val_D, bin_error_D = all_ABCDEdges[ch][nj].get("nEventsD", "{:.3f}".format(disc1), "{:.3f}".format(disc2), process+var)

                            tempHist.SetBinContent(nj-first_jet + 1, bin_val_A)
                            tempHist.SetBinError(nj-first_jet + 1, bin_error_A)

                            tempHist.SetBinContent(nj-first_jet + 7, bin_val_B)
                            tempHist.SetBinError(nj-first_jet + 7, bin_error_B)

                            tempHist.SetBinContent(nj-first_jet + 13, bin_val_C)
                            tempHist.SetBinError(nj-first_jet + 13, bin_error_C)

                            tempHist.SetBinContent(nj-first_jet + 19, bin_val_D)
                            tempHist.SetBinError(nj-first_jet + 19, bin_error_D)

                            # Using average weight and raw event count for MC stat uncertainty
                            # Histogram bin definitions:
                            #   1: Raw event count A
                            #   2: Raw event count B
                            #   3: Raw event count C
                            #   4: Raw event count D
                            #   5: Average weight
                            # For QCD, the binning has to be a little bit different:
                            #   1: Raw event count A (QCDCR in Data)
                            #   2: Raw event count B
                            #   3: Raw event count C
                            #   4: Raw event count D
                            #   5: TF for the A region (SR/CR)
                            #   5: TF for the B region
                            #   5: TF for the C region
                            #   5: TF for the D region
                            if "QCD" in process:
                                
                                nEvents_A = all_ABCDEdges[ch][nj].get("nEventsA", "{:.3f}".format(disc1), "{:.3f}".format(disc2), process)
                                nEvents_B = all_ABCDEdges[ch][nj].get("nEventsB", "{:.3f}".format(disc1), "{:.3f}".format(disc2), process)
                                nEvents_C = all_ABCDEdges[ch][nj].get("nEventsC", "{:.3f}".format(disc1), "{:.3f}".format(disc2), process)
                                nEvents_D = all_ABCDEdges[ch][nj].get("nEventsD", "{:.3f}".format(disc1), "{:.3f}".format(disc2), process)

                                tf_A = all_ABCDEdges[ch][nj].get("QCDTF_A", "{:.3f}".format(disc1), "{:.3f}".format(disc2), "QCD")
                                tf_B = all_ABCDEdges[ch][nj].get("QCDTF_B", "{:.3f}".format(disc1), "{:.3f}".format(disc2), "QCD")
                                tf_C = all_ABCDEdges[ch][nj].get("QCDTF_C", "{:.3f}".format(disc1), "{:.3f}".format(disc2), "QCD")
                                tf_D = all_ABCDEdges[ch][nj].get("QCDTF_D", "{:.3f}".format(disc1), "{:.3f}".format(disc2), "QCD")

                                tempWeight.SetBinContent(1, nEntries_A)
                                tempWeight.SetBinContent(2, nEntries_B)
                                tempWeight.SetBinContent(3, nEntries_C)
                                tempWeight.SetBinContent(4, nEntries_D)
                                tempWeight.SetBinContent(5, tf_A)
                                tempWeight.SetBinContent(6, tf_B)
                                tempWeight.SetBinContent(7, tf_C)
                                tempWeight.SetBinContent(8, tf_D)

                            else:
                                weight = all_ABCDEdges[ch][nj].get("weight", "{:.3f}".format(disc1), "{:.3f}".format(disc2), process+var)
                                nEntries_A = all_ABCDEdges[ch][nj].get("nEntriesA", "{:.3f}".format(disc1), "{:.3f}".format(disc2), process+var)
                                nEntries_B = all_ABCDEdges[ch][nj].get("nEntriesB", "{:.3f}".format(disc1), "{:.3f}".format(disc2), process+var)
                                nEntries_C = all_ABCDEdges[ch][nj].get("nEntriesC", "{:.3f}".format(disc1), "{:.3f}".format(disc2), process+var)
                                nEntries_D = all_ABCDEdges[ch][nj].get("nEntriesD", "{:.3f}".format(disc1), "{:.3f}".format(disc2), process+var)

                                tempWeight.SetBinContent(1, nEntries_A)
                                tempWeight.SetBinContent(2, nEntries_B)
                                tempWeight.SetBinContent(3, nEntries_C)
                                tempWeight.SetBinContent(4, nEntries_D)
                                tempWeight.SetBinContent(5, weight)

                            if process is "TT":
                                
                                mcCorr, mcCorrUnc = all_ABCDEdges[ch][nj].get("closureCorr", "{:.3f}".format(disc1), "{:.3f}".format(disc2), process)
                                if str(nj) in ttSys[ch]["({}, {})".format(disc1,disc2)].keys():
                                    ncSys = ttSys[ch]["({}, {})".format(disc1,disc2)][str(nj)]
                                else:
                                    ncSys = ttSys[ch]["({}, {})".format(disc1,disc2)]["max"]

                                tempMCcorr.SetBinContent(nj-first_jet + 1, mcCorr)
                                tempMCcorr.SetBinError(nj-first_jet + 1, mcCorrUnc)

                                tempClosureSys.SetBinContent(nj-first_jet + 1, ncSys)
                            
                            if "TT_" not in process:
                                qcdcr_bin_val_A, qcdcr_bin_error_A = all_ABCDEdges[ch][nj].get("nEventsA", "{:.3f}".format(disc1), "{:.3f}".format(disc2), process+"_QCDCR")
                                qcdcr_bin_val_B, qcdcr_bin_error_B = all_ABCDEdges[ch][nj].get("nEventsB", "{:.3f}".format(disc1), "{:.3f}".format(disc2), process+"_QCDCR")
                                qcdcr_bin_val_C, qcdcr_bin_error_C = all_ABCDEdges[ch][nj].get("nEventsC", "{:.3f}".format(disc1), "{:.3f}".format(disc2), process+"_QCDCR")
                                qcdcr_bin_val_D, qcdcr_bin_error_D = all_ABCDEdges[ch][nj].get("nEventsD", "{:.3f}".format(disc1), "{:.3f}".format(disc2), process+"_QCDCR")

                                tempQCDCRHist.SetBinContent(nj-first_jet + 1, qcdcr_bin_val_A)
                                tempQCDCRHist.SetBinError(nj-first_jet + 1, qcdcr_bin_error_A)

                                tempQCDCRHist.SetBinContent(nj-first_jet + 7, qcdcr_bin_val_B)
                                tempQCDCRHist.SetBinError(nj-first_jet + 7, qcdcr_bin_error_B)

                                tempQCDCRHist.SetBinContent(nj-first_jet + 13, qcdcr_bin_val_C)
                                tempQCDCRHist.SetBinError(nj-first_jet + 13, qcdcr_bin_error_C)

                                tempQCDCRHist.SetBinContent(nj-first_jet + 19, qcdcr_bin_val_D)
                                tempQCDCRHist.SetBinError(nj-first_jet + 19, qcdcr_bin_error_D)


                        if process is "TT" and var is "":
                            tempMCcorr.Write()
                            tempClosureSys.Write()

                        if ("TT_" not in process and var is "") or (process is "TTX" or process is "BG_OTHER" or signal in process or process is "QCD"):
                            tempQCDCRHist.Write()
                            tempWeight.Write()

                            print("Writing for {} {}".format(process, var))
                            tempHist.Write()

        procFile.Close()


    for d1 in disc_range:
        
        for d2 in disc_range:

            disc1 = d1 * step + start
            disc2 = d2 * step + start


            for ch in njets.keys():
               
                be_str = "{}_{}".format(int(100*disc1), int(100*disc2))
                print("Getting QCDCR TF for {}".format(be_str))
                get_QCDCR("./inputsAll/", "./inputsAll/", year, ch, signal, be_str)
   
def write_to_root_file(dc_histos):

    if not os.path.isdir("inputsAll"):
        os.makedirs("inputsAll")

    for proc in dc_histos.keys():
        
        procFile = ROOT.TFile.Open("inputsAll/Run2UL_{}.root".format(proc), "w")

        for histo in dc_histos[proc]:
            
            procFile.WriteObject(histo, histo.GetName())

        procFile.Close()
 

# Run the module and create all data cards for the appropriate choices of bin edges
def main():

    parser = ArgumentParser()

    parser.add_argument("--allMass",    action="store_true",      default=False,       help="Make inputs for all masses (specified signal) for Run2UL")
    parser.add_argument("--sig",        action="store", type=str, default="RPV",       help="Signal model (RPV, StealthSYY)")
    parser.add_argument("--mass",       action="store", type=str, default="300",       help="Signal mass (300-1400)")
    parser.add_argument("--year",       action="store", type=str, default="Run2UL",    help="Year prefix for inputs")
    parser.add_argument("--path",       action="store", type=str, default="/uscms/home/bcrossma/nobackup/analysis/CMSSW_10_2_13/src/CombineFits/DataCardProducer/inputs_v3.3/",    help="Input path")
    parser.add_argument("--output",   action="store", type=str, default="/uscms/home/bcrossma/nobackup/analysis/CMSSW_10_2_13/src/CombineFits/DataCardProducer/input_scan/",    help="Output path")
    parser.add_argument("--min",   action="store", type=float, default=0.1,    help="Bin edge minimum")
    parser.add_argument("--max",   action="store", type=float, default=1.0,    help="Bin edge maximum")
    parser.add_argument("--step",   action="store", type=float, default=0.1,    help="Bin edge step")

    args = parser.parse_args()

    files = {
        "TT"             : ROOT.TFile.Open(args.path + "/" + args.year + "_TT.root"            ), # TT.root file includes also fsrUp/Down, isrUp/Down JECup/down, JERup/down histograms anymore
        "TT_fsrDown"     : ROOT.TFile.Open(args.path + "/" + args.year + "_TT.root"            ),
        "TT_fsrUp"       : ROOT.TFile.Open(args.path + "/" + args.year + "_TT.root"            ),
        "TT_isrDown"     : ROOT.TFile.Open(args.path + "/" + args.year + "_TT.root"            ),
        "TT_isrUp"       : ROOT.TFile.Open(args.path + "/" + args.year + "_TT.root"            ),
        "TT_erdON"       : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_erdON.root"      ),
        "TT_hdampDOWN"   : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_hdampDOWN.root"  ),
        "TT_hdampUP"     : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_hdampUP.root"    ),
        "TT_TuneCP5down" : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_TuneCP5down.root"),
        "TT_TuneCP5up"   : ROOT.TFile.Open(args.path + "/" + args.year + "_TT_TuneCP5up.root"  ),
        "TT_JECdown"     : ROOT.TFile.Open(args.path + "/" + args.year + "_TT.root"            ),
        "TT_JECup"       : ROOT.TFile.Open(args.path + "/" + args.year + "_TT.root"            ),
        "TT_JERdown"     : ROOT.TFile.Open(args.path + "/" + args.year + "_TT.root"            ),
        "TT_JERup"       : ROOT.TFile.Open(args.path + "/" + args.year + "_TT.root"            ),
        #"NonTT"          : ROOT.TFile.Open(args.path + "/" + args.year + "_Non_TT.root"        ),
        "QCD"            : ROOT.TFile.Open(args.path + "/" + args.year + "_QCD.root"           ),
        "TTX"            : ROOT.TFile.Open(args.path + "/" + args.year + "_TTX.root"           ),
        "BG_OTHER"       : ROOT.TFile.Open(args.path + "/" + args.year + "_BG_OTHER.root"      ),
        "Data"           : ROOT.TFile.Open(args.path + "/" + args.year + "_Data.root"          ),
    }
        
    Sig = args.sig + "_" + args.mass

    if not args.allMass: 

        if "SYY" in args.sig:
            files[Sig] = ROOT.TFile.Open(args.path + "/" + args.year + "_Stealth%s_2t6j_mStop-%s.root"%(args.sig, args.mass))
        else:
            files[Sig] = ROOT.TFile.Open(args.path + "/" + args.year + "_%s_2t6j_mStop-%s.root"%(args.sig, args.mass))

    else:   
    
        if "SYY" in args.sig:
            files[Sig] = ROOT.TFile.Open(args.path + "/" + args.year + "_Stealth%s_2t6j_mStop-%s.root"%(args.sig, args.mass))
        else:
            files[Sig] = ROOT.TFile.Open(args.path + "/" + args.year + "_%s_2t6j_mStop-%s.root"%(args.sig, args.mass))
        for mass in range(400, 1000, 200):
            if "SYY" in args.sig:
                files["{}_{}".format(args.sig, mass)] = ROOT.TFile.Open(args.path + "/" + args.year + "_Stealth%s_2t6j_mStop-%s.root"%(args.sig, mass))
            else:
                files["{}_{}".format(args.sig, mass)] = ROOT.TFile.Open(args.path + "/" + args.year + "_%s_2t6j_mStop-%s.root"%(args.sig, mass))
        

    print(files)

    # Variations that are the same for all channels
    var_list = ["", "JECup", "JECdown", "JERup", "JERdown", "btgUp", "btgDown", "pdfUp", "pdfDown", "prfUp", "prfDown", "puUp", "puDown", "sclUp", "sclDown", "isrUp", "isrDown", "fsrUp", "fsrDown"]

    var_dict = {
        "0l": var_list + ["ttgUp", "ttgDown", "jetUp", "jetDown"],
        "1l": var_list + ["lepUp", "lepDown"],
        "2l": var_list + ["lepUp", "lepDown"],
    }

    njets = {"0l": range(8, 14), "1l": range(7, 13), "2l": range(6,12)}

    # ---------------------------------------------------------------
    # make regionis list for adding all edges to DoubleDisCo cfg file
    # ---------------------------------------------------------------
    regions = ["ABCD",
               "Val_BD",
               "Val_CD",
               "Val_D", 
    ]

    # ------------------------------------------
    # initialize the dictionaries of any regions
    # ------------------------------------------
    translator = {"ABCD"   : {"A" : "A",  "B" : "B",  "C" : "C",  "D" : "D" },
                  "Val_BD" : {"A" : "b",  "B" : "E",  "C" : "d",  "D" : "F" },
                  "Val_CD" : {"A" : "c",  "B" : "di", "C" : "G",  "D" : "H" },
                  "Val_D"  : {"A" : "dA", "B" : "dB", "C" : "dC", "D" : "dD"},
    }
    ttSys = {"0l": {}, "1l": {}, "2l": {}}
    QCDCRInfo = {"0l": {}, "1l": {}, "2l": {}}
    for d1 in np.arange(args.min, args.max, args.step):
        for d2 in np.arange(args.min, args.max, args.step):
            for channel in ["0l", "1l", "2l"]:
                # ---------------------
                # get the 2D histograms
                # --------------------- 
                histName = "h_DoubleDisCo_%s_disc1_disc2_%s_Njets${NJET}_ABCD"%(args.sig,channel)
                print("Getting TT non-closure sys for ({},{})".format(d1, d2))

                ttSys[channel]["({}, {})".format(d1, d2)], QCDCRInfo[channel]["({:.3f},{:.3f})".format(d1,d2)] = BryansHack(files, channel, Sig, args.mass, histName, regions, translator, d1, d2)

                maxSys = -999.9
                for val in ttSys[channel]["({}, {})".format(d1, d2)].values():
                    if abs(val - 1) > abs(maxSys - 1) or maxSys == -999.9:
                        maxSys = val 
                
                ttSys[channel]["({}, {})".format(d1, d2)]["max"] = maxSys
                    
    print(QCDCRInfo["0l"]["(0.600,0.600)"])
    print(QCDCRInfo["1l"]["(0.600,0.600)"])

    print("Files Loaded")
    all_ABCDEdges = get_edge_info(files, Sig, njets, var_dict, args.step, args.min, args.max, QCDCRInfo)

    make_ABCD_hists(all_ABCDEdges, files, Sig, args.step, args.min, args.max, njets, var_dict, ttSys, args.path, args.output, args.year)
    #write_to_root_file(dc_histos)

if __name__ == "__main__":
    main()    

