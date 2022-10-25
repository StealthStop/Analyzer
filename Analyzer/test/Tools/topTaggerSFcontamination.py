#! /usr/bin/env python

import ROOT, random, os, argparse, string, copy, math
ROOT.PyConfig.IgnoreCommandLineOptions = True

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat("")
ROOT.gStyle.SetPaintTextFormat("3.2f")
ROOT.gStyle.SetFrameLineWidth(2)
ROOT.gStyle.SetEndErrorSize(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH2.SetDefaultSumw2()

path = "/uscms/home/jhiltb/nobackup/susy/ZeroAndTwoLep/CMSSW_10_2_9/src/Analyzer/Analyzer/test/condor/2016_TopTaggerSF"

ttFile  = ROOT.TFile.Open(path+"/2016_TT.root")
qcdFile  = ROOT.TFile.Open(path+"/2016_QCD.root")
bgFile  = ROOT.TFile.Open(path+"/2016_BG.root")

sgFile1 = ROOT.TFile.Open(path+"/2016_RPV_2t6j_mStop-350.root")
sgFile2 = ROOT.TFile.Open(path+"/2016_RPV_2t6j_mStop-550.root")
sgFile3 = ROOT.TFile.Open(path+"/2016_RPV_2t6j_mStop-850.root")
sgFile4 = ROOT.TFile.Open(path+"/2016_RPV_2t6j_mStop-1150.root")

selections = ["tt", "qcd"]

Njets = ["7", "8", "9", "10", "11", "12incl"]

histo = "HT"

counts = {}

def makeTable(name, selection="tt"):

    procStr = r"\ttjets"
    if selection == "qcd":
        procStr = "QCD multi-jet"

    table = open("%s.tex"%(name), "w")
    table.write(r"\def\arraystretch{1.6}")
    table.write("\n")
    table.write(r"    \begin{tabular}{|c|c|c|c|c|c|}")
    table.write("\n")
    table.write(r"        \hline")
    table.write("\n")
    table.write(r"        \mboldcol{UMNMaroon}{\njets} & \mboldcol{UMNMaroon}{%s} & \mboldcol{UMNMaroon}{RPV $\mstop=350$ GeV} & \mboldcol{UMNMaroon}{RPV $\mstop=550$ GeV} & \mboldcol{UMNMaroon}{RPV $\mstop=850$ GeV} & \mboldcol{UMNMaroon}{RPV $\mstop=1150$ GeV}\\\hline"%(procStr))
    table.write("\n")

    return table

def calcRatioUnc(A, B, Aunc=-1.0, Bunc=-1.0, mode=1):

    if Aunc == -1.0:
        Aunc = math.sqrt(A)
    if Bunc == -1.0:
        Bunc = math.sqrt(B)

    if mode == 1:
        if A != 0.0 and B != 0.0:
            return math.sqrt((Aunc/B)**2.0 + (A*Bunc/B**2.0)**2.0)
        else:
            return -1.0 
    elif mode == 2:
        if A != 0.0 and B != 0.0:
            return math.sqrt((B*Aunc**2.0/(A + B)**2.0)**2.0 + (A*Bunc/(A + B)**2.0)**2.0)
        else:
            return -1.0

for selection in selections:
    for Njet in Njets:
    
        h_tt  = ttFile.Get("Njet"+Njet+"_"+selection+"/"+histo+"_Njet"+Njet+"_"+selection)
        h_qcd = qcdFile.Get("Njet"+Njet+"_"+selection+"/"+histo+"_Njet"+Njet+"_"+selection)
        h_bg  = bgFile.Get("Njet"+Njet+"_"+selection+"/"+histo+"_Njet"+Njet+"_"+selection)
        h_sg1 = sgFile1.Get("Njet"+Njet+"_"+selection+"/"+histo+"_Njet"+Njet+"_"+selection)
        h_sg2 = sgFile2.Get("Njet"+Njet+"_"+selection+"/"+histo+"_Njet"+Njet+"_"+selection)
        h_sg3 = sgFile3.Get("Njet"+Njet+"_"+selection+"/"+histo+"_Njet"+Njet+"_"+selection)
        h_sg4 = sgFile4.Get("Njet"+Njet+"_"+selection+"/"+histo+"_Njet"+Njet+"_"+selection)

        if selection not in counts:
            counts[selection] = {}

        if Njet not in counts[selection]:
            counts[selection][Njet] = {}

        ttUnc  = ROOT.Double(0.0)
        qcdUnc = ROOT.Double(0.0)
        bgUnc  = ROOT.Double(0.0)
        sg1Unc = ROOT.Double(0.0)
        sg2Unc = ROOT.Double(0.0)
        sg3Unc = ROOT.Double(0.0)
        sg4Unc = ROOT.Double(0.0)

        counts[selection][Njet]["TT"]  = h_tt.IntegralAndError(0, -1, ttUnc)
        counts[selection][Njet]["BG"]  = h_bg.IntegralAndError(0, -1, bgUnc)
        counts[selection][Njet]["QCD"] = h_qcd.IntegralAndError(0, -1, qcdUnc)
        counts[selection][Njet]["SG1"] = h_sg1.IntegralAndError(0, -1, sg1Unc)
        counts[selection][Njet]["SG2"] = h_sg2.IntegralAndError(0, -1, sg2Unc)
        counts[selection][Njet]["SG3"] = h_sg3.IntegralAndError(0, -1, sg3Unc)
        counts[selection][Njet]["SG4"] = h_sg4.IntegralAndError(0, -1, sg4Unc)

        counts[selection][Njet]["TTunc"]  = ttUnc
        counts[selection][Njet]["BGunc"]  = bgUnc
        counts[selection][Njet]["QCDunc"] = qcdUnc
        counts[selection][Njet]["SG1unc"] = sg1Unc
        counts[selection][Njet]["SG2unc"] = sg2Unc
        counts[selection][Njet]["SG3unc"] = sg3Unc
        counts[selection][Njet]["SG4unc"] = sg4Unc

tableTT  = makeTable("topTaggerSFcontamination_TT")
tableQCD = makeTable("topTaggerSFcontamination_QCD", "qcd")
tableTT_evts = makeTable("topTaggerSFcontamination_TT_evts")
tableQCD_evts = makeTable("topTaggerSFcontamination_QCD_evts", "qcd")

for selection, dict1 in counts.items():
    for Njet, dict2 in dict1.items():
    
        tt  = dict2["TT"]
        qcd = dict2["QCD"]
        bg  = dict2["BG"]
        sg1 = dict2["SG1"]
        sg2 = dict2["SG2"]
        sg3 = dict2["SG3"]
        sg4 = dict2["SG4"]

        ttUnc  = dict2["TTunc"]
        qcdUnc = dict2["QCDunc"]
        bgUnc  = dict2["BGunc"]

        sg1Unc = dict2["SG1unc"]
        sg2Unc = dict2["SG2unc"]
        sg3Unc = dict2["SG3unc"]
        sg4Unc = dict2["SG4unc"]

        ttFrac  = "%.3f"%(tt / bg)
        qcdFrac = "%.3f"%(qcd / bg)
        sg1Frac = "%.3f"%(sg1 / (sg1 + bg))
        sg2Frac = "%.3f"%(sg2 / (sg2 + bg))
        sg3Frac = "%.3f"%(sg3 / (sg3 + bg))
        sg4Frac = "%.3f"%(sg4 / (sg4 + bg))

        ttFracUnc  = ("%.3f"%(calcRatioUnc(tt,  bg, ttUnc, bgUnc))).replace("-1.000", r"\text{N/A}")
        qcdFracUnc = ("%.3f"%(calcRatioUnc(qcd, bg, ttUnc, bgUnc))).replace("-1.000", r"\text{N/A}")
        sg1FracUnc = ("%.3f"%(calcRatioUnc(sg1, bg, sg1Unc, bgUnc, 2))).replace("-1.000", r"\text{N/A}")
        sg2FracUnc = ("%.3f"%(calcRatioUnc(sg2, bg, sg2Unc, bgUnc, 2))).replace("-1.000", r"\text{N/A}")
        sg3FracUnc = ("%.3f"%(calcRatioUnc(sg3, bg, sg3Unc, bgUnc, 2))).replace("-1.000", r"\text{N/A}")
        sg4FracUnc = ("%.3f"%(calcRatioUnc(sg4, bg, sg4Unc, bgUnc, 2))).replace("-1.000", r"\text{N/A}")

        if selection == "tt":
            if "incl" not in Njet:
                tableTT.write(r"        \mboldcol{UMNMaroon}{$\njets=%s$}    & $%s \pm %s$ & $%s \pm %s$ & $%s \pm %s$ & $%s \pm %s$ & $%s \pm %s$\\\hline"%(Njet, ttFrac, ttFracUnc, sg1Frac, sg1FracUnc, sg2Frac, sg2FracUnc, sg3Frac, sg3FracUnc, sg4Frac, sg4FracUnc))
                tableTT_evts.write(r"        \mboldcol{UMNMaroon}{$\njets=%s$}    & $%.1f$ & $%.1f$ & $%.1f$ & $%.1f$ & $%.1f$\\\hline"%(Njet, tt, sg1, sg2, sg3, sg4))
            else:
                tableTT.write(r"        \mboldcol{UMNMaroon}{$\njets\geq%s$} & $%s \pm %s$ & $%s \pm %s$ & $%s \pm %s$ & $%s \pm %s$ & $%s \pm %s$\\\hline"%(Njet.replace("incl", ""), ttFrac, ttFracUnc, sg1Frac, sg1FracUnc, sg2Frac, sg2FracUnc, sg3Frac, sg3FracUnc, sg4Frac, sg4FracUnc))
                tableTT_evts.write(r"        \mboldcol{UMNMaroon}{$\njets\geq%s$}    & $%.1f$ & $%.1f$ & $%.1f$ & $%.1f$ & $%.1f$\\\hline"%(Njet.replace("incl", ""), tt, sg1, sg2, sg3, sg4))
              
            tableTT.write("\n")
            tableTT_evts.write("\n")

        else:
            if "incl" not in Njet:
                tableQCD.write(r"        \mboldcol{UMNMaroon}{$\njets=%s$}    & $%s \pm %s$ & $%s \pm %s$ & $%s \pm %s$ & $%s \pm %s$ & $%s \pm %s$\\\hline"%(Njet, qcdFrac, qcdFracUnc, sg1Frac, sg1FracUnc, sg2Frac, sg2FracUnc, sg3Frac, sg3FracUnc, sg4Frac, sg4FracUnc))
                tableQCD_evts.write(r"        \mboldcol{UMNMaroon}{$\njets=%s$}    & $%.1f$ & $%.1f$ & $%.1f$ & $%.1f$ & $%.1f$\\\hline"%(Njet, qcd, sg1, sg2, sg3, sg4))
            else:
                tableQCD.write(r"        \mboldcol{UMNMaroon}{$\njets\geq%s$} & $%s \pm %s$ & $%s \pm %s$ & $%s \pm %s$ & $%s \pm %s$ & $%s \pm %s$\\\hline"%(Njet.replace("incl", ""), qcdFrac, qcdFracUnc, sg1Frac, sg1FracUnc, sg2Frac, sg2FracUnc, sg3Frac, sg3FracUnc, sg4Frac, sg4FracUnc))
                tableQCD_evts.write(r"        \mboldcol{UMNMaroon}{$\njets\geq%s$}    & $%.1f$ & $%.1f$ & $%.1f$ & $%.1f$ & $%.1f$\\\hline"%(Njet.replace("incl", ""), qcd, sg1, sg2, sg3, sg4))

            tableQCD.write("\n")
            tableQCD_evts.write("\n")


tableTT.write(r"    \end{tabular}")
tableTT_evts.write(r"    \end{tabular}")
tableQCD.write(r"    \end{tabular}")
tableQCD_evts.write(r"    \end{tabular}")

tableTT.close()
tableTT_evts.close()
tableQCD.close()
tableQCD_evts.close()
