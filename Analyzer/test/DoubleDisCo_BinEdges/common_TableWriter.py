# -------------------------------------
# Base class for writing latex tables
# The writeHeader and writeLine methods
# get overridden by the derived classes
# -------------------------------------
class TableWriter:

    def __init__(self, tablesPath, channel, year, tag, model):

        self.f = open("%s/%s_%s_%s_%s.tex" %(tablesPath, year, tag, model, channel), "w")
        self.writeHeader()

    def writeHeader(self):
        print(self.f)

    def writeLine(self, **kwargs):
        print(kwargs)

    def writeClose(self):
        self.f.write("    \end{tabular}")
        self.f.write("\n")
        self.f.write("}")
        self.f.write("\n")
        self.f.write("\n")
        self.f.close()


# ----------------------------------------------------
# Derived table writer class for the ABCD fracs table
# Over ride the writeHeader and writeLine method
# ----------------------------------------------------
class ABCDfracsTable(TableWriter):

    def writeHeader(self):
        self.f.write("\\resizebox{\linewidth}{!}{%")
        self.f.write("\n")
        self.f.write("    \def\\arraystretch{0.6}")
        self.f.write("\n")
        self.f.write("    \\begin{tabular}{| c | c | c | c | c | c | c | c | c |}")
        self.f.write("\n")
        self.f.write("        \hline")
        self.f.write("\n")
        self.f.write("        \\textbf{\\njets} & \\multicolumn{2}{c|}{\\textbf{In Region A}} & \\multicolumn{2}{c|}{\\textbf{In Region B}} & \\multicolumn{2}{c|}{\\textbf{In Region C}} & \\multicolumn{2}{c|}{\\textbf{In Region D}} \\\\")
        self.f.write("\n")
        self.f.write("        \hline")
        self.f.write("\n")
        self.f.write("                          & \\textbf{\\ttjets Frac.} & \\textbf{Sig. Frac.} & \\textbf{\\ttjets Frac.} & \\textbf{Sig. Frac.} & \\textbf{\\ttjets Frac.} & \\textbf{Sig. Frac.} & \\textbf{\\ttjets Frac.} & \\textbf{Sig. Frac.} \\\\")
        self.f.write("\n")
        self.f.write("        \hline")
        self.f.write("\n")

    def writeLine(self, **kwargs):
        self.f.write("        %s & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f\\\\" %(kwargs["njet"], kwargs["finalTTfracs"]["A"][0], kwargs["finalSigFracs"]["A"][0], kwargs["finalTTfracs"]["B"][0], kwargs["finalSigFracs"]["B"][0], kwargs["finalTTfracs"]["C"][0], kwargs["finalSigFracs"]["C"][0], kwargs["finalTTfracs"]["D"][0], kwargs["finalSigFracs"]["D"][0]))
        self.f.write("\n")
        self.f.write("        \hline")
        self.f.write("\n")

# ----------------------------------------------------
# Derived table writer class for the Val fracs table
# Over ride the writeHeader and writeLine method
# ----------------------------------------------------
class ValFracsTable(TableWriter):

    def writeHeader(self):
        self.f.write("\\resizebox{\linewidth}{!}{%")
        self.f.write("\n")
        self.f.write("    \def\\arraystretch{0.6}")
        self.f.write("\n")
        self.f.write("    \\begin{tabular}{| c | c | c | c | c | c | c | c | c | c | c | c | c |}")
        self.f.write("\n")
        self.f.write("        \hline")
        self.f.write("\n")
        self.f.write("        \\textbf{\\njets} & \\multicolumn{3}{c|}{\\textbf{In Region A'}} & \\multicolumn{3}{c|}{\\textbf{In Region B'}} & \\multicolumn{3}{c|}{\\textbf{In Region C'}} & \\multicolumn{3}{c|}{\\textbf{In Region D'}} \\\\")
        self.f.write("\n")
        self.f.write("        \hline")
        self.f.write("\n")
        self.f.write("                          & \\textbf{Val I} & \\textbf{Val II} & \\textbf{Val III} & \\textbf{Val I} & \\textbf{Val II} & \\textbf{Val III} & \\textbf{Val I} & \\textbf{Val II} & \\textbf{Val III} & \\textbf{Val I} & \\textbf{Val II} & \\textbf{Val III} \\\\")
        self.f.write("\n")
        self.f.write("        \hline")
        self.f.write("\n")

    def writeLine(self, **kwargs):
        self.f.write("        %s & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f\\\\" %(kwargs["njet"], kwargs["finalSigFracs"]["Val_BD"]["A"][0], kwargs["finalSigFracs"]["Val_CD"]["A"][0], kwargs["finalSigFracs"]["Val_D"]["A"][0], kwargs["finalSigFracs"]["Val_BD"]["B"][0], kwargs["finalSigFracs"]["Val_CD"]["B"][0], kwargs["finalSigFracs"]["Val_D"]["B"][0], kwargs["finalSigFracs"]["Val_BD"]["C"][0], kwargs["finalSigFracs"]["Val_CD"]["C"][0], kwargs["finalSigFracs"]["Val_D"]["C"][0], kwargs["finalSigFracs"]["Val_BD"]["D"][0], kwargs["finalSigFracs"]["Val_CD"]["D"][0], kwargs["finalSigFracs"]["Val_D"]["D"][0]))
        self.f.write("\n")
        self.f.write("        \hline")
        self.f.write("\n")


# --------------------------------------------------------------------
# Derived table writer class for the ABCD and validation regions table
# Over ride the writeHeader and writeLine method
# --------------------------------------------------------------------
class SignalFractionsAllRegionsTable(TableWriter):

    def __init__(self, tablesPath, channel, year, tag, model, mass):

        self.f = open("%s/%s_%s_%s_%s_%s.tex" %(tablesPath, year, tag, model, mass, channel), "w")
        self.writeHeader()

        self.startedABCD = False
        self.startedBDEF = False
        self.startedCDGH = False
        self.startedSubD = False

        self.ABCDstr = ""
        self.BDEFstr = ""
        self.CDGHstr = ""
        self.SubDstr = ""

    def writeHeader(self):
        self.f.write("\\resizebox{0.8\linewidth}{!}{%")
        self.f.write("\n")
        self.f.write("    \def\\arraystretch{0.6}")
        self.f.write("\n")
        self.f.write("    \\begin{tabular}{| c | c | c | c | c |}")
        self.f.write("\n")
        self.f.write("        \hline")
        self.f.write("\n")

    def writeLine(self, **kwargs):

        Njet = kwargs["njet"]
        if "incl" in Njet:
            Njet = "$\geq%s$"%(Njet.replace("incl", ""))

        region = kwargs["region"]

        if region == "ABCD":
            if not self.startedABCD:
                self.ABCDstr += "        \scriptsize \\textcolor{ttjetscol}{NJets} & \scriptsize \\textcolor{click}{sigFracA} & \scriptsize \\textcolor{click}{sigFracB} & \scriptsize \\textcolor{click}{sigFracC} & \scriptsize \\textcolor{click}{sigFracD} \\\\"
                self.startedABCD = True

            self.ABCDstr += "\n"
            self.ABCDstr += "        \hline"
            self.ABCDstr += "\n"
            self.ABCDstr += "        \\tiny \\textcolor{ttjetscol}{%s} & \\tiny \\textcolor{click}{%.3f} & \\tiny \\textcolor{click}{%.3f} & \\tiny \\textcolor{click}{%.3f} & \\tiny \\textcolor{click}{%.3f} \\\\" %(Njet, kwargs["finalSigFracs"]["A"][0], kwargs["finalSigFracs"]["B"][0], kwargs["finalSigFracs"]["C"][0], kwargs["finalSigFracs"]["D"][0])
            self.ABCDstr += "\n"
            self.ABCDstr += "        \hline"
            self.ABCDstr += "\n"

        elif region == "Val_BD":
            if not self.startedBDEF: 
                self.BDEFstr += "        \scriptsize \\textcolor{ttjetscol}{NJets} & \scriptsize \\textcolor{rpvcol}{sigFracB'} & \scriptsize \\textcolor{rpvcol}{sigFracD'} & \scriptsize \\textcolor{rpvcol}{sigFracE} & \scriptsize \\textcolor{rpvcol}{sigFracF} \\\\"
                self.startedBDEF = True

            self.BDEFstr += "\n"
            self.BDEFstr += "        \hline"
            self.BDEFstr += "\n"
            self.BDEFstr += "        \\tiny \\textcolor{ttjetscol}{%s} & \\tiny \\textcolor{rpvcol}{%.3f} & \\tiny \\textcolor{rpvcol}{%.3f} & \\tiny \\textcolor{rpvcol}{%.3f} & \\tiny \\textcolor{rpvcol}{%.3f} \\\\" %(Njet, kwargs["finalSigFracs"]["A"][0], kwargs["finalSigFracs"]["B"][0], kwargs["finalSigFracs"]["C"][0], kwargs["finalSigFracs"]["D"][0])
            self.BDEFstr += "\n"
            self.BDEFstr += "        \hline"
            self.BDEFstr += "\n"

        elif region == "Val_CD":
            if not self.startedCDGH:
                self.CDGHstr += "        \scriptsize \\textcolor{ttjetscol}{NJets} & \scriptsize \\textcolor{disc}{sigFracC'} & \scriptsize \\textcolor{disc}{sigFracD'} & \scriptsize \\textcolor{disc}{sigFracG} & \scriptsize \\textcolor{disc}{sigFracH} \\\\"
                self.startedCDGH = True

            self.CDGHstr += "\n"
            self.CDGHstr += "        \hline"
            self.CDGHstr += "\n"
            self.CDGHstr += "        \\tiny \\textcolor{ttjetscol}{%s} & \\tiny \\textcolor{disc}{%.3f} & \\tiny \\textcolor{disc}{%.3f} & \\tiny \\textcolor{disc}{%.3f} & \\tiny \\textcolor{disc}{%.3f} \\\\" %(Njet, kwargs["finalSigFracs"]["A"][0], kwargs["finalSigFracs"]["B"][0], kwargs["finalSigFracs"]["C"][0], kwargs["finalSigFracs"]["D"][0])
            self.CDGHstr += "\n"
            self.CDGHstr += "        \hline"
            self.CDGHstr += "\n"

        elif region == "Val_D":
            if not self.startedSubD:
                self.SubDstr += "        \\scriptsize \\textcolor{ttjetscol}{NJets} & \\scriptsize \\textcolor{subDivD}{sigFrac\\_dA} & \\scriptsize \\textcolor{subDivD}{sigFrac\\_dB} & \\scriptsize \\textcolor{subDivD}{SigFrac\\_dC} & \\scriptsize \\textcolor{subDivD}{SigFrag\\_dD} \\\\"
                self.startedSubD = True

            self.SubDstr += "\n"
            self.SubDstr += "        \hline"
            self.SubDstr += "\n"
            self.SubDstr += "        \\tiny \\textcolor{ttjetscol}{%s} & \\tiny \\textcolor{subDivD}{%.3f} & \\tiny \\textcolor{subDivD}{%.3f} & \\tiny \\textcolor{subDivD}{%.3f} & \\tiny \\textcolor{subDivD}{%.3f} \\\\" %(Njet, kwargs["finalSigFracs"]["A"][0], kwargs["finalSigFracs"]["B"][0], kwargs["finalSigFracs"]["C"][0], kwargs["finalSigFracs"]["D"][0])
            self.SubDstr += "\n"
            self.SubDstr += "        \hline"
            self.SubDstr += "\n"

    def writeClose(self):
        self.f.write(self.ABCDstr)
        self.f.write(self.BDEFstr)
        self.f.write(self.CDGHstr)
        self.f.write(self.SubDstr)
        self.f.write("    \end{tabular}")
        self.f.write("\n")
        self.f.write("}")
        self.f.write("\n")
        self.f.write("\n")
        self.f.close()

# -------------------------------------------------------
# Derived table writer class for the optimized ABCD edges
# Over ride the writeHeader and writeLine method
# -------------------------------------------------------
class Optimized_ABCDedges(TableWriter):
   
    def writeHeader(self):
        self.f.write("\\resizebox{\linewidth}{!}{%")
        self.f.write("\n")
        self.f.write("    \def\\arraystretch{0.6}")
        self.f.write("\n")
        self.f.write("    \\begin{tabular}{| c | c | c | c | c | c | c | c |}")
        self.f.write("\n")
        self.f.write("        \hline")
        self.f.write("\n")
        self.f.write("        \\textbf{ABCD edges} & \\textbf{Significance} & \\textbf{\\njets} & \\textbf{Non-closure} & \\textbf{Pull} & \\textbf{Sig. Fracs B} & \\textbf{Sig. Fracs C} & \\textbf{Sig. Fracs D} \\\\")
        self.f.write("\n")
        self.f.write("        \hline")
        self.f.write("\n")

    def writeLine(self, **kwargs):
        tempDict   = kwargs["Best_Parameters"]
        njets      = kwargs["Njets"]
        maxOptions = kwargs["numEdgeChoices"]

        sortedKeys = sorted(tempDict.keys(), reverse=True)
        
        iOption = 0
        maxOptions = 4
        for key in sortedKeys:

            if iOption > maxOptions:
                break

            value = tempDict[key]

            firstWrite = False
            for njet in njets:
                Njet = njet.replace("incl", "")
                if "incl" in njet:
                    Njet = "$\geq%s$"%(Njet)
                if not firstWrite:
                    self.f.write("        \\multirow{%d}{*}{%s} & \\multirow{%d}{*}{%.3f} & %s & %.3f & %.3f & %.3f & %.3f & %.3f \\\\"%(len(njets), value["ABCDedges"], len(njets), value["Significance"], Njet, value["nonClosure_njet%s"%njet], value["nonCllosurePull_njet%s"%njet], value["sigFracB_njet%s"%njet], value["sigFracC_njet%s"%njet], value["sigFracD_njet%s"%njet]))
                    firstWrite = True
                else:
                    self.f.write("                                            &                        & %s & %.3f & %.3f & %.3f & %.3f & %.3f \\\\"%(Njet, value["nonClosure_njet%s"%njet], value["nonCllosurePull_njet%s"%njet], value["sigFracB_njet%s"%njet], value["sigFracC_njet%s"%njet], value["sigFracD_njet%s"%njet]))

                self.f.write("\n")

            self.f.write("        \hline")
            self.f.write("\n")

            iOption += 1

# ----------------------------------------------------------------------------------
# Derived table writer class for the tt sys with the maxi value of MC corrected data
# Over ride the writeHeader and writeLine method
# ----------------------------------------------------------------------------------
class maximumCorrectedData_ttSyst(TableWriter):

    def writeHeader(self):
        self.f.write("\\centering")
        self.f.write("\n")
        self.f.write("\\resizebox{0.6\linewidth}{!}{%")
        self.f.write("\n")
        self.f.write("    \def\\arraystretch{0.7}")
        self.f.write("\n")
        self.f.write("    \\begin{tabular}{| c | c | c |}")
        self.f.write("\n")
        self.f.write("        \hline")
        self.f.write("\n")
        self.f.write("        \\textbf{\\njets} & \\textbf{maxi Corrected Data Closure} & \\textbf{\\ttbar syst.} \\\\")
        self.f.write("\n")
        self.f.write("        \hline")
        self.f.write("\n")

    def writeLine(self, **kwargs):

        njetToWrite = kwargs["njet"]
        if "incl" in njetToWrite:
            njetToWrite = "$\\geq%s$"%(njetToWrite.partition("incl")[0])
        
        self.f.write("        %s & %s & %.3f \\\\" %(njetToWrite, kwargs["maxCorrData"], kwargs["ttSyst"]))
        self.f.write("\n")
        self.f.write("        \hline")
        self.f.write("\n")


 
