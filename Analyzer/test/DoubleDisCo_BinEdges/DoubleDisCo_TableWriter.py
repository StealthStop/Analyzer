# Base class for writing latex tables
# The writeHeader and writeLine methods
# get overridden by the derived classes
class TableWriter:

    def __init__(self, channel, year, tag, model, mass):

        self.f = open("tables/%s/%s_%s_%s_%s_%s.tex" %(channel, year, tag, model, mass, channel), "w")

    def writeHeader(self):
        print(self.f)

    def writeLine(self, **kwargs):
        print(kwargs)

    def writeClose(self):
        self.f.write("\end{tabular}")
        self.f.write("\n")
        self.f.write("}")
        self.f.write("\n")
        self.f.write("\n")
        self.f.close()

# Derived table writer class for the bin edges table
# Over ride the writeHeader and writeLine method
class BinEdgesTable(TableWriter):

    def writeHeader(self):
        self.f.write("\\resizebox{\linewidth}{!}{%")
        self.f.write("\n")
        self.f.write("\\begin{tabular}{| c | c | c | c | c | c | c |}")
        self.f.write("\n")
        self.f.write("\hline")
        self.f.write("\n")
        self.f.write("\\textcolor{ttjetscol}{NJets} &  \multicolumn{2}{c|} {\\textcolor{click}{ABCD}} & \multicolumn{2}{c|} {\\textcolor{rpvcol}{B'D'EF}} & \multicolumn{2}{c|} {\\textcolor{disc}{C'D'GH}} \\\\")
        self.f.write("\n")
        self.f.write("\hline")
        self.f.write("\n")
        self.f.write("& \scriptsize \\textcolor{click}{disc1 edges} & \scriptsize \\textcolor{click}{disc2 edges} & \scriptsize \\textcolor{rpvcol}{disc1 edges} & \scriptsize \\textcolor{click}{disc2 edges} & \scriptsize \\textcolor{click}{disc1 edges} & \scriptsize \\textcolor{disc}{disc2 edges} \\\\")
        self.f.write("\n")
        self.f.write("\hline")
        self.f.write("\n")

    def writeLine(self, **kwargs):
        self.f.write("    \\tiny \\textcolor{ttjetscol}{%s}  & \\tiny \\textcolor{click}{%s} & \\tiny \\textcolor{click}{%s} & \\tiny \\textcolor{rpvcol}{%s} & \\tiny \\textcolor{click}{%s} & \\tiny \\textcolor{click}{%s} & \\tiny \\textcolor{disc}{%s} \\\\" %(kwargs["njet"], kwargs["finalDiscs"]["ABCD"][0], kwargs["finalDiscs"]["ABCD"][1], kwargs["finalDiscs"]["Val_bdEF"][0], kwargs["finalDiscs"]["Val_bdEF"][1], kwargs["finalDiscs"]["Val_cdiGH"][0], kwargs["finalDiscs"]["Val_cdiGH"][1]))
        self.f.write("\n")
        self.f.write("\hline")
        self.f.write("\n")

# Derived table writer class for the ABCD events table
# Over ride the writeHeader and writeLine method
class ABCDeventsTable(TableWriter):

    def writeHeader(self):
        self.f.write("\\resizebox{\linewidth}{!}{%")
        self.f.write("\n")
        self.f.write("\\begin{tabular}{| c | c | c | c | c | c | c |}")
        self.f.write("\n")
        self.f.write("\hline")
        self.f.write("\n")
        self.f.write("\scriptsize \\textcolor{ttjetscol}{NJets} & \scriptsize \\textcolor{click}{sigFracA} & \scriptsize \\textcolor{click}{sigFracB} & \scriptsize \\textcolor{click}{sigFracC} & \scriptsize \\textcolor{click}{sigFracD} & \scriptsize \\textcolor{click}{nBkgEvents(A+C)} & \scriptsize \\textcolor{click}{nBkgEvents(A+B)} \\\\")
        self.f.write("\n")
        self.f.write("\hline")
        self.f.write("\n")

    def writeLine(self, **kwargs):
        self.f.write("\\tiny \\textcolor{ttjetscol}{%s} & \\tiny \\textcolor{click}{%.2f} & \\tiny \\textcolor{click}{%.2f} & \\tiny \\textcolor{click}{%.2f} & \\tiny \\textcolor{click}{%.2f} & \\tiny \\textcolor{click}{%.2f} & \\tiny \\textcolor{click}{%.2f} \\\\" %(kwargs["njet"], kwargs["finalSigFracs"]["A"][0], kwargs["finalSigFracs"]["B"][0], kwargs["finalSigFracs"]["C"][0], kwargs["finalSigFracs"]["D"][0], kwargs["nEvents_AC"], kwargs["nEvents_AB"]))
        self.f.write("\n")
        self.f.write("\hline")
        self.f.write("\n")

# Derived table writer class for the bdEF events table
# Over ride the writeHeader and writeLine method
class BDEFeventsTable(TableWriter):

    def writeHeader(self):
        self.f.write("\\resizebox{\linewidth}{!}{%")
        self.f.write("\n")
        self.f.write("\\begin{tabular}{| c | c | c | c | c | c | c |}")
        self.f.write("\n")
        self.f.write("\hline")
        self.f.write("\n")
        self.f.write("\scriptsize \\textcolor{ttjetscol}{NJets} & \scriptsize \\textcolor{rpvcol}{sigFracB'} & \scriptsize \\textcolor{rpvcol}{sigFracD'} & \scriptsize \\textcolor{rpvcol}{sigFracE} & \scriptsize \\textcolor{rpvcol}{sigFracF} & \scriptsize \\textcolor{rpvcol}{nBkgEvents(B'+D')} & \scriptsize \\textcolor{rpvcol}{nBkgEvents(E+F)} \\\\")
        self.f.write("\n")
        self.f.write("\hline")
        self.f.write("\n")

    def writeLine(self, **kwargs):
        self.f.write("\\tiny \\textcolor{ttjetscol}{%s} & \\tiny \\textcolor{rpvcol}{%.2f} & \\tiny \\textcolor{rpvcol}{%.2f} & \\tiny \\textcolor{rpvcol}{%.2f} & \\tiny \\textcolor{rpvcol}{%.2f} & \\tiny \\textcolor{rpvcol}{%.2f} & \\tiny \\textcolor{rpvcol}{%.2f} \\\\" %(kwargs["njet"], kwargs["finalSigFracs"]["A"][0], kwargs["finalSigFracs"]["B"][0], kwargs["finalSigFracs"]["C"][0], kwargs["finalSigFracs"]["D"][0], kwargs["nEvents_AC"], kwargs["nEvents_AB"]))
        self.f.write("\n")
        self.f.write("\hline")
        self.f.write("\n")

# Derived table writer class for the cdiGH events
# Over ride the writeHeader and writeLine method
class CDGHeventsTable(TableWriter):

    def writeHeader(self):
        self.f.write("\\resizebox{\linewidth}{!}{%")
        self.f.write("\n")
        self.f.write("\\begin{tabular}{| c | c | c | c | c | c | c |}")
        self.f.write("\n")
        self.f.write("\hline")
        self.f.write("\n")
        self.f.write("\scriptsize \\textcolor{ttjetscol}{NJets} & \scriptsize \\textcolor{disc}{sigFracC'} & \scriptsize \\textcolor{disc}{sigFracD'} & \scriptsize \\textcolor{disc}{sigFracG} & \scriptsize \\textcolor{disc}{sigFracH} & \scriptsize \\textcolor{disc}{nBkgEvents(C'+D')} & \scriptsize \\textcolor{disc}{nBkgEvents(G+H)} \\\\")
        self.f.write("\n")
        self.f.write("\hline")
        self.f.write("\n") 

    def writeLine(self, **kwargs):
        self.f.write("\\tiny \\textcolor{ttjetscol}{%s} & \\tiny \\textcolor{disc}{%.2f} & \\tiny \\textcolor{disc}{%.2f} & \\tiny \\textcolor{disc}{%.2f} & \\tiny \\textcolor{disc}{%.2f} & \\tiny \\textcolor{disc}{%.2f} & \\tiny \\textcolor{disc}{%.2f} \\\\" %(kwargs["njet"], kwargs["finalSigFracs"]["A"][0], kwargs["finalSigFracs"]["B"][0], kwargs["finalSigFracs"]["C"][0], kwargs["finalSigFracs"]["D"][0], kwargs["nEvents_AC"], kwargs["nEvents_AB"]))
        self.f.write("\n")
        self.f.write("\hline")
        self.f.write("\n")

# Derived table writer class for the cdiGH events
# Over ride the writeHeader and writeLine method
class SubDivDeventsTable(TableWriter):

    def writeHeader(self):
        d.write("\\resizebox{\linewidth}{!}{%")
        d.write("\n")
        d.write("\def\\arraystretch{0.6}")
        d.write("\n")
        d.write("\\begin{tabular}{| c | c | c | c | c | c | c | c | c | c | c |}")
        d.write("\n")
        d.write("\hline")
        d.write("\n")
        d.write("\\textcolor{ttjetscol}{NJets} &  \multicolumn{2}{c|} {\\textcolor{click}{Edges}} & \multicolumn{4}{c|} {\\textcolor{massReg}{SigFracs}} & \multicolumn{4}{c|} {\\textcolor{disc}{nEvents(Sig+Bkg)}} \\\\")
        d.write("\n")
        d.write("\hline")
        d.write("\n")
        d.write("& \scriptsize \\textcolor{click}{disc1} & \scriptsize \\textcolor{click}{disc2} & \scriptsize \\textcolor{massReg}{sigFrac\_dA} & \scriptsize \\textcolor{massReg}{sigFrac\_dB} & \scriptsize \\textcolor{massReg}{SigFrac\_dC} & \scriptsize \\textcolor{massReg}{SigFrac\_dD} & \scriptsize \\textcolor{disc}{in region dA} & \scriptsize \\textcolor{disc}{in region dB} & \scriptsize \\textcolor{disc}{in region dC} & \scriptsize \\textcolor{disc}{in region dD} \\\\")
        d.write("\n")
        d.write("\hline")
        d.write("\n")

    def writeLine(self, **kwargs):
        self.f.write("\\tiny \\textcolor{ttjetscol}{%s} &  \\tiny \\textcolor{click}{%s} & \\tiny \\textcolor{click}{%s} & \\tiny \\textcolor{massReg}{%.3f} & \\tiny \\textcolor{massReg}{%.3f} & \\tiny \\textcolor{massReg}{%.3f} & \\tiny \\textcolor{massReg}{%.3f} & \\tiny \\textcolor{disc}{%.2f} & \\tiny \\textcolor{disc}{%.2f} & \\tiny \\textcolor{disc}{%.2f} & \\tiny \\textcolor{disc}{%.2f} \\\\" %(kwargs["njet"], kwargs["finalDiscs"][0], kwargs["finalDiscs"][1], kwargs["finalSigFracs"]["A"][0], kwargs["finalSigFracs"]["B"][0], kwargs["finalSigFracs"]["C"][0], kwargs["finalSigFracs"]["D"][0], kwargs["final_nTot_Sig"]["A"][0]+kwargs["final_nTot_Bkg"]["A"][0], kwargs["final_nTot_Sig"]["B"][0]+kwargs["final_nTot_Bkg"]["B"][0], kwargs["final_nTot_Sig"]["C"][0]+kwargs["final_nTot_Bkg"]["C"][0], kwargs["final_nTot_Sig"]["D"][0]+kwargs["final_nTot_Bkg"]["D"][0]) )
        self.f.write("\n")
        self.f.write("\hline")
        self.f.write("\n")
