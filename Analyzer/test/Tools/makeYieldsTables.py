import os
import re
import copy
import ROOT
import argparse

# Sort a list of strings and include natural sorting e.g. 1, 2, 3, ..., 11, 12...
def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)

class Yields:

    def __init__(self, outputDir, njets, channel):

        # Maps a sample, Njet bin to a number of weighted event counts
        self.data = {}

        # Map each sample to the total number of weighted event counts (for sorting in the table)
        # Keep track separately for signal and background
        self.bkgSamples = {}
        self.sigSamples = {}

        # List of njet bins to loop over
        self.njets   = njets

        self.channel = channel

        self.outputDir = outputDir

        # Amount of spaces in LaTeX table to square everything up
        self.njetFlush   = 33 
        self.sampleFlush = 35

        # A histogram that is the summation of all background processes
        # Useful for calculating percentage yields
        self.totalHisto = None

        self.niceNames = {"TT"              : "$t\\bar{t}$ + jets",    
                          "QCD"             : "QCD multijet",      
                          "WJets"           : "W + jets",          
                          "DYJetsToLL_M-50" : "DY + jets",         
                          "Diboson"         : "Diboson",           
                          "Triboson"        : "Triboson",          
                          "TTZToLLNuNu_M-10": "$t\\bar{t}Z(\\ell\\nu\\ell\\nu)$ + jets",          
                          "TTZToQQ"         : "$t\\bar{t}Z(qq)$ + jets",          
                          "ttHJetToNonbb"   : "$t\\bar{t}H(\\text{non-}bb)$ + jets",   
                          "ttHJetTobb"      : "$t\\bar{t}H(bb)$ + jets",   
                          "TTWJetsToLNu"    : "$t\\bar{t}W(\\ell\\nu)$ + jets",   
                          "TTWJetsToQQ"     : "$t\\bar{t}W(qq)$ + jets",   
                          "TTTT"            : "$t\\bar{t}t\\bar{t}$",   
                          "TTWW"            : "$t\\bar{t} + WW$",      
                          "TTWZ"            : "$t\\bar{t} + WZ$",      
                          "TTTW"            : "$t\\bar{t} + tW$",      
                          "TTZH"            : "$t\\bar{t} + ZH$",      
                          "TTWH"            : "$t\\bar{t} + WH$",      
                          "TTZZ"            : "$t\\bar{t} + ZZ$",      
                          "TTHH"            : "$t\\bar{t} + HH$",      
                          "TTTJ"            : "$t\\bar{t} + t$ + jets"
        }


    # Loop over an Njets histogram and store per bin contents
    # in the self.data dictionary
    def processHisto(self, sample, inclusiveBin, histo=None):

        if histo == None:
            histo = self.totalHisto

        # Store the square of the error to making adding simple
        # Just take a square root in the end when using the error
        sampleTotal       = 0.0
        errorSquaredTotal = 0.0
        for bin in range(1, histo.GetNbinsX()+1):

            content = histo.GetBinContent(bin)
            error   = histo.GetBinError(bin)

            sampleTotal       += content
            errorSquaredTotal += error**2.0

            # Want to make the Njet = 11/12/13 bin inclusive (depending on channel)
            # So call any Njet > 11/12/13 as "11/12/13", respectively so that they get added together
            njetBin = bin - 1
            if njetBin >= inclusiveBin:
                njetBin = inclusiveBin

            key = sample + "_Njet" + str(njetBin)

            if key not in self.data:
                self.data[key] = [content, errorSquaredTotal]
            else:
                self.data[key][0] += content
                self.data[key][1] += error**2.0

        # Map each sample to the total number of weighted event counts
        # This is used to sort the samples by importance later on
        if "mStop" in sample:
            if sample not in self.sigSamples:
                self.sigSamples[sample]  = [sampleTotal, errorSquaredTotal]
            else:
                self.sigSamples[sample][0] += sampleTotal
                self.sigSamples[sample][1] += errorSquaredTotal
        else:
            if sample not in self.bkgSamples:
                self.bkgSamples[sample]  = [sampleTotal, errorSquaredTotal]
            else:
                self.bkgSamples[sample][0] += sampleTotal
                self.bkgSamples[sample][1] += errorSquaredTotal

            if sample != "AllBkg":
                # Construct the total histogram only for all backgrounds
                if self.totalHisto != None:
                    self.totalHisto.Add(histo)
                else:
                    self.totalHisto = copy.deepcopy(histo.Clone("aclone"))

    # Write out generic first part of table with header information
    def writeHeader(self):
        rowHeader = "| c "
        njetRow = ""
        for njet in self.njets:
            rowHeader += "| c " 

            operation = "="
            if "incl" in njet:
                operation = "\\geq"

            njetRow += ("& \\textbf{$\\mathbf{\\njets%s%s}$} "%(operation, njet.replace("incl", ""))).ljust(self.njetFlush)

        rowHeader += "|| c "
        njetRow += ("& \\textbf{$\\mathbf{\\njets\\geq%s}$} "%(self.njets[0])).ljust(self.njetFlush)
            
        indent = "    "

        header = ""
        header += "\\resizebox{\linewidth}{!}{%"
        header += "\n"
        header += "%s\def\\arraystretch{1.0}"%(indent)
        header += "\n"
        header += "%s\\begin{tabular}{%s|}"%(indent,rowHeader)
        header += "\n"
        header += "%s\hline"%(2*indent)
        header += "\n"

        temp = "\\textbf{Process}".ljust(self.sampleFlush)

        header += "%s%s %s \\\\"%(2*indent, temp, njetRow)
        header += "\n"
        header += "%s\hline"%(2*indent)
        header += "\n"

        return header

    def writeFooter(self):
        footer = ""
        footer += "    \end{tabular}"
        footer += "\n"
        footer += "}"
        footer += "\n"
        footer += "\n"

        return footer

    def makeCell(self, pair, summedPair, doFractions=False):
        quantity    = pair[0]
        uncertainty = pair[1]**0.5
        if doFractions:
            summedQuantity    = summedPair[0]
            summedUncertainty = summedPair[1]**0.5

            # Put uncertainty on percentage in range 0 to 100
            if quantity == 0.0:
                quantity += 10e-10

            uncertainty = 100 * (quantity/summedQuantity) * ((summedUncertainty / summedQuantity)**2.0 + (uncertainty / quantity)**2.0)**0.5

            quantity /= summedQuantity
            quantity *= 100.0

        aCell = ""
        if quantity < 1.0:
            aCell = (" & $%.3f \\pm %.3f$"%(quantity, uncertainty)).ljust(self.njetFlush)
        elif quantity < 10.0:
            aCell = (" & $%.2f \\pm %.2f$"%(quantity, uncertainty)).ljust(self.njetFlush)
        elif quantity > 100.0:
            aCell = (" & $%.1f \\pm %.1f$"%(quantity, uncertainty)).ljust(self.njetFlush)
        else:
            aCell = (" & $%s \\pm %s$"%(float("%.3g"%(quantity)), float("%.1f"%(uncertainty)))).ljust(self.njetFlush)
        
        isNegligible = doFractions and quantity < 0.1

        return aCell, isNegligible

    def makeTable(self, tableName, doFractions=False):

        # Create the output directory if it does not exist
        if not os.path.isdir(self.outputDir):
            os.makedirs(self.outputDir) 

        table = ""
        table += self.writeHeader()

        for keyVal in sorted(self.bkgSamples.items(), key=lambda t:t[1][0], reverse=True):

            sample = keyVal[0]
            if "AllBkg" in sample:
                continue

            newSample = self.niceNames[sample].replace("_", "\\_")
            if "mStop" in sample:
                chunks = sample.split("_")
                mass   = chunks[-1].split("-")[-1]
                signal = chunks[0]

                newSample = "%s $\\mathbf{m_{\\tilde{t}}=%s}$~GeV"%(signal, mass)

            aRow = ("\\textbf{%s}"%(newSample)).ljust(self.sampleFlush)
            aRow = 2*"    " + aRow

            isNegligible = True
            for njet in self.njets: 

                key = sample + "_Njet" + njet.replace("incl", "")

                if key in self.data:
                    aCell, njetIsNegligible = self.makeCell(self.data[key], self.data[key.replace(sample, "AllBkg")], doFractions)
                    aRow += aCell
                    isNegligible &= njetIsNegligible
                else:
                    aRow += " & N/A".ljust(self.njetFlush)
                
            aRow += self.makeCell(self.bkgSamples[sample],  self.bkgSamples["AllBkg"], doFractions)[0]

            if doFractions and isNegligible:
                aRow = aRow.replace(newSample, "\\textcolor{black}{%s}"%(newSample))

            table += aRow
            table += " \\\\\\hline \n"

        if not doFractions:

          for sample in natural_sort(self.sigSamples.keys()):

             newSample = sample
             if "mStop" in sample:
                 chunks = sample.split("_")
                 mass   = chunks[-1].split("-")[-1]
                 signal = chunks[0]

                 newSample = "%s $\\mathbf{m_{\\tilde{t}}=%s}$~GeV"%(signal, mass)

             aRow = ("\\textbf{%s}"%(newSample)).ljust(self.sampleFlush)
             aRow = 2*"    " + aRow

             for njet in self.njets: 
                 key = sample + "_Njet" + njet.replace("incl", "")

                 aRow += self.makeCell(self.data[key], None, doFractions)[0]

             aRow += self.makeCell(self.sigSamples[sample], None, doFractions)[0]

             table += aRow
             table += " \\\\\\hline \n"

        table += self.writeFooter()

        f = open("%s/%s_%s.tex"%(self.outputDir, tableName, self.channel), "w")
        f.write(table)
        f.close()

if __name__ == "__main__":

    usage = "usage: %tableYields [options]"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("--channel",   dest="channel",   help="which channel to process",     required=True)
    parser.add_argument("--inputDir",  dest="inputDir",  help="directory for input ROOT",     required=True)
    parser.add_argument("--outputDir", dest="outputDir", help="where to put tex table files", required=True)
    parser.add_argument("--QCDCR",     dest="QCDCR",     help="do for QCD CR",                default=False, action="store_true")
    parser.add_argument("--year",      dest="year",      help="which year to process",        default="Run2")
    args = parser.parse_args()

    extraStr = ""
    if args.QCDCR:
        extraStr = "_QCDCR"

    njets = ["7", "8", "9", "10", "11", "12incl"]
    inclusiveBin = 12
    if args.channel == "0l":
        njets = ["8", "9", "10", "11", "12", "13incl"]
        inclusiveBin = 13
    if args.channel == "2l":
        njets = ["6", "7", "8", "9", "10", "11incl"]
        inclusiveBin = 11

    theYields = Yields(args.outputDir, njets, args.channel)

    # Loop over all ROOT files in the input directory to read a histogram
    for rootFile in os.listdir(args.inputDir):

        chunks = rootFile.partition("_")

        year   = chunks[0]
        sample = chunks[-1].rstrip(".root")

        # Skip any Data ROOT files as well as those for TTX and BG_OTHER
        if "Data" in sample or "TTX" in sample or "BG_OTHER" in sample:
            continue

        # Only put a few signal
        if "mStop" in sample and (not any(mass in sample for mass in ["350", "550", "850", "1150"])):
            continue

        # If considering a single year, then pay attention to the year
        # in the name of the ROOT files
        if args.year != "Run2" and args.year != year:
            continue

        # Don't use Run2 files that are present, we will
        # combine individual year files ourselves
        if args.year == "Run2" and args.year == year:
            continue

        f = ROOT.TFile.Open(args.inputDir + "/" + rootFile, "READ")
        h = f.Get("h_Njets_%s%s_ABCD"%(args.channel, extraStr))

        # Extract per-bin information from the histogram
        # and store in an internal dictionary
        theYields.processHisto(sample, inclusiveBin, h)

    # Now that all backgrounds have been read in
    # Calculate and store the same information for total background
    theYields.processHisto("AllBkg")

    # Make per-Njets yields and fractions LaTeX tables
    theYields.makeTable("%s_Yields"%(args.year))
    theYields.makeTable("%s_Fractions"%(args.year), doFractions=True)
