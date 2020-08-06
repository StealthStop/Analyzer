import ROOT
import numpy as np
import math
from optparse import OptionParser
import root_numpy as rnp
import glob
from ROOT import TFile, gROOT, gStyle

parser = OptionParser()
parser.add_option("-o", "--output", action="store", type="string", dest="filename",
    default="datacard.txt", help="Name of desired output file")
parser.add_option("-s", "--signal", action="store", type="string", dest="signalName",
    default="RPV350", help="Name of signal root file to use as input")
parser.add_option("-b", "--binSize", action="store", type="int", dest="binSize",
     default="100", help="Desired size of bins in GeV (square bins)")

(options, args) = parser.parse_args()

def my_range(start, end, step):
    while start <= end:
        yield start
        start += step

def checkBins(binVals):
    x = 0
    for b in binVals[0]:
        if (binVals[0][x] < 0.1 and binVals[1][x] < 0.1 and binVals[2][x] < 0.1):
            del binVals[0][x]
            del binVals[1][x]
            del binVals[2][x]
        x += 1
    return binVals
    

def writeDataCard(binVals):
    processNames = [options.signalName, "QCD", "TT"]

    with open("Card%s.txt" % options.signalName, 'w') as f:
        f.write("Datacard for 2016 %s" % options.signalName)
        f.write( "\n" )
        f.write( "imax %d number of bins\n" % len(binVals[0]) )
        f.write( "jmax 2 number of processes minus 1\n" )
        f.write( "kmax 1 number of nuisance parameters\n" )
        f.write( "\n" )
        f.write( "-------------------------------------------------------------------------------------------------------------------------------------------\n" )
        f.write( "\n" )
        bins        = "bin         "
        observation = "observation "
        for i in range(len(binVals[0])):
            bins        += "D{0: <4}".format(i)
            observation += "{0: <5}".format(round(binVals[0][i],1))
        f.write( bins+"\n" )
        f.write( observation+"\n" )
        f.write( "-------------------------------------------------------------------------------------------------------------------------------------------\n" )
        f.write( "# background rate taken from events in bin\n" )
        bins      = "bin              "
        process1 = "process          "
        process2 = "process          "
        rate     = "rate             "
        f.write( "" )
        for i in range(len(binVals[0])):
            for j in range(3):
                bins      += "{0: <16} ".format("D"+str(i))
                process1 += "{0: <16} ".format(processNames[j])
                process2 += "{0: <16} ".format(j)
                rate     += "{0: <16} ".format(binVals[j][i])
        f.write( bins+"\n" )
        f.write( process1+"\n" )
        f.write( process2+"\n" )
        f.write( rate+"\n" )
        f.write( "-------------------------------------------------------------------------------------------------------------------------------------------\n" )
        f.write( "# Normal uncertainties in the signal region\n" )
        lumiSys = "lumi_13TeV      lnN  "
        for i in range(len(binVals[0])):
            for j in range(3):
                if(j == 0):
                    lumiSys += "{0: <5} ".format("1.05")
                else:
                    lumiSys += "{0: <5} ".format("-")
        f.write( lumiSys+"\n" )
        f.write( "-------------------------------------------------------------------------------------------------------------------------------------------\n" )

def main():
    #gROOT.SetBatch(True)
    gStyle.SetOptStat(0)

    # Create list of histograms that we will need to pull data from
    flist = []
    hlist = []

    for f in glob.glob("/uscms/home/bcrossma/nobackup/CMSSW_10_2_9/src/Analyzer/Analyzer/test/condor/StopMass_PtRank/output-files/rootFiles/*.root"):
        if(f.find(options.signalName) != -1):
            flist = [ROOT.TFile.Open(f)] + flist
        elif (f.find("TT") != -1 or f.find("QCD") != -1):
            flist.append(ROOT.TFile.Open(f))

#    print flist

    for f in flist:
        hlist.append(f.Get("h_stopMasses_diffVSavg_0l_HT500_ge2b_ge6j_ge2t_ge1dRbjets"))


    # Define array of bin boundaries
    binVals = [[],[],[]]
    binInt = 0
    i = 0

    for x in my_range(1, 150, 30):
        for y in my_range(1, 150, 30):
            i = 0
            temp_binVals = []
            for h in hlist:
                binInt = h.Integral(x, x+10, y, y+10)
                if (round(binInt,1) < 0.1):
                    break
                else:
                    temp_binVals.append(round(binInt,1))
                i += 1
            j = 0
            for b in temp_binVals:
                if (len(temp_binVals) == len(flist)):
                    binVals[j].append(b)
                j += 1
    
#    binVals = checkBins(binVals)

    print binVals
#    print "TTbar: %d, QCD: %d, Signal: %d" % (len(binVals[0]), len(binVals[1]), len(binVals[2]))

    print len(binVals[0])
    writeDataCard(binVals)

if __name__ == "__main__":
    main()

