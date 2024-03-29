import os, ROOT, argparse, shutil
from datetime import datetime

# The Splitter class is designed to split a signal file
# which could contain any stop mass point for a given even
# and produces N output ROOT files for the N mass points
# present in the original file
class Splitter:

    # inputFile : complete path to an input signal ROOT file
    # ttreePath : full path inside of ROOT file to get TreeMaker TTree
    def __init__(self, inputFile, ttreePath):

        self.inputFile = inputFile
        self.ttreePath = ttreePath
        self.stopMasses = list(range(300,1450,50))
        self.ctaus = [""]
        if "300to1500" in inputFile:
            self.stopMasses = [300,500,700,900,1100,1300,1500]
            self.ctaus = [0.01,0.1,1,10,100,1000]

        self.splitFile()

    # Get a nice little time stamp for printing
    def timeStamp(self):
        return datetime.now().strftime("%Y/%m/%d %H:%M:%S")

    # Open up the input ROOT file and tree and use the CopyTree
    # method with a cut on the variable SignalParameters, which
    # contains the mother (stop) mass for an event
    def splitFile(self):

        f = ROOT.TFile.Open(self.inputFile, "READ")
        t = f.Get(self.ttreePath)

        chunks   = list(filter(len, self.inputFile.split("/")))
        oldName  = chunks[-2]
        fileStub = chunks[-1]

        # "Dumb" loop over all possible masses each time
        # Not every mass will be found in a given input file
        # So just delete output files that contain no events
        for stopMass in self.stopMasses:

            if "300to1400" in self.inputFile:
                newName = oldName.replace("300to1400", str(stopMass))

                os.makedirs(newName)

                outFile = ROOT.TFile("%s/%s"%(newName,fileStub), "RECREATE")
                newTree = t.CopyTree("SignalParameters[0]==%d"%(stopMass))

                if newTree.GetEntries() != 0:
                    outFile.Write()
                    outFile.Close()
                    print("%s [INFO]: Created file: %s/%s"%(self.timeStamp(), newName, fileStub))
                else:
                    outFile.Close()
                    shutil.rmtree(newName)
                    print("%s [INFO]: Removed empty file: %s"%(self.timeStamp(), newName))
    
            else:
                for lspMass in [100, stopMass - 225]:
                    for ctau in self.ctaus:

                        newName = oldName.replace("300to1500", str(stopMass)).replace("lowandhigh", str(lspMass)).replace("0p01to1000", str(ctau).replace(".", "p"))

                        os.makedirs(newName)
        
                        outFile = ROOT.TFile("%s/%s"%(newName,fileStub), "RECREATE")
                        newTree = t.CopyTree("SignalParameters[0]==%d&&SignalParameters[1]==%d&&SignalParameters[2]>=%f&&SignalParameters[2]<=%f"%(stopMass,lspMass,0.5*ctau,1.5*ctau))

                        if newTree.GetEntries() != 0:
                            outFile.Write()
                            outFile.Close()
                            print("%s [INFO]: Created file: %s/%s"%(self.timeStamp(), newName, fileStub))
                        else:
                            outFile.Close()
                            shutil.rmtree(newName)
                            print("%s [INFO]: Removed empty file: %s"%(self.timeStamp(), newName))

# User can choose the path to the input ROOT file and the path to the TTree name from the command line
# Example call:
# python runSplitSignal.py --inputFile root://cmseos.fnal.gov///store/user/semrat/SusyRA2Analysis2015/Run2ProductionV20/Summer20UL16APV.RPV_2t6j_mStop-300to1400_mN1-100_TuneCP5_13TeV-madgraphMLM-pythia8_96_RA2AnalysisTree.root
#                          --ttreePath "TreeMaker2/PreSelection"
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--inputFile", dest="inputFile", help="Input file to split", type=str, required=True)
    parser.add_argument("--ttreePath", dest="ttreePath", help="TTree name to read" , type=str, default="TreeMaker2/PreSelection")
    args = parser.parse_args()

    theSplitter = Splitter(args.inputFile, args.ttreePath)
