import os, ROOT, argparse
from datetime import datetime

class Splitter:

    def __init__(self, inputFile, ttreeName):

        self.inputFile = inputFile
        self.ttreeName = ttreeName
        self.masses = list(range(300, 1450, 50))

        self.splitFile()

    def timeStamp(self):
        return datetime.now().strftime("%Y/%m/%d %H:%M:%S")

    def splitFile(self):

        f = ROOT.TFile.Open(self.inputFile, "READ")
        t = f.Get(self.ttreeName)

        oldName = self.inputFile.rpartition("/")[-1]

        for mass in self.masses:

            newName = oldName.replace("300to1400", str(mass))
        
            outFile = ROOT.TFile(newName, "RECREATE")
            newTree = t.CopyTree("SignalParameters[0]==%d"%(mass))

            if newTree.GetEntries() != 0:
                outFile.Write()
                outFile.Close()
                print("%s [INFO]: Created file: %s"%(self.timeStamp(), newName))
            else:
                outFile.Close()
                os.remove(newName)
                print("%s [INFO]: Removed empty file: %s"%(self.timeStamp(), newName))

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--inputFile", dest="inputFile", help="Input file to split", type=str, required=True)
    parser.add_argument("--ttreeName", dest="ttreeName", help="TTree name to read" , type=str, default="TreeMaker2/PreSelection")
    args = parser.parse_args()

    theSplitter = Splitter(args.inputFile, args.ttreeName)
