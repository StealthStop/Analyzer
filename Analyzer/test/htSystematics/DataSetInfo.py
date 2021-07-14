import ROOT

class DataSetInfo:
    def __init__(self, basedir, fileName, label, sys=None, processName=None, process=None, rate=None, lumiSys=None, scale=-1.0):
        self.basedir = basedir
        self.fileName = fileName
        self.label = label
        try: 
            self.file = ROOT.TFile.Open(basedir+fileName)
        except:
            pass
        self.sys = sys
        self.processName = processName
        self.process = process
        self.rate = rate
        self.lumiSys = lumiSys
        self.scale = scale

    def getHisto(self, name, scale=-1):
        histo = self.file.Get(name)
        if(self.scale != -1.0):
            #print self.scale
            histo.Scale(self.scale)
        return histo

    def getFile(self):
        return self.file

    def __del__(self):
        self.file.Close()
