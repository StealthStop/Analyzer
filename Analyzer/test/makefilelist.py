import os, re, argparse
from collections import defaultdict, OrderedDict

# Usage: Populates a directory "filelists_Kevin_<versioning>" with a text file for each year + sample combination.
#        Contents of the text file are simply a listing of paths to all relevant ROOT files for the year + sample
#        The "filelists_Kevin_<versioning>" should be placed in /eos/uscms/store/user/lpcsusyhad/StealthStop/
#        where it will be referenced by sampleSets.cfg. Likewise, a sampleSets_<versioning>.cfg is constructed

class FileLister:

    def __init__(self, production, tag):

        self.production   = production
        self.tag          = tag

        self.filesDir     = "/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2Production%s/"%(production)
        self.ttreeDir     = "TreeMaker2/PreSelection"
        self.tempFileName = "tmp.txt"
        self.fileListsDir = "filelists_Kevin_%s/"%(production)
        self.eosPath      = "/eos/uscms/store/user/lpcsusyhad/StealthStop/"

        self.treeMakerDir = "/uscms/home/jhiltb/nobackup/susy/ZeroAndTwoLep/CMSSW_10_6_29_patch1/src/TreeMaker/WeightProducer/python"

        # Dictionary for holding information about samples
        # Info pulled from TreeMaker
        self.auxInfo = {}

        # Will hold all ROOT files for each sample for each year
        self.samples = OrderedDict()
        self.samples["2016"]    = defaultdict(list)
        self.samples["2016APV"] = defaultdict(list)
        self.samples["2017"]    = defaultdict(list)
        self.samples["2018"]    = defaultdict(list)

        if not os.path.isdir(self.fileListsDir):
            os.mkdir(self.fileListsDir)

        # ls the directory containing all TreeMaker ntuple files
        # and pipe to temp output txt file
        command = os.system("eos root://cmseos.fnal.gov ls %s > tmp.txt"%(self.filesDir))
        self.collectFiles("tmp.txt")

        # If the TreeMaker area exists, then use it to read in event counts and xsections
        # These functions will populate self.auxInfo
        self.getEventCounts()
        self.getXsecInfo()

        # Write out the txt file with corresponding list of files for a sample
        # Also writes the sample set file
        self.writeFileLists()

    # Take a sample name e.g. "2016_QCD_Pt_120to170_TuneCP5_13TeV_pythia8"
    # and return everything before the "TuneCP5" => "2016_QCD_Pt_120to170."
    # Special case if "TuneCP5" is modified e.g "TuneCP5down", then include
    # the TuneCP5down in the final name
    def makeNiceName(self, oldname):
        chunks = oldname.split("_")
        endpoint = 0
        for chunk in chunks:
            if "TuneCP5" in chunk:
                if chunk != "TuneCP5":
                    endpoint += 1
                break
        
            endpoint += 1
    
        return "_".join(chunks[0:endpoint]).replace("_NLO", "")
    
    # Search dictionary self.auxInfo for all entries with keys containing
    # "name" and add "payload" subdictionary with matching key. 
    # Do not mix inclusive sample names with those for particular HT or Pt bins
    def addInfo(self, name, payload):
        
        for sample in self.auxInfo.keys():
            if name in sample:
                if ("_HT" in sample or "_Pt" in sample) and name != sample:
                    continue 
    
                self.auxInfo[sample].update(payload) 
    
    # Given a sample string e.g. "TTTW", determine which string
    # in the samplesGroup list contains the substring
    def findGroup(self, sampleGroups, aSample):

        chunks = aSample.split("_")
        uniqueStr = chunks[1]

        if len(chunks) > 2:
            if "HT" in chunks[2]:
                uniqueStr += "_HT"
            elif "Pt" in chunks[2]:
                uniqueStr += "_Pt"
        
        bestGroupMatch = "NULL"
        for sampleGroup in sampleGroups:
        
            samples = sampleGroup.split("|")

            for sample in samples:
    
                if uniqueStr == sample:
                    return sampleGroup
                    
                elif sample in uniqueStr:
                    bestGroupMatch = sampleGroup
    
        return bestGroupMatch
    
    def naturalSort(self, unsortedList): 

        convert      = lambda text: int(text) if text.isdigit() else text.lower() 
        alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)] 

        return sorted(unsortedList, key=alphanum_key)
    
    # Extract event numbers from MCSample files in the TreeMaker area
    def getEventCounts(self):
    
        for year in ["16", "16APV", "17", "18"]:
            if os.path.exists(self.treeMakerDir + "/MCSamples_Summer20UL%s.py"%(year)):
                countsFile = open(self.treeMakerDir + "/MCSamples_Summer20UL%s.py"%(year))
                
                lines = countsFile.readlines()
                countsFile.close()
        
                for line in lines:
                    if "MCSample(" not in line:
                        continue
        
                    temp = line.replace("'", "").replace(",", "")
                    temp = temp[temp.find("(")+1:temp.find(")")]
                    chunks = temp.split(" ")
        
                    process = chunks[0]
        
                    ntot = -1.0; ndiff = -1.0
        
                    if len(chunks) == 5:
                        ntot = chunks[-1]
                    elif len(chunks) > 5:
                        ntot = chunks[-2]
                        ndiff = chunks[-1]
        
                    name = self.makeNiceName(process)
                    if name not in self.auxInfo.keys():
                        self.auxInfo[name] = {}
        
                    if ndiff == -1.0:
                        self.auxInfo[name]["20%s_npos"%(year)] = int(ntot)
                        self.auxInfo[name]["20%s_nneg"%(year)] = 0
                    else:
                        self.auxInfo[name]["20%s_npos"%(year)] = int((float(ntot) + float(ndiff))/2.0)
                        self.auxInfo[name]["20%s_nneg"%(year)] = int((float(ntot) - float(ndiff))/2.0)
    
    # If there is a MCSamplesValues, use it to get info on xsec, kfactor, etc
    def getXsecInfo(self):
    
        if os.path.exists(self.treeMakerDir + "/MCSampleValues.py"):
            xsecFile = open(self.treeMakerDir + "/MCSampleValues.py")
            lines = xsecFile.readlines()
            xsecFile.close()
            
            seenDict = False
            process = ""
            xsec    = "-1.0" 
            brf     = "1.0"
            kfactor = "1.0"
            for line in lines:
            
                if "values_dict =" in line:
                    seenDict = True
            
                if not seenDict:
                    continue
            
                if ":" in line and "{" in line:
                    
                    if process != "":
                        self.addInfo(process, {"xsec" : xsec, "brf" : brf, "kfactor" : kfactor})
                        process = ""; xsec = "-1.0"; brf = "1.0"; kfactor = "1.0"
            
                    process = line.split("\"")[1] 
            
                if "XS_" in line:
                    xsec = line.split("=")[1].split(",")[0]
            
                if "BR_" in line:
                    brf = line.split("=")[1].split(",")[0] 
            
                if "kFactor_" in line:
                    kfactor = line.split("=")[1].split(",")[0]
   
    # Determine the common sample names amongst the files and 
    # group all these files under the common name 
    def collectFiles(self, tempFilePath):
        with open(tempFilePath, 'r') as tempfile:
            for line in tempfile:
        
                # Each line will be simple ROOT file name e.g. Summer20UL18.ttHJetTobb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8_9_RA2AnalysisTree.root
                if ".root" not in line or "Fast" in line:
                    continue
        
                era = ""
                if   ("UL2016" in line and "HIPM" in line) or "UL16APV" in line: 
                    era = "2016APV"
                elif "UL2016" in line or "UL16" in line: 
                    era = "2016"
                elif "UL2017" in line or "UL17" in line: 
                    era = "2017"
                elif "UL2018" in line or "UL18" in line: 
                    era = "2018"
                else:
                    continue
        
                # With example above, drop the "Summer20UL18" and ".root", 
                # grab string before the last two "_" i.e. "ttHJetTobb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8"
                shortline = line.split(".")[1].rpartition("_")[0].rpartition("_")[0]
        
                # Now we construct shortline as "2018_ttHJetTobb" and remove extra things
                shortline = era + "_" + shortline.replace("_ext1", "").replace("_ext2", "").replace("_ext3", "").replace("_backup", "")
                
                # Use "2018_ttHJetTobb" as key to list of all ttHJetTobb files for 2018
                newline = "root://cmseos.fnal.gov/" + self.filesDir + line
                self.samples[era][shortline].append(newline)
    
    # Write out a txt file for each sample with the list of corresponding files
    def writeFileLists(self):
    
        sampleGroups = ["Data", "TTJets", "DYJetsToLL", "TTToSemiLeptonic", "TTTo2L2Nu", "TTToHadronic", 
                        "WJetsToLNu", "QCD_HT", "QCD_Pt", "TTTT|TTTJ|TTTW|TTTZ|TTTH|TTWW|TTWZ|TTZZ|TTHH|TTWH|TTZH",
                        "WWW|WWG|WWZ|WZZ|ZZZ|WZG", "WW|WZ|ZZ", "TTZTo", "TTWJets", "WWTo|ZZTo|WZTo", "ttHJet" 
        ]

        for era, sampleLists in self.samples.items():
        
            finalDict = {sampleGroup : {} for sampleGroup in sampleGroups}
        
            # For some sanity, sort the samples alphabetically
            theSamples = sampleLists.keys()
            theSamples.sort()
            for sample in theSamples:
                newfile = open("filelists_Kevin_%s/"%(self.tag) + sample + ".txt", 'w')
        
                xsec       = "-1.0,"
                nposevents = "0,"
                nnegevents = "0,"
                kfactor    = "1.0"
        
                # Everything before "TuneCP5" will be included in name
                # "ttHJetTobb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8" ==> "ttHJetTobb_M125"
                name    = self.makeNiceName(sample)
                auxName = name.partition("_")[-1]
        
                # Do not care about these data sets
                if "SinglePhoton" in name or "MET" in name or "MHT" in name:
                    continue
        
                # For the inclusive W + jets and DY + jets samples add "Incl" string to name
                if (name.startswith("WJetsToLNu") or name.startswith("DYJetsToLL")) and "HT" not in name:
                    name += "_Incl"
        
                # For any data set, insert "Data" string into name
                name = name.replace("SingleMuon", "Data_SingleMuon").replace("SingleElectron", "Data_SingleElectron") \
                           .replace("JetHT", "Data_JetHT").replace("EGamma", "Data_SingleElectron")
        
                # Special case to get "erdON" into name
                if "erdON" in sample:
                    name += "_erdON"

                # Restrict to madgraph DY sample
                if "DYJetsToLL" in name in "NLO" in name:
                    continue
        
                # Special case in some boson samples to remove extra string
                name = name.replace("_4F", "").replace("_4f", "").replace("_M125", "")
        
                for f in sampleLists[sample]:
                    newfile.write(f)
                newfile.close()
        
                # Get any xsec, kfactor info that was retrieved from TreeMaker
                if auxName in self.auxInfo.keys():
                    xsec       = str(float(self.auxInfo[auxName]["xsec"]) * eval(self.auxInfo[auxName]["brf"])) + ","
                    kfactor    = str(self.auxInfo[auxName]["kfactor"])
        
                    if "%s_npos"%(era) in self.auxInfo[auxName]:
                        nposevents = str(self.auxInfo[auxName]["%s_npos"%(era)]) + ","
                        nnegevents = str(self.auxInfo[auxName]["%s_nneg"%(era)]) + ","
                       
                sampleGroup = self.findGroup(sampleGroups, name)

                if sampleGroup == "NULL":
                    print("No group found for sample with name \"%s\"---skipping !"%(name))
                    continue

                name   += ","
                sample += ".txt,"
        
                finalDict[sampleGroup][name] = "%s %s, %s %s, %s %s %s %s\n"%(name.ljust(40), self.eosPath + "/" +self.fileListsDir, sample.ljust(85), self.ttreeDir, xsec.rjust(14), nposevents.rjust(12), nnegevents.rjust(12), kfactor.rjust(6))
    
            self.writeSampleSet(finalDict)

    # Append lines to the sampleSets file for each sample
    # Similar samples will be grouped together
    def writeSampleSet(self, dictionary):

        sampleSet = open("sampleSets_%s.cfg"%(self.tag), "a")

        sortedGroups = dictionary.keys()
        sortedGroups.sort()
        for sampleGroup in sortedGroups:
        
            samples = self.naturalSort(dictionary[sampleGroup].keys())
        
            writeNewline = False
            for sample in samples:
                writeNewline = True
                sampleSet.write(dictionary[sampleGroup][sample])
        
            if writeNewline:
                sampleSet.write("\n")
        
        sampleSet.close()

# Run the script
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--prod", dest="prod", help="Unique tag for output", type=str, default="V20")
    parser.add_argument("--tag" , dest="tag" , help="Path to PU file"      , type=str, default="UL_v1")
    args = parser.parse_args()
    
    theLister = FileLister(args.prod, args.tag)
