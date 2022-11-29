import os, re, argparse, subprocess
from collections import defaultdict, OrderedDict

# Usage: Populates a directory "filelists_Kevin_<versioning>" with a text file for each year + sample combination.
#        Contents of the text file are simply a listing of paths to all relevant ROOT files for the year + sample
#        The "filelists_Kevin_<versioning>" should be placed in /eos/uscms/store/user/lpcsusyhad/StealthStop/
#        where it will be referenced by sampleSets.cfg. Likewise, a sampleSets_<versioning>.cfg is constructed

class FileLister:

    def __init__(self, production, tag, doSkim=False):

        self.production   = production
        self.tag          = tag
        self.doSkim       = doSkim

        # Whether we are generating file lists for a TreeMaker ntuple data set or a skimmed data set
        # created by our MakeMiniTree, we need to adjust some input paths, tree names, output folder names, etc
        if not doSkim:
            self.filesDir     = "/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2Production%s"%(production)
            self.ttreePath    = "TreeMaker2/PreSelection"
            self.ttreePathSig = "PreSelection"
            self.tempFileName = "tmp.txt"
            self.workingDir   = "filelists_Kevin_%s/"%(production)
        else:
            self.filesDir     = "/store/user/lpcsusystealth/StealthStop/Skims"
            self.ttreePath    = "PreSelection"
            self.ttreePathSig = "PreSelection"
            self.tempFileName = "tmp.txt"
            self.workingDir   = "filelists_Kevin_%s_skim/"%(production)

        self.fileListPath = "filelists"

        # This works as long as the user jhiltb maintains a TreeMaker setup here and it is up-to-date
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

        self.sigMasses = list(range(300, 1450, 50))

        self.wroteHeader = False

        if not os.path.isdir(self.workingDir):
            os.mkdir(self.workingDir)

        # Get listing of all available ROOT files
        self.collectFiles()

        # If the TreeMaker area exists, then use it to read in event counts and xsections
        # These functions will populate self.auxInfo
        self.getEventCounts()
        self.getXsecInfo()

        # Write out the txt file with corresponding list of files for a sample
        # Also writes the sample set file
        self.writeFileLists()

    # Wrapper to list things in an EOS path
    # Specifically, find all files that were modified within the last 30 days
    def listFilesEOS(self, path):
        proc = subprocess.Popen(["eos", "root://cmseos.fnal.gov", "find", "-mtime", "-30", path], stdout=subprocess.PIPE)
        payload = proc.stdout.readlines()
        lines = []
        for line in payload:
            if ".root" not in line:
                continue
            temp = str(line.strip())
            fileName = temp.split(" ")[0].split("/")[-1]
            lines.append(fileName)    

        return self.naturalSort(lines)

    # Specific wrapper to summon eos command to do a simple ls
    # Used for listing directories without needing to know modification time
    def listDirsEOS(self, path):
        proc = subprocess.Popen(["eos", "root://cmseos.fnal.gov", "ls", path], stdout=subprocess.PIPE)
        payload = proc.stdout.readlines()
        lines = [str(line.strip()) for line in payload]

        return lines

    # Take a sample name e.g. "2016_QCD_Pt_120to170_TuneCP5_13TeV_pythia8"
    # and return everything before the "TuneCP5" => "2016_QCD_Pt_120to170."
    # Special case if "TuneCP5" is modified e.g "TuneCP5down", then include
    # the TuneCP5down in the final name
    def makeNiceName(self, oldname):
        chunks = oldname.split("_")
        endpoint = 0
        extra = ""
        for chunk in chunks:
            if "TuneCP5" in chunk:
                if chunk != "TuneCP5":
                    endpoint += 1
                break

            if "erdON" in oldname:
                extra = "_erdON" 
        
            endpoint += 1
    
        return "_".join(chunks[0:endpoint]).replace("_NLO", "") + extra
    
    # Search dictionary self.auxInfo for all entries with keys containing
    # "name" and add "payload" subdictionary with matching key. 
    # Do not mix inclusive sample names with those for particular HT or Pt bins
    def addInfo(self, name, payload):
        
        for sample in self.auxInfo.keys():
            if name in sample:
                if ("_HT" in sample or "_Pt" in sample) and name != sample:
                    continue 

                if ("_4F" in sample or "_4f" in sample) and name != sample:
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
    
    # Function taken from a stack overflow post which sorts a list naturally i.e.
    # ["1_file.root", "11_file.root", "12_file.root", "2_file.root"] ===>
    # ["1_file.root", "2_file.root", "11_file.root", "12_file.root"]
    def naturalSort(self, unsortedList): 

        convert      = lambda text: int(text) if text.isdigit() else text.lower() 
        alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)] 

        return sorted(unsortedList, key=alphanum_key)
    
    # Extract event count numbers from MCSample files in the TreeMaker area
    # The MCSample files report NposEvts+NnegEvts and NposEvts-NnegEvts
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
        
                    name  = self.makeNiceName(process)
                    names = []

                    if ("RPV" in name or "SYY" in name or "SHH" in name) and "300to1400" in name:
                        names = [name.replace("300to1400", str(mass)) for mass in self.sigMasses]
                    else:
                        names = [name]
                    
                    for aName in names:
                        if aName not in self.auxInfo.keys():
                            self.auxInfo[aName] = {}
        
                        if ndiff == -1.0:
                            self.auxInfo[aName]["20%s_npos"%(year)] = int(ntot)
                            self.auxInfo[aName]["20%s_nneg"%(year)] = 0
                        else:
                            self.auxInfo[aName]["20%s_npos"%(year)] = int((float(ntot) + float(ndiff))/2.0)
                            self.auxInfo[aName]["20%s_nneg"%(year)] = int((float(ntot) - float(ndiff))/2.0)
    
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
            
                if (":" in line and "{" in line) or ("}" in line and "," not in line and seenDict):
                    
                    if process != "":
                        self.addInfo(process, {"xsec" : xsec, "brf" : brf, "kfactor" : kfactor})
                        process = ""; xsec = "-1.0"; brf = "1.0"; kfactor = "1.0"
            
                    if "}" not in line:
                        process = line.split("\"")[1] 
            
                if "XS_" in line:
                    xsec = line.split("=")[1].split(",")[0]
            
                if "BR_" in line:
                    brf = line.split("=")[1].split(",")[0] 
            
                if "kFactor_" in line:
                    kfactor = line.split("=")[1].split(",")[0]
   
    # Determine the common sample names amongst the files and 
    # group all these files under the common name 
    def collectFiles(self):

        allFileList = []

        extraSubDir = ""
        if self.doSkim:
            extraSubDir = "/output-files/"

        # Will get a list of folders corresponding to data and MC eras
        # i.e. Summer20UL16 or Run2017D-UL2017-v2
        eraDirs = self.listDirsEOS("%s/*"%(self.filesDir))
        for eraDir in eraDirs:

            # Skip any ROOT file at this level
            if ".root" in eraDir: continue

            # Within each eraDir, get the list of folders, each corresponding
            # to a single unique sample i.e. ZZZ_TuneCP5_13TeV-amcatnlo-pythia8
            sampleDirs = self.listDirsEOS("%s/%s%s"%(self.filesDir, eraDir, extraSubDir))
            for sampleDir in sampleDirs:

                if "300to1400" in sampleDir:
                    continue

                fileList = self.listFilesEOS("%s/%s%s/%s"%(self.filesDir, eraDir, extraSubDir, sampleDir))

                # Finally within each sampleDir is a set of ROOT files
                for aFile in fileList:
        
                    if ".root" not in aFile or "Fast" in aFile:
                        continue
        
                    era = ""
                    if   ("UL2016" in eraDir and "HIPM" in eraDir) or "UL16APV" in eraDir or "preVFP" in eraDir: 
                        era = "2016APV"
                    elif "UL2016" in eraDir or "UL16" in eraDir or "postVFP" in eraDir: 
                        era = "2016"
                    elif "UL2017" in eraDir or "UL17" in eraDir or "2017" in eraDir: 
                        era = "2017"
                    elif "UL2018" in eraDir or "UL18" in eraDir or "2018" in eraDir: 
                        era = "2018"
                    else:
                        continue
        
                    newName = sampleDir.replace("_ext1", "").replace("_ext2", "").replace("_ext3", "").replace("_backup", "")
                    newline = "root://cmseos.fnal.gov/" + self.filesDir + "/" + eraDir + "%s/"%(extraSubDir) + sampleDir + "/" + aFile + "\n"

                    # For skims, we will extra the newName per file
                    # This then will aline the naming with how it is setup in sampleSet.cfg
                    # Thus, if skims are for 2016_TT, it will still separate into TTToSemiLeptonic, TTToHadronic, TTTo2L2Nu
                    temp = None
                    if "MyAnalysis" in aFile:

                        # Start with aFile = /some/path/MyAnalysis_2016postVFP_Data_JetHT_420.root
                        # Now temp1 = 2016postVFP_Data_JetHT_420.root
                        temp = aFile.rpartition("MyAnalysis_")[-1]

                        # Now split off year to get temp1 = JetHT_420.root
                        temp = temp.partition("_")[-1]

                        # Finally split off anything related to job number
                        temp = temp.rpartition("_")[0]

                    else:
                        temp = newName
 
                    self.samples[era][era + "_" + temp].append(newline)

    
    # Write out a txt file for each sample with the list of corresponding files
    def writeFileLists(self):
    
        sampleGroups = ["Data", "TTJets", "DYJetsToLL", "TTToSemiLeptonic", "TTTo2L2Nu", "TTToHadronic", 
                        "WJetsToLNu|WJetsToQQ", "QCD_HT", "TTTT|TTTJ|TTTW|TTTZ|TTTH|TTWW|TTWZ|TTZZ|TTHH|TTWH|TTZH", "ST|tZq",
                        "WWW|WWG|WWZ|WZZ|ZZZ|WZG", "WW|WZ|ZZ", "TTZTo", "TTWJets", "WWTo|ZZTo|WZTo", "ttHJet", "RPV|StealthSYY|StealthSHH"
        ]

        for era, sampleLists in self.samples.items():
        
            finalDict = {sampleGroup : {} for sampleGroup in sampleGroups}
        
            # For some sanity, sort the samples alphabetically
            theSamples = sampleLists.keys()
            theSamples.sort()
            for sample in theSamples:

                # Some samples we do not care about at this time
                if "GJets" in sample or "TGamma" in sample or "QCD_Pt" in sample or "ZJetsToNuNu" in sample  or \
                   "MET" in sample or "HTMHT" in sample or "SinglePhoton" in sample or "genMET" in sample:
                    continue
            
                extraText = ""
                if self.doSkim:
                    extraText = "_skim"

                newfile = open(self.workingDir + sample + extraText + ".txt", 'w')
        
                xsec       = "-1,"
                nposevents = "-1,"
                nnegevents = "-1,"
                kfactor    = "1.0"
        
                # Everything before "TuneCP5" will be included in name
                # "ttHJetTobb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8" ==> "ttHJetTobb_M125"
                name    = self.makeNiceName(sample)
                auxName = name.partition("_")[-1]

                # For the inclusive W + jets and DY + jets samples add "Incl" string to name
                if ("_WJetsToLNu" in name or "_DYJetsToLL" in name or "_TTJets" in name) and "HT" not in name:
                    name += "_Incl"
        
                # For any data set, insert "Data" string into name
                if "Data" not in name:
                    name = name.replace("SingleMuon", "Data_SingleMuon").replace("SingleElectron", "Data_SingleElectron") \
                               .replace("JetHT", "Data_JetHT").replace("EGamma", "Data_SingleElectron")

                # Use the other alias for the two 2016s...the naming conventions are terrible...
                name = name.replace("2016_", "2016postVFP_").replace("2016APV_", "2016preVFP_")
        
                # Restrict to madgraph DY sample
                if "DYJetsToLL" in name and "NLO" in name:
                    continue
        
                # Special case in some boson and single top samples to remove extra string
                name = name.replace("_4F", "").replace("_4f", "").replace("_M125", "").replace("_5f_InclusiveDecays", "_Incl").replace("_5f_inclusiveDecays", "_Incl")
        
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

                # For skim sample sets, append the "_skim" suffix
                if self.doSkim:
                    name   += "_skim"
                    sample += "_skim"
                
                name   += ","
                sample += ".txt,"

                isSignal = "mStop" in name
                ttreePath = self.ttreePath
                if isSignal:
                    ttreePath = self.ttreePathSig

                    name = name.replace("_mSo-100", "") \
                               .replace("_mN1-100", "")

                if self.doSkim:
                    # Here we cannot get the number of events for the skims, because they were made outside TreeMaker
                    # So just put in the true number of events (based on dictionary in TreeMaker)
                    # We will run nEvts.py in any case to fix the numbers in the end
                    finalDict[sampleGroup][name] = "%s %s, %s %s, %s %s %s %s %s %s\n"%(name.ljust(45), self.fileListPath, sample.ljust(90), ttreePath.rjust(len(self.ttreePath)), xsec.rjust(14), nposevents.rjust(12), nnegevents.rjust(12), nposevents.rjust(12), nnegevents.rjust(12), kfactor.rjust(6))
                else:
                    finalDict[sampleGroup][name] = "%s %s, %s %s, %s %s %s %s\n"%(name.ljust(45), self.fileListPath, sample.ljust(90), ttreePath.rjust(len(self.ttreePath)), xsec.rjust(14), nposevents.rjust(12), nnegevents.rjust(12), kfactor.rjust(6))
   
            self.writeSampleSet(finalDict)

    # Append lines to the sampleSets file for each sample
    # Similar samples will be grouped together
    def writeSampleSet(self, dictionary):

        sampleSet = open("sampleSets_%s.cfg"%(self.tag), "a")

        if not self.wroteHeader:
            sampleSet.write("%s %s %s %s %s %s %s %s\n\n"%("# Sample name,".ljust(45), "filelists,".rjust(len(self.fileListPath)+1), "Sample_file_list.txt,".ljust(90), "TTree Name,".rjust(len(self.ttreePath)+1), "Xsec,".rjust(14), "+evt cnts,".rjust(12), "-evt cnts,".rjust(12), "kfact".rjust(6)))
            self.wroteHeader = True

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
    parser.add_argument("--prod", dest="prod", help="Unique tag for output", type=str,            default="V20")
    parser.add_argument("--tag" , dest="tag" , help="Path to PU file"      , type=str,            default="UL_v2")
    parser.add_argument("--skim", dest="skim", help="Make for skim"        , action="store_true", default=False)
    args = parser.parse_args()
    
    theLister = FileLister(args.prod, args.tag, args.skim)
