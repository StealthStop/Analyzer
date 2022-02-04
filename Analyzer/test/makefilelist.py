import os
from collections import defaultdict, OrderedDict

# Take a sample name e.g. "2016_QCD_Pt_120to170_TuneCP5_13TeV_pythia8"
# and return everything before the "TuneCP5" => "2016_QCD_Pt_120to170."
# Special case if "TuneCP5" is modified e.g "TuneCP5down"
def makeNiceName(oldname):
    chunks = oldname.split("_")
    endpoint = 0
    for chunk in chunks:
        if "TuneCP5" in chunk:
            if chunk != "TuneCP5":
                endpoint += 1
            break
    
        endpoint += 1

    return "_".join(chunks[0:endpoint]).replace("_NLO", "")

# Search dictionary "d" for all entries with keys containing
# "name" and add "payload" subdictionary". Do not mix inclusive sample
# names with those for particular HT or Pt bins
def addInfo(d, name, payload):
    
    for k in d.keys():
        if name in k:
            if ("_HT" in k or "_Pt" in k) and name != k:
                continue 

            d[k].update(payload) 

    return d

# Usage: Populates a directory "filelists_Kevin_<versioning>" with a text file for each year + sample combination.
#        Contents of the text file are simply a listing of paths to all relevant ROOT files for the year + sample
#        The "filelists_Kevin_<versioning>" should be placed in /eos/uscms/store/user/lpcsusyhad/StealthStop/
#        where it will be referenced by sampleSets.cfg. Likewise, a sampleSets_<versioning>.cfg is constructed

treemakerdir = "/uscms/home/jhiltb/nobackup/susy/ZeroAndTwoLep/CMSSW_10_6_29_patch1/src/TreeMaker/WeightProducer/python"

prod         = "V20"
basedir      = "/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2Production%s/"%(prod)
filelistsdir = "/eos/uscms/store/user/lpcsusyhad/StealthStop/filelists_Kevin_%s/"%(prod)
ttreedir     = "TreeMaker2/PreSelection"
tempfilename = "tmp.txt"
fdir         = "filelists_Kevin_%s/"%prod

if not os.path.isdir(fdir):
    os.mkdir(fdir)

auxiliary = {}

# Extract event numbers from MCSample files in the TreeMaker area
for year in ["16", "16APV", "17", "18"]:
    if os.path.exists(treemakerdir + "/MCSamples_Summer20UL%s.py"%(year)):
        countsFile = open(treemakerdir + "/MCSamples_Summer20UL%s.py"%(year))
        
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

            name = makeNiceName(process)
            if name not in auxiliary.keys():
                auxiliary[name] = {}

            if ndiff == -1.0:
                auxiliary[name]["20%s_npos"%(year)] = float(ntot)
                auxiliary[name]["20%s_nneg"%(year)] = 0.0
            else:
                auxiliary[name]["20%s_npos"%(year)] = (float(ntot) + float(ndiff))/2.0
                auxiliary[name]["20%s_nneg"%(year)] = (float(ntot) - float(ndiff))/2.0

# If there is a MCSamplesValues, use it to get info on xsec, kfactor, etc
if os.path.exists(treemakerdir + "/MCSampleValues.py"):
    xsecFile = open(treemakerdir + "/MCSampleValues.py")
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
                auxiliary = addInfo(auxiliary, process, {"xsec" : xsec, "brf" : brf, "kfactor" : kfactor})
                process = ""; xsec = "-1.0"; brf = "1.0"; kfactor = "1.0"
    
            process = line.split("\"")[1] 
    
        if "XS_" in line:
            xsec = line.split("=")[1].split(",")[0]
    
        if "BR_" in line:
            brf = line.split("=")[1].split(",")[0] 
    
        if "kFactor_" in line:
            kfactor = line.split("=")[1].split(",")[0]

samples = OrderedDict()
samples["2016"]    = defaultdict(list)
samples["2016APV"] = defaultdict(list)
samples["2017"]    = defaultdict(list)
samples["2018"]    = defaultdict(list)

command = os.system("eos root://cmseos.fnal.gov ls %s > %s" % (basedir, tempfilename))
with open(tempfilename, 'r') as tempfile:
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
        newline = "root://cmseos.fnal.gov/" + basedir + line
        samples[era][shortline].append(newline)

sampleSet = open("sampleSets_%s.txt"%(prod), "w")

# Write out a "2018_ttHJetTobb.txt" file with the list of corresponding files
for era, sampleLists in samples.iteritems():

    # For some sanity, sort the samples alphabetically
    theSamples = sampleLists.keys()
    theSamples.sort()
    for sample in theSamples:
        newfile = open("filelists_Kevin_%s/"%(prod) + sample + ".txt", 'w')

        xsec       = "-1.0,"
        nposevents = "-1.0,"
        nnegevents = "-1.0,"
        kfactor    = "-1.0"

        # Everything before "TuneCP5" will be included in name
        # "ttHJetTobb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8" ==> "ttHJetTobb_M125"
        name    = makeNiceName(sample)
        auxName = name.partition("_")[-1]

        name   += ","
        sample += ".txt,"

        # Get any xsec, kfactor info that was retrieved from TreeMaker
        if auxName in auxiliary.keys():
            xsec       = str(float(auxiliary[auxName]["xsec"]) * eval(auxiliary[auxName]["brf"])) + ","
            kfactor    = str(auxiliary[auxName]["kfactor"])

            if "%s_npos"%(era) in auxiliary[auxName]:
                nposevents = str(auxiliary[auxName]["%s_npos"%(era)]) + ","
                nnegevents = str(auxiliary[auxName]["%s_nneg"%(era)]) + ","
               
        sampleSet.write("%s %s, %s %s, %s %s %s %s\n"%(name.ljust(40), filelistsdir, sample.ljust(85), ttreedir, xsec.rjust(14), nposevents.rjust(12), nnegevents.rjust(12), kfactor.rjust(6)))

        for f in sampleLists[sample]:
            newfile.write(f)
        newfile.close()

    sampleSet.write("\n")

sampleSet.close()
