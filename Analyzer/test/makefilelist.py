import os
from collections import defaultdict, OrderedDict

# Usage: Populates a directory "filelists_Kevin_<versioning>" with a text file for each year + sample combination.
#        Contents of the text file are simply a listing of paths to all relevant ROOT files for the year + sample
#        The "filelists_Kevin_<versioning>" should be placed in /eos/uscms/store/user/lpcsusyhad/StealthStop/
#        where it will be referenced by sampleSets.cfg

treemakerdir = "/uscms/home/jhiltb/nobackup/susy/ZeroAndTwoLep/CMSSW_10_6_29_patch1/src/TreeMaker/WeightProducer/python"

prod         = "V20"
basedir      = "/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2Production%s/"%(prod)
filelistsdir = "/eos/uscms/store/user/lpcsusyhad/StealthStop/filelists_Kevin_%s/"%(prod)
ttreedir     = "TreeMaker2/PreSelection"
tempfilename = "tmp.txt"
fdir         = "filelists_Kevin_%s/"%prod

if not os.path.isdir(fdir):
    os.mkdir(fdir)

samples = OrderedDict()
samples["2016"]    = defaultdict(list)
samples["2016APV"] = defaultdict(list)
samples["2017"]    = defaultdict(list)
samples["2018"]    = defaultdict(list)

command = os.system("eos root://cmseos.fnal.gov ls %s > %s" % (basedir, tempfilename))
with open(tempfilename, 'r') as tempfile:
    for line in tempfile:

        # Each line will be simple ROOT file name e.g. Summer20UL18.ttHJetTobb_M125_TuneCP5_13TeV_amcatnloFXFX_madspin_pythia8_9_RA2AnalysisTree.root
        if not ".root" in line: continue
        if     "Fast"  in line: continue

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

sampleSet = open("sampleSet_%s.txt"%(prod), "w")

# Write out a "2018_ttHJetTobb.txt" file with the list of corresponding files
for era, sampleLists in samples.iteritems():

    # For some sanity, sort the samples alphabetically
    theSamples = sampleLists.keys()
    theSamples.sort()
    for sample in theSamples:
        newfile = open("filelists_Kevin_%s/"%(prod) + sample + ".txt", 'w')

        xsec     = "-1.0,"
        nevents  = "-1.0,"
        nnevents = "-1.0,"
        kfactor  = "-1.0,"

        chunks = sample.split("_")

        endpoint = 0
        for chunk in chunks:
            if "TuneCP5" in chunk:
                if chunk != "TuneCP5":
                    endpoint += 1
                break
    
            endpoint += 1

        name = "_".join(chunks[0:endpoint]) + ","
        sample += ".txt,"
                
        sampleSet.write("%s %s, %s %s, %s %s %s %s\n"%(name.ljust(40), filelistsdir, sample.ljust(85), ttreedir, xsec.ljust(10), nevents.ljust(10), nnevents.ljust(10), kfactor.ljust(10)))

        for f in sampleLists[sample]:
            newfile.write(f)
        newfile.close()

    sampleSet.write("\n")

sampleSet.close()
