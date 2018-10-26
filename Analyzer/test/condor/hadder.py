import sys, os
from os import system, environ
sys.path = [environ["CMSSW_BASE"] + "/src/SusyAnaTools/Tools/condor/",] + sys.path

from samples import SampleCollection
import optparse
from glob import glob
import datetime
import shutil

# Parse command line arguments
parser = optparse.OptionParser("usage: %prog [options]\n")
parser.add_option('-d', dest='datasets', type='string', default='',          help="Lists of datasets, comma separated")
parser.add_option('-H', dest='outDir',   type='string', default='rootfiles', help="Give the output directory name")
parser.add_option('-o', action='store_true',                                 help="Overwrite output directory")
options, args = parser.parse_args()

# Checks if user specified a dataset(s)
datasets = []
if options.datasets:
    datasets = options.datasets.split(',')
else:
    print "Failed: No dataset specified"
    exit(0)

# Check if output directory exits and makes it if not
outDir = options.outDir
overwrite = options.o
if os.path.exists(outDir):
    if overwrite: 
        print "Overwriting output directory"
        shutil.rmtree(outDir)
        os.makedirs(outDir)
    else:
        print "Failed: Output directory %s already exits" % ('"'+outDir+'"')
        exit(0)

else:
    os.makedirs(outDir) 

sc = SampleCollection("../sampleSets.cfg", "../sampleCollections.cfg")
scl = sc.sampleCollectionList()
for sampleCollection in scl:
    sl = sc.sampleList(sampleCollection)
    if sampleCollection in datasets:
        directory = sampleCollection
        rootName = sampleCollection
        files = ""
        print "-----------------------------------------------------------"
        print sampleCollection
        print "-----------------------------------------------------------"
        
        # copy signal root files
        if sampleCollection == "AllSignal":
            for sample in sl:
                command = "cp output-files/%s/MyAnalysis_%s_0.root %s/%s.root" % (directory, sample[1], outDir, sample[1])
                print command
                system(command)

        # hadd other condor jobs
        else:
            for sample in sl:
                files += " " + " ".join(glob("output-files/%s/MyAnalysis_%s_*.root" % (directory, sample[1])))
            command = "hadd %s/%s.root %s" % (outDir, sampleCollection, files)
            #print command
            system(command)

## Hack to make the BG_noTT.root file
#print outDir
#print overwrite
