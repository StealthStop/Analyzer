import sys, os
from os import system, environ
import subprocess
import optparse
from glob import glob
import datetime
import shutil
import ROOT

def red(string):
     CRED = "\033[91m"
     CEND = "\033[0m"
     return CRED + str(string) + CEND

def main():
    # Parse command line arguments
    parser = optparse.OptionParser("usage: %prog [options]\n")
    parser.add_option('-H', dest='outDir',   type='string', default='rootfiles',    help="Can pass in the output directory name")
    parser.add_option('-p', dest='inPath',   type='string', default='output-files', help="Can pass in the input directory name")
    parser.add_option('-o', action='store_true',                                    help="Overwrite output directory")
    options, args = parser.parse_args()

    # Get input directory path
    inPath = options.inPath
            
    # Check if output directory exits and makes it if not
    outDir = options.outDir
    overwrite = options.o
    if os.path.exists(outDir):
        if overwrite: 
            print red("Warning: Overwriting output directory")
            shutil.rmtree(outDir)
            os.makedirs(outDir)
        else:
            print red("Error: Output directory %s already exits" % ('"'+outDir+'"'))
            exit(0)    
    else:
        os.makedirs(outDir) 

    # Figure out the files to hadd
    files = glob("%s/20*_*.root" % (inPath))
    sampleSet = set((f.split("/")[1]).split("_", 1)[1] for f in files)
    samplesToHadd = dict((s,[]) for s in sampleSet)
    for s in sampleSet:
        #print s
        for i, f in enumerate(files):
            if s in f:
                samplesToHadd[s].append(f)

    # Loop over all samples and run hadd command
    for key in samplesToHadd:
        #print key, len(samplesToHadd[key])
        command = "hadd %s/Run2_%s" % (outDir,key)
        for f in samplesToHadd[key]:
            command += " %s" % (f)
        print "-----------------------------------------------------------"
        print command
        print "-----------------------------------------------------------"
        system(command)

if __name__ == "__main__":
    main()
