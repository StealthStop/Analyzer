import os
import subprocess

from samples import SampleCollection
import argparse
import glob
import shutil
import ROOT

def red(string):
     CRED = "\033[91m"
     CEND = "\033[0m"
     return CRED + str(string) + CEND

def cyan(string):
     CRED = "\033[96m"
     CEND = "\033[0m"
     return CRED + str(string) + CEND

def orange(string):
     CRED = "\033[93m"
     CEND = "\033[0m"
     return CRED + str(string) + CEND

def eosCompatible(query):

    compatible = True
    if "[" in query or "]" in query or "?" in query:
        print cyan("Running in EOS compatiblity mode with \"ls\" rather than \"eosls\"")
        compatible = False

    return compatible

def printError(log):
    if len(log) > 0:
         print red("------------------------------------------------------------------------------------------------")
         print red("There was some jobs that didn't match the epected number of events, see summary below")
         for l in log:
              print red(l)
         print red("------------------------------------------------------------------------------------------------")

def checkNumEvents(nEvents, rootFile, sampleCollection, log):
    try:
         f = ROOT.TFile.Open(rootFile)
         f.cd()
         try:
              h = f.Get("EventCounter")
              nNeg = h.GetBinContent(1)
              nPos = h.GetBinContent(2)
              diff = nEvents-(nPos-nNeg)
              if abs(diff) > 5.0:
                   message = "Error: Sample: "+sampleCollection+" Expected nEvents:  "+str(nEvents)+" EventCounter nEvents: "+str(nPos-nNeg)+" = "+str(nPos)+" "+str(-nNeg)
                   log.append(message)
                   print red("----------------------------------------------------------------------------------------------------------")
                   print red("Num events in \"EventCounter\" doesn't match the number in \"sampleSet.cfg\"")
                   print red(message)
                   print red("----------------------------------------------------------------------------------------------------------")
         except:
              print red("Error: Problem opening and reading from histogram \"EventCounter\"")
              pass
         f.Close()
    except:
         print red("Error: Can't open rootFile: %s" % rootFile)
         pass

    return log

def listFiles(inPath, query):

    files = None
    if eosCompatible(query):
        stub = query.split("/MyAnalysis")[0]
        commandList = ["eos", "root://cmseos.fnal.gov", "ls", query]

        proc = subprocess.Popen(commandList, stdout=subprocess.PIPE)
        payload = proc.stdout.readlines()
        files   = ["/eos/uscms" + stub + "/" + file.rstrip() for file in payload]
    else:
        files = glob.glob("/eos/uscms" + query)

    return files

def getDataSets(inPath):
    l = glob.glob(inPath+"/*")
    print "-------------------------------------------------------" 
    print cyan("Automatically hadding all directory names in input path")
    print "-------------------------------------------------------\n" 
    return list(path.split("/")[-1] for path in l)

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser("usage: %prog [options]\n")
    parser.add_argument('-d', dest='datasets', type=str,            default = '',             help = "Lists of datasets, comma separated")
    parser.add_argument('-H', dest='outDir',   type=str,            default = 'rootfiles',    help = "Can pass in the output directory name")
    parser.add_argument('-p', dest='inPath',   type=str,            default = 'output-files', help = "Can pass in the input directory name")
    parser.add_argument('-y', dest='year',     type=str,            default = '',             help = "Can pass in the year for this data")
    parser.add_argument('-f',                  action='store_true',                           help = "ok to write in existing folder")
    parser.add_argument('-v',                  action='store_true',                           help = "make hadd verbose")
    parser.add_argument('-o',                  action='store_true',                           help = "Overwrite output directory")
    parser.add_argument('-m',                  action='store_true',                           help = "multiprocess, restricted to 4")
    parser.add_argument('--noHadd',            action='store_true',                           help = "Dont hadd the the root files")
    parser.add_argument('--haddOther',         action='store_true',                           help = "Do the hack to make BG_OTHER.root")
    parser.add_argument('--haddData',          action='store_true',                           help = "Do the hack to make Data.root")
    parser.add_argument('--haddAll',           action='store_true',                           help = "Do the hack to hadd All_Bg, All_Signal, All_Data")
    options = parser.parse_args()

    # Get input directory path
    inPath = options.inPath

    userName = os.getenv("USER")

    stubPath = "/store/user/%s/StealthStop"%(userName)
    eosPath = "root://cmseos.fnal.gov//%s/%s"%(stubPath, inPath)
        
    # Checks if user specified a dataset(s)
    datasets = []
    if options.datasets:
        datasets = options.datasets.split(',')
    else:
        datasets = getDataSets("/eos/uscms" + stubPath + "/" + inPath)

    haddArgs = "-v 0"
    if options.v:
        haddArgs = "-v 99"
    if options.m:
        print cyan("Running hadd multithreaded with 4 threads")
        haddArgs += " -j 4"
    
    # Check if output directory exits and makes it if not
    outDir = options.outDir
    overwrite = options.o
    useSameDir = options.f
    if os.path.exists(outDir):
        if overwrite: 
            print red("Warning: Overwriting output directory !")
            shutil.rmtree(outDir)
            os.makedirs(outDir)
            pass
        elif not useSameDir:
            print red("Error: Output directory %s already exits ! Use option -f to OK writing into existing folder." % ('"'+outDir+'"'))
            exit(0)    
        else:
            print orange("User OK'd writing hadded output to existing folder: \"%s\" !"%(outDir))
    else:
        os.makedirs(outDir) 

    # Here we are in the nominal case of hadder.py
    # Combine files within one input directory
    # Loop over all sample options to find files to hadd
    log = []
    sc = SampleCollection("../sampleSets.cfg", "../sampleCollections.cfg")
    scl = sc.sampleCollectionList()
    for sampleCollection in scl:
        sl = sc.sampleList(sampleCollection)
        if sampleCollection in datasets:
            directory = sampleCollection
            files = ""
            print "-----------------------------------------------------------"
            print sampleCollection
            print "-----------------------------------------------------------"
            
            # hadd individual samples within a collection
            # this includes the AllSignal and the AllTT
            sampleSetsToHadd = ["2016preVFP_AllSignal", "2016postVFP_AllSignal", "2017_AllSignal", "2018_AllSignal",
                                "2016preVFP_AllTT",     "2016postVFP_AllTT",     "2017_AllTT",     "2018_AllTT",
                                "2016preVFP_RPV",       "2016postVFP_RPV",       "2017_RPV",       "2018_RPV",
                                "2016preVFP_StealthSYY","2016postVFP_StealthSYY","2017_StealthSYY","2018_StealthSYY",
                                "2016preVFP_TTX",       "2016postVFP_TTX",       "2017_TTX",       "2018_TTX"
                                ]

            if sampleCollection in sampleSetsToHadd:
                individualSamples = ""
                for sample in sl:
                    files = None

                    # For TT and its variations choose one of the three channels e.g. SemiLeptonic
                    # And make wildcard string to hadd together all three channels files together
                    # This is specifically here when choosing the AllTT collection: it should not
                    # be hadded together, rather the individual variations should be hadded separately
                    tempSample = sample[1]
                    cleanName = sample[1]

                    # Use the SemiLepton name to make a wildcard string to pick up all three channel names for hadding together
                    # tempSample will contain the wildcards and cleanName goes into the output filename
                    if "_TTTo" in tempSample and "SemiLep" in tempSample:
                        tempSample = sample[1].replace("_TTToSemiLeptonic", "_TTTo*[cu]")
                        cleanName = sample[1].replace("ToSemiLeptonic", "")
                    # Skip the other two channels for TT since they already were hadded in with the SemiLep file
                    elif "_TTTo" in tempSample:
                        continue

                    query = "%s/%s/%s/MyAnalysis_%s_[0-9]*.root"%(stubPath, inPath, directory, tempSample)

                    files = " ".join(listFiles(inPath, query))

                    outfile = "%s/%s.root" % (outDir, cleanName)
                    command = "hadd %s %s/%s.root %s" % (haddArgs, outDir, cleanName, files)
                    if not options.noHadd: subprocess.call(command.split(" "))
                    log = checkNumEvents(nEvents=float(sample[2]), rootFile=outfile, sampleCollection=cleanName, log=log)

                    individualSamples += " %s/%s.root"%(outDir, cleanName)

                # After hadding all the little guys individually, put them together into the big sampleCollection (if desired)
                if "TTX" in sampleCollection:
                    command = "hadd %s %s/%s.root%s"%(haddArgs, outDir, sampleCollection, individualSamples)
                    if not options.noHadd: subprocess.call(command.split(" "))
    
            # hadd samples that were not produced under and a broad collection name
            # I.e. something like 2017_QCD
            else:
                nEvents=0.0
                files = []
                for sample in sl:
                    query = "%s/%s/%s/MyAnalysis_%s_*.root" % (stubPath, inPath, directory, sample[1])
                    files += listFiles(inPath, query)
                    nEvents+=float(sample[2])
    
                outfile = "%s/%s.root" % (outDir,sampleCollection)
                command = "hadd %s %s %s" % (haddArgs, outfile, " ".join(files))
                try:
                    if not options.noHadd: 
                        process = subprocess.Popen(command, shell=True)
                        process.wait()
                except:
                    print red("Warning: Too many files to hadd, using the exception setup")
                    query = "%s/%s/%s/MyAnalysis_*"%(stubPath, inPath, sampleCollection)
                    files = " ".join(listFiles(inPath, query))
                    
                    command = "hadd %s %s/%s.root %s" % (haddArgs, outDir, sampleCollection, files)
                    if not options.noHadd: os.system(command)
                    pass
    
                log = checkNumEvents(nEvents=nEvents, rootFile=outfile, sampleCollection=sampleCollection, log=log)
   
    if options.haddAll:
        for year in ['2016preVFP', '2016postVFP', '2017', '2018']:
            query_TT = "%s/%s/%s_AllBg/MyAnalysis_%s*TTTo*.root" % (stubPath, inPath, year, year)
            files_TT = " ".join(listFiles(inPath, query_TT)) 
            command = "hadd %s %s/%s_TT.root %s" % (haddArgs, outDir, year, files_TT)
            os.system(command)

            query_QCD = "%s/%s/%s_AllBg/MyAnalysis_%s*QCD*.root" % (stubPath, inPath, year, year)
            files_QCD = " ".join(listFiles(inPath, query_QCD))
            command = "hadd %s %s/%s_QCD.root %s" % (haddArgs, outDir, year, files_QCD)
            os.system(command)

            query_TTX1 = "%s/%s/%s_AllBg/MyAnalysis_%s*TT[WZ][TJ]*.root" % (stubPath, inPath, year, year)
            query_TTX2 = "%s/%s/%s_AllBg/MyAnalysis_%s*ttH*.root" % (stubPath,year, year)
            files_TTX = " ".join(listFiles(inPath, query_TTX1)) + " " + " ".join(listFiles(inPath, query_TTX2))
            command = "hadd %s %s/%s_TTX.root %s" % (haddArgs, outDir, year, files_TTX)
            os.system(command)

            query_DY = "%s/%s/%s_AllBg/MyAnalysis_%s*DYJets*.root" % (stubPath, inPath, year, year)
            query_2B = "%s/%s/%s_AllBg/MyAnalysis_%s_??_*.root" % (stubPath, inPath, year, year)
            query_3B = "%s/%s/%s_AllBg/MyAnalysis_%s_???_*.root" % (stubPath, inPath, year, year)
            query_Wj = "%s/%s/%s_AllBg/MyAnalysis_%s_WJets*.root" % (stubPath, inPath, year, year)

            files_DY = " ".join(listFiles(inPath, query_DY))
            files_2B = " ".join(listFiles(inPath, query_2B))
            files_3B = " ".join(listFiles(inPath, query_3B))
            files_Wj = " ".join(listFiles(inPath, query_Wj))

            command = "hadd %s %s/%s_BG_OTHER.root %s" % (haddArgs, outDir, year, files_DY + " " + files_2B + " " + files_3B + " " + files_Wj)
            os.system(command)

            query_Data = "%s/%s/%s_Data/MyAnalysis_%s*.root" % (stubPath, inPath, year, year)
            files_Data = " ".join(listFiles(inPath, query_Data))
            command = "hadd %s %s/%s_Data.root %s" % (haddArgs, outDir, year, files_Data)
            os.system(command)

    #Print log of hadd at the end
    printError(log)    

    if options.haddOther:
        # Hack to make the BG_OTHER.root file
        other_2016preVFP  = ["2016preVFP_Triboson",  "2016preVFP_Diboson.root",  "2016preVFP_DYJetsToLL_M-50.root",  "2016preVFP_WJets.root"]
        other_2016postVFP = ["2016postVFP_Triboson", "2016postVFP_Diboson.root", "2016postVFP_DYJetsToLL_M-50.root", "2016postVFP_WJets.root"]
        other_2017        = ["2017_Triboson",        "2017_Diboson.root",        "2017_DYJetsToLL_M-50.root",        "2017_WJets.root"]
        other_2018        = ["2018_Triboson",        "2018_Diboson.root",        "2018_DYJetsToLL_M-50.root",        "2018_WJets.root"]
        other = other_2016preVFP + other_2016postVFP + other_2017 + other_2018

        files = ""
        for sampleCollection in scl:
            sl = sc.sampleList(sampleCollection)
            if sampleCollection in datasets and sampleCollection in other: 
                directory = sampleCollection
                files += " %s/%s.root " % (outDir, directory)
        if options.year:
            command = "hadd %s %s/%s_BG_OTHER.root %s" % (haddArgs, outDir, options.year, files)
        else:
            command = "hadd %s %s/BG_OTHER.root %s" % (haddArgs, outDir, files)
        print "-----------------------------------------------------------"
        print command
        print "-----------------------------------------------------------"
        os.system(command)
        
    if options.haddData:
        # Hack to make the Data.root file (hadd all the data together)
        dataFiles = ["Data_SingleMuon.root",             "Data_SingleElectron.root",             "Data_JetHT.root",
                     "2016preVFP_Data_SingleMuon.root",  "2016preVFP_Data_SingleElectron.root",  "2016preVFP_Data_JetHT.root",
                     "2016postVFP_Data_SingleMuon.root", "2016postVFP_Data_SingleElectron.root", "2016postVFP_Data_JetHT.root",
                     "2017_Data_SingleMuon.root",        "2017_Data_SingleElectron.root",        "2017_Data_JetHT.root",
                     "2018_Data_SingleMuon.root",        "2018_Data_SingleElectron.root",        "2018_Data_JetHT.root"]
        if options.year:
            command = "hadd %s %s/%s_Data.root " % (haddArgs, outDir,options.year)
        else:
            command = "hadd %s %s/Data.root " % (haddArgs, outDir)
        for f in dataFiles:
            if os.path.exists(outDir+"/"+f):
                command += " %s/%s" % (outDir, f)
        print "-----------------------------------------------------------"
        print command
        print "-----------------------------------------------------------"
        os.system(command)

if __name__ == "__main__":
    main()
