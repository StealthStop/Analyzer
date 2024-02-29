import sys, os
from os import system, environ

from samples import SampleCollection
import optparse 
import subprocess

def red(string):
     CRED = "\033[91m"
     CEND = "\033[0m"
     return CRED + str(string) + CEND

def removeCopies(x):
  return list(dict.fromkeys(x))

def makeExeAndFriendsTarball(filestoTransfer, fname, path):
    system("mkdir -p %s" % fname)
    for fn in removeCopies(filestoTransfer):
        system("cd %s; ln -s %s" % (fname, fn))
        
    tarallinputs = "tar czf %s/%s.tar.gz %s --dereference"% (path, fname, fname)
    system(tarallinputs)
    system("rm -r %s" % fname)

def getTopTaggerTrainingFile(topTaggerFile):
    name = ""
    with file(topTaggerFile) as meowttcfgFile:
        for line in meowttcfgFile:
            if "modelFile" in line:
                name = line.split("=")[1].strip().strip("\"")
                break
    return name

def main():
    repo = "Analyzer/Analyzer"    
    # Parse command line arguments
    parser = optparse.OptionParser("usage: %prog [options]\n")    
    parser.add_option ('-n',        dest='numfile',  type='int',                         default = 10,            help="number of files per job")
    parser.add_option ('-d',        dest='datasets', type='string',                      default = '',            help="List of datasets, comma separated")
    parser.add_option ('-l',        dest='dataCollections',         action='store_true', default = False,         help="List all datacollections")
    parser.add_option ('-L',        dest='dataCollectionslong',     action='store_true', default = False,         help="List all datacollections and sub collections")
    parser.add_option ('-c',        dest='noSubmit',                action='store_true', default = False,         help="Do not submit jobs.  Only create condor_submit.txt.")
    parser.add_option ('-s',        dest='fastMode',                action='store_true', default = False,         help="Run Analyzer in fast mode")
    parser.add_option ('-u',        dest='userOverride',  type='string',                 default = '',            help="Override username with something else")
    parser.add_option ('--output',  dest='outPath',  type='string',                      default = '.',           help="Name of directory where output of each condor job goes")
    parser.add_option ('--analyze', dest='analyze',                                      default = 'Analyze1Lep', help="AnalyzeBackground, AnalyzeEventSelection, Analyze0Lep, Analyze1Lep, MakeNJetDists")    
    options, args = parser.parse_args()

    filesPerJobSkim = {"TT"                            :    10,
                       "QCD_HT50to100"                 :   200,
                       "QCD_HT100to200"                :   100,
                       "QCD_HT200to300"                :    50,
                       "QCD_HT300to500"                :    20,
                       "QCD_HT500to700"                :     5,
                       "QCD_HT700to1000"               :     5,
                       "QCD_HT1000to1500"              :     5,
                       "QCD_HT1500to2000"              :     5,
                       "QCD_HT2000toInf"               :     3,
                       "WJetsToQQ_HT-200to400"         :   200,
                       "WJetsToQQ_HT-400to600"         :     8,
                       "WJetsToQQ_HT-600to800"         :     3,
                       "WJetsToQQ_HT-800toInf"         :     2,
                       "WJetsToLNu_Incl"               :   200,
                       "WJetsToLNu_HT-70To100"         :   200,
                       "WJetsToLNu_HT-100To200"        :   200,
                       "WJetsToLNu_HT-200To400"        :   100,
                       "WJetsToLNu_HT-400To600"        :    10,
                       "WJetsToLNu_HT-600To800"        :     8,
                       "WJetsToLNu_HT-800To1200"       :     5,
                       "WJetsToLNu_HT-1200To2500"      :     3,
                       "WJetsToLNu_HT-2500ToInf"       :     3,
                       "DYJetsToLL_M-50_Incl"          :   200,
                       "DYJetsToLL_M-50_HT-70to100"    :   200,
                       "DYJetsToLL_M-50_HT-100to200"   :   200,
                       "DYJetsToLL_M-50_HT-200to400"   :    50,
                       "DYJetsToLL_M-50_HT-400to600"   :     8,
                       "DYJetsToLL_M-50_HT-600to800"   :     5,
                       "DYJetsToLL_M-50_HT-800to1200"  :     4,
                       "DYJetsToLL_M-50_HT-1200to2500" :     3,
                       "DYJetsToLL_M-50_HT-2500toInf"  :     3,
                       "Diboson"                       :    90,
                       "Triboson"                      :    10,
                       "ST"                            :    20,
                       "TTX"                           :     2,
                       "JetHT"                         :     5,
                       "SingleMuon"                    :    75,
                       "SingleElectron"                :    50,
    }

    filesPerJobSkim = {"Data_SingleMuon" : 40,
                       "Data_JetHT"      : 10,
                       "TT"              : 10,
                       "TTX"             : 4,
                       "ST"              : 6,
                       "Diboson"         : 20,
                       "Triboson"        : 20,
                       "WJets"           : 20,
                       "DYJetsToLL_M-50" : 20,
    }

    srcDir   = environ["CMSSW_BASE"] + "/src"
    testDir  = environ["CMSSW_BASE"] + "/src/%s/test"%(repo) 
    userName = environ["USER"]

    if options.userOverride != "":
        userName = options.userOverride

    hostName = environ["HOSTNAME"]
    cmsConnect = "uscms.org" in hostName

    redirector = "root://cmseos.fnal.gov/"
    workingDir = options.outPath
    eosDir     = "%s/store/user/%s/StealthStop/%s"%(redirector, userName, options.outPath)

    # Prepare the list of files to transfer
    mvaFileName2016preVFP  = getTopTaggerTrainingFile(environ["CMSSW_BASE"] + "/src/%s/test/TopTaggerCfg_2016preVFP.cfg" % repo)
    mvaFileName2016postVFP = getTopTaggerTrainingFile(environ["CMSSW_BASE"] + "/src/%s/test/TopTaggerCfg_2016postVFP.cfg" % repo)
    mvaFileName2017        = getTopTaggerTrainingFile(environ["CMSSW_BASE"] + "/src/%s/test/TopTaggerCfg_2017.cfg" % repo)
    mvaFileName2018        = getTopTaggerTrainingFile(environ["CMSSW_BASE"] + "/src/%s/test/TopTaggerCfg_2018.cfg" % repo)

    filestoTransfer = [testDir + "/MyAnalysis", 
                       testDir + "/%s" % (mvaFileName2016preVFP),
                       testDir + "/%s" % (mvaFileName2016postVFP),
                       testDir + "/%s" % (mvaFileName2017),
                       testDir + "/%s" % (mvaFileName2018),
                       testDir + "/TopTaggerCfg_2016preVFP.cfg",
                       testDir + "/TopTaggerCfg_2016postVFP.cfg",
                       testDir + "/TopTaggerCfg_2017.cfg",
                       testDir + "/TopTaggerCfg_2018.cfg",
                       srcDir  + "/TopTagger/TopTagger/test/libTopTagger.so",
                       testDir + "/sampleSets.cfg",
                       testDir + "/sampleCollections.cfg",
                       testDir + "/Keras_Tensorflow_DoubleDisCo_Reg_0l_RPV_Run2.cfg",
                       testDir + "/Keras_Tensorflow_DoubleDisCo_Reg_1l_RPV_Run2.cfg",
                       testDir + "/Keras_Tensorflow_DoubleDisCo_Reg_2l_RPV_Run2.cfg",
                       testDir + "/Keras_Tensorflow_DoubleDisCo_Reg_0l_SYY_Run2.cfg",
                       testDir + "/Keras_Tensorflow_DoubleDisCo_Reg_1l_SYY_Run2.cfg",
                       testDir + "/Keras_Tensorflow_DoubleDisCo_Reg_2l_SYY_Run2.cfg",
                       testDir + "/Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_0l_RPV_Run2.cfg",
                       testDir + "/Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_1l_RPV_Run2.cfg",
                       testDir + "/Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_2l_RPV_Run2.cfg",
                       testDir + "/Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_0l_SYY_Run2.cfg",
                       testDir + "/Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_1l_SYY_Run2.cfg",
                       testDir + "/Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_2l_SYY_Run2.cfg",
                       testDir + "/keras_frozen_DoubleDisCo_Reg_0l_RPV_Run2.pb",
                       testDir + "/keras_frozen_DoubleDisCo_Reg_1l_RPV_Run2.pb",
                       testDir + "/keras_frozen_DoubleDisCo_Reg_2l_RPV_Run2.pb",
                       testDir + "/keras_frozen_DoubleDisCo_Reg_0l_SYY_Run2.pb",
                       testDir + "/keras_frozen_DoubleDisCo_Reg_1l_SYY_Run2.pb",
                       testDir + "/keras_frozen_DoubleDisCo_Reg_2l_SYY_Run2.pb",
                       testDir + "/allInOne_BTagEff_UL.root",
                       testDir + "/allInOne_SFMean_UL.root",
                       testDir + "/allInOne_leptonicSF_UL.root",
                       testDir + "/allInOne_hadronicSF_UL.root",
                       testDir + "/allInOne_TopTagEffandSF_UL.root",
                       testDir + "/wp_deepJet_106XUL16postVFP_v3.csv",
                       testDir + "/wp_deepJet_106XUL16preVFP_v2.csv",
                       testDir + "/wp_deepJet_106XUL17_v3.csv",
                       testDir + "/wp_deepJet_106XUL18_v2.csv",
                       testDir + "/filelists",
                       # Holding a place for using the reshaping btag SFs if we need them later
                       #testDir + "/reshaping_deepJet_106XUL16postVFP_v3.csv",
                       #testDir + "/reshaping_deepJet_106XUL16preVFP_v2.csv",
                       #testDir + "/reshaping_deepJet_106XUL17_v3.csv",
                       #testDir + "/reshaping_deepJet_106XUL18_v2.csv",
                       ]
    
    print "--------------Files to Transfer-----------------"
    for i in filestoTransfer:    
        print i
    print "------------------------------------------------"
    
    sc = SampleCollection(testDir + "/sampleSets.cfg", testDir + "/sampleCollections.cfg")
    if options.dataCollections or options.dataCollectionslong:
        scl = sc.sampleCollectionList()
        for sampleCollection in scl:
            sl = sc.sampleList(sampleCollection)
            print sampleCollection
            if options.dataCollectionslong:
                sys.stdout.write("\t")
                for sample in sl:
                    sys.stdout.write("%s  "%sample[1])
                print ""
                print ""
        exit(0)
    
    datasets = []
    if options.datasets:
        datasets = options.datasets.split(',')
    else:
        print "No dataset specified"
        exit(0)
    
    fileParts = []
    fileParts.append("Universe             = vanilla\n")
    fileParts.append("Executable           = run_Analyzer_condor.sh\n")
    fileParts.append("Transfer_Input_Files = %s/%s.tar.gz, %s/exestuff.tar.gz\n" % (options.outPath,environ["CMSSW_VERSION"],options.outPath))
    fileParts.append("Request_Memory       = 2.5 Gb\n")
    fileParts.append("x509userproxy        = $ENV(X509_USER_PROXY)\n\n")

    nFilesPerJob = options.numfile
    numberOfJobs = 0
    for ds in datasets:
        ds = ds.strip()

        stubDir = "output-files/%s"%(ds)
        logsDir = "log-files/%s"%(ds)
        # create the directory
        if not os.path.isdir("%s/%s" %(workingDir, logsDir)):
            system('mkdir -p %s/%s' %(workingDir, logsDir))
        else:
            print red("Job directory \"%s/%s\" already exists and cannot proceed safely ! Exiting..."%(workingDir, logsDir))
            exit(0)
   
        dataSetName = ds.partition("_")[-1]
        for s, n, e in sc.sampleList(ds):

            # When running skim jobs with the MiniTreeMaker analyzer,
            # The number of files per job is custom tuned based on skimming efficiency
            # Custom numbers for each main collection are in filesPerJobSkim
            if options.analyze == "MakeMiniTree" or options.analyze == "MakeTopTagSFTree":
                proc = n.partition("_")[-1]
        
                if dataSetName in filesPerJobSkim:
                    nFilesPerJob = filesPerJobSkim[dataSetName]
                elif proc in filesPerJobSkim:
                    nFilesPerJob = filesPerJobSkim[proc]
                else:
                    nFilesPerJob = options.numfile

            print "SampleSet:", n, ", nEvents:", e
            f = open(environ["CMSSW_BASE"] + "/src/Analyzer/Analyzer/test/" + s)
            if not f == None:
                count = 0
                for l in f:
                    if '.root' in l:
                        count = count + 1
                for startFileNum in xrange(0, count, nFilesPerJob):
                    numberOfJobs+=1

                    fileParts.append("Arguments = %s %i %i %s %s %s %s %d %d\n"%(n, nFilesPerJob, startFileNum, s, options.analyze, environ["CMSSW_VERSION"], eosDir + "/" + stubDir, cmsConnect, options.fastMode))
                    fileParts.append("Output    = %s/%s/MyAnalysis_%s_%i.stdout\n"%(workingDir, logsDir, n, startFileNum))
                    fileParts.append("Error     = %s/%s/MyAnalysis_%s_%i.stderr\n"%(workingDir, logsDir, n, startFileNum))
                    fileParts.append("Log       = %s/%s/MyAnalysis_%s_%i.log\n"%(workingDir,    logsDir, n, startFileNum))
                    fileParts.append("Queue\n\n")
    
                f.close()
    
    fout = open("condor_submit.txt", "w")
    fout.write(''.join(fileParts))
    fout.close()

    if not options.dataCollections and not options.dataCollectionslong:
        makeExeAndFriendsTarball(filestoTransfer, "exestuff", options.outPath)
        system("tar --exclude-caches-all --exclude-vcs -zcf %s/${CMSSW_VERSION}.tar.gz -C ${CMSSW_BASE}/.. ${CMSSW_VERSION} --exclude=src --exclude=tmp" % options.outPath)
        
    if not options.noSubmit: 
        system("echo 'condor_submit condor_submit.txt'")
        system('condor_submit condor_submit.txt')
    else:
        print "------------------------------------------"
        print "Number of Jobs:", numberOfJobs
        print "------------------------------------------"

if __name__ == "__main__":
    main()
