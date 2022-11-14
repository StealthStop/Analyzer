from samples import SampleSet
import optparse 
import os
import fnmatch
from os import system, environ

def makeExeAndFriendsTarrball(filestoTransfer, path):
    system("mkdir -p exestuff")
    for fn in filestoTransfer:
        system("cd exestuff; ln -s %s" % fn)
    
    tarallinputs = "tar czf %s/exestuff.tar.gz exestuff --dereference"%(path)
    system(tarallinputs)
    system("rm -r exestuff")

parser = optparse.OptionParser("usage: %prog [options]\n")
parser.add_option ('--noSubmit',   dest='noSubmit',   action='store_true', default = False,         help="Do not submit jobs. Only create condor_submit.txt")
parser.add_option ('--outPath',    dest='outPath',    type='string',       default = 'nEvtsOutput', help="Name of directory where output of each condor job goes")
parser.add_option ('--sampleSets', dest='sampleSets', type='string',       default = "sampleSets",  help="Sample sets config file")
parser.add_option ('--wildcard',   dest='wildcard',   type='string',       default = "*",           help="Wildcard expression for picking only some sample sets")

options, args = parser.parse_args()
sampleSets = options.sampleSets

repo = "Analyzer/Analyzer"    

testDir = environ["CMSSW_BASE"] + "/src/%s/test"%(repo)

#Here is the configuration for the Data/MC validation of the TopTagger 
filestoTransfer  = [
                         testDir + "/condor/nEvts.py",
                         testDir + "/%s.cfg"%(sampleSets),
                         testDir + "/condor/samples.py",
                         testDir + "/filelists",
                         testDir + "/obj/samplesModule.so"
                   ]

# make directory for condor submission
if not os.path.isdir("%s/output-files" % (options.outPath)):
    os.makedirs("%s/output-files" % (options.outPath))

if not os.path.isdir("%s/log-files" % (options.outPath)):
    os.makedirs("%s/log-files" % (options.outPath))

makeExeAndFriendsTarrball(filestoTransfer, options.outPath)
system("tar --exclude-caches-all --exclude-vcs -zcf %s/${CMSSW_VERSION}.tar.gz -C ${CMSSW_BASE}/.. ${CMSSW_VERSION} --exclude=src --exclude=tmp"%(options.outPath))

# Write out condor submit file
fout = open("condor_submit.txt", "w")
fout.write("universe              = vanilla\n")
fout.write("Executable            = goNEvts.sh\n")
fout.write("Requirements          = OpSys == \"LINUX\" && (Arch != \"DUMMY\")\n")
fout.write("Should_Transfer_Files = YES\n")
fout.write("WhenToTransferOutput  = ON_EXIT\n")
fout.write("Transfer_Input_Files  = %s/%s.tar.gz, %s/exestuff.tar.gz\n" % (options.outPath,environ["CMSSW_VERSION"],options.outPath))
fout.write("x509userproxy         = $ENV(X509_USER_PROXY)\n\n")

ss = SampleSet(environ["CMSSW_BASE"] + "/src/%s/test/%s.cfg"%(repo,sampleSets))
for ds in ss.sampleSetList():
    dsn = ds[0]

    if not fnmatch.fnmatch(dsn, options.wildcard):
        continue

    fout.write("Arguments = %s %s %s.cfg\n"%(dsn, environ["CMSSW_VERSION"], sampleSets))
    fout.write("transfer_output_remaps = \"output_%s.txt = %s/output-files/output_%s.txt\"\n"%(dsn, options.outPath, dsn))
    fout.write("Output    = %s/log-files/nEvts_%s.stdout\n"%(options.outPath, dsn))
    fout.write("Error     = %s/log-files/nEvts_%s.stderr\n"%(options.outPath, dsn))
    fout.write("Log       = %s/log-files/nEvts_%s.log\n"%(options.outPath, dsn))
    fout.write("Queue\n\n")

fout.close()

if not options.noSubmit: 
    system('mkdir -p logs')
    system("echo 'condor_submit condor_submit.txt'")
    system('condor_submit condor_submit.txt')

print "Sample sets file: %s"%(sampleSets)
print "Submission directory: %s"%(options.outPath)
