import os
import glob
import optparse 
import subprocess

from samples import SampleCollection

def red(string):
     CRED = "\033[91m"
     CEND = "\033[0m"
     return CRED + str(string) + CEND

# List things in a specified path (no wildcards) and handle
# whether path is to EOS or just somewhere locally
# If on EOS, then use the proper xrdfs ls command to get the list
# Process the return and format into a clean list of strings
def listPath(path, onEOS=False):

    # Depending on using xrdfs ls or bash ls, the position of the file sizes will be different
    proc = None
    pos = None
    if onEOS:
        proc = subprocess.Popen(["xrdfs", "root://cmseos.fnal.gov", "ls", "-l", path], stdout=subprocess.PIPE)
        pos = 3
    else:
        proc = subprocess.Popen(["ls", "-l", path], stdout=subprocess.PIPE)
        pos = 4

    lines = []
    fileList = proc.stdout.readlines()
    for line in fileList:

        chunks = None
        if ".root" in line:
            raw = line.rstrip("\n").split(" ")
            chunks = filter(str.strip, raw)
            size = chunks[pos]

            if int(size) < 512:
                continue

        # Ensure that each element in the list has its full path
        if not onEOS:
            line = path + "/" + chunks[-1]
        lines.append(line)

    return lines

# Use the list of log files for a particular job folder and the
# corresponding ROOT files in the output directory on EOS to determine
# which jobs did not complete, and also the number of files per job that was run
def getMissingJobs(rootFiles, logFiles):

    # Subfunction to take a list of files and extract the job ID
    def extractJobs(listOfFiles):

        jobIDs = {}
        for f in listOfFiles:

            filename = os.path.basename(f).rpartition(".")[0]
            temp   = filename.partition("_")[-1]
            chunks = temp.rpartition("_")

            sampleset = str(chunks[0])
            jobid     = int(chunks[-1])

            if sampleset not in jobIDs:
                jobIDs[sampleset] = []

            jobIDs[sampleset].append(jobid)

        sortedJobIDs = {}
        for sample, idList in jobIDs.items():
            sortedJobIDs[sample] = sorted(idList)

        return sortedJobIDs

    logDir = os.path.dirname(logFiles[0])
    proc = subprocess.Popen("grep -l TNet %s/*.stderr"%(logDir), stdout=subprocess.PIPE, shell=True)
    errFiles = []
    fileList = proc.stdout.readlines()
    for line in fileList:
        errFiles.append(line.rstrip("\n"))

    finishedJobs = extractJobs(rootFiles)
    totalJobs    = extractJobs(logFiles)
    errorJobs    = extractJobs(errFiles)
      
    nFilesPerJob = {}
    missingJobs = {}
    for sample, allJobIDs in totalJobs.items():
        if len(allJobIDs) > 1:
            nFilesPerJob[sample] = abs(totalJobs[sample][0] - totalJobs[sample][1])
        else:
            nFilesPerJob[sample] = 999999

        # If none of the jobs for a sample ran
        # then put the entire list of job IDs
        if sample not in finishedJobs:
            missingJobs[sample] = allJobIDs
            continue
            
        finishedJobIDs = finishedJobs[sample]
        errorJobIDs    = []
        if sample in errorJobs:
            errorJobIDs = errorJobs[sample]
        for jobID in allJobIDs:
            if jobID not in finishedJobIDs or jobID in errorJobIDs: 

                if sample not in missingJobs:
                    missingJobs[sample] = []

                missingJobs[sample].append(jobID)

    return nFilesPerJob, missingJobs

def main():
    parser = optparse.OptionParser("usage: %prog [options]\n")    
    parser.add_option ('-c',        dest='noSubmit',                action='store_true', default = False,         help="Do not submit jobs.  Only create condor_submit.txt.")
    parser.add_option ('-s',        dest='fastMode',                action='store_true', default = False,         help="Run Analyzer in fast mode")
    parser.add_option ('-u',        dest='userOverride',  type='string',                 default = '',            help="Override username with something else")
    parser.add_option ('--jobdir',  dest='jobdir',        type='string',                 default = '.',           help="Name of directory where output of each condor job goes")
    options, args = parser.parse_args()

    userName = os.environ["USER"]

    if options.userOverride != "":
        userName = options.userOverride

    hostName = os.environ["HOSTNAME"]
    cmsConnect = "uscms.org" in hostName

    redirector = "root://cmseos.fnal.gov/"
    testDir    = os.environ["CMSSW_BASE"] + "/src/Analyzer/Analyzer/test"
    workingDir = options.jobdir
    eosDir     = "/store/user/%s/StealthStop/%s"%(userName, options.jobdir)

    jobSubdirs = [jobSubdir.rstrip("\n").split(" ")[-1] for jobSubdir in listPath(eosDir + "/output-files", onEOS=True)]

    sc = SampleCollection(testDir + "/sampleSets.cfg", testDir + "/sampleCollections.cfg")

    fileParts = []
    fileParts.append("Universe             = vanilla\n")
    fileParts.append("Executable           = run_Analyzer_condor.sh\n")
    fileParts.append("Transfer_Input_Files = %s/%s.tar.gz, %s/exestuff.tar.gz\n" % (options.jobdir,os.environ["CMSSW_VERSION"],options.jobdir))
    fileParts.append("Request_Memory       = 2.5 Gb\n")
    fileParts.append("x509userproxy        = $ENV(X509_USER_PROXY)\n\n")

    nJob = 0
    for jobSubdir in jobSubdirs:
        
        # The last part of the job subdir path is actually the corresponding sample collection name
        collection = jobSubdir.split("/")[-1]

        stubDir = "output-files/%s"%(collection)
        logsDir = "log-files/%s"%(collection)

        rootFiles = listPath(jobSubdir, onEOS=True)
        logFiles  = glob.glob(workingDir + "/" + logsDir + "/*.log")

        outFile = logFiles[0].replace(".log", ".stdout")
        analyzer = None
        with open(outFile) as f:
            for line in f:
                if "Running the" == line[0:11]:
                    analyzer = line.split(" ")[2]
                    break

        nFilesPerJob, missingJobs = getMissingJobs(rootFiles, logFiles)

        for sample, missingJobIDs in missingJobs.items():

            # Need to search the sample list to get the correct
            # filelist to go with that sample
            filelist = ""
            for f, s, _ in sc.sampleList(collection):
                if sample == s:
                    filelist = f
                    break

            for missingJobID in missingJobIDs:

                fileParts.append("Arguments = %s %i %i %s %s %s %s %d %d\n"%(sample, nFilesPerJob[sample], missingJobID, filelist, analyzer, os.environ["CMSSW_VERSION"], redirector + eosDir + "/" + stubDir, cmsConnect, options.fastMode))
                fileParts.append("Output    = %s/%s/MyAnalysis_%s_%i.stdout\n"%(workingDir, logsDir, sample, missingJobID))
                fileParts.append("Error     = %s/%s/MyAnalysis_%s_%i.stderr\n"%(workingDir, logsDir, sample, missingJobID))
                fileParts.append("Log       = %s/%s/MyAnalysis_%s_%i.log\n"%(workingDir,    logsDir, sample, missingJobID))
                fileParts.append("Queue\n\n")
                nJob += 1
            
    fout = open("condor_submit.txt", "w")
    fout.write(''.join(fileParts))
    fout.close()

    # No need to try and submit if there are no missing jobs
    if not options.noSubmit and nJob > 0: 
        os.system("echo 'condor_submit condor_submit.txt'")
        os.system('condor_submit condor_submit.txt')
    else:
        print "------------------------------------------"
        print "Number of Jobs:", nJob
        print "------------------------------------------"

if __name__ == "__main__":
    main()
