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

    proc = None
    if onEOS:
        proc = subprocess.Popen(["xrdfs", "root://cmseos.fnal.gov", "ls", path], stdout=subprocess.PIPE)
    else:
        proc = subprocess.Popen(["ls", path], stdout=subprocess.PIPE)

    payload = proc.stdout.readlines()

    lines = []
    for line in payload:
        line = str(line).rstrip("\n")

        # Check the size, if clearly a broken file (< 512 bytes) do not count as "done"
        if ".root" in line:

            proc = subprocess.Popen(["xrdfs", "root://cmseos.fnal.gov", "stat", line], stdout=subprocess.PIPE)
            info = proc.stdout.readlines()
            size = None
            for stat in info:
                if "Size" in stat:
                    size = stat.partition("Size:")[-1].replace(" ", "").rstrip()
                    break
                else:
                    continue
                    
            if int(size) < 512:
                continue

        # Ensure that each element in the list has its full path
        if not onEOS:
            line = path + "/" + line
        lines.append(line)

    return lines

# Use the list of log files for a particular job folder and the
# corresponding ROOT files in the output directory on EOS to determine
# which jobs did not complete, and also the number of files per job that was run
def getMissingJobs(rootFiles, logFiles):

    # Subfunction to take a list of files and extract the job ID
    def extractJobs(listOfFiles, checkError = False):

        jobIDs = {}
        for f in listOfFiles:

            filename = os.path.basename(f).rpartition(".")[0]
            temp   = filename.partition("_")[-1]
            chunks = temp.rpartition("_")

            sampleset = str(chunks[0])
            jobid     = int(chunks[-1])

            if sampleset not in jobIDs:
                jobIDs[sampleset] = []

            if checkError:
                try:
                    subprocess.check_output("grep \"TNet\" %s"%(os.path.abspath(f).replace(".log", ".stderr")), shell=True)
                except:
                    jobIDs[sampleset].append(jobid)
            else:
                jobIDs[sampleset].append(jobid)

        sortedJobIDs = {}
        for sample, idList in jobIDs.items():
            sortedJobIDs[sample] = sorted(idList)

        return sortedJobIDs

    finishedJobs = extractJobs(rootFiles)
    totalJobs    = extractJobs(logFiles)
    errFreeJobs  = extractJobs(logFiles, checkError=True)
      
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
        errFreeJobIDs  = errFreeJobs[sample]
        for jobID in allJobIDs:
            if jobID not in finishedJobIDs or jobID not in errFreeJobIDs: 

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
    parser.add_option ('--analyze', dest='analyze',                                      default = 'Analyze1Lep', help="AnalyzeBackground, AnalyzeEventSelection, Analyze0Lep, Analyze1Lep, MakeNJetDists")    
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

    jobSubdirs = listPath(eosDir + "/output-files", onEOS=True)

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

                fileParts.append("Arguments = %s %i %i %s %s %s %s %d %d\n"%(sample, nFilesPerJob[sample], missingJobID, filelist, options.analyze, os.environ["CMSSW_VERSION"], redirector + eosDir + "/" + stubDir, cmsConnect, options.fastMode))
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
