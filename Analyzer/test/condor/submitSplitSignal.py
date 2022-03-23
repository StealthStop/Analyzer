import ROOT, subprocess, os, argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--eosPath"  , dest="eosPath"  , help="Path to files on EOS"   , type=str,  required=True                     )
    parser.add_argument("--outPath"  , dest="outPath"  , help="Output path for jobs"   , type=str,  default="SplitJobs"               )
    parser.add_argument("--treeName" , dest="treeName" , help="TTree name to read"     , type=str,  default="TreeMaker2/PreSelection" )
    parser.add_argument("--wildcard" , dest="wildcard" , help="Wildcards for picking"  , type=str,  default="*300to1400*"             )
    parser.add_argument("--noSubmit" , dest="noSubmit" , help="Do not submit to condor",            default=False, action="store_true")
    args = parser.parse_args()

    USER = os.getenv("USER")
    outputDir = "/eos/uscms/store/user/%s/%s"%(USER, args.outPath)

    # Make directory for condor submission logs and output
    if not os.path.isdir("%s/log-files" % (args.outPath)):
        os.makedirs("%s/log-files" % (args.outPath))
    
    if not os.path.isdir("/eos/uscms/store/user/%s/%s"%(USER, args.outPath)):
        os.makedirs("/eos/uscms/store/user/%s/%s"%(USER, args.outPath))

    outputDir = outputDir.replace("/eos/uscms/", "root://cmseos.fnal.gov///")

    # Write out condor submit file
    fout = open("condor_submit.txt", "w")
    fout.write("universe              = vanilla\n")
    fout.write("Executable            = runSplitSignal.sh\n")
    fout.write("Requirements          = OpSys == \"LINUX\" && (Arch != \"DUMMY\")\n")
    fout.write("Transfer_Input_Files  = runSplitSignal.py\n")
    fout.write("Should_Transfer_Files = YES\n")
    fout.write("WhenToTransferOutput  = ON_EXIT\n")
    fout.write("x509userproxy         = $ENV(X509_USER_PROXY)\n\n")

    # Use EOS commands to get the list of files at the eosPath that contain the wildcard phrase
    proc     = subprocess.Popen(["eos", "root://cmseos.fnal.gov", "ls", "%s/%s.root"%(args.eosPath, args.wildcard)], stdout=subprocess.PIPE)
    eosls    = proc.stdout.readlines()
    fileList = [line.decode("utf8").strip() for line in eosls]

    iFile = 0
    for file in fileList:
        fout.write("Arguments = %s/%s %s %s\n"%(args.eosPath.replace("/eos/uscms/", "root://cmseos.fnal.gov///"), file, outputDir, args.treeName))
        fout.write("Output    = %s/log-files/sigSplit_%s.stdout\n"%(args.outPath, iFile))
        fout.write("Error     = %s/log-files/sigSplit_%s.stderr\n"%(args.outPath, iFile))
        fout.write("Log       = %s/log-files/sigSplit_%s.log\n"%(args.outPath, iFile))
        fout.write("Queue\n\n")
    
    fout.close()

    if not args.noSubmit: 
        os.system('mkdir -p logs')
        os.system("echo 'condor_submit condor_submit.txt'")
        os.system('condor_submit condor_submit.txt')
