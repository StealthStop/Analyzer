import ROOT
import subprocess as sb
import argparse
import os
import glob

backgrounds = ["TTToHadronic", "TTToSemiLeptonic", "TTTo2L2Nu"]
variations  = ["", "TuneCP5up", "TuneCP5down", "hdampUP", "hdampDOWN", "erdON"]

models = ["RPV", "StealthSYY"]
masses = list(range(300, 1450, 50))

splits = ["Train", "Val", "Test"]

if __name__ == "__main__":

    usage = "usage: %prog [options]"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("--outputdir",     dest="outputdir", help="Directory to store hadded output",      required=True)
    parser.add_argument("--inputdir",      dest="inputdir",  help="Directory of unhadded input",           required=True)
    parser.add_argument("--year",          dest="year",      help="Which year to hadd",                    required=True)
    parser.add_argument("--dryrun",        dest="dryrun",    help="Print what will happen, don't do it",   action="store_true", default=False)
    parser.add_argument("--verbose", "-v", dest="verbose",   help="Commands are printed before execution", action="store_true", default=False)
    args = parser.parse_args()

    if not os.path.exists(args.outputdir):
        os.makedirs(args.outputdir)

    inputdir = ""
    if args.inputdir[-1] == "/":
        inputdir = args.inputdir[:-1]
    else:
        inputdir = args.inputdir

    # Hadd the background files together
    for background in backgrounds:
        for variation in variations:
            for split in splits:

                varStr = ""
                if variation != "":
                    varStr = "_%s"%(variation)

                files = glob.glob("%s/%s*/*_%s%s_?_%s.root"   %(inputdir,args.year,background,varStr,split)) \
                      + glob.glob("%s/%s*/*_%s%s_??_%s.root"  %(inputdir,args.year,background,varStr,split)) \
                      + glob.glob("%s/%s*/*_%s%s_???_%s.root" %(inputdir,args.year,background,varStr,split)) \
                      + glob.glob("%s/%s*/*_%s%s_????_%s.root"%(inputdir,args.year,background,varStr,split))

                nFilesPerChunk = 20
                nChunks = len(files) / nFilesPerChunk + 1
                if nChunks == 0:
                    nChunks += 1

                for i in range(0, nChunks):

                    filesChunk = []
                    for j in range(i*nFilesPerChunk, (i+1)*nFilesPerChunk):
                        if j >= len(files):
                            break
                        filesChunk.append(files[j])

                    command = ["hadd", "-f", "%s/MyAnalysis_%s_%s%s_%s_%d.root"%(args.outputdir,args.year,background,varStr,split,i)] + filesChunk

                    if args.verbose or args.dryrun: print("Executing command: \"%s\""%(" ".join(command)))
                    if not args.dryrun:
                        p = sb.Popen(command)
                        p.wait()

                tempFiles = glob.glob("%s/MyAnalysis_%s_%s%s_%s_*.root"%(args.outputdir,args.year,background,varStr,split))
                command = ["hadd", "-f", "%s/MyAnalysis_%s_%s%s_%s.root"%(args.outputdir,args.year,background,varStr,split)] + tempFiles

                if args.verbose or args.dryrun: print("Executing command: \"%s\""%(" ".join(command)))
                if not args.dryrun:
                    p = sb.Popen(command)
                    p.wait()

                for tempFile in tempFiles:
                    os.remove(tempFile) 

    for model in models:
        for mass in masses:
            for split in splits:

                command = ["hadd", "-f", "%s/MyAnalysis_%s_%s_2t6j_mStop-%s_%s.root"%(args.outputdir,args.year,model,mass,split)] \
                        + glob.glob("%s/%s*/*_%s_2t6j_mStop-%s_*_%s.root"%(inputdir,args.year,model,mass,split))

                if args.verbose or args.dryrun: print("Executing command: \"%s\""%(" ".join(command)))
                if not args.dryrun:
                    p = sb.Popen(command)
                    p.wait()
