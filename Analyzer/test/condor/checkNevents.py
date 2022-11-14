import argparse, glob

# Run the script
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--sampleSet", dest="sampleSet", help="Path to sample set file",   type=str, default="../sampleSets.cfg")
    parser.add_argument("--nEvtsDir",  dest="nEvtsDir",  help="Directory to nEvts output", type=str, default="./nEvtsOutput")

    args = parser.parse_args()

    # Hold positive and negative event counts
    # measured by nEvts.py
    obsCounts  = {}
    
    # Get all the text files corresponding to each sample
    sampleOutputs = glob.glob(args.nEvtsDir + "/output-files/*")

    # From the nEvt job output files, save the observed positive and negative
    # event counts
    for sampleFile in sampleOutputs:

        # Text filenames in form: "output_2018_TTToHadronic.txt"
        sampleName = sampleFile.split("/")[-1].split(".txt")[0].split("output_")[-1]
        tmp = open(sampleFile, "r")
        lines = tmp.readlines()
        tmp.close()

        # Parse nEvts.py output containing the string
        # Positive weights: 95170542, Negative weights: 0
        for line in lines:
            chunks = line.split(",") 
   
            npos = chunks[-2].split(": ")[-1].rstrip()
            nneg = chunks[-1].split(": ")[-1].rstrip()
    
            obsCounts[sampleName] = {"npos" : npos, "nneg" : nneg}
         
    # Read in the sampleSets.cfg file
    tmp = open(args.sampleSet, "r")
    sampleSetLines = tmp.readlines()
    tmp.close()

    # Also, prepare a new sampleSet file with any corrections needed
    newSampleSet = open(args.sampleSet.replace(".cfg", "_new.cfg"), "w")

    # From the original sampleSet file, get the expected
    # number of positive and negative event counts for each sample
    for line in sampleSetLines:

        # Write comment lines and newlines without processing
        if "#" in line or line == "\n":
            newSampleSet.write(line)
            continue
    
        chunks = (line.replace(" ", "")).split(",")
        sampleName = chunks[0]
        predPos    = chunks[-3].rstrip()
        predNeg    = chunks[-2].rstrip()

        # If sample was not run over with nEvts.py
        # just put back the line in the new sampleSet file
        if sampleName not in obsCounts:
            newSampleSet.write(line)
            continue

        obsPos = obsCounts[sampleName]["npos"]
        obsNeg = obsCounts[sampleName]["nneg"]

        posDiffLen = len(predPos) - len(obsPos)
        negDiffLen = len(predNeg) - len(obsNeg)

        # Make corrections to the line if the nEvts do not match up
        # Do some fancy string stuff to preserve the nice formatting in sampleSet.cfg
        newLine = line
        if int(obsNeg) != int(predNeg):
            print("%s: Expected %s negative weight events but measured %s"%(sampleName.ljust(40), predNeg, obsNeg))
            old = ""
            new = ""
            if negDiffLen >= 0:
                old = " "+predNeg+","
                new = " "*(negDiffLen+1)+obsNeg+","
            else:
                old = " "*abs(negDiffLen)+predNeg+","
                new = obsNeg+","

            newLine = new.join(newLine.rsplit(old, 1))

        if int(obsPos) != int(predPos):
            print("%s: Expected %s positive weight events but measured %s"%(sampleName.ljust(40), predPos, obsPos))
            old = ""
            new = ""
            if posDiffLen >= 0:
                old = " "+predPos+","
                new = " "*(posDiffLen+1)+obsPos+","
            else:
                old = " "*abs(posDiffLen)+predPos+","
                new = obsPos+","

            newLine = new.join(newLine.rsplit(old, 1))

        newSampleSet.write(newLine)
    newSampleSet.close()
