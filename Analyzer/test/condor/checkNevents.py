import argparse, glob

# Run the script
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--sampleSet", dest="sampleSet", help="Path to sample set file",  type=str, default="../sampleSets.cfg")
    parser.add_argument("--nEvtDir",   dest="nEvtDir",   help="Directory to nEvt output", type=str, default="./nEvtsOutput")

    args = parser.parse_args()

    # Hold positive and negative event counts
    # measured by nEvts.py
    obsCounts  = {}
    
    # Get all the text files corresponding to each sample
    sampleOutputs = glob.glob(args.nEvtDir + "/output-files/*")

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
   
            npos = chunks[-2].split(": ")[-1]
            nneg = chunks[-1].split(": ")[-1]
    
            obsCounts[sampleName] = {"npos" : int(npos), "nneg" : int(nneg)}
         
    # Read in the sampleSets.cfg file
    tmp = open(args.sampleSet, "r")
    sampleSetLines = tmp.readlines()
    tmp.close()

    # Also, prepare a new sampleSet file with any corrections needed
    newSampleSet = open(args.sampleSet + ".new", "w")

    # From the original sampleSet file, get the expected
    # number of positive and negative event counts for each sample
    for line in sampleSetLines:

        # Write comment lines and newlines without processing
        if "#" in line or line == "\n":
            newSampleSet.write(line)
            continue
    
        chunks = (line.replace(" ", "")).split(",")
        sampleName = chunks[0]
        predPos    = chunks[-3]
        predNeg    = chunks[-2]

        # If sample was not run over with nEvts.py
        # just put back the line in the new sampleSet file
        if sampleName not in obsCounts:
            newSampleSet.write(line)
            continue

        obsPos = obsCounts[sampleName]["npos"]
        obsNeg = obsCounts[sampleName]["nneg"]

        posDiffLen = len(predPos) - len(str(obsPos))
        negDiffLen = len(predNeg) - len(str(obsNeg))

        # Make corrections to the line if the nEvts do not match up
        # Do some fancy string stuff to preserve the nice formatting in sampleSet.cfg
        newLine = line
        if obsPos != int(predPos):
            print("%s: Expected %s positive weight events but measured %d"%(sampleName.ljust(40), npos, obsCounts[sampleName]["npos"]))
            if posDiffLen >= 0:
                newLine = newLine.replace(" "+str(predPos)+",", " "*(posDiffLen+1)+str(obsPos)+",")
            else:
                newLine = newLine.replace(" "*abs(posDiffLen)+str(predPos)+",", str(obsPos)+",")

        if obsNeg != int(nneg):
            print("%s: Expected %s negative weight events but measured %d"%(sampleName.ljust(40), nneg, obsCounts[sampleName]["nneg"]))
            if negDiffLen >= 0:
                newLine = newLine.replace(" "+str(predNeg)+",", " "*(negDiffLen+1)+str(obsNeg)+",")
            else:
                newLine = newLine.replace(" "*abs(negDiffLen)+str(predNeg)+",", str(obsNeg)+",")
       
        newSampleSet.write(newLine)
    newSampleSet.close()
