import argparse, glob, re

def naturalSort(unsortedList):
    convert      = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    
    return sorted(unsortedList, key=alphanum_key)

# Run the script
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--sampleSet", dest="sampleSet", help="Path to sample set file",  type=str, default="../sampleSets.cfg")
    parser.add_argument("--nEvtDir",   dest="nEvtDir",   help="Directory to nEvt output", type=str, default="./nEvtsOutput")

    args = parser.parse_args()

    obsCounts  = {}
    predCounts = {}
    
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

        for line in lines:
            chunks = line.split(",") 
   
            npos = chunks[-2].split(": ")[-1]
            nneg = chunks[-1].split(": ")[-1]
    
            obsCounts[sampleName] = {"npos" : int(npos), "nneg" : int(nneg)}
         
    tmp = open(args.sampleSet, "r")
    lines = tmp.readlines()
    tmp.close()

    # From the original sampleSet file, get the expected
    # number of positive and negative event counts for each sample
    for line in lines:

        if "#" in line or line == "\n":
            continue
    
        chunks = (line.replace(" ", "")).split(",")

        sampleName = chunks[0]

        npos = chunks[-3]
        nneg = chunks[-2]

        predCounts[sampleName] = {"npos" : int(npos), "nneg" : int(nneg)}

    # Do the comparison between observed and expected
    for sample in naturalSort(obsCounts.keys()):

        if sample in predCounts:

            obsPos = obsCounts[sample]["npos"]
            obsNeg = obsCounts[sample]["nneg"]
    
            predPos = predCounts[sample]["npos"]
            predNeg = predCounts[sample]["nneg"]

            if obsPos != predPos:
                print("%s: Expected %d positive weight events but measured %d"%(sample.ljust(40), predCounts[sample]["npos"], obsCounts[sample]["npos"]))

            if obsNeg != predNeg:
                print("%s: Expected %d negative weight events but measured %d"%(sample.ljust(40), predCounts[sample]["nneg"], obsCounts[sample]["nneg"]))
