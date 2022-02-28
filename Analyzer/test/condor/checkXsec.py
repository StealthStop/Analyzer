import argparse, glob

# Run the script
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--xsecDir",   dest="xsecDir",   help="Directory to xsec ana output", type=str, default="./xsec_temp_2016")

    args = parser.parse_args()

    # Get all the text files corresponding to each sample
    sampleOutputs = glob.glob(args.xsecDir + "/log-files/*.stdout")

    # From the xsec analyzer job output files
    for sampleFile in sampleOutputs:

        # Text filenames in form: "MyAnalysis_2018_TTToHadronic.stdout"
        sampleName = sampleFile.split("/")[-1].split(".stdout")[0].split("MyAnalysis_")[-1]

        tmp = open(sampleFile, "r")
        lines = tmp.readlines()
        tmp.close()

        for line in lines:

            if "XSEC INFO" in line:

                chunks = line.split(", ") 
   
                TMxsec = float(chunks[-2].split(":")[-1])
                SSxsec = float(chunks[-1].split(":")[-1])
    
                if TMxsec != SSxsec:
                    print("{}: TreeMaker xsec is {:.6e}, but stealth stop xsec is {:.6e}".format(sampleName.ljust(35), TMxsec, SSxsec))
