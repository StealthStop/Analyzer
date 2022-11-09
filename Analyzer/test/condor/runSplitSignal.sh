#! /bin/bash

INPUTFILE=$1 # full path to a ROOT file on EOS and properly uses the xrd address
OUTPUTDIR=$2 # full path to user's EOS area where output ROOT files will be copied
TREENAME=$3  # name of the TTree inside the input ROOT file to read

# Get a ROOT without using CMSSW
source /cvmfs/sft.cern.ch/lcg/views/LCG_99/x86_64-centos7-gcc10-opt/setup.sh

python runSplitSignal.py --inputFile $INPUTFILE --ttreePath $TREENAME

# Copy output files to EOS and remove from job node
for FOLDER in *mStop* 
do
    for FILE in $FOLDER/*
    do
        xrdcp -f $FILE $OUTPUTDIR/$FOLDER/
    done
done

rm -rf *mStop*
