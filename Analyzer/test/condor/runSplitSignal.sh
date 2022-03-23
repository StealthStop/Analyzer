#! /bin/bash

INPUTFILE=$1
OUTPUTDIR=$2
TREENAME=$3

source /cvmfs/sft.cern.ch/lcg/views/LCG_99/x86_64-centos7-gcc10-opt/setup.sh

python runSplitSignal.py --inputFile $INPUTFILE --ttreeName $TREENAME

for FILE in *.root
do
    xrdcp -r $FILE $OUTPUTDIR
    rm $FILE
done
