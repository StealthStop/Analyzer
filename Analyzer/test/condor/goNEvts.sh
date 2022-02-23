#!/bin/bash

dataset=$1
CMSSW_VERSION=$2
sampleset=$3

_PWD=${PWD}

printenv

echo ""

ls 

echo ""

export PATH=${PATH}:/cvmfs/cms.cern.ch/common
export CMS_PATH=/cvmfs/cms.cern.ch

#move the files from /home to /var/tmp
if [ ! -f ${_CONDOR_SCRATCH_DIR}/${CMSSW_VERSION}.tar.gz ]; then
    mv ${_PWD}/${CMSSW_VERSION}.tar.gz ${_CONDOR_SCRATCH_DIR}
fi
if [ ! -f ${_CONDOR_SCRATCH_DIR}/exestuff.tar.gz ]; then
    mv ${_PWD}/exestuff.tar.gz ${_CONDOR_SCRATCH_DIR}
fi

cd ${_CONDOR_SCRATCH_DIR}

#get the release setup and in place
tar -xzf ${CMSSW_VERSION}.tar.gz
cd ${CMSSW_VERSION} 
mkdir -p src
cd src
scram b ProjectRename
eval `scramv1 runtime -sh`

#set up local code
tar -xzf ${_CONDOR_SCRATCH_DIR}/exestuff.tar.gz
cd exestuff

ls -lrht

mkdir -p ../Analyzer/Analyzer/test/obj
mv samplesModule.so ../Analyzer/Analyzer/test/obj

pwd
ls -lhrt

python nEvts.py -s ${sampleset} -d "^$1$" > output_${dataset}.txt

ls -lhrt

mv output_$1.txt ${_CONDOR_SCRATCH_DIR}

rm -r ${_CONDOR_SCRATCH_DIR}/$2
