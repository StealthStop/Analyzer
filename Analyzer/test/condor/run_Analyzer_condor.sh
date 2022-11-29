#!/bin/bash

#--- You should copy this to a scratch disk where you will manage your job submission and output
#    and prepare the two tar files as described below.

##--- See this web page: http://uscms.org/uscms_at_work/computing/setup/condor_worker_node.shtml

##--- See also make_condor_jdl_files.c which makes jdl files that go with this executable script.

dataset=$1
nfiles=$2
startfile=$3
filelist=$4
analyzer=$5
CMSSW_VERSION=$6
eosPath=$7
cmsConnect=$8
fastMode=$9
base_dir=`pwd`
printf "\n\n base dir is $base_dir\n\n"

source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc820

printf "\n\n ls output\n"
ls -l

#--- You need to make a tar file of your CMSSW directory that contains the
#      TopTagger stuff and your compiled code from the Exploration directory.
#      You may need to change some symbolic links that reference full path
#      names to relative path names in order for the untarred file to work
#      on the batch node.  I forgot which specific changes I made.
#      Once you have done that, make the tar file with something like this
#
#        tar -cvf cmssw-toptagger.tar CMSSW_8_0_28
#

printf "\n\n unpacking CMSSW tar file.\n\n"
tar -xf ${CMSSW_VERSION}.tar.gz

printf "\n\n ls output\n"
ls -l

printf "\n\n changing to CMSSW_X_X_X/ dir\n"
cd ${CMSSW_VERSION}/
mkdir -p src
cd src
scram b ProjectRename
eval `scramv1 runtime -sh`

printf "\n\n ls output\n"
ls -l

printf "\n\n output of uname -s : "
uname -s
printf "\n\n"

cp ${base_dir}/exestuff.tar.gz .
tar xzvf exestuff.tar.gz
cd exestuff/

export LD_LIBRARY_PATH=${PWD}:${LD_LIBRARY_PATH}

printf "\n\n LD_LIBRARY_PATH: ${LD_LIBRARY_PATH}\n\n"
printf "\n\n ls output\n"
ls -l

fastModeStr=""
if [[ $fastMode -eq 1 ]]
then
    fastModeStr="-s"
fi

printf "\n\n Attempting to run MyAnalysis executable.\n\n"
./MyAnalysis -A ${analyzer} --condor -D ${dataset} -N ${nfiles} -M ${startfile} ${fastModeStr}

printf "\n\n ls output\n"
ls -l

HOST=$(hostname)

for outfile in MyAnalysis*.root;
do
    if [[ ( "$HOST" == *"fnal.gov"* ) && ( "$cmsConnect" -eq 1 ) ]]; then
        export WEBDAV_ENDPOINT="davs://cmseos.fnal.gov:9000/eos/uscms/store/user/"
        export OUTDIR=${WEBDAV_ENDPOINT}${eosPath#root://cmseos.fnal.gov//store/user/}
        echo "env -i X509_USER_PROXY=${X509_USER_PROXY} gfal-copy -f -p $outfile $OUTDIR/$outfile"
        env -i X509_USER_PROXY=${X509_USER_PROXY} gfal-copy -f -p $outfile $OUTDIR/$outfile
    else
        echo "xrdcp -f $outfile $eosPath/$outfile"
        xrdcp -f $outfile $eosPath/$outfile
    fi
done

rm -rf MyAnalysis*.root
cd ${base_dir}

printf "\n\n ls output\n"
ls -l
