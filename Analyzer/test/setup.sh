#!/bin/bash
cmsenv

SCRIPTSDIR=${CMSSW_BASE}/src/Framework/Framework/scripts
if [[ $PATH != *"${SCRIPTSDIR}"* ]]
then
    export PATH=${SCRIPTSDIR}:${PATH}
    echo "adding ${SCRIPTSDIR} to the path"
fi

# set up the top tagger libraries etc
echo "|-----------------------------|"
echo "|      Set up top tagger      |"
echo "|-----------------------------|"
echo ""
source ${CMSSW_BASE}/src/TopTagger/TopTagger/test/taggerSetup.sh
echo "sourced taggerSetup"

# get most up to date samples config
if [ ! -f sampleSets.cfg ] 
then
    echo ""
    echo "|--------------------------------------|"
    echo "|     Soft linking the sample configs  |"
    echo "|--------------------------------------|"
    getSamplesCfg.sh
fi

# Check needed SF files and ensure local file matches with EOS file
# Only copy down from EOS if necessary
sfPath="root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/ScaleFactorHistograms/FullRun2_UL"
sfFiles=("wp_deepCSV_106XUL16preVFP_v2.csv" "wp_deepCSV_106XUL16postVFP_v3.csv" "wp_deepCSV_106XUL17_v3.csv" "wp_deepCSV_106XUL18_v2.csv" "allInOne_leptonSF_UL.root" "allInOne_hadronicSF_UL.root" "allInOne_BTagEff_UL.root" "allInOne_SFMean_UL.root")
for i in ${!sfFiles[@]};
do
    sfFile=${sfFiles[$i]}

    if [[ ! -f $sfFile ]]
    then
        echo "Copying SF file: $sfFile"
        xrdcp -f $sfPath/$sfFile .
    fi

    localChkSum=`xrdadler32 $sfPath/$sfFile`
    eosChkSum=`xrdadler32 $sfFile` 

    localChkSum=${localChkSum% *}
    eosChkSum=${eosChkSum% *}

    if [[ $localChkSum != $eosChkSum ]]
    then
        echo "Updating SF file: $sfFile"
        xrdcp -f $sfPath/$sfFile .
    fi 
done

# Check repos for updates
if [[ "$1" == "-s" ]] 
then
    echo ""
    echo "|--------------------------------------|"
    echo "|      Checking repos for updates      |"
    echo "|--------------------------------------|"
    echo "If it asks for your password too many times you can do something like the following:"
    echo "         eval `ssh-agent`  "
    echo "         ssh-add ~/.ssh/id_rsa"
    status.sh
fi
