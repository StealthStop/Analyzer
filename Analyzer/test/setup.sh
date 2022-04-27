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

# get scale factor root files (Should fix this)
if [ ! -f allInOne_BTagEff_UL.root ] 
then
    echo ""
    echo "|--------------------------------------|"
    echo "|  Copying scale factor files          |"
    echo "|--------------------------------------|"
    xrdcp -f root://cmseos.fnal.gov//store/user/bcrossma/StealthStop/ScaleFactorHistograms/FullRun2/reshaping_deepCSV_106XUL16postVFP_v3.csv .
    xrdcp -f root://cmseos.fnal.gov//store/user/bcrossma/StealthStop/ScaleFactorHistograms/FullRun2/reshaping_deepCSV_106XUL16preVFP_v2.csv .
    xrdcp -f root://cmseos.fnal.gov//store/user/bcrossma/StealthStop/ScaleFactorHistograms/FullRun2/reshaping_deepCSV_106XUL17_v3.csv .
    xrdcp -f root://cmseos.fnal.gov//store/user/bcrossma/StealthStop/ScaleFactorHistograms/FullRun2/reshaping_deepCSV_106XUL18_v2.csv .
    xrdcp -f root://cmseos.fnal.gov//store/user/bcrossma/StealthStop/ScaleFactorHistograms/FullRun2/allInOne_leptonSF_UL.root .
    xrdcp -f root://cmseos.fnal.gov//store/user/bcrossma/StealthStop/ScaleFactorHistograms/FullRun2/allInOne_BTagEff_UL.root .
    xrdcp -f root://cmseos.fnal.gov//store/user/bcrossma/StealthStop/ScaleFactorHistograms/FullRun2/allInOne_SFMean_UL.root .
fi

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
