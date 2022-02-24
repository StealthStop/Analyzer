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
if [ ! -f allInOne_BTagEff.root ] 
then
    echo ""
    echo "|--------------------------------------|"
    echo "|  Copying scale factor files          |"
    echo "|--------------------------------------|"
    xrdcp -f root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/ScaleFactorHistograms/SUS-19-004_Final/CSVv2_Moriond17_B_H.csv .

    xrdcp -f root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/ScaleFactorHistograms/FullRun2/DeepCSV_102XSF_WP_V1.csv .
    xrdcp -f root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/ScaleFactorHistograms/FullRun2/DeepCSV_2016LegacySF_WP_V1.csv .
    xrdcp -f root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/ScaleFactorHistograms/FullRun2/DeepCSV_94XSF_WP_V4_B_F.csv .
    xrdcp -f root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/ScaleFactorHistograms/FullRun2/allInOne_leptonSF_2016.root .
    xrdcp -f root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/ScaleFactorHistograms/FullRun2/allInOne_leptonSF_2017.root .
    xrdcp -f root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/ScaleFactorHistograms/FullRun2/allInOne_leptonSF_2018.root .
    xrdcp -f root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/ScaleFactorHistograms/FullRun2/allInOne_BTagEff.root .
    xrdcp -f root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/ScaleFactorHistograms/FullRun2/allInOne_SFMean.root .
    xrdcp -f root://cmseos.fnal.gov//store/user/lpcsusyhad/StealthStop/ScaleFactorHistograms/FullRun2/L1prefiring_jetpt_2017BtoF.root .
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
