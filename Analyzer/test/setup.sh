#!/bin/bash
cmsenv

redirector="root://cmseos.fnal.gov/"

SCRIPTSDIR=${CMSSW_BASE}/src/Framework/Framework/scripts
if [[ $PATH != *"${SCRIPTSDIR}"* ]]
then
    export PATH=${SCRIPTSDIR}:${PATH}
    printf "Adding ${SCRIPTSDIR} to the path\n"
fi

PYTHONDIR=${CMSSW_BASE}/src/Analyzer/Analyzer/test/Plotters/General
if [[ $PYTHONPATH != *"${PYTHONDIR}"* ]]
then
    export PYTHONPATH=${PYTHONDIR}:${PYTHONPATH}
fi

# set up the top tagger libraries etc
printf "Set up top tagger\n"
printf ""
source ${CMSSW_BASE}/src/TopTagger/TopTagger/test/taggerSetup.sh
printf "%s> sourced taggerSetup\n" "--"

# get most up to date samples config
if [ ! -f sampleSets.cfg ] 
then
    printf "\n"
    printf "Soft linking the sample configs\n"
    getSamplesCfg.sh
fi

# Check needed SF files and ensure local file matches with EOS file
# Only copy down from EOS if necessary
printf "Checking for updated SF ROOT files\n"
sfPath="$redirector/store/user/lpcsusystealth/StealthStop/ScaleFactorHistograms/FullRun2_UL"
sfFiles=("wp_deepJet_106XUL16preVFP_v2.csv" "wp_deepJet_106XUL16postVFP_v3.csv" "wp_deepJet_106XUL17_v3.csv" "wp_deepJet_106XUL18_v2.csv" "allInOne_leptonicSF_UL.root" "allInOne_hadronicSF_UL.root" "allInOne_BTagEff_UL.root" "allInOne_SFMean_UL.root" "allInOne_TopTagEffandSF_UL.root")
for i in ${!sfFiles[@]};
do
    sfFile=${sfFiles[$i]}

    if [[ ! -f $sfFile ]]
    then
        printf "%s> Copying SF file: $sfFile\n" "--"
        xrdcp -f $sfPath/$sfFile .
    fi

    localChkSum=`xrdadler32 $sfPath/$sfFile`
    eosChkSum=`xrdadler32 $sfFile` 

    localChkSum=${localChkSum% *}
    eosChkSum=${eosChkSum% *}

    if [[ $localChkSum != $eosChkSum ]]
    then
        printf "%s> Updating SF file: $sfFile\n" "--"
        xrdcp -f $sfPath/$sfFile .
    fi 
done

fileListPath="/store/user/lpcsusystealth/StealthStop"
fileLists="filelists_Kevin_V20_2"
localFileLists="filelists"
if [[ ! -d $localFileLists ]]
then
    printf "First-time copying of file lists\n"

    mkdir $localFileLists
    xrdcp -s -r --parallel 4 --streams 2 $redirector/$fileListPath/$fileLists $localFileLists
else

    printf "Checking for updated file lists\n"
    
    localLastUpdate=`find $localFileLists -type f -exec stat -c "%Y %n" {} + | sort -n | tail -n 1 | cut -d ' ' -f 1` &> /dev/null
    eosLastUpdate=`find /eos/uscms/$fileListPath/$fileLists -type f -exec stat -c "%Y %n" {} + | sort -n | tail -n 1 | cut -d ' ' -f 1` &> /dev/null
    
    wait
    
    if [[ $localLastUpdate -lt $eosLastUpdate ]]
    then
        printf "%s> Scanning for detected filelist changes\n" "--"
    
        for file in `eosls $fileListPath/$fileLists/`
        do
        
            localChkSum=`xrdadler32 $localFileLists/$file` && localChkSum=${localChkSum% *} &> /dev/null
            eosChkSum=`xrdadler32 $redirector/$fileListPath/$fileLists/$file` && eosChkSum=${eosChkSum% *} &> /dev/null
        
            wait
           
            if [[ $localChkSum != $eosChkSum ]]
            then
                printf "    %s> Updating filelist: $file\n" "--"
                xrdcp -f -s $redirector/$fileListPath/$fileLists/$file $localFileLists/
            fi 
        done
    fi
fi

# Check repos for updates
if [[ "$1" == "-s" ]] 
then
    printf "\n"
    printf "Checking repos for updates\n"
    printf "If it asks for your password too many times you can do something like the following:\n"
    printf "         eval `ssh-agent`  \n"
    printf "         ssh-add ~/.ssh/id_rsa\n"
    status.sh
fi
