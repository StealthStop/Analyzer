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
noUpdates=1
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
        noUpdates=0
    fi 
done
if [[ $noUpdates == 1 ]]
then
    printf "%s> All SF files already up-to-date !\n" "--"
fi

# Check for any changes to file lists
# For efficiency, only scan files if EOS file list dir has a newer mod time than the local folder
# or if there are more files in the EOS file list dir than in the local folder
# If either of these are true, then the individual file lists files are scanned and chk summed
fileListPath="/store/user/lpcsusystealth/StealthStop"
fileLists="filelists_Kevin_V20_2"
localFileLists="filelists"
noUpdates=1
if [[ ! -d $localFileLists ]]
then
    printf "First-time copying of file lists\n"

    mkdir $localFileLists
    xrdcp -s -r --parallel 4 --streams 2 $redirector/$fileListPath/$fileLists $localFileLists
else

    printf "Checking for updated file lists\n"
    
    eosFileCount=`xrdfs root://cmseos.fnal.gov ls -l $fileListPath/$fileLists | wc -l` &> /dev/null
    localFileCount=`ls -1 $localFileLists | wc -l` &> /dev/null

    wait
    
    eosLastUpdate=`xrdfs root://cmseos.fnal.gov ls -l $fileListPath | grep $fileLists | awk '{print $2, $3}' | xargs -I {} date -d "{}" +"%s" --utc` &> /dev/null
    localLastUpdate=`ls -ld --full-time $localFileLists | awk '{print $6, $7, $8}' | xargs -I {} date -d "{}" +"%s"` &> /dev/null
    
    wait

    if [[ $localLastUpdate -lt $eosLastUpdate ]] || [[ $eosFileCount -gt $localFileCount ]]
    then
        noUpdates=0
        printf "%s> Changes detected, scanning individual file lists\n" "--"
    
        for file in `xrdfs root://cmseos.fnal.gov ls $fileListPath/$fileLists/`
        do
        
            localChkSum=`xrdadler32 $localFileLists/${file##*/}` && localChkSum=${localChkSum% *} &> /dev/null
            eosChkSum=`xrdadler32 $redirector/$file` && eosChkSum=${eosChkSum% *} &> /dev/null
        
            wait
           
            if [[ $localChkSum != $eosChkSum ]]
            then
                printf "    %s> Updating filelist: $(basename $file)\n" "--"
                xrdcp -f -s $redirector/$file $localFileLists/
            fi 
        done
        if [[ $eosFileCount -lt $localFileCount ]]
        then
            printf "    (Note: some file lists on EOS appear to have been removed, though no local files are removed)\n"
            touch -m $localFileLists
        fi
    fi
fi
if [[ $noUpdates == 1 ]]
then
    printf "%s> All file lists up-to-date !\n" "--"
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
