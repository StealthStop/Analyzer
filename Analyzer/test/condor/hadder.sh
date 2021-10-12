#!/bin/bash

dir=${CMSSW_BASE}/src/Analyzer/Analyzer/test/condor/output-files
finaldir=${CMSSW_BASE}/src/Analyzer/Analyzer/test/condor/final

function strip {
    local STRING=${1#$"$2"}
    echo ${STRING%$"$2"}
}

# Make sure directory ends with "/*"
if [[ $dir != */ ]]
then
    d="$dir/*"
else
    d="$dir*"
fi

echo "${dir}"

# Loop over output-files
for f in $d
do
    # Only interested in directories
    [ -d "${f}" ] || continue

    # Get the name of the sample collection
    sc=$(strip "$f" "$dir/")
    
    if [ "${sc}" = "AllSignal" ] 
    then 
        continue
    fi

    #echo "hadd ${sc}.root *.root"
    hadd -fk ${finaldir}/${sc}.root ${f}/*.root
    

    #if [ "${sc}" = "Rare" ]
    #then
    #    #echo "hadd ${sc}.root *.root"
    #    hadd ${f}/${sc}.root ${f}/*.root
    #else
    #    #echo "hadd ${sc}.root *${sc}*.root"
    #    hadd ${f}/${sc}.root ${f}/*${sc}*.root
    #fi

done
