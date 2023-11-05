#!/bin/bash

voms-proxy-init --voms cms --hours 168


command=$1

#YEARS=("2016preVFP" "2016postVFP" "2017" "2018")
YEARS=("2016preVFP" "2016postVFP" "2017" "2018")
DATE=("MaxSig_Fix_11_05_23")


for year in ${YEARS[@]}; do

    # -------------------
    # top tagger resolved
    # -------------------
    if [ $command == "resolved" ] || [ $command == "all" ]; then
        echo "running ResolvedTopTagger_Analyzer :---------------------------------" 
        python condorSubmit.py --analyze ResolvedTopTagger_Analyzer -d ${year}_TT,${year}_QCD  -n 20 --output ${year}_ResolvedTopTagger_fakeRateEfficiency_withDeepCSV_$DATE
    fi

    # -----------
    # trigger SFs
    # -----------
    if [ $command == "trigger" ] || [ $command == "all" ]; then
        echo "running HadTriggers_Analyzer :---------------------------------"
        python condorSubmit.py --analyze HadTriggers_Analyzer -d ${year}_TT,${year}_AllSignal,${year}_Data_SingleMuon -n 40 --output ${year}_HadronicTriggerEfficiencySF_$DATE

        echo "running AnalyzeLepTrigger :---------------------------------"
        #python condorSubmit.py --analyze AnalyzeLepTrigger -d ${year}_QCD,${year}_AllSignal,${year}_Data_SingleMuon,${year}_Data_SingleElectron -n 40 --output ${year}_LeptonicTriggerEfficiencySF_$DATE
        #python condorSubmit.py --analyze AnalyzeLepTrigger -d ${year}_TT,${year}_Data_SingleMuon,${year}_Data_SingleElectron -n 40 --output ${year}_LeptonicTriggerEfficiencySF_$DATE
    fi
 
    # ---------
    # HEM issue
    # ---------
    if [ $command == "hem" ] || [ $command == "all" ]; then
        echo "running HEM_Analyzer :---------------------------------"
        if [ ${year} == "2018" ]; then
            python condorSubmit.py --analyze HEM_Analyzer -d ${year}_Data_JetHT,${year}_Data_SingleElectron,${year}_Data_SingleMuon -n 20 --output ${year}_HEMissue_0l_1l_2l_$DATE 
        fi
    fi

    # --------------------------------
    # data/MC for checking trigger SFs
    # --------------------------------
    #if [ $command == "Semra" ] || [ $command == "all" ]; then
    #    echo "running Semra_Analyzer :---------------------------------"
    #    python condorSubmit.py --analyze Semra_Analyzer -d ${year}_TT,${year}_QCD,${year}_TTX,${year}_DYJetsToLL_M-50,${year}_Diboson,${year}_ST,${year}_Triboson,${year}_WJets,${year}_Data -n 20 --output ${year}_dataMC_triggersSFsWith-1b-ge2b_0L_1L_$DATE
    #    python condorSubmit.py --analyze Semra_Analyzer -d ${year}_TT,${year}_QCD,${year}_TTX,${year}_DYJetsToLL_M-50,${year}_Diboson,${year}_ST,${year}_Triboson,${year}_WJets,${year}_Data -n 20 --output ${year}_dataMC_triggersSFsWith-2b-3b-ge4b_0L_1L_$DATE
    #fi

    # -----------------
    # validation inputs
    # -----------------
    if [ $command == "disco" ] || [ $command == "all" ]; then
        echo "running AnalyzeDoubleDisCo :---------------------------------"
        python condorSubmit.py --analyze AnalyzeDoubleDisCo -d ${year}_TT_skim,${year}_TT_erdON_skim,${year}_TT_hdampUP_skim,${year}_TT_hdampDOWN_skim,${year}_TT_TuneCP5up_skim,${year}_TT_TuneCP5down_skim,${year}_QCD_skim,${year}_TTX_skim,${year}_DYJetsToLL_M-50_skim,${year}_Diboson_skim,${year}_ST_skim,${year}_Triboson_skim,${year}_WJets_skim,${year}_AllSignal,${year}_Data_skim -n 4 --output ${year}_DisCo_outputs_0l_1l_2l_$DATE -s

    fi


    # -------------
    # make minitree
    # -------------
    if [ $command == "minitree" ] || [ $command == "all" ]; then
        echo "running MakeMiniTree :---------------------------------"
        python condorSubmit.py --analyze MakeMiniTree -d ${year}_TT -n 10 -u lpcsusystealth --output Skims/${year} 

    fi

done
