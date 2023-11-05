#!/bin/bash

command=$1

YEARS=("2016preVFP" "2016postVFP" "2017" "2018")
#YEARS=("2018")
DATE=("MassExclusion_Fix_11_3_23")
#samples=("Data" "QCD" "TTX" "BG_OTHER" "TT" "TT_erdON" "TT_hdampDOWN" "TT_hdampUP" "TT_TuneCP5down" "TT_TuneCP5up") 
samples=("TT_skim" "TT_erdON_skim" "TT_hdampUP_skim" "TT_hdampDOWN_skim" "TT_TuneCP5up_skim" "TT_TuneCP5down_skim" "QCD_skim" "TTX_skim" "BG_OTHER_skim" "Data_skim")


for year in ${YEARS[@]}; do

    # -------------------
    # top tagger resolved
    # -------------------
    if [ $command == "resolved" ] || [ $command == "all" ]; then
        echo "hadding for ResolvedTopTagger study :---------------------------------" 
        python hadder.py -d ${year}_TT,${year}_QCD -H hadd_ResolvedTopTagger_fakeRateEfficiency_$DATE/ -f -m -p ${year}_ResolvedTopTagger_fakeRateEfficiency_$DATE/output-files 
    fi

    # -----------
    # trigger SFs
    # -----------
    if [ $command == "trigger" ] || [ $command == "all" ]; then
        echo "hadding for hadronic trigger SFs :---------------------------------"
        python hadder.py -d ${year}_TT,${year}_AllSignal,${year}_Data_SingleMuon -H Thesis_AN_FullStatusTalk/hadd_Run2UL_HadronicTriggerEfficiencySF/ -f -m -p ${year}_HadronicTriggerEfficiencySF_$DATE/output-files

        echo "hadding for leptonic trigger SFs :---------------------------------"
        #python hadder.py -d ${year}_QCD,${year}_AllSignal,${year}_Data_SingleMuon,${year}_Data_SingleElectron -H Thesis_AN_FullStatusTalk/hadd_Run2UL_LeptonicTriggerEfficiencySF/ -f -m -p ${year}_LeptonicTriggerEfficiencySF_$DATE/output-files
        #python hadder.py -d ${year}_TT,${year}_Data_SingleMuon,${year}_Data_SingleElectron -H Thesis_AN_FullStatusTalk/hadd_Run2UL_LeptonicTriggerEfficiencySF/ -f -m -p ${year}_LeptonicTriggerEfficiencySF_$DATE/output-files
    fi

    # ---------------
    # HEM issue plots
    # ---------------
    if [ $command == "hem" ] || [ $command == "all" ]; then
        echo "hadding for HEM isuue :---------------------------------"
        if [ ${year} == "2018" ]; then
            python hadder.py -d ${year}_Data_JetHT,${year}_Data_SingleElectron,${year}_Data_SingleMuon -H hadd_${year}_HEMissue_0l_1l_2l_$DATE -p ${year}_HEMissue_0l_1l_2l_$DATE/output-files
        fi
    fi

    # --------------------------------
    # data/MC for checking trigger SFs
    # --------------------------------
    #if [ $command == "Semra" ] || [ $command == "all" ]; then
    #    echo "hadding for data-MC tests :---------------------------------"
    #    python hadder.py -y ${year} -d ${year}_TT,${year}_QCD,${year}_TTX,${year}_DYJetsToLL_M-50,${year}_Diboson,${year}_ST,${year}_Triboson,${year}_WJets,${year}_Data -H TriggerSFsTest_dataMC_2022/hadd_${year}_dataMC_triggersSFsWith-1b-ge2b_0L_1L_$DATE -p ${year}_dataMC_triggersSFsWith-1b-ge2b_0L_1L_$DATE/output-files --haddOther
    #    python hadder.py -y ${year} -d ${year}_TT,${year}_QCD,${year}_TTX,${year}_DYJetsToLL_M-50,${year}_Diboson,${year}_ST,${year}_Triboson,${year}_WJets,${year}_Data -H TriggerSFsTest_dataMC_2022/hadd_${year}_dataMC_triggersSFsWith-2b-3b-ge4b_0L_1L_$DATE -p ${year}_dataMC_triggersSFsWith-2b-3b-ge4b_0L_1L_$DATE/output-files --haddOther
    #fi

    # -----------------
    # validation inputs
    # -----------------
    if [ $command == "disco" ] || [ $command == "all" ]; then
        echo "hadding for validation/higgs inputs :---------------------------------"
        python hadder.py -y ${year} -d ${year}_TT_skim,${year}_TT_erdON_skim,${year}_TT_hdampUP_skim,${year}_TT_hdampDOWN_skim,${year}_TT_TuneCP5up_skim,${year}_TT_TuneCP5down_skim,${year}_QCD_skim,${year}_TTX_skim,${year}_DYJetsToLL_M-50_skim,${year}_Diboson_skim,${year}_ST_skim,${year}_Triboson_skim,${year}_WJets_skim,${year}_AllSignal,${year}_Data_skim -H DisCo_outputs_0l_1l_2l_$DATE/ -f -m -p ${year}_DisCo_outputs_0l_1l_2l_$DATE/output-files --haddOther
        #python hadder.py -y ${year} -d ${year}_AllSignal -H DisCo_outputs_0l_1l_2l_$DATE/ -f -m -p ${year}_DisCo_outputs_0l_1l_2l_$DATE/output-files --haddOther > ${year}_haddDisCo.log
    fi

 
done

# -------------------------
# Merge all years as Run2UL
# -------------------------
if [ $command == "merge" ] || [ $command == "all" ]; then
    echo "merging all years as Run2UL :---------------------------------"

    cd DisCo_outputs_0l_1l_2l_${DATE}/
    
    # background
    for sample in ${samples[@]}; do
        hadd Run2UL_${sample}.root 2016preVFP_${sample}.root 2016postVFP_${sample}.root 2017_${sample}.root 2018_${sample}.root
    done
    
    
    # signal
    for mass in {300..1400..50}; do
        hadd Run2UL_RPV_2t6j_mStop-${mass}.root 2016preVFP_RPV_2t6j_mStop-${mass}.root 2016postVFP_RPV_2t6j_mStop-${mass}.root 2017_RPV_2t6j_mStop-${mass}.root 2018_RPV_2t6j_mStop-${mass}.root
        hadd Run2UL_StealthSYY_2t6j_mStop-${mass}.root 2016preVFP_StealthSYY_2t6j_mStop-${mass}.root 2016postVFP_StealthSYY_2t6j_mStop-${mass}.root 2017_StealthSYY_2t6j_mStop-${mass}.root 2018_StealthSYY_2t6j_mStop-${mass}.root
    done

fi

