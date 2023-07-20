###########################################
# Making all closure correction ratio plots
###########################################

command=$1

YEARS=("Run2UL") #"2016preVFP" "2016postVFP" "2017" "2018")
MODELS=("RPV" "SYY")
MASSES=("550")
CHANNELS=("0l" "1l" "2l")


# ------------------------------------------
# MC-based systematics
#   -- Closure Correction Ratio = ttvar / tt
# -----------------------------------------
for MODEL in ${MODELS[@]}; do

        for MASS in ${MASSES[@]}; do

            #################################
            if [ $command == "MassExclusion" ] || [ $command == "All" ]; then

                echo "running for RPV MassExclusion : ---------------------------------"
                python plotting_MCcorrRatio.py --sig RPV --outpath Run2UL_MassExclusion_RPV --mass ${MASS} --div MassExclusion
                
                echo "running for SYY MassExclusion : ---------------------------------" 
                python plotting_MCcorrRatio.py --sig SYY --outpath Run2UL_MassExclusion_StealthSYY --mass ${MASS} --div MassExclusion

            fi

            #################################
            if [ $command == "MaxSign" ] || [ $command == "All" ]; then

                echo "running for RPV MaxSign : ---------------------------------"
                python plotting_MCcorrRatio.py --sig RPV --outpath Run2UL_MaxSign_RPV --mass ${MASS} --div MaxSign

                echo "running for SYY MaxSign : ---------------------------------" 
                python plotting_MCcorrRatio.py --sig SYY --outpath Run2UL_MaxSign_StealthSYY --mass ${MASS} --div MaxSign

            fi
    done
done
