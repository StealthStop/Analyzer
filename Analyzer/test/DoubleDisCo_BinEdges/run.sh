#!/bin/bash

# ----------------------------------------------------------------------
# How run this script:
# ./run.sh All      for all commands
# ./run.sh BinEdges for this specific command
#
#
# Optimized_BinEdges Class makes only table based on optimization metric
#   -- then you can get optimized ABCD bin edges based on that table
#   -- for now we use non-optimized ABCD bin edges (0.6, 0.6) 
# 
# BinEdges Class makes all standard plots
#   -- 1D closure per-njets
#   -- 2D non-closure, pull, signal fractions, significance
#
# MCcorrectionFactor_TT Calss makes plots for only TT
#   -- 1D plots as a function of baoundary value
#
# MCcorrectionFactor_TTvar Class makes
#   -- MC-based systematics plots
#   -- Data-based systematics plots
#   -- Input root files, including syst., for Higgs Combine
# 
# "plotting_MCcorrRatio.py" makes 
#   -- Closure Correction Ratio plots for MC-based syst.
# ----------------------------------------------------------------------


#DATE=("11.10.2022")
command=$1
command2=$2 # channel option on the command line

YEARS=("Run2UL") #"2016preVFP" "2016postVFP" "2017" "2018")
MODELS=("RPV" "SYY") 
#MODELS=("SYY")
MASSES=("550")
CHANNELS=("0l" "1l" "2l")
DISC1s=("0.6") 
DISC2s=("0.6") 
#CHANNELS=("1l")

for YEAR in ${YEARS[@]}; do

    for MODEL in ${MODELS[@]}; do

        for MASS in ${MASSES[@]}; do

            for CHANNEL in ${CHANNELS[@]}; do

                for DISC1 in ${DISC1s[@]}; do 

                    for DISC2 in ${DISC2s[@]}; do
                   
                        echo "year: ${YEAR}"

                        #################################
                        if [ $command == "Optimized_BinEdges" ] || [ $command == "All" ]; then
                             echo "running for Optimized_BinEdges:---------------------------------"

                             python run_DoubleDisCo_Validation.py --run Optimized_BinEdges --year ${YEAR} --channel 0l --sig ${MODEL} --mass ${MASS} --njets "8" "9" "10" "11" "12" "13incl" 
                             python run_DoubleDisCo_Validation.py --run Optimized_BinEdges --year ${YEAR} --channel 1l --sig ${MODEL} --mass ${MASS} --njets "7" "8" "9" "10" "11" "12incl"
                             python run_DoubleDisCo_Validation.py --run Optimized_BinEdges --year ${YEAR} --channel 2l --sig ${MODEL} --mass ${MASS} --njets "6" "7" "8" "9" "10" "11incl"
 
                        fi
                        
                        #################################
                        if [ $command == "BinEdges" ] || [ $command == "All" ]; then 
                            echo "running for BinEdges:---------------------------------"

                            python run_DoubleDisCo_Validation.py --run BinEdges --year ${YEAR} --channel 0l --sig ${MODEL} --mass ${MASS} --disc1edge ${DISC1} --disc2edge ${DISC2} --plotVars2D --plotDisc1VsDisc2 --njets "8" "9" "10" "11" "12" "13incl"
                            python run_DoubleDisCo_Validation.py --run BinEdges --year ${YEAR} --channel 1l --sig ${MODEL} --mass ${MASS} --disc1edge ${DISC1} --disc2edge ${DISC2} --plotVars2D --plotDisc1VsDisc2 --njets "7" "8" "9" "10" "11" "12incl"
                            python run_DoubleDisCo_Validation.py --run BinEdges --year ${YEAR} --channel 2l --sig ${MODEL} --mass ${MASS} --disc1edge ${DISC1} --disc2edge ${DISC2} --plotVars2D --plotDisc1VsDisc2 --njets "6" "7" "8" "9" "10" "11incl" 

                        fi
                                                    
                        #################################
                        if [ $command == "MCcorrectionFactor_TT" ] || [ $command == "All" ]; then 
                            echo "running for MCcorrectionFactor_TT:---------------------------------"

                            python run_DoubleDisCo_Validation.py --run MCcorrectionFactor_TT --year ${YEAR} --channel 0l --sig ${MODEL} --mass ${MASS} --disc1edge ${DISC1} --disc2edge ${DISC2} --plotVarVsBoundary --fastMode --njets "8" "9" "10" "11" "12" "13incl"
                            python run_DoubleDisCo_Validation.py --run MCcorrectionFactor_TT --year ${YEAR} --channel 1l --sig ${MODEL} --mass ${MASS} --disc1edge ${DISC1} --disc2edge ${DISC2} --plotVarVsBoundary --fastMode --njets "7" "8" "9" "10" "11" "12incl"
                            python run_DoubleDisCo_Validation.py --run MCcorrectionFactor_TT --year ${YEAR} --channel 2l --sig ${MODEL} --mass ${MASS} --disc1edge ${DISC1} --disc2edge ${DISC2} --plotVarVsBoundary --fastMode --njets "6" "7" "8" "9" "10" "11incl"

                        fi

                        #################################
                        if [ $command == "MCcorrectionFactor_TTvar" ] || [ $command == "All" ]; then 
                            echo "running for MCcorrectionFactor_TTvar:---------------------------------"

                            python run_DoubleDisCo_Validation.py --run MCcorrectionFactor_TTvar --year ${YEAR} --channel 0l --sig ${MODEL} --mass ${MASS} --disc1edge ${DISC1} --disc2edge ${DISC2} --plotVarVsBoundary --fastMode --njets "8" "9" "10" "11" "12" "13incl"
                            python run_DoubleDisCo_Validation.py --run MCcorrectionFactor_TTvar --year ${YEAR} --channel 1l --sig ${MODEL} --mass ${MASS} --disc1edge ${DISC1} --disc2edge ${DISC2} --plotVarVsBoundary --fastMode --njets "7" "8" "9" "10" "11" "12incl"
                            python run_DoubleDisCo_Validation.py --run MCcorrectionFactor_TTvar --year ${YEAR} --channel 2l --sig ${MODEL} --mass ${MASS} --disc1edge ${DISC1} --disc2edge ${DISC2} --plotVarVsBoundary --fastMode --njets "6" "7" "8" "9" "10" "11incl"

                        fi

                    done 
                done
            done                
        done
    done
done


# ------------------------------------------
# MC-based systematics
#   -- Closure Correction Ratio = ttvar / tt
# -----------------------------------------
for MODEL in ${MODELS[@]}; do

        for MASS in ${MASSES[@]}; do

            if [ $command == "mc" ] || [ $command == "All" ]; then
            
                python plotting_MCcorrRatio.py --sig ${MODEL} --mass ${MASS}
            
            fi
    done
done
