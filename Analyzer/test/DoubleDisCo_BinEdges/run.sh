#!/bin/bash

# ----------------------------------------------------------------------
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
YEARS=("Run2UL" "2016" "2017" "2018")
#MODELS=("RPV")
CHANNELS=("1l")
DISC1s=("0.6")
DISC2s=("0.6")

    for YEAR in ${YEARS[@]}; do

    #for MODEL in ${MODELS[@]}; do

        for CHANNEL in ${CHANNELS[@]}; do

            for DISC1 in ${DISC1s[@]}; do 

                for DISC2 in ${DISC2s[@]}; do
                
                    echo "year: ${YEAR}"

                    echo "running for Optimized_BinEdges:---------------------------------"
                    python run_DoubleDisCo_Validation.py --run Optimized_BinEdges --year ${YEAR} --channel ${CHANNEL} --numEdgeChoices 6

                    echo "running for BinEdges:---------------------------------"
                    python run_DoubleDisCo_Validation.py --run BinEdges --year ${YEAR} --channel ${CHANNEL} --disc1edge ${DISC1} --disc2edge ${DISC2} --plotVars2D 

                    echo "running for MCcorrectionFactor_TT:---------------------------------"
                    python run_DoubleDisCo_Validation.py --run MCcorrectionFactor_TT --year ${YEAR} --channel ${CHANNEL} --disc1edge ${DISC1} --disc2edge ${DISC2} --plotVarVsBoundary --fastMode

                    echo "running for MCcorrectionFactor_TTvar:---------------------------------"
                    python run_DoubleDisCo_Validation.py --run MCcorrectionFactor_TTvar --year ${YEAR} --channel ${CHANNEL} --disc1edge ${DISC1} --disc2edge ${DISC2} --plotVarVsBoundary --fastMode                    
                done
            done                
        done
    #done
    done

# Closure Correction Ratio 
python plotting_MCcorrRatio.py
