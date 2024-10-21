##########################################################################
# Used optimized ABCD bin edges by Max Significance (Higgs combine method)
########################################################################## 

command=$1

YEARS=("Run2UL") #"2016preVFP" "2016postVFP" "2017" "2018")
MODELS=("SYY")
#MASSES=("300" "350" "400" "450" "500" "550" "600" "650" "700" "750" "800" "850" "900" "950" "1000" "1050" "1100" "1150" "1200" "1250" "1300" "1350" "1400")
MASSES=("550")
CHANNELS=("0l" "1l" "2l")
OUTPATH="../condor/DisCo_outputs_0l_1l_2l_MaxSign_2_29_24"


for YEAR in ${YEARS[@]}; do

    for MODEL in ${MODELS[@]}; do

        for MASS in ${MASSES[@]}; do

                echo "year: ${YEAR}"


                #################################
                if [ $command == "BinEdges" ] || [ $command == "All" ]; then
                    echo "running for BinEdges:---------------------------------"

                    python run_DoubleDisCo_Validation.py --run BinEdges --year ${YEAR} --outpath Run2UL_MaxSign_2_29_24_StealthSYY --channel 0l --sig ${MODEL} --mass ${MASS} --disc1edge 0.76 --disc2edge 0.70 --plotVars2D --plotDisc1VsDisc2 --njets "8" "9" "10" "11" "12incl"  --path ${OUTPATH}
                    #python run_DoubleDisCo_Validation.py --run BinEdges --year ${YEAR} --outpath Run2UL_MaxSign_2_29_24_StealthSYY --channel 1l --sig ${MODEL} --mass ${MASS} --disc1edge 0.44 --disc2edge 0.42 --plotVars2D --plotDisc1VsDisc2 --njets "7" "8" "9" "10" "11incl"  --path ${OUTPATH}
                    python run_DoubleDisCo_Validation.py --run BinEdges --year ${YEAR} --outpath Run2UL_MaxSign_2_29_24_StealthSYY --channel 1l --sig ${MODEL} --mass ${MASS} --disc1edge 0.68 --disc2edge 0.48 --plotVars2D --plotDisc1VsDisc2 --njets "7" "8" "9" "10" "11incl"  --path ${OUTPATH}
                    python run_DoubleDisCo_Validation.py --run BinEdges --year ${YEAR} --outpath Run2UL_MaxSign_2_29_24_StealthSYY --channel 2l --sig ${MODEL} --mass ${MASS} --disc1edge 0.40 --disc2edge 0.42 --plotVars2D --plotDisc1VsDisc2 --njets "6" "7" "8" "9" "10incl"  --path ${OUTPATH}

                fi

                #################################
                if [ $command == "MCcorrectionFactor_TT" ] || [ $command == "All" ]; then
                    echo "running for MCcorrectionFactor_TT:---------------------------------"

                    #python run_DoubleDisCo_Validation.py --run MCcorrectionFactor_TT --year ${YEAR} --outpath Run2UL_MaxSign_2_29_24_StealthSYY --channel 0l --sig ${MODEL} --mass ${MASS} --disc1edge 0.76 --disc2edge 0.70 --plotVarVsBoundary --njets "8" "9" "10" "11" "12incl" --path ${OUTPATH}  --fastMode
                    python run_DoubleDisCo_Validation.py --run MCcorrectionFactor_TT --year ${YEAR} --outpath Run2UL_MaxSign_2_29_24_StealthSYY --channel 1l --sig ${MODEL} --mass ${MASS} --disc1edge 0.44 --disc2edge 0.42 --plotVarVsBoundary --njets "7" "8" "9" "10" "11incl" --path ${OUTPATH}  --fastMode
                    #python run_DoubleDisCo_Validation.py --run MCcorrectionFactor_TT --year ${YEAR} --outpath Run2UL_MaxSign_2_29_24_StealthSYY --channel 1l --sig ${MODEL} --mass ${MASS} --disc1edge 0.68 --disc2edge 0.48 --plotVarVsBoundary --njets "7" "8" "9" "10" "11incl" --path ${OUTPATH} --fastMode
                    #python run_DoubleDisCo_Validation.py --run MCcorrectionFactor_TT --year ${YEAR} --outpath Run2UL_MaxSign_2_29_24_StealthSYY --channel 2l --sig ${MODEL} --mass ${MASS} --disc1edge 0.40 --disc2edge 0.42 --plotVarVsBoundary --njets "6" "7" "8" "9" "10incl" --path ${OUTPATH}  --fastMode

                fi

                #################################
                if [ $command == "MCcorrectionFactor_TTvar" ] || [ $command == "All" ]; then
                    echo "running for MCcorrectionFactor_TTvar:---------------------------------"

                    python run_DoubleDisCo_Validation.py --run MCcorrectionFactor_TTvar --year ${YEAR} --outpath Run2UL_MaxSign_2_29_24_StealthSYY --channel 0l --sig ${MODEL} --mass ${MASS} --disc1edge 0.76 --disc2edge 0.70 --plotVarVsBoundary --fastMode --njets "8" "9" "10" "11" "12incl" --path ${OUTPATH}
                    #python run_DoubleDisCo_Validation.py --run MCcorrectionFactor_TTvar --year ${YEAR} --outpath Run2UL_MaxSign_2_29_24_StealthSYY --channel 1l --sig ${MODEL} --mass ${MASS} --disc1edge 0.44 --disc2edge 0.42 --plotVarVsBoundary --fastMode --njets "7" "8" "9" "10" "11incl" --path ${OUTPATH}
                    python run_DoubleDisCo_Validation.py --run MCcorrectionFactor_TTvar --year ${YEAR} --outpath Run2UL_MaxSign_2_29_24_StealthSYY --channel 1l --sig ${MODEL} --mass ${MASS} --disc1edge 0.68 --disc2edge 0.48 --plotVarVsBoundary --fastMode --njets "7" "8" "9" "10" "11incl" --path ${OUTPATH}
                    python run_DoubleDisCo_Validation.py --run MCcorrectionFactor_TTvar --year ${YEAR} --outpath Run2UL_MaxSign_2_29_24_StealthSYY --channel 2l --sig ${MODEL} --mass ${MASS} --disc1edge 0.40 --disc2edge 0.42 --plotVarVsBoundary --fastMode --njets "6" "7" "8" "9" "10incl" --path ${OUTPATH}

                fi

            done
        done
    done

