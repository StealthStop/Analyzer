#ifndef MakeNJetDists_h
#define MakeNJetDists_h

#include "Analyzer/Analyzer/include/AnalyzeBase.h"
#include "Analyzer/Analyzer/include/Histo.h"

#include "NTupleReader/include/NTupleReader.h"

#include "Framework/Framework/include/RunTopTagger.h"
#include "Framework/Framework/include/Muon.h"
#include "Framework/Framework/include/Electron.h"
#include "Framework/Framework/include/Photon.h"
#include "Framework/Framework/include/Jet.h"
#include "Framework/Framework/include/BJet.h"
#include "Framework/Framework/include/CommonVariables.h"
#include "Framework/Framework/include/FatJetCombine.h"
#include "Framework/Framework/include/Baseline.h"
#include "Framework/Framework/include/MakeMVAVariables.h"
#include "Framework/Framework/include/DeepEventShape.h"
#include "Framework/Framework/include/BTagCorrector.h"
#include "Framework/Framework/include/ScaleFactors.h"

#include <TH1D.h>
#include <TH2D.h>
#include <iostream>

class MakeNJetDists : public AnalyzeBase
{
private:
    std::vector<std::pair<std::string, std::string>> myVarSuffixPairs;

public:
    void InitHistos(const std::string& runtype)
    {
        TH1::SetDefaultSumw2();
        TH2::SetDefaultSumw2();

        my_Histos.emplace_back(new Histo1D("EventCounter", 2,  -1.1, 1.1, "eventCounter", {}, {}));

        std::vector<std::string> weightVec;
        if( runtype == "MC" )
        {
            myVarSuffixPairs = {{"",""}, {"JECup","_JECUp"}, {"JECdown","_JECDown"}, {"JERup","_JERUp"}, {"JERdown","_JERDown"}, {"mpTScaled","_mpTScaled"}};
            weightVec = {"Lumi", "Weight"};

            //--------------------------------------------------------------------------------
            // Plots that are made only once and for MC only
            //--------------------------------------------------------------------------------

            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_btgUp",   6, 0.0,  6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Up", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor"}));
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_btgDown", 6, 0.0,  6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Down", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor"}));

            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_lepUp",   6, 0.0,  6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF_Up", "totGoodMuonSF_Up", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor"}));
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_lepDown", 6, 0.0,  6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF_Down", "totGoodMuonSF_Down", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor"}));

            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_isrUp",   6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_ISRUp"}));
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_isrDown", 6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_ISRDown"}));

            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_fsrUp",   6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_FSRUp"}));
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_fsrDown", 6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_FSRDown"}));

            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_isr2Up",   6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_ISRUp_2"}));
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_isr2Down", 6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_ISRDown_2"}));

            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_fsr2Up",   6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_FSRUp_2"}));
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_fsr2Down", 6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_FSRDown_2"}));

            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_pdfUp",   6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PDFweightUp"}));
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_pdfDown", 6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PDFweightDown"}));

            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_htUp",   6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htScaleUp", "puWeightCorr", "prefiringScaleFactor"}));
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_htDown", 6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htScaleDown", "puWeightCorr", "prefiringScaleFactor"}));

            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_puUp",   6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puSysUpCorr", "prefiringScaleFactor"}));
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_puDown", 6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puSysDownCorr", "prefiringScaleFactor"}));

            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_sclUp",   6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "scaleWeightUp"}));
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_sclDown", 6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "scaleWeightDown"}));
       
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_prfUp",   6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactorUp"}));
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_prfDown", 6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactorDown"}));

            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_noHT", 6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good"}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "puWeightCorr", "prefiringScaleFactor"}));
       
            for(int i = 0; i < 4; i++)
            {
                std::string index = std::to_string(i+1);

                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_btgUp",   6, 0.0,  6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Up", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor"}));
                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_btgDown", 6, 0.0,  6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Down", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor"}));

                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_lepUp",   6, 0.0,  6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF_Up", "totGoodMuonSF_Up", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor"}));
                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_lepDown", 6, 0.0,  6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF_Down", "totGoodMuonSF_Down", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor"}));

                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_isrUp",   6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_ISRUp"}));
                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_isrDown", 6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_ISRDown"}));
            
                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_fsrUp",   6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_FSRUp"}));
                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_fsrDown", 6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_FSRDown"}));
            
                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_isr2Up",   6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_ISRUp_2"}));
                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_isr2Down", 6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_ISRDown_2"}));
        
                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_fsr2Up",   6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_FSRUp_2"}));
                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_fsr2Down", 6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PSweight_FSRDown_2"}));

                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_pdfUp",   6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PDFweightUp"}));
                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_pdfDown", 6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "PDFweightDown"}));

                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_htUp",   6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htScaleUp","puWeightCorr", "prefiringScaleFactor"}));
                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_htDown", 6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htScaleDown", "puWeightCorr", "prefiringScaleFactor"}));
                
                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_puUp",   6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight","puSysUpCorr", "prefiringScaleFactor"}));
                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_puDown", 6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puSysDownCorr", "prefiringScaleFactor"}));

                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_sclUp",   6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "scaleWeightUp"}));
                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_sclDown", 6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactor", "scaleWeightDown"}));

                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_prfUp",   6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactorUp"}));
                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_prfDown", 6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "htDerivedweight", "puWeightCorr", "prefiringScaleFactorDown"}));

                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+"_noHT", 6,  0.0, 6.0, "NGoodJets_pt30_inclusive_shift", {"passBaseline1l_Good", "deepESM_bin"+index}, {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central", "totGoodElectronSF", "totGoodMuonSF", "puWeightCorr", "prefiringScaleFactor"}));
            }            
        }
        else
        {
            myVarSuffixPairs = {{"",""}};
            weightVec = {};
        }    

        //--------------------------------------------------------------------------------
        // Plots that are made with JEC/R variation and Data
        //--------------------------------------------------------------------------------
        for(const auto& pair : myVarSuffixPairs)
        {
            const std::string& s = pair.first;
            const std::string& n = pair.second;

            std::vector<std::string> weightVecNoHT, weightVecAll;
            if( runtype == "MC" )
            {
                weightVecNoHT = {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central"+s, "totGoodElectronSF"+s, "totGoodMuonSF"+s,                      "puWeightCorr"+s, "prefiringScaleFactor"+s};
                weightVecAll  = {"Lumi", "Weight", "bTagSF_EventWeightSimple_Central"+s, "totGoodElectronSF"+s, "totGoodMuonSF"+s, "htDerivedweight"+s, "puWeightCorr"+s, "prefiringScaleFactor"+s};
            }
            else
            {
                weightVecNoHT = {};
                weightVecAll = {};                
            }

            //-----------------------------------------------------------------
            // NJet plots
            //-----------------------------------------------------------------       
            my_Histos.emplace_back(new Histo1D("h_njetsShifted"+n,  6, 0.0,  6.0, "NGoodJets_pt30_inclusive_shift"+s, {}, weightVec));
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l"+n+"noHT", 6, 0.0,  6.0, "NGoodJets_pt30_inclusive_shift"+s, {"passBaseline1l_Good"+s}, weightVecNoHT));
            my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l"+n,        6, 0.0,  6.0, "NGoodJets_pt30_inclusive_shift"+s, {"passBaseline1l_Good"+s}, weightVecAll));
            for(int i = 0; i < 4; i++)
            {
                std::string index = std::to_string(i+1);
                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+n+"noHT", 6, 0.0, 6.0, "NGoodJets_pt30_inclusive_shift"+s, {"passBaseline1l_Good"+s,"deepESM_bin"+index+s}, weightVecNoHT));
                my_Histos.emplace_back(new Histo1D("h_njetsShifted_pt30_1l_D"+index+n,        6, 0.0, 6.0, "NGoodJets_pt30_inclusive_shift"+s, {"passBaseline1l_Good"+s,"deepESM_bin"+index+s}, weightVecAll));
            }
        }
    }//END of init histos

    void Loop(NTupleReader& tr, double, int maxevents, bool)
    {
        const auto& runtype = tr.getVar<std::string>("runtype");
        const auto& filetag = tr.getVar<std::string>("filetag");
        const auto& runYear = tr.getVar<std::string>("runYear");
        const auto& DeepESMCfg = tr.getVar<std::string>("DeepESMCfg");
        const auto& ModelFile = tr.getVar<std::string>("DeepESMModel");
        const auto& bjetFileName = tr.getVar<std::string>("bjetFileName");
        const auto& bjetCSVFileName = tr.getVar<std::string>("bjetCSVFileName");
        const auto& leptonFileName = tr.getVar<std::string>("leptonFileName");
        const auto& meanFileName = tr.getVar<std::string>("meanFileName");
        const auto& TopTaggerCfg    = tr.getVar<std::string>("TopTaggerCfg");

        //-------------------------------------
        //-- Initialize histograms to be filled
        //-------------------------------------
        InitHistos(runtype);

        for(const auto& pair : myVarSuffixPairs)
        {
            const std::string& myVarSuffix = pair.first;
            if(myVarSuffix == "") continue;
            Jet                 jet(myVarSuffix);
            BJet                bjet(myVarSuffix);
            Muon                muon(myVarSuffix);
            Photon              photon(myVarSuffix);
            Baseline            baseline(myVarSuffix);
            Electron            electron(myVarSuffix);
            ScaleFactors        scaleFactors( runYear, leptonFileName, meanFileName, myVarSuffix);
            RunTopTagger        topTagger(TopTaggerCfg, myVarSuffix);
            BTagCorrector       bTagCorrector(bjetFileName, "", bjetCSVFileName, filetag);
            DeepEventShape      deepEventShape(DeepESMCfg, ModelFile, "Info", true, myVarSuffix);
            CommonVariables     commonVariables(myVarSuffix);
            MakeMVAVariables    makeMVAVariables(false, myVarSuffix);
            FatJetCombine       fatJetCombine(myVarSuffix);
            bTagCorrector.SetVarNames("GenParticles_PdgId", "Jets"+myVarSuffix, "GoodJets_pt30"+myVarSuffix, "Jets"+myVarSuffix+"_bJetTagDeepCSVtotb", "Jets"+myVarSuffix+"_partonFlavor", myVarSuffix);
  
            // Remember, order matters here !
            // Follow what is done in Config.h
            tr.registerFunction(muon);
            tr.registerFunction(electron);
            tr.registerFunction(photon);
            tr.registerFunction(jet);
            tr.registerFunction(bjet);
            tr.registerFunction(topTagger);
            tr.registerFunction(commonVariables);
            tr.registerFunction(fatJetCombine);
            tr.registerFunction(baseline);
            tr.registerFunction(makeMVAVariables);
            tr.registerFunction(deepEventShape);
            tr.registerFunction(bTagCorrector);
            tr.registerFunction(scaleFactors);

        }

        while( tr.getNextEvent() )
        {
            //------------------------------------
            //-- Print Event Number
            //------------------------------------
            const bool breakLoop = printEventNum(maxevents, tr.getEvtNum());
            if(breakLoop) break;

            //-----------------------------------
            //-- Fill Histograms Below
            //-----------------------------------
            Fill(tr);
        }//END of while tr.getNextEvent loop   
    }//END of function
};

#endif
