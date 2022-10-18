#define AnalyzeDoubleDisCo_cxx
#include "Analyzer/Analyzer/include/AnalyzeDoubleDisCo.h"
#include "NTupleReader/include/NTupleReader.h"
#include "Framework/Framework/include/Utility.h"

#include "Framework/Framework/include/Muon.h"
#include "Framework/Framework/include/Electron.h"
#include "Framework/Framework/include/Photon.h"
#include "Framework/Framework/include/Jet.h"
#include "Framework/Framework/include/BJet.h"
#include "Framework/Framework/include/CommonVariables.h"
#include "Framework/Framework/include/RunTopTagger.h"
#include "Framework/Framework/include/FatJetCombine.h"
#include "Framework/Framework/include/MakeMVAVariables.h"
#include "Framework/Framework/include/Baseline.h"
#include "Framework/Framework/include/BTagCorrector.h"
#include "Framework/Framework/include/ScaleFactors.h"
#include "Framework/Framework/include/StopJets.h"
#include "Framework/Framework/include/StopGenMatch.h"
#include "Framework/Framework/include/MakeStopHemispheres.h"
#include "Framework/Framework/include/DeepEventShape.h"

#include <TH1D.h>
#include <TH2D.h>
#include <iostream>

AnalyzeDoubleDisCo::AnalyzeDoubleDisCo() : initHistos(false)
{

    histInfos = {
        {"h_DoubleDisCo_disc1",          80,    0,    1},
        {"h_DoubleDisCo_disc2",          80,    0,    1},
        {"h_DoubleDisCo_massReg",       150,    0, 1500},
        {"fwm2_top6",                    50,    0,    1},
        {"fwm3_top6",                    50,    0,    1},
        {"fwm4_top6",                    50,    0,    1},
        {"fwm5_top6",                    50,    0,    1},
        {"jmt_ev0_top6",                 50,    0,    1},
        {"jmt_ev1_top6",                 50,    0,    1},
        {"jmt_ev2_top6",                 50,    0,    1},
        {"Stop1_pt_cm_OldSeed",          72,    0, 1500},
        {"Stop1_eta_cm_OldSeed",         80,   -6,    6},
        {"Stop1_phi_cm_OldSeed",         64,   -4,    4},
        {"Stop1_mass_cm_OldSeed",        72,    0, 1500},
        {"Stop2_pt_cm_OldSeed",          72,    0, 1500},
        {"Stop2_eta_cm_OldSeed",         80,   -6,    6},
        {"Stop2_phi_cm_OldSeed",         64,   -4,    4},
        {"Stop2_mass_cm_OldSeed",        72,    0, 1500},
        {"Stop1_pt_cm_TopSeed",          72,    0, 1500},
        {"Stop1_eta_cm_TopSeed",         80,   -6,    6},
        {"Stop1_phi_cm_TopSeed",         64,   -4,    4},
        {"Stop1_mass_cm_TopSeed",        72,    0, 1500},
        {"Stop2_pt_cm_TopSeed",          72,    0, 1500},
        {"Stop2_eta_cm_TopSeed",         80,   -6,    6},
        {"Stop2_phi_cm_TopSeed",         64,   -4,    4},
        {"Stop2_mass_cm_TopSeed",        72,    0, 1500},
        {"h_njets",                      20,    0, 20.0},
        {"h_njets_11incl",               20, -0.5, 19.5},
        {"h_njets_12incl",               24, -0.5, 23.5},
        {"h_njets_13incl",               28, -0.5, 27.5},
        {"h_njets_14incl",               32, -0.5, 31.5},
        {"h_ht",                        500,    0, 5000},
        {"Stop1_mass_PtRank_matched",   360,    0, 1500},
        {"Stop2_mass_PtRank_matched",   360,    0, 1500},
        {"Stop1_mass_MassRank_matched", 360,    0, 1500},
        {"Stop2_mass_MassRank_matched", 360,    0, 1500},
        {"Stop_mass_average_matched",   360,    0, 1500},
        {"h_jPt",                       200,    0, 2000},
        {"h_jPhi",                      200,   -4,    4},
        {"h_jEta",                      200,   -6,    6},
        {"h_jM",                        200,    0,  200},
        {"h_lPt",                       200,    0, 2000},
        {"h_lIso",                      360,    0, 1500},
        {"h_lPhi",                      200,   -4,    4},
        {"h_lEta",                      200,   -6,    6},
        {"h_lCharge",                     2,   -1,    1},
        {"h_lMiniIso",                 4000,    0,    5},
        {"h_ePt",                       200,    0, 2000},
        {"h_eIso",                      360,    0, 1500},
        {"h_ePhi",                      200,   -4,    4},
        {"h_eEta",                      200,   -6,    6},
        {"h_eCharge",                     2,   -1,    1},
        {"h_eMiniIso",                  400,    0,    5},
        {"h_mPt",                       360,    0, 1500},
        {"h_mIso",                      360,    0, 1500},
        {"h_mPhi",                      200,   -4,    4},
        {"h_mEta",                      200,   -6,    6},
        {"h_mCharge",                     2,   -1,    1},
        {"h_mMiniIso",                 4000,    0,    5},
    };

    hist2DInfos = {
        {"h_DoubleDisCo_disc1_disc2", 100, 0, 1, 100, 0, 1}, 
        {"h_cm_pt_jetRank", 150, 0, 1500, 10, 0, 10},
        {"h_cm_ptrHT_jetRank", 150, 0, 1, 10, 0, 10},
    };

    njets = {"Incl", "7", "8", "9", "10", "11", "12", "12incl"};

    my_var_suffix = {""};
    //my_var_suffix = {"", "JECup", "JECdown", "JERup", "JERdown"};


}

void AnalyzeDoubleDisCo::makeSubregions(const std::vector<std::vector<std::string>>& regionVec)
{

    // Translate the generic "A", "B", "C", "D" names into unique names based on the region.
    // E.g. region "fixedbdEF" would yield "b", "d", "E", "F"
    // Naming conventions agreed upon with Validation code
    for (auto& regions : regionVec)
    {
        for (auto& region : regions)
        { 

            if (subRegionsMap.find(region) != subRegionsMap.end())
                continue;

            if (region.find("subDiv") == std::string::npos) {

                int offset;
                if (region.find("fixed") != std::string::npos) {
                    offset = 5;
                }
                else {
                    offset = 0;
                }

                // Attention to subregion names for bdEF
                if (region[0+offset] == 'b' or region[0+offset] == 'B') {
                    subRegionsMap[region].push_back(std::string(1, region[0+offset]));
                    subRegionsMap[region].push_back(std::string(1, region[2+offset]));
                    subRegionsMap[region].push_back(std::string(1, region[1+offset]));
                    subRegionsMap[region].push_back(std::string(1, region[3+offset]));
                } else {
                    subRegionsMap[region].push_back(std::string(1, region[0+offset]));
                    subRegionsMap[region].push_back(std::string(1, region[1+offset]));
                    subRegionsMap[region].push_back(std::string(1, region[2+offset]));
                    subRegionsMap[region].push_back(std::string(1, region[3+offset]));
                }

            } else if (region.find("subDiv") != std::string::npos) {
                subRegionsMap[region].push_back("dA");
                subRegionsMap[region].push_back("dB");
                subRegionsMap[region].push_back("dC");
                subRegionsMap[region].push_back("dD");
            } else {
                subRegionsMap[region].push_back("A");
                subRegionsMap[region].push_back("B");
                subRegionsMap[region].push_back("C");
                subRegionsMap[region].push_back("D");
            }
        }
    }
    subRegionsMap["Incl"].push_back("A");
    subRegionsMap["Incl"].push_back("B");
    subRegionsMap["Incl"].push_back("C");
    subRegionsMap["Incl"].push_back("D");
}

void AnalyzeDoubleDisCo::Preinit(unsigned int nNNJets)
{
    for(unsigned int i = 1; i <= nNNJets ; i++)
    {

        histInfos.push_back({"Jet_cm_ptrHT_"   + std::to_string(i), 150,  0,    1});
        histInfos.push_back({"Jet_cm_pt_"      + std::to_string(i), 150,  0, 1500});
        histInfos.push_back({"Jet_cm_eta_"     + std::to_string(i), 100, -6,    6});
        histInfos.push_back({"Jet_cm_phi_"     + std::to_string(i),  80, -4,    4});
        histInfos.push_back({"Jet_cm_m_"       + std::to_string(i), 150,  0,  300});
        histInfos.push_back({"Jet_cm_E_"       + std::to_string(i), 150,  0, 1500});
        histInfos.push_back({"Jet_cm_flavb_"   + std::to_string(i),  80,  0,    1});
        histInfos.push_back({"Jet_cm_flavc_"   + std::to_string(i),  80,  0,    1});
        histInfos.push_back({"Jet_cm_flavg_"   + std::to_string(i),  80,  0,    1});
        histInfos.push_back({"Jet_cm_flavq_"   + std::to_string(i),  80,  0,    1});
        histInfos.push_back({"Jet_cm_flavuds_" + std::to_string(i),  80,  0,    1});
    }
}

void AnalyzeDoubleDisCo::InitHistos(const std::map<std::string, bool>& cutMap, const std::vector<std::vector<std::string>>& regionsVec)
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    // Generates a map of region to constituent subregions
    makeSubregions(regionsVec);

    my_histos.emplace("EventCounter", std::make_shared<TH1D>("EventCounter", "EventCounter", 2, -1.1, 1.1) );

    // mycut.first : cut string, mycut.second : boolean if cut passed
    for (auto& mycut : cutMap)
    {
        for (const auto& hInfo : histInfos)
        { 
            // Njet string, can also be "Incl"
            for (const auto& Njet : njets)
            {
                std::string njetStr = "";
                // For 1D histos, don't make any where we exclude all but one njets bin
                if (Njet != "Incl")
                    njetStr = "_Njets" + Njet;

                // regions : a vector of region string names for 0L or 1L
                for (const auto& regionPair : subRegionsMap)
                {
                    std::string region = regionPair.first; 
                    // Only make mass reg histos for per-njets scenarios not inclusive
                    if (hInfo.name.find("massReg") != std::string::npos and (Njet == "Incl" or region == "Incl"))
                        continue;

                    // Only make njets histos for njets-inclusive scenarios
                    if (hInfo.name.find("njets") != std::string::npos and Njet != "Incl")
                        continue;

                    // For Combine njets histos, skip those for the inclusive region
                    if (hInfo.name.find("njets") != std::string::npos and hInfo.name.find("incl") != std::string::npos and region == "Incl")
                        continue;

                    if (hInfo.name.find("DoubleDisCo") == std::string::npos and Njet != "Incl")
                        continue;

                    if (hInfo.name.find("disc") != std::string::npos and Njet == "Incl")
                        continue;

                    std::string regionStr = "";
                    if (region != "Incl")
                        regionStr = "_" + region;

                    std::string name = hInfo.name + mycut.first + njetStr + regionStr;
                    my_histos.emplace(name, std::make_shared<TH1D>((name).c_str(),(name).c_str(), hInfo.nBins, hInfo.low, hInfo.high));
                }
            }
        }

        for(const auto& h2dInfo : hist2DInfos)
        {
            for (const auto& Njet : njets)
            {
                // For the disc1 vs disc2 plots, no need for njets inclusive one
                if (Njet == "Incl" && h2dInfo.name.find("disco") != std::string::npos)
                    continue;

                std::string njetStr = "";
                // For 2D histos, don't make any where we exclude all but one njets bin
                if (Njet != "Incl")
                    njetStr = "_Njets" + Njet;

                for (const auto& regionPair : subRegionsMap)
                {
                    std::string region = regionPair.first;
                    std::vector<std::string> subregions = regionPair.second;

                    std::string regionStr = "";
                    if (region != "Incl")
                        regionStr = "_" + region;

                    for (const auto& subregion : subregions)
                    {
                        std::string name = h2dInfo.name + mycut.first + njetStr + regionStr + "_" + subregion;
                        my_2d_histos.emplace(name, std::make_shared<TH2D>((name).c_str(),(name).c_str(), h2dInfo.nBinsX, h2dInfo.lowX, h2dInfo.highX, h2dInfo.nBinsY, h2dInfo.lowY, h2dInfo.highY));
                    }
                    std::string name = h2dInfo.name + mycut.first + njetStr + regionStr;
                    my_2d_histos.emplace(name, std::make_shared<TH2D>((name).c_str(),(name).c_str(), h2dInfo.nBinsX, h2dInfo.lowX, h2dInfo.highX, h2dInfo.nBinsY, h2dInfo.lowY, h2dInfo.highY));
                }
            }
        }
    }
}

void AnalyzeDoubleDisCo::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    const auto& dataset                           = tr.getVar<std::string>("dataset"                          ); // to check the collection name for FSR/ISR, JEC/JER 
    const auto& filetag                           = tr.getVar<std::string>("filetag"                          );
    const auto& runYear                           = tr.getVar<std::string>("runYear"                          );
    const auto& bjetFileName                      = tr.getVar<std::string>("bjetFileName"                     );
    const auto& bjetCSVFileName                   = tr.getVar<std::string>("bjetCSVFileName"                  );
    const auto& leptonFileName                    = tr.getVar<std::string>("leptonFileName"                   );
    const auto& hadronicFileName                  = tr.getVar<std::string>("hadronicFileName"                 );
    const auto& meanFileName                      = tr.getVar<std::string>("meanFileName"                     );
    const auto& TopTaggerCfg                      = tr.getVar<std::string>("TopTaggerCfg"                     );
    const auto& DoubleDisCo_Cfg_0l_RPV            = tr.getVar<std::string>("DoubleDisCo_Cfg_0l_RPV"           );
    const auto& DoubleDisCo_Model_0l_RPV          = tr.getVar<std::string>("DoubleDisCo_Model_0l_RPV"         );
    const auto& DoubleDisCo_Cfg_NonIsoMuon_0l_RPV = tr.getVar<std::string>("DoubleDisCo_Cfg_NonIsoMuon_0l_RPV");
    const auto& DoubleDisCo_Cfg_1l_RPV            = tr.getVar<std::string>("DoubleDisCo_Cfg_1l_RPV"           );  
    const auto& DoubleDisCo_Model_1l_RPV          = tr.getVar<std::string>("DoubleDisCo_Model_1l_RPV"         );    
    const auto& DoubleDisCo_Cfg_NonIsoMuon_1l_RPV = tr.getVar<std::string>("DoubleDisCo_Cfg_NonIsoMuon_1l_RPV");

    // to check the collection name for JEC/JER 
    if (dataset.find("JECup") != std::string::npos)
        my_var_suffix = {"JECup"};

    else if (dataset.find("JECdown") != std::string::npos)
        my_var_suffix = {"JECdown"};
        
    else if (dataset.find("JERup") != std::string::npos)
        my_var_suffix = {"JERup"};

    else if (dataset.find("JERdown") != std::string::npos)
        my_var_suffix = {"JERdown"};


    for(const auto& myVarSuffix : my_var_suffix)
    {

        if (myVarSuffix == "") continue;
        Jet                 jet(myVarSuffix);
        BJet                bjet(myVarSuffix);
        Muon                muon(myVarSuffix);
        Photon              photon(myVarSuffix);
        Baseline            baseline(myVarSuffix);
        Electron            electron(myVarSuffix);
        StopJets            stopJets(myVarSuffix);
        RunTopTagger        topTagger(TopTaggerCfg, myVarSuffix);
        ScaleFactors        scaleFactors( runYear, leptonFileName, hadronicFileName, meanFileName, myVarSuffix);
        StopGenMatch        stopGenMatch(myVarSuffix);
        FatJetCombine       fatJetCombine(myVarSuffix);
        BTagCorrector       bTagCorrector(bjetFileName, "", bjetCSVFileName, filetag);
        DeepEventShape      neuralNetwork0L(DoubleDisCo_Cfg_0l_RPV, DoubleDisCo_Model_0l_RPV, "Info", true, myVarSuffix);
        DeepEventShape      neuralNetwork0L_NonIsoMuon(DoubleDisCo_Cfg_NonIsoMuon_0l_RPV, DoubleDisCo_Model_0l_RPV, "Info", true, myVarSuffix);
        DeepEventShape      neuralNetwork1L(DoubleDisCo_Cfg_1l_RPV, DoubleDisCo_Model_1l_RPV, "Info", true, myVarSuffix); 
        DeepEventShape      neuralNetwork1L_NonIsoMuon(DoubleDisCo_Cfg_NonIsoMuon_1l_RPV, DoubleDisCo_Model_1l_RPV, "Info", true, myVarSuffix);
        CommonVariables     commonVariables(myVarSuffix);
        MakeMVAVariables    makeMVAVariables0L_NonIsoMuon(false, myVarSuffix, "GoodJets_pt30",       false, true, 7,  2, "_0l");
        MakeMVAVariables    makeMVAVariables1L_NonIsoMuon(false, myVarSuffix, "NonIsoMuonJets_pt30", false, true, 7,  2, "_1l");
        MakeMVAVariables    makeMVAVariables0L(false, myVarSuffix, "GoodJets_pt30", false, true, 7, 2, "_0l");
        MakeMVAVariables    makeMVAVariables1L(false, myVarSuffix, "GoodJets_pt30", false, true, 7, 2, "_1l");
        MakeStopHemispheres stopHemispheres_OldSeed("Jets",     "GoodJets_pt20", "NGoodJets_pt20", "_OldSeed", myVarSuffix, Hemisphere::InvMassSeed);
        MakeStopHemispheres stopHemispheres_TopSeed("StopJets", "GoodStopJets",  "NGoodStopJets",  "_TopSeed", myVarSuffix, Hemisphere::TopSeed    );
        MakeStopHemispheres stopHemispheres_OldSeed_NonIsoMuon("Jets",     "NonIsoMuonJets_pt20", "NNonIsoMuonJets_pt30", "_OldSeed_NonIsoMuon", myVarSuffix, Hemisphere::InvMassSeed);
        MakeStopHemispheres stopHemispheres_TopSeed_NonIsoMuon("StopJets", "GoodStopJets",        "NGoodStopJets",        "_TopSeed_NonIsoMuon", myVarSuffix, Hemisphere::InvMassSeed);
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
        tr.registerFunction(baseline);
        tr.registerFunction(fatJetCombine);
        tr.registerFunction(makeMVAVariables0L_NonIsoMuon);
        tr.registerFunction(makeMVAVariables1L_NonIsoMuon);
        tr.registerFunction(makeMVAVariables0L);
        tr.registerFunction(makeMVAVariables1L);
        tr.registerFunction(stopJets);
        tr.registerFunction(stopHemispheres_OldSeed);
        tr.registerFunction(stopHemispheres_TopSeed);
        tr.registerFunction(stopHemispheres_OldSeed_NonIsoMuon);
        tr.registerFunction(stopHemispheres_TopSeed_NonIsoMuon);
        tr.registerFunction(bTagCorrector);
        tr.registerFunction(scaleFactors);
        tr.registerFunction(stopGenMatch);
        tr.registerFunction(neuralNetwork0L);
        tr.registerFunction(neuralNetwork0L_NonIsoMuon);
        tr.registerFunction(neuralNetwork1L);
        tr.registerFunction(neuralNetwork1L_NonIsoMuon);

    }

    while( tr.getNextEvent() )
    {
        if (maxevents != -1 && tr.getEvtNum() >= maxevents)
            break;        

        if (tr.getEvtNum() % 1000 == 0)
            printf("  Event %i\n", tr.getEvtNum() );

        for(const auto& myVarSuffix : my_var_suffix)
        {

            const auto& runtype                           = tr.getVar<std::string>("runtype");    
            const auto& Jets                              = tr.getVec<utility::LorentzVector>("Jets"                       +myVarSuffix);
            const auto& GoodJets_pt30                     = tr.getVec<bool>("GoodJets_pt30"                                +myVarSuffix);
            const auto& NGoodJets_pt30                    = tr.getVar<int>("NGoodJets_pt30"                                +myVarSuffix);
            const auto& NNonIsoMuonJets_pt30              = tr.getVar<int>("NNonIsoMuonJets_pt30"                          +myVarSuffix);
            const auto& HT_trigger_pt30                   = tr.getVar<double>("HT_trigger_pt30"                            +myVarSuffix);
            const auto& HT_NonIsoMuon_pt30                = tr.getVar<double>("HT_NonIsoMuon_pt30"                         +myVarSuffix);

            const auto& passBaseline0l_Good               = tr.getVar<bool>("passBaseline0l_good"                          +myVarSuffix);
            const auto& passBaseline1l_Good               = tr.getVar<bool>("passBaseline1l_Good"                          +myVarSuffix);
            const auto& passBaseline1l_HT500_Good         = tr.getVar<bool>("passBaseline1l_HT500_Good"                    +myVarSuffix);
            const auto& passBaseline1l_HT700_Good         = tr.getVar<bool>("passBaseline1l_HT700_Good"                    +myVarSuffix);
            const auto& passBaseline0l_NonIsoMuon         = tr.getVar<bool>("pass_qcdCR"                                   +myVarSuffix);
            const auto& passBaseline1l_NonIsoMuon         = tr.getVar<bool>("passBaseline1l_NonIsoMuon"                    +myVarSuffix);

            const auto& DoubleDisCo_massReg_0l            = tr.getVar<double>("DoubleDisCo_massReg_0l_RPV"                 +myVarSuffix);
            const auto& DoubleDisCo_massReg_1l            = tr.getVar<double>("DoubleDisCo_massReg_1l_RPV"                 +myVarSuffix);
            const auto& DoubleDisCo_massReg_NonIsoMuon_0l = tr.getVar<double>("DoubleDisCo_massReg_NonIsoMuon_0l_RPV"      +myVarSuffix);
            const auto& DoubleDisCo_massReg_NonIsoMuon_1l = tr.getVar<double>("DoubleDisCo_massReg_NonIsoMuon_1l_RPV"      +myVarSuffix);
            const auto& DoubleDisCo_disc1_0l              = tr.getVar<double>("DoubleDisCo_disc1_0l_RPV"                   +myVarSuffix);
            const auto& DoubleDisCo_disc2_0l              = tr.getVar<double>("DoubleDisCo_disc2_0l_RPV"                   +myVarSuffix);
            const auto& DoubleDisCo_disc1_1l              = tr.getVar<double>("DoubleDisCo_disc1_1l_RPV"                   +myVarSuffix);
            const auto& DoubleDisCo_disc2_1l              = tr.getVar<double>("DoubleDisCo_disc2_1l_RPV"                   +myVarSuffix);
            const auto& DoubleDisCo_disc1_NonIsoMuon_0l   = tr.getVar<double>("DoubleDisCo_disc1_NonIsoMuon_0l_RPV"        +myVarSuffix);
            const auto& DoubleDisCo_disc2_NonIsoMuon_0l   = tr.getVar<double>("DoubleDisCo_disc2_NonIsoMuon_0l_RPV"        +myVarSuffix);
            const auto& DoubleDisCo_disc1_NonIsoMuon_1l   = tr.getVar<double>("DoubleDisCo_disc1_NonIsoMuon_1l_RPV"        +myVarSuffix);
            const auto& DoubleDisCo_disc2_NonIsoMuon_1l   = tr.getVar<double>("DoubleDisCo_disc2_NonIsoMuon_1l_RPV"        +myVarSuffix);

            const auto& fwm2_top6_0l                      = tr.getVar<double>("fwm2_top6_0l"                               +myVarSuffix);
            const auto& fwm3_top6_0l                      = tr.getVar<double>("fwm3_top6_0l"                               +myVarSuffix);
            const auto& fwm4_top6_0l                      = tr.getVar<double>("fwm4_top6_0l"                               +myVarSuffix);
            const auto& fwm5_top6_0l                      = tr.getVar<double>("fwm5_top6_0l"                               +myVarSuffix);
            const auto& jmt_ev0_top6_0l                   = tr.getVar<double>("jmt_ev0_top6_0l"                            +myVarSuffix);
            const auto& jmt_ev1_top6_0l                   = tr.getVar<double>("jmt_ev1_top6_0l"                            +myVarSuffix);
            const auto& jmt_ev2_top6_0l                   = tr.getVar<double>("jmt_ev2_top6_0l"                            +myVarSuffix);
            const auto& fwm2_top6_1l                      = tr.getVar<double>("fwm2_top6_1l"                               +myVarSuffix);
            const auto& fwm3_top6_1l                      = tr.getVar<double>("fwm3_top6_1l"                               +myVarSuffix);
            const auto& fwm4_top6_1l                      = tr.getVar<double>("fwm4_top6_1l"                               +myVarSuffix);
            const auto& fwm5_top6_1l                      = tr.getVar<double>("fwm5_top6_1l"                               +myVarSuffix);
            const auto& jmt_ev0_top6_1l                   = tr.getVar<double>("jmt_ev0_top6_1l"                            +myVarSuffix);
            const auto& jmt_ev1_top6_1l                   = tr.getVar<double>("jmt_ev1_top6_1l"                            +myVarSuffix);
            const auto& jmt_ev2_top6_1l                   = tr.getVar<double>("jmt_ev2_top6_1l"                            +myVarSuffix);
            const auto& NonIsoMuons_fwm2_top6_1l          = tr.getVar<double>("NonIsoMuons_fwm2_top6_1l"                   +myVarSuffix);
            const auto& NonIsoMuons_fwm3_top6_1l          = tr.getVar<double>("NonIsoMuons_fwm3_top6_1l"                   +myVarSuffix);
            const auto& NonIsoMuons_fwm4_top6_1l          = tr.getVar<double>("NonIsoMuons_fwm4_top6_1l"                   +myVarSuffix);
            const auto& NonIsoMuons_fwm5_top6_1l          = tr.getVar<double>("NonIsoMuons_fwm5_top6_1l"                   +myVarSuffix);
            const auto& NonIsoMuons_jmt_ev0_top6_1l       = tr.getVar<double>("NonIsoMuons_jmt_ev0_top6_1l"                +myVarSuffix);
            const auto& NonIsoMuons_jmt_ev1_top6_1l       = tr.getVar<double>("NonIsoMuons_jmt_ev1_top6_1l"                +myVarSuffix);
            const auto& NonIsoMuons_jmt_ev2_top6_1l       = tr.getVar<double>("NonIsoMuons_jmt_ev2_top6_1l"                +myVarSuffix);

            const auto& Jets_cm_top6_0l                   = tr.getVec<utility::LorentzVector>("Jets_cm_top6_0l"            +myVarSuffix);
            const auto& Jets_cm_top6_1l                   = tr.getVec<utility::LorentzVector>("Jets_cm_top6_1l"            +myVarSuffix);
            const auto& NonIsoMuons_Jets_cm_top6_1l       = tr.getVec<utility::LorentzVector>("NonIsoMuons_Jets_cm_top6_1l"+myVarSuffix);
            const auto& nMVAJets_0l                       = tr.getVar<unsigned int>("nMVAJets_0l"                          +myVarSuffix);
            const auto& nMVAJets_1l                       = tr.getVar<unsigned int>("nMVAJets_1l"                          +myVarSuffix);

            const auto& eventCounter                      = tr.getVar<int>("eventCounter");
            const auto& Stop1_pt_cm_OldSeed               = tr.getVar<double>("Stop1_pt_cm_OldSeed"                        +myVarSuffix);
            const auto& Stop1_eta_cm_OldSeed              = tr.getVar<double>("Stop1_eta_cm_OldSeed"                       +myVarSuffix);
            const auto& Stop1_phi_cm_OldSeed              = tr.getVar<double>("Stop1_phi_cm_OldSeed"                       +myVarSuffix);
            const auto& Stop1_mass_cm_OldSeed             = tr.getVar<double>("Stop1_mass_cm_OldSeed"                      +myVarSuffix);
            const auto& Stop2_pt_cm_OldSeed               = tr.getVar<double>("Stop2_pt_cm_OldSeed"                        +myVarSuffix);
            const auto& Stop2_eta_cm_OldSeed              = tr.getVar<double>("Stop2_eta_cm_OldSeed"                       +myVarSuffix);
            const auto& Stop2_phi_cm_OldSeed              = tr.getVar<double>("Stop2_phi_cm_OldSeed"                       +myVarSuffix);
            const auto& Stop2_mass_cm_OldSeed             = tr.getVar<double>("Stop2_mass_cm_OldSeed"                      +myVarSuffix);
            const auto& Stop1_pt_cm_TopSeed               = tr.getVar<double>("Stop1_pt_cm_TopSeed"                        +myVarSuffix);
            const auto& Stop1_eta_cm_TopSeed              = tr.getVar<double>("Stop1_eta_cm_TopSeed"                       +myVarSuffix);
            const auto& Stop1_phi_cm_TopSeed              = tr.getVar<double>("Stop1_phi_cm_TopSeed"                       +myVarSuffix);
            const auto& Stop1_mass_cm_TopSeed             = tr.getVar<double>("Stop1_mass_cm_TopSeed"                      +myVarSuffix);
            const auto& Stop2_pt_cm_TopSeed               = tr.getVar<double>("Stop2_pt_cm_TopSeed"                        +myVarSuffix);
            const auto& Stop2_eta_cm_TopSeed              = tr.getVar<double>("Stop2_eta_cm_TopSeed"                       +myVarSuffix);
            const auto& Stop2_phi_cm_TopSeed              = tr.getVar<double>("Stop2_phi_cm_TopSeed"                       +myVarSuffix);
            const auto& Stop2_mass_cm_TopSeed             = tr.getVar<double>("Stop2_mass_cm_TopSeed"                      +myVarSuffix);

            const auto& Stop1_pt_cm_OldSeed_NonIsoMuon    = tr.getVar<double>("Stop1_pt_cm_OldSeed_NonIsoMuon"             +myVarSuffix);
            const auto& Stop1_eta_cm_OldSeed_NonIsoMuon   = tr.getVar<double>("Stop1_eta_cm_OldSeed_NonIsoMuon"            +myVarSuffix);
            const auto& Stop1_phi_cm_OldSeed_NonIsoMuon   = tr.getVar<double>("Stop1_phi_cm_OldSeed_NonIsoMuon"            +myVarSuffix);
            const auto& Stop1_mass_cm_OldSeed_NonIsoMuon  = tr.getVar<double>("Stop1_mass_cm_OldSeed_NonIsoMuon"           +myVarSuffix);
            const auto& Stop2_pt_cm_OldSeed_NonIsoMuon    = tr.getVar<double>("Stop2_pt_cm_OldSeed_NonIsoMuon"             +myVarSuffix);
            const auto& Stop2_eta_cm_OldSeed_NonIsoMuon   = tr.getVar<double>("Stop2_eta_cm_OldSeed_NonIsoMuon"            +myVarSuffix);
            const auto& Stop2_phi_cm_OldSeed_NonIsoMuon   = tr.getVar<double>("Stop2_phi_cm_OldSeed_NonIsoMuon"            +myVarSuffix);
            const auto& Stop2_mass_cm_OldSeed_NonIsoMuon  = tr.getVar<double>("Stop2_mass_cm_OldSeed_NonIsoMuon"           +myVarSuffix);
            const auto& Stop1_pt_cm_TopSeed_NonIsoMuon    = tr.getVar<double>("Stop1_pt_cm_TopSeed_NonIsoMuon"             +myVarSuffix);
            const auto& Stop1_eta_cm_TopSeed_NonIsoMuon   = tr.getVar<double>("Stop1_eta_cm_TopSeed_NonIsoMuon"            +myVarSuffix);
            const auto& Stop1_phi_cm_TopSeed_NonIsoMuon   = tr.getVar<double>("Stop1_phi_cm_TopSeed_NonIsoMuon"            +myVarSuffix);
            const auto& Stop1_mass_cm_TopSeed_NonIsoMuon  = tr.getVar<double>("Stop1_mass_cm_TopSeed_NonIsoMuon"           +myVarSuffix);
            const auto& Stop2_pt_cm_TopSeed_NonIsoMuon    = tr.getVar<double>("Stop2_pt_cm_TopSeed_NonIsoMuon"             +myVarSuffix);
            const auto& Stop2_eta_cm_TopSeed_NonIsoMuon   = tr.getVar<double>("Stop2_eta_cm_TopSeed_NonIsoMuon"            +myVarSuffix);
            const auto& Stop2_phi_cm_TopSeed_NonIsoMuon   = tr.getVar<double>("Stop2_phi_cm_TopSeed_NonIsoMuon"            +myVarSuffix);
            const auto& Stop2_mass_cm_TopSeed_NonIsoMuon  = tr.getVar<double>("Stop2_mass_cm_TopSeed_NonIsoMuon"           +myVarSuffix);

            const auto& GoodLeptons                       = tr.getVec<std::pair<std::string, utility::LorentzVector>>("GoodLeptons"+myVarSuffix);
            const auto& GoodLeptonsCharge                 = tr.getVec<int>("GoodLeptonsCharge"                                     +myVarSuffix);
            const auto& GoodLeptonsMiniIso                = tr.getVec<double>("GoodLeptonsMiniIso"                                 +myVarSuffix);

            const auto& regions_0l                        = tr.getVec<std::string>("regions_0l_RPV");
            const auto& regions_1l                        = tr.getVec<std::string>("regions_1l_RPV");

            const auto& Stop1_mass_PtRank_matched         = runtype == "Data" ? -999.0 : tr.getVar<float>("stop1_ptrank_mass"+myVarSuffix);
            const auto& Stop2_mass_PtRank_matched         = runtype == "Data" ? -999.0 : tr.getVar<float>("stop2_ptrank_mass"+myVarSuffix);
            const auto& Stop1_mass_MassRank_matched       = runtype == "Data" ? -999.0 : tr.getVar<float>("stop1_mrank_mass" +myVarSuffix);
            const auto& Stop2_mass_MassRank_matched       = runtype == "Data" ? -999.0 : tr.getVar<float>("stop2_mrank_mass" +myVarSuffix);
            const auto& Stop_mass_average_matched         = runtype == "Data" ? -999.0 : tr.getVar<double>("stop_avemass"    +myVarSuffix);

            std::map<std::string, std::vector<bool> > DoubleDisCo_passRegions_0l;
            std::map<std::string, std::vector<bool> > DoubleDisCo_passRegions_1l;
            std::map<std::string, std::vector<bool> > DoubleDisCo_passRegions_NonIsoMuon_0l;
            std::map<std::string, std::vector<bool> > DoubleDisCo_passRegions_NonIsoMuon_1l;

            for (const auto region : regions_0l) {
                DoubleDisCo_passRegions_0l[region]            = tr.getVec<bool>("DoubleDisCo_"+region+"_0l_RPV"+myVarSuffix); 
                DoubleDisCo_passRegions_NonIsoMuon_0l[region] = tr.getVec<bool>("DoubleDisCo_"+region+"_NonIsoMuon_0l_RPV"+myVarSuffix); 
            }
            for (const auto region : regions_1l) {
                DoubleDisCo_passRegions_1l[region]            = tr.getVec<bool>("DoubleDisCo_"+region+"_1l_RPV"+myVarSuffix); 
                DoubleDisCo_passRegions_NonIsoMuon_1l[region] = tr.getVec<bool>("DoubleDisCo_"+region+"_NonIsoMuon_1l_RPV"+myVarSuffix); 
            }

            std::vector<double> Jets_flavb_0l;   std::vector<double> Jets_flavb_1l;
            std::vector<double> Jets_flavc_0l;   std::vector<double> Jets_flavc_1l;
            std::vector<double> Jets_flavg_0l;   std::vector<double> Jets_flavg_1l;
            std::vector<double> Jets_flavuds_0l; std::vector<double> Jets_flavuds_1l;
            std::vector<double> Jets_flavq_0l;   std::vector<double> Jets_flavq_1l;
            std::vector<double> JetNonIsoMuons_flavb_1l;
            std::vector<double> JetNonIsoMuons_flavc_1l;
            std::vector<double> JetNonIsoMuons_flavg_1l;
            std::vector<double> JetNonIsoMuons_flavuds_1l;
            std::vector<double> JetNonIsoMuons_flavq_1l;

            // Here assume number cm jets is same in CR and SR selection
            for (unsigned int iJet = 1; iJet <= nMVAJets_0l; iJet++) {
                Jets_flavb_0l.push_back(tr.getVar<double>("Jet_flavb_"    +std::to_string(iJet)+"_0l"+myVarSuffix));
                Jets_flavc_0l.push_back(tr.getVar<double>("Jet_flavc_"    +std::to_string(iJet)+"_0l"+myVarSuffix));
                Jets_flavg_0l.push_back(tr.getVar<double>("Jet_flavg_"    +std::to_string(iJet)+"_0l"+myVarSuffix));
                Jets_flavuds_0l.push_back(tr.getVar<double>("Jet_flavuds_"+std::to_string(iJet)+"_0l"+myVarSuffix));
                Jets_flavq_0l.push_back(tr.getVar<double>("Jet_flavq_"    +std::to_string(iJet)+"_0l"+myVarSuffix));
            }

            for (unsigned int iJet = 1; iJet <= nMVAJets_1l; iJet++) {
                Jets_flavb_1l.push_back(tr.getVar<double>("Jet_flavb_"    +std::to_string(iJet)+"_1l"+myVarSuffix));
                Jets_flavc_1l.push_back(tr.getVar<double>("Jet_flavc_"    +std::to_string(iJet)+"_1l"+myVarSuffix));
                Jets_flavg_1l.push_back(tr.getVar<double>("Jet_flavg_"    +std::to_string(iJet)+"_1l"+myVarSuffix));
                Jets_flavuds_1l.push_back(tr.getVar<double>("Jet_flavuds_"+std::to_string(iJet)+"_1l"+myVarSuffix));
                Jets_flavq_1l.push_back(tr.getVar<double>("Jet_flavq_"    +std::to_string(iJet)+"_1l"+myVarSuffix));

                JetNonIsoMuons_flavb_1l.push_back(  tr.getVar<double>("JetNonIsoMuons_flavb_"  +std::to_string(iJet)+"_1l"+myVarSuffix));
                JetNonIsoMuons_flavc_1l.push_back(  tr.getVar<double>("JetNonIsoMuons_flavc_"  +std::to_string(iJet)+"_1l"+myVarSuffix));
                JetNonIsoMuons_flavg_1l.push_back(  tr.getVar<double>("JetNonIsoMuons_flavg_"  +std::to_string(iJet)+"_1l"+myVarSuffix));
                JetNonIsoMuons_flavuds_1l.push_back(tr.getVar<double>("JetNonIsoMuons_flavuds_"+std::to_string(iJet)+"_1l"+myVarSuffix));
                JetNonIsoMuons_flavq_1l.push_back(  tr.getVar<double>("JetNonIsoMuons_flavq_"  +std::to_string(iJet)+"_1l"+myVarSuffix));
            }

            // Put 0L and 1L version of variables into vector
            // 0th position is for 0L and 1st position is for 1L for convenience in the event loop
            // For 0L things in the CR, nominal version of jets and derived quantities are used, take note !
            std::vector<int>                                        NGoodJets                     {NGoodJets_pt30,                        NGoodJets_pt30                       };
            std::vector<int>                                        NNonIsoMuonJets               {NGoodJets_pt30,                        NNonIsoMuonJets_pt30                 };
            std::vector<double>                                     HT_trigger                    {HT_trigger_pt30,                       HT_trigger_pt30                      };
            std::vector<double>                                     HT_NonIsoMuon                 {HT_trigger_pt30,                       HT_NonIsoMuon_pt30                   };
            std::vector<double>                                     DoubleDisCo_massReg           {DoubleDisCo_massReg_0l,                DoubleDisCo_massReg_1l               };
            std::vector<double>                                     DoubleDisCo_disc1             {DoubleDisCo_disc1_0l,                  DoubleDisCo_disc1_1l                 };
            std::vector<double>                                     DoubleDisCo_disc2             {DoubleDisCo_disc2_0l,                  DoubleDisCo_disc2_1l                 };
            std::vector<double>                                     DoubleDisCo_QCDCR_massReg     {DoubleDisCo_massReg_NonIsoMuon_0l,     DoubleDisCo_massReg_NonIsoMuon_1l    };
            std::vector<double>                                     DoubleDisCo_QCDCR_disc1       {DoubleDisCo_disc1_NonIsoMuon_0l,       DoubleDisCo_disc1_NonIsoMuon_1l      };
            std::vector<double>                                     DoubleDisCo_QCDCR_disc2       {DoubleDisCo_disc2_NonIsoMuon_0l,       DoubleDisCo_disc2_NonIsoMuon_1l      };
            std::vector<double>                                     fwm2_top6                     {fwm2_top6_0l,                          fwm2_top6_1l                         };
            std::vector<double>                                     fwm3_top6                     {fwm3_top6_0l,                          fwm3_top6_1l                         };
            std::vector<double>                                     fwm4_top6                     {fwm4_top6_0l,                          fwm4_top6_1l                         };
            std::vector<double>                                     fwm5_top6                     {fwm5_top6_0l,                          fwm5_top6_1l                         };
            std::vector<double>                                     jmt_ev0_top6                  {jmt_ev0_top6_0l,                       jmt_ev0_top6_1l                      };
            std::vector<double>                                     jmt_ev1_top6                  {jmt_ev1_top6_0l,                       jmt_ev1_top6_1l                      };
            std::vector<double>                                     jmt_ev2_top6                  {jmt_ev2_top6_0l,                       jmt_ev2_top6_1l                      };
            std::vector<double>                                     fwm2_top6_QCDCR               {fwm2_top6_0l,                          NonIsoMuons_fwm2_top6_1l             };
            std::vector<double>                                     fwm3_top6_QCDCR               {fwm3_top6_0l,                          NonIsoMuons_fwm3_top6_1l             };
            std::vector<double>                                     fwm4_top6_QCDCR               {fwm4_top6_0l,                          NonIsoMuons_fwm4_top6_1l             };
            std::vector<double>                                     fwm5_top6_QCDCR               {fwm5_top6_0l,                          NonIsoMuons_fwm5_top6_1l             };
            std::vector<double>                                     jmt_ev0_top6_QCDCR            {jmt_ev0_top6_0l,                       NonIsoMuons_jmt_ev0_top6_1l          };
            std::vector<double>                                     jmt_ev1_top6_QCDCR            {jmt_ev1_top6_0l,                       NonIsoMuons_jmt_ev1_top6_1l          };
            std::vector<double>                                     jmt_ev2_top6_QCDCR            {jmt_ev2_top6_0l,                       NonIsoMuons_jmt_ev2_top6_1l          };
            std::vector<unsigned int>                               nMVAJets                      {nMVAJets_0l,                           nMVAJets_1l                          };
            std::vector<std::vector<double> >                       Jets_flavb                    {Jets_flavb_0l,                         Jets_flavb_1l                        };
            std::vector<std::vector<double> >                       Jets_flavc                    {Jets_flavc_0l,                         Jets_flavc_1l                        };
            std::vector<std::vector<double> >                       Jets_flavg                    {Jets_flavg_0l,                         Jets_flavg_1l                        };
            std::vector<std::vector<double> >                       Jets_flavuds                  {Jets_flavuds_0l,                       Jets_flavuds_1l                      };
            std::vector<std::vector<double> >                       Jets_flavq                    {Jets_flavq_0l,                         Jets_flavq_1l                        };
            std::vector<std::vector<double> >                       Jets_flavb_QCDCR              {Jets_flavb_0l,                         JetNonIsoMuons_flavb_1l              };
            std::vector<std::vector<double> >                       Jets_flavc_QCDCR              {Jets_flavc_0l,                         JetNonIsoMuons_flavc_1l              };
            std::vector<std::vector<double> >                       Jets_flavg_QCDCR              {Jets_flavg_0l,                         JetNonIsoMuons_flavg_1l              };
            std::vector<std::vector<double> >                       Jets_flavuds_QCDCR            {Jets_flavuds_0l,                       JetNonIsoMuons_flavuds_1l            };
            std::vector<std::vector<double> >                       Jets_flavq_QCDCR              {Jets_flavq_0l,                         JetNonIsoMuons_flavq_1l              };
            std::vector<std::vector<std::string> >                  regions                       {regions_0l,                            regions_1l                           };
            std::vector<std::vector<utility::LorentzVector> >       Jets_cm_top6                  {Jets_cm_top6_0l,                       Jets_cm_top6_1l                      };
            std::vector<std::vector<utility::LorentzVector> >       Jets_cm_top6_QCDCR            {Jets_cm_top6_0l,                       NonIsoMuons_Jets_cm_top6_1l          };
            std::vector<std::map<std::string, std::vector<bool> > > DoubleDisCo_passRegions       {DoubleDisCo_passRegions_0l,            DoubleDisCo_passRegions_1l           }; 
            std::vector<std::map<std::string, std::vector<bool> > > DoubleDisCo_passRegions_QCDCR {DoubleDisCo_passRegions_NonIsoMuon_0l, DoubleDisCo_passRegions_NonIsoMuon_1l}; 

            // ------------------------
            // -- Define weights
            // ------------------------
            double eventweight = 1.0, pileupWeight = 1.0, prefiringScaleFactor = 1.0, bTagWeight = 1.0, leptonweight = 1.0, hadronicWeight = 1.0;
            double weight0L             = 1.0, weight1L             = 1.0;
            double weight0L_NonIsoMuon  = 1.0, weight1L_NonIsoMuon  = 1.0;
            
            if(runtype == "MC" )
            {
                // common event weights
                const auto& Weight           = tr.getVar<float>("Weight"    );
                const auto& lumi             = tr.getVar<double>("FinalLumi");
                eventweight                  = lumi * Weight;
                pileupWeight                 = tr.getVar<double>("puWeightCorr"                    +myVarSuffix);
                prefiringScaleFactor         = tr.getVar<double>("prefiringScaleFactor"            +myVarSuffix);
                bTagWeight                   = tr.getVar<double>("bTagSF_EventWeightSimple_Central"+myVarSuffix);
                // leptonic event weights
                const auto& eleLepWeight     = tr.getVar<double>("totGoodElectronSF"               +myVarSuffix);
                const auto& muLepWeight      = tr.getVar<double>("totGoodMuonSF"                   +myVarSuffix);
                const auto& muNonIso         = tr.getVar<double>("totNonIsoMuonSF"                 +myVarSuffix);
                const auto& topPtScaleFactor = tr.getVar<double>("topPtScaleFactor"                +myVarSuffix);
                leptonweight                 = eleLepWeight * muLepWeight;
                // hadronic event weights
                hadronicWeight = tr.getVar<double>("jetTrigSF"+myVarSuffix);
                
                weight0L             *= eventweight * pileupWeight * prefiringScaleFactor * bTagWeight * hadronicWeight;                
                weight1L             *= eventweight * pileupWeight * prefiringScaleFactor * bTagWeight * leptonweight; 
                weight0L_NonIsoMuon  *= eventweight * pileupWeight * prefiringScaleFactor * muNonIso;
                weight1L_NonIsoMuon  *= eventweight * pileupWeight * prefiringScaleFactor * muNonIso;

                // tt variations have unique samples "except" FSR/ISR
                // FSR/ISR are event weight based variation 
                // So, calculate the event weights for FSR/ISR to make unique samples
                const auto& PSweight_FSRUp   = tr.getVar<double>("PSweight_FSRUp"  );
                const auto& PSweight_FSRDown = tr.getVar<double>("PSweight_FSRDown");
                const auto& PSweight_ISRUp   = tr.getVar<double>("PSweight_ISRUp"  );
                const auto& PSweight_ISRDown = tr.getVar<double>("PSweight_ISRDown");
 
                if (dataset.find("fsrUP") != std::string::npos)
                {
                    weight0L            *= PSweight_FSRUp;  
                    weight1L            *= PSweight_FSRUp;
                    weight1L_NonIsoMuon *= PSweight_FSRUp;
                    weight0L_NonIsoMuon *= PSweight_FSRUp;
                }                 
   
                else if (dataset.find("fsrDOWN") != std::string::npos)
                {
                    weight0L            *= PSweight_FSRDown; 
                    weight1L            *= PSweight_FSRDown;
                    weight1L_NonIsoMuon *= PSweight_FSRDown;
                    weight0L_NonIsoMuon *= PSweight_FSRDown;
                }               

                else if (dataset.find("isrUP") != std::string::npos)
                {
                    weight0L            *= PSweight_ISRUp; 
                    weight1L            *= PSweight_ISRUp;
                    weight1L_NonIsoMuon *= PSweight_ISRUp;
                    weight0L_NonIsoMuon *= PSweight_ISRUp;
                }

                else if (dataset.find("isrDOWN") != std::string::npos)
                {
                    weight0L            *= PSweight_ISRDown; 
                    weight1L            *= PSweight_ISRDown;
                    weight1L_NonIsoMuon *= PSweight_ISRDown;
                    weight0L_NonIsoMuon *= PSweight_ISRDown;
                } 

            }

            std::vector<double> weight       {weight0L,            weight1L           };
            std::vector<double> weight_QCDCR {weight0L_NonIsoMuon, weight1L_NonIsoMuon};

            const std::map<std::string, bool> cut_map
            {
                {"_0l"               , passBaseline0l_Good},
                //{"_0l_QCDCR"         , passBaseline0l_NonIsoMuon},                         
                {"_1l"               , passBaseline1l_Good},                         
                //{"_1l_QCDCR"         , passBaseline1l_NonIsoMuon}, 
            };

            std::map<std::string, bool>               njetsMap;
            std::map<std::string, std::vector<bool> > ABCDmap;

            if(!initHistos)
            {
                Preinit(nMVAJets_0l > nMVAJets_1l ? nMVAJets_0l : nMVAJets_1l);
                InitHistos(cut_map, regions);
                initHistos = true;
            }

            // Fill once per event
            my_histos["EventCounter"]->Fill(eventCounter);

            for(auto& kv : cut_map)
            {

                // Extract "0" or "1" from cut string e.g. _1l_QCDCR
                int channel = 1;
                if (kv.first.size() > 0 and kv.first.substr(2,1) == "l")
                    channel = std::stoi(kv.first.substr(1,1));

                bool isQCD          = kv.first.find("QCDCR")    != std::string::npos;
                bool isTopReweight  = kv.first.find("top")      != std::string::npos;
                unsigned int nJets = !isQCD ? Jets_cm_top6[channel].size() : Jets_cm_top6_QCDCR[channel].size();

                // For 0L, we always use the NGoodJets case
                njetsMap = {{"Incl",     true},
                    {"7",      (!isQCD or channel==0) ? NGoodJets[channel]==7  : NNonIsoMuonJets[channel]==7},
                    {"8",      (!isQCD or channel==0) ? NGoodJets[channel]==8  : NNonIsoMuonJets[channel]==8},
                    {"9",      (!isQCD or channel==0) ? NGoodJets[channel]==9  : NNonIsoMuonJets[channel]==9},
                    {"10",     (!isQCD or channel==0) ? NGoodJets[channel]==10 : NNonIsoMuonJets[channel]==10},
                    {"11",     (!isQCD or channel==0) ? NGoodJets[channel]==11 : NNonIsoMuonJets[channel]==11},
                    {"12",     (!isQCD or channel==0) ? NGoodJets[channel]==12 : NNonIsoMuonJets[channel]==12},
                    {"12incl", (!isQCD or channel==0) ? NGoodJets[channel]>=12 : NNonIsoMuonJets[channel]>=12},

                };

                ABCDmap = { {"Incl", {true}} };
                for (const auto& region : regions[channel])
                    ABCDmap[region] = !isQCD ? DoubleDisCo_passRegions[channel][region] : DoubleDisCo_passRegions_QCDCR[channel][region];

                double w  = !isQCD ? weight[channel]     : weight_QCDCR[channel];
                double ht = !isQCD ? HT_trigger[channel] : HT_NonIsoMuon[channel];

                std::string name;
                for (auto& njetsPass : njetsMap)
                {
                    std::string njets      = njetsPass.first;
                    bool        inNjetsBin = njetsPass.second;
                    std::string njetsStr = "";
                    if (njets != "Incl")  
                        njetsStr = "_Njets" + njets;

                    for (auto& regionPass : ABCDmap)
                    {
                        std::string       region      = regionPass.first;
                        std::vector<bool> inRegionBin = regionPass.second;
                        bool inRegion = false;
                        for (bool pass : inRegionBin) {
                            inRegion |= pass;
                        }

                        std::string regionStr = "";
                        if (region != "Incl")
                            regionStr = "_" + region;

                        if (kv.second and inNjetsBin and inRegion)
                        {
                            name = kv.first + njetsStr + regionStr;
                            if (njets == "Incl")
                            {

                                my_histos["h_njets"      + name]->Fill(!isQCD ? NGoodJets[channel]    : NNonIsoMuonJets[channel],    w);
                                my_histos["fwm2_top6"    + name]->Fill(!isQCD ? fwm2_top6[channel]    : fwm2_top6_QCDCR[channel],    w);
                                my_histos["fwm3_top6"    + name]->Fill(!isQCD ? fwm3_top6[channel]    : fwm3_top6_QCDCR[channel],    w);
                                my_histos["fwm4_top6"    + name]->Fill(!isQCD ? fwm4_top6[channel]    : fwm4_top6_QCDCR[channel],    w);
                                my_histos["fwm5_top6"    + name]->Fill(!isQCD ? fwm5_top6[channel]    : fwm5_top6_QCDCR[channel],    w);
                                my_histos["jmt_ev0_top6" + name]->Fill(!isQCD ? jmt_ev0_top6[channel] : jmt_ev0_top6_QCDCR[channel], w);
                                my_histos["jmt_ev1_top6" + name]->Fill(!isQCD ? jmt_ev1_top6[channel] : jmt_ev1_top6_QCDCR[channel], w);
                                my_histos["jmt_ev2_top6" + name]->Fill(!isQCD ? jmt_ev2_top6[channel] : jmt_ev2_top6_QCDCR[channel], w);

                                // Plots of stop 4-vector are made with pt-ranked stops
                                if (Stop1_pt_cm_OldSeed  > Stop2_pt_cm_OldSeed)
                                {
                                    my_histos["Stop1_pt_cm_OldSeed"   + name]->Fill((!isQCD or channel==0) ? Stop1_pt_cm_OldSeed   : Stop1_pt_cm_OldSeed_NonIsoMuon,   w);
                                    my_histos["Stop1_eta_cm_OldSeed"  + name]->Fill((!isQCD or channel==0) ? Stop1_eta_cm_OldSeed  : Stop1_eta_cm_OldSeed_NonIsoMuon,  w);
                                    my_histos["Stop1_phi_cm_OldSeed"  + name]->Fill((!isQCD or channel==0) ? Stop1_phi_cm_OldSeed  : Stop1_phi_cm_OldSeed_NonIsoMuon,  w);
                                    my_histos["Stop1_mass_cm_OldSeed" + name]->Fill((!isQCD or channel==0) ? Stop1_mass_cm_OldSeed : Stop1_mass_cm_OldSeed_NonIsoMuon, w);

                                    my_histos["Stop2_pt_cm_OldSeed"   + name]->Fill((!isQCD or channel==0) ? Stop2_pt_cm_OldSeed   : Stop2_pt_cm_OldSeed_NonIsoMuon,   w);
                                    my_histos["Stop2_eta_cm_OldSeed"  + name]->Fill((!isQCD or channel==0) ? Stop2_eta_cm_OldSeed  : Stop2_eta_cm_OldSeed_NonIsoMuon,  w);
                                    my_histos["Stop2_phi_cm_OldSeed"  + name]->Fill((!isQCD or channel==0) ? Stop2_phi_cm_OldSeed  : Stop2_phi_cm_OldSeed_NonIsoMuon,  w);
                                    my_histos["Stop2_mass_cm_OldSeed" + name]->Fill((!isQCD or channel==0) ? Stop2_mass_cm_OldSeed : Stop2_mass_cm_OldSeed_NonIsoMuon, w);
                                } else
                                {
                                    my_histos["Stop1_pt_cm_OldSeed"   + name]->Fill((!isQCD or channel==0) ? Stop2_pt_cm_OldSeed   : Stop2_pt_cm_OldSeed_NonIsoMuon,   w);
                                    my_histos["Stop1_eta_cm_OldSeed"  + name]->Fill((!isQCD or channel==0) ? Stop2_eta_cm_OldSeed  : Stop2_eta_cm_OldSeed_NonIsoMuon,  w);
                                    my_histos["Stop1_phi_cm_OldSeed"  + name]->Fill((!isQCD or channel==0) ? Stop2_phi_cm_OldSeed  : Stop2_phi_cm_OldSeed_NonIsoMuon,  w);
                                    my_histos["Stop1_mass_cm_OldSeed" + name]->Fill((!isQCD or channel==0) ? Stop2_mass_cm_OldSeed : Stop2_mass_cm_OldSeed_NonIsoMuon, w);

                                    my_histos["Stop2_pt_cm_OldSeed"   + name]->Fill((!isQCD or channel==0) ? Stop1_pt_cm_OldSeed   : Stop1_pt_cm_OldSeed_NonIsoMuon,   w);
                                    my_histos["Stop2_eta_cm_OldSeed"  + name]->Fill((!isQCD or channel==0) ? Stop1_eta_cm_OldSeed  : Stop1_eta_cm_OldSeed_NonIsoMuon,  w);
                                    my_histos["Stop2_phi_cm_OldSeed"  + name]->Fill((!isQCD or channel==0) ? Stop1_phi_cm_OldSeed  : Stop1_phi_cm_OldSeed_NonIsoMuon,  w);
                                    my_histos["Stop2_mass_cm_OldSeed" + name]->Fill((!isQCD or channel==0) ? Stop1_mass_cm_OldSeed : Stop1_mass_cm_OldSeed_NonIsoMuon, w);
                                }

                                if (Stop1_pt_cm_TopSeed  > Stop2_pt_cm_TopSeed)
                                {
                                    my_histos["Stop1_pt_cm_TopSeed"   + name]->Fill((!isQCD or channel==0) ? Stop1_pt_cm_TopSeed   : Stop1_pt_cm_TopSeed_NonIsoMuon,   w);
                                    my_histos["Stop1_eta_cm_TopSeed"  + name]->Fill((!isQCD or channel==0) ? Stop1_eta_cm_TopSeed  : Stop1_eta_cm_TopSeed_NonIsoMuon,  w);
                                    my_histos["Stop1_phi_cm_TopSeed"  + name]->Fill((!isQCD or channel==0) ? Stop1_phi_cm_TopSeed  : Stop1_phi_cm_TopSeed_NonIsoMuon,  w);
                                    my_histos["Stop1_mass_cm_TopSeed" + name]->Fill((!isQCD or channel==0) ? Stop1_mass_cm_TopSeed : Stop1_mass_cm_TopSeed_NonIsoMuon, w);

                                    my_histos["Stop2_pt_cm_TopSeed"   + name]->Fill((!isQCD or channel==0) ? Stop2_pt_cm_TopSeed   : Stop2_pt_cm_TopSeed_NonIsoMuon,   w);
                                    my_histos["Stop2_eta_cm_TopSeed"  + name]->Fill((!isQCD or channel==0) ? Stop2_eta_cm_TopSeed  : Stop2_eta_cm_TopSeed_NonIsoMuon,  w);
                                    my_histos["Stop2_phi_cm_TopSeed"  + name]->Fill((!isQCD or channel==0) ? Stop2_phi_cm_TopSeed  : Stop2_phi_cm_TopSeed_NonIsoMuon,  w);
                                    my_histos["Stop2_mass_cm_TopSeed" + name]->Fill((!isQCD or channel==0) ? Stop2_mass_cm_TopSeed : Stop2_mass_cm_TopSeed_NonIsoMuon, w);
                                } else
                                {
                                    my_histos["Stop1_pt_cm_TopSeed"   + name]->Fill((!isQCD or channel==0) ? Stop2_pt_cm_TopSeed   : Stop2_pt_cm_TopSeed_NonIsoMuon,   w);
                                    my_histos["Stop1_eta_cm_TopSeed"  + name]->Fill((!isQCD or channel==0) ? Stop2_eta_cm_TopSeed  : Stop2_eta_cm_TopSeed_NonIsoMuon,  w);
                                    my_histos["Stop1_phi_cm_TopSeed"  + name]->Fill((!isQCD or channel==0) ? Stop2_phi_cm_TopSeed  : Stop2_phi_cm_TopSeed_NonIsoMuon,  w);
                                    my_histos["Stop1_mass_cm_TopSeed" + name]->Fill((!isQCD or channel==0) ? Stop2_mass_cm_TopSeed : Stop2_mass_cm_TopSeed_NonIsoMuon, w);

                                    my_histos["Stop2_pt_cm_TopSeed"   + name]->Fill((!isQCD or channel==0) ? Stop1_pt_cm_TopSeed   : Stop1_pt_cm_TopSeed_NonIsoMuon,   w);
                                    my_histos["Stop2_eta_cm_TopSeed"  + name]->Fill((!isQCD or channel==0) ? Stop1_eta_cm_TopSeed  : Stop1_eta_cm_TopSeed_NonIsoMuon,  w);
                                    my_histos["Stop2_phi_cm_TopSeed"  + name]->Fill((!isQCD or channel==0) ? Stop1_phi_cm_TopSeed  : Stop1_phi_cm_TopSeed_NonIsoMuon,  w);
                                    my_histos["Stop2_mass_cm_TopSeed" + name]->Fill((!isQCD or channel==0) ? Stop1_mass_cm_TopSeed : Stop1_mass_cm_TopSeed_NonIsoMuon, w);
                                }

                                my_histos["Stop1_mass_PtRank_matched"   + name]->Fill(Stop1_mass_PtRank_matched,   w);
                                my_histos["Stop2_mass_PtRank_matched"   + name]->Fill(Stop2_mass_PtRank_matched,   w);
                                my_histos["Stop1_mass_MassRank_matched" + name]->Fill(Stop1_mass_MassRank_matched, w);
                                my_histos["Stop2_mass_MassRank_matched" + name]->Fill(Stop2_mass_MassRank_matched, w);
                                my_histos["Stop_mass_average_matched"   + name]->Fill(Stop_mass_average_matched,   w);

                                for(unsigned int i = 1; i <= nJets; i++)
                                {
                                    double pt = 0.0, eta = 0.0, phi = 0.0, m = 0.0, E = 0.0;

                                    if (!isQCD)
                                    {
                                        pt  = static_cast<double>(Jets_cm_top6[channel].at(i-1).Pt() );
                                        eta = static_cast<double>(Jets_cm_top6[channel].at(i-1).Eta());
                                        phi = static_cast<double>(Jets_cm_top6[channel].at(i-1).Phi());
                                        m   = static_cast<double>(Jets_cm_top6[channel].at(i-1).M()  );
                                        E   = static_cast<double>(Jets_cm_top6[channel].at(i-1).E()  );
                                    } else
                                    {
                                        pt  = static_cast<double>(Jets_cm_top6_QCDCR[channel].at(i-1).Pt() );
                                        eta = static_cast<double>(Jets_cm_top6_QCDCR[channel].at(i-1).Eta());
                                        phi = static_cast<double>(Jets_cm_top6_QCDCR[channel].at(i-1).Phi());
                                        m   = static_cast<double>(Jets_cm_top6_QCDCR[channel].at(i-1).M()  );
                                        E   = static_cast<double>(Jets_cm_top6_QCDCR[channel].at(i-1).E()  );
                                    } 

                                    my_histos["Jet_cm_ptrHT_" + std::to_string(i) + name]->Fill(pt/ht, w);
                                    my_histos["Jet_cm_pt_"    + std::to_string(i) + name]->Fill(pt,    w);
                                    my_histos["Jet_cm_eta_"   + std::to_string(i) + name]->Fill(eta,   w);
                                    my_histos["Jet_cm_phi_"   + std::to_string(i) + name]->Fill(phi,   w);
                                    my_histos["Jet_cm_m_"     + std::to_string(i) + name]->Fill(m,     w);
                                    my_histos["Jet_cm_E_"     + std::to_string(i) + name]->Fill(E,     w);

                                    if (!isQCD)
                                    {
                                        my_histos["Jet_cm_flavb_"   + std::to_string(i) + name]->Fill(Jets_flavb[channel].at(i-1),   w);
                                        my_histos["Jet_cm_flavc_"   + std::to_string(i) + name]->Fill(Jets_flavc[channel].at(i-1),   w);
                                        my_histos["Jet_cm_flavg_"   + std::to_string(i) + name]->Fill(Jets_flavg[channel].at(i-1),   w);
                                        my_histos["Jet_cm_flavq_"   + std::to_string(i) + name]->Fill(Jets_flavq[channel].at(i-1),   w);
                                        my_histos["Jet_cm_flavuds_" + std::to_string(i) + name]->Fill(Jets_flavuds[channel].at(i-1), w);
                                    } else
                                    {
                                        my_histos["Jet_cm_flavb_"   + std::to_string(i) + name]->Fill(Jets_flavb_QCDCR[channel].at(i-1),   w);
                                        my_histos["Jet_cm_flavc_"   + std::to_string(i) + name]->Fill(Jets_flavc_QCDCR[channel].at(i-1),   w);
                                        my_histos["Jet_cm_flavg_"   + std::to_string(i) + name]->Fill(Jets_flavg_QCDCR[channel].at(i-1),   w);
                                        my_histos["Jet_cm_flavq_"   + std::to_string(i) + name]->Fill(Jets_flavq_QCDCR[channel].at(i-1),   w);
                                        my_histos["Jet_cm_flavuds_" + std::to_string(i) + name]->Fill(Jets_flavuds_QCDCR[channel].at(i-1), w);
                                    }
                                }

                                my_histos["h_ht"                  + name]->Fill(ht, w);
                                   for(std::size_t j = 0 ; j < std::size(GoodLeptons); ++j){
                                    auto& type = GoodLeptons[j].first;
                                    auto& lvec = GoodLeptons[j].second; 
                                    auto& charge = GoodLeptonsCharge[j];
                                    auto&  iso = GoodLeptonsMiniIso[j];

                                    my_histos["h_lPt"      + name]->Fill(lvec.Pt(),  w);
                                    my_histos["h_lPhi"     + name]->Fill(lvec.Phi(), w);
                                    my_histos["h_lEta"     + name]->Fill(lvec.Eta(), w);
                                    my_histos["h_lCharge"  + name]->Fill(charge ,    w);
                                    my_histos["h_lMiniIso" + name]->Fill(iso ,       w);
                                    if(type == 'e'){
                                    my_histos["h_ePt"  + name]->Fill(lvec.Pt() , w);
                                    my_histos["h_ePhi" + name]->Fill(lvec.Phi(), w);
                                    my_histos["h_eEta" + name]->Fill(lvec.Eta(), w);
                                    my_histos["h_eCharge" + name]->Fill(charge,  w);
                                    my_histos["h_eMiniIso" + name]->Fill(iso,    w);
                                    } else if (type == 'm') {
                                    my_histos["h_mPt"  + name]->Fill(lvec.Pt(),  w);
                                    my_histos["h_mPhi" + name]->Fill(lvec.Phi(), w);
                                    my_histos["h_mEta" + name]->Fill(lvec.Eta(), w);
                                    my_histos["h_mCharge" + name]->Fill(charge,  w);
                                    my_histos["h_mMiniIso" + name]->Fill(iso,    w);
                                            }
                                       }
                                
                                for(unsigned int j = 0; j < Jets.size(); j++) 
                                {
                                    if(!GoodJets_pt30[j]) continue;
                                    my_histos["h_jPt"  + name]->Fill(Jets.at(j).Pt(),  w); 
                                    my_histos["h_jEta" + name]->Fill(Jets.at(j).Eta(), w); 
                                    my_histos["h_jPhi" + name]->Fill(Jets.at(j).Phi(), w); 
                                    my_histos["h_jM"   + name]->Fill(Jets.at(j).M(),   w); 
                                    
                                }
                            }

                            for(unsigned int i = 1; i <= nJets; i++)
                            {
                                double pt = 0.0, eta = 0.0, phi = 0.0, m = 0.0, E = 0.0;
                                
                                if (!isQCD)
                                {
                                    pt  = static_cast<double>(Jets_cm_top6[channel].at(i-1).Pt());
                                    //eta = static_cast<double>(Jets_cm_top6[channel].at(i-1).Eta());
                                    //phi = static_cast<double>(Jets_cm_top6[channel].at(i-1).Phi());
                                    //m   = static_cast<double>(Jets_cm_top6[channel].at(i-1).M());
                                    //E   = static_cast<double>(Jets_cm_top6[channel].at(i-1).E());
                                } else
                                {
                                    pt  = static_cast<double>(Jets_cm_top6_QCDCR[channel].at(i-1).Pt());
                                    //eta = static_cast<double>(Jets_cm_top6_QCDCR[channel].at(i-1).Eta());
                                    //phi = static_cast<double>(Jets_cm_top6_QCDCR[channel].at(i-1).Phi());
                                    //m   = static_cast<double>(Jets_cm_top6_QCDCR[channel].at(i-1).M());
                                    //E   = static_cast<double>(Jets_cm_top6_QCDCR[channel].at(i-1).E());
                                } 
                                
                                my_2d_histos["h_cm_pt_jetRank"    + name]->Fill(pt,    i, w);
                                my_2d_histos["h_cm_ptrHT_jetRank" + name]->Fill(pt/ht, i, w);
                            }

                            if (region != "Incl" and njets != "Incl")
                                my_histos["h_DoubleDisCo_massReg" + name]->Fill(!isQCD ? DoubleDisCo_massReg[channel] : DoubleDisCo_QCDCR_massReg[channel], w);

                            // if plotting disco, no need to make plots when cutting on
                            if (njets != "Incl")
                            {
                                if (!isQCD)
                                {
                                    my_histos["h_DoubleDisCo_disc1"          + name]->Fill(DoubleDisCo_disc1[channel], w);
                                    my_histos["h_DoubleDisCo_disc2"          + name]->Fill(DoubleDisCo_disc2[channel], w);
                                } else
                                {
                                    my_histos["h_DoubleDisCo_disc1"          + name]->Fill(DoubleDisCo_QCDCR_disc1[channel], w);
                                    my_histos["h_DoubleDisCo_disc2"          + name]->Fill(DoubleDisCo_QCDCR_disc2[channel], w);
                                }
                            }

                            if (njets != "Incl")
                            {
                                if (!isQCD)
                                    my_2d_histos["h_DoubleDisCo_disc1_disc2" + name]->Fill(DoubleDisCo_disc1[channel], DoubleDisCo_disc2[channel], w);
                                else
                                    my_2d_histos["h_DoubleDisCo_disc1_disc2" + name]->Fill(DoubleDisCo_QCDCR_disc1[channel], DoubleDisCo_QCDCR_disc2[channel], w);
                            }
                        }

                        for (unsigned int iSubRegion = 0; iSubRegion < inRegionBin.size(); iSubRegion++)
                        {
                            if (kv.second and inNjetsBin and inRegion and inRegionBin[iSubRegion])
                            {
                                name = kv.first + njetsStr + regionStr;
                                int shift = inRegionBin[iSubRegion] ? iSubRegion : 100;

                                // if plotting Njets, don't care about individual Njet cuts
                                if (njets == "Incl" and region != "Incl")
                                {
                                    if (!isQCD or channel==0)
                                    {
                                        my_histos["h_njets_12incl" + name]->Fill(NGoodJets[channel]>=12       ? 12-nMVAJets[channel]+shift*6 : NGoodJets[channel]-nMVAJets[channel]+shift*6, w);

                                    } else
                                    {
                                        my_histos["h_njets_12incl" + name]->Fill(NNonIsoMuonJets[channel]>=12 ? 12-nMVAJets[channel]+shift*6 : NNonIsoMuonJets[channel]-nMVAJets[channel]+shift*6, w);

                                    }
                                }
                                if (njets != "Incl")
                                {
                                    name = kv.first + njetsStr + regionStr + "_" + subRegionsMap[region][iSubRegion];
                                    if (!isQCD)
                                        my_2d_histos["h_DoubleDisCo_disc1_disc2" + name]->Fill(DoubleDisCo_disc1[channel], DoubleDisCo_disc2[channel], w);
                                    else
                                        my_2d_histos["h_DoubleDisCo_disc1_disc2" + name]->Fill(DoubleDisCo_QCDCR_disc1[channel], DoubleDisCo_QCDCR_disc2[channel], w);

                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void AnalyzeDoubleDisCo::WriteHistos(TFile* outfile)
{
    outfile->cd();

    for(const auto& p : my_histos) 
        p.second->Write();

    for(const auto& p : my_2d_histos) 
        p.second->Write();
}
