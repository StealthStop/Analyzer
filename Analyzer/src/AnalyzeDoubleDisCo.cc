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
        {"h_njets_12incl",               24, -0.5, 23.5},
        {"h_ht",                        500,    0, 5000},
        {"h_dRbjets",                   320,    0,    6},
        {"h_mbl",                       320,    0,  320},
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
        {"h_DoubleDisCo_disc1_disc2", 100, 0,    1, 100, 0,  1}, 
        {"h_cm_pt_jetRank",           150, 0, 1500,  10, 0, 10},
        {"h_cm_ptrHT_jetRank",        150, 0,    1,  10, 0, 10},
    };

    njets = {"Incl", "7", "8", "9", "10", "11", "12incl"};

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

            // Attention to subregion names for bdEF or BD
            if (region.find("bd") != std::string::npos or region.find("BD") != std::string::npos) {
                subRegionsMap[region].push_back("B");
                subRegionsMap[region].push_back("E");
                subRegionsMap[region].push_back("D");
                subRegionsMap[region].push_back("F");
            // Subregion names for cdiGH or CDGH
            } else if (region.find("cd") != std::string::npos or region.find("CD") != std::string::npos) {
                subRegionsMap[region].push_back("C");
                subRegionsMap[region].push_back("D");
                subRegionsMap[region].push_back("G");
                subRegionsMap[region].push_back("H");
            // Subregion names for subDivD or D
            } else if (region.find("subDiv") != std::string::npos or region.find("_D") != std::string::npos) {
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
        MakeMVAVariables    makeMVAVariables0L_NonIsoMuon(false, myVarSuffix, "GoodJets_pt30",        false, true, 7,  1, "_0l");
        MakeMVAVariables    makeMVAVariables1L_NonIsoMuon(false, myVarSuffix, "NonIsoMuonJets_pt30",  false, true, 7,  1, "_1l");
        MakeMVAVariables    makeMVAVariables2L_NonIsoMuon(false, myVarSuffix, "NonIsoMuonJets_pt30",  false, true, 7,  1, "_2l");
        MakeMVAVariables    makeMVAVariables0L(false, myVarSuffix, "GoodJets_pt30", false, true, 7, 0, "_0l");
        MakeMVAVariables    makeMVAVariables1L(false, myVarSuffix, "GoodJets_pt30", false, true, 7, 1, "_1l");
        MakeMVAVariables    makeMVAVariables2L(false, myVarSuffix, "GoodJets_pt30", false, true, 6, 2, "_2l");
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
        tr.registerFunction(makeMVAVariables2L_NonIsoMuon);
        tr.registerFunction(makeMVAVariables0L);
        tr.registerFunction(makeMVAVariables1L);
        tr.registerFunction(makeMVAVariables2L);
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

        std::vector<std::string> channels = {"0", "1", "2"};

        for(const auto& myVarSuffix : my_var_suffix)
        {
            const auto& runtype      = tr.getVar<std::string>("runtype");    
            const auto& eventCounter = tr.getVar<int>("eventCounter");

            // Put 0L and 1L version of variables into vector
            // 0th position is for 0L and 1st position is for 1L for convenience in the event loop
            // For 0L things in the CR, nominal version of jets and derived quantities are used, take note !
            std::vector<std::vector<utility::LorentzVector> >       Jets_CR                        ;
            std::vector<bool>                                       Baseline_CR                    ;
            std::vector<double>                                     weight_CR                      ;
            std::vector<std::vector<bool> >                         GoodJets_CR                    ;
            std::vector<int>                                        NGoodJets_CR                   ;
            std::vector<double>                                     Mbl_CR                         ;
            std::vector<double>                                     dRbjets_CR                     ;
            std::vector<double>                                     HT_pt30_CR                     ;
            std::vector<double>                                     DoubleDisCo_massReg_CR         ;
            std::vector<double>                                     DoubleDisCo_disc1_CR           ;
            std::vector<double>                                     DoubleDisCo_disc2_CR           ;
            std::vector<double>                                     fwm2_top6_CR                   ;
            std::vector<double>                                     fwm3_top6_CR                   ;
            std::vector<double>                                     fwm4_top6_CR                   ;
            std::vector<double>                                     fwm5_top6_CR                   ;
            std::vector<double>                                     jmt_ev0_top6_CR                ;
            std::vector<double>                                     jmt_ev1_top6_CR                ;
            std::vector<double>                                     jmt_ev2_top6_CR                ;
            std::vector<unsigned int>                               nMVAJets_CR                    ;
            std::vector<std::vector<double> >                       Jets_flavb_CR                  ;
            std::vector<std::vector<double> >                       Jets_flavc_CR                  ;
            std::vector<std::vector<double> >                       Jets_flavg_CR                  ;
            std::vector<std::vector<double> >                       Jets_flavuds_CR                ;
            std::vector<std::vector<double> >                       Jets_flavq_CR                  ;
            std::vector<std::vector<std::string> >                  regions_CR                     ;
            std::vector<std::vector<utility::LorentzVector> >       Jets_cm_top6_CR                ;
            std::vector<std::map<std::string, std::vector<bool> > > DoubleDisCo_passRegions_CR     ; 

            std::vector<double>                                     Stop1_pt_cm_OldSeed_CR         ;
            std::vector<double>                                     Stop1_eta_cm_OldSeed_CR        ;
            std::vector<double>                                     Stop1_phi_cm_OldSeed_CR        ;
            std::vector<double>                                     Stop1_mass_cm_OldSeed_CR       ;
            std::vector<double>                                     Stop2_pt_cm_OldSeed_CR         ;
            std::vector<double>                                     Stop2_eta_cm_OldSeed_CR        ;
            std::vector<double>                                     Stop2_phi_cm_OldSeed_CR        ;
            std::vector<double>                                     Stop2_mass_cm_OldSeed_CR       ;
            std::vector<double>                                     Stop1_pt_cm_TopSeed_CR         ;
            std::vector<double>                                     Stop1_eta_cm_TopSeed_CR        ;
            std::vector<double>                                     Stop1_phi_cm_TopSeed_CR        ;
            std::vector<double>                                     Stop1_mass_cm_TopSeed_CR       ;
            std::vector<double>                                     Stop2_pt_cm_TopSeed_CR         ;
            std::vector<double>                                     Stop2_eta_cm_TopSeed_CR        ;
            std::vector<double>                                     Stop2_phi_cm_TopSeed_CR        ;
            std::vector<double>                                     Stop2_mass_cm_TopSeed_CR       ;

            std::vector<std::vector<std::pair<std::string, utility::LorentzVector> > > GoodLeptons_CR ;
            std::vector<std::vector<int> > GoodLeptonsCharge_CR                                    ;
            std::vector<std::vector<double> > GoodLeptonsMiniIso_CR                                ;

            std::vector<double> Stop1_mass_PtRank_matched_CR                                       ;
            std::vector<double> Stop2_mass_PtRank_matched_CR                                       ;
            std::vector<double> Stop1_mass_MassRank_matched_CR                                     ;
            std::vector<double> Stop2_mass_MassRank_matched_CR                                     ;
            std::vector<double> Stop_mass_average_matched_CR                                       ;

            std::vector<std::vector<utility::LorentzVector> >       Jets                           ;
            std::vector<bool>                                       Baseline                       ;
            std::vector<double>                                     weight                         ;
            std::vector<std::vector<bool> >                         GoodJets                       ;
            std::vector<int>                                        NGoodJets                      ;
            std::vector<double>                                     HT_pt30                        ;
            std::vector<double>                                     Mbl                            ;
            std::vector<double>                                     dRbjets                        ;
            std::vector<double>                                     DoubleDisCo_massReg            ;
            std::vector<double>                                     DoubleDisCo_disc1              ;
            std::vector<double>                                     DoubleDisCo_disc2              ;
            std::vector<double>                                     fwm2_top6                      ;
            std::vector<double>                                     fwm3_top6                      ;
            std::vector<double>                                     fwm4_top6                      ;
            std::vector<double>                                     fwm5_top6                      ;
            std::vector<double>                                     jmt_ev0_top6                   ;
            std::vector<double>                                     jmt_ev1_top6                   ;
            std::vector<double>                                     jmt_ev2_top6                   ;
            std::vector<unsigned int>                               nMVAJets                       ;
            std::vector<std::vector<double> >                       Jets_flavb                     ;
            std::vector<std::vector<double> >                       Jets_flavc                     ;
            std::vector<std::vector<double> >                       Jets_flavg                     ;
            std::vector<std::vector<double> >                       Jets_flavuds                   ;
            std::vector<std::vector<double> >                       Jets_flavq                     ;
            std::vector<std::vector<std::string> >                  regions                        ;
            std::vector<std::vector<utility::LorentzVector> >       Jets_cm_top6                   ;
            std::vector<std::map<std::string, std::vector<bool> > > DoubleDisCo_passRegions        ; 

            std::vector<double>                                     Stop1_pt_cm_OldSeed            ;
            std::vector<double>                                     Stop1_eta_cm_OldSeed           ;
            std::vector<double>                                     Stop1_phi_cm_OldSeed           ;
            std::vector<double>                                     Stop1_mass_cm_OldSeed          ;
            std::vector<double>                                     Stop2_pt_cm_OldSeed            ;
            std::vector<double>                                     Stop2_eta_cm_OldSeed           ;
            std::vector<double>                                     Stop2_phi_cm_OldSeed           ;
            std::vector<double>                                     Stop2_mass_cm_OldSeed          ;
            std::vector<double>                                     Stop1_pt_cm_TopSeed            ;
            std::vector<double>                                     Stop1_eta_cm_TopSeed           ;
            std::vector<double>                                     Stop1_phi_cm_TopSeed           ;
            std::vector<double>                                     Stop1_mass_cm_TopSeed          ;
            std::vector<double>                                     Stop2_pt_cm_TopSeed            ;
            std::vector<double>                                     Stop2_eta_cm_TopSeed           ;
            std::vector<double>                                     Stop2_phi_cm_TopSeed           ;
            std::vector<double>                                     Stop2_mass_cm_TopSeed          ;

            std::vector<std::vector<std::pair<std::string, utility::LorentzVector> > > GoodLeptons ;
            std::vector<std::vector<int> > GoodLeptonsCharge                                       ;
            std::vector<std::vector<double> > GoodLeptonsMiniIso                                   ;

            std::vector<double> Stop1_mass_PtRank_matched                                          ;
            std::vector<double> Stop2_mass_PtRank_matched                                          ;
            std::vector<double> Stop1_mass_MassRank_matched                                        ;
            std::vector<double> Stop2_mass_MassRank_matched                                        ;
            std::vector<double> Stop_mass_average_matched                                          ;

            for (auto& channel : channels)
            {

                // Collection names for 0L and 2L do  not change when in the QCDCR, because of the jet collection used...
                std::string mvaName  = "NonIsoMuons_";
                std::string flavName = "NonIsoMuons";
                std::string jetsName = "NonIsoMuon";
                std::string htName   = "NonIsoMuon";
                if (channel != "1")
                {
                    mvaName  = "";
                    flavName = "";
                    jetsName = "Good";
                    htName   = "trigger";
                }

                std::string tag = "_" + channel + "l";

                Jets.push_back(tr.getVec<utility::LorentzVector>("Jets"                + myVarSuffix));
                GoodJets.push_back(tr.getVec<bool>("GoodJets_pt30"                     + myVarSuffix));
                NGoodJets.push_back(tr.getVar<int>("NGoodJets_pt30"                    + myVarSuffix));
                HT_pt30_CR.push_back(tr.getVar<double>("HT_trigger_pt30"               + myVarSuffix));
                Baseline.push_back(tr.getVar<bool>("passBaseline" + channel + "l_Good" + myVarSuffix));
                Mbl.push_back(tr.getVar<double>("Mbl"                                  + myVarSuffix));
                dRbjets.push_back(tr.getVar<double>("dR_bjets"                         + myVarSuffix));

                if (channel != "2")
                {
                    DoubleDisCo_massReg.push_back(tr.getVar<double>("DoubleDisCo_massReg_" + channel + "l_RPV" + myVarSuffix));
                    DoubleDisCo_disc1.push_back(tr.getVar<double>("DoubleDisCo_disc1_" + channel + "l_RPV"     + myVarSuffix));
                    DoubleDisCo_disc2.push_back(tr.getVar<double>("DoubleDisCo_disc2_" + channel + "l_RPV"     + myVarSuffix));
                } else
                {
                    DoubleDisCo_massReg.push_back(tr.getVar<double>("DoubleDisCo_massReg_1l_RPV" + myVarSuffix));
                    DoubleDisCo_disc1.push_back(tr.getVar<double>("DoubleDisCo_disc1_1l_RPV"     + myVarSuffix));
                    DoubleDisCo_disc2.push_back(tr.getVar<double>("DoubleDisCo_disc2_1l_RPV"     + myVarSuffix));
                }

                fwm2_top6.push_back(tr.getVar<double>("fwm2_top6" + tag          + myVarSuffix));
                fwm3_top6.push_back(tr.getVar<double>("fwm3_top6" + tag          + myVarSuffix));
                fwm4_top6.push_back(tr.getVar<double>("fwm4_top6" + tag          + myVarSuffix));
                fwm5_top6.push_back(tr.getVar<double>("fwm5_top6" + tag          + myVarSuffix));
                jmt_ev0_top6.push_back(tr.getVar<double>("jmt_ev0_top6" + tag    + myVarSuffix));
                jmt_ev1_top6.push_back(tr.getVar<double>("jmt_ev1_top6" + tag    + myVarSuffix));
                jmt_ev2_top6.push_back(tr.getVar<double>("jmt_ev2_top6" + tag    + myVarSuffix));
                Jets_cm_top6.push_back(tr.getVec<utility::LorentzVector>("Jets_cm_top6" + tag            +myVarSuffix));
                nMVAJets.push_back(tr.getVar<unsigned int>("nMVAJets" + tag      + myVarSuffix));
                Stop1_pt_cm_OldSeed.push_back(tr.getVar<double>("Stop1_pt_cm_OldSeed"       + myVarSuffix));
                Stop1_eta_cm_OldSeed.push_back(tr.getVar<double>("Stop1_eta_cm_OldSeed"     + myVarSuffix));
                Stop1_phi_cm_OldSeed.push_back(tr.getVar<double>("Stop1_phi_cm_OldSeed"     + myVarSuffix));
                Stop1_mass_cm_OldSeed.push_back(tr.getVar<double>("Stop1_mass_cm_OldSeed"   + myVarSuffix));
                Stop2_pt_cm_OldSeed.push_back(tr.getVar<double>("Stop2_pt_cm_OldSeed"       + myVarSuffix));
                Stop2_eta_cm_OldSeed.push_back(tr.getVar<double>("Stop2_eta_cm_OldSeed"     + myVarSuffix));
                Stop2_phi_cm_OldSeed.push_back(tr.getVar<double>("Stop2_phi_cm_OldSeed"     + myVarSuffix));
                Stop2_mass_cm_OldSeed.push_back(tr.getVar<double>("Stop2_mass_cm_OldSeed"   + myVarSuffix));
                Stop1_pt_cm_TopSeed.push_back(tr.getVar<double>("Stop1_pt_cm_TopSeed"       + myVarSuffix));
                Stop1_eta_cm_TopSeed.push_back(tr.getVar<double>("Stop1_eta_cm_TopSeed"     + myVarSuffix));
                Stop1_phi_cm_TopSeed.push_back(tr.getVar<double>("Stop1_phi_cm_TopSeed"     + myVarSuffix));
                Stop1_mass_cm_TopSeed.push_back(tr.getVar<double>("Stop1_mass_cm_TopSeed"   + myVarSuffix));
                Stop2_pt_cm_TopSeed.push_back(tr.getVar<double>("Stop2_pt_cm_TopSeed"       + myVarSuffix));
                Stop2_eta_cm_TopSeed.push_back(tr.getVar<double>("Stop2_eta_cm_TopSeed"     + myVarSuffix));
                Stop2_phi_cm_TopSeed.push_back(tr.getVar<double>("Stop2_phi_cm_TopSeed"     + myVarSuffix));
                Stop2_mass_cm_TopSeed.push_back(tr.getVar<double>("Stop2_mass_cm_TopSeed"   + myVarSuffix));

                GoodLeptons.push_back(tr.getVec<std::pair<std::string, utility::LorentzVector>>("GoodLeptons" + myVarSuffix));
                GoodLeptonsCharge.push_back(tr.getVec<int>("GoodLeptonsCharge"                                + myVarSuffix));
                GoodLeptonsMiniIso.push_back(tr.getVec<double>("GoodLeptonsMiniIso"                           + myVarSuffix));

                Stop1_mass_PtRank_matched.push_back(  runtype == "Data" ? -999.0 : tr.getVar<float>("stop1_ptrank_mass"+myVarSuffix));
                Stop2_mass_PtRank_matched.push_back(  runtype == "Data" ? -999.0 : tr.getVar<float>("stop2_ptrank_mass"+myVarSuffix));
                Stop1_mass_MassRank_matched.push_back(runtype == "Data" ? -999.0 : tr.getVar<float>("stop1_mrank_mass" +myVarSuffix));
                Stop2_mass_MassRank_matched.push_back(runtype == "Data" ? -999.0 : tr.getVar<float>("stop2_mrank_mass" +myVarSuffix));
                Stop_mass_average_matched.push_back(  runtype == "Data" ? -999.0 : tr.getVar<double>("stop_avemass"    +myVarSuffix));

                std::vector<double> tempJets_flavb;  
                std::vector<double> tempJets_flavc;  
                std::vector<double> tempJets_flavg;  
                std::vector<double> tempJets_flavuds;
                std::vector<double> tempJets_flavq;  
                // Here assume number cm jets is same in CR and SR selection
                for (unsigned int iJet = 1; iJet <= nMVAJets.back(); iJet++) {
                    tempJets_flavb.push_back(tr.getVar<double>("Jet_flavb_"    +std::to_string(iJet)+ tag+myVarSuffix));
                    tempJets_flavc.push_back(tr.getVar<double>("Jet_flavc_"    +std::to_string(iJet)+ tag+myVarSuffix));
                    tempJets_flavg.push_back(tr.getVar<double>("Jet_flavg_"    +std::to_string(iJet)+ tag+myVarSuffix));
                    tempJets_flavuds.push_back(tr.getVar<double>("Jet_flavuds_"+std::to_string(iJet)+ tag+myVarSuffix));
                    tempJets_flavq.push_back(tr.getVar<double>("Jet_flavq_"    +std::to_string(iJet)+ tag+myVarSuffix));
                }

                Jets_flavb.push_back(tempJets_flavb);
                Jets_flavc.push_back(tempJets_flavc);
                Jets_flavg.push_back(tempJets_flavg);
                Jets_flavuds.push_back(tempJets_flavuds);
                Jets_flavq.push_back(tempJets_flavq);

                std::map<std::string, std::vector<bool> > tempRegionMap;
                std::map<std::string, std::vector<bool> > tempRegionMap_CR;

                if (channel != "2")
                { 
                    regions.push_back(tr.getVec<std::string>("regions_" + channel + "l_RPV"));
                } else
                {
                    regions.push_back(tr.getVec<std::string>("regions_1l_RPV"));
                }
                for (const std::string& region : regions.back())
                {
                    if (channel != "2")
                    {
                        tempRegionMap[region]    = tr.getVec<bool>("DoubleDisCo_"+region+"_" + channel + "l_RPV"+myVarSuffix); 
                        tempRegionMap_CR[region] = tr.getVec<bool>("DoubleDisCo_"+region+"_NonIsoMuon_" + channel + "l_RPV"+myVarSuffix); 
                    } else
                    {
                        tempRegionMap[region]    = tr.getVec<bool>("DoubleDisCo_"+region+"_1l_RPV"+myVarSuffix); 
                        tempRegionMap_CR[region] = tr.getVec<bool>("DoubleDisCo_"+region+"_NonIsoMuon_1l_RPV"+myVarSuffix); 
                    }
                }

                DoubleDisCo_passRegions.push_back(tempRegionMap); 
                DoubleDisCo_passRegions_CR.push_back(tempRegionMap_CR); 

                Jets_CR.push_back(tr.getVec<utility::LorentzVector>("Jets"         + myVarSuffix));
                GoodJets_CR.push_back(tr.getVec<bool>("NonIsoMuonJets_pt30"        + myVarSuffix));
                NGoodJets_CR.push_back(tr.getVar<int>("N" + jetsName + "Jets_pt30" + myVarSuffix));
                HT_pt30_CR.push_back(tr.getVar<double>("HT_" + htName + "_pt30"    + myVarSuffix));
                Baseline_CR.push_back(tr.getVar<bool>("pass_qcdCR"                 + myVarSuffix));
                Mbl_CR.push_back(tr.getVar<double>("Mbl"                           + myVarSuffix));
                dRbjets_CR.push_back(tr.getVar<double>("dR_bjets"                  + myVarSuffix));

                if (channel != "2")
                {
                    DoubleDisCo_massReg_CR.push_back(tr.getVar<double>("DoubleDisCo_massReg_NonIsoMuon_" + channel + "l_RPV"  + myVarSuffix));
                    DoubleDisCo_disc1_CR.push_back(tr.getVar<double>("DoubleDisCo_disc1_NonIsoMuon_" + channel + "l_RPV"      + myVarSuffix));
                    DoubleDisCo_disc2_CR.push_back(tr.getVar<double>("DoubleDisCo_disc2_NonIsoMuon_" + channel + "l_RPV"      + myVarSuffix));
                } else
                {
                    DoubleDisCo_massReg_CR.push_back(tr.getVar<double>("DoubleDisCo_massReg_NonIsoMuon_1l_RPV"    + myVarSuffix));
                    DoubleDisCo_disc1_CR.push_back(tr.getVar<double>("DoubleDisCo_disc1_NonIsoMuon_1l_RPV"        + myVarSuffix));
                    DoubleDisCo_disc2_CR.push_back(tr.getVar<double>("DoubleDisCo_disc2_NonIsoMuon_1l_RPV"        + myVarSuffix));
                }

                fwm2_top6_CR.push_back(tr.getVar<double>(mvaName + "fwm2_top6" + tag         + myVarSuffix));
                fwm3_top6_CR.push_back(tr.getVar<double>(mvaName + "fwm3_top6" + tag         + myVarSuffix));
                fwm4_top6_CR.push_back(tr.getVar<double>(mvaName + "fwm4_top6" + tag         + myVarSuffix));
                fwm5_top6_CR.push_back(tr.getVar<double>(mvaName + "fwm5_top6" + tag         + myVarSuffix));
                jmt_ev0_top6_CR.push_back(tr.getVar<double>(mvaName + "jmt_ev0_top6" + tag   + myVarSuffix));
                jmt_ev1_top6_CR.push_back(tr.getVar<double>(mvaName + "jmt_ev1_top6" + tag   + myVarSuffix));
                jmt_ev2_top6_CR.push_back(tr.getVar<double>(mvaName + "jmt_ev2_top6" + tag   + myVarSuffix));
                Jets_cm_top6_CR.push_back(tr.getVec<utility::LorentzVector>(mvaName + "Jets_cm_top6" + tag+myVarSuffix));
                nMVAJets_CR.push_back(tr.getVar<unsigned int>("nMVAJets" + tag                 + myVarSuffix));
                Stop1_pt_cm_OldSeed_CR.push_back(tr.getVar<double>("Stop1_pt_cm_OldSeed_NonIsoMuon"       + myVarSuffix));
                Stop1_eta_cm_OldSeed_CR.push_back(tr.getVar<double>("Stop1_eta_cm_OldSeed_NonIsoMuon"     + myVarSuffix));
                Stop1_phi_cm_OldSeed_CR.push_back(tr.getVar<double>("Stop1_phi_cm_OldSeed_NonIsoMuon"     + myVarSuffix));
                Stop1_mass_cm_OldSeed_CR.push_back(tr.getVar<double>("Stop1_mass_cm_OldSeed_NonIsoMuon"   + myVarSuffix));
                Stop2_pt_cm_OldSeed_CR.push_back(tr.getVar<double>("Stop2_pt_cm_OldSeed_NonIsoMuon"       + myVarSuffix));
                Stop2_eta_cm_OldSeed_CR.push_back(tr.getVar<double>("Stop2_eta_cm_OldSeed_NonIsoMuon"     + myVarSuffix));
                Stop2_phi_cm_OldSeed_CR.push_back(tr.getVar<double>("Stop2_phi_cm_OldSeed_NonIsoMuon"     + myVarSuffix));
                Stop2_mass_cm_OldSeed_CR.push_back(tr.getVar<double>("Stop2_mass_cm_OldSeed_NonIsoMuon"   + myVarSuffix));
                Stop1_pt_cm_TopSeed_CR.push_back(tr.getVar<double>("Stop1_pt_cm_TopSeed_NonIsoMuon"       + myVarSuffix));
                Stop1_eta_cm_TopSeed_CR.push_back(tr.getVar<double>("Stop1_eta_cm_TopSeed_NonIsoMuon"     + myVarSuffix));
                Stop1_phi_cm_TopSeed_CR.push_back(tr.getVar<double>("Stop1_phi_cm_TopSeed_NonIsoMuon"     + myVarSuffix));
                Stop1_mass_cm_TopSeed_CR.push_back(tr.getVar<double>("Stop1_mass_cm_TopSeed_NonIsoMuon"   + myVarSuffix));
                Stop2_pt_cm_TopSeed_CR.push_back(tr.getVar<double>("Stop2_pt_cm_TopSeed_NonIsoMuon"       + myVarSuffix));
                Stop2_eta_cm_TopSeed_CR.push_back(tr.getVar<double>("Stop2_eta_cm_TopSeed_NonIsoMuon"     + myVarSuffix));
                Stop2_phi_cm_TopSeed_CR.push_back(tr.getVar<double>("Stop2_phi_cm_TopSeed_NonIsoMuon"     + myVarSuffix));
                Stop2_mass_cm_TopSeed_CR.push_back(tr.getVar<double>("Stop2_mass_cm_TopSeed_NonIsoMuon"   + myVarSuffix));

                GoodLeptons_CR.push_back(tr.getVec<std::pair<std::string, utility::LorentzVector>>("GoodNonIsoMuons" + myVarSuffix));
                GoodLeptonsCharge_CR.push_back(tr.getVec<int>("GoodNonIsoMuonsCharge"                                + myVarSuffix));
                GoodLeptonsMiniIso_CR.push_back(tr.getVec<double>("GoodNonIsoMuonsMiniIso"                           + myVarSuffix));
               
                Stop1_mass_PtRank_matched_CR.push_back(  runtype == "Data" ? -999.0 : tr.getVar<float>("stop1_ptrank_mass" + myVarSuffix));
                Stop2_mass_PtRank_matched_CR.push_back(  runtype == "Data" ? -999.0 : tr.getVar<float>("stop2_ptrank_mass" + myVarSuffix));
                Stop1_mass_MassRank_matched_CR.push_back(runtype == "Data" ? -999.0 : tr.getVar<float>("stop1_mrank_mass"  + myVarSuffix));
                Stop2_mass_MassRank_matched_CR.push_back(runtype == "Data" ? -999.0 : tr.getVar<float>("stop2_mrank_mass"  + myVarSuffix));
                Stop_mass_average_matched_CR.push_back(  runtype == "Data" ? -999.0 : tr.getVar<double>("stop_avemass"     + myVarSuffix));

                std::vector<double> tempJets_flavb_CR;
                std::vector<double> tempJets_flavc_CR;
                std::vector<double> tempJets_flavg_CR;
                std::vector<double> tempJets_flavuds_CR;
                std::vector<double> tempJets_flavq_CR;

                for (unsigned int iJet = 1; iJet <= nMVAJets_CR.back(); iJet++) {
                    tempJets_flavb_CR.push_back(  tr.getVar<double>("Jet" + flavName + "_flavb_"  +std::to_string(iJet)+ tag +myVarSuffix));
                    tempJets_flavc_CR.push_back(  tr.getVar<double>("Jet" + flavName + "_flavc_"  +std::to_string(iJet)+ tag +myVarSuffix));
                    tempJets_flavg_CR.push_back(  tr.getVar<double>("Jet" + flavName + "_flavg_"  +std::to_string(iJet)+ tag +myVarSuffix));
                    tempJets_flavuds_CR.push_back(tr.getVar<double>("Jet" + flavName + "_flavuds_"+std::to_string(iJet)+ tag +myVarSuffix));
                    tempJets_flavq_CR.push_back(  tr.getVar<double>("Jet" + flavName + "_flavq_"  +std::to_string(iJet)+ tag +myVarSuffix));
                }

                Jets_flavb_CR.push_back(tempJets_flavb_CR);
                Jets_flavc_CR.push_back(tempJets_flavc_CR);
                Jets_flavg_CR.push_back(tempJets_flavg_CR);
                Jets_flavuds_CR.push_back(tempJets_flavuds_CR);
                Jets_flavq_CR.push_back(tempJets_flavq_CR);

                double theWeight = 1.0;
                if(runtype == "MC" )
                {
                    auto totalWeight = tr.getVar<double>("TotalWeight" + myVarSuffix);
                    theWeight *= totalWeight;
                    const auto& datasetName = tr.getVar<std::string>("dataset" + myVarSuffix);

                    // tt variations have unique samples "except" FSR/ISR
                    // FSR/ISR are event weight based variation 
                    // So, calculate the event weights for FSR/ISR to make unique samples
                    const auto& PSweight_FSRUp   = tr.getVar<double>("PSweight_FSRUp"  );
                    const auto& PSweight_FSRDown = tr.getVar<double>("PSweight_FSRDown");
                    const auto& PSweight_ISRUp   = tr.getVar<double>("PSweight_ISRUp"  );
                    const auto& PSweight_ISRDown = tr.getVar<double>("PSweight_ISRDown");
 
                    if (datasetName.find("fsrUP") != std::string::npos)
                        totalWeight *= PSweight_FSRUp;  

                    else if (datasetName.find("fsrDOWN") != std::string::npos)
                        totalWeight *= PSweight_FSRDown; 

                    else if (datasetName.find("isrUP") != std::string::npos)
                        totalWeight *= PSweight_ISRUp; 

                    else if (datasetName.find("isrDOWN") != std::string::npos)
                        totalWeight *= PSweight_ISRDown; 
                }

                weight.push_back(theWeight);
                weight_CR.push_back(theWeight);
            }

            const std::map<std::string, bool> cut_map
            {
                {"_0l"               , Baseline[0]},
                {"_1l"               , Baseline[1]},                         
                {"_2l"               , Baseline[2]}, 
                {"_QCDCR"            , Baseline_CR[0]}, 
            };

            // Let us assume that nMVAJets is the same for all channels
            // As well as definition of regions
            // Not a bad assumption at the moment...
            if(!initHistos)
            {
                Preinit(nMVAJets.front() > nMVAJets.back() ? nMVAJets.front() : nMVAJets.back());
                InitHistos(cut_map, regions);
                initHistos = true;
            }

            // Fill once per event
            my_histos["EventCounter"]->Fill(eventCounter);

            std::map<std::string, bool>               njetsMap;
            std::map<std::string, std::vector<bool> > ABCDmap;

            for(auto& kv : cut_map)
            {

                bool isQCD;
                // Extract "0", "1", "2", or "Q" from cut string e.g. _1l
                int channel = 1;
                if (kv.first.size() > 0 and kv.first.substr(2,1) == "l")
                {
                    std::string chunk = kv.first.substr(1,1);
                    
                    // One control region for the moment, so pick any of three channels
                    if (chunk == "Q")
                    {
                        channel = 1;
                        isQCD = true;
                    }
                    else
                    {
                        channel = std::stoi(chunk);
                    }
                }

                unsigned int nJets = !isQCD ? Jets_cm_top6[channel].size() : Jets_cm_top6_CR[channel].size();

                // For 0L, we always use the NGoodJets case
                njetsMap = {
                    {"Incl",                                                        true},
                    {"7",      !isQCD ? NGoodJets[channel]==7   : NGoodJets_CR[channel]==7},
                    {"8",      !isQCD ? NGoodJets[channel]==8   : NGoodJets_CR[channel]==8},
                    {"9",      !isQCD ? NGoodJets[channel]==9   : NGoodJets_CR[channel]==9},
                    {"10",     !isQCD ? NGoodJets[channel]==10  : NGoodJets_CR[channel]==10},
                    {"11",     !isQCD ? NGoodJets[channel]==11  : NGoodJets_CR[channel]==11},
                    {"12incl", !isQCD ? NGoodJets[channel]>=12  : NGoodJets_CR[channel]>=12},
                };

                ABCDmap = {};
                for (const auto& region : regions[channel])
                    ABCDmap[region] = !isQCD ? DoubleDisCo_passRegions[channel][region] : DoubleDisCo_passRegions_CR[channel][region];

                double w  = !isQCD ? weight[channel]  : weight_CR[channel];
                double ht = !isQCD ? HT_pt30[channel] : HT_pt30_CR[channel];

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

                                my_histos["h_njets"      + name]->Fill(!isQCD ?    NGoodJets[channel] : NGoodJets_CR[channel],    w);
                                my_histos["fwm2_top6"    + name]->Fill(!isQCD ?    fwm2_top6[channel] : fwm2_top6_CR[channel],    w);
                                my_histos["fwm3_top6"    + name]->Fill(!isQCD ?    fwm3_top6[channel] : fwm3_top6_CR[channel],    w);
                                my_histos["fwm4_top6"    + name]->Fill(!isQCD ?    fwm4_top6[channel] : fwm4_top6_CR[channel],    w);
                                my_histos["fwm5_top6"    + name]->Fill(!isQCD ?    fwm5_top6[channel] : fwm5_top6_CR[channel],    w);
                                my_histos["jmt_ev0_top6" + name]->Fill(!isQCD ? jmt_ev0_top6[channel] : jmt_ev0_top6_CR[channel], w);
                                my_histos["jmt_ev1_top6" + name]->Fill(!isQCD ? jmt_ev1_top6[channel] : jmt_ev1_top6_CR[channel], w);
                                my_histos["jmt_ev2_top6" + name]->Fill(!isQCD ? jmt_ev2_top6[channel] : jmt_ev2_top6_CR[channel], w);

                                my_histos["h_mbl"        + name]->Fill(!isQCD ? Mbl[channel] : Mbl_CR[channel], w);
                                my_histos["h_dRbjets"    + name]->Fill(!isQCD ? dRbjets[channel] : dRbjets_CR[channel], w);

                                // Plots of stop 4-vector are made with pt-ranked stops
                                if (Stop1_pt_cm_OldSeed  > Stop2_pt_cm_OldSeed)
                                {
                                    my_histos["Stop1_pt_cm_OldSeed"   + name]->Fill(!isQCD ? Stop1_pt_cm_OldSeed[channel]   : Stop1_pt_cm_OldSeed_CR[channel],   w);
                                    my_histos["Stop1_eta_cm_OldSeed"  + name]->Fill(!isQCD ? Stop1_eta_cm_OldSeed[channel]  : Stop1_eta_cm_OldSeed_CR[channel],  w);
                                    my_histos["Stop1_phi_cm_OldSeed"  + name]->Fill(!isQCD ? Stop1_phi_cm_OldSeed[channel]  : Stop1_phi_cm_OldSeed_CR[channel],  w);
                                    my_histos["Stop1_mass_cm_OldSeed" + name]->Fill(!isQCD ? Stop1_mass_cm_OldSeed[channel] : Stop1_mass_cm_OldSeed_CR[channel], w);

                                    my_histos["Stop2_pt_cm_OldSeed"   + name]->Fill(!isQCD ? Stop2_pt_cm_OldSeed[channel]   : Stop2_pt_cm_OldSeed_CR[channel],   w);
                                    my_histos["Stop2_eta_cm_OldSeed"  + name]->Fill(!isQCD ? Stop2_eta_cm_OldSeed[channel]  : Stop2_eta_cm_OldSeed_CR[channel],  w);
                                    my_histos["Stop2_phi_cm_OldSeed"  + name]->Fill(!isQCD ? Stop2_phi_cm_OldSeed[channel]  : Stop2_phi_cm_OldSeed_CR[channel],  w);
                                    my_histos["Stop2_mass_cm_OldSeed" + name]->Fill(!isQCD ? Stop2_mass_cm_OldSeed[channel] : Stop2_mass_cm_OldSeed_CR[channel], w);
                                } else
                                {
                                    my_histos["Stop1_pt_cm_OldSeed"   + name]->Fill(!isQCD ? Stop2_pt_cm_OldSeed[channel]   : Stop2_pt_cm_OldSeed_CR[channel],   w);
                                    my_histos["Stop1_eta_cm_OldSeed"  + name]->Fill(!isQCD ? Stop2_eta_cm_OldSeed[channel]  : Stop2_eta_cm_OldSeed_CR[channel],  w);
                                    my_histos["Stop1_phi_cm_OldSeed"  + name]->Fill(!isQCD ? Stop2_phi_cm_OldSeed[channel]  : Stop2_phi_cm_OldSeed_CR[channel],  w);
                                    my_histos["Stop1_mass_cm_OldSeed" + name]->Fill(!isQCD ? Stop2_mass_cm_OldSeed[channel] : Stop2_mass_cm_OldSeed_CR[channel], w);

                                    my_histos["Stop2_pt_cm_OldSeed"   + name]->Fill(!isQCD ? Stop1_pt_cm_OldSeed[channel]   : Stop1_pt_cm_OldSeed_CR[channel],   w);
                                    my_histos["Stop2_eta_cm_OldSeed"  + name]->Fill(!isQCD ? Stop1_eta_cm_OldSeed[channel]  : Stop1_eta_cm_OldSeed_CR[channel],  w);
                                    my_histos["Stop2_phi_cm_OldSeed"  + name]->Fill(!isQCD ? Stop1_phi_cm_OldSeed[channel]  : Stop1_phi_cm_OldSeed_CR[channel],  w);
                                    my_histos["Stop2_mass_cm_OldSeed" + name]->Fill(!isQCD ? Stop1_mass_cm_OldSeed[channel] : Stop1_mass_cm_OldSeed_CR[channel], w);
                                }

                                if (Stop1_pt_cm_TopSeed  > Stop2_pt_cm_TopSeed)
                                {
                                    my_histos["Stop1_pt_cm_TopSeed"   + name]->Fill(!isQCD ? Stop1_pt_cm_TopSeed[channel]   : Stop1_pt_cm_TopSeed_CR[channel],   w);
                                    my_histos["Stop1_eta_cm_TopSeed"  + name]->Fill(!isQCD ? Stop1_eta_cm_TopSeed[channel]  : Stop1_eta_cm_TopSeed_CR[channel],  w);
                                    my_histos["Stop1_phi_cm_TopSeed"  + name]->Fill(!isQCD ? Stop1_phi_cm_TopSeed[channel]  : Stop1_phi_cm_TopSeed_CR[channel],  w);
                                    my_histos["Stop1_mass_cm_TopSeed" + name]->Fill(!isQCD ? Stop1_mass_cm_TopSeed[channel] : Stop1_mass_cm_TopSeed_CR[channel], w);

                                    my_histos["Stop2_pt_cm_TopSeed"   + name]->Fill(!isQCD ? Stop2_pt_cm_TopSeed[channel]   : Stop2_pt_cm_TopSeed_CR[channel],   w);
                                    my_histos["Stop2_eta_cm_TopSeed"  + name]->Fill(!isQCD ? Stop2_eta_cm_TopSeed[channel]  : Stop2_eta_cm_TopSeed_CR[channel],  w);
                                    my_histos["Stop2_phi_cm_TopSeed"  + name]->Fill(!isQCD ? Stop2_phi_cm_TopSeed[channel]  : Stop2_phi_cm_TopSeed_CR[channel],  w);
                                    my_histos["Stop2_mass_cm_TopSeed" + name]->Fill(!isQCD ? Stop2_mass_cm_TopSeed[channel] : Stop2_mass_cm_TopSeed_CR[channel], w);
                                } else
                                {
                                    my_histos["Stop1_pt_cm_TopSeed"   + name]->Fill(!isQCD ? Stop2_pt_cm_TopSeed[channel]   : Stop2_pt_cm_TopSeed_CR[channel],   w);
                                    my_histos["Stop1_eta_cm_TopSeed"  + name]->Fill(!isQCD ? Stop2_eta_cm_TopSeed[channel]  : Stop2_eta_cm_TopSeed_CR[channel],  w);
                                    my_histos["Stop1_phi_cm_TopSeed"  + name]->Fill(!isQCD ? Stop2_phi_cm_TopSeed[channel]  : Stop2_phi_cm_TopSeed_CR[channel],  w);
                                    my_histos["Stop1_mass_cm_TopSeed" + name]->Fill(!isQCD ? Stop2_mass_cm_TopSeed[channel] : Stop2_mass_cm_TopSeed_CR[channel], w);

                                    my_histos["Stop2_pt_cm_TopSeed"   + name]->Fill(!isQCD ? Stop1_pt_cm_TopSeed[channel]   : Stop1_pt_cm_TopSeed_CR[channel],   w);
                                    my_histos["Stop2_eta_cm_TopSeed"  + name]->Fill(!isQCD ? Stop1_eta_cm_TopSeed[channel]  : Stop1_eta_cm_TopSeed_CR[channel],  w);
                                    my_histos["Stop2_phi_cm_TopSeed"  + name]->Fill(!isQCD ? Stop1_phi_cm_TopSeed[channel]  : Stop1_phi_cm_TopSeed_CR[channel],  w);
                                    my_histos["Stop2_mass_cm_TopSeed" + name]->Fill(!isQCD ? Stop1_mass_cm_TopSeed[channel] : Stop1_mass_cm_TopSeed_CR[channel], w);
                                }

                                my_histos["Stop1_mass_PtRank_matched"   + name]->Fill(!isQCD ? Stop1_mass_PtRank_matched[channel]   : Stop1_mass_PtRank_matched_CR[channel],   w);
                                my_histos["Stop2_mass_PtRank_matched"   + name]->Fill(!isQCD ? Stop2_mass_PtRank_matched[channel]   : Stop2_mass_PtRank_matched_CR[channel],   w);
                                my_histos["Stop1_mass_MassRank_matched" + name]->Fill(!isQCD ? Stop1_mass_MassRank_matched[channel] : Stop1_mass_MassRank_matched_CR[channel], w);
                                my_histos["Stop2_mass_MassRank_matched" + name]->Fill(!isQCD ? Stop2_mass_MassRank_matched[channel] : Stop2_mass_MassRank_matched_CR[channel], w);
                                my_histos["Stop_mass_average_matched"   + name]->Fill(!isQCD ? Stop_mass_average_matched[channel]   : Stop_mass_average_matched_CR[channel],   w);

                                auto& cmJets = Jets_cm_top6[channel];
                                auto& jetFlavb   = Jets_flavb[channel];
                                auto& jetFlavc   = Jets_flavc[channel];
                                auto& jetFlavg   = Jets_flavg[channel];
                                auto& jetFlavq   = Jets_flavq[channel];
                                auto& jetFlavuds = Jets_flavuds[channel];
                                if (isQCD)
                                {
                                    cmJets     = Jets_cm_top6_CR[channel];
                                    jetFlavb   = Jets_flavb_CR[channel];
                                    jetFlavc   = Jets_flavc_CR[channel];
                                    jetFlavg   = Jets_flavg_CR[channel];
                                    jetFlavq   = Jets_flavq_CR[channel];
                                    jetFlavuds = Jets_flavuds_CR[channel];
                                }
                                
                                for(unsigned int i = 1; i <= nJets; i++)
                                {
                                    my_histos["Jet_cm_ptrHT_" + std::to_string(i) + name]->Fill(cmJets.at(i-1).Pt()/ht, w);
                                    my_histos["Jet_cm_pt_"    + std::to_string(i) + name]->Fill(cmJets.at(i-1).Pt(),    w);
                                    my_histos["Jet_cm_eta_"   + std::to_string(i) + name]->Fill(cmJets.at(i-1).Eta(),   w);
                                    my_histos["Jet_cm_phi_"   + std::to_string(i) + name]->Fill(cmJets.at(i-1).Phi(),   w);
                                    my_histos["Jet_cm_m_"     + std::to_string(i) + name]->Fill(cmJets.at(i-1).M(),     w);
                                    my_histos["Jet_cm_E_"     + std::to_string(i) + name]->Fill(cmJets.at(i-1).E(),     w);

                                    my_histos["Jet_cm_flavb_"   + std::to_string(i) + name]->Fill(jetFlavb.at(i-1),   w);
                                    my_histos["Jet_cm_flavc_"   + std::to_string(i) + name]->Fill(jetFlavc.at(i-1),   w);
                                    my_histos["Jet_cm_flavg_"   + std::to_string(i) + name]->Fill(jetFlavg.at(i-1),   w);
                                    my_histos["Jet_cm_flavq_"   + std::to_string(i) + name]->Fill(jetFlavq.at(i-1),   w);
                                    my_histos["Jet_cm_flavuds_" + std::to_string(i) + name]->Fill(jetFlavuds.at(i-1), w);
                                }

                                my_histos["h_ht" + name]->Fill(ht, w);

                                auto& Leptons  = GoodLeptons[channel];
                                auto& LeptonsC = GoodLeptonsCharge[channel];
                                auto& LeptonsI = GoodLeptonsMiniIso[channel];
                                if (isQCD)
                                {
                                   Leptons  = GoodLeptons_CR[channel]; 
                                   LeptonsC = GoodLeptonsCharge_CR[channel]; 
                                   LeptonsI = GoodLeptonsMiniIso_CR[channel]; 
                                }

                                for(std::size_t j = 0 ; j < Leptons.size(); ++j)
                                {
                                    std::string type            = Leptons[j].first;
                                    utility::LorentzVector lvec = Leptons[j].second; 
                                    double charge               = LeptonsC[j];
                                    double iso                  = LeptonsI[j];

                                    my_histos["h_lPt"      + name]->Fill(lvec.Pt(),  w);
                                    my_histos["h_lPhi"     + name]->Fill(lvec.Phi(), w);
                                    my_histos["h_lEta"     + name]->Fill(lvec.Eta(), w);
                                    my_histos["h_lCharge"  + name]->Fill(charge ,    w);
                                    my_histos["h_lMiniIso" + name]->Fill(iso ,       w);

                                    if(type == 'e'){
                                        my_histos["h_ePt"      + name]->Fill(lvec.Pt() , w);
                                        my_histos["h_ePhi"     + name]->Fill(lvec.Phi(), w);
                                        my_histos["h_eEta"     + name]->Fill(lvec.Eta(), w);
                                        my_histos["h_eCharge"  + name]->Fill(charge,     w);
                                        my_histos["h_eMiniIso" + name]->Fill(iso,        w);
                                    } else if (type == 'm' or type == 'n') {
                                        my_histos["h_mPt"      + name]->Fill(lvec.Pt(),  w);
                                        my_histos["h_mPhi"     + name]->Fill(lvec.Phi(), w);
                                        my_histos["h_mEta"     + name]->Fill(lvec.Eta(), w);
                                        my_histos["h_mCharge"  + name]->Fill(charge,     w);
                                        my_histos["h_mMiniIso" + name]->Fill(iso,        w);
                                    }
                                }
                                auto& theJets     = Jets[channel];
                                auto& theGoodJets = GoodJets[channel];
                                if (isQCD)
                                {
                                    theJets     = Jets_CR[channel];
                                    theGoodJets = GoodJets_CR[channel];
                                }

                                for(unsigned int j = 0; j < theJets.size(); j++) 
                                {
                                    if(!theGoodJets[j]) continue;
                                    my_histos["h_jPt"  + name]->Fill(theJets.at(j).Pt(),  w); 
                                    my_histos["h_jEta" + name]->Fill(theJets.at(j).Eta(), w); 
                                    my_histos["h_jPhi" + name]->Fill(theJets.at(j).Phi(), w); 
                                    my_histos["h_jM"   + name]->Fill(theJets.at(j).M(),   w); 
                                }
                            }

                            auto& cmJets = Jets_cm_top6[channel];
                            if (isQCD)
                                cmJets = Jets_cm_top6_CR[channel];

                            for(unsigned int i = 1; i <= nJets; i++)
                            {
                                double pt  = static_cast<double>(cmJets.at(i-1).Pt());
                                my_2d_histos["h_cm_pt_jetRank"    + name]->Fill(pt,    i, w);
                                my_2d_histos["h_cm_ptrHT_jetRank" + name]->Fill(pt/ht, i, w);
                            }

                            if (region != "Incl" and njets != "Incl")
                                my_histos["h_DoubleDisCo_massReg" + name]->Fill(!isQCD ? DoubleDisCo_massReg[channel] : DoubleDisCo_massReg_CR[channel], w);

                            auto& disc1 = DoubleDisCo_disc1[channel];
                            auto& disc2 = DoubleDisCo_disc2[channel];
                            if (isQCD)
                            {
                                disc1 = DoubleDisCo_disc1_CR[channel];
                                disc2 = DoubleDisCo_disc2_CR[channel];
                            }

                            // if plotting disco, no need to make plots when cutting on
                            if (njets != "Incl")
                            {
                                my_histos["h_DoubleDisCo_disc1" + name]->Fill(disc1, w);
                                my_histos["h_DoubleDisCo_disc2" + name]->Fill(disc2, w);
                                my_2d_histos["h_DoubleDisCo_disc1_disc2" + name]->Fill(disc1, disc2, w);
                            }
                        }

                        for (unsigned int iSubRegion = 0; iSubRegion < inRegionBin.size(); iSubRegion++)
                        {
                            if (kv.second and inNjetsBin and inRegion and inRegionBin[iSubRegion])
                            {
                                name = kv.first + njetsStr + regionStr;
                                int shift = inRegionBin[iSubRegion] ? iSubRegion : 100;

                                auto& theGoodJets = NGoodJets[channel];
                                auto& theMVAjets  = nMVAJets[channel]; 
                                if (isQCD)
                                {
                                    theGoodJets = NGoodJets_CR[channel];
                                    theMVAjets  = nMVAJets_CR[channel];
                                }

                                // if plotting Njets, don't care about individual Njet cuts
                                if (njets == "Incl" and region != "Incl")
                                    my_histos["h_njets_12incl" + name]->Fill(theGoodJets>=12 ? 12-theMVAjets+shift*6 : theGoodJets-theMVAjets+shift*6, w);

                                if (njets != "Incl")
                                {
                                    auto& disc1 = DoubleDisCo_disc1[channel];
                                    auto& disc2 = DoubleDisCo_disc2[channel];
                                    if (isQCD)
                                    {
                                        disc1 = DoubleDisCo_disc1_CR[channel];
                                        disc2 = DoubleDisCo_disc2_CR[channel];
                                    }

                                    name = kv.first + njetsStr + regionStr + "_" + subRegionsMap[region][iSubRegion];
                                    my_2d_histos["h_DoubleDisCo_disc1_disc2" + name]->Fill(disc1, disc2, w);

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
