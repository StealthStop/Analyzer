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
        {"h_DoubleDisCo_massReg",       180,    0, 1800},
        {"fwm2_top6",                    50,    0,    1},
        {"fwm3_top6",                    50,    0,    1},
        {"fwm4_top6",                    50,    0,    1},
        {"fwm5_top6",                    50,    0,    1},
        {"jmt_ev0_top6",                 50,    0,    1},
        {"jmt_ev1_top6",                 50,    0,    1},
        {"jmt_ev2_top6",                 50,    0,    1},
        {"Stop1_pt_cm_OldSeed",         180,    0, 1800},
        {"Stop1_eta_cm_OldSeed",         80,   -6,    6},
        {"Stop1_phi_cm_OldSeed",         64,   -4,    4},
        {"Stop1_mass_cm_OldSeed",       180,    0, 1800},
        {"Stop2_pt_cm_OldSeed",         180,    0, 1800},
        {"Stop2_eta_cm_OldSeed",         80,   -6,    6},
        {"Stop2_phi_cm_OldSeed",         64,   -4,    4},
        {"Stop2_mass_cm_OldSeed",       180,    0, 1800},
        {"Stop1_pt_cm_TopSeed",         180,    0, 1800},
        {"Stop1_eta_cm_TopSeed",         80,   -6,    6},
        {"Stop1_phi_cm_TopSeed",         64,   -4,    4},
        {"Stop1_mass_cm_TopSeed",       180,    0, 1800},
        {"Stop2_pt_cm_TopSeed",         180,    0, 1800},
        {"Stop2_eta_cm_TopSeed",         80,   -6,    6},
        {"Stop2_phi_cm_TopSeed",         64,   -4,    4},
        {"Stop2_mass_cm_TopSeed",       180,    0, 1800},
        {"h_njets",                      21, -0.5, 20.5},
        {"h_nbjets",                     21, -0.5, 20.5},
        {"h_ntops",                       7, -0.5,  6.5},
        {"h_njets_12incl",               24, -0.5, 23.5},
        {"h_ht",                        720,    0, 7200},
        {"h_dRbjets",                   180,    0,    6},
        {"h_mbl",                       180,    0,  360},
        {"Stop1_mass_PtRank_matched",   180,    0, 1800},
        {"Stop2_mass_PtRank_matched",   180,    0, 1800},
        {"Stop1_mass_MassRank_matched", 180,    0, 1800},
        {"Stop2_mass_MassRank_matched", 180,    0, 1800},
        {"Stop_mass_average_matched",   180,    0, 1800},
    };

    hist2DInfos = {
        {"h_DoubleDisCo_disc1_disc2", 100,    0,    1, 100,    0,   1}, 
        {"h_cm_pt_jetRank",           150,    0, 1500,  10,    0,  10},
        {"h_cm_ptrHT_jetRank",        150,    0,    1,  10,    0,  10},
        {"h_nRtops_vs_nMtops",          7, -0.5,  6.5,   7, -0.5, 6.5},
    };

    channels      = {"0", "1", "2"};
    njets         = {"Incl", "6", "7", "8", "9", "10", "11", "12incl"};
    ttvars        = {"nom", "fsrUp", "fsrDown", "isrUp", "isrDown"   };
    my_var_suffix = {""                                              };

    my_histos.emplace("EventCounter", std::make_shared<TH1D>("EventCounter", "EventCounter", 2, -1.1, 1.1) );
}

void AnalyzeDoubleDisCo::Debug(const std::string& message, int line = 0)
{
    if (debug)
        std::cout << "L" << line << ": " << message << std::endl;
}

void AnalyzeDoubleDisCo::makeSubregions(const std::vector<std::vector<std::string>>& regionVec)
{

    Debug("Creating map of subregions and regions", __LINE__);

    // ------------------------------------------------------------------------------------
    // Translate the generic "A", "B", "C", "D" names into unique names based on the region
    // In each SubDivision of BD, CD, D : subregions are the same notation like A,B,C,D
    //  -- E.g., in SubDiv-BD: "B" : A' / "D" : C' / "E" : B' / "F" : D' 
    //  -- Naming conventions agreed upon with Validation code
    // ------------------------------------------------------------------------------------
    for (auto& regions : regionVec)
    {
        for (auto& region : regions)
        { 
            if (subRegionsMap.find(region) != subRegionsMap.end())
                continue;

            // SubDiv-BD: subregions (B,D,E,F) in BD regions
            if (region.find("bd") != std::string::npos or region.find("_BD") != std::string::npos) {
                subRegionsMap[region].push_back("B");
                subRegionsMap[region].push_back("E");
                subRegionsMap[region].push_back("D");
                subRegionsMap[region].push_back("F");
            // SubDiv-CD: subregions (C,D,G,H) in CD regions
            } else if (region.find("cd") != std::string::npos or region.find("_CD") != std::string::npos) {
                subRegionsMap[region].push_back("C");
                subRegionsMap[region].push_back("D");
                subRegionsMap[region].push_back("G");
                subRegionsMap[region].push_back("H");
            // SubDiv-D: subregions (dA, dB, dC, dD) in D region
            } else if (region.find("subDiv") != std::string::npos or region.find("_D") != std::string::npos) {
                subRegionsMap[region].push_back("dA");
                subRegionsMap[region].push_back("dB");
                subRegionsMap[region].push_back("dC");
                subRegionsMap[region].push_back("dD");
            // ABCD region
            } else {
                subRegionsMap[region].push_back("A");
                subRegionsMap[region].push_back("B");
                subRegionsMap[region].push_back("C");
                subRegionsMap[region].push_back("D");
            }
        }
    }
}

void AnalyzeDoubleDisCo::Preinit(unsigned int nNNJets, unsigned int nLeptons)
{

    Debug("Doing pre-initialization to define jet and lepton histograms", __LINE__);
    for(unsigned int i = 1; i <= nNNJets ; i++)
    {
        histInfos.push_back({"Jet_cm_pt_"      + std::to_string(i), 180,  0,  1800});
        histInfos.push_back({"Jet_cm_ptrHT_"   + std::to_string(i), 180,  0,     1});
        histInfos.push_back({"Jet_cm_eta_"     + std::to_string(i),  80, -6,     6});
        histInfos.push_back({"Jet_cm_phi_"     + std::to_string(i),  64, -4,     4});
        histInfos.push_back({"Jet_cm_m_"       + std::to_string(i), 180,  0,   360});
        histInfos.push_back({"Jet_cm_E_"       + std::to_string(i), 180,  0,  1800});

        histInfos.push_back({"Jet_cm_flavb_"   + std::to_string(i),  80,  0,     1});
        histInfos.push_back({"Jet_cm_flavc_"   + std::to_string(i),  80,  0,     1});
        histInfos.push_back({"Jet_cm_flavg_"   + std::to_string(i),  80,  0,     1});
        histInfos.push_back({"Jet_cm_flavq_"   + std::to_string(i),  80,  0,     1});
        histInfos.push_back({"Jet_cm_flavuds_" + std::to_string(i),  80,  0,     1});

        histInfos.push_back({"h_jPt_"          + std::to_string(i), 180,  0,  1800});
        histInfos.push_back({"h_jPtrHT_"       + std::to_string(i), 180,  0,  1800});
        histInfos.push_back({"h_jEta_"         + std::to_string(i),  80, -6,     6});
        histInfos.push_back({"h_jPhi_"         + std::to_string(i),  64, -4,     4});
        histInfos.push_back({"h_jM_"           + std::to_string(i), 180,  0,   360});
        histInfos.push_back({"h_jE_"           + std::to_string(i), 180,  0,  1800});

        histInfos.push_back({"h_jFlavb_"       + std::to_string(i),  80,  0,     1});
        histInfos.push_back({"h_jFlavc_"       + std::to_string(i),  80,  0,     1});
        histInfos.push_back({"h_jFlavg_"       + std::to_string(i),  80,  0,     1});
        histInfos.push_back({"h_jFlavq_"       + std::to_string(i),  80,  0,     1});
        histInfos.push_back({"h_jFlavuds_"     + std::to_string(i),  80,  0,     1});
    }

    // Make one-less-jet combined jet histogram
    histInfos.push_back({"h_combined" + std::to_string(nNNJets-1) + "thJet_pt",    180,  0, 1800});
    histInfos.push_back({"h_combined" + std::to_string(nNNJets-1) + "thJet_ptrHT", 180,  0,    1});
    histInfos.push_back({"h_combined" + std::to_string(nNNJets-1) + "thJet_eta",    80, -6,    6});
    histInfos.push_back({"h_combined" + std::to_string(nNNJets-1) + "thJet_phi",    64, -4,    4});
    histInfos.push_back({"h_combined" + std::to_string(nNNJets-1) + "thJet_m",     180,  0,  360});
    histInfos.push_back({"h_combined" + std::to_string(nNNJets-1) + "thJet_E",     180,  0, 1800});

    histInfos.push_back({"h_combined" + std::to_string(nNNJets) + "thJet_pt",    180,  0, 1800});
    histInfos.push_back({"h_combined" + std::to_string(nNNJets) + "thJet_ptrHT", 180,  0,    1});
    histInfos.push_back({"h_combined" + std::to_string(nNNJets) + "thJet_eta",    80, -6,    6});
    histInfos.push_back({"h_combined" + std::to_string(nNNJets) + "thJet_phi",    64, -4,    4});
    histInfos.push_back({"h_combined" + std::to_string(nNNJets) + "thJet_m",     180,  0,  360});
    histInfos.push_back({"h_combined" + std::to_string(nNNJets) + "thJet_E",     180,  0, 1800});

    for(unsigned int i = 1; i <= nLeptons; i++)
    {
        histInfos.push_back({"h_lPt_"      + std::to_string(i),  180,    0, 1800});
        histInfos.push_back({"h_lEta_"     + std::to_string(i),   80,   -6,    6});
        histInfos.push_back({"h_lPhi_"     + std::to_string(i),   64,   -4,    4});
        histInfos.push_back({"h_lCharge_"  + std::to_string(i),    2,   -1,    1});
        histInfos.push_back({"h_lMiniIso_" + std::to_string(i), 2880,    0,    5});

        histInfos.push_back({"h_ePt_"      + std::to_string(i),  180,    0, 1800});
        histInfos.push_back({"h_eEta_"     + std::to_string(i),   80,   -6,    6});
        histInfos.push_back({"h_ePhi_"     + std::to_string(i),   64,   -4,    4});
        histInfos.push_back({"h_eCharge_"  + std::to_string(i),    2,   -1,    1});
        histInfos.push_back({"h_eMiniIso_" + std::to_string(i), 2880,    0,    5});

        histInfos.push_back({"h_mPt_"      + std::to_string(i),  180,    0, 1800});
        histInfos.push_back({"h_mPhi_"     + std::to_string(i),   80,   -4,    4});
        histInfos.push_back({"h_mEta_"     + std::to_string(i),   64,   -6,    6});
        histInfos.push_back({"h_mCharge_"  + std::to_string(i),    2,   -1,    1});
        histInfos.push_back({"h_mMiniIso_" + std::to_string(i), 2880,    0,    5});
    }
}

void AnalyzeDoubleDisCo::InitHistos(const std::map<std::string, bool>& cutMap, const std::vector<std::vector<std::string>>& regionsVec)
{

    Debug("Initializing all histograms", __LINE__);
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    // Generates a map of region to constituent subregions
    makeSubregions(regionsVec);

    // --------------------------------------------------------------
    // loop over the cutmap
    // mycut.first : cut string, mycut.second : boolean if cut passed
    // --------------------------------------------------------------
    for (auto& mycut : cutMap)
    {
        // ----------------------
        // loop over 1D hist info
        // ----------------------
        for (const auto& hInfo : histInfos)
        { 
            // -------------------------------
            // loop over njets
            // Njet string, can also be "Incl"
            // -------------------------------
            for (const auto& Njet : njets)
            {
                std::string njetStr = "";
                // For 1D histos, don't make any where we exclude all but one njets bin
                if (Njet != "Incl")
                    njetStr = "_Njets" + Njet;

                // --------------------------------------------------------------------------            
                // loop over the tt and ttvars fsrUp/Down, isrUp/Down
                // to make indivudual histograms with the label tt and ttvars in TT root file
                // --------------------------------------------------------------------------
                for (const auto& ttvar : ttvars)
                {
                    std::string ttvarStr = "";
                    if (ttvar != "nom")
                        ttvarStr = "_" + ttvar; 
                
                    // ------------------------------------------------------
                    // loop over the regions    
                    // regions : a vector of region string names for 0L or 1L
                    // ------------------------------------------------------
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

                        if (hInfo.name.find("njets") == std::string::npos and region != "ABCD")
                            continue;

                        if (mycut.first.find("QCDCR_") != std::string::npos and ttvar != "nom")
                            continue; 

                        std::string name = hInfo.name + mycut.first + njetStr + regionStr + ttvarStr;
                        my_histos.emplace(name, std::make_shared<TH1D>((name).c_str(),(name).c_str(), hInfo.nBins, hInfo.low, hInfo.high));

                        Debug("Initializing histogram with name: " + name, __LINE__);
                    } 
                } 
            } 
        } 

        // ----------------------
        // loop over 2D hist info
        // ----------------------
        for(const auto& h2dInfo : hist2DInfos)
        {
            // ---------------
            // loop over njets
            // ---------------
            for (const auto& Njet : njets)
            {
                // For the disc1 vs disc2 plots, no need for njets inclusive one
                if (Njet == "Incl" && h2dInfo.name.find("DisCo") != std::string::npos)
                    continue;

                std::string njetStr = "";
                // For 2D histos, don't make any where we exclude all but one njets bin
                if (Njet != "Incl")
                    njetStr = "_Njets" + Njet;

                // --------------------------------------------------------------------------            
                // loop over the tt and ttvars fsrUp/Down, isrUp/Down
                // to make indivudual histograms with the label tt and ttvars in TT root file
                // --------------------------------------------------------------------------
                for (const auto& ttvar : ttvars)
                {   
                    std::string ttvarStr = "";
                    if (ttvar != "nom")
                        ttvarStr = "_" + ttvar;

                    // ---------------------
                    // loop over the regions
                    // ---------------------
                    for (const auto& regionPair : subRegionsMap)
                    {
                        std::string region = regionPair.first;
                        std::vector<std::string> subregions = regionPair.second;

                        std::string regionStr = "";
                        if (region != "Incl")
                            regionStr = "_" + region;

                        if (region != "ABCD")
                            continue;

                        if (ttvar != "nom" and h2dInfo.name.find("DoubleDisCo") == std::string::npos)
                            continue;

                        if (mycut.first.find("QCDCR_") != std::string::npos and ttvar != "nom")
                            continue; 

                        std::string name = h2dInfo.name + mycut.first + njetStr + regionStr + ttvarStr;
                        my_2d_histos.emplace(name, std::make_shared<TH2D>((name).c_str(),(name).c_str(), h2dInfo.nBinsX, h2dInfo.lowX, h2dInfo.highX, h2dInfo.nBinsY, h2dInfo.lowY, h2dInfo.highY));

                        Debug("Initializing histogram with name: " + name, __LINE__);
                    } 
                } 
            } 
        } 
    } 
}

void AnalyzeDoubleDisCo::Loop(NTupleReader& tr, double, int maxevents, bool isQuiet)
{
    debug = !isQuiet;

    Debug("Entering the Loop function", __LINE__);

    const auto& dataset                           = tr.getVar<std::string>("dataset"                          ); // to check the collection name for JEC/JER 
    const auto& filetag                           = tr.getVar<std::string>("filetag"                          );
    const auto& runtype                           = tr.getVar<std::string>("runtype"                          );    
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
    const auto& DoubleDisCo_Cfg_2l_RPV            = tr.getVar<std::string>("DoubleDisCo_Cfg_2l_RPV"           );  
    const auto& DoubleDisCo_Model_2l_RPV          = tr.getVar<std::string>("DoubleDisCo_Model_2l_RPV"         );    
    const auto& DoubleDisCo_Cfg_NonIsoMuon_2l_RPV = tr.getVar<std::string>("DoubleDisCo_Cfg_NonIsoMuon_2l_RPV");

    // -------------------------------------
    // check the collection name for JEC/JER
    // ------------------------------------- 
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

        Jet                 jet(myVarSuffix);
        BJet                bjet(myVarSuffix);
        Muon                muon(myVarSuffix);
        Photon              photon(myVarSuffix);
        Baseline            baseline(myVarSuffix);
        Electron            electron(myVarSuffix);
        StopJets            stopJets(myVarSuffix);
        RunTopTagger        topTagger(TopTaggerCfg, myVarSuffix);
        StopGenMatch        stopGenMatch(myVarSuffix);
        CommonVariables     commonVariables(myVarSuffix);
        MakeMVAVariables    makeMVAVariables(                false,  myVarSuffix,        "GoodJets_pt30",       false, true, 7, 2, ""                                  );
        MakeMVAVariables    makeMVAVariables_NonIsoMuon(     false,  myVarSuffix,        "NonIsoMuonJets_pt30", false, true, 7, 2, ""                                  );
        // 0l
        // note that if we make the inputs to the NN, we use just the GoodJets_pt30 collection to derive things
        // but, if we define the QCD CR selection, we use the NonIsoMuonJets_pt30 collection 
        DeepEventShape      neuralNetwork0L(           DoubleDisCo_Cfg_0l_RPV,            DoubleDisCo_Model_0l_RPV, "Info", true, myVarSuffix                          );
        DeepEventShape      neuralNetwork0L_NonIsoMuon(DoubleDisCo_Cfg_NonIsoMuon_0l_RPV, DoubleDisCo_Model_0l_RPV, "Info", true, myVarSuffix                          );
        MakeStopHemispheres stopHemispheres_TopSeed(           "StopJets", "GoodStopJets", "NGoodStopJets", "_TopSeed",            myVarSuffix, Hemisphere::TopSeed    );
        MakeStopHemispheres stopHemispheres_TopSeed_NonIsoMuon("StopJets", "GoodStopJets", "NGoodStopJets", "_TopSeed_NonIsoMuon", myVarSuffix, Hemisphere::InvMassSeed);
        // 1l
        DeepEventShape      neuralNetwork1L(           DoubleDisCo_Cfg_1l_RPV,            DoubleDisCo_Model_1l_RPV, "Info", true, myVarSuffix                                    ); 
        DeepEventShape      neuralNetwork1L_NonIsoMuon(DoubleDisCo_Cfg_NonIsoMuon_1l_RPV, DoubleDisCo_Model_1l_RPV, "Info", true, myVarSuffix                                    );
        MakeStopHemispheres stopHemispheres_OldSeed(           "Jets", "GoodJets_pt20",       "NGoodJets_pt20",       "_OldSeed",            myVarSuffix, Hemisphere::InvMassSeed);
        MakeStopHemispheres stopHemispheres_OldSeed_NonIsoMuon("Jets", "NonIsoMuonJets_pt20", "NNonIsoMuonJets_pt30", "_OldSeed_NonIsoMuon", myVarSuffix, Hemisphere::InvMassSeed);        
        // 2l
        DeepEventShape      neuralNetwork2L(           DoubleDisCo_Cfg_2l_RPV,            DoubleDisCo_Model_2l_RPV, "Info", true, myVarSuffix); 
        DeepEventShape      neuralNetwork2L_NonIsoMuon(DoubleDisCo_Cfg_NonIsoMuon_2l_RPV, DoubleDisCo_Model_2l_RPV, "Info", true, myVarSuffix);
        
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
        tr.registerFunction(makeMVAVariables);
        tr.registerFunction(makeMVAVariables_NonIsoMuon);
        tr.registerFunction(stopJets);
        tr.registerFunction(stopHemispheres_TopSeed);
        tr.registerFunction(stopHemispheres_OldSeed);
        tr.registerFunction(stopHemispheres_TopSeed_NonIsoMuon);
        tr.registerFunction(stopHemispheres_OldSeed_NonIsoMuon);

        if (runtype == "MC")
        {
            ScaleFactors        scaleFactors(runYear, leptonFileName, hadronicFileName, meanFileName, myVarSuffix);
            BTagCorrector       bTagCorrector(bjetFileName, "", bjetCSVFileName, "", filetag);
            bTagCorrector.SetVarNames("GenParticles_PdgId", "Jets"+myVarSuffix, "GoodJets_pt30"+myVarSuffix, "Jets"+myVarSuffix+"_bJetTagDeepCSVtotb", "Jets"+myVarSuffix+"_partonFlavor", myVarSuffix);

            tr.registerFunction(bTagCorrector);
            tr.registerFunction(scaleFactors);
            tr.registerFunction(stopGenMatch);
       }

        tr.registerFunction(neuralNetwork0L);
        tr.registerFunction(neuralNetwork0L_NonIsoMuon);
        tr.registerFunction(neuralNetwork1L);
        tr.registerFunction(neuralNetwork1L_NonIsoMuon);
        tr.registerFunction(neuralNetwork2L);
        tr.registerFunction(neuralNetwork2L_NonIsoMuon);
    }

    Debug("Initialized modules to run", __LINE__);

    // ---------------
    // main event loop
    // ---------------
    while( tr.getNextEvent() )
    {
        if (maxevents != -1 && tr.getEvtNum() > maxevents)
            break;        

        if (tr.getEvtNum() % 1000 == 0)
            printf("  Event %i\n", tr.getEvtNum() );

        // Fill once per event
        const auto& eventCounter = tr.getVar<int>("eventCounter"); 
        my_histos["EventCounter"]->Fill(eventCounter);

        Debug("Initializing variables for signal and control regions", __LINE__);

        for(const auto& myVarSuffix : my_var_suffix)
        {
            // If an event is not interesting for any channel selection (see explanation and definition
            // of lostCauseEvent boolean in Baseline.h), then move on here and do not waste any more 
            // time on looping over histos or getting vars.
            // This only matters if the user specifies -s on the command line
            // N.B. We already filled the EventCounter histogram for our due-diligence of counting up
            // every single event.
            const auto& lostCauseEvent = tr.getVar<bool>("lostCauseEvent" + myVarSuffix);
            const auto& fastMode       = tr.getVar<bool>("fastMode");

            if (lostCauseEvent and fastMode)
                break;

            // Put 0L and 1L version of variables into vector
            // 0th position is for 0L and 1st position is for 1L for convenience in the event loop
            // For 0L things in the CR, nominal version of jets and derived quantities are used, take note !
            std::vector<std::vector<utility::LorentzVector> >       Jets_CR                        ;
            std::vector<std::vector<utility::LorentzVector> >       Jets_cm_top6_CR                ;

            std::vector<bool>                                       Baseline_CR                    ;
            std::vector<double>                                     weight_CR                      ;
            std::vector<double>                                     weight_CR_fsrUp                ;
            std::vector<double>                                     weight_CR_fsrDown              ;
            std::vector<double>                                     weight_CR_isrUp                ;
            std::vector<double>                                     weight_CR_isrDown              ;
            std::vector<std::vector<bool> >                         GoodJets_CR                    ;
            std::vector<int>                                        NGoodJets_CR                   ;
            std::vector<int>                                        NGoodBJets_CR                  ;
            std::vector<int>                                        NTops_CR                       ;
            std::vector<int>                                        NRTops_CR                      ;
            std::vector<int>                                        NMTops_CR                      ;
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
            std::vector<unsigned int>                               nMVALeptons_CR                 ;
            std::vector<std::vector<double> >                       Jets_cm_flavb_CR               ;
            std::vector<std::vector<double> >                       Jets_cm_flavc_CR               ;
            std::vector<std::vector<double> >                       Jets_cm_flavg_CR               ;
            std::vector<std::vector<double> >                       Jets_cm_flavuds_CR             ;
            std::vector<std::vector<double> >                       Jets_cm_flavq_CR               ;

            std::vector<std::vector<float> >                        Jets_flavb_CR                  ;
            std::vector<std::vector<float> >                        Jets_flavc_CR                  ;
            std::vector<std::vector<float> >                        Jets_flavg_CR                  ;
            std::vector<std::vector<float> >                        Jets_flavuds_CR                ;
            std::vector<std::vector<float> >                        Jets_flavq_CR                  ;

            std::vector<std::vector<std::string> >                  regions_CR                     ;
            std::vector<double>                                     combinedNthJetPt_CR            ;
            std::vector<double>                                     combinedNthJetPtrHT_CR         ;
            std::vector<double>                                     combinedNthJetEta_CR           ;
            std::vector<double>                                     combinedNthJetPhi_CR           ;
            std::vector<double>                                     combinedNthJetM_CR             ;
            std::vector<double>                                     combinedNthJetE_CR             ;

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

            // Setup all same set of variables for the signal region selections i.e. no _CR
            std::vector<std::vector<utility::LorentzVector> >       Jets                           ;
            std::vector<std::vector<utility::LorentzVector> >       Jets_cm_top6                   ;

            std::vector<bool>                                       Baseline                       ;
            std::vector<bool>                                       Baseline_blind                 ;
            std::vector<double>                                     weight                         ;
            std::vector<double>                                     weight_fsrUp                   ;
            std::vector<double>                                     weight_fsrDown                 ;
            std::vector<double>                                     weight_isrUp                   ;
            std::vector<double>                                     weight_isrDown                 ;            
            std::vector<std::vector<bool> >                         GoodJets                       ;
            std::vector<int>                                        NGoodJets                      ;
            std::vector<int>                                        NGoodBJets                     ;
            std::vector<int>                                        NTops                          ;
            std::vector<int>                                        NRTops                         ;
            std::vector<int>                                        NMTops                         ;
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
            std::vector<unsigned int>                               nMVALeptons                    ;
            std::vector<std::vector<double> >                       Jets_cm_flavb                  ;
            std::vector<std::vector<double> >                       Jets_cm_flavc                  ;
            std::vector<std::vector<double> >                       Jets_cm_flavg                  ;
            std::vector<std::vector<double> >                       Jets_cm_flavuds                ;
            std::vector<std::vector<double> >                       Jets_cm_flavq                  ;

            std::vector<std::vector<float> >                        Jets_flavb                     ;
            std::vector<std::vector<float> >                        Jets_flavc                     ;
            std::vector<std::vector<float> >                        Jets_flavg                     ;
            std::vector<std::vector<float> >                        Jets_flavuds                   ;
            std::vector<std::vector<float> >                        Jets_flavq                     ;

            std::vector<std::vector<std::string> >                  regions                        ;
            std::vector<double>                                     combinedNthJetPt               ;
            std::vector<double>                                     combinedNthJetPtrHT            ;
            std::vector<double>                                     combinedNthJetEta              ;
            std::vector<double>                                     combinedNthJetPhi              ;
            std::vector<double>                                     combinedNthJetM                ;
            std::vector<double>                                     combinedNthJetE                ;

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

            Debug("Filling vector of variables for signal and control regions", __LINE__);

            // ------------------
            // loop over channels
            // ------------------
            for (auto& channel : channels)
            {
                // Collection names for 0L, 1L, 2L do  not change when in the QCDCR, because of the jet collection used...
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

                Jets.push_back(tr.getVec<utility::LorentzVector>("Jets"                       + myVarSuffix));
                Jets_cm_top6.push_back(tr.getVec<utility::LorentzVector>("Jets_cm_top6" + myVarSuffix));

                GoodJets.push_back(tr.getVec<bool>("GoodJets_pt30"                                 + myVarSuffix));
                NGoodJets.push_back(tr.getVar<int>("NGoodJets_pt30"                                + myVarSuffix));
                NGoodBJets.push_back(tr.getVar<int>("NGoodBJets_pt30"                              + myVarSuffix));
                NTops.push_back(tr.getVar<int>("ntops"                                             + myVarSuffix));
                NRTops.push_back(tr.getVar<int>("ntops_3jet"                                       + myVarSuffix));
                NMTops.push_back(tr.getVar<int>("ntops_1jet"                                       + myVarSuffix));
                HT_pt30.push_back(tr.getVar<double>("HT_trigger_pt30"                              + myVarSuffix));
                Baseline.push_back(tr.getVar<bool>("passBaseline" + channel + "l_Good"             + myVarSuffix));
                Baseline_blind.push_back(tr.getVar<bool>("passBaseline" + channel + "l_Good_blind" + myVarSuffix));
                Mbl.push_back(tr.getVar<double>("Mbl"                                              + myVarSuffix));
                dRbjets.push_back(tr.getVar<double>("dR_bjets"                                     + myVarSuffix));

                DoubleDisCo_massReg.push_back(tr.getVar<double>("DoubleDisCo_massReg_" + channel + "l_RPV" + myVarSuffix));
                DoubleDisCo_disc1.push_back(tr.getVar<double>("DoubleDisCo_disc1_"     + channel + "l_RPV" + myVarSuffix));
                DoubleDisCo_disc2.push_back(tr.getVar<double>("DoubleDisCo_disc2_"     + channel + "l_RPV" + myVarSuffix));

                fwm2_top6.push_back(tr.getVar<double>("fwm2_top6"          + myVarSuffix));
                fwm3_top6.push_back(tr.getVar<double>("fwm3_top6"          + myVarSuffix));
                fwm4_top6.push_back(tr.getVar<double>("fwm4_top6"          + myVarSuffix));
                fwm5_top6.push_back(tr.getVar<double>("fwm5_top6"          + myVarSuffix));
                jmt_ev0_top6.push_back(tr.getVar<double>("jmt_ev0_top6"    + myVarSuffix));
                jmt_ev1_top6.push_back(tr.getVar<double>("jmt_ev1_top6"    + myVarSuffix));
                jmt_ev2_top6.push_back(tr.getVar<double>("jmt_ev2_top6"    + myVarSuffix));
                nMVAJets.push_back(tr.getVar<unsigned int>("nMVAJets"      + myVarSuffix));
                nMVALeptons.push_back(tr.getVar<unsigned int>("nMVALeptons"      + myVarSuffix));

                combinedNthJetPt.push_back(tr.getVar<double>("combined"    + std::to_string(nMVAJets.back()) + "thToLastJet" + flavName + "_pt_cm" + myVarSuffix));
                combinedNthJetPtrHT.push_back(tr.getVar<double>("combined" + std::to_string(nMVAJets.back()) + "thToLastJet" + flavName + "_ptrHT_cm" + myVarSuffix));
                combinedNthJetEta.push_back(tr.getVar<double>("combined"   + std::to_string(nMVAJets.back()) + "thToLastJet" + flavName + "_eta_cm" + myVarSuffix));
                combinedNthJetPhi.push_back(tr.getVar<double>("combined"   + std::to_string(nMVAJets.back()) + "thToLastJet" + flavName + "_phi_cm" + myVarSuffix));
                combinedNthJetM.push_back(tr.getVar<double>("combined"     + std::to_string(nMVAJets.back()) + "thToLastJet" + flavName + "_m_cm" + myVarSuffix));
                combinedNthJetE.push_back(tr.getVar<double>("combined"     + std::to_string(nMVAJets.back()) + "thToLastJet" + flavName + "_E_cm" + myVarSuffix));

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

                // ---------------------------------------------------------
                // Here assume number cm jets is same in CR and SR selection
                // ---------------------------------------------------------
                std::vector<double> tempJets_cm_flavb;  
                std::vector<double> tempJets_cm_flavc;  
                std::vector<double> tempJets_cm_flavg;  
                std::vector<double> tempJets_cm_flavuds;
                std::vector<double> tempJets_cm_flavq;  
                
                for (unsigned int iJet = 1; iJet <= nMVAJets.back(); iJet++) {
                    tempJets_cm_flavb.push_back(tr.getVar<double>("Jet_flavb_"     + std::to_string(iJet) + myVarSuffix));
                    tempJets_cm_flavc.push_back(tr.getVar<double>("Jet_flavc_"     + std::to_string(iJet) + myVarSuffix));
                    tempJets_cm_flavg.push_back(tr.getVar<double>("Jet_flavg_"     + std::to_string(iJet) + myVarSuffix));
                    tempJets_cm_flavuds.push_back(tr.getVar<double>("Jet_flavuds_" + std::to_string(iJet) + myVarSuffix));
                    tempJets_cm_flavq.push_back(tr.getVar<double>("Jet_flavq_"     + std::to_string(iJet) + myVarSuffix));
                }

                Jets_cm_flavb.push_back(tempJets_cm_flavb);
                Jets_cm_flavc.push_back(tempJets_cm_flavc);
                Jets_cm_flavg.push_back(tempJets_cm_flavg);
                Jets_cm_flavuds.push_back(tempJets_cm_flavuds);
                Jets_cm_flavq.push_back(tempJets_cm_flavq);

                Jets_flavb.push_back(tr.getVec<float>("Jets" + myVarSuffix + "_bJetTagDeepFlavourtotb"));
                Jets_flavc.push_back(tr.getVec<float>("Jets" + myVarSuffix + "_bJetTagDeepFlavourprobc"));
                Jets_flavg.push_back(tr.getVec<float>("Jets" + myVarSuffix + "_bJetTagDeepFlavourprobg"));
                Jets_flavuds.push_back(tr.getVec<float>("Jets" + myVarSuffix + "_bJetTagDeepFlavourprobuds"));
                Jets_flavq.push_back(tr.getVec<float>("Jets" + myVarSuffix + "_bJetTagDeepFlavourtotq"));

                std::map<std::string, std::vector<bool> > tempRegionMap;
                std::map<std::string, std::vector<bool> > tempRegionMap_CR;

                regions.push_back(tr.getVec<std::string>("regions_" + channel + "l_RPV"));
                for (const std::string& region : regions.back())
                {
                    tempRegionMap[region]    = tr.getVec<bool>("DoubleDisCo_" + region + "_"            + channel + "l_RPV" + myVarSuffix); 
                    tempRegionMap_CR[region] = tr.getVec<bool>("DoubleDisCo_" + region + "_NonIsoMuon_" + channel + "l_RPV" + myVarSuffix); 
                }

                // --------------------------------------
                // Now fill up the CR vector of variables
                // --------------------------------------
                DoubleDisCo_passRegions.push_back(tempRegionMap); 
                DoubleDisCo_passRegions_CR.push_back(tempRegionMap_CR); 

                Jets_CR.push_back(tr.getVec<utility::LorentzVector>("Jets"         + myVarSuffix));
                Jets_cm_top6_CR.push_back(tr.getVec<utility::LorentzVector>(mvaName + "Jets_cm_top6" + myVarSuffix));

                GoodJets_CR.push_back(tr.getVec<bool>(jetsName + "Jets_pt30"        + myVarSuffix));
                NGoodJets_CR.push_back(tr.getVar<int>("N" + jetsName + "Jets_pt30" + myVarSuffix));
                NGoodBJets_CR.push_back(tr.getVar<int>("NGoodBJets_pt30" + myVarSuffix));
                NTops_CR.push_back(tr.getVar<int>("ntops" + myVarSuffix));
                NRTops_CR.push_back(tr.getVar<int>("ntops_3jet" + myVarSuffix));
                NMTops_CR.push_back(tr.getVar<int>("ntops_1jet" + myVarSuffix));

                HT_pt30_CR.push_back(tr.getVar<double>("HT_" + htName + "_pt30"    + myVarSuffix));
                Baseline_CR.push_back(tr.getVar<bool>("pass_qcdCR"                 + myVarSuffix));
                Baseline_CR.push_back(tr.getVar<bool>("pass_qcdCR_1b"              + myVarSuffix));
                Baseline_CR.push_back(tr.getVar<bool>("pass_qcdCR_1t"              + myVarSuffix));
                Baseline_CR.push_back(tr.getVar<bool>("pass_qcdCR_1b_1t"           + myVarSuffix));
                Baseline_CR.push_back(tr.getVar<bool>("pass_qcdCR_2b"              + myVarSuffix));
                Mbl_CR.push_back(tr.getVar<double>("Mbl"                           + myVarSuffix));
                dRbjets_CR.push_back(tr.getVar<double>("dR_bjets"                  + myVarSuffix));

                DoubleDisCo_massReg_CR.push_back(tr.getVar<double>("DoubleDisCo_massReg_NonIsoMuon_" + channel + "l_RPV"  + myVarSuffix));
                DoubleDisCo_disc1_CR.push_back(tr.getVar<double>("DoubleDisCo_disc1_NonIsoMuon_" + channel + "l_RPV"      + myVarSuffix));
                DoubleDisCo_disc2_CR.push_back(tr.getVar<double>("DoubleDisCo_disc2_NonIsoMuon_" + channel + "l_RPV"      + myVarSuffix));

                fwm2_top6_CR.push_back(tr.getVar<double>(mvaName + "fwm2_top6"         + myVarSuffix));
                fwm3_top6_CR.push_back(tr.getVar<double>(mvaName + "fwm3_top6"         + myVarSuffix));
                fwm4_top6_CR.push_back(tr.getVar<double>(mvaName + "fwm4_top6"         + myVarSuffix));
                fwm5_top6_CR.push_back(tr.getVar<double>(mvaName + "fwm5_top6"         + myVarSuffix));
                jmt_ev0_top6_CR.push_back(tr.getVar<double>(mvaName + "jmt_ev0_top6"   + myVarSuffix));
                jmt_ev1_top6_CR.push_back(tr.getVar<double>(mvaName + "jmt_ev1_top6"   + myVarSuffix));
                jmt_ev2_top6_CR.push_back(tr.getVar<double>(mvaName + "jmt_ev2_top6"   + myVarSuffix));
                nMVAJets_CR.push_back(tr.getVar<unsigned int>(mvaName + "nMVAJets"                  + myVarSuffix));
                nMVALeptons_CR.push_back(tr.getVar<unsigned int>(mvaName + "nMVALeptons"            + myVarSuffix));

                combinedNthJetPt_CR.push_back(tr.getVar<double>("combined"    + std::to_string(nMVAJets_CR.back()) + "thToLastJet" + flavName + "_pt_cm" + myVarSuffix));
                combinedNthJetPtrHT_CR.push_back(tr.getVar<double>("combined" + std::to_string(nMVAJets_CR.back()) + "thToLastJet" + flavName + "_ptrHT_cm" + myVarSuffix));
                combinedNthJetEta_CR.push_back(tr.getVar<double>("combined"   + std::to_string(nMVAJets_CR.back()) + "thToLastJet" + flavName + "_eta_cm" + myVarSuffix));
                combinedNthJetPhi_CR.push_back(tr.getVar<double>("combined"   + std::to_string(nMVAJets_CR.back()) + "thToLastJet" + flavName + "_phi_cm" + myVarSuffix));
                combinedNthJetM_CR.push_back(tr.getVar<double>("combined"     + std::to_string(nMVAJets_CR.back()) + "thToLastJet" + flavName + "_m_cm" + myVarSuffix));
                combinedNthJetE_CR.push_back(tr.getVar<double>("combined"     + std::to_string(nMVAJets_CR.back()) + "thToLastJet" + flavName + "_E_cm" + myVarSuffix));

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

                std::vector<double> tempJets_cm_flavb_CR;
                std::vector<double> tempJets_cm_flavc_CR;
                std::vector<double> tempJets_cm_flavg_CR;
                std::vector<double> tempJets_cm_flavuds_CR;
                std::vector<double> tempJets_cm_flavq_CR;

                for (unsigned int iJet = 1; iJet <= nMVAJets_CR.back(); iJet++) {
                    tempJets_cm_flavb_CR.push_back(  tr.getVar<double>("Jet" + flavName + "_flavb_"  +std::to_string(iJet)+myVarSuffix));
                    tempJets_cm_flavc_CR.push_back(  tr.getVar<double>("Jet" + flavName + "_flavc_"  +std::to_string(iJet)+myVarSuffix));
                    tempJets_cm_flavg_CR.push_back(  tr.getVar<double>("Jet" + flavName + "_flavg_"  +std::to_string(iJet)+myVarSuffix));
                    tempJets_cm_flavuds_CR.push_back(tr.getVar<double>("Jet" + flavName + "_flavuds_"+std::to_string(iJet)+myVarSuffix));
                    tempJets_cm_flavq_CR.push_back(  tr.getVar<double>("Jet" + flavName + "_flavq_"  +std::to_string(iJet)+myVarSuffix));
                }

                Jets_cm_flavb_CR.push_back(tempJets_cm_flavb_CR);
                Jets_cm_flavc_CR.push_back(tempJets_cm_flavc_CR);
                Jets_cm_flavg_CR.push_back(tempJets_cm_flavg_CR);
                Jets_cm_flavuds_CR.push_back(tempJets_cm_flavuds_CR);
                Jets_cm_flavq_CR.push_back(tempJets_cm_flavq_CR);

                Jets_flavb_CR.push_back(tr.getVec<float>("Jets" + myVarSuffix + "_bJetTagDeepFlavourtotb"));
                Jets_flavc_CR.push_back(tr.getVec<float>("Jets" + myVarSuffix + "_bJetTagDeepFlavourprobc"));
                Jets_flavg_CR.push_back(tr.getVec<float>("Jets" + myVarSuffix + "_bJetTagDeepFlavourprobg"));
                Jets_flavuds_CR.push_back(tr.getVec<float>("Jets" + myVarSuffix + "_bJetTagDeepFlavourprobuds"));
                Jets_flavq_CR.push_back(tr.getVec<float>("Jets" + myVarSuffix + "_bJetTagDeepFlavourtotq"));

                // ---------------------------------------------------------------------------
                // Calculate the event weights for tt and ttvars fsrUp/Down, isrUp/Down
                //  -- fsr/isr are event weight based variations
                //  -- calculate event weights for them
                //  -- make histograms for tt and ttvars fsrUp/Down, isrUp/Down in TT root file
                // ----------------------------------------------------------------------------
                double theWeight = 1.0;
                double theWeight_fsrUp = 1.0, theWeight_fsrDown = 1.0;
                double theWeight_isrUp = 1.0, theWeight_isrDown = 1.0;
                if(runtype == "MC" )
                {
                    // tt event weight
                    auto totalWeight = tr.getVar<double>("TotalWeight" + myVarSuffix);
                    theWeight *= totalWeight;
                    // ttvars fsrUp/Down, isrUp/Down event weights
                    const auto& PSweight_fsrUp   = tr.getVar<double>("PSweight_FSRUp"   + myVarSuffix);
                    const auto& PSweight_fsrDown = tr.getVar<double>("PSweight_FSRDown" + myVarSuffix);
                    const auto& PSweight_isrUp   = tr.getVar<double>("PSweight_ISRUp"   + myVarSuffix);
                    const auto& PSweight_isrDown = tr.getVar<double>("PSweight_ISRDown" + myVarSuffix);
                    theWeight_fsrUp   *= totalWeight * PSweight_fsrUp;
                    theWeight_fsrDown *= totalWeight * PSweight_fsrDown;
                    theWeight_isrUp   *= totalWeight * PSweight_isrUp;
                    theWeight_isrDown *= totalWeight * PSweight_isrDown;
                }
                
                // baseline
                weight.push_back(theWeight);
                weight_fsrUp.push_back(theWeight_fsrUp);
                weight_fsrDown.push_back(theWeight_fsrDown);
                weight_isrUp.push_back(theWeight_isrUp);
                weight_isrDown.push_back(theWeight_isrDown);
                // qcd cr
                weight_CR.push_back(theWeight);
                weight_CR_fsrUp.push_back(theWeight_fsrUp);
                weight_CR_fsrDown.push_back(theWeight_fsrDown);
                weight_CR_isrUp.push_back(theWeight_isrUp);
                weight_CR_isrDown.push_back(theWeight_isrDown);
            }

            Debug("Preparing cut map", __LINE__);

            // Define subregions of the 2D disco plane for studying populations
            bool leftEdge   = DoubleDisCo_disc1[1] < 0.1;
            bool bottomEdge = DoubleDisCo_disc2[1] < 0.1;
            bool central    = !leftEdge and !bottomEdge and DoubleDisCo_disc1[1] < 0.8 and DoubleDisCo_disc2[1] < 0.8;
            bool outside    = !central and !leftEdge and !bottomEdge;

            const std::map<std::string, bool> cut_map
            {
                {"_0l"               , Baseline[0]},
                {"_1l"               , Baseline[1]},                         
                {"_1l_leftEdge"      , Baseline[1] && leftEdge},                         
                {"_1l_bottomEdge"    , Baseline[1] && bottomEdge},                         
                {"_1l_central"       , Baseline[1] && central},                         
                {"_1l_outside"       , Baseline[1] && outside},                         
                {"_2l"               , Baseline[2]}, 
                {"_0l_blind"         , Baseline_blind[0]},
                {"_1l_blind"         , Baseline_blind[1]},                         
                {"_2l_blind"         , Baseline_blind[2]}, 
                {"_0l_QCDCR"         , Baseline_CR[0]}, 
                {"_0l_QCDCR_1b"      , Baseline_CR[1]}, 
                {"_0l_QCDCR_1t"      , Baseline_CR[2]}, 
                {"_0l_QCDCR_1b_1t"   , Baseline_CR[3]}, 
                {"_0l_QCDCR_2b"      , Baseline_CR[4]},   
                {"_1l_QCDCR"         , Baseline_CR[0]}, 
                {"_2l_QCDCR"         , Baseline_CR[0]}, 

            };

            // Put assume 7 jets and 2 leptons for making the histograms
            if(!initHistos)
            {
                Debug("Initializing the histograms", __LINE__);

                Preinit(7, 2);
                InitHistos(cut_map, regions);
                initHistos = true;
            }

            std::map<std::string, bool>               njetsMap;
            std::map<std::string, std::vector<bool> > ABCDmap;

            for(auto& kv : cut_map)
            {
                Debug("Top of cut loop for cut: " + kv.first, __LINE__);

                bool isQCD = false;
                // Extract "0", "1", "2", or "Q" from cut string e.g. _1l
                int channel = 1;
                if (kv.first.size() > 0)
                {
                    // Take first character after assumed underscore
                    std::string chunk = kv.first.substr(1,1);
                    channel = std::stoi(chunk);
                   
                    // One control region for the moment, so pick any of three channels
                    if (kv.first.find("QCDCR") != std::string::npos)
                    {
                        isQCD = true;
                    }
                }

                unsigned int nJets = !isQCD ? Jets_cm_top6[channel].size() : Jets_cm_top6_CR[channel].size();

                // For 0L, we always use the NGoodJets case
                njetsMap = {
                    {"Incl",                                                         true },
                    {"6",      !isQCD ? NGoodJets[channel]==6  : NGoodJets_CR[channel]==6 },
                    {"7",      !isQCD ? NGoodJets[channel]==7  : NGoodJets_CR[channel]==7 },
                    {"8",      !isQCD ? NGoodJets[channel]==8  : NGoodJets_CR[channel]==8 },
                    {"9",      !isQCD ? NGoodJets[channel]==9  : NGoodJets_CR[channel]==9 },
                    {"10",     !isQCD ? NGoodJets[channel]==10 : NGoodJets_CR[channel]==10},
                    {"11",     !isQCD ? NGoodJets[channel]==11 : NGoodJets_CR[channel]==11},
                    {"12incl", !isQCD ? NGoodJets[channel]>=12 : NGoodJets_CR[channel]>=12},
                };

                ABCDmap = {};
                for (const auto& region : regions[channel])
                    ABCDmap[region] = !isQCD ? DoubleDisCo_passRegions[channel][region] : DoubleDisCo_passRegions_CR[channel][region];

                double ht = !isQCD ? HT_pt30[channel] : HT_pt30_CR[channel];

                // ---------------
                // loop over njets
                // ---------------
                std::string name;
                for (auto& njetsPass : njetsMap)
                {
                    Debug("Top of njet loop for njets: " + njetsPass.first, __LINE__);

                    std::string njets      = njetsPass.first;
                    bool        inNjetsBin = njetsPass.second;
                    std::string njetsStr = "";
                    if (njets != "Incl")  
                        njetsStr = "_Njets" + njets;
    
                    // --------------------------------------------------------------------------            
                    // loop over the tt and ttvars fsrUp/Down, isrUp/Down
                    // to make indivudual histograms with the label tt and ttvars in TT root file
                    // --------------------------------------------------------------------------
                    for (const auto& ttvar : ttvars)
                    {   
                        Debug("Top of ttvar loop for ttvar: " + ttvar, __LINE__);

                        std::string ttvarStr = "";
                        if (ttvar != "nom")
                            ttvarStr = "_" + ttvar; 

                        if (kv.first.find("QCDCR_") != std::string::npos and ttvar != "nom")
                            continue;

                        double w = 1.0;

                        // get the event weights for fsrUp/Down, isrUp/Down
                        if (ttvar == "nom")
                        {
                            w  = !isQCD ? weight[channel] : weight_CR[channel];
                        }
                        
                        else if (ttvar == "fsrUp")
                        {
                            w = !isQCD ? weight_fsrUp[channel] : weight_CR_fsrUp[channel];
                        }

                        else if (ttvar == "fsrDown")
                        {
                            w = !isQCD ? weight_fsrDown[channel] : weight_CR_fsrDown[channel];
                        }
                        
                        else if (ttvar == "isrUp")
                        {
                            w = !isQCD ? weight_isrUp[channel] : weight_CR_isrUp[channel];
                        }

                        else if (ttvar == "isrDown")
                        {
                            w = !isQCD ? weight_isrDown[channel] : weight_CR_isrDown[channel];
                        }

                        // -----------------
                        // loop over regions
                        // -----------------
                        for (auto& regionPass : ABCDmap)
                        {
                            Debug("Top of region loop for region: " + regionPass.first, __LINE__);

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

                                name = kv.first + njetsStr + regionStr + ttvarStr;
                                if (njets == "Incl" and region == "ABCD")
                                {

                                    Debug("Filling event variable histograms with name: " + name, __LINE__);
                                    my_histos["h_njets"      + name]->Fill(!isQCD ?    NGoodJets[channel] : NGoodJets_CR[channel],    w);
                                    my_histos["h_nbjets"     + name]->Fill(!isQCD ?   NGoodBJets[channel] : NGoodBJets_CR[channel],   w);
                                    my_histos["h_ntops"      + name]->Fill(!isQCD ?        NTops[channel] : NTops_CR[channel],        w);
                                    my_histos["fwm2_top6"    + name]->Fill(!isQCD ?    fwm2_top6[channel] : fwm2_top6_CR[channel],    w);
                                    my_histos["fwm3_top6"    + name]->Fill(!isQCD ?    fwm3_top6[channel] : fwm3_top6_CR[channel],    w);
                                    my_histos["fwm4_top6"    + name]->Fill(!isQCD ?    fwm4_top6[channel] : fwm4_top6_CR[channel],    w);
                                    my_histos["fwm5_top6"    + name]->Fill(!isQCD ?    fwm5_top6[channel] : fwm5_top6_CR[channel],    w);
                                    my_histos["jmt_ev0_top6" + name]->Fill(!isQCD ? jmt_ev0_top6[channel] : jmt_ev0_top6_CR[channel], w);
                                    my_histos["jmt_ev1_top6" + name]->Fill(!isQCD ? jmt_ev1_top6[channel] : jmt_ev1_top6_CR[channel], w);
                                    my_histos["jmt_ev2_top6" + name]->Fill(!isQCD ? jmt_ev2_top6[channel] : jmt_ev2_top6_CR[channel], w);

                                    my_histos["h_mbl"        + name]->Fill(!isQCD ? Mbl[channel] : Mbl_CR[channel], w);
                                    my_histos["h_dRbjets"    + name]->Fill(!isQCD ? dRbjets[channel] : dRbjets_CR[channel], w);

                                    Debug("Filling stop variable histograms with name: " + name, __LINE__);

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
                                    
                                    Debug("Filling jet variable histograms with name: " + name, __LINE__);

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

                                    auto& theGoodLeptons  = GoodLeptons[channel];
                                    auto& theGoodLeptonsC = GoodLeptonsCharge[channel];
                                    auto& theGoodLeptonsI = GoodLeptonsMiniIso[channel];
                                    if (isQCD)
                                    {
                                       theGoodLeptons  = GoodLeptons_CR[channel]; 
                                       theGoodLeptonsC = GoodLeptonsCharge_CR[channel]; 
                                       theGoodLeptonsI = GoodLeptonsMiniIso_CR[channel]; 
                                    }

                                    Debug("Filling lepton variable histograms with name: " + name, __LINE__);

                                    for(std::size_t j = 0 ; j < theGoodLeptons.size(); ++j)
                                    {
                                        std::string type            = theGoodLeptons[j].first;
                                        utility::LorentzVector lvec = theGoodLeptons[j].second; 
                                        double charge               = theGoodLeptonsC[j];
                                        double iso                  = theGoodLeptonsI[j];

                                        my_histos["h_lPt_"      + std::to_string(j+1) + name]->Fill(lvec.Pt(),  w);
                                        my_histos["h_lPhi_"     + std::to_string(j+1) + name]->Fill(lvec.Phi(), w);
                                        my_histos["h_lEta_"     + std::to_string(j+1) + name]->Fill(lvec.Eta(), w);
                                        my_histos["h_lCharge_"  + std::to_string(j+1) + name]->Fill(charge ,    w);
                                        my_histos["h_lMiniIso_" + std::to_string(j+1) + name]->Fill(iso ,       w);

                                        if(type == 'e'){
                                            my_histos["h_ePt_"      + std::to_string(j+1) + name]->Fill(lvec.Pt() , w);
                                            my_histos["h_ePhi_"     + std::to_string(j+1) + name]->Fill(lvec.Phi(), w);
                                            my_histos["h_eEta_"     + std::to_string(j+1) + name]->Fill(lvec.Eta(), w);
                                            my_histos["h_eCharge_"  + std::to_string(j+1) + name]->Fill(charge,     w);
                                            my_histos["h_eMiniIso_" + std::to_string(j+1) + name]->Fill(iso,        w);
                                        } else if (type == 'm' or type == 'n') {
                                            my_histos["h_mPt_"      + std::to_string(j+1) + name]->Fill(lvec.Pt(),  w);
                                            my_histos["h_mPhi_"     + std::to_string(j+1) + name]->Fill(lvec.Phi(), w);
                                            my_histos["h_mEta_"     + std::to_string(j+1) + name]->Fill(lvec.Eta(), w);
                                            my_histos["h_mCharge_"  + std::to_string(j+1) + name]->Fill(charge,     w);
                                            my_histos["h_mMiniIso_" + std::to_string(j+1) + name]->Fill(iso,        w);
                                        }
                                    }

                                    auto& theJets         = Jets[channel];
                                    auto& theGoodJets     = GoodJets[channel];
                                    auto& theJetsFlavb    = Jets_flavb[channel];
                                    auto& theJetsFlavc    = Jets_flavc[channel];
                                    auto& theJetsFlavg    = Jets_flavg[channel];
                                    auto& theJetsFlavuds  = Jets_flavuds[channel];
                                    auto& theJetsFlavq    = Jets_flavq[channel];

                                    auto& nMVAjets        = nMVAJets[channel]; 
                                    auto& theCombJetPt    = combinedNthJetPt[channel];
                                    auto& theCombJetPtrHT = combinedNthJetPtrHT[channel];
                                    auto& theCombJetEta   = combinedNthJetEta[channel];
                                    auto& theCombJetPhi   = combinedNthJetPhi[channel];
                                    auto& theCombJetE     = combinedNthJetE[channel];
                                    auto& theCombJetM     = combinedNthJetM[channel];

                                    if (isQCD)
                                    {
                                        theJets         = Jets_CR[channel];
                                        theGoodJets     = GoodJets_CR[channel];
                                        nMVAjets        = nMVAJets_CR[channel];
                                        theCombJetPt    = combinedNthJetPt_CR[channel];
                                        theCombJetPtrHT = combinedNthJetPtrHT_CR[channel];
                                        theCombJetEta   = combinedNthJetEta_CR[channel];
                                        theCombJetPhi   = combinedNthJetPhi_CR[channel];
                                        theCombJetE     = combinedNthJetE_CR[channel];
                                        theCombJetM     = combinedNthJetM_CR[channel];
                                        theJetsFlavb    = Jets_flavb_CR[channel];
                                        theJetsFlavc    = Jets_flavc_CR[channel];
                                        theJetsFlavg    = Jets_flavg_CR[channel];
                                        theJetsFlavuds  = Jets_flavuds_CR[channel];
                                        theJetsFlavq    = Jets_flavq_CR[channel];
                                    }

                                    Debug("Filling non-cm jet variable histograms with name: " + name, __LINE__);

                                    unsigned int iGoodJet = 1;
                                    for(unsigned int j = 0; j < theJets.size(); j++) 
                                    {

                                        if(!theGoodJets[j]) continue;
                                        my_histos["h_jPt_"  + std::to_string(iGoodJet) + name]->Fill(theJets.at(j).Pt(),  w); 
                                        my_histos["h_jEta_" + std::to_string(iGoodJet) + name]->Fill(theJets.at(j).Eta(), w); 
                                        my_histos["h_jPhi_" + std::to_string(iGoodJet) + name]->Fill(theJets.at(j).Phi(), w); 
                                        my_histos["h_jM_"   + std::to_string(iGoodJet) + name]->Fill(theJets.at(j).M(),   w); 
                                        my_histos["h_jE_"   + std::to_string(iGoodJet) + name]->Fill(theJets.at(j).E(),   w); 

                                        my_histos["h_jFlavb_"   + std::to_string(iGoodJet) + name]->Fill(theJetsFlavb.at(j),   w); 
                                        my_histos["h_jFlavc_"   + std::to_string(iGoodJet) + name]->Fill(theJetsFlavc.at(j),   w); 
                                        my_histos["h_jFlavg_"   + std::to_string(iGoodJet) + name]->Fill(theJetsFlavg.at(j),   w); 
                                        my_histos["h_jFlavuds_" + std::to_string(iGoodJet) + name]->Fill(theJetsFlavuds.at(j), w); 
                                        my_histos["h_jFlavq_"   + std::to_string(iGoodJet) + name]->Fill(theJetsFlavq.at(j),   w); 

                                        iGoodJet++;

                                        if (iGoodJet > nMVAjets)
                                            break;
                                    }                        

                                    Debug("Filling combined jet variable histograms with name: " + name, __LINE__);

                                    my_histos["h_combined" + std::to_string(nMVAjets) + "thJet_pt"    + name]->Fill(theCombJetPt,    w); 
                                    my_histos["h_combined" + std::to_string(nMVAjets) + "thJet_ptrHT" + name]->Fill(theCombJetPtrHT, w); 
                                    my_histos["h_combined" + std::to_string(nMVAjets) + "thJet_eta"   + name]->Fill(theCombJetEta,   w); 
                                    my_histos["h_combined" + std::to_string(nMVAjets) + "thJet_phi"   + name]->Fill(theCombJetPhi,   w); 
                                    my_histos["h_combined" + std::to_string(nMVAjets) + "thJet_m"     + name]->Fill(theCombJetM,     w); 
                                    my_histos["h_combined" + std::to_string(nMVAjets) + "thJet_E"     + name]->Fill(theCombJetE,     w); 
                                }                            

                                if (region == "ABCD" and ttvar == "nom")
                                    my_2d_histos["h_nRtops_vs_nMtops"      + name]->Fill(!isQCD ? NRTops[channel] : NRTops_CR[channel], !isQCD ? NMTops[channel] : NMTops_CR[channel], w);

                                auto& cmJets = Jets_cm_top6[channel];
                                if (isQCD)
                                    cmJets = Jets_cm_top6_CR[channel];

                                Debug("Filling 2D pt HT variable histograms with name: " + name, __LINE__);

                                if (region == "ABCD" and ttvar == "nom")
                                {
                                    for(unsigned int i = 1; i <= nJets; i++)
                                    {
                                        double pt  = static_cast<double>(cmJets.at(i-1).Pt());
                                        my_2d_histos["h_cm_pt_jetRank"    + name]->Fill(pt,    i, w);
                                        my_2d_histos["h_cm_ptrHT_jetRank" + name]->Fill(pt/ht, i, w);
                                    }
                                }

                                Debug("Filling 2D disc1 disc2 histograms with name: " + name, __LINE__);

                                if (region == "ABCD" and njets != "Incl")
                                    my_histos["h_DoubleDisCo_massReg" + name]->Fill(!isQCD ? DoubleDisCo_massReg[channel] : DoubleDisCo_massReg_CR[channel], w);

                                auto& disc1 = DoubleDisCo_disc1[channel];
                                auto& disc2 = DoubleDisCo_disc2[channel];
                                if (isQCD)
                                {
                                    disc1 = DoubleDisCo_disc1_CR[channel];
                                    disc2 = DoubleDisCo_disc2_CR[channel];
                                }
                                Debug("Filling 1D discriminants histograms with name: " + name, __LINE__);
                                
                                // -----------------------------------------------------------
                                // make inputs for Validation Study
                                //  -- it is a 2D histogram of whole ABCD region for per-njets 
                                // -----------------------------------------------------------
                                if (njets != "Incl" and region == "ABCD")
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
                                    Debug("Filling 1D njets histograms with name: " + name, __LINE__);

                                    name = kv.first + njetsStr + regionStr;
                                    int shift = inRegionBin[iSubRegion] ? iSubRegion : 100;

                                    auto& theGoodJets = NGoodJets[channel];
                                    auto& theMVAjets  = nMVAJets[channel]; 
                                    if (isQCD)
                                    {
                                        theGoodJets = NGoodJets_CR[channel];
                                        theMVAjets  = nMVAJets_CR[channel];
                                    }

                                    // ---------------------------------------------------------
                                    // make inputs for Higgs Combine
                                    //  -- it is 1D histogram, including 4 groups A, B, C, D
                                    //  -- each of A, B, C, D has 7,8,9,10,11, 12incl Njets bins
                                    // ---------------------------------------------------------
                                    if (njets == "Incl")
                                        my_histos["h_njets_12incl" + name]->Fill(theGoodJets>=12 ? 12-theMVAjets+shift*6 : theGoodJets-theMVAjets+shift*6, w);

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
