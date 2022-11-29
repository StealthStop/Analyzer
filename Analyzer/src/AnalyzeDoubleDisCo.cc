#define AnalyzeDoubleDisCo_cxx
#include "Analyzer/Analyzer/include/AnalyzeDoubleDisCo.h"
#include "NTupleReader/include/NTupleReader.h"
#include "Framework/Framework/include/Utility.h"

#include "Framework/Framework/include/Muon.h"
#include "Framework/Framework/include/Electron.h"
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
        {"h_DoubleDisCo_disc1",            80,    0,    1},
        {"h_DoubleDisCo_disc2",            80,    0,    1},
        {"h_DoubleDisCo_massReg",         180,    0, 1800},
        {"h_FWM2_top6",                    50,    0,    1},
        {"h_FWM3_top6",                    50,    0,    1},
        {"h_FWM4_top6",                    50,    0,    1},
        {"h_FWM5_top6",                    50,    0,    1},
        {"h_JMT_ev0_top6",                 50,    0,    1},
        {"h_JMT_ev1_top6",                 50,    0,    1},
        {"h_JMT_ev2_top6",                 50,    0,    1},
        {"h_Stop1_Pt_cm_OldSeed",         360,    0, 1800},
        {"h_Stop1_Eta_cm_OldSeed",         80,   -6,    6},
        {"h_Stop1_Phi_cm_OldSeed",         64,   -4,    4},
        {"h_Stop1_Mass_cm_OldSeed",       180,    0, 1800},
        {"h_Stop2_Pt_cm_OldSeed",         360,    0, 1800},
        {"h_Stop2_Eta_cm_OldSeed",         80,   -6,    6},
        {"h_Stop2_Phi_cm_OldSeed",         64,   -4,    4},
        {"h_Stop2_Mass_cm_OldSeed",       180,    0, 1800},
        {"h_Stop1_Pt_cm_TopSeed",         360,    0, 1800},
        {"h_Stop1_Eta_cm_TopSeed",         80,   -6,    6},
        {"h_Stop1_Phi_cm_TopSeed",         64,   -4,    4},
        {"h_Stop1_Mass_cm_TopSeed",       180,    0, 1800},
        {"h_Stop2_Pt_cm_TopSeed",         360,    0, 1800},
        {"h_Stop2_Eta_cm_TopSeed",         80,   -6,    6},
        {"h_Stop2_Phi_cm_TopSeed",         64,   -4,    4},
        {"h_Stop2_Mass_cm_TopSeed",       180,    0, 1800},
        {"h_Top1_Pt_cm",                  360,    0, 1800},
        {"h_Top1_Eta_cm",                  80,   -6,    6},
        {"h_Top1_Phi_cm",                  64,   -4,    4},
        {"h_Top1_Mass_cm",                180,    0, 1800},
        {"h_Top2_Pt_cm",                  360,    0, 1800},
        {"h_Top2_Eta_cm",                  80,   -6,    6},
        {"h_Top2_Phi_cm",                  64,   -4,    4},
        {"h_Top2_Mass_cm",                180,    0, 1800},
        {"h_Top1_Pt",                     360,    0, 1800},
        {"h_Top1_Eta",                     80,   -6,    6},
        {"h_Top1_Phi",                     64,   -4,    4},
        {"h_Top1_Mass",                   180,    0, 1800},
        {"h_Top2_Pt",                     360,    0, 1800},
        {"h_Top2_Eta",                     80,   -6,    6},
        {"h_Top2_Phi",                     64,   -4,    4},
        {"h_Top2_Mass",                   180,    0, 1800},
        {"h_Njets",                        21, -0.5, 20.5},
        {"h_Nbjets",                       21, -0.5, 20.5},
        {"h_Ntops",                         7, -0.5,  6.5},
        {"h_njets_11incl",                 24, -0.5, 23.5},
        {"h_njets_12incl",                 24, -0.5, 23.5},
        {"h_njets_13incl",                 24, -0.5, 23.5},
        {"h_HT",                          720,    0, 7200},
        {"h_dRbjets",                     180,    0,    6},
        {"h_Mbl",                         180,    0,  360},
        {"h_Mll",                         180,    0,  360},
        {"h_Stop1_mass_PtRank_matched",   180,    0, 1800},
        {"h_Stop2_mass_PtRank_matched",   180,    0, 1800},
        {"h_Stop1_mass_MassRank_matched", 180,    0, 1800},
        {"h_Stop2_mass_MassRank_matched", 180,    0, 1800},
        {"h_Stop_mass_average_matched",   180,    0, 1800},
    };

    hist2DInfos = {
        {"h_DoubleDisCo_disc1_disc2", 100,    0,    1, 100,    0,   1}, 
        {"h_pt_jetRank_cm",           150,    0, 1500,  10,    0,  10},
        {"h_ptrHT_jetRank_cm",        150,    0,    1,  10,    0,  10},
        {"h_nRtops_vs_nMtops",          7, -0.5,  6.5,   7, -0.5, 6.5},
    };

    channels = {"0", "1", "2"};
    njets    = {"Incl", "6", "7", "8", "9", "10", "11", "11incl", "12incl", "13incl"};
    systvars = {"nom", "fsrUp", "fsrDown", "isrUp", "isrDown"   
                       "pdfUp", "pdfDown", "sclUp", "sclDown",
                       "prfUp", "prfDown", "btgUp", "btgDown",
                       "jetUp", "jetDown", "lepUp", "lepDown",
                       "puUp",  "puDown"   
    };
    jecvars  = {"", "JECup", "JECdown", "JERup", "JERdown"      };

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
        histInfos.push_back({"h_Jet" + std::to_string(i) + "_cm_Pt"      , 360,  0,  1800});
        histInfos.push_back({"h_Jet" + std::to_string(i) + "_cm_PtrHT"   , 180,  0,     1});
        histInfos.push_back({"h_Jet" + std::to_string(i) + "_cm_Eta"     ,  80, -6,     6});
        histInfos.push_back({"h_Jet" + std::to_string(i) + "_cm_Phi"     ,  64, -4,     4});
        histInfos.push_back({"h_Jet" + std::to_string(i) + "_cm_Mass"    , 180,  0,   360});
        histInfos.push_back({"h_Jet" + std::to_string(i) + "_cm_Energy"  , 360,  0,  1800});

        histInfos.push_back({"h_Jet" + std::to_string(i) + "_cm_Flavb"   ,  80,  0,     1});
        histInfos.push_back({"h_Jet" + std::to_string(i) + "_cm_Flavc"   ,  80,  0,     1});
        histInfos.push_back({"h_Jet" + std::to_string(i) + "_cm_Flavg"   ,  80,  0,     1});
        histInfos.push_back({"h_Jet" + std::to_string(i) + "_cm_Flavq"   ,  80,  0,     1});
        histInfos.push_back({"h_Jet" + std::to_string(i) + "_cm_Flavuds" ,  80,  0,     1});

        histInfos.push_back({"h_Jet" + std::to_string(i) + "_Pt"         , 360,  0,  1800});
        histInfos.push_back({"h_Jet" + std::to_string(i) + "_PtrHT"      , 180,  0,  1800});
        histInfos.push_back({"h_Jet" + std::to_string(i) + "_Eta"        ,  80, -6,     6});
        histInfos.push_back({"h_Jet" + std::to_string(i) + "_Phi"        ,  64, -4,     4});
        histInfos.push_back({"h_Jet" + std::to_string(i) + "_Mass"       , 180,  0,   360});
        histInfos.push_back({"h_Jet" + std::to_string(i) + "_Energy"     , 360,  0,  1800});

        histInfos.push_back({"h_Jet" + std::to_string(i) + "_Flavb"      ,  80,  0,     1});
        histInfos.push_back({"h_Jet" + std::to_string(i) + "_Flavc"      ,  80,  0,     1});
        histInfos.push_back({"h_Jet" + std::to_string(i) + "_Flavg"      ,  80,  0,     1});
        histInfos.push_back({"h_Jet" + std::to_string(i) + "_Flavq"      ,  80,  0,     1});
        histInfos.push_back({"h_Jet" + std::to_string(i) + "_Flavuds"    ,  80,  0,     1});
    }

    // Make one-less-jet combined jet histogram
    histInfos.push_back({"h_combined" + std::to_string(nNNJets-1) + "thJet_Pt",     360,  0, 1800});
    histInfos.push_back({"h_combined" + std::to_string(nNNJets-1) + "thJet_PtrHT",  180,  0,    1});
    histInfos.push_back({"h_combined" + std::to_string(nNNJets-1) + "thJet_Eta",     80, -6,    6});
    histInfos.push_back({"h_combined" + std::to_string(nNNJets-1) + "thJet_Phi",     64, -4,    4});
    histInfos.push_back({"h_combined" + std::to_string(nNNJets-1) + "thJet_Mass",   180,  0,  360});
    histInfos.push_back({"h_combined" + std::to_string(nNNJets-1) + "thJet_Energy", 360,  0, 1800});

    histInfos.push_back({"h_combined" + std::to_string(nNNJets) + "thJet_Pt",     360,  0, 1800});
    histInfos.push_back({"h_combined" + std::to_string(nNNJets) + "thJet_PtrHT",  180,  0,    1});
    histInfos.push_back({"h_combined" + std::to_string(nNNJets) + "thJet_Eta",     80, -6,    6});
    histInfos.push_back({"h_combined" + std::to_string(nNNJets) + "thJet_Phi",     64, -4,    4});
    histInfos.push_back({"h_combined" + std::to_string(nNNJets) + "thJet_Mass",   180,  0,  360});
    histInfos.push_back({"h_combined" + std::to_string(nNNJets) + "thJet_Energy", 360,  0, 1800});

    for(unsigned int i = 1; i <= nLeptons; i++)
    {
        histInfos.push_back({"h_Lepton" + std::to_string(i) + "_cm_Pt"    ,  360,    0, 1800});
        histInfos.push_back({"h_Lepton" + std::to_string(i) + "_cm_Eta"   ,   80,   -6,    6});
        histInfos.push_back({"h_Lepton" + std::to_string(i) + "_cm_Phi"   ,   64,   -4,    4});

        histInfos.push_back({"h_Lepton" + std::to_string(i) + "_Pt"       ,  360,    0, 1800});
        histInfos.push_back({"h_Lepton" + std::to_string(i) + "_Eta"      ,   80,   -6,    6});
        histInfos.push_back({"h_Lepton" + std::to_string(i) + "_Phi"      ,   64,   -4,    4});
        histInfos.push_back({"h_Lepton" + std::to_string(i) + "_Charge"   ,    2,   -1,    1});
        histInfos.push_back({"h_Lepton" + std::to_string(i) + "_MiniIso"  , 1440,    0,    5});

        histInfos.push_back({"h_Electron" + std::to_string(i) + "_Pt"     ,  360,    0, 1800});
        histInfos.push_back({"h_Electron" + std::to_string(i) + "_Eta"    ,   80,   -6,    6});
        histInfos.push_back({"h_Electron" + std::to_string(i) + "_Phi"    ,   64,   -4,    4});
        histInfos.push_back({"h_Electron" + std::to_string(i) + "_Charge" ,    2,   -1,    1});
        histInfos.push_back({"h_Electron" + std::to_string(i) + "_MiniIso", 1440,    0,    5});

        histInfos.push_back({"h_Muon" + std::to_string(i) + "_Pt"         ,  360,    0, 1800});
        histInfos.push_back({"h_Muon" + std::to_string(i) + "_Phi"        ,   80,   -4,    4});
        histInfos.push_back({"h_Muon" + std::to_string(i) + "_Eta"        ,   64,   -6,    6});
        histInfos.push_back({"h_Muon" + std::to_string(i) + "_Charge"     ,    2,   -1,    1});
        histInfos.push_back({"h_Muon" + std::to_string(i) + "_MiniIso"    , 1440,    0,    5});
    }
}

void AnalyzeDoubleDisCo::InitHistos(const std::map<std::string, bool>& cutMap, const std::vector<std::vector<std::string>>& regionsVec, const std::string& runtype)
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
                // loop over the systvars fsrUp/Down, isrUp/Down, etc.
                // to make indivudual histograms with the label systvars in root file
                // --------------------------------------------------------------------------
                for (const auto& ttvar : systvars)
                {
                    // Variations irrelevant for data
                    if (ttvar != "nom" and runtype == "Data")
                        continue;

                    std::string ttvarStr = "";
                    if (ttvar != "nom")
                        ttvarStr = "_" + ttvar; 
                
                    // -------------------------------------------------------------
                    // loop over jecvars to pick up different JEC and JER variations 
                    // -------------------------------------------------------------
                    for (const auto& jecvar : jecvars)
                    {
                        // Variations irrelevant for data
                        if (jecvar != "" and runtype == "Data")
                            continue;

                        // We don't do double variations
                        if (jecvar != "" and ttvar != "nom")
                            continue;

                        std::string jecStr = "";
                        if (jecvar != "")
                            jecStr = "_" + jecvar;

                        // ------------------------------------------------------
                        // loop over the regions    
                        // regions : a vector of region string names for 0L or 1L
                        // ------------------------------------------------------
                        for (const auto& regionPair : subRegionsMap)
                        {
                            std::string region = regionPair.first; 

                            std::string regionStr = "";
                            if (region != "Incl")
                                regionStr = "_" + region;

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

                            if (hInfo.name.find("njets") == std::string::npos and region != "ABCD")
                                continue;

                            if (mycut.first.find("QCDCR_") != std::string::npos and ttvar != "nom")
                                continue; 

                            if ((ttvar != "nom" or jecvar != "") and hInfo.name.find("njets") == std::string::npos and hInfo.name.find("DoubleDisCo") == std::string::npos)
                                continue;

                            std::string name = hInfo.name + mycut.first + njetStr + regionStr + ttvarStr + jecStr;
                            my_histos.emplace(name, std::make_shared<TH1D>((name).c_str(),(name).c_str(), hInfo.nBins, hInfo.low, hInfo.high));

                            Debug("Initializing histogram with name: " + name, __LINE__);
                        } 
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
                // loop over the systvars fsrUp/Down, isrUp/Down
                // to make indivudual histograms with the label systvars in root file
                // --------------------------------------------------------------------------
                for (const auto& ttvar : systvars)
                {   
                    // Variations irrelevant for data
                    if (ttvar != "nom" and runtype == "Data")
                        continue;

                    std::string ttvarStr = "";
                    if (ttvar != "nom")
                        ttvarStr = "_" + ttvar;

                    // -------------------------------------------------------------
                    // loop over jecvars to pick up different JEC and JER variations 
                    // -------------------------------------------------------------
                    for (const auto& jecvar : jecvars)
                    {
                       // Variations irrelevant for data
                       if (jecvar != "" and runtype == "Data")
                           continue;

                       // We don't do double variations
                        if (jecvar != "" and ttvar != "nom")
                            continue;

                        std::string jecStr = "";
                        if (jecvar != "")
                            jecStr = "_" + jecvar;

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

                            if ((ttvar != "nom" or jecvar != "") and h2dInfo.name.find("DoubleDisCo") == std::string::npos)
                                continue;

                            std::string name = h2dInfo.name + mycut.first + njetStr + regionStr + ttvarStr + jecStr;
                            my_2d_histos.emplace(name, std::make_shared<TH2D>((name).c_str(),(name).c_str(), h2dInfo.nBinsX, h2dInfo.lowX, h2dInfo.highX, h2dInfo.nBinsY, h2dInfo.lowY, h2dInfo.highY));

                            Debug("Initializing histogram with name: " + name, __LINE__);
                        } 
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

    const auto& filetag                           = tr.getVar<std::string>("filetag"                          );
    const auto& runtype                           = tr.getVar<std::string>("runtype"                          );    
    const auto& runYear                           = tr.getVar<std::string>("runYear"                          );
    const auto& btagEffFileName                   = tr.getVar<std::string>("btagEffFileName"                  );
    const auto& bjetTagFileName                   = tr.getVar<std::string>("bjetTagFileName"                  );
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

    for(const auto& jecvar : jecvars)
    {
        // Cannot do JEC and JER variations for data
        if (jecvar != "" and runtype == "Data")
            continue;

        Jet                 jet(jecvar);
        BJet                bjet(jecvar);
        Muon                muon(jecvar);
        Baseline            baseline(jecvar);
        Electron            electron(jecvar);
        StopJets            stopJets(jecvar);
        RunTopTagger        topTagger(TopTaggerCfg, jecvar);
        StopGenMatch        stopGenMatch(jecvar);
        CommonVariables     commonVariables(jecvar);
        MakeMVAVariables    makeMVAVariables(                false,  jecvar,        "GoodJets_pt30",       false, true, 7, 2, ""                                  );
        MakeMVAVariables    makeMVAVariables_NonIsoMuon(     false,  jecvar,        "NonIsoMuonJets_pt30", false, true, 7, 2, ""                                  );
        // 0l
        // note that if we make the inputs to the NN, we use just the GoodJets_pt30 collection to derive things
        // but, if we define the QCD CR selection, we use the NonIsoMuonJets_pt30 collection 
        DeepEventShape      neuralNetwork0L(           DoubleDisCo_Cfg_0l_RPV,            DoubleDisCo_Model_0l_RPV, "Info", true, jecvar                          );
        DeepEventShape      neuralNetwork0L_NonIsoMuon(DoubleDisCo_Cfg_NonIsoMuon_0l_RPV, DoubleDisCo_Model_0l_RPV, "Info", true, jecvar                          );
        MakeStopHemispheres stopHemispheres_TopSeed(           "StopJets", "GoodStopJets", "NGoodStopJets", "_TopSeed",            jecvar, Hemisphere::TopSeed    );
        MakeStopHemispheres stopHemispheres_TopSeed_NonIsoMuon("StopJets", "GoodStopJets", "NGoodStopJets", "_TopSeed_NonIsoMuon", jecvar, Hemisphere::InvMassSeed);
        // 1l
        DeepEventShape      neuralNetwork1L(           DoubleDisCo_Cfg_1l_RPV,            DoubleDisCo_Model_1l_RPV, "Info", true, jecvar                                    ); 
        DeepEventShape      neuralNetwork1L_NonIsoMuon(DoubleDisCo_Cfg_NonIsoMuon_1l_RPV, DoubleDisCo_Model_1l_RPV, "Info", true, jecvar                                    );
        MakeStopHemispheres stopHemispheres_OldSeed(           "Jets", "GoodJets_pt20",       "NGoodJets_pt20",       "_OldSeed",            jecvar, Hemisphere::InvMassSeed);
        MakeStopHemispheres stopHemispheres_OldSeed_NonIsoMuon("Jets", "NonIsoMuonJets_pt20", "NNonIsoMuonJets_pt30", "_OldSeed_NonIsoMuon", jecvar, Hemisphere::InvMassSeed);        
        // 2l
        DeepEventShape      neuralNetwork2L(           DoubleDisCo_Cfg_2l_RPV,            DoubleDisCo_Model_2l_RPV, "Info", true, jecvar); 
        DeepEventShape      neuralNetwork2L_NonIsoMuon(DoubleDisCo_Cfg_NonIsoMuon_2l_RPV, DoubleDisCo_Model_2l_RPV, "Info", true, jecvar);
        
        // Remember, order matters here !
        // Follow what is done in Config.h
        tr.registerFunction(muon);
        tr.registerFunction(electron);
        tr.registerFunction(jet);
        tr.registerFunction(bjet);
        tr.registerFunction(commonVariables);
        tr.registerFunction(topTagger);
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
            ScaleFactors        scaleFactors(runYear, leptonFileName, hadronicFileName, meanFileName, jecvar);
            BTagCorrector       bTagCorrector(btagEffFileName, "", bjetTagFileName, "", filetag);
            bTagCorrector.SetVarNames("GenParticles_PdgId", "Jets"+jecvar, "GoodJets_pt30"+jecvar, "Jets"+jecvar+"_bJetTagDeepFlavourtotb", "Jets"+jecvar+"_partonFlavor", jecvar);

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

        for(const auto& jecvar : jecvars)
        {
            // Cannot do JEC and JER variations for data
            if (jecvar != "" and runtype == "Data")
                continue;

            std::string jecStr = "";
            if (jecvar != "")
                jecStr = "_" + jecvar;

            // If an event is not interesting for any channel selection (see explanation and definition
            // of lostCauseEvent boolean in Baseline.h), then move on here and do not waste any more 
            // time on looping over histos or getting vars.
            // This only matters if the user specifies -s on the command line
            // N.B. We already filled the EventCounter histogram for our due-diligence of counting up
            // every single event.
            const auto& lostCauseEvent = tr.getVar<bool>("lostCauseEvent" + jecvar);
            const auto& fastMode       = tr.getVar<bool>("fastMode");

            if (lostCauseEvent and fastMode)
                break;

            // Put 0L and 1L version of variables into vector
            // 0th position is for 0L and 1st position is for 1L for convenience in the event loop
            // For 0L things in the CR, nominal version of jets and derived quantities are used, take note !
            std::vector<std::vector<utility::LorentzVector> >       Jets_CR                        ;
            std::vector<std::vector<utility::LorentzVector> >       Jets_cm_top6_CR                ;
            std::vector<std::vector<utility::LorentzVector> >       Tops_CR                        ;
            std::vector<double>                                     Top1_pt_cm_CR                  ;
            std::vector<double>                                     Top1_eta_cm_CR                 ;
            std::vector<double>                                     Top1_mass_cm_CR                ;
            std::vector<double>                                     Top1_phi_cm_CR                 ;
            std::vector<double>                                     Top2_pt_cm_CR                  ;
            std::vector<double>                                     Top2_eta_cm_CR                 ;
            std::vector<double>                                     Top2_mass_cm_CR                ;
            std::vector<double>                                     Top2_phi_cm_CR                 ;

            std::vector<bool>                                       Baseline_CR                    ;
            std::vector<double>                                     weight_CR                      ;
            std::vector<std::vector<bool> >                         GoodJets_CR                    ;
            std::vector<int>                                        NGoodJets_CR                   ;
            std::vector<int>                                        NGoodBJets_CR                  ;
            std::vector<int>                                        NTops_CR                       ;
            std::vector<int>                                        NRTops_CR                      ;
            std::vector<int>                                        NMTops_CR                      ;
            std::vector<double>                                     Mbl_CR                         ;
            std::vector<double>                                     Mll_CR                         ;
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

            std::vector<double>                                     Lepton1_pt_cm_CR               ;
            std::vector<double>                                     Lepton1_eta_cm_CR              ;
            std::vector<double>                                     Lepton1_phi_cm_CR              ;
            std::vector<double>                                     Lepton2_pt_cm_CR               ;
            std::vector<double>                                     Lepton2_eta_cm_CR              ;
            std::vector<double>                                     Lepton2_phi_cm_CR              ;

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
            std::vector<std::vector<utility::LorentzVector> >       Tops                           ;
            std::vector<double>                                     Top1_pt_cm                     ;
            std::vector<double>                                     Top1_eta_cm                    ;
            std::vector<double>                                     Top1_mass_cm                   ;
            std::vector<double>                                     Top1_phi_cm                    ;
            std::vector<double>                                     Top2_pt_cm                     ;
            std::vector<double>                                     Top2_eta_cm                    ;
            std::vector<double>                                     Top2_mass_cm                   ;
            std::vector<double>                                     Top2_phi_cm                    ;

            std::vector<bool>                                       Baseline                       ;
            std::vector<bool>                                       Baseline_blind                 ;
            std::vector<double>                                     weight                         ;
            std::vector<double>                                     weight_fsrUp                   ;
            std::vector<double>                                     weight_fsrDown                 ;
            std::vector<double>                                     weight_isrUp                   ;
            std::vector<double>                                     weight_isrDown                 ;            
            std::vector<double>                                     weight_sclUp                   ;
            std::vector<double>                                     weight_sclDown                 ;
            std::vector<double>                                     weight_pdfUp                   ;
            std::vector<double>                                     weight_pdfDown                 ;            
            std::vector<double>                                     weight_prfUp                   ;
            std::vector<double>                                     weight_prfDown                 ;
            std::vector<double>                                     weight_btgUp                   ;
            std::vector<double>                                     weight_btgDown                 ;            
            std::vector<double>                                     weight_jetUp                   ;
            std::vector<double>                                     weight_jetDown                 ;
            std::vector<double>                                     weight_lepUp                   ;
            std::vector<double>                                     weight_lepDown                 ;            
            std::vector<double>                                     weight_puUp                   ;
            std::vector<double>                                     weight_puDown                 ;            
            std::vector<std::vector<bool> >                         GoodJets                       ;
            std::vector<int>                                        NGoodJets                      ;
            std::vector<int>                                        NGoodBJets                     ;
            std::vector<int>                                        NTops                          ;
            std::vector<int>                                        NRTops                         ;
            std::vector<int>                                        NMTops                         ;
            std::vector<double>                                     HT_pt30                        ;
            std::vector<double>                                     Mbl                            ;
            std::vector<double>                                     Mll                            ;
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

            std::vector<double>                                     Lepton1_pt_cm                  ;
            std::vector<double>                                     Lepton1_eta_cm                 ;
            std::vector<double>                                     Lepton1_phi_cm                 ;
            std::vector<double>                                     Lepton2_pt_cm                  ;
            std::vector<double>                                     Lepton2_eta_cm                 ;
            std::vector<double>                                     Lepton2_phi_cm                 ;

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

                Jets.push_back(tr.getVec<utility::LorentzVector>("Jets"                       + jecvar));
                Jets_cm_top6.push_back(tr.getVec<utility::LorentzVector>("Jets_cm_top6"       + jecvar));
                Tops.push_back(tr.getVec<utility::LorentzVector>("topsLV"                     + jecvar));
                Top1_pt_cm.push_back(tr.getVar<double>("top1_pt_cm"                           + jecvar));
                Top1_eta_cm.push_back(tr.getVar<double>("top1_eta_cm"                         + jecvar));
                Top1_phi_cm.push_back(tr.getVar<double>("top1_phi_cm"                         + jecvar));
                Top1_mass_cm.push_back(tr.getVar<double>("top1_mass_cm"                       + jecvar));
                Top2_pt_cm.push_back(tr.getVar<double>("top2_pt_cm"                           + jecvar));
                Top2_eta_cm.push_back(tr.getVar<double>("top2_eta_cm"                         + jecvar));
                Top2_phi_cm.push_back(tr.getVar<double>("top2_phi_cm"                         + jecvar));
                Top2_mass_cm.push_back(tr.getVar<double>("top2_mass_cm"                       + jecvar));

                GoodJets.push_back(tr.getVec<bool>("GoodJets_pt30"                                 + jecvar));
                NGoodJets.push_back(tr.getVar<int>("NGoodJets_pt30"                                + jecvar));
                NGoodBJets.push_back(tr.getVar<int>("NGoodBJets_pt30"                              + jecvar));
                NTops.push_back(tr.getVar<int>("ntops"                                             + jecvar));
                NRTops.push_back(tr.getVar<int>("ntops_3jet"                                       + jecvar));
                NMTops.push_back(tr.getVar<int>("ntops_1jet"                                       + jecvar));
                HT_pt30.push_back(tr.getVar<double>("HT_trigger_pt30"                              + jecvar));
                Baseline.push_back(tr.getVar<bool>("passBaseline" + channel + "l_Good"             + jecvar));
                Baseline_blind.push_back(tr.getVar<bool>("passBaseline" + channel + "l_Good_blind" + jecvar));
                Mbl.push_back(tr.getVar<double>("Mbl"                                              + jecvar));
                Mll.push_back(tr.getVar<double>("mll"                                              + jecvar));
                dRbjets.push_back(tr.getVar<double>("dR_bjets"                                     + jecvar));

                DoubleDisCo_massReg.push_back(tr.getVar<double>("DoubleDisCo_massReg_" + channel + "l_RPV" + jecvar));
                DoubleDisCo_disc1.push_back(tr.getVar<double>("DoubleDisCo_disc1_"     + channel + "l_RPV" + jecvar));
                DoubleDisCo_disc2.push_back(tr.getVar<double>("DoubleDisCo_disc2_"     + channel + "l_RPV" + jecvar));

                fwm2_top6.push_back(tr.getVar<double>("fwm2_top6"          + jecvar));
                fwm3_top6.push_back(tr.getVar<double>("fwm3_top6"          + jecvar));
                fwm4_top6.push_back(tr.getVar<double>("fwm4_top6"          + jecvar));
                fwm5_top6.push_back(tr.getVar<double>("fwm5_top6"          + jecvar));
                jmt_ev0_top6.push_back(tr.getVar<double>("jmt_ev0_top6"    + jecvar));
                jmt_ev1_top6.push_back(tr.getVar<double>("jmt_ev1_top6"    + jecvar));
                jmt_ev2_top6.push_back(tr.getVar<double>("jmt_ev2_top6"    + jecvar));
                nMVAJets.push_back(tr.getVar<unsigned int>("nMVAJets"      + jecvar));
                nMVALeptons.push_back(tr.getVar<unsigned int>("nMVALeptons"      + jecvar));

                combinedNthJetPt.push_back(tr.getVar<double>("combined"    + std::to_string(nMVAJets.back()) + "thToLastJet" + flavName + "_pt_cm" + jecvar));
                combinedNthJetPtrHT.push_back(tr.getVar<double>("combined" + std::to_string(nMVAJets.back()) + "thToLastJet" + flavName + "_ptrHT_cm" + jecvar));
                combinedNthJetEta.push_back(tr.getVar<double>("combined"   + std::to_string(nMVAJets.back()) + "thToLastJet" + flavName + "_eta_cm" + jecvar));
                combinedNthJetPhi.push_back(tr.getVar<double>("combined"   + std::to_string(nMVAJets.back()) + "thToLastJet" + flavName + "_phi_cm" + jecvar));
                combinedNthJetM.push_back(tr.getVar<double>("combined"     + std::to_string(nMVAJets.back()) + "thToLastJet" + flavName + "_m_cm" + jecvar));
                combinedNthJetE.push_back(tr.getVar<double>("combined"     + std::to_string(nMVAJets.back()) + "thToLastJet" + flavName + "_E_cm" + jecvar));

                Stop1_pt_cm_OldSeed.push_back(tr.getVar<double>("Stop1_pt_cm_OldSeed"       + jecvar));
                Stop1_eta_cm_OldSeed.push_back(tr.getVar<double>("Stop1_eta_cm_OldSeed"     + jecvar));
                Stop1_phi_cm_OldSeed.push_back(tr.getVar<double>("Stop1_phi_cm_OldSeed"     + jecvar));
                Stop1_mass_cm_OldSeed.push_back(tr.getVar<double>("Stop1_mass_cm_OldSeed"   + jecvar));
                Stop2_pt_cm_OldSeed.push_back(tr.getVar<double>("Stop2_pt_cm_OldSeed"       + jecvar));
                Stop2_eta_cm_OldSeed.push_back(tr.getVar<double>("Stop2_eta_cm_OldSeed"     + jecvar));
                Stop2_phi_cm_OldSeed.push_back(tr.getVar<double>("Stop2_phi_cm_OldSeed"     + jecvar));
                Stop2_mass_cm_OldSeed.push_back(tr.getVar<double>("Stop2_mass_cm_OldSeed"   + jecvar));
                Stop1_pt_cm_TopSeed.push_back(tr.getVar<double>("Stop1_pt_cm_TopSeed"       + jecvar));
                Stop1_eta_cm_TopSeed.push_back(tr.getVar<double>("Stop1_eta_cm_TopSeed"     + jecvar));
                Stop1_phi_cm_TopSeed.push_back(tr.getVar<double>("Stop1_phi_cm_TopSeed"     + jecvar));
                Stop1_mass_cm_TopSeed.push_back(tr.getVar<double>("Stop1_mass_cm_TopSeed"   + jecvar));
                Stop2_pt_cm_TopSeed.push_back(tr.getVar<double>("Stop2_pt_cm_TopSeed"       + jecvar));
                Stop2_eta_cm_TopSeed.push_back(tr.getVar<double>("Stop2_eta_cm_TopSeed"     + jecvar));
                Stop2_phi_cm_TopSeed.push_back(tr.getVar<double>("Stop2_phi_cm_TopSeed"     + jecvar));
                Stop2_mass_cm_TopSeed.push_back(tr.getVar<double>("Stop2_mass_cm_TopSeed"   + jecvar));

                Lepton1_pt_cm.push_back(tr.getVar<double>("GoodLeptons_pt_1"       + jecvar));
                Lepton1_eta_cm.push_back(tr.getVar<double>("GoodLeptons_eta_1"     + jecvar));
                Lepton1_phi_cm.push_back(tr.getVar<double>("GoodLeptons_phi_1"     + jecvar));
                Lepton2_pt_cm.push_back(tr.getVar<double>("GoodLeptons_pt_2"       + jecvar));
                Lepton2_eta_cm.push_back(tr.getVar<double>("GoodLeptons_eta_2"     + jecvar));
                Lepton2_phi_cm.push_back(tr.getVar<double>("GoodLeptons_phi_2"     + jecvar));

                GoodLeptons.push_back(tr.getVec<std::pair<std::string, utility::LorentzVector>>("GoodLeptons" + jecvar));
                GoodLeptonsCharge.push_back(tr.getVec<int>("GoodLeptonsCharge"                                + jecvar));
                GoodLeptonsMiniIso.push_back(tr.getVec<double>("GoodLeptonsMiniIso"                           + jecvar));

                Stop1_mass_PtRank_matched.push_back(  runtype == "Data" ? -999.0 : tr.getVar<float>("stop1_ptrank_mass"+jecvar));
                Stop2_mass_PtRank_matched.push_back(  runtype == "Data" ? -999.0 : tr.getVar<float>("stop2_ptrank_mass"+jecvar));
                Stop1_mass_MassRank_matched.push_back(runtype == "Data" ? -999.0 : tr.getVar<float>("stop1_mrank_mass" +jecvar));
                Stop2_mass_MassRank_matched.push_back(runtype == "Data" ? -999.0 : tr.getVar<float>("stop2_mrank_mass" +jecvar));
                Stop_mass_average_matched.push_back(  runtype == "Data" ? -999.0 : tr.getVar<double>("stop_avemass"    +jecvar));

                // ---------------------------------------------------------
                // Here assume number cm jets is same in CR and SR selection
                // ---------------------------------------------------------
                std::vector<double> tempJets_cm_flavb;  
                std::vector<double> tempJets_cm_flavc;  
                std::vector<double> tempJets_cm_flavg;  
                std::vector<double> tempJets_cm_flavuds;
                std::vector<double> tempJets_cm_flavq;  
                
                for (unsigned int iJet = 1; iJet <= nMVAJets.back(); iJet++) {
                    tempJets_cm_flavb.push_back(tr.getVar<double>("Jet_flavb_"     + std::to_string(iJet) + jecvar));
                    tempJets_cm_flavc.push_back(tr.getVar<double>("Jet_flavc_"     + std::to_string(iJet) + jecvar));
                    tempJets_cm_flavg.push_back(tr.getVar<double>("Jet_flavg_"     + std::to_string(iJet) + jecvar));
                    tempJets_cm_flavuds.push_back(tr.getVar<double>("Jet_flavuds_" + std::to_string(iJet) + jecvar));
                    tempJets_cm_flavq.push_back(tr.getVar<double>("Jet_flavq_"     + std::to_string(iJet) + jecvar));
                }

                Jets_cm_flavb.push_back(tempJets_cm_flavb);
                Jets_cm_flavc.push_back(tempJets_cm_flavc);
                Jets_cm_flavg.push_back(tempJets_cm_flavg);
                Jets_cm_flavuds.push_back(tempJets_cm_flavuds);
                Jets_cm_flavq.push_back(tempJets_cm_flavq);

                Jets_flavb.push_back(tr.getVec<float>("Jets" + jecvar + "_bJetTagDeepFlavourtotb"));
                Jets_flavc.push_back(tr.getVec<float>("Jets" + jecvar + "_bJetTagDeepFlavourprobc"));
                Jets_flavg.push_back(tr.getVec<float>("Jets" + jecvar + "_bJetTagDeepFlavourprobg"));
                Jets_flavuds.push_back(tr.getVec<float>("Jets" + jecvar + "_bJetTagDeepFlavourprobuds"));
                Jets_flavq.push_back(tr.getVec<float>("Jets" + jecvar + "_bJetTagDeepFlavourtotq"));

                std::map<std::string, std::vector<bool> > tempRegionMap;
                std::map<std::string, std::vector<bool> > tempRegionMap_CR;

                regions.push_back(tr.getVec<std::string>("regions_" + channel + "l_RPV"));
                for (const std::string& region : regions.back())
                {
                    tempRegionMap[region]    = tr.getVec<bool>("DoubleDisCo_" + region + "_"            + channel + "l_RPV" + jecvar); 
                    tempRegionMap_CR[region] = tr.getVec<bool>("DoubleDisCo_" + region + "_NonIsoMuon_" + channel + "l_RPV" + jecvar); 
                }

                // --------------------------------------
                // Now fill up the CR vector of variables
                // --------------------------------------
                DoubleDisCo_passRegions.push_back(tempRegionMap); 
                DoubleDisCo_passRegions_CR.push_back(tempRegionMap_CR); 

                Jets_CR.push_back(tr.getVec<utility::LorentzVector>("Jets"         + jecvar));
                Jets_cm_top6_CR.push_back(tr.getVec<utility::LorentzVector>(mvaName + "Jets_cm_top6" + jecvar));
                Tops_CR.push_back(tr.getVec<utility::LorentzVector>("topsLV"                         + jecvar));
                Top1_pt_cm_CR.push_back(tr.getVar<double>("top1_pt_cm"                               + jecvar));
                Top1_eta_cm_CR.push_back(tr.getVar<double>("top1_eta_cm"                             + jecvar));
                Top1_phi_cm_CR.push_back(tr.getVar<double>("top1_phi_cm"                             + jecvar));
                Top1_mass_cm_CR.push_back(tr.getVar<double>("top1_mass_cm"                           + jecvar));
                Top2_pt_cm_CR.push_back(tr.getVar<double>("top2_pt_cm"                               + jecvar));
                Top2_eta_cm_CR.push_back(tr.getVar<double>("top2_eta_cm"                             + jecvar));
                Top2_phi_cm_CR.push_back(tr.getVar<double>("top2_phi_cm"                             + jecvar));
                Top2_mass_cm_CR.push_back(tr.getVar<double>("top2_mass_cm"                           + jecvar));

                GoodJets_CR.push_back(tr.getVec<bool>(jetsName + "Jets_pt30"        + jecvar));
                NGoodJets_CR.push_back(tr.getVar<int>("N" + jetsName + "Jets_pt30" + jecvar));
                NGoodBJets_CR.push_back(tr.getVar<int>("NGoodBJets_pt30" + jecvar));
                NTops_CR.push_back(tr.getVar<int>("ntops" + jecvar));
                NRTops_CR.push_back(tr.getVar<int>("ntops_3jet" + jecvar));
                NMTops_CR.push_back(tr.getVar<int>("ntops_1jet" + jecvar));

                HT_pt30_CR.push_back(tr.getVar<double>("HT_" + htName + "_pt30"    + jecvar));
                Baseline_CR.push_back(tr.getVar<bool>("pass_qcdCR"                 + jecvar));
                Baseline_CR.push_back(tr.getVar<bool>("pass_qcdCR_1b"              + jecvar));
                Baseline_CR.push_back(tr.getVar<bool>("pass_qcdCR_1t"              + jecvar));
                Baseline_CR.push_back(tr.getVar<bool>("pass_qcdCR_1b_1t"           + jecvar));
                Baseline_CR.push_back(tr.getVar<bool>("pass_qcdCR_2b"              + jecvar));
                Baseline_CR.push_back(tr.getVar<bool>("pass_qcdCR_45"              + jecvar));
                Baseline_CR.push_back(tr.getVar<bool>("pass_qcdCR_45_1b"           + jecvar));
                Baseline_CR.push_back(tr.getVar<bool>("pass_qcdCR_45_2b"           + jecvar));
                Baseline_CR.push_back(tr.getVar<bool>("pass_qcdCR_45_1t"           + jecvar));
                Baseline_CR.push_back(tr.getVar<bool>("pass_qcdCR_45_2t"           + jecvar));
                Baseline_CR.push_back(tr.getVar<bool>("pass_qcdCR_45_1b_1t"        + jecvar));
                Baseline_CR.push_back(tr.getVar<bool>("pass_qcdCR_45_2b_1t"        + jecvar));
                Baseline_CR.push_back(tr.getVar<bool>("pass_qcdCR_all"             + jecvar));

                Mbl_CR.push_back(tr.getVar<double>("Mbl"                           + jecvar));
                Mll_CR.push_back(tr.getVar<double>("mll"                           + jecvar));
                dRbjets_CR.push_back(tr.getVar<double>("dR_bjets"                  + jecvar));

                DoubleDisCo_massReg_CR.push_back(tr.getVar<double>("DoubleDisCo_massReg_NonIsoMuon_" + channel + "l_RPV"  + jecvar));
                DoubleDisCo_disc1_CR.push_back(tr.getVar<double>("DoubleDisCo_disc1_NonIsoMuon_" + channel + "l_RPV"      + jecvar));
                DoubleDisCo_disc2_CR.push_back(tr.getVar<double>("DoubleDisCo_disc2_NonIsoMuon_" + channel + "l_RPV"      + jecvar));

                fwm2_top6_CR.push_back(tr.getVar<double>(mvaName + "fwm2_top6"         + jecvar));
                fwm3_top6_CR.push_back(tr.getVar<double>(mvaName + "fwm3_top6"         + jecvar));
                fwm4_top6_CR.push_back(tr.getVar<double>(mvaName + "fwm4_top6"         + jecvar));
                fwm5_top6_CR.push_back(tr.getVar<double>(mvaName + "fwm5_top6"         + jecvar));
                jmt_ev0_top6_CR.push_back(tr.getVar<double>(mvaName + "jmt_ev0_top6"   + jecvar));
                jmt_ev1_top6_CR.push_back(tr.getVar<double>(mvaName + "jmt_ev1_top6"   + jecvar));
                jmt_ev2_top6_CR.push_back(tr.getVar<double>(mvaName + "jmt_ev2_top6"   + jecvar));
                nMVAJets_CR.push_back(tr.getVar<unsigned int>(mvaName + "nMVAJets"                  + jecvar));
                nMVALeptons_CR.push_back(tr.getVar<unsigned int>(mvaName + "nMVALeptons"            + jecvar));

                combinedNthJetPt_CR.push_back(tr.getVar<double>("combined"    + std::to_string(nMVAJets_CR.back()) + "thToLastJet" + flavName + "_pt_cm" + jecvar));
                combinedNthJetPtrHT_CR.push_back(tr.getVar<double>("combined" + std::to_string(nMVAJets_CR.back()) + "thToLastJet" + flavName + "_ptrHT_cm" + jecvar));
                combinedNthJetEta_CR.push_back(tr.getVar<double>("combined"   + std::to_string(nMVAJets_CR.back()) + "thToLastJet" + flavName + "_eta_cm" + jecvar));
                combinedNthJetPhi_CR.push_back(tr.getVar<double>("combined"   + std::to_string(nMVAJets_CR.back()) + "thToLastJet" + flavName + "_phi_cm" + jecvar));
                combinedNthJetM_CR.push_back(tr.getVar<double>("combined"     + std::to_string(nMVAJets_CR.back()) + "thToLastJet" + flavName + "_m_cm" + jecvar));
                combinedNthJetE_CR.push_back(tr.getVar<double>("combined"     + std::to_string(nMVAJets_CR.back()) + "thToLastJet" + flavName + "_E_cm" + jecvar));

                Stop1_pt_cm_OldSeed_CR.push_back(tr.getVar<double>("Stop1_pt_cm_OldSeed_NonIsoMuon"       + jecvar));
                Stop1_eta_cm_OldSeed_CR.push_back(tr.getVar<double>("Stop1_eta_cm_OldSeed_NonIsoMuon"     + jecvar));
                Stop1_phi_cm_OldSeed_CR.push_back(tr.getVar<double>("Stop1_phi_cm_OldSeed_NonIsoMuon"     + jecvar));
                Stop1_mass_cm_OldSeed_CR.push_back(tr.getVar<double>("Stop1_mass_cm_OldSeed_NonIsoMuon"   + jecvar));
                Stop2_pt_cm_OldSeed_CR.push_back(tr.getVar<double>("Stop2_pt_cm_OldSeed_NonIsoMuon"       + jecvar));
                Stop2_eta_cm_OldSeed_CR.push_back(tr.getVar<double>("Stop2_eta_cm_OldSeed_NonIsoMuon"     + jecvar));
                Stop2_phi_cm_OldSeed_CR.push_back(tr.getVar<double>("Stop2_phi_cm_OldSeed_NonIsoMuon"     + jecvar));
                Stop2_mass_cm_OldSeed_CR.push_back(tr.getVar<double>("Stop2_mass_cm_OldSeed_NonIsoMuon"   + jecvar));
                Stop1_pt_cm_TopSeed_CR.push_back(tr.getVar<double>("Stop1_pt_cm_TopSeed_NonIsoMuon"       + jecvar));
                Stop1_eta_cm_TopSeed_CR.push_back(tr.getVar<double>("Stop1_eta_cm_TopSeed_NonIsoMuon"     + jecvar));
                Stop1_phi_cm_TopSeed_CR.push_back(tr.getVar<double>("Stop1_phi_cm_TopSeed_NonIsoMuon"     + jecvar));
                Stop1_mass_cm_TopSeed_CR.push_back(tr.getVar<double>("Stop1_mass_cm_TopSeed_NonIsoMuon"   + jecvar));
                Stop2_pt_cm_TopSeed_CR.push_back(tr.getVar<double>("Stop2_pt_cm_TopSeed_NonIsoMuon"       + jecvar));
                Stop2_eta_cm_TopSeed_CR.push_back(tr.getVar<double>("Stop2_eta_cm_TopSeed_NonIsoMuon"     + jecvar));
                Stop2_phi_cm_TopSeed_CR.push_back(tr.getVar<double>("Stop2_phi_cm_TopSeed_NonIsoMuon"     + jecvar));
                Stop2_mass_cm_TopSeed_CR.push_back(tr.getVar<double>("Stop2_mass_cm_TopSeed_NonIsoMuon"   + jecvar));

                Lepton1_pt_cm_CR.push_back(tr.getVar<double>("GoodNonIsoMuons_pt_1"       + jecvar));
                Lepton1_eta_cm_CR.push_back(tr.getVar<double>("GoodNonIsoMuons_eta_1"     + jecvar));
                Lepton1_phi_cm_CR.push_back(tr.getVar<double>("GoodNonIsoMuons_phi_1"     + jecvar));
                Lepton2_pt_cm_CR.push_back(tr.getVar<double>("GoodNonIsoMuons_pt_2"       + jecvar));
                Lepton2_eta_cm_CR.push_back(tr.getVar<double>("GoodNonIsoMuons_eta_2"     + jecvar));
                Lepton2_phi_cm_CR.push_back(tr.getVar<double>("GoodNonIsoMuons_phi_2"     + jecvar));

                GoodLeptons_CR.push_back(tr.getVec<std::pair<std::string, utility::LorentzVector>>("GoodNonIsoMuons" + jecvar));
                GoodLeptonsCharge_CR.push_back(tr.getVec<int>("GoodNonIsoMuonsCharge"                                + jecvar));
                GoodLeptonsMiniIso_CR.push_back(tr.getVec<double>("GoodNonIsoMuonsMiniIso"                           + jecvar));
               
                Stop1_mass_PtRank_matched_CR.push_back(  runtype == "Data" ? -999.0 : tr.getVar<float>("stop1_ptrank_mass" + jecvar));
                Stop2_mass_PtRank_matched_CR.push_back(  runtype == "Data" ? -999.0 : tr.getVar<float>("stop2_ptrank_mass" + jecvar));
                Stop1_mass_MassRank_matched_CR.push_back(runtype == "Data" ? -999.0 : tr.getVar<float>("stop1_mrank_mass"  + jecvar));
                Stop2_mass_MassRank_matched_CR.push_back(runtype == "Data" ? -999.0 : tr.getVar<float>("stop2_mrank_mass"  + jecvar));
                Stop_mass_average_matched_CR.push_back(  runtype == "Data" ? -999.0 : tr.getVar<double>("stop_avemass"     + jecvar));

                std::vector<double> tempJets_cm_flavb_CR;
                std::vector<double> tempJets_cm_flavc_CR;
                std::vector<double> tempJets_cm_flavg_CR;
                std::vector<double> tempJets_cm_flavuds_CR;
                std::vector<double> tempJets_cm_flavq_CR;

                for (unsigned int iJet = 1; iJet <= nMVAJets_CR.back(); iJet++) {
                    tempJets_cm_flavb_CR.push_back(  tr.getVar<double>("Jet" + flavName + "_flavb_"  +std::to_string(iJet)+jecvar));
                    tempJets_cm_flavc_CR.push_back(  tr.getVar<double>("Jet" + flavName + "_flavc_"  +std::to_string(iJet)+jecvar));
                    tempJets_cm_flavg_CR.push_back(  tr.getVar<double>("Jet" + flavName + "_flavg_"  +std::to_string(iJet)+jecvar));
                    tempJets_cm_flavuds_CR.push_back(tr.getVar<double>("Jet" + flavName + "_flavuds_"+std::to_string(iJet)+jecvar));
                    tempJets_cm_flavq_CR.push_back(  tr.getVar<double>("Jet" + flavName + "_flavq_"  +std::to_string(iJet)+jecvar));
                }

                Jets_cm_flavb_CR.push_back(tempJets_cm_flavb_CR);
                Jets_cm_flavc_CR.push_back(tempJets_cm_flavc_CR);
                Jets_cm_flavg_CR.push_back(tempJets_cm_flavg_CR);
                Jets_cm_flavuds_CR.push_back(tempJets_cm_flavuds_CR);
                Jets_cm_flavq_CR.push_back(tempJets_cm_flavq_CR);

                Jets_flavb_CR.push_back(tr.getVec<float>("Jets" + jecvar + "_bJetTagDeepFlavourtotb"));
                Jets_flavc_CR.push_back(tr.getVec<float>("Jets" + jecvar + "_bJetTagDeepFlavourprobc"));
                Jets_flavg_CR.push_back(tr.getVec<float>("Jets" + jecvar + "_bJetTagDeepFlavourprobg"));
                Jets_flavuds_CR.push_back(tr.getVec<float>("Jets" + jecvar + "_bJetTagDeepFlavourprobuds"));
                Jets_flavq_CR.push_back(tr.getVec<float>("Jets" + jecvar + "_bJetTagDeepFlavourtotq"));

                // ---------------------------------------------------------------------------
                // Calculate the event weights for systvars fsrUp/Down, isrUp/Down, etc
                //  -- fsr/isr are event weight based variations
                //  -- calculate event weights for them
                //  -- make histograms for systvars fsrUp/Down, isrUp/Down in root file
                // ----------------------------------------------------------------------------
                double theWeight       = 1.0;
                double theWeightQCDCR  = 1.0;
                double theWeight_fsrUp = 1.0, theWeight_fsrDown = 1.0;
                double theWeight_isrUp = 1.0, theWeight_isrDown = 1.0;
                double theWeight_sclUp = 1.0, theWeight_sclDown = 1.0;
                double theWeight_pdfUp = 1.0, theWeight_pdfDown = 1.0;
                double theWeight_prfUp = 1.0, theWeight_prfDown = 1.0;
                double theWeight_lepUp = 1.0, theWeight_lepDown = 1.0;
                double theWeight_jetUp = 1.0, theWeight_jetDown = 1.0;
                double theWeight_btgUp = 1.0, theWeight_btgDown = 1.0;
                double theWeight_puUp  = 1.0, theWeight_puDown  = 1.0;
                if(runtype == "MC" )
                {
                    theWeight         = tr.getVar<double>("TotalWeight_" + channel + "l"         + jecvar);
                    theWeightQCDCR    = tr.getVar<double>("TotalWeight_" + channel + "l_QCDCR"   + jecvar);
                    theWeight_fsrUp   = tr.getVar<double>("TotalWeight_" + channel + "l_FSRup"   + jecvar);
                    theWeight_fsrDown = tr.getVar<double>("TotalWeight_" + channel + "l_FSRdown" + jecvar);
                    theWeight_isrUp   = tr.getVar<double>("TotalWeight_" + channel + "l_ISRup"   + jecvar);
                    theWeight_isrDown = tr.getVar<double>("TotalWeight_" + channel + "l_ISRdown" + jecvar);
                    theWeight_sclUp   = tr.getVar<double>("TotalWeight_" + channel + "l_SclUp"   + jecvar);
                    theWeight_sclDown = tr.getVar<double>("TotalWeight_" + channel + "l_SclDown" + jecvar);
                    theWeight_pdfUp   = tr.getVar<double>("TotalWeight_" + channel + "l_PDFup"   + jecvar);
                    theWeight_pdfDown = tr.getVar<double>("TotalWeight_" + channel + "l_PDFdown" + jecvar);
                    theWeight_prfUp   = tr.getVar<double>("TotalWeight_" + channel + "l_PrfUp"   + jecvar);
                    theWeight_prfDown = tr.getVar<double>("TotalWeight_" + channel + "l_PrfDown" + jecvar);
                    theWeight_btgUp   = tr.getVar<double>("TotalWeight_" + channel + "l_BtgUp"   + jecvar);
                    theWeight_btgDown = tr.getVar<double>("TotalWeight_" + channel + "l_BtgDown" + jecvar);
                    theWeight_puUp    = tr.getVar<double>("TotalWeight_" + channel + "l_PUup"    + jecvar);
                    theWeight_puDown  = tr.getVar<double>("TotalWeight_" + channel + "l_PUdown"  + jecvar);

                    if (channel != "0")
                    {
                        theWeight_lepUp   = tr.getVar<double>("TotalWeight_" + channel + "l_LepUp"   + jecvar);
                        theWeight_lepDown = tr.getVar<double>("TotalWeight_" + channel + "l_LepDown" + jecvar);
                    }
                    else
                    {
                        theWeight_jetUp   = tr.getVar<double>("TotalWeight_" + channel + "l_JetUp"   + jecvar);
                        theWeight_jetDown = tr.getVar<double>("TotalWeight_" + channel + "l_JetDown" + jecvar);
                    }
                }
                
                // baseline
                weight.push_back(theWeight);
                weight_fsrUp.push_back(theWeight_fsrUp);
                weight_fsrDown.push_back(theWeight_fsrDown);
                weight_isrUp.push_back(theWeight_isrUp);
                weight_isrDown.push_back(theWeight_isrDown);
                weight_sclUp.push_back(theWeight_sclUp);
                weight_sclDown.push_back(theWeight_sclDown);
                weight_pdfUp.push_back(theWeight_pdfUp);
                weight_pdfDown.push_back(theWeight_pdfDown);
                weight_prfUp.push_back(theWeight_prfUp);
                weight_prfDown.push_back(theWeight_prfDown);
                weight_lepUp.push_back(theWeight_lepUp);
                weight_lepDown.push_back(theWeight_lepDown);
                weight_jetUp.push_back(theWeight_jetUp);
                weight_jetDown.push_back(theWeight_jetDown);
                weight_btgUp.push_back(theWeight_btgUp);
                weight_btgDown.push_back(theWeight_btgDown);
                weight_puUp.push_back(theWeight_puUp);
                weight_puDown.push_back(theWeight_puDown);

                // qcd cr
                weight_CR.push_back(theWeightQCDCR);
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
                {"_0l_QCDCR_45"      , Baseline_CR[5]},
                {"_0l_QCDCR_45_1b"   , Baseline_CR[6]},
                {"_0l_QCDCR_45_2b"   , Baseline_CR[7]},
                {"_0l_QCDCR_45_1t"   , Baseline_CR[8]},
                {"_0l_QCDCR_45_2t"   , Baseline_CR[9]},
                {"_0l_QCDCR_45_1b_1t", Baseline_CR[10]},
                {"_0l_QCDCR_45_2b_1t", Baseline_CR[11]},
                {"_0l_QCDCR_all"     , Baseline_CR[12]},
                {"_1l_QCDCR"         , Baseline_CR[0]}, 
                {"_2l_QCDCR"         , Baseline_CR[0]}, 
            };

            // Put assume 7 jets and 2 leptons for making the histograms
            if(!initHistos)
            {
                Debug("Initializing the histograms", __LINE__);

                Preinit(7, 2);
                InitHistos(cut_map, regions, runtype);
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
                    {"11incl", !isQCD ? NGoodJets[channel]>=11 : NGoodJets_CR[channel]>=11},
                    {"12incl", !isQCD ? NGoodJets[channel]>=12 : NGoodJets_CR[channel]>=12},
                    {"13incl", !isQCD ? NGoodJets[channel]>=13 : NGoodJets_CR[channel]>=13},
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
                    // loop over the tt and systvars fsrUp/Down, isrUp/Down
                    // to make indivudual histograms with the label tt and systvars in TT root file
                    // --------------------------------------------------------------------------
                    for (const auto& ttvar : systvars)
                    {   
                        Debug("Top of ttvar loop for ttvar: " + ttvar, __LINE__);

                        std::string ttvarStr = "";
                        if (ttvar != "nom")
                            ttvarStr = "_" + ttvar; 

                        if (ttvar != "nom" and runtype == "Data")
                            continue;

                        if (ttvar != "nom" and jecvar != "")
                            continue;

                        if (kv.first.find("QCDCR_") != std::string::npos and ttvar != "nom")
                            continue;

                        double w = 1.0;

                        // get the event weights for fsrUp/Down, isrUp/Down, etc.
                        if (ttvar == "nom")
                        {
                            w  = !isQCD ? weight[channel] : weight_CR[channel];
                        }
                        else if (ttvar == "fsrUp")
                        {
                            w = weight_fsrUp[channel];
                        }
                        else if (ttvar == "fsrDown")
                        {
                            w = weight_fsrDown[channel];
                        }
                        else if (ttvar == "isrUp")
                        {
                            w =weight_isrUp[channel];
                        }
                        else if (ttvar == "isrDown")
                        {
                            w = weight_isrDown[channel];
                        }
                        else if (ttvar == "pdfUp")
                        {
                            w = weight_pdfUp[channel];
                        }
                        else if (ttvar == "pdfDown")
                        {
                            w = weight_pdfDown[channel];
                        }
                        else if (ttvar == "prfUp")
                        {
                            w =weight_prfUp[channel];
                        }
                        else if (ttvar == "prfDown")
                        {
                            w = weight_prfDown[channel];
                        }
                        else if (ttvar == "btgUp")
                        {
                            w = weight_btgUp[channel];
                        }
                        else if (ttvar == "btgDown")
                        {
                            w = weight_btgDown[channel];
                        }
                        else if (ttvar == "jetUp")
                        {
                            w =weight_jetUp[channel];
                        }
                        else if (ttvar == "jetDown")
                        {
                            w = weight_jetDown[channel];
                        }
                        else if (ttvar == "lepUp")
                        {
                            w =weight_lepUp[channel];
                        }
                        else if (ttvar == "lepDown")
                        {
                            w = weight_lepDown[channel];
                        }
                        else if (ttvar == "puUp")
                        {
                            w =weight_puUp[channel];
                        }
                        else if (ttvar == "puDown")
                        {
                            w = weight_puDown[channel];
                        }
                        else if (ttvar == "sclUp")
                        {
                            w =weight_sclUp[channel];
                        }
                        else if (ttvar == "sclDown")
                        {
                            w = weight_sclDown[channel];
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

                            name = kv.first + njetsStr + regionStr + ttvarStr + jecStr;

                            if (kv.second and inNjetsBin and inRegion)
                            {
                                if (njets == "Incl" and region == "ABCD" and jecvar == "" and ttvar == "nom")
                                {

                                    Debug("Filling event variable histograms with name: " + name, __LINE__);
                                    // ----------------------------------------------------------------------------------------------------
                                    // Filling event-level variables
                                    // ----------------------------------------------------------------------------------------------------
                                    my_histos["h_Njets"        + name]->Fill(!isQCD ?    NGoodJets[channel] : NGoodJets_CR[channel],    w);
                                    my_histos["h_Nbjets"       + name]->Fill(!isQCD ?   NGoodBJets[channel] : NGoodBJets_CR[channel],   w);
                                    my_histos["h_Ntops"        + name]->Fill(!isQCD ?        NTops[channel] : NTops_CR[channel],        w);
                                    my_histos["h_FWM2_top6"    + name]->Fill(!isQCD ?    fwm2_top6[channel] : fwm2_top6_CR[channel],    w);
                                    my_histos["h_FWM3_top6"    + name]->Fill(!isQCD ?    fwm3_top6[channel] : fwm3_top6_CR[channel],    w);
                                    my_histos["h_FWM4_top6"    + name]->Fill(!isQCD ?    fwm4_top6[channel] : fwm4_top6_CR[channel],    w);
                                    my_histos["h_FWM5_top6"    + name]->Fill(!isQCD ?    fwm5_top6[channel] : fwm5_top6_CR[channel],    w);
                                    my_histos["h_JMT_ev0_top6" + name]->Fill(!isQCD ? jmt_ev0_top6[channel] : jmt_ev0_top6_CR[channel], w);
                                    my_histos["h_JMT_ev1_top6" + name]->Fill(!isQCD ? jmt_ev1_top6[channel] : jmt_ev1_top6_CR[channel], w);
                                    my_histos["h_JMT_ev2_top6" + name]->Fill(!isQCD ? jmt_ev2_top6[channel] : jmt_ev2_top6_CR[channel], w);
                                    my_histos["h_Mbl"          + name]->Fill(!isQCD ?          Mbl[channel] : Mbl_CR[channel],          w);
                                    my_histos["h_Mll"          + name]->Fill(!isQCD ?          Mll[channel] : Mll_CR[channel],          w);
                                    my_histos["h_dRbjets"      + name]->Fill(!isQCD ?      dRbjets[channel] : dRbjets_CR[channel],      w);
                                    my_histos["h_HT"           + name]->Fill(!isQCD ?      HT_pt30[channel] : HT_pt30_CR[channel],      w);

                                    Debug("Filling stop variable histograms with name: " + name, __LINE__);

                                    // ----------------------------------------------------------------------------------------------------
                                    // Filling top 4-vector and stop 4-vector variables, sort objects by transverse momentum
                                    // ----------------------------------------------------------------------------------------------------
                                    unsigned int leadIndex    = 0;
                                    unsigned int subleadIndex = 1;
                                    if ( (!isQCD and Top2_pt_cm[channel]    > Top1_pt_cm[channel]) or \
                                         (isQCD  and Top2_pt_cm_CR[channel] > Top1_pt_cm_CR[channel]))
                                    {
                                        leadIndex    = 1;
                                        subleadIndex = 0;
                                    }
                                    my_histos["h_Top" + std::to_string(leadIndex+1) + "_Pt_cm"   + name]->Fill(!isQCD ? Top1_pt_cm[channel]   : Top1_pt_cm_CR[channel],   w);
                                    my_histos["h_Top" + std::to_string(leadIndex+1) + "_Eta_cm"  + name]->Fill(!isQCD ? Top1_eta_cm[channel]  : Top1_eta_cm_CR[channel],  w);
                                    my_histos["h_Top" + std::to_string(leadIndex+1) + "_Phi_cm"  + name]->Fill(!isQCD ? Top1_phi_cm[channel]  : Top1_phi_cm_CR[channel],  w);
                                    my_histos["h_Top" + std::to_string(leadIndex+1) + "_Mass_cm" + name]->Fill(!isQCD ? Top1_mass_cm[channel] : Top1_mass_cm_CR[channel], w);

                                    my_histos["h_Top" + std::to_string(subleadIndex+1) + "_Pt_cm"   + name]->Fill(!isQCD ? Top2_pt_cm[channel]   : Top2_pt_cm_CR[channel],   w);
                                    my_histos["h_Top" + std::to_string(subleadIndex+1) + "_Eta_cm"  + name]->Fill(!isQCD ? Top2_eta_cm[channel]  : Top2_eta_cm_CR[channel],  w);
                                    my_histos["h_Top" + std::to_string(subleadIndex+1) + "_Phi_cm"  + name]->Fill(!isQCD ? Top2_phi_cm[channel]  : Top2_phi_cm_CR[channel],  w);
                                    my_histos["h_Top" + std::to_string(subleadIndex+1) + "_Mass_cm" + name]->Fill(!isQCD ? Top2_mass_cm[channel] : Top2_mass_cm_CR[channel], w);

                                    if ((!isQCD and Tops[channel].size() >= 2) or (isQCD and Tops_CR[channel].size() >= 2))
                                    {
                                        leadIndex    = 0;
                                        subleadIndex = 1;
                                        if ( (!isQCD and Tops[channel].at(1).Pt()    > Tops[channel].at(0).Pt()) or \
                                             (isQCD  and Tops_CR[channel].at(1).Pt() > Tops_CR[channel].at(0).Pt()))
                                        {
                                            leadIndex    = 1;
                                            subleadIndex = 0;
                                        }
                                        my_histos["h_Top" + std::to_string(leadIndex+1) + "_Pt"   + name]->Fill(!isQCD ? Tops[channel].at(leadIndex).Pt()  : Tops_CR[channel].at(leadIndex).Pt(),  w);
                                        my_histos["h_Top" + std::to_string(leadIndex+1) + "_Eta"  + name]->Fill(!isQCD ? Tops[channel].at(leadIndex).Eta() : Tops_CR[channel].at(leadIndex).Eta(), w);
                                        my_histos["h_Top" + std::to_string(leadIndex+1) + "_Phi"  + name]->Fill(!isQCD ? Tops[channel].at(leadIndex).Phi() : Tops_CR[channel].at(leadIndex).Phi(), w);
                                        my_histos["h_Top" + std::to_string(leadIndex+1) + "_Mass" + name]->Fill(!isQCD ? Tops[channel].at(leadIndex).M()   : Tops_CR[channel].at(leadIndex).M(),   w);

                                        my_histos["h_Top" + std::to_string(subleadIndex+1) + "_Pt"   + name]->Fill(!isQCD ? Tops[channel].at(subleadIndex).Pt()  : Tops_CR[channel].at(subleadIndex).Pt(),  w);
                                        my_histos["h_Top" + std::to_string(subleadIndex+1) + "_Eta"  + name]->Fill(!isQCD ? Tops[channel].at(subleadIndex).Eta() : Tops_CR[channel].at(subleadIndex).Eta(), w);
                                        my_histos["h_Top" + std::to_string(subleadIndex+1) + "_Phi"  + name]->Fill(!isQCD ? Tops[channel].at(subleadIndex).Phi() : Tops_CR[channel].at(subleadIndex).Phi(), w);
                                        my_histos["h_Top" + std::to_string(subleadIndex+1) + "_Mass" + name]->Fill(!isQCD ? Tops[channel].at(subleadIndex).M()   : Tops_CR[channel].at(subleadIndex).M(),   w);
                                    }

                                    leadIndex    = 0;
                                    subleadIndex = 1;
                                    if ( (!isQCD and Stop2_pt_cm_OldSeed[channel]    > Stop1_pt_cm_OldSeed[channel]) or \
                                         (isQCD  and Stop2_pt_cm_OldSeed_CR[channel] > Stop1_pt_cm_OldSeed_CR[channel]))
                                    {
                                        leadIndex    = 1;
                                        subleadIndex = 0; 
                                    }
                                    my_histos["h_Stop" + std::to_string(leadIndex+1) + "_Pt_cm_OldSeed"   + name]->Fill(!isQCD ? Stop1_pt_cm_OldSeed[channel]   : Stop1_pt_cm_OldSeed_CR[channel],   w);
                                    my_histos["h_Stop" + std::to_string(leadIndex+1) + "_Eta_cm_OldSeed"  + name]->Fill(!isQCD ? Stop1_eta_cm_OldSeed[channel]  : Stop1_eta_cm_OldSeed_CR[channel],  w);
                                    my_histos["h_Stop" + std::to_string(leadIndex+1) + "_Phi_cm_OldSeed"  + name]->Fill(!isQCD ? Stop1_phi_cm_OldSeed[channel]  : Stop1_phi_cm_OldSeed_CR[channel],  w);
                                    my_histos["h_Stop" + std::to_string(leadIndex+1) + "_Mass_cm_OldSeed" + name]->Fill(!isQCD ? Stop1_mass_cm_OldSeed[channel] : Stop1_mass_cm_OldSeed_CR[channel], w);

                                    my_histos["h_Stop" + std::to_string(subleadIndex+1) + "_Pt_cm_OldSeed"   + name]->Fill(!isQCD ? Stop2_pt_cm_OldSeed[channel]   : Stop2_pt_cm_OldSeed_CR[channel],   w);
                                    my_histos["h_Stop" + std::to_string(subleadIndex+1) + "_Eta_cm_OldSeed"  + name]->Fill(!isQCD ? Stop2_eta_cm_OldSeed[channel]  : Stop2_eta_cm_OldSeed_CR[channel],  w);
                                    my_histos["h_Stop" + std::to_string(subleadIndex+1) + "_Phi_cm_OldSeed"  + name]->Fill(!isQCD ? Stop2_phi_cm_OldSeed[channel]  : Stop2_phi_cm_OldSeed_CR[channel],  w);
                                    my_histos["h_Stop" + std::to_string(subleadIndex+1) + "_Mass_cm_OldSeed" + name]->Fill(!isQCD ? Stop2_mass_cm_OldSeed[channel] : Stop2_mass_cm_OldSeed_CR[channel], w);

                                    leadIndex    = 0;
                                    subleadIndex = 1;
                                    if ( (!isQCD and Stop2_pt_cm_TopSeed[channel]    > Stop1_pt_cm_TopSeed[channel]) or \
                                         (isQCD  and Stop2_pt_cm_TopSeed_CR[channel] > Stop1_pt_cm_TopSeed_CR[channel]))
                                    {
                                        leadIndex    = 1;
                                        subleadIndex = 0; 
                                    }
                                    my_histos["h_Stop" + std::to_string(leadIndex+1) + "_Pt_cm_TopSeed"   + name]->Fill(!isQCD ? Stop1_pt_cm_TopSeed[channel]   : Stop1_pt_cm_TopSeed_CR[channel],   w);
                                    my_histos["h_Stop" + std::to_string(leadIndex+1) + "_Eta_cm_TopSeed"  + name]->Fill(!isQCD ? Stop1_eta_cm_TopSeed[channel]  : Stop1_eta_cm_TopSeed_CR[channel],  w);
                                    my_histos["h_Stop" + std::to_string(leadIndex+1) + "_Phi_cm_TopSeed"  + name]->Fill(!isQCD ? Stop1_phi_cm_TopSeed[channel]  : Stop1_phi_cm_TopSeed_CR[channel],  w);
                                    my_histos["h_Stop" + std::to_string(leadIndex+1) + "_Mass_cm_TopSeed" + name]->Fill(!isQCD ? Stop1_mass_cm_TopSeed[channel] : Stop1_mass_cm_TopSeed_CR[channel], w);

                                    my_histos["h_Stop" + std::to_string(subleadIndex+1) + "_Pt_cm_TopSeed"   + name]->Fill(!isQCD ? Stop2_pt_cm_TopSeed[channel]   : Stop2_pt_cm_TopSeed_CR[channel],   w);
                                    my_histos["h_Stop" + std::to_string(subleadIndex+1) + "_Eta_cm_TopSeed"  + name]->Fill(!isQCD ? Stop2_eta_cm_TopSeed[channel]  : Stop2_eta_cm_TopSeed_CR[channel],  w);
                                    my_histos["h_Stop" + std::to_string(subleadIndex+1) + "_Phi_cm_TopSeed"  + name]->Fill(!isQCD ? Stop2_phi_cm_TopSeed[channel]  : Stop2_phi_cm_TopSeed_CR[channel],  w);
                                    my_histos["h_Stop" + std::to_string(subleadIndex+1) + "_Mass_cm_TopSeed" + name]->Fill(!isQCD ? Stop2_mass_cm_TopSeed[channel] : Stop2_mass_cm_TopSeed_CR[channel], w);

                                    my_histos["h_Stop1_mass_PtRank_matched"   + name]->Fill(!isQCD ? Stop1_mass_PtRank_matched[channel]   : Stop1_mass_PtRank_matched_CR[channel],   w);
                                    my_histos["h_Stop2_mass_PtRank_matched"   + name]->Fill(!isQCD ? Stop2_mass_PtRank_matched[channel]   : Stop2_mass_PtRank_matched_CR[channel],   w);
                                    my_histos["h_Stop1_mass_MassRank_matched" + name]->Fill(!isQCD ? Stop1_mass_MassRank_matched[channel] : Stop1_mass_MassRank_matched_CR[channel], w);
                                    my_histos["h_Stop2_mass_MassRank_matched" + name]->Fill(!isQCD ? Stop2_mass_MassRank_matched[channel] : Stop2_mass_MassRank_matched_CR[channel], w);
                                    my_histos["h_Stop_mass_average_matched"   + name]->Fill(!isQCD ? Stop_mass_average_matched[channel]   : Stop_mass_average_matched_CR[channel],   w);

                                    // -------------------------------------------------------------------------------------------------
                                    // Fill per jet variables including 4-vector and flavor information
                                    // -------------------------------------------------------------------------------------------------
                                    auto& cmJets     = Jets_cm_top6[channel];
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
                                        my_histos["h_Jet" + std::to_string(i) + "_cm_PtrHT"   + name]->Fill(cmJets.at(i-1).Pt()/ht, w);
                                        my_histos["h_Jet" + std::to_string(i) + "_cm_Pt"      + name]->Fill(cmJets.at(i-1).Pt(),    w);
                                        my_histos["h_Jet" + std::to_string(i) + "_cm_Eta"     + name]->Fill(cmJets.at(i-1).Eta(),   w);
                                        my_histos["h_Jet" + std::to_string(i) + "_cm_Phi"     + name]->Fill(cmJets.at(i-1).Phi(),   w);
                                        my_histos["h_Jet" + std::to_string(i) + "_cm_Mass"    + name]->Fill(cmJets.at(i-1).M(),     w);
                                        my_histos["h_Jet" + std::to_string(i) + "_cm_Energy"  + name]->Fill(cmJets.at(i-1).E(),     w);

                                        my_histos["h_Jet" + std::to_string(i) + "_cm_Flavb"   + name]->Fill(jetFlavb.at(i-1),   w);
                                        my_histos["h_Jet" + std::to_string(i) + "_cm_Flavc"   + name]->Fill(jetFlavc.at(i-1),   w);
                                        my_histos["h_Jet" + std::to_string(i) + "_cm_Flavg"   + name]->Fill(jetFlavg.at(i-1),   w);
                                        my_histos["h_Jet" + std::to_string(i) + "_cm_Flavq"   + name]->Fill(jetFlavq.at(i-1),   w);
                                        my_histos["h_Jet" + std::to_string(i) + "_cm_Flavuds" + name]->Fill(jetFlavuds.at(i-1), w);
                                    }

                                    Debug("Filling lepton variable histograms with name: " + name, __LINE__);

                                    // -------------------------------------------------------------------------------------------------
                                    // Fill per lepton (cm) variables including 4-vector and flavor information
                                    // -------------------------------------------------------------------------------------------------
                                    leadIndex    = 0;
                                    subleadIndex = 1;
                                    if ( (!isQCD and Lepton2_pt_cm[channel]    > Lepton1_pt_cm[channel]) or \
                                         (isQCD  and Lepton2_pt_cm_CR[channel] > Lepton1_pt_cm_CR[channel]))
                                    {
                                        leadIndex    = 1;
                                        subleadIndex = 0; 
                                    }

                                    my_histos["h_Lepton" + std::to_string(leadIndex+1) + "_cm_Pt"     + name]->Fill(!isQCD ? Lepton1_pt_cm[channel]  : Lepton1_pt_cm_CR[channel], w);
                                    my_histos["h_Lepton" + std::to_string(leadIndex+1) + "_cm_Phi"    + name]->Fill(!isQCD ? Lepton1_phi_cm[channel] : Lepton1_pt_cm_CR[channel], w);
                                    my_histos["h_Lepton" + std::to_string(leadIndex+1) + "_cm_Eta"    + name]->Fill(!isQCD ? Lepton1_eta_cm[channel] : Lepton1_pt_cm_CR[channel], w);

                                    my_histos["h_Lepton" + std::to_string(subleadIndex+1) + "_cm_Pt"  + name]->Fill(!isQCD ? Lepton2_pt_cm[channel]  : Lepton2_pt_cm_CR[channel], w);
                                    my_histos["h_Lepton" + std::to_string(subleadIndex+1) + "_cm_Phi" + name]->Fill(!isQCD ? Lepton2_phi_cm[channel] : Lepton2_pt_cm_CR[channel], w);
                                    my_histos["h_Lepton" + std::to_string(subleadIndex+1) + "_cm_Eta" + name]->Fill(!isQCD ? Lepton2_eta_cm[channel] : Lepton2_pt_cm_CR[channel], w);

                                    auto& theGoodLeptons  = GoodLeptons[channel];
                                    auto& theGoodLeptonsC = GoodLeptonsCharge[channel];
                                    auto& theGoodLeptonsI = GoodLeptonsMiniIso[channel];
                                    if (isQCD)
                                    {
                                       theGoodLeptons  = GoodLeptons_CR[channel]; 
                                       theGoodLeptonsC = GoodLeptonsCharge_CR[channel]; 
                                       theGoodLeptonsI = GoodLeptonsMiniIso_CR[channel]; 
                                    }

                                    if (theGoodLeptons.size() == 1)
                                    {
                                        std::string type1            = theGoodLeptons[0].first;
                                        utility::LorentzVector lvec1 = theGoodLeptons[0].second; 
                                        double charge1               = theGoodLeptonsC[0];
                                        double iso1                  = theGoodLeptonsI[0];

                                        my_histos["h_Lepton1_Pt"      + name]->Fill(lvec1.Pt(),  w);
                                        my_histos["h_Lepton1_Phi"     + name]->Fill(lvec1.Phi(), w);
                                        my_histos["h_Lepton1_Eta"     + name]->Fill(lvec1.Eta(), w);
                                        my_histos["h_Lepton1_Charge"  + name]->Fill(charge1,     w);
                                        my_histos["h_Lepton1_MiniIso" + name]->Fill(iso1,        w);

                                        if(type1 == 'e'){
                                            my_histos["h_Electron1_Pt"      + name]->Fill(lvec1.Pt() , w);
                                            my_histos["h_Electron1_Phi"     + name]->Fill(lvec1.Phi(), w);
                                            my_histos["h_Electron1_Eta"     + name]->Fill(lvec1.Eta(), w);
                                            my_histos["h_Electron1_Charge"  + name]->Fill(charge1,     w);
                                            my_histos["h_Electron1_MiniIso" + name]->Fill(iso1,        w);
                                        } else if (type1 == 'm' or type1 == 'n') {
                                            my_histos["h_Muon1_Pt"      + name]->Fill(lvec1.Pt(),  w);
                                            my_histos["h_Muon1_Phi"     + name]->Fill(lvec1.Phi(), w);
                                            my_histos["h_Muon1_Eta"     + name]->Fill(lvec1.Eta(), w);
                                            my_histos["h_Muon1_Charge"  + name]->Fill(charge1,     w);
                                            my_histos["h_Muon1_MiniIso" + name]->Fill(iso1,        w);
                                        }

                                    } else if (theGoodLeptons.size() == 2)
                                    {
                                        unsigned int leadIndex    = 0;
                                        unsigned int subleadIndex = 1;
                                        if (theGoodLeptons[1].second.Pt() > theGoodLeptons[0].second.Pt())
                                        {
                                            leadIndex    = 1;
                                            subleadIndex = 0;
                                        }

                                        std::string type1            = theGoodLeptons[leadIndex].first;
                                        utility::LorentzVector lvec1 = theGoodLeptons[leadIndex].second; 
                                        double charge1               = theGoodLeptonsC[leadIndex];
                                        double iso1                  = theGoodLeptonsI[leadIndex];

                                        std::string type2            = theGoodLeptons[subleadIndex].first;
                                        utility::LorentzVector lvec2 = theGoodLeptons[subleadIndex].second; 
                                        double charge2               = theGoodLeptonsC[subleadIndex];
                                        double iso2                  = theGoodLeptonsI[subleadIndex];

                                        my_histos["h_Lepton1_Pt"      + name]->Fill(lvec1.Pt(),  w);
                                        my_histos["h_Lepton1_Phi"     + name]->Fill(lvec1.Phi(), w);
                                        my_histos["h_Lepton1_Eta"     + name]->Fill(lvec1.Eta(), w);
                                        my_histos["h_Lepton1_Charge"  + name]->Fill(charge1,     w);
                                        my_histos["h_Lepton1_MiniIso" + name]->Fill(iso1,        w);

                                        if(type2 == 'e'){
                                            my_histos["h_Electron1_Pt"      + name]->Fill(lvec1.Pt() , w);
                                            my_histos["h_Electron1_Phi"     + name]->Fill(lvec1.Phi(), w);
                                            my_histos["h_Electron1_Eta"     + name]->Fill(lvec1.Eta(), w);
                                            my_histos["h_Electron1_Charge"  + name]->Fill(charge1,     w);
                                            my_histos["h_Electron1_MiniIso" + name]->Fill(iso1,        w);
                                        } else if (type2 == 'm' or type2 == 'n') {
                                            my_histos["h_Muon1_Pt"      + name]->Fill(lvec1.Pt(),  w);
                                            my_histos["h_Muon1_Phi"     + name]->Fill(lvec1.Phi(), w);
                                            my_histos["h_Muon1_Eta"     + name]->Fill(lvec1.Eta(), w);
                                            my_histos["h_Muon1_Charge"  + name]->Fill(charge1,     w);
                                            my_histos["h_Muon1_MiniIso" + name]->Fill(iso1,        w);
                                        }

                                        my_histos["h_Lepton2_Pt"      + name]->Fill(lvec2.Pt(),  w);
                                        my_histos["h_Lepton2_Phi"     + name]->Fill(lvec2.Phi(), w);
                                        my_histos["h_Lepton2_Eta"     + name]->Fill(lvec2.Eta(), w);
                                        my_histos["h_Lepton2_Charge"  + name]->Fill(charge2,     w);
                                        my_histos["h_Lepton2_MiniIso" + name]->Fill(iso2,        w);

                                        if(type2 == 'e'){
                                            my_histos["h_Electron2_Pt"      + name]->Fill(lvec2.Pt() , w);
                                            my_histos["h_Electron2_Phi"     + name]->Fill(lvec2.Phi(), w);
                                            my_histos["h_Electron2_Eta"     + name]->Fill(lvec2.Eta(), w);
                                            my_histos["h_Electron2_Charge"  + name]->Fill(charge2,     w);
                                            my_histos["h_Electron2_MiniIso" + name]->Fill(iso2,        w);
                                        } else if (type2 == 'm' or type2 == 'n') {
                                            my_histos["h_Muon2_Pt"      + name]->Fill(lvec2.Pt(),  w);
                                            my_histos["h_Muon2_Phi"     + name]->Fill(lvec2.Phi(), w);
                                            my_histos["h_Muon2_Eta"     + name]->Fill(lvec2.Eta(), w);
                                            my_histos["h_Muon2_Charge"  + name]->Fill(charge2,     w);
                                            my_histos["h_Muon2_MiniIso" + name]->Fill(iso2,        w);
                                        }

                                    }

                                    // -------------------------------------------------------------------------------------------------
                                    // Fill per jet (non-cm) variables including 4-vector and flavor informatiom
                                    // -------------------------------------------------------------------------------------------------
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
                                        my_histos["h_Jet" + std::to_string(iGoodJet) + "_Pt"     + name]->Fill(theJets.at(j).Pt(),  w); 
                                        my_histos["h_Jet" + std::to_string(iGoodJet) + "_Eta"    + name]->Fill(theJets.at(j).Eta(), w); 
                                        my_histos["h_Jet" + std::to_string(iGoodJet) + "_Phi"    + name]->Fill(theJets.at(j).Phi(), w); 
                                        my_histos["h_Jet" + std::to_string(iGoodJet) + "_Mass"   + name]->Fill(theJets.at(j).M(),   w); 
                                        my_histos["h_Jet" + std::to_string(iGoodJet) + "_Energy" + name]->Fill(theJets.at(j).E(),   w); 

                                        my_histos["h_Jet" + std::to_string(iGoodJet) + "_Flavb"   + name]->Fill(theJetsFlavb.at(j),   w); 
                                        my_histos["h_Jet" + std::to_string(iGoodJet) + "_Flavc"   + name]->Fill(theJetsFlavc.at(j),   w); 
                                        my_histos["h_Jet" + std::to_string(iGoodJet) + "_Flavg"   + name]->Fill(theJetsFlavg.at(j),   w); 
                                        my_histos["h_Jet" + std::to_string(iGoodJet) + "_Flavuds" + name]->Fill(theJetsFlavuds.at(j), w); 
                                        my_histos["h_Jet" + std::to_string(iGoodJet) + "_Flavq"   + name]->Fill(theJetsFlavq.at(j),   w); 

                                        iGoodJet++;

                                        if (iGoodJet > nMVAjets)
                                            break;
                                    }                        

                                    Debug("Filling combined jet variable histograms with name: " + name, __LINE__);

                                    my_histos["h_combined" + std::to_string(nMVAjets) + "thJet_Pt"     + name]->Fill(theCombJetPt,    w); 
                                    my_histos["h_combined" + std::to_string(nMVAjets) + "thJet_PtrHT"  + name]->Fill(theCombJetPtrHT, w); 
                                    my_histos["h_combined" + std::to_string(nMVAjets) + "thJet_Eta"    + name]->Fill(theCombJetEta,   w); 
                                    my_histos["h_combined" + std::to_string(nMVAjets) + "thJet_Phi"    + name]->Fill(theCombJetPhi,   w); 
                                    my_histos["h_combined" + std::to_string(nMVAjets) + "thJet_Mass"   + name]->Fill(theCombJetM,     w); 
                                    my_histos["h_combined" + std::to_string(nMVAjets) + "thJet_Energy" + name]->Fill(theCombJetE,     w); 
                                }                            

                                // -------------------------------------------------------------------------------------------------
                                // Fill 2D event-level varaibles
                                // -------------------------------------------------------------------------------------------------
                                if (region == "ABCD" and ttvar == "nom" and jecvar == "")
                                    my_2d_histos["h_nRtops_vs_nMtops"      + name]->Fill(!isQCD ? NRTops[channel] : NRTops_CR[channel], !isQCD ? NMTops[channel] : NMTops_CR[channel], w);

                                auto& cmJets = Jets_cm_top6[channel];
                                if (isQCD)
                                    cmJets = Jets_cm_top6_CR[channel];

                                Debug("Filling 2D pt HT variable histograms with name: " + name, __LINE__);

                                if (region == "ABCD" and ttvar == "nom" and jecvar == "")
                                {
                                    for(unsigned int i = 1; i <= nJets; i++)
                                    {
                                        double pt  = static_cast<double>(cmJets.at(i-1).Pt());
                                        my_2d_histos["h_pt_jetRank_cm"    + name]->Fill(pt,    i, w);
                                        my_2d_histos["h_ptrHT_jetRank_cm" + name]->Fill(pt/ht, i, w);
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

                                    int shift = inRegionBin[iSubRegion] ? iSubRegion : 100;

                                    auto& theGoodJets = NGoodJets[channel];
                                    if (isQCD)
                                    {
                                        theGoodJets = NGoodJets_CR[channel];
                                    }

                                    // ---------------------------------------------------------
                                    // make inputs for Higgs Combine
                                    //  -- it is 1D histogram, including 4 groups A, B, C, D
                                    //  -- each of A, B, C, D has 7,8,9,10,11, 12incl Njets bins
                                    // ---------------------------------------------------------
                                    if (njets == "Incl")
                                    {
                                        my_histos["h_njets_11incl" + name]->Fill(theGoodJets>=11 ? 11 - 6 + shift*6 : theGoodJets - 6 + shift*6, w);
                                        my_histos["h_njets_12incl" + name]->Fill(theGoodJets>=12 ? 12 - 7 + shift*6 : theGoodJets - 7 + shift*6, w);
                                        my_histos["h_njets_13incl" + name]->Fill(theGoodJets>=13 ? 13 - 8 + shift*6 : theGoodJets - 8 + shift*6, w);
                                    }
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
