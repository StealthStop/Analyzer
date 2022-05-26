#define AnalyzeDataMC_cxx
#include "Analyzer/Analyzer/include/AnalyzeDataMC.h"
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

#include <TH1D.h>
#include <TH2D.h>
#include <iostream>

AnalyzeDataMC::AnalyzeDataMC() : initHistos(false)
{

    histInfos = {
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
        {"h_lMiniIso",                  200,    0,  0.2},
        {"h_ePt",                       200,    0, 2000},
        {"h_ePhi",                      200,   -4,    4},
        {"h_eEta",                      200,   -6,    6},
        {"h_eCharge",                     2,   -1,    1},
        {"h_eMiniIso",                  200,    0,  0.2},
        {"h_mPt",                       360,    0, 1500},
        {"h_mPhi",                      200,   -4,    4},
        {"h_mEta",                      200,   -6,    6},
        {"h_mCharge",                     2,   -1,    1},
        {"h_mMiniIso",                  200,    0,  0.2},
    };

    hist2DInfos = {
        {"h_cm_ptrHT_jetRank", 100, 0, 1, 10, 0, 10},
        {"h_cm_pt_jetRank", 150, 0, 1500, 10, 0, 10},
        {"h_cm_ptrHT_fwm2_top6", 100, 0, 1, 50, 0, 1},
        {"h_cm_ptrHT_fwm3_top6", 100, 0, 1, 50, 0, 1},
        {"h_cm_ptrHT_fwm4_top6", 100, 0, 1, 50, 0, 1},
        {"h_cm_ptrHT_fwm5_top6", 100, 0, 1, 50, 0, 1},
        {"h_cm_ptrHT_jmt_ev0_top6", 100, 0, 1, 50, 0, 1},
        {"h_cm_ptrHT_jmt_ev1_top6", 100, 0, 1, 50, 0, 1},
        {"h_cm_ptrHT_jmt_ev2_top6", 100, 0, 1, 50, 0, 1},
    };

    njets = {"Incl", "7", "8", "9", "10", "11Incl"};

    my_var_suffix = {""};
    //my_var_suffix = {"JECup"};
    //my_var_suffix = {"JECdown"};
    //my_var_suffix = {"JERup"};
    //my_var_suffix = {"JERdown"};
    //my_var_suffix = {"", "JECup", "JECdown", "JERup", "JERdown"};


}

void AnalyzeDataMC::Preinit(unsigned int nNNJets)
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

void AnalyzeDataMC::InitHistos(const std::map<std::string, bool>& cutMap)
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

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

                std::string name = hInfo.name + mycut.first + njetStr;
                my_histos.emplace(name, std::make_shared<TH1D>((name).c_str(),(name).c_str(), hInfo.nBins, hInfo.low, hInfo.high));
            }
        }

        for(const auto& h2dInfo : hist2DInfos)
        {
            for (const auto& Njet : njets)
            {
                if (Njet == "Incl")
                    continue;

                std::string njetStr = "";
                // For 2D histos, don't make any where we exclude all but one njets bin
                if (Njet != "Incl")
                    njetStr = "_Njets" + Njet;

                std::string name = h2dInfo.name + mycut.first + njetStr;
                my_2d_histos.emplace(name, std::make_shared<TH2D>((name).c_str(),(name).c_str(), h2dInfo.nBinsX, h2dInfo.lowX, h2dInfo.highX, h2dInfo.nBinsY, h2dInfo.lowY, h2dInfo.highY));
            }
            std::string name = h2dInfo.name + mycut.first;
            my_2d_histos.emplace(name, std::make_shared<TH2D>((name).c_str(),(name).c_str(), h2dInfo.nBinsX, h2dInfo.lowX, h2dInfo.highX, h2dInfo.nBinsY, h2dInfo.lowY, h2dInfo.highY));
        }
    }
}

void AnalyzeDataMC::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    const auto& filetag         = tr.getVar<std::string>("filetag");
    const auto& runYear         = tr.getVar<std::string>("runYear");
    const auto& bjetFileName    = tr.getVar<std::string>("bjetFileName");
    const auto& bjetCSVFileName = tr.getVar<std::string>("bjetCSVFileName");
    const auto& leptonFileName  = tr.getVar<std::string>("leptonFileName");
    const auto& meanFileName    = tr.getVar<std::string>("meanFileName");
    const auto& TopTaggerCfg    = tr.getVar<std::string>("TopTaggerCfg");

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
        ScaleFactors        scaleFactors( runYear, leptonFileName, meanFileName, myVarSuffix);
        StopGenMatch        stopGenMatch(myVarSuffix);
        FatJetCombine       fatJetCombine(myVarSuffix);
        BTagCorrector       bTagCorrector(bjetFileName, "", bjetCSVFileName, filetag);
        CommonVariables     commonVariables(myVarSuffix);
        MakeMVAVariables    makeMVAVariables0L_NonIsoMuon(false, myVarSuffix, "GoodJets_pt30",       false, true, 7,  2, "_0l");
        MakeMVAVariables    makeMVAVariables1L_NonIsoMuon(false, myVarSuffix, "NonIsoMuonJets_pt30", false, true, 7,  2, "_1l");
        MakeMVAVariables    makeMVAVariables0L(false, myVarSuffix, "GoodJets_pt30", false, true, 7, 2, "_0l");
        MakeMVAVariables    makeMVAVariables1L(false, myVarSuffix, "GoodJets_pt30", false, true, 7, 2, "_1l");
        MakeStopHemispheres stopHemispheres_OldSeed("Jets",     "GoodJets_pt20", "NGoodJets_pt20", "_OldSeed", myVarSuffix, Hemisphere::InvMassSeed);
        MakeStopHemispheres stopHemispheres_TopSeed("StopJets", "GoodStopJets",  "NGoodStopJets",  "_TopSeed", myVarSuffix, Hemisphere::TopSeed);
        MakeStopHemispheres stopHemispheres_OldSeed_NonIsoMuon("Jets",     "NonIsoMuonJets_pt20",     "NNonIsoMuonJets_pt30",     "_OldSeed_NonIsoMuon", myVarSuffix, Hemisphere::InvMassSeed);
        MakeStopHemispheres stopHemispheres_TopSeed_NonIsoMuon("StopJets",     "GoodStopJets",     "NGoodStopJets",     "_TopSeed_NonIsoMuon", myVarSuffix, Hemisphere::InvMassSeed);
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
            const auto& Jets                              = tr.getVec<utility::LorentzVector>("Jets"+myVarSuffix);
            const auto& GoodJets_pt30                     = tr.getVec<bool>("GoodJets_pt30"+myVarSuffix);
            const auto& NGoodJets_pt30                    = tr.getVar<int>("NGoodJets_pt30"+myVarSuffix);
            const auto& NNonIsoMuonJets_pt30              = tr.getVar<int>("NNonIsoMuonJets_pt30"+myVarSuffix);
            const auto& HT_trigger_pt30                   = tr.getVar<double>("HT_trigger_pt30"+myVarSuffix);
            const auto& HT_NonIsoMuon_pt30                = tr.getVar<double>("HT_NonIsoMuon_pt30"+myVarSuffix);

            const auto& passBaseline0l_Good               = tr.getVar<bool>("passBaseline0l_good"+myVarSuffix);
            const auto& passBaseline1l_Good               = tr.getVar<bool>("passBaseline1l_Good"+myVarSuffix);
            const auto& passBaseline0l_Good_blind         = tr.getVar<bool>("passBaseline0l_good_blind"+myVarSuffix);
            const auto& passBaseline1l_Good_blind         = tr.getVar<bool>("passBaseline1l_Good_blind"+myVarSuffix);
            const auto& passBaseline1l_HT500_Good         = tr.getVar<bool>("passBaseline1l_HT500_Good"+myVarSuffix);
            const auto& passBaseline1l_HT700_Good         = tr.getVar<bool>("passBaseline1l_HT700_Good"+myVarSuffix);
            const auto& passBaseline0l_NonIsoMuon         = tr.getVar<bool>("pass_qcdCR"+myVarSuffix);
            const auto& passBaseline1l_NonIsoMuon         = tr.getVar<bool>("passBaseline1l_NonIsoMuon"+myVarSuffix);

            const auto& fwm2_top6_0l                      = tr.getVar<double>("fwm2_top6_0l"+myVarSuffix);
            const auto& fwm3_top6_0l                      = tr.getVar<double>("fwm3_top6_0l"+myVarSuffix);
            const auto& fwm4_top6_0l                      = tr.getVar<double>("fwm4_top6_0l"+myVarSuffix);
            const auto& fwm5_top6_0l                      = tr.getVar<double>("fwm5_top6_0l"+myVarSuffix);
            const auto& jmt_ev0_top6_0l                   = tr.getVar<double>("jmt_ev0_top6_0l"+myVarSuffix);
            const auto& jmt_ev1_top6_0l                   = tr.getVar<double>("jmt_ev1_top6_0l"+myVarSuffix);
            const auto& jmt_ev2_top6_0l                   = tr.getVar<double>("jmt_ev2_top6_0l"+myVarSuffix);
            const auto& fwm2_top6_1l                      = tr.getVar<double>("fwm2_top6_1l"+myVarSuffix);
            const auto& fwm3_top6_1l                      = tr.getVar<double>("fwm3_top6_1l"+myVarSuffix);
            const auto& fwm4_top6_1l                      = tr.getVar<double>("fwm4_top6_1l"+myVarSuffix);
            const auto& fwm5_top6_1l                      = tr.getVar<double>("fwm5_top6_1l"+myVarSuffix);
            const auto& jmt_ev0_top6_1l                   = tr.getVar<double>("jmt_ev0_top6_1l"+myVarSuffix);
            const auto& jmt_ev1_top6_1l                   = tr.getVar<double>("jmt_ev1_top6_1l"+myVarSuffix);
            const auto& jmt_ev2_top6_1l                   = tr.getVar<double>("jmt_ev2_top6_1l"+myVarSuffix);
            const auto& NonIsoMuons_fwm2_top6_1l          = tr.getVar<double>("NonIsoMuons_fwm2_top6_1l"+myVarSuffix);
            const auto& NonIsoMuons_fwm3_top6_1l          = tr.getVar<double>("NonIsoMuons_fwm3_top6_1l"+myVarSuffix);
            const auto& NonIsoMuons_fwm4_top6_1l          = tr.getVar<double>("NonIsoMuons_fwm4_top6_1l"+myVarSuffix);
            const auto& NonIsoMuons_fwm5_top6_1l          = tr.getVar<double>("NonIsoMuons_fwm5_top6_1l"+myVarSuffix);
            const auto& NonIsoMuons_jmt_ev0_top6_1l       = tr.getVar<double>("NonIsoMuons_jmt_ev0_top6_1l"+myVarSuffix);
            const auto& NonIsoMuons_jmt_ev1_top6_1l       = tr.getVar<double>("NonIsoMuons_jmt_ev1_top6_1l"+myVarSuffix);
            const auto& NonIsoMuons_jmt_ev2_top6_1l       = tr.getVar<double>("NonIsoMuons_jmt_ev2_top6_1l"+myVarSuffix);

            const auto& Jets_cm_top6_0l                   = tr.getVec<utility::LorentzVector>("Jets_cm_top6_0l"+myVarSuffix);
            const auto& Jets_cm_top6_1l                   = tr.getVec<utility::LorentzVector>("Jets_cm_top6_1l"+myVarSuffix);
            const auto& NonIsoMuons_Jets_cm_top6_1l       = tr.getVec<utility::LorentzVector>("NonIsoMuons_Jets_cm_top6_1l"+myVarSuffix);
            const auto& nMVAJets_0l                       = tr.getVar<unsigned int>("nMVAJets_0l"+myVarSuffix);
            const auto& nMVAJets_1l                       = tr.getVar<unsigned int>("nMVAJets_1l"+myVarSuffix);

            const auto& eventCounter                      = tr.getVar<int>("eventCounter");
            const auto& Stop1_pt_cm_OldSeed               = tr.getVar<double>("Stop1_pt_cm_OldSeed"+myVarSuffix);
            const auto& Stop1_eta_cm_OldSeed              = tr.getVar<double>("Stop1_eta_cm_OldSeed"+myVarSuffix);
            const auto& Stop1_phi_cm_OldSeed              = tr.getVar<double>("Stop1_phi_cm_OldSeed"+myVarSuffix);
            const auto& Stop1_mass_cm_OldSeed             = tr.getVar<double>("Stop1_mass_cm_OldSeed"+myVarSuffix);
            const auto& Stop2_pt_cm_OldSeed               = tr.getVar<double>("Stop2_pt_cm_OldSeed"+myVarSuffix);
            const auto& Stop2_eta_cm_OldSeed              = tr.getVar<double>("Stop2_eta_cm_OldSeed"+myVarSuffix);
            const auto& Stop2_phi_cm_OldSeed              = tr.getVar<double>("Stop2_phi_cm_OldSeed"+myVarSuffix);
            const auto& Stop2_mass_cm_OldSeed             = tr.getVar<double>("Stop2_mass_cm_OldSeed"+myVarSuffix);
            const auto& Stop1_pt_cm_TopSeed               = tr.getVar<double>("Stop1_pt_cm_TopSeed"+myVarSuffix);
            const auto& Stop1_eta_cm_TopSeed              = tr.getVar<double>("Stop1_eta_cm_TopSeed"+myVarSuffix);
            const auto& Stop1_phi_cm_TopSeed              = tr.getVar<double>("Stop1_phi_cm_TopSeed"+myVarSuffix);
            const auto& Stop1_mass_cm_TopSeed             = tr.getVar<double>("Stop1_mass_cm_TopSeed"+myVarSuffix);
            const auto& Stop2_pt_cm_TopSeed               = tr.getVar<double>("Stop2_pt_cm_TopSeed"+myVarSuffix);
            const auto& Stop2_eta_cm_TopSeed              = tr.getVar<double>("Stop2_eta_cm_TopSeed"+myVarSuffix);
            const auto& Stop2_phi_cm_TopSeed              = tr.getVar<double>("Stop2_phi_cm_TopSeed"+myVarSuffix);
            const auto& Stop2_mass_cm_TopSeed             = tr.getVar<double>("Stop2_mass_cm_TopSeed"+myVarSuffix);

            const auto& Stop1_pt_cm_OldSeed_NonIsoMuon    = tr.getVar<double>("Stop1_pt_cm_OldSeed_NonIsoMuon"+myVarSuffix);
            const auto& Stop1_eta_cm_OldSeed_NonIsoMuon   = tr.getVar<double>("Stop1_eta_cm_OldSeed_NonIsoMuon"+myVarSuffix);
            const auto& Stop1_phi_cm_OldSeed_NonIsoMuon   = tr.getVar<double>("Stop1_phi_cm_OldSeed_NonIsoMuon"+myVarSuffix);
            const auto& Stop1_mass_cm_OldSeed_NonIsoMuon  = tr.getVar<double>("Stop1_mass_cm_OldSeed_NonIsoMuon"+myVarSuffix);
            const auto& Stop2_pt_cm_OldSeed_NonIsoMuon    = tr.getVar<double>("Stop2_pt_cm_OldSeed_NonIsoMuon"+myVarSuffix);
            const auto& Stop2_eta_cm_OldSeed_NonIsoMuon   = tr.getVar<double>("Stop2_eta_cm_OldSeed_NonIsoMuon"+myVarSuffix);
            const auto& Stop2_phi_cm_OldSeed_NonIsoMuon   = tr.getVar<double>("Stop2_phi_cm_OldSeed_NonIsoMuon"+myVarSuffix);
            const auto& Stop2_mass_cm_OldSeed_NonIsoMuon  = tr.getVar<double>("Stop2_mass_cm_OldSeed_NonIsoMuon"+myVarSuffix);
            const auto& Stop1_pt_cm_TopSeed_NonIsoMuon    = tr.getVar<double>("Stop1_pt_cm_TopSeed_NonIsoMuon"+myVarSuffix);
            const auto& Stop1_eta_cm_TopSeed_NonIsoMuon   = tr.getVar<double>("Stop1_eta_cm_TopSeed_NonIsoMuon"+myVarSuffix);
            const auto& Stop1_phi_cm_TopSeed_NonIsoMuon   = tr.getVar<double>("Stop1_phi_cm_TopSeed_NonIsoMuon"+myVarSuffix);
            const auto& Stop1_mass_cm_TopSeed_NonIsoMuon  = tr.getVar<double>("Stop1_mass_cm_TopSeed_NonIsoMuon"+myVarSuffix);
            const auto& Stop2_pt_cm_TopSeed_NonIsoMuon    = tr.getVar<double>("Stop2_pt_cm_TopSeed_NonIsoMuon"+myVarSuffix);
            const auto& Stop2_eta_cm_TopSeed_NonIsoMuon   = tr.getVar<double>("Stop2_eta_cm_TopSeed_NonIsoMuon"+myVarSuffix);
            const auto& Stop2_phi_cm_TopSeed_NonIsoMuon   = tr.getVar<double>("Stop2_phi_cm_TopSeed_NonIsoMuon"+myVarSuffix);
            const auto& Stop2_mass_cm_TopSeed_NonIsoMuon  = tr.getVar<double>("Stop2_mass_cm_TopSeed_NonIsoMuon"+myVarSuffix);

            const auto& GoodLeptons                       = tr.getVec<std::pair<std::string, utility::LorentzVector>>("GoodLeptons"+myVarSuffix);
            const auto& GoodLeptonsCharge                 = tr.getVec<int>("GoodLeptonsCharge"+myVarSuffix);
            const auto& GoodLeptonsMiniIso                = tr.getVec<double>("GoodLeptonsMiniIso"+myVarSuffix);

            const auto& Stop1_mass_PtRank_matched         = runtype == "Data" ? -999.0 : tr.getVar<float>("stop1_ptrank_mass"+myVarSuffix);
            const auto& Stop2_mass_PtRank_matched         = runtype == "Data" ? -999.0 : tr.getVar<float>("stop2_ptrank_mass"+myVarSuffix);
            const auto& Stop1_mass_MassRank_matched       = runtype == "Data" ? -999.0 : tr.getVar<float>("stop1_mrank_mass"+myVarSuffix);
            const auto& Stop2_mass_MassRank_matched       = runtype == "Data" ? -999.0 : tr.getVar<float>("stop2_mrank_mass"+myVarSuffix);
            const auto& Stop_mass_average_matched         = runtype == "Data" ? -999.0 : tr.getVar<double>("stop_avemass"+myVarSuffix);

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
                Jets_flavb_0l.push_back(tr.getVar<double>("Jet_flavb_"+std::to_string(iJet)+"_0l"+myVarSuffix));
                Jets_flavc_0l.push_back(tr.getVar<double>("Jet_flavc_"+std::to_string(iJet)+"_0l"+myVarSuffix));
                Jets_flavg_0l.push_back(tr.getVar<double>("Jet_flavg_"+std::to_string(iJet)+"_0l"+myVarSuffix));
                Jets_flavuds_0l.push_back(tr.getVar<double>("Jet_flavuds_"+std::to_string(iJet)+"_0l"+myVarSuffix));
                Jets_flavq_0l.push_back(tr.getVar<double>("Jet_flavq_"+std::to_string(iJet)+"_0l"+myVarSuffix));
            }

            for (unsigned int iJet = 1; iJet <= nMVAJets_1l; iJet++) {
                Jets_flavb_1l.push_back(tr.getVar<double>("Jet_flavb_"+std::to_string(iJet)+"_1l"+myVarSuffix));
                Jets_flavc_1l.push_back(tr.getVar<double>("Jet_flavc_"+std::to_string(iJet)+"_1l"+myVarSuffix));
                Jets_flavg_1l.push_back(tr.getVar<double>("Jet_flavg_"+std::to_string(iJet)+"_1l"+myVarSuffix));
                Jets_flavuds_1l.push_back(tr.getVar<double>("Jet_flavuds_"+std::to_string(iJet)+"_1l"+myVarSuffix));
                Jets_flavq_1l.push_back(tr.getVar<double>("Jet_flavq_"+std::to_string(iJet)+"_1l"+myVarSuffix));

                JetNonIsoMuons_flavb_1l.push_back(  tr.getVar<double>("JetNonIsoMuons_flavb_"+std::to_string(iJet)+"_1l"+myVarSuffix));
                JetNonIsoMuons_flavc_1l.push_back(  tr.getVar<double>("JetNonIsoMuons_flavc_"+std::to_string(iJet)+"_1l"+myVarSuffix));
                JetNonIsoMuons_flavg_1l.push_back(  tr.getVar<double>("JetNonIsoMuons_flavg_"+std::to_string(iJet)+"_1l"+myVarSuffix));
                JetNonIsoMuons_flavuds_1l.push_back(tr.getVar<double>("JetNonIsoMuons_flavuds_"+std::to_string(iJet)+"_1l"+myVarSuffix));
                JetNonIsoMuons_flavq_1l.push_back(  tr.getVar<double>("JetNonIsoMuons_flavq_"+std::to_string(iJet)+"_1l"+myVarSuffix));
            }

            // Put 0L and 1L version of variables into vector
            // 0th position is for 0L and 1st position is for 1L for convenience in the event loop
            // For 0L things in the CR, nominal version of jets and derived quantities are used, take note !
            std::vector<int>                                        NGoodJets                     {NGoodJets_pt30,                        NGoodJets_pt30};
            std::vector<int>                                        NNonIsoMuonJets               {NGoodJets_pt30,                        NNonIsoMuonJets_pt30};
            std::vector<double>                                     HT_trigger                    {HT_trigger_pt30,                       HT_trigger_pt30};
            std::vector<double>                                     HT_NonIsoMuon                 {HT_trigger_pt30,                       HT_NonIsoMuon_pt30};
            std::vector<double>                                     fwm2_top6                     {fwm2_top6_0l,                          fwm2_top6_1l};
            std::vector<double>                                     fwm3_top6                     {fwm3_top6_0l,                          fwm3_top6_1l};
            std::vector<double>                                     fwm4_top6                     {fwm4_top6_0l,                          fwm4_top6_1l};
            std::vector<double>                                     fwm5_top6                     {fwm5_top6_0l,                          fwm5_top6_1l};
            std::vector<double>                                     jmt_ev0_top6                  {jmt_ev0_top6_0l,                       jmt_ev0_top6_1l};
            std::vector<double>                                     jmt_ev1_top6                  {jmt_ev1_top6_0l,                       jmt_ev1_top6_1l};
            std::vector<double>                                     jmt_ev2_top6                  {jmt_ev2_top6_0l,                       jmt_ev2_top6_1l};
            std::vector<double>                                     fwm2_top6_QCDCR               {fwm2_top6_0l,                          NonIsoMuons_fwm2_top6_1l};
            std::vector<double>                                     fwm3_top6_QCDCR               {fwm3_top6_0l,                          NonIsoMuons_fwm3_top6_1l};
            std::vector<double>                                     fwm4_top6_QCDCR               {fwm4_top6_0l,                          NonIsoMuons_fwm4_top6_1l};
            std::vector<double>                                     fwm5_top6_QCDCR               {fwm5_top6_0l,                          NonIsoMuons_fwm5_top6_1l};
            std::vector<double>                                     jmt_ev0_top6_QCDCR            {jmt_ev0_top6_0l,                       NonIsoMuons_jmt_ev0_top6_1l};
            std::vector<double>                                     jmt_ev1_top6_QCDCR            {jmt_ev1_top6_0l,                       NonIsoMuons_jmt_ev1_top6_1l};
            std::vector<double>                                     jmt_ev2_top6_QCDCR            {jmt_ev2_top6_0l,                       NonIsoMuons_jmt_ev2_top6_1l};
            std::vector<unsigned int>                               nMVAJets                      {nMVAJets_0l,                           nMVAJets_1l};
            std::vector<std::vector<double> >                       Jets_flavb                    {Jets_flavb_0l,                         Jets_flavb_1l};
            std::vector<std::vector<double> >                       Jets_flavc                    {Jets_flavc_0l,                         Jets_flavc_1l};
            std::vector<std::vector<double> >                       Jets_flavg                    {Jets_flavg_0l,                         Jets_flavg_1l};
            std::vector<std::vector<double> >                       Jets_flavuds                  {Jets_flavuds_0l,                       Jets_flavuds_1l};
            std::vector<std::vector<double> >                       Jets_flavq                    {Jets_flavq_0l,                         Jets_flavq_1l};
            std::vector<std::vector<double> >                       Jets_flavb_QCDCR              {Jets_flavb_0l,                         JetNonIsoMuons_flavb_1l};
            std::vector<std::vector<double> >                       Jets_flavc_QCDCR              {Jets_flavc_0l,                         JetNonIsoMuons_flavc_1l};
            std::vector<std::vector<double> >                       Jets_flavg_QCDCR              {Jets_flavg_0l,                         JetNonIsoMuons_flavg_1l};
            std::vector<std::vector<double> >                       Jets_flavuds_QCDCR            {Jets_flavuds_0l,                       JetNonIsoMuons_flavuds_1l};
            std::vector<std::vector<double> >                       Jets_flavq_QCDCR              {Jets_flavq_0l,                         JetNonIsoMuons_flavq_1l};
            std::vector<std::vector<utility::LorentzVector> >               Jets_cm_top6                  {Jets_cm_top6_0l,                       Jets_cm_top6_1l};
            std::vector<std::vector<utility::LorentzVector> >               Jets_cm_top6_QCDCR            {Jets_cm_top6_0l,                       NonIsoMuons_Jets_cm_top6_1l};

            // ------------------------
            // -- Define weight
            // ------------------------
            double eventweight = 1.0, leptonweight = 1.0, bTagWeight = 1.0, prefiringScaleFactor = 1.0, pileupWeight = 1.0, htDerivedweight = 1.0;
            double weight0L             = 1.0, weight1L            = 1.0;
            double weight0LNoB          = 1.0, weight1LNoB         = 1.0;
            double weight0L_NonIsoMuon  = 1.0, weight1L_NonIsoMuon = 1.0;
            double weight0L_topReweight = 1.0, weight1L_topReweight = 1.0;
            double weight0L_NonIsoMuon_topReweight = 1.0, weight1L_NonIsoMuon_topReweight = 1.0;
            if(runtype == "MC" )
            {
                // Define Lumi weight
                const auto& Weight = tr.getVar<float>("Weight");
                const auto& lumi   = tr.getVar<double>("Lumi");
                eventweight        = lumi * Weight;
                
                const auto& eleLepWeight     = tr.getVar<double>("totGoodElectronSF"+myVarSuffix);
                const auto& muLepWeight      = tr.getVar<double>("totGoodMuonSF"+myVarSuffix);
                const auto& muNonIso         = tr.getVar<double>("totNonIsoMuonSF"+myVarSuffix);
                const auto& topPtScaleFactor = tr.getVar<double>("topPtScaleFactor"+myVarSuffix);
                leptonweight             = eleLepWeight * muLepWeight;

                pileupWeight         = tr.getVar<double>("puWeightCorr"+myVarSuffix);
                bTagWeight           = tr.getVar<double>("bTagSF_EventWeightSimple_Central"+myVarSuffix);
                htDerivedweight      = tr.getVar<double>("htDerivedweight"+myVarSuffix);
                prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor"+myVarSuffix);

                weight1L             *= eventweight * leptonweight * bTagWeight * prefiringScaleFactor * pileupWeight; // * htDerivedweight;
                weight0L             *= eventweight *                bTagWeight * prefiringScaleFactor * pileupWeight; // * htDerivedweight;

                weight1LNoB          *= eventweight * leptonweight * prefiringScaleFactor * pileupWeight; // * htDerivedweight;
                weight0LNoB          *= eventweight *                prefiringScaleFactor * pileupWeight; // * htDerivedweight;

                weight1L_topReweight *= eventweight * leptonweight * bTagWeight * prefiringScaleFactor * pileupWeight * topPtScaleFactor; // * htDerivedweight;
                weight0L_topReweight *= eventweight *                bTagWeight * prefiringScaleFactor * pileupWeight * topPtScaleFactor; // * htDerivedweight;

                weight1L_NonIsoMuon  *= eventweight * muNonIso                  * prefiringScaleFactor * pileupWeight;
                weight0L_NonIsoMuon  *= eventweight * muNonIso                  * prefiringScaleFactor * pileupWeight;

                weight1L_NonIsoMuon_topReweight  *= eventweight * muNonIso                  * prefiringScaleFactor * pileupWeight * topPtScaleFactor;
                weight0L_NonIsoMuon_topReweight  *= eventweight * muNonIso                  * prefiringScaleFactor * pileupWeight * topPtScaleFactor;
            }

            std::vector<double> weight                      {weight0L,                         weight1L};
            std::vector<double> weight_NoB                  {weight0LNoB,                      weight1LNoB};
            std::vector<double> weight_QCDCR                {weight0L_NonIsoMuon,              weight1L_NonIsoMuon};
            std::vector<double> weight_topReweight          {weight0L_topReweight,             weight1L_topReweight};
            std::vector<double> weight_QCDCR_topReweight    {weight0L_NonIsoMuon_topReweight,  weight1L_NonIsoMuon_topReweight};

            const std::map<std::string, bool> cut_map
            {
                {"_1l"               , passBaseline1l_Good},                         
                {"_1l_blind"         , passBaseline1l_Good_blind},                         
                //{"_1l_NoB"           , passBaseline1l_Good},                         
                {"_1l_topPtReweight" , passBaseline1l_Good},                         
                //{"_1l_HT500"         , passBaseline1l_HT500_Good},                         
                //{"_1l_HT700"         , passBaseline1l_HT700_Good},                         
                {"_0l"               , passBaseline0l_Good},                         
                {"_0l_blind"         , passBaseline0l_Good_blind},                         
                //{"_0l_NoB"           , passBaseline0l_Good},                         
                {"_0l_topPtReweight" , passBaseline0l_Good},                         
                {"_1l_QCDCR"         , passBaseline1l_NonIsoMuon},                         
                {"_1l_QCDCR_topPtReweight"         , passBaseline1l_NonIsoMuon},                         
                //{"_0l_QCDCR"         , passBaseline0l_NonIsoMuon},                         
            };

            std::map<std::string, bool>               njetsMap;

            if(!initHistos)
            {
                Preinit(nMVAJets_0l > nMVAJets_1l ? nMVAJets_0l : nMVAJets_1l);
                InitHistos(cut_map);
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
                bool isNoB          = kv.first.find("NoB")      != std::string::npos;
                unsigned int nJets = !isQCD ? Jets_cm_top6[channel].size() : Jets_cm_top6_QCDCR[channel].size();

                // For 0L, we always use the NGoodJets case
                njetsMap = {
                    {"Incl",     true},
                    {"7",      (!isQCD or channel==0) ? NGoodJets[channel]==7  : NNonIsoMuonJets[channel]==7},
                    {"8",      (!isQCD or channel==0) ? NGoodJets[channel]==8  : NNonIsoMuonJets[channel]==8},
                    {"9",      (!isQCD or channel==0) ? NGoodJets[channel]==9  : NNonIsoMuonJets[channel]==9},
                    {"10",     (!isQCD or channel==0) ? NGoodJets[channel]==10  : NNonIsoMuonJets[channel]==10},
                    {"11Incl", (!isQCD or channel==0) ? NGoodJets[channel]>=11  : NNonIsoMuonJets[channel]>=11},
                };

                double w  = !isQCD ? weight[channel]        : weight_QCDCR[channel];
                double ht = !isQCD ? HT_trigger[channel]    : HT_NonIsoMuon[channel];

                w  = !isTopReweight ? w                     : weight_topReweight[channel];
                w  = !(isTopReweight && isQCD) ? w          : weight_QCDCR_topReweight[channel];
                w  = !isNoB ? w                             : weight_NoB[channel];

                std::string name;
                for (auto& njetsPass : njetsMap)
                {
                    std::string njets      = njetsPass.first;
                    bool        inNjetsBin = njetsPass.second;
                    std::string njetsStr = "";
                    if (njets != "Incl")  
                        njetsStr = "_Njets" + njets;

                    if (kv.second and inNjetsBin)
                    {
                        name = kv.first + njetsStr;
                        if (njets == "Incl")
                        {

                            my_histos["h_njets"      + name]->Fill(!isQCD ? NGoodJets[channel]    : NNonIsoMuonJets[channel], w);
                            my_histos["fwm2_top6"    + name]->Fill(!isQCD ? fwm2_top6[channel]    : fwm2_top6_QCDCR[channel], w);
                            my_histos["fwm3_top6"    + name]->Fill(!isQCD ? fwm3_top6[channel]    : fwm3_top6_QCDCR[channel], w);
                            my_histos["fwm4_top6"    + name]->Fill(!isQCD ? fwm4_top6[channel]    : fwm4_top6_QCDCR[channel], w);
                            my_histos["fwm5_top6"    + name]->Fill(!isQCD ? fwm5_top6[channel]    : fwm5_top6_QCDCR[channel], w);
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
                                    pt  = static_cast<double>(Jets_cm_top6[channel].at(i-1).Pt());
                                    eta = static_cast<double>(Jets_cm_top6[channel].at(i-1).Eta());
                                    phi = static_cast<double>(Jets_cm_top6[channel].at(i-1).Phi());
                                    m   = static_cast<double>(Jets_cm_top6[channel].at(i-1).M());
                                    E   = static_cast<double>(Jets_cm_top6[channel].at(i-1).E());
                                } else
                                {
                                    pt  = static_cast<double>(Jets_cm_top6_QCDCR[channel].at(i-1).Pt());
                                    eta = static_cast<double>(Jets_cm_top6_QCDCR[channel].at(i-1).Eta());
                                    phi = static_cast<double>(Jets_cm_top6_QCDCR[channel].at(i-1).Phi());
                                    m   = static_cast<double>(Jets_cm_top6_QCDCR[channel].at(i-1).M());
                                    E   = static_cast<double>(Jets_cm_top6_QCDCR[channel].at(i-1).E());
                                } 

                                my_histos["Jet_cm_ptrHT_"  + std::to_string(i) + name]->Fill(pt/ht, w);
                                my_histos["Jet_cm_pt_"     + std::to_string(i) + name]->Fill(pt, w);
                                my_histos["Jet_cm_eta_"    + std::to_string(i) + name]->Fill(eta, w);
                                my_histos["Jet_cm_phi_"    + std::to_string(i) + name]->Fill(phi, w);
                                my_histos["Jet_cm_m_"      + std::to_string(i) + name]->Fill(m, w);
                                my_histos["Jet_cm_E_"      + std::to_string(i) + name]->Fill(E, w);
                                my_2d_histos["h_cm_ptrHT_fwm2_top6" + name]->Fill(pt/ht, !isQCD ? fwm2_top6[channel]    : fwm2_top6_QCDCR[channel], w);
                                my_2d_histos["h_cm_ptrHT_fwm3_top6" + name]->Fill(pt/ht, !isQCD ? fwm3_top6[channel]    : fwm3_top6_QCDCR[channel], w);
                                my_2d_histos["h_cm_ptrHT_fwm4_top6" + name]->Fill(pt/ht, !isQCD ? fwm4_top6[channel]    : fwm4_top6_QCDCR[channel], w);
                                my_2d_histos["h_cm_ptrHT_fwm5_top6" + name]->Fill(pt/ht, !isQCD ? fwm5_top6[channel]    : fwm5_top6_QCDCR[channel], w);
                                my_2d_histos["h_cm_ptrHT_jmt_ev0_top6" + name]->Fill(pt/ht, !isQCD ? jmt_ev0_top6[channel] : jmt_ev0_top6_QCDCR[channel], w);
                                my_2d_histos["h_cm_ptrHT_jmt_ev1_top6" + name]->Fill(pt/ht, !isQCD ? jmt_ev1_top6[channel] : jmt_ev1_top6_QCDCR[channel], w);
                                my_2d_histos["h_cm_ptrHT_jmt_ev2_top6" + name]->Fill(pt/ht, !isQCD ? jmt_ev2_top6[channel] : jmt_ev2_top6_QCDCR[channel], w);
                                

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
                            for(std::size_t j = 0 ; j < std::size(GoodLeptons); ++j)
                            {
                                auto& type = GoodLeptons[j].first;
                                auto& lvec = GoodLeptons[j].second; 
                                auto& charge = GoodLeptonsCharge[j];
                                auto&  iso = GoodLeptonsMiniIso[j];

                                my_histos["h_lPt"  + name]->Fill(  lvec.Pt()  ,w);
                                my_histos["h_lPhi" + name]->Fill(  lvec.Phi() ,w);
                                my_histos["h_lEta" + name]->Fill(  lvec.Eta() ,w);
                                my_histos["h_lCharge" + name]->Fill(  charge ,w);
                                my_histos["h_lMiniIso" + name]->Fill(  iso ,w);

                                if(type == 'e')
                                {
                                    my_histos["h_ePt"  + name]->Fill(  lvec.Pt()  ,w);
                                    my_histos["h_ePhi" + name]->Fill(  lvec.Phi() ,w);
                                    my_histos["h_eEta" + name]->Fill(  lvec.Eta() ,w);
                                    my_histos["h_eCharge" + name]->Fill(  charge ,w);
                                    my_histos["h_eMiniIso" + name]->Fill(  iso ,w);
                                } else if (type == 'm') 
                                {
                                    my_histos["h_mPt"  + name]->Fill(  lvec.Pt()  ,w);
                                    my_histos["h_mPhi" + name]->Fill(  lvec.Phi() ,w);
                                    my_histos["h_mEta" + name]->Fill(  lvec.Eta() ,w);
                                    my_histos["h_mCharge" + name]->Fill(  charge ,w);
                                    my_histos["h_mMiniIso" + name]->Fill(  iso ,w);
                                }
                            }
                            
                            for(unsigned int j = 0; j < Jets.size(); j++) 
                            {
                                if(!GoodJets_pt30[j]) continue;
                                my_histos["h_jPt" + name]->Fill(Jets.at(j).Pt(), w); 
                                my_histos["h_jEta" + name]->Fill(Jets.at(j).Eta(), w); 
                                my_histos["h_jPhi" + name]->Fill(Jets.at(j).Phi(), w); 
                                my_histos["h_jM" + name]->Fill(Jets.at(j).M(), w); 
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
                            my_2d_histos["h_cm_pt_jetRank"    + name]->Fill(pt, i, w);
                            my_2d_histos["h_cm_ptrHT_jetRank" + name]->Fill(pt/ht, i, w);
                            my_2d_histos["h_cm_ptrHT_fwm2_top6" + name]->Fill(pt/ht, !isQCD ? fwm2_top6[channel]    : fwm2_top6_QCDCR[channel], w);
                            my_2d_histos["h_cm_ptrHT_fwm3_top6" + name]->Fill(pt/ht, !isQCD ? fwm3_top6[channel]    : fwm3_top6_QCDCR[channel], w);
                            my_2d_histos["h_cm_ptrHT_fwm4_top6" + name]->Fill(pt/ht, !isQCD ? fwm4_top6[channel]    : fwm4_top6_QCDCR[channel], w);
                            my_2d_histos["h_cm_ptrHT_fwm5_top6" + name]->Fill(pt/ht, !isQCD ? fwm5_top6[channel]    : fwm5_top6_QCDCR[channel], w);
                            my_2d_histos["h_cm_ptrHT_jmt_ev0_top6" + name]->Fill(pt/ht, !isQCD ? jmt_ev0_top6[channel] : jmt_ev0_top6_QCDCR[channel], w);
                            my_2d_histos["h_cm_ptrHT_jmt_ev1_top6" + name]->Fill(pt/ht, !isQCD ? jmt_ev1_top6[channel] : jmt_ev1_top6_QCDCR[channel], w);
                            my_2d_histos["h_cm_ptrHT_jmt_ev2_top6" + name]->Fill(pt/ht, !isQCD ? jmt_ev2_top6[channel] : jmt_ev2_top6_QCDCR[channel], w);
                        }
                    }
                }
            }
        }
    }
}

void AnalyzeDataMC::WriteHistos(TFile* outfile)
{
    outfile->cd();

    for(const auto& p : my_histos) 
        p.second->Write();

    for(const auto& p : my_2d_histos) 
        p.second->Write();
}
