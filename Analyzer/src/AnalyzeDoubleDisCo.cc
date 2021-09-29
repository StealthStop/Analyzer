#define AnalyzeDoubleDisCo_cxx
#include "Analyzer/Analyzer/include/AnalyzeDoubleDisCo.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "Framework/Framework/include/Utility.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TLorentzVector.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

AnalyzeDoubleDisCo::AnalyzeDoubleDisCo() : initHistos(false)
{

    histInfos = {
        {"h_DoubleDisCo_disc1",    80,    0,    1},
        {"h_DoubleDisCo_disc2",    80,    0,    1},
        {"h_DoubleDisCo_massReg", 150,    0, 1500},
        {"fwm2_top6",              50,    0,    1},
        {"fwm3_top6",              50,    0,    1},
        {"fwm4_top6",              50,    0,    1},
        {"fwm5_top6",              50,    0,    1},
        {"jmt_ev0_top6",           50,    0,    1},
        {"jmt_ev1_top6",           50,    0,    1},
        {"jmt_ev2_top6",           50,    0,    1},
        {"Stop1_pt_cm_OldSeed",    72,    0, 1500},
        {"Stop1_eta_cm_OldSeed",   80,   -6,    6},
        {"Stop1_phi_cm_OldSeed",   64,   -4,    4},
        {"Stop1_mass_cm_OldSeed",  72,    0, 1500},
        {"Stop2_pt_cm_OldSeed",    72,    0, 1500},
        {"Stop2_eta_cm_OldSeed",   80,   -6,    6},
        {"Stop2_phi_cm_OldSeed",   64,   -4,    4},
        {"Stop2_mass_cm_OldSeed",  72,    0, 1500},
        {"h_MET_phi",              80,   -4,    4},
        {"h_MET_pt" ,             150,    0,  150},
        {"h_njets",                21, -0.5, 20.5},
        {"h_njets_11incl",         20, -0.5, 19.5},
        {"h_njets_12incl",         24, -0.5, 23.5},
        {"h_ntops",                 6, -0.5,  5.5},
        {"h_ht",                  500,    0, 5000},
    };

    hist2DInfos = {
        {"h_DoubleDisCo_disc1_disc2",  100,    0,    1, 100,     0,     1}, 
    };

    abcds = {"", "A", "B", "C", "D"};
    njets = {"Incl", "7", "8", "9", "10", "11", "11incl", "12", "12incl"};

}

void AnalyzeDoubleDisCo::Preinit(unsigned int nNNJets)
{
    for(unsigned int i = 1; i <= nNNJets ; i++)
    {
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

void AnalyzeDoubleDisCo::InitHistos(const std::map<std::string, bool>& cutMap)
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    my_histos.emplace("EventCounter", std::make_shared<TH1D>("EventCounter", "EventCounter", 2, -1.1, 1.1) );

    for(auto& mycut : cutMap)
    {
        for(const auto& hInfo : histInfos)
        { 
            for(const auto& Njet : njets)
            {
                std::string njetStr = "";
                if (Njet != "Incl") njetStr = "_Njets" + Njet;

                // For 1D njets histos, don't make any where we exclude all but one njets bin
                if (hInfo.name.find("njets") != std::string::npos and Njet != "Incl")
                    continue;

                for(const auto& region : abcds)
                {

                    // For discriminant plots, don't make for trivially subselection of ABCD region
                    if (hInfo.name.find("DisCo_disc") != std::string::npos and region != "")
                        continue;

                    std::string regionStr = "";
                    if (region != "") regionStr = "_" + region;

                    std::string name = hInfo.name+mycut.first+njetStr+regionStr;
                    my_histos.emplace(name, std::make_shared<TH1D>((name).c_str(),(name).c_str(), hInfo.nBins, hInfo.low, hInfo.high));
                }
            }
        }

        for(const auto& h2dInfo : hist2DInfos)
        {
            for(const auto& Njet : njets)
            {
                std::string njetStr = "";
                if (Njet != "Incl") njetStr = "_Njets" + Njet;

                for(const auto& region : abcds)
                {

                    // For discriminant plots, don't make for trivially subselection of ABCD region
                    if (h2dInfo.name.find("DisCo_disc") != std::string::npos and region != "")
                        continue;

                    std::string regionStr = "";
                    if (region != "") regionStr = "_" + region;

                    std::string name = h2dInfo.name+mycut.first+njetStr+regionStr;
                    my_2d_histos.emplace(name, std::make_shared<TH2D>((name).c_str(),(name).c_str(), h2dInfo.nBinsX, h2dInfo.lowX, h2dInfo.highX, h2dInfo.nBinsY, h2dInfo.lowY, h2dInfo.highY));
                } 
            }
        }
    }
}

void AnalyzeDoubleDisCo::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    while( tr.getNextEvent() )
    {
        const auto& ntops                             = tr.getVar<int>("ntops");
        const auto& runtype                           = tr.getVar<std::string>("runtype");     
        const auto& NGoodJets_pt30                    = tr.getVar<int>("NGoodJets_pt30");
        const auto& NGoodJets_pt45                    = tr.getVar<int>("NGoodJets_pt45");
        const auto& NNonIsoMuonJets_pt30              = tr.getVar<int>("NNonIsoMuonJets_pt30");
        const auto& NNonIsoMuonJets_pt45              = tr.getVar<int>("NNonIsoMuonJets_pt45");
        const auto& HT_trigger_pt30                   = tr.getVar<double>("HT_trigger_pt30");
        const auto& HT_trigger_pt45                   = tr.getVar<double>("HT_trigger_pt45");
        const auto& HT_NonIsoMuon_pt30                = tr.getVar<double>("HT_NonIsoMuon_pt30");
        const auto& HT_NonIsoMuon_pt45                = tr.getVar<double>("HT_NonIsoMuon_pt45");

        const auto& passBaseline0l_Good               = tr.getVar<bool>("passBaseline0l_Good");
        const auto& passBaseline1l_Good               = tr.getVar<bool>("passBaseline1l_Good");
        const auto& passBaseline0l_NonIsoMuon         = tr.getVar<bool>("passBaseline0l_NonIsoMuon");
        const auto& passBaseline1l_NonIsoMuon         = tr.getVar<bool>("passBaseline1l_NonIsoMuon");

        const auto& DoubleDisCo_binA_0l               = tr.getVar<bool>("DoubleDisCo_binA_0l");
        const auto& DoubleDisCo_binB_0l               = tr.getVar<bool>("DoubleDisCo_binB_0l");
        const auto& DoubleDisCo_binC_0l               = tr.getVar<bool>("DoubleDisCo_binC_0l");
        const auto& DoubleDisCo_binD_0l               = tr.getVar<bool>("DoubleDisCo_binD_0l");
        const auto& DoubleDisCo_binA_1l               = tr.getVar<bool>("DoubleDisCo_binA_1l");
        const auto& DoubleDisCo_binB_1l               = tr.getVar<bool>("DoubleDisCo_binB_1l");
        const auto& DoubleDisCo_binC_1l               = tr.getVar<bool>("DoubleDisCo_binC_1l");
        const auto& DoubleDisCo_binD_1l               = tr.getVar<bool>("DoubleDisCo_binD_1l");
        const auto& DoubleDisCo_binA_NonIsoMuon_0l    = tr.getVar<bool>("DoubleDisCo_binA_NonIsoMuon_0l");
        const auto& DoubleDisCo_binB_NonIsoMuon_0l    = tr.getVar<bool>("DoubleDisCo_binB_NonIsoMuon_0l");
        const auto& DoubleDisCo_binC_NonIsoMuon_0l    = tr.getVar<bool>("DoubleDisCo_binC_NonIsoMuon_0l");
        const auto& DoubleDisCo_binD_NonIsoMuon_0l    = tr.getVar<bool>("DoubleDisCo_binD_NonIsoMuon_0l");
        const auto& DoubleDisCo_binA_NonIsoMuon_1l    = tr.getVar<bool>("DoubleDisCo_binA_NonIsoMuon_1l");
        const auto& DoubleDisCo_binB_NonIsoMuon_1l    = tr.getVar<bool>("DoubleDisCo_binB_NonIsoMuon_1l");
        const auto& DoubleDisCo_binC_NonIsoMuon_1l    = tr.getVar<bool>("DoubleDisCo_binC_NonIsoMuon_1l");
        const auto& DoubleDisCo_binD_NonIsoMuon_1l    = tr.getVar<bool>("DoubleDisCo_binD_NonIsoMuon_1l");
        const auto& DoubleDisCo_massReg_0l            = tr.getVar<double>("DoubleDisCo_massReg_0l");
        const auto& DoubleDisCo_massReg_1l            = tr.getVar<double>("DoubleDisCo_massReg_1l");
        const auto& DoubleDisCo_massReg_NonIsoMuon_0l = tr.getVar<double>("DoubleDisCo_massReg_NonIsoMuon_0l");
        const auto& DoubleDisCo_massReg_NonIsoMuon_1l = tr.getVar<double>("DoubleDisCo_massReg_NonIsoMuon_1l");
        const auto& DoubleDisCo_disc1_0l              = tr.getVar<double>("DoubleDisCo_disc1_0l");
        const auto& DoubleDisCo_disc2_0l              = tr.getVar<double>("DoubleDisCo_disc2_0l");
        const auto& DoubleDisCo_disc1_1l              = tr.getVar<double>("DoubleDisCo_disc1_1l");
        const auto& DoubleDisCo_disc2_1l              = tr.getVar<double>("DoubleDisCo_disc2_1l");
        const auto& DoubleDisCo_disc1_NonIsoMuon_0l   = tr.getVar<double>("DoubleDisCo_disc1_NonIsoMuon_0l");
        const auto& DoubleDisCo_disc2_NonIsoMuon_0l   = tr.getVar<double>("DoubleDisCo_disc2_NonIsoMuon_0l");
        const auto& DoubleDisCo_disc1_NonIsoMuon_1l   = tr.getVar<double>("DoubleDisCo_disc1_NonIsoMuon_1l");
        const auto& DoubleDisCo_disc2_NonIsoMuon_1l   = tr.getVar<double>("DoubleDisCo_disc2_NonIsoMuon_1l");

        const auto& fwm2_top6_0l                      = tr.getVar<double>("fwm2_top6_0l");
        const auto& fwm3_top6_0l                      = tr.getVar<double>("fwm3_top6_0l");
        const auto& fwm4_top6_0l                      = tr.getVar<double>("fwm4_top6_0l");
        const auto& fwm5_top6_0l                      = tr.getVar<double>("fwm5_top6_0l");
        const auto& jmt_ev0_top6_0l                   = tr.getVar<double>("jmt_ev0_top6_0l");
        const auto& jmt_ev1_top6_0l                   = tr.getVar<double>("jmt_ev1_top6_0l");
        const auto& jmt_ev2_top6_0l                   = tr.getVar<double>("jmt_ev2_top6_0l");
        const auto& fwm2_top6_1l                      = tr.getVar<double>("fwm2_top6_1l");
        const auto& fwm3_top6_1l                      = tr.getVar<double>("fwm3_top6_1l");
        const auto& fwm4_top6_1l                      = tr.getVar<double>("fwm4_top6_1l");
        const auto& fwm5_top6_1l                      = tr.getVar<double>("fwm5_top6_1l");
        const auto& jmt_ev0_top6_1l                   = tr.getVar<double>("jmt_ev0_top6_1l");
        const auto& jmt_ev1_top6_1l                   = tr.getVar<double>("jmt_ev1_top6_1l");
        const auto& jmt_ev2_top6_1l                   = tr.getVar<double>("jmt_ev2_top6_1l");
        const auto& NonIsoMuons_fwm2_top6_0l          = tr.getVar<double>("NonIsoMuons_fwm2_top6_0l");
        const auto& NonIsoMuons_fwm3_top6_0l          = tr.getVar<double>("NonIsoMuons_fwm3_top6_0l");
        const auto& NonIsoMuons_fwm4_top6_0l          = tr.getVar<double>("NonIsoMuons_fwm4_top6_0l");
        const auto& NonIsoMuons_fwm5_top6_0l          = tr.getVar<double>("NonIsoMuons_fwm5_top6_0l");
        const auto& NonIsoMuons_jmt_ev0_top6_0l       = tr.getVar<double>("NonIsoMuons_jmt_ev0_top6_0l");
        const auto& NonIsoMuons_jmt_ev1_top6_0l       = tr.getVar<double>("NonIsoMuons_jmt_ev1_top6_0l");
        const auto& NonIsoMuons_jmt_ev2_top6_0l       = tr.getVar<double>("NonIsoMuons_jmt_ev2_top6_0l");
        const auto& NonIsoMuons_fwm2_top6_1l          = tr.getVar<double>("NonIsoMuons_fwm2_top6_1l");
        const auto& NonIsoMuons_fwm3_top6_1l          = tr.getVar<double>("NonIsoMuons_fwm3_top6_1l");
        const auto& NonIsoMuons_fwm4_top6_1l          = tr.getVar<double>("NonIsoMuons_fwm4_top6_1l");
        const auto& NonIsoMuons_fwm5_top6_1l          = tr.getVar<double>("NonIsoMuons_fwm5_top6_1l");
        const auto& NonIsoMuons_jmt_ev0_top6_1l       = tr.getVar<double>("NonIsoMuons_jmt_ev0_top6_1l");
        const auto& NonIsoMuons_jmt_ev1_top6_1l       = tr.getVar<double>("NonIsoMuons_jmt_ev1_top6_1l");
        const auto& NonIsoMuons_jmt_ev2_top6_1l       = tr.getVar<double>("NonIsoMuons_jmt_ev2_top6_1l");

        const auto& Jets_cm_top6_0l                   = tr.getVec<TLorentzVector>("Jets_cm_top6_0l");
        const auto& Jets_cm_top6_1l                   = tr.getVec<TLorentzVector>("Jets_cm_top6_1l");
        const auto& NonIsoMuons_Jets_cm_top6_0l       = tr.getVec<TLorentzVector>("NonIsoMuons_Jets_cm_top6_0l");
        const auto& NonIsoMuons_Jets_cm_top6_1l       = tr.getVec<TLorentzVector>("NonIsoMuons_Jets_cm_top6_1l");
        const auto& nMVAJets_0l                       = tr.getVar<unsigned int>("nMVAJets_0l");
        const auto& nMVAJets_1l                       = tr.getVar<unsigned int>("nMVAJets_1l");

        const auto& eventCounter                      = tr.getVar<int>("eventCounter");
        const auto& Stop1_pt_cm_OldSeed               = tr.getVar<double>("Stop1_pt_cm_OldSeed");
        const auto& Stop1_eta_cm_OldSeed              = tr.getVar<double>("Stop1_eta_cm_OldSeed");
        const auto& Stop1_phi_cm_OldSeed              = tr.getVar<double>("Stop1_phi_cm_OldSeed");
        const auto& Stop1_mass_cm_OldSeed             = tr.getVar<double>("Stop1_mass_cm_OldSeed");
        const auto& Stop2_pt_cm_OldSeed               = tr.getVar<double>("Stop2_pt_cm_OldSeed");
        const auto& Stop2_eta_cm_OldSeed              = tr.getVar<double>("Stop2_eta_cm_OldSeed");
        const auto& Stop2_phi_cm_OldSeed              = tr.getVar<double>("Stop2_phi_cm_OldSeed");
        const auto& Stop2_mass_cm_OldSeed             = tr.getVar<double>("Stop2_mass_cm_OldSeed");

        const auto& Stop1_pt_cm_NonIsoMuon_OldSeed    = tr.getVar<double>("Stop1_pt_cm_NonIsoMuon_OldSeed");
        const auto& Stop1_eta_cm_NonIsoMuon_OldSeed   = tr.getVar<double>("Stop1_eta_cm_NonIsoMuon_OldSeed");
        const auto& Stop1_phi_cm_NonIsoMuon_OldSeed   = tr.getVar<double>("Stop1_phi_cm_NonIsoMuon_OldSeed");
        const auto& Stop1_mass_cm_NonIsoMuon_OldSeed  = tr.getVar<double>("Stop1_mass_cm_NonIsoMuon_OldSeed");
        const auto& Stop2_pt_cm_NonIsoMuon_OldSeed    = tr.getVar<double>("Stop2_pt_cm_NonIsoMuon_OldSeed");
        const auto& Stop2_eta_cm_NonIsoMuon_OldSeed   = tr.getVar<double>("Stop2_eta_cm_NonIsoMuon_OldSeed");
        const auto& Stop2_phi_cm_NonIsoMuon_OldSeed   = tr.getVar<double>("Stop2_phi_cm_NonIsoMuon_OldSeed");
        const auto& Stop2_mass_cm_NonIsoMuon_OldSeed  = tr.getVar<double>("Stop2_mass_cm_NonIsoMuon_OldSeed");

        std::vector<unsigned int> nMVAJets{nMVAJets_0l, nMVAJets_1l};
        std::vector<int>          NGoodJets{NGoodJets_pt45, NGoodJets_pt30};
        std::vector<int>          NNonIsoMuonJets{NNonIsoMuonJets_pt45, NNonIsoMuonJets_pt30};
        std::vector<double>       HT_trigger{HT_trigger_pt45, HT_trigger_pt30};
        std::vector<double>       HT_NonIsoMuon{HT_NonIsoMuon_pt45, HT_NonIsoMuon_pt30};
        std::vector<bool>         DoubleDisCo_binA{DoubleDisCo_binA_0l,    DoubleDisCo_binA_1l};
        std::vector<bool>         DoubleDisCo_binB{DoubleDisCo_binB_0l,    DoubleDisCo_binB_1l};
        std::vector<bool>         DoubleDisCo_binC{DoubleDisCo_binC_0l,    DoubleDisCo_binC_1l};
        std::vector<bool>         DoubleDisCo_binD{DoubleDisCo_binD_0l,    DoubleDisCo_binD_1l};
        std::vector<double>       DoubleDisCo_massReg{DoubleDisCo_massReg_0l, DoubleDisCo_massReg_1l};
        std::vector<double>       DoubleDisCo_disc1{DoubleDisCo_disc1_0l, DoubleDisCo_disc1_1l};
        std::vector<double>       DoubleDisCo_disc2{DoubleDisCo_disc2_0l, DoubleDisCo_disc2_1l};
        std::vector<bool>         DoubleDisCo_QCDCR_binA{DoubleDisCo_binA_NonIsoMuon_0l,    DoubleDisCo_binA_NonIsoMuon_1l};
        std::vector<bool>         DoubleDisCo_QCDCR_binB{DoubleDisCo_binB_NonIsoMuon_0l,    DoubleDisCo_binB_NonIsoMuon_1l};
        std::vector<bool>         DoubleDisCo_QCDCR_binC{DoubleDisCo_binC_NonIsoMuon_0l,    DoubleDisCo_binC_NonIsoMuon_1l};
        std::vector<bool>         DoubleDisCo_QCDCR_binD{DoubleDisCo_binD_NonIsoMuon_0l,    DoubleDisCo_binD_NonIsoMuon_1l};
        std::vector<double>       DoubleDisCo_QCDCR_massReg{DoubleDisCo_massReg_NonIsoMuon_0l, DoubleDisCo_massReg_NonIsoMuon_1l};
        std::vector<double>       DoubleDisCo_QCDCR_disc1{DoubleDisCo_disc1_NonIsoMuon_0l, DoubleDisCo_disc1_NonIsoMuon_1l};
        std::vector<double>       DoubleDisCo_QCDCR_disc2{DoubleDisCo_disc2_NonIsoMuon_0l, DoubleDisCo_disc2_NonIsoMuon_1l};

        std::vector<double>       fwm2_top6{fwm2_top6_0l,       fwm2_top6_1l};
        std::vector<double>       fwm3_top6{fwm3_top6_0l,       fwm3_top6_1l};
        std::vector<double>       fwm4_top6{fwm4_top6_0l,       fwm4_top6_1l};
        std::vector<double>       fwm5_top6{fwm5_top6_0l,       fwm5_top6_1l};
        std::vector<double>       jmt_ev0_top6{jmt_ev0_top6_0l, jmt_ev0_top6_1l};
        std::vector<double>       jmt_ev1_top6{jmt_ev1_top6_0l, jmt_ev1_top6_1l};
        std::vector<double>       jmt_ev2_top6{jmt_ev2_top6_0l, jmt_ev2_top6_1l};

        std::vector<double>       fwm2_top6_QCDCR{NonIsoMuons_fwm2_top6_0l,       NonIsoMuons_fwm2_top6_1l};
        std::vector<double>       fwm3_top6_QCDCR{NonIsoMuons_fwm3_top6_0l,       NonIsoMuons_fwm3_top6_1l};
        std::vector<double>       fwm4_top6_QCDCR{NonIsoMuons_fwm4_top6_0l,       NonIsoMuons_fwm4_top6_1l};
        std::vector<double>       fwm5_top6_QCDCR{NonIsoMuons_fwm5_top6_0l,       NonIsoMuons_fwm5_top6_1l};
        std::vector<double>       jmt_ev0_top6_QCDCR{NonIsoMuons_jmt_ev0_top6_0l, NonIsoMuons_jmt_ev0_top6_1l};
        std::vector<double>       jmt_ev1_top6_QCDCR{NonIsoMuons_jmt_ev1_top6_0l, NonIsoMuons_jmt_ev1_top6_1l};
        std::vector<double>       jmt_ev2_top6_QCDCR{NonIsoMuons_jmt_ev2_top6_0l, NonIsoMuons_jmt_ev2_top6_1l};

        std::vector<std::vector<TLorentzVector>> Jets_cm_top6{Jets_cm_top6_0l, Jets_cm_top6_1l};
        std::vector<std::vector<TLorentzVector>> Jets_cm_top6_QCDCR{NonIsoMuons_Jets_cm_top6_0l, NonIsoMuons_Jets_cm_top6_1l};

        const auto& met                    = tr.getVar<double>("MET");
        const auto& metPhi                 = tr.getVar<double>("METPhi");

        std::vector<double> Jets_flavb_0l;   std::vector<double> Jets_flavb_1l;
        std::vector<double> Jets_flavc_0l;   std::vector<double> Jets_flavc_1l;
        std::vector<double> Jets_flavg_0l;   std::vector<double> Jets_flavg_1l;
        std::vector<double> Jets_flavuds_0l; std::vector<double> Jets_flavuds_1l;
        std::vector<double> Jets_flavq_0l;   std::vector<double> Jets_flavq_1l;

        std::vector<double> JetNonIsoMuons_flavb_0l;   std::vector<double> JetNonIsoMuons_flavb_1l;
        std::vector<double> JetNonIsoMuons_flavc_0l;   std::vector<double> JetNonIsoMuons_flavc_1l;
        std::vector<double> JetNonIsoMuons_flavg_0l;   std::vector<double> JetNonIsoMuons_flavg_1l;
        std::vector<double> JetNonIsoMuons_flavuds_0l; std::vector<double> JetNonIsoMuons_flavuds_1l;
        std::vector<double> JetNonIsoMuons_flavq_0l;   std::vector<double> JetNonIsoMuons_flavq_1l;

        // Here assume number cm jets is same in CR and SR selection
        for (unsigned int iJet = 1; iJet <= nMVAJets_0l; iJet++) {
            Jets_flavb_0l.push_back(tr.getVar<double>("Jet_flavb_"+std::to_string(iJet)+"_0l"));
            Jets_flavc_0l.push_back(tr.getVar<double>("Jet_flavc_"+std::to_string(iJet)+"_0l"));
            Jets_flavg_0l.push_back(tr.getVar<double>("Jet_flavg_"+std::to_string(iJet)+"_0l"));
            Jets_flavuds_0l.push_back(tr.getVar<double>("Jet_flavuds_"+std::to_string(iJet)+"_0l"));
            Jets_flavq_0l.push_back(tr.getVar<double>("Jet_flavq_"+std::to_string(iJet)+"_0l"));

            JetNonIsoMuons_flavb_0l.push_back(tr.getVar<double>("JetNonIsoMuons_flavb_"+std::to_string(iJet)+"_0l"));
            JetNonIsoMuons_flavc_0l.push_back(tr.getVar<double>("JetNonIsoMuons_flavc_"+std::to_string(iJet)+"_0l"));
            JetNonIsoMuons_flavg_0l.push_back(tr.getVar<double>("JetNonIsoMuons_flavg_"+std::to_string(iJet)+"_0l"));
            JetNonIsoMuons_flavuds_0l.push_back(tr.getVar<double>("JetNonIsoMuons_flavuds_"+std::to_string(iJet)+"_0l"));
            JetNonIsoMuons_flavq_0l.push_back(tr.getVar<double>("JetNonIsoMuons_flavq_"+std::to_string(iJet)+"_0l"));
        }

        for (unsigned int iJet = 1; iJet <= nMVAJets_1l; iJet++) {
            Jets_flavb_1l.push_back(tr.getVar<double>("Jet_flavb_"+std::to_string(iJet)+"_1l"));
            Jets_flavc_1l.push_back(tr.getVar<double>("Jet_flavc_"+std::to_string(iJet)+"_1l"));
            Jets_flavg_1l.push_back(tr.getVar<double>("Jet_flavg_"+std::to_string(iJet)+"_1l"));
            Jets_flavuds_1l.push_back(tr.getVar<double>("Jet_flavuds_"+std::to_string(iJet)+"_1l"));
            Jets_flavq_1l.push_back(tr.getVar<double>("Jet_flavq_"+std::to_string(iJet)+"_1l"));

            JetNonIsoMuons_flavb_1l.push_back(tr.getVar<double>("JetNonIsoMuons_flavb_"+std::to_string(iJet)+"_1l"));
            JetNonIsoMuons_flavc_1l.push_back(tr.getVar<double>("JetNonIsoMuons_flavc_"+std::to_string(iJet)+"_1l"));
            JetNonIsoMuons_flavg_1l.push_back(tr.getVar<double>("JetNonIsoMuons_flavg_"+std::to_string(iJet)+"_1l"));
            JetNonIsoMuons_flavuds_1l.push_back(tr.getVar<double>("JetNonIsoMuons_flavuds_"+std::to_string(iJet)+"_1l"));
            JetNonIsoMuons_flavq_1l.push_back(tr.getVar<double>("JetNonIsoMuons_flavq_"+std::to_string(iJet)+"_1l"));
        }

        std::vector<std::vector<double>> Jets_flavb{Jets_flavb_0l,     Jets_flavb_1l};
        std::vector<std::vector<double>> Jets_flavc{Jets_flavc_0l,     Jets_flavc_1l};
        std::vector<std::vector<double>> Jets_flavg{Jets_flavg_0l,     Jets_flavg_1l};
        std::vector<std::vector<double>> Jets_flavuds{Jets_flavuds_0l, Jets_flavuds_1l};
        std::vector<std::vector<double>> Jets_flavq{Jets_flavq_0l,     Jets_flavq_1l};

        std::vector<std::vector<double>> Jets_flavb_QCDCR{JetNonIsoMuons_flavb_0l,     JetNonIsoMuons_flavb_1l};
        std::vector<std::vector<double>> Jets_flavc_QCDCR{JetNonIsoMuons_flavc_0l,     JetNonIsoMuons_flavc_1l};
        std::vector<std::vector<double>> Jets_flavg_QCDCR{JetNonIsoMuons_flavg_0l,     JetNonIsoMuons_flavg_1l};
        std::vector<std::vector<double>> Jets_flavuds_QCDCR{JetNonIsoMuons_flavuds_0l, JetNonIsoMuons_flavuds_1l};
        std::vector<std::vector<double>> Jets_flavq_QCDCR{JetNonIsoMuons_flavq_0l,     JetNonIsoMuons_flavq_1l};

        // ------------------------
        // -- Print event number
        // ------------------------       
        if(maxevents != -1 && tr.getEvtNum() >= maxevents) break;        
        if(tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() );

        // ------------------------
        // -- Define weight
        // ------------------------
        double weight1L=1.0, eventweight=1.0, leptonweight=1.0, bTagWeight=1.0, prefiringScaleFactor=1.0, pileupWeight=1.0, htDerivedweight=1.0;
        double weight0L=1.0, weight0L_NonIsoMuon=1.0, weight1L_NonIsoMuon=1.0;
        if(runtype == "MC")
        {
            // Define Lumi weight
            const auto& Weight  = tr.getVar<double>("Weight");
            const auto& lumi = tr.getVar<double>("Lumi");
            eventweight = lumi*Weight;
           
            const auto& eleLepWeight = tr.getVar<double>("totGoodElectronSF");
            const auto& muLepWeight  = tr.getVar<double>("totGoodMuonSF");
            const auto& muNonIso     = tr.getVar<double>("totNonIsoMuonSF");
            leptonweight = eleLepWeight*muLepWeight;
          
            pileupWeight = tr.getVar<double>("puWeightCorr");
            bTagWeight   = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            htDerivedweight = tr.getVar<double>("htDerivedweight");
            prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor");
            
            weight1L            *= eventweight*leptonweight*bTagWeight*prefiringScaleFactor*pileupWeight*htDerivedweight;
            weight0L            *= eventweight*bTagWeight*prefiringScaleFactor*pileupWeight;
            weight1L_NonIsoMuon *= eventweight*muNonIso*prefiringScaleFactor*pileupWeight;
            weight0L_NonIsoMuon *= eventweight*muNonIso*prefiringScaleFactor*pileupWeight;

        }

        std::vector<double> weight{weight0L, weight1L};
        std::vector<double> weight_QCDCR{weight0L_NonIsoMuon, weight1L_NonIsoMuon};

        const std::map<std::string, bool> cut_map 
        {
            {"_1l"       , passBaseline1l_Good},                         
            {"_0l"       , passBaseline0l_Good},                         
            {"_1l_QCDCR" , passBaseline1l_NonIsoMuon},                         
            {"_0l_QCDCR" , passBaseline0l_NonIsoMuon},                         

        };

        // Initialize Histograms
        if(!initHistos)
        {
            Preinit(nMVAJets_0l > nMVAJets_1l ? nMVAJets_0l : nMVAJets_1l);
            InitHistos(cut_map);
            initHistos = true;
        }

        my_histos["EventCounter"]->Fill(eventCounter);

        std::map<std::string, bool> njetsMap;
        std::map<std::string, bool> ABCDmap;
        for(auto& kv : cut_map)
        {

            // Extract "0" or "1" from cut string e.g. _1l_QCDCR
            int channel = 1;
            if (kv.first.size() > 0 and kv.first.substr(2,1) == "l")
                channel = std::stoi(kv.first.substr(1,1));

            bool isQCD = kv.first.find("QCDCR") != std::string::npos;

            njetsMap = {{"Incl",     true},
                          {"7",      !isQCD ? NGoodJets[channel]==7  : NNonIsoMuonJets[channel]==7},
                          {"8",      !isQCD ? NGoodJets[channel]==8  : NNonIsoMuonJets[channel]==8},
                          {"9",      !isQCD ? NGoodJets[channel]==9  : NNonIsoMuonJets[channel]==9},
                          {"10",     !isQCD ? NGoodJets[channel]==10 : NNonIsoMuonJets[channel]==10},
                          {"11",     !isQCD ? NGoodJets[channel]==11 : NNonIsoMuonJets[channel]==11},
                          {"11incl", !isQCD ? NGoodJets[channel]>=11 : NNonIsoMuonJets[channel]>=11},
                          {"12",     !isQCD ? NGoodJets[channel]==12 : NNonIsoMuonJets[channel]==12},
                          {"12incl", !isQCD ? NGoodJets[channel]>=12 : NNonIsoMuonJets[channel]>=12}
            };

            bool passBinA = !isQCD ? DoubleDisCo_binA[channel] : DoubleDisCo_QCDCR_binA[channel];
            bool passBinB = !isQCD ? DoubleDisCo_binB[channel] : DoubleDisCo_QCDCR_binB[channel];
            bool passBinC = !isQCD ? DoubleDisCo_binC[channel] : DoubleDisCo_QCDCR_binC[channel];
            bool passBinD = !isQCD ? DoubleDisCo_binD[channel] : DoubleDisCo_QCDCR_binD[channel];

            int shift = -1;
            if (passBinA) shift = 0;
            if (passBinB) shift = 1;
            if (passBinC) shift = 2;
            if (passBinD) shift = 3;

            ABCDmap = {{"",  true},
                       {"A", passBinA},
                       {"B", passBinB},
                       {"C", passBinC},
                       {"D", passBinD}};

            for(auto& ipass : njetsMap)
            {
                std::string njetStr = "";
                if (ipass.first != "Incl") njetStr = "_Njets" + ipass.first;

                for(auto& jpass : ABCDmap)
                {

                    std::string regionStr = "";
                    if (jpass.first != "") regionStr = "_" + jpass.first;

                    if(kv.second and ipass.second and jpass.second)
                    {

                        std::string name = kv.first+njetStr+regionStr;
                        double w = !isQCD ? weight[channel] : weight_QCDCR[channel];

                        // if plotting Njets, don't care about individual Njet cuts
                        if (ipass.first == "Incl") {
                            my_histos["h_njets"                   + name]->Fill(!isQCD ? NGoodJets[channel] : NNonIsoMuonJets[channel], w);
                            if (jpass.first == "") {
                                if (!isQCD) {
                                    my_histos["h_njets_11incl" + name]->Fill(NGoodJets[channel]>=11       ? 11-nMVAJets[channel]+shift*5 : NGoodJets[channel]-nMVAJets[channel]+shift*5, w);
                                    my_histos["h_njets_12incl" + name]->Fill(NGoodJets[channel]>=12       ? 12-nMVAJets[channel]+shift*6 : NGoodJets[channel]-nMVAJets[channel]+shift*6, w);
                                } else {
                                    my_histos["h_njets_11incl" + name]->Fill(NNonIsoMuonJets[channel]>=11 ? 11-nMVAJets[channel]+shift*5 : NNonIsoMuonJets[channel]-nMVAJets[channel]+shift*5, w);
                                    my_histos["h_njets_12incl" + name]->Fill(NNonIsoMuonJets[channel]>=12 ? 12-nMVAJets[channel]+shift*6 : NNonIsoMuonJets[channel]-nMVAJets[channel]+shift*6, w);
               
                                }
                            }
                        }

                        // if plotting disco, no need to make plots when cutting on
                        if (jpass.first == "") {
                            if (!isQCD) {
                                my_histos["h_DoubleDisCo_disc1"          + name]->Fill(DoubleDisCo_disc1[channel], w);
                                my_histos["h_DoubleDisCo_disc2"          + name]->Fill(DoubleDisCo_disc2[channel], w);
                                my_2d_histos["h_DoubleDisCo_disc1_disc2" + name]->Fill(DoubleDisCo_disc1[channel], DoubleDisCo_disc2[channel], w);
                            } else {
                                my_histos["h_DoubleDisCo_disc1"          + name]->Fill(DoubleDisCo_QCDCR_disc1[channel], w);
                                my_histos["h_DoubleDisCo_disc2"          + name]->Fill(DoubleDisCo_QCDCR_disc2[channel], w);
                                my_2d_histos["h_DoubleDisCo_disc1_disc2" + name]->Fill(DoubleDisCo_QCDCR_disc1[channel], DoubleDisCo_QCDCR_disc2[channel], w);
                            }
                        }

                        my_histos["fwm2_top6"               + name]->Fill(!isQCD ? fwm2_top6[channel]    : fwm2_top6_QCDCR[channel], w);
                        my_histos["fwm3_top6"               + name]->Fill(!isQCD ? fwm3_top6[channel]    : fwm3_top6_QCDCR[channel], w);
                        my_histos["fwm4_top6"               + name]->Fill(!isQCD ? fwm4_top6[channel]    : fwm4_top6_QCDCR[channel], w);
                        my_histos["fwm5_top6"               + name]->Fill(!isQCD ? fwm5_top6[channel]    : fwm5_top6_QCDCR[channel], w);
                        my_histos["jmt_ev0_top6"            + name]->Fill(!isQCD ? jmt_ev0_top6[channel] : jmt_ev0_top6_QCDCR[channel], w);
                        my_histos["jmt_ev1_top6"            + name]->Fill(!isQCD ? jmt_ev1_top6[channel] : jmt_ev1_top6_QCDCR[channel], w);
                        my_histos["jmt_ev2_top6"            + name]->Fill(!isQCD ? jmt_ev2_top6[channel] : jmt_ev2_top6_QCDCR[channel], w);

                        if (Stop1_pt_cm_OldSeed  > Stop2_pt_cm_OldSeed) {
                            my_histos["Stop1_pt_cm_OldSeed"     + name]->Fill(!isQCD ? Stop1_pt_cm_OldSeed   : Stop1_pt_cm_NonIsoMuon_OldSeed,   w);
                            my_histos["Stop1_eta_cm_OldSeed"    + name]->Fill(!isQCD ? Stop1_eta_cm_OldSeed  : Stop1_eta_cm_NonIsoMuon_OldSeed,  w);
                            my_histos["Stop1_phi_cm_OldSeed"    + name]->Fill(!isQCD ? Stop1_phi_cm_OldSeed  : Stop1_phi_cm_NonIsoMuon_OldSeed,  w);
                            my_histos["Stop1_mass_cm_OldSeed"   + name]->Fill(!isQCD ? Stop1_mass_cm_OldSeed : Stop1_mass_cm_NonIsoMuon_OldSeed, w);
                    
                            my_histos["Stop2_pt_cm_OldSeed"     + name]->Fill(!isQCD ? Stop2_pt_cm_OldSeed   : Stop2_pt_cm_NonIsoMuon_OldSeed,   w);
                            my_histos["Stop2_eta_cm_OldSeed"    + name]->Fill(!isQCD ? Stop2_eta_cm_OldSeed  : Stop2_eta_cm_NonIsoMuon_OldSeed,  w);
                            my_histos["Stop2_phi_cm_OldSeed"    + name]->Fill(!isQCD ? Stop2_phi_cm_OldSeed  : Stop2_phi_cm_NonIsoMuon_OldSeed,  w);
                            my_histos["Stop2_mass_cm_OldSeed"   + name]->Fill(!isQCD ? Stop2_mass_cm_OldSeed : Stop2_mass_cm_NonIsoMuon_OldSeed, w);
                        } else {
                            my_histos["Stop1_pt_cm_OldSeed"     + name]->Fill(!isQCD ? Stop2_pt_cm_OldSeed   : Stop2_pt_cm_NonIsoMuon_OldSeed,   w);
                            my_histos["Stop1_eta_cm_OldSeed"    + name]->Fill(!isQCD ? Stop2_eta_cm_OldSeed  : Stop2_eta_cm_NonIsoMuon_OldSeed,  w);
                            my_histos["Stop1_phi_cm_OldSeed"    + name]->Fill(!isQCD ? Stop2_phi_cm_OldSeed  : Stop2_phi_cm_NonIsoMuon_OldSeed,  w);
                            my_histos["Stop1_mass_cm_OldSeed"   + name]->Fill(!isQCD ? Stop2_mass_cm_OldSeed : Stop2_mass_cm_NonIsoMuon_OldSeed, w);
                    
                            my_histos["Stop2_pt_cm_OldSeed"     + name]->Fill(!isQCD ? Stop1_pt_cm_OldSeed   : Stop1_pt_cm_NonIsoMuon_OldSeed,   w);
                            my_histos["Stop2_eta_cm_OldSeed"    + name]->Fill(!isQCD ? Stop1_eta_cm_OldSeed  : Stop1_eta_cm_NonIsoMuon_OldSeed,  w);
                            my_histos["Stop2_phi_cm_OldSeed"    + name]->Fill(!isQCD ? Stop1_phi_cm_OldSeed  : Stop1_phi_cm_NonIsoMuon_OldSeed,  w);
                            my_histos["Stop2_mass_cm_OldSeed"   + name]->Fill(!isQCD ? Stop1_mass_cm_OldSeed : Stop1_mass_cm_NonIsoMuon_OldSeed, w);
                        }

                        unsigned int nJets = !isQCD ? Jets_cm_top6[channel].size() : Jets_cm_top6_QCDCR[channel].size();
                        for(unsigned int i = 1; i <= nJets; i++)
                        {
                            double pt=0.0, eta=0.0, phi=0.0, m=0.0, E=0.0;

                            if (!isQCD) {
                                pt  = static_cast<double>(Jets_cm_top6[channel].at(i-1).Pt());
                                eta = static_cast<double>(Jets_cm_top6[channel].at(i-1).Eta());
                                phi = static_cast<double>(Jets_cm_top6[channel].at(i-1).Phi());
                                m   = static_cast<double>(Jets_cm_top6[channel].at(i-1).M());
                                E   = static_cast<double>(Jets_cm_top6[channel].at(i-1).E());
                            } else {
                                pt  = static_cast<double>(Jets_cm_top6_QCDCR[channel].at(i-1).Pt());
                                eta = static_cast<double>(Jets_cm_top6_QCDCR[channel].at(i-1).Eta());
                                phi = static_cast<double>(Jets_cm_top6_QCDCR[channel].at(i-1).Phi());
                                m   = static_cast<double>(Jets_cm_top6_QCDCR[channel].at(i-1).M());
                                E   = static_cast<double>(Jets_cm_top6_QCDCR[channel].at(i-1).E());
                            } 
   
                            my_histos["Jet_cm_pt_"           + std::to_string(i) + name]->Fill(pt, w);
                            my_histos["Jet_cm_eta_"          + std::to_string(i) + name]->Fill(eta, w);
                            my_histos["Jet_cm_phi_"          + std::to_string(i) + name]->Fill(phi, w);
                            my_histos["Jet_cm_m_"            + std::to_string(i) + name]->Fill(m, w);
                            my_histos["Jet_cm_E_"            + std::to_string(i) + name]->Fill(E, w);

                            if (!isQCD) {
                                my_histos["Jet_cm_flavb_"   + std::to_string(i) + name]->Fill(Jets_flavb[channel].at(i-1), w);
                                my_histos["Jet_cm_flavc_"   + std::to_string(i) + name]->Fill(Jets_flavc[channel].at(i-1), w);
                                my_histos["Jet_cm_flavg_"   + std::to_string(i) + name]->Fill(Jets_flavg[channel].at(i-1), w);
                                my_histos["Jet_cm_flavq_"   + std::to_string(i) + name]->Fill(Jets_flavq[channel].at(i-1), w);
                                my_histos["Jet_cm_flavuds_" + std::to_string(i) + name]->Fill(Jets_flavuds[channel].at(i-1), w);
                            } else {
                                my_histos["Jet_cm_flavb_"   + std::to_string(i) + name]->Fill(Jets_flavb_QCDCR[channel].at(i-1), w);
                                my_histos["Jet_cm_flavc_"   + std::to_string(i) + name]->Fill(Jets_flavc_QCDCR[channel].at(i-1), w);
                                my_histos["Jet_cm_flavg_"   + std::to_string(i) + name]->Fill(Jets_flavg_QCDCR[channel].at(i-1), w);
                                my_histos["Jet_cm_flavq_"   + std::to_string(i) + name]->Fill(Jets_flavq_QCDCR[channel].at(i-1), w);
                                my_histos["Jet_cm_flavuds_" + std::to_string(i) + name]->Fill(Jets_flavuds_QCDCR[channel].at(i-1), w);
                            }
                        }

                        my_histos["h_ntops"                   + name]->Fill(ntops, w);
                        my_histos["h_ht"                      + name]->Fill(!isQCD ? HT_trigger[channel] : HT_NonIsoMuon[channel], w);
                        my_histos["h_MET_phi"                 + name]->Fill(metPhi, w);
                        my_histos["h_MET_pt"                  + name]->Fill(met, w);
                        my_histos["h_DoubleDisCo_massReg"     + name]->Fill(!isQCD ? DoubleDisCo_massReg[channel] : DoubleDisCo_QCDCR_massReg[channel], w);

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
    {
        p.second->Write();
    }
    
    for(const auto& p : my_2d_histos) 
    {
        p.second->Write();
    }
}
