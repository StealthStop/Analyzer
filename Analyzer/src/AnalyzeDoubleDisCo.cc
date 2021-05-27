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
        {"h_DoubleDisCo_disc1",   80,  0, 1},
        {"h_DoubleDisCo_disc2",   80,  0, 1},
        {"h_DoubleDisCo_massReg", 150, 0, 1500},
        {"h_DoubleDisCo_disc1",   80,  0, 1},
        {"h_DoubleDisCo_disc2",   80,  0, 1},
        {"h_DoubleDisCo_massReg", 150, 0, 1500},
        {"fwm2_top6",             50,  0, 1},
        {"fwm3_top6",             50,  0, 1},
        {"fwm4_top6",             50,  0, 1},
        {"fwm5_top6",             50,  0, 1},
        {"fwm6_top6",             50,  0, 1},
        {"fwm7_top6",             50,  0, 1},
        {"fwm8_top6",             50,  0, 1},
        {"fwm9_top6",             50,  0, 1},
        {"fwm10_top6",            50,  0, 1},
        {"jmt_ev0_top6",          50,  0, 1},
        {"jmt_ev1_top6",          50,  0, 1},
        {"jmt_ev2_top6",          50,  0, 1},
        {"GoodLeptons_pt_1",      150, 0, 1500},
        {"GoodLeptons_eta_1",     100, -6, 6},
        {"GoodLeptons_phi_1",     80,  -4, 4},
        {"GoodLeptons_m_1",       20,  0, 200},
        {"Stop1_pt_cm_OldSeed",      72, 0, 1500},
        {"Stop1_eta_cm_OldSeed",     80, -6, 6},
        {"Stop1_phi_cm_OldSeed",     64,  -4, 4},
        {"Stop1_mass_cm_OldSeed",       72,  0, 1500},
        {"Stop2_pt_cm_OldSeed",      72, 0, 1500},
        {"Stop2_eta_cm_OldSeed",     80, -6, 6},
        {"Stop2_phi_cm_OldSeed",     64,  -4, 4},
        {"Stop2_mass_cm_OldSeed",       72,  0, 1500},
        {    "h_njets",           20,   0.0,   20.0},
        {"blind_njets",           20,   0.0,   20.0},
        {    "h_ntops",           10,   0.0,   10.0},
        {"blind_ntops",           10,   0.0,   10.0},
        {    "h_nb",              10,   0.0,   10.0},
        {"blind_nb",              10,   0.0,   10.0},
        {    "h_ht",              500,   0.0, 5000.0},
        {"blind_ht",              500,   0.0, 5000.0},
        {    "h_mbl",             300,   0.0,  300.0},
        {"blind_mbl",             300,   0.0,  300.0},
        {    "h_lPt",             150,   0.0, 1500.0},
        {"blind_lPt",             150,   0.0, 1500.0},
        {    "h_lEta",            100,  -6.0,    6.0},
        {"blind_lEta",            100,  -6.0,    6.0},
        {    "h_lPhi",            80,  -4.0,    4.0},
        {"blind_lPhi",            80,  -4.0,    4.0},
        {    "h_jPt",             150,   0.0, 1500.0},
        {"blind_jPt",             150,   0.0, 1500.0},
        {    "h_jEta",            100,  -6.0,    6.0},
        {"blind_jEta",            100,  -6.0,    6.0},
        {    "h_jPhi",            80,  -4.0,    4.0},
        {"blind_jPhi",            80,  -4.0,    4.0},
    };

    hist2DInfos = {
        {"h_DoubleDisCo_disc1_disc2",  80,    0,    1, 80,     0,     1}, 
        {"h_DoubleDisCo_disc1_disc2",  80,    0,    1, 80,     0,     1}, 
        {    "h_lEta_lPhi",         100, -6.0,  6.0, 80,  -4.0,   4.0},
        {"blind_lEta_lPhi",         100, -6.0,  6.0, 80,  -4.0,   4.0},
        {    "h_jEta_jPhi",         100, -6.0,  6.0, 80,  -4.0,   4.0},
        {"blind_jEta_jPhi",         100, -6.0,  6.0, 80,  -4.0,   4.0},
    };

    abcds = {"", "A", "B", "C", "D"};
    njets = {"Incl", "6", "7", "8", "9", "10", "11", "11incl", "12", "12incl", "13", "14"};

}

void AnalyzeDoubleDisCo::Preinit(unsigned int nNNJets)
{
    for(unsigned int i = 1; i <= nNNJets ; i++)
    {
        histInfos.push_back({"Jet_cm_pt_"           + std::to_string(i), 150, 0, 1500});
        histInfos.push_back({"Jet_cm_eta_"          + std::to_string(i), 100, -6, 6});
        histInfos.push_back({"Jet_cm_phi_"          + std::to_string(i), 80, -4, 4});
        histInfos.push_back({"Jet_cm_m_"            + std::to_string(i), 20, 0, 200});
        histInfos.push_back({"Jet_cm_ptD_"          + std::to_string(i), 80, 0, 1});
        histInfos.push_back({"Jet_cm_axismajor_"    + std::to_string(i), 80, 0, 0.4});
        histInfos.push_back({"Jet_cm_axisminor_"    + std::to_string(i), 80, 0, 0.4});
        histInfos.push_back({"Jet_cm_multiplicity_" + std::to_string(i), 121, -0.5, 120.5});
        histInfos.push_back({"Jet_cm_flavb_"        + std::to_string(i), 80, 0, 1});
        histInfos.push_back({"Jet_cm_flavc_"        + std::to_string(i), 80, 0, 1});
        histInfos.push_back({"Jet_cm_flavg_"        + std::to_string(i), 80, 0, 1});
        histInfos.push_back({"Jet_cm_flavuds_"        + std::to_string(i), 80, 0, 1});

    }
}

void AnalyzeDoubleDisCo::InitHistos(const std::map<std::string, bool>& cutMap)
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    my_histos.emplace("EventCounter", std::make_shared<TH1D>("EventCounter","EventCounter", 2, -1.1, 1.1 ) );

    for(auto& mycut : cutMap)
    {
        for(const auto& hInfo : histInfos)
        { 
            for(const auto& Njet : njets)
            {
                std::string njetStr = "";
                if (Njet != "Incl") njetStr = "_Njets" + Njet;

                for(const auto& region : abcds)
                {
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
        const auto& ntops                  = tr.getVar<int>("ntops");
        const auto& runtype                = tr.getVar<std::string>("runtype");     

        const auto& Jets                   = tr.getVec<TLorentzVector>("Jets");
        const auto& GoodJets_pt30          = tr.getVec<bool>("GoodJets_pt30");
        const auto& GoodJets_pt45          = tr.getVec<bool>("GoodJets_pt45");
        const auto& NGoodJets_pt30         = tr.getVar<int>("NGoodJets_pt30");
        const auto& NGoodJets_pt45         = tr.getVar<int>("NGoodJets_pt45");
        const auto& NGoodBJets_pt30        = tr.getVar<int>("NGoodBJets_pt30");
        const auto& NGoodBJets_pt45        = tr.getVar<int>("NGoodBJets_pt45");

        std::vector<std::vector<bool>> GoodJets{GoodJets_pt45, GoodJets_pt30};

        const std::vector<int> NGoodJets{NGoodJets_pt45, NGoodJets_pt30};
        const std::vector<int> NGoodBJets{NGoodBJets_pt45, NGoodBJets_pt30};

        const auto& GoodLeptons            = tr.getVec<std::pair<std::string, TLorentzVector>>("GoodLeptons");

        const auto& HT_trigger_pt30        = tr.getVar<double>("HT_trigger_pt30");
        const auto& HT_trigger_pt45        = tr.getVar<double>("HT_trigger_pt45");

        const std::vector<double> HT_trigger{HT_trigger_pt45, HT_trigger_pt30};

        const auto& correct2018Split       = tr.getVar<bool>("correct2018Split");
        const auto& passTrigger            = tr.getVar<bool>("passTrigger");
        const auto& passTriggerMC          = tr.getVar<bool>("passTriggerMC");
        const auto& passMETFilters         = tr.getVar<bool>("passMETFilters");
        const auto& passMadHT              = tr.getVar<bool>("passMadHT");
        const auto& passHEMVeto            = tr.getVar<bool>("passHEMVeto");

        const auto& Mbl                    = tr.getVar<double>("Mbl");
        const auto& passBlind              = tr.getVar<bool>("passBlindLep_Good");            
        const auto& passBaseline0l_Good    = tr.getVar<bool>("passBaseline0l_Good");
        const auto& passBaseline1l_Good    = tr.getVar<bool>("passBaseline1l_Good");

        const auto& DoubleDisCo_binA_0l    = tr.getVar<bool>("DoubleDisCo_binA_0l");
        const auto& DoubleDisCo_binB_0l    = tr.getVar<bool>("DoubleDisCo_binB_0l");
        const auto& DoubleDisCo_binC_0l    = tr.getVar<bool>("DoubleDisCo_binC_0l");
        const auto& DoubleDisCo_binD_0l    = tr.getVar<bool>("DoubleDisCo_binD_0l");
        const auto& DoubleDisCo_binA_1l    = tr.getVar<bool>("DoubleDisCo_binA_1l");
        const auto& DoubleDisCo_binB_1l    = tr.getVar<bool>("DoubleDisCo_binB_1l");
        const auto& DoubleDisCo_binC_1l    = tr.getVar<bool>("DoubleDisCo_binC_1l");
        const auto& DoubleDisCo_binD_1l    = tr.getVar<bool>("DoubleDisCo_binD_1l");
        const auto& DoubleDisCo_massReg_0l = tr.getVar<double>("DoubleDisCo_massReg_0l");
        const auto& DoubleDisCo_massReg_1l = tr.getVar<double>("DoubleDisCo_massReg_1l");
        const auto& DoubleDisCo_disc1_0l   = tr.getVar<double>("DoubleDisCo_disc1_0l");
        const auto& DoubleDisCo_disc2_0l   = tr.getVar<double>("DoubleDisCo_disc2_0l");
        const auto& DoubleDisCo_disc1_1l   = tr.getVar<double>("DoubleDisCo_disc1_1l");
        const auto& DoubleDisCo_disc2_1l   = tr.getVar<double>("DoubleDisCo_disc2_1l");

        const std::vector<bool>   DoubleDisCo_binA{DoubleDisCo_binA_0l,    DoubleDisCo_binA_1l};
        const std::vector<bool>   DoubleDisCo_binB{DoubleDisCo_binB_0l,    DoubleDisCo_binB_1l};
        const std::vector<bool>   DoubleDisCo_binC{DoubleDisCo_binC_0l,    DoubleDisCo_binC_1l};
        const std::vector<bool>   DoubleDisCo_binD{DoubleDisCo_binD_0l,    DoubleDisCo_binD_1l};
        const std::vector<double> DoubleDisCo_massReg{DoubleDisCo_massReg_0l, DoubleDisCo_massReg_1l};
        const std::vector<double> DoubleDisCo_disc1{DoubleDisCo_disc1_0l, DoubleDisCo_disc1_1l};
        const std::vector<double> DoubleDisCo_disc2{DoubleDisCo_disc2_0l, DoubleDisCo_disc2_1l};

        const auto& fwm2_top6_0l              = tr.getVar<double>("fwm2_top6_0l");
        const auto& fwm3_top6_0l              = tr.getVar<double>("fwm3_top6_0l");
        const auto& fwm4_top6_0l              = tr.getVar<double>("fwm4_top6_0l");
        const auto& fwm5_top6_0l              = tr.getVar<double>("fwm5_top6_0l");
        const auto& fwm6_top6_0l              = tr.getVar<double>("fwm6_top6_0l");
        const auto& fwm7_top6_0l              = tr.getVar<double>("fwm7_top6_0l");
        const auto& fwm8_top6_0l              = tr.getVar<double>("fwm8_top6_0l");
        const auto& fwm9_top6_0l              = tr.getVar<double>("fwm9_top6_0l");
        const auto& fwm10_top6_0l             = tr.getVar<double>("fwm10_top6_0l");
        const auto& jmt_ev0_top6_0l           = tr.getVar<double>("jmt_ev0_top6_0l");
        const auto& jmt_ev1_top6_0l           = tr.getVar<double>("jmt_ev1_top6_0l");
        const auto& jmt_ev2_top6_0l           = tr.getVar<double>("jmt_ev2_top6_0l");
        const auto& fwm2_top6_1l              = tr.getVar<double>("fwm2_top6_1l");
        const auto& fwm3_top6_1l              = tr.getVar<double>("fwm3_top6_1l");
        const auto& fwm4_top6_1l              = tr.getVar<double>("fwm4_top6_1l");
        const auto& fwm5_top6_1l              = tr.getVar<double>("fwm5_top6_1l");
        const auto& fwm6_top6_1l              = tr.getVar<double>("fwm6_top6_1l");
        const auto& fwm7_top6_1l              = tr.getVar<double>("fwm7_top6_1l");
        const auto& fwm8_top6_1l              = tr.getVar<double>("fwm8_top6_1l");
        const auto& fwm9_top6_1l              = tr.getVar<double>("fwm9_top6_1l");
        const auto& fwm10_top6_1l             = tr.getVar<double>("fwm10_top6_1l");
        const auto& jmt_ev0_top6_1l           = tr.getVar<double>("jmt_ev0_top6_1l");
        const auto& jmt_ev1_top6_1l           = tr.getVar<double>("jmt_ev1_top6_1l");
        const auto& jmt_ev2_top6_1l           = tr.getVar<double>("jmt_ev2_top6_1l");

        const std::vector<double> fwm2_top6{fwm2_top6_0l,    fwm2_top6_1l};
        const std::vector<double> fwm3_top6{fwm3_top6_0l,    fwm3_top6_1l};
        const std::vector<double> fwm4_top6{fwm4_top6_0l,    fwm4_top6_1l};
        const std::vector<double> fwm5_top6{fwm5_top6_0l,    fwm5_top6_1l};
        const std::vector<double> fwm6_top6{fwm6_top6_0l,    fwm6_top6_1l};
        const std::vector<double> fwm7_top6{fwm7_top6_0l,    fwm7_top6_1l};
        const std::vector<double> fwm8_top6{fwm8_top6_0l,    fwm8_top6_1l};
        const std::vector<double> fwm9_top6{fwm9_top6_0l,    fwm9_top6_1l};
        const std::vector<double> fwm10_top6{fwm10_top6_0l,   fwm10_top6_1l};
        const std::vector<double> jmt_ev0_top6{jmt_ev0_top6_0l, jmt_ev0_top6_1l};
        const std::vector<double> jmt_ev1_top6{jmt_ev1_top6_0l, jmt_ev1_top6_1l};
        const std::vector<double> jmt_ev2_top6{jmt_ev2_top6_0l, jmt_ev2_top6_1l};

        const auto& Jets_cm_top6_0l           = tr.getVec<TLorentzVector>("Jets_cm_top6_0l");
        const auto& Jets_cm_top6_1l           = tr.getVec<TLorentzVector>("Jets_cm_top6_1l");

        std::vector<std::vector<TLorentzVector>> Jets_cm_top6{Jets_cm_top6_0l, Jets_cm_top6_1l};

        const auto& eventCounter           = tr.getVar<int>("eventCounter");
        const auto& nMVAJets               = tr.getVar<unsigned int>("nMVAJets");

        const auto& Stop1_pt_cm_OldSeed    = tr.getVar<double>("Stop1_pt_cm_OldSeed");
        const auto& Stop1_eta_cm_OldSeed   = tr.getVar<double>("Stop1_eta_cm_OldSeed");
        const auto& Stop1_phi_cm_OldSeed   = tr.getVar<double>("Stop1_phi_cm_OldSeed");
        const auto& Stop1_mass_cm_OldSeed  = tr.getVar<double>("Stop1_mass_cm_OldSeed");
        const auto& Stop2_pt_cm_OldSeed    = tr.getVar<double>("Stop2_pt_cm_OldSeed");
        const auto& Stop2_eta_cm_OldSeed   = tr.getVar<double>("Stop2_eta_cm_OldSeed");
        const auto& Stop2_phi_cm_OldSeed   = tr.getVar<double>("Stop2_phi_cm_OldSeed");
        const auto& Stop2_mass_cm_OldSeed  = tr.getVar<double>("Stop2_mass_cm_OldSeed");

        // ------------------------
        // -- Print event number
        // ------------------------       
        if(maxevents != -1 && tr.getEvtNum() >= maxevents) break;        
        if(tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() );

        // ------------------------
        // -- Define weight
        // ------------------------
        double weight=1.0, eventweight=1.0, leptonweight=1.0, bTagWeight=1.0, prefiringScaleFactor=1.0, pileupWeight=1.0, htDerivedweight=1.0;
        if(runtype == "MC")
        {
            // Define Lumi weight
            const auto& Weight  = tr.getVar<double>("Weight");
            const auto& lumi = tr.getVar<double>("Lumi");
            eventweight = lumi*Weight;
           
            const auto& eleLepWeight = tr.getVar<double>("totGoodElectronSF");
            const auto& muLepWeight  = tr.getVar<double>("totGoodMuonSF");
            leptonweight = eleLepWeight*muLepWeight;
          
            pileupWeight = tr.getVar<double>("puWeightCorr");
            bTagWeight   = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            htDerivedweight = tr.getVar<double>("htDerivedweight");
            prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor");
            
            weight *= eventweight*leptonweight*bTagWeight*prefiringScaleFactor*pileupWeight*htDerivedweight;
        }

        bool pass_general    = passTriggerMC && passTrigger && passMadHT && passBlind && passMETFilters && passHEMVeto && correct2018Split;

        const std::map<std::string, bool> cut_map 
        {
            {"_1l"                                   , pass_general && passBaseline1l_Good},                         
            {"_0l"                                   , pass_general && passBaseline0l_Good},                         
        };
    
        // Initialize Histograms
        if(!initHistos)
        {
            Preinit(nMVAJets);
            InitHistos(cut_map);
            initHistos = true;
        }

        my_histos["EventCounter"]->Fill(eventCounter);

        std::map<std::string, bool> njetsMap;
        std::map<std::string, bool> ABCDmap;
        for(auto& kv : cut_map)
        {

            int channel = std::stoi(kv.first.substr(1,1));

            njetsMap = {{"Incl",   true},
                          {"6",      NGoodJets[channel]==6},
                          {"7",      NGoodJets[channel]==7},
                          {"8",      NGoodJets[channel]==8},
                          {"9",      NGoodJets[channel]==9},
                          {"10",     NGoodJets[channel]==10},
                          {"11",     NGoodJets[channel]==11},
                          {"11incl", NGoodJets[channel]>=11},
                          {"12",     NGoodJets[channel]==12},
                          {"12incl", NGoodJets[channel]>=12},
                          {"13",     NGoodJets[channel]==13},
                          {"14",     NGoodJets[channel]==14}};

            ABCDmap = {{"",  true},
                       {"A", DoubleDisCo_binA[channel]},
                       {"B", DoubleDisCo_binB[channel]},
                       {"C", DoubleDisCo_binC[channel]},
                       {"D", DoubleDisCo_binD[channel]}};

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
                        double w = weight;

                        my_histos["fwm2_top6"               + name]->Fill(fwm2_top6[channel], w);
                        my_histos["fwm3_top6"               + name]->Fill(fwm3_top6[channel], w);
                        my_histos["fwm4_top6"               + name]->Fill(fwm4_top6[channel], w);
                        my_histos["fwm5_top6"               + name]->Fill(fwm5_top6[channel], w);
                        my_histos["fwm6_top6"               + name]->Fill(fwm6_top6[channel], w);
                        my_histos["fwm7_top6"               + name]->Fill(fwm7_top6[channel], w);
                        my_histos["fwm8_top6"               + name]->Fill(fwm8_top6[channel], w);
                        my_histos["fwm9_top6"               + name]->Fill(fwm9_top6[channel], w);
                        my_histos["fwm10_top6"              + name]->Fill(fwm10_top6[channel], w);
                        my_histos["jmt_ev0_top6"            + name]->Fill(jmt_ev0_top6[channel], w);
                        my_histos["jmt_ev1_top6"            + name]->Fill(jmt_ev1_top6[channel], w);
                        my_histos["jmt_ev2_top6"            + name]->Fill(jmt_ev2_top6[channel], w);

                        if (channel == 1)
                        {
                            my_histos["GoodLeptons_pt_1"        + name]->Fill(GoodLeptons.at(0).second.Pt(), w);
                            my_histos["GoodLeptons_eta_1"       + name]->Fill(GoodLeptons.at(0).second.Eta(), w);
                            my_histos["GoodLeptons_phi_1"       + name]->Fill(GoodLeptons.at(0).second.Phi(), w);
                            my_histos["GoodLeptons_m_1"         + name]->Fill(GoodLeptons.at(0).second.M(), w);
                        }

                        my_histos["Stop1_pt_cm_OldSeed"     + name]->Fill(Stop1_pt_cm_OldSeed, w);
                        my_histos["Stop1_eta_cm_OldSeed"    + name]->Fill(Stop1_eta_cm_OldSeed, w);
                        my_histos["Stop1_phi_cm_OldSeed"    + name]->Fill(Stop1_phi_cm_OldSeed, w);
                        my_histos["Stop1_mass_cm_OldSeed"   + name]->Fill(Stop1_mass_cm_OldSeed, w);
                        my_histos["Stop2_pt_cm_OldSeed"     + name]->Fill(Stop2_pt_cm_OldSeed, w);
                        my_histos["Stop2_eta_cm_OldSeed"    + name]->Fill(Stop2_eta_cm_OldSeed, w);
                        my_histos["Stop2_phi_cm_OldSeed"    + name]->Fill(Stop2_phi_cm_OldSeed, w);
                        my_histos["Stop2_mass_cm_OldSeed"   + name]->Fill(Stop2_mass_cm_OldSeed, w);

                        unsigned int nJets = Jets_cm_top6[channel].size();
                        unsigned int iVec  = 0;

                        for(unsigned int i = 1; i <= nMVAJets; i++)
                        {

                            double pt           = (iVec < nJets) ? static_cast<double>(Jets_cm_top6[channel].at(iVec).Pt())  : 0.0;
                            double eta          = (iVec < nJets) ? static_cast<double>(Jets_cm_top6[channel].at(iVec).Eta()) : 0.0;
                            double phi          = (iVec < nJets) ? static_cast<double>(Jets_cm_top6[channel].at(iVec).Phi()) : 0.0;
                            double m            = (iVec < nJets) ? static_cast<double>(Jets_cm_top6[channel].at(iVec).M())   : 0.0;
   
                            double axismajor    = (iVec < nJets) ? tr.getVar<double>("Jet_axismajor_"    + std::to_string(i) + "_" + std::to_string(channel) + "l") : 0.0;
                            double axisminor    = (iVec < nJets) ? tr.getVar<double>("Jet_axisminor_"    + std::to_string(i) + "_" + std::to_string(channel) + "l") : 0.0;
                            double ptD          = (iVec < nJets) ? tr.getVar<double>("Jet_ptD_"          + std::to_string(i) + "_" + std::to_string(channel) + "l") : 0.0;
                            double multiplicity = (iVec < nJets) ? tr.getVar<double>("Jet_multiplicity_" + std::to_string(i) + "_" + std::to_string(channel) + "l") : 0.0;

                            my_histos["Jet_cm_pt_"           + std::to_string(i) + name]->Fill(pt, w);
                            my_histos["Jet_cm_eta_"          + std::to_string(i) + name]->Fill(eta, w);
                            my_histos["Jet_cm_phi_"          + std::to_string(i) + name]->Fill(phi, w);
                            my_histos["Jet_cm_m_"            + std::to_string(i) + name]->Fill(m, w);
                            my_histos["Jet_cm_axismajor_"    + std::to_string(i) + name]->Fill(axismajor, w);
                            my_histos["Jet_cm_axisminor_"    + std::to_string(i) + name]->Fill(axisminor, w);
                            my_histos["Jet_cm_ptD_"          + std::to_string(i) + name]->Fill(ptD, w);
                            my_histos["Jet_cm_multiplicity_" + std::to_string(i) + name]->Fill(multiplicity, w);

                            iVec++;
                        }

                        my_histos["h_DoubleDisCo_disc1"       + name]->Fill(DoubleDisCo_disc1[channel], w);
                        my_histos["h_DoubleDisCo_disc2"       + name]->Fill(DoubleDisCo_disc2[channel], w);
                        my_histos["h_DoubleDisCo_massReg"     + name]->Fill(DoubleDisCo_massReg[channel], w);
                        my_histos["h_njets"                   + name]->Fill(NGoodJets[channel], w);
                        my_histos["h_ntops"                   + name]->Fill(ntops, w);
                        my_histos["h_nb"                      + name]->Fill(NGoodBJets_pt30, w);
                        my_histos["h_ht"                      + name]->Fill(HT_trigger[channel], w);
                        my_histos["h_mbl"                     + name]->Fill(Mbl, w);
                        my_2d_histos["h_DoubleDisCo_disc1_disc2" + name]->Fill(DoubleDisCo_disc1[channel], DoubleDisCo_disc2[channel], w);

                        if (channel == 1)
                        {
                            for(const auto& l : GoodLeptons)
                            {
                                my_histos["h_lPt"          + name]->Fill(l.second.Pt(), w);
                                my_histos["h_lEta"         + name]->Fill(l.second.Eta(), w);
                                my_histos["h_lPhi"         + name]->Fill(l.second.Phi(), w);
                                my_2d_histos["h_lEta_lPhi" + name]->Fill(l.second.Eta(), l.second.Phi(), w);
                            }
                        }

                        for(unsigned int j = 0; j < Jets.size(); j++)
                        {
                            if(!GoodJets[channel][j]) continue;

                            my_histos["h_jPt"          + name]->Fill(Jets.at(j).Pt(), w);
                            my_histos["h_jEta"         + name]->Fill(Jets.at(j).Eta(), w);
                            my_histos["h_jPhi"         + name]->Fill(Jets.at(j).Phi(), w);
                            my_2d_histos["h_jEta_jPhi" + name]->Fill(Jets.at(j).Eta(), Jets.at(j).Phi(), w);

                        }

                        if( NGoodJets[channel] <= 8 )
                        {
                            my_histos["blind_njets" + name]->Fill(NGoodJets[channel], w);
                            my_histos["blind_ntops" + name]->Fill(ntops, w);
                            my_histos["blind_nb"    + name]->Fill(NGoodBJets[channel], w);
                            my_histos["blind_ht"    + name]->Fill(HT_trigger[channel], w);
                            my_histos["blind_mbl"   + name]->Fill(Mbl, w);

                            if (channel == 1)
                            {
                                for(const auto l : GoodLeptons)
                                {
                                    my_histos["blind_lPt"          + name]->Fill(l.second.Pt(), w);
                                    my_histos["blind_lEta"         + name]->Fill(l.second.Eta(), w);
                                    my_histos["blind_lPhi"         + name]->Fill(l.second.Phi(), w);
                                    my_2d_histos["blind_lEta_lPhi" + name]->Fill(l.second.Eta(), l.second.Phi(), w);
                                }
                            }

                            for(unsigned int j = 0; j < Jets.size(); j++)
                            {
                                if(!GoodJets[channel][j]) continue;

                                my_histos["blind_jPt"          + name]->Fill(Jets.at(j).Pt(), w);
                                my_histos["blind_jEta"         + name]->Fill(Jets.at(j).Eta(), w);
                                my_histos["blind_jPhi"         + name]->Fill(Jets.at(j).Phi(), w);
                                my_2d_histos["blind_jEta_jPhi" + name]->Fill(Jets.at(j).Eta(), Jets.at(j).Phi(), w);

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
    {
        p.second->Write();
    }
    
    for(const auto& p : my_2d_histos) 
    {
        p.second->Write();
    }
}
