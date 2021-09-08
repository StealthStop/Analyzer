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
        {"fwm2_top6",             50,  0, 1},
        {"fwm3_top6",             50,  0, 1},
        {"fwm4_top6",             50,  0, 1},
        {"fwm5_top6",             50,  0, 1},
        {"jmt_ev0_top6",          50,  0, 1},
        {"jmt_ev1_top6",          50,  0, 1},
        {"jmt_ev2_top6",          50,  0, 1},
        {"Stop1_pt_cm_OldSeed",      72, 0, 1500},
        {"Stop1_eta_cm_OldSeed",     80, -6, 6},
        {"Stop1_phi_cm_OldSeed",     64,  -4, 4},
        {"Stop1_mass_cm_OldSeed",       72,  0, 1500},
        {"Stop2_pt_cm_OldSeed",      72, 0, 1500},
        {"Stop2_eta_cm_OldSeed",     80, -6, 6},
        {"Stop2_phi_cm_OldSeed",     64,  -4, 4},
        {"Stop2_mass_cm_OldSeed",       72,  0, 1500},
        {"h_MET_phi", 80, -4, 4},
        {"h_MET_pt" , 150, 0, 150},
        {    "h_njets",           20,   0.0,   20.0},
        {    "h_ntops",           10,   0.0,   10.0},
        {    "h_ht",              500,   0.0, 5000.0},
    };

    hist2DInfos = {
        {"h_DoubleDisCo_disc1_disc2",  100,    0,    1, 100,     0,     1}, 
    };

    abcds = {"", "A", "B", "C", "D"};
    njets = {"Incl", "6", "7", "8", "9", "10", "11", "11incl", "12", "12incl"};

}

void AnalyzeDoubleDisCo::Preinit(unsigned int nNNJets)
{
    for(unsigned int i = 1; i <= nNNJets ; i++)
    {
        histInfos.push_back({"Jet_cm_pt_"           + std::to_string(i), 150, 0, 1500});
        histInfos.push_back({"Jet_cm_eta_"          + std::to_string(i), 100, -6, 6});
        histInfos.push_back({"Jet_cm_phi_"          + std::to_string(i), 80, -4, 4});
        histInfos.push_back({"Jet_cm_m_"            + std::to_string(i), 20, 0, 200});
        histInfos.push_back({"Jet_cm_flavb_"        + std::to_string(i), 80, 0, 1});
        histInfos.push_back({"Jet_cm_flavc_"        + std::to_string(i), 80, 0, 1});
        histInfos.push_back({"Jet_cm_flavg_"        + std::to_string(i), 80, 0, 1});
        histInfos.push_back({"Jet_cm_flavq_"        + std::to_string(i), 80, 0, 1});
        histInfos.push_back({"Jet_cm_flavuds_"      + std::to_string(i), 80, 0, 1});

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

        const auto& GoodJets_pt30          = tr.getVec<bool>("GoodJets_pt30");
        const auto& GoodJets_pt45          = tr.getVec<bool>("GoodJets_pt45");
        const auto& NGoodJets_pt30         = tr.getVar<int>("NGoodJets_pt30");
        const auto& NGoodJets_pt45         = tr.getVar<int>("NGoodJets_pt45");
        const auto& NGoodBJets_pt30        = tr.getVar<int>("NGoodBJets_pt30");
        const auto& NGoodBJets_pt45        = tr.getVar<int>("NGoodBJets_pt45");

        std::vector<std::vector<bool>> GoodJets{GoodJets_pt45, GoodJets_pt30};

        const std::vector<int> NGoodJets{NGoodJets_pt45, NGoodJets_pt30};
        const std::vector<int> NGoodBJets{NGoodBJets_pt45, NGoodBJets_pt30};

        const auto& HT_trigger_pt30        = tr.getVar<double>("HT_trigger_pt30");
        const auto& HT_trigger_pt45        = tr.getVar<double>("HT_trigger_pt45");

        const std::vector<double> HT_trigger{HT_trigger_pt45, HT_trigger_pt30};

        const auto& correct2018Split       = tr.getVar<bool>("correct2018Split");
        const auto& passTrigger            = tr.getVar<bool>("passTrigger");
        const auto& passTriggerMC          = tr.getVar<bool>("passTriggerMC");
        const auto& passMETFilters         = tr.getVar<bool>("passMETFilters");
        const auto& passMadHT              = tr.getVar<bool>("passMadHT");
        const auto& passHEMVeto            = tr.getVar<bool>("passHEMVeto");

        const auto& passBlind              = tr.getVar<bool>("passBlindLep_Good");            
        const auto& passBaseline0l_Good    = tr.getVar<bool>("passBaseline0l_Good_Loose");
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
        const auto& DoubleDisCo_disc1_NonIsoMuon_0l   = tr.getVar<double>("DoubleDisCo_disc1_NonIsoMuon_0l");
        const auto& DoubleDisCo_disc2_NonIsoMuon_0l   = tr.getVar<double>("DoubleDisCo_disc2_NonIsoMuon_0l");
        const auto& DoubleDisCo_disc1_NonIsoMuon_1l   = tr.getVar<double>("DoubleDisCo_disc1_NonIsoMuon_1l");
        const auto& DoubleDisCo_disc2_NonIsoMuon_1l   = tr.getVar<double>("DoubleDisCo_disc2_NonIsoMuon_1l");

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
        const auto& jmt_ev0_top6_0l           = tr.getVar<double>("jmt_ev0_top6_0l");
        const auto& jmt_ev1_top6_0l           = tr.getVar<double>("jmt_ev1_top6_0l");
        const auto& jmt_ev2_top6_0l           = tr.getVar<double>("jmt_ev2_top6_0l");
        const auto& fwm2_top6_1l              = tr.getVar<double>("fwm2_top6_1l");
        const auto& fwm3_top6_1l              = tr.getVar<double>("fwm3_top6_1l");
        const auto& fwm4_top6_1l              = tr.getVar<double>("fwm4_top6_1l");
        const auto& fwm5_top6_1l              = tr.getVar<double>("fwm5_top6_1l");
        const auto& jmt_ev0_top6_1l           = tr.getVar<double>("jmt_ev0_top6_1l");
        const auto& jmt_ev1_top6_1l           = tr.getVar<double>("jmt_ev1_top6_1l");
        const auto& jmt_ev2_top6_1l           = tr.getVar<double>("jmt_ev2_top6_1l");

        const std::vector<double> fwm2_top6{fwm2_top6_0l,    fwm2_top6_1l};
        const std::vector<double> fwm3_top6{fwm3_top6_0l,    fwm3_top6_1l};
        const std::vector<double> fwm4_top6{fwm4_top6_0l,    fwm4_top6_1l};
        const std::vector<double> fwm5_top6{fwm5_top6_0l,    fwm5_top6_1l};
        const std::vector<double> jmt_ev0_top6{jmt_ev0_top6_0l, jmt_ev0_top6_1l};
        const std::vector<double> jmt_ev1_top6{jmt_ev1_top6_0l, jmt_ev1_top6_1l};
        const std::vector<double> jmt_ev2_top6{jmt_ev2_top6_0l, jmt_ev2_top6_1l};

        const auto& Jets_cm_top6_0l           = tr.getVec<TLorentzVector>("Jets_cm_top6_0l");
        const auto& Jets_cm_top6_1l           = tr.getVec<TLorentzVector>("Jets_cm_top6_1l");

        std::vector<std::vector<TLorentzVector>> Jets_cm_top6{Jets_cm_top6_0l, Jets_cm_top6_1l};

        const auto& eventCounter           = tr.getVar<int>("eventCounter");
        const auto& nMVAJets_0l            = tr.getVar<unsigned int>("nMVAJets_0l");
        const auto& nMVAJets_1l            = tr.getVar<unsigned int>("nMVAJets_1l");

        std::vector<unsigned int> nMVAJets{nMVAJets_0l, nMVAJets_1l};

        const auto& Stop1_pt_cm_OldSeed    = tr.getVar<double>("Stop1_pt_cm_OldSeed");
        const auto& Stop1_eta_cm_OldSeed   = tr.getVar<double>("Stop1_eta_cm_OldSeed");
        const auto& Stop1_phi_cm_OldSeed   = tr.getVar<double>("Stop1_phi_cm_OldSeed");
        const auto& Stop1_mass_cm_OldSeed  = tr.getVar<double>("Stop1_mass_cm_OldSeed");
        const auto& Stop2_pt_cm_OldSeed    = tr.getVar<double>("Stop2_pt_cm_OldSeed");
        const auto& Stop2_eta_cm_OldSeed   = tr.getVar<double>("Stop2_eta_cm_OldSeed");
        const auto& Stop2_phi_cm_OldSeed   = tr.getVar<double>("Stop2_phi_cm_OldSeed");
        const auto& Stop2_mass_cm_OldSeed  = tr.getVar<double>("Stop2_mass_cm_OldSeed");

        const auto& met                    = tr.getVar<double>("MET");
        const auto& metPhi                 = tr.getVar<double>("METPhi");

        const auto& blindABCDdata          = tr.getVar<bool>("blind") && tr.getVar<std::string>("runtype")!="MC";

        std::vector<double> Jets_flavb_0l;   std::vector<double> Jets_flavb_1l;
        std::vector<double> Jets_flavc_0l;   std::vector<double> Jets_flavc_1l;
        std::vector<double> Jets_flavg_0l;   std::vector<double> Jets_flavg_1l;
        std::vector<double> Jets_flavuds_0l; std::vector<double> Jets_flavuds_1l;
        std::vector<double> Jets_flavq_0l;   std::vector<double> Jets_flavq_1l;
        for (unsigned int iJet = 1; iJet < nMVAJets_0l+1; iJet++) {
            Jets_flavb_0l.push_back(tr.getVar<double>("Jet_flavb_"+std::to_string(iJet)+"_0l"));
            Jets_flavc_0l.push_back(tr.getVar<double>("Jet_flavc_"+std::to_string(iJet)+"_0l"));
            Jets_flavg_0l.push_back(tr.getVar<double>("Jet_flavg_"+std::to_string(iJet)+"_0l"));
            Jets_flavuds_0l.push_back(tr.getVar<double>("Jet_flavuds_"+std::to_string(iJet)+"_0l"));
            Jets_flavq_0l.push_back(tr.getVar<double>("Jet_flavq_"+std::to_string(iJet)+"_0l"));
        }

        for (unsigned int iJet = 1; iJet < nMVAJets_1l+1; iJet++) {
            Jets_flavb_1l.push_back(tr.getVar<double>("Jet_flavb_"+std::to_string(iJet)+"_1l"));
            Jets_flavc_1l.push_back(tr.getVar<double>("Jet_flavc_"+std::to_string(iJet)+"_1l"));
            Jets_flavg_1l.push_back(tr.getVar<double>("Jet_flavg_"+std::to_string(iJet)+"_1l"));
            Jets_flavuds_1l.push_back(tr.getVar<double>("Jet_flavuds_"+std::to_string(iJet)+"_1l"));
            Jets_flavq_1l.push_back(tr.getVar<double>("Jet_flavq_"+std::to_string(iJet)+"_1l"));
        }

        std::vector<std::vector<double>> Jets_flavb{Jets_flavb_0l, Jets_flavb_1l};
        std::vector<std::vector<double>> Jets_flavc{Jets_flavc_0l, Jets_flavc_1l};
        std::vector<std::vector<double>> Jets_flavg{Jets_flavg_0l, Jets_flavg_1l};
        std::vector<std::vector<double>> Jets_flavuds{Jets_flavuds_0l, Jets_flavuds_1l};
        std::vector<std::vector<double>> Jets_flavq{Jets_flavq_0l, Jets_flavq_1l};

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
            Preinit(nMVAJets_0l > nMVAJets_1l ? nMVAJets_0l : nMVAJets_1l);
            InitHistos(cut_map);
            initHistos = true;
        }

        my_histos["EventCounter"]->Fill(eventCounter);

        std::map<std::string, bool> njetsMap;
        std::map<std::string, bool> ABCDmap;
        for(auto& kv : cut_map)
        {

            int channel = std::stoi(kv.first.substr(1,1));

            njetsMap = {{"Incl",     true},
                          {"6",      NGoodJets[channel]==6},
                          {"7",      NGoodJets[channel]==7},
                          {"8",      NGoodJets[channel]==8},
                          {"9",      NGoodJets[channel]==9},
                          {"10",     NGoodJets[channel]==10},
                          {"11",     NGoodJets[channel]==11},
                          {"11incl", NGoodJets[channel]>=11},
                          {"12",     NGoodJets[channel]==12},
                          {"12incl", NGoodJets[channel]>=12}
            };

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

                    if (blindABCDdata && jpass.first!="")
                        continue;
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
                        my_histos["jmt_ev0_top6"            + name]->Fill(jmt_ev0_top6[channel], w);
                        my_histos["jmt_ev1_top6"            + name]->Fill(jmt_ev1_top6[channel], w);
                        my_histos["jmt_ev2_top6"            + name]->Fill(jmt_ev2_top6[channel], w);

                        my_histos["Stop1_pt_cm_OldSeed"     + name]->Fill(Stop1_pt_cm_OldSeed, w);
                        my_histos["Stop1_eta_cm_OldSeed"    + name]->Fill(Stop1_eta_cm_OldSeed, w);
                        my_histos["Stop1_phi_cm_OldSeed"    + name]->Fill(Stop1_phi_cm_OldSeed, w);
                        my_histos["Stop1_mass_cm_OldSeed"   + name]->Fill(Stop1_mass_cm_OldSeed, w);
                        my_histos["Stop2_pt_cm_OldSeed"     + name]->Fill(Stop2_pt_cm_OldSeed, w);
                        my_histos["Stop2_eta_cm_OldSeed"    + name]->Fill(Stop2_eta_cm_OldSeed, w);
                        my_histos["Stop2_phi_cm_OldSeed"    + name]->Fill(Stop2_phi_cm_OldSeed, w);
                        my_histos["Stop2_mass_cm_OldSeed"   + name]->Fill(Stop2_mass_cm_OldSeed, w);

                        unsigned int nJets = Jets_cm_top6[channel].size();
                        for(unsigned int i = 1; i <= nMVAJets[channel]; i++)
                        {
                            double pt           = (i-1 < nJets) ? static_cast<double>(Jets_cm_top6[channel].at(i-1).Pt())  : 0.0;
                            double eta          = (i-1 < nJets) ? static_cast<double>(Jets_cm_top6[channel].at(i-1).Eta()) : 0.0;
                            double phi          = (i-1 < nJets) ? static_cast<double>(Jets_cm_top6[channel].at(i-1).Phi()) : 0.0;
                            double m            = (i-1 < nJets) ? static_cast<double>(Jets_cm_top6[channel].at(i-1).M())   : 0.0;
   
                            my_histos["Jet_cm_pt_"           + std::to_string(i) + name]->Fill(pt, w);
                            my_histos["Jet_cm_eta_"          + std::to_string(i) + name]->Fill(eta, w);
                            my_histos["Jet_cm_phi_"          + std::to_string(i) + name]->Fill(phi, w);
                            my_histos["Jet_cm_m_"            + std::to_string(i) + name]->Fill(m, w);

                            my_histos["Jet_cm_flavb_"   + std::to_string(i) + name]->Fill(Jets_flavb[channel].at(i-1), w);
                            my_histos["Jet_cm_flavc_"   + std::to_string(i) + name]->Fill(Jets_flavc[channel].at(i-1), w);
                            my_histos["Jet_cm_flavg_"   + std::to_string(i) + name]->Fill(Jets_flavg[channel].at(i-1), w);
                            my_histos["Jet_cm_flavq_"   + std::to_string(i) + name]->Fill(Jets_flavq[channel].at(i-1), w);
                            my_histos["Jet_cm_flavuds_" + std::to_string(i) + name]->Fill(Jets_flavuds[channel].at(i-1), w);
                        }

                        if (! blindABCDdata) {
                            my_histos["h_DoubleDisCo_disc1"       + name]->Fill(DoubleDisCo_disc1[channel], w);
                            my_histos["h_DoubleDisCo_disc2"       + name]->Fill(DoubleDisCo_disc2[channel], w);
                            my_histos["h_DoubleDisCo_massReg"     + name]->Fill(DoubleDisCo_massReg[channel], w);
                            my_2d_histos["h_DoubleDisCo_disc1_disc2" + name]->Fill(DoubleDisCo_disc1[channel], DoubleDisCo_disc2[channel], w);
                        }
                        my_histos["h_njets"                   + name]->Fill(NGoodJets[channel], w);
                        my_histos["h_ntops"                   + name]->Fill(ntops, w);
                        my_histos["h_ht"                      + name]->Fill(HT_trigger[channel], w);
                        my_histos["h_MET_phi"                 + name]->Fill(metPhi, w);
                        my_histos["h_MET_pt"                  + name]->Fill(met, w);

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
