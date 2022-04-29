#define HEM_Analyzer_cxx
#include "Analyzer/Analyzer/include/HEM_Analyzer.h"
#include "Framework/Framework/include/Utility.h"
#include "NTupleReader/include/NTupleReader.h"

#include <iostream>

HEM_Analyzer::HEM_Analyzer() : inithisto(false) // define inithisto variable
{
}

// -------------------
// -- Define histos
// -------------------
void HEM_Analyzer::InitHistos(const std::map<std::string, bool>& cutmap) // define variable map
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    
    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ;
    
    for (const auto& cutVar : cutmap) 
    {
        // ----------------
        // before HEM issue
        // ----------------
        // 1D histograms
        my_histos.emplace( "h_njets_"        +cutVar.first, std::make_shared<TH1D> ( ("h_njets_"        +cutVar.first).c_str(), ("h_njets_"        +cutVar.first).c_str(), 13, 6.5, 19.5     ) );
        my_histos.emplace( "h_nbjets_"       +cutVar.first, std::make_shared<TH1D> ( ("h_nbjets_"       +cutVar.first).c_str(), ("h_nbjets_"       +cutVar.first).c_str(), 15, -0.5, 14.5    ) );
        my_histos.emplace( "h_met_"          +cutVar.first, std::make_shared<TH1D> ( ("h_met_"          +cutVar.first).c_str(), ("h_met_"          +cutVar.first).c_str(), 720, 0, 1500      ) );
        my_histos.emplace( "h_ht_"           +cutVar.first, std::make_shared<TH1D> ( ("h_ht_"           +cutVar.first).c_str(), ("h_ht_"           +cutVar.first).c_str(), 720, 300, 5000    ) );
        my_histos.emplace( "h_jetPt_"        +cutVar.first, std::make_shared<TH1D> ( ("h_jetPt_"        +cutVar.first).c_str(), ("h_jetPt_"        +cutVar.first).c_str(), 1440, 0, 5000     ) );
        my_histos.emplace( "h_jetPtMax_"     +cutVar.first, std::make_shared<TH1D> ( ("h_jetPtMax_"     +cutVar.first).c_str(), ("h_jetPtMax_"     +cutVar.first).c_str(), 1440, 0, 5000     ) );
        my_histos.emplace( "h_lvMET_cm_mass_"+cutVar.first, std::make_shared<TH1D> ( ("h_lvMET_cm_mass_"+cutVar.first).c_str(), ("h_lvMET_cm_mass_"+cutVar.first).c_str(), 1440, 0, 1000     ) );
        my_histos.emplace( "h_lvMET_cm_eta_" +cutVar.first, std::make_shared<TH1D> ( ("h_lvMET_cm_eta_" +cutVar.first).c_str(), ("h_lvMET_cm_eta_" +cutVar.first).c_str(), 1440, -2.5, 2.5   ) );
        my_histos.emplace( "h_lvMET_cm_phi_" +cutVar.first, std::make_shared<TH1D> ( ("h_lvMET_cm_phi_" +cutVar.first).c_str(), ("h_lvMET_cm_phi_" +cutVar.first).c_str(), 1440, -3.14, 3.14 ) );
        my_histos.emplace( "h_lvMET_cm_pt_"  +cutVar.first, std::make_shared<TH1D> ( ("h_lvMET_cm_pt_"  +cutVar.first).c_str(), ("h_lvMET_cm_pt_"  +cutVar.first).c_str(), 1440, 0, 1000     ) );
        my_histos.emplace( "h_fwm2_top6_"    +cutVar.first, std::make_shared<TH1D> ( ("h_fwm2_top6_"    +cutVar.first).c_str(), ("h_fwm2_top6_"    +cutVar.first).c_str(), 1440, 0, 1        ) );
        my_histos.emplace( "h_fwm3_top6_"    +cutVar.first, std::make_shared<TH1D> ( ("h_fwm3_top6_"    +cutVar.first).c_str(), ("h_fwm3_top6_"    +cutVar.first).c_str(), 1440, 0, 1        ) );
        my_histos.emplace( "h_fwm4_top6_"    +cutVar.first, std::make_shared<TH1D> ( ("h_fwm4_top6_"    +cutVar.first).c_str(), ("h_fwm4_top6_"    +cutVar.first).c_str(), 1440, 0, 1        ) );
        my_histos.emplace( "h_fwm5_top6_"    +cutVar.first, std::make_shared<TH1D> ( ("h_fwm5_top6_"    +cutVar.first).c_str(), ("h_fwm5_top6_"    +cutVar.first).c_str(), 1440, 0, 1        ) );
        my_histos.emplace( "h_fwm6_top6_"    +cutVar.first, std::make_shared<TH1D> ( ("h_fwm6_top6_"    +cutVar.first).c_str(), ("h_fwm6_top6_"    +cutVar.first).c_str(), 1440, 0, 1        ) );
        my_histos.emplace( "h_fwm7_top6_"    +cutVar.first, std::make_shared<TH1D> ( ("h_fwm7_top6_"    +cutVar.first).c_str(), ("h_fwm7_top6_"    +cutVar.first).c_str(), 1440, 0, 1        ) );
        my_histos.emplace( "h_fwm8_top6_"    +cutVar.first, std::make_shared<TH1D> ( ("h_fwm8_top6_"    +cutVar.first).c_str(), ("h_fwm8_top6_"    +cutVar.first).c_str(), 1440, 0, 1        ) );
        my_histos.emplace( "h_fwm9_top6_"    +cutVar.first, std::make_shared<TH1D> ( ("h_fwm9_top6_"    +cutVar.first).c_str(), ("h_fwm9_top6_"    +cutVar.first).c_str(), 1440, 0, 1        ) );
        my_histos.emplace( "h_fwm10_top6_"   +cutVar.first, std::make_shared<TH1D> ( ("h_fwm10_top6_"   +cutVar.first).c_str(), ("h_fwm10_top6_"   +cutVar.first).c_str(), 1440, 0, 1        ) );
        my_histos.emplace( "h_jmt_ev0_top6_" +cutVar.first, std::make_shared<TH1D> ( ("h_jmt_ev0_top6_" +cutVar.first).c_str(), ("h_jmt_ev0_top6_" +cutVar.first).c_str(), 1440, 0, 1        ) );
        my_histos.emplace( "h_jmt_ev1_top6_" +cutVar.first, std::make_shared<TH1D> ( ("h_jmt_ev1_top6_" +cutVar.first).c_str(), ("h_jmt_ev1_top6_" +cutVar.first).c_str(), 1440, 0, 1        ) );
        my_histos.emplace( "h_jmt_ev2_top6_" +cutVar.first, std::make_shared<TH1D> ( ("h_jmt_ev2_top6_" +cutVar.first).c_str(), ("h_jmt_ev2_top6_" +cutVar.first).c_str(), 1440, 0, 1        ) );
        my_histos.emplace( "h_event_beta_z_" +cutVar.first, std::make_shared<TH1D> ( ("h_event_beta_z_" +cutVar.first).c_str(), ("h_event_beta_z_" +cutVar.first).c_str(), 1440, -1, 1       ) );
        // 2D histograms
        my_2d_histos.emplace( "h_jet_EtaVsPhi_"     +cutVar.first, std::make_shared<TH2D>( ("h_jet_EtaVsPhi_"     +cutVar.first).c_str(), ("h_jet_EtaVsPhi_"     +cutVar.first).c_str(), 1440, -2.5, 2.5, 1440, -3.14, 3.14 ) );
        my_2d_histos.emplace( "h_bjet_EtaVsPhi_"    +cutVar.first, std::make_shared<TH2D>( ("h_bjet_EtaVsPhi_"    +cutVar.first).c_str(), ("h_bjet_EtaVsPhi_"    +cutVar.first).c_str(), 1440, -2.5, 2.5, 1440, -3.14, 3.14 ) );
        my_2d_histos.emplace( "h_electron_EtaVsPhi_"+cutVar.first, std::make_shared<TH2D>( ("h_electron_EtaVsPhi_"+cutVar.first).c_str(), ("h_electron_EtaVsPhi_"+cutVar.first).c_str(), 1440, -2.5, 2.5, 1440, -3.14, 3.14 ) );
        my_2d_histos.emplace( "h_muon_EtaVsPhi_"    +cutVar.first, std::make_shared<TH2D>( ("h_muon_EtaVsPhi_"    +cutVar.first).c_str(), ("h_muon_EtaVsPhi_"    +cutVar.first).c_str(), 1440, -2.5, 2.5, 1440, -3.14, 3.14 ) );

        // ---------------
        // after HEM issue
        // ---------------
        // 1D histograms
        my_histos.emplace( "h_njets_HEM_"        +cutVar.first, std::make_shared<TH1D> ( ("h_njets_HEM_"        +cutVar.first).c_str(), ("h_njets_HEM_"        +cutVar.first).c_str(), 13, 6.5, 19.5     ) );
        my_histos.emplace( "h_nbjets_HEM_"       +cutVar.first, std::make_shared<TH1D> ( ("h_nbjets_HEM_"       +cutVar.first).c_str(), ("h_nbjets_HEM_"       +cutVar.first).c_str(), 15, -0.5, 14.5    ) );
        my_histos.emplace( "h_met_HEM_"          +cutVar.first, std::make_shared<TH1D> ( ("h_met_HEM_"          +cutVar.first).c_str(), ("h_met_HEM_"          +cutVar.first).c_str(), 720, 0, 1500      ) );
        my_histos.emplace( "h_ht_HEM_"           +cutVar.first, std::make_shared<TH1D> ( ("h_ht_HEM_"           +cutVar.first).c_str(), ("h_ht_HEM_"           +cutVar.first).c_str(), 720, 300, 5000    ) );
        my_histos.emplace( "h_jetPt_HEM_"        +cutVar.first, std::make_shared<TH1D> ( ("h_jetPt_HEM_"        +cutVar.first).c_str(), ("h_jetPt_HEM_"        +cutVar.first).c_str(), 1440, 0, 5000     ) );
        my_histos.emplace( "h_jetPtMax_HEM_"     +cutVar.first, std::make_shared<TH1D> ( ("h_jetPtMax_HEM_"     +cutVar.first).c_str(), ("h_jetPtMax_HEM_"     +cutVar.first).c_str(), 1440, 0, 5000     ) );
        my_histos.emplace( "h_lvMET_cm_mass_HEM_"+cutVar.first, std::make_shared<TH1D> ( ("h_lvMET_cm_mass_HEM_"+cutVar.first).c_str(), ("h_lvMET_cm_mass_HEM_"+cutVar.first).c_str(), 1440, 0, 1000     ) );
        my_histos.emplace( "h_lvMET_cm_eta_HEM_" +cutVar.first, std::make_shared<TH1D> ( ("h_lvMET_cm_eta_HEM_" +cutVar.first).c_str(), ("h_lvMET_cm_eta_HEM_" +cutVar.first).c_str(), 1440, -2.5, 2.5   ) );
        my_histos.emplace( "h_lvMET_cm_phi_HEM_" +cutVar.first, std::make_shared<TH1D> ( ("h_lvMET_cm_phi_HEM_" +cutVar.first).c_str(), ("h_lvMET_cm_phi_HEM_" +cutVar.first).c_str(), 1440, -3.14, 3.14 ) );
        my_histos.emplace( "h_lvMET_cm_pt_HEM_"  +cutVar.first, std::make_shared<TH1D> ( ("h_lvMET_cm_pt_HEM_"  +cutVar.first).c_str(), ("h_lvMET_cm_pt_HEM_"  +cutVar.first).c_str(), 1440, 0, 1000     ) );
        my_histos.emplace( "h_fwm2_top6_HEM_"    +cutVar.first, std::make_shared<TH1D> ( ("h_fwm2_top6_HEM_"    +cutVar.first).c_str(), ("h_fwm2_top6_HEM_"    +cutVar.first).c_str(), 1440, 0, 1        ) );
        my_histos.emplace( "h_fwm3_top6_HEM_"    +cutVar.first, std::make_shared<TH1D> ( ("h_fwm3_top6_HEM_"    +cutVar.first).c_str(), ("h_fwm3_top6_HEM_"    +cutVar.first).c_str(), 1440, 0, 1        ) );
        my_histos.emplace( "h_fwm4_top6_HEM_"    +cutVar.first, std::make_shared<TH1D> ( ("h_fwm4_top6_HEM_"    +cutVar.first).c_str(), ("h_fwm4_top6_HEM_"    +cutVar.first).c_str(), 1440, 0, 1        ) );
        my_histos.emplace( "h_fwm5_top6_HEM_"    +cutVar.first, std::make_shared<TH1D> ( ("h_fwm5_top6_HEM_"    +cutVar.first).c_str(), ("h_fwm5_top6_HEM_"    +cutVar.first).c_str(), 1440, 0, 1        ) );
        my_histos.emplace( "h_fwm6_top6_HEM_"    +cutVar.first, std::make_shared<TH1D> ( ("h_fwm6_top6_HEM_"    +cutVar.first).c_str(), ("h_fwm6_top6_HEM_"    +cutVar.first).c_str(), 1440, 0, 1        ) );
        my_histos.emplace( "h_fwm7_top6_HEM_"    +cutVar.first, std::make_shared<TH1D> ( ("h_fwm7_top6_HEM_"    +cutVar.first).c_str(), ("h_fwm7_top6_HEM_"    +cutVar.first).c_str(), 1440, 0, 1        ) );
        my_histos.emplace( "h_fwm8_top6_HEM_"    +cutVar.first, std::make_shared<TH1D> ( ("h_fwm8_top6_HEM_"    +cutVar.first).c_str(), ("h_fwm8_top6_HEM_"    +cutVar.first).c_str(), 1440, 0, 1        ) );
        my_histos.emplace( "h_fwm9_top6_HEM_"    +cutVar.first, std::make_shared<TH1D> ( ("h_fwm9_top6_HEM_"    +cutVar.first).c_str(), ("h_fwm9_top6_HEM_"    +cutVar.first).c_str(), 1440, 0, 1        ) );
        my_histos.emplace( "h_fwm10_top6_HEM_"   +cutVar.first, std::make_shared<TH1D> ( ("h_fwm10_top6_HEM_"   +cutVar.first).c_str(), ("h_fwm10_top6_HEM_"   +cutVar.first).c_str(), 1440, 0, 1        ) );
        my_histos.emplace( "h_jmt_ev0_top6_HEM_" +cutVar.first, std::make_shared<TH1D> ( ("h_jmt_ev0_top6_HEM_" +cutVar.first).c_str(), ("h_jmt_ev0_top6_HEM_" +cutVar.first).c_str(), 1440, 0, 1        ) );
        my_histos.emplace( "h_jmt_ev1_top6_HEM_" +cutVar.first, std::make_shared<TH1D> ( ("h_jmt_ev1_top6_HEM_" +cutVar.first).c_str(), ("h_jmt_ev1_top6_HEM_" +cutVar.first).c_str(), 1440, 0, 1        ) );
        my_histos.emplace( "h_jmt_ev2_top6_HEM_" +cutVar.first, std::make_shared<TH1D> ( ("h_jmt_ev2_top6_HEM_" +cutVar.first).c_str(), ("h_jmt_ev2_top6_HEM_" +cutVar.first).c_str(), 1440, 0, 1        ) );
        my_histos.emplace( "h_event_beta_z_HEM_" +cutVar.first, std::make_shared<TH1D> ( ("h_event_beta_z_HEM_" +cutVar.first).c_str(), ("h_event_beta_z_HEM_" +cutVar.first).c_str(), 1440, -1, 1       ) ); 
        // 2D histograms
        my_2d_histos.emplace( "h_jet_EtaVsPhi_HEM_"     +cutVar.first, std::make_shared<TH2D>( ("h_jet_EtaVsPhi_HEM_"     +cutVar.first).c_str(), ("h_jet_EtaVsPhi_HEM_"     +cutVar.first).c_str(), 1440, -2.5, 2.5, 1440, -3.14, 3.14 ) );
        my_2d_histos.emplace( "h_bjet_EtaVsPhi_HEM_"    +cutVar.first, std::make_shared<TH2D>( ("h_bjet_EtaVsPhi_HEM_"    +cutVar.first).c_str(), ("h_bjet_EtaVsPhi_HEM_"     +cutVar.first).c_str(), 1440, -2.5, 2.5, 1440, -3.14, 3.14 ) );
        my_2d_histos.emplace( "h_electron_EtaVsPhi_HEM_"+cutVar.first, std::make_shared<TH2D>( ("h_electron_EtaVsPhi_HEM_"+cutVar.first).c_str(), ("h_electron_EtaVsPhi_HEM_"+cutVar.first).c_str(), 1440, -2.5, 2.5, 1440, -3.14, 3.14 ) );
        my_2d_histos.emplace( "h_muon_EtaVsPhi_HEM_"    +cutVar.first, std::make_shared<TH2D>( ("h_muon_EtaVsPhi_HEM_"    +cutVar.first).c_str(), ("h_muon_EtaVsPhi_HEM_"    +cutVar.first).c_str(), 1440, -2.5, 2.5, 1440, -3.14, 3.14 ) );

    }
}

// ---------------------------------------------
// -- Put everything you want to do per event 
// ---------------------------------------------
void HEM_Analyzer::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    while( tr.getNextEvent() )
    {
        
        const auto& eventCounter = tr.getVar<int>("eventCounter");

        // ------------------------
        // -- Print Event Number 
        // ------------------------
        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() % 1000 == 0) printf( " Event %i\n", tr.getEvtNum() );

        const auto& runtype                   = tr.getVar<std::string>("runtype");     
        const auto& RunNum                    = tr.getVar<unsigned int>("RunNum");
        const auto& Jets                      = tr.getVec<utility::LorentzVector>("Jets");
        const auto& GoodLeptons               = tr.getVec<std::pair<std::string, utility::LorentzVector>>("GoodLeptons");
        const auto& GoodJets_pt30             = tr.getVec<bool>("GoodJets_pt30");
        const auto& GoodBJets_pt30            = tr.getVec<bool>("GoodBJets_pt30");
        const auto& NGoodLeptons              = tr.getVar<int>("NGoodLeptons");
        const auto& NGoodJets_pt30            = tr.getVar<int>("NGoodJets_pt30");
        const auto& NGoodBJets_pt30           = tr.getVar<int>("NGoodBJets_pt30");
        const auto& HT_trigger_pt30           = tr.getVar<double>("HT_trigger_pt30");
        const auto& met                       = tr.getVar<float>("MET");
        const auto& lvMET_cm_m                = tr.getVar<double>("lvMET_cm_m");
        const auto& lvMET_cm_eta              = tr.getVar<double>("lvMET_cm_eta");
        const auto& lvMET_cm_phi              = tr.getVar<double>("lvMET_cm_phi");
        const auto& lvMET_cm_pt               = tr.getVar<double>("lvMET_cm_pt");
        const auto& fwm2_top6                 = tr.getVar<double>("fwm2_top6");
        const auto& fwm3_top6                 = tr.getVar<double>("fwm3_top6");
        const auto& fwm4_top6                 = tr.getVar<double>("fwm4_top6");
        const auto& fwm5_top6                 = tr.getVar<double>("fwm5_top6");
        const auto& fwm6_top6                 = tr.getVar<double>("fwm6_top6");
        const auto& fwm7_top6                 = tr.getVar<double>("fwm7_top6");
        const auto& fwm8_top6                 = tr.getVar<double>("fwm8_top6");
        const auto& fwm9_top6                 = tr.getVar<double>("fwm9_top6");
        const auto& fwm10_top6                = tr.getVar<double>("fwm10_top6");
        const auto& jmt_ev0_top6              = tr.getVar<double>("jmt_ev0_top6");
        const auto& jmt_ev1_top6              = tr.getVar<double>("jmt_ev1_top6");
        const auto& jmt_ev2_top6              = tr.getVar<double>("jmt_ev2_top6");
        const auto& event_beta_z              = tr.getVar<double>("event_beta_z");
        const auto& passMadHT                 = tr.getVar<bool>("passMadHT");
        const auto& passBaseline_0l           = tr.getVar<bool>("passBaseline0l_good");
        const auto& passBaseline_1l           = tr.getVar<bool>("passBaseline1l_Good");
        const auto& passBaseline_0l_noHEMveto = tr.getVar<bool>("passBaseline0l_good_noHEMveto");
        const auto& passBaseline_1l_noHEMveto = tr.getVar<bool>("passBaseline1l_Good_noHEMveto");

        // -------------------
        // -- Define weight
        // -------------------
        double weight               = 1.0;
        double eventweight          = 1.0;
        double bTagScaleFactor      = 1.0;
        double puScaleFactor        = 1.0;
        double prefiringScaleFactor = 1.0;

        double leptonScaleFactor    = 1.0;
        double htDerivedScaleFactor = 1.0;
        
        if(runtype == "MC")
        {
            eventweight          = tr.getVar<float>("LumiXsec");
            bTagScaleFactor      = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            puScaleFactor        = tr.getVar<double>("puWeightCorr");
            prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor");

            // Define lepton weight
            if(NGoodLeptons == 1)
            {
                const auto& eleLepWeight = tr.getVar<double>("totGoodElectronSF");
                const auto& muLepWeight  = tr.getVar<double>("totGoodMuonSF");
                leptonScaleFactor = (GoodLeptons[0].first == "e") ? eleLepWeight : muLepWeight;
            }
    
            weight *= eventweight*leptonScaleFactor*bTagScaleFactor*prefiringScaleFactor*puScaleFactor;
        }

        // -------------------------------------------------
        // -- Make cuts and fill histograms here & cutmap
        // -------------------------------------------------
        const std::map<std::string, bool>& cutmap
        {
            {"passBaseline_0l",           passBaseline_0l},
            {"passBaseline_1l",           passBaseline_1l}, 
            {"passBaseline_0l_noHEMveto", passBaseline_0l_noHEMveto},
            {"passBaseline_1l_noHEMveto", passBaseline_1l_noHEMveto}, 
        };

        if (!inithisto) 
        {
            InitHistos(cutmap);
            inithisto = true;
        }

        my_histos["EventCounter"]->Fill( eventCounter );

        // --------------------------------
        // -- Fill the cutmap histograms
        // --------------------------------     
        for (const auto& cutVar: cutmap) 
        {    
            if (cutVar.second) 
            {
                double jetPtMax = 0.0;
                
                // ----------------                
                // before HEM issue
                // ----------------
                if ( RunNum < 319077 ) 
                {
                    for (unsigned int goodJet = 0; goodJet < Jets.size(); goodJet++) 
                    {
                        if (!GoodJets_pt30[goodJet]) { continue; }

                        my_histos["h_jetPt_"          +cutVar.first]->Fill( Jets[goodJet].Pt(), weight                       );
                        my_2d_histos["h_jet_EtaVsPhi_"+cutVar.first]->Fill( Jets[goodJet].Eta(), Jets[goodJet].Phi(), weight ); 

                        if (GoodBJets_pt30[goodJet])
                            my_2d_histos["h_bjet_EtaVsPhi_"+cutVar.first]->Fill( Jets[goodJet].Eta(), Jets[goodJet].Phi(), weight ); 


                        if (Jets[goodJet].Pt() > jetPtMax) 
                        { jetPtMax = Jets[goodJet].Pt(); }
                    }
        
                    if (cutVar.first.find("passBaseline_1l") != std::string::npos) 
                    {    
                        for (unsigned int goodLep = 0; goodLep < GoodLeptons.size(); goodLep++) 
                        {
                            if (GoodLeptons[goodLep].first == "e") 
                            {
                                my_2d_histos["h_electron_EtaVsPhi_"+cutVar.first]->Fill(GoodLeptons[goodLep].second.Eta(), GoodLeptons[goodLep].second.Phi(), weight);
                            } 
                            else 
                            {
                                my_2d_histos["h_muon_EtaVsPhi_"+cutVar.first]->Fill(GoodLeptons[goodLep].second.Eta(), GoodLeptons[goodLep].second.Phi(), weight);
                            }
                        }        
                    }          
 
                    my_histos["h_njets_"        +cutVar.first]->Fill( NGoodJets_pt30, weight  );   
                    my_histos["h_nbjets_"       +cutVar.first]->Fill( NGoodBJets_pt30, weight ); 
                    my_histos["h_met_"          +cutVar.first]->Fill( met, weight             );
                    my_histos["h_ht_"           +cutVar.first]->Fill( HT_trigger_pt30, weight );      
                    my_histos["h_jetPtMax_"     +cutVar.first]->Fill( jetPtMax, weight        ); // same for 0l & 1l
                    my_histos["h_lvMET_cm_mass_"+cutVar.first]->Fill( lvMET_cm_m, weight      );                    
                    my_histos["h_lvMET_cm_eta_" +cutVar.first]->Fill( lvMET_cm_eta, weight    );
                    my_histos["h_lvMET_cm_phi_" +cutVar.first]->Fill( lvMET_cm_phi, weight    );
                    my_histos["h_lvMET_cm_pt_"  +cutVar.first]->Fill( lvMET_cm_pt, weight     );
                    my_histos["h_fwm2_top6_"    +cutVar.first]->Fill( fwm2_top6, weight       );
                    my_histos["h_fwm3_top6_"    +cutVar.first]->Fill( fwm3_top6, weight       );
                    my_histos["h_fwm4_top6_"    +cutVar.first]->Fill( fwm4_top6, weight       );
                    my_histos["h_fwm5_top6_"    +cutVar.first]->Fill( fwm5_top6, weight       );
                    my_histos["h_fwm6_top6_"    +cutVar.first]->Fill( fwm6_top6, weight       );
                    my_histos["h_fwm7_top6_"    +cutVar.first]->Fill( fwm7_top6, weight       );
                    my_histos["h_fwm8_top6_"    +cutVar.first]->Fill( fwm8_top6, weight       );
                    my_histos["h_fwm9_top6_"    +cutVar.first]->Fill( fwm9_top6, weight       );
                    my_histos["h_fwm10_top6_"   +cutVar.first]->Fill( fwm10_top6, weight      );
                    my_histos["h_jmt_ev0_top6_" +cutVar.first]->Fill( jmt_ev0_top6, weight    );
                    my_histos["h_jmt_ev1_top6_" +cutVar.first]->Fill( jmt_ev1_top6, weight    );
                    my_histos["h_jmt_ev2_top6_" +cutVar.first]->Fill( jmt_ev2_top6, weight    );
                    my_histos["h_event_beta_z_" +cutVar.first]->Fill( event_beta_z, weight    );

                }

                // --------------- 
                // after HEM issue
                // --------------- 
                else 
                {
                    for (unsigned int goodJet = 0; goodJet < Jets.size(); goodJet++)
                    {

                        if (!GoodJets_pt30[goodJet]) { continue; }

                        my_histos["h_jetPt_HEM_"          +cutVar.first]->Fill( Jets[goodJet].Pt(), weight               );
                        my_2d_histos["h_jet_EtaVsPhi_HEM_"+cutVar.first]->Fill( Jets[goodJet].Eta(), Jets[goodJet].Phi(), weight ); 

                        if (GoodBJets_pt30[goodJet])
                            my_2d_histos["h_bjet_EtaVsPhi_HEM_"+cutVar.first]->Fill( Jets[goodJet].Eta(), Jets[goodJet].Phi(), weight ); 


                        if (Jets[goodJet].Pt() > jetPtMax)
                        { jetPtMax = Jets[goodJet].Pt(); }
                    }

                    if (cutVar.first.find("passBaseline_1l") != std::string::npos)
                    {
                        for (unsigned int goodLep = 0; goodLep < GoodLeptons.size(); goodLep++)
                        {
                            if (GoodLeptons[goodLep].first == "e")
                            {
                                my_2d_histos["h_electron_EtaVsPhi_HEM_"+cutVar.first]->Fill(GoodLeptons[goodLep].second.Eta(), GoodLeptons[goodLep].second.Phi(), weight);
                            }
                            else
                            {
                                my_2d_histos["h_muon_EtaVsPhi_HEM_"+cutVar.first]->Fill(GoodLeptons[goodLep].second.Eta(), GoodLeptons[goodLep].second.Phi(), weight);
                            }
                        }
                    }

                    my_histos["h_njets_HEM_"        +cutVar.first]->Fill( NGoodJets_pt30, weight  );
                    my_histos["h_nbjets_HEM_"       +cutVar.first]->Fill( NGoodBJets_pt30, weight );
                    my_histos["h_met_HEM_"          +cutVar.first]->Fill( met, weight             );
                    my_histos["h_ht_HEM_"           +cutVar.first]->Fill( HT_trigger_pt30, weight );
                    my_histos["h_jetPtMax_HEM_"     +cutVar.first]->Fill( jetPtMax, weight        ); // same for 0l & 1l                    
                    my_histos["h_lvMET_cm_mass_HEM_"+cutVar.first]->Fill( lvMET_cm_m, weight      );
                    my_histos["h_lvMET_cm_eta_HEM_" +cutVar.first]->Fill( lvMET_cm_eta, weight    );
                    my_histos["h_lvMET_cm_phi_HEM_" +cutVar.first]->Fill( lvMET_cm_phi, weight    );
                    my_histos["h_lvMET_cm_pt_HEM_"  +cutVar.first]->Fill( lvMET_cm_pt, weight     );
                    my_histos["h_fwm2_top6_HEM_"    +cutVar.first]->Fill( fwm2_top6, weight       );
                    my_histos["h_fwm3_top6_HEM_"    +cutVar.first]->Fill( fwm3_top6, weight       );
                    my_histos["h_fwm4_top6_HEM_"    +cutVar.first]->Fill( fwm4_top6, weight       );
                    my_histos["h_fwm5_top6_HEM_"    +cutVar.first]->Fill( fwm5_top6, weight       );
                    my_histos["h_fwm6_top6_HEM_"    +cutVar.first]->Fill( fwm6_top6, weight       );
                    my_histos["h_fwm7_top6_HEM_"    +cutVar.first]->Fill( fwm7_top6, weight       );
                    my_histos["h_fwm8_top6_HEM_"    +cutVar.first]->Fill( fwm8_top6, weight       );
                    my_histos["h_fwm9_top6_HEM_"    +cutVar.first]->Fill( fwm9_top6, weight       );
                    my_histos["h_fwm10_top6_HEM_"   +cutVar.first]->Fill( fwm10_top6, weight      );
                    my_histos["h_jmt_ev0_top6_HEM_" +cutVar.first]->Fill( jmt_ev0_top6, weight    );
                    my_histos["h_jmt_ev1_top6_HEM_" +cutVar.first]->Fill( jmt_ev1_top6, weight    );
                    my_histos["h_jmt_ev2_top6_HEM_" +cutVar.first]->Fill( jmt_ev2_top6, weight    );
                    my_histos["h_event_beta_z_HEM_" +cutVar.first]->Fill( event_beta_z, weight    );
                }

            }
        }
    } 
}

void HEM_Analyzer::WriteHistos(TFile* outfile)
{
    outfile->cd();

    for (const auto &p : my_histos) {
        p.second->Write();
    }
 
    for (const auto &p : my_2d_histos) {
        p.second->Write();
    }
}
