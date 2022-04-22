#define AnalyzeHEM_cxx
#include "Analyzer/Analyzer/include/AnalyzeHEM.h"
#include "Framework/Framework/include/Utility.h"
#include "NTupleReader/include/NTupleReader.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TRandom3.h>
#include <iostream>
#include <TFile.h>

AnalyzeHEM::AnalyzeHEM()
{
    InitHistos();
}

//Define all your histograms here. 
void AnalyzeHEM::InitHistos()
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    //This event counter histogram is necessary so that we know that all the condor jobs ran successfully. If not, when you use the hadder script, you will see a discrepancy in red as the files are being hadded.
    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ;
    my_histos.emplace( "h_njets", std::make_shared<TH1D>( "h_njets", "h_njets", 13, 6.5, 19.5 ) ) ;
    my_histos.emplace( "h_ngoodleptons", std::make_shared<TH1D>( "h_ngoodleptons", "h_ngoodleptons", 5, 0.5, 5.5 ) ) ;
    my_histos.emplace( "h_met", std::make_shared<TH1D>( "h_met", "h_met", 720, 0, 1500) ) ;
    my_histos.emplace( "h_ht", std::make_shared<TH1D>( "h_ht", "h_ht", 720, 300, 5000) ) ;
    my_histos.emplace( "h_nb", std::make_shared<TH1D>( "h_nb", "h_nb", 15, -0.5, 14.5) ) ;
    my_histos.emplace( "h_jetptmax", std::make_shared<TH1D>( "h_jetptmax", "h_jetptmax", 1440, 0, 5000) ) ;
    my_histos.emplace( "h_jetpt", std::make_shared<TH1D>( "h_jetpt", "h_jetpt", 1440, 0, 5000) ) ;
    my_histos.emplace( "h_lvMET_cm_pt",std::make_shared<TH1D>( "h_lvMET_cm_pt", "h_lvMET_cm_pt", 1440, 0, 1000) ) ;
    my_histos.emplace( "h_lvMET_cm_eta",std::make_shared<TH1D>( "h_lvMET_cm_eta", "h_lvMET_cm_eta", 1440, -2.5, 2.5) ) ;
    my_histos.emplace( "h_lvMET_cm_phi",std::make_shared<TH1D>( "h_lvMET_cm_phi", "h_lvMET_cm_phi", 1440, -3.14, 3.14) ) ;
    my_histos.emplace( "h_lvMET_cm_m",std::make_shared<TH1D>( "h_lvMET_cm_m", "h_lvMET_cm_m", 1440, 0, 1000) ) ;
    my_histos.emplace( "h_fwm2_top6",std::make_shared<TH1D>( "h_fwm2_top6", "h_fwm2_top6", 1440, 0, 1) ) ;
    my_histos.emplace( "h_fwm3_top6",std::make_shared<TH1D>( "h_fwm3_top6", "h_fwm3_top6", 1440, 0, 1) ) ;
    my_histos.emplace( "h_fwm4_top6",std::make_shared<TH1D>( "h_fwm4_top6", "h_fwm4_top6", 1440, 0, 1) ) ;
    my_histos.emplace( "h_fwm5_top6",std::make_shared<TH1D>( "h_fwm5_top6", "h_fwm5_top6", 1440, 0, 1) ) ;
    my_histos.emplace( "h_fwm6_top6",std::make_shared<TH1D>( "h_fwm6_top6", "h_fwm6_top6", 1440, 0, 1) ) ;
    my_histos.emplace( "h_fwm7_top6",std::make_shared<TH1D>( "h_fwm7_top6", "h_fwm7_top6", 1440, 0, 1) ) ;
    my_histos.emplace( "h_fwm8_top6",std::make_shared<TH1D>( "h_fwm8_top6", "h_fwm8_top6", 1440, 0, 1) ) ;
    my_histos.emplace( "h_fwm9_top6",std::make_shared<TH1D>( "h_fwm9_top6", "h_fwm9_top6", 1440, 0, 1) ) ;
    my_histos.emplace( "h_fwm10_top6",std::make_shared<TH1D>( "h_fwm10_top6", "h_fwm10_top6", 1440, 0, 1) ) ;
    my_histos.emplace( "h_jmt_ev0_top6",std::make_shared<TH1D>( "h_jmt_ev0_top6", "h_jmt_ev0_top6", 1440, 0, 1) ) ;
    my_histos.emplace( "h_jmt_ev1_top6",std::make_shared<TH1D>( "h_jmt_ev1_top6", "h_jmt_ev1_top6", 1440, 0, 1) ) ;
    my_histos.emplace( "h_jmt_ev2_top6",std::make_shared<TH1D>( "h_jmt_ev2_top6", "h_jmt_ev2_top6", 1440, 0, 1) ) ;
    my_histos.emplace( "h_event_beta_z", std::make_shared<TH1D>( "h_beta_z", "h_beta_z", 1440, -1, 1) ) ;

    my_histos.emplace( "h_njets_HEM", std::make_shared<TH1D>( "h_njets_HEM", "h_njets_HEM", 13, 6.5, 19.5) ) ;
    my_histos.emplace( "h_ngoodleptons_HEM", std::make_shared<TH1D>( "h_ngoodleptons_HEM", "h_ngoodleptons_HEM", 5, 0.5, 5.5 ) ) ;
    my_histos.emplace( "h_met_HEM", std::make_shared<TH1D>( "h_met_HEM", "h_met_HEM", 720, 0, 1500) ) ;
    my_histos.emplace( "h_ht_HEM", std::make_shared<TH1D>( "h_ht_HEM", "h_ht_HEM", 720, 300, 5000) ) ;
    my_histos.emplace( "h_nb_HEM", std::make_shared<TH1D>( "h_nb_HEM", "h_nb_HEM", 15, -0.5, 14.5) ) ;
    my_histos.emplace( "h_jetptmax_HEM", std::make_shared<TH1D>( "h_jetptmax_HEM", "h_jetptmax_HEM", 1440, 0, 5000) ) ;
    my_histos.emplace( "h_jetpt_HEM", std::make_shared<TH1D>( "h_jetpt_HEM", "h_jetpt_HEM", 1440, 0, 5000) ) ;
    my_histos.emplace( "h_lvMET_cm_pt_HEM",std::make_shared<TH1D>( "h_lvMET_cm_pt_HEM", "h_lvMET_cm_pt_HEM", 1440, 0, 1000) ) ;
    my_histos.emplace( "h_lvMET_cm_eta_HEM",std::make_shared<TH1D>( "h_lvMET_cm_eta_HEM", "h_lvMET_cm_eta_HEM", 1440, -2.5, 2.5) ) ;
    my_histos.emplace( "h_lvMET_cm_phi_HEM",std::make_shared<TH1D>( "h_lvMET_cm_phi_HEM", "h_lvMET_cm_phi_HEM", 1440, -3.14, 3.14) ) ;
    my_histos.emplace( "h_lvMET_cm_m_HEM",std::make_shared<TH1D>( "h_lvMET_cm_m_HEM", "h_lvMET_cm_m_HEM", 1440, 0, 1000) ) ;
    my_histos.emplace( "h_fwm2_top6_HEM",std::make_shared<TH1D>( "h_fwm2_top6_HEM", "h_fwm2_top6_HEM", 1440, 0, 1) ) ;
    my_histos.emplace( "h_fwm3_top6_HEM",std::make_shared<TH1D>( "h_fwm3_top6_HEM", "h_fwm3_top6_HEM", 1440, 0, 1) ) ;
    my_histos.emplace( "h_fwm4_top6_HEM",std::make_shared<TH1D>( "h_fwm4_top6_HEM", "h_fwm4_top6_HEM", 1440, 0, 1) ) ;
    my_histos.emplace( "h_fwm5_top6_HEM",std::make_shared<TH1D>( "h_fwm5_top6_HEM", "h_fwm5_top6_HEM", 1440, 0, 1) ) ;
    my_histos.emplace( "h_fwm6_top6_HEM",std::make_shared<TH1D>( "h_fwm6_top6_HEM", "h_fwm6_top6_HEM", 1440, 0, 1) ) ;
    my_histos.emplace( "h_fwm7_top6_HEM",std::make_shared<TH1D>( "h_fwm7_top6_HEM", "h_fwm7_top6_HEM", 1440, 0, 1) ) ;
    my_histos.emplace( "h_fwm8_top6_HEM",std::make_shared<TH1D>( "h_fwm8_top6_HEM", "h_fwm8_top6_HEM", 1440, 0, 1) ) ;
    my_histos.emplace( "h_fwm9_top6_HEM",std::make_shared<TH1D>( "h_fwm9_top6_HEM", "h_fwm9_top6_HEM", 1440, 0, 1) ) ;
    my_histos.emplace( "h_fwm10_top6_HEM",std::make_shared<TH1D>( "h_fwm10_top6_HEM", "h_fwm10_top6_HEM", 1440, 0, 1) ) ;
    my_histos.emplace( "h_jmt_ev0_top6_HEM",std::make_shared<TH1D>( "h_jmt_ev0_top6_HEM", "h_jmt_ev0_top6_HEM", 1440, 0, 1) ) ;
    my_histos.emplace( "h_jmt_ev1_top6_HEM",std::make_shared<TH1D>( "h_jmt_ev1_top6_HEM", "h_jmt_ev1_top6_HEM", 1440, 0, 1) ) ;
    my_histos.emplace( "h_jmt_ev2_top6_HEM",std::make_shared<TH1D>( "h_jmt_ev2_top6_HEM", "h_jmt_ev2_top6_HEM", 1440, 0, 1) ) ;
    my_histos.emplace( "h_event_beta_z_HEM", std::make_shared<TH1D>( "h_beta_z_HEM", "h_beta_z_HEM", 1440, -1, 1) ) ;

    my_2d_histos.emplace( "h_electron_etaphi", std::make_shared<TH2D>( "h_electron_etaphi", "h_electron_etaphi", 1440, -2.5, 2.5, 1440, -3.14, 3.14) );
    my_2d_histos.emplace( "h_muon_etaphi", std::make_shared<TH2D>( "h_muon_etaphi", "h_muon_etaphi", 1440, -2.5, 2.5, 1440, -3.14, 3.14) );

    my_2d_histos.emplace( "h_electron_etaphi_HEM", std::make_shared<TH2D>( "h_electron_etaphi_HEM", "h_electron_etaphi_HEM", 1440, -2.5, 2.5, 1440, -3.14, 3.14) );
    my_2d_histos.emplace( "h_muon_etaphi_HEM", std::make_shared<TH2D>( "h_muon_etaphi_HEM", "h_muon_etaphi_HEM", 1440, -2.5, 2.5, 1440, -3.14, 3.14) );

}

//Put everything you want to do per event here.
void AnalyzeHEM::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    while( tr.getNextEvent() )
    {

        //This is added to count the number of events- do not change the next two lines.
        const auto& eventCounter        = tr.getVar<int>("eventCounter");
        my_histos["EventCounter"]->Fill( eventCounter );

        //--------------------------------------------------
        //-- Print Event Number 
        //--------------------------------------------------
        
        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() & (10000 == 0) ) printf( " Event %i\n", tr.getEvtNum() );
       
        const auto& runtype             = tr.getVar<std::string>("runtype");     
        const auto& GoodLeptons         = tr.getVec<std::pair<std::string, utility::LorentzVector>>("GoodLeptons");
        const auto& Jets                = tr.getVec<utility::LorentzVector>("Jets");
        const auto& GoodJets_pt30       = tr.getVec<bool>("GoodJets_pt30");
        const auto& NGoodLeptons        = tr.getVar<int>("NGoodLeptons");
        const auto& NGoodBJets_pt30     = tr.getVar<int>("NGoodBJets_pt30");
        const auto& HT_trigger_pt30     = tr.getVar<double>("HT_trigger_pt30");
        const auto& NGoodJets_pt30      = tr.getVar<int>("NGoodJets_pt30");
        const auto& RunNum              = tr.getVar<unsigned int>("RunNum");
        const auto& met                 = tr.getVar<float>("MET");
        const auto& lvMET_cm_pt         = tr.getVar<double>("lvMET_cm_pt");
        const auto& lvMET_cm_eta        = tr.getVar<double>("lvMET_cm_eta");
        const auto& lvMET_cm_phi        = tr.getVar<double>("lvMET_cm_phi");
        const auto& lvMET_cm_m          = tr.getVar<double>("lvMET_cm_m");
        const auto& fwm2_top6           = tr.getVar<double>("fwm2_top6");
        const auto& fwm3_top6           = tr.getVar<double>("fwm3_top6");
        const auto& fwm4_top6           = tr.getVar<double>("fwm4_top6");
        const auto& fwm5_top6           = tr.getVar<double>("fwm5_top6");
        const auto& fwm6_top6           = tr.getVar<double>("fwm6_top6");
        const auto& fwm7_top6           = tr.getVar<double>("fwm7_top6");
        const auto& fwm8_top6           = tr.getVar<double>("fwm8_top6");
        const auto& fwm9_top6           = tr.getVar<double>("fwm9_top6");
        const auto& fwm10_top6          = tr.getVar<double>("fwm10_top6");
        const auto& jmt_ev0_top6        = tr.getVar<double>("jmt_ev0_top6");
        const auto& jmt_ev1_top6        = tr.getVar<double>("jmt_ev1_top6");
        const auto& jmt_ev2_top6        = tr.getVar<double>("jmt_ev2_top6");
        const auto& event_beta_z        = tr.getVar<double>("event_beta_z");
        const auto& passMadHT           = tr.getVar<bool>("passMadHT");
        const auto& passBaseline_1l     = tr.getVar<bool>("passBaseline1l_Good");
        const auto& passBaseline_0l     = tr.getVar<bool>("passBaseline0l_good");
     
        if(maxevents != -1 && tr.getEvtNum() >= maxevents) break;        
        if ( tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() ) ;

        // ------------------------
        // -- Define weight
        // ------------------------
        double weight               = 1.0;
        double eventweight          = 1.0;
        double leptonScaleFactor    = 1.0;
        double bTagScaleFactor      = 1.0;
        double htDerivedScaleFactor = 1.0;
        double prefiringScaleFactor = 1.0;
        double puScaleFactor        = 1.0;
        
        if(runtype == "MC")
        {
            if( !passMadHT ) continue; //Make sure not to double count DY events
            // Define Lumi weight
            const auto& Weight = tr.getVar<float>("Weight");
            const auto& lumi   = tr.getVar<double>("Lumi");
            eventweight        = lumi*Weight;
            
            // Define lepton weight
            //if(NGoodLeptons == 1)
            //{
            //    const auto& eleLepWeight = tr.getVar<double>("totGoodElectronSF");
            //    const auto& muLepWeight  = tr.getVar<double>("totGoodMuonSF");
            //    leptonScaleFactor        = (GoodLeptons[0].first == "e") ? eleLepWeight : muLepWeight;
            //}
            
            //PileupWeight         = tr.getVar<float>("_PUweightFactor");
            //bTagScaleFactor      = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            //htDerivedScaleFactor = tr.getVar<double>("htDerivedweight");
            //prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor");
            //puScaleFactor        = tr.getVar<double>("puWeightCorr");
            
            //weight *= eventweight*leptonScaleFactor*bTagScaleFactor*htDerivedScaleFactor*prefiringScaleFactor*puScaleFactor;
            weight *= eventweight;
        }

        //Make cuts and fill histograms here
        if( passBaseline_1l ) {

            double jetPtMax = 0.0;
            if ( RunNum < 319077 ) {
               
                for (unsigned int goodJet = 0; goodJet < Jets.size(); goodJet++) {
                    if (!GoodJets_pt30[goodJet]) { continue; }
                        my_histos["h_jetpt"]->Fill(Jets[goodJet].Pt(), weight);
                    if (Jets[goodJet].Pt() > jetPtMax) { jetPtMax = Jets[goodJet].Pt(); }
                }

                for (unsigned int goodLep = 0; goodLep < GoodLeptons.size(); goodLep++) {
                    if (GoodLeptons[goodLep].first == "e") {
                        my_2d_histos["h_electron_etaphi"]->Fill(GoodLeptons[goodLep].second.Eta(), GoodLeptons[goodLep].second.Phi());
                    } else {
                        my_2d_histos["h_muon_etaphi"]->Fill(GoodLeptons[goodLep].second.Eta(), GoodLeptons[goodLep].second.Phi());
                    }
                }

                my_histos["h_njets"]->Fill( NGoodJets_pt30, weight); 
                my_histos["h_ngoodleptons"]->Fill( NGoodLeptons, weight);
                my_histos["h_met"]->Fill(met, weight);
                my_histos["h_ht"]->Fill(HT_trigger_pt30, weight);
                my_histos["h_nb"]->Fill(NGoodBJets_pt30, weight);
                my_histos["h_jetptmax"]->Fill(jetPtMax, weight);
                my_histos["h_lvMET_cm_pt"]->Fill(lvMET_cm_pt, weight);
                my_histos["h_lvMET_cm_eta"]->Fill(lvMET_cm_eta, weight);
                my_histos["h_lvMET_cm_phi"]->Fill(lvMET_cm_phi, weight);
                my_histos["h_lvMET_cm_m"]->Fill(lvMET_cm_m, weight);
                my_histos["h_fwm2_top6"]->Fill(fwm2_top6, weight);
                my_histos["h_fwm3_top6"]->Fill(fwm3_top6, weight);
                my_histos["h_fwm4_top6"]->Fill(fwm4_top6, weight);
                my_histos["h_fwm5_top6"]->Fill(fwm5_top6, weight);
                my_histos["h_fwm6_top6"]->Fill(fwm6_top6, weight);
                my_histos["h_fwm7_top6"]->Fill(fwm7_top6, weight);
                my_histos["h_fwm8_top6"]->Fill(fwm8_top6, weight);
                my_histos["h_fwm9_top6"]->Fill(fwm9_top6, weight);
                my_histos["h_fwm10_top6"]->Fill(fwm10_top6, weight);
                my_histos["h_jmt_ev0_top6"]->Fill(jmt_ev0_top6, weight);
                my_histos["h_jmt_ev1_top6"]->Fill(jmt_ev1_top6, weight);
                my_histos["h_jmt_ev2_top6"]->Fill(jmt_ev2_top6, weight);
                my_histos["h_event_beta_z"]->Fill(event_beta_z, weight);

            } else {
                for (unsigned int goodJet = 0; goodJet < Jets.size(); goodJet++) {
                    if (!GoodJets_pt30[goodJet]) { continue; }
                        my_histos["h_jetpt_HEM"]->Fill(Jets[goodJet].Pt(), weight);
                    if (Jets[goodJet].Pt() > jetPtMax) { jetPtMax = Jets[goodJet].Pt(); }
                }
               
                for (unsigned int goodLep = 0; goodLep < GoodLeptons.size(); goodLep++) {
                    if (GoodLeptons[goodLep].first == "e") {
                        my_2d_histos["h_electron_etaphi_HEM"]->Fill(GoodLeptons[goodLep].second.Eta(), GoodLeptons[goodLep].second.Phi());
                    } else {
                        my_2d_histos["h_muon_etaphi_HEM"]->Fill(GoodLeptons[goodLep].second.Eta(), GoodLeptons[goodLep].second.Phi()); 
                    }
                }

                my_histos["h_njets_HEM"]->Fill( NGoodJets_pt30, weight); 
                my_histos["h_ngoodleptons_HEM"]->Fill( NGoodLeptons, weight);
                my_histos["h_met_HEM"]->Fill(met, weight);
                my_histos["h_ht_HEM"]->Fill(HT_trigger_pt30, weight);
                my_histos["h_nb_HEM"]->Fill(NGoodBJets_pt30, weight);
                my_histos["h_jetptmax_HEM"]->Fill(jetPtMax, weight);
                my_histos["h_lvMET_cm_pt_HEM"]->Fill(lvMET_cm_pt, weight);
                my_histos["h_lvMET_cm_eta_HEM"]->Fill(lvMET_cm_eta, weight);
                my_histos["h_lvMET_cm_phi_HEM"]->Fill(lvMET_cm_phi, weight);
                my_histos["h_lvMET_cm_m_HEM"]->Fill(lvMET_cm_m, weight);
                my_histos["h_fwm2_top6_HEM"]->Fill(fwm2_top6, weight);
                my_histos["h_fwm3_top6_HEM"]->Fill(fwm3_top6, weight);
                my_histos["h_fwm4_top6_HEM"]->Fill(fwm4_top6, weight);
                my_histos["h_fwm5_top6_HEM"]->Fill(fwm5_top6, weight);
                my_histos["h_fwm6_top6_HEM"]->Fill(fwm6_top6, weight);
                my_histos["h_fwm7_top6_HEM"]->Fill(fwm7_top6, weight);
                my_histos["h_fwm8_top6_HEM"]->Fill(fwm8_top6, weight);
                my_histos["h_fwm9_top6_HEM"]->Fill(fwm9_top6, weight);
                my_histos["h_fwm10_top6_HEM"]->Fill(fwm10_top6, weight);
                my_histos["h_jmt_ev0_top6_HEM"]->Fill(jmt_ev0_top6, weight);
                my_histos["h_jmt_ev1_top6_HEM"]->Fill(jmt_ev1_top6, weight);
                my_histos["h_jmt_ev2_top6_HEM"]->Fill(jmt_ev2_top6, weight);
                my_histos["h_event_beta_z_HEM"]->Fill(event_beta_z, weight);

            }
        }
    } 
}

void AnalyzeHEM::WriteHistos(TFile* outfile)
{
    outfile->cd();

    for (const auto &p : my_histos) {
        p.second->Write();
    }
    
    for (const auto &p : my_2d_histos) {
        p.second->Write();
    }
}
