#define AnalyzeSignalModels_cxx
#include "Analyzer/Analyzer/include/AnalyzeSignalModels.h"
#include "Framework/Framework/include/Utility.h"
#include "SusyAnaTools/Tools/NTupleReader.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TRandom3.h>
#include <iostream>
#include <TFile.h>

AnalyzeSignalModels::AnalyzeSignalModels()
{
    InitHistos();
}

//Define all your histograms here. 
void AnalyzeSignalModels::InitHistos()
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    //This event counter histogram is necessary so that we know that all the condor jobs ran successfully. If not, when you use the hadder script, you will see a discrepancy in red as the files are being hadded.
    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ;

    //Define 1D histograms
    my_2d_histos.emplace( "h_njets", std::make_shared<TH2D>( "h_njets", "h_njets", 13, 6.5, 19.5, 3, -0.5, 2.5 ) ) ;
    my_2d_histos.emplace( "h_met", std::make_shared<TH2D>( "h_met", "h_met", 720, 0, 1500, 3, -0.5, 2.5 ) ) ;
    my_2d_histos.emplace( "h_ht", std::make_shared<TH2D>( "h_ht", "h_ht", 720, 300, 5000, 3, -0.5, 2.5 ) ) ;
    my_2d_histos.emplace( "h_nb", std::make_shared<TH2D>( "h_nb", "h_nb", 15, -0.5, 14.5, 3, -0.5, 2.5 ) ) ;
    my_2d_histos.emplace( "h_jetptmax", std::make_shared<TH2D>( "h_jetptmax", "h_jetptmax", 1440, 0, 5000, 3, -0.5, 2.5 ) ) ;
    my_2d_histos.emplace( "h_jetpt", std::make_shared<TH2D>( "h_jetpt", "h_jetpt", 1440, 0, 5000, 3, -0.5, 2.5 ) ) ;

    my_2d_histos.emplace( "h_top_nlino_dphi", std::make_shared<TH2D>( "h_top_nlino_dphi", "h_top_nlino_dphi", 720, 0, 3.14, 3, -0.5, 2.5 ) ) ;
    my_2d_histos.emplace( "h_top_nlino_deta", std::make_shared<TH2D>( "h_top_nlino_deta", "h_top_nlino_deta", 720, 0, 3.14, 3, -0.5, 2.5 ) ) ;
    my_2d_histos.emplace( "h_atop_nlino_dphi", std::make_shared<TH2D>( "h_atop_nlino_dphi", "h_atop_nlino_dphi", 720, 0, 3.14, 3, -0.5, 2.5 ) ) ;
    my_2d_histos.emplace( "h_atop_nlino_deta", std::make_shared<TH2D>( "h_atop_nlino_deta", "h_atop_nlino_deta", 720, 0, 3.14, 3, -0.5, 2.5 ) ) ;
    my_2d_histos.emplace( "h_atop_top_dphi", std::make_shared<TH2D>( "h_atop_top_dphi", "h_atop_top_dphi", 720, 0, 3.14, 3, -0.5, 2.5 ) ) ;
    my_2d_histos.emplace( "h_atop_top_deta", std::make_shared<TH2D>( "h_atop_top_deta", "h_atop_top_deta", 720, 0, 3.14, 3, -0.5, 2.5 ) ) ;
    my_3d_histos.emplace( "h_toppt_top_nlino_dphi", std::make_shared<TH3D>( "h_toppt_top_nlino_dphi", "h_toppt_top_nlino_dphi", 120, 0, 1000, 120, 0, 3.14, 3, -0.5, 2.5 ) ) ;

    //Define TEfficiencies if you are doing trigger studies (for proper error bars) or cut flow charts.
    my_efficiencies.emplace("event_sel_weight", std::make_shared<TEfficiency>("event_sel_weight","event_sel_weight",9,0,9));
}

//Put everything you want to do per event here.
void AnalyzeSignalModels::Loop(NTupleReader& tr, double, int maxevents, bool)
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
        
        const auto& runtype         = tr.getVar<std::string>("runtype");     
        const auto& GoodJets        = tr.getVec<bool>("GoodJets");
        const auto& Jets            = tr.getVec<utility::LorentzVector>("Jets");

        const auto& GenParticles    = tr.getVec<utility::LorentzVector>("GenParticles");
        const auto& GenParticlesPDG = tr.getVec<int>("GenParticles_PdgId");
        const auto& GenParentsPDG = tr.getVec<int>("GenParticles_ParentId");

        //const auto& JetID           = tr.getVar<bool>("JetID");
        const auto& NGoodLeptons    = tr.getVar<int>("NGoodLeptons");
        const auto& NGoodBJets_pt30 = tr.getVar<int>("NGoodBJets_pt30");
        //const auto& Mbl             = tr.getVar<double>("Mbl");
        const auto& HT_trigger_pt30 = tr.getVar<double>("HT_trigger_pt30");
        const auto& NGoodJets       = tr.getVar<int>("NGoodJets");
        //const auto& NGoodJets_pt30  = tr.getVar<int>("NGoodJets_pt30");
        const auto& passMadHT       = tr.getVar<bool>("passMadHT");
        const auto& met             = tr.getVar<double>("MET");
       
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
            const auto& Weight  = tr.getVar<double>("Weight");
            const auto& lumi = tr.getVar<double>("Lumi");
            eventweight = lumi*Weight;
            
            // Define lepton weight
            //if(NGoodLeptons == 1)
            //{
            //    const auto& eleLepWeight = tr.getVar<double>("totGoodElectronSF");
            //    const auto& muLepWeight  = tr.getVar<double>("totGoodMuonSF");
            //    leptonScaleFactor = (GoodLeptons[0].first == "e") ? eleLepWeight : muLepWeight;
            //}
            
            //PileupWeight = tr.getVar<double>("_PUweightFactor");
            //bTagScaleFactor   = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            //htDerivedScaleFactor = tr.getVar<double>("htDerivedweight");
            //prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor");
            //puScaleFactor = tr.getVar<double>("puWeightCorr");
            
            weight *= eventweight*leptonScaleFactor*bTagScaleFactor*htDerivedScaleFactor*prefiringScaleFactor*puScaleFactor;
        }

        utility::LorentzVector stop; utility::LorentzVector astop; utility::LorentzVector top; utility::LorentzVector atop;
        utility::LorentzVector nlino; utility::LorentzVector anlino;
        for (unsigned int genPart = 0; genPart < GenParticles.size(); genPart++) {
            
            if (GenParticlesPDG[genPart] == 1000006) {
                stop = GenParticles[genPart];
            }
                
            else if (GenParticlesPDG[genPart] == -1000006) {
                astop = GenParticles[genPart];
            }

            else if (GenParticlesPDG[genPart] == 6) {
                top = GenParticles[genPart];
            }

            else if (GenParticlesPDG[genPart] == -6) {
                atop = GenParticles[genPart];
            }

            else if (GenParticlesPDG[genPart] == 1000022) {
                if (GenParentsPDG[genPart] == 1000006) {
                    nlino = GenParticles[genPart];
                }
                
                else if (GenParentsPDG[genPart] == -1000006) {
                    anlino = GenParticles[genPart];
                }
            }
        }

        if (nlino.E() != 0) {
            my_2d_histos["h_top_nlino_dphi"]->Fill(utility::DeltaPhi(top, nlino), NGoodLeptons, weight);
            my_2d_histos["h_top_nlino_deta"]->Fill(std::abs(top.Eta() - nlino.Eta()), NGoodLeptons, weight);
            my_2d_histos["h_atop_nlino_dphi"]->Fill(utility::DeltaPhi(atop, anlino), NGoodLeptons, weight);
            my_2d_histos["h_atop_nlino_deta"]->Fill(std::abs(atop.Eta() - anlino.Eta()), NGoodLeptons, weight);
            my_3d_histos["h_toppt_top_nlino_dphi"]->Fill(stop.Pt(), utility::DeltaPhi(top, nlino), NGoodLeptons, weight) ;

        }
        my_2d_histos["h_atop_top_dphi"]->Fill(utility::DeltaPhi(atop, top), NGoodLeptons, weight);
        my_2d_histos["h_atop_top_deta"]->Fill(std::abs(atop.Eta() - top.Eta()), NGoodLeptons, weight);

        double jetPtMax = 0.0;
        for (unsigned int goodJet = 0; goodJet < GoodJets.size(); goodJet++) {
            if (!GoodJets[goodJet]) { continue; }
            my_2d_histos["h_jetpt"]->Fill(Jets[goodJet].Pt(), NGoodLeptons, weight);
            if (Jets[goodJet].Pt() > jetPtMax) { jetPtMax = Jets[goodJet].Pt(); }
        }
        
        //Make cuts and fill histograms here
        //if( passBaseline ) {
            my_2d_histos["h_njets"]->Fill( NGoodJets, NGoodLeptons, weight); 
            my_2d_histos["h_met"]->Fill(met, NGoodLeptons, weight);
            my_2d_histos["h_ht"]->Fill(HT_trigger_pt30, NGoodLeptons, weight);
            my_2d_histos["h_nb"]->Fill(NGoodBJets_pt30, NGoodLeptons, weight);
            my_2d_histos["h_jetptmax"]->Fill(jetPtMax, NGoodLeptons, weight);
        //}
    } 
}

void AnalyzeSignalModels::WriteHistos(TFile* outfile)
{
    outfile->cd();

    for (const auto &p : my_histos) {
        p.second->Write();
    }
    
    for (const auto &p : my_2d_histos) {
        p.second->Write();
    }

    for (const auto &p : my_3d_histos) {
        p.second->Write();
    }
    
    for (const auto &p : my_efficiencies) {
        p.second->Write();
    }
    
}
