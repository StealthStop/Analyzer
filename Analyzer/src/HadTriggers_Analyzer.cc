#define HadTriggers_Analyzer_cxx
#include "Analyzer/Analyzer/include/HadTriggers_Analyzer.h"
#include "Framework/Framework/include/Utility.h"
#include "NTupleReader/include/NTupleReader.h"

#include <iostream>

HadTriggers_Analyzer::HadTriggers_Analyzer()
{
    InitHistos();
}

void HadTriggers_Analyzer::InitHistos()
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ;

    std::vector<std::string> effTags      { "den", "num"                                                                              };
    std::vector<std::string> combTrigTags { "jet"                                                                                     }; 
    std::vector<std::string> trigTags     { "trig"                                                                                    };
    std::vector<std::string> nBJetCutTags { "1bjetCut", "ge1bjetCut", "2bjetCut", "ge2bjetCut", "3bjetCut", "ge3bjetCut", "ge4bjetCut"};

    const int htbins   = 6;
    const int ptbins   = 4;
    const int bjetbins = 3;
    double htbinEdges[htbins + 1]     = {500, 550, 600, 650, 700, 800, 1000};
    double ptbinEdges[ptbins + 1 ]    = {45, 50, 55, 60, 90                };
    double bjetbinEdges[bjetbins + 1] = {1.5, 2.5, 3.5, 8                  };

    // ---------------------------------------------------------
    // latest triggers and preselctions with pt45 with nbjet cut
    // --------------------------------------------------------- 
    for( std::string effTag : effTags )
    {
        for( std::string combTrigTag : combTrigTags )
        {
            for( std::string trigTag : trigTags )
            {
                for( std::string nBJetCutTag : nBJetCutTags )
                {
                        // 1D - Efficiency
                        my_histos.emplace( "h_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_wJetHtBin", std::make_shared<TH1D>(("h_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_wJetHtBin").c_str(), ("h_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_wJetHtBin").c_str(), htbins, htbinEdges ) );
                        my_histos.emplace( "h_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_w6thJetPtBin", std::make_shared<TH1D>(("h_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_w6thJetPtBin").c_str(), ("h_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_w6thJetPtBin").c_str(), ptbins, ptbinEdges ) );
                        // 2D - Scale Factor
                        my_2d_histos.emplace( "h2_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_wJetHt6thJetPtBin", std::make_shared<TH2D>( ( "h2_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_wJetHt6thJetPtBin" ).c_str(), ( "h2_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_wJetHt6thJetPtBin" ).c_str(), htbins, htbinEdges, ptbins, ptbinEdges ) );

                }
            }
        }
    }   

} //


void HadTriggers_Analyzer::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    while( tr.getNextEvent() )
    {
        const auto& eventCounter = tr.getVar<int>("eventCounter");
        my_histos["EventCounter"]->Fill( eventCounter );

        const auto& runtype                = tr.getVar<std::string>("runtype");
        const auto& filetag                = tr.getVar<std::string>("filetag");
        const auto& Jets                   = tr.getVec<utility::LorentzVector>("Jets");
        const auto& GoodJets_pt45          = tr.getVec<bool>("GoodJets_pt45");
        const auto& NGoodBJets_pt45        = tr.getVar<int>("NGoodBJets_pt45");
        const auto& HT_trigger_pt45        = tr.getVar<double>("HT_trigger_pt45");
        const auto& passBaseline0l_trigEff = tr.getVar<bool>("passBaseline0l_trigEff");
        const auto& passTriggerMuon        = tr.getVar<bool>("passTriggerMuon");
        const auto& passTriggerAllHad      = tr.getVar<bool>("passTriggerAllHad");

        bool pass_1bjetCut   = NGoodBJets_pt45 == 1;
        bool pass_2bjetCut   = NGoodBJets_pt45 == 2;
        bool pass_ge2bjetCut = NGoodBJets_pt45 >= 2;
        bool pass_3bjetCut   = NGoodBJets_pt45 == 3;
        bool pass_ge3bjetCut = NGoodBJets_pt45 >= 3;
        bool pass_ge4bjetCut = NGoodBJets_pt45 >= 4;

        // -----------------------------------------
        // get the 6th jet pt for reference analysis
        // ----------------------------------------- 
        int njetspt45 = 0;
        double SixthJetPt45 = 0.0;
        for (unsigned int j = 0; j < Jets.size(); j++)
        {
            if (!GoodJets_pt45[j]) continue;
            njetspt45++;

            if (njetspt45 == 6)
            {
                SixthJetPt45 = Jets.at(j).Pt();
                break;
            }
        }
        
        // ------------------
        // Print Event Number
        // ------------------
        if(maxevents != -1 && tr.getEvtNum() >= maxevents) break;        
        //if( tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() ) ;

        // ----------------------
        // Print list of triggers
        // ----------------------
        //const auto& TriggerNames = tr.getVec<std::string>("TriggerNames");
        //if( tr.getEvtNum() == 1 ) printTriggerList(TriggerNames); 

        // -------------
        // Define weight
        // -------------
        double weight               = 1.0;
        double eventweight          = 1.0;
        double bTagScaleFactor      = 1.0;
        double prefiringScaleFactor = 1.0;
        double puScaleFactor        = 1.0;
        double topPtScaleFactor     = 1.0;
        if(runtype == "MC")
        {
            // Define Lumi weight
            const auto& Weight   = tr.getVar<float>("Weight");
            const auto& lumi     = tr.getVar<double>("FinalLumi");
            eventweight          = lumi*Weight;

            bTagScaleFactor      = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor");
            puScaleFactor        = tr.getVar<double>("puWeightCorr");
            topPtScaleFactor     = tr.getVar<double>("topPtScaleFactor");

            weight *= eventweight*bTagScaleFactor*prefiringScaleFactor*puScaleFactor*topPtScaleFactor;
        }
        
        // -------------------------------------------------------
        // Jet Trigger Efficiency on the SingleMuon Dataset and MC
        // -------------------------------------------------------
        if ( (filetag.find("Data_SingleMuon") != std::string::npos || runtype == "MC") )
        {
            const std::map<std::string, bool> cut_map_combHadMuTriggers
            {
                { "jet_trig_1bjetCut"  , passBaseline0l_trigEff && passTriggerMuon && pass_1bjetCut   },
                { "jet_trig_ge1bjetCut", passBaseline0l_trigEff && passTriggerMuon                    }, // preselection requires ge1b
                { "jet_trig_2bjetCut"  , passBaseline0l_trigEff && passTriggerMuon && pass_2bjetCut   },
                { "jet_trig_ge2bjetCut", passBaseline0l_trigEff && passTriggerMuon && pass_ge2bjetCut },
                { "jet_trig_3bjetCut"  , passBaseline0l_trigEff && passTriggerMuon && pass_3bjetCut   },
                { "jet_trig_ge3bjetCut", passBaseline0l_trigEff && passTriggerMuon && pass_ge3bjetCut },
                { "jet_trig_ge4bjetCut", passBaseline0l_trigEff && passTriggerMuon && pass_ge4bjetCut },

            };

            fillHistos(cut_map_combHadMuTriggers, passTriggerAllHad, HT_trigger_pt45, SixthJetPt45, weight);
        }
    }
}

void HadTriggers_Analyzer::WriteHistos(TFile* outfile)
{
    outfile->cd();

    for (const auto &p : my_histos) 
    {
        p.second->Write();
    }
    
    for (const auto &p : my_2d_histos) 
    {
        p.second->Write();
    }
}

bool HadTriggers_Analyzer::containsGoodHadron( const std::vector<utility::LorentzVector>& hadrons, const std::vector<bool>& goodHadrons, double ptThreshold, double etaSelection) 
{ 
    // Require a good hadron in JetHT data
    for( unsigned int h = 0; h < hadrons.size(); h++ ) 
    {
        if( !goodHadrons.at(h) ) continue; 
    
        utility::LorentzVector myHadron = hadrons.at(h);
    
        if( myHadron.Pt() >= ptThreshold && std::fabs( myHadron.Eta() ) < etaSelection ) return true;
    }

    return false;
}


void HadTriggers_Analyzer::fillHistos( const std::map<std::string, bool>& cutMap, bool passTriggerAllHad, double HT, double pt, double weight )
{ 
    for( auto& kv : cutMap )
    {
        if( kv.second )
        {
            my_histos["h_den_"+kv.first+"_wJetHtBin"]->Fill( HT, weight );
            my_histos["h_den_"+kv.first+"_w6thJetPtBin"]->Fill( pt, weight );
            my_2d_histos["h2_den_"+kv.first+"_wJetHt6thJetPtBin"]->Fill( HT, pt, weight );
            
            if( passTriggerAllHad )
            {
                my_histos["h_num_"+kv.first+"_wJetHtBin"]->Fill( HT, weight );
                my_histos["h_num_"+kv.first+"_w6thJetPtBin"]->Fill( pt, weight );
                my_2d_histos["h2_num_"+kv.first+"_wJetHt6thJetPtBin"]->Fill( HT, pt, weight );
            }

        }
    }
}


void HadTriggers_Analyzer::printTriggerList( const std::vector<std::string>& TriggerNames )
{
    for( unsigned int i = 0; i < TriggerNames.size(); i++ ) 
    {
        std::string myString = TriggerNames.at(i);
        printf("%s\n", myString.c_str());
    }
}
