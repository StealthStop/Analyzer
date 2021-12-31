#define HadTriggers_Analyzer_cxx
#include "Analyzer/Analyzer/include/HadTriggers_Analyzer.h"
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

HadTriggers_Analyzer::HadTriggers_Analyzer()
{
    InitHistos();
}

void HadTriggers_Analyzer::InitHistos()
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ;

    std::vector<std::string> effTags      { "denominator", "numerator"                                            };
    std::vector<std::string> combTrigTags { "CombHadIsoMu"                                                        }; 
    std::vector<std::string> trigTags     { "trig", "noTrig"                                                      };
    std::vector<std::string> nBJetCutTags { "ge2bjetCut", "2bjetCut", "3bjetCut", "ge4bjetCut"                    };
    std::vector<std::string> nJetCutTags  { "ge6jetCut", "6jetCut", "7jetCut", "8jetCut", "9jetCut", "ge10jetCut" };
    std::vector<std::string> ptTags       { "pt45"                                                                }; 

    const int htbins   = 9;
    const int ptbins   = 6;
    const int njetbins = 9;
    const int bjetbins = 3;
    double htbinEdges[htbins + 1]     = {500, 550, 600, 650, 700, 800, 1000, 1500, 2000, 2500};
    double ptbinEdges[ptbins + 1 ]    = {45, 50, 55, 60, 70, 120, 200};
    double njetbinEdges[njetbins + 1] = {6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    double bjetbinEdges[bjetbins + 1] = {1.5, 2.5, 3.5, 8};

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
                    for( std::string ptTag : ptTags )
                    {
                        my_histos.emplace( "h_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_"+ptTag+"_HT", std::make_shared<TH1D>(("h_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_"+ptTag+"_HT").c_str(), ("h_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_"+ptTag+"_HT").c_str(), htbins, htbinEdges ) );
                        my_histos.emplace( "h_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_"+ptTag+"_6thJetPt", std::make_shared<TH1D>(("h_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_"+ptTag+"_6thJetPt").c_str(), ("h_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_"+ptTag+"_6thJetPt").c_str(), ptbins, ptbinEdges ) );
                        my_histos.emplace( "h_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_"+ptTag+"_NJet", std::make_shared<TH1D>(("h_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_"+ptTag+"_NJet").c_str(), ("h_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_"+ptTag+"_NJet").c_str(), njetbins, njetbinEdges ) );
                        my_histos.emplace( "h_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_"+ptTag+"_NBJet", std::make_shared<TH1D>(("h_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_"+ptTag+"_NBJet").c_str(), ("h_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_"+ptTag+"_NBJet").c_str(), bjetbins, bjetbinEdges ) );

                        my_2d_histos.emplace( "h_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_"+ptTag+"_HTvs6thJetPt", std::make_shared<TH2D>( ( "h_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_"+ptTag+"_HTvs6thJetPt" ).c_str(), ( "h_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_"+ptTag+"_HTvs6thJetPt" ).c_str(), htbins, htbinEdges, ptbins, ptbinEdges ) );

                    }
                }
            }
        }
    }   

    // --------------------------------------------------------
    // latest triggers and preselctions with pt45 with njet cut
    // --------------------------------------------------------
    for( std::string effTag : effTags )
    {
        for( std::string combTrigTag : combTrigTags )
        {
            for( std::string trigTag : trigTags )
            {
                for( std::string nJetCutTag : nJetCutTags )
                {
                    for( std::string ptTag : ptTags )
                    {
                        my_histos.emplace( "h_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nJetCutTag+"_"+ptTag+"_HT", std::make_shared<TH1D>(("h_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nJetCutTag+"_"+ptTag+"_HT").c_str(), ("h_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nJetCutTag+"_"+ptTag+"_HT").c_str(), htbins, htbinEdges ) );
                        my_histos.emplace( "h_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nJetCutTag+"_"+ptTag+"_6thJetPt", std::make_shared<TH1D>(("h_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nJetCutTag+"_"+ptTag+"_6thJetPt").c_str(), ("h_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nJetCutTag+"_"+ptTag+"_6thJetPt").c_str(), ptbins, ptbinEdges ) );
                        my_histos.emplace( "h_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nJetCutTag+"_"+ptTag+"_NJet", std::make_shared<TH1D>(("h_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nJetCutTag+"_"+ptTag+"_NJet").c_str(), ("h_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nJetCutTag+"_"+ptTag+"_NJet").c_str(), njetbins, njetbinEdges ) );                        
                        my_histos.emplace( "h_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nJetCutTag+"_"+ptTag+"_NBJet", std::make_shared<TH1D>(("h_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nJetCutTag+"_"+ptTag+"_NBJet").c_str(), ("h_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nJetCutTag+"_"+ptTag+"_NBJet").c_str(), bjetbins, bjetbinEdges ) );

                        my_2d_histos.emplace( "h_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nJetCutTag+"_"+ptTag+"_HTvs6thJetPt", std::make_shared<TH2D>( ( "h_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nJetCutTag+"_"+ptTag+"_HTvs6thJetPt" ).c_str(), ( "h_"+effTag+"_"+combTrigTag+"_"+trigTag+"_"+nJetCutTag+"_"+ptTag+"_HTvs6thJetPt" ).c_str(), htbins, htbinEdges, ptbins, ptbinEdges ) );

                    }
                }
            }
        }
    }

} //


void HadTriggers_Analyzer::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    while( tr.getNextEvent() )
    {
        const auto& eventCounter        = tr.getVar<int>("eventCounter");
        my_histos["EventCounter"]->Fill( eventCounter );

        const auto& runtype               = tr.getVar<std::string>("runtype");
        const auto& filetag               = tr.getVar<std::string>("filetag");
        const auto& Jets                  = tr.getVec<utility::LorentzVector>("Jets");
        const auto& GoodJets_pt45         = tr.getVec<bool>("GoodJets_pt45");
        const auto& NGoodJets_pt45        = tr.getVar<int>("NGoodJets_pt45");
        const auto& NGoodBJets_pt45       = tr.getVar<int>("NGoodBJets_pt45");
        const auto& HT_trigger_pt45       = tr.getVar<double>("HT_trigger_pt45");
        const auto& passBaseline0l_pt45   = tr.getVar<bool>("passBaseline0l_pt45");
        const auto& passTriggerMuonsRefAN = tr.getVar<bool>("passTriggerMuonsRefAN");
        const auto& passTriggerAllHad     = tr.getVar<bool>("passTriggerAllHad");

        bool pass_2bjetCut   = NGoodBJets_pt45 == 2;
        bool pass_3bjetCut   = NGoodBJets_pt45 == 3;
        bool pass_ge4bjetCut = NGoodBJets_pt45 >= 4;
        bool pass_6jetCut    = NGoodJets_pt45 == 6;
        bool pass_7jetCut    = NGoodJets_pt45 == 7;
        bool pass_8jetCut    = NGoodJets_pt45 == 8;
        bool pass_9jetCut    = NGoodJets_pt45 == 9;
        bool pass_ge10jetCut = NGoodJets_pt45 >= 10;

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
        
        // ------------------------
        // -- Print Event Number 
        // ------------------------
        if(maxevents != -1 && tr.getEvtNum() >= maxevents) break;        
        //if( tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() ) ;

        // ----------------------------
        // -- Print list of triggers 
        // ----------------------------
        //const auto& TriggerNames = tr.getVec<std::string>("TriggerNames");
        //if( tr.getEvtNum() == 1 ) printTriggerList(TriggerNames); 

        // -------------------
        // -- Define weight
        // -------------------
        double weight               = 1.0;
        double eventweight          = 1.0;
        double bTagScaleFactor      = 1.0;
        double prefiringScaleFactor = 1.0;
        double puScaleFactor        = 1.0;
        if(runtype == "MC")
        {
            // Define Lumi weight
            const auto& Weight   = tr.getVar<float>("Weight");
            const auto& lumi     = tr.getVar<double>("Lumi");
            eventweight          = lumi*Weight;

            bTagScaleFactor      = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor");
            puScaleFactor        = tr.getVar<double>("puWeightCorr");

            weight *= eventweight*bTagScaleFactor*prefiringScaleFactor*puScaleFactor;
        }
        
        // ---------------------------------------------------------
        // -- Trigger Efficiency on the SingleMuon Dataset and MC
        // ---------------------------------------------------------
        if ( (filetag.find("Data_SingleMuon") != std::string::npos || runtype == "MC") )
        {
            // ---------------------------------------------------------
            // latest triggers and preselctions with pt45 with nbjet cut
            // ---------------------------------------------------------
            const std::map<std::string, bool> cut_map_combHadMuTriggers
            {
                { "CombHadIsoMu_trig_ge2bjetCut_pt45",   passBaseline0l_pt45 && passTriggerMuonsRefAN && passTriggerAllHad                    },
                { "CombHadIsoMu_trig_2bjetCut_pt45",     passBaseline0l_pt45 && passTriggerMuonsRefAN && passTriggerAllHad && pass_2bjetCut   },
                { "CombHadIsoMu_trig_3bjetCut_pt45",     passBaseline0l_pt45 && passTriggerMuonsRefAN && passTriggerAllHad && pass_3bjetCut   },
                { "CombHadIsoMu_trig_ge4bjetCut_pt45",   passBaseline0l_pt45 && passTriggerMuonsRefAN && passTriggerAllHad && pass_ge4bjetCut },

                { "CombHadIsoMu_noTrig_ge2bjetCut_pt45", passBaseline0l_pt45 && passTriggerMuonsRefAN                    },
                { "CombHadIsoMu_noTrig_2bjetCut_pt45",   passBaseline0l_pt45 && passTriggerMuonsRefAN && pass_2bjetCut   },
                { "CombHadIsoMu_noTrig_3bjetCut_pt45",   passBaseline0l_pt45 && passTriggerMuonsRefAN && pass_3bjetCut   },
                { "CombHadIsoMu_noTrig_ge4bjetCut_pt45", passBaseline0l_pt45 && passTriggerMuonsRefAN && pass_ge4bjetCut },
            };
            fillHistosRefAN(cut_map_combHadMuTriggers, passTriggerAllHad, HT_trigger_pt45, SixthJetPt45, NGoodJets_pt45, NGoodBJets_pt45, weight);

            // --------------------------------------------------------
            // latest triggers and preselctions with pt45 with njet cut
            // --------------------------------------------------------
            const std::map<std::string, bool> cut_map2_combHadMuTriggers
            {
                { "CombHadIsoMu_trig_ge6jetCut_pt45",    passBaseline0l_pt45 && passTriggerMuonsRefAN && passTriggerAllHad                    },    
                { "CombHadIsoMu_trig_6jetCut_pt45",      passBaseline0l_pt45 && passTriggerMuonsRefAN && passTriggerAllHad && pass_6jetCut    },
                { "CombHadIsoMu_trig_7jetCut_pt45",      passBaseline0l_pt45 && passTriggerMuonsRefAN && passTriggerAllHad && pass_7jetCut    },
                { "CombHadIsoMu_trig_8jetCut_pt45",      passBaseline0l_pt45 && passTriggerMuonsRefAN && passTriggerAllHad && pass_8jetCut    },
                { "CombHadIsoMu_trig_9jetCut_pt45",      passBaseline0l_pt45 && passTriggerMuonsRefAN && passTriggerAllHad && pass_9jetCut    },
                { "CombHadIsoMu_trig_ge10jetCut_pt45",   passBaseline0l_pt45 && passTriggerMuonsRefAN && passTriggerAllHad && pass_ge10jetCut },

                { "CombHadIsoMu_noTrig_ge6jetCut_pt45",  passBaseline0l_pt45 && passTriggerMuonsRefAN                    },
                { "CombHadIsoMu_noTrig_6jetCut_pt45",    passBaseline0l_pt45 && passTriggerMuonsRefAN && pass_6jetCut    },
                { "CombHadIsoMu_noTrig_7jetCut_pt45",    passBaseline0l_pt45 && passTriggerMuonsRefAN && pass_7jetCut    },
                { "CombHadIsoMu_noTrig_8jetCut_pt45",    passBaseline0l_pt45 && passTriggerMuonsRefAN && pass_8jetCut    },
                { "CombHadIsoMu_noTrig_9jetCut_pt45",    passBaseline0l_pt45 && passTriggerMuonsRefAN && pass_9jetCut    },
                { "CombHadIsoMu_noTrig_ge10jetCut_pt45", passBaseline0l_pt45 && passTriggerMuonsRefAN && pass_ge10jetCut },

            };
            fillHistosRefAN(cut_map2_combHadMuTriggers, passTriggerAllHad, HT_trigger_pt45, SixthJetPt45, NGoodJets_pt45, NGoodBJets_pt45, weight); 

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
    
    for (const auto &p : my_efficiencies) 
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

void HadTriggers_Analyzer::fillHistos( const std::map<std::string, bool>& cutMap, bool passTriggerAllHad, double HT, int njet, int nbjet, double weight ) 
{
    for( auto& kv : cutMap ) 
    {
        if( kv.second ) 
        {
            my_histos["h_denominator_"+kv.first+"_HT"]->Fill( HT, weight );
            my_histos["h_denominator_"+kv.first+"_ht5000"]->Fill( HT, weight );
            my_histos["h_denominator_"+kv.first+"_NJet"]->Fill( njet, weight );
            my_histos["h_denominator_"+kv.first+"_NBJet"]->Fill( nbjet, weight );
            my_2d_histos["h_denominator_"+kv.first+"_NJetVsHT"]->Fill( njet, HT, weight );
            my_2d_histos["h_denominator_"+kv.first+"_NJetVsHt"]->Fill( njet, HT, weight );
            my_2d_histos["h_denominator_"+kv.first+"_NJetVsNBJet"]->Fill( njet, nbjet, weight );

            if( passTriggerAllHad ) 
            {
                my_histos["h_numerator_"+kv.first+"_HT"]->Fill( HT, weight );
                my_histos["h_numerator_"+kv.first+"_ht5000"]->Fill( HT, weight );
                my_histos["h_numerator_"+kv.first+"_NJet"]->Fill( njet, weight );
                my_histos["h_numerator_"+kv.first+"_NBJet"]->Fill( nbjet, weight );
                my_2d_histos["h_numerator_"+kv.first+"_NJetVsHT"]->Fill( njet, HT, weight );
                my_2d_histos["h_numerator_"+kv.first+"_NJetVsHt"]->Fill( njet, HT, weight );
                my_2d_histos["h_numerator_"+kv.first+"_NJetVsNBJet"]->Fill( njet, nbjet, weight );
            }
        }
    }
}

// function for the reference analysis histos
void HadTriggers_Analyzer::fillHistosRefAN( const std::map<std::string, bool>& cutMap, bool passTriggerRefAN, double HT, double pt, int njet, int nbjet, double weight )
{ 
    for( auto& kv : cutMap )
    {
        if( kv.second )
        {
            my_histos["h_denominator_"+kv.first+"_HT"]->Fill( HT, weight );
            my_histos["h_denominator_"+kv.first+"_6thJetPt"]->Fill( pt, weight );
            my_histos["h_denominator_"+kv.first+"_NJet"]->Fill( njet, weight );
            my_histos["h_denominator_"+kv.first+"_NBJet"]->Fill( nbjet, weight );
            my_2d_histos["h_denominator_"+kv.first+"_HTvs6thJetPt"]->Fill( HT, pt, weight );
            
            if( passTriggerRefAN )
            {
                my_histos["h_numerator_"+kv.first+"_HT"]->Fill( HT, weight );
                my_histos["h_numerator_"+kv.first+"_6thJetPt"]->Fill( pt, weight );
                my_histos["h_numerator_"+kv.first+"_NJet"]->Fill( njet, weight );
                my_histos["h_numerator_"+kv.first+"_NBJet"]->Fill( nbjet, weight );
                my_2d_histos["h_numerator_"+kv.first+"_HTvs6thJetPt"]->Fill( HT, pt, weight );
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
