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

    std::vector<std::string> effTags  { "denominator", "numerator" }; 
    std::vector<std::string> hadTags  { "had", "had_IsoMu"         };
    std::vector<std::string> trigTags { "trig", "noTrig"           };
    std::vector<std::string> ptTags   { "pt45"                     }; // label for GoodJets_pt45 & GoodBJets_pt45 & HT_trigger_pt45 
 
    const int nHTbins   = 7;
    const int nhtBins   = 13;
    const int nJetBins  = 9;
    const int nBJetBins = 5;
    double HTbinEdges[nHTbins + 1]      = {400, 500, 600, 700, 800, 900, 1000, 1100};
    double htBinEdges[nhtBins + 1 ]     = {0, 200, 400, 600, 800, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000};
    double njetBinEdges[nJetBins + 1]   = {6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    double nbjetBinEdges[nBJetBins + 1] = {0, 1, 2, 3, 4, 5};

    for( std::string effTag : effTags ) 
    {
        for( std::string hadTag : hadTags ) 
        {
            for( std::string trigTag : trigTags )
            {
                    for( std::string ptTag : ptTags ) 
                    {
                        my_histos.emplace( "h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_HT", std::make_shared<TH1D>(("h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_HT").c_str(), ("h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_HT").c_str(), nHTbins, HTbinEdges ) );
                        my_histos.emplace( "h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_ht5000", std::make_shared<TH1D>(("h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_ht5000").c_str(), ("h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_ht5000").c_str(), nhtBins, htBinEdges ) );
                        my_histos.emplace( "h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_NJet", std::make_shared<TH1D>(("h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_NJet").c_str(), ("h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_NJet").c_str(), nJetBins, njetBinEdges ) );
                        my_histos.emplace( "h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_NBJet", std::make_shared<TH1D>(("h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_NBJet").c_str(), ("h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_NBJet").c_str(), nBJetBins, nbjetBinEdges ) );                                

                        my_2d_histos.emplace( "h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_NJetVsHT", std::make_shared<TH2D>( ( "h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_NJetVsHT" ).c_str(), ( "h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_NJetVsHT" ).c_str(), nJetBins, njetBinEdges, nHTbins, HTbinEdges ) );
                        my_2d_histos.emplace( "h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_NJetVsHt", std::make_shared<TH2D>( ( "h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_NJetVsHt" ).c_str(), ( "h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_NJetVsHt" ).c_str(), nJetBins, njetBinEdges, nhtBins, htBinEdges ) );
                        my_2d_histos.emplace( "h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_NJetVsNBJet", std::make_shared<TH2D>( ( "h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_NJetVsNBJet" ).c_str(), ( "h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_NJetVsNBJet" ).c_str(), nJetBins, njetBinEdges, nBJetBins, nbjetBinEdges ) );

                }
            }
        }
    }
}

void HadTriggers_Analyzer::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    while( tr.getNextEvent() )
    {
        const auto& eventCounter        = tr.getVar<int>("eventCounter");
        my_histos["EventCounter"]->Fill( eventCounter );

        const auto& runtype                  = tr.getVar<std::string>("runtype");
        const auto& filetag                  = tr.getVar<std::string>("filetag");
        const auto& NGoodJets_pt45           = tr.getVar<int>("NGoodJets_pt45");
        const auto& NGoodBJets_pt45          = tr.getVar<int>("NGoodBJets_pt45");
        const auto& HT_trigger_pt45          = tr.getVar<double>("HT_trigger_pt45");
        const auto& passTriggerAllHad        = tr.getVar<bool>("passTriggerAllHad");     
        //const auto& passIsoMuTrigger         = tr.getVar<bool>("passIsoMuTrigger"); 
        const auto& passBaseline0l_hadTrig   = tr.getVar<bool>("passBaseline0l_hadTrig");
        const auto& passBaseline0l_hadMuTrig = tr.getVar<bool>("passBaseline0l_hadMuTrig"); 
       
        // ------------------------
        // -- Print Event Number 
        // ------------------------
        if(maxevents != -1 && tr.getEvtNum() >= maxevents) break;        
        //if( tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() ) ;

        // ----------------------------
        // -- Print list of triggers 
        // ----------------------------
        const auto& TriggerNames = tr.getVec<std::string>("TriggerNames");
        if( tr.getEvtNum() == 1 ) printTriggerList(TriggerNames); 

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
            const auto& Weight   = tr.getVar<double>("Weight");
            const auto& lumi     = tr.getVar<double>("Lumi");
            eventweight          = lumi*Weight;

            bTagScaleFactor      = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor");
            puScaleFactor        = tr.getVar<double>("puWeightCorr");

            weight *= eventweight*bTagScaleFactor*prefiringScaleFactor*puScaleFactor;
        }

        // ----------------------------------------------------
        // -- Trigger Efficiency on the JetHT Dataset and MC
        // ----------------------------------------------------
        
        if( (filetag.find("Data_JetHT") != std::string::npos || runtype == "MC") ) 
        {
            const std::map<std::string, bool> cut_map_hadTriggers 
            {
                { "had_trig_pt45",   passBaseline0l_hadTrig && passTriggerAllHad },                                               
                { "had_noTrig_pt45", passBaseline0l_hadTrig }, 
            };
        
            fillHistos(cut_map_hadTriggers, passTriggerAllHad, HT_trigger_pt45, NGoodJets_pt45, NGoodBJets_pt45, weight);
        }

        if ( (filetag.find("Data_SingleMuon") != std::string::npos || runtype == "MC") )
        {
            const std::map<std::string, bool> cut_map_hadMuTriggers
            {   
                { "had_IsoMu_trig_pt45",   passBaseline0l_hadMuTrig && passTriggerAllHad},
                { "had_IsoMu_noTrig_pt45", passBaseline0l_hadMuTrig },  
            };

            fillHistos(cut_map_hadMuTriggers, passTriggerAllHad, HT_trigger_pt45, NGoodJets_pt45, NGoodBJets_pt45, weight);
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

bool HadTriggers_Analyzer::containsGoodHadron( const std::vector<TLorentzVector>& hadrons, const std::vector<bool>& goodHadrons, double ptThreshold, double etaSelection) 
{ 
    // Require a good hadron in JetHT data
    for( unsigned int h = 0; h < hadrons.size(); h++ ) 
    {
        if( !goodHadrons.at(h) ) continue; 
    
        TLorentzVector myHadron = hadrons.at(h);
    
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

void HadTriggers_Analyzer::printTriggerList( const std::vector<std::string>& TriggerNames )
{
    for( unsigned int i = 0; i < TriggerNames.size(); i++ ) 
    {
        std::string myString = TriggerNames.at(i);
        printf("%s\n", myString.c_str());
    }
}
