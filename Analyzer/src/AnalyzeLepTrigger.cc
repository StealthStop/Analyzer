#define AnalyzeLepTrigger_cxx
#include "Analyzer/Analyzer/include/AnalyzeLepTrigger.h"
#include "Framework/Framework/include/Utility.h"
#include "NTupleReader/include/NTupleReader.h"

#include <iostream>

AnalyzeLepTrigger::AnalyzeLepTrigger()
{
    InitHistos();
}

//Define all your histograms here. 
void AnalyzeLepTrigger::InitHistos()
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    //This event counter histogram is necessary so that we know that all the condor jobs ran successfully. If not, when you use the hadder script, you will see a discrepancy in red as the files are being hadded.
    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ;

    // Define some strings that are used for different scenarios that we want to calculate trigger efficiencies for
    std::vector<std::string> effTags    { "den", "num"                                                   }; // eff = den / num
    std::vector<std::string> lepTags    { "el", "mu"                                                     }; // Electron, muon
    std::vector<std::string> ptTags     { "pt40"                                                         }; // Pt threshold 
    std::vector<std::string> trigTags   { "trig"                                                         }; 
    std::vector<std::string> nJetCutTags{ "ge1jetCut", "ge2jetCut", "ge3jetCut", "ge4jetCut", "ge5jetCut"}; 
 
    // Define binning for the histograms
    const Int_t nPtBins = 5;
    Double_t ptBinEdges[ nPtBins + 1 ] = { 30.0, 45.0, 60.0, 85.0, 120, 200 };
    const Int_t nEtaBins = 4;
    Double_t etaBinEdges[ nEtaBins + 1 ] = { -2.4, -1.4, 0, 1.4, 2.4 };

    for( std::string effTag : effTags ) 
    {
        for( std::string lepTag : lepTags ) 
        {
            for( std::string ptTag : ptTags ) 
            {
                for( std::string trigTag : trigTags ) 
                {
                    for( std::string nJetCutTag : nJetCutTags ) 
                    {
                        // 1D - Efficiency  
                        my_histos.emplace( "h_"+effTag+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"_wLepPtBin", std::make_shared<TH1D>( ( "h_"+effTag+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"_wLepPtBin" ).c_str(), ( "h_"+effTag+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"_wLepPtBin" ).c_str(), nPtBins, ptBinEdges ) );
                        my_histos.emplace( "h_"+effTag+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"_wLepEtaBin", std::make_shared<TH1D>( ( "h_"+effTag+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"_wLepEtaBin" ).c_str(), ( "h_"+effTag+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"_wLepEtaBin" ).c_str(), nEtaBins, etaBinEdges ) );

                        // 2D - Scale Factor
                        my_2d_histos.emplace( "h2_"+effTag+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"_wLepPtLepEtaBin", std::make_shared<TH2D>( ( "h2_"+effTag+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"_wLepPtLepEtaBin" ).c_str(), ( "h2_"+effTag+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"_wLepPtLepEtaBin" ).c_str(), nPtBins, ptBinEdges, nEtaBins, etaBinEdges ) );
                    }
                }
            }
        }
    }
}

//Put everything you want to do per event here.
void AnalyzeLepTrigger::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    while( tr.getNextEvent() )
    {
        //This is added to count the number of events- do not change the next two lines.
        const auto& eventCounter        = tr.getVar<int>("eventCounter");
        my_histos["EventCounter"]->Fill( eventCounter );

        //Define useful variables here
        const auto& runtype                = tr.getVar<std::string>("runtype");
        const auto& filetag                = tr.getVar<std::string>("filetag");
        const auto& etaCut                 = tr.getVar<double>("etaCut");
        const auto& GoodLeptons            = tr.getVec<std::pair<std::string, utility::LorentzVector>>("GoodLeptons");
        const auto& NGoodLeptons           = tr.getVar<int>("NGoodLeptons");
        const auto& NGoodJets_pt30         = tr.getVar<int>("NGoodJets_pt30");
        const auto& Muons                  = tr.getVec<utility::LorentzVector>("Muons");
        const auto& Electrons              = tr.getVec<utility::LorentzVector>("Electrons");
        const auto& NGoodMuons             = tr.getVar<int>("NGoodMuons");
        const auto& NGoodElectrons         = tr.getVar<int>("NGoodElectrons");
        const auto& passMadHT              = tr.getVar<bool>("passMadHT");
        const auto& GoodMuons              = tr.getVec<bool>("GoodMuons");
        const auto& GoodElectrons          = tr.getVec<bool>("GoodElectrons");
        const auto& passBaseline1l_trigEff = tr.getVar<bool>("passBaseline1l_trigEff");      

        bool passMuonTriggers      = tr.getVar<bool>("passTriggerMuon");
        bool passElectronTriggers  = tr.getVar<bool>("passTriggerElectron");

        bool pass_ge1JetCut = ( NGoodJets_pt30 >= 1 );
        bool pass_ge2JetCut = ( NGoodJets_pt30 >= 2 );
        bool pass_ge3JetCut = ( NGoodJets_pt30 >= 3 );
        bool pass_ge4JetCut = ( NGoodJets_pt30 >= 4 );
        bool pass_ge5JetCut = ( NGoodJets_pt30 >= 5 );

        // ------------------
        // Print Event Number
        // ------------------
        if(maxevents != -1 && tr.getEvtNum() >= maxevents) break;        
        if( tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() ) ;

        // ------------------------------------------------------
        //  Print list of triggers (only if you want to see them)
        // ------------------------------------------------------
        //const auto& TriggerNames        = tr.getVec<std::string>("TriggerNames");
        //if( tr.getEvtNum() == 1 ) printTriggerList(TriggerNames); 

        // ----------------
        // Define theweight
        // ----------------
        double theweight         = 1.0;
        double leptonScaleFactor = 1.0;
        if(runtype == "MC")
        {
            // Define Lumi weight
            const auto& Weight   = tr.getVar<float>("Weight");
            const auto& lumi     = tr.getVar<double>("FinalLumi");
            const auto& puWeight = tr.getVar<double>("puWeightCorr");

            // Define lepton weight
            if(NGoodLeptons == 1)
            {
                const auto& eleLepWeight = tr.getVar<double>("noTrigGoodElectronSF");
                const auto& muLepWeight  = tr.getVar<double>("noTrigGoodMuonSF");
                leptonScaleFactor = (GoodLeptons[0].first == "e") ? eleLepWeight : muLepWeight;
            }

            const auto& bTagScaleFactor      = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            const auto& prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor");
            const auto& topPtScaleFactor     = tr.getVar<double>("topPtScaleFactor");

            theweight = lumi*Weight*puWeight*leptonScaleFactor*bTagScaleFactor*prefiringScaleFactor*topPtScaleFactor;
        }


        // --------------------------------------------------------------------
        // Do the Electron Trigger Efficiency on the Single Muon Dataset and MC
        // --------------------------------------------------------------------
        if( (filetag.find("SingleMuon") != std::string::npos || runtype == "MC") ) 
        {
            if ( NGoodMuons >= 1 ) 
            {
                bool foundMuonPt40 = containsGoodLepton(Muons, GoodMuons, 40, etaCut);
                
                // Look at the first good electron
                int myGoodElectronIndex = goodLeptonIndex(Electrons, GoodElectrons);
       
                if( myGoodElectronIndex != -1 ) 
                {
                    const std::map<std::string, bool> cut_map_elTriggers 
                    {  
                        { "el_pt40_trig_ge1jetCut",   passBaseline1l_trigEff && foundMuonPt40 && passMuonTriggers                   }, // preselection requires ge1jetCut 
                        { "el_pt40_trig_ge2jetCut",   passBaseline1l_trigEff && foundMuonPt40 && passMuonTriggers && pass_ge2JetCut }, 
                        { "el_pt40_trig_ge3jetCut",   passBaseline1l_trigEff && foundMuonPt40 && passMuonTriggers && pass_ge3JetCut }, 
                        { "el_pt40_trig_ge4jetCut",   passBaseline1l_trigEff && foundMuonPt40 && passMuonTriggers && pass_ge4JetCut }, 
                        { "el_pt40_trig_ge5jetCut",   passBaseline1l_trigEff && foundMuonPt40 && passMuonTriggers && pass_ge5JetCut }, 
    
                    };
    
                    fillHistos(cut_map_elTriggers, passElectronTriggers, Electrons.at( myGoodElectronIndex ), theweight);
                }
            }
        } // 
        
        // -------------------------------------------------------------
        // Muon Trigger Efficiency on the Single Electron Dataset and MC
        // -------------------------------------------------------------
        
        if( (filetag.find("SingleElectron") != std::string::npos || runtype == "MC") ) 
        {
            if ( NGoodElectrons >= 1 ) 
            { 
                bool foundElectronPt40 = containsGoodLepton(Electrons, GoodElectrons, 40, etaCut);
                
                // Look at the first good muon
                int myGoodMuonIndex = goodLeptonIndex(Muons, GoodMuons);
                
                if( myGoodMuonIndex != -1 ) 
                {
                    const std::map<std::string, bool> cut_map_muTriggers 
                    {
                        { "mu_pt40_trig_ge1jetCut",  passBaseline1l_trigEff && foundElectronPt40 && passElectronTriggers                   }, // ge1jetCut
                        { "mu_pt40_trig_ge2jetCut",  passBaseline1l_trigEff && foundElectronPt40 && passElectronTriggers && pass_ge2JetCut }, 
                        { "mu_pt40_trig_ge3jetCut",  passBaseline1l_trigEff && foundElectronPt40 && passElectronTriggers && pass_ge3JetCut }, 
                        { "mu_pt40_trig_ge4jetCut",  passBaseline1l_trigEff && foundElectronPt40 && passElectronTriggers && pass_ge4JetCut }, 
                        { "mu_pt40_trig_ge5jetCut",  passBaseline1l_trigEff && foundElectronPt40 && passElectronTriggers && pass_ge5JetCut }, 
                       
                    };
                
                    fillHistos(cut_map_muTriggers, passMuonTriggers, Muons.at( myGoodMuonIndex ), theweight);
                }
            }
        } //
    } 
}

void AnalyzeLepTrigger::WriteHistos(TFile* outfile)
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

bool AnalyzeLepTrigger::containsGoodLepton( const std::vector<utility::LorentzVector>& leptons, const std::vector<bool>& goodLeptons, double ptThreshold, double etaSelection) 
{     
    // Require a good muon in the single muon dataset
    for( unsigned int iLep = 0; iLep < leptons.size(); ++iLep ) 
    {
        if( !goodLeptons.at( iLep ) ) continue; 
    
        utility::LorentzVector myLepton = leptons.at( iLep );
    
        if( myLepton.Pt() >= ptThreshold && std::fabs( myLepton.Eta() ) < etaSelection ) return true;
    }

    return false;
}

int AnalyzeLepTrigger::goodLeptonIndex( const std::vector<utility::LorentzVector>& leptons, const std::vector<bool>& goodLeptons) 
{
    for( unsigned int iLep = 0; iLep < leptons.size(); ++iLep ) 
    {
        if( !goodLeptons.at( iLep ) ) continue;

        return iLep;
    }

    return -1;
}

void AnalyzeLepTrigger::fillHistos( const std::map<std::string, bool>& cutMap, bool passLeptonTriggers, const utility::LorentzVector& lepton, double theWeight ) 
{
    for( auto& kv : cutMap ) 
    {
        if( kv.second ) 
        {
            my_histos["h_den_"+kv.first+"_wLepPtBin"]->Fill( lepton.Pt(), theWeight );
            my_histos["h_den_"+kv.first+"_wLepEtaBin"]->Fill( lepton.Eta(), theWeight );
            my_2d_histos["h2_den_"+kv.first+"_wLepPtLepEtaBin"]->Fill( lepton.Pt(), lepton.Eta(), theWeight );

            if( passLeptonTriggers ) 
            {
                my_histos["h_num_"+kv.first+"_wLepPtBin"]->Fill( lepton.Pt(), theWeight );
                my_histos["h_num_"+kv.first+"_wLepEtaBin"]->Fill( lepton.Eta(), theWeight );
                my_2d_histos["h2_num_"+kv.first+"_wLepPtLepEtaBin"]->Fill( lepton.Pt(), lepton.Eta(), theWeight );
            }
        }
    }
}

// Use this function to print out the entire trigger list in an ntuple (useful when trying to figure out which triggers to use)
void AnalyzeLepTrigger::printTriggerList( const std::vector<std::string>& TriggerNames ) 
{
    for( unsigned int i = 0; i < TriggerNames.size(); i++ ) 
    {
        std::string myString = TriggerNames.at(i);
        printf("%s\n", myString.c_str());
    }
}
