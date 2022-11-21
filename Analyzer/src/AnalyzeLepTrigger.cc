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

    //Define some strings that are used for different scenarios that we want to calculate trigger efficiencies for
    std::vector<std::string> effTags        { "den", "num" }; 
    std::vector<std::string> lepTags        { "el", "mu" }; //Electron, muon, non-iso muon
    std::vector<std::string> ptTags         { "pt35", "pt40" }; //Pt threshold of other lepton
    std::vector<std::string> trigTags       { "trig", "noTrig" }; //Require other lepton trigger? Used for studies, but conclusion was that it does not affect statistics much so it's okay
    std::vector<std::string> nJetCutTags    { "1", "2", "3", "4", "5", "6" }; //Cutting at 6 jets really killed statistics, but the 1 jet cut region is not representative of our signal region, so look at all intermediate steps as well
 
    //Define binning for the histograms
    const Int_t nPtBins = 5;
    Double_t ptBinEdges[ nPtBins + 1 ] = { 30.0, 45.0, 60.0, 85.0, 120, 200 };
    const Int_t nEtaBins = 4;
    Double_t etaBinEdges[ nEtaBins + 1 ] = { -2.4, -1.4, 0, 1.4, 2.4 };

    for( std::string effTag : effTags ) {
        for( std::string lepTag : lepTags ) {
            for( std::string ptTag : ptTags ) {
                if( lepTag == "mu" && ptTag == "pt35" ) continue;
                for( std::string trigTag : trigTags ) {
                    for( std::string nJetCutTag : nJetCutTags ) {
                        
                        //Define 1D histograms - use these to get event yields ( a little more convenient than when extracting it from the TEfficiency )
                        my_histos.emplace( "h_trig_"+effTag+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"jCut_wLepPtBin", std::make_shared<TH1D>( ( "h_trig_"+effTag+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"jCut_wLepPtBin" ).c_str(), ( "h_trig_"+effTag+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"jCut_wLepPtBin" ).c_str(), nPtBins, ptBinEdges ) );
                        my_histos.emplace( "h_trig_"+effTag+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"jCut_wLepEtaBin", std::make_shared<TH1D>( ( "h_trig_"+effTag+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"jCut_wLepEtaBin" ).c_str(), ( "h_trig_"+effTag+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"jCut_wLepEtaBin" ).c_str(), nEtaBins, etaBinEdges ) );

                        //Define 2D histograms - use these for the final plots ( i.e., it is easier to get these values and use these to create TEfficiencies later on )
                        my_2d_histos.emplace( "h2_trig_"+effTag+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"jCut_wLepPtLepEtaBin", std::make_shared<TH2D>( ( "h2_trig_"+effTag+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"jCut_wLepPtLepEtaBin" ).c_str(), ( "h2_trig_"+effTag+"_"+lepTag+"_"+ptTag+"_"+trigTag+"_"+nJetCutTag+"jCut_wLepPtLepEtaBin" ).c_str(), nPtBins, ptBinEdges, nEtaBins, etaBinEdges ) );
                    }//End of nJetCutTags loop
                }//End of trigTags loop
            }//End of ptTags loop
        }//End of lepTags loop
    }//End of effTags loop
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
        const auto& runtype             = tr.getVar<std::string>("runtype");
        const auto& filetag             = tr.getVar<std::string>("filetag");
        const auto& etaCut              = tr.getVar<double>("etaCut");
        const auto& GoodLeptons         = tr.getVec<std::pair<std::string, utility::LorentzVector>>("GoodLeptons");

        const auto& NGoodLeptons        = tr.getVar<int>("NGoodLeptons");
        const auto& NGoodJets_pt30      = tr.getVar<int>("NGoodJets_pt30");

        const auto& Muons               = tr.getVec<utility::LorentzVector>("Muons");
        const auto& Electrons           = tr.getVec<utility::LorentzVector>("Electrons");

        const auto& NGoodMuons          = tr.getVar<int>("NGoodMuons");
        const auto& NGoodElectrons      = tr.getVar<int>("NGoodElectrons");
        
        const auto& passMadHT           = tr.getVar<bool>("passMadHT");
        const auto& passBaseline        = tr.getVar<bool>("passBaselineGoodOffline1l");

        const auto& GoodMuons           = tr.getVec<bool>("GoodMuons");
        const auto& GoodElectrons       = tr.getVec<bool>("GoodElectrons");
       
        //--------------------------------------------------
        //-- Print Event Number 
        //--------------------------------------------------
        if(maxevents != -1 && tr.getEvtNum() >= maxevents) break;        
        if( tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() ) ;

        //--------------------------------------------------
        //-- Print list of triggers (only if you want to see them)
        //--------------------------------------------------
        //const auto& TriggerNames        = tr.getVec<std::string>("TriggerNames");
        //if( tr.getEvtNum() == 1 ) printTriggerList(TriggerNames); 

        // ------------------------
        // -- Define theweight
        // ------------------------
        double theweight            = 1.0;
        double leptonScaleFactor    = 1.0;
        if(runtype == "MC")
        {
            if( !passMadHT ) continue; //Make sure not to double count DY events
            // Define Lumi weight
            const auto& Weight  = tr.getVar<float>("Weight");
            const auto& lumi = tr.getVar<double>("FinalLumi");
            const auto& puWeight = tr.getVar<double>("puWeightCorr");

            // Define lepton weight
            if(NGoodLeptons == 1)
            {
                const auto& eleLepWeight = tr.getVar<double>("noTrigGoodElectronSF");
                const auto& muLepWeight  = tr.getVar<double>("noTrigGoodMuonSF");
                leptonScaleFactor = (GoodLeptons[0].first == "e") ? eleLepWeight : muLepWeight;
            }

            const auto& bTagScaleFactor   = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            const auto& prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor");

            const auto& topPtScaleFactor  = tr.getVar<double>("topPtScaleFactor");

            theweight = lumi*Weight*puWeight*leptonScaleFactor*bTagScaleFactor*prefiringScaleFactor*topPtScaleFactor;
        }

        bool passMuonTriggers       = tr.getVar<bool>("passTriggerMuon");
        bool passElectronTriggers   = tr.getVar<bool>("passTriggerElectron");
        
        bool pass6JetCut            = ( NGoodJets_pt30 >= 6 );
        bool pass5JetCut            = ( NGoodJets_pt30 >= 5 );
        bool pass4JetCut            = ( NGoodJets_pt30 >= 4 );
        bool pass3JetCut            = ( NGoodJets_pt30 >= 3 );
        bool pass2JetCut            = ( NGoodJets_pt30 >= 2 );
        bool pass1JetCut            = ( NGoodJets_pt30 >= 1 ); //A little redundant because we require one good b jet but for symmetry reasons

        //-----------------------------------------------------------------------
        //-- Do the Electron Trigger Efficiency on the Single Muon Dataset and MC
        //-----------------------------------------------------------------------
        
        if( (filetag.find("SingleMuon") != std::string::npos || runtype == "MC") ) {

            if ( NGoodMuons >= 1 ) {
                bool foundMuonPt35 = containsGoodLepton(Muons, GoodMuons, 35, etaCut);
                bool foundMuonPt40 = containsGoodLepton(Muons, GoodMuons, 40, etaCut);
                
                //Look at the first good electron
                int myGoodElectronIndex = goodLeptonIndex(Electrons, GoodElectrons);
       
                if( myGoodElectronIndex != -1 ) {
                    const std::map<std::string, bool> cut_map_elTriggers {
                        { "el_pt35_trig_1jCut",   passBaseline && foundMuonPt35 && passMuonTriggers && pass1JetCut }, 
                        { "el_pt35_trig_2jCut",   passBaseline && foundMuonPt35 && passMuonTriggers && pass2JetCut }, 
                        { "el_pt35_trig_3jCut",   passBaseline && foundMuonPt35 && passMuonTriggers && pass3JetCut }, 
                        { "el_pt35_trig_4jCut",   passBaseline && foundMuonPt35 && passMuonTriggers && pass4JetCut }, 
                        { "el_pt35_trig_5jCut",   passBaseline && foundMuonPt35 && passMuonTriggers && pass5JetCut }, 
                        { "el_pt35_trig_6jCut",   passBaseline && foundMuonPt35 && passMuonTriggers && pass6JetCut }, 
                        
                        { "el_pt35_noTrig_1jCut", passBaseline && foundMuonPt35 && pass1JetCut }, 
                        { "el_pt35_noTrig_2jCut", passBaseline && foundMuonPt35 && pass2JetCut }, 
                        { "el_pt35_noTrig_3jCut", passBaseline && foundMuonPt35 && pass3JetCut }, 
                        { "el_pt35_noTrig_4jCut", passBaseline && foundMuonPt35 && pass4JetCut }, 
                        { "el_pt35_noTrig_5jCut", passBaseline && foundMuonPt35 && pass5JetCut }, 
                        { "el_pt35_noTrig_6jCut", passBaseline && foundMuonPt35 && pass6JetCut }, 
                        
                        { "el_pt40_trig_1jCut",   passBaseline && foundMuonPt40 && passMuonTriggers && pass1JetCut }, 
                        { "el_pt40_trig_2jCut",   passBaseline && foundMuonPt40 && passMuonTriggers && pass2JetCut }, 
                        { "el_pt40_trig_3jCut",   passBaseline && foundMuonPt40 && passMuonTriggers && pass3JetCut }, 
                        { "el_pt40_trig_4jCut",   passBaseline && foundMuonPt40 && passMuonTriggers && pass4JetCut }, 
                        { "el_pt40_trig_5jCut",   passBaseline && foundMuonPt40 && passMuonTriggers && pass5JetCut }, 
                        { "el_pt40_trig_6jCut",   passBaseline && foundMuonPt40 && passMuonTriggers && pass6JetCut }, 
                        
                        { "el_pt40_noTrig_1jCut", passBaseline && foundMuonPt40 && pass1JetCut }, 
                        { "el_pt40_noTrig_2jCut", passBaseline && foundMuonPt40 && pass2JetCut }, 
                        { "el_pt40_noTrig_3jCut", passBaseline && foundMuonPt40 && pass3JetCut }, 
                        { "el_pt40_noTrig_4jCut", passBaseline && foundMuonPt40 && pass4JetCut }, 
                        { "el_pt40_noTrig_5jCut", passBaseline && foundMuonPt40 && pass5JetCut }, 
                        { "el_pt40_noTrig_6jCut", passBaseline && foundMuonPt40 && pass6JetCut } 
                    };
    
                    fillHistos(cut_map_elTriggers, passElectronTriggers, Electrons.at( myGoodElectronIndex ), theweight);
                }
            }
        }
        
        //-----------------------------------------------------------------------
        //-- Do the Muon Trigger Efficiency on the Single Electron Dataset and MC
        //-----------------------------------------------------------------------
        
        if( (filetag.find("SingleElectron") != std::string::npos || runtype == "MC") ) {

            if ( NGoodElectrons >= 1 ) { 
                bool foundElectronPt40 = containsGoodLepton(Electrons, GoodElectrons, 40, etaCut);
                
                //Look at the first good muon
                int myGoodMuonIndex = goodLeptonIndex(Muons, GoodMuons);
                
                if( myGoodMuonIndex != -1 ) {
        
                    const std::map<std::string, bool> cut_map_muTriggers {
                        
                        { "mu_pt40_trig_1jCut",  passBaseline && foundElectronPt40 && passElectronTriggers && pass1JetCut }, 
                        { "mu_pt40_trig_2jCut",  passBaseline && foundElectronPt40 && passElectronTriggers && pass2JetCut }, 
                        { "mu_pt40_trig_3jCut",  passBaseline && foundElectronPt40 && passElectronTriggers && pass3JetCut }, 
                        { "mu_pt40_trig_4jCut",  passBaseline && foundElectronPt40 && passElectronTriggers && pass4JetCut }, 
                        { "mu_pt40_trig_5jCut",  passBaseline && foundElectronPt40 && passElectronTriggers && pass5JetCut }, 
                        { "mu_pt40_trig_6jCut",  passBaseline && foundElectronPt40 && passElectronTriggers && pass6JetCut }, 
                        
                        { "mu_pt40_noTrig_1jCut", passBaseline && foundElectronPt40 && pass1JetCut }, 
                        { "mu_pt40_noTrig_2jCut", passBaseline && foundElectronPt40 && pass2JetCut }, 
                        { "mu_pt40_noTrig_3jCut", passBaseline && foundElectronPt40 && pass3JetCut }, 
                        { "mu_pt40_noTrig_4jCut", passBaseline && foundElectronPt40 && pass4JetCut }, 
                        { "mu_pt40_noTrig_5jCut", passBaseline && foundElectronPt40 && pass5JetCut }, 
                        { "mu_pt40_noTrig_6jCut", passBaseline && foundElectronPt40 && pass6JetCut }, 
                        
                    };
                
                    fillHistos(cut_map_muTriggers, passMuonTriggers, Muons.at( myGoodMuonIndex ), theweight);
                }
            }
        }
    } 
}

void AnalyzeLepTrigger::WriteHistos(TFile* outfile)
{
    outfile->cd();

    for (const auto &p : my_histos) {
        p.second->Write();
    }
    
    for (const auto &p : my_2d_histos) {
        p.second->Write();
    }
}

bool AnalyzeLepTrigger::containsGoodLepton( const std::vector<utility::LorentzVector>& leptons, const std::vector<bool>& goodLeptons, double ptThreshold, double etaSelection) { 
    //Require a good muon in the single muon dataset
    for( unsigned int iLep = 0; iLep < leptons.size(); ++iLep ) {
        if( !goodLeptons.at( iLep ) ) continue; 
    
        utility::LorentzVector myLepton = leptons.at( iLep );
    
        if( myLepton.Pt() >= ptThreshold && std::fabs( myLepton.Eta() ) < etaSelection ) return true;
    }

    return false;
}

int AnalyzeLepTrigger::goodLeptonIndex( const std::vector<utility::LorentzVector>& leptons, const std::vector<bool>& goodLeptons) {
    for( unsigned int iLep = 0; iLep < leptons.size(); ++iLep ) {
        if( !goodLeptons.at( iLep ) ) continue;

        return iLep;
    }

    return -1;
}

void AnalyzeLepTrigger::fillHistos( const std::map<std::string, bool>& cutMap, bool passLeptonTriggers, const utility::LorentzVector& lepton, double theWeight ) {
    for( auto& kv : cutMap ) {
        if( kv.second ) {
            my_histos["h_trig_den_"+kv.first+"_wLepPtBin"]->Fill( lepton.Pt(), theWeight );
            my_histos["h_trig_den_"+kv.first+"_wLepEtaBin"]->Fill( lepton.Eta(), theWeight );
            my_2d_histos["h2_trig_den_"+kv.first+"_wLepPtLepEtaBin"]->Fill( lepton.Pt(), lepton.Eta(), theWeight );

            if( passLeptonTriggers ) {
                my_histos["h_trig_num_"+kv.first+"_wLepPtBin"]->Fill( lepton.Pt(), theWeight );
                my_histos["h_trig_num_"+kv.first+"_wLepEtaBin"]->Fill( lepton.Eta(), theWeight );
                my_2d_histos["h2_trig_num_"+kv.first+"_wLepPtLepEtaBin"]->Fill( lepton.Pt(), lepton.Eta(), theWeight );
            }
        }
    }
}

//Use this function to print out the entire trigger list in an ntuple (useful when trying to figure out which triggers to use)
void AnalyzeLepTrigger::printTriggerList( const std::vector<std::string>& TriggerNames ) {
    for( unsigned int i = 0; i < TriggerNames.size(); i++ ) {
        std::string myString = TriggerNames.at(i);
        printf("%s\n", myString.c_str());
    }
}
