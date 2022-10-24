#define Semra_Analyzer_cxx
#include "Analyzer/Analyzer/include/Semra_Analyzer.h"
#include "Framework/Framework/include/Utility.h"
#include "NTupleReader/include/NTupleReader.h"

#include <iostream>

Semra_Analyzer::Semra_Analyzer() : inithisto(false) // define inithisto variable
{
}

// -------------------
// -- Define histos
// -------------------
void Semra_Analyzer::InitHistos(const std::map<std::string, bool>& cutmap) // define variable map
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    
    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ;
    
    for (const auto& cutVar : cutmap) 
    {  
        my_histos.emplace( "h_njets_"+cutVar.first,         std::make_shared<TH1D> ( ("h_njets_"+cutVar.first).c_str(),         ("h_njets_"+cutVar.first).c_str(),         20, 0, 20     ) );
        my_histos.emplace( "h_nbjets_"+cutVar.first,        std::make_shared<TH1D> ( ("h_nbjets_"+cutVar.first).c_str(),        ("h_nbjets_"+cutVar.first).c_str(),        20, 0, 20     ) );
        
        my_histos.emplace( "h_ntops_"+cutVar.first,         std::make_shared<TH1D> ( ("h_ntops_"+cutVar.first).c_str(),         ("h_ntops_"+cutVar.first).c_str(),         10, 0, 10     ) );
        my_histos.emplace( "h_nRtops_"+cutVar.first,        std::make_shared<TH1D> ( ("h_nRtops_"+cutVar.first).c_str(),        ("h_nRtops_"+cutVar.first).c_str(),        10, 0, 10     ) );
        my_histos.emplace( "h_nMtops_"+cutVar.first,        std::make_shared<TH1D> ( ("h_nMtops_"+cutVar.first).c_str(),        ("h_nMtops_"+cutVar.first).c_str(),        10, 0, 10     ) );
        my_histos.emplace( "h_topsMass_"+cutVar.first,      std::make_shared<TH1D> ( ("h_topsMass_"+cutVar.first).c_str(),      ("h_topsMass_"+cutVar.first).c_str(),      1000, 0, 500  ) );
        my_histos.emplace( "h_topsPt_"+cutVar.first,        std::make_shared<TH1D> ( ("h_topsPt_"+cutVar.first).c_str(),        ("h_topsPt_"+cutVar.first).c_str(),        1000, 0, 2000 ) );
        my_histos.emplace( "h_6thJetPt_"+cutVar.first,      std::make_shared<TH1D> ( ("h_6thJetPt_"+cutVar.first).c_str(),      ("h_6thJetPt_"+cutVar.first).c_str(),      1000, 45, 150 ) );
        my_histos.emplace( "h_ht_"+cutVar.first,            std::make_shared<TH1D> ( ("h_ht_"+cutVar.first).c_str(),            ("h_ht_"+cutVar.first).c_str(),            60, 0, 3000   ) );
        my_histos.emplace( "h_dR_bjets_"+cutVar.first,      std::make_shared<TH1D> ( ("h_dR_bjets_"+cutVar.first).c_str(),      ("h_dR_bjets_"+cutVar.first).c_str(),      50, 0, 10     ) );

        //my_histos.emplace( "h_jetsMass_"+cutVar.first,      std::make_shared<TH1D> ( ("h_jetsMass_"+cutVar.first).c_str(),      ("h_jetsMass_"+cutVar.first).c_str(),      1000, 0, 500  ) );
        //my_histos.emplace( "h_jetsPt_"+cutVar.first,        std::make_shared<TH1D> ( ("h_jetsPt_"+cutVar.first).c_str(),        ("h_jetsPt_"+cutVar.first).c_str(),        1000, 0, 2000 ) );
        my_histos.emplace( "h_jetsPhi_"+cutVar.first,       std::make_shared<TH1D> ( ("h_jetsPhi_"+cutVar.first).c_str(),       ("h_jetsPhi_"+cutVar.first).c_str(),       80, -4, 4     ) );        
        my_histos.emplace( "h_jetsEta_"+cutVar.first,       std::make_shared<TH1D> ( ("h_jetsEta_"+cutVar.first).c_str(),       ("h_jetsEta_"+cutVar.first).c_str(),       100, -6, 6    ) );
        my_histos.emplace( "h_mbl_"+cutVar.first,           std::make_shared<TH1D> ( ("h_mbl_"+cutVar.first).c_str(),           ("h_mbl_"+cutVar.first).c_str(),           1000, 0, 300  ) );
        //my_histos.emplace( "h_bjetsMass_"+cutVar.first,     std::make_shared<TH1D> ( ("h_bjetsMass_"+cutVar.first).c_str(),     ("h_bjetsMass_"+cutVar.first).c_str(),     1000, 0, 500  ) );
        //my_histos.emplace( "h_bjetsEta_"+cutVar.first,      std::make_shared<TH1D> ( ("h_bjetsEta_"+cutVar.first).c_str(),      ("h_bjetsEta_"+cutVar.first).c_str(),      100, -6, 6    ) );
        //my_histos.emplace( "h_bjetsPhi_"+cutVar.first,      std::make_shared<TH1D> ( ("h_bjetsPhi_"+cutVar.first).c_str(),      ("h_bjetsPhi_"+cutVar.first).c_str(),      80, -4, 4     ) ); 
        //my_histos.emplace( "h_bjetsPt_"+cutVar.first,       std::make_shared<TH1D> ( ("h_bjetsPt_"+cutVar.first).c_str(),       ("h_bjetsPt_"+cutVar.first).c_str(),       1000, 0, 2000 ) );

        // for cut optimization of dR_bjets cut
        my_2d_histos.emplace( "h_njets_dR_bjets_"+cutVar.first, std::make_shared<TH2D>( ("h_njets_dR_bjets_"+cutVar.first).c_str(), ("h_njets_dR_bjets_"+cutVar.first).c_str(), 1000, 0, 10, 20, 0, 20 ) );   
        // 2d ntops plot  
        my_2d_histos.emplace( "h_nRtops_vs_nMtops_"+cutVar.first, std::make_shared<TH2D>( ("h_nRtops_vs_nMtops_"+cutVar.first).c_str(), ("h_nRtops_vs_nMtops_"+cutVar.first).c_str(), 12, -0.5, 11.5, 12, -0.5, 11.5 ) );
    }

}

// ---------------------------------------------
// -- Put everything you want to do per event 
// ---------------------------------------------
void Semra_Analyzer::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    while( tr.getNextEvent() )
    {
        
        const auto& eventCounter = tr.getVar<int>("eventCounter");
        
        // -------------------
        //  Print Event Number
        // -------------------
        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() & (10000 == 0) ) printf( " Event %i\n", tr.getEvtNum() );

        // General variables
        const auto& runtype               = tr.getVar<std::string>("runtype");     
        const auto& Jets                  = tr.getVec<utility::LorentzVector>("Jets");
        const auto& NGoodJets_pt30        = tr.getVar<int>("NGoodJets_pt30");
        const auto& NGoodBJets_pt30       = tr.getVar<int>("NGoodBJets_pt30");
        // Top variables
        const auto& ntops                 = tr.getVar<int>("ntops");
        const auto& ntops_1jet            = tr.getVar<int>("ntops_1jet"); // merged
        const auto& ntops_3jet            = tr.getVar<int>("ntops_3jet"); // resolved 
        const auto& topsMass              = tr.getVec<double>("topsMass");
        const auto& topsPt                = tr.getVec<double>("topsPt");
        // 0-lepton baseline selections
        const auto& passBaseline0l_Good   = tr.getVar<bool>("passBaseline0l_Good");
        const auto& GoodJets_pt45         = tr.getVec<bool>("GoodJets_pt45");
        const auto& HT_trigger_pt30       = tr.getVar<double>("HT_trigger_pt30");
        const auto& dR_bjets              = tr.getVar<double>("dR_bjets");
        //const bool ZeroNonIsoMuon         = NNonIsoMuons == 0;
        //const bool pass_ge7j_pt30         = NGoodJets_pt30 >= 7;
        //const bool pass_ge2t              = ntops >= 2;
        //const bool pass_ge2b_pt30         = NGoodBJets_pt30 >= 2;
        //const bool pass_ge1dRbjets        = dR_bjets >= 1.0;       
        // 1-lepton baseline selections
        const auto& passBaseline1l_Good   = tr.getVar<bool>("passBaseline1l_Good");
        const auto& GoodJets_pt30         = tr.getVec<bool>("GoodJets_pt30");
        const auto& Mbl                   = tr.getVar<double>("Mbl");

        // -------------------
        // -- Define weight
        // -------------------
        double eventweight = 1.0, pileupWeight = 1.0, prefiringScaleFactor = 1.0, bTagWeight = 1.0;
        double weight_0l = 1.0, leptonicWeight = 1.0;
        double weight_1l = 1.0, hadronicWeight = 1.0;
        
        if(runtype == "MC")
        {
            // general event weight
            const auto& Weight   = tr.getVar<float>("Weight"    );
            const auto& lumi     = tr.getVar<double>("FinalLumi");
            eventweight          = lumi * Weight;
            pileupWeight         = tr.getVar<double>("puWeightCorr"                    );
            prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor"            );
            bTagWeight           = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            // leptonic event weight            
            const auto& eleLepWeight = tr.getVar<double>("totGoodElectronSF");
            const auto& muLepWeight  = tr.getVar<double>("totGoodMuonSF"    );
            leptonicWeight           = eleLepWeight * muLepWeight;           
            // hadronic event weight
            hadronicWeight       = tr.getVar<double>("jetTrigSF"                       );
       
            weight_0l *= eventweight * pileupWeight * prefiringScaleFactor * bTagWeight * hadronicWeight;
            weight_1l *= eventweight * pileupWeight * prefiringScaleFactor * bTagWeight * leptonicWeight;
        }

        // -------------------------------
        // get the 6th jet pt for 0-lepton
        // -------------------------------
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
 
        // -------------------------------------------
        // Make cuts and fill histograms here & cutmap
        // -------------------------------------------
        const std::map<std::string, bool>& cutmap
        {
            // -------------------
            // Baseline selections
            // -------------------
            {"passBaseline0l_Good", passBaseline0l_Good},
            {"passBaseline1l_Good", passBaseline1l_Good},

        };

        if (!inithisto) 
        {
            InitHistos(cutmap);
            inithisto = true;
        }

        my_histos["EventCounter"]->Fill( eventCounter );

        // --------------------------
        // Fill the cutmap histograms
        // --------------------------     
        for (const auto& cutVar: cutmap) 
        {    
            if (cutVar.second) 
            {
                // --------------------------
                // get variables for 0-lepton
                // --------------------------
                if (cutVar.first == "passBaseline0l_Good")
                {
                    my_histos["h_njets_"+cutVar.first]->Fill( NGoodJets_pt30, weight_0l );
                    my_histos["h_nbjets_"+cutVar.first]->Fill( NGoodBJets_pt30, weight_0l );

                    // get top variables
                    my_histos["h_ntops_"+cutVar.first]->Fill( ntops, weight_0l );
                    my_histos["h_nRtops_"+cutVar.first]->Fill( ntops_3jet, weight_0l);
                    my_histos["h_nMtops_"+cutVar.first]->Fill( ntops_1jet, weight_0l);

                    for (unsigned int itops = 0; itops < topsMass.size(); itops++)
                    {
                        my_histos["h_topsMass_"+cutVar.first]->Fill( topsMass.at(itops), weight_0l );
                    }

                    for (unsigned int itops = 0; itops < topsPt.size(); itops++)
                    {
                        my_histos["h_topsPt_"+cutVar.first]->Fill( topsPt.at(itops), weight_0l );
                    }
               
                    // get SF variables 
                    my_histos["h_6thJetPt_"+cutVar.first]->Fill( SixthJetPt45, weight_0l );
                    my_histos["h_ht_"+cutVar.first]->Fill( HT_trigger_pt30, weight_0l );
                  
                    // get others
                    my_histos["h_dR_bjets_"+cutVar.first]->Fill( dR_bjets, weight_0l ); 
                
                    // 2d histograms
                    my_2d_histos["h_njets_dR_bjets_"+cutVar.first]->Fill( dR_bjets, NGoodJets_pt30, weight_0l );
                    my_2d_histos["h_njets_dR_bjets_"+cutVar.first]->GetXaxis()->SetTitle("#DeltaR_{bjets}");
                    my_2d_histos["h_njets_dR_bjets_"+cutVar.first]->GetYaxis()->SetTitle("N_{J}");
                    my_2d_histos["h_nRtops_vs_nMtops_"+cutVar.first]->Fill( ntops_1jet, ntops_3jet, weight_0l );
                    my_2d_histos["h_nRtops_vs_nMtops_"+cutVar.first]->GetXaxis()->SetTitle("N_{MergedTops}");
                    my_2d_histos["h_nRtops_vs_nMtops_"+cutVar.first]->GetYaxis()->SetTitle("N_{ResolvedTops}"); 

                }
        
                // --------------------------
                // get variables for 1-lepton
                // --------------------------
                else 
                {
                    my_histos["h_njets_"+cutVar.first]->Fill( NGoodJets_pt30, weight_1l );
                    my_histos["h_nbjets_"+cutVar.first]->Fill( NGoodBJets_pt30, weight_1l );
                    
                    // get SF variables
                    for (unsigned int j = 0; j < Jets.size(); j++)
                    {
                        if(!GoodJets_pt30[j]) continue;
                        
                        my_histos["h_jetsEta_"+cutVar.first]->Fill( Jets.at(j).Eta(), weight_1l );
                        my_histos["h_jetsPhi_"+cutVar.first]->Fill( Jets.at(j).Phi(), weight_1l );
                    }
 
                    // get others
                    my_histos["h_ht_"+cutVar.first]->Fill( HT_trigger_pt30, weight_1l );
                    my_histos["h_mbl_"+cutVar.first]->Fill( Mbl, weight_1l );

                }
            
            }
       } 
    }
}

void Semra_Analyzer::WriteHistos(TFile* outfile)
{
    outfile->cd();

    for (const auto &p : my_histos) {
        p.second->Write();
    }
 
    for (const auto &p : my_2d_histos) {
        p.second->Write();
    }
}
