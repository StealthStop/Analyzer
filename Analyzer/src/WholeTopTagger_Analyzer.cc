#define WholeTopTagger_Analyzer_cxx
#include "Analyzer/Analyzer/include/WholeTopTagger_Analyzer.h"
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
#include <TDirectory.h>
#include <TH1F.h>

WholeTopTagger_Analyzer::WholeTopTagger_Analyzer() : inithisto(false) 
{
}

// -------------------
// -- Define histos
// -------------------
void WholeTopTagger_Analyzer::InitHistos(const std::map<std::string, bool>& cutmap) 
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
   
    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ; 
    
    for (const auto& cutVar : cutmap) 
    {  
        // -------------
        // 1D Histograms
        // -------------
        // Top variables
        my_histos.emplace( "h_ntops_"+cutVar.first,    std::make_shared<TH1D> ( ("h_ntops_"+cutVar.first).c_str(),    ("h_ntops_"+cutVar.first).c_str(),    10, 0, 10     ) );
        my_histos.emplace( "h_nRtops_"+cutVar.first,   std::make_shared<TH1D> ( ("h_nRtops_"+cutVar.first).c_str(),   ("h_nRtops_"+cutVar.first).c_str(),   10, 0, 10     ) );
        my_histos.emplace( "h_nMtops_"+cutVar.first,   std::make_shared<TH1D> ( ("h_nMtops_"+cutVar.first).c_str(),   ("h_nMtops_"+cutVar.first).c_str(),   10, 0, 10     ) );
        
        my_histos.emplace( "h_topsMass_"+cutVar.first, std::make_shared<TH1D> ( ("h_topsMass_"+cutVar.first).c_str(), ("h_topsMass_"+cutVar.first).c_str(), 1000, 0, 500  ) );
        my_histos.emplace( "h_topsEta_"+cutVar.first,  std::make_shared<TH1D> ( ("h_topsEta_"+cutVar.first).c_str(),  ("h_topsEta_"+cutVar.first).c_str(),  100, -6, 6    ) );
        my_histos.emplace( "h_topsPhi_"+cutVar.first,  std::make_shared<TH1D> ( ("h_topsPhi_"+cutVar.first).c_str(),  ("h_topsPhi_"+cutVar.first).c_str(),  80, -4, 4     ) );
        my_histos.emplace( "h_topsPt_"+cutVar.first,   std::make_shared<TH1D> ( ("h_topsPt_"+cutVar.first).c_str(),   ("h_topsPt_"+cutVar.first).c_str(),   1000, 0, 2000 ) );

        // Gen Match Top variables - Efficiency
        my_histos.emplace( "h_genMatch_TopMass_R_"+cutVar.first, std::make_shared<TH1D> ( ("h_genMatch_TopMass_R_"+cutVar.first).c_str(), ("h_genMatch_TopMass_R_"+cutVar.first).c_str(), 1000, 0, 500  ) );
        my_histos.emplace( "h_genMatch_TopEta_R_"+cutVar.first,  std::make_shared<TH1D> ( ("h_genMatch_TopEta_R_"+cutVar.first).c_str(),  ("h_genMatch_TopEta_R_"+cutVar.first).c_str(),  100, -6, 6    ) );
        my_histos.emplace( "h_genMatch_TopPhi_R_"+cutVar.first,  std::make_shared<TH1D> ( ("h_genMatch_TopPhi_R_"+cutVar.first).c_str(),  ("h_genMatch_TopPhi_R_"+cutVar.first).c_str(),  80, -4, 4     ) );
        my_histos.emplace( "h_genMatch_TopPt_R_"+cutVar.first,   std::make_shared<TH1D> ( ("h_genMatch_TopPt_R_"+cutVar.first).c_str(),   ("h_genMatch_TopPt_R_"+cutVar.first).c_str(),   1000, 0, 2000 ) );

        my_histos.emplace( "h_genMatch_TopMass_M_"+cutVar.first, std::make_shared<TH1D> ( ("h_genMatch_TopMass_M_"+cutVar.first).c_str(), ("h_genMatch_TopMass_M_"+cutVar.first).c_str(), 1000, 0, 500  ) );
        my_histos.emplace( "h_genMatch_TopEta_M_"+cutVar.first,  std::make_shared<TH1D> ( ("h_genMatch_TopEta_M_"+cutVar.first).c_str(),  ("h_genMatch_TopEta_M_"+cutVar.first).c_str(),  100, -6, 6    ) );
        my_histos.emplace( "h_genMatch_TopPhi_M_"+cutVar.first,  std::make_shared<TH1D> ( ("h_genMatch_TopPhi_M_"+cutVar.first).c_str(),  ("h_genMatch_TopPhi_M_"+cutVar.first).c_str(),  80, -4, 4     ) );
        my_histos.emplace( "h_genMatch_TopPt_M_"+cutVar.first,   std::make_shared<TH1D> ( ("h_genMatch_TopPt_M_"+cutVar.first).c_str(),   ("h_genMatch_TopPt_M_"+cutVar.first).c_str(),   1000, 0, 2000 ) );

        my_histos.emplace( "h_genMatch_TopMass_C_"+cutVar.first, std::make_shared<TH1D> ( ("h_genMatch_TopMass_C_"+cutVar.first).c_str(), ("h_genMatch_TopMass_C_"+cutVar.first).c_str(), 1000, 0, 500  ) );
        my_histos.emplace( "h_genMatch_TopEta_C_"+cutVar.first,  std::make_shared<TH1D> ( ("h_genMatch_TopEta_C_"+cutVar.first).c_str(),  ("h_genMatch_TopEta_C-"+cutVar.first).c_str(),  100, -6, 6    ) );
        my_histos.emplace( "h_genMatch_TopPhi_C_"+cutVar.first,  std::make_shared<TH1D> ( ("h_genMatch_TopPhi_C_"+cutVar.first).c_str(),  ("h_genMatch_TopPhi_C_"+cutVar.first).c_str(),  80, -4, 4     ) );
        my_histos.emplace( "h_genMatch_TopPt_C_"+cutVar.first,   std::make_shared<TH1D> ( ("h_genMatch_TopPt_C_"+cutVar.first).c_str(),   ("h_genMatch_TopPt_C_"+cutVar.first).c_str(),   1000, 0, 2000 ) );

        // Not Gen Match Top variables - Fake Rate
        my_histos.emplace( "h_fakeRate_njets_R_"+cutVar.first, std::make_shared<TH1D> ( ("h_fakeRate_njets_R_"+cutVar.first).c_str(), ("h_fakeRate_njets_R_"+cutVar.first).c_str(), 21, -0.5, 20.5) );
        my_histos.emplace( "h_fakeRate_njets_M_"+cutVar.first, std::make_shared<TH1D> ( ("h_fakeRate_njets_M_"+cutVar.first).c_str(), ("h_fakeRate_njets_M_"+cutVar.first).c_str(), 21, -0.5, 20.5) );
        my_histos.emplace( "h_fakeRate_njets_C_"+cutVar.first, std::make_shared<TH1D> ( ("h_fakeRate_njets_C_"+cutVar.first).c_str(), ("h_fakeRate_njets_C_"+cutVar.first).c_str(), 21, -0.5, 20.5) );    

        my_histos.emplace( "h_fakeRate_nbjets_R_"+cutVar.first, std::make_shared<TH1D> ( ("h_fakeRate_nbjets_R_"+cutVar.first).c_str(), ("h_fakeRate_nbjets_R_"+cutVar.first).c_str(), 21, -0.5, 20.5) );
        my_histos.emplace( "h_fakeRate_nbjets_M_"+cutVar.first, std::make_shared<TH1D> ( ("h_fakeRate_nbjets_M_"+cutVar.first).c_str(), ("h_fakeRate_nbjets_M_"+cutVar.first).c_str(), 21, -0.5, 20.5) );
        my_histos.emplace( "h_fakeRate_nbjets_C_"+cutVar.first, std::make_shared<TH1D> ( ("h_fakeRate_nbjets_C_"+cutVar.first).c_str(), ("h_fakeRate_nbjets_C_"+cutVar.first).c_str(), 21, -0.5, 20.5) );
        
        // Gen Match & Not Gen Match Top Discriminators - ROC
        my_histos.emplace( "h_genMatch_TopDisc_"+cutVar.first, std::make_shared<TH1D>    ( ("h_genMatch_TopDisc_"+cutVar.first).c_str(),    ("h_genMatch_TopDisc_"+cutVar.first).c_str(),    1000, 0, 1) );
        my_histos.emplace( "h_notGenMatch_TopDisc_"+cutVar.first, std::make_shared<TH1D> ( ("h_notGenMatch_TopDisc_"+cutVar.first).c_str(), ("h_notGenMatch_TopDisc_"+cutVar.first).c_str(), 1000, 0, 1) );

        // -------------
        // 2D Histograms
        // -------------
        my_2d_histos.emplace( "h_nRtops_vs_nMtops_"+cutVar.first, std::make_shared<TH2D>( ("h_nRtops_vs_nMtops_"+cutVar.first).c_str(), ("h_nRtops_vs_nMtops_"+cutVar.first).c_str(), 12, -0.5, 11.5, 12, -0.5, 11.5 ) );
    }

}

// ---------------------------------------------
// -- Put everything you want to do per event 
// ---------------------------------------------
void WholeTopTagger_Analyzer::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    while( tr.getNextEvent() )
    {
        
        const auto& eventCounter = tr.getVar<int>("eventCounter");
        
        //-------------------------
        // -- Print Event Number 
        //-------------------------
        
        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() & (10000 == 0) ) printf( " Event %i\n", tr.getEvtNum() );

        // General variables
        const auto& runtype            = tr.getVar<std::string>("runtype");     
        const auto& Jets               = tr.getVec<TLorentzVector>("Jets");
        // Top variables
        const auto& ntops              = tr.getVar<int>("ntops");
        const auto& ntops_3jet         = tr.getVar<int>("ntops_3jet"); // resolved 
        const auto& ntops_1jet         = tr.getVar<int>("ntops_1jet"); // merged
        const auto& topsMass           = tr.getVec<double>("topsMass");
        const auto& topsEta            = tr.getVec<double>("topsEta");
        const auto& topsPhi            = tr.getVec<double>("topsPhi");  
        const auto& topsPt             = tr.getVec<double>("topsPt");
        const auto& genTopMatchMass_R  = tr.getVar<double>("genTopMatchMass_R");
        const auto& genTopMatchEta_R   = tr.getVar<double>("genTopMatchEta_R");
        const auto& genTopMatchPhi_R   = tr.getVar<double>("genTopMatchPhi_R");
        const auto& genTopMatchPt_R    = tr.getVar<double>("genTopMatchPt_R");
        const auto& genTopMatchMass_M  = tr.getVar<double>("genTopMatchMass_M");
        const auto& genTopMatchEta_M   = tr.getVar<double>("genTopMatchEta_M");
        const auto& genTopMatchPhi_M   = tr.getVar<double>("genTopMatchPhi_M");
        const auto& genTopMatchPt_M    = tr.getVar<double>("genTopMatchPt_M");
        const auto& genTopMatchMass_C  = tr.getVar<double>("genTopMatchMass_C");
        const auto& genTopMatchEta_C   = tr.getVar<double>("genTopMatchEta_C");
        const auto& genTopMatchPhi_C   = tr.getVar<double>("genTopMatchPhi_C");
        const auto& genTopMatchPt_C    = tr.getVar<double>("genTopMatchPt_C");
        const auto& isFakeRate_R       = tr.getVar<bool>("isFakeRate_R");
        const auto& isFakeRate_M       = tr.getVar<bool>("isFakeRate_M");
        const auto& isFakeRate_C       = tr.getVar<bool>("isFakeRate_C");
        const auto& topDiscGenMatch    = tr.getVar<double>("topDiscGenMatch");
        const auto& topDiscNotGenMatch = tr.getVar<double>("topDiscNotGenMatch");
        // New Baseline selection
        const bool passBaseline0l_pre  = tr.getVar<bool>("passBaseline0l_pre");
        const auto& NNonIsoMuons       = tr.getVar<int>("NNonIsoMuons");
        const auto& NGoodJets_pt30     = tr.getVar<int>("NGoodJets_pt30");
        const auto& NGoodBJets_pt30    = tr.getVar<int>("NGoodBJets_pt30");
        const auto& dR_bjets           = tr.getVar<double>("dR_bjets");
        const bool pass_newTopTagger   = passBaseline0l_pre
                                        && NNonIsoMuons == 0   && NGoodBJets_pt30 >= 2
                                        && NGoodJets_pt30 >= 7 && dR_bjets >= 1.0;
        const bool pass_7j_pt30        = NGoodJets_pt30 == 7;
        const bool pass_8j_pt30        = NGoodJets_pt30 == 8;
        const bool pass_9j_pt30        = NGoodJets_pt30 == 9;
        const bool pass_10j_pt30       = NGoodJets_pt30 == 10;
        const bool pass_11j_pt30       = NGoodJets_pt30 == 11;
        const bool pass_ge12j_pt30     = NGoodJets_pt30 >= 12;

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
            const auto& lumi     = tr.getVar<double>("Lumi");
            const auto& Weight   = tr.getVar<double>("Weight");
            eventweight          = lumi*Weight;
            puScaleFactor        = tr.getVar<double>("puWeightCorr");
        
            weight *= eventweight * puScaleFactor;
        }

        // -------------------------------------------------
        // -- Make cuts and fill histograms here & cutmap
        // -------------------------------------------------
        const std::map<std::string, bool>& cutmap
        {
            // New baseline selection
            {"0l_HT500_0NonIsoMuon_ge7j_ge1dRbjets",  pass_newTopTagger                   },
            {"0l_HT500_0NonIsoMuon_7j_ge1dRbjets",    pass_newTopTagger && pass_7j_pt30   },
            {"0l_HT500_0NonIsoMuon_8j_ge1dRbjets",    pass_newTopTagger && pass_8j_pt30   },
            {"0l_HT500_0NonIsoMuon_9j_ge1dRbjets",    pass_newTopTagger && pass_9j_pt30   },
            {"0l_HT500_0NonIsoMuon_10j_ge1dRbjets",   pass_newTopTagger && pass_10j_pt30  },
            {"0l_HT500_0NonIsoMuon_11j_ge1dRbjets",   pass_newTopTagger && pass_11j_pt30  },
            {"0l_HT500_0NonIsoMuon_ge12j_ge1dRbjets", pass_newTopTagger && pass_ge12j_pt30},           
 
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
                // -------------
                // 1D Histograms
                // ------------- 
                // Top variables
                my_histos["h_ntops_"+cutVar.first]->Fill( ntops, weight );
                my_histos["h_nRtops_"+cutVar.first]->Fill( ntops_3jet, weight);
                my_histos["h_nMtops_"+cutVar.first]->Fill( ntops_1jet, weight);
                
                for (unsigned int itops = 0; itops < topsMass.size(); itops++) 
                {
                    my_histos["h_topsMass_"+cutVar.first]->Fill( topsMass.at(itops), weight );
                }
        
                for (unsigned int itops = 0; itops < topsEta.size(); itops++) 
                {
                    my_histos["h_topsEta_"+cutVar.first]->Fill( topsEta.at(itops), weight );
                }
       
                for (unsigned int itops = 0; itops < topsPhi.size(); itops++) 
                {
                    my_histos["h_topsPhi_"+cutVar.first]->Fill( topsPhi.at(itops), weight );
                }
 
                for (unsigned int itops = 0; itops < topsPt.size(); itops++) 
                {
                    my_histos["h_topsPt_"+cutVar.first]->Fill( topsPt.at(itops), weight );
                }      

                // Gen Match Top variables - Efficiency
                my_histos["h_genMatch_TopMass_R_"+cutVar.first]->Fill( genTopMatchMass_R, weight );
                my_histos["h_genMatch_TopEta_R_"+cutVar.first]->Fill( genTopMatchEta_R, weight );
                my_histos["h_genMatch_TopPhi_R_"+cutVar.first]->Fill( genTopMatchPhi_R, weight );
                my_histos["h_genMatch_TopPt_R_"+cutVar.first]->Fill( genTopMatchPt_R, weight );

                my_histos["h_genMatch_TopMass_M_"+cutVar.first]->Fill( genTopMatchMass_M, weight );
                my_histos["h_genMatch_TopEta_M_"+cutVar.first]->Fill( genTopMatchEta_M, weight );
                my_histos["h_genMatch_TopPhi_M_"+cutVar.first]->Fill( genTopMatchPhi_M, weight );
                my_histos["h_genMatch_TopPt_M_"+cutVar.first]->Fill( genTopMatchPt_M, weight );

                my_histos["h_genMatch_TopMass_C_"+cutVar.first]->Fill( genTopMatchMass_C, weight );
                my_histos["h_genMatch_TopEta_C_"+cutVar.first]->Fill( genTopMatchEta_C, weight );
                my_histos["h_genMatch_TopPhi_C_"+cutVar.first]->Fill( genTopMatchPhi_C, weight );
                my_histos["h_genMatch_TopPt_C_"+cutVar.first]->Fill( genTopMatchPt_C, weight );

                // Not Gen Match Top variables - Fake Rate
                if (isFakeRate_R)
                {
                    my_histos["h_fakeRate_njets_R_"+cutVar.first]->Fill( NGoodJets_pt30, weight );
                    my_histos["h_fakeRate_nbjets_R_"+cutVar.first]->Fill( NGoodBJets_pt30, weight );
                }

                if (isFakeRate_M)
                {
                    my_histos["h_fakeRate_njets_M_"+cutVar.first]->Fill( NGoodJets_pt30, weight );
                    my_histos["h_fakeRate_nbjets_M_"+cutVar.first]->Fill( NGoodBJets_pt30, weight ); 
                }

                if (isFakeRate_C)
                {
                    my_histos["h_fakeRate_njets_C_"+cutVar.first]->Fill( NGoodJets_pt30, weight );
                    my_histos["h_fakeRate_nbjets_C_"+cutVar.first]->Fill( NGoodBJets_pt30, weight );
                }            

                // Gen Match & Not Gen Match Top Discriminators - ROC
                my_histos["h_genMatch_TopDisc_"+cutVar.first]->Fill( topDiscGenMatch, weight );
                my_histos["h_notGenMatch_TopDisc_"+cutVar.first]->Fill( topDiscNotGenMatch, weight );

                // -------------
                // 2D Histograms
                // ------------- 
                my_2d_histos["h_nRtops_vs_nMtops_"+cutVar.first]->Fill(ntops_1jet, ntops_3jet, weight);
                my_2d_histos["h_nRtops_vs_nMtops_"+cutVar.first]->GetXaxis()->SetTitle("N_{MergedTops}");
                my_2d_histos["h_nRtops_vs_nMtops_"+cutVar.first]->GetYaxis()->SetTitle("N_{ResolvedTops}");
            
            }
            
        }
    } 
}

void WholeTopTagger_Analyzer::WriteHistos(TFile* outfile)
{
    outfile->cd();

    for (const auto &p : my_histos) {
        p.second->Write();
    }
 
    for (const auto &p : my_2d_histos) {
        p.second->Write();
    }

    for (const auto &p : my_efficiencies) {
        p.second->Write();
    }    
}
