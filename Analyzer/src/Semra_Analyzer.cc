#define Semra_Analyzer_cxx
#include "Analyzer/Analyzer/include/Semra_Analyzer.h"
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
#include <TDirectory.h>
#include <TH1F.h>

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
        // get jet and bjet variables
        my_histos.emplace( "h_njets_"+cutVar.first,         std::make_shared<TH1D> ( ("h_njets_"+cutVar.first).c_str(),         ("h_njets_"+cutVar.first).c_str(),         20, 0, 20     ) );
        my_histos.emplace( "h_jetsMass_"+cutVar.first,      std::make_shared<TH1D> ( ("h_jetsMass_"+cutVar.first).c_str(),      ("h_jetsMass_"+cutVar.first).c_str(),      1000, 0, 500  ) );
        my_histos.emplace( "h_jetsEta_"+cutVar.first,       std::make_shared<TH1D> ( ("h_jetsEta_"+cutVar.first).c_str(),       ("h_jetsEta_"+cutVar.first).c_str(),       100, -6, 6    ) );
        my_histos.emplace( "h_jetsPhi_"+cutVar.first,       std::make_shared<TH1D> ( ("h_jetsPhi_"+cutVar.first).c_str(),       ("h_jetsPhi_"+cutVar.first).c_str(),       80, -4, 4     ) );        
        my_histos.emplace( "h_jetsPt_"+cutVar.first,        std::make_shared<TH1D> ( ("h_jetsPt_"+cutVar.first).c_str(),        ("h_jetsPt_"+cutVar.first).c_str(),        1000, 0, 2000 ) );

        my_histos.emplace( "h_nbjets_"+cutVar.first,        std::make_shared<TH1D> ( ("h_nbjets_"+cutVar.first).c_str(),        ("h_nbjets_"+cutVar.first).c_str(),        20, 0, 20     ) );
        my_histos.emplace( "h_bjetsMass_"+cutVar.first,     std::make_shared<TH1D> ( ("h_bjetsMass_"+cutVar.first).c_str(),     ("h_bjetsMass_"+cutVar.first).c_str(),     1000, 0, 500  ) );
        my_histos.emplace( "h_bjetsEta_"+cutVar.first,      std::make_shared<TH1D> ( ("h_bjetsEta_"+cutVar.first).c_str(),      ("h_bjetsEta_"+cutVar.first).c_str(),      100, -6, 6    ) );
        my_histos.emplace( "h_bjetsPhi_"+cutVar.first,      std::make_shared<TH1D> ( ("h_bjetsPhi_"+cutVar.first).c_str(),      ("h_bjetsPhi_"+cutVar.first).c_str(),      80, -4, 4     ) ); 
        my_histos.emplace( "h_bjetsPt_"+cutVar.first,       std::make_shared<TH1D> ( ("h_bjetsPt_"+cutVar.first).c_str(),       ("h_bjetsPt_"+cutVar.first).c_str(),       1000, 0, 2000 ) );

        // get the top object (actual top jets from top utility::LorentzVector) jets' Mass, Eta, Phi, Pt  
        my_histos.emplace( "h_ntops_"+cutVar.first,         std::make_shared<TH1D> ( ("h_ntops_"+cutVar.first).c_str(),         ("h_ntops_"+cutVar.first).c_str(),         10, 0, 10     ) );
        my_histos.emplace( "h_nRtops_"+cutVar.first,        std::make_shared<TH1D> ( ("h_nRtops_"+cutVar.first).c_str(),        ("h_nRtops_"+cutVar.first).c_str(),        10, 0, 10     ) );
        my_histos.emplace( "h_nMtops_"+cutVar.first,        std::make_shared<TH1D> ( ("h_nMtops_"+cutVar.first).c_str(),        ("h_nMtops_"+cutVar.first).c_str(),        10, 0, 10     ) );

        my_histos.emplace( "h_topsMass_"+cutVar.first,      std::make_shared<TH1D> ( ("h_topsMass_"+cutVar.first).c_str(),      ("h_topsMass_"+cutVar.first).c_str(),      1000, 0, 500  ) );
        my_histos.emplace( "h_topsEta_"+cutVar.first,       std::make_shared<TH1D> ( ("h_topsEta_"+cutVar.first).c_str(),       ("h_topsEta_"+cutVar.first).c_str(),       100, -6, 6    ) );
        my_histos.emplace( "h_topsPhi_"+cutVar.first,       std::make_shared<TH1D> ( ("h_topsPhi_"+cutVar.first).c_str(),       ("h_topsPhi_"+cutVar.first).c_str(),       80, -4, 4     ) );
        my_histos.emplace( "h_topsPt_"+cutVar.first,        std::make_shared<TH1D> ( ("h_topsPt_"+cutVar.first).c_str(),        ("h_topsPt_"+cutVar.first).c_str(),        1000, 0, 2000 ) );

        my_histos.emplace( "h_bestTopMass_"+cutVar.first,   std::make_shared<TH1D> ( ("h_bestTopMass_"+cutVar.first).c_str(),   ("h_bestTopMass_"+cutVar.first).c_str(),   1000, 0, 500  ) );
        my_histos.emplace( "h_bestTopEta_"+cutVar.first,    std::make_shared<TH1D> ( ("h_bestTopEta_"+cutVar.first).c_str(),    ("h_bestTopEta_"+cutVar.first).c_str(),    100, -6, 6    ) );
        my_histos.emplace( "h_bestTopPhi_"+cutVar.first,    std::make_shared<TH1D> ( ("h_bestTopPhi_"+cutVar.first).c_str(),    ("h_bestTopPhi_"+cutVar.first).c_str(),    80, -4, 4     ) );
        my_histos.emplace( "h_bestTopPt_"+cutVar.first,     std::make_shared<TH1D> ( ("h_bestTopPt_"+cutVar.first).c_str(),     ("h_bestTopPt_"+cutVar.first).c_str(),     1000, 0, 2000 ) );
        
        // get other variables
        my_histos.emplace( "h_ht_"+cutVar.first,            std::make_shared<TH1D> ( ("h_ht_"+cutVar.first).c_str(),            ("h_ht_"+cutVar.first).c_str(),            60, 0, 3000   ) );
        my_histos.emplace( "h_met_"+cutVar.first,           std::make_shared<TH1D> ( ("h_met_"+cutVar.first).c_str(),           ("h_met_"+cutVar.first).c_str(),           200, 0, 2000  ) );
        my_histos.emplace( "h_dR_bjets_"+cutVar.first,      std::make_shared<TH1D> ( ("h_dR_bjets_"+cutVar.first).c_str(),      ("h_dR_bjets_"+cutVar.first).c_str(),      50, 0, 10     ) );
        my_histos.emplace( "h_dR_top1_top2_"+cutVar.first,  std::make_shared<TH1D> ( ("h_dR_top1_top2_"+cutVar.first).c_str(),  ("h_dR_top1_top2_"+cutVar.first).c_str(),  50, 0, 10     ) );
        my_histos.emplace( "h_dR_tops_bjets_"+cutVar.first, std::make_shared<TH1D> ( ("h_dR_tops_bjets_"+cutVar.first).c_str(), ("h_dR_tops_bjets_"+cutVar.first).c_str(), 50, 0, 10     ) );

        // get variables for QCD CR
        my_histos.emplace( "h_DoubleDisCo_disc1_"+cutVar.first,   std::make_shared<TH1D> ( ("h_DoubleDisCo_disc1_"+cutVar.first).c_str(),   ("h_DoubleDisCo_disc1_"+cutVar.first).c_str(),   100, 0, 1    ) );
        my_histos.emplace( "h_DoubleDisCo_disc2_"+cutVar.first,   std::make_shared<TH1D> ( ("h_DoubleDisCo_disc2_"+cutVar.first).c_str(),   ("h_DoubleDisCo_disc2_"+cutVar.first).c_str(),   100, 0, 1    ) );
        my_histos.emplace( "h_DoubleDisCo_massReg_"+cutVar.first, std::make_shared<TH1D> ( ("h_DoubleDisCo_massReg_"+cutVar.first).c_str(), ("h_DoubleDisCo_massReg_"+cutVar.first).c_str(), 150, 0, 1500 ) );

        // for cut optimization of dR_bjets cut
        my_2d_histos.emplace( "h_njets_dR_bjets_"+cutVar.first, std::make_shared<TH2D>( ("h_njets_dR_bjets_"+cutVar.first).c_str(), ("h_njets_dR_bjets_"+cutVar.first).c_str(), 1000, 0, 10, 20, 0, 20 ) );                    
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
        
        //-------------------------
        // -- Print Event Number 
        //-------------------------
        
        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() & (10000 == 0) ) printf( " Event %i\n", tr.getEvtNum() );

        // General variables
        const auto& runtype               = tr.getVar<std::string>("runtype");     
        const auto& Jets                  = tr.getVec<utility::LorentzVector>("Jets");
        const auto& MET                   = tr.getVar<float>("MET");
        const auto& HT_trigger_pt30       = tr.getVar<float>("HT_trigger_pt30");
        const auto& GoodJets_pt30         = tr.getVec<bool>("GoodJets_pt30");
        const auto& NGoodJets_pt30        = tr.getVar<int>("NGoodJets_pt30");
        const auto& GoodBJets_pt30        = tr.getVec<bool>("GoodBJets_pt30");
        const auto& NGoodBJets_pt30       = tr.getVar<int>("NGoodBJets_pt30"); 
        // Top variables
        const auto& ntops                 = tr.getVar<int>("ntops");
        const auto& ntops_1jet            = tr.getVar<int>("ntops_1jet"); // merged
        const auto& ntops_2jet            = tr.getVar<int>("ntops_2jet"); // medium
        const auto& ntops_3jet            = tr.getVar<int>("ntops_3jet"); // resolved 
        const auto& topsMass              = tr.getVec<float>("topsMass");
        const auto& topsEta               = tr.getVec<float>("topsEta");
        const auto& topsPhi               = tr.getVec<float>("topsPhi");  
        const auto& topsPt                = tr.getVec<float>("topsPt");
        const auto& topsLV                = tr.getVec<utility::LorentzVector>("topsLV");
        const auto& bestTopMass           = tr.getVar<float>("bestTopMass");
        const auto& bestTopEta            = tr.getVar<float>("bestTopEta");
        const auto& bestTopPhi            = tr.getVar<float>("bestTopPhi");
        const auto& bestTopPt             = tr.getVar<float>("bestTopPt");
        const auto& dR_bjets              = tr.getVar<float>("dR_bjets");
        const auto& dR_top1_top2          = tr.getVar<float>("dR_top1_top2");
        // Baseline selection
        const auto& passBaseline0l_pre    = tr.getVar<bool>("passBaseline0l_pre");
        const auto& passBaseline0l_good   = tr.getVar<bool>("passBaseline0l_good");
        const auto& NNonIsoMuons          = tr.getVar<int>("NNonIsoMuons");
        const bool ZeroNonIsoMuon         = NNonIsoMuons == 0;
        const bool pass_HT500_pt30        = HT_trigger_pt30 > 500;
        const bool pass_ge7j_pt30         = NGoodJets_pt30 >= 7;
        const bool pass_ge2b_pt30         = NGoodBJets_pt30 >= 2;
        const bool pass_ge1t              = ntops >= 1;
        const bool pass_ge1tR             = ntops >= 1 && ntops_1jet == 0 && ntops_2jet == 0;
        const bool pass_ge1tM             = ntops >= 1 && ntops_3jet == 0 && ntops_2jet == 0;
        const bool pass_ge2t              = ntops >= 2;
        const bool pass_ge2tR             = ntops >= 2 && ntops_1jet == 0 && ntops_2jet == 0;
        const bool pass_ge2tM             = ntops >= 2 && ntops_3jet == 0 && ntops_2jet == 0;
        const bool pass_ge2tRM            = ntops >= 2 && ntops_3jet >= 1 && ntops_1jet >= 1 && ntops_2jet == 0;
        const bool pass_ge1dRbjets        = dR_bjets >= 1.0;       
        // QCD CR
        const bool pass_qcdCR                            = tr.getVar<bool>("pass_qcdCR"); // 1l qcd cr
        const auto NNonIsoMuonJets_pt30                  = tr.getVar<int>("NNonIsoMuonJets_pt30"); 
        const auto DoubleDisCo_disc1_NonIsoMuon_0l_RPV   = tr.getVar<float>("DoubleDisCo_disc1_NonIsoMuon_0l_RPV");
        const auto DoubleDisCo_disc2_NonIsoMuon_0l_RPV   = tr.getVar<float>("DoubleDisCo_disc2_NonIsoMuon_0l_RPV");
        const auto DoubleDisCo_massReg_NonIsoMuon_0l_RPV = tr.getVar<float>("DoubleDisCo_massReg_NonIsoMuon_0l_RPV");
        const auto DoubleDisCo_disc1_0l_RPV              = tr.getVar<float>("DoubleDisCo_disc1_0l_RPV");
        const auto DoubleDisCo_disc2_0l_RPV              = tr.getVar<float>("DoubleDisCo_disc2_0l_RPV");
        const auto DoubleDisCo_massReg_0l_RPV            = tr.getVar<float>("DoubleDisCo_massReg_0l_RPV");
        const bool pass_7j_pt30           = NGoodJets_pt30 == 7;
        const bool pass_8j_pt30           = NGoodJets_pt30 == 8;
        const bool pass_9j_pt30           = NGoodJets_pt30 == 9;
        const bool pass_10j_pt30          = NGoodJets_pt30 == 10;
        const bool pass_ge11j_pt30        = NGoodJets_pt30 >= 11;

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
            const auto& Weight   = tr.getVar<float>("Weight");
            eventweight          = lumi*Weight;
            puScaleFactor        = tr.getVar<double>("puWeightCorr");
        
            weight *= eventweight * puScaleFactor;
        }

        // ---------------------------------------------
        // -- Calculate DeltaR between tops and bjets
        // ---------------------------------------------
        std::vector<utility::LorentzVector> bjets;
        
        for(unsigned int ijet = 0; ijet < Jets.size(); ijet++)
        {
            if(!GoodBJets_pt30[ijet]) continue;
            bjets.push_back(Jets.at(ijet));        
        }
        
        std::vector<double> dR_top_bjet;
        
        for (unsigned int t = 0; t < topsLV.size(); t++) 
        {
            for (unsigned int b = 0; b < bjets.size(); b++) 
            {
                double deltaR = utility::DeltaR(topsLV.at(t), bjets.at(b));
                dR_top_bjet.push_back(deltaR);
            }
        }
 

        // -------------------------------------------------
        // -- Make cuts and fill histograms here & cutmap
        // -------------------------------------------------
        const std::map<std::string, bool>& cutmap
        {
            // -------------------
            // Baseline selections
            // -------------------
            // without NonIsoMuon cut
            {"0l_HT500_ge7j_ge2t_ge1dRbjets",  passBaseline0l_pre && pass_ge7j_pt30  && pass_ge2t && pass_ge1dRbjets },
            {"0l_HT500_7j_ge2t_ge1dRbjets",    passBaseline0l_pre && pass_7j_pt30    && pass_ge2t && pass_ge1dRbjets },
            {"0l_HT500_8j_ge2t_ge1dRbjets",    passBaseline0l_pre && pass_8j_pt30    && pass_ge2t && pass_ge1dRbjets },
            {"0l_HT500_9j_ge2t_ge1dRbjets",    passBaseline0l_pre && pass_9j_pt30    && pass_ge2t && pass_ge1dRbjets },
            {"0l_HT500_10j_ge2t_ge1dRbjets",   passBaseline0l_pre && pass_10j_pt30   && pass_ge2t && pass_ge1dRbjets },
            {"0l_HT500_ge11j_ge2t_ge1dRbjets", passBaseline0l_pre && pass_ge11j_pt30 && pass_ge2t && pass_ge1dRbjets },

            // with NonIsoMuon cut
            {"0l_HT500_0NonIsoMuon_ge7j_ge2t_ge1dRbjets",  passBaseline0l_pre && passBaseline0l_good },
            {"0l_HT500_0NonIsoMuon_7j_ge2t_ge1dRbjets",    passBaseline0l_pre && passBaseline0l_good },
            {"0l_HT500_0NonIsoMuon_8j_ge2t_ge1dRbjets",    passBaseline0l_pre && passBaseline0l_good },
            {"0l_HT500_0NonIsoMuon_9j_ge2t_ge1dRbjets",    passBaseline0l_pre && passBaseline0l_good },
            {"0l_HT500_0NonIsoMuon_10j_ge2t_ge1dRbjets",   passBaseline0l_pre && passBaseline0l_good },
            {"0l_HT500_0NonIsoMuon_ge11j_ge2t_ge1dRbjets", passBaseline0l_pre && passBaseline0l_good },
            
            // -----------------
            // QCD CR selections
            // -----------------
            {"qcdCR_0l_HT300_1NonIsoMuon_ge7NonIsoMuonJet",  pass_qcdCR                    }, 
            {"qcdCR_0l_HT300_1NonIsoMuon_7NonIsoMuonJet",    pass_qcdCR && pass_7j_pt30    },
            {"qcdCR_0l_HT300_1NonIsoMuon_8NonIsoMuonJet",    pass_qcdCR && pass_8j_pt30    },
            {"qcdCR_0l_HT300_1NonIsoMuon_9NonIsoMuonJet",    pass_qcdCR && pass_9j_pt30    },
            {"qcdCR_0l_HT300_1NonIsoMuon_10NonIsoMuonJet",   pass_qcdCR && pass_10j_pt30   },
            {"qcdCR_0l_HT300_1NonIsoMuon_ge11NonIsoMuonJet", pass_qcdCR && pass_ge11j_pt30 },

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
                //my_histos["h_njets_"+cutVar.first]->Fill( NGoodJets_pt30, weight );
                my_histos["h_nbjets_"+cutVar.first]->Fill( NGoodBJets_pt30, weight );
                my_histos["h_ntops_"+cutVar.first]->Fill( ntops, weight );
                my_histos["h_nRtops_"+cutVar.first]->Fill( ntops_3jet, weight);
                my_histos["h_nMtops_"+cutVar.first]->Fill( ntops_1jet, weight);
 
                // -----------------------------
                // -- jets & bjets mass & pT
                // -----------------------------
                for(unsigned int ijet = 0; ijet < Jets.size(); ijet++) 
                {
                    if(!GoodJets_pt30[ijet]) continue;
                    my_histos["h_jetsMass_"+cutVar.first]->Fill(Jets.at(ijet).M(), weight);
                    my_histos["h_jetsEta_"+cutVar.first]->Fill(Jets.at(ijet).Eta(), weight);
                    my_histos["h_jetsPhi_"+cutVar.first]->Fill(Jets.at(ijet).Phi(), weight);
                    my_histos["h_jetsPt_"+cutVar.first]->Fill(Jets.at(ijet).Pt(), weight);
 
                    if(!GoodBJets_pt30[ijet]) continue;
                    const utility::LorentzVector& bjet = Jets.at(ijet);                     
                    my_histos["h_bjetsMass_"+cutVar.first]->Fill(bjet.M(), weight);
                    my_histos["h_bjetsEta_"+cutVar.first]->Fill(bjet.Eta(), weight);
                    my_histos["h_bjetsPhi_"+cutVar.first]->Fill(bjet.Phi(), weight);
                    my_histos["h_bjetsPt_"+cutVar.first]->Fill(bjet.Pt(), weight);
                }
        
                // --------------------------------------
                // -- get top jets' Mass, Eta, Phi, Pt 
                // --------------------------------------
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

                my_histos["h_bestTopMass_"+cutVar.first]->Fill( bestTopMass, weight );
                my_histos["h_bestTopEta_"+cutVar.first]->Fill( bestTopEta, weight );
                my_histos["h_bestTopPhi_"+cutVar.first]->Fill( bestTopPhi, weight );
                my_histos["h_bestTopPt_"+cutVar.first]->Fill( bestTopPt, weight );

                my_histos["h_ht_"+cutVar.first]->Fill( HT_trigger_pt30, weight );
                my_histos["h_met_"+cutVar.first]->Fill( MET, weight );
                my_histos["h_dR_bjets_"+cutVar.first]->Fill( dR_bjets, weight );
                my_histos["h_dR_top1_top2_"+cutVar.first]->Fill( dR_top1_top2, weight );
    
                // ---------------------------------
                // -- deltaR between top and bjet
                // ---------------------------------
                for (unsigned int idR = 0; idR < dR_top_bjet.size(); idR++) 
                {
                    my_histos["h_dR_tops_bjets_"+cutVar.first]->Fill( dR_top_bjet.at(idR), weight );        
                }
         
                my_2d_histos["h_njets_dR_bjets_"+cutVar.first]->Fill( dR_bjets, NGoodJets_pt30, weight );
                my_2d_histos["h_njets_dR_bjets_"+cutVar.first]->GetXaxis()->SetTitle("#DeltaR_{bjets}");
                my_2d_histos["h_njets_dR_bjets_"+cutVar.first]->GetYaxis()->SetTitle("N_{J}");

                // ------------------------------------------
                // -- get the njets with QCD CR selections 
                // ------------------------------------------
                if ( cutVar.first.find("qcdCR") != std::string::npos &&  cutVar.first.find("1NonIsoMuon") != std::string::npos )
                {
                    my_histos["h_njets_"+cutVar.first]->Fill( NNonIsoMuonJets_pt30, weight );   
                    my_histos["h_DoubleDisCo_disc1_"+cutVar.first]->Fill( DoubleDisCo_disc1_NonIsoMuon_0l_RPV, weight );
                    my_histos["h_DoubleDisCo_disc2_"+cutVar.first]->Fill( DoubleDisCo_disc2_NonIsoMuon_0l_RPV, weight );
                    my_histos["h_DoubleDisCo_massReg_"+cutVar.first]->Fill( DoubleDisCo_massReg_NonIsoMuon_0l_RPV, weight );          
                }

                else
                {
                    my_histos["h_njets_"+cutVar.first]->Fill( NGoodJets_pt30, weight );
                    my_histos["h_DoubleDisCo_disc1_"+cutVar.first]->Fill( DoubleDisCo_disc1_0l_RPV, weight );
                    my_histos["h_DoubleDisCo_disc2_"+cutVar.first]->Fill( DoubleDisCo_disc2_0l_RPV, weight );
                    my_histos["h_DoubleDisCo_massReg_"+cutVar.first]->Fill( DoubleDisCo_massReg_0l_RPV, weight );
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

    for (const auto &p : my_efficiencies) {
        p.second->Write();
    }    
}
