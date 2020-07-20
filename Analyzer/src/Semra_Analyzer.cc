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
        // --------------------------------------
        // -- variables for baseline selection
        // --------------------------------------
        my_histos.emplace( "h_ntops_"+cutVar.first, std::make_shared<TH1D> ( ("h_ntops_"+cutVar.first).c_str(), ("h_ntops_"+cutVar.first).c_str(), 10, 0, 10 ) );
        my_histos.emplace( "h_njets_"+cutVar.first, std::make_shared<TH1D> ( ("h_njets_"+cutVar.first).c_str(), ("h_njets_"+cutVar.first).c_str(), 20, 0, 20 ) );
        my_histos.emplace( "h_nbjets_"+cutVar.first, std::make_shared<TH1D> ( ("h_nbjets_"+cutVar.first).c_str(), ("h_nbjets_"+cutVar.first).c_str(), 20, 0, 20 ) );
        my_histos.emplace( "h_ht_"+cutVar.first, std::make_shared<TH1D> ( ("h_ht_"+cutVar.first).c_str(), ("h_ht_"+cutVar.first).c_str(), 60, 0, 3000 ) );
        my_histos.emplace( "h_met_"+cutVar.first, std::make_shared<TH1D> ( ("h_met_"+cutVar.first).c_str(), ("h_met_"+cutVar.first).c_str(), 200, 0, 2000 ) ) ;
        //my_histos.emplace( "h_jetsPt_"+cutVar.first, std::make_shared<TH1D> ( ("h_jetsPt_"+cutVar.first).c_str(), ("h_jetsPt_"+cutVar.first).c_str(), 1000, 0, 2000 ) );
        //my_histos.emplace( "h_jetsMass_"+cutVar.first, std::make_shared<TH1D> ( ("h_jetsMass_"+cutVar.first).c_str(), ("h_jetsMass_"+cutVar.first).c_str(), 1000, 0, 500) );
        //my_histos.emplace( "h_bjetsPt_"+cutVar.first, std::make_shared<TH1D> ( ("h_bjetsPt_"+cutVar.first).c_str(), ("h_bjetsPt_"+cutVar.first).c_str(), 1000, 0, 2000 ) );
        //my_histos.emplace( "h_bjetsMass_"+cutVar.first, std::make_shared<TH1D> ( ("h_bjetsMass_"+cutVar.first).c_str(), ("h_bjetsMass_"+cutVar.first).c_str(), 1000, 0, 500 ) );
        
        // get the top object (actual top jets from top TLorentzVector) jets' Mass, Eta, Phi, Pt  
        my_histos.emplace( "h_topsMass_"+cutVar.first, std::make_shared<TH1D> ( ("h_topsMass_"+cutVar.first).c_str(), ("h_topsMass_"+cutVar.first).c_str(), 1000, 0, 500 ) );
        my_histos.emplace( "h_topsEta_"+cutVar.first, std::make_shared<TH1D> ( ("h_topsEta_"+cutVar.first).c_str(), ("h_topsEta_"+cutVar.first).c_str(), 100, -6, 6 ) );
        my_histos.emplace( "h_topsPhi_"+cutVar.first, std::make_shared<TH1D> ( ("h_topsPhi_"+cutVar.first).c_str(), ("h_topsPhi_"+cutVar.first).c_str(), 80, -4, 4 ) );
        my_histos.emplace( "h_topsPt_"+cutVar.first, std::make_shared<TH1D> ( ("h_topsPt_"+cutVar.first).c_str(), ("h_topsPt_"+cutVar.first).c_str(), 1000, 0, 2000 ) );

        // get the each resolved jets' Mass, Eta, Phi, Pt 
        my_histos.emplace( "h_resolvedMass_"+cutVar.first, std::make_shared<TH1D> ( ("h_resolvedMass_"+cutVar.first).c_str(), ("h_resolvedMass_"+cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_resolvedEta_"+cutVar.first, std::make_shared<TH1D> ( ("h_resolvedEta_"+cutVar.first).c_str(), ("h-resolvedEta_"+cutVar.first).c_str(), 100, -6, 6 ) );
        my_histos.emplace( "h_resolvedPhi_"+cutVar.first, std::make_shared<TH1D> ( ("h_resolvedPhi_"+cutVar.first).c_str(), ("h_resolvedPhi_"+cutVar.first).c_str(), 80, -4, 4 ) );
        my_histos.emplace( "h_resolvedPt_"+cutVar.first, std::make_shared<TH1D> ( ("h_resolvedPt_"+cutVar.first).c_str(), ("h_resolvedPt_"+cutVar.first).c_str(), 100, 0, 1000 ) );

        my_histos.emplace( "h_bestTopMass_"+cutVar.first, std::make_shared<TH1D> ( ("h_bestTopMass_"+cutVar.first).c_str(), ("h_bestTopMass_"+cutVar.first).c_str(), 1000, 0, 500 ) );
        my_histos.emplace( "h_bestTopEta_"+cutVar.first, std::make_shared<TH1D> ( ("h_bestTopEta_"+cutVar.first).c_str(), ("h_bestTopEta_"+cutVar.first).c_str(), 100, -6, 6 ) );
        my_histos.emplace( "h_bestTopPt_"+cutVar.first, std::make_shared<TH1D> ( ("h_bestTopPt_"+cutVar.first).c_str(), ("h_bestTopPt_"+cutVar.first).c_str(), 1000, 0, 2000 ) );
        
        my_histos.emplace( "h_dR_bjets_"+cutVar.first, std::make_shared<TH1D> ( ("h_dR_bjets_"+cutVar.first).c_str(), ("h_dR_bjets_"+cutVar.first).c_str(), 50, 0, 10 ) );
        my_histos.emplace( "h_dR_top1_top2_"+cutVar.first, std::make_shared<TH1D> ( ("h_dR_top1_top2_"+cutVar.first).c_str(), ("h_dR_top1_top2_"+cutVar.first).c_str(), 50, 0, 10 ) );
        my_histos.emplace( "h_dR_tops_bjets_"+cutVar.first, std::make_shared<TH1D> ( ("h_dR_tops_bjets_"+cutVar.first).c_str(), ("h_dR_tops_bjets_"+cutVar.first).c_str(), 50, 0, 10 ) );

        //my_2d_histos.emplace( "h_njets_MVA_"+cutVar.first, std::make_shared<TH2D>( ("h_njets_MVA_"+cutVar.first).c_str(), ("h_njets_MVA_"+cutVar.first).c_str(), 8, 7, 15, 50, 0, 1.0 ) );
        //my_2d_histos.emplace( "h_njets_dR_bjets_"+cutVar.first, std::make_shared<TH2D>( ("h_njets_dR_bjets_"+cutVar.first).c_str(), ("h_njets_dR_bjets_"+cutVar.first).c_str(), 1000, 0, 10, 20, 0, 20 ) ); // for cut optimization of dR_bjets cut                   

        // --------------------------
        // -- stop MT2 hemispheres  
        // --------------------------
        // pt rank 
        //my_histos.emplace( "h_stop1Mass_PtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1Mass_PtRank_"+cutVar.first).c_str(), ("h_stop1Mass_PtRank_"+cutVar.first).c_str(), 500, 0, 1500) );
        //my_histos.emplace( "h_stop1Eta_PtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1Eta_PtRank_"+cutVar.first).c_str(), ("h_stop1Eta_PtRank_"+cutVar.first).c_str(), 100, -6, 6 ) );
        //my_histos.emplace( "h_stop1Phi_PtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1Phi_PtRank_"+cutVar.first).c_str(), ("h_stop1Phi_PtRank_"+cutVar.first).c_str(), 80, -4, 4 ) );
        //my_histos.emplace( "h_stop1Pt_PtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1Pt_PtRank_"+cutVar.first).c_str(), ("h_stop1Pt_PtRank_"+cutVar.first).c_str(), 100, 0, 1000 ) ); 
        //my_histos.emplace( "h_stop2Mass_PtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2Mass_PtRank_"+cutVar.first).c_str(), ("h_stop2Mass_PtRank_"+cutVar.first).c_str(), 500, 0, 1500) );
        //my_histos.emplace( "h_stop2Eta_PtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2Eta_PtRank_"+cutVar.first).c_str(), ("h_stop2Eta_PtRank_"+cutVar.first).c_str(), 100, -6, 6 ) );
        //my_histos.emplace( "h_stop2Phi_PtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2Phi_PtRank_"+cutVar.first).c_str(), ("h_stop2Phi_PtRank_"+cutVar.first).c_str(), 80, -4, 4 ) );
        //my_histos.emplace( "h_stop2Pt_PtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2Pt_PtRank_"+cutVar.first).c_str(), ("h_stop2Pt_PtRank_"+cutVar.first).c_str(), 100, 0, 1000 ) );
        // mass rank
        //my_histos.emplace( "h_stop1Mass_MassRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1Mass_MassRank_"+cutVar.first).c_str(), ("h_stop1Mass_MassRank_"+cutVar.first).c_str(), 500, 0, 1500) );
        //my_histos.emplace( "h_stop1Eta_MassRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1Eta_MassRank_"+cutVar.first).c_str(), ("h_stop1Eta_MassRank_"+cutVar.first).c_str(), 100, -6, 6 ) );
        //my_histos.emplace( "h_stop1Phi_MassRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1Phi_MassRank_"+cutVar.first).c_str(), ("h_stop1Phi_MassRank_"+cutVar.first).c_str(), 80, -4, 4 ) );
        //my_histos.emplace( "h_stop1Pt_MassRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1Pt_MassRank_"+cutVar.first).c_str(), ("h_stop1Pt_MassRank_"+cutVar.first).c_str(), 100, 0, 1000 ) );
        //my_histos.emplace( "h_stop2Mass_MassRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2Mass_MassRank_"+cutVar.first).c_str(), ("h_stop2Mass_MassRank_"+cutVar.first).c_str(), 500, 0, 1500) );
        //my_histos.emplace( "h_stop2Eta_MassRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2Eta_MassRank_"+cutVar.first).c_str(), ("h_stop2Eta_MassRank_"+cutVar.first).c_str(), 100, -6, 6 ) );
        //my_histos.emplace( "h_stop2Phi_MassRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2Phi_MassRank_"+cutVar.first).c_str(), ("h_stop2Phi_MassRank_"+cutVar.first).c_str(), 80, -4, 4 ) );
        //my_histos.emplace( "h_stop2Pt_MassRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2Pt_MassRank_"+cutVar.first).c_str(), ("h_stop2Pt_MassRank_"+cutVar.first).c_str(), 100, 0, 1000 ) );
        // scalarPt rank
        //my_histos.emplace( "h_stop1Mass_ScalarPtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1Mass_ScalarPtRank_"+cutVar.first).c_str(), ("h_stop1Mass_ScalarPtRank_"+cutVar.first).c_str(), 500, 0, 1500) );
        //my_histos.emplace( "h_stop1Eta_ScalarPtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1Eta_ScalarPtRank_"+cutVar.first).c_str(), ("h_stop1Eta_ScalarPtRank_"+cutVar.first).c_str(), 100, -6, 6 ) );
        //my_histos.emplace( "h_stop1Phi_ScalarPtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1Phi_ScalarPtRank_"+cutVar.first).c_str(), ("h_stop1Phi_ScalarPtRank_"+cutVar.first).c_str(), 80, -4, 4 ) );
        //my_histos.emplace( "h_stop1Pt_ScalarPtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1Pt_ScalarPtRank_"+cutVar.first).c_str(), ("h_stop1Pt_ScalarPtRank_"+cutVar.first).c_str(), 100, 0, 1000 ) );
        //my_histos.emplace( "h_stop2Mass_ScalarPtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2Mass_ScalarPtRank_"+cutVar.first).c_str(), ("h_stop2Mass_ScalarPtRank_"+cutVar.first).c_str(), 500, 0, 1500) );
        //my_histos.emplace( "h_stop2Eta_ScalarPtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2Eta_ScalarPtRank_"+cutVar.first).c_str(), ("h_stop2Eta_ScalarPtRank_"+cutVar.first).c_str(), 100, -6, 6 ) );
        //my_histos.emplace( "h_stop2Phi_ScalarPtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2Phi_ScalarPtRank_"+cutVar.first).c_str(), ("h_stop2Phi_ScalarPtRank_"+cutVar.first).c_str(), 80, -4, 4 ) );
        //my_histos.emplace( "h_stop2Pt_ScalarPtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2Pt_ScalarPtRank_"+cutVar.first).c_str(), ("h_stop2Pt_ScalarPtRank_"+cutVar.first).c_str(), 100, 0, 1000 ) );       
        
        //my_histos.emplace( "h_stop1ScalarPt_ScalarPtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1ScalarPt_ScalarPtRank_"+cutVar.first).c_str(), ("h_stop1ScalarPt_ScalarPtRank_"+cutVar.first).c_str(), 100, 0, 1000 ) );
        //my_histos.emplace( "h_stop2ScalarPt_ScalarPtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2ScalarPt_ScalarPtRank_"+cutVar.first).c_str(), ("h_stop2ScalarPt_ScalarPtRank_"+cutVar.first).c_str(), 100, 0, 1000 ) );
        
        //my_histos.emplace( "h_MT2_"+cutVar.first, std::make_shared<TH1D> ( ("h_MT2_"+cutVar.first).c_str(), ("h_MT2_"+cutVar.first).c_str(), 500, 0, 1500) );
        //my_histos.emplace( "h_dR_stop1stop2_"+cutVar.first, std::make_shared<TH1D> ( ("h_dR_stop1stop2_"+cutVar.first).c_str(), ("h_dR_stop1stop2_"+cutVar.first).c_str(), 50, 0, 10 ) );
        //my_histos.emplace( "h_dPhi_stop1stop2_"+cutVar.first, std::make_shared<TH1D> ( ("h_dPhi_stop1stop2_"+cutVar.first).c_str(), ("h_dPhi_stop1stop2_"+cutVar.first).c_str(), 50, 0, 10 ) );
        //my_histos.emplace( "h_difference_stopMasses_"+cutVar.first, std::make_shared<TH1D> ( ("h_difference_stopMasses_"+cutVar.first).c_str(), ("h_difference_stopMasses_"+cutVar.first).c_str(), 500, 0, 1500) );
        //my_histos.emplace( "h_average_stopMasses_"+cutVar.first, std::make_shared<TH1D> ( ("h_average_stopMasses_"+cutVar.first).c_str(), ("h_average_stopMasses_"+cutVar.first).c_str(), 500, 0, 1500) );
        //my_histos.emplace( "h_relativeDiff_stopMasses_"+cutVar.first, std::make_shared<TH1D> ( ("h_relativeDiff_stopMasses_"+cutVar.first).c_str(), ("h_relativeDiff_stopMasses_"+cutVar.first).c_str(), 500, -1500, 1500) );

        // stop1VSstop2 Mass, Eta, Phi, Pt       
        // pt rank
        //my_2d_histos.emplace( "h_Mass_stop1vsstop2_PtRank_"+cutVar.first, std::make_shared<TH2D>( ("h_Mass_stop1vsstop2_PtRank_"+cutVar.first).c_str(), ("h_Mass_stop1vsstop2_PtRank_"+cutVar.first).c_str(), 500, 0, 1500, 500, 0, 1500 ) );
        //my_2d_histos.emplace( "h_Eta_stop1vsstop2_PtRank_"+cutVar.first, std::make_shared<TH2D>( ("h_Eta_stop1vsstop2_PtRank_"+cutVar.first).c_str(), ("h_Eta_stop1vsstop2_PtRank_"+cutVar.first).c_str(), 100, -6, 6, 100, -6, 6 ) );
        //my_2d_histos.emplace( "h_Phi_stop1vsstop2_PtRank_"+cutVar.first, std::make_shared<TH2D>( ("h_Phi_stop1vsstop2_PtRank_"+cutVar.first).c_str(), ("h_Phi_stop1vsstop2_PtRank_"+cutVar.first).c_str(), 80, -4, 4, 80, -4, 4 ) );
        //my_2d_histos.emplace( "h_Pt_stop1vsstop2_PtRank_"+cutVar.first, std::make_shared<TH2D>( ("h_Pt_stop1vsstop2_PtRank_"+cutVar.first).c_str(), ("h_Pt_stop1vsstop2_PtRank_"+cutVar.first).c_str(), 100, 0, 1000, 100, 0, 1000) );
        // mass rank
        //my_2d_histos.emplace( "h_Mass_stop1vsstop2_MassRank_"+cutVar.first, std::make_shared<TH2D>( ("h_Mass_stop1vsstop2_MassRank_"+cutVar.first).c_str(), ("h_Mass_stop1vsstop2_MassRank_"+cutVar.first).c_str(), 500, 0, 1500, 500, 0, 1500 ) );
        //my_2d_histos.emplace( "h_Eta_stop1vsstop2_MassRank_"+cutVar.first, std::make_shared<TH2D>( ("h_Eta_stop1vsstop2_MassRank_"+cutVar.first).c_str(), ("h_Eta_stop1vsstop2_MassRank_"+cutVar.first).c_str(), 100, -6, 6, 100, -6, 6 ) );
        //my_2d_histos.emplace( "h_Phi_stop1vsstop2_MassRank_"+cutVar.first, std::make_shared<TH2D>( ("h_Phi_stop1vsstop2_MassRank_"+cutVar.first).c_str(), ("h_Phi_stop1vsstop2_MassRank_"+cutVar.first).c_str(), 80, -4, 4, 80, -4, 4 ) );
        //my_2d_histos.emplace( "h_Pt_stop1vsstop2_MassRank_"+cutVar.first, std::make_shared<TH2D>( ("h_Pt_stop1vsstop2_MassRank_"+cutVar.first).c_str(), ("h_Pt_stop1vsstop2_MassRank_"+cutVar.first).c_str(), 100, 0, 1000, 100, 0, 1000) );
        // scalarPt rank
        //my_2d_histos.emplace( "h_Mass_stop1vsstop2_ScalarPtRank_"+cutVar.first, std::make_shared<TH2D>( ("h_Mass_stop1vsstop2_ScalarPtRank_"+cutVar.first).c_str(), ("h_Mass_stop1vsstop2_ScalarPtRank_"+cutVar.first).c_str(), 500, 0, 1500, 500, 0, 1500 ) );    
        //my_2d_histos.emplace( "h_Eta_stop1vsstop2_ScalarPtRank_"+cutVar.first, std::make_shared<TH2D>( ("h_Eta_stop1vsstop2_ScalarPtRank_"+cutVar.first).c_str(), ("h_Eta_stop1vsstop2_ScalarPtRank_"+cutVar.first).c_str(), 100, -6, 6, 100, -6, 6 ) );
        //my_2d_histos.emplace( "h_Phi_stop1vsstop2_ScalarPtRank_"+cutVar.first, std::make_shared<TH2D>( ("h_Phi_stop1vsstop2_ScalarPtRank_"+cutVar.first).c_str(), ("h_Phi_stop1vsstop2_ScalarPtRank_"+cutVar.first).c_str(), 80, -4, 4, 80, -4, 4 ) );
        //my_2d_histos.emplace( "h_Pt_stop1vsstop2_ScalarPtRank_"+cutVar.first, std::make_shared<TH2D>( ("h_Pt_stop1vsstop2_ScalarPtRank_"+cutVar.first).c_str(), ("h_Pt_stop1vsstop2_ScalarPtRank_"+cutVar.first).c_str(), 100, 0, 1000, 100, 0, 1000) ); 
        
        // NJetsVSstops
        //my_2d_histos.emplace( "h_Mass_NJetsVSstop1_PtRank_"+cutVar.first, std::make_shared<TH2D> ( ("h_Mass_NJetsVSstop1_PtRank_"+cutVar.first).c_str(), ("h_Mass_NJetsVSstop1_PtRank_"+cutVar.first).c_str(), 20, 0, 20, 500, 0, 1500 ) );
        //my_2d_histos.emplace( "h_Mass_NJetsVSstop2_PtRank_"+cutVar.first, std::make_shared<TH2D> ( ("h_Mass_NJetsVSstop2_PtRank_"+cutVar.first).c_str(), ("h_Mass_NJetsVSstop2_PtRank_"+cutVar.first).c_str(), 20, 0, 20, 500, 0, 1500 ) );
        //my_2d_histos.emplace( "h_Mass_NJetsVSstop1_MassRank_"+cutVar.first, std::make_shared<TH2D> ( ("h_Mass_NJetsVSstop1_MassRank_"+cutVar.first).c_str(), ("h_Mass_NJetsVSstop1_MassRank_"+cutVar.first).c_str(), 20, 0, 20, 500, 0, 1500 ) );
        //my_2d_histos.emplace( "h_Mass_NJetsVSstop2_MassRank_"+cutVar.first, std::make_shared<TH2D> ( ("h_Mass_NJetsVSstop2_MassRank_"+cutVar.first).c_str(), ("h_Mass_NJetsVSstop2_MassRank_"+cutVar.first).c_str(), 20, 0, 20, 500, 0, 1500 ) );       
        //my_2d_histos.emplace( "h_Mass_NJetsVSstop1_ScalarPtRank_"+cutVar.first, std::make_shared<TH2D> ( ("h_Mass_NJetsVSstop1_ScalarPtRank_"+cutVar.first).c_str(), ("h_Mass_NJetsVSstop1_ScalarPtRank_"+cutVar.first).c_str(), 20, 0, 20, 500, 0, 1500 ) );
        //my_2d_histos.emplace( "h_Mass_NJetsVSstop2_ScalarPtRank_"+cutVar.first, std::make_shared<TH2D> ( ("h_Mass_NJetsVSstop2_ScalarPtRank_"+cutVar.first).c_str(), ("h_Mass_NJetsVSstop2_ScalarPtRank_"+cutVar.first).c_str(), 20, 0, 20, 500, 0, 1500 ) );
        //my_2d_histos.emplace( "h_stopMasses_diffVSavg_"+cutVar.first, std::make_shared<TH2D>( ("h_stopMasses_diffVSavg_"+cutVar.first).c_str(), ( "h_stopMasses_diffVSavg"+cutVar.first).c_str(), 150, -1500, 1500, 150, 0, 1500 ) );

        // stop TaggedTop Pt combinations 
        //my_2d_histos.emplace( "h_ScalarPt_stop1vsstop2_ScalarPtRank_"+cutVar.first, std::make_shared<TH2D>( ("h_ScalarPt_stop1vsstop2_ScalarPtRank_"+cutVar.first).c_str(), ("h_ScalarPt_stop1vsstop2_ScalarPtRank_"+cutVar.first).c_str(), 100, 0, 1000, 100, 0, 1000) );
        //my_2d_histos.emplace( "h_Pt1_PtRankVsScalarPtRank_"+cutVar.first, std::make_shared<TH2D>( ("h_Pt1_PtRankVsScalarPtRank_"+cutVar.first).c_str(), ("h_Pt1_PtRankVsScalarPtRank_"+cutVar.first).c_str(), 100, 0, 1000, 100, 0, 1000) );
        //my_2d_histos.emplace( "h_Pt2_PtRankVsScalarPtRank_"+cutVar.first, std::make_shared<TH2D>( ("h_Pt2_PtRankVsScalarPtRank_"+cutVar.first).c_str(), ("h_Pt2_PtRankVsScalarPtRank_"+cutVar.first).c_str(), 100, 0, 1000, 100, 0, 1000) );
        //my_2d_histos.emplace( "h_Pt1_PtRankVsScalarPt1_ScalarPtRank_"+cutVar.first, std::make_shared<TH2D>( ("h_Pt1_PtRankVsScalarPt1_ScalarPtRank_"+cutVar.first).c_str(), ("h_Pt1_PtRankVsScalarPt1_ScalarPtRank_"+cutVar.first).c_str(), 100, 0, 1000, 100, 0, 1000) );
        //my_2d_histos.emplace( "h_Pt2_PtRankVsScalarPt2_ScalarPtRank"+cutVar.first, std::make_shared<TH2D>( ("h_Pt2_PtRankVsScalarPt2_ScalarPtRank_"+cutVar.first).c_str(), ("h_Pt2_PtRankVsScalarPt2_ScalarPtRank_"+cutVar.first).c_str(), 100, 0, 1000, 100, 0, 1000) );

        // stops MassVsPt
        // pt rank
        //my_2d_histos.emplace( "h_stop1_MassVsPt_PtRank_"+cutVar.first, std::make_shared<TH2D>( ("h_stop1_MassVsPt_PtRank_"+cutVar.first).c_str(), ("h_stop1_MassVsPt_PtRank_"+cutVar.first).c_str(), 500, 0, 1500, 100, 0, 1000 ) );
        //my_2d_histos.emplace( "h_stop2_MassVsPt_PtRank_"+cutVar.first, std::make_shared<TH2D>( ("h_stop2_MassVsPt_PtRank_"+cutVar.first).c_str(), ("h_stop2_MassVsPt_PtRank_"+cutVar.first).c_str(), 500, 0, 1500, 100, 0, 1000 ) );
        // mask rank
        //my_2d_histos.emplace( "h_stop1_MassVsPt_MassRank_"+cutVar.first, std::make_shared<TH2D>( ("h_stop1_MassVsPt_MassRank_"+cutVar.first).c_str(), ("h_stop1_MassVsPt_MassRank_"+cutVar.first).c_str(), 500, 0, 1500, 100, 0, 1000 ) );
        //my_2d_histos.emplace( "h_stop2_MassVsPt_MassRank_"+cutVar.first, std::make_shared<TH2D>( ("h_stop2_MassVsPt_MassRank_"+cutVar.first).c_str(), ("h_stop2_MassVsPt_MassRank_"+cutVar.first).c_str(), 500, 0, 1500, 100, 0, 1000 ) );
        // scalarPt rank
        //my_2d_histos.emplace( "h_stop1_MassVsPt_ScalarPtRank_"+cutVar.first, std::make_shared<TH2D>( ("h_stop1_MassVsPt_ScalarPtRank_"+cutVar.first).c_str(), ("h_stop1_MassVsPt_ScalarPtRank_"+cutVar.first).c_str(), 500, 0, 1500, 100, 0, 1000 ) );
        //my_2d_histos.emplace( "h_stop2_MassVsPt_ScalarPtRank_"+cutVar.first, std::make_shared<TH2D>( ("h_stop2_MassVsPt_ScalarPtRank_"+cutVar.first).c_str(), ("h_stop2_MassVsPt_ScalarPtRank_"+cutVar.first).c_str(), 500, 0, 1500, 100, 0, 1000 ) );
        //my_2d_histos.emplace( "h_stop1_MassVsScalarPt_ScalarPtRank_"+cutVar.first, std::make_shared<TH2D>( ("h_stop1_MassVsScalarPt_ScalarPtRank_"+cutVar.first).c_str(), ("h_stop1_MassVsScalarPt_ScalarPtRank_"+cutVar.first).c_str(), 500, 0, 1500, 100, 0, 1000 ) );
        //my_2d_histos.emplace( "h_stop2_MassVsScalarPt_ScalarPtRank_"+cutVar.first, std::make_shared<TH2D>( ("h_stop2_MassVsScalarPt_ScalarPtRank_"+cutVar.first).c_str(), ("h_stop2_MassVsScalarPt_ScalarPtRank_"+cutVar.first).c_str(), 500, 0, 1500, 100, 0, 1000 ) );

        // MT2
        //my_2d_histos.emplace( "h_NJetsVsMT2_"+cutVar.first, std::make_shared<TH2D> ( ("h_NJetsVsMT2_"+cutVar.first).c_str(), ("h_NJetsVsMT2_"+cutVar.first).c_str(), 20, 0, 20, 500, 0, 1500 ) );
        //my_2d_histos.emplace( "h_Mass_MT2vsstop1_PtRank_"+cutVar.first, std::make_shared<TH2D> ( ("h_Mass_MT2vsstop1_PtRank_"+cutVar.first).c_str(), ("h_Mass_MT2vsstop1_PtRank_"+cutVar.first).c_str(), 500, 0, 1500, 500, 0, 1500 ) );
        //my_2d_histos.emplace( "h_Mass_MT2vsstop2_PtRank_"+cutVar.first, std::make_shared<TH2D> ( ("h_Mass_MT2vsstop2_PtRank_"+cutVar.first).c_str(), ("h_Mass_MT2vsstop2_PtRank_"+cutVar.first).c_str(), 500, 0, 1500, 500, 0, 1500 ) );
        //my_2d_histos.emplace( "h_Mass_MT2vsstop1_MassRank_"+cutVar.first, std::make_shared<TH2D> ( ("h_Mass_MT2vsstop1_MassRank_"+cutVar.first).c_str(), ("h_Mass_MT2vsstop1_MassRank_"+cutVar.first).c_str(), 500, 0, 1500, 500, 0, 1500 ) );
        //my_2d_histos.emplace( "h_Mass_MT2vsstop2_MassRank_"+cutVar.first, std::make_shared<TH2D> ( ("h_Mass_MT2vsstop2_MassRank_"+cutVar.first).c_str(), ("h_Mass_MT2vsstop2_MassRank_"+cutVar.first).c_str(), 500, 0, 1500, 500, 0, 1500 ) ); 
        //my_2d_histos.emplace( "h_Mass_MT2vsstop1_ScalarPtRank_"+cutVar.first, std::make_shared<TH2D> ( ("h_Mass_MT2vsstop1_ScalarPtRank_"+cutVar.first).c_str(), ("h_Mass_MT2vsstop1_ScalarPtRank_"+cutVar.first).c_str(), 500, 0, 1500, 500, 0, 1500 ) );
        //my_2d_histos.emplace( "h_Mass_MT2vsstop2_ScalarPtRank_"+cutVar.first, std::make_shared<TH2D> ( ("h_Mass_MT2vsstop2_ScalarPtRank_"+cutVar.first).c_str(), ("h_Mass_MT2vsstop2_ScalarPtRank_"+cutVar.first).c_str(), 500, 0, 1500, 500, 0, 1500 ) );
 
        // --------------------
        // -- stop gen level
        // --------------------
        // Reco
        //my_histos.emplace( "h_StopRecoMT2_"+cutVar.first, std::make_shared<TH1D> ( ("h_StopRecoMT2_"+cutVar.first).c_str(), ("h_StopRecoMT2_"+cutVar.first).c_str(), 500, 0, 1500) ); 
        //my_histos.emplace( "h_Stop1RecoMass_"+cutVar.first, std::make_shared<TH1D> ( ("h_Stop1RecoMass_"+cutVar.first).c_str(), ("h_Stop1RecoMass_"+cutVar.first).c_str(), 500, 0, 1500) );
        //my_histos.emplace( "h_Stop2RecoMass_"+cutVar.first, std::make_shared<TH1D> ( ("h_Stop2RecoMass_"+cutVar.first).c_str(), ("h_Stop2RecoMass_"+cutVar.first).c_str(), 500, 0, 1500) );
        //my_histos.emplace( "h_Nlino1RecoMass_"+cutVar.first, std::make_shared<TH1D> ( ("h_Nlino1RecoMass_"+cutVar.first).c_str(), ("h_Nlino1RecoMass_"+cutVar.first).c_str(), 500, 0, 1500) );
        //my_histos.emplace( "h_Nlino2RecoMass_"+cutVar.first, std::make_shared<TH1D> ( ("h_Nlino2RecoMass_"+cutVar.first).c_str(), ("h_Nlino2RecoMass_"+cutVar.first).c_str(), 500, 0, 1500) );
        //my_histos.emplace( "h_Single1RecoMass_"+cutVar.first, std::make_shared<TH1D> ( ("h_Single1RecoMass_"+cutVar.first).c_str(), ("h_Single1RecoMass_"+cutVar.first).c_str(), 500, 0, 1500) );
        //my_histos.emplace( "h_Single2RecoMass_"+cutVar.first, std::make_shared<TH1D> ( ("h_Single2RecoMass_"+cutVar.first).c_str(), ("h_Single2RecoMass_"+cutVar.first).c_str(), 500, 0, 1500) );
        // Gen
        //my_histos.emplace( "h_StopGenMT2_"+cutVar.first, std::make_shared<TH1D> ( ("h_StopGenMT2_"+cutVar.first).c_str(), ("h_StopGenMT2_"+cutVar.first).c_str(), 500, 0, 1500) );
        //my_histos.emplace( "h_Stop1GenMass_"+cutVar.first, std::make_shared<TH1D> ( ("h_Stop1GenMass_"+cutVar.first).c_str(), ("h_Stop1GenMass_"+cutVar.first).c_str(), 500, 0, 1500) );
        //my_histos.emplace( "h_Stop2GenMass_"+cutVar.first, std::make_shared<TH1D> ( ("h_Stop2GenMass_"+cutVar.first).c_str(), ("h_Stop2GenMass_"+cutVar.first).c_str(), 500, 0, 1500) );
        //my_histos.emplace( "h_Nlino1GenMass_"+cutVar.first, std::make_shared<TH1D> ( ("h_Nlino1GenMass_"+cutVar.first).c_str(), ("h_Nlino1GenMass_"+cutVar.first).c_str(), 500, 0, 1500) );
        //my_histos.emplace( "h_Nlino2GenMass_"+cutVar.first, std::make_shared<TH1D> ( ("h_Nlino2GenMass_"+cutVar.first).c_str(), ("h_Nlino2GenMass_"+cutVar.first).c_str(), 500, 0, 1500) );
        //my_histos.emplace( "h_Single1GenMass_"+cutVar.first, std::make_shared<TH1D> ( ("h_Single1GenMass_"+cutVar.first).c_str(), ("h_Single1GenMass_"+cutVar.first).c_str(), 500, 0, 1500) );
        //my_histos.emplace( "h_Single2GenMass_"+cutVar.first, std::make_shared<TH1D> ( ("h_Single2GenMass_"+cutVar.first).c_str(), ("h_Single2GenMass_"+cutVar.first).c_str(), 500, 0, 1500) );

    }

    // cut flow absolute numbers 
    //my_histos.emplace( "h_cutFlow_absolute", std::make_shared<TH1D>("h_cutFlow_absolute", "h_cutFlow_absolute", 9,0,9));    
}

// ---------------------------------------------
// -- Put everything you want to do per event 
// ---------------------------------------------
void Semra_Analyzer::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    while( tr.getNextEvent() )
    {
        const auto& eventCounter    = tr.getVar<int>("eventCounter");
        
        //-------------------------
        // -- Print Event Number 
        //-------------------------
        
        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() & (10000 == 0) ) printf( " Event %i\n", tr.getEvtNum() );

        const auto& runtype         = tr.getVar<std::string>("runtype");     
        const auto& JetID           = tr.getVar<bool>("JetID");
        const auto& Jets            = tr.getVec<TLorentzVector>("Jets");
        const auto& NGoodLeptons    = tr.getVar<int>("NGoodLeptons");
        const auto& MET             = tr.getVar<double>("MET");
        const auto& GoodJets_pt45   = tr.getVec<bool>("GoodJets_pt45");
        const auto& GoodBJets_pt45  = tr.getVec<bool>("GoodBJets_pt45");
        const auto& HT_trigger_pt45 = tr.getVar<double>("HT_trigger_pt45");
        const auto& NGoodJets_pt45  = tr.getVar<int>("NGoodJets_pt45");
        const auto& NGoodBJets_pt45 = tr.getVar<int>("NGoodBJets_pt45");
        //const auto& deepESM_val     = tr.getVar<double>("deepESM_val");
        const auto& dR_bjets        = tr.getVar<double>("dR_bjets");               
        const auto& dR_top1_top2    = tr.getVar<double>("dR_top1_top2");
        const auto& topsLV          = tr.getVec<TLorentzVector>("topsLV");
 
        // ------------------------------
        // -- Define Top Tag variables
        // ------------------------------
        const auto& ntops          = tr.getVar<int>("ntops");
        const auto& ntops_1jet     = tr.getVar<int>("ntops_1jet");
        const auto& ntops_2jet     = tr.getVar<int>("ntops_2jet");
        const auto& ntops_3jet     = tr.getVar<int>("ntops_3jet");
        const auto& topsMass       = tr.getVec<double>("topsMass");
        const auto& topsEta        = tr.getVec<double>("topsEta");
        const auto& topsPhi        = tr.getVec<double>("topsPhi");  
        const auto& topsPt         = tr.getVec<double>("topsPt");
        const auto& resolvedMass   = tr.getVar<double>("resolvedMass");
        const auto& resolvedEta    = tr.getVar<double>("resolvedEta");
        const auto& resolvedPhi    = tr.getVar<double>("resolvedPhi");
        const auto& resolvedPt     = tr.getVar<double>("resolvedPt");
        const auto& bestTopMass    = tr.getVar<double>("bestTopMass");
        const auto& bestTopEta     = tr.getVar<double>("bestTopEta");
        const auto& bestTopPt      = tr.getVar<double>("bestTopPt");
       
        const auto& passMadHT      = tr.getVar<bool>("passMadHT");
        const auto& passBaseline0l = tr.getVar<bool>("passBaseline0l_Good");
        const auto& passMETFilters = tr.getVar<bool>("passMETFilters");
        const bool pass_general    = JetID && passMETFilters && passMadHT;
        const bool pass_0l         = NGoodLeptons==0;  
        const bool pass_HT500      = HT_trigger_pt45 > 500;
        const bool pass_ge2b       = NGoodBJets_pt45 >= 2;
        const bool pass_ge6j       = NGoodJets_pt45 >= 6;
        const bool pass_ge2t       = ntops >= 2;
        const bool pass_ge2t1j     = ntops >= 2 && ntops_3jet == 0 && ntops_2jet==0;
        const bool pass_ge2t3j     = ntops >= 2 && ntops_1jet == 0 && ntops_2jet==0;
        const bool pass_ge2t1j3j   = ntops >= 2 && ntops_1jet >= 1 && ntops_3jet >= 1 && ntops_2jet==0;
        const bool pass_ge1dRbjets = dR_bjets >= 1.0;       
    
        // ------------------------------------------
        // -- MT2 or Stealth hemispheres variables
        // ------------------------------------------
        //const auto& stop1_PtRank               = tr.getVar<TLorentzVector>("stop1_PtRank_TaggedTop");
        //const auto& stop2_PtRank               = tr.getVar<TLorentzVector>("stop2_PtRank_TaggedTop");
        //const auto& stop1_MassRank             = tr.getVar<TLorentzVector>("stop1_MassRank_TaggedTop");
        //const auto& stop2_MassRank             = tr.getVar<TLorentzVector>("stop2_MassRank_TaggedTop");
        //const auto& stop1_ScalarPtRank         = tr.getVar<TLorentzVector>("stop1_ScalarPtRank_TaggedTop");
        //const auto& stop2_ScalarPtRank         = tr.getVar<TLorentzVector>("stop2_ScalarPtRank_TaggedTop");
        //const auto& stop1ScalarPt_ScalarPtRank = tr.getVar<double>("stop1ScalarPt_ScalarPtRank_TaggedTop");
        //const auto& stop2ScalarPt_ScalarPtRank = tr.getVar<double>("stop2ScalarPt_ScalarPtRank_TaggedTop");
        //const auto& MT2                        = tr.getVar<double>("MT2_TaggedTop"); 
        //const auto& dR_stop1stop2              = tr.getVar<double>("dR_stop1stop2_TaggedTop");
        //const auto& dPhi_stop1stop2            = tr.getVar<double>("dPhi_stop1stop2_TaggedTop");
        //const auto& difference_stopMasses      = tr.getVar<double>("difference_stopMasses_TaggedTop");
        //const auto& average_stopMasses         = tr.getVar<double>("average_stopMasses_TaggedTop");
        //const auto& relativeDiff_stopMasses    = tr.getVar<double>("relativeDiff_stopMasses_TaggedTop");

        // --------------------
        // -- stop gen level
        // --------------------
        //const auto& StopRecoMT2     = tr.getVar<double>("GM_StopMT2"); // Reco Sum List
        //const auto& Stop1RecoMass   = tr.getVar<double>("GM_Stop1Mass"); // Reco 
        //const auto& Stop2RecoMass   = tr.getVar<double>("GM_Stop2Mass"); // Reco 
        //const auto& Nlino1RecoMass  = tr.getVar<double>("GM_Nlino1Mass"); // Reco 
        //const auto& Nlino2RecoMass  = tr.getVar<double>("GM_Nlino2Mass"); // Reco
        //const auto& Single1RecoMass = tr.getVar<double>("GM_Single1Mass"); // Reco 
        //const auto& Single2RecoMass = tr.getVar<double>("GM_Single2Mass"); // Reco
        //const auto& StopGenMT2      = tr.getVar<double>("GM_StopGenMT2"); // Gen Sum List
        //const auto& Stop1GenMass    = tr.getVar<double>("GM_Stop1GenMass"); // Gen 
        //const auto& Stop2GenMass    = tr.getVar<double>("GM_Stop2GenMass"); // Gen 
        //const auto& Nlino1GenMass   = tr.getVar<double>("GM_Nlino1GenMass"); // Gen
        //const auto& Nlino2GenMass   = tr.getVar<double>("GM_Nlino2GenMass"); // Gen
        //const auto& Single1GenMass  = tr.getVar<double>("GM_Single1GenMass"); // Gen 
        //const auto& Single2GenMass  = tr.getVar<double>("GM_Single2GenMass"); // Gen

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
            const auto& Weight      = tr.getVar<double>("Weight");
            const auto& lumi        = tr.getVar<double>("Lumi");
            eventweight             = lumi*Weight;
        
            bTagScaleFactor         = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            prefiringScaleFactor    = tr.getVar<double>("prefiringScaleFactor");
            puScaleFactor           = tr.getVar<double>("puWeightCorr");
        
            weight *= eventweight*bTagScaleFactor*prefiringScaleFactor*puScaleFactor;
        }

        // ---------------------------------------------
        // -- Calculate DeltaR between tops and bjets
        // ---------------------------------------------
        std::vector<TLorentzVector> bjets;
        for(unsigned int ijet = 0; ijet < Jets.size(); ijet++) {
            if(!GoodBJets_pt45[ijet]) continue;
            bjets.push_back(Jets.at(ijet));        
        }
        std::vector<double> dR_top_bjet;
        for (unsigned int t = 0; t < topsLV.size(); t++) {
            for (unsigned int b = 0; b < bjets.size(); b++) {
                double deltaR = topsLV.at(t).DeltaR(bjets.at(b));
                dR_top_bjet.push_back(deltaR);
            }
        }
  
        // -----------------------------------
        // -- get MT2 hemispheres variables
        // ------------------------------------
        //double stop1Mass_PtRank       = 0.0, stop1Eta_PtRank       = 0.0, stop1Phi_PtRank       = 0.0, stop1Pt_PtRank       = 0.0;
        //double stop2Mass_PtRank       = 0.0, stop2Eta_PtRank       = 0.0, stop2Phi_PtRank       = 0.0, stop2Pt_PtRank       = 0.0;
        //double stop1Mass_MassRank     = 0.0, stop1Eta_MassRank     = 0.0, stop1Phi_MassRank     = 0.0, stop1Pt_MassRank     = 0.0;
        //double stop2Mass_MassRank     = 0.0, stop2Eta_MassRank     = 0.0, stop2Phi_MassRank     = 0.0, stop2Pt_MassRank     = 0.0;
        //double stop1Mass_ScalarPtRank = 0.0, stop1Eta_ScalarPtRank = 0.0, stop1Phi_ScalarPtRank = 0.0, stop1Pt_ScalarPtRank = 0.0;
        //double stop2Mass_ScalarPtRank = 0.0, stop2Eta_ScalarPtRank = 0.0, stop2Phi_ScalarPtRank = 0.0, stop2Pt_ScalarPtRank = 0.0;

        //stop1Mass_PtRank       = stop1_PtRank.M();
        //stop1Eta_PtRank        = stop1_PtRank.Eta();
        //stop1Phi_PtRank        = stop1_PtRank.Phi();       
        //stop1Pt_PtRank         = stop1_PtRank.Pt();
        //stop2Mass_PtRank       = stop2_PtRank.M();
        //stop2Eta_PtRank        = stop2_PtRank.Eta();
        //stop2Phi_PtRank        = stop2_PtRank.Phi(); 
        //stop2Pt_PtRank         = stop2_PtRank.Pt();

        //stop1Mass_MassRank     = stop1_MassRank.M();
        //stop1Eta_MassRank      = stop1_MassRank.Eta();
        //stop1Phi_MassRank      = stop1_MassRank.Phi(); 
        //stop1Pt_MassRank       = stop1_MassRank.Pt();
        //stop2Mass_MassRank     = stop2_MassRank.M();
        //stop2Eta_MassRank      = stop2_MassRank.Eta();
        //stop2Phi_MassRank      = stop2_MassRank.Phi();
        //stop2Pt_MassRank       = stop2_MassRank.Pt();

        //stop1Mass_ScalarPtRank = stop1_ScalarPtRank.M();
        //stop1Eta_ScalarPtRank  = stop1_ScalarPtRank.Eta();
        //stop1Phi_ScalarPtRank  = stop1_ScalarPtRank.Phi();
        //stop1Pt_ScalarPtRank   = stop1_ScalarPtRank.Pt();
        //stop2Mass_ScalarPtRank = stop2_ScalarPtRank.M();
        //stop2Eta_ScalarPtRank  = stop2_ScalarPtRank.Eta();
        //stop2Phi_ScalarPtRank  = stop2_ScalarPtRank.Phi();
        //stop2Pt_ScalarPtRank   = stop2_ScalarPtRank.Pt();

        // -------------------------------------------------
        // -- Make cuts and fill histograms here & cutmap
        // -------------------------------------------------
        const std::map<std::string, bool>& cutmap
        {
            //{"",                                       true                                                            },
            //{"0l",                                     pass_general && pass_0l                                         },
            //{"0l_HT500",                               pass_general && pass_0l && pass_HT500                           },           
            //{"0l_HT500_ge2b",                          pass_general && pass_0l && pass_HT500 && pass_ge2b              },     
            //{"0l_HT500_ge2b_ge6j",                     pass_general && pass_0l && pass_HT500 && pass_ge2b && pass_ge6j },
            
            // >= 2 tops
            //{"0l_HT500_ge2b_ge6j_ge2t",                passBaseline0l && pass_ge2t     },
            //{"0l_HT500_ge2b_ge6j_ge2t1j",              passBaseline0l && pass_ge2t1j   },
            //{"0l_HT500_ge2b_ge6j_ge2t3j",              passBaseline0l && pass_ge2t3j   },
            //{"0l_HT500_ge2b_ge6j_ge2t1j3j",            passBaseline0l && pass_ge2t1j3j }, 
            
            // dR_bjets >= 1
            {"0l_HT500_ge2b_ge6j_ge2t_ge1dRbjets",     passBaseline0l && pass_ge2t && pass_ge1dRbjets     },
            {"0l_HT500_ge2b_ge6j_ge2t1j_ge1dRbjets",   passBaseline0l && pass_ge2t1j && pass_ge1dRbjets   },
            {"0l_HT500_ge2b_ge6j_ge2t3j_ge1dRbjets",   passBaseline0l && pass_ge2t3j && pass_ge1dRbjets   },
            {"0l_HT500_ge2b_ge6j_ge2t1j3j_ge1dRbjets", passBaseline0l && pass_ge2t1j3j && pass_ge1dRbjets },
       
            // pt Rank stopMasses cuts
            //{"baseline_0l_stopMassesg200_PtRank",          passBaseline0l && pass_ge2t && pass_ge1dRbjets && stop1Mass_PtRank > 200 && stop2Mass_PtRank > 200     },
            //{"baseline_0l_ge2t1j_stopMassesg200_PtRank",   passBaseline0l && pass_ge2t1j && pass_ge1dRbjets && stop1Mass_PtRank > 200 && stop2Mass_PtRank > 200   },
            //{"baseline_0l_ge2t3j_stopMassesg200_PtRank",   passBaseline0l && pass_ge2t3j && pass_ge1dRbjets  && stop1Mass_PtRank > 200 && stop2Mass_PtRank > 200  },
            //{"baseline_0l_ge2t1j3j_stopMassesg200_PtRank", passBaseline0l && pass_ge2t1j3j && pass_ge1dRbjets && stop1Mass_PtRank > 200 && stop2Mass_PtRank > 200 },

            //{"baseline_0l_stopMassesle200_PtRank",          passBaseline0l && pass_ge2t && pass_ge1dRbjets && stop1Mass_PtRank < 200 && stop2Mass_PtRank < 200     },
            //{"baseline_0l_ge2t1j_stopMassesle200_PtRank",   passBaseline0l && pass_ge2t1j && pass_ge1dRbjets && stop1Mass_PtRank < 200 && stop2Mass_PtRank < 200   },
            //{"baseline_0l_ge2t3j_stopMassesle200_PtRank",   passBaseline0l && pass_ge2t3j && pass_ge1dRbjets  && stop1Mass_PtRank < 200 && stop2Mass_PtRank < 200  },
            //{"baseline_0l_ge2t1j3j_stopMassesle200_PtRank", passBaseline0l && pass_ge2t1j3j && pass_ge1dRbjets && stop1Mass_PtRank < 200 && stop2Mass_PtRank < 200 },

            // mass Rank stopMasses cuts
            //{"baseline_0l_stopMassesg200_MassRank",          passBaseline0l && pass_ge2t && pass_ge1dRbjets && stop1Mass_MassRank > 200 && stop2Mass_MassRank > 200     },
            //{"baseline_0l_ge2t1j_stopMassesg200_MassRank",   passBaseline0l && pass_ge2t1j && pass_ge1dRbjets && stop1Mass_MassRank > 200 && stop2Mass_MassRank > 200   },
            //{"baseline_0l_ge2t3j_stopMassesg200_MassRank",   passBaseline0l && pass_ge2t3j && pass_ge1dRbjets  && stop1Mass_MassRank > 200 && stop2Mass_MassRank > 200  },
            //{"baseline_0l_ge2t1j3j_stopMassesg200_MassRank", passBaseline0l && pass_ge2t1j3j && pass_ge1dRbjets && stop1Mass_MassRank > 200 && stop2Mass_MassRank > 200 },

            //{"baseline_0l_stopMassesle200_MassRank",          passBaseline0l && pass_ge2t && pass_ge1dRbjets && stop1Mass_MassRank < 200 && stop2Mass_MassRank < 200     },
            //{"baseline_0l_ge2t1j_stopMassesle200_MassRank",   passBaseline0l && pass_ge2t1j && pass_ge1dRbjets && stop1Mass_MassRank < 200 && stop2Mass_MassRank < 200   },
            //{"baseline_0l_ge2t3j_stopMassesle200_MassRank",   passBaseline0l && pass_ge2t3j && pass_ge1dRbjets  && stop1Mass_MassRank < 200 && stop2Mass_MassRank < 200  },
            //{"baseline_0l_ge2t1j3j_stopMassesle200_MassRank", passBaseline0l && pass_ge2t1j3j && pass_ge1dRbjets && stop1Mass_MassRank < 200 && stop2Mass_MassRank < 200 },
            
            // scalar pt Rank stopMasses cuts
            //{"baseline_0l_stopMassesg200_ScalarPtRank",          passBaseline0l && pass_ge2t && pass_ge1dRbjets && stop1Mass_ScalarPtRank > 200 && stop2Mass_ScalarPtRank > 200     },
            //{"baseline_0l_ge2t1j_stopMassesg200_ScalarPtRank",   passBaseline0l && pass_ge2t1j && pass_ge1dRbjets && stop1Mass_ScalarPtRank > 200 && stop2Mass_ScalarPtRank > 200   },
            //{"baseline_0l_ge2t3j_stopMassesg200_ScalarPtRank",   passBaseline0l && pass_ge2t3j && pass_ge1dRbjets  && stop1Mass_ScalarPtRank > 200 && stop2Mass_ScalarPtRank > 200  },
            //{"baseline_0l_ge2t1j3j_stopMassesg200_ScalarPtRank", passBaseline0l && pass_ge2t1j3j && pass_ge1dRbjets && stop1Mass_ScalarPtRank > 200 && stop2Mass_ScalarPtRank > 200 },

            //{"baseline_0l_stopMassesle200_ScalarPtRank",          passBaseline0l && pass_ge2t && pass_ge1dRbjets && stop1Mass_ScalarPtRank < 200 && stop2Mass_ScalarPtRank < 200     },
            //{"baseline_0l_ge2t1j_stopMassesle200_ScalarPtRank",   passBaseline0l && pass_ge2t1j && pass_ge1dRbjets && stop1Mass_ScalarPtRank < 200 && stop2Mass_ScalarPtRank < 200   },
            //{"baseline_0l_ge2t3j_stopMassesle200_ScalarPtRank",   passBaseline0l && pass_ge2t3j && pass_ge1dRbjets  && stop1Mass_ScalarPtRank < 200 && stop2Mass_ScalarPtRank < 200  },
            //{"baseline_0l_ge2t1j3j_stopMassesle200_ScalarPtRank", passBaseline0l && pass_ge2t1j3j && pass_ge1dRbjets && stop1Mass_ScalarPtRank < 200 && stop2Mass_ScalarPtRank < 200 },

            // NJets cuts for stop hemispheres
            //{"baseline_0l_Njet7",  passBaseline0l && pass_ge2t && pass_ge1dRbjets && NGoodJets_pt45 == 7  },
            //{"baseline_0l_Njet8",  passBaseline0l && pass_ge2t && pass_ge1dRbjets && NGoodJets_pt45 == 8  },
            //{"baseline_0l_Njet9",  passBaseline0l && pass_ge2t && pass_ge1dRbjets && NGoodJets_pt45 == 9  },
            //{"baseline_0l_Njet10", passBaseline0l && pass_ge2t && pass_ge1dRbjets && NGoodJets_pt45 == 10 },
            //{"baseline_0l_Njet11", passBaseline0l && pass_ge2t && pass_ge1dRbjets && NGoodJets_pt45 == 11 },
            //{"baseline_0l_Njet12", passBaseline0l && pass_ge2t && pass_ge1dRbjets && NGoodJets_pt45 == 12 },
            //{"baseline_0l_Njet13", passBaseline0l && pass_ge2t && pass_ge1dRbjets && NGoodJets_pt45 == 13 },
            //{"baseline_0l_Njet14", passBaseline0l && pass_ge2t && pass_ge1dRbjets && NGoodJets_pt45 == 14 },
            //{"baseline_0l_Njet15", passBaseline0l && pass_ge2t && pass_ge1dRbjets && NGoodJets_pt45 == 15 },
    };

    if (!inithisto) {
        InitHistos(cutmap);
        inithisto = true;
    }

    my_histos["EventCounter"]->Fill( eventCounter );

    // --------------------------------
    // -- Fill the cutmap histograms
    // --------------------------------     
    for (const auto& cutVar: cutmap) {        
        if (cutVar.second) {
            my_histos["h_ntops_"+cutVar.first]->Fill( ntops, weight );
            my_histos["h_njets_"+cutVar.first]->Fill( NGoodJets_pt45, weight );
            my_histos["h_nbjets_"+cutVar.first]->Fill( NGoodBJets_pt45, weight );
            my_histos["h_ht_"+cutVar.first]->Fill( HT_trigger_pt45, weight );
            my_histos["h_met_"+cutVar.first]->Fill( MET, weight );
            
            // -----------------------------
            // -- jets & bjets mass & pT
            // -----------------------------
            //for(int ijet = 0; ijet < Jets.size(); ijet++) {
            //    if(!GoodJets_pt45[ijet]) continue;
            //    my_histos["h_jetsPt_"+cutVar.first]->Fill(Jets.at(ijet).Pt(), weight);
            //    my_histos["h_jetsMass_"+cutVar.first]->Fill(Jets.at(ijet).M(), weight);                 
    
            //    if(!GoodBJets_pt45[ijet]) continue;
            //    const TLorentzVector& bjet = Jets.at(ijet);                     
            //    my_histos["h_bjetsPt_"+cutVar.first]->Fill(bjet.Pt(), weight);
            //    my_histos["h_bjetsMass_"+cutVar.first]->Fill(bjet.M(), weight);
            //}
        
            // --------------------------------------
            // -- get top jets' Mass, Eta, Phi, Pt 
            // --------------------------------------
            for (unsigned int itops = 0; itops < topsMass.size(); itops++) 
            {
                my_histos["h_topsMass_"+cutVar.first]->Fill( topsMass.at(itops), weight );
            }
        
            for (unsigned int itops = 0; itops < topsEta.size(); itops++) {
                my_histos["h_topsEta_"+cutVar.first]->Fill( topsEta.at(itops), weight );
            }
       
            for (unsigned int itops = 0; itops < topsPhi.size(); itops++) {
                my_histos["h_topsPhi_"+cutVar.first]->Fill( topsPhi.at(itops), weight );
            }
 
            for (unsigned int itops = 0; itops < topsPt.size(); itops++) {
                my_histos["h_topsPt_"+cutVar.first]->Fill( topsPt.at(itops), weight );
            }      

            // get the resolved jets' Mass, Eta, Phi, Pt
            my_histos["h_resolvedMass_"+cutVar.first]->Fill( resolvedMass, weight );
            my_histos["h_resolvedEta_"+cutVar.first]->Fill( resolvedEta, weight );
            my_histos["h_resolvedPhi_"+cutVar.first]->Fill( resolvedPhi, weight );
            my_histos["h_resolvedPt_"+cutVar.first]->Fill( resolvedPt, weight ); 
     
            my_histos["h_bestTopMass_"+cutVar.first]->Fill( bestTopMass, weight );
            my_histos["h_bestTopEta_"+cutVar.first]->Fill( bestTopEta, weight );
            my_histos["h_bestTopPt_"+cutVar.first]->Fill( bestTopPt, weight );
            my_histos["h_dR_bjets_"+cutVar.first]->Fill( dR_bjets, weight );
            my_histos["h_dR_top1_top2_"+cutVar.first]->Fill( dR_top1_top2, weight );
    
            // ---------------------------------
            // -- deltaR between top and bjet
            // ---------------------------------
            for (unsigned int idR = 0; idR < dR_top_bjet.size(); idR++) {
                my_histos["h_dR_tops_bjets_"+cutVar.first]->Fill( dR_top_bjet.at(idR), weight );        
            }
         
            //my_2d_histos["h_njets_MVA_"+cutVar.first]->Fill( NGoodJets_pt45, deepESM_val, weight );
            //my_2d_histos["h_njets_MVA_"+cutVar.first]->GetXaxis()->SetTitle("N_{J}");
            //my_2d_histos["h_njets_MVA_"+cutVar.first]->GetYaxis()->SetTitle("MVA");
            //my_2d_histos["h_njets_dR_bjets_"+cutVar.first]->Fill( dR_bjets, NGoodJets_pt45, weight );
            //my_2d_histos["h_njets_dR_bjets_"+cutVar.first]->GetXaxis()->SetTitle("#DeltaR_{bjets}");
            //my_2d_histos["h_njets_dR_bjets_"+cutVar.first]->GetYaxis()->SetTitle("N_{J}");

            // --------------------------
            // -- stop MT2 hemispheres 
            // --------------------------
            // ptRank
            //my_histos["h_stop1Mass_PtRank_"+cutVar.first]->Fill( stop1Mass_PtRank, weight );
            //my_histos["h_stop1Eta_PtRank_"+cutVar.first]->Fill( stop1Eta_PtRank, weight );
            //my_histos["h_stop1Phi_PtRank_"+cutVar.first]->Fill( stop1Phi_PtRank, weight );
            //my_histos["h_stop1Pt_PtRank_"+cutVar.first]->Fill( stop1Pt_PtRank, weight );
            //my_histos["h_stop2Mass_PtRank_"+cutVar.first]->Fill( stop2Mass_PtRank, weight );
            //my_histos["h_stop2Eta_PtRank_"+cutVar.first]->Fill( stop2Eta_PtRank, weight );
            //my_histos["h_stop2Phi_PtRank_"+cutVar.first]->Fill( stop2Phi_PtRank, weight );
            //my_histos["h_stop2Pt_PtRank_"+cutVar.first]->Fill( stop2Pt_PtRank, weight );
            // massRank
            //my_histos["h_stop1Mass_MassRank_"+cutVar.first]->Fill( stop1Mass_MassRank, weight );
            //my_histos["h_stop1Eta_MassRank_"+cutVar.first]->Fill( stop1Eta_MassRank, weight );
            //my_histos["h_stop1Phi_MassRank_"+cutVar.first]->Fill( stop1Phi_MassRank, weight );
            //my_histos["h_stop1Pt_MassRank_"+cutVar.first]->Fill( stop1Pt_MassRank, weight );
            //my_histos["h_stop2Mass_MassRank_"+cutVar.first]->Fill( stop2Mass_MassRank, weight );
            //my_histos["h_stop2Eta_MassRank_"+cutVar.first]->Fill( stop2Eta_MassRank, weight );
            //my_histos["h_stop2Phi_MassRank_"+cutVar.first]->Fill( stop2Phi_MassRank, weight );
            //my_histos["h_stop2Pt_MassRank_"+cutVar.first]->Fill( stop2Pt_MassRank, weight );
            // scalarPtRank
            //my_histos["h_stop1Mass_ScalarPtRank_"+cutVar.first]->Fill( stop1Mass_ScalarPtRank, weight );
            //my_histos["h_stop1Eta_ScalarPtRank_"+cutVar.first]->Fill( stop1Eta_ScalarPtRank, weight );
            //my_histos["h_stop1Phi_ScalarPtRank_"+cutVar.first]->Fill( stop1Phi_ScalarPtRank, weight );
            //my_histos["h_stop1Pt_ScalarPtRank_"+cutVar.first]->Fill( stop1Pt_ScalarPtRank, weight );
            //my_histos["h_stop2Mass_ScalarPtRank_"+cutVar.first]->Fill( stop2Mass_ScalarPtRank, weight );
            //my_histos["h_stop2Eta_ScalarPtRank_"+cutVar.first]->Fill( stop2Eta_ScalarPtRank, weight );
            //my_histos["h_stop2Phi_ScalarPtRank_"+cutVar.first]->Fill( stop2Phi_ScalarPtRank, weight );
            //my_histos["h_stop2Pt_ScalarPtRank_"+cutVar.first]->Fill( stop2Pt_ScalarPtRank, weight );
            //my_histos["h_stop1ScalarPt_ScalarPtRank_"+cutVar.first]->Fill( stop1ScalarPt_ScalarPtRank, weight );
            //my_histos["h_stop2ScalarPt_ScalarPtRank_"+cutVar.first]->Fill( stop2ScalarPt_ScalarPtRank, weight );
            //my_histos["h_MT2_"+cutVar.first]->Fill( MT2, weight );
            //my_histos["h_dR_stop1stop2_"+cutVar.first]->Fill( dR_stop1stop2, weight );
            //my_histos["h_dPhi_stop1stop2_"+cutVar.first]->Fill( dPhi_stop1stop2, weight );
            //my_histos["h_difference_stopMasses_"+cutVar.first]->Fill( difference_stopMasses, weight );
            //my_histos["h_average_stopMasses_"+cutVar.first]->Fill( average_stopMasses, weight );
            //my_histos["h_relativeDiff_stopMasses_"+cutVar.first]->Fill( relativeDiff_stopMasses, weight );
            
            // stop1VSstop2 Mass, Eta, Phi, Pt
            //my_2d_histos["h_Mass_stop1vsstop2_PtRank_"+cutVar.first]->Fill(stop1Mass_PtRank, stop2Mass_PtRank, weight);
            //my_2d_histos["h_Mass_stop1vsstop2_PtRank_"+cutVar.first]->GetXaxis()->SetTitle("Pt Rank M_{#tildet}_{1}");
            //my_2d_histos["h_Mass_stop1vsstop2_PtRank_"+cutVar.first]->GetYaxis()->SetTitle("Pt Rank M_{#tildet}_{2}");
            //my_2d_histos["h_Eta_stop1vsstop2_PtRank_"+cutVar.first]->Fill(stop1Eta_PtRank, stop2Eta_PtRank, weight);
            //my_2d_histos["h_Eta_stop1vsstop2_PtRank_"+cutVar.first]->GetXaxis()->SetTitle("Pt Rank #eta_{#tildet}_{1}");
            //my_2d_histos["h_Eta_stop1vsstop2_PtRank_"+cutVar.first]->GetYaxis()->SetTitle("Pt Rank #eta_{#tildet}_{2}");
            //my_2d_histos["h_Phi_stop1vsstop2_PtRank_"+cutVar.first]->Fill(stop1Phi_PtRank, stop2Phi_PtRank, weight);
            //my_2d_histos["h_Phi_stop1vsstop2_PtRank_"+cutVar.first]->GetXaxis()->SetTitle("Pt Rank #phi_{#tildet}_{1}");
            //my_2d_histos["h_Phi_stop1vsstop2_PtRank_"+cutVar.first]->GetYaxis()->SetTitle("Pt Rank #phi_{#tildet}_{2}");
            //my_2d_histos["h_Pt_stop1vsstop2_PtRank_"+cutVar.first]->Fill(stop1Pt_PtRank, stop2Pt_PtRank, weight);
            //my_2d_histos["h_Pt_stop1vsstop2_PtRank_"+cutVar.first]->GetXaxis()->SetTitle("Pt Rank pT_{#tildet}_{1}");
            //my_2d_histos["h_Pt_stop1vsstop2_PtRank_"+cutVar.first]->GetYaxis()->SetTitle("Pt Rank pT_{#tildet}_{2}");
            //my_2d_histos["h_Mass_stop1vsstop2_MassRank_"+cutVar.first]->Fill(stop1Mass_MassRank, stop2Mass_MassRank, weight);
            //my_2d_histos["h_Mass_stop1vsstop2_MassRank_"+cutVar.first]->GetXaxis()->SetTitle("Mass Rank M_{#tildet}_{1}");
            //my_2d_histos["h_Mass_stop1vsstop2_MassRank_"+cutVar.first]->GetYaxis()->SetTitle("Mass Rank M_{#tildet}_{2}");
            //my_2d_histos["h_Eta_stop1vsstop2_MassRank_"+cutVar.first]->Fill(stop1Eta_MassRank, stop2Eta_MassRank, weight);
            //my_2d_histos["h_Eta_stop1vsstop2_MassRank_"+cutVar.first]->GetXaxis()->SetTitle("Mass Rank #eta_{#tildet}_{1}");
            //my_2d_histos["h_Eta_stop1vsstop2_MassRank_"+cutVar.first]->GetYaxis()->SetTitle("Mass Rank #eta_{#tildet}_{2}");
            //my_2d_histos["h_Phi_stop1vsstop2_MassRank_"+cutVar.first]->Fill(stop1Phi_MassRank, stop2Phi_MassRank, weight);
            //my_2d_histos["h_Phi_stop1vsstop2_MassRank_"+cutVar.first]->GetXaxis()->SetTitle("Mass Rank #phi_{#tildet}_{1}");
            //my_2d_histos["h_Phi_stop1vsstop2_MassRank_"+cutVar.first]->GetYaxis()->SetTitle("Mass Rank #phi_{#tildet}_{2}");
            //my_2d_histos["h_Pt_stop1vsstop2_MassRank_"+cutVar.first]->Fill(stop1Pt_MassRank, stop2Pt_MassRank, weight);
            //my_2d_histos["h_Pt_stop1vsstop2_MassRank_"+cutVar.first]->GetXaxis()->SetTitle("Mass Rank pT_{#tildet}_{1}");
            //my_2d_histos["h_Pt_stop1vsstop2_MassRank_"+cutVar.first]->GetYaxis()->SetTitle("Mass Rank pT_{#tildet}_{2}");
            //my_2d_histos["h_Mass_stop1vsstop2_ScalarPtRank_"+cutVar.first]->Fill(stop1Mass_ScalarPtRank, stop2Mass_ScalarPtRank, weight);
            //my_2d_histos["h_Mass_stop1vsstop2_ScalarPtRank_"+cutVar.first]->GetXaxis()->SetTitle("ScalarPt Rank M_{#tildet}_{1}");
            //my_2d_histos["h_Mass_stop1vsstop2_ScalarPtRank_"+cutVar.first]->GetYaxis()->SetTitle("ScalarPt Rank M_{#tildet}_{2}");
            //my_2d_histos["h_Eta_stop1vsstop2_ScalarPtRank_"+cutVar.first]->Fill(stop1Eta_ScalarPtRank, stop2Eta_ScalarPtRank, weight);
            //my_2d_histos["h_Eta_stop1vsstop2_ScalarPtRank_"+cutVar.first]->GetXaxis()->SetTitle("ScalarPt Rank #eta_{#tildet}_{1}");
            //my_2d_histos["h_Eta_stop1vsstop2_ScalarPtRank_"+cutVar.first]->GetYaxis()->SetTitle("ScalarPt Rank #eta_{#tildet}_{2}");
            //my_2d_histos["h_Phi_stop1vsstop2_ScalarPtRank_"+cutVar.first]->Fill(stop1Phi_ScalarPtRank, stop2Phi_ScalarPtRank, weight);
            //my_2d_histos["h_Phi_stop1vsstop2_ScalarPtRank_"+cutVar.first]->GetXaxis()->SetTitle("ScalarPt Rank #phi_{#tildet}_{1}");
            //my_2d_histos["h_Phi_stop1vsstop2_ScalarPtRank_"+cutVar.first]->GetYaxis()->SetTitle("ScalarPt Rank #phi_{#tildet}_{2}");
            //my_2d_histos["h_Pt_stop1vsstop2_ScalarPtRank_"+cutVar.first]->Fill(stop1Pt_ScalarPtRank, stop2Pt_ScalarPtRank, weight);
            //my_2d_histos["h_Pt_stop1vsstop2_ScalarPtRank_"+cutVar.first]->GetXaxis()->SetTitle("ScalarPt Rank pT_{#tildet}_{1}");
            //my_2d_histos["h_Pt_stop1vsstop2_ScalarPtRank_"+cutVar.first]->GetYaxis()->SetTitle("ScalarPt Rank pT_{#tildet}_{2}");
            
            // NJetsVSstops
            //my_2d_histos["h_Mass_NJetsVSstop1_PtRank_"+cutVar.first]->Fill( NGoodJets_pt45, stop1Mass_PtRank, weight );
            //my_2d_histos["h_Mass_NJetsVSstop1_PtRank_"+cutVar.first]->GetXaxis()->SetTitle("N_{J}");
            //my_2d_histos["h_Mass_NJetsVSstop1_PtRank_"+cutVar.first]->GetYaxis()->SetTitle("Pt Rank M_{#tildet}_{1}"); 
            //my_2d_histos["h_Mass_NJetsVSstop2_PtRank_"+cutVar.first]->Fill( NGoodJets_pt45, stop2Mass_PtRank, weight );
            //my_2d_histos["h_Mass_NJetsVSstop2_PtRank_"+cutVar.first]->GetXaxis()->SetTitle("N_{J}");
            //my_2d_histos["h_Mass_NJetsVSstop2_PtRank_"+cutVar.first]->GetYaxis()->SetTitle("Pt Rank M_{#tildet}_{2}");               
            //my_2d_histos["h_Mass_NJetsVSstop1_MassRank_"+cutVar.first]->Fill( NGoodJets_pt45, stop1Mass_MassRank, weight );
            //my_2d_histos["h_Mass_NJetsVSstop1_MassRank_"+cutVar.first]->GetXaxis()->SetTitle("N_{J}");
            //my_2d_histos["h_Mass_NJetsVSstop1_MassRank_"+cutVar.first]->GetYaxis()->SetTitle("Mass Rank M_{#tildet}_{1}");
            //my_2d_histos["h_Mass_NJetsVSstop2_MassRank_"+cutVar.first]->Fill( NGoodJets_pt45, stop2Mass_MassRank, weight );
            //my_2d_histos["h_Mass_NJetsVSstop2_MassRank_"+cutVar.first]->GetXaxis()->SetTitle("N_{J}");
            //my_2d_histos["h_Mass_NJetsVSstop2_MassRank_"+cutVar.first]->GetYaxis()->SetTitle("Mass Rank M_{#tildet}_{2}");
            //my_2d_histos["h_Mass_NJetsVSstop1_ScalarPtRank_"+cutVar.first]->Fill( NGoodJets_pt45, stop1Mass_ScalarPtRank, weight );
            //my_2d_histos["h_Mass_NJetsVSstop1_ScalarPtRank_"+cutVar.first]->GetXaxis()->SetTitle("N_{J}");
            //my_2d_histos["h_Mass_NJetsVSstop1_ScalarPtRank_"+cutVar.first]->GetYaxis()->SetTitle("ScalarPt Rank M_{#tildet}_{1}");
            //my_2d_histos["h_Mass_NJetsVSstop2_ScalarPtRank_"+cutVar.first]->Fill( NGoodJets_pt45, stop2Mass_ScalarPtRank, weight );
            //my_2d_histos["h_Mass_NJetsVSstop2_ScalarPtRank_"+cutVar.first]->GetXaxis()->SetTitle("N_{J}");
            //my_2d_histos["h_Mass_NJetsVSstop2_ScalarPtRank_"+cutVar.first]->GetYaxis()->SetTitle("ScalarPt Rank M_{#tildet}_{2}");
            //my_2d_histos["h_stopMasses_diffVSavg_"+cutVar.first]->Fill( difference_stopMasses, average_stopMasses, weight);
            //my_2d_histos["h_stopMasses_diffVSavg_"+cutVar.first]->GetXaxis()->SetTitle("difference");
            //my_2d_histos["h_stopMasses_diffVSavg_"+cutVar.first]->GetYaxis()->SetTitle("average");
            
            // stop All Pt combinations
            //my_2d_histos["h_ScalarPt_stop1vsstop2_ScalarPtRank_"+cutVar.first]->Fill( stop1ScalarPt_ScalarPtRank, stop2ScalarPt_ScalarPtRank, weight );
            //my_2d_histos["h_ScalarPt_stop1vsstop2_ScalarPtRank_"+cutVar.first]->GetXaxis()->SetTitle("ScalarPt Rank Scalar pT_{#tildet}_{1}");
            //my_2d_histos["h_ScalarPt_stop1vsstop2_ScalarPtRank_"+cutVar.first]->GetYaxis()->SetTitle("ScalarPt Rank Scalar pT_{#tildet}_{2}");
            //my_2d_histos["h_Pt1_PtRankVsScalarPtRank_"+cutVar.first]->Fill( stop1Pt_PtRank, stop1Pt_ScalarPtRank, weight );
            //my_2d_histos["h_Pt1_PtRankVsScalarPtRank_"+cutVar.first]->GetXaxis()->SetTitle("Pt Rank pT_{#tildet}_{1}");
            //my_2d_histos["h_Pt1_PtRankVsScalarPtRank_"+cutVar.first]->GetYaxis()->SetTitle("ScalarPt Rank pT_{#tildet}_{1}");
            //my_2d_histos["h_Pt2_PtRankVsScalarPtRank_"+cutVar.first]->Fill( stop2Pt_PtRank, stop2Pt_ScalarPtRank, weight );
            //my_2d_histos["h_Pt2_PtRankVsScalarPtRank_"+cutVar.first]->GetXaxis()->SetTitle("Pt Rank pT_{#tildet}_{2}");
            //my_2d_histos["h_Pt2_PtRankVsScalarPtRank_"+cutVar.first]->GetYaxis()->SetTitle("ScalarPt Rank pT_{#tildet}_{2}");
            //my_2d_histos["h_Pt1_PtRankVsScalarPt1_ScalarPtRank_"+cutVar.first]->Fill( stop1Pt_PtRank, stop1ScalarPt_ScalarPtRank, weight );
            //my_2d_histos["h_Pt1_PtRankVsScalarPt1_ScalarPtRank_"+cutVar.first]->GetXaxis()->SetTitle("Pt Rank pT_{#tildet}_{1}");
            //my_2d_histos["h_Pt1_PtRankVsScalarPt1_ScalarPtRank_"+cutVar.first]->GetYaxis()->SetTitle("ScalarPt Rank Scalar pT_{#tildet}_{1}");
            //my_2d_histos["h_Pt2_PtRankVsScalarPt2_ScalarPtRank"+cutVar.first]->Fill( stop2Pt_PtRank, stop2ScalarPt_ScalarPtRank, weight );
            //my_2d_histos["h_Pt2_PtRankVsScalarPt2_ScalarPtRank"+cutVar.first]->GetXaxis()->SetTitle("Pt Rank pT_{#tildet}_{2}");
            //my_2d_histos["h_Pt2_PtRankVsScalarPt2_ScalarPtRank"+cutVar.first]->GetYaxis()->SetTitle("ScalarPt Rank Scalar pT_{#tildet}_{2}");
            
            // stops MassVsPt
            //my_2d_histos["h_stop1_MassVsPt_PtRank_"+cutVar.first]->Fill(stop1Mass_PtRank, stop1Pt_PtRank, weight);
            //my_2d_histos["h_stop1_MassVsPt_PtRank_"+cutVar.first]->GetXaxis()->SetTitle("Pt Rank M_{#tildet}_{1}");
            //my_2d_histos["h_stop1_MassVsPt_PtRank_"+cutVar.first]->GetYaxis()->SetTitle("Pt Rank pT_{#tildet}_{1}");
            //my_2d_histos["h_stop2_MassVsPt_PtRank_"+cutVar.first]->Fill(stop2Mass_PtRank, stop2Pt_PtRank, weight);
            //my_2d_histos["h_stop2_MassVsPt_PtRank_"+cutVar.first]->GetXaxis()->SetTitle("Pt Rank M_{#tildet}_{2}");
            //my_2d_histos["h_stop2_MassVsPt_PtRank_"+cutVar.first]->GetYaxis()->SetTitle("Pt Rank pT_{#tildet}_{2}");
            //my_2d_histos["h_stop1_MassVsPt_MassRank_"+cutVar.first]->Fill(stop1Mass_MassRank, stop1Pt_MassRank, weight);
            //my_2d_histos["h_stop1_MassVsPt_MassRank_"+cutVar.first]->GetXaxis()->SetTitle("Mass Rank M_{#tildet}_{1}");
            //my_2d_histos["h_stop1_MassVsPt_MassRank_"+cutVar.first]->GetYaxis()->SetTitle("Mass Rank pT_{#tildet}_{1}");
            //my_2d_histos["h_stop2_MassVsPt_MassRank_"+cutVar.first]->Fill(stop2Mass_MassRank, stop2Pt_MassRank, weight);
            //my_2d_histos["h_stop2_MassVsPt_MassRank_"+cutVar.first]->GetXaxis()->SetTitle("Mass Rank M_{#tildet}_{2}");
            //my_2d_histos["h_stop2_MassVsPt_MassRank_"+cutVar.first]->GetYaxis()->SetTitle("Mass Rank pT_{#tildet}_{2}");
            //my_2d_histos["h_stop1_MassVsPt_ScalarPtRank_"+cutVar.first]->Fill(stop1Mass_ScalarPtRank, stop1Pt_ScalarPtRank, weight);
            //my_2d_histos["h_stop1_MassVsPt_ScalarPtRank_"+cutVar.first]->GetXaxis()->SetTitle("ScalarPt Rank M_{#tildet}_{1}");
            //my_2d_histos["h_stop1_MassVsPt_ScalarPtRank_"+cutVar.first]->GetYaxis()->SetTitle("ScalarPt Rank pT_{#tildet}_{1}");
            //my_2d_histos["h_stop2_MassVsPt_ScalarPtRank_"+cutVar.first]->Fill(stop2Mass_ScalarPtRank, stop2Pt_ScalarPtRank, weight);
            //my_2d_histos["h_stop2_MassVsPt_ScalarPtRank_"+cutVar.first]->GetXaxis()->SetTitle("ScalarPt Rank M_{#tildet}_{2}");
            //my_2d_histos["h_stop2_MassVsPt_ScalarPtRank_"+cutVar.first]->GetYaxis()->SetTitle("ScalarPt Rank pT_{#tildet}_{2}");
            //my_2d_histos["h_stop1_MassVsScalarPt_ScalarPtRank_"+cutVar.first]->Fill(stop1Mass_ScalarPtRank, stop1ScalarPt_ScalarPtRank, weight);
            //my_2d_histos["h_stop1_MassVsScalarPt_ScalarPtRank_"+cutVar.first]->GetXaxis()->SetTitle("ScalarPt Rank M_{#tildet}_{1}");
            //my_2d_histos["h_stop1_MassVsScalarPt_ScalarPtRank_"+cutVar.first]->GetYaxis()->SetTitle("ScalarPt Rank Scalar pT_{#tildet}_{1}");
            //my_2d_histos["h_stop2_MassVsScalarPt_ScalarPtRank_"+cutVar.first]->Fill(stop2Mass_ScalarPtRank, stop1ScalarPt_ScalarPtRank, weight);
            //my_2d_histos["h_stop2_MassVsScalarPt_ScalarPtRank_"+cutVar.first]->GetXaxis()->SetTitle("ScalarPt Rank M_{#tildet}_{2}");
            //my_2d_histos["h_stop2_MassVsScalarPt_ScalarPtRank_"+cutVar.first]->GetYaxis()->SetTitle("ScalarPt Rank Scalar pT_{#tildet}_{2}");

            // MT2
            //my_2d_histos["h_NJetsVsMT2_"+cutVar.first]->Fill( NGoodJets_pt45, MT2, weight );
            //my_2d_histos["h_NJetsVsMT2_"+cutVar.first]->GetXaxis()->SetTitle("N_{J}");
            //my_2d_histos["h_NJetsVsMT2_"+cutVar.first]->GetYaxis()->SetTitle("MT2");
            //my_2d_histos["h_Mass_MT2vsstop1_PtRank_"+cutVar.first]->Fill( MT2, stop1Mass_PtRank, weight );
            //my_2d_histos["h_Mass_MT2vsstop1_PtRank_"+cutVar.first]->GetXaxis()->SetTitle("MT2");
            //my_2d_histos["h_Mass_MT2vsstop1_PtRank_"+cutVar.first]->GetYaxis()->SetTitle("Pt Rank M_{#tildet}_{1} [GeV]");
            //my_2d_histos["h_Mass_MT2vsstop2_PtRank_"+cutVar.first]->Fill( MT2, stop2Mass_PtRank, weight );
            //my_2d_histos["h_Mass_MT2vsstop2_PtRank_"+cutVar.first]->GetXaxis()->SetTitle("MT2");
            //my_2d_histos["h_Mass_MT2vsstop2_PtRank_"+cutVar.first]->GetYaxis()->SetTitle("Pt Rank M_{#tildet}_{2} [GeV]");
            //my_2d_histos["h_Mass_MT2vsstop1_MassRank_"+cutVar.first]->Fill( MT2, stop1Mass_MassRank, weight );
            //my_2d_histos["h_Mass_MT2vsstop1_MassRank_"+cutVar.first]->GetXaxis()->SetTitle("MT2");
            //my_2d_histos["h_Mass_MT2vsstop1_MassRank_"+cutVar.first]->GetYaxis()->SetTitle("Mass Rank M_{#tildet}_{1} [GeV]");
            //my_2d_histos["h_Mass_MT2vsstop2_MassRank_"+cutVar.first]->Fill( MT2, stop2Mass_MassRank, weight );
            //my_2d_histos["h_Mass_MT2vsstop2_MassRank_"+cutVar.first]->GetXaxis()->SetTitle("MT2");
            //my_2d_histos["h_Mass_MT2vsstop2_MassRank_"+cutVar.first]->GetYaxis()->SetTitle("Mass Rank M_{#tildet}_{2} [GeV]");
            //my_2d_histos["h_Mass_MT2vsstop1_ScalarPtRank_"+cutVar.first]->Fill( MT2, stop1Mass_ScalarPtRank, weight );
            //my_2d_histos["h_Mass_MT2vsstop1_ScalarPtRank_"+cutVar.first]->GetXaxis()->SetTitle("MT2");
            //my_2d_histos["h_Mass_MT2vsstop1_ScalarPtRank_"+cutVar.first]->GetYaxis()->SetTitle("ScalarPt Rank M_{#tildet}_{1} [GeV]");
            //my_2d_histos["h_Mass_MT2vsstop2_ScalarPtRank_"+cutVar.first]->Fill( MT2, stop2Mass_ScalarPtRank, weight );
            //my_2d_histos["h_Mass_MT2vsstop2_ScalarPtRank_"+cutVar.first]->GetXaxis()->SetTitle("MT2");
            //my_2d_histos["h_Mass_MT2vsstop2_ScalarPtRank_"+cutVar.first]->GetYaxis()->SetTitle("ScalarPt Rank M_{#tildet}_{2} [GeV]");

            // --------------------
            // -- stop gen level
            // --------------------
            // Reco
            //my_histos["h_StopRecoMT2_"+cutVar.first]->Fill( StopRecoMT2, weight );    
            //my_histos["h_Stop1RecoMass_"+cutVar.first]->Fill( Stop1RecoMass, weight );
            //my_histos["h_Stop2RecoMass_"+cutVar.first]->Fill( Stop2RecoMass, weight );
            //my_histos["h_Nlino1RecoMass_"+cutVar.first]->Fill( Nlino1RecoMass, weight );
            //my_histos["h_Nlino2RecoMass_"+cutVar.first]->Fill( Nlino2RecoMass, weight );
            //my_histos["h_Single1RecoMass_"+cutVar.first]->Fill( Single1RecoMass, weight );
            //my_histos["h_Single2RecoMass_"+cutVar.first]->Fill( Single2RecoMass, weight );
            // Gen
            //my_histos["h_StopGenMT2_"+cutVar.first]->Fill( StopGenMT2, weight );
            //my_histos["h_Stop1GenMass_"+cutVar.first]->Fill( Stop1GenMass, weight );
            //my_histos["h_Stop2GenMass_"+cutVar.first]->Fill( Stop2GenMass, weight );
            //my_histos["h_Nlino1GenMass_"+cutVar.first]->Fill( Nlino1GenMass, weight );
            //my_histos["h_Nlino2GenMass_"+cutVar.first]->Fill( Nlino2GenMass, weight );
            //my_histos["h_Single1GenMass_"+cutVar.first]->Fill( Single1GenMass, weight );
            //my_histos["h_Single2RecoMass_"+cutVar.first]->Fill( Single2GenMass, weight );

            }
        }

        // -------------------------------
        // -- Cut flow absolute numbers
        // -------------------------------
        //if(true) my_histos["h_cutFlow_absolute"]->Fill(0.5, weight);
        //if(true && pass_general) my_histos["h_cutFlow_absolute"]->Fill(1.5, weight);  
        //if(true && pass_general && pass_0l) my_histos["h_cutFlow_absolute"]->Fill(2.5, weight);
        //if(true && pass_general && pass_0l && pass_HT500) my_histos["h_cutFlow_absolute"]->Fill(3.5, weight);
        //if(true && pass_general && pass_0l && pass_HT500 && pass_ge2b) my_histos["h_cutFlow_absolute"]->Fill(4.5, weight);
        //if(true && pass_general && pass_0l && pass_HT500 && pass_ge2b && pass_ge6j) my_histos["h_cutFlow_absolute"]->Fill(5.5, weight);
        //if(true && pass_general && pass_0l && pass_HT500 && pass_ge2b && pass_ge6j && pass_ge2t) my_histos["h_cutFlow_absolute"]->Fill(6.5, weight);
        //if(true && pass_general && pass_0l && pass_HT500 && pass_ge2b && pass_ge6j && pass_ge2t && pass_ge1dRbjets) my_histos["h_cutFlow_absolute"]->Fill(7.5, weight);        
    
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
