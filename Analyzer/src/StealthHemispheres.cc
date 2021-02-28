#define StealthHemispheres_cxx
#include "Analyzer/Analyzer/include/StealthHemispheres.h"
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

StealthHemispheres::StealthHemispheres() : inithisto(false) 
{
}

// -------------------
// -- Define histos
// -------------------
void StealthHemispheres::InitHistos(const std::map<std::string, bool>& cutmap) 
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    
    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ;
    
    for (const auto& cutVar : cutmap) 
    {  
        // ---------------------------
        // -- Make Stop Hemispheres  
        // ---------------------------
        // pt rank 
        my_histos.emplace( "h_stop1Mass_PtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1Mass_PtRank_"+cutVar.first).c_str(), ("h_stop1Mass_PtRank_"+cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_stop1Eta_PtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1Eta_PtRank_"+cutVar.first).c_str(), ("h_stop1Eta_PtRank_"+cutVar.first).c_str(), 100, -6, 6 ) );
        my_histos.emplace( "h_stop1Phi_PtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1Phi_PtRank_"+cutVar.first).c_str(), ("h_stop1Phi_PtRank_"+cutVar.first).c_str(), 80, -4, 4 ) );
        my_histos.emplace( "h_stop1Pt_PtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1Pt_PtRank_"+cutVar.first).c_str(), ("h_stop1Pt_PtRank_"+cutVar.first).c_str(), 100, 0, 1000 ) ); 
        my_histos.emplace( "h_stop2Mass_PtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2Mass_PtRank_"+cutVar.first).c_str(), ("h_stop2Mass_PtRank_"+cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_stop2Eta_PtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2Eta_PtRank_"+cutVar.first).c_str(), ("h_stop2Eta_PtRank_"+cutVar.first).c_str(), 100, -6, 6 ) );
        my_histos.emplace( "h_stop2Phi_PtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2Phi_PtRank_"+cutVar.first).c_str(), ("h_stop2Phi_PtRank_"+cutVar.first).c_str(), 80, -4, 4 ) );
        my_histos.emplace( "h_stop2Pt_PtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2Pt_PtRank_"+cutVar.first).c_str(), ("h_stop2Pt_PtRank_"+cutVar.first).c_str(), 100, 0, 1000 ) );
        // mass rank
        my_histos.emplace( "h_stop1Mass_MassRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1Mass_MassRank_"+cutVar.first).c_str(), ("h_stop1Mass_MassRank_"+cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_stop1Eta_MassRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1Eta_MassRank_"+cutVar.first).c_str(), ("h_stop1Eta_MassRank_"+cutVar.first).c_str(), 100, -6, 6 ) );
        my_histos.emplace( "h_stop1Phi_MassRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1Phi_MassRank_"+cutVar.first).c_str(), ("h_stop1Phi_MassRank_"+cutVar.first).c_str(), 80, -4, 4 ) );
        my_histos.emplace( "h_stop1Pt_MassRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1Pt_MassRank_"+cutVar.first).c_str(), ("h_stop1Pt_MassRank_"+cutVar.first).c_str(), 100, 0, 1000 ) );
        my_histos.emplace( "h_stop2Mass_MassRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2Mass_MassRank_"+cutVar.first).c_str(), ("h_stop2Mass_MassRank_"+cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_stop2Eta_MassRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2Eta_MassRank_"+cutVar.first).c_str(), ("h_stop2Eta_MassRank_"+cutVar.first).c_str(), 100, -6, 6 ) );
        my_histos.emplace( "h_stop2Phi_MassRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2Phi_MassRank_"+cutVar.first).c_str(), ("h_stop2Phi_MassRank_"+cutVar.first).c_str(), 80, -4, 4 ) );
        my_histos.emplace( "h_stop2Pt_MassRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2Pt_MassRank_"+cutVar.first).c_str(), ("h_stop2Pt_MassRank_"+cutVar.first).c_str(), 100, 0, 1000 ) );
        // scalarPt rank
        my_histos.emplace( "h_stop1Mass_ScalarPtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1Mass_ScalarPtRank_"+cutVar.first).c_str(), ("h_stop1Mass_ScalarPtRank_"+cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_stop1Eta_ScalarPtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1Eta_ScalarPtRank_"+cutVar.first).c_str(), ("h_stop1Eta_ScalarPtRank_"+cutVar.first).c_str(), 100, -6, 6 ) );
        my_histos.emplace( "h_stop1Phi_ScalarPtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1Phi_ScalarPtRank_"+cutVar.first).c_str(), ("h_stop1Phi_ScalarPtRank_"+cutVar.first).c_str(), 80, -4, 4 ) );
        my_histos.emplace( "h_stop1Pt_ScalarPtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1Pt_ScalarPtRank_"+cutVar.first).c_str(), ("h_stop1Pt_ScalarPtRank_"+cutVar.first).c_str(), 100, 0, 1000 ) );
        my_histos.emplace( "h_stop2Mass_ScalarPtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2Mass_ScalarPtRank_"+cutVar.first).c_str(), ("h_stop2Mass_ScalarPtRank_"+cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_stop2Eta_ScalarPtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2Eta_ScalarPtRank_"+cutVar.first).c_str(), ("h_stop2Eta_ScalarPtRank_"+cutVar.first).c_str(), 100, -6, 6 ) );
        my_histos.emplace( "h_stop2Phi_ScalarPtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2Phi_ScalarPtRank_"+cutVar.first).c_str(), ("h_stop2Phi_ScalarPtRank_"+cutVar.first).c_str(), 80, -4, 4 ) );
        my_histos.emplace( "h_stop2Pt_ScalarPtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2Pt_ScalarPtRank_"+cutVar.first).c_str(), ("h_stop2Pt_ScalarPtRank_"+cutVar.first).c_str(), 100, 0, 1000 ) );       
        
        my_histos.emplace( "h_stop1ScalarPt_ScalarPtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1ScalarPt_ScalarPtRank_"+cutVar.first).c_str(), ("h_stop1ScalarPt_ScalarPtRank_"+cutVar.first).c_str(), 100, 0, 1000 ) );
        my_histos.emplace( "h_stop2ScalarPt_ScalarPtRank_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2ScalarPt_ScalarPtRank_"+cutVar.first).c_str(), ("h_stop2ScalarPt_ScalarPtRank_"+cutVar.first).c_str(), 100, 0, 1000 ) );
        
        my_histos.emplace( "h_MT2_"+cutVar.first, std::make_shared<TH1D> ( ("h_MT2_"+cutVar.first).c_str(), ("h_MT2_"+cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_dR_stop1stop2_"+cutVar.first, std::make_shared<TH1D> ( ("h_dR_stop1stop2_"+cutVar.first).c_str(), ("h_dR_stop1stop2_"+cutVar.first).c_str(), 50, 0, 10 ) );
        my_histos.emplace( "h_dPhi_stop1stop2_"+cutVar.first, std::make_shared<TH1D> ( ("h_dPhi_stop1stop2_"+cutVar.first).c_str(), ("h_dPhi_stop1stop2_"+cutVar.first).c_str(), 50, 0, 10 ) );
        my_histos.emplace( "h_difference_stopMasses_"+cutVar.first, std::make_shared<TH1D> ( ("h_difference_stopMasses_"+cutVar.first).c_str(), ("h_difference_stopMasses_"+cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_average_stopMasses_"+cutVar.first, std::make_shared<TH1D> ( ("h_average_stopMasses_"+cutVar.first).c_str(), ("h_average_stopMasses_"+cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_relativeDiff_stopMasses_"+cutVar.first, std::make_shared<TH1D> ( ("h_relativeDiff_stopMasses_"+cutVar.first).c_str(), ("h_relativeDiff_stopMasses_"+cutVar.first).c_str(), 500, -1500, 1500) );

        // stop1VSstop2 Mass, Eta, Phi, Pt       
        // pt rank
        my_2d_histos.emplace( "h_Mass_stop1vsstop2_PtRank_"+cutVar.first, std::make_shared<TH2D>( ("h_Mass_stop1vsstop2_PtRank_"+cutVar.first).c_str(), ("h_Mass_stop1vsstop2_PtRank_"+cutVar.first).c_str(), 500, 0, 1500, 500, 0, 1500 ) );
        my_2d_histos.emplace( "h_Eta_stop1vsstop2_PtRank_"+cutVar.first, std::make_shared<TH2D>( ("h_Eta_stop1vsstop2_PtRank_"+cutVar.first).c_str(), ("h_Eta_stop1vsstop2_PtRank_"+cutVar.first).c_str(), 100, -6, 6, 100, -6, 6 ) );
        my_2d_histos.emplace( "h_Phi_stop1vsstop2_PtRank_"+cutVar.first, std::make_shared<TH2D>( ("h_Phi_stop1vsstop2_PtRank_"+cutVar.first).c_str(), ("h_Phi_stop1vsstop2_PtRank_"+cutVar.first).c_str(), 80, -4, 4, 80, -4, 4 ) );
        my_2d_histos.emplace( "h_Pt_stop1vsstop2_PtRank_"+cutVar.first, std::make_shared<TH2D>( ("h_Pt_stop1vsstop2_PtRank_"+cutVar.first).c_str(), ("h_Pt_stop1vsstop2_PtRank_"+cutVar.first).c_str(), 100, 0, 1000, 100, 0, 1000) );
        // mass rank
        my_2d_histos.emplace( "h_Mass_stop1vsstop2_MassRank_"+cutVar.first, std::make_shared<TH2D>( ("h_Mass_stop1vsstop2_MassRank_"+cutVar.first).c_str(), ("h_Mass_stop1vsstop2_MassRank_"+cutVar.first).c_str(), 500, 0, 1500, 500, 0, 1500 ) );
        my_2d_histos.emplace( "h_Eta_stop1vsstop2_MassRank_"+cutVar.first, std::make_shared<TH2D>( ("h_Eta_stop1vsstop2_MassRank_"+cutVar.first).c_str(), ("h_Eta_stop1vsstop2_MassRank_"+cutVar.first).c_str(), 100, -6, 6, 100, -6, 6 ) );
        my_2d_histos.emplace( "h_Phi_stop1vsstop2_MassRank_"+cutVar.first, std::make_shared<TH2D>( ("h_Phi_stop1vsstop2_MassRank_"+cutVar.first).c_str(), ("h_Phi_stop1vsstop2_MassRank_"+cutVar.first).c_str(), 80, -4, 4, 80, -4, 4 ) );
        my_2d_histos.emplace( "h_Pt_stop1vsstop2_MassRank_"+cutVar.first, std::make_shared<TH2D>( ("h_Pt_stop1vsstop2_MassRank_"+cutVar.first).c_str(), ("h_Pt_stop1vsstop2_MassRank_"+cutVar.first).c_str(), 100, 0, 1000, 100, 0, 1000) );
        // scalarPt rank
        my_2d_histos.emplace( "h_Mass_stop1vsstop2_ScalarPtRank_"+cutVar.first, std::make_shared<TH2D>( ("h_Mass_stop1vsstop2_ScalarPtRank_"+cutVar.first).c_str(), ("h_Mass_stop1vsstop2_ScalarPtRank_"+cutVar.first).c_str(), 500, 0, 1500, 500, 0, 1500 ) );    
        my_2d_histos.emplace( "h_Eta_stop1vsstop2_ScalarPtRank_"+cutVar.first, std::make_shared<TH2D>( ("h_Eta_stop1vsstop2_ScalarPtRank_"+cutVar.first).c_str(), ("h_Eta_stop1vsstop2_ScalarPtRank_"+cutVar.first).c_str(), 100, -6, 6, 100, -6, 6 ) );
        my_2d_histos.emplace( "h_Phi_stop1vsstop2_ScalarPtRank_"+cutVar.first, std::make_shared<TH2D>( ("h_Phi_stop1vsstop2_ScalarPtRank_"+cutVar.first).c_str(), ("h_Phi_stop1vsstop2_ScalarPtRank_"+cutVar.first).c_str(), 80, -4, 4, 80, -4, 4 ) );
        my_2d_histos.emplace( "h_Pt_stop1vsstop2_ScalarPtRank_"+cutVar.first, std::make_shared<TH2D>( ("h_Pt_stop1vsstop2_ScalarPtRank_"+cutVar.first).c_str(), ("h_Pt_stop1vsstop2_ScalarPtRank_"+cutVar.first).c_str(), 100, 0, 1000, 100, 0, 1000) ); 
        
        // NJetsVSstops
        my_2d_histos.emplace( "h_Mass_NJetsVSstop1_PtRank_"+cutVar.first, std::make_shared<TH2D> ( ("h_Mass_NJetsVSstop1_PtRank_"+cutVar.first).c_str(), ("h_Mass_NJetsVSstop1_PtRank_"+cutVar.first).c_str(), 20, 0, 20, 500, 0, 1500 ) );
        my_2d_histos.emplace( "h_Mass_NJetsVSstop2_PtRank_"+cutVar.first, std::make_shared<TH2D> ( ("h_Mass_NJetsVSstop2_PtRank_"+cutVar.first).c_str(), ("h_Mass_NJetsVSstop2_PtRank_"+cutVar.first).c_str(), 20, 0, 20, 500, 0, 1500 ) );
        my_2d_histos.emplace( "h_Mass_NJetsVSstop1_MassRank_"+cutVar.first, std::make_shared<TH2D> ( ("h_Mass_NJetsVSstop1_MassRank_"+cutVar.first).c_str(), ("h_Mass_NJetsVSstop1_MassRank_"+cutVar.first).c_str(), 20, 0, 20, 500, 0, 1500 ) );
        my_2d_histos.emplace( "h_Mass_NJetsVSstop2_MassRank_"+cutVar.first, std::make_shared<TH2D> ( ("h_Mass_NJetsVSstop2_MassRank_"+cutVar.first).c_str(), ("h_Mass_NJetsVSstop2_MassRank_"+cutVar.first).c_str(), 20, 0, 20, 500, 0, 1500 ) );       
        my_2d_histos.emplace( "h_Mass_NJetsVSstop1_ScalarPtRank_"+cutVar.first, std::make_shared<TH2D> ( ("h_Mass_NJetsVSstop1_ScalarPtRank_"+cutVar.first).c_str(), ("h_Mass_NJetsVSstop1_ScalarPtRank_"+cutVar.first).c_str(), 20, 0, 20, 500, 0, 1500 ) );
        my_2d_histos.emplace( "h_Mass_NJetsVSstop2_ScalarPtRank_"+cutVar.first, std::make_shared<TH2D> ( ("h_Mass_NJetsVSstop2_ScalarPtRank_"+cutVar.first).c_str(), ("h_Mass_NJetsVSstop2_ScalarPtRank_"+cutVar.first).c_str(), 20, 0, 20, 500, 0, 1500 ) );
        my_2d_histos.emplace( "h_stopMasses_diffVSavg_"+cutVar.first, std::make_shared<TH2D>( ("h_stopMasses_diffVSavg_"+cutVar.first).c_str(), ( "h_stopMasses_diffVSavg_"+cutVar.first).c_str(), 150, 0, 1500, 150, 0, 1500 ) ); // 150, -1500, 1500

        // stop Pt combinations 
        my_2d_histos.emplace( "h_ScalarPt_stop1vsstop2_ScalarPtRank_"+cutVar.first, std::make_shared<TH2D>( ("h_ScalarPt_stop1vsstop2_ScalarPtRank_"+cutVar.first).c_str(), ("h_ScalarPt_stop1vsstop2_ScalarPtRank_"+cutVar.first).c_str(), 100, 0, 1000, 100, 0, 1000) );
        my_2d_histos.emplace( "h_Pt1_PtRankVsScalarPtRank_"+cutVar.first, std::make_shared<TH2D>( ("h_Pt1_PtRankVsScalarPtRank_"+cutVar.first).c_str(), ("h_Pt1_PtRankVsScalarPtRank_"+cutVar.first).c_str(), 100, 0, 1000, 100, 0, 1000) );
        my_2d_histos.emplace( "h_Pt2_PtRankVsScalarPtRank_"+cutVar.first, std::make_shared<TH2D>( ("h_Pt2_PtRankVsScalarPtRank_"+cutVar.first).c_str(), ("h_Pt2_PtRankVsScalarPtRank_"+cutVar.first).c_str(), 100, 0, 1000, 100, 0, 1000) );
        my_2d_histos.emplace( "h_Pt1_PtRankVsScalarPt1_ScalarPtRank_"+cutVar.first, std::make_shared<TH2D>( ("h_Pt1_PtRankVsScalarPt1_ScalarPtRank_"+cutVar.first).c_str(), ("h_Pt1_PtRankVsScalarPt1_ScalarPtRank_"+cutVar.first).c_str(), 100, 0, 1000, 100, 0, 1000) );
        my_2d_histos.emplace( "h_Pt2_PtRankVsScalarPt2_ScalarPtRank"+cutVar.first, std::make_shared<TH2D>( ("h_Pt2_PtRankVsScalarPt2_ScalarPtRank_"+cutVar.first).c_str(), ("h_Pt2_PtRankVsScalarPt2_ScalarPtRank_"+cutVar.first).c_str(), 100, 0, 1000, 100, 0, 1000) );

        // stops MassVsPt
        // pt rank
        my_2d_histos.emplace( "h_stop1_MassVsPt_PtRank_"+cutVar.first, std::make_shared<TH2D>( ("h_stop1_MassVsPt_PtRank_"+cutVar.first).c_str(), ("h_stop1_MassVsPt_PtRank_"+cutVar.first).c_str(), 500, 0, 1500, 100, 0, 1000 ) );
        my_2d_histos.emplace( "h_stop2_MassVsPt_PtRank_"+cutVar.first, std::make_shared<TH2D>( ("h_stop2_MassVsPt_PtRank_"+cutVar.first).c_str(), ("h_stop2_MassVsPt_PtRank_"+cutVar.first).c_str(), 500, 0, 1500, 100, 0, 1000 ) );
        // mask rank
        my_2d_histos.emplace( "h_stop1_MassVsPt_MassRank_"+cutVar.first, std::make_shared<TH2D>( ("h_stop1_MassVsPt_MassRank_"+cutVar.first).c_str(), ("h_stop1_MassVsPt_MassRank_"+cutVar.first).c_str(), 500, 0, 1500, 100, 0, 1000 ) );
        my_2d_histos.emplace( "h_stop2_MassVsPt_MassRank_"+cutVar.first, std::make_shared<TH2D>( ("h_stop2_MassVsPt_MassRank_"+cutVar.first).c_str(), ("h_stop2_MassVsPt_MassRank_"+cutVar.first).c_str(), 500, 0, 1500, 100, 0, 1000 ) );
        // scalarPt rank
        my_2d_histos.emplace( "h_stop1_MassVsPt_ScalarPtRank_"+cutVar.first, std::make_shared<TH2D>( ("h_stop1_MassVsPt_ScalarPtRank_"+cutVar.first).c_str(), ("h_stop1_MassVsPt_ScalarPtRank_"+cutVar.first).c_str(), 500, 0, 1500, 100, 0, 1000 ) );
        my_2d_histos.emplace( "h_stop2_MassVsPt_ScalarPtRank_"+cutVar.first, std::make_shared<TH2D>( ("h_stop2_MassVsPt_ScalarPtRank_"+cutVar.first).c_str(), ("h_stop2_MassVsPt_ScalarPtRank_"+cutVar.first).c_str(), 500, 0, 1500, 100, 0, 1000 ) );
        my_2d_histos.emplace( "h_stop1_MassVsScalarPt_ScalarPtRank_"+cutVar.first, std::make_shared<TH2D>( ("h_stop1_MassVsScalarPt_ScalarPtRank_"+cutVar.first).c_str(), ("h_stop1_MassVsScalarPt_ScalarPtRank_"+cutVar.first).c_str(), 500, 0, 1500, 100, 0, 1000 ) );
        my_2d_histos.emplace( "h_stop2_MassVsScalarPt_ScalarPtRank_"+cutVar.first, std::make_shared<TH2D>( ("h_stop2_MassVsScalarPt_ScalarPtRank_"+cutVar.first).c_str(), ("h_stop2_MassVsScalarPt_ScalarPtRank_"+cutVar.first).c_str(), 500, 0, 1500, 100, 0, 1000 ) );

        // MT2
        my_2d_histos.emplace( "h_NJetsVsMT2_"+cutVar.first, std::make_shared<TH2D> ( ("h_NJetsVsMT2_"+cutVar.first).c_str(), ("h_NJetsVsMT2_"+cutVar.first).c_str(), 20, 0, 20, 500, 0, 1500 ) );
        my_2d_histos.emplace( "h_Mass_MT2vsstop1_PtRank_"+cutVar.first, std::make_shared<TH2D> ( ("h_Mass_MT2vsstop1_PtRank_"+cutVar.first).c_str(), ("h_Mass_MT2vsstop1_PtRank_"+cutVar.first).c_str(), 500, 0, 1500, 500, 0, 1500 ) );
        my_2d_histos.emplace( "h_Mass_MT2vsstop2_PtRank_"+cutVar.first, std::make_shared<TH2D> ( ("h_Mass_MT2vsstop2_PtRank_"+cutVar.first).c_str(), ("h_Mass_MT2vsstop2_PtRank_"+cutVar.first).c_str(), 500, 0, 1500, 500, 0, 1500 ) );
        my_2d_histos.emplace( "h_Mass_MT2vsstop1_MassRank_"+cutVar.first, std::make_shared<TH2D> ( ("h_Mass_MT2vsstop1_MassRank_"+cutVar.first).c_str(), ("h_Mass_MT2vsstop1_MassRank_"+cutVar.first).c_str(), 500, 0, 1500, 500, 0, 1500 ) );
        my_2d_histos.emplace( "h_Mass_MT2vsstop2_MassRank_"+cutVar.first, std::make_shared<TH2D> ( ("h_Mass_MT2vsstop2_MassRank_"+cutVar.first).c_str(), ("h_Mass_MT2vsstop2_MassRank_"+cutVar.first).c_str(), 500, 0, 1500, 500, 0, 1500 ) ); 
        my_2d_histos.emplace( "h_Mass_MT2vsstop1_ScalarPtRank_"+cutVar.first, std::make_shared<TH2D> ( ("h_Mass_MT2vsstop1_ScalarPtRank_"+cutVar.first).c_str(), ("h_Mass_MT2vsstop1_ScalarPtRank_"+cutVar.first).c_str(), 500, 0, 1500, 500, 0, 1500 ) );
        my_2d_histos.emplace( "h_Mass_MT2vsstop2_ScalarPtRank_"+cutVar.first, std::make_shared<TH2D> ( ("h_Mass_MT2vsstop2_ScalarPtRank_"+cutVar.first).c_str(), ("h_Mass_MT2vsstop2_ScalarPtRank_"+cutVar.first).c_str(), 500, 0, 1500, 500, 0, 1500 ) );
    
    }
}

// ---------------------------------------------
// -- Put everything you want to do per event 
// ---------------------------------------------
void StealthHemispheres::Loop(NTupleReader& tr, double, int maxevents, bool)
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
        const auto& NGoodLeptons    = tr.getVar<int>("NGoodLeptons");
        const auto& HT_trigger_pt45 = tr.getVar<double>("HT_trigger_pt45");
        const auto& NGoodJets_pt45  = tr.getVar<int>("NGoodJets_pt45");
        const auto& NGoodBJets_pt45 = tr.getVar<int>("NGoodBJets_pt45");
        const auto& dR_bjets        = tr.getVar<double>("dR_bjets");               
        const auto& ntops           = tr.getVar<int>("ntops");
        const auto& ntops_1jet      = tr.getVar<int>("ntops_1jet"); // merged
        const auto& ntops_2jet      = tr.getVar<int>("ntops_2jet");
        const auto& ntops_3jet      = tr.getVar<int>("ntops_3jet"); // resolved 
        const auto& passMadHT       = tr.getVar<bool>("passMadHT");
        const auto& passBaseline0l  = tr.getVar<bool>("passBaseline0l_Good");
        const auto& passMETFilters  = tr.getVar<bool>("passMETFilters");
        const bool pass_general     = JetID && passMETFilters && passMadHT;
        const bool pass_0l          = NGoodLeptons == 0;  
        const bool pass_HT500       = HT_trigger_pt45 > 500;
        const bool pass_ge2b        = NGoodBJets_pt45 >= 2;
        const bool pass_ge2t        = ntops >= 2;
        const bool pass_ge2t1j      = ntops >= 2 && ntops_3jet == 0 && ntops_2jet == 0;
        const bool pass_ge2t3j      = ntops >= 2 && ntops_1jet == 0 && ntops_2jet == 0;
        const bool pass_ge2t1j3j    = ntops >= 2 && ntops_1jet >= 1 && ntops_3jet >= 1 && ntops_2jet == 0;
        const bool pass_ge1dRbjets  = dR_bjets >= 1.0;       
    
        // -------------------------------------
        // -- Make Stop Hemispheres variables
        // -------------------------------------
        const auto& Stop1_PtRank               = tr.getVar<TLorentzVector>("Stop1_PtRank_OldSeed_maskedISR");
        const auto& Stop2_PtRank               = tr.getVar<TLorentzVector>("Stop2_PtRank_OldSeed_maskedISR");
        const auto& Stop1_MassRank             = tr.getVar<TLorentzVector>("Stop1_MassRank_OldSeed_maskedISR");
        const auto& Stop2_MassRank             = tr.getVar<TLorentzVector>("Stop2_MassRank_OldSeed_maskedISR");
        const auto& Stop1_ScalarPtRank         = tr.getVar<TLorentzVector>("Stop1_ScalarPtRank_OldSeed_maskedISR");
        const auto& Stop2_ScalarPtRank         = tr.getVar<TLorentzVector>("Stop2_ScalarPtRank_OldSeed_maskedISR");
        const auto& Stop1ScalarPt_ScalarPtRank = tr.getVar<double>("Stop1ScalarPt_ScalarPtRank_OldSeed_maskedISR");
        const auto& Stop2ScalarPt_ScalarPtRank = tr.getVar<double>("Stop2ScalarPt_ScalarPtRank_OldSeed_maskedISR");
        const auto& MT2                        = tr.getVar<double>("MT2_OldSeed_maskedISR"); 
        const auto& dR_Stop1Stop2              = tr.getVar<double>("dR_Stop1Stop2_OldSeed_maskedISR");
        const auto& dPhi_Stop1Stop2            = tr.getVar<double>("dPhi_Stop1Stop2_OldSeed_maskedISR");
        const auto& difference_stopMasses      = tr.getVar<double>("difference_stopMasses_OldSeed_maskedISR");
        const auto& average_stopMasses         = tr.getVar<double>("average_stopMasses_OldSeed_maskedISR");
        const auto& relativeDiff_stopMasses    = tr.getVar<double>("relativeDiff_stopMasses_OldSeed_maskedISR");

        double stop1Mass_PtRank       = 0.0, stop1Eta_PtRank       = 0.0, stop1Phi_PtRank       = 0.0, stop1Pt_PtRank       = 0.0;
        double stop2Mass_PtRank       = 0.0, stop2Eta_PtRank       = 0.0, stop2Phi_PtRank       = 0.0, stop2Pt_PtRank       = 0.0;
        double stop1Mass_MassRank     = 0.0, stop1Eta_MassRank     = 0.0, stop1Phi_MassRank     = 0.0, stop1Pt_MassRank     = 0.0;
        double stop2Mass_MassRank     = 0.0, stop2Eta_MassRank     = 0.0, stop2Phi_MassRank     = 0.0, stop2Pt_MassRank     = 0.0;
        double stop1Mass_ScalarPtRank = 0.0, stop1Eta_ScalarPtRank = 0.0, stop1Phi_ScalarPtRank = 0.0, stop1Pt_ScalarPtRank = 0.0;
        double stop2Mass_ScalarPtRank = 0.0, stop2Eta_ScalarPtRank = 0.0, stop2Phi_ScalarPtRank = 0.0, stop2Pt_ScalarPtRank = 0.0;        
        
        stop1Mass_PtRank       = Stop1_PtRank.M();
        stop1Eta_PtRank        = Stop1_PtRank.Eta();
        stop1Phi_PtRank        = Stop1_PtRank.Phi();       
        stop1Pt_PtRank         = Stop1_PtRank.Pt();
        stop2Mass_PtRank       = Stop2_PtRank.M();
        stop2Eta_PtRank        = Stop2_PtRank.Eta();
        stop2Phi_PtRank        = Stop2_PtRank.Phi(); 
        stop2Pt_PtRank         = Stop2_PtRank.Pt();

        stop1Mass_MassRank     = Stop1_MassRank.M();
        stop1Eta_MassRank      = Stop1_MassRank.Eta();
        stop1Phi_MassRank      = Stop1_MassRank.Phi(); 
        stop1Pt_MassRank       = Stop1_MassRank.Pt();
        stop2Mass_MassRank     = Stop2_MassRank.M();
        stop2Eta_MassRank      = Stop2_MassRank.Eta();
        stop2Phi_MassRank      = Stop2_MassRank.Phi();
        stop2Pt_MassRank       = Stop2_MassRank.Pt();

        stop1Mass_ScalarPtRank = Stop1_ScalarPtRank.M();
        stop1Eta_ScalarPtRank  = Stop1_ScalarPtRank.Eta();
        stop1Phi_ScalarPtRank  = Stop1_ScalarPtRank.Phi();
        stop1Pt_ScalarPtRank   = Stop1_ScalarPtRank.Pt();
        stop2Mass_ScalarPtRank = Stop2_ScalarPtRank.M();
        stop2Eta_ScalarPtRank  = Stop2_ScalarPtRank.Eta();
        stop2Phi_ScalarPtRank  = Stop2_ScalarPtRank.Phi();
        stop2Pt_ScalarPtRank   = Stop2_ScalarPtRank.Pt();        

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

        // -------------------------------------------------
        // -- Make cuts and fill histograms here & cutmap
        // -------------------------------------------------
        const std::map<std::string, bool>& cutmap
        {
            {"",                                       true},
           
            // baseline 
            {"0l_HT500_ge2b_ge6j_ge2t_ge1dRbjets",     passBaseline0l && pass_ge2t && pass_ge1dRbjets    },
            {"0l_HT500_ge2b_ge6j_ge2t1j_ge1dRbjets",   passBaseline0l && pass_ge2t1j && pass_ge1dRbjets  },
            {"0l_HT500_ge2b_ge6j_ge2t3j_ge1dRbjets",   passBaseline0l && pass_ge2t3j && pass_ge1dRbjets  },
            {"0l_HT500_ge2b_ge6j_ge2t1j3j_ge1dRbjets", passBaseline0l && pass_ge2t1j3j && pass_ge1dRbjets},
       
            // NJets cuts for stop hemispheres
            {"baseline_0l_Njet6",  pass_general && pass_0l && pass_HT500 && pass_ge2b && pass_ge2t && pass_ge1dRbjets && NGoodJets_pt45 == 6},
            {"baseline_0l_Njet7",  passBaseline0l && pass_ge2t && pass_ge1dRbjets && NGoodJets_pt45 == 7 },
            {"baseline_0l_Njet8",  passBaseline0l && pass_ge2t && pass_ge1dRbjets && NGoodJets_pt45 == 8 },
            {"baseline_0l_Njet9",  passBaseline0l && pass_ge2t && pass_ge1dRbjets && NGoodJets_pt45 == 9 },
            {"baseline_0l_Njet10", passBaseline0l && pass_ge2t && pass_ge1dRbjets && NGoodJets_pt45 == 10},
            {"baseline_0l_Njet11", passBaseline0l && pass_ge2t && pass_ge1dRbjets && NGoodJets_pt45 == 11},
            {"baseline_0l_Njet12", passBaseline0l && pass_ge2t && pass_ge1dRbjets && NGoodJets_pt45 == 12},
            {"baseline_0l_Njet13", passBaseline0l && pass_ge2t && pass_ge1dRbjets && NGoodJets_pt45 == 13},
            {"baseline_0l_Njet14", passBaseline0l && pass_ge2t && pass_ge1dRbjets && NGoodJets_pt45 == 14},
            {"baseline_0l_Njet15", passBaseline0l && pass_ge2t && pass_ge1dRbjets && NGoodJets_pt45 == 15},
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
                // ---------------------------
                // -- Make Stop Hemispheres  
                // ---------------------------
                // 1D - ptRank
                my_histos["h_stop1Mass_PtRank_"+cutVar.first]->Fill( stop1Mass_PtRank, weight );
                my_histos["h_stop1Eta_PtRank_"+cutVar.first]->Fill( stop1Eta_PtRank, weight );
                my_histos["h_stop1Phi_PtRank_"+cutVar.first]->Fill( stop1Phi_PtRank, weight );
                my_histos["h_stop1Pt_PtRank_"+cutVar.first]->Fill( stop1Pt_PtRank, weight );
                my_histos["h_stop2Mass_PtRank_"+cutVar.first]->Fill( stop2Mass_PtRank, weight );
                my_histos["h_stop2Eta_PtRank_"+cutVar.first]->Fill( stop2Eta_PtRank, weight );
                my_histos["h_stop2Phi_PtRank_"+cutVar.first]->Fill( stop2Phi_PtRank, weight );
                my_histos["h_stop2Pt_PtRank_"+cutVar.first]->Fill( stop2Pt_PtRank, weight );
                // 1D - massRank
                my_histos["h_stop1Mass_MassRank_"+cutVar.first]->Fill( stop1Mass_MassRank, weight );
                my_histos["h_stop1Eta_MassRank_"+cutVar.first]->Fill( stop1Eta_MassRank, weight );
                my_histos["h_stop1Phi_MassRank_"+cutVar.first]->Fill( stop1Phi_MassRank, weight );
                my_histos["h_stop1Pt_MassRank_"+cutVar.first]->Fill( stop1Pt_MassRank, weight );
                my_histos["h_stop2Mass_MassRank_"+cutVar.first]->Fill( stop2Mass_MassRank, weight );
                my_histos["h_stop2Eta_MassRank_"+cutVar.first]->Fill( stop2Eta_MassRank, weight );
                my_histos["h_stop2Phi_MassRank_"+cutVar.first]->Fill( stop2Phi_MassRank, weight );
                my_histos["h_stop2Pt_MassRank_"+cutVar.first]->Fill( stop2Pt_MassRank, weight );
                // 1D - scalarPtRank
                my_histos["h_stop1Mass_ScalarPtRank_"+cutVar.first]->Fill( stop1Mass_ScalarPtRank, weight );
                my_histos["h_stop1Eta_ScalarPtRank_"+cutVar.first]->Fill( stop1Eta_ScalarPtRank, weight );
                my_histos["h_stop1Phi_ScalarPtRank_"+cutVar.first]->Fill( stop1Phi_ScalarPtRank, weight );
                my_histos["h_stop1Pt_ScalarPtRank_"+cutVar.first]->Fill( stop1Pt_ScalarPtRank, weight );
                my_histos["h_stop2Mass_ScalarPtRank_"+cutVar.first]->Fill( stop2Mass_ScalarPtRank, weight );
                my_histos["h_stop2Eta_ScalarPtRank_"+cutVar.first]->Fill( stop2Eta_ScalarPtRank, weight );
                my_histos["h_stop2Phi_ScalarPtRank_"+cutVar.first]->Fill( stop2Phi_ScalarPtRank, weight );
                my_histos["h_stop2Pt_ScalarPtRank_"+cutVar.first]->Fill( stop2Pt_ScalarPtRank, weight );
                my_histos["h_stop1ScalarPt_ScalarPtRank_"+cutVar.first]->Fill( Stop1ScalarPt_ScalarPtRank, weight );
                my_histos["h_stop2ScalarPt_ScalarPtRank_"+cutVar.first]->Fill( Stop2ScalarPt_ScalarPtRank, weight );
                // 1D - others
                my_histos["h_MT2_"+cutVar.first]->Fill( MT2, weight );
                my_histos["h_dR_stop1stop2_"+cutVar.first]->Fill( dR_Stop1Stop2, weight );
                my_histos["h_dPhi_stop1stop2_"+cutVar.first]->Fill( dPhi_Stop1Stop2, weight );
                my_histos["h_difference_stopMasses_"+cutVar.first]->Fill( difference_stopMasses, weight );
                my_histos["h_average_stopMasses_"+cutVar.first]->Fill( average_stopMasses, weight );
                my_histos["h_relativeDiff_stopMasses_"+cutVar.first]->Fill( relativeDiff_stopMasses, weight );
                // 2D - stop1VSstop2 Mass, Eta, Phi, Pt
                my_2d_histos["h_Mass_stop1vsstop2_PtRank_"+cutVar.first]->Fill(stop1Mass_PtRank, stop2Mass_PtRank, weight);
                my_2d_histos["h_Mass_stop1vsstop2_PtRank_"+cutVar.first]->GetXaxis()->SetTitle("Pt Rank M_{#tildet}_{1}");
                my_2d_histos["h_Mass_stop1vsstop2_PtRank_"+cutVar.first]->GetYaxis()->SetTitle("Pt Rank M_{#tildet}_{2}");
                my_2d_histos["h_Eta_stop1vsstop2_PtRank_"+cutVar.first]->Fill(stop1Eta_PtRank, stop2Eta_PtRank, weight);
                my_2d_histos["h_Eta_stop1vsstop2_PtRank_"+cutVar.first]->GetXaxis()->SetTitle("Pt Rank #eta_{#tildet}_{1}");
                my_2d_histos["h_Eta_stop1vsstop2_PtRank_"+cutVar.first]->GetYaxis()->SetTitle("Pt Rank #eta_{#tildet}_{2}");
                my_2d_histos["h_Phi_stop1vsstop2_PtRank_"+cutVar.first]->Fill(stop1Phi_PtRank, stop2Phi_PtRank, weight);
                my_2d_histos["h_Phi_stop1vsstop2_PtRank_"+cutVar.first]->GetXaxis()->SetTitle("Pt Rank #phi_{#tildet}_{1}");
                my_2d_histos["h_Phi_stop1vsstop2_PtRank_"+cutVar.first]->GetYaxis()->SetTitle("Pt Rank #phi_{#tildet}_{2}");
                my_2d_histos["h_Pt_stop1vsstop2_PtRank_"+cutVar.first]->Fill(stop1Pt_PtRank, stop2Pt_PtRank, weight);
                my_2d_histos["h_Pt_stop1vsstop2_PtRank_"+cutVar.first]->GetXaxis()->SetTitle("Pt Rank pT_{#tildet}_{1}");
                my_2d_histos["h_Pt_stop1vsstop2_PtRank_"+cutVar.first]->GetYaxis()->SetTitle("Pt Rank pT_{#tildet}_{2}");
                my_2d_histos["h_Mass_stop1vsstop2_MassRank_"+cutVar.first]->Fill(stop1Mass_MassRank, stop2Mass_MassRank, weight);
                my_2d_histos["h_Mass_stop1vsstop2_MassRank_"+cutVar.first]->GetXaxis()->SetTitle("Mass Rank M_{#tildet}_{1}");
                my_2d_histos["h_Mass_stop1vsstop2_MassRank_"+cutVar.first]->GetYaxis()->SetTitle("Mass Rank M_{#tildet}_{2}");
                my_2d_histos["h_Eta_stop1vsstop2_MassRank_"+cutVar.first]->Fill(stop1Eta_MassRank, stop2Eta_MassRank, weight);
                my_2d_histos["h_Eta_stop1vsstop2_MassRank_"+cutVar.first]->GetXaxis()->SetTitle("Mass Rank #eta_{#tildet}_{1}");
                my_2d_histos["h_Eta_stop1vsstop2_MassRank_"+cutVar.first]->GetYaxis()->SetTitle("Mass Rank #eta_{#tildet}_{2}");
                my_2d_histos["h_Phi_stop1vsstop2_MassRank_"+cutVar.first]->Fill(stop1Phi_MassRank, stop2Phi_MassRank, weight);
                my_2d_histos["h_Phi_stop1vsstop2_MassRank_"+cutVar.first]->GetXaxis()->SetTitle("Mass Rank #phi_{#tildet}_{1}");
                my_2d_histos["h_Phi_stop1vsstop2_MassRank_"+cutVar.first]->GetYaxis()->SetTitle("Mass Rank #phi_{#tildet}_{2}");
                my_2d_histos["h_Pt_stop1vsstop2_MassRank_"+cutVar.first]->Fill(stop1Pt_MassRank, stop2Pt_MassRank, weight);
                my_2d_histos["h_Pt_stop1vsstop2_MassRank_"+cutVar.first]->GetXaxis()->SetTitle("Mass Rank pT_{#tildet}_{1}");
                my_2d_histos["h_Pt_stop1vsstop2_MassRank_"+cutVar.first]->GetYaxis()->SetTitle("Mass Rank pT_{#tildet}_{2}");
                my_2d_histos["h_Mass_stop1vsstop2_ScalarPtRank_"+cutVar.first]->Fill(stop1Mass_ScalarPtRank, stop2Mass_ScalarPtRank, weight);
                my_2d_histos["h_Mass_stop1vsstop2_ScalarPtRank_"+cutVar.first]->GetXaxis()->SetTitle("ScalarPt Rank M_{#tildet}_{1}");
                my_2d_histos["h_Mass_stop1vsstop2_ScalarPtRank_"+cutVar.first]->GetYaxis()->SetTitle("ScalarPt Rank M_{#tildet}_{2}");
                my_2d_histos["h_Eta_stop1vsstop2_ScalarPtRank_"+cutVar.first]->Fill(stop1Eta_ScalarPtRank, stop2Eta_ScalarPtRank, weight);
                my_2d_histos["h_Eta_stop1vsstop2_ScalarPtRank_"+cutVar.first]->GetXaxis()->SetTitle("ScalarPt Rank #eta_{#tildet}_{1}");
                my_2d_histos["h_Eta_stop1vsstop2_ScalarPtRank_"+cutVar.first]->GetYaxis()->SetTitle("ScalarPt Rank #eta_{#tildet}_{2}");
                my_2d_histos["h_Phi_stop1vsstop2_ScalarPtRank_"+cutVar.first]->Fill(stop1Phi_ScalarPtRank, stop2Phi_ScalarPtRank, weight);
                my_2d_histos["h_Phi_stop1vsstop2_ScalarPtRank_"+cutVar.first]->GetXaxis()->SetTitle("ScalarPt Rank #phi_{#tildet}_{1}");
                my_2d_histos["h_Phi_stop1vsstop2_ScalarPtRank_"+cutVar.first]->GetYaxis()->SetTitle("ScalarPt Rank #phi_{#tildet}_{2}");
                my_2d_histos["h_Pt_stop1vsstop2_ScalarPtRank_"+cutVar.first]->Fill(stop1Pt_ScalarPtRank, stop2Pt_ScalarPtRank, weight);
                my_2d_histos["h_Pt_stop1vsstop2_ScalarPtRank_"+cutVar.first]->GetXaxis()->SetTitle("ScalarPt Rank pT_{#tildet}_{1}");
                my_2d_histos["h_Pt_stop1vsstop2_ScalarPtRank_"+cutVar.first]->GetYaxis()->SetTitle("ScalarPt Rank pT_{#tildet}_{2}");
                // 2D - NJetsVSstops
                my_2d_histos["h_Mass_NJetsVSstop1_PtRank_"+cutVar.first]->Fill( NGoodJets_pt45, stop1Mass_PtRank, weight );
                my_2d_histos["h_Mass_NJetsVSstop1_PtRank_"+cutVar.first]->GetXaxis()->SetTitle("N_{J}");
                my_2d_histos["h_Mass_NJetsVSstop1_PtRank_"+cutVar.first]->GetYaxis()->SetTitle("Pt Rank M_{#tildet}_{1}"); 
                my_2d_histos["h_Mass_NJetsVSstop2_PtRank_"+cutVar.first]->Fill( NGoodJets_pt45, stop2Mass_PtRank, weight );
                my_2d_histos["h_Mass_NJetsVSstop2_PtRank_"+cutVar.first]->GetXaxis()->SetTitle("N_{J}");
                my_2d_histos["h_Mass_NJetsVSstop2_PtRank_"+cutVar.first]->GetYaxis()->SetTitle("Pt Rank M_{#tildet}_{2}");               
                my_2d_histos["h_Mass_NJetsVSstop1_MassRank_"+cutVar.first]->Fill( NGoodJets_pt45, stop1Mass_MassRank, weight );
                my_2d_histos["h_Mass_NJetsVSstop1_MassRank_"+cutVar.first]->GetXaxis()->SetTitle("N_{J}");
                my_2d_histos["h_Mass_NJetsVSstop1_MassRank_"+cutVar.first]->GetYaxis()->SetTitle("Mass Rank M_{#tildet}_{1}");
                my_2d_histos["h_Mass_NJetsVSstop2_MassRank_"+cutVar.first]->Fill( NGoodJets_pt45, stop2Mass_MassRank, weight );
                my_2d_histos["h_Mass_NJetsVSstop2_MassRank_"+cutVar.first]->GetXaxis()->SetTitle("N_{J}");
                my_2d_histos["h_Mass_NJetsVSstop2_MassRank_"+cutVar.first]->GetYaxis()->SetTitle("Mass Rank M_{#tildet}_{2}");
                my_2d_histos["h_Mass_NJetsVSstop1_ScalarPtRank_"+cutVar.first]->Fill( NGoodJets_pt45, stop1Mass_ScalarPtRank, weight );
                my_2d_histos["h_Mass_NJetsVSstop1_ScalarPtRank_"+cutVar.first]->GetXaxis()->SetTitle("N_{J}");
                my_2d_histos["h_Mass_NJetsVSstop1_ScalarPtRank_"+cutVar.first]->GetYaxis()->SetTitle("ScalarPt Rank M_{#tildet}_{1}");
                my_2d_histos["h_Mass_NJetsVSstop2_ScalarPtRank_"+cutVar.first]->Fill( NGoodJets_pt45, stop2Mass_ScalarPtRank, weight );
                my_2d_histos["h_Mass_NJetsVSstop2_ScalarPtRank_"+cutVar.first]->GetXaxis()->SetTitle("N_{J}");
                my_2d_histos["h_Mass_NJetsVSstop2_ScalarPtRank_"+cutVar.first]->GetYaxis()->SetTitle("ScalarPt Rank M_{#tildet}_{2}");
                my_2d_histos["h_stopMasses_diffVSavg_"+cutVar.first]->Fill( difference_stopMasses, average_stopMasses, weight);
                my_2d_histos["h_stopMasses_diffVSavg_"+cutVar.first]->GetXaxis()->SetTitle("difference");
                my_2d_histos["h_stopMasses_diffVSavg_"+cutVar.first]->GetYaxis()->SetTitle("average");
                // 2D - stop All Pt combinations
                my_2d_histos["h_ScalarPt_stop1vsstop2_ScalarPtRank_"+cutVar.first]->Fill( Stop1ScalarPt_ScalarPtRank, Stop2ScalarPt_ScalarPtRank, weight );
                my_2d_histos["h_ScalarPt_stop1vsstop2_ScalarPtRank_"+cutVar.first]->GetXaxis()->SetTitle("ScalarPt Rank Scalar pT_{#tildet}_{1}");
                my_2d_histos["h_ScalarPt_stop1vsstop2_ScalarPtRank_"+cutVar.first]->GetYaxis()->SetTitle("ScalarPt Rank Scalar pT_{#tildet}_{2}");
                my_2d_histos["h_Pt1_PtRankVsScalarPtRank_"+cutVar.first]->Fill( stop1Pt_PtRank, stop1Pt_ScalarPtRank, weight );
                my_2d_histos["h_Pt1_PtRankVsScalarPtRank_"+cutVar.first]->GetXaxis()->SetTitle("Pt Rank pT_{#tildet}_{1}");
                my_2d_histos["h_Pt1_PtRankVsScalarPtRank_"+cutVar.first]->GetYaxis()->SetTitle("ScalarPt Rank pT_{#tildet}_{1}");
                my_2d_histos["h_Pt2_PtRankVsScalarPtRank_"+cutVar.first]->Fill( stop2Pt_PtRank, stop2Pt_ScalarPtRank, weight );
                my_2d_histos["h_Pt2_PtRankVsScalarPtRank_"+cutVar.first]->GetXaxis()->SetTitle("Pt Rank pT_{#tildet}_{2}");
                my_2d_histos["h_Pt2_PtRankVsScalarPtRank_"+cutVar.first]->GetYaxis()->SetTitle("ScalarPt Rank pT_{#tildet}_{2}");
                my_2d_histos["h_Pt1_PtRankVsScalarPt1_ScalarPtRank_"+cutVar.first]->Fill( stop1Pt_PtRank, Stop1ScalarPt_ScalarPtRank, weight );
                my_2d_histos["h_Pt1_PtRankVsScalarPt1_ScalarPtRank_"+cutVar.first]->GetXaxis()->SetTitle("Pt Rank pT_{#tildet}_{1}");
                my_2d_histos["h_Pt1_PtRankVsScalarPt1_ScalarPtRank_"+cutVar.first]->GetYaxis()->SetTitle("ScalarPt Rank Scalar pT_{#tildet}_{1}");
                my_2d_histos["h_Pt2_PtRankVsScalarPt2_ScalarPtRank"+cutVar.first]->Fill( stop2Pt_PtRank, Stop2ScalarPt_ScalarPtRank, weight );
                my_2d_histos["h_Pt2_PtRankVsScalarPt2_ScalarPtRank"+cutVar.first]->GetXaxis()->SetTitle("Pt Rank pT_{#tildet}_{2}");
                my_2d_histos["h_Pt2_PtRankVsScalarPt2_ScalarPtRank"+cutVar.first]->GetYaxis()->SetTitle("ScalarPt Rank Scalar pT_{#tildet}_{2}");
                // 2D - stops MassVsPt
                my_2d_histos["h_stop1_MassVsPt_PtRank_"+cutVar.first]->Fill(stop1Mass_PtRank, stop1Pt_PtRank, weight);
                my_2d_histos["h_stop1_MassVsPt_PtRank_"+cutVar.first]->GetXaxis()->SetTitle("Pt Rank M_{#tildet}_{1}");
                my_2d_histos["h_stop1_MassVsPt_PtRank_"+cutVar.first]->GetYaxis()->SetTitle("Pt Rank pT_{#tildet}_{1}");
                my_2d_histos["h_stop2_MassVsPt_PtRank_"+cutVar.first]->Fill(stop2Mass_PtRank, stop2Pt_PtRank, weight);
                my_2d_histos["h_stop2_MassVsPt_PtRank_"+cutVar.first]->GetXaxis()->SetTitle("Pt Rank M_{#tildet}_{2}");
                my_2d_histos["h_stop2_MassVsPt_PtRank_"+cutVar.first]->GetYaxis()->SetTitle("Pt Rank pT_{#tildet}_{2}");
                my_2d_histos["h_stop1_MassVsPt_MassRank_"+cutVar.first]->Fill(stop1Mass_MassRank, stop1Pt_MassRank, weight);
                my_2d_histos["h_stop1_MassVsPt_MassRank_"+cutVar.first]->GetXaxis()->SetTitle("Mass Rank M_{#tildet}_{1}");
                my_2d_histos["h_stop1_MassVsPt_MassRank_"+cutVar.first]->GetYaxis()->SetTitle("Mass Rank pT_{#tildet}_{1}");
                my_2d_histos["h_stop2_MassVsPt_MassRank_"+cutVar.first]->Fill(stop2Mass_MassRank, stop2Pt_MassRank, weight);
                my_2d_histos["h_stop2_MassVsPt_MassRank_"+cutVar.first]->GetXaxis()->SetTitle("Mass Rank M_{#tildet}_{2}");
                my_2d_histos["h_stop2_MassVsPt_MassRank_"+cutVar.first]->GetYaxis()->SetTitle("Mass Rank pT_{#tildet}_{2}");
                my_2d_histos["h_stop1_MassVsPt_ScalarPtRank_"+cutVar.first]->Fill(stop1Mass_ScalarPtRank, stop1Pt_ScalarPtRank, weight);
                my_2d_histos["h_stop1_MassVsPt_ScalarPtRank_"+cutVar.first]->GetXaxis()->SetTitle("ScalarPt Rank M_{#tildet}_{1}");
                my_2d_histos["h_stop1_MassVsPt_ScalarPtRank_"+cutVar.first]->GetYaxis()->SetTitle("ScalarPt Rank pT_{#tildet}_{1}");
                my_2d_histos["h_stop2_MassVsPt_ScalarPtRank_"+cutVar.first]->Fill(stop2Mass_ScalarPtRank, stop2Pt_ScalarPtRank, weight);
                my_2d_histos["h_stop2_MassVsPt_ScalarPtRank_"+cutVar.first]->GetXaxis()->SetTitle("ScalarPt Rank M_{#tildet}_{2}");
                my_2d_histos["h_stop2_MassVsPt_ScalarPtRank_"+cutVar.first]->GetYaxis()->SetTitle("ScalarPt Rank pT_{#tildet}_{2}");
                my_2d_histos["h_stop1_MassVsScalarPt_ScalarPtRank_"+cutVar.first]->Fill(stop1Mass_ScalarPtRank, Stop1ScalarPt_ScalarPtRank, weight);
                my_2d_histos["h_stop1_MassVsScalarPt_ScalarPtRank_"+cutVar.first]->GetXaxis()->SetTitle("ScalarPt Rank M_{#tildet}_{1}");
                my_2d_histos["h_stop1_MassVsScalarPt_ScalarPtRank_"+cutVar.first]->GetYaxis()->SetTitle("ScalarPt Rank Scalar pT_{#tildet}_{1}");
                my_2d_histos["h_stop2_MassVsScalarPt_ScalarPtRank_"+cutVar.first]->Fill(stop2Mass_ScalarPtRank, Stop1ScalarPt_ScalarPtRank, weight);
                my_2d_histos["h_stop2_MassVsScalarPt_ScalarPtRank_"+cutVar.first]->GetXaxis()->SetTitle("ScalarPt Rank M_{#tildet}_{2}");
                my_2d_histos["h_stop2_MassVsScalarPt_ScalarPtRank_"+cutVar.first]->GetYaxis()->SetTitle("ScalarPt Rank Scalar pT_{#tildet}_{2}");
                // 2D - MT2
                my_2d_histos["h_NJetsVsMT2_"+cutVar.first]->Fill( NGoodJets_pt45, MT2, weight );
                my_2d_histos["h_NJetsVsMT2_"+cutVar.first]->GetXaxis()->SetTitle("N_{J}");
                my_2d_histos["h_NJetsVsMT2_"+cutVar.first]->GetYaxis()->SetTitle("MT2");
                my_2d_histos["h_Mass_MT2vsstop1_PtRank_"+cutVar.first]->Fill( MT2, stop1Mass_PtRank, weight );
                my_2d_histos["h_Mass_MT2vsstop1_PtRank_"+cutVar.first]->GetXaxis()->SetTitle("MT2");
                my_2d_histos["h_Mass_MT2vsstop1_PtRank_"+cutVar.first]->GetYaxis()->SetTitle("Pt Rank M_{#tildet}_{1} [GeV]");
                my_2d_histos["h_Mass_MT2vsstop2_PtRank_"+cutVar.first]->Fill( MT2, stop2Mass_PtRank, weight );
                my_2d_histos["h_Mass_MT2vsstop2_PtRank_"+cutVar.first]->GetXaxis()->SetTitle("MT2");
                my_2d_histos["h_Mass_MT2vsstop2_PtRank_"+cutVar.first]->GetYaxis()->SetTitle("Pt Rank M_{#tildet}_{2} [GeV]");
                my_2d_histos["h_Mass_MT2vsstop1_MassRank_"+cutVar.first]->Fill( MT2, stop1Mass_MassRank, weight );
                my_2d_histos["h_Mass_MT2vsstop1_MassRank_"+cutVar.first]->GetXaxis()->SetTitle("MT2");
                my_2d_histos["h_Mass_MT2vsstop1_MassRank_"+cutVar.first]->GetYaxis()->SetTitle("Mass Rank M_{#tildet}_{1} [GeV]");
                my_2d_histos["h_Mass_MT2vsstop2_MassRank_"+cutVar.first]->Fill( MT2, stop2Mass_MassRank, weight );
                my_2d_histos["h_Mass_MT2vsstop2_MassRank_"+cutVar.first]->GetXaxis()->SetTitle("MT2");
                my_2d_histos["h_Mass_MT2vsstop2_MassRank_"+cutVar.first]->GetYaxis()->SetTitle("Mass Rank M_{#tildet}_{2} [GeV]");
                my_2d_histos["h_Mass_MT2vsstop1_ScalarPtRank_"+cutVar.first]->Fill( MT2, stop1Mass_ScalarPtRank, weight );
                my_2d_histos["h_Mass_MT2vsstop1_ScalarPtRank_"+cutVar.first]->GetXaxis()->SetTitle("MT2");
                my_2d_histos["h_Mass_MT2vsstop1_ScalarPtRank_"+cutVar.first]->GetYaxis()->SetTitle("ScalarPt Rank M_{#tildet}_{1} [GeV]");
                my_2d_histos["h_Mass_MT2vsstop2_ScalarPtRank_"+cutVar.first]->Fill( MT2, stop2Mass_ScalarPtRank, weight );
                my_2d_histos["h_Mass_MT2vsstop2_ScalarPtRank_"+cutVar.first]->GetXaxis()->SetTitle("MT2");
                my_2d_histos["h_Mass_MT2vsstop2_ScalarPtRank_"+cutVar.first]->GetYaxis()->SetTitle("ScalarPt Rank M_{#tildet}_{2} [GeV]");

            } 
        } 
    } 
}

void StealthHemispheres::WriteHistos(TFile* outfile)
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
