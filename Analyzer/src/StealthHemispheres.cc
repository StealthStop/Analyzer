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
        // without any rank
        my_histos.emplace( "h_stop1Mass_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop1Mass_"+cutVar.first).c_str(), ("h_stop1Mass_"+cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_stop1Eta_"+cutVar.first,  std::make_shared<TH1D> ( ("h_stop1Eta_"+cutVar.first).c_str(),  ("h_stop1Eta_"+cutVar.first).c_str(), 100, -6, 6 ) );
        my_histos.emplace( "h_stop1Phi_"+cutVar.first,  std::make_shared<TH1D> ( ("h_stop1Phi_"+cutVar.first).c_str(),  ("h_stop1Phi_"+cutVar.first).c_str(), 80, -4, 4 ) );
        my_histos.emplace( "h_stop1Pt_"+cutVar.first,   std::make_shared<TH1D> ( ("h_stop1Pt_"+cutVar.first).c_str(),   ("h_stop1Pt_"+cutVar.first).c_str(), 100, 0, 1000 ) );
        my_histos.emplace( "h_stop2Mass_"+cutVar.first, std::make_shared<TH1D> ( ("h_stop2Mass_"+cutVar.first).c_str(), ("h_stop2Mass_"+cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_stop2Eta_"+cutVar.first,  std::make_shared<TH1D> ( ("h_stop2Eta_"+cutVar.first).c_str(),  ("h_stop2Eta_"+cutVar.first).c_str(), 100, -6, 6 ) );
        my_histos.emplace( "h_stop2Phi_"+cutVar.first,  std::make_shared<TH1D> ( ("h_stop2Phi_"+cutVar.first).c_str(),  ("h_stop2Phi_"+cutVar.first).c_str(), 80, -4, 4 ) );
        my_histos.emplace( "h_stop2Pt_"+cutVar.first,   std::make_shared<TH1D> ( ("h_stop2Pt_"+cutVar.first).c_str(),   ("h_stop2Pt_"+cutVar.first).c_str(), 100, 0, 1000 ) );
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
        // without any rank
        my_2d_histos.emplace( "h_Mass_stop1vsstop2_"+cutVar.first, std::make_shared<TH2D>( ("h_Mass_stop1vsstop2_"+cutVar.first).c_str(), ("h_Mass_stop1vsstop2_"+cutVar.first).c_str(), 500, 0, 1500, 500, 0, 1500 ) );
        my_2d_histos.emplace( "h_Eta_stop1vsstop2_"+cutVar.first,  std::make_shared<TH2D>( ("h_Eta_stop1vsstop2_"+cutVar.first).c_str(),  ("h_Eta_stop1vsstop2_"+cutVar.first).c_str(),  100, -6, 6, 100, -6, 6 ) );
        my_2d_histos.emplace( "h_Phi_stop1vsstop2_"+cutVar.first,  std::make_shared<TH2D>( ("h_Phi_stop1vsstop2_"+cutVar.first).c_str(),  ("h_Phi_stop1vsstop2_"+cutVar.first).c_str(),  80, -4, 4, 80, -4, 4 ) );
        my_2d_histos.emplace( "h_Pt_stop1vsstop2_"+cutVar.first,   std::make_shared<TH2D>( ("h_Pt_stop1vsstop2_"+cutVar.first).c_str(),   ("h_Pt_stop1vsstop2_"+cutVar.first).c_str(),   100, 0, 1000, 100, 0, 1000) );    
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
        my_2d_histos.emplace( "h_Mass_NJetsVSstop1_"+cutVar.first, std::make_shared<TH2D> ( ("h_Mass_NJetsVSstop1_"+cutVar.first).c_str(), ("h_Mass_NJetsVSstop1_"+cutVar.first).c_str(), 20, 0, 20, 500, 0, 1500 ) );
        my_2d_histos.emplace( "h_Mass_NJetsVSstop2_"+cutVar.first, std::make_shared<TH2D> ( ("h_Mass_NJetsVSstop2_"+cutVar.first).c_str(), ("h_Mass_NJetsVSstop2_"+cutVar.first).c_str(), 20, 0, 20, 500, 0, 1500 ) );
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
        // without any rank
        my_2d_histos.emplace( "h_stop1_MassVsPt_"+cutVar.first, std::make_shared<TH2D>( ("h_stop1_MassVsPt_"+cutVar.first).c_str(), ("h_stop1_MassVsPt_"+cutVar.first).c_str(), 500, 0, 1500, 100, 0, 1000 ) );
        my_2d_histos.emplace( "h_stop2_MassVsPt_"+cutVar.first, std::make_shared<TH2D>( ("h_stop2_MassVsPt_"+cutVar.first).c_str(), ("h_stop2_MassVsPt_"+cutVar.first).c_str(), 500, 0, 1500, 100, 0, 1000 ) );
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

        // General variables
        const auto& runtype               = tr.getVar<std::string>("runtype");
        // Old Baseline selection
        const auto& passBaseline0l_old    = tr.getVar<bool>("passBaseline0l_old");
        const auto& NGoodJets_pt45        = tr.getVar<int>("NGoodJets_pt45");
        const bool pass_ge7j_pt45         = NGoodJets_pt45 >= 7;
        const bool pass_6j_pt45           = NGoodJets_pt45 == 6;
        const bool pass_7j_pt45           = NGoodJets_pt45 == 7;
        const bool pass_8j_pt45           = NGoodJets_pt45 == 8;
        const bool pass_9j_pt45           = NGoodJets_pt45 == 9;
        const bool pass_10j_pt45          = NGoodJets_pt45 == 10;
        const bool pass_11j_pt45          = NGoodJets_pt45 == 11;
        const bool pass_ge12j_pt45        = NGoodJets_pt45 >= 12;
        // New Baseline selection 
        const auto& passBaseline0l_pre    = tr.getVar<bool>("passBaseline0l_pre");
        const auto& passBaseline0l_good   = tr.getVar<bool>("passBaseline0l_good");
        const auto& NGoodJets_pt30        = tr.getVar<int>("NGoodJets_pt30");
        const bool pass_7j_pt30           = NGoodJets_pt30 == 7;
        const bool pass_8j_pt30           = NGoodJets_pt30 == 8;
        const bool pass_9j_pt30           = NGoodJets_pt30 == 9;
        const bool pass_10j_pt30          = NGoodJets_pt30 == 10;
        const bool pass_11j_pt30          = NGoodJets_pt30 == 11;
        const bool pass_ge12j_pt30        = NGoodJets_pt30 >= 12;

        // -------------------------------------
        // -- Make Stop Hemispheres variables
        // -------------------------------------
        const auto& Stop1                      = tr.getVar<TLorentzVector>("Stop1_OldSeed");
        const auto& Stop2                      = tr.getVar<TLorentzVector>("Stop2_OldSeed");
        const auto& Stop1_PtRank               = tr.getVar<TLorentzVector>("Stop1_PtRank_OldSeed");
        const auto& Stop2_PtRank               = tr.getVar<TLorentzVector>("Stop2_PtRank_OldSeed");
        const auto& Stop1_MassRank             = tr.getVar<TLorentzVector>("Stop1_MassRank_OldSeed");
        const auto& Stop2_MassRank             = tr.getVar<TLorentzVector>("Stop2_MassRank_OldSeed");
        const auto& Stop1_ScalarPtRank         = tr.getVar<TLorentzVector>("Stop1_ScalarPtRank_OldSeed");
        const auto& Stop2_ScalarPtRank         = tr.getVar<TLorentzVector>("Stop2_ScalarPtRank_OldSeed");
        const auto& Stop1ScalarPt_ScalarPtRank = tr.getVar<double>("Stop1ScalarPt_ScalarPtRank_OldSeed");
        const auto& Stop2ScalarPt_ScalarPtRank = tr.getVar<double>("Stop2ScalarPt_ScalarPtRank_OldSeed");
        const auto& MT2                        = tr.getVar<double>("MT2_OldSeed"); 
        const auto& dR_Stop1Stop2              = tr.getVar<double>("dR_Stop1Stop2_OldSeed");
        const auto& dPhi_Stop1Stop2            = tr.getVar<double>("dPhi_Stop1Stop2_OldSeed");
        const auto& difference_stopMasses      = tr.getVar<double>("difference_stopMasses_OldSeed");
        const auto& average_stopMasses         = tr.getVar<double>("average_stopMasses_OldSeed");
        const auto& relativeDiff_stopMasses    = tr.getVar<double>("relativeDiff_stopMasses_OldSeed");

        double stop1Mass              = 0.0, stop1Eta              = 0.0, stop1Phi              = 0.0, stop1Pt              = 0.0;
        double stop2Mass              = 0.0, stop2Eta              = 0.0, stop2Phi              = 0.0, stop2Pt              = 0.0;
        double stop1Mass_PtRank       = 0.0, stop1Eta_PtRank       = 0.0, stop1Phi_PtRank       = 0.0, stop1Pt_PtRank       = 0.0;
        double stop2Mass_PtRank       = 0.0, stop2Eta_PtRank       = 0.0, stop2Phi_PtRank       = 0.0, stop2Pt_PtRank       = 0.0;
        double stop1Mass_MassRank     = 0.0, stop1Eta_MassRank     = 0.0, stop1Phi_MassRank     = 0.0, stop1Pt_MassRank     = 0.0;
        double stop2Mass_MassRank     = 0.0, stop2Eta_MassRank     = 0.0, stop2Phi_MassRank     = 0.0, stop2Pt_MassRank     = 0.0;
        double stop1Mass_ScalarPtRank = 0.0, stop1Eta_ScalarPtRank = 0.0, stop1Phi_ScalarPtRank = 0.0, stop1Pt_ScalarPtRank = 0.0;
        double stop2Mass_ScalarPtRank = 0.0, stop2Eta_ScalarPtRank = 0.0, stop2Phi_ScalarPtRank = 0.0, stop2Pt_ScalarPtRank = 0.0;        
        
        stop1Mass              = Stop1.M();
        stop1Eta               = Stop1.Eta();
        stop1Phi               = Stop1.Phi();
        stop1Pt                = Stop1.Pt();
        stop2Mass              = Stop2.M();
        stop2Eta               = Stop2.Eta();
        stop2Phi               = Stop2.Phi();
        stop2Pt                = Stop2.Pt();
    
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
            // Old baseline selection
            {"0l_HT500_ge6j_ge2b_ge2t_ge1dRbjets",  passBaseline0l_old                   },
            {"0l_HT500_ge7j_ge2b_ge2t_ge1dRbjets",  passBaseline0l_old && pass_ge7j_pt45 },
            {"0l_HT500_6j_ge2b_ge2t_ge1dRbjets",    passBaseline0l_old && pass_6j_pt45   },
            {"0l_HT500_7j_ge2b_ge2t_ge1dRbjets",    passBaseline0l_old && pass_7j_pt45   },
            {"0l_HT500_8j_ge2b_ge2t_ge1dRbjets",    passBaseline0l_old && pass_8j_pt45   },
            {"0l_HT500_9j_ge2b_ge2t_ge1dRbjets",    passBaseline0l_old && pass_9j_pt45   },
            {"0l_HT500_10j_ge2b_ge2t_ge1dRbjets",   passBaseline0l_old && pass_10j_pt45  },
            {"0l_HT500_11j_ge2b_ge2t_ge1dRbjets",   passBaseline0l_old && pass_11j_pt45  },
            {"0l_HT500_ge12j_ge2b_ge2t_ge1dRbjets", passBaseline0l_old && pass_ge12j_pt45},

            // New baseline selection
            {"0l_HT500_0NonIsoMuon_ge7j_ge2t_ge1dRbjets",  passBaseline0l_pre && passBaseline0l_good                   },
            {"0l_HT500_0NonIsoMuon_7j_ge2t_ge1dRbjets",    passBaseline0l_pre && passBaseline0l_good && pass_7j_pt30   },
            {"0l_HT500_0NonIsoMuon_8j_ge2t_ge1dRbjets",    passBaseline0l_pre && passBaseline0l_good && pass_8j_pt30   },
            {"0l_HT500_0NonIsoMuon_9j_ge2t_ge1dRbjets",    passBaseline0l_pre && passBaseline0l_good && pass_9j_pt30   },
            {"0l_HT500_0NonIsoMuon_10j_ge2t_ge1dRbjets",   passBaseline0l_pre && passBaseline0l_good && pass_10j_pt30  },
            {"0l_HT500_0NonIsoMuon_11j_ge2t_ge1dRbjets",   passBaseline0l_pre && passBaseline0l_good && pass_11j_pt30  },
            {"0l_HT500_0NonIsoMuon_ge12j_ge2t_ge1dRbjets", passBaseline0l_pre && passBaseline0l_good && pass_ge12j_pt30},
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
                // 1D - without any rank
                my_histos["h_stop1Mass_"+cutVar.first]->Fill( stop1Mass, weight );
                my_histos["h_stop1Eta_"+cutVar.first]->Fill( stop1Eta, weight );
                my_histos["h_stop1Phi_"+cutVar.first]->Fill( stop1Phi, weight );
                my_histos["h_stop1Pt_"+cutVar.first]->Fill( stop1Pt, weight );
                my_histos["h_stop2Mass_"+cutVar.first]->Fill( stop2Mass, weight );
                my_histos["h_stop2Eta_"+cutVar.first]->Fill( stop2Eta, weight );
                my_histos["h_stop2Phi_"+cutVar.first]->Fill( stop2Phi, weight );
                my_histos["h_stop2Pt_"+cutVar.first]->Fill( stop2Pt, weight );
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
                my_2d_histos["h_Mass_stop1vsstop2_"+cutVar.first]->Fill(stop1Mass, stop2Mass, weight);
                my_2d_histos["h_Mass_stop1vsstop2_"+cutVar.first]->GetXaxis()->SetTitle("M_{#tildet}_{1}");
                my_2d_histos["h_Mass_stop1vsstop2_"+cutVar.first]->GetYaxis()->SetTitle("M_{#tildet}_{2}");
                my_2d_histos["h_Eta_stop1vsstop2_"+cutVar.first]->Fill(stop1Eta, stop2Eta, weight);
                my_2d_histos["h_Eta_stop1vsstop2_"+cutVar.first]->GetXaxis()->SetTitle("#eta_{#tildet}_{1}");
                my_2d_histos["h_Eta_stop1vsstop2_"+cutVar.first]->GetYaxis()->SetTitle("#eta_{#tildet}_{2}");
                my_2d_histos["h_Phi_stop1vsstop2_"+cutVar.first]->Fill(stop1Phi, stop2Phi, weight);
                my_2d_histos["h_Phi_stop1vsstop2_"+cutVar.first]->GetXaxis()->SetTitle("#phi_{#tildet}_{1}");
                my_2d_histos["h_Phi_stop1vsstop2_"+cutVar.first]->GetYaxis()->SetTitle("#phi_{#tildet}_{2}");
                my_2d_histos["h_Pt_stop1vsstop2_"+cutVar.first]->Fill(stop1Pt, stop2Pt, weight);
                my_2d_histos["h_Pt_stop1vsstop2_"+cutVar.first]->GetXaxis()->SetTitle("pT_{#tildet}_{1}");
                my_2d_histos["h_Pt_stop1vsstop2_"+cutVar.first]->GetYaxis()->SetTitle("pT_{#tildet}_{2}");
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
                my_2d_histos["h_Mass_NJetsVSstop1_"+cutVar.first]->Fill( NGoodJets_pt45, stop1Mass, weight );
                my_2d_histos["h_Mass_NJetsVSstop1_"+cutVar.first]->GetXaxis()->SetTitle("N_{J}");
                my_2d_histos["h_Mass_NJetsVSstop1_"+cutVar.first]->GetYaxis()->SetTitle("M_{#tildet}_{1}");
                my_2d_histos["h_Mass_NJetsVSstop2_"+cutVar.first]->Fill( NGoodJets_pt45, stop2Mass, weight );
                my_2d_histos["h_Mass_NJetsVSstop2_"+cutVar.first]->GetXaxis()->SetTitle("N_{J}");
                my_2d_histos["h_Mass_NJetsVSstop2_"+cutVar.first]->GetYaxis()->SetTitle("M_{#tildet}_{2}");
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
                my_2d_histos["h_stop1_MassVsPt_"+cutVar.first]->Fill(stop1Mass, stop1Pt, weight);
                my_2d_histos["h_stop1_MassVsPt_"+cutVar.first]->GetXaxis()->SetTitle("M_{#tildet}_{1}");
                my_2d_histos["h_stop1_MassVsPt_"+cutVar.first]->GetYaxis()->SetTitle("pT_{#tildet}_{1}");
                my_2d_histos["h_stop2_MassVsPt_"+cutVar.first]->Fill(stop2Mass, stop2Pt, weight);
                my_2d_histos["h_stop2_MassVsPt_"+cutVar.first]->GetXaxis()->SetTitle("M_{#tildet}_{2}");
                my_2d_histos["h_stop2_MassVsPt_"+cutVar.first]->GetYaxis()->SetTitle("pT_{#tildet}_{2}");
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
