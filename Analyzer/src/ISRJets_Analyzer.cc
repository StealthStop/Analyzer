#define ISRJets_Analyzer_cxx
#include "Analyzer/Analyzer/include/ISRJets_Analyzer.h"
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

ISRJets_Analyzer::ISRJets_Analyzer() : inithisto(false) 
{
}

// -------------------
// -- Define histos
// -------------------
void ISRJets_Analyzer::InitHistos(const std::map<std::string, bool>& cutmap) 
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    
    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ;
    
    for (const auto& cutVar : cutmap) 
    {   
        // --------------------------------
        // -- ISR gen matching variables
        // --------------------------------
        // ISR gen matching deltaR
        my_histos.emplace( "h_GM_ISRmatching_allDR_"+cutVar.first, std::make_shared<TH1D> ( ("h_GM_ISRmatching_allDR_"+cutVar.first).c_str(), ("h_GM_ISRmatching_allDR_"+cutVar.first).c_str(), 1000, 0, 10 ) );
        my_histos.emplace( "h_GM_ISRmatching_bestDR_"+cutVar.first, std::make_shared<TH1D> ( ("h_GM_ISRmatching_bestDR_"+cutVar.first).c_str(), ("h_GM_ISRmatching_bestDR_"+cutVar.first).c_str(), 1000, 0, 10 ) );
        my_histos.emplace( "h_GM_ISRmatching_CutOnDR_DR_"+cutVar.first, std::make_shared<TH1D> ( ("h_GM_ISRmatching_CutOnDR_DR_"+cutVar.first).c_str(), ("h_GM_ISRmatching_CutOnDR_DR_"+cutVar.first).c_str(), 1000, 0, 10 ) );
        my_histos.emplace( "h_GM_ISRmatching_CutOnPtRatio_DR_"+cutVar.first, std::make_shared<TH1D> ( ("h_GM_ISRmatching_CutOnPtRatio_DR_"+cutVar.first).c_str(), ("h_GM_ISRmatching_CutOnPtRatio_DR_"+cutVar.first).c_str(), 1000, 0, 10 ) );        
        // ISR gen matching pt ratio
        my_histos.emplace( "h_GM_ISRmatching_allPtRatio_"+cutVar.first, std::make_shared<TH1D> ( ("h_GM_ISRmatching_allPtRatio_"+cutVar.first).c_str(), ("h_GM_ISRmatching_allPtRatio_"+cutVar.first).c_str(), 1000, 0, 10 ) );
        my_histos.emplace( "h_GM_ISRmatching_bestPtRatio_"+cutVar.first, std::make_shared<TH1D> ( ("h_GM_ISRmatching_bestPtRatio_"+cutVar.first).c_str(), ("h_GM_ISRmatching_bestPtRatio_"+cutVar.first).c_str(), 1000, 0, 10 ) );
        my_histos.emplace( "h_GM_ISRmatching_CutOnDR_PtRatio_"+cutVar.first, std::make_shared<TH1D> ( ("h_GM_ISRmatching_CutOnDR_PtRatio_"+cutVar.first).c_str(), ("h_GM_ISRmatching_CutOnDR_PtRatio_"+cutVar.first).c_str(), 1000, 0, 10 ) );
        my_histos.emplace( "h_GM_ISRmatching_CutOnPtRatio_PtRatio_"+cutVar.first, std::make_shared<TH1D> ( ("h_GM_ISRmatching_CutOnPtRatio_PtRatio_"+cutVar.first).c_str(), ("h_GM_ISRmatching_CutOnPtRatio_PtRatio_"+cutVar.first).c_str(), 1000, 0, 10 ) ); 

        // -----------------------------------------
        // -- Get the Gen ISR kinematic variables
        // -----------------------------------------
        my_histos.emplace( "h_nGenISR_"+cutVar.first, std::make_shared<TH1D> ( ("h_nGenISR_"+cutVar.first).c_str(), ("h_nGenISR_"+cutVar.first).c_str(), 20, 0, 20 ) );
        my_histos.emplace( "h_GenISR_Mass_"+cutVar.first, std::make_shared<TH1D> ( ("h_GenISR_Mass_"+cutVar.first).c_str(), ("h_GenISR_Mass_"+cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_GenISR_Pt_"+cutVar.first, std::make_shared<TH1D> ( ("h_GenISR_Pt_"+cutVar.first).c_str(), ("h_GenISR_Pt_"+cutVar.first).c_str(), 100, 0, 1000 ) );
        my_histos.emplace( "h_GenISR_Phi_"+cutVar.first, std::make_shared<TH1D> ( ("h_GenISR_Phi_"+cutVar.first).c_str(), ("h_GenISR_Phi_"+cutVar.first).c_str(), 80, -4, 4 ) );
        my_histos.emplace( "h_GenISR_Eta_"+cutVar.first, std::make_shared<TH1D> ( ("h_GenISR_Eta_"+cutVar.first).c_str(), ("h_GenISR_Eta_"+cutVar.first).c_str(), 100, -6, 6 ) );        

        my_histos.emplace( "h_nGenISR_qg_"+cutVar.first, std::make_shared<TH1D> ( ("h_nGenISR_qg_"+cutVar.first).c_str(), ("h_nGenISR_qg_"+cutVar.first).c_str(), 20, 0, 20 ) );
        my_histos.emplace( "h_GenISR_qg_Mass_"+cutVar.first, std::make_shared<TH1D> ( ("h_GenISR_qg_Mass_"+cutVar.first).c_str(), ("h_GenISR_Mass_"+cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_GenISR_qg_Pt_"+cutVar.first, std::make_shared<TH1D> ( ("h_GenISR_qg_Pt_"+cutVar.first).c_str(), ("h_GenISR_qg_Pt_"+cutVar.first).c_str(), 100, 0, 1000 ) );
        my_histos.emplace( "h_GenISR_qg_Phi_"+cutVar.first, std::make_shared<TH1D> ( ("h_GenISR_qg_Phi_"+cutVar.first).c_str(), ("h_GenISR_qg_Phi_"+cutVar.first).c_str(), 80, -4, 4 ) );
        my_histos.emplace( "h_GenISR_qg_Eta_"+cutVar.first, std::make_shared<TH1D> ( ("h_GenISR_qg_Eta_"+cutVar.first).c_str(), ("h_GenISR_qg_Eta_"+cutVar.first).c_str(), 100, -6, 6 ) ); 

        my_histos.emplace( "h_nGenISR_gq_"+cutVar.first, std::make_shared<TH1D> ( ("h_nGenISR_gq_"+cutVar.first).c_str(), ("h_nGenISR_gq_"+cutVar.first).c_str(), 20, 0, 20 ) );
        my_histos.emplace( "h_GenISR_gq_Mass_"+cutVar.first, std::make_shared<TH1D> ( ("h_GenISR_gq_Mass_"+cutVar.first).c_str(), ("h_GenISR_gq_Mass_"+cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_GenISR_gq_Pt_"+cutVar.first, std::make_shared<TH1D> ( ("h_GenISR_gq_Pt_"+cutVar.first).c_str(), ("h_GenISR_gq_Pt_"+cutVar.first).c_str(), 100, 0, 1000 ) );
        my_histos.emplace( "h_GenISR_gq_Phi_"+cutVar.first, std::make_shared<TH1D> ( ("h_GenISR_gq_Phi_"+cutVar.first).c_str(), ("h_GenISR_gq_Phi_"+cutVar.first).c_str(), 80, -4, 4 ) );
        my_histos.emplace( "h_GenISR_gq_Eta_"+cutVar.first, std::make_shared<TH1D> ( ("h_GenISR_gq_Eta_"+cutVar.first).c_str(), ("h_GenISR_gq_Eta_"+cutVar.first).c_str(), 100, -6, 6 ) );

        my_histos.emplace( "h_nGenISR_gg_"+cutVar.first, std::make_shared<TH1D> ( ("h_nGenISR_gg_"+cutVar.first).c_str(), ("h_nGenISR_gg_"+cutVar.first).c_str(), 20, 0, 20 ) );
        my_histos.emplace( "h_GenISR_gg_Mass_"+cutVar.first, std::make_shared<TH1D> ( ("h_GenISR_gg_Mass_"+cutVar.first).c_str(), ("h_GenISR_gg_Mass_"+cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_GenISR_gg_Pt_"+cutVar.first, std::make_shared<TH1D> ( ("h_GenISR_gg_Pt_"+cutVar.first).c_str(), ("h_GenISR_gg_Pt_"+cutVar.first).c_str(), 100, 0, 1000 ) );
        my_histos.emplace( "h_GenISR_gg_Phi_"+cutVar.first, std::make_shared<TH1D> ( ("h_GenISR_gg_Phi_"+cutVar.first).c_str(), ("h_GenISR_gg_Phi_"+cutVar.first).c_str(), 80, -4, 4 ) );
        my_histos.emplace( "h_GenISR_gg_Eta_"+cutVar.first, std::make_shared<TH1D> ( ("h_GenISR_gg_Eta_"+cutVar.first).c_str(), ("h_GenISR_gg_Eta_"+cutVar.first).c_str(), 100, -6, 6 ) );

        // ---------------------------------
        // -- Testing Jets_ISRMask filter
        // ---------------------------------
        my_histos.emplace( "h_Jets_ISRMask_Mass_"+cutVar.first, std::make_shared<TH1D> ( ("h_Jets_ISRMask_Mass_"+cutVar.first).c_str(), ("h_Jets_ISRMask_Mass_"+cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_Jets_ISRMask_Pt_"+cutVar.first, std::make_shared<TH1D> ( ("h_Jets_ISRMask_Pt_"+cutVar.first).c_str(), ("h_Jets_ISRMask_Pt_"+cutVar.first).c_str(), 100, 0, 1000 ) );
        my_histos.emplace( "h_Jets_ISRMask_Phi_"+cutVar.first, std::make_shared<TH1D> ( ("h_Jets_ISRMask_Phi_"+cutVar.first).c_str(), ("h_Jets_ISRMask_Phi_"+cutVar.first).c_str(), 80, -4, 4 ) );
        my_histos.emplace( "h_Jets_ISRMask_Eta_"+cutVar.first, std::make_shared<TH1D> ( ("h_Jets_ISRMask_Eta_"+cutVar.first).c_str(), ("h_Jets_ISRMask_Eta_"+cutVar.first).c_str(), 100, -6, 6 ) );

        // -----------------------------------------
        // -- ISR & NonISR & Other Jets variables
        // -----------------------------------------
        // --------------------------------------------------------------------
        // ISR Jets - matches cutting on dr - by using truth information of ISR
        // --------------------------------------------------------------------
        my_histos.emplace( "h_nISRJets_"+cutVar.first, std::make_shared<TH1D> ( ("h_nISRJets_"+cutVar.first).c_str(), ("h_nISRJets_"+cutVar.first).c_str(), 20, 0, 20 ) );
        my_histos.emplace( "h_ISRJets_Mass_"+cutVar.first, std::make_shared<TH1D> ( ("h_ISRJets_Mass_"+cutVar.first).c_str(), ("h_ISRJets_Mass_"+cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_ISRJets_Pt_"+cutVar.first, std::make_shared<TH1D> ( ("h_ISRJets_Pt_"+cutVar.first).c_str(), ("h_ISRJets_Pt_"+cutVar.first).c_str(), 100, 0, 1000 ) );
        my_histos.emplace( "h_ISRJets_Phi_"+cutVar.first, std::make_shared<TH1D> ( ("h_ISRJets_Phi_"+cutVar.first).c_str(), ("h_ISRJets_Phi_"+cutVar.first).c_str(), 80, -4, 4 ) );
        my_histos.emplace( "h_ISRJets_Eta_"+cutVar.first, std::make_shared<TH1D> ( ("h_ISRJets_Eta_"+cutVar.first).c_str(), ("h_ISRJets_Eta_"+cutVar.first).c_str(), 100, -6, 6 ) );
        my_histos.emplace( "h_ISRJets_axismajor_"+cutVar.first, std::make_shared<TH1D> ( ("h_ISRJets_axismajor_"+cutVar.first).c_str(), ("h_ISRJets_axismajor_"+cutVar.first).c_str(), 100, 0, 1 ) );
        my_histos.emplace( "h_ISRJets_axisminor_"+cutVar.first, std::make_shared<TH1D> ( ("h_ISRJets_axisminor_"+cutVar.first).c_str(), ("h_ISRJets_axisminor_"+cutVar.first).c_str(), 100, 0, 1 ) );
        my_histos.emplace( "h_ISRJets_chargedHadronEnergyFraction_"+cutVar.first, std::make_shared<TH1D> ( ("h_ISRJets_chargedHadronEnergyFraction_"+cutVar.first).c_str(), ("h_ISRJets_chargedHadronEnergyFraction_"+cutVar.first).c_str(), 100, 0, 1 ) );
        my_histos.emplace( "h_ISRJets_neutralHadronEnergyFraction_"+cutVar.first, std::make_shared<TH1D> ( ("h_ISRJets_neutralHadronEnergyFraction_"+cutVar.first).c_str(), ("h_ISRJets_neutralHadronEnergyFraction_"+cutVar.first).c_str(), 100, 0, 1 ) );
        my_histos.emplace( "h_ISRJets_ptD_"+cutVar.first, std::make_shared<TH1D> ( ("h_ISRJets_ptD_"+cutVar.first).c_str(), ("h_ISRJets_ptD_"+cutVar.first).c_str(), 100, 0, 1 ) );  
        my_histos.emplace( "h_ISRJets_chargedHadronMultiplicity_"+cutVar.first, std::make_shared<TH1D> ( ("h_ISRJets_chargedHadronMultiplicity_"+cutVar.first).c_str(), ("h_ISRJets_chargedHadronMultiplicity_"+cutVar.first).c_str(), 1000, 0, 100 ) );
        my_histos.emplace( "h_ISRJets_neutralHadronMultiplicity_"+cutVar.first, std::make_shared<TH1D> ( ("h_ISRJets_neutralHadronMultiplicity_"+cutVar.first).c_str(), ("h_ISRJets_neutralHadronMultiplicity_"+cutVar.first).c_str(), 1000, 0, 100 ) );
        my_histos.emplace( "h_ISRJets_chargedMultiplicity_"+cutVar.first, std::make_shared<TH1D> ( ("h_ISRJets_chargedMultiplicity_"+cutVar.first).c_str(), ("h_ISRJets_chargedMultiplicity_"+cutVar.first).c_str(), 1000, 0, 100 ) );
        my_histos.emplace( "h_ISRJets_neutralMultiplicity_"+cutVar.first, std::make_shared<TH1D> ( ("h_ISRJets_neutralMultiplicity_"+cutVar.first).c_str(), ("h_ISRJets_neutralMultiplicity_"+cutVar.first).c_str(), 1000, 0, 100 ) );
        my_histos.emplace( "h_ISRJets_qgLikelihood_"+cutVar.first, std::make_shared<TH1D> ( ("h_ISRJets_qgLikelihood_"+cutVar.first).c_str(), ("h_ISRJets_qgLikelihood_"+cutVar.first).c_str(), 100, 0, 1 ) );
        my_histos.emplace( "h_ISRJets_bJetTagDeepCSVtotb_"+cutVar.first, std::make_shared<TH1D> ( ("h_ISRJets_bJetTagDeepCSVtotb_"+cutVar.first).c_str(), ("h_ISRJets_bJetTagDeepCSVtotb_"+cutVar.first).c_str(), 100, 0, 1 ) );

        // --------------------------------
        // NonISR Jets - by using TreeMaker
        // --------------------------------
        my_histos.emplace( "h_nNonISRJets_"+cutVar.first, std::make_shared<TH1D> ( ("h_nNonISRJets_"+cutVar.first).c_str(), ("h_nNonISRJets_"+cutVar.first).c_str(), 20, 0, 20 ) );
        my_histos.emplace( "h_NonISRJets_Mass_"+cutVar.first, std::make_shared<TH1D> ( ("h_NonISRJets_Mass_"+cutVar.first).c_str(), ("h_NonISRJets_Mass_"+cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_NonISRJets_Pt_"+cutVar.first, std::make_shared<TH1D> ( ("h_NonISRJets_Pt_"+cutVar.first).c_str(), ("h_NonISRJets_Pt_"+cutVar.first).c_str(), 100, 0, 1000 ) );
        my_histos.emplace( "h_NonISRJets_Phi_"+cutVar.first, std::make_shared<TH1D> ( ("h_NonISRJets_Phi_"+cutVar.first).c_str(), ("h_NonISRJets_Phi_"+cutVar.first).c_str(), 80, -4, 4 ) );
        my_histos.emplace( "h_NonISRJets_Eta_"+cutVar.first, std::make_shared<TH1D> ( ("h_NonISRJets_Eta_"+cutVar.first).c_str(), ("h_NonISRJets_Eta_"+cutVar.first).c_str(), 100, -6, 6 ) );
        my_histos.emplace( "h_NonISRJets_axismajor_"+cutVar.first, std::make_shared<TH1D> ( ("h_NonISRJets_axismajor_"+cutVar.first).c_str(), ("h_NonISRJets_axismajor_"+cutVar.first).c_str(), 100, 0, 1 ) );
        my_histos.emplace( "h_NonISRJets_axisminor_"+cutVar.first, std::make_shared<TH1D> ( ("h_NonISRJets_axisminor_"+cutVar.first).c_str(), ("h_NonISRJets_axisminor_"+cutVar.first).c_str(), 100, 0, 1 ) );
        my_histos.emplace( "h_NonISRJets_ptD_"+cutVar.first, std::make_shared<TH1D> ( ("h_NonISRJets_ptD_"+cutVar.first).c_str(), ("h_NonISRJets_ptD_"+cutVar.first).c_str(), 100, 0, 1 ) ); 
        my_histos.emplace( "h_NonISRJets_chargedHadronEnergyFraction_"+cutVar.first, std::make_shared<TH1D> ( ("h_NonISRJets_chargedHadronEnergyFraction_"+cutVar.first).c_str(), ("h_NonISRJets_chargedHadronEnergyFraction_"+cutVar.first).c_str(), 100, 0, 1 ) );
        my_histos.emplace( "h_NonISRJets_neutralHadronEnergyFraction_"+cutVar.first, std::make_shared<TH1D> ( ("h_NonISRJets_neutralHadronEnergyFraction_"+cutVar.first).c_str(), ("h_NonISRJets_neutralHadronEnergyFraction_"+cutVar.first).c_str(), 100, 0, 1 ) );
        my_histos.emplace( "h_NonISRJets_chargedHadronMultiplicity_"+cutVar.first, std::make_shared<TH1D> ( ("h_NonISRJets_chargedHadronMultiplicity_"+cutVar.first).c_str(), ("h_NonISRJets_chargedHadronMultiplicity_"+cutVar.first).c_str(), 1000, 0, 100 ) );
        my_histos.emplace( "h_NonISRJets_neutralHadronMultiplicity_"+cutVar.first, std::make_shared<TH1D> ( ("h_NonISRJets_neutralHadronMultiplicity_"+cutVar.first).c_str(), ("h_NonISRJets_neutralHadronMultiplicity_"+cutVar.first).c_str(), 1000, 0, 100 ) );
        my_histos.emplace( "h_NonISRJets_chargedMultiplicity_"+cutVar.first, std::make_shared<TH1D> ( ("h_NonISRJets_chargedMultiplicity_"+cutVar.first).c_str(), ("h_NonISRJets_chargedMultiplicity_"+cutVar.first).c_str(), 1000, 0, 100 ) );
        my_histos.emplace( "h_NonISRJets_neutralMultiplicity_"+cutVar.first, std::make_shared<TH1D> ( ("h_NonISRJets_neutralMultiplicity_"+cutVar.first).c_str(), ("h_NonISRJets_neutralMultiplicity_"+cutVar.first).c_str(), 1000, 0, 100 ) );
        my_histos.emplace( "h_NonISRJets_qgLikelihood_"+cutVar.first, std::make_shared<TH1D> ( ("h_NonISRJets_qgLikelihood_"+cutVar.first).c_str(), ("h_NonISRJets_qgLikelihood_"+cutVar.first).c_str(), 100, 0, 1 ) );      
        my_histos.emplace( "h_NonISRJets_bJetTagDeepCSVtotb_"+cutVar.first, std::make_shared<TH1D> ( ("h_NonISRJets_bJetTagDeepCSVtotb_"+cutVar.first).c_str(), ("h_NonISRJets_bJetTagDeepCSVtotb_"+cutVar.first).c_str(), 100, 0, 1 ) );

        // ---------------------------------
        // Other Jets - Not ISR + not NonISR 
        // ---------------------------------
        my_histos.emplace( "h_nOtherJets_"+cutVar.first, std::make_shared<TH1D> ( ("h_nOtherJets_"+cutVar.first).c_str(), ("h_nOtherJets_"+cutVar.first).c_str(), 20, 0, 20 ) );
        my_histos.emplace( "h_OtherJets_Mass_"+cutVar.first, std::make_shared<TH1D> ( ("h_OtherJets_Mass_"+cutVar.first).c_str(), ("h_OtherJets_Mass_"+cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_OtherJets_Pt_"+cutVar.first, std::make_shared<TH1D> ( ("h_OtherJets_Pt_"+cutVar.first).c_str(), ("h_OtherJets_Pt_"+cutVar.first).c_str(), 100, 0, 1000 ) );
        my_histos.emplace( "h_OtherJets_Phi_"+cutVar.first, std::make_shared<TH1D> ( ("h_OtherJets_Phi_"+cutVar.first).c_str(), ("h_OtherJets_Phi_"+cutVar.first).c_str(), 80, -4, 4 ) );
        my_histos.emplace( "h_OtherJets_Eta_"+cutVar.first, std::make_shared<TH1D> ( ("h_OtherJets_Eta_"+cutVar.first).c_str(), ("h_OtherJets_Eta_"+cutVar.first).c_str(), 100, -6, 6 ) );
        my_histos.emplace( "h_OtherJets_axismajor_"+cutVar.first, std::make_shared<TH1D> ( ("h_OtherJets_axismajor_"+cutVar.first).c_str(), ("h_OtherJets_axismajor_"+cutVar.first).c_str(), 100, 0, 1 ) );
        my_histos.emplace( "h_OtherJets_axisminor_"+cutVar.first, std::make_shared<TH1D> ( ("h_OtherJets_axisminor_"+cutVar.first).c_str(), ("h_OtherJets_axisminor_"+cutVar.first).c_str(), 100, 0, 1 ) );
        my_histos.emplace( "h_OtherJets_ptD_"+cutVar.first, std::make_shared<TH1D> ( ("h_OtherJets_ptD_"+cutVar.first).c_str(), ("h_OtherJets_ptD_"+cutVar.first).c_str(), 100, 0, 1 ) ); 
        my_histos.emplace( "h_OtherJets_chargedHadronEnergyFraction_"+cutVar.first, std::make_shared<TH1D> ( ("h_OtherJets_chargedHadronEnergyFraction_"+cutVar.first).c_str(), ("h_OtherJets_chargedHadronEnergyFraction_"+cutVar.first).c_str(), 100, 0, 1 ) );
        my_histos.emplace( "h_OtherJets_neutralHadronEnergyFraction_"+cutVar.first, std::make_shared<TH1D> ( ("h_OtherJets_neutralHadronEnergyFraction_"+cutVar.first).c_str(), ("h_OtherJets_neutralHadronEnergyFraction_"+cutVar.first).c_str(), 100, 0, 1 ) );
        my_histos.emplace( "h_OtherJets_chargedHadronMultiplicity_"+cutVar.first, std::make_shared<TH1D> ( ("h_OtherJets_chargedHadronMultiplicity_"+cutVar.first).c_str(), ("h_OtherJets_chargedHadronMultiplicity_"+cutVar.first).c_str(), 1000, 0, 100 ) );
        my_histos.emplace( "h_OtherJets_neutralHadronMultiplicity_"+cutVar.first, std::make_shared<TH1D> ( ("h_OtherJets_neutralHadronMultiplicity_"+cutVar.first).c_str(), ("h_OtherJets_neutralHadronMultiplicity_"+cutVar.first).c_str(), 1000, 0, 100 ) );
        my_histos.emplace( "h_OtherJets_chargedMultiplicity_"+cutVar.first, std::make_shared<TH1D> ( ("h_OtherJets_chargedMultiplicity_"+cutVar.first).c_str(), ("h_OtherJets_chargedMultiplicity_"+cutVar.first).c_str(), 1000, 0, 100 ) );
        my_histos.emplace( "h_OtherJets_neutralMultiplicity_"+cutVar.first, std::make_shared<TH1D> ( ("h_OtherJets_neutralMultiplicity_"+cutVar.first).c_str(), ("h_OtherJets_neutralMultiplicity_"+cutVar.first).c_str(), 1000, 0, 100 ) ); 
        my_histos.emplace( "h_OtherJets_qgLikelihood_"+cutVar.first, std::make_shared<TH1D> ( ("h_OtherJets_qgLikelihood_"+cutVar.first).c_str(), ("h_OtherJets_qgLikelihood_"+cutVar.first).c_str(), 100, 0, 1 ) ); 
        my_histos.emplace( "h_OtherJets_bJetTagDeepCSVtotb_"+cutVar.first, std::make_shared<TH1D> ( ("h_OtherJets_bJetTagDeepCSVtotb_"+cutVar.first).c_str(), ("h_OtherJets_bJetTagDeepCSVtotb_"+cutVar.first).c_str(), 100, 0, 1 ) );

        // -----------------------------------
        // DeltaEta between GenISR and RecoISR
        // -----------------------------------
        my_histos.emplace( "h_dEta_RecoISR_GenISR_"+cutVar.first, std::make_shared<TH1D> ( ("h_dEta_RecoISR_GenISR_"+cutVar.first).c_str(), ("h_dEta_RecoISR_GenISR_"+cutVar.first).c_str(), 100, -2, 2 ) );   

        // ----------------------------------- 
        // DeltaR between ISR and Non ISR jets
        // -----------------------------------
        my_histos.emplace( "h_dR_ISR_NonISR_"+cutVar.first, std::make_shared<TH1D> ( ("h_dR_ISR_NonISR_"+cutVar.first).c_str(), ("h_dR_ISR_NonISR_"+cutVar.first).c_str(), 50, 0, 10 ) );
        my_histos.emplace( "h_dR_top_ISR_"+cutVar.first, std::make_shared<TH1D> ( ("h_dR_top_ISR_"+cutVar.first).c_str(), ("h_dR_top_ISR_"+cutVar.first).c_str(), 50, 0, 10 ) );
        my_histos.emplace( "h_dR_top_NonISR_"+cutVar.first, std::make_shared<TH1D> ( ("h_dR_top_NonISR_"+cutVar.first).c_str(), ("h_dR_top_NonISR_"+cutVar.first).c_str(), 50, 0, 10 ) );
        my_histos.emplace( "h_dR_bjet_ISR_"+cutVar.first, std::make_shared<TH1D> ( ("h_dR_bjet_ISR_"+cutVar.first).c_str(), ("h_dR_bjet_ISR_"+cutVar.first).c_str(), 50, 0, 10 ) );
        my_histos.emplace( "h_dR_bjet_NonISR_"+cutVar.first, std::make_shared<TH1D> ( ("h_dR_bjet_NonISR_"+cutVar.first).c_str(), ("h_dR_bjet_NonISR_"+cutVar.first).c_str(), 50, 0, 10 ) );

        // -------------------
        // -- 2D histograms
        // -------------------
        // ------------
        // Gen Matching
        // ------------
        my_2d_histos.emplace( "h_GM_all_PtRvsDR_"+cutVar.first, std::make_shared<TH2D>( ("h_GM_all_PtRvsDR_"+cutVar.first).c_str(), ("h_GM_all_PtRvsDR_"+cutVar.first).c_str(), 1000, 0, 10, 1000, 0, 10) );
        my_2d_histos.emplace( "h_GM_best_PtRvsDR_"+cutVar.first, std::make_shared<TH2D>( ("h_GM_best_PtRvsDR_"+cutVar.first).c_str(), ("h_GM_best_PtRvsDR_"+cutVar.first).c_str(), 1000, 0, 10, 1000, 0, 10) );
        my_2d_histos.emplace( "h_GM_CutOnDR_PtRvsDR_"+cutVar.first, std::make_shared<TH2D>( ("h_GM_CutOnDR_PtRvsDR_"+cutVar.first).c_str(), ("h_GM_CutOnDR_PtRvsDR_"+cutVar.first).c_str(), 1000, 0, 10, 1000, 0, 10) );
        my_2d_histos.emplace( "h_GM_CutOnPtRatio_PtRvsDR_"+cutVar.first, std::make_shared<TH2D>( ("h_GM_CutOnPtRatio_PtRvsDR_"+cutVar.first).c_str(), ("h_GM_CutOnPtRatio_PtRvsDR_"+cutVar.first).c_str(), 1000, 0, 10, 1000, 0, 10) );
        
        // --------------------------------------------
        // ISR Jets - by using truth information of ISR
        // --------------------------------------------
        my_2d_histos.emplace( "h_ISRJets_EtaVsPhi_"+cutVar.first, std::make_shared<TH2D> ( ("h_ISRJets_EtaVsPhi_"+cutVar.first).c_str(), ("h_ISRJets_EtaVsPhi_"+cutVar.first).c_str(), 100, -6, 6, 80, -4, 4 ) );
        my_2d_histos.emplace( "h_ISRJets_MassVsPt_"+cutVar.first, std::make_shared<TH2D>( ("h_ISRJets_MassVsPt_"+cutVar.first).c_str(), ("h_ISRJets_MassVsPt_"+cutVar.first).c_str(), 500, 0, 1500, 100, 0, 1000) );
        // Mass
        my_2d_histos.emplace( "h_ISRJets_MassVsChargedHadronEnergyFraction_"+cutVar.first, std::make_shared<TH2D>( ("h_ISRJets_MassVsChargedHadronEnergyFraction_"+cutVar.first).c_str(), ("h_ISRJets_MassVsChargedHadronEnergyFraction_"+cutVar.first).c_str(), 500, 0, 1500, 100, 0, 1) );
        my_2d_histos.emplace( "h_ISRJets_MassVsNeutralHadronEnergyFraction_"+cutVar.first, std::make_shared<TH2D>( ("h_ISRJets_MassVsNeutralHadronEnergyFraction_"+cutVar.first).c_str(), ("h_ISRJets_MassVsNeutralHadronEnergyFraction_"+cutVar.first).c_str(), 500, 0, 1500, 100, 0, 1) );
        my_2d_histos.emplace( "h_ISRJets_MassVsChargedHadronMultiplicity_"+cutVar.first, std::make_shared<TH2D>( ("h_ISRJets_MassVsChargedHadronMultiplicity_"+cutVar.first).c_str(), ("h_ISRJets_MassVsChargedHadronMultiplicity_"+cutVar.first).c_str(), 500, 0, 1500, 1000, 0, 100) );
        my_2d_histos.emplace( "h_ISRJets_MassVsNeutralHadronMultiplicity_"+cutVar.first, std::make_shared<TH2D>( ("h_ISRJets_MassVsNeutralHadronMultiplicity_"+cutVar.first).c_str(), ("h_ISRJets_MassVsNeutralHadronMultiplicity_"+cutVar.first).c_str(), 500, 0, 1500, 1000, 0, 100) );
        my_2d_histos.emplace( "h_ISRJets_MassVsChargedMultiplicity_"+cutVar.first, std::make_shared<TH2D>( ("h_ISRJets_MassVsChargedMultiplicity_"+cutVar.first).c_str(), ("h_ISRJets_MassVsChargedMultiplicity_"+cutVar.first).c_str(), 500, 0, 1500, 1000, 0, 100) );
        my_2d_histos.emplace( "h_ISRJets_MassVsNeutralMultiplicity_"+cutVar.first, std::make_shared<TH2D>( ("h_ISRJets_MassVsNeutralMultiplicity_"+cutVar.first).c_str(), ("h_ISRJets_MassVsNeutralMultiplicity_"+cutVar.first).c_str(), 500, 0, 1500, 1000, 0, 100) );
        my_2d_histos.emplace( "h_ISRJets_MassVSqgLikelihood_"+cutVar.first, std::make_shared<TH2D>( ("h_ISRJets_MassVSqgLikelihood_"+cutVar.first).c_str(), ("h_ISRJets_MassVSqgLikelihood_"+cutVar.first).c_str(), 500, 0, 1500, 100, 0, 1) );
        my_2d_histos.emplace( "h_ISRJets_MassVSbJetTagDeepCSVtotb_"+cutVar.first, std::make_shared<TH2D>( ("h_ISRJets_MassVSbJetTagDeepCSVtotb_"+cutVar.first).c_str(), ("h_ISRJets_MassVSbJetTagDeepCSVtotb_"+cutVar.first).c_str(), 500, 0, 1500, 100, 0, 1) );
        // Pt 
        my_2d_histos.emplace( "h_ISRJets_PtVsChargedHadronEnergyFraction_"+cutVar.first, std::make_shared<TH2D>( ("h_ISRJets_PtVsChargedHadronEnergyFraction_"+cutVar.first).c_str(), ("h_ISRJets_PtVsChargedHadronEnergyFraction_"+cutVar.first).c_str(), 100, 0, 1000, 100, 0, 1) );
        my_2d_histos.emplace( "h_ISRJets_PtVsNeutralHadronEnergyFraction_"+cutVar.first, std::make_shared<TH2D>( ("h_ISRJets_PtVsNeutralHadronEnergyFraction_"+cutVar.first).c_str(), ("h_ISRJets_PtVsNeutralHadronEnergyFraction_"+cutVar.first).c_str(), 100, 0, 1000, 100, 0, 1) );
        my_2d_histos.emplace( "h_ISRJets_PtVsChargedHadronMultiplicity_"+cutVar.first, std::make_shared<TH2D>( ("h_ISRJets_PtVsChargedHadronMultiplicity_"+cutVar.first).c_str(), ("h_ISRJets_PtVsChargedHadronMultiplicity_"+cutVar.first).c_str(), 100, 0, 1000, 1000, 0, 100) );
        my_2d_histos.emplace( "h_ISRJets_PtVsNeutralHadronMultiplicity_"+cutVar.first, std::make_shared<TH2D>( ("h_ISRJets_PtVsNeutralHadronMultiplicity_"+cutVar.first).c_str(), ("h_ISRJets_PtVsNeutralHadronMultiplicity_"+cutVar.first).c_str(), 100, 0, 1000, 1000, 0, 100) );
        my_2d_histos.emplace( "h_ISRJets_PtVsChargedMultiplicity_"+cutVar.first, std::make_shared<TH2D>( ("h_ISRJets_PtVsChargedMultiplicity_"+cutVar.first).c_str(), ("h_ISRJets_PtVsChargedMultiplicity_"+cutVar.first).c_str(), 100, 0, 1000, 1000, 0, 100) );       
        my_2d_histos.emplace( "h_ISRJets_PtVsNeutralMultiplicity_"+cutVar.first, std::make_shared<TH2D>( ("h_ISRJets_PtVsNeutralMultiplicity_"+cutVar.first).c_str(), ("h_ISRJets_PtVsNeutralMultiplicity_"+cutVar.first).c_str(), 100, 0, 1000, 1000, 0, 100) );
        my_2d_histos.emplace( "h_ISRJets_PtVSqgLikelihood_"+cutVar.first, std::make_shared<TH2D>( ("h_ISRJets_PtVSqgLikelihood_"+cutVar.first).c_str(), ("h_ISRJets_PtVSqgLikelihood_"+cutVar.first).c_str(), 100, 0, 1000, 100, 0, 1) );
        my_2d_histos.emplace( "h_ISRJets_PtVSbJetTagDeepCSVtotb_"+cutVar.first, std::make_shared<TH2D>( ("h_ISRJets_PtVSbJetTagDeepCSVtotb_"+cutVar.first).c_str(), ("h_ISRJets_PtVSbJetTagDeepCSVtotb_"+cutVar.first).c_str(), 100, 0, 1000, 100, 0, 1) );

        // --------------------------------
        // NonISR Jets - by using TreeMaker
        // -------------------------------- 
        my_2d_histos.emplace( "h_NonISRJets_EtaVsPhi_"+cutVar.first, std::make_shared<TH2D>( ("h_NonISRJets_EtaVsPhi_"+cutVar.first).c_str(), ("h_NonISRJets_EtaVsPhi_"+cutVar.first).c_str(), 100, -6, 6, 80, -4, 4) ); 
        my_2d_histos.emplace( "h_NonISRJets_MassVsPt_"+cutVar.first, std::make_shared<TH2D>( ("h_NonISRJets_MassVsPt_"+cutVar.first).c_str(), ("h_NonISRJets_MassVsPt_"+cutVar.first).c_str(), 500, 0, 1500, 100, 0, 1000) );
        // Mass
        my_2d_histos.emplace( "h_NonISRJets_MassVsChargedHadronEnergyFraction_"+cutVar.first, std::make_shared<TH2D>( ("h_NonISRJets_MassVsChargedHadronEnergyFraction_"+cutVar.first).c_str(), ("h_NonISRJets_MassVsChargedHadronEnergyFraction_"+cutVar.first).c_str(), 500, 0, 1500, 100, 0, 1) );
        my_2d_histos.emplace( "h_NonISRJets_MassVsNeutralHadronEnergyFraction_"+cutVar.first, std::make_shared<TH2D>( ("h_NonISRJets_MassVsNeutralHadronEnergyFraction_"+cutVar.first).c_str(), ("h_NonISRJets_MassVsNeutralHadronEnergyFraction_"+cutVar.first).c_str(), 500, 0, 1500, 100, 0, 1) );
        my_2d_histos.emplace( "h_NonISRJets_MassVsChargedHadronMultiplicity_"+cutVar.first, std::make_shared<TH2D>( ("h_NonISRJets_MassVsChargedHadronMultiplicity_"+cutVar.first).c_str(), ("h_NonISRJets_MassVsChargedHadronMultiplicity_"+cutVar.first).c_str(), 500, 0, 1500, 1000, 0, 100) );
        my_2d_histos.emplace( "h_NonISRJets_MassVsNeutralHadronMultiplicity_"+cutVar.first, std::make_shared<TH2D>( ("h_NonISRJets_MassVsNeutralHadronMultiplicity_"+cutVar.first).c_str(), ("h_NonISRJets_MassVsNeutralHadronMultiplicity_"+cutVar.first).c_str(), 500, 0, 1500, 1000, 0, 100) );        
        my_2d_histos.emplace( "h_NonISRJets_MassVsChargedMultiplicity_"+cutVar.first, std::make_shared<TH2D>( ("h_NonISRJets_MassVsChargedMultiplicity_"+cutVar.first).c_str(), ("h_NonISRJets_MassVsChargedMultiplicity_"+cutVar.first).c_str(), 500, 0, 1500, 1000, 0, 100) );
        my_2d_histos.emplace( "h_NonISRJets_MassVsNeutralMultiplicity_"+cutVar.first, std::make_shared<TH2D>( ("h_NonISRJets_MassVsNeutralMultiplicity_"+cutVar.first).c_str(), ("h_NonISRJets_MassVsNeutralMultiplicity_"+cutVar.first).c_str(), 500, 0, 1500, 1000, 0, 100) );        
        my_2d_histos.emplace( "h_NonISRJets_MassVSqgLikelihood_"+cutVar.first, std::make_shared<TH2D>( ("h_NonISRJets_MassVSqgLikelihood_"+cutVar.first).c_str(), ("h_NonISRJets_MassVSqgLikelihood_"+cutVar.first).c_str(), 500, 0, 1500, 100, 0, 1) );      
        my_2d_histos.emplace( "h_NonISRJets_MassVSbJetTagDeepCSVtotb_"+cutVar.first, std::make_shared<TH2D>( ("h_NonISRJets_MassVSbJetTagDeepCSVtotb_"+cutVar.first).c_str(), ("h_NonISRJets_MassVSbJetTagDeepCSVtotb_"+cutVar.first).c_str(), 500, 0, 1500, 100, 0, 1) );
        // Pt
        my_2d_histos.emplace( "h_NonISRJets_PtVsChargedHadronEnergyFraction_"+cutVar.first, std::make_shared<TH2D>( ("h_NonISRJets_PtVsChargedHadronEnergyFraction_"+cutVar.first).c_str(), ("h_NonISRJets_PtVsChargedHadronEnergyFraction_"+cutVar.first).c_str(), 100, 0, 1000, 100, 0, 1) );
        my_2d_histos.emplace( "h_NonISRJets_PtVsNeutralHadronEnergyFraction_"+cutVar.first, std::make_shared<TH2D>( ("h_NonISRJets_PtVsNeutralHadronEnergyFraction_"+cutVar.first).c_str(), ("h_NonISRJets_PtVsNeutralHadronEnergyFraction_"+cutVar.first).c_str(), 100, 0, 1000, 100, 0, 1) );
        my_2d_histos.emplace( "h_NonISRJets_PtVsChargedHadronMultiplicity_"+cutVar.first, std::make_shared<TH2D>( ("h_NonISRJets_PtVsChargedHadronMultiplicity_"+cutVar.first).c_str(), ("h_NonISRJets_PtVsChargedHadronMultiplicity_"+cutVar.first).c_str(), 100, 0, 1000, 1000, 0, 100) );
        my_2d_histos.emplace( "h_NonISRJets_PtVsNeutralHadronMultiplicity_"+cutVar.first, std::make_shared<TH2D>( ("h_NonISRJets_PtVsNeutralHadronMultiplicity_"+cutVar.first).c_str(), ("h_NonISRJets_PtVsNeutralHadronMultiplicity_"+cutVar.first).c_str(), 100, 0, 1000, 1000, 0, 100) );
        my_2d_histos.emplace( "h_NonISRJets_PtVsChargedMultiplicity_"+cutVar.first, std::make_shared<TH2D>( ("h_NonISRJets_PtVsChargedMultiplicity_"+cutVar.first).c_str(), ("h_NonISRJets_PtVsChargedMultiplicity_"+cutVar.first).c_str(), 100, 0, 1000, 1000, 0, 100) ); 
        my_2d_histos.emplace( "h_NonISRJets_PtVsNeutralMultiplicity_"+cutVar.first, std::make_shared<TH2D>( ("h_NonISRJets_PtVsNeutralMultiplicity_"+cutVar.first).c_str(), ("h_NonISRJets_PtVsNeutralMultiplicity_"+cutVar.first).c_str(), 100, 0, 1000, 1000, 0, 100) );
        my_2d_histos.emplace( "h_NonISRJets_PtVSqgLikelihood_"+cutVar.first, std::make_shared<TH2D>( ("h_NonISRJets_PtVSqgLikelihood_"+cutVar.first).c_str(), ("h_NonISRJets_PtVSqgLikelihood_"+cutVar.first).c_str(), 100, 0, 1000, 100, 0, 1) );
        my_2d_histos.emplace( "h_NonISRJets_PtVSbJetTagDeepCSVtotb_"+cutVar.first, std::make_shared<TH2D>( ("h_NonISRJets_PtVSbJetTagDeepCSVtotb_"+cutVar.first).c_str(), ("h_NonISRJets_PtVSbJetTagDeepCSVtotb_"+cutVar.first).c_str(), 100, 0, 1000, 100, 0, 1) );

        // ---------------------------------
        // Other Jets - not ISR + not NonISR
        // ---------------------------------        
        my_2d_histos.emplace( "h_OtherJets_EtaVsPhi_"+cutVar.first, std::make_shared<TH2D>( ("h_OtherJets_EtaVsPhi_"+cutVar.first).c_str(), ("h_OtherJets_EtaVsPhi_"+cutVar.first).c_str(), 100, -6, 6, 80, -4, 4) );
        my_2d_histos.emplace( "h_OtherJets_MassVsPt_"+cutVar.first, std::make_shared<TH2D>( ("h_OtherJets_MassVsPt_"+cutVar.first).c_str(), ("h_OtherJets_MassVsPt_"+cutVar.first).c_str(), 500, 0, 1500, 100, 0, 1000) );
        // Mass
        my_2d_histos.emplace( "h_OtherJets_MassVsChargedHadronEnergyFraction_"+cutVar.first, std::make_shared<TH2D>( ("h_OtherJets_MassVsChargedHadronEnergyFraction_"+cutVar.first).c_str(), ("h_OtherJets_MassVsChargedHadronEnergyFraction_"+cutVar.first).c_str(), 500, 0, 1500, 100, 0, 1) );
        my_2d_histos.emplace( "h_OtherJets_MassVsNeutralHadronEnergyFraction_"+cutVar.first, std::make_shared<TH2D>( ("h_OtherJets_MassVsNeutralHadronEnergyFraction_"+cutVar.first).c_str(), ("h_OtherJets_MassVsNeutralHadronEnergyFraction_"+cutVar.first).c_str(), 500, 0, 1500, 100, 0, 1) );
        my_2d_histos.emplace( "h_OtherJets_MassVsChargedHadronMultiplicity_"+cutVar.first, std::make_shared<TH2D>( ("h_OtherJets_MassVsChargedHadronMultiplicity_"+cutVar.first).c_str(), ("h_OtherJets_MassVsChargedHadronMultiplicity_"+cutVar.first).c_str(), 500, 0, 1500, 1000, 0, 100) );
        my_2d_histos.emplace( "h_OtherJets_MassVsNeutralHadronMultiplicity_"+cutVar.first, std::make_shared<TH2D>( ("h_OtherJets_MassVsNeutralHadronMultiplicity_"+cutVar.first).c_str(), ("h_OtherJets_MassVsNeutralHadronMultiplicity_"+cutVar.first).c_str(), 500, 0, 1500, 1000, 0, 100) );
        my_2d_histos.emplace( "h_OtherJets_MassVsChargedMultiplicity_"+cutVar.first, std::make_shared<TH2D>( ("h_OtherJets_MassVsChargedMultiplicity_"+cutVar.first).c_str(), ("h_OtherJets_MassVsChargedMultiplicity_"+cutVar.first).c_str(), 500, 0, 1500, 1000, 0, 100) );
        my_2d_histos.emplace( "h_OtherJets_MassVsNeutralMultiplicity_"+cutVar.first, std::make_shared<TH2D>( ("h_OtherJets_MassVsNeutralMultiplicity_"+cutVar.first).c_str(), ("h_OtherJets_MassVsNeutralMultiplicity_"+cutVar.first).c_str(), 500, 0, 1500, 1000, 0, 100) );
        my_2d_histos.emplace( "h_OtherJets_MassVSqgLikelihood_"+cutVar.first, std::make_shared<TH2D>( ("h_OtherJets_MassVSqgLikelihood_"+cutVar.first).c_str(), ("h_OtherJets_MassVSqgLikelihood_"+cutVar.first).c_str(), 500, 0, 1500, 100, 0, 1) );
        my_2d_histos.emplace( "h_OtherJets_MassVSbJetTagDeepCSVtotb_"+cutVar.first, std::make_shared<TH2D>( ("h_OtherJets_MassVSbJetTagDeepCSVtotb_"+cutVar.first).c_str(), ("h_OtherJets_MassVSbJetTagDeepCSVtotb_"+cutVar.first).c_str(), 500, 0, 1500, 100, 0, 1) );
        // Pt
        my_2d_histos.emplace( "h_OtherJets_PtVsChargedHadronEnergyFraction_"+cutVar.first, std::make_shared<TH2D>( ("h_OtherJets_PtVsChargedHadronEnergyFraction_"+cutVar.first).c_str(), ("h_OtherJets_PtVsChargedHadronEnergyFraction_"+cutVar.first).c_str(), 100, 0, 1000, 100, 0, 1) );
        my_2d_histos.emplace( "h_OtherJets_PtVsNeutralHadronEnergyFraction_"+cutVar.first, std::make_shared<TH2D>( ("h_OtherJets_PtVsNeutralHadronEnergyFraction_"+cutVar.first).c_str(), ("h_OtherJets_PtVsNeutralHadronEnergyFraction_"+cutVar.first).c_str(), 100, 0, 1000, 100, 0, 1) );
        my_2d_histos.emplace( "h_OtherJets_PtVsChargedHadronMultiplicity_"+cutVar.first, std::make_shared<TH2D>( ("h_OtherJets_PtVsChargedHadronMultiplicity_"+cutVar.first).c_str(), ("h_OtherJets_PtVsChargedHadronMultiplicity_"+cutVar.first).c_str(), 100, 0, 1000, 1000, 0, 100) );
        my_2d_histos.emplace( "h_OtherJets_PtVsNeutralHadronMultiplicity_"+cutVar.first, std::make_shared<TH2D>( ("h_OtherJets_PtVsNeutralHadronMultiplicity_"+cutVar.first).c_str(), ("h_OtherJets_PtVsNeutralHadronMultiplicity_"+cutVar.first).c_str(), 100, 0, 1000, 1000, 0, 100) );        
        my_2d_histos.emplace( "h_OtherJets_PtVsChargedMultiplicity_"+cutVar.first, std::make_shared<TH2D>( ("h_OtherJets_PtVsChargedMultiplicity_"+cutVar.first).c_str(), ("h_OtherJets_PtVsChargedMultiplicity_"+cutVar.first).c_str(), 100, 0, 1000, 1000, 0, 100) );
        my_2d_histos.emplace( "h_OtherJets_PtVsNeutralMultiplicity_"+cutVar.first, std::make_shared<TH2D>( ("h_OtherJets_PtVsNeutralMultiplicity_"+cutVar.first).c_str(), ("h_OtherJets_PtVsNeutralMultiplicity_"+cutVar.first).c_str(), 100, 0, 1000, 1000, 0, 100) );        
        my_2d_histos.emplace( "h_OtherJets_PtVSqgLikelihood_"+cutVar.first, std::make_shared<TH2D>( ("h_OtherJets_PtVSqgLikelihood_"+cutVar.first).c_str(), ("h_OtherJets_PtVSqgLikelihood_"+cutVar.first).c_str(), 100, 0, 1000, 100, 0, 1) );
        my_2d_histos.emplace( "h_OtherJets_PtVSbJetTagDeepCSVtotb_"+cutVar.first, std::make_shared<TH2D>( ("h_OtherJets_PtVSbJetTagDeepCSVtotb_"+cutVar.first).c_str(), ("h_OtherJets_PtVSbJetTagDeepCSVtotb_"+cutVar.first).c_str(), 100, 0, 1000, 100, 0, 1) );

    }
}

// ---------------------------------------------
// -- Put everything you want to do per event 
// ---------------------------------------------
void ISRJets_Analyzer::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    while( tr.getNextEvent() )
    {
        const auto& eventCounter    = tr.getVar<int>("eventCounter");
        
        //-------------------------
        // -- Print Event Number 
        //-------------------------
        
        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() & (10000 == 0) ) printf( " Event %i\n", tr.getEvtNum() );

        const auto& runtype               = tr.getVar<std::string>("runtype");     
        const auto& Jets                  = tr.getVec<TLorentzVector>("Jets");
        const auto& GoodJets_pt20         = tr.getVec<bool>("GoodJets_pt20");
        const auto& GoodBJets_pt45        = tr.getVec<bool>("GoodBJets_pt45");
        //const auto& GenParticles          = tr.getVec<TLorentzVector>("GenParticles");
        const auto& passBaseline0l        = tr.getVar<bool>("passBaseline0l_Good");
        const auto& ntops                 = tr.getVar<int>("ntops");
        const auto& topsLV                = tr.getVec<TLorentzVector>("topsLV");
        const auto& dR_bjets              = tr.getVar<double>("dR_bjets");
        const bool  pass_ge2t             = ntops >= 2;
        const bool  pass_ge1dRbjets       = dR_bjets >= 1.0;       

        // ----------------------
        // -- ISR gen matching
        // ----------------------
        const auto& GM_ISRmatching_allDR                    = tr.getVec<double>("GM_ISRmatching_allDR");
        const auto& GM_ISRmatching_bestDR                   = tr.getVec<double>("GM_ISRmatching_bestDR");
        const auto& GM_ISRmatching_justCutOnDR_DR           = tr.getVec<double>("GM_ISRmatching_justCutOnDR_DR");
        const auto& GM_ISRmatching_justCutOnPtRatio_DR      = tr.getVec<double>("GM_ISRmatching_justCutOnPtRatio_DR");
        const auto& GM_ISRmatching_allPtRatio               = tr.getVec<double>("GM_ISRmatching_allPtRatio");
        const auto& GM_ISRmatching_bestPtRatio              = tr.getVec<double>("GM_ISRmatching_bestPtRatio");
        const auto& GM_ISRmatching_justCutOnDR_PtRatio      = tr.getVec<double>("GM_ISRmatching_justCutOnDR_PtRatio");
        const auto& GM_ISRmatching_justCutOnPtRatio_PtRatio = tr.getVec<double>("GM_ISRmatching_justCutOnPtRatio_PtRatio");
        const auto& nGenISR                                 = tr.getVar<int>("nGenISR");
        const auto& nGenISR_qg                              = tr.getVar<int>("nGenISR_qg");
        const auto& nGenISR_gq                              = tr.getVar<int>("nGenISR_gq");
        const auto& nGenISR_gg                              = tr.getVar<int>("nGenISR_gg");
        const auto& GenISR                                  = tr.getVec<TLorentzVector>("GenISR");
        const auto& GenISR_qg                               = tr.getVec<TLorentzVector>("GenISR_qg");
        const auto& GenISR_gq                               = tr.getVec<TLorentzVector>("GenISR_gq");
        const auto& GenISR_gg                               = tr.getVec<TLorentzVector>("GenISR_gg");
        const auto& dEta_RecoISR_GenISR                     = tr.getVec<double>("dEta_RecoISR_GenISR");
        const auto& ISRmatched_dr_ptr                       = tr.getVec<bool>("ISRmatched_dr_ptr");
        const auto& NonISRmatched_dr_ptr                    = tr.getVec<bool>("NonISRmatched_dr_ptr"); 
        // AK4 variables from ntuple
        const auto& Jets_axismajor                          = tr.getVec<double>("Jets_axismajor");
        const auto& Jets_axisminor                          = tr.getVec<double>("Jets_axisminor");
        const auto& Jets_ptD                                = tr.getVec<double>("Jets_ptD");
        const auto& Jets_chargedHadronEnergyFraction        = tr.getVec<double>("Jets_chargedHadronEnergyFraction");
        const auto& Jets_neutralHadronEnergyFraction        = tr.getVec<double>("Jets_neutralHadronEnergyFraction");
        const auto& Jets_chargedHadronMultiplicity          = tr.getVec<int>("Jets_chargedHadronMultiplicity");
        const auto& Jets_neutralHadronMultiplicity          = tr.getVec<int>("Jets_neutralHadronMultiplicity");
        const auto& Jets_chargedMultiplicity                = tr.getVec<int>("Jets_chargedMultiplicity");
        const auto& Jets_neutralMultiplicity                = tr.getVec<int>("Jets_neutralMultiplicity");
        const auto& Jets_qgLikelihood                       = tr.getVec<double>("Jets_qgLikelihood");
        const auto& Jets_bJetTagDeepCSVtotb                 = tr.getVec<double>("Jets_bJetTagDeepCSVtotb");

        const auto& Jets_ISRMask                            = tr.getVec<bool>("Jets_ISRMask");

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
            //weight *= eventweight;
        }

        // --------------------------------------------
        // DeltaR between each ISR and each NonISR jets
        // --------------------------------------------
        std::vector<TLorentzVector> ISRJets;
        std::vector<TLorentzVector> NonISRJets;
        for(unsigned int j = 0; j < Jets.size(); j++) 
        {
            if(!GoodJets_pt20[j]) continue;
            
            if (ISRmatched_dr_ptr[j]) 
            {   
                ISRJets.push_back(Jets.at(j));
            }

            if (NonISRmatched_dr_ptr[j])
            {   
                NonISRJets.push_back(Jets.at(j));
            }
        }
        
        std::vector<double> dR_ISR_NonISR;
        for (unsigned int i = 0; i < ISRJets.size(); i++)
        {
            double tempDR = 999; 
            
            for (unsigned int n = 0; n < NonISRJets.size(); n++)
            {
                double deltaR = ISRJets.at(i).DeltaR(NonISRJets.at(n));
                if (deltaR < tempDR) tempDR = deltaR; 
            }
            
            dR_ISR_NonISR.push_back(tempDR);
        }

        // -------------------------------------------
        // DeltaR between closet top and each ISR jets
        // -------------------------------------------
        std::vector<double> dR_top_ISR; 
        for (unsigned int j = 0; j < ISRJets.size(); j++)
        {   
            double tempDR = 999; // keep track of smallest dR
         
            for (unsigned int t = 0; t < topsLV.size(); t++)
            {   
                double deltaR = topsLV.at(t).DeltaR(ISRJets.at(j));
                if (deltaR < tempDR) tempDR = deltaR; // set temp variable to newest smaller dR
            }
            
            dR_top_ISR.push_back(tempDR);
        } 

        // -----------------------------------------------
        // DeltaR between closest top and each NonISR jets
        // -----------------------------------------------
        std::vector<double> dR_top_NonISR;
        for (unsigned int j = 0; j < NonISRJets.size(); j++)
        {
            double tempDR = 999; // keep track of smallest dR
            
            for (unsigned int t = 0; t < topsLV.size(); t++)
            {
                double deltaR = topsLV.at(t).DeltaR(NonISRJets.at(j));
                if (deltaR < tempDR) tempDR = deltaR; // set temp variable to newest smaller dR
            }
            
            dR_top_NonISR.push_back(tempDR);
        }

        // ---------------------------------------------
        // DeltaR between closest bjet and each ISR jets  
        // ---------------------------------------------
        std::vector<TLorentzVector> bjets;
        for(unsigned int ijet = 0; ijet < Jets.size(); ijet++) 
        {
            if(!GoodBJets_pt45[ijet]) continue;
            bjets.push_back(Jets.at(ijet));        
        } 
    
        std::vector<double> dR_bjet_ISR;
        for (unsigned int j = 0; j < ISRJets.size(); j++)
        {
            double tempDR = 999; 

            for (unsigned int t = 0; t < bjets.size(); t++)
            {
                double deltaR = bjets.at(t).DeltaR(ISRJets.at(j));
                if (deltaR < tempDR) tempDR = deltaR; 
            }

            dR_bjet_ISR.push_back(tempDR);
        }

        // ------------------------------------------------
        // DeltaR between closest bjet and each NonISR jets
        // ------------------------------------------------
        std::vector<double> dR_bjet_NonISR;
        for (unsigned int j = 0; j < NonISRJets.size(); j++)
        {
            double tempDR = 999; 

            for (unsigned int t = 0; t < bjets.size(); t++)
            {
                double deltaR = bjets.at(t).DeltaR(NonISRJets.at(j));
                if (deltaR < tempDR) tempDR = deltaR;
            }

            dR_bjet_NonISR.push_back(tempDR);
        }

        // -------------------------------------------------
        // -- Make cuts and fill histograms here & cutmap
        // -------------------------------------------------
        const std::map<std::string, bool>& cutmap
        {
            {"", true},
            {"0l_baseline", passBaseline0l && pass_ge2t && pass_ge1dRbjets},
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
                // -------------------------------
                // -- ISR gen matching variables
                // -------------------------------
                // --------------------
                // all possible matches
                // --------------------
                for (double dr = 0.0; dr < GM_ISRmatching_allDR.size(); dr++)
                {
                    my_histos["h_GM_ISRmatching_allDR_"+cutVar.first]->Fill( GM_ISRmatching_allDR.at(dr), weight );
                    my_histos["h_GM_ISRmatching_allDR_"+cutVar.first]->GetXaxis()->SetTitle("ISRmatching all #DeltaR");
                    my_histos["h_GM_ISRmatching_allDR_"+cutVar.first]->GetYaxis()->SetTitle("Events");

                    // 2D histogram - dr vs pt ratio
                    my_2d_histos["h_GM_all_PtRvsDR_"+cutVar.first]->Fill(GM_ISRmatching_allPtRatio.at(dr), GM_ISRmatching_allDR.at(dr), weight);
                    my_2d_histos["h_GM_all_PtRvsDR_"+cutVar.first]->GetXaxis()->SetTitle("ISRmatching all (Pt^{Reco}/Pt^{Gen})");
                    my_2d_histos["h_GM_all_PtRvsDR_"+cutVar.first]->GetYaxis()->SetTitle("ISRmatching all #DeltaR");
                }

                for (double pt = 0.0; pt < GM_ISRmatching_allPtRatio.size(); pt++)
                {
                    my_histos["h_GM_ISRmatching_allPtRatio_"+cutVar.first]->Fill( GM_ISRmatching_allPtRatio.at(pt), weight );
                    my_histos["h_GM_ISRmatching_allPtRatio_"+cutVar.first]->GetXaxis()->SetTitle("ISRmatching all (Pt^{Reco}/Pt^{Gen})");
                    my_histos["h_GM_ISRmatching_allPtRatio_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                }

                // -------------------------------------------
                // matches with cutting on deltaR and pt ratio
                // -------------------------------------------
                for (double dr = 0.0; dr < GM_ISRmatching_bestDR.size(); dr++)
                {
                    my_histos["h_GM_ISRmatching_bestDR_"+cutVar.first]->Fill( GM_ISRmatching_bestDR.at(dr), weight );
                    my_histos["h_GM_ISRmatching_bestDR_"+cutVar.first]->GetXaxis()->SetTitle("ISRmatching best #DeltaR");
                    my_histos["h_GM_ISRmatching_bestDR_"+cutVar.first]->GetYaxis()->SetTitle("Events");

                    // 2D histogram - dr vs pt ratio
                    my_2d_histos["h_GM_best_PtRvsDR_"+cutVar.first]->Fill(GM_ISRmatching_bestPtRatio.at(dr), GM_ISRmatching_bestDR.at(dr), weight);
                    my_2d_histos["h_GM_best_PtRvsDR_"+cutVar.first]->GetXaxis()->SetTitle("ISRmatching best (Pt^{Reco}/Pt^{Gen})");
                    my_2d_histos["h_GM_best_PtRvsDR_"+cutVar.first]->GetYaxis()->SetTitle("ISRmatching best #DeltaR");
                }

                for (double pt = 0.0; pt < GM_ISRmatching_bestPtRatio.size(); pt++)
                {
                    my_histos["h_GM_ISRmatching_bestPtRatio_"+cutVar.first]->Fill( GM_ISRmatching_bestPtRatio.at(pt), weight );
                    my_histos["h_GM_ISRmatching_bestPtRatio_"+cutVar.first]->GetXaxis()->SetTitle("ISRmatching best (Pt^{Reco}/Pt^{Gen})");
                    my_histos["h_GM_ISRmatching_bestPtRatio_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                }

                // ------------------------------
                // matches with cutting on deltaR
                // ------------------------------
                for (double dr = 0.0; dr < GM_ISRmatching_justCutOnDR_DR.size(); dr++)
                {
                    my_histos["h_GM_ISRmatching_CutOnDR_DR_"+cutVar.first]->Fill( GM_ISRmatching_justCutOnDR_DR.at(dr), weight );
                    my_histos["h_GM_ISRmatching_CutOnDR_DR_"+cutVar.first]->GetXaxis()->SetTitle("ISRmatching CutOnDR #DeltaR");
                    my_histos["h_GM_ISRmatching_CutOnDR_DR_"+cutVar.first]->GetYaxis()->SetTitle("Events");

                    // 2D histogram - dr vs pt ratio
                    my_2d_histos["h_GM_CutOnDR_PtRvsDR_"+cutVar.first]->Fill(GM_ISRmatching_justCutOnDR_PtRatio.at(dr), GM_ISRmatching_justCutOnDR_DR.at(dr), weight);
                    my_2d_histos["h_GM_CutOnDR_PtRvsDR_"+cutVar.first]->GetXaxis()->SetTitle("ISRmatching CutOnDR (Pt^{Reco}/Pt^{Gen})");
                    my_2d_histos["h_GM_CutOnDR_PtRvsDR_"+cutVar.first]->GetYaxis()->SetTitle("ISRmatching CutOnDR #DeltaR");   
                }

                for (double pt = 0.0; pt < GM_ISRmatching_justCutOnDR_PtRatio.size(); pt++)
                {
                    my_histos["h_GM_ISRmatching_CutOnDR_PtRatio_"+cutVar.first]->Fill( GM_ISRmatching_justCutOnDR_PtRatio.at(pt), weight );
                    my_histos["h_GM_ISRmatching_CutOnDR_PtRatio_"+cutVar.first]->GetXaxis()->SetTitle("ISRmatching CutOnDR (Pt^{Reco}/Pt^{Gen})");
                    my_histos["h_GM_ISRmatching_CutOnDR_PtRatio_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                }

                // --------------------------------
                // matches with cutting on pt ratio
                // --------------------------------
                for (double dr = 0.0; dr < GM_ISRmatching_justCutOnPtRatio_DR.size(); dr++)
                {
                    my_histos["h_GM_ISRmatching_CutOnPtRatio_DR_"+cutVar.first]->Fill( GM_ISRmatching_justCutOnPtRatio_DR.at(dr), weight );
                    my_histos["h_GM_ISRmatching_CutOnPtRatio_DR_"+cutVar.first]->GetXaxis()->SetTitle("ISRmatching CutOnPtRatio #DeltaR");
                    my_histos["h_GM_ISRmatching_CutOnPtRatio_DR_"+cutVar.first]->GetYaxis()->SetTitle("Events");                
                    // 2D histogram - dr vs pt ratio
                    my_2d_histos["h_GM_CutOnPtRatio_PtRvsDR_"+cutVar.first]->Fill(GM_ISRmatching_justCutOnPtRatio_PtRatio.at(dr), GM_ISRmatching_justCutOnPtRatio_DR.at(dr), weight);
                    my_2d_histos["h_GM_CutOnPtRatio_PtRvsDR_"+cutVar.first]->GetXaxis()->SetTitle("ISRmatching CutOnPtRatio (Pt^{Reco}/Pt^{Gen})");
                    my_2d_histos["h_GM_CutOnPtRatio_PtRvsDR_"+cutVar.first]->GetYaxis()->SetTitle("ISRmatching CutOnPtRatio #DeltaR");
                }

                for (double pt = 0.0; pt < GM_ISRmatching_justCutOnPtRatio_PtRatio.size(); pt++)
                {
                    my_histos["h_GM_ISRmatching_CutOnPtRatio_PtRatio_"+cutVar.first]->Fill( GM_ISRmatching_justCutOnPtRatio_PtRatio.at(pt), weight );
                    my_histos["h_GM_ISRmatching_CutOnPtRatio_PtRatio_"+cutVar.first]->GetXaxis()->SetTitle("ISRmatching CutOnPtRatio (Pt^{Reco}/Pt^{Gen})");
                    my_histos["h_GM_ISRmatching_CutOnPtRatio_PtRatio_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                }            

                // -----------------------------------------
                // -- Get the Gen ISR kinematic variables
                // -----------------------------------------
                my_histos["h_nGenISR_"+cutVar.first]->Fill( nGenISR, weight );
                my_histos["h_nGenISR_"+cutVar.first]->GetXaxis()->SetTitle("nGenISR");
                my_histos["h_nGenISR_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                for (unsigned int i = 0;  i < GenISR.size(); i ++)
                {
                    my_histos["h_GenISR_Mass_"+cutVar.first]->Fill( GenISR.at(i).M(), weight );
                    my_histos["h_GenISR_Mass_"+cutVar.first]->GetXaxis()->SetTitle("GenISR Mass");
                    my_histos["h_GenISR_Mass_"+cutVar.first]->GetYaxis()->SetTitle("Events"); 
                    my_histos["h_GenISR_Pt_"+cutVar.first]->Fill( GenISR.at(i).Pt(), weight );
                    my_histos["h_GenISR_Pt_"+cutVar.first]->GetXaxis()->SetTitle("GenISR Pt");
                    my_histos["h_GenISR_Pt_"+cutVar.first]->GetYaxis()->SetTitle("Events"); 
                    my_histos["h_GenISR_Phi_"+cutVar.first]->Fill( GenISR.at(i).Phi(), weight );
                    my_histos["h_GenISR_Phi_"+cutVar.first]->GetXaxis()->SetTitle("GenISR Phi");
                    my_histos["h_GenISR_Phi_"+cutVar.first]->GetYaxis()->SetTitle("Events"); 
                    my_histos["h_GenISR_Eta_"+cutVar.first]->Fill( GenISR.at(i).Eta(), weight );
                    my_histos["h_GenISR_Eta_"+cutVar.first]->GetXaxis()->SetTitle("GenISR Eta");
                    my_histos["h_GenISR_Eta_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                }
                
                my_histos["h_nGenISR_qg_"+cutVar.first]->Fill( nGenISR_qg, weight );
                my_histos["h_nGenISR_qg_"+cutVar.first]->GetXaxis()->SetTitle("nGenISR_qg");
                my_histos["h_nGenISR_qg_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                for (unsigned int i = 0;  i < GenISR_qg.size(); i ++)
                {
                    my_histos["h_GenISR_qg_Mass_"+cutVar.first]->Fill( GenISR_qg.at(i).M(), weight );
                    my_histos["h_GenISR_qg_Mass_"+cutVar.first]->GetXaxis()->SetTitle("GenISR_qg Mass");
                    my_histos["h_GenISR_qg_Mass_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                    my_histos["h_GenISR_qg_Pt_"+cutVar.first]->Fill( GenISR_qg.at(i).Pt(), weight );
                    my_histos["h_GenISR_qg_Pt_"+cutVar.first]->GetXaxis()->SetTitle("GenISR_qg Pt");
                    my_histos["h_GenISR_qg_Pt_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                    my_histos["h_GenISR_qg_Phi_"+cutVar.first]->Fill( GenISR_qg.at(i).Phi(), weight );
                    my_histos["h_GenISR_qg_Phi_"+cutVar.first]->GetXaxis()->SetTitle("GenISR_qg Phi");
                    my_histos["h_GenISR_qg_Phi_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                    my_histos["h_GenISR_qg_Eta_"+cutVar.first]->Fill( GenISR_qg.at(i).Eta(), weight );
                    my_histos["h_GenISR_qg_Eta_"+cutVar.first]->GetXaxis()->SetTitle("GenISR_qg Eta");
                    my_histos["h_GenISR_qg_Eta_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                }

                my_histos["h_nGenISR_gq_"+cutVar.first]->Fill( nGenISR_gq, weight );
                my_histos["h_nGenISR_gq_"+cutVar.first]->GetXaxis()->SetTitle("nGenISR_gq");
                my_histos["h_nGenISR_gq_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                for (unsigned int i = 0;  i < GenISR_gq.size(); i ++)
                {
                    my_histos["h_GenISR_gq_Mass_"+cutVar.first]->Fill( GenISR_gq.at(i).M(), weight );
                    my_histos["h_GenISR_gq_Mass_"+cutVar.first]->GetXaxis()->SetTitle("GenISR_gq Mass");
                    my_histos["h_GenISR_gq_Mass_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                    my_histos["h_GenISR_gq_Pt_"+cutVar.first]->Fill( GenISR_gq.at(i).Pt(), weight );
                    my_histos["h_GenISR_gq_Pt_"+cutVar.first]->GetXaxis()->SetTitle("GenISR_gq Pt");
                    my_histos["h_GenISR_gq_Pt_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                    my_histos["h_GenISR_gq_Phi_"+cutVar.first]->Fill( GenISR_gq.at(i).Phi(), weight );
                    my_histos["h_GenISR_gq_Phi_"+cutVar.first]->GetXaxis()->SetTitle("GenISR_gq Phi");
                    my_histos["h_GenISR_gq_Phi_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                    my_histos["h_GenISR_gq_Eta_"+cutVar.first]->Fill( GenISR_gq.at(i).Eta(), weight );
                    my_histos["h_GenISR_gq_Eta_"+cutVar.first]->GetXaxis()->SetTitle("GenISR_gq Eta");
                    my_histos["h_GenISR_gq_Eta_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                }

                my_histos["h_nGenISR_gg_"+cutVar.first]->Fill( nGenISR_gg, weight );
                my_histos["h_nGenISR_gg_"+cutVar.first]->GetXaxis()->SetTitle("nGenISR_gg");
                my_histos["h_nGenISR_gg_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                for (unsigned int i = 0;  i < GenISR_gg.size(); i ++)
                {
                    my_histos["h_GenISR_gg_Mass_"+cutVar.first]->Fill( GenISR_gg.at(i).M(), weight );
                    my_histos["h_GenISR_gg_Mass_"+cutVar.first]->GetXaxis()->SetTitle("GenISR_gg Mass");
                    my_histos["h_GenISR_gg_Mass_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                    my_histos["h_GenISR_gg_Pt_"+cutVar.first]->Fill( GenISR_gg.at(i).Pt(), weight );
                    my_histos["h_GenISR_gg_Pt_"+cutVar.first]->GetXaxis()->SetTitle("GenISR_gg Pt");
                    my_histos["h_GenISR_gg_Pt_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                    my_histos["h_GenISR_gg_Phi_"+cutVar.first]->Fill( GenISR_gg.at(i).Phi(), weight );
                    my_histos["h_GenISR_gg_Phi_"+cutVar.first]->GetXaxis()->SetTitle("GenISR_gg Phi");
                    my_histos["h_GenISR_gg_Phi_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                    my_histos["h_GenISR_gg_Eta_"+cutVar.first]->Fill( GenISR_gg.at(i).Eta(), weight );
                    my_histos["h_GenISR_gg_Eta_"+cutVar.first]->GetXaxis()->SetTitle("GenISR_gg Eta");
                    my_histos["h_GenISR_gg_Eta_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                }

                // ---------------------------------
                // -- Testing Jets_ISRMask filter
                // ---------------------------------
                for (unsigned int j = 0; j < Jets.size(); j++)
                {
                    if (!GoodJets_pt20[j]) continue;

                    if (Jets_ISRMask[j])
                    {
                        my_histos["h_Jets_ISRMask_Mass_"+cutVar.first]->Fill( Jets[j].M(), weight );
                        my_histos["h_Jets_ISRMask_Mass_"+cutVar.first]->GetXaxis()->SetTitle("Jets_ISRMask Mass");
                        my_histos["h_Jets_ISRMask_Mass_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_Jets_ISRMask_Pt_"+cutVar.first]->Fill( Jets[j].Pt(), weight );
                        my_histos["h_Jets_ISRMask_Pt_"+cutVar.first]->GetXaxis()->SetTitle("Jets_ISRMask Pt");
                        my_histos["h_Jets_ISRMask_Pt_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_Jets_ISRMask_Phi_"+cutVar.first]->Fill( Jets[j].Phi(), weight );
                        my_histos["h_Jets_ISRMask_Phi_"+cutVar.first]->GetXaxis()->SetTitle("Jets_ISRMask #phi");
                        my_histos["h_Jets_ISRMask_Phi_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_Jets_ISRMask_Eta_"+cutVar.first]->Fill( Jets[j].Eta(), weight );
                        my_histos["h_Jets_ISRMask_Eta_"+cutVar.first]->GetXaxis()->SetTitle("Jets_ISRMask #eta");
                        my_histos["h_Jets_ISRMask_Eta_"+cutVar.first]->GetYaxis()->SetTitle("Events");   
                    }
                }


                // -----------------------------------------
                // -- ISR & NonISR & Other Jets variables
                // -----------------------------------------
                int nISRJets    = 0;
                int nNonISRJets = 0;
                int nOtherJets  = 0;

                for (unsigned int j = 0; j < Jets.size(); j++)
                {
                    if (!GoodJets_pt20[j]) continue;
                
                    // --------------------------------------------------------------------
                    // ISR Jets - matches cutting on dr - by using truth information of ISR
                    // -------------------------------------------------------------------- 
                    if (ISRmatched_dr_ptr[j])
                    {
                        // 1D
                        my_histos["h_ISRJets_Mass_"+cutVar.first]->Fill( Jets[j].M(), weight );
                        my_histos["h_ISRJets_Mass_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets Mass");
                        my_histos["h_ISRJets_Mass_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_ISRJets_Pt_"+cutVar.first]->Fill( Jets[j].Pt(), weight );
                        my_histos["h_ISRJets_Pt_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets Pt");
                        my_histos["h_ISRJets_Pt_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_ISRJets_Phi_"+cutVar.first]->Fill( Jets[j].Phi(), weight );
                        my_histos["h_ISRJets_Phi_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets #phi");
                        my_histos["h_ISRJets_Phi_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_ISRJets_Eta_"+cutVar.first]->Fill( Jets[j].Eta(), weight );
                        my_histos["h_ISRJets_Eta_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets #eta");
                        my_histos["h_ISRJets_Eta_"+cutVar.first]->GetYaxis()->SetTitle("Events");        
                        my_histos["h_ISRJets_axismajor_"+cutVar.first]->Fill( Jets_axismajor[j], weight );
                        my_histos["h_ISRJets_axismajor_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets axismajor");
                        my_histos["h_ISRJets_axismajor_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_ISRJets_axisminor_"+cutVar.first]->Fill( Jets_axisminor[j], weight );
                        my_histos["h_ISRJets_axisminor_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets axisminor");
                        my_histos["h_ISRJets_axisminor_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_ISRJets_chargedHadronEnergyFraction_"+cutVar.first]->Fill( Jets_chargedHadronEnergyFraction[j], weight );
                        my_histos["h_ISRJets_chargedHadronEnergyFraction_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets Energy fraction");
                        my_histos["h_ISRJets_chargedHadronEnergyFraction_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_ISRJets_neutralHadronEnergyFraction_"+cutVar.first]->Fill( Jets_neutralHadronEnergyFraction[j], weight );
                        my_histos["h_ISRJets_neutralHadronEnergyFraction_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets Energy fraction");
                        my_histos["h_ISRJets_neutralHadronEnergyFraction_"+cutVar.first]->GetYaxis()->SetTitle("Events");   
                        my_histos["h_ISRJets_ptD_"+cutVar.first]->Fill( Jets_ptD[j], weight );
                        my_histos["h_ISRJets_ptD_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets ptD");
                        my_histos["h_ISRJets_ptD_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_ISRJets_chargedHadronMultiplicity_"+cutVar.first]->Fill( Jets_chargedHadronMultiplicity[j], weight );
                        my_histos["h_ISRJets_chargedHadronMultiplicity_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets CH Multiplicity");
                        my_histos["h_ISRJets_chargedHadronMultiplicity_"+cutVar.first]->GetYaxis()->SetTitle("Event");
                        my_histos["h_ISRJets_neutralHadronMultiplicity_"+cutVar.first]->Fill( Jets_neutralHadronMultiplicity[j], weight );
                        my_histos["h_ISRJets_neutralHadronMultiplicity_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets NH Multiplicity");
                        my_histos["h_ISRJets_neutralHadronMultiplicity_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_ISRJets_chargedMultiplicity_"+cutVar.first]->Fill( Jets_chargedMultiplicity[j], weight );
                        my_histos["h_ISRJets_chargedMultiplicity_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets Charged Multiplicity");
                        my_histos["h_ISRJets_chargedMultiplicity_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_ISRJets_neutralMultiplicity_"+cutVar.first]->Fill( Jets_neutralMultiplicity[j], weight );
                        my_histos["h_ISRJets_neutralMultiplicity_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets Neural Multiplicity");
                        my_histos["h_ISRJets_neutralMultiplicity_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_ISRJets_qgLikelihood_"+cutVar.first]->Fill( Jets_qgLikelihood[j], weight );
                        my_histos["h_ISRJets_qgLikelihood_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets qgLikelihood");
                        my_histos["h_ISRJets_qgLikelihood_"+cutVar.first]->GetYaxis()->SetTitle("Events");    
                        my_histos["h_ISRJets_bJetTagDeepCSVtotb_"+cutVar.first]->Fill( Jets_bJetTagDeepCSVtotb[j], weight ); 
                        my_histos["h_ISRJets_bJetTagDeepCSVtotb_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets bJetTagDeepCSVtotb");
                        my_histos["h_ISRJets_bJetTagDeepCSVtotb_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        // 2D 
                        my_2d_histos["h_ISRJets_EtaVsPhi_"+cutVar.first]->Fill( Jets[j].Eta(), Jets[j].Phi(), weight );
                        my_2d_histos["h_ISRJets_EtaVsPhi_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets #eta");
                        my_2d_histos["h_ISRJets_EtaVsPhi_"+cutVar.first]->GetYaxis()->SetTitle("ISRJets #phi");
                        my_2d_histos["h_ISRJets_MassVsPt_"+cutVar.first]->Fill( Jets[j].M(), Jets[j].Pt(), weight );
                        my_2d_histos["h_ISRJets_MassVsPt_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets Mass");
                        my_2d_histos["h_ISRJets_MassVsPt_"+cutVar.first]->GetYaxis()->SetTitle("ISRJets Pt");
                        my_2d_histos["h_ISRJets_MassVsChargedHadronEnergyFraction_"+cutVar.first]->Fill( Jets[j].M(), Jets_chargedHadronEnergyFraction[j], weight );
                        my_2d_histos["h_ISRJets_MassVsChargedHadronEnergyFraction_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets Mass");
                        my_2d_histos["h_ISRJets_MassVsChargedHadronEnergyFraction_"+cutVar.first]->GetYaxis()->SetTitle("ISRJets ChargedHadronEnergyFraction");
                        my_2d_histos["h_ISRJets_MassVsNeutralHadronEnergyFraction_"+cutVar.first]->Fill( Jets[j].M(), Jets_neutralHadronEnergyFraction[j], weight );
                        my_2d_histos["h_ISRJets_MassVsNeutralHadronEnergyFraction_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets Mass");
                        my_2d_histos["h_ISRJets_MassVsNeutralHadronEnergyFraction_"+cutVar.first]->GetYaxis()->SetTitle("ISRJets NeutralHadronEnergyFraction");
                        my_2d_histos["h_ISRJets_MassVsChargedHadronMultiplicity_"+cutVar.first]->Fill( Jets[j].M(), Jets_chargedHadronMultiplicity[j], weight );
                        my_2d_histos["h_ISRJets_MassVsChargedHadronMultiplicity_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets Mass");
                        my_2d_histos["h_ISRJets_MassVsChargedHadronMultiplicity_"+cutVar.first]->GetYaxis()->SetTitle("ISRJets ChargedHadronMultiplicity");
                        my_2d_histos["h_ISRJets_MassVsNeutralHadronMultiplicity_"+cutVar.first]->Fill( Jets[j].M(), Jets_neutralHadronMultiplicity[j], weight );
                        my_2d_histos["h_ISRJets_MassVsNeutralHadronMultiplicity_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets Mass");
                        my_2d_histos["h_ISRJets_MassVsNeutralHadronMultiplicity_"+cutVar.first]->GetYaxis()->SetTitle("ISRJets NeutralHadronMultiplicity");                      
                        my_2d_histos["h_ISRJets_MassVsChargedMultiplicity_"+cutVar.first]->Fill( Jets[j].M(), Jets_chargedMultiplicity[j], weight );
                        my_2d_histos["h_ISRJets_MassVsChargedMultiplicity_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets Mass");
                        my_2d_histos["h_ISRJets_MassVsChargedMultiplicity_"+cutVar.first]->GetYaxis()->SetTitle("ISRJets ChargedMultiplicity");
                        my_2d_histos["h_ISRJets_MassVsNeutralMultiplicity_"+cutVar.first]->Fill( Jets[j].M(), Jets_neutralMultiplicity[j], weight );
                        my_2d_histos["h_ISRJets_MassVsNeutralMultiplicity_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets Mass");
                        my_2d_histos["h_ISRJets_MassVsNeutralMultiplicity_"+cutVar.first]->GetYaxis()->SetTitle("ISRJets NeutralMultiplicity");
                        my_2d_histos["h_ISRJets_MassVSqgLikelihood_"+cutVar.first]->Fill( Jets[j].M(), Jets_qgLikelihood[j], weight );
                        my_2d_histos["h_ISRJets_MassVSqgLikelihood_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets Mass");
                        my_2d_histos["h_ISRJets_MassVSqgLikelihood_"+cutVar.first]->GetYaxis()->SetTitle("ISRJets qgLikelihood");
                        my_2d_histos["h_ISRJets_MassVSbJetTagDeepCSVtotb_"+cutVar.first]->Fill( Jets[j].M(), Jets_bJetTagDeepCSVtotb[j], weight );
                        my_2d_histos["h_ISRJets_MassVSbJetTagDeepCSVtotb_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets Mass");
                        my_2d_histos["h_ISRJets_MassVSbJetTagDeepCSVtotb_"+cutVar.first]->GetYaxis()->SetTitle("ISRJets bJetTagDeepCSVtotb");
                        my_2d_histos["h_ISRJets_PtVsChargedHadronEnergyFraction_"+cutVar.first]->Fill( Jets[j].Pt(), Jets_chargedHadronEnergyFraction[j], weight );
                        my_2d_histos["h_ISRJets_PtVsChargedHadronEnergyFraction_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets Pt");
                        my_2d_histos["h_ISRJets_PtVsChargedHadronEnergyFraction_"+cutVar.first]->GetYaxis()->SetTitle("ISRJets ChargedHadronEnergyFraction");
                        my_2d_histos["h_ISRJets_PtVsNeutralHadronEnergyFraction_"+cutVar.first]->Fill( Jets[j].Pt(), Jets_neutralHadronEnergyFraction[j], weight );
                        my_2d_histos["h_ISRJets_PtVsNeutralHadronEnergyFraction_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets Pt");
                        my_2d_histos["h_ISRJets_PtVsNeutralHadronEnergyFraction_"+cutVar.first]->GetYaxis()->SetTitle("ISRJets NeutralHadronEnergyFraction");
                        my_2d_histos["h_ISRJets_PtVsChargedHadronMultiplicity_"+cutVar.first]->Fill( Jets[j].Pt(), Jets_chargedHadronMultiplicity[j], weight );
                        my_2d_histos["h_ISRJets_PtVsChargedHadronMultiplicity_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets Pt");
                        my_2d_histos["h_ISRJets_PtVsChargedHadronMultiplicity_"+cutVar.first]->GetYaxis()->SetTitle("ISRJets ChargedHadronMultiplicity");
                        my_2d_histos["h_ISRJets_PtVsNeutralHadronMultiplicity_"+cutVar.first]->Fill( Jets[j].Pt(), Jets_neutralHadronMultiplicity[j], weight );
                        my_2d_histos["h_ISRJets_PtVsNeutralHadronMultiplicity_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets Pt");
                        my_2d_histos["h_ISRJets_PtVsNeutralHadronMultiplicity_"+cutVar.first]->GetYaxis()->SetTitle("ISRJets NeutralHadronMultiplicity"); 
                        my_2d_histos["h_ISRJets_PtVsChargedMultiplicity_"+cutVar.first]->Fill( Jets[j].Pt(), Jets_chargedMultiplicity[j], weight );
                        my_2d_histos["h_ISRJets_PtVsChargedMultiplicity_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets Pt");
                        my_2d_histos["h_ISRJets_PtVsChargedMultiplicity_"+cutVar.first]->GetYaxis()->SetTitle("ISRJets ChargedMultiplicity");
                        my_2d_histos["h_ISRJets_PtVsNeutralMultiplicity_"+cutVar.first]->Fill( Jets[j].Pt(), Jets_neutralMultiplicity[j], weight );
                        my_2d_histos["h_ISRJets_PtVsNeutralMultiplicity_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets Pt");
                        my_2d_histos["h_ISRJets_PtVsNeutralMultiplicity_"+cutVar.first]->GetYaxis()->SetTitle("ISRJets NeutralMultiplicity");
                        my_2d_histos["h_ISRJets_PtVSqgLikelihood_"+cutVar.first]->Fill( Jets[j].Pt(), Jets_qgLikelihood[j], weight );
                        my_2d_histos["h_ISRJets_PtVSqgLikelihood_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets Pt");
                        my_2d_histos["h_ISRJets_PtVSqgLikelihood_"+cutVar.first]->GetYaxis()->SetTitle("ISRJets qgLikelihood");
                        my_2d_histos["h_ISRJets_PtVSbJetTagDeepCSVtotb_"+cutVar.first]->Fill( Jets[j].Pt(), Jets_bJetTagDeepCSVtotb[j], weight );
                        my_2d_histos["h_ISRJets_PtVSbJetTagDeepCSVtotb_"+cutVar.first]->GetXaxis()->SetTitle("ISRJets Pt");
                        my_2d_histos["h_ISRJets_PtVSbJetTagDeepCSVtotb_"+cutVar.first]->GetYaxis()->SetTitle("ISRJets bJetTagDeepCSVtotb");
                        nISRJets++;
                    }                    
 
                    // --------------------------------
                    // NonISR Jets - by using TreeMaker
                    // -------------------------------- 
                    if (NonISRmatched_dr_ptr[j])
                    {
                        // 1D
                        my_histos["h_NonISRJets_Mass_"+cutVar.first]->Fill( Jets[j].M(), weight );
                        my_histos["h_NonISRJets_Mass_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets Mass");
                        my_histos["h_NonISRJets_Mass_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_NonISRJets_Pt_"+cutVar.first]->Fill( Jets[j].Pt(), weight );
                        my_histos["h_NonISRJets_Pt_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets Pt");
                        my_histos["h_NonISRJets_Pt_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_NonISRJets_Phi_"+cutVar.first]->Fill( Jets[j].Phi(), weight );
                        my_histos["h_NonISRJets_Phi_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets #phi");
                        my_histos["h_NonISRJets_Phi_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_NonISRJets_Eta_"+cutVar.first]->Fill( Jets[j].Eta(), weight );
                        my_histos["h_NonISRJets_Eta_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets #eta");
                        my_histos["h_NonISRJets_Eta_"+cutVar.first]->GetYaxis()->SetTitle("Events"); 
                        my_histos["h_NonISRJets_axismajor_"+cutVar.first]->Fill( Jets_axismajor[j], weight );
                        my_histos["h_NonISRJets_axismajor_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets axismajor");
                        my_histos["h_NonISRJets_axismajor_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_NonISRJets_axisminor_"+cutVar.first]->Fill( Jets_axisminor[j], weight );
                        my_histos["h_NonISRJets_axisminor_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets axisminor");
                        my_histos["h_NonISRJets_axisminor_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_NonISRJets_chargedHadronEnergyFraction_"+cutVar.first]->Fill( Jets_chargedHadronEnergyFraction[j], weight );
                        my_histos["h_NonISRJets_chargedHadronEnergyFraction_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets Energy fraction");
                        my_histos["h_NonISRJets_chargedHadronEnergyFraction_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_NonISRJets_neutralHadronEnergyFraction_"+cutVar.first]->Fill( Jets_neutralHadronEnergyFraction[j], weight );
                        my_histos["h_NonISRJets_neutralHadronEnergyFraction_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets Energy fraction");
                        my_histos["h_NonISRJets_neutralHadronEnergyFraction_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_NonISRJets_ptD_"+cutVar.first]->Fill( Jets_ptD[j], weight );
                        my_histos["h_NonISRJets_ptD_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets ptD");
                        my_histos["h_NonISRJets_ptD_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_NonISRJets_chargedHadronMultiplicity_"+cutVar.first]->Fill( Jets_chargedHadronMultiplicity[j], weight );
                        my_histos["h_NonISRJets_chargedHadronMultiplicity_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets CH Multiplicity");
                        my_histos["h_NonISRJets_chargedHadronMultiplicity_"+cutVar.first]->GetYaxis()->SetTitle("Event");
                        my_histos["h_NonISRJets_neutralHadronMultiplicity_"+cutVar.first]->Fill( Jets_neutralHadronMultiplicity[j], weight );
                        my_histos["h_NonISRJets_neutralHadronMultiplicity_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets NH Multiplicity");
                        my_histos["h_NonISRJets_neutralHadronMultiplicity_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_NonISRJets_chargedMultiplicity_"+cutVar.first]->Fill( Jets_chargedMultiplicity[j], weight );
                        my_histos["h_NonISRJets_chargedMultiplicity_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets Charged Multiplicity");
                        my_histos["h_NonISRJets_chargedMultiplicity_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_NonISRJets_neutralMultiplicity_"+cutVar.first]->Fill( Jets_neutralMultiplicity[j], weight );
                        my_histos["h_NonISRJets_neutralMultiplicity_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets Neural Multiplicity");
                        my_histos["h_NonISRJets_neutralMultiplicity_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_NonISRJets_qgLikelihood_"+cutVar.first]->Fill( Jets_qgLikelihood[j], weight );
                        my_histos["h_NonISRJets_qgLikelihood_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets qgLikelihood");
                        my_histos["h_NonISRJets_qgLikelihood_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_NonISRJets_bJetTagDeepCSVtotb_"+cutVar.first]->Fill( Jets_bJetTagDeepCSVtotb[j], weight ); 
                        my_histos["h_NonISRJets_bJetTagDeepCSVtotb_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets bJetTagDeepCSVtotb");
                        my_histos["h_NonISRJets_bJetTagDeepCSVtotb_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        // 2D
                        my_2d_histos["h_NonISRJets_EtaVsPhi_"+cutVar.first]->Fill( Jets[j].Eta(), Jets[j].Phi(), weight );
                        my_2d_histos["h_NonISRJets_EtaVsPhi_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets #eta");
                        my_2d_histos["h_NonISRJets_EtaVsPhi_"+cutVar.first]->GetYaxis()->SetTitle("NonISRJets #phi");
                        my_2d_histos["h_NonISRJets_MassVsPt_"+cutVar.first]->Fill( Jets[j].M(), Jets[j].Pt(), weight );
                        my_2d_histos["h_NonISRJets_MassVsPt_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets Mass");
                        my_2d_histos["h_NonISRJets_MassVsPt_"+cutVar.first]->GetYaxis()->SetTitle("NonISRJets Pt");
                        my_2d_histos["h_NonISRJets_MassVsChargedHadronEnergyFraction_"+cutVar.first]->Fill( Jets[j].M(), Jets_chargedHadronEnergyFraction[j], weight );
                        my_2d_histos["h_NonISRJets_MassVsChargedHadronEnergyFraction_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets Mass");
                        my_2d_histos["h_NonISRJets_MassVsChargedHadronEnergyFraction_"+cutVar.first]->GetYaxis()->SetTitle("NonISRJets ChargedHadronEnergyFraction");
                        my_2d_histos["h_NonISRJets_MassVsNeutralHadronEnergyFraction_"+cutVar.first]->Fill( Jets[j].M(), Jets_neutralHadronEnergyFraction[j], weight );
                        my_2d_histos["h_NonISRJets_MassVsNeutralHadronEnergyFraction_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets Mass");
                        my_2d_histos["h_NonISRJets_MassVsNeutralHadronEnergyFraction_"+cutVar.first]->GetYaxis()->SetTitle("NonISRJets NeutralHadronEnergyFraction");
                        my_2d_histos["h_NonISRJets_MassVsChargedHadronMultiplicity_"+cutVar.first]->Fill( Jets[j].M(), Jets_chargedHadronMultiplicity[j], weight );
                        my_2d_histos["h_NonISRJets_MassVsChargedHadronMultiplicity_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets Mass");
                        my_2d_histos["h_NonISRJets_MassVsChargedHadronMultiplicity_"+cutVar.first]->GetYaxis()->SetTitle("NonISRJets ChargedHadronMultiplicity");
                        my_2d_histos["h_NonISRJets_MassVsNeutralHadronMultiplicity_"+cutVar.first]->Fill( Jets[j].M(), Jets_neutralHadronMultiplicity[j], weight );
                        my_2d_histos["h_NonISRJets_MassVsNeutralHadronMultiplicity_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets Mass");
                        my_2d_histos["h_NonISRJets_MassVsNeutralHadronMultiplicity_"+cutVar.first]->GetYaxis()->SetTitle("NonISRJets NeutralHadronMultiplicity");
                        my_2d_histos["h_NonISRJets_MassVsChargedMultiplicity_"+cutVar.first]->Fill( Jets[j].M(), Jets_chargedMultiplicity[j], weight );
                        my_2d_histos["h_NonISRJets_MassVsChargedMultiplicity_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets Mass");
                        my_2d_histos["h_NonISRJets_MassVsChargedMultiplicity_"+cutVar.first]->GetYaxis()->SetTitle("NonISRJets ChargedMultiplicity");
                        my_2d_histos["h_NonISRJets_MassVsNeutralMultiplicity_"+cutVar.first]->Fill( Jets[j].M(), Jets_neutralMultiplicity[j], weight );
                        my_2d_histos["h_NonISRJets_MassVsNeutralMultiplicity_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets Mass");
                        my_2d_histos["h_NonISRJets_MassVsNeutralMultiplicity_"+cutVar.first]->GetYaxis()->SetTitle("NonISRJets NeutralMultiplicity");
                        my_2d_histos["h_NonISRJets_MassVSqgLikelihood_"+cutVar.first]->Fill( Jets[j].M(), Jets_qgLikelihood[j], weight );
                        my_2d_histos["h_NonISRJets_MassVSqgLikelihood_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets Mass");
                        my_2d_histos["h_NonISRJets_MassVSqgLikelihood_"+cutVar.first]->GetYaxis()->SetTitle("NonISRJets qgLikelihood");
                        my_2d_histos["h_NonISRJets_MassVSbJetTagDeepCSVtotb_"+cutVar.first]->Fill( Jets[j].M(), Jets_bJetTagDeepCSVtotb[j], weight );
                        my_2d_histos["h_NonISRJets_MassVSbJetTagDeepCSVtotb_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets Mass");
                        my_2d_histos["h_NonISRJets_MassVSbJetTagDeepCSVtotb_"+cutVar.first]->GetYaxis()->SetTitle("NonISRJets bJetTagDeepCSVtotb");
                        my_2d_histos["h_NonISRJets_PtVsChargedHadronEnergyFraction_"+cutVar.first]->Fill( Jets[j].Pt(), Jets_chargedHadronEnergyFraction[j], weight );
                        my_2d_histos["h_NonISRJets_PtVsChargedHadronEnergyFraction_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets Pt");
                        my_2d_histos["h_NonISRJets_PtVsChargedHadronEnergyFraction_"+cutVar.first]->GetYaxis()->SetTitle("NonISRJets ChargedHadronEnergyFraction");
                        my_2d_histos["h_NonISRJets_PtVsNeutralHadronEnergyFraction_"+cutVar.first]->Fill( Jets[j].Pt(), Jets_neutralHadronEnergyFraction[j], weight );
                        my_2d_histos["h_NonISRJets_PtVsNeutralHadronEnergyFraction_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets Pt");
                        my_2d_histos["h_NonISRJets_PtVsNeutralHadronEnergyFraction_"+cutVar.first]->GetYaxis()->SetTitle("NonISRJets NeutralHadronEnergyFraction");
                        my_2d_histos["h_NonISRJets_PtVsChargedHadronMultiplicity_"+cutVar.first]->Fill( Jets[j].Pt(), Jets_chargedHadronMultiplicity[j], weight );
                        my_2d_histos["h_NonISRJets_PtVsChargedHadronMultiplicity_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets Pt");
                        my_2d_histos["h_NonISRJets_PtVsChargedHadronMultiplicity_"+cutVar.first]->GetYaxis()->SetTitle("NonISRJets ChargedHadronMultiplicity");
                        my_2d_histos["h_NonISRJets_PtVsNeutralHadronMultiplicity_"+cutVar.first]->Fill( Jets[j].Pt(), Jets_neutralHadronMultiplicity[j], weight );
                        my_2d_histos["h_NonISRJets_PtVsNeutralHadronMultiplicity_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets Pt");
                        my_2d_histos["h_NonISRJets_PtVsNeutralHadronMultiplicity_"+cutVar.first]->GetYaxis()->SetTitle("NonISRJets NeutralHadronMultiplicity");
                        my_2d_histos["h_NonISRJets_PtVsChargedMultiplicity_"+cutVar.first]->Fill( Jets[j].Pt(), Jets_chargedMultiplicity[j], weight );
                        my_2d_histos["h_NonISRJets_PtVsChargedMultiplicity_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets Pt");
                        my_2d_histos["h_NonISRJets_PtVsChargedMultiplicity_"+cutVar.first]->GetYaxis()->SetTitle("NonISRJets ChargedMultiplicity");
                        my_2d_histos["h_NonISRJets_PtVsNeutralMultiplicity_"+cutVar.first]->Fill( Jets[j].Pt(), Jets_neutralMultiplicity[j], weight );
                        my_2d_histos["h_NonISRJets_PtVsNeutralMultiplicity_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets Pt");
                        my_2d_histos["h_NonISRJets_PtVsNeutralMultiplicity_"+cutVar.first]->GetYaxis()->SetTitle("NonISRJets NeutralMultiplicity");
                        my_2d_histos["h_NonISRJets_PtVSqgLikelihood_"+cutVar.first]->Fill( Jets[j].Pt(), Jets_qgLikelihood[j], weight );
                        my_2d_histos["h_NonISRJets_PtVSqgLikelihood_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets Pt");
                        my_2d_histos["h_NonISRJets_PtVSqgLikelihood_"+cutVar.first]->GetYaxis()->SetTitle("NonISRJets qgLikelihood");
                        my_2d_histos["h_NonISRJets_PtVSbJetTagDeepCSVtotb_"+cutVar.first]->Fill( Jets[j].Pt(), Jets_bJetTagDeepCSVtotb[j], weight );
                        my_2d_histos["h_NonISRJets_PtVSbJetTagDeepCSVtotb_"+cutVar.first]->GetXaxis()->SetTitle("NonISRJets Pt");
                        my_2d_histos["h_NonISRJets_PtVSbJetTagDeepCSVtotb_"+cutVar.first]->GetYaxis()->SetTitle("NonISRJets bJetTagDeepCSVtotb");
                        nNonISRJets++;
                    }

                    // ---------------------------------
                    // Other Jets - not ISR + not NonISR
                    // ---------------------------------
                    //if (OtherJets[j])
                    if ( (!(ISRmatched_dr_ptr[j])) && (!(NonISRmatched_dr_ptr[j])) )
                    {
                        // 1D 
                        my_histos["h_OtherJets_Mass_"+cutVar.first]->Fill( Jets[j].M(), weight );
                        my_histos["h_OtherJets_Mass_"+cutVar.first]->GetXaxis()->SetTitle("OtherJets Mass");
                        my_histos["h_OtherJets_Mass_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_OtherJets_Pt_"+cutVar.first]->Fill( Jets[j].Pt(), weight );
                        my_histos["h_OtherJets_Pt_"+cutVar.first]->GetXaxis()->SetTitle("OtherJets Pt");
                        my_histos["h_OtherJets_Pt_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_OtherJets_Phi_"+cutVar.first]->Fill( Jets[j].Phi(), weight );
                        my_histos["h_OtherJets_Phi_"+cutVar.first]->GetXaxis()->SetTitle("OtherJets #phi");
                        my_histos["h_OtherJets_Phi_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_OtherJets_Eta_"+cutVar.first]->Fill( Jets[j].Eta(), weight );
                        my_histos["h_OtherJets_Eta_"+cutVar.first]->GetXaxis()->SetTitle("OtherJets #eta");
                        my_histos["h_OtherJets_axismajor_"+cutVar.first]->Fill( Jets_axismajor[j], weight );
                        my_histos["h_OtherJets_axismajor_"+cutVar.first]->GetXaxis()->SetTitle("OtherJets axismajor");
                        my_histos["h_OtherJets_axismajor_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_OtherJets_axisminor_"+cutVar.first]->Fill( Jets_axisminor[j], weight );
                        my_histos["h_OtherJets_axisminor_"+cutVar.first]->GetXaxis()->SetTitle("OtherJets axisminor");
                        my_histos["h_OtherJets_axisminor_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_OtherJets_chargedHadronEnergyFraction_"+cutVar.first]->Fill( Jets_chargedHadronEnergyFraction[j], weight );
                        my_histos["h_OtherJets_chargedHadronEnergyFraction_"+cutVar.first]->GetXaxis()->SetTitle("OtherJets Energy fraction");
                        my_histos["h_OtherJets_chargedHadronEnergyFraction_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_OtherJets_neutralHadronEnergyFraction_"+cutVar.first]->Fill( Jets_neutralHadronEnergyFraction[j], weight );
                        my_histos["h_OtherJets_neutralHadronEnergyFraction_"+cutVar.first]->GetXaxis()->SetTitle("OtherJets Energy fraction");
                        my_histos["h_OtherJets_neutralHadronEnergyFraction_"+cutVar.first]->GetYaxis()->SetTitle("Events");       
                        my_histos["h_OtherJets_ptD_"+cutVar.first]->Fill( Jets_ptD[j], weight );
                        my_histos["h_OtherJets_ptD_"+cutVar.first]->GetXaxis()->SetTitle("OtherJets ptD");
                        my_histos["h_OtherJets_ptD_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_OtherJets_chargedHadronMultiplicity_"+cutVar.first]->Fill( Jets_chargedHadronMultiplicity[j], weight );
                        my_histos["h_OtherJets_chargedHadronMultiplicity_"+cutVar.first]->GetXaxis()->SetTitle("OtherJets CH Multiplicity");
                        my_histos["h_OtherJets_chargedHadronMultiplicity_"+cutVar.first]->GetYaxis()->SetTitle("Event");
                        my_histos["h_OtherJets_neutralHadronMultiplicity_"+cutVar.first]->Fill( Jets_neutralHadronMultiplicity[j], weight );
                        my_histos["h_OtherJets_neutralHadronMultiplicity_"+cutVar.first]->GetXaxis()->SetTitle("OtherJets NH Multiplicity");
                        my_histos["h_OtherJets_neutralHadronMultiplicity_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_OtherJets_chargedMultiplicity_"+cutVar.first]->Fill( Jets_chargedMultiplicity[j], weight );
                        my_histos["h_OtherJets_chargedMultiplicity_"+cutVar.first]->GetXaxis()->SetTitle("OtherJets Charged Multiplicity");
                        my_histos["h_OtherJets_chargedMultiplicity_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_OtherJets_neutralMultiplicity_"+cutVar.first]->Fill( Jets_neutralMultiplicity[j], weight );
                        my_histos["h_OtherJets_neutralMultiplicity_"+cutVar.first]->GetXaxis()->SetTitle("OtherJets Neural Multiplicity");
                        my_histos["h_OtherJets_neutralMultiplicity_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                        my_histos["h_OtherJets_qgLikelihood_"+cutVar.first]->Fill( Jets_qgLikelihood[j], weight );
                        my_histos["h_OtherJets_qgLikelihood_"+cutVar.first]->GetXaxis()->SetTitle("OtherJets qgLikelihood");
                        my_histos["h_OtherJets_qgLikelihood_"+cutVar.first]->GetYaxis()->SetTitle("Events");  
                        my_histos["h_OtherJets_bJetTagDeepCSVtotb_"+cutVar.first]->Fill( Jets_bJetTagDeepCSVtotb[j], weight ); 
                        my_histos["h_OtherJets_bJetTagDeepCSVtotb_"+cutVar.first]->GetXaxis()->SetTitle("OtherJets bJetTagDeepCSVtotb");
                        my_histos["h_OtherJets_bJetTagDeepCSVtotb_"+cutVar.first]->GetYaxis()->SetTitle("Events"); 
                        // 2D
                        my_2d_histos["h_OtherJets_EtaVsPhi_"+cutVar.first]->Fill( Jets[j].Eta(), Jets[j].Phi(), weight );
                        my_2d_histos["h_OtherJets_EtaVsPhi_"+cutVar.first]->GetXaxis()->SetTitle("OtherJets #eta");
                        my_2d_histos["h_OtherJets_EtaVsPhi_"+cutVar.first]->GetYaxis()->SetTitle("OtherJets #phi");
                        my_2d_histos["h_OtherJets_MassVsPt_"+cutVar.first]->Fill( Jets[j].M(), Jets[j].Pt(), weight );
                        my_2d_histos["h_OtherJets_MassVsPt_"+cutVar.first]->GetXaxis()->SetTitle("OtherJets Mass");
                        my_2d_histos["h_OtherJets_MassVsPt_"+cutVar.first]->GetYaxis()->SetTitle("OtherJets Pt");
                        my_2d_histos["h_OtherJets_MassVsChargedHadronEnergyFraction_"+cutVar.first]->Fill( Jets[j].M(), Jets_chargedHadronEnergyFraction[j], weight );
                        my_2d_histos["h_OtherJets_MassVsChargedHadronEnergyFraction_"+cutVar.first]->GetXaxis()->SetTitle("OtherJets Mass");
                        my_2d_histos["h_OtherJets_MassVsChargedHadronEnergyFraction_"+cutVar.first]->GetYaxis()->SetTitle("OtherJets ChargedHadronEnergyFraction");
                        my_2d_histos["h_OtherJets_MassVsNeutralHadronEnergyFraction_"+cutVar.first]->Fill( Jets[j].M(), Jets_neutralHadronEnergyFraction[j], weight );
                        my_2d_histos["h_OtherJets_MassVsNeutralHadronEnergyFraction_"+cutVar.first]->GetXaxis()->SetTitle("OtherJets Mass");
                        my_2d_histos["h_OtherJets_MassVsNeutralHadronEnergyFraction_"+cutVar.first]->GetYaxis()->SetTitle("OtherJets NeutralHadronEnergyFraction");
                        my_2d_histos["h_OtherJets_MassVsChargedHadronMultiplicity_"+cutVar.first]->Fill( Jets[j].M(), Jets_chargedHadronMultiplicity[j], weight );
                        my_2d_histos["h_OtherJets_MassVsChargedHadronMultiplicity_"+cutVar.first]->GetXaxis()->SetTitle("OtherJets Mass");
                        my_2d_histos["h_OtherJets_MassVsChargedHadronMultiplicity_"+cutVar.first]->GetYaxis()->SetTitle("OtherJets ChargedHadronMultiplicity");
                        my_2d_histos["h_OtherJets_MassVsNeutralHadronMultiplicity_"+cutVar.first]->Fill( Jets[j].M(), Jets_neutralHadronMultiplicity[j], weight );
                        my_2d_histos["h_OtherJets_MassVsNeutralHadronMultiplicity_"+cutVar.first]->GetXaxis()->SetTitle("OtherJets Mass");
                        my_2d_histos["h_OtherJets_MassVsNeutralHadronMultiplicity_"+cutVar.first]->GetYaxis()->SetTitle("OtherJets NeutralHadronMultiplicity");
                        my_2d_histos["h_OtherJets_MassVsChargedMultiplicity_"+cutVar.first]->Fill( Jets[j].M(), Jets_chargedMultiplicity[j], weight );
                        my_2d_histos["h_OtherJets_MassVsChargedMultiplicity_"+cutVar.first]->GetXaxis()->SetTitle("OtherJets Mass");
                        my_2d_histos["h_OtherJets_MassVsChargedMultiplicity_"+cutVar.first]->GetYaxis()->SetTitle("OtherJets ChargedMultiplicity");
                        my_2d_histos["h_OtherJets_MassVsNeutralMultiplicity_"+cutVar.first]->Fill( Jets[j].M(), Jets_neutralMultiplicity[j], weight );
                        my_2d_histos["h_OtherJets_MassVsNeutralMultiplicity_"+cutVar.first]->GetXaxis()->SetTitle("OtherJets Mass");
                        my_2d_histos["h_OtherJets_MassVsNeutralMultiplicity_"+cutVar.first]->GetYaxis()->SetTitle("OtherJets NeutralMultiplicity");
                        my_2d_histos["h_OtherJets_MassVSqgLikelihood_"+cutVar.first]->Fill( Jets[j].M(), Jets_qgLikelihood[j], weight );
                        my_2d_histos["h_OtherJets_MassVSqgLikelihood_"+cutVar.first]->GetXaxis()->SetTitle("OtherJets Mass");
                        my_2d_histos["h_OtherJets_MassVSqgLikelihood_"+cutVar.first]->GetYaxis()->SetTitle("OtherJets qgLikelihood");
                        my_2d_histos["h_OtherJets_MassVSbJetTagDeepCSVtotb_"+cutVar.first]->Fill( Jets[j].M(), Jets_bJetTagDeepCSVtotb[j], weight );
                        my_2d_histos["h_OtherJets_MassVSbJetTagDeepCSVtotb_"+cutVar.first]->GetXaxis()->SetTitle("OtherJets Mass");
                        my_2d_histos["h_OtherJets_MassVSbJetTagDeepCSVtotb_"+cutVar.first]->GetYaxis()->SetTitle("OtherJets bJetTagDeepCSVtotb");
                        my_2d_histos["h_OtherJets_PtVsChargedHadronEnergyFraction_"+cutVar.first]->Fill( Jets[j].Pt(), Jets_chargedHadronEnergyFraction[j], weight );
                        my_2d_histos["h_OtherJets_PtVsChargedHadronEnergyFraction_"+cutVar.first]->GetXaxis()->SetTitle("OtherJets Pt");
                        my_2d_histos["h_OtherJets_PtVsChargedHadronEnergyFraction_"+cutVar.first]->GetYaxis()->SetTitle("OtherJets ChargedHadronEnergyFraction");
                        my_2d_histos["h_OtherJets_PtVsNeutralHadronEnergyFraction_"+cutVar.first]->Fill( Jets[j].Pt(), Jets_neutralHadronEnergyFraction[j], weight );
                        my_2d_histos["h_OtherJets_PtVsNeutralHadronEnergyFraction_"+cutVar.first]->GetXaxis()->SetTitle("OtherJets Pt");
                        my_2d_histos["h_OtherJets_PtVsNeutralHadronEnergyFraction_"+cutVar.first]->GetYaxis()->SetTitle("OtherJets NeutralHadronEnergyFraction");
                        my_2d_histos["h_OtherJets_PtVsChargedHadronMultiplicity_"+cutVar.first]->Fill( Jets[j].Pt(), Jets_chargedHadronMultiplicity[j], weight );
                        my_2d_histos["h_OtherJets_PtVsChargedHadronMultiplicity_"+cutVar.first]->GetXaxis()->SetTitle("OtherJets Pt");
                        my_2d_histos["h_OtherJets_PtVsChargedHadronMultiplicity_"+cutVar.first]->GetYaxis()->SetTitle("OtherJets ChargedHadronMultiplicity");
                        my_2d_histos["h_OtherJets_PtVsNeutralHadronMultiplicity_"+cutVar.first]->Fill( Jets[j].Pt(), Jets_neutralHadronMultiplicity[j], weight );
                        my_2d_histos["h_OtherJets_PtVsNeutralHadronMultiplicity_"+cutVar.first]->GetXaxis()->SetTitle("OtherJets Pt");
                        my_2d_histos["h_OtherJets_PtVsNeutralHadronMultiplicity_"+cutVar.first]->GetYaxis()->SetTitle("OtherJets NeutralHadronMultiplicity");
                        my_2d_histos["h_OtherJets_PtVsChargedMultiplicity_"+cutVar.first]->Fill( Jets[j].Pt(), Jets_chargedMultiplicity[j], weight );
                        my_2d_histos["h_OtherJets_PtVsChargedMultiplicity_"+cutVar.first]->GetXaxis()->SetTitle("OtherJets Pt");
                        my_2d_histos["h_OtherJets_PtVsChargedMultiplicity_"+cutVar.first]->GetYaxis()->SetTitle("OtherJets ChargedMultiplicity");
                        my_2d_histos["h_OtherJets_PtVsNeutralMultiplicity_"+cutVar.first]->Fill( Jets[j].Pt(), Jets_neutralMultiplicity[j], weight );
                        my_2d_histos["h_OtherJets_PtVsNeutralMultiplicity_"+cutVar.first]->GetXaxis()->SetTitle("OtherJets Pt");
                        my_2d_histos["h_OtherJets_PtVsNeutralMultiplicity_"+cutVar.first]->GetYaxis()->SetTitle("OtherJets NeutralMultiplicity");  
                        my_2d_histos["h_OtherJets_PtVSqgLikelihood_"+cutVar.first]->Fill( Jets[j].Pt(), Jets_qgLikelihood[j], weight );
                        my_2d_histos["h_OtherJets_PtVSqgLikelihood_"+cutVar.first]->GetXaxis()->SetTitle("OtherJets Pt");
                        my_2d_histos["h_OtherJets_PtVSqgLikelihood_"+cutVar.first]->GetYaxis()->SetTitle("OtherJets qgLikelihood");
                        my_2d_histos["h_OtherJets_PtVSbJetTagDeepCSVtotb_"+cutVar.first]->Fill( Jets[j].Pt(), Jets_bJetTagDeepCSVtotb[j], weight );
                        my_2d_histos["h_OtherJets_PtVSbJetTagDeepCSVtotb_"+cutVar.first]->GetXaxis()->SetTitle("OtherJets Pt");
                        my_2d_histos["h_OtherJets_PtVSbJetTagDeepCSVtotb_"+cutVar.first]->GetYaxis()->SetTitle("OtherJets bJetTagDeepCSVtotb"); 
                        nOtherJets++;
                    }
                }
                my_histos["h_nISRJets_"+cutVar.first]->Fill( nISRJets, weight );
                my_histos["h_nISRJets_"+cutVar.first]->GetXaxis()->SetTitle("nISRJets");
                my_histos["h_nISRJets_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                my_histos["h_nNonISRJets_"+cutVar.first]->Fill( nNonISRJets, weight );
                my_histos["h_nNonISRJets_"+cutVar.first]->GetXaxis()->SetTitle("nNonISRJets");
                my_histos["h_nNonISRJets_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                my_histos["h_nOtherJets_"+cutVar.first]->Fill( nOtherJets, weight );
                my_histos["h_nOtherJets_"+cutVar.first]->GetXaxis()->SetTitle("nOtherJets");
                my_histos["h_nOtherJets_"+cutVar.first]->GetYaxis()->SetTitle("Events");

                // -----------------------------------
                // DeltaEta between RecoISR and GenISR
                // -----------------------------------
                for (unsigned int i = 0; i < dEta_RecoISR_GenISR.size(); i++)
                {
                    my_histos["h_dEta_RecoISR_GenISR_"+cutVar.first]->Fill( dEta_RecoISR_GenISR.at(i), weight );     
                    my_histos["h_dEta_RecoISR_GenISR_"+cutVar.first]->GetXaxis()->SetTitle("#DeltaEta(RecoISR,GenISR)");
                    my_histos["h_dEta_RecoISR_GenISR_"+cutVar.first]->GetYaxis()->SetTitle("Events"); 
                }

                // --------------------------------------------
                // DeltaR between each ISR and each NonISR jets
                // -------------------------------------------- 
                for (unsigned int i = 0; i < dR_ISR_NonISR.size(); i++)
                {
                    my_histos["h_dR_ISR_NonISR_"+cutVar.first]->Fill( dR_ISR_NonISR.at(i), weight );           
                    my_histos["h_dR_ISR_NonISR_"+cutVar.first]->GetXaxis()->SetTitle("#DeltaR(ISR,NonISR)");
                    my_histos["h_dR_ISR_NonISR_"+cutVar.first]->GetYaxis()->SetTitle("Events"); 
                }

                // --------------------------------------------
                // DeltaR between closest top and each ISR jets
                // --------------------------------------------
                for (unsigned int i = 0; i < dR_top_ISR.size(); i++)
                {
                    my_histos["h_dR_top_ISR_"+cutVar.first]->Fill( dR_top_ISR.at(i), weight );
                    my_histos["h_dR_top_ISR_"+cutVar.first]->GetXaxis()->SetTitle("#DeltaR(top,ISR)");
                    my_histos["h_dR_top_ISR_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                }
           
                // -----------------------------------------------
                // DeltaR between closest top and each NonISR jets
                // -----------------------------------------------
                for (unsigned int i = 0; i < dR_top_NonISR.size(); i++)
                {   
                    my_histos["h_dR_top_NonISR_"+cutVar.first]->Fill( dR_top_NonISR.at(i), weight );
                    my_histos["h_dR_top_NonISR_"+cutVar.first]->GetXaxis()->SetTitle("#DeltaR(top,NonISR)");
                    my_histos["h_dR_top_NonISR_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                }   

                // ---------------------------------------------
                // DeltaR between closest bjet and each ISR jets
                // ---------------------------------------------
                for (unsigned int i = 0; i < dR_bjet_ISR.size(); i++)
                {
                    my_histos["h_dR_bjet_ISR_"+cutVar.first]->Fill( dR_bjet_ISR.at(i), weight );
                    my_histos["h_dR_bjet_ISR_"+cutVar.first]->GetXaxis()->SetTitle("#DeltaR(bjet,ISR)");
                    my_histos["h_dR_bjet_ISR_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                }                 

                // ------------------------------------------------
                // DeltaR between closest bjet and each NonISR jets
                // ------------------------------------------------                
                for (unsigned int i = 0; i < dR_bjet_NonISR.size(); i++)  
                {
                    my_histos["h_dR_bjet_NonISR_"+cutVar.first]->Fill( dR_bjet_NonISR.at(i), weight );
                    my_histos["h_dR_bjet_NonISR_"+cutVar.first]->GetXaxis()->SetTitle("#DeltaR(bjet,NonISR)");
                    my_histos["h_dR_bjet_NonISR_"+cutVar.first]->GetYaxis()->SetTitle("Events");
                }  
            }
        }
    }    
}

void ISRJets_Analyzer::WriteHistos(TFile* outfile)
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
