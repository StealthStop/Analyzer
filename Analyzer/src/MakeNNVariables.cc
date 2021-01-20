#define MakeNNVariables_cxx
#include "Analyzer/Analyzer/include/MakeNNVariables.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/MiniTupleMaker.h"
#include "Framework/Framework/include/Utility.h" 

#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TRandom3.h>
#include <iostream>
#include <stdio.h> 
//#include <fstream>
//#include <cstdio>

MakeNNVariables::MakeNNVariables()
{
    InitHistos();
}

void MakeNNVariables::InitHistos()
{
    my_histos.emplace( "EventCounterTrain_0l", std::make_shared<TH1D>( "EventCounterTrain_0l", "EventCounterTrain_0l", 2, -1.1, 1.1 ) ); 
    my_histos.emplace( "EventCounterTrain_1l", std::make_shared<TH1D>( "EventCounterTrain_1l", "EventCounterTrain_1l", 2, -1.1, 1.1 ) );
    my_histos.emplace( "EventCounterTrain_2l", std::make_shared<TH1D>( "EventCounterTrain_2l", "EventCounterTrain_2l", 2, -1.1, 1.1 ) );
    my_histos.emplace( "EventCounterTest_0l",  std::make_shared<TH1D>( "EventCounterTest_0l",  "EventCounterTest_0l",  2, -1.1, 1.1 ) ); 
    my_histos.emplace( "EventCounterTest_1l",  std::make_shared<TH1D>( "EventCounterTest_1l",  "EventCounterTest_1l",  2, -1.1, 1.1 ) );
    my_histos.emplace( "EventCounterTest_2l",  std::make_shared<TH1D>( "EventCounterTest_2l",  "EventCounterTest_2l",  2, -1.1, 1.1 ) );
    my_histos.emplace( "EventCounterVal_0l",   std::make_shared<TH1D>( "EventCounterVal_0l",   "EventCounterVal_0l",   2, -1.1, 1.1 ) ); 
    my_histos.emplace( "EventCounterVal_1l",   std::make_shared<TH1D>( "EventCounterVal_1l",   "EventCounterVal_1l",   2, -1.1, 1.1 ) ); 
    my_histos.emplace( "EventCounterVal_2l",   std::make_shared<TH1D>( "EventCounterVal_2l",   "EventCounterVal_2l",   2, -1.1, 1.1 ) ); 

}//END of init histos

void MakeNNVariables::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    int count_0l = 0, numPassTrain_0l = 0, numPassTrain_1l = 0, numPassTrain_2l = 0; 
    int count_1l = 0, numPassTest_0l  = 0, numPassTest_1l  = 0, numPassTest_2l  = 0; 
    int count_2l = 0, numPassVal_0l   = 0, numPassVal_1l   = 0, numPassVal_2l   = 0;

    while( tr.getNextEvent() )
    {
        const auto& eventCounter        = tr.getVar<int>("eventCounter");
        const auto& isSignal            = tr.getVar<bool>("isSignal");
        const auto& passBaseline0l      = tr.getVar<bool>("passBaseline0l_Good"); 
        const auto& passBaseline1l      = tr.getVar<bool>("passBaseline1l_Good");
        const auto& passBaseline2l_pt20 = tr.getVar<bool>("passBaseline2l_pt20");
        const auto& filetag             = tr.getVar<std::string>("filetag");

        auto& mass = tr.createDerivedVar<double>("mass", 0.0);
        if(!isSignal)
        {
            mass = 173.0;
        }
        else
        {
            for(unsigned int m = 300; m < 1500; m+=50)
            {
                mass = (filetag.find(std::to_string(m)) != std::string::npos) ? m : mass;
            }
        }
       
        //------------------------------------
        //-- Print Event Number
        //------------------------------------
        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() % 1000 == 0 ) printf( " Event %i\n", tr.getEvtNum() );        

        //-----------------------------------
        //  Initialize the tree
        //-----------------------------------       
        std::set<std::string> varGeneral = 
        {
            "Lumi",   
            "mass",   
            "isSignal",
            "Weight",  
            "totalEventWeight",
            "HT_trigger_pt30",
            "HT_trigger_pt45",
            "NGoodJets_pt20_double",
            "NGoodJets_pt30_double",
            "NGoodJets_pt45_double",
            "NGoodBJets_pt30_double",
            "NGoodBJets_pt45_double",          
        };

        std::set<std::string> varEventShape = 
        {
            "fwm2_top6",    "fwm3_top6",    "fwm4_top6",   "fwm5_top6",
            "fwm6_top6",    "fwm7_top6",    "fwm8_top6",   "fwm9_top6", "fwm10_top6",
            "jmt_ev0_top6", "jmt_ev1_top6", "jmt_ev2_top6",

        };

        std::set<std::string> varJets = 
        {
            "Jet_m_1",            "Jet_m_2",            "Jet_m_3",            "Jet_m_4",             "Jet_m_5",             "Jet_m_6",
            "Jet_m_7",            "Jet_m_8",            "Jet_m_9",            "Jet_m_10",            "Jet_m_11",            "Jet_m_12",
            "Jet_eta_1",          "Jet_eta_2",          "Jet_eta_3",          "Jet_eta_4",           "Jet_eta_5",           "Jet_eta_6",  
            "Jet_eta_7",          "Jet_eta_8",          "Jet_eta_9",          "Jet_eta_10",          "Jet_eta_11",          "Jet_eta_12",
            "Jet_phi_1",          "Jet_phi_2",          "Jet_phi_3",          "Jet_phi_4",           "Jet_phi_5",           "Jet_phi_6",  
            "Jet_phi_7",          "Jet_phi_8",          "Jet_phi_9",          "Jet_phi_10",          "Jet_phi_11",          "Jet_phi_12",
            "Jet_pt_1",           "Jet_pt_2",           "Jet_pt_3",           "Jet_pt_4",            "Jet_pt_5",            "Jet_pt_6", 
            "Jet_pt_7",           "Jet_pt_8",           "Jet_pt_9",           "Jet_pt_10",           "Jet_pt_11",           "Jet_pt_12",
            "Jet_flavb_1",        "Jet_flavb_2",        "Jet_flavb_3",        "Jet_flavb_4",         "Jet_flavb_5",         "Jet_flavb_6", 
            "Jet_flavb_7",        "Jet_flavb_8",        "Jet_flavb_9",        "Jet_flavb_10",        "Jet_flavb_11",        "Jet_flavb_12",
            "Jet_flavg_1",        "Jet_flavg_2",        "Jet_flavg_3",        "Jet_flavg_4",         "Jet_flavg_5",         "Jet_flavg_6", 
            "Jet_flavg_7",        "Jet_flavg_8",        "Jet_flavg_9",        "Jet_flavg_10",        "Jet_flavg_11",        "Jet_flavg_12",
            "Jet_flavq_1",        "Jet_flavq_2",        "Jet_flavq_3",        "Jet_flavq_4",         "Jet_flavq_5",         "Jet_flavq_6", 
            "Jet_flavq_7",        "Jet_flavq_8",        "Jet_flavq_9",        "Jet_flavq_10",        "Jet_flavq_11",        "Jet_flavq_12",
            "Jet_ptD_1",          "Jet_ptD_2",          "Jet_ptD_3",          "Jet_ptD_4",           "Jet_ptD_5",           "Jet_ptD_6",  
            "Jet_ptD_7",          "Jet_ptD_8",          "Jet_ptD_9",          "Jet_ptD_10",          "Jet_ptD_11",          "Jet_ptD_12",
            "Jet_nEF_1",          "Jet_nEF_2",          "Jet_nEF_3",          "Jet_nEF_4",           "Jet_nEF_5",           "Jet_nEF_6",  
            "Jet_nEF_7",          "Jet_nEF_8",          "Jet_nEF_9",          "Jet_nEF_10",          "Jet_nEF_11",          "Jet_nEF_12",
            "Jet_cEF_1",          "Jet_cEF_2",          "Jet_cEF_3",          "Jet_cEF_4",           "Jet_cEF_5",           "Jet_cEF_6",  
            "Jet_cEF_7",          "Jet_cEF_8",          "Jet_cEF_9",          "Jet_cEF_10",          "Jet_cEF_11",          "Jet_cEF_12",
            "Jet_nHF_1",          "Jet_nHF_2",          "Jet_nHF_3",          "Jet_nHF_4",           "Jet_nHF_5",           "Jet_nHF_6",  
            "Jet_nHF_7",          "Jet_nHF_8",          "Jet_nHF_9",          "Jet_nHF_10",          "Jet_nHF_11",          "Jet_nHF_12",
            "Jet_cHF_1",          "Jet_cHF_2",          "Jet_cHF_3",          "Jet_cHF_4",           "Jet_cHF_5",           "Jet_cHF_6",  
            "Jet_cHF_7",          "Jet_cHF_8",          "Jet_cHF_9",          "Jet_cHF_10",          "Jet_cHF_11",          "Jet_cHF_12",
            "Jet_axismajor_1",    "Jet_axismajor_2",    "Jet_axismajor_3",    "Jet_axismajor_4",     "Jet_axismajor_5",     "Jet_axismajor_6", 
            "Jet_axismajor_7",    "Jet_axismajor_8",    "Jet_axismajor_9",    "Jet_axismajor_10",    "Jet_axismajor_11",    "Jet_axismajor_12",
            "Jet_axisminor_1",    "Jet_axisminor_2",    "Jet_axisminor_3",    "Jet_axisminor_4",     "Jet_axisminor_5",     "Jet_axisminor_6", 
            "Jet_axisminor_7",    "Jet_axisminor_8",    "Jet_axisminor_9",    "Jet_axisminor_10",    "Jet_axisminor_11",    "Jet_axisminor_12",
            "Jet_multiplicity_1", "Jet_multiplicity_2", "Jet_multiplicity_3", "Jet_multiplicity_4",  "Jet_multiplicity_5",  "Jet_multiplicity_6", 
            "Jet_multiplicity_7", "Jet_multiplicity_8", "Jet_multiplicity_9", "Jet_multiplicity_10", "Jet_multiplicity_11", "Jet_multiplicity_12",
        };            

        std::set<std::string> varJetsAK8 =
        {
            "JetsAK8_m_1",              "JetsAK8_m_2",              "JetsAK8_m_3",              "JetsAK8_m_4",              "JetsAK8_m_5",
            "JetsAK8_eta_1",            "JetsAK8_eta_2",            "JetsAK8_eta_3",            "JetsAK8_eta_4",            "JetsAK8_eta_5",
            "JetsAK8_phi_1",            "JetsAK8_phi_2",            "JetsAK8_phi_3",            "JetsAK8_phi_4",            "JetsAK8_phi_5",
            "JetsAK8_pt_1",             "JetsAK8_pt_2",             "JetsAK8_pt_3",             "JetsAK8_pt_4",             "JetsAK8_pt_5",            
            "JetsAK8_SDM_1",            "JetsAK8_SDM_2",            "JetsAK8_SDM_3",            "JetsAK8_SDM_4",            "JetsAK8_SDM_5",
            "JetsAK8_Pruned_1",         "JetsAK8_Pruned_2",         "JetsAK8_Pruned_3",         "JetsAK8_Pruned_4",         "JetsAK8_Pruned_5",
            "JetsAK8_Tau1_1",           "JetsAK8_Tau1_2",           "JetsAK8_Tau1_3",           "JetsAK8_Tau1_4",           "JetsAK8_Tau1_5",
            "JetsAK8_Tau2_1",           "JetsAK8_Tau2_2",           "JetsAK8_Tau2_3",           "JetsAK8_Tau2_4",           "JetsAK8_Tau2_5",
            "JetsAK8_Tau3_1",           "JetsAK8_Tau3_2",           "JetsAK8_Tau3_3",           "JetsAK8_Tau3_4",           "JetsAK8_Tau3_5",
            "JetsAK8_axismajor_1",      "JetsAK8_axismajor_2",      "JetsAK8_axismajor_3",      "JetsAK8_axismajor_4",      "JetsAK8_axismajor_5",
            "JetsAK8_axisminor_1",      "JetsAK8_axisminor_2",      "JetsAK8_axisminor_3",      "JetsAK8_axisminor_4",      "JetsAK8_axisminor_5",
            "JetsAK8_nsubjets_1",       "JetsAK8_nsubjets_2",       "JetsAK8_nsubjets_3",       "JetsAK8_nsubjets_4",       "JetsAK8_nsubjets_5",
            "JetsAK8_tDiscriminator_1", "JetsAK8_tDiscriminator_2", "JetsAK8_tDiscriminator_3", "JetsAK8_tDiscriminator_4", "JetsAK8_tDiscriminator_5",
            "JetsAK8_wDiscriminator_1", "JetsAK8_wDiscriminator_2", "JetsAK8_wDiscriminator_3", "JetsAK8_wDiscriminator_4", "JetsAK8_wDiscriminator_5",
            "JetsAK8_hDiscriminator_1", "JetsAK8_hDiscriminator_2", "JetsAK8_hDiscriminator_3", "JetsAK8_hDiscriminator_4", "JetsAK8_hDiscriminator_5",
            "JetsAK8_multiplicity_1",   "JetsAK8_multiplicity_2",   "JetsAK8_multiplicity_3",   "JetsAK8_multiplicity_4",   "JetsAK8_multiplicity_5",
        };

        std::set<std::string> varLeptonic =
        {
            "Mbl",
            "deepESM_val", // legacy from 1l analysis
            "lvMET_cm_m",
            "lvMET_cm_eta",
            "lvMET_cm_phi",
            "lvMET_cm_pt",
            "GoodLeptons_m_1",   "GoodLeptons_m_2",
            "GoodLeptons_eta_1", "GoodLeptons_eta_2",
            "GoodLeptons_phi_1", "GoodLeptons_phi_2",
            "GoodLeptons_pt_1",  "GoodLeptons_pt_2",
            "GoodLeptons_miniIso_1",  "GoodLeptons_miniIso_2",
        };

        std::set<std::string> varOldSeed =
        {            
            // without ISR cheating filter
            "Mass_stop1_PtRank_cm_OldSeed",                 "Mass_stop2_PtRank_cm_OldSeed", 
            "Pt_stop1_PtRank_cm_OldSeed",                   "Pt_stop2_PtRank_cm_OldSeed",
            "Phi_stop1_PtRank_cm_OldSeed",                  "Phi_stop2_PtRank_cm_OldSeed",  
            "Eta_stop1_PtRank_cm_OldSeed",                  "Eta_stop2_PtRank_cm_OldSeed",
            "Mass_stop1_MassRank_cm_OldSeed",               "Mass_stop2_MassRank_cm_OldSeed", 
            "Pt_stop1_MassRank_cm_OldSeed",                 "Pt_stop2_MassRank_cm_OldSeed",
            "Phi_stop1_MassRank_cm_OldSeed",                "Phi_stop2_MassRank_cm_OldSeed",  
            "Eta_stop1_MassRank_cm_OldSeed",                "Eta_stop2_MassRank_cm_OldSeed",
            "Mass_stop1_ScalarPtRank_cm_OldSeed",           "Mass_stop2_ScalarPtRank_cm_OldSeed", 
            "Pt_stop1_ScalarPtRank_cm_OldSeed",             "Pt_stop2_ScalarPtRank_cm_OldSeed",
            "Phi_stop1_ScalarPtRank_cm_OldSeed",            "Phi_stop2_ScalarPtRank_cm_OldSeed",  
            "Eta_stop1_ScalarPtRank_cm_OldSeed",            "Eta_stop2_ScalarPtRank_cm_OldSeed",
            // with ISR cheating filter
            //"Mass_stop1_PtRank_cm_OldSeed_maskedISR",       "Mass_stop2_PtRank_cm_OldSeed_maskedISR", 
            //"Pt_stop1_PtRank_cm_OldSeed_maskedISR",         "Pt_stop2_PtRank_cm_OldSeed_maskedISR",
            //"Phi_stop1_PtRank_cm_OldSeed_maskedISR",        "Phi_stop2_PtRank_cm_OldSeed_maskedISR",
            //"Eta_stop1_PtRank_cm_OldSeed_maskedISR",        "Eta_stop2_PtRank_cm_OldSeed_maskedISR",
            //"Mass_stop1_MassRank_cm_OldSeed_maskedISR",     "Mass_stop2_MassRank_cm_OldSeed_maskedISR",
            //"Pt_stop1_MassRank_cm_OldSeed_maskedISR",       "Pt_stop2_MassRank_cm_OldSeed_maskedISR",
            //"Phi_stop1_MassRank_cm_OldSeed_maskedISR",      "Phi_stop2_MassRank_cm_OldSeed_maskedISR",
            //"Eta_stop1_MassRank_cm_OldSeed_maskedISR",      "Eta_stop2_MassRank_cm_OldSeed_maskedISR",
            //"Mass_stop1_ScalarPtRank_cm_OldSeed_maskedISR", "Mass_stop2_ScalarPtRank_cm_OldSeed_maskedISR",
            //"Pt_stop1_ScalarPtRank_cm_OldSeed_maskedISR",   "Pt_stop2_ScalarPtRank_cm_OldSeed_maskedISR",
            //"Phi_stop1_ScalarPtRank_cm_OldSeed_maskedISR",  "Phi_stop2_ScalarPtRank_cm_OldSeed_maskedISR",
            //"Eta_stop1_ScalarPtRank_cm_OldSeed_maskedISR",  "Eta_stop2_ScalarPtRank_cm_OldSeed_maskedISR",
        };        

        std::set<std::string> varTopSeed =
        {
            // without ISR cheating filter - Only for 0l case
            "Mass_stop1_PtRank_cm_TopSeed",                 "Mass_stop2_PtRank_cm_TopSeed",
            "Pt_stop1_PtRank_cm_TopSeed",                   "Pt_stop2_PtRank_cm_TopSeed",
            "Phi_stop1_PtRank_cm_TopSeed",                  "Phi_stop2_PtRank_cm_TopSeed",  
            "Eta_stop1_PtRank_cm_TopSeed",                  "Eta_stop2_PtRank_cm_TopSeed",
            "Mass_stop1_MassRank_cm_TopSeed",               "Mass_stop2_MassRank_cm_TopSeed", 
            "Pt_stop1_MassRank_cm_TopSeed",                 "Pt_stop2_MassRank_cm_TopSeed",
            "Phi_stop1_MassRank_cm_TopSeed",                "Phi_stop2_MassRank_cm_TopSeed",  
            "Eta_stop1_MassRank_cm_TopSeed",                "Eta_stop2_MassRank_cm_TopSeed",
            "Mass_stop1_ScalarPtRank_cm_TopSeed",           "Mass_stop2_ScalarPtRank_cm_TopSeed", 
            "Pt_stop1_ScalarPtRank_cm_TopSeed",             "Pt_stop2_ScalarPtRank_cm_TopSeed",
            "Phi_stop1_ScalarPtRank_cm_TopSeed",            "Phi_stop2_ScalarPtRank_cm_TopSeed",  
            "Eta_stop1_ScalarPtRank_cm_TopSeed",            "Eta_stop2_ScalarPtRank_cm_TopSeed",
            // with ISR cheating filter - Only for 0l case
            //"Mass_stop1_PtRank_cm_TopSeed_maskedISR",       "Mass_stop2_PtRank_cm_TopSeed_maskedISR", 
            //"Pt_stop1_PtRank_cm_TopSeed_maskedISR",         "Pt_stop2_PtRank_cm_TopSeed_maskedISR",
            //"Phi_stop1_PtRank_cm_TopSeed_maskedISR",        "Phi_stop2_PtRank_cm_TopSeed_maskedISR",
            //"Eta_stop1_PtRank_cm_TopSeed_maskedISR",        "Eta_stop2_PtRank_cm_TopSeed_maskedISR",
            //"Mass_stop1_MassRank_cm_TopSeed_maskedISR",     "Mass_stop2_MassRank_cm_TopSeed_maskedISR",
            //"Pt_stop1_MassRank_cm_TopSeed_maskedISR",       "Pt_stop2_MassRank_cm_TopSeed_maskedISR",
            //"Phi_stop1_MassRank_cm_TopSeed_maskedISR",      "Phi_stop2_MassRank_cm_TopSeed_maskedISR",
            //"Eta_stop1_MassRank_cm_TopSeed_maskedISR",      "Eta_stop2_MassRank_cm_TopSeed_maskedISR",
            //"Mass_stop1_ScalarPtRank_cm_TopSeed_maskedISR", "Mass_stop2_ScalarPtRank_cm_TopSeed_maskedISR",
            //"Pt_stop1_ScalarPtRank_cm_TopSeed_maskedISR",   "Pt_stop2_ScalarPtRank_cm_TopSeed_maskedISR",
            //"Phi_stop1_ScalarPtRank_cm_TopSeed_maskedISR",  "Phi_stop2_ScalarPtRank_cm_TopSeed_maskedISR",
            //"Eta_stop1_ScalarPtRank_cm_TopSeed_maskedISR",  "Eta_stop2_ScalarPtRank_cm_TopSeed_maskedISR",            
        };

        if( tr.isFirstEvent() ) 
        {
            std::string myTreeName_0l = "myMiniTree_0l";
            std::string myTreeName_1l = "myMiniTree_1l";
            std::string myTreeName_2l = "myMiniTree_2l";

            // ------------
            // for 0 lepton
            // ------------
            my_histos["EventCounterTrain_0l"]->Fill( eventCounter );
            myTreeTrain_0l      = new TTree( (myTreeName_0l).c_str() , (myTreeName_0l).c_str() );
            myMiniTupleTrain_0l = new MiniTupleMaker( myTreeTrain_0l );
            myMiniTupleTrain_0l->setTupleVars(varGeneral);
            myMiniTupleTrain_0l->setTupleVars(varEventShape); 
            myMiniTupleTrain_0l->setTupleVars(varJets);
            myMiniTupleTrain_0l->setTupleVars(varJetsAK8);
            myMiniTupleTrain_0l->setTupleVars(varOldSeed);
            myMiniTupleTrain_0l->setTupleVars(varTopSeed);
            myMiniTupleTrain_0l->initBranches(tr);

            my_histos["EventCounterTest_0l"]->Fill( eventCounter );
            myTreeTest_0l      = new TTree( (myTreeName_0l).c_str() , (myTreeName_0l).c_str() );
            myMiniTupleTest_0l = new MiniTupleMaker( myTreeTest_0l );
            myMiniTupleTest_0l->setTupleVars(varGeneral);
            myMiniTupleTest_0l->setTupleVars(varEventShape);      
            myMiniTupleTest_0l->setTupleVars(varJets);
            myMiniTupleTest_0l->setTupleVars(varJetsAK8);
            myMiniTupleTest_0l->setTupleVars(varOldSeed);
            myMiniTupleTest_0l->setTupleVars(varTopSeed);
            myMiniTupleTest_0l->initBranches(tr);

            my_histos["EventCounterVal_0l"]->Fill( eventCounter );
            myTreeVal_0l      = new TTree( (myTreeName_0l).c_str() , (myTreeName_0l).c_str() );
            myMiniTupleVal_0l = new MiniTupleMaker( myTreeVal_0l ); 
            myMiniTupleVal_0l->setTupleVars(varGeneral);
            myMiniTupleVal_0l->setTupleVars(varEventShape); 
            myMiniTupleVal_0l->setTupleVars(varJets);
            myMiniTupleVal_0l->setTupleVars(varJetsAK8);
            myMiniTupleVal_0l->setTupleVars(varOldSeed);
            myMiniTupleVal_0l->setTupleVars(varTopSeed);
            myMiniTupleVal_0l->initBranches(tr);

            // ------------
            // for 1 lepton
            // ------------
            my_histos["EventCounterTrain_1l"]->Fill( eventCounter );
            myTreeTrain_1l      = new TTree( (myTreeName_1l).c_str() , (myTreeName_1l).c_str() );
            myMiniTupleTrain_1l = new MiniTupleMaker( myTreeTrain_1l );
            myMiniTupleTrain_1l->setTupleVars(varGeneral);
            myMiniTupleTrain_1l->setTupleVars(varEventShape);      
            myMiniTupleTrain_1l->setTupleVars(varJets);
            myMiniTupleTrain_1l->setTupleVars(varJetsAK8);
            myMiniTupleTrain_1l->setTupleVars(varLeptonic);
            myMiniTupleTrain_1l->setTupleVars(varOldSeed);
            myMiniTupleTrain_1l->initBranches(tr);

            my_histos["EventCounterTest_1l"]->Fill( eventCounter );
            myTreeTest_1l      = new TTree( (myTreeName_1l).c_str() , (myTreeName_1l).c_str() );
            myMiniTupleTest_1l = new MiniTupleMaker( myTreeTest_1l );
            myMiniTupleTest_1l->setTupleVars(varGeneral);
            myMiniTupleTest_1l->setTupleVars(varEventShape); 
            myMiniTupleTest_1l->setTupleVars(varJets);
            myMiniTupleTest_1l->setTupleVars(varJetsAK8);
            myMiniTupleTest_1l->setTupleVars(varLeptonic);
            myMiniTupleTest_1l->setTupleVars(varOldSeed);
            myMiniTupleTest_1l->initBranches(tr);

            my_histos["EventCounterVal_1l"]->Fill( eventCounter );
            myTreeVal_1l      = new TTree( (myTreeName_1l).c_str() , (myTreeName_1l).c_str() );
            myMiniTupleVal_1l = new MiniTupleMaker( myTreeVal_1l );
            myMiniTupleVal_1l->setTupleVars(varGeneral);
            myMiniTupleVal_1l->setTupleVars(varEventShape);
            myMiniTupleVal_1l->setTupleVars(varJets);
            myMiniTupleVal_1l->setTupleVars(varJetsAK8);
            myMiniTupleVal_1l->setTupleVars(varLeptonic);
            myMiniTupleVal_1l->setTupleVars(varOldSeed);
            myMiniTupleVal_1l->initBranches(tr);
            
            // ------------
            // for 2 lepton
            // ------------
            my_histos["EventCounterTrain_2l"]->Fill( eventCounter );
            myTreeTrain_2l      = new TTree( (myTreeName_2l).c_str() , (myTreeName_2l).c_str() );
            myMiniTupleTrain_2l = new MiniTupleMaker( myTreeTrain_2l );
            myMiniTupleTrain_2l->setTupleVars(varGeneral);
            myMiniTupleTrain_2l->setTupleVars(varEventShape);
            myMiniTupleTrain_2l->setTupleVars(varJets);
            myMiniTupleTrain_2l->setTupleVars(varJetsAK8);
            myMiniTupleTrain_2l->setTupleVars(varLeptonic);
            myMiniTupleTrain_2l->initBranches(tr);
            
            my_histos["EventCounterTest_2l"]->Fill( eventCounter );
            myTreeTest_2l      = new TTree( (myTreeName_2l).c_str() , (myTreeName_2l).c_str() );
            myMiniTupleTest_2l = new MiniTupleMaker( myTreeTest_2l );
            myMiniTupleTest_2l->setTupleVars(varGeneral);
            myMiniTupleTest_2l->setTupleVars(varEventShape);
            myMiniTupleTest_2l->setTupleVars(varJets);
            myMiniTupleTest_2l->setTupleVars(varJetsAK8);
            myMiniTupleTest_2l->setTupleVars(varLeptonic);
            myMiniTupleTest_2l->initBranches(tr);
            
            my_histos["EventCounterVal_2l"]->Fill( eventCounter );
            myTreeVal_2l      = new TTree( (myTreeName_2l).c_str() , (myTreeName_2l).c_str() );
            myMiniTupleVal_2l = new MiniTupleMaker( myTreeVal_2l );
            myMiniTupleVal_2l->setTupleVars(varGeneral);
            myMiniTupleVal_2l->setTupleVars(varEventShape);
            myMiniTupleVal_2l->setTupleVars(varJets);
            myMiniTupleVal_2l->setTupleVars(varJetsAK8);
            myMiniTupleVal_2l->setTupleVars(varLeptonic);
            myMiniTupleVal_2l->initBranches(tr);

        }
        
        //-----------------------------------
        //-- Fill Histograms Below
        //-----------------------------------
        // for 0 lepton 
        if ( passBaseline0l )
        {
            int mod = count_0l % 10;
            if(mod < 8)
            {
                myMiniTupleTrain_0l->fill();
                numPassTrain_0l++;
            }
            else if(mod == 8)
            {
                myMiniTupleTest_0l->fill();
                numPassTest_0l++;
            }
            else
            {
                myMiniTupleVal_0l->fill();
                numPassVal_0l++;
            }
            count_0l++;
        }

        // for 1 lepton
        if( passBaseline1l ) 
        {
            int mod = count_1l % 10;
            if(mod < 8)
            {
                myMiniTupleTrain_1l->fill();
                numPassTrain_1l++;
            }
            else if(mod == 8)
            {
                myMiniTupleTest_1l->fill();
                numPassTest_1l++;
            }
            else
            {
                myMiniTupleVal_1l->fill();
                numPassVal_1l++;
            }
            count_1l++;
        }

        // for 2 lepton 
        if( passBaseline2l_pt20 )
        {
            int mod = count_2l % 10;
            if(mod < 8)
            {
                myMiniTupleTrain_2l->fill();
                numPassTrain_2l++;
            }
            else if(mod == 8)
            {
                myMiniTupleTest_2l->fill();
                numPassTest_2l++;
            }
            else
            {
                myMiniTupleVal_2l->fill();
                numPassVal_2l++;
            }
            count_2l++;
        }

    }//END of while tr.getNextEvent loop   
    std::cout << "Total_0l: " << count_0l << "   Train_0l: " << numPassTrain_0l << "   Test_0l: " << numPassTest_0l << "   Val_0l: " << numPassVal_0l << std::endl;
    std::cout << "Total_1l: " << count_1l << "   Train_1l: " << numPassTrain_1l << "   Test_1l: " << numPassTest_1l << "   Val_1l: " << numPassVal_1l << std::endl;
    std::cout << "Total_2l: " << count_2l << "   Train_2l: " << numPassTrain_2l << "   Test_2l: " << numPassTest_2l << "   Val_2l: " << numPassVal_2l << std::endl;

}//END of function
      
void MakeNNVariables::WriteHistos( TFile* outfile ) 
{
    const auto& outFileName = std::string(outfile->GetName());
    const auto& name = utility::split("first", outFileName, ".");

    TFile* outfileTrain = TFile::Open((name+"_Train.root").c_str(), "RECREATE");
    outfileTrain->cd();
    myTreeTrain_0l->Write();
    myTreeTrain_1l->Write();
    myTreeTrain_2l->Write();
    my_histos["EventCounterTrain_0l"]->Write();
    my_histos["EventCounterTrain_1l"]->Write();
    my_histos["EventCounterTrain_2l"]->Write();
    delete myTreeTrain_0l; 
    delete myTreeTrain_1l;
    delete myTreeTrain_2l;   
    delete myMiniTupleTrain_0l;
    delete myMiniTupleTrain_1l;
    delete myMiniTupleTrain_2l;
    outfileTrain->Close();

    TFile* outfileTest = TFile::Open((name+"_Test.root").c_str(), "RECREATE");
    outfileTest->cd();
    myTreeTest_0l->Write();
    myTreeTest_1l->Write();
    myTreeTest_2l->Write();
    my_histos["EventCounterTest_0l"]->Write();
    my_histos["EventCounterTest_1l"]->Write();
    my_histos["EventCounterTest_2l"]->Write();
    delete myTreeTest_0l;
    delete myTreeTest_1l;
    delete myTreeTest_2l;    
    delete myMiniTupleTest_0l;
    delete myMiniTupleTest_1l;
    delete myMiniTupleTest_2l;
    outfileTest->Close();

    TFile* outfileVal = TFile::Open((name+"_Val.root").c_str(), "RECREATE");
    outfileVal->cd();
    myTreeVal_0l->Write();
    myTreeVal_1l->Write();
    myTreeVal_2l->Write();
    my_histos["EventCounterVal_0l"]->Write();
    my_histos["EventCounterVal_1l"]->Write();
    my_histos["EventCounterVal_2l"]->Write();
    delete myTreeVal_0l;
    delete myTreeVal_1l; 
    delete myTreeVal_2l;   
    delete myMiniTupleVal_0l;
    delete myMiniTupleVal_1l;
    delete myMiniTupleVal_2l;
    outfileVal->Close();

    remove(outFileName.c_str());
}

