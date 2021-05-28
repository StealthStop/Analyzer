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
    my_histos.emplace( "EventCounterTest_0l",  std::make_shared<TH1D>( "EventCounterTest_0l",  "EventCounterTest_0l",  2, -1.1, 1.1 ) );
    my_histos.emplace( "EventCounterVal_0l",   std::make_shared<TH1D>( "EventCounterVal_0l",   "EventCounterVal_0l",   2, -1.1, 1.1 ) );
    
    my_histos.emplace( "EventCounterTrain_1l", std::make_shared<TH1D>( "EventCounterTrain_1l", "EventCounterTrain_1l", 2, -1.1, 1.1 ) );
    my_histos.emplace( "EventCounterTest_1l",  std::make_shared<TH1D>( "EventCounterTest_1l",  "EventCounterTest_1l",  2, -1.1, 1.1 ) );
    my_histos.emplace( "EventCounterVal_1l",   std::make_shared<TH1D>( "EventCounterVal_1l",   "EventCounterVal_1l",   2, -1.1, 1.1 ) ); 
    
    //my_histos.emplace( "EventCounterTrain_2l", std::make_shared<TH1D>( "EventCounterTrain_2l", "EventCounterTrain_2l", 2, -1.1, 1.1 ) ); 
    //my_histos.emplace( "EventCounterTest_2l",  std::make_shared<TH1D>( "EventCounterTest_2l",  "EventCounterTest_2l",  2, -1.1, 1.1 ) );
    //my_histos.emplace( "EventCounterVal_2l",   std::make_shared<TH1D>( "EventCounterVal_2l",   "EventCounterVal_2l",   2, -1.1, 1.1 ) ); 

}//END of init histos

void MakeNNVariables::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    int count_0l = 0, numPassTrain_0l = 0, numPassTest_0l = 0, numPassVal_0l = 0;
    int count_1l = 0, numPassTrain_1l = 0, numPassTest_1l = 0, numPassVal_1l = 0;
    //int count_2l = 0, numPassTrain_2l = 0, numPassTest_2l = 0, numPassVal_2l = 0;


    while( tr.getNextEvent() )
    {
        const auto& isSignal            = tr.getVar<bool>("isSignal");
        const auto& filetag             = tr.getVar<std::string>("filetag");
        const auto& eventCounter        = tr.getVar<int>("eventCounter");
        const auto& passBaseline0l_Good = tr.getVar<bool>("passBaseline0l_Good"); 
        const auto& passBaseline1l      = tr.getVar<bool>("passBaseline1l_Good");
        //const auto& passBaseline2l_pt20 = tr.getVar<bool>("passBaseline2l_pt20");

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
            "stop1_ptrank_mass",
            "stop2_ptrank_mass",
            "stop1_mrank_mass",
            "stop2_mrank_mass",
            "stop_avemass",
        };

        std::set<std::string> varLeptonic =
        {
            "Mbl",
            "lvMET_cm_m",
            "lvMET_cm_eta",
            "lvMET_cm_phi",
            "lvMET_cm_pt",
            "GoodLeptons_m_1",       "GoodLeptons_m_2",
            "GoodLeptons_eta_1",     "GoodLeptons_eta_2",
            "GoodLeptons_phi_1",     "GoodLeptons_phi_2",
            "GoodLeptons_pt_1",      "GoodLeptons_pt_2",
        };

        std::set<std::string> varOldSeed =
        {
            "MT2_cm_OldSeed",
            "dR_Stop1Stop2_cm_OldSeed",
            "dPhi_Stop1Stop2_cm_OldSeed",
            "Stop1_mass_cm_OldSeed",     "Stop2_mass_cm_OldSeed",
            "Stop1_pt_cm_OldSeed",       "Stop2_pt_cm_OldSeed",
            "Stop1_phi_cm_OldSeed",      "Stop2_phi_cm_OldSeed",
            "Stop1_eta_cm_OldSeed",      "Stop2_eta_cm_OldSeed",
            "Stop1_scalarPt_cm_OldSeed", "Stop2_scalarPt_cm_OldSeed",
        };

        std::set<std::string> varTopSeed =
        {
            "MT2_cm_TopSeed",
            "dR_Stop1Stop2_cm_TopSeed",
            "dPhi_Stop1Stop2_cm_TopSeed",
            "Stop1_mass_cm_TopSeed",     "Stop2_mass_cm_TopSeed",
            "Stop1_pt_cm_TopSeed",       "Stop2_pt_cm_TopSeed",
            "Stop1_phi_cm_TopSeed",      "Stop2_phi_cm_TopSeed",
            "Stop1_eta_cm_TopSeed",      "Stop2_eta_cm_TopSeed",
            "Stop1_scalarPt_cm_TopSeed", "Stop2_scalarPt_cm_TopSeed",
        };

        // -----------------------------------------------
        // get the jet variables separately for 0l, 1l, 2l
        // -----------------------------------------------
        std::vector<std::string> labels = {"_0l", "_1l"}; // for now no 2l variables

        for (std::string label : labels)
        {
            std::set<std::string> varEventShape = 
            {
                "fwm2_top6"+label,    "fwm3_top6"+label,    "fwm4_top6"+label,   "fwm5_top6"+label,
                "fwm6_top6"+label,    "fwm7_top6"+label,    "fwm8_top6"+label,   "fwm9_top6"+label, "fwm10_top6"+label,
                "jmt_ev0_top6"+label, "jmt_ev1_top6"+label, "jmt_ev2_top6"+label,

            };

            std::set<std::string> varJets = 
            {
                "Jet_m_1"+label,            "Jet_m_2"+label,            "Jet_m_3"+label,            "Jet_m_4"+label,             "Jet_m_5"+label,             "Jet_m_6"+label,
                "Jet_m_7"+label,            "Jet_m_8"+label,            "Jet_m_9"+label,            "Jet_m_10"+label,            "Jet_m_11"+label,            "Jet_m_12"+label,
                "Jet_eta_1"+label,          "Jet_eta_2"+label,          "Jet_eta_3"+label,          "Jet_eta_4"+label,           "Jet_eta_5"+label,           "Jet_eta_6"+label,  
                "Jet_eta_7"+label,          "Jet_eta_8"+label,          "Jet_eta_9"+label,          "Jet_eta_10"+label,          "Jet_eta_11"+label,          "Jet_eta_12"+label,
                "Jet_phi_1"+label,          "Jet_phi_2"+label,          "Jet_phi_3"+label,          "Jet_phi_4"+label,           "Jet_phi_5"+label,           "Jet_phi_6"+label,  
                "Jet_phi_7"+label,          "Jet_phi_8"+label,          "Jet_phi_9"+label,          "Jet_phi_10"+label,          "Jet_phi_11"+label,          "Jet_phi_12"+label,
                "Jet_pt_1"+label,           "Jet_pt_2"+label,           "Jet_pt_3"+label,           "Jet_pt_4"+label,            "Jet_pt_5"+label,            "Jet_pt_6"+label, 
                "Jet_pt_7"+label,           "Jet_pt_8"+label,           "Jet_pt_9"+label,           "Jet_pt_10"+label,           "Jet_pt_11"+label,           "Jet_pt_12"+label,
                "Jet_flavb_1"+label,        "Jet_flavb_2"+label,        "Jet_flavb_3"+label,        "Jet_flavb_4"+label,         "Jet_flavb_5"+label,         "Jet_flavb_6"+label, 
                "Jet_flavb_7"+label,        "Jet_flavb_8"+label,        "Jet_flavb_9"+label,        "Jet_flavb_10"+label,        "Jet_flavb_11"+label,        "Jet_flavb_12"+label,
                "Jet_flavg_1"+label,        "Jet_flavg_2"+label,        "Jet_flavg_3"+label,        "Jet_flavg_4"+label,         "Jet_flavg_5"+label,         "Jet_flavg_6"+label, 
                "Jet_flavg_7"+label,        "Jet_flavg_8"+label,        "Jet_flavg_9"+label,        "Jet_flavg_10"+label,        "Jet_flavg_11"+label,        "Jet_flavg_12"+label,
                "Jet_flavc_1"+label,        "Jet_flavc_2"+label,        "Jet_flavc_3"+label,        "Jet_flavc_4"+label,         "Jet_flavc_5"+label,         "Jet_flavc_6"+label, 
                "Jet_flavc_7"+label,        "Jet_flavc_8"+label,        "Jet_flavc_9"+label,        "Jet_flavc_10"+label,        "Jet_flavc_11"+label,        "Jet_flavc_12"+label,
                "Jet_flavuds_1"+label,      "Jet_flavuds_2"+label,      "Jet_flavuds_3"+label,      "Jet_flavuds_4"+label,       "Jet_flavuds_5"+label,       "Jet_flavuds_6"+label, 
                "Jet_flavuds_7"+label,      "Jet_flavuds_8"+label,      "Jet_flavuds_9"+label,      "Jet_flavuds_10"+label,      "Jet_flavuds_11"+label,      "Jet_flavuds_12"+label,
                "Jet_flavq_1"+label,        "Jet_flavq_2"+label,        "Jet_flavq_3"+label,        "Jet_flavq_4"+label,         "Jet_flavq_5"+label,         "Jet_flavq_6"+label, 
                "Jet_flavq_7"+label,        "Jet_flavq_8"+label,        "Jet_flavq_9"+label,        "Jet_flavq_10"+label,        "Jet_flavq_11"+label,        "Jet_flavq_12"+label,
                "Jet_ptD_1"+label,          "Jet_ptD_2"+label,          "Jet_ptD_3"+label,          "Jet_ptD_4"+label,           "Jet_ptD_5"+label,           "Jet_ptD_6"+label,  
                "Jet_ptD_7"+label,          "Jet_ptD_8"+label,          "Jet_ptD_9"+label,          "Jet_ptD_10"+label,          "Jet_ptD_11"+label,          "Jet_ptD_12"+label,
                "Jet_nEF_1"+label,          "Jet_nEF_2"+label,          "Jet_nEF_3"+label,          "Jet_nEF_4"+label,           "Jet_nEF_5"+label,           "Jet_nEF_6"+label,  
                "Jet_nEF_7"+label,          "Jet_nEF_8"+label,          "Jet_nEF_9"+label,          "Jet_nEF_10"+label,          "Jet_nEF_11"+label,          "Jet_nEF_12"+label,
                "Jet_cEF_1"+label,          "Jet_cEF_2"+label,          "Jet_cEF_3"+label,          "Jet_cEF_4"+label,           "Jet_cEF_5"+label,           "Jet_cEF_6"+label,  
                "Jet_cEF_7"+label,          "Jet_cEF_8"+label,          "Jet_cEF_9"+label,          "Jet_cEF_10"+label,          "Jet_cEF_11"+label,          "Jet_cEF_12"+label,
                "Jet_nHF_1"+label,          "Jet_nHF_2"+label,          "Jet_nHF_3"+label,          "Jet_nHF_4"+label,           "Jet_nHF_5"+label,           "Jet_nHF_6"+label,  
                "Jet_nHF_7"+label,          "Jet_nHF_8"+label,          "Jet_nHF_9"+label,          "Jet_nHF_10"+label,          "Jet_nHF_11"+label,          "Jet_nHF_12"+label,
                "Jet_cHF_1"+label,          "Jet_cHF_2"+label,          "Jet_cHF_3"+label,          "Jet_cHF_4"+label,           "Jet_cHF_5"+label,           "Jet_cHF_6"+label,  
                "Jet_cHF_7"+label,          "Jet_cHF_8"+label,          "Jet_cHF_9"+label,          "Jet_cHF_10"+label,          "Jet_cHF_11"+label,          "Jet_cHF_12"+label,
                "Jet_axismajor_1"+label,    "Jet_axismajor_2"+label,    "Jet_axismajor_3"+label,    "Jet_axismajor_4"+label,     "Jet_axismajor_5"+label,     "Jet_axismajor_6"+label, 
                "Jet_axismajor_7"+label,    "Jet_axismajor_8"+label,    "Jet_axismajor_9"+label,    "Jet_axismajor_10"+label,    "Jet_axismajor_11"+label,    "Jet_axismajor_12"+label,
                "Jet_axisminor_1"+label,    "Jet_axisminor_2"+label,    "Jet_axisminor_3"+label,    "Jet_axisminor_4"+label,     "Jet_axisminor_5"+label,     "Jet_axisminor_6"+label, 
                "Jet_axisminor_7"+label,    "Jet_axisminor_8"+label,    "Jet_axisminor_9"+label,    "Jet_axisminor_10"+label,    "Jet_axisminor_11"+label,    "Jet_axisminor_12"+label,
                "Jet_multiplicity_1"+label, "Jet_multiplicity_2"+label, "Jet_multiplicity_3"+label, "Jet_multiplicity_4"+label,  "Jet_multiplicity_5"+label,  "Jet_multiplicity_6"+label, 
                "Jet_multiplicity_7"+label, "Jet_multiplicity_8"+label, "Jet_multiplicity_9"+label, "Jet_multiplicity_10"+label, "Jet_multiplicity_11"+label, "Jet_multiplicity_12"+label,
            };            

            std::set<std::string> varJetsAK8 =
            {
                "JetsAK8_m_1"+label,              "JetsAK8_m_2"+label,              "JetsAK8_m_3"+label,              "JetsAK8_m_4"+label,              "JetsAK8_m_5"+label,
                "JetsAK8_eta_1"+label,            "JetsAK8_eta_2"+label,            "JetsAK8_eta_3"+label,            "JetsAK8_eta_4"+label,            "JetsAK8_eta_5"+label,
                "JetsAK8_phi_1"+label,            "JetsAK8_phi_2"+label,            "JetsAK8_phi_3"+label,            "JetsAK8_phi_4"+label,            "JetsAK8_phi_5"+label,
                "JetsAK8_pt_1"+label,             "JetsAK8_pt_2"+label,             "JetsAK8_pt_3"+label,             "JetsAK8_pt_4"+label,             "JetsAK8_pt_5"+label,            
                "JetsAK8_SDM_1"+label,            "JetsAK8_SDM_2"+label,            "JetsAK8_SDM_3"+label,            "JetsAK8_SDM_4"+label,            "JetsAK8_SDM_5"+label,
                "JetsAK8_Pruned_1"+label,         "JetsAK8_Pruned_2"+label,         "JetsAK8_Pruned_3"+label,         "JetsAK8_Pruned_4"+label,         "JetsAK8_Pruned_5"+label,
                "JetsAK8_Tau1_1"+label,           "JetsAK8_Tau1_2"+label,           "JetsAK8_Tau1_3"+label,           "JetsAK8_Tau1_4"+label,           "JetsAK8_Tau1_5"+label,
                "JetsAK8_Tau2_1"+label,           "JetsAK8_Tau2_2"+label,           "JetsAK8_Tau2_3"+label,           "JetsAK8_Tau2_4"+label,           "JetsAK8_Tau2_5"+label,
                "JetsAK8_Tau3_1"+label,           "JetsAK8_Tau3_2"+label,           "JetsAK8_Tau3_3"+label,           "JetsAK8_Tau3_4"+label,           "JetsAK8_Tau3_5"+label,
                "JetsAK8_axismajor_1"+label,      "JetsAK8_axismajor_2"+label,      "JetsAK8_axismajor_3"+label,      "JetsAK8_axismajor_4"+label,      "JetsAK8_axismajor_5"+label,
                "JetsAK8_axisminor_1"+label,      "JetsAK8_axisminor_2"+label,      "JetsAK8_axisminor_3"+label,      "JetsAK8_axisminor_4"+label,      "JetsAK8_axisminor_5"+label,
                "JetsAK8_nsubjets_1"+label,       "JetsAK8_nsubjets_2"+label,       "JetsAK8_nsubjets_3"+label,       "JetsAK8_nsubjets_4"+label,       "JetsAK8_nsubjets_5"+label,
                "JetsAK8_tDiscriminator_1"+label, "JetsAK8_tDiscriminator_2"+label, "JetsAK8_tDiscriminator_3"+label, "JetsAK8_tDiscriminator_4"+label, "JetsAK8_tDiscriminator_5"+label,
                "JetsAK8_wDiscriminator_1"+label, "JetsAK8_wDiscriminator_2"+label, "JetsAK8_wDiscriminator_3"+label, "JetsAK8_wDiscriminator_4"+label, "JetsAK8_wDiscriminator_5"+label,
                "JetsAK8_hDiscriminator_1"+label, "JetsAK8_hDiscriminator_2"+label, "JetsAK8_hDiscriminator_3"+label, "JetsAK8_hDiscriminator_4"+label, "JetsAK8_hDiscriminator_5"+label,
                "JetsAK8_multiplicity_1"+label,   "JetsAK8_multiplicity_2"+label,   "JetsAK8_multiplicity_3"+label,   "JetsAK8_multiplicity_4"+label,   "JetsAK8_multiplicity_5"+label,
            };

            if( tr.isFirstEvent() ) 
            {
                std::string myTreeName_0l = "myMiniTree_0l";
                std::string myTreeName_1l = "myMiniTree_1l";
                //std::string myTreeName_2l = "myMiniTree_2l";

                if (label == "_0l")
                {
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
                }

                if (label == "_1l")
                {
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
                }           

                //else
                //{
                //  my_histos["EventCounterTrain_2l"]->Fill( eventCounter );
                //  myTreeTrain_2l      = new TTree( (myTreeName_2l).c_str() , (myTreeName_2l).c_str() );
                //  myMiniTupleTrain_2l = new MiniTupleMaker( myTreeTrain_2l );
                //  myMiniTupleTrain_2l->setTupleVars(varGeneral);
                //  myMiniTupleTrain_2l->setTupleVars(varEventShape);
                //  myMiniTupleTrain_2l->setTupleVars(varJets);
                //  myMiniTupleTrain_2l->setTupleVars(varJetsAK8);
                //  myMiniTupleTrain_2l->setTupleVars(varLeptonic);
                //  myMiniTupleTrain_2l->initBranches(tr);
                //  
                //  my_histos["EventCounterTest_2l"]->Fill( eventCounter );
                //  myTreeTest_2l      = new TTree( (myTreeName_2l).c_str() , (myTreeName_2l).c_str() );
                //  myMiniTupleTest_2l = new MiniTupleMaker( myTreeTest_2l );
                //  myMiniTupleTest_2l->setTupleVars(varGeneral);
                //  myMiniTupleTest_2l->setTupleVars(varEventShape);
                //  myMiniTupleTest_2l->setTupleVars(varJets);
                //  myMiniTupleTest_2l->setTupleVars(varJetsAK8);
                //  myMiniTupleTest_2l->setTupleVars(varLeptonic);
                //  myMiniTupleTest_2l->initBranches(tr);
                //  
                //  my_histos["EventCounterVal_2l"]->Fill( eventCounter );
                //  myTreeVal_2l      = new TTree( (myTreeName_2l).c_str() , (myTreeName_2l).c_str() );
                //  myMiniTupleVal_2l = new MiniTupleMaker( myTreeVal_2l );
                //  myMiniTupleVal_2l->setTupleVars(varGeneral);
                //  myMiniTupleVal_2l->setTupleVars(varEventShape);
                //  myMiniTupleVal_2l->setTupleVars(varJets);
                //  myMiniTupleVal_2l->setTupleVars(varJetsAK8);
                //  myMiniTupleVal_2l->setTupleVars(varLeptonic);
                //  myMiniTupleVal_2l->initBranches(tr);
                //} 

            } //  
        } // jet varaible loop
        
        //-----------------------------------
        //-- Fill Histograms Below
        //-----------------------------------
        // for 0 lepton 
        if ( passBaseline0l_Good )
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
        //if( passBaseline2l_pt20 )
        //{
        //    int mod = count_2l % 10;
        //    if(mod < 8)
        //    {
        //        myMiniTupleTrain_2l->fill();
        //        numPassTrain_2l++;
        //    }
        //    else if(mod == 8)
        //    {
        //        myMiniTupleTest_2l->fill();
        //        numPassTest_2l++;
        //    }
        //    else
        //    {
        //        myMiniTupleVal_2l->fill();
        //        numPassVal_2l++;
        //    }
        //    count_2l++;
        //}

    }//END of while tr.getNextEvent loop   
    std::cout << "Total_0l: " << count_0l << "   Train_0l: " << numPassTrain_0l << "   Test_0l: " << numPassTest_0l << "   Val_0l: " << numPassVal_0l << std::endl;
    std::cout << "Total_1l: " << count_1l << "   Train_1l: " << numPassTrain_1l << "   Test_1l: " << numPassTest_1l << "   Val_1l: " << numPassVal_1l << std::endl;
    //std::cout << "Total_2l: " << count_2l << "   Train_2l: " << numPassTrain_2l << "   Test_2l: " << numPassTest_2l << "   Val_2l: " << numPassVal_2l << std::endl;

}//END of function
      
void MakeNNVariables::WriteHistos( TFile* outfile ) 
{
    const auto& outFileName = std::string(outfile->GetName());
    const auto& name = utility::split("first", outFileName, ".");

    TFile* outfileTrain = TFile::Open((name+"_Train.root").c_str(), "RECREATE");
    outfileTrain->cd();
    myTreeTrain_0l->Write();
    my_histos["EventCounterTrain_0l"]->Write();
    myTreeTrain_1l->Write();
    my_histos["EventCounterTrain_1l"]->Write();
    //myTreeTrain_2l->Write();
    //my_histos["EventCounterTrain_2l"]->Write();
    delete myTreeTrain_0l; 
    delete myMiniTupleTrain_0l;
    delete myTreeTrain_1l;
    delete myMiniTupleTrain_1l;
    //delete myTreeTrain_2l;   
    //delete myMiniTupleTrain_2l;
    outfileTrain->Close();

    TFile* outfileTest = TFile::Open((name+"_Test.root").c_str(), "RECREATE");
    outfileTest->cd();
    myTreeTest_0l->Write();
    my_histos["EventCounterTest_0l"]->Write();
    myTreeTest_1l->Write();
    my_histos["EventCounterTest_1l"]->Write();
    //myTreeTest_2l->Write();
    //my_histos["EventCounterTest_2l"]->Write();
    delete myTreeTest_0l;
    delete myMiniTupleTest_0l;
    delete myTreeTest_1l;
    delete myMiniTupleTest_1l;
    //delete myTreeTest_2l;    
    //delete myMiniTupleTest_2l;
    outfileTest->Close();

    TFile* outfileVal = TFile::Open((name+"_Val.root").c_str(), "RECREATE");
    outfileVal->cd();
    myTreeVal_0l->Write();
    my_histos["EventCounterVal_0l"]->Write();
    myTreeVal_1l->Write();
    my_histos["EventCounterVal_1l"]->Write();
    //myTreeVal_2l->Write();
    //my_histos["EventCounterVal_2l"]->Write();
    delete myTreeVal_0l;
    delete myMiniTupleVal_0l;
    delete myTreeVal_1l; 
    delete myMiniTupleVal_1l;
    //delete myTreeVal_2l;   
    //delete myMiniTupleVal_2l;
    outfileVal->Close();

    remove(outFileName.c_str());
}

