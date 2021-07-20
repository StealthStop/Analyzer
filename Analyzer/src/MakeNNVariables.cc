#define MakeNNVariables_cxx
#include "Analyzer/Analyzer/include/MakeNNVariables.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/MiniTupleMaker.h"
#include "Framework/Framework/include/Utility.h" 

#include "Framework/Framework/include/Muon.h"
#include "Framework/Framework/include/Electron.h"
#include "Framework/Framework/include/Photon.h"
#include "Framework/Framework/include/Jet.h"
#include "Framework/Framework/include/BJet.h"
#include "Framework/Framework/include/CommonVariables.h"
#include "Framework/Framework/include/RunTopTagger.h"
#include "Framework/Framework/include/FatJetCombine.h"
#include "Framework/Framework/include/MakeMVAVariables.h"
#include "Framework/Framework/include/Baseline.h"
#include "Framework/Framework/include/BTagCorrector.h"
#include "Framework/Framework/include/ScaleFactors.h"
#include "Framework/Framework/include/StopJets.h"
#include "Framework/Framework/include/StopGenMatch.h"
#include "Framework/Framework/include/MakeStopHemispheres.h"

#include <iostream>
#include <stdio.h> 

MakeNNVariables::MakeNNVariables()
{
    InitHistos();
    my_var_suffix = {"", "JECup", "JECdown", "JERup", "JERdown"};
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
    std::map<std::string, std::map<std::string, std::map<std::string, int> > > counts;

    const auto& filetag = tr.getVar<std::string>("filetag");
    const auto& runYear = tr.getVar<std::string>("runYear");
    const auto& bjetFileName = tr.getVar<std::string>("bjetFileName");
    const auto& bjetCSVFileName = tr.getVar<std::string>("bjetCSVFileName");
    const auto& leptonFileName = tr.getVar<std::string>("leptonFileName");
    const auto& puFileName = tr.getVar<std::string>("puFileName");
    const auto& meanFileName = tr.getVar<std::string>("meanFileName");
    const auto& TopTaggerCfg = tr.getVar<std::string>("TopTaggerCfg");

    for(const auto& myVarSuffix : my_var_suffix)
    {
        if(myVarSuffix == "") continue;
        Muon muon(myVarSuffix);
        Electron electron(myVarSuffix);
        Photon photon(myVarSuffix);
        Jet jet(myVarSuffix);
        BJet bjet(myVarSuffix);
        RunTopTagger topTagger(TopTaggerCfg, myVarSuffix);
        CommonVariables commonVariables(myVarSuffix);
        FatJetCombine fatJetCombine(myVarSuffix);
        MakeMVAVariables makeMVAVariables0L(false, myVarSuffix, "GoodJets_pt45", false, true, 12, 2, "_0l");
        MakeMVAVariables makeMVAVariables1L(false, myVarSuffix, "GoodJets_pt30", false, true, 12, 2, "_1l");
        Baseline baseline(myVarSuffix);
        BTagCorrector bTagCorrector(bjetFileName, "", bjetCSVFileName, filetag);
        bTagCorrector.SetVarNames("GenParticles_PdgId", "Jets"+myVarSuffix, "GoodJets_pt30"+myVarSuffix, "Jets"+myVarSuffix+"_bJetTagDeepCSVtotb", "Jets"+myVarSuffix+"_partonFlavor", myVarSuffix);
        ScaleFactors scaleFactors( runYear, leptonFileName, puFileName, meanFileName, myVarSuffix);
        StopJets stopJets(myVarSuffix);
        StopGenMatch stopGenMatch(myVarSuffix);
        MakeStopHemispheres stopHemispheres("Jets"+myVarSuffix, "GoodJets_pt20"+myVarSuffix, "NGoodJets_pt20"+myVarSuffix, "_OldSeed"+myVarSuffix, Hemisphere::InvMassSeed);
  
        tr.registerFunction(muon);
        tr.registerFunction(electron);
        tr.registerFunction(photon);
        tr.registerFunction(jet);
        tr.registerFunction(bjet);
        tr.registerFunction(topTagger);
        tr.registerFunction(commonVariables);
        tr.registerFunction(baseline);
        tr.registerFunction(fatJetCombine);
        tr.registerFunction(makeMVAVariables0L);
        tr.registerFunction(makeMVAVariables1L);
        tr.registerFunction(stopJets);
        tr.registerFunction(stopHemispheres);
        tr.registerFunction(bTagCorrector);
        tr.registerFunction(scaleFactors);
        tr.registerFunction(stopGenMatch);
    }

    while( tr.getNextEvent() )
    {
        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() % 1000 == 0 ) printf( " Event %i\n", tr.getEvtNum() );        


        for(const auto& myVarSuffix : my_var_suffix)
        {
            const auto& isSignal        = tr.getVar<bool>("isSignal");
            const auto& filetag         = tr.getVar<std::string>("filetag");
            const auto& eventCounter    = tr.getVar<int>("eventCounter");

            std::map<std::string, bool> baselines;
            baselines["_0l"] = tr.getVar<bool>("passBaseline0l_Good"+myVarSuffix); 
            baselines["_1l"] = tr.getVar<bool>("passBaseline1l_Good"+myVarSuffix);
            //baselines["_2l"] = tr.getVar<bool>("passBaseline2l_pt20");

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

            auto& model = tr.createDerivedVar<int>("model", 0);
            if(filetag.find("RPV") != std::string::npos)
            {
                model = 100;
            } else if (filetag.find("SYY") != std::string::npos) 
            {
                model = 101;
            } else if (filetag.find("SHH") != std::string::npos)
            {
                model = 102;
            }
       
            //------------------------------------
            //-- Print Event Number
            //------------------------------------
           
            std::vector<std::string> labels = {"_0l", "_1l"}; // for now no 2l variables
            std::vector<std::string> splits = {"count", "Train", "Test", "Val"};

            if( tr.isFirstEvent() ) 
            {

                for (const auto& label : labels)
                {
                    std::string myTreeName = "myMiniTree"+label+myVarSuffix;

                    for (const auto& split : splits)
                    {
                        counts[split][label][myVarSuffix] = 0;

                        if (split != "count")
                        {
                            myTree[split][label][myVarSuffix]      = new TTree( (myTreeName).c_str() , (myTreeName).c_str() );
                            myMiniTuple[split][label][myVarSuffix] = new MiniTupleMaker( myTree[split][label][myVarSuffix] );
                        }
                    }
                }

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
                    "HT_trigger_pt30"+myVarSuffix,
                    "HT_trigger_pt45"+myVarSuffix,
                    "NGoodJets_pt20_double"+myVarSuffix,
                    "NGoodJets_pt30_double"+myVarSuffix,
                    "NGoodJets_pt45_double"+myVarSuffix,
                    "NGoodBJets_pt30_double"+myVarSuffix,
                    "NGoodBJets_pt45_double"+myVarSuffix,          
                    "stop1_ptrank_mass"+myVarSuffix,
                    "stop2_ptrank_mass"+myVarSuffix,
                    "stop1_mrank_mass"+myVarSuffix,
                    "stop2_mrank_mass"+myVarSuffix,
                    "stop_avemass"+myVarSuffix,
                };

                std::set<std::string> varLeptonic =
                {
                    "Mbl"+myVarSuffix,
                    "lvMET_cm_m"+myVarSuffix,
                    "lvMET_cm_eta"+myVarSuffix,
                    "lvMET_cm_phi"+myVarSuffix,
                    "lvMET_cm_pt"+myVarSuffix,
                    "GoodLeptons_m_1"+myVarSuffix,       "GoodLeptons_m_2"+myVarSuffix,
                    "GoodLeptons_eta_1"+myVarSuffix,     "GoodLeptons_eta_2"+myVarSuffix,
                    "GoodLeptons_phi_1"+myVarSuffix,     "GoodLeptons_phi_2"+myVarSuffix,
                    "GoodLeptons_pt_1"+myVarSuffix,      "GoodLeptons_pt_2"+myVarSuffix,
                };

                std::set<std::string> varOldSeed =
                {
                    "MT2_cm_OldSeed"+myVarSuffix,
                    "dR_Stop1Stop2_cm_OldSeed"+myVarSuffix,
                    "dPhi_Stop1Stop2_cm_OldSeed"+myVarSuffix,
                    "Stop1_mass_cm_OldSeed"+myVarSuffix,     "Stop2_mass_cm_OldSeed"+myVarSuffix,
                    "Stop1_pt_cm_OldSeed"+myVarSuffix,       "Stop2_pt_cm_OldSeed"+myVarSuffix,
                    "Stop1_phi_cm_OldSeed"+myVarSuffix,      "Stop2_phi_cm_OldSeed"+myVarSuffix,
                    "Stop1_eta_cm_OldSeed"+myVarSuffix,      "Stop2_eta_cm_OldSeed"+myVarSuffix,
                    "Stop1_scalarPt_cm_OldSeed"+myVarSuffix, "Stop2_scalarPt_cm_OldSeed"+myVarSuffix,
                };

                // -----------------------------------------------
                // get the jet variables separately for 0l, 1l, 2l
                // -----------------------------------------------

                for (std::string label : labels)
                {
                    std::set<std::string> varEventShape = 
                    {
                        "fwm2_top6"+label+myVarSuffix,    "fwm3_top6"+label+myVarSuffix,    "fwm4_top6"+label+myVarSuffix,   "fwm5_top6"+label+myVarSuffix,
                        "fwm6_top6"+label+myVarSuffix,    "fwm7_top6"+label+myVarSuffix,    "fwm8_top6"+label+myVarSuffix,   "fwm9_top6"+label+myVarSuffix, "fwm10_top6"+label+myVarSuffix,
                        "jmt_ev0_top6"+label+myVarSuffix, "jmt_ev1_top6"+label+myVarSuffix, "jmt_ev2_top6"+label+myVarSuffix,

                    };

                    std::set<std::string> varJets = 
                    {
                        "Jet_m_1"+label+myVarSuffix,            "Jet_m_2"+label+myVarSuffix,            "Jet_m_3"+label+myVarSuffix,            "Jet_m_4"+label+myVarSuffix,             "Jet_m_5"+label+myVarSuffix,             "Jet_m_6"+label+myVarSuffix,
                        "Jet_m_7"+label+myVarSuffix,            "Jet_m_8"+label+myVarSuffix,            "Jet_m_9"+label+myVarSuffix,            "Jet_m_10"+label+myVarSuffix,            "Jet_m_11"+label+myVarSuffix,            "Jet_m_12"+label+myVarSuffix,
                        "Jet_eta_1"+label+myVarSuffix,          "Jet_eta_2"+label+myVarSuffix,          "Jet_eta_3"+label+myVarSuffix,          "Jet_eta_4"+label+myVarSuffix,           "Jet_eta_5"+label+myVarSuffix,           "Jet_eta_6"+label+myVarSuffix,  
                        "Jet_eta_7"+label+myVarSuffix,          "Jet_eta_8"+label+myVarSuffix,          "Jet_eta_9"+label+myVarSuffix,          "Jet_eta_10"+label+myVarSuffix,          "Jet_eta_11"+label+myVarSuffix,          "Jet_eta_12"+label+myVarSuffix,
                        "Jet_phi_1"+label+myVarSuffix,          "Jet_phi_2"+label+myVarSuffix,          "Jet_phi_3"+label+myVarSuffix,          "Jet_phi_4"+label+myVarSuffix,           "Jet_phi_5"+label+myVarSuffix,           "Jet_phi_6"+label+myVarSuffix,  
                        "Jet_phi_7"+label+myVarSuffix,          "Jet_phi_8"+label+myVarSuffix,          "Jet_phi_9"+label+myVarSuffix,          "Jet_phi_10"+label+myVarSuffix,          "Jet_phi_11"+label+myVarSuffix,          "Jet_phi_12"+label+myVarSuffix,
                        "Jet_pt_1"+label+myVarSuffix,           "Jet_pt_2"+label+myVarSuffix,           "Jet_pt_3"+label+myVarSuffix,           "Jet_pt_4"+label+myVarSuffix,            "Jet_pt_5"+label+myVarSuffix,            "Jet_pt_6"+label+myVarSuffix, 
                        "Jet_pt_7"+label+myVarSuffix,           "Jet_pt_8"+label+myVarSuffix,           "Jet_pt_9"+label+myVarSuffix,           "Jet_pt_10"+label+myVarSuffix,           "Jet_pt_11"+label+myVarSuffix,           "Jet_pt_12"+label+myVarSuffix,
                        "Jet_flavb_1"+label+myVarSuffix,        "Jet_flavb_2"+label+myVarSuffix,        "Jet_flavb_3"+label+myVarSuffix,        "Jet_flavb_4"+label+myVarSuffix,         "Jet_flavb_5"+label+myVarSuffix,         "Jet_flavb_6"+label+myVarSuffix, 
                        "Jet_flavb_7"+label+myVarSuffix,        "Jet_flavb_8"+label+myVarSuffix,        "Jet_flavb_9"+label+myVarSuffix,        "Jet_flavb_10"+label+myVarSuffix,        "Jet_flavb_11"+label+myVarSuffix,        "Jet_flavb_12"+label+myVarSuffix,
                        "Jet_flavg_1"+label+myVarSuffix,        "Jet_flavg_2"+label+myVarSuffix,        "Jet_flavg_3"+label+myVarSuffix,        "Jet_flavg_4"+label+myVarSuffix,         "Jet_flavg_5"+label+myVarSuffix,         "Jet_flavg_6"+label+myVarSuffix, 
                        "Jet_flavg_7"+label+myVarSuffix,        "Jet_flavg_8"+label+myVarSuffix,        "Jet_flavg_9"+label+myVarSuffix,        "Jet_flavg_10"+label+myVarSuffix,        "Jet_flavg_11"+label+myVarSuffix,        "Jet_flavg_12"+label+myVarSuffix,
                        "Jet_flavc_1"+label+myVarSuffix,        "Jet_flavc_2"+label+myVarSuffix,        "Jet_flavc_3"+label+myVarSuffix,        "Jet_flavc_4"+label+myVarSuffix,         "Jet_flavc_5"+label+myVarSuffix,         "Jet_flavc_6"+label+myVarSuffix, 
                        "Jet_flavc_7"+label+myVarSuffix,        "Jet_flavc_8"+label+myVarSuffix,        "Jet_flavc_9"+label+myVarSuffix,        "Jet_flavc_10"+label+myVarSuffix,        "Jet_flavc_11"+label+myVarSuffix,        "Jet_flavc_12"+label+myVarSuffix,
                        "Jet_flavuds_1"+label+myVarSuffix,      "Jet_flavuds_2"+label+myVarSuffix,      "Jet_flavuds_3"+label+myVarSuffix,      "Jet_flavuds_4"+label+myVarSuffix,       "Jet_flavuds_5"+label+myVarSuffix,       "Jet_flavuds_6"+label+myVarSuffix, 
                        "Jet_flavuds_7"+label+myVarSuffix,      "Jet_flavuds_8"+label+myVarSuffix,      "Jet_flavuds_9"+label+myVarSuffix,      "Jet_flavuds_10"+label+myVarSuffix,      "Jet_flavuds_11"+label+myVarSuffix,      "Jet_flavuds_12"+label+myVarSuffix,
                        "Jet_flavq_1"+label+myVarSuffix,        "Jet_flavq_2"+label+myVarSuffix,        "Jet_flavq_3"+label+myVarSuffix,        "Jet_flavq_4"+label+myVarSuffix,         "Jet_flavq_5"+label+myVarSuffix,         "Jet_flavq_6"+label+myVarSuffix, 
                        "Jet_flavq_7"+label+myVarSuffix,        "Jet_flavq_8"+label+myVarSuffix,        "Jet_flavq_9"+label+myVarSuffix,        "Jet_flavq_10"+label+myVarSuffix,        "Jet_flavq_11"+label+myVarSuffix,        "Jet_flavq_12"+label+myVarSuffix,
                        "Jet_ptD_1"+label+myVarSuffix,          "Jet_ptD_2"+label+myVarSuffix,          "Jet_ptD_3"+label+myVarSuffix,          "Jet_ptD_4"+label+myVarSuffix,           "Jet_ptD_5"+label+myVarSuffix,           "Jet_ptD_6"+label+myVarSuffix,  
                        "Jet_ptD_7"+label+myVarSuffix,          "Jet_ptD_8"+label+myVarSuffix,          "Jet_ptD_9"+label+myVarSuffix,          "Jet_ptD_10"+label+myVarSuffix,          "Jet_ptD_11"+label+myVarSuffix,          "Jet_ptD_12"+label+myVarSuffix,
                        "Jet_nEF_1"+label+myVarSuffix,          "Jet_nEF_2"+label+myVarSuffix,          "Jet_nEF_3"+label+myVarSuffix,          "Jet_nEF_4"+label+myVarSuffix,           "Jet_nEF_5"+label+myVarSuffix,           "Jet_nEF_6"+label+myVarSuffix,  
                        "Jet_nEF_7"+label+myVarSuffix,          "Jet_nEF_8"+label+myVarSuffix,          "Jet_nEF_9"+label+myVarSuffix,          "Jet_nEF_10"+label+myVarSuffix,          "Jet_nEF_11"+label+myVarSuffix,          "Jet_nEF_12"+label+myVarSuffix,
                        "Jet_cEF_1"+label+myVarSuffix,          "Jet_cEF_2"+label+myVarSuffix,          "Jet_cEF_3"+label+myVarSuffix,          "Jet_cEF_4"+label+myVarSuffix,           "Jet_cEF_5"+label+myVarSuffix,           "Jet_cEF_6"+label+myVarSuffix,  
                        "Jet_cEF_7"+label+myVarSuffix,          "Jet_cEF_8"+label+myVarSuffix,          "Jet_cEF_9"+label+myVarSuffix,          "Jet_cEF_10"+label+myVarSuffix,          "Jet_cEF_11"+label+myVarSuffix,          "Jet_cEF_12"+label+myVarSuffix,
                        "Jet_nHF_1"+label+myVarSuffix,          "Jet_nHF_2"+label+myVarSuffix,          "Jet_nHF_3"+label+myVarSuffix,          "Jet_nHF_4"+label+myVarSuffix,           "Jet_nHF_5"+label+myVarSuffix,           "Jet_nHF_6"+label+myVarSuffix,  
                        "Jet_nHF_7"+label+myVarSuffix,          "Jet_nHF_8"+label+myVarSuffix,          "Jet_nHF_9"+label+myVarSuffix,          "Jet_nHF_10"+label+myVarSuffix,          "Jet_nHF_11"+label+myVarSuffix,          "Jet_nHF_12"+label+myVarSuffix,
                        "Jet_cHF_1"+label+myVarSuffix,          "Jet_cHF_2"+label+myVarSuffix,          "Jet_cHF_3"+label+myVarSuffix,          "Jet_cHF_4"+label+myVarSuffix,           "Jet_cHF_5"+label+myVarSuffix,           "Jet_cHF_6"+label+myVarSuffix,  
                        "Jet_cHF_7"+label+myVarSuffix,          "Jet_cHF_8"+label+myVarSuffix,          "Jet_cHF_9"+label+myVarSuffix,          "Jet_cHF_10"+label+myVarSuffix,          "Jet_cHF_11"+label+myVarSuffix,          "Jet_cHF_12"+label+myVarSuffix,
                        "Jet_axismajor_1"+label+myVarSuffix,    "Jet_axismajor_2"+label+myVarSuffix,    "Jet_axismajor_3"+label+myVarSuffix,    "Jet_axismajor_4"+label+myVarSuffix,     "Jet_axismajor_5"+label+myVarSuffix,     "Jet_axismajor_6"+label+myVarSuffix, 
                        "Jet_axismajor_7"+label+myVarSuffix,    "Jet_axismajor_8"+label+myVarSuffix,    "Jet_axismajor_9"+label+myVarSuffix,    "Jet_axismajor_10"+label+myVarSuffix,    "Jet_axismajor_11"+label+myVarSuffix,    "Jet_axismajor_12"+label+myVarSuffix,
                        "Jet_axisminor_1"+label+myVarSuffix,    "Jet_axisminor_2"+label+myVarSuffix,    "Jet_axisminor_3"+label+myVarSuffix,    "Jet_axisminor_4"+label+myVarSuffix,     "Jet_axisminor_5"+label+myVarSuffix,     "Jet_axisminor_6"+label+myVarSuffix, 
                        "Jet_axisminor_7"+label+myVarSuffix,    "Jet_axisminor_8"+label+myVarSuffix,    "Jet_axisminor_9"+label+myVarSuffix,    "Jet_axisminor_10"+label+myVarSuffix,    "Jet_axisminor_11"+label+myVarSuffix,    "Jet_axisminor_12"+label+myVarSuffix,
                        "Jet_multiplicity_1"+label+myVarSuffix, "Jet_multiplicity_2"+label+myVarSuffix, "Jet_multiplicity_3"+label+myVarSuffix, "Jet_multiplicity_4"+label+myVarSuffix,  "Jet_multiplicity_5"+label+myVarSuffix,  "Jet_multiplicity_6"+label+myVarSuffix, 
                        "Jet_multiplicity_7"+label+myVarSuffix, "Jet_multiplicity_8"+label+myVarSuffix, "Jet_multiplicity_9"+label+myVarSuffix, "Jet_multiplicity_10"+label+myVarSuffix, "Jet_multiplicity_11"+label+myVarSuffix, "Jet_multiplicity_12"+label+myVarSuffix,
                    };            

                    std::set<std::string> varJetsAK8 =
                    {
                        "JetsAK8_m_1"+label+myVarSuffix,              "JetsAK8_m_2"+label+myVarSuffix,              "JetsAK8_m_3"+label+myVarSuffix,              "JetsAK8_m_4"+label+myVarSuffix,              "JetsAK8_m_5"+label+myVarSuffix,
                        "JetsAK8_eta_1"+label+myVarSuffix,            "JetsAK8_eta_2"+label+myVarSuffix,            "JetsAK8_eta_3"+label+myVarSuffix,            "JetsAK8_eta_4"+label+myVarSuffix,            "JetsAK8_eta_5"+label+myVarSuffix,
                        "JetsAK8_phi_1"+label+myVarSuffix,            "JetsAK8_phi_2"+label+myVarSuffix,            "JetsAK8_phi_3"+label+myVarSuffix,            "JetsAK8_phi_4"+label+myVarSuffix,            "JetsAK8_phi_5"+label+myVarSuffix,
                        "JetsAK8_pt_1"+label+myVarSuffix,             "JetsAK8_pt_2"+label+myVarSuffix,             "JetsAK8_pt_3"+label+myVarSuffix,             "JetsAK8_pt_4"+label+myVarSuffix,             "JetsAK8_pt_5"+label+myVarSuffix,            
                        "JetsAK8_SDM_1"+label+myVarSuffix,            "JetsAK8_SDM_2"+label+myVarSuffix,            "JetsAK8_SDM_3"+label+myVarSuffix,            "JetsAK8_SDM_4"+label+myVarSuffix,            "JetsAK8_SDM_5"+label+myVarSuffix,
                        "JetsAK8_Pruned_1"+label+myVarSuffix,         "JetsAK8_Pruned_2"+label+myVarSuffix,         "JetsAK8_Pruned_3"+label+myVarSuffix,         "JetsAK8_Pruned_4"+label+myVarSuffix,         "JetsAK8_Pruned_5"+label+myVarSuffix,
                        "JetsAK8_Tau1_1"+label+myVarSuffix,           "JetsAK8_Tau1_2"+label+myVarSuffix,           "JetsAK8_Tau1_3"+label+myVarSuffix,           "JetsAK8_Tau1_4"+label+myVarSuffix,           "JetsAK8_Tau1_5"+label+myVarSuffix,
                        "JetsAK8_Tau2_1"+label+myVarSuffix,           "JetsAK8_Tau2_2"+label+myVarSuffix,           "JetsAK8_Tau2_3"+label+myVarSuffix,           "JetsAK8_Tau2_4"+label+myVarSuffix,           "JetsAK8_Tau2_5"+label+myVarSuffix,
                        "JetsAK8_Tau3_1"+label+myVarSuffix,           "JetsAK8_Tau3_2"+label+myVarSuffix,           "JetsAK8_Tau3_3"+label+myVarSuffix,           "JetsAK8_Tau3_4"+label+myVarSuffix,           "JetsAK8_Tau3_5"+label+myVarSuffix,
                        "JetsAK8_axismajor_1"+label+myVarSuffix,      "JetsAK8_axismajor_2"+label+myVarSuffix,      "JetsAK8_axismajor_3"+label+myVarSuffix,      "JetsAK8_axismajor_4"+label+myVarSuffix,      "JetsAK8_axismajor_5"+label+myVarSuffix,
                        "JetsAK8_axisminor_1"+label+myVarSuffix,      "JetsAK8_axisminor_2"+label+myVarSuffix,      "JetsAK8_axisminor_3"+label+myVarSuffix,      "JetsAK8_axisminor_4"+label+myVarSuffix,      "JetsAK8_axisminor_5"+label+myVarSuffix,
                        "JetsAK8_nsubjets_1"+label+myVarSuffix,       "JetsAK8_nsubjets_2"+label+myVarSuffix,       "JetsAK8_nsubjets_3"+label+myVarSuffix,       "JetsAK8_nsubjets_4"+label+myVarSuffix,       "JetsAK8_nsubjets_5"+label+myVarSuffix,
                        "JetsAK8_tDiscriminator_1"+label+myVarSuffix, "JetsAK8_tDiscriminator_2"+label+myVarSuffix, "JetsAK8_tDiscriminator_3"+label+myVarSuffix, "JetsAK8_tDiscriminator_4"+label+myVarSuffix, "JetsAK8_tDiscriminator_5"+label+myVarSuffix,
                        "JetsAK8_wDiscriminator_1"+label+myVarSuffix, "JetsAK8_wDiscriminator_2"+label+myVarSuffix, "JetsAK8_wDiscriminator_3"+label+myVarSuffix, "JetsAK8_wDiscriminator_4"+label+myVarSuffix, "JetsAK8_wDiscriminator_5"+label+myVarSuffix,
                        "JetsAK8_hDiscriminator_1"+label+myVarSuffix, "JetsAK8_hDiscriminator_2"+label+myVarSuffix, "JetsAK8_hDiscriminator_3"+label+myVarSuffix, "JetsAK8_hDiscriminator_4"+label+myVarSuffix, "JetsAK8_hDiscriminator_5"+label+myVarSuffix,
                        "JetsAK8_multiplicity_1"+label+myVarSuffix,   "JetsAK8_multiplicity_2"+label+myVarSuffix,   "JetsAK8_multiplicity_3"+label+myVarSuffix,   "JetsAK8_multiplicity_4"+label+myVarSuffix,   "JetsAK8_multiplicity_5"+label+myVarSuffix,
                    };

                    for (const auto& split : myTree)
                    {
                        my_histos["EventCounter"+split.first+label]->Fill( eventCounter );
                        myMiniTuple[split.first][label][myVarSuffix]->setTupleVars(varGeneral);
                        myMiniTuple[split.first][label][myVarSuffix]->setTupleVars(varEventShape); 
                        myMiniTuple[split.first][label][myVarSuffix]->setTupleVars(varJets);
                        myMiniTuple[split.first][label][myVarSuffix]->setTupleVars(varJetsAK8);
                        myMiniTuple[split.first][label][myVarSuffix]->setTupleVars(varOldSeed);
                        myMiniTuple[split.first][label][myVarSuffix]->initBranches(tr);
                    }
                }
            }
            
            //-----------------------------------
            //-- Fill Histograms Below
            //-----------------------------------
            for (const auto& selection : baselines)
            {
                if (selection.second)
                {
                    int mod = counts["count"][selection.first][myVarSuffix] % 10;
                    if(mod < 8)
                    {
                        myMiniTuple["Train"][selection.first][myVarSuffix]->fill();
                        counts["Train"][selection.first][myVarSuffix]++;
                    }
                    else if(mod == 8)
                    {
                        myMiniTuple["Test"][selection.first][myVarSuffix]->fill();
                        counts["Test"][selection.first][myVarSuffix]++;
                    }
                    else
                    {
                        myMiniTuple["Val"][selection.first][myVarSuffix]->fill();
                        counts["Val"][selection.first][myVarSuffix]++;
                    }
                    counts["count"][selection.first][myVarSuffix]++;
                }
            }
        }
    }

    for (const auto& train : counts)
    {
        for (const auto& channel : train.second)
        {
            for (const auto& suffix : channel.second)
            {
                std::cout << train.first+channel.first+suffix.first+": " << suffix.second << std::endl;
            }
        }
    }

}//END of function
      
void MakeNNVariables::WriteHistos( TFile* outfile ) 
{
    const auto& outFileName = std::string(outfile->GetName());
    const auto& name = utility::split("first", outFileName, ".");

    for (const auto& split : myTree)
    {

        TFile* outfileTrain = TFile::Open((name+"_"+split.first+".root").c_str(), "RECREATE");
        outfileTrain->cd();

        for (auto& channel : split.second)
        { 
            for (auto& suffix : channel.second)
            {
                suffix.second->Write();

                delete suffix.second;
                delete myMiniTuple[split.first][channel.first][suffix.first];
            }
        }

        for (const auto& histo : my_histos)
        {
            if (histo.first == split.first)
                histo.second->Write();
        }

        outfileTrain->Close();
    }

    remove(outFileName.c_str());
}
