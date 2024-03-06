#define MakeAnaSkimTree_cxx
#include "Analyzer/Analyzer/include/MakeAnaSkimTree.h"
#include "NTupleReader/include/NTupleReader.h"
#include "Framework/Framework/include/MiniTupleMaker.h"

#include "Framework/Framework/include/Muon.h"
#include "Framework/Framework/include/Electron.h"
#include "Framework/Framework/include/Jet.h"
#include "Framework/Framework/include/BJet.h"
#include "Framework/Framework/include/CommonVariables.h"
#include "Framework/Framework/include/RunTopTagger.h"
#include "Framework/Framework/include/MakeMVAVariables.h"
#include "Framework/Framework/include/Baseline.h"
#include "Framework/Framework/include/BTagCorrector.h"
#include "Framework/Framework/include/ScaleFactors.h"
#include "Framework/Framework/include/StopJets.h"
#include "Framework/Framework/include/StopGenMatch.h"
#include "Framework/Framework/include/MakeStopHemispheres.h"
#include "Framework/Framework/include/DeepEventShape.h"

#include <iostream>

MakeAnaSkimTree::MakeAnaSkimTree()
{
    InitHistos();

    jecvars  = {"", "JECup", "JECdown", "JERup", "JERdown"};

    ttvars = {"erdON", "hdamp", "TuneCP5"};

    for (const auto& jecvar : jecvars)
    {
        treeInits[jecvar] = false;
    }
}

void MakeAnaSkimTree::InitHistos()
{
    eventCounter = std::make_shared<TH1D>("EventCounter", "EventCounter", 2, -1.1, 1.1);
}

void MakeAnaSkimTree::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    const auto& filetag                           = tr.getVar<std::string>("filetag"                          );
    const auto& runtype                           = tr.getVar<std::string>("runtype"                          );
    const auto& runYear                           = tr.getVar<std::string>("runYear"                          );
    const auto& btagEffFileName                   = tr.getVar<std::string>("btagEffFileName"                  );
    const auto& bjetTagFileName                   = tr.getVar<std::string>("bjetTagFileName"                  );
    const auto& leptonFileName                    = tr.getVar<std::string>("leptonFileName"                   );
    const auto& hadronicFileName                  = tr.getVar<std::string>("hadronicFileName"                 );
    const auto& toptaggerFileName                 = tr.getVar<std::string>("toptaggerFileName"                );
    const auto& meanFileName                      = tr.getVar<std::string>("meanFileName"                     );
    const auto& TopTaggerCfg                      = tr.getVar<std::string>("TopTaggerCfg"                     );
    const auto& DoubleDisCo_Cfg_0l_RPV            = tr.getVar<std::string>("DoubleDisCo_Cfg_0l_RPV"           );
    const auto& DoubleDisCo_Model_0l_RPV          = tr.getVar<std::string>("DoubleDisCo_Model_0l_RPV"         );
    const auto& DoubleDisCo_Cfg_NonIsoMuon_0l_RPV = tr.getVar<std::string>("DoubleDisCo_Cfg_NonIsoMuon_0l_RPV");
    const auto& DoubleDisCo_Cfg_1l_RPV            = tr.getVar<std::string>("DoubleDisCo_Cfg_1l_RPV"           );
    const auto& DoubleDisCo_Model_1l_RPV          = tr.getVar<std::string>("DoubleDisCo_Model_1l_RPV"         );
    const auto& DoubleDisCo_Cfg_NonIsoMuon_1l_RPV = tr.getVar<std::string>("DoubleDisCo_Cfg_NonIsoMuon_1l_RPV");
    const auto& DoubleDisCo_Cfg_2l_RPV            = tr.getVar<std::string>("DoubleDisCo_Cfg_2l_RPV"           );
    const auto& DoubleDisCo_Model_2l_RPV          = tr.getVar<std::string>("DoubleDisCo_Model_2l_RPV"         );
    const auto& DoubleDisCo_Cfg_NonIsoMuon_2l_RPV = tr.getVar<std::string>("DoubleDisCo_Cfg_NonIsoMuon_2l_RPV");
    const auto& DoubleDisCo_Cfg_0l_SYY            = tr.getVar<std::string>("DoubleDisCo_Cfg_0l_SYY"           );
    const auto& DoubleDisCo_Model_0l_SYY          = tr.getVar<std::string>("DoubleDisCo_Model_0l_SYY"         );
    const auto& DoubleDisCo_Cfg_NonIsoMuon_0l_SYY = tr.getVar<std::string>("DoubleDisCo_Cfg_NonIsoMuon_0l_SYY");
    const auto& DoubleDisCo_Cfg_1l_SYY            = tr.getVar<std::string>("DoubleDisCo_Cfg_1l_SYY"           );
    const auto& DoubleDisCo_Model_1l_SYY          = tr.getVar<std::string>("DoubleDisCo_Model_1l_SYY"         );
    const auto& DoubleDisCo_Cfg_NonIsoMuon_1l_SYY = tr.getVar<std::string>("DoubleDisCo_Cfg_NonIsoMuon_1l_SYY");
    const auto& DoubleDisCo_Cfg_2l_SYY            = tr.getVar<std::string>("DoubleDisCo_Cfg_2l_SYY"           );
    const auto& DoubleDisCo_Model_2l_SYY          = tr.getVar<std::string>("DoubleDisCo_Model_2l_SYY"         );
    const auto& DoubleDisCo_Cfg_NonIsoMuon_2l_SYY = tr.getVar<std::string>("DoubleDisCo_Cfg_NonIsoMuon_2l_SYY");

    bool runningOnTTvar = false;
    for (const auto& ttvar : ttvars)
    {
        if ( filetag.find(ttvar) != std::string::npos )
        {
            runningOnTTvar = true;
            break;
        }
    }

    // Cannot do JEC and JER variations for data
    // and no need for JEC/JER on top of TT var
    if (runtype == "Data" or runningOnTTvar)
        jecvars = {""};

    for(const auto& jecvar : jecvars)
    {

        Jet                 jet(jecvar);
        BJet                bjet(jecvar);
        Muon                muon(jecvar);
        Baseline            baseline(jecvar);
        Electron            electron(jecvar);
        StopJets            stopJets(jecvar);
        RunTopTagger        topTagger(TopTaggerCfg, jecvar);
        StopGenMatch        stopGenMatch(jecvar);
        CommonVariables     commonVariables(jecvar);
        MakeMVAVariables    makeMVAVariables(                false,  jecvar,        "GoodJets_pt30",       false, true, 7, 2, ""                                  );
        MakeMVAVariables    makeMVAVariables_NonIsoMuon(     false,  jecvar,        "NonIsoMuonJets_pt30", false, true, 7, 2, ""                                  );
        // 0l
        // note that if we make the inputs to the NN, we use just the GoodJets_pt30 collection to derive things
        // but, if we define the QCD CR selection, we use the NonIsoMuonJets_pt30 collection
        DeepEventShape      neuralNetwork0L_RPV(           DoubleDisCo_Cfg_0l_RPV,            DoubleDisCo_Model_0l_RPV, "Info", true, jecvar                          );
        DeepEventShape      neuralNetwork0L_NonIsoMuon_RPV(DoubleDisCo_Cfg_NonIsoMuon_0l_RPV, DoubleDisCo_Model_0l_RPV, "Info", true, jecvar                          );
        DeepEventShape      neuralNetwork0L_SYY(           DoubleDisCo_Cfg_0l_SYY,            DoubleDisCo_Model_0l_SYY, "Info", true, jecvar                          );
        DeepEventShape      neuralNetwork0L_NonIsoMuon_SYY(DoubleDisCo_Cfg_NonIsoMuon_0l_SYY, DoubleDisCo_Model_0l_SYY, "Info", true, jecvar                          );
        MakeStopHemispheres stopHemispheres_TopSeed(           "StopJets", "GoodStopJets", "NGoodStopJets", "_TopSeed",            jecvar, Hemisphere::TopSeed    );
        MakeStopHemispheres stopHemispheres_TopSeed_NonIsoMuon("StopJets", "GoodStopJets", "NGoodStopJets", "_TopSeed_NonIsoMuon", jecvar, Hemisphere::InvMassSeed);
        // 1l
        DeepEventShape      neuralNetwork1L_RPV(           DoubleDisCo_Cfg_1l_RPV,            DoubleDisCo_Model_1l_RPV, "Info", true, jecvar                                    );
        DeepEventShape      neuralNetwork1L_NonIsoMuon_RPV(DoubleDisCo_Cfg_NonIsoMuon_1l_RPV, DoubleDisCo_Model_1l_RPV, "Info", true, jecvar                                    );
        DeepEventShape      neuralNetwork1L_SYY(           DoubleDisCo_Cfg_1l_SYY,            DoubleDisCo_Model_1l_SYY, "Info", true, jecvar                                    );
        DeepEventShape      neuralNetwork1L_NonIsoMuon_SYY(DoubleDisCo_Cfg_NonIsoMuon_1l_SYY, DoubleDisCo_Model_1l_SYY, "Info", true, jecvar                                    );
        MakeStopHemispheres stopHemispheres_OldSeed(           "Jets", "GoodJets_pt20",       "NGoodJets_pt20",       "_OldSeed",            jecvar, Hemisphere::InvMassSeed);
        MakeStopHemispheres stopHemispheres_OldSeed_NonIsoMuon("Jets", "NonIsoMuonJets_pt20", "NNonIsoMuonJets_pt30", "_OldSeed_NonIsoMuon", jecvar, Hemisphere::InvMassSeed);
        // 2l
        DeepEventShape      neuralNetwork2L_RPV(           DoubleDisCo_Cfg_2l_RPV,            DoubleDisCo_Model_2l_RPV, "Info", true, jecvar);
        DeepEventShape      neuralNetwork2L_NonIsoMuon_RPV(DoubleDisCo_Cfg_NonIsoMuon_2l_RPV, DoubleDisCo_Model_2l_RPV, "Info", true, jecvar);
        DeepEventShape      neuralNetwork2L_SYY(           DoubleDisCo_Cfg_2l_SYY,            DoubleDisCo_Model_2l_SYY, "Info", true, jecvar);
        DeepEventShape      neuralNetwork2L_NonIsoMuon_SYY(DoubleDisCo_Cfg_NonIsoMuon_2l_SYY, DoubleDisCo_Model_2l_SYY, "Info", true, jecvar);

        // Remember, order matters here !
        // Follow what is done in Config.h
        tr.registerFunction(muon);
        tr.registerFunction(electron);
        tr.registerFunction(jet);
        tr.registerFunction(bjet);
        tr.registerFunction(commonVariables);
        tr.registerFunction(topTagger);
        tr.registerFunction(baseline);
        tr.registerFunction(makeMVAVariables);
        tr.registerFunction(makeMVAVariables_NonIsoMuon);
        tr.registerFunction(stopJets);
        tr.registerFunction(stopHemispheres_TopSeed);
        tr.registerFunction(stopHemispheres_OldSeed);
        tr.registerFunction(stopHemispheres_TopSeed_NonIsoMuon);
        tr.registerFunction(stopHemispheres_OldSeed_NonIsoMuon);

        if (runtype == "MC")
        {
            ScaleFactors        scaleFactors(runYear, leptonFileName, hadronicFileName, toptaggerFileName, meanFileName, filetag, jecvar);
            BTagCorrector       bTagCorrector(btagEffFileName, "", bjetTagFileName, "", filetag);
            bTagCorrector.SetVarNames("GenParticles_PdgId", "Jets"+jecvar, "GoodJets_pt30"+jecvar, "Jets"+jecvar+"_bJetTagDeepFlavourtotb", "Jets"+jecvar+"_partonFlavor", jecvar);

            tr.registerFunction(bTagCorrector);
            tr.registerFunction(scaleFactors);
            tr.registerFunction(stopGenMatch);
       }

        tr.registerFunction(neuralNetwork0L_RPV);
        tr.registerFunction(neuralNetwork0L_NonIsoMuon_RPV);
        tr.registerFunction(neuralNetwork0L_SYY);
        tr.registerFunction(neuralNetwork0L_NonIsoMuon_SYY);

        tr.registerFunction(neuralNetwork1L_RPV);
        tr.registerFunction(neuralNetwork1L_NonIsoMuon_RPV);
        tr.registerFunction(neuralNetwork1L_SYY);
        tr.registerFunction(neuralNetwork1L_NonIsoMuon_SYY);

        tr.registerFunction(neuralNetwork2L_RPV);
        tr.registerFunction(neuralNetwork2L_NonIsoMuon_RPV);
        tr.registerFunction(neuralNetwork2L_SYY);
        tr.registerFunction(neuralNetwork2L_NonIsoMuon_SYY);
    }

    while( tr.getNextEvent() )
    {
        //------------------------------------
        //-- Print Event Number
        //------------------------------------
        if( maxevents != -1 && tr.getEvtNum() > maxevents )
            break;
        if( tr.getEvtNum() % 1000 == 0 )
            printf( " Event %i\n", tr.getEvtNum() );

        const auto& evtCounter = tr.getVar<int>("eventCounter");
        eventCounter->Fill(evtCounter);

        const auto& runtype = tr.getVar<std::string>("runtype");
        const auto& runYear = tr.getVar<std::string>("runYear"); 

        auto& year = tr.createDerivedVar<int>("year", 0);
        if      (runYear == "2016preVFP")
            year = 160;
        else if (runYear == "2016postVFP")
            year = 161;
        else if (runYear == "2017")
            year = 17;
        else if (runYear == "2018")
            year = 18;

        for(const auto& jecvar : jecvars)
        {
            if( !treeInits[jecvar] ) {
                std::set<std::string> variables = {
                    "year",
                    "NGoodJets_pt30"+jecvar,
                    "NNonIsoMuonJets_pt30"+jecvar,
                    "DoubleDisCo_disc1_0l_RPV"+jecvar,
                    "DoubleDisCo_disc2_0l_RPV"+jecvar,
                    "DoubleDisCo_disc1_0l_SYY"+jecvar,
                    "DoubleDisCo_disc2_0l_SYY"+jecvar,
                    "DoubleDisCo_disc1_1l_RPV"+jecvar,
                    "DoubleDisCo_disc2_1l_RPV"+jecvar,
                    "DoubleDisCo_disc1_1l_SYY"+jecvar,
                    "DoubleDisCo_disc2_1l_SYY"+jecvar,
                    "DoubleDisCo_disc1_2l_RPV"+jecvar,
                    "DoubleDisCo_disc2_2l_RPV"+jecvar,
                    "DoubleDisCo_disc1_2l_SYY"+jecvar,
                    "DoubleDisCo_disc2_2l_SYY"+jecvar,
                    "DoubleDisCo_massReg_0l_SYY"+jecvar,
                    "DoubleDisCo_massReg_0l_RPV"+jecvar,
                    "DoubleDisCo_massReg_1l_SYY"+jecvar,
                    "DoubleDisCo_massReg_1l_RPV"+jecvar,
                    "DoubleDisCo_massReg_2l_SYY"+jecvar,
                    "DoubleDisCo_massReg_2l_RPV"+jecvar,
                    "DoubleDisCo_disc1_NonIsoMuon_0l_RPV"+jecvar,
                    "DoubleDisCo_disc2_NonIsoMuon_0l_RPV"+jecvar,
                    "DoubleDisCo_disc1_NonIsoMuon_0l_SYY"+jecvar,
                    "DoubleDisCo_disc2_NonIsoMuon_0l_SYY"+jecvar,
                    "DoubleDisCo_disc1_NonIsoMuon_1l_RPV"+jecvar,
                    "DoubleDisCo_disc2_NonIsoMuon_1l_RPV"+jecvar,
                    "DoubleDisCo_disc1_NonIsoMuon_1l_SYY"+jecvar,
                    "DoubleDisCo_disc2_NonIsoMuon_1l_SYY"+jecvar,
                    "DoubleDisCo_disc1_NonIsoMuon_2l_RPV"+jecvar,
                    "DoubleDisCo_disc2_NonIsoMuon_2l_RPV"+jecvar,
                    "DoubleDisCo_disc1_NonIsoMuon_2l_SYY"+jecvar,
                    "DoubleDisCo_disc2_NonIsoMuon_2l_SYY"+jecvar,
                    "DoubleDisCo_massReg_NonIsoMuon_0l_SYY"+jecvar,
                    "DoubleDisCo_massReg_NonIsoMuon_0l_RPV"+jecvar,
                    "DoubleDisCo_massReg_NonIsoMuon_1l_SYY"+jecvar,
                    "DoubleDisCo_massReg_NonIsoMuon_1l_RPV"+jecvar,
                    "DoubleDisCo_massReg_NonIsoMuon_2l_SYY"+jecvar,
                    "DoubleDisCo_massReg_NonIsoMuon_2l_RPV"+jecvar
                };

                if (runtype == "Data")
                {
                    variables.insert("RunNum");
                    variables.insert("EvtNum");
                }

                if (jecvar == "")
                {
                    if (runtype == "MC")
                    {
                        variables.insert("FinalLumi");
                        variables.insert("Weight");
                    }

                    variables.insert("combined6thToLastJet_pt_cm");
                    variables.insert("combined6thToLastJet_eta_cm");
                    variables.insert("combined6thToLastJet_phi_cm");
                    variables.insert("combined6thToLastJet_m_cm");
                    variables.insert("combined6thToLastJet_E_cm");
                    variables.insert("combined7thToLastJet_pt_cm");
                    variables.insert("combined7thToLastJet_eta_cm");
                    variables.insert("combined7thToLastJet_phi_cm");
                    variables.insert("combined7thToLastJet_m_cm");
                    variables.insert("combined7thToLastJet_E_cm");
                    variables.insert("combined8thToLastJet_pt_cm");
                    variables.insert("combined8thToLastJet_eta_cm");
                    variables.insert("combined8thToLastJet_phi_cm");
                    variables.insert("combined8thToLastJet_m_cm");
                    variables.insert("combined8thToLastJet_E_cm");
                    variables.insert("combined6thToLastJetNonIsoMuons_pt_cm");
                    variables.insert("combined6thToLastJetNonIsoMuons_eta_cm");
                    variables.insert("combined6thToLastJetNonIsoMuons_phi_cm");
                    variables.insert("combined6thToLastJetNonIsoMuons_m_cm");
                    variables.insert("combined6thToLastJetNonIsoMuons_E_cm");
                    variables.insert("combined7thToLastJetNonIsoMuons_pt_cm");
                    variables.insert("combined7thToLastJetNonIsoMuons_eta_cm");
                    variables.insert("combined7thToLastJetNonIsoMuons_phi_cm");
                    variables.insert("combined7thToLastJetNonIsoMuons_m_cm");
                    variables.insert("combined7thToLastJetNonIsoMuons_E_cm");
                    variables.insert("combined8thToLastJetNonIsoMuons_pt_cm");
                    variables.insert("combined8thToLastJetNonIsoMuons_eta_cm");
                    variables.insert("combined8thToLastJetNonIsoMuons_phi_cm");
                    variables.insert("combined8thToLastJetNonIsoMuons_m_cm");
                    variables.insert("combined8thToLastJetNonIsoMuons_E_cm");
                    variables.insert("dR_bjets");
                    variables.insert("event_beta_z");
                    variables.insert("event_phi_rotate");
                    variables.insert("fixedGridRhoFastjetAll"); 
                    variables.insert("fwm2_top6");
                    variables.insert("fwm3_top6");
                    variables.insert("fwm4_top6");
                    variables.insert("fwm5_top6");
                    variables.insert("NonIsoMuons_fwm2_top6");
                    variables.insert("NonIsoMuons_fwm3_top6");
                    variables.insert("NonIsoMuons_fwm4_top6");
                    variables.insert("NonIsoMuons_fwm5_top6");
                    variables.insert("GoodLeptons_m_1");    variables.insert("GoodLeptons_m_2");    variables.insert("GoodNonIsoMuons_m_1");   
                    variables.insert("GoodLeptons_eta_1");  variables.insert("GoodLeptons_eta_2");  variables.insert("GoodNonIsoMuons_eta_1"); 
                    variables.insert("GoodLeptons_phi_1");  variables.insert("GoodLeptons_phi_2");  variables.insert("GoodNonIsoMuons_phi_1"); 
                    variables.insert("GoodLeptons_pt_1");   variables.insert("GoodLeptons_pt_2");   variables.insert("GoodNonIsoMuons_pt_1");  
                    variables.insert("GoodLeptons_flav_1"); variables.insert("GoodLeptons_flav_2"); variables.insert("GoodNonIsoMuons_flav_1");
                    variables.insert("GoodLeptons_iso_1");  variables.insert("GoodLeptons_iso_2");  variables.insert("GoodNonIsoMuons_iso_1"); 
                    variables.insert("HT_trigger_pt30");
                    variables.insert("HT_NonIsoMuon_pt30");
                    variables.insert("HT_trigger_pt45");
                    variables.insert("jmt_ev0_top6");
                    variables.insert("jmt_ev1_top6");
                    variables.insert("jmt_ev2_top6");
                    variables.insert("NonIsoMuons_jmt_ev0_top6");
                    variables.insert("NonIsoMuons_jmt_ev1_top6");
                    variables.insert("NonIsoMuons_jmt_ev2_top6");
                    variables.insert("Jet_m_1");       variables.insert("Jet_m_2");       variables.insert("Jet_m_3");       variables.insert("Jet_m_4");       variables.insert("Jet_m_5");       variables.insert("Jet_m_6");       variables.insert("Jet_m_7");
                    variables.insert("Jet_E_1");       variables.insert("Jet_E_2");       variables.insert("Jet_E_3");       variables.insert("Jet_E_4");       variables.insert("Jet_E_5");       variables.insert("Jet_E_6");       variables.insert("Jet_E_7");
                    variables.insert("Jet_eta_1");     variables.insert("Jet_eta_2");     variables.insert("Jet_eta_3");     variables.insert("Jet_eta_4");     variables.insert("Jet_eta_5");     variables.insert("Jet_eta_6");     variables.insert("Jet_eta_7");
                    variables.insert("Jet_phi_1");     variables.insert("Jet_phi_2");     variables.insert("Jet_phi_3");     variables.insert("Jet_phi_4");     variables.insert("Jet_phi_5");     variables.insert("Jet_phi_6");     variables.insert("Jet_phi_7");
                    variables.insert("Jet_pt_1");      variables.insert("Jet_pt_2");      variables.insert("Jet_pt_3");      variables.insert("Jet_pt_4");      variables.insert("Jet_pt_5");      variables.insert("Jet_pt_6");      variables.insert("Jet_pt_7");
                    variables.insert("Jet_flavb_1");   variables.insert("Jet_flavb_2");   variables.insert("Jet_flavb_3");   variables.insert("Jet_flavb_4");   variables.insert("Jet_flavb_5");   variables.insert("Jet_flavb_6");   variables.insert("Jet_flavb_7");
                    variables.insert("Jet_flavg_1");   variables.insert("Jet_flavg_2");   variables.insert("Jet_flavg_3");   variables.insert("Jet_flavg_4");   variables.insert("Jet_flavg_5");   variables.insert("Jet_flavg_6");   variables.insert("Jet_flavg_7");
                    variables.insert("Jet_flavc_1");   variables.insert("Jet_flavc_2");   variables.insert("Jet_flavc_3");   variables.insert("Jet_flavc_4");   variables.insert("Jet_flavc_5");   variables.insert("Jet_flavc_6");   variables.insert("Jet_flavc_7");
                    variables.insert("Jet_flavuds_1"); variables.insert("Jet_flavuds_2"); variables.insert("Jet_flavuds_3"); variables.insert("Jet_flavuds_4"); variables.insert("Jet_flavuds_5"); variables.insert("Jet_flavuds_6"); variables.insert("Jet_flavuds_7");
                    variables.insert("Jet_flavq_1");   variables.insert("Jet_flavq_2");   variables.insert("Jet_flavq_3");   variables.insert("Jet_flavq_4");   variables.insert("Jet_flavq_5");   variables.insert("Jet_flavq_6");   variables.insert("Jet_flavq_7");
                    variables.insert("JetNonIsoMuons_m_1");       variables.insert("JetNonIsoMuons_m_2");       variables.insert("JetNonIsoMuons_m_3");       variables.insert("JetNonIsoMuons_m_4");       variables.insert("JetNonIsoMuons_m_5");       variables.insert("JetNonIsoMuons_m_6");       variables.insert("JetNonIsoMuons_m_7");
                    variables.insert("JetNonIsoMuons_E_1");       variables.insert("JetNonIsoMuons_E_2");       variables.insert("JetNonIsoMuons_E_3");       variables.insert("JetNonIsoMuons_E_4");       variables.insert("JetNonIsoMuons_E_5");       variables.insert("JetNonIsoMuons_E_6");       variables.insert("JetNonIsoMuons_E_7");
                    variables.insert("JetNonIsoMuons_eta_1");     variables.insert("JetNonIsoMuons_eta_2");     variables.insert("JetNonIsoMuons_eta_3");     variables.insert("JetNonIsoMuons_eta_4");     variables.insert("JetNonIsoMuons_eta_5");     variables.insert("JetNonIsoMuons_eta_6");     variables.insert("JetNonIsoMuons_eta_7");
                    variables.insert("JetNonIsoMuons_phi_1");     variables.insert("JetNonIsoMuons_phi_2");     variables.insert("JetNonIsoMuons_phi_3");     variables.insert("JetNonIsoMuons_phi_4");     variables.insert("JetNonIsoMuons_phi_5");     variables.insert("JetNonIsoMuons_phi_6");     variables.insert("JetNonIsoMuons_phi_7");
                    variables.insert("JetNonIsoMuons_pt_1");      variables.insert("JetNonIsoMuons_pt_2");      variables.insert("JetNonIsoMuons_pt_3");      variables.insert("JetNonIsoMuons_pt_4");      variables.insert("JetNonIsoMuons_pt_5");      variables.insert("JetNonIsoMuons_pt_6");      variables.insert("JetNonIsoMuons_pt_7");
                    variables.insert("JetNonIsoMuons_flavb_1");   variables.insert("JetNonIsoMuons_flavb_2");   variables.insert("JetNonIsoMuons_flavb_3");   variables.insert("JetNonIsoMuons_flavb_4");   variables.insert("JetNonIsoMuons_flavb_5");   variables.insert("JetNonIsoMuons_flavb_6");   variables.insert("JetNonIsoMuons_flavb_7");
                    variables.insert("JetNonIsoMuons_flavg_1");   variables.insert("JetNonIsoMuons_flavg_2");   variables.insert("JetNonIsoMuons_flavg_3");   variables.insert("JetNonIsoMuons_flavg_4");   variables.insert("JetNonIsoMuons_flavg_5");   variables.insert("JetNonIsoMuons_flavg_6");   variables.insert("JetNonIsoMuons_flavg_7");
                    variables.insert("JetNonIsoMuons_flavc_1");   variables.insert("JetNonIsoMuons_flavc_2");   variables.insert("JetNonIsoMuons_flavc_3");   variables.insert("JetNonIsoMuons_flavc_4");   variables.insert("JetNonIsoMuons_flavc_5");   variables.insert("JetNonIsoMuons_flavc_6");   variables.insert("JetNonIsoMuons_flavc_7");
                    variables.insert("JetNonIsoMuons_flavuds_1"); variables.insert("JetNonIsoMuons_flavuds_2"); variables.insert("JetNonIsoMuons_flavuds_3"); variables.insert("JetNonIsoMuons_flavuds_4"); variables.insert("JetNonIsoMuons_flavuds_5"); variables.insert("JetNonIsoMuons_flavuds_6"); variables.insert("JetNonIsoMuons_flavuds_7");
                    variables.insert("JetNonIsoMuons_flavq_1");   variables.insert("JetNonIsoMuons_flavq_2");   variables.insert("JetNonIsoMuons_flavq_3");   variables.insert("JetNonIsoMuons_flavq_4");   variables.insert("JetNonIsoMuons_flavq_5");   variables.insert("JetNonIsoMuons_flavq_6");   variables.insert("JetNonIsoMuons_flavq_7");
                    variables.insert("Mbl");
                    variables.insert("Mbb");
                    variables.insert("MET");
                    variables.insert("METPhi");
                    variables.insert("mll");
                    variables.insert("NGoodBJets_pt30");
                    variables.insert("NNonIsoMuonBJets_pt30");
                    variables.insert("NGoodJets_pt45");
                    variables.insert("NGoodBJets_pt45");
                    variables.insert("ntops_1jet");
                    variables.insert("ntops_3jet");
                    variables.insert("NVtx");
                    variables.insert("passBaseline0l_Good");
                    variables.insert("passBaseline1l_Good");
                    variables.insert("passBaseline2l_Good");
                    variables.insert("pass_qcdCR_0l");
                    variables.insert("pass_qcdCR_1l");
                    variables.insert("pass_qcdCR_2l");
                    variables.insert("Stop1_mass_cm_OldSeed"); variables.insert("Stop2_mass_cm_OldSeed");
                    variables.insert("Stop1_pt_cm_OldSeed");   variables.insert("Stop2_pt_cm_OldSeed");
                    variables.insert("Stop1_phi_cm_OldSeed");  variables.insert("Stop2_phi_cm_OldSeed");
                    variables.insert("Stop1_eta_cm_OldSeed");  variables.insert("Stop2_eta_cm_OldSeed");
                    variables.insert("Stop1_mass_cm_OldSeed_NonIsoMuon"); variables.insert("Stop2_mass_cm_OldSeed_NonIsoMuon");
                    variables.insert("Stop1_pt_cm_OldSeed_NonIsoMuon");   variables.insert("Stop2_pt_cm_OldSeed_NonIsoMuon");
                    variables.insert("Stop1_phi_cm_OldSeed_NonIsoMuon");  variables.insert("Stop2_phi_cm_OldSeed_NonIsoMuon");
                    variables.insert("Stop1_eta_cm_OldSeed_NonIsoMuon");  variables.insert("Stop2_eta_cm_OldSeed_NonIsoMuon");
                    variables.insert("Stop1_mass_cm_TopSeed"); variables.insert("Stop2_mass_cm_TopSeed");
                    variables.insert("Stop1_pt_cm_TopSeed");   variables.insert("Stop2_pt_cm_TopSeed");
                    variables.insert("Stop1_phi_cm_TopSeed");  variables.insert("Stop2_phi_cm_TopSeed");
                    variables.insert("Stop1_eta_cm_TopSeed");  variables.insert("Stop2_eta_cm_TopSeed");
                    variables.insert("top1_pt_cm");            variables.insert("top2_pt_cm");
                    variables.insert("top1_eta_cm");           variables.insert("top2_eta_cm");
                    variables.insert("top1_phi_cm");           variables.insert("top2_phi_cm");
                    variables.insert("top1_mass_cm");          variables.insert("top2_mass_cm");
                }

                if ( runtype == "MC" )
                {
                    variables.insert("TotalWeight_0l"+jecvar);
                    variables.insert("TotalWeight_1l"+jecvar);
                    variables.insert("TotalWeight_2l"+jecvar);
                    variables.insert("TotalWeight_QCDCR"+jecvar);
                    if (jecvar == "")
                    {
                        variables.insert("stop1_ptrank_mass");
                        variables.insert("stop2_ptrank_mass");

                        if (! runningOnTTvar)
                        {
                            variables.insert("jetTrigSF"+jecvar);
                            variables.insert("totGoodElectronSF"+jecvar);
                            variables.insert("totGoodMuonSF"+jecvar);
                            variables.insert("totNonIsoMuonSF"+jecvar);
                            variables.insert("topTaggerScaleFactor"+jecvar);
                            variables.insert("bTagSF_EventWeightSimple_Central"+jecvar);
                            variables.insert("scaleWeightUp");
                            variables.insert("scaleWeightDown");
                            variables.insert("PSweight_ISRUp");
                            variables.insert("PSweight_ISRDown");
                            variables.insert("PSweight_FSRUp");
                            variables.insert("PSweight_FSRDown");
                            variables.insert("PDFweightUp");
                            variables.insert("PDFweightDown");
                            variables.insert("puSysUpCorr");
                            variables.insert("puSysDownCorr");
                            variables.insert("puWeightCorr");
                            variables.insert("prefiringScaleFactor");
                            variables.insert("prefiringScaleFactorUp");
                            variables.insert("prefiringScaleFactorDown");
                            variables.insert("jetTrigSF_Up");
                            variables.insert("jetTrigSF_Down");
                            variables.insert("totGoodElectronSF_Up");
                            variables.insert("totGoodElectronSF_Down");
                            variables.insert("totGoodMuonSF_Up");
                            variables.insert("totGoodMuonSF_Down");
                            variables.insert("totNonIsoMuonSF_Up");
                            variables.insert("totNonIsoMuonSF_Down");
                            variables.insert("topTaggerScaleFactorUp");
                            variables.insert("topTaggerScaleFactorDown");
                            variables.insert("bTagSF_EventWeightSimple_Up");
                            variables.insert("bTagSF_EventWeightSimple_Down");
                        }
                    }
                }   

                std::string myTreeName = "AnaSkim"+jecvar;
                myTrees[jecvar] = new TTree( (myTreeName).c_str() , (myTreeName).c_str() );
                myAnaSkimTuples[jecvar] = new MiniTupleMaker( myTrees[jecvar] );
                myAnaSkimTuples[jecvar]->setTupleVars(variables);
                myAnaSkimTuples[jecvar]->initBranches(tr);

                treeInits[jecvar] = true;
            }

            // If an event is not interesting for any channel selection (see explanation and definition
            // of lostCauseEvent boolean in Baseline.h), then move on here and do not waste any more
            // time on looping over histos or getting vars.
            // This only matters if the user specifies -s on the command line
            // N.B. We already filled the EventCounter histogram for our due-diligence of counting up
            // every single event.
            const auto& lostCauseEvent = tr.getVar<bool>("lostCauseEvent" + jecvar);
            const auto& fastMode       = tr.getVar<bool>("fastMode");

            if (lostCauseEvent and fastMode)
                continue;

            const auto& passBaseline0l = tr.getVar<bool>("passBaseline0l_Good" + jecvar);
            const auto& passBaseline1l = tr.getVar<bool>("passBaseline1l_Good" + jecvar);
            const auto& passBaseline2l = tr.getVar<bool>("passBaseline2l_Good" + jecvar);
            const auto& passQCDCR0l    = tr.getVar<bool>("pass_qcdCR_0l" + jecvar);
            const auto& passQCDCR1l    = tr.getVar<bool>("pass_qcdCR_1l" + jecvar);
            const auto& passQCDCR2l    = tr.getVar<bool>("pass_qcdCR_2l" + jecvar);

            //-----------------------------------
            //-- Fill Histograms Below
            //-----------------------------------
            if( passBaseline0l or passBaseline1l or passBaseline2l or passQCDCR0l or passQCDCR1l or passQCDCR2l ) {
                myAnaSkimTuples[jecvar]->fill(tr);
            }
        }
    } 
}
      
void MakeAnaSkimTree::WriteHistos( TFile* outfile ) 
{
    outfile->cd();

    eventCounter->Write();

    for (const auto& jecvar : jecvars)
    {
        myTrees[jecvar]->Write();

        delete myTrees[jecvar];    
        delete myAnaSkimTuples[jecvar];
    }
}
