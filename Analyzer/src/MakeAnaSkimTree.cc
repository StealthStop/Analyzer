#define MakeAnaSkimTree_cxx
#include "Analyzer/Analyzer/include/MakeAnaSkimTree.h"
#include "NTupleReader/include/NTupleReader.h"
#include "Framework/Framework/include/MiniTupleMaker.h"

#include <iostream>

MakeAnaSkimTree::MakeAnaSkimTree()
{
    InitHistos();

    treeInit = false;
}

void MakeAnaSkimTree::InitHistos()
{
    eventCounter = std::make_shared<TH1D>("EventCounter", "EventCounter", 2, -1.1, 1.1);
}

void MakeAnaSkimTree::Loop(NTupleReader& tr, double, int maxevents, bool)
{
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

        const auto& lostCauseEvent = tr.getVar<bool>("lostCauseEvent");
        const auto& fastMode       = tr.getVar<bool>("fastMode");

        if (lostCauseEvent and fastMode){
            continue;
        }

        const auto& runtype        = tr.getVar<std::string>("runtype");
        const auto& passBaseline0l = tr.getVar<bool>("passBaseline0l_Good");
        const auto& passBaseline1l = tr.getVar<bool>("passBaseline1l_Good");
        const auto& passBaseline2l = tr.getVar<bool>("passBaseline2l_Good");
        const auto& passQCDCR0l    = tr.getVar<bool>("pass_qcdCR_0l");
        const auto& passQCDCR1l    = tr.getVar<bool>("pass_qcdCR_1l");

        if( !treeInit ) {
            std::set<std::string> variables = {
                "combined6thToLastJet_pt_cm",
                "combined6thToLastJet_eta_cm",
                "combined6thToLastJet_phi_cm",
                "combined6thToLastJet_m_cm",
                "combined6thToLastJet_E_cm",
                "combined7thToLastJet_pt_cm",
                "combined7thToLastJet_eta_cm",
                "combined7thToLastJet_phi_cm",
                "combined7thToLastJet_m_cm",
                "combined7thToLastJet_E_cm",
                "combined8thToLastJet_pt_cm",
                "combined8thToLastJet_eta_cm",
                "combined8thToLastJet_phi_cm",
                "combined8thToLastJet_m_cm",
                "combined8thToLastJet_E_cm",
                "combined7thToLastJetNonIsoMuons_pt_cm",
                "combined7thToLastJetNonIsoMuons_eta_cm",
                "combined7thToLastJetNonIsoMuons_phi_cm",
                "combined7thToLastJetNonIsoMuons_m_cm",
                "combined7thToLastJetNonIsoMuons_E_cm",
                "DoubleDisCo_disc1_0l_RPV",
                "DoubleDisCo_disc2_0l_RPV",
                "DoubleDisCo_disc1_0l_SYY",
                "DoubleDisCo_disc2_0l_SYY",
                "DoubleDisCo_disc1_1l_RPV",
                "DoubleDisCo_disc2_1l_RPV",
                "DoubleDisCo_disc1_1l_SYY",
                "DoubleDisCo_disc2_1l_SYY",
                "DoubleDisCo_disc1_2l_RPV",
                "DoubleDisCo_disc2_2l_RPV",
                "DoubleDisCo_disc1_2l_SYY",
                "DoubleDisCo_disc2_2l_SYY",
                "DoubleDisCo_massReg_0l_SYY",
                "DoubleDisCo_massReg_0l_RPV",
                "DoubleDisCo_massReg_1l_SYY",
                "DoubleDisCo_massReg_1l_RPV",
                "DoubleDisCo_massReg_2l_SYY",
                "DoubleDisCo_massReg_2l_RPV",
                "DoubleDisCo_disc1_NonIsoMuon_0l_RPV",
                "DoubleDisCo_disc2_NonIsoMuon_0l_RPV",
                "DoubleDisCo_disc1_NonIsoMuon_0l_SYY",
                "DoubleDisCo_disc2_NonIsoMuon_0l_SYY",
                "DoubleDisCo_disc1_NonIsoMuon_1l_RPV",
                "DoubleDisCo_disc2_NonIsoMuon_1l_RPV",
                "DoubleDisCo_disc1_NonIsoMuon_1l_SYY",
                "DoubleDisCo_disc2_NonIsoMuon_1l_SYY",
                "DoubleDisCo_massReg_NonIsoMuon_0l_SYY",
                "DoubleDisCo_massReg_NonIsoMuon_0l_RPV",
                "DoubleDisCo_massReg_NonIsoMuon_1l_SYY",
                "DoubleDisCo_massReg_NonIsoMuon_1l_RPV",
                "dR_bjets",
                "event_beta_z",
                "event_phi_rotate",
                "fixedGridRhoFastjetAll", 
                "fwm2_top6",
                "fwm3_top6",
                "fwm4_top6",
                "fwm5_top6",
                "NonIsoMuons_fwm2_top6",
                "NonIsoMuons_fwm3_top6",
                "NonIsoMuons_fwm4_top6",
                "NonIsoMuons_fwm5_top6",
                "GoodLeptons_m_1",    "GoodLeptons_m_2",    "GoodNonIsoMuons_m_1",   
                "GoodLeptons_eta_1",  "GoodLeptons_eta_2",  "GoodNonIsoMuons_eta_1", 
                "GoodLeptons_phi_1",  "GoodLeptons_phi_2",  "GoodNonIsoMuons_phi_1", 
                "GoodLeptons_pt_1",   "GoodLeptons_pt_2",   "GoodNonIsoMuons_pt_1",  
                "GoodLeptons_flav_1", "GoodLeptons_flav_2", "GoodNonIsoMuons_flav_1",
                "GoodLeptons_iso_1",  "GoodLeptons_iso_2",  "GoodNonIsoMuons_iso_1", 
                "HT_trigger_pt30",
                "HT_NonIsoMuon_pt30",
                "HT_trigger_pt45",
                "jmt_ev0_top6",
                "jmt_ev1_top6",
                "jmt_ev2_top6",
                "NonIsoMuons_jmt_ev0_top6",
                "NonIsoMuons_jmt_ev1_top6",
                "NonIsoMuons_jmt_ev2_top6",
                "Jet_m_1",       "Jet_m_2",       "Jet_m_3",       "Jet_m_4",       "Jet_m_5",       "Jet_m_6",       "Jet_m_7",
                "Jet_E_1",       "Jet_E_2",       "Jet_E_3",       "Jet_E_4",       "Jet_E_5",       "Jet_E_6",       "Jet_E_7",
                "Jet_eta_1",     "Jet_eta_2",     "Jet_eta_3",     "Jet_eta_4",     "Jet_eta_5",     "Jet_eta_6",     "Jet_eta_7",
                "Jet_phi_1",     "Jet_phi_2",     "Jet_phi_3",     "Jet_phi_4",     "Jet_phi_5",     "Jet_phi_6",     "Jet_phi_7",
                "Jet_pt_1",      "Jet_pt_2",      "Jet_pt_3",      "Jet_pt_4",      "Jet_pt_5",      "Jet_pt_6",      "Jet_pt_7",
                "Jet_flavb_1",   "Jet_flavb_2",   "Jet_flavb_3",   "Jet_flavb_4",   "Jet_flavb_5",   "Jet_flavb_6",   "Jet_flavb_7",
                "Jet_flavg_1",   "Jet_flavg_2",   "Jet_flavg_3",   "Jet_flavg_4",   "Jet_flavg_5",   "Jet_flavg_6",   "Jet_flavg_7",
                "Jet_flavc_1",   "Jet_flavc_2",   "Jet_flavc_3",   "Jet_flavc_4",   "Jet_flavc_5",   "Jet_flavc_6",   "Jet_flavc_7",
                "Jet_flavuds_1", "Jet_flavuds_2", "Jet_flavuds_3", "Jet_flavuds_4", "Jet_flavuds_5", "Jet_flavuds_6", "Jet_flavuds_7",
                "Jet_flavq_1",   "Jet_flavq_2",   "Jet_flavq_3",   "Jet_flavq_4",   "Jet_flavq_5",   "Jet_flavq_6",   "Jet_flavq_7",
                "JetNonIsoMuons_m_1",       "JetNonIsoMuons_m_2",       "JetNonIsoMuons_m_3",       "JetNonIsoMuons_m_4",       "JetNonIsoMuons_m_5",       "JetNonIsoMuons_m_6",       "JetNonIsoMuons_m_7",
                "JetNonIsoMuons_E_1",       "JetNonIsoMuons_E_2",       "JetNonIsoMuons_E_3",       "JetNonIsoMuons_E_4",       "JetNonIsoMuons_E_5",       "JetNonIsoMuons_E_6",       "JetNonIsoMuons_E_7",
                "JetNonIsoMuons_eta_1",     "JetNonIsoMuons_eta_2",     "JetNonIsoMuons_eta_3",     "JetNonIsoMuons_eta_4",     "JetNonIsoMuons_eta_5",     "JetNonIsoMuons_eta_6",     "JetNonIsoMuons_eta_7",
                "JetNonIsoMuons_phi_1",     "JetNonIsoMuons_phi_2",     "JetNonIsoMuons_phi_3",     "JetNonIsoMuons_phi_4",     "JetNonIsoMuons_phi_5",     "JetNonIsoMuons_phi_6",     "JetNonIsoMuons_phi_7",
                "JetNonIsoMuons_pt_1",      "JetNonIsoMuons_pt_2",      "JetNonIsoMuons_pt_3",      "JetNonIsoMuons_pt_4",      "JetNonIsoMuons_pt_5",      "JetNonIsoMuons_pt_6",      "JetNonIsoMuons_pt_7",
                "JetNonIsoMuons_flavb_1",   "JetNonIsoMuons_flavb_2",   "JetNonIsoMuons_flavb_3",   "JetNonIsoMuons_flavb_4",   "JetNonIsoMuons_flavb_5",   "JetNonIsoMuons_flavb_6",   "JetNonIsoMuons_flavb_7",
                "JetNonIsoMuons_flavg_1",   "JetNonIsoMuons_flavg_2",   "JetNonIsoMuons_flavg_3",   "JetNonIsoMuons_flavg_4",   "JetNonIsoMuons_flavg_5",   "JetNonIsoMuons_flavg_6",   "JetNonIsoMuons_flavg_7",
                "JetNonIsoMuons_flavc_1",   "JetNonIsoMuons_flavc_2",   "JetNonIsoMuons_flavc_3",   "JetNonIsoMuons_flavc_4",   "JetNonIsoMuons_flavc_5",   "JetNonIsoMuons_flavc_6",   "JetNonIsoMuons_flavc_7",
                "JetNonIsoMuons_flavuds_1", "JetNonIsoMuons_flavuds_2", "JetNonIsoMuons_flavuds_3", "JetNonIsoMuons_flavuds_4", "JetNonIsoMuons_flavuds_5", "JetNonIsoMuons_flavuds_6", "JetNonIsoMuons_flavuds_7",
                "JetNonIsoMuons_flavq_1",   "JetNonIsoMuons_flavq_2",   "JetNonIsoMuons_flavq_3",   "JetNonIsoMuons_flavq_4",   "JetNonIsoMuons_flavq_5",   "JetNonIsoMuons_flavq_6",   "JetNonIsoMuons_flavq_7",
                "Mbl",
                "Mbb",
                "MET",
                "METPhi",
                "mll",
                "NGoodJets_pt30",
                "NGoodBJets_pt30",
                "NNonIsoMuonJets_pt30",
                "NNonIsoMuonBJets_pt30",
                "NGoodJets_pt45",
                "NGoodBJets_pt45",
                "ntops_1jet",
                "ntops_3jet",
                "NVtx",
                "passBaseline0l_Good",
                "passBaseline1l_Good",
                "passBaseline2l_Good",
                "pass_qcdCR_0l",
                "pass_qcdCR_1l",
                "Stop1_mass_cm_OldSeed", "Stop2_mass_cm_OldSeed",
                "Stop1_pt_cm_OldSeed",   "Stop2_pt_cm_OldSeed",
                "Stop1_phi_cm_OldSeed",  "Stop2_phi_cm_OldSeed",
                "Stop1_eta_cm_OldSeed",  "Stop2_eta_cm_OldSeed",
                "Stop1_mass_cm_OldSeed_NonIsoMuon", "Stop2_mass_cm_OldSeed_NonIsoMuon",
                "Stop1_pt_cm_OldSeed_NonIsoMuon",   "Stop2_pt_cm_OldSeed_NonIsoMuon",
                "Stop1_phi_cm_OldSeed_NonIsoMuon",  "Stop2_phi_cm_OldSeed_NonIsoMuon",
                "Stop1_eta_cm_OldSeed_NonIsoMuon",  "Stop2_eta_cm_OldSeed_NonIsoMuon",
                "Stop1_mass_cm_TopSeed", "Stop2_mass_cm_TopSeed",
                "Stop1_pt_cm_TopSeed",   "Stop2_pt_cm_TopSeed",
                "Stop1_phi_cm_TopSeed",  "Stop2_phi_cm_TopSeed",
                "Stop1_eta_cm_TopSeed",  "Stop2_eta_cm_TopSeed",
                "top1_pt_cm",            "top2_pt_cm",
                "top1_eta_cm",           "top2_eta_cm",
                "top1_phi_cm",           "top2_phi_cm",
                "top1_mass_cm",          "top2_mass_cm",
            };

            if ( runtype == "MC" )
            {
                variables.insert("stop1_ptrank_mass");
                variables.insert("stop2_ptrank_mass");
                variables.insert("TotalWeight_0l");
                variables.insert("TotalWeight_1l");
                variables.insert("TotalWeight_2l");
                variables.insert("TotalWeight_QCDCR");
                variables.insert("scaleWeightUp");
                variables.insert("scaleWeightDown");
                variables.insert("PSweight_ISRUp");
                variables.insert("PSweight_ISRDown");
                variables.insert("PSweight_FSRUp");
                variables.insert("PSweight_FSRDown");
                variables.insert("PDFweightUp");
                variables.insert("PDFweightDown");
                variables.insert("jetTrigSF");
                variables.insert("jetTrigSF_Up");
                variables.insert("jetTrigSF_Down");
                variables.insert("totGoodElectronSF");
                variables.insert("totGoodElectronSF_Up");
                variables.insert("totGoodElectronSF_Down");
                variables.insert("totGoodMuonSF");
                variables.insert("totGoodMuonSF_Up");
                variables.insert("totGoodMuonSF_Down");
                variables.insert("totNonIsoMuonSF");
                variables.insert("totNonIsoMuonSF_Up");
                variables.insert("totNonIsoMuonSF_Down");
                variables.insert("topTaggerScaleFactor");
                variables.insert("topTaggerScaleFactorUp");
                variables.insert("topTaggerScaleFactorDown");
                variables.insert("prefiringScaleFactor");
                variables.insert("prefiringScaleFactorUp");
                variables.insert("prefiringScaleFactorDown");
                variables.insert("bTagSF_EventWeightSimple_Central");
                variables.insert("bTagSF_EventWeightSimple_Up");
                variables.insert("bTagSF_EventWeightSimple_Down");
                variables.insert("puWeightCorr");
                variables.insert("puSysUpCorr");
                variables.insert("puSysDownCorr");
            }   

            std::string myTreeName = "AnaSkim";
            myTree = new TTree( (myTreeName).c_str() , (myTreeName).c_str() );
            myAnaSkimTuple = new MiniTupleMaker( myTree );
            myAnaSkimTuple->setTupleVars(variables);
            myAnaSkimTuple->initBranches(tr);

            treeInit = true;
        }
        
        //-----------------------------------
        //-- Fill Histograms Below
        //-----------------------------------
        if( passBaseline0l or passBaseline1l or passBaseline2l or passQCDCR0l or passQCDCR1l ) {
            myAnaSkimTuple->fill(tr);
        }
    } 
}
      
void MakeAnaSkimTree::WriteHistos( TFile* outfile ) 
{
    outfile->cd();
    myTree->Write();
    eventCounter->Write();

    delete myTree;    
    delete myAnaSkimTuple;

}
