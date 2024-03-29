#define MakeMiniTree_cxx
#include "Analyzer/Analyzer/include/MakeMiniTree.h"
#include "NTupleReader/include/NTupleReader.h"
#include "Framework/Framework/include/MiniTupleMaker.h"

#include <iostream>

MakeMiniTree::MakeMiniTree()
{
    InitHistos();
}

void MakeMiniTree::InitHistos()
{
    eventCounter = std::make_shared<TH1D>("EventCounter", "EventCounter", 2, -1.1, 1.1);
}

void MakeMiniTree::Loop(NTupleReader& tr, double, int maxevents, bool)
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

        const auto& runtype         = tr.getVar<std::string>("runtype");

        const auto& NGoodJets_pt30  = tr.getVar<int>("NGoodJets_pt30");
        const auto& HT_trigger_pt30 = tr.getVar<double>("HT_trigger_pt30");

        //-----------------------------------
        //  Initialize the tree
        //-----------------------------------       
        std::set<std::string> variables = {
            // Event-level
            "RunNum",
            "TriggerPass",
            "fixedGridRhoFastjetAll",
            "NVtx",
    
            // Electrons
            "Electrons",
            "Electrons_passIso",
            "Electrons_iso",
            "Electrons_charge",
            "Electrons_tightID",

            // Muons
            "Muons",
            "Muons_passIso",
            "Muons_iso",
            "Muons_charge",
            "Muons_mediumID",

            // Photons
            "Photons",
            "Photons_fullID",
    
            // MET
            "MET",
            "METPhi",

            // Jets
            "Jets",
            "Jets_pileupId",
            "Jets_bDiscriminatorCSV" ,
            "Jets_bJetTagDeepCSVprobb",
            "Jets_bJetTagDeepCSVprobbb",
            "Jets_bJetTagDeepFlavourprobb",
            "Jets_bJetTagDeepFlavourprobbb",
            "Jets_bJetTagDeepFlavourproblepb",
            "Jets_bJetTagDeepFlavourprobc",
            "Jets_bJetTagDeepFlavourprobuds",
            "Jets_bJetTagDeepFlavourprobg",
            "Jets_bJetTagDeepCSVprobc",
            "Jets_bJetTagDeepCSVprobudsg",
            "Jets_qgLikelihood",
            "Jets_muonEnergyFraction",
            "Jets_hfHadronEnergyFraction",
            "Jets_hfEMEnergyFraction",
            "Jets_photonEnergyFraction",
            "Jets_electronEnergyFraction",
            "Jets_chargedHadronMultiplicity",
            "Jets_neutralHadronMultiplicity",
            "Jets_photonMultiplicity",
            "Jets_electronMultiplicity",
            "Jets_muonMultiplicity",
            "Jets_ID",
            "JetID",
            "Jets_partonFlavor",
            "Jets_ptD",
            "Jets_axismajor",
            "Jets_axisminor",
            "Jets_multiplicity",
            "Jets_neutralEmEnergyFraction",
            "Jets_chargedEmEnergyFraction",
            "Jets_neutralHadronEnergyFraction",
            "Jets_chargedHadronEnergyFraction",

            // JetsAK8
            "JetsAK8",
            "JetsAK8_NsubjettinessTau1",
            "JetsAK8_NsubjettinessTau2",
            "JetsAK8_NsubjettinessTau3",
            "JetsAK8_softDropMass",
            "JetsAK8_axismajor",
            "JetsAK8_axisminor",
            "JetsAK8_subjets",
            "JetsAK8_subjetsCounts",
            "JetsAK8_DeepTagTvsQCD",
            "JetsAK8_DeepTagWvsQCD",
            "JetsAK8_DeepTagHbbvsQCD",
            "JetsAK8_multiplicity",

            // Event filters
            "globalSuperTightHalo2016Filter",
            "PrimaryVertexFilter",
            "BadPFMuonFilter",
            "EcalDeadCellTriggerPrimitiveFilter",
            "HBHEIsoNoiseFilter",
            "HBHENoiseFilter",
        };

        if( runtype == "MC" )
        {
            variables.insert("Weight");
            variables.insert("madHT");
            variables.insert("GenElectrons");
            variables.insert("GenMuons");
            variables.insert("GenTaus");
            variables.insert("GenMET");
            variables.insert("GenMETPhi");
            variables.insert("GenParticles");
            variables.insert("GenParticles_PdgId");
            variables.insert("GenParticles_ParentId");
            variables.insert("GenParticles_ParentIdx");
            variables.insert("GenParticles_Status");
            variables.insert("ScaleWeights");
            variables.insert("PSweights");
            variables.insert("PDFweights");
            variables.insert("puWeight");
            variables.insert("puSysUp");
            variables.insert("puSysDown");
            variables.insert("NonPrefiringProb");
            variables.insert("NonPrefiringProbUp");
            variables.insert("NonPrefiringProbDown");
            variables.insert("Jets_origIndex");
            variables.insert("JetsJECup_origIndex");
            variables.insert("JetsJECdown_origIndex");
            variables.insert("JetsJERup_origIndex");
            variables.insert("JetsJERdown_origIndex");
            variables.insert("Jets_jerFactor");
            variables.insert("Jets_jecUnc");
            variables.insert("JetsJECup_jerFactor");
            variables.insert("JetsJECdown_jerFactor");
            variables.insert("Jets_jerFactorUp");
            variables.insert("Jets_jerFactorDown");
            variables.insert("JetIDJECdown");
            variables.insert("JetIDJECup");
            variables.insert("JetIDJERdown");
            variables.insert("JetIDJERup");
            variables.insert("JetsAK8_origIndex");
            variables.insert("JetsAK8JECup_origIndex");
            variables.insert("JetsAK8JECdown_origIndex");
            variables.insert("JetsAK8JERup_origIndex");
            variables.insert("JetsAK8JERdown_origIndex");
            variables.insert("JetsAK8_jerFactor");
            variables.insert("JetsAK8_jecUnc");
            variables.insert("JetsAK8JECup_jerFactor");
            variables.insert("JetsAK8JECdown_jerFactor");
            variables.insert("JetsAK8_jerFactorUp");
            variables.insert("JetsAK8_jerFactorDown");
        }

        if( tr.isFirstEvent() ) {
            std::string myTreeName = "PreSelection";
            myTree = new TTree( (myTreeName).c_str() , (myTreeName).c_str() );
            myMiniTuple = new MiniTupleMaker( myTree );
            myMiniTuple->setTupleVars(variables);
            myMiniTuple->initBranches(tr);
        }
        
        //-----------------------------------
        //-- Fill Histograms Below
        //-----------------------------------
        if( NGoodJets_pt30 >= 5 and HT_trigger_pt30 > 400.0 ) {
            myMiniTuple->fill(tr);
        }
    } 
}
      
void MakeMiniTree::WriteHistos( TFile* outfile ) 
{
    outfile->cd();
    myTree->Write();
    eventCounter->Write();

    delete myTree;    
    delete myMiniTuple;

}
