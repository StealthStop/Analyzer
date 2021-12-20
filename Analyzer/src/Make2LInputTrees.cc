#define Make2LInputTrees_cxx
#include "Analyzer/Analyzer/include/Make2LInputTrees.h"
#include "SusyAnaTools/Tools/NTupleReader.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TRandom3.h>
#include <iostream>

#include "SusyAnaTools/Tools/MiniTupleMaker.h"

Make2LInputTrees::Make2LInputTrees()
{
    InitHistos();
}


void Make2LInputTrees::InitHistos()
{
    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ); 
}//END of init histos

void Make2LInputTrees::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    while( tr.getNextEvent() )
    {
        const auto& eventCounter        = tr.getVar<int>("eventCounter");
        my_histos["EventCounter"]->Fill( eventCounter );

        const auto& runtype             = tr.getVar<std::string>("runtype");
//        const auto& GoodLeptons         = tr.getVec<std::pair<std::string, utility::LorentzVector>>("GoodLeptons");
        const auto& GoodLeptonsCharge   = tr.getVec<int>("GoodLeptonsCharge");
//        const auto& Jets                = tr.getVec<utility::LorentzVector>("Jets");
        const auto& NGoodLeptons        = tr.getVar<int>("NGoodLeptons");
//        const auto& onZ                 = tr.getVar<bool>("onZ");
        const auto& JetID               = tr.getVar<bool>("JetID");
        const auto& NGoodBJets_pt30     = tr.getVar<int>("NGoodBJets_pt30");
//        const auto& HT_trigger_pt30     = tr.getVar<double>("HT_trigger_pt30");
//        const auto& NGoodJets_pt30      = tr.getVar<int>("NGoodJets_pt30");

        const auto& TwoLep_Mbl1              = tr.getVar<double>("TwoLep_Mbl1");
        const auto& TwoLep_Mbl2              = tr.getVar<double>("TwoLep_Mbl2");
//        const auto& GoodBJets_pt30           = tr.getVec<bool>("GoodBJets_pt30");



        bool pass_2l = NGoodLeptons==2;
        bool pass_2l_opc =  pass_2l ? GoodLeptonsCharge[0]!=GoodLeptonsCharge[1] : false;
        bool baseline_2l =  JetID && pass_2l_opc && NGoodBJets_pt30 >= 1 && 25 < TwoLep_Mbl1 && 25 < TwoLep_Mbl2;

        //------------------------------------
        //-- Print Event Number
        //------------------------------------

        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() % 10000 == 0 ) printf( " Event %i\n", tr.getEvtNum() );

        //-----------------------------------
        //  Initialize the tree
        //-----------------------------------
        
        std::set<std::string> variables = {
//            "FirstComboCandidates",
//            "SecComboCandidates",
//            "FirstStopMassSums",
//            "SecStopMassSums",
//            "StopMassDiffs",
//            "StopMT2s",

//            "GM_Stop1",
//            "GM_Stop2",
//            "GM_Stop1Gen",
//            "GM_Stop2Gen",
//            "GM_Nlino1",
//            "GM_Nlino2",
//            "GM_Nlino1Gen",
//           "GM_Nlino2Gen",
///            "GM_Single1",
//            "GM_Single2",
//            "GM_Single1Gen",
//            "GM_Single2Gen",
//            "GM_StopMT2",
//            "GM_StopGenMT2",
//            "fracGenMatched",
            
            "GoodJetsMass",
            "GoodJetsPt",
            "GoodJetsEta",
            "GoodJetsPhi",
            "GoodJetsPx",
            "GoodJetsPy",
            "GoodJetsPz",
            "GoodJetsE",
            "GoodLeptonsMass",
            "GoodLeptonsPt",
            "GoodLeptonsEta",
            "GoodLeptonsPhi",
            "GoodLeptonsPx",
            "GoodLeptonsPy",
            "GoodLeptonsPz",
            "GoodLeptonsE",

            "NGoodJets_pt30",
            "NGoodBJets_pt30",
            "TwoLep_Mbl1",
            "TwoLep_Mbl2",
            "puWeightCorr",
            "Weight",
            "totalEventWeight",
            "NGoodBJets_pt45"
            
        };
        if( tr.isFirstEvent() ) {
            std::string myTreeName = "myMiniTree";
            myTree = new TTree( (myTreeName).c_str() , (myTreeName).c_str() );
            myMiniTuple = new MiniTupleMaker( myTree );
            myMiniTuple->setTupleVars(variables);
            myMiniTuple->initBranches(tr);
        }
        
        if( runtype == "MC" ) {
            //        const auto& passTriggerMC       = tr.getVar<bool>("passTriggerMC");
            const auto& passMadHT           = tr.getVar<bool>("passMadHT");
            
            if( !passMadHT ) continue; 
        }

        //-----------------------------------
        //-- Fill Histograms Below
        //-----------------------------------
        if( baseline_2l) {
            myMiniTuple->fill();
        }

    }//END of while tr.getNextEvent loop   
}//END of function
      
void Make2LInputTrees::WriteHistos( TFile* outfile ) 
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

    myTree->Write();
    delete myTree;    
    delete myMiniTuple;

}

