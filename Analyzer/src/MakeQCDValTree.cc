#define MakeQCDValTree_cxx
#include "Analyzer/Analyzer/include/MakeQCDValTree.h"
#include "NTupleReader/include/NTupleReader.h"
#include "Framework/Framework/include/MiniTupleMaker.h"

#include <iostream>

MakeQCDValTree::MakeQCDValTree()
{
    InitHistos();

    treeInit = false;
}

void MakeQCDValTree::InitHistos()
{
    eventCounter = std::make_shared<TH1D>("EventCounter", "EventCounter", 2, -1.1, 1.1);
}

void MakeQCDValTree::Loop(NTupleReader& tr, double, int maxevents, bool)
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

        //if (lostCauseEvent and fastMode){
        //    continue;
        //}

        const auto& runtype         = tr.getVar<std::string>("runtype");
        const auto& filetag         = tr.getVar<std::string>("filetag");

        // General requirements to always pass
        const auto& JetID               = tr.getVar<bool>("JetID");
        const auto& passMadHT           = tr.getVar<bool>("passMadHT");
        const auto& passMETFilters      = tr.getVar<bool>("passMETFilters");
        const auto& passNonIsoTrigger   = tr.getVar<bool>("passNonIsoTrigger");
        const auto& passNonIsoTriggerMC = tr.getVar<bool>("passNonIsoTriggerMC");
        const auto& passElectronHEMveto = tr.getVar<bool>("passElectronHEMveto");

        bool passGeneral = JetID               && 
                           passMadHT           &&  
                           passMETFilters      &&  
                           passNonIsoTrigger   &&  
                           passNonIsoTriggerMC &&
                           passElectronHEMveto &&
                           (runtype != "Data" || filetag.find("Data_SingleMuon") != std::string::npos);


        // Relevant lepton quantities
        const auto& NGoodMuons             = tr.getVar<int>("NGoodMuons");
        const auto& NGoodElectrons         = tr.getVar<int>("NGoodElectrons");
        const auto& NNonIsoMuons           = tr.getVar<int>("NNonIsoMuons");

        bool passLeptonReqs = NNonIsoMuons == 1 &&
                              NGoodMuons == 0   &&
                              NGoodElectrons == 0;

        const auto& NGoodJets_pt30       = tr.getVar<int>("NGoodJets_pt30");

        if( !treeInit ) {
            //-----------------------------------
            //  Initialize the tree
            //-----------------------------------       
            std::set<std::string> variables = {
                "NGoodJets_pt30",
                "NGoodBJets_pt30",
                "NNonIsoMuonJets_pt30",
                "NGoodJets_pt45",
                "NGoodBJets_pt45",
                "HT_trigger_pt30",
                "HT_NonIsoMuon_pt30",
                "Mbb",
                "dR_bjets",
                "nimMbb",
                "dR_nimbjets",
                "DoubleDisCo_disc1_NonIsoMuon_0l_RPV",
                "DoubleDisCo_disc2_NonIsoMuon_0l_RPV",
                "DoubleDisCo_disc1_NonIsoMuon_0l_SYY",
                "DoubleDisCo_disc2_NonIsoMuon_0l_SYY",
                "DoubleDisCo_disc1_NonIsoMuon_1l_RPV",
                "DoubleDisCo_disc2_NonIsoMuon_1l_RPV",
                "DoubleDisCo_disc1_NonIsoMuon_1l_SYY",
                "DoubleDisCo_disc2_NonIsoMuon_1l_SYY",
                "DoubleDisCo_disc1_0l_RPV",
                "DoubleDisCo_disc2_0l_RPV",
                "DoubleDisCo_disc1_0l_SYY",
                "DoubleDisCo_disc2_0l_SYY",
                "DoubleDisCo_disc1_1l_RPV",
                "DoubleDisCo_disc2_1l_RPV",
                "DoubleDisCo_disc1_1l_SYY",
                "DoubleDisCo_disc2_1l_SYY",
            };

            if ( runtype == "MC" )
            {
                variables.insert("TotalWeight_QCDCR");
            }   
            else
            {
                variables.insert("Weight");
            }

            std::string myTreeName = "PreSelection";
            myTree = new TTree( (myTreeName).c_str() , (myTreeName).c_str() );
            myQCDValTuple = new MiniTupleMaker( myTree );
            myQCDValTuple->setTupleVars(variables);
            myQCDValTuple->initBranches(tr);

            treeInit = true;
        }
        
        //-----------------------------------
        //-- Fill Histograms Below
        //-----------------------------------
        // Requirements not on NonIsoMuon jets or HT as those will be more restrictive with a non iso muon present
        if( passGeneral and passLeptonReqs and NGoodJets_pt30 >= 7 ) {
            myQCDValTuple->fill(tr);
        }
    } 
}
      
void MakeQCDValTree::WriteHistos( TFile* outfile ) 
{
    outfile->cd();
    myTree->Write();
    eventCounter->Write();

    delete myTree;    
    delete myQCDValTuple;

}
