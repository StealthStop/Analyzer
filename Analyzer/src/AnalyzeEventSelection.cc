#define AnalyzeEventSelection_cxx
#include "Analyzer/Analyzer/include/AnalyzeEventSelection.h"
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

AnalyzeEventSelection::AnalyzeEventSelection()
{
    InitHistos();
}

void AnalyzeEventSelection::InitHistos()
{
    TH1::SetDefaultSumw2();

    // Cut flows
    my_efficiencies.emplace("event_sel", std::make_shared<TEfficiency>("event_sel","event_sel",8,0,8));
    my_efficiencies.emplace("event_sel_weight", std::make_shared<TEfficiency>("event_sel_weight","event_sel_weight",9,0,9));
}

void AnalyzeEventSelection::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    while( tr.getNextEvent() )
    {
        const auto& runtype              = tr.getVar<std::string>("runtype");     
        const auto& GoodLeptons          = tr.getVec<std::pair<std::string, utility::LorentzVector>>("GoodLeptons");

        const auto& JetID                = tr.getVar<bool>("JetID");
        const auto& NGoodLeptons         = tr.getVar<int>("NGoodLeptons");
        const auto& passTriggerMC        = tr.getVar<bool>("passTriggerMC");
        const auto& NGoodBJets_pt30      = tr.getVar<int>("NGoodBJets_pt30");
        const auto& Mbl                  = tr.getVar<double>("Mbl");
        const auto& HT_trigger_pt30      = tr.getVar<double>("HT_trigger_pt30");
        const auto& NGoodJets_pt30       = tr.getVar<int>("NGoodJets_pt30");
       
        if(maxevents != -1 && tr.getEvtNum() >= maxevents) break;        
        if ( tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() ) ;

        // ------------------------
        // -- Define weight
        // ------------------------
        double weight = 1.0;
        double eventweight = 1.0;
        double leptonweight = 1.0;
        double bTagWeight = 1.0;
        double htDerivedweight = 1.0;
        double prefiringScaleFactor = 1.0;
        if(runtype == "MC")
        {
            // Define Lumi weight
            const auto& Weight  = tr.getVar<double>("Weight");
            const auto& lumi = tr.getVar<double>("Lumi");
            eventweight = lumi*Weight;
            
            // Define lepton weight
            if(NGoodLeptons == 1)
            {
                const auto& eleLepWeight = tr.getVar<double>("totGoodElectronSF");
                const auto& muLepWeight  = tr.getVar<double>("totGoodMuonSF");
                leptonweight = (GoodLeptons[0].first == "e") ? eleLepWeight : muLepWeight;
            }
            
            //PileupWeight = tr.getVar<double>("_PUweightFactor");
            bTagWeight   = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            htDerivedweight = tr.getVar<double>("htDerivedweight");
            prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor");
            
            weight *= eventweight*leptonweight*bTagWeight*htDerivedweight*prefiringScaleFactor;
        }
        
        // Fill event selection efficiencies
        my_efficiencies["event_sel"]->Fill(true,0);
        my_efficiencies["event_sel"]->Fill(true && JetID,1);
        my_efficiencies["event_sel"]->Fill(true && JetID && NGoodLeptons == 1,2);
        my_efficiencies["event_sel"]->Fill(true && JetID && NGoodLeptons == 1 && passTriggerMC,3);
        my_efficiencies["event_sel"]->Fill(true && JetID && NGoodLeptons == 1 && passTriggerMC && NGoodBJets_pt30 >= 1,4);
        my_efficiencies["event_sel"]->Fill(true && JetID && NGoodLeptons == 1 && passTriggerMC && NGoodBJets_pt30 >= 1 && 50 < Mbl && Mbl < 250,5);
        my_efficiencies["event_sel"]->Fill(true && JetID && NGoodLeptons == 1 && passTriggerMC && NGoodBJets_pt30 >= 1 && 50 < Mbl && Mbl < 250 && HT_trigger_pt30 > 300,6);
        my_efficiencies["event_sel"]->Fill(true && JetID && NGoodLeptons == 1 && passTriggerMC && NGoodBJets_pt30 >= 1 && 50 < Mbl && Mbl < 250 && HT_trigger_pt30 > 300 && NGoodJets_pt30 >= 7,7);

        my_efficiencies["event_sel_weight"]->FillWeighted(true,eventweight,0);
        my_efficiencies["event_sel_weight"]->FillWeighted(true && JetID,eventweight,1);
        my_efficiencies["event_sel_weight"]->FillWeighted(true && JetID && NGoodLeptons == 1,eventweight,2);
        my_efficiencies["event_sel_weight"]->FillWeighted(true && JetID && NGoodLeptons == 1 && passTriggerMC,eventweight,3);
        my_efficiencies["event_sel_weight"]->FillWeighted(true && JetID && NGoodLeptons == 1 && passTriggerMC && NGoodBJets_pt30 >= 1,eventweight,4);
        my_efficiencies["event_sel_weight"]->FillWeighted(true && JetID && NGoodLeptons == 1 && passTriggerMC && NGoodBJets_pt30 >= 1 && 50 < Mbl && Mbl < 250,eventweight,5);
        my_efficiencies["event_sel_weight"]->FillWeighted(true && JetID && NGoodLeptons == 1 && passTriggerMC && NGoodBJets_pt30 >= 1 && 50 < Mbl && Mbl < 250 && HT_trigger_pt30 > 300,eventweight,6);
        my_efficiencies["event_sel_weight"]->FillWeighted(true && JetID && NGoodLeptons == 1 && passTriggerMC && NGoodBJets_pt30 >= 1 && 50 < Mbl && Mbl < 250 && HT_trigger_pt30 > 300 && NGoodJets_pt30 >= 7,eventweight,7);
        my_efficiencies["event_sel_weight"]->FillWeighted(true && JetID && NGoodLeptons == 1 && passTriggerMC && NGoodBJets_pt30 >= 1 && 50 < Mbl && Mbl < 250 && HT_trigger_pt30 > 300 && NGoodJets_pt30 >= 7,weight,8);
    } 
}

void AnalyzeEventSelection::WriteHistos(TFile* outfile)
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
