#define ResolvedTopTagger_Analyzer_cxx
#include "Analyzer/Analyzer/include/ResolvedTopTagger_Analyzer.h"
#include "Framework/Framework/include/SetUpTopTagger.h"
#include "Framework/Framework/include/Utility.h"
#include "NTupleReader/include/NTupleReader.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TRandom3.h>
#include <iostream>

ResolvedTopTagger_Analyzer::ResolvedTopTagger_Analyzer() : hists_old("histos_old"), histNjet6_old("Njet6_old"), histNjet7_old("Njet7_old"), histNjet8_old("Njet8_old"), histNjet9_old("Njet9_old"), 
                                       histNjet10_old("Njet10_old"), histNjet11_old("Njet11_old"), histNjet12_old("Njet12_old"), histNjet12inc_old("Njet12inc_old"), 
                                       hists_new("histos_new"), histNjet7_new("Njet7_new"), histNjet8_new("Njet8_new"), histNjet9_new("Njet9_new"), 
                                       histNjet10_new("Njet10_new"), histNjet11_new("Njet11_new"), histNjet12_new("Njet12_new"), histNjet12inc_new("Njet12inc_new")
{
    InitHistos();
}

void ResolvedTopTagger_Analyzer::InitHistos()
{
    TH1::SetDefaultSumw2();
    my_histos.emplace("EventCounter", std::make_shared<TH1D>("EventCounter","EventCounter", 2, -1.1, 1.1 ) );
}

void ResolvedTopTagger_Analyzer::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    TRandom3 rand(123);
   
    while(tr.getNextEvent())
    {
        const auto& eventCounter           = tr.getVar<int>("eventCounter");
        my_histos["EventCounter"]->Fill(eventCounter);

        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() & (10000 == 0) ) printf( " Event %i\n", tr.getEvtNum() );        

        const auto& runtype          = tr.getVar<std::string>("runtype");
        const auto& Jets             = tr.getVec<utility::LorentzVector>("Jets");
        // old resolved selection based on old baseline
        const auto& JetID            = tr.getVar<bool>("JetID");        
        const auto& passMETFilters   = tr.getVar<bool>("passMETFilters");
        const auto& passMadHT        = tr.getVar<bool>("passMadHT");
        const auto& passTriggerHadMC = tr.getVar<bool>("passTriggerHadMC");
        const auto& NGoodLeptons     = tr.getVar<int>("NGoodLeptons");
        const auto& HT_trigger_pt45  = tr.getVar<double>("HT_trigger_pt45");
        const auto& NGoodBJets_pt45  = tr.getVar<int>("NGoodBJets_pt45");
        const auto& NGoodJets_pt45   = tr.getVar<int>("NGoodJets_pt45");
        const auto& dR_bjets         = tr.getVar<float>("dR_bjets");
        const auto& GoodJets_pt45    = tr.getVec<bool>("GoodJets_pt45");
        const bool pass_oldResolved  = JetID && passMETFilters && passMadHT && passTriggerHadMC 
                                      && NGoodLeptons==0       && HT_trigger_pt45 > 500 
                                      && NGoodBJets_pt45 >= 2  && NGoodJets_pt45 >= 6
                                      && dR_bjets >= 1.0;
        // new resolved selection based on new baseline
        const bool passBaseline0l_pre = tr.getVar<bool>("passBaseline0l_pre");
        const auto& NNonIsoMuons      = tr.getVar<int>("NNonIsoMuons");
        const auto& NGoodJets_pt30    = tr.getVar<int>("NGoodJets_pt30");
        const auto& NGoodBJets_pt30   = tr.getVar<int>("NGoodBJets_pt30");
        const bool pass_newResolved   = passBaseline0l_pre     
                                       && NNonIsoMuons == 0   && NGoodBJets_pt30 >= 2 
                                       && NGoodJets_pt30 >= 7 && dR_bjets >= 1.0;

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
            const auto& Weight = tr.getVar<float>("Weight");
            const auto& lumi   = tr.getVar<double>("Lumi");
            eventweight        = lumi*Weight;

            bTagScaleFactor      = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor");
            puScaleFactor        = tr.getVar<double>("puWeightCorr");

            //weight *= eventweight*bTagScaleFactor*prefiringScaleFactor*puScaleFactor;
            weight *= eventweight * puScaleFactor;
        }

        // ---------------------------------------------------------------
        // -- Create variables for setStealthStopVar() in HistoContainer.h
        // ---------------------------------------------------------------
        auto& goodjets_pt45 = tr.createDerivedVec<utility::LorentzVector>("GoodJets_pt45_tlv");
        
        for (unsigned int i = 0; i < Jets.size(); ++i )
        {
            if (!GoodJets_pt45[i]) continue;
            goodjets_pt45.emplace_back(Jets.at(i));            
        }

        // -----------------------------------------
        // -- Fill the histograms 
        // -----------------------------------------
        // ----------------
        // -- baseline cuts
        // ----------------
        std::vector<std::pair<std::string, bool>> old_resolved =
        {   
            {"pass_oldResolved", pass_oldResolved},
        };
        hists_old.fillWithCutFlow(old_resolved, tr, weight, &rand);

        std::vector<std::pair<std::string, bool>> new_resolved =
        {
            {"pass_newResolved", pass_newResolved},
        };
        hists_new.fillWithCutFlow(new_resolved, tr, weight, &rand);

        // -------------------------------------------------
        // for MVA score distribution as a function of Njets
        // -------------------------------------------------
        // --------------------------
        // baseline cuts + Njets == 6
        // --------------------------        
        std::vector<std::pair<std::string, bool>> Njets6_old =
        {
            {"pass_oldResolved", pass_oldResolved   },
            {"Njet6_old"       , NGoodJets_pt45 == 6},
        };
        histNjet6_old.fillWithCutFlow(Njets6_old, tr, weight, &rand);

        // --------------------------
        // baseline cuts + Njets == 7
        // --------------------------
        std::vector<std::pair<std::string, bool>> Njets7_old =
        {
            {"pass_oldResolved", pass_oldResolved   },
            {"Njet7_old"       , NGoodJets_pt45 == 7},
        };
        histNjet7_old.fillWithCutFlow(Njets7_old, tr, weight, &rand);

        std::vector<std::pair<std::string, bool>> Njets7_new =
        {
            {"pass_newResolved", pass_newResolved   },
            {"Njet7_new"       , NGoodJets_pt30 == 7},
        };
        histNjet7_new.fillWithCutFlow(Njets7_new, tr, weight, &rand);

        // --------------------------
        // baseline cuts + Njets == 8
        // --------------------------
        std::vector<std::pair<std::string, bool>> Njets8_old =
        {
            {"pass_oldResolved", pass_oldResolved   },
            {"Njet8_old"       , NGoodJets_pt45 == 8},
        };
        histNjet8_old.fillWithCutFlow(Njets8_old, tr, weight, &rand);

        std::vector<std::pair<std::string, bool>> Njets8_new =
        {
            {"pass_newResolved", pass_newResolved   },
            {"Njet8_new"       , NGoodJets_pt30 == 8},
        };
        histNjet8_new.fillWithCutFlow(Njets8_new, tr, weight, &rand);

        // --------------------------
        // baseline cuts + Njets == 9
        // --------------------------
        std::vector<std::pair<std::string, bool>> Njets9_old =
        {
            {"pass_oldResolved", pass_oldResolved   },
            {"Njet9_old"    ,    NGoodJets_pt45 == 9},
        };
        histNjet9_old.fillWithCutFlow(Njets9_old, tr, weight, &rand);        
       
        std::vector<std::pair<std::string, bool>> Njets9_new =
        {
            {"pass_newResolved", pass_newResolved   },
            {"Njet9_new"       , NGoodJets_pt30 == 9},
        };
        histNjet9_new.fillWithCutFlow(Njets9_new, tr, weight, &rand);
 
        // ---------------------------
        // baseline cuts + Njets == 10
        // ---------------------------
        std::vector<std::pair<std::string, bool>> Njets10_old =
        {   
            {"pass_oldResolved", pass_oldResolved    },
            {"Njet10_old"      , NGoodJets_pt45 == 10},
        };
        histNjet10_old.fillWithCutFlow(Njets10_old, tr, weight, &rand);

        std::vector<std::pair<std::string, bool>> Njets10_new =
        {
            {"pass_newResolved" , pass_newResolved    },
            {"Njet10_new"       , NGoodJets_pt30 == 10},
        };
        histNjet10_new.fillWithCutFlow(Njets10_new, tr, weight, &rand);

        // ---------------------------
        // baseline cuts + Njets == 11
        // ---------------------------
        std::vector<std::pair<std::string, bool>> Njets11_old =
        {
            {"pass_oldResolved", pass_oldResolved    },
            {"Njet11_old"      , NGoodJets_pt45 == 11},
        };
        histNjet11_old.fillWithCutFlow(Njets11_old, tr, weight, &rand);

        std::vector<std::pair<std::string, bool>> Njets11_new =
        {
            {"pass_newResolved" , pass_newResolved    },
            {"Njet11_new"       , NGoodJets_pt30 == 11},
        };
        histNjet11_new.fillWithCutFlow(Njets11_new, tr, weight, &rand);

        // ---------------------------
        // baseline cuts + Njets == 12
        // ---------------------------
        std::vector<std::pair<std::string, bool>> Njets12_old =
        {
            {"pass_oldResolved", pass_oldResolved    },
            {"Njet12_old"      , NGoodJets_pt45 == 12},
        };
        histNjet12_old.fillWithCutFlow(Njets12_old, tr, weight, &rand);

        std::vector<std::pair<std::string, bool>> Njets12_new =
        {
            {"pass_newResolved" , pass_newResolved    },
            {"Njet12_new"       , NGoodJets_pt30 == 12},
        };
        histNjet12_new.fillWithCutFlow(Njets12_new, tr, weight, &rand);        

        // ---------------------------
        // baseline cuts + Njets >= 12
        // ---------------------------
        std::vector<std::pair<std::string, bool>> Njets12inc_old =
        {
            {"pass_oldResolved", pass_oldResolved    },
            {"Njet12inc_old"   , NGoodJets_pt45 >= 12},
        };
        histNjet12inc_old.fillWithCutFlow(Njets12inc_old, tr, weight, &rand);

        std::vector<std::pair<std::string, bool>> Njets12inc_new =
        {
            {"pass_newResolved" , pass_newResolved    },
            {"Njet12inc_new"    , NGoodJets_pt30 >= 12},
        };
        histNjet12inc_new.fillWithCutFlow(Njets12inc_new, tr, weight, &rand);

    }
}

void ResolvedTopTagger_Analyzer::WriteHistos(TFile* outfile)
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
   
    hists_old.save(outfile);
    histNjet6_old.save(outfile);
    histNjet7_old.save(outfile);
    histNjet8_old.save(outfile);
    histNjet9_old.save(outfile);
    histNjet10_old.save(outfile);
    histNjet11_old.save(outfile);
    histNjet12_old.save(outfile);
    histNjet12inc_old.save(outfile);

    hists_new.save(outfile);
    histNjet7_new.save(outfile);
    histNjet8_new.save(outfile);
    histNjet9_new.save(outfile);
    histNjet10_new.save(outfile);
    histNjet11_new.save(outfile);
    histNjet12_new.save(outfile);
    histNjet12inc_new.save(outfile); 
    //outfile->Write();
    outfile->Close();
}
