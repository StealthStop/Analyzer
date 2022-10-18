#define ResolvedTopTagger_Analyzer_cxx
#include "Analyzer/Analyzer/include/ResolvedTopTagger_Analyzer.h"
#include "Framework/Framework/include/SetUpTopTagger.h"
#include "Framework/Framework/include/Utility.h"
#include "NTupleReader/include/NTupleReader.h"

#include <TH1D.h>
#include <TH2D.h>
#include <iostream>

ResolvedTopTagger_Analyzer::ResolvedTopTagger_Analyzer() : hists("histos"), histNjet7("Njet7"), histNjet8("Njet8"), histNjet9("Njet9"), 
                                                           histNjet10("Njet10"), histNjet11("Njet11"), histNjet12("Njet12"), histNjet12inc("Njet12inc")
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
        const auto& JetID               = tr.getVar<bool>("JetID");        
        const auto& passMETFilters      = tr.getVar<bool>("passMETFilters");
        const auto& passMadHT           = tr.getVar<bool>("passMadHT");
        const auto& passTriggerHadMC    = tr.getVar<bool>("passTriggerHadMC");
        const auto& NGoodLeptons        = tr.getVar<int>("NGoodLeptons");
        const auto& HT_trigger_pt45     = tr.getVar<double>("HT_trigger_pt45");
        const auto& NGoodBJets_pt45     = tr.getVar<int>("NGoodBJets_pt45");
        const auto& NGoodJets_pt45      = tr.getVar<int>("NGoodJets_pt45");
        const auto& dR_bjets_old        = tr.getVar<double>("dR_bjets_old");
        const auto& dR_bjets            = tr.getVar<double>("dR_bjets");
        const auto& GoodJets_pt45       = tr.getVec<bool>("GoodJets_pt45");
        const auto& passElectronHEMveto = tr.getVar<bool>("passElectronHEMveto");
        const bool pass_oldResolved     = JetID && passMETFilters && passMadHT && passTriggerHadMC 
                                         && NGoodLeptons==0       && HT_trigger_pt45 > 500 
                                         && NGoodBJets_pt45 >= 2  && NGoodJets_pt45 >= 6
                                         && dR_bjets_old >= 1.0;
        // new resolved selection based on new baseline
        const bool passBaseline0l_pre = tr.getVar<bool>("passBaseline0l_pre");
        const auto& NNonIsoMuons      = tr.getVar<int>("NNonIsoMuons");
        const auto& NGoodJets_pt30    = tr.getVar<int>("NGoodJets_pt30");
        const auto& NGoodBJets_pt30   = tr.getVar<int>("NGoodBJets_pt30");
        const bool pass_newResolved   = passBaseline0l_pre     
                                       && NNonIsoMuons == 0   && NGoodBJets_pt30 >= 2 
                                       && NGoodJets_pt30 >= 7 && dR_bjets >= 1.0
                                       && passElectronHEMveto;

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
            const auto& lumi   = tr.getVar<double>("FinalLumi");
            eventweight        = lumi*Weight;

            bTagScaleFactor      = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor");
            puScaleFactor        = tr.getVar<double>("puWeightCorr");

            weight *= eventweight * bTagScaleFactor * prefiringScaleFactor * puScaleFactor;
        }

        // ---------------------------------------------------------------
        // -- Create variables for setStealthStopVar() in HistoContainer.h
        // ---------------------------------------------------------------
        auto& goodjets_pt45 = tr.createDerivedVec<TLorentzVector>("GoodJets_pt45_tlv");
        
        for (unsigned int i = 0; i < Jets.size(); ++i )
        {
            if (!GoodJets_pt45[i]) continue;
            //goodjets_pt45.emplace_back(Jets.at(i));            
            goodjets_pt45.emplace_back(utility::convertLV<TLorentzVector, utility::LorentzVector>(Jets.at(i)));
        }

        // -----------------------------------------
        // -- Fill the histograms 
        // -----------------------------------------
        // ----------------
        // -- baseline cuts
        // ----------------
        std::vector<std::pair<std::string, bool>> new_resolved =
        {
            {"pass_newResolved", pass_newResolved},
        };
        hists.fillWithCutFlow(new_resolved, tr, weight, &rand);

        // -------------------------------------------------
        // for MVA score distribution as a function of Njets
        // -------------------------------------------------
        // --------------------------
        // baseline cuts + Njets == 7
        // --------------------------
        std::vector<std::pair<std::string, bool>> Njets7 =
        {
            {"pass_newResolved", pass_newResolved   },
            {"Njet7"           , NGoodJets_pt30 == 7},
        };
        histNjet7.fillWithCutFlow(Njets7, tr, weight, &rand);

        // --------------------------
        // baseline cuts + Njets == 8
        // --------------------------
        std::vector<std::pair<std::string, bool>> Njets8 =
        {
            {"pass_newResolved", pass_newResolved   },
            {"Njet8"           , NGoodJets_pt30 == 8},
        };
        histNjet8.fillWithCutFlow(Njets8, tr, weight, &rand);

        // --------------------------
        // baseline cuts + Njets == 9
        // --------------------------
        std::vector<std::pair<std::string, bool>> Njets9 =
        {
            {"pass_newResolved", pass_newResolved   },
            {"Njet9"           , NGoodJets_pt30 == 9},
        };
        histNjet9.fillWithCutFlow(Njets9, tr, weight, &rand);
 
        // ---------------------------
        // baseline cuts + Njets == 10
        // ---------------------------
        std::vector<std::pair<std::string, bool>> Njets10 =
        {
            {"pass_newResolved" , pass_newResolved    },
            {"Njet10"           , NGoodJets_pt30 == 10},
        };
        histNjet10.fillWithCutFlow(Njets10, tr, weight, &rand);

        // ---------------------------
        // baseline cuts + Njets == 11
        // ---------------------------
        std::vector<std::pair<std::string, bool>> Njets11 =
        {
            {"pass_newResolved" , pass_newResolved    },
            {"Njet11"           , NGoodJets_pt30 == 11},
        };
        histNjet11.fillWithCutFlow(Njets11, tr, weight, &rand);

        // ---------------------------
        // baseline cuts + Njets == 12
        // ---------------------------
        std::vector<std::pair<std::string, bool>> Njets12 =
        {
            {"pass_newResolved" , pass_newResolved    },
            {"Njet12"           , NGoodJets_pt30 == 12},
        };
        histNjet12.fillWithCutFlow(Njets12, tr, weight, &rand);        

        // ---------------------------
        // baseline cuts + Njets >= 12
        // ---------------------------

        std::vector<std::pair<std::string, bool>> Njets12inc =
        {
            {"pass_newResolved" , pass_newResolved    },
            {"Njet12inc"        , NGoodJets_pt30 >= 12},
        };
        histNjet12inc.fillWithCutFlow(Njets12inc, tr, weight, &rand);

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
    
    hists.save(outfile);
    histNjet7.save(outfile);
    histNjet8.save(outfile);
    histNjet9.save(outfile);
    histNjet10.save(outfile);
    histNjet11.save(outfile);
    histNjet12.save(outfile);
    histNjet12inc.save(outfile); 
    //outfile->Write();
    outfile->Close();
}
