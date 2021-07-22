#define AnalyzeQCDCR_cxx
#include "Analyzer/Analyzer/include/AnalyzeQCDCR.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "Framework/Framework/include/Utility.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TProfile2D.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TRandom3.h>
#include <iostream>

AnalyzeQCDCR::AnalyzeQCDCR() : initHistos(false)
{
}

void AnalyzeQCDCR::InitHistos(const std::map<std::string, bool>& cutMap, const std::vector<TH1DInfo>& histInfos, 
                             const std::vector<TH2DInfo>& hist2DInfos)
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    my_histos.emplace("EventCounter", std::make_shared<TH1D>("EventCounter","EventCounter", 2, -1.1, 1.1 ) );

    for(auto& mycut : cutMap)
    {
        for(const auto& hInfo : histInfos)
        { 
            my_histos.emplace(hInfo.name+mycut.first, 
                              std::make_shared<TH1D>((hInfo.name+mycut.first).c_str(),(hInfo.name+mycut.first).c_str(), hInfo.nBins, hInfo.low, hInfo.high));
        }

        for(const auto& h2dInfo : hist2DInfos)
        {
            my_2d_histos.emplace(h2dInfo.name+mycut.first, 
                                 std::make_shared<TH2D>((h2dInfo.name+mycut.first).c_str(),(h2dInfo.name+mycut.first).c_str(), 
                                                        h2dInfo.nBinsX, h2dInfo.lowX, h2dInfo.highX, h2dInfo.nBinsY, h2dInfo.lowY, h2dInfo.highY));
        }
    }

    //my_histos.emplace( "h_cutFlow", std::make_shared<TH1D>("h_cutFlow", "h_cutFlow", 9,0,9));    
}

void AnalyzeQCDCR::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    while( tr.getNextEvent() )
    {

        const auto& runtype                           = tr.getVar<std::string>("runtype");     
        const auto& Jets                              = tr.getVec<TLorentzVector>("Jets");
        const auto& GoodJets_pt30                     = tr.getVec<bool>("GoodJets_pt30");
        const auto& GoodJets_pt45                     = tr.getVec<bool>("GoodJets_pt45");
        const auto& NGoodJets_pt30                    = tr.getVar<int>("NGoodJets_pt30");
        const auto& NGoodJets_pt45                    = tr.getVar<int>("NGoodJets_pt45");
        const auto& NonIsoMuonJets_pt30               = tr.getVec<bool>("NonIsoMuonJets_pt30");
        const auto& NonIsoMuonJets_pt45               = tr.getVec<bool>("NonIsoMuonJets_pt45");
        const auto& NNonIsoMuonJets_pt30              = tr.getVar<int>("NNonIsoMuonJets_pt30");
        const auto& NNonIsoMuonJets_pt45              = tr.getVar<int>("NNonIsoMuonJets_pt45");
        const auto& NGoodBJets_pt30                   = tr.getVar<int>("NGoodBJets_pt30");
        const auto& NGoodBJets_pt45                   = tr.getVar<int>("NGoodBJets_pt45");
        const auto& ntops                             = tr.getVar<int>("ntops");
        const auto& GoodLeptons                       = tr.getVec<std::pair<std::string, TLorentzVector>>("GoodLeptons");
        const auto& GoodNonIsoMuons                   = tr.getVec<std::pair<std::string, TLorentzVector>>("GoodNonIsoMuons");
        const auto& HT_trigger_pt30                   = tr.getVar<double>("HT_trigger_pt30");
        const auto& HT_trigger_pt45                   = tr.getVar<double>("HT_trigger_pt45");
        const auto& HT_NonIsoMuon_pt30                = tr.getVar<double>("HT_NonIsoMuon_pt30");
        const auto& HT_NonIsoMuon_pt45                = tr.getVar<double>("HT_NonIsoMuon_pt45");
        const auto& correct2018Split                  = tr.getVar<bool>("correct2018Split");
        const auto& passTrigger                       = tr.getVar<bool>("passTrigger");
        const auto& passTriggerMC                     = tr.getVar<bool>("passTriggerMC");
        const auto& passMETFilters                    = tr.getVar<bool>("passMETFilters");
        const auto& passMadHT                         = tr.getVar<bool>("passMadHT");
        const auto& passBaseline0l_Good               = tr.getVar<bool>("passBaseline0l_Good");
        const auto& passBaseline1l_Good               = tr.getVar<bool>("passBaseline1l_Good");
        const auto& passBaseline1l_NonIsoMuon         = tr.getVar<bool>("passBaseline1l_NonIsoMuon");
        const auto& passHEMVeto                       = tr.getVar<bool>("passHEMVeto");
        const auto& Mbl                               = tr.getVar<double>("Mbl");
        const auto& MblVec                            = tr.getVec<double>("MblVec");
        const auto& passBlind                         = tr.getVar<bool>("passBlindLep_Good");            
        const auto& eventCounter                      = tr.getVar<int>("eventCounter");
        const auto& DoubleDisCo_disc1_NonIsoMuon_0l   = tr.getVar<double>("DoubleDisCo_disc1_NonIsoMuon_0l");
        const auto& DoubleDisCo_disc2_NonIsoMuon_0l   = tr.getVar<double>("DoubleDisCo_disc2_NonIsoMuon_0l");
        const auto& DoubleDisCo_disc1_NonIsoMuon_1l   = tr.getVar<double>("DoubleDisCo_disc1_NonIsoMuon_1l");
        const auto& DoubleDisCo_disc2_NonIsoMuon_1l   = tr.getVar<double>("DoubleDisCo_disc2_NonIsoMuon_1l");

        // ------------------------
        // -- Print event number
        // ------------------------       
        if(maxevents != -1 && tr.getEvtNum() >= maxevents) break;        
        if(tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() );

        // ------------------------
        // -- Define weight
        // ------------------------
        double weight=1.0, weightNoHT=1.0, weightQCDCR=1.0, weightNoBTag=1.0;
        double eventweight=1.0, leptonweight=1.0, bTagWeight=1.0, prefiringScaleFactor=1.0, pileupWeight=1.0, htDerivedweight=1.0;
        double weightNoLepton=1.0;
        if(runtype == "MC")
        {
            // Define Lumi weight
            const auto& Weight  = tr.getVar<double>("Weight");
            const auto& lumi = tr.getVar<double>("Lumi");
            eventweight = lumi*Weight;
            
            const auto& eleLepWeight = tr.getVar<double>("totGoodElectronSF");
            const auto& muLepWeight  = tr.getVar<double>("totGoodMuonSF");
            const auto& muNonIso     = tr.getVar<double>("totNonIsoMuonSF");
            leptonweight = eleLepWeight*muLepWeight;
            
            pileupWeight = tr.getVar<double>("puWeightCorr");
            bTagWeight   = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            htDerivedweight = tr.getVar<double>("htDerivedweight");
            prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor");
            
            weightQCDCR *= eventweight*muNonIso*prefiringScaleFactor*pileupWeight;
            weightNoHT *= eventweight*leptonweight*bTagWeight*prefiringScaleFactor*pileupWeight;
            weightNoLepton *= eventweight*bTagWeight*prefiringScaleFactor*pileupWeight*htDerivedweight;
            weightNoBTag *= eventweight*leptonweight*prefiringScaleFactor*pileupWeight*htDerivedweight;
            weight *= eventweight*leptonweight*bTagWeight*prefiringScaleFactor*pileupWeight*htDerivedweight;
        }

        // -------------------------------
        // -- Define cuts
        // -------------------------------
        bool pass_general    = passTriggerMC && passTrigger && passMadHT && passBlind && passMETFilters && passHEMVeto && correct2018Split;
        bool pass_0btag_pt30 = NGoodBJets_pt30 == 0;
        bool pass_0btag_pt45 = NGoodBJets_pt45 == 0;
        bool pass_0tops      = ntops == 0;       
 
        // -------------------
        // --- Fill Histos ---
        // -------------------                        
        const std::map<std::string, bool> cut_map 
        {
            {"_test"                        , true                                                                                     },
            {"_0l"                          , pass_general && passBaseline0l_Good                                                      },
            {"_1l"                          , pass_general && passBaseline1l_Good                                                      },
            {"_passQCDCR_0tops"             , passBaseline1l_NonIsoMuon && pass_0tops                                                  },
            {"_passQCDCR_0btag"             , passBaseline1l_NonIsoMuon && pass_0btag_pt45                                             },
            {"_passQCDCR_0tops_0btag"       , passBaseline1l_NonIsoMuon && pass_0btag_pt45 && pass_0tops                               },
            {"_passQCDCR_pt30"              , passBaseline1l_NonIsoMuon                                                                },
            {"_passQCDCR_pt30_0tops"        , passBaseline1l_NonIsoMuon && pass_0tops                                                  },
            {"_passQCDCR_pt30_0btag"        , passBaseline1l_NonIsoMuon && pass_0btag_pt45                                             },
            {"_passQCDCR_pt30_0tops_0btag"  , passBaseline1l_NonIsoMuon && pass_0btag_pt45 && pass_0tops                               },
            {"_passQCDCR_pt45"              , passBaseline1l_NonIsoMuon                                                                },
            {"_passQCDCR_pt45_0tops"        , passBaseline1l_NonIsoMuon && pass_0tops                                                  },
            {"_passQCDCR_pt45_0btag"        , passBaseline1l_NonIsoMuon && pass_0btag_pt45                                             },
            {"_passQCDCR_pt45_0tops_0btag"  , passBaseline1l_NonIsoMuon && pass_0btag_pt45 && pass_0tops                               },
            {"_passQCDCR_0l"                , pass_general && passBaseline1l_NonIsoMuon && passBaseline0l_Good                                         },
            {"_passQCDCR_0l_0tops"          , pass_general && passBaseline1l_NonIsoMuon && passBaseline0l_Good && pass_0tops                           },
            {"_passQCDCR_0l_0btag"          , pass_general && passBaseline1l_NonIsoMuon && passBaseline0l_Good && pass_0btag_pt45                      },
            {"_passQCDCR_0l_0tops_0btag"    , pass_general && passBaseline1l_NonIsoMuon && passBaseline0l_Good && pass_0btag_pt45 && pass_0tops        },
            {"_passQCDCR_1l"                , pass_general && passBaseline1l_NonIsoMuon && passBaseline1l_Good                                         },
            {"_passQCDCR_1l_0tops"          , pass_general && passBaseline1l_NonIsoMuon && passBaseline1l_Good && pass_0tops                           },
            {"_passQCDCR_1l_0btag"          , pass_general && passBaseline1l_NonIsoMuon && passBaseline1l_Good && pass_0btag_pt30                      },
            {"_passQCDCR_1l_0tops_0btag"    , pass_general && passBaseline1l_NonIsoMuon && passBaseline1l_Good && pass_0btag_pt30 && pass_0tops        },
        };

        std::vector<TH1DInfo> histInfos = {
            {"h_njets",                 20,   0.0,   20.0},
            {"h_njetsQCDCR",            20,   0.0,   20.0},
            {"h_disc1",                 50,   0.0,    1.0},
            {"h_disc2",                 50,   0.0,    1.0},
            {"h_nb",                    10,   0.0,   10.0},
            {"h_ht",                   500,   0.0, 5000.0},
            {"h_htQCDCR",              500,   0.0, 5000.0},
            {"h_mbl",                  300,   0.0,  300.0},
            {"h_lPt",                  200,   0.0, 2000.0},
            {"h_lEta",                 200,  -6.0,    6.0},
            {"h_lPhi",                 200,  -4.0,    4.0},
            {"h_isomPt",               200,   0.0, 2000.0},
            {"h_isomEta",              200,  -6.0,    6.0},
            {"h_isomPhi",              200,  -4.0,    4.0},
            {"h_jPt",                  200,   0.0, 2000.0},
            {"h_jEta",                 200,  -6.0,    6.0},
            {"h_jPhi",                 200,  -4.0,    4.0},
            {"h_allMbl",               300,   0.0,  300.0},            
        };

        std::vector<TH2DInfo> hist2DInfos = {
            {"h_disc1_disc2",            50, 0.0,  1.0, 50, 0.0, 1.0},
            {"h_njetsQCDCR_disc1",       20, 0.0, 20.0, 50, 0.0, 1.0},
            {"h_njetsQCDCR_disc2",       20, 0.0, 20.0, 50, 0.0, 1.0},
        };

        // Initialize Histograms
        if(!initHistos)
        {
            InitHistos(cut_map, histInfos, hist2DInfos);
            initHistos = true;
        }

        my_histos["EventCounter"]->Fill(eventCounter);

        for(auto& kv : cut_map)
        {
            if(kv.second)
            {
                double w = weight;
                if(kv.first.find("50to110mt")    != std::string::npos || 
                   kv.first.find("htCorr")       != std::string::npos || 
                   kv.first.find("noHTWeight")   != std::string::npos) w = weightNoHT;
                if(kv.first.find("passQCDCR")    != std::string::npos) w = weightQCDCR;
                if(kv.first.find("noLepWeight")  != std::string::npos) w = weightNoLepton;
                if(kv.first.find("noBTagWeight") != std::string::npos) w = weightNoBTag;
                if(kv.first.find("_1l") != std::string::npos || kv.first.find("pt30") != std::string::npos)
                {
                    my_histos["h_njets"               +kv.first]->Fill(NGoodJets_pt30, w);
                    my_histos["h_njetsQCDCR"          +kv.first]->Fill(NNonIsoMuonJets_pt30, w);
                    my_histos["h_nb"                  +kv.first]->Fill(NGoodBJets_pt30, w);
                    my_histos["h_ht"                  +kv.first]->Fill(HT_trigger_pt30, w);
                    my_histos["h_htQCDCR"             +kv.first]->Fill(HT_NonIsoMuon_pt30, w);
                    my_histos["h_disc1"               +kv.first]->Fill(DoubleDisCo_disc1_NonIsoMuon_1l, w);
                    my_histos["h_disc2"               +kv.first]->Fill(DoubleDisCo_disc2_NonIsoMuon_1l, w);
                    if(kv.first.find("_passQCDCR") != std::string::npos)
                    {
                        my_2d_histos["h_njetsQCDCR_disc1"      +kv.first]->Fill(NNonIsoMuonJets_pt30, DoubleDisCo_disc1_NonIsoMuon_1l, w);
                        my_2d_histos["h_njetsQCDCR_disc2"      +kv.first]->Fill(NNonIsoMuonJets_pt30, DoubleDisCo_disc2_NonIsoMuon_1l, w);
                        my_2d_histos["h_disc1_disc2"           +kv.first]->Fill(DoubleDisCo_disc1_NonIsoMuon_1l, DoubleDisCo_disc2_NonIsoMuon_1l, w);
                    }
                }
                else if(kv.first.find("_0l") != std::string::npos || kv.first.find("pt45") != std::string::npos)
                {
                    my_histos["h_njets"               +kv.first]->Fill(NGoodJets_pt45, w);
                    my_histos["h_njetsQCDCR"          +kv.first]->Fill(NNonIsoMuonJets_pt45, w);
                    my_histos["h_nb"                  +kv.first]->Fill(NGoodBJets_pt45, w);
                    my_histos["h_ht"                  +kv.first]->Fill(HT_trigger_pt45, w);
                    my_histos["h_htQCDCR"             +kv.first]->Fill(HT_NonIsoMuon_pt45, w);
                    my_histos["h_disc1"               +kv.first]->Fill(DoubleDisCo_disc1_NonIsoMuon_0l, w);
                    my_histos["h_disc2"               +kv.first]->Fill(DoubleDisCo_disc2_NonIsoMuon_0l, w);
                    if(kv.first.find("_passQCDCR") != std::string::npos)
                    {
                        my_2d_histos["h_njetsQCDCR_disc1"      +kv.first]->Fill(NNonIsoMuonJets_pt45, DoubleDisCo_disc1_NonIsoMuon_0l, w);
                        my_2d_histos["h_njetsQCDCR_disc2"      +kv.first]->Fill(NNonIsoMuonJets_pt45, DoubleDisCo_disc2_NonIsoMuon_0l, w);
                        my_2d_histos["h_disc1_disc2"           +kv.first]->Fill(DoubleDisCo_disc1_NonIsoMuon_0l, DoubleDisCo_disc2_NonIsoMuon_0l, w);
                    }
                }
                my_histos["h_mbl"                 +kv.first]->Fill(Mbl, w);
                for(const auto& l : GoodLeptons)
                {
                    my_histos["h_lPt"+kv.first]->Fill(l.second.Pt(), w);
                    my_histos["h_lEta"+kv.first]->Fill(l.second.Eta(), w);
                    my_histos["h_lPhi"+kv.first]->Fill(l.second.Phi(), w);
                }
                for(const auto& isoMuon : GoodNonIsoMuons)
                {
                    my_histos["h_isomPt"+kv.first]->Fill(isoMuon.second.Pt(), w);
                    my_histos["h_isomEta"+kv.first]->Fill(isoMuon.second.Eta(), w);
                    my_histos["h_isomPhi"+kv.first]->Fill(isoMuon.second.Phi(), w);
                }
                for(const auto& mbl : MblVec)
                {
                    my_histos["h_allMbl"+kv.first]->Fill(mbl, w);
                }
                for(unsigned int j = 0; j < Jets.size(); j++)
                {
                    if(kv.first.find("passQCDCR") != std::string::npos)
                    {
                        if(kv.first.find("_0l") != std::string::npos || kv.first.find("pt45") != std::string::npos)
                        { 
                            if(!NonIsoMuonJets_pt45[j]) continue;
                            my_histos["h_jPt"+kv.first]->Fill(Jets.at(j).Pt(), w);
                            my_histos["h_jEta"+kv.first]->Fill(Jets.at(j).Eta(), w);
                            my_histos["h_jPhi"+kv.first]->Fill(Jets.at(j).Phi(), w);
                        }
                        else if(kv.first.find("_1l") != std::string::npos || kv.first.find("pt30") != std::string::npos)
                        { 
                            if(!NonIsoMuonJets_pt30[j]) continue;
                            my_histos["h_jPt"+kv.first]->Fill(Jets.at(j).Pt(), w);
                            my_histos["h_jEta"+kv.first]->Fill(Jets.at(j).Eta(), w);
                            my_histos["h_jPhi"+kv.first]->Fill(Jets.at(j).Phi(), w);
                        }
                    }
                    else
                    {
                        if(kv.first.find("_0l") != std::string::npos || kv.first.find("pt45") != std::string::npos)
                        { 
                            if(!GoodJets_pt45[j]) continue;
                            my_histos["h_jPt"+kv.first]->Fill(Jets.at(j).Pt(), w);
                            my_histos["h_jEta"+kv.first]->Fill(Jets.at(j).Eta(), w);
                            my_histos["h_jPhi"+kv.first]->Fill(Jets.at(j).Phi(), w);
                        }
                        else if(kv.first.find("_1l") != std::string::npos || kv.first.find("pt30") != std::string::npos)
                        { 
                            if(!GoodJets_pt30[j]) continue;
                            my_histos["h_jPt"+kv.first]->Fill(Jets.at(j).Pt(), w);
                            my_histos["h_jEta"+kv.first]->Fill(Jets.at(j).Eta(), w);
                            my_histos["h_jPhi"+kv.first]->Fill(Jets.at(j).Phi(), w);
                        }
                    }
                }
            }
        }

        // ------------
        // -- Cut flow
        // ------------
        /*
        if(true) my_histos["h_cutFlow"]->Fill(0.5, weight);
        if(true && pass_general) my_histos["h_cutFlow"]->Fill(1.5, weight);  
        if(true && pass_general && pass_1l) my_histos["h_cutFlow"]->Fill(2.5, weight);
        if(true && pass_general && pass_1l && pass_ht) my_histos["h_cutFlow"]->Fill(3.5, weight);
        if(true && pass_general && pass_1l && pass_ht && JetID) my_histos["h_cutFlow"]->Fill(4.5, weight);
        if(true && pass_general && pass_1l && pass_ht && JetID && pass_1btag_pt30) my_histos["h_cutFlow"]->Fill(5.5, weight);
        if(true && pass_general && pass_1l && pass_ht && JetID && pass_1btag_pt30 && pass_MBL) my_histos["h_cutFlow"]->Fill(6.5, weight);
        if(true && pass_general && pass_1l && pass_ht && JetID && pass_1btag_pt30 && pass_MBL && pass_njet_pt30) my_histos["h_cutFlow"]->Fill(7.5, weight);
        if(true && pass_general && pass_1l && pass_ht && JetID && pass_1btag_pt30 && pass_MBL && pass_njet_pt30 && passHEMVeto) my_histos["h_cutFlow"]->Fill(8.5, weight);   
        */
    } // end of event loop
}

void AnalyzeQCDCR::WriteHistos(TFile* outfile)
{
    outfile->cd();
    
    for(const auto& p : my_histos) 
    {
        p.second->Write();
    }
   
    for(const auto& p : my_2d_histos) 
    {
        p.second->Write();
    }
}
