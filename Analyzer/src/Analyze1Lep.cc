#define Analyze1Lep_cxx
#include "Analyzer/Analyzer/include/Analyze1Lep.h"
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

Analyze1Lep::Analyze1Lep() : initHistos(false)
{
}

void Analyze1Lep::InitHistos(const std::map<std::string, bool>& cutMap, const std::vector<TH1DInfo>& histInfos, 
                             const std::vector<TH2DInfo>& hist2DInfos,  const std::vector<TH2DProfileInfo>& hist2DProfileInfos,
                             const unsigned int nMVAJets)
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    my_histos.emplace("EventCounter", std::make_shared<TH1D>("EventCounter","EventCounter", 2, -1.1, 1.1 ) );

    // Declare all your histograms here, that way we can fill them for multiple chains
    my_histos.emplace("fwm2_top6_1l_ge7j_ge1b", std::make_shared<TH1D>("fwm2_top6_1l_ge7j_ge1b","fwm2_top6_1l_ge7j_ge1b", 50, 0, 1 ) );
    my_histos.emplace("fwm3_top6_1l_ge7j_ge1b", std::make_shared<TH1D>("fwm3_top6_1l_ge7j_ge1b","fwm3_top6_1l_ge7j_ge1b", 50, 0, 1 ) );
    my_histos.emplace("fwm4_top6_1l_ge7j_ge1b", std::make_shared<TH1D>("fwm4_top6_1l_ge7j_ge1b","fwm4_top6_1l_ge7j_ge1b", 50, 0, 1 ) );
    my_histos.emplace("fwm5_top6_1l_ge7j_ge1b", std::make_shared<TH1D>("fwm5_top6_1l_ge7j_ge1b","fwm5_top6_1l_ge7j_ge1b", 50, 0, 1 ) );
    my_histos.emplace("fwm6_top6_1l_ge7j_ge1b", std::make_shared<TH1D>("fwm6_top6_1l_ge7j_ge1b","fwm6_top6_1l_ge7j_ge1b", 50, 0, 1 ) );
    my_histos.emplace("fwm7_top6_1l_ge7j_ge1b", std::make_shared<TH1D>("fwm7_top6_1l_ge7j_ge1b","fwm7_top6_1l_ge7j_ge1b", 50, 0, 1 ) );
    my_histos.emplace("fwm8_top6_1l_ge7j_ge1b", std::make_shared<TH1D>("fwm8_top6_1l_ge7j_ge1b","fwm8_top6_1l_ge7j_ge1b", 50, 0, 1 ) );
    my_histos.emplace("fwm9_top6_1l_ge7j_ge1b", std::make_shared<TH1D>("fwm9_top6_1l_ge7j_ge1b","fwm9_top6_1l_ge7j_ge1b", 50, 0, 1 ) );
    my_histos.emplace("fwm10_top6_1l_ge7j_ge1b", std::make_shared<TH1D>("fwm10_top6_1l_ge7j_ge1b","fwm10_top6_1l_ge7j_ge1b", 50, 0, 1 ) );
    my_histos.emplace("jmt_ev0_top6_1l_ge7j_ge1b", std::make_shared<TH1D>("jmt_ev0_top6_1l_ge7j_ge1b","jmt_ev0_top6_1l_ge7j_ge1b", 50, 0, 1 ) );
    my_histos.emplace("jmt_ev1_top6_1l_ge7j_ge1b", std::make_shared<TH1D>("jmt_ev1_top6_1l_ge7j_ge1b","jmt_ev1_top6_1l_ge7j_ge1b", 50, 0, 1 ) );
    my_histos.emplace("jmt_ev2_top6_1l_ge7j_ge1b", std::make_shared<TH1D>("jmt_ev2_top6_1l_ge7j_ge1b","jmt_ev2_top6_1l_ge7j_ge1b", 50, 0, 1 ) );
    my_histos.emplace("GoodLeptons_pt_1", std::make_shared<TH1D>("GoodLeptons_pt_1","GoodLeptons_pt_1", 150, 0, 1500 ) );
    my_histos.emplace("GoodLeptons_eta_1", std::make_shared<TH1D>("GoodLeptons_eta_1","GoodLeptons_eta_1", 100, -6, 6 ) );
    my_histos.emplace("GoodLeptons_phi_1", std::make_shared<TH1D>("GoodLeptons_phi_1","GoodLeptons_phi_1", 80, -4, 4 ) );
    my_histos.emplace("GoodLeptons_m_1", std::make_shared<TH1D>("GoodLeptons_m_1","GoodLeptons_m_1", 20, 0, 200 ) );

    for(unsigned int i = 1; i <= nMVAJets ; i++)
    {
        my_histos.emplace("Jet_cm_pt_"+std::to_string(i)+"_1l_ge7j_ge1b",  std::make_shared<TH1D>(("Jet_cm_pt_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(),("Jet_cm_pt_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(), 150, 0, 1500 ));
        my_histos.emplace("Jet_cm_eta_"+std::to_string(i)+"_1l_ge7j_ge1b", std::make_shared<TH1D>(("Jet_cm_eta_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(),("Jet_cm_eta_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(), 100, -6, 6 ));
        my_histos.emplace("Jet_cm_phi_"+std::to_string(i)+"_1l_ge7j_ge1b", std::make_shared<TH1D>(("Jet_cm_phi_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(),("Jet_cm_phi_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(), 80, -4, 4 ));
        my_histos.emplace("Jet_cm_m_"+std::to_string(i)+"_1l_ge7j_ge1b",   std::make_shared<TH1D>(("Jet_cm_m_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(),("Jet_cm_m_"+std::to_string(i)+"_1l_ge7j_ge1b").c_str(), 20, 0, 200 ));

        my_2d_histos.emplace("Jet_cm_pt_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_pt_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b").c_str(),("Jet_cm_pt_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b").c_str(), 200, 0.0, 1.0, 300, 0, 1500));
        my_2d_histos.emplace("Jet_cm_eta_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_eta_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b").c_str(),("Jet_cm_eta_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b").c_str(), 200, 0.0, 1.0, 100, -6, 6));
        my_2d_histos.emplace("Jet_cm_phi_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_phi_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b").c_str(),("Jet_cm_phi_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b").c_str(), 200, 0.0, 1.0, 80, -4, 4));
        my_2d_histos.emplace("Jet_cm_m_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_m_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b").c_str(),("Jet_cm_m_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b").c_str(), 200, 0.0, 1.0, 40, 0, 200));

        my_2d_histos.emplace("Jet_cm_pt_"+std::to_string(i)+"_njets_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_pt_"+std::to_string(i)+"_njets_1l_ge7j_ge1b").c_str(),("Jet_cm_pt_"+std::to_string(i)+"_njets_1l_ge7j_ge1b").c_str(), 20, 0, 20, 300, 0, 1500));
        my_2d_histos.emplace("Jet_cm_eta_"+std::to_string(i)+"_njets_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_eta_"+std::to_string(i)+"_njets_1l_ge7j_ge1b").c_str(),("Jet_cm_eta_"+std::to_string(i)+"_njets_1l_ge7j_ge1b").c_str(), 20, 0, 20, 100, -6, 6));
        my_2d_histos.emplace("Jet_cm_phi_"+std::to_string(i)+"_njets_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_phi_"+std::to_string(i)+"_njets_1l_ge7j_ge1b").c_str(),("Jet_cm_phi_"+std::to_string(i)+"_njets_1l_ge7j_ge1b").c_str(), 20, 0, 20, 80, -4, 4));
        my_2d_histos.emplace("Jet_cm_m_"+std::to_string(i)+"_njets_1l_ge7j_ge1b", std::make_shared<TH2D>(("Jet_cm_m_"+std::to_string(i)+"_njets_1l_ge7j_ge1b").c_str(),("Jet_cm_m_"+std::to_string(i)+"_njets_1l_ge7j_ge1b").c_str(), 20, 0, 20, 40, 0, 200));
    }

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

        for(const auto& h2dProfile : hist2DProfileInfos)
        {
            my_2d_tp_histos.emplace(h2dProfile.name+mycut.first,
                                    std::make_shared<TProfile2D>((h2dProfile.name+mycut.first).c_str(),(h2dProfile.name+mycut.first).c_str(), 
                                                                 h2dProfile.nBinsX, h2dProfile.lowX, h2dProfile.highX, h2dProfile.nBinsY, h2dProfile.lowY, h2dProfile.highY, h2dProfile.lowZ, h2dProfile.highZ));
        }
    }

    my_histos.emplace( "h_cutFlow", std::make_shared<TH1D>("h_cutFlow", "h_cutFlow", 9,0,9));    
}

void Analyze1Lep::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    while( tr.getNextEvent() )
    {

        const auto& MET                       = tr.getVar<float>("MET");
        const auto& METPhi                    = tr.getVar<float>("METPhi");
        const auto& runtype                   = tr.getVar<std::string>("runtype");     
        const auto& filetag                   = tr.getVar<std::string>("filetag");
        const auto& Jets                      = tr.getVec<utility::LorentzVector>("Jets");
        const auto& GoodJets_pt30             = tr.getVec<bool>("GoodJets_pt30");
        const auto& NJet                      = tr.getVar<int>("NJets");
        const auto& NGoodJets_pt30            = tr.getVar<int>("NGoodJets_pt30");
        const auto& NNonIsoMuonJets_pt30      = tr.getVar<int>("NNonIsoMuonJets_pt30");
        const auto& NGoodBJets_pt30           = tr.getVar<int>("NGoodBJets_pt30");
        const auto& NGoodMuons                = tr.getVar<int>("NGoodMuons");
        const auto& NGoodElectrons            = tr.getVar<int>("NGoodElectrons");
        const auto& NGoodLeptons              = tr.getVar<int>("NGoodLeptons");
        const auto& GoodLeptons               = tr.getVec<std::pair<std::string, utility::LorentzVector>>("GoodLeptons");
        const auto& GoodNonIsoMuons           = tr.getVec<std::pair<std::string, utility::LorentzVector>>("GoodNonIsoMuons");
        const auto& HT_trigger_pt30           = tr.getVar<double>("HT_trigger_pt30");
        const auto& HT_NonIsoMuon_pt30        = tr.getVar<double>("HT_NonIsoMuon_pt30");
        const auto& JetID                     = tr.getVar<bool>("JetID");
        const auto& correct2018Split          = tr.getVar<bool>("correct2018Split");
        const auto& passTrigger               = tr.getVar<bool>("passTrigger");
        const auto& passTriggerMC             = tr.getVar<bool>("passTriggerMC");
        const auto& passMETFilters            = tr.getVar<bool>("passMETFilters");
        const auto& passMadHT                 = tr.getVar<bool>("passMadHT");
        const auto& passBaseline1l_Good       = tr.getVar<bool>("passBaseline1l_Good");
        const auto& passBaseline1e1m          = tr.getVar<bool>("passBaseline1e1m_Good");
        const auto& passBaselineGoodOffline1l = tr.getVar<bool>("passBaselineGoodOffline1l");
        const auto& passBaseline1l_NonIsoMuon = tr.getVar<bool>("passBaseline1l_NonIsoMuon");
        const auto& passHEMVeto               = tr.getVar<bool>("passHEMVeto");
        const auto& Mbl                       = tr.getVar<double>("Mbl");
        const auto& MblVec                    = tr.getVec<double>("MblVec");
        const auto& passBlind                 = tr.getVar<bool>("passBlindLep_Good");            
        const auto& deepESM_val               = tr.getVar<double>("deepESM_val");
        const auto& deepESM_valNonIsoMuon     = tr.getVar<double>("deepESM_valNonIsoMuon");
        const auto& deepESM_bin1              = tr.getVar<bool>("deepESM_bin1");
        const auto& deepESM_bin2              = tr.getVar<bool>("deepESM_bin2");
        const auto& deepESM_bin3              = tr.getVar<bool>("deepESM_bin3");
        const auto& deepESM_bin4              = tr.getVar<bool>("deepESM_bin4");
        const auto& deepESM_binNum            = tr.getVar<int>("deepESM_binNum");
        const auto& fwm2_top6                 = tr.getVar<double>("fwm2_top6");
        const auto& fwm3_top6                 = tr.getVar<double>("fwm3_top6");
        const auto& fwm4_top6                 = tr.getVar<double>("fwm4_top6");
        const auto& fwm5_top6                 = tr.getVar<double>("fwm5_top6");
        const auto& fwm6_top6                 = tr.getVar<double>("fwm6_top6");
        const auto& fwm7_top6                 = tr.getVar<double>("fwm7_top6");
        const auto& fwm8_top6                 = tr.getVar<double>("fwm8_top6");
        const auto& fwm9_top6                 = tr.getVar<double>("fwm9_top6");
        const auto& fwm10_top6                = tr.getVar<double>("fwm10_top6");
        const auto& jmt_ev0_top6              = tr.getVar<double>("jmt_ev0_top6");
        const auto& jmt_ev1_top6              = tr.getVar<double>("jmt_ev1_top6");
        const auto& jmt_ev2_top6              = tr.getVar<double>("jmt_ev2_top6");
        const auto& Jets_cm_top6              = tr.getVec<utility::LorentzVector>("Jets_cm_top6");
        const auto& eventCounter              = tr.getVar<int>("eventCounter");
        const auto& nMVAJets                  = tr.getVar<unsigned int>("nMVAJets");

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
            const auto& Weight  = tr.getVar<float>("Weight");
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

        int NGenJets = 0, NGenJets_pt30 = 0;
        if(runtype == "MC")
        {
            //Get the numer of GenJets 
            const auto& GenJets = tr.getVec<utility::LorentzVector>("GenJets");           
            for(const auto& lv : GenJets)
            {
                NGenJets++;
                if(lv.Pt() > 30)
                {
                    NGenJets_pt30++;
                }
            }
        }

        // -------------------------------
        // -- Define cuts
        // -------------------------------
        bool pass_general    = passTriggerMC && passTrigger && passMadHT && passBlind && passMETFilters && passHEMVeto && correct2018Split;
        //bool pass_0l         = NGoodLeptons == 0;
        bool pass_1l         = NGoodLeptons == 1;
        bool pass_ht         = HT_trigger_pt30 > 300;
        bool pass_MBL        = (50 < Mbl && Mbl < 250);
        //bool pass_1e_1m      = (NGoodLeptons == 2) ? GoodLeptons[0].first != GoodLeptons[1].first : false;
        bool pass_njet_pt30  = NGoodJets_pt30 >= 7;
        bool pass_1btag_pt30 = NGoodBJets_pt30 >= 1;
        bool pass_2btag_pt30 = NGoodBJets_pt30 >= 2;
        bool pass_1e         = false;
        bool pass_1m         = false;
        bool pass_lBarrel    = false;
        if(pass_1l)
        { 
            if(GoodLeptons[0].first == "e") pass_1e = true;
            if(GoodLeptons[0].first == "m") pass_1m = true;
            pass_lBarrel = abs( GoodLeptons[0].second.Eta() ) <= 1.2;
        }
        //bool pass_5to6njet_pt30 = (NGoodJets_pt30 == 5 || NGoodJets_pt30 == 6);
        
        bool passBaseline1l_AllJets = passBaselineGoodOffline1l &&
                                      passTrigger               &&
                                      passTriggerMC             &&
                                      passBlind                 &&
            (
                ((runtype != "Data" || filetag.find("Data_SingleMuon")     != std::string::npos) && NGoodMuons == 1     && NGoodElectrons == 0)
                                                                                     ||
                ((runtype != "Data" || filetag.find("Data_SingleElectron") != std::string::npos) && NGoodElectrons == 1 && NGoodMuons == 0)
            );
        
        // ------------------------------------------------
        // --  Temporary Home of the W+Jets Control Region
        // ------------------------------------------------
        double mT = -1;
        bool pass_mT = false, pass_Mbl_all = false;
        if(pass_1l)
        {
            utility::LorentzVector metLV;
            metLV.SetPt(MET); metLV.SetEta(0.0); metLV.SetPhi(METPhi); metLV.SetE(MET);
            mT = utility::calcMT(GoodLeptons[0].second, metLV);
            pass_mT = mT > 50 && mT < 110;
        
            for(unsigned int i = 0; i < Jets.size(); ++i)
            {
                if(!GoodJets_pt30[i]) continue;
                double mbl = (GoodLeptons[0].second+Jets[i]).M();
                if(mbl > 50 && mbl < 250)
                {
                    pass_Mbl_all = true;
                    break;
                }
            }            
        }

        const auto& NGoodBJets_pt30_loose = tr.getVar<int>("NGoodBJets_pt30_loose");        
        bool pass_0b_loose = NGoodBJets_pt30_loose == 0;

        bool passBaseline1l_WCR = JetID                   &&
                                  passMadHT               &&
                                  passTrigger             &&
                                  passTriggerMC           &&
                                  HT_trigger_pt30 > 300   &&
                                  pass_mT                 && 
                                  pass_0b_loose           && 
                                  !pass_Mbl_all           && 
                                  MET>30                  &&
            (
                ((runtype != "Data" || filetag.find("Data_SingleMuon")     != std::string::npos) && NGoodMuons == 1     && NGoodElectrons == 0)
                                                                                     ||
                ((runtype != "Data" || filetag.find("Data_SingleElectron") != std::string::npos) && NGoodElectrons == 1 && NGoodMuons == 0)
            );

        bool evenEvent = tr.getEvtNum() % 2 == 0;

        // -------------------
        // --- Fill Histos ---
        // -------------------                        
        const std::map<std::string, bool> cut_map_1l 
        {
            {""                                      , pass_general                                                                             },
            {"_1l_HT300"                             , pass_general && pass_1l && pass_ht                                                       },
            {"_1l_HT300_ge7j"                        , pass_general && pass_1l && pass_ht && pass_njet_pt30  && JetID                           },
            {"_1l_HT300_ge1b"                        , pass_general && pass_1l && pass_ht && pass_1btag_pt30 && JetID                           },
            {"_1l_HT300_ge2b"                        , pass_general && pass_1l && pass_ht && pass_2btag_pt30 && JetID                           },
            {"_1l_HT300_ge7j_ge1b"                   , pass_general && pass_1l && pass_ht && pass_njet_pt30  && pass_1btag_pt30 && JetID        },
            {"_1l_HT300_ge1b_Mbl"                    , pass_general && passBaseline1l_AllJets                                                   },
            {"_1l_HT300_ge4j_ge1b_Mbl"               , pass_general && passBaseline1l_AllJets && NGoodJets_pt30 >= 4                            },
            {"_1e_HT300_ge4j_ge1b_Mbl"               , pass_general && passBaseline1l_AllJets && NGoodJets_pt30 >= 4 && pass_1e                 },
            {"_1m_HT300_ge4j_ge1b_Mbl"               , pass_general && passBaseline1l_AllJets && NGoodJets_pt30 >= 4 && pass_1m                 },
            {"_1l_HT300_ge7j_ge1b_Mbl"               , pass_general && passBaseline1l_Good                                                      },                         
            {"_1e_HT300_ge7j_ge1b_Mbl"               , pass_general && passBaseline1l_Good && pass_1e                                           },                         
            {"_1m_HT300_ge7j_ge1b_Mbl"               , pass_general && passBaseline1l_Good && pass_1m                                           },   
            {"_1l_HT300_ge7j_ge1b_Mbl_posWeight"     , pass_general && passBaseline1l_Good && weight >  0.0                                     },                         
            {"_1l_HT300_ge7j_ge1b_Mbl_negWeight"     , pass_general && passBaseline1l_Good && weight <= 0.0                                     },                         
            {"_1l_HT300_ge7j_ge1b_Mbl_noHTWeight"    , pass_general && passBaseline1l_Good                                                      }, 
            {"_1l_HT300_ge7j_ge1b_Mbl_noLepWeight"   , pass_general && passBaseline1l_Good                                                      },                         
            {"_1l_HT300_ge7j_ge1b_Mbl_noBTagWeight"  , pass_general && passBaseline1l_Good                                                      },                         
            {"_1l_HT300_ge7j_ge2b_Mbl"               , pass_general && passBaseline1l_Good && pass_2btag_pt30                                   },
            {"_1l_HT300_ge7j_ge1b_Mbl_lBarrel"       , pass_general && passBaseline1l_Good && pass_lBarrel                                      },
            {"_1e_HT300_ge7j_ge1b_Mbl_lBarrel"       , pass_general && passBaseline1l_Good && pass_lBarrel && pass_1e                           },
            {"_1m_HT300_ge7j_ge1b_Mbl_lBarrel"       , pass_general && passBaseline1l_Good && pass_lBarrel && pass_1m                           },
            {"_1l_HT300_ge7j_ge1b_Mbl_lEndCap"       , pass_general && passBaseline1l_Good && !pass_lBarrel                                     },
            {"_1e_HT300_ge7j_ge1b_Mbl_lEndCap"       , pass_general && passBaseline1l_Good && !pass_lBarrel && pass_1e                          },
            {"_1m_HT300_ge7j_ge1b_Mbl_lEndCap"       , pass_general && passBaseline1l_Good && !pass_lBarrel && pass_1m                          },
            {"_1l_HT300_ge7j_ge1b_Mbl_d1"            , pass_general && passBaseline1l_Good && deepESM_bin1                                      },                         
            {"_1l_HT300_ge7j_ge1b_Mbl_d2"            , pass_general && passBaseline1l_Good && deepESM_bin2                                      },                         
            {"_1l_HT300_ge7j_ge1b_Mbl_d3"            , pass_general && passBaseline1l_Good && deepESM_bin3                                      },                         
            {"_1l_HT300_ge7j_ge1b_Mbl_d4"            , pass_general && passBaseline1l_Good && deepESM_bin4                                      },                         
            {"_1l_HT300_1j_ge1b_Mbl"                 , pass_general && passBaseline1l_AllJets && NGoodJets_pt30 == 1                            },
            {"_1l_HT300_2j_ge1b_Mbl"                 , pass_general && passBaseline1l_AllJets && NGoodJets_pt30 == 2                            },
            {"_1l_HT300_3j_ge1b_Mbl"                 , pass_general && passBaseline1l_AllJets && NGoodJets_pt30 == 3                            },
            {"_1l_HT300_4j_ge1b_Mbl"                 , pass_general && passBaseline1l_AllJets && NGoodJets_pt30 == 4                            },
            {"_1l_HT300_5j_ge1b_Mbl"                 , pass_general && passBaseline1l_AllJets && NGoodJets_pt30 == 5                            },
            {"_1l_HT300_6j_ge1b_Mbl"                 , pass_general && passBaseline1l_AllJets && NGoodJets_pt30 == 6                            },
            {"_1l_HT300_7j_ge1b_Mbl"                 , pass_general && passBaseline1l_Good && NGoodJets_pt30 == 7                               },
            {"_1l_HT300_8j_ge1b_Mbl"                 , pass_general && passBaseline1l_Good && NGoodJets_pt30 == 8                               },
            {"_1l_HT300_9j_ge1b_Mbl"                 , pass_general && passBaseline1l_Good && NGoodJets_pt30 == 9                               },
            {"_1l_HT300_10j_ge1b_Mbl"                , pass_general && passBaseline1l_Good && NGoodJets_pt30 == 10                              },
            {"_1l_HT300_11j_ge1b_Mbl"                , pass_general && passBaseline1l_Good && NGoodJets_pt30 == 11                              },
            {"_1l_HT300_12j_ge1b_Mbl"                , pass_general && passBaseline1l_Good && NGoodJets_pt30 == 12                              },
            {"_1l_HT300_13j_ge1b_Mbl"                , pass_general && passBaseline1l_Good && NGoodJets_pt30 == 13                              },
            {"_1l_HT300_14j_ge1b_Mbl"                , pass_general && passBaseline1l_Good && NGoodJets_pt30 == 14                              },
            {"_1l_HT300_15j_ge1b_Mbl"                , pass_general && passBaseline1l_Good && NGoodJets_pt30 == 15                              },
            {"_1l_HT300_5j_ge1b_Mbl_htCorr"          , pass_general && passBaseline1l_AllJets && NGoodJets_pt30 == 5                            },
            {"_1l_HT300_6j_ge1b_Mbl_htCorr"          , pass_general && passBaseline1l_AllJets && NGoodJets_pt30 == 6                            },
            {"_1l_HT300_7j_ge1b_Mbl_htCorr"          , pass_general && passBaseline1l_Good && NGoodJets_pt30 == 7                               },
            {"_1l_HT300_8j_ge1b_Mbl_htCorr"          , pass_general && passBaseline1l_Good && NGoodJets_pt30 == 8                               },
            {"_1l_HT300_ge7j_ge1b_Mbl_htCorr"        , pass_general && passBaseline1l_Good                                                      },

            {"_1l_0b_ge300ht_50to110mt_ge30MET"      , pass_general && passBaseline1l_WCR                                                       },
            {"_1l_0b_ge300ht_50to110mt_ge30MET_even" , pass_general && passBaseline1l_WCR && evenEvent                                          },
            {"_1l_0b_ge300ht_50to110mt_ge30MET_odd"  , pass_general && passBaseline1l_WCR && !evenEvent                                         },
            {"_1e_1m_ge2b_le5j"                      , pass_general && passBaseline1e1m                                                         },
            {"_passQCDCR"                            , passBaseline1l_NonIsoMuon                                                                },
        };

        std::vector<TH1DInfo> histInfos = {
            {    "h_njets",             20,   0.0,   20.0},
            {"blind_njets",             20,   0.0,   20.0},
            {    "h_njetsQCDCR",        20,   0.0,   20.0},
            {    "h_ngjets",            20,   0.0,   20.0},
            {    "h_ngjets_pt30",       20,   0.0,   20.0},
            {    "h_nb",                10,   0.0,   10.0},
            {"blind_nb",                10,   0.0,   10.0},
            {    "h_deepESM",          200,   0.0,    1.0},
            {"blind_deepESM",          200,   0.0,    1.0},
            {    "h_deepESMQCDCR",     200,   0.0,    1.0},
            {    "h_deepESMMerged",      4,   0.5,    4.5},
            {"blind_deepESMMerged",      4,   0.5,    4.5},
            {    "h_ht",               500,   0.0, 5000.0},
            {    "h_htQCDCR",          500,   0.0, 5000.0},
            {"blind_ht",               500,   0.0, 5000.0},
            {    "h_mbl",              300,   0.0,  300.0},
            {"blind_mbl",              300,   0.0,  300.0},
            {    "h_lPt",              200,   0.0, 2000.0},
            {"blind_lPt",              200,   0.0, 2000.0},
            {    "h_lEta",             200,  -6.0,    6.0},
            {"blind_lEta",             200,  -6.0,    6.0},
            {    "h_lPhi",             200,  -4.0,    4.0},
            {"blind_lPhi",             200,  -4.0,    4.0},
            {    "h_isomPt",           200,   0.0, 2000.0},
            {    "h_isomEta",          200,  -6.0,    6.0},
            {    "h_isomPhi",          200,  -4.0,    4.0},
            {    "h_jPt",              200,   0.0, 2000.0},
            {"blind_jPt",              200,   0.0, 2000.0},
            {    "h_jEta",             200,  -6.0,    6.0},
            {"blind_jEta",             200,  -6.0,    6.0},
            {    "h_jPhi",             200,  -4.0,    4.0},
            {"blind_jPhi",             200,  -4.0,    4.0},
            {    "h_allMbl",           300,   0.0,  300.0},            
            {"blind_allMbl",           300,   0.0,  300.0},
            {"h_weight",               200,  -5.0,    5.0},
            {"h_leptonweight",         200,  -5.0,    5.0},
            {"h_pileupWeight",         200,  -5.0,   20.0},
            {"h_bTagWeight",           200,  -5.0,   20.0},
            {"h_htDerivedweight",      200,  -5.0,    5.0},
            {"h_prefiringScaleFactor", 200,  -5.0,    5.0},
            {"h_eventweight",          200,  -5.0,   50.0},            
        };

        std::vector<TH2DInfo> hist2DInfos = {
            {    "h_njets_deepESM", 15,    0,   15, 100,   0.0,   1.0},
            {"blind_njets_deepESM", 15,    0,   15, 100,   0.0,   1.0},
            {    "h_njets_mbl",     15,    0,   15, 100,   0.0, 300.0},
            {"blind_njets_mbl",     15,    0,   15, 100,   0.0, 300.0},
            {    "h_ht_deepESM",   100,    0, 3000, 100,   0.0,   1.0},
            {"blind_ht_deepESM",   100,    0, 3000, 100,   0.0,   1.0},
            {    "h_lEta_lPhi",    100, -6.0,  6.0, 100,  -3.2,   3.2},
            {"blind_lEta_lPhi",    100, -6.0,  6.0, 100,  -3.2,   3.2},
            {    "h_jEta_jPhi",    100, -6.0,  6.0, 100,  -3.2,   3.2},
            {"blind_jEta_jPhi",    100, -6.0,  6.0, 100,  -3.2,   3.2},
            {    "h_lEta_nb",      100, -6.0,  6.0,  10,   0.0,  10.0},
            {"blind_lEta_nb",      100, -6.0,  6.0,  10,   0.0,  10.0},
        };

        std::vector<TH2DProfileInfo> hist2DProfileInfos = {
            {"h_njets_deepESMMerged_preFireSF", 15,0,15, 4,0.5,4.5, 0,1.0},
        };

        // Initialize Histograms
        if(!initHistos)
        {
            InitHistos(cut_map_1l, histInfos, hist2DInfos, hist2DProfileInfos, nMVAJets);
            initHistos = true;
        }

        my_histos["EventCounter"]->Fill(eventCounter);

        for(auto& kv : cut_map_1l)
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
                my_histos["h_njets"               +kv.first]->Fill(NGoodJets_pt30, w);
                my_histos["h_ngjets"              +kv.first]->Fill(NGenJets, eventweight);
                my_histos["h_ngjets_pt30"         +kv.first]->Fill(NGenJets_pt30, eventweight);
                my_histos["h_njetsQCDCR"          +kv.first]->Fill(NNonIsoMuonJets_pt30, w);
                my_histos["h_nb"                  +kv.first]->Fill(NGoodBJets_pt30, w);
                my_histos["h_deepESM"             +kv.first]->Fill(deepESM_val, w);
                my_histos["h_deepESMQCDCR"        +kv.first]->Fill(deepESM_valNonIsoMuon, w);
                my_histos["h_deepESMMerged"       +kv.first]->Fill(deepESM_binNum, w);
                my_histos["h_ht"                  +kv.first]->Fill(HT_trigger_pt30, w);
                my_histos["h_htQCDCR"             +kv.first]->Fill(HT_NonIsoMuon_pt30, w);
                my_histos["h_mbl"                 +kv.first]->Fill(Mbl, w);
                my_histos["h_weight"              +kv.first]->Fill(weight, w);
                my_histos["h_leptonweight"        +kv.first]->Fill(leptonweight, w);
                my_histos["h_pileupWeight"        +kv.first]->Fill(pileupWeight, w);
                my_histos["h_bTagWeight"          +kv.first]->Fill(bTagWeight, w);
                my_histos["h_htDerivedweight"     +kv.first]->Fill(htDerivedweight, w);
                my_histos["h_prefiringScaleFactor"+kv.first]->Fill(prefiringScaleFactor, w);
                my_histos["h_eventweight"         +kv.first]->Fill(eventweight, w);
                for(const auto& l : GoodLeptons)
                {
                    my_histos["h_lPt"+kv.first]->Fill(l.second.Pt(), w);
                    my_histos["h_lEta"+kv.first]->Fill(l.second.Eta(), w);
                    my_histos["h_lPhi"+kv.first]->Fill(l.second.Phi(), w);
                    my_2d_histos["h_lEta_lPhi"+kv.first]->Fill(l.second.Eta(), l.second.Phi(), w);
                    my_2d_histos["h_lEta_nb"+kv.first]->Fill(l.second.Eta(), NGoodBJets_pt30, w);
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
                    if(!GoodJets_pt30[j]) continue;
                    my_histos["h_jPt"+kv.first]->Fill(Jets.at(j).Pt(), w);
                    my_histos["h_jEta"+kv.first]->Fill(Jets.at(j).Eta(), w);
                    my_histos["h_jPhi"+kv.first]->Fill(Jets.at(j).Phi(), w);
                    my_2d_histos["h_jEta_jPhi"+kv.first]->Fill(Jets.at(j).Eta(), Jets.at(j).Phi(), w);
                }
                my_2d_histos["h_njets_deepESM"+kv.first]->Fill(NGoodJets_pt30, deepESM_val, w);
                my_2d_histos["h_njets_mbl"+kv.first]->Fill(NGoodJets_pt30, Mbl, w);
                my_2d_histos["h_ht_deepESM"+kv.first]->Fill(HT_trigger_pt30, deepESM_val, w);
                my_2d_tp_histos["h_njets_deepESMMerged_preFireSF"+kv.first]->Fill(NJet, deepESM_binNum, prefiringScaleFactor, w);

                if( NGoodJets_pt30 <= 8 )
                {
                    my_histos["blind_njets"         +kv.first]->Fill(NGoodJets_pt30, w);
                    my_histos["blind_nb"            +kv.first]->Fill(NGoodBJets_pt30, w);
                    my_histos["blind_deepESM"       +kv.first]->Fill(deepESM_val, w);
                    my_histos["blind_deepESMMerged" +kv.first]->Fill(deepESM_binNum, w);
                    my_histos["blind_ht"            +kv.first]->Fill(HT_trigger_pt30, w);
                    my_histos["blind_mbl"           +kv.first]->Fill(Mbl, w);
                    for(const auto l : GoodLeptons)
                    {
                        my_histos["blind_lPt"+kv.first]->Fill(l.second.Pt(), w);
                        my_histos["blind_lEta"+kv.first]->Fill(l.second.Eta(), w);
                        my_histos["blind_lPhi"+kv.first]->Fill(l.second.Phi(), w);
                        my_2d_histos["blind_lEta_lPhi"+kv.first]->Fill(l.second.Eta(), l.second.Phi(), w);
                        my_2d_histos["blind_lEta_nb"+kv.first]->Fill(l.second.Eta(), NGoodBJets_pt30, w);
                    }
                    for(const auto& mbl : MblVec)
                    {
                        my_histos["blind_allMbl"+kv.first]->Fill(mbl, w);
                    }
                    for(unsigned int j = 0; j < Jets.size(); j++)
                    {
                        if(!GoodJets_pt30[j]) continue;
                        my_histos["blind_jPt"+kv.first]->Fill(Jets.at(j).Pt(), w);
                        my_histos["blind_jEta"+kv.first]->Fill(Jets.at(j).Eta(), w);
                        my_histos["blind_jPhi"+kv.first]->Fill(Jets.at(j).Phi(), w);
                        my_2d_histos["blind_jEta_jPhi"+kv.first]->Fill(Jets.at(j).Eta(), Jets.at(j).Phi(), w);
                    }
                    my_2d_histos["blind_njets_deepESM"+kv.first]->Fill(NGoodJets_pt30, deepESM_val, w);
                    my_2d_histos["blind_njets_mbl"+kv.first]->Fill(NGoodJets_pt30, Mbl, w);
                    my_2d_histos["blind_ht_deepESM"+kv.first]->Fill(HT_trigger_pt30, deepESM_val, w);
                }
            }
        }

        // Fill histos for deepESM training
        if(passBaseline1l_Good)
        {
            my_histos["fwm2_top6_1l_ge7j_ge1b"]->Fill(fwm2_top6, weight);
            my_histos["fwm3_top6_1l_ge7j_ge1b"]->Fill(fwm3_top6, weight);
            my_histos["fwm4_top6_1l_ge7j_ge1b"]->Fill(fwm4_top6, weight);
            my_histos["fwm5_top6_1l_ge7j_ge1b"]->Fill(fwm5_top6, weight);
            my_histos["fwm6_top6_1l_ge7j_ge1b"]->Fill(fwm6_top6, weight);
            my_histos["fwm7_top6_1l_ge7j_ge1b"]->Fill(fwm7_top6, weight);
            my_histos["fwm8_top6_1l_ge7j_ge1b"]->Fill(fwm8_top6, weight);
            my_histos["fwm9_top6_1l_ge7j_ge1b"]->Fill(fwm9_top6, weight);
            my_histos["fwm10_top6_1l_ge7j_ge1b"]->Fill(fwm10_top6, weight);
            my_histos["jmt_ev0_top6_1l_ge7j_ge1b"]->Fill(jmt_ev0_top6, weight);
            my_histos["jmt_ev1_top6_1l_ge7j_ge1b"]->Fill(jmt_ev1_top6, weight);
            my_histos["jmt_ev2_top6_1l_ge7j_ge1b"]->Fill(jmt_ev2_top6, weight);
            my_histos["GoodLeptons_pt_1"]->Fill(GoodLeptons.at(0).second.Pt(), weight);
            my_histos["GoodLeptons_eta_1"]->Fill(GoodLeptons.at(0).second.Eta(), weight);
            my_histos["GoodLeptons_phi_1"]->Fill(GoodLeptons.at(0).second.Phi(), weight);
            my_histos["GoodLeptons_m_1"]->Fill(GoodLeptons.at(0).second.M(), weight);

            unsigned int nJets = Jets_cm_top6.size();
            unsigned int iVec  = 0;
            for(unsigned int i = 1; i <= nMVAJets; i++)
            {
                double pt  = (iVec < nJets) ? static_cast<double>(Jets_cm_top6.at(iVec).Pt())  : 0.0;
                double eta = (iVec < nJets) ? static_cast<double>(Jets_cm_top6.at(iVec).Eta()) : 0.0;
                double phi = (iVec < nJets) ? static_cast<double>(Jets_cm_top6.at(iVec).Phi()) : 0.0;
                double m   = (iVec < nJets) ? static_cast<double>(Jets_cm_top6.at(iVec).M())   : 0.0;

                my_histos["Jet_cm_pt_"+std::to_string(i)+"_1l_ge7j_ge1b"]->Fill(pt, weight);
                my_histos["Jet_cm_eta_"+std::to_string(i)+"_1l_ge7j_ge1b"]->Fill(eta, weight);
                my_histos["Jet_cm_phi_"+std::to_string(i)+"_1l_ge7j_ge1b"]->Fill(phi, weight);
                my_histos["Jet_cm_m_"+std::to_string(i)+"_1l_ge7j_ge1b"]->Fill(m, weight);
                my_2d_histos["Jet_cm_pt_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b"]->Fill(deepESM_val, pt, weight);
                my_2d_histos["Jet_cm_eta_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b"]->Fill(deepESM_val, eta, weight);
                my_2d_histos["Jet_cm_phi_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b"]->Fill(deepESM_val, phi, weight);
                my_2d_histos["Jet_cm_m_"+std::to_string(i)+"_deepESM_1l_ge7j_ge1b"]->Fill(deepESM_val, m, weight);
                my_2d_histos["Jet_cm_pt_"+std::to_string(i)+"_njets_1l_ge7j_ge1b"]->Fill(NGoodJets_pt30, pt, weight);
                my_2d_histos["Jet_cm_eta_"+std::to_string(i)+"_njets_1l_ge7j_ge1b"]->Fill(NGoodJets_pt30, eta, weight);
                my_2d_histos["Jet_cm_phi_"+std::to_string(i)+"_njets_1l_ge7j_ge1b"]->Fill(NGoodJets_pt30, phi, weight);
                my_2d_histos["Jet_cm_m_"+std::to_string(i)+"_njets_1l_ge7j_ge1b"]->Fill(NGoodJets_pt30, m, weight);

                iVec++;
            }
        }

        // ------------
        // -- Cut flow
        // ------------
        if(true) my_histos["h_cutFlow"]->Fill(0.5, weight);
        if(true && pass_general) my_histos["h_cutFlow"]->Fill(1.5, weight);  
        if(true && pass_general && pass_1l) my_histos["h_cutFlow"]->Fill(2.5, weight);
        if(true && pass_general && pass_1l && pass_ht) my_histos["h_cutFlow"]->Fill(3.5, weight);
        if(true && pass_general && pass_1l && pass_ht && JetID) my_histos["h_cutFlow"]->Fill(4.5, weight);
        if(true && pass_general && pass_1l && pass_ht && JetID && pass_1btag_pt30) my_histos["h_cutFlow"]->Fill(5.5, weight);
        if(true && pass_general && pass_1l && pass_ht && JetID && pass_1btag_pt30 && pass_MBL) my_histos["h_cutFlow"]->Fill(6.5, weight);
        if(true && pass_general && pass_1l && pass_ht && JetID && pass_1btag_pt30 && pass_MBL && pass_njet_pt30) my_histos["h_cutFlow"]->Fill(7.5, weight);
        if(true && pass_general && pass_1l && pass_ht && JetID && pass_1btag_pt30 && pass_MBL && pass_njet_pt30 && passHEMVeto) my_histos["h_cutFlow"]->Fill(8.5, weight);   
    } // end of event loop
}

void Analyze1Lep::WriteHistos(TFile* outfile)
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

    for(const auto& p : my_tp_histos) 
    {
        p.second->Write();
    }

    for(const auto& p : my_2d_tp_histos)
    {
        p.second->Write();
    }
    
    for(const auto& p : my_efficiencies) 
    {
        p.second->Write();
    }    
}
