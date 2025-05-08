#define MakeCutflow_cxx
#include "Analyzer/Analyzer/include/MakeCutflow.h"
#include "Framework/Framework/include/Utility.h"
#include "NTupleReader/include/NTupleReader.h"

#include <iostream>

MakeCutflow::MakeCutflow()
{
    InitHistos();
}

//Define all your histograms here. 
void MakeCutflow::InitHistos()
{
    TH1::SetDefaultSumw2();

    //This event counter histogram is necessary so that we know that all the condor jobs ran successfully. If not, when you use the hadder script, you will see a discrepancy in red as the files are being hadded.
    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ;

    //Define 1D histograms
    my_histos.emplace( "h_njets_0l", std::make_shared<TH1D>( "h_njets_0l", "h_njets_0l", 8, 8, 16 ) ) ;
    my_histos.emplace( "h_njets_1l", std::make_shared<TH1D>( "h_njets_1l", "h_njets_1l", 8, 7, 15 ) ) ;
    my_histos.emplace( "h_njets_2l", std::make_shared<TH1D>( "h_njets_2l", "h_njets_2l", 8, 6, 14 ) ) ;

    my_histos.emplace("h_cutflow_0l", std::make_shared<TH1D>("h_cutflow_0l", "h_cutflow_0l", 15, 0, 15));
    my_histos.emplace("h_cutflow_1l", std::make_shared<TH1D>("h_cutflow_1l", "h_cutflow_1l", 12, 0, 12));
    my_histos.emplace("h_cutflow_2l", std::make_shared<TH1D>("h_cutflow_2l", "h_cutflow_2l", 13, 0, 13));
}

//Put everything you want to do per event here.
void MakeCutflow::Loop(NTupleReader& tr, double, int maxevents, bool)
{

    while( tr.getNextEvent() )
    {
        //This is added to count the number of events- do not change the next two lines.
        const auto& eventCounter        = tr.getVar<int>("eventCounter");
        my_histos["EventCounter"]->Fill( eventCounter );

        //--------------------------------------------------
        //-- Print Event Number 
        //--------------------------------------------------
        
        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() & (10000 == 0) ) printf( " Event %i\n", tr.getEvtNum() );
        
        const auto& runtype             = tr.getVar<std::string>("runtype");     
        const auto& filetag             = tr.getVar<std::string>("filetag");
        const auto& GoodLeptons         = tr.getVec<std::pair<std::string, utility::LorentzVector>>("GoodLeptons");

        const auto& JetID               = tr.getVar<bool>("JetID");
        const auto& passMETFilters      = tr.getVar<bool>("passMETFilters");
        const auto& NGoodLeptons        = tr.getVar<int>("NGoodLeptons");
        const auto& passTrigger         = tr.getVar<bool>("passTrigger");
        const auto& passTriggerMC       = tr.getVar<bool>("passTriggerMC");
        const auto& passElectronHEMveto = tr.getVar<bool>("passElectronHEMveto");
        const auto& passTriggerHadMC    = tr.getVar<bool>("passTriggerHadMC");

        // get variables for 0-Lepton
        const auto& NGoodJets_pt45         = tr.getVar<int>("NGoodJets_pt45"    );
        const auto& NGoodBJets_pt45        = tr.getVar<int>("NGoodBJets_pt45"   );
        const auto& dR_bjets               = tr.getVar<double>("dR_bjets"       );
        const auto& ntops                  = tr.getVar<int>("ntops"             );
        // get variables for 1-Lepton
        const auto& Mbl                    = tr.getVar<double>("Mbl"            );
        // get variables for 2-Lepton
        const auto& NGoodMuons             = tr.getVar<int>("NGoodMuons"        );
        const auto& NGoodElectrons         = tr.getVar<int>("NGoodElectrons"    );
        const auto& GoodLeptonsCharge      = tr.getVec<int>("GoodLeptonsCharge" );
        const auto& onZ                    = tr.getVar<bool>("onZ"              );
        const auto& NNonIsoMuons           = tr.getVar<int>("NNonIsoMuons"      );
        const auto& NGoodBJets_pt30     = tr.getVar<int>("NGoodBJets_pt30");
        const auto& HT_trigger_pt30     = tr.getVar<double>("HT_trigger_pt30");
        const auto& NGoodJets_pt30      = tr.getVar<int>("NGoodJets_pt30");
        
        const auto& passMadHT           = tr.getVar<bool>("passMadHT");

        const auto& passBaseline0l        = tr.getVar<bool>("passBaseline0l_Good");
        const auto& passBaseline1l        = tr.getVar<bool>("passBaseline1l_Good");
        const auto& passBaseline2l        = tr.getVar<bool>("passBaseline2l_Good");
      
        if(maxevents != -1 && tr.getEvtNum() >= maxevents) break;        
        if ( tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() ) ;

        double weight_0l = tr.getVar<double>("TotalWeight_0l");
        double weight_1l = tr.getVar<double>("TotalWeight_1l");
        double weight_2l = tr.getVar<double>("TotalWeight_2l");
 
        if (passBaseline0l) {
            my_histos["h_njets_0l"]->Fill(NGoodJets_pt30, weight_0l);
        } 
        if (passBaseline1l) {
            my_histos["h_njets_1l"]->Fill(NGoodJets_pt30, weight_1l);
        }
        if (passBaseline2l) {
            my_histos["h_njets_2l"]->Fill(NGoodJets_pt30, weight_2l);
        }
       
        // All baseline selection variables
        // General Selections (all channels)
        bool pass_noniso = NNonIsoMuons == 0;
        bool pass_ht = HT_trigger_pt30 > 500;

        // 0l specific cuts
        bool pass_0l = NGoodLeptons == 0;
        bool pass_njets_0l_pt30 = NGoodJets_pt30 >= 8;
        bool pass_njets_0l_pt45 = NGoodJets_pt45 >= 6;
        bool pass_nbjets_0l_pt30 = NGoodBJets_pt30 >= 2;
        bool pass_nbjets_0l_pt45 = NGoodBJets_pt45 >= 1;
        bool pass_ntops = ntops >= 2;
        bool pass_dr_bjets = dR_bjets >= 1.0;
      
        my_histos["h_cutflow_0l"]->GetXaxis()->SetBinLabel(1,"Inclusive");
        my_histos["h_cutflow_0l"]->GetXaxis()->SetBinLabel(2,"Trigger");
        my_histos["h_cutflow_0l"]->GetXaxis()->SetBinLabel(3,"JetID");
        my_histos["h_cutflow_0l"]->GetXaxis()->SetBinLabel(4,"E_{T}^{Miss} Filter");
        my_histos["h_cutflow_0l"]->GetXaxis()->SetBinLabel(5,"HEM Electron Veto");
        my_histos["h_cutflow_0l"]->GetXaxis()->SetBinLabel(6,"H_{T} > 500 GeV");
        my_histos["h_cutflow_0l"]->GetXaxis()->SetBinLabel(7,"N_{Leptons} = 0");
        my_histos["h_cutflow_0l"]->GetXaxis()->SetBinLabel(8,"N_{Jets} #geq 6 (p_{T} > 45 GeV)");
        my_histos["h_cutflow_0l"]->GetXaxis()->SetBinLabel(9,"N_{Jets} #geq 8 (p_{T} > 30 GeV");
        my_histos["h_cutflow_0l"]->GetXaxis()->SetBinLabel(10,"N_b #geq 1 (p_{T} > 45 GeV)");
        my_histos["h_cutflow_0l"]->GetXaxis()->SetBinLabel(11,"N_b #geq 2 (p_{T} > 30 GeV)");
        my_histos["h_cutflow_0l"]->GetXaxis()->SetBinLabel(12,"N_{t} #geq 2");
        my_histos["h_cutflow_0l"]->GetXaxis()->SetBinLabel(13,"#Delta R_{b} #geq 1.0");
        my_histos["h_cutflow_0l"]->GetXaxis()->SetBinLabel(14,"N_{muon}^{non-iso} = 0");
        my_histos["h_cutflow_0l"]->GetXaxis()->SetBinLabel(15,"Madgraph HT Filter");
 
        my_histos["h_cutflow_0l"]->Fill(1, weight_0l);
        if (passTrigger && passTriggerHadMC) {
            my_histos["h_cutflow_0l"]->Fill(2, weight_0l);
            if (JetID) {
                my_histos["h_cutflow_0l"]->Fill(3, weight_0l);
                if (passMETFilters) {
                    my_histos["h_cutflow_0l"]->Fill(4, weight_0l);
                    if (passElectronHEMveto) {
                        my_histos["h_cutflow_0l"]->Fill(5, weight_0l);
                        if (pass_0l) {        
                            my_histos["h_cutflow_0l"]->Fill(6, weight_0l);
                            if (pass_ht) {
                                my_histos["h_cutflow_0l"]->Fill(7, weight_0l);
                                if (pass_njets_0l_pt45) {
                                    my_histos["h_cutflow_0l"]->Fill(8, weight_0l);
                                    if (pass_njets_0l_pt30) {
                                        my_histos["h_cutflow_0l"]->Fill(9, weight_0l);
                                        if (pass_nbjets_0l_pt45) {
                                            my_histos["h_cutflow_0l"]->Fill(10, weight_0l);
                                            if (pass_nbjets_0l_pt30) {
                                                my_histos["h_cutflow_0l"]->Fill(11, weight_0l);
                                                if (pass_ntops) {
                                                    my_histos["h_cutflow_0l"]->Fill(12, weight_0l);
                                                    if (pass_dr_bjets) {
                                                        my_histos["h_cutflow_0l"]->Fill(13, weight_0l);
                                                        if (pass_noniso) {
                                                            my_histos["h_cutflow_0l"]->Fill(14, weight_0l);
                                                            if (passMadHT) {
                                                                my_histos["h_cutflow_0l"]->Fill(15, weight_0l);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // 1l specific cuts
        bool pass_1l = NGoodLeptons == 1;
        bool pass_njets_1l_pt30 = NGoodJets_pt30 >= 7;
        bool pass_nbjets_1l_pt30 = NGoodBJets_pt30 >= 1;
        bool pass_mbl = 50 < Mbl && Mbl < 250;

        my_histos["h_cutflow_1l"]->GetXaxis()->SetBinLabel(1,"Inclusive");
        my_histos["h_cutflow_1l"]->GetXaxis()->SetBinLabel(2,"Trigger");
        my_histos["h_cutflow_1l"]->GetXaxis()->SetBinLabel(3,"JetID");
        my_histos["h_cutflow_1l"]->GetXaxis()->SetBinLabel(4,"E_{T}^{Miss} Filter");
        my_histos["h_cutflow_1l"]->GetXaxis()->SetBinLabel(5,"HEM Electron Veto");
        my_histos["h_cutflow_1l"]->GetXaxis()->SetBinLabel(6,"H_{T} > 500 GeV");
        my_histos["h_cutflow_1l"]->GetXaxis()->SetBinLabel(7,"N_{Leptons} = 1");
        my_histos["h_cutflow_1l"]->GetXaxis()->SetBinLabel(8,"N_{Jets} #geq 7 (p_{T} > 30 GeV");
        my_histos["h_cutflow_1l"]->GetXaxis()->SetBinLabel(9,"N_b #geq 1 (p_{T} > 30 GeV)");
        my_histos["h_cutflow_1l"]->GetXaxis()->SetBinLabel(10,"N_{muon}^{non-iso} = 0");
        my_histos["h_cutflow_1l"]->GetXaxis()->SetBinLabel(11,"Madgraph HT Filter");
        my_histos["h_cutflow_1l"]->GetXaxis()->SetBinLabel(12,"M_{bl}");

        my_histos["h_cutflow_1l"]->Fill(1, weight_1l);
        if (passTrigger && passTriggerMC) {
            my_histos["h_cutflow_1l"]->Fill(2, weight_1l);
            if (JetID) {
                my_histos["h_cutflow_1l"]->Fill(3, weight_1l);
                if (passMETFilters) {
                    my_histos["h_cutflow_1l"]->Fill(4, weight_1l);
                    if (passElectronHEMveto) {
                        my_histos["h_cutflow_1l"]->Fill(5, weight_1l);
                        if (pass_1l) {        
                            my_histos["h_cutflow_1l"]->Fill(6, weight_1l);
                            if (pass_ht) {
                                my_histos["h_cutflow_1l"]->Fill(7, weight_1l);
                                if (pass_njets_1l_pt30) {
                                    my_histos["h_cutflow_1l"]->Fill(8, weight_1l);
                                    if (pass_nbjets_1l_pt30) {
                                        my_histos["h_cutflow_1l"]->Fill(9, weight_1l);
                                        if (pass_noniso) {
                                            my_histos["h_cutflow_1l"]->Fill(10, weight_1l);
                                            if (passMadHT) {
                                                my_histos["h_cutflow_1l"]->Fill(11, weight_1l);
                                                if (pass_mbl) {
                                                    my_histos["h_cutflow_1l"]->Fill(12, weight_1l);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // 2l specific cuts
        bool pass_2l = NGoodLeptons == 2;
        bool pass_njets_2l_pt30 = NGoodJets_pt30 >= 6;
        bool pass_nbjets_2l_pt30 = NGoodBJets_pt30 >= 1;
        bool pass_opp_charge = false;

        my_histos["h_cutflow_2l"]->GetXaxis()->SetBinLabel(1,"Inclusive");
        my_histos["h_cutflow_2l"]->GetXaxis()->SetBinLabel(2,"Trigger");
        my_histos["h_cutflow_2l"]->GetXaxis()->SetBinLabel(3,"JetID");
        my_histos["h_cutflow_2l"]->GetXaxis()->SetBinLabel(4,"E_{T}^{Miss} Filter");
        my_histos["h_cutflow_2l"]->GetXaxis()->SetBinLabel(5,"HEM Electron Veto");
        my_histos["h_cutflow_2l"]->GetXaxis()->SetBinLabel(6,"H_{T} > 500 GeV");
        my_histos["h_cutflow_2l"]->GetXaxis()->SetBinLabel(7,"N_{Leptons} = 1");
        my_histos["h_cutflow_2l"]->GetXaxis()->SetBinLabel(8,"N_{Jets} #geq 7 (p_{T} > 30 GeV");
        my_histos["h_cutflow_2l"]->GetXaxis()->SetBinLabel(9,"N_b #geq 1 (p_{T} > 30 GeV)");
        my_histos["h_cutflow_2l"]->GetXaxis()->SetBinLabel(10,"N_{muon}^{non-iso} = 0");
        my_histos["h_cutflow_2l"]->GetXaxis()->SetBinLabel(11,"Madgraph HT Filter");
        my_histos["h_cutflow_2l"]->GetXaxis()->SetBinLabel(12,"Off Z Resonance");
        my_histos["h_cutflow_2l"]->GetXaxis()->SetBinLabel(13,"Opposite Charge Requirement");

        my_histos["h_cutflow_2l"]->Fill(1, weight_2l);
        if (passTrigger && passTriggerMC) {
            my_histos["h_cutflow_2l"]->Fill(2, weight_2l);
            if (JetID) {
                my_histos["h_cutflow_2l"]->Fill(3, weight_2l);
                if (passMETFilters) {
                    my_histos["h_cutflow_2l"]->Fill(4, weight_2l);
                    if (passElectronHEMveto) {
                        my_histos["h_cutflow_2l"]->Fill(5, weight_2l);
                        if (pass_2l) {        
                            pass_opp_charge = GoodLeptonsCharge[0]!=GoodLeptonsCharge[1];
                            my_histos["h_cutflow_2l"]->Fill(6, weight_2l);
                            if (pass_ht) {
                                my_histos["h_cutflow_2l"]->Fill(7, weight_2l);
                                if (pass_njets_2l_pt30) {
                                    my_histos["h_cutflow_2l"]->Fill(8, weight_2l);
                                    if (pass_nbjets_2l_pt30) {
                                        my_histos["h_cutflow_2l"]->Fill(9, weight_2l);
                                        if (pass_noniso) {
                                            my_histos["h_cutflow_2l"]->Fill(10, weight_2l);
                                            if (passMadHT) {
                                                my_histos["h_cutflow_2l"]->Fill(11, weight_2l);
                                                if (!onZ) {
                                                    my_histos["h_cutflow_2l"]->Fill(12, weight_2l);
                                                    if (pass_opp_charge) {
                                                        my_histos["h_cutflow_2l"]->Fill(13, weight_2l);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

    } 
}

void MakeCutflow::WriteHistos(TFile* outfile)
{
    outfile->cd();

    for (const auto &p : my_histos) {
        p.second->Write();
    }
    
}
