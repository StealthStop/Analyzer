#define Analyze0Lep_cxx
#include "Analyzer/Analyzer/include/Analyze0Lep.h"
#include "SusyAnaTools/Tools/NTupleReader.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TRandom3.h>
#include <iostream>

Analyze0Lep::Analyze0Lep() : initHistos(false)
{
}

void Analyze0Lep::InitHistos(const std::map<std::string, bool>& cutMap)
{
    TH1::SetDefaultSumw2();

    // Declare all your histograms here, that way we can fill them for multiple chains
    my_histos.emplace("h_met",      std::make_shared<TH1D>("h_met",     "h_met",      20,  0,    200  ) );
    my_histos.emplace("h_ht",       std::make_shared<TH1D>("h_ht",      "h_ht",       60,  0,   3000  ) );
    my_histos.emplace("h_njets",    std::make_shared<TH1D>("h_njets",   "h_njets",    20,  0,     20  ) );
    my_histos.emplace("h_nb",       std::make_shared<TH1D>("h_nb",      "h_nb",       10,  0,     10  ) );
    my_histos.emplace("h_ntops",    std::make_shared<TH1D>("h_ntops",   "h_ntops",    10,  0,     10  ) );
    my_histos.emplace("h_ntops_j1", std::make_shared<TH1D>("h_ntops_j1","h_ntops_j1", 10,  0,     10  ) );
    my_histos.emplace("h_ntops_j2", std::make_shared<TH1D>("h_ntops_j2","h_ntops_j2", 10,  0,     10  ) );
    my_histos.emplace("h_ntops_j3", std::make_shared<TH1D>("h_ntops_j3","h_ntops_j3", 10,  0,     10  ) );

    my_histos.emplace("blind_met",      std::make_shared<TH1D>("blind_met",     "blind_met",      20,  0,    200  ) );
    my_histos.emplace("blind_ht",       std::make_shared<TH1D>("blind_ht",      "blind_ht",       60,  0,   3000  ) );
    my_histos.emplace("blind_njets",    std::make_shared<TH1D>("blind_njets",   "blind_njets",    20,  0,     20  ) );
    my_histos.emplace("blind_nb",       std::make_shared<TH1D>("blind_nb",      "blind_nb",       10,  0,     10  ) );
    my_histos.emplace("blind_ntops",    std::make_shared<TH1D>("blind_ntops",   "blind_ntops",    10,  0,     10  ) );
    my_histos.emplace("blind_ntops_j1", std::make_shared<TH1D>("blind_ntops_j1","blind_ntops_j1", 10,  0,     10  ) );
    my_histos.emplace("blind_ntops_j2", std::make_shared<TH1D>("blind_ntops_j2","blind_ntops_j2", 10,  0,     10  ) );
    my_histos.emplace("blind_ntops_j3", std::make_shared<TH1D>("blind_ntops_j3","blind_ntops_j3", 10,  0,     10  ) );

    for(auto& mycut : cutMap)
    {
        my_histos.emplace("h_njets_0l_"+mycut.first, std::make_shared<TH1D>(("h_njets_0l_"+mycut.first).c_str(),("h_njets_0l_"+mycut.first).c_str(), 20, 0, 20));
        my_histos.emplace("h_ntops_0l_"+mycut.first, std::make_shared<TH1D>(("h_ntops_0l_"+mycut.first).c_str(),("h_ntops_0l_"+mycut.first).c_str(), 10, 0, 10));
        my_histos.emplace("h_nb_0l_"+mycut.first,    std::make_shared<TH1D>(("h_nb_0l_"+mycut.first).c_str(),("h_nb_0l_"+mycut.first).c_str(), 10, 0, 10));
        my_histos.emplace("h_HT_0l_"+mycut.first,    std::make_shared<TH1D>(("h_HT_0l_"+mycut.first).c_str(),("h_HT_0l_"+mycut.first).c_str(), 60, 0, 3000));
    

        my_histos.emplace("blind_njets_0l_"+mycut.first, std::make_shared<TH1D>(("blind_njets_0l_"+mycut.first).c_str(),("blind_njets_0l_"+mycut.first).c_str(), 20, 0, 20));
        my_histos.emplace("blind_ntops_0l_"+mycut.first, std::make_shared<TH1D>(("blind_ntops_0l_"+mycut.first).c_str(),("blind_ntops_0l_"+mycut.first).c_str(), 10, 0, 10));
        my_histos.emplace("blind_nb_0l_"+mycut.first,    std::make_shared<TH1D>(("blind_nb_0l_"+mycut.first).c_str(),("blind_nb_0l_"+mycut.first).c_str(), 10, 0, 10));
        my_histos.emplace("blind_HT_0l_"+mycut.first,    std::make_shared<TH1D>(("blind_HT_0l_"+mycut.first).c_str(),("blind_HT_0l_"+mycut.first).c_str(), 60, 0, 3000));
    
    }

    // Cut flows
    my_efficiencies.emplace("event_sel",       std::make_shared<TEfficiency>("event_sel","Event selection efficiency wrt previous cut;Cut;#epsilon",8,0,8));
    my_efficiencies.emplace("event_sel_total", std::make_shared<TEfficiency>("event_sel_total","Total event selection efficiency;Cut;#epsilon",8,0,8));

    my_efficiencies.emplace("event_sel_0l",       std::make_shared<TEfficiency>("event_sel_0l","0 lepton event selection efficiency wrt previous cut;Cut;#epsilon",8,0,8));
    my_efficiencies.emplace("event_sel_total_0l", std::make_shared<TEfficiency>("event_sel_total_0l","Total 0 lepton event selection efficiency;Cut;#epsilon",8,0,8));

    my_efficiencies.emplace("event_sel_1l",       std::make_shared<TEfficiency>("event_sel_1l","1 lepton event selection efficiency wrt previous cut;Cut;#epsilon",8,0,8));
    my_efficiencies.emplace("event_sel_total_1l", std::make_shared<TEfficiency>("event_sel_total_1l","Total 1 lepton event selection efficiency;Cut;#epsilon",8,0,8));

    my_efficiencies.emplace("event_sel_2l",       std::make_shared<TEfficiency>("event_sel_2l","2 lepton event selection efficiency wrt previous cut;Cut;#epsilon",8,0,8));
    my_efficiencies.emplace("event_sel_total_2l", std::make_shared<TEfficiency>("event_sel_total_2l","Total 2 lepton event selection efficiency;Cut;#epsilon",8,0,8));

}

void Analyze0Lep::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    while( tr.getNextEvent() )
    {
        const auto& MET                = tr.getVar<double>("MET");
        const auto& HT_trigger         = tr.getVar<double>("HT_trigger");
        const auto& ntops              = tr.getVar<int>("ntops");
        const auto& ntops_3jet         = tr.getVar<int>("ntops_3jet");
        const auto& ntops_2jet         = tr.getVar<int>("ntops_2jet");
        const auto& ntops_1jet         = tr.getVar<int>("ntops_1jet");
        const auto& runtype            = tr.getVar<std::string>("runtype");     
        const auto& NJets_pt30         = tr.getVar<int>("NGoodJets_pt30");
        const auto& NJets_pt45         = tr.getVar<int>("NGoodJets_pt45");
        const auto& NBJets             = tr.getVar<int>("NGoodBJets");
        const auto& NBJets_pt45        = tr.getVar<int>("NGoodBJets_pt45");
        const auto& NGoodLeptons       = tr.getVar<int>("NGoodLeptons");
        const auto& passBaseline0l     = tr.getVar<bool>("passBaseline0l_Good");
        const auto& passBlindHad       = tr.getVar<bool>("passBlindHad_Good");
        const auto& passTrigger        = tr.getVar<bool>("passTrigger");
        const auto& passMadHT          = tr.getVar<bool>("passMadHT");

        // ------------------------
        // -- Print event number
        // -----------------------        

        if(maxevents != -1 && tr.getEvtNum() >= maxevents) break;        
        if ( tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() ) ;

        // ------------------------
        // -- Define weight
        // -----------------------

        double eventweight = 1.0;        
        // Weight from samples.cc
        //eventweight = weight;

        if(runtype == "MC"){
            const auto& Weight  = tr.getVar<double>("Weight");
            double lumi = 35900; // Lumi for 2016
            // Weight from NTuples            
            eventweight = lumi*Weight;
        }

        // -------------------------------
        // -- Define cuts
        // -------------------------------

        // Global cuts
        if ( !(passBlindHad && passTrigger && passMadHT) ) continue;
        
        bool pass_0l              = NGoodLeptons==0;
        bool pass_njet_pt45       = NJets_pt45>=6;
        bool pass_HT_trigger      = HT_trigger > 500;
        bool pass_njet_pt45_1btag = NBJets_pt45 >= 1;
        bool pass_njet_pt45_2btag = NBJets_pt45 >= 2;

        // ---------------------------
        // --    1 Top selection
        // ---------------------------
        //exactly selections
        bool pass_1t    = ntops==1;
        bool pass_1t1 = pass_1t && ntops_1jet==1;
        bool pass_1t2 = pass_1t && ntops_2jet==1;
        bool pass_1t3 = pass_1t && ntops_3jet==1;
        bool pass_1t2or3 = pass_1t2 or pass_1t3;

        //at least selections
        bool pass_ge1t    = ntops>=1;

        bool pass_ge1t1 = pass_ge1t && ntops_1jet>=1;
        bool pass_ge1t2 = pass_ge1t && ntops_2jet>=1;
        bool pass_ge1t3 = pass_ge1t && ntops_3jet>=1;

        // ---------------------------
        // --    2 Top selection
        // ---------------------------

        // exactly selections
        bool pass_2t    = ntops==2;

        bool pass_2t11 = pass_2t && ntops_1jet==2;
        bool pass_2t12 = pass_2t && ntops_1jet==1 && ntops_2jet==1;
        bool pass_2t13 = pass_2t && ntops_1jet==1 && ntops_3jet==1;
        bool pass_2t22 = pass_2t && ntops_2jet==2;
        bool pass_2t23 = pass_2t && ntops_2jet==1 && ntops_3jet==1;
        bool pass_2t33 = pass_2t && ntops_3jet==2;

        // at least selections
        bool pass_ge2t    = ntops>=2;

        bool pass_ge2t11 = pass_ge2t && ntops_1jet>=2;
        bool pass_ge2t12 = pass_ge2t && ntops_1jet>=1 && ntops_2jet>=1;
        bool pass_ge2t13 = pass_ge2t && ntops_1jet>=1 && ntops_3jet>=1;
        bool pass_ge2t22 = pass_ge2t && ntops_2jet>=2;
        bool pass_ge2t23 = pass_ge2t && ntops_2jet>=1 && ntops_3jet>=1;
        bool pass_ge2t33 = pass_ge2t && ntops_3jet>=2;

        //Merge or no Merge selections
        bool pass_Merge = ntops_1jet >= 1, pass_noMerge = ntops_1jet == 0;
        bool pass_ge2t11or12or13 = pass_Merge   && (pass_ge2t11 or pass_ge2t12 or pass_ge2t13);
        bool pass_ge2t22or23or33 = pass_noMerge && (pass_ge2t22 or pass_ge2t23 or pass_ge2t33);

        // -------------------
        // --- Fill Histos ---
        // -------------------                        
        const std::map<std::string, bool> cut_map_0l 
        {
            {""                                  , pass_0l                                                             },
            {"ge6j"                              , pass_0l && pass_njet_pt45                                           },
            {"HT500"                             , pass_0l && pass_HT_trigger                                          },
            {"ge2b"                              , pass_0l && pass_njet_pt45_2btag                                     },
            {"1t"                                , pass_0l && pass_1t                                                  },
            {"2t"                                , pass_0l && pass_2t                                                  },
            {"ge1t"                              , pass_0l && pass_ge1t                                                },
            {"ge2t"                              , pass_0l && pass_ge2t                                                },
            {"ge6j_HT500"                        , pass_0l && pass_njet_pt45 && pass_HT_trigger                        },
            {"ge6j_HT500_ge1b"                   , pass_0l && pass_njet_pt45 && pass_HT_trigger && pass_njet_pt45_1btag},
            {"ge6j_HT500_ge2b"                   , passBaseline0l                                                      },                         
            {"6j_HT500_ge2b"                     , passBaseline0l && NJets_pt45 == 6                                   },
            {"7j_HT500_ge2b"                     , passBaseline0l && NJets_pt45 == 7                                   },
            {"8j_HT500_ge2b"                     , passBaseline0l && NJets_pt45 == 8                                   },
            {"9j_HT500_ge2b"                     , passBaseline0l && NJets_pt45 == 9                                   },
            {"10j_HT500_ge2b"                    , passBaseline0l && NJets_pt45 == 10                                  },
            {"11j_HT500_ge2b"                    , passBaseline0l && NJets_pt45 == 11                                  },
            {"12j_HT500_ge2b"                    , passBaseline0l && NJets_pt45 == 12                                  },
            {"13j_HT500_ge2b"                    , passBaseline0l && NJets_pt45 == 13                                  },
            {"14j_HT500_ge2b"                    , passBaseline0l && NJets_pt45 == 14                                  },
            {"15j_HT500_ge2b"                    , passBaseline0l && NJets_pt45 == 15                                  },
                                                 
            {"ge6j_HT500_ge2b_1t"                , passBaseline0l && pass_1t                                           },
            {"ge6j_HT500_ge2b_ge1t"              , passBaseline0l && pass_ge1t                                         },
                                                 
            {"ge6j_HT500_ge2b_1t1"               , passBaseline0l && pass_1t1                                          },
            {"ge6j_HT500_ge2b_1t2"               , passBaseline0l && pass_1t2                                          },
            {"ge6j_HT500_ge2b_1t3"               , passBaseline0l && pass_1t3                                          },
            {"ge6j_HT500_ge2b_1t2or3"            , passBaseline0l && pass_1t2or3                                       },
            {"ge6j_HT500_ge2b_ge1t1"             , passBaseline0l && pass_ge1t1                                        },
            {"ge6j_HT500_ge2b_ge1t2"             , passBaseline0l && pass_ge1t2                                        },
            {"ge6j_HT500_ge2b_ge1t3"             , passBaseline0l && pass_ge1t3                                        },
                                                 
            {"ge6j_HT500_ge2b_2t"                , passBaseline0l && pass_2t                                           },
            {"ge6j_HT500_ge2b_ge2t"              , passBaseline0l && pass_ge2t                                         },
                                                 
            {"ge6j_HT500_ge2b_2t11"              , passBaseline0l && pass_2t11                                         },
            {"ge6j_HT500_ge2b_2t12"              , passBaseline0l && pass_2t12                                         },
            {"ge6j_HT500_ge2b_2t13"              , passBaseline0l && pass_2t13                                         },
            {"ge6j_HT500_ge2b_2t22"              , passBaseline0l && pass_2t22                                         },
            {"ge6j_HT500_ge2b_2t23"              , passBaseline0l && pass_2t23                                         },
            {"ge6j_HT500_ge2b_2t33"              , passBaseline0l && pass_2t33                                         },
            {"ge6j_HT500_ge2b_ge2t11"            , passBaseline0l && pass_ge2t11                                       },
            {"ge6j_HT500_ge2b_ge2t12"            , passBaseline0l && pass_ge2t12                                       },
            {"ge6j_HT500_ge2b_ge2t13"            , passBaseline0l && pass_ge2t13                                       },
            {"ge6j_HT500_ge2b_ge2t11or12or13"    , passBaseline0l && pass_ge2t11or12or13                               },
            {"ge6j_HT500_ge2b_ge2t22"            , passBaseline0l && pass_ge2t22                                       },
            {"ge6j_HT500_ge2b_ge2t23"            , passBaseline0l && pass_ge2t23                                       },
            {"ge6j_HT500_ge2b_ge2t33"            , passBaseline0l && pass_ge2t33                                       },
            {"ge6j_HT500_ge2b_ge2t22or23or33"    , passBaseline0l && pass_ge2t22or23or33                               },
        };

        // Initialize Histograms
        if(!initHistos)
        {
            InitHistos(cut_map_0l);
            initHistos = true;
        }

        for(auto& kv : cut_map_0l)
        {
            if(kv.second)
            {
                my_histos["h_njets_0l_"   +kv.first]->Fill(NJets_pt30, eventweight);
                my_histos["h_ntops_0l_"   +kv.first]->Fill(ntops, eventweight);
                my_histos["h_nb_0l_"      +kv.first]->Fill(NBJets, eventweight);
                my_histos["h_HT_0l_"      +kv.first]->Fill(HT_trigger, eventweight);

                if ( NJets_pt30 < 9 )
                {
                    my_histos["blind_njets_0l_"   +kv.first]->Fill(NJets_pt30, eventweight);
                    my_histos["blind_ntops_0l_"   +kv.first]->Fill(ntops, eventweight);
                    my_histos["blind_nb_0l_"      +kv.first]->Fill(NBJets, eventweight);
                    my_histos["blind_HT_0l_"      +kv.first]->Fill(HT_trigger, eventweight);
                }
            }
        }

        // Fill event selection efficiencies
        my_efficiencies["event_sel_total"]->Fill(true,0);
        my_efficiencies["event_sel_total"]->Fill(HT_trigger>500,1);
        my_efficiencies["event_sel_total"]->Fill(HT_trigger>500 && NJets_pt45>=6 ,2);
        my_efficiencies["event_sel_total"]->Fill(HT_trigger>500 && NJets_pt45>=6 && NBJets_pt45>0 ,3);
        my_efficiencies["event_sel_total"]->Fill(HT_trigger>500 && NJets_pt45>=6 && NBJets_pt45>0 && ntops>0 ,4);
        my_efficiencies["event_sel_total"]->Fill(HT_trigger>500 && NJets_pt45>=6 && NBJets_pt45>0 && ntops>0 && NBJets_pt45>1 ,5);
        my_efficiencies["event_sel_total"]->Fill(HT_trigger>500 && NJets_pt45>=6 && NBJets_pt45>0 && ntops>0 && NBJets_pt45>1 && ntops>1 ,6);
        my_efficiencies["event_sel_total"]->Fill(HT_trigger>500 && NJets_pt45>=6 && NBJets_pt45>0 && ntops>0 && NBJets_pt45>1 && ntops>1 && NJets_pt30>=8 ,7);
        
        my_efficiencies["event_sel"]->Fill(true,0);
        my_efficiencies["event_sel"]->Fill(HT_trigger>500,1);
        if(HT_trigger>500)
        {
            my_efficiencies["event_sel"]->Fill(NJets_pt45>=6,2);
            if (NJets_pt45>=6)
            {
                my_efficiencies["event_sel"]->Fill(NBJets_pt45>0,3);
                if (NBJets_pt45>0)
                {
                    my_efficiencies["event_sel"]->Fill(ntops>0,4);
                    if (ntops>0)
                    {
                        my_efficiencies["event_sel"]->Fill(NBJets_pt45>1,5);
                        if (NBJets_pt45>1)
                        {
                            my_efficiencies["event_sel"]->Fill(ntops>1,6);
                            if (ntops>1)
                            {
                                my_efficiencies["event_sel"]->Fill(NJets_pt30>=8,7);
                            }
                        }
                    }
                }
            }
        }
        
        // No local cuts applied here
        my_histos["h_met"     ]->Fill(MET, eventweight);
        my_histos["h_ht"      ]->Fill(HT_trigger, eventweight);
        my_histos["h_njets"   ]->Fill(NJets_pt30, eventweight);
        my_histos["h_nb"      ]->Fill(NBJets, eventweight);
        my_histos["h_ntops"   ]->Fill(ntops, eventweight);        
        my_histos["h_ntops_j1"]->Fill(ntops_1jet, eventweight);
        my_histos["h_ntops_j2"]->Fill(ntops_2jet, eventweight);
        my_histos["h_ntops_j3"]->Fill(ntops_3jet, eventweight);        

        // Do the blind cuts, for now, for data MC plots
        if ( NJets_pt30 < 9 )
        {
            my_histos["blind_met"     ]->Fill(MET, eventweight);
            my_histos["blind_ht"      ]->Fill(HT_trigger, eventweight);
            my_histos["blind_njets"   ]->Fill(NJets_pt30, eventweight);
            my_histos["blind_nb"      ]->Fill(NBJets, eventweight);
            my_histos["blind_ntops"   ]->Fill(ntops, eventweight);        
            my_histos["blind_ntops_j1"]->Fill(ntops_1jet, eventweight);
            my_histos["blind_ntops_j2"]->Fill(ntops_2jet, eventweight);
            my_histos["blind_ntops_j3"]->Fill(ntops_3jet, eventweight);        
        }

    } // end of event loop

}

void Analyze0Lep::WriteHistos(TFile* outfile)
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
