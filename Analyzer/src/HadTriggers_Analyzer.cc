#define HadTriggers_Analyzer_cxx
#include "Analyzer/Analyzer/include/HadTriggers_Analyzer.h"
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

HadTriggers_Analyzer::HadTriggers_Analyzer()
{
    InitHistos();
}

void HadTriggers_Analyzer::InitHistos()
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ;

    // ----------------
    // for Our Analysis
    // ----------------
    std::vector<std::string> effTags   { "denominator", "numerator" }; 
    std::vector<std::string> hadTags   { "had", "had_IsoMu"         };
    std::vector<std::string> trigTags  { "trig", "noTrig"           };
    std::vector<std::string> ptTags    { "pt45"                     }; // label for GoodJets_pt45 & GoodBJets_pt45 & HT_trigger_pt45 
 
    const int nHTbins   = 7;
    const int nhtBins   = 13;
    const int nJetBins  = 9;
    const int nBJetBins = 5;
    double HTbinEdges[nHTbins + 1]      = {400, 500, 600, 700, 800, 900, 1000, 1100};
    double htBinEdges[nhtBins + 1 ]     = {0, 200, 400, 600, 800, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000};
    double njetBinEdges[nJetBins + 1]   = {6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    double nbjetBinEdges[nBJetBins + 1] = {0, 1, 2, 3, 4, 5};

    for( std::string effTag : effTags ) 
    {
        for( std::string hadTag : hadTags ) 
        {
            for( std::string trigTag : trigTags )
            {
                    for( std::string ptTag : ptTags ) 
                    {
                        my_histos.emplace( "h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_HT", std::make_shared<TH1D>(("h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_HT").c_str(), ("h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_HT").c_str(), nHTbins, HTbinEdges ) );
                        my_histos.emplace( "h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_ht5000", std::make_shared<TH1D>(("h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_ht5000").c_str(), ("h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_ht5000").c_str(), nhtBins, htBinEdges ) );
                        my_histos.emplace( "h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_NJet", std::make_shared<TH1D>(("h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_NJet").c_str(), ("h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_NJet").c_str(), nJetBins, njetBinEdges ) );
                        my_histos.emplace( "h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_NBJet", std::make_shared<TH1D>(("h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_NBJet").c_str(), ("h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_NBJet").c_str(), nBJetBins, nbjetBinEdges ) );                                

                        my_2d_histos.emplace( "h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_NJetVsHT", std::make_shared<TH2D>( ( "h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_NJetVsHT" ).c_str(), ( "h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_NJetVsHT" ).c_str(), nJetBins, njetBinEdges, nHTbins, HTbinEdges ) );
                        my_2d_histos.emplace( "h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_NJetVsHt", std::make_shared<TH2D>( ( "h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_NJetVsHt" ).c_str(), ( "h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_NJetVsHt" ).c_str(), nJetBins, njetBinEdges, nhtBins, htBinEdges ) );
                        my_2d_histos.emplace( "h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_NJetVsNBJet", std::make_shared<TH2D>( ( "h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_NJetVsNBJet" ).c_str(), ( "h_"+effTag+"_"+hadTag+"_"+trigTag+"_"+ptTag+"_NJetVsNBJet" ).c_str(), nJetBins, njetBinEdges, nBJetBins, nbjetBinEdges ) );

                }
            }
        }
    }

    // ------------------------------------
    // for Reference Analysis - AN-2016/411
    // ------------------------------------
    std::vector<std::string> refANTags    { "refAnHad", "refAnHadIsoMu" };
    std::vector<std::string> nBJetCutTags { "ge2bjetCut", "2bjetCut", "3bjetCut", "ge4bjetCut" };
    const int htbins   = 9;
    const int ptbins   = 7;
    const int bjetbins = 3;
    double htbinEdges[htbins + 1]     = {500, 550, 600, 650, 700, 800, 1000, 1500, 2000, 2500};
    double ptbinEdges[ptbins + 1 ]    = {40, 45, 50, 55, 60, 70, 120, 200};
    double bjetbinEdges[bjetbins + 1] = {1.5, 2.5, 3.5, 8};   

    for( std::string effTag : effTags )
    {
        for( std::string refANTag : refANTags )
        {
            for( std::string trigTag : trigTags ) 
            {
                for( std::string nBJetCutTag : nBJetCutTags )
                {
                    my_histos.emplace( "h_"+effTag+"_"+refANTag+"_"+trigTag+"_"+nBJetCutTag+"_HT", std::make_shared<TH1D>(("h_"+effTag+"_"+refANTag+"_"+trigTag+"_"+nBJetCutTag+"_HT").c_str(), ("h_"+effTag+"_"+refANTag+"_"+trigTag+"_"+nBJetCutTag+"_HT").c_str(), htbins, htbinEdges ) );
                    my_histos.emplace( "h_"+effTag+"_"+refANTag+"_"+trigTag+"_"+nBJetCutTag+"_6thJetPt", std::make_shared<TH1D>(("h_"+effTag+"_"+refANTag+"_"+trigTag+"_"+nBJetCutTag+"_6thJetPt").c_str(), ("h_"+effTag+"_"+refANTag+"_"+trigTag+"_"+nBJetCutTag+"_6thJetPt").c_str(), ptbins, ptbinEdges ) );
                    my_histos.emplace( "h_"+effTag+"_"+refANTag+"_"+trigTag+"_"+nBJetCutTag+"_NBJet", std::make_shared<TH1D>(("h_"+effTag+"_"+refANTag+"_"+trigTag+"_"+nBJetCutTag+"_NBJet").c_str(), ("h_"+effTag+"_"+refANTag+"_"+trigTag+"_"+nBJetCutTag+"_NBJet").c_str(), bjetbins, bjetbinEdges ) );

                    my_2d_histos.emplace( "h_"+effTag+"_"+refANTag+"_"+trigTag+"_"+nBJetCutTag+"_HTvs6thJetPt", std::make_shared<TH2D>( ( "h_"+effTag+"_"+refANTag+"_"+trigTag+"_"+nBJetCutTag+"_HTvs6thJetPt" ).c_str(), ( "h_"+effTag+"_"+refANTag+"_"+trigTag+"_"+nBJetCutTag+"_HTvs6thJetPt" ).c_str(), htbins, htbinEdges, ptbins, ptbinEdges ) );

                }        
            }
        }
    }
   
    // --------------------------------------------------
    // their triggers + their preselections with our pt45
    // --------------------------------------------------
    std::vector<std::string> theirTrigTags { "TheirHadIsoMu" };

    for( std::string effTag : effTags )
    {
        for( std::string theirTrigTag : theirTrigTags )
        {
            for( std::string trigTag : trigTags )
            {
                for( std::string nBJetCutTag : nBJetCutTags )
                {
                    my_histos.emplace( "h_"+effTag+"_"+theirTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_HT", std::make_shared<TH1D>(("h_"+effTag+"_"+theirTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_HT").c_str(), ("h_"+effTag+"_"+theirTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_HT").c_str(), htbins, htbinEdges ) );
                    my_histos.emplace( "h_"+effTag+"_"+theirTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_6thJetPt", std::make_shared<TH1D>(("h_"+effTag+"_"+theirTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_6thJetPt").c_str(), ("h_"+effTag+"_"+theirTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_6thJetPt").c_str(), ptbins, ptbinEdges ) );
                    my_histos.emplace( "h_"+effTag+"_"+theirTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_NBJet", std::make_shared<TH1D>(("h_"+effTag+"_"+theirTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_NBJet").c_str(), ("h_"+effTag+"_"+theirTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_NBJet").c_str(), bjetbins, bjetbinEdges ) );

                    my_2d_histos.emplace( "h_"+effTag+"_"+theirTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_HTvs6thJetPt", std::make_shared<TH2D>( ( "h_"+effTag+"_"+theirTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_HTvs6thJetPt" ).c_str(), ( "h_"+effTag+"_"+theirTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_HTvs6thJetPt" ).c_str(), htbins, htbinEdges, ptbins, ptbinEdges ) );

                }
            }
        }
    }
 
    // -----------------------------------------------
    // our triggers + reference analysis preselections
    // -----------------------------------------------
    std::vector<std::string> ourTrigTags { "OurHadIsoMu" };
    
    for( std::string effTag : effTags )
    {
        for( std::string ourTrigTag : ourTrigTags ) 
        {
            for( std::string trigTag : trigTags )
            {
                for( std::string nBJetCutTag : nBJetCutTags )
                {
                    my_histos.emplace( "h_"+effTag+"_"+ourTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_HT", std::make_shared<TH1D>(("h_"+effTag+"_"+ourTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_HT").c_str(), ("h_"+effTag+"_"+ourTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_HT").c_str(), htbins, htbinEdges ) );
                    my_histos.emplace( "h_"+effTag+"_"+ourTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_6thJetPt", std::make_shared<TH1D>(("h_"+effTag+"_"+ourTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_6thJetPt").c_str(), ("h_"+effTag+"_"+ourTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_6thJetPt").c_str(), ptbins, ptbinEdges ) );
                    my_histos.emplace( "h_"+effTag+"_"+ourTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_NBJet", std::make_shared<TH1D>(("h_"+effTag+"_"+ourTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_NBJet").c_str(), ("h_"+effTag+"_"+ourTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_NBJet").c_str(), bjetbins, bjetbinEdges ) );
                    
                    my_2d_histos.emplace( "h_"+effTag+"_"+ourTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_HTvs6thJetPt", std::make_shared<TH2D>( ( "h_"+effTag+"_"+ourTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_HTvs6thJetPt" ).c_str(), ( "h_"+effTag+"_"+ourTrigTag+"_"+trigTag+"_"+nBJetCutTag+"_HTvs6thJetPt" ).c_str(), htbins, htbinEdges, ptbins, ptbinEdges ) );
 
                }          
            }
        }
    }

} //


void HadTriggers_Analyzer::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    while( tr.getNextEvent() )
    {
        const auto& eventCounter        = tr.getVar<int>("eventCounter");
        my_histos["EventCounter"]->Fill( eventCounter );

        const auto& runtype                   = tr.getVar<std::string>("runtype");
        const auto& filetag                   = tr.getVar<std::string>("filetag");
        // for our analysis
        const auto& NGoodJets_pt45            = tr.getVar<int>("NGoodJets_pt45");
        const auto& NGoodBJets_pt45           = tr.getVar<int>("NGoodBJets_pt45");
        const auto& HT_trigger_pt45           = tr.getVar<double>("HT_trigger_pt45");
        const auto& passBaseline0l_hadTrig    = tr.getVar<bool>("passBaseline0l_hadTrig");
        const auto& passBaseline0l_hadMuTrig  = tr.getVar<bool>("passBaseline0l_hadMuTrig");
        const auto& passTriggerAllHad         = tr.getVar<bool>("passTriggerAllHad");
        // for reference analysis AN-2016/411 
        const auto& Jets                      = tr.getVec<TLorentzVector>("Jets");
        const auto& GoodJets_pt40             = tr.getVec<bool>("GoodJets_pt40");
        const auto& NGoodJets_pt40            = tr.getVar<int>("NGoodJets_pt40");
        const auto& NGoodBJets_pt30           = tr.getVar<int>("NGoodBJets_pt30");
        const auto& HT_trigger_pt30           = tr.getVar<double>("HT_trigger_pt30");
        const auto& passBaseline0l_refAN      = tr.getVar<bool>("passBaseline0l_refAN");
        const auto& passTriggerMuonsRefAN     = tr.getVar<bool>("passTriggerMuonsRefAN");
        const auto& passTriggerRefAN          = tr.getVar<bool>("passTriggerRefAN");
        // their triggers + their preselections with our pt45
        const auto& GoodJets_pt45             = tr.getVec<bool>("GoodJets_pt45");
        const auto& passBaseline0l_refAN_pt45 = tr.getVar<bool>("passBaseline0l_refAN_pt45");

        bool pass_2bjetCut   = NGoodBJets_pt30 == 2; 
        bool pass_3bjetCut   = NGoodBJets_pt30 == 3;
        bool pass_ge4bjetCut = NGoodBJets_pt30 >= 4;

        // -----------------------------------------
        // get the 6th jet pt for reference analysis
        // ----------------------------------------- 
        int njets = 0;
        double SixthJetPt = 0.0;

        for (unsigned int j = 0; j < Jets.size(); j++)
        {
            if (!GoodJets_pt40[j]) continue;
            njets++;

            if (njets == 6)
            {
                SixthJetPt = Jets.at(j).Pt();
                break;
            }
        }

        int njetspt45 = 0;
        double SixthJetPt45 = 0.0;
        for (unsigned int j = 0; j < Jets.size(); j++)
        {
            if (!GoodJets_pt45[j]) continue;
            njetspt45++;

            if (njetspt45 == 6)
            {
                SixthJetPt45 = Jets.at(j).Pt();
                break;
            }
        }

        // ------------------------
        // -- Print Event Number 
        // ------------------------
        if(maxevents != -1 && tr.getEvtNum() >= maxevents) break;        
        //if( tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() ) ;

        // ----------------------------
        // -- Print list of triggers 
        // ----------------------------
        //const auto& TriggerNames = tr.getVec<std::string>("TriggerNames");
        //if( tr.getEvtNum() == 1 ) printTriggerList(TriggerNames); 

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
            const auto& Weight   = tr.getVar<double>("Weight");
            const auto& lumi     = tr.getVar<double>("Lumi");
            eventweight          = lumi*Weight;

            bTagScaleFactor      = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor");
            puScaleFactor        = tr.getVar<double>("puWeightCorr");

            weight *= eventweight*bTagScaleFactor*prefiringScaleFactor*puScaleFactor;
        }

        // ----------------------------------------------------
        // -- Trigger Efficiency on the JetHT Dataset and MC
        // ----------------------------------------------------
        
        if( (filetag.find("Data_JetHT") != std::string::npos || runtype == "MC") ) 
        {
            // -------------------------
            // for our analysis triggers
            // -------------------------
            const std::map<std::string, bool> cut_map_hadTriggers 
            {
                { "had_trig_pt45",   passBaseline0l_hadTrig && passTriggerAllHad },                                               
                { "had_noTrig_pt45", passBaseline0l_hadTrig }, 
            };
            fillHistos(cut_map_hadTriggers, passTriggerAllHad, HT_trigger_pt45, NGoodJets_pt45, NGoodBJets_pt45, weight);

            // -------------------------------
            // for reference analysis triggers
            // -------------------------------
            const std::map<std::string, bool> cut_map_refAnTriggers
            {  
                { "refAnHad_trig_ge2bjetCut", passBaseline0l_refAN && passTriggerRefAN                    }, 
                { "refAnHad_trig_2bjetCut",   passBaseline0l_refAN && passTriggerRefAN && pass_2bjetCut   },
                { "refAnHad_trig_3bjetCut",   passBaseline0l_refAN && passTriggerRefAN && pass_3bjetCut   },
                { "refAnHad_trig_ge4bjetCut", passBaseline0l_refAN && passTriggerRefAN && pass_ge4bjetCut },

                { "refAnHad_noTrig_ge2bjetCut", passBaseline0l_refAN                    },
                { "refAnHad_noTrig_2bjetCut",   passBaseline0l_refAN && pass_2bjetCut   },
                { "refAnHad_noTrig_3bjetCut",   passBaseline0l_refAN && pass_3bjetCut   },
                { "refAnHad_noTrig_ge4bjetCut", passBaseline0l_refAN && pass_ge4bjetCut },
            };
            fillHistosRefAN(cut_map_refAnTriggers, passTriggerRefAN, HT_trigger_pt30, SixthJetPt, NGoodBJets_pt30, weight);

        }
        
        // ---------------------------------------------------------
        // -- Trigger Efficiency on the SingleMuon Dataset and MC
        // ---------------------------------------------------------
        if ( (filetag.find("Data_SingleMuon") != std::string::npos || runtype == "MC") )
        {
            // -------------------------
            // for our analysis triggers
            // -------------------------
            const std::map<std::string, bool> cut_map_hadMuTriggers
            {   
                { "had_IsoMu_trig_pt45",   passBaseline0l_hadMuTrig && passTriggerAllHad},
                { "had_IsoMu_noTrig_pt45", passBaseline0l_hadMuTrig },  
            };
            fillHistos(cut_map_hadMuTriggers, passTriggerAllHad, HT_trigger_pt45, NGoodJets_pt45, NGoodBJets_pt45, weight);

            // -------------------------------
            // for reference analysis triggers
            // -------------------------------
            const std::map<std::string, bool> cut_map_refAnMuTriggers
            {
                { "refAnHadIsoMu_trig_ge2bjetCut", passBaseline0l_refAN && passTriggerMuonsRefAN && passTriggerRefAN                    },
                { "refAnHadIsoMu_trig_2bjetCut",   passBaseline0l_refAN && passTriggerMuonsRefAN && passTriggerRefAN && pass_2bjetCut   },
                { "refAnHadIsoMu_trig_3bjetCut",   passBaseline0l_refAN && passTriggerMuonsRefAN && passTriggerRefAN && pass_3bjetCut   },
                { "refAnHadIsoMu_trig_ge4bjetCut", passBaseline0l_refAN && passTriggerMuonsRefAN && passTriggerRefAN && pass_ge4bjetCut },

                { "refAnHadIsoMu_noTrig_ge2bjetCut", passBaseline0l_refAN && passTriggerMuonsRefAN                    },
                { "refAnHadIsoMu_noTrig_2bjetCut",   passBaseline0l_refAN && passTriggerMuonsRefAN && pass_2bjetCut   },
                { "refAnHadIsoMu_noTrig_3bjetCut",   passBaseline0l_refAN && passTriggerMuonsRefAN && pass_3bjetCut   },
                { "refAnHadIsoMu_noTrig_ge4bjetCut", passBaseline0l_refAN && passTriggerMuonsRefAN && pass_ge4bjetCut },
            };
            fillHistosRefAN(cut_map_refAnMuTriggers, passTriggerRefAN, HT_trigger_pt30, SixthJetPt, NGoodBJets_pt30, weight);
            
            // --------------------------------------------------
            // their triggers + their preselections with our pt45
            // --------------------------------------------------
            const std::map<std::string, bool> cut_map_theirTriggers
            {   
                { "TheirHadIsoMu_trig_ge2bjetCut",   passBaseline0l_refAN_pt45 && passTriggerMuonsRefAN && passTriggerRefAN                    },
                { "TheirHadIsoMu_trig_2bjetCut",     passBaseline0l_refAN_pt45 && passTriggerMuonsRefAN && passTriggerRefAN && pass_2bjetCut   },
                { "TheirHadIsoMu_trig_3bjetCut",     passBaseline0l_refAN_pt45 && passTriggerMuonsRefAN && passTriggerRefAN && pass_3bjetCut   },
                { "TheirHadIsoMu_trig_ge4bjetCut",   passBaseline0l_refAN_pt45 && passTriggerMuonsRefAN && passTriggerRefAN && pass_ge4bjetCut },
                
                { "TheirHadIsoMu_noTrig_ge2bjetCut", passBaseline0l_refAN_pt45 && passTriggerMuonsRefAN                    },
                { "TheirHadIsoMu_noTrig_2bjetCut",   passBaseline0l_refAN_pt45 && passTriggerMuonsRefAN && pass_2bjetCut   },
                { "TheirHadIsoMu_noTrig_3bjetCut",   passBaseline0l_refAN_pt45 && passTriggerMuonsRefAN && pass_3bjetCut   },
                { "TheirHadIsoMu_noTrig_ge4bjetCut", passBaseline0l_refAN_pt45 && passTriggerMuonsRefAN && pass_ge4bjetCut },
            };
            fillHistosRefAN(cut_map_theirTriggers, passTriggerRefAN, HT_trigger_pt45, SixthJetPt45, NGoodBJets_pt45, weight);
            
            // -----------------------------------------------
            // our triggers + reference analysis preselections
            // -----------------------------------------------
            const std::map<std::string, bool> cut_map_ourTriggers
            {
                { "OurHadIsoMu_trig_ge2bjetCut",   passBaseline0l_refAN && passTriggerMuonsRefAN && passTriggerAllHad                    },
                { "OurHadIsoMu_trig_2bjetCut",     passBaseline0l_refAN && passTriggerMuonsRefAN && passTriggerAllHad && pass_2bjetCut   },
                { "OurHadIsoMu_trig_3bjetCut",     passBaseline0l_refAN && passTriggerMuonsRefAN && passTriggerAllHad && pass_3bjetCut   },
                { "OurHadIsoMu_trig_ge4bjetCut",   passBaseline0l_refAN && passTriggerMuonsRefAN && passTriggerAllHad && pass_ge4bjetCut },

                { "OurHadIsoMu_noTrig_ge2bjetCut", passBaseline0l_refAN && passTriggerMuonsRefAN                    },
                { "OurHadIsoMu_noTrig_2bjetCut",   passBaseline0l_refAN && passTriggerMuonsRefAN && pass_2bjetCut   },
                { "OurHadIsoMu_noTrig_3bjetCut",   passBaseline0l_refAN && passTriggerMuonsRefAN && pass_3bjetCut   },
                { "OurHadIsoMu_noTrig_ge4bjetCut", passBaseline0l_refAN && passTriggerMuonsRefAN && pass_ge4bjetCut },
            };
            fillHistosRefAN(cut_map_ourTriggers, passTriggerAllHad, HT_trigger_pt30, SixthJetPt, NGoodBJets_pt30, weight);

        }
    }
}

void HadTriggers_Analyzer::WriteHistos(TFile* outfile)
{
    outfile->cd();

    for (const auto &p : my_histos) 
    {
        p.second->Write();
    }
    
    for (const auto &p : my_2d_histos) 
    {
        p.second->Write();
    }
    
    for (const auto &p : my_efficiencies) 
    {
        p.second->Write();
    }
    
}

bool HadTriggers_Analyzer::containsGoodHadron( const std::vector<TLorentzVector>& hadrons, const std::vector<bool>& goodHadrons, double ptThreshold, double etaSelection) 
{ 
    // Require a good hadron in JetHT data
    for( unsigned int h = 0; h < hadrons.size(); h++ ) 
    {
        if( !goodHadrons.at(h) ) continue; 
    
        TLorentzVector myHadron = hadrons.at(h);
    
        if( myHadron.Pt() >= ptThreshold && std::fabs( myHadron.Eta() ) < etaSelection ) return true;
    }

    return false;
}

void HadTriggers_Analyzer::fillHistos( const std::map<std::string, bool>& cutMap, bool passTriggerAllHad, double HT, int njet, int nbjet, double weight ) 
{
    for( auto& kv : cutMap ) 
    {
        if( kv.second ) 
        {
            my_histos["h_denominator_"+kv.first+"_HT"]->Fill( HT, weight );
            my_histos["h_denominator_"+kv.first+"_ht5000"]->Fill( HT, weight );
            my_histos["h_denominator_"+kv.first+"_NJet"]->Fill( njet, weight );
            my_histos["h_denominator_"+kv.first+"_NBJet"]->Fill( nbjet, weight );
            my_2d_histos["h_denominator_"+kv.first+"_NJetVsHT"]->Fill( njet, HT, weight );
            my_2d_histos["h_denominator_"+kv.first+"_NJetVsHt"]->Fill( njet, HT, weight );
            my_2d_histos["h_denominator_"+kv.first+"_NJetVsNBJet"]->Fill( njet, nbjet, weight );

            if( passTriggerAllHad ) 
            {
                my_histos["h_numerator_"+kv.first+"_HT"]->Fill( HT, weight );
                my_histos["h_numerator_"+kv.first+"_ht5000"]->Fill( HT, weight );
                my_histos["h_numerator_"+kv.first+"_NJet"]->Fill( njet, weight );
                my_histos["h_numerator_"+kv.first+"_NBJet"]->Fill( nbjet, weight );
                my_2d_histos["h_numerator_"+kv.first+"_NJetVsHT"]->Fill( njet, HT, weight );
                my_2d_histos["h_numerator_"+kv.first+"_NJetVsHt"]->Fill( njet, HT, weight );
                my_2d_histos["h_numerator_"+kv.first+"_NJetVsNBJet"]->Fill( njet, nbjet, weight );
            }
        }
    }
}

// function for the reference analysis histos
void HadTriggers_Analyzer::fillHistosRefAN( const std::map<std::string, bool>& cutMap, bool passTriggerRefAN, double HT, double pt, int nbjet, double weight )
{ 
    for( auto& kv : cutMap )
    {
        if( kv.second )
        {
            my_histos["h_denominator_"+kv.first+"_HT"]->Fill( HT, weight );
            my_histos["h_denominator_"+kv.first+"_6thJetPt"]->Fill( pt, weight );
            my_histos["h_denominator_"+kv.first+"_NBJet"]->Fill( nbjet, weight );
            my_2d_histos["h_denominator_"+kv.first+"_HTvs6thJetPt"]->Fill( HT, pt, weight );
            
            if( passTriggerRefAN )
            {
                my_histos["h_numerator_"+kv.first+"_HT"]->Fill( HT, weight );
                my_histos["h_numerator_"+kv.first+"_6thJetPt"]->Fill( pt, weight );
                my_histos["h_numerator_"+kv.first+"_NBJet"]->Fill( nbjet, weight );
                my_2d_histos["h_numerator_"+kv.first+"_HTvs6thJetPt"]->Fill( HT, pt, weight );
            }

        }
    }
}


void HadTriggers_Analyzer::printTriggerList( const std::vector<std::string>& TriggerNames )
{
    for( unsigned int i = 0; i < TriggerNames.size(); i++ ) 
    {
        std::string myString = TriggerNames.at(i);
        printf("%s\n", myString.c_str());
    }
}
