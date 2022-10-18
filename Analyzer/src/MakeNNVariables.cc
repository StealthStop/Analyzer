#define MakeNNVariables_cxx
#include "Analyzer/Analyzer/include/MakeNNVariables.h"
#include "NTupleReader/include/NTupleReader.h"

#include "Framework/Framework/include/MiniTupleMaker.h"
#include "Framework/Framework/include/Utility.h" 
#include "Framework/Framework/include/Muon.h"
#include "Framework/Framework/include/Electron.h"
#include "Framework/Framework/include/Photon.h"
#include "Framework/Framework/include/Jet.h"
#include "Framework/Framework/include/BJet.h"
#include "Framework/Framework/include/CommonVariables.h"
#include "Framework/Framework/include/RunTopTagger.h"
#include "Framework/Framework/include/FatJetCombine.h"
#include "Framework/Framework/include/MakeMVAVariables.h"
#include "Framework/Framework/include/Baseline.h"
#include "Framework/Framework/include/BTagCorrector.h"
#include "Framework/Framework/include/ScaleFactors.h"
#include "Framework/Framework/include/StopJets.h"
#include "Framework/Framework/include/StopGenMatch.h"
#include "Framework/Framework/include/MakeStopHemispheres.h"

#include <iostream>
#include <stdio.h> 

MakeNNVariables::MakeNNVariables()
{
    //my_channels   = {"0l", "1l", "2l"};
    my_channels   = {"2l"};
    my_splits     = {"count", "Train", "Test", "Val"};
    my_var_suffix = {"", "JECup", "JECdown", "JERup", "JERdown"};
}

void MakeNNVariables::Loop(NTupleReader& tr, double, int maxevents, bool)
{

    const auto& TopTaggerCfg    = tr.getVar<std::string>("TopTaggerCfg");
    for(const auto& myVarSuffix : my_var_suffix)
    {
        if (myVarSuffix == "") continue;
        Jet                 jet(myVarSuffix);
        BJet                bjet(myVarSuffix);
        Muon                muon(myVarSuffix);
        Photon              photon(myVarSuffix);
        Baseline            baseline(myVarSuffix);
        Electron            electron(myVarSuffix);
        StopJets            stopJets(myVarSuffix);
        RunTopTagger        topTagger(TopTaggerCfg, myVarSuffix);
        StopGenMatch        stopGenMatch(myVarSuffix);
        FatJetCombine       fatJetCombine(myVarSuffix);
        CommonVariables     commonVariables(myVarSuffix);
        MakeMVAVariables    makeMVAVariables(false, myVarSuffix, "GoodJets_pt30", false, true, 7, 2, "");
        MakeStopHemispheres stopHemispheres_OldSeed("Jets",     "GoodJets_pt20", "NGoodJets_pt20", "_OldSeed", myVarSuffix, Hemisphere::InvMassSeed);
        MakeStopHemispheres stopHemispheres_TopSeed("StopJets", "GoodStopJets",  "NGoodStopJets",  "_TopSeed", myVarSuffix, Hemisphere::TopSeed);

        // Remember, order matters here !
        // Follow what is done in Config.h
        tr.registerFunction(muon);
        tr.registerFunction(electron);
        tr.registerFunction(photon);
        tr.registerFunction(jet);
        tr.registerFunction(bjet);
        tr.registerFunction(topTagger);
        tr.registerFunction(commonVariables);
        tr.registerFunction(baseline);
        tr.registerFunction(fatJetCombine);
        tr.registerFunction(makeMVAVariables);
        tr.registerFunction(stopJets);
        tr.registerFunction(stopHemispheres_OldSeed);
        tr.registerFunction(stopHemispheres_TopSeed);
        tr.registerFunction(stopGenMatch);
    }

    while( tr.getNextEvent() )
    {
        if ( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if ( tr.getEvtNum() % 1000 == 0 ) printf( " Event %i\n", tr.getEvtNum() );        

        for(const auto& myVarSuffix : my_var_suffix)
        {
            const auto& isSignal     = tr.getVar<bool>("isSignal");
            const auto& filetag      = tr.getVar<std::string>("filetag");

            std::map<std::string, bool> baselines;
            //baselines["0l"] = tr.getVar<bool>("passBaseline0l_good"+myVarSuffix); 
            //baselines["1l"] = tr.getVar<bool>("passBaseline1l_Good"+myVarSuffix);
            baselines["2l"] = tr.getVar<bool>("passBaseline2l_Good"+myVarSuffix);

            // Add a branch containing the mass for the stop
            // In the case of signal, use the top mass
            auto& mass = tr.createDerivedVar<double>("mass", 0.0);
            if(!isSignal)
            {
                mass = 173.0;
            }
            else
            {
                for(unsigned int m = 300; m < 1500; m+=50)
                {
                    mass = (filetag.find(std::to_string(m)) != std::string::npos) ? m : mass;
                }
            }

            // Set a model field for uniquely identifying signal models and background
            // For background, use the field to denote which variation e.g. erdOn, JECup etc.
            // Signal will live in the 100 block
            // Background will live in the 0+ block
            // For JECup, add 10
            // For JECdown, add 20
            // For JERup, add 30
            // For JERdown, add 40
            auto& model = tr.createDerivedVar<int>("model", 0);
            if (isSignal) {
                if(filetag.find("RPV") != std::string::npos)
                {
                    model = 100;
                } else if (filetag.find("SYY") != std::string::npos) 
                {
                    model = 101;
                } else if (filetag.find("SHH") != std::string::npos)
                {
                    model = 102;
                }
            } else {
                if (filetag.find("TTJets") != std::string::npos)
                {
                    model = 1;
                } else if (filetag.find("erdOn") != std::string::npos)
                {   
                    model = 2;
                } else if (filetag.find("hdampUp") != std::string::npos)
                {
                    model = 3;
                } else if (filetag.find("hdampDown") != std::string::npos)
                {
                    model = 4;
                } else if (filetag.find("underlyingEvtUp") != std::string::npos)
                {
                    model = 5;
                } else if (filetag.find("underlyingEvtDown") != std::string::npos)
                {
                    model = 6;
                } else if (filetag.find("fsrUp") != std::string::npos)
                {
                    model = 7;
                } else if (filetag.find("fsrDown") != std::string::npos)
                {
                    model = 8;
                } else if (filetag.find("isrUp") != std::string::npos)
                {
                    model = 9;
                } else if (filetag.find("isrDown") != std::string::npos)
                {
                    model = 10;
                } 
            }

            // Put JEC/JER variation in a different number block by adding
            if (myVarSuffix == "JECup")
            {
                model += 10;
            } else if (myVarSuffix == "JECdown")
            {
                model += 20;
            } else if (myVarSuffix == "JERup")
            {
                model += 30;
            } else if (myVarSuffix == "JERdown")
            {
                model += 40;
            }

            if( tr.isFirstEvent() ) 
            {
                for (const auto& channel : my_channels)
                {
                    std::string myTreeName = "myMiniTree_"+channel+myVarSuffix;

                    for (const auto& split : my_splits)
                    {
                        my_counts[split][channel][myVarSuffix] = 0;

                        if (split != "count")
                        {
                            myTree[split][channel][myVarSuffix]      = new TTree( (myTreeName).c_str() , (myTreeName).c_str() );
                            myTree[split][channel][myVarSuffix]->SetDirectory(0);
                            myMiniTuple[split][channel][myVarSuffix] = new MiniTupleMaker( myTree[split][channel][myVarSuffix] );
                        }
                    }
                }

                //-----------------------------------
                //  Initialize the tree
                //-----------------------------------       
                std::set<std::string> varGeneral = 
                {
                    "FinalLumi",   
                    "mass",   
                    "model",
                    "Weight",  
                    "stop1_ptrank_mass"+myVarSuffix,
                    "stop2_ptrank_mass"+myVarSuffix,
                };

                std::set<std::string> varLeptonic =
                {
                    "Mbl"+myVarSuffix,
                    "lvMET_cm_m"+myVarSuffix,
                    "lvMET_cm_eta"+myVarSuffix,
                    "lvMET_cm_phi"+myVarSuffix,
                    "lvMET_cm_pt"+myVarSuffix,
                    "GoodLeptons_m_1"+myVarSuffix,       "GoodLeptons_m_2"+myVarSuffix,
                    "GoodLeptons_eta_1"+myVarSuffix,     "GoodLeptons_eta_2"+myVarSuffix,
                    "GoodLeptons_phi_1"+myVarSuffix,     "GoodLeptons_phi_2"+myVarSuffix,
                    "GoodLeptons_pt_1"+myVarSuffix,      "GoodLeptons_pt_2"+myVarSuffix,
                };

                std::set<std::string> varOldSeed =
                {
                    "dR_Stop1Stop2_cm_OldSeed"+myVarSuffix,
                    "dPhi_Stop1Stop2_cm_OldSeed"+myVarSuffix,
                    "Stop1_mass_cm_OldSeed"+myVarSuffix,     "Stop2_mass_cm_OldSeed"+myVarSuffix,
                    "Stop1_pt_cm_OldSeed"+myVarSuffix,       "Stop2_pt_cm_OldSeed"+myVarSuffix,
                    "Stop1_phi_cm_OldSeed"+myVarSuffix,      "Stop2_phi_cm_OldSeed"+myVarSuffix,
                    "Stop1_eta_cm_OldSeed"+myVarSuffix,      "Stop2_eta_cm_OldSeed"+myVarSuffix,
                    "Stop1_scalarPt_cm_OldSeed"+myVarSuffix, "Stop2_scalarPt_cm_OldSeed"+myVarSuffix,
                };

                std::set<std::string> varTopSeed =
                {
                    "dR_Stop1Stop2_cm_TopSeed"+myVarSuffix,
                    "dPhi_Stop1Stop2_cm_TopSeed"+myVarSuffix,
                    "Stop1_mass_cm_TopSeed"+myVarSuffix,     "Stop2_mass_cm_TopSeed"+myVarSuffix,
                    "Stop1_pt_cm_TopSeed"+myVarSuffix,       "Stop2_pt_cm_TopSeed"+myVarSuffix,
                    "Stop1_phi_cm_TopSeed"+myVarSuffix,      "Stop2_phi_cm_TopSeed"+myVarSuffix,
                    "Stop1_eta_cm_TopSeed"+myVarSuffix,      "Stop2_eta_cm_TopSeed"+myVarSuffix,
                    "Stop1_scalarPt_cm_TopSeed"+myVarSuffix, "Stop2_scalarPt_cm_TopSeed"+myVarSuffix,
                };

                std::set<std::string> varTops = 
                {
                    "top1_pt_cm"+myVarSuffix,   "top2_pt_cm"+myVarSuffix,
                    "top1_eta_cm"+myVarSuffix,  "top2_eta_cm"+myVarSuffix,
                    "top1_phi_cm"+myVarSuffix,  "top2_phi_cm"+myVarSuffix,
                    "top1_mass_cm"+myVarSuffix, "top2_mass_cm"+myVarSuffix,
                };

                std::set<std::string> var7toLastJet =
                {
                    "combined7thToLastJet_pt_cm"+myVarSuffix,
                    "combined7thToLastJet_eta_cm"+myVarSuffix,
                    "combined7thToLastJet_phi_cm"+myVarSuffix,
                    "combined7thToLastJet_m_cm"+myVarSuffix,
                    "combined7thToLastJet_E_cm"+myVarSuffix,
                };

                // -----------------------------------------------
                // get the jet variables separately for 0l, 1l, 2l
                // -----------------------------------------------

                std::string ptCut = "pt30";
                
                std::set<std::string> varChannelSpecific =
                {
                    "HT_trigger_"+ptCut+myVarSuffix,
                    "NGoodJets_"+ptCut+"_double"+myVarSuffix,
                };

                std::set<std::string> varEventShape = 
                {
                    "fwm2_top6"+myVarSuffix,    "fwm3_top6"+myVarSuffix,    "fwm4_top6"+myVarSuffix,   "fwm5_top6"+myVarSuffix,
                    "jmt_ev0_top6"+myVarSuffix, "jmt_ev1_top6"+myVarSuffix, "jmt_ev2_top6"+myVarSuffix,

                };

                std::set<std::string> varJets = 
                {
                    "Jet_m_1"+myVarSuffix,         "Jet_m_2"+myVarSuffix,         "Jet_m_3"+myVarSuffix,         "Jet_m_4"+myVarSuffix,         "Jet_m_5"+myVarSuffix,         "Jet_m_6"+myVarSuffix,         "Jet_m_7"+myVarSuffix,
                    "Jet_E_1"+myVarSuffix,         "Jet_E_2"+myVarSuffix,         "Jet_E_3"+myVarSuffix,         "Jet_E_4"+myVarSuffix,         "Jet_E_5"+myVarSuffix,         "Jet_E_6"+myVarSuffix,         "Jet_E_7"+myVarSuffix,
                    "Jet_eta_1"+myVarSuffix,       "Jet_eta_2"+myVarSuffix,       "Jet_eta_3"+myVarSuffix,       "Jet_eta_4"+myVarSuffix,       "Jet_eta_5"+myVarSuffix,       "Jet_eta_6"+myVarSuffix,       "Jet_eta_7"+myVarSuffix,
                    "Jet_phi_1"+myVarSuffix,       "Jet_phi_2"+myVarSuffix,       "Jet_phi_3"+myVarSuffix,       "Jet_phi_4"+myVarSuffix,       "Jet_phi_5"+myVarSuffix,       "Jet_phi_6"+myVarSuffix,       "Jet_phi_7"+myVarSuffix,
                    "Jet_pt_1"+myVarSuffix,        "Jet_pt_2"+myVarSuffix,        "Jet_pt_3"+myVarSuffix,        "Jet_pt_4"+myVarSuffix,        "Jet_pt_5"+myVarSuffix,        "Jet_pt_6"+myVarSuffix,        "Jet_pt_7"+myVarSuffix,
                    "Jet_flavb_1"+myVarSuffix,     "Jet_flavb_2"+myVarSuffix,     "Jet_flavb_3"+myVarSuffix,     "Jet_flavb_4"+myVarSuffix,     "Jet_flavb_5"+myVarSuffix,     "Jet_flavb_6"+myVarSuffix,     "Jet_flavb_7"+myVarSuffix,
                    "Jet_flavg_1"+myVarSuffix,     "Jet_flavg_2"+myVarSuffix,     "Jet_flavg_3"+myVarSuffix,     "Jet_flavg_4"+myVarSuffix,     "Jet_flavg_5"+myVarSuffix,     "Jet_flavg_6"+myVarSuffix,     "Jet_flavg_7"+myVarSuffix,
                    "Jet_flavc_1"+myVarSuffix,     "Jet_flavc_2"+myVarSuffix,     "Jet_flavc_3"+myVarSuffix,     "Jet_flavc_4"+myVarSuffix,     "Jet_flavc_5"+myVarSuffix,     "Jet_flavc_6"+myVarSuffix,     "Jet_flavc_7"+myVarSuffix,
                    "Jet_flavuds_1"+myVarSuffix,   "Jet_flavuds_2"+myVarSuffix,   "Jet_flavuds_3"+myVarSuffix,   "Jet_flavuds_4"+myVarSuffix,   "Jet_flavuds_5"+myVarSuffix,   "Jet_flavuds_6"+myVarSuffix,   "Jet_flavuds_7"+myVarSuffix,
                    "Jet_flavq_1"+myVarSuffix,     "Jet_flavq_2"+myVarSuffix,     "Jet_flavq_3"+myVarSuffix,     "Jet_flavq_4"+myVarSuffix,     "Jet_flavq_5"+myVarSuffix,     "Jet_flavq_6"+myVarSuffix,     "Jet_flavq_7"+myVarSuffix,
                    "Jet_CSVb_1"+myVarSuffix,      "Jet_CSVb_2"+myVarSuffix,      "Jet_CSVb_3"+myVarSuffix,      "Jet_CSVb_4"+myVarSuffix,      "Jet_CSVb_5"+myVarSuffix,      "Jet_CSVb_6"+myVarSuffix,      "Jet_CSVb_7"+myVarSuffix,
                };            

                for (std::string channel : my_channels)
                {

                    for (const auto& split : myTree)
                    {
                        myMiniTuple[split.first][channel][myVarSuffix]->setTupleVars(varGeneral);
                        myMiniTuple[split.first][channel][myVarSuffix]->setTupleVars(varEventShape); 
                        myMiniTuple[split.first][channel][myVarSuffix]->setTupleVars(varChannelSpecific);
                        myMiniTuple[split.first][channel][myVarSuffix]->setTupleVars(varJets);
                        myMiniTuple[split.first][channel][myVarSuffix]->setTupleVars(varOldSeed);
                        myMiniTuple[split.first][channel][myVarSuffix]->setTupleVars(var7toLastJet);
                        myMiniTuple[split.first][channel][myVarSuffix]->setTupleVars(varLeptonic);
                
                        if (channel == "0l") 
                        {
                            myMiniTuple[split.first][channel][myVarSuffix]->setTupleVars(varTopSeed);
                            myMiniTuple[split.first][channel][myVarSuffix]->setTupleVars(varTops);
                        }

                        myMiniTuple[split.first][channel][myVarSuffix]->initBranches(tr);
                    }
                }
            }
            
            //-----------------------------------
            //-- Fill Histograms Below
            //-----------------------------------
            for (const auto& selection : baselines)
            {
                if (selection.second)
                {
                    int mod = my_counts["count"][selection.first][myVarSuffix] % 10;
                    if(mod < 8)
                    {
                        myMiniTuple["Train"][selection.first][myVarSuffix]->fill();
                        my_counts["Train"][selection.first][myVarSuffix]++;
                    }
                    else if(mod == 8)
                    {
                        myMiniTuple["Test"][selection.first][myVarSuffix]->fill();
                        my_counts["Test"][selection.first][myVarSuffix]++;
                    }
                    else
                    {
                        myMiniTuple["Val"][selection.first][myVarSuffix]->fill();
                        my_counts["Val"][selection.first][myVarSuffix]++;
                    }
                    my_counts["count"][selection.first][myVarSuffix]++;
                }
            }
        }
    }

    for (const auto& train : my_counts)
    {
        for (const auto& channel : train.second)
        {
            for (const auto& suffix : channel.second)
            {
                std::cout << train.first+channel.first+suffix.first+": " << suffix.second << std::endl;
            }
        }
    }

}
      
void MakeNNVariables::WriteHistos( TFile* outfile ) 
{
    const auto& outFileName = std::string(outfile->GetName());
    const auto& name = utility::split("first", outFileName, ".");

    for (const auto& split : myTree)
    {
        TFile* theOutfile = TFile::Open((name+"_"+split.first+".root").c_str(), "RECREATE");
        theOutfile->cd();

        for (auto& channel : split.second)
        { 
            for (auto& suffix : channel.second)
            {
                suffix.second->Write();

                delete suffix.second;
                delete myMiniTuple[split.first][channel.first][suffix.first];
            }
        }

        theOutfile->Close();
        delete theOutfile;
    }

    remove(outFileName.c_str());
}
