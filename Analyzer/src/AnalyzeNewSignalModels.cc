#define AnalyzeNewSignalModels_cxx
#include "Analyzer/Analyzer/include/AnalyzeNewSignalModels.h"
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

AnalyzeNewSignalModels::AnalyzeNewSignalModels()
{
    InitHistos();
}

//Define all your histograms here. 
void AnalyzeNewSignalModels::InitHistos()
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    //This event counter histogram is necessary so that we know that all the condor jobs ran successfully. If not, when you use the hadder script, you will see a discrepancy in red as the files are being hadded.
    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ;

    //1D histograms of njet distributions for two different categories of NonIsoMuons 
    //-----------------------------------
    //          Njets histos
    //-----------------------------------
    
    my_histos.emplace ( "h_RPV_2t6j_mStop-300_Njets", std::make_shared<TH1D>("h_RPV_mStop-300_Njets", "Njets for RPV stop mass 300", 15 ,0 , 15) );
    my_histos.emplace ( "h_RPV_2t6j_mStop-350_Njets", std::make_shared<TH1D>("h_RPV_mStop-350_Njets", "Njets for RPV stop mass 350", 15 ,0 , 15) );
    my_histos.emplace ( "h_RPV_2t6j_mStop-800_Njets", std::make_shared<TH1D>("h_RPV_mStop-800_Njets", "Njets for RPV stop mass 800", 15 ,0 , 15) );
    my_histos.emplace ( "h_RPV_2t6j_mStop-850_Njets", std::make_shared<TH1D>("h_RPV_mStop-850_Njets", "Njets for RPV stop mass 850", 15 ,0 , 15) ); 
    my_histos.emplace ( "h_RPV_2t6j_mStop-1300_Njets", std::make_shared<TH1D>("h_RPV_mStop-1300_Njets", "Njets for RPV stop mass 1300", 15 ,0 , 15) );
    my_histos.emplace ( "h_RPV_2t6j_mStop-1350_Njets", std::make_shared<TH1D>("h_RPV_mStop-1350_Njets", "Njets for RPV stop mass 1350", 15 ,0 , 15) );

    my_histos.emplace ( "h_StealthSYY_2t6j_mStop-300_Njets", std::make_shared<TH1D>("h_StealthSYY_2t6j_mStop-300_Njets", "Njets for Stealth SYY stop mass 300", 15 ,0 , 15) );
    my_histos.emplace ( "h_StealthSYY_2t6j_mStop-350_Njets", std::make_shared<TH1D>("h_StealthSYY_2t6j_mStop-350_Njets", "Njets for Stealth SYY stop mass 350", 15 ,0 , 15) );
    my_histos.emplace ( "h_StealthSYY_2t6j_mStop-800_Njets", std::make_shared<TH1D>("h_StealthSYY_2t6j_mStop-800_Njets", "Njets for Stealth SYY stop mass 800", 15 ,0 , 15) );
    my_histos.emplace ( "h_StealthSYY_2t6j_mStop-850_Njets", std::make_shared<TH1D>("h_StealthSYY_2t6j_mStop-850_Njets", "Njets for Stealth SYY stop mass 850", 15 ,0 , 15) ); 
    my_histos.emplace ( "h_StealthSYY_2t6j_mStop-1300_Njets", std::make_shared<TH1D>("h_StealthSYY_2t6j_mStop-1300_Njets", "Njets for Stealth SYY stop mass 1300", 15 ,0 , 15) );
    my_histos.emplace ( "h_StealthSYY_2t6j_mStop-1350_Njets", std::make_shared<TH1D>("h_StealthSYY_2t6j_mStop-1350_Njets", "Njets for Stealth SYY stop mass 1350", 15 ,0 , 15) );


    my_histos.emplace ( "h_StealthSHH_2t4b_mStop-300_Njets", std::make_shared<TH1D>("h_StealthSHH_2t4b_mStop-300_Njets", "Njets for Stealth SHH stop mass 300", 15 ,0 , 15) );
    my_histos.emplace ( "h_StealthSHH_2t4b_mStop-350_Njets", std::make_shared<TH1D>("h_StealthSHH_2t4b_mStop-350_Njets", "Njets for Stealth SHH stop mass 350", 15 ,0 , 15) );
    my_histos.emplace ( "h_StealthSHH_2t4b_mStop-800_Njets", std::make_shared<TH1D>("h_StealthSHH_2t4b_mStop-800_Njets", "Njets for Stealth SHH stop mass 800", 15 ,0 , 15) );
    my_histos.emplace ( "h_StealthSHH_2t4b_mStop-850_Njets", std::make_shared<TH1D>("h_StealthSHH_2t4b_mStop-850_Njets", "Njets for Stealth SHH stop mass 850", 15 ,0 , 15) ); 
    my_histos.emplace ( "h_StealthSHH_2t4b_mStop-1300_Njets", std::make_shared<TH1D>("h_StealthSHH_2t4b_mStop-1300_Njets", "Njets for Stealth SHH stop mass 1300", 15 ,0 , 15) );
    my_histos.emplace ( "h_StealthSHH_2t4b_mStop-1350_Njets", std::make_shared<TH1D>("h_StealthSHH_2t4b_mStop-1350_Njets", "Njets for Stealth SHH stop mass 1350", 15 ,0 , 15) );

    //-----------------------------------
    //           MVA histos
    //-----------------------------------
        
    my_histos.emplace ( "h_RPV_2t6j_mStop-300_MVA", std::make_shared<TH1D>("h_RPV_mStop-300_MVA", "MVA for RPV stop mass 300", 20 , 0 , 1) );
    my_histos.emplace ( "h_RPV_2t6j_mStop-350_MVA", std::make_shared<TH1D>("h_RPV_mStop-350_MVA", "MVA for RPV stop mass 350", 20 , 0 , 1) );
    my_histos.emplace ( "h_RPV_2t6j_mStop-800_MVA", std::make_shared<TH1D>("h_RPV_mStop-800_MVA", "MVA for RPV stop mass 800", 20 , 0 , 1) );
    my_histos.emplace ( "h_RPV_2t6j_mStop-850_MVA", std::make_shared<TH1D>("h_RPV_mStop-850_MVA", "MVA for RPV stop mass 850", 20 , 0 , 1) ); 
    my_histos.emplace ( "h_RPV_2t6j_mStop-1300_MVA", std::make_shared<TH1D>("h_RPV_mStop-1300_MVA", "MVA for RPV stop mass 1300", 20 , 0 , 1) );
    my_histos.emplace ( "h_RPV_2t6j_mStop-1350_MVA", std::make_shared<TH1D>("h_RPV_mStop-1350_MVA", "MVA for RPV stop mass 1350", 20 , 0 , 1) );

    my_histos.emplace ( "h_StealthSYY_2t6j_mStop-300_MVA", std::make_shared<TH1D>("h_StealthSYY_2t6j_mStop-300_MVA", "MVA for Stealth SYY stop mass 300", 20 , 0 , 1) );
    my_histos.emplace ( "h_StealthSYY_2t6j_mStop-350_MVA", std::make_shared<TH1D>("h_StealthSYY_2t6j_mStop-350_MVA", "MVA for Stealth SYY stop mass 350", 20 , 0 , 1) );
    my_histos.emplace ( "h_StealthSYY_2t6j_mStop-800_MVA", std::make_shared<TH1D>("h_StealthSYY_2t6j_mStop-800_MVA", "MVA for Stealth SYY stop mass 800", 20 , 0 , 1) );
    my_histos.emplace ( "h_StealthSYY_2t6j_mStop-850_MVA", std::make_shared<TH1D>("h_StealthSYY_2t6j_mStop-850_MVA", "MVA for Stealth SYY stop mass 850", 20 , 0 , 1) ); 
    my_histos.emplace ( "h_StealthSYY_2t6j_mStop-1300_MVA", std::make_shared<TH1D>("h_StealthSYY_2t6j_mStop-1300_MVA", "MVA for Stealth SYY stop mass 1300", 20 , 0 , 1) );
    my_histos.emplace ( "h_StealthSYY_2t6j_mStop-1350_MVA", std::make_shared<TH1D>("h_StealthSYY_2t6j_mStop-1350_MVA", "MVA for Stealth SYY stop mass 1350", 20 , 0 , 1) );


    my_histos.emplace ( "h_StealthSHH_2t4b_mStop-300_MVA", std::make_shared<TH1D>("h_StealthSHH_2t4b_mStop-300_MVA", "MVA for Stealth SHH stop mass 300", 20 , 0 , 1) );
    my_histos.emplace ( "h_StealthSHH_2t4b_mStop-350_MVA", std::make_shared<TH1D>("h_StealthSHH_2t4b_mStop-350_MVA", "MVA for Stealth SHH stop mass 350", 20 , 0 , 1) );
    my_histos.emplace ( "h_StealthSHH_2t4b_mStop-800_MVA", std::make_shared<TH1D>("h_StealthSHH_2t4b_mStop-800_MVA", "MVA for Stealth SHH stop mass 800", 20 , 0 , 1) );
    my_histos.emplace ( "h_StealthSHH_2t4b_mStop-850_MVA", std::make_shared<TH1D>("h_StealthSHH_2t4b_mStop-850_MVA", "MVA for Stealth SHH stop mass 850", 20 , 0 , 1) ); 
    my_histos.emplace ( "h_StealthSHH_2t4b_mStop-1300_MVA", std::make_shared<TH1D>("h_StealthSHH_2t4b_mStop-1300_MVA", "MVA for Stealth SHH stop mass 1300", 20 , 0 , 1) );
    my_histos.emplace ( "h_StealthSHH_2t4b_mStop-1350_MVA", std::make_shared<TH1D>("h_StealthSHH_2t4b_mStop-1350_MVA", "MVA for Stealth SHH stop mass 1350", 20 , 0 , 1) );
    

    //-----------------------------------
    //     2D Njets vs MVA Histos
    //-----------------------------------


    my_2d_histos.emplace ( "h_RPV_2t6j_mStop-300_Njets_MVA", std::make_shared<TH2D>("h_RPV_2t6j_mStop-300_Njets_MVA", "h_RPV_2t6j_mStop-300_Njets_MVA", 15, 0, 15, 20, 0, 1) );
    my_2d_histos.emplace ( "h_RPV_2t6j_mStop-350_Njets_MVA", std::make_shared<TH2D>("h_RPV_2t6j_mStop-350_Njets_MVA", "h_RPV_2t6j_mStop-350_Njets_MVA", 15, 0, 15, 20, 0, 1) );
    my_2d_histos.emplace ( "h_RPV_2t6j_mStop-800_Njets_MVA", std::make_shared<TH2D>("h_RPV_2t6j_mStop-800_Njets_MVA", "h_RPV_2t6j_mStop-800_Njets_MVA", 15, 0, 15, 20, 0, 1) );
    my_2d_histos.emplace ( "h_RPV_2t6j_mStop-850_Njets_MVA", std::make_shared<TH2D>("h_RPV_2t6j_mStop-850_Njets_MVA", "h_RPV_2t6j_mStop-850_Njets_MVA", 15, 0, 15, 20, 0, 1) );
    my_2d_histos.emplace ( "h_RPV_2t6j_mStop-1300_Njets_MVA", std::make_shared<TH2D>("h_RPV_2t6j_mStop-1300_Njets_MVA", "h_RPV_2t6j_mStop-1300_Njets_MVA", 15, 0, 15, 20, 0, 1) );
    my_2d_histos.emplace ( "h_RPV_2t6j_mStop-1350_Njets_MVA", std::make_shared<TH2D>("h_RPV_2t6j_mStop-1350_Njets_MVA", "h_RPV_2t6j_mStop-1350_Njets_MVA", 15, 0, 15, 20, 0, 1) );

    my_2d_histos.emplace ( "h_StealthSYY_2t6j_mStop-300_Njets_MVA", std::make_shared<TH2D>("h_StealthSYY_2t6j_mStop-300_Njets_MVA", "h_StealthSYY_2t6j_mStop-300_Njets_MVA", 15, 0, 15, 20, 0, 1) );
    my_2d_histos.emplace ( "h_StealthSYY_2t6j_mStop-350_Njets_MVA", std::make_shared<TH2D>("h_StealthSYY_2t6j_mStop-350_Njets_MVA", "h_StealthSYY_2t6j_mStop-350_Njets_MVA", 15, 0, 15, 20, 0, 1) );
    my_2d_histos.emplace ( "h_StealthSYY_2t6j_mStop-800_Njets_MVA", std::make_shared<TH2D>("h_StealthSYY_2t6j_mStop-800_Njets_MVA", "h_StealthSYY_2t6j_mStop-800_Njets_MVA", 15, 0, 15, 20, 0, 1) );
    my_2d_histos.emplace ( "h_StealthSYY_2t6j_mStop-850_Njets_MVA", std::make_shared<TH2D>("h_StealthSYY_2t6j_mStop-850_Njets_MVA", "h_StealthSYY_2t6j_mStop-850_Njets_MVA", 15, 0, 15, 20, 0, 1) );
    my_2d_histos.emplace ( "h_StealthSYY_2t6j_mStop-1300_Njets_MVA", std::make_shared<TH2D>("h_StealthSYY_2t6j_mStop-1300_Njets_MVA", "h_StealthSYY_2t6j_mStop-1300_Njets_MVA", 15, 0, 15, 20, 0, 1) );
    my_2d_histos.emplace ( "h_StealthSYY_2t6j_mStop-1350_Njets_MVA", std::make_shared<TH2D>("h_StealthSYY_2t6j_mStop-1350_Njets_MVA", "h_StealthSYY_2t6j_mStop-1350_Njets_MVA", 15, 0, 15, 20, 0, 1) );    

    my_2d_histos.emplace ( "h_StealthSHH_2t4b_mStop-300_Njets_MVA", std::make_shared<TH2D>("h_StealthSHH_2t4b_mStop-300_Njets_MVA", "h_StealthSHH_2t4b_mStop-300_Njets_MVA", 15, 0, 15, 20, 0, 1) );
    my_2d_histos.emplace ( "h_StealthSHH_2t4b_mStop-350_Njets_MVA", std::make_shared<TH2D>("h_StealthSHH_2t4b_mStop-350_Njets_MVA", "h_StealthSHH_2t4b_mStop-350_Njets_MVA", 15, 0, 15, 20, 0, 1) );
    my_2d_histos.emplace ( "h_StealthSHH_2t4b_mStop-800_Njets_MVA", std::make_shared<TH2D>("h_StealthSHH_2t4b_mStop-800_Njets_MVA", "h_StealthSHH_2t4b_mStop-800_Njets_MVA", 15, 0, 15, 20, 0, 1) );
    my_2d_histos.emplace ( "h_StealthSHH_2t4b_mStop-850_Njets_MVA", std::make_shared<TH2D>("h_StealthSHH_2t4b_mStop-850_Njets_MVA", "h_StealthSHH_2t4b_mStop-850_Njets_MVA", 15, 0, 15, 20, 0, 1) );
    my_2d_histos.emplace ( "h_StealthSHH_2t4b_mStop-1300_Njets_MVA", std::make_shared<TH2D>("h_StealthSHH_2t4b_mStop-1300_Njets_MVA", "h_StealthSHH_2t4b_mStop-1300_Njets_MVA", 15, 0, 15, 20, 0, 1) );
    my_2d_histos.emplace ( "h_StealthSHH_2t4b_mStop-1350_Njets_MVA", std::make_shared<TH2D>("h_StealthSHH_2t4b_mStop-1350_Njets_MVA", "h_StealthSHH_2t4b_mStop-1350_Njets_MVA", 15, 0, 15, 20, 0, 1) );

}

//Put everything you want to do per event here.
void AnalyzeNewSignalModels::Loop(NTupleReader& tr, double, int maxevents, bool)
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

	const auto& NGoodJets_pt30      = tr.getVar<int>("NGoodJets_pt30");
	const auto& deepESM_val         = tr.getVar<double>("deepESM_val");
	const auto& GoodLeptons		= tr.getVec<std::pair<std::string, TLorentzVector>>("GoodLeptons");	
	const auto& NGoodLeptons	= tr.getVar<int>("NGoodLeptons");

        const auto& passMadHT           = tr.getVar<bool>("passMadHT");
	const auto& passBaseline   	= tr.getVar<bool>("passBaseline1l_Good");
        
       
        if(maxevents != -1 && tr.getEvtNum() >= maxevents) break;        
        if ( tr.getEvtNum() % 1000 == 0 ) printf("  Event %i\n", tr.getEvtNum() ) ;

        // ------------------------
        // -- Define weight
        // ------------------------
        double weight               = 1.0;
        double eventweight          = 1.0;
        double bTagWeight	    = 1.0;
        double htDerivedweight	    = 1.0;
        double prefiringScaleFactor = 1.0;
        double pileupWeight         = 1.0;        
	double leptonScaleFactor    = 1.0;

        if(runtype == "MC")
        {
            if( !passMadHT ) continue; //Make sure not to double count DY events
            // Define Lumi weight
            const auto& Weight  = tr.getVar<double>("Weight");
            const auto& lumi = tr.getVar<double>("Lumi");
            eventweight = lumi*Weight;
            
            if(NGoodLeptons == 1)
            {
                const auto& eleLepWeight = tr.getVar<double>("totGoodElectronSF");
                const auto& muLepWeight  = tr.getVar<double>("totGoodMuonSF");
                leptonScaleFactor = (GoodLeptons[0].first == "e") ? eleLepWeight : muLepWeight;
            }
	
            //PileupWeight = tr.getVar<double>("_PUweightFactor");
            bTagWeight   = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            htDerivedweight = tr.getVar<double>("htDerivedweight");
            prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor");
	    pileupWeight = tr.getVar<double>("puWeightCorr");            

            weight *= eventweight*bTagWeight*prefiringScaleFactor*pileupWeight*htDerivedweight*leptonScaleFactor;
        }
        
        //Signal region  cut
	if( !passBaseline ) continue;

	//Make cuts based on the mass point mass value
	std::string suffix = filetag.substr(filetag.find("-") + 1);

	if( filetag.find("RPV") != std::string::npos ){
	  my_histos["h_RPV_2t6j_mStop-" + suffix + "_Njets"]-> Fill( NGoodJets_pt30, weight);
	  my_histos["h_RPV_2t6j_mStop-" + suffix + "_MVA"]-> Fill( deepESM_val, weight);
	  my_2d_histos["h_RPV_2t6j_mStop-" + suffix + "_Njets_MVA"]-> Fill( NGoodJets_pt30, deepESM_val, weight);
	}

	if( filetag.find("StealthSYY") != std::string::npos ){
	  my_histos["h_StealthSYY_2t6j_mStop-" + suffix + "_Njets"]-> Fill( NGoodJets_pt30, weight);
	  my_histos["h_StealthSYY_2t6j_mStop-" + suffix + "_MVA"]-> Fill( deepESM_val, weight);
	  my_2d_histos["h_StealthSYY_2t6j_mStop-" + suffix + "_Njets_MVA"]-> Fill( NGoodJets_pt30, deepESM_val, weight);
	}

	if( filetag.find("StealthSHH") != std::string::npos ){
	  my_histos["h_StealthSHH_2t4b_mStop-" + suffix + "_Njets"]-> Fill( NGoodJets_pt30, weight);
	  my_histos["h_StealthSHH_2t4b_mStop-" + suffix + "_MVA"]-> Fill( deepESM_val, weight);
	  my_2d_histos["h_StealthSHH_2t4b_mStop-" + suffix + "_Njets_MVA"]-> Fill( NGoodJets_pt30, deepESM_val, weight);
        }
	
    }//end while loop
}

void AnalyzeNewSignalModels::WriteHistos(TFile* outfile)
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
