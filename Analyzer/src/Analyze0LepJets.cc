#define Analyze0LepJets_cxx
#include "Analyzer/Analyzer/include/Analyze0LepJets.h"
#include "SusyAnaTools/Tools/NTupleReader.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TRandom3.h>
#include <iostream>

Analyze0LepJets::Analyze0LepJets()
{
    InitHistos();
}

void Analyze0LepJets::InitHistos()
{
    TH1::SetDefaultSumw2();

    // Declare all your histograms here, that way we can fill them for multiple chains

	my_histos.emplace("EventCounter", std::make_shared<TH1D>("EventCounter","EventCounter", 2, -1.1, 1.1 ) );

	my_histos.emplace("h_jPt", std::make_shared<TH1D>("h_jPt","h_jPt", 150, 0.0, 1500.0));
	my_histos.emplace("h_jEta", std::make_shared<TH1D>("h_jEta","h_jEta", 200, -6.0, 6.0));
	my_histos.emplace("h_jPhi", std::make_shared<TH1D>("h_jPhi","h_jPhi", 200, -4.0, 4.0));
	my_histos.emplace("h_jMass", std::make_shared<TH1D>("h_jMass","h_jMass", 200, 0.0, 200.0));

	for( unsigned int i = 1; i <= 7; i++ ){

		my_histos.emplace("h_cm_jPt_jet" + std::to_string(i), std::make_shared<TH1D>(("h_cm_jPt_jet" + std::to_string(i)).c_str(), ("h_cm_jPt_jet" + std::to_string(i)).c_str(), 150, 0.0, 1500.0) );
		my_histos.emplace("h_cm_jEta_jet" + std::to_string(i), std::make_shared<TH1D>(("h_cm_jEta_jet" + std::to_string(i)).c_str(), ("h_cm_jEta_jet" + std::to_string(i)).c_str(), 200, -6.0, 6.0) );
		my_histos.emplace("h_cm_jPhi_jet" + std::to_string(i), std::make_shared<TH1D>(("h_cm_jPhi_jet" + std::to_string(i)).c_str(), ("h_cm_jPhi_jet" + std::to_string(i)).c_str(), 200, -4.0, 4.0) );
		my_histos.emplace("h_cm_jMass_jet" + std::to_string(i), std::make_shared<TH1D>(("h_cm_jMass_jet" + std::to_string(i)).c_str(), ("h_cm_jMass_jet" + std::to_string(i)).c_str(), 200, 0.0, 200.0) );
		
		my_histos.emplace("h_lab_jPt_jet" + std::to_string(i), std::make_shared<TH1D>(("h_lab_jPt_jet" + std::to_string(i)).c_str(), ("h_lab_jPt_jet" + std::to_string(i)).c_str(), 150, 0.0, 1500.0) );
		my_histos.emplace("h_lab_jEta_jet" + std::to_string(i), std::make_shared<TH1D>(("h_lab_jEta_jet" + std::to_string(i)).c_str(), ("h_lab_jEta_jet" + std::to_string(i)).c_str(), 200, -6.0, 6.0) );
		my_histos.emplace("h_lab_jPhi_jet" + std::to_string(i), std::make_shared<TH1D>(("h_lab_jPhi_jet" + std::to_string(i)).c_str(), ("h_lab_jPhi_jet" + std::to_string(i)).c_str(), 200, -4.0, 4.0) );
		my_histos.emplace("h_lab_jMass_jet" + std::to_string(i), std::make_shared<TH1D>(("h_lab_jMass_jet" + std::to_string(i)).c_str(), ("h_lab_jMass_jet" + std::to_string(i)).c_str(), 200, 0.0, 200.0) );

	}

}

void Analyze0LepJets::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    while( tr.getNextEvent() )
    {
        const auto& runtype            = tr.getVar<std::string>("runtype");     
	const auto& passBaseline0l     = tr.getVar<bool>("passBaseline0l_Good");
	const auto& Jets 	       = tr.getVec<TLorentzVector>("Jets");
	const auto& GoodJets_pt45      = tr.getVec<bool>("GoodJets_pt45");	
	const auto& Jets_cm_top6       = tr.getVec<TLorentzVector>("Jets_cm_top6");
	const auto& Jets_top6	       = tr.getVec<TLorentzVector>("Jets_top6");
	const auto& eventCounter       = tr.getVar<int>("eventCounter");	
		
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

        // -------------------
        // --- Fill Histos ---
        // -------------------                        

        // Initialize Histograms
        if(!initHistos)
        {
            InitHistos();
            initHistos = true;
        }


	// Fill histograms

	my_histos["EventCounter"	] -> Fill(eventCounter);
	
	if(!(passBaseline0l)) continue;
	
	for(unsigned int j = 0; j < Jets.size(); j++ ) {
	if( !(GoodJets_pt45[j]) ) continue;
		my_histos["h_jPt"	] -> Fill(Jets[j].Pt(), eventweight);
		my_histos["h_jEta"	] -> Fill(Jets[j].Eta(), eventweight);
		my_histos["h_jPhi"	] -> Fill(Jets[j].Phi(), eventweight);
		my_histos["h_jMass"	] -> Fill(Jets[j].M(), eventweight);
	}
	 
	for( unsigned int j = 0; j < Jets_cm_top6.size(); j++ ) {
		my_histos["h_cm_jPt_jet" + std::to_string(j+1)	] -> Fill(Jets_cm_top6[j].Pt(), eventweight);
		my_histos["h_cm_jEta_jet" + std::to_string(j+1)	] -> Fill(Jets_cm_top6[j].Eta(), eventweight);
		my_histos["h_cm_jPhi_jet" + std::to_string(j+1)	] -> Fill(Jets_cm_top6[j].Phi(), eventweight);
		my_histos["h_cm_jMass_jet" + std::to_string(j+1)] -> Fill(Jets_cm_top6[j].M(), eventweight);
	}       
	
	for( unsigned int j = 0; j < Jets_top6.size(); j++ ) {
		my_histos["h_lab_jPt_jet" + std::to_string(j+1)	] -> Fill(Jets_top6[j].Pt(), eventweight);
		my_histos["h_lab_jEta_jet" + std::to_string(j+1)] -> Fill(Jets_top6[j].Eta(), eventweight);
		my_histos["h_lab_jPhi_jet" + std::to_string(j+1)] -> Fill(Jets_top6[j].Phi(), eventweight);
		my_histos["h_lab_jMass_jet" + std::to_string(j+1)] -> Fill(Jets_top6[j].M(), eventweight);
	}       

    } // end of event loop

}

void Analyze0LepJets::WriteHistos(TFile* outfile)
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
