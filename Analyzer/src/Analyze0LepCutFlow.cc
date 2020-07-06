#define Analyze0LepCutFlow_cxx
#include "Analyzer/Analyzer/include/Analyze0LepCutFlow.h"
#include "SusyAnaTools/Tools/NTupleReader.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TRandom3.h>
#include <iostream>

Analyze0LepCutFlow::Analyze0LepCutFlow() : initHistos(false)
{
}

void Analyze0LepCutFlow::InitHistos( std::map<std::string, bool> cutMap )
{
    TH1::SetDefaultSumw2();

    // Declare all your histograms here, that way we can fill them for multiple chains

	for(auto& mycut : cutMap) {

		my_histos.emplace("h_njets_0l_" + mycut.first, std::make_shared<TH1D>( ("h_njets_0l_" + mycut.first).c_str(), ("h_njets_0l_" + mycut.first).c_str(), 20, 0, 20) );
		my_histos.emplace("h_ntops_0l_" + mycut.first, std::make_shared<TH1D>( ("h_ntops_0l_" + mycut.first).c_str(), ("h_ntops_0l_" + mycut.first).c_str(), 10, 0, 10) );
		my_histos.emplace("h_nb_0l_" + mycut.first, std::make_shared<TH1D>( ("h_nb_0l_" + mycut.first).c_str(), ("h_nb_0l_" + mycut.first).c_str(), 10, 0, 10) );
		my_histos.emplace("h_HT_0l_" + mycut.first, std::make_shared<TH1D>( ("h_HT_0l_" + mycut.first).c_str(), ("h_HT_0l_" + mycut.first).c_str(), 60, 0, 3000) );


	my_histos.emplace("h_jPt_0l_" + mycut.first, std::make_shared<TH1D>( ("h_jPt_0l_" + mycut.first).c_str(),("h_jPt_0l_" + mycut.first).c_str(), 150, 0.0, 1500.0));
	my_histos.emplace("h_jEta_0l_" + mycut.first, std::make_shared<TH1D>( ("h_jEta_0l_" + mycut.first).c_str(),("h_jEta_0l_" + mycut.first).c_str(), 200, -6.0, 6.0));
	my_histos.emplace("h_jPhi_0l_" + mycut.first, std::make_shared<TH1D>( ("h_jPhi_0l_" + mycut.first).c_str(),("h_jPhi_0l_" + mycut.first).c_str(), 200, -4.0, 4.0));
	my_histos.emplace("h_jMass_0l_" + mycut.first, std::make_shared<TH1D>( ("h_jMass_0l_" + mycut.first).c_str(),("h_jMass_0l_" + mycut.first).c_str(), 200, 0.0, 200.0));

	my_histos.emplace("h_dM_stops_0l_" + mycut.first, std::make_shared<TH1D>( ("h_dM_stops_0l_" + mycut.first).c_str(), ("h_dM_stops_0l_" + mycut.first).c_str(), 100, 0, 500) );
	my_histos.emplace("h_dR_stops_0l_" + mycut.first, std::make_shared<TH1D>( ("h_dR_stops_0l_" + mycut.first).c_str(), ("h_dR_stops_0l_" + mycut.first).c_str(), 40, 0, 4) );
	my_histos.emplace("h_dPhi_stops_0l_" + mycut.first, std::make_shared<TH1D>( ("h_dPhi_stops_0l_" + mycut.first).c_str(), ("h_dPhi_stops_0l_" + mycut.first).c_str(), 40, -4, 4) );
	my_histos.emplace("h_dEta_stops_0l_" + mycut.first, std::make_shared<TH1D>( ("h_dEta_stops_0l_" + mycut.first).c_str(), ("h_dEta_stops_0l_" + mycut.first).c_str(), 40, -4, 4) );

	my_histos.emplace("h_dR_bjets_0l_" + mycut.first, std::make_shared<TH1D>( ("h_dR_bjets_0l_" + mycut.first).c_str(), ("h_dR_bjets_0l_" + mycut.first).c_str(), 40, 0, 4) );
	my_histos.emplace("h_dPhi_bjets_0l_" + mycut.first, std::make_shared<TH1D>( ("h_dPhi_bjets_0l_" + mycut.first).c_str(), ("h_dPhi_bjets_0l_" + mycut.first).c_str(), 40, -4, 4) );
	my_histos.emplace("h_dEta_bjets_0l_" + mycut.first, std::make_shared<TH1D>( ("h_dEta_bjets_0l_" + mycut.first).c_str(), ("h_dEta_bjets_0l_" + mycut.first).c_str(), 40, -4, 4) );

    my_histos.emplace("h_AvgMass_stops_0l_" + mycut.first, std::make_shared<TH1D>( ("h_AvgMass_stops_0l_" + mycut.first ).c_str(), ("h_AvgMass_stops_0l_" + mycut.first ).c_str(), 500, 0, 1500 ) );
    my_histos.emplace("h_WeightedMassDiff_stops_0l_" + mycut.first, std::make_shared<TH1D>( ("h_WeightedMassDiff_stops_0l_" + mycut.first ).c_str(), ("h_WeightedMassDiff_stops_0l_" + mycut.first ).c_str(), 500, 0, 5 ) );

    my_2d_histos.emplace("h_Weighted_Avg_MassStop_0l_" + mycut.first, std::make_shared<TH2D>(("h_Weighted_Avg_MassStop_0l_" + mycut.first).c_str(), ("h_Weighted_Avg_MassStop_0l_" + mycut.first).c_str(), 500, 0, 5, 500, 0, 1500 ) );
	
	}

	my_histos.emplace( "h_event_sel_cutFlow",	std::make_shared<TH1D>( "h_event_sel_cutFlow",	"h_event_sel_cutFlow",	10,	0,	10));
	my_histos.emplace("h_mstop_cuts", std::make_shared<TH1D>("h_mstop_cuts", "h_mstop_cuts", 8, 0, 8) );

    my_histos.emplace("h_ntops",    std::make_shared<TH1D>("h_ntops",  "h_ntops",   10,  0,     10  ) );
    my_histos.emplace("h_ntops_j1", std::make_shared<TH1D>("h_ntops_j1","h_ntops_j1", 10,  0,     10  ) );
    my_histos.emplace("h_ntops_j2", std::make_shared<TH1D>("h_ntops_j2","h_ntops_j2", 10,  0,     10  ) );
    my_histos.emplace("h_ntops_j3", std::make_shared<TH1D>("h_ntops_j3","h_ntops_j3", 10,  0,     10  ) );
    my_histos.emplace("h_ntops_all", std::make_shared<TH1D>("h_tops_all","h_tops_all;ntops;All/merged/between/resolved", 4, 0, 4));

    my_histos.emplace("h_1j_vs_3j", std::make_shared<TH1D>("h_1j_vs_3j", "h_1j_vs_3j", 3, 0, 3 ) );
    my_histos.emplace("h_topsPt", std::make_shared<TH1D>( "h_topsPt", "h_topsPt",100, 0, 3000 ) );

    my_histos.emplace("h_Mass_stop1", std::make_shared<TH1D>( "h_Mass_stop1", "h_Mass_stop1", 500, 0, 1500 ) );
    my_histos.emplace("h_Mass_stop2", std::make_shared<TH1D>( "h_Mass_stop2", "h_Mass_stop2", 500, 0, 1500 ) );

//2D histos

    my_2d_histos.emplace("h_Mass_stop1stop2", std::make_shared<TH2D>("h_Mass_stop1stop2", "h_Mass_stop1stop2", 500, 0, 1500, 500, 0, 1500 ) );

	

    my_efficiencies.emplace("event_sel_total", std::make_shared<TEfficiency>("event_sel_total","Total event selection efficiency;Cut;#epsilon",8,0,8));

}

void Analyze0LepCutFlow::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    while( tr.getNextEvent() )
    {
        const auto& filetag	       = tr.getVar<std::string>("filetag");
	const auto& HT_trigger         = tr.getVar<double>("HT_trigger");
        const auto& ntops              = tr.getVar<int>("ntops");
        const auto& ntops_3jet         = tr.getVar<int>("ntops_3jet");
        const auto& ntops_2jet         = tr.getVar<int>("ntops_2jet");
        const auto& ntops_1jet         = tr.getVar<int>("ntops_1jet");
        const auto& runtype            = tr.getVar<std::string>("runtype");     
        const auto& NJets_pt45         = tr.getVar<int>("NGoodJets_pt45");
        const auto& NBJets_pt45        = tr.getVar<int>("NGoodBJets_pt45");
        const auto& NGoodLeptons       = tr.getVar<int>("NGoodLeptons");
	const auto& Jets 	       = tr.getVec<TLorentzVector>("Jets");
	const auto& GoodJets_pt45      = tr.getVec<bool>("GoodJets_pt45");	
	const auto& GoodBJets_pt45     = tr.getVec<bool>("GoodBJets_pt45");
	const auto& passBaseline0l     = tr.getVar<bool>("passBaseline0l_Good");
	const auto& dR_bjets	       = tr.getVar<double>("dR_bjets");
	const auto& topsPt	       = tr.getVec<double>("topsPt");
        const auto& stop1_PtRank       = tr.getVar<TLorentzVector>("stop1_PtRank_0l");
        const auto& stop2_PtRank       = tr.getVar<TLorentzVector>("stop2_PtRank_0l");

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
            double lumi_16 = 35900; // Lumi for 2016
	    double lumi_17 = 41500; // Lumi for 2017
            // Weight from NTuples
            if ( filetag.find("2016") != std::string::npos )
            	eventweight = lumi_16*Weight;
	    else if ( filetag.find("2017") != std::string::npos )
		eventweight = lumi_17*Weight;
        }

        // -------------------------------
        // -- Define cuts
        // -------------------------------

        // Global cuts
        bool pass_0l              	= NGoodLeptons==0;
        bool pass_njet_pt45    		= NJets_pt45>=6;
        bool pass_HT_trigger 	 	= HT_trigger > 500;
        bool pass_njet_pt45_2btag 	= NBJets_pt45 >= 2;
	bool pass_ge2tops	 	= ntops >= 2;
	bool pass_dR_bjets	 	= dR_bjets >= 1;
        bool pass_2tops 	 	= ntops == 2;
 	bool pass_2_1j		  	= ntops_1jet == 2;
	bool pass_2_3j			= ntops_3jet == 2;

	// Stop mass cuts 

	bool pass_mstop1_100	  	= stop1_PtRank.M() >= 100; 
	bool pass_mstop1_200	  	= stop1_PtRank.M() >= 200; 
	bool pass_mstop1_250	  	= stop1_PtRank.M() >= 250; 
	bool pass_mstop1_300 		= stop1_PtRank.M() >= 300; 
	bool pass_mstop1_350 		= stop1_PtRank.M() >= 350; 
	bool pass_mstop1_400 		= stop1_PtRank.M() >= 400;
	bool pass_mstop1_500		= stop1_PtRank.M() >= 500;

	bool pass_mstop2_100	  	= stop2_PtRank.M() >= 100; 
	bool pass_mstop2_200	  	= stop2_PtRank.M() >= 200; 
	bool pass_mstop2_250	  	= stop2_PtRank.M() >= 250; 
	bool pass_mstop2_300 		= stop2_PtRank.M() >= 300; 
	bool pass_mstop2_350 		= stop2_PtRank.M() >= 350; 
	bool pass_mstop2_400 		= stop2_PtRank.M() >= 400;
	bool pass_mstop2_500		= stop2_PtRank.M() >= 500;
	
	// Combined stop mass cuts

	bool pass_m12ge100		= pass_mstop1_100 && pass_mstop2_100;
	bool pass_m12ge200		= pass_mstop1_200 && pass_mstop2_200;
	bool pass_m12ge250		= pass_mstop1_250 && pass_mstop2_250;
	bool pass_m12ge300		= pass_mstop1_300 && pass_mstop2_300;
	bool pass_m12ge350		= pass_mstop1_350 && pass_mstop2_350;
	bool pass_m12ge400		= pass_mstop1_400 && pass_mstop2_400;
	bool pass_m12ge500		= pass_mstop1_500 && pass_mstop2_500;

	//-------------------
	//--- Stop Vars -----
	//-------------------

	double stop1_PtRank_Mass 	= stop1_PtRank.M();
	double stop1_PtRank_Pt		= stop1_PtRank.Pt();
	double stop1_PtRank_Phi 	= stop1_PtRank.Phi();
	double stop1_PtRank_Eta 	= stop1_PtRank.Eta();
	
	double stop2_PtRank_Mass 	= stop2_PtRank.M();
	double stop2_PtRank_Pt 		= stop2_PtRank.Pt();
	double stop2_PtRank_Phi 	= stop2_PtRank.Phi();
	double stop2_PtRank_Eta 	= stop2_PtRank.Eta();
		
	double dM_stop1_stop2 		= abs ( stop1_PtRank_Mass - stop2_PtRank_Mass );
	double dR_stop1_stop2		= stop1_PtRank.DeltaR( stop2_PtRank );
	double dPhi_stop1_stop2		= stop1_PtRank.DeltaPhi( stop2_PtRank );
	double dEta_stop1_stop2		= stop1_PtRank_Eta - stop2_PtRank_Eta;
	double Avg_stop1_stop2		= (stop1_PtRank_Mass + stop2_PtRank_Mass) / 2.0;
	double s12_WeightedDiff		= dM_stop1_stop2 / Avg_stop1_stop2;

	// -------------------
        // --- Fill Histos ---
        // -------------------                        

	const std::map<std::string, bool> cut_map_0l
	{
		{"",			passBaseline0l && pass_ge2tops && pass_dR_bjets },
		{"m12ge100",		passBaseline0l && pass_ge2tops && pass_dR_bjets && pass_m12ge100 },
		{"m12ge200",		passBaseline0l && pass_ge2tops && pass_dR_bjets && pass_m12ge200 },
		{"m12ge250",		passBaseline0l && pass_ge2tops && pass_dR_bjets && pass_m12ge250 },
		{"m12ge300",		passBaseline0l && pass_ge2tops && pass_dR_bjets && pass_m12ge300 },
		{"m12ge350",		passBaseline0l && pass_ge2tops && pass_dR_bjets && pass_m12ge350 },
		{"m12ge400",		passBaseline0l && pass_ge2tops && pass_dR_bjets && pass_m12ge400 },
		{"m12ge500",		passBaseline0l && pass_ge2tops && pass_dR_bjets && pass_m12ge500 },
	};

        // Initialize Histograms
        if(!initHistos)
        {
            InitHistos( cut_map_0l );
            initHistos = true;
        }


	my_histos["h_event_sel_cutFlow"]->Fill(0.5, eventweight);
	if( pass_0l ){
		my_histos["h_event_sel_cutFlow"]->Fill(1.5 ,eventweight);
		if( pass_HT_trigger ){
			my_histos["h_event_sel_cutFlow"]->Fill(2.5, eventweight);
			if( pass_njet_pt45 ){
				my_histos["h_event_sel_cutFlow"]->Fill(3.5, eventweight);
				if( pass_ge2tops ){
					my_histos["h_event_sel_cutFlow"]->Fill(4.5, eventweight);
					if( pass_njet_pt45_2btag ){
						my_histos["h_event_sel_cutFlow"]->Fill(5.5, eventweight);
						if( pass_dR_bjets && passBaseline0l ){
							my_histos["h_event_sel_cutFlow"]->Fill(6.5, eventweight);
						}
					}
				}
			}
		}
	}

	my_efficiencies["event_sel_total"]->Fill(true, 0);
	my_efficiencies["event_sel_total"]->Fill(pass_0l, 1);
	my_efficiencies["event_sel_total"]->Fill(pass_0l && pass_HT_trigger, 2);
	my_efficiencies["event_sel_total"]->Fill(pass_0l && pass_HT_trigger && pass_njet_pt45, 3);
	my_efficiencies["event_sel_total"]->Fill(pass_0l && pass_HT_trigger && pass_njet_pt45 && pass_njet_pt45_2btag, 4);
	my_efficiencies["event_sel_total"]->Fill(pass_0l && pass_HT_trigger && pass_njet_pt45 && pass_njet_pt45_2btag && pass_ge2tops, 5);


	// Fill ntops histos with only events which pass the 0 lep baseline
        
	if( !( passBaseline0l && pass_ge2tops && pass_dR_bjets ) ) continue;

	if( pass_ge2tops ) {

		my_histos["h_ntops"   ]->Fill(ntops, eventweight);        
		my_histos["h_ntops_j1"]->Fill(ntops_1jet, eventweight);
		my_histos["h_ntops_j2"]->Fill(ntops_2jet, eventweight);
		my_histos["h_ntops_j3"]->Fill(ntops_3jet, eventweight);        

		my_histos["h_ntops_all"]->Fill(0.5, ntops);
		my_histos["h_ntops_all"]->Fill(1.5, ntops_1jet);
		my_histos["h_ntops_all"]->Fill(2.5, ntops_2jet);
		my_histos["h_ntops_all"]->Fill(3.5, ntops_3jet);

	}

	if( pass_2tops ) {
		if( pass_2_1j )
			my_histos["h_1j_vs_3j"] -> Fill(0.5, eventweight);
		else if( pass_2_3j )
			my_histos["h_1j_vs_3j"] -> Fill(2.5, eventweight);
		else
			my_histos["h_1j_vs_3j"] -> Fill(1.5, eventweight);
	}

	for( auto& t : topsPt ) {

		my_histos["h_topsPt"] -> Fill(t, eventweight);

	}
	
	my_histos["h_mstop_cuts"] -> Fill( 0.5, eventweight );
	
	if( pass_m12ge100 ) {
		my_histos["h_mstop_cuts"] -> Fill( 1.5, eventweight );
	}
	if( pass_m12ge200 ) {
		my_histos["h_mstop_cuts"] -> Fill( 2.5, eventweight );
	}
	if( pass_m12ge250 ) {
		my_histos["h_mstop_cuts"] -> Fill( 3.5, eventweight );
	}
	if( pass_m12ge300 ) {
		my_histos["h_mstop_cuts"] -> Fill( 4.5, eventweight );
	}
	if( pass_m12ge350 ) {
		my_histos["h_mstop_cuts"] -> Fill( 5.5, eventweight );
	}
	if( pass_m12ge400 ) {
		my_histos["h_mstop_cuts"] -> Fill( 6.5, eventweight );
	}
	if( pass_m12ge500 ) {
		my_histos["h_mstop_cuts"] -> Fill( 7.5, eventweight );
	}

	my_histos["h_Mass_stop1"] -> Fill( stop1_PtRank.M(), eventweight );
	my_histos["h_Mass_stop2"] -> Fill( stop2_PtRank.M(), eventweight );

	my_2d_histos["h_Mass_stop1stop2"] -> Fill( stop1_PtRank.M(), stop2_PtRank.M(), eventweight );


	std::vector<TLorentzVector> bjets;
	
	for( unsigned int iJet = 0; iJet < Jets.size(); iJet++ ) {
				
		if( !GoodBJets_pt45[iJet] ) continue;
		bjets.push_back( Jets[iJet] );

	}

	for( auto& kv : cut_map_0l ) {

		if( kv.second ) {

			my_histos["h_njets_0l_"	+ kv.first] -> Fill(NJets_pt45, eventweight);
			my_histos["h_nb_0l_"	+ kv.first] -> Fill(NBJets_pt45, eventweight);
			my_histos["h_ntops_0l_"	+ kv.first] -> Fill(ntops, eventweight);
			my_histos["h_HT_0l_"	+ kv.first] -> Fill(HT_trigger, eventweight);
			my_2d_histos["h_Weighted_Avg_MassStop_0l_" + kv.first] -> Fill(s12_WeightedDiff, Avg_stop1_stop2, eventweight);
	
			for(unsigned int j = 0; j < Jets.size(); j++ ) {
			if( !(GoodJets_pt45[j]) ) continue;
				my_histos["h_jPt_0l_" 	+ kv.first] -> Fill(Jets[j].Pt(), eventweight);
				my_histos["h_jEta_0l_" 	+ kv.first] -> Fill(Jets[j].Eta(), eventweight);
				my_histos["h_jPhi_0l_" 	+ kv.first] -> Fill(Jets[j].Phi(), eventweight);
				my_histos["h_jMass_0l_"	+ kv.first] -> Fill(Jets[j].M(), eventweight);
			}

			
			my_histos["h_dM_stops_0l_" + kv.first] -> Fill(dM_stop1_stop2, eventweight);
			my_histos["h_dR_stops_0l_" + kv.first] -> Fill(dR_stop1_stop2, eventweight);
			my_histos["h_dPhi_stops_0l_" + kv.first] -> Fill(dPhi_stop1_stop2, eventweight);
			my_histos["h_dEta_stops_0l_" + kv.first] -> Fill(dEta_stop1_stop2, eventweight);
			my_histos["h_AvgMass_stops_0l_" + kv.first] -> Fill(Avg_stop1_stop2, eventweight);
			my_histos["h_WeightedMassDiff_stops_0l_" + kv.first] -> Fill(s12_WeightedDiff, eventweight);
			
			for( unsigned int i = 0; i < bjets.size(); i++) {
				for(unsigned int j = i + 1; j < bjets.size(); j++ ) {
						
						my_histos["h_dR_bjets_0l_" + kv.first] -> Fill(bjets[i].DeltaR(bjets[j]), eventweight);
						my_histos["h_dPhi_bjets_0l_" + kv.first] -> Fill(bjets[i].DeltaPhi(bjets[j]), eventweight);
						my_histos["h_dEta_bjets_0l_" + kv.first] -> Fill(bjets[i].Eta() - bjets[j].Eta(), eventweight);
						
					
		
				}	

			}
		
		}

	}

    
    } // end of event loop

}

void Analyze0LepCutFlow::WriteHistos(TFile* outfile)
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
