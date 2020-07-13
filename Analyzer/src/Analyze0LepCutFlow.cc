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
	
    my_2d_histos.emplace("h_Mass_stop1stop2_0l_" + mycut.first, std::make_shared<TH2D>( ("h_Mass_stop1stop2_0l_" + mycut.first ).c_str(), ( "h_Mass_stop1stop2_0l_" + mycut.first ).c_str(), 500, 0, 1500, 500, 0, 1500 ) );
	
	my_histos.emplace("h_dM_stops_TaggedTop_" + mycut.first, std::make_shared<TH1D>( ("h_dM_stops_TaggedTop_" + mycut.first).c_str(), ("h_dM_stops_TaggedTop_" + mycut.first).c_str(), 100, 0, 500) );
	my_histos.emplace("h_dR_stops_TaggedTop_" + mycut.first, std::make_shared<TH1D>( ("h_dR_stops_TaggedTop_" + mycut.first).c_str(), ("h_dR_stops_TaggedTop_" + mycut.first).c_str(), 40, 0, 4) );
	my_histos.emplace("h_dPhi_stops_TaggedTop_" + mycut.first, std::make_shared<TH1D>( ("h_dPhi_stops_TaggedTop_" + mycut.first).c_str(), ("h_dPhi_stops_TaggedTop_" + mycut.first).c_str(), 40, -4, 4) );
	my_histos.emplace("h_dEta_stops_TaggedTop_" + mycut.first, std::make_shared<TH1D>( ("h_dEta_stops_TaggedTop_" + mycut.first).c_str(), ("h_dEta_stops_TaggedTop_" + mycut.first).c_str(), 40, -4, 4) );

    my_histos.emplace("h_AvgMass_stops_TaggedTop_" + mycut.first, std::make_shared<TH1D>( ("h_AvgMass_stops_TaggedTop_" + mycut.first ).c_str(), ("h_AvgMass_stops_TaggedTop_" + mycut.first ).c_str(), 500, 0, 1500 ) );
    my_histos.emplace("h_WeightedMassDiff_stops_TaggedTop_" + mycut.first, std::make_shared<TH1D>( ("h_WeightedMassDiff_stops_TaggedTop_" + mycut.first ).c_str(), ("h_WeightedMassDiff_stops_TaggedTop_" + mycut.first ).c_str(), 500, 0, 5 ) );

    my_2d_histos.emplace("h_Weighted_Avg_MassStop_TaggedTop_" + mycut.first, std::make_shared<TH2D>(("h_Weighted_Avg_MassStop_TaggedTop_" + mycut.first).c_str(), ("h_Weighted_Avg_MassStop_TaggedTop_" + mycut.first).c_str(), 500, 0, 5, 500, 0, 1500 ) );
	
    my_2d_histos.emplace("h_Mass_stop1stop2_TaggedTop_" + mycut.first, std::make_shared<TH2D>( ("h_Mass_stop1stop2_TaggedTop_" + mycut.first ).c_str(), ( "h_Mass_stop1stop2_TaggedTop_" + mycut.first ).c_str(), 500, 0, 1500, 500, 0, 1500 ) );
	
    }

	my_histos.emplace( "h_event_sel_cutFlow",	std::make_shared<TH1D>( "h_event_sel_cutFlow",	"h_event_sel_cutFlow",	10,	0,	10));
	my_histos.emplace("h_mstop_cuts_0l", std::make_shared<TH1D>("h_mstop_cuts_0l", "h_mstop_cuts_0l", 8, 0, 8) );
	my_histos.emplace("h_avg_diff_cuts_0l", std::make_shared<TH1D>("h_avg_diff_cuts_0l", "h_avg_diff_cuts_0l", 8, 0, 8 ) );
	my_histos.emplace("h_mstop_cuts_TaggedTop", std::make_shared<TH1D>("h_mstop_cuts_TaggedTop", "h_mstop_cuts_TaggedTop", 8, 0, 8) );
	my_histos.emplace("h_avg_diff_cuts_TaggedTop", std::make_shared<TH1D>("h_avg_diff_cuts_TaggedTop", "h_avg_diff_cuts_TaggedTop", 8, 0, 8 ) );
    my_histos.emplace("h_ntops",    std::make_shared<TH1D>("h_ntops",  "h_ntops",   10,  0,     10  ) );
    my_histos.emplace("h_ntops_j1", std::make_shared<TH1D>("h_ntops_j1","h_ntops_j1", 10,  0,     10  ) );
    my_histos.emplace("h_ntops_j2", std::make_shared<TH1D>("h_ntops_j2","h_ntops_j2", 10,  0,     10  ) );
    my_histos.emplace("h_ntops_j3", std::make_shared<TH1D>("h_ntops_j3","h_ntops_j3", 10,  0,     10  ) );
    my_histos.emplace("h_ntops_all", std::make_shared<TH1D>("h_tops_all","h_tops_all;ntops;All/merged/between/resolved", 4, 0, 4));

    my_histos.emplace("h_1j_vs_3j", std::make_shared<TH1D>("h_1j_vs_3j", "h_1j_vs_3j", 3, 0, 3 ) );
    my_histos.emplace("h_topsPt", std::make_shared<TH1D>( "h_topsPt", "h_topsPt",100, 0, 3000 ) );

    my_histos.emplace("h_Mass_stop1_0l", std::make_shared<TH1D>( "h_Mass_stop1_0l", "h_Mass_stop1_0l", 500, 0, 1500 ) );
    my_histos.emplace("h_Mass_stop2_0l", std::make_shared<TH1D>( "h_Mass_stop2_0l", "h_Mass_stop2_0l", 500, 0, 1500 ) );
    my_histos.emplace("h_Mass_stop1_TaggedTop", std::make_shared<TH1D>( "h_Mass_stop1_TaggedTop", "h_Mass_stop1_TaggedTop", 500, 0, 1500 ) );
    my_histos.emplace("h_Mass_stop2_TaggedTop", std::make_shared<TH1D>( "h_Mass_stop2_TaggedTop", "h_Mass_stop2_TaggedTop", 500, 0, 1500 ) );

    my_histos.emplace("h_nISR", std::make_shared<TH1D>( "h_nISR", "h_nISR", 20, 0, 20) );

//2D histos

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
        
        const auto& GenParticles    = tr.getVec<TLorentzVector>("GenParticles");
        const auto& GenParticlesPDG = tr.getVec<int>("GenParticles_PdgId");
        const auto& GenParticles_Status     = tr.getVec<int>("GenParticles_Status");

        const auto& stop1_PtRank_0l         = tr.getVar<TLorentzVector>("stop1_PtRank_0l");
        const auto& stop2_PtRank_0l         = tr.getVar<TLorentzVector>("stop2_PtRank_0l");

        const auto& stop1_PtRank_TaggedTop  = tr.getVar<TLorentzVector>("stop1_PtRank_TaggedTop");
        const auto& stop2_PtRank_TaggedTop  = tr.getVar<TLorentzVector>("stop2_PtRank_TaggedTop");

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

        int ISRCount = 0;

        for( unsigned int i = 0; i < GenParticles.size(); i++ ) {

            if ( ( (GenParticles_Status.at(i) >= 20 && GenParticles_Status.at(i) <=29)
                || (GenParticles_Status.at(i) >= 40 && GenParticles_Status.at(i) <=49) )
                && GenParticlesPDG.at(i) == 24 ) {

                std::cout << "Found ISR" << std::endl;
                ISRCount++;

                }

        }

        my_histos["h_nISR"] -> Fill(ISRCount, eventweight);
        
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

        bool pass_mstop1_100_0l	  	= stop1_PtRank_0l.M() >= 100; 
        bool pass_mstop1_200_0l	  	= stop1_PtRank_0l.M() >= 200; 
        bool pass_mstop1_250_0l	  	= stop1_PtRank_0l.M() >= 250; 
        bool pass_mstop1_300_0l 	= stop1_PtRank_0l.M() >= 300; 
        bool pass_mstop1_350_0l 	= stop1_PtRank_0l.M() >= 350; 
        bool pass_mstop1_400_0l 	= stop1_PtRank_0l.M() >= 400;
        bool pass_mstop1_500_0l		= stop1_PtRank_0l.M() >= 500;

        bool pass_mstop2_100_0l	  	= stop2_PtRank_0l.M() >= 100; 
        bool pass_mstop2_200_0l	  	= stop2_PtRank_0l.M() >= 200; 
        bool pass_mstop2_250_0l	  	= stop2_PtRank_0l.M() >= 250; 
        bool pass_mstop2_300_0l 	= stop2_PtRank_0l.M() >= 300; 
        bool pass_mstop2_350_0l 	= stop2_PtRank_0l.M() >= 350; 
        bool pass_mstop2_400_0l 	= stop2_PtRank_0l.M() >= 400;
        bool pass_mstop2_500_0l		= stop2_PtRank_0l.M() >= 500;
        
        bool pass_mstop1_100_TaggedTop	  	= stop1_PtRank_TaggedTop.M() >= 100; 
        bool pass_mstop1_200_TaggedTop	  	= stop1_PtRank_TaggedTop.M() >= 200; 
        bool pass_mstop1_250_TaggedTop	  	= stop1_PtRank_TaggedTop.M() >= 250; 
        bool pass_mstop1_300_TaggedTop 	= stop1_PtRank_TaggedTop.M() >= 300; 
        bool pass_mstop1_350_TaggedTop 	= stop1_PtRank_TaggedTop.M() >= 350; 
        bool pass_mstop1_400_TaggedTop 	= stop1_PtRank_TaggedTop.M() >= 400;
        bool pass_mstop1_500_TaggedTop		= stop1_PtRank_TaggedTop.M() >= 500;

        bool pass_mstop2_100_TaggedTop	  	= stop2_PtRank_TaggedTop.M() >= 100; 
        bool pass_mstop2_200_TaggedTop	  	= stop2_PtRank_TaggedTop.M() >= 200; 
        bool pass_mstop2_250_TaggedTop	  	= stop2_PtRank_TaggedTop.M() >= 250; 
        bool pass_mstop2_300_TaggedTop 	= stop2_PtRank_TaggedTop.M() >= 300; 
        bool pass_mstop2_350_TaggedTop 	= stop2_PtRank_TaggedTop.M() >= 350; 
        bool pass_mstop2_400_TaggedTop 	= stop2_PtRank_TaggedTop.M() >= 400;
        bool pass_mstop2_500_TaggedTop		= stop2_PtRank_TaggedTop.M() >= 500;
        
        // Combined stop mass cuts

        bool pass_m12ge100_0l		= pass_mstop1_100_0l && pass_mstop2_100_0l;
        bool pass_m12ge200_0l		= pass_mstop1_200_0l && pass_mstop2_200_0l;
        bool pass_m12ge250_0l		= pass_mstop1_250_0l && pass_mstop2_250_0l;
        bool pass_m12ge300_0l		= pass_mstop1_300_0l && pass_mstop2_300_0l;
        bool pass_m12ge350_0l		= pass_mstop1_350_0l && pass_mstop2_350_0l;
        bool pass_m12ge400_0l		= pass_mstop1_400_0l && pass_mstop2_400_0l;
        bool pass_m12ge500_0l		= pass_mstop1_500_0l && pass_mstop2_500_0l;

        bool pass_m12ge100_TaggedTop		= pass_mstop1_100_TaggedTop && pass_mstop2_100_TaggedTop;
        bool pass_m12ge200_TaggedTop		= pass_mstop1_200_TaggedTop && pass_mstop2_200_TaggedTop;
        bool pass_m12ge250_TaggedTop		= pass_mstop1_250_TaggedTop && pass_mstop2_250_TaggedTop;
        bool pass_m12ge300_TaggedTop		= pass_mstop1_300_TaggedTop && pass_mstop2_300_TaggedTop;
        bool pass_m12ge350_TaggedTop		= pass_mstop1_350_TaggedTop && pass_mstop2_350_TaggedTop;
        bool pass_m12ge400_TaggedTop		= pass_mstop1_400_TaggedTop && pass_mstop2_400_TaggedTop;
        bool pass_m12ge500_TaggedTop		= pass_mstop1_500_TaggedTop && pass_mstop2_500_TaggedTop;
        
        //----------------------
        //--- Stop Vars 0l -----
        //----------------------

        double stop1_PtRank_0l_Mass 	= stop1_PtRank_0l.M();
        double stop1_PtRank_0l_Pt		= stop1_PtRank_0l.Pt();
        double stop1_PtRank_0l_Phi 	= stop1_PtRank_0l.Phi();
        double stop1_PtRank_0l_Eta 	= stop1_PtRank_0l.Eta();
        
        double stop2_PtRank_0l_Mass 	= stop2_PtRank_0l.M();
        double stop2_PtRank_0l_Pt 		= stop2_PtRank_0l.Pt();
        double stop2_PtRank_0l_Phi 	= stop2_PtRank_0l.Phi();
        double stop2_PtRank_0l_Eta 	= stop2_PtRank_0l.Eta();
            
        double dM_stop1_stop2_0l 		= abs ( stop1_PtRank_0l_Mass - stop2_PtRank_0l_Mass );
        double dR_stop1_stop2_0l		= stop1_PtRank_0l.DeltaR( stop2_PtRank_0l );
        double dPhi_stop1_stop2_0l		= stop1_PtRank_0l.DeltaPhi( stop2_PtRank_0l );
        double dEta_stop1_stop2_0l		= stop1_PtRank_0l_Eta - stop2_PtRank_0l_Eta;
        double Avg_stop1_stop2_0l		= (stop1_PtRank_0l_Mass + stop2_PtRank_0l_Mass) / 2.0;
        double s12_WeightedDiff_0l		= dM_stop1_stop2_0l / Avg_stop1_stop2_0l;

        //----------------------------
        //--- Stop Vars Tagged Top ---
        //----------------------------
         
        double stop1_PtRank_TaggedTop_Mass 	= stop1_PtRank_TaggedTop.M();
        double stop1_PtRank_TaggedTop_Pt		= stop1_PtRank_TaggedTop.Pt();
        double stop1_PtRank_TaggedTop_Phi 	= stop1_PtRank_TaggedTop.Phi();
        double stop1_PtRank_TaggedTop_Eta 	= stop1_PtRank_TaggedTop.Eta();
        
        double stop2_PtRank_TaggedTop_Mass 	= stop2_PtRank_TaggedTop.M();
        double stop2_PtRank_TaggedTop_Pt 		= stop2_PtRank_TaggedTop.Pt();
        double stop2_PtRank_TaggedTop_Phi 	= stop2_PtRank_TaggedTop.Phi();
        double stop2_PtRank_TaggedTop_Eta 	= stop2_PtRank_TaggedTop.Eta();
            
        double dM_stop1_stop2_TaggedTop 		= abs ( stop1_PtRank_TaggedTop_Mass - stop2_PtRank_TaggedTop_Mass );
        double dR_stop1_stop2_TaggedTop		= stop1_PtRank_TaggedTop.DeltaR( stop2_PtRank_TaggedTop );
        double dPhi_stop1_stop2_TaggedTop		= stop1_PtRank_TaggedTop.DeltaPhi( stop2_PtRank_TaggedTop );
        double dEta_stop1_stop2_TaggedTop		= stop1_PtRank_TaggedTop_Eta - stop2_PtRank_TaggedTop_Eta;
        double Avg_stop1_stop2_TaggedTop		= (stop1_PtRank_TaggedTop_Mass + stop2_PtRank_TaggedTop_Mass) / 2.0;
        double s12_WeightedDiff_TaggedTop		= dM_stop1_stop2_TaggedTop / Avg_stop1_stop2_TaggedTop;

        //-----------------------
        //--- Other stop cuts ---
        //-----------------------
        
        bool pass_diffle1_0l			    = s12_WeightedDiff_0l <= 1;
        bool pass_diff_0to1_avg_200to300_0l	= Avg_stop1_stop2_0l >= ( s12_WeightedDiff_0l * 100 ) + 200;	
        bool pass_diff_0to1_avg_200to350_0l	= Avg_stop1_stop2_0l >= ( s12_WeightedDiff_0l * 150 ) + 200;	
        bool pass_diff_0to1_avg_200to400_0l	= Avg_stop1_stop2_0l >= ( s12_WeightedDiff_0l * 200 ) + 200;

        bool pass_diffle1_TaggedTop			    = s12_WeightedDiff_TaggedTop <= 1;
        bool pass_diff_0to1_avg_200to300_TaggedTop	= Avg_stop1_stop2_TaggedTop >= ( s12_WeightedDiff_TaggedTop * 100 ) + 200;	
        bool pass_diff_0to1_avg_200to350_TaggedTop	= Avg_stop1_stop2_TaggedTop >= ( s12_WeightedDiff_TaggedTop * 150 ) + 200;	
        bool pass_diff_0to1_avg_200to400_TaggedTop	= Avg_stop1_stop2_TaggedTop >= ( s12_WeightedDiff_TaggedTop * 200 ) + 200;
        
        //------------------
        //--- Njets cuts ---
        //------------------
        
//        bool pass_nj7 = NJets_pt45 == 7; 
//        bool pass_nj8 = NJets_pt45 == 8; 
//        bool pass_nj9 = NJets_pt45 == 9; 
//        bool pass_nj10 = NJets_pt45 == 10; 
//        bool pass_nj11 = NJets_pt45 == 11; 
//        bool pass_nj12 = NJets_pt45 == 12; 
//        bool pass_nj13 = NJets_pt45 == 13;
//        bool pass_nj14 = NJets_pt45 == 14;
        
        //------------------
        //--- Multi Cuts ---
        //------------------
        
        bool pass_0l_Full 		= passBaseline0l && pass_ge2tops && pass_dR_bjets;
        bool pass_0l_m12ge200		= pass_0l_Full && pass_m12ge200_0l;	
        bool pass_0l_m12ge250		= pass_0l_Full && pass_m12ge250_0l;	
        bool pass_0l_m12ge300		= pass_0l_Full && pass_m12ge300_0l;	
        bool pass_0l_diff_avg300	= pass_0l_Full && pass_diff_0to1_avg_200to300_0l && pass_diffle1_0l;
        bool pass_0l_diff_avg350	= pass_0l_Full && pass_diff_0to1_avg_200to350_0l && pass_diffle1_0l;
        bool pass_0l_diff_avg400	= pass_0l_Full && pass_diff_0to1_avg_200to400_0l && pass_diffle1_0l;

        bool pass_TaggedTop_m12ge200		= pass_0l_Full && pass_m12ge200_TaggedTop;	
        bool pass_TaggedTop_m12ge250		= pass_0l_Full && pass_m12ge250_TaggedTop;	
        bool pass_TaggedTop_m12ge300		= pass_0l_Full && pass_m12ge300_TaggedTop;	
        bool pass_TaggedTop_diff_avg300	= pass_0l_Full && pass_diff_0to1_avg_200to300_TaggedTop && pass_diffle1_TaggedTop;
        bool pass_TaggedTop_diff_avg350	= pass_0l_Full && pass_diff_0to1_avg_200to350_TaggedTop && pass_diffle1_TaggedTop;
        bool pass_TaggedTop_diff_avg400	= pass_0l_Full && pass_diff_0to1_avg_200to400_TaggedTop && pass_diffle1_TaggedTop;

        // -------------------
            // --- Fill Histos ---
            // -------------------                        

        const std::map<std::string, bool> cut_map_0l
        {
            {"",				pass_0l_Full },
        //    {"m12ge100",			pass_0l_Full && pass_m12ge100_0l },
            {"m12ge200_0l",			pass_0l_Full && pass_0l_m12ge200 },
            {"m12ge250_0l",			pass_0l_Full && pass_0l_m12ge250 },
            {"m12ge300_0l",			pass_0l_Full && pass_0l_m12ge300 },
            {"m12ge200_TaggedTop",			pass_0l_Full && pass_TaggedTop_m12ge200 },
            {"m12ge250_TaggedTop",			pass_0l_Full && pass_TaggedTop_m12ge250 },
            {"m12ge300_TaggedTop",			pass_0l_Full && pass_TaggedTop_m12ge300 },
        //    {"m12ge350",			pass_0l_Full && pass_m12ge350_0l },
        //    {"m12ge400",			pass_0l_Full && pass_m12ge400_0l },
        //    {"m12ge500",			pass_0l_Full && pass_m12ge500_0l },
            {"diff_0to1_avg_200to300_0l",	pass_0l_Full && pass_diff_0to1_avg_200to300_0l && pass_diffle1_0l },
            {"diff_0to1_avg_200to350_0l",	pass_0l_Full && pass_diff_0to1_avg_200to350_0l && pass_diffle1_0l },
            {"diff_0to1_avg_200to400_0l",	pass_0l_Full && pass_diff_0to1_avg_200to400_0l && pass_diffle1_0l },
            {"diff_0to1_avg_200to300_TaggedTop",	pass_0l_Full && pass_diff_0to1_avg_200to300_TaggedTop && pass_diffle1_TaggedTop },
            {"diff_0to1_avg_200to350_TaggedTop",	pass_0l_Full && pass_diff_0to1_avg_200to350_TaggedTop && pass_diffle1_TaggedTop },
            {"diff_0to1_avg_200to400_TaggedTop",	pass_0l_Full && pass_diff_0to1_avg_200to400_TaggedTop && pass_diffle1_TaggedTop },
//            {"m12ge200_nj7", 		pass_0l_m12ge200 && pass_nj7 },	
//            {"m12ge200_nj8", 		pass_0l_m12ge200 && pass_nj8 },	
//            {"m12ge200_nj9", 		pass_0l_m12ge200 && pass_nj9 },	
//            {"m12ge200_nj10", 		pass_0l_m12ge200 && pass_nj10 },	
//            {"m12ge200_nj11", 		pass_0l_m12ge200 && pass_nj11 },	
//            {"m12ge200_nj12", 		pass_0l_m12ge200 && pass_nj12 },	
//            {"m12ge200_nj13", 		pass_0l_m12ge200 && pass_nj13 },	
//            {"m12ge200_nj14", 		pass_0l_m12ge200 && pass_nj14 },	
//            {"m12ge250_nj7", 		pass_0l_m12ge250 && pass_nj7 },	
//            {"m12ge250_nj8", 		pass_0l_m12ge250 && pass_nj8 },	
//            {"m12ge250_nj9", 		pass_0l_m12ge250 && pass_nj9 },	
//            {"m12ge250_nj10", 		pass_0l_m12ge250 && pass_nj10 },	
//            {"m12ge250_nj11", 		pass_0l_m12ge250 && pass_nj11 },	
//            {"m12ge250_nj12", 		pass_0l_m12ge250 && pass_nj12 },	
//            {"m12ge250_nj13", 		pass_0l_m12ge250 && pass_nj13 },	
//            {"m12ge250_nj14", 		pass_0l_m12ge250 && pass_nj14 },	
//            {"m12ge300_nj7", 		pass_0l_m12ge300 && pass_nj7 },	
//            {"m12ge300_nj8", 		pass_0l_m12ge300 && pass_nj8 },	
//            {"m12ge300_nj9", 		pass_0l_m12ge300 && pass_nj9 },	
//            {"m12ge300_nj10", 		pass_0l_m12ge300 && pass_nj10 },	
//            {"m12ge300_nj11", 		pass_0l_m12ge300 && pass_nj11 },	
//            {"m12ge300_nj12", 		pass_0l_m12ge300 && pass_nj12 },	
//            {"m12ge300_nj13", 		pass_0l_m12ge300 && pass_nj13 },	
//            {"m12ge300_nj14", 		pass_0l_m12ge300 && pass_nj14 },	
//            {"diff_0to1_avg_200to300_nj7",	 pass_0l_diff_avg300 && pass_nj7 },
//            {"diff_0to1_avg_200to300_nj8",	 pass_0l_diff_avg300 && pass_nj8 },
//            {"diff_0to1_avg_200to300_nj9",	 pass_0l_diff_avg300 && pass_nj9 },
//            {"diff_0to1_avg_200to300_nj10",	 pass_0l_diff_avg300 && pass_nj10 },
//            {"diff_0to1_avg_200to300_nj11",	 pass_0l_diff_avg300 && pass_nj11 },
//            {"diff_0to1_avg_200to300_nj12",	 pass_0l_diff_avg300 && pass_nj12 },
//            {"diff_0to1_avg_200to300_nj13",	 pass_0l_diff_avg300 && pass_nj13 },
//            {"diff_0to1_avg_200to300_nj14",	 pass_0l_diff_avg300 && pass_nj14 },
//            {"diff_0to1_avg_200to350_nj7",	 pass_0l_diff_avg350 && pass_nj7 },
//            {"diff_0to1_avg_200to350_nj8",	 pass_0l_diff_avg350 && pass_nj8 },
//            {"diff_0to1_avg_200to350_nj9",	 pass_0l_diff_avg350 && pass_nj9 },
//            {"diff_0to1_avg_200to350_nj10",	 pass_0l_diff_avg350 && pass_nj10 },
//            {"diff_0to1_avg_200to350_nj11",	 pass_0l_diff_avg350 && pass_nj11 },
//            {"diff_0to1_avg_200to350_nj12",	 pass_0l_diff_avg350 && pass_nj12 },
//            {"diff_0to1_avg_200to350_nj13",	 pass_0l_diff_avg350 && pass_nj13 },
//            {"diff_0to1_avg_200to350_nj14",	 pass_0l_diff_avg350 && pass_nj14 },
//            {"diff_0to1_avg_200to400_nj7",	 pass_0l_diff_avg400 && pass_nj7 },
//            {"diff_0to1_avg_200to400_nj8",	 pass_0l_diff_avg400 && pass_nj8 },
//            {"diff_0to1_avg_200to400_nj9",	 pass_0l_diff_avg400 && pass_nj9 },
//            {"diff_0to1_avg_200to400_nj10",	 pass_0l_diff_avg400 && pass_nj10 },
//            {"diff_0to1_avg_200to400_nj11",	 pass_0l_diff_avg400 && pass_nj11 },
//            {"diff_0to1_avg_200to400_nj12",	 pass_0l_diff_avg400 && pass_nj12 },
//            {"diff_0to1_avg_200to400_nj13",	 pass_0l_diff_avg400 && pass_nj13 },
//            {"diff_0to1_avg_200to400_nj14",	 pass_0l_diff_avg400 && pass_nj14 },
        
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
        
        my_histos["h_mstop_cuts_0l"] -> Fill( 0.5, eventweight );
        
        if( pass_m12ge100_0l ) {
            my_histos["h_mstop_cuts_0l"] -> Fill( 1.5, eventweight );
        }
        if( pass_m12ge200_0l ) {
            my_histos["h_mstop_cuts_0l"] -> Fill( 2.5, eventweight );
        }
        if( pass_m12ge250_0l ) {
            my_histos["h_mstop_cuts_0l"] -> Fill( 3.5, eventweight );
        }
        if( pass_m12ge300_0l ) {
            my_histos["h_mstop_cuts_0l"] -> Fill( 4.5, eventweight );
        }
        if( pass_m12ge350_0l ) {
            my_histos["h_mstop_cuts_0l"] -> Fill( 5.5, eventweight );
        }
        if( pass_m12ge400_0l ) {
            my_histos["h_mstop_cuts_0l"] -> Fill( 6.5, eventweight );
        }
        if( pass_m12ge500_0l ) {
            my_histos["h_mstop_cuts_0l"] -> Fill( 7.5, eventweight );
        }

        my_histos["h_avg_diff_cuts_0l"] -> Fill( 0.5, eventweight );

        if( pass_diffle1_0l ) {
            my_histos["h_avg_diff_cuts_0l"] -> Fill( 1.5, eventweight );
            if( pass_diff_0to1_avg_200to300_0l ) {
                my_histos["h_avg_diff_cuts_0l"] -> Fill( 2.5, eventweight ); 
                if( pass_diff_0to1_avg_200to350_0l ) {
                    my_histos["h_avg_diff_cuts_0l"] -> Fill( 3.5, eventweight );
                    if( pass_diff_0to1_avg_200to400_0l ) {
                        my_histos["h_avg_diff_cuts_0l"] -> Fill( 4.5, eventweight );
                    }
                }
            }
        }
        my_histos["h_Mass_stop1_0l"] -> Fill( stop1_PtRank_0l.M(), eventweight );
        my_histos["h_Mass_stop2_0l"] -> Fill( stop2_PtRank_0l.M(), eventweight );

        my_histos["h_mstop_cuts_TaggedTop"] -> Fill( 0.5, eventweight );
        
        if( pass_m12ge100_TaggedTop ) {
            my_histos["h_mstop_cuts_TaggedTop"] -> Fill( 1.5, eventweight );
        }
        if( pass_m12ge200_TaggedTop ) {
            my_histos["h_mstop_cuts_TaggedTop"] -> Fill( 2.5, eventweight );
        }
        if( pass_m12ge250_TaggedTop ) {
            my_histos["h_mstop_cuts_TaggedTop"] -> Fill( 3.5, eventweight );
        }
        if( pass_m12ge300_TaggedTop ) {
            my_histos["h_mstop_cuts_TaggedTop"] -> Fill( 4.5, eventweight );
        }
        if( pass_m12ge350_TaggedTop ) {
            my_histos["h_mstop_cuts_TaggedTop"] -> Fill( 5.5, eventweight );
        }
        if( pass_m12ge400_TaggedTop ) {
            my_histos["h_mstop_cuts_TaggedTop"] -> Fill( 6.5, eventweight );
        }
        if( pass_m12ge500_TaggedTop ) {
            my_histos["h_mstop_cuts_TaggedTop"] -> Fill( 7.5, eventweight );
        }

        my_histos["h_avg_diff_cuts_TaggedTop"] -> Fill( 0.5, eventweight );

        if( pass_diffle1_TaggedTop ) {
            my_histos["h_avg_diff_cuts_TaggedTop"] -> Fill( 1.5, eventweight );
            if( pass_diff_0to1_avg_200to300_TaggedTop ) {
                my_histos["h_avg_diff_cuts_TaggedTop"] -> Fill( 2.5, eventweight ); 
                if( pass_diff_0to1_avg_200to350_TaggedTop ) {
                    my_histos["h_avg_diff_cuts_TaggedTop"] -> Fill( 3.5, eventweight );
                    if( pass_diff_0to1_avg_200to400_TaggedTop ) {
                        my_histos["h_avg_diff_cuts_TaggedTop"] -> Fill( 4.5, eventweight );
                    }
                }
            }
        }
        my_histos["h_Mass_stop1_TaggedTop"] -> Fill( stop1_PtRank_0l.M(), eventweight );
        my_histos["h_Mass_stop2_TaggedTop"] -> Fill( stop2_PtRank_0l.M(), eventweight );
        
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
                my_2d_histos["h_Weighted_Avg_MassStop_0l_" + kv.first] -> Fill(s12_WeightedDiff_0l, Avg_stop1_stop2_0l, eventweight); 
                my_2d_histos["h_Weighted_Avg_MassStop_TaggedTop_" + kv.first] -> Fill(s12_WeightedDiff_TaggedTop, Avg_stop1_stop2_TaggedTop, eventweight);
                
                for(unsigned int j = 0; j < Jets.size(); j++ ) {
                    if( !(GoodJets_pt45[j]) ) continue;
                    my_histos["h_jPt_0l_" 	+ kv.first] -> Fill(Jets[j].Pt(), eventweight);
                    my_histos["h_jEta_0l_" 	+ kv.first] -> Fill(Jets[j].Eta(), eventweight);
                    my_histos["h_jPhi_0l_" 	+ kv.first] -> Fill(Jets[j].Phi(), eventweight);
                    my_histos["h_jMass_0l_"	+ kv.first] -> Fill(Jets[j].M(), eventweight);
                }

                my_histos["h_dM_stops_0l_" + kv.first] -> Fill(dM_stop1_stop2_0l, eventweight);
                my_histos["h_dR_stops_0l_" + kv.first] -> Fill(dR_stop1_stop2_0l, eventweight);
                my_histos["h_dPhi_stops_0l_" + kv.first] -> Fill(dPhi_stop1_stop2_0l, eventweight);
                my_histos["h_dEta_stops_0l_" + kv.first] -> Fill(dEta_stop1_stop2_0l, eventweight);
                my_histos["h_AvgMass_stops_0l_" + kv.first] -> Fill(Avg_stop1_stop2_0l, eventweight);
                my_histos["h_WeightedMassDiff_stops_0l_" + kv.first] -> Fill(s12_WeightedDiff_0l, eventweight);
                
                my_histos["h_dM_stops_TaggedTop_" + kv.first] -> Fill(dM_stop1_stop2_TaggedTop, eventweight);
                my_histos["h_dR_stops_TaggedTop_" + kv.first] -> Fill(dR_stop1_stop2_TaggedTop, eventweight);
                my_histos["h_dPhi_stops_TaggedTop_" + kv.first] -> Fill(dPhi_stop1_stop2_TaggedTop, eventweight);
                my_histos["h_dEta_stops_TaggedTop_" + kv.first] -> Fill(dEta_stop1_stop2_TaggedTop, eventweight);
                my_histos["h_AvgMass_stops_TaggedTop_" + kv.first] -> Fill(Avg_stop1_stop2_TaggedTop, eventweight);
                my_histos["h_WeightedMassDiff_stops_TaggedTop_" + kv.first] -> Fill(s12_WeightedDiff_TaggedTop, eventweight);
                
                my_2d_histos["h_Mass_stop1stop2_0l_" + kv.first] -> Fill( stop1_PtRank_0l.M(), stop2_PtRank_0l.M(), eventweight );
                my_2d_histos["h_Mass_stop1stop2_TaggedTop_" + kv.first] -> Fill( stop1_PtRank_TaggedTop.M(), stop2_PtRank_TaggedTop.M(), eventweight );

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
