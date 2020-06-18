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

Analyze0LepCutFlow::Analyze0LepCutFlow()
{
    InitHistos();
}

void Analyze0LepCutFlow::InitHistos()
{
    TH1::SetDefaultSumw2();

    // Declare all your histograms here, that way we can fill them for multiple chains
    my_histos.emplace("h_event_sel_cutFlow",	std::make_shared<TH1D>("h_event_sel_cutFlow",	"h_event_sel_cutFlow",	10,	0,	10));
    my_histos.emplace("h_ntops",    std::make_shared<TH1D>("h_ntops",  "h_ntops",   10,  0,     10  ) );
    my_histos.emplace("h_ntops_j1", std::make_shared<TH1D>("h_ntops_j1","h_ntops_j1", 10,  0,     10  ) );
    my_histos.emplace("h_ntops_j2", std::make_shared<TH1D>("h_ntops_j2","h_ntops_j2", 10,  0,     10  ) );
    my_histos.emplace("h_ntops_j3", std::make_shared<TH1D>("h_ntops_j3","h_ntops_j3", 10,  0,     10  ) );
    my_histos.emplace("h_ntops_all", std::make_shared<TH1D>("h_tops_all","h_tops_all;ntops;All/merged/between/resolved", 4, 0, 4));

    my_histos.emplace("h_1j_vs_3j", std::make_shared<TH1D>("h_1j_vs_3j", "h_1j_vs_3j", 3, 0, 3 ) );

    my_efficiencies.emplace("event_sel_total", std::make_shared<TEfficiency>("event_sel_total","Total event selection efficiency;Cut;#epsilon",8,0,8));

}

void Analyze0LepCutFlow::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    while( tr.getNextEvent() )
    {
        const auto& HT_trigger         = tr.getVar<double>("HT_trigger");
        const auto& ntops              = tr.getVar<int>("ntops");
        const auto& ntops_3jet         = tr.getVar<int>("ntops_3jet");
        const auto& ntops_2jet         = tr.getVar<int>("ntops_2jet");
        const auto& ntops_1jet         = tr.getVar<int>("ntops_1jet");
        const auto& runtype            = tr.getVar<std::string>("runtype");     
        const auto& NJets_pt45         = tr.getVar<int>("NGoodJets_pt45");
        const auto& NBJets_pt45        = tr.getVar<int>("NGoodBJets_pt45");
        const auto& NGoodLeptons       = tr.getVar<int>("NGoodLeptons");
	const auto& passBaseline0l     = tr.getVar<bool>("passBaseline0l_Good");
	const auto& dR_bjets	       = tr.getVar<double>("dR_bjets");

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
        bool pass_0l              = NGoodLeptons==0;
        bool pass_njet_pt45       = NJets_pt45>=6;
        bool pass_HT_trigger      = HT_trigger > 500;
        bool pass_njet_pt45_2btag = NBJets_pt45 >= 2;
	bool pass_ge2tops	  = ntops >= 2;
	bool pass_dR_bjets	  = dR_bjets >= 1;
        bool pass_2tops 	  = ntops == 2;
 	bool pass_2_1j		  = ntops_1jet == 2;
	bool pass_2_3j		  = ntops_3jet == 2;

	// -------------------
        // --- Fill Histos ---
        // -------------------                        

        // Initialize Histograms
        if(!initHistos)
        {
            InitHistos();
            initHistos = true;
        }


	my_histos["h_event_sel_cutFlow"]->Fill(0.5, eventweight);
	if( pass_0l ){
		my_histos["h_event_sel_cutFlow"]->Fill(1.5 ,eventweight);
		if( pass_HT_trigger ){
			my_histos["h_event_sel_cutFlow"]->Fill(2.5, eventweight);
			if( pass_njet_pt45 ){
				my_histos["h_event_sel_cutFlow"]->Fill(3.5, eventweight);
				if( pass_njet_pt45_2btag ){
					my_histos["h_event_sel_cutFlow"]->Fill(4.5, eventweight);
					if( pass_ge2tops ){
						my_histos["h_event_sel_cutFlow"]->Fill(5.5, eventweight);
						if( pass_dR_bjets ){
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
        
	if( !( passBaseline0l && pass_dR_bjets ) ) continue;

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
