#define AnalyzeNewSignalJets_cxx
#include "Analyzer/Analyzer/include/AnalyzeNewSignalJets.h"
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

AnalyzeNewSignalJets::AnalyzeNewSignalJets()
{
    InitHistos();
}

//Define all your histograms here. 
void AnalyzeNewSignalJets::InitHistos()
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    //This event counter histogram is necessary so that we know that all the condor jobs ran successfully. If not, when you use the hadder script, you will see a discrepancy in red as the files are being hadded.
    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ;

    //1D histograms of njet distributions for two different categories of NonIsoMuons 
    my_histos.emplace ( "h_passBaseline", std::make_shared<TH1D>( "h_passBaseline", "h_passBaseline", 1, 0, 1) );
    my_histos.emplace( "h_neutralinos", std::make_shared<TH1D>( "h_neutralinos", "h_neutralinos", 1, 0, 1) );
    my_histos.emplace( "h_neutralinoDaughters", std::make_shared<TH1D>( "h_neutralinoDaughters", "h_neutralinoDaughters", 17, -8, 8) );
    
    my_histos.emplace( "h_jetDeltaR_RPV_mStop-300", std::make_shared<TH1D>( "h_jetDeltaR_RPV_mStop-300", "h_jetDeltaR_RPV_mStop-300", 40, 0, 4) );
    my_histos.emplace( "h_jetDeltaR_RPV_mStop-350", std::make_shared<TH1D>( "h_jetDeltaR_RPV_mStop-350", "h_jetDeltaR_RPV_mStop-350", 40, 0, 4) );
    my_histos.emplace( "h_jetDeltaR_RPV_mStop-800", std::make_shared<TH1D>( "h_jetDeltaR_RPV_mStop-800", "h_jetDeltaR_RPV_mStop-800", 40, 0, 4) );
    my_histos.emplace( "h_jetDeltaR_RPV_mStop-850", std::make_shared<TH1D>( "h_jetDeltaR_RPV_mStop-850", "h_jetDeltaR_RPV_mStop-850", 40, 0, 4) );
    my_histos.emplace( "h_jetDeltaR_RPV_mStop-1300", std::make_shared<TH1D>( "h_jetDeltaR_RPV_mStop-1300", "h_jetDeltaR_RPV_mStop-1300", 40, 0, 4) );
    my_histos.emplace( "h_jetDeltaR_RPV_mStop-1350", std::make_shared<TH1D>( "h_jetDeltaR_RPV_mStop-1350", "h_jetDeltaR_RPV_mStop-1350", 40, 0, 4) );

     my_histos.emplace( "h_jetDeltaR_SHH_mStop-300", std::make_shared<TH1D>( "h_jetDeltaR_SHH_mStop-300", "h_jetDeltaR_SHH_mStop-300", 40, 0, 4) );
    my_histos.emplace( "h_jetDeltaR_SHH_mStop-350", std::make_shared<TH1D>( "h_jetDeltaR_SHH_mStop-350", "h_jetDeltaR_SHH_mStop-350", 40, 0, 4) );
    my_histos.emplace( "h_jetDeltaR_SHH_mStop-800", std::make_shared<TH1D>( "h_jetDeltaR_SHH_mStop-800", "h_jetDeltaR_SHH_mStop-800", 40, 0, 4) );
    my_histos.emplace( "h_jetDeltaR_SHH_mStop-850", std::make_shared<TH1D>( "h_jetDeltaR_SHH_mStop-850", "h_jetDeltaR_SHH_mStop-850", 40, 0, 4) );
    my_histos.emplace( "h_jetDeltaR_SHH_mStop-1300", std::make_shared<TH1D>( "h_jetDeltaR_SHH_mStop-1300", "h_jetDeltaR_SHH_mStop-1300", 40, 0, 4) );
    my_histos.emplace( "h_jetDeltaR_SHH_mStop-1350", std::make_shared<TH1D>( "h_jetDeltaR_SHH_mStop-1350", "h_jetDeltaR_SHH_mStop-1350", 40, 0, 4) );

     my_histos.emplace( "h_jetDeltaR_SYY_mStop-300", std::make_shared<TH1D>( "h_jetDeltaR_SYY_mStop-300", "h_jetDeltaR_SYY_mStop-300", 40, 0, 4) );
    my_histos.emplace( "h_jetDeltaR_SYY_mStop-350", std::make_shared<TH1D>( "h_jetDeltaR_SYY_mStop-350", "h_jetDeltaR_SYY_mStop-350", 40, 0, 4) );
    my_histos.emplace( "h_jetDeltaR_SYY_mStop-800", std::make_shared<TH1D>( "h_jetDeltaR_SYY_mStop-800", "h_jetDeltaR_SYY_mStop-800", 40, 0, 4) );
    my_histos.emplace( "h_jetDeltaR_SYY_mStop-850", std::make_shared<TH1D>( "h_jetDeltaR_SYY_mStop-850", "h_jetDeltaR_SYY_mStop-850", 40, 0, 4) );
    my_histos.emplace( "h_jetDeltaR_SYY_mStop-1300", std::make_shared<TH1D>( "h_jetDeltaR_SYY_mStop-1300", "h_jetDeltaR_SYY_mStop-1300", 40, 0, 4) );
    my_histos.emplace( "h_jetDeltaR_SYY_mStop-1350", std::make_shared<TH1D>( "h_jetDeltaR_SYY_mStop-1350", "h_jetDeltaR_SYY_mStop-1350", 40, 0, 4) );
   
}

//Put everything you want to do per event here.
void AnalyzeNewSignalJets::Loop(NTupleReader& tr, double, int maxevents, bool)
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

       	const auto& GenParticles        = tr.getVec<TLorentzVector>("GenParticles");
       	const auto& GenParticlesPDG     = tr.getVec<int>("GenParticles_PdgId");
       	const auto& GenParticlesParent  = tr.getVec<int>("GenParticles_ParentId");
	const auto& GenParticlesParentIdx = tr.getVec<int>("GenParticles_ParentIdx");

	const auto& passBaseline        = tr.getVar<bool>("passBaseline1l_Good");
        const auto& passMadHT           = tr.getVar<bool>("passMadHT");

	const auto& GoodLeptons		= tr.getVec<std::pair<std::string, TLorentzVector>>("GoodLeptons");	
	const auto& NGoodLeptons	= tr.getVar<int>("NGoodLeptons");
	
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

	my_histos["h_passBaseline"]->Fill( 1, weight );

	TLorentzVector first;
	int firstMomidx = 0;

	std::string suffix = filetag.substr( filetag.find("-") + 1 );
	
	//Find jets from neutralino decay
	for( unsigned int gpi = 0; gpi < GenParticles.size(); gpi++ )
	  {
	    TLorentzVector lv( GenParticles.at(gpi) );
	    int pdgid = abs( GenParticlesPDG.at(gpi) ) ;

      	    int momid = abs( GenParticlesParent.at(gpi) );
	    int momidx = GenParticlesParentIdx.at(gpi);
	    
	    
				   
	    if ( pdgid == 1000022 && filetag.find("StealthSHH") != std::string::npos )
	      {
		my_histos["h_neutralinos"]->Fill( 0.5, weight );
	      }

	    if ( momid == 1000022 )
	      {
		my_histos["h_neutralinoDaughters"]->Fill( pdgid, weight );
		
	      }
	    if ( momid == 1000022 && first.Pt() == 0 )
	      {
		first = lv;
		firstMomidx = momidx;
	      }
	    else if ( momid == 1000022 && momidx == firstMomidx )
	      {
		if ( filetag.find("RPV") != std::string::npos )
		  my_histos["h_jetDeltaR_RPV_mStop-" + suffix]->Fill( lv.DeltaR( first ), weight );
		if ( filetag.find("StealthSHH") != std::string::npos )
		  my_histos["h_jetDeltaR_SHH_mStop-" + suffix]->Fill( lv.DeltaR( first ), weight );
		if ( filetag.find("StealthSYY") != std::string::npos )
		  my_histos["h_jetDeltaR_SYY_mStop-" + suffix]->Fill( lv.DeltaR( first ), weight );
	      }
	    
	  }
    }//end while loop
}

void AnalyzeNewSignalJets::WriteHistos(TFile* outfile)
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
