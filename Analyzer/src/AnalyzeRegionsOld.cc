#define AnalyzeRegionsOld_cxx 
#include "Analyzer/Analyzer/include/AnalyzeRegionsOld.h" 
#include "Framework/Framework/include/Utility.h" 
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "fastjet/ClusterSequence.hh"

#include <TH1D.h> 
#include <TH2D.h>
#include <TStyle.h> 
#include <TCanvas.h> 
#include <TEfficiency.h> 
#include <TRandom3.h> 
#include <iostream> 
#include <TFile.h>

AnalyzeRegionsOld::AnalyzeRegionsOld()
{
	InitHistos();
}

//Define all your histograms here. 
void AnalyzeRegionsOld::InitHistos() 
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    //This event counter histogram is necessary so that we know that all the condor jobs ran successfully. If not, when you use the hadder script, you will see a discrepancy in red as the files are being hadded.
    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ;

    //Define 1D histograms
    my_histos.emplace( "h_nGoodLeps_CR1", std::make_shared<TH1D>( "h_nGoodLeps_CR1", "h_nGoodLeps_CR1", 10, 0, 10 ) ) ;
    my_histos.emplace( "h_nGoodLeps_CR2", std::make_shared<TH1D>( "h_nGoodLeps_CR2", "h_nGoodLeps_CR2", 10, 0, 10 ) ) ;
    my_histos.emplace( "h_nGoodLeps_CR3", std::make_shared<TH1D>( "h_nGoodLeps_CR3", "h_nGoodLeps_CR3", 10, 0, 10 ) ) ;
    my_histos.emplace( "h_nGoodLeps_SR1", std::make_shared<TH1D>( "h_nGoodLeps_SR1", "h_nGoodLeps_SR1", 10, 0, 10 ) ) ;
    my_histos.emplace( "h_nGoodLeps_SR2", std::make_shared<TH1D>( "h_nGoodLeps_SR2", "h_nGoodLeps_SR2", 10, 0, 10 ) ) ;
    my_histos.emplace( "h_nGoodLeps_SR3", std::make_shared<TH1D>( "h_nGoodLeps_SR3", "h_nGoodLeps_SR3", 10, 0, 10 ) ) ;
    my_histos.emplace( "h_nGoodLeps_SR4", std::make_shared<TH1D>( "h_nGoodLeps_SR4", "h_nGoodLeps_SR4", 10, 0, 10 ) ) ;
    my_histos.emplace( "h_nGoodLeps_SR5", std::make_shared<TH1D>( "h_nGoodLeps_SR5", "h_nGoodLeps_SR5", 10, 0, 10 ) ) ;
    my_histos.emplace( "h_nGoodJets_CR1", std::make_shared<TH1D>( "h_nGoodJets_CR1", "h_nGoodJets_CR1", 30, 0, 30 ) ) ;
    my_histos.emplace( "h_nGoodJets_CR2", std::make_shared<TH1D>( "h_nGoodJets_CR2", "h_nGoodJets_CR2", 30, 0, 30 ) ) ;
    my_histos.emplace( "h_nGoodJets_CR3", std::make_shared<TH1D>( "h_nGoodJets_CR3", "h_nGoodJets_CR3", 30, 0, 30 ) ) ;
    my_histos.emplace( "h_nGoodJets_SR1", std::make_shared<TH1D>( "h_nGoodJets_SR1", "h_nGoodJets_SR1", 30, 0, 30 ) ) ;
    my_histos.emplace( "h_nGoodJets_SR2", std::make_shared<TH1D>( "h_nGoodJets_SR2", "h_nGoodJets_SR2", 30, 0, 30 ) ) ;
    my_histos.emplace( "h_nGoodJets_SR3", std::make_shared<TH1D>( "h_nGoodJets_SR3", "h_nGoodJets_SR3", 30, 0, 30 ) ) ;
    my_histos.emplace( "h_nGoodJets_SR4", std::make_shared<TH1D>( "h_nGoodJets_SR4", "h_nGoodJets_SR4", 30, 0, 30 ) ) ;
    my_histos.emplace( "h_nGoodJets_SR5", std::make_shared<TH1D>( "h_nGoodJets_SR5", "h_nGoodJets_SR5", 30, 0, 30 ) ) ;
    my_histos.emplace( "h_nGoodbJets_CR1", std::make_shared<TH1D>( "h_nGoodbJets_CR1", "h_nGoodbJets_CR1", 15, 0, 15 ) ) ;
    my_histos.emplace( "h_nGoodbJets_CR2", std::make_shared<TH1D>( "h_nGoodbJets_CR2", "h_nGoodbJets_CR2", 15, 0, 15 ) ) ;
    my_histos.emplace( "h_nGoodbJets_CR3", std::make_shared<TH1D>( "h_nGoodbJets_CR3", "h_nGoodbJets_CR3", 15, 0, 15 ) ) ;
    my_histos.emplace( "h_nGoodbJets_SR1", std::make_shared<TH1D>( "h_nGoodbJets_SR1", "h_nGoodbJets_SR1", 15, 0, 15 ) ) ;
    my_histos.emplace( "h_nGoodbJets_SR2", std::make_shared<TH1D>( "h_nGoodbJets_SR2", "h_nGoodbJets_SR2", 15, 0, 15 ) ) ;
    my_histos.emplace( "h_nGoodbJets_SR3", std::make_shared<TH1D>( "h_nGoodbJets_SR3", "h_nGoodbJets_SR3", 15, 0, 15 ) ) ;
    my_histos.emplace( "h_nGoodbJets_SR4", std::make_shared<TH1D>( "h_nGoodbJets_SR4", "h_nGoodbJets_SR4", 15, 0, 15 ) ) ;
    my_histos.emplace( "h_nGoodbJets_SR5", std::make_shared<TH1D>( "h_nGoodbJets_SR5", "h_nGoodbJets_SR5", 15, 0, 15 ) ) ;
    my_histos.emplace( "h_HT_CR1", std::make_shared<TH1D>( "h_HT_CR1", "h_HT_CR1", 100, 0, 3000 ) ) ;
    my_histos.emplace( "h_HT_CR2", std::make_shared<TH1D>( "h_HT_CR2", "h_HT_CR2", 100, 0, 3000 ) ) ;
    my_histos.emplace( "h_HT_CR3", std::make_shared<TH1D>( "h_HT_CR3", "h_HT_CR3", 100, 0, 3000 ) ) ;
    my_histos.emplace( "h_HT_SR1", std::make_shared<TH1D>( "h_HT_SR1", "h_HT_SR1", 100, 0, 3000 ) ) ;
    my_histos.emplace( "h_HT_SR2", std::make_shared<TH1D>( "h_HT_SR2", "h_HT_SR2", 100, 0, 3000 ) ) ;
    my_histos.emplace( "h_HT_SR3", std::make_shared<TH1D>( "h_HT_SR3", "h_HT_SR3", 100, 0, 3000 ) ) ;
    my_histos.emplace( "h_HT_SR4", std::make_shared<TH1D>( "h_HT_SR4", "h_HT_SR4", 100, 0, 3000 ) ) ;
    my_histos.emplace( "h_HT_SR5", std::make_shared<TH1D>( "h_HT_SR5", "h_HT_SR5", 100, 0, 3000 ) ) ;
    my_histos.emplace( "h_MJ_CR1", std::make_shared<TH1D>( "h_MJ_CR1", "h_MJ_CR1", 50, 0, 2000 ) ) ;
    my_histos.emplace( "h_MJ_CR2", std::make_shared<TH1D>( "h_MJ_CR2", "h_MJ_CR2", 50, 0, 2000 ) ) ;
    my_histos.emplace( "h_MJ_CR3", std::make_shared<TH1D>( "h_MJ_CR3", "h_MJ_CR3", 50, 0, 2000 ) ) ;
    my_histos.emplace( "h_MJ_SR1", std::make_shared<TH1D>( "h_MJ_SR1", "h_MJ_SR1", 50, 0, 2000 ) ) ;
    my_histos.emplace( "h_MJ_SR2", std::make_shared<TH1D>( "h_MJ_SR2", "h_MJ_SR2", 50, 0, 2000 ) ) ;
    my_histos.emplace( "h_MJ_SR3", std::make_shared<TH1D>( "h_MJ_SR3", "h_MJ_SR3", 50, 0, 2000 ) ) ;
    my_histos.emplace( "h_MJ_SR4", std::make_shared<TH1D>( "h_MJ_SR4", "h_MJ_SR4", 50, 0, 2000 ) ) ;
    my_histos.emplace( "h_MJ_SR5", std::make_shared<TH1D>( "h_MJ_SR5", "h_MJ_SR5", 50, 0, 2000 ) ) ;
    my_histos.emplace( "h_largeRJets_CR1", std::make_shared<TH1D>( "h_largeRJets_CR1", "h_largeRJets_CR1", 14, 0, 14 ) ) ;
    my_histos.emplace( "h_largeRJets_CR2", std::make_shared<TH1D>( "h_largeRJets_CR2", "h_largeRJets_CR2", 14, 0, 14 ) ) ;
    my_histos.emplace( "h_largeRJets_CR3", std::make_shared<TH1D>( "h_largeRJets_CR3", "h_largeRJets_CR3", 14, 0, 14 ) ) ;
    my_histos.emplace( "h_largeRJets_SR1", std::make_shared<TH1D>( "h_largeRJets_SR1", "h_largeRJets_SR1", 14, 0, 14 ) ) ;
    my_histos.emplace( "h_largeRJets_SR2", std::make_shared<TH1D>( "h_largeRJets_SR2", "h_largeRJets_SR2", 14, 0, 14 ) ) ;
    my_histos.emplace( "h_largeRJets_SR3", std::make_shared<TH1D>( "h_largeRJets_SR3", "h_largeRJets_SR3", 14, 0, 14 ) ) ;
    my_histos.emplace( "h_largeRJets_SR4", std::make_shared<TH1D>( "h_largeRJets_SR4", "h_largeRJets_SR4", 14, 0, 14 ) ) ;
    my_histos.emplace( "h_largeRJets_SR5", std::make_shared<TH1D>( "h_largeRJets_SR5", "h_largeRJets_SR5", 14, 0, 14 ) ) ;
    my_histos.emplace( "Events_by_Region", std::make_shared<TH1D>( "Events_by_Region", "Events_by_Region", 8, 0, 8 ) ) ;
}

//Put everything you want to do per event here. 
void AnalyzeRegionsOld::Loop(NTupleReader& tr, double, int maxevents, bool) 
{
    while( tr.getNextEvent() )
    {
        //This is added to count the number of events- do not change the next two lines.
        const auto& eventCounter = tr.getVar<int>("eventCounter");
        my_histos["EventCounter"]->Fill( eventCounter );

        //--------------------------------------------------
        //-- Print Event Number
        //--------------------------------------------------
        
        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() & (10000 == 0) ) printf( " Event %i\n", tr.getEvtNum() );
        
        const auto& runtype = tr.getVar<std::string>("runtype");

        const auto& Jets_ID = tr.getVec<bool>("Jets_ID");
	const auto& allElectrons = tr.getVec<TLorentzVector>("Electrons");
	const auto& allMuons = tr.getVec<TLorentzVector>("Muons");
	const auto& Jets = tr.getVec<TLorentzVector>("Jets");
        const auto& allElectrons_passIso = tr.getVec<bool>("Electrons_passIso");
	const auto& allElectrons_tightID = tr.getVec<bool>("Electrons_tightID");
	const auto& allMuons_passIso = tr.getVec<bool>("Muons_passIso");
        const auto& allMuons_tightID = tr.getVec<bool>("Muons_tightID");
	const auto& Jets_bDiscriminatorCSV = tr.getVec<double>("Jets_bDiscriminatorCSV");
        const auto& TriggerNames         = tr.getVec<std::string>("TriggerNames");
	const auto& TriggerPass          = tr.getVec<int>("TriggerPass");

	std::vector<bool> GoodElectrons( allElectrons.size(), false );
	std::vector<bool> GoodMuons( allMuons.size() , false );
        std::vector<bool> GoodLeptons( ( allElectrons.size() + allMuons.size() ) , false );
	std::vector<bool> goodJets( Jets.size(), false );
	std::vector<bool> goodbJets( Jets.size(), false );
	std::vector<fastjet::PseudoJet> particles;	

	bool OneGoodLepton = false;
	bool fourplusJets = false;
        bool oneplusbJets = false;
        bool passHT = false;
	bool passTriggers = false;
	bool passMJ = false;

	bool CRone = false;
	bool CRtwo = false;
	bool CRthree = false;
	bool SRone = false;
	bool SRtwo = false;
	bool SRthree = false;
	bool SRfour = false;
	bool SRfive = false;

	bool allcut = false;

	const auto& passMadHT = tr.getVar<bool>("passMadHT");
	
        if(maxevents != -1 && tr.getEvtNum() >= maxevents) break;
        if ( tr.getEvtNum() % 1000 == 0 ) printf(" Event %i\n", tr.getEvtNum() ) ;
	
        // ------------------------
        // -- Define weight
        // ------------------------
        double weight = 1.0;
        double eventweight = 1.0;
        double leptonScaleFactor = 1.0;
        double prefiringScaleFactor = 1.0;
        double puScaleFactor = 1.0;
        
        if(runtype == "MC")
        {
            if( !passMadHT ) continue; //Make sure not to double count DY events
            // Define Lumi weight
            const auto& Weight = tr.getVar<double>("Weight");
            const auto& lumi = tr.getVar<double>("Lumi");
            eventweight = lumi*Weight;
            
            //PileupWeight = tr.getVar<double>("_PUweightFactor");
            prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor");
            puScaleFactor = tr.getVar<double>("puWeightCorr");
            
            weight *= eventweight*leptonScaleFactor*prefiringScaleFactor*puScaleFactor;
        }
        
        //Make cuts and fill histograms here
	int lepgood = 0;

	for(unsigned int j = 0; j < allElectrons.size(); j++)
	{
		if( allElectrons.at(j).Pt() > 20 and abs( allElectrons.at(j).Eta() ) < 2.5 and allElectrons_passIso.at(j) and allElectrons_tightID.at(j) )
		{
			GoodElectrons[j] = true;
			GoodLeptons[j] = true;
			lepgood++;
		}
	}

	for(unsigned int j = 0; j < allMuons.size(); j++)
	{
                if( allMuons.at(j).Pt() > 20 and abs( allMuons.at(j).Eta() ) < 2.4 and allMuons_passIso.at(j) and allMuons_tightID.at(j) )
		{
                	GoodMuons[j] = true;
                	GoodLeptons[(allElectrons.size() +j)] = true;
			lepgood++;
		}
        }

	if( lepgood == 1 ) OneGoodLepton = true;
	
	int Jetsgood = 0;
	int bJetsgood = 0;
	double eventHT = 0;
	
	for(unsigned int j = 0; j < Jets.size(); j++)
	{
		if( Jets.at(j).Pt() > 30 and abs( Jets.at(j).Eta() ) <= 2.4 and Jets_ID.at(j) ) goodJets[j] = true;
	} 

	for(unsigned int es = 0; es < allElectrons.size(); es++)
	{
		if( GoodElectrons[es] == false ) continue;
		double tempdR = 1.0;
		int jClose = 0;
		TLorentzVector myElec = allElectrons.at(es);
		for(unsigned int j = 0; j < Jets.size(); j++)
		{
			TLorentzVector myJet = Jets.at(j);
			double JetdR = myElec.DeltaR(myJet);
			if( JetdR < tempdR )
			{
				tempdR = JetdR;
				jClose = j;
			}
		}		
		if( tempdR < 0.4 ) goodJets[jClose] = false;
	}
		
	for(unsigned int mus = 0; mus < allMuons.size(); mus++)
	{
		if( GoodMuons[mus] == false ) continue;
		double tempdR = 1.0;
                int jClose = 0;
                TLorentzVector myMuon = allMuons.at(mus);
		for(unsigned int j = 0; j < Jets.size(); j++)
		{
			TLorentzVector myJet = Jets.at(j);
			double JetdR = myMuon.DeltaR(myJet);
			if( JetdR < tempdR )
			{
				tempdR = JetdR;
				jClose = j;
			}
		}
		if( tempdR < 0.4 ) goodJets[jClose] = false;
        }			
	
	for( unsigned int j = 0; j < Jets.size(); j++)
	{
		if( goodJets[j] )
		{
			Jetsgood++;
        		eventHT = eventHT + Jets.at(j).Pt();
        		if( Jets_bDiscriminatorCSV.at(j) > 0.8484 )
        		{
          			bJetsgood++;
                		goodbJets[j] = true;
        		}
			particles.push_back( fastjet::PseudoJet( Jets.at(j).Px() , Jets.at(j).Py() , Jets.at(j).Pz() , Jets.at(j).E() ) ); 
		}
	}
	
	std::cout << "Jetsgood " << Jetsgood << std::endl;	

	if( Jetsgood > 3 ) fourplusJets = true;
	if( bJetsgood > 0 ) oneplusbJets = true;
	if( eventHT > 1200 ) passHT = true;
	
	int tone = 0;
	int ttwo = 0;
	for(unsigned int j = 0; j < TriggerNames.size(); j++)
	{
		if( TriggerNames.at(j) == "HLT_PFHT900_v" ) tone = j;
		if( TriggerNames.at(j) == "HLT_PFJet450_v" ) ttwo = j;
	} 
	if( TriggerPass[tone] == 1 or TriggerPass[ttwo] == 1) passTriggers = true;	
	
	if( OneGoodLepton )
        {
                for(unsigned int j = 0; j < GoodLeptons.size(); j++)
                {
                        if( GoodLeptons[j] )
                        {
                                if( j >= allElectrons.size() )
                                {
                                        unsigned int k = j - allElectrons.size();
                                        particles.push_back( fastjet::PseudoJet( allMuons.at(k).Px() , allMuons.at(k).Py() , allMuons.at(k).Pz() , allMuons.at(k).E() ) );
                                }
                                else
                                {
                                        particles.push_back( fastjet::PseudoJet( allElectrons.at(j).Px() , allElectrons.at(j).Py() , allElectrons.at(j).Pz() , allElectrons.at(j).E() ) );
                                }
                        }
                }
        }	
        
	double R = 1.2;
	fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
	fastjet::ClusterSequence cs(particles, jet_def);
	std::vector<fastjet::PseudoJet> largeRjets = fastjet::sorted_by_pt(cs.inclusive_jets());  
		
	double MJ = 0.0;
	for(unsigned int j = 0; j < largeRjets.size(); j++)
	{
		MJ = MJ + largeRjets.at(j).m();
	}
	
	if( MJ > 500 ) passMJ = true;	
	
	if( OneGoodLepton and fourplusJets and oneplusbJets and passHT and passMJ and passTriggers ) allcut = true;

	if( MJ > 500 and MJ <= 800 and allcut )
	{
		if( Jetsgood == 4 || Jetsgood == 5 ) CRone = true;
		if( Jetsgood == 6 || Jetsgood == 7 ) CRtwo = true;
		if( Jetsgood >= 8 ) SRone = true;
	}
	if( MJ > 800 and  MJ <= 1000 and allcut )
	{
		if( Jetsgood == 4 || Jetsgood == 5 ) CRthree = true;
		if( Jetsgood == 6 || Jetsgood == 7 ) SRtwo = true;
		if( Jetsgood >=	8 ) SRthree = true;
	}
	if( MJ > 1000 and allcut )
	{
		if( Jetsgood == 4 || Jetsgood == 5 ) CRthree = true;
                if( Jetsgood == 6 || Jetsgood == 7 ) SRfour = true;
                if( Jetsgood >= 8 ) SRfive = true;
	}
	
	if( CRone )
	{
		my_histos["h_nGoodLeps_CR1"]->Fill( lepgood , weight );
		my_histos["h_nGoodJets_CR1"]->Fill( Jetsgood , weight );
		my_histos["h_nGoodbJets_CR1"]->Fill( bJetsgood , weight );
		my_histos["h_HT_CR1"]->Fill( eventHT , weight );
		my_histos["h_MJ_CR1"]->Fill( MJ , weight );
		my_histos["h_largeRJets_CR1"]->Fill( largeRjets.size() , weight );
		my_histos["Events_by_Region"]->Fill( 0.0 , weight );
	}

	if( CRtwo )
        {
                my_histos["h_nGoodLeps_CR2"]->Fill( lepgood , weight );
                my_histos["h_nGoodJets_CR2"]->Fill( Jetsgood , weight );
                my_histos["h_nGoodbJets_CR2"]->Fill( bJetsgood , weight );
                my_histos["h_HT_CR2"]->Fill( eventHT , weight );
                my_histos["h_MJ_CR2"]->Fill( MJ , weight );
                my_histos["h_largeRJets_CR2"]->Fill( largeRjets.size() , weight );
                my_histos["Events_by_Region"]->Fill( 1.0 , weight );
        }

	if( CRthree )
	{
		my_histos["h_nGoodLeps_CR3"]->Fill( lepgood , weight );
                my_histos["h_nGoodJets_CR3"]->Fill( Jetsgood , weight );
                my_histos["h_nGoodbJets_CR3"]->Fill( bJetsgood , weight );
                my_histos["h_HT_CR3"]->Fill( eventHT , weight );
                my_histos["h_MJ_CR3"]->Fill( MJ , weight );
                my_histos["h_largeRJets_CR3"]->Fill( largeRjets.size() , weight );
                my_histos["Events_by_Region"]->Fill( 2.0 , weight );
	}
		
	if( SRone )
        {
                my_histos["h_nGoodLeps_SR1"]->Fill( lepgood , weight );
                my_histos["h_nGoodJets_SR1"]->Fill( Jetsgood , weight );
                my_histos["h_nGoodbJets_SR1"]->Fill( bJetsgood , weight );
                my_histos["h_HT_SR1"]->Fill( eventHT , weight );
                my_histos["h_MJ_SR1"]->Fill( MJ , weight );
                my_histos["h_largeRJets_SR1"]->Fill( largeRjets.size() , weight );
                my_histos["Events_by_Region"]->Fill( 3.0 , weight );
	}

	if( SRtwo )
        {
                my_histos["h_nGoodLeps_SR2"]->Fill( lepgood , weight );
                my_histos["h_nGoodJets_SR2"]->Fill( Jetsgood , weight );
                my_histos["h_nGoodbJets_SR2"]->Fill( bJetsgood , weight );
                my_histos["h_HT_SR2"]->Fill( eventHT , weight );
                my_histos["h_MJ_SR2"]->Fill( MJ , weight );
                my_histos["h_largeRJets_SR2"]->Fill( largeRjets.size() , weight );
                my_histos["Events_by_Region"]->Fill( 4.0 , weight );
        }
	
	if( SRthree )
        {
                my_histos["h_nGoodLeps_SR3"]->Fill( lepgood , weight );
                my_histos["h_nGoodJets_SR3"]->Fill( Jetsgood , weight );
                my_histos["h_nGoodbJets_SR3"]->Fill( bJetsgood , weight );
                my_histos["h_HT_SR3"]->Fill( eventHT , weight );
                my_histos["h_MJ_SR3"]->Fill( MJ , weight );
                my_histos["h_largeRJets_SR3"]->Fill( largeRjets.size() , weight );
                my_histos["Events_by_Region"]->Fill( 5.0 , weight );
        }	

	if( SRfour )
        {
                my_histos["h_nGoodLeps_SR4"]->Fill( lepgood , weight );
                my_histos["h_nGoodJets_SR4"]->Fill( Jetsgood , weight );
                my_histos["h_nGoodbJets_SR4"]->Fill( bJetsgood , weight );
                my_histos["h_HT_SR4"]->Fill( eventHT , weight );
                my_histos["h_MJ_SR4"]->Fill( MJ , weight );
                my_histos["h_largeRJets_SR4"]->Fill( largeRjets.size() , weight );
                my_histos["Events_by_Region"]->Fill( 6.0 , weight );
        }

	if( SRfive )
        {
                my_histos["h_nGoodLeps_SR5"]->Fill( lepgood , weight );
                my_histos["h_nGoodJets_SR5"]->Fill( Jetsgood , weight );
                my_histos["h_nGoodbJets_SR5"]->Fill( bJetsgood , weight );
                my_histos["h_HT_SR5"]->Fill( eventHT , weight );
                my_histos["h_MJ_SR5"]->Fill( MJ , weight );
                my_histos["h_largeRJets_SR5"]->Fill( largeRjets.size() , weight );
                my_histos["Events_by_Region"]->Fill( 7.0 , weight );
        }
    } 
}

void AnalyzeRegionsOld::WriteHistos(TFile* outfile) 
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

