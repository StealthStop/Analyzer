#define AnalyzeWorkspace_cxx 
#include "Analyzer/Analyzer/include/AnalyzeWorkspace.h" 
#include "Framework/Framework/include/Utility.h" 
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "fastjet/ClusterSequence.hh"

#include <TH1D.h> 
#include <TH2D.h>
#include <TProfile.h>
#include <TProfile2D.h> 
#include <TStyle.h> 
#include <TCanvas.h> 
#include <TEfficiency.h> 
#include <TRandom3.h> 
#include <iostream> 
#include <TFile.h>

AnalyzeWorkspace::AnalyzeWorkspace() : inithistos(false)
{
}

//Define all your histograms here. 
void AnalyzeWorkspace::InitHistos(const std::map<std::string, bool>& cutMap, const std::vector<TH1DInfo>& histInfos) 
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    //This event counter histogram is necessary so that we know that all the condor jobs ran successfully. If not, when you use the hadder script, you will see a discrepancy in red as the files are being hadded.
    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ;

    //Define 1D histograms
    my_histos.emplace( "Cutflow", std::make_shared<TH1D>( "Cutflow", "Cutflow", 7, 0, 7 ) ) ;
    my_histos.emplace( "Events_by_Region", std::make_shared<TH1D>( "Events_by_Region", "Events_by_Region", 32, 0, 32 ) ) ;	

    for(auto& mycut : cutMap)
    {
        for(const auto& hInfo : histInfos)
        { 
            my_histos.emplace(hInfo.name+mycut.first, std::make_shared<TH1D>((hInfo.name+mycut.first).c_str(),(hInfo.name+mycut.first).c_str(), hInfo.nBins, hInfo.low, hInfo.high));
        }

    }    
}

//Put everything you want to do per event here. 
void AnalyzeWorkspace::Loop(NTupleReader& tr, double, int maxevents, bool) 
{
    while( tr.getNextEvent() )
    {
	const auto& eventCounter = tr.getVar<int>("eventCounter");

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

	bool onebJet = false;
	bool twobJets = false;
	bool threebJets = false;
	bool fourplusbJets = false;

	const auto& passMadHT = tr.getVar<bool>("passMadHT");
	
	if(maxevents != -1 && tr.getEvtNum() >= maxevents) break;
        if ( tr.getEvtNum() % 1000 == 0 ) printf(" Event %i\n", tr.getEvtNum() ) ;
	
        // ------------------------
        // -- Define weight
        // ------------------------
        double weight = 1.0;
        double eventweight = 1.0;
        double leptonScaleFactor = 1.0;
        //double prefiringScaleFactor = 1.0;
        //double puScaleFactor = 1.0;
        
        if(runtype == "MC")
        {
            if( !passMadHT ) continue; //Make sure not to double count DY events
            // Define Lumi weight
            const auto& Weight = tr.getVar<double>("Weight");
            const auto& lumi = tr.getVar<double>("Lumi");
	    std::cout << "lumi: " + std::to_string(lumi) << std::endl;
            eventweight = lumi*Weight;
            
            //PileupWeight = tr.getVar<double>("_PUweightFactor");
            //prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor");
            //puScaleFactor = tr.getVar<double>("puWeightCorr");
            
            weight *= eventweight*leptonScaleFactor;
            //weight *= eventweight*leptonScaleFactor*prefiringScaleFactor*puScaleFactor;
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
                if( Jetsgood >= 8 ) SRthree = true;
        }
	if( MJ > 1000 and allcut )
        {
                if( Jetsgood == 4 || Jetsgood == 5 ) CRthree = true;
                if( Jetsgood == 6 || Jetsgood == 7 ) SRfour = true;
                if( Jetsgood >= 8 ) SRfive = true;
        }

	if( bJetsgood == 1 ) onebJet = true;
	if( bJetsgood == 2 ) twobJets = true;
	if( bJetsgood == 3 ) threebJets = true;
	if( bJetsgood >= 4 ) fourplusbJets = true;

	bool placehold = true;

	const std::map<std::string, bool> cut_map
	{
		{"" , placehold },
		{"_OneGoodLepton" , OneGoodLepton },
		{"_OneGoodLepton_fourplusJets" , OneGoodLepton and fourplusJets },
		{"_OneGoodLepton_fourplusJets_oneplusbJets" , OneGoodLepton and fourplusJets and oneplusbJets },
		{"_OneGoodLepton_fourplusJets_oneplusbJets_passHT" , OneGoodLepton and fourplusJets and oneplusbJets and passHT },
		{"_OneGoodLepton_fourplusJets_oneplusbJets_passHT_passMJ" , OneGoodLepton and fourplusJets and oneplusbJets and passHT and passMJ },
		{"_OneGoodLepton_fourplusJets_oneplusbJets_passHT_passMJ_passTriggers" , passTriggers and OneGoodLepton and fourplusJets and oneplusbJets and passHT and passMJ },
		{"_CR1" , passTriggers and OneGoodLepton and fourplusJets and oneplusbJets and passHT and passMJ and CRone },
		{"_CR2" , passTriggers and OneGoodLepton and fourplusJets and oneplusbJets and passHT and passMJ and CRtwo },
		{"_CR3" , passTriggers and OneGoodLepton and fourplusJets and oneplusbJets and passHT and passMJ and CRthree },
		{"_SR1" , passTriggers and OneGoodLepton and fourplusJets and oneplusbJets and passHT and passMJ and SRone },
		{"_SR2" , passTriggers and OneGoodLepton and fourplusJets and oneplusbJets and passHT and passMJ and SRtwo },		
		{"_SR3" , passTriggers and OneGoodLepton and fourplusJets and oneplusbJets and passHT and passMJ and SRthree },
		{"_SR4" , passTriggers and OneGoodLepton and fourplusJets and oneplusbJets and passHT and passMJ and SRfour },
		{"_SR5" , passTriggers and OneGoodLepton and fourplusJets and oneplusbJets and passHT and passMJ and SRfive },
	};
	
	std::vector<TH1DInfo> histInfos = {
		{ "h_nGoodLeps", 10, 0, 10},
        	{ "h_nGoodJets", 30, 0, 30},
		{ "h_nGoodbJets", 15, 0, 15},
		{ "h_HT", 100, 0, 3000},
		{ "h_MJ", 50, 0, 2000},
		{ "h_largeRJets", 14, 0, 14},
        };
	
	if( inithistos == false )
	{
		InitHistos( cut_map, histInfos );
		inithistos = true;
	}
	
	for(auto& kv : cut_map)
	{
		if( kv.second )
		{
			my_histos["h_nGoodLeps" + kv.first]->Fill( lepgood , weight );
			my_histos["h_nGoodJets" + kv.first]->Fill( Jetsgood , weight );
			my_histos["h_nGoodbJets" + kv.first]->Fill( bJetsgood , weight );
			my_histos["h_HT" + kv.first]->Fill( eventHT , weight );
			my_histos["h_MJ" + kv.first]->Fill( MJ , weight );
			my_histos["h_largeRJets" + kv.first]->Fill( largeRjets.size() , weight );
			
		}
	}
	
        my_histos["EventCounter"]->Fill( eventCounter );

	my_histos["Cutflow"]->Fill( 0.5 , weight );
        if( OneGoodLepton ) my_histos["Cutflow"]->Fill( 1.5 , weight );
        if( OneGoodLepton and fourplusJets ) my_histos["Cutflow"]->Fill( 2.5 , weight );
        if( OneGoodLepton and fourplusJets and oneplusbJets ) my_histos["Cutflow"]->Fill( 3.5 , weight );
        if( OneGoodLepton and fourplusJets and oneplusbJets and passHT ) my_histos["Cutflow"]->Fill( 4.5 , weight );
        if( OneGoodLepton and fourplusJets and oneplusbJets and passHT and passMJ ) my_histos["Cutflow"]->Fill( 5.5 , weight );
        if( allcut ) my_histos["Cutflow"]->Fill( 6.5 , weight );

	if( allcut and CRone and onebJet ) my_histos["Events_by_Region"]->Fill( 0.5 , weight );
	if( allcut and CRone and twobJets ) my_histos["Events_by_Region"]->Fill( 1.5 , weight );
	if( allcut and CRone and threebJets ) my_histos["Events_by_Region"]->Fill( 2.5 , weight );
	if( allcut and CRone and fourplusbJets ) my_histos["Events_by_Region"]->Fill( 3.5 , weight );	
	if( allcut and CRtwo and onebJet ) my_histos["Events_by_Region"]->Fill( 4.5 , weight );
        if( allcut and CRtwo and twobJets ) my_histos["Events_by_Region"]->Fill( 5.5 , weight );
        if( allcut and CRtwo and threebJets ) my_histos["Events_by_Region"]->Fill( 6.5 , weight );
        if( allcut and CRtwo and fourplusbJets ) my_histos["Events_by_Region"]->Fill( 7.5 , weight );
	if( allcut and CRthree and onebJet ) my_histos["Events_by_Region"]->Fill( 8.5 , weight );
        if( allcut and CRthree and twobJets ) my_histos["Events_by_Region"]->Fill( 9.5 , weight );
        if( allcut and CRthree and threebJets ) my_histos["Events_by_Region"]->Fill( 10.5 , weight );
        if( allcut and CRthree and fourplusbJets ) my_histos["Events_by_Region"]->Fill( 11.5 , weight );
	if( allcut and SRone and onebJet ) my_histos["Events_by_Region"]->Fill( 12.5 , weight );
        if( allcut and SRone and twobJets ) my_histos["Events_by_Region"]->Fill( 13.5 , weight );
        if( allcut and SRone and threebJets ) my_histos["Events_by_Region"]->Fill( 14.5 , weight );
        if( allcut and SRone and fourplusbJets ) my_histos["Events_by_Region"]->Fill( 15.5 , weight );
	if( allcut and SRtwo and onebJet ) my_histos["Events_by_Region"]->Fill( 16.5 , weight );
        if( allcut and SRtwo and twobJets ) my_histos["Events_by_Region"]->Fill( 17.5 , weight );
        if( allcut and SRtwo and threebJets ) my_histos["Events_by_Region"]->Fill( 18.5 , weight );
        if( allcut and SRtwo and fourplusbJets ) my_histos["Events_by_Region"]->Fill( 19.5 , weight );
	if( allcut and SRthree and onebJet ) my_histos["Events_by_Region"]->Fill( 20.5 , weight );
        if( allcut and SRthree and twobJets ) my_histos["Events_by_Region"]->Fill( 21.5 , weight );
        if( allcut and SRthree and threebJets ) my_histos["Events_by_Region"]->Fill( 22.5 , weight );
        if( allcut and SRthree and fourplusbJets ) my_histos["Events_by_Region"]->Fill( 23.5 , weight );
	if( allcut and SRfour and onebJet ) my_histos["Events_by_Region"]->Fill( 24.5 , weight );
        if( allcut and SRfour and twobJets ) my_histos["Events_by_Region"]->Fill( 25.5 , weight );
        if( allcut and SRfour and threebJets ) my_histos["Events_by_Region"]->Fill( 26.5 , weight );
        if( allcut and SRfour and fourplusbJets ) my_histos["Events_by_Region"]->Fill( 27.5 , weight );
	if( allcut and SRfive and onebJet ) my_histos["Events_by_Region"]->Fill( 28.5 , weight );
        if( allcut and SRfive and twobJets ) my_histos["Events_by_Region"]->Fill( 29.5 , weight );
        if( allcut and SRfive and threebJets ) my_histos["Events_by_Region"]->Fill( 30.5 , weight );
        if( allcut and SRfive and fourplusbJets ) my_histos["Events_by_Region"]->Fill( 31.5 , weight );
	
    } 
}

void AnalyzeWorkspace::WriteHistos(TFile* outfile) 
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

