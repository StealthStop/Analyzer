#define CalculateBTagSF_cxx
#include "Analyzer/Analyzer/include/CalculateBTagSF.h"
#include "Framework/Framework/include/Utility.h"
#include "NTupleReader/include/NTupleReader.h"

#include <iostream>

// Calcualte b-tag eff. needed for the scale factor 
// Info from this Twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation
CalculateBTagSF::CalculateBTagSF() : initHistos(false)
{
}

void CalculateBTagSF::InitHistos(const std::string& histoFileTag)
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    const std::vector<double> ptBins  = { 20, 30, 45, 60, 75, 90, 120, 160, 210, 260, 320, 400, 500, 700, 1000 };
    const std::vector<double> etaBins = { -2.4, -1.4, 0.0, 1.4, 2.4 };
    const int nPtBins = ptBins.size() - 1;
    const int nEtaBins = etaBins.size() - 1;

    // Declare all your histograms here, that way we can fill them for multiple chains
    my_histos.emplace("EventCounter", std::make_shared<TH1D>("EventCounter","EventCounter", 2, -1.1, 1.1 ) );

    my_2d_histos.emplace( "n_eff_b_"+histoFileTag,     std::make_shared<TH2D>( ("n_eff_b_"+histoFileTag ).c_str(),     ( "n_eff_b_Efficiency_"+histoFileTag ).c_str(),     nPtBins, ptBins.data(), nEtaBins, etaBins.data() ) );
    my_2d_histos.emplace( "n_eff_c_"+histoFileTag,     std::make_shared<TH2D>( ("n_eff_c_"+histoFileTag ).c_str(),     ( "n_eff_c_Efficiency_"+histoFileTag ).c_str(),     nPtBins, ptBins.data(), nEtaBins, etaBins.data() ) );
    my_2d_histos.emplace( "n_eff_udsg_"+histoFileTag,  std::make_shared<TH2D>( ("n_eff_udsg_"+histoFileTag ).c_str(),  ( "n_eff_udsg_Efficiency_"+histoFileTag ).c_str(),  nPtBins, ptBins.data(), nEtaBins, etaBins.data() ) );
    my_2d_histos.emplace( "d_eff_b_"+histoFileTag,     std::make_shared<TH2D>( ("d_eff_b_"+histoFileTag ).c_str(),     ( "d_eff_b_Efficiency_"+histoFileTag ).c_str(),     nPtBins, ptBins.data(), nEtaBins, etaBins.data() ) );
    my_2d_histos.emplace( "d_eff_c_"+histoFileTag,     std::make_shared<TH2D>( ("d_eff_c_"+histoFileTag ).c_str(),     ( "d_eff_c_Efficiency_"+histoFileTag ).c_str(),     nPtBins, ptBins.data(), nEtaBins, etaBins.data() ) );
    my_2d_histos.emplace( "d_eff_udsg_"+histoFileTag,  std::make_shared<TH2D>( ("d_eff_udsg_"+histoFileTag ).c_str(),  ( "d_eff_udsg_Efficiency_"+histoFileTag ).c_str(),  nPtBins, ptBins.data(), nEtaBins, etaBins.data() ) );

    my_2d_histos["n_eff_b_"+histoFileTag]->GetXaxis()->SetTitle( "p_{T} [GeV]" );
    my_2d_histos["n_eff_b_"+histoFileTag]->GetYaxis()->SetTitle( "#eta" );
    my_2d_histos["n_eff_c_"+histoFileTag]->GetXaxis()->SetTitle( "p_{T} [GeV]" );
    my_2d_histos["n_eff_c_"+histoFileTag]->GetYaxis()->SetTitle( "#eta" );
    my_2d_histos["n_eff_udsg_"+histoFileTag]->GetXaxis()->SetTitle( "p_{T} [GeV]" );
    my_2d_histos["n_eff_udsg_"+histoFileTag]->GetYaxis()->SetTitle( "#eta" );
}//END of init histos

void CalculateBTagSF::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    while( tr.getNextEvent() )
    {
        const auto& filetag      = tr.getVar<std::string>("filetag");
        const auto& passMadHT    = tr.getVar<bool>("passMadHT");
        const auto& eventCounter = tr.getVar<int>("eventCounter");

        //-----------------------------------
        //-- Initialize Histograms
        //----------------------------------
        if( !initHistos ) 
        {
            InitHistos(filetag);
            initHistos = true;
        }

        //------------------------------------
        //-- Print Event Number
        //------------------------------------
        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() % 1000 == 0 ) printf( " Event %i\n", tr.getEvtNum() );

        //------------------------------------
        //-- Fill the histo
        //------------------------------------
        my_histos["EventCounter"]->Fill(eventCounter);
        
        if( passMadHT )
        {
            const auto& Jets           = tr.getVec<utility::LorentzVector>("Jets");
            const auto& GoodJets_pt30  = tr.getVec<bool>("GoodJets_pt30");
            const auto& recoJetsBtag   = tr.getVec<float>("Jets_bJetTagDeepCSVtotb");
            const auto& recoJetsFlavor = tr.getVec<int>("Jets_hadronFlavor");
            const auto& deepCSV_WP_medium = tr.getVar<double>("deepCSV_WP_medium");

            const auto& Weight = tr.getVar<float>("Weight");
            const auto& Lumi   = tr.getVar<double>("FinalLumi");
            const double eventweight = Lumi*Weight;
            
            for( unsigned int ij = 0; ij < Jets.size(); ++ij ) 
            {                
                if(!GoodJets_pt30[ij]) continue;
                const int myJetFlavor = std::abs( recoJetsFlavor.at( ij ) );
                const double myJetDeepCSV = recoJetsBtag.at( ij );
    
                if( myJetFlavor == 5 ) 
                {
                    my_2d_histos["d_eff_b_"+filetag]->Fill( Jets.at( ij ).Pt(), Jets.at( ij ).Eta(), eventweight );
                    if( myJetDeepCSV > deepCSV_WP_medium ) my_2d_histos["n_eff_b_"+filetag]->Fill( Jets.at( ij ).Pt(), Jets.at( ij ).Eta(), eventweight ); 
                }
                else if( myJetFlavor == 4 ) 
                {
                    my_2d_histos["d_eff_c_"+filetag]->Fill( Jets.at( ij ).Pt(), Jets.at( ij ).Eta(), eventweight );
                    if( myJetDeepCSV > deepCSV_WP_medium ) my_2d_histos["n_eff_c_"+filetag]->Fill( Jets.at( ij ).Pt(), Jets.at( ij ).Eta(), eventweight ); 
                }
                else if( myJetFlavor < 4 || myJetFlavor == 21 ) 
                {
                    my_2d_histos["d_eff_udsg_"+filetag]->Fill( Jets.at( ij ).Pt(), Jets.at( ij ).Eta(), eventweight );
                    if( myJetDeepCSV > deepCSV_WP_medium ) my_2d_histos["n_eff_udsg_"+filetag]->Fill( Jets.at( ij ).Pt(), Jets.at( ij ).Eta(), eventweight ); 
                }
            }
        }
    }//END of while tr.getNextEvent loop   
}//END of function
      
void CalculateBTagSF::WriteHistos( TFile* outfile ) 
{
    outfile->cd();

    for( const auto& p : my_histos ) 
    {
        p.second->SetDirectory(outfile);
        p.second->Write();
    }
    
    for( const auto& p : my_2d_histos ) 
    {
        p.second->SetDirectory(outfile);
        p.second->Write();
    }
}
