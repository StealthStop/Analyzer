#define CalculateTopTagSF_cxx
#include "Analyzer/Analyzer/include/CalculateTopTagSF.h"
#include "Framework/Framework/include/Utility.h"
#include "NTupleReader/include/NTupleReader.h"

#include "TopTagger/TopTagger/interface/TopTaggerResults.h"
#include "TopTagger/TopTagger/interface/TopObject.h"

#include <iostream>

// Calcualte top-tag eff. needed for the final scale factor a la b tagging
// Info from this Twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation (method 1a)
CalculateTopTagSF::CalculateTopTagSF() : initHistos(false)
{
}

void CalculateTopTagSF::InitHistos(const std::string& histoFileTag)
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    // Pt bin edges for resolved case are directly from AN-18-273 v2, Fig. 7
    const std::vector<double> ptBinsRes  = { 0, 150, 250, 300, 350, 400, 450, 500, 600, 1000 };
    const std::vector<double> ptBinsMrg  = { 400, 450, 500, 600, 800, 1000 };
    const std::vector<double> etaBins    = { -2.4, -1.4, 0.0, 1.4, 2.4 };
    const int nPtBinsRes = ptBinsRes.size() - 1;
    const int nPtBinsMrg = ptBinsMrg.size() - 1;
    const int nEtaBins = etaBins.size() - 1;

    // Declare all your histograms here, that way we can fill them for multiple chains
    my_histos.emplace("EventCounter", std::make_shared<TH1D>("EventCounter","EventCounter", 2, -1.1, 1.1 ) );

    my_2d_histos.emplace( "n_eff_res_"+histoFileTag, std::make_shared<TH2D>( ( "n_eff_res_"+histoFileTag ).c_str(), ( "n_eff_res_Efficiency_"+histoFileTag ).c_str(), nPtBinsRes, ptBinsRes.data(), nEtaBins, etaBins.data() ) );
    my_2d_histos.emplace( "n_eff_mrg_"+histoFileTag, std::make_shared<TH2D>( ( "n_eff_mrg_"+histoFileTag ).c_str(), ( "n_eff_mrg_Efficiency_"+histoFileTag ).c_str(), nPtBinsMrg, ptBinsMrg.data(), nEtaBins, etaBins.data() ) );
    my_2d_histos.emplace( "d_eff_res_"+histoFileTag, std::make_shared<TH2D>( ( "d_eff_res_"+histoFileTag ).c_str(), ( "d_eff_res_Efficiency_"+histoFileTag ).c_str(), nPtBinsRes, ptBinsRes.data(), nEtaBins, etaBins.data() ) );
    my_2d_histos.emplace( "d_eff_mrg_"+histoFileTag, std::make_shared<TH2D>( ( "d_eff_mrg_"+histoFileTag ).c_str(), ( "d_eff_mrg_Efficiency_"+histoFileTag ).c_str(), nPtBinsMrg, ptBinsMrg.data(), nEtaBins, etaBins.data() ) );

    my_2d_histos.emplace( "n_mis_res_"+histoFileTag, std::make_shared<TH2D>( ( "n_mis_res_"+histoFileTag ).c_str(), ( "n_mis_res_Mistag_"+histoFileTag ).c_str(),     nPtBinsRes, ptBinsRes.data(), nEtaBins, etaBins.data() ) );
    my_2d_histos.emplace( "n_mis_mrg_"+histoFileTag, std::make_shared<TH2D>( ( "n_mis_mrg_"+histoFileTag ).c_str(), ( "n_mis_mrg_Mistag_"+histoFileTag ).c_str(),     nPtBinsMrg, ptBinsMrg.data(), nEtaBins, etaBins.data() ) );
    my_2d_histos.emplace( "d_mis_res_"+histoFileTag, std::make_shared<TH2D>( ( "d_mis_res_"+histoFileTag ).c_str(), ( "d_mis_res_Mistag_"+histoFileTag ).c_str(),     nPtBinsRes, ptBinsRes.data(), nEtaBins, etaBins.data() ) );
    my_2d_histos.emplace( "d_mis_mrg_"+histoFileTag, std::make_shared<TH2D>( ( "d_mis_mrg_"+histoFileTag ).c_str(), ( "d_mis_mrg_Mistag_"+histoFileTag ).c_str(),     nPtBinsMrg, ptBinsMrg.data(), nEtaBins, etaBins.data() ) );

    my_2d_histos["n_eff_res_"+histoFileTag]->GetXaxis()->SetTitle( "p_{T} [GeV]" );
    my_2d_histos["n_eff_res_"+histoFileTag]->GetYaxis()->SetTitle( "#eta" );
    my_2d_histos["n_eff_mrg_"+histoFileTag]->GetXaxis()->SetTitle( "p_{T} [GeV]" );
    my_2d_histos["n_eff_mrg_"+histoFileTag]->GetYaxis()->SetTitle( "#eta" );
    my_2d_histos["n_mis_res_"+histoFileTag]->GetXaxis()->SetTitle( "p_{T} [GeV]" );
    my_2d_histos["n_mis_res_"+histoFileTag]->GetYaxis()->SetTitle( "#eta" );
    my_2d_histos["n_mis_mrg_"+histoFileTag]->GetXaxis()->SetTitle( "p_{T} [GeV]" );
    my_2d_histos["n_mis_mrg_"+histoFileTag]->GetYaxis()->SetTitle( "#eta" );
}

void CalculateTopTagSF::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    while( tr.getNextEvent() )
    {
        const auto& filetag      = tr.getVar<std::string>("filetag");
        const auto& runyear      = tr.getVar<std::string>("runYear");
        const auto& passMadHT    = tr.getVar<bool>("passMadHT");
        const auto& eventCounter = tr.getVar<int>("eventCounter");

        // Distinguish if a top candidate is resolved or merged
        double resolvedWP = 0.95;
        double mergedWP   = 0.937;
        if ( runyear.find("2016") == std::string::npos )
            mergedWP = 0.895;

        //-----------------------------------
        //-- Initialize Histograms
        //-----------------------------------
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
            const auto* topTagRes = tr.getVar<TopTaggerResults*>("ttr");

            const auto& Weight = tr.getVar<float>("Weight");
            const auto& Lumi   = tr.getVar<double>("FinalLumi");
            const double eventweight = Lumi * Weight;

            const auto& tops = topTagRes->getTops();
            for( const auto& top : tops ) 
            {                
                const auto* genTop = top->getBestGenTopMatch();

                // Ensure that the last pt bin is inclusive
                double thePt = top->P().Pt();
                if ( thePt > 1000.0 )
                    thePt = 999.0;

                // Ensure the edge eta bins are inclusive
                // Not necessary due to internal top tagger requirement of a top
                // being within |eta| < 2.0, but just put here for symmetry
                double theEta = top->P().Eta();
                if      ( theEta < -2.4 )
                    theEta = -2.3;
                else if ( theEta > 2.4 )
                    theEta = 2.3;
 
                // Efficiency when dealing with an actual top
                if ( genTop )
                {
                    if      ( top->getType() == TopObject::RESOLVED_TOP )
                    {
                        my_2d_histos["d_eff_res_"+filetag]->Fill( thePt, theEta, eventweight );
                        if( top->getDiscriminator() > resolvedWP )
                            my_2d_histos["n_eff_res_"+filetag]->Fill( thePt, theEta, eventweight ); 
                    }
                    else if ( top->getType() == TopObject::MERGED_TOP )
                    {
                        my_2d_histos["d_eff_mrg_"+filetag]->Fill( thePt, theEta, eventweight );
                        if( top->getDiscriminator() > mergedWP )
                            my_2d_histos["n_eff_mrg_"+filetag]->Fill( thePt, theEta, eventweight ); 
                    }
                }
                // Mistag when dealing with a fake top i.e. no GEN top present
                else
                {
                    if      ( top->getType() == TopObject::RESOLVED_TOP )
                    {
                        my_2d_histos["d_mis_res_"+filetag]->Fill( thePt, theEta, eventweight );
                        if( top->getDiscriminator() > resolvedWP )
                            my_2d_histos["n_mis_res_"+filetag]->Fill( thePt, theEta, eventweight ); 
                    }
                    else if ( top->getType() == TopObject::MERGED_TOP )
                    {
                        my_2d_histos["d_mis_mrg_"+filetag]->Fill( thePt, theEta, eventweight );
                        if( top->getDiscriminator() > mergedWP )
                            my_2d_histos["n_mis_mrg_"+filetag]->Fill( thePt, theEta, eventweight ); 
                    }
                }
            }
        }
    }
}
      
void CalculateTopTagSF::WriteHistos( TFile* outfile ) 
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
