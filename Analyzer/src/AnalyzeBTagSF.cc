#define AnalyzeBTagSF_cxx
#include "Analyzer/Analyzer/include/AnalyzeBTagSF.h"
#include "Framework/Framework/include/Utility.h"
#include "SusyAnaTools/Tools/NTupleReader.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TRandom3.h>
#include <iostream>

AnalyzeBTagSF::AnalyzeBTagSF()
{
    InitHistos();
}


void AnalyzeBTagSF::InitHistos()
{
    
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    // Declare all your histograms here, that way we can fill them for multiple chains

    my_histos.emplace( "h_btagSF",          std::make_shared<TH1D>( "h_btagSF",         "h_btagSF",         60, 0.0, 2.0 ) );
    my_histos.emplace( "h_btagSF_u",        std::make_shared<TH1D>( "h_btagSF_u",       "h_btagSF_u",       60, 0.0, 2.0 ) );
    my_histos.emplace( "h_btagSF_d",        std::make_shared<TH1D>( "h_btagSF_d",       "h_btagSF_d",       60, 0.0, 2.0 ) );

    my_histos.emplace( "h_nb",              std::make_shared<TH1D>( "h_nb",             "h_nb",             10, 0.0, 10.0 ) );
    my_histos.emplace( "h_nb_btagSF",       std::make_shared<TH1D>( "h_nb_btagSF",      "h_nb_btagSF",      10, 0.0, 10.0 ) );
    my_histos.emplace( "h_nb_btagSF_u",     std::make_shared<TH1D>( "h_nb_btagSF_u",    "h_nb_btagSF_u",    10, 0.0, 10.0 ) );
    my_histos.emplace( "h_nb_btagSF_d",     std::make_shared<TH1D>( "h_nb_btagSF_d",    "h_nb_btagSF_d",    10, 0.0, 10.0 ) );
    
    my_histos.emplace( "h_CSV",             std::make_shared<TH1D>( "h_CSV",            "h_CSV",            60, 0.0, 1.0 ) );
    my_histos.emplace( "h_CSV_btagSF",      std::make_shared<TH1D>( "h_CSV_btagSF",     "h_CSV_btagSF",     60, 0.0, 1.0 ) );
    my_histos.emplace( "h_CSV_btagSF_u",    std::make_shared<TH1D>( "h_CSV_btagSF_u",   "h_CSV_btagSF_u",   60, 0.0, 1.0 ) );
    my_histos.emplace( "h_CSV_btagSF_d",    std::make_shared<TH1D>( "h_CSV_btagSF_d",   "h_CSV_btagSF_d",   60, 0.0, 1.0 ) );
    
    my_histos.emplace( "h_mistagSF_u",      std::make_shared<TH1D>( "h_mistagSF_u",     "h_mistagSF_u",     60, 0.0, 2.0 ) );
    my_histos.emplace( "h_mistagSF_d",      std::make_shared<TH1D>( "h_mistagSF_d",     "h_mistagSF_d",     60, 0.0, 2.0 ) );
    
    my_histos.emplace( "h_puSF",            std::make_shared<TH1D>( "h_puSF",           "h_puSF",           60, 0.0, 2.0 ) );
    //my_histos.emplace( "h_puSF_JP",         std::make_shared<TH1D>( "h_puSF_JP",        "h_puSF_JP",        60, 0.0, 2.0 ) );

    my_histos.emplace( "h_pdfSF",           std::make_shared<TH1D>( "h_pdfSF",          "h_pdfSF",          60, 0.0, 2.0 ) );
    my_histos.emplace( "h_pdfSF_u",         std::make_shared<TH1D>( "h_pdfSF_u",        "h_pdfSF_u",        60, 0.0, 2.0 ) );
    my_histos.emplace( "h_pdfSF_d",         std::make_shared<TH1D>( "h_pdfSF_d",        "h_pdfSF_d",        60, 0.0, 2.0 ) );

    my_histos.emplace( "h_scaleSF",         std::make_shared<TH1D>( "h_scaleSF",        "h_scaleSF",        60, 0.0, 2.0 ) );
    my_histos.emplace( "h_scaleSF_u",       std::make_shared<TH1D>( "h_scaleSF_u",      "h_scaleSF_u",      60, 0.0, 2.0 ) );
    my_histos.emplace( "h_scaleSF_d",       std::make_shared<TH1D>( "h_scaleSF_d",      "h_scaleSF_d",      60, 0.0, 2.0 ) );

    my_histos.emplace( "h_nVtx_nPUSF",      std::make_shared<TH1D>( "h_nVtx_nPUSF",     "h_nVtx_nPUSF",     50, 0.0, 50.0) ); 
    my_histos.emplace( "h_nVtx_wPUSF",      std::make_shared<TH1D>( "h_nVtx_wPUSF",     "h_nVtx_wPUSF",     50, 0.0, 50.0) );
    //my_histos.emplace( "h_nVtx_wPUSF_JP",   std::make_shared<TH1D>( "h_nVtx_wPUSF_JP",  "h_nVtx_wPUSF_JP",  50, 0.0, 50.0) );

    my_histos.emplace( "h_nj_scaleSF",      std::make_shared<TH1D>( "h_nj_scaleSF",     "h_nj_scaleSF",     20, 0.0, 20.0) );
    my_histos.emplace( "h_nj_scaleSF_u",    std::make_shared<TH1D>( "h_nj_scaleSF_u",   "h_nj_scaleSF_u",   20, 0.0, 20.0) );
    my_histos.emplace( "h_nj_scaleSF_d",    std::make_shared<TH1D>( "h_nj_scaleSF_d",   "h_nj_scaleSF_d",   20, 0.0, 20.0) );
    
    my_histos.emplace( "h_nj_pdfSF",        std::make_shared<TH1D>( "h_nj_pdfSF",       "h_nj_pdfSF",       20, 0.0, 20.0) );
    my_histos.emplace( "h_nj_pdfSF_u",      std::make_shared<TH1D>( "h_nj_pdfSF_u",     "h_nj_pdfSF_u",     20, 0.0, 20.0) );
    my_histos.emplace( "h_nj_pdfSF_d",      std::make_shared<TH1D>( "h_nj_pdfSF_d",     "h_nj_pdfSF_d",     20, 0.0, 20.0) );

    my_histos.emplace( "h_eleSF",           std::make_shared<TH1D>( "h_eleSF",          "h_eleSF",          60, -1.2, 2.0 ) );
    my_histos.emplace( "h_muSF",            std::make_shared<TH1D>( "h_muSF",           "h_muSF",           60, -1.2, 2.0 ) );
/*    std::vector<std::string> genericVectorOfTags    { "myTag" };

    for( std::string nJetsTag : njets_fisherDist_tags ) {
        my_histos.emplace( "h_fishDist_"+nJetsTag+"j_"+jetPtTag+"_"+fishBinTag, std::make_shared<TH1D>( ( "h_fishDist_"+nJetsTag+"j_"+jetPtTag+"_"+fishBinTag ).c_str(), ( "h_fishDist_"+nJetsTag+"j_"+jetPtTag+"_"+fishBinTag ).c_str(), 120, -0.6, 0.6 ) );
    }//END of njets_fisherDist_tags*/

}//END of init histos

void AnalyzeBTagSF::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    while( tr.getNextEvent() )
    {
        const auto& NJets_pt30          = tr.getVar<int>("NGoodJets_pt30");
        const auto& NBJets_pt30         = tr.getVar<int>("NGoodBJets_pt30");
        const auto& Jets_CSV            = tr.getVec<float>("Jets_bDiscriminatorCSV");

        const auto& passMadHT           = tr.getVar<bool>("passMadHT");
        
        const auto& bTagSF              = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
        const auto& bTagSF_u            = tr.getVar<double>("bTagSF_EventWeightSimple_Up");
        const auto& bTagSF_d            = tr.getVar<double>("bTagSF_EventWeightSimple_Down");
        
        const auto& bMisTagSF_u         = tr.getVar<float>("mistagSF_EventWeightSimple_Up");
        const auto& bMisTagSF_d         = tr.getVar<float>("mistagSF_EventWeightSimple_Down");
        
        const auto& puSF                = tr.getVar<float>("puWeight");
        //const auto& puSF_JP             = tr.getVar<float>("_PUweightFactor");
        const auto& ntru_PV             = tr.getVar<float>("TrueNumInteractions");

        const auto& scaleSF             = tr.getVar<float>("scaleWeightNom");
        const auto& scaleSF_u           = tr.getVar<float>("scaleWeightUp");
        const auto& scaleSF_d           = tr.getVar<float>("scaleWeightDown");

        const auto& pdfSF               = tr.getVar<float>("PDFweightNom");
        const auto& pdfSF_u             = tr.getVar<float>("PDFweightUp");
        const auto& pdfSF_d             = tr.getVar<float>("PDFweightDown");

        const auto& electronSF          = tr.getVar<float>("leadGoodElectronSF");
        const auto& muonSF              = tr.getVar<float>("leadGoodMuonSF");


        //------------------------------------
        //-- Print Event Number
        //------------------------------------

        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() % 1000 == 0 ) printf( " Event %i\n", tr.getEvtNum() );

        //------------------------------------
        //-- Get the Proper Event Weight
        //------------------------------------

        double eventweight          = 1.0; //For data
        double PUaddweight          = 1.0;//, PU_JPaddweight           = 1.0;
        double scaleSFaddweight     = 1.0, scaleSFaddweight_u       = 1.0, scaleSFaddweight_d       = 1.0;
        double pdfSFaddweight       = 1.0, pdfSFaddweight_u         = 1.0, pdfSFaddweight_d         = 1.0;
        double bTagSFaddweight      = 1.0, bTagSFaddweight_u        = 1.0, bTagSFaddweight_d        = 1.0;
        double Lumi = 35900;

        if( passMadHT )
        {
            const auto& Weight = tr.getVar<float>("Weight");
            
            eventweight         = Lumi*Weight;
            PUaddweight         = Lumi*Weight*puSF;
            //PU_JPaddweight      = Lumi*Weight*puSF_JP;

            scaleSFaddweight    = Lumi*Weight*scaleSF;
            scaleSFaddweight_u  = Lumi*Weight*scaleSF_u;
            scaleSFaddweight_d  = Lumi*Weight*scaleSF_d;

            pdfSFaddweight      = Lumi*Weight*pdfSF;
            pdfSFaddweight_u    = Lumi*Weight*pdfSF_u;
            pdfSFaddweight_d    = Lumi*Weight*pdfSF_d;

            bTagSFaddweight     = Lumi*Weight*bTagSF;
            bTagSFaddweight_u   = Lumi*Weight*bTagSF_u;
            bTagSFaddweight_d   = Lumi*Weight*bTagSF_d;

        }

        else continue;
 
        my_histos["h_muSF"]->Fill(muonSF);
        my_histos["h_eleSF"]->Fill(electronSF);
        //std::cout<<ntru_PV<<std::endl;
        my_histos["h_btagSF"]->Fill(bTagSF);
        my_histos["h_btagSF_u"]->Fill(bTagSF_u);
        my_histos["h_btagSF_d"]->Fill(bTagSF_d);
        
        my_histos["h_mistagSF_u"]->Fill(bMisTagSF_u);
        my_histos["h_mistagSF_d"]->Fill(bMisTagSF_d);

        my_histos["h_nb"]->Fill(NBJets_pt30, eventweight);
        my_histos["h_nb_btagSF"]->Fill(NBJets_pt30, bTagSFaddweight);
        my_histos["h_nb_btagSF_u"]->Fill(NBJets_pt30, bTagSFaddweight_u);
        my_histos["h_nb_btagSF_d"]->Fill(NBJets_pt30, bTagSFaddweight_d);
       
        if(Jets_CSV.size() != 0) {
            my_histos["h_CSV"]->Fill(Jets_CSV.at(0), eventweight);
            my_histos["h_CSV_btagSF"]->Fill(Jets_CSV.at(0), bTagSFaddweight);
            my_histos["h_CSV_btagSF_u"]->Fill(Jets_CSV.at(0), bTagSFaddweight_u);
            my_histos["h_CSV_btagSF_d"]->Fill(Jets_CSV.at(0), bTagSFaddweight_d);
        }
        my_histos["h_puSF"]->Fill(puSF);
        //my_histos["h_puSF_JP"]->Fill(puSF_JP);

        my_histos["h_nVtx_nPUSF"]->Fill(ntru_PV,eventweight);
        my_histos["h_nVtx_wPUSF"]->Fill(ntru_PV,PUaddweight);
        //my_histos["h_nVtx_wPUSF_JP"]->Fill(ntru_PV,PU_JPaddweight);

        my_histos["h_scaleSF"]->Fill(scaleSF);
        my_histos["h_scaleSF_u"]->Fill(scaleSF_u);
        my_histos["h_scaleSF_d"]->Fill(scaleSF_d);

        if( !std::isfinite(NJets_pt30) ) std::cout<<"NJets_pt30 is not finite"<<std::endl;
        if( !std::isfinite(scaleSFaddweight) ) std::cout<<"Scale SF add weight is not finite"<<std::endl;
        if( !std::isfinite(pdfSFaddweight) ) std::cout<<"PDF SF add weight is not finite"<<std::endl;
        my_histos["h_nj_scaleSF"]->Fill(NJets_pt30, scaleSFaddweight);
        my_histos["h_nj_scaleSF_u"]->Fill(NJets_pt30, scaleSFaddweight_u);
        my_histos["h_nj_scaleSF_d"]->Fill(NJets_pt30, scaleSFaddweight_d);

        my_histos["h_pdfSF"]->Fill(pdfSF);
        my_histos["h_pdfSF_u"]->Fill(pdfSF_u);
        my_histos["h_pdfSF_d"]->Fill(pdfSF_d);

        my_histos["h_nj_pdfSF"]->Fill(NJets_pt30, pdfSFaddweight);
        my_histos["h_nj_pdfSF_u"]->Fill(NJets_pt30, pdfSFaddweight_u);
        my_histos["h_nj_pdfSF_d"]->Fill(NJets_pt30, pdfSFaddweight_d);

    }//END of while tr.getNextEvent loop   
}//END of function
      
void AnalyzeBTagSF::WriteHistos( TFile* outfile ) 
{
    outfile->cd();

    for (const auto &p : my_histos) {
        p.second->SetDirectory(outfile);
        p.second->Write();
    }
    
    for (const auto &p : my_2d_histos) {
        p.second->SetDirectory(outfile);
        p.second->Write();
    }
    
    for (const auto &p : my_efficiencies) {
        p.second->SetDirectory(outfile);
        p.second->Write();
    }   
}

