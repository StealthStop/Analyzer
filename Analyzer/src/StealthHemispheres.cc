#include "Analyzer/Analyzer/include/StealthHemispheres.h"
#include "SusyAnaTools/Tools/NTupleReader.h"

#include "TopTagger/TopTagger/interface/TopTagger.h"
#include "TopTagger/TopTagger/interface/TopTaggerResults.h"
#include "TopTagger/TopTagger/interface/TopTaggerUtilities.h"
#include "TopTagger/CfgParser/interface/TTException.h"
#include "Framework/Framework/include/SetUpTopTagger.h"


#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include "TLorentzVector.h"
#include "TMath.h"
#include "TRandom.h"

StealthHemispheres::StealthHemispheres() : inithisto(false)
{
}

void StealthHemispheres::InitHistos(const std::map<std::string, bool>& cutmap) 
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ;

    std::vector<std::string> sufs = {"_TaggedTop", "_0l", "_0l_cm"};
    for (const auto& suf    : sufs  ) 
    {
    for (const auto& cutVar : cutmap)
    {
        // 1 Lepton Case
        my_histos.emplace( "h_MT2"+suf + "_" + cutVar.first, std::make_shared<TH1D> ( ("h_MT2" + suf + "_" + cutVar.first).c_str(), ("h_MT2" + suf + "_" + cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_stop1Mass"+suf + "_" + cutVar.first, std::make_shared<TH1D> ( ("h_stop1Mass" + suf + "_" + cutVar.first).c_str(), ("h_stop1Mass" + suf + "_" + cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_stop1Eta"+suf + "_" + cutVar.first, std::make_shared<TH1D> ( ("h_stop1Eta" + suf + "_" + cutVar.first).c_str(), ("h_stop1Eta" + suf + "_" + cutVar.first).c_str(), 100, -6, 6 ) );
        my_histos.emplace( "h_stop1Phi"+suf + "_" + cutVar.first, std::make_shared<TH1D> ( ("h_stop1Phi" + suf + "_" + cutVar.first).c_str(), ("h_stop1Phi" + suf + "_" + cutVar.first).c_str(), 80, -4, 4 ) );
        my_histos.emplace( "h_stop1Pt"+suf + "_" + cutVar.first, std::make_shared<TH1D> ( ("h_stop1Pt" + suf + "_" + cutVar.first).c_str(), ("h_stop1Pt" + suf + "_" + cutVar.first).c_str(), 100, 0, 1000 ) );
        my_histos.emplace( "h_stop2Mass"+suf + "_" + cutVar.first, std::make_shared<TH1D> ( ("h_stop2Mass" + suf + "_" + cutVar.first).c_str(), ("h_stop2Mass" + suf + "_" + cutVar.first).c_str(), 500, 0, 1500) );
        my_histos.emplace( "h_stop2Eta"+suf + "_" + cutVar.first, std::make_shared<TH1D> ( ("h_stop2Eta" + suf + "_" + cutVar.first).c_str(), ("h_stop2Eta" + suf + "_" + cutVar.first).c_str(), 100, -6, 6 ) );
        my_histos.emplace( "h_stop2Phi"+suf + "_" + cutVar.first, std::make_shared<TH1D> ( ("h_stop2Phi" + suf + "_" + cutVar.first).c_str(), ("h_stop2Phi" + suf + "_" + cutVar.first).c_str(), 80, -4, 4 ) );
        my_histos.emplace( "h_stop2Pt"+suf + "_" + cutVar.first, std::make_shared<TH1D> ( ("h_stop2Pt" + suf + "_" + cutVar.first).c_str(), ("h_stop2Pt" + suf + "_" + cutVar.first).c_str(), 100, 0, 1000 ) );
        my_histos.emplace( "h_seed1Pt"+suf + "_" + cutVar.first, std::make_shared<TH1D> ( ("h_seed1Pt" + suf + "_" + cutVar.first).c_str(), ("h_seed1Pt" + suf + "_" + cutVar.first).c_str(), 100, 0, 1000 ) );
        my_histos.emplace( "h_seed1Eta"+suf + "_" + cutVar.first, std::make_shared<TH1D> ( ("h_seed1Eta" + suf + "_" + cutVar.first).c_str(), ("h_seed1Eta" + suf + "_" + cutVar.first).c_str(), 100, -6, 6 ) );
        my_histos.emplace( "h_seed1Phi"+suf + "_" + cutVar.first, std::make_shared<TH1D> ( ("h_seed1Phi" + suf + "_" + cutVar.first).c_str(), ("h_seed1Phi" + suf + "_" + cutVar.first).c_str(), 80, -4, 4 ) );
        my_histos.emplace( "h_seed1Mass"+suf + "_" + cutVar.first, std::make_shared<TH1D> ( ("h_seed1Mass" + suf + "_" + cutVar.first).c_str(), ("h_seed1Mass" + suf + "_" + cutVar.first).c_str(), 100, 0, 1000 ) );
        my_histos.emplace( "h_seed2Pt"+suf + "_" + cutVar.first, std::make_shared<TH1D> ( ("h_seed2Pt" + suf + "_" + cutVar.first).c_str(), ("h_seed2Pt" + suf + "_" + cutVar.first).c_str(), 100, 0, 1000 ) );
        my_histos.emplace( "h_seed2Eta"+suf + "_" + cutVar.first, std::make_shared<TH1D> ( ("h_seed2Eta" + suf + "_" + cutVar.first).c_str(), ("h_seed2Eta" + suf + "_" + cutVar.first).c_str(), 100, -6, 6 ) );
        my_histos.emplace( "h_seed2Phi"+suf + "_" + cutVar.first, std::make_shared<TH1D> ( ("h_seed2Phi" + suf + "_" + cutVar.first).c_str(), ("h_seed2Phi" + suf + "_" + cutVar.first).c_str(), 100, -4, 4 ) );
        my_histos.emplace( "h_seed2Mass"+suf + "_" + cutVar.first, std::make_shared<TH1D> ( ("h_seed2Mass" + suf + "_" + cutVar.first).c_str(), ("h_seed2Mass" + suf + "_" + cutVar.first).c_str(), 100, 0, 1000 ) );
        my_histos.emplace( "h_dR_seed1seed2"+suf + "_" + cutVar.first, std::make_shared<TH1D> ( ("h_dR_seed1seed2" + suf + "_" + cutVar.first).c_str(), ("h_dR_seed1seed2" + suf + "_" + cutVar.first).c_str(), 50, 0, 10 ) );
        my_histos.emplace( "h_dR_stop1stop2"+suf + "_" + cutVar.first, std::make_shared<TH1D> ( ("h_dR_stop1stop2" + suf + "_" + cutVar.first).c_str(), ("h_dR_stop1stop2" + suf + "_" + cutVar.first).c_str(), 50, 0, 10 ) );
        my_histos.emplace( "h_dPhi_stop1stop2"+suf + "_" + cutVar.first, std::make_shared<TH1D> ( ("h_dPhi_stop1stop2" + suf + "_" + cutVar.first).c_str(), ("h_dPhi_stop1stop2" + suf + "_" + cutVar.first).c_str(), 50, 0, 10 ) );
        my_2d_histos.emplace( "h_Mass_stop1vsstop2"+suf + "_" + cutVar.first, std::make_shared<TH2D>( ("h_Mass_stop1vsstop2" + suf + "_" + cutVar.first).c_str(), ("h_Mass_stop1vsstop2" + suf + "_" + cutVar.first).c_str(), 500, 0, 1500, 500, 0, 1500 ) );
        my_2d_histos.emplace( "h_Eta_stop1vsstop2"+suf + "_" + cutVar.first, std::make_shared<TH2D>( ("h_Eta_stop1vsstop2" + suf + "_" + cutVar.first).c_str(), ("h_Eta_stop1vsstop2" + suf + "_" + cutVar.first).c_str(), 100, -6, 6, 100, -6, 6 ) );
        my_2d_histos.emplace( "h_Phi_stop1vsstop2"+suf + "_" + cutVar.first, std::make_shared<TH2D>( ("h_Phi_stop1vsstop2" + suf + "_" + cutVar.first).c_str(), ("h_Phi_stop1vsstop2" + suf + "_" + cutVar.first).c_str(), 80, -4, 4, 80, -4, 4 ) );
        my_2d_histos.emplace( "h_Pt_stop1vsstop2"+suf + "_" + cutVar.first, std::make_shared<TH2D>( ("h_Pt_stop1vsstop2" + suf + "_" + cutVar.first).c_str(), ("h_Pt_stop1vsstop2" + suf + "_" + cutVar.first).c_str(), 100, 0, 1000, 100, 0, 1000 ) );
        my_2d_histos.emplace( "h_MT2VsNJets"+suf + "_" + cutVar.first, std::make_shared<TH2D> ( ("h_MT2VsNJets" + suf + "_" + cutVar.first).c_str(), ("h_MT2VsNJets" + suf + "_" + cutVar.first).c_str(), 500, 0, 1500, 20, 0, 20 ) );
        my_2d_histos.emplace( "h_MT2Vsstop1Mass"+suf + "_" + cutVar.first, std::make_shared<TH2D> ( ("h_MT2Vsstop1Mass" + suf + "_" + cutVar.first).c_str(), ("h_MT2Vsstop1Mass" + suf + "_" + cutVar.first).c_str(), 500, 0, 1500, 500, 0, 1500 ) );
        my_2d_histos.emplace( "h_MT2Vsstop2Mass"+suf + "_" + cutVar.first, std::make_shared<TH2D> ( ("h_MT2Vsstop2Mass" + suf + "_" + cutVar.first).c_str(), ("h_MT2Vsstop2Mass" + suf + "_" + cutVar.first).c_str(), 500, 0, 1500, 500, 0, 1500 ) );
        my_2d_histos.emplace( "h_stop1MassVsNJets"+suf + "_" + cutVar.first, std::make_shared<TH2D> ( ("h_stop1MassVsNJets" + suf + "_" + cutVar.first).c_str(), ("h_stop1MassVsNJets" + suf + "_" + cutVar.first).c_str(), 500, 0, 1500, 20, 0, 20 ) );
        my_2d_histos.emplace( "h_stop2MassVsNJets"+suf + "_" + cutVar.first, std::make_shared<TH2D> ( ("h_stop2MassVsNJets" + suf + "_" + cutVar.first).c_str(), ("h_stop2MassVsNJets" + suf + "_" + cutVar.first).c_str(), 500, 0, 1500, 20, 0, 20 ) );

    }
    }
}

void StealthHemispheres::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    while( tr.getNextEvent() )
    {

        const auto& eventCounter    = tr.getVar<int>("eventCounter");
        
        // ------------------------
        // -- Print event number   
        // ------------------------     
        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() & (10000 == 0) ) printf( " Event %i\n", tr.getEvtNum() );

        const auto& runtype             = tr.getVar<std::string>("runtype");
        const auto& filetag             = tr.getVar<std::string>("filetag");
        const auto& NGoodJets_pt45      = tr.getVar<int>("NGoodJets_pt45");
        const auto& passBaseline0l_Good = tr.getVar<bool>("passBaseline0l_Good");
        const auto& passBaseline1l_Good = tr.getVar<bool>("passBaseline1l_Good");
        const auto& passBaseline2l_Good = tr.getVar<bool>("passBaseline2l_Good");       
        const auto& ntops               = tr.getVar<int>("ntops");       
        const auto& dR_bjets            = tr.getVar<double>("dR_bjets"); 
        const auto& isSignal            = tr.getVar<bool>("isSignal");
        const auto& StopMass            = isSignal ? std::stoi(filetag.substr(filetag.find("-") + 1)) : 0; 

        // -------------------------------
        // -- MT2 hemispheres variables
        // -------------------------------
        
        std::vector<std::string> sufs = {"_TaggedTop", "_0l", "_0l_cm"};

        for ( auto& suf : sufs ) {
        
         auto MT2                = tr.getVar<double>("MT2" + suf);
        
         auto stop1              = tr.getVar<TLorentzVector>("stop1_PtRank" + suf);
         auto stop2              = tr.getVar<TLorentzVector>("stop2_PtRank" + suf);
         auto Jets               = tr.getVec<TLorentzVector>("Jets");        
         auto StopJets           = tr.getVec<TLorentzVector>("StopJets"); 

         auto stop1Mass          = stop1.M();
         auto stop1Eta           = stop1.Eta();
         auto stop1Phi           = stop1.Phi();
         auto stop1Pt            = stop1.Pt();
         auto stop2Mass          = stop2.M();
         auto stop2Eta           = stop2.Eta();
         auto stop2Phi           = stop2.Phi();
         auto stop2Pt            = stop2.M();
         auto dR_stop1stop2      = tr.getVar<double>("dR_stop1stop2" + suf);
         auto dPhi_stop1stop2    = tr.getVar<double>("dPhi_stop1stop2" + suf);
         auto seed1              = tr.getVar<TLorentzVector>("seed1" + suf);
         auto seed2              = tr.getVar<TLorentzVector>("seed2" + suf);
         auto dR_seed1seed2      = seed1.DeltaR(seed2);
        
        // -------------------
        // -- Define weight
        // -------------------
        double weight               = 1.0;
        double eventweight          = 1.0;
        double bTagScaleFactor      = 1.0;
        double prefiringScaleFactor = 1.0;
        double puScaleFactor        = 1.0;
        if(runtype == "MC")
        {
            //Define Lumi weight
            const auto& Weight      = tr.getVar<double>("Weight");
            const auto& lumi        = tr.getVar<double>("Lumi");

            eventweight             = lumi*Weight;

            bTagScaleFactor         = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            prefiringScaleFactor    = tr.getVar<double>("prefiringScaleFactor");
            puScaleFactor           = tr.getVar<double>("puWeightCorr");
            
            weight *= eventweight*bTagScaleFactor*prefiringScaleFactor*puScaleFactor;
        }
                                                                                                                         
        // -------------------------------------------------
        // -- Make cuts and fill histograms here & cutmap
        // -------------------------------------------------
        const std::map<std::string, bool>& cutmap
        {
            {"baseline_0l",         passBaseline0l_Good && ntops >= 2 && dR_bjets >= 1},
            {"baseline_0l_dRle2",   passBaseline0l_Good && ntops >= 2 && dR_bjets >= 1 && dR_seed1seed2 <= 2},
            {"baseline_0l_dRge2",   passBaseline0l_Good && ntops >= 2 && dR_bjets >= 1 && dR_seed1seed2 >= 2},
            {"baseline_0l_largeStop_Signal",         passBaseline0l_Good && ntops >= 2 && dR_bjets >= 1 && (stop1Mass >= StopMass + 400 || stop2Mass >= StopMass + 400) && isSignal},
            {"baseline_0l_largeStop_Background",         passBaseline0l_Good && ntops >= 2 && dR_bjets >= 1 && (stop1Mass >= 1000 || stop2Mass >= 1000) && !isSignal}
        };

        if (!inithisto) 
        {
            InitHistos(cutmap);
            inithisto = true;
        }
       
        my_histos["EventCounter"]->Fill( eventCounter );
 
        for (const auto& cutVar: cutmap) 
        {
            if (cutVar.second) 
            {

                my_histos["h_MT2" + suf + "_" + cutVar.first]->Fill( MT2, weight );
                my_histos["h_stop1Mass" + suf + "_" + cutVar.first]->Fill( stop1Mass, weight );
                my_histos["h_stop1Eta" + suf + "_" + cutVar.first]->Fill( stop1Eta, weight );
                my_histos["h_stop1Phi" + suf + "_" + cutVar.first]->Fill( stop1Phi, weight );
                my_histos["h_stop1Pt" + suf + "_" + cutVar.first]->Fill( stop1Pt, weight );
                my_histos["h_stop2Mass" + suf + "_" + cutVar.first]->Fill( stop2Mass, weight );
                my_histos["h_stop2Eta" + suf + "_" + cutVar.first]->Fill( stop2Eta, weight );
                my_histos["h_stop2Phi" + suf + "_" + cutVar.first]->Fill( stop2Phi, weight );
                my_histos["h_stop2Pt" + suf + "_" + cutVar.first]->Fill( stop2Pt, weight );
               
                my_histos["h_seed1Pt" + suf + "_" + cutVar.first]->Fill( seed1.Pt(), weight );
                my_histos["h_seed1Eta" + suf + "_" + cutVar.first]->Fill( seed1.Eta(), weight );
                my_histos["h_seed1Phi" + suf + "_" + cutVar.first]->Fill( seed1.Phi(), weight );
                my_histos["h_seed1Mass" + suf + "_" + cutVar.first]->Fill( seed1.M(), weight );
                my_histos["h_seed2Pt" + suf + "_" + cutVar.first]->Fill( seed2.Pt(), weight );
                my_histos["h_seed2Eta" + suf + "_" + cutVar.first]->Fill( seed2.Eta(), weight );
                my_histos["h_seed2Phi" + suf + "_" + cutVar.first]->Fill( seed2.Phi(), weight );
                my_histos["h_seed2Mass" + suf + "_" + cutVar.first]->Fill( seed2.M(), weight );
                my_histos["h_dR_seed1seed2" + suf + "_" + cutVar.first]->Fill( dR_seed1seed2, weight );
                
                my_histos["h_dR_stop1stop2" + suf + "_" + cutVar.first]->Fill( dR_stop1stop2, weight );
                my_histos["h_dPhi_stop1stop2" + suf + "_" + cutVar.first]->Fill( dPhi_stop1stop2, weight );
                my_2d_histos["h_Mass_stop1vsstop2" + suf + "_" + cutVar.first]->Fill( stop1Mass, stop2Mass, weight );
                my_2d_histos["h_Mass_stop1vsstop2" + suf + "_" + cutVar.first]->GetXaxis()->SetTitle("M_{#tildet}_{1} [GeV]");
                my_2d_histos["h_Mass_stop1vsstop2" + suf + "_" + cutVar.first]->GetYaxis()->SetTitle("M_{#tildet}_{2} [GeV]");
                my_2d_histos["h_Eta_stop1vsstop2" + suf + "_" + cutVar.first]->Fill( stop1Eta, stop2Eta, weight );
                my_2d_histos["h_Eta_stop1vsstop2" + suf + "_" + cutVar.first]->GetXaxis()->SetTitle("#eta_{#tildet}_{1}");
                my_2d_histos["h_Eta_stop1vsstop2" + suf + "_" + cutVar.first]->GetYaxis()->SetTitle("#eta_{#tildet}_{2}");
                my_2d_histos["h_Phi_stop1vsstop2" + suf + "_" + cutVar.first]->Fill( stop1Phi, stop2Phi, weight );
                my_2d_histos["h_Phi_stop1vsstop2" + suf + "_" + cutVar.first]->GetXaxis()->SetTitle("#phi_{#tildet}_{1}");
                my_2d_histos["h_Phi_stop1vsstop2" + suf + "_" + cutVar.first]->GetYaxis()->SetTitle("#phi_{#tildet}_{2}");
                my_2d_histos["h_Pt_stop1vsstop2" + suf + "_" + cutVar.first]->Fill( stop1Pt, stop2Pt, weight );
                my_2d_histos["h_Pt_stop1vsstop2" + suf + "_" + cutVar.first]->GetXaxis()->SetTitle("pT_{#tildet}_{1}");
                my_2d_histos["h_Pt_stop1vsstop2" + suf + "_" + cutVar.first]->GetYaxis()->SetTitle("pT_{#tildet}_{2}");
                my_2d_histos["h_MT2VsNJets" + suf + "_" + cutVar.first]->Fill( MT2,  NGoodJets_pt45, weight );
                my_2d_histos["h_MT2VsNJets" + suf + "_" + cutVar.first]->GetXaxis()->SetTitle("MT2");
                my_2d_histos["h_MT2VsNJets" + suf + "_" + cutVar.first]->GetYaxis()->SetTitle("N_{J}");
                my_2d_histos["h_MT2Vsstop1Mass" + suf + "_" + cutVar.first]->Fill( MT2, stop1Mass, weight );
                my_2d_histos["h_MT2Vsstop1Mass" + suf + "_" + cutVar.first]->GetXaxis()->SetTitle("MT2");
                my_2d_histos["h_MT2Vsstop1Mass" + suf + "_" + cutVar.first]->GetYaxis()->SetTitle("M_{#tildet}_{1} [GeV]");
                my_2d_histos["h_MT2Vsstop2Mass" + suf + "_" + cutVar.first]->Fill( MT2, stop2Mass, weight );
                my_2d_histos["h_MT2Vsstop2Mass" + suf + "_" + cutVar.first]->GetXaxis()->SetTitle("MT2");
                my_2d_histos["h_MT2Vsstop2Mass" + suf + "_" + cutVar.first]->GetYaxis()->SetTitle("M_{#tildet}_{2} [GeV]");
                my_2d_histos["h_stop1MassVsNJets" + suf + "_" + cutVar.first]->Fill( stop1Mass, NGoodJets_pt45, weight );
                my_2d_histos["h_stop1MassVsNJets" + suf + "_" + cutVar.first]->GetXaxis()->SetTitle("M_{#tildet}_{1}");
                my_2d_histos["h_stop1MassVsNJets" + suf + "_" + cutVar.first]->GetYaxis()->SetTitle("N_{J}");
                my_2d_histos["h_stop2MassVsNJets" + suf + "_" + cutVar.first]->Fill( stop2Mass, NGoodJets_pt45, weight );
                my_2d_histos["h_stop2MassVsNJets" + suf + "_" + cutVar.first]->GetXaxis()->SetTitle("M_{#tildet}_{2}");
                my_2d_histos["h_stop2MassVsNJets" + suf + "_" + cutVar.first]->GetYaxis()->SetTitle("N_{J}");
         
            }
        }
        }
    } 
} 

void StealthHemispheres::WriteHistos(TFile* outfile)
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
