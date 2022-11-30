#define AnalyzeGenStop_cxx
#include "Analyzer/Analyzer/include/AnalyzeGenStop.h"
#include "Framework/Framework/include/Utility.h"
#include "NTupleReader/include/NTupleReader.h"

#include <iostream>

AnalyzeGenStop::AnalyzeGenStop() : initHistos(false)

{
}

//Define all your histograms here. 
void AnalyzeGenStop::InitHistos(const std::map<std::string, bool>& cutMap, const std::vector<TH1DInfo>& histInfos, const std::vector<TH2DInfo>& hist2DInfos)
{

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    //This event counter histogram is necessary so that we know that all the condor jobs ran successfully. If not, when you use the hadder script, you will see a discrepancy in red as the files are being hadded.
    my_histos.emplace( "EventCounter", std::make_shared<TH1D>( "EventCounter", "EventCounter", 2, -1.1, 1.1 ) ) ;

    for(auto& mycut : cutMap)
    {
        for(const auto& hInfo : histInfos)
        { 
            my_histos.emplace(hInfo.name+mycut.first, 
                              std::make_shared<TH1D>((hInfo.name+mycut.first).c_str(),(hInfo.name+mycut.first).c_str(), hInfo.nBins, hInfo.low, hInfo.high));
        }

        for(const auto& h2dInfo : hist2DInfos)
        {
            my_2d_histos.emplace(h2dInfo.name+mycut.first, 
                                 std::make_shared<TH2D>((h2dInfo.name+mycut.first).c_str(),(h2dInfo.name+mycut.first).c_str(), 
                                                        h2dInfo.nBinsX, h2dInfo.lowX, h2dInfo.highX, h2dInfo.nBinsY, h2dInfo.lowY, h2dInfo.highY));
        }
    }
}

//Put everything you want to do per event here.
void AnalyzeGenStop::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    while( tr.getNextEvent() )
    {
        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() & (1000 == 0) ) printf( " Event %i\n", tr.getEvtNum() );

        const auto& eventCounter        = tr.getVar<int>("eventCounter");

        const auto& runtype             = tr.getVar<std::string>("runtype");     

        const auto& NGoodLeptons        = tr.getVar<int>("NGoodLeptons");
        const auto& GoodLeptons         = tr.getVec<std::pair<std::string, utility::LorentzVector>>("GoodLeptons");

        const auto& recoStopMass1       = tr.getVar<utility::LorentzVector>("GM_Stop1");
        const auto& recoStopMass2       = tr.getVar<utility::LorentzVector>("GM_Stop2");

        const auto& Stop1pdgs       = tr.getVec<int>("GM_Stop1_pdgs");
        const auto& Stop2pdgs       = tr.getVec<int>("GM_Stop2_pdgs");

        const auto& Stop1mom       = tr.getVec<int>("GM_Stop1_mom");
        const auto& Stop2mom       = tr.getVec<int>("GM_Stop2_mom");

        const auto& Stop1genetas       = tr.getVec<double>("GM_Stop1_genetas");
        const auto& Stop2genetas       = tr.getVec<double>("GM_Stop2_genetas");

        const auto& Stop1recetas       = tr.getVec<double>("GM_Stop1_recetas");
        const auto& Stop2recetas       = tr.getVec<double>("GM_Stop2_recetas");

        const auto& Stop1genpts       = tr.getVec<double>("GM_Stop1_genpts");
        const auto& Stop2genpts       = tr.getVec<double>("GM_Stop2_genpts");

        const auto& Stop1recphis       = tr.getVec<double>("GM_Stop1_recphis");
        const auto& Stop2recphis       = tr.getVec<double>("GM_Stop2_recphis");

        const auto& Stop1genphis       = tr.getVec<double>("GM_Stop1_genphis");
        const auto& Stop2genphis       = tr.getVec<double>("GM_Stop2_genphis");

        const auto& Stop1DR       = tr.getVec<double>("GM_Stop1_DR");
        const auto& Stop2DR       = tr.getVec<double>("GM_Stop2_DR");

        const auto& Stop1PT       = tr.getVec<double>("GM_Stop1_PT");
        const auto& Stop2PT       = tr.getVec<double>("GM_Stop2_PT");

        const auto& Stop1MassPtRank    = tr.getVar<float>("stop1_ptrank_mass");
        const auto& Stop2MassPtRank    = tr.getVar<float>("stop2_ptrank_mass");
        const auto& Stop1MassMassRank  = tr.getVar<float>("stop1_mrank_mass");
        const auto& Stop2MassMassRank  = tr.getVar<float>("stop2_mrank_mass");
        const auto& StopMass           = tr.getVar<double>("stop_avemass");

        const auto& passBaseline1l      = tr.getVar<bool>("passBaseline1l_Good");
        const auto& passBaseline0l      = tr.getVar<bool>("passBaseline0l_Good");

        const auto& passTrigger               = tr.getVar<bool>("passTrigger");
        const auto& passTriggerMC             = tr.getVar<bool>("passTriggerMC");
        const auto& passMETFilters            = tr.getVar<bool>("passMETFilters");
        const auto& passMadHT                 = tr.getVar<bool>("passMadHT");
        const auto& passHEMVeto               = tr.getVar<bool>("passElectronHEMveto");
        bool pass_general    = passTriggerMC && passTrigger && passMadHT && passMETFilters && passHEMVeto;

        // ------------------------
        // -- Define weight
        // ------------------------
        double weight               = 1.0;
        double eventweight          = 1.0;
        double leptonScaleFactor    = 1.0;
        double bTagScaleFactor      = 1.0;
        double prefiringScaleFactor = 1.0;
        double puScaleFactor        = 1.0;
        
        if(runtype == "MC")
        {
            if( !passMadHT ) continue; //Make sure not to double count DY events
            // Define Lumi weight
            const auto& Weight  = tr.getVar<float>("Weight");
            const auto& lumi = tr.getVar<double>("FinalLumi");
            eventweight = lumi*Weight;
            
            // Define lepton weight
            if(NGoodLeptons == 1)
            {
                const auto& eleLepWeight = tr.getVar<double>("totGoodElectronSF");
                const auto& muLepWeight  = tr.getVar<double>("totGoodMuonSF");
                leptonScaleFactor = (GoodLeptons[0].first == "e") ? eleLepWeight : muLepWeight;
            }
            
            bTagScaleFactor   = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor");
            puScaleFactor = tr.getVar<double>("puWeightCorr");
            
                    weight *= eventweight*leptonScaleFactor*bTagScaleFactor*prefiringScaleFactor*puScaleFactor;
        }

        //This is added to count the number of events- do not change the next two lines.
        std::vector<TH1DInfo> histInfos = {
            {"h_stop1_Wdau_Gen_Match_dR", 720, 0.0, 5.0},
            {"h_stop2_Wdau_Gen_Match_dR", 720, 0.0, 5.0},
            {"h_stop1_Ndau_Gen_Match_dR", 720, 0.0, 5.0},
            {"h_stop2_Ndau_Gen_Match_dR", 720, 0.0, 5.0},
            {"h_stop1_tdau_Gen_Match_dR", 720, 0.0, 5.0},
            {"h_stop2_tdau_Gen_Match_dR", 720, 0.0, 5.0},
            {"h_stop1_Gen_Match_dR", 720, 0.0, 5.0},
            {"h_stop2_Gen_Match_dR", 720, 0.0, 5.0},
            {"h_stop1_Wdau_Gen_Match_PtRatio", 720, 0.0, 10.0},
            {"h_stop2_Wdau_Gen_Match_PtRatio", 720, 0.0, 10.0},
            {"h_stop1_Ndau_Gen_Match_PtRatio", 720, 0.0, 10.0},
            {"h_stop2_Ndau_Gen_Match_PtRatio", 720, 0.0, 10.0},
            {"h_stop1_tdau_Gen_Match_PtRatio", 720, 0.0, 10.0},
            {"h_stop2_tdau_Gen_Match_PtRatio", 720, 0.0, 10.0},
            {"h_stop1_Gen_Match_PtRatio", 720, 0.0, 10.0},
            {"h_stop2_Gen_Match_PtRatio", 720, 0.0, 10.0},
            {"h_stop1_ptrank_mass", 5000, 0, 5000},
            {"h_stop2_ptrank_mass", 5000, 0, 5000},
            {"h_stop1_mrank_mass", 5000, 0, 5000},
            {"h_stop2_mrank_mass", 5000, 0, 5000},
            {"h_stop_avemass", 5000, 0, 5000},
        };

        std::vector<TH2DInfo> hist2DInfos = {
            { "h_stop1_Wdau_Gen_Match_dR_PtRatio", 720, 0.0, 5.0, 720, 0.0, 5.0},
            { "h_stop2_Wdau_Gen_Match_dR_PtRatio", 720, 0.0, 5.0, 720, 0.0, 5.0},
            { "h_stop1_Ndau_Gen_Match_dR_PtRatio", 720, 0.0, 5.0, 720, 0.0, 5.0},
            { "h_stop2_Ndau_Gen_Match_dR_PtRatio", 720, 0.0, 5.0, 720, 0.0, 5.0},
            { "h_stop1_tdau_Gen_Match_dR_PtRatio", 720, 0.0, 5.0, 720, 0.0, 5.0},
            { "h_stop2_tdau_Gen_Match_dR_PtRatio", 720, 0.0, 5.0, 720, 0.0, 5.0},
            { "h_stop1_Gen_Match_dR_PtRatio", 720, 0.0, 5.0, 720, 0.0, 5.0},
            { "h_stop2_Gen_Match_dR_PtRatio", 720, 0.0, 5.0, 720, 0.0, 5.0},
            { "h_stop1_ptrank_mass_Gen_Match_dR", 5000, 0.0, 5000.0, 720, 0.0, 5.0},
            { "h_stop2_ptrank_mass_Gen_Match_dR", 5000, 0.0, 5000.0, 720, 0.0, 5.0},
            { "h_stop1_mrank_mass_Gen_Match_dR", 5000, 0.0, 5000.0, 720, 0.0, 5.0},
            { "h_stop2_mrank_mass_Gen_Match_dR", 5000, 0.0, 5000.0, 720, 0.0, 5.0},
            { "h_stop1_ptrank_mass_Gen_Match_PtRatio", 5000, 0.0, 5000.0, 720, 0.0, 10.0},
            { "h_stop2_ptrank_mass_Gen_Match_PtRatio", 5000, 0.0, 5000.0, 720, 0.0, 10.0},
            { "h_stop1_mrank_mass_Gen_Match_PtRatio", 5000, 0.0, 5000.0, 720, 0.0, 10.0},
            { "h_stop2_mrank_mass_Gen_Match_PtRatio", 5000, 0.0, 5000.0, 720, 0.0, 10.0},
            { "h_stop1_gendau_eta_phi", 720, -6.0, 6.0, 720, -4.0, 4.0},
            { "h_stop2_gendau_eta_phi", 720, -6.0, 6.0, 720, -4.0, 4.0},
            { "h_stop1_recdau_eta_phi", 720, -6.0, 6.0, 720, -4.0, 4.0},
            { "h_stop2_recdau_eta_phi", 720, -6.0, 6.0, 720, -4.0, 4.0},

        };

        const std::map<std::string, bool> cut_map 
        {
            {"_0L"               , pass_general && passBaseline0l                                                      },                         
            {"_1L"               , pass_general && passBaseline1l                                                      },                         
        };

        // Initialize Histograms
        if(!initHistos)
        {
            InitHistos(cut_map, histInfos, hist2DInfos);
            initHistos = true;

        }

        my_histos["EventCounter"]->Fill( eventCounter );

        for(auto& kv : cut_map)
        {
            if(kv.second)
            {
                double w = weight;
                for (unsigned int p = 0; p < Stop1genpts.size(); p++)
                {
                    if (abs(Stop1mom[p]) == 24)
                    {
                        my_histos["h_stop1_Wdau_Gen_Match_dR"+kv.first]->Fill(Stop1DR[p], w);
                        my_histos["h_stop1_Wdau_Gen_Match_PtRatio"+kv.first]->Fill(Stop1PT[p], w);
                        my_2d_histos["h_stop1_Wdau_Gen_Match_dR_PtRatio"+kv.first]->Fill(Stop1DR[p], Stop1PT[p], w);
                    } else if (abs(Stop1mom[p]) == 1000022)
                    {
                        my_histos["h_stop1_Ndau_Gen_Match_dR"+kv.first]->Fill(Stop1DR[p], w);
                        my_histos["h_stop1_Ndau_Gen_Match_PtRatio"+kv.first]->Fill(Stop1PT[p], w);
                        my_2d_histos["h_stop1_Ndau_Gen_Match_dR_PtRatio"+kv.first]->Fill(Stop1DR[p], Stop1PT[p], w);
                    } else if (abs(Stop1mom[p]) == 6)
                    {
                        my_histos["h_stop1_tdau_Gen_Match_dR"+kv.first]->Fill(Stop1DR[p], w);
                        my_histos["h_stop1_tdau_Gen_Match_PtRatio"+kv.first]->Fill(Stop1PT[p], w);
                        my_2d_histos["h_stop1_tdau_Gen_Match_dR_PtRatio"+kv.first]->Fill(Stop1DR[p], Stop1PT[p], w);
                    }

                    my_histos["h_stop1_Gen_Match_dR"+kv.first]->Fill(Stop1DR[p], w);
                    my_histos["h_stop1_Gen_Match_PtRatio"+kv.first]->Fill(Stop1PT[p], w);
                    my_2d_histos["h_stop1_Gen_Match_dR_PtRatio"+kv.first]->Fill(Stop1DR[p], Stop1PT[p], w);
                    my_2d_histos["h_stop1_ptrank_mass_Gen_Match_PtRatio"+kv.first]->Fill(Stop1MassPtRank, Stop1PT[p], w);
                    my_2d_histos["h_stop1_ptrank_mass_Gen_Match_dR"+kv.first]->Fill(Stop1MassPtRank, Stop1DR[p], w);
                    my_2d_histos["h_stop1_mrank_mass_Gen_Match_PtRatio"+kv.first]->Fill(Stop1MassMassRank, Stop1PT[p], w);
                    my_2d_histos["h_stop1_mrank_mass_Gen_Match_dR"+kv.first]->Fill(Stop1MassMassRank, Stop1DR[p], w);

                    //if (Stop2MassPtRank < 300.0)
                    //{
                        my_2d_histos["h_stop1_gendau_eta_phi"+kv.first]->Fill(Stop1genetas[p], Stop1genphis[p], 5.0);
                        my_2d_histos["h_stop1_recdau_eta_phi"+kv.first]->Fill(Stop1recetas[p], Stop1recphis[p], 10.0);
                    //}
                }

                for (unsigned int p = 0; p < Stop2genpts.size(); p++)
                {
                    if (abs(Stop2mom[p]) == 24)
                    {
                        my_histos["h_stop2_Wdau_Gen_Match_dR"+kv.first]->Fill(Stop2DR[p], w);
                        my_histos["h_stop2_Wdau_Gen_Match_PtRatio"+kv.first]->Fill(Stop2PT[p], w);
                        my_2d_histos["h_stop2_Wdau_Gen_Match_dR_PtRatio"+kv.first]->Fill(Stop2DR[p], Stop2PT[p], w);
                    } else if (abs(Stop2mom[p]) == 1000022)
                    {
                        my_histos["h_stop2_Ndau_Gen_Match_dR"+kv.first]->Fill(Stop2DR[p], w);
                        my_histos["h_stop2_Ndau_Gen_Match_PtRatio"+kv.first]->Fill(Stop2PT[p], w);
                        my_2d_histos["h_stop2_Ndau_Gen_Match_dR_PtRatio"+kv.first]->Fill(Stop2DR[p], Stop2PT[p], w);
                    } else if (abs(Stop2mom[p]) == 6)
                    {
                        my_histos["h_stop2_tdau_Gen_Match_dR"+kv.first]->Fill(Stop2DR[p], w);
                        my_histos["h_stop2_tdau_Gen_Match_PtRatio"+kv.first]->Fill(Stop2PT[p], w);
                        my_2d_histos["h_stop2_tdau_Gen_Match_dR_PtRatio"+kv.first]->Fill(Stop2DR[p], Stop2PT[p], w);
                    }
                    my_histos["h_stop2_Gen_Match_dR"+kv.first]->Fill(Stop2DR[p], w);
                    my_histos["h_stop2_Gen_Match_PtRatio"+kv.first]->Fill(Stop2PT[p], w);
                    my_2d_histos["h_stop2_Gen_Match_dR_PtRatio"+kv.first]->Fill(Stop2DR[p], Stop2PT[p], w);
                    my_2d_histos["h_stop2_ptrank_mass_Gen_Match_PtRatio"+kv.first]->Fill(Stop2MassPtRank, Stop2PT[p], w);
                    my_2d_histos["h_stop2_ptrank_mass_Gen_Match_dR"+kv.first]->Fill(Stop2MassPtRank, Stop2DR[p], w);
                    my_2d_histos["h_stop2_mrank_mass_Gen_Match_PtRatio"+kv.first]->Fill(Stop2MassMassRank, Stop2PT[p], w);
                    my_2d_histos["h_stop2_mrank_mass_Gen_Match_dR"+kv.first]->Fill(Stop2MassMassRank, Stop2DR[p], w);
                
                    //if (Stop2MassPtRank < 300.0)
                    //{
                        my_2d_histos["h_stop2_gendau_eta_phi"+kv.first]->Fill(Stop2genetas[p], Stop2genphis[p], 5.0);
                        my_2d_histos["h_stop2_recdau_eta_phi"+kv.first]->Fill(Stop2recetas[p], Stop2recphis[p], 10.0);
                    //}

                }

                my_histos["h_stop1_ptrank_mass" + kv.first]->Fill(Stop1MassPtRank, w);
                my_histos["h_stop2_ptrank_mass" + kv.first]->Fill(Stop2MassPtRank, w);
                my_histos["h_stop1_mrank_mass" + kv.first]->Fill(Stop1MassMassRank, w);
                my_histos["h_stop2_mrank_mass" + kv.first]->Fill(Stop2MassMassRank, w);
                my_histos["h_stop_avemass" + kv.first]->Fill(StopMass, w);
            }
        }
        
        //Make cuts and fill histograms here
        if ( passBaseline0l and pass_general ) {

            std::cout << "STOP1 PARTICLES:" << std::endl;
            for (unsigned int p=0; p<Stop1pdgs.size(); p++) {
                std::cout << "    " << Stop1mom[p] << " --> " << Stop1pdgs[p] << " PT: " << Stop1genpts[p] << " ETA: " << Stop1genetas[p] << " Gen_Match_dR: " << Stop1DR[p] << " dPT: " << Stop1PT[p] << std::endl;
            }
            std::cout << "STOP2 PARTICLES:" << std::endl;
            for (unsigned int p=0; p<Stop2pdgs.size(); p++) {
                std::cout << "    " << Stop2mom[p] << " --> " << Stop2pdgs[p] << " PT: " << Stop2genpts[p] << " ETA: " << Stop2genetas[p] << " Gen_Match_dR: " << Stop2DR[p] << " dPT: " << Stop2PT[p] << std::endl;
            }

            std::cout << "STOP1: M: " << recoStopMass1.M() << " ETA: " << recoStopMass1.Eta() << " PT: " << recoStopMass1.Pt() << std::endl;
            std::cout << "STOP2: M: " << recoStopMass2.M() << " ETA: " << recoStopMass2.Eta() << " PT: " << recoStopMass2.Pt() << std::endl;
            std::cout << std::endl;
        }
    } 
}

void AnalyzeGenStop::WriteHistos(TFile* outfile)
{
    outfile->cd();

    for (const auto &p : my_histos) {
        p.second->Write();
    }
    
    for (const auto &p : my_2d_histos) {
        p.second->Write();
    }
}
