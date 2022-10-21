#define TopTaggerSF_Analyzer_cxx
#include "Analyzer/Analyzer/include/TopTaggerSF_Analyzer.h"
#include "Framework/Framework/include/SetUpTopTagger.h"
#include "Framework/Framework/include/Utility.h"
#include "NTupleReader/include/NTupleReader.h"

#include <iostream>

TopTaggerSF_Analyzer::TopTaggerSF_Analyzer()
{

    initHistos = false;

    hist1Dinfos = {
        {"h_nElectrons",     6,   -0.5,    5.5},
        {"h_nMuons",         6,   -0.5,    5.5},
        {"h_nJets",          21,  -0.5,   20.5},
        {"h_nBJets",         21,  -0.5,   20.5},
        {"h_nBJetsLoose",    21,  -0.5,   20.5},
        {"h_ht",            720,   0.0, 4000.0},
        {"h_met",           360,   0.0, 1500.0},
        {"h_jet1METdPhi",    80,  -4.0,    4.0},
        {"h_jet2METdPhi",    80,  -4.0,    4.0},
        {"h_jet3METdPhi",    80,  -4.0,    4.0},
        {"h_muonBjetDR",    360,   0.0,   10.0},
        {"h_muonBjetMass",  720,   0.0, 1000.0},
        {"h_muonMETdPhi",    80,  -4.0,    4.0},
        {"h_muonMETmT",     720,   0.0, 1000.0},
        {"h_topPt",         360,   0.0, 1500.0}
    };

    hist2Dinfos = {
        {"h_topPt_vs_HT",    360, 0.0, 1000.0, 360, 0.0, 4000.0}, 
        {"h_topPt_vs_nJets", 360, 0.0, 1000.0, 13, -0.5, 12.5}, 
    };

    njets = {"7incl", "7", "8", "9", "10", "11", "12incl"};
    
    tags = {"num", "den"};

    my_var_suffix = {""};
}

void TopTaggerSF_Analyzer::InitHistos(const std::map<std::string, bool>& cutMap)
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    my_1D_histos.emplace("EventCounter", std::make_shared<TH1D>("EventCounter","EventCounter", 2, -1.1, 1.1 ) );

    for(auto& cutVar : cutMap)
    {
        for (auto& njet : njets)
        {
            for (auto& tag : tags)
            {
                for(const auto& h1Dinfo : hist1Dinfos)
                { 
                    std::string name = h1Dinfo.name + cutVar.first + "_" + tag + "_Njet" + njet;
                    my_1D_histos.emplace(name, std::make_shared<TH1D>((name).c_str(),(name).c_str(), h1Dinfo.nBins, h1Dinfo.low, h1Dinfo.high));
                }
                for(const auto& h2Dinfo : hist2Dinfos)
                { 
                    std::string name = h2Dinfo.name + cutVar.first + "_" + tag + "_Njet" + njet;
                    my_2D_histos.emplace(name, std::make_shared<TH2D>((name).c_str(),(name).c_str(), h2Dinfo.nBinsX, h2Dinfo.lowX, h2Dinfo.highX, h2Dinfo.nBinsY, h2Dinfo.lowY, h2Dinfo.highY));
                }
            }
        }
    }
}

void TopTaggerSF_Analyzer::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    while (tr.getNextEvent())
    {
        const auto& eventCounter           = tr.getVar<int>("eventCounter");

        if ( maxevents != -1 && tr.getEvtNum() >= maxevents )   
            break;

        if ( tr.getEvtNum() % 1000 == 0)
            printf( " Event %i\n", tr.getEvtNum() );        

        const auto& runtype        = tr.getVar<std::string>("runtype");

        const auto& bestTopPt   = tr.getVar<double>("bestTopPt");
        const auto& bestTopDisc = tr.getVar<double>("bestTopDisc");

        // Event-level booleans (triggers, filters, etc.)
        const auto& passMETFilters  = tr.getVar<bool>("passMETFilters");
        const auto& passMadHT       = tr.getVar<bool>("passMadHT");
        const auto& passMuonTrigger = tr.getVar<bool>("passTriggerMuon");
        const auto& passHadTrigger  = tr.getVar<bool>("passTriggerAllHad");

        // Lepton quantities
        const auto& Muons          = tr.getVec<utility::LorentzVector>("Muons");
        const auto& GoodMuons      = tr.getVec<bool>("GoodMuons");
        const auto& NGoodMuons     = tr.getVar<int>("NGoodMuons");
        const auto& NGoodElectrons = tr.getVar<int>("NGoodElectrons");

        // Jet quantities
        const auto& Jets                  = tr.getVec<utility::LorentzVector>("Jets");
        const auto& JetID                 = tr.getVar<bool>("JetID");        
        const auto& GoodJets_pt30         = tr.getVec<bool>("GoodJets_pt30");
        const auto& NGoodJets_pt30        = tr.getVar<int>("NGoodJets_pt30");
        const auto& NGoodBJets_pt30       = tr.getVar<int>("NGoodBJets_pt30");
        const auto& GoodBJets_pt30_loose  = tr.getVec<bool>("GoodBJets_pt30_loose");
        const auto& GoodBJets_pt30        = tr.getVec<bool>("GoodBJets_pt30");
        const auto& NGoodBJets_pt30_loose = tr.getVar<int>("NBJets_pt30_loose");

        // Event-level quantities
        const auto& HT_trigger_pt30 = tr.getVar<double>("HT_trigger_pt30");
        const auto& HT              = tr.getVar<float>("HT");
        const auto& MET             = tr.getVar<float>("MET"); 
        const auto& METPhi          = tr.getVar<float>("METPhi");

        // Get the 4-vec for the MET
        utility::LorentzVector lvMET;
        lvMET.SetPt(MET); lvMET.SetEta(0.0); lvMET.SetPhi(METPhi); lvMET.SetE(MET);

        bool pass_ZeroMuons     = NGoodMuons     == 0;
        bool pass_ZeroElectrons = NGoodElectrons == 0;

        bool pass_SingleMuon    = NGoodMuons            == 1;
        bool pass_ExtraLepVeto  = NGoodElectrons        == 0;
        bool pass_4jets         = NGoodJets_pt30        >= 4;
        bool pass_7jets         = NGoodJets_pt30        >= 7;
        bool pass_1bjet         = NGoodBJets_pt30       >= 1;

        // The medium b jet (above) also passes the loose working point
        // So if we want an additional b jet that is loose we will have at least 2
        bool pass_1bjetLoose    = NGoodBJets_pt30_loose >= 2;
        bool pass_HT200         = HT                     > 200.0;
        bool pass_HT500         = HT_trigger_pt30        > 500.0;
        bool pass_HT1000        = HT                     > 1000.0;
        bool pass_MET50         = MET                    > 50.0;

        // Here pt ordered jet collection is assumed...
        // Calculate dPhi between MET and leading three jets
        unsigned int jetOrder = 1;
        bool pass_Jet1METdPhi  = false; double Jet1METdPhi = -999.0;
        bool pass_Jet2METdPhi  = false; double Jet2METdPhi = -999.0;
        bool pass_Jet3METdPhi  = false; double Jet3METdPhi = -999.0;
        bool pass_JetMETdPhi   = false;

        bool pass_MuonBjetdR     = false; double MuonBjetdR   = 9999.0;
        bool pass_MuonBjetMass   = false; double MuonBjetMass = -1.0;
        bool pass_MuonBjetdRMass = false;
        auto& goodjets_pt30 = tr.createDerivedVec<TLorentzVector>("GoodJets_pt30_tlv");
        for (unsigned int iJet = 0; iJet < Jets.size(); iJet++)
        {
            if (GoodJets_pt30.at(iJet))
            {
                goodjets_pt30.emplace_back(utility::convertLV<TLorentzVector>(Jets.at(iJet)));

                if (jetOrder == 1)
                {
                    Jet1METdPhi      = utility::DeltaPhi(Jets.at(iJet), lvMET);
                    pass_Jet1METdPhi = fabs(Jet1METdPhi) > 0.5;
                    jetOrder++;
                }
                else if (jetOrder == 2)
                {
                    Jet2METdPhi      = utility::DeltaPhi(Jets.at(iJet), lvMET);
                    pass_Jet2METdPhi = fabs(Jet2METdPhi) > 0.5;
                    jetOrder++;
                }
                else if (jetOrder == 3)
                {
                    Jet3METdPhi      = utility::DeltaPhi(Jets.at(iJet), lvMET);
                    pass_Jet3METdPhi = fabs(Jet3METdPhi) > 0.3;
                    jetOrder++;
                }
            } 

            // Check for another b jet that is loose but not also the medium one
            if (GoodBJets_pt30_loose.at(iJet) && !GoodBJets_pt30.at(iJet))
            {
                for (unsigned int iMuon = 0; iMuon < Muons.size(); iMuon++)
                {
                    if (!GoodMuons.at(iMuon))
                        continue;

                    double tempMuonBjetdR = utility::DeltaR(Jets.at(iJet), Muons.at(iMuon));

                    // At least one loose b jet be within a dR of 1.5
                    pass_MuonBjetdR |= (tempMuonBjetdR < 1.5);
    
                    auto combinedMuonBjet   = Jets.at(iJet) + Muons.at(iMuon);
                    double tempMuonBjetMass = combinedMuonBjet.M();
    
                    if (tempMuonBjetdR < MuonBjetdR)
                    {
                        MuonBjetdR   = tempMuonBjetdR;
                        MuonBjetMass = tempMuonBjetMass;
                    }

                    // At least one loose b jet + muon combination to give inv mass inside window around top---the tag
                    pass_MuonBjetMass |= (tempMuonBjetMass >= 30.0 && tempMuonBjetMass <= 180.0);

                    // Require that the loose b jet satisfies both criteria at same time
                    // E.g. one loose b cannot satisfy dR while another loose b satisfies invariant mass 
                    pass_MuonBjetdRMass |= ((tempMuonBjetdR < 1.5) && pass_MuonBjetMass);
                }
            }
        }

        pass_JetMETdPhi = (pass_Jet1METdPhi && pass_Jet2METdPhi && pass_Jet3METdPhi);

        // Put the good muon together with MET to calculate dPhi and transverse mass 
        bool pass_MuonMETdPhi = false; double MuonMETdPhi = 9999.0;
        bool pass_MuonMETmT   = false; double MuonMETmT   = 9999.0;
        for (unsigned int iMuon = 0; iMuon < Muons.size(); iMuon++)
        {
            if (!GoodMuons.at(iMuon))
                continue;

            double tempMuonMETdPhi  = utility::DeltaPhi(Muons.at(iMuon), lvMET);
            double tempMuonMETmT    = utility::calcMT(Muons.at(iMuon),   lvMET);

            if (fabs(tempMuonMETdPhi) < fabs(MuonMETdPhi))
                MuonMETdPhi = tempMuonMETdPhi;
            if (fabs(tempMuonMETmT) < fabs(MuonMETmT))
                MuonMETmT = tempMuonMETmT;

            pass_MuonMETdPhi |= (fabs(tempMuonMETdPhi) < 0.8);
            pass_MuonMETmT   |= (tempMuonMETmT         < 100.0);
        }

        // -------------------
        // -- Define weight
        // -------------------
        double weightQCD              = 1.0;
        double weightTTbar            = 1.0;
        double eventweight            = 1.0;
        double bTagScaleFactor        = 1.0;
        double prefiringScaleFactor   = 1.0;
        double puScaleFactor          = 1.0;
        double totGoodMuonScaleFactor = 1.0;
        if(runtype == "MC")
        {
            // Define Lumi weight
            const auto& Weight = tr.getVar<float>("Weight");
            const auto& lumi   = tr.getVar<double>("FinalLumi");
            eventweight        = lumi*Weight;

            bTagScaleFactor        = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            prefiringScaleFactor   = tr.getVar<double>("prefiringScaleFactor");
            puScaleFactor          = tr.getVar<double>("puWeightCorr");
            totGoodMuonScaleFactor = tr.getVar<double>("totGoodMuonSF");

            weightTTbar *= eventweight * totGoodMuonScaleFactor * puScaleFactor * bTagScaleFactor * prefiringScaleFactor;
            weightQCD   *= eventweight /* jetTrigSF ? */        * puScaleFactor * bTagScaleFactor * prefiringScaleFactor;
        }

        // To be used with Single-Muon-triggered data and measuring tagging efficiency SF
        bool pass_SemiLepTTbarCR = JetID                                                    &&
                                   passMETFilters                                           &&
                                   passMadHT                                                &&
                                   passMuonTrigger                                          &&
                                   pass_SingleMuon                                          &&
                                   pass_7jets                                               &&
                                   pass_ExtraLepVeto                                        &&
                                   pass_1bjet                                               &&
                                   //pass_Jet1METdPhi && pass_Jet2METdPhi && pass_Jet3METdPhi &&
                                   pass_1bjetLoose                                          &&
                                   pass_MuonBjetdRMass                                      &&
                                   pass_MuonMETdPhi                                         &&
                                   pass_MuonMETmT                                           &&
                                   pass_HT500                                               ; 
                                   //pass_MET50                                               ;

        bool pass_SemiLepTTbarCR_Joe = JetID                                                    &&
                                       passMETFilters                                           &&
                                       passMadHT                                                &&
                                       passMuonTrigger                                          &&
                                       pass_SingleMuon                                          &&
                                       pass_4jets                                               &&
                                       pass_ExtraLepVeto                                        &&
                                       pass_1bjet                                               &&
                                       pass_Jet1METdPhi && pass_Jet2METdPhi && pass_Jet3METdPhi &&
                                       pass_1bjetLoose                                          &&
                                       pass_MuonBjetdRMass                                      &&
                                       pass_MuonMETdPhi                                         &&
                                       pass_MuonMETmT                                           &&
                                       pass_HT200                                               &&
                                       pass_MET50                                               ;

        // For use with JetHT-triggered data and measuring mistag SF
        bool pass_QCDCR = JetID              && 
                          passMETFilters     &&
                          passMadHT          &&
                          pass_7jets         &&
                          pass_ZeroMuons     &&
                          pass_ZeroElectrons &&
                          pass_HT500;

        bool pass_QCDCR_Joe = JetID              && 
                              passMETFilters     &&
                              passMadHT          &&
                              pass_4jets         &&
                              pass_ZeroMuons     &&
                              pass_ZeroElectrons &&
                              pass_HT1000;

        const std::map<std::string, bool> cut_map
        {
            {                                                                   "", true},
            {                                                               "_1mu", passMuonTrigger && pass_SingleMuon},
            {                                                            "_1mu_7j", passMuonTrigger && pass_SingleMuon && pass_7jets},
            {                                                        "_1mu_7j_0el", passMuonTrigger && pass_SingleMuon && pass_7jets && pass_ExtraLepVeto},
            {                                                     "_1mu_7j_0el_1b", passMuonTrigger && pass_SingleMuon && pass_7jets && pass_ExtraLepVeto && pass_1bjet},
            {                                            "_1mu_7j_0el_1b_jMETdPhi", passMuonTrigger && pass_SingleMuon && pass_7jets && pass_ExtraLepVeto && pass_1bjet && pass_JetMETdPhi},
            {                                        "_1mu_7j_0el_1b_jMETdPhi_1bl", passMuonTrigger && pass_SingleMuon && pass_7jets && pass_ExtraLepVeto && pass_1bjet && pass_JetMETdPhi && pass_1bjetLoose},
            {                             "_1mu_7j_0el_1b_jMETdPhi_1bl_mubdR_mubM", passMuonTrigger && pass_SingleMuon && pass_7jets && pass_ExtraLepVeto && pass_1bjet && pass_JetMETdPhi && pass_1bjetLoose && pass_MuonBjetdRMass},
            {                    "_1mu_7j_0el_1b_jMETdPhi_1bl_mubdR_mubM_muMETphi", passMuonTrigger && pass_SingleMuon && pass_7jets && pass_ExtraLepVeto && pass_1bjet && pass_JetMETdPhi && pass_1bjetLoose && pass_MuonBjetdRMass && pass_MuonMETdPhi},
            {            "_1mu_7j_0el_1b_jMETdPhi_1bl_mubdR_mubM_muMETphi_muMETmT", passMuonTrigger && pass_SingleMuon && pass_7jets && pass_ExtraLepVeto && pass_1bjet && pass_JetMETdPhi && pass_1bjetLoose && pass_MuonBjetdRMass && pass_MuonMETdPhi && pass_MuonMETmT},
            {      "_1mu_7j_0el_1b_jMETdPhi_1bl_mubdR_mubM_muMETphi_muMETmT_HT500", passMuonTrigger && pass_SingleMuon && pass_7jets && pass_ExtraLepVeto && pass_1bjet && pass_JetMETdPhi && pass_1bjetLoose && pass_MuonBjetdRMass && pass_MuonMETdPhi && pass_MuonMETmT && pass_HT500},
            {"_1mu_7j_0el_1b_jMETdPhi_1bl_mubdR_mubM_muMETphi_muMETmT_HT500_MET50", passMuonTrigger && pass_SingleMuon && pass_7jets && pass_ExtraLepVeto && pass_1bjet && pass_JetMETdPhi && pass_1bjetLoose && pass_MuonBjetdRMass && pass_MuonMETdPhi && pass_MuonMETmT && pass_HT500 && pass_MET50},

            {"_QCDCR"                                                           , pass_QCDCR},
            {"_QCDCR_Joe"                                                       , pass_QCDCR_Joe},
            {"_TTbarCR"                                                         , pass_SemiLepTTbarCR},
            {"_TTbarCR_Joe"                                                     , pass_SemiLepTTbarCR_Joe},
       };

        // Initialize Histograms
        if(!initHistos)
        {
            InitHistos(cut_map);
            initHistos = true;
        }

        my_1D_histos["EventCounter"]->Fill(eventCounter);

        std::map<std::string, bool> njetsMap =
        {
            {"7",      NGoodJets_pt30 == 7},
            {"7incl",  NGoodJets_pt30 >= 7},
            {"8",      NGoodJets_pt30 == 8},
            {"9",      NGoodJets_pt30 == 9},
            {"10",     NGoodJets_pt30 == 10},
            {"11",     NGoodJets_pt30 == 11},
            {"12incl", NGoodJets_pt30 >= 12},
        };

        std::map<std::string, bool> numDenMap =
        {
            {"num", bestTopDisc > 0.95},
            {"den", bestTopDisc > 0.0},
        };

        for (auto& kv : cut_map)
        {
            for (auto& njetCut : njetsMap)
            {
                for (auto& numDen : numDenMap)
                {
                    if (kv.second && njetCut.second && numDen.second)
                    {

                        double weight = weightQCD;
                        if (kv.first.find("QCD") == std::string::npos)
                            weight = weightTTbar;

                        my_1D_histos["h_nElectrons"   + kv.first + "_" + numDen.first + "_Njet" + njetCut.first]->Fill(NGoodElectrons,        weight);
                        my_1D_histos["h_nMuons"       + kv.first + "_" + numDen.first + "_Njet" + njetCut.first]->Fill(NGoodMuons,            weight);
                        my_1D_histos["h_nJets"        + kv.first + "_" + numDen.first + "_Njet" + njetCut.first]->Fill(NGoodJets_pt30,        weight);
                        my_1D_histos["h_nBJets"       + kv.first + "_" + numDen.first + "_Njet" + njetCut.first]->Fill(NGoodBJets_pt30,       weight);
                        my_1D_histos["h_nBJetsLoose"  + kv.first + "_" + numDen.first + "_Njet" + njetCut.first]->Fill(NGoodBJets_pt30_loose, weight);
                        my_1D_histos["h_ht"           + kv.first + "_" + numDen.first + "_Njet" + njetCut.first]->Fill(HT_trigger_pt30,       weight);
                        my_1D_histos["h_met"          + kv.first + "_" + numDen.first + "_Njet" + njetCut.first]->Fill(MET,                   weight);
                        my_1D_histos["h_jet1METdPhi"  + kv.first + "_" + numDen.first + "_Njet" + njetCut.first]->Fill(Jet1METdPhi,           weight);
                        my_1D_histos["h_jet2METdPhi"  + kv.first + "_" + numDen.first + "_Njet" + njetCut.first]->Fill(Jet2METdPhi,           weight);
                        my_1D_histos["h_jet3METdPhi"  + kv.first + "_" + numDen.first + "_Njet" + njetCut.first]->Fill(Jet3METdPhi,           weight);
                        my_1D_histos["h_muonBjetDR"   + kv.first + "_" + numDen.first + "_Njet" + njetCut.first]->Fill(MuonBjetdR,            weight);
                        my_1D_histos["h_muonBjetMass" + kv.first + "_" + numDen.first + "_Njet" + njetCut.first]->Fill(MuonBjetMass,          weight);
                        my_1D_histos["h_muonMETdPhi"  + kv.first + "_" + numDen.first + "_Njet" + njetCut.first]->Fill(MuonMETdPhi,           weight);
                        my_1D_histos["h_muonMETmT"    + kv.first + "_" + numDen.first + "_Njet" + njetCut.first]->Fill(MuonMETmT,             weight);
                        my_1D_histos["h_topPt"        + kv.first + "_" + numDen.first + "_Njet" + njetCut.first]->Fill(bestTopPt,             weight);

                        my_2D_histos["h_topPt_vs_HT"    + kv.first + "_" + numDen.first + "_Njet" + njetCut.first]->Fill(bestTopPt, HT_trigger_pt30, weight);

                        int correctNjets = NGoodJets_pt30;
                        if (correctNjets >= 12)
                            correctNjets = 12;

                        my_2D_histos["h_topPt_vs_nJets" + kv.first + "_" + numDen.first + "_Njet" + njetCut.first]->Fill(bestTopPt, correctNjets, weight);

                    }
                }
            }
        }
    }
}

void TopTaggerSF_Analyzer::WriteHistos(TFile* outfile)
{
    outfile->cd();

    for (const auto &p : my_1D_histos) {
        p.second->Write();
    }
    
    for (const auto &p : my_2D_histos) {
        p.second->Write();
    }

    outfile->Close();
}
