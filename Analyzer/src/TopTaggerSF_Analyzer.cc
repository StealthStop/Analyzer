#define TopTaggerSF_Analyzer_cxx
#include "Analyzer/Analyzer/include/TopTaggerSF_Analyzer.h"
#include "Framework/Framework/include/Utility.h"
#include "NTupleReader/include/NTupleReader.h"

#include <iostream>

TopTaggerSF_Analyzer::TopTaggerSF_Analyzer()
{
    initHistos = false;

    hist1Dinfos = {
        {"nElectrons",     6,   -0.5,    5.5},
        {"nMuons",         6,   -0.5,    5.5},
        {"nJets",          22,  -0.5,   21.5},
        {"nBJets",         22,  -0.5,   21.5},
        {"nBJetsLoose",    22,  -0.5,   21.5},
        {"ht",            720,   0.0, 7200.0},
        {"met",           360,   0.0, 1440.0},
        {"jet1METdPhi",    80,  -4.0,    4.0},
        {"jet2METdPhi",    80,  -4.0,    4.0},
        {"jet3METdPhi",    80,  -4.0,    4.0},
        {"muonBjetDR",    360,   0.0,   10.0},
        {"muonBjetMass",  720,   0.0, 1440.0},
        {"muonMETdPhi",    80,  -4.0,    4.0},
        {"muonMETmT",     720,   0.0, 1440.0},
        {"topPt",         360,   0.0, 1440.0},

        // Override the binning in InitHistos for variable bin widths
        {"Binned_topPt", -1, -1.0, -1.0},
        {"Binned_nJets", -1, -1.0, -1.0},
        {"Binned_HT",    -1, -1.0, -1.0},
        {"Binned_Disc1", -1, -1.0, -1.0},
        {"Binned_Disc2", -1, -1.0, -1.0},
    };

    hist2Dinfos = {

        // We will override the binning inside of InitHistos in order to variable bin widths
        {"Binned_topPt_vs_nJets", -1, -1.0, -1.0, -1, -1.0, -1.0},

        {"topPt_vs_nJets",    120, 0.0, 1500.0,   5, 5.5,   10.5},
        {"topPt_vs_HT",       120, 0.0, 1500.0, 120, 0.0, 3000.0},
    };

    njets = {"", "6", "7", "8", "9incl"};
    
    // For setting histograms of events passing numerator and denominator selections
    tags = {"", "numR", "numM", "denR", "denM"};
}

Double_t* TopTaggerSF_Analyzer::vecToDouble_t(const std::vector<double>& input)
{
    unsigned int length = input.size();
    
    Double_t* output = new Double_t[length];

    for (unsigned int i = 0; i < length; i++)
        output[i] = input.at(i);

    return output;
}

void TopTaggerSF_Analyzer::InitHistos(const std::map<std::string, bool>& cutMap)
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    my_1D_histos.emplace("EventCounter", std::make_shared<TH1D>("EventCounter","EventCounter", 2, -1.1, 1.1 ) );

    //Define binning for the histograms
    const Int_t nPtBins_Njet6_R = 6;
    std::vector<double> ptBinEdges_Njet6_R = { 0.0, 100.0, 150.0, 200.0, 300.0, 400.0, 1500.0 };

    const Int_t nPtBins_Njet7_R = 6;
    std::vector<double> ptBinEdges_Njet7_R = { 0.0, 100.0, 150.0, 200.0, 300.0, 400.0, 1500.0 };

    const Int_t nPtBins_Njet8_R = 4;
    std::vector<double> ptBinEdges_Njet8_R = { 0.0, 120.0, 200.0, 300.0, 1500.0 };

    const Int_t nPtBins_Njet9_R = 3;
    std::vector<double> ptBinEdges_Njet9_R = { 0.0, 120.0, 200.0, 1500.0 };

    const Int_t nPtBins_Njet6_M = 2;
    std::vector<double> ptBinEdges_Njet6_M = { 400.0, 600.0, 1500.0 };
                                          
    const Int_t nPtBins_Njet7_M = 2;      
    std::vector<double> ptBinEdges_Njet7_M = { 400.0, 600.0, 1500.0 };
                                          
    const Int_t nPtBins_Njet8_M = 1;      
    std::vector<double> ptBinEdges_Njet8_M = { 400.0, 1500.0 };
                                          
    const Int_t nPtBins_Njet9_M = 1;      
    std::vector<double> ptBinEdges_Njet9_M = { 400.0, 1500.0 };

    const Int_t nJetBins = 4;
    Double_t njetBinEdges[ nJetBins + 1 ] = { 5.5, 6.5, 7.5, 8.5, 9.5 };

    const Int_t nHTBins = 4;
    Double_t htBinEdges[ nHTBins + 1 ] = { 200.0, 450.0, 600.0, 850.0, 3000.0 };

    const Int_t nDiscBins = 7;
    Double_t discBinEdges[ nDiscBins + 1 ] = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0 };

    for(auto& cutVar : cutMap)
    {
        for (auto& njet : njets)
        {
            for (auto& tag : tags)
            {
                std::string njetStr = "";
                if (njet != "")
                    njetStr = "_Njet" + njet;

                std::string tagStr = "";
                if (tag != "")
                    tagStr = "_" + tag;

                Int_t nPtBins;
                std::vector<double> ptBinEdges;
                // Choose binning for resolved case
                if (tagStr.find("R") != std::string::npos)
                {
                    if (njetStr.find("6") != std::string::npos)
                    {
                        nPtBins    = nPtBins_Njet6_R;
                        ptBinEdges = ptBinEdges_Njet6_R; 
                    }
                    else if (njetStr.find("7") != std::string::npos)
                    {
                        nPtBins    = nPtBins_Njet7_R;
                        ptBinEdges = ptBinEdges_Njet7_R; 
                    }
                    else if (njetStr.find("8") != std::string::npos)
                    {
                        nPtBins    = nPtBins_Njet8_R;
                        ptBinEdges = ptBinEdges_Njet8_R; 
                    }
                    else if (njetStr.find("9") != std::string::npos)
                    {
                        nPtBins    = nPtBins_Njet9_R;
                        ptBinEdges = ptBinEdges_Njet9_R; 
                    }
                    else
                    {
                        nPtBins    = nPtBins_Njet6_R;
                        ptBinEdges = ptBinEdges_Njet6_R; 
                    }
                }
                // Choose the binning for the merged case
                else if (tagStr.find("M") != std::string::npos)
                {
                    if (njetStr.find("6") != std::string::npos)
                    {
                        nPtBins    = nPtBins_Njet6_M;
                        ptBinEdges = ptBinEdges_Njet6_M; 
                    }
                    else if (njetStr.find("7") != std::string::npos)
                    {
                        nPtBins    = nPtBins_Njet7_M;
                        ptBinEdges = ptBinEdges_Njet7_M; 
                    }
                    else if (njetStr.find("8") != std::string::npos)
                    {
                        nPtBins    = nPtBins_Njet8_M;
                        ptBinEdges = ptBinEdges_Njet8_M; 
                    }
                    else if (njetStr.find("9") != std::string::npos)
                    {
                        nPtBins    = nPtBins_Njet9_M;
                        ptBinEdges = ptBinEdges_Njet9_M; 
                    }
                    else
                    {
                        nPtBins    = nPtBins_Njet6_M;
                        ptBinEdges = ptBinEdges_Njet6_M; 
                    }
                }

                for(const auto& h1Dinfo : hist1Dinfos)
                { 
                    std::string name = h1Dinfo.name + cutVar.first + tagStr + njetStr;

                    // Only make the 1D specially-binned histograms for specifically numerator and denominator
                    // and only for the fully formed CRs (not partial selections)
                    if      (name.find("Binned") != std::string::npos and tag != "" and name.find("CR") != std::string::npos)
                    {
                        if      (name.find("topPt") != std::string::npos)
                            my_1D_histos.emplace(name, std::make_shared<TH1D>(name.c_str(), name.c_str(), nPtBins,   vecToDouble_t(ptBinEdges)));
                        else if (name.find("nJets") != std::string::npos)
                            my_1D_histos.emplace(name, std::make_shared<TH1D>(name.c_str(), name.c_str(), nJetBins,  njetBinEdges));
                        else if (name.find("HT") != std::string::npos)
                            my_1D_histos.emplace(name, std::make_shared<TH1D>(name.c_str(), name.c_str(), nHTBins,   htBinEdges));
                        else if (name.find("Disc") != std::string::npos)
                            my_1D_histos.emplace(name, std::make_shared<TH1D>(name.c_str(), name.c_str(), nDiscBins, discBinEdges));
                    }
            
                    // For a 1D histo of the selection variables, make for any version except for specific numerator/denominator
                    else if (name.find("Binned") == std::string::npos and tag == "" and njet == "")
                        my_1D_histos.emplace(name, std::make_shared<TH1D>(name.c_str(), name.c_str(), h1Dinfo.nBins, h1Dinfo.low, h1Dinfo.high));
                }
                for(const auto& h2Dinfo : hist2Dinfos)
                { 
                    std::string name = h2Dinfo.name + cutVar.first + tagStr + njetStr;

                    // No 2D histos are made for any partial selection, just the full CR selections
                    if (name.find("CR") == std::string::npos)
                        continue;

                    if      (njet == "" and name.find("Binned") == std::string::npos)
                    {
                        if (name.find("HT") != std::string::npos)
                            my_2D_histos.emplace(name, std::make_shared<TH2D>(name.c_str(), name.c_str(), h2Dinfo.nBinsX, h2Dinfo.lowX, h2Dinfo.highX, h2Dinfo.nBinsY, h2Dinfo.lowY, h2Dinfo.highY));
                        else if (name.find("nJets") != std::string::npos)
                            my_2D_histos.emplace(name, std::make_shared<TH2D>(name.c_str(), name.c_str(), h2Dinfo.nBinsX, h2Dinfo.lowX, h2Dinfo.highX, h2Dinfo.nBinsY, h2Dinfo.lowY, h2Dinfo.highY));
                    }
                    else if (njet == "" and tag != "" and name.find("Binned") != std::string::npos)
                    {
                        if (name.find("nJets") != std::string::npos)
                            my_2D_histos.emplace(name, std::make_shared<TH2D>(name.c_str(), name.c_str(), nPtBins, vecToDouble_t(ptBinEdges), nJetBins,  njetBinEdges));
                    }
                }
            }
        }
    }
}

void TopTaggerSF_Analyzer::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    while (tr.getNextEvent())
    {
        const auto& eventCounter = tr.getVar<int>("eventCounter");

        if ( maxevents != -1 && tr.getEvtNum() > maxevents )   
            break;

        if ( tr.getEvtNum() % 1000 == 0)
            printf( " Event %i\n", tr.getEvtNum() );        

        const auto& runtype = tr.getVar<std::string>("runtype");

        // Top canddiate pt and disc
        const auto& bestTopPt     = tr.getVar<double>("bestTopPt");
        const auto& bestTopDisc   = tr.getVar<double>("bestTopDisc");
        const auto& bestTopNconst = tr.getVar<int>("bestTopNconst");

        const auto& resolvedTop_WP = tr.getVar<double>("resolvedTop_WP");
        const auto& mergedTop_WP   = tr.getVar<double>("mergedTop_WP");

        // Distinguish if the best top candidate is resolved or merged
        bool pass_TightResolvedTop = bestTopDisc > resolvedTop_WP and bestTopNconst == 3;
        bool pass_AnyResolvedTop   = bestTopDisc > 0.0            and bestTopNconst == 3;
       
        bool pass_TightMergedTop   = bestTopDisc > mergedTop_WP and bestTopNconst == 1;
        bool pass_AnyMergedTop     = bestTopDisc > 0.0          and bestTopNconst == 1;

        // Event-level booleans (triggers, filters, etc.)
        const auto& pass_METFilters  = tr.getVar<bool>("passMETFilters");
        const auto& pass_MadHT       = tr.getVar<bool>("passMadHT");
        const auto& pass_MuonTrigger = tr.getVar<bool>("passTriggerMuon");
        const auto& pass_QCDTrigger  = tr.getVar<bool>("passTriggerAllHad");
        const auto& pass_HEMveto     = tr.getVar<bool>("passElectronHEMveto");
        const auto& pass_JetID       = tr.getVar<bool>("JetID");        

        const auto& disc1            = tr.getVar<double>("DoubleDisCo_disc1_0l_RPV");
        const auto& disc2            = tr.getVar<double>("DoubleDisCo_disc2_0l_RPV");

        bool pass_DataQuality = pass_METFilters && pass_HEMveto && pass_JetID;

        // Lepton quantities
        const auto& Muons          = tr.getVec<utility::LorentzVector>("Muons");
        const auto& GoodMuons      = tr.getVec<bool>("GoodMuons");
        const auto& NGoodMuons     = tr.getVar<int>("NGoodMuons");
        const auto& NGoodElectrons = tr.getVar<int>("NGoodElectrons");

        // Jet quantities
        const auto& Jets                  = tr.getVec<utility::LorentzVector>("Jets");
        const auto& GoodJets_pt30         = tr.getVec<bool>("GoodJets_pt30");
        const auto& NGoodJets_pt30        = tr.getVar<int>("NGoodJets_pt30");
        const auto& NGoodBJets_pt30       = tr.getVar<int>("NGoodBJets_pt30");
        const auto& GoodBJets_pt30_loose  = tr.getVec<bool>("GoodBJets_pt30_loose");
        const auto& GoodBJets_pt30        = tr.getVec<bool>("GoodBJets_pt30");
        const auto& NGoodBJets_pt30_loose = tr.getVar<int>("NBJets_pt30_loose");

        // Event-level quantities
        const auto& HT_trigger_pt30 = tr.getVar<double>("HT_trigger_pt30");
        const auto& MET             = tr.getVar<float>("MET"); 
        const auto& METPhi          = tr.getVar<float>("METPhi");

        // Get the 4-vec for the MET
        utility::LorentzVector lvMET;
        lvMET.SetPt(MET); lvMET.SetEta(0.0); lvMET.SetPhi(METPhi); lvMET.SetE(MET);

        bool pass_ZeroMuons     = NGoodMuons     == 0;
        bool pass_ZeroElectrons = NGoodElectrons == 0;

        bool pass_SingleMuon    = NGoodMuons            == 1;
        bool pass_ExtraLepVeto  = NGoodElectrons        == 0;
        bool pass_ge6jets       = NGoodJets_pt30        >= 6;
        bool pass_ge9jets       = NGoodJets_pt30        >= 9;
        bool pass_1bjet         = NGoodBJets_pt30       >= 1;

        // The medium b jet (above) also passes the loose working point
        // So if we want an additional b jet that is loose we will have at least 2 loose b
        bool pass_1bjetLoose    = NGoodBJets_pt30_loose >= 2;
        bool pass_HT200         = HT_trigger_pt30        > 200.0;
        bool pass_HT500         = HT_trigger_pt30        > 500.0;
        bool pass_MET50         = MET                    > 50.0;

        bool pass_HT_QCDCR = pass_HT500;

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
        for (unsigned int iJet = 0; iJet < Jets.size(); iJet++)
        {
            if (GoodJets_pt30.at(iJet))
            {
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

            // Check for another b jet that is loose but not also just the medium one
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
    
                    // Save the dR and mass for the loose b + muon combination with smallest dR
                    if (tempMuonBjetdR < MuonBjetdR)
                    {
                        MuonBjetdR   = tempMuonBjetdR;
                        MuonBjetMass = tempMuonBjetMass;
                    }

                    // At least one loose b jet + muon combination to give inv mass inside window around top---the tag
                    pass_MuonBjetMass |= (tempMuonBjetMass > 30.0 && tempMuonBjetMass < 180.0);

                    // Require that the loose b jet satisfies both criteria at same time
                    // E.g. one loose b cannot satisfy dR while another loose b satisfies invariant mass 
                    pass_MuonBjetdRMass |= (pass_MuonBjetdR && pass_MuonBjetMass);
                }
            }
        }

        pass_JetMETdPhi = pass_Jet1METdPhi && pass_Jet2METdPhi && pass_Jet3METdPhi;

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
        if(runtype == "MC")
        {
            // Define Lumi weight
            const auto& Weight = tr.getVar<float>("Weight");
            const auto& lumi   = tr.getVar<double>("FinalLumi");
            double eventweight = lumi * Weight;

            double bTagScaleFactor        = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            double prefiringScaleFactor   = tr.getVar<double>("prefiringScaleFactor");
            double puScaleFactor          = tr.getVar<double>("puWeightCorr");
            double topPtScaleFactor       = tr.getVar<double>("topPtScaleFactor");

            weightTTbar *= eventweight * puScaleFactor * bTagScaleFactor * prefiringScaleFactor * topPtScaleFactor;
            weightQCD   *= eventweight * puScaleFactor * bTagScaleFactor * prefiringScaleFactor * topPtScaleFactor;
        }
    
        // To be used with Single-Muon-triggered data and measuring tagging efficiency SF
        bool pass_SemiLepTTbarCR = pass_DataQuality                                         &&
                                   pass_MadHT                                               &&
                                   pass_MuonTrigger                                         &&
                                   pass_SingleMuon                                          &&
                                   pass_ge6jets                                             &&
                                   pass_ExtraLepVeto                                        &&
                                   pass_1bjet                                               &&
                                   //pass_Jet1METdPhi && pass_Jet2METdPhi && pass_Jet3METdPhi &&
                                   pass_1bjetLoose                                          &&
                                   pass_MuonBjetdRMass                                      &&
                                   pass_MuonMETdPhi                                         &&
                                   pass_MuonMETmT                                           &&
                                   pass_HT200                                                ; 
                                   //pass_MET50                                               ;

        // For use with JetHT-triggered data and measuring mistag SF
        bool pass_QCDCR = pass_DataQuality   && 
                          pass_QCDTrigger    &&
                          pass_MadHT         &&
                          pass_ge6jets       &&
                          pass_ZeroMuons     &&
                          pass_ZeroElectrons &&
                          pass_HT_QCDCR       ;

        const std::map<std::string, bool> cut_map
        {
            {                                                                          "", pass_DataQuality && pass_MadHT},
            {                                                                      "_1mu", pass_DataQuality && pass_MadHT && pass_MuonTrigger && pass_SingleMuon},
            {                                                                 "_1mu_ge6j", pass_DataQuality && pass_MadHT && pass_MuonTrigger && pass_SingleMuon && pass_ge6jets},
            {                                                             "_1mu_ge6j_0el", pass_DataQuality && pass_MadHT && pass_MuonTrigger && pass_SingleMuon && pass_ge6jets && pass_ExtraLepVeto},
            {                                                          "_1mu_ge6j_0el_1b", pass_DataQuality && pass_MadHT && pass_MuonTrigger && pass_SingleMuon && pass_ge6jets && pass_ExtraLepVeto && pass_1bjet},
            {                                                 "_1mu_ge6j_0el_1b_jMETdPhi", pass_DataQuality && pass_MadHT && pass_MuonTrigger && pass_SingleMuon && pass_ge6jets && pass_ExtraLepVeto && pass_1bjet && pass_JetMETdPhi},
            {                                             "_1mu_ge6j_0el_1b_jMETdPhi_1bl", pass_DataQuality && pass_MadHT && pass_MuonTrigger && pass_SingleMuon && pass_ge6jets && pass_ExtraLepVeto && pass_1bjet && pass_JetMETdPhi && pass_1bjetLoose},
            {                                  "_1mu_ge6j_0el_1b_jMETdPhi_1bl_mubdR_mubM", pass_DataQuality && pass_MadHT && pass_MuonTrigger && pass_SingleMuon && pass_ge6jets && pass_ExtraLepVeto && pass_1bjet && pass_JetMETdPhi && pass_1bjetLoose && pass_MuonBjetdRMass},
            {                         "_1mu_ge6j_0el_1b_jMETdPhi_1bl_mubdR_mubM_muMETphi", pass_DataQuality && pass_MadHT && pass_MuonTrigger && pass_SingleMuon && pass_ge6jets && pass_ExtraLepVeto && pass_1bjet && pass_JetMETdPhi && pass_1bjetLoose && pass_MuonBjetdRMass && pass_MuonMETdPhi},
            {                 "_1mu_ge6j_0el_1b_jMETdPhi_1bl_mubdR_mubM_muMETphi_muMETmT", pass_DataQuality && pass_MadHT && pass_MuonTrigger && pass_SingleMuon && pass_ge6jets && pass_ExtraLepVeto && pass_1bjet && pass_JetMETdPhi && pass_1bjetLoose && pass_MuonBjetdRMass && pass_MuonMETdPhi && pass_MuonMETmT},
            {      "_1mu_ge6j_0el_1b_jMETdPhi_1bl_mubdR_mubM_muMETphi_muMETmT_HT200_pt30", pass_DataQuality && pass_MadHT && pass_MuonTrigger && pass_SingleMuon && pass_ge6jets && pass_ExtraLepVeto && pass_1bjet && pass_JetMETdPhi && pass_1bjetLoose && pass_MuonBjetdRMass && pass_MuonMETdPhi && pass_MuonMETmT && pass_HT200},
            {"_1mu_ge6j_0el_1b_jMETdPhi_1bl_mubdR_mubM_muMETphi_muMETmT_HT200_pt30_MET50", pass_DataQuality && pass_MadHT && pass_MuonTrigger && pass_SingleMuon && pass_ge6jets && pass_ExtraLepVeto && pass_1bjet && pass_JetMETdPhi && pass_1bjetLoose && pass_MuonBjetdRMass && pass_MuonMETdPhi && pass_MuonMETmT && pass_HT200 && pass_MET50},

            {"_QCDCR"                                                           , pass_QCDCR},
            {"_TTbarCR"                                                         , pass_SemiLepTTbarCR},
            {"_QCDCR_ge9j"                                                      , pass_QCDCR          && pass_ge9jets},
            {"_TTbarCR_ge9j"                                                    , pass_SemiLepTTbarCR && pass_ge9jets},

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
            {"",       true},
            {"6",      NGoodJets_pt30 == 6},
            {"7",      NGoodJets_pt30 == 7},
            {"8",      NGoodJets_pt30 == 8},
            {"9incl",  NGoodJets_pt30 >= 9},
        };

        std::map<std::string, bool> numDenMap =
        {
            {"",    true},
            {"numR", pass_AnyResolvedTop && pass_TightResolvedTop},
            {"denR", pass_AnyResolvedTop},
            {"numM", pass_AnyMergedTop && pass_TightMergedTop},
            {"denM", pass_AnyMergedTop},
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

                        std::string tagStr = "";
                        if (numDen.first != "")
                            tagStr = "_" + numDen.first;

                        std::string njetStr = "";
                        if (njetCut.first != "")
                            njetStr = "_Njets" + njetCut.first;

                        std::string selection = kv.first + tagStr + njetStr;

                        if (njetCut.first == "" and numDen.first == "")
                        {
                            my_1D_histos["nElectrons"   + selection]->Fill(NGoodElectrons,        weight);
                            my_1D_histos["nMuons"       + selection]->Fill(NGoodMuons,            weight);
                            my_1D_histos["nJets"        + selection]->Fill(NGoodJets_pt30,        weight);
                            my_1D_histos["nBJets"       + selection]->Fill(NGoodBJets_pt30,       weight);
                            my_1D_histos["nBJetsLoose"  + selection]->Fill(NGoodBJets_pt30_loose, weight);
                            my_1D_histos["ht"           + selection]->Fill(HT_trigger_pt30,       weight);
                            my_1D_histos["met"          + selection]->Fill(MET,                   weight);
                            my_1D_histos["jet1METdPhi"  + selection]->Fill(Jet1METdPhi,           weight);
                            my_1D_histos["jet2METdPhi"  + selection]->Fill(Jet2METdPhi,           weight);
                            my_1D_histos["jet3METdPhi"  + selection]->Fill(Jet3METdPhi,           weight);
                            my_1D_histos["muonBjetDR"   + selection]->Fill(MuonBjetdR,            weight);
                            my_1D_histos["muonBjetMass" + selection]->Fill(MuonBjetMass,          weight);
                            my_1D_histos["muonMETdPhi"  + selection]->Fill(MuonMETdPhi,           weight);
                            my_1D_histos["muonMETmT"    + selection]->Fill(MuonMETmT,             weight);
                            my_1D_histos["topPt"        + selection]->Fill(bestTopPt,             weight);
                        }

                        // Make versions of variables if they go into overflow bin
                        // Such that they make it into the previous bin
                        double correctTopPt = bestTopPt;
                        if (correctTopPt >= 1500.0)
                            correctTopPt = 1499.0;

                        double correctHT = HT_trigger_pt30;
                        if (correctHT >= 3000.0)
                            correctHT = 2999.0;

                        int correctNjets = NGoodJets_pt30;
                        if (correctNjets >= 10)
                            correctNjets = 10;

                        if (njetCut.first == "" and selection.find("CR") != std::string::npos)
                        {
                            my_2D_histos["topPt_vs_nJets" + selection]->Fill(correctTopPt, correctNjets, weight);
                            my_2D_histos["topPt_vs_HT"    + selection]->Fill(correctTopPt, correctHT,    weight);
                        }

                        if (numDen.first != "" and selection.find("CR") != std::string::npos)
                        {
                            my_1D_histos["Binned_topPt" + selection]->Fill(correctTopPt, weight);
                            my_1D_histos["Binned_nJets" + selection]->Fill(correctNjets, weight);
                            my_1D_histos["Binned_HT"    + selection]->Fill(correctHT,    weight);
                            my_1D_histos["Binned_Disc1" + selection]->Fill(disc1,        weight);
                            my_1D_histos["Binned_Disc2" + selection]->Fill(disc2,        weight);

                            if (njetCut.first == "")
                            {
                                my_2D_histos["Binned_topPt_vs_nJets" + selection]->Fill(correctTopPt, correctNjets, weight);
                            }
                        }
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
