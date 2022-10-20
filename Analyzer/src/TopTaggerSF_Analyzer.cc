#define TopTaggerSF_Analyzer_cxx
#include "Analyzer/Analyzer/include/TopTaggerSF_Analyzer.h"
#include "Framework/Framework/include/SetUpTopTagger.h"
#include "Framework/Framework/include/Utility.h"
#include "NTupleReader/include/NTupleReader.h"

#include <iostream>

TopTaggerSF_Analyzer::TopTaggerSF_Analyzer() : ttHistsNjetIncl("NjetIncl_tt"), ttHistsNjetIncl_Joe("NjetIncl_tt_joe"), ttHistsNjet7("Njet7_tt"), ttHistsNjet8("Njet8_tt"), ttHistsNjet9("Njet9_tt"), ttHistsNjet10("Njet10_tt"), ttHistsNjet11("Njet11_tt"), ttHistsNjet12incl("Njet12incl_tt"), qcdHistsNjetIncl("NjetIncl_qcd"), qcdHistsNjetIncl_Joe("NjetIncl_qcd_joe"), qcdHistsNjet7("Njet7_qcd"), qcdHistsNjet8("Njet8_qcd"), qcdHistsNjet9("Njet9_qcd"), qcdHistsNjet10("Njet10_qcd"), qcdHistsNjet11("Njet11_qcd"), qcdHistsNjet12incl("Njet12incl_qcd")
{

    initHistos = false;

    histInfos = {
        {"h_nElectrons",           6,   -0.5,    5.5},
        {"h_nMuons",               6,   -0.5,    5.5},
        {"h_nJets",                21,  -0.5,   20.5},
        {"h_nBJets",               21,  -0.5,   20.5},
        {"h_nBJetsLoose",          21,  -0.5,   20.5},
        {"h_ht",                  500,     0,   5000},
        {"h_met",                 200,     0,   1500},
        {"h_jet1METdPhi",          80,  -4.0,    4.0},
        {"h_jet2METdPhi",          80,  -4.0,    4.0},
        {"h_jet3METdPhi",          80,  -4.0,    4.0},
        {"h_muonBjetDR",          200,   0.0,   10.0},
        {"h_muonBjetMass",        500,   0.0, 1000.0},
        {"h_muonMETdPhi",          80,  -4.0,    4.0},
        {"h_muonMETmT",          1000,   0.0, 1000.0},
    };
}

void TopTaggerSF_Analyzer::InitHistos(const std::map<std::string, bool>& cutMap)
{
    TH1::SetDefaultSumw2();

    my_histos.emplace("EventCounter", std::make_shared<TH1D>("EventCounter","EventCounter", 2, -1.1, 1.1 ) );

    for(auto& mycut : cutMap)
    {
        for(const auto& hInfo : histInfos)
        { 
            my_histos.emplace(hInfo.name+mycut.first, std::make_shared<TH1D>((hInfo.name+mycut.first).c_str(),(hInfo.name+mycut.first).c_str(), hInfo.nBins, hInfo.low, hInfo.high));
        }
    }
}

void TopTaggerSF_Analyzer::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    TRandom3 rand(123);
   
    while (tr.getNextEvent())
    {
        const auto& eventCounter           = tr.getVar<int>("eventCounter");

        if ( maxevents != -1 && tr.getEvtNum() >= maxevents )   
            break;

        if ( tr.getEvtNum() % 1000 == 0)
            printf( " Event %i\n", tr.getEvtNum() );        

        const auto& runtype          = tr.getVar<std::string>("runtype");

        // Event booleans
        const auto& passMETFilters   = tr.getVar<bool>("passMETFilters");
        const auto& passMadHT        = tr.getVar<bool>("passMadHT");

        // Lepton quantities
        const auto& Muons            = tr.getVec<utility::LorentzVector>("Muons");
        const auto& GoodMuons        = tr.getVec<bool>("GoodMuons");
        const auto& NGoodMuons       = tr.getVar<int>("NGoodMuons");
        const auto& NGoodElectrons   = tr.getVar<int>("NGoodElectrons");

        // Jet quantities
        const auto& Jets                  = tr.getVec<utility::LorentzVector>("Jets");
        const auto& JetID                 = tr.getVar<bool>("JetID");        
        const auto& GoodJets_pt30         = tr.getVec<bool>("GoodJets_pt30");
        const auto& NGoodJets_pt30        = tr.getVar<int>("NGoodJets_pt30");
        const auto& NGoodBJets_pt30       = tr.getVar<int>("NGoodBJets_pt30");
        const auto& GoodBJets_pt30_loose  = tr.getVec<bool>("GoodBJets_pt30_loose");
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
        bool pass_1bjetLoose    = NGoodBJets_pt30_loose >= 1;
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

            if (GoodBJets_pt30_loose.at(iJet))
            {
                for (unsigned int iMuon = 0; iMuon < Muons.size(); iMuon++)
                {
                    if (!GoodMuons.at(iMuon))
                        continue;

                    double tempMuonBjetdR = utility::DeltaR(Jets.at(iJet), Muons.at(iMuon));

                    // At least one loose Bjet be within a dR of 1.5
                    pass_MuonBjetdR |= (tempMuonBjetdR < 1.5);
    
                    auto combinedMuonBjet = Jets.at(iJet) + Muons.at(iMuon);
                    double tempMuonBjetMass = combinedMuonBjet.M();
    
                    if (tempMuonBjetdR < MuonBjetdR)
                    {
                        MuonBjetdR = tempMuonBjetdR;
                        MuonBjetMass = tempMuonBjetMass;
                    }

                    // At least one loose Bjet + muon combination to give inv mass inside window around top---the tag
                    pass_MuonBjetMass |= (tempMuonBjetMass >= 30.0 && tempMuonBjetMass <= 180.0);

                    // Require that the loose b jet satisfies both criteria at same time
                    // E.g. one loose b cannot satisfy dR while another loose b satisfies invariant mass 
                    pass_MuonBjetdRMass |= ((tempMuonBjetdR < 1.5) && (tempMuonBjetMass >= 30.0 && tempMuonBjetMass <= 180.0));
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

            auto combinedMuonMET = Muons.at(iMuon) + lvMET;

            double tempMuonMETdPhi  = utility::DeltaPhi(Muons.at(iMuon), lvMET);
            double tempMuonMETmT    = utility::calcMT(Muons.at(iMuon), lvMET);

            if (fabs(tempMuonMETdPhi) < fabs(MuonMETdPhi))
                MuonMETdPhi = tempMuonMETdPhi;
            if (fabs(tempMuonMETmT) < fabs(MuonMETmT))
                MuonMETmT = tempMuonMETmT;

            pass_MuonMETdPhi |= (fabs(tempMuonMETdPhi) < 0.8);
            pass_MuonMETmT   |= (tempMuonMETmT < 100.0);
        }

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
            // Define Lumi weight
            const auto& Weight = tr.getVar<float>("Weight");
            const auto& lumi   = tr.getVar<double>("FinalLumi");
            eventweight        = lumi*Weight;

            bTagScaleFactor      = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor");
            puScaleFactor        = tr.getVar<double>("puWeightCorr");

            weight *= eventweight*puScaleFactor;
        }

        // To be used with Single-Muon-triggered data and measuring tagging efficiency SF
        bool pass_SemiLepTTbarCR = JetID                                                    &&
                                   passMETFilters                                           &&
                                   passMadHT                                                &&
                                   pass_SingleMuon                                          &&
                                   pass_7jets                                               &&
                                   pass_ExtraLepVeto                                        &&
                                   pass_1bjet                                               &&
                                   //pass_Jet1METdPhi && pass_Jet2METdPhi && pass_Jet3METdPhi &&
                                   pass_1bjetLoose                                          &&
                                   pass_MuonBjetdRMass                                      &&
                                   pass_MuonMETdPhi                                         &&
                                   pass_MuonMETmT                                           &&
                                   pass_HT500                                               &&
                                   pass_MET50                                               ;

        bool pass_SemiLepTTbarCR_Joe = JetID                                                    &&
                                       passMETFilters                                           &&
                                       passMadHT                                                &&
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
            {                                                                  "", true},
            {                                                              "_1mu", pass_SingleMuon},
            {                                                           "_1mu_7j", pass_SingleMuon && pass_7jets},
            {                                                       "_1mu_7j_0el", pass_SingleMuon && pass_7jets && pass_ExtraLepVeto},
            {                                                    "_1mu_7j_0el_1b", pass_SingleMuon && pass_7jets && pass_ExtraLepVeto && pass_1bjet},
            {                                            "_1mu_7j_0el_1b",         pass_SingleMuon && pass_7jets && pass_ExtraLepVeto && pass_1bjet /*&& pass_JetMETdPhi*/},
            {                                        "_1mu_7j_0el_1b_1bl",         pass_SingleMuon && pass_7jets && pass_ExtraLepVeto && pass_1bjet /*&& pass_JetMETdPhi*/ && pass_1bjetLoose},
            {                             "_1mu_7j_0el_1b_1bl_mubdR_mubM",         pass_SingleMuon && pass_7jets && pass_ExtraLepVeto && pass_1bjet /*&& pass_JetMETdPhi*/ && pass_1bjetLoose && pass_MuonBjetdRMass},
            {                    "_1mu_7j_0el_1b_1bl_mubdR_mubM_muMETphi",         pass_SingleMuon && pass_7jets && pass_ExtraLepVeto && pass_1bjet /*&& pass_JetMETdPhi*/ && pass_1bjetLoose && pass_MuonBjetdRMass && pass_MuonMETdPhi},
            {            "_1mu_7j_0el_1b_1bl_mubdR_mubM_muMETphi_muMETmT",         pass_SingleMuon && pass_7jets && pass_ExtraLepVeto && pass_1bjet /*&& pass_JetMETdPhi*/ && pass_1bjetLoose && pass_MuonBjetdRMass && pass_MuonMETdPhi && pass_MuonMETmT},
            {      "_1mu_7j_0el_1b_1bl_mubdR_mubM_muMETphi_muMETmT_HT500",         pass_SingleMuon && pass_7jets && pass_ExtraLepVeto && pass_1bjet /*&& pass_JetMETdPhi*/ && pass_1bjetLoose && pass_MuonBjetdRMass && pass_MuonMETdPhi && pass_MuonMETmT && pass_HT500},
            {"_1mu_7j_0el_1b_1bl_mubdR_mubM_muMETphi_muMETmT_HT500_MET50",         pass_SingleMuon && pass_7jets && pass_ExtraLepVeto && pass_1bjet /*&& pass_JetMETdPhi*/ && pass_1bjetLoose && pass_MuonBjetdRMass && pass_MuonMETdPhi && pass_MuonMETmT && pass_HT500 && pass_MET50},
        };

        // Initialize Histograms
        if(!initHistos)
        {
            InitHistos(cut_map);
            initHistos = true;
        }

        my_histos["EventCounter"]->Fill(eventCounter);

        for(auto& kv : cut_map)
        {
            if(kv.second)
            {
                my_histos["h_nElectrons"   + kv.first]->Fill(NGoodElectrons,        weight);
                my_histos["h_nMuons"       + kv.first]->Fill(NGoodMuons,            weight);
                my_histos["h_nJets"        + kv.first]->Fill(NGoodJets_pt30,        weight);
                my_histos["h_nBJets"       + kv.first]->Fill(NGoodBJets_pt30,       weight);
                my_histos["h_nBJetsLoose"  + kv.first]->Fill(NGoodBJets_pt30_loose, weight);
                my_histos["h_ht"           + kv.first]->Fill(HT_trigger_pt30,       weight);
                my_histos["h_met"          + kv.first]->Fill(MET,                   weight);
                my_histos["h_jet1METdPhi"  + kv.first]->Fill(Jet1METdPhi,           weight);
                my_histos["h_jet2METdPhi"  + kv.first]->Fill(Jet2METdPhi,           weight);
                my_histos["h_jet3METdPhi"  + kv.first]->Fill(Jet3METdPhi,           weight);
                my_histos["h_muonBjetDR"   + kv.first]->Fill(MuonBjetdR,            weight);
                my_histos["h_muonBjetMass" + kv.first]->Fill(MuonBjetMass,          weight);
                my_histos["h_muonMETdPhi"  + kv.first]->Fill(MuonMETdPhi,           weight);
                my_histos["h_muonMETmT"    + kv.first]->Fill(MuonMETmT,             weight);
            }
        }

        // -------------------------------------------------
        // -- Fill the histograms for the ttbar CR selection
        // -------------------------------------------------
        std::vector<std::pair<std::string, bool>> ttbarCuts =
        {   
            {"pass_SemiLepTTbarCR", pass_SemiLepTTbarCR},
        };
        ttHistsNjetIncl.fillWithCutFlow(ttbarCuts, tr, weight, &rand);

        std::vector<std::pair<std::string, bool>> ttbarCuts_Joe =
        {   
            {"pass_SemiLepTTbarCR_Joe", pass_SemiLepTTbarCR_Joe},
        };
        ttHistsNjetIncl_Joe.fillWithCutFlow(ttbarCuts_Joe, tr, weight, &rand);

        // baseline cuts + ttNjets == 7
        std::vector<std::pair<std::string, bool>> ttNjets7 =
        {
            {"pass_SemiLepTTbarCR", pass_SemiLepTTbarCR},
            {"Njet7"              , NGoodJets_pt30 == 7},
        };
        ttHistsNjet7.fillWithCutFlow(ttNjets7, tr, weight, &rand);

        // baseline cuts + ttNjets == 8
        std::vector<std::pair<std::string, bool>> ttNjets8 =
        {
            {"pass_SemiLepTTbarCR", pass_SemiLepTTbarCR},
            {"Njet8"              , NGoodJets_pt30 == 8},
        };
        ttHistsNjet8.fillWithCutFlow(ttNjets8, tr, weight, &rand);

        // baseline cuts + ttNjets == 9
        std::vector<std::pair<std::string, bool>> ttNjets9 =
        {
            {"pass_SemiLepTTbarCR", pass_SemiLepTTbarCR},
            {"Njet9"              , NGoodJets_pt30 == 9},
        };
        ttHistsNjet9.fillWithCutFlow(ttNjets9, tr, weight, &rand);        
        
        // baseline cuts + ttNjets == 10
        std::vector<std::pair<std::string, bool>> ttNjets10 =
        {   
            {"pass_SemiLepTTbarCR", pass_SemiLepTTbarCR},
            {"Njet10"             , NGoodJets_pt30 == 10},
        };
        ttHistsNjet10.fillWithCutFlow(ttNjets10, tr, weight, &rand);

        // baseline cuts + ttNjets == 11
        std::vector<std::pair<std::string, bool>> ttNjets11 =
        {
            {"pass_SemiLepTTbarCR", pass_SemiLepTTbarCR},
            {"Njet11"             , NGoodJets_pt30 == 11},
        };
        ttHistsNjet11.fillWithCutFlow(ttNjets11, tr, weight, &rand);

        // baseline cuts + ttNjets == 12
        std::vector<std::pair<std::string, bool>> ttNjets12incl =
        {
            {"pass_SemiLepTTbarCR", pass_SemiLepTTbarCR},
            {"Njet12incl"         , NGoodJets_pt30 >= 12},
        };
        ttHistsNjet12incl.fillWithCutFlow(ttNjets12incl, tr, weight, &rand);


        // -------------------------------------------------
        // -- Fill the histograms for the QCD CR selection
        // -------------------------------------------------
        std::vector<std::pair<std::string, bool>> qcdCuts_Joe =
        {   
            {"pass_QCDCR_Joe", pass_QCDCR_Joe},
        };
        qcdHistsNjetIncl_Joe.fillWithCutFlow(qcdCuts_Joe, tr, weight, &rand);

        std::vector<std::pair<std::string, bool>> qcdCuts =
        {   
            {"pass_QCDCR", pass_QCDCR},
        };
        qcdHistsNjetIncl.fillWithCutFlow(qcdCuts, tr, weight, &rand);

        // baseline cuts + Njets == 7
        std::vector<std::pair<std::string, bool>> qcdNjets7 =
        {
            {"pass_QCDCR", pass_QCDCR},
            {"Njet7"     , NGoodJets_pt30 == 7},
        };
        qcdHistsNjet7.fillWithCutFlow(qcdNjets7, tr, weight, &rand);

        // baseline cuts + qcdNjets == 8
        std::vector<std::pair<std::string, bool>> qcdNjets8 =
        {
            {"pass_QCDCR", pass_QCDCR},
            {"Njet8"     , NGoodJets_pt30 == 8},
        };
        qcdHistsNjet8.fillWithCutFlow(qcdNjets8, tr, weight, &rand);

        // baseline cuts + qcdNjets == 9
        std::vector<std::pair<std::string, bool>> qcdNjets9 =
        {
            {"pass_QCDCR", pass_QCDCR},
            {"Njet9"     , NGoodJets_pt30 == 9},
        };
        qcdHistsNjet9.fillWithCutFlow(qcdNjets9, tr, weight, &rand);        
        
        // baseline cuts + qcdNjets == 10
        std::vector<std::pair<std::string, bool>> qcdNjets10 =
        {   
            {"pass_QCDCR", pass_QCDCR},
            {"Njet10"    , NGoodJets_pt30 == 10},
        };
        qcdHistsNjet10.fillWithCutFlow(qcdNjets10, tr, weight, &rand);

        // baseline cuts + qcdNjets == 11
        std::vector<std::pair<std::string, bool>> qcdNjets11 =
        {
            {"pass_QCDCR", pass_QCDCR},
            {"Njet11"    , NGoodJets_pt30 == 11},
        };
        qcdHistsNjet11.fillWithCutFlow(qcdNjets11, tr, weight, &rand);

        // baseline cuts + qcdNjets == 12
        std::vector<std::pair<std::string, bool>> qcdNjets12incl =
        {
            {"pass_QCDCR", pass_QCDCR},
            {"Njet12incl", NGoodJets_pt30 >= 12},
        };
        qcdHistsNjet12incl.fillWithCutFlow(qcdNjets12incl, tr, weight, &rand);
    }
}

void TopTaggerSF_Analyzer::WriteHistos(TFile* outfile)
{
    outfile->cd();

    for (const auto &p : my_histos) {
        p.second->Write();
    }
    
    ttHistsNjetIncl_Joe.save(outfile);
    ttHistsNjetIncl.save(outfile);
    ttHistsNjet7.save(outfile);
    ttHistsNjet8.save(outfile);
    ttHistsNjet9.save(outfile);
    ttHistsNjet10.save(outfile);
    ttHistsNjet11.save(outfile);
    ttHistsNjet12incl.save(outfile);

    qcdHistsNjetIncl_Joe.save(outfile);
    qcdHistsNjetIncl.save(outfile);
    qcdHistsNjet7.save(outfile);
    qcdHistsNjet8.save(outfile);
    qcdHistsNjet9.save(outfile);
    qcdHistsNjet10.save(outfile);
    qcdHistsNjet11.save(outfile);
    qcdHistsNjet12incl.save(outfile);

    outfile->Close();
}
