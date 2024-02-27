#ifndef Confg_h
#define Confg_h

#include "NTupleReader/include/NTupleReader.h"

#include "Framework/Framework/include/PrepNTupleVars.h"
#include "Framework/Framework/include/RunTopTagger.h"
#include "Framework/Framework/include/Muon.h"
#include "Framework/Framework/include/Electron.h"
#include "Framework/Framework/include/Photon.h"
#include "Framework/Framework/include/Jet.h"
#include "Framework/Framework/include/BJet.h"
#include "Framework/Framework/include/CommonVariables.h"
#include "Framework/Framework/include/MakeMVAVariables.h"
#include "Framework/Framework/include/Baseline.h"
#include "Framework/Framework/include/DeepEventShape.h"
#include "Framework/Framework/include/BTagCorrector.h"
#include "Framework/Framework/include/ScaleFactors.h"
#include "Framework/Framework/include/StopGenMatch.h"
#include "Framework/Framework/include/MegaJetCombine.h"
#include "Framework/Framework/include/MakeStopHemispheres.h"
#include "Framework/Framework/include/StopJets.h"
#include "Framework/Framework/include/ISRJets.h"

class Config
{
private:
    void registerModules(NTupleReader& tr, const std::vector<std::string>&& modules) const
    {
        const auto& runtype                           = tr.getVar<std::string>("runtype"                          );
        const auto& runYear                           = tr.getVar<std::string>("runYear"                          );
        const auto& DoubleDisCo_Cfg_0l_RPV            = tr.getVar<std::string>("DoubleDisCo_Cfg_0l_RPV"           );
        const auto& DoubleDisCo_Model_0l_RPV          = tr.getVar<std::string>("DoubleDisCo_Model_0l_RPV"         );
        const auto& DoubleDisCo_Cfg_NonIsoMuon_0l_RPV = tr.getVar<std::string>("DoubleDisCo_Cfg_NonIsoMuon_0l_RPV");
        const auto& DoubleDisCo_Cfg_1l_RPV            = tr.getVar<std::string>("DoubleDisCo_Cfg_1l_RPV"           );
        const auto& DoubleDisCo_Model_1l_RPV          = tr.getVar<std::string>("DoubleDisCo_Model_1l_RPV"         );
        const auto& DoubleDisCo_Cfg_NonIsoMuon_1l_RPV = tr.getVar<std::string>("DoubleDisCo_Cfg_NonIsoMuon_1l_RPV");
        const auto& DoubleDisCo_Cfg_0l_SYY            = tr.getVar<std::string>("DoubleDisCo_Cfg_0l_SYY"           );
        const auto& DoubleDisCo_Model_0l_SYY          = tr.getVar<std::string>("DoubleDisCo_Model_0l_SYY"         );
        const auto& DoubleDisCo_Cfg_NonIsoMuon_0l_SYY = tr.getVar<std::string>("DoubleDisCo_Cfg_NonIsoMuon_0l_SYY");
        const auto& DoubleDisCo_Cfg_1l_SYY            = tr.getVar<std::string>("DoubleDisCo_Cfg_1l_SYY"           );
        const auto& DoubleDisCo_Model_1l_SYY          = tr.getVar<std::string>("DoubleDisCo_Model_1l_SYY"         );
        const auto& DoubleDisCo_Cfg_NonIsoMuon_1l_SYY = tr.getVar<std::string>("DoubleDisCo_Cfg_NonIsoMuon_1l_SYY");
        const auto& DoubleDisCo_Cfg_2l_RPV            = tr.getVar<std::string>("DoubleDisCo_Cfg_2l_RPV"           );
        const auto& DoubleDisCo_Model_2l_RPV          = tr.getVar<std::string>("DoubleDisCo_Model_2l_RPV"         );
        const auto& DoubleDisCo_Cfg_NonIsoMuon_2l_RPV = tr.getVar<std::string>("DoubleDisCo_Cfg_NonIsoMuon_2l_RPV");
        const auto& DoubleDisCo_Cfg_2l_SYY            = tr.getVar<std::string>("DoubleDisCo_Cfg_2l_SYY"           );
        const auto& DoubleDisCo_Model_2l_SYY          = tr.getVar<std::string>("DoubleDisCo_Model_2l_SYY"         );
        const auto& DoubleDisCo_Cfg_NonIsoMuon_2l_SYY = tr.getVar<std::string>("DoubleDisCo_Cfg_NonIsoMuon_2l_SYY");
        const auto& leptonFileName                    = tr.getVar<std::string>("leptonFileName"                    );
        const auto& hadronicFileName                  = tr.getVar<std::string>("hadronicFileName"                  );
        const auto& toptaggerFileName                 = tr.getVar<std::string>("toptaggerFileName"                 );
        const auto& btagEffFileName                   = tr.getVar<std::string>("btagEffFileName"                   );
        const auto& bjetTagFileName                   = tr.getVar<std::string>("bjetTagFileName"                   );
        //const auto& bjetTagFileNameReshape            = tr.getVar<std::string>("bjetTagFileNameReshape"            );
        const auto& filetag                           = tr.getVar<std::string>("filetag"                           );
        const auto& meanFileName                      = tr.getVar<std::string>("meanFileName"                      );
        const auto& TopTaggerCfg                      = tr.getVar<std::string>("TopTaggerCfg"                      );

        for(const auto& module : modules)
        {
            if     (module=="PrepNTupleVars")                        tr.emplaceModule<PrepNTupleVars>();
            else if(module=="RunTopTagger")                          tr.emplaceModule<RunTopTagger>(TopTaggerCfg);
            else if(module=="Muon")                                  tr.emplaceModule<Muon>();
            else if(module=="Electron")                              tr.emplaceModule<Electron>();
            else if(module=="Photon")                                tr.emplaceModule<Photon>();
            else if(module=="Jet")                                   tr.emplaceModule<Jet>();
            else if(module=="JetAK8")                                tr.emplaceModule<JetAK8>();
            else if(module=="BJet")                                  tr.emplaceModule<BJet>();
            else if(module=="CommonVariables")                       tr.emplaceModule<CommonVariables>();
            else if(module=="MakeMVAVariables")                      tr.emplaceModule<MakeMVAVariables>(false, "", "GoodJets_pt30",        false, true, 7, 2, "");
            else if(module=="MakeMVAVariables_NonIsoMuon")           tr.emplaceModule<MakeMVAVariables>(false, "", "NonIsoMuonJets_pt30",  false, true, 7, 2, "");
            else if(module=="MakeMVAVariables_0l")                   tr.emplaceModule<MakeMVAVariables>(false, "", "GoodJets_pt30",        false, true, 7, 0, "_0l");
            else if(module=="MakeMVAVariables_NonIsoMuon_0l")        tr.emplaceModule<MakeMVAVariables>(false, "", "GoodJets_pt30",        false, true, 7, 0, "_0l");
            else if(module=="MakeMVAVariables_1l")                   tr.emplaceModule<MakeMVAVariables>(false, "", "GoodJets_pt30",        false, true, 7, 1, "_1l");
            else if(module=="MakeMVAVariables_NonIsoMuon_1l")        tr.emplaceModule<MakeMVAVariables>(false, "", "NonIsoMuonJets_pt30",  false, true, 7, 1, "_1l");
            else if(module=="MakeMVAVariables_2l")                   tr.emplaceModule<MakeMVAVariables>(false, "", "GoodJets_pt30",        false, true, 6, 2, "_2l");
            else if(module=="MakeMVAVariables_NonIsoMuon_2l")        tr.emplaceModule<MakeMVAVariables>(false, "", "NonIsoMuonJets_pt30",  false, true, 7, 2, "_2l");
            else if(module=="Baseline")                              tr.emplaceModule<Baseline>();
            else if(module=="StopGenMatch")                          tr.emplaceModule<StopGenMatch>();
            else if(module=="MegaJetCombine")                        tr.emplaceModule<MegaJetCombine>();
            else if(module=="MakeStopHemispheres_All")               tr.emplaceModule<MakeStopHemispheres>("Jets",     "AllJets",                 "NJets",                    "_All",                "", Hemisphere::InvMassSeed);
            else if(module=="MakeStopHemispheres_OldSeed")           tr.emplaceModule<MakeStopHemispheres>("Jets",     "GoodJets_pt20",           "NGoodJets_pt20",           "_OldSeed",            "", Hemisphere::InvMassSeed);
            else if(module=="MakeStopHemispheres_OldSeed_maskedISR") tr.emplaceModule<MakeStopHemispheres>("Jets",     "GoodJets_pt20_maskedISR", "NGoodJets_pt20_maskedISR", "_OldSeed_maskedISR",  "", Hemisphere::InvMassSeed);
            else if(module=="MakeStopHemispheres_OldSeed_NonIsoMuon")tr.emplaceModule<MakeStopHemispheres>("Jets",     "NonIsoMuonJets_pt20",     "NNonIsoMuonJets_pt30",     "_OldSeed_NonIsoMuon", "", Hemisphere::InvMassSeed);
            else if(module=="MakeStopHemispheres_TopSeed")           tr.emplaceModule<MakeStopHemispheres>("StopJets", "GoodStopJets",            "NGoodStopJets",            "_TopSeed",            "", Hemisphere::TopSeed);
            else if(module=="MakeStopHemispheres_TopSeed_maskedISR") tr.emplaceModule<MakeStopHemispheres>("StopJets", "GoodStopJets_maskedISR",  "NGoodStopJets_maskedISR",  "_TopSeed_maskedISR",  "", Hemisphere::TopSeed);
            else if(module=="MakeStopHemispheres_TopSeed_NonIsoMuon")tr.emplaceModule<MakeStopHemispheres>("StopJets", "GoodStopJets",            "NGoodStopJets",            "_TopSeed_NonIsoMuon", "", Hemisphere::TopSeed);
            else if(module=="StopJets")                              tr.emplaceModule<StopJets>();
            else if(module=="ISRJets")                               tr.emplaceModule<ISRJets>();
            else if(module=="DoubleDisCo_0l_RPV")                    tr.emplaceModule<DeepEventShape>(DoubleDisCo_Cfg_0l_RPV, DoubleDisCo_Model_0l_RPV);
            else if(module=="DoubleDisCo_NonIsoMuon_0l_RPV")         tr.emplaceModule<DeepEventShape>(DoubleDisCo_Cfg_NonIsoMuon_0l_RPV, DoubleDisCo_Model_0l_RPV);
            else if(module=="DoubleDisCo_1l_RPV")                    tr.emplaceModule<DeepEventShape>(DoubleDisCo_Cfg_1l_RPV, DoubleDisCo_Model_1l_RPV);
            else if(module=="DoubleDisCo_NonIsoMuon_1l_RPV")         tr.emplaceModule<DeepEventShape>(DoubleDisCo_Cfg_NonIsoMuon_1l_RPV, DoubleDisCo_Model_1l_RPV);
            else if(module=="DoubleDisCo_2l_RPV")                    tr.emplaceModule<DeepEventShape>(DoubleDisCo_Cfg_2l_RPV, DoubleDisCo_Model_2l_RPV);
            else if(module=="DoubleDisCo_NonIsoMuon_2l_RPV")         tr.emplaceModule<DeepEventShape>(DoubleDisCo_Cfg_NonIsoMuon_2l_RPV, DoubleDisCo_Model_2l_RPV);
            else if(module=="DoubleDisCo_0l_SYY")                    tr.emplaceModule<DeepEventShape>(DoubleDisCo_Cfg_0l_SYY, DoubleDisCo_Model_0l_SYY);
            else if(module=="DoubleDisCo_NonIsoMuon_0l_SYY")         tr.emplaceModule<DeepEventShape>(DoubleDisCo_Cfg_NonIsoMuon_0l_SYY, DoubleDisCo_Model_0l_SYY);
            else if(module=="DoubleDisCo_1l_SYY")                    tr.emplaceModule<DeepEventShape>(DoubleDisCo_Cfg_1l_SYY, DoubleDisCo_Model_1l_SYY);
            else if(module=="DoubleDisCo_NonIsoMuon_1l_SYY")         tr.emplaceModule<DeepEventShape>(DoubleDisCo_Cfg_NonIsoMuon_1l_SYY, DoubleDisCo_Model_1l_SYY);
            else if(module=="DoubleDisCo_2l_SYY")                    tr.emplaceModule<DeepEventShape>(DoubleDisCo_Cfg_2l_SYY, DoubleDisCo_Model_2l_SYY);
            else if(module=="DoubleDisCo_NonIsoMuon_2l_SYY")         tr.emplaceModule<DeepEventShape>(DoubleDisCo_Cfg_NonIsoMuon_2l_SYY, DoubleDisCo_Model_2l_SYY);

            if(runtype == "MC")
            {
                if     (module=="ScaleFactors")  tr.emplaceModule<ScaleFactors>(runYear, leptonFileName, hadronicFileName, toptaggerFileName, meanFileName, filetag);
                else if(module=="BTagCorrector")
                {
                    std::string bjetTagFileNameReshape = "";
                    auto& bTagCorrector = tr.emplaceModule<BTagCorrector>(btagEffFileName, "", bjetTagFileName, bjetTagFileNameReshape, filetag);
                    bTagCorrector.SetVarNames("GenParticles_PdgId", "Jets", "GoodJets_pt30", "Jets_bJetTagDeepFlavourtotb", "Jets_partonFlavor");
                }
            }
        }
    }

public:
    Config()
    {
    }

    void setUp(NTupleReader& tr) const
    {
        //Get and make needed info
        const auto& filetag  = tr.getVar<std::string>("filetag");
        const auto& analyzer = tr.getVar<std::string>("analyzer");
        const bool isSignal  = (filetag.find("_stop") != std::string::npos || filetag.find("_mStop") != std::string::npos || filetag.find("VLQ_2t4b") != std::string::npos) ? true : false;

        std::string runYear;
        std::string bjetTagFileName, bjetTagFileNameReshape, TopTaggerCfg;
 
        double Lumi=0.0, Lumi_postHEM=-1.0, Lumi_preHEM=-1.0;
        double deepCSV_WP_loose=0.0, deepCSV_WP_medium=0.0, deepCSV_WP_tight=0.0, deepFlavour_WP_loose=0.0, deepFlavour_WP_medium=0.0, deepFlavour_WP_tight=0.0;
        double resolvedTop_WP = 0.0, mergedTop_WP = 0.0;
        bool blind = true;

        std::string DoubleDisCo_Cfg_0l_RPV            = "Keras_Tensorflow_DoubleDisCo_Reg_0l_RPV_Run2_MaxSig.cfg";           
        std::string DoubleDisCo_Model_0l_RPV          = "keras_frozen_DoubleDisCo_Reg_0l_RPV_Run2_MaxSig.pb";
        std::string DoubleDisCo_Cfg_NonIsoMuon_0l_RPV = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_0l_RPV_Run2_MaxSig.cfg"; 
        std::string DoubleDisCo_Cfg_1l_RPV            = "Keras_Tensorflow_DoubleDisCo_Reg_1l_RPV_Run2_MaxSig.cfg";
        std::string DoubleDisCo_Model_1l_RPV          = "keras_frozen_DoubleDisCo_Reg_1l_RPV_Run2_MaxSig.pb";
        std::string DoubleDisCo_Cfg_NonIsoMuon_1l_RPV = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_1l_RPV_Run2_MaxSig.cfg";
        std::string DoubleDisCo_Cfg_2l_RPV            = "Keras_Tensorflow_DoubleDisCo_Reg_2l_RPV_Run2_MaxSig.cfg";
        std::string DoubleDisCo_Model_2l_RPV          = "keras_frozen_DoubleDisCo_Reg_2l_RPV_Run2_MaxSig.pb";
        std::string DoubleDisCo_Cfg_NonIsoMuon_2l_RPV = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_2l_RPV_Run2_MaxSig.cfg";
        std::string DoubleDisCo_Cfg_0l_SYY            = "Keras_Tensorflow_DoubleDisCo_Reg_0l_SYY_Run2_MaxSig.cfg";           
        std::string DoubleDisCo_Model_0l_SYY          = "keras_frozen_DoubleDisCo_Reg_0l_SYY_Run2_MaxSig.pb";
        std::string DoubleDisCo_Cfg_NonIsoMuon_0l_SYY = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_0l_SYY_Run2_MaxSig.cfg";
        std::string DoubleDisCo_Cfg_1l_SYY            = "Keras_Tensorflow_DoubleDisCo_Reg_1l_SYY_Run2_MaxSig.cfg";
        std::string DoubleDisCo_Model_1l_SYY          = "keras_frozen_DoubleDisCo_Reg_1l_SYY_Run2_MaxSig.pb";
        std::string DoubleDisCo_Cfg_NonIsoMuon_1l_SYY = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_1l_SYY_Run2_MaxSig.cfg";
        std::string DoubleDisCo_Cfg_2l_SYY            = "Keras_Tensorflow_DoubleDisCo_Reg_2l_SYY_Run2_MaxSig.cfg";
        std::string DoubleDisCo_Model_2l_SYY          = "keras_frozen_DoubleDisCo_Reg_2l_SYY_Run2_MaxSig.pb";
        std::string DoubleDisCo_Cfg_NonIsoMuon_2l_SYY = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_2l_SYY_Run2_MaxSig.cfg";
        // SF root files
        std::string leptonFileName                    = "allInOne_leptonicSF_UL.root";
        std::string hadronicFileName                  = "allInOne_hadronicSF_UL.root";
        std::string btagEffFileName                   = "allInOne_BTagEff_UL.root";
        std::string meanFileName                      = "allInOne_SFMean_UL.root";
        std::string toptaggerFileName                 = "allInOne_TopTagEffandSF_UL_new.root";

        // Determines if an analyzer is compatible with fastmode
        // If it is not, then the fastMode flag is essentially neutralized
        bool fastModeCompatible = false;

        if(filetag.find("2016preVFP") != std::string::npos)
        {
            runYear                           = "2016preVFP";
            Lumi                              = 19520.0;
            // For b-tagging working point definitions, see link: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL16preVFP
            deepCSV_WP_loose                  = 0.2027;
            deepCSV_WP_medium                 = 0.6001;
            deepCSV_WP_tight                  = 0.8819;
            deepFlavour_WP_loose              = 0.0508;
            deepFlavour_WP_medium             = 0.2598;
            deepFlavour_WP_tight              = 0.6502;           
            resolvedTop_WP                    = 0.95;
            mergedTop_WP                      = 0.937;
            bjetTagFileName                   = "wp_deepJet_106XUL16preVFP_v2.csv";
            //bjetTagFileNameReshape            = "reshaping_deepJet_106XUL16preVFP_v2.csv";
            blind                             = true;
            TopTaggerCfg                      = "TopTaggerCfg_2016preVFP.cfg";
        }
        else if(filetag.find("2016postVFP") != std::string::npos)
        {
            runYear                           = "2016postVFP";
            Lumi                              = 16810.0;
            // For b-tagging working point definitions, see link: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL16postVFP
            deepCSV_WP_loose                  = 0.1918;
            deepCSV_WP_medium                 = 0.5847;
            deepCSV_WP_tight                  = 0.8767;
            deepFlavour_WP_loose              = 0.0480;
            deepFlavour_WP_medium             = 0.2489;
            deepFlavour_WP_tight              = 0.6377;           
            resolvedTop_WP                    = 0.95;
            mergedTop_WP                      = 0.937;
            bjetTagFileName                   = "wp_deepJet_106XUL16postVFP_v3.csv";
            //bjetTagFileNameReshape            = "reshaping_deepJet_106XUL16postVFP_v3.csv";
            blind                             = true;
            TopTaggerCfg                      = "TopTaggerCfg_2016postVFP.cfg";
        }
        else if(filetag.find("2017") != std::string::npos)
        {
            runYear                           = "2017";
            Lumi                              = 41480.0;
            // For b-tagging working point definitions, see link: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL17
            deepCSV_WP_loose                  = 0.1355;
            deepCSV_WP_medium                 = 0.4506;
            deepCSV_WP_tight                  = 0.7738;
            deepFlavour_WP_loose              = 0.0532;
            deepFlavour_WP_medium             = 0.3040;
            deepFlavour_WP_tight              = 0.7476;
            resolvedTop_WP                    = 0.95;
            mergedTop_WP                      = 0.895;
            bjetTagFileName                   = "wp_deepJet_106XUL17_v3.csv";
            //bjetTagFileNameReshape            = "reshaping_deepJet_106XUL17_v3.csv";
            blind                             = true;
            TopTaggerCfg                      = "TopTaggerCfg_2017.cfg";
        }
        else if(filetag.find("2018") != std::string::npos)
        {
            runYear                           = "2018";
            Lumi                              = 59830.0;
            Lumi_preHEM                       = 21071.0;
            Lumi_postHEM                      = 38654.0;
            // For b-tagging working point definitions, see link: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL18
            deepCSV_WP_loose                  = 0.1208;
            deepCSV_WP_medium                 = 0.4168;
            deepCSV_WP_tight                  = 0.7665;
            deepFlavour_WP_loose              = 0.0490;
            deepFlavour_WP_medium             = 0.2783;
            deepFlavour_WP_tight              = 0.7100;
            resolvedTop_WP                    = 0.95;
            mergedTop_WP                      = 0.895;
            bjetTagFileName                   = "wp_deepJet_106XUL18_v2.csv";
            //bjetTagFileNameReshape            = "reshaping_deepJet_106XUL18_v2.csv";
            blind                             = true;
            TopTaggerCfg                      = "TopTaggerCfg_2018.cfg";
        }

        // When making skims for top tag SF measurement,
        // need to consider full discriminant range of best tops
        // so neutralize the working point in this case
        if (analyzer == "MakeTopTagSFTree")
        {
            resolvedTop_WP = 0.0;
            mergedTop_WP   = 0.0;
        }

        tr.registerDerivedVar("runYear",                           runYear                          );
        tr.registerDerivedVar("Lumi",                              Lumi                             );
        tr.registerDerivedVar("Lumi_preHEM",                       Lumi_preHEM                      );
        tr.registerDerivedVar("Lumi_postHEM",                      Lumi_postHEM                     );
        tr.registerDerivedVar("deepCSV_WP_loose",                  deepCSV_WP_loose                 );
        tr.registerDerivedVar("deepCSV_WP_medium",                 deepCSV_WP_medium                );
        tr.registerDerivedVar("deepCSV_WP_tight",                  deepCSV_WP_tight                 );
        tr.registerDerivedVar("deepFlavour_WP_loose",              deepFlavour_WP_loose             );
        tr.registerDerivedVar("deepFlavour_WP_medium",             deepFlavour_WP_medium            );
        tr.registerDerivedVar("deepFlavour_WP_tight",              deepFlavour_WP_tight             );
        tr.registerDerivedVar("resolvedTop_WP",                    resolvedTop_WP                   );
        tr.registerDerivedVar("mergedTop_WP",                      mergedTop_WP                     );
        tr.registerDerivedVar("isSignal",                          isSignal                         );
        tr.registerDerivedVar("DoubleDisCo_Cfg_0l_RPV",            DoubleDisCo_Cfg_0l_RPV           );
        tr.registerDerivedVar("DoubleDisCo_Model_0l_RPV",          DoubleDisCo_Model_0l_RPV         );
        tr.registerDerivedVar("DoubleDisCo_Cfg_NonIsoMuon_0l_RPV", DoubleDisCo_Cfg_NonIsoMuon_0l_RPV);
        tr.registerDerivedVar("DoubleDisCo_Cfg_1l_RPV",            DoubleDisCo_Cfg_1l_RPV           );
        tr.registerDerivedVar("DoubleDisCo_Model_1l_RPV",          DoubleDisCo_Model_1l_RPV         );
        tr.registerDerivedVar("DoubleDisCo_Cfg_NonIsoMuon_1l_RPV", DoubleDisCo_Cfg_NonIsoMuon_1l_RPV);
        tr.registerDerivedVar("DoubleDisCo_Cfg_2l_RPV",            DoubleDisCo_Cfg_2l_RPV           );
        tr.registerDerivedVar("DoubleDisCo_Model_2l_RPV",          DoubleDisCo_Model_2l_RPV         );
        tr.registerDerivedVar("DoubleDisCo_Cfg_NonIsoMuon_2l_RPV", DoubleDisCo_Cfg_NonIsoMuon_2l_RPV);
        tr.registerDerivedVar("DoubleDisCo_Cfg_0l_SYY",            DoubleDisCo_Cfg_0l_SYY           );
        tr.registerDerivedVar("DoubleDisCo_Model_0l_SYY",          DoubleDisCo_Model_0l_SYY         );
        tr.registerDerivedVar("DoubleDisCo_Cfg_NonIsoMuon_0l_SYY", DoubleDisCo_Cfg_NonIsoMuon_0l_SYY);
        tr.registerDerivedVar("DoubleDisCo_Cfg_1l_SYY",            DoubleDisCo_Cfg_1l_SYY           );
        tr.registerDerivedVar("DoubleDisCo_Model_1l_SYY",          DoubleDisCo_Model_1l_SYY         );
        tr.registerDerivedVar("DoubleDisCo_Cfg_NonIsoMuon_1l_SYY", DoubleDisCo_Cfg_NonIsoMuon_1l_SYY);
        tr.registerDerivedVar("DoubleDisCo_Cfg_2l_SYY",            DoubleDisCo_Cfg_2l_SYY           );
        tr.registerDerivedVar("DoubleDisCo_Model_2l_SYY",          DoubleDisCo_Model_2l_SYY         );
        tr.registerDerivedVar("DoubleDisCo_Cfg_NonIsoMuon_2l_SYY", DoubleDisCo_Cfg_NonIsoMuon_2l_SYY);
        tr.registerDerivedVar("leptonFileName",                    leptonFileName                   );
        tr.registerDerivedVar("hadronicFileName",                  hadronicFileName                 );
        tr.registerDerivedVar("toptaggerFileName",                 toptaggerFileName                );
        tr.registerDerivedVar("btagEffFileName",                   btagEffFileName                  );
        tr.registerDerivedVar("bjetTagFileName",                   bjetTagFileName                  );
        //tr.registerDerivedVar("bjetTagFileNameReshape",            bjetTagFileNameReshape           );
        tr.registerDerivedVar("meanFileName",                      meanFileName                     );
        tr.registerDerivedVar("etaCut",                            2.4                              );
        tr.registerDerivedVar("blind",                             blind                            );
        tr.registerDerivedVar("TopTaggerCfg",                      TopTaggerCfg                     );

        // Register Modules that are needed for each Analyzer
        if(analyzer=="MakeNJetDists") // for legacy 1l
        {
            const std::vector<std::string> modulesList = {
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Jet",
                "BJet",
                "CommonVariables",
                "RunTopTagger",
                "Baseline",
                "MakeMVAVariables",
                "DeepEventShape",
                "BTagCorrector",
                "ScaleFactors"
            };
            registerModules(tr, std::move(modulesList));
        }
        else if(analyzer=="AnalyzeDoubleDisCo")
        {
            fastModeCompatible = true;
            const std::vector<std::string> modulesList = {
                "PrepNTupleVars",
                // All other necessary modules instantiated explicitly in AnalyzeDoubleDisCo !!!
            };
            registerModules(tr, std::move(modulesList));
        }
        else if(analyzer=="MakeNNVariables")
        {
            fastModeCompatible = true;
            const std::vector<std::string> modulesList = {
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Jet",
                "BJet",
                "CommonVariables",
                "RunTopTagger",
                "Baseline",
                "MakeMVAVariables",
                "StopJets",
                "MakeStopHemispheres_OldSeed",
                "MakeStopHemispheres_TopSeed",
                "StopGenMatch",
            };
            registerModules(tr, std::move(modulesList));
        }
        else if(analyzer=="AnalyzeLepTrigger" || analyzer=="HadTriggers_Analyzer")
        {
            const std::vector<std::string> modulesList = {
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Jet",
                "BJet",
                "CommonVariables",
                "RunTopTagger",
                "Baseline",
                "BTagCorrector",
                "ScaleFactors"
            };
            registerModules(tr, std::move(modulesList));
        }
        else if(analyzer=="CalculateBTagSF")
        {
            const std::vector<std::string> modulesList = {
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Jet",
                "BJet",
                "CommonVariables",
            };
            registerModules(tr, std::move(modulesList));
        }
        else if(analyzer=="CalculateTopTagSF")
        {
            const std::vector<std::string> modulesList = {
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Jet",
                "BJet",
                "CommonVariables",
                "RunTopTagger"
            };
            registerModules(tr, std::move(modulesList));
        }
        else if(analyzer=="CalculateSFMean")
        {
            const std::vector<std::string> modulesList = {
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Jet",
                "BJet",
                "CommonVariables",
                "ScaleFactors",
            };
            registerModules(tr, std::move(modulesList));
        }
        else if(analyzer=="HEM_Analyzer")
        {
            const std::vector<std::string> modulesList = {
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Jet",
                "BJet",
                "CommonVariables",
                "RunTopTagger",
                "Baseline",
                "BTagCorrector",
                "ScaleFactors",
                "MakeMVAVariables",
            };
            registerModules(tr, std::move(modulesList));
        }
        else if(analyzer=="TopTaggerSF_Analyzer")
        {
            const std::vector<std::string> modulesList = {
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Jet",
                "BJet",
                "CommonVariables",
                "RunTopTagger",
                "Baseline",
                "BTagCorrector",
                "ScaleFactors",
                "StopJets",
                "MakeMVAVariables",
                "MakeStopHemispheres_TopSeed",
                "DoubleDisCo_0l_RPV"
            };
            registerModules(tr, std::move(modulesList));
        }
        else if(analyzer=="ResolvedTopTagger_Analyzer")
        {
            const std::vector<std::string> modulesList = {
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Jet",
                "BJet",
                "CommonVariables",
                "RunTopTagger",
                "Baseline",
                "BTagCorrector",
                "ScaleFactors"
            };
            registerModules(tr, std::move(modulesList));
        }
        else if(analyzer=="Semra_Analyzer")
        {
            const std::vector<std::string> modulesList = {
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Jet",
                "BJet",
                "CommonVariables",
                "RunTopTagger",
                "Baseline",
                "BTagCorrector",
                "ScaleFactors",
                "StopJets",
                "MakeMVAVariables_0l",
                "MakeStopHemispheres_OldSeed",
                "MakeStopHemispheres_OldSeed_NonIsoMuon",
                "MakeStopHemispheres_TopSeed",
                "StopGenMatch",
            };
            registerModules(tr, std::move(modulesList));
        }
        else if(analyzer=="StealthHemispheres" || analyzer=="ISRJets_Analyzer")
        {
            const std::vector<std::string> modulesList = {
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Jet",
                "BJet",
                "CommonVariables",
                "RunTopTagger",
                "Baseline",
                "MakeMVAVariables_0l",
                "ISRJets",
                "StopJets",
                "MakeStopHemispheres_All",
                "MakeStopHemispheres_OldSeed",
                //"MakeStopHemispheres_OldSeed_maskedISR",
                "MakeStopHemispheres_TopSeed",
                //"MakeStopHemispheres_TopSeed_maskedISR",
                "StopGenMatch",
                "BTagCorrector",
                "ScaleFactors",
            };
            registerModules(tr, std::move(modulesList));
        }
        else if(analyzer=="AnalyzeXsec")
        {
            const std::vector<std::string> modulesList = {
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Jet",
                "BJet",
                "CommonVariables",
            };
            registerModules(tr, std::move(modulesList));
        }
        else if(analyzer=="MakeMiniTree")
        {
            const std::vector<std::string> modulesList = {
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Photon",
                "Jet",
                "BJet",
                "CommonVariables",
            };
            registerModules(tr, std::move(modulesList));
        }
        else if(analyzer=="MakeQCDValTree")
        {
            fastModeCompatible = true;
            const std::vector<std::string> modulesList = {
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Photon",
                "Jet",
                "BJet",
                "CommonVariables",
                "RunTopTagger",
                "Baseline",
                "MakeMVAVariables",
                "MakeMVAVariables_NonIsoMuon",
                "StopJets",
                "MakeStopHemispheres_TopSeed",
                "MakeStopHemispheres_OldSeed",
                "MakeStopHemispheres_TopSeed_NonIsoMuon",
                "MakeStopHemispheres_OldSeed_NonIsoMuon",
                "BTagCorrector",
                "ScaleFactors",
                "StopGenMatch",
                "DoubleDisCo_0l_RPV",
                "DoubleDisCo_NonIsoMuon_0l_RPV",
                "DoubleDisCo_0l_SYY",
                "DoubleDisCo_NonIsoMuon_0l_SYY",
                "DoubleDisCo_1l_RPV",
                "DoubleDisCo_NonIsoMuon_1l_RPV",
                "DoubleDisCo_1l_SYY",
                "DoubleDisCo_NonIsoMuon_1l_SYY",
            };
            registerModules(tr, std::move(modulesList));
        }
        else if(analyzer=="MakeTopTagSFTree")
        {
            fastModeCompatible = true;
            const std::vector<std::string> modulesList = {
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Photon",
                "Jet",
                "BJet",
                "CommonVariables",
                "RunTopTagger",
                "Baseline",
                "BTagCorrector",
                "ScaleFactors",
            };
            registerModules(tr, std::move(modulesList));
        }

        else if(analyzer=="MakeAnaSkimTree")
        {
            fastModeCompatible = true;
            const std::vector<std::string> modulesList = {
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Photon",
                "Jet",
                "BJet",
                "CommonVariables",
                "RunTopTagger",
                "Baseline",
                "MakeMVAVariables",
                "MakeMVAVariables_NonIsoMuon",
                "StopJets",
                "MakeStopHemispheres_TopSeed",
                "MakeStopHemispheres_OldSeed",
                "MakeStopHemispheres_OldSeed_NonIsoMuon",
                "BTagCorrector",
                "ScaleFactors",
                "StopGenMatch",
                "DoubleDisCo_0l_RPV",
                "DoubleDisCo_0l_SYY",
                "DoubleDisCo_1l_RPV",
                "DoubleDisCo_1l_SYY",
                "DoubleDisCo_2l_RPV",
                "DoubleDisCo_2l_SYY",
                "DoubleDisCo_NonIsoMuon_0l_RPV",
                "DoubleDisCo_NonIsoMuon_0l_SYY",
                "DoubleDisCo_NonIsoMuon_1l_RPV",
                "DoubleDisCo_NonIsoMuon_1l_SYY",
            };
            registerModules(tr, std::move(modulesList));
        }
        else
        {
            const std::vector<std::string> modulesList = {
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Jet",
                "BJet",
                "CommonVariables",
                "RunTopTagger",
                "Baseline",
                "BTagCorrector",
                "ScaleFactors",
                "StopGenMatch"
            };
            registerModules(tr, std::move(modulesList));
        }

        const auto& fastMode = tr.getVar<bool>("fastMode");
        if (fastMode and !fastModeCompatible)
        {
            std::cerr << utility::color("Error: Analyzer \"" + analyzer + "\" is not compatible with fast mode !!! Exiting...", "red") << std::endl;
            exit(-1);
        }
    }
};

#endif
