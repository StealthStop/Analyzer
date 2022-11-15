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
        const auto& DoubleDisCo_Cfg_2l_RPV            = tr.getVar<std::string>("DoubleDisCo_Cfg_2l_RPV"   ); 
        const auto& DoubleDisCo_Model_2l_RPV          = tr.getVar<std::string>("DoubleDisCo_Model_2l_RPV" );
        const auto& DoubleDisCo_Cfg_NonIsoMuon_2l_RPV = tr.getVar<std::string>("DoubleDisCo_Cfg_NonIsoMuon_2l_RPV");
        const auto& DoubleDisCo_Cfg_2l_SYY            = tr.getVar<std::string>("DoubleDisCo_Cfg_2l_SYY"   ); 
        const auto& DoubleDisCo_Model_2l_SYY          = tr.getVar<std::string>("DoubleDisCo_Model_2l_SYY" );
        const auto& DoubleDisCo_Cfg_NonIsoMuon_2l_SYY = tr.getVar<std::string>("DoubleDisCo_Cfg_NonIsoMuon_2l_SYY");
        const auto& leptonFileName            = tr.getVar<std::string>("leptonFileName"              );
        const auto& hadronicFileName          = tr.getVar<std::string>("hadronicFileName"            );
        const auto& bjetFileName              = tr.getVar<std::string>("bjetFileName"                );
        const auto& bjetCSVFileName           = tr.getVar<std::string>("bjetCSVFileName"             );
        //const auto& bjetCSVFileNameReshape    = tr.getVar<std::string>("bjetCSVFileNameReshape"      );
        const auto& filetag                   = tr.getVar<std::string>("filetag"                     );
        const auto& meanFileName              = tr.getVar<std::string>("meanFileName"                );
        const auto& TopTaggerCfg              = tr.getVar<std::string>("TopTaggerCfg"                );
        const auto& TopTaggerCfg_ResolvedOnly = tr.getVar<std::string>("TopTaggerCfg_ResolvedOnly"   );

 
        for(const auto& module : modules)
        {
            if     (module=="PrepNTupleVars")                        tr.emplaceModule<PrepNTupleVars>();
            else if(module=="RunTopTagger")                          tr.emplaceModule<RunTopTagger>(TopTaggerCfg);
            else if(module=="RunTopTagger_ResolvedOnly")             tr.emplaceModule<RunTopTagger>(TopTaggerCfg_ResolvedOnly);
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
                if     (module=="ScaleFactors")  tr.emplaceModule<ScaleFactors>(runYear, leptonFileName, hadronicFileName, meanFileName);
                else if(module=="BTagCorrector")
                {
                    std::string bjetCSVFileNameReshape = "";
                    auto& bTagCorrector = tr.emplaceModule<BTagCorrector>(bjetFileName, "", bjetCSVFileName, bjetCSVFileNameReshape, filetag);
                    bTagCorrector.SetVarNames("GenParticles_PdgId", "Jets", "GoodJets_pt30", "Jets_bJetTagDeepCSVtotb", "Jets_partonFlavor");
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
        std::string DoubleDisCo_Cfg_0l_RPV, DoubleDisCo_Model_0l_RPV, DoubleDisCo_Cfg_NonIsoMuon_0l_RPV; 
        std::string DoubleDisCo_Cfg_1l_RPV, DoubleDisCo_Model_1l_RPV, DoubleDisCo_Cfg_NonIsoMuon_1l_RPV; 
        std::string DoubleDisCo_Cfg_2l_RPV, DoubleDisCo_Model_2l_RPV, DoubleDisCo_Cfg_NonIsoMuon_2l_RPV; 
        std::string DoubleDisCo_Cfg_0l_SYY, DoubleDisCo_Model_0l_SYY, DoubleDisCo_Cfg_NonIsoMuon_0l_SYY; 
        std::string DoubleDisCo_Cfg_1l_SYY, DoubleDisCo_Model_1l_SYY, DoubleDisCo_Cfg_NonIsoMuon_1l_SYY; 
        std::string DoubleDisCo_Cfg_2l_SYY, DoubleDisCo_Model_2l_SYY, DoubleDisCo_Cfg_NonIsoMuon_2l_SYY; 
        std::string leptonFileName, hadronicFileName, bjetFileName, bjetCSVFileName, bjetCSVFileNameReshape, meanFileName, TopTaggerCfg, TopTaggerCfg_ResolvedOnly;
 
        double Lumi=0.0, Lumi_postHEM=-1.0, Lumi_preHEM=-1.0;
        double deepCSV_WP_loose=0.0, deepCSV_WP_medium=0.0, deepCSV_WP_tight=0.0;
        bool blind = true;

        // Determines if an analyzer is compatible with fastmode
        // If it is not, then the fastMode flag is essentially neutralized
        bool fastModeCompatible = false;

        if(filetag.find("2016preVFP") != std::string::npos)
        {
            runYear                           = "2016preVFP";
            Lumi                              = 19520.0;
            deepCSV_WP_loose                  = 0.2027;
            deepCSV_WP_medium                 = 0.6001;
            deepCSV_WP_tight                  = 0.8819;           
            DoubleDisCo_Cfg_0l_RPV            = "Keras_Tensorflow_DoubleDisCo_Reg_0l_RPV_Run2.cfg";           
            DoubleDisCo_Model_0l_RPV          = "keras_frozen_DoubleDisCo_Reg_0l_RPV_Run2.pb";
            DoubleDisCo_Cfg_NonIsoMuon_0l_RPV = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_0l_RPV_Run2.cfg"; 
            DoubleDisCo_Cfg_1l_RPV            = "Keras_Tensorflow_DoubleDisCo_Reg_1l_RPV_Run2.cfg";
            DoubleDisCo_Model_1l_RPV          = "keras_frozen_DoubleDisCo_Reg_1l_RPV_Run2.pb";
            DoubleDisCo_Cfg_NonIsoMuon_1l_RPV = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_1l_RPV_Run2.cfg";
            DoubleDisCo_Cfg_2l_RPV            = "Keras_Tensorflow_DoubleDisCo_Reg_2l_RPV_Run2.cfg";
            DoubleDisCo_Model_2l_RPV          = "keras_frozen_DoubleDisCo_Reg_2l_RPV_Run2.pb";
            DoubleDisCo_Cfg_NonIsoMuon_2l_RPV = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_2l_RPV_Run2.cfg";

            // Use RPV config for now --- switch to SYY with dedicated training
            DoubleDisCo_Cfg_0l_SYY            = "Keras_Tensorflow_DoubleDisCo_Reg_0l_RPV_Run2.cfg";           
            DoubleDisCo_Model_0l_SYY          = "keras_frozen_DoubleDisCo_Reg_0l_RPV_Run2.pb";
            DoubleDisCo_Cfg_NonIsoMuon_0l_SYY = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_0l_RPV_Run2.cfg";
            DoubleDisCo_Cfg_1l_SYY            = "Keras_Tensorflow_DoubleDisCo_Reg_1l_RPV_Run2.cfg";
            DoubleDisCo_Model_1l_SYY          = "keras_frozen_DoubleDisCo_Reg_1l_RPV_Run2.pb";
            DoubleDisCo_Cfg_NonIsoMuon_1l_SYY = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_1l_RPV_Run2.cfg";
            DoubleDisCo_Cfg_2l_SYY            = "Keras_Tensorflow_DoubleDisCo_Reg_2l_RPV_Run2.cfg";
            DoubleDisCo_Model_2l_SYY          = "keras_frozen_DoubleDisCo_Reg_2l_RPV_Run2.pb";
            DoubleDisCo_Cfg_NonIsoMuon_2l_SYY = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_2l_RPV_Run2.cfg";

            leptonFileName                    = "allInOne_leptonSF_UL.root";
            hadronicFileName                  = "allInOne_hadronicSF_UL.root";
            bjetFileName                      = "allInOne_BTagEff_UL.root";
            bjetCSVFileName                   = "wp_deepCSV_106XUL16preVFP_v2.csv";
            //bjetCSVFileNameReshape            = "reshaping_deepJet_106XUL16preVFP_v2.csv";
            meanFileName                      = "allInOne_SFMean_UL.root";
            blind                             = true;
            TopTaggerCfg                      = "TopTaggerCfg_2016preVFP.cfg";
            TopTaggerCfg_ResolvedOnly         = "TopTaggerCfg_ResolvedOnly_2016preVFP.cfg";

        }

        else if(filetag.find("2016postVFP") != std::string::npos)
        {
            runYear                           = "2016postVFP";
            Lumi                              = 16810.0;
            deepCSV_WP_loose                  = 0.1918;
            deepCSV_WP_medium                 = 0.5847;
            deepCSV_WP_tight                  = 0.8767;           
            DoubleDisCo_Cfg_0l_RPV            = "Keras_Tensorflow_DoubleDisCo_Reg_0l_RPV_Run2.cfg";           
            DoubleDisCo_Model_0l_RPV          = "keras_frozen_DoubleDisCo_Reg_0l_RPV_Run2.pb";
            DoubleDisCo_Cfg_NonIsoMuon_0l_RPV = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_0l_RPV_Run2.cfg";
            DoubleDisCo_Cfg_1l_RPV            = "Keras_Tensorflow_DoubleDisCo_Reg_1l_RPV_Run2.cfg";
            DoubleDisCo_Model_1l_RPV          = "keras_frozen_DoubleDisCo_Reg_1l_RPV_Run2.pb";
            DoubleDisCo_Cfg_NonIsoMuon_1l_RPV = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_1l_RPV_Run2.cfg";
            DoubleDisCo_Cfg_2l_RPV            = "Keras_Tensorflow_DoubleDisCo_Reg_2l_RPV_Run2.cfg";
            DoubleDisCo_Model_2l_RPV          = "keras_frozen_DoubleDisCo_Reg_2l_RPV_Run2.pb";
            DoubleDisCo_Cfg_NonIsoMuon_2l_RPV = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_2l_RPV_Run2.cfg";

            // SYY config for now --- switch to SYY with dedicated training
            DoubleDisCo_Cfg_0l_SYY            = "Keras_Tensorflow_DoubleDisCo_Reg_0l_RPV_Run2.cfg";           
            DoubleDisCo_Model_0l_SYY          = "keras_frozen_DoubleDisCo_Reg_0l_RPV_Run2.pb";
            DoubleDisCo_Cfg_NonIsoMuon_0l_SYY = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_0l_RPV_Run2.cfg";
            DoubleDisCo_Cfg_1l_SYY            = "Keras_Tensorflow_DoubleDisCo_Reg_1l_RPV_Run2.cfg";
            DoubleDisCo_Model_1l_SYY          = "keras_frozen_DoubleDisCo_Reg_1l_RPV_Run2.pb";
            DoubleDisCo_Cfg_NonIsoMuon_1l_SYY = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_1l_RPV_Run2.cfg";
            DoubleDisCo_Cfg_2l_SYY            = "Keras_Tensorflow_DoubleDisCo_Reg_2l_RPV_Run2.cfg";
            DoubleDisCo_Model_2l_SYY          = "keras_frozen_DoubleDisCo_Reg_2l_RPV_Run2.pb";
            DoubleDisCo_Cfg_NonIsoMuon_2l_SYY = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_2l_RPV_Run2.cfg";

            // SF root files
            leptonFileName                    = "allInOne_leptonSF_UL.root";
            hadronicFileName                  = "allInOne_hadronicSF_UL.root";
            bjetFileName                      = "allInOne_BTagEff_UL.root";
            bjetCSVFileName                   = "wp_deepCSV_106XUL16postVFP_v3.csv";
            //bjetCSVFileNameReshape            = "reshaping_deepJet_106XUL16postVFP_v3.csv";
            meanFileName                      = "allInOne_SFMean_UL.root";
            blind                             = true;
            TopTaggerCfg                      = "TopTaggerCfg_2016postVFP.cfg";
            TopTaggerCfg_ResolvedOnly         = "TopTaggerCfg_ResolvedOnly_2016postVFP.cfg";

        }
        else if(filetag.find("2017") != std::string::npos)
        { 
            runYear                           = "2017";
            Lumi                              = 41480.0;
            deepCSV_WP_loose                  = 0.1355;
            deepCSV_WP_medium                 = 0.4506;       
            deepCSV_WP_tight                  = 0.7738;
            DoubleDisCo_Cfg_0l_RPV            = "Keras_Tensorflow_DoubleDisCo_Reg_0l_RPV_Run2.cfg";           
            DoubleDisCo_Model_0l_RPV          = "keras_frozen_DoubleDisCo_Reg_0l_RPV_Run2.pb";
            DoubleDisCo_Cfg_NonIsoMuon_0l_RPV = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_0l_RPV_Run2.cfg";
            DoubleDisCo_Cfg_1l_RPV            = "Keras_Tensorflow_DoubleDisCo_Reg_1l_RPV_Run2.cfg";
            DoubleDisCo_Model_1l_RPV          = "keras_frozen_DoubleDisCo_Reg_1l_RPV_Run2.pb";
            DoubleDisCo_Cfg_NonIsoMuon_1l_RPV = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_1l_RPV_Run2.cfg";
            DoubleDisCo_Cfg_2l_RPV            = "Keras_Tensorflow_DoubleDisCo_Reg_2l_RPV_Run2.cfg";
            DoubleDisCo_Model_2l_RPV          = "keras_frozen_DoubleDisCo_Reg_2l_RPV_Run2.pb";
            DoubleDisCo_Cfg_NonIsoMuon_2l_RPV = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_2l_RPV_Run2.cfg";

            // Use RPV config for now --- switch to SYY with dedicated training
            DoubleDisCo_Cfg_0l_SYY            = "Keras_Tensorflow_DoubleDisCo_Reg_0l_RPV_Run2.cfg";           
            DoubleDisCo_Model_0l_SYY          = "keras_frozen_DoubleDisCo_Reg_0l_RPV_Run2.pb";
            DoubleDisCo_Cfg_NonIsoMuon_0l_SYY = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_0l_RPV_Run2.cfg";
            DoubleDisCo_Cfg_1l_SYY            = "Keras_Tensorflow_DoubleDisCo_Reg_1l_RPV_Run2.cfg";
            DoubleDisCo_Model_1l_SYY          = "keras_frozen_DoubleDisCo_Reg_1l_RPV_Run2.pb";
            DoubleDisCo_Cfg_NonIsoMuon_1l_SYY = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_1l_RPV_Run2.cfg";
            DoubleDisCo_Cfg_2l_SYY            = "Keras_Tensorflow_DoubleDisCo_Reg_2l_RPV_Run2.cfg";
            DoubleDisCo_Model_2l_SYY          = "keras_frozen_DoubleDisCo_Reg_2l_RPV_Run2.pb";
            DoubleDisCo_Cfg_NonIsoMuon_2l_SYY = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_2l_RPV_Run2.cfg";

            // SF root files
            leptonFileName                    = "allInOne_leptonSF_UL.root";
            hadronicFileName                  = "allInOne_hadronicSF_UL.root";
            bjetFileName                      = "allInOne_BTagEff_UL.root";
            bjetCSVFileName                   = "wp_deepCSV_106XUL17_v3.csv";
            //bjetCSVFileNameReshape            = "reshaping_deepJet_106XUL17_v3.csv";
            meanFileName                      = "allInOne_SFMean_UL.root";
            blind                             = true;
            TopTaggerCfg                      = "TopTaggerCfg_2017.cfg";
            TopTaggerCfg_ResolvedOnly         = "TopTaggerCfg_ResolvedOnly_2017.cfg";

        }
        else if(filetag.find("2018") != std::string::npos) 
        {
            runYear                           = "2018";
            Lumi                              = 59830.0;
            Lumi_preHEM                       = 21071.0;
            Lumi_postHEM                      = 38654.0;
            deepCSV_WP_loose                  = 0.1208;
            deepCSV_WP_medium                 = 0.4168;       
            deepCSV_WP_tight                  = 0.7665;
            DoubleDisCo_Cfg_0l_RPV            = "Keras_Tensorflow_DoubleDisCo_Reg_0l_RPV_Run2.cfg";           
            DoubleDisCo_Model_0l_RPV          = "keras_frozen_DoubleDisCo_Reg_0l_RPV_Run2.pb";
            DoubleDisCo_Cfg_NonIsoMuon_0l_RPV = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_0l_RPV_Run2.cfg";
            DoubleDisCo_Cfg_1l_RPV            = "Keras_Tensorflow_DoubleDisCo_Reg_1l_RPV_Run2.cfg";
            DoubleDisCo_Model_1l_RPV          = "keras_frozen_DoubleDisCo_Reg_1l_RPV_Run2.pb";
            DoubleDisCo_Cfg_NonIsoMuon_1l_RPV = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_1l_RPV_Run2.cfg";
            DoubleDisCo_Cfg_2l_RPV            = "Keras_Tensorflow_DoubleDisCo_Reg_2l_RPV_Run2.cfg";
            DoubleDisCo_Model_2l_RPV          = "keras_frozen_DoubleDisCo_Reg_2l_RPV_Run2.pb";
            DoubleDisCo_Cfg_NonIsoMuon_2l_RPV = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_2l_RPV_Run2.cfg";

            // Use RPV config for now --- switch to SYY with dedicated training
            DoubleDisCo_Cfg_0l_SYY            = "Keras_Tensorflow_DoubleDisCo_Reg_0l_RPV_Run2.cfg";           
            DoubleDisCo_Model_0l_SYY          = "keras_frozen_DoubleDisCo_Reg_0l_RPV_Run2.pb";
            DoubleDisCo_Cfg_NonIsoMuon_0l_SYY = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_0l_RPV_Run2.cfg";
            DoubleDisCo_Cfg_1l_SYY            = "Keras_Tensorflow_DoubleDisCo_Reg_1l_RPV_Run2.cfg";
            DoubleDisCo_Model_1l_SYY          = "keras_frozen_DoubleDisCo_Reg_1l_RPV_Run2.pb";
            DoubleDisCo_Cfg_NonIsoMuon_1l_SYY = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_1l_RPV_Run2.cfg";
            DoubleDisCo_Cfg_2l_SYY            = "Keras_Tensorflow_DoubleDisCo_Reg_2l_RPV_Run2.cfg";
            DoubleDisCo_Model_2l_SYY          = "keras_frozen_DoubleDisCo_Reg_2l_RPV_Run2.pb";
            DoubleDisCo_Cfg_NonIsoMuon_2l_SYY = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_2l_RPV_Run2.cfg";

            // SF root files
            leptonFileName                    = "allInOne_leptonSF_UL.root";
            hadronicFileName                  = "allInOne_hadronicSF_UL.root";
            bjetFileName                      = "allInOne_BTagEff_UL.root";
            bjetCSVFileName                   = "wp_deepCSV_106XUL18_v2.csv";
            //bjetCSVFileNameReshape            = "reshaping_deepJet_106XUL18_v2.csv";
            meanFileName                      = "allInOne_SFMean_UL.root";
            blind                             = true;
            TopTaggerCfg                      = "TopTaggerCfg_2018.cfg";
            TopTaggerCfg_ResolvedOnly         = "TopTaggerCfg_ResolvedOnly_2018.cfg";

        }

        tr.registerDerivedVar("runYear",                           runYear                          );
        tr.registerDerivedVar("Lumi",                              Lumi                             );
        tr.registerDerivedVar("Lumi_preHEM",                       Lumi_preHEM                      );
        tr.registerDerivedVar("Lumi_postHEM",                      Lumi_postHEM                     );
        tr.registerDerivedVar("deepCSV_WP_loose",                  deepCSV_WP_loose                 );
        tr.registerDerivedVar("deepCSV_WP_medium",                 deepCSV_WP_medium                );
        tr.registerDerivedVar("deepCSV_WP_tight",                  deepCSV_WP_tight                 );
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
        tr.registerDerivedVar("bjetFileName",                      bjetFileName                     );        
        tr.registerDerivedVar("bjetCSVFileName",                   bjetCSVFileName                  );        
        //tr.registerDerivedVar("bjetCSVFileNameReshape",            bjetCSVFileNameReshape           );        
        tr.registerDerivedVar("meanFileName",                      meanFileName                     );        
        tr.registerDerivedVar("etaCut",                            2.4                              ); 
        tr.registerDerivedVar("blind",                             blind                            );
        tr.registerDerivedVar("TopTaggerCfg",                      TopTaggerCfg                     );
        tr.registerDerivedVar("TopTaggerCfg_ResolvedOnly",         TopTaggerCfg_ResolvedOnly        );

        // Register Modules that are needed for each Analyzer
        if(analyzer=="MakeNJetDists") // for legacy 1l
        {
            const std::vector<std::string> modulesList = {
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Jet",
                "BJet",
                "RunTopTagger",
                "CommonVariables",
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
            const std::vector<std::string> modulesList = {
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Jet",
                "BJet",
                "RunTopTagger",
                "CommonVariables",
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
                "RunTopTagger",
                "Baseline",
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
                "RunTopTagger",
                "Baseline",
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
                "RunTopTagger",
                "CommonVariables",
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
                "RunTopTagger",
                "CommonVariables",
                "Baseline",
                "BTagCorrector",
                "ScaleFactors"
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
                "RunTopTagger_ResolvedOnly",
                "CommonVariables",
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
                "RunTopTagger",
                "CommonVariables",
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
                "RunTopTagger",
                "CommonVariables",
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
        else
        {
            const std::vector<std::string> modulesList = {
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Jet",
                "BJet",
                "RunTopTagger",
                "CommonVariables",
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
