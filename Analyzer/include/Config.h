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
#include "Framework/Framework/include/PartialUnBlinding.h"
#include "Framework/Framework/include/StopGenMatch.h"
#include "Framework/Framework/include/MegaJetCombine.h"
#include "Framework/Framework/include/TrainingNTupleVars.h"
#include "Framework/Framework/include/MakeStopHemispheres.h"
#include "Framework/Framework/include/StopJets.h"
#include "Framework/Framework/include/ISRJets.h"
#include "Framework/Framework/include/FatJetCombine.h"

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

        //const auto& DoubleDisCo_Cfg_2l    = tr.getVar<std::string>("DoubleDisCo_Cfg_2l"   ); 
        //const auto& DoubleDisCo_Model_2l  = tr.getVar<std::string>("DoubleDisCo_Model_2l" );
        //const auto& DoubleDisCo_Cfg_NonIsoMuon_2l = tr.getVar<std::string>("DoubleDisCo_Cfg_NonIsoMuon_2l");
        const auto& leptonFileName        = tr.getVar<std::string>("leptonFileName"       );
        const auto& bjetFileName          = tr.getVar<std::string>("bjetFileName"         );
        const auto& bjetCSVFileName       = tr.getVar<std::string>("bjetCSVFileName"      );
        const auto& filetag               = tr.getVar<std::string>("filetag"              );
        const auto& meanFileName          = tr.getVar<std::string>("meanFileName"         );
        const auto& TopTaggerCfg          = tr.getVar<std::string>("TopTaggerCfg"         );
 
        for(const auto& module : modules)
        {
            if     (module=="PartialUnBlinding")                     tr.emplaceModule<PartialUnBlinding>();
            else if(module=="PrepNTupleVars")                        tr.emplaceModule<PrepNTupleVars>();
            else if(module=="RunTopTagger")                          tr.emplaceModule<RunTopTagger>(TopTaggerCfg);
            else if(module=="Muon")                                  tr.emplaceModule<Muon>();
            else if(module=="Electron")                              tr.emplaceModule<Electron>();
            else if(module=="Photon")                                tr.emplaceModule<Photon>();
            else if(module=="Jet")                                   tr.emplaceModule<Jet>();
            else if(module=="BJet")                                  tr.emplaceModule<BJet>();
            else if(module=="CommonVariables")                       tr.emplaceModule<CommonVariables>();
            else if(module=="MakeMVAVariables")                      tr.emplaceModule<MakeMVAVariables>(false, "", "GoodJets_pt30",                  false, true, 12, 2, "");
            else if(module=="MakeMVAVariables_NonIsoMuon")           tr.emplaceModule<MakeMVAVariables>(false, "", "NonIsoMuonJets_pt30",            false, true, 12, 2, "");
            else if(module=="MakeMVAVariables_0l_old")               tr.emplaceModule<MakeMVAVariables>(false, "", "GoodJets_pt45",                  false, true, 7,  2, "_0l");
            else if(module=="MakeMVAVariables_0l")                   tr.emplaceModule<MakeMVAVariables>(false, "", "GoodJets_pt30",                  false, true, 7,  2, "_0l");
            else if(module=="MakeMVAVariables_NonIsoMuon_0l")        tr.emplaceModule<MakeMVAVariables>(false, "", "GoodJets_pt30",                  false, true, 7,  2, "_0l");
            else if(module=="MakeMVAVariables_1l")                   tr.emplaceModule<MakeMVAVariables>(false, "", "GoodJets_pt30",                  false, true, 7,  2, "_1l");
            else if(module=="MakeMVAVariables_NonIsoMuon_1l")        tr.emplaceModule<MakeMVAVariables>(false, "", "NonIsoMuonJets_pt30",            false, true, 7,  2, "_1l");
            else if(module=="MakeMVAVariables_2l")                   tr.emplaceModule<MakeMVAVariables>(false, "", "GoodJets_pt30_GoodLeptons_pt20", false, true, 7,  2, "_2l");
            else if(module=="Baseline")                              tr.emplaceModule<Baseline>();
            else if(module=="StopGenMatch")                          tr.emplaceModule<StopGenMatch>();
            else if(module=="MegaJetCombine")                        tr.emplaceModule<MegaJetCombine>();
            else if(module=="FatJetCombine")                         tr.emplaceModule<FatJetCombine>();
            else if(module=="TrainingNTupleVars")                    tr.emplaceModule<TrainingNTupleVars>();
            else if(module=="MakeStopHemispheres_All")               tr.emplaceModule<MakeStopHemispheres>("Jets",     "AllJets",                 "NJets",                    "_All",                "", Hemisphere::InvMassSeed);
            else if(module=="MakeStopHemispheres_OldSeed")           tr.emplaceModule<MakeStopHemispheres>("Jets",     "GoodJets_pt20",           "NGoodJets_pt20",           "_OldSeed",            "", Hemisphere::InvMassSeed);
            else if(module=="MakeStopHemispheres_OldSeed_maskedISR") tr.emplaceModule<MakeStopHemispheres>("Jets",     "GoodJets_pt20_maskedISR", "NGoodJets_pt20_maskedISR", "_OldSeed_maskedISR",  "", Hemisphere::InvMassSeed);
            else if(module=="MakeStopHemispheres_OldSeed_NonIsoMuon")tr.emplaceModule<MakeStopHemispheres>("Jets",     "NonIsoMuonJets_pt20",     "NNonIsoMuonJets_pt30",     "_OldSeed_NonIsoMuon", "", Hemisphere::InvMassSeed);
            else if(module=="MakeStopHemispheres_TopSeed")           tr.emplaceModule<MakeStopHemispheres>("StopJets", "GoodStopJets",            "NGoodStopJets",            "_TopSeed",            "", Hemisphere::TopSeed);
            else if(module=="MakeStopHemispheres_TopSeed_maskedISR") tr.emplaceModule<MakeStopHemispheres>("StopJets", "GoodStopJets_maskedISR",  "NGoodStopJets_maskedISR",  "_TopSeed_maskedISR",  "", Hemisphere::TopSeed);
            else if(module=="StopJets")                              tr.emplaceModule<StopJets>();
            else if(module=="ISRJets")                               tr.emplaceModule<ISRJets>();
            else if(module=="DoubleDisCo_0l_RPV")                    tr.emplaceModule<DeepEventShape>(DoubleDisCo_Cfg_0l_RPV, DoubleDisCo_Model_0l_RPV);
            else if(module=="DoubleDisCo_NonIsoMuon_0l_RPV")         tr.emplaceModule<DeepEventShape>(DoubleDisCo_Cfg_NonIsoMuon_0l_RPV, DoubleDisCo_Model_0l_RPV);
            else if(module=="DoubleDisCo_1l_RPV")                    tr.emplaceModule<DeepEventShape>(DoubleDisCo_Cfg_1l_RPV, DoubleDisCo_Model_1l_RPV); 
            else if(module=="DoubleDisCo_NonIsoMuon_1l_RPV")         tr.emplaceModule<DeepEventShape>(DoubleDisCo_Cfg_NonIsoMuon_1l_RPV, DoubleDisCo_Model_1l_RPV);
            else if(module=="DoubleDisCo_0l_SYY")                    tr.emplaceModule<DeepEventShape>(DoubleDisCo_Cfg_0l_SYY, DoubleDisCo_Model_0l_SYY);
            else if(module=="DoubleDisCo_NonIsoMuon_0l_SYY")         tr.emplaceModule<DeepEventShape>(DoubleDisCo_Cfg_NonIsoMuon_0l_SYY, DoubleDisCo_Model_0l_SYY);
            else if(module=="DoubleDisCo_1l_SYY")                    tr.emplaceModule<DeepEventShape>(DoubleDisCo_Cfg_1l_SYY, DoubleDisCo_Model_1l_SYY); 
            else if(module=="DoubleDisCo_NonIsoMuon_1l_SYY")         tr.emplaceModule<DeepEventShape>(DoubleDisCo_Cfg_NonIsoMuon_1l_SYY, DoubleDisCo_Model_1l_SYY);

            //else if(module=="DoubleDisCo_2l")                        tr.emplaceModule<DeepEventShape>(DoubleDisCo_Cfg_2l, DoubleDisCo_Model_2l);
 
            if(runtype == "MC")
            {
                if     (module=="ScaleFactors")  tr.emplaceModule<ScaleFactors>(runYear, leptonFileName, meanFileName);
                else if(module=="BTagCorrector")
                {
                    auto& bTagCorrector = tr.emplaceModule<BTagCorrector>(bjetFileName, "", bjetCSVFileName, filetag);
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
        std::string DoubleDisCo_Cfg_0l_SYY, DoubleDisCo_Model_0l_SYY, DoubleDisCo_Cfg_NonIsoMuon_0l_SYY; 
        std::string DoubleDisCo_Cfg_1l_SYY, DoubleDisCo_Model_1l_SYY, DoubleDisCo_Cfg_NonIsoMuon_1l_SYY; 
        //std::string DoubleDisCo_Cfg_2l, DoubleDisCo_Model_2l, DoubleDisCo_Cfg_NonIsoMuon_2l;      
        std::string leptonFileName, bjetFileName, bjetCSVFileName, meanFileName, TopTaggerCfg;
 
        double Lumi=0.0, deepCSV_WP_loose=0.0, deepCSV_WP_medium=0.0, deepCSV_WP_tight=0.0;
        bool blind = true;

        if(filetag.find("2016preVFP") != std::string::npos)
        {
            runYear                           = "2016preVFP";
            Lumi                              = 19520.0;
            deepCSV_WP_loose                  = 0.2027;
            deepCSV_WP_medium                 = 0.6001;
            deepCSV_WP_tight                  = 0.8819;            
            DoubleDisCo_Cfg_0l_RPV            = "Keras_Tensorflow_DoubleDisCo_Reg_0l_RPV_2016.cfg";           
            DoubleDisCo_Model_0l_RPV          = "keras_frozen_DoubleDisCo_Reg_0l_RPV_2016.pb";
            DoubleDisCo_Cfg_NonIsoMuon_0l_RPV = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_0l_RPV_2016.cfg";
            DoubleDisCo_Cfg_1l_RPV            = "Keras_Tensorflow_DoubleDisCo_Reg_1l_RPV_2016.cfg";
            DoubleDisCo_Model_1l_RPV          = "keras_frozen_DoubleDisCo_Reg_1l_RPV_2016.pb";
            DoubleDisCo_Cfg_NonIsoMuon_1l_RPV = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_1l_RPV_2016.cfg";

            // Use RPV config for now --- switch to SYY with dedicated training
            DoubleDisCo_Cfg_0l_SYY            = "Keras_Tensorflow_DoubleDisCo_Reg_0l_RPV_2016.cfg";           
            DoubleDisCo_Model_0l_SYY          = "keras_frozen_DoubleDisCo_Reg_0l_RPV_2016.pb";
            DoubleDisCo_Cfg_NonIsoMuon_0l_SYY = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_0l_RPV_2016.cfg";
            DoubleDisCo_Cfg_1l_SYY            = "Keras_Tensorflow_DoubleDisCo_Reg_1l_RPV_2016.cfg";
            DoubleDisCo_Model_1l_SYY          = "keras_frozen_DoubleDisCo_Reg_1l_RPV_2016.pb";
            DoubleDisCo_Cfg_NonIsoMuon_1l_SYY = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_1l_RPV_2016.cfg";

            //DoubleDisCo_Cfg_2l    = "Keras_Tensorflow_DoubleDisCo_Reg_2l_2016.cfg";
            //DoubleDisCo_Model_2l  = "keras_frozen_DoubleDisCo_Reg_2l_2016.pb";
            //DoubleDisCo_Cfg_NonIsoMuon_2l = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_2l_2016.cfg";
            leptonFileName        = "allInOne_leptonSF_UL.root";
            bjetFileName          = "allInOne_BTagEff_UL.root";
            bjetCSVFileName       = "wp_deepCSV_106XUL16preVFP_v2.csv";
            meanFileName          = "allInOne_SFMean_UL.root";
            blind                 = false;
            TopTaggerCfg          = "TopTaggerCfg_2016.cfg";
        }

        else if(filetag.find("2016postVFP") != std::string::npos)
        {
            runYear                           = "2016postVFP";
            Lumi                              = 16810.0;
            deepCSV_WP_loose                  = 0.1918;
            deepCSV_WP_medium                 = 0.5847;
            deepCSV_WP_tight                  = 0.8767;            
            DoubleDisCo_Cfg_0l_RPV            = "Keras_Tensorflow_DoubleDisCo_Reg_0l_RPV_2016.cfg";           
            DoubleDisCo_Model_0l_RPV          = "keras_frozen_DoubleDisCo_Reg_0l_RPV_2016.pb";
            DoubleDisCo_Cfg_NonIsoMuon_0l_RPV = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_0l_RPV_2016.cfg";
            DoubleDisCo_Cfg_1l_RPV            = "Keras_Tensorflow_DoubleDisCo_Reg_1l_RPV_2016.cfg";
            DoubleDisCo_Model_1l_RPV          = "keras_frozen_DoubleDisCo_Reg_1l_RPV_2016.pb";
            DoubleDisCo_Cfg_NonIsoMuon_1l_RPV = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_1l_RPV_2016.cfg";

            // Use RPV config for now --- switch to SYY with dedicated training
            DoubleDisCo_Cfg_0l_SYY            = "Keras_Tensorflow_DoubleDisCo_Reg_0l_RPV_2016.cfg";           
            DoubleDisCo_Model_0l_SYY          = "keras_frozen_DoubleDisCo_Reg_0l_RPV_2016.pb";
            DoubleDisCo_Cfg_NonIsoMuon_0l_SYY = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_0l_RPV_2016.cfg";
            DoubleDisCo_Cfg_1l_SYY            = "Keras_Tensorflow_DoubleDisCo_Reg_1l_RPV_2016.cfg";
            DoubleDisCo_Model_1l_SYY          = "keras_frozen_DoubleDisCo_Reg_1l_RPV_2016.pb";
            DoubleDisCo_Cfg_NonIsoMuon_1l_SYY = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_1l_RPV_2016.cfg";

            //DoubleDisCo_Cfg_2l    = "Keras_Tensorflow_DoubleDisCo_Reg_2l_2016.cfg";
            //DoubleDisCo_Model_2l  = "keras_frozen_DoubleDisCo_Reg_2l_2016.pb";
            //DoubleDisCo_Cfg_NonIsoMuon_2l = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_2l_2016.cfg";
            leptonFileName        = "allInOne_leptonSF_UL.root";
            bjetFileName          = "allInOne_BTagEff_UL.root";
            bjetCSVFileName       = "wp_deepCSV_106XUL16postVFP_v3.csv";
            meanFileName          = "allInOne_SFMean_UL.root";
            blind                 = false;
            TopTaggerCfg          = "TopTaggerCfg_2016.cfg";
        }
        else if(filetag.find("2017") != std::string::npos)
        { 
            runYear               = "2017";
            Lumi                  = 41480.0;
            deepCSV_WP_loose      = 0.1355;
            deepCSV_WP_medium     = 0.4506;       
            deepCSV_WP_tight      = 0.7738;

            // Switch to proper year's training when done for UL
            DoubleDisCo_Cfg_0l_RPV            = "Keras_Tensorflow_DoubleDisCo_Reg_0l_RPV_2016.cfg";           
            DoubleDisCo_Model_0l_RPV          = "keras_frozen_DoubleDisCo_Reg_0l_RPV_2016.pb";
            DoubleDisCo_Cfg_NonIsoMuon_0l_RPV = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_0l_RPV_2016.cfg";
            DoubleDisCo_Cfg_1l_RPV            = "Keras_Tensorflow_DoubleDisCo_Reg_1l_RPV_2016.cfg";
            DoubleDisCo_Model_1l_RPV          = "keras_frozen_DoubleDisCo_Reg_1l_RPV_2016.pb";
            DoubleDisCo_Cfg_NonIsoMuon_1l_RPV = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_1l_RPV_2016.cfg";

            // Use RPV config for now --- switch to SYY with dedicated training
            DoubleDisCo_Cfg_0l_SYY            = "Keras_Tensorflow_DoubleDisCo_Reg_0l_RPV_2016.cfg";           
            DoubleDisCo_Model_0l_SYY          = "keras_frozen_DoubleDisCo_Reg_0l_RPV_2016.pb";
            DoubleDisCo_Cfg_NonIsoMuon_0l_SYY = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_0l_RPV_2016.cfg";
            DoubleDisCo_Cfg_1l_SYY            = "Keras_Tensorflow_DoubleDisCo_Reg_1l_RPV_2016.cfg";
            DoubleDisCo_Model_1l_SYY          = "keras_frozen_DoubleDisCo_Reg_1l_RPV_2016.pb";
            DoubleDisCo_Cfg_NonIsoMuon_1l_SYY = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_1l_RPV_2016.cfg";

            //DoubleDisCo_Cfg_2l    = "Keras_Tensorflow_DoubleDisCo_Reg_2l_2016.cfg";
            //DoubleDisCo_Model_2l  = "keras_frozen_DoubleDisCo_Reg_2l_2016.pb";
            //DoubleDisCo_Cfg_NonIsoMuon_2l = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_2l_2016.cfg";
            leptonFileName        = "allInOne_leptonSF_UL.root";
            bjetFileName          = "allInOne_BTagEff_UL.root";
            bjetCSVFileName       = "wp_deepCSV_106XUL17_v3.csv";
            meanFileName          = "allInOne_SFMean_UL.root";
            blind                 = false;
            TopTaggerCfg          = "TopTaggerCfg_2017.cfg";
        }
        else if(filetag.find("2018pre") != std::string::npos) 
        {
            runYear               = "2018pre";
            Lumi                  = 21071.0;
            deepCSV_WP_loose      = 0.1208;
            deepCSV_WP_medium     = 0.4168;       
            deepCSV_WP_tight      = 0.7665;

            // Switch to proper year's training when done for UL
            DoubleDisCo_Cfg_0l_RPV            = "Keras_Tensorflow_DoubleDisCo_Reg_0l_RPV_2016.cfg";           
            DoubleDisCo_Model_0l_RPV          = "keras_frozen_DoubleDisCo_Reg_0l_RPV_2016.pb";
            DoubleDisCo_Cfg_NonIsoMuon_0l_RPV = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_0l_RPV_2016.cfg";
            DoubleDisCo_Cfg_1l_RPV            = "Keras_Tensorflow_DoubleDisCo_Reg_1l_RPV_2016.cfg";
            DoubleDisCo_Model_1l_RPV          = "keras_frozen_DoubleDisCo_Reg_1l_RPV_2016.pb";
            DoubleDisCo_Cfg_NonIsoMuon_1l_RPV = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_1l_RPV_2016.cfg";

            // Use RPV config for now --- switch to SYY with dedicated training
            DoubleDisCo_Cfg_0l_SYY            = "Keras_Tensorflow_DoubleDisCo_Reg_0l_RPV_2016.cfg";           
            DoubleDisCo_Model_0l_SYY          = "keras_frozen_DoubleDisCo_Reg_0l_RPV_2016.pb";
            DoubleDisCo_Cfg_NonIsoMuon_0l_SYY = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_0l_RPV_2016.cfg";
            DoubleDisCo_Cfg_1l_SYY            = "Keras_Tensorflow_DoubleDisCo_Reg_1l_RPV_2016.cfg";
            DoubleDisCo_Model_1l_SYY          = "keras_frozen_DoubleDisCo_Reg_1l_RPV_2016.pb";
            DoubleDisCo_Cfg_NonIsoMuon_1l_SYY = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_1l_RPV_2016.cfg";

            //DoubleDisCo_Cfg_2l    = "Keras_Tensorflow_DoubleDisCo_Reg_2l_2016.cfg";
            //DoubleDisCo_Model_2l  = "keras_frozen_DoubleDisCo_Reg_2l_2016.pb";
            //DoubleDisCo_Cfg_NonIsoMuon_2l = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_2l_2016.cfg";
            leptonFileName        = "allInOne_leptonSF_UL.root";
            bjetFileName          = "allInOne_BTagEff_UL.root";
            bjetCSVFileName       = "wp_deepCSV_106XUL18_v2.csv";
            meanFileName          = "allInOne_SFMean_UL.root";
            blind                 = false;
            TopTaggerCfg          = "TopTaggerCfg_2018.cfg";
        }
        else if(filetag.find("2018post") != std::string::npos) 
        {
            runYear               = "2018post";
            Lumi                  = 38654.0;
            deepCSV_WP_loose      = 0.1208;
            deepCSV_WP_medium     = 0.4168;       
            deepCSV_WP_tight      = 0.7665;

            // Switch to proper year's training when done for UL
            DoubleDisCo_Cfg_0l_RPV            = "Keras_Tensorflow_DoubleDisCo_Reg_0l_RPV_2016.cfg";           
            DoubleDisCo_Model_0l_RPV          = "keras_frozen_DoubleDisCo_Reg_0l_RPV_2016.pb";
            DoubleDisCo_Cfg_NonIsoMuon_0l_RPV = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_0l_RPV_2016.cfg";
            DoubleDisCo_Cfg_1l_RPV            = "Keras_Tensorflow_DoubleDisCo_Reg_1l_RPV_2016.cfg";
            DoubleDisCo_Model_1l_RPV          = "keras_frozen_DoubleDisCo_Reg_1l_RPV_2016.pb";
            DoubleDisCo_Cfg_NonIsoMuon_1l_RPV = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_1l_RPV_2016.cfg";

            // Use RPV config for now --- switch to SYY with dedicated training
            DoubleDisCo_Cfg_0l_SYY            = "Keras_Tensorflow_DoubleDisCo_Reg_0l_RPV_2016.cfg";           
            DoubleDisCo_Model_0l_SYY          = "keras_frozen_DoubleDisCo_Reg_0l_RPV_2016.pb";
            DoubleDisCo_Cfg_NonIsoMuon_0l_SYY = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_0l_RPV_2016.cfg";
            DoubleDisCo_Cfg_1l_SYY            = "Keras_Tensorflow_DoubleDisCo_Reg_1l_RPV_2016.cfg";
            DoubleDisCo_Model_1l_SYY          = "keras_frozen_DoubleDisCo_Reg_1l_RPV_2016.pb";
            DoubleDisCo_Cfg_NonIsoMuon_1l_SYY = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_1l_RPV_2016.cfg";

            //DoubleDisCo_Cfg_2l    = "Keras_Tensorflow_DoubleDisCo_Reg_2l_2016.cfg";
            //DoubleDisCo_Model_2l  = "keras_frozen_DoubleDisCo_Reg_2l_2016.pb";
            //DoubleDisCo_Cfg_NonIsoMuon_2l = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_2l_2016.cfg";
            leptonFileName        = "allInOne_leptonSF_UL.root";
            bjetFileName          = "allInOne_BTagEff_UL.root";
            bjetCSVFileName       = "wp_deepCSV_106XUL18_v2.csv";
            meanFileName          = "allInOne_SFMean_UL.root";
            blind                 = false;
            TopTaggerCfg          = "TopTaggerCfg_2018.cfg";
        }
        else if(filetag.find("2018") != std::string::npos) 
        {
            runYear               = "2018";
            Lumi                  = 59830.0;
            deepCSV_WP_loose      = 0.1208;
            deepCSV_WP_medium     = 0.4168;       
            deepCSV_WP_tight      = 0.7665;

            // Switch to proper year's training when done for UL
            DoubleDisCo_Cfg_0l_RPV            = "Keras_Tensorflow_DoubleDisCo_Reg_0l_RPV_2016.cfg";           
            DoubleDisCo_Model_0l_RPV          = "keras_frozen_DoubleDisCo_Reg_0l_RPV_2016.pb";
            DoubleDisCo_Cfg_NonIsoMuon_0l_RPV = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_0l_RPV_2016.cfg";
            DoubleDisCo_Cfg_1l_RPV            = "Keras_Tensorflow_DoubleDisCo_Reg_1l_RPV_2016.cfg";
            DoubleDisCo_Model_1l_RPV          = "keras_frozen_DoubleDisCo_Reg_1l_RPV_2016.pb";
            DoubleDisCo_Cfg_NonIsoMuon_1l_RPV = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_1l_RPV_2016.cfg";

            // Use RPV config for now --- switch to SYY with dedicated training
            DoubleDisCo_Cfg_0l_SYY            = "Keras_Tensorflow_DoubleDisCo_Reg_0l_RPV_2016.cfg";           
            DoubleDisCo_Model_0l_SYY          = "keras_frozen_DoubleDisCo_Reg_0l_RPV_2016.pb";
            DoubleDisCo_Cfg_NonIsoMuon_0l_SYY = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_0l_RPV_2016.cfg";
            DoubleDisCo_Cfg_1l_SYY            = "Keras_Tensorflow_DoubleDisCo_Reg_1l_RPV_2016.cfg";
            DoubleDisCo_Model_1l_SYY          = "keras_frozen_DoubleDisCo_Reg_1l_RPV_2016.pb";
            DoubleDisCo_Cfg_NonIsoMuon_1l_SYY = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_1l_RPV_2016.cfg";

            //DoubleDisCo_Cfg_2l    = "Keras_Tensorflow_DoubleDisCo_Reg_2l_2016.cfg";
            //DoubleDisCo_Model_2l  = "keras_frozen_DoubleDisCo_Reg_2l_2016.pb";
            //DoubleDisCo_Cfg_NonIsoMuon_2l = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_2l_2016.cfg";
            leptonFileName        = "allInOne_leptonSF_UL.root";
            bjetFileName          = "allInOne_BTagEff_UL.root";
            bjetCSVFileName       = "wp_deepCSV_106XUL18_v2.csv";
            meanFileName          = "allInOne_SFMean_UL.root";
            blind                 = false;
            TopTaggerCfg          = "TopTaggerCfg_2018.cfg";
        }

        tr.registerDerivedVar("runYear",                           runYear                          );
        tr.registerDerivedVar("Lumi",                              Lumi                             );
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
        tr.registerDerivedVar("DoubleDisCo_Cfg_0l_SYY",            DoubleDisCo_Cfg_0l_SYY           );
        tr.registerDerivedVar("DoubleDisCo_Model_0l_SYY",          DoubleDisCo_Model_0l_SYY         );
        tr.registerDerivedVar("DoubleDisCo_Cfg_NonIsoMuon_0l_SYY", DoubleDisCo_Cfg_NonIsoMuon_0l_SYY);
        tr.registerDerivedVar("DoubleDisCo_Cfg_1l_SYY",            DoubleDisCo_Cfg_1l_SYY           );
        tr.registerDerivedVar("DoubleDisCo_Model_1l_SYY",          DoubleDisCo_Model_1l_SYY         );
        tr.registerDerivedVar("DoubleDisCo_Cfg_NonIsoMuon_1l_SYY", DoubleDisCo_Cfg_NonIsoMuon_1l_SYY);

        //tr.registerDerivedVar("DoubleDisCo_Cfg_2l",    DoubleDisCo_Cfg_2l   );
        //tr.registerDerivedVar("DoubleDisCo_Model_2l",  DoubleDisCo_Model_2l );
        //tr.registerDerivedVar("DoubleDisCo_Cfg_NonIsoMuon_2l",    DoubleDisCo_Cfg_NonIsoMuon_2l   );
        tr.registerDerivedVar("leptonFileName",        leptonFileName       );        
        tr.registerDerivedVar("bjetFileName",          bjetFileName         );        
        tr.registerDerivedVar("bjetCSVFileName",       bjetCSVFileName      );        
        tr.registerDerivedVar("meanFileName",          meanFileName         );        
        tr.registerDerivedVar("etaCut",                2.4                  ); 
        tr.registerDerivedVar("blind",                 blind                );
        tr.registerDerivedVar("TopTaggerCfg",          TopTaggerCfg         );

        // Register Modules that are needed for each Analyzer
        if(analyzer=="MakeNJetDists") // for legacy 1l
        {
            const std::vector<std::string> modulesList = {
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Photon",
                "Jet",
                "BJet",
                "RunTopTagger",
                "CommonVariables",
                "FatJetCombine",
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
            const std::vector<std::string> modulesList = {
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Photon",
                "Jet",
                "BJet",
                "RunTopTagger",
                "CommonVariables",
                "Baseline",
                "FatJetCombine",
                "MakeMVAVariables_0l",
                "MakeMVAVariables_1l",
                "MakeMVAVariables_NonIsoMuon_0l",
                "MakeMVAVariables_NonIsoMuon_1l",
                "StopJets",
                "MakeStopHemispheres_OldSeed",
                "MakeStopHemispheres_OldSeed_NonIsoMuon",
                "BTagCorrector",
                "ScaleFactors",
                "StopGenMatch",
                "DoubleDisCo_0l_RPV",
                "DoubleDisCo_1l_RPV",
                "DoubleDisCo_NonIsoMuon_0l_RPV",
                "DoubleDisCo_NonIsoMuon_1l_RPV",
                "DoubleDisCo_0l_SYY",
                "DoubleDisCo_1l_SYY",
                "DoubleDisCo_NonIsoMuon_0l_SYY",
                "DoubleDisCo_NonIsoMuon_1l_SYY",
            };
            registerModules(tr, std::move(modulesList));
        }
        else if(analyzer=="MakeNNVariables")
        {
            const std::vector<std::string> modulesList = {
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Photon",
                "Jet",
                "BJet",
                "RunTopTagger",
                "CommonVariables",
                "Baseline",
                "FatJetCombine",
                "MakeMVAVariables_0l_old",
                //"MakeMVAVariables_0l",
                //"MakeMVAVariables_1l",
                //"MakeMVAVariables_2l",
                //"MakeMVAVariables_NonIsoMuon_0l",
                //"MakeMVAVariables_NonIsoMuon_1l",
                "StopJets",
                "MakeStopHemispheres_OldSeed",
                //"MakeStopHemispheres_OldSeed_NonIsoMuon",
                "MakeStopHemispheres_TopSeed",
                "BTagCorrector",
                "ScaleFactors",
                "StopGenMatch",
            };
            registerModules(tr, std::move(modulesList));
        }

        else if(analyzer=="AnalyzeLepTrigger" || analyzer=="HadTriggers_Analyzer" || analyzer=="CalculateBTagSF" || analyzer=="CalculateSFMean")
        {
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
                "ScaleFactors"
            };
            registerModules(tr, std::move(modulesList));
        }
        else if(analyzer=="Semra_Analyzer" || analyzer=="TopTagger_Analyzer" || analyzer=="AnalyzeTopTagger")
        {   
            const std::vector<std::string> modulesList = {
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Photon",
                "Jet",
                "BJet",
                "RunTopTagger",
                "CommonVariables",
                "FatJetCombine",
                "Baseline",
                "BTagCorrector",
                "ScaleFactors",
                "MakeMVAVariables_0l",
                "MakeMVAVariables_1l",
                "MakeMVAVariables_2l",
                "MakeMVAVariables_NonIsoMuon_0l",
                "MakeMVAVariables_NonIsoMuon_1l",
                "StopJets",
                "MakeStopHemispheres_OldSeed",
                "MakeStopHemispheres_OldSeed_NonIsoMuon",
                "MakeStopHemispheres_TopSeed",
                "StopGenMatch",
                //"DoubleDisCo_0l_RPV",
                //"DoubleDisCo_1l_RPV",
                //"DoubleDisCo_NonIsoMuon_0l_RPV",
                //"DoubleDisCo_NonIsoMuon_1l_RPV",
            };
            registerModules(tr, std::move(modulesList));
        } 
        else if(analyzer=="StealthHemispheres" || analyzer=="ISRJets_Analyzer")
        {
            const std::vector<std::string> modulesList = {
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Photon",
                "Jet",
                "BJet",
                "RunTopTagger",
                "CommonVariables",
                "FatJetCombine",
                "Baseline",
                "MakeMVAVariables_0l",
                //"MakeMVAVariables_1l",
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
        else if(analyzer=="TwoLepAnalyzer" || analyzer=="Make2LInputTrees")
        {   
            const std::vector<std::string> modulesList = {
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Photon",
                "Jet",
                "BJet",
                "CommonVariables",
                "FatJetCombine",
                "Baseline",
                "MakeMVAVariables_2l",
                "StopGenMatch",
                "BTagCorrector",
                "ScaleFactors",
                "TrainingNTupleVars",
            };
            registerModules(tr, std::move(modulesList));
        }
        else if(analyzer=="AnalyzeXsec")
        {   
            const std::vector<std::string> modulesList = {
                "PrepNTupleVars"
            };
            registerModules(tr, std::move(modulesList));
        }

        else
        {
            const std::vector<std::string> modulesList = {
                "PrepNTupleVars",
                "Muon",
                "Electron",
                "Photon",
                "Jet",
                "BJet",
                "RunTopTagger",
                "CommonVariables",
                "FatJetCombine",
                "Baseline",
                "BTagCorrector",
                "ScaleFactors"
            };
            registerModules(tr, std::move(modulesList));
        }
    }
};

#endif
