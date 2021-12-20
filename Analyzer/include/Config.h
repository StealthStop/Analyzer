#ifndef Confg_h
#define Confg_h

#include "SusyAnaTools/Tools/NTupleReader.h"

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
        const auto& DeepESMCfg                        = tr.getVar<std::string>("DeepESMCfg"                       );
        const auto& DeepESMModel                      = tr.getVar<std::string>("DeepESMModel"                     );
        const auto& DeepESMCfg_NonIsoMuon             = tr.getVar<std::string>("DeepESMCfg_NonIsoMuon"            );
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
            else if(module=="DeepEventShape")                        tr.emplaceModule<DeepEventShape>(DeepESMCfg, DeepESMModel);
            else if(module=="DeepEventShape_NonIsoMuon")             tr.emplaceModule<DeepEventShape>(DeepESMCfg_NonIsoMuon, DeepESMModel);
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

        std::string runYear, DeepESMCfg_NonIsoMuon, DeepESMCfg, DeepESMModel;
        std::string DoubleDisCo_Cfg_0l_RPV, DoubleDisCo_Model_0l_RPV, DoubleDisCo_Cfg_NonIsoMuon_0l_RPV; 
        std::string DoubleDisCo_Cfg_1l_RPV, DoubleDisCo_Model_1l_RPV, DoubleDisCo_Cfg_NonIsoMuon_1l_RPV; 
        std::string DoubleDisCo_Cfg_0l_SYY, DoubleDisCo_Model_0l_SYY, DoubleDisCo_Cfg_NonIsoMuon_0l_SYY; 
        std::string DoubleDisCo_Cfg_1l_SYY, DoubleDisCo_Model_1l_SYY, DoubleDisCo_Cfg_NonIsoMuon_1l_SYY; 
        //std::string DoubleDisCo_Cfg_2l, DoubleDisCo_Model_2l, DoubleDisCo_Cfg_NonIsoMuon_2l;      
        std::string leptonFileName, bjetFileName, bjetCSVFileName, meanFileName, TopTaggerCfg;
 
        double Lumi=0.0, deepCSV_WP_loose=0.0, deepCSV_WP_medium=0.0, deepCSV_WP_tight=0.0;
        bool blind = true;

        if(filetag.find("2016") != std::string::npos)
        {
            runYear                           = "2016";
            Lumi                              = 36330.0;
            deepCSV_WP_loose                  = 0.1918;
            deepCSV_WP_medium                 = 0.5847;
            deepCSV_WP_tight                  = 0.8767;            
            DeepESMCfg                        = "DeepEventShape_2016.cfg";
            DeepESMModel                      = "keras_frozen_2016.pb";
            DeepESMCfg_NonIsoMuon             = "DeepEventShape_NonIsoMuon_2016.cfg";
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
            leptonFileName        = "allInOne_leptonSF_2016.root";
            bjetFileName          = "allInOne_BTagEff.root";
            bjetCSVFileName       = "DeepCSV_2016LegacySF_WP_V1.csv";
            meanFileName          = "allInOne_SFMean.root";
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
            DeepESMCfg            = "DeepEventShape_2017.cfg";
            DeepESMModel          = "keras_frozen_2017.pb";
            DeepESMCfg_NonIsoMuon = "DeepEventShape_NonIsoMuon_2017.cfg";
            //DoubleDisCo_Cfg_0l    = "Keras_Tensorflow_DoubleDisCo_Reg_0l_2017.cfg";    
            //DoubleDisCo_Model_0l  = "keras_frozen_DoubleDisCo_Reg_0l_2017.pb";
            //DoubleDisCo_Cfg_NonIsoMuon_0l = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_0l_2017.cfg";
            //DoubleDisCo_Cfg_1l    = "Keras_Tensorflow_DoubleDisCo_Reg_1l_2017.cfg";
            //DoubleDisCo_Model_1l  = "keras_frozen_DoubleDisCo_Reg_1l_2017.pb";
            //DoubleDisCo_Cfg_NonIsoMuon_1l = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_1l_2017.cfg";
            //DoubleDisCo_Cfg_2l    = "Keras_Tensorflow_DoubleDisCo_Reg_2l_2017.cfg";
            //DoubleDisCo_Model_2l  = "keras_frozen_DoubleDisCo_Reg_2l_2017.pb";
            //DoubleDisCo_Cfg_NonIsoMuon_2l = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_2l_2017.cfg";
            leptonFileName        = "allInOne_leptonSF_2017.root";
            bjetFileName          = "allInOne_BTagEff.root";
            bjetCSVFileName       = "DeepCSV_94XSF_WP_V4_B_F.csv";
            meanFileName          = "allInOne_SFMean.root";
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
            DeepESMCfg            = "DeepEventShape_2018pre.cfg";
            DeepESMModel          = "keras_frozen_2018pre.pb";
            DeepESMCfg_NonIsoMuon = "DeepEventShape_NonIsoMuon_2018pre.cfg";
            //DoubleDisCo_Cfg_0l    = "Keras_Tensorflow_DoubleDisCo_Reg_0l_2018pre.cfg";    
            //DoubleDisCo_Model_0l  = "keras_frozen_DoubleDisCo_Reg_0l_2018pre.pb";
            //DoubleDisCo_Cfg_NonIsoMuon_0l = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_0l_2018pre.cfg";
            //DoubleDisCo_Cfg_1l    = "Keras_Tensorflow_DoubleDisCo_Reg_1l_2018pre.cfg";
            //DoubleDisCo_Model_1l  = "keras_frozen_DoubleDisCo_Reg_1l_2018pre.pb";
            //DoubleDisCo_Cfg_NonIsoMuon_1l = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_1l_2018pre.cfg";
            //DoubleDisCo_Cfg_2l    = "Keras_Tensorflow_DoubleDisCo_Reg_2l_2018pre.cfg";
            //DoubleDisCo_Model_2l  = "keras_frozen_DoubleDisCo_Reg_2l_2018pre.pb";
            //DoubleDisCo_Cfg_NonIsoMuon_2l = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_2l_2018pre.cfg";
            leptonFileName        = "allInOne_leptonSF_2018.root";
            bjetFileName          = "allInOne_BTagEff.root";
            bjetCSVFileName       = "DeepCSV_102XSF_WP_V1.csv";
            meanFileName          = "allInOne_SFMean.root";
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
            DeepESMCfg            = "DeepEventShape_2018post.cfg";
            DeepESMModel          = "keras_frozen_2018post.pb";
            DeepESMCfg_NonIsoMuon = "DeepEventShape_NonIsoMuon_2018post.cfg";
            //DoubleDisCo_Cfg_0l    = "Keras_Tensorflow_DoubleDisCo_Reg_0l_2018post.cfg";    
            //DoubleDisCo_Model_0l  = "eras_frozen_DoubleDisCo_Reg_0l_2018post.pb";
            //DoubleDisCo_Cfg_NonIsoMuon_0l = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_0l_2018post.cfg";
            //DoubleDisCo_Cfg_1l    = "Keras_Tensorflow_DoubleDisCo_Reg_1l_2018post.cfg";
            //DoubleDisCo_Model_1l  = "keras_frozen_DoubleDisCo_Reg_1l_2018pre.pb";
            //DoubleDisCo_Cfg_NonIsoMuon_1l = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_1l_2018post.cfg";
            //DoubleDisCo_Cfg_2l    = "Keras_Tensorflow_DoubleDisCo_Reg_2l_2018post.cfg";
            //DoubleDisCo_Model_2l  = "keras_frozen_DoubleDisCo_Reg_2l_2018post.pb";
            //DoubleDisCo_Cfg_NonIsoMuon_2l = "Keras_Tensorflow_NonIsoMuon_DoubleDisCo_Reg_2l_2018post.cfg";
            leptonFileName        = "allInOne_leptonSF_2018.root";
            bjetFileName          = "allInOne_BTagEff.root";
            bjetCSVFileName       = "DeepCSV_102XSF_WP_V1.csv";
            meanFileName          = "allInOne_SFMean.root";
            blind                 = false;
            TopTaggerCfg          = "TopTaggerCfg_2018.cfg";
        }

        tr.registerDerivedVar("runYear",                           runYear              );
        tr.registerDerivedVar("Lumi",                              Lumi                 );
        tr.registerDerivedVar("deepCSV_WP_loose",                  deepCSV_WP_loose     );
        tr.registerDerivedVar("deepCSV_WP_medium",                 deepCSV_WP_medium    );
        tr.registerDerivedVar("deepCSV_WP_tight",                  deepCSV_WP_tight     );
        tr.registerDerivedVar("isSignal",                          isSignal             );
        tr.registerDerivedVar("DeepESMCfg",                        DeepESMCfg           );
        tr.registerDerivedVar("DeepESMCfg_NonIsoMuon",             DeepESMCfg_NonIsoMuon);
        tr.registerDerivedVar("DeepESMModel",                      DeepESMModel         );        
        tr.registerDerivedVar("DoubleDisCo_Cfg_0l_RPV",            DoubleDisCo_Cfg_0l_RPV   );
        tr.registerDerivedVar("DoubleDisCo_Model_0l_RPV",          DoubleDisCo_Model_0l_RPV );
        tr.registerDerivedVar("DoubleDisCo_Cfg_NonIsoMuon_0l_RPV", DoubleDisCo_Cfg_NonIsoMuon_0l_RPV   );
        tr.registerDerivedVar("DoubleDisCo_Cfg_1l_RPV",            DoubleDisCo_Cfg_1l_RPV   );
        tr.registerDerivedVar("DoubleDisCo_Model_1l_RPV",          DoubleDisCo_Model_1l_RPV );
        tr.registerDerivedVar("DoubleDisCo_Cfg_NonIsoMuon_1l_RPV", DoubleDisCo_Cfg_NonIsoMuon_1l_RPV   );
        tr.registerDerivedVar("DoubleDisCo_Cfg_0l_SYY",            DoubleDisCo_Cfg_0l_SYY   );
        tr.registerDerivedVar("DoubleDisCo_Model_0l_SYY",          DoubleDisCo_Model_0l_SYY );
        tr.registerDerivedVar("DoubleDisCo_Cfg_NonIsoMuon_0l_SYY", DoubleDisCo_Cfg_NonIsoMuon_0l_SYY   );
        tr.registerDerivedVar("DoubleDisCo_Cfg_1l_SYY",            DoubleDisCo_Cfg_1l_SYY   );
        tr.registerDerivedVar("DoubleDisCo_Model_1l_SYY",          DoubleDisCo_Model_1l_SYY );
        tr.registerDerivedVar("DoubleDisCo_Cfg_NonIsoMuon_1l_SYY", DoubleDisCo_Cfg_NonIsoMuon_1l_SYY   );

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
                "MakeMVAVariables_0l",
                "MakeMVAVariables_1l",
                "MakeMVAVariables_2l",
                "MakeMVAVariables_NonIsoMuon_0l",
                "MakeMVAVariables_NonIsoMuon_1l",
                "StopJets",
                "MakeStopHemispheres_OldSeed",
                "MakeStopHemispheres_OldSeed_NonIsoMuon",
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
                "Baseline",
                "BTagCorrector",
                "ScaleFactors"
            };
            registerModules(tr, std::move(modulesList));
        }
        else if(analyzer=="Semra_Analyzer" || analyzer=="AnalyzeTopTagger")
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
                "DoubleDisCo_0l_RPV",
                "DoubleDisCo_1l_RPV",
                "DoubleDisCo_NonIsoMuon_0l_RPV",
                "DoubleDisCo_NonIsoMuon_1l_RPV",
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
        else if(analyzer=="Analyze1Lep")
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
                "MakeMVAVariables",
                "MakeMVAVariables_NonIsoMuon",
                "RunTopTagger",
                "Baseline",
                "DeepEventShape",
                "DeepEventShape_NonIsoMuon",
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
