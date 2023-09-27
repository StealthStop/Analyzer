#define MakeTopTagSFTree_cxx

#include "Analyzer/Analyzer/include/MakeTopTagSFTree.h"
#include "NTupleReader/include/NTupleReader.h"
#include "Framework/Framework/include/Utility.h"
#include "Framework/Framework/include/MiniTupleMaker.h"

#include "Framework/Framework/include/Muon.h"
#include "Framework/Framework/include/Electron.h"
#include "Framework/Framework/include/Photon.h"
#include "Framework/Framework/include/Jet.h"
#include "Framework/Framework/include/BJet.h"
#include "Framework/Framework/include/CommonVariables.h"
#include "Framework/Framework/include/RunTopTagger.h"
#include "Framework/Framework/include/Baseline.h"
#include "Framework/Framework/include/BTagCorrector.h"
#include "Framework/Framework/include/ScaleFactors.h"

#include <iostream>

MakeTopTagSFTree::MakeTopTagSFTree()
{
    eventCounter = std::make_shared<TH1D>("EventCounter", "EventCounter", 2, -1.1, 1.1);
    jecvars  = {"", "JECup", "JECdown", "JERup", "JERdown"};

    for (const auto& jecvar : jecvars)
    {
        treeInit[jecvar] = false;
    }
}

void MakeTopTagSFTree::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    const auto& filetag            = tr.getVar<std::string>("filetag"          );
    const auto& runtype            = tr.getVar<std::string>("runtype"          );
    const auto& runYear            = tr.getVar<std::string>("runYear"          );
    const auto& btagEffFileName    = tr.getVar<std::string>("btagEffFileName"  );
    const auto& bjetTagFileName    = tr.getVar<std::string>("bjetTagFileName"  );
    const auto& leptonFileName     = tr.getVar<std::string>("leptonFileName"   );
    const auto& hadronicFileName   = tr.getVar<std::string>("hadronicFileName" );
    const auto& toptaggerFileName  = tr.getVar<std::string>("toptaggerFileName");
    const auto& meanFileName       = tr.getVar<std::string>("meanFileName"     );
    const auto& TopTaggerCfg       = tr.getVar<std::string>("TopTaggerCfg"     );

    for(const auto& jecvar : jecvars)
    {
        // Cannot do JEC and JER variations for data
        if (runtype == "Data" or jecvar == "")
            continue;

        Jet             jet(jecvar);
        BJet            bjet(jecvar);
        Muon            muon(jecvar);
        Photon          photon(jecvar);
        Electron        electron(jecvar);
        Baseline        baseline(jecvar);
        RunTopTagger    topTagger(TopTaggerCfg, jecvar);
        CommonVariables commonVariables(jecvar);

        // Remember, order matters here !
        // Follow what is done in Config.h
        tr.registerFunction(muon);
        tr.registerFunction(electron);
        tr.registerFunction(photon);
        tr.registerFunction(jet);
        tr.registerFunction(bjet);
        tr.registerFunction(commonVariables);
        tr.registerFunction(topTagger);
        tr.registerFunction(baseline);

        if (runtype == "MC")
        {
            ScaleFactors  scaleFactors(runYear, leptonFileName, hadronicFileName, toptaggerFileName, meanFileName, filetag, jecvar);
            BTagCorrector bTagCorrector(btagEffFileName, "", bjetTagFileName, "", filetag);
            bTagCorrector.SetVarNames("GenParticles_PdgId", "Jets"+jecvar, "GoodJets_pt30"+jecvar, "Jets"+jecvar+"_bJetTagDeepFlavourtotb", "Jets"+jecvar+"_partonFlavor", jecvar);

            tr.registerFunction(bTagCorrector);
            tr.registerFunction(scaleFactors);
        }
    }

    while( tr.getNextEvent() )
    {
        if( maxevents != -1 && tr.getEvtNum() > maxevents )
            break;
        if( tr.getEvtNum() % 1000 == 0 )
            printf( " Event %i\n", tr.getEvtNum() );

        const auto& evtCounter = tr.getVar<int>("eventCounter");
        eventCounter->Fill(evtCounter);

        for(const auto& jecvar : jecvars)
        {
            // Cannot do JEC and JER variations for data
            if (jecvar != "" and runtype == "Data")
                continue;

            // If an event is not interesting for any channel selection (see explanation and definition
            // of lostCauseEvent boolean in Baseline.h), then move on here and do not waste any more
            // time on looping over histos or getting vars.
            // This only matters if the user specifies -s on the command line
            // N.B. We already filled the EventCounter histogram for our due-diligence of counting up
            // every single event.
            const auto& lostCauseEvent = tr.getVar<bool>("lostCauseEvent" + jecvar);
            const auto& fastMode       = tr.getVar<bool>("fastMode");

            if (lostCauseEvent and fastMode)
                continue;

            // Event-level booleans (triggers, filters, etc.)
            const auto& pass_METFilters  = tr.getVar<bool>("passMETFilters");
            const auto& pass_MadHT       = tr.getVar<bool>("passMadHT");
            const auto& pass_MuonTrigger = tr.getVar<bool>("passTriggerMuon");
            const auto& pass_QCDTrigger  = tr.getVar<bool>("passTriggerQCD");
            const auto& pass_HEMveto     = tr.getVar<bool>("passElectronHEMveto");
            const auto& pass_JetID       = tr.getVar<bool>("JetID");        

            bool pass_DataQuality = pass_METFilters && pass_HEMveto && pass_JetID;

            // Lepton quantities
            const auto& Muons          = tr.getVec<utility::LorentzVector>("Muons");
            const auto& GoodMuons      = tr.getVec<bool>("GoodMuons" + jecvar);
            const auto& NGoodMuons     = tr.getVar<int>("NGoodMuons" + jecvar);
            const auto& NNonIsoMuons   = tr.getVar<int>("NNonIsoMuons" + jecvar);
            const auto& NGoodElectrons = tr.getVar<int>("NGoodElectrons" + jecvar);

            // Jet quantities
            const auto& Jets                  = tr.getVec<utility::LorentzVector>("Jets");
            const auto& NGoodJets_pt30        = tr.getVar<int>("NGoodJets_pt30" + jecvar);
            const auto& NGoodBJets_pt30       = tr.getVar<int>("NGoodBJets_pt30" + jecvar);
            const auto& GoodBJets_pt30_loose  = tr.getVec<bool>("GoodBJets_pt30_loose" + jecvar);
            const auto& GoodBJets_pt30        = tr.getVec<bool>("GoodBJets_pt30" + jecvar);
            const auto& NGoodBJets_pt30_loose = tr.getVar<int>("NBJets_pt30_loose" + jecvar);

            // Event-level quantities
            const auto& HT_trigger_pt30 = tr.getVar<double>("HT_trigger_pt30" + jecvar);
            const auto& METPhi          = tr.getVar<float>("METPhi");
            const auto& MET             = tr.getVar<float>("MET");

            // Get the 4-vec for the MET
            utility::LorentzVector lvMET;
            lvMET.SetPt(MET); lvMET.SetEta(0.0); lvMET.SetPhi(METPhi); lvMET.SetE(MET);

            bool pass_ZeroMuons     = NGoodMuons      == 0;
            bool pass_ZeroElectrons = NGoodElectrons  == 0;

            bool pass_SingleMuon    = NGoodMuons      == 1;
            bool pass_ExtraLepVeto  = NGoodElectrons  == 0 and NNonIsoMuons == 0;
            bool pass_ge4jets       = NGoodJets_pt30  >= 4;
            bool pass_1bjet         = NGoodBJets_pt30 >= 1;

            // The medium b jet (above) also passes the loose working point
            // So if we want an additional b jet that is loose we will have at least 2 loose b
            bool pass_1bjetLoose    = NGoodBJets_pt30_loose >= 2;
            bool pass_HT200         = HT_trigger_pt30        > 200.0;
            bool pass_HT1000        = HT_trigger_pt30        > 1000.0;
            bool pass_HT1500        = HT_trigger_pt30        > 1500.0;

            bool pass_HT = runYear.find("2016") != std::string::npos ? pass_HT1000 : pass_HT1500;

            bool pass_MuonBjetdR     = false; double MuonBjetdR   = 9999.0;
            bool pass_MuonBjetMass   = false; double MuonBjetMass = -1.0;
            bool pass_MuonBjetdRMass = false;
            for (unsigned int iJet = 0; iJet < Jets.size(); iJet++)
            {
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

            tr.registerDerivedVar("MuonBjetdR" + jecvar,   MuonBjetdR);
            tr.registerDerivedVar("MuonBjetMass" + jecvar, MuonBjetMass);
            tr.registerDerivedVar("MuonMETdPhi" + jecvar,  MuonMETdPhi);
            tr.registerDerivedVar("MuonMETmT" + jecvar,    MuonMETmT);

            bool pass_toSaveTTCR = pass_DataQuality                                         &&
                                   pass_MadHT                                               &&
                                   pass_MuonTrigger                                         &&
                                   pass_SingleMuon                                          &&
                                   pass_ExtraLepVeto                                        &&
                                   pass_ge4jets                                             &&
                                   pass_1bjet                                               &&
                                   pass_1bjetLoose                                          &&
                                   pass_HT200                                                ; 

            // To be used with Single-Muon-triggered data and measuring tagging efficiency SF
            bool pass_SemiLepTTbarCR = pass_toSaveTTCR      && 
                                       pass_MuonBjetdRMass  &&
                                       pass_MuonMETdPhi     &&
                                       pass_MuonMETmT       ;


            bool pass_toSaveQCDCR = pass_DataQuality   && 
                                    pass_QCDTrigger    &&
                                    pass_MadHT         &&
                                    pass_ge4jets       &&
                                    pass_ZeroMuons     &&
                                    pass_ZeroElectrons &&
                                    pass_HT             ;

            // For use with JetHT-triggered data and measuring mistag SF
            bool pass_QCDCR = pass_toSaveQCDCR;

            tr.registerDerivedVar("pass_TTCR" + jecvar,  pass_SemiLepTTbarCR);
            tr.registerDerivedVar("pass_QCDCR" + jecvar, pass_QCDCR);

            tr.registerDerivedVar("pass_preTTCR" + jecvar,  pass_toSaveTTCR);
            tr.registerDerivedVar("pass_preQCDCR" + jecvar, pass_toSaveQCDCR);

            double weightQCD   = 1.0;
            double weightTTbar = 1.0;
            if(runtype == "MC")
            {
                // Define Lumi weight
                const auto& Weight = tr.getVar<float>("Weight");
                const auto& lumi   = tr.getVar<double>("FinalLumi");
                double eventweight = lumi * Weight;

                double puScaleFactor        = tr.getVar<double>("puWeightCorr");
                double topPtScaleFactor     = tr.getVar<double>("topPtScaleFactor");
                double prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor");
                double bTagScaleFactor      = tr.getVar<double>("bTagSF_EventWeightSimple_Central" + jecvar);
                double totGoodMuonSF        = tr.getVar<double>("totGoodMuonSF" + jecvar);

                weightTTbar *= eventweight * puScaleFactor * prefiringScaleFactor * topPtScaleFactor * bTagScaleFactor * totGoodMuonSF;
                weightQCD   *= eventweight * puScaleFactor * prefiringScaleFactor * topPtScaleFactor;
            }

            tr.registerDerivedVar("weightTTbar" + jecvar, weightTTbar);
            tr.registerDerivedVar("weightQCD" + jecvar,   weightQCD);

            if( !treeInit[jecvar] ) {
                //-----------------------------------
                //  Initialize the tree
                //-----------------------------------       
                std::set<std::string> variables = {
                    "bestTopPt" + jecvar,
                    "bestTopEta" + jecvar,
                    "bestTopPhi" + jecvar,
                    "bestTopMass" + jecvar,
                    "bestTopDisc" + jecvar,
                    "bestTopNconst" + jecvar,
                    "bestTopMassGenMatch" + jecvar,
                    "bestRTopPt" + jecvar,
                    "bestRTopEta" + jecvar,
                    "bestRTopPhi" + jecvar,
                    "bestRTopMass" + jecvar,
                    "bestRTopDisc" + jecvar,
                    "bestRTopMassGenMatch" + jecvar,
                    "bestMTopPt" + jecvar,
                    "bestMTopEta" + jecvar,
                    "bestMTopPhi" + jecvar,
                    "bestMTopMass" + jecvar,
                    "bestMTopDisc" + jecvar,
                    "bestMTopMassGenMatch" + jecvar,
                    "NGoodJets_pt30" + jecvar,
                    "NGoodBJets_pt30" + jecvar,
                    "NGoodBJets_pt30_loose" + jecvar,
                    "HT_trigger_pt30" + jecvar,
                    "MuonBjetdR" + jecvar,
                    "MuonBjetMass" + jecvar,
                    "MuonMETdPhi" + jecvar,   
                    "MuonMETmT" + jecvar,     
                    "pass_TTCR" + jecvar, 
                    "pass_QCDCR" + jecvar,
                    "pass_preTTCR" + jecvar, 
                    "pass_preQCDCR" + jecvar,
                };

                if (runtype == "MC")
                {
                    variables.insert("weightTTbar" + jecvar);
                    variables.insert("weightQCD" + jecvar);
                    variables.insert("totGoodMuonSF" + jecvar);
                    variables.insert("totGoodMuonSF_Up" + jecvar);
                    variables.insert("totGoodMuonSF_Down" + jecvar);
                    variables.insert("bTagSF_EventWeightSimple_Central" + jecvar);
                    variables.insert("bTagSF_EventWeightSimple_Up" + jecvar);
                    variables.insert("bTagSF_EventWeightSimple_Down" + jecvar);

                    if (jecvar == "")
                    {
                        variables.insert("MET");
                        variables.insert("scaleWeightUp");
                        variables.insert("scaleWeightDown");
                        variables.insert("PSweight_ISRUp");
                        variables.insert("PSweight_ISRDown");
                        variables.insert("PSweight_FSRUp");
                        variables.insert("PSweight_FSRDown");
                        variables.insert("PDFweightUp");
                        variables.insert("PDFweightDown");
                        variables.insert("puWeightCorr");
                        variables.insert("puSysUpCorr");
                        variables.insert("puSysDownCorr");
                        variables.insert("prefiringScaleFactor");
                        variables.insert("prefiringScaleFactorUp");
                        variables.insert("prefiringScaleFactorDown");
                    }
                }   
                else if (runtype == "Data")
                {
                    variables.insert("Weight");
                }

                std::string myTreeName = "TopTagSFSkim" + jecvar;
                myTree[jecvar] = new TTree( (myTreeName).c_str() , (myTreeName).c_str() );
                myTopTagSFTuple[jecvar] = new MiniTupleMaker( myTree[jecvar] );
                myTopTagSFTuple[jecvar]->setTupleVars(variables);
                myTopTagSFTuple[jecvar]->initBranches(tr);

                treeInit[jecvar] = true;
            }

            // Requirements not on NonIsoMuon jets or HT as those will be more restrictive with a non iso muon present
            if( pass_toSaveQCDCR or pass_toSaveTTCR ) {
                myTopTagSFTuple[jecvar]->fill(tr);
            }
        }
    } 
}
      
void MakeTopTagSFTree::WriteHistos( TFile* outfile ) 
{
    outfile->cd();

    eventCounter->Write();

    for (auto& iter : myTree)
    {
        iter.second->Write();

        delete iter.second;    
        delete myTopTagSFTuple[iter.first];
    }
}
