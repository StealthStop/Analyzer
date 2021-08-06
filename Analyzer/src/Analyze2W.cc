#define Analyze2W_cxx

#include "Analyzer/Analyzer/include/Analyze2W.h"
#include "Framework/Framework/include/Utility.h"
#include "SusyAnaTools/Tools/NTupleReader.h"

#include <TCanvas.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <utility>
#include <vector>

constexpr int W_PDGID = 24;
constexpr int T_PDGID = 6;
constexpr int STOP_PDGID = 1000006;

Analyze2W::Analyze2W() { InitHistos(); }

auto createNewHistogram(auto &myhistos, std::string name, int v1, int v2,
                        int v3, std::vector<std::string> props) {
  name = std::accumulate(props.begin(), props.end(), name,
                         [](std::string &x, auto y) {
                           x += "_";
                           return x += y;
                         });
  const char *x = name.c_str();
  auto ret = std::make_shared<TH1D>(x, x, v1, v2, v3);
  myhistos.emplace(name, std::move(ret));
  return ret;
}

auto createNewHistogram(auto &myhistos, std::string name, int v1, int v2,
                        int v3) {
  return createNewHistogram(myhistos, name, v1, v2, v3, {});
}

// Define all your histograms here.
void Analyze2W::InitHistos() {
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  // This event counter histogram is necessary so that we know that all the
  // condor jobs ran successfully. If not, when you use the hadder script, you
  // will see a discrepancy in red as the files are being hadded.
  //    createNewHistogram(my_histos, "EventCounter", 2, -1.1, 1.1  ) ;

  createNewHistogram(my_histos, "h_njets_l0", 15, 0, 15);
  createNewHistogram(my_histos, "h_njets_l1", 15, 0, 15);
  createNewHistogram(my_histos, "h_njets_l2", 15, 0, 15);
  createNewHistogram(my_histos, "h_njets_lINF", 15, 0, 15);
  createNewHistogram(my_histos, "h_met", 200, 0, 2000);
  // Define 1D histograms
}

// Put everything you want to do per event here.
void Analyze2W::Loop(NTupleReader &tr, double, int maxevents, bool) {
  while (tr.getNextEvent()) {
      std::cout << "Analyzing event";
    // This is added to count the number of events- do not change the next two
    // lines.
    //        const auto& eventCounter        = tr.getVar<int>("eventCounter");
    //        my_histos["EventCounter"]->Fill( eventCounter );

    //--------------------------------------------------
    //-- Print Event Number
    //--------------------------------------------------

    if (maxevents != -1 && tr.getEvtNum() >= maxevents)
      break;
    if (tr.getEvtNum() & (10000 == 0))
      printf(" Event %i\n", tr.getEvtNum());

#define makeVar(name, type) const type &name = tr.getVar<type>(#name)
#define makeVec(name, type)                                                    \
  const std::vector<type> &name = tr.getVec<type>(#name)

    // clang-format off
    makeVar(runtype                  , std::string    );
    makeVar(NJets                     , int            );
    makeVar(JetID                    , bool           );
    makeVar(NGoodLeptons             , int            );
    /*
    makeVec(Jets                     , TLorentzVector );
    makeVar(MET                      , double         );
    makeVar(METPhi                   , double         );
    makeVec(GenElectrons, TLorentzVector);
    makeVec(GenMuons, TLorentzVector);
    makeVec(GenParticles, TLorentzVector );
    makeVec(GenParticles_Charge, int );
    makeVec(GenParticles_ParentId, int );
    makeVec(GenParticles_ParentIdx, int );
    makeVec(GenParticles_PdgId, int );
    makeVec(GenParticles_Status, int );
    */

     // makeVar(passMadHT                , bool           );
    // clang-format on

    // ------------------------
    // -- Define weight
    // ------------------------
    double weight = 1.0;
    double eventweight = 1.0;
    double leptonScaleFactor = 1.0;
    double bTagScaleFactor = 1.0;
    double htDerivedScaleFactor = 1.0;
    double prefiringScaleFactor = 1.0;
    double puScaleFactor = 1.0;

    if (runtype == "MC") {
    //  if (!passMadHT)
     //   continue; // Make sure not to double count DY events
      // Define Lumi weight
      const auto &Weight = tr.getVar<double>("Weight");
      const auto &lumi = tr.getVar<double>("Lumi");
      eventweight = lumi * Weight;

      // Define lepton weight

      // PileupWeight = tr.getVar<double>("_PUweightFactor");
      // bTagScaleFactor   =
      // tr.getVar<double>("bTagSF_EventWeightSimple_Central");
     //  htDerivedScaleFactor = tr.getVar<double>("htDerivedweight");
     //  prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor");
      // puScaleFactor = tr.getVar<double>("puWeightCorr");

      weight *= eventweight * leptonScaleFactor * bTagScaleFactor *
                htDerivedScaleFactor * prefiringScaleFactor * puScaleFactor;
    }
    switch (NGoodLeptons) {
    case 0:
      my_histos["h_njets_l0"]->Fill(NJets, weight);
      break;
    case 1:
      my_histos["h_njets_l1"]->Fill(NJets, weight);
      break;
    case 2:
      my_histos["h_njets_l2"]->Fill(NJets, weight);
      break;
    default:
      my_histos["h_njets_lINF"]->Fill(NJets, weight);
      break;
    }
  }
}

void Analyze2W::WriteHistos(TFile *outfile) {
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

















