#define Analyze2W_cxx

#define SHOW_PROGRESS 1
#define PER_EVENT_LOG 0
#define IMAGE_PRODUCTION 1

#if (SHOW_PROGRESS && PER_EVENT_LOG)
#error "Cannot have both per-event logging and a progress bar"
#endif

#include "Analyzer/Analyzer/include/Analyze2W.h"
#include "Framework/Framework/include/Utility.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <utility>
#include <vector>
#include "fastjet/ClusterSequence.hh"
#include "2W_Utils.hpp"

/*
   To Add:
   - Individual Jet Pts up to some number
   - Individual Lep Pts up to some number
   - N Bjets (look at diff working points)
   - Angle between leptons (2 Lep case)
   - Angle between W
   - Fraction of W above 200GeV pt

   - TTBar background
   - QCD  for 0 lepton
   - W and Jets
   - W and bosons

   2 Lep:
   - DY + Jets

   Basic cuts:
   - 4 jets

   Take a gander at W taggers

   Use WP: 0.918, check twiki for PT cut for jets, check for eta

   What does clean AK8 mean. What is getting cleaned out.

   Plot of HT.
   Check how produced: HT_trigger_pt30 (just sums the pt of Jets)
   Study the paper, take a look at the cuts they make, make plots.
   Take a look at the cutflow.

   N-1 plots. N cuts, apply all but one, look at distribution of missing one.

 */


class HTCut: public Cut{ 
    public:
        HTCut(): Cut("HTCut", false, {"HT<=900", "HT>900"}){}
     void compute(const SliceData& data)  const override{
            return data->Ht>900;
         }
};



const int W_PDGID = 24, E_PDGID = 11, M_PDGID = 13, T_PDGID = 6,
          STOP_PDGID = 1000006;
const int JETCOUNT = 10;
const int LEP_COUNT = 10;
const double W_WP = 0.918;
const double WTAG_PT = 200.0f;

static std::vector<std::string> cuts;

Analyze2W::Analyze2W() {
  cuts = getAllNames<std::vector<std::vector<std::string>>>(
      {{"l0", "l1", "l2", "lINF"},
       // {"njets<=4", "njets>4"},
       {"HT<=900", "HT>900"},
       {"LeadingJetsPass", "NotLeadingJetsPass"},
       {"SigPassed", "NotSigPassed"}}
       );
  InitHistos();
}

auto getHistName(int num_leps, int NJets, float HT,
                 const std::vector<TLorentzVector> &Jets,
                 const std::vector<double> nsr21, 
                 const std::vector<double> nsr43,
                const std::vector<double> nsr42
                 ) {
  auto lepsel = [num_leps]() {
    const std::map<int, std::string> vals = {
        {0, "l0"}, {1, "l1"}, {2, "l2"}, {99999, "lINF"}};
    const auto x =
        std::lower_bound(vals.begin(), vals.end(), num_leps,
                         [](auto &&x, int y) { return x.first < y; });
    return std::string("_") + x->second;
  };
  auto njetsel = [NJets]() {
    return std::string("_") + ((NJets > 4) ? "njets>4" : "njets<=4");
  };
  auto htsel = [HT]() {
    return std::string("_") + ((HT > 900) ? "HT>900" : "HT<=900");
  };
auto firstLargestPt =
    std::max_element(Jets.begin(), Jets.end(),
                     [](auto &x, auto &y) { return x.Pt() < y.Pt(); });
auto secondLargestPt = std::max_element(
    Jets.begin(), Jets.end(),
    [&firstLargestPt](auto &x, auto &y) {
      return (x != *firstLargestPt) && (x.Pt() < y.Pt());
    });
    std::size_t j1_index =  std::distance(Jets.begin(), firstLargestPt);
    std::size_t j2_index =  std::distance(Jets.begin(), secondLargestPt);
  auto bothjets = [&]() {
        if (std::size(Jets) < 2) {
          return std::string("_NotLeadingJetsPass");
        }
          auto jetpass = [](const TLorentzVector &j,double t21){
            return j.Pt() > 0.4 && std::abs(j.Eta()) > 2 && t21 < 0.75;
          };
        double largestTau21 = nsr21[j1_index];
        double secondTau21 = nsr21[j2_index];
        return (jetpass(*firstLargestPt, largestTau21) &&
                jetpass(*secondLargestPt, secondTau21))
                   ? std::string("_LeadingJetsPass")
                   : std::string("_NotLeadingJetsPass");
      };

  auto subjetiness = [&](){
  if(bothjets() == std::string("_LeadingJetsPass")){
    return (nsr43[j1_index] < 0.8 
           && nsr43[j2_index] <0.8 
           && nsr42[j1_index] < 0.5 
           && nsr42[j2_index] < 0.5
           && std::abs(firstLargestPt->Mag() - secondLargestPt->Mag())/(firstLargestPt->Mag() + secondLargestPt->Mag()) < 0.1
           && std::abs(firstLargestPt->Eta() - secondLargestPt->Eta()) < 1.0
            )? std::string("_SigPassed") : std::string("_NotSigPassed");
  } else {return std::string("_NotSigPassed");}
  };

  std::array<std::string, 4>
          temp{lepsel(), htsel(), bothjets(), subjetiness()};
  return combine(temp);
}

std::vector<double> computeRatio(std::vector<double> v1,
                                 const std::vector<double> &v2) {
  for (int i = 0; i < std::size(v1); ++i) {
    v1[i] /= v2[i];
  }
  return v1;
}

void Analyze2W::InitHistos() {
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  createNewHistograms(my_histos, "EventCounter", 2, -1.1, 1.1, cuts,
                      "Event Counter");
  for (int i = 0; i < JETCOUNT; ++i) {
    createNewHistograms(my_histos, std::string("JetPt_") + std::to_string(i),
                        50, 0, 500, cuts, "Jet PT (GoodJets_pt20)", "Events",
                        true, true);
  }
  createNewHistograms(my_histos, "njets", 15, 0, 15, cuts, "NJets");
  createNewHistograms(my_histos, "nAK8jets", 15, 0, 15, cuts, "NAK8Jets");

  createNewHistograms(my_histos, "nCA12Jets", 15, 0, 15, cuts, "NAK8Jets");

  createNewHistograms(my_histos, "nbjets", 5, 0, 5, cuts, "NBJets");
  createNewHistograms(my_histos, "nwjets", 5, 0, 5, cuts, "NWJets_Reco");
  createNewHistograms(my_histos, "MT2", 100, 0, 1000, cuts, "MT2_All", "Events",
                      true, true);
  createNewHistograms(my_histos, "met", 200, 0, 2000, cuts, "MET");
  createNewHistograms(my_histos, "WPhi", 100, -3.2, 3.2, cuts, "W Phi");
  createNewHistograms(my_histos, "WEta", 100, -3.2, 3.2, cuts, "W Eta");

  createNewHistograms(my_histos, "NSubJettiness4", 100, 0, 1, cuts, "NSJ 4");
  createNewHistograms(my_histos, "NSubJettiness3", 100, 0, 1, cuts, "NSJ 3");
  createNewHistograms(my_histos, "NSubJettiness2", 100, 0, 1, cuts, "NSJ 2");
  createNewHistograms(my_histos, "NSubJettiness1", 100, 0, 1, cuts, "NSJ 1");

  createNewHistograms(my_histos, "NSubRatio43", 100, 0, 2, cuts, "NSR 42");
  createNewHistograms(my_histos, "NSubRatio42", 100, 0, 2, cuts, "NSR 43");
  createNewHistograms(my_histos, "NSubRatio21", 100, 0, 2, cuts, "NSR 21");

  createNewHistograms(my_histos, "WE", 200, 0, 2000, cuts, "W Energy");
  createNewHistograms(my_histos, "WP", 200, 0, 2000, cuts, "W P");
  createNewHistograms(my_histos, "WPt", 200, 0, 2000, cuts, "W PT");
  createNewHistograms(my_histos, "JetPt", 50, 0, 500, cuts,
                      "jet pt (goodjets_pt20)", "", true, true);

  createNewHistograms(my_histos, "CA12Pt", 50, 0, 500, cuts,
                      "CA12Pt ", "Events", true, true);

  createNewHistograms(my_histos, "HT_pt30", 200, 0, 2000, cuts, "HT_pt30",
                      "Events", true, true);

  createNewHistograms(my_histos, "Gen_Lep_Angle", 100, 0, 4, cuts, "Angle");
  createNewHistograms(my_histos, "Gen_W_Angle", 100, 0, 4, cuts, "Angle");
}

void Analyze2W::Loop(NTupleReader &tr, double, int maxevents, bool) {
#if SHOW_PROGRESS
  ProgressBar prog(maxevents, "Processing Event ");
#endif
  while (tr.getNextEvent()) {
#if SHOW_PROGRESS
    (++prog).display();
#endif
#if PER_EVENT_LOG
    std::cout << "--------------------------------------------\n";
    std::cout << "Event: " << tr.getEvtNum() << "\n";
#endif
    if (maxevents != -1 && tr.getEvtNum() >= maxevents) {
      std::cout << std::endl;
      break;
    }
    const auto &eventCounter = tr.getVar<int>("eventCounter");

#if !SHOW_PROGRESS && !PER_EVENT_LOG
    if (tr.getEvtNum() & (10000 == 0))
      printf(" Event %i\n", tr.getEvtNum());
#endif

#define makeVec(name, t)                                                       \
  const std::vector<argument_type<void(t)>::type> &name =                      \
      tr.getVec<argument_type<void(t)>::type>(#name)

#define makeVarIdenticalC(name, type) type name = tr.getVar<type>(#name)
#define makeVarCustomC(name, type, myname) type myname = tr.getVar<type>(#name)
#define GET_MACRO(_1, _2, _3, NAME, ...) NAME
#define makeVarC(...)                                                          \
  GET_MACRO(__VA_ARGS__, makeVarCustomC, makeVarIdenticalC)(__VA_ARGS__)

#define makeVarIdentical(name, type) const type &name = tr.getVar<type>(#name)
#define makeVarCustom(name, type, myname)                                      \
  const type &myname = tr.getVar<type>(#name)
#define GET_MACRO(_1, _2, _3, NAME, ...) NAME
#define makeVar(...)                                                           \
  GET_MACRO(__VA_ARGS__, makeVarCustom, makeVarIdentical)(__VA_ARGS__)

    makeVar(runtype, std::string);
    makeVar(NJets, int);
    makeVec(Jets, TLorentzVector);
    makeVec(JetsAK8Clean, TLorentzVector);

    makeVec(JetsAK8, TLorentzVector);

    makeVec(GenJetsAK8, TLorentzVector);
    makeVec(GoodJets_pt20, bool);
    makeVar(NGoodBJets_pt30, int);

    makeVec(JetsCA12, TLorentzVector);
    makeVec(JetsCA12_NsubjettinessTau4, double);
    makeVec(JetsCA12_NsubjettinessTau3, double);
    makeVec(JetsCA12_NsubjettinessTau2, double);
    makeVec(JetsCA12_NsubjettinessTau1, double);

#define MAKE_RATIO(x,y) auto nsr##x##y = computeRatio(JetsCA12_NsubjettinessTau##x,JetsCA12_NsubjettinessTau##y)
    MAKE_RATIO(4,3);
    MAKE_RATIO(4,2);
    MAKE_RATIO(2,1);


    makeVar(MET, double);
    makeVar(MT2, double);
    makeVar(HT_trigger_pt30, double);

    makeVec(GenParticles, TLorentzVector);
    makeVec(GenParticles_PdgId, int);
    makeVec(GenParticles_ParentId, int);
    makeVec(GenParticles_Status, int);

    makeVec(JetsAK8_wDiscriminatorDeep, double);

    makeVec(GoodLeptons, (std::pair<std::string, TLorentzVector>));

#if PER_EVENT_LOG && 0
    std::cout << "NBJets: " << NGoodBJets_pt30 << "\n";
    std::cout << "Clean AK8 Jets: " << JetsAK8Clean << "\n";
    std::cout << "AK8 Jets: " << JetsAK8 << "\n";
    std::cout << "GenAK8 Jets: " << GenJetsAK8 << "\n";
    std::cout << "Jets: " << Jets;
#endif



/*
    int nAK8jets_wtaggable =
        std::count_if(JetsAK8.begin(), JetsAK8.end(),
                      [](const auto &v) { return v.Pt() >= WTAG_PT; });
    int nwjets = std::count_if(std::begin(JetsAK8_wDiscriminatorDeep),
                               std::end(JetsAK8_wDiscriminatorDeep),
                               [](double x) { return x >= W_WP; });
                               */

    // std::sort(GoodLeptons.begin(), GoodLeptons.end(), [](auto x,
    // auto y){return x.second.Pt() < y.second.Pt();});

    double weight = 1.0, eventweight = 1.0, leptonScaleFactor = 1.0,
           bTagScaleFactor = 1.0, htDerivedScaleFactor = 1.0,
           prefiringScaleFactor = 1.0, puScaleFactor = 1.0;

    if (runtype == "MC") {
      const auto &Weight = tr.getVar<double>("Weight");
      const auto &lumi = tr.getVar<double>("Lumi");
      eventweight = lumi * Weight;
      weight *= eventweight * leptonScaleFactor * bTagScaleFactor *
                htDerivedScaleFactor * prefiringScaleFactor * puScaleFactor;
    }

    std::vector<TLorentzVector> Ws, Leps;
    auto isHard = [](int x) {
      return x == 1 || (std::abs(x) < 60 && std::abs(x) > 20);
    };
    for (std::size_t i = 0; i < GenParticles.size(); ++i) { // clang-format off
      switch (GenParticles_PdgId[i]) {
      case W_PDGID: case -W_PDGID:
        if (!isHard(GenParticles_Status[i])) continue;
        Ws.push_back(GenParticles[i]);
        break;
      case E_PDGID: case M_PDGID: case T_PDGID:
      case -E_PDGID: case -M_PDGID: case -T_PDGID:
        Leps.push_back(GenParticles[i]);
        break;
      } // clang-format on
    }
    int num_leps = std::size(Leps);
    auto cut = getHistName(num_leps, NJets, HT_trigger_pt30, JetsAK8 , nsr21, nsr43, nsr42);
#define Fill(table, var)  fillHistos(my_histos, table, var, weight, cut)
#define FillU(table, var) fillHistos(my_histos, table, var, weight, cut)
    Fill("EventCounter", eventCounter);
    Fill("njets", NJets);
    Fill("nCA12Jets", std::size(JetsCA12));
    Fill("nbjets", NGoodBJets_pt30);
    // Fill("nwjets", nwjets);
    Fill("HT_pt30", HT_trigger_pt30);
    Fill("met", MET);
    Fill("MT2", MT2);
    for (const auto &v : Ws) {
      Fill("WEta", v.Eta());
      Fill("WE", v.E());
      Fill("WPhi", v.Phi());
      Fill("WP", v.P());
      Fill("WPt", v.Pt());
    }
    auto computeAngle = [](const auto &pair) {
      return pair[0].Vect().Angle(pair[1].Vect());
    };
    for (std::size_t i = 0; i < std::size(JetsCA12); ++i) {
      Fill("NSubJettiness3", JetsCA12_NsubjettinessTau3[i]);
      Fill("NSubJettiness2", JetsCA12_NsubjettinessTau2[i]);
      Fill("NSubJettiness1", JetsCA12_NsubjettinessTau1[i]);
      Fill("NSubJettiness4", JetsCA12_NsubjettinessTau4[i]);
      Fill("NSubRatio42", nsr42[i]);
      Fill("NSubRatio43", nsr43[i]);
      Fill("NSubRatio21", nsr21[i]);
    }
    if (std::size(Ws) == 2)
      Fill("Gen_W_Angle", computeAngle(Ws));
    if (std::size(Leps) == 2)
      FillU("Gen_Lep_Angle", computeAngle(Leps));
    for (std::size_t i = 0; i < std::size(JetsCA12); ++i) {
        Fill("CA12Pt", JetsCA12[i].Pt());
    }

#if PER_EVENT_LOG
    std::cout << "\n";
#endif
  }
}

void Analyze2W::WriteHistos(TFile *outfile) {
  outfile->cd();
  for (const auto &cut : cuts) {
    std::string name = "WPt_" + cut;
    my_histos[name].mods.push_back([](TH1D *h, TCanvas *c) {
      double frac = (h->Integral(20, 99999)) / (h->Integral(0, 99999));
      addTextToStats(
          h, c, (std::string("Wpt200/W = ") + std::to_string(frac)).c_str());
    });
  }

#if SHOW_PROGRESS
  ProgressBar prog(std::size(my_histos), "Producing Histogram: ");
#endif

  for (const auto &p : my_histos) {
#if SHOW_PROGRESS
    (++prog).display(p.first);
#endif

    p.second.Write();

#if IMAGE_PRODUCTION
    if (p.second.GetEntries())
      p.second.DrawImage("Output");
#endif

  }
  for (const auto &p : my_2d_histos)
    p.second.Write();
//   for (const auto &p : my_efficiencies)
//     p.second.Write();
}
