#define Analyze2W_cxx

#include "Analyzer/Analyzer/include/Analyze2W.h"
#include "Framework/Framework/include/Utility.h"
#include "SusyAnaTools/Tools/NTupleReader.h"

#include <TCanvas.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TPaveStats.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <utility>
#include <vector>

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

 */

using Hist1 = Histogram<TH1D>;
using Hist2 = Histogram<TH2D>;

const int W_PDGID = 24, E_PDGID = 11, M_PDGID = 13, T_PDGID = 6,
          STOP_PDGID = 1000006;
const int JETCOUNT = 10;
const int LEP_COUNT = 10;
const float W_WP = 0.918;

template <typename T> struct argument_type;
template <typename T, typename U> struct argument_type<T(U)> {
  typedef U type;
};

template <typename U, typename T, typename NameExtractor>
void fillHistos(U &histos, const std::string &name, T val, double weight,
                NameExtractor ext, bool uncut = true, bool cut = true) {
  if (cut) {
    const std::string temp = ext(name);
    histos[temp].Fill(val, weight);
  }
  if (uncut)
    histos[name].Fill(val, weight);
}

std::vector<std::string>
getAllNames(std::vector<std::vector<std::string>> vals) {
  std::vector<std::string> ret;
  if (vals.size() == 1) {
    return vals[0];
  } else {
    auto current = vals.back();
    vals.pop_back();
    for (const auto &s : current) {
      std::vector<std::string> next = getAllNames(vals);
      for (std::string n : next) {
        ret.push_back(n + "_" + s);
      }
    }
    return ret;
  }
}

static std::vector<std::string> cuts;
Analyze2W::Analyze2W() {
  cuts = getAllNames({{"l0", "l1", "l2", "lINF"}, {"njets<=4", "njets>4"}});
  InitHistos();
}

template <typename T>
void createNewHistogram(T &myhistos, std::string name, int v1, int v2, int v3,
                        std::string xlabel = "", std::string ylabel = "",
                        bool logaxis = false) {
  std::string basename = /*std::string("h_") + */ name;
  std::string newname =
      basename + ';' + std::move(xlabel) + ';' + std::move(ylabel);
  const char *x = basename.c_str();
  auto ret = Hist1(x, x, v1, v2, v3);
  ret.log_y = logaxis;
  ret.SetTitle(newname.c_str());
  myhistos.emplace(basename, ret);
}

template <typename T>
void createNewHistograms(T &myhistos, const std::string &name, int v1, int v2,
                         int v3, const std::vector<std::string> &cuts,
                         std::string xlabel = "", std::string ylabel = "Events",
                         bool uncut = true, bool logaxis = false) {
  for (const std::string &cut : cuts) {
    std::string newname = name + "_";
    newname += cut;
    createNewHistogram(myhistos, newname, v1, v2, v3, xlabel, ylabel, logaxis);
  }
  if (uncut)
    createNewHistogram(myhistos, name, v1, v2, v3, xlabel, ylabel, logaxis);
}

void Analyze2W::InitHistos() {
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  createNewHistogram(my_histos, "EventCounter", 2, -1.1, 1.1);
  for (int i = 0; i < JETCOUNT; ++i) {
    createNewHistograms(my_histos, std::string("JetPt_") + std::to_string(i),
                        50, 0, 500, cuts, "Jet PT (GoodJets_pt20)", "Events",
                        true, true);
  }
  createNewHistograms(my_histos, "njets", 15, 0, 15, cuts, "NJets");
  createNewHistograms(my_histos, "nbjets", 5, 0, 5, cuts, "NBJets");
  createNewHistograms(my_histos, "MT2", 100, 0, 1000, cuts, "MT2_All", "Events",
                      true, true);
  createNewHistograms(my_histos, "met", 200, 0, 2000, cuts, "MET");
  createNewHistograms(my_histos, "WPhi", 100, -3.2, 3.2, cuts, "W Phi");
  createNewHistograms(my_histos, "WEta", 100, -3.2, 3.2, cuts, "W Eta");
  createNewHistograms(my_histos, "WE", 200, 0, 2000, cuts, "W Energy");
  createNewHistograms(my_histos, "WP", 200, 0, 2000, cuts, "W P");
  createNewHistograms(my_histos, "WPt", 200, 0, 2000, cuts, "W PT");
  createNewHistograms(my_histos, "JetPt", 50, 0, 500, cuts,
                      "Jet PT (GoodJets_pt20)", "", true, true);
  createNewHistograms(my_histos, "LepPt", 100, 0, 1000, cuts,
                      "Lepton PT (GoodLeptons Pt 10)");
  for (int i = 0; i < LEP_COUNT; ++i) {
    createNewHistograms(my_histos, std::string("LepPt_") + std::to_string(i),
                        50, 0, 500, cuts, "Lepton Pt (Goodleptons Pt 10)",
                        "Events", true, true);
  }

  createNewHistogram(my_histos, "Gen_Lep_Angle", 100, 0, 4, "Angle");
  createNewHistograms(my_histos, "Gen_W_Angle", 100, 0, 4, cuts, "Angle");
}

void Analyze2W::Loop(NTupleReader &tr, double, int maxevents, bool) {
  float progress = 0.0;
  while (tr.getNextEvent()) {

#if 1
    float barWidth = 40;

    progress = float(tr.getEvtNum()) / maxevents;
    std::cout << "[";

    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
      if (i < pos)
        std::cout << "=";
      else if (i == pos)
        std::cout << ">";
      else
        std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << "% Processing Event "
              << tr.getEvtNum() << "\r";
    std::cout.flush();
    if (maxevents != -1 && tr.getEvtNum() >= maxevents) {
      std::cout << std::endl;
      break;
    }
#endif
    const auto &eventCounter = tr.getVar<int>("eventCounter");
    my_histos["EventCounter"].Fill(eventCounter);

    // if (tr.getEvtNum() & (10000 == 0))
    //   printf(" Event %i\n", tr.getEvtNum());

#define makeVar(name, type) const type &name = tr.getVar<type>(#name)
#define makeVec(name, t)                                                       \
  const std::vector<argument_type<void(t)>::type> &name =                      \
      tr.getVec<argument_type<void(t)>::type>(#name)

    makeVar(runtype, std::string);
    makeVar(NJets, int);
    makeVar(JetID, bool);
    makeVar(NGoodLeptons, int);
    makeVec(Jets, TLorentzVector);
    makeVec(GoodJets_pt20, bool);
    makeVar(NGoodBJets_pt30, int);

    makeVar(MET, double);
    makeVar(METPhi, double);
    makeVar(MT2, double);

    makeVec(Muons, TLorentzVector);
    makeVec(Electrons, TLorentzVector);
    makeVec(GenElectrons, TLorentzVector);
    makeVec(GenMuons, TLorentzVector);
    makeVec(GenParticles, TLorentzVector);
    makeVec(GenParticles_PdgId, int);
    makeVec(GenParticles_Status, int);
    // makeVec(GenParticles_ParentId, int);
    // makeVec(GenParticles_ParentIdx, int);
    makeVec(GoodLeptons, (std::pair<std::string, TLorentzVector>));

    // std::sort(Jets.begin(), Jets.end(), [](auto x, auto y){return x.Pt() <
    // y.Pt();}); std::sort(GoodLeptons.begin(), GoodLeptons.end(), [](auto x,
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
    int found = 0;
    int w1s = -9999, w2s = -9999;
    auto isHard = [](int x) {
      return x == 1 || (std::abs(x) < 60 && std::abs(x) > 20);
    };
    for (int i = 0; i < GenParticles.size(); ++i) {
      switch (GenParticles_PdgId[i]) {
      case W_PDGID:
      case -W_PDGID:
        if (w1s > -100)
          w2s = GenParticles_Status[i];
        else
          w1s = GenParticles_Status[i];
        if (!isHard(GenParticles_Status[i]))
          continue;
        Ws.push_back(GenParticles[i]);
        break;
      case E_PDGID: case M_PDGID: case T_PDGID:
      case -E_PDGID: case -M_PDGID: case -T_PDGID:
        Leps.push_back(GenParticles[i]);
        break;
      }
    }

    int num_leps = std::size(Leps);
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

    auto getHistName = [&](std::string base) {
      // base = std::string("h_") + std::move(base);
      base += lepsel();
      base += njetsel();
      return base;
    };

#define Fill(table, var)                                                       \
  fillHistos(my_histos, table, var, weight, getHistName, true, true)
#define FillU(table, var)                                                      \
  fillHistos(my_histos, table, var, weight, getHistName, true, false)
    Fill("njets", NJets);
    Fill("nbjets", NGoodBJets_pt30);
    Fill("met", MET);
    Fill("MT2", MT2);
    for (const auto &v : Ws) {
      Fill("WEta", v.Eta());
      Fill("WE", v.E());
      Fill("WPhi", v.Phi());
      Fill("WP", v.P());
      Fill("WPt", v.Pt());
    }

    if (std::size(Ws) == 2) {
      Fill("Gen_W_Angle", Ws[0].Vect().Angle(Ws[1].Vect()));
    } else {
    }
    if (std::size(Leps) == 2)
      FillU("Gen_Lep_Angle", Leps[0].Vect().Angle(Leps[1].Vect()));

    for (int i = 0; i < std::size(Jets); ++i) {
      if (GoodJets_pt20[i])
        Fill("JetPt", Jets[i].Pt());
      if (i < JETCOUNT)
        Fill(std::string("JetPt_") + std::to_string(i), Jets[i].Pt());
    }

    for (int i = 0; i < std::size(GoodLeptons); ++i) {
      const auto tlv = GoodLeptons[i].second;
      const double pt = tlv.Pt();
      if (pt >= 10)
        Fill("LepPt", pt);
      if (i < LEP_COUNT)
        Fill(std::string("LepPt_") + std::to_string(i), pt);
    }
  }
}

void Analyze2W::WriteHistos(TFile *outfile) {
  outfile->cd();
  for (const auto &cut : cuts) {
    std::string name = "WPt_" + cut;
    my_histos[name].mods.push_back([](TH1D *h, TCanvas *c) {
      double frac = (h->Integral(20, 99999)) / (h->Integral(0, 99999));
      TPaveStats *st = static_cast<TPaveStats *>(c->GetPrimitive("stats"));
      st->SetName("mystats");
      TList *listOfLines = st->GetListOfLines();
      TLatex *t = new TLatex(
          0, 0, (std::string("Wpt200/W = ") + std::to_string(frac)).c_str());
      t->SetTextFont(42);
      t->SetTextColor(kBlack);
      t->SetTextSize(0.03);
      listOfLines->Add(t);
      h->SetStats(0);
    });
  }
  for (const auto &p : my_histos) {
    p.second.Write();
    p.second.DrawImage("IMAGES");
  }
  for (const auto &p : my_2d_histos)
    p.second.Write();
  for (const auto &p : my_efficiencies)
    p.second.Write();
}
