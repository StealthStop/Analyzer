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
#include <map>
#include <numeric>
#include <utility>
#include <vector>

const int W_PDGID = 24;
const int T_PDGID = 6;
const int STOP_PDGID = 1000006;

template <typename T> struct argument_type;
template <typename T, typename U> struct argument_type<T(U)> {
  typedef U type;
};

template <typename T, typename NameExtractor>
void fillHistos(auto &histos, const std::string &name, T val, double weight,
                NameExtractor ext, bool uncut) {
  const std::string temp = ext(name);
  histos[temp]->Fill(val, weight);
  if (uncut)
    histos[std::string("h_") + name]->Fill(val, weight);
}

Analyze2W::Analyze2W() { InitHistos(); }

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

void createNewHistogram(auto &myhistos, std::string name, int v1, int v2,
                        int v3) {
  const std::string newname = std::string("h_") + name;
  const char *x = newname.c_str();
  auto ret = std::make_shared<TH1D>(x, x, v1, v2, v3);
  myhistos.emplace(newname, std::move(ret));
}

void createNewHistograms(auto &myhistos, const std::string &name, int v1,
                         int v2, int v3, const std::vector<std::string> &cuts,
                         bool uncut = true) {

  for (const std::string &cut : cuts) {
    std::string newname = name + "_";
    newname += cut;
    createNewHistogram(myhistos, newname, v1, v2, v3);
  }
  if (uncut)
    createNewHistogram(myhistos, name, v1, v2, v3);
}

void Analyze2W::InitHistos() {
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  std::vector<std::string> cuts = getAllNames({{"l0", "l1", "l2", "lINF"}});

  createNewHistogram(my_histos, "EventCounter", 2, -1.1, 1.1);

  createNewHistograms(my_histos, "njets", 15, 0, 15, cuts);
  createNewHistograms(my_histos, "MT2", 120, 0, 600, cuts);
  createNewHistograms(my_histos, "met", 200, 0, 2000, cuts);
  createNewHistograms(my_histos, "WPhi", 100, -3.2, 3.2, cuts);
  createNewHistograms(my_histos, "WEta", 100, -3.2, 3.2, cuts);
  createNewHistograms(my_histos, "WE", 200, 0, 2000, cuts);
  createNewHistograms(my_histos, "WP", 200, 0, 2000, cuts);
  createNewHistograms(my_histos, "WPt", 200, 0, 2000, cuts);
}

void Analyze2W::Loop(NTupleReader &tr, double, int maxevents, bool) {
  float progress = 0.0;
  while (tr.getNextEvent()) {
    float barWidth = 80;
    progress = float(tr.getEvtNum()) / maxevents;
    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
      if (i < pos) std::cout << "=";
      else if (i == pos) std::cout << ">";
      else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << "% Processing Event "
              << tr.getEvtNum() << "\r";
    std::cout.flush();
    if (maxevents != -1 && tr.getEvtNum() >= maxevents) {
      std::cout << std::endl;
      break;
    }
    const auto &eventCounter = tr.getVar<int>("eventCounter");
    my_histos["h_EventCounter"]->Fill(eventCounter);

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
    makeVar(MET, double);
    makeVar(METPhi, double);
    makeVar(MT2, double);
    makeVec(GenElectrons, TLorentzVector);
    makeVec(GenMuons, TLorentzVector);
    makeVec(GenParticles, TLorentzVector);
    makeVec(GenParticles_ParentId, int);
    makeVec(GenParticles_ParentIdx, int);
    makeVec(GenParticles_PdgId, int);
    makeVec(GenParticles_Status, int);
    makeVec(GoodLeptons, (std::pair<std::string, TLorentzVector>));

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

    double W1Phi, W2Phi, W1Eta, W2Eta, W1E, W2E, W1P, W2P, W2Pt, W1Pt;
    int found = 0;
    for (int i = 0; i < GenParticles.size(); ++i) {
      if (GenParticles_PdgId[i] == W_PDGID) {
        const TLorentzVector &v = GenParticles[i];
        W1P = v.P();
        W1Pt = v.Pt();
        W1Phi = v.Phi();
        W1Eta = v.Eta();
        W1E = v.Energy();
        ++found;
      } else if (GenParticles_PdgId[i] == -W_PDGID) {
        const TLorentzVector &v = GenParticles[i];
        W2Phi = v.Phi();
        W2P = v.P();
        W2Pt = v.Pt();
        W2Eta = v.Eta();
        W2E = v.Energy();
        ++found;
      }
    }

    int num_leps = GenElectrons.size() + GenMuons.size();

    auto lepsel = [num_leps]() {
      const std::map<int, std::string> vals = {
          {0, "l0"}, {1, "l1"}, {2, "l2"}, {99999, "linf"}};
      const auto x =
          std::lower_bound(vals.begin(), vals.end(), num_leps,
                           [](auto &&x, int y) { return x.first < y; });
      return x->second;
    };

    auto getHistName = [&](std::string base) {
      base = std::string("h_") + std::move(base);
      base += std::string("_") + lepsel();
      return base;
    };

#define Fill(table, var)                                                       \
  fillHistos(my_histos, table, var, weight, getHistName, true)
    Fill("njets", NJets);
    Fill("met", MET);
    Fill("MT2", MT2);
    Fill("WEta", W1Eta);
    Fill("WE", W1E);
    Fill("WPhi", W1Phi);
    Fill("WP", W1P);
    Fill("WPt", W1Pt);
    Fill("WEta", W2Eta);
    Fill("WE", W2E);
    Fill("WPhi", W2Phi);
    Fill("WP", W2P);
    Fill("WPt", W2Pt);
  }
}

void Analyze2W::WriteHistos(TFile *outfile) {
  outfile->cd();
  for (const auto &p : my_histos)
    p.second->Write();
  for (const auto &p : my_2d_histos)
    p.second->Write();
  for (const auto &p : my_efficiencies)
    p.second->Write();
}
