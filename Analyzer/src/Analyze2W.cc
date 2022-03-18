#define Analyze2W_cxx 
#define SHOW_PROGRESS 0
#define PER_EVENT_LOG 0
#define IMAGE_PRODUCTION 0

#if (SHOW_PROGRESS && PER_EVENT_LOG)
#error "Cannot have both per-event logging and a progress bar"
#endif

#include "Analyzer/Analyzer/include/Analyze2W.h"
#include "2W_Utils.hpp"
#include "Framework/Framework/include/Utility.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "fastjet/ClusterSequence.hh"
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
#include <fstream>

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

struct SliceData {
    double Ht;
    double genHt;
    int NJets;
    int n_gen_leps;
    const std::vector<TLorentzVector> &Jets;
    const std::vector<double> &nsr21;
    const std::vector<double> &nsr42;
    const std::vector<double> &nsr43;
    std::size_t j1_index=0;
    std::size_t j2_index=1;
    bool has_2_jets=false;
};

std::ostream& operator<<(std::ostream &os, const TLorentzVector &v) {
    os << "(" << v.Pt() << ", " << v.Eta() << ", " << v.Phi() << "," << v.E() << ")";
    return os;
}

    template<class ForwardIt, class Compare>
ForwardIt max_element_second(ForwardIt first, ForwardIt last, 
        Compare comp, ForwardIt otherthan)
{
    if (first == last) return last;

    ForwardIt largest = first;
    ++first;
    for (; first != last; ++first) {
        if (comp(*largest, *first) && first != otherthan) {
            largest = first;
        }
    }
    return largest;
}

class LeadingJetGen : public Generator {
    public:
        LeadingJetGen() : Generator("LeadingJetGen"){}
        void calculate(SliceData &data) override {
            auto &Jets = data.Jets;
            /*
               DEBUG("Jets: " << data.Jets);
               auto firstLargestPt =
               std::max_element(Jets.begin(), Jets.end(),
               [](auto &x, auto &y) { return x.Pt() < y.Pt(); });
               auto secondLargestPt = max_element_second(
               Jets.begin(), Jets.end(), [&firstLargestPt](auto &x, auto &y) {
               return  (x.Pt() < y.Pt());
               }, firstLargestPt);
               data.j1_index = std::distance(Jets.begin(), firstLargestPt);
               data.j2_index = std::distance(Jets.begin(), secondLargestPt);
             */
            data.has_2_jets = std::size(Jets) > 1;
        }
};


class HTCut : public Cut {
    public:
        HTCut() : Cut("HTCut", false, {"HT<=900", "HT>900"}) {}
        void calculate(const SliceData &data) override {
            //            DEBUG("Running cut HT");
            if (data.Ht > 900) {
                passed = true;
                value = possible_values[1];
            } else {
                passed = false;
                value = possible_values[0];
            }
        }
};

class GenHTCut : public Cut {
    public:
        GenHTCut() : Cut("GenHTCut", false, {"HT>700", "HT<=700"}) {}
        void calculate(const SliceData &data) override {
            if (data.genHt > 700) {
                passed = true;
                value = possible_values[0];
            } else {
                passed = false;
                value = possible_values[1];
            }
        }
};

class GenLepCut : public Cut {
    public:
        GenLepCut() : Cut("GenLepCut", false, {"0Lep", "1Lep","2Lep"}) {}
        void calculate(const SliceData &data) override {
            switch(data.n_gen_leps){
                case 0:
                    passed = true;
                    value = possible_values[0];
                    break;
                case 1:
                    passed = false;
                    value = possible_values[1];
                    break;
                case 2:
                    passed = false;
                    value = possible_values[2];
                    break;
            }
        }
};

class SelectionCut : public Cut {
    public:
        SelectionCut()
            : Cut("SelectionCut", false, {"Selection", "NotSelection"}) {
            }
        void calculate(const SliceData &data) override {
            if (!data.has_2_jets) {
                value = possible_values[1];
                passed = false;
                return;
            }
            auto& Jets = data.Jets;
            auto& nsr21 = data.nsr21;
            std::size_t j1_index = data.j1_index;
            std::size_t j2_index = data.j2_index;

            auto& firstLargestPt = Jets[j1_index];
            auto& secondLargestPt = Jets[j2_index];
            auto jetpass = [](const TLorentzVector &j, double t21) {
                return j.Pt() > 400 && std::abs(j.Eta()) < 2.0 && t21 < 0.75;
            };
            double largestTau21 = nsr21[j1_index];
            double secondTau21 = nsr21[j2_index];
            if ((jetpass(firstLargestPt, largestTau21) &&
                        jetpass(secondLargestPt, secondTau21))&& data.Ht > 900) {
                value = possible_values[0];
                passed = true;
            } else {
                value = possible_values[1];
                passed = false;
            }
        }
};



class TauCut : public Cut {
    public:
        TauCut() : Cut("TauCut", false, {"TauPass", "TauFail"}) {}
        void calculate(const SliceData &data) override {
            auto &Jets = data.Jets;
            if (!data.has_2_jets) {
                value = possible_values[1];
                passed = false;
                return;
            }
            auto &nsr21 = data.nsr21;
            auto &nsr43 = data.nsr43;
            auto &nsr42 = data.nsr42;
            std::size_t j1_index = data.j1_index;
            std::size_t j2_index = data.j2_index;
            if (nsr43[j1_index] < 0.8 && nsr43[j2_index] < 0.8 &&
                    nsr42[j1_index] < 0.5 && nsr42[j2_index] < 0.5 ){
                passed = true;
                value = possible_values[0];
            } else {
                passed = false;
                value = possible_values[1];
            }
        }
};

class MassRatioCut : public Cut {
    public:
        MassRatioCut() : Cut("MassRatioCut", false, {"APass", "AFail"}) {}
        void calculate(const SliceData &data) override {
            auto &Jets = data.Jets;
            if (!data.has_2_jets) {
                value = possible_values[1];
                return;
            }
            auto &nsr21 = data.nsr21;
            auto &nsr43 = data.nsr43;
            auto &nsr42 = data.nsr42;
            std::size_t j1_index = data.j1_index;
            std::size_t j2_index = data.j2_index;
            auto& firstLargestPt = Jets[j1_index];
            auto& secondLargestPt = Jets[j2_index];

            if ( std::abs(firstLargestPt.M() - secondLargestPt.M()) /
                    (firstLargestPt.M() + secondLargestPt.M()) <
                    0.1 ){
                passed = true;
                value = possible_values[0];
            } else {
                passed = false;
                value = possible_values[1];
            }
        }
};

class EtaCut : public Cut {
    public:
        EtaCut() : Cut("EtaCut", false, {"EtaPass", "NotSig"}) {}
        void calculate(const SliceData &data) override {
            auto &Jets = data.Jets;
            std::size_t j1_index = data.j1_index;
            std::size_t j2_index = data.j2_index;
            auto& firstLargestPt = Jets[j1_index];
            auto& secondLargestPt = Jets[j2_index];
            if ( data.has_2_jets && std::abs(firstLargestPt.Eta() - secondLargestPt.Eta()) < 1.0) {
                passed = true;
                value = possible_values[0];
            } else {
                passed = false;
                value = possible_values[1];
            }
        }
};

const int W_PDGID = 24, E_PDGID = 11, M_PDGID = 13, T_PDGID = 15,
      STOP_PDGID = 1000006;
const int JETCOUNT = 10;
const int LEP_COUNT = 10;
const double W_WP = 0.918;
const double WTAG_PT = 200.0f;

Analyze2W::Analyze2W() {
    my_histos.generators.push_back(std::make_unique<LeadingJetGen>());

    my_histos.cuts.push_back(std::make_unique<GenHTCut>());
    my_histos.cuts.push_back(std::make_unique<GenLepCut>());
    my_histos.cuts.push_back(std::make_unique<SelectionCut>());
    my_histos.cuts.push_back(std::make_unique<MassRatioCut>());
    my_histos.cuts.push_back(std::make_unique<TauCut>());
    my_histos.cuts.push_back(std::make_unique<EtaCut>());
    my_histos.constructChains({
            {"GenLepCut", "GenHTCut"},
            {"GenLepCut", "SelectionCut"},
            {"GenLepCut", "SelectionCut","EtaCut"},
            {"GenLepCut", "SelectionCut","EtaCut", "MassRatioCut"},
            {"GenLepCut", "SelectionCut","EtaCut", "MassRatioCut", "TauCut"},
            });
    InitHistos();
    my_histos.printHistos();
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

    my_histos.createNewHistogram("EventCounter", 2, -1.1, 1.1, "Event Counter");
    my_histos.createNewHistogram("njets", 15, 0, 15, "NJets");
    my_histos.createNewHistogram("nAK8jets", 15, 0, 15, "NAK8Jets");

    my_histos.createNewHistogram("nCA12Jets", 15, 0, 15, "NAK8Jets");

    my_histos.createNewHistogram("nbjets", 5, 0, 5, "NBJets");
    my_histos.createNewHistogram("nwjets", 5, 0, 5, "NWJets_Reco");
    my_histos.createNewHistogram("MT2", 100, 0, 1000, "MT2_All", "Events", true);
    my_histos.createNewHistogram("met", 200, 0, 2000, "MET");

    my_histos.createNewHistogram("NSubJettiness4", 100, 0, 1, "NSJ 4");
    my_histos.createNewHistogram("NSubJettiness3", 100, 0, 1, "NSJ 3");
    my_histos.createNewHistogram("NSubJettiness2", 100, 0, 1, "NSJ 2");
    my_histos.createNewHistogram("NSubJettiness1", 100, 0, 1, "NSJ 1");

    my_histos.createNewHistogram("NSubRatio43", 100, 0, 2, "NSR 43");
    my_histos.createNewHistogram("NSubRatio42", 100, 0, 2, "NSR 42");
    my_histos.createNewHistogram("NSubRatio21", 100, 0, 2, "NSR 21");

    my_histos.createNewHistogram("WE", 200, 0, 2000, "W Energy");
    my_histos.createNewHistogram("WP", 200, 0, 2000, "W P");
    my_histos.createNewHistogram("WPt", 200, 0, 2000, "W PT");
    my_histos.createNewHistogram("WPhi", 100, -3.2, 3.2, "W Phi");
    my_histos.createNewHistogram("WEta", 100, -3.2, 3.2, "W Eta");


    my_histos.createNewHistogram("StopE", 200, 0, 2000, "Stop Energy");
    my_histos.createNewHistogram("StopP", 200, 0, 2000, "Stop P");
    my_histos.createNewHistogram("StopPt", 200, 0, 2000, "Stop PT");
    my_histos.createNewHistogram("StopPhi", 100, -3.2, 3.2, "Stop Phi");
    my_histos.createNewHistogram("StopEta", 100, -3.2, 3.2, "Stop Eta");

    my_histos.createNewHistogram("JetPt", 50, 0, 500, "jet pt (goodjets_pt20)",
            "", true);

    my_histos.createNewHistogram("CA12Pt", 50, 0, 500, "CA12Pt ", "Events", true);

    my_histos.createNewHistogram("HT_pt30", 200, 0, 2000, "HT_pt30", "Events",
            true);

    my_histos.createNewHistogram("Gen_Lep_Angle", 100, 0, 4, "Angle");
    my_histos.createNewHistogram("Gen_W_Angle", 100, 0, 4, "Angle");
}

void Analyze2W::Loop(NTupleReader &tr, double, int maxevents, bool) {
    filetag = tr.getVar<std::string>("filetag");
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

#define MAKE_RATIO(x, y)                                                       \
        auto nsr##x##y =                                                             \
        computeRatio(JetsCA12_NsubjettinessTau##x, JetsCA12_NsubjettinessTau##y)
        MAKE_RATIO(4, 3);
        MAKE_RATIO(4, 2);
        MAKE_RATIO(2, 1);

        makeVar(MET, double);
        makeVar(MT2, double);
        makeVar(HT_trigger_pt30, double);
        makeVar(GenHT, double);

        makeVec(GenParticles, TLorentzVector);
        makeVec(GenParticles_PdgId, int);
        makeVec(GenParticles_ParentId, int);
        makeVec(GenParticles_ParentIdx, int);
        makeVec(GenParticles_Status, int);
        makeVec(JetsAK8_wDiscriminatorDeep, double);

        makeVec(GoodLeptons, (std::pair<std::string, TLorentzVector>));

        std::vector<TLorentzVector> Ws, Leps, Stops, Quarks;
        auto isHard = [](int x) {
            return x == 1 || (std::abs(x) < 60 && std::abs(x) > 20);
        };
        for (std::size_t i = 0; i < GenParticles.size(); ++i) { // clang-format off
#if 0
            std::cout << "      "<< i << ": " << GenParticles_PdgId[i]  << "(" << GenParticles_Status[i] << ") "
                << " -- " << GenParticles_ParentIdx[i] << "(" << GenParticles_ParentId[i] << ")" << std::endl;
#endif
            switch (GenParticles_PdgId[i]) {
                case W_PDGID: case -W_PDGID:
                    if (!isHard(GenParticles_Status[i])) continue;
                    Ws.push_back(GenParticles[i]);
                   // if(GenParticles_ParentId[i] >= 0) {
                   // if (std::abs(GenParticles_PdgId[GenParticles_ParentId[i]]) == STOP_PDGID)
                   //
                   // }
                    break;
                case E_PDGID: case M_PDGID: case T_PDGID:
                case -E_PDGID: case -M_PDGID: case -T_PDGID:
                    Leps.push_back(GenParticles[i]);
                    if( std::abs(GenParticles_PdgId[i]) == 15)
#if PER_EVENT_LOG
                                std::cout << "      "<< i << ": " << GenParticles_PdgId[i]  << "(" << GenParticles_Status[i] << ") " 
                                                    << " -- " << GenParticles_ParentIdx[i] << "(" << GenParticles_ParentId[i] << ")" << std::endl;
#endif
                    break;
                case STOP_PDGID: case -STOP_PDGID:
                    if(std::abs(GenParticles_Status[i])==22)
                        Stops.push_back(GenParticles[i]);
                        break;
                case 1: case -1: case 2: case -2: case 3: case -3:
                case 4: case -4: case 5: case -5: case 6: case -6:
                case 7: case -7: case 8: case -8: 
                        if(std::abs(GenParticles_ParentId[i]) == 24){
                            Quarks.push_back(GenParticles[i]);
                        }
                        break;
            } // clang-format on
        }
        int num_leps = std::size(Leps);
//if(num_leps > 0) continue;


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

        SliceData data{HT_trigger_pt30, GenHT, NJets, num_leps, JetsCA12, nsr21, nsr42, nsr43};
        my_histos.processCuts(data);

#define Fill(table, var) my_histos.fill(table, weight, var)

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
        for (const auto &v : Ws) {
            Fill("WEta", v.Eta());
            Fill("WE", v.E());
            Fill("WPhi", v.Phi());
            Fill("WP", v.P());
            Fill("WPt", v.Pt());
        }
        for (const auto &v : Stops) {
            Fill("StopEta", v.Eta());
            Fill("StopE", v.E());
            Fill("StopPhi", v.Phi());
            Fill("StopP", v.P());
            Fill("StopPt", v.Pt());
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
        for (std::size_t i = 0; i < std::size(JetsCA12); ++i) {
            Fill("CA12Pt", JetsCA12[i].Pt());
        }

#if PER_EVENT_LOG
        std::cout << "\n";
#endif
    }
}

void Analyze2W::WriteHistos(TFile *outfile) {
    std::string outputDir = std::string("Output/")+filetag;
    struct stat info;
    if( stat( outputDir.c_str(), &info ) != 0 ) system((std::string("mkdir -p ")+outputDir).c_str());
    outfile = TFile::Open((outputDir+"/histos.root").c_str(), "RECREATE");
    outfile->cd();
    for (auto &h : my_histos) {
        if (h.first.find("WPt") == std::string::npos)
            continue;
        my_histos[h.first].mods.push_back([](TH1D *h, TCanvas *c) {
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
            p.second.DrawImage("Output"+filetag);
#endif
    }
    for (const auto &p : my_2d_histos)
        p.second.Write();
    //   for (const auto &p : my_efficiencies)
    std::ofstream cutfile;
    cutfile.open("Output/cutflow_out.txt", std::ios_base::app);
    cutfile << filetag << "    ";
    cutfile << my_histos.makeCutflow("EventCounter",true)<< '\n';
    //     p.second.Write();
}
