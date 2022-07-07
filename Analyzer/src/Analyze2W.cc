#define Analyze2W_cxx #pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#include <TPaveStats.h>
#pragma GCC diagnostic pop

#if (SHOW_PROGRESS && PER_EVENT_LOG)
#error "Cannot have both per-event logging and a progress bar"
#endif

#include "NTupleReader/include/NTupleReader.h"
#include "Analyzer/Analyzer/include/Analyze2W.h"
#include "Framework/Framework/include/Utility.h"

#include "fastjet/ClusterSequence.hh"
#include "Framework/Framework/include/Jet.h"

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


namespace detail {
    template <class BidirIter>
        void rotate_discontinuous(
                BidirIter first1, BidirIter last1,
                typename std::iterator_traits<BidirIter>::difference_type d1,
                BidirIter first2, BidirIter last2,
                typename std::iterator_traits<BidirIter>::difference_type d2) {
            using std::swap;
            if (d1 <= d2)
                std::rotate(first2, std::swap_ranges(first1, last1, first2), last2);
            else {
                BidirIter i1 = last1;
                while (first2 != last2) swap(*--i1, *--last2);
                std::rotate(first1, i1, last1);
            }
        }
    template <class BidirIter, class Function>
        bool combine_discontinuous(
                BidirIter first1, BidirIter last1,
                typename std::iterator_traits<BidirIter>::difference_type d1,
                BidirIter first2, BidirIter last2,
                typename std::iterator_traits<BidirIter>::difference_type d2, Function &f,
                typename std::iterator_traits<BidirIter>::difference_type d = 0) {
            typedef typename std::iterator_traits<BidirIter>::difference_type D;
            using std::swap;
            if (d1 == 0 || d2 == 0) return f();
            if (d1 == 1) {
                for (BidirIter i2 = first2; i2 != last2; ++i2) {
                    if (f()) return true;
                    swap(*first1, *i2);
                }
            } else {
                BidirIter f1p = std::next(first1);
                BidirIter i2 = first2;
                for (D d22 = d2; i2 != last2; ++i2, --d22) {
                    if (combine_discontinuous(f1p, last1, d1 - 1, i2, last2, d22, f,
                                d + 1))
                        return true;
                    swap(*first1, *i2);
                }
            }
            if (f()) return true;
            if (d != 0)
                rotate_discontinuous(first1, last1, d1, std::next(first2), last2,
                        d2 - 1);
            else
                rotate_discontinuous(first1, last1, d1, first2, last2, d2);
            return false;
        }

    // Creates a functor with no arguments which calls f_(first_, last_).
    //   Also has a variant that takes two It and ignores them.
    template <class Function, class It>
        class bound_range {
            Function f_;
            It first_;
            It last_;

            public:
            bound_range(Function f, It first, It last)
                : f_(f), first_(first), last_(last) {}

            bool operator()() { return f_(first_, last_); }

            bool operator()(It, It) { return f_(first_, last_); }
        };
};  // namespace detail
template <class BidirIter, class Function>
Function for_each_combination(BidirIter first, BidirIter mid, BidirIter last,
        Function f) {
    detail::bound_range<Function &, BidirIter> wfunc(f, first, mid);
    detail::combine_discontinuous(first, mid, std::distance(first, mid), mid,
            last, std::distance(mid, last), wfunc);
    return std::move(f);
}
template <class It>
unsigned display(It begin, It end) {
    unsigned r = 0;
    if (begin != end) {
        std::cout << *begin;
        ++r;
        for (++begin; begin != end; ++begin) {
            std::cout << ", " << *begin;
            ++r;
        }
    }
    return r;
}




const double WTAG_MEDIUM = 0.960;
const double W_JET_MATCH_DR = 0.2;


template< class... >
using void_t = void;

struct nonesuch {
    nonesuch() = delete;
    ~nonesuch() = delete;
    nonesuch(nonesuch const&) = delete;
    void operator=(nonesuch const&) = delete;
};


template <typename T> struct argument_type;
template <typename T, typename U> struct argument_type<T(U)> {
    typedef U type;
};

template <typename U, typename T>
void fillHistos(U &histos, const std::string &name, T val, double weight,
        const std::string& append, bool uncut = true, bool cut = true) {
    if (cut) {
        histos[name+append].Fill(val, weight);
    }
    if (uncut)
        histos[name].Fill(val, weight);
}

template <typename U, typename T, std::size_t N>
void fillHistos(U &histos, const std::string &name, T val, double weight,
        const std::array<std::string,N>& append){
    for(std::size_t i = 0 ; i < N; ++i){ 
        histos[name+append[i]].Fill(val, weight);
    }
}


/*
   template <typename T> typename T::value_type getAllNames(T vals) {
   typename T::value_type ret;
   ret.reserve(std::accumulate(vals.begin(), vals.end(), 1, [](auto x, const auto& y){return x*(std::size(y)+1);}));
   for(int i = 0 ; i < std::size(vals); ++i){
   int s = ret.size();
   for (int j = 0 ; j < s; ++j){
   for(const auto& y: vals[i]){
   ret.push_back(ret[j] + "_" +  y);
   }
   }
   for(const auto& y: vals[i]){
   ret.push_back( y);
   }
   }
   return ret;
   }
 */


template<typename T, typename U>
T constexpr power(T base, U exponent) {
    return exponent == 0 ? 1 : base * pow(base, exponent - 1);
}

template <typename T, std::size_t N>
auto combine(const std::array<T,N>& vals) {
    std::array<T, power(static_cast<std::size_t>(2),N)> ret;
    ret[0] = "";
    for(int i = 0; i < N; ++i){
        int cur = power(2,i);
        for(int j = 0 ; j < cur; ++j){
            ret[cur + j]=ret[j] + vals[i];
        }
    }
    return ret;
}
void addTextToStats(TH1D *h, TCanvas *c, const std::string& s){
    TPaveStats *st = static_cast<TPaveStats *>(c->GetPrimitive("stats"));
    st->SetName("mystats");
    TList *listOfLines = st->GetListOfLines();
    TLatex *t = new TLatex(0, 0, s.c_str());
    t->SetTextFont(42);
    t->SetTextColor(kBlack);
    t->SetTextSize(0.03);
    listOfLines->Add(t);
    h->SetStats(0);
}





template <typename T>
void createNewHistogram(T &myhistos, std::string name, int v1, int v2, int v3,
        std::string xlabel = "", std::string ylabel = "",
        bool logaxis = false) {
    std::string basename = /*std::string("h_") + */ name;
    std::string newname =
        basename + ';' + std::move(xlabel) + ';' + std::move(ylabel);
    const char *x = basename.c_str();
    auto ret = Histogram<TH1D>(x, x, v1, v2, v3);
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


struct SliceData {
    double Ht;
    double genHt;
    int NJets;
    int n_gen_leps;
    const std::vector<utility::LorentzVector> &Jets;
    const std::vector<float> &nsr21;
    const std::vector<float> &nsr42;
    const std::vector<float> &nsr43;
    std::size_t j1_index=0;
    std::size_t j2_index=1;
    bool has_2_jets=false;
    double mass_ratio = 0;
};

std::ostream& operator<<(std::ostream &os, const utility::LorentzVector &v) {
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
            data.has_2_jets = std::size(Jets) > 1;
        }
};


class HTCut : public Cut {
    public:
        HTCut() : Cut("HTCut", false, {"HT<=900", "HT>900"}) {}
        void calculate(SliceData &data) override {
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
        void calculate(SliceData &data) override {
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
        void calculate(SliceData &data) override {
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
        void calculate(SliceData &data) override {
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
            auto jetpass = [](const utility::LorentzVector &j, double t21) {
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
        void calculate(SliceData &data) override {
            //auto &Jets = data.Jets;
            if (!data.has_2_jets) {
                value = possible_values[1];
                passed = false;
                return;
            }
            //auto &nsr21 = data.nsr21;
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
        void calculate(SliceData &data) override {
            auto &Jets = data.Jets;
            if (!data.has_2_jets) {
                value = possible_values[1];
                return;
            }
            //auto &nsr21 = data.nsr21;
            //auto &nsr43 = data.nsr43;
            //auto &nsr42 = data.nsr42;
            std::size_t j1_index = data.j1_index;
            std::size_t j2_index = data.j2_index;
            auto& firstLargestPt = Jets[j1_index];
            auto& secondLargestPt = Jets[j2_index];
            auto mass_ratio = std::abs(firstLargestPt.M() - secondLargestPt.M()) /
                (firstLargestPt.M() + secondLargestPt.M());
            data.mass_ratio=mass_ratio;

            if (mass_ratio< 0.1){
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
        void calculate(SliceData &data) override {
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
const std::size_t MAX_JETS=6;

Analyze2W::Analyze2W() {
    my_histos.generators.push_back(std::make_unique<LeadingJetGen>());

    my_histos.cuts.push_back(std::make_unique<GenHTCut>());
    my_histos.cuts.push_back(std::make_unique<GenLepCut>());
    my_histos.cuts.push_back(std::make_unique<SelectionCut>());
    my_histos.cuts.push_back(std::make_unique<MassRatioCut>());
    my_histos.cuts.push_back(std::make_unique<TauCut>());
    my_histos.cuts.push_back(std::make_unique<EtaCut>());
    my_histos.constructChains({
            {"GenLepCut"}
            //{"GenLepCut", "GenHTCut"},
            //{"GenLepCut", "SelectionCut"},
            //{"GenLepCut", "SelectionCut","EtaCut"},
            //{"GenLepCut", "SelectionCut","EtaCut", "MassRatioCut"},
            //{"GenLepCut", "SelectionCut","EtaCut", "MassRatioCut", "TauCut"},
            });
    InitHistos();
    //my_histos.printHistos();
}

std::vector<float> computeRatio(std::vector<float> v1,
        const std::vector<float> &v2) {
    for (std::size_t i = 0; i < std::size(v1); ++i) {
        v1[i] /= v2[i];
    }
    return v1;
}

void Analyze2W::InitHistos() {
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    my_histos.createNewHistogram("EventCounter", 2, -1.1, 1.1, "Event Counter");

    my_histos.createNewHistogram("NJets_pt20", 15, 0, 15, "NJets PT20");
    my_histos.createNewHistogram("NJets_pt30", 15, 0, 15, "NJets PT30");

    my_histos.createNewHistogram("HT_pt30", 200, 0, 2000, "HT_pt30", "Events",
            true);

    // Number of reco b jets
    my_histos.createNewHistogram("nbjets_loose", 5, 0, 5, "NBjets_loose");
    my_histos.createNewHistogram("DeepAK8TagW_medium_wp", 5, 0, 5, "DeepAK8TagW_medium_wp");
    my_histos.createNewHistogram("DeepAK8TagW_HasNearGenW", 5, 0, 5, "DeepAK8TagW_HasNearGenW");
    my_histos.createNewHistogram("DeepAK8TagW_MatchedEnergyRatio", 100,0,3, "DeepAK8TagW_MatchedEnergyRatio");


    my_histos.createNewHistogram("RecoStopMass", 150,200,1400, "RecoStopMass");
    my_histos.createNewHistogram("RecoStopMassImbalance", 100,0, 500 , "RecoStopMassImbalance");
    my_histos.createNewHistogram("RecoSBottomMass", 120,50,650, "RecoSBottomMass");
    my_histos.createNewHistogram("RecoSBottomMassImbalance", 100,0, 500 , "RecoSBottomMassImbalance");




        // CA12 jet information
        my_histos.createNewHistogram("nCA12Jets", 15, 0, 15, "nCA12Jets");
    my_histos.createNewHistogram("CA12Pt", 50, 0, 500, "CA12Pt ", "Events", true);

    // Number of reco w jets
    //my_histos.createNewHistogram("nwjets", 5, 0, 5, "NWJets_Reco");


    // MET
    my_histos.createNewHistogram("met", 200, 0, 2000, "MET");

    // Nsubjetiness of CA12 jets
    my_histos.createNewHistogram("CA12_NSubJettiness4", 100, 0, 1, "CA12_NSJ 4");
    my_histos.createNewHistogram("CA12_NSubJettiness3", 100, 0, 1, "CA12_NSJ 3");
    my_histos.createNewHistogram("CA12_NSubJettiness2", 100, 0, 1, "CA12_NSJ 2");
    my_histos.createNewHistogram("CA12_NSubJettiness1", 100, 0, 1, "CA12_NSJ 1");
    my_histos.createNewHistogram("CA12_NSubRatio43", 100, 0, 2, "CA12_NSR 43");
    my_histos.createNewHistogram("CA12_NSubRatio42", 100, 0, 2, "CA12_NSR 42");
    my_histos.createNewHistogram("CA12_NSubRatio21", 100, 0, 2, "CA12_NSR 21");

    my_histos.createNewHistogram("AK8_NSubJettiness3", 100, 0, 1, "AK8_NSJ 3");
    my_histos.createNewHistogram("AK8_NSubJettiness2", 100, 0, 1, "AK8_NSJ 2");
    my_histos.createNewHistogram("AK8_NSubJettiness1", 100, 0, 1, "AK8_NSJ 1");
    my_histos.createNewHistogram("AK8_NSubRatio21", 100, 0, 2, "AK8_NSR 21");
    my_histos.createNewHistogram("AK8_NSubRatio31", 100, 0, 2, "AK8_NSR 31");
    my_histos.createNewHistogram("AK8_NSubRatio32", 100, 0, 2, "AK8_NSR 32");

    // Gen W information
    my_histos.createNewHistogram("WE", 200, 0, 2000, "W Energy");
    my_histos.createNewHistogram("WP", 200, 0, 2000, "W P");
    my_histos.createNewHistogram("WPt", 200, 0, 2000, "W PT");
    my_histos.createNewHistogram("WPhi", 100, -3.2, 3.2, "W Phi");
    my_histos.createNewHistogram("WEta", 100, -3.2, 3.2, "W Eta");


    // Gen Stop information
    my_histos.createNewHistogram("StopE", 200, 0, 2000, "Stop Energy");
    my_histos.createNewHistogram("StopP", 200, 0, 2000, "Stop P");
    my_histos.createNewHistogram("StopPt", 200, 0, 2000, "Stop PT");
    my_histos.createNewHistogram("StopPhi", 100, -3.2, 3.2, "Stop Phi");
    my_histos.createNewHistogram("StopEta", 100, -3.2, 3.2, "Stop Eta");

    // Leptons
    my_histos.createNewHistogram("NGoodLeptons", 100, 0, 10, "N Good Leptons");
    for(const std::string& s: {"l", "e", "m"}){
        my_histos.createNewHistogram(s + "_Pt" ,   360,    0, 1500 , s + "_Pt" );
        my_histos.createNewHistogram(s + "_Phi" ,   200,    -4, 4, s + "_Phi" );
        my_histos.createNewHistogram(s + "_Eta" ,   200,    -6, 6, s + "_Eta" );
        my_histos.createNewHistogram(s + "_Charge" , 2, -1, 1, s + "_Charge" );
        my_histos.createNewHistogram(s + "_MiniIso" , 100, 0, 1, s + "_MiniIso" );
    }

    // Angles between Ws and Leptons
    //my_histos.createNewHistogram("Gen_Lep_Angle", 100, 0, 4, "Angle");
    my_histos.createNewHistogram("Gen_W_Angle", 100, 0, 4, "Angle");

    // Mass difference of two leading CA12 jets
    my_histos.createNewHistogram("mass_ratio", 100, 0, 1, "Mass Difference Ratio");

    // Jet kinematic information
    for(std::size_t i = 0; i < 6; ++i){
        my_histos.createNewHistogram("Jet_" + std::to_string(i) + "_E", 200, 0, 2000,"Jet_" + std::to_string(i) + "_E" );
        my_histos.createNewHistogram("Jet_" + std::to_string(i) + "_P", 200, 0, 2000, "Jet_" + std::to_string(i) + "_P");
        my_histos.createNewHistogram("Jet_" + std::to_string(i) + "_Pt", 200, 0, 2000, "Jet_" + std::to_string(i) + "_Pt");
        my_histos.createNewHistogram("Jet_" + std::to_string(i) + "_Phi", 100, -3.2, 3.2, "Jet_" + std::to_string(i) + "_Phi");
        my_histos.createNewHistogram("Jet_" + std::to_string(i) + "_Eta", 100, -3.2, 3.2, "Jet_" + std::to_string(i) + "_Eta");
    }

}

void Analyze2W::Loop(NTupleReader &tr, double, int maxevents, bool) {
    filetag = tr.getVar<std::string>("filetag");
    while (tr.getNextEvent()) {
        if (maxevents != -1 && tr.getEvtNum() >= maxevents) {
            std::cout << std::endl;
            break;
        }
        const auto &eventCounter = tr.getVar<int>("eventCounter");

        if (tr.getEvtNum() & (1000 == 0))
            printf(" Event %i\n", tr.getEvtNum());

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
        GET_MACRO(__VA_ARGS__, makeVarCustom, makeVarIdentical, void)(__VA_ARGS__)

        makeVar(runtype, std::string);

        makeVar(NJets_pt20, int);
        makeVar(NJets_pt30, int);
        //makeVec(Jets, utility::LorentzVector);

        makeVec(JetsAK8, utility::LorentzVector);
        makeVec(Jets, utility::LorentzVector); // AK4

        makeVec(JetsAK8_DeepTagWvsQCD, float);

        makeVar(NBJets_loose, int);

        makeVec(JetsCA12, utility::LorentzVector);
        makeVec(JetsCA12_NsubjettinessTau4, float);
        makeVec(JetsCA12_NsubjettinessTau3, float);
        makeVec(JetsCA12_NsubjettinessTau2, float);
        makeVec(JetsCA12_NsubjettinessTau1, float);
        makeVec(JetsAK8_NsubjettinessTau3, float);
        makeVec(JetsAK8_NsubjettinessTau2, float);
        makeVec(JetsAK8_NsubjettinessTau1, float);


#define MAKE_RATIO(name, x, y)                                                       \
        auto nsr##name##_##x##y =                                                             \
        computeRatio(Jets##name##_NsubjettinessTau##x, Jets##name##_NsubjettinessTau##y)
        MAKE_RATIO(CA12,4, 3);
        MAKE_RATIO(CA12,4, 2);
        MAKE_RATIO(CA12,2, 1);

        MAKE_RATIO(AK8,3, 2);
        MAKE_RATIO(AK8,3, 1);
        MAKE_RATIO(AK8,2, 1);

        makeVar(HT_trigger_pt30, double);
        makeVar(MET, float);
        makeVar(GenHT, float);

        makeVec(GenParticles, utility::LorentzVector);
        makeVec(GenParticles_PdgId, int);
        makeVec(GenParticles_ParentId, int);
        makeVec(GenParticles_ParentIdx, int);
        makeVec(GenParticles_Status, int);


        makeVec(GoodLeptons, (std::pair<std::string, utility::LorentzVector>));
        makeVec(GoodLeptonsCharge, int);
        makeVec(GoodLeptonsMiniIso, double);


        std::vector<utility::LorentzVector> Ws, Leps, Stops, Quarks;
        auto isHard = [](int x) {
            return x == 1 || (std::abs(x) < 60 && std::abs(x) > 20);
        };

        for (std::size_t i = 0; i < GenParticles.size(); ++i) { // clang-format off
            switch (GenParticles_PdgId[i]) {
                case W_PDGID: case -W_PDGID:
                    if (!isHard(GenParticles_Status[i])) continue;
                    if(GenParticles_ParentIdx[i] >= 0) {
                        if (GenParticles_ParentId[i] == STOP_PDGID || 
                                GenParticles_ParentId[i] == -STOP_PDGID)
                            Ws.push_back(GenParticles[i]);
                    }
                    break;
                case E_PDGID: case M_PDGID: case T_PDGID:
                case -E_PDGID: case -M_PDGID: case -T_PDGID:
                    if(std::abs(GenParticles_ParentId[i]) == W_PDGID) Leps.push_back(GenParticles[i]);
                    if( std::abs(GenParticles_PdgId[i]) == 15)
                        break;
                case STOP_PDGID: case -STOP_PDGID:
                    if(std::abs(GenParticles_Status[i])==62)
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

        // W tagging
        int nw_deep_tag = 0; //std::count_if(JetsAK8_DeepTagWvsQCD.begin(), JetsAK8_DeepTagWvsQCD.end(),[&](auto&& x){return x>WTAG_MEDIUM});
        int nw_deep_tagjet_has_gen_w = 0; 
        std::vector<int> deep_w_jet_tags;
        std::vector<utility::LorentzVector> deep_w_jets;
        std::vector<std::pair<std::size_t, std::size_t>> matched_w_jets;
        for(std::size_t i = 0 ; i < std::size(JetsAK8); ++i){
            if(JetsAK8_DeepTagWvsQCD[i] > WTAG_MEDIUM) {
                ++nw_deep_tag;
                deep_w_jets.push_back(JetsAK8[i]);
                for(std::size_t j = 0 ; j < std::size(Ws); ++j){
                    if(ROOT::Math::VectorUtil::DeltaR(Ws[j],JetsAK8[i])<W_JET_MATCH_DR){
                        matched_w_jets.push_back({i,j});
                    }
                }
            }
        }

        std::vector<utility::LorentzVector> ak4_exclusive;

        if(nw_deep_tag==2){
            for(const auto& jet : Jets ){
                bool is_close = false;
                for(const auto& wjet : deep_w_jets){
                    if( ROOT::Math::VectorUtil::DeltaR(jet,wjet)<0.8){
                        is_close = true;
                    }
                }
                if(!is_close){
                    ak4_exclusive.push_back(jet);
                }
                if(std::size(ak4_exclusive) > 4) break;
            }
        }
        std::array<int, 4> jets_idx_combos = {1,2,3,4};


        float imbalance=1000000.0f, reco_stop_mass = 0, reco_sbottom_mass=0, reco_sbottom_imbalance;
        auto callinner = [&](auto&& start, auto&& end){
                    float reco_stop_mass_test_1  = (deep_w_jets[0] + ak4_exclusive[jets_idx_combos[0]] + ak4_exclusive[jets_idx_combos[1]]).M();
                    float reco_stop_mass_test_2  = (deep_w_jets[1] + ak4_exclusive[jets_idx_combos[2]] + ak4_exclusive[jets_idx_combos[3]]).M();
                    float test_imbalance=std::abs(reco_stop_mass_test_1 - reco_stop_mass_test_2)/2 ;
                    if(test_imbalance < imbalance){
                        imbalance = test_imbalance;
                        reco_stop_mass = (reco_stop_mass_test_1 + reco_stop_mass_test_2)/2;
                    float reco_sbottom_mass_1  = (ak4_exclusive[jets_idx_combos[0]] + ak4_exclusive[jets_idx_combos[1]]).M();
                    float reco_sbottom_mass_2  = (ak4_exclusive[jets_idx_combos[2]] + ak4_exclusive[jets_idx_combos[3]]).M();
                    reco_sbottom_imbalance = std::abs(reco_sbottom_mass_1 - reco_sbottom_mass_2) / 2.0f;
                    reco_sbottom_mass = (reco_sbottom_mass_1 + reco_sbottom_mass_2) / 2.0f;
                    }
                    return false;
        };

        auto callouter =// [&deep_w_jets, &ak4_exclusive, &imbalance, &reco_stop_mass]
            [&](auto&& start, auto&& end){ for_each_combination(jets_idx_combos.begin(), jets_idx_combos.begin()+2, jets_idx_combos.end(), callinner); return false;};
        if(std::size(ak4_exclusive) >= 4 && nw_deep_tag >= 2){
            for_each_combination(deep_w_jets.begin(), deep_w_jets.begin()+2, deep_w_jets.end(), callouter);
        }


        double weight = 1.0, eventweight = 1.0, leptonScaleFactor = 1.0,
               bTagScaleFactor = 1.0, htDerivedScaleFactor = 1.0,
               prefiringScaleFactor = 1.0, puScaleFactor = 1.0;

        if (runtype == "MC") {
            const auto &Weight = tr.getVar<float>("Weight");
            const auto &lumi = tr.getVar<double>("Lumi");
            eventweight = lumi * Weight;
            weight *= eventweight * leptonScaleFactor * bTagScaleFactor *
                htDerivedScaleFactor * prefiringScaleFactor * puScaleFactor;
             weight = 1.0;
        }

        SliceData data{HT_trigger_pt30, GenHT, NJets_pt20, num_leps, JetsCA12, nsrCA12_21, nsrCA12_42, nsrCA12_43};
        my_histos.processCuts(data);

#define Fill(table, var) my_histos.fill(table, weight, var)

        Fill("EventCounter", eventCounter);

        if(num_leps==0 && nw_deep_tag >= 2 && reco_stop_mass > 0.0f ){
            Fill("RecoStopMass", reco_stop_mass);
            Fill("RecoStopMassImbalance", imbalance);
            Fill("RecoSBottomMass", reco_sbottom_mass);
            Fill("RecoSBottomMassImbalance", reco_sbottom_imbalance);
        }


        Fill("NJets_pt20", NJets_pt20);
        Fill("NJets_pt30", NJets_pt30);
        Fill("nCA12Jets", std::size(JetsCA12));
        Fill("nbjets_loose", NBJets_loose);
        // Fill("nwjets", nwjets);
        Fill("HT_pt30", HT_trigger_pt30);
        Fill("met", MET);
        Fill("mass_ratio", data.mass_ratio);
        Fill("DeepAK8TagW_medium_wp", nw_deep_tag);
        Fill("DeepAK8TagW_HasNearGenW", std::size(matched_w_jets));
        for(const auto& pair: matched_w_jets){
            Fill("DeepAK8TagW_MatchedEnergyRatio", Ws[pair.first].E()/JetsAK8[pair.second].E());
        }

        //Fill("MT2", MT2);
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
            return ROOT::Math::VectorUtil::Angle(pair[0].Vect(), pair[1].Vect());
            //return pair[0].Vect().Angle(pair[1].Vect());
        };

        for (std::size_t i = 0; i < std::size(JetsCA12); ++i) {
            Fill("CA12_NSubJettiness3", JetsCA12_NsubjettinessTau3[i]);
            Fill("CA12_NSubJettiness2", JetsCA12_NsubjettinessTau2[i]);
            Fill("CA12_NSubJettiness1", JetsCA12_NsubjettinessTau1[i]);
            Fill("CA12_NSubJettiness4", JetsCA12_NsubjettinessTau4[i]);
            Fill("CA12_NSubRatio42", nsrCA12_42[i]);
            Fill("CA12_NSubRatio43", nsrCA12_43[i]);
            Fill("CA12_NSubRatio21", nsrCA12_21[i]);
        }

        for (std::size_t i = 0; i < std::size(JetsAK8); ++i) {
            Fill("AK8_NSubJettiness3", JetsAK8_NsubjettinessTau3[i]);
            Fill("AK8_NSubJettiness2", JetsAK8_NsubjettinessTau2[i]);
            Fill("AK8_NSubJettiness1", JetsAK8_NsubjettinessTau1[i]);
            Fill("AK8_NSubRatio32", nsrAK8_32[i]);
            Fill("AK8_NSubRatio31", nsrAK8_31[i]);
            Fill("AK8_NSubRatio21", nsrAK8_21[i]);
        }

        Fill("NGoodLeptons", std::size(GoodLeptons));
        for(std::size_t i = 0 ; i < std::size(GoodLeptons); ++i){
            const auto& lepton =  GoodLeptons[i];
            const auto& v = lepton.second;
            Fill("l_Pt" , v.Pt() );
            Fill("l_Phi" , v.Phi());
            Fill("l_Eta" , v.Eta());
            Fill("l_Charge" , GoodLeptonsCharge[i]);
            Fill("l_MiniIso", GoodLeptonsMiniIso[i]);
            if(lepton.first == "e"){
                Fill("e_Pt" , v.Pt() );
                Fill("e_Phi" , v.Phi());
                Fill("e_Eta" , v.Eta());
                Fill("e_Charge" , GoodLeptonsCharge[i]);
                Fill("e_MiniIso", GoodLeptonsMiniIso[i]);
            } else if (lepton.first == "m") {
                Fill("m_Pt" , v.Pt() );
                Fill("m_Phi" , v.Phi());
                Fill("m_Eta" , v.Eta());
                Fill("m_Charge" , GoodLeptonsCharge[i]);
                Fill("m_MiniIso", GoodLeptonsMiniIso[i]);

            }
        }



        for (std::size_t i = 0; i < std::min(MAX_JETS, std::size(Jets)); ++i) {
            auto& jet = Jets[i];
            Fill("Jet_" + std::to_string(i) + "_E", jet.E());
            Fill("Jet_" + std::to_string(i) + "_P", jet.P());
            Fill("Jet_" + std::to_string(i) + "_Pt",jet.Pt());
            Fill("Jet_" + std::to_string(i) + "_Phi",jet.Phi());
            Fill("Jet_" + std::to_string(i) + "_Eta",jet.Eta());
        }

        if (std::size(Ws) == 2) Fill("Gen_W_Angle", computeAngle(Ws));

        for (std::size_t i = 0; i < std::size(JetsCA12); ++i) {
            Fill("CA12Pt", JetsCA12[i].Pt());
        }

    }
}

void Analyze2W::WriteHistos(TFile *outfile) {
    outfile->cd();
    for (const auto &p : my_histos) {
        p.second.Write();
    }
    for (const auto &p : my_2d_histos) p.second.Write();

    std::ofstream cutfile;
    cutfile.open("Output/cutflow_out.txt", std::ios_base::app);
    cutfile << filetag << "    ";
    cutfile << my_histos.makeCutflow("EventCounter",true)<< '\n';
}

