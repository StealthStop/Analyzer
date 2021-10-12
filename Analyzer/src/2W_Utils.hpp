#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#include <TPaveStats.h>
#pragma GCC diagnostic pop
#include "Analyzer/Analyzer/include/Analyze2W.h"

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
#include <type_traits>

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


std::ostream& operator<<(std::ostream &os, const TLorentzVector &v) {
    os << "(" << v.P() << ", " << v.Eta() << ", " << v.Phi() << "," << v.E() << ")";
    return os;
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

template <typename T> class ProgressBar {
    private:
        float barWidth = 40;

    public:
        T max_items = 1, current_items = 0;
        std::string text = "";
        ProgressBar(T m, std::string s) : max_items{m}, text{std::move(s)} {}
        ProgressBar<T> &operator++() {
            T next = current_items + 1;
            current_items = (next <= max_items) ? next : max_items;
            return *this;
        }
        void display(const std::string& item = "") {
            std::cout << "[";
            float progress = float(current_items) / float(max_items);
            int pos = barWidth * progress;
            for (int i = 0; i < barWidth; ++i) {
                if (i < pos)
                    std::cout << "=";
                else if (i == pos)
                    std::cout << ">";
                else
                    std::cout << " ";
            }
            std::cout << "] "
                << int(progress * 100.0) << "% "
                << text << ((item == "" ) ? std::to_string(current_items) : item)
                << "\r";
            std::cout.flush();
        }
};
