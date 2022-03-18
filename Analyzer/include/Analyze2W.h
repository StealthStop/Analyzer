#ifndef Analyse2W_h
#define Analyse2W_h

#include <TCanvas.h>
#include <TEfficiency.h>
#include <TError.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>

#include <algorithm>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>

/*
constexpr inline auto string_hash(const char *s) {
    unsigned long long hash{}, c{};
    for (auto p = s; *p; ++p, ++c) hash += *p << c;
    return hash;
}

constexpr inline auto operator"" _sh(const char *s, size_t) {
        return string_hash(s);
}
*/


#define DO_DEBUG false
#define DEBUG(x)                                                               \
  do {                                                                         \
    if (DO_DEBUG) {                                                            \
      std::cout << __FILE__ << ":" << __LINE__ << ": " << x << std::endl;      \
    }                                                                          \
  } while (false)

template <typename T> class Histogram {
public:
  using HT = std::shared_ptr<T>;

private:
  HT histogram;

public:
  std::vector<std::function<void(T *, TCanvas *)>> mods;
  bool log_y = false;
#define FDef(func)                                                             \
  template <typename... Args> void func(Args &&... args) {                     \
    histogram->func(std::forward<Args>(args)...);                              \
  }
#define FDefC(func)                                                            \
  template <typename... Args> void func(Args &&... args) const {               \
    histogram->func(std::forward<Args>(args)...);                              \
  }
  FDef(Fill) FDef(FillRandom) FDef(SetTitle) FDefC(Draw) FDefC(Print)
      FDefC(Write)
#undef FDef
#undef FDefC
          template <typename... Args>
          double GetEntries(Args &&... args) const {
    return histogram->GetEntries(std::forward<Args>(args)...);
  }
  template <typename... Args> double Integral(Args &&... args) const {
    return histogram->Integral(std::forward<Args>(args)...);
  }
  Histogram(const char *c1, const char *c2, double v1, double v2, double v3)
      : histogram{std::make_shared<T>(c1, c2, v1, v2, v3)} {}
  Histogram() {}

  void DrawImage(std::string outputDir = ".",
                 const std::string type = "pdf") const {
    gErrorIgnoreLevel = // kPrint,
                        // kInfo,
        kWarning        //, kError, kBreak, kSysError, kFatal
        ;
    if (histogram == nullptr)
      return;
    auto c = std::make_unique<TCanvas>();
    if (log_y)
      c->SetLogy();
    Draw();
    c->Update();
    for (auto &f : mods)
      f(histogram.get(), c.get());
    c->Modified();

    struct stat info;

    outputDir = std::string("Output/")+outputDir;
    if( stat( outputDir.c_str(), &info ) != 0 ) system((std::string("mkdir -p ")+outputDir).c_str());
    c->Print((outputDir + "/" + histogram->GetName() + '.' + type).c_str());
  }
};

struct SliceData;

class Cut {
public:
  Cut(std::string name, bool destructive,
      std::vector<std::string> possible_values)
      : name{name}, destructive{destructive}, possible_values{
                                                  std::move(possible_values)} {}
  std::string name;
  bool destructive = false;
  virtual void calculate(const SliceData &data) = 0;
  const std::vector<std::string> possible_values;

  bool passed = false;
  std::string value = "";
  int pass_idx = 0;

  virtual ~Cut(){};
  friend std::ostream &operator<<(std::ostream &os, Cut *&c) {
    return os << c->name;
  }
};

class Generator {
public:
  Generator(std::string name) : name{name} {}
  std::string name;
  virtual void calculate(SliceData &data) = 0;
  virtual ~Generator(){};
};

/*
template <typename T> typename T::value_type getAllNames(T vals) {
  typename T::value_type ret;
  ret.reserve(
      std::accumulate(vals.begin(), vals.end(), 1, [](auto x, const auto &y) {
        return x * (std::size(y) + 1);
      }));
  for (int i = 0; i < std::size(vals); ++i) {
    int s = ret.size();
    for (int j = 0; j < s; ++j) {
      for (const auto &y : vals[i]) {
        ret.push_back(ret[j] + "_" + y);
      }
    }
    for (const auto &y : vals[i])
      ret.push_back(y);
  }
  return ret;
}
*/

template <typename T> typename T::value_type getAllNames(T vals) {
  typename T::value_type ret;
  std::vector<typename T::value_type::iterator> it;
  const std::size_t K = std::size(vals);
  for(auto& x : vals) it.push_back(x.begin());
  while (it[0] != vals[0].end()) {
        // process the pointed-to elements
        std::string temp = "";
        for(std::size_t z=0;z<std::size(it);++z){
            if (z) temp += "_";
            temp += *it[z];
        }
        ret.push_back(temp);
        // the following increments the "odometer" by 1
        ++it[K-1];
        for (int i = K-1; (i > 0) && (it[i] == vals[i].end()); --i) {
                it[i] = vals[i].begin();
            ++it[i-1];
            }
        }
  return ret;
}
template <typename T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &v) {
  bool first = true;
  for (const T &val : v) {
    if (!first) {
      os << ", ";
    }
    os << val;
    first = false;
  }
  return os;
}

struct Chain {
  Chain(auto &&name, auto &&c) : name{name}, chain{c} {}
  Chain(auto &&x) : Chain("", x) {}
  std::string name;
  std::vector<Cut *> chain;
  auto begin() { return chain.begin(); }
  auto end() { return chain.end(); }
  auto begin() const { return chain.cbegin(); }
  auto end() const { return chain.cend(); }
  std::size_t size() const { return std::size(chain); }
  std::string getPassedName() const {
    std::string current_name = "";
    for (const auto &cut : chain) {
      current_name += "_";
      current_name += cut->possible_values[cut->pass_idx];
    }
    return current_name;
  }
};

template<typename T>
std::string to_string (const std::vector<T>& v){
        std::string ret;
        std::size_t len = std::size(v);
        for(std::size_t i = 0 ; i< len; ++i){
            if constexpr (std::is_same_v<T,std::string>){
           ret += v[i];
           } else{
           ret += std::to_string(v[i]);
           }
           if(i < len-1) ret += ", ";
        }
        return ret;
    };

template <typename H> class HistogramManager {
public:
  std::vector<std::unique_ptr<Cut>> cuts;
  std::vector<Chain> chains;
  std::unordered_map<std::string, Histogram<H>> my_histos;
  std::vector<std::unique_ptr<Generator>> generators;

  void constructChains(const std::vector<std::vector<std::string>> &vals) {
    DEBUG("Constructing Chains");
    chains.clear();
    for (const auto &chain : vals) {
      std::vector<Cut *> new_chain;
      for (const auto &name : chain) {
        auto found = std::find_if(cuts.begin(), cuts.end(),
                                  [&](auto &x) { return x->name == name; });
        new_chain.push_back(found->get());
      }
      DEBUG(new_chain);
      chains.push_back(new_chain);
    }
  }

  std::vector<std::vector<std::string>>
  computePossibleFromChain(const Chain &c) {
    std::vector<std::vector<std::string>> ret;
    for (const auto &cut : c) {
      std::vector<std::string> possible_values = cut->possible_values;
      if (std::size(possible_values) && !cut->destructive) {
        ret.push_back(possible_values);
      }
    }
    return ret;
  }

  template <typename... Args>
  void createNewHistogram(std::string name, Args &&... args) {
    createAndAdd(name, std::forward<Args>(args)...);
    for (const auto &chain : chains) {
      std::vector<std::string> cuts =
          getAllNames(computePossibleFromChain(chain));
      for (const std::string &cut : cuts) {
        std::string newname = name + "_";
        newname += cut;
        createAndAdd(newname, std::forward<Args>(args)...);
      }
    }
  }

  void createAndAdd(std::string name, int v1, int v2, int v3,
                    std::string xlabel = "", std::string ylabel = "",
                    bool logaxis = false) {
    std::string basename = /*std::string("h_") + */ name;
    std::string newname =
        basename + ';' + std::move(xlabel) + ';' + std::move(ylabel);
    const char *x = basename.c_str();
    auto ret = Histogram<H>(x, x, v1, v2, v3);
    ret.log_y = logaxis;
    ret.SetTitle(newname.c_str());
    my_histos.emplace(basename, ret);
    DEBUG("Added new histogram with name " << basename);
  }

  auto &operator[](const std::string &s) {
    DEBUG("Attempting to fetch histogram " << s);
    return my_histos[s];
  }

  void processCuts(SliceData &d) {
    for (auto &g : generators)
      g->calculate(d);
    for (auto &cut : cuts)
      cut->calculate(d);
  }

  template <typename... T>
  void fill(const std::string &name, double weight, T... fill) {
    // DEBUG("Starting to fill variable" << name);
    std::vector<std::string> valid_names;
    for (const auto &chain : chains) {
      std::string current_name = "";
      bool fail = false;
      for (const auto &cut : chain) {
        if (cut->destructive && !cut->passed) {
          fail = true;
          break;
        } else if (!cut->destructive) {
          current_name += "_";
          current_name += cut->value;
        }
      }
      if (!fail)
        valid_names.push_back(current_name);
    }
    valid_names.push_back("");
    for (const std::string &s : valid_names) {
      // DEBUG("Attempting to fill histogram with name" << name+s << " with
      // value " << (fill << ...));
      my_histos.at(name + s).Fill(fill..., weight);
    }
  }
  std::string makeCutflow(const std::string &name, bool do_entries = false) {
    std::vector<std::string> names = {name};
    auto getval = [do_entries](auto &&x)->double {
      if (do_entries) {
        return x.GetEntries();
      } else {
        return x.Integral();
      }
    };
    std::vector<double> vals = {getval(my_histos[name])};
    for (const auto &chain : chains) {
      std::string hname = name + chain.getPassedName();
      names.push_back(hname);
      vals.push_back(getval(my_histos[hname]));
    }
    return  to_string(names) + "\n" + to_string(vals);
  };
  auto begin() { return my_histos.begin(); }
  auto end() { return my_histos.end(); }
  auto cbegin() const { return my_histos.cbegin(); }
  auto cend() const { return my_histos.cend(); }
  auto size() const { return my_histos.size(); }
  auto printHistos() {
    for (const auto &x : my_histos) {
      std::cout << x.first << ", ";
    }
    std::cout << "\n";
  }
};

class NTupleReader;

class Analyze2W {
private:
    std::string filetag;
public:
  //    std::unordered_map<std::string, Histogram<TH1D>> my_histos;
  //    std::unordered_map<std::string, Histogram<TH2D>> my_2d_histos;

  HistogramManager<TH1D> my_histos;
  HistogramManager<TH2D> my_2d_histos;

  Analyze2W();
  ~Analyze2W(){};

  void Loop(NTupleReader &tr, double weight, int maxevents = -1,
            bool isQuiet = false);
  void InitHistos();
  void WriteHistos(TFile *outfile);
};

#endif
