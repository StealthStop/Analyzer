#ifndef Analyse2W_h
#define Analyse2W_h

#include <TEfficiency.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TCanvas.h>


#include <map>
#include <functional>
#include <vector>
#include <string>
#include <unordered_map>


template<typename T>
class Histogram{
    public:
        using HT =  std::shared_ptr<T>;

    private:
        HT histogram;

    public:
        std::vector<std::function<void(T*, TCanvas*)>> mods;
        bool log_y = false;
#define FDef(func) template<typename... Args> \
        void func(Args&&... args){ histogram->func(std::forward<Args>(args)...); }
#define FDefC(func) template<typename... Args> \
        void func(Args&&... args) const { histogram->func(std::forward<Args>(args)...); }
        FDef(Fill) FDef(FillRandom) FDef(SetTitle)
        FDefC(Draw) FDefC(Print) FDefC(Write)
#undef FDef

        /*
        template<typename... Args>
            Histogram(Args&& ... args): histogram{std::make_shared<T>(std::forward<Args>(args)...)} {}

            */
            Histogram(const char* c1, const char* c2, double v1, double v2, double v3): histogram{std::make_shared<T>(c1,c2,v1,v2,v3)} {}
            Histogram(){}

        void DrawImage(const std::string& outputDir=".", const std::string type = "pdf") const {
            if(histogram == nullptr) return;
            auto c = std::make_unique<TCanvas>();
            if(log_y) c->SetLogy();
            Draw();
            c->Update();
            for(auto& f: mods) f(histogram.get(), c.get());
            c->Modified();
            c->Print((outputDir + "/" +  histogram->GetName() + '.' + "pdf").c_str());
        }
};



class NTupleReader;

class Analyze2W {
private:
   public:
    std::unordered_map<std::string, Histogram<TH1D>> my_histos;
    std::unordered_map<std::string, Histogram<TH2D>> my_2d_histos;
    std::unordered_map<std::string, Histogram<TEfficiency>> my_efficiencies;

    Analyze2W();
    ~Analyze2W(){};

    void Loop(NTupleReader& tr, double weight, int maxevents = -1,
              bool isQuiet = false);
    void InitHistos();
    void WriteHistos(TFile* outfile);
};


#endif
