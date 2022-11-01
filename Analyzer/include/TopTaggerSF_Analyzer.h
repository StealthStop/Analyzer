#ifndef TopTaggerSF_Analyzer_h
#define TopTaggerSF_Analyzer_h

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>

#include <map>
#include <string>

class NTupleReader;

class TopTaggerSF_Analyzer
{

private:
    class TH1Dinfo
    {
    public:
        std::string name;
        int nBins;
        double low;
        double high;        
    };

    class TH2Dinfo
    {
    public:
        std::string name;
        int nBinsX;
        double lowX;
        double highX;        
        int nBinsY;
        double lowY;
        double highY;        
    };

public:
    std::map<std::string, std::shared_ptr<TH1D>> my_1D_histos;
    std::map<std::string, std::shared_ptr<TH2D>> my_2D_histos;

    std::vector<TH1Dinfo> hist1Dinfos;
    std::vector<TH2Dinfo> hist2Dinfos;

    std::vector<std::string> njets;
    std::vector<std::string> tags;
    std::vector<std::string> my_var_suffix;

    bool initHistos;

    TopTaggerSF_Analyzer();
    ~TopTaggerSF_Analyzer(){};

    void Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
    void InitHistos(const std::map<std::string, bool>& cutMap);
    void WriteHistos(TFile* outfile);

    Double_t vecToDouble_t(const std::vector<double>& input);
};

#endif
