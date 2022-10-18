#ifndef AnalyzeGenStop_h
#define AnalyzeGenStop_h

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>

#include <map>
#include <string>

class NTupleReader;

class AnalyzeGenStop 
{
public:
    class TH1DInfo
    {
    public:
        std::string name;
        int nBins;
        double low;
        double high;        
    };

    class TH2DInfo
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

    std::map<std::string, std::shared_ptr<TH1D>>  my_histos;
    std::map<std::string, std::shared_ptr<TH2D>>  my_2d_histos;
    bool initHistos;
   
    AnalyzeGenStop();
    ~AnalyzeGenStop(){};
    
    void Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
    void InitHistos(const std::map<std::string, bool>& cutMap, const std::vector<TH1DInfo>& histInfos, const std::vector<TH2DInfo>& hist2DInfos);
    void WriteHistos(TFile* outfile);
    
};

#endif
