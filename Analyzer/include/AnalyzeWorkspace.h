#ifndef AnalyzeWorkspace_h
#define AnalyzeWorkspace_h

#include <TH1D.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TProfile2D.h>
#include <TTree.h>
#include <TFile.h>

#include <map>
#include <string>

class NTupleReader;

class AnalyzeWorkspace 
{
private:
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

    class TH2DProfileInfo
    {
    public:
        std::string name;
        int nBinsX;
        double lowX;
        double highX;        
        int nBinsY;
        double lowY;
        double highY;        
        double lowZ;
        double highZ;
    };

    bool inithistos;

public:
    std::map<std::string, std::shared_ptr<TH1D>>  my_histos;
    std::map<std::string, std::shared_ptr<TH2D>>  my_2d_histos;
    std::map<std::string, std::shared_ptr<TEfficiency>>  my_efficiencies;    

    AnalyzeWorkspace();
    ~AnalyzeWorkspace(){};
    
    void Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
    void InitHistos(const std::map<std::string, bool>& cutMap, const std::vector<TH1DInfo>& histInfos);
    void WriteHistos(TFile* outfile);
    
};

#endif
