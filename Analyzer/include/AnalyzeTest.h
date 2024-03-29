#ifndef AnalyzeTest_h
#define AnalyzeTest_h

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TEfficiency.h>

#include <map>
#include <string>

class NTupleReader;

class AnalyzeTest 
{
public:
    std::map<std::string, std::shared_ptr<TH1D>>  my_histos;
    std::map<std::string, std::shared_ptr<TH2D>>  my_2d_histos;
    std::map<std::string, std::shared_ptr<TEfficiency>>  my_efficiencies;
    
    AnalyzeTest();
    ~AnalyzeTest(){};
    
    void Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
    void InitHistos();
    void WriteHistos(TFile* outfile);
    
};

#endif
