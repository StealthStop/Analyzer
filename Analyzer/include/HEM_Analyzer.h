#ifndef HEM_Analyzer_h
#define HEM_Analyzer_h

#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>

#include <map>
#include <string>

class NTupleReader;

class HEM_Analyzer 
{
public:
    std::map<std::string, std::shared_ptr<TH1D>>  my_histos;
    std::map<std::string, std::shared_ptr<TH2D>>  my_2d_histos;
    bool inithisto;   
 
    HEM_Analyzer();
    ~HEM_Analyzer(){};
    
    void Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
    void InitHistos(const std::map<std::string, bool>& cutmap);
    void WriteHistos(TFile* outfile);   
    
};

#endif
