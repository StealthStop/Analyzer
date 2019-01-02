#ifndef AnalyzeNjetsMinusOneCSFillDijetHists_h
#define AnalyzeNjetsMinusOneCSFillDijetHists_h

#include <TH1D.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TTree.h>
#include <TFile.h>

#include <map>
#include <string>

class NTupleReader;

class AnalyzeNjetsMinusOneCSFillDijetHists 
{
public:
    std::map<std::string, std::shared_ptr<TH1D>>  my_histos;
    std::map<std::string, std::shared_ptr<TH2D>>  my_2d_histos;
    bool initHistos;
    
    AnalyzeNjetsMinusOneCSFillDijetHists();
    ~AnalyzeNjetsMinusOneCSFillDijetHists(){};
    
    void Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
    void InitHistos();
    void WriteHistos(TFile* outfile);
private:
    
};

#endif
