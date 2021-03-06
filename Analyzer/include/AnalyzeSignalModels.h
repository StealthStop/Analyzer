#ifndef AnalyzeSignalModels_h
#define AnalyzeSignalModels_h

#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TEfficiency.h>
#include <TTree.h>

#include <map>
#include <string>

class NTupleReader;

class AnalyzeSignalModels 
{
public:
    std::map<std::string, std::shared_ptr<TH1D>>  my_histos;
    std::map<std::string, std::shared_ptr<TH2D>>  my_2d_histos;
    std::map<std::string, std::shared_ptr<TH3D>>  my_3d_histos;
    std::map<std::string, std::shared_ptr<TEfficiency>>  my_efficiencies;
    
    AnalyzeSignalModels();
    ~AnalyzeSignalModels(){};
    
    void Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
    void InitHistos();
    void WriteHistos(TFile* outfile);
    
};

#endif

