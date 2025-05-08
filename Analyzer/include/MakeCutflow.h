#ifndef MakeCutflow_h
#define MakeCutflow_h

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TEfficiency.h>

#include <map>
#include <string>

class NTupleReader;

class MakeCutflow 
{
public:
    std::map<std::string, std::shared_ptr<TH1D>>  my_histos;
    
    MakeCutflow();
    ~MakeCutflow(){};
    
    void Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
    void InitHistos();
    void WriteHistos(TFile* outfile);
    
};

#endif
