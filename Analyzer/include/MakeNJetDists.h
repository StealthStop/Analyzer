#ifndef MakeNJetDists_h
#define MakeNJetDists_h

#include <TH1D.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TTree.h>

#include <map>
#include <string>

#include "Analyzer/Analyzer/include/Histo.h"
class NTupleReader;

class MakeNJetDists
{
public:
    std::vector<std::unique_ptr<Histo_Base>> my_Histos;
    std::vector<std::string> myVarSuffixs;
    
    MakeNJetDists();
    ~MakeNJetDists(){};
    
    void Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
    void InitHistos();
    void WriteHistos(TFile* outfile); 
};

#endif
