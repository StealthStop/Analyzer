#ifndef AnalyzeXsec_h
#define AnalyzeXsec_h

#include "TFile.h"

class NTupleReader;

class AnalyzeXsec 
{
public:
    AnalyzeXsec();
    ~AnalyzeXsec(){};
    
    void Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
    void WriteHistos(TFile* outfile);
};

#endif
