#ifndef HadTriggers_Analyzer_h
#define HadTriggers_Analyzer_h

#include "Framework/Framework/include/Utility.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>

#include <map>
#include <string>

class NTupleReader;

class HadTriggers_Analyzer 
{
public:
    std::map<std::string, std::shared_ptr<TH1D>>  my_histos;
    std::map<std::string, std::shared_ptr<TH2D>>  my_2d_histos;
    
    HadTriggers_Analyzer();
    ~HadTriggers_Analyzer(){};
    
    void Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
    void InitHistos();
    void WriteHistos(TFile* outfile);
    void printTriggerList( const std::vector<std::string>& TriggerNames );
    
    bool containsGoodHadron( const std::vector<utility::LorentzVector>& hadrons, const std::vector<bool>& goodHadrons, double ptThreshold, double etaSelection);
    void fillHistos( const std::map<std::string, bool>& cutMap, bool passTriggerAllHad, double pt, double HT, int njet, int nbjet, double weight );
};

#endif
