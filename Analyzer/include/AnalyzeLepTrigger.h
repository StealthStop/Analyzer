#ifndef AnalyzeLepTrigger_h
#define AnalyzeLepTrigger_h

#include <TH1D.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TTree.h>
#include <TLorentzVector.h>

#include <map>
#include <string>

class NTupleReader;

class AnalyzeLepTrigger 
{
public:
    std::map<std::string, std::shared_ptr<TH1D>>  my_histos;
    std::map<std::string, std::shared_ptr<TH2D>>  my_2d_histos;
    std::map<std::string, std::shared_ptr<TEfficiency>>  my_efficiencies;
    
    AnalyzeLepTrigger();
    ~AnalyzeLepTrigger(){};
    
    void Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
    void InitHistos();
    void WriteHistos(TFile* outfile);
    void printTriggerList( const std::vector<std::string>& TriggerNames );
    
    bool containsGoodLepton( const std::vector<TLorentzVector>& leptons, const std::vector<bool>& goodLeptons, double ptThreshold, double etaSelection);
    int goodLeptonIndex( const std::vector<TLorentzVector>& leptons, const std::vector<bool>& goodLeptons );
    void fillHistos( const std::map<std::string, bool>& cutMap, bool passLeptonTriggers, const TLorentzVector& lepton, double theWeight );
};

#endif
