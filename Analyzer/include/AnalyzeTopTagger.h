#ifndef AnalyzeTopTagger_h
#define AnalyzeTopTagger_h

#include <TH1D.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TTree.h>

#include <map>
#include <string>
#include "TopTaggerTools/Tools/include/HistoContainer.h"

class NTupleReader;

class AnalyzeTopTagger
{
private:
   HistoContainer<NTupleReader> hists_old, histNjet6_old, histNjet7_old, histNjet8_old, histNjet9_old,  histNjet10_old, histNjet11_old, histNjet12_old, histNjet12inc_old, 
                                hists_new, histNjet7_new, histNjet8_new, histNjet9_new, histNjet10_new, histNjet11_new, histNjet12_new, histNjet12inc_new; 
public:
   std::map<std::string, std::shared_ptr<TH1D>>  my_histos;
   std::map<std::string, std::shared_ptr<TH2D>>  my_2d_histos;
   std::map<std::string, std::shared_ptr<TEfficiency>>  my_efficiencies;
    
   AnalyzeTopTagger();
   ~AnalyzeTopTagger(){};

   void Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
   void InitHistos();
   void WriteHistos(TFile* outfile);

};

#endif
