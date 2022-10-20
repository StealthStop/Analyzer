#ifndef ResolvedTopTagger_Analyzer_h
#define ResolvedTopTagger_Analyzer_h

#include <TH1D.h>
#include <TH2D.h>

#include <map>
#include <string>
#include "TopTaggerTools/Tools/include/HistoContainer.h"

class NTupleReader;

class ResolvedTopTagger_Analyzer
{
private:
   HistoContainer<NTupleReader> hists, histNjet7, histNjet8, histNjet9, histNjet10, histNjet11, histNjet12, histNjet12inc; 
public:
   std::map<std::string, std::shared_ptr<TH1D>>  my_histos;
   std::map<std::string, std::shared_ptr<TH2D>>  my_2d_histos;
    
   ResolvedTopTagger_Analyzer();
   ~ResolvedTopTagger_Analyzer(){};

   void Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
   void InitHistos();
   void WriteHistos(TFile* outfile);

};

#endif
