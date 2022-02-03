#ifndef TopTaggerSF_Analyzer_h
#define TopTaggerSF_Analyzer_h

#include <TH1D.h>
#include <TTree.h>

#include <map>
#include <string>
#include "TopTaggerTools/Tools/include/HistoContainer.h"

class NTupleReader;

class TopTaggerSF_Analyzer
{

private:

   HistoContainer<NTupleReader> ttHistsNjetIncl,  ttHistsNjetIncl_Joe,  ttHistsNjet7,  ttHistsNjet8,  ttHistsNjet9,  ttHistsNjet10,  ttHistsNjet11,  ttHistsNjet12incl;
   HistoContainer<NTupleReader> qcdHistsNjetIncl, qcdHistsNjetIncl_Joe, qcdHistsNjet7, qcdHistsNjet8, qcdHistsNjet9, qcdHistsNjet10, qcdHistsNjet11, qcdHistsNjet12incl;

   class TH1DInfo
   {
   public:
       std::string name;
       int nBins;
       double low;
       double high;        
   };

   std::vector<TH1DInfo> histInfos;
   bool initHistos;

public:

   std::map<std::string, std::shared_ptr<TH1D>> my_histos;
    
   TopTaggerSF_Analyzer();
   ~TopTaggerSF_Analyzer(){};

   void Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
   void InitHistos(const std::map<std::string, bool>& cutMap);
   void WriteHistos(TFile* outfile);

};

#endif
