#ifndef AnalyzeBTagSF_h
#define AnalyzeBTagSF_h

#include <TH1D.h>
#include <TFile.h>

#include <map>
#include <string>

class NTupleReader;

class AnalyzeBTagSF{

public :
   std::map<std::string, std::shared_ptr<TH1D>>  my_histos;

   AnalyzeBTagSF();
   ~AnalyzeBTagSF(){};

   void     Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
   virtual void     InitHistos();
   virtual void     WriteHistos(TFile* outfile);
};

#endif
