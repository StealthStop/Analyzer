#ifndef CalculateTopTagSF_h
#define CalculateTopTagSF_h

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>

#include <map>
#include <string>

class NTupleReader;

class CalculateTopTagSF{

public :
   std::map<std::string, std::shared_ptr<TH1D>>  my_histos;
   std::map<std::string, std::shared_ptr<TH2D>>  my_2d_histos;
   bool initHistos;

   CalculateTopTagSF();
   ~CalculateTopTagSF(){};

   void     Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
   virtual void     InitHistos(const std::string& histoFileTag);
   virtual void     WriteHistos(TFile* outfile);
};

#endif
