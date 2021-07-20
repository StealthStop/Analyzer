#ifndef MakeNNVariables_h
#define MakeNNVariables_h

#include <TH1D.h>
#include <TTree.h>

#include <map>
#include <string>

class NTupleReader;
class MiniTupleMaker;

class MakeNNVariables{

public :
   std::map<std::string, std::shared_ptr<TH1D>>  my_histos;

   std::vector<std::string> my_var_suffix;

   MakeNNVariables();
   ~MakeNNVariables(){};

   void     Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
   void     InitHistos();
   void     WriteHistos(TFile* outfile); 

   std::map<std::string, std::map<std::string, std::map<std::string, MiniTupleMaker*> > > myMiniTuple;
   std::map<std::string, std::map<std::string, std::map<std::string, TTree*> > > myTree;

};

#endif
