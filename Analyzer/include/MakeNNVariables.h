#ifndef MakeNNVariables_h
#define MakeNNVariables_h

#include <TTree.h>

#include <map>
#include <string>

class NTupleReader;
class MiniTupleMaker;

class MakeNNVariables{

public :

   std::vector<std::string> my_var_suffix;
   std::vector<std::string> my_channels;
   std::vector<std::string> my_splits;

   std::map<std::string, std::map<std::string, std::map<std::string, int> > > my_counts;

   MakeNNVariables();
   ~MakeNNVariables(){};

   void Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
   void WriteHistos(TFile* outfile); 

   std::map<std::string, std::map<std::string, std::map<std::string, MiniTupleMaker*> > > myMiniTuple;
   std::map<std::string, std::map<std::string, std::map<std::string, TTree*> > > myTree;

   std::map<std::string, bool> treeInit;

};

#endif
