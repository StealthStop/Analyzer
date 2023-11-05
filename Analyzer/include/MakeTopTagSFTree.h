#ifndef MakeTopTagSFTree_h
#define MakeTopTagSFTree_h

#include <TTree.h>
#include <TH1D.h>

#include <map>
#include <string>

class NTupleReader;

class MiniTupleMaker;

class MakeTopTagSFTree{

public :

   std::shared_ptr<TH1D> eventCounter;

   MakeTopTagSFTree();
   ~MakeTopTagSFTree(){};

   void Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
   void WriteHistos(TFile* outfile); 

   std::map<std::string, MiniTupleMaker*> myTopTagSFTuple;
   std::map<std::string, TTree*>          myTree;

   std::map<std::string, bool> treeInit;

   std::vector<std::string> jecvars;
};

#endif
