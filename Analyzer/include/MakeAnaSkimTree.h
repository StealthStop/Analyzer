#ifndef MakeAnaSkimTree_h
#define MakeAnaSkimTree_h

#include <TTree.h>
#include <TH1D.h>

#include <map>
#include <string>

class NTupleReader;

class MiniTupleMaker;

class MakeAnaSkimTree{

public :

   std::shared_ptr<TH1D> eventCounter;

   MakeAnaSkimTree();
   ~MakeAnaSkimTree(){};

   void Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
   void InitHistos();
   void WriteHistos(TFile* outfile); 

   std::map<std::string, MiniTupleMaker*> myAnaSkimTuples;
   std::map<std::string, TTree*>          myTrees;

   std::map<std::string, bool> treeInits;
   std::vector<std::string> jecvars;
   std::vector<std::string> ttvars;
};

#endif
