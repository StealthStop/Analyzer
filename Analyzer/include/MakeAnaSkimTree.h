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

   MiniTupleMaker *myAnaSkimTuple;
   TTree          *myTree;

   bool treeInit;
};

#endif
