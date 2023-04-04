#ifndef MakeQCDValTree_h
#define MakeQCDValTree_h

#include <TTree.h>
#include <TH1D.h>

#include <map>
#include <string>

class NTupleReader;

class MiniTupleMaker;

class MakeQCDValTree{

public :

   std::shared_ptr<TH1D> eventCounter;

   MakeQCDValTree();
   ~MakeQCDValTree(){};

   void Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
   void InitHistos();
   void WriteHistos(TFile* outfile); 

   MiniTupleMaker *myQCDValTuple;
   TTree          *myTree;

   bool treeInit;
};

#endif
