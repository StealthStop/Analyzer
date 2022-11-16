#ifndef MakeMiniTree_h
#define MakeMiniTree_h

#include <TTree.h>
#include <TH1D.h>

#include <map>
#include <string>

class NTupleReader;

class MiniTupleMaker;

class MakeMiniTree{

public :

   std::shared_ptr<TH1D> eventCounter;

   MakeMiniTree();
   ~MakeMiniTree(){};

   void Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
   void InitHistos();
   void WriteHistos(TFile* outfile); 

   MiniTupleMaker *myMiniTuple;
   TTree          *myTree;
};

#endif
