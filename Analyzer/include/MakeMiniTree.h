#ifndef MakeMiniTree_h
#define MakeMiniTree_h

#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>

#include <map>
#include <string>

class NTupleReader;

class MiniTupleMaker;

class MakeMiniTree{

public :
   std::map<std::string, std::shared_ptr<TH1D>>  my_histos;
   std::map<std::string, std::shared_ptr<TH2D>>  my_2d_histos;

   MakeMiniTree();
   ~MakeMiniTree(){};

   void     Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
   void     InitHistos();
   void     WriteHistos(TFile* outfile); 

   MiniTupleMaker *myMiniTuple;
   TTree          *myTree;
};

#endif
