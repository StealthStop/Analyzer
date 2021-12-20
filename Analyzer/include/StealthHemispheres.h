#ifndef StealthHemispheres_H
#define StealthHemispheres_H

#include "TObject.h"
#include "TVector3.h"
#include "TFile.h"
#include "TFile.h"
#include "TEfficiency.h"
#include "TH1D.h"
#include "TH2D.h"

#include <map>

class NTupleReader;

class StealthHemispheres 
{
private:
    bool inithisto;
    
public:
    std::map<std::string, std::shared_ptr<TH1D>>  my_histos;
    std::map<std::string, std::shared_ptr<TH2D>>  my_2d_histos;
    std::map<std::string, std::shared_ptr<TEfficiency>>  my_efficiencies;

    void Loop(NTupleReader& tr, double weight, int maxevents = -1, bool isQuiet = false);
    void InitHistos(const std::map<std::string, bool>& cutmap);
    void WriteHistos(TFile* outfile);
    StealthHemispheres();
};

#endif
