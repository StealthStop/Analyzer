#define AnalyzeXsec_cxx
#include "Analyzer/Analyzer/include/AnalyzeXsec.h"
#include "Framework/Framework/include/Utility.h"
#include "NTupleReader/include/NTupleReader.h"

#include <iostream>

AnalyzeXsec::AnalyzeXsec()
{
}

void AnalyzeXsec::WriteHistos(TFile* outfile)
{
    // Do not need to do anything here
}

//Put everything you want to do per event here.
void AnalyzeXsec::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    while( tr.getNextEvent() )
    {
        if (maxevents != -1 && tr.getEvtNum() >= maxevents)
            break;        

        const auto& filetag         = tr.getVar<std::string>("filetag");     
        const auto& TreeMakerWeight = tr.getVar<float>("Weight");
        const auto& SampleSetWeight = tr.getVar<double>("weight");

        if (tr.getEvtNum() % 1000 == 0)
            printf("  Event %i\n", tr.getEvtNum() );

        std::cout << "XSEC INFO: " << filetag << ", TM:" << TreeMakerWeight << ", SS:" << SampleSetWeight << std::endl;
        break;
    }
}
