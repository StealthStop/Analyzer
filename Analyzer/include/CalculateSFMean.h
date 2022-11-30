#ifndef CalculateSFMean_h
#define CalculateSFMean_h

#include "Analyzer/Analyzer/include/AnalyzeBase.h"
#include "Analyzer/Analyzer/include/Histo.h"
class NTupleReader;

class CalculateSFMean : public AnalyzeBase
{
public:
    void InitHistos(const std::string& filetag)
    {
        TH1::SetDefaultSumw2();
        TH2::SetDefaultSumw2();

        my_Histos.emplace_back(new Histo1D("EventCounter", 2,  -1.1, 1.1, "eventCounter", {}, {}));
        my_Histos.emplace_back(new Histo1D(filetag+"_sclUp",        100, 0, 5, "scaleWeightUpUncor",            {"passMadHT"}, {"Lumi", "Weight"}));
        my_Histos.emplace_back(new Histo1D(filetag+"_sclDown",      100, 0, 5, "scaleWeightDownUncor",          {"passMadHT"}, {"Lumi", "Weight"}));
        my_Histos.emplace_back(new Histo1D(filetag+"_pdf_Up",       100, 0, 5, "PDFweightUpUncor",              {"passMadHT"}, {"Lumi", "Weight"}));
        my_Histos.emplace_back(new Histo1D(filetag+"_pdf_Down",     100, 0, 5, "PDFweightDownUncor",            {"passMadHT"}, {"Lumi", "Weight"}));    
        my_Histos.emplace_back(new Histo1D(filetag+"_pu",           100, 0, 5, "puWeightUnCorr",                {"passMadHT"}, {"Lumi", "Weight"}));    
        my_Histos.emplace_back(new Histo1D(filetag+"_pu_Up",        100, 0, 5, "puSysUpUnCorr",                 {"passMadHT"}, {"Lumi", "Weight"}));    
        my_Histos.emplace_back(new Histo1D(filetag+"_pu_Down",      100, 0, 5, "puSysDownUnCorr",               {"passMadHT"}, {"Lumi", "Weight"}));    
    }

    void Loop(NTupleReader& tr, double, int maxevents, bool)
    {
        const auto& runtype = tr.getVar<std::string>("runtype");
        const auto& filetag = tr.getVar<std::string>("filetag");

        //-------------------------------------
        //-- Initialize histograms to be filled
        //-------------------------------------
        InitHistos(filetag);

        while( tr.getNextEvent() )
        {
            //------------------------------------
            //-- Print Event Number
            //------------------------------------
            const bool breakLoop = printEventNum(maxevents, tr.getEvtNum());
            if(breakLoop) break;

            //-----------------------------------
            //-- Code unique to this Analyzer
            //-----------------------------------
            if( runtype != "MC" ) 
            {
                std::cerr<<"Please run over an MC file since these scale factors should not be applied to data!!"<<std::endl;
                break;
            }
               
            //-----------------------------------
            //-- Fill Histograms Below
            //-----------------------------------
            Fill(tr);

        }//END of while tr.getNextEvent loop   
    }//END of function
};

#endif
