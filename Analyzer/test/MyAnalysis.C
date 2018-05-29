#include "SusyAnaTools/Tools/samples.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "TopTagger/CfgParser/include/TTException.h"

#include "Analyzer/Analyzer/include/AnalyzeBackground.h"
#include "Analyzer/Analyzer/include/AnalyzeTopTagger.h"
#include "Analyzer/Analyzer/include/AnalyzeEventSelection.h"
#include "Analyzer/Analyzer/include/AnalyzeEventShape.h"
#include "Analyzer/Analyzer/include/Analyze0Lep.h"
#include "Analyzer/Analyzer/include/AnalyzeStealthTopTagger.h"
#include "Analyzer/Analyzer/include/AnalyzeSHuHdEventSelection.h"
#include "Analyzer/Analyzer/include/AnalyzeDataMC.h"
#include "Analyzer/Analyzer/include/AnalyzeTest.h"
#include "Analyzer/Analyzer/include/AnalyzeHadTrigger.h"
#include "Analyzer/Analyzer/include/AnalyzeTemp.h"
#include "Analyzer/Analyzer/include/AnalyzeCleanJet.h"
#include "SusyAnaTools/Tools/NTupleReader.h"

#include "Framework/Framework/include/RunTopTagger.h"
#include "Framework/Framework/include/RunFisher.h"
#include "Framework/Framework/include/Muon.h"
#include "Framework/Framework/include/Electron.h"
#include "Framework/Framework/include/Jet.h"
#include "Framework/Framework/include/BJet.h"
#include "Framework/Framework/include/CommonVariables.h"
#include "Framework/Framework/include/Baseline.h"

#include "TH1D.h"
#include "TFile.h"
#include "TChain.h"

#include <iostream>
#include <getopt.h>
#include <string>

template<typename Analyze> void run(std::set<AnaSamples::FileSummary> vvf, 
                                    int startFile, int nFiles, int maxEvts, 
                                    bool isSkim, TFile* outfile)
{
    std::cout << "Initializing..." << std::endl;
    Analyze a;
    for(const AnaSamples::FileSummary& file : vvf)
    {
        // Define what is needed per sample set
        std::cout << "Running over sample " << file.tag << std::endl;
        TChain* ch = new TChain( (file.treePath).c_str() );
        file.addFilesToChain(ch, startFile, nFiles);
        NTupleReader tr(ch);
        double weight = file.getWeight(); // not used currently
        std::string runtype = (file.tag.find("Data") != std::string::npos) ? "Data" : "MC";
        std::cout << "Starting loop (in run)" << std::endl;
        printf( "runtype: %s fileWeight: %f nFiles: %i startFile: %i maxEvts: %i \n",runtype.c_str(),weight,nFiles,startFile,maxEvts ); fflush( stdout );
        tr.registerDerivedVar<std::string>("runtype",runtype);
        tr.registerDerivedVar<std::string>("filetag",file.tag);
        tr.registerDerivedVar<double>("etaCut",2.4);
        tr.registerDerivedVar<bool>("blind",true);

        // Define classes/functions that add variables on the fly
        std::shared_ptr<RunTopTagger> rtt;
        if ( !isSkim ) rtt = std::make_shared<RunTopTagger>();
        RunFisher runFisher;
        //RunFisher runFisher("test");
        Muon muon;
        Electron electron;
        Jet jet;
        BJet bjet;
        CommonVariables commonVariables;
        Baseline baseline;

        // Register classes/functions that add variables on the fly
        if ( !isSkim ) tr.registerFunction( std::move(*rtt) );
        tr.registerFunction( std::move(muon) );
        tr.registerFunction( std::move(electron) );
        tr.registerFunction( std::move(jet) );
        tr.registerFunction( std::move(runFisher) );
        tr.registerFunction( std::move(bjet) );
        tr.registerFunction( std::move(commonVariables) );
        tr.registerFunction( std::move(baseline) );

        // Loop over all of the events and fill histos
        a.Loop(tr, weight, maxEvts);

        // Cleaning up dynamic memory
        delete ch;
            
    }
    std::cout << "Writing histograms..." << std::endl;
    a.WriteHistos(outfile);
}

std::set<AnaSamples::FileSummary> setFS(const std::string& dataSets, const bool& isCondor)
{
    AnaSamples::SampleSet        ss("sampleSets.cfg", isCondor);
    AnaSamples::SampleCollection sc("sampleCollections.cfg", ss);

    std::map<std::string, std::vector<AnaSamples::FileSummary>> fileMap;
    if(ss[dataSets] != ss.null())
    {
        fileMap[dataSets] = {ss[dataSets]};
        for(const auto& colls : ss[dataSets].getCollections())
        {
            fileMap[colls] = {ss[dataSets]};
        }
    }
    else if(sc[dataSets] != sc.null())
    {
        fileMap[dataSets] = {sc[dataSets]};
        int i = 0;
        for(const auto& fs : sc[dataSets])
        {
            fileMap[sc.getSampleLabels(dataSets)[i++]].push_back(fs);
        }
    }
    std::set<AnaSamples::FileSummary> vvf;
    for(auto& fsVec : fileMap) for(auto& fs : fsVec.second) vvf.insert(fs);

    return vvf;
}

int main(int argc, char *argv[])
{
    int opt, option_index = 0;
    bool doBackground = false, doTopTagger = false, doEventSelection = false, 
         doEventShape = false, do0Lep = false, doStealthTT = false, doSHuHd = false, 
         doDataMC = false, doTest = false, doHadTrigger = false, doTemp = false,
         doCleanJet = false;
    bool runOnCondor = false;
    bool isSkim = false;
    std::string histFile = "", dataSets = "";
    int nFiles = -1, startFile = 0, maxEvts = -1;

    static struct option long_options[] = {
        {"doBackground",       no_argument, 0, 'b'},
        {"doTopTagger",        no_argument, 0, 't'},
        {"doEventSelection",   no_argument, 0, 's'},
        {"doEventShape",       no_argument, 0, 'p'},
        {"do0Lep",             no_argument, 0, 'z'},
        {"doStealthTT",        no_argument, 0, 'x'},
        {"doSHuHd",            no_argument, 0, 'k'},
        {"doDataMC",           no_argument, 0, 'm'},
        {"doTest",             no_argument, 0, 'a'},
        {"doHadTrigger",       no_argument, 0, 'T'},
        {"doTemp",             no_argument, 0, 'e'},
        {"doCleanJet",         no_argument, 0, 'f'},
        {"condor",             no_argument, 0, 'c'},
        {"histFile",     required_argument, 0, 'H'},
        {"dataSets",     required_argument, 0, 'D'},
        {"numFiles",     required_argument, 0, 'N'},
        {"startFile",    required_argument, 0, 'M'},
        {"numEvts",      required_argument, 0, 'E'},
    };

    while((opt = getopt_long(argc, argv, "btspzxkmaTefcH:D:N:M:E:", long_options, &option_index)) != -1)
    {
        switch(opt)
        {
            case 'b': doBackground     = true;              break;
            case 't': doTopTagger      = true;              break;
            case 's': doEventSelection = true;              break;
            case 'p': doEventShape     = true;              break;
            case 'z': do0Lep           = true;              break;
            case 'x': doStealthTT      = true;              break;
            case 'k': doSHuHd          = true;              break;
            case 'm': doDataMC         = true;              break;
            case 'a': doTest           = true;              break;
            case 'T': doHadTrigger     = true;              break;
            case 'e': doTemp           = true;              break;
            case 'f': doCleanJet       = true;              break;
            case 'c': runOnCondor      = true;              break;
            case 'H': histFile         = optarg;            break;
            case 'D': dataSets         = optarg;            break;
            case 'N': nFiles           = int(atoi(optarg)); break;
            case 'M': startFile        = int(atoi(optarg)); break;
            case 'E': maxEvts          = int(atoi(optarg)); break;
        }
    }

    if(dataSets.find("skim") != std::string::npos) isSkim = true;  

    if(runOnCondor)
    {
        char thistFile[128];
        sprintf(thistFile, "MyAnalysis_%s_%d.root", dataSets.c_str(), startFile);
        histFile = thistFile;
    }

    std::set<AnaSamples::FileSummary> vvf = setFS(dataSets, runOnCondor); 
    TFile* outfile = TFile::Open(histFile.c_str(), "RECREATE");
    
    try
    {
        if(doBackground)
        {
            printf("\n\n running AnalyzeBackground\n\n" ) ;
            run<AnalyzeBackground>(vvf,startFile,nFiles,maxEvts,isSkim,outfile);
        }
        else if(doTopTagger)
        {
            printf("\n\n running AnalyzeTopTagger\n\n" ) ;
            run<AnalyzeTopTagger>(vvf,startFile,nFiles,maxEvts,isSkim,outfile);
        }
        else if(doEventSelection)
        {
            printf("\n\n running AnalyzeEventSelection\n\n") ;
            run<AnalyzeEventSelection>(vvf,startFile,nFiles,maxEvts,isSkim,outfile);
        }
        else if(doEventShape)
        {
            printf("\n\n running AnalyzeEventShape\n\n") ;
            run<AnalyzeEventShape>(vvf,startFile,nFiles,maxEvts,isSkim,outfile);
        }
        else if(do0Lep)
        {
            printf("\n\n running Analyze0Lep\n\n") ;
            run<Analyze0Lep>(vvf,startFile,nFiles,maxEvts,isSkim,outfile);
        }
        else if(doStealthTT)
        {
            printf("\n\n running AnalyzeStealthTopTagger\n\n") ;
            run<AnalyzeStealthTopTagger>(vvf,startFile,nFiles,maxEvts,isSkim,outfile);
        }
        else if(doSHuHd)
        {
            printf("\n\n running AnalyzeSHuHdEventSelection\n\n") ;
            run<AnalyzeSHuHdEventSelection>(vvf,startFile,nFiles,maxEvts,isSkim,outfile);
        }
        else if(doDataMC)
        {
            printf("\n\n running AnalyzeDataMC\n\n") ;
            run<AnalyzeDataMC>(vvf,startFile,nFiles,maxEvts,isSkim,outfile);
        }
        else if(doTest)
        {
            printf("\n\n running AnalyzeTest\n\n") ;
            run<AnalyzeTest>(vvf,startFile,nFiles,maxEvts,isSkim,outfile);
        }
        else if(doTemp)
        {
            printf("\n\n running AnalyzeTemp\n\n") ;
            run<AnalyzeTemp>(vvf,startFile,nFiles,maxEvts,isSkim,outfile);
        }
        else if(doCleanJet)
        {
            printf("\n\n running AnalyzeCleanJet\n\n") ;
            run<AnalyzeCleanJet>(vvf,startFile,nFiles,maxEvts,isSkim,outfile);
        }
        else if(doHadTrigger)
        {
            printf("\n\n runnning AnalyzeHadTrigger\n\n") ;
            run<AnalyzeHadTrigger>(vvf,startFile,nFiles,maxEvts,isSkim,outfile);
        }
    
        outfile->Close();
    }
    catch(const std::string e)
    {
        std::cout << e << std::endl;
        return 0;
    }
    catch(const TTException e)
    {
        std::cout << e << std::endl;
        return 0;
    }
    catch(const SATException e)
    {
        std::cout << e << std::endl;
        return 0;
    }

    return 0;
}
