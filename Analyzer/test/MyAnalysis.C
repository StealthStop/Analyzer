#include "NTupleReader/include/NTupleReader.h"

#include "Framework/Framework/include/samples.h"
#include "Framework/Framework/include/MiniTupleMaker.h"
#include "Framework/Framework/include/Utility.h"

#include "TopTagger/CfgParser/interface/TTException.h"

#include "Analyzer/Analyzer/include/AnalyzeDoubleDisCo.h"
#include "Analyzer/Analyzer/include/AnalyzeTest.h"
#include "Analyzer/Analyzer/include/AnalyzeLepTrigger.h"
#include "Analyzer/Analyzer/include/MakeNJetDists.h"
#include "Analyzer/Analyzer/include/MakeMiniTree.h"
#include "Analyzer/Analyzer/include/MakeQCDValTree.h"
#include "Analyzer/Analyzer/include/MakeTopTagSFTree.h"
#include "Analyzer/Analyzer/include/MakeAnaSkimTree.h"
#include "Analyzer/Analyzer/include/CalculateBTagSF.h"
#include "Analyzer/Analyzer/include/CalculateTopTagSF.h"
#include "Analyzer/Analyzer/include/CalculateSFMean.h"
#include "Analyzer/Analyzer/include/Config.h"
#include "Analyzer/Analyzer/include/Semra_Analyzer.h"
#include "Analyzer/Analyzer/include/ResolvedTopTagger_Analyzer.h"
#include "Analyzer/Analyzer/include/HEM_Analyzer.h"
#include "Analyzer/Analyzer/include/TopTaggerSF_Analyzer.h"
#include "Analyzer/Analyzer/include/ISRJets_Analyzer.h"
#include "Analyzer/Analyzer/include/HadTriggers_Analyzer.h"
#include "Analyzer/Analyzer/include/StealthHemispheres.h"
#include "Analyzer/Analyzer/include/AnalyzeTemplate.h"
#include "Analyzer/Analyzer/include/MakeNNVariables.h"
#include "Analyzer/Analyzer/include/AnalyzeGenStop.h"
#include "Analyzer/Analyzer/include/AnalyzeXsec.h"

#include "TH1D.h"
#include "TFile.h"
#include "TChain.h"
#include "TGraph.h"

#include <iostream>
#include <getopt.h>
#include <string>
#include <functional>
#include <unistd.h>

const std::string getFullPath(const std::string& file)
{
    char buf[512];
    int count = readlink(file.c_str(), buf, sizeof(buf));
    if(count >= 0)
    {
        buf[count] = '\0';
        return std::string(buf);
    }
    else
    {
        std::cout<<"Could not get full path of "<<file<<std::endl;
        return std::string();
    }
}

template<typename Analyze> void run(const std::set<AnaSamples::FileSummary>& vvf, const std::string dataSets,
                                    const int startFile, const int nFiles, const int maxEvts, 
                                    TFile* const outfile, const bool isQuiet, const bool fastMode, const std::string& analyzer)
{
    std::cout << "Initializing..." << std::endl;
    Analyze a;
    for(const auto& file : vvf)
    {
        // Define what is needed per sample set
        std::cout << "Running over sample " << file.tag << std::endl;
        TChain* ch = new TChain( (file.treePath).c_str() );
        file.addFilesToChain(ch, startFile, nFiles);
        NTupleReader tr(ch, {"RunNum"});
        const std::string runtype = (file.tag.find("Data") != std::string::npos) ? "Data" : "MC";
        tr.registerDerivedVar("runtype",runtype);
        tr.registerDerivedVar("filetag",file.tag);
        tr.registerDerivedVar("analyzer",analyzer);
        tr.registerDerivedVar("dataset",dataSets);

        // Registerd event weight computed by FileSummary class, no sign information
        tr.registerDerivedVar("weightAbsVal", file.getWeight());

        // Determine if smart early-abort-event-module-pipeline is active
        tr.registerDerivedVar("fastMode", fastMode);

        printf( "runtype: %s nFiles: %i startFile: %i maxEvts: %i \n",runtype.c_str(),nFiles,startFile,maxEvts ); fflush( stdout );

        // Define classes/functions that add variables on the fly        
        Config c;
        c.setUp(tr);

        // Loop over all of the events and fill histos
        std::cout << "Starting event loop (in run)" << std::endl;
        a.Loop(tr, 1.0, maxEvts, isQuiet);

        // Cleaning up dynamic memory
        delete ch;            
    }
    std::cout << "Writing histograms..." << std::endl;
    a.WriteHistos(outfile);
}

std::set<AnaSamples::FileSummary> setFS(const std::string& dataSets)
{
    AnaSamples::SampleSet        ss("sampleSets.cfg");
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
    if(vvf.size() == 0)             std::cout<< utility::color("No samples for \""+std::string(dataSets)+"\" in the sampleSet.cfg or sampleCollection.cfg","red")  <<std::endl;
    else if(vvf.begin()->tag == "") std::cout<< utility::color("A filetag is empty, Check if all sampleSet(s) make sense for \""+std::string(dataSets)+"\"","red") <<std::endl;
    return vvf;
}

int main(int argc, char *argv[])
{
    int opt, option_index = 0;
    bool runOnCondor = false, isQuiet = true, fastMode = false;
    std::string histFile = "", dataSets = "", analyzer = "";
    int nFiles = -1, startFile = 0, maxEvts = -1;

    static struct option long_options[] = {
        {"condor",             no_argument, 0, 'c'},
        {"verbose",            no_argument, 0, 'v'},
        {"fastMode",           no_argument, 0, 's'},
        {"analyzer",     required_argument, 0, 'A'},
        {"histFile",     required_argument, 0, 'H'},
        {"dataSets",     required_argument, 0, 'D'},
        {"numFiles",     required_argument, 0, 'N'},
        {"startFile",    required_argument, 0, 'M'},
        {"numEvts",      required_argument, 0, 'E'},
    };

    // here is the options to run the codes / can add options
    while((opt = getopt_long(argc, argv, "cvsA:H:D:N:M:E:", long_options, &option_index)) != -1)
    {
        switch(opt)
        {
            case 'c': runOnCondor       = true;              break;
            case 'v': isQuiet           = false;             break;
            case 's': fastMode          = true;              break;
            case 'A': analyzer          = optarg;            break;
            case 'H': histFile          = optarg;            break;
            case 'D': dataSets          = optarg;            break;
            case 'N': nFiles            = int(atoi(optarg)); break;
            case 'M': startFile         = int(atoi(optarg)); break;
            case 'E': maxEvts           = int(atoi(optarg)); break;
        }
    }

    if(runOnCondor)
    {
        char thistFile[128];
        sprintf(thistFile, "MyAnalysis_%s_%d.root", dataSets.c_str(), startFile);
        histFile = thistFile;
    }

    std::set<AnaSamples::FileSummary> vvf = setFS(dataSets); 
    TFile* outfile = TFile::Open(histFile.c_str(), "RECREATE");

    std::vector<std::pair<std::string, std::function<void(const std::set<AnaSamples::FileSummary>&, const std::string, const int,const int,const int,TFile* const,const bool,const bool,const std::string&)>>> AnalyzerPairVec = {
        {"AnalyzeDoubleDisCo",         run<AnalyzeDoubleDisCo>        },
        {"AnalyzeLepTrigger",          run<AnalyzeLepTrigger>         },
        {"AnalyzeTest",                run<AnalyzeTest>               },
        {"CalculateBTagSF",            run<CalculateBTagSF>           },
        {"CalculateTopTagSF",          run<CalculateTopTagSF>         },
        {"CalculateSFMean",            run<CalculateSFMean>           },
        {"MakeMiniTree",               run<MakeMiniTree>              },
        {"MakeQCDValTree",             run<MakeQCDValTree>            },
        {"MakeTopTagSFTree",           run<MakeTopTagSFTree>          },
        {"MakeAnaSkimTree",            run<MakeAnaSkimTree>           },
        {"MakeNJetDists",              run<MakeNJetDists>             },
        {"Semra_Analyzer",             run<Semra_Analyzer>            },
        {"ResolvedTopTagger_Analyzer", run<ResolvedTopTagger_Analyzer>},
        {"HEM_Analyzer",               run<HEM_Analyzer>              },
        {"TopTaggerSF_Analyzer",       run<TopTaggerSF_Analyzer>      },
        {"ISRJets_Analyzer",           run<ISRJets_Analyzer>          },
        {"HadTriggers_Analyzer",       run<HadTriggers_Analyzer>      },
        {"StealthHemispheres",         run<StealthHemispheres>        },
        {"AnalyzeTemplate",            run<AnalyzeTemplate>           },
        {"MakeNNVariables",            run<MakeNNVariables>           },
        {"AnalyzeGenStop",             run<AnalyzeGenStop>            },
        {"AnalyzeXsec",                run<AnalyzeXsec>               },

    }; 

    try
    {
        bool foundAnalyzer = false;
        for(auto& pair : AnalyzerPairVec)
        {
            if(pair.first==analyzer) 
            {
                std::cout<<"Running the " << analyzer << " Analyzer" <<std::endl;
                pair.second(vvf,dataSets,startFile,nFiles,maxEvts,outfile,isQuiet,fastMode,analyzer); 
                foundAnalyzer = true;
            }
        }

        if (!foundAnalyzer)
        {
            std::cout << utility::color("ERROR: The analyzer \"" + analyzer + "\" is not an analyzer option! Please add it to the MyAnalysis.C list.", "red") << std::endl;        
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
    catch(const NTRException e)
    {
        std::cout << e << std::endl;
        return 0;
    }

    return 0;
}
