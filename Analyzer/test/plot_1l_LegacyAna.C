#include "Analyzer/Analyzer/test/plotterLegacyAna.h"
#include <iostream>
#include <getopt.h>

// This plotting script is setup to make Fig. 2 of the SUS-19-004 paper
// as well as the supplementary plots in Fig. 9 and 10 found on the twiki page
// https://twiki.cern.ch/twiki/bin/view/CMS/SUS19004SupplementaryMaterial

class FisherHolder
{
public:
    std::vector<std::string> cutNames_;
    std::string plotName_;
};

void setHistInfo(const std::string& path, std::vector<histInfo>& data, std::vector<histInfo>& bg, std::vector<histInfo>& sig, std::vector<histInfo>& syst, std::vector<histInfo>& rsyst, const std::string& year, int color = 0, const bool doQCD = false)
{
    //entry for data
    //this uses the initializer syntax to initialize the histInfo object
    //               leg entry root file                 draw options  draw color
    data = {
        {"Data" , path + "/"+year+"_Data.root", "PEX0", kBlack},
    };
    
    if (!doQCD) {
        bg = {
            {"t#bar{t} + X", path + "/"+year+"_TTX.root",   "hist", kOrange + 2  , 1.0, true, 0},
            {"QCD multijet", path + "/"+year+"_QCD.root",   "hist", kGreen + 1   , 1.0, true, 0},
            {"Other",        path + "/"+year+"_Other.root", "hist", kMagenta + 2 , 1.0, true, 0},        
            {"t#bar{t}",     path + "/"+year+"_TT.root",    "hist", kBlue - 6    , 1.0, true, 0},
        };
    } else {
        bg = {
            {"t#bar{t} + X", path + "/"+year+"_TTX.root",   "hist", kOrange + 2  , 1.0, true, 0},
            {"Other",        path + "/"+year+"_Other.root", "hist", kMagenta + 2 , 1.0, true, 0},        
            {"t#bar{t}",     path + "/"+year+"_TT.root",    "hist", kBlue - 6    , 1.0, true, 0},
            {"QCD multijet", path + "/"+year+"_QCD.root",   "hist", kGreen + 1   , 1.0, true, 0},
        };
    }

    //sig = {        
    //    {"RPV m_{#tilde{t}} = 350 GeV (#sigma_{#tilde{t} #bar{#tilde{t}}} #times 4)", path + "/"+year+"_RPV_2t6j_mStop-350.root",  "hist", kMagenta, printNEvents, 4.0, false, 2},
    //    {"RPV m_{#tilde{t}} = 850 GeV (#sigma_{#tilde{t} #bar{#tilde{t}}} #times 16)", path + "/"+year+"_RPV_2t6j_mStop-850.root",  "hist", kRed   , printNEvents, 16.0, false, 3},
    //    {"Stealth SYY m_{#tilde{t}} = 350 GeV (#sigma_{#tilde{t} #bar{#tilde{t}}} #times 4)", path + "/"+year+"_StealthSYY_2t6j_mStop-350.root",  "hist", kCyan + 1, printNEvents, 4.0, false, 1},        
    //    {"Stealth SYY m_{#tilde{t}} = 850 GeV (#sigma_{#tilde{t} #bar{#tilde{t}}} #times 16)", path + "/"+year+"_StealthSYY_2t6j_mStop-850.root",  "hist", kBlue, printNEvents, 16.0, false, 4},        
    //};

    sig = {        
        {"RPV m_{#tilde{t}} = 350 GeV)",         path + "/"+year+"_RPV_2t6j_mStop-350.root",        "hist", kMagenta,  1.0, false, 2},
        {"RPV m_{#tilde{t}} = 850 GeV)",         path + "/"+year+"_RPV_2t6j_mStop-850.root",        "hist", kRed,      1.0, false, 3},
        {"Stealth SYY m_{#tilde{t}} = 350 GeV)", path + "/"+year+"_StealthSYY_2t6j_mStop-350.root", "hist", kCyan + 1, 1.0, false, 1},        
        {"Stealth SYY m_{#tilde{t}} = 850 GeV)", path + "/"+year+"_StealthSYY_2t6j_mStop-850.root", "hist", kBlue,     1.0, false, 4},        
    };

    if (!doQCD) {
        syst = {
            {"SYST",   path + "/"+year+"_MC_Syst_wNormUncJan08.root", "hist", kBlack, 1.0, true, 0},
        };

        rsyst = {
            {"RSYST",  path + "/"+year+"_MC_Ratio_Syst_wNormUncJan08.root", "hist", kBlack, 1.0, true, 0},
        };
    }
}

int main(int argc, char *argv[])
{
    TH1::AddDirectory(false);

    int opt, option_index = 0;    
    std::string year = "2018post";
    std::string directory = "";
    int approved = 0;

    static struct option long_options[] = {
        {"year",      required_argument, 0, 'y'},
        {"directory", required_argument, 0, 't'},
        {"approved",  optional_argument, 0, 'a'},
    };

    while((opt = getopt_long(argc, argv, "y:t:a:", long_options, &option_index)) != -1)
    {
        switch(opt)
        {
            case 't': directory = optarg; break;
            case 'y': year = optarg; break;
            case 'a': approved = std::stoi(optarg); break;
        }
    }

    std::string path = "condor/hadded/" + directory + "/" + year;

    std::vector<histInfo> data, bg, sig, syst, rsyst;
    setHistInfo(path, data, bg, sig, syst, rsyst, year, 1);
    HistInfoCollection histInfoCollection(data, bg, sig, syst, rsyst);

    //make plotter object with the required sources for histograms specified
    Plotter plt(std::move(histInfoCollection), "PlotsForLegacyAna/" + directory + "/" +year);

    double lumi = 0.0;
    if      (year == "2016")     lumi = 35900.0;
    else if (year == "2017")     lumi = 41500.0;
    else if (year == "2018pre")  lumi = 21100.0;
    else if (year == "2018post") lumi = 38700.0;
    else if (year == "2020")     lumi = 101300.0;

    // --------------------
    // - Make stack plots
    // --------------------

    std::vector<std::string> mycuts_1l = {
        "_1l_HT300_ge7j_ge1b_Mbl",
    };

    for(std::string mycut : mycuts_1l)
    {
        plt.plotStack("h_deepESM"+mycut, "S_{NN}" ,            "Events",        false, 10, true, true, 999.9, -999.9, lumi, approved);
        plt.plotStack("h_njets"+mycut,   "N_{jets}" ,          "Events / bin",  true,  -1, true, true, 999.9, -999.9, lumi, approved);
        plt.plotStack("h_mbl"+mycut,     "M(l,b) [GeV]",       "Events / bin",  true,  10, true, true, 999.9, -999.9, lumi, approved);
        plt.plotStack("h_ht"+mycut,      "H_{T} [GeV]",        "Events / bin",  true,  10, true, true, 999.9, -999.9, lumi, approved);
        plt.plotStack("h_lPt"+mycut,     "Lepton p_{T} [GeV]", "Leptons / bin", true,  2,  true, true, 0,     1000,   lumi, approved);
        plt.plotStack("h_lEta"+mycut,    "Lepton #eta",        "Leptons / bin", false, 10, true, true, 999.9, -999.9, lumi, approved);
        plt.plotStack("h_lPhi"+mycut,    "Lepton #phi",        "Leptons / bin", false, 2,  true, true, 999.9, -999.9, lumi, approved);
        plt.plotStack("h_jPt"+mycut,     "Jet p_{T} [GeV]",    "Jets / bin",    true,  2,  true, true, 0,     1000,   lumi, approved);
    }

    plt.plotStack( "fwm2_top6_1l_ge7j_ge1b",    "FWM2", "Events / bin", false, 1, true, true, 999.0, -999.9, lumi, approved);
    plt.plotStack( "fwm3_top6_1l_ge7j_ge1b",    "FWM3", "Events / bin", false, 1, true, true, 999.0, -999.9, lumi, approved);
    plt.plotStack( "fwm4_top6_1l_ge7j_ge1b",    "FWM4", "Events / bin", false, 1, true, true, 999.0, -999.9, lumi, approved);
    plt.plotStack( "jmt_ev0_top6_1l_ge7j_ge1b", "JMT0", "Events / bin", false, 1, true, true, 999.0, -999.9, lumi, approved);
    plt.plotStack( "jmt_ev1_top6_1l_ge7j_ge1b", "JMT1", "Events / bin", false, 1, true, true, 999.0, -999.9, lumi, approved);
    plt.plotStack( "jmt_ev2_top6_1l_ge7j_ge1b", "JMT2", "Events / bin", false, 1, true, true, 999.0, -999.9, lumi, approved);
    for(unsigned int i = 0; i < 7; i++)
    {
        std::string order = "";
        switch (i)
        {
            case 0: order = "Leading"; break;
            case 1: order = "Subleading"; break;
            case 2: order = "Third"; break;
            case 3: order = "Fourth"; break;
            case 4: order = "Fifth"; break;
            case 5: order = "Sixth"; break;
            case 6: order = "Seventh"; break;
        }
 
        plt.plotStack("Jet_cm_pt_"+std::to_string(i+1)+"_1l_ge7j_ge1b",  order + " Jet" + " p_{T} [GeV]", "Events / 30 GeV", true,   3, true, true, 0, 1500, lumi, approved);
        plt.plotStack("Jet_cm_eta_"+std::to_string(i+1)+"_1l_ge7j_ge1b", order + " Jet" + " #eta",        "Events / bin",    false,  1, true, true, -6, 6, lumi, approved);
        plt.plotStack("Jet_cm_phi_"+std::to_string(i+1)+"_1l_ge7j_ge1b", order + " Jet" + " #phi",        "Events / bin",    false,  1, true, true, -4, 4, lumi, approved);
        plt.plotStack("Jet_cm_m_"+std::to_string(i+1)  +"_1l_ge7j_ge1b", order + " Jet" + " Mass [GeV]",  "Events / 5 GeV",  false,  5, true, true,  0, 200, lumi, approved);
    }
}
