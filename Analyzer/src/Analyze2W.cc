#define Analyze2W_cxx
#include "Analyzer/Analyzer/include/Analyze2W.h"
#include "Framework/Framework/include/Utility.h"
#include "SusyAnaTools/Tools/NTupleReader.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TRandom3.h>
#include <iostream>
#include <TFile.h>
#include <algorithm>
#include <vector>
#include <utility>
#include <numeric>



template<typename... Args>
void getHistogramTarget(Args&&... args){
  return (args() + ... );
}


/*
  Fill("Variable", value, tr);
  
*/

Analyze2W::Analyze2W()
{
    InitHistos();

}

void createNewHistogram(auto& myhistos, std::string name, int v1, int v2, int v3, std::vector<std::string> props){
  name = std::accumulate(props.begin(),props.end(),name, [](std::string& x , auto y){ x += "_"; return x += y;});
    const char* x = name.c_str();
    myhistos.emplace(name,std::make_shared<TH1D>(x,x,v1,v2,v3));
}

void createNewHistogram(auto& myhistos, std::string name, int v1, int v2, int v3){
    createNewHistogram(myhistos, name, v1,v2,v3,{});
}

//Define all your histograms here.
void Analyze2W::InitHistos()
{
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    //This event counter histogram is necessary so that we know that all the condor jobs ran successfully. If not, when you use the hadder script, you will see a discrepancy in red as the files are being hadded.
    createNewHistogram(my_histos, "EventCounter", 2, -1.1, 1.1  ) ;

    const int lepton_min = 0 ;
    const int lepton_max = 2 ;
    for(int i = lepton_min; i <= lepton_max; ++i){
        std::string lep = "lep" + std::to_string(i);
        createNewHistogram(my_histos, "h_njets", 8, 7, 15,  {lep} ) ;
        createNewHistogram(my_histos, "h_met", 200,0, 2000, {lep} ) ;
    }
    //Define 1D histograms
}

//Put everything you want to do per event here.
void Analyze2W::Loop(NTupleReader& tr, double, int maxevents, bool)
{
    std::cout << "Analyzing event";
    const int lepton_min = 0 ;
    const int lepton_max = 2 ;
    while( tr.getNextEvent() )
    {
        //This is added to count the number of events- do not change the next two lines.
        const auto& eventCounter        = tr.getVar<int>("eventCounter");
        my_histos["EventCounter"]->Fill( eventCounter );

        //--------------------------------------------------
        //-- Print Event Number 
        //--------------------------------------------------

        if( maxevents != -1 && tr.getEvtNum() >= maxevents ) break;
        if( tr.getEvtNum() & (10000 == 0) ) printf( " Event %i\n", tr.getEvtNum() );


        #define makeVar(name, type) const type& name = tr.getVar<type>(#name)
        #define makeVec(name, type) const std::vector<type>& name = tr.getVec<type>(#name)

        makeVar(runtype         ,std::string);     
        makeVar(JetID           ,bool);
        makeVec(Jets            ,TLorentzVector);
        makeVar(NJet                      ,int);
        makeVar(NGoodLeptons    ,int);
        makeVar(MET             ,double);
        makeVar(METPhi          ,double);
        makeVec(gen_particles_pdg ,int);
        makeVar(gen_particles ,TLorentzVector);
        makeVar(gen_particles_parent_idx ,int);
        makeVar(gen_particles_parent_id ,int);
        makeVar(gen_particles_status ,int);

        constexpr int W_PDGID    = 24;
        constexpr int T_PDGID    = 6;
        constexpr int STOP_PDGID = 1000006;

        makeVar(passMadHT           ,bool);
        // ------------------------
        // -- Define weight
        // ------------------------
        double weight               = 1.0;
        double eventweight          = 1.0;
        double leptonScaleFactor    = 1.0;
        double bTagScaleFactor      = 1.0;
        double htDerivedScaleFactor = 1.0;
        double prefiringScaleFactor = 1.0;
        double puScaleFactor        = 1.0;

        if(runtype == "MC")
        {
            if( !passMadHT ) continue; //Make sure not to double count DY events
            // Define Lumi weight
            const auto& Weight  = tr.getVar<double>("Weight");
            const auto& lumi = tr.getVar<double>("Lumi");
            eventweight = lumi*Weight;

            // Define lepton weight

            //PileupWeight = tr.getVar<double>("_PUweightFactor");
            // bTagScaleFactor   = tr.getVar<double>("bTagSF_EventWeightSimple_Central");
            htDerivedScaleFactor = tr.getVar<double>("htDerivedweight");
            prefiringScaleFactor = tr.getVar<double>("prefiringScaleFactor");
            puScaleFactor = tr.getVar<double>("puWeightCorr");

            weight *= eventweight*leptonScaleFactor*bTagScaleFactor*htDerivedScaleFactor*prefiringScaleFactor*puScaleFactor;
        }

        if(NGoodLeptons <= 2){
            std::string name = "h_njets_lep" + std::to_string(NGoodLeptons);
            my_histos[name]->Fill(NJet, weight);
        }else {
          std::string name = "h_njets_lep-overflow";
            my_histos[name]->Fill(NJet, weight);
        }
    }
}

void Analyze2W::WriteHistos(TFile* outfile)
{
    outfile->cd();

    for (const auto &p : my_histos) {
        p.second->Write();
    }

    for (const auto &p : my_2d_histos) {
        p.second->Write();
    }

    for (const auto &p : my_efficiencies) {
        p.second->Write();
    }

}
