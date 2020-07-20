#include "TH1.h"
#include "THStack.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"

#include <memory>
#include <vector>
#include <string>
#include <cstdio>
#include <iostream>
#include <algorithm>

#include "Framework/Framework/include/Utility.h"

// ----------------------------------------------------
// -- keep the plot from overlapping with the legend 
// ----------------------------------------------------
void smartMax(const TH1 * const h, const TLegend* const l, const TPad* const p, double& gmin, double& gmax, double& gpThreshMax, const bool error)
{
    //const bool isLog  = p->GetLogy();
    double min    = 9e99;
    double max    = -9e99;
    double pThreshMax = -9e99;
    int threshold     = static_cast<int>(h->GetNbinsX()*(l->GetX1() - p->GetLeftMargin())/((1 - p->GetRightMargin()) - p->GetLeftMargin()));

    for(int i = 1; i <= h->GetNbinsX(); ++i) {

        double bin = 0.0;
        
        if(error) 
            bin = h->GetBinContent(i) + h->GetBinError(i);
        else      
            bin = h->GetBinContent(i);
        if(bin > max) 
            max = bin;
        else if(bin > 1e-10 && bin < min) 
            min = bin;
        if(i >= threshold && bin > pThreshMax) 
            pThreshMax = bin;
    }

    gpThreshMax = std::max(gpThreshMax, pThreshMax);
    gmax    = std::max(gmax, max);
    gmin    = std::min(gmin, min);

}

// -----------------------------------------------------
// -- CLASS1: hold TH1* with various helper functions 
// -----------------------------------------------------
class histInfo
{
public:
    std::string legName, legEntry, histFile, histName, drawOptions;
    int color, rebin;
    double nEvents;
    std::shared_ptr<TH1> h;
    bool drawHisto; // define for empty data plot

    // -----------------------------------------------------------
    // get histogram from file and configure its optional settings
    // -----------------------------------------------------------
    void retrieveHistogram()
    {
        if(drawHisto) // define for empty data plot
        {
            //Open the file for this histogram
            TFile *f = TFile::Open(histFile.c_str());

            if(!f) {
                printf("File \"%s\" could not be opened!!!\n", histFile.c_str());
                h = nullptr;
                return;
            }

            //get histogram & close file
            h.reset(static_cast<TH1*>(f->Get(histName.c_str())));
            f->Close();
            delete f;

            if(!h) {
                printf("Histogram \"%s\" could not be found in file \"%s\"!!!\n", histName.c_str(), histFile.c_str());
                return;
            }

            h->SetLineColor(color);
            h->SetLineWidth(3);
            h->SetMarkerColor(color);
            h->SetMarkerStyle(20);

            //Get number of events and save it as a member variable
            nEvents  = h->Integral();
            legEntry = legName + " ("+std::to_string(nEvents)+") ";     

            // rebin the histogram if desired
            if(rebin >0) 
                h->Rebin(rebin);
        } 
    }

    // -------------
    // help for axes
    // ------------- 
    void setupAxes()
    {
        h->SetStats(0);
        h->SetTitle(0);
        h->GetYaxis()->SetTitleOffset(1.2);
        h->GetXaxis()->SetTitleOffset(1.1);
        h->GetXaxis()->SetTitleSize(0.045);
        h->GetXaxis()->SetLabelSize(0.045);
        h->GetYaxis()->SetTitleSize(0.045);
        h->GetYaxis()->SetLabelSize(0.045);
        if(h->GetXaxis()->GetNdivisions() % 100 > 5) 
            h->GetXaxis()->SetNdivisions(6, 5, 0);
    }

    // -------------------------
    // wrapper to draw histogram
    // -------------------------
    void draw(const std::string& additionalOptions = "", bool noSame = false) const 
    {
        h->Draw(((noSame?"":"same " + drawOptions + " " + additionalOptions)).c_str());
    }

    // ------------------------------
    // add for signal plot line style
    // ------------------------------
    void setLineStyle(int style = 1)
    {
        h->SetLineStyle(style);
    }

    // ------------
    // setFillColor 
    // ------------
    void setFillColor(int newColor = -1) {
        if(newColor >= 0) 
            h->SetFillColor(newColor);
        else          
            h->SetFillColor(color);
    }

    // constructor & distructor
    histInfo(const std::string& legName, const std::string& histFile, const std::string& drawOptions, const int color, const bool drawHisto = true) : legName(legName), legEntry(legName), histFile(histFile), histName(""), drawOptions(drawOptions), color(color), rebin(-1), nEvents(-1), h(nullptr), drawHisto(drawHisto)
    {
    }

    histInfo(TH1* h) : legName(h->GetName()), legEntry(h->GetName()), histFile(""), histName(h->GetName()), drawOptions(""), color(0), rebin(0), nEvents(-1), h(h)
    {
    }

    ~histInfo()
    {
    }
};

// ---------------------
// -- CLASS2: Plotter 
// ---------------------
class Plotter
{

private:
    //entry for data
    histInfo data_;
    //vector summarizing background & signal histograms to include in the plot
    std::vector<histInfo> bgEntries_;
    std::vector<histInfo> sigEntries_;

    static bool compareNEvents(histInfo h1, histInfo h2) 
    {
        const double nE1 = h1.nEvents;
        const double nE2 = h2.nEvents;  
        return (nE1 > nE2);
    }
    
public:
    Plotter(histInfo&& data, std::vector<histInfo>&& bgEntries, std::vector<histInfo>&& sigEntries) : data_(data), bgEntries_(bgEntries), sigEntries_(sigEntries) {}

    // ---------------------
    // function for plotting
    // ---------------------   
    void plot(const std::string& histName, const std::string& xAxisLabel, const std::string& yAxisLabel = "Events", const bool isLogY = false, const std::string& cutlabel = "", const double xmin = 999.9, const double xmax = -999.9, int rebin = -1, double lumi = 35900) // lumi 2016 = 35900, lumi 2017 = 41500, lumi 2018 = 59740
    {
        TH1::AddDirectory(false);

        // ----------------------------------------
        // create the canvas & TLegend for the plot
        // ----------------------------------------
        TCanvas *c = new TCanvas("c1", "c1", 800, 800);
        c->cd();
        TLegend *leg = new TLegend(0.50, 0.56, 0.89, 0.88); 
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetLineWidth(1);
        leg->SetNColumns(1);
        leg->SetTextFont(42); 

        // get maximum from histos and fill TLegend
        double min  = 0.0;
        double max  = 0.0;
        double lmax = 0.0;

        // -------------------------------------
        // create the THStack for the background
        // ------------------------------------- 
        THStack *bgStack = new THStack();
        //Make seperate histogram from sum of BG histograms  
        TH1* hbgSum      = nullptr;

        // -----------------------------------------------------------------------
        // background / get new histogram & sort the bgEntries by number of events
        // -----------------------------------------------------------------------
        for(auto& h : bgEntries_)
        {
            h.histName = histName;
            h.rebin    = rebin;
            h.retrieveHistogram();
        }
       
        // sort the bgEntries by number of events
        std::sort(bgEntries_.begin(), bgEntries_.end(), compareNEvents);

        for(int iBG = bgEntries_.size() - 1; iBG >= 0; --iBG) 
        {
            bgStack->Add(bgEntries_[iBG].h.get(), bgEntries_[iBG].drawOptions.c_str());

            if(!hbgSum)
                hbgSum = static_cast<TH1*>(bgEntries_[iBG].h->Clone());
            else
            hbgSum->Add(bgEntries_[iBG].h.get());
        }

        // ------------------------
        // data / get new histogram 
        // ------------------------
        data_.histName = histName;
        data_.rebin    = rebin;
        data_.retrieveHistogram();
        
        if(data_.drawHisto) // define for empty data plot  
            leg->AddEntry(data_.h.get(), data_.legEntry.c_str(), data_.drawOptions.c_str());
        smartMax(hbgSum, leg, static_cast<TPad*>(gPad), min, max, lmax, true);
    
        // ----------
        // background 
        // ---------- 
        for(auto& entry : bgEntries_) 
        {
            entry.setFillColor();

            //add histograms to TLegend
            leg->AddEntry(entry.h.get(), entry.legEntry.c_str(), "F");
        }
        smartMax(hbgSum, leg, static_cast<TPad*>(gPad), min, max, lmax, false);

        // --------------------------
        // signal / get new histogram
        // --------------------------
        for(auto& entry : sigEntries_) 
        {
            entry.histName = histName;
            entry.rebin    = rebin;
            entry.retrieveHistogram();
            entry.setLineStyle(2);

            //add histograms to TLegend
            leg->AddEntry(entry.h.get(), entry.legEntry.c_str(), "L"); 
            smartMax(entry.h.get(), leg, static_cast<TPad*>(gPad), min, max, lmax, false);
        }
    
        // -----------------
        // set Canvas margin
        // ----------------- 
        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(0.06);
        gPad->SetTopMargin(0.08);
        gPad->SetBottomMargin(0.12);

        // -------------------------------------------
        // create a dummy histogram to act as the axes
        // -------------------------------------------
        histInfo dummy(new TH1D("dummy", "dummy", 1000, hbgSum->GetBinLowEdge(1), hbgSum->GetBinLowEdge(hbgSum->GetNbinsX()) + hbgSum->GetBinWidth(hbgSum->GetNbinsX())));
        dummy.setupAxes();
        dummy.h->GetYaxis()->SetTitle(yAxisLabel.c_str());
        dummy.h->GetXaxis()->SetTitle(xAxisLabel.c_str());
        
        // set x-axis range
        if(xmin < xmax)
            dummy.h->GetXaxis()->SetRangeUser(xmin, xmax);

        // set y-axis range 
        if(isLogY) {

            double locMin  = std::min(0.2, std::max(0.2, 0.05 * min));
            double legSpan = (log10(3*max) - log10(locMin)) * (leg->GetY1() - gPad->GetBottomMargin()) / ((1 - gPad->GetTopMargin()) - gPad->GetBottomMargin());
            double legMin  = legSpan + log10(locMin);

            if(log10(lmax) > legMin) {
                double scale = (log10(lmax) - log10(locMin)) / (legMin - log10(locMin));
                max = pow(max/locMin, scale)*locMin;
        }
        dummy.h->GetYaxis()->SetRangeUser(locMin, 10*max);
    
        } else {
        
            double locMin = 0.0;
            double legMin = (1.2*max - locMin) * (leg->GetY1() - gPad->GetBottomMargin()) / ((1 - gPad->GetTopMargin()) - gPad->GetBottomMargin());
            if(lmax > legMin) max *= (lmax - locMin)/(legMin - locMin);
                dummy.h->GetYaxis()->SetRangeUser(0.0, max*1.2);
        }


        dummy.draw(); 

        gPad->SetLogy(isLogY);

        bgStack->Draw("same");

        // ----
        // draw
        // ----
        // plot signal histograms
        for(const auto& entry : sigEntries_) 
        {
            entry.draw();
        }
    
        // plot data histogram
        if(data_.drawHisto) // define for empty data plot
            data_.draw();

        // plot legend
        leg->Draw("same");

        // draw dummy hist again to get axes on top of histograms
        dummy.draw("AXIS");

        // --------------------------------------------------
        // draw CMS & lumi lables & cut labels & significance
        // --------------------------------------------------
        char lumistamp[128];
        sprintf(lumistamp, "%.1f fb^{-1} (13 TeV)", lumi / 1000.0);
    
        TLatex mark;
        mark.SetNDC(true);
        
        mark.SetTextAlign(11);
        mark.SetTextSize(0.050);
        mark.SetTextFont(61);
        mark.DrawLatex(gPad->GetLeftMargin(), 1 - (gPad->GetTopMargin() - 0.017), "CMS"); 
        
        mark.SetTextSize(0.030); 
        mark.SetTextFont(52);
        mark.DrawLatex(gPad->GetLeftMargin() + 0.11, 1 - (gPad->GetTopMargin() - 0.017), "Preliminary"); 
        
        mark.SetTextFont(42); 
        mark.SetTextAlign(31);
        mark.DrawLatex(1 - gPad->GetRightMargin(), 1 - (gPad->GetTopMargin() - 0.017), lumistamp);

        mark.SetTextAlign(11);
        mark.SetTextFont(42);
        mark.SetTextSize(0.030);
        mark.DrawLatex(0.51, 0.89, cutlabel.c_str()); 

        // calculate significance bin by bin
        double sig = 0.0;
        for(int i = 0; i < sigEntries_.at(0).h->GetNbinsX(); i++)
        {
            const double totBG = hbgSum->GetBinContent(i);
            const double nSig = sigEntries_.at(0).h->GetBinContent(i);
            if(totBG > 1.0 && nSig > 1.0)
            {
                const double s = nSig / sqrt( totBG + pow ( 0.3*totBG, 2) ) ;
                sig = utility::addInQuad(sig, s);
            }
        }
        TLatex significance;  
        significance.SetNDC(true);
        significance.SetTextAlign(11);
        significance.SetTextFont(52);
        significance.SetTextSize(0.030);
        significance.DrawLatex(0.15, 0.89, ("Significance = "+std::to_string(sig)).c_str());

        // -----------------
        // save plots as pdf
        // -----------------
        //c->Print((histName + ".png").c_str());
        c->Print((histName + ".pdf").c_str());

        // -----------------------
        // clean up dynamic memory
        // -----------------------
        delete c;
        delete leg;
        delete bgStack;
        delete hbgSum;
    }
};

// -------------------
// -- Main function 
// -------------------
int main()
{
    // --------------
    // entry for data
    // --------------
    // 'leg entry'  'root file'     'draw options'  'draw color'
    histInfo data = {"Data", "condor/hadd_2016_wholeTT_WP98_GoodJetsMaskTest_16.07.2020/2016_BG_OTHER.root", "PEX0", kBlack, false};

    // ----------------------------------------
    // vector summarizing background histograms 
    // ----------------------------------------
    std::vector<histInfo> bgEntries = {

        {"t#bar{t}",        "condor/hadd_2016_wholeTT_WP98_GoodJetsMaskTest_16.07.2020/2016_TT.root",              "hist", kBlue - 6    },
        {"t#bar{t}+X",      "condor/hadd_2016_wholeTT_WP98_GoodJetsMaskTest_16.07.2020/2016_TTX.root",             "hist", kOrange + 2   },
        {"QCD",             "condor/hadd_2016_wholeTT_WP98_GoodJetsMaskTest_16.07.2020/2016_QCD.root",             "hist", kGreen + 1  },
        {"Other",           "condor/hadd_2016_wholeTT_WP98_GoodJetsMaskTest_16.07.2020/2016_Other.root",           "hist", kMagenta + 2  },
         
        //{"WJets",           "condor/hadd_2016_wholeTT_WP98_GoodJetsMaskTest_16.07.2020/2016_WJets.root",           "hist", kYellow + 1 },
        //{"DYJetsToLL_M-50", "condor/hadd_2016_wholeTT_WP98_GoodJetsMaskTest_16.07.2020/2016_DYJetsToLL_M-50.root", "hist", kOrange + 2 },
        //{"ST",              "condor/hadd_2016_wholeTT_WP98_GoodJetsMaskTest_16.07.2020/2016_ST.root",              "hist", kRed + 1    },
        //{"Diboson",         "condor/hadd_2016_wholeTT_WP98_GoodJetsMaskTest_16.07.2020/2016_Diboson.root",         "hist", kMagenta + 1},
        //{"Triboson",        "condor/hadd_2016_wholeTT_WP98_GoodJetsMaskTest_16.07.2020/2016_Triboson.root",        "hist", kGray       },

    
    };
    
    // ------------------------------------
    // vector summarizing signal histograms
    // ------------------------------------ 
    std::vector<histInfo> sigEntries = { 

        {"RPV m_{#tildet} = 350", "condor/hadd_2016_wholeTT_WP98_GoodJetsMaskTest_16.07.2020/2016_RPV_2t6j_mStop-350.root",        "hist", kCyan+1    }, 
        {"RPV m_{#tildet} = 550", "condor/hadd_2016_wholeTT_WP98_GoodJetsMaskTest_16.07.2020/2016_RPV_2t6j_mStop-550.root",        "hist", kMagenta+2 },
        {"RPV m_{#tildet} = 850", "condor/hadd_2016_wholeTT_WP98_GoodJetsMaskTest_16.07.2020/2016_RPV_2t6j_mStop-850.root",        "hist", kRed },
        //{"SYY m_{#tildet} = 900", "condor/hadd_2016_wholeTT_WP98_GoodJetsMaskTest_16.07.2020/2016_StealthSYY_2t6j_mStop-900.root", "hist", kGray },
    
    };

    // -------------------
    // make plotter object
    // ------------------- 
    Plotter plt(std::move(data), std::move(bgEntries), std::move(sigEntries));

    std::vector<std::string> cut {

        //"", "0l", 
        
        //"0l_HT500", "0l_HT500_ge2b", "0l_HT500_ge2b_ge6j", 
 
        //"0l_HT500_ge2b_ge6j_ge2t", "0l_HT500_ge2b_ge6j_ge2t1j", "0l_HT500_ge2b_ge6j_ge2t3j", "0l_HT500_ge2b_ge6j_ge2t1j3j",
        
        "0l_HT500_ge2b_ge6j_ge2t_ge1dRbjets", 
        //"0l_HT500_ge2b_ge6j_ge2t1j_ge1dRbjets", "0l_HT500_ge2b_ge6j_ge2t3j_ge1dRbjets", "0l_HT500_ge2b_ge6j_ge2t1j3j_ge1dRbjets",
        //"baseline_0l_Njet7", "baseline_0l_Njet8", "baseline_0l_Njet9", 
        //"baseline_0l_Njet10", "baseline_0l_Njet11", "baseline_0l_Njet12",
        //"baseline_0l_Njet13", "baseline_0l_Njet14", "baseline_0l_Njet15",
        
        //"baseline_0l_stopMassesg200_PtRank", "baseline_0l_stopMassesg200_ScalarPtRank",
        //"baseline_0l_stopMassesle200_PtRank", "baseline_0l_stopMassesle200_ScalarPtRank", 

    };

    for (const auto& cutlabel : cut) 
    {
        plt.plot( "h_ntops_"+cutlabel,        "N_{T}",              "Events", true, cutlabel );
        plt.plot( "h_njets_"+cutlabel,        "N_{J}",              "Events", true, cutlabel );
        plt.plot( "h_nbjets_"+cutlabel,       "N_{BJ}",             "Events", true, cutlabel );
        plt.plot( "h_ht_"+cutlabel,           "HT [GeV]",           "Events", true, cutlabel );
        plt.plot( "h_met_"+cutlabel,          "MET [GeV]",          "Events", true, cutlabel );

        plt.plot( "h_topsMass_"+cutlabel,      "M_{Tops} [GeV]",       "Events", true, cutlabel, 0, 50);
        plt.plot( "h_topsEta_"+cutlabel,       "#eta_{Tops}",          "Events", true, cutlabel, 5 );
        plt.plot( "h_topsPhi_"+cutlabel,       "#phi_{Tops}",          "Events", true, cutlabel, 5 );
        plt.plot( "h_topsPt_"+cutlabel,        "pT_{Tops} [GeV]",      "Events", true, cutlabel, 0, 500);

        plt.plot( "h_resolvedMass_"+cutlabel,  "resolved jet M",       "Events", true, cutlabel, 0, 50);
        plt.plot( "h_resolvedEta_"+cutlabel,   "resolved jet #eta",    "Events", true, cutlabel, 5 );
        plt.plot( "h_resolvedPhi_"+cutlabel,   "resolved jet #phi",    "Events", true, cutlabel, 5 );
        plt.plot( "h_resolvedPt_"+cutlabel,    "resolved jet pT",      "Events", true, cutlabel, 0, 500);

        plt.plot( "h_bestTopMass_"+cutlabel,   "M_{BestTop} [GeV]",    "Events", true, cutlabel );
        plt.plot( "h_bestTopEta_"+cutlabel,    "#eta_{BestTop}",       "Events", true, cutlabel );         
        plt.plot( "h_bestTopPt_"+cutlabel,     "pT_{BestTop} [GeV]",   "Events", true, cutlabel );

        plt.plot( "h_dR_bjets_"+cutlabel,      "#DeltaR_{bjets}",      "Events", false, cutlabel ); // true: for not log scale
        plt.plot( "h_dR_top1_top2_"+cutlabel,  "#DeltaR_{t1-t2}",      "Events", false, cutlabel );
        plt.plot( "h_dR_tops_bjets_"+cutlabel, "#DeltaR_{tops-bjets}", "Events", false, cutlabel ); 


        // stop MT2 hemispheres 
        //plt.plot( "h_stop1Mass_PtRank_"+cutlabel,           "Pt Rank M_{#tildet}_{1}",               "Events", true, cutlabel, 0, 1500, 5 );
        //plt.plot( "h_stop1Eta_PtRank_"+cutlabel,            "Pt Rank #eta_{#tildet}_{1}",            "Events", true, cutlabel, 5 );
        //plt.plot( "h_stop1Phi_PtRank_"+cutlabel,            "Pt Rank #phi_{#tildet}_{1}",            "Events", true, cutlabel, 5 );
        //plt.plot( "h_stop1Pt_PtRank_"+cutlabel,             "Pt Rank pT_{#tildet}_{1}",              "Events", true, cutlabel, 0, 1500, 5 );        
        //plt.plot( "h_stop2Mass_PtRank_"+cutlabel,           "Pt Rank M_{#tildet}_{2}",               "Events", true, cutlabel, 0, 1500, 5 );
        //plt.plot( "h_stop2Eta_PtRank_"+cutlabel,            "Pt Rank #eta_{#tildet}_{2}",            "Events", true, cutlabel, 5 );
        //plt.plot( "h_stop2Phi_PtRank_"+cutlabel,            "Pt Rank #phi_{#tildet}_{2}",            "Events", true, cutlabel, 5 );
        //plt.plot( "h_stop2Pt_PtRank_"+cutlabel,             "Pt Rank pT_{#tildet}_{2}",              "Events", true, cutlabel, 0, 1500, 5 );

        //plt.plot( "h_stop1Mass_MassRank_"+cutlabel,         "Mass Rank M_{#tildet}_{1}",             "Events", true, cutlabel, 0, 1500, 5 );
        //plt.plot( "h_stop1Eta_MassRank_"+cutlabel,          "Mass Rank #eta_{#tildet}_{1}",          "Events", true, cutlabel, 5 );
        //plt.plot( "h_stop1Phi_MassRank_"+cutlabel,          "Mass Rank #phi_{#tildet}_{1}",          "Events", true, cutlabel, 5 );
        //plt.plot( "h_stop1Pt_MassRank_"+cutlabel,           "Mass Rank pT_{#tildet}_{1}",            "Events", true, cutlabel, 0, 1500, 5 );
        //plt.plot( "h_stop2Mass_MassRank_"+cutlabel,         "Mass Rank M_{#tildet}_{2}",             "Events", true, cutlabel, 0, 1500, 5 );
        //plt.plot( "h_stop2Eta_MassRank_"+cutlabel,          "Mass Rank #eta_{#tildet}_{2}",          "Events", true, cutlabel, 5 );
        //plt.plot( "h_stop2Phi_MassRank_"+cutlabel,          "Mass Rank #phi_{#tildet}_{2}",          "Events", true, cutlabel, 5 );
        //plt.plot( "h_stop2Pt_MassRank_"+cutlabel,           "Mass Rank pT_{#tildet}_{2}",            "Events", true, cutlabel, 0, 1500, 5 );

        //plt.plot( "h_stop1Mass_ScalarPtRank_"+cutlabel,     "ScalarPt Rank M_{#tildet}_{1}",         "Events", true, cutlabel, 0, 1500, 5 );
        //plt.plot( "h_stop1Eta_ScalarPtRank_"+cutlabel,      "ScalarPt Rank #eta_{#tildet}_{1}",      "Events", true, cutlabel, 5 );
        //plt.plot( "h_stop1Phi_ScalarPtRank_"+cutlabel,      "ScalarPt Rank #phi_{#tildet}_{1}",      "Events", true, cutlabel, 5 );
        //plt.plot( "h_stop1Pt_ScalarPtRank_"+cutlabel,       "ScalarPt Rank pT_{#tildet}_{1}",        "Events", true, cutlabel, 0, 1500, 5 );
        //plt.plot( "h_stop2Mass_ScalarPtRank_"+cutlabel,     "ScalarPt Rank M_{#tildet}_{2}",         "Events", true, cutlabel, 0, 1500, 5 );
        //plt.plot( "h_stop2Eta_ScalarPtRank_"+cutlabel,      "ScalarPt Rank #eta_{#tildet}_{2}",      "Events", true, cutlabel, 5 );
        //plt.plot( "h_stop2Phi_ScalarPtRank_"+cutlabel,      "ScalarPt Rank #phi_{#tildet}_{2}",      "Events", true, cutlabel, 5 );
        //plt.plot( "h_stop2Pt_ScalarPtRank_"+cutlabel,       "ScalarPt Rank pT_{#tildet}_{2}",        "Events", true, cutlabel, 0, 1500, 5 );

        //plt.plot( "h_stop1ScalarPt_ScalarPtRank_"+cutlabel, "ScalarPt Rank Scalar pT_{#tildet}_{1}", "Events", true, cutlabel, 0, 1500, 5 );
        //plt.plot( "h_stop2ScalarPt_ScalarPtRank_"+cutlabel, "ScalarPt Rank Scalar pT_{#tildet}_{2}", "Events", true, cutlabel, 0, 1500, 5 );

        //plt.plot( "h_MT2_"+cutlabel,                        "MT2",                                   "Events", true, cutlabel, 0, 1500, 5 );
        //plt.plot( "h_dR_stop1stop2_"+cutlabel,              "#DeltaR",                               "Events", true, cutlabel, 5 );
        //plt.plot( "h_dPhi_stop1stop2_"+cutlabel,            "#Delta#phi",                            "Events", true, cutlabel, 0, 4 );
        //plt.plot( "h_difference_stopMasses_"+cutlabel,      "Pt Rank difference",                    "Events", true, cutlabel, 0, 1500, 5 );
        //plt.plot( "h_average_stopMasses_"+cutlabel,         "average",                               "Events", true, cutlabel, 0, 1500, 5 );
        //plt.plot( "h_relativeDiff_stopMasses_"+cutlabel,    "Pt Rank relative difference",           "Events", true, cutlabel, 0, 100,  5 );

    }
}
