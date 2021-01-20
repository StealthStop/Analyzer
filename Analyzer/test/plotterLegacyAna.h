#ifndef PLOTTERLEGACYANA_H
#define PLOTTERLEGACYANA_H

#include "TH1.h"
#include "THStack.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TF1.h"
#include "TString.h"
#include "TPaveText.h"
#include "HistInfoCollectionLegacyAna.h"
#include "Framework/Framework/include/Utility.h"

#include <memory>
#include <vector>
#include <string>
#include <cstdio>
#include <iostream>
#include <map>

class Plotter
{
private:
    //Collection of histInfo
    HistInfoCollection hc_;
    std::shared_ptr<TH1> hbgSum_;
    std::map< std::string, HistInfoCollection > mhc_;
    std::vector<std::shared_ptr<TH1>> hbgSumVec_;
    std::string outpath_;

public:
    Plotter(HistInfoCollection&& hc, const std::string& outpath = "outputPlots") : hc_(hc), outpath_(outpath){}
    Plotter(std::map< std::string, HistInfoCollection >&& mhc, const std::string& outpath = "outputPlots") : mhc_(mhc), outpath_(outpath){}

    void plotStack(const std::string& histName, const std::string& xAxisLabel, const std::string& yAxisLabel = "Events", const bool isLogY = false, int rebin = -1, const bool scale = false, const bool doFill = true, const double xmin = 999.9, const double xmax = -999.9, double lumi = 36100, bool supplementary = false, bool approved = true, double padDiv = 0.3)
    {
        //This is a magic incantation to disassociate opened histograms from their files so the files can be closed
        TH1::AddDirectory(false);

        double xTitleSize = 0.065; double yTitleSize = 0.065;
        double xLabelSize = 0.045;  double yLabelSize = 0.045;
        double xOffset = 2.0;      double yOffset = 1.1;
        double sf = (1.0-padDiv)/padDiv;
        //create the canvas for the plot
        TCanvas *c = new TCanvas("c1", "c1", 800, 800);
        //switch to the canvas to ensure it is the active object
        c->cd();

        // Upper plot will be in pad1: TPad(x1, y1, x2, y2)
        TPad *pad1 = new TPad("pad1", "pad1", 0.0, padDiv, 1.0, 1.0);
        pad1->SetBottomMargin(0.); // Upper and lower plot are joined
        pad1->SetLeftMargin(0.14); // Upper and lower plot are joined
        pad1->SetRightMargin(0.05);
        pad1->Draw();             // Draw the upper pad: pad1
        pad1->cd();               // pad1 becomes the current pad
        
        //Create TLegend
        TLegend *leg = nullptr;
        if (not supplementary) leg = new TLegend(0.17, 0.50, 0.35, 0.88);
        else leg = new TLegend(0.17, 0.46, 0.35, 0.84);
        leg->SetFillStyle(0);
        leg->SetColumnSeparation(0.15);
        leg->SetBorderSize(0);
        leg->SetLineWidth(1);
        leg->SetNColumns(1);
        leg->SetTextFont(42);
        leg->SetTextSize(0.040);

        TLegend *dataleg = nullptr;
        if (not supplementary) dataleg = new TLegend(0.26, 0.50, 0.40, 0.88);
        else dataleg = new TLegend(0.26, 0.46, 0.40, 0.84);
        dataleg->SetFillStyle(0);
        dataleg->SetColumnSeparation(0.15);
        dataleg->SetBorderSize(0);
        dataleg->SetLineWidth(1);
        dataleg->SetNColumns(1);
        dataleg->SetTextFont(42);
        dataleg->SetTextSize(0.040);

        TLegend *sigleg = nullptr;
        if (not supplementary) sigleg = new TLegend(0.40, 0.56, 0.84, 0.88);
        else sigleg = new TLegend(0.56, 0.56, 0.82, 0.85);
        sigleg->SetFillStyle(0);
        sigleg->SetColumnSeparation(0.15);
        sigleg->SetBorderSize(0);
        sigleg->SetLineWidth(1);
        sigleg->SetNColumns(1);
        sigleg->SetTextFont(42);
        sigleg->SetTextSize(0.040);

        // ------------------------
        // -  Setup plots
        // ------------------------
        THStack *bgStack = new THStack();
        hc_.setUpBG(histName, rebin, bgStack, hbgSum_, doFill, scale);

        hc_.setUpSignal(histName, rebin);

        hc_.setUpData(histName, rebin);

        hc_.setUpSyst(histName);

        //create a dummy histogram to act as the axes
        histInfo dummy(new TH1D("dummy", "dummy", 1000, hbgSum_->GetBinLowEdge(1), hbgSum_->GetBinLowEdge(hbgSum_->GetNbinsX()) + hbgSum_->GetBinWidth(hbgSum_->GetNbinsX())));        

        //draw dummy axes
        dummy.setupAxes(xOffset, yOffset, xTitleSize, yTitleSize, xLabelSize, yLabelSize);
        dummy.draw();

        //Switch to logY if desired
        gPad->SetLogy(isLogY);
        gPad->SetTicks(1,1);

        //get maximum from histos and fill TLegend
        double min = 0.0;
        double max = 0.0;
        double lmax = 0.0;

        // -----------------------
        // -  Plot Background
        // -----------------------
        TString legType = (doFill) ? "F" : "L";
        if (hc_.bgVec_.size() != 0) {
            for(auto& entry : hc_.bgVec_)
            {
                leg->AddEntry(entry.h.get(), entry.legEntry.c_str(), legType);
            }
            smartMax(hbgSum_.get(), leg, static_cast<TPad*>(gPad), min, max, lmax, false);
            bgStack->Draw("same");
        }

        // -------------------------
        // -   Plot Signal
        // -------------------------
        if (hc_.sigVec_.size() != 0) {
            for(const auto& entry : hc_.sigVec_)
            {
                sigleg->AddEntry(entry.h.get(), entry.legEntry.c_str(), "L");
                smartMax(entry.h.get(), leg, static_cast<TPad*>(gPad), min, max, lmax, false);
                entry.draw();
            }
        }

        //------------------------
        //-  Plot Data
        //------------------------
        if (hc_.dataVec_.size() != 0) {
            dataleg->AddEntry((TObject*)0, "", "");
            dataleg->AddEntry((TObject*)0, "", "");
            dataleg->AddEntry((TObject*)0, "", "");

            for(const auto& entry : hc_.dataVec_)
            {
                legType = (entry.drawOptions=="hist") ? "L" : entry.drawOptions;
                dataleg->AddEntry(entry.h.get(), entry.legEntry.c_str(), legType);
                smartMax(entry.h.get(), leg, static_cast<TPad*>(gPad), min, max, lmax, true);
                entry.draw();
                entry.h->Draw("E0 SAME");
            }

        }

        if (hc_.systVec_.size() != 0) {
            for(const auto& entry : hc_.systVec_)
            {
                smartMax(entry.h.get(), leg, static_cast<TPad*>(gPad), min, max, lmax, true);
                entry.ge->Draw("2SAME");
            }
        }


        //plot legend
        leg->Draw("same");
        dataleg->Draw("same");
        sigleg->Draw("same");

        TLatex significance;  
        significance.SetNDC(true);
        significance.SetTextAlign(11);
        significance.SetTextFont(52);
        significance.SetTextSize(0.030);

        //Draw CMS and lumi lables
        drawLables(lumi, supplementary, approved);
            
        //Draw dummy hist again to get axes on top of histograms
        setupDummy(dummy, sigleg, "", yAxisLabel, isLogY, xmin, xmax, min, max, lmax);
        dummy.setupAxes(xOffset, yOffset, xTitleSize, yTitleSize, xLabelSize, yLabelSize);

        dummy.draw("AXIS");
                        
        // lower plot will be in pad2
        c->cd();          // Go back to the main canvas before defining pad2
        TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1.0, padDiv);
        pad2->SetTopMargin(0.);
        pad2->SetBottomMargin(0.30);
        pad2->SetRightMargin(0.05);
        pad2->SetLeftMargin(0.14);
        pad2->SetGridy(); // Horizontal grid
        pad2->Draw();
        pad2->cd();       // pad2 becomes the current pad        
        gPad->SetTicks(1,1);
            
        //make ratio dummy
        histInfo ratioDummy(new TH1D("rdummy", "rdummy", 1000, hc_.bgVec_[0].h->GetBinLowEdge(1), 
                                     hc_.bgVec_[0].h->GetBinLowEdge(hc_.bgVec_[0].h->GetNbinsX()) + hc_.bgVec_[0].h->GetBinWidth(hc_.bgVec_[0].h->GetNbinsX())));
        ratioDummy.h->GetXaxis()->SetTitle(xAxisLabel.c_str());
        ratioDummy.h->GetYaxis()->SetTitle("Data / MC");

        ratioDummy.h->GetXaxis()->SetTickLength(0.1);
        ratioDummy.h->GetYaxis()->SetTickLength(0.045);
        ratioDummy.setupAxes(xOffset/sf, yOffset/sf, xTitleSize*sf, yTitleSize*sf, xLabelSize*sf, yLabelSize*sf);
        ratioDummy.h->GetYaxis()->SetNdivisions(3,5,0);
        ratioDummy.h->GetXaxis()->SetRangeUser(xmin, xmax);
        ratioDummy.h->GetYaxis()->SetRangeUser(0.0, 5.0);
        ratioDummy.h->SetStats(0);
        ratioDummy.h->SetMinimum(0.4);
        ratioDummy.h->SetMaximum(1.6);
            
        //Make ratio histogram for data / background.
        histInfo ratio((TH1*)hc_.dataVec_[0].h->Clone());
                
        ratio.drawOptions = "ep";
        ratio.color = kBlack;
            
        ratio.h->Divide(hbgSum_.get());
        ratio.h->SetMarkerStyle(20);
            
        ratioDummy.draw();
        ratio.draw("same");
            
        if (hc_.rsystVec_.size() != 0) {
            for(const auto& entry : hc_.rsystVec_)
            {
                entry.ge->Draw("2SAME");
            }
        }

        ratio.draw("same");
        ratio.h->Draw("E0 SAME");

        pad1->cd();

        //save new plot to file
        if (not approved) {
            c->Print((outpath_ + "/" + histName + "_prelim.pdf").c_str());
        } else {
            c->Print((outpath_ + "/" + histName + ".pdf").c_str());
        }

        //clean up dynamic memory
        delete c;
        delete leg;
        delete sigleg;
        delete bgStack;
    }

    //This is a helper function which will keep the plot from overlapping with the legend
    void smartMax(const TH1* const h, const TLegend* const l, const TPad* const p, double& gmin, double& gmax, double& gpThreshMax, const bool error)
    {
        const bool isLog = p->GetLogy();
        double min = 9e99;
        double max = -9e99;
        double pThreshMax = -9e99;
        int threshold = static_cast<int>(h->GetNbinsX()*(l->GetX1() - p->GetLeftMargin())/((1 - p->GetRightMargin()) - p->GetLeftMargin()));

        for(int i = 1; i <= h->GetNbinsX(); ++i)
        {
            double bin = 0.0;
            if(error) bin = h->GetBinContent(i) + h->GetBinError(i);
            else      bin = h->GetBinContent(i);
            if(bin > max) max = bin;
            else if(bin > 1e-10 && bin < min) min = bin;
            if(i >= threshold && bin > pThreshMax) pThreshMax = bin;
        }
        
        gpThreshMax = std::max(gpThreshMax, pThreshMax);
        if (isLog) gmax = 10*std::max(gmax, max);
        else gmax = 1.06*std::max(gmax, max);
        gmin = std::min(gmin, min);
    }

    void setupDummy(histInfo dummy, TLegend *leg, const std::string& xAxisLabel, const std::string& yAxisLabel, const bool isLogY, 
                    const double xmin, const double xmax, double min, double max, double lmax)
    {
        dummy.setupAxes(1.2, 1.4, 0.065, 0.065, 0.045, 0.045);
        dummy.h->GetYaxis()->SetTitle(yAxisLabel.c_str());
        dummy.h->GetXaxis()->SetTitle(xAxisLabel.c_str());
        dummy.h->SetTitle("");
        //Set the y-range of the histogram
        if(isLogY)
        {
            double locMin = std::min(0.2, std::max(0.2, 0.05 * min));
            double legSpan = (log10(3*max) - log10(locMin)) * (leg->GetY1() - gPad->GetBottomMargin()) / ((1 - gPad->GetTopMargin()) - gPad->GetBottomMargin());
            double legMin = legSpan + log10(locMin);
            if(log10(lmax) > legMin)
            {
                double scale = (log10(lmax) - log10(locMin)) / (legMin - log10(locMin));
                max = pow(max/locMin, scale)*locMin;
            }
            dummy.h->GetYaxis()->SetRangeUser(locMin, 10*max);
        }
        else
        {
            double locMin = 0.0;
            double legMin = (1.2*max - locMin) * (leg->GetY1() - gPad->GetBottomMargin()) / ((1 - gPad->GetTopMargin()) - gPad->GetBottomMargin());
            if(lmax > legMin) max *= (lmax - locMin)/(legMin - locMin);
            dummy.h->GetYaxis()->SetRangeUser(1.0, max*1.4);
        }
        //set x-axis range
        if(xmin < xmax) dummy.h->GetXaxis()->SetRangeUser(xmin, xmax);
    }

    void drawLables(double lumi, bool supplementary = false, bool approved = false)
    {
        //Draw CMS and lumi lables
        char lumistamp[128];
        if (lumi > 100000) sprintf(lumistamp, "%.0f fb^{-1} (13 TeV)", lumi / 1000.0);
        else sprintf(lumistamp, "%.1f fb^{-1} (13 TeV)", lumi / 1000.0);

        TLatex mark;
        mark.SetNDC(true);

        //Draw CMS mark
        mark.SetTextAlign(11);
        mark.SetTextSize(0.065);
        mark.SetTextFont(61);
        mark.DrawLatex(gPad->GetLeftMargin(), 1 - (gPad->GetTopMargin() - 0.017), "CMS");
        mark.SetTextFont(52);
        if (not approved and not supplementary) mark.DrawLatex(gPad->GetLeftMargin() + 0.1025, 1 - (gPad->GetTopMargin() - 0.017), "Preliminary");
        else if (supplementary) mark.DrawLatex(gPad->GetLeftMargin() + 0.1025, 1 - (gPad->GetTopMargin() - 0.017), "Supplementary");
        else mark.DrawLatex(gPad->GetLeftMargin() + 0.1025, 1 - (gPad->GetTopMargin() - 0.017), "");

        //Draw lumistamp
        mark.SetTextFont(42);
        mark.SetTextAlign(31);
        mark.DrawLatex(1 - gPad->GetRightMargin(), 1 - (gPad->GetTopMargin() - 0.017), lumistamp);        

        if (supplementary)
        {
            mark.SetTextSize(0.04);
            mark.DrawLatex(gPad->GetLeftMargin() + 0.24, 1 - (gPad->GetTopMargin() + 0.055), "arXiv XXXX.XXXX");
        }
    }    
};

#endif
