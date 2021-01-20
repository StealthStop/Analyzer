#ifndef HISTINFOCOLLECTIONLEGACYANA_H
#define HISTINFOCOLLECTIONLEGACYANA_H

#include "TH1.h"
#include "THStack.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"

#include <memory>
#include <vector>
#include <string>
#include <cstdio>
#include <iostream>
#include <map>

//Class to hold TH1* with various helper functions 
class histInfo
{
public:
    std::string legName, histFile, histName, drawOptions, legEntry;
    int color, rebin;
    double scale; 
    bool norm2Data;
    int lineStyle;
    std::shared_ptr<TH1> h;
    std::shared_ptr<TGraphErrors> ge;

    // Merge the njets >= 12 bins together
    void fixNJetsHisto(std::shared_ptr<TH1>& njetsHisto)
    {

        int nbinsX = njetsHisto->GetNbinsX();
        int bin12jet = njetsHisto->FindBin(12.2);
        for (int bin = bin12jet+1; bin < nbinsX; bin++) {

        
            double content = njetsHisto->GetBinContent(bin);
            njetsHisto->Fill(12.2, content);

            njetsHisto->SetBinContent(bin, 0.0);
            njetsHisto->SetBinError(bin, 0.0);

        }
    }

    //helper function to get histogram from file and configure its optional settings
    void retrieveHistogram()
    {
        //Open the file for this histogram
        TFile *f = TFile::Open(histFile.c_str());

        //check that the file was opened successfully
        if(!f)
        {
            printf("File \"%s\" could not be opened!!!\n", histFile.c_str());
            h = nullptr;
            return;
        }

        //get the histogram from the file
        h.reset(static_cast<TH1*>(f->Get(histName.c_str())));

        ge.reset(new TGraphErrors(static_cast<TH1*>(f->Get(histName.c_str()))));

        if (histName.find("njets") != std::string::npos) { fixNJetsHisto(h); }

        //with the histogram retrieved, close the file
        f->Close();
        delete f;

        //check that the histogram was retireved from the file successfully
        if(!h)
        {
            printf("Histogram \"%s\" could not be found in file \"%s\"!!!\n", histName.c_str(), histFile.c_str());
            return;
        }

        //set the histogram color
        h->SetLineColor(color);
        h->SetLineWidth(3);
        h->SetMarkerColor(color);
        h->SetMarkerStyle(20);
        h->Scale(scale);

        ge->SetLineWidth(0);
        ge->SetMarkerSize(0);
        ge->SetFillStyle(3254);
        ge->SetFillColor(kBlack);
        
        legEntry = legName;

        // rebin the histogram if desired
        if(rebin >0) h->Rebin(rebin);
    }

    //helper function for axes
    void setupAxes(double xOffset, double yOffset, double xTitle, double yTitle, double xLabel, double yLabel)
    {
        h->SetStats(0);
        h->SetTitle(0);
        h->GetXaxis()->SetTitleOffset(xOffset);
        h->GetYaxis()->SetTitleOffset(yOffset);
        h->GetXaxis()->SetTitleSize(xTitle);
        h->GetYaxis()->SetTitleSize(yTitle);
        h->GetXaxis()->SetLabelSize(xLabel);
        h->GetYaxis()->SetLabelSize(yLabel);
        if(h->GetXaxis()->GetNdivisions() % 100 > 5) h->GetXaxis()->SetNdivisions(6, 5, 0);
    }

    //helper function for pads
    void setupPad(double left, double right, double top, double bottom)
    {
        gPad->SetLeftMargin(left);
        gPad->SetRightMargin(right);
        gPad->SetTopMargin(top);
        gPad->SetBottomMargin(bottom);
        gPad->SetTicks(1,1);
    }

    //wrapper to draw histogram
    void draw(const std::string& additionalOptions = "", bool noSame = false) const
    {
        h->Draw(((noSame?"":"same " + drawOptions + " " + additionalOptions)).c_str());
    }

    void setFillColor(int newColor = -1)
    {
        if(newColor >= 0) h->SetFillColor(newColor);
        else              h->SetFillColor(color);
    }

    void setLineStyle(int style = 1)
    {
        h->SetLineStyle(style);
    }

    histInfo(const std::string& legName, const std::string& histFile, const std::string& drawOptions, const int color, const double scale = 1.0, const bool norm2Data = false, int lineStyle = 0) : legName(legName), histFile(histFile), histName(""), drawOptions(drawOptions), color(color), rebin(-1), scale(scale), norm2Data(norm2Data), lineStyle(lineStyle), h(nullptr)
    {
    }

    histInfo(TH1* h) : legName(h->GetName()), histFile(""), histName(h->GetName()), drawOptions(""), color(kWhite), rebin(0), scale(1.0), h(h)
    {
    }

    ~histInfo()
    {
    }
};

class HistInfoCollection
{
public:
    std::vector<histInfo> dataVec_, bgVec_, sigVec_, systVec_, rsystVec_;

    HistInfoCollection(std::vector<histInfo>&& dataVec, std::vector<histInfo>&& bgVec, std::vector<histInfo>&& sigVec, std::vector<histInfo>&& systVec, std::vector<histInfo>&& rsystVec) : dataVec_(dataVec), bgVec_(bgVec), sigVec_(sigVec), systVec_(systVec), rsystVec_(rsystVec) {}
    HistInfoCollection(std::vector<histInfo>& dataVec, std::vector<histInfo>& bgVec, std::vector<histInfo>& sigVec, std::vector<histInfo>& systVec, std::vector<histInfo>& rsystVec) : dataVec_(dataVec), bgVec_(bgVec), sigVec_(sigVec), systVec_(systVec), rsystVec_(rsystVec) {}
    HistInfoCollection() {};

    ~HistInfoCollection() {}

    void setUpBG(const std::string& histName, int rebin, THStack* bgStack, std::shared_ptr<TH1>& hbgSum, const bool& setFill = true, const bool& scale = false)
    {
        if (dataVec_.size() != 0 and bgVec_.size() != 0) {
            double dNum = 1, bNum = 1;
            if(scale)
            {
                dNum = 0;
                bNum = 0;
                for(auto& entry : dataVec_)
                {
                    entry.histName = histName;
                    entry.rebin = rebin;
                    entry.retrieveHistogram();            
                    dNum += entry.h->Integral(); 
                }
                for(auto& entry : bgVec_)
                {
                    entry.histName = histName;
                    entry.rebin = rebin;
                    entry.retrieveHistogram();
                    bNum += entry.h->Integral();
                } 
            }
            
            bool firstPass = true;
            for(auto& entry : bgVec_)
            {
                entry.histName = histName;
                entry.rebin = rebin;
                entry.retrieveHistogram();
                if (entry.norm2Data) entry.h->Scale(dNum/bNum);

                bgStack->Add(entry.h.get(), entry.drawOptions.c_str());
                if(firstPass) 
                {
                    hbgSum.reset( static_cast<TH1*>(entry.h->Clone()) );
                    firstPass = false;
                }
                else 
                {
                    hbgSum->Add(entry.h.get());
                }
                
                if(setFill)
                {
                    entry.setFillColor();
                }
            }
        }
    }

    void setUpSignal(const std::string& histName, int rebin)
    {
        if (sigVec_.size() != 0) {

            for(auto& entry : sigVec_)
            {
                entry.histName = histName;
                entry.rebin = rebin;
                entry.retrieveHistogram();
                entry.h->Scale(entry.scale);
                entry.setLineStyle(entry.lineStyle);    
            }
        }
    }

    void setUpData(const std::string& histName, int rebin)
    {
        if (dataVec_.size() != 0) {
            for(auto& entry : dataVec_)
            {
                entry.histName = histName;
                entry.rebin = rebin;
                entry.retrieveHistogram();
            }
        }
    }

    void setUpSyst(const std::string& histName)
    {
        if (systVec_.size() != 0) {

            for(auto& entry : systVec_)
            {
                entry.histName = histName;
                entry.retrieveHistogram();
            }

            for(auto& entry : rsystVec_)
            {
            entry.histName = histName;
            entry.retrieveHistogram();
            }
        }
    }
};

#endif
