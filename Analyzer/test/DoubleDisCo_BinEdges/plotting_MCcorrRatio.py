import ROOT
import argparse
import os

# -----------------------------
# get canvas for all histograms
# -----------------------------
def get_canvas(year, channel):

    canvas = ROOT.TCanvas(year + "_MC_correction_ratio_" + channel, year + "_MC_correction_ratio_" + channel, 800, 800)
    canvas.cd()
    ROOT.gPad.SetGrid()
    ROOT.gPad.SetTopMargin(0.1)
    ROOT.gPad.SetBottomMargin(0.1)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetRightMargin(0.1)
    ROOT.gPad.SetTicks()

    return canvas

# -----------------------------
# get canvas for each histogram
# -----------------------------
def get_canvas_eachTTvar(year, ttvar, channel):
    
    canvas = ROOT.TCanvas(year + "_MC_correction_ratio_" + ttvar + "_" + channel, year + "_MC_correction_ratio_" + ttvar + "_" + channel, 800, 800)
    canvas.cd()
    canvas.Divide(1,2)

    split      = 0.3
    upperSplit = 1.0
    lowerSplit = 1.0

    upperSplit = 1.0-split
    lowerSplit = split
    scale = upperSplit / lowerSplit

    canvas.cd(1)
    ROOT.gPad.SetPad(0.0, split, 1.0, 1.0)
    ROOT.gPad.SetTopMargin(0.1 / upperSplit) 
    ROOT.gPad.SetBottomMargin(0)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetRightMargin(0.1)
    ROOT.gPad.SetGrid()
    ROOT.gPad.SetTicks()

    canvas.cd(2)
    ROOT.gPad.SetPad(0.0, 0.0, 1.0, 0.3)
    ROOT.gPad.SetTopMargin(0)
    ROOT.gPad.SetBottomMargin(0.1 / lowerSplit)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetRightMargin(0.1)
    ROOT.gPad.SetGrid()
    ROOT.gPad.SetTicks()

    return canvas, scale

# -----------------------------------
# get legend to superimpose all ttvar
# -----------------------------------
def get_legend(textsize=0.02):

    legend = ROOT.TLegend(0.2, 0.6, 0.5, 0.85, "", "trNDC")
    legend.SetNColumns(2)
    legend.SetFillStyle(0)
    legend.SetTextSize(textsize)
    legend.SetLineWidth(0)

    return legend

# --------------------------------------
# get legend to superimpose ttvar and tt
# --------------------------------------
def get_legend_eachTTvar(textsize=0.02):

    legend = ROOT.TLegend(0.2, 0.78, 0.5, 0.88, "", "trNDC")
    legend.SetFillStyle(0)
    legend.SetTextSize(textsize)
    legend.SetLineWidth(0)

    return legend

# ----------------------
# add CMS logo to canvas
# ----------------------
def addCMSlogo(canvas, year, TopMargin, LeftMargin, RightMargin, SF=1.0):

    canvas.cd()
    mark = ROOT.TLatex()
    mark.SetNDC(True)
    mark.SetTextAlign(11)
    mark.SetTextSize(0.048)
    mark.SetTextFont(61)
    mark.DrawLatex(LeftMargin, (1 - (TopMargin - 0.01)*SF), "CMS")
    mark.SetTextFont(52)
    mark.SetTextSize(0.038)
    mark.DrawLatex(LeftMargin + 0.12, (1 - (TopMargin - 0.01)*SF), "Work in Progress")
    mark.SetTextAlign(31)
    mark.DrawLatex(1 - RightMargin, (1 - (TopMargin - 0.01)*SF), "%s (13 TeV)"%(year))


# ---------------------
# main part of plotting
# ---------------------
def main():

    usage  = "usage: %prog [options]"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("--sig",  dest="sig",  help="signal model RPV, SYY", default="RPV")
    parser.add_argument("--mass", dest="mass", help="signal mass",           default="550")
    args = parser.parse_args()

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)

    labels = [ #"0l_0.6_0.6",
               #"1l_0.6_0.6", 
               #"2l_0.6_0.6",
            
               # SYY optimized ABCD bin edges
               "0l_0.85_0.74",
               "1l_0.7_0.85",
               "2l_0.69_0.57", 
             ]

    path = "year_TT_TTvar_Syst_%s_%s_label.root"%(args.sig, args.mass)
   
 
    years = [
             #"2016preVFP" ,
             #"2016postVFP" ,
             #"2017" ,
             #"2018" ,
             "Run2UL",
            ]
    
    ttvars = {
             "TT",
             "TT_fsrUp",            
             "TT_fsrDown",          
             "TT_isrUp",            
             "TT_isrDown",          
             "TT_hdampUP",          
             "TT_hdampDOWN",        
             "TT_TuneCP5up",
             "TT_TuneCP5down",
             "TT_erdON",
             "TT_JECup",            
             "TT_JECdown",          
             "TT_JERup",            
             "TT_JERdown",
            }
    
    ttvar_colors = {
                    "TT"             : "#000000",
                    "TT_erdON"       : "#878787",
                    "TT_fsrUp"       : "#238443",
                    "TT_fsrDown"     : "#78c679",
                    "TT_isrUp"       : "#88419d",
                    "TT_isrDown"     : "#8c96c6",
                    "TT_hdampUP"     : "#0570b0",
                    "TT_hdampDOWN"   : "#74a9cf",
                    "TT_JECup"       : "#ae017e",
                    "TT_JECdown"     : "#f768a1",
                    "TT_JERup"       : "#8c510a",
                    "TT_JERdown"     : "#dfc27d",
                    "TT_TuneCP5up"   : "#d73027",
                    "TT_TuneCP5down" : "#fdae61",
                   }
    
    ttvar_names = {
                   "TT"             : "Nominal TT",
                   "TT_erdON"       : "Color Reconn.",
                   "TT_fsrUp"       : "FSR Up",
                   "TT_fsrDown"     : "FSR Down",
                   "TT_isrUp"       : "ISR Up",
                   "TT_isrDown"     : "ISR Down",
                   "TT_hdampUP"     : "ME-PS Up",
                   "TT_hdampDOWN"   : "ME-PS Down",
                   "TT_JECup"       : "JEC Up",
                   "TT_JECdown"     : "JEC Down",
                   "TT_JERup"       : "JER Up",
                   "TT_JERdown"     : "JER Down",
                   "TT_TuneCP5up"   : "Tune Up",
                   "TT_TuneCP5down" : "Tune Down",
                  }
    
    # -------------------
    # loop over the years
    # -------------------
    for year in years:

        # create output directory for each year
        outputPath = "%s_plots_MCcorrectionFactorRatio/"%(year)
        if not os.path.exists(outputPath):
            os.makedirs(outputPath)

        # ----------------------
        # loop over the channels
        # ----------------------
        for label in labels:

            # open root files
            file1 = path.replace("year", year).replace("label", label)
            f1    = ROOT.TFile.Open(file1, "READ")
    
            # get canvas for all histograms
            canvas = get_canvas(year, label)
            
            # get legend
            legend = get_legend()
  
            draw = False 

            MCcorr_TT = f1.Get("%s_MCcorr_TT_TT"%(year))
            MCcorr_TT.SetTitle("")
            MCcorr_TT.SetLineWidth(4)
            MCcorr_TT.SetLineColor(ROOT.kBlack)
            MCcorr_TT.GetYaxis().SetTitle("Closure Correction [TT]")

            # ----------------
            # loop over ttvars
            # ----------------
            for ttvar in ttvars:

                if ttvar == "TT":
                    continue

                histName = "%s_MCcorr_Ratio_MC_%s"%(year,ttvar)
                yTitle   = "Closure Correction Ratio [TTvar / TT]"
                tag      = "_MC_correction_ratio_"
                
                if ttvar == "TT":
                    histName = "%s_MCcorr_TT_TT"%(year)
                    yTitle = "MC Correction [TT]"
                    tag = "_MC_correction_"

                # -------------------------
                # canvas for all histograms
                # -------------------------
                hist = f1.Get(histName)
                hist.SetTitle("")
                hist.SetLineWidth(4)
                #hist.SetLineColor(ttvar_colors[ttvar])
                hist.SetLineColor(ROOT.TColor.GetColor(ttvar_colors[ttvar]))
                hist.GetXaxis().SetTitle("N_{jets}")
                hist.GetYaxis().SetTitle(yTitle)

                MCcorr_TTvar = f1.Get("%s_MCcorr_TTvar_%s"%(year,ttvar))
                MCcorr_TTvar.SetTitle("")
                MCcorr_TTvar.SetLineWidth(4)
                MCcorr_TTvar.SetLineColor(ROOT.TColor.GetColor(ttvar_colors[ttvar]))
                MCcorr_TTvar.GetXaxis().SetTitle("N_{jets}")
                MCcorr_TTvar.GetYaxis().SetTitle("MC Correction")

                globalScale = 1.10
                xLabelSize = 0.05; yLabelSize = 0.05
                xTitleSize = 0.05; yTitleSize = 0.05

                hist.GetXaxis().SetLabelSize(xLabelSize * 0.8)
                hist.GetYaxis().SetLabelSize(yLabelSize * 0.8)
                hist.GetXaxis().SetTitleSize(xTitleSize * 0.8)
                hist.GetYaxis().SetTitleSize(yTitleSize * 0.8)
                hist.GetYaxis().SetTitleOffset(1.4)
                hist.GetYaxis().SetTitleOffset(1.4)

                histBottomPanel = hist.Clone(hist.GetName() + "_clone")

                if ttvar != "TT":
                    canvas.cd()
                    if not draw:
                        hist.GetYaxis().SetRangeUser(0.0, 2.0)
                        hist.Draw("hist")
                        draw = True
                        
                    else:
                        hist.Draw("hist SAME")

                    legend.AddEntry(hist, ttvar_names[ttvar], "l")

                # ----------------------------------
                # canvas & legend for each histogram
                # ----------------------------------
                canvas_each, scale = get_canvas_eachTTvar(year, ttvar, label)
                legend_each = get_legend_eachTTvar(0.04)

                MCcorr_TT.GetXaxis().SetLabelSize(xLabelSize * globalScale)
                MCcorr_TT.GetYaxis().SetLabelSize(yLabelSize * globalScale)
                MCcorr_TT.GetXaxis().SetTitleSize(xTitleSize * globalScale)
                MCcorr_TT.GetYaxis().SetTitleSize(yTitleSize * globalScale)
                MCcorr_TT.GetXaxis().SetTitleOffset(1)
                MCcorr_TT.GetYaxis().SetTitleOffset(1)

                canvas_each.cd(1)

                end = MCcorr_TT.GetNbinsX()

                line = ROOT.TLine(0.0, 1.0, float(end), 1.0)
                line.SetLineColor(ROOT.kBlack)
                line.SetLineWidth(2)
                line.SetLineStyle(7)

                MCcorr_TT.Draw("hist E")
                MCcorr_TTvar.Draw("hist E SAME")
                line.Draw("same")

                legend_each.AddEntry(MCcorr_TT, "Nominal TT", "l")            
                legend_each.AddEntry(MCcorr_TTvar, ttvar_names[ttvar], "l")            
                addCMSlogo(canvas_each, year, TopMargin=0.1, LeftMargin=0.12, RightMargin=0.1, SF=1.0)
                legend_each.Draw("SAME")

                canvas_each.cd(2)
                histBottomPanel.GetXaxis().SetLabelSize(xLabelSize * scale * globalScale)
                histBottomPanel.GetYaxis().SetLabelSize(yLabelSize * scale * globalScale)
                histBottomPanel.GetXaxis().SetTitleSize(xTitleSize * scale * globalScale)
                histBottomPanel.GetYaxis().SetTitleSize(yTitleSize * scale * globalScale)
                histBottomPanel.GetYaxis().SetTitleOffset(0.3/0.7)

                histBottomPanel.GetYaxis().SetTitleOffset(histBottomPanel.GetYaxis().GetTitleOffset() / globalScale)

                histBottomPanel.GetYaxis().SetRangeUser(0.04,1.96) 
                histBottomPanel.GetYaxis().SetNdivisions(7)
                histBottomPanel.GetYaxis().SetTitle("Ratio [TTvar / TT]")
                histBottomPanel.Draw("hist E")

                line.Draw("same")
                canvas_each.SaveAs("%s_plots_MCcorrectionFactorRatio/"%(year) + year + "_" + args.sig + "_" + args.mass + tag + ttvar + "_" + label + ".pdf")
 
            # ---------------------------------------------
            # save canvas & legend including all histograms
            # ---------------------------------------------
            canvas.cd()
            addCMSlogo(canvas, year, TopMargin=0.1, LeftMargin=0.12, RightMargin=0.1, SF=1.0)    
            legend.Draw("SAME")
            canvas.SaveAs("%s_plots_MCcorrectionFactorRatio/"%(year) + year + "_" + args.sig + "_" + args.mass + "_MCcorr_Ratio_MC_" + label + ".pdf")     
         
        f1.Close()
 
if __name__ == "__main__":
    main()


         
