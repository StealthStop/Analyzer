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
    ROOT.gPad.SetTopMargin(0.06)
    ROOT.gPad.SetBottomMargin(0.1)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetRightMargin(0.03)
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
    ROOT.gPad.SetTopMargin(0.06 / upperSplit) 
    ROOT.gPad.SetBottomMargin(0)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetRightMargin(0.03)
    ROOT.gPad.SetGrid()
    ROOT.gPad.SetTicks()

    canvas.cd(2)
    ROOT.gPad.SetPad(0.0, 0.0, 1.0, 0.3)
    ROOT.gPad.SetTopMargin(0)
    ROOT.gPad.SetBottomMargin(0.1 / lowerSplit)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetRightMargin(0.03)
    ROOT.gPad.SetGrid()
    ROOT.gPad.SetTicks()

    return canvas, scale

# -----------------------------------
# get legend to superimpose all ttvar
# -----------------------------------
def get_legend(textsize=0.02):

    legend = ROOT.TLegend(0.40, 0.65, 0.88, 0.92, "", "trNDC")
    legend.SetNColumns(2)
    legend.SetFillStyle(0)
    legend.SetTextSize(textsize)
    legend.SetLineWidth(0)

    return legend

# --------------------------------------
# get legend to superimpose ttvar and tt
# --------------------------------------
def get_legend_eachTTvar(textsize=0.02):

    legend = ROOT.TLegend(0.36, 0.78, 0.71, 0.88, "", "trNDC")
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
    mark.DrawLatex(LeftMargin + 0.11, (1 - (TopMargin - 0.01)*SF), "Work in Progress")
    mark.SetTextAlign(31)
    mark.SetTextFont(42)
    mark.DrawLatex(1 - RightMargin, (1 - (TopMargin - 0.01)*SF), "%s (13 TeV)"%(year))

# ---------------------
# add extra information
# ---------------------
def addExtraInfo(canvas, LeftMargin, TopMargin, textsize, model, channel, div):

    canvas.cd(1)
    text = ROOT.TLatex()
    text.SetNDC(True)
    text.SetTextAlign(13)
    text.SetTextSize(textsize)
    text.SetTextFont(42)
    text.SetTextColor(ROOT.TColor.GetColor("#7C99D1"))

    name = ""
    if   channel == "0l":
        name = "0L"
    elif channel == "1l":
        name = "1L"
    elif channel == "2l":
        name = "2L"

    modelName = ""
    if "SYY" in model:
        modelName = "Stealth SYY"
    else:
        modelName = model
    
    divName = ""
    if "NonOptimized" in div:
        divName = "NonOptimized "
    elif "MassExclusion" in div:
        divName = "High Mass"
    elif "MaxSign" in div:
        divName = "Low Mass"
    
    text.DrawLatex(LeftMargin + 0.03, TopMargin - 0.04, modelName)
    text.DrawLatex(LeftMargin + 0.03, TopMargin - 0.08, name     )
    text.DrawLatex(LeftMargin + 0.03, TopMargin - 0.12, divName  )

# ------------------------------------
# relabel the njets bins in the x-axis
# ------------------------------------
def relabelNjetsBins(hist, channel, labelSize):

    if channel == "0l":
        njets = [8, 12]

    elif channel == "1l":
        njets = [7, 11]

    elif channel == "2l":
        njets = [6, 10]

    for bin in range(1, hist.GetNbinsX()+1):
        
        label = str(njets[0] + bin - 1)
        
        if label == str(njets[-1]):
            label = "#geq " + str(njets[0] + bin - 1)

        hist.GetXaxis().SetBinLabel(bin, label)
        hist.GetXaxis().SetLabelSize(labelSize)

# ---------------------
# main part of plotting
# ---------------------
def main():

    usage  = "usage: %prog [options]"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("--sig",     dest="sig",     help="signal model RPV, SYY",               default="RPV"                     )
    parser.add_argument("--mass",    dest="mass",    help="signal mass",                         default="550"                     )
    parser.add_argument("--div",     dest="div",     help="NonOptimized, MassEclusion, MaxSign", default="MassExclusion"           )
    parser.add_argument("--cat",   dest="cat",   help="1=Event Weight, 2=Independent Samples, 3=Jet Based", default=""           )
    parser.add_argument("--outpath", dest="outpath", help="output dir where group dirs",         default="Run2UL_MassExclusion_RPV") 
    
    args = parser.parse_args()

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)

    # Non-optimized edges (0.6, 0.6)
    if args.div == "NonOptimized" and args.sig == "RPV":
        labels = ["0l_0.6_0.6",
                  "1l_0.6_0.6",
                  "2l_0.6_0.6",
                 ]
    
    if args.div == "NonOptimized" and args.sig == "SYY":
        labels = ["0l_0.6_0.6",
                  "1l_0.6_0.6",
                  "2l_0.6_0.6",
                 ]

    # Bin edges for Mass Exclusion 
    if args.div == "MassExclusion" and args.sig == "RPV": 
        labels = ["0l_0.74_0.8", 
                  "1l_0.8_0.72", 
                  "2l_0.5_0.5",
                 ]

    if args.div == "MassExclusion" and args.sig == "SYY":
        labels = ["0l_0.54_0.56", 
                  "1l_0.68_0.82", 
                  "2l_0.48_0.48",
                 ]

    # Bin edges for Max Significance
    if args.div == "MaxSign" and args.sig == "RPV": 
        labels = ["0l_0.52_0.54", 
                  "1l_0.84_0.42", 
                  "2l_0.52_0.58",
                 ]

    if args.div == "MaxSign" and args.sig == "SYY": 
        labels = ["0l_0.76_0.7", 
                  "1l_0.44_0.42", 
                  "2l_0.4_0.42",
                 ]


    if args.sig == "SYY":
        args.outpath = args.outpath.replace("SYY", "StealthSYY")
    path = "%s/year_TT_TTvar_Syst_%s_%s_label.root"%(args.outpath, args.sig, args.mass)
   
    print(path)
 
    years = ["Run2UL"]
   
    if args.cat == "1":

        ttvars = {
                 "TT",
                 "TT_fsrUp",            
                 "TT_fsrDown",          
                 "TT_isrUp",            
                 "TT_isrDown",          
                }

    elif args.cat == "2": 

        ttvars = {
                 "TT",
                 "TT_hdampUP",          
                 "TT_hdampDOWN",        
                 "TT_TuneCP5up",
                 "TT_TuneCP5down",
                 "TT_erdON",
                }
    
    elif args.cat == "3":

        ttvars = {
                 "TT",
                 "TT_JECup",            
                 "TT_JECdown",          
                 "TT_JERup",            
                 "TT_JERdown",
                }

    else:
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
        outputPath = "%s/%s_plots_MCcorrectionFactorRatio/"%(args.outpath, year)
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
            legend = get_legend(0.028)
  
            draw = False 

            # get closure correction histogram to put plot including tt and all other ttvar
            MCcorr_TT = f1.Get("%s_MCcorr_TT_TT"%(year))
            MCcorr_TT.SetTitle("")
            MCcorr_TT.SetLineWidth(4)
            MCcorr_TT.SetLineColor(ROOT.kBlack)
            MCcorr_TT.GetYaxis().SetTitle("Closure Correction Ratio (Var/Nominal)")

            # get closure correction histogram to put comparion plot (each ttvar vs tt)
            # calculate statistical unc. on closure coerection
            MCcorr_TT_Unc = f1.Get("Run2UL_maximum_MCcorrectedData_Syst_All")
            #MCcorr_TT_Unc = MCcorr_TT.Clone("MCcorr_TT_Unc")
          
            unc_val = MCcorr_TT_Unc.GetBinContent(1)  
            for i in range(1, MCcorr_TT_Unc.GetNbinsX()+1):
                content        = MCcorr_TT.GetBinContent(i) 
                contentError   = MCcorr_TT.GetBinError(i)
                #closureCorrUnc = (contentError / content) 
                closureSyst = (abs(unc_val - 1) / content)
                
                MCcorr_TT_Unc.SetBinContent(i, 1.0)
                #MCcorr_TT_Unc.SetBinError(i, closureCorrUnc)
                MCcorr_TT_Unc.SetBinError(i, closureSyst)

            MCcorr_TT_Unc.SetTitle("")
            MCcorr_TT_Unc.SetLineWidth(0)
            MCcorr_TT_Unc.SetFillColorAlpha(ROOT.kGray, 0.8)
            MCcorr_TT_Unc.GetYaxis().SetTitle("#kappa_{var.} / #kappa_{nom.}")
            MCcorr_TT_Unc.GetXaxis().SetTitle("N_{jets}")

            relabelNjetsBins(MCcorr_TT_Unc, label.partition("_")[0], 1.1 * 0.05)

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
                hist.GetXaxis().SetTitle("N_{ jets}")
                hist.GetYaxis().SetTitle(yTitle)

                MCcorr_TTvar = f1.Get("%s_MCcorr_TTvar_%s"%(year,ttvar))
                MCcorr_TTvar.SetTitle("")
                MCcorr_TTvar.SetLineWidth(4)
                MCcorr_TTvar.SetLineColor(ROOT.TColor.GetColor(ttvar_colors[ttvar]))
                MCcorr_TTvar.GetXaxis().SetTitle("N_{ jets}")
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
                        MCcorr_TT_Unc.GetYaxis().SetRangeUser(0.0, 2.4)
                        MCcorr_TT_Unc.Draw("E2")
                        hist.Draw("hist SAME")
                        draw = True
        
                        legend.AddEntry(MCcorr_TT_Unc, "Stat. + Syst. Unc.", "F")
                        
                    else:
                        hist.Draw("hist SAME")

                    
                    legend.AddEntry(hist, ttvar_names[ttvar], "l")

                # ----------------------------------
                # canvas & legend for each histogram
                # ----------------------------------
                canvas_each, scale = get_canvas_eachTTvar(year, ttvar, label)
                legend_each = get_legend_eachTTvar(0.042)

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
                addCMSlogo(canvas_each, year, TopMargin=0.06, LeftMargin=0.12, RightMargin=0.03, SF=1.0)
                addExtraInfo(canvas_each, 0.12, 0.91, 0.050, args.sig, label.partition("_")[0], args.div)

                legend_each.Draw("SAME")

                canvas_each.cd(2)
                histBottomPanel.GetXaxis().SetLabelSize(xLabelSize * scale * globalScale)
                histBottomPanel.GetYaxis().SetLabelSize(yLabelSize * scale * globalScale)
                histBottomPanel.GetXaxis().SetTitleSize(xTitleSize * scale * globalScale)
                histBottomPanel.GetYaxis().SetTitleSize(yTitleSize * scale * globalScale)
                histBottomPanel.GetYaxis().SetTitleOffset(0.3/0.7)

                histBottomPanel.GetYaxis().SetTitleOffset(histBottomPanel.GetYaxis().GetTitleOffset() / globalScale)
            
                histBottomPanel.GetYaxis().SetRangeUser(0.55,1.45) 
                histBottomPanel.GetYaxis().SetNdivisions(5)
                histBottomPanel.GetYaxis().SetTitle("Ratio [TTvar / TT]")
                relabelNjetsBins(histBottomPanel, label.partition("_")[0], 1.1 * 0.18)
                histBottomPanel.Draw("hist E")

                line.Draw("same")

                # NonOptimized
                if args.div == "NonOptimized" and args.sig == "RPV":
                    canvas_each.SaveAs("%s/%s_plots_MCcorrectionFactorRatio/"%(args.outpath, year) + year + "_" + args.sig + "_" + args.mass + tag + ttvar + "_" + label + ".pdf")

                if args.div == "NonOptimized" and args.sig == "SYY":
                    canvas_each.SaveAs("%s/%s_plots_MCcorrectionFactorRatio/"%(args.outpath, year) + year + "_" + args.sig + "_" + args.mass + tag + ttvar + "_" + label + ".pdf")
        
 
                # MassExclusion
                if args.div == "MassExclusion" and args.sig == "RPV":
                    canvas_each.SaveAs("%s/%s_plots_MCcorrectionFactorRatio/"%(args.outpath, year) + year + "_" + args.sig + "_" + args.mass + tag + ttvar + "_" + label + ".pdf")
 
                if args.div == "MassExclusion" and args.sig == "SYY":
                    canvas_each.SaveAs("%s/%s_plots_MCcorrectionFactorRatio/"%(args.outpath, year) + year + "_" + args.sig + "_" + args.mass + tag + ttvar + "_" + label + ".pdf")

                
                # MaxSign
                if args.div == "MaxSign" and args.sig == "RPV":
                    canvas_each.SaveAs("%s/%s_plots_MCcorrectionFactorRatio/"%(args.outpath, year) + year + "_" + args.sig + "_" + args.mass + tag + ttvar + "_" + label + ".pdf")

                if args.div == "MaxSign" and args.sig == "SYY":
                    canvas_each.SaveAs("%s/%s_plots_MCcorrectionFactorRatio/"%(args.outpath, year) + year + "_" + args.sig + "_" + args.mass + tag + ttvar + "_" + label + ".pdf")


            # ---------------------------------------------
            # save canvas & legend including all histograms
            # ---------------------------------------------
            canvas.cd()
            addCMSlogo(canvas, year, TopMargin=0.06, LeftMargin=0.12, RightMargin=0.03, SF=1.0)    
            addExtraInfo(canvas, 0.12, 0.95, 0.035, args.sig, label.partition("_")[0], args.div)

            legend.Draw("SAME")

            # NonOptimized
            if args.div == "NonOptimized" and args.sig == "RPV":
                canvas.SaveAs("%s/%s_plots_MCcorrectionFactorRatio/"%(args.outpath, year) + year + "_" + args.sig + "_" + args.mass + "_MCcorr_Ratio_MC_" + label + ".pdf")

            if args.div == "NonOptimized" and args.sig == "SYY":
                canvas.SaveAs("%s/%s_plots_MCcorrectionFactorRatio/"%(args.outpath, year) + year + "_" + args.sig + "_" + args.mass + "_MCcorr_Ratio_MC_" + label + ".pdf")


            # MassExclusion
            if args.div == "MassExclusion" and args.sig == "RPV":
                canvas.SaveAs("%s/%s_plots_MCcorrectionFactorRatio/"%(args.outpath, year) + year + "_" + args.sig + "_" + args.mass + "_MCcorr_Ratio_MC_" + label + "_Cat" + args.cat + ".pdf")     
        
            if args.div == "MassExclusion" and args.sig == "SYY":
                canvas.SaveAs("%s/%s_plots_MCcorrectionFactorRatio/"%(args.outpath, year) + year + "_" + args.sig + "_" + args.mass + "_MCcorr_Ratio_MC_" + label + "_Cat" + args.cat +  ".pdf")


            # MaxSign
            if args.div == "MaxSign" and args.sig == "RPV":
                canvas.SaveAs("%s/%s_plots_MCcorrectionFactorRatio/"%(args.outpath, year) + year + "_" + args.sig + "_" + args.mass + "_MCcorr_Ratio_MC_" + label + "_Cat" + args.cat + ".pdf")

            if args.div == "MaxSign" and args.sig == "SYY":
                canvas.SaveAs("%s/%s_plots_MCcorrectionFactorRatio/"%(args.outpath, year) + year + "_" + args.sig + "_" + args.mass + "_MCcorr_Ratio_MC_" + label + "_Cat" + args.cat + ".pdf")

 
        f1.Close()
 
if __name__ == "__main__":
    main()


         
