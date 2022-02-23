import ROOT 

# This file defines various plot routines that are general enough to be useable for many small studies. 

# Define some global colors
mycolors = [ROOT.kMagenta+1, ROOT.kOrange+2, ROOT.kBlue-7, ROOT.kGreen+1, ROOT.kCyan+1, ROOT.kYellow+1, ROOT.kRed+1]
mysignalcolors = [ROOT.kViolet-6, ROOT.kAzure-6, ROOT.kPink-6, ROOT.kTeal-6]
deflumi = 35900
defcolors = [ROOT.kRed+2, ROOT.kBlue+1, ROOT.kGreen+2, ROOT.kMagenta+1, ROOT.kCyan, ROOT.kOrange+1, ROOT.kRed, ROOT.kGreen, ROOT.kBlack, ROOT.kGray+2]

ROOT.gROOT.SetBatch(True)

# Make shape comparison plot
def makeplot(histos, names, xlabel, plotname, plotdir="./", linear=True, lumi=deflumi, legendColumns=1, colors=mycolors, norm=True, drawstyle="default", 
             ymin=0, ymax=None, ylabel="A.U.", dropzeroes=True):
    """
histos: list of histograms to plot
names: corresponding histogram name to go in the legend
xlabel: xaxis title
plotname: name of the plot
plotdir: where to save the plot
"""
    legend = ROOT.TLegend(0.6, 0.75, 0.9, 0.9)
    legend.SetNColumns(legendColumns)

    max_y = 0
    print len(histos)
    for i,h in enumerate(histos):
        if dropzeroes and h.Integral() == 0: continue
        h.SetLineColor(colors[i])
        h.SetLineWidth(2)
        #h.SetFillColorAlpha(colors[i],0.3)
        #if names[i] == "TT" or names[i] == "ttbar":
        #    h.SetLineWidth(4)
        h.GetXaxis().SetTitle(xlabel)
        h.GetXaxis().SetTitleSize(0.055)
        h.GetXaxis().SetTitleOffset(0.85)
        h.GetYaxis().SetTitle(ylabel)
        h.GetYaxis().SetTitleSize(0.055)
        h.GetYaxis().SetTitleOffset(0.85)
        h.SetTitle("")
        if norm:
            h.Scale(1./h.Integral())
        legend.AddEntry(h, names[i], "l")
        if h.GetMaximum() > max_y:
            max_y = h.GetMaximum()

    canvas = ROOT.TCanvas("c_"+plotname,"c_"+plotname, 800, 800)
    firsthisto = None
    if drawstyle == "default":
        firsthisto = histos[0].DrawCopy("hist")
        histos[0].SetFillColorAlpha(colors[0],0.3)
        histos[0].Draw("E2 same")
        for i,h in enumerate(histos[1:]):
            #h.Draw("histsame")
            h.DrawCopy("hist same")
            h.SetFillColorAlpha(colors[i+1],0.3)
            h.Draw("E2 same")
    elif drawstyle == "comp":
        firsthisto = histos[0].DrawCopy("hist")
        histos[0].SetFillColorAlpha(colors[0],0.3)
        histos[0].Draw("E2 same")
        for i,h in enumerate(histos[1:]):
            h.Draw("histsame")
            #h.DrawCopy("hist same")
            #h.SetFillColorAlpha(colors[i+1],0.3)
            #h.Draw("E2 same")
    elif drawstyle == "lastP":
        firsthisto = histos[0]
        histos[0].Draw("hist")
        for i,h in enumerate(histos[1:-1]):
            h.Draw("histsame")
        histos[-1].SetMarkerStyle(4)
        histos[-1].SetMarkerColor(colors[len(histos)-1])
        histos[-1].Draw("PEX0same")
    else:
        firsthisto = histos[0]
        histos[0].Draw(drawstyle)
        for i,h in enumerate(histos[1:]):
            h.Draw(drawstyle+"same")       
    legend.Draw()
    if not linear:
        canvas.SetLogy()
        histos[0].GetYaxis().SetRangeUser(ymin if ymin else 0.000001, ymax if ymax else max_y*3)
    else:
        firsthisto.GetYaxis().SetRangeUser(ymin if ymin else 0, ymax if ymax else max_y*1.4)
        #histos[0].GetYaxis().SetRangeUser(ymin if ymin else 0, ymax if ymax else max_y*1.4)
    canvas.SaveAs(plotdir+plotname+".pdf")
    canvas.SaveAs(plotdir+plotname+".png")


# Make a stack plot with background, signal, and data
def makestackplot(histos_bg_, names_bg, histos_signal_, names_signal, histo_data_, name_data, xlabel, plotname, plotdir="./", linear=False, myrange=None, lumi=deflumi, legendColumns=1, colors=mycolors, signalcolors=mysignalcolors):
    """
histos_bg: histograms for the backgrounds that should be stacked
names_bg: names of the backgrounds that will go in the legend
histos_signal: histograms for the signals that should be overlaid
names_signal: names of the signals that will go in the legend
histo_data: histogram for the data that should be drawn
name_data: names of the data that will go in the legend
xlabel: name of the xaxis title
plotname: Name of the plot (and the stack)
    """
    histos_bg = [h.Clone() for h in histos_bg_]
    if histos_signal_:
        histos_signal = [h.Clone() for h in histos_signal_]
    else:
        histos_signal = None
    if histo_data_:
        histo_data = histo_data_.Clone()
    else:
        histo_data = None
    legend = None
    if legendColumns >= 2:
        legend = ROOT.TLegend(0.35, 0.75, 0.95, 0.85)
        legend.SetNColumns(3)
    else:
        legend = ROOT.TLegend(0.55, 0.7, 0.95, 0.85)
        legend.SetNColumns(legendColumns)
    legend.SetBorderSize(0)

    stack = ROOT.THStack(plotname,"")

    canvas = ROOT.TCanvas("c_"+plotname, "c_"+plotname, 800, 800)
    canvas.SetRightMargin(0.03)
    
    if histo_data:
        canvas.Divide(1, 2)    
        canvas.cd(1)
        ROOT.gPad.SetPad("p1", "p1", 0, 2.5 / 9.0, 1, 1, ROOT.kWhite, 0, 0)
        ROOT.gPad.SetBottomMargin(0.01)
        ROOT.gPad.SetLeftMargin(0.12)
        ROOT.gPad.SetRightMargin(0.03)
        ROOT.gPad.SetTopMargin(0.06 * (8.0 / 6.5))
    if not linear:
        ROOT.gPad.SetLogy()
        
    htotal_bg = None
    for i,h in enumerate(histos_bg):
        #h.Scale(lumi)
        h.SetLineColor(colors[i])
        h.SetMarkerColor(colors[i])
        h.SetFillColor(colors[i])
        stack.Add(h)
        legend.AddEntry(h, names_bg[i], "f")
        if htotal_bg is None:
            htotal_bg = h.Clone("htotal_bg")
        else:
            htotal_bg.Add(h)

    stack.Draw("hist")
    if not histo_data:
        stack.GetXaxis().SetTitle(xlabel)
    stack.GetYaxis().SetTitle("Events")
    #stack.GetYaxis().SetLabelSize(0.13)
    stack.GetYaxis().SetTitleSize(0.05)
    stack.GetXaxis().SetTitleSize(0.05)
    stack.GetYaxis().SetTitleOffset(0.9)
    stack.GetXaxis().SetTitleOffset(0.85)
    if histo_data:
        stack.GetYaxis().SetLabelSize(0.13*2.5/6.5)
        stack.GetYaxis().SetTitleSize(0.15*2.5/6.5)
        stack.GetYaxis().SetTitleOffset(0.95)
    stack.SetMinimum(1)
    if linear:
        stack.SetMaximum(1.3*stack.GetMaximum())
    else:
        stack.SetMaximum(50*stack.GetMaximum())

    if histos_signal:
        for i,h in enumerate(histos_signal):
            #h.Scale(lumi)
            h.SetLineColor(signalcolors[i])
            h.SetLineStyle(7)
            h.SetLineWidth(4)
            legend.AddEntry(h, names_signal[i], "l")
            h.Draw("histsame")

    if histo_data:
        histo_data.SetMarkerColor(ROOT.kBlack)
        histo_data.SetMarkerStyle(20)
        histo_data.SetMarkerSize(1.2)
        legend.AddEntry(histo_data, name_data, "ep")
        histo_data.Draw("PEX0same")

    if myrange is not None:
        stack.GetXaxis().SetRangeUser(myrange[0], myrange[1])

    legend.Draw()

    if histo_data:
        canvas.cd(2)
        ROOT.gPad.SetPad("p2", "p2", 0, 0, 1, 2.5 / 9.0, ROOT.kWhite, 0, 0)
        ROOT.gPad.SetLeftMargin(0.12)
        ROOT.gPad.SetRightMargin(0.03)
        ROOT.gPad.SetTopMargin(0.01)
        ROOT.gPad.SetBottomMargin(0.37)

        # Make ratio of Data over total stack
        myratio = histo_data.Clone("ratio_"+histo_data.GetName())
        myratio.SetTitle("")
        myratio.Divide(htotal_bg)
        myratio.GetXaxis().SetTitle(xlabel)
        myratio.GetYaxis().SetTitle("Data/MC")
        myratio.GetXaxis().SetLabelSize(0.15)
        myratio.GetXaxis().SetTitleSize(0.15)
        myratio.GetYaxis().SetLabelSize(0.13)
        myratio.GetYaxis().SetTitleSize(0.15)
        myratio.GetYaxis().SetTitleOffset(0.37)
        myratio.GetYaxis().SetRangeUser(0, 2.2)
        myratio.GetYaxis().SetNdivisions(4, 2, 0)
        myratio.Draw()
        canvas.Update()
        line = ROOT.TF1("line", "1", myratio.GetBinLowEdge(1), myratio.GetBinLowEdge(myratio.GetNbinsX()) + myratio.GetBinWidth(myratio.GetNbinsX()))
        line.SetLineColor(ROOT.kGray+2)
        line.SetLineStyle(ROOT.kDashed)
        line.Draw("same")

    if linear:
        canvas.SaveAs(plotdir+plotname+"_linear.pdf")
    else:
        canvas.SaveAs(plotdir+plotname+".pdf")
    




def makeplot2d(filename, histoname, variableX, variableY, label, plotname, dotext=True):

    files = ROOT.TFile.Open(filename)

    histo = files.Get(histoname)
    histo.GetXaxis().SetTitle(variableX)
    histo.GetYaxis().SetTitle(variableY)
    #histo.GetZaxis().SetTitle("Fraction of events")
    histo.SetTitle("")
    histo.SetMarkerSize(2)
    total = histo.GetEntries()
    #histo.Scale(1./total)

    canvas = ROOT.TCanvas("c","c")
    canvas.SetRightMargin(0.15)
    if dotext:
        histo.Draw("colztext")
        histo.GetZaxis().SetRangeUser(0.,1.)
    else:
        histo.Draw("colz")

    # also add profile
    prof = histo.ProfileX()
    prof.SetLineWidth(3)
    prof.SetLineColor(ROOT.kBlack)
    prof.Draw("same")

    text = ROOT.TText(0.1,0.95,label)
    text.SetNDC()
    text.SetTextColor(ROOT.kRed+2)
    text.Draw()

    canvas.SaveAs(plotdir+plotname+".pdf")

    closed = files.Close()

def makeplotEff2d(filename, histoname, label, plotname):

    files = ROOT.TFile.Open(filename)
    canvas = ROOT.TCanvas("c","c")
    canvas.SetRightMargin(0.15)
    histo = files.Get(histoname)
    ROOT.gStyle.SetPaintTextFormat(".3f")
    histo.SetMarkerSize(2)
    histo.Draw("colztext")
    ROOT.gPad.Update()
    histo.GetPaintedHistogram().GetXaxis().SetNdivisions(3)
    histo.GetPaintedHistogram().GetYaxis().SetNdivisions(3)
    histo.GetPaintedHistogram().GetZaxis().SetRangeUser(0.,1.)
    histo.GetPaintedHistogram().GetZaxis().SetTitle("fraction")

    text = ROOT.TText(0.1,0.95,label)
    text.SetNDC()
    text.SetTextColor(ROOT.kRed+2)
    text.Draw()

    canvas.SaveAs(plotdir+plotname+".pdf")

    closed = files.Close()

def makeplotEff2dQCD(filenames_qcd, histoname, label, plotname):
    files = [ROOT.TFile.Open(filename) for filename in filenames_qcd]
    histos = [f.Get(histoname) for f in files]

    for i,histo in enumerate(histos):
        label = filenames_qcd[i].rsplit("/")[-1].replace(".root","")
        histo.SetWeight(weights[label])
    total = histos[0].Clone()
    for h in histos[1:]:
        total.Add(h)

    canvas = ROOT.TCanvas("c","c")
    canvas.SetRightMargin(0.15)
    ROOT.gStyle.SetPaintTextFormat(".3f")

    total.SetMarkerSize(2)
    total.Draw("colztext")
    ROOT.gPad.Update()
    total.GetPaintedHistogram().GetXaxis().SetNdivisions(3)
    total.GetPaintedHistogram().GetYaxis().SetNdivisions(3)
    total.GetPaintedHistogram().GetZaxis().SetRangeUser(0.,1.)
    total.GetPaintedHistogram().GetZaxis().SetTitle("fraction")

    text = ROOT.TText(0.1,0.95,label)
    text.SetNDC()
    text.SetTextColor(ROOT.kRed+2)
    text.Draw()

    canvas.SaveAs(plotdir+plotname+".pdf")

    closed = [f.Close() for f in files]

def makeplot_multivar(filename, names, variables, plotname, normalized=False):

    files = ROOT.TFile.Open(filename)

    legend = ROOT.TLegend(0.6, 0.7, 0.9, 0.9)

    histos = [files.Get(variable) for variable in variables]
    colors = [ROOT.kRed+2, ROOT.kGreen+2, ROOT.kBlue+1, ROOT.kMagenta+1, ROOT.kGray+2]
    max_y = 0
    for i,h in enumerate(histos):
        if h.Integral() == 0: continue
        h.SetLineColor(colors[i])
        h.SetLineWidth(2)
        h.GetXaxis().SetTitle(plotname)
        if normalized:
            h.GetYaxis().SetTitle("A.U.")
        else:
            h.GetYaxis().SetTitle("Events")
        h.SetTitle("")
        if normalized:
            h.Scale(1./h.Integral())
        legend.AddEntry(h, names[i], "l")
        if h.GetMaximum() > max_y:
            max_y = h.GetMaximum()

    canvas = ROOT.TCanvas("c","c")
    histos[0].Draw("hist")
    for h in histos[1:]:
        h.Draw("histsame")
    legend.Draw()
    histos[0].GetYaxis().SetRangeUser(0,max_y*1.4)
    canvas.SaveAs(plotdir+plotname+".pdf")

    closed = files.Close()

from math import sqrt
def calcUnc(a,unc_a,b,unc_b):
    return  a/b * sqrt( (unc_a/a)**2+(unc_b/b)**2)

def makeNJplot(inputfile, labels, jettype, cuts, plotname):
    basename = "h_njets_" + jettype
    histos = [inputfile.Get(basename+"_"+cut) for cut in cuts]
    # scale histos
    for i,h in enumerate(histos):
        h.Scale(35900.)
        #h.SetName(h.GetName()+"_"+labels[i])
        #h.SetTitle(h.GetTitle()+"_"+labels[i])
        #h.Write()

    legend = ROOT.TLegend(0.7, 0.7, 0.98, 0.9)

    rj_histos = []
    for i,h in enumerate(histos):
        hist_rj = ROOT.TH1D(h.GetName().replace("h_njets","rj"), h.GetName().replace("h_njets","rj"),
                            h.GetNbinsX()-1, h.GetBinLowEdge(1), h.GetBinLowEdge(h.GetNbinsX()))

        for nj in range(h.GetNbinsX()-1):
            Njplus1 = h.GetBinContent(nj+2)
            Njplus1_unc = h.GetBinError(nj+2)
            Nj = h.GetBinContent(nj+1)
            Nj_unc = h.GetBinError(nj+1)
            rj = 0
            unc_rj = 0
            if Nj>0 and Njplus1>0:
                rj = Njplus1/Nj
                unc_rj = calcUnc(Njplus1, Njplus1_unc, Nj, Nj_unc)

            hist_rj.GetXaxis().SetBinLabel(nj+1, "%s/%s"%(nj+1+6, nj+6))
            hist_rj.SetBinContent(nj+1, rj)
            hist_rj.SetBinError(nj+1, unc_rj)

        legend.AddEntry(hist_rj, labels[i], "l")
        rj_histos.append(hist_rj)

    colors = [ROOT.kRed+2, ROOT.kGreen+2, ROOT.kBlue+1, ROOT.kMagenta+1, ROOT.kGray+2]
    for i,h in enumerate(rj_histos):
        h.SetLineColor(colors[i])
        h.SetLineWidth(2)
        h.GetXaxis().SetTitle("")
        h.GetXaxis().SetLabelSize(0.06)
        h.GetYaxis().SetTitle("R(j)")
        h.GetYaxis().SetTitleSize(0.06)
        h.GetYaxis().SetTitleOffset(1.15)
        h.SetTitle("")

    canvas = ROOT.TCanvas("c","c", 800,800)
    canvas.SetLeftMargin(0.15)
    canvas.SetRightMargin(0.02)
    rj_histos[0].Draw("")
    for h in rj_histos[1:]:
        h.Draw("same")
    legend.Draw()
    rj_histos[0].GetYaxis().SetRangeUser(0,0.5)
    canvas.SaveAs(plotdir+plotname+".pdf")


