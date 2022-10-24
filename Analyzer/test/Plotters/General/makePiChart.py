#!/bin/python
import ROOT
import math
import array
import numpy as np

ROOT.gROOT.SetBatch(True)

def main() :

    #vals        = array.array('f', [81255,1626,3865,6455] )
    vals        = array.array('f', [29079,358,1797,371] )

    colors          = array.array('i', [40,38,30,41] )
    labels          = ["t#bar{t}", "t#bar{t}+X", "QCD multijet", "Other" ]

    tpie            = ROOT.TPie( "tpie", "", 4, vals, colors )

    c1              = ROOT.TCanvas( "c1", "c1", 1400, 700 )
    c1.Divide(2,1)
    c1.SetFillStyle(4000)
    c1.SetFillColorAlpha(1, 0.)

    c1.cd(1)
    ROOT.gPad.SetFillStyle(4000)
    ROOT.gPad.SetFillColorAlpha(1, 0.)
    ROOT.gPad.SetBottomMargin(0.0)
    #ROOT.gPad.SetLeftMargin(0.0)
    #ROOT.gPad.SetRightMargin(0.0)
    tpie.SetRadius( .30 )
    tpie.SetLabelsOffset( .05 )
    tpie.SetEntryLabel( 0, "t#bar{t}" )
    tpie.SetEntryLabel( 2, "QCD" )
    tpie.SetEntryLabel( 1, "t#bar{t}+X" )
    tpie.SetEntryLabel( 3, "Other" )
    tpie.SetLabelFormat( "%txt (%perc)")
    tpie.SetTextSize(.025)
    #tpie.SetX(0.1)
    #tpie.SetY(0.5)

    l1              = tpie.MakeLegend()
    l1.SetY1(0.88)
    l1.SetY2(0.95)
    l1.SetX1(0.0)
    l1.SetX2(1.05)
    l1.SetTextSize(.05)
    l1.SetBorderSize(0)
    l1.SetFillStyle(4000)
    l1.SetFillColorAlpha(1, 0.)
    l1.SetNColumns(4)
    tpie.Draw("3d ")
    l1.Draw("SAME")

    c1.Update()
    c1.SaveAs("2016McComposition_3d.pdf")
  
    #tpie.Draw("nol <")
    #l1.SetX1(.80)
    #l1.SetY2(.86)
    #l1.Draw("SAME")

    #c1.Update()
    #c1.SaveAs("Run2McComposition_2d.pdf")

    c1.cd(2)
    ROOT.gPad.SetFillStyle(4000)
    ROOT.gPad.SetFillColorAlpha(1, 0.)
    ROOT.gPad.SetTopMargin(0.01)
    ROOT.gPad.SetLeftMargin(0.0)
    ROOT.gPad.SetRightMargin(0.0)

    #crvals      = array.array('f', [16424,131565,229,2886])
    crvals      = array.array('f', [58458,377328,779,8925])
    crcolors     = array.array('i', [ROOT.kBlue-6,ROOT.kGreen+1,ROOT.kOrange+2,ROOT.kMagenta+2] )
    crtpie            = ROOT.TPie( "crtpie", "", 4, crvals, crcolors )

    #c2              = ROOT.TCanvas( "c2", "c2", 700, 700 )
    #c2.SetFillStyle(4000)
    #c2.SetFillColorAlpha(1, 0.)

    crtpie.SetRadius( .30 )
    crtpie.SetLabelsOffset( .05 )
    crtpie.SetEntryLabel( 0, "t#bar{t}" )
    crtpie.SetEntryLabel( 1, "QCD" )
    crtpie.SetEntryLabel( 2, "t#bar{t}+X" )
    crtpie.SetEntryLabel( 3, "Other" )
    crtpie.SetLabelFormat( "%txt (%perc)")
    crtpie.SetTextSize(.025)
    #crtpie.SetX(0.9)
    #crtpie.SetY(0.5)

    #l1              = crtpie.MakeLegend()
    #l1.SetY1(.66)
    #l1.SetY2(.86)
    #l1.SetTextSize(.025)
    #l1.SetBorderSize(0)
    #l1.SetFillStyle(4000)
    #l1.SetFillColorAlpha(1, 0.)
    crtpie.Draw("3d SAME")
    #l1.Draw("SAME")

    c1.Update()
    c1.SaveAs("Run2McComposition_3d.pdf")
 
    #crtpie.Draw("nol <")
    #l1.SetX1(.80)
    #l1.SetY2(.86)
    #l1.Draw("SAME")

    #c2.Update()
    #c2.SaveAs("2016McComposition_CR_2d.pdf")

if __name__ == "__main__" :
    main()
