import os 
import sys
import ROOT
import argparse


# -------------------------------------------------
# make up the histograms such as rebin, title, etc.
# -------------------------------------------------
def do_Options(histo, histoName, theMap):

    is1D = "TH1" in histo.ClassName()

    for axis, options in theMap[histoName].iteritems():

        # x-axis
        if axis == "X":

            if "rebin" in options:
                if is1D: 
                    histo.Rebin(options["rebin"])
                else: 
                    histo.RebinX(options["rebin"])

            if "min" in options and "max" in options: 
                histo.GetXaxis().SetRangeUser(options["min"],options["max"])

            if "title" in options: 
                histo.GetXaxis().SetTitle(options["title"])

        #y-axis
        if axis == "Y":

            if "rebin" in options:
                if is1D: 
                    histo.Rebin(options["rebin"])
                else: 
                    histo.RebinY(options["rebin"])

            if "min" in options and "max" in options: 
                histo.GetYaxis().SetRangeUser(options["min"],options["max"])

            if "title" in options: 
                histo.GetYaxis().SetTitle(options["title"])

        #z-axis
        if axis == "Z":

            if "min" in options and "max" in options: 
                histo.GetZaxis().SetRangeUser(options["min"],options["max"])

# ---------------------
# update all font sizes
# ---------------------
def pretty_Histo(histo,magicFactor=1.0,magicFactor2=1.0):
   
    histo.GetXaxis().SetLabelSize(magicFactor*0.055)
    histo.GetXaxis().SetTitleSize(magicFactor*0.08)
    histo.GetXaxis().SetTitleOffset(0.7/magicFactor2)
 
    histo.GetYaxis().SetLabelSize(magicFactor*0.055)
    histo.GetYaxis().SetTitleSize(magicFactor*0.08)
    histo.GetYaxis().SetTitleOffset(0.6/magicFactor)

    histo.GetZaxis().SetLabelSize(magicFactor*0.055)
    histo.GetZaxis().SetTitleSize(magicFactor*0.06)

# --------------------------------------
# get the histograms
# make a dictionary for putting all them
#   -- with/without HEM labels
#   -- with 0l/1l labels
# --------------------------------------
def fill_Map(inRootDir, theMap):

    theMap["HEM"]   = {}
    theMap["NOHEM"] = {}

    for histoFile in os.listdir(inRootDir):

        if ".root" not in histoFile: continue
        
        histoFile = ROOT.TFile.Open(inRootDir + histoFile, "READ")
        
        for hkey in histoFile.GetListOfKeys():
            if "TH" not in hkey.GetClassName(): continue

            if hkey.GetName() == "EventCounter": continue

            keyName = ""
            if "HEM" in hkey.GetName(): 
                keyName = "HEM"
            else: 
                keyName = "NOHEM"

            name  = hkey.GetName()
            name  = name.replace("_HEM","")
            histo = hkey.ReadObj()
            histo.SetDirectory(0)
            histo.Sumw2()
           
            if name in theMap[keyName].keys(): 
                theMap[keyName][name].Add(histo)
            else: 
                theMap[keyName][name] = histo



if __name__ == '__main__':

    # --------------------------
    # make sure about statistics
    # --------------------------
    ROOT.TH1.SetDefaultSumw2()
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat("")
    ROOT.gStyle.SetLineWidth(2)
    ROOT.gStyle.SetPaintTextFormat("3.2f")
    ROOT.gStyle.SetFrameLineWidth(2)
    ROOT.gStyle.SetErrorX(0)

    # ---------------------------------------------------------------------------------------------------
    # command line options:
    #   -- python HEM_Plotter.py --inputDir 2018UL_HEMstudy_0l_1l --data JetHT --ratio          // for 0l
    #   -- python HEM_Plotter.py --inputDir 2018UL_HEMstudy_0l_1l --data SingleElectron --ratio // for 1l
    #   -- python HEM_Plotter.py --inputDir 2018UL_HEMstudy_0l_1l --data SingleMuon --ratio     // for 1l
    # ---------------------------------------------------------------------------------------------------
    usage  = "usage: %prog [options]"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("--inputDir", dest="inputDir", help="input directory",                   default="",     type=str)
    parser.add_argument("--data",     dest="data",     help="JetHT, SingleElectron, SingleMuon", default="NULL", type=str)
    parser.add_argument("--ratio",    dest="ratio",    help="Draw ratio", action="store_true",   default=False           )
    
    args = parser.parse_args()

    # --------------------
    # make histograms list
    # --------------------
    optionsMap = {
            "h_njets"           : {"X" : {"rebin" : 1,                          "title" : "N_{jets}"             }},
            "h_nbjets"          : {"X" : {"rebin" : 1,                          "title" : "N_{bjets}"            }},
            "h_ht"              : {"X" : {"rebin" : 5, "min" : 0, "max" : 2500, "title" : "H_{T} [GeV]"          }},
            "h_met"             : {"X" : {"rebin" : 4, "min" : 0, "max" : 600,  "title" : "MET [GeV]"            }},
            "h_jetPt"           : {"X" : {"rebin" : 5, "min" : 0, "max" : 1000, "title" : "Jet p_{T} [GeV]"      }},
            "h_jetPtMax"        : {"X" : {"rebin" : 5, "min" : 0, "max" : 1000, "title" : "Max Jet p_{T} [GeV]"  }},
            "h_lvMET_cm_mass"   : {"X" : {"rebin" : 1, "min" : 0, "max" : 1000, "title" : "Mass [GeV]"           }},
            "h_lvMET_cm_eta"    : {"X" : {"rebin" : 1, "min" : 0, "max" : 1000, "title" : "#eta [GeV]"           }},
            "h_lvMET_cm_phi"    : {"X" : {"rebin" : 1, "min" : 0, "max" : 1000, "title" : "#phi [GeV]"           }},
            "h_lvMET_cm_pt"     : {"X" : {"rebin" : 1, "min" : 0, "max" : 1000, "title" : "p_{T} [GeV]"          }},
            "h_fwm2_top6"       : {"X" : {"rebin" : 30                                                           }},
            "h_fwm3_top6"       : {"X" : {"rebin" : 30                                                           }},
            "h_fwm4_top6"       : {"X" : {"rebin" : 30                                                           }},
            "h_fwm5_top6"       : {"X" : {"rebin" : 30                                                           }},
            "h_fwm6_top6"       : {"X" : {"rebin" : 30                                                           }},
            "h_fwm7_top6"       : {"X" : {"rebin" : 30                                                           }},
            "h_fwm8_top6"       : {"X" : {"rebin" : 30                                                           }},
            "h_fwm9_top6"       : {"X" : {"rebin" : 30                                                           }},
            "h_fwm10_top6"      : {"X" : {"rebin" : 30                                                           }},
            "h_jmt_ev0_top6"    : {"X" : {"rebin" : 30                                                           }},
            "h_jmt_ev1_top6"    : {"X" : {"rebin" : 30                                                           }},
            "h_jmt_ev2_top6"    : {"X" : {"rebin" : 30                                                           }},
            "h_lvMET_cm_m"      : {"X" : {"rebin" : 30, "title" : "Mass [GeV]"                                   }},
            "h_lvMET_cm_phi"    : {"X" : {"rebin" : 30, "title" : "#phi"                                         }},
            "h_lvMET_cm_eta"    : {"X" : {"rebin" : 30, "title" : "#eta"                                         }},
            "h_lvMET_cm_pt"     : {"X" : {"rebin" : 30, "title" : "p_{T} [GeV]"                                  }},
            "h_beta_z"          : {"X" : {"rebin" : 30                                                           }},
            "h_jet_etaphi"      : {"X" : {"rebin" : 360, "title" : "#eta"}, "Y" : {"rebin" : 80, "title" : "#phi"}},
            "h_electron_etaphi" : {"X" : {"rebin" : 360, "title" : "#eta"}, "Y" : {"rebin" : 80, "title" : "#phi"}},
            "h_muon_etaphi"     : {"X" : {"rebin" : 360, "title" : "#eta"}, "Y" : {"rebin" : 80, "title" : "#phi"}},
    }

    # -----------------------------------
    # grab things for command line option
    # -----------------------------------
    ratio   = args.ratio
    XCANVAS = 2400
    YCANVAS = 2400

    inRootDir = args.inputDir
    outpath   = "./2018UL_HEM_Study_0l_1l"
    
    if not os.path.exists(outpath): 
        os.makedirs(outpath)

    mapPFAhistos = {}
    fill_Map(inRootDir, mapPFAhistos)

    # -------------------------
    # Save the final histograms
    # -------------------------
    for name in mapPFAhistos.values()[0].keys():

        if "0l" in name and "Single" in args.data: continue 

        if "1l" in name and "JetHT" in args.data: continue

        magicMargins = {"T" : 0.02625, "B" : 0.13375, "L" : 0.11, "R" : 0.12}

        # -------------
        # make 2D plots
        # ------------- 
        if "TH2" in mapPFAhistos.values()[0][name].ClassName():
            
            # -----------------------
            # 2D plots with HEM issue
            # -----------------------
            c1 = ROOT.TCanvas("%s_HEM"%(name), "%s_HEM"%(name), XCANVAS, int(0.8*YCANVAS))
            c1.cd()

            ROOT.gPad.SetTopMargin(magicMargins["T"])
            ROOT.gPad.SetBottomMargin(magicMargins["B"])
            ROOT.gPad.SetLeftMargin(magicMargins["L"])
            ROOT.gPad.SetRightMargin(magicMargins["R"])

            theName = name.replace("_HEM","")
            data1   = mapPFAhistos["HEM"][name]
            data2   = mapPFAhistos["NOHEM"][name]

            pretty_Histo(data1,0.875,0.8)
            pretty_Histo(data2,0.875,0.8)

            do_Options(data1, theName, optionsMap)
            do_Options(data2, theName, optionsMap)

            data1.SetTitle(""); data2.SetTitle("")
            data1.SetContour(255); data2.SetContour(255)
            data1.Draw("COLZ TEXT E")
            c1.SaveAs("%s/%s_HEM.pdf"%(outpath,name))

            # --------------------------
            # 2D plots without HEM issue
            # --------------------------
            c2 = ROOT.TCanvas("%s_NOHEM"%(name), "%s_NOHEM"%(name), XCANVAS, int(0.8*YCANVAS))
            c2.cd()

            ROOT.gPad.SetTopMargin(magicMargins["T"])
            ROOT.gPad.SetBottomMargin(magicMargins["B"])
            ROOT.gPad.SetLeftMargin(magicMargins["L"])
            ROOT.gPad.SetRightMargin(magicMargins["R"])

            data2.Draw("COLZ TEXT E")
            c2.SaveAs("%s/%s_NOHEM.pdf"%(outpath,name))

            # -----------------------
            # make ratio for 2D plots
            # -----------------------
            if data1.Integral() != 0:

                c1 = ROOT.TCanvas("%s_ratio"%(name), "%s_ratio"%(name), XCANVAS, int(0.8*YCANVAS))
                c1.cd()

                ROOT.gPad.SetTopMargin(magicMargins["T"])
                ROOT.gPad.SetBottomMargin(magicMargins["B"])
                ROOT.gPad.SetLeftMargin(magicMargins["L"])
                ROOT.gPad.SetRightMargin(magicMargins["R"])

                data1.Scale(1./data1.Integral())
                data2.Scale(1./data2.Integral())
                data1.Divide(data2)
                data1.GetZaxis().SetRangeUser(0.5,2.0)
                data1.Draw("COLZ TEXT E")
                c1.SaveAs("%s/%s_ratio.pdf"%(outpath,name))

        # -------------
        # make 1D plots
        # -------------
        elif "TH1" in mapPFAhistos.values()[0][name].ClassName():

            if ratio:
                XMin = 0;    XMax = 1; RatioXMin = 0; RatioXMax = 1 
                YMin = 0.30; YMax = 1; RatioYMin = 0; RatioYMax = 0.30
                PadFactor = (YMax-YMin) / (RatioYMax-RatioYMin)

                c1 = ROOT.TCanvas("%s"%(name), "%s"%(name), XCANVAS, YCANVAS) 
                c1.Divide(1,2)
               
                # -------------------------------------------
                # superimpose the histograms with/without HEM
                # -------------------------------------------
                c1.cd(1)
                ROOT.gPad.SetLogy()
                ROOT.gPad.SetLogz()
                ROOT.gPad.SetPad(XMin, YMin, XMax, YMax)
                ROOT.gPad.SetGridy(); ROOT.gPad.SetGridx()
                ROOT.gPad.SetTopMargin(0.03)
                ROOT.gPad.SetLeftMargin(0.11)
                ROOT.gPad.SetBottomMargin(0.02)
                ROOT.gPad.SetRightMargin(0.04)

                data1 = mapPFAhistos["HEM"][name]
                data2 = mapPFAhistos["NOHEM"][name]
                ratio = ROOT.TH1F("ratio_%s"%(name), "ratio_%s"%(name), data1.GetNbinsX(), data1.GetXaxis().GetXmin(), data1.GetXaxis().GetXmax()) 

                pretty_Histo(data1)
                pretty_Histo(data2)
                pretty_Histo(ratio,PadFactor)

                theName = data1.GetName().replace("_HEM", "")
                if theName in optionsMap:

                    do_Options(data1, theName, optionsMap)
                    do_Options(data2, theName, optionsMap)
                    do_Options(ratio, theName, optionsMap)

                data2.SetMarkerColor(ROOT.kBlack)
                data2.SetLineColor(ROOT.kBlack)
                data2.SetMarkerSize(3)
                data2.SetLineWidth(2)
                data2.SetMarkerStyle(20)
                data1.SetMarkerColor(ROOT.kRed)
                data1.SetLineColor(ROOT.kRed)
                data1.SetMarkerSize(3)
                data1.SetLineWidth(2)
                data1.SetMarkerStyle(20)
                data2.SetTitle("")
                data1.SetTitle("")

                data1.Scale(1./data1.Integral())
                data2.Scale(1./data2.Integral())
                #data1.GetYaxis().SetRangeUser(0.001,1.1*data1.GetMaximum())
                data1.GetXaxis().SetLabelSize(0)
                data1.Draw("EP")
                data2.Draw("EP SAME")

                # ---------------------------------------------------
                # put the ratio plot between with HEM and without HEM
                # ---------------------------------------------------
                c1.cd(2)

                ROOT.gPad.SetGridy()
                ROOT.gPad.SetTopMargin(0.10)
                ROOT.gPad.SetBottomMargin(0.30)
                ROOT.gPad.SetRightMargin(0.04)
                ROOT.gPad.SetLeftMargin(0.11)
                ROOT.gPad.SetPad(RatioXMin, RatioYMin, RatioXMax, RatioYMax)

                ratio.SetTitle("")
                ratio.Divide(data1,data2)

                ratio.GetYaxis().SetRangeUser(0.50,1.5)
                ratio.GetYaxis().SetNdivisions(-304)
                ratio.GetYaxis().SetTitle("HEM / Nominal")
                ratio.GetYaxis().SetTitleSize(ratio.GetYaxis().GetTitleSize()*0.6)
                ratio.GetXaxis().SetTitleSize(ratio.GetXaxis().GetTitleSize()*0.8)
                ratio.GetYaxis().SetTitleOffset(0.5)
                ratio.GetXaxis().SetTitleOffset(0.9)
                ratio.SetMarkerStyle(20); ratio.SetMarkerSize(3); ratio.SetMarkerColor(ROOT.kBlue+2)
                ratio.SetLineWidth(2); ratio.SetLineColor(ROOT.kBlue+2)

                ratio.Draw("EP")
                c1.SaveAs("%s/%s.pdf"%(outpath,name))



