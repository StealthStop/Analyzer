import os 
import sys
import ROOT
import argparse

# Sequence to run
# python HEM_Plotter.py --inputDir ./condor/HEMstudy1L_hadd --data SingleElectron --tag bothVeto --ratio
# python HEM_Plotter.py --inputDir ./condor/HEMstudy1L_hadd --data SingleMuon     --tag bothVeto --ratio
# python HEM_Plotter.py --inputDir ./condor/HEMstudy0L_hadd --data JetHT          --tag bothVeto --ratio
# 
# python HEM_Plotter.py --inputDir ./condor/HEMstudy1L_hadd --data SingleElectron --tag ePostVeto --ratio
# python HEM_Plotter.py --inputDir ./condor/HEMstudy1L_hadd --data SingleMuon     --tag ePostVeto --ratio
# python HEM_Plotter.py --inputDir ./condor/HEMstudy0L_hadd --data JetHT          --tag ePostVeto --ratio
# 
# python HEM_Plotter.py --inputDir ./condor/HEMstudy1L_hadd --data SingleElectron --tag bPostVeto --ratio
# python HEM_Plotter.py --inputDir ./condor/HEMstudy1L_hadd --data SingleMuon     --tag bPostVeto --ratio
# python HEM_Plotter.py --inputDir ./condor/HEMstudy0L_hadd --data JetHT          --tag bPostVeto --ratio
# 
# python HEM_Plotter.py --inputDir ./condor/HEMstudy1L_hadd --data SingleElectron --tag ebPostVeto --ratio
# python HEM_Plotter.py --inputDir ./condor/HEMstudy1L_hadd --data SingleMuon     --tag ebPostVeto --ratio
# python HEM_Plotter.py --inputDir ./condor/HEMstudy0L_hadd --data JetHT          --tag ebPostVeto --ratio

# -------------------------------------------------
# make up the histograms such as rebin, title, etc.
# -------------------------------------------------
def do_Options(histo, histoName, theMap):

    is1D   = "TH1" in histo.ClassName()
    doLogY = False

    optDict = None 
    for key, value in theMap.items():
        if key in histoName:
            optDict = value
            break

    if optDict == None:
        return False

    for axis, options in optDict.items():

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

        # y-axis
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

            if "logY" in options:
                doLogY = options["logY"]

        # z-axis
        if axis == "Z":

            if "min" in options and "max" in options: 
                histo.GetZaxis().SetRangeUser(options["min"],options["max"])       

    return doLogY

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
        
        histoFile = ROOT.TFile.Open(inRootDir + "/" + histoFile, "READ")
        
        for hkey in histoFile.GetListOfKeys():
            if "TH" not in hkey.GetClassName(): continue

            if hkey.GetName() == "EventCounter": continue

            keyName = ""
            if "HEM" in hkey.GetName(): 
                keyName = "HEM"
            elif "HEM" not in hkey.GetName(): 
                keyName = "NOHEM"
            else:
                continue

            name  = hkey.GetName()
            name  = name.replace("_HEM","").replace("veto", "")

            histo = hkey.ReadObj()
            histo.SetDirectory(0)
            histo.Sumw2()
           
            if name in theMap[keyName].keys(): 
                theMap[keyName][name].Add(histo)
            else: 
                theMap[keyName][name] = histo

# ------------------------------
# add the CMS logo to the canvas
# ------------------------------
def add_CMSlogo(canvas):

    TopMargin    = 0.06
    BottomMargin = 0.35
    RightMargin  = 0.05
    LeftMargin   = 0.10

    canvas.cd()
    
    mark = ROOT.TLatex()
    mark.SetNDC(True)

    mark.SetTextAlign(11)
    mark.SetTextFont(61)
    mark.SetTextSize(0.050)
    mark.DrawLatex(LeftMargin,        1 - (TopMargin - 0.015), "CMS"             )
    
    mark.SetTextFont(52)
    mark.SetTextSize(0.032)
    mark.DrawLatex(LeftMargin + 0.12, 1 - (TopMargin - 0.017), "Work in Progress")

    mark.SetTextAlign(31)
    mark.SetTextFont(42)
    mark.DrawLatex(1 - RightMargin,   1 - (TopMargin - 0.017), "2018 (13 TeV)"   )


if __name__ == '__main__':

    # --------------------------
    # make sure about statistics
    # --------------------------
    ROOT.TH1.SetDefaultSumw2()
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat("")
    ROOT.gStyle.SetLineWidth(2)
    ROOT.gStyle.SetPaintTextFormat("3.4f")
    ROOT.gStyle.SetFrameLineWidth(2)
    ROOT.gStyle.SetErrorX(0)

    # ------------------------------------------------------------------
    # command line options:
    #   -- python HEM_Plotter.py --data JetHT --ratio          // for 0l
    #   -- python HEM_Plotter.py --data SingleElectron --ratio // for 1l
    #   -- python HEM_Plotter.py --data SingleMuon --ratio     // for 1l
    # ------------------------------------------------------------------
    usage  = "usage: %prog [options]"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("--inputDir", dest="inputDir", help="input directory",                   default="condor/2_2018UL_HEM_study_2022/hadd_2018_HEM_issue_0l_1l_18.04.2022/", type=str)
    parser.add_argument("--data",     dest="data",     help="JetHT, SingleElectron, SingleMuon", default="NULL", type=str)
    parser.add_argument("--tag",      dest="tag",      help="Tag to put in folder name",         default="NULL", type=str)
    parser.add_argument("--ratio",    dest="ratio",    help="Draw ratio", action="store_true",   default=False           )
    
    args = parser.parse_args()

    # --------------------
    # make histograms list
    # --------------------
    optionsMap = {
            "h_njets"             : {  "Y" : {"logY" : True, },  "X" : {"rebin" : 1,  "min" : 6, "max" : 16,   "title" : "N_{jets}"                    }  },
            "h_nbjets"            : {  "Y" : {"logY" : True, },  "X" : {"rebin" : 1,  "min" : 0, "max" : 8,    "title" : "N_{bjets}"                   }  },
            "h_ht"                : {  "Y" : {"logY" : True, },  "X" : {"rebin" : 5,  "min" : 0, "max" : 2500, "title" : "H_{T} [GeV]"                 }  },
            "h_met"               : {  "Y" : {"logY" : True, },  "X" : {"rebin" : 4,  "min" : 0, "max" : 600,  "title" : "MET [GeV]"                   }  },
            "h_jetPt"             : {  "Y" : {"logY" : True, },  "X" : {"rebin" : 5,  "min" : 0, "max" : 1000, "title" : "Jet p_{T} [GeV]"             }  },
            "h_jetPtMax"          : {  "Y" : {"logY" : True, },  "X" : {"rebin" : 5,  "min" : 0, "max" : 1000, "title" : "Max Jet p_{T} [GeV]"         }  },
            "h_lvMET_cm_mass"     : {  "Y" : {"logY" : True, },  "X" : {"rebin" : 20, "min" : 0, "max" : 1000, "title" : "Mass [GeV]"                  }  },
            "h_lvMET_cm_eta"      : {  "Y" : {"logY" : True, },  "X" : {"rebin" : 20, "min" : 0, "max" : 1000, "title" : "#eta"                        }  },
            "h_lvMET_cm_phi"      : {  "Y" : {"logY" : True, },  "X" : {"rebin" : 20, "min" : 0, "max" : 1000, "title" : "#phi"                        }  },
            "h_lvMET_cm_pt"       : {  "Y" : {"logY" : True, },  "X" : {"rebin" : 20, "min" : 0, "max" : 1000, "title" : "p_{T} [GeV/c]"               }  },
            "h_fwm2_top6"         : {  "Y" : {"logY" : True, },  "X" : {"rebin" : 30,  "title" : "FWM 2"                                               }  },
            "h_fwm3_top6"         : {  "Y" : {"logY" : True, },  "X" : {"rebin" : 30,  "title" : "FWM 3"                                               }  },
            "h_fwm4_top6"         : {  "Y" : {"logY" : True, },  "X" : {"rebin" : 30,  "title" : "FWM 4"                                               }  },
            "h_fwm5_top6"         : {  "Y" : {"logY" : True, },  "X" : {"rebin" : 30,  "title" : "FWM 5"                                               }  },
            "h_fwm6_top6"         : {  "Y" : {"logY" : True, },  "X" : {"rebin" : 30,  "title" : "FWM 6"                                               }  },
            "h_fwm7_top6"         : {  "Y" : {"logY" : True, },  "X" : {"rebin" : 30,  "title" : "FWM 7"                                               }  },
            "h_fwm8_top6"         : {  "Y" : {"logY" : True, },  "X" : {"rebin" : 30,  "title" : "FWM 8"                                               }  },
            "h_fwm9_top6"         : {  "Y" : {"logY" : True, },  "X" : {"rebin" : 30,  "title" : "FWM 9"                                               }  },
            "h_fwm10_top6"        : {  "Y" : {"logY" : True, },  "X" : {"rebin" : 30,  "title" : "FWM 10"                                              }  },
            "h_jmt_ev0_top6"      : {  "Y" : {"logY" : True, },  "X" : {"rebin" : 30,  "title" : "JMT 0"                                               }  },
            "h_jmt_ev1_top6"      : {  "Y" : {"logY" : True, },  "X" : {"rebin" : 30,  "title" : "JMT 1"                                               }  },
            "h_jmt_ev2_top6"      : {  "Y" : {"logY" : True, },  "X" : {"rebin" : 30,  "title" : "JMT 2"                                               }  },
            "h_event_beta_z"      : {  "Y" : {"logY" : False,},  "X" : {"rebin" : 50,  "title" : "#beta_{z}"                                           }  },
            "h_jet_EtaVsPhi"      : {  "Y" : {"logY" : True, "rebin" : 12, "title" : "Jet #phi"},      "X" : {"rebin" : 12, "title" : "Jet #eta"       }  },
            "h_bjet_EtaVsPhi"     : {  "Y" : {"logY" : True, "rebin" : 12, "title" : "b Jet #phi"},    "X" : {"rebin" : 12, "title" : "b Jet #eta"     }  },
            "h_electron_EtaVsPhi" : {  "Y" : {"logY" : True, "rebin" : 12, "title" : "Electron #phi"}, "X" : {"rebin" : 12, "title" : "Electron #eta"  }  },
            "h_muon_EtaVsPhi"     : {  "Y" : {"logY" : True, "rebin" : 12, "title" : "Muon #phi"},     "X" : {"rebin" : 12, "title" : "Muon #eta"      }  },
    }

    # -----------------------------------
    # grab things for command line option
    # -----------------------------------
    ratio   = args.ratio
    XCANVAS = 2400
    YCANVAS = 2400

    inRootDir = args.inputDir
    outpath   = "./2018UL_HEM_Study_0l_1l_%s/%s/"%(args.tag,args.data)
    
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

        if ("h_electron_EtaVsPhi" in name or "h_muon_EtaVsPhi" in name) and "JetHT" in args.data: continue

        if "h_electron_EtaVsPhi" in name and ("SingleMuon" in args.data or "JetHT" in args.data): continue

        if "h_muon_EtaVsPhi" in name and ("SingleElectron" in args.data or "JetHT" in args.data): continue

        magicMargins_1D = {"T" : 0.085, "B" : 0.04, "L" : 0.11, "R" : 0.04}
        magicMargins_2D = {"T" : 0.08, "B" : 0.14, "L" : 0.11, "R" : 0.17}

        # -------------
        # make 2D plots
        # ------------- 
        if "TH2" in mapPFAhistos.values()[0][name].ClassName():

            maxRange = 0.025 
            option = "COLZ E TEXT"
            
            # -----------------------
            # 2D plots with HEM issue
            # -----------------------
            c1 = ROOT.TCanvas("%s_HEM"%(name), "%s_HEM"%(name), XCANVAS, YCANVAS)
            c1.cd()

            ROOT.gPad.SetTopMargin(magicMargins_2D["T"])
            ROOT.gPad.SetBottomMargin(magicMargins_2D["B"])
            ROOT.gPad.SetLeftMargin(magicMargins_2D["L"])
            ROOT.gPad.SetRightMargin(magicMargins_2D["R"])

            theName = name.replace("_HEM","")
            data1   = mapPFAhistos["HEM"][name]
            data2   = mapPFAhistos["NOHEM"][name]

            pretty_Histo(data1,0.875,0.8)
            pretty_Histo(data2,0.875,0.8)

            theName = data1.GetName()
            do_Options(data1, theName, optionsMap)
            do_Options(data2, theName, optionsMap)

            data1.SetTitle("")
            data2.SetTitle("")
            data1.SetContour(255)
            data2.SetContour(255)

            if data1.Integral() != 0.0:
                data1.Scale(1./data1.Integral())
            data1.Draw(option)
            data1.GetZaxis().SetRangeUser(0.0,maxRange)
            data1.GetZaxis().SetMaxDigits(2)

            data1.SetTitle("post HEM")
            #add_CMSlogo(c1)
            c1.SaveAs("%s/%s_HEM.pdf"%(outpath,name))

            # --------------------------
            # 2D plots without HEM issue
            # --------------------------
            c2 = ROOT.TCanvas("%s_NOHEM"%(name), "%s_NOHEM"%(name), XCANVAS, YCANVAS)
            c2.cd()

            ROOT.gPad.SetTopMargin(magicMargins_2D["T"])
            ROOT.gPad.SetBottomMargin(magicMargins_2D["B"])
            ROOT.gPad.SetLeftMargin(magicMargins_2D["L"])
            ROOT.gPad.SetRightMargin(magicMargins_2D["R"])

            if data2.Integral() != 0.0:
                data2.Scale(1./data2.Integral())
            data2.GetZaxis().SetRangeUser(0.0,maxRange)
            data2.GetZaxis().SetMaxDigits(2)

            data2.Draw(option)
            data2.SetTitle("pre HEM")
            #add_CMSlogo(c2)
            c2.SaveAs("%s/%s_NOHEM.pdf"%(outpath,name))

            # --------------------------
            # 2D plots combined 
            # --------------------------
            c3 = ROOT.TCanvas("%s_ALLHEM"%(name), "%s_ALLHEM"%(name), XCANVAS, YCANVAS)
            c3.cd()

            ROOT.gPad.SetTopMargin(magicMargins_2D["T"])
            ROOT.gPad.SetBottomMargin(magicMargins_2D["B"])
            ROOT.gPad.SetLeftMargin(magicMargins_2D["L"])
            ROOT.gPad.SetRightMargin(magicMargins_2D["R"])

            data2.Add(data1)
            if data2.Integral() != 0.0:
                data2.Scale(1./data2.Integral())
            data2.GetZaxis().SetRangeUser(0.0,maxRange)
            data2.GetZaxis().SetMaxDigits(2)

            data2.Draw(option)
            data2.SetTitle("all HEM")
            #add_CMSlogo(c2)
            c3.SaveAs("%s/%s_ALLHEM.pdf"%(outpath,name))

            # ---------------------------
            # make HEM ratio for 2D plots
            # ---------------------------
            if data2.Integral() != 0:

                c1 = ROOT.TCanvas("%s_ratio"%(name), "%s_ratio"%(name), XCANVAS, YCANVAS)
                c1.cd()

                ROOT.gPad.SetTopMargin(magicMargins_2D["T"])
                ROOT.gPad.SetBottomMargin(magicMargins_2D["B"])
                ROOT.gPad.SetLeftMargin(magicMargins_2D["L"])
                ROOT.gPad.SetRightMargin(magicMargins_2D["R"])

                if data1.Integral() > 0.0:
                    data1.Scale(1./data1.Integral())
                if data2.Integral() > 0.0:
                    data2.Scale(1./data2.Integral())

                data1.Divide(data2)
                data1.GetZaxis().SetRangeUser(0.5,2.0)
                data1.Draw(option)
                data1.SetTitle("HEM Ratio")
                c1.SaveAs("%s/%s_ratio.pdf"%(outpath,name))

        # -------------
        # make 1D plots
        # -------------
        elif "TH1" in mapPFAhistos.values()[0][name].ClassName():

            if ratio:
                XMin = 0;    XMax = 1; RatioXMin = 0; RatioXMax = 1 
                YMin = 0.30; YMax = 1; RatioYMin = 0; RatioYMax = 0.30
                PadFactor = (YMax-YMin) / (RatioYMax-RatioYMin)

                c1     = ROOT.TCanvas("%s"%(name), "%s"%(name), XCANVAS, YCANVAS) 
                legend = ROOT.TLegend(0.76, 0.8, 0.99, 0.9)
                c1.Divide(1,2)
               
                # -------------------------------------------
                # superimpose the histograms with/without HEM
                # -------------------------------------------
                c1.cd(1)
                ROOT.gPad.SetLogy()
                ROOT.gPad.SetLogz()
                ROOT.gPad.SetPad(XMin, YMin, XMax, YMax)
                ROOT.gPad.SetGridy()
                ROOT.gPad.SetGridx()
                ROOT.gPad.SetTopMargin(magicMargins_1D["T"])
                ROOT.gPad.SetBottomMargin(magicMargins_1D["B"])
                ROOT.gPad.SetLeftMargin(magicMargins_1D["L"])
                ROOT.gPad.SetRightMargin(magicMargins_1D["R"])

                data1 = mapPFAhistos["HEM"][name]
                data2 = mapPFAhistos["NOHEM"][name]
                ratio = ROOT.TH1F("ratio_%s"%(name), "ratio_%s"%(name), data1.GetNbinsX(), data1.GetXaxis().GetXmin(), data1.GetXaxis().GetXmax()) 

                pretty_Histo(data1)
                pretty_Histo(data2)
                pretty_Histo(ratio,PadFactor)

                theName = data1.GetName()
                do_Options(data1, theName, optionsMap)
                do_Options(data2, theName, optionsMap)
                do_Options(ratio, theName, optionsMap)

                data2.SetMarkerColor(38)
                data2.SetLineColor(38)
                data2.SetMarkerSize(7) # 5 
                data2.SetLineWidth(3)
                data2.SetMarkerStyle(29)
                data2.SetTitle("")
                data1.SetMarkerColor(46)
                data1.SetLineColor(46)
                data1.SetMarkerSize(7) # 5
                data1.SetLineWidth(3)
                data1.SetMarkerStyle(29)
                data1.SetTitle("")

                if data1.Integral() != 0.0:
                    data1.Scale(1./data1.Integral())
                if data2.Integral() != 0.0:
                    data2.Scale(1./data2.Integral())
                #data1.GetYaxis().SetRangeUser(0.001,1.1*data1.GetMaximum())
                data1.GetXaxis().SetLabelSize(0)

                legend.SetBorderSize(0)
                legend.SetFillStyle(0)
                legend.SetTextSize(0.04) # 0.035
                legend.AddEntry(data2, "pre HEM",  "l")
                legend.AddEntry(data1, "post HEM", "l") 
                
                data1.Draw("EP")
                data2.Draw("EP SAME")
                legend.Draw()
                # add CMS logo           
                add_CMSlogo(c1)            

                if data2.Integral() != 0.0:
                    # ---------------------------------------------------
                    # put the ratio plot between with HEM and without HEM
                    # ---------------------------------------------------
                    c1.cd(2)

                    ROOT.gPad.SetGridy()
                    ROOT.gPad.SetTopMargin(0.002)
                    ROOT.gPad.SetBottomMargin(0.35)
                    ROOT.gPad.SetLeftMargin(magicMargins_1D["L"])
                    ROOT.gPad.SetRightMargin(magicMargins_1D["R"])
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
        


