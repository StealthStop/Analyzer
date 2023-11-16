import ROOT
import argparse

def run(inputPath):

    tag = inputPath.split("/")[-1]

    channels = ["0l", "1l", "2l"]
    lines = [r"\begin{tabular}{c c c c c c c}", r"\hline", r"Year & channel & \multicolumn{5}{c}{\njets (\%)} \\", r"\hline"]

    for channel in channels:

        njetsChunks = []
        njetsRange = None
        if   channel == "0l":
            njetsRange = list(range(8, 13))
        elif channel == "1l":
            njetsRange = list(range(7, 12))
        elif channel == "2l":
            njetsRange = list(range(6, 11))

        for Njets in njetsRange:
            if Njets == njetsRange[-1]:
                njetsChunks.append((r"$\geq %d$"%(Njets)).ljust(16))
            else:
                njetsChunks.append((r"$= %d$"%(Njets)).ljust(16))
            
        lines.append(r"\multicolumn{7}{c}{} \\\hline")

        lines.append("     & " + " ".ljust(14) + " & %s"%(" & ".join(njetsChunks)) + r"\\")

        # For the table, can grab any output since the TFs inside will be the same
        f = None 
        model = ""
        try:
            f = ROOT.TFile.Open("%s/Total_Run2UL_RPV_%s_QCDCR_Prediction.root"%(inputPath, channel))
            model = "RPV"
        except:
            f = ROOT.TFile.Open("%s/Total_Run2UL_SYY_%s_QCDCR_Prediction.root"%(inputPath, channel))
            model = "SYY"
        
        h = f.Get("Run2UL_TF_%s_%sOver%s_%s_QCDCRABCD_perNjets"%(model, channel, model, channel))

        chunks = ["Run2"]
        if channel == "0l":
             chunks.append("fully-hadronic".ljust(14))
        elif channel == "1l":
            chunks.append("semi-leptonic".ljust(14))
        elif channel == "2l":
            chunks.append("fully-leptonic".ljust(14))
        
        for ibin in range(1, h.GetNbinsX() / 4 + 1):
                
            tf    = 100.0 * h.GetBinContent(ibin)
            tferr = 100.0 * h.GetBinError(ibin)
            if tf < 0.1 and tf != 0.0:
                chunks.append((r"%.2f (%.0f)"%(tf, 100.0*tferr/tf)).ljust(16))
            else:
                chunks.append((r"%.1f (%.0f)"%(tf, 100.0*tferr/tf)).ljust(16))

        lines.append(" & ".join(chunks) + r" \\") 
       
        lines.append(r"\hline")
        
    lines.append(r"\end{tabular}")

    table = "\n".join(lines)
    
    tfTable = open("QCD_TF_Table_%s.tex"%(tag), "w")
    
    tfTable.write(table)
    
    tfTable.close()

if __name__ == "__main__":

    usage = "usage: %tabQCDTFs [options]"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("--inputPath", dest="inputPath", help="Path to QCD CR root files", default="NULL",  required=True)
    args = parser.parse_args()

    run(args.inputPath)
