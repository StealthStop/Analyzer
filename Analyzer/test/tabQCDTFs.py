import ROOT

def run(massExcPath, maxSignPath):
    opts     = [massExcPath, maxSignPath]
    channels = ["0l", "1l", "2l"]
    models   = ["SYY", "RPV"]
    
    lines = [r"\begin{tabular}{c c c c c c}", r"\hline", r"Year & channel & TF A (\%) & TF B (\%) & TF C (\%) & TF D (\%) \\", r"\hline"]
    
    for model in models:
        for opt in opts:
            if opt == opts[-1]:
                lines.append(r"\multicolumn{6}{c}{%s (significance-optimized ABCD regions)} \\"%(model))
            else:
                lines.append(r"\multicolumn{6}{c}{%s (limit-optimized ABCD regions)} \\"%(model))
    
            lines.append(r"\hline")
               
            for channel in channels:
                f = ROOT.TFile.Open("%s/Total_Run2UL_%s_%s_QCDCR_Prediction.root"%(opt, model, channel))
    
                h = f.Get("Run2UL_TF_%s_%sOver%s_%s_QCDCRABCD"%(model, channel, model, channel))
    
                chunks = ["Run2"]
                if channel == "0l":
                     chunks.append("fully-hadronic".ljust(14))
                elif channel == "1l":
                    chunks.append("semi-leptonic".ljust(14))
                elif channel == "2l":
                    chunks.append("fully-leptonic".ljust(14))
    
                for ibin in range(1, h.GetNbinsX()+1):
                    tf    = 100.0 * h.GetBinContent(ibin)
                    tferr = 100.0 * h.GetBinError(ibin)
                    if tf < 0.1 and tf != 0.0:
                        chunks.append((r"%.2f $\pm$ %.2f"%(tf, tferr)).ljust(13))
                    else:
                        chunks.append((r"%.1f $\pm$ %.1f"%(tf, tferr)).ljust(13))
    
                lines.append(" & ".join(chunks) + r" \\") 
    
            lines.append(r"\hline")
    lines.append(r"\end{tabular}")
    
    table = "\n".join(lines)
    
    tfTable = open("QCD_TF_Table.tex", "w")
    
    tfTable.write(table)
    
    tfTable.close()

if __name__ == "__main__":

    usage = "usage: %tabQCDTFs [options]"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument("--massExcPath", dest="massExcPath", help="Path to mass exclusion QCD CR root files",   default="NULL", required=True)
    parser.add_argument("--maxSignPath", dest="maxSignPath", help="Path to max significance QCD CR root files", default="NULL", required=True)
    args = parser.parse_args()

    run(args.massExcPath, args.maxSignPath)
