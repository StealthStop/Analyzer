import numpy as np

import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.lines as ml

plt.rcParams["font.family"] = "Nimbus Sans"

def addABCD(ax, c1, c2):
    plt.text(c1 + (1.0-c1)/2.0,  c2 + (1.0-c2)/2.0, "A", fontsize="58", fontweight="bold", va="center", ha="center")
    plt.text(c1/2.0,             c2 + (1.0-c2)/2.0, "B", fontsize="58", fontweight="bold", va="center", ha="center")
    plt.text(c1 + (1.0-c1)/2.0,  c2/2.0,            "C", fontsize="58", fontweight="bold", va="center", ha="center")
    plt.text(c1/2.0,             c2/2.0,            "D", fontsize="58", fontweight="bold", va="center", ha="center")

    l1 = ml.Line2D([c1, c1], [0.0, 1.0], color="black", linewidth=5, linestyle="dashed")
    l2 = ml.Line2D([0.0, 1.0], [c2, c2], color="black", linewidth=5, linestyle="dashed")

    ax.add_line(l1)
    ax.add_line(l2)

def plotDisc1vsDisc2(disc1, disc2, bw, c1, c2, tag, drawABCD = True, nBins = 100):

    fig = plt.figure(figsize=(7.5,8))

    plt.hist2d(disc1, disc2, bins=[nBins, nBins], range=[[0, 1], [0, 1]], cmap=plt.cm.jet, weights=bw, cmin = bw.min())
    ax = plt.gca()

    ax.set_ylabel("Variable 2", fontsize=30, fontweight="bold")
    ax.set_xlabel("Variable 1", fontsize=30, fontweight="bold")

    ax.tick_params(axis='both', which='major', labelsize=24)

    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(0.2))
    ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.2))

    if drawABCD:
        addABCD(ax, c1, c2)

    fig.tight_layout()

    fig.savefig("2D_%s_Disc1VsDisc2.pdf"%(tag), dpi=fig.dpi)

    plt.close(fig)

def plotVarVsBinEdges(var, edges, c1, c2, tag, drawABCD = True, nBins = 100):

    binWidth = 1.0 / nBins

    fig = plt.figure(figsize=(7.5,8))

    plt.hist2d(edges[:,0], edges[:,1], bins=[nBins+1, nBins+1], range=[[-binWidth/2.0, 1.0+binWidth/2.0], [-binWidth/2.0, 1.0+binWidth/2.0]], cmap=plt.cm.jet, weights=var, cmin=10e-10, cmax=20.0, vmin=0.0, vmax=0.5)
    ax = plt.gca()

    ax.set_ylabel("Horizontal Bin Edge", fontsize=30, fontweight="bold")
    ax.set_xlabel("Vertical Bin Edge",   fontsize=30, fontweight="bold")

    ax.tick_params(axis='both', which='major', labelsize=24)

    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(0.2))
    ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.2))

    if drawABCD:
        addABCD(ax, c1, c2)

    fig.tight_layout()

    fig.savefig("%s_vs_Disc1Disc2.pdf"%(tag), dpi=fig.dpi)

    plt.close(fig)

def main():
    
    bg1s = []; bg2s = []; wbg = []
    sg1s = []; sg2s = []; wsg = []
    
    theEdge1 = 0.6
    theEdge2 = 0.6
    nBins    = 100
    drawABCD = True

    edges = np.arange(0.0, 1.0, 1.0 / nBins)

    for i in range(0, 1000000):
    
        bg1 = np.random.exponential(1.0)
        bg2 = np.random.exponential(1.0)
    
        sg1 = np.random.normal(0.95, 0.2) 
        sg2 = np.random.normal(0.95, 0.2) 
    
        if bg1 >= 0.0 and bg1 <= 1.0 and bg2 >= 0.0 and bg2 <= 1.0: 
            bg1s.append(bg1); bg2s.append(bg2); wbg.append(1.0)
    
        if sg1 >= 0.0 and sg1 <= 1.0 and sg2 >= 0.0 and sg2 <= 1.0 and i%35==0:
            sg1s.append(sg1); sg2s.append(sg2); wsg.append(1.0)
    
    bg1sArr = np.array(bg1s)
    bg2sArr = np.array(bg2s)
    sg1sArr = np.array(sg1s)
    sg2sArr = np.array(sg2s)
    wsgArr  = np.array(wsg)
    wbgArr  = np.array(wbg)
    
    dg1sArr = np.array(bg1s+sg1s)
    dg2sArr = np.array(bg2s+sg2s)
    wgArr   = np.array(wbg+wsg)
    
    plotDisc1vsDisc2(bg1sArr, bg2sArr, wbgArr, theEdge1, theEdge2, "BG", drawABCD, nBins)
    plotDisc1vsDisc2(sg1sArr, sg2sArr, wsgArr, theEdge1, theEdge2, "SG", drawABCD, nBins)
    plotDisc1vsDisc2(dg1sArr, dg2sArr, wgArr,  theEdge1, theEdge2, "SB", drawABCD, nBins)
    
    nonclosures = [] 
    nonclosuress = [] 
    theEdges = [] 
    
    for edge1 in edges:
        for edge2 in edges:
    
            N = float(np.count_nonzero(wbgArr))
            A = float(np.count_nonzero((bg1sArr>edge1)&(bg2sArr>edge2)))
            B = float(np.count_nonzero(bg2sArr>edge2)) - A
            C = float(np.count_nonzero(bg1sArr>edge1)) - A
            D = N - A - B - C
    
            Ns = float(np.count_nonzero(wsgArr))
            As = float(np.count_nonzero((sg1sArr>edge1)&(sg2sArr>edge2)))
            Bs = float(np.count_nonzero(sg2sArr>edge2)) - As
            Cs = float(np.count_nonzero(sg1sArr>edge1)) - As
            Ds = Ns - As - Bs - Cs
    
            if B > 0.0 and C > 0.0:
    
                nonclosures.append(abs(1. - (A*D)/(B*C)))
                nonclosuress.append(abs(1. - ((A+As)*(D+Ds))/((B+Bs)*(C+Cs))))
                theEdges.append([edge1, edge2])
    
    plotVarVsBinEdges(np.array(nonclosures),  np.array(theEdges), theEdge1, theEdge2, "BG_NonClosure", drawABCD, nBins)
    plotVarVsBinEdges(np.array(nonclosuress), np.array(theEdges), theEdge1, theEdge2, "SB_NonClosure", drawABCD, nBins)

if __name__ == "__main__":
    main()
