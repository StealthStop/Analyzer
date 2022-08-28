# Double DisCo ABCD Validation Framework

This framework is for performing the validation of the Double DisCo based ABCD method.
The code analyzes 2D background and signal histograms and computes quantities about the method.
Files beginning with `common_` contain core classes for performing fundamental calculations for any validation study, while
files starting with `DoubleDisCo_` are ``applications'' that perform different validation studies.

## `common` Modules

There are four core `common` modules: `common_Regions.py`, `common_Aggregator.py`, `common_Plotter.py`, `common_TableWriter.py`.

`common_Regions.py`: This module provides an `All_Regions` class, which takes in 2D histograms for background, signal processes as well as data.
The exact bounds of a 2D region to perform the ABCD method in are provided by the user, and the method is then performed for every choice of edge 
values within the bounded region to define the A, B, C, and D subregions.
With this scanning of all choices of edge values within the bounded 2D region, quantities such as ABCD closure, MC correction factors etc., are calculated.
Getter functions are provided to extract any of these quantities for any choice of the edges.

`common_Aggregtator.py`: This module supplies an `Aggregator` class to combine information from multiple instances of the `All_Regions` for extracting useful information in bulk,
i.e. for mutliple choices of the definition of a 2D region.

`common_Plotter.py`: The module supplies the `Common_Calculations_Plotters` class, which contains all custom plotting functions.
Quantities are aggregated from the `All_Regions` objects and curated for sending to this plotter class.

`common_TableWriter.py`: A module that provides a `TableWriter` class that has curated information from the `All_Regions` class passed to it for making custom LaTeX tables.

`common_HiggsCombineInputs.py`: A module that provides functions to make a root file including all sys for Higgs combine.

## `DoubleDisCo` Applications

Currently, there are four ``applications'' that perform specific validation studies: `DoubleDisCo_BinEdges.py`, `DoubleDisCo_MCcorrectionFactor_TT.py`, `DoubleDisCo_MCcorrectionFactor_TTvar.py`, `DoubleDisCo_Optimized_BinEdges.py`.

`DoubleDisCo_BinEdges.py`: This application performs a general validation where several quantities are plotted as a function of ABCD edges in both 1D and 2D.
Three validation regions (subdivision of BD, subdivision of CD, and subdivision of D) are defined with predefined boundary values, and quantities are plotted
individually for these validation boundaries.

`DoubleDisCo_Optimized_BinEdges.py`: This application specifically performs an optimization procedure to determine what the best choice of edges defining the proper ABCD regions.
The best choice is determined based on best signal sensitivy as well as minimal non-closure, and limited signal contamination in the B, C, and D regions.
The top `n` choices (requested by the user) are printed in a LaTeX table.

`DoubleDisCo_MCcorrectionFactor_TT.py`: This application plots quantities related a correction factor to be applied to the ABCD calculation.
Things are plotted as a function of the boundary value defining the validation regions. 

`DoubleDisCo_MCcorrectionFactor_TTvar.py`: The ``TTvar'' version of the `MCcorrectionFactor` application plots quantities related a correction factor to be applied to the ABCD calculation for several variations of the TT background e.g. TT_erdOn, TT_fsrUp, TT_fsrDown, etc.

## Performing Validation Studies

Any validation study is run with the main `run_DoubleDisCo_Validation.py` script
where the name of the application is passed to the `--run` argument.

```
usage: usage: %prog [options] [-h] --run RUN --year YEAR [--path PATH]
                              [--tt TT] [--nontt NONTT] [--ttVar TTVAR]
                              [--sig SIG] [--mass MASS] [--data DATA]
                              --channel CHANNEL [--disc1edge DISC1EDGE]
                              [--disc2edge DISC2EDGE] [--fastMode]
                              [--njets NJETS [NJETS ...]] [--plotVars1D]
                              [--plotVars2D] [--plotDisc1VsDisc2]
                              [--plotVarVsBoundary]
                              [--numEdgeChoices NUMEDGECHOICES]

optional arguments:
  -h, --help            show this help message and exit
  --run RUN             which code to run
  --year YEAR           which year
  --path PATH           Input dir with histos
  --tt TT               name of TT sample
  --nontt NONTT         name of NonTT sample
  --ttVar TTVAR         TT MCcorrectionFactor_TTvar (default no var)
  --sig SIG             signal model RPV, SYY
  --mass MASS           signal mass
  --data DATA           name of Data sample
  --channel CHANNEL     0l or 1l
  --disc1edge DISC1EDGE
                        fixed d1 edge
  --disc2edge DISC2EDGE
                        fixed d2 edge
  --fastMode            Fast mode, don't scan all choices
  --njets NJETS [NJETS ...]
                        which njet bins to run on
  --plotVars1D          Plot 1D var vs disc (slices)
  --plotVars2D          Plot var vs disc1 and disc2 (2D)
  --plotDisc1VsDisc2    Plot disc1 and disc2 (2D)
  --plotVarVsBoundary   Plot var vs boundary
  --numEdgeChoices NUMEDGECHOICES
                        number of edge choices to print
```
Various useful examples for running the different applications with command-line arguments.

```
# Running for Optimized_BinEdges
python run_DoubleDisCo_Validation.py --run Optimized_BinEdges --year 2016 --channel 0l
python run_DoubleDisCo_Validation.py --run Optimized_BinEdges --year 2016 --channel 1l

# Running for MCcorrectionFactor_TT
python run_DoubleDisCo_Validation.py --run MCcorrectionFactor_TT --year 2016 --channel 0l --disc1edge 0.69 --disc2edge 0.74 --fastMode
python run_DoubleDisCo_Validation.py --run MCcorrectionFactor_TT --year 2016 --channel 1l --disc1edge 0.59 --disc2edge 0.70 --fastMode

# Running for MCcorrectionFactor_TTvar
python run_DoubleDisCo_Validation.py --run MCcorrectionFactor_TTvar --year 2016 --channel 0l  --disc1edge 0.69 --disc2edge 0.74
python run_DoubleDisCo_Validation.py --run MCcorrectionFactor_TTvar --year 2016 --channel 1l  --disc1edge 0.59 --disc2edge 0.70

# Running for BinEdges
python run_DoubleDisCo_Validation.py --run BinEdges --year 2016 --channel 0l --disc1edge 0.6 --disc2edge 0.6
python run_DoubleDisCo_Validation.py --run BinEdges --year 2016 --channel 1l --disc1edge 0.6 --disc2edge 0.6

python run_DoubleDisCo_Validation.py --run BinEdges --year 2016 --channel 0l --disc1edge 0.69 --disc2edge 0.74
python run_DoubleDisCo_Validation.py --run BinEdges --year 2016 --channel 1l --disc1edge 0.59 --disc2edge 0.70
```
