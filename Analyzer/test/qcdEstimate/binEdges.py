#!/bin/python

#Create a file with the latest and greatest bin edges that can be easily ported to all files.
#Region is defined as:
#       SRBE - bin edges derived in TT MC in the SR
#       CRBE - bin edges derived in Data from the SingleMuon data set in the CR

def binEdgeDictReturn( region, year ) :
    
    binEdges = {}
   
    # These are the bin edges derived from the QCD control region using data from the Single Muon data set - Updated January 30, 2020

    if region == "CRBE" and year == "2018post" :
        binEdges[7]     = [ 0.364, 0.729, 0.842, 1.000 ]
        binEdges[8]     = [ 0.386, 0.771, 0.874, 1.000 ]
        binEdges[9]     = [ 0.405, 0.798, 0.895, 1.000 ]
        binEdges[10]    = [ 0.411, 0.826, 0.914, 1.000 ]
        binEdges[11]    = [ 0.411, 0.826, 0.914, 1.000 ]
        binEdges[12]    = [ 0.411, 0.826, 0.914, 1.000 ]
    
    if region == "CRBE" and year == "2018pre" :
        binEdges[7]     = [ 0.345, 0.713, 0.832, 1.000 ]
        binEdges[8]     = [ 0.363, 0.752, 0.863, 1.000 ]
        binEdges[9]     = [ 0.374, 0.779, 0.883, 1.000 ]
        binEdges[10]    = [ 0.386, 0.805, 0.900, 1.000 ]
        binEdges[11]    = [ 0.386, 0.805, 0.911, 1.000 ]
        binEdges[12]    = [ 0.386, 0.805, 0.911, 1.000 ]
    
    if region == "CRBE" and year == "2017" :
        binEdges[7]     = [ 0.391, 0.769, 0.869, 1.000 ]
        binEdges[8]     = [ 0.413, 0.807, 0.892, 1.000 ]
        binEdges[9]     = [ 0.431, 0.831, 0.913, 1.000 ]
        binEdges[10]    = [ 0.437, 0.862, 0.932, 1.000 ]
        binEdges[11]    = [ 0.444, 0.879, 0.942, 1.000 ]
        binEdges[12]    = [ 0.484, 0.900, 0.961, 1.000 ]
                
    if region == "CRBE" and year == "2016" :
        binEdges[7]     = [ 0.407, 0.768, 0.885, 1.000 ]
        binEdges[8]     = [ 0.428, 0.811, 0.907, 1.000 ]
        binEdges[9]     = [ 0.454, 0.847, 0.933, 1.000 ]
        binEdges[10]    = [ 0.466, 0.860, 0.937, 1.000 ]
        binEdges[11]    = [ 0.505, 0.892, 0.948, 1.000 ]
        binEdges[12]    = [ 0.494, 0.911, 0.960, 1.000 ]
            
    # These bin edges are derived from the signal region using TT bar MC - Updated January 29,2020
    if region == "SRBE" and year == "2018post" :
        binEdges[7]     = [ 0.364, 0.729, 0.842, 1.000 ]
        binEdges[8]     = [ 0.386, 0.771, 0.874, 1.000 ]
        binEdges[9]     = [ 0.405, 0.798, 0.895, 1.000 ]
        binEdges[10]    = [ 0.411, 0.826, 0.914, 1.000 ]
        binEdges[11]    = [ 0.411, 0.826, 0.914, 1.000 ]
        binEdges[12]    = [ 0.411, 0.826, 0.914, 1.000 ]

    if region == "SRBE" and year == "2018pre" :
        binEdges[7]     = [ 0.345, 0.713, 0.832, 1.000 ]
        binEdges[8]     = [ 0.363, 0.752, 0.863, 1.000 ]
        binEdges[9]     = [ 0.374, 0.779, 0.883, 1.000 ]
        binEdges[10]    = [ 0.386, 0.805, 0.900, 1.000 ]
        binEdges[11]    = [ 0.386, 0.805, 0.911, 1.000 ]
        binEdges[12]    = [ 0.386, 0.805, 0.911, 1.000 ]

    if region == "SRBE" and year == "2017" :
        binEdges[7]     = [ 0.350, 0.716, 0.833, 1.000 ]
        binEdges[8]     = [ 0.368, 0.753, 0.864, 1.000 ]
        binEdges[9]     = [ 0.381, 0.780, 0.885, 1.000 ]
        binEdges[10]    = [ 0.387, 0.793, 0.894, 1.000 ]
        binEdges[11]    = [ 0.402, 0.806, 0.898, 1.000 ]
        binEdges[12]    = [ 0.402, 0.806, 0.898, 1.000 ]
            
    if region == "SRBE" and year == "2016" :
        binEdges[7]     = [ 0.349, 0.680, 0.835, 1.000 ]
        binEdges[8]     = [ 0.358, 0.741, 0.878, 1.000 ]
        binEdges[9]     = [ 0.371, 0.787, 0.901, 1.000 ]
        binEdges[10]    = [ 0.386, 0.816, 0.917, 1.000 ]
        binEdges[11]    = [ 0.410, 0.848, 0.936, 1.000 ]
        binEdges[12]    = [ 0.410, 0.853, 0.936, 1.000 ]

    return binEdges
