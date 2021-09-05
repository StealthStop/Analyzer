# ------------------------
# Significance calculation
# ------------------------
def cal_Significance(nSigEvents, nBkgEvents, sys=0.3):
    if (nBkgEvents == 0.0):
        return 0.0

    significance = nSigEvents / ( nBkgEvents + (sys * nBkgEvents)**2.0 )**0.5
    return significance

# -------------------------
# Closure error calculation
# -------------------------
def cal_ClosureError(nEvents_A, nEvents_B, nEvents_C, nEvents_D, nEventsErr_A, nEventsErr_B, nEventsErr_C, nEventsErr_D):

    if nEvents_A == 0.0 or nEvents_D == 0.0:
        return -999.0, -999.0
    
    closureError = abs(1.0 - ( (nEvents_B * nEvents_C) / (nEvents_A * nEvents_D) ) )
    
    closureErrUnc = ( ( ( nEvents_C * nEventsErr_B ) / ( nEvents_A * nEvents_D) )**2.0 
                    + ( ( nEvents_B * nEventsErr_C ) / ( nEvents_A * nEvents_D) )**2.0 
                    + ( ( nEvents_B * nEvents_C * nEventsErr_A ) / ( nEvents_A**2.0 * nEvents_D ) )**2.0 
                    + ( ( nEvents_B * nEvents_C * nEventsErr_D ) / ( nEvents_A * nEvents_D**2.0 ) )**2.0 )**0.5

    return closureError, closureErrUnc

# ---------------------
# ABCD prediction for A
# ---------------------
def cal_SimpleABCD(nEvents_A, nEvents_B, nEvents_C, nEvents_D, nEventsErr_A, nEventsErr_B, nEventsErr_C, nEventsErr_D):

    if nEvents_D == 0.0:
        return -999.0, -999.0
    
    nPred_A = (nEvents_B * nEvents_C) / nEvents_D
    nPred_Aunc = ((nEvents_C * nEventsErr_B / nEvents_D)**2.0 + (nEventsErr_C * nEvents_B / nEvents_D)**2.0 + (nEvents_C * nEvents_B * nEventsErr_D / nEvents_D**2.0)**2.0)**0.5

    return nPred_A, nPred_Aunc
