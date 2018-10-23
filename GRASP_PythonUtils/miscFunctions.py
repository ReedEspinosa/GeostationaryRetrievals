#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

def angstrmIntrp(lmbd, tau, lmbdTrgt):
    tau = tau[lmbd.argsort()]
    lmbd.sort()
    frstInd = np.nonzero((lmbd - lmbdTrgt) < 0)[0][-1]
    alpha = angstrm(lmbd[frstInd:frstInd+2], tau[frstInd:frstInd+2])
    return tau[frstInd]*(lmbd[frstInd]/lmbdTrgt)**alpha
    
def angstrm(lmbd, tau):
    assert (lmbd.shape[0]==2 and tau.shape[0]==2), "Fitting of more than two values not supported"
    return -np.log(tau[0]/tau[1])/np.log(lmbd[0]/lmbd[1])