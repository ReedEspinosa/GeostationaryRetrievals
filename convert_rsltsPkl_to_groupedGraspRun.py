#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from os import path
import numpy as np
sys.path.append("/Users/wrespino/Synced/Local_Code_MacBook/GSFC-GRASP-Python-Interface")
from runGRASP import graspDB, graspRun, pixel

# paths
pklInputPath = '/Users/wrespino/Synced/RST_CAN-GRASP/AERONET_collocation/ABI16_ALM_TestFiles/ABI16_ALM_5lines_V1.pkl'
# pklInputPath = '/Users/wrespino/Synced/RST_CAN-GRASP/AERONET_collocation/ABI16/ABI16_ALM.pkl'
# pklInputPath = '/Users/wrespino/Synced/RST_CAN-GRASP/AERONET_collocation/ABI17/ABI17_ALM.pkl'

maxT = 60

# load the data
gDB_in = graspDB()
rslts = gDB_in.loadResults(pklInputPath)
del gDB_in

# (1) determine grouping of rslts
# HINT: We may want to condition siteIDs (i.e. make them integers)

def packPixel(rs, surf='ocean'):
    print('We still need to right this one...')
    if np.isnan(rslts[0]['meas_I_'+surf]).all(): return None
    land_prct = 100 if surf=='land' else 0 
    pix = pixel(rs['datetime'], 1, 1, s['longitude'], rs['latitude'], rs['masl'], land_prct)
    # Prep rs against pixel.populateFromRslt; careful of changes in calling scope (will reuse for land)
    pix.populateFromRslt(rs, dataStage=meas) # a cheap hack would be to store measurements we want (land/ocean) as the fits (actually any string will do)
    return pix

print('Sorting rslts by observations datetime...')
rslts = rslts[np.argsort([r['datetime'] for r in rslts])]
print('Looping over all AERONET sites...')
grObjs = []
grInd['ocean'] = -1
grInd['land'] = -1
for siteID in np.unique([r['AERO_siteID'] for r in rslts]):
    matchingRsltsInds = [r['AERO_siteID']==sitID for r in rslts] 
    for rslt in rslts[matchingRsltsInds]:
        for surfType in ['ocean', 'land']:
            surfPix = packPixel(rslt, surfType)
            if surfPix is not None:
                if grInd[surfType]<0 or len(grObjs[grInd[surfType]].pixels)==maxT:
                    grInd[surfType] = len(grObjs)
                    grObjs.append(graspRun(verbose=True))
                grObjs[grInd[surfType]].addPix(oceanPix)
        
# (3) combine (1) and (2) to pack them all into graspRun members of gDB

# (4) save results to a pkl... gDB does not have a method for this (if above is fast we could just run here...)
