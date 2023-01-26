#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# TODO: This seems to be working okay so far but 0.47 land and ocean radiances are suspiciously close for first few pixels

import sys
from os import path
import numpy as np
sys.path.append("/Users/wrespino/Synced/Local_Code_MacBook/GSFC-GRASP-Python-Interface")
from runGRASP import graspDB, graspRun, pixel

# paths
# pklInputPath = '/Users/wrespino/Synced/RST_CAN-GRASP/AERONET_collocation/ABI16_ALM_TestFiles/ABI16_ALM_5lines_V1.pkl'
pklInputPath = '/Users/wrespino/Synced/RST_CAN-GRASP/AERONET_collocation/ABI16/ABI16_ALM.pkl'
# pklInputPath = '/Users/wrespino/Synced/RST_CAN-GRASP/AERONET_collocation/ABI17/ABI17_ALM.pkl'

# maxPixPerSite = None # ALL Pixels
maxPixPerSite = 10

surfTypes = ['ocean', 'land'] # ALL Surfaces
# surfTypes = ['ocean']

# siteNames = None # ALL Sites
siteNames = ['Capo_Verde', 'Cape_San_Juan', 'Dakar_Belair']

maxT = 60 # The maximum number of pixels to pack together in a single run of GRASP
measLwrBnd = 1e-6 # Minimum measurement values passed to GRASP 

# load the data
gDB_in = graspDB()
rslts = gDB_in.loadResults(pklInputPath)
del gDB_in

# (1) determine grouping of rslts
# HINT: We may want to condition siteIDs (i.e. make them integers)

def packPixel(rs, surf='ocean'):
    if np.isnan(rs['meas_I_'+surf]).all(): return None # There was no data for this surface type
    # Create the pixel
    land_prct = 100 if surf=='land' else 0
    pix = pixel(rs['datetime'], 1, 1, rs['longitude'], rs['latitude'], rs['masl'], land_prct)
    # Update units to match GRASP UPDATE THESE VARIABLE NAMES TO COME FROM rs
    assert np.all(rs['fis']>=0), 'Negative phi detected in MODIS data but GRASP likes phi on the interval [0,360]' # FIXES: (1) comment this line out and accept less accuracy in GRASP RT OR (2) ensure phi<0 => phi=phi+360
    mu = np.cos(rs['sza']*np.pi/180) # [0] needed b/c sza might have dummyAng at end
    rs['meas_I_'+surf] = np.maximum(rs['meas_I_'+surf]*mu, measLwrBnd) # MODIS R=L/FO*pi/mu0; GRASP R=L/FO*pi w/ R>1e-6      
    # Adjust rs variable names and populate pixel observational data
    rs['meas_'+surf+'_I'] = rs.pop('meas_I_'+surf) # this is format populateFromRslt() can handle; should have just started with it in convert_YingxiFile2_to_rsltsPkl.py
    pix.populateFromRslt(rs, dataStage='meas_'+surf)
    return pix

def checkRsltConsistent(rsBase, rsNow, lenDispStr=0):
    vld = []
    vld.append(np.isclose(rsBase['latitude'], rsNow['latitude'], atol=0.5))
    vld.append(np.isclose(rsBase['longitude'], rsNow['longitude'], atol=0.5))
    vld.append(np.isclose(rsBase['masl'], rsNow['masl'], atol=500))
    vld.append(rsBase['AERO_siteName'] == rsNow['AERO_siteName'])
    if np.all(vld): return True
    frmt = '%s [%d]: (%.2f,%.2f) - %dm'
    frmtKeys = ['AERO_siteName','AERO_siteID','latitude','longitude','masl']
    baseTup = tuple(rsBase[key] for key in frmtKeys)
    nowTup = tuple(rsNow[key] for key in frmtKeys)
    preStr = '\b'*lenDispStr + 'AERONET site mismatch >> '
    print((preStr + frmt + ' || ' + frmt) % (baseTup + nowTup), flush=True)
    return False

print('Sorting rslts by observations datetime...')
rslts = rslts[np.argsort([r['datetime'] for r in rslts])]
print('Looping over all AERONET sites...')
dispString = 'Packing %5d pixels at %26s' # len is 50 characters
lenDispStr = 50
print(dispString % (0,''), end='')
grObjs = []
if siteNames is None: siteNames = np.unique([r['AERO_siteName'] for r in rslts])
for siteName in siteNames:
    matchingRsltsInds = [r['AERO_siteName']==siteName for r in rslts]
    siteName = rslts[matchingRsltsInds][0]['AERO_siteName']
    inconsistRslt = False
    print('\b'*lenDispStr + (dispString % (sum(matchingRsltsInds), siteName)), end='', flush=True)
    grInd = {'ocean':-1, 'land':-1} # trigger creation of new GRASP run specific to this site
    for rslt in rslts[matchingRsltsInds][0:maxPixPerSite]:
        inconsistRslt = inconsistRslt or not checkRsltConsistent(rslts[matchingRsltsInds][0], rslt, lenDispStr)
        for surfType in surfTypes:
            surfPix = packPixel(rslt, surfType)
            if surfPix is not None:
                if grInd[surfType]<0 or len(grObjs[grInd[surfType]].pixels)==maxT:
                    grInd[surfType] = len(grObjs)
                    grObjs.append(graspRun(verbose=True))
                grObjs[grInd[surfType]].addPix(surfPix)
    if inconsistRslt: print(dispString % (0,''), end='')
# summarize and sanity check results
Npix = sum(len(gr.pixels) for gr in grObjs)
print('\b'*lenDispStr + dispString.replace('king','ked ') % (Npix,('%d sites' % len(siteNames)))) # extra space after ked needed to preserve length
for gr in grObjs:
    for pix in gr.pixels:
        assert np.isclose(gr.pixels[0].land_prct, pix.land_prct, atol=1), 'Land and ocean pixels in same run!'
Nocean = sum([grObj.pixels[0].land_prct < 1 for grObj in grObjs])
Nland = sum([grObj.pixels[0].land_prct > 99 for grObj in grObjs])
print('%d land runs and %d ocean runs for a total of %d runs packed' % (Nland, Nocean, len(grObjs)))


# (3) combine (1) and (2) to pack them all into graspRun members of gDB

# (4) save results to a pkl... gDB does not have a method for this (if above is fast we could just run here...)
