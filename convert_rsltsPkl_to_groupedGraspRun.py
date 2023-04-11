#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" 
NOTES:
AOD shows some correlation but with large error on 53 pixel test
    Cape_San_Juan costVal is always >1 (goes to 4+); Verde and Dakar always <1

X may not equal sum of pixels at above sites when script states "Packed X pixels with land or ocean data..." at end

"""

# TODO: GRASP v1.1.2 will not take three lognormal modes and merging from way back is not trivial
# TODO: auto adjust YAML breaks with blue noise by itself (see GitHub issue #6 in GRASP-Python-interface repo)
# TODO: adjust grouping so that runs are roughly even in size (e.g., 61 pixels does not become runs with 30, 30 & 1 pixels)

import sys
from os import path
import numpy as np
sys.path.append("/Users/wrespino/Synced/Local_Code_MacBook/GSFC-GRASP-Python-Interface")
sys.path.append("/Users/wrespino/Synced/Local_Code_MacBook/GSFC-Retrieval-Simulators")
from runGRASP import graspDB, graspRun, pixel
from simulateRetrieval import simulation

### SCRIPT INPUTS ###
# pklInputPath = '/Users/wrespino/Synced/RST_CAN-GRASP/AERONET_collocation/ABI16_ALM_TestFiles/ABI16_ALM_5lines_V1.pkl'
pklInputPath = '/Users/wrespino/Synced/RST_CAN-GRASP/AERONET_collocation/ABI16/ABI16_ALM.pkl'
# pklInputPath = '/Users/wrespino/Synced/RST_CAN-GRASP/AERONET_collocation/ABI17/ABI17_ALM.pkl'

# maxPixPerSite = None # ALL Pixels
maxPixPerSite = 100

# surfTypes = ['ocean', 'land'] # ALL Surfaces
surfTypes = ['ocean', 'land']

# siteNames = None # ALL Sites
siteNames = ['Capo_Verde', 'Cape_San_Juan', 'Dakar_Belair'] 

maxT = 30 # The maximum number of pixels to pack together in a single run of GRASP
measLwrBnd = 1e-6 # Minimum measurement values passed to GRASP

# GRASP processing options
maxCPU = 8
binPathGRASP = '/Users/wrespino/Synced/Local_Code_MacBook/grasp_open/build/bin/grasp'
krnlPathGRASP = '/Users/wrespino/Synced/Local_Code_MacBook/grasp_open/src/retrieval/internal_files'
yamlPathGRASP = '/Users/wrespino/Synced/Local_Code_MacBook/MODAERO_Analysis/YAML_files/settings_BCK_GRASPv1.1.2_2lgnrm_Ocean.yml' # easily split into land vs ocean or any other breakdown below ~LN95
savePath = '/Users/wrespino/Synced/Working/ABI_initialTests/Test_threeSites_Ocean_V02.pkl'
rndGuess = False
showWorkingPath = False # True prints a large amount to console

### END INPUTS – Begin Processing ###

# define functions
def packPixel(rs, surf='ocean'):
    if np.isnan(rs['meas_I_'+surf]).all(): return None # There was no data for this surface type
    # Create the pixel
    land_prct = 100 if surf=='land' else 0
    pix = pixel(rs['datetime'], 1, 1, rs['longitude'], rs['latitude'], rs['masl'], land_prct)
    # Update units to match GRASP UPDATE THESE VARIABLE NAMES TO COME FROM rs
    assert np.all(rs['fis']>=0), 'Negative phi detected in MODIS data but GRASP likes phi on the interval [0,360]' # FIXES: (1) comment this line out and accept less accuracy in GRASP RT OR (2) ensure phi<0 => phi=phi+360
    mu = np.cos(rs['sza']*np.pi/180) # [0] needed b/c sza might have dummyAng at end
    rs['meas_I_'+surf] = np.maximum(rs['meas_I_'+surf]*mu, measLwrBnd) # MODIS R=L/FO*pi/mu0; GRASP R=L/FO*pi w/ R>1e-6      
    # Populate pixel with observational data
    pix.populateFromRslt(rs, dataStage='meas', endStr=surf)
    return pix

def checkRsltConsistent(rsBase, rsNow):
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
    print(('AERONET site mismatch >> '+frmt+' || '+frmt) % (baseTup + nowTup))
    return False

def calc_maxTnow(n, maxT):
    r = 0.5 # r≈0 -> last run often very short; r≈1 -> maxTnow << maxT for many cases
    while True:
        if n % maxT == 0 or n % maxT > r*maxT: return maxT
        maxT -= 1

print('Loading collocated data from %s' % (path.basename(pklInputPath)))
gDB_in = graspDB()
rslts = gDB_in.loadResults(pklInputPath)
del gDB_in
for var in ['n','k']: # ensure 2D arrays with size 2xNλ for consistency with rslt dict (AERONET retrieve same value for fine/coarse RI)
    if var in rslts[0] and (rslts[0][var].ndim==1 or rslts[0][var].shape[0]==1):
        for i, rs in enumerate(rslts):
            rslts[i][var] = np.tile(rslts[i][var], (2,1))

print('Sorting rslts by observations datetime...')
rslts = rslts[np.argsort([r['datetime'] for r in rslts])]
for i in range(len(rslts)): rslts[i]['rsltsID'] = i # add unique rsltIDs to add in tracking later

print('Looping over all AERONET sites...') # This is quite fast, so no point in saving results before calling GRASP
dispString = 'Packing %5d pixels with land and/or ocean data at %26s' # len is 82 characters

grObjs = []
if siteNames is None: siteNames = np.unique([r['AERO_siteName'] for r in rslts])
for siteName in siteNames:
    mtchInds = [r['AERO_siteName']==siteName for r in rslts]
    NpixSite = sum(mtchInds) if maxPixPerSite is None else min(sum(mtchInds), maxPixPerSite)
    print(dispString % (NpixSite, siteName))
    grInd = {'ocean':-1, 'land':-1} # trigger creation of new GRASP run specific to this site
    maxTsite = dict()
    for surf in surfTypes: 
        inds = [~np.isnan(rs['meas_I_'+surf]).all() for rs in rslts[mtchInds][0:maxPixPerSite]]
        maxTsite[surf] = calc_maxTnow(sum(inds), maxT) # this is a bit of an art, see function above
    for rslt in rslts[mtchInds][0:maxPixPerSite]: # better to sub-select here b/c mtchInds is large boolean list
        checkRsltConsistent(rslts[mtchInds][0], rslt)
        for surfType in surfTypes:
            surfPix = packPixel(rslt, surfType)
            if surfPix is not None:
                surfPix.pixID = rslt['rsltsID']
                if grInd[surfType]<0 or len(grObjs[grInd[surfType]].pixels)==maxTsite[surf]:
                    grInd[surfType] = len(grObjs)
                    grObjs.append(graspRun(pathYAML=yamlPathGRASP, releaseYAML=True, verbose=showWorkingPath))
                grObjs[grInd[surfType]].addPix(surfPix)

print('Summarizing and sanity checking resulting list of graspRun objects...')
Npix = sum(len(gr.pixels) for gr in grObjs)
dispStringMod = dispString.replace('king','ked').replace('and/or', 'or') 
print(dispStringMod % (Npix,('%d sites' % len(siteNames))) + ' '*10) # extra spaces needed to ensure overwrite
for gr in grObjs:
    for pix in gr.pixels:
        assert np.isclose(gr.pixels[0].land_prct, pix.land_prct, atol=1), 'Land and ocean pixels in same run!'
Nocean = sum([grObj.pixels[0].land_prct < 1 for grObj in grObjs])
Nland = sum([grObj.pixels[0].land_prct > 99 for grObj in grObjs])
print('%d land runs and %d ocean runs for a total of %d runs packed' % (Nland, Nocean, len(grObjs)))

print('Creating graspDB object and running grasp...', flush=True)
gDB = graspDB(grObjs)
simABI = simulation()
rsltBck, failRun = gDB.processData(maxCPU, binPathGRASP, False, krnlPathGRASP, rndGuess=rndGuess)
simABI.rsltBck = rsltBck
print('%d of %d pixels processed successfully' % (sum(~failRun), Npix))
# del gDB # helpful in large scale processing, but not debugging
triedPixIDs = [px.pixID for grObj in grObjs for px in grObj.pixels]
simABI.rsltFwd = rslts[triedPixIDs][~failRun]
simABI.spectralInterpFwdToBck() # line up fwd with back wavelengths
simABI.saveSim(savePath, verbose=True)





























