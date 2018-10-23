#!/usr/bin/env python3

import matplotlib.pyplot as plt
import datetime as dt
import numpy as np
import dill
import os
import sys
sys.path.append(os.path.join("..", "GRASP_PythonUtils"))
from miscFunctions import angstrmIntrp

#dillFN = '/Users/wrespino/Synced/Working/MODAERO_16binLoose.pkl'
dillFN = '/Users/wrespino/Synced/Working/MODAERO_OceanMODall_Uranus.pkl'

maxAOD = 2

siteID_plot = 106 # False for all sites (doesn't apply to AOD)
ttleNm = 'Ascension Island'
lmbdInd = 3 # [0.469, 0.555, 0.645, 0.8585, 1.24, 1.64, 2.13]
AODonly = True

dill.load_session(dillFN)
radii = rslts[0][0]['r'][0] # assume the same radii are used in every retrieval; analysis:ignore (it was loaded)
wvlngth = rslts[0][0]['lambda'] # assume the same lambda used in every retrieval; analysis:ignore (it was loaded)
Nmodes = rslts[0][0]['r'].shape[0] # analysis:ignore (it was loaded)
Npsd = radii.shape[0]
Nlmbd = wvlngth.shape[0] 
aodGRASP = np.array([])
rriGRASP = np.array([])
iriGRASP = np.array([])
h20Alb550 = np.array([])
h20AlbAll = np.array([]).reshape(0, Nlmbd) # index 1 will be 550nm
wndSpd = np.array([])
dateDays = np.array([]) 
dVdlnr = np.array([]).reshape(0, Npsd)
rv = np.array([]).reshape(0, Nmodes)
sigma = np.array([]).reshape(0, Nmodes)
sphGRASP = np.array([]).reshape(0, Nmodes)
# MODIS GRASP AODs
for rsltRun in rslts: # analysis:ignore (it was loaded)
    for rslt in rsltRun:
        aodGRASP = np.append(aodGRASP, rslt['aod'][lmbdInd])
        rriGRASP = np.append(rriGRASP, rslt['n'][lmbdInd]) # BUG: the following three are only fine mode values
        iriGRASP = np.append(iriGRASP, rslt['k'][lmbdInd]) #      fix needed here and in reading of GRASP output
        sphGRASP = np.vstack([sphGRASP, rslt['sph']])
        dVdlnr = np.vstack([dVdlnr, rslt['vol']@rslt['dVdlnr']]) # scale to um^3/um^2 & sum multiple modes
        h20AlbAll = np.vstack([h20AlbAll, rslt['wtrSurf'][0]])    
        wndSpd = np.append(wndSpd, (rslt['wtrSurf'][2][1]-0.003)/0.00512)
        tDelta = (rslt['datetime'] - dt.datetime(2000, 1, 1, 0, 0))
        dateDays = np.append(dateDays, tDelta.days + tDelta.seconds/86400) # days since start of 2000's
        rv = np.vstack([rv, rslt['rv']])
        sigma = np.vstack([sigma, rslt['sigma']])
nPixels = aodGRASP.shape[0]
lmbdVal = rslt['lambda'][lmbdInd]
lmbdStr = ' (%4.2f μm)' % lmbdVal

# AERONET AODs
aodAERO = np.array([])
siteID = np.array([])
for seg in DB.siteSegment: # analysis:ignore (it was loaded)
    unqDTs = np.unique(seg.mod_loc[:,-1])
    for unqDT in unqDTs:
        nowInd = np.nonzero(seg.mod_loc[:,-1] == unqDT)[0][0]
        if np.abs(lmbdVal - 0.555) < 0.01: # use values in file
            aodAERO = np.append(aodAERO, seg.aod[nowInd,-1])
        else:
            aodAERO_lmbd = angstrmIntrp(seg.AERO_LAMDA, seg.aod[nowInd,:], lmbdVal)
            aodAERO = np.append(aodAERO, aodAERO_lmbd)            
        siteID = np.append(siteID, seg.aero_loc[0])
assert (aodAERO.shape[0] == aodGRASP.shape[0]), 'rslts and DB.siteSegment contained a different number of pixels!'

# compare with AERONET
aodAERO = np.array(aodAERO)
aodGRASP = np.array(aodGRASP)
aodGRASPcln = aodGRASP[~np.isnan(aodAERO)]

clrVar = wndSpd[~np.isnan(aodAERO)] # WIND SPEED
#clrVar = h20AlbAll[~np.isnan(aodAERO), 0] # OCEAN ALBEDO, last ind. wavelength
clrVar[clrVar < np.percentile(clrVar,2)] = np.percentile(clrVar,2)
clrVar[clrVar > np.percentile(clrVar,98)] = np.percentile(clrVar,98)
#clrVar = siteID[~np.isnan(aodAERO)]
aodAERO = aodAERO[~np.isnan(aodAERO)]

line = plt.figure()
plt.scatter(aodAERO, aodGRASPcln, c=clrVar, marker='.')
plt.plot(np.r_[0, maxAOD], np.r_[0, maxAOD], 'k')
plt.xlabel("AOD AERONET" + lmbdStr)
plt.ylabel("AOD MODIS/GRASP" + lmbdStr)
plt.xlim([0, maxAOD])
plt.ylim([0, maxAOD])
Rcoef = np.corrcoef(aodAERO, aodGRASPcln)[0,1]
RMSE = np.sqrt(np.mean((aodAERO - aodGRASPcln)**2))
textstr = 'N=%d\nR=%.3f\nRMS=%.3f\n'%(len(aodAERO), Rcoef, RMSE)
plt.text(0.7*maxAOD, 0.03, textstr, fontsize=12)
plt.title('All Sites')

if AODonly: sys.exit()
# plot PSD
N = len(dVdlnr)
if siteID_plot:
    pltInd = siteID == siteID_plot
    Nplt = np.sum(pltInd)
else:
    pltInd = slice(None)
    Nplt = N
line = plt.figure()
ttleStr = ''
if (Nplt < 100):
    [plt.semilogx(radii, sngdV, color=[i/Nplt, 0.2, 1/(i+1)], linewidth=0.3) for i,sngdV in enumerate(dVdlnr)]
    ttleStr = ' (blue=old, red=recent, black=mean)'
plt.semilogx(radii, np.mean(dVdlnr[pltInd,:], axis=0), 'k')
plt.title(ttleNm + ttleStr)
plt.xlabel("radius (um)")
plt.ylabel("dV/dlnr (um^3/um^2)")

# plot RRI PDF
fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(10, 5))
ax[0].hist(rriGRASP[pltInd], 30)
ax[0].set_xlabel("RRI" + lmbdStr)
ax[0].set_ylabel("Frequency")
# plot IRI PDF
ax[1].hist(iriGRASP[pltInd], 60)
ax[1].set_title(ttleNm)
ax[1].set_xlabel("IRI" + lmbdStr)
# plot SPH PDF
ax[2].hist(sphGRASP[pltInd,0], 30) # ONLY FINE MODE
ax[2].set_xlabel("Spherical Fraction")

# plot H20 Albedo Time Series
fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(5, 10))
ax[0].plot(dateDays[pltInd], h20AlbAll[pltInd,1])
#plt.plot(dateDays, aodGRASP/3)
ax[0].set_ylabel("Ocean Albedo Term" + lmbdStr)
ax[0].set_title(ttleNm + " (fine mode only)")
# plot RRI Time Series
ax[1].plot(dateDays[pltInd], rriGRASP[pltInd])
ax[1].set_ylabel("Real Refractive Index" + lmbdStr)
# plot Rv and sigma Time Series
ax[2].plot(dateDays[pltInd], rv[pltInd,:], dateDays[pltInd], sigma[pltInd,:])
ax[2].set_yscale('log')
ax[2].set_xlabel("Days since 2000")
ax[2].set_ylabel("rv & σ's")


# plot H20 Albedo vs AOD
fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(10, 5))
ax[0].plot(h20AlbAll[:,1], aodGRASP, '.')
ax[0].set_ylabel("GRASP AOD" + lmbdStr)
ax[0].set_xlabel("Ocean Albedo Term")
# plot wind speed vs aod
ax[1].plot(wndSpd, aodGRASP, '.')
#ax[1].plot(wndSpdCln, aodAERO, '.')
ax[1].set_title('All Sites')
ax[1].set_ylabel("GRASP AOD" + lmbdStr)
ax[1].set_xlabel("GRASP Wind Speed (m/s)")
ax[1].legend(('GRASP','AERONET'))
# plot wind speed vs ocean albedo
ax[1].scatter(wndSpd, h20AlbAll[:,1], c=aodGRASP)
ax[1].set_ylabel("Ocean Albedo Term")
ax[1].set_xlabel("GRASP Wind Speed (m/s)")

# plot ocean albedo vs wavelength
line = plt.figure()
[plt.scatter(wvlngth, alb, color=[i/Nplt, 0.2, 1/(i +1)], s=1) for i,alb in enumerate(h20AlbAll[pltInd,:])]
plt.plot(wvlngth, np.mean(h20AlbAll[pltInd,:], axis=0), 'k')
plt.title(ttleNm)
plt.ylabel("Ocean Albedo Term")
plt.xlabel("Wavelength (um)")