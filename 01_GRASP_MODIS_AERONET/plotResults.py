#!/usr/bin/env python3

import matplotlib.pyplot as plt
import datetime as dt
import numpy as np
import dill

#dillFN = '/Users/wrespino/Synced/Working/MODAERO_16binLoose.pkl'
dillFN = '/Users/wrespino/Synced/Working/MODAERO_OceanMODall_Uranus.pkl'

maxAOD = 0.3
ttleNm = 'Nauru'

dill.load_session(dillFN)
radii = rslts[0][0]['r'][0] # assume the same radii are used in every retrieval; analysis:ignore (it was loaded)
wvlngth = rslts[0][0]['lambda'] # assume the same lambda used in every retrieval; analysis:ignore (it was loaded)
Nmodes = rslts[0][0]['r'].shape[0] # analysis:ignore (it was loaded)
Npsd = radii.shape[0]
Nlmbd = wvlngth.shape[0] 
aod550AERO = np.array([])
aod550GRASP = np.array([])
rri550GRASP = np.array([])
iri550GRASP = np.array([])
sph550GRASP = np.array([])
h20Alb550 = np.array([])
h20AlbAll = np.array([]).reshape(0, Nlmbd) # index 1 will be 550nm
wndSpd = np.array([])
dateDays = np.array([]) 
dVdlnr = np.array([]).reshape(0, Npsd)
rv = np.array([]).reshape(0, Nmodes)
sigma = np.array([]).reshape(0, Nmodes)
# MODIS GRASP AODs
for rsltRun in rslts: # analysis:ignore (it was loaded)
    for rslt in rsltRun:
        aod550GRASP = np.append(aod550GRASP, rslt['aod'][1])
        rri550GRASP = np.append(rri550GRASP, rslt['n'][1]) # BUG: the following three are only fine mode values
        iri550GRASP = np.append(iri550GRASP, rslt['k'][1]) #      fix needed here and in reading of GRASP output
        sph550GRASP = np.append(sph550GRASP, rslt['sph'])
        dVdlnr = np.vstack([dVdlnr, rslt['vol']@rslt['dVdlnr']]) # scale to um^3/um^2 & sum multiple modes
        h20AlbAll = np.vstack([h20AlbAll, rslt['wtrSurf'][0]])
        
        wndSpd = np.append(wndSpd, (rslt['wtrSurf'][2][1]-0.003)/0.00512)
        tDelta = (rslt['datetime'] - dt.datetime(2000, 1, 1, 0, 0))
        dateDays = np.append(dateDays, tDelta.days + tDelta.seconds/86400) # days since start of 2000's
        rv = np.vstack([rv, rslt['rv']])
        sigma = np.vstack([sigma, rslt['sigma']])
   
# AERONET AODs
for seg in DB.siteSegment: # analysis:ignore (it was loaded)
    unqDTs = np.unique(seg.mod_loc[:,-1])
    for unqDT in unqDTs:
        nowInd = np.nonzero(seg.mod_loc[:,-1] == unqDT)[0][0]
        aod550AERO = np.append(aod550AERO, seg.aod[nowInd,-1])
assert (aod550AERO.shape[0] == aod550GRASP.shape[0]), 'rslts and DB.siteSegment contained a different number of pixels!'

aod550AERO = np.array(aod550AERO)
aod550GRASP = np.array(aod550GRASP)
aod550GRASPcln = aod550GRASP[~np.isnan(aod550AERO)]
wndSpdCln = wndSpd[~np.isnan(aod550AERO)]
h20Alb550Cln = h20AlbAll[~np.isnan(aod550AERO),1]
aod550AERO = aod550AERO[~np.isnan(aod550AERO)]


# compare with AERONET
line = plt.figure()
#plt.plot(aod550AERO, aod550GRASPcln, '.')
plt.scatter(aod550AERO, aod550GRASPcln, c=wndSpdCln, marker='.')
#plt.scatter(aod550AERO, aod550GRASPcln, c=h20Alb550Cln)
plt.plot(np.r_[0, maxAOD], np.r_[0, maxAOD], 'k')
plt.xlabel("AOD AERONET (550nm)")
plt.ylabel("AOD MODIS/GRASP (550nm)")
plt.xlim([0, maxAOD])
plt.ylim([0, maxAOD])
Rcoef = np.corrcoef(aod550AERO, aod550GRASPcln)[0,1]
RMSE = np.sqrt(np.mean((aod550AERO - aod550GRASPcln)**2))
textstr = 'N=%d\nR=%.3f\nRMS=%.3f\n'%(len(aod550AERO), Rcoef, RMSE)
plt.text(0.7*maxAOD, 0.03, textstr, fontsize=12)
plt.title(ttleNm)

# plot PSD
line = plt.figure()
N = len(dVdlnr)
[plt.semilogx(radii, sngdV, color=[i/N, 0.2, 1/(i+1)], linewidth=0.3) for i,sngdV in enumerate(dVdlnr)]
plt.semilogx(radii, np.mean(dVdlnr, axis=0), 'k')
plt.title(ttleNm + ' (blue=old, red=recent)')
plt.xlabel("radius (um)")
plt.ylabel("dV/dlnr (um^3/um^2)")

# plot RRI PDF
line = plt.figure()
plt.hist(rri550GRASP, 30)
plt.title(ttleNm)
plt.xlabel("Real Refractive Index (550nm)")
plt.ylabel("Frequency")
# plot IRI PDF
line = plt.figure()
plt.hist(iri550GRASP, 60)
plt.title(ttleNm)
plt.xlabel("Imaginary Refractive Index (550nm)")
plt.ylabel("Frequency")
# plot SPH PDF
line = plt.figure()
plt.hist(sph550GRASP, 30)
plt.title(ttleNm)
plt.xlabel("Spherical Fraction")
plt.ylabel("Frequency")

# plot H20 Albedo Time Series
line = plt.figure()
plt.plot(dateDays, h20AlbAll[:,1])
#plt.plot(dateDays, aod550GRASP/3)
plt.title(ttleNm)
plt.xlabel("Days since 2000")
plt.ylabel("Ocean Albedo Term")
# plot RRI Time Series
line = plt.figure()
plt.plot(dateDays, rri550GRASP)
plt.title(ttleNm)
plt.xlabel("Days since 2000")
plt.ylabel("Real Refractive Index (550nm)")
# plot Rv and sigma Time Series
fig = plt.figure()
plt.plot(dateDays, rv, dateDays, sigma)
fig.axes[0].set_yscale('log')
plt.title(ttleNm)
plt.xlabel("Days since 2000")

# plot H20 Albedo vs AOD
line = plt.figure()
plt.plot(h20AlbAll[:,1], aod550GRASP, '.')
plt.title(ttleNm)
plt.ylabel("GRASP AOD (550nm)")
plt.xlabel("Ocean Albedo Term")
# plot wind speed vs aod
line = plt.figure()
plt.plot(wndSpd, aod550GRASP, '.')
plt.plot(wndSpdCln, aod550AERO, '.')
plt.title(ttleNm)
plt.ylabel("GRASP AOD (550nm)")
plt.xlabel("GRASP Wind Speed (m/s)")
plt.legend(('GRASP','AERONET'))
# plot wind speed vs ocean albedo
line = plt.figure()
plt.scatter(wndSpd, h20AlbAll[:,1], c=aod550GRASP)
plt.title(ttleNm)
plt.ylabel("Ocean Albedo Term")
plt.xlabel("GRASP Wind Speed (m/s)")

# plot ocean albedo vs wavelength
line = plt.figure()
[plt.scatter(wvlngth, alb, color=[i/N, 0.2, 1/(i+1)], linewidth=0.3) for i,alb in enumerate(h20AlbAll)]
plt.plot(wvlngth, np.mean(h20AlbAll, axis=0), 'k')
plt.title(ttleNm)
plt.ylabel("Ocean Albedo Term")
plt.xlabel("Wavelength (um)")