#!/usr/bin/env python3

import matplotlib.pyplot as plt
import datetime as dt
import numpy as np
import dill
import os
import sys
sys.path.append(os.path.join("..", "GRASP_PythonUtils"))
from miscFunctions import angstrmIntrp, angstrm
from scipy.stats import gaussian_kde

#dillFN = '/Users/wrespino/Synced/Working/MODAERO_retrievalPickles/MODAERO_GRASP_MODISonly_Terra.pkl'
#dillFN = '/Users/wrespino/Synced/Working/MODAERO_retrievalPickles/MODAERO_GRASP_AERONETincld_surfLoose.pkl'
#dillFN = '/Users/wrespino/Synced/Working/MODAERO_retrievalPickles/MODAERO_GRASP_MODISonly.pkl'
dillFN = '/Users/wrespino/Synced/Working/MODAERO_retrievalPickles/MODAERO_GRASP_MODISonly_Aqua.pkl'
#dillFN = '/Users/wrespino/Synced/Working/MODAERO_retrievalPickles/MODAERO_GRASP_MODISonly_singlePix.pkl'
#dillFN = '/Users/wrespino/Synced/Working/MODAERO_retrievalPickles/MODAERO_GRASP_MODISonly_noAeroTimeSmooth.pkl'

ttlTxt = 'MODIS AQUA' 

maxAOD = 2

siteID_plot = False # False for all sites (doesn't apply to AOD)
ttleNm = ''
lmbdInd = 2 # [0.469, 0.555, 0.645, 0.8585, 1.24, 1.64, 2.13]
#biasCrct = [1.85, 1.32, 1, 0.62]
biasCrct = [1, 1, 1, 1]
plotDT = False # AE IS ALWAYS GRASP (despite label)
AODonly = True
difVSclrPlot = False;
logLogAOD = False # plot AOD on log-log plot
plotMode = 1

dill.load_session(dillFN)
radii = rslts[0][0]['r'][0] # assume the same radii are used in every retrieval; analysis:ignore (it was loaded)
wvlngth = rslts[0][0]['lambda'] # assume the same lambda used in every retrieval; analysis:ignore (it was loaded)
Nmodes = rslts[0][0]['r'].shape[0] # analysis:ignore (it was loaded)
Npsd = radii.shape[0]
Nlmbd = wvlngth.shape[0] 
aodGRASP = np.array([])
aeGRASP = np.array([])
rriGRASP = np.array([])
iriGRASP = np.array([])
h20Alb550 = np.array([])
h20AlbAll = np.array([]).reshape(0, Nlmbd) # index 1 will be 550nm
wndSpd = np.array([])
dateDays = np.array([]) 
dVdlnr = np.array([]).reshape(0, Npsd)
volFMF = np.array([])
volTot = np.array([])
rv = np.array([]).reshape(0, Nmodes)
sigma = np.array([]).reshape(0, Nmodes)
sphGRASP = np.array([]).reshape(0, Nmodes)
# MODIS GRASP AODs
for rsltRun in rslts: # analysis:ignore (it was loaded)
    for rslt in rsltRun:
        aodGRASP = np.append(aodGRASP, rslt['aod'][lmbdInd])
        aeGRASP = np.append(aeGRASP, angstrm(wvlngth[[1,3]], rslt['aod'][[1,3]]))
        rriGRASP = np.append(rriGRASP, rslt['n'][plotMode,lmbdInd]) # BUG: these two are only fine mode values
        iriGRASP = np.append(iriGRASP, rslt['k'][plotMode,lmbdInd]) #      fix needed here and in reading of GRASP output
        sphGRASP = np.vstack([sphGRASP, rslt['sph']])
        dVdlnr = np.vstack([dVdlnr, rslt['vol']@rslt['dVdlnr']]) # scale to um^3/um^2 & sum multiple modes
        volFMF = np.append(volFMF, rslt['vol'][0]/(np.sum(rslt['vol']))) # technically this is "first mode fraction"
        volTot = np.append(volTot, np.sum(rslt['vol']))
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
aeAERO = np.array([])
aodDT = np.array([])
windNCEP = np.array([])
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
        grnInd = np.argmin(np.abs(seg.AERO_LAMDA-0.554))
        nirInd = np.argmin(np.abs(seg.AERO_LAMDA-0.87))
        aeInd = [grnInd, nirInd]
        aeAERO = np.append(aeAERO, angstrm(seg.AERO_LAMDA[aeInd], seg.aod[nowInd,aeInd]))
        siteID = np.append(siteID, seg.aero_loc[0])
        windNCEP = np.append(windNCEP, seg.metaData[nowInd,1])
        aodDT = np.append(aodDT, seg.modDT_aod[nowInd, lmbdInd])
assert (aodAERO.shape[0] == aodGRASP.shape[0]), 'rslts and DB.siteSegment contained a different number of pixels!'
aodDT[aodDT<0.00001] = 0.00001

# compare with AERONET
aodAERO = np.array(aodAERO)
if plotDT:
    aodGRASP = np.array(aodDT)
    rtrvlStr = 'Dark Target'
else:
    aodGRASP = np.array(aodGRASP)
    rtrvlStr = 'GRASP'
aodRtrvdCln = aodGRASP[~np.isnan(aodAERO)]/biasCrct[lmbdInd]
#aodRtrvdCln[aodRtrvdCln < 0.009] = 0.009 # HACK
aeNotNan = ~np.isnan(aeAERO)
aeRtrvdCln = aeGRASP[aeNotNan]
aeRtrvdCln = aeRtrvdCln[aodAERO[aeNotNan]>0.05]

# THESE ALL COLOR THE POINTS BY SOME AUX VARIABLE, TO COLOR BY DENSITY SEE BELOW
#clrVar = wndSpd[~np.isnan(aodAERO)] # WIND SPEED
#clrVar = h20AlbAll[~np.isnan(aodAERO), 3] # OCEAN ALBEDO, last ind. wavelength
#clrVar = volFMF[~np.isnan(aodAERO)] # FIRST MODE FRACTION (by volume)
#clrVar = sphGRASP[~np.isnan(aodAERO),1] # SPH FRACTION (mode select left)
#clrVar = rriGRASP[~np.isnan(aodAERO)] # RRI (Mode selected above)
#clrVar = iriGRASP[~np.isnan(aodAERO)] # IRI (Mode selected above)
#clrVar = rv[~np.isnan(aodAERO),1] # rv (mode select left)
#clrVar = sigma[~np.isnan(aodAERO),0] # sigma (mode select left)
#clrVar[clrVar < np.percentile(clrVar,2)] = np.percentile(clrVar,2)
#clrVar[clrVar > np.percentile(clrVar,98)] = np.percentile(clrVar,98)

aeAERO = aeAERO[aeNotNan]
aeAERO = aeAERO[aodAERO[aeNotNan]>0.05]
windNCEPcln = windNCEP[~np.isnan(aodAERO)]
wndSpdcln = wndSpd[~np.isnan(aodAERO)]
aodAERO = aodAERO[~np.isnan(aodAERO)]

#COLOR BY DENSITY
xy = np.log(np.vstack([aodAERO,aodRtrvdCln])) if logLogAOD else np.vstack([aodAERO,aodRtrvdCln])
clrVar = gaussian_kde(xy)(xy)

line = plt.figure()
ax = plt.scatter(aodAERO, aodRtrvdCln, c=clrVar, marker='.')
pax = plt.plot(np.r_[0, maxAOD], np.r_[0, maxAOD], 'k')
plt.xlabel("AOD AERONET" + lmbdStr)
plt.ylabel("AOD " + rtrvlStr + " " + lmbdStr)
plt.xlim([0.01, maxAOD])
plt.ylim([0.01, maxAOD])
Rcoef = np.corrcoef(aodAERO, aodRtrvdCln)[0,1]
RMSE = np.sqrt(np.mean((aodAERO - aodRtrvdCln)**2))
#fTau = np.mean(aodRtrvdCln/aodAERO)
#fTau = 2*np.mean((aodRtrvdCln-aodAERO)/(aodRtrvdCln+aodAERO))
#textstr = 'N=%d\nR=%.3f\nRMS=%.3f\nbias=%.3f%%\n'%(len(aodAERO), Rcoef, RMSE, 100*fTau)
fTau = np.mean((aodRtrvdCln-aodAERO))
textstr = 'N=%d\nR=%.3f\nRMS=%.3f\nbias=%.3f\n'%(len(aodAERO), Rcoef, RMSE, fTau)
if logLogAOD:
    plt.text(0.2*maxAOD, 0.015, textstr, fontsize=12)
    plt.yscale('log')
    plt.xscale('log')
else:
    plt.text(0.7*maxAOD, 0.03, textstr, fontsize=12)
plt.title(ttlTxt)

if difVSclrPlot:
    line = plt.figure()
    ax = plt.scatter(aodRtrvdCln/aodAERO, clrVar, marker='.')
    plt.xlabel("AOD: MODIS/GRASP - AERONET" + lmbdStr)
    plt.ylabel("Color Parameter")
    plt.xlim([0.3, 3])
    plt.title('All Sites')

# AOD diff plot
minBin = 5e-3
maxBin = 2
FS=14
binEdge = np.r_[np.logspace(np.log10(minBin),np.log10(maxBin), 13)];
binMid = np.exp(np.log(binEdge[1:]*binEdge[:-1])/2)
aodDif = aodRtrvdCln-aodAERO
aodRng = np.zeros([binEdge.shape[0]-1,3]) # lower (16%), mean, upper (84%)
for i in range(binEdge.shape[0]-1):
    aodDifNow = aodDif[np.logical_and(aodAERO > binEdge[i], aodAERO <= binEdge[i+1])]
    if aodDifNow.shape[0]==0:
        aodRng[i,:] = np.nan
    else:
        aodRng[i,0] = np.percentile(aodDifNow, 16)
        aodRng[i,1] = np.mean(aodDifNow)
        aodRng[i,2] = np.percentile(aodDifNow, 84)
plt.figure()
plt.xscale('log')
plt.plot([1e-5,maxAOD], [0,0], 'k')
x = np.logspace(-4, 1, 1000)
y = 0.03+0.1*x
#y = 0.05+0.15*x
plt.plot(x, y, '--', color=[0.5,0.5,0.5])
plt.plot(x, -y, '--', color=[0.5,0.5,0.5])
errBnds = np.abs(aodRng[:,1].reshape(-1,1)-aodRng[:,[0,2]]).T
plt.errorbar(binMid, aodRng[:,1], errBnds, ecolor='r', color='r', marker='s', linstyle=None)
plt.xlim([minBin, maxAOD])
plt.ylim([-0.55, 0.55])
plt.xlabel('AERONET AOD' + lmbdStr)
plt.ylabel('GRASP-AERONET AOD' + lmbdStr)
inEE = np.sum(np.abs(aodDif) < 0.03+0.1*aodAERO)/aodAERO.shape[0]
in2EE = np.sum(np.abs(aodDif) < 2*(0.03+0.1*aodAERO))/aodAERO.shape[0]
plt.text(0.01, -0.3, 'N=%d\nwithin 1xEE: %4.1f%%\nwithin 2xEE: %4.1f%%' % (aodAERO.shape[0],100*inEE,100*in2EE), FontSize=FS)
plt.text(0.02, 0.06, 'EE=0.03+0.1τ', FontSize=FS)


if AODonly: sys.exit()

# AE Plots
xy = np.vstack([aeAERO,aeRtrvdCln])
clrVar = gaussian_kde(xy)(xy)

line = plt.figure()
ax = plt.scatter(aeAERO, aeRtrvdCln, c=clrVar, marker='.')
pax = plt.plot(np.r_[0, 3], np.r_[0, 3], 'k')
plt.xlabel("AE AERONET")
plt.ylabel("AE " + rtrvlStr + " " + lmbdStr)
plt.xlim([-0.35, 2.8])
plt.ylim([-0.35, 2.8])
Rcoef = np.corrcoef(aeAERO, aeRtrvdCln)[0,1]
RMSE = np.sqrt(np.mean((aeAERO - aeRtrvdCln)**2))
#fTau = np.mean(aeRtrvdCln/aeAERO)
fTau = 2*np.mean((aeRtrvdCln-aeAERO)/(aeRtrvdCln+aeAERO))
textstr = 'N=%d\nR=%.3f\nRMS=%.3f\nbiasτ=%.3f%%\n'%(len(aodAERO), Rcoef, RMSE, 100*fTau)
plt.text(2, 0, textstr, fontsize=12)

# WIND STUFF
plt.figure()
#y = 2*(aodRtrvdCln-aodAERO)/(aodRtrvdCln+aodAERO)
#x = 2*(windNCEPcln-wndSpdcln)/(windNCEPcln+wndSpdcln)
#xy = np.vstack([x,y])
xy = np.vstack([windNCEP,wndSpd])
clrVar = gaussian_kde(xy)(xy)
ax = plt.scatter(windNCEP, wndSpd, c=clrVar, marker='.')
plt.xlabel('NCEP')
plt.ylabel('GRASP')
plt.title('Wind Speed (m/s)')
#ax = plt.scatter(x, y, c=clrVar, marker='.')
#plt.ylim([-2,2])
#plt.xlim([-2,2])
#plt.xlabel('2(vGRASP-vNCEP)/(vGRASP+vNCEP)')
#plt.ylabel('2(τGRASP-τAERONET)/(τGRASP+τAERONET)')
#plt.title('Wind Speed vs AOD Biases ' + lmbdStr)
#Rcoef = np.corrcoef(xy[0,:], xy[1,:])[0,1]
#RMSE = np.sqrt(np.mean((xy[0,:] - xy[1,:])**2))
#bias = np.mean(xy[1,:]-xy[0,:])
#textstr = 'N=%d\nR=%.3f\nRMS=%.3f\nbiasτ=%.3f%%\n'%(len(aodAERO), Rcoef, RMSE, bias)
#plt.text(6, 3, textstr, fontsize=12)



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
lowSphInd = sphGRASP[:,1] < 70 # BUG: this ignores plotting site (fine if slice(None))
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(7, 5))
ax[0].hist(rriGRASP[pltInd], 30)
ax[0].set_xlabel("RRI" + lmbdStr)
ax[0].set_ylabel("Frequency (a.u.)")
ax[0].hist(rriGRASP[lowSphInd], 30)
# plot IRI PDF
ax[1].hist(iriGRASP[pltInd], 60)
ax[1].set_title(ttleNm)
ax[1].set_xlabel("IRI" + lmbdStr)
ax[1].hist(iriGRASP[lowSphInd], 60)

# plot rv PDF
fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(10, 5))
ax[0].hist(rv[pltInd,plotMode], 30)
ax[0].set_xlabel("rv")
ax[0].set_ylabel("Frequency (a.u.)")
ax[0].hist(rv[lowSphInd,plotMode], 30)
# plot sigma PDF
ax[1].hist(sigma[pltInd,plotMode], 60)
ax[1].set_title(ttleNm)
ax[1].set_xlabel("σ")
ax[1].hist(sigma[lowSphInd,plotMode], 60)
# plot SPH PDF
ax[2].hist(sphGRASP[pltInd,1], 30) # ONLY COARSE MODE
ax[2].set_xlabel("Spherical Fraction")
ax[2].hist(sphGRASP[lowSphInd,1], 30)

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