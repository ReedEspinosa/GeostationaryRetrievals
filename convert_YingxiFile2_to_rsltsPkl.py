#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pickle
import re

csvFileIn = '/Users/wrespino/Synced/RST_CAN-GRASP/AERONET_collocation/ABI16_ALM/ABI16_ALM_5lines_V1.csv'

def poorMansAvg(trgt, val):
    if np.isnan(trgt) or np.isclose(trgt,val): return val
    return (trgt + val)/2

def AzmWvlVals2Ind(waveAzmVals, azms, wvls):
    waveAzmInd = np.empty(waveAzmVals.shape, dtype=np.int64)
    for i,pair in enumerate(waveAzmVals):
        azmInd, = np.isclose(azms, np.abs(pair[1])).nonzero()
        wvlInd, = np.isclose(wvls, np.abs(pair[2])).nonzero()
        assert len(azmInd)==1 and len(wvlInd)==1, 'Failure identifying index for azimuth %f and wavelength %f' % tuple(pair[1:3])
        waveAzmInd[i, :] = [pair[0], azmInd[0], wvlInd[0]]
    return waveAzmInd

azmPtrn = re.compile('([0-9]+)_Radiance_for_Azimuth_Angle_in_Degrees\(uW/cm\^2/sr/nm\)[ \-]+([\-0-9]+\.[0-9]+)')
scaPtrn = re.compile('([0-9]+)_Scattering_Angle_for_Nominal_Azimuth_Angle\(Degrees\)[ \-]+([\-0-9]+\.[0-9]+)')


with open(csvFileIn) as lines: 
    colNames = lines.readline().rstrip().split(',')

waveAzmValsRad = []; waveAzmValsSca = []
for i,cn in enumerate(colNames):
    mtch = azmPtrn.match(cn)
    if mtch is not None:
        waveAzmValsRad.append([i, mtch.group(2), mtch.group(1)]) # this is [i, azm, wave]
    mtch = scaPtrn.match(cn)
    if mtch is not None:
        waveAzmValsSca.append([i, mtch.group(2), mtch.group(1)])

waveAzmValsRad = np.asarray(waveAzmValsRad, dtype=np.float64)
waveAzmValsSca = np.asarray(waveAzmValsSca, dtype=np.float64)

azms = np.unique(waveAzmValsRad[:,1])
assert np.all(azms==np.unique(waveAzmValsSca[:,1])), 'The relative azimuths for AERONET radiance and scattering angle data do not match!'
wvls = np.unique(waveAzmValsRad[:,2])
assert np.all(wvls==np.unique(waveAzmValsSca[:,2])), 'The wavelengths for AERONET radiance and scattering angle data do not match!'

waveAzmIndRad = AzmWvlVals2Ind(waveAzmValsRad, azms, wvls)
waveAzmIndSca = AzmWvlVals2Ind(waveAzmValsSca, azms, wvls)

#indDict['latitude', 'longitude',...] then for each of these keys we add item in rslts[n][key]
# we need to code to run through colNames and swap all these to int key dicts (e.g., {3:'longitude', 4:'latitude',...})
fl2ToRsltScalar = {
    'longitude':'Latitude(Degrees)',
    'latitude':'Longitude(Degrees)',
    'masl':'Elevation(m)',
    'sph':'Sphericity_Factor(%)',
    }

fl2ToRsltSpectral = {
    'aod':'AOD_Extinction-Total',
    'ssa':'Single_Scattering_Albedo', # this will need to be repeated Nmode times
    'g':'Asymmetry_Factor-Total',
    'n':'Refractive_Index-Real_Part',
    'k':'Refractive_Index-Imaginary_Part',
    'LidarRatio':'Lidar_Ratio',
    'LidarDepol':'Depolarization_Ratio',
    }
    
fl2ToRsltMode = {
    'vol':'VolC-',
    'rEffMode':'REff-',
    'rv':'VMR-',
    'sigma':'Std-',
    }

# Need two more: one for anlges and one for radii


csvData = np.genfromtxt(csvFileIn, delimiter=',', skip_header=1)
rslts = []
for pxRow in csvData:
    rslts.append({
        'lambda':wvls,
        'fis':np.tile(azms, [len(wvls),1]).T,
        'meas_Isky':np.full((len(azms), len(wvls)), np.nan),
        'sca_ang':np.full((len(azms), len(wvls)), np.nan),
        })
    for pair in waveAzmIndRad: 
        newVal = poorMansAvg(rslts[-1]['meas_Isky'][pair[1], pair[2]], pxRow[pair[0]])
        rslts[-1]['meas_Isky'][pair[1], pair[2]] = newVal
    for pair in waveAzmIndSca: 
        newVal = poorMansAvg(rslts[-1]['sca_ang'][pair[1], pair[2]], np.abs(pxRow[pair[0]]))
        rslts[-1]['sca_ang'][pair[1], pair[2]] = newVal

# AERONET INVERSION PRODUCTS COME NEXT

# THEN ABI DT ocean [o] and land [l] RADIANCES (but not AOD retrievals)


pklFileOut = csvFileIn[:-3]+'pkl'
with open(pklFileOut, 'wb') as f:
    pickle.dump(rslts, f, pickle.HIGHEST_PROTOCOL)
print('Converted data saved to %s' % pklFileOut)








