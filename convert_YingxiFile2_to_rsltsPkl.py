#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# This script converts Yingxi's collocated File 2 into a list of rslt dicts saved in pickle format 
# Designed to work only with Almucantar scans (no Hybrid scans allowed)

import numpy as np
import pickle
import re

csvFileIn = '/Users/wrespino/Synced/RST_CAN-GRASP/AERONET_collocation/ABI16_ALM/ABI16_ALM_5lines_V1.csv'


## Function and Variable Definitions ##
def poorMansAvg(trgt, val):
    if np.isnan(trgt) or np.isclose(trgt,val): return val
    return (trgt + val)/2

def AzmWvlVals2Ind(waveAzmVals, azms, wvls):
    waveAzmInd = np.empty(waveAzmVals.shape, dtype=np.uint16) # will only break if more than 64k columns in CSV file...
    for i,triplet in enumerate(waveAzmVals):
        azmInd, = np.isclose(azms, np.abs(triplet[1])).nonzero()
        wvlInd, = np.isclose(wvls, np.abs(triplet[2])).nonzero()
        assert len(azmInd)==1 and len(wvlInd)==1, 'Failure identifying index for azimuth %f and wavelength %f' % tuple(triplet[1:3])
        waveAzmInd[i, :] = [triplet[0], azmInd[0], wvlInd[0]]
    return waveAzmInd

azmPtrn = re.compile('([0-9]+)_Radiance_for_Azimuth_Angle_in_Degrees\(uW/cm\^2/sr/nm\)[ \-]+([\-0-9]+\.[0-9]+)') # TD: Do we need both "-" toward end of match?
scaPtrn = re.compile('([0-9]+)_Scattering_Angle_for_Nominal_Azimuth_Angle\(Degrees\)[ \-]+([\-0-9]+\.[0-9]+)')      # - Here too; and need to updated below if |angle|


csv2rsltScalar = { # TODO: Add AERONET angles and measurement time here
    'latitude':'Latitude(Degrees)',
    'longitude':'Longitude(Degrees)',
    'masl':'Elevation(m)',
    'sph':'Sphericity_Factor(%)',
    'AERO_siteName':'AERONET_Site',
    'AERO_siteID':'Instrument_Number',
    'AERO_sza':'Solar_Zenith_Angle(Degrees)',
    }

csv2rsltSpectral = { # CSV header format: "Lidar_Ratio[440nm]"
    'aod':'AOD_Extinction-Total\[([0-9]+)nm\]',
    'ssa':'Single_Scattering_Albedo\[([0-9]+)nm\]', # this will need to be repeated Nmode times
    'g':'Asymmetry_Factor-Total\[([0-9]+)nm\]',
    'n':'Refractive_Index-Real_Part\[([0-9]+)nm\]',
    'k':'Refractive_Index-Imaginary_Part\[([0-9]+)nm\]',
    'LidarRatio':'Lidar_Ratio\[([0-9]+)nm\]',
    'LidarDepol':'Depolarization_Ratio\[([0-9]+)nm\]',
    }
    
csv2rsltMode = {# CSV header format: "VolC-F"
    'vol':'VolC-([FC])',
    'rEffMode':'REff-([FC])',
    'rv':'VMR-([FC])',
    'sigma':'Std-([FC])',
    }

# TODO: add an additional csv2rsltPSD for radii dependent variables
# TODO: We need to include ABI angles too... Are they scalars? Either way, arrays needed in rslts dicts
# TODO: We are using 'meas_Isky' below but convention for things like angles seems to be 'AERO_sza'. Make consistent?

## Parse Header Row ##
with open(csvFileIn) as lines:
    colNames = lines.readline().rstrip().split(',')

# parse sky-scan data first (valid wavelengths are determined here)
waveAzmValsRad = []; waveAzmValsSca = []
for i,cn in enumerate(colNames):
    mtch = azmPtrn.match(cn)
    if mtch is not None:
        waveAzmValsRad.append([i, mtch.group(2), mtch.group(1)]) # this is [i, azm, wave]
    mtch = scaPtrn.match(cn)
    if mtch is not None:
        waveAzmValsSca.append([i, mtch.group(2), mtch.group(1)])

waveAzmValsRad = np.asarray(waveAzmValsRad, dtype=np.float64) # N_angles x 3; rows are: [i_column, azimuth_angle (deg), wavelength (nm)]
waveAzmValsSca = np.asarray(waveAzmValsSca, dtype=np.float32) # N_angles x 3; rows are: [i_column, azimuth_angle (deg), wavelength (nm)]

azms = np.unique(waveAzmValsRad[:,1])
assert np.all(azms==np.unique(waveAzmValsSca[:,1])), 'The relative azimuths for AERONET radiance and scattering angle data do not match!'
wvls = np.unique(waveAzmValsRad[:,2])
assert np.all(wvls==np.unique(waveAzmValsSca[:,2])), 'The wavelengths for AERONET radiance and scattering angle data do not match!'

waveAzmIndRad = AzmWvlVals2Ind(waveAzmValsRad, azms, wvls) # N_angles x 3; rows are: [i_column, azms_ind, wvls_ind]
waveAzmIndSca = AzmWvlVals2Ind(waveAzmValsSca, azms, wvls) # N_angles x 3; rows are: [i_column, azms_ind, wvls_ind]

# parse the header names for the rest of the columns
indDict = {}; indDictSpctrl = {}; indDictMode = {}
for i,cn in enumerate(colNames):
    for key,value in csv2rsltScalar.items():
        if cn==value: indDict[i]=key # indDict[i_column] = rsltKey; i.e., csv file column i corresponds to rslts[key] 
    for key,regEx in csv2rsltSpectral.items():
        regMtch = re.match(regEx, cn)
        if regMtch is not None:
            wvlsInd, = np.nonzero(wvls==int(regMtch.group(1)))
            assert len(wvlsInd)==1, 'Error matching %s with expected wavelengths!' % cn
            indDictSpctrl[i] = (key, wvlsInd[0]) # indDictSpctrl[i_column] = (rsltKey, wvls_ind)
    for key,regEx in csv2rsltMode.items():
        regMtch = re.match(regEx, cn)
        if regMtch is not None: 
            indDictMode[i] = (key, 0) if regMtch.group(1)=='F' else (key, 1) # indDictMode[i_column] = (rsltKey, mod_ind)

## Loop Over All Data Rows ##
csvData = np.genfromtxt(csvFileIn, delimiter=',', skip_header=1) # csvData[Nrows, NColumns]
rslts = [] # PERFORMANCE: This might be slow... Do we need to preallocate?
for pxRow in csvData: # loop over rows of CSV file; pxRow[NColumns]
    # TODO: Add progress indicator for when we process tens of thousands of rows
    rslts.append({
        'lambda':wvls,
        'fis':np.tile(azms, [len(wvls),1]).T,
        'meas_Isky':np.full((len(azms), len(wvls)), np.nan),
        'sca_ang':np.full((len(azms), len(wvls)), np.nan),
        })
    for triplet in waveAzmIndRad: # triplet[i_column, azms_ind, wvls_ind] for a given radiance value
        newVal = poorMansAvg(rslts[-1]['meas_Isky'][triplet[1], triplet[2]], pxRow[triplet[0]])
        rslts[-1]['meas_Isky'][triplet[1], triplet[2]] = newVal
    for triplet in waveAzmIndSca: # triplet[i_column, azms_ind, wvls_ind] for a given scattering angle value
        newVal = poorMansAvg(rslts[-1]['sca_ang'][triplet[1], triplet[2]], np.abs(pxRow[triplet[0]]))
        rslts[-1]['sca_ang'][triplet[1], triplet[2]] = newVal
    for csvInd,rsltKey in indDict.items():
        rslts[-1][rsltKey] = pxRow[csvInd]
    for csvInd, keyIndTuple in indDictSpctrl.items(): # keyIndTuple = (rsltKey, wvls_ind)
        if keyIndTuple[0] not in rslts[-1]: # allocate the keyIndTuple[0] array first...
            rslts[-1][keyIndTuple[0]] = np.full(len(wvls), np.nan, dtype=np.float64) # wavelengths may be unset so NANs needed
        rslts[-1][keyIndTuple[0]][keyIndTuple[1]] = pxRow[csvInd]
    for csvInd, keyIndTuple in indDictMode.items(): # keyIndTuple = (rsltKey, mode_ind)
        if keyIndTuple[0] not in rslts[-1]: # allocate the keyIndTuple[0] array first...
            rslts[-1][keyIndTuple[0]] = np.full(2, np.nan, dtype=np.float64) # wavelengths may be unset so NANs needed
        rslts[-1][keyIndTuple[0]][keyIndTuple[1]] = pxRow[csvInd]
    
# THEN ABI DT ocean [o] and land [l] RADIANCES (but not AOD retrievals.. Why not? L2 DT products useful before.)


pklFileOut = csvFileIn[:-3]+'pkl'
with open(pklFileOut, 'wb') as f:
    pickle.dump(rslts, f, pickle.HIGHEST_PROTOCOL)
print('Converted data saved to %s' % pklFileOut)








