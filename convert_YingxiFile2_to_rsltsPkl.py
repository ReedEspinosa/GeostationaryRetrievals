#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# This script converts Yingxi's collocated File 2 into a list of rslt dicts saved in pickle format 
# Designed to work only with Almucantar scans (no Hybrid scans allowed)
# rslts keys generally match convention for ABI; AERONET data keys appended with "_sky"
#   Although in some cases ABI does not match convention either (e.g., "_land" or "_ocean")

# TODO: Confirm ABI wavelengths with Yingxi
# TODO: We need to include ABI angles too... Are they scalars? Either way, arrays needed in rslts dicts
# TODO: Read in the Dark Target products
# TODO: add an additional csv2rsltPSD for radii dependent variables


import numpy as np
import pickle
import re
from os import path

# Path to CSV file for conversions
csvFileIn = '/Users/wrespino/Synced/RST_CAN-GRASP/AERONET_collocation/ABI16_ALM_TestFiles/ABI16_ALM_5lines_V1.csv'

## Patterns for AERONET measurement data ##
radPtrn_sky = re.compile('([0-9]+)_Radiance_for_Azimuth_Angle_in_Degrees\(uW/cm\^2/sr/nm\)[ \-]+([\-0-9]+\.[0-9]+)') # TD: Do we need both "-" toward end of match?
scaPtrn = re.compile('([0-9]+)_Scattering_Angle_for_Nominal_Azimuth_Angle\(Degrees\)[ \-]+([\-0-9]+\.[0-9]+)')      # - Here too; and need to updated below if |angle|

## Patterns for ABI measurement data ##
radPtrn_landABI = re.compile('refl([0-9]+)')
radPtrn_landStdABI = re.compile('refstdl([0-9]+)')
radPtrn_oceanABI = re.compile('refo([0-9]+)')
radPtrn_oceanStdABI = re.compile('refstdo([0-9]+)')
waveLand_ABI = [0.47, 0.64, 2.2] # NEEDS YINGXI's CONFIRMATION
waveOcean_ABI = [0.47, 0.55, 0.64, 0.86, 1.37, 1.6, 2.2] # NEEDS YINGXI's CONFIRMATION
atolWave = 0.01 # wavelength difference ≤10 nm will be considered negligible

## Mapping of CSV headers to rslt dict keys (csv2rsltsX have form {rsltKey:CSV_HEADER_STRING}) ##
csv2rsltStr = { # CSV headers found only once with for string values
    'datetime':'times', # 'datetime' is baked in below, careful with changes
    'AERO_siteName':'AERONET_Site',
    'ABI_fileName':'abifile' # 'ABI_fileName' is baked in below, careful with changes
    }

csv2rsltScalar = { # CSV headers found only once with for scalar (numeric) values
    'latitude':'Latitude(Degrees)',
    'longitude':'Longitude(Degrees)',
    'masl':'Elevation(m)',
    'sph':'Sphericity_Factor(%)',
    'AERO_siteID':'Instrument_Number',
    'sza_sky':'Solar_Zenith_Angle(Degrees)',
    }

csv2rsltSpectral = { # CSV header format: "Lidar_Ratio[440nm]"
    'aod':'AOD_Extinction-Total\[([0-9]+)nm\]', # TODO: Should aod from AERONET should be called aod_sky? Or aod_dt?
    'ssa':'Single_Scattering_Albedo\[([0-9]+)nm\]', # ssaMode would need to be repeated Nmode times...
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


## Function and Variable Definitions ##
def poorMansAvg(trgt, val):
    if np.isnan(trgt) or np.isclose(trgt,val): return val
    return (trgt + val)/2

def AzmWvlVals2Ind(waveAzmVals, wvls, azms=None):
    # waveAzmVals should be list with [i_column, wavelength (μm), scattering_angle (deg)]
    waveAzmInd = np.empty(waveAzmVals.shape, dtype=np.uint16) # will only break if more than 64k columns in CSV file...
    for i,triplet in enumerate(waveAzmVals):
        wvlInd, = np.isclose(wvls, triplet[1], atol=atolWave).nonzero()
        assert len(wvlInd)==1, 'Failure identifying wavelength index for wavelength %f' % (triplet[1])
        if azms is not None: 
            azmInd, = np.isclose(azms, np.abs(triplet[2])).nonzero() # relative azimuth so radiance at azm and -azm are (ideally) optically identical in 1D RT
            assert len(azmInd)==1, 'Failure identifying azimuth index for wavelength %f and azimuth %f' % tuple(triplet[1:3])
        waveAzmInd[i, :] = [triplet[0], wvlInd[0]] if azms is None else [triplet[0], wvlInd[0], azmInd[0]]
    return waveAzmInd


## Parse Header Row ##
with open(csvFileIn) as lines:
    colNames = lines.readline().rstrip().split(',')

# parse sky-scan data first (valid wavelengths are determined here)
waveAzmValsRad = []; waveAzmValsSca = []
waveValsRadLandABI = []; waveValsStdLandABI = []; waveValsRadOceanABI = []; waveValsStdOceanABI = [];
for i,cn in enumerate(colNames):
    mtch = radPtrn_sky.match(cn) # AERONET sky-scan radiances
    if mtch is not None:
        waveAzmValsRad.append([i, mtch.group(1), mtch.group(2)]) # this is [i, wave, azm]
    mtch = scaPtrn.match(cn) # AERONET sky-scan scattering angles
    if mtch is not None:
        waveAzmValsSca.append([i, mtch.group(1), mtch.group(2)])
    mtch = radPtrn_landABI.match(cn) # ABI land radiances
    if mtch is not None:
        waveValsRadLandABI.append([i, mtch.group(1)])
    mtch = radPtrn_landStdABI.match(cn) # ABI land radiances STD
    if mtch is not None:
        waveValsStdLandABI.append([i, mtch.group(1)])
    mtch = radPtrn_oceanABI.match(cn) # ABI ocean radiances
    if mtch is not None:
        waveValsRadOceanABI.append([i, mtch.group(1)])
    mtch = radPtrn_oceanStdABI.match(cn) # ABI ocean radiances STD
    if mtch is not None:
        waveValsStdOceanABI.append([i, mtch.group(1)])

waveAzmValsRad = np.asarray(waveAzmValsRad, dtype=np.float64) # N_angles x 3; rows are: [i_column, wavelength (nm), azimuth_angle (deg)]
waveAzmValsRad[:,1] = waveAzmValsRad[:,1]/1e3 # Convert nm to μm; [i_column, wavelength (μm), azimuth_angle (deg)]
waveAzmValsSca = np.asarray(waveAzmValsSca, dtype=np.float32) # N_angles x 3; rows are: [i_column, wavelength (nm), scattering_angle (deg)]
waveAzmValsSca[:,1] = waveAzmValsSca[:,1]/1e3 # Convert nm to μm; [i_column, wavelength (μm), scattering_angle (deg)]

# define the set of valid wavelengths and sky-scan azimuth angles
azms_sky = np.unique(waveAzmValsRad[:,2])
assert np.all(azms_sky==np.unique(waveAzmValsSca[:,2])), 'The relative azimuths for AERONET radiance and scattering angle data do not match!'
wvls = np.unique(waveAzmValsRad[:,1])
assert np.all(wvls==np.unique(waveAzmValsSca[:,1])), 'The wavelengths for AERONET radiance and scattering angle data do not match!'
for wvl in (waveLand_ABI + waveOcean_ABI): # Note this is list concatenation, not addition
    if not np.isclose(wvl, wvls, atol=atolWave).any(): wvls = np.r_[wvls, wvl]
wvls = np.sort(wvls)

# match CSV columns with wavelengths and angles
waveAzmIndRad = AzmWvlVals2Ind(waveAzmValsRad, wvls, azms_sky) # N_angles*N_wvls x 3; rows are: [i_column, wvls_ind, azms_ind]
waveAzmIndSca = AzmWvlVals2Ind(waveAzmValsSca, wvls, azms_sky) # N_angles*N_wvls x 3; rows are: [i_column, wvls_ind, azms_ind]
waveAzmIndRadLandABI = AzmWvlVals2Ind(waveValsRadLandABI, wvls) # N_wvls x 2; rows are: [i_column, wvls_ind]
waveAzmIndStdLandABI = AzmWvlVals2Ind(waveValsStdLandABI, wvls) # N_wvls x 2; rows are: [i_column, wvls_ind]
waveAzmIndRadOceanABI = AzmWvlVals2Ind(waveValsRadOceanABI, wvls) # N_wvls x 2; rows are: [i_column, wvls_ind]
waveAzmIndStdOceanABI = AzmWvlVals2Ind(waveValsStdOceanABI, wvls) # N_wvls x 2; rows are: [i_column, wvls_ind]

# parse the header names for the rest of the columns
indDictScalar = {}; indDictStr = {}; indDictSpctrl = {}; indDictMode = {}
for i,cn in enumerate(colNames):
    for key,value in csv2rsltScalar.items():
        if cn==value: indDictScalar[i]=key # indDictScalar[i_column] = rsltKey; i.e., csv file column i corresponds to rslts[key] 
    for key,value in csv2rsltStr.items():
        if cn==value: indDictStr[i]=key # indDictStr[i_column] = rsltKey; i.e., csv file column i corresponds to rslts[key] 
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

# find location of misc. special format columns
spcInd = dict()
spcInd['time'] = np.nonzero(colNames[0]=='times')[0]
spcInd['siteName'] = np.nonzero(colNames[0]=='AERONET_Site')[0]
spcInd['abifile'] = np.nonzero(colNames[0]=='abifile')[0]
for key,vals in spcInd.items():
    assert len(vals)==1, 'Error deteremining %s column while parsing CSV header row!' % key 

## Loop Over All Data Rows ##
#  parse both string and numerical data and create list of rslt dicts
strColumnInds = tuple(colInd for colInd in indDictStr.keys()) # Python 3.7+ => dict iteration order guaranteed to be in order of insertion
csvStrData = np.genfromtxt(csvFileIn, delimiter=',', skip_header=1, usecols=strColumnInds, dtype=str) # csvStrData[Nrows, len(csv2rsltStr)]
csvNumData = np.genfromtxt(csvFileIn, delimiter=',', skip_header=1) # csvNumData[Nrows, NColumnsTotal]
rslts = []
dispString = '%07d/%07d rows processed' # len is 30 characters
print(dispString % (0,len(csvNumData)))
for i,(numRow, strRow) in enumerate(zip(csvNumData, csvStrData)): # loop over rows of CSV file; numRow[NColumns]
    if np.mod(i,2)==0: print('\b'*30 + (dispString % (i, len(csvNumData)))) 
    rslts.append({
        'lambda':wvls,
        'fis_sky':np.tile(azms_sky, [len(wvls),1]).T,
        'meas_I_sky':np.full((len(azms_sky), len(wvls)), np.nan),
        'meas_I_land':np.full((len(wvls)), np.nan),
        'meas_I_ocean':np.full((len(wvls)), np.nan),
        'meas_I_landStd':np.full((len(wvls)), np.nan),
        'meas_I_oceanStd':np.full((len(wvls)), np.nan),
        'sca_ang_sky':np.full((len(azms_sky), len(wvls)), np.nan),
        })
    for jStr, key in enumerate(indDictStr.values()): # jStr is column in csvStrData of rslt[key]; Python 3.7+ => dict iteration order guaranteed to be in order of insertion
        if key == 'datetime':
            rslts[-1]['datetime'] = datetime.strptime(strRow[jStr], '%Y-%m-%d %H:%M:%S')
        elif key == 'ABI_fileName':
            rslts[-1]['ABI_fileName'] = path.basename(strRow[jStr])
        else:
            rslts[-1][key] = strRow[jStr]
    datetime.strptime(ts, '%Y-%m-%d %H:%M:%S')
    for triplet in waveAzmIndRad: # triplet[i_column, wvls_ind, azms_ind] for a given radiance value
        newVal = poorMansAvg(rslts[-1]['meas_I_sky'][triplet[2], triplet[1]], numRow[triplet[0]])
        rslts[-1]['meas_I_sky'][triplet[2], triplet[1]] = newVal
    for triplet in waveAzmIndSca: # triplet[i_column, wvls_ind, azms_ind] for a given scattering angle value
        newVal = poorMansAvg(rslts[-1]['sca_ang_sky'][triplet[2], triplet[1]], np.abs(numRow[triplet[0]])) # absolute values to eliminate negative scattering angles
        rslts[-1]['sca_ang_sky'][triplet[2], triplet[1]] = newVal
    for pair, pairStd in zip(waveAzmIndRadOceanABI, waveAzmIndStdOceanABI): # pair[i_column, wvls_ind] for a given ABI OCEAN radiance value
        rslts[-1]['meas_I_ocean'][pair[1]] = numRow[pair[0]]
        rslts[-1]['meas_I_oceanStd'][pairStd[1]] = numRow[pairStd[0]]
    for pair, pairStd in zip(waveAzmIndRadLandABI, waveAzmIndStdLandABI): # pair[i_column, wvls_ind] for a given ABI LAND radiance value
        rslts[-1]['meas_I_land'][pair[1]] = numRow[pair[0]]
        rslts[-1]['meas_I_landStd'][pairStd[1]] = numRow[pairStd[0]]
    for csvInd,rsltKey in indDictScalar.items():
        rslts[-1][rsltKey] = numRow[csvInd]
    for csvInd, keyIndTuple in indDictSpctrl.items(): # keyIndTuple = (rsltKey, wvls_ind)
        if keyIndTuple[0] not in rslts[-1]: # allocate the keyIndTuple[0] array first...
            rslts[-1][keyIndTuple[0]] = np.full(len(wvls), np.nan, dtype=np.float64) # wavelengths may be unset so NANs needed
        rslts[-1][keyIndTuple[0]][keyIndTuple[1]] = numRow[csvInd]
    for csvInd, keyIndTuple in indDictMode.items(): # keyIndTuple = (rsltKey, mode_ind)
        if keyIndTuple[0] not in rslts[-1]: # allocate the keyIndTuple[0] array first...
            rslts[-1][keyIndTuple[0]] = np.full(2, np.nan, dtype=np.float64) # wavelengths may be unset so NANs needed
        rslts[-1][keyIndTuple[0]][keyIndTuple[1]] = numRow[csvInd]

pklFileOut = csvFileIn[:-3]+'pkl'
with open(pklFileOut, 'wb') as f:
    pickle.dump(rslts, f, pickle.HIGHEST_PROTOCOL)
print('Converted data saved to %s' % pklFileOut)








