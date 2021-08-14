#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import warnings
import numpy as np
from read_MODAERO import modaeroDB
sys.path.append(os.path.join("..", "GRASP_scripts")) # relative PATH to runGRASP class
from runGRASP import graspDB

# DISCOVER
#basePath = '/discover/nobackup/wrespino/synced'
#binPathGRASP = '/discover/nobackup/wrespino/grasp_open/build/bin/grasp'
#maxCPUs = 6 # max simaltanous instances of GRASP
# NASA MacBook Air
basePath = '/Users/wrespino/Synced'
# binPathGRASP = '/usr/local/bin/grasp'
binPathGRASP = '/Users/wrespino/Synced/Local_Code_MacBook/grasp_open/build/bin/grasp'
maxCPUs = 2 # max simaltanous instances of GRASP

#siteID 9=Wallops, 210=Nauru, 106=Ascension Island, 919=Fukue(Japan), 58=Cape Verde, Honolulu=97, Amr. Somoa=1185, Midway=278, Stennis=NOLA
# NOTE: these are above 600m (bit us on 946): 518,  268,  196,  514,  430, 1197,  946,   33,  806,  316, 666 
siteID_global = [210, 9, 106, 919, 58, 97, 1185, 278, 243, 908, 809, 210, 103, 1101, 739, 383, 916, 475, 436, 667] # globaly representative samplw
siteID_dust = [987, 352, 119, 392, 769, 793, 1031, 187, 962, 469, 471, 470, 543, 461, 467, 906, 1108, 1109, 1114, 507, 256] # DUST SITES
#siteID = np.unique(np.r_[siteID_global, siteID_dust])
#siteID = np.r_[9, 161, 649, 1126, 1149, 518, 510, 268, 196, 514, 59, 1, 430, 221, 1150, 1197, 946, 33, 1077, 1107, 337, 806, 1102, 176, 339, 989, 658, 532, 1156, 746, 961, 981, 748, 316, 666, 485, 940, 922, 285, 126] # LAND SITES
#siteID = siteID[0::4]
siteID = np.r_[9]
# siteID = np.r_[77,1,285,514,961,221,746]

orbHght = 713 # km
saveNPZ = False # False to load pre-generated NPZ file
incldAERO = False # Include AERONET AOD as input to retrieval

npzPath = os.path.join(basePath, 'CAN-GRASP/AERONET_collocation/ABI16/collocation_average_1920jj_oceanv10_20km.avg.npz')
dirPath = os.path.join(basePath, 'CAN-GRASP/AERONET_collocation/ABI16/collocation_average_1920jj_oceanv10_20km.avg') 
fileType = 'ocean' # 'land' or 'ocean'

savePath = os.path.join(basePath, 'Working/ABIAERO_retrievalPickles/Test_Wallops_ABI16_V0.pkl')  # Save results here
# pathYAML = os.path.join(basePath, 'Local_Code_MacBook/MODAERO_Analysis/YAML_files/settings_modis_3lgnrm_optimalLand_maxBlue.yml') # YAML file to use (V<1.0)
pathYAML = '/Users/wrespino/Synced/Local_Code_MacBook/MADCAP_Analysis/ACCP_ArchitectureAndCanonicalCases/settings_BCK_POLAR_3modes.yml' # YAML file to use
dirGRASP = False
#dirGRASP = os.path.join(basePath, 'Working/retreivalSandbox') # working path for GRASP, good for debuging
maxNtPerSeg=20 # how many pixels per grasp run? (DICOVER Test use 4 -> 45 grasp runs)
# maxNtPerSeg=200000 # HACK001 [2/2]

# --- CODE ---
if dirGRASP and maxCPUs > 1: # need seperate directories to run parallel
    warnings.warn('System temp directories must be used in parallel mode. Ignoring dirGRASP setting...')
    dirGRASP = False

if saveNPZ:
    DB = modaeroDB(None, fileType, maxNtPerSeg, instrument='ABI', verbose=1)
else: 
    DB = modaeroDB(npzPath, fileType, maxNtPerSeg, instrument='ABI', verbose=1)

if DB.aod.shape[0]==0:
    DB.readDIR(dirPath)
    DB.sortData()
    if saveNPZ:
        DB.saveData(npzPath)
    else:
        warnings.warn('Satellite and AERONET data was loaded from text files but not saved.')

# siteID_low = np.unique(DB.aero_loc[DB.aero_loc[:,1]<20,0])
# siteID = np.unique(np.r_[siteID_global, siteID_dust, siteID_low])

DB.groupData(siteID)

grspObjs = DB.graspPackData(pathYAML, incldAERO, orbHght, dirGRASP, releaseYAML=True)

gDB = graspDB(grspObjs)

gDB.processData(maxCPUs, binPathGRASP, savePath)
