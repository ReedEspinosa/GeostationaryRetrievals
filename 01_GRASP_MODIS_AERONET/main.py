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
binPathGRASP = '/usr/local/bin/grasp'
maxCPUs = 1 # max simaltanous instances of GRASP
# Personel MacBook Pro
#basePath = '/Volumes/MagHDD/LACO'
#binPathGRASP = '/usr/local/bin/grasp'
# tsingtau
#basePath = '/home/wrespino/data/synced/'
#binPathGRASP = '/home/wrespino/data/grasp_open/build/bin/grasp'
# Uranus
#basePath = '/home/respinosa/ReedWorking/'
#binPathGRASP = '/home/respinosa/ReedWorking/Local_Code_MacBook/grasp_open/build/bin/grasp'
#maxCPUs = 22 # max simaltanous instances of GRASP

#siteID 9=Wallops, 210=Nauru, 106=Ascension Island, 919=Fukue(Japan), 58=Cape Verde, Honolulu=97, Amr. Somoa=1185, Midway=278, Stennis=NOLA
# NOTE: these are above 600m (bit us on 946): 518,  268,  196,  514,  430, 1197,  946,   33,  806,  316, 666 
siteID_global = [210, 9, 106, 919, 58, 97, 1185, 278, 243, 908, 809, 210, 103, 1101, 739, 383, 916, 475, 436, 667] # globaly representative samplw
siteID_dust = [987, 352, 119, 392, 769, 793, 1031, 187, 962, 469, 471, 470, 543, 461, 467, 906, 1108, 1109, 1114, 507, 256] # DUST SITES
#siteID = np.unique(np.r_[siteID_global, siteID_dust])
#siteID = np.r_[9, 161, 649, 1126, 1149, 518, 510, 268, 196, 514, 59, 1, 430, 221, 1150, 1197, 946, 33, 1077, 1107, 337, 806, 1102, 176, 339, 989, 658, 532, 1156, 746, 961, 981, 748, 316, 666, 485, 940, 922, 285, 126] # LAND SITES
#siteID = siteID[0::4]
#siteID = np.r_[9]
siteID = np.r_[77,1,285,514,961,221,746]

orbHght = 713 # km
saveNPZ = False
incldAERO = False # Include AERONET AOD as input to retrieval

#npzPath = os.path.join(basePath, 'Remote_Sensing_Projects/GRASP_MODIS/Land_MYD04_L2_C61_V2c.out.npz') # Use AQUA land Data
#dirPath = os.path.join(basePath, 'Remote_Sensing_Projects/GRASP_MODIS/landV2/MYD*.outv2') 
#npzPath = os.path.join(basePath, 'Remote_Sensing_Projects/GRASP_MODIS/Land_MOD04_L2_C61_V2c.out.npz') # Use Terra land Data
#dirPath = os.path.join(basePath, 'Remote_Sensing_Projects/GRASP_MODIS/landV2/MOD*.outv2') 
#npzPath = os.path.join(basePath, 'Remote_Sensing_Projects/GRASP_MODIS/Ocean_MYD04_L2_C61_V2c.out.npz') # Use AQUA ocean Data
#dirPath = os.path.join(basePath, 'Remote_Sensing_Projects/GRASP_MODIS/oceanV2/MYD*.out') 
npzPath = os.path.join(basePath, 'Remote_Sensing_Projects/GRASP_MODIS/Ocean_MOD04_L2_C61_V2c.out.npz') # Use Terra ocean Data
dirPath = os.path.join(basePath, 'Remote_Sensing_Projects/GRASP_MODIS/oceanV2/MOD*.out') 
fileType = 'ocean' # 'land' or 'ocean'

savePath = os.path.join(basePath, 'Working/MODAERO_retrievalPickles/TerraLand_Select8Sites_YAMLa74e24c_maxBlue_3lgnrm_Nt300_V2c.pkl')  # Save results here
pathYAML = os.path.join(basePath, 'Local_Code_MacBook/MODAERO_Analysis/YAML_files/settings_modis_3lgnrm_optimalLand_maxBlue.yml') # YAML file to use
dirGRASP = False
#dirGRASP = os.path.join(basePath, 'Working/retreivalSandbox') # working path for GRASP, good for debuging
maxNtPerSeg=300 # how many pixels per grasp run? (DICOVER Test use 4 -> 45 grasp runs)


# --- CODE ---
if dirGRASP and maxCPUs > 1: # need seperate directories to run parallel
    warnings.warn('System temp directories must be used in parallel mode. Ignoring dirGRASP setting...')
    dirGRASP = False

DB = modaeroDB(None, fileType, maxNtPerSeg) if saveNPZ else modaeroDB(npzPath, fileType, maxNtPerSeg)

if DB.aod.shape[0]==0:
    DB.readDIR(dirPath)
    DB.sortData()
    if saveNPZ:
        DB.saveData(npzPath)
    else:
        warnings.warn('MODIS/AERONET data was loaded from text files but not saved.')

siteID_low = np.unique(DB.aero_loc[DB.aero_loc[:,1]<20,0])
siteID = np.unique(np.r_[siteID_global, siteID_dust, siteID_low])

DB.groupData(siteID)

grspObjs = DB.graspPackData(pathYAML, incldAERO, orbHght, dirGRASP)
sys.exit()
gDB = graspDB(grspObjs)

gDB.processData(maxCPUs, binPathGRASP, savePath)
