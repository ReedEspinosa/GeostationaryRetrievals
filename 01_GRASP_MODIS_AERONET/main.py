#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import warnings
from read_MODAERO import modaeroDB
sys.path.append(os.path.join("..", "GRASP_PythonUtils"))
from runGRASP import graspDB

# NASA MacBook Air
#basePath = '/Users/wrespino/Synced'
#binPathGRASP = '/usr/local/bin/grasp'
#maxCPUs = 2 # max simaltanous instances of GRASP
# Personel MacBook Pro
#basePath = '/Volumes/MagHDD/LACO'
#binPathGRASP = '/usr/local/bin/grasp'
# tsingtau
#basePath = '/home/wrespino/data/synced/'
#binPathGRASP = '/home/wrespino/data/grasp_open/build/bin/grasp'
# Uranus
basePath = '/home/respinosa/ReedWorking/'
binPathGRASP = '/home/respinosa/ReedWorking/Local_Code_MacBook/grasp_open/build/bin/grasp'
maxCPUs = 12 # max simaltanous instances of GRASP

orbHght = 713 # km
#siteID 9=Wallops, 210=Nauru, 106=Ascension Island, 919=Fukue(Japan), 58=Cape Verde, Honolulu=97, Amr. Somoa=1185, Midway=278, Stennis=NOLA
#siteID = [210, 9, 106, 919, 58, 97, 1185, 278, 243] # can be int or list (not numpy array)
#siteID = [97] # can be int or list (not numpy array)
siteID = [210, 9, 106, 919, 58, 97, 1185, 278, 243, 179, 908, 809, 210, 1146, 103, 1101, 535, 739, 383, 840, 916, 475, 436, 667, 175] # can be int or list (not numpy array)
saveNPZ = False
incldAERO = False # Include AERONET AOD as input to retrieval

#savePath = os.path.join(basePath, 'Working/MODAERO_retrievalPickles/MODAERO_GRASP_MODISonly_Aqua_FullList.pkl')  # Save results here
savePath = os.path.join(basePath, 'Working/GRASP_testCases/testUranusFull.pkl') # Save results here
dirPath = os.path.join(basePath, 'Remote_Sensing_Projects/GRASP_MODIS/oceanV2/MOD*.out') 
#npzPath = os.path.join(basePath, 'Remote_Sensing_Projects/GRASP_MODIS/Ocean_MYD04_L2_C61_V2.out.npz') # Use AQUA Data
npzPath = os.path.join(basePath, 'Remote_Sensing_Projects/GRASP_MODIS/Ocean_MOD04_L2_C61_V2.out.npz') # Use Terra Data
pathYAML = os.path.join(basePath, 'Remote_Sensing_Analysis/GRASP_PythonUtils/settings_modis_2lgnrm_optimalNauru.yml') # YAML file to use
#pathYAML = os.path.join(basePath, 'Remote_Sensing_Analysis/GRASP_PythonUtils/settings_modis_2lgnrm_optimalNauru_noAeroTimeSmth.yml') # YAML file to use
dirGRASP = os.path.join(basePath, 'Remote_Sensing_Analysis/GRASP_testCases') # working path for GRASP

# --- CODE ---
if maxCPUs > 1: # need seperate directories to run parallel
    warnings.warn('System temp directories must be used in parallel mode. Ignoring dirGRASP setting...')
    dirGRASP = False

DB = modaeroDB(npzPath)

if DB.aod.shape[0]==0:
    DB.readDIR(dirPath)
    DB.sortData()
    if saveNPZ:
        DB.saveData(npzPath)
    else:
        warnings.warn('MODIS/AERONET data was loaded from text files but not saved.')

DB.groupData(siteID)

grspObjs = DB.graspPackData(pathYAML, incldAERO, orbHght, dirGRASP)
gDB = graspDB(grspObjs)
gDB.processData(maxCPUs, binPathGRASP, savePath)
