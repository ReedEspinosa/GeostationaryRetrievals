#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import dill
from read_MODAERO import modaeroDB
sys.path.append(os.path.join("..", "GRASP_PythonUtils"))
from runGRASP import graspDB

# NASA MacBook Air
#basePath = '/Users/wrespino/Synced'
#binPathGRASP = '/usr/local/bin/grasp'
# Personel MacBook Pro
#basePath = '/Volumes/MagHDD/LACO'
#binPathGRASP = '/usr/local/bin/grasp'
# tsingtau
#basePath = '/home/wrespino/data/synced/'
#binPathGRASP = '/home/wrespino/data/grasp_open/build/bin/grasp'
# Uranus
basePath = '/home/respinosa/ReedWorking/'
binPathGRASP = '/home/respinosa/ReedWorking/Local_Code_MacBook/grasp_open/build/bin/grasp'

maxCPUs = 12 # (12) max simaltanous instances of GRASP
orbHght = 713 # km
#siteID 9=Wallops, 210=Nauru, 106=Ascension Island, 919=Azores, 58=Cape Verde, Honolulu=97, Amr. Somoa=1185, Midway=278, San Cristobal=243
siteID = [210, 9, 106, 919, 58, 97, 1185, 278, 243] # can be int or list (not numpy array)
#siteID = [106, 210, 278] # can be int or list (not numpy array)
saveNPZ = False
incldAERO = False # Include AERONET AOD as input to retrieval

#savePath = os.path.join(basePath, 'Working/MODAERO_retrievalPickles/MODAERO_GRASP_AERONETincld.pkl')
#savePath = os.path.join(basePath, 'Working/MODAERO_retrievalPickles/MODAERO_GRASP_MODISonly_SinglePix.pkl')
#savePath = os.path.join(basePath, 'Working/MODAERO_retrievalPickles/MODAERO_GRASP_MODISonly_Aqua.pkl')
savePath = os.path.join(basePath, 'Working/MODAERO_retrievalPickles/MODAERO_GRASP_MODISonly_Terra.pkl')
dirPath = os.path.join(basePath, 'Remote_Sensing_Projects/GRASP_MODIS/oceanV2/MOD*.out')
#npzPath = os.path.join(basePath, 'Remote_Sensing_Projects/GRASP_MODIS/Ocean_MYD04_L2_C61_V2.out.npz')
npzPath = os.path.join(basePath, 'Remote_Sensing_Projects/GRASP_MODIS/Ocean_MOD04_L2_C61_V2.out.npz')
#pathYAML = os.path.join(basePath, 'Remote_Sensing_Analysis/GRASP_PythonUtils/settings_modisAERONET_2lgnrm_7lamda_loose.yml')
pathYAML = os.path.join(basePath, 'Remote_Sensing_Analysis/GRASP_PythonUtils/settings_modis_2lgnrm_optimalNauru.yml')
#pathYAML = os.path.join(basePath, 'Remote_Sensing_Analysis/GRASP_PythonUtils/settings_modis_2lgnrm_optimalNauru_noAeroTimeSmth.yml')
dirGRASP = os.path.dirname(pathYAML)

# --- CODE ---
if maxCPUs > 1: # need seperate directories to run parallel
    dirGRASP = False

DB = modaeroDB(npzPath)

if DB.aod.shape[0]==0:
    DB.readDIR(dirPath)
    DB.sortData()
    if saveNPZ:
        DB.saveData(npzPath)

DB.groupData(siteID)

grspObjs = DB.graspPackData(pathYAML, incldAERO, orbHght, dirGRASP)
gDB = graspDB(grspObjs)
rslts = gDB.processData(maxCPUs, binPathGRASP)
dill.dump_session(savePath)
