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

maxCPUs = 12 # max simaltanous instances of GRASP
orbHght = 713 # km
#siteID 9=Wallops, 210=Nauru, 106=Ascension Island, 919=Azores, 58=Cape Verde, Honolulu=97, Amr. Somoa=1185, Midway=278, San Cristobal=243
siteID = [210, 9, 106, 919, 58, 97, 1185, 278, 243] # can be int or list (not numpy array)
#siteID = [106, 210, 278] # can be int or list (not numpy array)
saveNPZ = False

savePath = pathPath = os.path.join(basePath, 'Working/MODAERO_retrievalPickles/MODAERO_GRASPref_updatedYAMLX.pkl')
dirPath = pathPath = os.path.join(basePath, 'Remote_Sensing_Projects/GRASP_MODIS/ocean/MOD*.out')
npzPath = pathPath = os.path.join(basePath, 'Remote_Sensing_Projects/GRASP_MODIS/Ocean_MOD04_L2_C61.out.npz')
pathYAML = pathPath = os.path.join(basePath, 'Remote_Sensing_Analysis/GRASP_PythonUtils/settings_modis_2lgnrm_optimalNauru.yml')
dirGRASP = os.path.dirname(pathYAML)

# --- CODE ---
if maxCPUs > 1: # need seperate directories to run parallel
    dirGRASP = False

DB = modaeroDB(npzPath)

if DB.aod.shape[0]==0:
    DB.readDIR(dirPath)
    if saveNPZ:
        DB.saveData(npzPath)

DB.sorted = False # we shoud be able to remove once the NPZ file is rebuilt
DB.groupData(siteID)

grspObjs = DB.graspPackData(pathYAML, orbHght, dirGRASP)
gDB = graspDB(grspObjs)
rslts = gDB.processData(maxCPUs, binPathGRASP)
dill.dump_session(savePath)
