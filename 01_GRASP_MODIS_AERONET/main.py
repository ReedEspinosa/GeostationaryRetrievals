#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import dill
from read_MODAERO import modaeroDB
sys.path.append(os.path.join("..", "GRASP_PythonUtils"))
from runGRASP import graspDB

# NASA MacBook Air
basePath = '/Users/wrespino/Synced'
binPathGRASP = '/usr/local/bin/grasp'
# Personel MacBook Pro
#basePath = '/Volumes/MagHDD/LACO'
#binPathGRASP = '/usr/local/bin/grasp'
# tsingtau
#basePath = '/home/wrespino/data/synced/'
#binPathGRASP = '/home/wrespino/data/grasp_open/build/bin/grasp'
# Uranus
#basePath = '/home/respinosa/ReedWorking/'
#binPathGRASP = '/home/respinosa/grasp_open/build/bin/grasp'

maxCPUs = 2
dirGRASP = os.path.dirname(pathYAML)
orbHght = 713 #km
siteID = 9 #Nauru
segIndGRASP = slice(None)

savePath = pathPath = os.path.join(basePath, 'Working/MODAERO_16binLoose.pkl')
dirPath = pathPath = os.path.join(basePath, 'Remote_Sensing_Projects/GRASP_MODIS/ocean/MOD*.out')
npzPath = pathPath = os.path.join(basePath, 'Remote_Sensing_Projects/GRASP_MODIS/Ocean_MOD04_L2_C61.out.npz')
pathYAML = pathPath = os.path.join(basePath, 'Remote_Sensing_Projects/GRASP_MODIS/settings_modis.yml')

# --- CODE ---
DB = modaeroDB(npzPath)

if DB.aod.shape[0]==0:
    DB.readDIR(dirPath)
#    DB.saveData(npzPath)

DB.sorted = False # we can remove this once we rebuild NPZ file
DB.groupData(siteID)

#graspObjs = DB.graspPackData(pathYAML, dirGRASP, orbHght, segIndGRASP)
grspObjs = DB.graspPackData(pathYAML)
gDB = graspDB(grspObjs)
rslts = gDB.processData(maxCPUs, binPathGRASP)
dill.dump_session(savePath)
