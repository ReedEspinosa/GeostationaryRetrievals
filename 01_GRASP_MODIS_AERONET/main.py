#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from read_MODAERO import modaeroDB
#import matplotlib.pyplot as plt

dirPath = '/Users/wrespino/Synced/Remote_Sensing_Projects/GRASP_MODIS/ocean/MOD*.out'
npzPath = '/Users/wrespino/Synced/Remote_Sensing_Projects/GRASP_MODIS/Ocean_MOD04_L2_C61.out.npz'

DB = modaeroDB(npzPath)

if DB.aod.shape[0]==0:
    DB.readDIR(dirPath)
    DB.saveData(npzPath)

#scatAng = DB.geom[:,4]
#
## plot results
#line = plt.figure()
#plt.hist(scatAng, bins=200)


