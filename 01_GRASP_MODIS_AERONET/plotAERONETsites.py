#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from read_MODAERO import modaeroDB
import os
os.environ["PROJ_LIB"] = "/Users/wrespino/anaconda3/share/proj" # fix for "KeyError: 'PROJ_LIB'" bug
from mpl_toolkits.basemap import Basemap

# NASA MacBook Air
basePath = '/Users/wrespino/Synced'
binPathGRASP = '/usr/local/bin/grasp'

#siteID 9=Wallops, 210=Nauru, 106=Ascension Island, 919=Fukue(Japan), 58=Cape Verde, Honolulu=97, Amr. Somoa=1185, Midway=278, Stennis=NOLA
#siteID = [210, 9, 106, 919, 58, 97, 1185, 278, 243, 179, 908, 809, 210, 1146, 103, 1101, 535, 739, 383, 840, 916, 475, 436, 667, 175] # can be int or list (not numpy array)
#siteID = [210, 9, 106, 919, 58, 97, 1185, 278, 243] # can be int or list (not numpy array)
#siteID = [987, 352, 119, 392, 769, 793, 1031, 187, 962, 469, 471, 470, 543, 461, 467, 906, 1108, 1109, 1114, 507, 256] # DUST SITES
#siteID = [9]
#siteID = np.unique(DB.aero_loc[:,0][DB.aero_loc[:,1]<40]) # find sites below 40m, only used after 1st run
siteID = np.r_[9, 161, 649, 1126, 1149, 518, 510, 268, 196, 514, 59, 1, 430, 221, 1150, 1197, 946, 33, 1077, 1107, 337, 806, 1102, 176, 339, 989, 658, 532, 1156, 746, 961, 981, 748, 316, 666, 485, 940, 922, 285, 126] # LAND SITES
# SEE "HINT" BELOW TO PLOT ALL SITES

npzPath = os.path.join(basePath, 'Remote_Sensing_Projects/GRASP_MODIS/Land_MOD04_L2_C61_V2c.out.npz')

DB = modaeroDB(npzPath, 'land')
#DB.groupData(siteID)

#siteID = np.unique(DB.aero_loc[:,0]) # HINT plot all sites
siteInd = [np.nonzero(DB.aero_loc[:,0]==ID)[0][0] for ID in siteID]
elev = DB.aero_loc[siteInd,1]
lat = DB.aero_loc[siteInd,2]
lon = DB.aero_loc[siteInd,3]

plt.figure(figsize=(12, 6))
m = Basemap(projection='robin', resolution='h', lat_0=0, lon_0=0)
m.bluemarble(scale=1);
#m.shadedrelief(scale=0.2)
x, y = m(lon, lat)
#plt.scatter(x, y, c='r', s=30, facecolors='none', edgecolors='r', cmap='plasma')
plt.plot(x, y, c='r', markersize=8, marker='o', mfc='none', linestyle='none')
[plt.text(x0,y0,'%d'%ID, color='r') for x0,y0,ID in zip(x,y,siteID)]
#plt.title('AERONET Sites')
#cbar = plt.colorbar()
##cbar.set_label("Elevation (m)", FontSize=14)
#cbar.set_label("80th percentile NCEP wind speed (m/s)", FontSize=14)