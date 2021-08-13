import numpy as np
from read_MODAERO import modaeroDB as mdb

# load modis terra collocation data
x = mdb('/Users/wrespino/Synced/GRASP_MODIS/Land_MOD04_L2_C61_V2c.out.npz')

# get indices at NASA_Ames
savePath = '/Users/wrespino/Synced/Proposals/ROSES_TASNPP_Yingxi_2020/retrievalSimulation/NASA_Ames_MOD_angles-SZA-VZA-PHI.txt'
ind = np.where(np.logical_and(np.isclose(x.aero_loc[:,-1],-122, rtol=0.01),np.isclose(x.aero_loc[:,-2],37.5, rtol=0.01)))[0]
Nind = len(ind)
ind = ind[round(0.11*Nind):round(0.8*Nind)] # these are very consecutive for Ames site

# pull first set of angles for each day (would be better to pick closest to site angles but little difference)
days = []
geomOut = []
for modloc,geom in zip(x.mod_loc[ind], x.geom[ind]):
    if len(days)==0 or not days[-1]==np.floor(modloc[-1]):
        days.append(np.floor(modloc[-1]))
        print('new day – %f' % (days[-1] - days[0]))
        sza = geom[0]
        thtv = geom[2]
        phi = geom[1] - geom[3]
        geomOut.append([sza, thtv, phi])
geomOut = np.asarray(geomOut)

# save the three angles (columns: [sza, thtv, phi]) on each day (rows) to a text file
fmt='%5.1f'
with open(savePath, "w") as fid:
    np.savetxt(fid, geomOut, fmt=fmt)
    