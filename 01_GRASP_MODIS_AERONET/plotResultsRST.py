import matplotlib.pyplot as plt
import numpy as np
import sys
import os
sys.path.append(os.path.join("..", "..", "GRASP_Scripts"))
from importlib import reload
import runGRASP
reload(runGRASP)
from runGRASP import graspDB
from miscFunctions import angstrm
os.environ["PROJ_LIB"] = "/Users/wrespino/anaconda3/share/proj" # fix for "KeyError: 'PROJ_LIB'" bug
from mpl_toolkits.basemap import Basemap


#rsltsFile = '/Users/wrespino/Synced/Working/MODAERO_retrievalPickles/Terra_GLOBALandDUSTsites_YAML1b95f26_3lgnrm_Nt120_V2c.pkl' # A [OLD] new king! (and his challengers were defeated)
rsltsFile = '/Users/wrespino/Synced/Working/MODAERO_retrievalPickles/Terra_GLOBALandDUSTandLT20Msites_YAML1b95f26_3lgnrm_Nt120_V2c.pkl'
#dict_keys(['datetime', 'longitude', 'latitude', 'r', 'dVdlnr', 'rv', 'sigma', 'vol', 'sph', 'lambda', 'aod', 'aodMode', 'ssa', 'ssaMode', 'n', 'k', 'albedo', 'wtrSurf', 'brdf', 'sca_ang', 'meas_I', 'fit_I', 'aodAERO', 'aodDT', 'metaData', 'AEROloc'])
angLvarNm = ['lambda','lambda','lambdaDT',]
angTvarNm = ['aodAERO', 'aod', 'aodDT']
angInd = [[1,3], [1,3], [1,3]]
#cVarNm = 'brdf'
cVarNm= False
EEfunc = lambda x: 0.03+0.1*x

siteTune = [210, 9, 106, 919, 58, 97, 1185, 278, 243, 987, 352, 119, 392, 769, 793, 1031, 187, 962, 469, 471, 470, 543, 461, 467, 906, 1108, 1109, 1114, 507, 256]

filterSites = True
frbdnSites = np.r_[946,179,1146,535,840, 175, 436, 778, 1185, 480,  771, 511, 804,  982, 174] # need to uncomment/comment below to activate
#keepSites = np.r_[961] #9,518,77,1,285,514,946,961
#9
#518, N=2
#77, 1,285,514,
#946, too low aero aod, GRASP and DT
#961


gDB = graspDB()
gDB.loadResults(rsltsFile)

for j in range(2): # note we go 0-2
    for i,rslt in enumerate(gDB.rslts):
        lmbd = rslt[angLvarNm[j]][angInd[j]]
        tau = rslt[angTvarNm[j]][angInd[j]]
        gDB.rslts[i]['angstrm%d'%j] = angstrm(lmbd, tau) 

if filterSites:
    siteID = np.array([rslt['AEROloc'][0] for rslt in gDB.rslts])
    vldInd = ~np.array([np.sum(a == frbdnSites, dtype=bool) for a in siteID])
else:
    vldInd = slice(None)

scatAODLim = 2.6
scatAELim = 2.6
diffAODLim = 0.49
aodAEROlbl = '$\\tau_{AERONET}$'
aodGRASPlbl = '$\\tau_{DTGRASP}$'
aodDTlbl = '$\\tau_{DarkTarget}$'
lblFS = 14

# Generate Plots
plt.figure(figsize=(9,5.3))

axHnd = plt.subplot(2,3,1)
gDB.scatterPlot('aodAERO', 'aod', 1, 1, False, 0, one2oneScale=True, clnLayout=False,
                    logScl=False, customAx=axHnd, Rstats=True, rsltInds=vldInd, FS=12)
axHnd.set_ylim(0,scatAODLim), axHnd.set_xlim(0,scatAODLim)
axHnd.set_xlabel(aodAEROlbl, FontSize=lblFS)
axHnd.set_ylabel(aodGRASPlbl, FontSize=lblFS)
axHnd.set_xticks(np.arange(0,scatAODLim,0.5))
axHnd.set_yticks(np.arange(0,scatAODLim,0.5))

axHnd = plt.subplot(2,3,4)
gDB.scatterPlot('aodAERO', 'aodDT', 1, 1, False, 0, one2oneScale=True, clnLayout=False, 
                    logScl=False, customAx=axHnd, Rstats=True, rsltInds=vldInd, FS=12)
axHnd.set_ylim(0,2.5), axHnd.set_xlim(0,2.5)
axHnd.set_ylim(0,scatAODLim), axHnd.set_xlim(0,scatAODLim)
axHnd.set_xlabel(aodAEROlbl, FontSize=lblFS)
axHnd.set_ylabel(aodDTlbl, FontSize=lblFS)
axHnd.set_xticks(np.arange(0,scatAODLim,0.5))
axHnd.set_yticks(np.arange(0,scatAODLim,0.5))


axHnd = plt.subplot(2,3,2)
gDB.diffPlot('aodAERO', 'aod', 1, 1, lambdaFuncEE=EEfunc, FS=12, clnLayout=False, rsltInds=vldInd)
axHnd.set_ylim(-diffAODLim,diffAODLim)
axHnd.set_xlabel(aodAEROlbl, FontSize=lblFS)
axHnd.set_ylabel(aodGRASPlbl+' $-$ '+aodAEROlbl, FontSize=lblFS)

axHnd = plt.subplot(2,3,5)
gDB.diffPlot('aodAERO', 'aodDT', 1, 1, lambdaFuncEE=EEfunc, FS=12, clnLayout=False, rsltInds=vldInd)
axHnd.set_ylim(-diffAODLim,diffAODLim)
axHnd.set_xlabel(aodAEROlbl, FontSize=lblFS)
axHnd.set_ylabel(aodDTlbl+' $-$ '+aodAEROlbl, FontSize=lblFS)

axHnd = plt.subplot(2,3,3)
gDB.scatterPlot('angstrm0', 'angstrm1', 0, 0, False, 0, one2oneScale=True, clnLayout=False, 
                    logScl=False, customAx=axHnd, Rstats=True, rsltInds=vldInd, FS=12)
axHnd.set_ylim(0,scatAELim), axHnd.set_xlim(0,scatAELim)
axHnd.set_xlabel('$AE_{AERONET}$', FontSize=lblFS)
axHnd.set_ylabel('$AE_{DTGRASP}$', FontSize=lblFS)
axHnd.set_xticks(np.arange(0,scatAELim,0.5))
axHnd.set_yticks(np.arange(0,scatAELim,0.5))

axHnd = plt.subplot(2,3,6)
lon = []
lat = []
tuneInd = np.array([],dtype='bool')
for st in np.unique(siteID):
    stInd = np.nonzero(siteID==st)[0]
    if stInd.shape[0] > 0:
        lon = np.r_[lon, gDB.rslts[stInd[0]]['longitude']]
        lat = np.r_[lat, gDB.rslts[stInd[0]]['latitude']]
        tuneInd = np.r_[tuneInd, (st in siteTune)]
m = Basemap(projection='merc', llcrnrlon=-170,llcrnrlat=-62,urcrnrlon=191,urcrnrlat=84)
#m = Basemap(projection='merc', llcrnrlon=-170,llcrnrlat=-80,urcrnrlon=191,urcrnrlat=85)
#m.drawcoastlines(color=[0.5,0.5,0.5])
m.drawmapboundary(color='black')
m.fillcontinents(color=[0.7,0.8,0.7],lake_color='white')
x, y = m(lon, lat)
plt.plot(x, y, c='b', markersize=2, marker='.', linestyle='none')
plt.plot(x[tuneInd], y[tuneInd], c='r', markersize=2, marker='.', linestyle='none')

plt.tight_layout()

pos1 = axHnd.get_position() # get the original position 
pos2 = [pos1.x0 - 0.062, pos1.y0 - 0.053,  pos1.width*1.3, pos1.height*1.3] 
axHnd.set_position(pos2) # set a new position
axHnd.set_xlabel('Tuning (red) & test (blue) sites', FontSize=12, labelpad=17)