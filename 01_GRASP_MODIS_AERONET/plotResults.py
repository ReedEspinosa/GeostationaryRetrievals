import matplotlib.pyplot as plt
import numpy as np
import sys
import os
sys.path.append(os.path.join("..", "GRASP_PythonUtils"))
from importlib import reload
import runGRASP
reload(runGRASP)
from runGRASP import graspDB

#rsltsFile = '/Users/wrespino/Synced/Working/MODAERO_retrievalPickles/TerraLand_Select8LandSites_YAML0e60019_maxBlue_3lgnrm_dlPar1_Nt300_V2c.pkl'
#rsltsFile = '/Users/wrespino/Synced/Working/MODAERO_retrievalPickles/TerraLand_Select8LandSites_YAML0e60019_maxBlue_3lgnrm_Nt300_AEROfrcd_V2c.pkl'
rsltsFile = '/Users/wrespino/Synced/Working/MODAERO_retrievalPickles/Terra_GLOBALandDUSTsites_YAML7b69bf3_3lgnrm_V2b.pkl' # A new king! (and his challengers were defeated)
#rsltsFile = '/Users/wrespino/Synced/Working/MODAERO_retrievalPickles_land/TerraLand_WallopsOnly_YAML4cae9bd_3lgnrm_AEROincld_V2b.pkl'
#dict_keys(['datetime', 'longitude', 'latitude', 'r', 'dVdlnr', 'rv', 'sigma', 'vol', 'sph', 'lambda', 'aod', 'aodMode', 'ssa', 'ssaMode', 'n', 'k', 'albedo', 'wtrSurf', 'brdf', 'sca_ang', 'meas_I', 'fit_I', 'aodAERO', 'aodDT', 'metaData', 'AEROloc'])
xVarNm = 'aodAERO' #ANY CHANCE THIS IS BEING INTERPOLATED WRONG NOW?
#yVarNm = 'aodDT'
yVarNm = 'aod'
#xVarNm = 'n'
#yVarNm = 'k'
#xVarNm = 'datetime'
#xVarNm = 'brdf'
#yInd = [1]
#cVarNm = 'vol'
cVarNm= False
EEfunc = lambda x: 0.03+0.1*x

filterSites = False
#frbdnSites = np.r_[179,1146,535,840,175] # need to uncomment/comment below to activate
keepSites = np.r_[9] #9,518,77,1,285,514,946,961


gDB = graspDB()
gDB.loadResults(rsltsFile)

if xVarNm=='relRes' or yVarNm=='relRes' or cVarNm=='relRes':
    for i,rslt in enumerate(gDB.rslts):
        relRes = (rslt['aod']-rslt['aodAERO'])/(2*rslt['aod']+rslt['aodAERO'])
        gDB.rslts[i]['relRes'] = np.minimum(relRes,1e9) # can lower if needed

## TESTING HACK!
#for i,rslt in enumerate(gDB.rslts):
#    exactEE = rslt['aodAERO']+EEfunc(rslt['aodAERO']-0.0001)
#    gDB.rslts[i]['aod'] = exactEE


if filterSites:
    siteID = [rslt['AEROloc'][0] for rslt in gDB.rslts]
#    vldInd = ~np.array([np.sum(a == frbdnSites, dtype=bool) for a in siteID])
    vldInd = np.array([np.sum(a == keepSites, dtype=bool) for a in siteID])
else:
    vldInd = slice(None)

# MANUAL SCAT PLOT
#plt.figure(figsize=(6,5))
#gDB.scatterPlot(xVarNm, yVarNm, 1, 1, cVarNm, 1, one2oneScale=False, 
#                    logScl=False, customAx=plt.gca(), Rstats=False, rsltInds=vldInd, pltLabel=os.path.basename(rsltsFile))

# DIFF PLOTS
plt.figure(figsize=(9,5))
#axHnd = plt.subplot(2,1,2)
Nplts=8
for i in range(0,Nplts):
    axHnd = plt.subplot(2,int(np.ceil(Nplts/2)),i+1)
    gDB.diffPlot(xVarNm, yVarNm, i, i, lambdaFuncEE=EEfunc, FS=12, rsltInds=vldInd,
                 pltLabel=os.path.basename(rsltsFile), clnLayout=(i==(Nplts-1)))

# SCAT PLOTS
plt.figure(figsize=(9,5))
#axHnd = plt.subplot(2,1,1)
Nplts=8
for i in range(0,Nplts):
    axHnd = plt.subplot(2,int(np.ceil(Nplts/2)),i+1)
    gDB.scatterPlot(xVarNm, yVarNm, i, i, cVarNm, 0, one2oneScale=True, pltLabel=os.path.basename(rsltsFile),
                    logScl=True, customAx=axHnd, Rstats=True, rsltInds=vldInd, clnLayout=(i==(Nplts-1)))
######    ind1= int(np.floor(i/2)) # 0, 0, 1
######    ind2= int(np.ceil((i)/2)+1) # 1 2 2
#    ind1 = [0,i]
#    ind2 = [0,i]
#    ind3 = [0]
#    gDB.scatterPlot(xVarNm, yVarNm, ind1, ind2, cVarNm, ind3, one2oneScale=False, 
#                    logScl=False, customAx=axHnd, Rstats=False, rsltInds=vldInd, clnLayout=(i==(Nplts-1)))

# HIST PLOTS
#plt.figure(figsize=(9,6))
#bot = np.array([])
#median = np.array([]) 
#top = np.array([])
#Nplts= 1
#for i in range(0,Nplts):
#    axHnd = plt.subplot(2,int(np.ceil(Nplts/2)),i+1)
#    index = np.r_[0,i] # THE LENGTH/1ST INDEX OF THIS MAY CHANGE FOR DIFFERENT VARS
#    gDB.histPlot(xVarNm, index, rsltInds=vldInd, customAx=axHnd)
#    vals = [rslt[xVarNm][tuple(index)] for rslt in gDB.rslts]
#    bot = np.r_[bot, np.percentile(vals,2)]
#    median = np.r_[median, np.percentile(vals,50)]
#    top = np.r_[top, np.percentile(vals,99.9)] 
#print('median:' + ', '.join(['%7.5f' %y for y in median]))
#print('bottom 2%:' + ', '.join(['%7.5f' %y for y in bot]))
#print('top 99.9%:' + ', '.join(['%7.5f' %y for y in top]))
#    
    
    
    
    