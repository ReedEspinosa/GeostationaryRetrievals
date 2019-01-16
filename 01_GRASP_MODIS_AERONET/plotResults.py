import matplotlib.pyplot as plt
import numpy as np
import sys
import os
sys.path.append(os.path.join("..", "GRASP_PythonUtils"))
from importlib import reload
import runGRASP
reload(runGRASP)
from runGRASP import graspDB

#rsltsFile = '/Users/wrespino/Synced/Working/MODAERO_retrievalPickles/Terra_GLOBALandDUSTsites_YAML7b69bf3_3lgnrm_V2b.pkl'
rsltsFile = '/Users/wrespino/Synced/Working/MODAERO_retrievalPickles/Terra_GLOBALandDUSTsites_YAMLce20bcb_3lgnrm_AEROincld_V2b.pkl'
#rsltsFile = '/Users/wrespino/Synced/Working/MODAERO_retrievalPickles/Terra_DUSTsites_02a0d00YAML_2lgnrm_V2b.pkl'
#rsltsFile = '/Users/wrespino/Synced/Working/MODAERO_retrievalPickles/Terra_20sites_b877e92YAML_newFineMode_V2b.pkl' # NEW REINING KING
#rsltsFile = '/Users/wrespino/Synced/Working/MODAERO_retrievalPickles/Terra_DUSTsites_02a0d00YAML_3lgnrm_AEROforced_V2b.pkl'
#dict_keys(['datetime', 'longitude', 'latitude', 'r', 'dVdlnr', 'rv', 'sigma', 'vol', 'sph', 'lambda', 'aod', 'aodMode', 'ssa', 'ssaMode', 'n', 'k', 'albedo', 'wtrSurf', 'brdf', 'sca_ang', 'meas_I', 'fit_I', 'aodAERO', 'aodDT', 'metaData', 'AEROloc'])
#xVarNm = 'aodAERO'
#yVarNm = 'aodDT'
#yVarNm = 'aod'
xVarNm = 'k'
yVarNm = 'rv'
#yInd = [1]
#cVarNm= 'relRes'
cVarNm= False
cInd= 2
EEfunc = lambda x: 0.03+0.1*x

filterSites = False
frbdnSites = np.r_[179,1146,535,840,175]

gDB = graspDB()
gDB.loadResults(rsltsFile)

for i,rslt in enumerate(gDB.rslts):
    relRes = (rslt['aod']-rslt['aodDT'])/(2*rslt['aod']+rslt['aodDT'])
    gDB.rslts[i]['relRes'] = relRes

if filterSites:
    siteID = [rslt['AEROloc'][0] for rslt in gDB.rslts]
    vldInd = ~np.array([np.sum(a == frbdnSites, dtype=bool) for a in siteID])
else:
    vldInd = slice(None)


# DIFF PLOTS
#plt.figure(figsize=(10,6))
##axHnd = plt.subplot(2,1,2)
#for i in range(0,6):
#    axHnd = plt.subplot(2,3,i+1)
#    gDB.diffPlot(xVarNm, yVarNm, i, i, lambdaFuncEE=EEfunc, FS=12, rsltInds=vldInd,
#                 pltLabel=os.path.basename(rsltsFile))

# SCAT PLOTS
plt.figure(figsize=(8,5))
#axHnd = plt.subplot(2,1,1)
for i in range(0,6):
    axHnd = plt.subplot(2,3,i+1)
#    gDB.scatterPlot(xVarNm, yVarNm, i, i, cVarNm, cInd, one2oneScale=True, pltLabel=os.path.basename(rsltsFile),
#                    logScl=True, customAx=axHnd, Rstats=True, rsltInds=vldInd)
    gDB.scatterPlot(xVarNm, yVarNm, [1,i], 1, cVarNm, 1, one2oneScale=False, 
                    logScl=False, customAx=axHnd, Rstats=False, rsltInds=vldInd)

# HIST PLOTS
#plt.figure(figsize=(9,6))
#bot = np.array([])
#median = np.array([]) 
#top = np.array([])
#for i in range(0,6):
#    axHnd = plt.subplot(2,3,i+1)
#    index = np.r_[i] # THE LENGTH/1ST INDEX OF THIS MAY CHANGE FOR DIFFERENT VARS
#    gDB.histPlot(xVarNm, index, rsltInds=vldInd, customAx=axHnd)
#    vals = [rslt[xVarNm][tuple(index)] for rslt in gDB.rslts]
#    bot = np.r_[bot, np.percentile(vals,25)]
#    median = np.r_[median, np.percentile(vals,50)]
#    top = np.r_[top, np.percentile(vals,75)] 
#print('median:' + ', '.join(['%7.5f' %y for y in median]))
#print('bottom 25%:' + ', '.join(['%7.5f' %y for y in bot]))
#print('top 75%:' + ', '.join(['%7.5f' %y for y in top]))
    
    
    
    
    