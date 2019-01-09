import matplotlib.pyplot as plt
import numpy as np
import sys
import os
sys.path.append(os.path.join("..", "GRASP_PythonUtils"))
from importlib import reload
import runGRASP
reload(runGRASP)
from runGRASP import graspDB


rsltsFile = '/Users/wrespino/Synced/Working/GRASP_testCases/Terra_20sites_bfdc446YAML_looseIRI_incldAero.pkl'
#dict_keys(['datetime', 'longitude', 'latitude', 'r', 'dVdlnr', 'rv', 'sigma', 'vol', 'sph', 'lambda', 'aod', 'aodMode', 'ssa', 'ssaMode', 'n', 'k', 'albedo', 'wtrSurf', 'brdf', 'sca_ang', 'meas_I', 'fit_I', 'aodAERO', 'aodDT', 'metaData', 'AEROloc'])
xVarNm = 'sigma'
#yVarNm = 'aod'
yVarNm = 'k'
xInd = 0
#yInd = [1]
#cVarNm= 'latitude'
cVarNm= False
cInd= 2
EEfunc = lambda x: 0.03+0.1*x

filterSites = True
frbdnSites = np.r_[179,1146,535,840,175]

gDB = graspDB()
gDB.loadResults(rsltsFile)

if filterSites:
    siteID = [rslt['AEROloc'][0] for rslt in gDB.rslts]
    vldInd = ~np.array([np.sum(a == frbdnSites, dtype=bool) for a in siteID])
else:
    vldInd = slice(None)


# DIFF PLOTS
#plt.figure()
#for i in range(0,4):
#    axHnd = plt.subplot(2,2,i+1)
#    gDB.diffPlot(xVarNm, yVarNm, i, i, lambdaFuncEE=EEfunc, FS=12, rsltInds=vldInd)

# SCAT PLOTS
plt.figure()
for i in range(0,4):
    axHnd = plt.subplot(2,2,i+1)
#    gDB.scatterPlot(xVarNm, yVarNm, i, i, cVarNm, cInd, one2oneScale=True, 
#                    logScl=True, customAx=axHnd, Rstats=True, rsltInds=vldInd)
    gDB.scatterPlot(xVarNm, yVarNm, [0], [0,i], cVarNm, cInd, one2oneScale=False, 
                    logScl=True, customAx=axHnd, Rstats=False, rsltInds=vldInd)

# HIST PLOTS
#plt.figure()
#bot = np.array([])
#top = np.array([])
#for i in range(0,4):
#    axHnd = plt.subplot(2,2,i+1)
#    index = np.r_[1,i] # THE LENGTH/1ST INDEX OF THIS MAY CHANGE FOR DIFFERENT VARS
#    gDB.histPlot(xVarNm, index, rsltInds=vldInd, customAx=axHnd)
#    vals = [rslt[xVarNm][tuple(index)] for rslt in gDB.rslts]
#    bot = np.r_[bot, np.percentile(vals,15)]
#    top = np.r_[top, np.percentile(vals,85)] 
#print('bottom 5%:' + ', '.join(['%7.5f' %y for y in bot]))
#print('top 95%:' + ', '.join(['%7.5f' %y for y in top]))
    
    
    
    
    