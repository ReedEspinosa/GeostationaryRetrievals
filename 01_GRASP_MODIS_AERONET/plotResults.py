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
#rsltsFile = '/Users/wrespino/Synced/Working/MODAERO_retrievalPickles/TerraLand_Select7Sites_YAMLc3db8ed_maxBlue_2lgnrm_Nt200_aeroWeakFrc_V2c.pkl' 
#rsltsFile = '/Users/wrespino/Synced/Working/MODAERO_retrievalPickles/TerraLand_Select7Sites_YAML6a60730_maxBlue_2lgnrm_Nt200_V2c.pkl'# LAND KING
#rsltsFile = '/Users/wrespino/Synced/Working/MODAERO_retrievalPickles/TerraLand_Select7Sites_YAMLd74dc51_maxBlue_2lgnrm_Nt200_V2c.pkl'
#rsltsFile = '/Users/wrespino/Synced/Working/MODAERO_retrievalPickles/Terra_GLOBALandDUSTsites_YAML1b95f26_3lgnrm_V2b.pkl' # the real ocean king!
#rsltsFile = '/Users/wrespino/Synced/Working/MODAERO_retrievalPickles/Terra_GLOBALandDUSTsites_YAML1b95f26_3lgnrm_Nt120_V2c.pkl' # A [OLD] new king! (and his challengers were defeated)
rsltsFile = '/Users/wrespino/Synced/Working/MODAERO_retrievalPickles/Terra_GLOBALandDUSTandLT20Msites_YAML1b95f26_3lgnrm_Nt120_V2c.pkl'
#dict_keys(['datetime', 'longitude', 'latitude', 'r', 'dVdlnr', 'rv', 'sigma', 'vol', 'sph', 'lambda', 'aod', 'aodMode', 'ssa', 'ssaMode', 'n', 'k', 'albedo', 'wtrSurf', 'brdf', 'sca_ang', 'meas_I', 'fit_I', 'aodAERO', 'aodDT', 'metaData', 'AEROloc'])
#xVarNm = 'aodAERO' #ANY CHANCE THIS IS BEING INTERPOLATED WRONG NOW?
#yVarNm = 'aodDT'
#yVarNm = 'aod'
#xVarNm = 'vol'
#xVarNm = 'angstrm0'
#yVarNm = 'angstrm1'
angLvarNm = ['lambda','lambda','lambdaDT',]
angTvarNm = ['aodAERO', 'aod', 'aodDT']
angInd = [[0,3], [0,3], [1,3]]
#cVarNm = 'brdf'
cVarNm= False
EEfunc = lambda x: 0.03+0.1*x

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

if xVarNm=='relRes' or yVarNm=='relRes' or cVarNm=='relRes':
    for i,rslt in enumerate(gDB.rslts):
        relRes = (rslt['aod']-rslt['aodAERO'])/(2*rslt['aod']+rslt['aodAERO'])
        gDB.rslts[i]['relRes'] = np.minimum(relRes,1e9) # can lower if needed

if xVarNm=='relResFit' or yVarNm=='relResFit' or cVarNm=='relResFit':
    for i,rslt in enumerate(gDB.rslts):
        relRes = (rslt['meas_I']-rslt['fit_I'])/(2*rslt['meas_I']+rslt['fit_I'])
        gDB.rslts[i]['relResFit'] = np.minimum(relRes,1e9) # can lower if needed

for j in range(3): # note we go 0-2
    if xVarNm=='angstrm%d'%j or yVarNm=='angstrm%d'%j or cVarNm=='angstrm%d'%j:
        for i,rslt in enumerate(gDB.rslts):
            lmbd = rslt[angLvarNm[j]][angInd[j]]
            tau = rslt[angTvarNm[j]][angInd[j]]
            gDB.rslts[i]['angstrm%d'%j] = angstrm(lmbd, tau) 

if filterSites:
    siteID = [rslt['AEROloc'][0] for rslt in gDB.rslts]
    vldInd = ~np.array([np.sum(a == frbdnSites, dtype=bool) for a in siteID])
#    vldInd = np.array([np.sum(a == keepSites, dtype=bool) for a in siteID])
else:
    vldInd = slice(None)

# MANUAL SCAT PLOT
#plt.figure(figsize=(6,5))
#gDB.scatterPlot(xVarNm, yVarNm, 0, 0, cVarNm, 3, one2oneScale=True, 
#                    logScl=False, customAx=plt.gca(), Rstats=True, rsltInds=vldInd, pltLabel=os.path.basename(rsltsFile))
#sys.exit()

# DIFF PLOTS
plt.figure(figsize=(9,5))
#axHnd = plt.subplot(2,1,2)
Nplts=6
for i in range(0,Nplts):
    axHnd = plt.subplot(2,int(np.ceil(Nplts/2)),i+1)
    gDB.diffPlot(xVarNm, yVarNm, i, i, lambdaFuncEE=EEfunc, FS=12, rsltInds=vldInd,
                 pltLabel=os.path.basename(rsltsFile), clnLayout=(i==(Nplts-1)))
##
## SCAT PLOTS
plt.figure(figsize=(9,5))
#axHnd = plt.subplot(2,1,1)
Nplts=6
for i in range(0,Nplts):
    axHnd = plt.subplot(2,int(np.ceil(Nplts/2)),i+1)
    gDB.scatterPlot(xVarNm, yVarNm, i, i, cVarNm, [0,i], one2oneScale=True, pltLabel=os.path.basename(rsltsFile),
                    logScl=True, customAx=axHnd, Rstats=True, rsltInds=vldInd, clnLayout=(i==(Nplts-1)))
#    ind1 = [0]
#    ind2 = [0,i]
#    ind3 = [0]
#    gDB.scatterPlot(xVarNm, yVarNm, ind1, ind2, cVarNm, ind3, one2oneScale=False, 
#                    logScl=False, customAx=axHnd, Rstats=False, rsltInds=vldInd, clnLayout=(i==(Nplts-1)))
#
# HIST PLOTS
#plt.figure(figsize=(9,6))
#bot = np.array([])
#median = np.array([]) 
#top = np.array([])
##clrVarVal = [rslt[cVarNm][3] for rslt in gDB.rslts]
##bestInd = np.abs(clrVarVal)<np.nanpercentile(np.abs(clrVarVal),20.0)
#Nplts= 9
#for i in range(0,Nplts):
#    axHnd = plt.subplot(2,int(np.ceil(Nplts/2)),i+1)
#    index = np.r_[0,i] # THE LENGTH/1ST INDEX OF THIS MAY CHANGE FOR DIFFERENT VARS
#    gDB.histPlot(xVarNm, index, rsltInds=vldInd, customAx=axHnd)
#    vals = np.array([rslt[xVarNm][tuple(index)] for rslt in gDB.rslts])
##    vals = vals[bestInd]
#    bot = np.r_[bot, np.percentile(vals,5)]
#    median = np.r_[median, np.percentile(vals,50)]
#    top = np.r_[top, np.percentile(vals,95)] 
#print('median:' + ', '.join(['%7.5f' %y for y in median]))
#print('bottom 5%:' + ', '.join(['%7.5f' %y for y in bot]))
#print('top 95%:' + ', '.join(['%7.5f' %y for y in top]))
#    
    
    
    
    