#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

#sdataFN = '/Volumes/MagHDD/LACO/Remote_Sensing_Projects/GRASP_MODIS/GRASP_modisOnly/modeoro_Nauru%d.txt'


Nfiles = 6
maxAOD = 0.7
Npsd = 16
ttleNm = 'Wallops'

aod550AERO = []
aod550GRASP = []
rri550GRASP = []
iri550GRASP = []
sph550GRASP = []
dVdlnr = np.array([]).reshape(0,Npsd)
for i in range(Nfiles):
    # MODIS GRASP AODs
#    grspRslt = graspObjs[0].readOutput(sdataFN % i)
    sdataFN = os.path.join(graspObjs[i].dirGRASP, graspObjs[i].findStream_FN())
    grspRslt = graspObjs[0].readOutput(sdataFN)
    aod550GRASP += [rslt['aod'][1] for rslt in grspRslt]
    rri550GRASP += [rslt['n'][1] for rslt in grspRslt]
    iri550GRASP += [rslt['k'][1] for rslt in grspRslt]
    sph550GRASP += [rslt['sph'] for rslt in grspRslt]
    for rslt in grspRslt:
        dVdlnr = np.vstack([dVdlnr, rslt['dVdlnr']])
    
    # AERONET AODs
    seg = DB.siteSegment[i]
    unqDTs = np.unique(seg.mod_loc[:,-1])
    for unqDT in unqDTs:
        nowInd = np.nonzero(seg.mod_loc[:,-1] == unqDT)[0]
        aod550AERO.append(seg.aod[nowInd[0],-1])

aod550AERO = np.array(aod550AERO)
aod550GRASP = np.array(aod550GRASP)
aod550GRASPcln = aod550GRASP[~np.isnan(aod550AERO)]
aod550AERO = aod550AERO[~np.isnan(aod550AERO)]

# compare with AERONET
line = plt.figure()
plt.plot(aod550AERO, aod550GRASPcln, '.')
plt.plot(np.r_[0, maxAOD], np.r_[0, maxAOD], 'k')
plt.xlabel("AOD AERONET (550nm)")
plt.ylabel("AOD MODIS/GRASP (550nm)")
plt.xlim([0, maxAOD])
plt.ylim([0, maxAOD])
Rcoef = np.corrcoef(aod550AERO, aod550GRASPcln)[0,1]
RMSE = np.sqrt(np.mean((aod550AERO - aod550GRASPcln)**2))
textstr = 'N=%d\nR=%.3f\nRMS=%.3f\n'%(len(aod550AERO), Rcoef, RMSE)
plt.text(0.7*maxAOD, 0.03, textstr, fontsize=12)
plt.title(ttleNm)


# plot PSD
radii = grspRslt[0]['r']
line = plt.figure()
N = len(dVdlnr)
[plt.semilogx(radii, sngdV, color=[i/N, 0.2, 1/(i+1)], linewidth=0.3) for i,sngdV in enumerate(dVdlnr)]
plt.semilogx(radii, np.mean(dVdlnr, axis=0), 'k')
plt.title(ttleNm + ' (blue=old, red=recent)')
plt.xlabel("radius (um)")
plt.ylabel("dV/dlnr (normalized)")

line = plt.figure()
plt.hist(rri550GRASP, 30)
plt.title(ttleNm)
plt.xlabel("Real Refractive Index (550nm)")
plt.ylabel("Frequency")

line = plt.figure()
plt.hist(iri550GRASP, 60)
plt.title(ttleNm)
plt.xlabel("Imaginary Refractive Index (550nm)")
plt.ylabel("Frequency")

line = plt.figure()
plt.hist(sph550GRASP, 30)
plt.title(ttleNm)
plt.xlabel("Spherical Fraction")
plt.ylabel("Frequency")


####################
#import sys
#import os
#sys.path.append(os.path.join("..", "GRASP_PythonUtils"))
#from runGRASP import graspRun

#def hello(a,b=8):
#    print(a)
#    print(b)
#    
#hello(99,90)

#import sys
#sys.path.insert(0, '../GRASP_PythonUtils/')
#from read_MODAERO import modaeroDB
#from .tester import bot
#
#x = bot()
#greeting('I dont think it will work')

#class top(object):
#    
#    def __init__(self):
#        self.botObj = bot()
#        
#    def giveBack(self,value):
#        self.botObj.setVal(value)
#        print(self.botObj.aod)
#        return self.botObj
#    
#class bot(object):
#    
#    def __init__(self):
#        self.aod = 1
#        
#    def setVal(self, value):
#        self.aod = value