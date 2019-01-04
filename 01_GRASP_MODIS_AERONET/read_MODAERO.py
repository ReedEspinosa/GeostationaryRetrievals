# -*- coding: utf-8 -*-

import numpy as np
import glob
import warnings
import os
import sys
import hashlib
from math import ceil
sys.path.append(os.path.join("..", "GRASP_PythonUtils"))
from runGRASP import graspRun, pixel
from miscFunctions import angstrmIntrp

class modaeroDB(object):
    MOD_LAMDA = np.array([0.469, 0.555, 0.645, 0.8585, 1.24, 1.64, 2.13, 0.412, 0.443, 0.745])
    RFLCT_IND = np.r_[117:127]    
    AERO_LAMDA = np.array([1.64, 1.02, 0.87, 0.865, 0.779, 0.675, 0.667, 0.62, 0.56, 0.555, 0.551, 0.532, 0.531, 0.51, 0.5, 0.49, 0.443, 0.44, 0.412, 0.4, 0.38, 0.34, 0.554])
    AOD_IND = np.r_[2:25] 
    MOD_LOC_IND = np.r_[0:2, 50:52] # year, day, LAT, LON, DATENUM* (*added after calling sortData())
    AERO_LOC_IND = np.r_[46:50] # SITE_ID, ELEV, LAT, LON
    GEOM_IND = np.r_[52:57] # SOL_ZEN, SOL_AZM, SEN_ZEN, SEN_ASM, SCAT_ANG
    META_IND = np.r_[142:144] # GLINT_ANG, WIND_SPEED
    MOD_AOD_IND = np.r_[71:78] # DT retrieved AOD "Average", wavelengths 1st seven values of MOD_LAMDA
    hashObj = hashlib.sha1(MOD_LAMDA.tostring()+RFLCT_IND.tostring()+AERO_LAMDA.tostring()
        +AOD_IND.tostring()+MOD_LOC_IND.tostring()+AERO_LOC_IND.tostring()+GEOM_IND.tostring())
    SET_HASH = np.frombuffer(hashObj.digest(), dtype='uint32')[0] # we convert the first 32 bits to unsigned int
    GRP_DAY_LMT = np.r_[120] # maximum nt in GRASP build, each site seg will have this many unique datenums

    
    def __init__(self, loadPath=None):
        if loadPath is None or not self.loadData(loadPath):            
            self.rflct = np.array([]).reshape(0,self.RFLCT_IND.shape[0])
            self.aod = np.array([]).reshape(0,self.AOD_IND.shape[0])        
            self.mod_loc = np.array([]).reshape(0,self.MOD_LOC_IND.shape[0])
            self.aero_loc = np.array([]).reshape(0,self.AERO_LOC_IND.shape[0])
            self.geom = np.array([]).reshape(0,self.GEOM_IND.shape[0])
            self.metaData = np.array([]).reshape(0,self.META_IND.shape[0])
            self.modDT_aod = np.array([]).reshape(0,self.MOD_AOD_IND.shape[0])
            self.sorted = False    
        self.siteSegment = [] # grouped data is after saving stage
        
    def readFile(self, filePath):
        warnings.filterwarnings('ignore', message="some errors were detected") # HDF load error lines will produce warnings
        fileData = np.genfromtxt(filePath, delimiter='  ', invalid_raise=False) 
        warnings.resetwarnings()
        warnings.filterwarnings('ignore', message="invalid value encountered in less_equal") # HDF has "-np.nan" sometimes 
        fileData[fileData <= -9999] = np.nan
        warnings.resetwarnings()
        self.rflct = np.block([[self.rflct], [fileData[:, self.RFLCT_IND]]])
        self.aod = np.block([[self.aod], [fileData[:, self.AOD_IND]]])
        self.mod_loc = np.block([[self.mod_loc], [fileData[:, self.MOD_LOC_IND]]])
        self.aero_loc = np.block([[self.aero_loc], [fileData[:, self.AERO_LOC_IND]]])
        self.geom = np.block([[self.geom], [fileData[:, self.GEOM_IND]]])
        self.metaData = np.block([[self.metaData], [fileData[:, self.META_IND]]])
        self.modDT_aod = np.block([[self.modDT_aod], [fileData[:, self.MOD_AOD_IND]]])
        self.sorted = False     

    def readDIR(self, dirPath):
        filePaths = glob.glob(dirPath)
        print('%3d files found...' % np.shape(filePaths)[0])
        for filePath in filePaths:
            print('  %s' % filePath)
            self.readFile(filePath)
            
    def sortData(self):
        # Works from 1900-2100 & made sortData() 5x faster than w/ strptime(); strptime() version in Commit 9849d2a
        datenums = np.array([730120 + (yr-2000)*365 + ceil((yr-2000)/4)+dy for yr,dy in self.mod_loc[:,:2]])
        srtInd = datenums.argsort();
        self.rflct = self.rflct[srtInd,:]
        self.aod = self.aod[srtInd,:]
        self.mod_loc = np.block([[self.mod_loc[srtInd,:], datenums[srtInd,None]]])
        self.aero_loc = self.aero_loc[srtInd,:]
        self.geom = self.geom[srtInd,:]
        self.metaData = self.metaData[srtInd,:]
        self.modDT_aod = self.modDT_aod[srtInd,:]    
        self.sorted = True
    
    def groupData(self, siteIDFrc=False):
        if not self.sorted:
            self.sortData()
        if siteIDFrc:
            if type(siteIDFrc) is not list:
                siteList = [siteIDFrc]
            else:
                siteList = siteIDFrc
        else:
            siteList = np.unique(self.aero_loc[:,0])
        datenumSep = np.diff(self.mod_loc[:, -1]) # below won't catch same site on consecutive orbits,last of orbit != 0
        datenumSep[np.abs(np.diff(self.aero_loc[:,0]))>0] = -1 # in case same time registers at two sites
        for siteID in siteList:
            nowInd = np.nonzero(self.aero_loc[:,0] == siteID)[0]
            assert (len(nowInd)>0), "No sites found with ID %d" % siteID
            AEROloc = self.aero_loc[nowInd[0],:]
            if not np.array_equiv(AEROloc, self.aero_loc[nowInd,:]):
                warnings.warn('At least two AERONET site LAT/LON/ELEV were not the same within a single segment!')
            self.siteSegment.append(aeroSite(AEROloc, self.MOD_LAMDA, self.AERO_LAMDA)) # Only one site per segment
            for i in nowInd:
                self.siteSegment[-1].addMeas(self.rflct[i,:], self.aod[i,:], self.mod_loc[i,:],
                                             self.geom[i,:], self.modDT_aod[i,:], self.metaData[i,:])
                numDays = len(np.unique(np.atleast_2d(self.siteSegment[-1].mod_loc)[:,4]))
                if (numDays == self.GRP_DAY_LMT) and (datenumSep[i] != 0) and (i != nowInd[-1]):
                    self.siteSegment.append(aeroSite(AEROloc, self.MOD_LAMDA, self.AERO_LAMDA)) # This segment is full, start a new one
            # code to prevent siteSegments with only a couple days would go here
        [seg.condenceAOD() for seg in self.siteSegment] # Remove unused AOD wavelengths
        return self.siteSegment
                             
    def graspPackData(self, pathYAML, incldAERO=False, orbHghtKM=713, dirGRASP=False): # HINT: THIS WILL CHANGE FOR ix>1, land and AERONET included
        lndPrct = 0 # b/c ocean only for now
        graspObjs = []
        measLwrBnd = 0.00001 # minimum allowed value for radiance and AOD
        lambdaUsed = slice(0,7) # only use first 7 lambda for now, must be in ascending order
        for seg in self.siteSegment:
            gObj = graspRun(pathYAML, orbHghtKM, dirGRASP)
            gObj.AUX_dict = []
            unqDTs = np.unique(seg.mod_loc[:,-1])
#            unqDTs = np.unique(seg.mod_loc[:,-1])[0:3] # HACK to make run faster
            for unqDT in unqDTs:
                nowInd = np.nonzero(seg.mod_loc[:,-1] == unqDT)[0]
                MODlon = np.mean(seg.mod_loc[nowInd, 3])
                MODlat = np.mean(seg.mod_loc[nowInd, 2])
                nowPix = pixel(unqDT, 1, 1, MODlon, MODlat, seg.aero_loc[1], lndPrct)
                aodAERO = np.zeros(lambdaUsed.stop)
                for i,wl in enumerate(seg.MOD_LAMDA[lambdaUsed]): 
                    aodAERO[i] = angstrmIntrp(seg.AERO_LAMDA, seg.aod[nowInd[0],:], wl)
                    aodAEROinpt = aodAERO[i] if incldAERO else np.nan
                    msTyp = np.r_[41] if np.isnan(aodAEROinpt) else np.r_[41, 12] # normalized radiances, AOD
                    nip = msTyp.shape[0]  
                    sza = np.mean(seg.geom[nowInd, 0])
                    mu = np.cos(sza*np.pi/180) # [0] needed b/c sza might have dummyAng at end
                    dummyAng = [] if np.isnan(aodAEROinpt) else 0
                    thtv = np.r_[np.mean(seg.geom[nowInd, 2]), dummyAng]
                    phi = np.r_[np.mean(seg.geom[nowInd, 1] - seg.geom[nowInd, 3]), dummyAng]
                    radiance = max(np.mean(seg.rflct[nowInd,i])*mu, measLwrBnd) # MODIS R=L/FO*pi/mu0; GRASP R=L/FO*pi w/ R>1e-6                      
                    aodAEROinpt =  [] if np.isnan(aodAEROinpt) else max(aodAEROinpt, measLwrBnd)
                    nowPix.addMeas(wl, msTyp, np.repeat(1, nip), sza, thtv, phi, np.r_[radiance, aodAEROinpt])
                gObj.addPix(nowPix)
                aodDT = np.mean(seg.modDT_aod[nowInd,:], axis=0)
                metaData = np.mean(seg.metaData[nowInd,:], axis=0) # glint angle, wind speed NCEP
                aero_loc = seg.aero_loc # siteID, elev, lat, lon
                gObj.AUX_dict.append({'aodAERO':aodAERO, 'aodDT':aodDT, 'metaData':metaData, 'AEROloc':aero_loc})
            graspObjs.append(gObj)
        return graspObjs
            
    def saveData(self, filePath): # only saves data after readDir and sortData, grouped data is not inlcuded
        if not os.path.isdir(os.path.dirname(filePath)):
            warnings.warn('The directory containing the NPZ file path specified does not exist!')
            return False
        np.savez_compressed(filePath, rflct=self.rflct, aod=self.aod, mod_loc=self.mod_loc,
                            aero_loc=self.aero_loc, geom=self.geom, set_hash=self.SET_HASH,
                            sort=self.sorted, modDT_aod=self.modDT_aod, metaData=self.metaData)
        return True
        
    def loadData(self, filePath):
        if not os.path.isfile(filePath):
            warnings.warn('The file '+filePath+' does not exist!', stacklevel=2)
            return False
        loaded = np.load(filePath)
        if not np.equal(self.SET_HASH,loaded['set_hash']): # we need equal() because loaded set_hash is numpy array
            warnings.warn('The setting hash of the loaded file does not match the current defaults! Discarding loaded data...', stacklevel=2)
            return False
        self.rflct = loaded['rflct']
        self.aod = loaded['aod']
        self.mod_loc = loaded['mod_loc']
        self.aero_loc = loaded['aero_loc']
        self.geom = loaded['geom']
        self.sorted = loaded['sort']
        self.modDT_aod = loaded['modDT_aod']
        self.metaData = loaded['metaData']
        return True
        
    
class aeroSite(object):
    def __init__(self, location, MOD_LAMDA, AERO_LAMDA):
        self.aero_loc = location
        self.MOD_LAMDA = MOD_LAMDA
        self.AERO_LAMDA = AERO_LAMDA
        self.rflct = np.array([]).reshape(0,MOD_LAMDA.shape[0])
        self.aod = np.array([]).reshape(0,AERO_LAMDA.shape[0])
        self.mod_loc = np.array([]) # we don't know Ncol for the rest
        self.geom = np.array([])
        self.metaData = np.array([])
        self.modDT_aod = np.array([])        
        self.Nmeas = 0;
               
    def addMeas(self, rflct, aod, mod_loc, geom, modDT_aod=False, metaData=False):
        self.rflct = np.vstack([self.rflct, rflct])
        self.aod = np.vstack([self.aod, aod])
        self.mod_loc = np.vstack([self.mod_loc, mod_loc]) if self.mod_loc.size else mod_loc
        self.geom = np.vstack([self.geom, geom]) if self.geom.size else geom
        self.metaData = np.vstack([self.metaData, metaData]) if self.metaData.size else metaData
        self.modDT_aod = np.vstack([self.modDT_aod, modDT_aod]) if self.modDT_aod.size else modDT_aod
        self.Nmeas += 1
              
    def condenceAOD(self): # HINT: we may want to "condence" AERONET AOD to MODIS wavelengths here
        emptyInd = np.isnan(self.aod).all(0)
        if np.sum(~emptyInd) < 5:
            warnings.warn('Less than 5 valid wavelengths found at AERONET site %d' % self.aero_loc[0])        
        self.aod = np.delete(self.aod, emptyInd.nonzero(), 1)
        self.AERO_LAMDA = np.delete(self.AERO_LAMDA, emptyInd.nonzero())
        
        
        
        
        
        
        
        
        
