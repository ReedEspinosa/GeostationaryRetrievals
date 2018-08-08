# -*- coding: utf-8 -*-

# ultimatly we want a list of objects with all times at a given aeronet station

import numpy as np
import glob
import warnings
import os.path
import hashlib
import datetime as dt

class modaeroDB(object):
    MOD_LAMDA = np.array([0.469, 0.555, 0.645, 0.8585, 1.24, 1.64, 2.13, 0.412, 0.443, 0.745])
    RFLCT_IND = np.r_[117:127]    
    AERO_LAMDA = np.array([1.64, 1.02, 0.87, 0.865, 0.779, 0.675, 0.667, 0.62, 0.56, 0.555, 0.551, 0.532, 0.531, 0.51, 0.5, 0.49, 0.443, 0.44, 0.412, 0.4, 0.38, 0.34, 0.554])
    AOD_IND = np.r_[2:25] 
    MOD_LOC_IND = np.r_[0:2, 50:52] # year, day, LAT, LON, DATENUM* (*added after calling sortData())
    AERO_LOC_IND = np.r_[46:50] # SITE_ID, ELEV, LAT, LON
    GEOM_IND = np.r_[52:57] # SOL_ZEN, SOL_AZM, SEN_ZEN, SEN_ASM, SCAT_ANG
    hashObj = hashlib.sha1(MOD_LAMDA.tostring()+RFLCT_IND.tostring()+AERO_LAMDA.tostring()
        +AOD_IND.tostring()+MOD_LOC_IND.tostring()+AERO_LOC_IND.tostring()+GEOM_IND.tostring())
    SET_HASH = np.frombuffer(hashObj.digest(), dtype='uint32')[0] # we convert the first 32 bits to unsigned int
    GRP_DAY_LMT = np.r_[30] # gaps at a site of this many days will trigger new aeroSite

    
    def __init__(self, loadPath=None):
        if loadPath is None or not self.loadData(loadPath):            
            self.rflct = np.array([]).reshape(0,self.RFLCT_IND.shape[0])
            self.aod = np.array([]).reshape(0,self.AOD_IND.shape[0])        
            self.mod_loc = np.array([]).reshape(0,self.MOD_LOC_IND.shape[0])
            self.aero_loc = np.array([]).reshape(0,self.AERO_LOC_IND.shape[0])
            self.geom = np.array([]).reshape(0,self.GEOM_IND.shape[0])
            self.sorted = False
                
        
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
        self.sorted = False
        

    def readDIR(self, dirPath):
        filePaths = glob.glob(dirPath)
        print('%3d files found...' % np.shape(filePaths)[0])
        for filePath in filePaths:
            print('  %s' % filePath)
            self.readFile(filePath)

            
    def sortData(self):
        datenums = np.array([dt.strptime(str(yr), "%Y").toordinal()+dy for yr,dy in self.mod_loc[:,:2]])
        srtInd = datenums.argsort();
        self.rflct = self.rflct[srtInd,:]
        self.aod = self.aod[srtInd,:]
        self.mod_loc = np.block([[self.mod_loc[srtInd,:], datenums[srtInd,None]]])
        self.aero_loc = self.aero_loc[srtInd,:]
        self.geom = self.geom[srtInd,:]
        self.sorted = True

    
    def groupData(self):
        siteSegment = []
        if not self.sorted:
            self.sortData()
        siteList = np.unique(self.aero_loc[:,0]);
        for siteID in siteList:
            nowInd = self.aero_loc[:,0] == siteID
            datenumSep = np.diff(self.mod_loc[nowInd, 4])
            jumpInd = np.nonzero(datenumSep > self.GRP_DAY_LMT)[0]
            # We should now loop through the above, the ind of the N aeroSite will be jumpInd[N]:(jumpInd[N+1]-1)
            # Each new aeroSite should then be appended to siteSegment (see commented portion of temp2.py)
            
                
    def saveData(self, filePath):
        np.savez_compressed(filePath, rflct=self.rflct, aod=self.aod, mod_loc=self.mod_loc, aero_loc=self.aero_loc, geom=self.geom, set_hash=self.SET_HASH)        
        
        
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
        del(loaded)
        return True
        
        
        

class aeroSite(object):
    def __init__(self, location):
        self.aero_loc = location
        
        
        
        
        
        
        
        
        
        
        
        
        
        