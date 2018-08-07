# -*- coding: utf-8 -*-

# ultimatly we want a list of objects with all times at a given aeronet station

from numpy import array, r_, block, genfromtxt, nan, shape, savez_compressed, load
import glob
import warnings
import os.path


class modaeroDB(object):
    MOD_LAMDA = array([0.469, 0.555, 0.645, 0.8585, 1.24, 1.64, 2.13, 0.412, 0.443, 0.745])
    RFLCT_IND = r_[117:127]    
    AERO_LAMDA = array([1.64, 1.02, 0.87, 0.865, 0.779, 0.675, 0.667, 0.62, 0.56, 0.555, 0.551, 0.532, 0.531, 0.51, 0.5, 0.49, 0.443, 0.44, 0.412, 0.4, 0.38, 0.34, 0.554])
    AOD_IND = r_[2:25] 
    MOD_LOC_IND = r_[0:2, 50:52] # year, day, LAT, LON
    AERO_LOC_IND = r_[46:50] # SITE_ID, ELEV, LAT, LON
    GEOM_IND = r_[52:57] # SOL_ZEN, SOL_AZM, SEN_ZEN, SEN_ASM, SCAT_ANG
    SET_HASH = hash(MOD_LAMDA.tostring()+RFLCT_IND.tostring()+AERO_LAMDA.tostring()
        +AOD_IND.tostring()+MOD_LOC_IND.tostring()+AERO_LOC_IND.tostring()+GEOM_IND.tostring())
    
    def __init__(self, loadPath=None):
        if loadPath is None or not self.loadData(loadPath):            
            self.rflct = array([]).reshape(0,self.RFLCT_IND.shape[0])
            self.aod = array([]).reshape(0,self.AOD_IND.shape[0])        
            self.mod_loc = array([]).reshape(0,self.MOD_LOC_IND.shape[0])
            self.aero_loc = array([]).reshape(0,self.AERO_LOC_IND.shape[0])
            self.geom = array([]).reshape(0,self.GEOM_IND.shape[0])                
                
        
    def readFile(self, filePath):
        warnings.filterwarnings('ignore', message="some errors were detected") # HDF load error lines will produce warnings
        fileData = genfromtxt(filePath, delimiter='  ', invalid_raise=False) 
        warnings.resetwarnings()
        
        warnings.filterwarnings('ignore', message="invalid value encountered in less_equal") # HDF has "-nan" sometimes 
        fileData[fileData <= -9999] = nan
        warnings.resetwarnings()
        
        self.rflct = block([[self.rflct], [fileData[:, self.RFLCT_IND]]])
        self.aod = block([[self.aod], [fileData[:, self.AOD_IND]]])
        self.mod_loc = block([[self.mod_loc], [fileData[:, self.MOD_LOC_IND]]])
        self.aero_loc = block([[self.aero_loc], [fileData[:, self.AERO_LOC_IND]]])
        self.geom = block([[self.geom], [fileData[:, self.GEOM_IND]]])


    def readDIR(self, dirPath):
        filePaths = glob.glob(dirPath)
        print('%3d files found...' % shape(filePaths)[0])
        for filePath in filePaths:
            print('  %s' % filePath)
            self.readFile(filePath)
            
                
    def saveData(self, filePath):
        savez_compressed(filePath, rflct=self.rflct, aod=self.aod, mod_loc=self.mod_loc, aero_loc=self.aero_loc, geom=self.geom, set_hash=self.SET_HASH)        
        
        
    def loadData(self, filePath):
        if not os.path.isfile(filePath):
            warnings.warn('The file '+filePath+' does not exist!')
            return False

        loaded = load(filePath)
        if self.SET_HASH != loaded['set_hash']:
            warnings.warn('The setting hash of the loaded file does not match the current defaults! Discarding loaded data...')
            return False

        self.rflct = loaded['rflct']
        self.aod = loaded['aod']
        self.mod_loc = loaded['mod_loc']
        self.aero_loc = loaded['aero_loc']
        self.geom = loaded['geom']
        del(loaded)
        return True
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        