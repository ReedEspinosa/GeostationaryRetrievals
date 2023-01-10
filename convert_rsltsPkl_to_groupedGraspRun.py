#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from os import path
import numpy as np
sys.path.append("/Users/wrespino/Synced/Local_Code_MacBook/GSFC-GRASP-Python-Interface")
from runGRASP import graspDB, graspRun, pixel

# paths
pklInputPath = '/Users/wrespino/Synced/RST_CAN-GRASP/AERONET_collocation/ABI16_ALM_TestFiles/ABI16_ALM_5lines_V1.pkl'
# pklInputPath = '/Users/wrespino/Synced/RST_CAN-GRASP/AERONET_collocation/ABI16/ABI16_ALM.pkl'
# pklInputPath = '/Users/wrespino/Synced/RST_CAN-GRASP/AERONET_collocation/ABI17/ABI17_ALM.pkl'

# load the data
gDB_in = graspDB()
rslts = gDB_in.loadResults(pklInputPath)
del gDB_in

# (1) determine grouping of rslts

# (2) convert each rslt to a pixel objects

# (3) combine (1) and (2) to pack them all into graspRun members of gDB

# (4) save results to a pkl... gDB does not have a method for this (if above is fast we could just run here...)
