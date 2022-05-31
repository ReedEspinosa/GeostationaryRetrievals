#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pickle
import re

csvFileIn = '/Users/wrespino/Synced/RST_CAN-GRASP/AERONET_collocation/ABI16_ALM/ABI16_ALM_5lines_V1.csv'

with open(csvFileIn) as lines: 
    colNames = lines.readline().rstrip().split(','))

csvData = np.genfromtxt(csvFileIn, delimiter=',', skip_header=1)



pklFileOut = csvFileIn[:-3]+'pkl'
with open(pklFileOut, 'wb') as f:
    pickle.dump(self.rslts, f, pickle.HIGHEST_PROTOCOL)
print('Converted data saved to %s' % pklFileOut)