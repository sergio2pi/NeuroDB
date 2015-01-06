'''
Created on Dec 22, 2014

@author: sergio
'''

import numpy as np
import ctypes
import numpy.ctypeslib as npct
import matplotlib.pyplot as plt

array_1d_double = npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array_1d_int = npct.ndpointer(dtype=np.int64, ndim=1, flags='CONTIGUOUS')
array_2d_double = npct.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')

libcd = npct.load_library("cfsfdp", "/home/sergio/iibm/workspace/NeuroDB/NeuroDB/cfunctions/cfsfdp")

libcd.get_centers_cluster_dp.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
libcd.get_centers_cluster_dp.restype = array_1d_double

a = libcd.get_centers_cluster_dp("dbname=demo host=192.168.2.2 user=postgres password=postgres", "54");

pass