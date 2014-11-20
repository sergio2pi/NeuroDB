'''
Created on Oct 21, 2014

@author: sergio
'''

import numpy as np
import ctypes
import numpy.ctypeslib as npct
import matplotlib.pyplot as plt

#cfsfd = ctypes.cdll.LoadLibrary('/home/sergio/iibm/sandbox/t.so')
#cfsfd.get_dc.restype = ctypes.c_float
#dc = cfsfd.get_dc("dbname=demo host=192.168.2.2 user=postgres password=postgres", "54")
#print dc

array_1d_double = npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array_1d_int = npct.ndpointer(dtype=np.int64, ndim=1, flags='CONTIGUOUS')

libcd = npct.load_library("cfsfdp", "/home/sergio/iibm/workspace/NeuroDB/NeuroDB/cfunctions/cfsfdp")

libcd.get_local_density.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_float, array_1d_double, ctypes.c_char_p]
libcd.get_local_density.restype = ctypes.c_int

libcd.get_distance_to_higher_density.argtypes = [ctypes.c_char_p, ctypes.c_char_p, array_1d_double, array_1d_double, ctypes.c_int]
libcd.get_distance_to_higher_density.restype = ctypes.c_int

libcd.get_dc.argtypes = []
libcd.get_dc.restype = ctypes.c_float

dc = libcd.get_dc("dbname=demo host=192.168.2.2 user=postgres password=postgres", "54")

local_density = np.empty(1026)
distance_to_higher_density = np.empty(1026)

libcd.get_local_density("dbname=demo host=192.168.2.2 user=postgres password=postgres", "54", dc, local_density, "cutoff")
libcd.get_distance_to_higher_density("dbname=demo host=192.168.2.2 user=postgres password=postgres", "54", local_density, distance_to_higher_density, len(local_density))
print "dc: ", dc
print local_density
print distance_to_higher_density
plt.plot(local_density, distance_to_higher_density, 'o')
plt.show()

#print dc("dbname=demo user=postgres password=postgres hostaddr=192.168.2.2 port=5432")
pass