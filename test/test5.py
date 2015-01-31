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

libcd.get_dc.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_float]
libcd.get_dc.restype = ctypes.c_float

libcd.get_cluster_dp.argtypes = [ctypes.c_char_p, ctypes.c_char_p, array_1d_double]
libcd.get_cluster_dp.restype = array_1d_double


dc = libcd.get_dc("dbname=demo host=192.168.2.2 user=postgres password=postgres", "54", 2.0)

local_density = np.empty(1026)
distance_to_higher_density = np.empty(1026)

print "dc: ", dc, type(dc)

libcd.get_local_density("dbname=demo host=192.168.2.2 user=postgres password=postgres", "54", dc, local_density, "gaussian")
libcd.get_distance_to_higher_density("dbname=demo host=192.168.2.2 user=postgres password=postgres", "54", local_density, distance_to_higher_density, len(local_density))

gamma = local_density*distance_to_higher_density

dp = distance_to_higher_density[local_density.argsort()]

dp2 = np.empty(1026)
for i in range(len(dp)):
    dp2[i] = i * dp[i]

#gamma2 = libcd.get_cluster_dp("dbname=demo host=192.168.2.2 user=postgres password=postgres", "54")

# plt.subplot(4,1,1)
# plt.plot(local_density, 'o')
# plt.subplot(4,1,2)
# plt.plot(distance_to_higher_density, 'o')
plt.subplot(2,1,1)
plt.plot(local_density, distance_to_higher_density, 'o')
plt.subplot(2,1,2)
plt.plot(dp2, 'o')
plt.show()

#print dc("dbname=demo user=postgres password=postgres hostaddr=192.168.2.2 port=5432")
pass