'''
Created on Dec 22, 2014

@author: sergio
'''

import numpy as np
import ctypes
import numpy.ctypeslib as npct
import matplotlib.pyplot as plt
import psycopg2

import neodb.core
username = 'postgres'
password = 'postgres'
host = '192.168.2.2'
dbname = 'demo'
url = 'postgresql://%s:%s@%s/%s'%(username, password, host, dbname)
dbconn = psycopg2.connect('dbname=%s user=%s password=%s host=%s'%(dbname, username, password, host))

array_1d_double = npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array_1d_int = npct.ndpointer(dtype=np.int64, ndim=1, flags='CONTIGUOUS')
array_2d_double = npct.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')

libcd = npct.load_library("cfsfdp", "/home/sergio/iibm/workspace/NeuroDB/NeuroDB/cfunctions/cfsfdp")

libcd.get_centers_cluster_dp.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, array_1d_double, ctypes.c_double]
libcd.get_centers_cluster_dp.restype = ctypes.c_int

libcd.get_cluster_dp.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_double, array_1d_double, ctypes.c_double]
libcd.get_cluster_dp.restype = ctypes.c_int

libcd.get_dc.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_float]
libcd.get_dc.restype = ctypes.c_float

connect = "dbname=demo host=192.168.2.2 user=postgres password=postgres"
id_block = "54"
channel = "3"

dc = libcd.get_dc(connect, id_block, channel, 2.0)
cluster1 = np.empty(1026)
cluster2 = np.empty(1026)
centers = np.empty(1026)

libcd.get_centers_cluster_dp(connect, id_block, channel, centers, dc)
libcd.get_cluster_dp(connect, id_block, channel,  centers[1], cluster1, dc)
libcd.get_cluster_dp(connect, id_block, channel,  centers[2], cluster2, dc)

plt.subplot(2,1,1)
for i in range(int(cluster1[0])):
    if i != 0:
        spikes = neodb.core.spikedb.get_from_db(dbconn, id_block = 54, channel = 3, id = int(cluster1[i]))
        signal = spikes[0].waveform
        plt.plot(signal)

plt.subplot(2,1,2)
for i in range(int(cluster2[0])):
    if i != 0:
        spikes = neodb.core.spikedb.get_from_db(dbconn, id_block = 54, channel = 3, id = int(cluster2[i]))
        signal = spikes[0].waveform
        plt.plot(signal)
plt.show()
    
pass