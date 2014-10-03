'''
Created on Dec 22, 2013

@author: Sergio Hinojosa Rojas
'''

# import ctypes
# import NeuroDB.config as config
# import os
# import re
# import numpy as np
# import neo.core
# from neo import io
# import quantities
# import datetime
# import argparse
# 
# import threading
# import shutil

import threading
import NeuroDB.config as config
import argparse
import neodb.core
import datetime
import neodb.dbutils
import os
import quantities
from neo import io
import utils
import threading
import Queue
import time
import re
import signalProcessor.spike as Spike
import numpy as np

NDB = None

def connect_db(user=config.DBUSER, password=config.DBPASSWORD, hostname=config.DBHOSTNAME, dbname=config.DBNAME):
    """connect_db(user=config.DBUSER, password=config.DBPASSWORD, hostname=config.DBHOSTNAME, dbname=config.DBNAME)
       
       Create connection NDB that establishes a real DBAPI connection to the database.
       Neurodb uses a database schema create by NeoDB
       
       Return
       Function returns object SQLAlchemy Connection
       """
    global NDB
    try:
        NDB = neodb.dbconnect(dbname, user, password, hostname)
    except argparse.ArgumentError:
        pass
    
    return NDB


#### PROJECT ####


def create_project(name, date, description = None, projectcode = None):
    """create_project(name, description = None, date, projectcode)

       Creates project in database. Project may contains several individuals.

       Parameters
       name: Name of project. Mandatory.
       projectcode: Code of project. Number to identify the project. Optional
       description: Description or info about the project. Optional
       date: Date of creation. Mandatory.

       Return
       Function returns the id that identified it in database. Each project have 
       unique id.
       
       USE:
       
        create_project("Project1", "Testing project", "10/10/2014")
       """
    global NDB
    
    if NDB == None:
        connect_db()
    
    project = neodb.core.Project(name, date, description, projectcode)
    id = project.save(NDB)
    
    print "Project created. Id: %s"%id
    return id

def find_project(name = None, date_from = None, date_end = None):
    global NDB
    
    if NDB == None:
        connect_db()
    
    projects =[]
    
    if name != None:
        projects = neodb.get_id(NDB, 'project', name = name)
    else:
        if date_from != None and date_end != None:
            projects = neodb.get_id(NDB, 'project',
                                       date_start = date_from,
                                       date_end = date_end)
        else:
            if date_from != None:
                projects = neodb.get_id(NDB, 'project',
                                           date_start = date_from)
                
            elif date_end != None:
                projects = neodb.get_id(NDB, 'project',
                                           date_end = date_end)
    
    return projects

def get_project(id):
    global NDB
    
    if NDB == None:
        connect_db()
    
    project = neodb.core.Project().get_from_db(NDB, id)
    
    return project


#### INDIVIDUAL #####


def add_session(id_project, id_individual, date, name, session_path,
                description = None, sample_rate = 14400, dtype = 'i2', 
                unit = 'mv', scale_factor = 1, nchannels = 25):
    """
    add_session(id_project, id_individual, date, name, session_path,
                description = None, sample_rate = 14400, dtype = 'i2', 
                unit = 'mv', scale_factor = 1)
    
    Saves the records of a session.
    
    Parameters:
    id_project (mandatory)
    id_individual (mandatory)
    date (mandatory)
    name (mandatory)
    session_path (mandatory)
    description
    sample_rate
    dtype
    unit
    scale_factor
    """
    #TODO: implement "scale factor" as the factor that adjusts the measured value 
    # for an electric interpretation with unit.
    # Add nchannel parameter to create a recordingchannels objects that gatter the
    # segments of analogsignals
    
    global NDB
    
    if NDB == None:
        connect_db()
        
    if type(date) == str:
        date = neodb.dbutils.get_datetimedate(date)
    elif type(date) != datetime.date:
        raise StandardError("Invalid date type. It must be 'datetime.date' or " + 
                             "string with format 'dd-mm-yyyy' or 'dd/mm/yyyy'")
    
    if not os.path.isdir(session_path):
        raise StandardError("Invalid path.")
    
    # create Session Block
    session = neodb.core.BlockDB(id_project = id_project,
                                 id_individual = id_individual, name = name, 
                                 rec_datetime = date, description = description)
    
    id_session = session.save(NDB)
    
    # sort files of session
#     carpeta = session_path
#     archivos = []
#     for nombre in os.listdir(carpeta):
#         direccion = os.path.join(carpeta, nombre)
#         if os.path.isfile(direccion):  # solo archivos
#             archivos.append((os.path.getmtime(direccion), nombre, direccion))
    
    
    carpeta = session_path
    archivos = []
    for nombre in os.listdir(carpeta):
        match = re.match('(^(.*)-(\d+)-(\d+)$)', nombre)
        if match:
            direccion = os.path.join(carpeta, nombre)
            if os.path.isfile(direccion):  # solo archivos
                ndate = match.groups()[1]
                index = int(match.groups()[3])
                archivos.append((index, os.path.getmtime(direccion), nombre, direccion))
                
    archivos.sort()
    
    # Create recordingchannels
    rchannels = []  
    for n in range(nchannels):
        rc = neodb.core.RecordingChannelDB(id_block = id_session, index = n)
        id_rc = rc.save(NDB)
        rchannels.append((n, id_rc))
    
    # Create segments
    tstart = 0.0
    tlength = 0.0
    
    for nindex, _, name, segmentfile_path in archivos:
        segmentdb = neodb.core.SegmentDB(id_session, name, file_origin = segmentfile_path)
        id_segment = segmentdb.save(NDB)
        
        segment = io.RawBinarySignalIO(filename = segmentfile_path).read_segment(sampling_rate = sample_rate,          
                                                                             dtype = dtype,
                                                                             nbchannel=nchannels, 
                                                                             rangemin = -16380, 
                                                                             rangemax = 16380)
        
        cindex = 0
        
        # Save analogsignals
        for channel in segment.analogsignals:
            analogsignaldb = neodb.core.AnalogSignalDB(id_segment = id_segment,
                                                       id_recordingchannel =  rchannels[cindex][1],
                                                       name = name, 
                                                       signal = channel.tolist(), 
                                                       channel_index = rchannels[cindex][0],
                                                       units = utils.get_quantitie(unit),
                                                       file_origin = segmentfile_path,
                                                       sampling_rate = sample_rate*quantities.Hz,
                                                       t_start = tstart,
                                                       index = nindex)
            analogsignaldb.save(NDB)
            cindex = cindex + 1
        
        tlength = len(channel)/sample_rate
        tstart = tstart + tlength


def create_individual(name, description, birth_date, picture_path):
    """
    create_individual(name, description, birth_date, picture_path)
    
    Creates a individual in DB.
    
    USE:
    
    create_individual("Bolt", "Raton albino macho", "10/10/2013", "/home/user/pictures/bolt.jpg")
    
    """
    global NDB
    
    if NDB == None:
        connect_db()
    
    individual = neodb.core.Individual(name, description, birth_date, picture_path)
    id = individual.save(NDB)
    
    print "Individual created. Id: %s"%id
    return id

def find_individual(name = None, birth_date_from = None, birth_date_end = None):
    global NDB
    
    if NDB == None:
        connect_db()
    
    individuals =[]
    
    if name != None:
        individuals = neodb.get_id(NDB, 'individual_vw', name = name)
    else:
        if birth_date_from != None and birth_date_end != None:
            individuals = neodb.get_id(NDB, 'individual_vw',
                                       date_start = birth_date_from,
                                       date_end = birth_date_end)
        else:
            if birth_date_from != None:
                individuals = neodb.get_id(NDB, 'individual_vw',
                                           date_start = birth_date_from)
                
            elif birth_date_end != None:
                individuals = neodb.get_id(NDB, 'individual_vw',
                                           date_end = birth_date_end)
    
    return individuals

def get_individual(id):
    global NDB
    
    if NDB == None:
        connect_db()
    
    individual = neodb.core.Individual().get_from_db(NDB, id)
    
    return individual




def get_spikes(id_block, channel):
    """
    Process spikes of signals from session. Before run this functions you must run
    """
    global NDB
    
    if NDB == None:
        connect_db()
    
    ansig = neodb.core.analogsignaldb.get_from_db3(NDB, 
                                                   recordingchannel = channel, 
                                                   id_block = id_block)
    
    if len(ansig) != 0:
        ndbdetector = Spike.Detector()
        spikes, index, thr = ndbdetector.get_spikes(ansig)
        
        for i in range(len(spikes)):
            t_spike = index[i]/float(ansig.sampling_rate)
            segment = neodb.core.analogsignaldb.get_from_db3(NDB, 
                                                   recordingchannel = channel, 
                                                   id_block = id_block,
                                                   t_start = t_spike)
            #[(id_seg_ansig, _)] = neodb.get_id(NDB, "analogsignal", )
            spike = neodb.core.SpikeDB(id_segment = segment.id_segment, 
                                       id_recordingchannel = segment.id_recordingchannel,
                                       waveform = segment.signal,
                                       time = t_spike,
                                       index = index[i],
                                       sampling_rate = segment.sampling_rate)
            id_spike = spike.save(NDB)
        
        print ansig.signal
    else:
        return None
    
# def eliminar(ndb):
#     
#     query = "delete from public.analogsignal where name = '030712_1_00'"
#     cursor.execute(query)
#     
#     query = "delete from public.segment where name = '030712_1_00'"
#     cursor.execute(query)
#         
#     query = "delete from public.block where name = 'session3'"
#     cursor.execute(query)
#         
#     ndb.commit()    

if __name__ == '__main__':
    ndb = connect_db()
    cursor = ndb.cursor()
    #neodb.config_server('192.168.2.2', "oss", "ohm", "demo", "postgres", "postgres")
    #create_project("Project1", "Testing project", "10/10/2014")
    #image = "/home/sergio/iibm/workspace2/NeuroDB/test/feliz.jpg"
    #save = save_session_live(5, 1, '01/01/2014', 'asdfgqh', '/home/sergio/iibm/sandbox')
    #save.start()
    #create_individual("individuo2", "asdfg", "19-02-2013", image)
    
    #print get_individual(8)
    #eliminar(ndb)
    #add_session(14, 6, datetime.date.today(), 'session3', '/home/sergio/sandbox/test', description = "test")
    
    #eliminar(ndb)
    #last_file('/home/sergio/iibm/sandbox')

    get_spikes(54, 3)
    

    pass