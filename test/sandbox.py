# funciones para el guardado de session online
class save_session_live():
    def __init__(self, id_project, id_individual, date, name, session_path,
                description = None, sample_rate = 14400, dtype = 'i2', 
                unit = 'mv', nbchannel = 25, scale_factor = 1):
        
        if type(date) == str:
            date = neodb.dbutils.get_datetimedate(date)
        elif type(date) != datetime.date:
            raise StandardError("Invalid date type. It must be 'datetime.date' or " + 
                                 "string with format 'dd-mm-yyyy' or 'dd/mm/yyyy'")
        
        if not os.path.isdir(session_path):
            raise StandardError("Invalid path.")
        
        self.in_queue =  Queue.Queue()
        self.out_queue =  Queue.Queue()
        
        self.save = save_session(id_project, id_individual, date, name,
                                 session_path, description, sample_rate, dtype, 
                                 unit, scale_factor, nbchannel, self.in_queue, self.out_queue)
        
    def start(self):
        self.save.start()
    
    def stop(self):
        self.in_queue.put('quit')
    
class save_segment(threading.Thread):
    def __init__(self, id_session, name, segmentfile_path, sample_rate, dtype, unit, nbchannel, index):
        threading.Thread.__init__(self)
        self.id_session = id_session
        self.name = name
        self.segmentfile_path = segmentfile_path
        self.sample_rate = sample_rate
        self.dtype = dtype
        self.index = index
        self.nbchannel = nbchannel
        self.unit = unit
    
    def run(self):
        segmentdb = neodb.core.SegmentDB(self.id_session, self.name, file_origin = self.segmentfile_path)
        id_segment = segmentdb.save(NDB)
        
        # Me aseguro que el archivo no se esta editando o que se haya terminado de copiar
        file = open(self.segmentfile_path)
        time.sleep(1)
        file.seek(0,2)
        while(file.read(1)):
            time.sleep(1)
        file.close()
        
        # Reading
        segment = io.RawBinarySignalIO(filename = self.segmentfile_path).read_segment(sampling_rate = self.sample_rate,          
                                                                             dtype = self.dtype,
                                                                             nbchannel = self.nbchannel, 
                                                                             rangemin = -16380, 
                                                                             rangemax = 16380)
        
        cindex = 0
        t_start = 0.0
        # Save analogsignals
        for channel in segment.analogsignals:
            analogsignaldb = neodb.core.AnalogSignalDB(id_segment = id_segment, 
                                                       name = self.name, 
                                                       signal = channel.tolist(), 
                                                       channel_index = cindex,
                                                       units = utils.get_quantitie(self.unit),
                                                       file_origin = self.segmentfile_path,
                                                       sampling_rate = self.sample_rate*quantities.Hz,
                                                       index = self.index,
                                                       t_start = t_start)
            id_analogsignal = analogsignaldb.save(NDB)
            cindex = cindex + 1
        
        tstart = len(channel)/self.sample_rate
        
        return id_analogsignal, tstart

class save_session(threading.Thread):  
    def __init__(self, id_project, id_individual, date, name, session_path,
                description, sample_rate, dtype, unit, scale_factor,
                nbchannel, in_queue, out_queue):
        global NDB
    
        if NDB == None:
            connect_db()
        threading.Thread.__init__(self)
        self.id_project = id_project
        self.id_individual = id_individual
        self.date = date
        self.name = name
        self.session_path = session_path
        self.description = description
        self.sample_rate = sample_rate
        self.dtype = dtype
        self.unit = unit
        self.scale_factor = scale_factor
        self.tstart = 0.0
        self.nbchannel = nbchannel
        
        self.in_queue = in_queue
        self.out_queue = out_queue
        
    def run(self):
        # create Session Block
        session = neodb.core.BlockDB(id_project = self.id_project,
                                     id_individual = self.id_individual, name = self.name, 
                                     rec_datetime = self.date, description = self.description)
        
        self.id_session = session.save(NDB)
        
        files = []
        msg = ''
        index = 0
        
        while True:
            name, segmentfile_path = utils.last_file(self.session_path)
            time.sleep(0.5)
            if (name != None) and (name not in files):
                files.append(name)
                sg = save_segment(self.id_session, self.name, segmentfile_path, self.sample_rate, self.dtype, self.unit, self.nbchannel, index)
                sg.start()
                index = index + 1
            try:
                msg = self.in_queue.get_nowait()
                if msg == 'quit':
                    break
            except:
                pass
        pass

def last_file(path):
    files = []
    
    while True:
        name, segmentfile_path = utils.last_file(path)
        if (name != None) and (name not in files):
            file = open(segmentfile_path)
            print "procesando: ", name, segmentfile_path
            time.sleep(1)
            file.seek(0,2)
            while(file.read(1)):
                time.sleep(1)
            
            print "terminado"
            file.close()
            files.append(name)
            print name