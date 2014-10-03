
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from neodb import db, core, neodb
import numpy as np


if __name__ == '__main__':
    # an Engine, which the Session will use for connection
    # resources
#     engine = create_engine('postgresql://postgres:postgres@192.168.2.2/demo')
#     
#     # create a configured "Session" class
#     Session = sessionmaker(bind=engine)
# 
#     # create a Session
#     session = Session()
    
#     engine, metadata, session = db.connect('postgres', 'postgres', '192.168.2.2', 'demo')
#     core.mapper(metadata)
#     query = session.query(core.AnalogSignalDB).filter(core.AnalogSignalDB.id == 1)
#     analogsignal = query.all()[0]
#     analogsignal.signal = np.fromstring(analogsignal.signal,dtype=np.int16)
#     
    ndb = neodb('postgres', 'postgres', '192.168.2.2', 'demo')
    kwargs = {"arg3": 3, "arg2": "two","arg1":5}
    k = {"arg3": 3}
    analogsignal = ndb.read("AnalogSignal", id=1)
    
    pass