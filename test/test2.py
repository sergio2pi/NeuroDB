import sys
import re
#sys.path.append("/home/sergio/iibm/workspace2/NeoDB")
#sys.path.append("/home/sergio/iibm/workspace2/NeuroDB/src")

import NeuroDB.neurodb as ndb
#id_individual = ndb.create_individual("Elektra","Raton albino hembra","10/02/2014","/home/sergio/Pictures/elektra.jpg")
#id_project = ndb.create_project("Elektra", "05/05/2014", "Testing project.")
id_session = ndb.add_session(17, 15, "05/05/2014", "Session 1", "/home/sergio/iibm/registro/elektra", "Description of session 1.")

# NDB = ndb.connect_db()
# cursor = NDB.cursor()
# 
# query = """select analogsignal.index
#            from analogsignal join segment
#            on id_segment = segment.id
#            and segment.id_block = 50
#            and channel_index = 0"""
#            
# cursor.execute(query)
# results = cursor.fetchall()
# 
# query = """select analogsignal.signal
#            from analogsignal join segment
#            on id_segment = segment.id
#            and segment.id_block = 50
#            and analogsignal.name = '030712_1_2'"""
# 
# cursor.execute(query)
# results = cursor.fetchall()
# 
# query = """select analogsignal.index,
#            analogsignal.signal,
#            from analogsignal join segment
#            on id_segment = segment.id
#            and segment.id_block = 50
#            and channel_index = 0"""
#            
# 
# arch = ['030712-1-2', '030712-1-36', '030712-1-30', '030712-1-28', '030712-1-45', '030712-1-32', '030712-1-18', '030712-1-13', '030712-1-29', '030712-1-10', '030712-1-46', '030712-1-26', '030712-1-20', '030712-1-6', '030712-1-33', '030712-1-43', '030712-1-25', '030712-1-24', '030712-1-34', '030712-1-16', '030712-1-48', '030712-1-35', '030712-1-44', '030712-1-47', '030712-1-7', '030712-1-11', '030712-1-37', '030712-1-0', '030712-1-14', '030712-1-38', '030712-1-49', '030712-1-27', '030712-1-17', '030712-1-15', '030712-1-19', '030712-1-3', '030712-1-21', '030712-1-4', '030712-1-39', '030712-1-8', '030712-1-12', '030712-1-40', '030712-1-41', '030712-1-31', '030712-1-42', '030712-1-50', '030712-1-1', '030712-1-23', '030712-1-22', '030712-1-5']
# 
# for name in arch:
#     match = re.match('(^(.*)-(\d+)-(\d+)$)', name)
#     if match:
#         query = "UPDATE analogsignal SET index=%s WHERE name='%s'"%(match.groups()[3],name)
#         cursor.execute(query)
#         ndb.commit()
# 
# 
# pass