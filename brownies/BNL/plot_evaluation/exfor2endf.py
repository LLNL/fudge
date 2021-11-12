from sqlalchemy import Column, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine
import os

'''
COLUMNS   DEFINITION                                                           1
=======   ==========                                                           1
  1- 48   EXFOR REACTION (WARNING...DO NOT GO BEYOND COLUMN 48).               1
 49- 53   PROJECTILE ZA (E.G. NEUTRON = 1, PROTON = 1001)                      1
 54- 56   MF NUMBER (USE ENDF/B CONVENTION FOR MF =1 TO 99, MF=100 TO 999      1
          TO TRANSLATE DATA NOT EQUIVALENT TO ENDF/B)                          1
 57- 60   MT NUMBER (USE ENDF/B CONVENTION FOR MT =1 TO 999, MT=1000 TO 9999   1
          TO TRANSLATE DATA NOT EQUIVALENT TO ENDF/B)                          1
 61- 63   OPERATION NUMBER (SEE, LIST BELOW FOR DEFINITIONS).                  1
'''

# Define the tables
Base = declarative_base()


class EXFOR2ENDF(Base):
    __tablename__ = 'MFMTMapper'
    id = Column(Integer, primary_key=True)
    exfor_reaction = Column(String(48))
    projectile_za = Column(Integer)
    endf_mf = Column(Integer)
    endf_mt = Column(Integer)
    operation_number = Column(Integer)


# Turn on a database
engine = create_engine('sqlite:///:memory:')
Base.metadata.create_all(engine)

# Make a session
from sqlalchemy.orm import sessionmaker

session = sessionmaker()
session.configure(bind=engine)


# ------------------------------------------------
# database loading utility
# ------------------------------------------------
def load_database(data_path=os.path.split(os.path.abspath(__file__))[0], verbose=False):
    """
    Database loading utility
    :param data_path: path to EXFOR14A.DAT file, defaults to the directory containing this file
    :param verbose: toggle verbose output
    :return: None
    """
    s = session()
    for line in open(data_path + os.sep + 'EXFOR14A.DAT').readlines():
        if len(line) > 78 and line[79] == '1':
            continue
        if verbose:
            print(line)
            print('  rxn:', '"' + line[:48] + '"')
            print('   za:', '"' + line[48:53] + '"')
            print('   mf:', '"' + line[53:56] + '"')
            print('   mt:', '"' + line[56:60] + '"')
            print('   op:', '"' + line[60:63] + '"')
        s.add(
            EXFOR2ENDF(
                exfor_reaction=line[:48].strip(),
                projectile_za=int(line[48:53]),
                endf_mf=int(line[53:56]),
                endf_mt=int(line[56:60]),
                operation_number=int(line[60:63])))
    s.commit()


load_database()


# ------------------------------------------------
# queries
# ------------------------------------------------
def get_exfor_reactions(za, mf=None, mt=None):
    s = session()
    if mf is None and mt is None:
        out = s.query(EXFOR2ENDF).filter_by(projectile_za=za)
    elif mf is None:
        out = s.query(EXFOR2ENDF).filter_by(projectile_za=za, endf_mt=mt)
    elif mt is None:
        out = s.query(EXFOR2ENDF).filter_by(projectile_za=za, endf_mf=mf)
    else:
        out = s.query(EXFOR2ENDF).filter_by(projectile_za=za, endf_mf=mf, endf_mt=mt)
    return [x.exfor_reaction for x in out]


def get_endf_mf_mt(reaction):
    s = session()
    out = s.query(EXFOR2ENDF).filter_by(exfor_reaction=reaction)
    return [(x.endf_mf, x.endf_mt) for x in out]

