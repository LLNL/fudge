from pqu import PQU
import os.path

ROOMTEMPERATURE    = PQU.PQU(20.,"degC")
ZERODEGREESCELCIUS = PQU.PQU(0.,"degC")
CADMIUMCUTTOFF     = PQU.PQU(0.5, 'eV')
SPEEDOFLIGHT       = PQU.PQU(299792458., 'm / s')
AVAGADROSNUMBER    = PQU.PQU(6.0221415e23,'1/mol')
INTERDIR = os.path.dirname(__file__)

# switches to turn on not-yet-implemented things
TURNONNEUTRONNSOURCES=False
TURNONFUNDIPPE=False

