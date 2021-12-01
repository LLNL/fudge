from pqu import PQU


def convert_units_to_energy(T, energyUnits='eV'):
    if type(T) == str:
        T = PQU.PQU(T)
    if not isinstance(T, PQU.PQU):
        raise TypeError("argument to convert_units_to_energy must be a PQU, got type=%s" % str(type(T)))
    if T.isEnergy():
        return T.inUnitsOf(energyUnits)
    elif T.isTemperature():
        return (T.inUnitsOf('K') * PQU.PQU(8.6173324e-5, 'eV/K')).inUnitsOf(energyUnits)
    #    elif T.isMass(): return
    #    elif T.isLength(): return
    #    elif T.isTime(): return
    else:
        raise Exception("Cannot convert %s to energy units" % str(T))

