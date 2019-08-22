try:
    from ..consts import *
    from ..utils import *
except ValueError:
    try:
        from BNL.inter.consts import *
        from BNL.inter.utils import *
    except:
        from site_packages.BNL.inter.consts import *
        from site_packages.BNL.inter.utils import *
import json


spectra=json.load(open(os.sep.join([INTERDIR,'spectra','spectra.json']),mode='r'))

def get_realpath(path):
    realpath = path
    if not os.path.exists(realpath):
        realpath = os.sep.join([INTERDIR, 'spectra', path])
        if not os.path.exists(realpath):
            raise IOError("Cannot find spectrum file at %s or %s" % (path, realpath))
    return realpath

def get_spectrum(k):
    if spectra[k]['format'] == 'KENO':
        return get_KENO_flux( get_realpath(spectra[k]['path']) )
    elif spectra[k]['format'] == 'ENDF-6':
        path,mat,mf,mt=spectra[k]['path'].split('/')
        return get_ENDF_flux( get_realpath(path), mat, mf, mt )
    else: # free format
        return get_free_format_flux( spectra, k )

def get_KENO_flux( file ):
    energies=[]
    energyUnits='eV'
    fluxes=[]
    fluxUnits='1/eV'
    for line in open(file,mode='r').readlines()[6:305]:
        sline=line.split()
        energies.append(float(sline[1]))
        fluxes.append(float(sline[2]))
    energies.append(1e-5)
    energies.reverse()
    fluxes.reverse()
    return grouped_values_to_XYs( energies, fluxes, domainUnit=energyUnits, rangeUnit=fluxUnits )

def get_ENDF_flux( file, mat, mf, mt ):
    import fudge.legacy.converting.ENDFToGNDS.endfFileToGNDSMisc as endfFileToGNDSMiscModule
    import xData.axes as axesModule
    import xData.XYs as XYsModule
    # Grab the data corresponding to the spectrum
    MF3Data = [ line.strip('\r\n') for line in open(file).readlines() if
                int(line[67:70])==int(mat) and
                int(line[70:72])==int(mf) and
                int(line[72:75])==int(mt)]
   # Define the axes
    spectrumAxes= axesModule.axes( rank = 2 )
    spectrumAxes[0] = axesModule.axis( 'spectrum', 0, '1/eV' )
    spectrumAxes[1] = axesModule.axis( 'energy_in', 1, 'eV' )
    # Now make the spectrum as an XYs1d
    dataLine, TAB1, spectrumRegions = endfFileToGNDSMiscModule.getTAB1Regions(1, MF3Data, allowInterpolation6=True,
                                                                                 logFile=None, #info.logs,
                                                                                 axes=spectrumAxes,
                                                                                 cls=XYsModule.XYs1d)
    if( len( spectrumRegions ) == 1 ) :         # Store as XYs1d.
        spectrumForm = XYsModule.XYs1d( data = spectrumRegions[0], label = None,
                axes = spectrumAxes, interpolation = spectrumRegions[0].interpolation )
    else :
        raise NotImplementedError("Fluxes with multiple regions for %s/%s/%s/%s"%(file, mat, mf, mt))
        spectrumForm = crossSectionModule.regions1d( label = None, axes = spectrumAxes )
        for region in spectrumRegions :
            if( len( region ) > 1 ) : spectrumForm.append( region )
    return( spectrumForm )

def get_user_defined_flux( file, key ):
    if not os.path.exists( file ): raise IOError( "Cannot find spectrum file %s"%file )
    return get_free_format_flux( json.load(open(file, mode='r')), key )

def get_free_format_flux( spectraDict, key ):
    energyUnit=str(spectraDict[key]['xunit'])
    fluxUnit=str(spectraDict[key]['yunit'])
    minE=1e-5
    energyConversion={'eV':1,'keV':1e3,'MeV':1e6}[energyUnit]
    fluxConversion={'1/eV':1,'1/keV':1e-3,'1/MeV':1e-6}[fluxUnit]
    energies=[]
    fluxes=[]
    format=spectraDict[key]['format'].split()
    try:
        iy=format.index('y')
    except ValueError:
        raise KeyError("Format specifier does not give 'y' for the flux column")
    if 'x' in format:
        ix=format.index('x')
        for row in spectraDict[key]['data']:
            energies.append(row[ix]*energyConversion)
            fluxes.append(row[iy]*fluxConversion)
        if energies[0]>energies[1]:
            energies.append(minE)
            energies.reverse()
            fluxes.reverse()
        else:
            energies.insert(0,minE)
    elif 'xmin' in format and 'xmax' in format:
        ixmin=format.index('xmin')
        ixmax=format.index('xmax')
        firstEmin=None
        for row in spectraDict[key]['data']:
            enMin=row[ixmin]*energyConversion
            enMax=row[ixmax]*energyConversion
            if firstEmin is None: firstEmin=enMin
            if enMin not in energies: energies.append(enMin)
            if enMax not in energies: energies.append(enMax)
            fluxes.append(row[iy]*fluxConversion)
        if firstEmin>energies[-1]:fluxes.reverse()
        energies.sort()
        #print energies
        #print fluxes
    else:
        raise KeyError("Format specifier does not give 'x' or both of 'xmin' and 'xmax' for the energy column(s)")
    if energies[0] < minE: energies[0]=minE
    if spectraDict[key]['fluxIntegratedOverGroup']:
        for i in range(len(fluxes)):
            fluxes[i]/=abs(energies[i+1]-energies[i])
    return grouped_values_to_XYs( energies, fluxes, domainUnit='eV', rangeUnit='1/eV' )


