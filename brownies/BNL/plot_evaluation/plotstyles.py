import os
import json
import fnmatch
import copy
from fudge.vis.matplotlib import defaultSymbols, nSymbols, defaultColors, nColors, defaultLineStyles, nLineStyles

# ---------------------------------------------------
# plot style resolution
# ---------------------------------------------------


DEFAULT_STYLE_DICT = {
    "plotStyles": {
        "crossSection": {
            "title": None,
            "showLegend": True,
            "legendX": 0.05,
            "legendY": 0.95,
            "referenceFrame": None,
            "xAxis": {
                "unit": "MeV",
                "min": 1e-11,
                "max": 20.0,
                "log": None,
                "label": "Incident energy"},
            "yAxis": {
                "unit": "b",
                "min": None,
                "max": None,
                "log": None,
                "label": "Cross section"}},
        "crossSectionUncertainty": {
            "title": None,
            "showLegend": True,
            "legendX": 0.05,
            "legendY": 0.95,
            "referenceFrame": None,
            "xAxis": {
                "unit": "MeV",
                "min": 1e-11,
                "max": 20.0,
                "log": None,
                "label": "Incident energy"},
            "yAxis": {
                "unit": "",
                "min": 0.0,
                "max": 5.0,
                "log": False,
                "label": "d(Cross section)/(Cross section)"}},
        "crossSectionIntegrals": {
            "title": None,
            "showLegend": True,
            "referenceFrame": None,
            "legendX": 0.05,
            "legendY": 0.95,
            "xAxis": {
                "unit": "MeV",
                "min": 1e-11,
                "max": 20.0,
                "log": True,
                "label": "Incident energy"},
            "yAxis": {
                "unit": "b",
                "min": 0.0,
                "max": 2.0,
                "log": False,
                "label": "C/E"}},
        "mubar": {
            "title": None,
            "showLegend": True,
            "legendX": 0.05,
            "legendY": 0.95,
            "referenceFrame": "lab",
            "xAxis": {
                "unit": "MeV",
                "min": 1e-11,
                "max": 20.0,
                "log": None,
                "label": "Incident energy"},
            "yAxis": {
                "unit": "",
                "min": -1.0,
                "max": 1.0,
                "log": False,
                "label": "<mu>"}},
        "momentumDeposit": {
            "title": None,
            "showLegend": False,
            "legendX": 0.05,
            "legendY": 0.95,
            "referenceFrame": "lab",
            "xAxis": {
                "unit": "MeV",
                "min": None,
                "max": None,
                "log": None,
                "label": "Incident energy"},
            "yAxis": {
                "unit": "MeV/c",
                "min": None,
                "max": None,
                "log": False,
                "label": "Outgoing forward momentum"}},
        "momentumBalance": {
            "title": None,
            "showLegend": False,
            "legendX": 0.05,
            "legendY": 0.95,
            "referenceFrame": "lab",
            "xAxis": {
                "unit": "MeV",
                "min": None,
                "max": None,
                "log": None,
                "label": "Incident energy"},
            "yAxis": {
                "unit": "",
                "min": None,
                "max": None,
                "log": False,
                "label": "Fraction of forward momentum"}},
        "fissionEnergyRelease": {
            "title": None,
            "showLegend": False,
            "legendX": 0.05,
            "legendY": 0.95,
            "referenceFrame": "lab",
            "xAxis": {
                "unit": "MeV",
                "min": None,
                "max": None,
                "log": None,
                "label": "Incident energy"},
            "yAxis": {
                "unit": "MeV",
                "min": 0.0,
                "max": 10.0,
                "log": False,
                "label": "Energy release"}},
        "energyDeposit": {
            "title": None,
            "showLegend": False,
            "legendX": 0.05,
            "legendY": 0.95,
            "referenceFrame": "lab",
            "xAxis": {
                "unit": "MeV",
                "min": None,
                "max": None,
                "log": None,
                "label": "Incident energy"},
            "yAxis": {
                "unit": "MeV",
                "min": 0.0,
                "max": 10.0,
                "log": False,
                "label": "Outgoing energy"}},
        "energyBalance": {
            "title": None,
            "showLegend": False,
            "legendX": 0.05,
            "legendY": 0.95,
            "referenceFrame": "lab",
            "xAxis": {
                "unit": "MeV",
                "min": None,
                "max": None,
                "log": None,
                "label": "Incident energy"},
            "yAxis": {
                "unit": "",
                "min": 0.0,
                "max": 2.0,
                "log": False,
                "label": "Fraction of available outgoing energy"}},
        "LegendreMoment": {
            "title": None,
            "showLegend": True,
            "legendX": 0.05,
            "legendY": 0.95,
            "referenceFrame": "centerOfMass",
            "xAxis": {
                "unit": "MeV",
                "min": 1e-11,
                "max": 20.0,
                "log": None,
                "label": "Incident energy"},
            "yAxis": {
                "unit": "",
                "min": -1.0,
                "max": 1.0,
                "log": False,
                "label": "Legendre moment"}},
        "energyTransfer": {
            "title": None,
            "showLegend": True,
            "legendX": 0.05,
            "legendY": 0.95,
            "referenceFrame": None,
            "xAxis": {
                "unit": "MeV",
                "min": None,
                "max": None,
                "log": None,
                "label": "Incident energy"},
            "yAxis": {
                "unit": "MeV",
                "min": None,
                "max": None,
                "log": False,
                "label": "Transferred energy"}},
        "formFactor": {
            "title": None,
            "showLegend": True,
            "legendX": 0.05,
            "legendY": 0.95,
            "referenceFrame": None,
            "xAxis": {
                "unit": "1/Ang",
                "min": None,
                "max": None,
                "log": None,
                "label": "electron recoil momentum, q"},
            "yAxis": {
                "unit": "",
                "min": None,
                "max": None,
                "log": False,
                "label": "Form factor"}},
        "anomolousScatteringFactor": {
            "title": None,
            "showLegend": True,
            "legendX": 0.05,
            "legendY": 0.95,
            "referenceFrame": None,
            "xAxis": {
                "unit": "MeV",
                "min": 1e-6,
                "max": 10.0,
                "log": None,
                "label": "Incident energy"},
            "yAxis": {
                "unit": "",
                "min": None,
                "max": None,
                "log": False,
                "label": "anomolous scattering factor"}},
        "nubar": {
            "title": None,
            "showLegend": True,
            "legendX": 0.05,
            "legendY": 0.95,
            "referenceFrame": None,
            "xAxis": {
                "unit": "MeV",
                "min": 1e-11,
                "max": 20.0,
                "log": None,
                "label": "Incident energy"},
            "yAxis": {
                "unit": "",
                "min": 0.0,
                "max": 5.0,
                "log": False,
                "label": "nubar"}}},
    "evaluation": {
        "ACEFile": {
            "filePattern": "*.ace.endf",
            "legend": "ACE, 293.6 K",
            "lineStyle": "-",
            "lineColor": "SandyBrown",
            "lineWidth": 1,
            "errorColor": "None"},
        #        "ENDFFile": {
        #            "filePattern":"*.endf",
        #            "legend":"ENDF, 0 K",
        #            "lineStyle":"-",
        #            "lineColor":"b",
        #            "lineWidth":3,
        #            "errorColor":"b" },
        "heatedPREPROFile": {
            "filePattern": "*.sigma1",
            "legend": "PREPRO, 293.6 K",
            "lineStyle": "-.",
            "lineColor": "k",
            "lineWidth": 2,
            "errorColor": "k"},
        "PREPROFile": {
            "filePattern": "*.recent",
            "legend": "PREPRO, 0 K",
            "lineStyle": "--",
            "lineColor": "k",
            "lineWidth": 2,
            "errorColor": "k"},
        "GNDSFile": {
            "filePattern": "*.gnds.xml",
            "legend": "GNDS, 0 K",
            "lineStyle": "-",
            "lineColor": "b",
            "lineWidth": 2,
            "errorColor": "b"},
        "AMPXFile": {
            "filePattern": "*.ampx*",
            "legend": "AMPX",
            "lineStyle": "-",
            "lineColor": "DarkGreen",
            "lineWidth": 1,
            "errorColor": "None"},
        "default": {
            "filePattern": "*",
            "legend": "ENDF, 0 K",
            "lineStyle": "-",
            "lineColor": "b",
            "lineWidth": 3,
            "errorColor": "b"}},
    "xyCurve": {
        "default": {
            "filePattern": "*.xy.dat",
            "legend": "XY Pairs",
            "symbol": "",
            "xUnit": "eV",
            "yUnit": "b",
            "dataType": "line",
            "lineStyle": "-",
            "lineColor": "k",
            "lineWidth": 1,
            "errorColor": "None"}},
    "xydyCurve": {
        "default": {
            "filePattern": "*.xydy.dat",
            "legend": "XYdY Triplets",
            "symbol": "",
            "xUnit": "eV",
            "yUnit": "b",
            "dataType": "line",
            "lineStyle": "-",
            "lineColor": "k",
            "lineWidth": 1,
            "errorColor": "k"}},
    "xdxydyCurve": {
        "default": {
            "filePattern": "*.xdxydy.dat",
            "legend": "XdXYdY Quadruplets",
            "symbol": "",
            "xUnit": "eV",
            "yUnit": "b",
            "dataType": "line",
            "lineStyle": "-",
            "lineColor": "k",
            "lineWidth": 1,
            "errorColor": "k"}},
    "C4Set": {
        "default": {
            "filePattern": "*",
            "legend": None,
            "symbol": "",
            "xUnit": "eV",
            "yUnit": "b",
            "dataType": "points",
            "lineStyle": "",
            "lineColor": "k",
            "lineWidth": 1,
            "errorColor": "k"}},
    "EXFORSet": {
        "default": {
            "filePattern": "*",
            "legend": None,
            "symbol": "",
            "xUnit": None,
            "yUnit": None,
            "dataType": "points",
            "lineStyle": "",
            "lineColor": "k",
            "lineWidth": 1,
            "errorColor": "k"}}}


# Plot symbol & color for cases where it is not set in the DEFAULT_STYLE_DICT
def getPlotSymbol(i):
    return defaultSymbols[i % nSymbols]


def getPlotColor(i, defaultScheme=True):
    if defaultScheme:
        return defaultColors[i % nColors]
    if type(i) == str and i.startswith('_') and i.endswith('_'):
        return '0.5'
    return '#' + str(hex(hash(i) % (256 * 256 * 256))).replace('0x', '').zfill(6)


def getLineStyle(i): return defaultLineStyles[i % nLineStyles]


# These define the allowed styling options
allowedObservables = list(DEFAULT_STYLE_DICT["plotStyles"].keys())
allowedSetKinds = ['xydyCurve', 'xyCurve', 'xdxydyCurve', 'C4Set', 'evaluation', 'EXFORSet']
allowedSetStyleKeys = ["filePattern", "legend", "symbol", "lineStyle", "lineColor", "lineWidth", "errorColor", "xUnit",
                       "yUnit", "dataType"]
allowedPlotStyleKeys = ["title", "showLegend", "legendX", "legendY", "referenceFrame", "xAxis", "yAxis"]
allowedAxisStyleKeys = ["unit", "min", "max", "log", "label"]


def getThisPlotStyle(userdefs=None, observable=None):
    """
    Use this to generate the plot style for plots of the observable
    
    userdefs: should be main one retrieved by __main__ routine
    observable: 
    """
    plotstyles = copy.copy(DEFAULT_STYLE_DICT["plotStyles"][observable])
    # Update, with checking
    if userdefs is not None and "plotStyles" in userdefs:
        if observable is not None and observable in userdefs["plotStyles"]:
            for axis in ['xAxis', 'yAxis']:
                if axis in userdefs["plotStyles"][observable]:
                    for k in userdefs["plotStyles"][observable][axis]:
                        if k in allowedAxisStyleKeys:
                            plotstyles[axis][k] = userdefs["plotStyles"][observable][axis][k]
                        else:
                            print("In axis '%s', ignored unknown user defined axisStyle key: '%s'" % (axis, k))
            for k in userdefs["plotStyles"][observable]:
                if k in ['xAxis', 'yAxis']:
                    continue  # skip axes, did them already
                if k in allowedPlotStyleKeys:
                    plotstyles[k] = userdefs["plotStyles"][observable][k]
                else:
                    print("Ignored unknown user defined plotStyle key: '%s'" % k)
    return plotstyles


def getThisSetStyle(userdefs, setkind, setname, verbose=False):
    """
    Use this to work out what style to use for a data set
    
    userdefs: should be main one retrieved by __main__ routine
    setkind: in allowedSetKinds list
    setname: name to key off of when attempting to resolve the style information.
             Usually the data file name, but could be the EXFOR/C4 entry.
    """
    if setkind not in allowedSetKinds:
        raise KeyError("Cannot resolve data kind '%s', so cannot look up style information" % setkind)

    # Set the global defaults first
    result = copy.copy(DEFAULT_STYLE_DICT[setkind]['default'])

    # Now search through for wildcard overrides in the default styles
    if verbose:
        print("\n...checking default styles...")
    for key in DEFAULT_STYLE_DICT[setkind]:
        if key == 'default':
            continue
        if verbose:
            print("-> does %s match filePattern '%s' of setkind '%s'/'%s'?" %
                  (setname, DEFAULT_STYLE_DICT[setkind][key]['filePattern'], setkind, key))
        if fnmatch.fnmatch(setname, DEFAULT_STYLE_DICT[setkind][key]['filePattern']):
            if verbose:
                print("    yes: %s matches filePattern '%s' of setkind '%s'/'%s'" %
                      (setname, DEFAULT_STYLE_DICT[setkind][key]['filePattern'], setkind, key))
            result.update(DEFAULT_STYLE_DICT[setkind][key])
        else:
            if verbose:
                print("    no...")

    # Final pass for user specific overrides
    if verbose:
        print("\n...checking user defined styles...")
    if setkind in userdefs:
        for key in userdefs[setkind]:
            if verbose:
                print("-> does %s match filePattern '%s' of setkind '%s'/'%s'?" %
                      (setname, userdefs[setkind][key]['filePattern'], setkind, key))
            if 'filePattern' in userdefs[setkind][key] and fnmatch.fnmatch(setname,
                                                                           userdefs[setkind][key]['filePattern']):
                if verbose:
                    print("    yes: %s matches filePattern '%s' of setkind '%s'/'%s'" %
                          (setname, userdefs[setkind][key]['filePattern'], setkind, key))
                for skey in userdefs[setkind][key]:
                    if skey not in allowedSetStyleKeys:
                        print("Ignored unknown user defined setStyle key: %s" % skey)
                result.update(userdefs[setkind][key])
            else:
                if verbose:
                    print("    no...")

    return result


def readUserStyles(styleFile):
    if styleFile is None:
        return {}
    if not os.path.exists(styleFile):
        raise IOError("Cannot file user defined style file '%s'" % styleFile)
    with open(styleFile) as jsonfile:
        return json.loads(jsonfile.read())
