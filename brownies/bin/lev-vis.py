#! /usr/bin/env python

import sys
import os
import argparse
import brownies.BNL.RIPL.level_density as LD
import brownies.BNL.RIPL.level_scheme as LS
from PoPs.chemicalElements .misc import symbolFromZ as elementSymbolFromZ
import numpy as np
import matplotlib.pyplot as plt
HOME=os.environ['HOME']
DESKTOP=HOME+'/Desktop'
ATLAS=DESKTOP+"/atlas/"
RIPL=DESKTOP+"/empire.trunk/RIPL/"


# --------------------------------------------------------------------
#   Command line
# --------------------------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser(description="Plot ENTIRE level schemes, as determined from RIPL and the Atlas")
    parser.add_argument('Z', type=int, help='Nucleus charge')
    parser.add_argument('A', type=int, help="Nuclear atomic number")
    parser.add_argument('-v', dest='verbose', default=False, action='store_true', help='Run verbosely')
    parser.add_argument('-q', dest='verbose', action='store_false', help='Run quietly')
    parser.add_argument('--RIPL', default=RIPL, help="Path to RIPL files")
    parser.add_argument('--ATLAS', default=ATLAS, help="Path to atlas project")
    parser.add_argument('--dEmax', type=float, default=2.0, help="Plot from Emin=0 to Emax=Esep+dEmax, dEmax in MeV")
    return parser.parse_args()


# -----------------------------------------------------------------------
# Set defaults
# -----------------------------------------------------------------------
args = parse_args()

Z, A=args.Z, args.A
elem=elementSymbolFromZ[Z]
sym='%s%i'%(elem, A)
symTarget='%s%i'%(elem, A-1)  # for n+target = isotope of interest

ld, D, levs=None, None, None

sys.path.append(args.ATLAS)
import atlas.io as aio

# -----------------------------------------------------------------------
# Get data
# -----------------------------------------------------------------------

# Get the level scheme for Z,A.  From RIPL-3.
fname=args.RIPL+"/levels/z%s.dat"%str(Z).zfill(3)
with open(fname, mode='r') as f:
    levMap=LS.readRIPLLevelScheme(f.read(), verbose=False)
levs=levMap[sym]
spinTarget = levMap[symTarget].levels[0].spin
parityTarget = levMap[symTarget].levels[0].parity
spin0=levMap[sym].levels[0].spin
parity0=levMap[sym].levels[0].parity

if args.verbose:
    print('target:', symTarget, spinTarget, parityTarget)
    print('compound:', sym, spin0, parity0)
    print(levs.name, levs.Sn)
    print(levs.levelReport())

# Get the HFB level density for Z,A.  From RIPL-3.
fname=args.RIPL+"/densities/total/level-densities-hfb/z%s.tab"%str(Z).zfill(3)
with open(fname, mode='r') as f:
    ld=LD.readHFBMLevelDensityTable(f.read(), Z, A, verbose=True)

# Get the mean level spacing for the resonance region for Z,A-1.
# The compound nucleus # for neutron + Z,A-1 is Z,A.  From RIPL-3.
# FIXME

# Get the resonances for Z,A-1.  The compound nucleus for neutron + Z,A-1 is Z,A.
# From the Atlas of Neutron Resonances, 6th edition.
try:
    res = aio.read_atlas(isotope=None, element=None, Z=Z, A=A-1, ZA=None, verbose=False)
except KeyError:
    res = None

if args.verbose and res is not None:
    for r in res.resonance_parameters:
        print('\t'.join([str(x) for x in r]))

icut = levs.lastLevelInCompleteScheme
Ecut = levs.levels[icut].energy.value
Esep = levs.Sn.value
Emin = 0.0
Emax = Esep+args.dEmax

# -----------------------------------------------------------------------
# Hi!
# -----------------------------------------------------------------------

for NAME in ["Angie", 'Nathaniel', "Mami"]:
    print('hi '+NAME)


# -----------------------------------------------------------------------
# Make plot
# -----------------------------------------------------------------------

# Set up axes and title
plt.title(levs.name)
plt.xlabel("$\Pi*J$")
plt.ylabel("$E^*$ (MeV)")

# Widget to get J & Pi, a common theme
def get_J_and_Pi(__lev, useNone=True):
    if __lev.spin is None:
        if useNone:
            J = None
        else:
            J = -11.33333
    else:
        J = float(lev.spin)
    if lev.parity is None:
        if useNone:
            Pi = None
        else:
            Pi = 1.0
    else:
        Pi = float(str(lev.parity))
    return J, Pi

# Plot discrete levels with complete information
if True:
    x, y, xerr = [], [], [] # for completely know levels
    for lev in levs.levels:
        J, Pi=get_J_and_Pi(lev, True)
        if J is None or Pi is None:
            pass
        else:
            if lev.energy.value > Emax:
                continue
            y.append(lev.energy.value)
            x.append(J*Pi)
            xerr.append(0.25)
    plt.errorbar(x=x, y=y, xerr=xerr, linewidth=0, elinewidth=2, color='k')

# Highlight ground state
if True:
    plt.errorbar(x=[float(spin0)*int(parity0)], y=[0.0], xerr=[0.25], linewidth=0, elinewidth=2, color='g')

# Highlight gamma transitions
if False:
    raise NotImplementedError("Highlight gamma transitions")

# Plot discrete levels missing either J or Pi
if True:
    x, y, xerr = [], [], [] # for completely know levels
    for lev in levs.levels:
        J, Pi=get_J_and_Pi(lev, True)
        if J is None or Pi is None:
            J, Pi=get_J_and_Pi(lev, False)
            if lev.energy.value > Emax:
                continue
            y.append(lev.energy.value)
            x.append(J*Pi)
            xerr.append(0.25)
        else:
            pass
    plt.errorbar(x=x, y=y, xerr=xerr, linewidth=0, elinewidth=2, color='blue')

# Highlight rotational bands
if False:
    raise NotImplementedError("Highlight rotational bands")

# Highlight vibrational "bands"
if False:
    raise NotImplementedError('Highlight vibrational "bands"')

# Plot Ecut
if True:
    plt.axhline(y=Ecut, color='black', alpha=0.25, linestyle=':')
    plt.text(11, Ecut+0.1, r'$E_{cut}$')

# Plot level density contour plot
if True:
    JPigrid = []
    for Pi in ld.spin_dep_level_density:
        for twoJ in ld.spin_dep_level_density[Pi].keys():
            if twoJ > 30: continue
            JPigrid.append(Pi*twoJ/2.0)
    JPigrid.sort()
    JPigrid = np.array(JPigrid)
    Exgrid = np.arange(Emin, Emax, 0.1) # np.arange(Ecut, 5.0+Esep, 0.1)
    X, Y = np.meshgrid(JPigrid, Exgrid)
    vevaluate = np.vectorize(ld.evaluate)
    Z = vevaluate(Y, np.abs(X), np.sign(X))
    CS = plt.contour(X, Y, Z, levels=[0, 1, 5]+[int(x) for x in np.logspace(1,3,14)])
    plt.clabel(CS, fontsize=9, inline=1)

# Plot Esep
if True:
    plt.axhline(y=Esep, color='black', alpha=0.25)
    plt.text(11, Esep+0.1, r'$E_{sep}$')

# Plot box of levels that can be excited by n+target reactions
if False:
    raise NotImplementedError("Plot box of levels that can be excited by n+target reactions")

# Plot resonances with known J
if True and res is not None:
    x, y, xerr, yerr = [], [], [], [] # for completely know levels
    for r in res.resonance_parameters:
        Er=r[0].value*1e-6
        if Er < 0.0:
            continue
        Er += Esep
        J, L = r[1].spin, r[2].spin
        if J is None or L is None:
            continue
        Pi = pow(-1, int(L)) * int(parityTarget)  # FIXME: check math!
        y.append(Er)
        x.append(float(J)*Pi)
        # FIXME: fuzzy band from width
        xerr.append(0.25)
        yerr.append(0.0)
    plt.errorbar(x=x, y=y, xerr=xerr, linewidth=0, elinewidth=2, color='red')

# Highlight primary gammas
if False:
    raise NotImplementedError("Highlight primary gammas")

# Plot resonances with unknown J
if False:
    raise NotImplementedError("Plot resonances with unknown J")

plt.show()
