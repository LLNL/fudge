# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from pqu import PQU as PQUModule

from fudge import sums as sumsModule

from . import gndsToENDF6 as gndsToENDF6Module


#
# summed cross section writing back to ENDF
#
def toENDF6( self, endfMFList, flags, targetInfo, verbosityIndent = '' ) :

    if self.ENDF_MT in range(851,871) : return   # for lumped-sum covariance, only write MF 33
    if 'omit' in targetInfo['ENDFconversionFlags'].get(self,""): return
    if flags['verbosity'] >= 10 : print( '%ssummed reaction: %s' % (verbosityIndent, self.label) )
    Q = gndsToENDF6Module.getForm( targetInfo['style'], self.Q )
    Q = PQUModule.PQU( Q.value, Q.axes[0].unit ).getValueAs( 'eV' )
    targetInfo['Q'] = Q                     # FIXME may need to account for ground state
    crossSection = self.crossSection
    targetInfo['EMin'], targetInfo['EMax'] = crossSection.domainMin, crossSection.domainMax
    if self.summands:
        level = min( [ max( [ product.getLevelAsFloat('eV') for product in reaction.link.ancestor.outputChannel ] )
            for reaction in self.summands ] )
    else :
        level = 0
    crossSection.toENDF6( self.ENDF_MT, endfMFList, targetInfo, level, LR=0 )

sumsModule.CrossSectionSum.toENDF6 = toENDF6

#
# summed multiplicities are written back in gndsToENDF6.gammasToENDF6_MF12_13. Just implement do-nothing routine here:
#
def toENDF6( self, endfMFList, flags, targetInfo, verbosityIndent = '' ) :
    pass

sumsModule.MultiplicitySum.toENDF6 = toENDF6
