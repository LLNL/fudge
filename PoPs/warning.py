# <<BEGIN-copyright>>
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Nuclear Data and Theory group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-683960.
# All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the disclaimer below.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the disclaimer (as noted below) in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of LLNS/LLNL nor the names of its contributors may be used
#       to endorse or promote products derived from this software without specific
#       prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
# Additional BSD Notice
# 
# 1. This notice is required to be provided under our contract with the U.S.
# Department of Energy (DOE). This work was produced at Lawrence Livermore
# National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
# 
# 2. Neither the United States Government nor Lawrence Livermore National Security,
# LLC nor any of their employees, makes any warranty, express or implied, or assumes
# any liability or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that its use
# would not infringe privately-owned rights.
# 
# 3. Also, reference herein to any specific commercial products, process, or services
# by trade name, trademark, manufacturer or otherwise does not necessarily constitute
# or imply its endorsement, recommendation, or favoring by the United States Government
# or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
# herein do not necessarily state or reflect those of the United States Government or
# Lawrence Livermore National Security, LLC, and shall not be used for advertising or
# product endorsement purposes.
# 
# <<END-copyright>>
"""
Store and report warnings and errors in a PoPs database.
PoPs.check() returns a nested list of warning objects:

    >>> warnings = PoPs.check()
    >>> print warnings

May include or exclude specific classes of warning using the filter command.
filter() returns a new context instance:

    >>> warnings2 = warnings.filter( exclude=[warning.unnormalizedGammas] )

Or, for easier searching you may wish to flatten the list (to get warnings alone without context messages):

    >>> flat = warnings.flatten()
"""

# FIXME context class and base warning class are both identical to stuff in fudge.gnd.warning. Move to external utility?
__metaclass__ = type


class context:
    """
    Store warnings in context. This class contains location information (reactionSuite, reaction, etc)
    plus a nested list of warnings or other context instances
    """

    def __init__( self, message='', warningList=None ):
        self.message = message
        self.warningList = warningList or []

    def __len__( self ):
        return len(self.warningList)

    def __getitem__( self, idx ):
        return self.warningList[idx]

    def __str__( self ):
        if len(self.warningList) == 0:
            return self.message + ": no problems encountered"
        return '\n'.join(self.toStringList())

    def __eq__( self, other ):
        return self.message == other.message and self.warningList == other.warningList

    def filter( self, include=None, exclude=None ):
        """Filter warning list to only include (or exclude) specific classes of warning. For example:

            >>> newWarnings = warnings.filter( exclude=[warning.discreteLevelsOutOfOrder] )

        Note that if both 'include' and 'exclude' lists are provided, exclude is ignored."""

        if include is None and exclude is None: return self
        newWarningList = []
        for warning in self.warningList:
            if isinstance(warning, context):
                newContext = warning.filter(include, exclude)
                if newContext: newWarningList.append(newContext)
            elif include is not None:
                if warning.__class__ in include:
                    newWarningList.append(warning)
            else:  # exclude is not None:
                if warning.__class__ not in exclude:
                    newWarningList.append(warning)
        return context(self.message, newWarningList)

    def flatten( self ):
        """From a nested hierarchy of warnings, get back a flat list for easier searching:

            >>> w = PoPs.check()
            >>> warningList = w.flatten()"""

        List = []
        for val in self.warningList:
            if isinstance(val, warning):
                List.append(val)
            else:
                List += val.flatten()
        return List

    def toStringList( self, indent='', dIndent='    ' ):
        """ Format warnings for printing. Returns a list of warning strings with indentation. """
        s = ['%s%s' % (indent, self.message)]
        for warning in self.warningList:
            s += warning.toStringList(indent + dIndent)
        return s


class warning:
    """
    General warning class. Contains link to problem object,
    xpath in case the object leaves memory,
    and information about the warning or error.
    """

    def __init__( self, obj=None ):
        self.obj = obj
        self.xpath = ''
        if hasattr(obj, 'toXLink'):
            self.xpath = obj.toXLink()

    def __str__( self ):
        return "Generic warning for %s" % self.xpath

    def __eq__( self, other ):
        return self.xpath == other.xpath

    def toStringList( self, indent='' ):
        return ['%sWARNING: %s' % (indent, self)]


#
# specific warning classes:
#

class NotImplemented(warning):
    def __init__( self, form, obj=None ):
        warning.__init__(self, obj)
        self.form = form

    def __str__( self ):
        return "Checking not yet implemented for %s type data" % self.form

    def __eq__( self, other ):
        return (self.form == other.form and self.xpath == other.xpath)

class discreteLevelsOutOfOrder( warning ):
    def __init__(self, lidx, obj=None):
        warning.__init__(self, obj)
        self.lidx = lidx

    def __str__(self):
        return "Discrete level %s is out of order" % self.lidx

    def __eq__(self, other):
        return (self.lidx == other.lidx)

class unnormalizedDecayProbabilities( warning ):
    def __init__(self, branchingSum, obj=None):
        warning.__init__(self, obj)
        self.branchingSum = branchingSum

    def __str__(self):
        return "Sum of decay probabilities = %s, should be 1.0!" % (self.branchingSum)

    def __eq__(self, other):
        return (self.xpath == other.xpath and self.branchingSum == other.branchingSum)

class AliasToNonExistentParticle(warning):
    def __init__( self, id, pid, obj=None ):
        warning.__init__(self, obj)
        self.id = id
        self.pid = pid

    def __str__(self):
        return "Alias '%s' points to non-existant particle '%s'" % (self.id, self.pid)

    def __eq__(self, other):
        return (self.id == other.id and self.pid == other.pid)