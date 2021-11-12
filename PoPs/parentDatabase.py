# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

from xData import link as linkModule

class parentDatabase( linkModule.link ):
    """
    This class is used to link a 'derived' PoPs database back to its 'parent' database.
    """

    moniker = 'parentDatabase'

    def __init__( self, link=None, root=None, path=None, label=None, relative=False, **attributes ):

        from PoPs import database as databaseModule     # avoid circular import by defining here
        if link is not None and not isinstance( link, databaseModule.database ):
            raise TypeError("parentDatabase must link to a PoPs.database instance, not '%s'" % type(link))

        linkModule.link.__init__( self, link, root, path, label, relative, **attributes )
