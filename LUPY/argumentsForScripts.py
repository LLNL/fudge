# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains common argparse options used by some of the FUDGE scripts.
"""

from argparse import ArgumentParser

from PoPs import IDs as IDsModule

from fudge import map as mapModule
from fudge import GNDS_file as GNDS_fileModule
from fudge import reactionSuite as reactionSuiteModule


class SingleProtareArguments :

    def __init__( self, parser ) :
        """
        :param parser:      An argparse.ArgumentParser instance which will be updated with parameters for accessing a profile from a file.
        """

        self.__parser = parser

        parser.add_argument( 'mapOrProtareFileName',                                    help = 'The path to a map or protare file. If it is a map file then option "--tid" (and maybe "--pid") is required.' )
        parser.add_argument( '--pid', action = 'store', default = IDsModule.neutron,    help = 'The PoPs id for the projectile. Only used if the argument "mapOrProtareFileName" is a map file; otherwise, it is ignored. Default is "%s".' % IDsModule.neutron )
        parser.add_argument( '--tid', action = 'store', default = None,                 help = 'The PoPs id for the target. Required if the argument "mapOrProtareFileName" is a map file; otherwise, it is ignored.' )
        parser.add_argument( '--noLazyParse', action = 'store_true',                    help = 'Disable lazy parsing (increases load time but sometimes useful).' )

    def protare(self, args, verbosity=0, lazyParsing=True):
        """
        Returns a protare (i.e., "reactionSuite") instance that has been read per the "mapOrProtareFileName", "--pid" and "--tid" parameters.
        """

        lazyParsing = lazyParsing and not args.noLazyParse

        mapOrProtare = GNDS_fileModule.read(args.mapOrProtareFileName, verbosity=verbosity, lazyParsing=lazyParsing)

        if( isinstance( mapOrProtare, reactionSuiteModule.ReactionSuite ) ) :
            return( mapOrProtare )
        elif( isinstance( mapOrProtare, mapModule.Map ) ) :
            if( args.tid is None ) : raise ValueError( '--tid option required for the map file.' )
            mapProtare = mapOrProtare.find( args.pid, args.tid )
            if( mapProtare is None ) : raise Exception( 'No match in map file for pid = "%s" and tid = "%s".' % ( args.pid, args.tid ) )
            return(mapProtare.protare(verbosity=verbosity, lazyParsing=lazyParsing))
        else :
            raise TypeError( 'File "%s" is not a map or protare file.' )
