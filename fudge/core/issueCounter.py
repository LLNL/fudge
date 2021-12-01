# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This modules contains the classes **IssueCounter** and **IssueCounters**.
"""

class IssueCounter :
    """
    This class is used to aid in keeping track of how many times a specific issue occurs and if a message should be
    logged for the issue. An issue can be an informational, warning or error.
    """

    def __init__( self, name, maximumReportsToLog = 3 ) :
        """
        Constructor for the **IssueCounter** class.

        :param name:                    The name of the issue. Must be unique within a list of *IssueCounter*.
        :param maximumReportsToLog:     The first *maximumReportsToLog* should be logged
        """

        self.__name = name
        self.__maximumReportsToLog = maximumReportsToLog
        self.__counter = 0

    @property
    def name( self ) :
        """Returns self's name."""

        return( self.__name )

    @property
    def maximumReportsToLog( self ) :
        """Returns self's maximumReportsToLog."""

        return( self.__maximumReportsToLog )

    @property
    def counter( self ) :
        """Returns self's counter."""

        return( self.__counter )

    def increment( self ) :
        """Increments the counter and returns the status state prior to the incrementation."""

        status = self.status( )
        self.__counter += 1

        return( status )

    def reset( self ) :
        """Reset the counter for *self* to 0."""

        self.__counter = 0

    def status( self ) :
        """Returns *True* if *counter* less than *maximumReportsToLog* and *False* otherwise."""

        return( self.__counter < self.__maximumReportsToLog )

class IssueCounters :
    """
    Contains a list of *IssueCounter* instances.
    """

    def __init__( self ) :
        """
        Constructor for the **IssueCounters** class.
        """

        self.__issueCounters = {}

    def __getitem__( self, name ) :
        """
        Returns the **IssueCounter** instance named *name*.

        :param name:                The name of the issue to return.
        """

        return( self.__issueCounters[name] )

    def __contains__( self, name ) :
        """
        Returns *True* if an **IssueCounter** instance named *name* is in self and *False* otherwise.

        :param name:                The name of the issue to check for.
        """

        return( name in self.__issueCounters )

    def add( self, issueCounter, overWrite = False ) :
        """
        Added *issueCounter* to self. If an **IssueCounter** with the same name exist and *overWrite* is *False*,
        *issueCounter* is not added. Otherwise, *issueCounter* is added.

        :param issueCounter:        The **IssueCounter** instance to add.
        :param overWrite:           Determines action to be taken if an instance with the same name as *issueCounter* is already in self.
        """

        if( not( isinstance( issueCounter, IssueCounter ) ) ) : raise TypeError( 'Invalid "IssueCounter" instance.' )
        if( ( issueCounter.name in self ) and not( overWrite ) ) : return
        self.__issueCounters[issueCounter.name] = issueCounter

    def get( self, name, maximumReportsToLog, add = False ) :
        """
        Returns the **IssueCounter** instance named *name* if it exists. If the named **IssueCounter** instance does not exist
        returns a newly constructed **IssueCounter** instance with parameters *name* and *maximumReportsToLog*. If *name* does not exist
        and *add* is *True*, the created **IssueCounter** instance is added to self.

        :param name:                    The name of the issue to create if one does not exist.
        :param maximumReportsToLog:     The *maximumReportsToLog* argument passed to a created **IssueCounter** instance if one does not exist.
        :param add:                     If true and *name* does not exist, adds the create **IssueCounter** instance to self.
        """

        try :
            return( self.__issueCounters[name] )
        except :
            pass

        issueCounter = IssueCounter( name, maximumReportsToLog )
        if( add ) : self.add( issueCounter )
        return( issueCounter )

    def keys( self ) :
        """
        Returns the list of *name*s in self.
        """

        return( self.__issueCounters.keys( ) )
