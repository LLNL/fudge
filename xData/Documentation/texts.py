# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
This module contains the GNDS documentation child node endfCompatible class.
"""

from .. import text as textModule

class Abstract( textModule.Text ) :
    """A class representing a GNDS documentation/abstract node."""

    moniker = 'abstract'

class Body( textModule.Text ) :
    """A class representing a GNDS documentation/body node."""

    moniker = 'body'

class EndfCompatible( textModule.Text ) :
    """A class representing a GNDS documentation/endfCompatible node."""

    moniker = 'endfCompatible'

class Title( textModule.Text ) :
    """A class representing a GNDS documentation/title node."""

    moniker = 'title'
