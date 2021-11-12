# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""Defines the supported version numbers for the GNDS format (not the same as the FUDGE version number)."""

version_1_10 = '1.10'
version_2_0_LLNL_4 = '2.0.LLNL_4'

allowed = ( version_1_10, version_2_0_LLNL_4 )

default = version_1_10      # version written by default in translation scripts
gnds2_0 = version_2_0_LLNL_4

MAJORVERSION, MINORVERSION = default.split( '.' )[:2]

