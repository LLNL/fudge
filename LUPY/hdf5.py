# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

"""
Use this module to import h5py.
If module not found, sets HDF5_present = False

Also guards against Python interpreter crash that may occur on if MPI environment variable is set before importing h5py.
"""

import os
PMI_FD = os.environ.get('PMI_FD')
if PMI_FD is not None:
    del os.environ['PMI_FD']

try:
    import h5py
    HDF5_present = True
except ImportError:
    h5py = None
    HDF5_present = False

if PMI_FD is not None:
    os.environ['PMI_FD'] = PMI_FD
