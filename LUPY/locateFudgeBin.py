# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os
import sys

import fudge as fudgeModule

'''
    Methods to find scripts or executables in the FUDGE bin folder since the location is dependent on teh way FUDGE is installed.
'''

def locateFileInBin( fileInBin, exceptionIfNotFound=False ):
    '''
       Find the full path to a file located in the FUDGE bin folder

       :param fileInBin:           The file for which the path is seeked.
       :param exceptionIfNotFound: Option to raise a FileNotFoundError exception if the requested path is not found. 
       :return:                    Full path to the requested file path. 
    '''

    filePath = os.path.join(os.path.abspath('./'), fileInBin)
    if not os.path.exists( filePath ):
        filePath = os.path.join(sys.prefix, 'bin', fileInBin)
        if not os.path.exists( filePath ):
            filePath = os.path.join(os.path.split(os.path.dirname(fudgeModule.__file__))[0], 'bin', fileInBin)
            if not os.path.exists(filePath):
                if exceptionIfNotFound:
                    raise FileNotFoundError(f'FUDGE bin folder file "{fileInBin}" not found.')
                else:
                    filePath = None

    return filePath

def locateMerced( exceptionIfNotFound=False ):
    '''
        Find the full path to the mercury executable.

       :param exceptionIfNotFound: Option to raise a FileNotFoundError exception if no "mercury" executable is found.
       :return:                    Full path of the "merced" executable. 
    '''

    return locateFileInBin( 'merced', exceptionIfNotFound )

