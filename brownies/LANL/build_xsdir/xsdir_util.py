#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os

def xsdir_entry(ace_file, relative=True):
    """
    Generate xsdir entry for specified ACE file. Use absolute path if relative == False.

    :return: string xsdir entry
    """
    fileName = ace_file
    if not relative:
        fileName = os.path.abspath(fileName)

    with open(ace_file) as fin:
        header = [next(fin) for idx in range(12)]

        offset = 0
        if header[0].startswith('2.'):
            # 'new-style' ACE file with extra header info. Read past the extra stuff
            offset = int(header[1].split()[-1])
            header += [next(fin) for idx in range(offset)]
            extraHeader = header[:offset]
            header = header[offset:]

    tableName = header[0][:10]
    atomicMassRatio = header[0][10:22]
    ZAID = header[1].split( )[0]
    if '_e' in ZAID:
        raise NotImplementedError('Meta-stable targets not currently supported!')

    ZA, evaluationID = tableName.strip().split('.')
    aceType = evaluationID[-1]     # e.g. 'c' for neutrons, 't' for tnsl, 'h' for protons, etc.

    ptable = ""
    if aceType == "c":
        urrData = int(header[10].split()[6])
        if urrData:
            ptable = " +\n          ptable"

    accessRoute = '0'
    fileType = 1
    address = offset + 1
    tableLength = int(header[6][:9])
    recordLength = 0
    numberOfEntriesPerRecord = 0
    temperature = header[0][22:34]
    xsdir_str = "%10s %s %s %s %d %d %d %d %d %s%s" % (
        tableName, atomicMassRatio, fileName, accessRoute, fileType, address, tableLength,
        recordLength, numberOfEntriesPerRecord, temperature, ptable)

    return xsdir_str

if __name__ == '__main__':
    import sys
    for aceFile in sys.argv[1:]:
        print(xsdir_entry(aceFile))
