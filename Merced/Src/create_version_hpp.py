# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

import os,glob
boilerPlate="""
/*
 * ******** get_transfer: calculate the transfer matrix *********
 * 
 * Copyright (c) , The Regents of the University of California. 
 * All rights reserved.
 * 
 * Produced at the Lawrence Livermore National Laboratory. 
 * Written by Gerald Hedstrom
 * 
 * This file is part of  v1.0  (UCRL-CODE-)
 * 
 * Please read the COPYING file for "Our Notice and GNU General 
 * Public License" in the root of this software distribution.  
 * 
 * This program is free software; you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License (as published by 
 * the Free Software Foundation) version 2, dated June 1991. 
 * 
 * This program is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms 
 * and conditions of the GNU General Public License for more details. 
 * 
 * You should have received a copy of the GNU General Public License along 
 * with this program; if not, write to the Free Software Foundation, Inc., 
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 
 * 
 * ******** get_transfer: calculate the transfer matrix *********
 */
"""

version_dot_hpp=boilerPlate+'#ifndef __NDFGEN_VERSION_HPP\n#define __NDFGEN_VERSION_HPP\n\n#define NDFGEN_VERSION "1.0"\n#define NDFGEN_VERSION_DATE "PUTDATETIMEHERE"\n#define NDFGEN_REVISION "PUTREVISIONHERE"\n\n#endif\n'

flist = [ fileName for fileName in glob.glob( '*pp' ) if fileName != 'version.hpp' ]
last_datetime = None
last_revision = 0
if os.path.exists("version.hpp"):
    for line in open("version.hpp",mode='r').readlines():
        if line.find("FETE_REVISION")!=-1: 
            last_revision = int(line.strip().replace('"','').split()[-1])
first_revision = last_revision
for f in flist:
    lines = os.popen('ident '+f,'r').readlines()
    for line in lines:
        if 'Date' in line: 
            dmy = line[:-1].split()[1]
            time = line[:-1].split()[2]
            this_datetime = (dmy, time)
            last_datetime = max(last_datetime,this_datetime)
        if 'Revision' in line:
            this_revision = line[:-1].split()[1]
            last_revision = max(last_revision,this_revision)
if last_revision > first_revision:
    open("version.hpp",mode="w").writelines(version_dot_hpp.replace('PUTDATETIMEHERE'," ".join(last_datetime)).replace('PUTREVISIONHERE',last_revision))
