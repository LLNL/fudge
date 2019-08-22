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

flist = filter(lambda x:x!='version.hpp', glob.glob('*pp'))
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
