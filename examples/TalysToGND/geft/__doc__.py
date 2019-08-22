# <<BEGIN-copyright>>
# Copyright (c) 2011, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Computational Nuclear Physics group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-494171 All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
# <<END-copyright>>

"""
Module for interacting with TALYS to create evaluations.
Can run TALYS creating inputs and extracting the data from the output files and exporting evaluation into ENDL format.

Simple examples:

1) If have a TALYS calculation that you want to import from directory 'talysDir':
>>> import fudge,geft
>>> ZA = 36077
>>> talysDir = 'talysDir'
>>> t = geft.talysZA(fudge.endlZA(ZA,1),talysDir=talysDir)
>>> a = geft.CompleteEvaluation(t)
>>> a.toEndl()
>>> t.process()
>>> t.save()

2) Run TALYS in default mode and export to ENDL:
>>> import fudge,geft
>>> ZA = 36077
>>> t = geft.talysZA(fudge.endlZA(ZA,1))
>>> t.run()
>>> a = geft.CompleteEvaluation(t)
>>> a.toEndl()
>>> t.process()
>>> t.save()

To use the Audi-Wapster 2003 mass evaluation include the following 2 lines before you import fudge.
>>> import bdfls
>>> b = bdfls.getDefaultBdfls( template = '/usr/gapps/data/nuclear/bdfls.Audi_etal.2003.12.22' )

The basic idea with geft is you make a talysZA object, which corresponds to a endlZA object
(i.e. unique object for a projectile-target combination).
This object can be used to create talys inputs, and run talys.
Then you can create an evaluation using CompleteEvaluation and passing it the talysZA object.
This will contain a list of all the data that will be put into ENDL.
The data in each item is found in .data e.g.
>>> a
'prints out a complete list of the data found in a'
>>> a[0]
(1, (n,total) integrated)
>>> a[0].data
'prints out endl2dmath object of the total cross section'
>>> a[0].data.plot()
'will produce a gnuplot of the data'

The CompleteEvaluation is a endlList object. This has some functions to help manage the list.
Type help(endlList) for more details, e.g. to get all the integrated cross sections
>>> a.findDatas(I=0)
'returns a list of all integrated cross sections'
>>> fudge.multiPlot([d.data for d in a.findDatas(I=0)])
'plots all the integrated cross sections on one graph'
This actually works better once the data is exported to ENDL format and you can use the commands in fudge
>>> a.toEndl()
>>> fudge.multiPlot(t.findDatas(I=0))

To create your own input file for the talys run with custumized input parameters you need to pass the talysZA
object a talysInput object when you setup the talysZA object.
>>> i = geft.talysInput.Input()
You can add keywords to the default run using the command addKeyword()
>>> i.addKeyword( maxrot=4, spherical='y' )
Then before you run() talys you can use the setup() command and pass it the input file
(run() automatically calls setup() if not already called)
>>> t.setup( input=i, endf=True )
The input lines for the default header, e.g. project, target element, mass, and the energy list are automatically
generated along with the required output flags needed to create ENDL files.
The endf=True keyword is so that the energy list goes from 10e-11 upto 20 MeV.
"""

