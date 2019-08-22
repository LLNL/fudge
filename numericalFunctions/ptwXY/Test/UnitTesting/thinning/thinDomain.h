/*
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
*/

double test1_input[] = {
  1.00000000000000000e+00,   8.41470984807896505e-01,
  1.00000000006000000e+00,   8.41470984840314684e-01,
  1.00000000012000001e+00,   8.41470984872732752e-01,
  1.00000000018000001e+00,   8.41470984905150932e-01,
  1.00000000024000002e+00,   8.41470984937569111e-01,
  1.00000000030000002e+00,   8.41470984969987179e-01,
  1.00000000036000003e+00,   8.41470985002405358e-01,
  1.00000000042000003e+00,   8.41470985034823538e-01,
  1.00000000048000004e+00,   8.41470985067241606e-01,
  1.00000000054000004e+00,   8.41470985099659785e-01 };
int n_test1_input = sizeof( test1_input ) / ( 2 * sizeof( test1_input[0] ) );

double test1_1output[] = {
  1.00000000000000000e+00,   8.41470984807896505e-01,
  1.00000000012000001e+00,   8.41470984872732752e-01,
  1.00000000024000002e+00,   8.41470984937569111e-01,
  1.00000000036000003e+00,   8.41470985002405358e-01,
  1.00000000054000004e+00,   8.41470985099659785e-01 };
int n_test1_1output = sizeof( test1_1output ) / ( 2 * sizeof( test1_1output[0] ) );

double test1_2output[] = {
  1.00000000000000000e+00,   8.41470984807896505e-01,
  1.00000000006000000e+00,   8.41470984840314684e-01,
  1.00000000012000001e+00,   8.41470984872732752e-01,
  1.00000000018000001e+00,   8.41470984905150932e-01,
  1.00000000024000002e+00,   8.41470984937569111e-01,
  1.00000000030000002e+00,   8.41470984969987179e-01,
  1.00000000036000003e+00,   8.41470985002405358e-01,
  1.00000000042000003e+00,   8.41470985034823538e-01,
  1.00000000048000004e+00,   8.41470985067241606e-01,
  1.00000000054000004e+00,   8.41470985099659785e-01 };
int n_test1_2output = sizeof( test1_2output ) / ( 2 * sizeof( test1_2output[0] ) );

double test1_3output[] = {
  1.00000000000000000e+00,   8.41470984807896505e-01,
  1.00000000054000004e+00,   8.41470985099659785e-01 };
int n_test1_3output = sizeof( test1_3output ) / ( 2 * sizeof( test1_3output[0] ) );

double test2_input[] = {
  1.00000000000000000e+10,   5.00000000000000000e-01,
  1.00000000005000000e+10,   1.00000000000000000e+00,
  1.00000000100000000e+10,   1.50000000000000000e+00 };
int n_test2_input = sizeof( test2_input ) / ( 2 * sizeof( test2_input[0] ) );

double test2_1output[] = {
  1.00000000000000000e+10,   5.00000000000000000e-01,
  1.00000000010000000e+10,   1.02631578947368429e+00,
  1.00000000100000000e+10,   1.50000000000000000e+00 };
int n_test2_1output = sizeof( test2_1output ) / ( 2 * sizeof( test2_1output[0] ) );

double test2_2output[] = {
  1.00000000000000000e+10,   5.00000000000000000e-01,
  1.00000000005000000e+10,   1.00000000000000000e+00,
  1.00000000100000000e+10,   1.50000000000000000e+00 };
int n_test2_2output = sizeof( test2_2output ) / ( 2 * sizeof( test2_2output[0] ) );

double test2_3output[] = {
  1.00000000000000000e+10,   5.00000000000000000e-01,
  1.00000000030000000e+10,   1.13157894736842102e+00,
  1.00000000100000000e+10,   1.50000000000000000e+00 };
int n_test2_3output = sizeof( test2_3output ) / ( 2 * sizeof( test2_3output[0] ) );

double test3_input[] = {
  1.00000000000000000e+10,   5.00000000000000000e-01,
  1.00000000095000000e+10,   1.00000000000000000e+00,
  1.00000000100000000e+10,   1.50000000000000000e+00 };
int n_test3_input = sizeof( test3_input ) / ( 2 * sizeof( test3_input[0] ) );

double test3_1output[] = {
  1.00000000000000000e+10,   5.00000000000000000e-01,
  1.00000000090000000e+10,   9.73684210526315819e-01,
  1.00000000100000000e+10,   1.50000000000000000e+00 };
int n_test3_1output = sizeof( test3_1output ) / ( 2 * sizeof( test3_1output[0] ) );

double test3_2output[] = {
  1.00000000000000000e+10,   5.00000000000000000e-01,
  1.00000000095000000e+10,   1.00000000000000000e+00,
  1.00000000100000000e+10,   1.50000000000000000e+00 };
int n_test3_2output = sizeof( test3_2output ) / ( 2 * sizeof( test3_2output[0] ) );

double test3_3output[] = {
  1.00000000000000000e+10,   5.00000000000000000e-01,
  1.00000000070000000e+10,   8.68421052631578982e-01,
  1.00000000100000000e+10,   1.50000000000000000e+00 };
int n_test3_3output = sizeof( test3_3output ) / ( 2 * sizeof( test3_3output[0] ) );
