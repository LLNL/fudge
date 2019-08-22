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

NUC_XSEC_ADJUST_FOR_HEATED_TARGET_DIR = ../Src
CFLAGS = -g -O0 -Wall
local_CFLAGS = -I$(NUC_XSEC_ADJUST_FOR_HEATED_TARGET_DIR) $(CFLAGS)

TARGETS = generatePeaks test_threshold_xsec test_one_over_v_xsec test_constant_xsec test_xsec_integrate_xn_qauss

.PHONY: all check clean realclean

check: test_constant_xsec test_one_over_v_xsec
	test_constant_xsec -check; if [ $$? -ne 0 ]; then echo 'ERROR in crossSectionAdjustForHeatedTarget: test_constant_xsec'; fi
	test_one_over_v_xsec -check; if [ $$? -ne 0 ]; then echo 'ERROR in crossSectionAdjustForHeatedTarget: test_one_over_v_xsec'; fi

all: $(TARGETS)

generatePeaks: generatePeaks.c
	$(CC) $(local_CFLAGS) -o generatePeaks generatePeaks.c -lm

test_threshold_xsec: test_threshold_xsec.c $(NUC_XSEC_ADJUST_FOR_HEATED_TARGET_DIR)/crossSectionAdjustForHeatedTarget.o
	$(CC) $(local_CFLAGS) -o test_threshold_xsec test_threshold_xsec.c $(NUC_XSEC_ADJUST_FOR_HEATED_TARGET_DIR)/crossSectionAdjustForHeatedTarget.o -lm

test_one_over_v_xsec: test_one_over_v_xsec.c $(NUC_XSEC_ADJUST_FOR_HEATED_TARGET_DIR)/crossSectionAdjustForHeatedTarget.o
	$(CC) $(local_CFLAGS) -o test_one_over_v_xsec test_one_over_v_xsec.c $(NUC_XSEC_ADJUST_FOR_HEATED_TARGET_DIR)/crossSectionAdjustForHeatedTarget.o -lm

test_constant_xsec: test_constant_xsec.c $(NUC_XSEC_ADJUST_FOR_HEATED_TARGET_DIR)/crossSectionAdjustForHeatedTarget.o
	$(CC) $(local_CFLAGS) -o test_constant_xsec test_constant_xsec.c $(NUC_XSEC_ADJUST_FOR_HEATED_TARGET_DIR)/crossSectionAdjustForHeatedTarget.o -lm

test_xsec_integrate_xn_qauss: test_xsec_integrate_xn_qauss.c $(NUC_XSEC_ADJUST_FOR_HEATED_TARGET_DIR)/crossSectionAdjustForHeatedTarget.o
	$(CC) $(local_CFLAGS) -o test_xsec_integrate_xn_qauss test_xsec_integrate_xn_qauss.c $(NUC_XSEC_ADJUST_FOR_HEATED_TARGET_DIR)/crossSectionAdjustForHeatedTarget.o -lm

$(NUC_XSEC_ADJUST_FOR_HEATED_TARGET_DIR)/crossSectionAdjustForHeatedTarget.o:
	cd $(NUC_XSEC_ADJUST_FOR_HEATED_TARGET_DIR); $(MAKE)

clean:
	rm -f *.o $(TARGETS)

realclean: clean
	rm -f test_constant_xsec.data
