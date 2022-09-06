# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

NUC_XSEC_ADJUST_FOR_HEATED_TARGET_DIR = ../Src
CFLAGS = -g -O0 -Wall
local_CFLAGS = -I$(NUC_XSEC_ADJUST_FOR_HEATED_TARGET_DIR) $(CFLAGS)

TARGETS = generatePeaks test_threshold_xsec test_one_over_v_xsec test_constant_xsec test_xsec_integrate_xn_qauss

.PHONY: all check clean realclean

check: test_constant_xsec test_one_over_v_xsec
	./test_constant_xsec -check; if [ $$? -ne 0 ]; then echo 'ERROR in crossSectionAdjustForHeatedTarget: test_constant_xsec'; fi
	./test_one_over_v_xsec -check; if [ $$? -ne 0 ]; then echo 'ERROR in crossSectionAdjustForHeatedTarget: test_one_over_v_xsec'; fi

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
