# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

SMR_PATH = ../../../../statusMessageReporting
NUMERICALFUNCTIONS_PATH = ../../../..

ptwXYCFLAGS = -lnfut_utilities -lptwXY_utilities \
			  -I$(NUMERICALFUNCTIONS_PATH)/include -L$(NUMERICALFUNCTIONS_PATH)/lib -lnumericalfunctions \
              -I$(SMR_PATH)/include -L$(SMR_PATH)/lib -lstatusMessageReporting
localCFLAGS = $(CFLAGS) $(ptwXYCFLAGS)

TARGETS = xArray clip trim dullEdges mutualify mod areDomainsMutual mergeClosePoints copy clone ysMappedToXs

.PHONY: default check clean realclean $(TARGETS)

default: $(TARGETS)

xArray: xArray.c
	$(CC) -o xArray xArray.c $(localCFLAGS) -lm

clip: clip.c
	$(CC) -o clip clip.c $(localCFLAGS) -lm

trim: trim.c
	$(CC) -o trim trim.c $(localCFLAGS) -lm

dullEdges: dullEdges.c
	$(CC) -o dullEdges dullEdges.c $(localCFLAGS) -lm

mutualify: mutualify.c
	$(CC) -o mutualify mutualify.c $(localCFLAGS) -lm

mod: mod.c
	$(CC) -o mod mod.c $(localCFLAGS) -lm

areDomainsMutual: areDomainsMutual.c
	$(CC) -o areDomainsMutual areDomainsMutual.c $(localCFLAGS) -lm

mergeClosePoints: mergeClosePoints.c
	$(CC) -o mergeClosePoints mergeClosePoints.c $(localCFLAGS) -lm

copy: copy.c
	$(CC) -o copy copy.c $(localCFLAGS) -lm

clone: clone.c
	$(CC) -o clone clone.c $(localCFLAGS) -lm

ysMappedToXs: ysMappedToXs.c
	$(CC) -o ysMappedToXs ysMappedToXs.c $(localCFLAGS) -lm

check: default
	for f in $(TARGETS); do ./$$f $(CHECKOPTIONS); done

clean:
	rm -f v

realclean: clean
	rm -fr $(TARGETS) *.dSYM
