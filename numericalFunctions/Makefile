# <<BEGIN-copyright>>
# <<END-copyright>>

SHELL = /bin/ksh

.PHONY: default purgeIncludeAndLib includeAndLib ptwC ptwX ptwXY check checke clean realclean doTarget

SOURCE_DIRS = ptwC ptwX ptwXY nf_Legendre nf_specialFunctions
SUB_DIRS = include lib $(SOURCE_DIRS) Doc

default: purgeIncludeAndLib $(SOURCE_DIRS)

purgeIncludeAndLib:
	cd include; $(MAKE) realclean
	cd lib; $(MAKE) realclean

includeAndLib:
	cd include; $(MAKE) NFUPDATING_PATH=$(NFUPDATING_PATH)
	cd lib; $(MAKE) NFUPDATING_PATH=$(NFUPDATING_PATH)

ptwC:
	cd ptwC; $(MAKE)
	$(MAKE) includeAndLib NFUPDATING_PATH=../ptwC

ptwX: ptwC
	cd ptwX; $(MAKE)
	$(MAKE) includeAndLib NFUPDATING_PATH=../ptwX

ptwXY: ptwX
	cd ptwXY; $(MAKE)
	$(MAKE) includeAndLib NFUPDATING_PATH=../ptwXY

nf_Legendre: ptwX ptwXY
	cd nf_Legendre; $(MAKE)
	$(MAKE) includeAndLib NFUPDATING_PATH=../nf_Legendre

nf_specialFunctions: ptwC
	cd nf_specialFunctions ; $(MAKE)
	$(MAKE) includeAndLib NFUPDATING_PATH=../nf_specialFunctions

python:
	python setup.py --quiet build

check:
	$(MAKE) doTarget MODE=check SUB_DIRS="$(SOURCE_DIRS)"

checke:
	$(MAKE) doTarget MODE=check SUB_DIRS="$(SOURCE_DIRS)" CHECKOPTIONS=-e

clean:
	rm -rf build
	$(MAKE) doTarget MODE=clean

realclean:
	rm -rf build *.pyc
	$(MAKE) doTarget MODE=realclean

doTarget:
	SAVED_PWD=`pwd`; \
	for directory in $(SUB_DIRS); do cd $$directory; $(MAKE) $(MODE); cd $$SAVED_PWD; done
