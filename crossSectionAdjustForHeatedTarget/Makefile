# <<BEGIN-copyright>>
# <<END-copyright>>

SHELL = /bin/ksh
dirList = Doc Src Test Python Misc include lib

.PHONY: lib all python check clean realclean cycleDirList

lib:
	cd Src; $(MAKE)
	cd include; $(MAKE) 
	cd lib; $(MAKE)

all:
	$(MAKE) dirList="Src Doc Python" cycleDirsTarget="all" cycleDirList

python:
	python setup.py build

check:
	cd Test; $(MAKE) check

clean:
	rm -f *.pyc
	$(MAKE) cycleDirsTarget="clean" cycleDirList

realclean:
	rm -fr *.pyc build
	$(MAKE) cycleDirsTarget="realclean" cycleDirList

cycleDirList:
	for d in $(dirList); do \
		cd $$d; \
		$(MAKE) $(cycleDirsTarget); \
		cd ..; \
	done
