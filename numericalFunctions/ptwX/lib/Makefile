# <<BEGIN-copyright>>
# <<END-copyright>>

SHELL = /bin/ksh

sourceDir = ../Src

.PHONY: default check clean realclean

default: libptwX.a

libptwX.a:
	cd $(sourceDir); $(MAKE)
	cp $(sourceDir)/libptwX.a .

check:

clean:
	rm -f *.a

realclean: clean
