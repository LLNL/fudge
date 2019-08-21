# <<BEGIN-copyright>>
# <<END-copyright>>

SHELL = /bin/ksh

sourceDir = ../Src

.PHONY: default clean realclean

default: libnf_utilities.a

libnf_utilities.a: $(objects)
	$(AR) -cr libnf_utilities.a $(sourceDir)/*.o

clean:
	rm -f *.a

realclean: clean
