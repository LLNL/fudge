# <<BEGIN-copyright>>
# <<END-copyright>>

SHELL = /bin/ksh

sourceDir = ../Src

.PHONY: clean realclean

libnf_Legendre.a: $(sourceDir)/*.c $(sourceDir)/*.h
	cd $(sourceDir); $(MAKE)
	$(AR) -cr libnf_Legendre.a $(sourceDir)/*.o

clean:
	rm -f *.a

realclean: clean
