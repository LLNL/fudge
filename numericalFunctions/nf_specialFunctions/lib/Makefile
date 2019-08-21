# <<BEGIN-copyright>>
# <<END-copyright>>

SHELL = /bin/ksh

sourceDir = ../Src

.PHONY: clean realclean

libnf_specialFunctions.a: $(sourceDir)/*.c $(sourceDir)/*.h
	cd $(sourceDir); $(MAKE)
	$(AR) -cr libnf_specialFunctions.a $(sourceDir)/*.o

clean:
	rm -f *.a

realclean: clean
