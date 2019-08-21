# <<BEGIN-copyright>>
# <<END-copyright>>

SHELL = /bin/ksh

sourceDir = ../Src

.PHONY: clean realclean

libptwXY.a: $(sourceDir)/*.c $(sourceDir)/*.h
	cd $(sourceDir); $(MAKE)
	$(AR) -cr libptwXY.a $(sourceDir)/*.o

clean:
	rm -f *.a

realclean: clean
