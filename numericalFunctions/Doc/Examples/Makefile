# <<BEGIN-copyright>>
# <<END-copyright>>

SHELL = /bin/ksh

NUMERICALFUNCTIONS_DIR = ../../
Libs = -lnumericalfunctions -L$(NUMERICALFUNCTIONS_DIR)/lib -lm

.PHONY: clean realclean

E1: E1.c
	$(CC) $(CFLAGS) -I$(NUMERICALFUNCTIONS_DIR)/include -o E1 E1.c $(Libs)

clean:
	rm -f E1

realclean: clean
